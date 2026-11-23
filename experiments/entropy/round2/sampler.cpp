#include <CommonCrypto/CommonDigest.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <numeric>
#include <openssl/evp.h>
#include <stdexcept>
#include <string>
#include <vector>
using U = uint64_t;
using W = __uint128_t;
using Byte = unsigned char;
const std::string key = "ent-c-R2-v1/prime-pool";
struct HashPool {
  CC_SHA256_CTX inner{}, outer{};
  HashPool() {
    Byte in[64], out[64];
    for (int j = 0; j < 64; j++) {
      Byte b = j < int(key.size()) ? key[j] : 0;
      in[j] = b ^ 0x36;
      out[j] = b ^ 0x5c;
    }
    CC_SHA256_Init(&inner);
    CC_SHA256_Update(&inner, in, 64);
    CC_SHA256_Init(&outer);
    CC_SHA256_Update(&outer, out, 64);
  }
  int operator()(U p) const {
    auto a = inner, b = outer;
    Byte d[32];
    CC_SHA256_Update(&a, &p, 8);
    CC_SHA256_Final(d, &a);
    CC_SHA256_Update(&b, d, 32);
    CC_SHA256_Final(d, &b);
    unsigned k = (unsigned(d[0]) << 8) | d[1];
    return k < 19661 ? 0 : (k < 39322 ? 1 : 2);
  }
};
struct Random {
  EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
  std::array<U, 2048> buf{};
  size_t pos = buf.size();
  explicit Random(std::string label) {
    Byte k[32], iv[16]{};
    CC_SHA256(label.data(), label.size(), k);
    if (!EVP_EncryptInit_ex(ctx, EVP_chacha20(), nullptr, k, iv))
      throw std::runtime_error("ChaCha init");
  }
  ~Random() { EVP_CIPHER_CTX_free(ctx); }
  U next() {
    if (pos == buf.size()) {
      static const std::array<Byte, 16384> zeros{};
      int len = 0;
      EVP_EncryptUpdate(ctx, reinterpret_cast<Byte *>(buf.data()), &len,
                        zeros.data(), zeros.size());
      pos = 0;
    }
    return buf[pos++];
  }
  U below(U n) {
    W x = W(next()) * n;
    if (U(x) < n) {
      U threshold = (-n) % n;
      while (U(x) < threshold)
        x = W(next()) * n;
    }
    return x >> 64;
  }
};
std::string tab(int n, int s) {
  return "data/table" + std::to_string(n) + "_" + std::to_string(s) + ".u32";
}
void table(int n) {
  U lo = U(3) << (n / 2 - 2), hi = U(1) << (n / 2);
  HashPool pool;
  std::vector<unsigned> base;
  for (unsigned i = 3; U(i) * i < hi; i += 2) {
    bool prime = true;
    for (auto p : base) {
      if (p * p > i)
        break;
      if (i % p == 0) {
        prime = false;
        break;
      }
    }
    if (prime)
      base.push_back(i);
  }
  std::array<std::ofstream, 3> files;
  for (int s = 0; s < 3; s++)
    files[s].open(tab(n, s), std::ios::binary);
  const U block = 1 << 19;
  std::vector<Byte> composite(block);
  for (U first = lo + 1; first < hi; first += 2 * block) {
    U cnt = std::min(block, (hi - first + 1) / 2);
    std::fill(composite.begin(), composite.end(), 0);
    for (U p : base) {
      U j = (p - first % p) % p;
      if (j % 2)
        j += p;
      for (j /= 2; j < cnt; j += p)
        composite[j] = 1;
    }
    std::array<std::vector<uint32_t>, 3> parts;
    for (U j = 0; j < cnt; j++)
      if (!composite[j]) {
        U p = first + 2 * j;
        parts[pool(p)].push_back(p);
      }
    for (int s = 0; s < 3; s++)
      files[s].write(reinterpret_cast<char *>(parts[s].data()),
                     parts[s].size() * 4);
  }
}
struct Prime {
  Random rng;
  HashPool pool;
  int stage, n;
  std::vector<uint32_t> primes;
  std::vector<unsigned> wheel;
  U lo, hi, first, last;
  int (*isprime)(unsigned long);
  Prime(int n_, int s, std::string seed) : rng(seed), stage(s), n(n_) {
    lo = U(3) << (n / 2 - 2);
    hi = U(1) << (n / 2);
    if (n <= 64) {
      std::ifstream f(tab(n, s), std::ios::binary | std::ios::ate);
      if (!f)
        throw std::runtime_error("prime table missing");
      primes.resize(f.tellg() / 4);
      f.seekg(0);
      f.read(reinterpret_cast<char *>(primes.data()), primes.size() * 4);
      return;
    }
    const char *lib = getenv("FLINTLIB");
    void *handle = dlopen(lib, RTLD_NOW | RTLD_LOCAL);
    if (!handle)
      throw std::runtime_error(dlerror());
    isprime = (int (*)(unsigned long))dlsym(handle, "n_is_prime");
    if (!isprime)
      throw std::runtime_error("FLINT primality missing");
    for (unsigned i = 1; i < 30030; i += 2)
      if (std::gcd(i, 30030u) == 1)
        wheel.push_back(i);
    first = lo / 30030;
    last = (hi - 1) / 30030;
  }
  U draw() {
    if (!primes.empty())
      return primes[rng.below(primes.size())];
    for (;;) {
      U x = (first + rng.below(last - first + 1)) * 30030 +
            wheel[rng.below(wheel.size())];
      if (x < lo || x >= hi || x % 17 == 0 || x % 19 == 0 || x % 23 == 0 ||
          x % 29 == 0 || x % 31 == 0 || x % 37 == 0 || x % 41 == 0 ||
          x % 43 == 0 || x % 47 == 0)
        continue;
      if (isprime(x) && pool(x) == stage)
        return x;
    }
  }
};
int main(int argc, char **argv) {
  if (argc == 3 && std::string(argv[1]) == "table") {
    table(std::stoi(argv[2]));
    return 0;
  }
  if (argc != 6 && argc != 7)
    return 2;
  bool equal8 = argc == 7 && std::string(argv[6]) == "equal8";
  int n = std::stoi(argv[1]), s = std::stoi(argv[2]);
  U m = std::stoull(argv[3]);
  Prime src(n, s,
            "ent-c-R2-v1/data/" + std::to_string(n) + "/" + std::to_string(s) +
                "/" + argv[4]);
  std::ofstream f(argv[5], std::ios::binary);
  std::vector<U> buf;
  buf.reserve(65536);
  while (m) {
    U p = src.draw(), q = src.draw();
    if (p == q || (equal8 && (p % 8 != q % 8)))
      continue;
    if (p > q)
      std::swap(p, q);
    buf.push_back(p);
    buf.push_back(q);
    --m;
    if (buf.size() == 65536 || !m) {
      f.write(reinterpret_cast<char *>(buf.data()), buf.size() * 8);
      buf.clear();
    }
  }
  return f ? 0 : 1;
}
