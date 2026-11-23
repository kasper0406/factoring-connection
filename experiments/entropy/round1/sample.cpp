#include <array>
#include <cstdint>
#include <fstream>
#include <gmpxx.h>
#include <iostream>
#include <random>
#include <vector>
using U = uint64_t;
using W = __uint128_t;

struct PrimeSource {
  int h;
  U low;
  std::vector<unsigned> primes;
  explicit PrimeSource(int bits) : h(bits), low(U(1) << (h - 1)) {
    if (h > 24)
      return;
    std::vector<bool> composite(low);
    for (U i = 3; i * i < 2 * low; i += 2)
      if (!composite[i / 2])
        for (U j = i * i; j < 2 * low; j += 2 * i)
          composite[j / 2] = true;
    for (U i = low + 1; i < 2 * low; i += 2)
      if (!composite[i / 2])
        primes.push_back(i);
  }
  U draw(std::mt19937_64 &rng) const {
    if (!primes.empty())
      return primes[std::uniform_int_distribution<size_t>(0, primes.size() -
                                                                 1)(rng)];
    mpz_class z;
    for (;;) {
      U x = low | (rng() & (low - 1)) | 1;
      mpz_set_ui(z.get_mpz_t(), x);
      int cert = mpz_probab_prime_p(z.get_mpz_t(), 25);
      if (cert == 2)
        return x;
      if (cert == 1) { // Deterministic for unsigned 64-bit inputs.
        U d = x - 1;
        int s = 0;
        while (!(d & 1)) {
          d >>= 1;
          ++s;
        }
        bool prime = true;
        for (U a :
             {2ul, 325ul, 9375ul, 28178ul, 450775ul, 9780504ul, 1795265022ul}) {
          if (a % x == 0)
            continue;
          mpz_class v(static_cast<unsigned long>(a));
          mpz_powm_ui(v.get_mpz_t(), v.get_mpz_t(), d, z.get_mpz_t());
          U t = v.get_ui();
          bool ok = t == 1 || t == x - 1;
          for (int j = 1; j < s && !ok; ++j) {
            t = W(t) * t % x;
            ok = t == x - 1;
          }
          if (!ok) {
            prime = false;
            break;
          }
        }
        if (prime)
          return x;
      }
    }
  }
  int pair(std::mt19937_64 &rng, U &p, U &q) const {
    for (;;) {
      p = draw(rng);
      q = draw(rng);
      W v = W(p) * q;
      if (p == q || !(v >> (2 * h - 1)))
        continue;
      if (p > q)
        std::swap(p, q);
      return int(v >> (2 * h - 8)) - 128;
    }
  }
};
int main(int argc, char **argv) {
  if (argc != 6)
    return 2;
  int n = std::stoi(argv[1]), count = std::stoi(argv[2]);
  PrimeSource primes(n / 2);
  std::mt19937_64 primary(std::stoull(argv[3])),
      reference(std::stoull(argv[4]));
  std::array<std::vector<int>, 128> bins;
  std::vector<std::array<U, 4>> rows(count);
  for (int i = 0; i < count; ++i) {
    auto &r = rows[i];
    bins[primes.pair(primary, r[0], r[1])].push_back(i);
  }
  std::array<size_t, 128> used{};
  int filled = 0;
  while (filled < count) {
    U p, q;
    int b = primes.pair(reference, p, q);
    if (used[b] == bins[b].size())
      continue;
    auto &r = rows[bins[b][used[b]++]];
    r[2] = p;
    r[3] = q;
    ++filled;
  }
  std::ofstream out(argv[5], std::ios::binary);
  out.write(reinterpret_cast<const char *>(rows.data()),
            rows.size() * sizeof(rows[0]));
  return out ? 0 : 1;
}
