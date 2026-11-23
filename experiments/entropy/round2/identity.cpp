#include "ring.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
namespace fs = std::filesystem;
struct Chunk {
  bool extra;
  std::string name;
};
std::vector<Chunk> files(const char *f) {
  std::ifstream in(f);
  int e;
  std::string x;
  std::vector<Chunk> out;
  while (in >> e >> x)
    out.push_back({bool(e), x});
  return out;
}
template <class T> void writevec(std::ofstream &f, const std::vector<T> &v) {
  f.write((const char *)v.data(), sizeof(T) * v.size());
}
std::vector<U> readu(std::string name) {
  std::ifstream f(name, std::ios::binary | std::ios::ate);
  if (!f)
    throw std::runtime_error(name);
  std::vector<U> a(f.tellg() / 8);
  f.seekg(0);
  f.read((char *)a.data(), a.size() * 8);
  return a;
}
using Pair = std::array<U, 2>;
struct CHash {
  size_t operator()(const Pair &a) const {
    return mix(a[0]) ^ mix(a[1] + 0x2654bc8478da3ULL);
  }
};
int main(int argc, char **argv) {
  if (argc != 4)
    return 2;
  int n = std::stoi(argv[1]);
  auto chunks = files(argv[2]);
  std::string prefix = argv[3];
  U lo = U(3) << (n / 2 - 2);
  constexpr unsigned B = 64;
  std::vector<U> repeats, badhash;
  U distinct[2]{}, maxfreq[2]{}, end[2]{}, overlap = 0;
  for (bool nodes : {true, false}) {
    std::array<std::ofstream, B> f;
    std::array<std::vector<U>, B> buf;
    for (int b = 0; b < B; ++b)
      f[b].open(prefix + "_" + std::to_string(b) + ".tmp", std::ios::binary);
    for (auto &c : chunks) {
      std::ifstream in(c.name + (nodes ? "" : ".hash"), std::ios::binary);
      U v;
      while (in.read((char *)&v, 8)) {
        unsigned b = nodes ? (v - lo) >> (n / 2 - 8) : v % 64;
        // Hash has 64 bits: retain its upper 58 bits plus the phase, implied
        // bucket carries the lower bits.
        U key = nodes ? (v << 1) | c.extra : ((v >> 6) << 1) | c.extra;
        buf[b].push_back(key);
        if (buf[b].size() >= 8192) {
          writevec(f[b], buf[b]);
          buf[b].clear();
        }
      }
    }
    for (int b = 0; b < B; ++b) {
      writevec(f[b], buf[b]);
      f[b].close();
      std::string file = prefix + "_" + std::to_string(b) + ".tmp";
      auto a = readu(file);
      std::sort(a.begin(), a.end());
      for (size_t i = 0; i < a.size();) {
        U key = a[i] >> 1;
        size_t j = i;
        U base = 0;
        while (j < a.size() && (a[j] >> 1) == key) {
          base += !(a[j] & 1);
          ++j;
        }
        U total = j - i;
        if (nodes) {
          distinct[0] += base > 0;
          ++distinct[1];
          maxfreq[0] = std::max(maxfreq[0], base);
          maxfreq[1] = std::max(maxfreq[1], total);
          end[0] += base;
          end[1] += total;
          overlap += base > 0 && total > base;
          if (total > 1)
            repeats.push_back(key);
        } else if (total > 1)
          badhash.push_back((key << 6) | b);
        i = j;
      }
      fs::remove(file);
    }
  }
  std::unordered_set<U> potential(badhash.begin(), badhash.end());
  std::unordered_map<Pair, std::array<U, 2>, CHash> duplicates;
  for (auto &c : chunks) {
    std::ifstream in(c.name, std::ios::binary);
    Pair row;
    while (in.read((char *)&row, 16))
      if (potential.count(CHash()(row))) {
        auto &cnt = duplicates[row];
        cnt[0] += !c.extra;
        ++cnt[1];
      }
  }
  U dup[2]{};
  std::vector<Pair> actual;
  for (auto &kv : duplicates) {
    for (int w = 0; w < 2; w++)
      if (kv.second[w])
        dup[w] += kv.second[w] - 1;
    if (kv.second[1] > 1)
      actual.push_back(kv.first);
  }
  std::sort(actual.begin(), actual.end());
  std::ofstream ef(prefix + ".repeated_pairs", std::ios::binary);
  writevec(ef, actual);
  std::ofstream of(prefix + ".repeated_primes", std::ios::binary);
  writevec(of, repeats);
  std::cout << "{\"n\":" << n << ",\"base_rows\":" << end[0] / 2
            << ",\"combined_rows\":" << end[1] / 2
            << ",\"base_distinct_primes\":" << distinct[0]
            << ",\"combined_distinct_primes\":" << distinct[1]
            << ",\"base_max_prime_frequency\":" << maxfreq[0]
            << ",\"combined_max_prime_frequency\":" << maxfreq[1]
            << ",\"base_repeat_product_draws\":" << dup[0]
            << ",\"combined_repeat_product_draws\":" << dup[1]
            << ",\"prime_overlap_base_supplement\":" << overlap << "}\n";
}
