#include "ring.hpp"
#include <fstream>
#include <iostream>
#include <unordered_set>
template <unsigned P, unsigned... Rest> U strip(U x) {
  while (x % P == 0)
    x /= P;
  if constexpr (sizeof...(Rest))
    return strip<Rest...>(x);
  else
    return x;
}
inline bool smooth(U v) {
  return strip<2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
               61>(v) == 1;
}
struct Pair {
  U p, q;
  bool operator==(const Pair &o) const { return p == o.p && q == o.q; }
};
struct Hash {
  size_t operator()(const Pair &x) const {
    return mix(x.p) ^ mix(x.q + 0x2654bc8478da3ULL);
  }
};
int main(int argc, char **v) {
  if (argc != 6)
    return 2;
  int n = std::stoi(v[1]);
  U lo = U(3) << (n / 2 - 2);
  W lower = W(lo) * lo, range = (W(1) << n) - lower;
  std::ifstream old(v[2], std::ios::binary);
  std::unordered_set<Pair, Hash> known;
  Pair x;
  while (old.read((char *)&x, sizeof x))
    if (x.p >= lo)
      known.insert(x);
  std::ifstream in(v[3], std::ios::binary);
  std::ofstream fl(v[4], std::ios::binary), hashes(v[5], std::ios::binary);
  while (in.read((char *)&x, sizeof x)) {
    W N = W(x.p) * x.q;
    U gap = x.q - x.p;
    uint8_t flags = 0;
    if (gap <= isqrtw(isqrtw(N << 32)))
      flags |= 1;
    for (auto ij :
         {std::pair{2, 2}, std::pair{3, 3}, std::pair{4, 3}, std::pair{4, 4},
          std::pair{5, 4}, std::pair{5, 5}, std::pair{6, 5}}) {
      W a = W(ij.first) * x.p, b = W(ij.second) * x.q,
        delta = a > b ? a - b : b - a;
      if (delta * delta < 16 * (a + b) - 64) {
        flags |= 2;
        break;
      }
    }
    if (smooth(x.p - 1) || smooth(x.q - 1))
      flags |= 4;
    if (smooth(x.p + 1) || smooth(x.q + 1))
      flags |= 8;
    if (known.count(x))
      flags |= 16;
    if ((N - lower) * 1000 < range || ((W(1) << n) - N) * 1000 < range)
      flags |= 32;
    fl.put(flags);
    const U h = Hash()(x);
    hashes.write((char *)&h, 8);
  }
  return in.eof() && fl && hashes ? 0 : 1;
}
