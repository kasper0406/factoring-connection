#include "ring.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
struct H {
  size_t operator()(U x) const { return mix(x); }
};
using I = uint32_t;
int main(int argc, char **argv) {
  if (argc != 4)
    return 2;
  int n = std::stoi(argv[1]);
  U lo = U(3) << (n / 2 - 2);
  std::ifstream input(std::string(argv[3]) + ".repeated_primes",
                      std::ios::binary);
  std::vector<U> bitmap;
  std::vector<I> prefix;
  std::unordered_map<U, I, H> rep;
  if (n <= 64)
    bitmap.resize((U(1) << (n / 2 - 3)) / 64);
  U p;
  I count = 0;
  while (input.read((char *)&p, 8)) {
    ++count;
    if (n <= 64) {
      U id = (p - lo) / 2;
      bitmap[id / 64] |= U(1) << (id % 64);
    } else
      rep[p] = count;
  }
  if (n <= 64) {
    prefix.resize(bitmap.size());
    I acc = 0;
    for (size_t i = 0; i < bitmap.size(); i++) {
      prefix[i] = acc;
      acc += __builtin_popcountll(bitmap[i]);
    }
  }
  auto id = [&](U p) -> I {
    if (n > 64) {
      auto it = rep.find(p);
      return it == rep.end() ? 0 : it->second;
    }
    U i = (p - lo) / 2, b = bitmap[i / 64], mask = U(1) << (i % 64);
    return (b & mask)
               ? 1 + prefix[i / 64] + __builtin_popcountll(b & (mask - 1))
               : 0;
  };
  std::ifstream list(argv[2]);
  int extra;
  std::string name;
  while (list >> extra >> name) {
    std::ifstream a(name, std::ios::binary);
    std::ofstream out(name + ".nodes", std::ios::binary);
    U pair[2];
    I i = 0;
    while (a.read((char *)pair, 16)) {
      I v[2] = {id(pair[0]), id(pair[1])};
      if (n <= 64)
        out.write((char *)v, 8);
      else
        for (I node : v)
          if (node) {
            I entry[2] = {i, node};
            out.write((char *)entry, 8);
          }
      ++i;
    }
  }
}
