#include <algorithm>
#include <cstdint>
#include <vector>
using U = uint64_t;
using I = int32_t;
extern "C" void matched(int words, int ntests, int ntarget, int nu,
                        const I *spec, const U *x, const U *y, uint32_t *counts,
                        uint32_t *marginal) {
  std::vector<U> plus(words), minus(words);
  for (int t = 0; t < ntests; ++t) {
    const I *row = spec + 4 * t;
    int op = row[0];
    unsigned ones = 0;
    for (int i = 0; i < words; ++i) {
      auto value = [&](int member) {
        const U *u = x + size_t(member) * nu * words + i;
        if (!op)
          return u[size_t(t) * words];
        U a = u[size_t(row[1]) * words], b = u[size_t(row[2]) * words];
        if (op == 1)
          return a ^ b;
        U c = u[size_t(row[3]) * words];
        return (a & b) | (a & c) | (b & c);
      };
      U a = value(0), b = value(1);
      plus[i] = a & ~b;
      minus[i] = b & ~a;
      ones += __builtin_popcountll(a) + __builtin_popcountll(b);
    }
    marginal[t] = ones;
    for (int j = 0; j < ntarget; ++j) {
      const U *a = y + size_t(2 * j) * words;
      const U *b = a + words;
      unsigned pos = 0, neg = 0;
      for (int k = 0; k < words; ++k) {
        pos += __builtin_popcountll((a[k] & plus[k]) | (b[k] & minus[k]));
        neg += __builtin_popcountll((a[k] & minus[k]) | (b[k] & plus[k]));
      }
      counts[(size_t(t) * ntarget + j) * 2] = pos;
      counts[(size_t(t) * ntarget + j) * 2 + 1] = neg;
    }
  }
}
