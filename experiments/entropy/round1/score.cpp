#include <cstdint>
#include <vector>
extern "C" void measure(int words, int targets, int views, int pairs,
                        const int32_t *tasks, const int32_t *ops,
                        const uint64_t *unary, const uint64_t *masks,
                        uint32_t *counts) {
  std::vector<uint64_t> tmp(words);
  int previous = -1;
  const uint64_t *values = nullptr;
  for (int k = 0; k < pairs; ++k) {
    int t = tasks[2 * k], y = tasks[2 * k + 1];
    if (t != previous) {
      int op = ops[3 * t], a = ops[3 * t + 1], b = ops[3 * t + 2];
      if (!op)
        values = unary + t * words;
      else {
        const uint64_t *left = unary + a * words;
        const uint64_t *right = unary + b * words;
        if (op == 1)
          for (int j = 0; j < words; ++j)
            tmp[j] = left[j] ^ right[j];
        if (op == 2)
          for (int j = 0; j < words; ++j)
            tmp[j] = left[j] & right[j];
        if (op == 3)
          for (int j = 0; j < words; ++j)
            tmp[j] = left[j] | right[j];
        values = tmp.data();
      }
      previous = t;
    }
    for (int v = 0; v < views; ++v)
      for (int label = 0; label < 2; ++label) {
        const uint64_t *m = masks + ((y * views + v) * 2 + label) * words;
        unsigned total = 0;
        for (int j = 0; j < words; ++j)
          total += __builtin_popcountll(values[j] & m[j]);
        counts[(k * views + v) * 2 + label] = total;
      }
  }
}
