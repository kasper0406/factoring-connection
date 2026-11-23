#include <cstdint>
using I = uint32_t;
extern "C" double accumulate(int64_t len, const I *nodes, const double *e,
                             double *sums, int sparse) {
  double diff = 0;
  if (sparse) {
    for (int64_t i = 0; i < len; i++) {
      I row = nodes[2 * i], node = nodes[2 * i + 1];
      double z = e[row];
      diff += 2 * sums[node] * z;
      sums[node] += z;
    }
  } else {
    for (int64_t i = 0; i < len; i++) {
      double z = e[i];
      if (!z)
        continue;
      for (int k = 0; k < 2; k++) {
        I node = nodes[2 * i + k];
        if (node) {
          diff += 2 * sums[node] * z;
          sums[node] += z;
        }
      }
    }
  }
  return diff;
}
