#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <gmpxx.h>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using U = uint64_t;
using W = __uint128_t;
using Z = mpz_class;
inline Z big(W w) {
  unsigned long l[2] = {static_cast<unsigned long>(w),
                        static_cast<unsigned long>(w >> 64)};
  Z a;
  mpz_import(a.get_mpz_t(), 2, -1, 8, 0, 0, l);
  return a;
}
inline W small(const Z &a) {
  return W(mpz_getlimbn(a.get_mpz_t(), 1)) << 64 |
         mpz_getlimbn(a.get_mpz_t(), 0);
}
inline Z sqroot(const Z &a) {
  Z v;
  mpz_sqrt(v.get_mpz_t(), a.get_mpz_t());
  return v;
}
inline U isqrtw(W a) {
  U x = static_cast<U>(std::sqrt((long double)a));
  while (W(x) * x > a)
    --x;
  while (W(x + 1) * (x + 1) <= a)
    ++x;
  return x;
}
inline bool square64(U a) {
  unsigned r = a & 15;
  if (r != 0 && r != 1 && r != 4 && r != 9)
    return false;
  U s = isqrtw(a);
  return s * s == a;
}
inline int jac(const Z &a, const Z &n) {
  return mpz_jacobi(a.get_mpz_t(), n.get_mpz_t());
}
inline U mix(U x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
struct Trace {
  U p, q;
  bool split = false;
  void see(W z) {
    if (!split && z >= p && ((z % p == 0) != (z % q == 0)))
      split = true;
  }
  void seez(const Z &z) {
    if (split || z == 0)
      return;
    if (mpz_size(z.get_mpz_t()) <= 2)
      see(small(z));
    else if (bool(mpz_divisible_ui_p(z.get_mpz_t(), p)) !=
             bool(mpz_divisible_ui_p(z.get_mpz_t(), q)))
      split = true;
  }
};
struct Ring {
  W N;
  Z nz;
  Trace &tr;
  bool singular = false;
  W watched(W x) {
    tr.see(x);
    return x;
  }
  W add(W a, W b) {
    W v = a + b;
    return watched(v >= N ? v - N : v);
  }
  W sub(W a, W b) { return watched(a >= b ? a - b : N - (b - a)); }
  W mul(W a, W b) {
    if ((N >> 64) == 0)
      return watched(W(U(a)) * U(b) % U(N));
    mp_limb_t aa[2] = {mp_limb_t(a), mp_limb_t(a >> 64)},
              bb[2] = {mp_limb_t(b), mp_limb_t(b >> 64)},
              den[2] = {mp_limb_t(N), mp_limb_t(N >> 64)}, prod[4], quot[3],
              rem[2];
    mpn_mul_n(prod, aa, bb, 2);
    mpn_tdiv_qr(quot, rem, 0, prod, 4, den, 2);
    return watched(W(rem[1]) << 64 | rem[0]);
  }
  W power(W x, W e) {
    tr.see(x);
    Z v = big(x), ex = big(e);
    mpz_powm(v.get_mpz_t(), v.get_mpz_t(), ex.get_mpz_t(), nz.get_mpz_t());
    return watched(small(v));
  }
  W inv(W x) {
    tr.see(x);
    Z a = big(x), v;
    if (!mpz_invert(v.get_mpz_t(), a.get_mpz_t(), nz.get_mpz_t())) {
      singular = true;
      return 0;
    }
    return watched(small(v));
  }
};
using Out = std::array<Z, 4>;
inline Out output(W a, W b = 0, W c = 0, W d = 0) {
  return {big(a), big(b), big(c), big(d)};
}
inline W lcmpow(int b) {
  Z m = 1;
  for (int j = 2; j <= b; ++j)
    mpz_lcm_ui(m.get_mpz_t(), m.get_mpz_t(), j);
  return small(m);
}
inline W M(int b) {
  static const W m16 = lcmpow(16), m32 = lcmpow(32), m64 = lcmpow(64);
  return b == 16 ? m16 : b == 32 ? m32 : m64;
}
inline unsigned bitlen(W x) {
  return x >> 64 ? 128 - __builtin_clzll(U(x >> 64))
                 : (x ? 64 - __builtin_clzll(U(x)) : 0);
}
