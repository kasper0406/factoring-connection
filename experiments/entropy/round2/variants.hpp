#pragma once
#include "geometry.hpp"
inline const std::vector<std::string> names = {
    "G1/summed_inner", "G1/summed_outer",
    "G1/shift_tower",  "G1/inverse",
    "G1/SLP16",        "G1/SLP32",
    "G1/SLP64",        "G1/SLP_Fermat",
    "G1/cubic",        "G1/quartic",
    "G1/AKS5",         "G1/asymmetric_Lucas",
    "G2/k1_pow2",      "G2/k3_pow2",
    "G2/k5_pow2",      "G2/k1_powN",
    "G3/cf1_16",       "G3/cf3_16",
    "G3/cf5_64",       "G3/cf7_64",
    "G3/Euclid",       "G3/infra_fourth",
    "G3/infra_sqrt",   "G4/Frob_t1_ell3",
    "G4/Frob_t3_ell5", "G4/37a1_psi3",
    "G4/Hurwitz_edge", "G5/Coppersmith3",
    "G5/Coppersmith4", "G5/GAR3",
    "G5/GAR4",         "G6/N",
    "G6/Nplus1",       "G6/Hasse_edge",
    "G6/SLP32",        "G6/Ntwist",
    "G6/CM1728",       "G6/nested",
    "G7/QR_dense",     "G7/digit_sums",
    "G7/automata",     "G7/carries",
    "G7/public_floor", "G6/37a1_N",
    "G6/37a1_Nplus1",  "G1/Redei",
    "G3/cf11_64",      "G3/cf13_16",
    "G5/Coppersmith5", "G5/GAR5",
    "G1/AKS7",         "G2/k7_pow2",
    "G1/inner_mixed"};
inline Out variant(int v, Ring &r) {
  W N = r.N;
  if (v == 43 || v == 44)
    return ec37(r, N + (v == 44));
  if (v == 45) {
    W x = r.power(2, N);
    W t = r.mul(r.add(r.mul(x, x), 3), r.inv(r.mul(2, x)));
    return output(r.power(r.add(t, 1), N), t, r.sub(t, 1), r.add(t, 1));
  }
  if (v == 46 || v == 47)
    return cf(r, v == 46 ? 11 : 13, v == 46 ? 64 : 16);
  if (v == 48 || v == 49)
    return lattice(r, v == 49, 5);
  if (v == 50) {
    Poly base{}, rel{};
    base[0] = 2;
    base[1] = 1;
    rel[0] = 1;
    auto a = polypow(r, base, rel, 7, N);
    return output(a[0], a[1], a[3], a[6]);
  }
  if (v == 51)
    return forms(r, 7, false);
  if (v == 52) {
    W x = r.power(5, N - 1), y = r.power(7, N + 1);
    W t = r.power(r.add(r.add(x, y), 3), N);
    return output(t, r.sub(t, 5), r.add(t, y), r.sub(t, x));
  }
  if (v <= 6) {
    W a = r.power(2, N), b = 0, z = 0, y = 0;
    switch (v) {
    case 0:
      b = r.power(3, (N - 1) / 2);
      z = r.power(r.add(r.add(a, b), 1), N);
      break;
    case 1:
      b = r.power(3, N);
      z = r.add(r.power(r.add(a, 1), N), r.power(r.sub(b, 1), N));
      break;
    case 2:
      b = r.power(3, N);
      y = r.power(r.add(a, 1), N);
      z = r.sub(r.power(r.add(y, 2), N), b);
      break;
    case 3:
      y = r.inv(r.add(a, 1));
      z = r.power(r.add(y, 2), N);
      break;
    default:
      y = r.power(r.add(a, 1), M(v == 4 ? 16 : v == 5 ? 32 : 64));
      z = r.power(r.add(y, 3), N);
      break;
    }
    return output(z, r.sub(z, 3), v >= 4 ? r.sub(y, 1) : r.add(z, 1),
                  r.sub(z, a));
  }
  if (v == 7) {
    W x = 2;
    for (int j = 0; j < 10; ++j)
      x = r.mul(x, x);
    return output(r.add(x, 3), r.sub(x, 1), r.add(x, 1), x);
  }
  if ((v >= 8 && v <= 11) || (v >= 23 && v <= 25)) {
    int d = v == 8 ? 3 : v == 9 || v == 25 ? 4 : v == 10 ? 5 : 2;
    Poly rel{}, base{};
    base[1] = 1;
    W e = N;
    if (v == 10) {
      rel[0] = 1;
      base[0] = 2;
    } else if (v == 11) {
      rel[0] = N - 3;
      rel[1] = 5;
    } else if (v == 23 || v == 24) {
      rel[0] = N - (v == 23 ? 3 : 5);
      rel[1] = v == 23 ? 1 : 3;
      e += v == 23 ? 1 : -1;
    } else if (v == 25) {
      rel[0] = r.inv(3);
      rel[1] = N - 1;
      rel[2] = 2;
    } else {
      rel[0] = rel[1] = 1;
    }
    auto a = polypow(r, base, rel, d, e);
    if (v == 23 || v == 24) {
      int ell = v == 23 ? 3 : 5;
      unsigned fiber = 0;
      for (int t = 0; t * t < 4 * ell; ++t)
        fiber = 3 * fiber + jac(Z(t * t - 4 * ell), r.nz) + 1;
      return output(a[0], a[1], fiber, r.sub(a[1], 1));
    }
    return output(a[0], a[1], d > 2 ? a[2] : r.add(a[0], a[1]),
                  d > 3 ? a[3] : r.sub(a[0], a[1]));
  }
  if (v >= 12 && v <= 15)
    return forms(r, v == 13 ? 3 : v == 14 ? 5 : 1, v == 15);
  if (v >= 16 && v <= 19)
    return cf(r, v == 16 ? 1 : v == 17 ? 3 : v == 18 ? 5 : 7, v < 18 ? 16 : 64);
  if (v == 20)
    return euclid(r);
  if (v == 21 || v == 22)
    return infrastructure(r, v == 22);
  if (v == 26) {
    U s = isqrtw(4 * N);
    unsigned jacc = 0;
    for (int j = 0; j < 8; ++j) {
      Z delta = big(W(s - j) * (s - j)) - 4 * r.nz;
      r.tr.seez(delta);
      jacc = 3 * jacc + jac(delta, r.nz) + 1;
    }
    Z delta = big(W(s) * s) - 4 * r.nz;
    return {delta, Z(jacc), big(s), big(N % 3840)};
  }
  if (v >= 27 && v <= 30)
    return lattice(r, v >= 29, v == 28 || v == 30 ? 4 : 3);
  if (v >= 31 && v <= 37)
    return ec(r, v - 31);
  if (v == 38) {
    Z values[] = {r.nz, r.nz * r.nz, (Z(1) << (2 * bitlen(N))) / r.nz};
    unsigned out[3]{};
    int primes[] = {7, 11, 19, 23, 31, 43, 47, 59};
    for (int j = 0; j < 3; ++j) {
      int len = mpz_sizeinbase(values[j].get_mpz_t(), 2);
      for (int k = 0; k < 8; ++k) {
        U mask = 0;
        int p = primes[k];
        for (int h = 1; h < p; ++h)
          mask |= U(1) << (h * h % p);
        int bit = 0;
        for (int h = 0; h < len; ++h)
          if ((mask >> (h % p)) & 1)
            bit ^= mpz_tstbit(values[j].get_mpz_t(), h);
        out[j] |= bit << k;
      }
    }
    return output(out[0], out[1], out[2], out[0] ^ out[2]);
  }
  if (v == 39) {
    Out out;
    for (int j = 0; j < 4; ++j) {
      unsigned b = (int[]){3, 5, 7, 3}[j];
      W a = N;
      int s = 0;
      while (a) {
        int d = a % b;
        a /= b;
        if (j == 3 && d == 2) {
          d = -1;
          ++a;
        }
        s += d;
      }
      out[j] = s;
    }
    return out;
  }
  if (v == 40) {
    unsigned a = 0, b = 1, c = 2;
    unsigned A = 1, B = 0, C = 0, D = 1;
    for (int i = bitlen(N) - 1; i >= 0; --i) {
      if ((N >> i) & 1) {
        unsigned t = a;
        a = b;
        b = c;
        c = t;
        C = (C + A) % 3;
        D = (D + B) % 3;
      } else {
        std::swap(a, b);
        A = (A + C) % 3;
        B = (B + D) % 3;
      }
    }
    return output(9 * a + 3 * b + c, 27 * A + 9 * B + 3 * C + D, A + D, a);
  }
  if (v == 41) {
    W a = r.power(2, N), b = r.power(3, N);
    unsigned longest = 0, run = 0, total = 0, carry = 0;
    for (unsigned j = 0; j < bitlen(N); ++j) {
      carry = (((a >> j) & 1) + ((b >> j) & 1) + carry) >= 2;
      run = carry ? run + 1 : 0;
      longest = std::max(run, longest);
      total += carry;
    }
    return output(longest, total, r.add(a, b), a ^ b);
  }
  U s = isqrtw(N);
  return output(N % 256, N % 3, N % 5, s);
}
inline uint32_t projection(const Z &z, Ring &r) {
  r.tr.seez(z);
  Z a = abs(z);
  U lo = a.get_ui();
  Z rep;
  mpz_mod(rep.get_mpz_t(), z.get_mpz_t(), r.nz.get_mpz_t());
  unsigned j = jac(z, r.nz) + 1;
  W w = small(rep);
  U frac = (w * 8) / r.N;
  return (lo & 63) | (__builtin_parity(lo & 255) << 6) |
         ((mpz_popcount(a.get_mpz_t()) & 1) << 7) | ((j == 0) << 8) |
         ((j == 1) << 9) | (frac << 10) |
         ((mpz_popcount(rep.get_mpz_t()) > bitlen(r.N) / 2) << 13) |
         ((z < 0) << 14) |
         ((__builtin_parityll(U(w) & 0xd5a65bb99d13a6adULL) ^
           __builtin_parityll(U(w >> 64) & 0xe1990aafa5babc65ULL))
          << 15);
}
