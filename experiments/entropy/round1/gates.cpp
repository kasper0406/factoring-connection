#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <gmpxx.h>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using Z = mpz_class;
using U = uint64_t;
using W = __uint128_t;
constexpr int G = 50;
Z big(W x) {
  U v[2] = {U(x), U(x >> 64)};
  Z z;
  mpz_import(z.get_mpz_t(), 2, -1, 8, 0, 0, v);
  return z;
}
W word(const Z &z) {
  auto a = z.get_mpz_t();
  W v = a->_mp_size ? a->_mp_d[0] : 0;
  if (a->_mp_size >= 2)
    v |= W(a->_mp_d[1]) << 64;
  return v;
}
Z powmod(const Z &x, const Z &e, const Z &n) {
  Z v;
  mpz_powm(v.get_mpz_t(), x.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
  return v;
}
Z root(const Z &x) {
  Z r;
  mpz_sqrt(r.get_mpz_t(), x.get_mpz_t());
  return r;
}

enum {
  GATE = 1,
  LADDER = 2,
  CF = 4,
  FERMAT = 8,
  LEHMAN = 16,
  PM1 = 32,
  PP1 = 64
};
struct Audit {
  U p, q;
  unsigned flags = 0;
  void check(const Z &v, unsigned kind) {
    if (!p || (flags & kind) || v == 0)
      return;
    if (v < 0) {
      Z t = -v;
      check(t, kind);
      return;
    }
    W w = word(v);
    if (w < p)
      return;
    if ((w % p == 0) != (w % q == 0))
      flags |= kind;
  }
};
// Exact Fermat/Lehman window membership from the divisors of kN.
// A^2-4kN=B^2 and a proper split imply A=i*p+j*q, i*j=k.
bool in_window(W u, W v, U width) {
  W delta = u > v ? u - v : v - u, sum = u + v;
  return delta * delta < W(2) * width * sum - W(width) * width;
}
const std::array<unsigned, 31> powers = {
    128, 81, 125, 49, 121, 13, 17, 19, 23, 29,  31,  37,  41,  43,  47, 53,
    59,  61, 67,  71, 73,  79, 83, 89, 97, 101, 103, 107, 109, 113, 127};
U mul(U x, U y, U p) { return W(x) * y % p; }
U smallpow(U x, U e, U p) {
  U v = 1;
  while (e) {
    if (e & 1)
      v = mul(v, x, p);
    e >>= 1;
    if (e)
      x = mul(x, x, p);
  }
  return v;
}
U lucasprime(U x, U e, U p) {
  U a = 2, b = x;
  for (int j = 63 - __builtin_clzll(e); j >= 0; --j) {
    U mid = (mul(a, b, p) + p - x) % p;
    if ((e >> j) & 1) {
      b = (mul(b, b, p) + p - 2) % p;
      a = mid;
    } else {
      a = (mul(a, a, p) + p - 2) % p;
      b = mid;
    }
  }
  return a;
}
int first_stage(U prime, bool plus) {
  // Order divides p-1 for base 2, p-(5|p) for roots of X²-3X+1.
  U order = prime - 1;
  if (plus && (prime % 5 == 2 || prime % 5 == 3))
    order = prime + 1;
  U smooth = 1;
  for (auto e : powers) {
    U d = std::gcd<U>(order, e);
    smooth *= d;
    order /= d;
  }
  U value = plus ? lucasprime(3, smooth, prime) : smallpow(2, smooth, prime);
  if (value != (plus ? 2 : 1))
    return powers.size();
  U a = plus ? 3 : 2;
  for (size_t j = 0; j < powers.size(); ++j) {
    a = plus ? lucasprime(a, powers[j], prime) : smallpow(a, powers[j], prime);
    if (a == (plus ? 2 : 1))
      return j;
  }
  throw std::logic_error("stage-order audit");
}
void split_methods(Audit &audit) {
  U p = audit.p, q = audit.q;
  if (!p)
    return;
  for (unsigned i = 1; i <= 32 && !(audit.flags & LEHMAN); ++i)
    for (unsigned j = 1; i * j <= 32; ++j)
      if (i * j >= 2 && in_window(W(i) * p, W(j) * q, 8)) {
        audit.flags |= LEHMAN;
        break;
      }
  // For odd nonsquare N, 2*ceil(sqrt(N)) is the nearest even integer above
  // 2sqrt(N). The exact 256-step condition uses a^2 vs N at a=(p+q)/2-256.
  W a = (W(p) + q) / 2;
  bool fermat = (a <= 256 || (a - 256) * (a - 256) < W(p) * q);
  if (fermat)
    audit.flags |= FERMAT;
  else
    audit.flags &= ~FERMAT;
  for (bool plus : {false, true})
    if (first_stage(p, plus) != first_stage(q, plus))
      audit.flags |= plus ? PP1 : PM1;
}
Z lucas(const Z &N, int P, int Q, Audit &audit) {
  Z a = 2, b = P, mid, sq;
  int sign = 1;
  for (int j = mpz_sizeinbase(N.get_mpz_t(), 2) - 1; j >= 0; --j) {
    mid = (a * b - P * sign) % N;
    if (mid < 0)
      mid += N;
    if (mpz_tstbit(N.get_mpz_t(), j)) {
      sq = (b * b - 2 * sign * Q) % N;
      if (sq < 0)
        sq += N;
      a.swap(mid);
      b.swap(sq);
      sign = Q;
    } else {
      sq = (a * a - 2 * sign) % N;
      if (sq < 0)
        sq += N;
      b.swap(mid);
      a.swap(sq);
      sign = 1;
    }
    audit.check(a, LADDER);
    audit.check(b, LADDER);
  }
  return a;
}
void evaluate(const Z &N, int n, Audit &au, uint16_t *out) {
  int idx = 0;
  auto add = [&](const Z &z, bool modular = false) {
    au.check(z, GATE);
    if (modular) {
      au.check(z - 1, GATE);
      au.check(z + 1, GATE);
    }
    if (out) {
      unsigned w = z.get_ui(), jacbit = 0;
      int j = mpz_jacobi(z.get_mpz_t(), N.get_mpz_t());
      jacbit = (unsigned(j < 0) << 11) | (unsigned(j == 0) << 12);
      out[idx] = (w & 255) | (__builtin_parity(w & 255) << 8) |
                 (__builtin_parity(w & 65535) << 9) |
                 ((mpz_popcount(z.get_mpz_t()) & 1) << 10) | jacbit;
    }
    ++idx;
  };
  add(N);
  Z e = (N - 1) / 2;
  au.check(e, GATE);
  au.check(N - 1, GATE);
  au.check(N + 1, GATE);
  std::array<Z, 2> inner;
  int xi = 0;
  for (int x : {2, 3, 5, 7}) {
    Z t = powmod(Z(x), e, N), minus = (t * t) % N, v = (minus * x) % N,
      plus = (v * x) % N;
    add(v, true);
    add(t, true);
    add(minus, true);
    add(plus, true);
    if (xi < 2)
      inner[xi] = v;
    ++xi;
  }
  for (auto pq : {std::pair{3, 1}, std::pair{5, 1}, std::pair{3, -1}})
    add(lucas(N, pq.first, pq.second, au), true);
  for (auto &x : inner)
    for (int c : {-1, 1}) {
      Z v = x + c;
      au.check(v, GATE);
      add(powmod(v, N, N), true);
    }
  for (int k : {3, 5, 7, 11, 13, 16, 31}) {
    add(N / k);
    add(N % k);
  }
  Z a0 = root(N);
  add(a0);
  Z fixed = root(N << 32);
  au.check(fixed, GATE);
  add(fixed % 65536);
  for (int d : {2, 3}) {
    fixed = root((N * N * d) << 32);
    // Both scaled root and integer part may exceed 128 bits only for larger,
    // unsupported n.
    au.check(fixed, GATE);
    au.check(fixed >> 16, GATE);
    add(fixed % 65536);
  }
  Z m = 0, den = 1, a = a0;
  for (int i = 0; i < 4; ++i) {
    m = den * a - m;
    au.check(m, CF);
    den = (N - m * m) / den;
    au.check(den, CF);
    if (den == 0)
      throw std::logic_error("square witness");
    Z num = a0 + m;
    au.check(num, CF);
    a = num / den;
    add(a);
  }
  Z numer = N;
  den = Z(1) << (n / 2);
  for (int i = -1; i < 4; ++i) {
    Z quot = 0, rem = 0;
    if (den != 0) {
      mpz_fdiv_qr(quot.get_mpz_t(), rem.get_mpz_t(), numer.get_mpz_t(),
                  den.get_mpz_t());
      au.check(rem, CF);
    }
    if (i < 0)
      au.check(quot, CF);
    else
      add(quot);
    numer.swap(den);
    den.swap(rem);
  }
  if (idx != G)
    throw std::logic_error("gate count");
  split_methods(au);
}
int main(int argc, char **argv) {
  if (argc != 6)
    return 2;
  int n = std::stoi(argv[1]);
  bool witness = std::string(argv[2]) == "witness";
  std::ifstream in(argv[3], std::ios::binary);
  std::ofstream features(argv[4], std::ios::binary),
      flags(argv[5], std::ios::binary);
  std::array<U, 4> r;
  while (in.read(reinterpret_cast<char *>(r.data()), sizeof(r))) {
    uint16_t f[G];
    uint8_t mask[2]{};
    Audit au{witness ? 0 : r[0], witness ? 0 : r[1]};
    Z N = big(witness ? (W(r[1]) << 64 | r[0]) : W(r[0]) * r[1]);
    evaluate(N, n, au, f);
    mask[0] = au.flags;
    if (!witness) {
      Audit ref{r[2], r[3]};
      evaluate(big(W(r[2]) * r[3]), n, ref, nullptr);
      mask[1] = ref.flags;
    }
    features.write(reinterpret_cast<char *>(f), sizeof(f));
    flags.write(reinterpret_cast<char *>(mask), 2);
  }
  return (in.eof() && features && flags) ? 0 : 1;
}
