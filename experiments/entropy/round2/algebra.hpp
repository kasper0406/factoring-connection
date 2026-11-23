#pragma once
#include "ring.hpp"
#include <dlfcn.h>
struct Flint {
  void *lib = dlopen(getenv("FLINTLIB"), RTLD_LOCAL | RTLD_NOW);
  template <class F> F fn(const char *name) {
    if (!lib)
      throw std::runtime_error(dlerror());
    auto p = dlsym(lib, name);
    if (!p)
      throw std::runtime_error(name);
    return reinterpret_cast<F>(p);
  }
  using F = long;
  struct Form {
    F a, b, c;
  };
  void (*set)(F *, mpz_srcptr) = fn<decltype(set)>("fmpz_set_mpz");
  void (*get)(mpz_ptr, const F *) = fn<decltype(get)>("fmpz_get_mpz");
  void (*clear)(F *) = fn<decltype(clear)>("fmpz_clear");
  void (*initf)(Form *) = fn<decltype(initf)>("qfb_init");
  void (*clearf)(Form *) = fn<decltype(clearf)>("qfb_clear");
  void (*dupl)(Form *, const Form *, const F *,
               const F *) = fn<decltype(dupl)>("qfb_nudupl");
  void (*comp)(Form *, const Form *, const Form *, const F *,
               const F *) = fn<decltype(comp)>("qfb_nucomp");
  void (*red)(Form *, const Form *,
              const F *) = fn<decltype(red)>("qfb_reduce");
  void (*mi)(void *, long, long) = fn<decltype(mi)>("fmpz_mat_init");
  void (*mc)(void *) = fn<decltype(mc)>("fmpz_mat_clear");
  F *(*entry)(void *, long, long) = fn<decltype(entry)>("fmpz_mat_entry");
  void (*li)(void *) = fn<decltype(li)>("fmpz_lll_context_init_default");
  void (*lll)(void *, void *, void *) = fn<decltype(lll)>("fmpz_lll");
  Z val(const F *z) {
    Z v;
    get(v.get_mpz_t(), z);
    return v;
  }
  void audit(const Form &f, Trace &t) {
    if (!t.split)
      for (auto x : {&f.a, &f.b, &f.c})
        t.seez(val(x));
  }
};
inline Flint &fl() {
  static Flint f;
  return f;
}
inline Out forms(Ring &r, int k, bool full) {
  auto &f = fl();
  Z D = -4 * k * r.nz;
  Flint::F disc = 0, L = 0;
  f.set(&disc, D.get_mpz_t());
  Z rt = sqroot(sqroot(-D / 4));
  f.set(&L, rt.get_mpz_t());
  Flint::Form base{}, cur{}, tmp{};
  f.initf(&base);
  f.initf(&cur);
  f.initf(&tmp);
  unsigned ell = 3, rr = 0;
  for (;; ell += 2) {
    if (!mpz_probab_prime_p(big(ell).get_mpz_t(), 10))
      continue;
    unsigned v = unsigned((ell - (r.N % ell) * k % ell) % ell);
    if (!v)
      continue;
    for (rr = 1; rr < ell && rr * rr % ell != v; ++rr) {
    }
    if (rr < ell)
      break;
  }
  if (2 * rr > ell)
    rr = ell - rr;
  Z c = (k * r.nz + rr * rr) / ell;
  f.set(&base.a, big(ell).get_mpz_t());
  f.set(&base.b, big(2 * rr).get_mpz_t());
  f.set(&base.c, c.get_mpz_t());
  W e = full ? r.N : (W(1) << ((bitlen(r.N) + 1) / 2));
  f.red(&cur, &base, &disc);
  f.audit(cur, r.tr);
  for (int j = int(bitlen(e)) - 2; j >= 0; --j) {
    f.dupl(&tmp, &cur, &disc, &L);
    f.audit(tmp, r.tr);
    f.red(&cur, &tmp, &disc);
    f.audit(cur, r.tr);
    if ((e >> j) & 1) {
      f.comp(&tmp, &cur, &base, &disc, &L);
      f.audit(tmp, r.tr);
      f.red(&cur, &tmp, &disc);
      f.audit(cur, r.tr);
    }
  }
  Z a = f.val(&cur.a), b = f.val(&cur.b), cc = f.val(&cur.c);
  Out out{a, b, cc, Z(4 * abs(b) <= a)};
  f.clearf(&tmp);
  f.clearf(&cur);
  f.clearf(&base);
  f.clear(&disc);
  f.clear(&L);
  return out;
}
inline Out lattice(Ring &r, bool gar, int dim) {
  auto &f = fl();
  alignas(16) unsigned char mat[256]{}, ctx[256]{};
  f.mi(mat, dim, dim);
  f.li(ctx);
  auto put = [&](int i, int j, const Z &v) {
    r.tr.seez(v);
    f.set(f.entry(mat, i, j), v.get_mpz_t());
  };
  Z root = sqroot(r.nz), x = sqroot(root);
  if (gar) {
    put(0, 0, r.nz);
    for (int i = 1; i < dim; ++i) {
      Z v = i == 1   ? big(r.power(2, r.N))
            : i == 2 ? big(r.power(3, r.N))
            : i == 3 ? root
                     : big(r.power(5, r.N));
      put(i, 0, v);
      put(i, i, x);
    }
  } else { // Coppersmith shift-polynomial basis N, x0+Xh, Xh(x0+Xh), ...
    put(0, 0, r.nz);
    Z scale = 1;
    for (int i = 1; i < dim; ++i) {
      put(i, i - 1, root * scale);
      scale *= x;
      put(i, i, scale);
    }
  }
  f.lll(mat, nullptr, ctx);
  Z norm = 0, l1 = 0, signs = 0, a = f.val(f.entry(mat, 0, 0));
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j) {
      Z v = f.val(f.entry(mat, i, j));
      r.tr.seez(v);
      if (!i) {
        norm += v * v;
        l1 += abs(v);
        if (v < 0)
          mpz_setbit(signs.get_mpz_t(), j);
      }
    }
  f.mc(mat);
  return {a, norm, l1, signs};
}
