#pragma once
#include "algebra.hpp"
struct Pari {
  using G = long *;
  void *lib = dlopen("/Volumes/git/factoring-connection/.venv/lib/python3.12/"
                     "site-packages/cypari2/.dylibs/libpari-gmp-tls.dylib",
                     RTLD_NOW | RTLD_LOCAL);
  template <class F> F fn(const char *s) {
    if (!lib)
      throw std::runtime_error(dlerror());
    auto p = dlsym(lib, s);
    if (!p)
      throw std::runtime_error(s);
    return reinterpret_cast<F>(p);
  }
  void (*init)(size_t, U) = fn<decltype(init)>("pari_init");
  U (*save)() = fn<decltype(save)>("get_avma");
  void (*restore)(U) = fn<decltype(restore)>("set_avma");
  G (*ui)(unsigned long) = fn<decltype(ui)>("utoi");
  G (*si)(long) = fn<decltype(si)>("stoi");
  G (*read)(const char *) = fn<decltype(read)>("strtoi");
  long (*itos)(G) = fn<decltype(itos)>("itos");
  G (*vec)(G, G, G) = fn<decltype(vec)>("mkvec3");
  void (*setup)(G, long, void *) = fn<decltype(setup)>("qfr_data_init");
  G (*five)(G, long) = fn<decltype(five)>("qfr_to_qfr5");
  G (*rho)(G, void *) = fn<decltype(rho)>("qfr5_rho");
  G (*comp)(G, G, void *) = fn<decltype(comp)>("qfr5_comp");
  G (*dist)(G, G, long) = fn<decltype(dist)>("qfr5_dist");
  int (*cmp)(G, G) = fn<decltype(cmp)>("gcmp");
  long (*expo)(G) = fn<decltype(expo)>("gexpo");
  int (*acmp)(G, G) = fn<decltype(acmp)>("abscmprr");
  G (*exp)(G, long) = fn<decltype(exp)>("gexp");
  G (*shift)(G, long) = fn<decltype(shift)>("shiftr");
  struct Dist {
    long e;
    G d;
  };
  Dist key(G v) {
    return {itos(reinterpret_cast<G>(v[4])) * (1L << 22),
            reinterpret_cast<G>(v[5])};
  }
  int compare(Dist a, Dist b) {
    if (a.e == b.e)
      return acmp(a.d, b.d);
    long h = a.e + expo(a.d) - b.e - expo(b.d);
    if (h)
      return h > 0 ? 1 : -1;
    return acmp(shift(a.d, a.e - b.e), b.d);
  }
  Pari() { init(16000000, 20000); }
  G distance(G v) {
    return dist(reinterpret_cast<G>(v[4]), reinterpret_cast<G>(v[5]), 5);
  }
  void audit(G v, Trace &t) {
    for (int i = 1; i <= 3 && !t.split; ++i)
      t.see(U(std::abs(itos(reinterpret_cast<G>(v[i])))));
  }
};
inline Out infrastructure(Ring &r, bool far) {
  static Pari f;
  auto av = f.save();
  alignas(16) unsigned char ctx[128]{};
  U s = isqrtw(r.N);
  auto D = f.read(Z(4 * r.nz).get_str().c_str());
  f.setup(D, 5, ctx);
  auto a = f.five(f.vec(f.ui(1), f.ui(2 * s), f.si(-long(r.N - W(s) * s))), 5);
  Pari::Dist target{0, f.exp(f.ui(2 * (far ? s : isqrtw(s))), 5)};
  auto cur = a;
  auto cd = f.key(a);
  std::vector<Pari::G> bases;
  auto one = f.rho(cur, ctx);
  f.audit(one, r.tr);
  for (int i = 0; i < 100; ++i) {
    auto d = f.key(one);
    if (f.compare(d, target) > 0)
      break;
    bases.push_back(one);
    auto next = f.comp(one, one, ctx);
    f.audit(next, r.tr);
    if (f.compare(f.key(next), d) <= 0)
      break;
    one = next;
  }
  for (auto it = bases.rbegin(); it != bases.rend(); ++it) {
    auto trial = f.comp(cur, *it, ctx);
    f.audit(trial, r.tr);
    auto td = f.key(trial);
    if (f.compare(td, cd) > 0 && f.compare(td, target) <= 0) {
      cur = trial;
      cd = td;
    }
  }
  bool bracket = false;
  for (int i = 0; i < 1024; ++i) {
    auto next = f.rho(cur, ctx);
    f.audit(next, r.tr);
    auto nd = f.key(next);
    if (f.compare(nd, target) > 0) {
      bracket = true;
      break;
    }
    cur = next;
  }
  if (!bracket)
    throw std::runtime_error("infrastructure distance bracket");
  Z u = f.itos(reinterpret_cast<Pari::G>(cur[1])),
    v = f.itos(reinterpret_cast<Pari::G>(cur[2])),
    w = f.itos(reinterpret_cast<Pari::G>(cur[3]));
  Out out{u, v, w, Z(2 * abs(u) < big(s))};
  f.restore(av);
  return out;
}
using Poly = std::array<W, 7>;
inline Poly polyprod(Ring &r, const Poly &a, const Poly &b, const Poly &rel,
                     int d) {
  std::array<W, 13> t{};
  for (int i = 0; i < d; ++i)
    for (int j = 0; j < d; ++j)
      if (a[i] && b[j])
        t[i + j] = r.add(t[i + j], r.mul(a[i], b[j]));
  for (int i = 2 * d - 2; i >= d; --i)
    for (int j = 0; j < d; ++j)
      if (rel[j])
        t[i - d + j] =
            r.add(t[i - d + j], rel[j] == 1 ? t[i] : r.mul(t[i], rel[j]));
  Poly v{};
  std::copy_n(t.begin(), d, v.begin());
  return v;
}
inline Poly polypow(Ring &r, Poly a, Poly rel, int d, W e) {
  Poly v{};
  v[0] = 1;
  for (int j = bitlen(e) - 1; j >= 0; --j) {
    v = polyprod(r, v, v, rel, d);
    if ((e >> j) & 1)
      v = polyprod(r, v, a, rel, d);
  }
  return v;
}
inline Out cf(Ring &r, int k, int len) {
  W D = r.N * k;
  U s = isqrtw(D), m = 0, q = 1, a = s;
  W old = 0, A = 1;
  unsigned even = 0, odd = 0;
  U mx = 0;
  int kr = 1;
  for (int j = 0; j < len; ++j) {
    W val = r.add(r.mul(a, A), old);
    old = A;
    A = val;
    m = q * a - m;
    q = (D - W(m) * m) / q;
    if (!q)
      throw std::runtime_error("square kN");
    r.tr.see(m);
    r.tr.see(q);
    r.tr.see(s + m);
    a = (s + m) / q;
    r.tr.see(a);
    mx = std::max(mx, a);
    U h = q;
    for (unsigned v : {2, 3, 5, 7})
      while (h % v == 0)
        h /= v;
    if (square64(h))
      ((j + 1) % 2 ? odd : even)++;
    if (square64(q)) {
      U t = isqrtw(q);
      r.add(A, t);
      r.sub(A, t);
    }
    static const int leg[11] = {0, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1};
    kr *= leg[q % 11];
  }
  return {big(q), big(mx), Z(even + 128 * odd), Z(kr)};
}
inline Out euclid(Ring &r) {
  W x = r.power(2, r.N), a = r.N, b = x, den0 = 0, den1 = 1, mid = 0;
  Z alt = 0;
  unsigned len = 0;
  while (b) {
    W q = a / b, y = a % b;
    r.tr.see(q);
    r.tr.see(y);
    W next = q * den1 + den0;
    r.tr.see(next);
    if (!mid && next >= isqrtw(r.N))
      mid = r.mul(next % r.N, y);
    if (len & 1)
      alt -= big(q);
    else
      alt += big(q);
    ++len;
    a = b;
    b = y;
    den0 = den1;
    den1 = next;
  }
  return {alt, big(mid), big(den0), Z(len)};
}
inline std::pair<W, W> ladder(Ring &r, W x, W scalar, W a24) {
  W X0 = 1, Z0 = 0, X1 = x, Z1 = 1;
  for (int j = bitlen(scalar) - 1; j >= 0; --j) {
    bool bit = (scalar >> j) & 1;
    if (bit) {
      std::swap(X0, X1);
      std::swap(Z0, Z1);
    }
    W A = r.add(X0, Z0), B = r.sub(X0, Z0), aa = r.mul(A, A), bb = r.mul(B, B),
      e = r.sub(aa, bb);
    W da = r.mul(r.sub(X1, Z1), A), cb = r.mul(r.add(X1, Z1), B);
    W l = r.add(da, cb), h = r.sub(da, cb);
    X1 = r.mul(l, l);
    Z1 = r.mul(x, r.mul(h, h));
    X0 = r.mul(aa, bb);
    Z0 = r.mul(e, r.add(bb, r.mul(a24, e)));
    if (bit) {
      std::swap(X0, X1);
      std::swap(Z0, Z1);
    }
  }
  return {X0, Z0};
}
inline Out ec(Ring &r, int mode) {
  W A = mode == 4 ? 4 * ((r.N >> 7) % 17) + 6 : mode == 5 ? 0 : 6;
  W a24 = r.mul(A + 2, r.inv(4)), x = mode == 1 ? 3 : 2;
  auto curvecheck = [&](W xx) {
    W xx2 = r.mul(xx, xx);
    W b = r.add(r.add(r.mul(xx2, xx), r.mul(A, xx2)), xx);
    r.tr.see(b);
    r.sub(r.mul(A, A), 4);
  };
  curvecheck(x);
  W scalar = mode == 1   ? r.N + 1
             : mode == 2 ? r.N + 1 - 2 * isqrtw(r.N)
             : mode == 3 ? M(32)
                         : r.N;
  auto [X, Zc] = ladder(r, x, scalar, a24);
  if (mode == 6) {
    x = r.add(r.mul(X, r.inv(Zc)), 1);
    curvecheck(x);
    std::tie(X, Zc) = ladder(r, x, r.N, a24);
  }
  W affine = r.mul(X, r.inv(Zc));
  return output(X, Zc, affine, r.sub(X, Zc));
}
struct Jacobian {
  W x, y, z;
};
inline Out ec37(Ring &r, W e) {
  W half = r.inv(2);
  Jacobian g{0, half, 1}, p = g;
  auto dbl = [&](Jacobian P) {
    auto [x, y, z] = P;
    W xx = r.mul(x, x), yy = r.mul(y, y), yyyy = r.mul(yy, yy),
      zz = r.mul(z, z);
    W t = r.add(x, yy);
    W S = r.mul(2, r.sub(r.sub(r.mul(t, t), xx), yyyy));
    W slope = r.sub(r.mul(3, xx), r.mul(zz, zz));
    W X = r.sub(r.mul(slope, slope), r.mul(2, S));
    W Y = r.sub(r.mul(slope, r.sub(S, X)), r.mul(8, yyyy));
    return Jacobian{X, Y, r.mul(2, r.mul(y, z))};
  };
  for (int j = int(bitlen(e)) - 2; j >= 0; --j) {
    p = dbl(p);
    if (!((e >> j) & 1))
      continue;
    if (!p.z) {
      p = g;
      continue;
    }
    W zz = r.mul(p.z, p.z), h = r.sub(0, p.x),
      diff = r.sub(r.mul(half, r.mul(p.z, zz)), p.y);
    if (h == 0) {
      p = diff == 0 ? dbl(p) : Jacobian{0, 1, 0};
      continue;
    }
    W hh = r.mul(h, h), I = r.mul(4, hh), J = r.mul(h, I), R = r.mul(2, diff),
      V = r.mul(p.x, I);
    W X = r.sub(r.sub(r.mul(R, R), J), r.mul(2, V)),
      Y = r.sub(r.mul(R, r.sub(V, X)), r.mul(2, r.mul(p.y, J)));
    W Zc = r.mul(2, r.mul(p.z, h));
    p = {X, Y, Zc};
  }
  W zi = r.inv(p.z);
  W affine = r.mul(p.x, r.mul(zi, zi));
  return output(p.x, p.y, p.z, affine);
}
