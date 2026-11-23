#include "variants.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
// The sampler's independent ChaCha streams are reused for orientation and
// planted bits.
#define main entropy_sampler_entry
#include "sampler.cpp"
#undef main
struct Ref {
  int v, slot, bit;
};
struct Lock {
  int op;
  std::vector<Ref> refs;
};
int main(int argc, char **argv) {
  if (argc != 7)
    return 2;
  int n = std::stoi(argv[1]);
  std::string key = argv[2];
  std::ifstream cfg(argv[3]), in(argv[4], std::ios::binary);
  std::ofstream out(argv[5], std::ios::binary);
  std::string audit_mode = argv[6];
  std::vector<Lock> locks;
  std::set<int> vs;
  bool needed[64][4]{};
  int op, k;
  while (cfg >> op >> k) {
    Lock l{op, {}};
    for (int i = 0; i < k; ++i) {
      Ref v;
      cfg >> v.v >> v.slot >> v.bit;
      l.refs.push_back(v);
      vs.insert(v.v);
      needed[v.v][v.slot] = true;
    }
    locks.push_back(l);
  }
  if (locks.size() > 32 || !in || !out)
    throw std::runtime_error("lock input");
  Random orientation("ent-c-R2-v1/confirmation/orientation/" + key),
      one("ent-c-R2-v1/confirmation/1n/" + key),
      mid("ent-c-R2-v1/confirmation/16n2/" + key),
      weak("ent-c-R2-v1/confirmation/1n2/" + key);
  unsigned index[3840]{}, j = 0;
  for (unsigned r = 0; r < 3840; ++r)
    if (std::gcd(r, 3840u) == 1)
      index[r] = j++;
  U pq[2];
  W lower = W(9) << (n - 4);
  uint16_t bits[64][4]{};
  bool hit[64]{};
  std::vector<uint32_t> buf;
  buf.reserve(65536);
  while (in.read((char *)pq, 16)) {
    U p = pq[0], q = pq[1];
    W N = W(p) * q;
    Z nz = big(N);
    uint32_t xb = 0, flags = 0;
    for (int v : vs) {
      Trace tr{p, q};
      Ring r{N, nz, tr};
      auto z = variant(v, r);
      for (int slot = 0; slot < 4; ++slot) {
        if (needed[v][slot]) {
          auto b = projection(z[slot], r);
          bits[v][slot] = r.singular ? 0 : b;
        } else
          tr.seez(z[slot]);
      }
      hit[v] = tr.split;
    }
    for (unsigned i = 0; i < locks.size(); ++i) {
      auto &l = locks[i];
      bool sum[3]{};
      bool flag = false;
      for (unsigned j = 0; j < l.refs.size(); ++j) {
        auto e = l.refs[j];
        sum[j] = (bits[e.v][e.slot] >> e.bit) & 1;
        flag |= hit[e.v];
      }
      bool bit = l.op == 0   ? sum[0]
                 : l.op == 1 ? sum[0] ^ sum[1]
                             : (sum[0] + sum[1] + sum[2] >= 2);
      xb |= uint32_t(bit) << i;
      flags |= uint32_t(flag) << i;
    }
    unsigned coin = orientation.next() & 1;
    unsigned success1 = one.below(2 * n) < U(n + 2),
             success16 = mid.below(2 * n * n) < U(n * n + 32),
             success2 = weak.below(2 * n * n) < U(n * n + 2);
    unsigned sen = __builtin_parityll(U(N) & 0xd5a65bb99d13a6adULL) ^
                   __builtin_parityll(U(N >> 64) & 0xe1990aafa5babc65ULL);
    unsigned c = coin | (((p >> (n / 4)) & 1) << 1) | (success1 << 2) |
                 (success16 << 3) | (success2 << 4) | (sen << 5) |
                 ((((coin ? q : p) >> 1) & 1) << 6) | (((N % 4) == 3) << 7);
    unsigned h = index[(p % 3840) * (q % 3840) % 3840] * 1024 +
                 unsigned(((N - lower) >> (n - 14)) / 7);
    buf.insert(buf.end(), {h, xb, flags, c});
    if (buf.size() == 65536) {
      out.write((char *)buf.data(), buf.size() * 4);
      buf.clear();
    }
  }
  out.write((char *)buf.data(), buf.size() * 4);
  return in.eof() && out ? 0 : 1;
}
