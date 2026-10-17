"""Track A (a): binary section transforms E/O on hash-consed arithmetic circuits.

A circuit is a DAG over R_N = Z/NZ in one variable.  Nodes are hash-consed
(structurally identical nodes are shared) and trivial identities are folded:
constants are reduced mod N, 0+a=a, 0*a=0, 1*a=a, const op const = const.

Section transform (one level).  For a circuit in variable Z, every node v is
mapped to a pair (e(v), o(v)) of nodes of a NEW circuit in variable Y = Z^2,
with v = e(v)(Z^2) + Z o(v)(Z^2):
    const c   -> (c, 0)
    Z         -> (0, 1)
    a + b     -> (e(a)+e(b), o(a)+o(b))
    a * b     -> (e(a)e(b) + Y o(a)o(b),  e(a)o(b) + o(a)e(b))
Only one section (bit 0 -> E, bit 1 -> O) is kept per level, and the new
circuit is pruned to the nodes reachable from that output.  Iterating k times
gives the circuit of S_{r,2^k}(Y) where r has binary digits = chosen bits
(least significant first).

The circuit for (1+X)^N (or F_N = (1+X)^N - 1 - X^N) is built by left-to-right
repeated squaring: O(log N) nodes.

Operational: uses only N and the chosen bits.  `run_path` reports per level the
number of unique reachable nodes, number of mul/add nodes, and an upper bound
on the degree in the current variable.
"""
from __future__ import annotations

import sys
from typing import Dict, List, Optional, Sequence, Tuple

CONST, VAR, ADD, MUL = 0, 1, 2, 3


class Circuit:
    """Hash-consed arithmetic circuit over Z/NZ in one variable."""

    def __init__(self, N: int):
        self.N = N
        self.nodes: List[Tuple] = []          # node id -> (op, a, b) / (CONST, c, None)
        self.table: Dict[Tuple, int] = {}
        self.deg: List[int] = []
        self.zero = self.const(0)
        self.one = self.const(1)
        self.var = self._mk((VAR, 0, None), 1)

    def _mk(self, key: Tuple, deg: int) -> int:
        nid = self.table.get(key)
        if nid is None:
            nid = len(self.nodes)
            self.nodes.append(key)
            self.table[key] = nid
            self.deg.append(deg)
        return nid

    def const(self, c: int) -> int:
        return self._mk((CONST, c % self.N, None), 0)

    def is_const(self, a: int) -> Optional[int]:
        node = self.nodes[a]
        return node[1] if node[0] == CONST else None

    def add(self, a: int, b: int) -> int:
        ca, cb = self.is_const(a), self.is_const(b)
        if ca == 0:
            return b
        if cb == 0:
            return a
        if ca is not None and cb is not None:
            return self.const(ca + cb)
        if a > b:
            a, b = b, a
        return self._mk((ADD, a, b), max(self.deg[a], self.deg[b]))

    def mul(self, a: int, b: int) -> int:
        ca, cb = self.is_const(a), self.is_const(b)
        if ca == 0 or cb == 0:
            return self.zero
        if ca == 1:
            return b
        if cb == 1:
            return a
        if ca is not None and cb is not None:
            return self.const(ca * cb)
        if a > b:
            a, b = b, a
        return self._mk((MUL, a, b), self.deg[a] + self.deg[b])

    def sub(self, a: int, b: int) -> int:
        return self.add(a, self.mul(self.const(-1), b))

    # ------------------------------------------------------------------
    def reachable(self, out: int) -> List[int]:
        seen = set()
        stack = [out]
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            op, a, b = self.nodes[v]
            if op in (ADD, MUL):
                stack.append(a)
                stack.append(b)
        return sorted(seen)

    def stats(self, out: int) -> Dict[str, int]:
        r = self.reachable(out)
        ops = [self.nodes[v][0] for v in r]
        return {"nodes": len(r), "mul": ops.count(MUL), "add": ops.count(ADD),
                "const": ops.count(CONST), "deg_bound": self.deg[out]}

    def evaluate(self, out: int, x: int) -> int:
        """Evaluate at a scalar (topological order = id order)."""
        N = self.N
        val: Dict[int, int] = {}
        for v in self.reachable(out):
            op, a, b = self.nodes[v]
            if op == CONST:
                val[v] = a
            elif op == VAR:
                val[v] = x % N
            elif op == ADD:
                val[v] = (val[a] + val[b]) % N
            else:
                val[v] = val[a] * val[b] % N
        return val[out]

    def expand(self, out: int) -> List[int]:
        """Dense coefficient list (small cases only, for testing)."""
        N = self.N
        val: Dict[int, List[int]] = {}
        for v in self.reachable(out):
            op, a, b = self.nodes[v]
            if op == CONST:
                val[v] = [a]
            elif op == VAR:
                val[v] = [0, 1]
            elif op == ADD:
                x, y = val[a], val[b]
                n = max(len(x), len(y))
                val[v] = [((x[i] if i < len(x) else 0) + (y[i] if i < len(y) else 0)) % N
                          for i in range(n)]
            else:
                x, y = val[a], val[b]
                z = [0] * (len(x) + len(y) - 1)
                for i, xi in enumerate(x):
                    if xi:
                        for j, yj in enumerate(y):
                            z[i + j] += xi * yj
                val[v] = [t % N for t in z]
        res = val[out]
        while len(res) > 1 and res[-1] == 0:
            res.pop()
        return res


def build_power(N: int, exponent: Optional[int] = None, subtract_endpoints: bool = False
                ) -> Tuple[Circuit, int]:
    """Circuit for (1+X)^e over Z/NZ by left-to-right repeated squaring
    (e = N by default).  With subtract_endpoints, builds F_N = (1+X)^N - 1 - X^N
    (X^N also by repeated squaring, sharing nothing with the first chain)."""
    e = N if exponent is None else exponent
    C = Circuit(N)
    lin = C.add(C.one, C.var)

    def power(base: int, k: int) -> int:
        r = None
        for bit in bin(k)[2:]:
            if r is not None:
                r = C.mul(r, r)
            if bit == "1":
                r = base if r is None else C.mul(r, base)
        return r

    out = power(lin, e)
    if subtract_endpoints:
        out = C.sub(C.sub(out, C.one), power(C.var, e))
    return C, out


def section_once(C: Circuit, out: int, bit: int, karatsuba: bool = False) -> Tuple[Circuit, int]:
    """Apply E (bit=0) or O (bit=1) to the circuit rooted at `out`.
    Returns a fresh circuit (variable Y = X^2) and its output node.
    karatsuba=True computes O(ab) = (e_a+o_a)(e_b+o_b) - e_a e_b - o_a o_b
    (3 products per mul gate instead of 4)."""
    D = Circuit(C.N)
    Y = D.var
    pair: Dict[int, Tuple[int, int]] = {}
    for v in C.reachable(out):  # ids are topologically ordered
        op, a, b = C.nodes[v]
        if op == CONST:
            pair[v] = (D.const(a), D.zero)
        elif op == VAR:
            pair[v] = (D.zero, D.one)
        elif op == ADD:
            ea, oa = pair[a]
            eb, ob = pair[b]
            pair[v] = (D.add(ea, eb), D.add(oa, ob))
        else:
            ea, oa = pair[a]
            eb, ob = pair[b]
            ee, oo = D.mul(ea, eb), D.mul(oa, ob)
            e = D.add(ee, D.mul(Y, oo))
            if karatsuba and a != b and ea != D.zero and oa != D.zero:
                o = D.sub(D.sub(D.mul(D.add(ea, oa), D.add(eb, ob)), ee), oo)
            elif a == b:
                o = D.mul(D.const(2), D.mul(ea, oa))
            else:
                o = D.add(D.mul(ea, ob), D.mul(oa, eb))
            pair[v] = (e, o)
    new_out = pair[out][bit]
    # compact: rebuild only reachable part so node ids/counts reflect the section
    return compact(D, new_out)


def compact(C: Circuit, out: int) -> Tuple[Circuit, int]:
    D = Circuit(C.N)
    m: Dict[int, int] = {}
    for v in C.reachable(out):
        op, a, b = C.nodes[v]
        if op == CONST:
            m[v] = D.const(a)
        elif op == VAR:
            m[v] = D.var
        elif op == ADD:
            m[v] = D.add(m[a], m[b])
        else:
            m[v] = D.mul(m[a], m[b])
    return D, m[out]


def run_path(N: int, bits: Sequence[int], subtract_endpoints: bool = True,
             node_cap: int = 3_000_000, karatsuba: bool = False) -> List[Dict[str, int]]:
    """Section the circuit of F_N (or (1+X)^N) along `bits` (LSB first).
    Returns per-level stats; stops early if the node count exceeds node_cap."""
    C, out = build_power(N, subtract_endpoints=subtract_endpoints)
    rows = [dict(level=0, **C.stats(out))]
    for k, bit in enumerate(bits, start=1):
        C, out = section_once(C, out, bit, karatsuba)
        st = C.stats(out)
        rows.append(dict(level=k, bit=bit, **st))
        if st["nodes"] > node_cap:
            break
    return rows


if __name__ == "__main__":  # pragma: no cover
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 143
    for row in run_path(N, [0] * 8):
        print(row)
