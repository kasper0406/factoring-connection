"""detGRH diagnostics (HARNESS / DIAGNOSTIC: uses the factorisation of N; never imported by core.py).

- instance generation by congruence class (blum / mixed / both1 / tri), SHA-seeded;
- Lemma MR exhaustive check (E1a);
- the full relation lattice Lambda = ker(ev: Z^FB -> (Z/N)^*) from discrete logarithms modulo each p | N,
  the 2-rank c of Lambda/Lambda_c, and dim W' = dim (Lambda cap 2Z^n)/2Lambda (E1c).
"""
from __future__ import annotations

import random
from itertools import product as iproduct
from math import gcd, prod
from typing import Dict, List, Optional, Sequence, Tuple

import flint
import gmpy2
import sympy

from .core import FactorBase, bits_of


# ----------------------------------------------------------------------------- instances

def _rng(cls: str, bits: int, seed: int) -> random.Random:
    return random.Random(f"detGRH|{cls}|bits={bits}|{seed}")


def _prime(rng: random.Random, bits: int, mod4: Optional[int]) -> int:
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if mod4 is not None and x % 4 != mod4:
            continue
        if gmpy2.is_prime(x, 50):
            return x


def make_instance(cls: str, bits: int, seed: int, min_factor: int = 0) -> Tuple[int, List[int]]:
    """N with exactly `bits` bits and distinct prime factors > min_factor, of the given class."""
    rng = _rng(cls, bits, seed)
    while True:
        if cls == "tri":
            b1 = -(-bits // 3)          # ceil(bits/3) for all three, keep only products of exactly `bits` bits
            ps = [_prime(rng, b1, None), _prime(rng, b1, None), _prime(rng, b1, None)]
        else:
            m1, m2 = {"blum": (3, 3), "mixed": (1, 3), "both1": (1, 1)}[cls]
            ps = [_prime(rng, (bits + 1) // 2, m1), _prime(rng, bits // 2, m2)]
        if len(set(ps)) < len(ps) or min(ps) <= min_factor:
            continue
        N = prod(ps)
        if N.bit_length() == bits:
            return N, sorted(ps)


# ----------------------------------------------------------------------------- Lemma MR (E1a)

def _cyclic_components(N: int) -> List[Tuple[int, int, int]]:
    """[(modulus p^e, generator lifted by CRT to Z/N (= g mod p^e, 1 elsewhere), order)] for odd N."""
    out = []
    for p, e in sympy.factorint(N).items():
        q = p ** e
        g = int(sympy.primitive_root(q))
        rest = N // q
        # CRT: x = g mod q, x = 1 mod rest
        x = (g * rest * pow(rest, -1, q) + 1 * q * pow(q, -1, rest)) % N if rest > 1 else g % N
        out.append((q, x, (p - 1) * p ** (e - 1)))
    return out


def mu2(N: int) -> List[int]:
    return [x for x in range(1, N) if x * x % N == 1]


def count_sections(N: int, brute: bool = True) -> int:
    """Number of homomorphisms s: (Z/N)^{*2} -> (Z/N)^* with s(x)^2 = x.  Squares decompose as the direct
    product of the cyclic groups <g_i^2> (order n_i/2); s is determined by s(g_i^2), which must be a root of
    g_i^2 (i.e. g_i * eps, eps in mu_2) of order dividing n_i/2.  If brute, every candidate section is also
    expanded on all of Q and s(x)^2 = x and multiplicativity are checked element by element."""
    comps = _cyclic_components(N)
    m2 = mu2(N)
    choices = []
    for (_, g, n) in comps:
        c = [g * e % N for e in m2 if pow(g * e % N, n // 2, N) == 1]
        choices.append(c)
    total = 0
    for imgs in iproduct(*choices):
        if brute:
            # expand s on Q = prod <g_i^2>
            qs = [(g * g % N, n // 2) for (_, g, n) in comps]
            table: Dict[int, int] = {}
            ok = True
            for exps in iproduct(*[range(o) for (_, o) in qs]):
                x, y = 1, 1
                for (q, _), a, im in zip(qs, exps, imgs):
                    x = x * pow(q, a, N) % N
                    y = y * pow(im, a, N) % N
                if x in table and table[x] != y:
                    ok = False
                    break
                table[x] = y
                if y * y % N != x:
                    ok = False
                    break
            if not ok:
                continue
        total += 1
    return total


def principal_root_multiplicative(N: int) -> bool:
    """For N with every p = 3 mod 4: the root of each square that is itself a square; check multiplicativity
    over all pairs of squares (brute force)."""
    units = [x for x in range(1, N) if gcd(x, N) == 1]
    sq = sorted({x * x % N for x in units})
    sqset = set(sq)
    pr: Dict[int, int] = {}
    for z in units:
        s = z * z % N
        if z in sqset:
            if s in pr and pr[s] != z:
                return False       # two square roots that are squares: not unique
            pr[s] = z
    if len(pr) != len(sq):
        return False
    for a in sq:
        for b in sq:
            if pr[a * b % N] != pr[a] * pr[b] % N:
                return False
    return True


def all_3mod4(N: int) -> bool:
    return all(p % 4 == 3 for p in sympy.factorint(N))


# ----------------------------------------------------------------------------- full lattice (E1c)

def dlog_table(p: int) -> Tuple[int, Dict[int, int]]:
    g = int(sympy.primitive_root(p))
    tab = {}
    x = 1
    for k in range(p - 1):
        tab[x] = k
        x = x * g % p
    return g, tab


def full_lattice(fb: FactorBase, factors: Sequence[int]) -> List[List[int]]:
    """Basis (rows) of Lambda = {e in Z^n : prod q^{e_q} = 1 mod N} for squarefree N = prod factors.
    Built as the kernel of Z^n -> prod Z/(p-1), e -> (sum_q e_q dlog_p(q)) via HNF."""
    n = fb.n
    cols = []    # one per prime p: dlog of each FB element mod p
    mods = []
    for p in factors:
        _, tab = dlog_table(p)
        col = [(p - 1) // 2] + [tab[q % p] for q in fb.primes]   # dlog(-1) = (p-1)/2
        cols.append(col)
        mods.append(p - 1)
    r = len(factors)
    rows = []
    for j in range(n):
        rows.append([cols[i][j] for i in range(r)] + [1 if k == j else 0 for k in range(n)])
    for i in range(r):
        rows.append([mods[i] if k == i else 0 for k in range(r)] + [0] * n)
    H = flint.fmpz_mat(rows).hnf().tolist()
    basis = []
    for row in H:
        row = [int(x) for x in row]
        if any(row) and all(x == 0 for x in row[:r]):
            basis.append(row[r:])
    assert len(basis) == n, (len(basis), n)
    return basis


def two_rank_quotient(Lam: List[List[int]], sub: List[List[int]]) -> int:
    """c = dim_F2 Lambda/(Lambda_c + 2 Lambda): coordinates of the generators of Lambda_c in the basis of
    Lambda (exact rational solve; integral since Lambda_c is inside Lambda), reduced mod 2."""
    n = len(Lam)
    Hs = [[int(x) for x in r] for r in flint.fmpz_mat(sub).hnf().tolist()]
    Hs = [r for r in Hs if any(r)]
    assert len(Hs) == n, "Lambda_c must have full rank"
    At = flint.fmpq_mat(flint.fmpz_mat(Lam).transpose())
    sol = At.solve(flint.fmpq_mat(flint.fmpz_mat(Hs).transpose()))   # columns = coordinates
    coords = []
    for j in range(n):
        ci = []
        for i in range(n):
            c = sol[i, j]
            assert c.q == 1, "Lambda_c not inside Lambda"
            ci.append(int(c.p))
        coords.append(ci)
    # F2 rank of the coordinate rows
    piv: Dict[int, int] = {}
    pivmask = 0
    rk = 0
    for c in coords:
        b = bits_of(c)
        t = b & pivmask
        while t:
            j = (t & -t).bit_length() - 1
            b ^= piv[j]
            t = b & pivmask
        if b:
            j = (b & -b).bit_length() - 1
            for jj in list(piv):
                if (piv[jj] >> j) & 1:
                    piv[jj] ^= b
            piv[j] = b
            pivmask |= 1 << j
            rk += 1
    return n - rk


def dim_W(Lam: List[List[int]]) -> int:
    """dim (Lambda cap 2Z^n)/2Lambda = n - rank_F2(Lambda basis mod 2)  (= 2-rank of G_B)."""
    n = len(Lam)
    piv: Dict[int, int] = {}
    pivmask = 0
    rk = 0
    for v in Lam:
        b = bits_of(v)
        t = b & pivmask
        while t:
            j = (t & -t).bit_length() - 1
            b ^= piv[j]
            t = b & pivmask
        if b:
            j = (b & -b).bit_length() - 1
            for jj in list(piv):
                if (piv[jj] >> j) & 1:
                    piv[jj] ^= b
            piv[j] = b
            pivmask |= 1 << j
            rk += 1
    return n - rk


def lattice_det(Lam: List[List[int]]) -> int:
    return abs(int(flint.fmpz_mat(Lam).det()))


def small_primes_generate(N: int, factors: Sequence[int], bound: int) -> bool:
    """Do the primes < bound (and -1) generate (Z/N)^*?  Index of their relation lattice = phi(N)."""
    ps = [q for q in sympy.primerange(2, bound)]
    fb = FactorBase(N=N, B=bound, primes=ps)
    fb.index = {p: i + 1 for i, p in enumerate(ps)}
    Lam = full_lattice(fb, factors)
    phi = prod(p - 1 for p in factors)
    return lattice_det(Lam) == phi
