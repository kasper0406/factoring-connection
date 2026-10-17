"""OPERATIONAL decoders for derived partial-information formats.

Every function receives N, public parameters and the format object only.  Each
reduces the format to a *coset list*  {(M_i, r_i)}: "p ≡ r_i (mod M_i) for some i",
and calls the affine decoder (coppersmith.decode_affine) on each element until a
divisor is verified.  Cost of a format = Σ_i decode(M_i)  (list tolerance §2).

Auxiliary moduli l are public primes chosen by the decoder/format, never derived from
factors.  Square roots modulo a *public prime* l (Tonelli–Shanks) do not involve N.
"""
from __future__ import annotations

import itertools
import math
import random
from math import gcd, isqrt, prod
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from src.common.metrics import Metrics, NullMetrics
from src.trackE.coppersmith import AffineHint, balanced_promise, decode_affine


# ---------------------------------------------------------------------------
# helpers on public primes
# ---------------------------------------------------------------------------

def small_primes(B: int) -> List[int]:
    sieve = bytearray([1]) * (B + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(B) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(B + 1) if sieve[i]]


def legendre(a: int, l: int) -> int:
    a %= l
    if a == 0:
        return 0
    return 1 if pow(a, (l - 1) // 2, l) == 1 else -1


def sqrt_mod_prime(a: int, l: int) -> Optional[int]:
    """Tonelli–Shanks for a public odd prime l; None if a is a non-residue."""
    a %= l
    if a == 0:
        return 0
    if pow(a, (l - 1) // 2, l) != 1:
        return None
    if l % 4 == 3:
        return pow(a, (l + 1) // 4, l)
    Q, S = l - 1, 0
    while Q % 2 == 0:
        Q //= 2; S += 1
    z = 2
    while pow(z, (l - 1) // 2, l) != l - 1:
        z += 1
    Mx, c, t, R = S, pow(z, Q, l), pow(a, Q, l), pow(a, (Q + 1) // 2, l)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % l; i += 1
        b = pow(c, 1 << (Mx - i - 1), l)
        Mx, c, t, R = i, b * b % l, t * b * b % l, R * b % l
    return R


def crt(residues: Sequence[int], moduli: Sequence[int]) -> Tuple[int, int]:
    r, M = 0, 1
    for a, m in zip(residues, moduli):
        # moduli pairwise coprime (public primes)
        t = ((a - r) * pow(M, -1, m)) % m
        r, M = r + M * t, M * m
    return r % M, M


# ---------------------------------------------------------------------------
# root cosets of  x^2 - s x + N  modulo l^k
# ---------------------------------------------------------------------------

def quad_root_cosets(s: int, N: int, l: int, k: int, max_nodes: int = 10 ** 6) -> List[Tuple[int, int]]:
    """All x mod l^k with x^2 - s x + N ≡ 0, as disjoint cosets (a, e): x ≡ a mod l^e.

    p-adic root tree: node (a, j) has f(a) ≡ 0 mod l^j; it is emitted as a whole
    coset once f(a + l^j z) ≡ 0 mod l^k identically in z.  Branching ≤ l per level
    (use for small public l only)."""
    lk = l ** k
    out: List[Tuple[int, int]] = []
    stack = [(0, 0)]
    nodes = 0
    while stack:
        a, j = stack.pop()
        nodes += 1
        if nodes > max_nodes:
            raise RuntimeError("root tree too large")
        lj = l ** j
        fa = (a * a - s * a + N) % lk
        dfa = (2 * a - s) * lj % lk
        if fa == 0 and dfa == 0 and (lj * lj) % lk == 0:
            out.append((a, j))
            continue
        if j == k:
            continue
        lj1 = lj * l
        for z in range(l):
            a2 = a + lj * z
            if (a2 * a2 - s * a2 + N) % lj1 == 0:
                stack.append((a2, j + 1))
    return sorted(out, key=lambda t: (-t[1], t[0]))


def quad_roots_mod_prime(s: int, N: int, l: int) -> List[int]:
    """Roots of x^2 - s x + N mod an odd public prime l (Tonelli–Shanks on s^2 - 4N)."""
    D = (s * s - 4 * N) % l
    t = sqrt_mod_prime(D, l)
    if t is None:
        return []
    inv2 = (l + 1) // 2
    return sorted({(s + t) * inv2 % l, (s - t) * inv2 % l})


# ---------------------------------------------------------------------------
# generic coset-list decoder
# ---------------------------------------------------------------------------

def decode_coset_list(N: int, cosets: Iterable[Tuple[int, int]], d: Optional[int] = None,
                      metrics: Optional[Metrics] = None, max_decodes: Optional[int] = None,
                      promise: Optional[Tuple[int, int]] = None) -> Optional[int]:
    """Try each (M, r) in order; returns the first verified proper divisor."""
    mt = metrics if metrics is not None else NullMetrics()
    L, U = promise or balanced_promise(N)
    for i, (M, r) in enumerate(cosets):
        if max_decodes is not None and i >= max_decodes:
            break
        sub = Metrics()
        g = decode_affine(N, AffineHint(M, r % M, L, U), d=d, metrics=sub)
        mt.inc("decodes")
        mt.merge(sub, prefix="dec_")
        mt.fields["decode_time_s"] = mt.fields.get("decode_time_s", 0.0) + sub.elapsed()
        if g is not None and 1 < g < N and N % g == 0:
            mt.set("success_index", i)
            mt.set("success_M_bits", M.bit_length())
            return g
    return None


# ---------------------------------------------------------------------------
# (i)  s = p + q mod 2^k   (equivalently phi(N) mod 2^k)
# ---------------------------------------------------------------------------

def sum_mod_2k_cosets(N: int, s: int, k: int) -> List[Tuple[int, int]]:
    """(M, r) list for format (i): all root cosets of x^2 - s x + N mod 2^k."""
    return [(1 << e, a) for a, e in quad_root_cosets(s % (1 << k), N, 2, k)]


def decode_sum_mod_2k(N: int, s: int, k: int, d: Optional[int] = None,
                      metrics: Optional[Metrics] = None) -> Optional[int]:
    mt = metrics if metrics is not None else NullMetrics()
    cos = sum_mod_2k_cosets(N, s, k)
    mt.set("n_cosets", len(cos))
    mt.set("coset_bits", [M.bit_length() - 1 for M, _ in cos])
    return decode_coset_list(N, cos, d=d, metrics=mt)


def decode_phi_mod_2k(N: int, phi_mod: int, k: int, **kw) -> Optional[int]:
    """phi(N) = N + 1 - s, so this is format (i) verbatim."""
    return decode_sum_mod_2k(N, (N + 1 - phi_mod) % (1 << k), k, **kw)


# ---------------------------------------------------------------------------
# (ii)  s mod M,  M = prod of distinct odd public primes, optional per-l labels
# ---------------------------------------------------------------------------

def sum_mod_M_candidates(N: int, s_res: Dict[int, int],
                         labels: Optional[Dict[int, int]] = None) -> Tuple[List[List[int]], List[int]]:
    """Per-prime root lists for s mod l; `labels[l]` = Legendre symbol of the target
    factor mod l.  A label selects a root iff the two roots have different symbols,
    which happens exactly when (N/l) = -1 (proof in notes/trackE_analysis.md §3)."""
    roots, mods = [], []
    for l, sl in s_res.items():
        rts = quad_roots_mod_prime(sl, N, l)
        if labels is not None and l in labels and len(rts) == 2:
            sel = [x for x in rts if legendre(x, l) == labels[l]]
            if len(sel) == 1:
                rts = sel
        roots.append(rts)
        mods.append(l)
    return roots, mods


def decode_sum_mod_M(N: int, s_res: Dict[int, int], labels: Optional[Dict[int, int]] = None,
                     d: Optional[int] = None, metrics: Optional[Metrics] = None,
                     max_decodes: Optional[int] = None, order_seed: Optional[int] = None) -> Optional[int]:
    """Enumerate all 2^{omega'} labelings (CRT), decode each (in seeded random order)."""
    mt = metrics if metrics is not None else NullMetrics()
    roots, mods = sum_mod_M_candidates(N, s_res, labels)
    T = prod(len(r) for r in roots)
    mt.set("list_size", T)
    mt.set("omega", len(mods))
    mt.set("omega_ambiguous", sum(1 for r in roots if len(r) == 2))
    combos = list(itertools.product(*roots)) if T <= 2 ** 20 else None
    if combos is None:
        mt.set("status", "list_too_large")
        return None
    if order_seed is not None:
        random.Random(order_seed).shuffle(combos)
    Mtot = prod(mods)
    cos = ((Mtot, crt(c, mods)[0]) for c in combos)
    return decode_coset_list(N, cos, d=d, metrics=mt, max_decodes=max_decodes)


# ---------------------------------------------------------------------------
# (iii) p mod l for many small primes (consistent labels), possibly with errors
# ---------------------------------------------------------------------------

def decode_residues(N: int, primes: Sequence[int], residues: Sequence[int], d: Optional[int] = None,
                    metrics: Optional[Metrics] = None) -> Optional[int]:
    r, M = crt(residues, primes)
    mt = metrics if metrics is not None else NullMetrics()
    mt.set("M_bits", M.bit_length())
    return decode_coset_list(N, [(M, r)], d=d, metrics=mt)


def decode_residues_noisy(N: int, primes: Sequence[int], residues: Sequence[int], need_bits: int,
                          trials: int, seed: int = 0, d: Optional[int] = None,
                          metrics: Optional[Metrics] = None) -> Optional[int]:
    """Random-subset decoder: repeatedly pick a random subset of the residues whose
    modulus has >= need_bits bits and decode.  Succeeds once a subset is error-free."""
    mt = metrics if metrics is not None else NullMetrics()
    rng = random.Random(seed)
    idx = list(range(len(primes)))
    cos = []
    for _ in range(trials):
        rng.shuffle(idx)
        sub, Mb = [], 1
        for i in idx:
            sub.append(i); Mb *= primes[i]
            if Mb.bit_length() >= need_bits:
                break
        r, M = crt([residues[i] for i in sub], [primes[i] for i in sub])
        cos.append((M, r))
    return decode_coset_list(N, cos, d=d, metrics=mt)


# ---------------------------------------------------------------------------
# (iv) Legendre bits
# ---------------------------------------------------------------------------

def legendre_coset_list_size(primes: Sequence[int]) -> int:
    """|{r mod prod l : (r/l) = given signs}| = prod (l-1)/2 (every sign pattern)."""
    return prod((l - 1) // 2 for l in primes)


def labels_resolve(N: int, l: int) -> bool:
    """A Legendre label of the target factor separates the two roots of s mod l
    iff (N/l) = -1."""
    return legendre(N, l) == -1


# ---------------------------------------------------------------------------
# (v) e = N mod (p-1)  (Track B/D D3 linear complexity)  -> Pollard p-1 exponent
# ---------------------------------------------------------------------------

def decode_N_mod_pm1(N: int, e: int, bases: Sequence[int] = (2, 3, 5, 7, 11, 13),
                     metrics: Optional[Metrics] = None) -> Optional[int]:
    """E = N - e is a multiple of p-1, so a^E ≡ 1 (mod p).  Miller–Rabin style descent
    over the 2-part of E separates p from q when q-1 also divides E."""
    mt = metrics if metrics is not None else NullMetrics()
    E = N - e
    if E <= 0:
        return None
    o, t = E, 0
    while o % 2 == 0:
        o //= 2; t += 1
    for a in bases:
        mt.inc("gcd")
        g = gcd(a, N)
        if 1 < g < N:
            return g
        x = pow(a, o, N)
        mt.inc("rn_powmod")
        for _ in range(t + 1):
            mt.inc("gcd")
            g = gcd(x - 1, N)
            if 1 < g < N:
                mt.set("base", a)
                return g
            x = x * x % N
    return None


# ---------------------------------------------------------------------------
# (vi) random known bits of p and q (Heninger–Shacham-style branch and prune) + (a)
# ---------------------------------------------------------------------------

def hs_branch(N: int, depth: int, known_p: Dict[int, int], known_q: Dict[int, int],
              max_width: int = 1 << 22, metrics: Optional[Metrics] = None) -> List[Tuple[int, int]]:
    """All (a, b) mod 2^depth with a b ≡ N (mod 2^depth) consistent with known bits."""
    mt = metrics if metrics is not None else NullMetrics()
    cur = [(1, 1)]
    if known_p.get(0, 1) != 1 or known_q.get(0, 1) != 1:
        return []
    for i in range(1, depth):
        nxt = []
        bit = 1 << i
        mod = bit << 1
        for a, b in cur:
            for x in ((known_p[i],) if i in known_p else (0, 1)):
                for y in ((known_q[i],) if i in known_q else (0, 1)):
                    a2, b2 = a + x * bit, b + y * bit
                    if (a2 * b2 - N) % mod == 0:
                        nxt.append((a2, b2))
        cur = nxt
        mt.hi("max_width", len(cur))
        if len(cur) > max_width:
            raise RuntimeError("width exceeded")
    return cur


def decode_random_bits(N: int, nbits: int, known_p: Dict[int, int], known_q: Dict[int, int],
                       depth: int, d: Optional[int] = None,
                       metrics: Optional[Metrics] = None, max_width: int = 1 << 20,
                       max_decodes: Optional[int] = None) -> Optional[int]:
    mt = metrics if metrics is not None else NullMetrics()
    cands = hs_branch(N, depth, known_p, known_q, max_width=max_width, metrics=mt)
    mt.set("width_final", len(cands))
    cos = [(1 << depth, a) for a, _ in cands]
    return decode_coset_list(N, cos, d=d, metrics=mt, max_decodes=max_decodes)
