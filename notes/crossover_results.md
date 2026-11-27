# crossover: hidden bits vs b-adic precision (2026-11-27)

**Verdict: no cheap middle band exists, in any base. L[1/3] is not broken.** Every count turns hidden at a fixed small
depth, and its first hidden bit is a quadratic character, or for KH the factor-order bit. Past that point a family
is either bounded (character type: saturates below 0.6 bits) or digit type (about log₂ b bits per digit). The
hidden digits are digits of p+q or q, i.e. the Frobenius of p in Q(ζ_{b^k}), and no family has a native evaluator
below about √N at that depth. H1: holds except for KH (its first hidden bit is the factor-order bit, not a
character). H2: holds. H3 (the hope): false.

Setup: prereg 9e6c3b9; experiments/crossover/crossover.py; results.json. N = pq with p, q in [2^13, 2^14).
200k samples, and 20k for S(N), which the modrec hull evaluator computes from N alone at 20 ms per N. Y = collision
conditional entropy given N mod b^{k+1} and five Jacobi symbols. The control is 0 everywhere. Entries at the
largest k in base 5 are biased low because the cells are sparse.

| Family | b=2, k=1..8 | b=3, k=1..5 | Type |
|---|---|---|---|
| r2, rQ2 | 0 until k=5 / k=4, then 0.41 flat | 0.41 flat | character (the Blum bit) |
| R23 (x²+xy+6y²) | 0, 0, 0.27, 0.34 flat | 0.34 flat | character (cubic class) |
| S(N) floor8, N-only | 0.21, 0.36, 0.58 flat (k≤5) | 0.42, 0.58 | character; costs N^{1/3} at every k |
| r4 = 8σ | 0 until k=5, then 0.41, 1.00, 1.66 | 0.41, 1.26, 2.35, 3.51, 4.27 | digit type (p+q) |
| r8 = 16σ₃ | 0 until k=6, then 0.41, 0.99 | 0.41, 0.41, 1.26, 1.70, 2.54 | digit type |
| τ | 0, 0, 0.41, 1.00, 1.68, 2.41, 3.16, 3.90 | identical to r4 (Eisenstein congruences) | digit type |
| KH = 2(N+q) | 0, 0, 1.00, 2.00, 2.99, 3.96, 4.88, 5.58 | 1.00, 2.58, 4.12, 5.33 | digit type, from the factor-order bit |

Base 2 only delays the start, because of the forced 2^ω and orbit factors. Odd bases turn hidden at k=1 with the
same kind of bit, (−3/p) or (5/p). Digit families need about n/4 hidden digits for Coppersmith (on low bits) to
finish. The cheapest known way to reach them: KH costs √N in class numbers; r4, r8 and τ need σ(N) or Frobenius at
p. Caveat: "public" means these particular invariants. The experiment locates the information; it is not a lower
bound on native algorithms.
