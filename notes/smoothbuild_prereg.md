# smoothbuild pre-registration (phase 24)

Agent `smoothbuild`, 2026-10-23. Written and committed before any experiment below was run. Code lives in
`src/smoothbuild/`; raw output goes to `experiments/raw/smoothbuild_E{1,2,3}.jsonl`, summary to
`experiments/analysis/smoothbuild_summary.json`. The codex blind round (`codex/smoothbuild_1.md`) was launched before
this file was written; its output has not been read.

Every experiment is a toy-scale check of a structural claim. None of them is evidence about L-scale costs by itself.
Factorisations computed with python-flint in the harness play the role of the "pre-factored" oracle (angle 1) and are
harness code, never imported by operational code.

## Claims under test

- **C1 (Galois recurrence law).** For an irreducible F ∈ Z[x] of degree d with Galois group G acting on its roots, the
  number of pairs of family members F(a), F(a′) (a ≠ a′ in an interval of length A) sharing a prime p > A, summed over
  p ∈ (P₁, P₂], is ≈ (r_G − 1)·C(A, 2)·Σ 1/p², where r_G = number of G-orbits on ordered pairs of roots (the rank).
  For independent uniform random integers the factor is 1. r_G − 1 ≤ d − 1, with equality iff G acts regularly
  (e.g. cyclotomic Φ_n, abelian). So structured "pre-factored" families make large primes recur at most d − 1 times
  more often than random integers.
- **C2 (known-prime relaxation for random sides; CGPT).** For uniform random integers ≤ X, the first index T at which
  some subsequence product is a perfect square (any primes allowed, i.e. full known factorisations) is a constant
  fraction of J₀(X) = min_y π(y)·X/Ψ(X, y): Croot–Granville–Pemantle–Tetali (Annals 2012) prove
  T ∈ [(π/4)e^{−γ}, e^{−γ}]·J₀ (1 + o(1)) with probability 1 − o(1). Allowing all primes beats the best single
  smoothness bound by at most a constant factor.
- **C3 (constructed side is free, tested side follows its size).** In residue-class construction (smooth Gauss
  splits of random quadratic-residue targets mod N), steering by free smooth multipliers or by placing a free
  factor with a 2-dimensional lattice changes the target, not the size distribution of the tested pieces; the
  per-test success rate equals that of size-matched random integers. Products of two smooth lists are flat mod N.

## E1: recurrence of large primes in pre-factored families (C1)

- Families (degree, Galois group, rank r): x⁴+x+1 (S4, 2), x⁴+8x+12 (A4, 2), x⁴−2 (D4, 3), Φ₅ (C4, 4), x⁴+1 (V4, 4),
  x⁶+x+1 (S6, 2), Φ₇ (C6, 6), Φ₉ = x⁶+x³+1 (C6, 6). Groups as identified by sympy `galois_group`; ranks by
  Burnside on explicit permutation groups.
- Members F(a), a ∈ [a₀ⱼ, a₀ⱼ + A), A = 2¹⁵, a₀ⱼ = 2¹⁶ + j·2¹⁵, j = 0..4 (5 replicates).
- Control R: A uniform random integers in [2⁶³, 2⁶⁴), seeds 0..4.
- Primes p ∈ (P₁, P₂] = (2¹⁶, 2²²]. c_p = number of members divisible by p: by root sieving for polynomial families
  (cross-checked against flint factorisations on a subsample in the tests), by flint factorisation for R.
- Statistic S = Σ_p C(c_p, 2). Chebotarev-level prediction E = (r − 1)·C(A, 2)·Σ_p 1/p² (families),
  E_R = C(A, 2)·Σ_p (⌈2⁶⁴/p⌉ − ⌈2⁶³/p⌉)²/2¹²⁶ (control). Σ computed in exact integer fixed point.
- **P1a.** For each polynomial family, mean over the 5 replicates of S/E ∈ [0.80, 1.25].
- **P1b.** For R, mean S/E_R ∈ [0.90, 1.10].
- **P1c.** Ratios of replicate means: S(Φ₅)/S(S4) and S(x⁴+1)/S(S4) each ∈ [2.2, 4.0]; S(D4)/S(S4) ∈ [1.5, 2.6];
  S(Φ₇)/S(S6) and S(Φ₉)/S(S6) each ∈ [3.5, 6.5]; S(A4)/S(S4) ∈ [0.75, 1.33].
- Descriptive (no prediction): the arrangement-level ratio S / ((A²/2)·Σ_p ρ_F(p)(ρ_F(p) − 1)/p²).

## E2: residue-class construction vs sampling at matched size (C3)

- N: `balanced_semiprime(bits, seed)`, bits ∈ {40, 48, 56, 64}, seeds 0..5. y = 2^{bits/4} (so u = 4 for the
  canonical residue). Targets t = r² mod N with r a uniform unit (seeded); the constructors see only (N, t).
- M0 (sampling, Dixon): test the canonical residue t. 20 000 targets per N.
- M1 (Gauss split): t ≡ u/v with |u|, |v| ≤ ~√N from the first extended-Euclid remainder below √N; success iff
  |u| and |v| are both y-smooth. Same 20 000 targets.
- M2 (steering by free multipliers = placing a free factor by a 2-D lattice): for 2 000 targets per N, Gauss-split
  t·s_j for the K = 16 multipliers s_j = the first 16 primes; record per-split and per-target success.
- M4 (flatness of smooth products): A = A′ = the 1 024 smallest y-smooth integers ≥ 2^{bits/2 − 2}; for 3 targets t
  per N count pairs (a, a′) with (a·a′·t⁻¹ mod N) ∈ [1, W], W = ⌈400·N/|A|²⌉ (expected ≈ 400 per target).
- Nulls: M0, 20 000 uniform integers in [1, N) per N; M1/M2, per split one uniform integer of the same bit length as
  |u| and one of the same bit length as |v|, both tested.
- **P2a.** M0 pooled observed/null ∈ [0.85, 1.18] at each size.
- **P2b.** M1 pooled observed/null ∈ [0.75, 1.25] at each size.
- **P2c.** M2 pooled per-split success / M1 pooled success ∈ [0.90, 1.10] at each size.
- **P2d.** M2 per-target success / (1 − (1 − p̂)^16) ∈ [0.90, 1.10] at each size, p̂ = M2 pooled per-split rate.
- **P2e.** M4 pooled hits / (|A|²·W/N summed) ∈ [0.85, 1.18] at each size.
- Descriptive: relations per smoothness test for M0, M1, M2 vs bits.

## E3: first square dependency with all primes allowed (C2)

- X ∈ {2³², 2⁴⁰, 2⁴⁸}; seeds 0..11. a_i uniform in [1, X]. Each a_i is fully factored (flint); incremental GF(2)
  elimination on the odd-exponent support over all primes. T = first index with a dependency.
- On the same sequence, for each y in the grid y_k = ⌊2^{k/2}⌋, k = 8..2·log₂X − 4, T_y = first index at which the
  y-smooth members have a dependency (∞ if not reached by the cap 6·Ĵ₀); T_best = min_y T_y.
- Ĵ₀ = min_k π(y_k)/f̂(y_k), f̂ = fraction of y-smooth integers among 200 000 independent uniform integers ≤ X
  (seed 1000).
- **P3a.** Median over seeds of T/Ĵ₀ ∈ [0.35, 0.80] at each X.
- **P3b.** Minimum over seeds of T/Ĵ₀ ≥ 0.25 at each X.
- **P3c.** Median over seeds of T_best/T ∈ [1.0, 2.0] at each X.
- Descriptive: largest prime in the first dependency's support, relative to the y attaining Ĵ₀.

## Resource limits

One heavy job at a time, ≤ 2 GB RSS, `memory_pressure` checked first. Any deviation from this file will be disclosed
in the results note as a deviation.
