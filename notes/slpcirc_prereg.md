# slpcirc pre-registration (2026-10-20)

Agent `slpcirc` (phase 17). Thread: non-relational factoring via short straight-line programs (SLPs).
This file is committed BEFORE any confirmatory run. Development/debugging uses only SLP size τ ≤ 4
(disclosed); the confirmatory runs below are E1 (proof check), E2 (census, τ ≤ 6) and E3 (sparse
dichotomy check).

## Objects

- Constant-free division-free SLPs over Z: start from {1, x}; each step applies +, −, × to two earlier
  values. τ(F) = least number of steps. All polynomials with τ ≤ 6 are enumerated exactly by breadth-first
  search over sets of computed values (development count: 133 742 polynomials incl. constants; max
  degree 64). Code: `src/slpcirc/enum_slp.py`.
- Window of primes W = all primes 67 ≤ p ≤ 251 (37 primes). For F and p, r(F,p) = number of distinct
  roots of F mod p in F_p (F mod p not identically zero as a polynomial; otherwise the pair is skipped
  and counted separately).

## Factor classification (fixed in advance)

Each irreducible factor g of F over Z (flint `fmpz_poly.factor`) gets the first matching type:
1. LIN: deg g = 1 (rational root; "explicit").
2. CYC: the primitive part of g(s·x + c) is cyclotomic for some s ∈ {±1, ±2, ±1/2}, c ∈ {−3,…,3}
   (coset/translate of a set of roots of unity; G_m, group order p − 1).
3. KUM: x^k ≡ constant mod g(x + c) for some 1 ≤ k ≤ 128, c ∈ {−3,…,3} (Kummer/binomial, G_m).
4. CHEB: for some a ∈ {1, −1}, κ ∈ {1, 1/2}, c ∈ {−2,…,2}, every irreducible factor of
   u^d·g(κ(u + a/u) + c) (cleared of denominators) is cyclotomic (Dickson/Chebyshev; group orders p ± 1).
5. LOW: deg g ≤ 4, none of the above ("explicit algebraic").
6. OTHER: everything else.
A root of F mod p is attributed to the first type (in this order) among the factors it is a root of.
"Group-or-explicit" = LIN ∪ CYC ∪ KUM ∪ CHEB ∪ LOW.

## E1 — proof check of Lemma S1 (sparse root bound on a cyclic group)

Claim to be checked (proved in notes/slpcirc_results.md): if G ≤ F_p^* is cyclic of order n, ℓ | n is
prime, F = Σ_{i ≤ t} c_i u^{e_i} with c_i ≠ 0 and the e_i pairwise distinct mod ℓ, then
#{u ∈ G : F(u) = 0} ≤ (n/ℓ)·2(⌈ℓ/K⌉ − 1) with K = ⌊(ℓ − 1)^{1/(t−1)}⌋.
Exhaustive for t = 2, 3 and all primes p ≤ 43, G = F_p^* (WLOG e_1 = 0, c_1 = 1); random 2·10^5
instances for t = 4, 5 with p ≤ 211. **Prediction: 0 violations** (a violation refutes the lemma).

## E2 — census of tiny SLPs (τ ≤ 6) over W

Pre-registered predictions:
- **H1 (dense zero sets are group-structured or explicit).** Call (F, p) dense if r(F,p) ≥ p/8. In at
  least 95% of dense pairs, ≥ 90% of the roots are group-or-explicit. Every dense pair failing this is
  listed and examined by hand.
- **H2 (no small-Frobenius anomaly for non-group factors).** For each distinct OTHER factor g with
  deg g ≥ 8 occurring in the census, let m(g) = mean over p ∈ W (p ∤ disc g · lc g) of #roots(g mod p).
  Prediction: the mean of m(g) over all such g lies in [0.7, 1.3] (Chebotarev: one root per factor on
  average), and no g has m(g) > 3. Every g with m(g) > 2 is listed and examined.
- **H3 (density ceiling).** The maximum over F of min_{p∈W} r(F,p)/p is attained by a polynomial whose
  dense roots are CYC/KUM/CHEB/LIN (not OTHER).
Failure of any prediction is reported as such; it is evidence, not proof, either way. Tiny SLPs have
degree ≤ 64, so this census only probes p up to ~4·deg; it cannot test asymptotic claims.

## E3 — sparse dichotomy check (Theorem SP, per-prime form)

For X = 2^14 and every prime p ∈ (X, 2X]: F ranges over (a) x^M − 1, M = lcm(1..B), B ∈ {50, 200, 1000};
(b) 200 seeded random trinomials and 4-nomials with exponents < 2^40 and coefficients in [−5, 5];
(c) x^{M}+x^{M'}−2 with M, M' smooth. Compute δ_p = r(F,p)/p via exponent reduction mod p − 1 and
gcd with x^{p−1} − 1 (plus x = 0). Check the implication
"δ_p ≥ λ (λ = 0.05) ⇒ every prime ℓ ≥ ℓ_0(t, λ − 1/p) dividing p − 1 divides Λ = ∏_{i<j}(e_i − e_j)".
**Prediction: 0 violations** (a violation refutes Theorem SP). Also report the fraction of p with
δ_p ≥ λ for each F (expected: smooth-(p−1) frequency for (a)/(c), ≈ 0 for (b)).

Resources: single core, < 2 GB RSS, `memory_pressure` checked first. Raw output:
experiments/raw/slpcirc_E{1,2,3}.jsonl.

## Addendum A (2026-10-20, before running E4; E1–E3 already run)

Motivated by the codex blind round (codex/slpcirc_1.out, open question: dense fibres of logarithmically
many quadratic composition layers with small fixed coefficients, outside group-induced families).

**E4 — fibres of quadratic compositions.** F = f_k ∘ … ∘ f_1, f_i(x) = x² + c_i, c_i ∈ {−2,…,2}
(SLP size 2k, constants from a fixed small set, N-independent). Fibre T(F,p) = {x ∈ F_p : F(x) = 0},
computed exactly by backward preimage propagation mod p. Primes: 60 primes drawn with seed 17 from
(2^12, 2^13]. (a) all 5^8 sequences with k = 8 (exhaustive, DFS sharing suffixes); (b) 10^4 seeded random
sequences with k = 39 (≈ 3 log2 p). Predictions:
- P4a: max over all k = 8 sequences of the median (over the 60 primes) of |T| is ≤ 16.
- P4b: for random k = 39 sequences, the mean of |T| over (sequence, prime) lies in [0.5, 2]
  (critical branching: each point has 1 + χ(s − c) preimages).
- P4c: no sequence in (a) or (b) has |T| ≥ p/(log2 p)² (≈ 25–30) at ≥ half of the 60 primes.
Top sequences are listed and checked for group structure (c ≡ 0 layers = powers, c = −2 = Chebyshev).
Raw: experiments/raw/slpcirc_E4.jsonl.
