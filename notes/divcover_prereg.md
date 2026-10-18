# divcover: pre-registration (frozen before the search runs)

Agent `divcover`, 2026-10-18. This file is committed before any run of the experiment
drivers `src/divcover/run_exact.py` (E1) and `src/divcover/run_structured.py` (E2).
Before writing it I ran only smoke tests of the exact C search (`src/divcover/gapsearch.c`):
- a cross-check against a Python brute force (n ∈ {8,12,16}, shapes (2),(3),(4),(2,2),(3,2),
  both modes; 30/30 agree);
- timing runs at n ∈ {30,40,60,80}, L ∈ {6,10,16} (rank 1), and n = 60, L ∈ {10,16} (rank 2);
- one off-grid function smoke test of E2 at n = 200: crt_design for r = 1,2,3 with L ≈ 14
  (all covers verified; log H within 1% of mass/(r+1)), and poly_family(t = 10, r = 3, K = 3).

The smoke tests showed exact optima within about a factor 3 in H of the random-model
prediction below, for example n = 60, L = 16, rank 1: H* = 1422 against a prediction of 646.
That observation fixed the threshold in §3, which is set far from it (a factor 2 in log H,
not 3 in H).

## 0. Question

Umans–Wang (UW, arXiv 2511.10851) Conj. 5.1: GAPs S, T of rank n^{o(1)}, |S|,|T| ≤ n^β,
with elements ≤ exp(n^α) and prefactored differences, such that S − T covers [1,n]. Thm 5.5
then factors N in Õ(N^{max(α,β)/2+o(1)}). He–Sahai (arXiv 2608.06681) refute rank 1.
Question for the search: **at small n, do rank-r GAP covers (r = 1, 2, 3) show any structure
that beats a "random-residue" model?** A structure gain would be the only empirical sign of
a route to max(α,β) < 2/5. Under the random model the rank helps only through an (r+1)-th
root of a CRT mass, which forces α → 1.

## 1. Objects (frozen)

- GAP Q = {b + Σ_k x_k c_k : 0 ≤ x_k < L_k}, b ≥ 1, c_k ≥ 1 (WLOG: flip coordinates), shape
  L_1 ≥ … ≥ L_r ≥ 2, box size L = Π L_k, height H = b + Σ (L_k − 1) c_k = max element.
  Any GAP Q of box size L splits as S − T with |S|·|T| = L (UW Prop. 3.4 for the long side),
  so β = log L / (2 log n).
- Modes: **prime** (every prime p ≤ n divides an element) and **full** (every d ∈ [2,n] divides
  an element). All elements are positive, so "nonzero" is automatic.
- H*_r(n, L, mode) = min height over all rank-r shapes with Π L_k = L. This is the exact
  minimum, from the exhaustive C search (`gapsearch.c`: exact residue masks per target, a
  CRT pre-filter on the two most restrictive coprime targets, and only the ordering symmetry
  c_1 ≤ c_2 when L_1 = L_2).
- **Random model (prime mode).** P_p(shape) = Pr_{b,c uniform mod p}[0 ∈ Q mod p] (exact
  enumeration). #GAPs of height ≤ H ≈ H^{r+1}/((r+1)! Π(L_k − 1)). H_pred is the H at which
  Σ_shapes #GAPs·Π_{p≤n} P_p = 1 (`core.random_model_log_height`).

## 2. Experiments (frozen grid)

**E1 (exact, small n).**
- n ∈ {20, 30, 40, 50, 60, 70, 80, 90, 100}; L ∈ {4, 6, 8, 12, 16, 24, 32} with L ≤ n/2;
  ranks r ∈ {1, 2, 3} (every shape with Π L_k = L); modes prime and full.
- A (n, L, r) cell runs iff its prime-mode H_pred ≤ cap_r, with cap_1 = 2·10^6, cap_2 = 5·10^3,
  cap_3 = 300. Search bound Hmax = ⌈20·H_pred⌉ (prime) and ⌈200·H_pred⌉ (full).
  The time limit is 300 s per C call. A timeout or "none ≤ Hmax" is recorded as a censored
  lower bound and is never imputed.
- Raw output: `experiments/raw/divcover_exact.jsonl`. Config: `experiments/configs/divcover_exact.json`.

**E2 (structured families, larger n; prime mode).**
- (a) **CRT/LLL design.** n ∈ {10^3, 10^4}; r ∈ {1,2,3,4}; balanced shapes with
  L ≈ n^{0.5}, n^{0.6}, n^{0.7}, n^{0.8}. Each prime p > L_1 is assigned to a box point
  greedily: largest prime first, to the point with the least log-mass, among points where
  the designed element stays positive. Then LLL is run on the (r+1)-dimensional lattice
  {(b,c) : b + ⟨x_p, c⟩ ≡ 0 mod p ∀p} and the GAP is re-anchored so that it is positive.
  Coverage is verified exactly. Recorded: α_eff = log log H / log n, β_eff = log L / (2 log n).
- (b) **Polynomial-value / cyclotomic GAP.** Generators t^0, …, t^r with box [0,K)^{r+1}, so
  that S − T ⊇ {f(t) : f ∈ Z[X], deg ≤ r, 0 ≤ coeffs < K} (signed version: shift by
  −(K−1)/2·Σ t^i, which contains t^k − 1, t^a − t^b, …).
  - n ∈ {10^3, 10^4}; t ∈ {3, 10, 2^16+1, 1000003, 10^9+7}; (r, K) ∈ {(4,3),(6,3),(8,3),(5,5)},
    signed boxes.
  - Measured: the fraction of primes p ∈ (box, n] that divide a nonzero element, against the
    random expectation Σ_p (1 − (1 − 1/p)^{#distinct nonzero elements}) / #primes.
- Raw: `experiments/raw/divcover_structured.jsonl`. Config: `experiments/configs/divcover_structured.json`.

## 3. Success / signal criteria (frozen)

- **E1 structure signal:** some cell with n ≥ 40 and L ≤ n/4, prime mode, any rank, has
  log H*_exact ≤ 0.5·log H_pred. Otherwise E1 is recorded as "consistent with the random
  model". The constant-factor deviation seen in the smoke tests counts as consistent.
- **E1 exponent fits (descriptive only, NOT evidence about asymptotics):** for α ∈ {1/3, 1/2},
  L*(n; α) is the least grid L (over ranks ≤ 3) with H* ≤ exp(n^α). Fit log L* against
  log n by least squares on the uncensored points. Report β_fit = slope/2 per rank.
  Pre-stated caveat: at n ≤ 100, exp(n^{1/3}) ≤ 100, so this regime is pre-asymptotic.
- **E2 success** (a candidate for theory): a verified cover at n ≥ 10^3 with α_eff < 0.4 and
  β_eff < 0.4. **E2(b) structure signal:** a covered fraction ≥ 2× the random expectation
  for some (t, r, K) at n = 10^4.
- **Overall.** Without a signal, the result is "no small-scale evidence for the Strong
  Prefactored Divisor Conjecture below 2/5". A signal triggers a theory follow-up, a
  GPT-6.1-sol review and an independent Opus critic before anything is reported.
