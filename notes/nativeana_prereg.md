# nativeana pre-registration (phase 29)

Agent `nativeana`, 2026-10-27. Written and committed before any experiment code exists or runs. The theory
(notes/nativeana_results.md, to be written) carries the main weight. These are small exact or numerical checks of
precision-versus-cost trade-offs. Toy sizes cannot show asymptotics, and any factoring success is a correctness
check only (rule 8).

## Common setup
- Instances: `src/common/instances.py::balanced_semiprime(n, seed)`, seeds 0, 1, 2, … in order. Where a congruence
  filter is stated (N ≡ 1 mod 4, so that D = −4N is fundamental), failing seeds are skipped and recorded.
- Operational estimators (`src/nativeana/estimators.py`) see N (or m) only. Ground truth (p, q, h(D) and Hurwitz
  class numbers from PARI via cypari2) lives in `src/nativeana/harness.py` and is never imported by operational code.
  PARI's class-number routines factor D internally, so they are harness-only.
- Arithmetic: χ_D(n) and all identities are exact (integers / Fractions). L-value approximations, least-squares fits
  and slopes use mpmath (50 digits) or float64 where stated; this is allowed for fits (rule 1).
- Raw output: `experiments/raw/nativeana_E{1,2,3}.jsonl`. Memory ≤ 2 GB; one job at a time; no background jobs.
- "Slope" means the least-squares slope of log(median over instances) against log(parameter) over the stated range.

## E1: precision versus cost for h(−4N)
- n ∈ {28, 36, 44}, 20 instances each, N ≡ 1 (mod 4). Q := |D| = 4N; χ = χ_{−4N}; h = (√Q/π)·L(1, χ).
- Estimators of h (each returns ĥ; error Δ = |ĥ − h|):
  - EP(x): truncated Euler product over primes ℓ ≤ x, x = 2^j, j = 4..20.
  - PS(T): sharp partial sum Σ_{n≤T} χ(n)/n, T = 2^j, j = 4..min(20, ⌊log₂(√Q/4)⌋).
  - SM(T): smoothed sum Σ_n χ(n)e^{−(n/T)²}/n (summed to n ≤ 6T), same T grid as PS.
  - AFE(T): the exact erfc formula h = Σ_n χ(n)[erfc(n√(π/Q)) + (√Q/(πn))e^{−πn²/Q}] truncated at n ≤ T, T = 2^j up to
    2^⌈log₂(8√Q)⌉.
  - RS(K): Srinivasan-type random summation with K samples (cutoff T* = 2^12·Q, dyadic blocks chosen uniformly, n
    uniform in its block; seeded), K = 2^j, j = 6..16.
- **P1.1** EP slope over x ∈ [2^6, 2^20] lies in [−0.65, −0.35] at each n.
- **P1.2** PS slope and **P1.3** SM slope over their T grids (T ≥ 2^6) lie in [−0.65, −0.35] at each n.
- **P1.4** AFE: Δ/h < 10^{−6} for every instance once T ≥ 4√Q. For T ≤ √Q/16 the median of Δ_AFE/Δ_PS lies in
  [1/3, 3].
- **P1.5** RS: slope over K lies in [−0.6, −0.4], and the median of Δ·√K/(√Q·B) lies in [0.05, 0.5], where B is the
  number of dyadic blocks.
- **P1.6** Model cost C = (cost) + √(2Δ), with cost = x (EP), T (PS, SM, AFE) or K (RS). Report log C_min / log N per
  estimator. Prediction for EP: in [0.15, 0.30] at every n. This uses the true (harness) error, so it is a
  retrospective score, not an algorithm.

## E2: optimal linear rules (numerical companion of the linear-model barrier)
- Q = 2^30. Family A: m uniform in [Q, 2Q] with m ≡ 1 (mod 4) (any m; D = −4m, possibly non-fundamental), drawn from
  `random.Random("nativeana|E2|A|<i>")`. Family B: `balanced_semiprime(31, seed)` with N ≡ 1 (mod 4) and N ∈ [Q, 2Q]
  (skips recorded). 6000 training and 3000 test members per family.
- Target y = L(1, χ_{−4m}) = π·h(−4m)/(2√m) (h from PARI, harness). Features χ_{−4m}(n) for odd n ≤ 2K (K features),
  K ∈ {2, 4, …, 512}.
- Rules, fitted by least squares on training data, scored by test RMS:
  - (i) constant weights (plus an intercept);
  - (ii) weights affine in u = (m − Q)/Q (2K + 2 parameters).
- Log-linear companion (Euler-product model): fit log L(1, χ_{−4m}) by least squares on χ_{−4m}(ℓ) and χ_{−4m}(ℓ)² for
  primes ℓ ≤ x, K = π(x) ∈ {2, …, ≈ 500}.
- **P2.1** Slope of log test-RMS against log K, over K ∈ [8, 512], lies in [−0.65, −0.35] for both families and both
  rules.
- **P2.2** Rule (ii) test RMS ≥ 0.9 × rule (i) test RMS at every K, in both families.
- **P2.3** Test RMS of family B / family A lies in [0.7, 1.4] at every K ≥ 8.
- **P2.4** Test RMS·√K lies in [0.1, 1.0] at every K ≥ 8 (heuristic constant ≈ 0.4).
- **P2.5** Log-linear companion: slope of log test-RMS(log L) against log K over K ∈ [8, ≈ 500] lies in
  [−0.65, −0.35], both families.

## E3: Kronecker–Hurwitz amortisation
- n ∈ {16, 20, 24}, 20 instances each, no congruence filter. H(n) from PARI `qfbhclassno` (harness), exact rationals.
- **P3.1** For every instance, exactly: Σ_{|t|<2√N} H(4N − t²) = 2(N + q); Σ (t² − N)H(4N − t²) = −2(1 + p³);
  Σ (t⁴ − 3t²N + N²)H(4N − t²) = −2(1 + p⁵) (Eichler–Selberg at weights 2, 4, 6, where S_k(1) = 0).
- **P3.2** Local densities (the ingredient of the amortisation proposition): for ℓ ∈ {3, 5, 7}, k ∈ {1, 2, 3} with
  ℓ^{4k} ≤ 2·10^7, and every a mod ℓ^k with v = v_ℓ(a) < k, the number of M ∈ M₂(Z/ℓ^k) with det M ≡ a equals
  ℓ^{3k}(1 − ℓ^{−2})·Σ_{j=0}^{v} ℓ^{−j}, by brute force.
- **P3.3 (exploratory)** C^∞ bump weight w(u) = exp(−1/(1 − u²)) on u = t/(2√N): S_w = Σ_t w·H(4N − t²) and
  A_w = (σ(N)/N)(1/π)Σ_t w·√(4N − t²). Report R_w = S_w − A_w, |R_w|/√N, and the R² of R_w against the hyperbolic
  prediction −2p·J_w(√(p/q)) with J_w(r) = (2/π)∫ w(u)√(1 − u²)/(1 − 2ur + r²) du. Prediction: median |R_w|/√N in
  [0.05, 20] at every n. The R² is reported, not predicted.

## Tests (not experiments)
- The Hankel window bound is checked in `tests/test_nativeana.py` on random and least-squares-optimised instances:
  a sum of T exponentials (nodes anywhere in the closed unit disc) differs from δ_N by ≥ 1 in ℓ¹ on the window
  |n − N| ≤ W whenever T ≤ W; and T = W + 1 equally spaced nodes on the unit circle reproduce δ_N exactly on that
  window, so the bound is tight.
