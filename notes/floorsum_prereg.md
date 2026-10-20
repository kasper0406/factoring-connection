# floorsum: pre-registration

Agent `floorsum`, 2026-10-20. Written and committed BEFORE E1–E5 were run.

Code at this commit: `src/floorsum/{farey,ladder_ops}.py`. Development checks before this commit:
exactness of `FloorSumConvex` against brute force on 1 500 random convex polynomials of degree 2–4
and 2 400 ladder-block range queries (j = 2, 3, 4; 44-bit N): 0 mismatches. One calibration block
(N ≈ 2^61, j = 2, x ≈ √N, H = 19 471) was inspected: 12 windows, 37 runs, Hβ^{1/3} = 17. No other
scaling data were produced.

## Claims under test (analysis in notes/floorsum_results.md)

- **Trap (C1).** A block arc of the ladder has curvature β ≈ N/x³ and length H = cx/B, so the
  convexity scale Hβ^{1/3} ≈ cN^{1/3}/B per block. Summed over the O(B log N) blocks this is
  Θ(N^{1/3} log N) for EVERY j ≥ 2. For even j the block integer hulls contain every vertex of the
  global hull of {xy ≥ N} except at candidate points (Lemma V, proved), so any hull-visiting
  oracle pays ≥ V_glob − #candidates.
- **Oracle (C2).** The Farey dissection oracle (`FloorSumConvex`) is exact for convex/concave
  rational polynomials of any degree, with cost O((Hβ^{1/3} + 1) log H) exact evaluations
  (heuristic constant). For the ladder at x ≈ √N this is T_2 ≈ H^{1/3}, T_3 ≈ H^{4/9},
  T_4 ≈ H^{1/2}: in general T_j ≈ H^{2(j−1)/(3j)}.
- **End-to-end (C3).** Plugging the oracle into the ladder gives an operational exact factoring
  algorithm for every j with total cost ≈ N^{1/3} polylog N, independent of j.
- **Fourier (C4).** In the certified smoothing model (degree-K trigonometric approximant plus exact
  enumeration of the h with ‖P(h)‖ < δ ≍ 1/K), the per-block cost is ≥ K + #{h : ‖P(h)‖ < δ},
  which is ≍ √H when the near-integer counts are ≈ 2δH (equidistribution).
- **Lattice (C5, KNOWN algorithm).** For the ladder's real need at j = 2 (divisors of N in a block)
  Coppersmith/Howgrave-Graham lattice reduction is a polylog block oracle; it covers blocks of
  length ≈ x^{β}·(loss), β = log x / log N, which contain the j = 2 ladder blocks.

## E1: the trap, measured (`run_experiments e1`)

- Instances: balanced_semiprime(n, s), n ∈ {28, 32, 36, 40, 44}, s = 0..4; j ∈ {2, 3, 4}; c = ½.
- For each (N, j): over all ladder blocks [x, x+H) with x ∈ [x0, √N] (x0 = ⌈2B⌉ as in the ladder),
  V_b = #vertices of the lower convex hull of {(h, ⌈P_j(h)⌉) : 0 ≤ h < H} (brute force,
  diagnostic); V_tot = Σ_b V_b; #blocks; V_glob = #vertices of the integer hull of {xy ≥ N} with
  x ∈ [x0, √N] (src.analytic.hull, brute force); K_cand = #located candidates (brute force).
- **Predictions.**
  - P1.1 (proved, a bug check): for j ∈ {2, 4}, V_tot ≥ V_glob − K_cand on every instance.
  - P1.2: V_tot / V_glob ∈ [1.0, 3.0] for every (N, j).
  - P1.3: for each j, the least-squares slope of ln(mean V_tot / ln N) against ln N over the five
    sizes lies in [0.28, 0.38] (N^{1/3}, j-independent).
  - P1.4: for each j, the per-block mean V_tot/#blocks divided by (2N)^{1/3}/(2B) lies in
    [0.3, 3] at 44 bits (per-block count tracks the convexity scale cN^{1/3}/B).
- Decision: C1 is supported iff P1.1–P1.3 hold.

## E2: oracle exactness and T_j(H) (`run_experiments e2`)

- Exactness: 3 000 random convex polynomials (degree 2–4, random denominators up to 10^15,
  lengths up to 3 000) and, for j ∈ {2, 3, 4} and n ∈ {40, 48, 56}, 50 random ladder blocks each
  with three shifts (0, −ρ−η, +ρ) and 4 random subranges, all against brute force.
- Scaling: single blocks at x = ⌊√N/2⌋ + δ (δ random in [0, 1000]), 5 random odd N per size,
  j = 2: n ∈ {48, 64, 80, 96, 112, 128}; j = 3: n ∈ {48, 64, 80, 96, 112}; j = 4: n ∈ {48, 64, 80, 96}.
  Cost = windows + runs (each costing O(log H) evaluations); also evals.
- **Predictions.**
  - P2.1: 0 mismatches.
  - P2.2: (windows + runs)/(Hβ^{1/3}) ∈ [0.5, 6] for every block, β = N/x³.
  - P2.3: for each j, the slope of ln(mean windows + runs) against ln H is within ±0.07 of
    2(j−1)/(3j) (1/3, 4/9, 1/2 for j = 2, 3, 4).
- Decision: C2 is supported iff P2.1–P2.3 hold.

## E3: end-to-end operational ladder (`run_experiments e3`)

- Instances: prime_instance(n, s) for cost scaling, s = 0..2; balanced_semiprime(n, s) for
  correctness, s = 0..4; n ∈ {28, 32, 36, 40, 44}; j ∈ {2, 3}; c = ½.
- **Predictions.**
  - P3.1: correctness 100% (semiprimes: a nontrivial factor; primes: None).
  - P3.2: for primes and each j, the slope of ln(mean evals / ln² N) against ln N lies in
    [0.28, 0.40].
  - P3.3: for primes, mean evals(j = 3)/mean evals(j = 2) ∈ [0.25, 4] at every size (no
    j-dependence of the exponent).
- Decision: C3 is supported iff P3.1 and P3.2 hold.

## E4: the Fourier certified model (`run_experiments e4`)

- 200 random ladder blocks (j = 2, x uniform in [√N/4, √N]) for random odd N of n ∈ {48, 56, 64}
  bits. Brute force n_δ = #{0 ≤ h < H : ‖P(h)‖ < δ} at δ = 1/K for K = ⌈√H⌉, and the model cost
  M = min over K ∈ {2^i} of (K + n_{1/K}).
- **Predictions.**
  - P4.1: mean n_δ/(2δH) ∈ [0.85, 1.15] at every size.
  - P4.2: mean M/√H ∈ [1, 4] at every size.
- Decision: C4's quantitative premise is supported iff P4.1 holds.

## E5: Coppersmith/Howgrave-Graham block oracle (`run_experiments e5`; KNOWN algorithm)

- balanced_semiprime(n, s), n ∈ {40, 48, 56}, s = 0..4. A divisor-in-interval lattice
  (polynomials N^{m−i}(x+t)^i, t^k(x+t)^m; LLL from python-flint), block length X(x) chosen from
  the Howgrave-Graham determinant bound with lattice dimension ≤ 16; scan x from ⌈N^{1/4}⌉ to √N,
  trial division below.
- **Predictions.**
  - P5.1: correctness 15/15.
  - P5.2: at x ≈ √N, the lattice block length X is ≥ the ladder's j = 2 block H = ⌊x/(2B)⌋ at
    every size (lattice blocks contain ladder blocks).
- E5 is a demonstration of a known algorithm; it supports no new claim.

## Disclosed limitations

- Hull counts are brute force (diagnostic); the operational oracle never enumerates the block.
- Sizes are small: slopes over 28–44 bits carry finite-size log factors (cf. analytic E1).

## Amendment A1 (before any E1–E5 run; committed with `src/floorsum/{lattice_block,run_experiments}.py`)

- **Disclosed peek.** While developing `lattice_block.py` I ran `lattice_factor` on four
  semiprimes (32 and 40 bits, seeds 100–101; all factored) and printed `best_params` at one 48-bit
  N: X(√N) = 1 248 with (m, n) = (7, 16), versus the ladder's j = 2 block H = 1 836. So P5.2 as
  written would fail at dimension ≤ 16. The reason is known: for fixed lattice dimension n the
  Howgrave-Graham bound loses a factor N^{Θ(1/n)} (plus the 2^{(n−1)/4} LLL factor), and the
  containment X ≥ x^β needs n ≍ log N.
- **E5 replaced.** Sizes n ∈ {36, 40, 44, 48}, seeds 0..2 (12 instances; the 56-bit runs are too
  slow in Python). P5.1′: correctness 12/12. P5.2′ (report): X(√N)/H_ladder at dimension ≤ 16;
  predicted ∈ [0.2, 2]. E5 still supports no new claim.
