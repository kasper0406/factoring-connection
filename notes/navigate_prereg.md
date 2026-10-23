# navigate: pre-registration (agent `navigate`, phase 25, 2026-10-23)

Written and committed before any registered statistic is computed. The only computation done before this file was
a resource calibration (`scratchpad/navigate/calib.py`): the *number* of signed y-smooth exponent vectors of weighted
ℓ1 size ≤ β ln N for y ∈ {13, 19}, N ≈ 2^18, 2^20 (to size E2). No count per class, no smoothness rate, no Fourier
coefficient was computed. The codex blind round (`codex/navigate_1.out`) suggested E2's "does conditional information
survive near the end" question; it was read before this file.

Target (smoothbuild §7, W11; "NAV"): given odd N and a uniform square a (root withheld), output y-smooth u, v with
u ≡ a·v (mod N), ln|u| + ln|v| ≤ 3 ln N. The experiments are toy-scale mechanism checks. They cannot test L-scale
asymptotics, and no outcome will be read as evidence for or against an L[1/3] algorithm.

## E1: toy navigation algorithms vs sampling (operational; sees only N, a, factor base, seed)

**Instances.** Balanced semiprimes from `src/common/instances.balanced_semiprime(bits, seed)`, bits ∈ {40, 48, 56},
seeds 0–5. Factor base FB = primes ≤ y = ⌊2^{bits/4}⌋ (so the Gauss pieces have u = ln√N/ln y ≈ 2 at every size).
100 targets per N: a = r² mod N, r uniform in [1, N − 1] (harness; gcd(r, N) checked; r is never passed on).

**Cost units.** A *test* is one smoothness test of one nonzero integer (the second piece of a split is tested only if the
first is smooth). A *split* is one rational reconstruction (`gauss_split`). Each (algorithm, target) runs until it has
used 400 tests.

**Valid output.** Integers u, v ≠ 0 with u ≡ a·v (mod N), |u| and |v| y-smooth, gcd(v, N) = 1, ln|u| + ln|v| ≤ 3 ln N.
Signs are allowed (−1 is a column, as in EQ1); the fraction with u, v > 0 is reported descriptively. Every output is
verified exactly. *Distinct* outputs: distinct reduced rationals |u|/|v| with sign, per target.

**Algorithms.**
- **S0 (sampling).** Repeat: w = ℓ₁ℓ₂ with ℓ₁, ℓ₂ independent uniform in FB; t = a·w mod N; (u′, v′) = gauss_split(t);
  test |u′|, then |v′|; if both smooth output (u′, w·v′).
- **G1 (height-steered walk).** State (t, w) = (a, 1). Step: draw c = 16 uniform FB primes ℓ; for each compute the split
  of t·ℓ and its score ln|u′| + ln|v′|; move to the minimiser (first on ties): t ← tℓ, w ← wℓ; test the chosen split as
  in S0 and output (u′, w·v′) if both pieces are smooth. Reset to (a, 1) when ln w + ln y > 1.9 ln N.
- **J2 (Monte Carlo counting walk, JVV-style with the obvious estimator).** State (t, w) = (a, 1). Round: draw c = 8
  uniform FB primes ℓ; for each, run m = 16 S0-trials on target t·ℓ (fresh 2-prime multiplier w′ each; every success
  is recorded as an output (u′, w·ℓ·w′·v′)); score = number of successes; move to the argmax (uniform among ties);
  same reset rule as G1. The run stops at 400 tests even mid-round.

**Null (harness).** For each split tested by S0 or G1, draw an independent coprime pair (U, V) uniform among integers of
the same bit lengths as |u′|, |v′| (resample until gcd = 1) and test both. R_alg := (alg's success rate per split
tested) / (null success rate on the same splits).

**Statistics.** Pooled per size (all N, all targets): Y_alg = valid outputs / tests; Y_alg^dist = distinct valid
outputs / tests; R_S0, R_G1; mean tested log-size (ln|u′| + ln|v′|) for S0 and G1.

**Predictions (per size, each must hold at all three sizes).**
- P1a: R_S0 ∈ [0.85, 1.15]. (Gauss pieces of random targets are smooth like coprime random pairs of their sizes;
  smoothbuild's post-hoc coprime null gave 1.00–1.01.)
- P1b: Y_J2 / Y_S0 ∈ [0.85, 1.15]. (Monte Carlo completion counting adds no yield per test.)
- P1c: Y_G1 / Y_S0 ≥ 1.15. (Steering on the cheap height signal lowers the tested size.)
- P1d: R_G1 ∈ [0.85, 1.15]. (G1's gain is fully explained by the tested size: no residual "navigation" effect.)
- Descriptive: Y^dist ratios, splits per output, positive-sign fraction, mean tested log-sizes.

## E2: exact completion counts and the information a counting oracle must supply (diagnostic, factor-aware)

**Instances.** balanced_semiprime(20, seed), seeds 0–3. FB = primes ≤ y, y ∈ {13, 19}. The harness uses p, q only to
compute discrete logarithms (for the Fourier transform and the subgroup H = ⟨FB⟩ ≤ G = (Z/N)^×).

**Counts.** Enumerate every signed exponent vector e over FB with weighted size |e| := Σ|e_i| ln p_i ≤ β_max ln N; its
class is ∏p_i^{e_i} mod N (this is the coprime pair (u, v) with u/v ≡ class). For budget β: C_β(t) = #{e : |e| ≤ β ln N,
class(e) = t}. Grids: y = 13: β ∈ {1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0}; y = 19: β ∈ {1.0, 1.25, 1.5, 1.75, 2.0}.
All statistics are over t ∈ H (classes outside H have C ≡ 0); μ_β := mean of C_β over H.

**Statistics.**
- h_β := fraction of t ∈ H with C_β(t) ≥ 1; Poisson reference 1 − e^{−μ_β}.
- D_β := Var/mean of C_β over H (index of dispersion).
- Fourier spread: Ĉ_β(χ) for all characters χ of G (FFT in log coordinates); exclude characters trivial on H; sort the
  rest by |Ĉ|²; k75_β := least number of characters carrying ≥ 75% of their total |Ĉ|² mass; report
  f75_β := k75_β / #(characters nontrivial on H). An i.i.d. noise vector gives f75 ≈ 0.39 (exponential order statistics).
- Squares vs non-squares: h_β and μ_β separately on H ∩ G² and H ∖ G².
- Descriptive: the exact-oracle JVV profile: along uniformly sampled solutions of the full budget (canonical order:
  largest prime first), the fraction σ of next letters (prime and side) with a nonzero completion count, as a function
  of the remaining budget.

**Predictions** (each over all (seed, y, β) cells with μ_β in the stated range).
- P2a: |h_β − (1 − e^{−μ_β})| ≤ 0.05 for μ_β ∈ [0.1, 5].
- P2b: D_β ∈ [0.8, 1.6] for μ_β ≥ 0.5.
- P2c: f75_β ∈ [0.2, 0.5] for μ_β ∈ [0.2, 10] (no sparse character description of the counts).
- P2d: |h_β(squares) − h_β(non-squares)| ≤ 0.05 for μ_β ∈ [0.1, 5] (in cells where H ⊄ G²).

## Analysis rules
- A prediction fails if any cell violates its band; failures are reported with the cell. Post-hoc analyses are labelled.
- Raw JSONL: `experiments/raw/navigate_E{1,2}.jsonl`; summary `experiments/analysis/navigate_summary.json`.
- Code: `src/navigate/core.py` (operational), `src/navigate/e1_walks.py`, `src/navigate/e2_counts.py` (harness),
  `src/navigate/analyze.py`; tests `tests/test_navigate.py`.
