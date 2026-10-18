# ecmavg pre-registration (frozen before any pre-registered run)

Agent `ecmavg`, 2026-10-18. This is codex r18 §5's "first concrete step", extended. The code is
frozen together with this file: `src/ecmavg/{curves,labels,experiments}.py`, tested by
`tests/test_ecmavg.py`. The tests check three things: the vectorised labels against an affine
group-law brute force, truncation against recomputation, and the operational prefilter against
the labels on 150 pairs. No pre-registered statistic has been computed. The only data touched
so far is test data: primes below 1400, plus 150 pairs at 12 bits that are used solely for the
operational-consistency test, whose statistics are not inspected.

## What is measured

- **Families (fixed, explicit, J_max = 64).**
  - `suyama`: σ = 6, …, 69. u = σ² − 5, v = 4σ, x₀ = u³, z₀ = v³,
    a24 = (v−u)³(3u+v)/(16u³v). The curve is Montgomery and the arithmetic is x-only.
  - `mont` (control): A = j + 3 (j = 0..63), x₀ = 2.
- **Cutoffs.** For b-bit primes (b = ⌈n/2⌉ for the operational N), B = ⌊2^{bβ}⌋ with
  β ∈ {1/2, 1/3, 1/4}. β = 1/4 is codex's B = N^{1/8}.
  - Stage-1 prime powers are l^e ≤ X_b = 2^b + 2^{⌈b/2⌉+1} + 1.
- **Exposure label (exact marked-point label).** k_p(j) = 0 if p divides the degeneracy integer
  D_j. Otherwise k_p(j) = i if Z first vanishes mod p after the i-th prime-power step (i.e. i is
  the index of P⁺(ord P_j mod p) among the primes ≤ B). Otherwise k_p(j) = INF (the order is
  not B-smooth). The exposure indicator is X_j(p) = [k_p(j) ≠ INF].
- **Statistics, over all unordered pairs of distinct primes in a sample.** Same bit length,
  hence p < q < 2p. The product bit length is not restricted (disclosed).
  - Zu_J = #{j < J : X_j(p) ≠ X_j(q)} is codex's Z_N (one gcd per curve).
  - Zs_J = #{j < J : k_p(j) ≠ k_q(j)} is the staged version (a gcd after every prime-power
    step). It is exactly what the operational prefilter separates.
  - Per J ∈ {8, 16, 32, 64}: E Z/J, Var Z/J and P(Z = 0) for both Zu and Zs.
  - The staged failure is split into "both primes unexposed by all J curves" and
    "same nonzero signature".
  - Single prime: W_p(J) = Σ_{j<J} X_j(p), with E W/J, Var W/J, P(W = 0), the mean
    off-diagonal covariance c̄ of (X_j)_j over primes, and the between-class share of
    Var W_J for classes p mod 120.
  - Signatures: the fraction of pairs with identical unstaged (resp. staged) 64-curve
    signatures, the number of distinct unstaged signatures, and the fraction of primes with
    the all-zero signature.
  - R = [Var(Z_64)/64]/[Var(Z_8)/8] and the log–log slope of Var(Z_J) against J, for Zu, Zs
    and W.
- **Samples (seed 20261018 + b).**
  - S1 is codex's experiment: 128 primes at each b ∈ {12, 14, 16}.
  - S2: all primes at b ∈ {12, 14, 16}.
  - S3 is the extension: 2048 sampled primes at b ∈ {18, 20, 22}.
- **Commands.**
  - `.venv/bin/python -m src.ecmavg.experiments labels`
  - `.venv/bin/python -m src.ecmavg.experiments analyze`
  - Outputs: `experiments/raw/ecmavg_labels.jsonl` and
    `experiments/analysis/ecmavg_summary.json`.

## Kill criteria and predictions (frozen)

The theory says Var Z_N ≤ J·N^{o(1)} (EC2) is what Chebyshev needs. Below, "toy sizes" means
b ≤ 22. No outcome at these sizes proves or refutes an asymptotic statement. A kill means
"refuted for this family at toy sizes" only.

- **K1 (literal EC2, codex's criterion "variance growing like J²").**
  - Statistic: R_u for `suyama` on S2 at b = 14 and b = 16, and on S3 at b = 22, for
    β = 1/3 and β = 1/4.
  - If R_u ≥ 4 in every one of these cells: KILL literal EC2 for this family at toy sizes.
  - If R_u ≤ 2 in every cell: CONSISTENT with EC2.
  - Otherwise: INCONCLUSIVE.
  - Independence would give R_u ≈ 1. Pure J² growth would give R_u ≈ 8.
- **K2 (route-level, at toy sizes).** For `suyama` on S2, b = 16, β = 1/4: if
  fail_s(64)/fail_s(8) ≥ 1/2, the staged failure does not decay with J. Verdict: KILL the
  route for this family at toy sizes.
- **Predictions (mine, before running).**
  - **Q1.** K1 is INCONCLUSIVE, with 1.1 ≤ R_u < 4 in most K1 cells. The motivation is a
    common-mode heuristic: fluctuations of the family-average exposure rate with p, e.g.
    through p mod small primes, add a J²·Var_p(ḡ) term. Its relative size is ~ 1/log B.
  - **Q2.** The mean off-diagonal covariance c̄ > 0 in every S2/S3 cell, for both families
    and every β.
  - **Q3.** K2 does not fire. fail_s(64)/fail_s(8) < 0.2 in every S2/S3 cell with β ∈ {1/2, 1/3}.
  - **Q4.** At J = 64, ≥ 90% of staged failures are "both unexposed" in every S2/S3 cell with
    β ∈ {1/3, 1/4}. Same-nonzero-signature coincidences are rare.
  - **Q5.** The `mont` control family shows the same qualitative picture: its R_u lies within
    a factor 2 of `suyama`'s in every K1 cell.
  - **Q6.** R_u − 1 for `suyama`, β = 1/3, changes by less than a factor 3 between b = 12 (S2)
    and b = 22 (S3). The common mode decays slowly, if at all, over this range.
- **Exploratory, not pre-registered:** anything else, e.g. within-class variances, other
  conditioning sets, and larger b.

## Deviations

Any deviation will be recorded here as a dated, disclosed addendum.
