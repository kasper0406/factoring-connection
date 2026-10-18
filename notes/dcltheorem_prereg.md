# dcltheorem: pre-registration (frozen before the pre-registered runs)

Agent `dcltheorem`, 2026-10-18. Committed together with `src/dcltheorem/*` and
`tests/test_dcltheorem.py` before any pre-registered run of
`src/dcltheorem/experiments.py`.

## 0. Disclosure of exploratory runs before this file

- Smoke runs of DCL-det (design C1) at n = 32, 40 (pair seeds 1000–1004) and at n = 40, 56, 64
  (seed 7). One run each of C1@80, C3@80, C2@64 and C2@72 (seed 9000). These seeds are outside
  the pre-registered ranges, **except seed 7 at n = 40, which lies in the `lemma` range**.
  - Seen: every output was correct; one mismatch in the harness check of Lemma 5/7 at n = 32.
    It was a harness bug: the closed form used the superset |δ| < m instead of the exact range
    δ ∈ [−i*, m−1−i*]. After the fix there were 0 mismatches.
  - Block-0 metric E: C1@64 ≈ 4.5·10⁶ (before the per-ℓ handle precomputation) and
    C1@80 ≈ 2.8·10⁷. List elements: C1@80 block 0 = 71 830.
- HH baseline (signsplit V0 with charged fmpz stage 3) on seed 3 at n = 40/56/64/72/80.
  E = 2.1·10⁴ / 2.6·10⁵ / 1.1·10⁶ / (exact match) / 1.2·10⁷, with list = HH model count.
- **So the finite-size constant is already known to favour HH:** DCL-det block 0 is ≈ 2.4× HH
  in E at 80 bits. What is not known is the slope over 48–80 bits.

## 1. Question

Is the deterministic, end-to-end DCL algorithm (DCL-det; notes/dcltheorem_results.md §2)
correct as the theorem says, and how does its real arithmetic cost scale against the
calibrated operational Harvey–Hittmeir (HH) algorithm? DCL-det uses base 2, four disjoint
prime bands with ℓ ≡ 3 (mod 4), roots by powering, a balanced exponent tree, per-giant
remainder-tree evaluation and gcd extraction. It sees only N.

## 2. Population and instances

- S_n = {(p,q) primes : p < q ≤ 2p, 2^{n−1} ≤ pq < 2^n}.
- Pairs are sampled uniformly from S_n by `harness.sample_pair(n, seed)`, rejection from a box.
- Note: this is the theorem's population, **not** `balanced_semiprime`, which uses top-2-bit
  primes and so q/p < 4/3.

## 3. Designs (frozen; `params.schedule_design`, n-only)

- **C1 (primary).** ℓ ≡ 3 (mod 4), K = 4 dyadic bands (2^kY, 2^{k+1}Y].
  - H = the smallest H with block-0 multiplicity λ₀ = |V|·Σ_{ℓ∈P₀} 1/ℓ ≥ 3.
  - Y minimises an n-only expected block-0 list cost on a 2^{1/4} grid.
  - This is the theorem's algorithm with the finite-size H-schedule λ₀ ≈ 3 instead of
    H = n·2^{εn/2}.
  - **Bridge to the theorem.** Theorem D′ allows any n-only (H, Y). The schedule
    corresponds to H ≍ √n. That is the polylog regime of the theorem, where the
    *large-sieve-proved* exceptional fraction is O(n²/H²) and does not tend to 0; the
    Poisson heuristic predicts e^{−λ}. The experiment therefore tests the algorithm and
    its cost exponent, not the proved exceptional rate.
- **C2.** As C1, with λ₀ ≥ ln(P/ln P) + 2, the heuristic full-coverage regime. Run only at
  n ≤ 64, 8 seeds.
- **C3.** As C1, but with all odd band primes, using Tonelli–Shanks and a charged
  least-nonresidue search.

## 4. Metric (frozen)

E counts n-bit modular-multiplication equivalents:

E = mulN + mul_small + poly_cost + ⌈log₂ n⌉·(gcdN + invN) + ⌈bigmul_bits/n⌉

- poly_cost = Σ (d+1)·⌈log₂(d+2)⌉ over every polynomial product and remainder, in the FFT
  model. Every square-and-multiply is charged, including the Bluestein chirp powers of the
  HH baseline (audit item 3).
- **DCL-det.** E_block0, E_first (the blocks up to and including the first that yields a
  factor), and E_full (all 4 blocks, the worst-case bound).
- **HH.** E of the as-run algorithm. Slopes use runs that reached stage 3 (status
  "collision"); exact-match early exits are reported but excluded.
- **List elements.** m + #giants for DCL and HH, as a secondary metric.

## 5. Experiments

- **lemma.** n ∈ {28, 32, 36, 40}, seeds 0–59, C1, all blocks recorded. Checks:
  - Lemmas 1–4 on every covering target;
  - the Lemma 5/7 equivalence (gcd = p ⇔ not FULL) on every covering target of every
    clean block;
  - "certified ⇒ success".
- **oper.** C1 and C3 at n ∈ {40, 48, 56, 64, 72, 80}, seeds 100–131; C2 at n ≤ 64, seeds
  100–107. Per instance:
  - per-block metrics, success, the class of the first hit (STRUCTURAL_P / STRUCTURAL_Q /
    ACCIDENTAL);
  - the null expectation of accidental hits (#giants·m·(1/ord_p γ + 1/ord_q γ));
  - the certificate and the lemma checks.
- **hh.** n ∈ {40, …, 80}, seeds 100–115, the same pairs. Parameters: the minimum over
  primorial M of signsplit `best_shared`.
- **exc.** C1 certificate components only (clean blocks, coverage, FULL counts) for 1000
  pairs per n ∈ {48, 64, 80, 96}.

## 6. Predictions and decision rules (frozen)

- **R-a (must hold; a violation refutes a lemma or exposes a bug).**
  - 0 incorrect outputs;
  - 0 lemma violations in `lemma` and `oper`;
  - every certified instance factors.
- **R-b (exponent, C1 vs HH).** Least-squares slope of log₂(median E) against n over
  n ∈ {48, …, 80}, for C1's E_first and for HH.
  - Predictions: C1 slope ∈ [0.16, 0.21]; HH ∈ [0.18, 0.23].
  - Verdict "DCL-det has a smaller measured exponent" if slope_C1 ≤ slope_HH − 0.01.
    "Not distinguishable" if |difference| < 0.01. "Larger" otherwise.
  - The same is reported for E_full and E_block0 (secondary).
- **R-c (constant).** Prediction: median E_first(C1)/E(HH) > 1 at every size ≤ 80 bits. The
  extrapolated crossover is reported but is not a claim.
- **R-d (C2).** Prediction: C2's slope is ≥ HH − 0.005, i.e. no advantage at heuristically
  full coverage, as for multimod R2.
- **R-e (exc).** Predictions:
  - certified fraction ≥ 0.99 at every n;
  - fraction of q with block 0 clean ≈ 1 − Σ_{ℓ∈P₀} 1/(ℓ−1), within ±0.02;
  - FULL covering targets ≤ 1% of covering targets.
- **Classes.** Reported, with no rule. The accidental rate is compared with its null.
