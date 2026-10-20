# almostall13: pre-registration of E1–E2

Agent `almostall13`, phase 14, 2026-10-20. Frozen before any run. Code is written after this file is committed and
frozen at the commit that first runs it. Neither experiment tests an asymptotic claim about factoring or proves
anything about relation supply. E1 is an exact check of a linear-algebra statement. E2 is a toy statistical comparison.

Hypothesis under test (to be proved or refuted in `notes/almostall13_results.md`). Averaging over the input N gives no
new free parameters to an NFS-type field beyond the d + 1 coefficients of F. Every descent statistic is a function of the
field pair (F, G) alone, so an average over (N, field) is an average over fields under a pushforward law.

Code (planned):
- `src/almostall13/rank.py`: E1, exact affine-hull dimensions over Q (python-flint / fraction-free elimination).
- `src/almostall13/harness_e2.py`, `src/almostall13/analyze_e2.py`: E2. DIAGNOSTIC/HARNESS code that reuses
  `src/nested/toy_e1.py` by import, uses only N from each instance, and never reads factor labels.
- Tests: `tests/test_almostall13.py`. Config: `experiments/configs/almostall13_e2.json`. Raw output:
  `experiments/raw/almostall13_e1.json`, `experiments/raw/almostall13_e2.jsonl`, summary
  `experiments/raw/almostall13_e2_summary.json`.

## E1. Exact parameter count under N-averaging (checks Prop NB)

**Setting.** Design (a‴) as implemented in `src/thmA/numfield.make_f`: F = F̂_{g,w} + G·R_c with F(g₀, g₁) = N·w.
- d ∈ {3, 5}.
- Points: s primitive points P_i = (a_i, b_i), drawn with seed `almostall13|E1|pts|d|s`, |a_i| ≤ 2^10,
  1 ≤ b_i ≤ 2^10, pairwise non-proportional; s ∈ {d − 1, d, d + 1, d + 2, 2d + 2}.
- N: 40-bit integers.
- g₀, g₁ ≈ 2^10, coprime.
- w: a prime in [2^5, 2^6].
- c ∈ [−2^12, 2^12]^d.

Regimes. Each draws K := 3(d + 3) + 10 samples, seeded.
- (A) (N, g, w) fixed; c varies.
- (B) (N, g) fixed; (w, c) vary.
- (C) g fixed; (N, w, c) vary, N uniform among 40-bit integers.
- (D) everything varies.

Statistics:
- a_F: the affine-hull dimension of {(F(P_i))_{i≤s}} ⊂ Q^s;
- a_FG: the same for the joint vectors {(F(P_i), G(P_i))_{i≤s}} ⊂ Q^{2s} (regime D only);
- an exact check that the integer kernel vector(s) λ of the s × (d + 1) evaluation matrix, computed from the points
  alone, satisfy λ·(F(P_i))_i = 0 for every sample in every regime (s > d + 1).

**Predictions.**
- R1: a_F = min(s, d) in (A); a_F = min(s, d + 1) in (B), (C) and (D).
- R2: a_F(C) = a_F(B) for every (d, s). Varying N adds no dimension once w varies.
- R3: a_FG(D) = min(s, d + 1) + min(s, 2).
- R4: for s > d + 1, every sample in every regime satisfies the kernel relations exactly (0 violations).

## E2. Toy one-step descendability: random N versus fixed N (checks the "pushforward" picture)

**Setting.** The toy of nested E1 (`src/nested/toy_e1.py`, unchanged, imported):
- d = 3, 44-bit balanced semiprimes (`src/common/instances.balanced_semiprime`), M_bits = 20;
- B = 2^10, Y = 2^14, S = 32;
- parents from rational ℓ ∈ [2^22, 2^22 + 2^16].

Two regimes, 72 fields each:
- FIX: 6 fixed N (instance seeds 101–106) × 12 fields each, drawn with field seeds `almostall13|FIX|<N-seed>|<j>`;
- RND: 72 distinct N (instance seeds 201–272), one field each, field seed `almostall13|RND|<N-seed>`.

Per field:
- n_u := 150 uniform algebraic targets and 150 uniform rational targets (uniform over the pools of `Uniform` in
  (B, Y], seeded), with their D1 values, giving p̂_alg(A) and p̂_rat(A);
- up to 80 parent nodes (`parent_node`, S₁ = S), with D1 for every child and one matched-uniform draw per child
  (same type, norm within a factor 1.25), as in nested E1.

**Statistics.**
- The between-field variance of p̂_alg: V_btw := sample variance − mean of p̂(1 − p̂)/(n_u − 1), and SD_btw := √V_btw
  (floored at 0).
- Jensen ratios over parents with ≥ 2 children:
  - J_pool := (Σ observed all-children-good)/(Σ ∏_children pooled-per-regime matched-uniform rate of the child's type);
  - J_field: the same with each field's own per-type matched-uniform rates.

**Predictions.**
- Q1: |mean_RND p̂_alg − mean_FIX p̂_alg| ≤ 0.05, and the 95% bootstrap CI of the difference lies in [−0.08, 0.08]. The
  bootstrap resamples fields in RND, and resamples N-clusters and then fields within clusters in FIX (2000 resamples).
- Q2: SD_btw(RND)/SD_btw(FIX) ∈ [0.67, 1.5]. N-averaging does not reduce the between-field heterogeneity.
- Q3: in FIX, a one-way ANOVA of p̂_alg across the 6 N-groups has p-value > 0.01. There is no N effect beyond the
  field. p is computed by a permutation test, 2000 permutations of fields across groups.
- Q4: J_pool > 1.2 in both regimes. J_field ∈ [0.75, 1.33] in both regimes. The Jensen gap of nested P3 persists under
  N-averaging and is again explained by per-field rates.

Descriptive only (no prediction):
- the same quantities for p̂_rat;
- the correlation of p̂_alg with the number of roots of f modulo primes ≤ 50;
- the run time.

**Interpretation, fixed in advance.** If Q1–Q4 hold, the toy is consistent with the claim that averaging over N acts
only through the field law, and leaves the field-level heterogeneity (the empirical face of nested §6) intact. If Q2
fails with a ratio < 0.67, N-averaging reduces heterogeneity in the toy. That would be a reason to re-examine Prop NB's
practical scope, though not its algebra. Toy sizes support nothing rigorous.

Deviations, if any, will be disclosed in the results note.
