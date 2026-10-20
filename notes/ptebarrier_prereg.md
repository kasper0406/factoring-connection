# ptebarrier pre-registration (agent `ptebarrier`, phase 16, 2026-10-20)

Written and committed before any run of the experiments below. The runs are exact checks of a
proved lemma and exploratory counts; they are not evidence for the lemma's proof, only a bug
detector for it.

## The statement under test (Lemma R, frozen here)

Let P, P' ∈ Z[X] be monic of degree M, split over Z, all roots in [−B, B] (B ≥ 1).
Let Q, Q' ∈ Z[X] \ {0} have degree ≤ j. Put G := Q·P − Q'·P'. Suppose deg G ≤ K − 1 and
‖G‖_∞ ≤ h. Put Δ := 2^{K−1}·√K·h·B^{K−1} and L := ⌊log₂ Δ⌋.

- (R0, j = 0, Q = Q' = 1) If M ≥ K + 2L + 3, then P = P'.
- (R, general j) If M ≥ K + 2L + 4j + 3, then QP = Q'P'.

Define h_R(B, M, K, j) := max{h ∈ Z_{≥1} : the hypothesis on M holds}, or 0 if there is none.

## E1: exhaustive check of R0

- Grid: roots in [−B, B]. (B, M_max) ∈ {(1, 30), (2, 22), (3, 16), (4, 13), (5, 11)}. M runs from 2 to M_max, and K ∈ {1, 2, 3, 4} with K < M.
- For each (B, M, K), enumerate all multisets of size M and group them by the coefficients of P in degrees ≥ K. Let h_min(B, M, K) be the minimum of ‖P − P'‖_∞ over pairs of distinct multisets in one group (∞ if no group has two members).
- **P1 (theorem check; pass/fail).** h_min(B, M, K) > h_R(B, M, K, 0) in every cell. Any violation means the proof is wrong and will be reported as such.
- **P2 (exploratory; not a test of the lemma).** We expect the bound to be far from tight: at the largest M of each (B, K) row, h_min should exceed h_R by at least a factor of 10 wherever h_R ≥ 1. We expect log₂ h_min to grow at least linearly in M; the lemma guarantees slope ≥ 1/2 asymptotically.

## E1b: exhaustive check of R with j = 1

- Grid: (B, M_max) ∈ {(1, 24), (2, 16), (3, 12)}. Q ranges over the 24 polynomials c₁X + c₀ with c₀, c₁ ∈ [−2, 2], not both 0. K ∈ {1, 2, 3}.
- For each (B, M, K), consider all pairs ((A, Q), (A', Q')) with QP_A ≠ Q'P_{A'} and deg(QP_A − Q'P_{A'}) ≤ K − 1. Let h_min be the minimum height of the difference.
- **P1b (pass/fail).** h_min > h_R(B, M, K, 1) in every cell.

## E2: fiber statistics (exploratory)

- Grid: B ∈ {2, 3, 4} with the M ranges of E1. s = M − k for k ∈ {1, 2, 3, 4}. Record the maximum and the mean size of nonempty fibers of (e₁, …, e_s), and the number of fibers of size ≥ 2.
- **P3 (sanity; pass/fail).** max fiber ≤ C(q + k − 1, k), where q = 2B + 1 (the injection bound).
- **P4 (sanity; pass/fail).** For k = 1 and every fiber, the members are pairwise disjoint as multisets.
- **P5 (exploratory).** Under the heuristic local-limit count, fibers with s ≥ √(2M) + 2 should be generically trivial. We expect nontrivial fibers to persist there for small M (structured PTE solutions), and we will report how many.

## E3: the Dickson example (exploratory)

- F = X⁶ − 6X⁴ + 9X² − X − 7. Count (a, b) ∈ [−200, 200]² with F + aX + b split over Z.
- Expectation: only the members (X² − 4)(X² − 1)², i.e. (a, b) = (1, 3), plus a handful of degenerate ones.

## Reproduction

`.venv/bin/python -m src.ptebarrier.run_experiment all`. Raw results go to `experiments/raw/ptebarrier_E*.jsonl` and the config to `experiments/configs/ptebarrier.json`. Memory stays below 2 GB; one job at a time.
