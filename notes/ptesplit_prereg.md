# ptesplit pre-registration (frozen before any E1/E2 run)

Agent `ptesplit`, 2026-10-20.
- Code: `src/ptesplit/{census,nccheck,run_experiment}.py`; tests `tests/test_ptesplit.py` (5 pass; the tests ran only
  tiny cells (M, R) ≤ (6, 30) as correctness checks).
- Config: `experiments/configs/ptesplit.json`.
- Commands: `.venv/bin/python -m src.ptesplit.run_experiment E1`, then `... E2` (one job at a time; check
  `memory_pressure` first; target RSS ≤ 2 GB).

Computed before freezing (analytic, not a run): box-model pair predictions, from `census.box_model_pairs`.

## E1: exhaustive PTE fibre census

**Cells.** (M, R) = (4,60), (5,60), (6,40), (7,30), (8,24), (9,20), (10,18): all multisets of size M over {0..R}
(equivalent to [−B, B] with R = 2B by translation). Counts range from 6.4·10⁵ to 1.3·10⁷.

**Statistics.** For k = 1, 2, 3 free coefficients (s = M − k prescribed power sums), we record:
- the fibres;
- the nontrivial (size ≥ 2) fibres;
- colliding pairs;
- the maximum fibre.

For k = 1 we also record the constants |∏A − ∏A′| (since P_A − P_{A′} is constant).

**Method.** Grouping uses two independent 64-bit hashes of (p₁..p_s). The false-merge probability is ≲ 10⁻²³ per
cell; this is disclosed, and a brute-force cross-check is in the tests.

Box-model predicted pairs (T²/2∏_{i≤s}(MR^i + 1)):

| (M, R) | k = 1 | k = 2 | k = 3 |
|---|---|---|---|
| (4,60) | 0.067 | 5.8e4 | 8.4e8 |
| (5,60) | 9e−8 | 5.8 | 6.3e6 |
| (6,40) | 5e−15 | 3e−6 | 49 |
| (7,30) | 4e−23 | 2e−13 | 4e−5 |
| (8,24) | 6e−32 | 2e−21 | 3e−12 |
| (9,20) | 2e−41 | 4e−30 | 4e−20 |
| (10,18) | 3e−52 | 6e−40 | 6e−29 |

**Predictions (frozen).**
- **P1 (theorem check).** No fibre exceeds the injection bound C(R + k, k).
- **P2 (theorem check, Kleiman 1975).** Every k = 1 constant is a nonzero multiple of (M − 1)!.
- **P3 (structure beats the box model at k = 1).** Nontrivial k = 1 pairs exist for M = 4, 5, 6, although the box model
  predicts < 0.1. None exist for M ≥ 7 in these ranges, since the smallest known ideal solutions of size ≥ 7 have
  entries ≫ 30.
- **P4 (structure at k = 2, 3).**
  - In every cell with box prediction ≥ 10⁻³, observed pairs ≥ 10 × the box prediction.
  - Nontrivial k = 3 fibres exist in every cell with M ≤ 8.
- **P5 (determinant-method scale).** The k = 1 maximum fibre is ≤ 2√(R + 1) + 2 in every cell. This is a numeric proxy for
  ptebarrier's B^{1/2} bound with c = 1, κ = 0.
- **P6 (fibres are far below the trivial bound).** For M ≥ 6 and k ∈ {2, 3}, the maximum fibre is ≤ 0.1 × C(R + k, k).
- **Kill / interest signal.** Any k = 1 fibre of size ≥ R/2 would contradict the B^{1/2+o(1)} bound at these sizes and
  would be investigated as a bug or a structured family.

## E2: norm compression on exact presentations

**Cells.** F ∈ {X³ − 2, X³ + X + 1, X² − 2, X² + 1, X⁴ + X + 1} (monic and irreducible; N := F(m) is implicit and irrelevant
to an exact identity). We take (M, B) ∈ {(4, 8), (5, 8), (6, 6)} with roots in [−B, B], and h ∈ {20, 200}.

**Members.** These are the multisets with 0 < ‖P_A mod F‖_∞ ≤ h. Statistics:
- members;
- distinct residuals r;
- distinct (r, non-unit part);
- the maximum number of non-unit pieces (|F(a)| ≥ 2);
- the fraction of members with ≤ 2 non-unit pieces.

**Predictions (frozen).**
- **P7 (theorem check).** ∏_j|F(a_j)| = |Res(F, r)| ≤ ‖r‖₁^d‖F‖₂^{d−1} for every member: 0 violations.
- **P8 (no new information beyond the residual).** In every cell, #distinct (r, non-unit part) ≤ 3 × #distinct r.
- **P9 (unit domination).** For d ≥ 3 and M ≥ 5, more than half of the members have ≤ 2 non-unit pieces.
  Equivalently, the "dense split family" consists of unit pieces (|F(a)| = 1) with multiplicities.
