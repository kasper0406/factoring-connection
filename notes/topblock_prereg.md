# topblock pre-registration (phase 18)

Agent `topblock`, 2026-10-20. Written and committed **before** any experiment below was run.
Code: `src/topblock/{core,lattices,run_experiments}.py`; tests `tests/test_topblock.py` (7 pass).
Command: `PYTHONPATH=. .venv/bin/python -m src.topblock.run_experiments <x>`; raw data
`experiments/raw/topblock_<x>.jsonl`; summary `experiments/analysis/topblock_summary.json`.
Instances: `balanced_semiprime(bits, seed)` from `src/common/instances.py`. Environment: fpylll 0.6.x
needed `cysignals`, which was installed into `.venv` (`uv pip install cysignals`).

## Hypotheses under test

- **H-EQ.** The ladder's top-block window sets are exactly near-divisor sets (identity (W)), and
  window counts reduce to divisor search with a provable call bound (reduction (iii)); a general
  divisors-in-interval finder applied to N − k answers window counts (reduction (iv)).
- **H-MIX (the main question).** Adding window congruences mod the known D = x^{j+1} (one window,
  two windows on adjacent blocks, or two Taylor degrees) to the Coppersmith divisor lattice does
  **not** extend its range beyond what a univariate lattice in (N, x+h)^M achieves. Theory (notes):
  the substitution r ↦ r(h) maps every mixed vector into D^{m2}(N, x+h)^{M} with M = m1 + T and degree
  ≤ E + T; vectors in its kernel are multiples of the curve relation F = (x+h)r + N(−h)^{j+1} and locate
  nothing. F·N^{m1−1}D^{m2−1} lies in the mixed lattice and is Howgrave-Graham-short up to
  X ≈ x·q^{−1/(j+1)}, so naive "short vector exists" counts overstate the range.
- **H-HM.** Splitting the unknown as (block index u, offset v), Herrmann–May style, gives no range
  beyond the univariate lattice; its short vectors are functions of z = Hu + v.
- **H-SVP.** At fixed dimension, exact SVP extends the univariate range by at most a bit over LLL.

## Scan protocol (X3–X6)

For each lattice, log₂X runs over 4.0, 4.5, 5.0, … For each point the planted offset h₀ is drawn
uniformly from [X/2, X) with a seed derived from (bits, seed, point), the same for every lattice
at that point. The block start is x = p − h₀; lattices see only (N, x, parameters). **X_max** is the
last success before the first two consecutive failures. Success means a divisor of N was recovered
from the reduced basis: integer roots of a Howgrave-Graham-good vector (U); of a good vector's
nontrivial univariate image G ≠ 0 (MIX); or from linear factors/resultants (HM). HG-good means
n‖v‖² < (x^{m1}∏D_b^{m2_b})², using the public bound p ≥ x.

## Experiments and predictions

**X1 (identity (W)).** 20/24/28-bit N, seeds 0–9, j ∈ {2,3,4}. For each, the top block containing p
plus 3 random top blocks, with the divisor window and 10 random rational windows; window set by
direct evaluation vs. by the near-divisor criterion.
- P1: 0 mismatches.

**X2 (reductions).** Reduction (iii): bisection with exact-divisor windows, brute-force oracle, on the
top block containing p; 40 bits (j = 2,3,4), 48 bits (j = 2,3), 44 bits (j = 4), seeds 0–9.
Reduction (iv): 32 bits, seeds 0–9, j = 2,3, 20 random (subinterval, window) queries each.
- P2.1: p found in 100% of blocks.
- P2.2: calls ≤ the proved bound C^w (notes §1.3) in 100% of blocks.
- P2.3: median calls/⌈log₂H⌉ ≤ 3.0.
- P2.4: (iv) counts equal oracle counts in 100% of queries.

**X3 (single-window mixed lattice).** 64/96/128 bits, seeds 0–2, j ∈ {2,3}. Configs (m1, m2, E):
j = 2: (1,1,8), (2,1,10), (2,2,10); j = 3: (1,1,9), (2,1,11), (2,2,12). That is 54 cells.
Comparators: U_img = U(M = m1+m2, n = E+m2+1), the transfer target; U_dim = best X_max over
U(m,n) ∈ {(3,8),(4,10),(5,12),(6,14),(7,16),(9,20),(11,24),(13,28),(15,32)} with n ≤ dim(MIX).
For m1 = 1 the scan continues in 1-bit steps up to log₂X = 0.45·bits, recording HG-good counts.
- P3.1: log₂X_max(MIX) ≤ log₂X_max(U_img) + 3 in every cell.
- P3.2: log₂X_max(MIX) ≤ log₂X_max(U_dim) + 3 in every cell.
- P3.3: in ≥ 80% of the 18 m1 = 1 cells, MIX has HG-good vectors at some log₂X ≥ log₂X_max(U_img) + 4.
  The total number of nontrivial (G ≠ 0) HG-good vectors at log₂X > log₂X_max(U_img) + 3 is 0.
- P3.4: for m1 = 1 cells, (largest scanned log₂X with an HG-good vector) − (largest scanned log₂X
  at which the F row is HG-short) ∈ [−3, +3].

**X4 (two windows).** 64/96 bits, seeds 0–2. Windows: adjacent blocks [(Δ = 0, j = 2), (Δ = X, j = 2)]
with (m1, E) ∈ {(1,7), (2,9)}; two degrees [(0, j = 2), (0, j = 3)] with (m1, E) ∈ {(1,8), (2,10)};
m2 = 1 each. That is 24 cells. U_img = U(m1+2, E+3).
- P4.1: log₂X_max(MIX2) ≤ log₂X_max(U_img) + 3 in every cell.
- P4.2: log₂X_max(MIX2) ≤ log₂X_max(U_dim) + 3 in every cell.

**X5 (HM split).** 64/96 bits, seeds 0–2, (m, E) ∈ {(2,4), (3,6)} (dimensions 15, 28);
H_s = 2^{⌊log₂X/2⌋}, U = ⌈X/H_s⌉. That is 12 cells.
- P5.1: log₂X_max(HM) ≤ log₂X_max(U_dim) + 1 in every cell.
- P5.2: ≥ 80% of all HG-good HM vectors (over all scans) have a linear factor b(H_s u + v) + c.

**X6 (exact SVP vs LLL, univariate).** 64/96 bits, seeds 0–2, (m, n) ∈ {(3,8),(4,10),(5,12),(6,14),(7,16)};
exact SVP by unpruned enumeration. That is 30 cells.
- P6.1: 0 ≤ log₂X_max(SVP-HG) − log₂X_max(LLL) ≤ 1.0 in every cell.
- Descriptive only (no prediction): log₂(λ₁/GH) beyond the LLL threshold.

## Decision rules

- **H-MIX rejected as a mechanism** (windows add no range) if P3.1, P3.2, P4.1 and P4.2 all hold. If
  any cell exceeds its comparator by more than 3 bits, the gap is re-measured at the next bit size up,
  and the result is reported as a candidate mechanism only if the gap grows with log N.
- **H-HM** supported if P5.1 and P5.2 hold. **H-SVP** supported if P6.1 holds.
- Toy successes are correctness checks only (AGENT_RULES §8).

## Disclosed development checks (before this commit)

All at sizes or seeds outside the grid:
- 48-bit seed 999 smoke runs of U(3,7), U(4,12), MIX(j = 2, m1 = 1, m2 = 1, E = 6), MIX2 adjacent
  (m1 = 1, E = 7), HM (2,4)/(3,6) and SVP n = 12.
- 96-bit seed 12345 SVP timing at n = 12, 14, 16 (three X points).
- 128-bit seed 12345 timing.

What the peeks showed:
- In the MIX smoke run, HG-good vectors persisted to log₂X = 15, all trivial, while U stopped at 10.
- In the HM run, all good vectors had the z-factor (which motivated the linear-factor recovery).
- λ₁ was far below GH beyond the threshold, so no λ₁/GH prediction is made.
- Exact SVP at n = 20 and 128 bits exceeded 9 GB RSS and was killed, so X6 is capped at n ≤ 16 and
  96 bits.
