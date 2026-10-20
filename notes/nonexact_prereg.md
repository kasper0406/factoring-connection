# nonexact: pre-registration (phase 17, O1)

Agent `nonexact`, 2026-10-20. Committed BEFORE any experiment below is run. Before this commit, only unit tests on
tiny instances (`tests/test_nonexact.py`), a count of fibre pairs per E3 cell (below) and a formula evaluation of the cost model on the E1 grid (no
lattice, no hit counting) were run. Config: `experiments/configs/nonexact.json`. Code: `src/nonexact/`.
Raw output: `experiments/raw/nonexact_E{1,1b,2,3}.jsonl`.

## What is being tested

The theory side (to be written up in `notes/nonexact_results.md`) argues as follows.
- Non-exact and evaluated compressions are **abundant**: there is no counting barrier.
- Their only obstruction is the cost of *finding* them.
- Lattice reduction imposes the congruence mod N linearly, but loses the split condition, which is nonlinear and
  thin.

The experiments measure the three lattice mechanisms the task names.

**E1, linearised lattice (Prop LIN).**
- *Setup.* View B with d = 1: μ is uniform mod N.
  - A random prefix multiset C of size M − s gives c = P_C(μ) mod N.
  - The unknown suffix S (an s-multiset in [−B, B]) is linearised through its elementary symmetric functions:
    c·P_S(μ) ≡ w₀ + Σ e_i w_i (mod N).
  - The member condition |c·P_S(μ) mod± N| ≤ V is then one linear congruence in (e₁..e_s, v).
  - The box is |e_i| ≤ C(s,i)B^i, |v| ≤ V.
- *Model (HEURISTIC: equidistribution of the coset).*
  - Box points per lattice: p = Vol_e·(2V+1)/N, with Vol_e = ∏(2C(s,i)B^i + 1).
  - Split hits per lattice: h = C(2B+s, s)·(2V+1)/N.
  - Cost per hit: max(1, p)/h (at least one unit per lattice setup, one per enumerated point).
  - Brute force: N/(2V+1) per hit.
- *Measured.*
  - Hits are counted over T lattices, by exact evaluation of all s-multisets.
  - Box points are counted exactly on a subset.
  - The LLL + exact Fincke–Pohst enumerator (`lin_lattice_enumerate`) runs on a subset where p ≤ 200, and is
    checked against the exhaustive box count.

**E2, multiplier lattices (Prop LAT).**
- *Setup.*
  - F = X − μ (d = 1), or F = (X − μ)(X − μ′) lifted mod N (d = 2), with μ, μ′ uniform.
  - P is a random split polynomial with M = 12 roots in [−16, 16].
  - L_P = {(Q, r) : deg Q ≤ j, deg r < k, QP ≡ r mod (N, F)}, of dimension n = j + 1 + k and determinant N^d.
  - L_P is reduced exactly with flint LLL.
- *Measured.*
  - ‖b₁‖ against the Gaussian heuristic GH = (N^d·Γ(n/2+1))^{1/n}/√π. Floats are used only for this reported
    ratio.
  - The residual height log₂‖r‖_∞/log₂ N of b₁.
  - The number of "useful" basis vectors: Q and r both nonzero, and both split over Z (constant × ∏(X − a),
    |a| ≤ 16).
  - Integer-root rates of Q, against a baseline of random polynomials with the same coordinate-wise heights.

**E3, inheritance (Prop IH, k ≥ 1).**
- *Setup.*
  - Fibre pairs: P ≠ P′ monic split, roots in [−R, R], with deg(P − P′) = k − 1 exactly (shared e₁..e_{M−k}).
  - For each pair, primes N > max(4H+1, 50) dividing the coefficients of P − P′ − δ for some δ with
    ‖δ‖ ≤ 2H, H = 3.
  - A monic F of degree d = max(k, 2) dividing P − r mod N, for some r with ‖r‖, ‖r − δ‖ ≤ H. This is built
    by factoring mod the prime N (DIAGNOSTIC harness; N prime is allowed because Prop IH does not need N
    composite).
  - Both P ≡ r and P′ ≡ r − δ (mod N, F) are verified.
- *Measured.* The integer relations of Prop IH at the non-common roots: u_b = ±D(b) ≡ v_b = ±δ(b) (mod N).

## Frozen predictions

**E1.** Cells: n ∈ {24, 32, 40} bits (balanced semiprimes, `balanced_semiprime(n, 1)`), B ∈ {8, 16},
s ∈ {1, 2, 3}, V = 255 (n ≤ 32) or 65535 (n = 40). T = ⌈80/h⌉ lattices; a cell needing T > 10⁶ is not run and is
reported as a prediction only. M = ⌈n·ln 2/ln(2B)⌉ + 2.
- **P1 (hit rate).** Observed/predicted total hits ∈ [0.6, 1.5] in every run cell. With ≥ 200 predicted hits the
  band is [0.8, 1.25].
- **P2 (box points).** On the box subset, Σ observed p / Σ predicted p ∈ [0.9, 1.1] when Σ predicted ≥ 1000, and
  ∈ [0.5, 2] when it lies between 30 and 1000.
- **P3 (cost per hit).** Measured (T + T·p̄)/hits ∈ [0.6, 1.5] × the model's max(1, p)/h. Here p̄ is the subset
  mean, or the model p where p < 1. The best s in each (n, B) row is the model's arg-min.
- **P4 (no excess over brute force beyond the suffix count).** In every run cell, brute-force cost per hit
  divided by measured cost per hit ≤ 1.5·C(2B+s, s). In words: the lattice saves at most the number of
  suffixes, a subexponential factor B^{s}/s!.
- **P5 (algorithm correctness and cost).** On the enumeration subset the LLL + Fincke–Pohst point sets equal the
  exhaustive box sets: 0 mismatches. Median enumeration nodes ≤ 30·(1 + p).
- **Kill signal (construction side).** Any run cell with observed hits > 3× predicted (and predicted ≥ 20)
  means structured excess. That would be a positive signal for O1 and would be investigated.

**E2.** n ∈ {32, 48, 64}; d = 1 with (j, k) ∈ {(2,1), (4,1), (8,1)}; d = 2 with (j, k) ∈ {(1,2), (3,2), (6,2),
(10,2)}; T = 100 split P per cell.
- **P6 (Minkowski / GH).** The median of ‖b₁‖₂/GH is in [0.6, 1.4] in every cell.
- **P7 (residual size).** The median of log₂‖r(b₁)‖_∞/log₂N is in [d/n − 0.12, d/n + 0.12], where n is the
  lattice dimension. If r(b₁) = 0, use the smallest nonzero r among the basis vectors.
- **P8 (no useful vectors).** The total number of useful vectors with deg Q ≥ 1 is ≤ 2 per cell. Their
  integer-root rate (Q has ≥ 1 integer root) is within a factor 3 of the random baseline when the baseline count
  is ≥ 10; otherwise both are ≤ 10 per cell.
- **Kill signal.** More than 5 useful vectors with deg Q ≥ 1 in any cell.

**E3.** Cells (k, M, R) = (1,4,8), (1,5,10), (1,6,11), (2,4,8), (2,5,8), (2,6,6), (3,5,8), (3,6,6). *Disclosed:* this
list was fixed after counting how many fibre pairs each candidate cell contains (33, 4, 8, 3256, 593, 6, 49308,
1116; the k = 1 cells (1,5,8) and (1,6,6) have none). Only that count was computed, no N, F or relation.

*Also disclosed:* `test_runner_smoke` runs E1–E3 end to end on a tiny non-registered config (n = 20 and 24 bits,
B ≤ 6, 3 lattices per E2 cell, 1–2 pairs per E3 cell), writing to a temp directory. Its printed output was seen:
E1 hit ratios 1.0 and 2.2 (5 and 11 hits), E2 ‖b₁‖/GH ≈ 1.27 at dimension 4, E3 0 violations. None of it
changed a prediction above. Up to 200 pairs per cell; H = 3.
- **P9 (Prop IH).** 0 violations of the following, for every constructed pair:
  - (a) u_b ≡ v_b (mod N);
  - (b) nontrivial ⇔ E(b) ≠ 0;
  - (c) |v_b| ≤ 2H·Σ_{l<k}R^l;
  - (d) every piece |b − a| ≤ 2R;
  - (e) the number of distinct nontrivial relations is ≥ the Prop IH lower bound
    (⌈(n_supp − (k−1))/deg D⌉, or 1 if deg D = 0).
- **P10.** For k = 1, every pair yields exactly one distinct nontrivial relation (u_b = ±C for all b).
- **Measured, not predicted.** The mean number of distinct nontrivial relations per pair for k = 2, 3.

## Resources

- One process at a time, after checking `memory_pressure`.
- Target ≤ 2 GB RSS, measured with `/usr/bin/time -l`.
- Expected run time: E1 ≈ 10–20 min, E1b ≈ 2 min, E2 ≈ 2 min, E3 ≈ 1 min.
