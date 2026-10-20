# lv5unif pre-registration (frozen before any run)

Agent `lv5unif`, phase 12, 2026-10-20. Written and committed before the code below is run.

## E1 — exact check of the dimension-free kernel bound (Lemma Kern-U)

**Claim under test (proved in notes/lv5unif_results.md, Lemma Kern-U).** Let d ≥ 3, P = (a, b), P′ = (a′, b′)
primitive with b, b′ ≥ 1 and Δ = ab′ − a′b ≠ 0, and Φ(c) = (Σ c_i a^{d−1−i}b^i, Σ c_i a′^{d−1−i}b′^i) on
C_H = [−H, H]^d ∩ Z^d. Then for every z,
  #{c ∈ C_H : Φ(c) = z} ≤ (⌊2H/(bb′)⌋ + 1)^{d−2}.

**Code (frozen with this file):** `src/lv5unif/kern.py` (`fibre_histogram_brute`, `fibre_count_lattice`, `triangular_bound`),
`src/lv5unif/e1_kern.py`. The fibre count is computed by brute force over C_H for (2H+1)^d ≤ 2·10⁶, and
otherwise by exact lattice enumeration over the kernel coordinates with *all d* box constraints checked (the
enumeration prunes on the triangular coordinates only to save time; it never uses the inequality being tested).
Honest scope: in lattice mode each triangular level has at most ⌊2H/(bb′)⌋ + 1 choices by construction, so P1 is an
*independent* test only in brute-force mode (small d, H) and through the brute-force/lattice cross-check on every
brute-force pair (which tests the kernel basis and the fibre parametrisation). P2 and P3 are informative in both modes.

**Grid.** d ∈ {3, 4, 5, 6, 7, 9}; 40 pairs (P, P′) per d drawn with a fixed seed (seed 20261020) from
V/2 < b, b′ ≤ V (the candidate range), |a|, |a′| ≤ V, gcd(a, b) = gcd(a′, b′) = 1, Δ ≠ 0, with V ∈ {3, 4, 6} and
H ∈ {V², 2V², dV²}. A drawn pair whose triangular bound exceeds 3·10⁶ is redrawn (the enumeration visits at most
that many leaves), and the number of redraws is reported; so large d is tested only at small H/(bb′).
For each pair: the central fibre z = 0 and 8 fibres z = Φ(c) for c uniform on C_H (seeded).

**Predictions.**
- P1 (correctness): 0 violations of the bound over all tested fibres. Any violation refutes the lemma as stated.
- P2 (the bound is not wildly loose where it matters): for pairs with b, b′ > V(1 − 1/d) and H = dV², the
  median over pairs of (central fibre count)/(triangular bound) is ≥ 0.05 for every d tested.
- P3 (dimension-free constant): over all pairs with H ≥ dV², the maximum of
  (2H + 1)^{2}(bb′)^{d−2}·max_z Pr[Φ(c) = z] is ≤ e^{1/2}·4 ≈ 6.6 (this is P1 rewritten with the constant
  of the lemma; reported separately to make the d-independence visible).

**Reporting.** Counts, the P1–P3 verdicts, and the per-d maxima/medians, in `experiments/raw/lv5unif_e1.jsonl` and
`experiments/raw/lv5unif_e1_summary.json`. E1 is a correctness check of an exact combinatorial lemma; it cannot test
any asymptotic claim.

## No other experiments

The τ-table, the parameter/margin checks and the cost constant are exact symbolic computations
(`src/lv5unif/params.py`), not prediction tests.
