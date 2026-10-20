# rankdegen pre-registration (agent rankdegen, phase 18, 2026-10-20)

Committed before E1–E3 are run. Code: `src/rankdegen/{fh,ledger,e1,e2,e3}.py`. Raw output:
`experiments/raw/rankdegen_E{1,2,3}.jsonl`. Memory budget ≤ 2 GB RSS, one job at a time.

**Disclosed prior computation (development, before this file).**
- The Φ-program values below were computed numerically while developing the theory (`ledger.phi_numeric`):
  Φ ≈ (r−1)/r for every type listed in P5.
- The frame structure of codex's α⁶ = 2 example (`codex/rankdegen_1.out`, eq. (2)) was computed once at N = 31
  to test `fh.py`: 40 of 220 unordered 3-minors are zero, 225 frames, LP max–min coverage 1/9 under uniform marginals.
- No box counts, yields or large-box counts have been computed.

## E1: Theorem FH and the linking exponent, exact (r = 3)

**Instances.**
- E1-C2 (codex rich-line type): f = g = x⁶ − 2, P = (1, x, x³), Q = (1, x + N, x³), with
  N ∈ {31, 31·41, 31·41·47} and the common root μ fixed by CRT from 2 (mod 31), 13 (mod 41), 16 (mod 47).
- E1-B (subfield type): f = g = x⁶ + x² + 1 (= h(x²), h = x³ + x + 1), P = (1, x², x⁴), Q = (1, x² + N, x⁴),
  N = h(14) = 2759 = 31·89 (subfield map x² ↦ 14). Only the subfield order is linked (Lemma SUB).
- E1-A (coincident type): f = g = Φ₅, P = Q = (1, x, x²), N = 31·41.
- E1-GP, E1-C1, E1-FP (seeds 0, 1, 2): N = pq with p, q seeded 10/11-bit primes (`rank3.e1.inst_A` generator,
  d = d' = 3), common root m.
  - GP: P = Q = (1, x, x²).
  - C1 (rank-deficient K side, ρ = 2): P = (1, x, 0), Q = (1, x, x² − m̃x), m̃ the balanced residue of m.
  - FP (both sides rank 2, fibre product): P = (1, x, 0), Q = (1, x, x − m̃).
- Box [−E, E]³ with E = 30, all instances.

**Distributions.** For each instance: pure-s distributions (uniform over frames with s K-rows, s = 0..4, when
nonempty), and the exact-LP max–min coverage distribution under uniform marginals (r+1)/D (when feasible).

**Checks.** For each distribution w with side weights (y_K, y_K′) and min cross-pair coverage t:
- the N-only bound log J ≤ log(2^r C_r) − ((r−1)/(r+1))·t·log N;
- for pure-s w, also the actual-minor bound log J ≤ log(2^r C_r) − Σ_S w_S log Φ_S/(r+1);
- count{a ∈ box : 0 < |W| ≤ Y₁, 0 < |W′| ≤ Y₂} ≤ J·Y₁^{y_K}Y₂^{y_K′} + r(d+d′)(2E+1)^{r−1}, on the grid
  Y_i = 2^k (all k from the smallest to the largest observed bit length).

**Predictions.**
- **P1 (FH):** 0 violations of the count check, for both bounds, over all instances, distributions and grid points.
  Also report how many checks are informative (count > error term).
- **P2 (linking valuation):** for every pure-s w with t > 0, certified Σ_S w_S log Φ_S ≥ (r−1)·t·log N.
  0 violations.
- **P3 (type A):** for E1-A, the LP max–min cross-pair coverage is exactly 0 (coincident cross pairs lie in no frame).
- **P4 (C2 LP):** for E1-C2, the LP value is 1/9 at every N (structure-determined).
- **P5 (Φ program):** `ledger.phi_numeric` returns (r−1)/r within 10⁻³ for: GP two-sided r = 2, 3, 4; C1 r = 3
  (ρ, ρ′) = (2, 3), (2, 2); r = 4 (3, 4), (3, 3), (2, 4), (2, 3).

## E2: degenerate-family yields vs rank 2 (toy sizes; measured, not asymptotic)

Setup as `rank3.e2` (reused by import): n = log₂N ∈ {40, 50, 60}, 3 seeds each (same seed strings), B = 2^12,
same rank-2 designs R2-rat and R2-nn and boxes (rank 2: |a| ≤ 256, 1 ≤ b ≤ 256; rank 3: |a|, |b| ≤ 32,
1 ≤ c ≤ 32; primitive candidates only). Degenerate rank-3 designs on the lattice pair (h₁, h₂):
- D-C1: P = (1, −x, 0) on h₁, Q = (1, −x, ν) on h₂, with ν the shortest (LLL) quadratic h with h(m) ≡ 0 (mod N)
  that is not a multiple of a linear one.
- D-FP: P = (1, −x, 0), Q = (1, −x, x − m̃).
- D-RK (rational side + rank-deficient K side): f_rat, g_rat (base m), P = (1, −x, 0), Q = (1, −x, N).
Yield = number of candidates with both norms B-smooth.

**Predictions.**
- **P6:** best degenerate yield ≤ 2 × best rank-2 yield, in 9/9 instances.
- **P7:** every degenerate design's median bit length of |F| ≥ (best rank-2 design's median) − 2, in 9/9.
- A lead signal would be any violation of P6.

## E3: lattice points in excess of the volume bound, by exceptional type (exact)

Instances: the nine (B) rank-3 (d, d′) = (3, 3) instances of `experiments/raw/rank3_E1.jsonl` (SNFS-like,
power basis), plus the three (A) instances of that cell. Box [−40, 40]³. main(X) = C₃|Δ_F|^{−1/120}X^{1/2}
(R1 with the exact Δ_F, `rank3.forms.count_bound_terms`). Points are classified as
- L: c = 0 (linear candidates, the rank-2 NFS sublattice);
- Rd: c ≠ 0 and a + bx + cx² has a rational root (b² − 4ac a perfect square; products of linear candidates);
- G: the rest (irreducible quadratic candidates).

Counts T(X) = L + Rd + G with 0 < |F| ≤ X on X = 2^k.

**Predictions.**
- **P8 (excess exists):** in every instance, some grid X has T(X) > main(X).
- **P9 (LP with exceptions, strongest form):** in every instance and at every grid X, G(X) ≤ 2·main(X) + 20.
- Reported without prediction: L/main, Rd/main, G/main as functions of X, and the largest X at which
  T(X) > main(X).
