# rank3 pre-registration (phase 17, 2026-10-20)

Agent `rank3`. Written and committed **before** any experiment below is run. Experiments are exact (integers,
python-flint `arb`/`acb` balls only to certify integer values, never as results). Memory target ≤ 2 GB RSS, one
heavy job at a time.

## Claims under test (stated in advance; proofs go in notes/rank3_results.md)

Setting: orders O ⊂ K, O′ ⊂ K′ (degrees d, d′, D = d + d′), ring maps φ, φ′ to Z/N, a linked rank-r module
w_j = (u_j, u′_j), φ(u_j) = φ′(u′_j), u_j, u′_j algebraic integers; ℓ_σ = (σ(u_j))_j, ℓ′_τ = (τ(u′_j))_j;
F = W·W′ = ∏ℓ_σ ∏ℓ′_τ; Δ_F = ∏ over ordered r-tuples of distinct factors of det (Bean–Thunder's invariant).

- **T1 (divisibility).** Δ_F ∈ Z, and v_p(Δ_F) ≥ P(r, D) := r(r−1)(D−2)!/(D−r)! for every p | N unramified
  in K and K′.
- **T2 (volume).** With Bean–Thunder (JLMS 54, 1996), V_F ≤ C_r |Δ_F|^{−(D−r)!/D!} ≤ C_r N^{−r(r−1)/(D(D−1))}.
- **T3 (box count).** #{a ∈ Z^r, |a_j| ≤ E : 0 < |F(a)| ≤ X} ≤ C_r|Δ_F|^{−(D−r)!/D!} X^{r/D} + r(D+1)(2E+1)^{r−1}.
- **T4 (rank-r Liouville).** 𝓜 := ∏_{i ≤ D} max_j |L_i(e_j)| ≥ (2^r/(r!·C_r))^{D/r}·|Δ_F|^{(D−r)!/(r·(D−1)!)}
  (≥ N^{(r−1)/(D−1)}·(2^r/(r!C_r))^{D/r} under T1).

## E1: theorem checks (exact)

Instances (seeded, `src/rank3/e1.py`, seed 20261020):
- (A) "GNFS-like two-algebraic": N = pq with p, q random primes of 10 and 11 bits (sympy.randprime, seeded; primality
  test used), m random mod N, f, g monic of degree d with root m mod N from small random top coefficients and a
  reduced constant term (as in a lattice/base-m construction). Power basis, r ∈ {3, 4} with r ≤ min(d, d′),
  (d, d′) ∈ {(3,3), (4,3), (4,4)}; 3 seeds each.
- (B) "SNFS-like": small random monic irreducible f, g (coefficients in [−3, 3]), N := the largest squarefree
  divisor of |Res(f, g)| coprime to disc(f)·disc(g) and to all primes ≤ D; skip if N < 1000. Power basis, same
  (d, d′, r) grid, 3 seeds.
- (C) non-power-basis linked modules: as (A) but u_j = P_j(α), u′_j = Q_j(β) with random small P_j and
  Q_j = P_j + s_j(x − m̃) (m̃ the centred lift of m), r = 3, (d, d′) ∈ {(3,3), (4,4)}, 3 seeds.
- Box counts for T3: E ∈ {3, 6} (r = 3), E = 3 (r = 4); X on a dyadic grid from min|F| to max|F| over the box.

Predictions (theorem checks: any violation means a bug or a false theorem and is reported as such):
- **P1.** T1 holds in every instance (0 violations), with Δ_F certified to be a nonzero integer.
- **P2.** T3 holds at every (instance, X) (0 violations), using C_r from the Beta-function formula, computed as a
  certified `arb` upper bound.
- **P3.** T4 holds in every instance (0 violations).
- **P4 (sharpness, may fail).** The main-term saturation s(X) = count(X)/(C_r|Δ_F|^{−(D−r)!/D!}X^{r/D}),
  maximised over grid points X where the main term exceeds the error term, satisfies s ≥ 10^{−3} in at least
  half of the (B) instances. (This says the volume bound with the exact Δ_F is not vacuous by many orders.)
- **P5 (GNFS vs SNFS, may fail).** Median over seeds of max_X s(X) is larger for (B) than for (A) in every
  (d, d′, r) cell where both have eligible X.

## E2: construction comparison at toy sizes (measured; NOT a test of asymptotics)

Sizes N = pq, n = log₂N ∈ {40, 50, 60}, 3 seeds each (p, q balanced, sympy.randprime seeded).
Designs (all candidates in a box, total T ≈ 2^17 points, zero and non-primitive candidates discarded):
- R2-rat: GNFS pair f (degree 3, base m = ⌈N^{1/4}⌉), g = x − m; candidates a − bα, |a| ≤ E, 1 ≤ b ≤ E.
- R2-nn: two nonlinear sides: f, g of degree 3 with common root m mod N from a reduced lattice basis of
  {h ∈ Z[x]_{≤3} : h(m) ≡ 0 mod N} (two shortest basis vectors that are irreducible and coprime; flint LLL);
  candidates a − bθ on both sides.
- R3-rat: power basis rank 3 on the R2-rat pair: a + bα + cα² ↦ (A(α), A(m)), box |a|,|b| ≤ E, 1 ≤ c ≤ E.
- R3-nn: power basis rank 3 on the R2-nn pair.
- R3-lll: rank-3 module spanned by the three shortest vectors of the linked lattice Λ = {(h(α), h(β)) : h ∈ Z[x]}
  (degree < 6) in a weighted Minkowski norm, found by integer LLL on a scaled embedding (scaling chosen to equalise
  the two sides' typical sizes; the choice is a heuristic search step, counts are exact).
Measure exactly: median log₂|F| over the box, and the number of candidates with both W and W′ B-smooth,
B = 2^12 (by repeated gcd with the primorial), per design.

Predictions:
- **P6.** In every (n, seed), yield(best rank-3 design) ≤ 2 × yield(best rank-2 design).
- **P7.** In every (n, seed), median log₂|F| of every rank-3 design ≥ median log₂|F| of the best rank-2 design
  minus 2 bits.
- **Kill/lead signal.** If at n = 60 the best rank-3 yield exceeds 2× the best rank-2 yield in ≥ 2 of 3 seeds and
  the ratio increases from n = 40 to n = 60, I flag it as a construction lead and analyse it (this would not by
  itself contradict an asymptotic theorem).

Resources: E1 a few minutes; E2 ≈ 15 × 2^17 resultant evaluations per design, ≤ 30 min, < 1 GB.
Raw data: `experiments/raw/rank3_E1.jsonl`, `experiments/raw/rank3_E2.jsonl`.
