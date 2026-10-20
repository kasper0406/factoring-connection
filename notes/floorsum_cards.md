# floorsum: candidate cards (BRIEF §4 template)

Agent `floorsum`. Verdicts point to `notes/floorsum_results.md`.

---
**F1. Farey-dissection exact oracle (Algorithm O).**
- **Mechanism.** Split each ladder block into windows, one per continued-fraction convergent p/q of the
  local slope; on each window g = qP − pu has O(1) integer crossings, and between crossings ⌊P⌋ is a
  linear floor function (Lemma R) summed by Euclid.
- **Target.** Exact Σ⌊P + s⌋ on subintervals, any degree, P″ of fixed sign.
- **Representation.** Integer numerator polynomial over x^{j+1}, from N and x.
- **Extraction.** Run decomposition with prefix sums; bisection of candidate counts as in Theorem L_j.
- **Why it exposes a factor.** Through Theorem L_j (divisors are candidates).
- **Known lemmas.** Continued fractions; Euclid floor sums.
- **New lemma.** Lemma R (one line) and Theorem O (cost).
- **Correctness.** PROVED; 0 mismatches in E2.
- **Cost.** O((1 + Hβ^{1/3}) log H) per block; N^{1/3} polylog in total for every j.
- **Randomness.** Deterministic.
- **Closest prior.** Hull walks (Sladkey; Har-Peled; Alcántara et al., whose Cor 1.4 already factors
  in O(N^{1/3} log² N)); codex's fixed-cell variant.
  Difference: generic degree, rational shifts and subinterval queries for the ladder.
- **Smallest experiment.** E2/E2b/E3, done.
- **Refutation criterion.** Cost slope ≠ 2(j−1)/(3j) on generic blocks. Not seen (0.337/0.442/0.506).
- **Repair.** Not needed.
- **Verdict.** Works, and is exactly the trap: N^{1/3} for every j.

---
**F2. Fourier smoothing + exact correction.**
- **Mechanism.** ⌊y⌋ = y − ½ − ψ(y); replace ψ by a degree-K trigonometric approximant (Vaaler); each
  frequency is an incomplete Weyl sum, computed in polylog time by Hiary for j = 2 (and for j = 3 only
  at low frequencies; nothing for j ≥ 4); correct exactly the h with P(h) near an integer.
- **Target.** Same as F1.
- **New lemma required.** K polylog with few corrections.
- **Correctness/cost.** Prop F (PROVED): the approximant is off by ≥ ¼ on measure ≥ 1/(51K), so
  ≥ 1/(51Kg) − O(K) corrections, cost ≥ 0.099/√g − 3 ≍ √H for ladder blocks (E4).
- **Closest prior.** Vaaler; Hiary's theta algorithm; barrier F of analytic.
- **Refutation criterion.** Near-integer counts ≪ δH on ladder blocks. Not seen (ratio 0.998–1.003).
- **Verdict.** REFUTED for the certified model as a route below N^{1/4}: the total is
  ≳ N^{1/4 + 1/(2(j+2))}, experimentally supported for j = 2 (it needs a small max gap g). For j ≥ 5
  that is below N^{1/3}; there the blockers are the missing Weyl-sum algorithms and the cost of
  locating E.

---
**F3. Quadratic Euclid / reciprocity.**
- **Mechanism.** Mimic the linear Euclid recursion: reduce coefficients, swap axes, recurse.
- **Obstruction.** The swap produces square-root floors, and the double swap is the identity (PROVED).
  Lattice moves preserve affine arc length (Prop A), so no sequence of them straightens the arc.
  Unrestricted quadratic floor sums are NP-hard (Manders–Adleman).
- **Verdict.** REFUTED for the equi-affine-move model; no reciprocity law found for incomplete sums.

---
**F4. Lattice (Coppersmith / Howgrave-Graham) block oracle.**
- **Mechanism.** Find divisors of N in [x, x + X) with a Howgrave-Graham lattice (LLL).
- **Target.** The ladder's real need (divisor location), not floor sums.
- **Correctness.** KNOWN; the determinant inequality is checked exactly; E5 12/12.
- **Cost.** Polylog per block; blocks of length ≈ N^{α²} at x = N^α contain the j = 2 ladder blocks
  (PROVED, (α − ½)² ≥ 0). Total KNOWN N^{1/4+o(1)}.
- **New lemma for j ≥ 3.** Beat Coppersmith's interval N^{α²} at α ≈ ½, i.e. factor from fewer than
  (1/4)log₂N known high bits. OPEN since 1996.
- **Verdict.** Solves j = 2 (known exponent). The j ≥ 3 version is the classical open problem.

---
**F5. Voronoi / Li-type de-smoothing for restricted divisor sums.**
- **Mechanism.** Compute a smoothed Σ_{d∈I}⌊N/d⌋ fast; correct the boundary exactly.
- **Obstruction.** The correction set is {d : ‖N/d‖ < δ}: divisors of N − k for small k, including
  k = 0. Li corrects by factoring nearby integers; for a restricted sum that includes N. Circular.
- **Verdict.** Not applicable to factoring; works only for full D(N).

---
**F6. Stern–Brocot chain amortisation (universal-Euclid style).**
- **Mechanism.** Aggregate consecutive Farey windows that form an arithmetic chain of fractions.
- **Obstruction.** Generic chains have O(1) length (long only near major arcs); constant-factor gain.
- **Verdict.** HEURISTIC no-go for an exponent; not implemented.

---
**F7. Lattice window oracle at top blocks (bivariate Coppersmith; orchestrator request).**
- **Mechanism.** {P_j(h)} ∈ [0, ρ] ⇔ r − g(h) ≡ 0 (mod x^{j+1}) with 0 ≤ r ≤ ρD; basic-strategy
  lattice h^u f^k D^{m−k}; good rows vanish on all roots in the box.
- **Target.** The ladder's positive-width windows on top blocks x ∈ [√N/2, √N].
- **Correctness.** Howgrave-Graham (PROVED for the box); root extraction needs two independent rows
  (E6c: planted roots recovered 12/12).
- **Cost.** Measured X_max ≈ xN^{−1/3} (ratio 0.69–1.75, E6b), i.e. Hiary's scale: cN^{1/3}/B calls
  per block, N^{1/3} total. HEURISTIC explanation: the nonlinear part can be absorbed into a linear
  congruence only while X·N x^{j−2}X² ≲ x^{j+1}.
- **Verdict.** No gain over the trap, in this construction. Major-arc blocks (slope ≈ integer) are
  resolved whole.
