# ptecons round 2 (ADVERSARIAL): attack these construction-side claims

You are an adversarial referee. Below are claims from a research note about whether dense,
N-dependent families of split polynomials can give sub-L[1/3] integer factoring. Try to BREAK each
claim: find counterexamples, gaps in proofs, mislabelled statuses (PROVED vs HEURISTIC), overstated
scope, missed constructions that evade them, and missing or wrong citations. You may read
`notes/ptecons_results.md` (draft, may be in progress), `notes/systeq_results.md` and
`src/ptecons/*.py` in this repository, and you may run small Python checks in the current directory.
Prior blind-round output (yours) is `codex/ptecons_1.out`.

Setting. Relation = (u, v), u ≡ v (mod N), u ≠ v. u = P(m) = ∏_{l=1}^M (m − a_l) for a split
polynomial P with |a_l| ≤ A; v = (P mod F)(m) or the centered residue of u mod N. Target (systeq §6):
(a) ≳ log N/log B pieces of size ≤ B^{O(1)}; (b) log|v| = o(x^{2/3}) (x = log N), better a power
saving x^γ with γ < 2/3; (c) |𝓕| ≥ L[α], α < 1/3, computable without the factors.

## Claims to attack

1. **AVG [HEURISTIC + measured].** For any family U of integers fixed independently of N, the
   expected number of u ∈ U (|u| > N) with |u mod± N| ≤ V, over N drawn from a wide range, equals
   Σ_u (2V+1)/N up to (1+o(1)). Hence structure without N-dependence gives no excess yield.
   Measured: G = obs/exp = 1.05, 1.07, 1.04, 0.96, 0.99 for generic sextics, symmetric clusters,
   PTE-class members, Dickson integer-split members, composition clusters (200 semiprimes N of
   24 bits, V = 1000). Question: is there a rigorous version (divisor-sum / Barban–Davenport–
   Halberstam type average over N), and are there structured U that provably DEVIATE?

2. **GAP-R [PROVED].** If real y_1..y_M satisfy e_i(y) = e_{i+1}(y) = 0 for some 1 ≤ i ≤ M−1, then at
   most i − 1 of the y_l are nonzero. Proof: R(Z) = ∏(Z + y_l) is real-rooted; differentiate
   M − i − 1 times; the result is Z²·(...) since its Z¹ and Z⁰ coefficients are multiples of e_i, e_{i+1};
   a multiple root of a derivative of a real-rooted polynomial is a root of the polynomial of
   multiplicity ≥ (order of derivative) + (its multiplicity), so 0 has multiplicity ≥ M − i + 1 in R.
   **LIN-R [PROVED].** If P_λ(X) = ∏(X − α_l − L_l(λ)) with α_l ∈ C and L_l REAL linear forms on R^s,
   and the coefficients of X^{M−2} and X^{M−3} are affine in λ, then at most one L_l is nonzero.
   (Degree-i homogeneous part of e_i(α + L) is e_i(L); apply GAP-R pointwise with i = 2.)
   Interpretation: the BGJT mechanism (roots b + c·a, c ∈ F_q, all middle coefficients vanish because
   Σ_{c∈F_q} c^j = 0) has no real-direction analogue; complex directions need vanishing power sums
   p_2..p_k of the direction multiset, e.g. roots of unity (μ_M), which pushes the pieces into
   cyclotomic fields where the Splitting Lemma (norms of pieces ≥ N^{φ(n)/n}) applies.
   Attack: are there complex-direction families over a FIXED small field (e.g. Z[i], Z[ω]) with
   e_2..e_k ≡ 0 for k growing, small house, and useful norms? (We believe vanishing p_1..p_k in Q(i)
   is possible for all k by a convexity argument, but with heights growing — quantify.)

3. **CUN [PROVED].** For n ≥ 3 and integer roots a_l, with r = P mod (X^n − 1):
   ∏_{|a_l|≥2} |a_l| ≤ ‖r‖_1 · (1 − 2^{−n})^{−M/(n−2)}.
   Proof: |P(ζ)| = |r(ζ)| ≤ ‖r‖_1 at the nonreal n-th roots of unity; ∏_{ζ nonreal} |ζ − a| ≥
   |a|^{n−2}(1 − |a|^{−n}) for |a| ≥ 2 and ≥ 1 for a ∈ {−1, 0, 1}. So for Cunningham presentations
   (N | m^n − 1, X^n ≡ 1) every compressed split polynomial has at most log₂‖r‖_1 + o(1) roots outside
   {−1, 0, 1}: "the constants are not units" (over F_q all nonzero constants have absolute value 1).
   General form LEDGER: for monic F, ∏_l |F(a_l)| = |Res(F, P mod F)| ≤ ‖r‖_1^d M(F)^{deg r}.

4. **CW [PROVED, from your round 1].** Ideal classes (P − Q = C ≠ 0, both split over Z, degree M)
   have log|C| ≥ ½ M log M − O(M); with M ≥ x/(C₀ log B) a common-offset cluster has residual
   log-size ≥ (x/(2C₀ log B)) log(x/log B), which exceeds x^{2/3} exactly when log B = x^{α}, α < 1/3.
   Checked on the known ideal solutions of sizes 2–8, 10, 12.

5. **Alignment capacity / HALF [HEURISTIC + measured].** "Digit matching" — choosing m ≈ (kN)^{1/M}
   and forcing the top coefficients e_1, e_2, … of a split P to cancel the base-m digits of kN —
   needs |e_2| ~ m, hence roots A ≳ m^{1/2}; then the bottom coefficient e_M ~ A^M ≳ m^{M/2} = √N, so
   typical residuals are ≥ N^{1/2}. Generic PTE counting can only match O(√M) top coefficients, giving
   v ≥ N^{1−O(M^{−1/2})}. Measured (pre-registered, 64/80-bit semiprimes, pool 20000): median
   λ = log|v|/log N: COS (M = 2) 0.592/0.574; generic quartic 0.595–0.606/0.565–0.583;
   symmetric quartic 0.633/0.606; symmetric sextic two-digit matching 0–4 members (infeasible);
   one-digit sextic 0.675/0.648. No higher-degree aligned variant beat COS.

6. **Overall verdict [assessment].** No construction meets (a)–(c). Every candidate fails by one
   of: N-independence (AVG), orbit/Galois collapse and root-support rank, field-degree growth
   (Splitting Lemma), archimedean non-units (CUN/LEDGER), real rigidity (LIN-R), ideal cluster width
   (CW), or alignment capacity (HALF). The precise remaining primitive is an output-sensitive
   enumerator of split polynomials in the N-lattice {P : P(m) ≡ small (mod N)} — equivalently
   smoothenum's "doubly/algebraic listing", which is open and has no candidate.

## Tasks
- For each claim: verdict (holds / holds with fix / false), the fix, and the correct label.
- Find any construction family that escapes ALL seven failure modes simultaneously, even at the
  level of a heuristic, or argue why the failure modes are exhaustive for "split-by-identity" families.
- Check the cost threshold: is "γ < 2/3 (power saving)" the right requirement, and does CW's
  threshold land exactly at α = 1/3 or only up to logs?
- Literature: anything that already states LIN-R/GAP-R, CUN, CW (e.g. Kleiman, Rees–Smyth,
  Caley on PTE constants; Borwein; Costello–Meyer–Naehrig twin-smooth PTE sieving).
Label every statement PROVED / KNOWN (cite) / HEURISTIC / UNVERIFIED. Be concise and concrete.
