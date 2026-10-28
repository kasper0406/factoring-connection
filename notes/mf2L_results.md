# mf2L: transferring the U-MF trick to the two-level chain — results

Agent `mf2L`, phase 31, 2026-10-28. **v2.1**: v0 + proof-critic fixes (v1) + codex blind round's Drappeau CJM route
(§4.4–4.5) + the repairs from the literature critic and codex adversarial round mf2L_2 (MAJOR REVISION, no FATAL);
§8 records the reviews; codex round 3 (ACCEPT WITH MINOR FIXES) re-checked v2.

Files:
- this note;
- code: `src/mf2L/params.py` (exact-rational exponent margins of the input substitution, 16 checks + informational
  comparisons), `src/mf2L/checks.py` (exact checks of the finite statements: Props INV, AVG, NID, the d-th power
  kernel, Lemma GS in exact cyclotomic arithmetic);
- tests: `tests/test_mf2L.py` (7 tests, exact arithmetic);
- codex: `codex/mf2L_1.{md,out}` (blind round, launched before anything here was written), `codex/mf2L_2.*`
  (adversarial round on v1), `codex/mf2L_3.*` (re-check of v2).
- No experiments were run (every finite claim is an exact identity, checked by the tests), so there is no
  pre-registration.

Labels as in thmAB/c12rat: **proved** (internal, unrefereed, from the cited published results), **sketch**,
**analysis** (an argument about a proof architecture, not a theorem about factoring), **open**, **UNVERIFIED**.
Notation: λ := log log N, D := log N/λ, L[α, c] = e^{(c+o(1))D^αλ}. Design 2L_DW(α) of c12alg §1.1, α ∈ [3/8, 2/5).
Ψ(z, y; χ) := Σ_{n≤z, P⁺(n)≤y}χ(n); u = log z/log y.

---

## 0. Verdict

**Plain language.** The U-MF trick itself does not carry over, but the goal is reached by another route. Theorem 2L
(rigorous L[α] for every fixed α ∈ (3/8, 2/5), and L[3/8] under the old constraint a > 2c₂) can be made to rest on
**published theorems with full proofs only**. As before, it also rests on the project's internal chain, which no human
has refereed. Drappeau's Lemme 5 and Harper's unpublished preprint are no longer used.
- **What the integer multiplier does.** It is already part of the two-level design (design DW). It turns every
  *one-point* estimate into an interval density: level-1 first moments, child means, parent counts, and U′. So Lemma
  8.5H is not needed. The remaining use of Theorem U, factoring values recursively, can call the BV-free FACTOR_MF.
- **What it cannot do, in this proof architecture.** The two-level proof bounds the dead parent–child incidences by a
  union bound. Each child must then fail with probability ≪ η ≍ D^{−(3−7α)}, so pair covariances are needed to
  *second order*. Two facts block an interval-density proof:
  - for every integral form F and points P₁, P₂, F(P₁) ≡ λ^dF(P₂) (mod det(P₁, P₂)); no multiplier, shift, extra
    rational side or MF family of fields changes this (Prop INV);
  - interval densities do not determine pair counts on such lattices (Prop NID).

  So characters modulo pair determinants enter the covariance bound. MF's independent algebraic sides do not remove
  this: per field, the event "parent accepted and all children good" is a joint event that nobody can bound (§3).
  These are obstructions to this architecture and to the MF2 design. They are not a proof that no other method could
  avoid characters.
- **The BV-type input needed is weak, and published.** The four character-error families need
  (WA-hyp): Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀}|û(χ)|² ≤ εU², with ε ≈ e^{−CD^θλ²} (the old proof had e^{−cD^α/λ²}). We prove it
  from two pieces.
  - (a) *Conductors above R₁ = e^{C_RD^θλ²}.* Harper's large-values bound (4.1) (Compositio 152 (2016), Prop. 2). It is
    published, and its proof does not use the preprint, unlike his Theorem 2. See Lemmas GS, LV, IND and Prop SUB.
    Alternatively (codex), the classical large sieve plus a log-free zero-density estimate, if the constant a in B is
    large (§4.5, sketch).
  - (b) *Small conductors.* A pointwise bound, Hyp SC: |Ψ(z, y; χ)| ≤ e^{−u^κ}Ψ for cond χ ≤ e^{u^κ}. This follows
    from **Drappeau, Canad. J. Math. 67 (2015), Props 2.9 and 2.13** (Lemma SC-D), which have full proofs, include
    the exceptional (Siegel) character, and do not use Siegel's theorem. Codex's blind round found these propositions.
    The literature critic and codex then checked them against the published PDF.
- **Theorem 2L-MF(α) [proved relative to the internal chain; §5].** For every fixed α ∈ (3/8, 2/5), and at α = 3/8 with
  a > 2c₂, FACTOR_2L (design DW, value factorisations by FACTOR_MF) is Las Vegas and factors every N in expected time
  L_N[α, c(α)]. The exponent and, via route (a), the constant c(α) are those of Theorem 2L-final. The bounds are not
  uniform as α → 2/5.
- **Minimal BV-type statement (Q5).** In isolated form: Hyp SC(κ) for some κ ∈ (κ(α), 1), plus any mean-square
  large-conductor bound with polylogarithmic loss. Here κ(α) = max(θ/α, θ_p/(3 − 7α)), so κ(3/8) = 2/3 and κ → 1 as
  α → 2/5. Lemme 5 implies both. So do Drappeau CJM and Harper (4.1), which is what we use.

**Caveats.**
- Everything internal (thmAB, c12rat, c12alg, uprime, pc3, multidescent, U-MF) remains unrefereed. 2L-final's own
  reviewers put it at ≈ 0.65 (window) and ≈ 0.52 (3/8). This note changes only the external inputs; `audit2L` audits
  the rest.
- Harper's (4.1) was read in the accepted manuscript, not the typeset journal version. Iwaniec–Kowalski theorem numbers
  are taken from secondary citations.
- The v2 repairs were re-checked by codex round 3 (ACCEPT WITH MINOR FIXES, applied in v2.1). No Opus critic re-read v2.

**Probabilities (ours, after review; reviewers' in §8).**
- Dependency map complete and correct: 88%. [proof critic 90%]
- Props INV, AVG (corrected phase) and NID: 97%.
- The substitution (GS, LV, IND, SUB, SC-D) after the v2 repairs: 88%. [codex mf2L_3: 92% for 2L-MF relative to the
  chain; proof critic: 85–90% on v0/v1]
- The literature facts (Drappeau CJM full proofs; (4.1) independent of the preprint): 90%. [literature critic: FULL
  PROOF for both]
- Theorem 2L-MF(α) correct relative to 2L-final's internal chain: ≈ 0.8. Unconditionally, including the chain:
  ≈ 0.5 (window) and ≈ 0.42 (3/8), i.e. 2L-final's probabilities times 0.8.

## 1. Dependency map (Q1)

Where 2L-final (thmAB §0, §12; c12alg §5; c12rat §9; pc3; uprime; multidescent) uses a BV / BDH / large-sieve input,
and whether an integer-multiplier or interval-density argument replaces it.

| # | Place | Input used in 2L-final | Kind of statement | Replaceable by OP-type / elementary? |
|---|---|---|---|---|
| 1 | level-1 first moment FM-U′ (algebraic factor) | under DW: Lemma DW(a)⁺ (pointwise); Theorem U's version used Lemma 8.5′ | one-point lower bound | **already replaced** (c12alg §2.7): uniform residue mod G(P) ⇒ interval density, exactly as Lemma OP |
| 2 | Lemma 8.5H (uprime, hybrid BV) | not used under DW (c12alg §2.7) | — | **already gone** |
| 3 | level-1 second moment LD*-U′ | Lemma W-DW (all residues), Kern-U, DS-U′, L2-U′ (classical large sieve, Montgomery–Vaughan) | upper bound, polylog loss | **BV-free**. The classical large sieve is a published theorem (not BV-type). MF would remove LD* altogether, but it is not needed for that |
| 4 | child one-point mean (Prop CLB) | DW(a)⁺ with root coset, Lemma FIB, Lemma U1 | one-point lower bound | **BV-free** (pointwise) |
| 5 | parent count first moment, Lemma η, Lemma W-P(P3) | Hildebrand at B₂, DW(a)⁺, FIB | one-point | **BV-free** |
| 6 | pair upper bounds (Lemma PUB, PU-A, BTS), thin pairs, skew (DOM, SPR, LS, TF, MULT) | elementary sieve-free bounds | upper bounds | **BV-free** |
| 7 | one-point cell moments (Lemma FO, OPV, OPV-A, LF-A) | Hildebrand short intervals, HT, dBT | first-order one-point | **BV-free** |
| 8 | **C₁ rational** (c12rat RAT-1, Lemma EN + Theorem WA + FAM) | **Drappeau Lemme 5** (conductors ≤ 4Q₂S₂², contains q) | second-order pair covariance | **No.** Characters mod qΔ′/δ are forced (Prop INV); replaced in §4 by Harper (4.1) + Hyp SC |
| 9 | **C₂ rational** (RAT-2, library 𝒜^J) | **Drappeau Lemme 5** (conductors ≤ 2V₁²) | same | **No**; §4 |
| 10 | **C₁ algebraic** (thmAB AB₂(i), Lemma EN-A + FAM-A + WA) | **Drappeau Lemme 5** (conductors ≤ 2S₂²) | same | **No**; §4 |
| 11 | **C₂ algebraic** (AB₂(ii), parents at B₂) | **Drappeau Lemme 5** (conductors ≤ 2V₁²/q) | same | **No**; §4 |
| 12 | recursive value factorisations (uprime §5: cofactors ≤ X_F by FACTOR_U) | Theorem U ⇒ Lemma 8.5′ ⇒ Lemme 5 | oracle | **Replaced**: call FACTOR_MF (Theorem U-MF), same exponent L_N[7/20]-scale = L[α]^{o(1)}; legitimate by Lemma OS (transcript law independent of the exact oracle) and Lemma W |
| 13 | post-processing (LV6) | Lagarias–Odlyzko, Stark, Landau | — | not BV |
| 14 | Prop UB / DC / DC_e, Lemma VD, Prop PB | elementary | reduction | — |

So in 2L-final, **Drappeau's Lemme 5 enters exactly through the four character families 8–11, plus the
FACTOR_U oracle (12)**. Item 12 is removable at no cost. Items 8–11 are the subject of §§2–5.

## 2. Why second-order pair estimates force characters (Q2, Q3)

**Prop INV (homogeneity is design-invariant) [proved; pc3 Lemma LAT(c)].** Let P₁, P₂ ∈ Z² be primitive and
non-proportional, Δ := |det(P₁, P₂)|, and λ the unit mod Δ with P₁ ≡ λP₂ (mod Δ). Then for **every** integral binary
form H of degree e, H(P₁) ≡ λ^eH(P₂) (mod Δ).
*Proof.* H(P₁) ≡ H(λP₂) = λ^eH(P₂), coefficientwise mod Δ. ∎ (`checks.homogeneity`: 400 random forms of degree
1–7 with coefficients up to 10⁶.)

*Consequences.*
- *Rational side (e = 1).* Whatever the law of g, and whatever further randomisation keeps the rational side a linear
  form, (G(P₁), G(P₂)) lies on {n₁ ≡ λn₂ (mod Δ)}. A multiple of G only rescales the condition. A shift G + N·(·) makes
  the values ≈ N, which destroys smoothness. So neither a random shift, nor a second multiplier, nor randomising g makes
  the rational residues mod Δ uniform: the obstruction is the residue *relation* between the two coordinates, not the
  residue of either one. Prop INV covers one form at two points. Using *different* rational forms G^{(i)} at
  different child points is not excluded by INV (proof critic #7), but each extra rational root shared with F costs
  about N^{1/d} in F's height (L[α]^{−Θ(1)} at child scale), and with r rational sides the same-side pairs (a fraction
  ≈ 1/r of all pairs) still need precision η. So the conclusion stands. [analysis]
- *Algebraic side (e = d).* The same for F = wF̂_N + G·R_c, for any further multiplier, and for any of the J_f
  independent MF sides: their pair laws all live on L′ = {x₁ ≡ λ^dx₂ (mod Δ/q)} (c12alg Cor TP-q).
- The integer multiplier w randomises the residue of F(P) **modulo G(P)**, which is a one-point statement. For a pair it
  makes (F(P₁), F(P₂)) uniform on L′ (DW-L), which is the best possible. It cannot touch the congruence between the
  coordinates.

**Prop AVG (where characters come from) [proved; exact].** For finitely supported a, b ≥ 0 and k ≥ 1, write A*, B* for
the sums over arguments coprime to k. Then Σ_{ι ∈ (Z/k)^×}#_{a,b}{x₁ ≡ ιx₂ (mod k)} = A*B*, where #_{a,b} weights pairs
by a(x₁)b(x₂) and x₂ is restricted to units. For a fixed ι the count is A*B*/φ(k) + φ(k)^{−1}Σ_{χ≠χ₀}χ̄(ι)â(χ)\overline{b̂(χ)} (orthogonality: 1[x₁ ≡ ιx₂] = φ(k)^{−1}Σ_χχ(x₁)χ̄(ιx₂) on units; v1 wrote χ(ι), codex mf2L_2).
*Proof.* Each pair of units (x₁, x₂) satisfies x₁ ≡ ιx₂ for exactly one ι. The second statement is orthogonality.
∎ (`checks.iota_average`, six moduli.)
So the character term vanishes if (sufficient, not necessary: codex mf2L_2) the pair's unit ι is equidistributed over
the pair family. Otherwise it must be bounded through the twisted sums â(χ), i.e. through the distribution of smooth numbers in residue classes.

**Remark (d-th powers) [proved].** On the algebraic side ι = λ^d. A character χ with χ^d = χ₀ has χ(ι) = 1 for every
pair, so averaging over pairs never cancels it. Nontrivial such characters exist mod k iff gcd(d, φ(k)) > 1
(`checks.dth_power_kernel`). For quadratic χ and odd d, χ(λ^d) = χ(λ) (`checks.quadratic_odd_power`). So characters
of order not dividing d, in particular the real ones, could in principle benefit from pair-averaging on both sides
(§7, open route); those of order dividing d cannot.

**Prop NID (interval densities do not determine pair counts) [proved; model].** Let k ∈ {3, 4, 6}, χ the real
non-principal character mod k, s(n) := 1[(n, k) = 1] and t(n) := s(n)(1 + εχ(n)). On every interval I of length ≥ k,
|Σ_It − Σ_Is| = ε|Σ_Iχ| ≤ ε, so the two weightings have the same interval densities up to ε/|I|. But with
c_r := #{n ≡ r (mod k)} in the support, Σ_{x₁≡x₂ (k)}t(x₁)t(x₂) − Σ_{x₁≡x₂ (k)}s(x₁)s(x₂) = Σ_{r unit}c_r²(2εχ(r) + ε²),
which is ε²Σ_rc_r² > 0 when the c_r are equal on units. (`checks.nid_example`: k = 3, 4, 6, ε = 1/2.)
*Consequence.* A covariance bound on L′ cannot follow from interval densities (Hildebrand–Tenenbaum short-interval
asymptotics) alone. It needs information on smooth numbers in residue classes modulo the divisors of Δ′. Q3's
answer is therefore **no**. What *can* be proved from interval densities is the principal-character part, which is
what Lemmas CF/CP and OPV-A already do.

**Precision required by this architecture (pc3 Prop DC, c12alg Prop PB, multidescent Lemma MP).** The dead-incidence
bound is W(A + √(AB)) with η := μ₁/W ≥ (3u₂)^{−1} ≍ D^{−(3−7α)} (c12alg Lemma η: the *averaged* number of children per
accepted parent is ≤ 3u₂, while the deterministic cap is D^{2−4α}; that η is not much larger is heuristic). Chebyshev-type union bounds therefore need A ≪ η and
√(AB) ≪ η, i.e. relative variances to precision D^{−(3−7α)}·λ^{−C} or better (at 3/8: D^{−3/8}). Upper bounds with polylog
loss (Lemma DS-type, the U-MF rational device) give only Paley–Zygmund success (log N)^{−C} per child and cannot make
the dead mass small. This is a statement about the union-bound/Chebyshev architecture, not about every conceivable
proof (codex mf2L_2).

## 3. Design changes (Q2, Q4) [analysis]

*MF at level 2.* Design MF2: shared g, J_f independent (w_j, c_j), first rationally good parent P₀(g, ℓ). Fields are
independent given g, so a per-field success probability p(g, ℓ) ≥ ρ₂·L^{−O(1)} would suffice, and losses
(log N)^{O(m)} from m children would be search work, not root density (m ≤ the cap D^{2−4α}; on incidence average
m ≤ 3u₂ ≍ D^{3−7α}, proved, and ≍ u₂ heuristically). But:
- algebraic children of field j must be descended in field j (Cor EX), so per field the event is the **joint** event
  "F_j(P₀) accepted and all m children of P₀ in field j good". Heuristically m ≍ u₂ ≍ D^{3−7α} > d ≍ D^{1−2α} (as
  3 − 7α > 1 − 2α for α < 2/5; c12alg Lemma η proves only the upper bound W/μ₁ ≤ 3u₂, so this typical count is
  [heuristic], codex mf2L_3), so it concerns more than d + 1 points of one form, and no rank argument makes the
  children independent (Lemma OD / Prop R3). Codex's countermodel: making exactly one uniformly chosen child dead per
  accepted parent gives per-child failure 1/m and joint success 0, so per-child information cannot bound the joint
  event.
- a lower bound for that joint event goes either through a union bound (per-child failure ≪ 1/m again), or
  through a joint estimate (LS_C). The union bound must decouple the parent from its children. pc3's Prop DC does that
  with the parent **count** Z_𝔮, which concentrates because many parents share a label. In MF2 the parent of (g, ℓ) is
  one indicator 1[acc(P₀) ∧ q | F(P₀)]. Its relative standard deviation is ≈ ρ₂^{−1/2}, which is L^{Θ(1)} at 3/8 and
  e^{Θ(D^{3−7α}λ)} ≫ D^C for α ∈ (3/8, 2/5), so Cauchy–Schwarz decoupling fails on the whole range. Summing over ℓ restores a parent count and with it C₁ and C₂, with the same
  character families. The other option is conditioning on the exact value F(P₀), which fixes the residues of F(P′)
  modulo det(P₀, P′)/q (Prop INV, three points). That needs smooth numbers in one *prescribed* class modulo moduli of
  size L[1 − α]^{O(1)}, which is harder.
- conclusion: MF removes the level-1 algebraic second moment, which is not a BV input anyway (§1 item 3). It does not
  remove C₁ or C₂ *in this architecture* (an MF2 proof would need a new joint lower bound p_joint ≥ Pr(A_j)e^{−Cmλ},
  which nobody has; codex mf2L_1 (9)). This agrees with multidescent §4's table ("Prop MF fields sharing g … does not export a child").

*Extra rational randomness* (§2, Prop INV): none available without breaking the root or the height.

## 4. The input substitution: Drappeau's Lemme 5 → Harper (4.1) + Hyp SC

All four families use Lemme 5 only through the library-energy hypothesis of c12rat's Theorem WA:

  (WA-hyp) Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀ (mod k)}|û_j(χ)|² ≤ εU_j² for every library sequence u_j,

followed by ε·𝔇 ≪ (budget), with 𝔇 ≤ e^{C_𝔇D^θλ²} (children, rational parents) or e^{C_𝔇D^{θ_p}λ²} (algebraic
parents) (thmAB Lemma FAM-A, c12rat Lemma FAM). The library sequences are:
- thmAB children/parents: u(m) = 1[m ≤ z, P⁺(m) ≤ y], y = B (children) or B₂ (parents);
- c12rat children: u(m) = 1[m ≤ z, P⁺(m) ≤ B, m P₀-free];
- c12rat parents: u = Σ_{ℓ∈J}1[ℓm, m ≤ z/ℓ, P⁺(m) ≤ B, m P₀-free];
with z ≥ x^{1/2} at the family's value scale x and Q ≤ x^{o(1)} (`params`: "Farey points").

### 4.1 Two lemmas on one sequence

Let a = (a_n)_{n≤z} with |a_n| ≤ 1, supported on y-smooth n, S(θ) := Σa_ne(nθ), Ψ := Ψ(z, y), â(χ) := Σa_nχ(n).

**Lemma GS (Gauss transfer) [proved; textbook].** For every k ≥ 2, Σ*_{χ mod k}|â(χ)|² ≤ Σ_{(b,k)=1}|S(b/k)|² (Σ* over
primitive χ).
*Proof.* For primitive χ, Σ_{b mod k}χ̄(b)e(bn/k) = χ(n)τ(χ̄) for all n, and |τ(χ̄)|² = k. So
|â(χ)|² = k^{−1}|Σ_{(b,k)=1}χ̄(b)S(b/k)|². Summing over all χ mod k and using orthogonality,
Σ_χ|Σ_bχ̄(b)S(b/k)|² = φ(k)Σ_{(b,k)=1}|S(b/k)|². Drop the imprimitive χ and use φ(k)/k ≤ 1. ∎
(`checks.gauss_transfer`: both identities exactly in Q(ζ_{kφ(k)}) for k = 5, 7, 9.)

**Fact H4.1 (Harper, Compositio Math. 152 (2016) 1121–1158, Proposition 2, bound (4.1); statement as in arXiv v1).** There is an absolute C such that: let
log^{21}z ≤ y ≤ z^{1/100} be large, (a_n)_{n≤z} complex with |a_n| ≤ 1 supported on y-smooth n, z^{−1/1000} < δ ≤ 1, and
θ₁, …, θ_R with ‖θ_r − θ_s‖ ≥ 1/z (r ≠ s) and |S(θ_r)| ≥ δΨ(z, y) for every r. Then
R ≪_ε δ^{−2−20(1−α(z,y))−ε}(log z)^{3+32(1−α(z,y))+ε}.
*Status (verified by the literature critic and by codex mf2L_2 in the accepted journal manuscript, WRAP eprint 107457;
same numbering as arXiv v1; the typeset version was not accessible).* Harper deduces his Theorem 2 (restriction) from
(4.1) **and** (4.2). Only (4.2) and the §5 applications use the major-arc estimates, which rest on §3.3 of
arXiv:1208.5992 (Appendix A.1). (4.1) is proved via his Props 3–4 from Theorem 1, Erdős–Turán, Bourgain's lemma,
Smooth Numbers Results 2–3 (de la Bretèche–Tenenbaum; proved in the paper) and Fouvry–Tenenbaum's minor-arc bound.
So **(4.1) has a full published proof**; Theorem 2 as stated does not.

**Lemma LV (large conductors) [proved, given Fact H4.1].** Let 1 − α(z, y) ≤ 1/40, R ≥ 2 and 8R² ≤ z^{1/1000}. Then
  Σ_{R<k≤2R}φ(k)^{−1}Σ*_{χ mod k}|â(χ)|² ≤ C·(log log 3R)·R^{−1/3}(log z)^{3}·Ψ(z, y)²
(in Fact H4.1's domain: log^{21}z ≤ y ≤ z^{1/100}, z large).
*Proof.* Let Θ be the Farey fractions b/k ∈ (0, 1), k ≤ 2R, (b, k) = 1. They are pairwise ≥ 1/(4R²) ≥ 1/z apart, and
#Θ ≤ 4R². |S| ≤ Ψ. Apply Fact H4.1 with ε = 1/4 to Θ_j := {θ ∈ Θ : |S(θ)| ≥ 2^{−j}Ψ} for 2^{−j} > z^{−1/1000}:
#Θ_j ≤ C2^{j(2.75)}(log z)^{4.05}, since 20γ ≤ 1/2 and 32γ ≤ 0.8. Hence
Σ_{θ∈Θ}|S(θ)|³ ≤ Σ_j(2^{1−j}Ψ)³#Θ_j + #Θ·z^{−3/1000}Ψ³ ≤ C(log z)^{4.05}Ψ³ + Ψ³.
By Hölder, Σ_{θ∈Θ}|S(θ)|² ≤ (#Θ)^{1/3}(Σ|S|³)^{2/3} ≤ CR^{2/3}(log z)^{2.7}Ψ². Lemma GS, summed over R < k ≤ 2R (each
b/k occurs once), bounds Σ_kΣ*_χ|â(χ)|² by the same quantity. Finally φ(k) ≥ k/(C log log k). ∎

The exponent 1/3 is not optimal (any p ∈ (2 + 20γ, 4) works). Only a power saving in R with a polylogarithmic loss
is needed. The essential point is that the loss is **polylogarithmic**, not 1/ρ as in the classical large sieve: that
closes the λ-gap (§0).

### 4.2 From primitive characters to the WA hypothesis

**Lemma IND (imprimitive characters) [proved].** Let u be one of the library sequences of §4 with y-smooth support
and endpoint z, Q ≤ z^{1/2002} with z large (so that 8Q² ≤ (z/Q)^{1/1000} and Lemma LV applies at every endpoint
z/e ≥ z/Q), normaliser U(z′) := Ψ(z′, y) (the *unfiltered* count, as in c12rat Lemma EN; for P₀-free u this is
larger than #supp u by ≤ (log P₀)^C), and put T(z′) := Σ_{1<k₀≤Q}φ(k₀)^{−1}Σ*_{χ mod k₀}|Ψ_u(z′; χ)|²/U(z′)², where
Ψ_u(z′; χ) := Σ_{m≤z′}u(m)χ(m) is the same kind of sequence with endpoint z′. Then
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀ (k)}|û(χ)|² ≤ Cλ^{C}(log Q)·U(z)²·sup_{z/Q≤z′≤z}T(z′).
*Proof.* A non-principal χ mod k is induced by a unique primitive χ* mod k₀, k₀ | k, k₀ > 1. Let k₁ be the product of
the primes dividing k but not k₀. Then û(χ) = Σ_{e|k₁}μ(e)χ*(e)V_e with V_e := Σ_mu(em)χ*(m). If e is not y-smooth,
V_e = 0. Otherwise V_e is the same kind of sequence with endpoint z/e (P₀-freeness: V_e = 0 if e has a prime ≤ P₀).
Cauchy–Schwarz with σ := 1 − 1/λ gives |û(χ)|² ≤ ∏_{p|k}(1 + p^{−σ})·Σ_{e|k₁}e^{σ}|V_e|². Here
∏_{p|k}(1 + p^{−σ}) ≤ λ^C: the product is largest when the primes are the first ω(k) ≤ log Q primes, and
p^{1/λ} ≤ e for p ≤ D. Next, Σ_{k≤Q: k₀e|k}φ(k)^{−1} ≤ C log Q/(φ(k₀)φ(e)), because (k₀, e) = 1. So the left side is
≤ λ^C log Q·Σ_e(e^{σ}/φ(e))Σ_{k₀}φ(k₀)^{−1}Σ*_{χ*}|V_e(χ*)|². This is ≤ λ^C log Q·Σ_e(e^{σ}/φ(e))U(z/e)²T(z/e). Fact dBT1
gives U(z/e) ≤ Ce^{−β}U(z) up to (log P₀)^C. Then Σ_e e^{σ−2β}/φ(e) ≤ ∏_p(1 + Cp^{−2+2γ−1/λ}) = O(1). ∎
**The J-library (codex mf2L_2, MAJOR #3).** 𝒜^J = {ℓm} is *not* supported on B-smooth integers, so Lemma IND does not
apply to it directly. Apply IND (and Prop SUB below) to each smooth cofactor sequence A_{z/ℓ} first, then sum the fibres
with Minkowski's inequality for the energy norm 𝓔_Q(·)^{1/2} := (Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀}|·̂(χ)|²)^{1/2}: since
û^J(χ) = Σ_ℓχ(ℓ)Â_{z/ℓ}(χ) and |χ(ℓ)| ≤ 1, 𝓔_Q(A^J)^{1/2} ≤ Σ_{ℓ∈J}𝓔_Q(A_{z/ℓ})^{1/2} ≤ √ε Σ_{ℓ∈J}Ψ(z/ℓ, B), which is
(WA-hyp) with c12rat's normaliser Σ_ℓΨ(z/ℓ, B). No factor #J and no assumption on (ℓ, k) is needed. The endpoints z/ℓ
are ≥ x^{1/3} because log ℓ = o(log x).

**Hyp SC(κ) (small conductors) [hypothesis; proved below as Lemma SC-D].** There is x₀ such that for all x ≥ x₀ and
y with (log x)^{1/2} ≤ log y, y ≤ x^{1/100} and u := log x/log y ≥ (log log x)⁵: for every primitive χ with
1 < cond χ ≤ exp(u^κ) and every x^{1/3} ≤ z′ ≤ x, |Ψ(z′, y; χ)| ≤ exp(−u^κ)Ψ(z′, y). (v1 allowed (log x)^{1/3} ≤ log y;
codex mf2L_2 MAJOR #1–2: on that larger domain log y ≥ u fails and Lemma SC-D does not apply. All four families satisfy
the restricted domain at every endpoint, `params` "SC-D regime": log y/√(log x) ≍ √λ for children and rational parents,
≍ D^{(11α−4)/2}√λ for algebraic parents.)
*Effectivity.* Theorem 2L-SC's conclusion (an expected-time bound) does not need x₀ to be effective: the algorithm uses
no analytic constant (alternation with trial division, as in Theorem U-MF), so an ineffective x₀ only makes the o(1)
in L_N[α, c(α) + o(1)] ineffective. If Hyp SC is proved via Siegel's theorem for a possible exceptional character, this
is what happens. An effective version is needed only for an effective N₀.

**Remark (Lemme 5 ⇒ Hyp SC(κ) for every κ < 1) [proved; proof critic #3].** Keep only the term k = k₀ = cond χ in
Drappeau's sum Σ_{k≤Q}φ(k)^{−1}Σ_{χ mod k, 1<cond χ≤x^{η_D}}|Ψ(x, y; χ)| ≪ Ψ(x, y)(H(u)^{−δ_D} + y^{−δ_D}). This gives
|Ψ(z′, y; χ)| ≤ Cφ(k₀)Ψ(z′, y)(H(u′)^{−δ_D} + y^{−δ_D}) at z′ (u′ ≥ u/3). For k₀ ≤ e^{u^κ} and κ < 1 this is
≤ e^{−u^κ}Ψ(z′, y) for large x, since H(u′)^{δ_D} = e^{δ_Du′/log²u′} ≫ e^{2u^κ}, and on Hyp SC's domain
log y ≥ (log x)^{1/2} ≥ u, so y^{−δ_D} ≤ e^{−δ_Du} ≤ e^{−2u^κ}. So Theorem 2L-SC contains
Theorem 2L-final.

**Prop SUB (substitution) [proved, given Fact H4.1 and Hyp SC(κ) with κ > κ(α)].** In each of the four families,
(WA-hyp) holds with ε ≤ Cλ^C(log x)^{4}(R₁^{−1/3} + R₁e^{−2u^κ}) for every 2 ≤ R₁ ≤ e^{u^κ}. Let 𝔇̃ := 𝔇/Σŵ be
the *normalised* family constant (c12rat's 𝔇 contains Σŵ ≍ μ₁·polylog; thmAB's is already normalised), with
log 𝔇̃ ≤ C_𝔇D^{θ}λ² (children, rational parents) or C_𝔇D^{θ_p}λ² (algebraic parents). Put
R₁ := e^{C_R D^θλ²} with C_R > 3(C_𝔇 + C_T + 2) chosen after C_T and C_𝔇 (children and rational parents; T₀ = e^{C_TD^θλ²}),
and R₁^p := e^{C_R D^{θ_p}λ²} with the **parent** threshold ϑ₀^p = e^{C_TD^{θ_p}λ²} in place of T₀ (algebraic parents;
proof critic #1: the child T₀ would make R₁ ≤ e^{u₂^κ} impossible for α ≥ 7/18). Then the normalised character totals are
O(e^{−(C_T+1)D^θλ²}) (resp. O(e^{−(C_T+1)D^{θ_p}λ²})).
*Downstream (codex mf2L_2 MAJOR #5: this is not a literal replacement of the old numerical bound e^{−cD^α/λ²}).* The new
totals are absorbed at each later use: thmAB §12(i): the child total is o(T₀^{−1}); thmAB §12(ii): the parent total is
smaller than every power of D and is absorbed into the polynomial parent term; c12rat §8.3: summing the dominating
weights costs a polynomial factor times μ₁ (Lemma DOM(ii)), and W ≥ (1 − o(1))μ₁ (DOM(iii)), so both rational totals
are absorbed into the polynomial RAT-1/RAT-2 terms after increasing C. One retains A ≪ D^{−2α}λ^C and
B ≪ (D^{−2α} + D^{−2(3−7α)})λ^C, hence D^{3−7α}(A + √(AB)) ≪ λ^C(D^{3−9α} + D^{−α}) = o(1) on [3/8, 2/5): the incidence
budget of c12alg Prop PB / pc3 Cor DC holds, including at 3/8.
*Proof.* Split the primitive conductors in T(z′) at R₁. For k₀ > R₁, sum Lemma LV over dyadic R ≥ R₁ (using
P₀-freeness and J-fibres as above; Fact H4.1 allows any weights |a_n| ≤ 1 on smooth n). This gives
Cλ^C(log x)^{3}R₁^{−1/3}. For k₀ ≤ R₁, Hyp SC (at κ with R₁ ≤ e^{u^κ}) gives Σ_{k₀≤R₁}φ(k₀)^{−1}·φ(k₀)e^{−2u^κ} ≤ R₁e^{−2u^κ}. For P₀-free
sequences, Möbius over d | P₀# turns Hyp SC into the same bound with a factor (log P₀)^C (as in c12rat Lemma EN; the
endpoints z′/d ≥ x^{1/3}). Lemma IND assembles the families. The choice of R₁ is admissible: log R₁ = O(D^θλ²) and
u^κ ≍ D^{ακ} with ακ > θ (resp. log R₁^p = O(D^{θ_p}λ²) and u₂^κ ≍ D^{(3−7α)κ} with (3 − 7α)κ > θ_p), so
R₁e^{−2u^κ}𝔇̃ = o(T₀^{−1}e^{−D^θλ²}). Hypotheses of Fact H4.1: y ≤ x^{1/100}, 8Q² ≤ x^{1/1000},
1 − α ≤ 1/40, log^{21}x ≤ y (`params`, all families). ∎
For the algebraic parents, θ is replaced by θ_p, u by u₂ ≍ D^{3−7α} and T₀ by ϑ₀^p throughout.

### 4.3 The constant κ(α) (`params.kappa`)


κ(α) := max over the families of (loss exponent)/(u-exponent) = max(θ/α, θ_p/(3 − 7α)). On [3/8, 2/5) this gives
θ = 2α − 1/2 (children, both rational and algebraic, and c12rat's parent family, whose h_* contains the child
threshold) and θ_p = max(1 − 2α, 2 − 5α) (algebraic parents). Values: κ(3/8) = 2/3; κ(19/50) = 12/17; κ(7/18) = 4/5;
κ(39/100) = 22/27; κ → 1 as α → 2/5. Hyp SC(κ) is needed for some fixed κ ∈ (κ(α), 1). The families' (x, y):
children x = e^{Θ(D^{2α}λ)}, y = B; algebraic parents x = e^{Θ(D^{2−3α}λ)}, y = B₂; rational parents x = e^{Θ(D^{2α}λ)},
y = B. In every case log y ≥ (log x)^{1/2} (strictly, at leading order), y ≤ x^{1/100} and u ≥ (log log x)⁵.

**What Hyp SC is.** The conductors exp(u^κ) are y^{o(1)} but far beyond the Siegel–Walfisz range exp(c√log y)
(`params.informational`; exception: the algebraic parents inside the open window). The saving e^{−u^κ} is far weaker
than the truth (≈ e^{−cu} for non-exceptional χ). The standard proof is a zero-free region for L(s, χ) at conductor
y^{o(1)}, with at most one exceptional real χ₁ (Landau–Page), plus a Rankin / saddle-point bound for the twisted Euler
product ∏_{p≤y}(1 − χ(p)p^{−s})^{−1}. Heuristically an exceptional zero *helps* (χ₁(p) = −1 for most small p makes the
twisted sum smaller), but this must be proved. This is exactly what §3.3 of Harper's arXiv:1208.5992 does (and why
Lemme 5 and Harper 2016's Appendix A.1 cite it).

### 4.4 Hyp SC from published theorems (Q5): Lemma SC-D

Source, found by codex's blind round (`codex/mf2L_1.out` §5.1) and checked against the published CJM PDF by the Opus
literature critic and by codex mf2L_2: S. Drappeau, *Sommes friables d'exponentielles et applications*, Canad. J.
Math. 67 (2015) 597–638 ("[D-CJM]"). CJM numbering differs from arXiv 1302.4318 and the author copy (CJM 2.9 = arXiv
2.1; 2.10 = 2.2, the one Harper 2016 cites; 2.13 = 2.5; (2.16) = (2.20)). Prop 2.13 is already in arXiv v1.
- **Exceptional characters, [D-CJM] (1.2).** There is b > 0 such that for all Q, T ≥ 2, L(s, χ) has no zero with
  σ ≥ 1 − b/log(QT), |τ| ≤ T for every χ mod q ≤ Q, except possibly for the characters induced by one primitive
  quadratic χ₁ (mod q₁), which then has a unique simple real zero β₁; χ₁ is called *Q-exceptional* (for the height T).
  This is the Landau–Page theorem (Iwaniec–Kowalski Thms 5.26, 5.28; Davenport Ch. 14; IK numbering from secondary
  citations, UNVERIFIED against the book).
- **[D-CJM] Prop 2.9** (near-verbatim). There are c₁, c₂ > 0 such that, for 4 ≤ T ≤ y^{c₂}, ε > 0, 2 ≤ q ≤ y^{c₂} and
  χ mod q non-principal with L(s, χ) ≠ 0 on {σ > 1 − ε, |τ| ≤ T}: if (log x)^{c₁} ≤ y ≤ x and z ∈ [x^{2/3}, x], then
  Ψ₀(z, y; χ, γ) ≪ z^αζ(α, y)(1 + |γ|z)((log T)x^{−c₂ε} + T^{−c₂}), with Ψ₀(z, y; χ, γ) := Σ_{n∈S(z,y)}e(nγ)χ(n) and
  α = α(x, y). *Hidden hypotheses* (its proof uses Lemma 2.7): ε log y ≥ 1/c₂, T ≥ y^{c₁ε}(log x)², and the zero-free
  region up to |τ| ≤ 2T.
- **[D-CJM] Prop 2.13** (near-verbatim). For (log x)^{c₁} ≤ y ≤ exp(log x/(log log x)⁴), Q ≤ y^{c₂/log₃x}, χ
  Q-exceptional of modulus q ≤ Q and |η|x ≤ x^{1/4}:
  W(x, y; q, η) ≪ [2^{ω(q/q₁)}√q₁Ψ(x, y)/φ(q)](q/q₁)^{1−α}·(1/((1 + |η|x)^{α+β₁−1}x^{1−β₁}H(u)^{δ}) + (1 + |η|x)R),
  R := y^{−c₂/log₃x} + e^{−c₂√(log x)} + (log x)^{3/2}x^{−c₂/log Q}. Its proof takes the height T′ := min{y^{c₂/log₃x},
  e^{√(log x)}}, so "Q-exceptional" is for that height. With q = q₁ and η = 0, W(x, y; q₁, 0) = τ(χ₁)Ψ(x, y; χ₁)/φ(q₁)
  and |τ(χ₁)| = √q₁, so the prefactor cancels and |Ψ(x, y; χ₁)| ≪ Ψ(x, y)(H(u)^{−δ} + R) (drop x^{β₁−1} ≤ 1).
- **Dependencies.** Props 2.9 and 2.13 have complete written proofs (pp. 609–628), resting on Lemmas 2.2 and 2.4–2.8
  (proved in the paper; Lemma 2.6 credits Harper's preprint but is re-proved in full "par souci de complétude") and on
  the log-free zero-density estimate (2.16), cited from Huxley 1974 and Jutila 1977 (published). The remarks citing
  [Har12b] enter no proof. Siegel's theorem is not used, so the constants appear effective. Harper 2016 still needs
  §3.3 of his preprint because his major arcs need the range y ≤ x and Q up to min{y^{c₂}, e^{c₂√log x}}, beyond
  Prop 2.13's y ≤ exp(log x/(log log x)⁴) and Q ≤ y^{c₂/log₃x}; there is no contradiction.

**Lemma SC-D [proved, from [D-CJM] Props 2.9, 2.13, (1.2) and Hildebrand–Tenenbaum Thm 1].** Hyp SC(κ) holds for every
fixed κ < 1 (on its domain: (log x)^{1/2} ≤ log y, y ≤ x^{1/100}, u ≥ (log log x)⁵, x ≥ x₀(κ)).
*Proof (v2; repaired after the literature critic and codex mf2L_2).* Fix z′ ∈ [x^{1/3}, x] and work at the scale z′
(Prop 2.9 with x := z := z′). Write L := log z′, Y := log y and u′ := L/Y ∈ [u/3, u]. On the domain,
Y ≥ (log x)^{1/2} ≥ √L, so u′ ≤ √L and u′ ≤ Y; and Y = o(L) because u′ → ∞. Put
  Q := R₁ := e^{u^κ},  T′ := min{y^{c₂/log₃z′}, e^{√L}},  so log T′ = min(c₂Y/log₃z′, √L) ≥ c₂u′/log₃z′.
Hence log Q ≤ (3u′)^κ = o(log T′), so log(2QT′) ≤ 2 log T′, Q ≤ y^{c₂/log₃z′} ≤ y^{c₂} and T′ ≤ y^{c₂}. Use (1.2) at
height 2T′ (Lemma 2.7 needs the region up to |τ| ≤ 2T) with Drappeau's b. If χ₁ is exceptional for (Q, 2T′), its zero
β₁ ≥ 1 − b/log(2QT′) > 1 − b/log(QT′) also lies in the height-T′ region, so by uniqueness χ₁ is also the Q-exceptional
character for height T′ in Prop 2.13's convention. Let χ be primitive with 1 < cond χ ≤ Q.
- *χ is not Q-exceptional, or χ = χ₁ with β₁ ≤ 1 − √c₂/log(2QT′).* Then L(s, χ) ≠ 0 on σ > 1 − ε, |τ| ≤ 2T′, with
  ε := min(b, √c₂)/log(2QT′) ≍ 1/log T′. Lemma 2.7's side conditions:
  (i) ε log y ≍ Y/log T′ = max{log₃z′/c₂, Y/√L} ≥ log₃z′/c₂ → ∞ (codex mf2L_3);
  (ii) T′ ≥ y^{c₁ε}(log x)²: c₁εY ≍ Y/log T′, and Y = o((log T′)²) because (log T′)² = min(c₂²Y²/(log₃z′)², L) and both
  Y²/(log₃)² and L are ≫ Y; so c₁εY + 2 log log x = o(log T′).
  Prop 2.9 (γ = 0) gives |Ψ(z′, y; χ)| ≪ z′^αζ(α, y)((log T′)z′^{−c₂ε} + T′^{−c₂}), where z′^{−c₂ε} = e^{−cL/log T′} ≤ e^{−c√L}
  ≤ e^{−cu′} (as log T′ ≤ √L) and T′^{−c₂} ≤ e^{−c₂²u′/log₃z′}. By Hildebrand–Tenenbaum Thm 1,
  z′^αζ(α, y) ≍ Ψ(z′, y)α√(2πφ₂) ≪ Ψ(z′, y)(log x)². Since u^κ = o(u/log₃x), the bound is ≤ e^{−u^κ}Ψ(z′, y) for large x.
- *χ = χ₁ Q-exceptional with β₁ > 1 − √c₂/log(2QT′).* Prop 2.13 at (z′, y), with q = q₁ ≤ Q and η = 0, gives
  |Ψ(z′, y; χ₁)| ≪ Ψ(z′, y)(H(u′)^{−δ} + R) (prefactor cancelled as above). Its hypotheses: (log z′)^{c₁} ≤ y;
  y ≤ exp(L/(log L)⁴) since u′ ≥ u/3 ≥ (log log x)⁵/3; Q ≤ y^{c₂/log₃z′}. Each term is ≤ e^{−2u^κ}:
  H(u′)^{−δ} = e^{−δu′/log²u′}; y^{−c₂/log₃z′} = e^{−c₂Y/log₃z′} with Y ≥ u′; e^{−c₂√L} with √L ≥ u′; and
  (log x)^{3/2}z′^{−c₂/log Q} = (log x)^{3/2}e^{−c₂L/u^κ} with L/u^κ ≥ u′^{1−κ}Y/3 ≥ u′/3. ∎

So **Hyp SC is a consequence of published theorems with full proofs.** Since [D-CJM] does not use Siegel's theorem,
x₀ appears effective (not re-derived here).

### 4.5 A second route for large conductors (codex), avoiding Harper 2016 [sketch]

Codex's energy lemma (`codex/mf2L_1.out` (11)) replaces Lemma LV by: pointwise bounds for conductors up to
R_* = e^{cu/log²u}; for R_* < r ≤ R₀ := y^{c₀}, the log-free zero-density estimate (2.16) bounds the number of characters
without a zero-free region by O(R^{1/10}), and the others are handled by Prop 2.9; for r > R₀, the classical
multiplicative large sieve (Iwaniec–Kowalski Thm 7.13: Σ_{q≤Q}(q/φ(q))Σ*|Σa_nχ(n)|² ≤ (N + Q² − 1)Σ|a_n|²) costs a
factor 1/(ρR). This needs the density restriction ρ ≥ y^{−c₀/2}, i.e. u log u ≤ ε₀ log y at every library scale, which
holds once the constant a in B = L[α, a] is large (child scales: log(1/ρ)/log B = O(1/a²); algebraic parents:
u₂ log u₂/log B₂ = O(D^{4−11α}) → 0). Enlarging a changes c(α), not α. Caveats (codex mf2L_2): Lemma SC-D at *fixed*
κ does not reach R_*, so route B needs the quantitative small-conductor consequence of the same proof (available: the
exceptional term only gives e^{−cu/log²u}, so R_* = e^{cu/log²u} with c small, which must be matched against the
zero-density range); and the full energy bound must keep the large-sieve term Q/(xρ) and the induced-character
prefactor. Codex judges these harmless in our families. Route A needs no condition on a; route B needs no Harper.

## 5. Theorems 2L-SC(α) and 2L-MF(α)

**Theorem 2L-SC(α) [proved conditional on Hyp SC(κ) for some κ ∈ (κ(α), 1), on Fact H4.1, and on the internal chain
as in 2L-final].** Let α ∈ (3/8, 2/5) be fixed, or α = 3/8 with a > 2c₂. Then FACTOR_2L in design 2L_DW(α), with
value factorisations by FACTOR_MF, is Las Vegas and factors every N in expected time L_N[α, c(α)], where c(α) is the
constant of Theorem 2L-final.
*Proof.* Theorem 2L-final's proof (thmAB §12, c12alg §5, c12rat §9) uses Drappeau's Lemme 5 only in items 8–11 and 12
of §1. In items 8–11 it enters only through (WA-hyp) with ε = ε_lib (c12rat Cor WA′, thmAB Lemma FAM-A). Prop SUB
supplies (WA-hyp) with ε𝔇̃ = o(T₀^{−1}e^{−D^θλ²}), and the parent analogue. Every later use of the character total is
additive: an absorption into T₀^{−1} (thmAB §12(i)), into the parent budget (thmAB §12(ii)), or into the additive
term of RAT-1/RAT-2 (c12rat §8.3). The weaker bound does each of these jobs, because it is ≪ D^{−2α} and C_R is chosen
after C_T. Item 12: by Lemma OS the law of the outer transcript is the same for every oracle that returns canonical
factorisations, and by Lemma W the cost is the oracle's expected cost. FACTOR_MF costs
L_{X_F}[2/5, c_MF + o(1)] = e^{O(D^{(2−3α)·2/5}λ)} = L_N[α]^{o(1)} on values ≤ X_F, because (2 − 3α)·2/5 < α for
α > 4/11. The cost constant is unchanged, since no character term enters the ledger. ∎

**Theorem 2L-MF(α) [proved relative to the internal chain].** Under the hypotheses of Theorem 2L-SC(α) with Hyp SC
removed, the same conclusion holds, with only published external inputs (§6).
*Proof.* Lemma SC-D gives Hyp SC(κ) for every κ < 1 on its domain, which contains every endpoint of all four
families (`params`, SC-D regime; codex mf2L_2 checked the subsidiary endpoints after induction, Möbius and J-division). ∎
*Variant [not proved; sketch].* §4.5's route B in place of Fact H4.1 should give the same for a ≥ a₀(α, d₀, κ, b₂),
with c(α) evaluated at that a. It needs the quantitative small-conductor range R_* = e^{cu/log²u} (the exceptional term
gives e^{−cu/log²u}, not e^{−cu/log₃x}) and the complete energy bound; codex mf2L_3 judges it incomplete as written.

## 6. Published inputs of Theorem 2L-MF(α)

Those of 2L-final (c12alg §5, c12rat §9.1, thmAB §14, uprime §10, lv5unif §14) **minus** Drappeau's Lemme 5,
Pascadi's Lemma 3.2, Harper's arXiv:1208.5992 and Lemma 8.5′/8.5H (via FACTOR_U); **plus**:

| Input | Statement used | Source | Checked |
|---|---|---|---|
| Fact H4.1 | large-values bound (4.1) for 1/x-spaced points, any weights |a_n| ≤ 1 on smooth n | Harper, Compositio Math. 152 (2016) 1121–1158, Prop. 2, (4.1) | accepted manuscript (lit. critic, codex); full proof, no 1208.5992 |
| [D-CJM] 2.9 | twisted smooth sum under a zero-free region, q, T ≤ y^{c₂} | Drappeau, Canad. J. Math. 67 (2015) 597–638, Prop. 2.9 (= arXiv 1302.4318 Prop 2.1) | published PDF (lit. critic, codex); full proof |
| [D-CJM] 2.13, (1.2) | exceptional character, Q ≤ y^{c₂/log₃x}; definition of Q-exceptional | ibid., Prop. 2.13 (= arXiv Prop 2.5), eq. (1.2) | published PDF; full proof |
| Landau–Page | at most one exceptional real primitive character below conductor Q, height T | Iwaniec–Kowalski Thms 5.26, 5.28; Davenport Ch. 14 (as used in [D-CJM] (1.2)) | IK numbering from secondary sources |
| log-free zero density | [D-CJM] (2.16), used inside Prop 2.9's proof | Huxley 1974 (1.1); Jutila 1977 (1.8) | cited by Drappeau |
| Gauss sums, Hölder, Möbius | textbook | — | — |
| (route B only) [D-CJM] (2.16); Iwaniec–Kowalski Thm 7.13 | log-free zero density; multiplicative large sieve | ibid.; IK, *Analytic Number Theory* (2004) | codex |
| Theorem U-MF's inputs (for FACTOR_MF) | HT 1986 Thms 1, 2(i), 3; CEP; Rosser–Schoenfeld; Stark Lemma 8; Landau; Lagarias–Odlyzko Thm 1.3; AKS; Cantor–Zassenhaus; Kaltofen–Saunders | notes/mfU_theorem.md §1.2 | there |

Still used from 2L-final: Hildebrand, JNT 22 (1986) Thms 1, 3; Hildebrand–Tenenbaum 1986 Thms 1–3; de la
Bretèche–Tenenbaum, Ramanujan J. 9 (2005) Thm 2.4(i); Montgomery–Vaughan's large sieve (L2-U′); Tenenbaum's Ψ_q
theorem; Alladi; the Hildebrand–Tenenbaum survey (JTNB 1993) Cor 2.4, Thm 5.1; Trudgian 2016; Lagarias–Odlyzko, Stark,
Landau (LV6); AKS, Kaltofen–Saunders, Cantor–Zassenhaus, Bernstein.

## 7. Open routes and remaining risk

0. **Route B** (§4.5) is incomplete. It is not needed, since route A suffices.
1. **The internal chain.** Theorem 2L-MF is only as good as thmAB, c12rat (v2 repairs), c12alg, uprime, pc3 and
   multidescent. The parallel agent `audit2L` audits that.
2. **Reduce the losses 𝔇 to polylogarithmic (proof critic).** This would make the needed small-conductor range
   (log x)^{O(1)} (a Siegel–Walfisz theorem for smooth numbers), and so make the external input much more robust. The
   losses come from the thinness threshold T₀, the filter spread (N^{η₀/4})^γ and the box spread (SPR), all of exponent
   D^θλ². This would mean rewriting thmAB §§8–10 and c12rat §§6–8.
3. **Pair-phase averaging** (Prop AVG). On the rational side ι = λ, so small-conductor characters could cancel over
   pairs. Codex's exact example (determinant 41 in [−5, 5]²: only two slopes occur) shows that fixed-determinant
   averaging alone does not equidistribute λ. It is not needed now.
4. **Below 3/8.** Nothing here changes multidescent §3.5: depth ≥ 3 needs precision below the heuristic D^{−2α}
   covariance floor. The input substitution is exponent-neutral.

## 8. Reviews

**Opus proof critic (read-only, on v0): ACCEPT WITH FIXES, no FATAL.** It checked the dependency map by grepping the
whole chain (complete as far as it could verify), Lemmas GS, LV and IND line by line, κ(α), Fact H4.1's hypotheses in all
four families, and the FACTOR_MF substitution. Probabilities: 2L-SC correct relative to 2L-final ≈ 85% after fix 1;
map complete 90%; GS/LV/IND/SUB 90%. Findings, all applied in v1:

| # | Sev. | Finding | Action |
|---|---|---|---|
| 1 | MINOR (breaks α ≥ 7/18 as written) | Prop SUB used the child T₀ in the algebraic-parent R₁, so R₁ ≤ e^{u₂^κ} fails for α ≥ 7/18 | R₁^p := e^{C_RD^{θ_p}λ²} with the parent threshold ϑ₀^p |
| 2 | MINOR | Lemma LV log exponents | corrected (4.05, 2.7; the stated (log z)³ survives) |
| 3 | MINOR | Lemme 5 ⇒ Hyp SC asserted, not proved | Remark in §4.2 |
| 4 | MINOR | Hyp SC not uniformly stated; effectivity | uniform statement + effectivity paragraph |
| 5 | MINOR | two definitions of R₁; c12rat's 𝔇 is unnormalised | R₁ := e^{C_RD^θλ²}, normalised 𝔇̃ |
| 6 | MINOR | classical LS starts at ρ^{−1}, not ρ^{−1/2} | corrected |
| 7 | MINOR | "rational side is a multiple of G" false as stated | replaced by the cost argument (several rational sides still leave same-side pairs) |
| 8 | MINOR | §3 stated only at 3/8; child count loose | restated on the window; m ≍ D^α |
| 9 | MINOR | Q ≤ z^{1/2001}; sign in an exponent | fixed |
| 10 | MINOR | Q5 unanswered | §4.4 (Lemma SC-D) |
| 11 | MINOR | framing: trading a by-reference lemma for a hypothesis | §0 caveat (now moot if Lemma SC-D stands) |

**Codex blind round `codex/mf2L_1` (xhigh, launched before v0 was written).** It found the same dependency map, Prop
INV, the d-th power obstruction, a periodic-set version of Prop NID, and the failure of MF for the per-field joint
event (its countermodel: one uniformly chosen dead child per accepted parent gives per-child failure 1/m but joint
success 0). It also found what we had missed: **Drappeau CJM 2015 Props 2.9 and 2.13** (with full proofs, checked
against the published PDF) as the small-conductor input, and an energy lemma (route B, §4.5) using the log-free
zero-density estimate (2.16) and the classical large sieve under a density restriction (a large). Its literature table:
Harper 2016 Appendix A.1 and Pascadi 2025 Lemma 3.2 still rely on 1208.5992; Harper JLMS 2025 Cor 1 uses the preprint's
Thm 1; Harper JNT 2012, Soundararajan 2008, Granville 1993 (I, II), Fouvry–Tenenbaum 1991/1996 and de la
Bretèche–Drappeau 2020 do not give the needed family estimate (ranges or error terms). Its assessment: "a
published-input analytic repair below 2/5 is realistic"; it did not certify the internal chain. It did not use
Harper's (4.1) (it read only Smooth Numbers Result 3 and Appendix A.1 of Harper 2016).

**Opus literature critic (read-only, on v1): both (SC) and (LC) available from published full proofs, with two fixes.**
- Drappeau CJM: Props 2.9 and 2.13 are FULL PROOF (pp. 609–628; Lemma 2.6 re-proves Harper's lemma; the [Har12b]
  remarks enter no proof; no Siegel). It reported the numbering correspondence, the exact statements and the hidden
  hypotheses of Prop 2.9 (Lemma 2.7: ε log y ≥ 1/c₂, T ≥ y^{c₁ε}(log x)², region to height 2T). It also explained why
  Harper 2016 still cites the preprint: wider ranges than Prop 2.13's.
- Harper (4.1): FULL PROOF, same numbering in the accepted manuscript. (4.2), Theorem 2 and the §5 applications rely on
  Appendix A.1, hence on the preprint.
- Required fixes, applied in v2: one height T′ = min{y^{c₂/log₃x}, e^{√log x}} throughout Lemma SC-D; height 2T′; the
  Hildebrand–Tenenbaum normalisation.
- Iwaniec–Kowalski numbering is UNVERIFIED (book not accessible). Soundararajan, Granville, Harper JNT 2012,
  Fouvry–Tenenbaum and de la Bretèche–Tenenbaum give at most o(1) savings or classical-range bounds, too weak for (SC).
  So [D-CJM] is the source to use.

**Codex adversarial round `codex/mf2L_2` (xhigh, on v1): MAJOR REVISION, no FATAL; 85% that the repaired substitution is
correct; "after the specified repairs, Theorem 2L-MF follows relative to the stipulated internal chain for every fixed
3/8 ≤ α < 2/5"** (route A with no extra condition on a; at 3/8, a > 2c₂). It re-verified [D-CJM] Props 2.9 and 2.13
(including the Gauss-sum prefactor cancellation) and (4.1) independently. All findings were applied in v2:

| # | Sev. | Finding | Action |
|---|---|---|---|
| 1–2 | MAJOR | Hyp SC's domain (log y ≥ (log x)^{1/3}) was larger than SC-D's; there "log y ≥ u" and the Remark fail | domain restricted to log y ≥ (log x)^{1/2}, u ≥ (log₂x)⁵; all families satisfy it at every endpoint |
| 3 | MAJOR | the J-library is not smooth-supported; IND cannot apply to it | IND on each fibre, then Minkowski in the energy norm (no #J) |
| 4 | MAJOR | normalisers and endpoint domains implicit | unfiltered EN normalisers Ψ(v, B), Σ_ℓΨ(v/ℓ, B); endpoint ranges stated |
| 5 | MAJOR | weaker saving is not a literal replacement of the old e^{−cD^α/λ²} bounds | explicit downstream absorption (Prop SUB), C_R > 3(C_𝔇 + C_T + 2); incidence budget rechecked incl. 3/8 |
| 6 | MAJOR | INV/NID do not prove that every design needs BV | verdict restricted to this architecture and MF2 |
| 7 | MAJOR | "m ≍ D^α", "η ≍ D^{−α}" inaccurate in the window | cap D^{2−4α}, averaged count u₂ ≍ D^{3−7α}, η ≳ D^{−(3−7α)}; rank obstruction kept (3 − 7α > 1 − 2α) |
| 8 | MINOR | Prop AVG needs χ̄(ι) (counterexample k = 5, ι = 2) | fixed; exact test added |
| 9 | MINOR | "vanishes exactly if"; "nontrivial" | fixed |
| 10 | MINOR | log log R at R = 2; Q ≤ z^{1/2002} | fixed |
| 11 | MINOR | Prop 2.13 is a W-bound; constants must be compatible | specialisation written out; Drappeau's b used |
| 12 | MINOR | route B needs the quantitative small-conductor range and the omitted terms | §4.5 marked sketch with these caveats |

**Codex round 3 `codex/mf2L_3` (high, on v2): ACCEPT WITH MINOR FIXES; 92% that Theorem 2L-MF(α) is correct relative
to the internal chain; "route A supports Theorem 2L-MF(α) … route B remains incomplete".** Ten of the twelve round-2
findings PASS, including the rewritten Lemma SC-D, which it checked against the published Props 2.9 and 2.13 and (1.2).
Three minor findings, applied in v2.1:
- the typical child count m ≍ u₂ is heuristic (Lemma η gives only W/μ₁ ≤ 3u₂);
- route B is incomplete, so the §5 variant is now labelled a sketch;
- SC-D condition (i) should read Y/log T′ = max{log₃z′/c₂, Y/√L}.

## 9. Reproduction

`.venv/bin/python -m pytest -q tests/test_mf2L.py` (7 tests, < 1 s); `.venv/bin/python -m src.mf2L.params`.
