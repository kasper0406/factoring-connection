# mf2L: transferring the U-MF trick to the two-level chain — results

Agent `mf2L`, phase 31, 2026-10-28. **v1**: v0 + proof-critic fixes + codex blind round's Drappeau CJM route (§4.4–4.5); §8 records the reviews.

Files:
- this note;
- code: `src/mf2L/params.py` (exact-rational exponent margins of the input substitution, 16 checks + informational
  comparisons), `src/mf2L/checks.py` (exact checks of the finite statements: Props INV, AVG, NID, the d-th power
  kernel, Lemma GS in exact cyclotomic arithmetic);
- tests: `tests/test_mf2L.py` (6 tests, exact arithmetic);
- codex: `codex/mf2L_1.{md,out}` (blind round, launched before anything here was written), `codex/mf2L_2.*`
  (adversarial round on v0).
- No experiments were run (every finite claim is an exact identity, checked by the tests), so there is no
  pre-registration.

Labels as in thmAB/c12rat: **proved** (internal, unrefereed, from the cited published results), **sketch**,
**analysis** (an argument about a proof architecture, not a theorem about factoring), **open**, **UNVERIFIED**.
Notation: λ := log log N, D := log N/λ, L[α, c] = e^{(c+o(1))D^αλ}. Design 2L_DW(α) of c12alg §1.1, α ∈ [3/8, 2/5).
Ψ(z, y; χ) := Σ_{n≤z, P⁺(n)≤y}χ(n); u = log z/log y.

---

## 0. Verdict

**Plain language.** The U-MF trick itself does not carry over, but the goal is reached by a different route.
Theorem 2L (rigorous L[α] for every α ∈ (3/8, 2/5), and L[3/8] under the same constraint as before) can be made to rest
on **published theorems with full proofs only**. This assumes the project's internal chain, which is unrefereed as
before, and two literature facts that are being double-checked (§8).
- **What the integer multiplier does.** It is already part of the two-level design (design DW). It turns every
  *one-point* estimate into an interval density: level-1 first moments, child means, parent counts, and U′. So
  Lemma 8.5H is not needed. The one remaining use of Theorem U (factoring values recursively) can call the BV-free
  FACTOR_MF instead.
- **What it cannot do.** The two-level proof needs each child to fail with probability below about 1/(number of
  children), so it needs pair covariances to *second order*. For every integral form F and points P₁, P₂,
  F(P₁) ≡ λ^dF(P₂) (mod det(P₁, P₂)). No multiplier, shift, extra rational side or MF family of fields changes this
  (Prop INV). Interval densities do not determine pair counts on such lattices (Prop NID). Characters modulo pair
  determinants therefore enter unavoidably, and MF's independent algebraic sides do not help at level 2 (§3). **This
  is where a BV-type input is unavoidable.**
- **But the BV-type input needed is weak, and published.** The four character-error families need
  (WA-hyp): Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀}|û(χ)|² ≤ εU², with ε ≈ e^{−CD^θλ²}, not e^{−cD^α/λ²}. We prove it from:
  - (a) large conductors (> R₁ = e^{C_RD^θλ²}): Harper's **published** large-values bound (Compositio 152 (2016),
    Prop. 2, (4.1)), which does not use his preprint (Lemmas GS, LV, IND, Prop SUB). Alternatively (codex's route), the
    classical large sieve plus a log-free zero-density estimate, if the constant a in B is taken large (§4.5);
  - (b) small conductors: a pointwise bound **Hyp SC** (|Ψ(z, y; χ)| ≤ e^{−u^κ}Ψ for cond χ ≤ e^{u^κ}). This follows
    from **Drappeau, Canad. J. Math. 67 (2015), Props 2.9 and 2.13** together with Landau–Page (Lemma SC-D). Those
    propositions were found by codex's blind round, and they include the exceptional character.
- **Theorem 2L-MF(α) [proved relative to the internal chain, given the literature facts of §4.4; §5].** For every
  fixed α ∈ (3/8, 2/5), and at α = 3/8 with a > 2c₂, FACTOR_2L (design DW, value factorisations by FACTOR_MF) is Las
  Vegas and factors every N in expected time L_N[α, c(α)]. Drappeau's Lemme 5, Pascadi's Lemma 3.2 and Harper's
  arXiv:1208.5992 are not used. The exponent is that of 2L-final. With route (a) the constant c(α) is unchanged; with
  route (b) a must be large, which changes c(α).
- **Minimal BV-type statement (Q5)**, if one wants to isolate it: (WA-hyp) at the family's moduli, or, more primitively,
  Hyp SC(κ) for some κ ∈ (κ(α), 1) plus any mean-square large-conductor bound with polylogarithmic loss. Here
  κ(3/8) = 2/3 and κ → 1 as α → 2/5. Lemme 5 implies both, and so do Drappeau CJM + Harper (4.1).

**Caveats.**
- Everything internal (thmAB, c12rat, c12alg, uprime, pc3, multidescent, U-MF) remains unrefereed. 2L-final's own
  reviewers put it at ≈ 0.65 (window) and ≈ 0.52 (3/8). This note changes only the external inputs.
- The literature facts on which the "published only" claim rests are: Drappeau CJM Props 2.9 and 2.13 (statements,
  the notion "R-exceptional", and no dependence on 1208.5992), and Harper (4.1)'s independence of the preprint in the
  published version. §8 records the literature critic's verdict on them.

**Probabilities (ours, before the final reviews; §8 has the reviewers').**
- Dependency map complete and correct: 85%. [proof critic: 90%]
- Props INV, AVG, NID and the d-th power remark: 97%.
- Lemmas GS, LV, IND and Prop SUB, given (4.1) as quoted: 88%. [proof critic: 90% after fixes]
- Lemma SC-D, given Drappeau's statements as quoted: 85%.
- The literature facts (Drappeau CJM full proofs; (4.1) BV-free): 70%.
- Theorem 2L-MF(α) correct relative to 2L-final's internal chain: ≈ 60% (window), ≈ 55% at 3/8.

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
by a(x₁)b(x₂) and x₂ is restricted to units. For a fixed ι the count is A*B*/φ(k) + φ(k)^{−1}Σ_{χ≠χ₀}χ(ι)â(χ)\overline{b̂(χ)}.
*Proof.* Each pair of units (x₁, x₂) satisfies x₁ ≡ ιx₂ for exactly one ι. The second statement is orthogonality.
∎ (`checks.iota_average`, six moduli.)
So the character term vanishes exactly if the pair's unit ι is equidistributed over the pair family. Otherwise it must
be bounded through the twisted sums â(χ), i.e. through the distribution of smooth numbers in residue classes.

**Remark (d-th powers) [proved].** On the algebraic side ι = λ^d. A character χ with χ^d = χ₀ has χ(ι) = 1 for every
pair, so averaging over pairs never cancels it. Such characters exist mod k iff gcd(d, φ(k)) > 1
(`checks.dth_power_kernel`). For quadratic χ and odd d, χ(λ^d) = χ(λ) (`checks.quadratic_odd_power`). So the
exceptional (real) character is the one character for which pair-averaging *could* help on both sides (§7, open
route).

**Prop NID (interval densities do not determine pair counts) [proved; model].** Let k ∈ {3, 4, 6}, χ the real
non-principal character mod k, s(n) := 1[(n, k) = 1] and t(n) := s(n)(1 + εχ(n)). On every interval I of length ≥ k,
|Σ_It − Σ_Is| = ε|Σ_Iχ| ≤ ε, so the two weightings have the same interval densities up to ε/|I|. But with
c_r := #{n ≡ r (mod k)} in the support, Σ_{x₁≡x₂ (k)}t(x₁)t(x₂) − Σ_{x₁≡x₂ (k)}s(x₁)s(x₂) = Σ_{r unit}c_r²(2εχ(r) + ε²),
which is ε²Σ_rc_r² > 0 when the c_r are equal on units. (`checks.nid_example`: k = 3, 4, 6, ε = 1/2.)
*Consequence.* A covariance bound on L′ cannot follow from interval densities (Hildebrand–Tenenbaum short-interval
asymptotics) alone. It needs information on smooth numbers in residue classes modulo the divisors of Δ′. Q3's
answer is therefore **no**. What *can* be proved from interval densities is the principal-character part, which is
what Lemmas CF/CP and OPV-A already do.

**Necessity of the precision (multidescent Lemma MP, pc3).** Any moment-based proof that a child fails with
probability < η needs moments accurate to relative error η ≍ D^{−α}λ^{−C}. So upper bounds with polylog loss
(Lemma DS-type, the U-MF rational device) cannot suffice at level 2.

## 3. Design changes (Q2, Q4) [analysis]

*MF at level 2.* Design MF2: shared g, J_f independent (w_j, c_j), first rationally good parent P₀(g, ℓ). Fields are
independent given g, so a per-field success probability p(g, ℓ) ≥ ρ₂·L^{−O(1)} would suffice, and losses
(log N)^{O(m)} from m ≍ D^α children would be search work, not root density. But:
- algebraic children of field j must be descended in field j (Cor EX), so per field the event is the **joint** event
  "F_j(P₀) accepted and all m ≍ D^α children of P₀ in field j good". It concerns ≥ m > d + 1 points of one form, so no
  rank argument makes the children independent (Lemma OD / Prop R3).
- a lower bound for that joint event goes either through a union bound (precision η ≍ 1/m again, Lemma MP), or
  through a joint estimate (LS_C). The union bound must decouple the parent from its children. pc3's Prop DC does that
  with the parent **count** Z_𝔮, which concentrates because many parents share a label. In MF2 the parent of (g, ℓ) is
  one indicator 1[acc(P₀) ∧ q | F(P₀)]. Its relative standard deviation is ≈ ρ₂^{−1/2}, which is L^{Θ(1)} at 3/8 and
  e^{Θ(D^{3−7α}λ)} ≫ D^C for α ∈ (3/8, 2/5), so Cauchy–Schwarz decoupling fails on the whole range. Summing over ℓ restores a parent count and with it C₁ and C₂, with the same
  character families. The other option is conditioning on the exact value F(P₀), which fixes the residues of F(P′)
  modulo det(P₀, P′)/q (Prop INV, three points). That needs smooth numbers in one *prescribed* class modulo moduli of
  size L[1 − α]^{O(1)}, which is harder.
- conclusion: MF removes the level-1 algebraic second moment, which is not a BV input anyway (§1 item 3). It does not
  remove C₁ or C₂. This agrees with multidescent §4's table ("Prop MF fields sharing g … does not export a child").

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
*Status.* Harper deduces his Theorem 2 (restriction) from (4.1) **and** (4.2). Only (4.2) uses major-arc estimates,
which rest on §3.3 of arXiv:1208.5992 (his Appendix A.1). (4.1) is proved from his Theorem 1, Erdős–Turán,
Fouvry–Tenenbaum's exponential-sum bound and Bourgain's lemma (web-fetch of v1, this session; the published version is
to be checked by the literature critic, §8).

**Lemma LV (large conductors) [proved, given Fact H4.1].** Let 1 − α(z, y) ≤ 1/40, R ≥ 2 and 8R² ≤ z^{1/1000}. Then
  Σ_{R<k≤2R}φ(k)^{−1}Σ*_{χ mod k}|â(χ)|² ≤ C·(log log R)·R^{−1/3}(log z)^{3}·Ψ(z, y)².
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
and endpoint z, Q ≤ z^{1/2001} (so that Lemma LV applies at every endpoint z/e ≥ z/Q), and put T(z′) := Σ_{1<k₀≤Q}φ(k₀)^{−1}Σ*_{χ mod k₀}|Ψ_u(z′; χ)|²/U(z′)², where
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
For the parent library 𝒜^J, |Σ_ℓχ(ℓ)A_{z/ℓ}(χ)|² ≤ #J·Σ_ℓ|A_{z/ℓ}(χ)|² and (Σ_ℓU(z/ℓ))² ≥ #J Σ_ℓU(z/ℓ)²/4 (U(z/ℓ)
varies by ≤ 2 over ℓ ∈ (Y, 2Y]). So the J-library reduces to its fibres with no loss.

**Hyp SC(κ) (small conductors) [hypothesis].** There is x₀ such that for all x ≥ x₀ and y with
(log x)^{1/3} ≤ log y and y ≤ x^{1/100}, putting u := log x/log y: for every primitive χ with 1 < cond χ ≤ exp(u^κ)
and every x^{1/3} ≤ z′ ≤ x, |Ψ(z′, y; χ)| ≤ exp(−u^κ)Ψ(z′, y).
*Effectivity.* Theorem 2L-SC's conclusion (an expected-time bound) does not need x₀ to be effective: the algorithm uses
no analytic constant (alternation with trial division, as in Theorem U-MF), so an ineffective x₀ only makes the o(1)
in L_N[α, c(α) + o(1)] ineffective. If Hyp SC is proved via Siegel's theorem for a possible exceptional character, this
is what happens. An effective version is needed only for an effective N₀.

**Remark (Lemme 5 ⇒ Hyp SC(κ) for every κ < 1) [proved; proof critic #3].** Keep only the term k = k₀ = cond χ in
Drappeau's sum Σ_{k≤Q}φ(k)^{−1}Σ_{χ mod k, 1<cond χ≤x^{η_D}}|Ψ(x, y; χ)| ≪ Ψ(x, y)(H(u)^{−δ_D} + y^{−δ_D}). This gives
|Ψ(z′, y; χ)| ≤ Cφ(k₀)Ψ(z′, y)(H(u′)^{−δ_D} + y^{−δ_D}) at z′ (u′ ≥ u/3). For k₀ ≤ e^{u^κ} and κ < 1 this is
≤ e^{−u^κ}Ψ(z′, y) for large x, since H(u′)^{δ_D} = e^{δ_Du′/log²u′} ≫ e^{2u^κ}, and log y ≥ (log x)^{1/3} ≥ u^{1/2}·(log y)^{1/2}
gives log y ≥ u, so y^{−δ_D} ≤ e^{−δ_Du} ≤ e^{−2u^κ}. So Theorem 2L-SC contains
Theorem 2L-final.

**Prop SUB (substitution) [proved, given Fact H4.1 and Hyp SC(κ) with κ > κ(α)].** In each of the four families,
(WA-hyp) holds with ε ≤ Cλ^C(log x)^{4}(R₁^{−1/3} + R₁e^{−2u^κ}) for every 2 ≤ R₁ ≤ e^{u^κ}. Let 𝔇̃ := 𝔇/Σŵ be
the *normalised* family constant (c12rat's 𝔇 contains Σŵ ≍ μ₁·polylog; thmAB's is already normalised), with
log 𝔇̃ ≤ C_𝔇D^{θ}λ² (children, rational parents) or C_𝔇D^{θ_p}λ² (algebraic parents). Put
R₁ := e^{C_R D^θλ²} with C_R := 6(C_𝔇 + C_T) chosen after C_T and C_𝔇 (children and rational parents; T₀ = e^{C_TD^θλ²}),
and R₁^p := e^{C_R D^{θ_p}λ²} with the **parent** threshold ϑ₀^p = e^{C_TD^{θ_p}λ²} in place of T₀ (algebraic parents;
proof critic #1: the child T₀ would make R₁ ≤ e^{u₂^κ} impossible for α ≥ 7/18). Then ε𝔇̃ ≤ T₀^{−1}·e^{−D^θλ²}
(resp. ≤ (ϑ₀^p)^{−1}e^{−D^{θ_p}λ²}), which thmAB §12 absorbs into T₀^{−1} (resp. into the parent budget, §11/§12(ii)) and
c12rat §8.3 into its additive character term.
*Proof.* Split the primitive conductors in T(z′) at R₁. For k₀ > R₁, sum Lemma LV over dyadic R ≥ R₁ (using
P₀-freeness and J-fibres as above; Fact H4.1 allows any weights |a_n| ≤ 1 on smooth n). This gives
Cλ^C(log x)^{3}R₁^{−1/3}. For k₀ ≤ R₁, Hyp SC gives Σ_{k₀≤R₁}φ(k₀)^{−1}·φ(k₀)e^{−2u^κ} ≤ R₁e^{−2u^κ}. For P₀-free
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
y = B. In every case log y ≥ (log x)^{1/3} and y ≤ x^{1/100}.

**What Hyp SC is.** The conductors exp(u^κ) are y^{o(1)} but far beyond the Siegel–Walfisz range exp(c√log y)
(`params.informational`; exception: the algebraic parents inside the open window). The saving e^{−u^κ} is far weaker
than the truth (≈ e^{−cu} for non-exceptional χ). The standard proof is a zero-free region for L(s, χ) at conductor
y^{o(1)}, with at most one exceptional real χ₁ (Landau–Page), plus a Rankin / saddle-point bound for the twisted Euler
product ∏_{p≤y}(1 − χ(p)p^{−s})^{−1}. Heuristically an exceptional zero *helps* (χ₁(p) = −1 for most small p makes the
twisted sum smaller), but this must be proved. This is exactly what §3.3 of Harper's arXiv:1208.5992 does (and why
Lemme 5 and Harper 2016's Appendix A.1 cite it).

### 4.4 Hyp SC from published theorems (Q5): Lemma SC-D

Codex's blind round (`codex/mf2L_1.out` §5.1) pointed to a published source with full proofs: S. Drappeau, *Sommes
friables d'exponentielles et applications*, Canad. J. Math. 67 (2015) 597–638 ("[D-CJM]"), Propositions 2.9 and 2.13
(published numbering). Statements as reported by codex, to be confirmed by the literature critic (§8):
- **[D-CJM] Prop 2.9.** For q, T ≤ y^c, T ≥ 4 and (log x)^C ≤ y ≤ x, a non-principal χ mod q whose L-function has no
  zero in Re s > 1 − ε, |Im s| ≤ T satisfies |Ψ(x, y; χ)| ≪ x^βζ(β, y)((log T)x^{−cε} + T^{−c}).
- **[D-CJM] Prop 2.13.** For a primitive R-exceptional character of conductor ≤ R, y ≤ exp(log x/(log log x)⁴) and
  R ≤ y^{c/log₃x}: |Ψ(x, y; χ)| ≪ Ψ(x, y)(H(u)^{−c} + 𝓡), 𝓡 := y^{−c/log₃x} + e^{−c√(log x)} + (log x)^{3/2}x^{−c/log R}.
- **Landau–Page** (Davenport, *Multiplicative Number Theory*, Ch. 14; exact reference in §8): for some absolute c_L,
  among the primitive characters of conductor ≤ R, at most one (real) χ₁ has L(s, χ) = 0 with σ > 1 − c_L/log(RT),
  |t| ≤ T.

**Lemma SC-D [proved, given the three statements above].** Hyp SC(κ) holds for every fixed κ < 1 in the regime
(log x)^{1/2} ≤ log y, y ≤ x^{1/100}, u ≥ (log log x)⁵, x ≥ x₀(κ). (All four families are in this regime:
`params`, "SC-D regime".)
*Proof.* Let R₁ := e^{u^κ}, T := y^{c_T} with c_T ≤ c, χ primitive of conductor r ∈ (1, R₁], z′ ∈ [x^{1/3}, x],
u′ := log z′/log y ∈ [u/3, u]. Since u ≤ (log x)^{1/2} ≤ log y, we have log R₁ = u^κ = o(log y/log₃x), so R₁ ≤ y^{c/log₃x}
and R₁T ≤ y^{2c_T}.
- *χ is not the exceptional character for (R₁, T).* By Landau–Page, L(s, χ) ≠ 0 for σ > 1 − ε, |t| ≤ T, with
  ε := c_L/(2c_T log y). Prop 2.9 at (z′, y) gives |Ψ(z′, y; χ)| ≪ z′^βζ(β, y)((log T)e^{−cc_Lu′/(2c_T)} + y^{−cc_T}).
  By Hildebrand–Tenenbaum Thm 1, z′^βζ(β, y) ≤ (log x)^CΨ(z′, y). As u^κ = o(u) and log y ≥ u, the bracket times
  (log x)^C is ≤ e^{−u^κ} for large x.
- *χ = χ₁, with conductor ≤ R₁.* Prop 2.13 at (z′, y) with R = R₁ applies: u′ ≥ u/3 ≥ (log log z′)⁴ and
  R₁ ≤ y^{c/log₃z′}. Each error term is ≤ e^{−2u^κ}: H(u′)^{−c} = e^{−cu′/log²u′}; y^{−c/log₃x} = e^{−c log y/log₃x} and
  log y ≥ u; e^{−c√(log z′)} with √(log x) ≥ u; and z′^{−c/log R₁} = e^{−c log z′/u^κ} with log x/u^κ ≥ u^{1−κ}log y ≥ u. ∎

So **Hyp SC is a consequence of published theorems with full proofs**, *provided* the literature critic confirms (i) the
statements above, (ii) that "R-exceptional" in Prop 2.13 covers the Landau–Page exceptional character for height
T = y^{c_T}, and (iii) that the proofs of Props 2.9 and 2.13 do not use Harper's arXiv:1208.5992. Item (iii) matters
because Harper 2016, Appendix A.1, still cites §3.3 of the preprint for the exceptional-character term of his
major-arc estimate (in a different range).

### 4.5 A second route for large conductors (codex), avoiding Harper 2016

Codex's energy lemma (`codex/mf2L_1.out` (11)) replaces Lemma LV by: pointwise bounds up to R_* = e^{cu/log²u}
(Lemma SC-D); for R_* < r ≤ R₀ := y^{c₀}, [D-CJM]'s log-free zero-density estimate (its eq. (2.16)) bounds the number of
characters without a zero-free region by O(R^{1/10}), and the others are handled by Prop 2.9; for r > R₀, the classical
multiplicative large sieve (Iwaniec–Kowalski Thm 7.13) costs a factor 1/(ρR). This needs the *density restriction*
ρ = Ψ(x, y)/x ≥ y^{−c₀/2}, i.e. u log u ≤ ε₀ log y, at every library scale. That holds once the constant a in
B = L[α, a] is large (child scales: log(1/ρ)/log B = O(1/a²); algebraic parents: u₂ log u₂/log B₂ = O(D^{4−11α}) → 0).
Enlarging a changes c(α), not α. [sketch: codex's proof as written in its §5.1; not re-derived line by line here.] Route A
(Harper (4.1)) needs no condition on a. Route B needs no Harper. Either suffices with Lemma SC-D.

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

**Theorem 2L-MF(α) [proved relative to the internal chain and to the literature facts of §4.4].** Under the hypotheses
of Theorem 2L-SC(α) with Hyp SC removed, the same conclusion holds.
*Proof.* Lemma SC-D gives Hyp SC(κ) for every κ < 1 in all four families. ∎
*Variant.* Using §4.5's route in place of Fact H4.1 gives the same, provided a ≥ a₀(α, d₀, κ, b₂). Then c(α) is
the constant of 2L-final evaluated at that a.

## 6. Published inputs of Theorem 2L-MF(α)

Those of 2L-final (c12alg §5, c12rat §9.1, thmAB §14, uprime §10, lv5unif §14) **minus** Drappeau's Lemme 5,
Pascadi's Lemma 3.2, Harper's arXiv:1208.5992 and Lemma 8.5′/8.5H (via FACTOR_U); **plus**:

| Input | Statement used | Source | Checked |
|---|---|---|---|
| Fact H4.1 | large-values bound (4.1) for 1/x-spaced points, any weights |a_n| ≤ 1 on smooth n | Harper, Compositio Math. 152 (2016) 1121–1158, Prop. 2, (4.1) | arXiv v1 (web-fetch); published: §8 |
| [D-CJM] 2.9 | twisted smooth sum under a zero-free region, q, T ≤ y^c | Drappeau, Canad. J. Math. 67 (2015) 597–638, Prop. 2.9 | codex (published PDF); §8 |
| [D-CJM] 2.13 | exceptional character, R ≤ y^{c/log₃x} | ibid., Prop. 2.13 | codex (published PDF); §8 |
| Landau–Page | at most one exceptional real primitive character below conductor R, height T | Davenport, *Multiplicative Number Theory*, Ch. 14 | §8 |
| Gauss sums, Hölder, Möbius | textbook | — | — |
| (route B only) [D-CJM] (2.16); Iwaniec–Kowalski Thm 7.13 | log-free zero density; multiplicative large sieve | ibid.; IK, *Analytic Number Theory* (2004) | codex |
| Theorem U-MF's inputs (for FACTOR_MF) | HT 1986 Thms 1, 2(i), 3; CEP; Rosser–Schoenfeld; Stark Lemma 8; Landau; Lagarias–Odlyzko Thm 1.3; AKS; Cantor–Zassenhaus; Kaltofen–Saunders | notes/mfU_theorem.md §1.2 | there |

Still used from 2L-final: Hildebrand, JNT 22 (1986) Thms 1, 3; Hildebrand–Tenenbaum 1986 Thms 1–3; de la
Bretèche–Tenenbaum, Ramanujan J. 9 (2005) Thm 2.4(i); Montgomery–Vaughan's large sieve (L2-U′); Tenenbaum's Ψ_q
theorem; Alladi; the Hildebrand–Tenenbaum survey (JTNB 1993) Cor 2.4, Thm 5.1; Trudgian 2016; Lagarias–Odlyzko, Stark,
Landau (LV6); AKS, Kaltofen–Saunders, Cantor–Zassenhaus, Bernstein.

## 7. Open routes and remaining risk

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

**Opus literature critic:** [pending]

**Codex adversarial round `codex/mf2L_2` (on v1):** [pending]

## 9. Reproduction

`.venv/bin/python -m pytest -q tests/test_mf2L.py` (6 tests, < 1 s); `.venv/bin/python -m src.mf2L.params`.
