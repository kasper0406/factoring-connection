# mf2L: transferring the U-MF trick to the two-level chain — results

Agent `mf2L`, phase 31, 2026-10-28. **v0** (pre-review; §8 records the reviews).

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

**Plain language.** The transfer does not go through as hoped. Theorem 2L-final (L[α] for α ∈ (3/8, 2/5), and
L[3/8]) still needs one estimate of Bombieri–Vinogradov type, and we could not get it from published theorems alone.
But we narrowed it down a lot. 

- **What already transferred.** The integer multiplier is already part of the two-level design (design DW). It makes
  every *one-point* estimate an interval density: level-1 first moments, the child means, the parent counts and U′.
  Lemma 8.5H is not needed. The one remaining use of Theorem U (to factor values recursively) can call the BV-free
  FACTOR_MF instead.
- **What cannot transfer.** The two-level proof needs each child to fail with probability below about 1/(number of
  children). So it needs pair covariances to *second order*, not just upper bounds. Two values of the same form at two
  points P₁, P₂ always satisfy F(P₁) ≡ λ^d F(P₂) (mod det(P₁, P₂)), for every integral form F. No multiplier,
  shift, second rational side or MF-style family of fields can change this (Prop INV). So pair counts live on a
  lattice. Interval densities do not determine them (Prop NID), and characters modulo the pair determinants enter
  exactly. The MF design (independent algebraic sides) does not help at level 2. Children must be descended in their
  parent's field, and a single parent indicator has no concentration to decouple against (§3).
- **What the BV input reduces to (new).** The four character-error families (rational and algebraic, children and
  parents) need only:
  - (a) for conductors above R₁ = e^{O(D^θλ²)}: a mean-square bound, which we derive from Harper's **published**
    large-values estimate (Compositio 152 (2016), Proposition 2, bound (4.1)). By Harper's own account, that bound does not
    use his unpublished preprint (Lemmas GS, LV);
  - (b) for conductors 1 < k ≤ R₁ = exp(u^{κ}): a pointwise bound |Ψ(z, y; χ)| ≤ e^{−u^κ}Ψ(z, y)
    (**Hyp SC(κ)**, §4.3), for a fixed κ = κ(α) ∈ [2/3, 1). This is a Siegel–Walfisz-type statement for smooth
    numbers, in a conductor range beyond the classical one, and it needs a treatment of a possible exceptional (Siegel)
    character.

  Harper's (4.1) bound with polylogarithmic loss closes the "λ-gap" between pointwise bounds (which reach only conductors
  e^{O(u)}) and the classical large sieve (which starts at ρ^{−1/2} = e^{(u log u)/2}).
- **Theorem 2L-SC(α) [proved conditional on Hyp SC(κ(α)) and on the internal chain; §5].** For every fixed
  α ∈ (3/8, 2/5), and at α = 3/8 under a > 2c₂: a Las Vegas algorithm factors every N in expected time L_N[α, c(α)].
  c(α) is the constant of 2L-final, unchanged. Its published inputs are listed in §6 and exclude Drappeau's Lemme 5 and
  Harper's arXiv:1208.5992. **We do not have Theorem 2L-MF**, i.e. 2L with only published inputs. That would need Hyp SC
  from a published source with a full proof. We found none, and codex found none (§4.4).
- **Minimal remaining input (Q5).** Hyp SC(κ): for primitive χ with 1 < cond χ ≤ exp(u^κ) and z in the ranges of §4.3,
  |Ψ(z, y; χ)| ≤ exp(−u^κ)Ψ(z, y). Here κ(3/8) = 2/3, and κ(α) → 1 as α → 2/5 (`params.kappa`). It is much weaker
  than Drappeau's Lemme 5 (pointwise, conductors x^{o(1)}, saving e^{−o(u)}). Its natural proof (zero-free regions for
  L(s, χ) at conductor y^{o(1)}, the Landau–Page exceptional character, a saddle-point bound) is the content of §3.3 of
  Harper's preprint, which is not published.

**Probabilities (ours, before review).**
- The dependency map (§1) is complete and correct: 80%.
- Prop INV, AVG and NID, and the d-th power remark: 97% (elementary).
- The substitution (§4: Lemmas GS, LV, IND, Prop SUB) is correct given Harper's (4.1) as quoted: 80%.
- (4.1) is BV-free and stated as quoted in the published version: 75% (v1 read; published text not yet seen by us).
- Theorem 2L-SC(α) is correct relative to the internal chain and 2L-final: 55% (2L-final's own 0.65 times the above).
- Hyp SC is provable with known methods: 85%. Hyp SC is already in a published paper with full proof: 15%.

---

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
- *Rational side (e = 1).* Whatever the law of g, and whatever further randomisation of the rational side keeps it a
  linear form, (G(P₁), G(P₂)) lies on {n₁ ≡ λn₂ (mod Δ)}. A rational side with root m mod N and coefficients
  ≪ N^{1/2} is an integer multiple of G (the lattice {(u, v) : u ≡ mv (mod N)} has covolume N). A multiple of G only
  rescales the lattice condition. A shift G + N·(·) makes the values ≈ N, which destroys smoothness. Neither the
  "random shift", nor the "second multiplier", nor "randomising g" makes rational residues mod Δ uniform. The
  residue *relation* between the two coordinates is the obstruction, not the residue of either coordinate.
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
  "F_j(P₀) accepted and all m children of P₀ in field j good". It concerns ≈ d^{3/2} > d + 1 points of one form, so no
  rank argument makes the children independent (Lemma OD / Prop R3).
- a lower bound for that joint event goes either through a union bound (precision η ≍ 1/m again, Lemma MP), or
  through a joint estimate (LS_C). The union bound must decouple the parent from its children. pc3's Prop DC does that
  with the parent **count** Z_𝔮, which concentrates because many parents share a label. In MF2 the parent of (g, ℓ) is
  one indicator 1[acc(P₀) ∧ q | F(P₀)]. Its relative standard deviation is ≈ ρ₂^{−1/2} = L^{Θ(1)} at 3/8, so
  Cauchy–Schwarz decoupling fails. Summing over ℓ restores a parent count and with it C₁ and C₂, with the same
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
#Θ_j ≤ C2^{j(2.75)}(log z)^{3.3}, since 20γ ≤ 1/2 and 32γ ≤ 0.8. Hence
Σ_{θ∈Θ}|S(θ)|³ ≤ Σ_j(2^{1−j}Ψ)³#Θ_j + #Θ·z^{−3/1000}Ψ³ ≤ C(log z)^{3.3}Ψ³ + Ψ³.
By Hölder, Σ_{θ∈Θ}|S(θ)|² ≤ (#Θ)^{1/3}(Σ|S|³)^{2/3} ≤ CR^{2/3}(log z)^{2.2}Ψ². Lemma GS, summed over R < k ≤ 2R (each
b/k occurs once), bounds Σ_kΣ*_χ|â(χ)|² by the same quantity. Finally φ(k) ≥ k/(C log log k). ∎

The exponent 1/3 is not optimal (any p ∈ (2 + 20γ, 4) works). Only a power saving in R with a polylogarithmic loss
is needed. The essential point is that the loss is **polylogarithmic**, not 1/ρ as in the classical large sieve: that
closes the λ-gap (§0).

### 4.2 From primitive characters to the WA hypothesis

**Lemma IND (imprimitive characters) [proved].** Let u be one of the library sequences of §4 with y-smooth support
and endpoint z, Q ≤ z^{1/2000}, and put T(z′) := Σ_{1<k₀≤Q}φ(k₀)^{−1}Σ*_{χ mod k₀}|Ψ_u(z′; χ)|²/U(z′)², where
Ψ_u(z′; χ) := Σ_{m≤z′}u(m)χ(m) is the same kind of sequence with endpoint z′. Then
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀ (k)}|û(χ)|² ≤ Cλ^{C}(log Q)·U(z)²·sup_{z/Q≤z′≤z}T(z′).
*Proof.* A non-principal χ mod k is induced by a unique primitive χ* mod k₀, k₀ | k, k₀ > 1. Let k₁ be the product of
the primes dividing k but not k₀. Then û(χ) = Σ_{e|k₁}μ(e)χ*(e)V_e with V_e := Σ_mu(em)χ*(m). If e is not y-smooth,
V_e = 0. Otherwise V_e is the same kind of sequence with endpoint z/e (P₀-freeness: V_e = 0 if e has a prime ≤ P₀).
Cauchy–Schwarz with σ := 1 − 1/λ gives |û(χ)|² ≤ ∏_{p|k}(1 + p^{−σ})·Σ_{e|k₁}e^{σ}|V_e|². Here
∏_{p|k}(1 + p^{−σ}) ≤ λ^C: the product is largest when the primes are the first ω(k) ≤ log Q primes, and
p^{1/λ} ≤ e for p ≤ D. Next, Σ_{k≤Q: k₀e|k}φ(k)^{−1} ≤ C log Q/(φ(k₀)φ(e)), because (k₀, e) = 1. So the left side is
≤ λ^C log Q·Σ_e(e^{σ}/φ(e))Σ_{k₀}φ(k₀)^{−1}Σ*_{χ*}|V_e(χ*)|². This is ≤ λ^C log Q·Σ_e(e^{σ}/φ(e))U(z/e)²T(z/e). Fact dBT1
gives U(z/e) ≤ Ce^{−β}U(z) up to (log P₀)^C. Then Σ_e e^{σ−2β}/φ(e) ≤ ∏_p(1 + Cp^{−2+2γ+1/λ}) = O(1). ∎
For the parent library 𝒜^J, |Σ_ℓχ(ℓ)A_{z/ℓ}(χ)|² ≤ #J·Σ_ℓ|A_{z/ℓ}(χ)|² and (Σ_ℓU(z/ℓ))² ≥ #J Σ_ℓU(z/ℓ)²/4 (U(z/ℓ)
varies by ≤ 2 over ℓ ∈ (Y, 2Y]). So the J-library reduces to its fibres with no loss.

**Hyp SC(κ) (small conductors) [hypothesis].** For the family's (x, y) (§4.3 lists them), every primitive χ with
1 < cond χ ≤ exp(u^κ), u := log x/log y, and every x^{1/3} ≤ z′ ≤ x: |Ψ(z′, y; χ)| ≤ exp(−u^κ)Ψ(z′, y).

**Prop SUB (substitution) [proved, given Fact H4.1 and Hyp SC(κ) with κ > κ(α)].** In each of the four families,
(WA-hyp) holds with ε ≤ Cλ^C(log x)^{4}(R₁^{−1/3} + R₁e^{−2u^κ}) for every 2 ≤ R₁ ≤ e^{u^κ}. With
R₁ := (𝔇T₀)⁶ (T₀ the thin-pair threshold e^{C_TD^θλ²}), ε𝔇 ≤ T₀^{−1}·(𝔇T₀)^{−1}λ^C(log x)^4 + o(T₀^{−1}).
*Proof.* Split the primitive conductors in T(z′) at R₁. For k₀ > R₁, sum Lemma LV over dyadic R ≥ R₁ (using
P₀-freeness and J-fibres as above; Fact H4.1 allows any weights |a_n| ≤ 1 on smooth n). This gives
Cλ^C(log x)^{3}R₁^{−1/3}. For k₀ ≤ R₁, Hyp SC gives Σ_{k₀≤R₁}φ(k₀)^{−1}·φ(k₀)e^{−2u^κ} ≤ R₁e^{−2u^κ}. For P₀-free
sequences, Möbius over d | P₀# turns Hyp SC into the same bound with a factor (log P₀)^C (as in c12rat Lemma EN; the
endpoints z′/d ≥ x^{1/3}). Lemma IND assembles the families. The choice of R₁ is admissible: log R₁ = O(D^θλ²) and
u^κ ≍ D^{ακ} with ακ > θ, and then R₁e^{−2u^κ}𝔇 = o(T₀^{−1}). Hypotheses of Fact H4.1: y ≤ x^{1/100}, 8Q² ≤ x^{1/1000},
1 − α ≤ 1/40, log^{21}x ≤ y (`params`, all families). ∎
For the algebraic parents, θ is replaced by θ_p and u by u₂ ≍ D^{3−7α}.

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

### 4.4 Is Hyp SC published? (Q5)

[To be completed from codex/mf2L_1 and the literature critic; §8.]

## 5. Theorem 2L-SC(α)

**Theorem 2L-SC(α) [proved conditional on Hyp SC(κ) for some κ ∈ (κ(α), 1), and on the internal chain as in
2L-final].** Let α ∈ (3/8, 2/5) be fixed, or α = 3/8 with a > 2c₂. Then FACTOR_2L in design 2L_DW(α), with value
factorisations by FACTOR_MF, is Las Vegas and factors every N in expected time L_N[α, c(α)], with the constant c(α)
of Theorem 2L-final.
*Proof.* Theorem 2L-final's proof (thmAB §12, c12alg §5, c12rat §9) uses Drappeau's Lemme 5 only in items 8–11 and 12
of §1. In 8–11 it enters only through (WA-hyp) with ε = ε_lib (c12rat Cor WA′, thmAB Lemma FAM-A). Prop SUB supplies
(WA-hyp) with ε𝔇 = o(T₀^{−1}). Every later use of the character total is an absorption into T₀^{−1} or into
e^{−cD^α/λ²}μ₁, and o(T₀^{−1})μ₁ ≪ D^{−2α}μ₁ does the same job (thmAB §12, c12rat §8.3: the character totals enter
only additively). Item 12: by Lemma OS the law of the outer transcript is the same for every oracle that returns
canonical factorisations, and by Lemma W the cost is the oracle's expected cost. FACTOR_MF (Theorem U-MF) costs
L_{X_F}[2/5, c_MF + o(1)] = e^{O(D^{(2−3α)·2/5}λ)} = L_N[α]^{o(1)} on values ≤ X_F, because (2 − 3α)·2/5 < α for
α > 4/11. The cost constant is unchanged: none of the character terms enters the ledger. ∎

**Theorem 2L-MF(α)** (the same, with only published inputs) is **not proved**. It follows from Theorem 2L-SC(α) if
Hyp SC(κ(α) + ε) is a published theorem with full proof.

## 6. Published inputs of Theorem 2L-SC(α)

Those of 2L-final (c12alg §5, c12rat §9.1, thmAB, uprime §10, lv5unif §14) **minus** Drappeau Lemme 5, Pascadi Lemma 3.2,
Harper arXiv:1208.5992, Lemma 8.5′/8.5H (via FACTOR_U), **plus**:
- Harper, Compositio Math. 152 (2016) 1121–1158, Proposition 2, bound (4.1) (arXiv:1408.1662v1 read; published
  numbering to be confirmed);
- Theorem U-MF's inputs for the oracle (notes/mfU_theorem.md §1.2: Hildebrand–Tenenbaum 1986 Thms 1, 2(i), 3; CEP;
  Rosser–Schoenfeld; Stark Lemma 8; Landau; Lagarias–Odlyzko Thm 1.3; AKS; Cantor–Zassenhaus; Kaltofen–Saunders);
- Gauss sums (textbook), Hölder.
Still used from 2L-final: Hildebrand, JNT 22 (1986) Thms 1, 3; Hildebrand–Tenenbaum 1986; de la Bretèche–Tenenbaum,
Ramanujan J. 9 (2005) Thm 2.4(i); Montgomery–Vaughan large sieve (L2-U′); Tenenbaum's Ψ_q theorem; Alladi; the
Hildebrand–Tenenbaum survey (JTNB 1993) Cor 2.4, Thm 5.1; Trudgian 2016.

## 7. Open routes

1. **Prove Hyp SC.** A self-contained proof would need zero-free regions, Landau–Page, and a saddle-point bound for
   twisted smooth sums at conductor y^{o(1)} with saving e^{−u^κ}. It is "routine" for an analytic number theorist,
   but it is real work, and it is exactly the unpublished part of Harper's preprint.
2. **Reduce the losses 𝔇 to polylogarithmic.** Then R₁ = (log x)^{O(1)}, and Hyp SC would shrink to a
   Siegel–Walfisz theorem for smooth numbers (conductors ≤ (log x)^A, possibly ineffective). The losses come from the
   thinness threshold T₀, the filter spread (N^{η₀/4})^γ and the box spread (SPR), all exponents D^θλ². This would
   mean rewriting thmAB §§8–10 and c12rat §§6–8.
3. **Pair-phase averaging** for small conductors (Prop AVG). On the rational side ι = λ, so every small-conductor
   character could in principle cancel over the pair family; on the algebraic side only characters of order ∤ d can,
   and real characters always can (d odd). This needs equidistribution of λ mod k₀ over pairs with given cells, which
   is not established, and the characters of order | d would remain.

## 8. Reviews

[pending]

## 9. Reproduction

`.venv/bin/python -m pytest -q tests/test_mf2L.py` (6 tests, < 1 s); `.venv/bin/python -m src.mf2L.params`.
