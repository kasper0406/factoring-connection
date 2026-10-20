# resbarrier: the generalized resultant barrier and its escapes (route 2)

Agent `resbarrier`, phase 16, 2026-10-20.
- Code: `src/resbarrier/forms.py` (exact binary-form tools, Lemma RES/C2/KL checkers),
  `src/resbarrier/ledger.py` (exact LPs, constants, exponent floors), `src/resbarrier/e1.py` (experiment E1).
- Tests: `tests/test_resbarrier.py`.
- Pre-registration: `notes/resbarrier_prereg.md` (commit `81d9f71`, before E1 ran). Raw data:
  `experiments/raw/resbarrier_E1.jsonl`.
- GPT-6.1-sol: blind round `codex/resbarrier_1.{md,out}` (§7); adversarial round `codex/resbarrier_2.{md,out}` (§8).
- Supplementary post-hoc data: `experiments/raw/resbarrier_E1_supp.jsonl` (§5).
- Literature: one research subagent (Kirchner 2018 read in full from mathnet.ru; Bean 1994; Thunder 2001;
  Stewart 2024; Schmidt 1993; Schnorr refutations). Reviews: two read-only Opus critics (§8).

Labels: **proved** (proof here, or a cited published theorem), **proved (model)** (an exact consequence of a stated
heuristic model), **heuristic** (named assumption), **KNOWN** (cited), **experimentally supported**, **open**.

---

## 0. Verdict

**Plain language.**
- No escape we examined gives a factoring exponent below 1/3, not even heuristically. This is within relation
  methods of the linked-order model (§1), for one N, without precomputation.
- The resultant wall is more general than "two polynomials". It is a statement about the *discriminant of the
  order the candidates live in*:
  - methods that extract factors through ring maps (all NFS variants) need that discriminant divisible by N²;
  - methods that use ramification (Kirchner-type class-group methods) need it divisible by N;
  - a classical area theorem (Bean 1994), combined with a lattice-point lemma (Davenport/Thunder), then caps how
    many candidates in a box can have small norms.
- This rank-2 barrier holds for every polynomial pair, for special N too, and for any selection inside a box
  (or special-q lattice) of side at most about the factor-base size B. It uses only the usual smoothness heuristic, with no counting heuristic for the polynomials.
- The counting heuristic itself ("G") is a theorem for almost all N, by an elementary union bound.
- One gap is not closed rigorously: families with ≥ 3 sieved coefficients and no rational side. There a
  counting-and-independence heuristic (G_r) closes it. A hypothetical family violating G_r would formally give
  L[1/4], but only under a further typical-size heuristic. Nobody knows how to build one, and G_r predicts
  that none exists for typical N.

**Precisely.** The model is single N, RAM, sparse linear algebra B², heuristic H ("an examined candidate's norm
is B-smooth as often as a random integer of the same size"), and box hypothesis H5. Batch NFS and the
factorization factory are outside it.

1. **Theorem E (elimination certificate) [proved; found blind by codex too].** For rings, commutative or not,
   with maps to Z/N agreeing on a common source, every nonzero integer in the sum of the kernels is a multiple
   of N. Res(f, g) is the one-variable instance. Without a sampling hypothesis there is *no* norm lower bound:
   unit powers, unipotents and Z[ε]/ε² have norm ±1.
2. **Lemma ND [proved, for extraction comparing two square roots of one algebraic element].** If
   gcd([O_A : R], N) = 1, then every split the relations can produce is already produced by evaluating the
   primitive idempotents, in polynomial time. Extraction through two different components needs N | Res, hence
   N² | disc. Ramification-based extraction needs N | disc, or a free gcd.
3. **Lemma C2 + Corollary RB2 [count: corollary of Bean 1994 Thm 1 + Thunder 2001 Lemma 14; cost: proved
   (model)].** #{box points with 0 < |F| ≤ X} ≤ 16·N^{−k/(D(D−1))}X^{2/D} + 2(D+1)(2E+1) when N^k | disc F.
   - With the yield integration of §2.3(e), α ≥ 1/3 and the constant is ≥ (16k/9)^{1/3}.
   - Scope: box or lattice side ≤ B^{1−o(1)} (H5), rank 2, D ≥ 3.
   - k = 2 (any NFS pair, any N): 1.5263, so SNFS is optimal among single-pair NFS *in this model*.
   - k = 1 (ramified): 1.2114, which is Kirchner's special-n constant.
   - Special-q / p-adic sublattices gain ≤ D²Q^{2/D} = L^{o(1)}.
   - E1 (pre-registered): 0 violations; SNFS/GNFS saturation ratio 59 (predicted ≈ 64).
4. **Theorem GAE [proved, elementary union bound; no prior statement found by us].** For all but N₀^{1−ε+o(1)} of
   the N ∈ [N₀, 2N₀), every polynomial pair (heights restricted as in §2.4) with ∏(2H_i+1)^{d_i+1} < N₀^{1−ε}
   has no common root mod N. The same holds for squarefree f with N | disc f.
   - So heurnfs's heuristic G holds for almost all N, at the level of coefficient heights. With full-box
     sampling (codex Thm 2), 1.9230 is then the single-pair floor [proved (model)].
   - The height part of Kirchner's "for most n" 1.5263 is necessary but not sufficient.
5. **Lemma KL (rank r with a rational side) [proved].** N ≤ C(d+r−1, r−1)·H(W₂)·(√(r−1)‖v‖^{1/(r−1)})^d. So
   the height floor is min{1, (r−1)/d}, and in the size model r = 2 is optimal.
6. **Rank r ≥ 3 with no rational side [open].**
   - 2-plane resultants give a family height ≥ x/d_max, *independent of r*. This permits the aggregate ledger
     x/d + (D/r)C. The floor 1/4 holds only under a typical-size heuristic.
   - G_r (Schmidt's counts of fixed-degree algebraic points plus independence of the r − 1 linking conditions)
     restores r = 2 optimality and α ≥ 1/3, and predicts that aggregate families are absent for typical N.

**Escapes (table §3).** None gives α < 1/3 in the forms we could analyse.
- MNFS changes only constants. The RES-only floor of *symmetric* MNFS is 1.71143, the known MexTNFS constant; mixed configurations reach the single-pair 1.5263.
- Multivariate and arithmetic-surface designs reduce to item 6.
- Non-commutative rings: structural obstructions, but not a proved exclusion.
- Other heights and p-adic smallness are sublattices.
- p-adic/adelic: no ring maps Z_p → Z/N.
- Kirchner: rank 2 with k = 1.
- Schnorr: D = 2 with |uw| ≥ N − 1; the lattice-shortness step fails (KNOWN refutations).
- Units, class groups, Stickelberger: no information without a DLP, a class group of discriminant ≥ N, or roots of
  unity mod p.

**Best mechanism.** None survives. The one conditional mechanism is the aggregate rank-r family (§4.2), which
needs a dependent-linking structure that nobody has exhibited.

## 1. The model

A **relation method in the linked-order model (LOM)** consists of the following.
- An étale Q-algebra A = K₁ × ⋯ × K_V (number fields of degrees d_i, D = Σd_i) and an order R ⊆ O_A generated by
  the candidate module.
- Candidates: elements θ(a) = Σ_{j=1}^r a_j w_j ∈ R, with integer coefficient vectors a examined by the
  algorithm. "Rank r" means r generators. NFS has r = 2, w = (1, −α) in each component; t-dimensional sieving
  has r = t.
- Smoothness integers: n_i(a) = (normalised) N_{K_i}(θ(a)_i). Their product W(a) is a decomposable form of
  degree D in r variables.
- Extraction, of one of two kinds:
  - **(ring-map)** ring maps ψ_i: R_i → Z/N on overorders R ⊆ R_i ⊆ O_A, agreeing on R. NFS, MNFS and TNFS use
    ψ_i = φ_i ∘ pr_i.
  - **(ramified)** ambiguous ideals / genus theory at primes ramified in some K_i. This covers class-group
    methods, including Seysen, Lenstra–Pomerance and Kirchner 2018.
- Cost: at least the number of examined candidates, plus B^ω for linear algebra (ω = 2), with B the factor-base
  size.

**Heuristic H.** Each examined candidate is B-smooth with probability ρ(log|W(a)|/log B)·L^{o(1)}. This is
heurnfs's H in per-candidate form. It excludes output-sensitive enumeration of smooth values (heurnfs lever (a)),
which in any case changes only constants (1.6386).

## 2. The generalized barrier: statements and proofs

### 2.1 Theorem E (elimination certificate) [proved]

Let S be a unital ring with maps ψ_i: S → R_i and ring maps φ_i: R_i → A, where Z/N ↪ A, and suppose the
composites φ_i ∘ ψ_i agree. Put J = Σ_i ker ψ_i (two-sided). Then J ∩ Z ⊆ NZ.

*Proof.* The common composite kills J, so d ∈ J ∩ Z gives d·1_A = 0, hence N | d. ∎

(codex resbarrier_1 Thm 1, independently.)

- J ∩ Z ≠ 0 iff (S/J) ⊗ Q = 0: the rings share no characteristic-0 point.
- One variable: J ∋ Res(f, g), so N | Res. Multivariate source: for every g ∈ J₂ with g ∉ J₁ ⊗ Q,
  |N_{K₁}(g(β))| ≥ N(J₂R₁) ≥ N. The ideal J₂R₁ lies in ker φ₁, which has index N.

**Theorem E alone gives no norm bound.** Three examples show why.
- f = x² − Nx − 1 and g = x² + Nx − 1 share the root 1 mod N, yet every x^j has norm ±1 on both sides (codex).
- Unipotent matrices in M_n(Z) have det 1.
- Z[ε]/ε² gives norm 1.

So every norm barrier needs a **sampling hypothesis**: candidates of bounded height in fixed lattice coordinates,
or equivalently computable cheaply without discrete logs (§3, units).

### 2.2 Lemma ND (non-étaleness is necessary for ring-map extraction) [proved, in the stated extraction form]

**Extraction form.** Relations give an element S ∈ R that is a square in A. The algorithm obtains two square roots
γ, γ′ ∈ A of the *same algebraic element* S and outputs gcd(ψ₁(γ) − ψ₂(γ′), N), where ψ₁, ψ₂ are ring maps to
Z/N from overorders R ⊆ R_i ⊆ O_A that agree on R. NFS, MNFS and TNFS have this form.

This restriction matters. Comparing roots of elements that are equal only mod N is a different scheme: in R = Z
with N = 15, the roots 1 and 4 of 1 ≡ 16 already give a factor. Dixon/QS-type schemes of that kind are outside the
lemma. They are covered by C2 through their own forms; QS uses x² − N, whose disc is 4N.

Let c = [O_A : R]. Then c·O_A ⊆ R, so R[1/c] = O_A[1/c].
- (i) **Agreement.** If gcd(c, N) = 1, every map from any intermediate overorder extends uniquely to R[1/c]
  by ψ(z) = ψ(cz)/ψ(c). So ψ₁ = ψ₂ on O_A. This covers non-maximal overorders and overorders that are not
  localisations (codex resbarrier_2).
- (ii) **Square roots.** γ² = S ∈ R makes γ integral, so γ ∈ O_A. Two roots of S differ by ε = Σ ±e_i, the
  e_i being the primitive idempotents. Then ψ₂(γ′) = ψ₁(ε)ψ₁(γ), and the values ψ₁(e_i) ∈ Z/N are computable
  without relations. So every split the relations can give is already given by gcd(ψ₁(e_i), N).
- (iii) **Partial overlap.** If 1 < gcd(c, N) < N, then gcd(c, N) is a factor. Here c is computable, e.g. as a
  multiple of |Res| times the order indices.

**Conclusion.** Nontrivial extraction needs rad(N) | c. For extraction through two *different* components at
p (the NFS diagram), the stronger condition (i′) holds: two maps through K₁ and K₂ that agree on R mod p must send
the idempotent e₁ to both 1 and 0, unless the prime of R below them lies on both components. For R = Z[θ] with
θ = (α₁, α₂) that means p | Res(f₁, f₂).

In general c = |Res(f₁, f₂)|·[O_{K₁} : Z[α₁]]·[O_{K₂} : Z[α₂]] (monic case). The draft's "c = |Res|" was wrong;
codex's counterexample is f₁ = x² − 5, f₂ = x − 1, with |Res| = 4 and c = 8. It is the Res part that the
two-component diagram forces. Hence disc(F₁F₂) = disc F₁·disc F₂·Res² is divisible by N² (k = 2).

For a general rank-2 module with generators u, v (critic A), argue directly. Both F_i vanish mod N on the line
au + bv ∈ ker φ, because the norm of an element of ker φ_i is divisible by N. So N | Res(F₁, F₂) whenever (u, v)
is primitive mod N and the F_i are coprime.

Lemma ND is proved only in this extraction form. Its ramified clause below is a reduction, not a theorem about
class-group methods.

**Ramified extraction (k = 1) [reduction, proved under its own hypothesis].** Ambiguous-ideal / genus-theory
extraction reveals only primes that ramify in some K_i, so those primes divide disc O_K, which divides the
discriminant of the integral norm form of any order in K. If they are only some of the primes of N, then
gcd(disc F, N) is a proper factor, computed for free. Otherwise N | disc F. This is a statement about
ramification-based extraction, not about every method that uses class groups.

Test: `test_lemma_nd_etale_gives_trivial_idempotent`.

### 2.3 Lemma C2 (rank-2 count) [proved; a corollary of known results] and Corollary RB2 [proved (model)]

**Lemma C2.** Let F ∈ Z[a, b] be a binary form of degree D ≥ 3 with disc F ≠ 0. For E ≥ 1 and X > 0:

  #{(a, b) ∈ Z², |a| ≤ E, 0 ≤ b ≤ E : 0 < |F(a, b)| ≤ X} ≤ A_F·X^{2/D} + 2(D+1)(2E+1),

where A_F = area{(x, y) ∈ R² : |F| ≤ 1} ≤ 3B(1/3,1/3)·|disc F|^{−1/(D(D−1))}, with 3B(1/3,1/3) = 15.8997 < 16.

**Sources.**
- The area bound is **Bean, Compositio Math. 92 (1994) 115–131, Theorem 1** (= Theorem 2 of the announcement,
  Bull. AMS 31 (1994) 204–207). It holds for complex forms of degree ≥ 3 with D_F ≠ 0, including reducible
  squarefree ones. Equality holds for cubics that are, up to a complex scalar, GL₂(R)-equivalent to xy(x − y).
- The lattice-point step is **Thunder 2001 (Annals 153) Lemma 14** for n = 2 (error 2D(2E+1)), i.e. Davenport's
  lemma (J. London Math. Soc. 26 (1951)).
- So C2 is a corollary of known results. Only its use as an NFS barrier is ours ("no prior statement found by
  us").

*Proof of the count (Davenport/Thunder argument).*
- Rows: count the closed sublevel set {|F| ≤ X} and discard zeros at the end. For fixed b, F(·, b) is either
  ≡ 0 (possible: F = ab(a − b) at b = 0; those points are discarded) or nonconstant. In the second case
  {a : |F(a, b)| ≤ X} ∩ [−E, E] is a union of ≤ D + 1 closed intervals. Each holds ≤ length + 1 integers.
- Columns: by Fubini, Σ_b length_b = ∫_{|a|≤E} #{b ∈ {0..E} : |F(a, b)| ≤ X} da. The b-sets are again ≤ D + 1
  intervals, so this is ≤ area + (D+1)(2E+1).
- Area: area{|F| ≤ X} = X^{2/D}A_F by homogeneity. ∎

**Consequences.**
- (a) **Two-map NFS:** |disc(F₁F₂)| ≥ Res² ≥ N² (k = 2). Checked in `disc_product_floor`.
- (b) **Ramified:** N | disc F (k = 1), with the §2.2 reduction.
- (c) **Sublattices (special-q, prescribed p-adic valuation), corrected.** Take a q-lattice M (det M = q) and
  candidates with cofactor |F(Ma′)|/q ≤ Y. They are exactly the points of that lattice with |F| ≤ qY. Their
  number is ≤ A_F(qY)^{2/D}/q + O(D·E_q), since the area scales by 1/det M.
  - Equivalently, G = F∘M/q has A_G = q^{2/D−1}A_F and disc G = disc F·q^{(D−1)(D−2)}. The draft's "main term/q"
    is right only at the original threshold.
  - Sum over the ≤ D lattices of each q ≤ Q (overlaps only over-count, which is harmless for an upper bound):
    #{special-q candidates with cofactor ≤ Y} ≤ A_F·Y^{2/D}·D·Σ_{q≤Q} q^{2/D−1} ≤ A_F·(QY)^{2/D}·D².
  - So the special-q gain is ≤ D²Q^{2/D}. With Q ≤ B = L[1/3, β] and D ≍ (x/log x)^{1/3} this is
    exp(O(log x)) = L^{o(1)}.
- (d) **Nonlinear binary parametrisations.** Suppose the family is homogeneous of degree e in (a, b), linked
  coefficientwise, so that mod each p | N the forms share a degree-e factor, with nonzero resultant. Then the
  Sylvester nullity is ≥ e, so N^e | Res (N squarefree), and the degree becomes De. This is worse. It is not a
  claim about arbitrary parametrisations.
- (e) **Selection inside the box, and the yield step.** C2 bounds *all* box points, so the sizes of any examined
  subset are controlled. The draft omitted the step from count to relation yield; here it is.
  - Let count(Y) := #{box points with 0 < |F| ≤ Y} ≤ 32·N^{−k/(D(D−1))}Y^{2/D} whenever this exceeds the error
    term.
  - Under H, the expected yield of any examined set of size e^C is
    ≤ Σ_{dyadic Y} min(n_Y, count(Y))·ρ(u_Y) ≤ (log₂ Y_max)·max_Y min(e^C, count(Y))·ρ(u_Y),
    with u_Y = log Y/log B.
  - A yield ≥ B therefore needs some Y with λ := log min(e^C, count(Y)) satisfying λ − u_Y log u_Y ≥ b − log x
    and λ ≤ C, where log Y ≥ k x/(2(D−1)) + (D/2)(λ − log 32).
  - This is the same ledger with λ ≤ C in place of C, so the optimisation below is unchanged.
  - **Hypothesis H5 (corrected after critic A, M1):** box side E ≤ B/(16(D+1)) for each box or q-lattice.
    - The error term 2(D+1)(2E+1) does not depend on X. If it could exceed B, those points might in principle
      supply the relations. Under H5 it is ≤ B/4, so a yield ≥ B needs Σ_Y min(n_Y, main(Y))·ρ(u_Y) ≥ B/2.
      The argument above then runs with main(Y) in place of count(Y).
    - In L-units H5 is e ≤ β. The GNFS and SNFS optima have e = β, so the floors are unchanged (as limits
      e = β − ε).
    - Algorithms whose boxes or sublattices have side > B^{1+o(1)} are not covered. That includes sparse
      selection from very large boxes. (The draft's condition E ≤ e^C/(8(D+1)) only showed that most
      candidates cannot all be small; it did not bound the smooth yield.) The uniform counting theorems carry Thunder's 𝓗(F)^{D−2}
    (𝓗 = ∏‖L_i‖, a Mahler-type height; Thunder 2001 Thm 3, Stewart 2024 Cor. 4).

**Corollary RB2 [proved (model): H + H5 (box/lattice side ≤ B^{1−o(1)}) + rank 2 + D ≥ 3 + linear algebra B²,
single N, RAM, no precomputation].**
- The ledger log X ≥ k x/(2(D−1)) + (D/2)·C is the GNFS program with κ = k/2: heurnfs writes the norm-size
  coefficient as κ/δ + δe with C = 2e and D = δ(x/log x)^{1/3}. So cost ≥ L[1/3, (16k/9)^{1/3}]
  (`ledger.rank2_constant`).
- Exponent: α ≥ 1/3 (`alpha_floor(1/2, 1/2)`). This holds for every rank-2 method of §2.2 form (two-map, or
  ramified with the reduction) and every N, special or not, with no counting heuristic for the polynomials.
- **Phrasing.** "SNFS is optimal among single-pair NFS" means: within this model, (32/9)^{1/3} is a floor for
  every single pair and every N, and SNFS-form N attain it. It does not mean the floor is attainable for every N. Batch NFS and the factorization
  factory (precomputation; area–time) are outside it.

| k | architecture | floor (ω = 2) | attained by |
|---|---|---|---|
| 1 | ramified / ambiguous-ideal (Kirchner) | (16/9)^{1/3} = 1.2114 | Kirchner's "very special n" (unconstructed; he also assumes ω = 2) |
| 2 | two-map NFS pair, all N | (32/9)^{1/3} = 1.5263 | SNFS |
| 2 | *symmetric* MNFS with V ≥ 2 equal-degree algebraic sides, RES only (kr → 1; ka ≥ 1/2, typical-size model) | ((36+16√6)/15)^{1/3} = 1.71143 | — (KNOWN constant, see below) |
| G holds (GAE, almost all N; full-box sampling + codex Thm 2) | NFS pair | 1.9230 | GNFS |
| G holds for disc (GAE-ram, almost all N; same caveat) | ramified | (32/9)^{1/3} = 1.5263 | Kirchner's "for most n" hope (unconstructed) |

**The MNFS row.**
- The exact optimum c³ = min_t 32(1/2+t)/(9t²(2−t)²) = (36+16√6)/15 at t = β_a/β = √(2/3) was found by codex
  resbarrier_2 (`ledger.mnfs_res_only_closed`).
- It **coincides with the known medium-prime MexTNFS constant 1.7114** (Kim–Barbulescu 2016; Sarkar–Singh 2016;
  Barbulescu–Pierrot 2014 / Pierrot 2015 for MNFS). It is the same Coppersmith-type program with the algebraic
  height term halved. Not new.
- Pairwise RES constrains only all but one side. Critic A ran mixed configurations (one side at t, the rest at
  1 − t): t = 0 gives 1.52629, t = 0.25 gives 1.64414, t = 0.5 gives 1.71143. So the RES-only floor over this
  family is the single-pair 1.5263. That is numerically supported, not proved.
- Under GAE the height assumptions of Coppersmith's MNFS program hold for almost all N. Its optimum is 1.9019.

### 2.4 Theorem GAE (heuristic G holds for almost all N) [proved]

Fix ε > 0 and A ≥ 1, and let D_max = (log N₀)^{1−δ} for some δ > 0. Let P be the set of pairs of coprime
integer polynomials (degrees d₁, d₂ ≤ D_max, heights H₁, H₂) with H₁^{d₂}H₂^{d₁} ≤ N₀^{A} and
∏(2H_i+1)^{d_i+1} ≤ N₀^{1−ε}. Then #{N ∈ [N₀, 2N₀) : some pair in P has a common root mod N} ≤ N₀^{1−ε+o(1)}.
(That pairs with H₁^{d₂}H₂^{d₁} > N₀^{A} are irrelevant at L[1/3] is a heuristic remark, not part of the
theorem.)

*Proof.* A common root gives N | Res ≠ 0. Since log|Res| ≤ D log(D+1) + A log N₀ = O(log N₀), each pair covers
≤ τ(Res) ≤ exp(O(log N₀/log log N₀)) = N₀^{o(1)} values of N. Summing over degree pairs and dyadic heights adds a
factor (log N₀)^{O(1)}. ∎

- The same argument with N | disc f proves **GAE-ram**: for almost all N, every squarefree f (disc f ≠ 0) with
  N | disc f has (d+1)h ≥ 1 − ε. Here |disc f| ≤ (2d)^{2d}H^{2d−2} ≤ N₀^{O(1)} under the same restrictions.
  (Without disc f ≠ 0, f = x² would defeat it.)
- GAE is a *necessary-height* theorem, not an existence theorem at the threshold.
- We found no explicit statement (Fried–Gaudry–Heninger–Thomé 2017 call SNFS inputs a "relatively rare set").
  It formalises the folklore counting heuristic that base-m heights N^{1/(d+1)} are about optimal. It is an
  elementary union bound, and we claim no novelty.
- **GAE bounds coefficient heights, not discriminants (critic A, M2).** The small-value count is governed by
  A_F, which is SL₂(R)-invariant, so heights cannot bound it for selected sub-families. The GAE-based floors
  (1.9230 for a single pair; the heights behind Coppersmith's 1.9019; 1.5263 for Kirchner's "for most n") are
  therefore **proved (model) only under full-box sampling, using codex resbarrier_1 Thm 2**. That theorem says
  most box points have log|FG| ≥ log M(f) + log M(g) + D log A − O(D); it has not been independently reviewed
  here. H says nothing about which sub-families have small norms, so it cannot fill this gap. A
  discriminant-version of GAE would need counts of forms by discriminant, which is hard for D ≥ 4.

### 2.5 Lemma KL (rank-r Liouville with a rational side) [proved]

Let W₁(a) = Σ_{j≤r} a_j v_j with v ∈ Z^r primitive (the rational side). Let W₂ be an integral form of degree d₂
in r variables, nonzero on Z^r \ {0} ∩ ker W₁ (for a norm form this needs the w_j to be Q-independent, so
r ≤ d₂), with W₂ ≡ 0 (mod N) on ker W₁ ∩ Z^r. This linking holds when
W₂ = c·N_{K₂}(Σ a_j w_j) with φ₂(w_j) = v_j mod N and gcd(c·denominators, N) = 1. Then

  N ≤ C(d₂+r−1, r−1)·H(W₂)·(√(r−1)·‖v‖₂^{1/(r−1)})^{d₂}.

*Proof.* ker W₁ ∩ Z^r is a rank-(r−1) lattice of covolume ‖v‖₂. By Minkowski it has a nonzero a with
|a|_∞ ≤ √(r−1)‖v‖^{1/(r−1)}. Then W₂(a) is a nonzero multiple of N, and |W₂(a)| ≤ #monomials·H(W₂)·|a|_∞^{d₂}. ∎

- In heights: h₂ + (d₂/(r−1))h₁ ≥ 1 − o(1), so min(h₁ + h₂) = min{1, (r−1)/d₂} (`ledger.kl_min`); the
  combinatorial factors are O(d log r)/x.
- Hypotheses that are indispensable: W₂ ≠ 0 on the nonzero points of the integer kernel; defined reduction maps;
  denominators coprime to N. The resultant form of the r = 2 case is standard in polynomial selection (Coxon,
  arXiv 1412.6011, Lemma 3.1; Bai's thesis §2.1.2).
- r = 2 is the resultant bound.
- With typical size H₁H₂E^{d₂+1}, the ledger is (r−1)x/d + (d+1)C/r. Its minimum over d is 2√((r−1)xC/r),
  increasing in r. So **with a rational side, higher-dimensional sieving is worse than r = 2 [proved (model):
  KL is rigorous, the typical-size ledger is the model part]**, and α ≥ 1/3.
- Test: `test_lemma_kl_rank3_power_basis` (exact, on rank-3 power-basis families of GNFS pairs).

### 2.6 Rank r ≥ 3 with two nonlinear sides [open; heuristic closure only]

- **Rigorous: 2-plane resultants.** Choose a rank-2 sublattice P of the candidate lattice, spanned by short
  vectors, on which the restricted forms W_i|_P have no common factor in characteristic 0 and a nonzero
  resultant. Such a P need not be a coordinate plane, since restrictions can share factors. A generic plane
  with small entries avoids the codimension-2 locus {W₁ = W₂ = 0} at a height cost r^{O(d)} = N^{o(1)}. We also
  need a coordinate j with u_j a unit mod N: a nonzero non-unit u_j splits N, and u ≡ 0 is trivial. On P the forms have a
  common root mod N, so Lemma RES gives d₂h₁(P) + d₁h₂(P) ≥ 1. If P is spanned by vectors of size O(1) in the
  family's coordinates, then h_i(P) ≤ h_i + O(d log r)/x, and the family height is ≥ x/d_max − O(d log r),
  independent of r.
- **What this permits.**
  - Minkowski gives no constraint below full rank. Rank-r sublattices of Λ = O₁ ×_{Z/N} O₂ on which the linking
    functional vanishes identically have no index constraint; NFS's (1, α) lattice is one.
  - So the displayed rigorous inequalities do not exclude the aggregate ledger log X ≈ x/d + (D/r)C. This says
    nothing about existence, constructibility, or other arithmetic constraints.
- **The 1/4 balance is conditional [heuristic, relabelled after critic A, M3].** It holds *if* log X ≳ x/d + D. That would follow from a random-walk size
  heuristic with balanced embedding profiles: the relevant scale is (Σ_j|σ(w_j)|²)^{1/2}, which can be much
  smaller than √r·max_j|σ(w_j)| for unbalanced profiles. Under that assumption log X ≳ √x and α ≥ 1/4
  (`alpha_floor(1/2, 0)`, = heurnfs `alpha_floor_aggregate`). It also hides the comparison between H(W) and
  the typical |W(a)|: for r-variable forms the coefficient-to-Mahler factor can be e^{O(rD)}. Without these
  heuristics we have no positive floor for this sub-case.
- **Heuristic G_r.**
  - Schmidt 1993 (Monatsh. Math. 115): c·X^{e(max(e,n)+1)} ≤ #{points of degree e in P^n, H ≤ X} ≤ C·X^{e(e+n)},
    for X large. With n = r − 1 and M = X^d this gives M^{max(d+1, r)} ≪ # ≪ M^{d+r−1}. These are
    fixed-dimension estimates. Applying them at growing d, r and to norm-form heights is an *assumption*.
  - The second assumption is that the r − 1 linking conditions (a common linear factor ℓ ∈ P^{r−1}(Z/N)) are
    equidistributed and independent mod N.
  - Under both assumptions, at total height h, the expected number of linked pairs is N^{e(h)}, with
    max(d+1, r)·h − (r−1) ≤ e(h) ≤ (d+r−1)·h − (r−1). (A lower count cannot bound the number of families from
    above, so the *upper* exponent is the one that matters for nonexistence.)
  - At the aggregate height h ≈ 1/d both ends are negative for 3 ≤ r ≤ d. So aggregate families are
    heuristically absent for typical N.
  - The LP `ledger.gr_min` (lower count) gives height floor (r−1)/(d+1). With the upper count the floor is
    (r−1)/(d+r−1). Both grow with r, so r = 2 is optimal and α ≥ 1/3, *conditionally*.
  - Power-basis t-dimensional sieving (heurnfs B2: 2(t−1)x/D) matches the lower-count floor, so it is
    G_r-optimal at leading order.
- **Why r = 2 is special (structural remark).**
  - Every binary form is decomposable over Q̄. So for r = 2, "norm form with a prescribed root mod N" is a
    linear condition: a lattice, solvable by reduction.
  - For r ≥ 3 decomposable forms are a thin (Chow-variety) subset, and lattice reduction returns
    non-decomposable forms, which carry no multiplicative structure.
  - The constructions known to us for r ≥ 3 are ring-generated modules (powers), which pay geometric heights.
    This is an observation, not a lower bound on constructions.

### 2.7 The hypotheses, and which are escape routes

| # | hypothesis | used by | escape? |
|---|---|---|---|
| H1 | extraction via ring maps or ramification in orders of an étale algebra (commutative, abelian) | ND | non-commutative: §3 (d) |
| H2 | candidates of bounded height in fixed lattice coordinates (sampling) | C2, KL | units / compact representations: need DLP mod N or class groups (§3 (i)) |
| H3 | smoothness integers are archimedean norms | C2 | product formula; p-adic = sublattice (§3 (e), (f)) |
| H4 | rank r = 2 | C2 | rank ≥ 3: KL (rational side), otherwise only heuristic G_r (§2.6) |
| H5 | box/lattice side ≤ B/(16(D+1)) (e ≤ β) | C2 (e) | sparse selection from boxes larger than B: needs an X-dependent error term (uniform Thue counting with Thunder's 𝓗(F)^{D−2}); not covered |
| H6 | heuristic H (smoothness of examined candidates) | RB2 | output-sensitive enumeration: constants only (heurnfs) |
| H7 | cost ≥ #examined + B² | RB2 | ω < 2: constants only (heurnfs B4) |
| H8 | for almost-all statements: N generic | GAE | special N: SNFS 1.526, Kirchner 1.211 (rank 2 floors) |

## 3. Escape table

"Reduces to" means: in the stated form, the escape lands in a case already covered by §2. It is not a lower
bound for every algorithm that uses that object.

| Escape | What it actually faces | α < 1/3? | Status |
|---|---|---|---|
| (a) more rings (MNFS) | Each used pair is a rank-2 system (C2, k = 2). V sides multiply the relation count by ≤ V ≤ T. | No. Symmetric MNFS ≥ 1.71143 under RES only (= known MexTNFS constant; mixed configurations reach 1.5263); 1.9019 is Coppersmith's optimum under the GAE-valid heights | proved (model) |
| (b) multivariate source Z[x₁..x_n]/I; varieties with a point mod N | Theorem E gives N(J₂R₁) ≥ N. The candidate module is a rank-r linked module (§2.6). | With a rational side: no (KL, model). Without one: no under G_r; rigorously open; floor 1/4 only under a typical-size heuristic | proved (model) / heuristic |
| (c) function fields over number fields; arithmetic surfaces | For specialisation-based designs (candidate functions restricted to horizontal curves Z_i, i.e. number fields), this is a rank-r linked module (H⁰(L) restricted to Z_i). Hilbert irreducibility makes horizontal divisors irreducible for generic candidates, so smoothness is carried by the vertical (norm) part. | Not for this design class: reduces to (b). Arbitrary adaptive arithmetic-surface algorithms are not covered. | heuristic sketch |
| (d) non-commutative: M_n, quaternion orders, Z[G] | Observations, not an exclusion. (1) For a full matrix representation, det φ(u) ≡ nrd u (mod N) is an integer identity, so the reduced-norm character carries no linking. (2) Indefinite Eichler orders of level N: the norm-one units surject onto the Borel diagonal (strong approximation; Voight, *Quaternion Algebras* Thm 28.2.11), so the diagonal character is not defined on classes. (3) Definite ones: mass φ(D_B)ψ(N)/12 ≥ N/12 (Eichler mass formula). That is not a runtime lower bound, but a factor base would have to resolve ≍ N classes. (4) Z[G] ⊗ Q splits by Wedderburn into simple algebras over number fields. Integral gluing and non-liftable modular representations are not analysed; the latter are group-order methods (route 3). | Not by any mechanism we could construct | open; withdrawn as a proved exclusion (codex resbarrier_2 #8) |
| (e) other multiplicative heights | Smoothness concerns integers. The product formula equates the finite and archimedean parts of the norm. A prescribed S-part (p-adic smallness) is a sublattice: C2(c), gain ≤ D²Q^{2/D} = L^{o(1)}. | No, in the sublattice form | proved (model) |
| (f) p-adic / adelic relations | Hom(Z_p, Z/N) = ∅ for p ∤ N (N is a unit versus 0). Globalising a local identity restores log|U| ≥ sL (codex (10)). Arbitrary non-global adelic identities are not analysed. | No, in the global form | proved (scoped) |
| (g) Kirchner 2018 | Ramified, k = 1: rank-2 floor 1.2114, equal to his special-n constant. GAE-ram: 1.5263 for almost all N, equal to his "for most n" hope. | No (exponent 1/3 in all rank-2 versions) | §4.1 |
| (h) Schnorr 2021 | D = 2 rational–rational linking, u − vN = w. For nonzero u, v, w, |uw| ≥ N − 1 (attained by u = N − 1, w = −1, v = 1). Charging smoothness of both u and w gives an L[1/2]-type balance at best (KNOWN: Dixon/CFRAC regime). His step from a short lattice vector to |u − vN| = p_n^{O(1)} fails: shortness gives only log W ≤ L + log v + log R₀ − L/(n+1) (codex (14)), not the needed |δ| ≲ B/(Nv) (codex (12)–(13)). Ducas (SchnorrGate, n = 47, 1000 trials): 353 / 65 / 0 / 0 relations at 10 / 20 / 40 / 400 bits. Aboumrad–Widdows–Kaushik (arXiv 2308.07804, Table 1): unique relations per lattice fall from 8.26% at 15 bits to 0.00% at 60 bits. The quantum-assisted variant (Yan et al., arXiv 2212.12372) was refuted by Khattar–Yosri (arXiv 2307.09651). | No | KNOWN + proved (|uw| ≥ N − 1) |
| (i) units, class groups, Stickelberger | Evaluating known unit powers mod N needs no DLP, but those relations carry only r₁ + r₂ square classes (Dirichlet). Matching a *prescribed* image needs a DLP mod N, which is ≥ factoring (Bach 1984). Class groups of Q(√−N): norms ≥ √N, giving L[1/2] (KNOWN: Schnorr–Lenstra 1984, Seysen 1987, Lenstra–Pomerance 1992). Gauss sums: a ring map Z[ζ_m] → Z/N needs primitive m-th roots of unity mod every p \| N (m \| p − 1). Without one, mapping to Z[ζ_m]/N makes the resultant constraint N^{φ(m)} (heurnfs B3). That is a height cost, not the quotient's cardinality as runtime. | No | proved / KNOWN |

## 4. Most promising escapes, developed

### 4.1 Kirchner's fields: exactly what is needed

The literature subagent read the full text (mathnet.ru/eng/pdma386), and literature critic B verified the quotes.
- **Method.** Class-group relations in a Galois extension L/K of prime degree p ("unlikely to be useful for any
  p > 2"). The goal is to factor Δ_{L/K}, at a cost measured in N(Δ_{L/K}). A nontrivial strongly ambiguous 𝔞
  gives gcd(Δ, 𝔞).
- **Height.** L = Q[X]/f with height(f) ≈ N(Δ)^{1/2/(d−1)}, the reason given in a parenthetical: "(the
  discriminant is a homogeneous polynomial of degree 2d − 2 in the coefficients of f)".
- **Special case.** L(1/3, (16/9)^{1/3}) "for very special n (so special that we do not know how to find them, or
  if they even exist)". This assumes n^{2+o(1)} linear algebra.
- **"Average-case factorization?".** "by relaxing the height of f to ≈ n^{1/d}, we may hope to factor the integer
  n in time which is only L_n(1/3, (32/9)^{1/3}) for most n. Of course, we also need the density of such f with a
  suitable subfield to be not too tiny, and a fast way to generate such a f."

**Our reading.**
- (1) **[proved (model)]** The special height N^{1/(2(d−1))} is the Vandermonde–Hadamard Liouville bound
  aL ≤ (2d−2)log M(f) + d log d (codex (11)). Lemma C2 with k = 1 then gives the floor (16/9)^{1/3}. His special
  case sits exactly on our rank-2 barrier. In the model (H5: box side ≤ B^{1−o(1)}), no rank-2 ramified method
  goes below it.
- (2) **[proved, necessary condition]** GAE-ram: for almost all N, any squarefree f with N | disc f has height
  ≥ N^{(1−ε)/(d+1)}. So his "for most n" relaxation to ≈ n^{1/d} is necessary up to d versus d + 1, giving 1.5263.
  GAE does not show that such f exist at the threshold, and certainly not within the thin family with "a suitable
  subfield".
- (3) **[assessment]** Even granting existence, a construction would at minimum have to find a small solution of
  disc(f) ≡ 0 (mod N) (a single congruence of degree 2d − 2 in d + 1 unknowns, total size ≈ N), plus the
  ramification, Galois-subfield and principalisation structure. We know no method below exhaustive search.
  - Linearising by prescribing a double root (f ≡ (x − m)²h mod N) imposes two linear conditions, so the generic
    lattice height is N^{2/(d+1)}, κ = 2, back to 1.923.
  - That is not a proof that every construction pays this cost.
- (4) **Exponent.** All versions are rank 2 with N | disc, so Corollary RB2 gives α ≥ 1/3 within the model. A
  rank-r ramified variant would face §2.6 again.

### 4.2 The aggregate rank-r family (the one place the barrier is only heuristic)

**Mechanism, if it existed (conditional).**
- Two degree-d fields K₁, K₂ with r ≪ d linked elements (φ₁(w_{1j}) = φ₂(w_{2j})), balanced embedding profiles,
  and family height H₁H₂ ≈ N^{1/d}.
- Candidates: Σ a_j(w_{1j}, w_{2j}) with a_j ∈ {0, ±1}, so C = r log 3 and log W ≈ x/d + d log r.
- Take d ≍ √(x/log r) ≍ √(x/log x) and r ≍ x^{1/4}(log x)^{3/4} (critic A #7, codex #7). Then log W ≍ √(x log x). With B = L[1/4] smoothing,
  u ≍ x^{1/4}, and the formal balance is **L[1/4]** (corrected parameters from codex resbarrier_2 #7; the draft's
  "r ≈ d" was inconsistent).
- Requirements: existence, construction, extraction maps (§2.2 form), relation rank.

**Status.**
- Under G_r (both Schmidt bounds), such families are absent for typical N: the expected number is N^{e} with
  e ≤ 1 − (r−1)(1 − 1/d) < 0 for r ≥ 3 (`gr_expected_log_count_upper`).
- No construction is known. Decomposability is nonlinear (§2.6), and the only ring-generated modules known to us
  (powers) have geometric heights.
- "Search costs ≥ N^{r−2}" in the draft was not a lower bound and is withdrawn.

**What would revive it.** A family of N, or a structure such as Galois or CM symmetry, that makes the r − 1
linking conditions *dependent* without the geometric height growth of powers. We found none.
- Cyclotomic fields have flat small profiles, but need primitive roots of unity mod p.
- "N is a norm" (O_K/(β) ≅ Z/N) gives a full-rank linked lattice with generic minima, which is useless.
- Fields of small discriminant usually have ring maps to Z/N, but computing them is root-finding mod N.

This is **open as a theorem**, excluded only under G_r. It is not a candidate algorithm.

## 5. Experiment E1 (pre-registered, exact)

Command: `.venv/bin/python -m src.resbarrier.e1` (≈ 2 s). Raw output: `experiments/raw/resbarrier_E1.jsonl` (36
instances).

- **P1 (theorem check): PASS.** 0 violations of the Lemma C2 bound with the N² floor, and 0 with the exact
  |disc(F₁F₂)| floor. All 36 pairs have N | Res and satisfy Hadamard.
- **P2 (saturation): PASS.** At n = 72, d = 3, median s_max(SNFS)/s_max(GNFS) = **59.4**. The pre-registered
  threshold was 8, and the predicted ratio ≈ 2^6 = 64.

| n | d | s_max GNFS (N² floor) | s_max SNFS | GNFS with exact-disc floor |
|---|---|---|---|---|
| 40 | 2 / 3 / 4 | 0.0081 / 0.017 / 0.028 | 0.268 / 0.093 / 0.157 | 0.19 / 0.21 / 0.18 |
| 56 | 2 / 3 / 4 | 0.0029 / 0.0081 / 0.014 | 0.268 / 0.127 / 0.157 | 0.21 / 0.21 / 0.17 |
| 72 | 2 / 3 / 4 | 0.00071 / 0.0028 / 0.0090 | 0.268 / 0.167 / 0.157 | 0.20 / 0.21 / 0.17 |

**Disclosed deviations and post-hoc check (after critic A).**
- The X grid starts at ⌊log₂ min |F|⌋ rather than the pre-registered ⌈·⌉. This is harmless: it adds at most
  one grid point.
- 6 of the 9 pre-registered SNFS instances have the reducible f = x^d − 1 (dk = e). That includes the P2
  instance (e = 72, d = 3). The pre-registration allowed this, since it only required disc ≠ 0.
- *Post-hoc supplementary run* (not pre-registered): N = 2^e + 3, f = x^d + 3·2^{dk−e}, irreducible by flint
  (`supplementary_irreducible_snfs`; `experiments/raw/resbarrier_E1_supp.jsonl`). Result: 0 violations in all
  9 cells. s_max = 0.071–0.133, which is 0.133 at n = 72, d = 3, i.e. 47× the GNFS median. The saturation
  conclusion does not rest on the degenerate pairs.

Reading.
- With the true discriminant the bound is saturated to a constant ≈ 0.2 for every pair. Lemma C2 is sharp up to
  a shape constant.
- GNFS sits below the N² floor by its discriminant excess N^{(2d−2)/(d+1)}, which is the factor-2 (κ = 2 versus 1)
  gap that heuristic G and Theorem GAE explain.
- SNFS saturates the N² floor.
- This is a correctness and sharpness check. It says nothing about speed.

## 6. Relation to earlier project notes

- heurnfs B1 used heuristic G. Here the exponent part (α ≥ 1/3) is proved without G (C2), and G itself is proved
  for almost all N (GAE).
- heurnfs B2 (t-dimensional sieving) is G_r-optimal (§2.6). With a rational side it is rigorous in the size
  model (KL).
- heurnfs B3 (towers) is consistent: mapping into Z[ζ]/N costs N^{φ(m)}.
- snfsrig Lemma RES is the one-variable case of Theorem E, combined with ND. escapeR2 Lemma TI concerns
  per-field parity, which is orthogonal: it constrains descent, not norm sizes.

## 7. GPT-6.1-sol blind round (codex/resbarrier_1.out, before our derivations)

**Agreements.** It independently found:
- Theorem E, in the same generality (non-commutative S);
- the determinantal/Fitting version;
- a full-box typical-norm theorem (its Thm 2: log|FG| ≥ (1 − o(1))√(2LC) for most box points), which
  complements our count (C2 bounds *all* box points, so it covers selection);
- the counterexamples showing that a sampling model is needed (unit powers, unipotents, Z[ε]/ε²);
- the Vandermonde bound (11) for N | disc;
- the precise Schnorr gap (12)–(14);
- the verdict that no listed escape gives α < 1/3, and that a theorem for all architectures is out of reach.

**Where we go further.**
- Lemma ND (why N | disc is necessary).
- Using the (known) Bean + Davenport/Thunder count as an NFS barrier with selection inside boxes of side ≤ B
  (C2 (e)), and its sublattice invariance.
- The k-indexed floors, including SNFS optimality in the model. The symmetric MNFS-RES floor turned out to be the
  known MexTNFS constant 1.7114.
- GAE.
- KL for rank r with a rational side.
- The Schmidt-count G_r.
- The 1/4 floor.
- The exact reading of Kirchner, whose full paper codex could not retrieve.

**Disagreement.** Codex ranked "ramification/class-group methods incl. Kirchner" as the top near-miss, because
"the resultant theorem does not exclude it". Lemma C2 with k = 1 does exclude exponent < 1/3 for that
architecture at rank 2.

## 8. Reviews

**GPT-6.1-sol adversarial round (codex/resbarrier_2.out): "substantial revision required".** It confirmed
Theorem E, the C2 count, GAE with its restrictions, and KL. No sub-L[1/3] mechanism emerged. Fixes applied:
- ND restricted to extraction from two square roots of one algebraic element (its counterexample was Z, N = 15).
  c = |Res|·indices, not |Res| (counterexample f₁ = x² − 5, f₂ = x − 1). Overorders are handled by R[1/c] =
  O_A[1/c].
- C2(c) restated with thresholds. The special-q gain is ≤ D²Q^{2/D} = L^{o(1)}.
- C2(d) hypotheses (homogeneous, coefficientwise linking, nonzero resultant) made explicit.
- The count-to-yield step was added.
- The ramified clause is now a reduction under its own hypothesis.
- GAE-ram requires disc f ≠ 0. GAE is a necessary-height theorem, not existence.
- §2.6: plane choice qualified. Schmidt bounds are labelled as fixed-dimension inputs, the independence of the
  r − 1 conditions is labelled an assumption, and the exponent is given as a range.
- §4.1 "construction, not existence" was withdrawn (existence is not established either). §4.2 parameters were
  corrected, and the claim "search ≥ N^{r−2}" was withdrawn.
- §3: (c) and (d) downgraded; (h) |uw| ≥ N − 1; (i) units and cyclotomic wording corrected.
- MNFS: the exact optimum (36+16√6)/15 was found by codex.

**Opus critic A (proofs and heuristics).** It reproduced E1 byte for byte and passed all tests. Three MAJOR
issues, all fixed:
- M1: the X-independent error term could in principle supply the B relations. H5 is now box side ≤ B/(16(D+1))
  (e ≤ β). The floors survive as limits, and larger boxes are not covered.
- M2: GAE bounds heights, not discriminants. The GAE-based floors are relabelled proved (model) under full-box
  sampling plus codex Thm 2.
- M3: the 1/4 floor was mislabelled and is now heuristic.

Minor fixes:
- the KL hypothesis r ≤ d₂ (Q-independence);
- the direct N | Res(F₁, F₂) argument for general rank-2 modules;
- a generic plane in §2.6 and the unit-coordinate remark;
- d ≍ √(x/log r) in §4.2;
- MNFS: the 1.7114 floor holds for *symmetric* MNFS only, while mixed configurations reach 1.5263;
- E1 disclosures plus a post-hoc irreducible-SNFS run;
- a real test of the ND idempotent mechanism (`test_lemma_nd_idempotent_mechanism`).

**Opus critic B (literature).** All fixes applied:
- Bean is Compositio *Theorem 1* (= BAMS Thm 2), with the "up to a complex scalar" equality clause.
- The C2 count is Thunder 2001 Lemma 14 / Davenport 1951, so C2 is presented as a corollary of known results.
- The Kirchner quotes were corrected ("with a suitable subfield"; "for most n"; the parenthetical).
- 1.7114 is the known MexTNFS constant (Kim–Barbulescu 2016; Sarkar–Singh 2016), so the novelty claim was
  removed.
- Eichler mass φ(D_B)ψ(N)/12 and strong approximation are cited to Voight GTM 288.
- Bach 1984 (DLP ⇒ factoring); Schnorr–Lenstra 1984 / Seysen 1987 / Lenstra–Pomerance 1992 (class groups,
  L[1/2]); Coxon arXiv 1412.6011 / Bai (resultant bound in polynomial selection); Yan et al. 2022 and
  Khattar–Yosri 2023 (Schnorr variant); Aboumrad–Widdows–Kaushik Table 1.
- "Within LOM, single N, no precomputation (Batch NFS, factory excluded)" was added to the verdict.
- No prior statement was found of (a) a C2-type NFS lower bound, (b) GAE, or (c) ND. ND's closest sources are
  Stevenhagen 2008 (NFS as Z × Z[α] → diagonal of Z/n × Z/n) and Buchmann–Lenstra 1994. We claim no novelty.
- Further citations added: BGGM 2015; Gaudry–Grémy–Videau 2016 (higher-dimensional sieving); Kleinjung 2006;
  Lee–Venkatesan 2018; Bernstein–Lange 2014/921; Coppersmith 1993.

**Not done.** A second codex round on the revised note was not run. The revisions only narrow claims and add
scope labels.

## 9. Reproduction

```bash
.venv/bin/python -m src.resbarrier.e1                 # E1 (pre-registered), ~2 s
.venv/bin/python -m pytest -q tests/test_resbarrier.py
.venv/bin/python -c "from src.resbarrier.e1 import supplementary_irreducible_snfs as f; f()"   # post-hoc E1 supplement
.venv/bin/python -c "from src.resbarrier.ledger import mnfs_floors; print(mnfs_floors())"   # 1.7114 / 1.9019
```
