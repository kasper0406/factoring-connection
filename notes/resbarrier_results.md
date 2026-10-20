# resbarrier: the generalized resultant barrier and its escapes (route 2)

Agent `resbarrier`, phase 16, 2026-10-20.
- Code: `src/resbarrier/forms.py` (exact binary-form tools, Lemma RES/C2/KL checkers),
  `src/resbarrier/ledger.py` (exact LPs, constants, exponent floors), `src/resbarrier/e1.py` (experiment E1).
- Tests: `tests/test_resbarrier.py`.
- Pre-registration: `notes/resbarrier_prereg.md` (commit `81d9f71`, before E1 ran). Raw data:
  `experiments/raw/resbarrier_E1.jsonl`.
- GPT-6.1-sol: blind round `codex/resbarrier_1.{md,out}` (§7); adversarial round `codex/resbarrier_2.{md,out}` (§8).
- Literature: one research subagent (Kirchner 2018 read in full from mathnet.ru; Bean 1994; Thunder 2001;
  Stewart 2024; Schmidt 1993; Schnorr refutations). Reviews: two read-only Opus critics (§8).

Labels: **proved** (proof here, or a cited published theorem), **proved (model)** (an exact consequence of a stated
heuristic model), **heuristic** (named assumption), **KNOWN** (cited), **experimentally supported**, **open**.

---

## 0. Verdict

**Plain language.** No escape gives a factoring exponent below 1/3, not even heuristically. The resultant wall is
more general than "two polynomials": it is a statement about the *discriminant of the order that the candidates
live in*. Any relation method that extracts factors through ring maps or ramification needs that discriminant to
be divisible by N, and then a classical area theorem (Mahler–Bean) limits how many candidates can have small
norms. This holds for every polynomial pair and every way of choosing candidates from a box or lattice, for
special N too, with no counting heuristic for the polynomials. The only gap we cannot close rigorously is
higher-rank families (sieving over ≥ 3 coefficients) with no rational side. There a counting heuristic closes
it, and even if that heuristic failed maximally the floor would be 1/4, not lower.

**Precisely** (single N, RAM model, sparse linear algebra exponent ω = 2, heuristic H = "an examined candidate's
norm is B-smooth as often as a random integer of the same size"):

1. **Theorem E (elimination certificate) [proved; also found blind by codex].** For any rings, commutative or
   not, with maps to Z/N that agree on a common source, every nonzero integer in the sum of the kernels is a
   multiple of N. Res(f, g) is the one-variable instance. This is the only universal content of "Res ≥ N".
   Without a sampling model there is *no* norm lower bound: unit powers and unipotent matrices have norm ±1
   (codex counterexamples, §2.1).
2. **Lemma ND (non-étaleness is necessary) [proved].** Suppose an order R in an étale algebra A extracts factors
   through ring maps to Z/N. If gcd([O_A : R], N) = 1, every congruence of squares it yields is y ≡ ψ(ε)x with
   ψ(ε) computable in polynomial time without any relations. So relations are useless unless N | [O_A : R],
   which gives N² | disc(R). The ramified (class-group) architecture instead needs N | disc(O_K).
3. **Lemma C2 + Corollary RB2 (rank-2 discriminant barrier) [proved count; proved (model) cost].** Let the
   smoothed integers be values of a binary form F of degree D ≥ 3 with N^k | disc F. Then, in any box of side E,
   #{|F| ≤ X} ≤ 16·N^{−k/(D(D−1))}·X^{2/D} + 2(D+1)(2E+1).
   - Under H this forces α ≥ 1/3 and a constant ≥ (16k/9)^{1/3}.
   - k = 2 (NFS pairs, any polynomials): 1.5263, i.e. **SNFS is optimal among single-pair NFS for every N**.
   - k = 1 (ramified / class-group architectures): 1.2114, exactly Kirchner's special-N constant.
   - It covers special-q and p-adic cofactors (sublattices) and nonlinear binary parametrisations.
   - Experiment E1: 0 violations, and the SNFS/GNFS saturation ratio was 59 against a predicted ≈ 64.
4. **Theorem GAE (G for almost all N) [proved, elementary].** For all but N₀^{1−ε+o(1)} integers N ∈ [N₀, 2N₀],
   no pair of polynomials of total count height ∏(2H_i+1)^{d_i+1} < N₀^{1−ε} has a common root mod N. The same holds
   for single polynomials with N | disc f. So heuristic G of heurnfs is a theorem for almost all N.
   - Combined with heurnfs B1, GNFS 1.9230 is the H-model floor of a single two-map pair for almost all N.
     GAE also makes rigorous the height assumptions behind Coppersmith's MNFS program, whose optimum is
     1.9019.
   - Kirchner's general-N hope, 1.5263, is the almost-all floor of his architecture. It requires constructing f
     with N | disc f at the counting threshold: a nonlinear congruence of degree 2d − 2, which no lattice method
     reaches (§4.1).
5. **Lemma KL (rank-r Liouville with a rational side) [proved].** Let the rational side be Σ a_j v_j and the
   other side a degree-d norm form W₂ of rank r. Then
   N ≤ C(d+r−1, r−1)·H(W₂)·(√(r−1)‖v‖^{1/(r−1)})^d.
   The height floor is (r−1)x/d, so higher-dimensional sieving with a rational side is rigorously worse than
   r = 2 (α ≥ 1/3 in the size model).
6. **Rank r ≥ 3 with no rational side [open rigorously; heuristic closure].**
   - Only the 2-plane resultants are rigorous. They give family height ≥ x/d_max *independent of r*, so they
     permit the "aggregate ledger" log X ≈ x/d + (D/r)C. Under a random-walk size heuristic that ledger has
     floor **α ≥ 1/4**.
   - Heuristic G_r closes the gap: count rank-r norm forms by Schmidt's bounds for algebraic points of fixed
     degree, and treat the r − 1 linking conditions as independent. Then the height floor grows like
     (r−1)x/(d+1), r = 2 is optimal, power-basis sieving is G_r-optimal, and α ≥ 1/3.
   - G_r also predicts that aggregate families do not exist for typical N: the expected number is
     N^{(d+1)/d − (r−1)} ≪ 1.

**Escape verdicts (§3 table).**
- MNFS: constants only.
- Multivariate presentations and 2-dimensional schemes reduce to rank-r linked modules: item 6.
- Non-commutative rings: the abelian quotient is the reduced norm, which gives trivial relations; otherwise
  strong approximation kills the linking (indefinite), or the class number is ≍ N (definite Eichler orders).
- Other heights: the product formula; p-adic smallness is a sublattice, covered by C2.
- p-adic/adelic rings: no ring maps Z_p → Z/N for p ∤ N.
- Kirchner: §4.1.
- Schnorr: D = 2 rational–rational linking. The norm product is ≥ N, so the best possible is L[1/2]. The
  claimed |u − vN| = p_n^{O(1)} does not follow from lattice shortness (§3).
- Units, class groups and Stickelberger: norm-free but information-free, unless one solves a discrete log mod N,
  or computes a class group of discriminant ≥ N (L[1/2]), or uses a ring map Z[ζ_m] → Z/N (needs m | p − 1).

**Best mechanism.** No mechanism survives. The best near-miss is the aggregate rank-r family of item 6, which
would give L[1/4] if it existed and could be constructed (§4.2). G_r says it does not exist for typical N, and no
construction is known even for special N. This is not a candidate algorithm. It marks the one place where the
barrier is a heuristic rather than a theorem.

---

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

### 2.2 Lemma ND (non-étaleness is necessary for ring-map extraction) [proved]

Let R ⊆ O_A and c = [O_A : R]. Let the ring maps ψ₁, ψ₂ (from overorders) agree on R.

- (i) If p ∤ c, then R_(p) = O_{A,(p)}, so the two maps agree mod p. But a map through pr_i kills the idempotent
  e_j for j ≠ i. So two maps through different components both reduce mod p only if p | c.
- (ii) Two square roots of the same S ∈ A^{×2} differ by ε = Σ ±e_i. When gcd(c, N) = 1, the e_i lie in
  R[1/c], so ψ(ε) ∈ Z/N is computable directly. Then y = ψ(ε)x, and gcd(ψ(ε) ∓ 1, N) gives every split that
  any relation could give.
- (iii) If 1 < gcd(c, N) < N, then gcd(c, N) is already a factor.

So nontrivial ring-map extraction needs N | c (N squarefree; in general every p | N divides c). Since
disc R = disc O_A · c², this gives **N² | disc R**. For R = Z[x]/(f₁f₂), c = |Res(f₁, f₂)|, which recovers Lemma RES.

- A single field with a single map yields only ±γ, which is trivial.
- Ramified extraction therefore needs p | disc O_K for the primes it reveals, so N | disc R (k = 1).
- Test: `test_lemma_nd_etale_gives_trivial_idempotent`.

### 2.3 Lemma C2 (rank-2 count) [proved] and Corollary RB2 [proved (model)]

**Lemma C2.** Let F ∈ Z[a, b] be a binary form of degree D ≥ 3 with disc F ≠ 0. For E ≥ 1 and X > 0:

  #{(a, b) ∈ Z², |a| ≤ E, 0 ≤ b ≤ E : 0 < |F(a, b)| ≤ X} ≤ A_F·X^{2/D} + 2(D+1)(2E+1),

where A_F = area{|F| ≤ 1} ≤ 3B(1/3,1/3)·|disc F|^{−1/(D(D−1))}, with 3B(1/3,1/3) = 15.8997 < 16.

The area bound is **Bean, Compositio Math. 92 (1994) 115–131, Theorem 2**, verified by the literature subagent.
Equality holds for cubics GL₂(R)-equivalent to xy(x − y).

*Proof of the count.*
- Rows: for fixed b, {a ∈ R : |F(a, b)| ≤ X} is a union of ≤ D intervals. Its endpoints are among the ≤ 2D
  roots of F(·, b) = ±X; F(·, b) is nonconstant because disc F ≠ 0. An interval of length ℓ holds ≤ ℓ + 1
  integers. So the count is ≤ Σ_{b=0}^{E} len_b + D(E+1).
- Columns: Σ_b len_b = ∫_{|a|≤E} #{b ∈ {0..E} : |F(a, b)| ≤ X} da. For each a the b-set is again ≤ D intervals,
  so this is ≤ area({|F| ≤ X} ∩ box) + D·2E.
- Area: area{|F| ≤ X} = X^{2/D}A_F by homogeneity. ∎

**Consequences used.**
- (a) **Two-map NFS:** disc(F₁F₂) = disc F₁·disc F₂·Res(F₁,F₂)², and the discriminants are nonzero integers, so
  |disc| ≥ N² (k = 2). Checked in `disc_product_floor`.
- (b) **Ramified:** N | disc (k = 1).
- (c) **Sublattices** (special-q, prescribed p-adic valuation, S-unit cofactors). For M ∈ M₂(Z) with det M = q
  and G = F∘M/s integral, disc G = disc F · q^{D(D−1)}/s^{2(D−1)}. For s = q this is disc F·q^{(D−1)(D−2)}, and
  the main term becomes (main term of F)/q. Summing over all q-lattices re-partitions the plane at a loss of
  O(D log log Q). So special-q and p-adic smallness gain nothing beyond L^{o(1)}.
- (d) **Nonlinear binary parametrisations** a ↦ θ(a, b) of degree e. The forms share a degree-e factor mod p, so
  N^e | Res and N^{2e} | disc(F₁F₂). The degree is De, which is worse.
- (e) **Arbitrary selection inside the box.** The bound counts *all* box points, so any subset is covered. To
  find e^C candidates with |F| ≤ X it suffices that E ≤ e^C/(8(D+1)): then the error term is ≤ e^C/2, so
  e^C ≤ 2·16·N^{−k/(D(D−1))}X^{2/D}, i.e.
  **log X ≥ k·x/(2(D−1)) + (D/2)(C − log 32).**
  - Boxes much larger than the budget are not covered by the elementary lemma. The uniform counting theorems
    carry a factor H(F)^{D−2} (Thunder 2001 Thm 3; Stewart 2024 Cor. 4), so we do not claim them.

**Corollary RB2 [proved (model)].**
- The ledger above is the GNFS program with κ = k/2, since the norm-size coefficient is κ/δ with
  D = δ(x/log x)^{1/3}. Under H with LA B², cost ≥ L[1/3, (16k/9)^{1/3}] (`ledger.rank2_constant`).
- **Exponent.** Smoothing L[α] needs C ≳ (log X/log B)·log(…). With log X ≳ √(k x C), the generic ledger, this
  gives α ≥ 1/3 (`alpha_floor(1/2, 1/2)`). This holds for every two-map or ramified rank-2 method, every N
  (special or not), and needs no counting heuristic for the polynomials.

| k | architecture | floor (ω = 2) | attained by |
|---|---|---|---|
| 1 | ramified / class group (Kirchner) | (16/9)^{1/3} = 1.2114 | Kirchner's "very special n" (unconstructed) |
| 2 | two-map NFS pair, all N | (32/9)^{1/3} = 1.5263 | SNFS |
| 2 | MNFS with V ≥ 2 algebraic sides (Res only: kr → 1, and ka ≥ 1/2 from pairwise algebraic resultants; typical-size model) | 1.7114 (numeric) | — |
| G holds (GAE, almost all N) | NFS pair / MNFS | 1.9230 / 1.9019 | GNFS / MNFS |
| G holds for disc (GAE-ram, almost all N) | ramified | (32/9)^{1/3} = 1.5263 | Kirchner's general-N hope (unconstructed) |

The MNFS row (no prior statement found by us): under Lemma RES alone, MNFS with V ≥ 2 sides cannot go below
1.7114; the overall RES-only floor stays the single-pair 1.5263. The value comes from heurnfs's program
`mnfs(kr=1, ka=0.5)`.

### 2.4 Theorem GAE (heuristic G holds for almost all N) [proved]

Fix ε > 0 and A ≥ 1, and let D_max = (log N₀)^{1−δ} for some δ > 0. Let P be the set of pairs of coprime
integer polynomials (degrees d₁, d₂ ≤ D_max, heights H₁, H₂) with H₁^{d₂}H₂^{d₁} ≤ N₀^{A} and
∏(2H_i+1)^{d_i+1} ≤ N₀^{1−ε}. Then #{N ∈ [N₀, 2N₀) : some pair in P has a common root mod N} ≤ N₀^{1−ε+o(1)}.
(Pairs with H₁^{d₂}H₂^{d₁} > N₀^{A} have typical norms ≥ N₀^{A/D}·… and are irrelevant at L[1/3].)

*Proof.* A common root gives N | Res ≠ 0. Since log|Res| ≤ D log(D+1) + A log N₀ = O(log N₀), each pair covers
≤ τ(Res) ≤ exp(O(log N₀/log log N₀)) = N₀^{o(1)} values of N. Summing over degree pairs and dyadic heights adds a
factor (log N₀)^{O(1)}. ∎

- The same argument with N | disc f (|disc f| ≤ (2d)^{2d}H^{2d−2} ≤ N₀^{O(1)} under the same restrictions) proves
  **GAE-ram**: for almost all N, every
  f with N | disc f has (d+1)h ≥ 1 − ε.
- The literature subagent found no statement of this anywhere (Fried–Gaudry–Heninger–Thomé 2017 call SNFS inputs
  a "relatively rare set"). We claim no novelty; it is elementary.
- GAE only bounds polynomial heights. To turn heights into norm sizes we use the full-box typical-size theorem
  (codex resbarrier_1 Thm 2: most box points have log|FG| ≥ log M(f) + log M(g) + D log A − O(D)), or H for
  sub-families.

### 2.5 Lemma KL (rank-r Liouville with a rational side) [proved]

Let W₁(a) = Σ_{j≤r} a_j v_j with v ∈ Z^r primitive (the rational side). Let W₂ be an integral form of degree d₂
in r variables, nonzero on Z^r \ {0} ∩ ker W₁, with W₂ ≡ 0 (mod N) on ker W₁ ∩ Z^r. This linking holds when
W₂ = c·N_{K₂}(Σ a_j w_j) with φ₂(w_j) = v_j mod N and gcd(c·denominators, N) = 1. Then

  N ≤ C(d₂+r−1, r−1)·H(W₂)·(√(r−1)·‖v‖₂^{1/(r−1)})^{d₂}.

*Proof.* ker W₁ ∩ Z^r is a rank-(r−1) lattice of covolume ‖v‖₂. By Minkowski it has a nonzero a with
|a|_∞ ≤ √(r−1)‖v‖^{1/(r−1)}. Then W₂(a) is a nonzero multiple of N, and |W₂(a)| ≤ #monomials·H(W₂)·|a|_∞^{d₂}. ∎

- In heights: h₂ + (d₂/(r−1))h₁ ≥ 1 − o(1), so min(h₁ + h₂) = (r−1)/d₂ (`ledger.kl_min`).
- r = 2 is the resultant bound.
- With typical size H₁H₂E^{d₂+1}, the ledger is (r−1)x/d + (d+1)C/r. Its minimum over d is 2√((r−1)xC/r),
  increasing in r. So **with a rational side, higher-dimensional sieving is rigorously (in the size model) worse
  than r = 2**, and α ≥ 1/3.
- Test: `test_lemma_kl_rank3_power_basis` (exact, on rank-3 power-basis families of GNFS pairs).

### 2.6 Rank r ≥ 3 with two nonlinear sides [open rigorously]

- **Rigorous (2-plane resultants).** Take any 2-plane P spanned by basis vectors. Restricted to P the forms are
  binary with a common root mod N, so d₂h₁ + d₁h₂ ≥ 1 (Lemma RES on W_i|_P, whose heights are ≤ H(W_i)). The
  family height is therefore ≥ x/d_max, *with no dependence on r*.
- **Minkowski does not help below full rank.** The linked lattice Λ = O₁ ×_{Z/N} O₂ has covolume ≍ N. For r < D,
  rank-r sublattices on which the linking functional vanishes identically have no index constraint. NFS's (1, α)
  lattice is one.
- **So the rigorous constraints permit the aggregate ledger log X ≈ x/d + (D/r)C**, with r up to min d_i.
  - Under the random-walk size heuristic, |Σ a_jσ(w_j)| ≈ √r·max_j|σ(w_j)| for most a, so there is a
    D·(½)log r term. With that, log X ≳ x/d + D, so log X ≳ √x and **α ≥ 1/4** (`alpha_floor(1/2, 0)`; this
    matches heurnfs `alpha_floor_aggregate`).
- **Heuristic G_r.** The number of rank-r norm forms of degree d with Mahler height ≤ M lies between
  M^{max(d+1, r)} and M^{d+r−1}. This is Schmidt 1993 for points of degree d in P^{r−1}, converted from absolute
  height X to M = X^d. Linking needs r − 1 conditions mod N: the common linear factor ℓ ∈ P^{r−1}(Z/N).
  - The lower count gives Σ max(r, d_i+1)h_i ≥ r − 1, i.e. h ≥ (r−1)/(d+1) (`ledger.gr_min`). The upper count
    gives (r−1)/(d+r−1).
  - In both cases the x-coefficient grows with r, so r = 2 is optimal and α ≥ 1/3.
  - Power-basis t-dimensional sieving (heurnfs B2: 2(t−1)x/D) has the same x-coefficient, so it is
    G_r-optimal.
  - The expected number of aggregate families (h_total ≈ 1/d) is N^{(d+1)/d − (r−1)} ≪ 1 for r ≥ 3
    (`gr_expected_log_count`).
- **Why r = 2 is special (structural remark).** Every binary form is decomposable over Q̄, so for r = 2 "norm
  form with a prescribed root mod N" is a *linear* condition: a lattice, solvable by reduction. For r ≥ 3,
  decomposable forms are a thin (Chow-variety) subset. Lattice reduction returns non-decomposable forms, which
  carry no multiplicative structure. Rank r ≥ 3 constructions are therefore limited to ring-generated modules
  (powers), which pay geometric heights.

### 2.7 The hypotheses, and which are escape routes

| # | hypothesis | used by | escape? |
|---|---|---|---|
| H1 | extraction via ring maps or ramification in orders of an étale algebra (commutative, abelian) | ND | non-commutative: §3 (d) |
| H2 | candidates of bounded height in fixed lattice coordinates (sampling) | C2, KL | units / compact representations: need DLP mod N or class groups (§3 (i)) |
| H3 | smoothness integers are archimedean norms | C2 | product formula; p-adic = sublattice (§3 (e), (f)) |
| H4 | rank r = 2 | C2 | rank ≥ 3: KL (rational side), otherwise only heuristic G_r (§2.6) |
| H5 | box of side ≤ e^C/(8(D+1)) | C2 (e) | sparse selection from huge boxes: needs uniform Thue counting; heuristic H covers it |
| H6 | heuristic H (smoothness of examined candidates) | RB2 | output-sensitive enumeration: constants only (heurnfs) |
| H7 | cost ≥ #examined + B² | RB2 | ω < 2: constants only (heurnfs B4) |
| H8 | for almost-all statements: N generic | GAE | special N: SNFS 1.526, Kirchner 1.211 (rank 2 floors) |

## 3. Escape table

| Escape | Norm lower bound it actually faces | α < 1/3? | Status |
|---|---|---|---|
| (a) more rings (MNFS) | each used pair is a rank-2 system (C2, k = 2); V sides multiply relations by ≤ V ≤ T | no; constant ≥ 1.7114 (RES only), 1.9019 for almost all N | proved (model) |
| (b) multivariate source Z[x₁..x_n]/I, varieties with a point mod N | the candidate module is a rank-r linked module (§2.6); Theorem E gives N(J₂R₁) ≥ N | rational side: no (KL). No rational side: no under G_r; floor 1/4 rigorously-in-size even if G_r fails | proved (KL) / heuristic (G_r) |
| (c) function fields over number fields, arithmetic surfaces | horizontal divisors over Q are irreducible for almost all candidates (Hilbert irreducibility), so smoothness lives on the vertical part, i.e. norms of restrictions to horizontal curves Z_i = number fields; H⁰(L) restricted to Z_i is a rank-r linked module | no; reduces to (b) | proved reduction (structural); heuristic for the count |
| (d) non-commutative: M_n, quaternion orders, Z[G] | det/nrd is the abelianisation: det φ(u) ≡ nrd u (mod N) is an integer identity, so the relations are trivial. Eichler orders of level N carry a Borel character, but indefinite ones have unit groups surjecting onto (Z/N)^× (strong approximation; the character is not defined on classes) and definite ones have class number ≍ N/12 (Eichler mass formula). Z[G] reduces by Wedderburn to number fields or simple algebras. Non-liftable modular representations are group-order methods (route 3) | no | proved (structural), standard facts |
| (e) other multiplicative heights | smoothness concerns the norm; by the product formula the finite part equals the archimedean part. Allowing an S-part (p-adic smallness) is a sublattice: C2 (c) | no | proved |
| (f) p-adic / adelic relations | Hom(Z_p, Z/N) = ∅ for p ∤ N (N invertible versus 0); globalising local identities restores (10) of codex, log|U| ≥ sL | no | proved |
| (g) Kirchner 2018 | ramified architecture, k = 1: rigorous floor 1.2114 (matches his special-N constant); GAE-ram: 1.5263 for almost all N (matches his general hope) | no (exponent 1/3 in all versions) | §4.1 |
| (h) Schnorr 2021 lattices | D = 2 rational–rational linking (u ≡ w with u − vN = w). Linked pairs with |u|·|w| ≤ X number ≍ X log X/N, so the smoothed integers multiply to ≥ N: L[1/2] at best even with free smooth u, v. His step "short vector ⇒ |u − vN| = p_n^{O(1)}" fails: lattice shortness gives only log W ≤ L + log v + log R₀ − L/(n+1) (codex (14)), not the |δ| ≲ B/(Nv) precision required (codex (12)–(13)). Ducas: 0 relations in 1000 trials at 40 and 400 bits (n = 47); Aboumrad–Widdows–Kaushik: unique relations per lattice fall from 8.26% (15 bits) to ≈ 0% (60 bits) | no | KNOWN (refutations) + proved (our D = 2 count) |
| (i) units, class groups, Stickelberger | norm ±1 / smooth by construction, but no information. Unit-twisted candidates need a DLP in (Z/N)^×, which is ≥ factoring (Bach). Cl(Q(√−N)) relations have norms ≥ √N: L[1/2], KNOWN (Seysen, Lenstra–Pomerance). Gauss sums need Z[ζ_m] → Z/N, i.e. m-th roots of unity mod N (m | p − 1), else Z[ζ_m]/N costs N^{φ(m)} (heurnfs B3) | no | proved / KNOWN |

## 4. Most promising escapes, developed

### 4.1 Kirchner's fields: exactly what is needed and why nobody can build them

Read in full from mathnet.ru/eng/pdma386 by the literature subagent.
- The method is class-group relations in a Galois extension L/K of prime degree p, in practice p = 2.
- The number factored is N(Δ_{L/K}). A nontrivial strongly ambiguous ideal 𝔞 gives gcd(Δ, 𝔞).
- L = Q[X]/f, with height(f) ≈ N^{1/(2(d−1))}, "because the discriminant is a homogeneous polynomial of degree
  2d − 2 in the coefficients".
- This yields L(1/3, (16/9)^{1/3}) "for very special n (so special that we do not know how to find them, or if
  they even exist)", and possibly (32/9)^{1/3} for general n "by relaxing the height of f to ≈ n^{1/d}". That
  requires "the density of such f … not too tiny, and a fast way to generate such a f."

**Our reading [proved / model].**
- (1) Kirchner's special height N^{1/(2(d−1))} is exactly the Vandermonde–Hadamard Liouville bound
  aL ≤ (2d−2)log M(f) + d log d (codex (11)). Lemma C2 with k = 1 then gives the floor (16/9)^{1/3}. His special
  case sits exactly on our rigorous barrier, and no ramified method can go below it, for any N.
- (2) For general N, GAE-ram shows that f with N | disc f exist only at (d+1)h ≥ 1 for almost all N. This is
  height n^{1/(d+1)}, Kirchner's relaxed height up to d versus d + 1, giving 1.5263. Exponent 1/3 in every case.
- (3) **Construction is the obstacle, not existence.** At the counting threshold about one f per H^{d+1}/N
  exists, but finding one means solving disc(f) ≡ 0 (mod N): one congruence of degree 2d − 2 in d + 1 unknowns,
  with solutions of total size N^{1}.
  - Brute force costs ≈ N/τ.
  - Coppersmith-type methods reach only roots of total size ≈ N^{1/(2d−2)}.
  - Linearising by prescribing the double root (f ≡ (x−m)²h mod N) imposes 2 linear conditions, so the height
    is N^{2/(d+1)}, κ = 2, and the result is GNFS 1.923: no gain.
  - So "why nobody can construct them" is precise: a nonlinear small-solution problem mod N, of the same type as
    computing ring maps for fields not built from N.
- (4) **Different exponent?** No. Every version is rank 2 with N | disc, so Corollary RB2 gives α ≥ 1/3. A
  rank-r ramified variant would face §2.6 again.

### 4.2 The aggregate rank-r family (the only place the barrier is heuristic)

**Mechanism, if it existed.**
- Two degree-d fields K₁, K₂ with r ≈ d elements w_{ij} ∈ K_i, linked (φ₁(w_{1j}) = φ₂(w_{2j})), with family
  height H₁H₂ ≈ N^{1/d}.
- Candidates: Σ a_j(w_{1j}, w_{2j}) with a_j ∈ {0, ±1} (C = r log 3), so log W ≈ x/d + d log r.
- With d ≈ √x and r ≈ x^{α′}: log X ≈ √x log x. B = L[1/4] smoothing gives u ≈ x^{1/4}, and the cost
  L[1/4, O(1)] needs C ≈ x^{1/4}(log x)^{3/4}, i.e. r ≈ x^{1/4}. This would be an **L[1/4] factoring algorithm**.
- Requirements: (i) existence; (ii) construction; (iii) relation rank. (iii) is plausible by the usual heuristics.

**Kill (heuristic G_r).**
- Existence: the expected number of linked pairs at that height is N^{(d+1)/d − (r−1)}. That is ≪ 1 already at
  r = 3, under either Schmidt count.
- Construction: decomposability is nonlinear (§2.6), so no lattice method applies. Search costs ≥ N^{r−2}.

**What would revive it.** A family of N, or a structure such as Galois or CM symmetry, for which the r − 1
linking conditions are dependent: one condition implies the others. The ring structure provides this only for
powers of a single element, which costs geometric heights. We found no other mechanism.
- Cyclotomic fields have flat small profiles, but their ring maps to Z/N need m | p − 1.
- "N is a norm" (O_K/(β) ≅ Z/N) gives a full-rank linked lattice with generic minima N^{1/(2d)}, which is useless.
- Fields of small discriminant almost always have ring maps to Z/N (probability ≈ (1 − 1/e)^ω per field), but
  computing them is root-finding mod N.

Status: **open as a theorem**. It is excluded under G_r and needs a dependent-linking structure that nobody has
exhibited.

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
- The count with any selection (C2 (e)) and its sublattice invariance.
- The k-indexed floors, including SNFS optimality and the MNFS-RES 1.7114.
- GAE.
- KL for rank r with a rational side.
- The Schmidt-count G_r.
- The 1/4 floor.
- The exact reading of Kirchner, whose full paper codex could not retrieve.

**Disagreement.** Codex ranked "ramification/class-group methods incl. Kirchner" as the top near-miss, because
"the resultant theorem does not exclude it". Lemma C2 with k = 1 does exclude exponent < 1/3 for that
architecture at rank 2.

## 8. Reviews

(pending: codex adversarial round resbarrier_2; Opus critic A (proofs/heuristics); Opus critic B (literature).)

## 9. Reproduction

```bash
.venv/bin/python -m src.resbarrier.e1                 # E1 (pre-registered), ~2 s
.venv/bin/python -m pytest -q tests/test_resbarrier.py
.venv/bin/python -c "from src.resbarrier.ledger import mnfs_floors; print(mnfs_floors())"   # 1.7114 / 1.9019
```
