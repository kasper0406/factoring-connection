# nonexact: O1, non-exact and evaluated compressions modulo N

Agent `nonexact`, phase 17, 2026-10-20. DRAFT for review: codex adversarial round and Opus critics are pending
(§8).

**Files.**
- Code: `src/nonexact/`.
  - `core.py`: exact helpers (reductions mod (N, F), flint LLL, kernel lattices, exact Fincke–Pohst enumeration,
    the Prop LIN counts, the Prop IH relation extractor).
  - `run_experiment.py`: the pre-registered runs E1–E3.
  - `posthoc.py`: disclosed post-hoc diagnostics.
- Tests: `tests/test_nonexact.py` (15).
- Pre-registration: `notes/nonexact_prereg.md` (commit e4e203e, before any registered run).
- Config: `experiments/configs/nonexact.json`.
- Raw results:
  - pre-registered: `experiments/raw/nonexact_E{1,1b,2,3}.jsonl`;
  - post hoc: `nonexact_E1_distinct_posthoc.jsonl`, `nonexact_E1_exhaustive_posthoc.jsonl`,
    `nonexact_E2_baseline2_posthoc.jsonl`, `nonexact_E3_reasons_posthoc.jsonl`.
- Codex (GPT-6.1-sol): blind round `codex/nonexact_1.{md,out}`; adversarial round `codex/nonexact_2.{md,out}`.

**Labels.** PROVED (proof here) · PROVED (model) (exact consequence of a named model) · KNOWN (cited) ·
HEURISTIC (assumption named) · measured · open · assessment. Novelty is never claimed.

**Notation.** As in `notes/ptesplit_formulation.md` and ptebarrier §1.1.
- x = log N; B = L[β] with β < 1/3; M ≍ x/log B pieces; ℓ := log V = O(x^{2β}(log x)^{1−2β}) is the
  small-side target (T-v).
- P_A = ∏_{a∈A}(X − a) for a multiset A ⊂ [−B, B]. The presentation is F̃ ∈ Z[X] monic of degree d.
- I_N := (N, F̃) ⊂ Z[X].
  - View A: φ(X) = m ∈ Z with N | F̃(m), and J_m := (N, X − m) = ker(Z[X] → Z/N, X ↦ m).
  - View B: φ(X) = μ ∈ Z/N, and J_μ := (N, X − μ̃).
- W := 2V + 1 is the number of admissible small values.

---

## 0. Verdict

**Plain language.**
- We find **neither a barrier nor a construction**. We do prove that O1 is not a separate escape.
- Once the reduction is allowed to hold only modulo N, or only in value at m, the polynomial vocabulary drops out:
  - presentations, multipliers, small remainders and PTE fibres no longer constrain anything;
  - the problem *becomes* the general one: list many split-by-construction products congruent to something
    small modulo N, faster than random search.
- Two consequences follow.
  - **No counting barrier is possible.** Such compressions exist in abundance; we prove a coincidence count
    unconditionally. Any barrier would have to be a lower bound on *finding* them: a complexity statement about a
    natural class of factoring algorithms. The rigidity methods that closed the exact case do not supply that.
  - **Lattice reduction does not find them cheaply.**
    - It imposes the congruence mod N linearly, but loses the split condition, which is thin.
    - With the last s roots linearised, the cost per relation is (N/V)·exp(−(1+o(1))√(2 log(N/V)·log B)), which
      is N^{1−o(1)}: exponential. This is a model statement, confirmed at toy scale.
    - Multiplier lattices reach Minkowski heights N^{d/n} but return generic, non-split vectors.

**Proved (main items).**
- **Prop EQ (collapse).** Evaluated O1 is *exactly* the listing problem Π (§1.3): every evaluated member is an
  integer relation, and every relation is a member on F̃ = X − m. The non-exact ring version is contained in
  Π (through J_μ ⊇ I_N).
- **Prop BRK (where R and W break).**
  - (i) mod (N, F̃). Lemma EX forces every non-exact difference to have height ≥ N − 2H. Lemma R then needs
    M ≥ 2 log₂ N ≈ 2.885x, while members need only M ≍ x/log B. The archimedean step 4 of Lemma R is vacuous.
  - (ii) Evaluated. The defect ideal J_m contains X − m, of degree 1 < k. So Lemma EX's "height ≥ N" fails
    and value-smallness gives no Mahler bound (step 3).
- **Prop AB (abundance, unconditional).** Any S distinct split products contain ≥ S²/N − S ordered pairs that are
  congruent mod N. These are doubly-split relations, abundant as soon as S > 2N.
- **Prop CS.** Rigorous upper counts of members in a lattice coset are stuck at the trivial bound (one root
  determined). Beating it needs power cancellation in character sums over intervals of length p^{o(1)}, far
  below Burgess's range. The determinant method does not apply (Prop BP-NO).
- **Prop IH (inheritance is relation-valued, all k).**
  - This settles the gap ptebarrier left open ("defects are injective, but the conditions are not
    independent"), in the PTE setting with d ≥ k.
  - Any inherited compatibility of two non-exact fibre members yields ≥ 1 congruence
    ∏_{a∈A}(b − a) ≡ v (mod N) of the target type. For k ≥ 2 it yields ≥ ⌈(M − 2k + 3)/(2(k − 1))⌉ distinct
    ones, with pieces ≤ 2B and |v| ≤ 2H·Σ_{l<k}B^l.
  - Codex proved the existence part independently.
- **Prop NF′.** N-independent fibres inherit compatibility for at most an e^{−x+o(x)} fraction of balanced
  semiprimes.
- **Prop COH.** A family whose congruences follow linearly from e + 1 seed congruences is an instance of O2: split
  members of a fixed rank-(e+1) polynomial lattice, with an N-aligned basis. Rank-2 pencils with linear pieces
  collapse to rank 1: two seeds force μ to be a small rational mod N, and every member is then the same
  relation.
- **Prop LIN (model).** The linearised-lattice cost per relation is
  ≥ max(N/(T_s W), D_s/T_s) ≥ (N/W)·exp(−(1+o(1))√(2 ln(N/W)·ln(2B+s))).
- **Prop LAT / VA.** The lattice determinant (Smith form), Minkowski and Henk counts. In View A, multiplier
  compressions are points of the 2-dimensional value lattice of determinant N (the Dixon/CFRAC ledger).

**Experiments (pre-registered; failures disclosed).**
- **E1 (linearised lattice).**
  - 9 of 15 cells were run; 6 needed more than 10⁶ lattices and are reported as predictions.
  - Hit ratios 0.75–1.20 in 7 cells. **P1 failed in 2 cells** (0.54, and 4.94, which formally triggered the
    kill signal).
  - Post-hoc diagnosis: hits are repeat detections of a few full multisets, because prefixes are drawn with
    replacement from a small space. The 395 hits are 4 distinct members.
  - Exhaustive member counts are 31 against 28.0, 564 against 584.6, and 4 against 0.84 expected. The last has
    Poisson p ≈ 0.01; we found no structure in it.
  - The LLL + Fincke–Pohst enumerator matched the exhaustive box counts in 450/450 lattices.
  - The measured best s equals the model's arg-min in every (n, B) row.
  - The lattice beat brute force by at most the suffix count C(2B+s, s), as predicted.
- **E2 (multiplier lattices, 21 cells, 2100 lattices).**
  - ‖b₁‖/GH has medians 1.01–1.15. Residual exponents sit at d/n (within 0.075 below).
  - 0 useful vectors with deg Q ≥ 1 for d = 2. One for d = 1, a lucky Dixon-type vector.
  - **P8's root-rate comparison failed in 3 cells.** The LLL multipliers have *fewer* integer roots than the
    registered baseline. Post hoc: that baseline was miscalibrated; against a uniform-height baseline the
    rates agree within ~1.5×.
- **E3 (Prop IH, 246 constructed pairs).**
  - Clauses (a)–(d): 0 violations. k = 1 gives exactly one relation per pair.
  - **The registered distinct-count bound failed in 31 k = 2 pairs.** It ignored the sign identification
    D(b′) = −D(b). The corrected bound, with the factor 2 now in Prop IH, has 0 violations.

**Assessment.** Probability < 1% that O1 yields sub-L[1/3]. Codex (blind): about 5%. The difference is about
how likely an unknown output-sensitive lister for Π is. Both of us find that O1 *is* that lister problem.

**Residual open statement (§6).** For fixed β < 1/3, is there a uniform classical algorithm Π_β of the following
kind?
- It outputs, from N alone and in time L[β], at least L[β] congruences ∏_i c_i ≡ v (mod N).
- The pieces c_i are integers ≤ L[β], or values m − a at one computed m.
- 0 < |v| ≤ exp(O(x^{2β}(log x)^{1−2β})), or v is itself a product of pieces.
- The exponent rows have rank ≥ (1 − o(1))·#columns.

O1 is equivalent to Π_β (Prop EQ). Every mechanism examined here either reduces to it, reduces to O2
(coherent families), or pays an exponential search (Prop LIN; the query model).

---

## 1. Formalisation (Task 1)

### 1.1 The two versions

A member is (A, h, Q, r) as in ptebarrier §1.1. Its relation is φ(h)φ(P_A) = φ(r).

- **(E) Exact.** hP_A − r ∈ F̃Z[X]. Closed in regime (ptebarrier Theorem W, Corollary P, Prop S).
- **(O1-i) Non-exact.** hP_A − r ∈ I_N \ F̃Z[X].
  - By Lemma EX (F̃ monic), red_F̃(hP_A − r) = N·D_A with 0 ≠ D_A ∈ Z[X]_{<d}.
  - Equivalently hP_A − F̃Q_A − r = N·T_A for some Q_A, T_A ∈ Z[X] (codex eq. 5). Bounds on the certificate Q_A
    are *extra* restrictions; they are not part of the congruence.
  - The congruence is d scalar conditions mod N: d − k exact vanishings and k windows.
- **(O1-ii) Evaluated.** hP_A − r ∈ J_m (View A) with |r(m)| ≤ V, and no height condition on r. It is
  non-trivial when h(m)P_A(m) ≠ r(m).
  - J_m = {g : g(m) ≡ 0 (mod N)}, since g = (X − m)s + g(m) (codex). This is **one** scalar condition.

**Inclusion [PROVED].** N | F̃(m) gives I_N ⊆ J_m, strictly for d > 1, since the quotients have N^d and N
elements. So exact ⊂ (i) ⊂ (ii), as sets of relations at m. The View B analogue uses J_μ.

### 1.2 Lemma COL and Prop EQ [PROVED, elementary]

**Lemma COL.**
- (a) An evaluated member on any presentation (F̃, m) is exactly an integer congruence
  u ≡ v (mod N) with u = h(m)P_A(m), v = r(m), |v| ≤ V. The data F̃, Q and deg/height of r impose nothing:
  replacing r by the constant v preserves the condition (codex).
- (b) Conversely, any multiset A and integer v with P_A(m) ≡ v (mod N) is an evaluated member on F̃ = X − m with
  h = 1 and r = v. At m = 0 and F̃ = X, this is P_A(0) = ∏(−a) ≡ v (codex eq. 51).

*Proof.* J_m is the kernel of evaluation at m, followed by reduction mod N. ∎

**The listing problem Π(N; 𝔛, V).** The factor base is 𝔛 = {x_a}, a set of |𝔛| ≤ L[β] residues computable from N:
small integers, m − a, or μ − a. Task: output many distinct multisets A with ∏_{a∈A} x_a ≡ v (mod N), where
v ∈ Z, 0 < |v| ≤ V (the "one-sided" form), or v = ∏_{a′∈A′} x_{a′} with A′ ≠ A (the "doubly split" form).

**Prop EQ.**
- (i) The class of evaluated O1 families on some presentation equals the class of families of one-sided Π
  outputs with 𝔛 = {m − a} (Lemma COL).
- (ii) Non-exact families with an evaluation map are contained in Π with 𝔛 = {μ − a}: I_N ⊆ J_μ. The small side
  r(μ) must still be factored, either by r splitting over small roots (the doubly split form) or through a second
  order.
- (iii) Pieces of size ≤ B^{O(1)} can be re-encoded as roots at m = 0, changing B to a fixed power of itself
  (codex).

**Scope.**
- The presentation must be *computed*, and its cost charged. An N-dependent F̃ is not free advice (codex).
- The containment Π ⊆ O1-(i) needs d = 1. With a prescribed d ≥ 2 and bounded certificates, Π is not
  automatically contained in O1-(i) (codex).
- So "O1 ≡ Π" is exact for (ii). For (i) it is "O1-(i) ⊆ Π, with equality at d = 1".

**Consequence.** For (ii), the polynomial structure of a family can only matter through how its members are
*generated*, never through what they are. A barrier for O1 is therefore a barrier for Π.

**Useful rows (codex).** Multiplying both sides of u ≡ v by the same smooth t gives new pairs with the same
exponent-difference row. So conditions (a)–(c) of systeq §6 must be read with a rank requirement. Π_β (§0)
includes one.

### 1.3 Which step of Lemma R / Theorem W breaks (Prop BRK) [PROVED]

**(i) Non-exact.**
- *Lemma R input.* In a PTE fibre with d ≥ k, take two non-exact members whose compatibility is not exact (E ≠ 0
  below). Then G = P − P′ = (r − r′) + N·E with 0 ≠ E ∈ Z[X]_{<k}, so ‖G‖_∞ ≥ N − 2H.
- *Threshold.* Lemma R's threshold M ≥ k + 2⌊log₂(2^{k−1}√k·h·B^{k−1})⌋ + 3 then needs
  M ≥ 2 log₂ N − O(1) = 2x/ln 2 − O(1). Members need only M·log(2B) ≥ x − 1 (T-u), and the cost range has
  M ≍ x/log B = o(x). So the regime is never reached.
- *Step 4.* Lemma R's step 4 ("|U′(a)| ≤ Δ, hence at most log₂Δ distant roots") is vacuous as soon as
  log₂Δ ≥ deg U′. That already holds when M ≤ log₂ N.
- *Theorem W.* Eliminating F̃ from two certificates gives (codex eq. 10)
  Q′hP − Qh′P′ = Q′r − Qr′ + N(Q′T − QT′).
  - The added term has degree up to d + j₂ and height ≥ N after any cancellation.
  - The certificate quotients Q_A are forced; they are not small.
  - Both inputs of Lemma R fail.
- *Explicit failure of W's conclusion.* ptebarrier's N = 25 example: F = P + 26, with P and P′ both ≡ −1 mod
  (25, F). It has two non-exact members with different Ψ. (It is outside the L-range; the point is the
  mechanism.)
- *What survives.*
  - Corollary P's step "Q_A = Q_{A′}" survives, by degree.
  - Prop NX survives: in regime, equal defects force equal members.

**(ii) Evaluated.**
- *The defect ideal.* J_m ∩ Z[X]_{<k} contains X − m. So, unlike (N, F̃) with deg F̃ ≥ k, its nonzero elements
  need not have height ≥ N: Lemma EX fails exactly because deg(X − m) = 1 < k.
- *Values carry no height.* A small value gives no height: r = C(X − m) + v has r(m) = v for every C (codex
  eq. 12). For fibre pairs, G = P − P′ is constrained only through G(m) ≡ v − v′ (mod N). That gives no bound on
  M(G), so step 3 (Mahler measure) has no input.
- *k = 1 is the exception.* Then G is a constant and (ii) coincides with (i): NX1 applies.

**Quantitative summary.** With log H + log H_Q = O(ℓ), the exact case sits in Lemma R's regime, since ℓ/M → 0
(codex eq. 1). Replacing the height by N^{1−o(1)} adds 2x/ln 2 to the threshold, against M = o(x).

---

## 2. Barrier attempts (Task 2)

### 2.1 No counting barrier: abundance

**Prop AB (coincidences) [PROVED].**
- *Statement.* Let 𝒰 be a set of S distinct objects with residues ρ: 𝒰 → Z/N. Examples: distinct integers
  u ↦ u mod N (View A); multisets A ↦ P_A(μ) mod N (View B, formal products). Then
  #{(u, u′) ∈ 𝒰², u ≠ u′, ρ(u) = ρ(u′)} ≥ S²/N − S.
- *Proof.* With n_c = #ρ^{−1}(c), Cauchy–Schwarz gives Σ_c n_c² ≥ (Σ n_c)²/N = S²/N; subtract the diagonal. ∎
- *Consequence.* Distinct split products are plentiful.
  - Take products of exactly M primes from [B^{1/4}, B]: they are distinct and number N^{c−o(1)} when
    Mb = cx (codex eq. 31).
  - Once S ≥ 2N·L, such a catalogue contains ≥ S(2L − 1) congruent pairs u ≡ u′.
  - These are *doubly split* relations, smooth by identity on both sides (ptecons row 2b). In View B they are
    genuine multiplicative relations among the symbols μ − a.
- *Not claimed.* Prop AB counts coincidences, not independent rows.

**Conditional one-sided abundance [PROVED conditional; codex].** Assume square-root cancellation
|Σ_{p∈[B^{1/4},B]} χ(p)| ≤ B^{1/2}x^{C₀} for every nonprincipal χ mod N (hypothesis CS). Then for c > 2 the
catalogue above has (1 + o(1))S·W/N members with |v| ≤ V. Unconditional versions at modulus N and
y = N^{o(1)} are not available. For example, Soundararajan's smooth-numbers-in-progressions theorem needs modulus
≤ y^{4√e−ε} (codex).

**Measured (post hoc, E1).** Exhaustive member counts in whole M-multiset spaces (View B, random μ):

| (n, B, M) | space | members | expected S·W/N |
|---|---|---|---|
| (24, 8, 8) | 735 471 | 31 | 28.0 |
| (24, 16, 7) | 15 380 937 | 564 | 584.6 |
| (32, 8, 10) | 5 311 735 | 4 | 0.84 |

**Reading.**
- Every statement that bounds the members of a large enough family (S ≫ N) below ≈ S/N is false.
- So the target's obstruction is not existence. It is the cost of *finding* members when the family's density
  of members is ≈ W/N.
- The only cheap regime would be over-dense subfamilies ("conspiracies") built from N, which is the construction
  problem itself.

### 2.2 The lattice view (Prop LAT, Prop VA)

**Prop LAT [PROVED; determinant formula from codex eq. 14–16].**
- *Setting.* Fix P. Let L_P = {(Q, r) : deg Q ≤ j, deg r < k ≤ d, QP − r ∈ I_N}, of dimension n = j + 1 + k.
- *Determinant.* det L_P = ∏_i N/gcd(N, δ_i), where the δ_i are the Smith invariants of
  C_P = [P mod F̃, …, X^jP mod F̃ | −E_k]. Equivalently det L_P = p^{rank_p}·q^{rank_q}.
  - If j + 1 ≥ d and P is a unit in (Z/N)[X]/(F̃), then det L_P = N^d. `test_multiplier_lattice_index_and_membership`
    checks this, and so does E2 (2100/2100 lattices).
  - If n ≤ d and C_P has full column rank mod p, then the box with heights < p contains only 0. This is a
    genuine fixed-P barrier under a rank hypothesis (codex).
- *Minkowski.* A nonzero (Q, r) exists with ‖·‖_∞ ≤ N^{d/n} [KNOWN].
- *Henk 2002, Thm 1.5 [KNOWN].* #(L_P ∩ K) < 2^{n−1}∏_i⌊1 + 2/λ_i(K)⌋ for symmetric convex K.
  - With K = [−R, R]^n, R ≥ λ_n, and Minkowski II (∏λ_i ≥ det/n! for the cube), this is
    ≤ 2^{n−1}3^n n!·R^n/N^d.
  - So it equals the Gaussian count up to e^{O(n log n)}.
- *Barroero–Widmer (IMRN 2014, Thm 1.3) [KNOWN].* This gives the matching asymptotic with an explicit
  error in the minima (codex).
- *What it does not give.*
  - The index gives no lower bound on individual minima: Z^{n−d} × NZ^d has determinant N^d and many short
    vectors (codex).
  - Compressions with a multiplier exist for *every* P at height N^{d/n}. The multiplier Q is generic (E2), not
    split by construction.
  - If Q is required to be split, QP is a split polynomial of degree M + j, and we are back at j = 0.

**Prop VA (View A value ledger) [PROVED].**
- *Statement.* In View A, (Q, r) ↦ (Q(m), r(m)) maps L_P into Λ_u = {(q, v) : qu ≡ v (mod N)},
  u = P(m), a 2-dimensional lattice of determinant N.
  - The kernel contains every (Q, r) with Q(m) = r(m) = 0: trivial relations, e.g. multiples of X − m.
  - Λ_u is the whole value content of L_P: any polynomial "compression" with Q(m) ≠ 0 is a point of Λ_u.
  - The exact average is Σ_{u mod N}#{(q, v) ∈ Λ_u : 0 < |q| ≤ Q₀, |v| ≤ V} = 2Q₀(2V + 1). So on average the
    residual product |q||v| ≈ N: the Dixon/CFRAC ledger.
- *Consequence.* In View A, polynomial multipliers add nothing beyond the 2-dimensional value lattice.

### 2.3 Charging the split structure

**Prop CS (character sums) [PROVED; the reading is an assessment].**
- *Setting.* Let p | N, I = {μ − a : |a| ≤ B} ⊂ F_p^* with n = |I| < p, and
  R_M(c) = #{(x₁..x_M) ∈ I^M : ∏x_i ≡ c (mod p)}.
- *Statement.* |R_M(c) − n^M/(p−1)| ≤ n·S_max^{M−2}, with S_max = max_{χ≠χ₀}|Σ_{x∈I}χ(x)|.
- *Proof.* Use orthogonality, and (p−1)^{−1}Σ_χ|S_χ|² = n. ∎
- *Trivial case.* With S_max ≤ n this gives R_M(c) ≤ n^M/(p−1) + n^{M−1}: the "last root determined" bound. So
  members are a fraction ≤ W/n + W/(p−1) of all tuples.
- *Reading.*
  - To reach the heuristic n^M/p even up to a factor 2 one needs S_max ≤ n·(n/p)^{1/(M−2)}. For M ≍ x/log B,
    that is S_max ≲ n^{1/2+o(1)}: **square-root cancellation over intervals of length L[β] = p^{o(1)}**.
  - Larger M only weakens this to a power saving.
  - Burgess (1962) gives power saving only for length ≥ p^{1/4+ε} [KNOWN].
  - Assessment: no rigorous count of split members of a lattice coset beyond the trivial one is within reach.
    This is the "Davenport-type bound" for the joint count.

**Prop BP-NO (the determinant method does not apply) [PROVED: the reductions; assessment: inapplicability].**
- (a) *Coefficient space.* The coefficient vectors of integer-split polynomials are Zariski-dense in A^M.
  - The Jacobian of (a) ↦ (e) is the Vandermonde, so the map is dominant.
  - A polynomial vanishing on all images pulls back to one vanishing on all integer tuples (codex).
  - So no proper variety contains them, and Bombieri–Pila or dimension-growth bounds (BP 1989; CCDN 2020)
    reduce to the trivial dimension count.
- (b) *Root space.* The congruence ∏(μ̃ − a_j) ≡ v (mod N) is not an algebraic condition over Q.
  - Clearing N gives Y: Nt + v − ∏(μ̃ − a_j) = 0 in A^{M+2}.
  - Y is linear in (t, v) and ruled by the lines (a, v + Nλ, t − λ).
  - Its integral points in a box are in bijection with the congruence solutions. The determinant method, which
    exploits curvature of a fixed variety over Q, cannot see the congruence.
  - Heath-Brown's p-adic variant (Annals 2002) still counts points of a variety over Q.
- (c) *Counting mod p in a box.* Counting points of a hypersurface over F_p in a box of side 2B + 1 = p^{o(1)}
  needs exponential sums. Completion needs side > p^{1/2}; Burgess-type arguments need side > p^{1/4}. Both are
  out of range.
- (d) *Where the determinant method does apply: the exact part of structured families.*
  - Members lying on a fixed irreducible curve over Q of degree D number ≤ c·D²B^{1/D}(log B)^κ
    (Binyamini–Cluckers–Novikov), independently of N. This is ptebarrier's k = 1 argument.
  - So a structured family of size L[β] must live on a variety of dimension ≥ 1 with many integer points, and its
    N-alignment must then hold on a dense set. That is Prop COH below.

### 2.4 Inheritance

**Prop IH (inherited compatibility is relation-valued, all k) [PROVED].**
- *Hypotheses.* P ≠ P′ are monic of degree M, split over Z, with roots in [−B, B]. D := P − P′ has deg D < k.
  δ ∈ Z[X]_{<k} satisfies D ≡ δ coefficientwise (mod N) and D ≠ δ. Put E := (D − δ)/N ∈ Z[X]_{<k} \ {0}.
- *(a) Congruences.* For every b ∈ Z, D(b) ≡ δ(b) (mod N). At a root b of P′ that is not a root of P:
  u_b := P(b) = ∏_{a∈A}(b − a) ≠ 0, and u_b − δ(b) = N·E(b). Symmetrically, at roots of P,
  u_b := P′(b) ≡ −δ(b).
- *(b) Non-triviality.* u_b ≠ v_b ⇔ E(b) ≠ 0.
  - Values in supp P ∩ supp P′ are roots of D, hence of E (because |δ(b)| < N).
  - The signed root-multiplicity measure of the coprime cofactors has vanishing moments 0..M−k. By Vandermonde
    its support has ≥ M − k + 2 points.
  - So there are ≥ M − 2k + 3 "good" points b: non-common, with E(b) ≠ 0.
- *(c) Distinct relations, counted up to the sign identification (u, v) ~ (−u, −v).*
  - deg D = 0: exactly one (u_b = ±C). This is ptebarrier's NX1.
  - deg D ≥ 1: at least ⌈(M − 2k + 3)/(2·deg D)⌉, since |D(b)| = t has ≤ 2 deg D solutions.
- *(d) Sizes.*
  - Every piece satisfies |b − a| ≤ 2B and |v_b| ≤ ‖δ‖·Σ_{l<k}B^l.
  - If 2B < min p | N, then u_b is a unit mod N, so v_b ≠ 0 (codex).
- *Application.* Let members P ≡ r and P′ ≡ r′ mod (N, F̃) be in one PTE fibre, with F̃ monic and d ≥ k.
  Lemma EX gives D − (r − r′) ∈ N·Z[X]_{<d}, so the hypotheses hold with δ = r − r′, ‖δ‖ ≤ 2H. In Lemma R's
  regime, E ≠ 0 for every pair of distinct members (Prop NX).
- *Reading.*
  - Codex's counterexample showed that inherited compatibility can make congruence conditions dependent.
    Prop IH shows that this dependence is never free: it *is* a target relation, computable in polynomial time
    from the pair, with pieces ≤ 2B and small side ≤ 2H·Σ B^l. So, for PTE-type non-exact families, inheritance
    is an instance of Π, not a way around it.
  - Codex proved the existence part (≥ 1 relation, under M > 3k) independently.
  - The bound in (c) is corrected after E3 (§5). The pre-registered version divided by deg D and failed in
    31 pairs.
- *Not covered.*
  - Arbitrary varying multipliers.
  - Identities P − P′ − δ = F̃S + NU whose left side has degree ≥ d. Evaluating at a root keeps F̃(b)S(b)
    (codex).
  - Evaluated (ii) with k ≥ 2. There only D(m) ≡ δ (mod N) is known; by Lagrange this is a k-term *linear*
    congruence among split values, which is ptecons' NC problem.

**Prop NF′ (N-independent inheritance) [PROVED; ptesplit Prop NF applied to differences].**
- *Setting.* Let Π₀ be a set of fibre pairs fixed independently of N, with coefficient bit-sizes ≤ T. Put
  T* = ⌈log₂(2^T + 4H)⌉.
- *Statement.* The balanced semiprimes N = pq ∈ [X, 2X] (p, q > √X/2) for which some pair in Π₀ is non-exactly
  compatible with ‖δ‖ ≤ 2H number at most |Π₀|·k·(4H + 1)·(T*/log₂(√X/2))².
- *Proof.* Compatibility makes N divide a nonzero integer c_l − δ_l (lowest l with c_l ≠ δ_l), which has
  ≤ T*/log₂(√X/2) prime factors above √X/2. ∎
- *Consequence.* For |Π₀|, H, T = e^{o(x)} this is a fraction e^{−x+o(x)} of all N. Inheritance must be built
  from N.

**Prop COH (coherent families are O2 instances) [PROVED, definitional; Remark PEN PROVED].**
- *Statement.* Let G₀..G_e ∈ Z[X] with G_i ≡ ρ_i (mod I), where I = I_N or J_m, and ρ_i ∈ Z[X]_{<k} small.
  Then every P_c = Σc_iG_i satisfies P_c ≡ ρ_c := Σc_iρ_i, with ‖ρ_c‖ ≤ (e+1)‖c‖_∞max‖ρ_i‖.
- *The family.* The members are the c with P_c monic and split with small roots. So e + 1 paid seed congruences
  make every split member of the lattice Z𝒢 (small c) free.
- *What is left.* How dense split polynomials are in a fixed low-rank lattice of polynomials. That is the
  N-free part of ptebarrier's O2 (the missing multi-term rigidity lemma). In addition the basis must be
  N-aligned, which costs e + 1 seeds.
- *Remark PEN (rank 2 with linear pieces).* Take 𝒢 = {T·(X − c), T·(X − c′)}, T split, c ≠ c′ small. Its
  members P_t = T·L_t with L_t = (1+t)X − (c′ + tc) are split for every t (non-monic linear pieces, allowed in
  ptesplit F1).
  - Suppose both seeds are non-exact members: T(μ)(μ − c) ≡ ρ and T(μ)(μ − c′) ≡ ρ′ (mod N), with ρ, ρ′ small
    and ρ − ρ′ a unit.
  - Subtracting gives T(μ)(c′ − c) ≡ ρ − ρ′. Hence (ρ − ρ′)μ ≡ c′ρ − cρ′ (mod N): μ ≡ a₀/b₀ with
    a₀ = c′ρ − cρ′ and b₀ = ρ − ρ′ small.
  - Then every factor-base symbol is a small rational, bμ − a ≡ (ba₀ − ab₀)/b₀. In particular
    L_t(μ) ≡ (c′ − c)(ρ′ + tρ)/(ρ − ρ′), an exact identity of rationals.
  - So member t's relation T(μ)L_t(μ) ≡ ρ′ + tρ is, after cancelling the unit ρ′ + tρ, the single relation
    T(μ)(c′ − c) ≡ ρ − ρ′ for every t. **The pencil has rank 1** (codex's useful-rows caveat).
  - For uniform μ, seeds with |ρ|, |ρ′| ≤ V′ and |c|, |c′| ≤ B exist with probability
    ≤ (4(B+1)V′ + 1)(4V′ + 1)/N.
  - The same holds in View A, with m in place of μ.
- *Rank ≥ 3.* Not covered here. It is the split-side-charged linked-module problem (resbarrier §2.6; ptesplit
  §2.4 joint statement).

### 2.5 What a barrier would have to be

Abundance (2.1) rules out barriers that depend only on the *set* of relations. So a barrier for O1 must bound
the cost of an algorithm class.
- **Query model [PROVED in model; codex].**
  - The model: each candidate's residue is revealed by a query and is uniform and independent, conditional on
    the past.
  - Then R outputs need ≥ RN/W queries in expectation. With W = N^{o(1)} this excludes the target.
  - It does not apply to algorithms that exploit the arithmetic of μ − a.
- **Generic groups [KNOWN].** Exact-target product search in a generic cyclic group reduces discrete logarithms
  (Wagner, CRYPTO 2002, Thm 2, Cor 1). So generic-group speed-ups do not transfer.
- **Linear-lattice model.** Prop LIN, §3.
- An unconditional barrier against all algorithms would be a lower bound for factoring-type relation collection.
  It is not expected from counting or rigidity (assessment).

---

## 3. Construction attempt: lattice reduction (Task 3)

### 3.1 Prop LIN (linearise s roots) [PROVED (model); HEURISTIC H_lin]

**Setting.**
- A prefix multiset C of size M − s gives c = P_C(μ), assumed a unit.
- The suffix S (an s-multiset in [−B, B]) is linearised through its elementary symmetric functions, which makes
  the member condition one linear congruence in (e₁..e_s, v) ∈ Z^{s+1}:
  c·P_S(μ) ≡ w₀ + Σ_i e_iw_i ≡ v (mod N).
- The coset has determinant N.
- The box is |e_i| ≤ C(s,i)B^i, |v| ≤ V, of volume D_s·W with D_s = ∏_i(2C(s,i)B^i + 1). The number of genuine
  suffixes is T_s = C(2B + s, s).

**Model.**
- H_lin: coset points and split points in the box are as frequent as for a uniform affine residue.
- Cost: ≥ 1 per lattice (setup), plus 1 per enumerated point.
- Expected box points p = D_sW/N; split hits h = T_sW/N per lattice.

**Statement.**
- cost per relation ≥ max(1, p)/h = max(N/(T_sW), D_s/T_s).
- Since D_s/T_s ≥ B^{s(s−1)/2}e^{−s²/(2B)} and T_s ≤ (2B + s)^s, minimising over s gives
  cost per relation ≥ (N/W)·exp(−(1 + o(1))√(2 ln(N/W)·ln(2B + s))).
- For W = e^{o(x)} and B = L[β] the saving √(2x ln B) = Θ(x^{(1+β)/2}(log x)^{(1−β)/2}) = o(x). So the cost per
  relation is **N^{1−o(1)}**.
- At the optimum each lattice yields h ≈ B^{−s(s−1)/2} ≪ 1 relation.

**Proof.** Exact arithmetic of the bounds. Optimise max(y − sb, s(s−1)b′/2) over s, with y = ln(N/W),
b = ln(2B + s), b′ = ln B. ∎

**Remarks.**
- *Localising the suffix* to an interval of length ℓ replaces B by ℓ in the formula. The saving √(2y ln ℓ) is
  then smaller.
- *Other linear coordinates* (e.g. Newton) change D_s by e^{O(s²)} = N^{o(1)} at the optimum: the image of the
  cube has volume B^{s(s+1)/2}·c_s with c_s = e^{−Θ(s²)}.
- *Regions not covered.* General non-convex or thin enumeration regions are not covered.
- *Codex's caveat.* In View A with small m, the relaxed box is degenerate: R_s(m) takes ≤ (|m| + B)^s values.
  So H_lin can fail there, and the search is over integers rather than residues.
- *Codex's analysis (eq. 41–48)* gives the same cost and the same conclusion: "small s has inadequate split
  entropy, large s has enormous rejection".

### 3.2 Other lattice mechanisms

- **Meet in the middle [PROVED / HEURISTIC; codex eq. 49–50].** Matching C·u·v mod N ∈ 𝒯 over two product lists
  costs Õ(min{LR, L + WR, R + WL}) + output. The work per relation is ≈ N/(LW) + N/R, so N^{1/2−o(1)}
  amortised. Wagner-type many-list speed-ups do not transfer to generic multiplicative groups (Wagner,
  Thm 2).
- **Coppersmith small roots [KNOWN; codex].** This handles a repeated-root tail C(m − a)^s ≡ v (mod N) for a
  *specified* v with |a| < ½N^{1/s−ε}. It has only 2B + 1 candidates, and is not a lister for split tails.
- **Multiplier lattices [Prop LAT; E2].** LLL returns n independent compressions of one P at height ≈ N^{d/n}.
  In View A they all project to Λ_u (Prop VA). In View B the multipliers are generic: 0 useful vectors with
  deg Q ≥ 1 at d = 2 in 1200 lattices.

---

## 4. (Reserved)

## 5. Experiments against the pre-registration (commit e4e203e)

Resources: one process at a time. Peak RSS: E1 73 MB, E2 50 MB, E3 271 MB. The post-hoc exhaustive count used
1.53 GB (at most). Times: E1 170 s, E2 6 s, E3 3 s.

### 5.1 E1: linearised lattice (`nonexact_E1.jsonl`, `nonexact_E1b.jsonl`)

Balanced semiprimes `balanced_semiprime(n, 1)`; μ uniform. V = 255 for n ≤ 32 and V = 65535 for n = 40.

| n | B | s | T | hits | pred | obs/pred | box pts (subset) obs/pred | cost/hit meas/model | brute/meas | C(2B+s,s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 24 | 8 | 1 | 123 819 | 43 | 80 | **0.54** | — | **1.86** | 9.1 | 17 |
| 24 | 8 | 2 | 13 758 | 69 | 80 | 0.86 | 0/1.6 | 1.16 | 132 | 153 |
| 24 | 8 | 3 | 2 173 | 74 | 80 | 0.92 | 7308/7349 | 1.08 | 1.2 | 969 |
| 24 | 16 | 1 | 63 786 | 65 | 80 | 0.81 | — | 1.23 | 26.8 | 33 |
| 24 | 16 | 2 | 3 753 | 75 | 80 | 0.94 | 14/12.7 | 1.18 | 376 | 561 |
| 24 | 16 | 3 | 322 | 74 | 80 | 0.92 | 464270/464240 | 1.08 | 0.1 | 6545 |
| 32 | 8 | 3 | 523 691 | 395 | 80 | **4.94** | 35/30.5 | **0.23** | 1367 | 969 |
| 32 | 16 | 2 | 904 556 | 92 | 80 | 1.15 | 0/0.05 | 0.87 | 645 | 561 |
| 32 | 16 | 3 | 77 534 | 96 | 80 | 1.20 | 1909/1926 | 0.83 | 41 | 6545 |
| 40 | 8 | 3 | 498 103 | 60 | 80 | 0.75 | 27/32.1 | 1.12 | 269 | 969 |
| 40 | 16 | 2 | 860 359 | 80 | 80 | 1.00 | 0/0.2 | 1.00 | 561 | 561 |
| 40 | 16 | 3 | 73 746 | 82 | 80 | 1.02 | 2018/2025 | 0.97 | 33 | 6545 |

Not run (T > 10⁶; predictions in the JSONL): n = 32 with s = 1 (both B) and (B, s) = (8, 2); the same three for
n = 40. Their model costs per hit are 1.9·10⁵–9.1·10⁷, against brute force 6.3·10⁶ and 1.5·10⁹.

**Evaluation.**
- **P1 (hit rate): FAILED in 2 cells**, (24, 8, 1) at 0.54 and (32, 8, 3) at 4.94. It held in the other 7.
- **Kill signal: formally triggered** by (32, 8, 3).
- **Post-hoc diagnosis (disclosed; `nonexact_E1_distinct_posthoc.jsonl`).**
  - A replay of the identical RNG stream reproduces every hit count. It shows that the hits are repeat detections
    of a few full multisets.
  - Prefixes are drawn with replacement from a small space: 17 symbols, C(23, 7) = 245 157 prefixes for B = 8.
  - In (32, 8, 3) the 395 hits are **4 distinct members**, with maximum multiplicity 205.
  - In (24, 8, 1), 43 hits are 19 members, from 28 expected in the whole space.
  - Cells with large spaces (B = 16) have multiplicities ≤ 7.
  - The registered hit count is therefore not an independent count, and P1's band was the wrong statistic for
    B = 8. This is a design flaw.
- **Exhaustive check (post hoc; `nonexact_E1_exhaustive_posthoc.jsonl`).** Counting all members of the three
  feasible spaces (table in §2.1) gives 31/28.0, 564/584.6 and 4/0.84.
  - The last is a Poisson upper tail with p ≈ 0.01. With 12 cells examined the look-elsewhere probability is
    ≈ 0.1.
  - Its 4 members share ≤ 4 roots pairwise. No pair of symbols μ − a, μ − a′ has a rational ratio of height
    < 20000.
  - We found no structure and read it as a fluctuation. The kill signal is not taken as a positive signal.
- **P2 (box points): held** in all 7 cells in its range: 0.99, 1.00, 0.99, 1.00 (≥ 1000 predicted); 1.15, 0.84
  (30–1000).
- **P3 (cost per hit within [0.6, 1.5]× model): FAILED in the same 2 cells** (1.86 and 0.23), for the same
  reason. Held in 7.
  - Best s per (n, B) row equals the model's arg-min in every row: s = 2 at n = 24, s = 2 at (32, 16), s = 3 at
    (40, 16).
- **P4 (brute/meas ≤ 1.5·C(2B+s, s)): held** in all 9 cells. The lattice saves at most the number of suffixes, a
  subexponential factor.
- **P5 (enumerator): held.** 450 lattices in 9 cells, **0 mismatches** between LLL + exact Fincke–Pohst and the
  exhaustive count. Median nodes/(1 + p) is 0.94–3.45, against ≤ 30 registered.

**Reading [measured].**
- Where the space is large enough for independent sampling, the Prop LIN model predicts hits, box points and
  cost per hit within Poisson error.
- Measured cost per relation at the best s:
  - 70 at n = 24 (brute 2.6·10⁴);
  - 9.8·10³ at n = 32 (brute 6.3·10⁶);
  - 2.3·10⁵ at n = 40 (brute 1.5·10⁹).
- The gain is the suffix count, and the cost per relation still grows like N^{1−o(1)}/W.
- Per-lattice yield at the optimum: 0.02, 10⁻⁴ and 4·10⁻⁶ relations per lattice.

### 5.2 E2: multiplier lattices (`nonexact_E2.jsonl`)

Setup: n ∈ {32, 48, 64}, 100 random split P (M = 12, roots in [−16, 16]) per cell, 21 cells. det L_P = N^d was
verified in 2100/2100 lattices.

**Evaluation.**
- **P6 (median ‖b₁‖/GH ∈ [0.6, 1.4]): held**, 1.01–1.15 in every cell.
- **P7 (median log₂‖r(b₁)‖/log₂N within ±0.12 of d/n): held** in every cell. The observed values are 0.03–0.075
  below d/n, as expected: the r-coordinates are a subset of b₁, and the GH carries the √(n/2πe) factor.
- **P8, first part (useful vectors with deg Q ≥ 1, ≤ 2 per cell): held.** 0 in all 12 d = 2 cells, and 1 in
  (32, d = 1, j = 2).
  - That vector is Q = c(X − a)(X − b) with r constant. Its relation c·(μ − a)(μ − b)P(μ) ≡ r has two integers
    c, r ≈ N^{1/4} to be made smooth: a Dixon-type ledger, consistent with chance.
- **P8, second part (integer-root rate within 3× of baseline): FAILED** in (48, d1, j4), (48, d2, j10) and
  (64, d2, j10), at factors 4.2, 5.4 and 5.5. It also missed narrowly in (32, d1, j4), at factor 3.13.
  - The failures are in the direction of *fewer* roots than baseline.
  - Post hoc (`nonexact_E2_baseline2_posthoc.jsonl`): the registered baseline resampled each coefficient within
    its own magnitude. That overweights zero and small coefficients, and hence roots 0 and ±1.
  - Against a uniform-height baseline (every coefficient uniform in [−‖Q‖, ‖Q‖]) the observed rates agree within
    ~1.5×. Examples: 171/1000 against 0.138; 84/1300 against 0.044; 10/1300 against 0.0088.
- **Kill signal (> 5 useful vectors): not triggered.**

**Reading [measured].** Lattice reduction finds non-exact compressions of a given split P at exactly the
Minkowski/Gaussian height N^{d/n}, n independent ones per lattice. The multipliers behave like random
polynomials of that height: they are not split, and they carry the N-sized ledger.

### 5.3 E3: inheritance (`nonexact_E3.jsonl`)

Setup: 8 cells, 246 constructed pairs. N is a prime dividing D − δ; F̃ is a monic factor of P − r mod N, of
degree max(k, 2) (diagnostic harness).

| (k, M, R) | fibre pairs | tested | mean distinct relations | registered bound (mean) | violations as registered | corrected |
|---|---|---|---|---|---|---|
| (1,4,8) | 33 | 23 | 1.00 | 1.00 | 0 | 0 |
| (1,5,10) | 4 | 4 | 1.00 | 1.00 | 0 | 0 |
| (1,6,11) | 8 | 8 | 1.00 | 1.00 | 0 | 0 |
| (2,4,8) | 3256 | 61 | 4.95 | 5.03 | **15** | 0 |
| (2,5,8) | 593 | 93 | 6.27 | 6.08 | **10** | 0 |
| (2,6,6) | 6 | 6 | 3.67 | 6.33 | **6** | 0 |
| (3,5,8) | 49308 | 21 | 5.24 | 2.14 | 0 | 0 |
| (3,6,6) | 1116 | 30 | 5.17 | 2.87 | 0 | 0 |

**Evaluation.**
- **P9 (a)–(d): held**: 0 violations of membership, congruence, the non-triviality criterion, the |v| bound and
  the piece bound (`nonexact_E3_reasons_posthoc.jsonl`).
- **P9 (e): FAILED in 31 k = 2 pairs.**
  - The registered bound ⌈(n_supp − (k−1))/deg D⌉ forgot that relations are counted up to sign. D(b′) = −D(b)
    makes (D(b), δ(b)) and (D(b′), δ(b′)) the same relation.
  - The corrected bound, with denominator 2·deg D, is a theorem (Prop IH(c)) and has 0 violations. It was found
    post hoc.
- **P10 (k = 1 gives exactly one relation): held**, 35/35 pairs.

---

## 6. Verdict and residual open statement (Task 4)

**What is settled.**
- O1 is **not a separate escape**. Evaluated compressions are exactly the outputs of the listing problem Π (Prop
  EQ). Non-exact ring compressions are a sub-case (with equality at d = 1).
- The exact-case tools fail for a precise reason (Prop BRK). The archimedean rigidity of Lemma R needs a
  difference of height below 2^{M/2}, and a non-exact difference has height ≥ N ≫ 2^{M/2} in the L-range.
- **No counting barrier exists** for O1 (Prop AB). Rigorous counting beyond the trivial bound needs short
  character-sum cancellation far beyond Burgess (Prop CS). The determinant method does not see the congruence
  (Prop BP-NO).
- The "dependent conditions" gap is closed for PTE-type families with d ≥ k: inheritance is itself a target
  relation (Prop IH), and is absent for N-independent fibres (Prop NF′).
- Coherent families reduce to O2 (Prop COH). Rank-2 pencils with linear pieces have rank 1 (Remark PEN).
- Lattice reduction in the linearised form costs N^{1−o(1)} per relation (Prop LIN, model; E1). Multiplier
  lattices give Minkowski-size generic multipliers (Prop LAT, VA; E2).

**Residual open statement Π_β.** Fix β < 1/3 and put ℓ = x^{2β}(log x)^{1−2β}. Is there a uniform classical
algorithm with all of the following properties?
- Its input is N alone and its running time is L[β]^{1+o(1)}.
- It outputs L[β]^{1+o(1)} congruences ∏_i c_i ≡ v (mod N).
- Each piece c_i lies in a factor base of size L[β]: integers ≤ L[β], or values m − a (|a| ≤ L[β]) at a computed
  m.
- Either 0 < |v| ≤ e^{O(ℓ)}, or v is a product of pieces.
- The exponent rows have rank ≥ (1 − o(1))·#columns.

**Mechanism classes and how each fares.**
- Incoherent members pay one congruence each: N/W per member in the query model, N^{1−o(1)} in the
  linear-lattice model.
- Coherent families reduce to O2.
- PTE inheritance yields Π-outputs.

Π_β is ptecons' open item 3 and smoothenum's algebraic/doubly-smooth listing primitive. Its doubly split form
would give heuristic polynomial time (smoothenum T7/T8). That strength is itself evidence against existence
(assessment).

**Probability.** < 1% that O1 yields sub-L[1/3] (ours); about 5% (codex, blind). ptebarrier, ptesplit and ptecons
all stood at < 1%.

---

## 7. Literature (to be checked by the literature critic)

- **Lattice-point counting.**
  - Minkowski's second theorem.
  - Henk, "Successive minima and lattice points", Rend. Circ. Mat. Palermo Suppl. 70 (2002) 377–384, Thm 1.5.
  - Barroero–Widmer, IMRN 2014(18) 4932–4957, Thm 1.3.
  - Davenport, J. London Math. Soc. 26 (1951).
- **Character sums.**
  - Burgess, Proc. LMS 12 (1962) / Mathematika 4 (1957).
  - Soundararajan, "The distribution of smooth numbers in arithmetic progressions", Anatomy of Integers, CRM
    Proc. 46 (2008) 115–128, Thm 1.
- **Determinant method.**
  - Bombieri–Pila, Duke 59 (1989), Thm 5.
  - Heath-Brown, Ann. Math. 155 (2002).
  - Castryck–Cluckers–Dittmann–Nguyen, ANT 14 (2020).
  - Binyamini–Cluckers–Novikov, arXiv 2302.04209.
- **Lattice factoring and small roots.**
  - Schnorr, EUROCRYPT '91; Adleman 1995 (unpublished).
  - Ducas 2021 experiments.
  - Coppersmith, J. Cryptology 10 (1997) 233–260, Thm 1.
  - LLL, Math. Ann. 261 (1982).
- **Generalized birthday.** Wagner, CRYPTO 2002, LNCS 2442, 288–304, Thm 2, Cor 1.
- **Kleiman, Rees–Smyth.** As in ptebarrier §10.
- **Linear sieve.** Schroeppel (unpublished; Pomerance 1982).

## 8. Reviews

- **Codex blind round (`nonexact_1.out`).** Same verdict, with no barrier and no construction. Its contributions,
  all adopted and credited above:
  - the Smith-form determinant;
  - Henk Thm 1.5 and Barroero–Widmer;
  - conditional abundance under CS;
  - the query-model bound;
  - the Wagner obstruction;
  - the existence form of Prop IH;
  - Zariski density;
  - the View A degeneracy of the relaxed box;
  - the MITM and Coppersmith remarks;
  - the containment (51);
  - the useful-rows caveat;
  - the residual statement.

  Its probability is about 5%.
- Codex adversarial round: pending.
- Opus proof critic: pending.
- Opus literature critic: pending.

## 9. Reproduction

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_nonexact.py
PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment E1   # ~3 min
PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment E2   # ~6 s
PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment E3   # ~3 s
PYTHONPATH=. .venv/bin/python -m src.nonexact.posthoc             # replay: distinct members (~2 min)
PYTHONPATH=. .venv/bin/python -m src.nonexact.posthoc exhaustive  # exhaustive member counts (1.5 GB peak)
# E2 baseline2 and E3 reasons: e2_baseline2(load_cfg()), e3_reasons(load_cfg()) in src/nonexact/posthoc.py
```
