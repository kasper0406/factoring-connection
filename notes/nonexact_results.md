# nonexact: O1, non-exact and evaluated compressions modulo N

Agent `nonexact`, phase 17, 2026-10-20. **Revised** after three reviews, all points applied (§8):
- codex nonexact_2: MAJOR REVISION;
- the Opus proof critic: ACCEPT WITH FIXES;
- the Opus literature critic: OK WITH FIXES.

**Files.**
- Code: `src/nonexact/`.
  - `core.py`: exact helpers.
    - Reductions mod (N, F̃), flint LLL, kernel lattices and exact Fincke–Pohst enumeration.
    - The Prop LIN counts and the Prop IH extractor.
    - Prop IH helpers in two versions. `ih_distinct_lower_bound` is the registered one, refuted and kept for
      history. `ih_hypothesis`, `ih_good_points_bound` and `ih_distinct_bound` are the corrected ones.
  - `run_experiment.py`: the pre-registered runs E1–E3.
  - `posthoc.py`: disclosed post-hoc diagnostics and rescoring.
- Tests: `tests/test_nonexact.py`, 20 tests, including the reviewers' counterexamples as regressions.
- Pre-registration: `notes/nonexact_prereg.md`, commit e4e203e, made before any registered run.
- Config: `experiments/configs/nonexact.json`.
- Raw results:
  - pre-registered: `experiments/raw/nonexact_E{1,1b,2,3}.jsonl`;
  - post hoc: `nonexact_E1_distinct_posthoc.jsonl`, `nonexact_E1_exhaustive_posthoc.jsonl`,
    `nonexact_E2_baseline2_posthoc.jsonl`, `nonexact_E3_reasons_posthoc.jsonl` and `nonexact_rescore_posthoc.json`.
- Codex (GPT-6.1-sol): blind round `codex/nonexact_1.{md,out}`; adversarial round `codex/nonexact_2.{md,out}`.

**Labels.**
- PROVED.
- PROVED under GRH.
- PROVED (model): an exact consequence of a named model.
- KNOWN: cited.
- HEURISTIC: the assumption is named.
- measured · open · assessment.

Novelty is never claimed.

**Notation**, as in `notes/ptesplit_formulation.md` and ptebarrier §1.1.
- x = log N; B = L[β] with 0 < β < 1/3; M ≍ x/log B.
- ℓ := log V = O(x^{2β}(log x)^{1−2β}), and W := 2V + 1.
- P_A = ∏_{a∈A}(X − a). The presentation F̃ is monic of degree d, and I_N := (N, F̃).
- View A: φ(X) = m ∈ Z with N | F̃(m), and J_m := (N, X − m).
- View B: φ(X) = μ ∈ Z/N.

---

## 0. Verdict

**Plain language.**
- We find **neither a barrier nor a construction**.
- Instead we show four things.
  - **Reformulation.** The reduction may hold only modulo N, or only in value at m. Then the polynomial
    apparatus stops constraining the relations: presentation, small remainder and PTE fibre alike.
    - The evaluated version, with admissible (split) multipliers, is a product-listing problem Π. Π asks for many
      products of factor-base pieces that are congruent to something small modulo N.
    - The polynomial structure can matter only through how members are generated.
  - **No counting barrier for natural catalogues.**
    - Unconditionally, any S distinct split products contain at least S²/N − S congruent pairs. Each such pair is
      a relation that is split on both sides.
    - Under GRH, natural catalogues of exponential size contain (1+o(1))·S·W/N one-sided members.
    - So the rigidity and counting methods that closed the exact case cannot close O1 for such catalogues. Any
      obstruction must concern the cost of *finding* the members.
  - **In the models we can analyse, lattice reduction is exponential.**
    - Linearising the last s roots and enumerating the whole box costs at least
      (N/W)·exp(−(1+o(1))√(2 ln(N/W)·ln(2B+s))) = N^{1−o(1)} per relation. This is a model statement, and the
      toy-scale measurements match the model.
    - Multiplier lattices reach Minkowski heights, but their multipliers are generic and not split.
  - **Inheritance is computably relation-valued.**
    - Take two non-exact members of one PTE fibre with d ≥ k, compatible with each other. Their compatibility
      converts into congruences of the target type (Prop IH, under an explicit smallness hypothesis).
    - This answers ptebarrier's *conversion* question. It does **not** prove that inherited conditions are
      independent.
- **Not excluded.** No known construction of any of the following exists, but none is ruled out either:
  - targeted (non-exhaustive) or batched lattice searches;
  - nonlinear N-aligned families;
  - higher-rank coherent modules (O2, resbarrier §2.6);
  - arithmetic that uses the additive structure of the pieces μ − a;
  - preprocessing shared across prefixes.

**Proved (main items).**
- **Lemma COL / Prop EQ.**
  - Evaluated O1 with h = 1, or with h split and absorbed into A, equals one-sided Π with 𝔛 = {m − a}.
  - Non-exact O1 is contained in Π for k = 1.
  - For k ≥ 2 it is contained in Π with V′ = kH|m|^{k−1} (View A), or in the doubly-split form when r splits
    (View B).
  - Non-exact members with general, non-split r are a two-sided, NFS-like variant. They are not covered.
- **Prop BRK.**
  - In (i), Lemma EX forces non-exact differences of height at least N − 2H. Lemma R then needs
    M ≥ 2 log₂N ≈ 2.885x, while M ≍ x/log B = o(x).
  - In (ii), J_m contains X − m, so a small value gives no height bound.
- **Prop AB** [PROVED] and **AB′** [PROVED under GRH]: the abundance statements above.
- **Prop IH** [PROVED, assuming (H): ‖δ‖Σ_{l<k}B^l < N].
  - A compatible pair has at least M − 2k + 3 good points.
  - It has at least ⌈(M − 2k + 3)/(2 deg D)⌉ distinct relations up to sign, with pieces ≤ 2B and
    |v| ≤ ‖δ‖ΣB^l.
  - Without (H) there are counterexamples, including one with no relation at all.
- **Prop NF′.** N-independent fibres inherit compatibility for at most an e^{−x+o(x)} fraction of balanced
  semiprimes.
- **Prop COH.** Linearly coherent families are instances of O2.
- **Remark PEN.** For the pencils T(X − c), T(X − c′), two seeds force μ to be a small rational mod N. In View A
  there is then a dichotomy:
  - either the pencil has rank 1 in the prime ledger;
  - or its pieces violate target (a).

  The draft's unconditional "rank 1" is withdrawn.
- **Prop LIN** (exhaustive-enumeration model), **Prop LAT** (actual determinant D, Henk upper count) and **Prop VA**
  (View A value lattice).

**Experiments (pre-registered; the failures are listed here in full).**
- **E1, linearised lattice.** 12 of 18 cells were run.
  - P1 failed in 2 cells, at 0.54 and 4.94, and held in 10. The **kill signal was triggered** at (32, 8, 3).
  - Under the frozen statistic, P3 failed in 4 cells and held in 8. P2, P4 and P5 held.
  - LLL + Fincke–Pohst matched the exhaustive counts in 450/450 lattices.
  - **(32, 8, 3).** Given the catalogue, the 395 detections are explained: the sampler's weights predict 396.7.
    The catalogue itself has 4 members against 0.84 expected; that excess is unexplained, and no structure was
    found.
  - **(24, 8, 1)** is **unexplained**: 43 detections against a weighted expectation of 70.7. It is not a bug, and
    five fresh seeds give 61–83.
- **E2, multiplier lattices.** 2100 lattices.
  - P6 and P7 held.
  - There were 0 useful vectors at d = 2.
  - The P8 root-rate comparison **failed in 5 cells**, all with fewer roots than the registered baseline. Post
    hoc, that baseline turned out to be miscalibrated.
- **E3, Prop IH.** 246 pairs.
  - P9(a)–(d) and P10 held.
  - The registered distinct-count bound **failed in 31 pairs**.
  - The corrected IH has 0 violations in the 201 pairs that satisfy (H) (post hoc).

**Assessment.** Probability below 1% that O1 yields sub-L[1/3] (ours; subjective). Codex's blind estimate is
about 5%. The experiments do not calibrate either figure.

**Residual open statement.** Π_β with its admissibility conditions (§6). O1 is reformulated as Π_β (evaluated)
or contained in it (non-exact, with the scope above).

---

## 1. Formalisation (Task 1)

### 1.1 The two versions

A member is (A, h, Q, r), with relation φ(h)φ(P_A) = φ(r).

- **(E) Exact.** hP_A − r ∈ F̃Z[X]. This case is closed in its regime by ptebarrier Theorem W, Corollary P and
  Prop S.
- **(O1-i) Non-exact.** hP_A − r ∈ I_N \ F̃Z[X].
  - Equivalently, hP_A − F̃Q − r = N·T for some certificate (Q, T) (codex, eq. 5).
  - Certificates are not unique: (Q, T) ↦ (Q + NS, T − F̃S).
  - The canonical Euclidean quotient gives red_F̃(hP_A − r) = N·D_A with 0 ≠ D_A ∈ Z[X]_{<d} (Lemma EX).
  - The member is d scalar conditions mod N. Bounds on Q are extra restrictions.
- **(O1-ii) Evaluated.** hP_A − r ∈ J_m with |r(m)| ≤ V.
  - There is no height condition.
  - It is non-trivial when h(m)P_A(m) ≠ r(m).
  - J_m = {g : g(m) ≡ 0 mod N}, so this is a single scalar condition.

**Relation between the two [PROVED].** I_N ⊆ J_m, strictly when d > 1.
- So a non-exact member satisfies the *congruence* of (ii).
- It is an evaluated member only if, in addition, |r(m)| ≤ V. In general |r(m)| ≤ kH|m|^{k−1}, and |m| may be as
  large as N^{1/d}.
- The converse fails. Codex's example: N = 77, m = 10, F̃ = X² − 23, P = (X−1)(X−2), r = −5. Here P(10) ≡ −5, but
  red_F̃(P − r) = 30 − 3X.
- (i) excludes exact members by definition. So "exact ⊂ (i)" holds only for the sets of relations.

### 1.2 Lemma COL, the problem Π, Prop EQ

**Lemma COL [PROVED, elementary].**
- (a) An evaluated member is exactly an integer congruence h(m)P_A(m) ≡ v (mod N) with |v| ≤ V. Replacing r by
  the constant v keeps the condition.
- (b) Conversely, any A and v with P_A(m) ≡ v is a member on F̃ = X − m with h = 1, r = v. At m = 0 this is
  ∏(−a) ≡ v (codex, eq. 51).

**The multiplier must be charged.** If h does not split over the pieces, then h(m) is an arbitrary integer.
- Example: N = 77, m = 0, A = {−2}, h = 39. Then 39·2 ≡ 1 (mod 77), but 2 ≢ 1.
- Such a member is a point of the value lattice with an unfactored multiplier (Prop VA), not an output of Π.

**Π(N; 𝔛, V).**
- 𝔛 = {x_a} is a set of at most L[β] units mod N computable from N: small integers, m − a, or μ − a. A non-unit is
  output as a factor.
- The outputs are distinct multisets A with ∏_{a∈A}x_a ≡ v (mod N). There are two forms:
  - one-sided: 0 < |v| ≤ V, with the product ≠ v as integers;
  - doubly split: v = ∏_{a′∈A′}x_{a′} with A′ ≠ A and the two sides unequal as integers. This excludes exact
    integer identities such as 2·6 = 3·4.

**Prop EQ [PROVED].**
- (i) Evaluated O1 with h = 1, or with h split over 𝔛 and absorbed into A, equals one-sided Π with 𝔛 = {m − a}.
- (ii) Non-exact O1 is contained in:
  - one-sided Π with V′ = kH·max|m|^{k−1} (View A);
  - doubly split Π when r splits over 𝔛 (View B).
  - For k = 1 no extra condition is needed.
  - General non-split r is not covered.
- (iii) Pieces of size at most B^{O(1)} can be re-encoded at m = 0 (codex).

**Scope.**
- The cost of the presentation must be charged (codex).
- Π ⊆ (i) needs d = 1.
- **Π is a listing problem, not a sufficient factoring primitive.** Factoring also needs a column convention,
  fruitful F₂-dependencies, smoothness tests and linear algebra (§6).

### 1.3 Which step of Lemma R / Theorem W breaks (Prop BRK) [PROVED]

**(i) Non-exact.**
- Take canonical defects and a fibre with d ≥ k. Then G = P − P′ = (r − r′) + N·E with 0 ≠ E ∈ Z[X]_{<k}, so
  ‖G‖ ≥ N − 2H.
- Lemma R's threshold then needs M ≥ 2 log₂N − O(1) = 2x/ln 2 − O(1). But members need only
  M log(2B) ≥ x − 1, and M ≍ x/log B = o(x).
- Step 4 is vacuous as soon as log₂Δ ≥ deg U′. This already holds for M ≤ log₂N.
- **Theorem W.** With canonical reduced defects T, T′, the identity is (codex, eq. 10)
  Q′hP − Qh′P′ = Q′r − Qr′ + N(Q′T − QT′).
  - When Q′T − QT′ ≠ 0, the right-hand side has height at least N in its low coefficients.
  - Canonical quotients can be small; codex gives F̃ = P − r + N, with Q = 1. What fails is the height of the
    right-hand side.
  - ptebarrier's N = 25 example gives two non-exact members with different Ψ. This illustrates the mechanism
    only; it is outside the L-range.
- **What survives.** The equality of canonical quotients in Corollary P, and Prop NX.

**(ii) Evaluated.**
- J_m ∩ Z[X]_{<k} contains X − m, of height |m| < N. Lemma EX fails exactly because deg(X − m) = 1 < k.
- r = C(X − m) + v has value v for every C.
- So fibre differences are constrained only by G(m), which gives no Mahler bound.
- For k = 1, *pairwise fibre compatibility* is a congruence on a constant, and NX1 applies.

**Scale (codex).** ℓ/M = O(x^{3β−1}(log x)^{2−3β}) → 0.

---

## 2. Barrier attempts (Task 2)

### 2.1 Abundance

**Prop AB [PROVED; standard Cauchy–Schwarz].**
- *Statement.* S distinct objects with residues mod N give at least S²/N − S ordered congruent pairs.
- *Catalogue.* Products of exactly M primes from [B^{1/4}, B] are distinct, and there are N^{c−o(1)} of them when
  Mb = cx (codex, eq. 31).
- *Consequence.* With S ≥ 2NL there are at least S(2L − 1) pairs of distinct integers, each split on both sides.
- *What it does not give.* One-sided small residues; independent or fruitful rows; abundance inside restricted
  families. A family can avoid the window altogether.

**Prop AB′ [PROVED under GRH; codex's argument, rechecked].**
- *Hypothesis CS.* |Σ_{p∈[B^{1/4},B]}χ(p)| ≤ B^{1/2}x^{C₀} for every nonprincipal χ mod N.
  - Under GRH, Σ_{p≤y}χ(p) ≪ y^{1/2}log²(Ny), so CS holds with C₀ = 2 (Iwaniec–Kowalski Ch. 5; Montgomery–Vaughan
    Ch. 13; literature critic).
- *Statement.* For c > 2 and V = N^{o(1)}, the catalogue has (1+o(1))·S·W/φ(N) members with 0 < |v| ≤ V and
  u ≠ v.
- *Proof.*
  - Use the generating function exp(Σ_j z^j/j·Σ_pχ(p)^j). It is O(1) on |z| = B^{−1/2}x^{−C₀−2}.
  - Hence |H_M(χ)| ≤ e^{Mb/2+O(M log x)}.
  - Since S ≥ e^{Mb−O(M log x)}, the relative error per class is φ(N)·N^{−c/2+o(1)} = o(1). ∎
- *Unconditional status.* The available unconditional theorems need q ≤ y^{4√e−ε}: Soundararajan 2008, Thm 1, and
  Harper 2012.

**Measured, post hoc (E1).** Exhaustive member counts: 31 vs 28.0, 564 vs 584.6, and 4 vs 0.84 (§5.1).

**Reading.**
- Natural catalogues have no scarcity obstruction: unconditionally in the doubly split form, under GRH in the
  one-sided form.
- So a barrier must bound the cost of finding members.
- This does not exclude counting barriers for specific restricted families.

### 2.2 Lattices

**Prop LAT [PROVED; determinant from codex, eqs. 14–16].**
- *Setting.* L_P = {(Q, r) : deg Q ≤ j, deg r < k ≤ d, QP − r ∈ I_N}, of dimension n = j + 1 + k.
- *Determinant.* D = ∏N/gcd(N, δ_i), where the δ_i are the Smith invariants. For N = pq this is
  p^{rank_p}q^{rank_q}.
  - D = N^d if j + 1 ≥ d and P is a unit. E2 confirms this in 2100/2100 lattices.
  - D can be smaller. Codex: P = F̃ = X², N = 77, j = 0, k = 1 gives D = 77.
  - Full column rank mod p, with n ≤ d, excludes every nonzero vector of height below p (codex).
- *Minkowski [KNOWN].* There is a nonzero coefficient vector of height at most D^{1/n}. It is not necessarily a
  non-exact relation or an admissible multiplier.
- *Henk 2002, Thm 1.5 [KNOWN].* #(L ∩ K) < 2^{n−1}∏⌊2/λ_i + 1⌋.
  - With λ_i relative to [−1, 1]^n, R ≥ λ_n and Minkowski II, this gives
    #(L_P ∩ [−R, R]^n) < 2^{n−1}3^n n!·R^n/D.
  - This is an **upper bound only**.
- *Barroero–Widmer (IMRN 2014, Thm 1.3).* Gives an asymptotic in terms of Euclidean minima. The constants are not
  explicit, and growing dimension is uncontrolled.
- *Short vectors.* The determinant gives no lower bound on individual minima; Z^{n−d} × NZ^d is an example.
  Compressions with multipliers exist for every P, but the multipliers are generic.

**Prop VA [PROVED].**
- *Statement.* (Q, r) ↦ (Q(m), r(m)) maps L_P into Λ_u = {(q, v) : qu ≡ v}, u = P(m), of determinant N.
  The kernel is trivial vectors, e.g. multiples of X − m.
- *Exact average for unit q* (for example, 0 < |q| < min p | N):
  Σ_{u mod N}#{(q, v) ∈ Λ_u : 0 < |q| ≤ Q₀, |v| ≤ V} = 2Q₀(2V+1).
  - For general q the weight is gcd(q, N)(2⌊V/gcd⌋ + 1) (codex).
  - This is average scarcity, the Dixon/CFRAC ledger (Dixon 1981; Morrison–Brillhart 1975).
  - It does not show that polynomial structure cannot help to locate points.
  - Vallée (1991) is the closest prior use of determinant-N two-dimensional lattices.

### 2.3 Charging the split structure

**Prop CS [PROVED inequality; the readings are assessments].**
- *Setting.* I = {μ − a : |a| ≤ B} ⊂ F_p^*, n = |I| < p, and R_M(c) = #{x ∈ I^M : ∏x_i ≡ c (mod p)}.
- *Statement.* For c ≠ 0 and M ≥ 2: |R_M(c) − n^M/(p−1)| ≤ n·S_max^{M−2}, where
  S_max = max_{χ≠χ₀}|Σ_{x∈I}χ(x)|.
- *What this route needs.* It certifies the heuristic count only if S_max ≤ n(n/p)^{1/(M−2)}.
  - With M ∼ c·x/log B and p = N^{1/2+o(1)}, this is S_max ≲ n^{1−1/(2c)+o(1)}.
  - That is square-root cancellation at c = 1, and power saving in general.
  - It is what this argument needs, not a necessary condition.
- *Known limits.*
  - Burgess (1957, 1962, 1963) needs length at least p^{1/4+ε}, or q^{1/4+ε} for cube-free q.
  - Under GRH, initial segments are handled only when log H/log log q → ∞ (Granville–Soundararajan 2001). As far
    as we know, nothing covers arbitrary shifts at length q^{o(1)}.
- *Integer-lifting bounds beat the trivial count.*
  - Cilleruelo–Garaev (GAFA 2011) give R_M(c) ≪ n^{M−2+o(1)} for n < p^{1/4}, and n^{M−3+o(1)} for n ≪ p^{1/8}.
  - BGKS (2013, 2014) handle ν variables for h < p^{1/(ν²−1)}.
  - The draft's "only the trivial bound" is **withdrawn**.
- *Reading.*
  - These bounds save at most n^{O(1)}.
  - In the target regime W = e^{Θ(ℓ)} exceeds every fixed power of n, because ℓ/log B → ∞. So member-fraction
    bounds of the form W·n^{−O(1)} are vacuous there.
  - No known method counts split members in the target window rigorously.
  - Rigorous upper counts would in any case bound over-density, not search cost.

**Prop BP-NO [assessment: the cited theorems do not directly give the needed bound].**
- (a) *Zariski density.* Integer-split coefficient vectors are Zariski-dense as the root bound varies, since the
  map has Vandermonde Jacobian (codex).
  - In a box, this holds only against degrees below about (2B+1)/M.
  - So BP 1989 and CCDN 2020, which count points on a specified variety, do not apply as stated.
- (b) *Lifting the congruence.* The lift Y: Nt + v − ∏(μ̃ − a_j) = 0 is ruled by lines, and its bounded integral
  points are exactly the congruence solutions.
  - This does not show that geometric methods are useless here. Modular small-box methods combine divisibility
    of determinants with archimedean bounds (Chang et al., arXiv 1111.1543).
  - Heath-Brown 2002 counts points of varieties over Q.
- (c) *Structured families.* A structured family on a fixed irreducible curve of degree D in A^n has at most
  c(n)D²B^{1/D}(log B)^κ members (BCN Thm 2, eq. 27).
  - Here c(n) depends on n = M, so the bound applies only after projecting to bounded dimension (literature
    critic).
- **Withdrawn:** "large structured families must be linearly coherent".

### 2.4 Inheritance

**Prop IH [PROVED; hypothesis (H) added after review].**
- *Hypotheses.*
  - P ≠ P′ are monic, split, of degree M, with roots in [−B, B].
  - D = P − P′ has degree < k.
  - D ≡ δ coefficientwise mod N, with δ ∈ Z[X]_{<k} and D ≠ δ.
  - **(H)** ‖δ‖_∞Σ_{l<k}B^l < N.
  - Put E = (D − δ)/N ≠ 0.
- *(a) [no (H) needed].* At a root b of exactly one of P, P′: u_b = ±∏(b − a) ≡ v_b = ±δ(b), and
  u_b − v_b = ±N·E(b). So u_b ≠ v_b if and only if E(b) ≠ 0. This step is tautological.
- *(b) [needs (H)].* There are at least M − 2k + 3 good points.
  - The signed cofactor measure has vanishing moments 0..M−k, so by Vandermonde it has support of size at least
    M − k + 2.
  - Support points that are roots of both polynomials have δ(b) ≡ 0, hence δ(b) = 0 by (H), hence E(b) = 0.
  - E has at most k − 1 roots.
  - Without (H), only M − 3k + 4 good points are safe.
- *(c) [needs (H)].* Distinct relations, counting (u, v) and (−u, −v) as one:
  - exactly one if deg D = 0, because (H) forces δ to be constant;
  - at least ⌈(M − 2k + 3)/(2 deg D)⌉ otherwise, since |D(b)| = t has at most 2 deg D solutions.
  - Existence therefore needs M ≥ 2k − 1. Codex's example P = X⁵, P′ = (X−1)⁵, N = 5 has no relation.
- *(d) Sizes.* Pieces are at most 2B, and |v_b| ≤ ‖δ‖ΣB^l. If 2B < min p | N, then v_b ≠ 0.
- *Application.*
  - For fibre members with F̃ monic and d ≥ k, Lemma EX gives δ = r − r′, with ‖δ‖ ≤ 2H.
  - (H) holds when 2HΣB^l < N.
  - The relations have target size when log‖δ‖ + (k−1)log B + log k = O(ℓ).
  - In regime, E ≠ 0 (Prop NX).
- *Counterexamples without (H)* (regression tests):
  - critic: k = 2, A = (−8,−8,−5,−5), A′ = (−8,−7,−7,−4), N = 6, δ = 2 − 2X;
  - critic: k = 3, N = 4, with 0 good points;
  - codex: k = 4, (1,1,1,1,2,3) vs (0,1,2,2,2,2), N = 101, where both relations are trivial.
- *Reading.*
  - IH answers the *conversion* question: an inherited compatibility is computably a target relation, i.e. an
    instance of Π.
  - It does **not** give independence across a family or a per-member search cost. "Closes the gap" is
    withdrawn.
  - Codex proved existence independently for M > 3k.
- *Not covered.*
  - Varying multipliers.
  - Identities of degree at least d.
  - Evaluated (ii) with k ≥ 2, where only D(m) ≡ δ(m) holds. That is ptecons' NC problem.

**Prop NF′ [PROVED].**
- *Statement.* Let Π₀ be N-independent, with bit-sizes at most T and T* = ⌈log₂(2^T + 4H)⌉. At most
  |Π₀|·k(4H+1)(T*/log₂(√X/2))² balanced semiprimes in [X, 2X] admit non-exact compatibility with ‖δ‖ ≤ 2H.
- *Scope.* This is an almost-all statement for explicitly bounded catalogues only.

**Prop COH [PROVED, definitional].**
- *Statement.* Let G_i ≡ ρ_i (mod I) with ρ_i small. Then P_c = Σc_iG_i ≡ Σc_iρ_i.
- *Consequence.* The family is the set of split members of a polynomial lattice of rank at most e + 1. This is
  ptebarrier's O2.
- *Scope.* The cost of the e + 1 seeds is not bounded here. COH classifies these stipulated linear families
  only.

**Remark PEN [PROVED; the draft's rank claim corrected].**
- *Setting.* Seeds T(μ)(μ − c) ≡ ρ and T(μ)(μ − c′) ≡ ρ′, with ρ − ρ′ a unit.
- *Alignment.* They force (ρ − ρ′)μ ≡ c′ρ − cρ′, i.e. μ ≡ a₀/b₀ (mod N).
  - For uniform μ this has probability at most (4(B+1)V′ + 1)(4V′ + 1)/N.
- *View A dichotomy.*
  - (1) If m = a₀/b₀ exactly, then L_t(m) = (c′ − c)(ρ′ + tρ)/(ρ − ρ′) as integers. Every member then has the
    seed's prime-exponent row, i.e. rank 1 (tested).
  - (2) Otherwise |m| ≥ (N − |a₀|)/|b₀|, and the pieces have size ≳ N/V′, outside target (a).
    Codex's example N = 77, μ = 75 is of this kind; its rows are independent.
- *View B.* Only the alignment is claimed. "Rank 1" is **withdrawn**.

### 2.5 What a barrier would have to be

- **Query model** [PROVED in the model; codex]. R outputs need at least R·N/W queries. This does not cover
  arithmetic algorithms.
- **Generic groups** [KNOWN]. Wagner, CRYPTO 2002, Thm 2 (credited there to W. Dai) and Cor. 1; Nechaev 1994;
  Shoup 1997. A k-sum solver gives discrete logarithms, and generic algorithms need Ω(√p).
- **The linear-lattice model.** Prop LIN.
- A barrier against all algorithms would be a lower bound for a natural class of factoring algorithms (assessment).

---

## 3. Construction attempt: lattice reduction (Task 3)

### 3.1 Prop LIN [PROVED (model)]

**Setting.**
- Take a prefix C of M − s roots, so c = P_C(μ) is a unit.
- Linearise the suffix: c·P_S(μ) ≡ w₀ + Σe_iw_i ≡ v (mod N). This is an affine lattice of determinant N.
- The box is |e_i| ≤ C(s,i)B^i and |v| ≤ V. Put D_s = ∏(2C(s,i)B^i + 1) and T_s = C(2B+s, s).

**Model.**
- Each prefix lattice is enumerated **exhaustively** at cost 1 + #points.
- **H_lin:** uniform-residue frequencies.

**Statement.**
- The cost per relation is at least max(N/(T_sW), D_s/T_s).
- Using D_s/T_s ≥ B^{s(s−1)/2}e^{−s²/(2B)} and T_s ≤ (2B+s)^s, this is at least
  (N/W)·exp(−(1+o(1))√(2 ln(N/W)·ln(2B+s))).
- The bound is optimised at s ∼ √(2 ln(N/W)/ln B). It equals N^{1−o(1)} for W = e^{o(x)}.
- Each lattice yields about B^{−s(s−1)/2} relations.

**Proof.** The two inequalities are elementary; the proof critic also checked them numerically for B < 200 and
s < 40. Then optimise over s. ∎

**Scope (codex).** The model needs exhaustive enumeration. It does **not** exclude:
- targeted CVP searches;
- partial enumeration that exploits split structure;
- thin regions or exceptional prefixes;
- batching or hybrid strategies.

A change of coordinates changes D_s by a factor e^{O(s²)}; this is absorbed in the 1 + o(1). In View A with small
m, H_lin can fail.

**Comparison with Schnorr and Adleman.**
- Schnorr (EUROCRYPT '91; ePrint 2021/232) and Adleman (1995) linearise *logarithmically*, in a prime-number
  lattice.
- LIN linearises *additively*, in elementary symmetric coordinates.
- Both impose the congruence linearly and lose the fine condition.
- The Schnorr line fails in practice: Ducas 2021; Aboumrad–Widdows–Kaushik (arXiv 2308.07804).

### 3.2 Other mechanisms

- **Meet in the middle** (codex). Amortised cost N^{1/2−o(1)} per relation.
- **Coppersmith** (J. Cryptology 1997, Thm 1). Finds a specified-v repeated-root tail; it is not a lister.
- **Boneh's CRT decoding** (STOC 2000; JCSS 2002). Gives integers in a short interval with a large smooth *part*.
  smoothenum's Prop BC limits its gain to L^{o(1)}.
- **Multiplier lattices** (Prop LAT, E2).
  - Reaching height e^{O(ℓ)} needs n ≥ d·x/ℓ, i.e. multipliers of degree ≍ x^{1−2β}(log x)^{2β−1}.
  - None of the 9 600 d = 2 basis vectors was useful (n ≤ 13).

---

## 4. (Reserved)

## 5. Experiments against the pre-registration (commit e4e203e)

**Resources.**
- Peak RSS: E1 73 MB, E2 50 MB, E3 271 MB, and 1.53 GB for the post-hoc exhaustive count.
- Run times: 170 s, 6 s and 3 s.
- The full repository suite: 893 passed, plus one transient collection error in another agent's
  `test_coverage38.py`, which passes on its own.

### 5.1 E1: linearised lattice

**Setup.**
- 18 cells, of which 12 were run.
- V = 255 for n ≤ 32 and 65535 for n = 40.
- The 6 cells not run have model costs per hit of 3.9·10⁴ to 3.7·10⁵. Brute force costs 6.3·10⁶ at n = 32 and
  6.0·10⁶ at n = 40.

| n | B | s | T | hits | pred | obs/pred | box obs/pred | P3 frozen ÷ model | brute/meas | C(2B+s,s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 24 | 8 | 1 | 123 819 | 43 | 80 | **0.54** | — | **1.86** | 9.1 | 17 |
| 24 | 8 | 2 | 13 758 | 69 | 80 | 0.86 | 0/1.6 | 1.35 | 132 | 153 |
| 24 | 8 | 3 | 2 173 | 74 | 80 | 0.92 | 7308/7349 | 1.08 | 1.2 | 969 |
| 24 | 16 | 1 | 63 786 | 65 | 80 | 0.81 | — | 1.23 | 26.8 | 33 |
| 24 | 16 | 2 | 3 753 | 75 | 80 | 0.94 | 14/12.7 | **2.02** | 376 | 561 |
| 24 | 16 | 3 | 322 | 74 | 80 | 0.92 | 464270/464240 | 1.08 | 0.1 | 6545 |
| 32 | 8 | 3 | 523 691 | 395 | 80 | **4.94** | 35/30.5 | **0.30** | 1367 | 969 |
| 32 | 16 | 2 | 904 556 | 92 | 80 | 1.15 | 0/0.05 | 0.87 | 645 | 561 |
| 32 | 16 | 3 | 77 534 | 96 | 80 | 1.20 | 1909/1926 | 0.83 | 41 | 6545 |
| 40 | 8 | 3 | 498 103 | 60 | 80 | 0.75 | 27/32.1 | **1.54** | 269 | 969 |
| 40 | 16 | 2 | 860 359 | 80 | 80 | 1.00 | 0/0.2 | 1.01 | 561 | 561 |
| 40 | 16 | 3 | 73 746 | 82 | 80 | 1.02 | 2018/2025 | 0.98 | 33 | 6545 |

**Evaluation.**
- **P1 failed in 2 cells, held in 10.**
- **The kill signal was triggered** at (32, 8, 3) and stays recorded as triggered.
- **P2 held** in the 6 cells within its range.
- **P3, under the frozen statistic T(1 + p̄)/hits, failed in 4 cells and held in 8.**
  - The draft had scored the runner's T·max(1, p̄)/hits instead; codex caught this.
  - In the four rows with more than one s, the best measured s equals the model's arg-min (s = 2).
- **P4 held** in 12/12. **P5 held**: 450/450 lattices with 0 mismatches.

**Post-hoc diagnosis (disclosed).**
- *(32, 8, 3).* The 395 hits are 4 distinct members, one of them found 205 times.
  - Weighting the members by the sampler predicts 396.7 detections (codex), so the count is explained *given the
    catalogue*.
  - The catalogue itself is 4 vs 0.84 (Poisson p ≈ 0.011 in an independent-residue model, given as an
    illustration). It is reported, not explained.
  - The 4 members share at most 4 roots pairwise. No pair of symbols has a small rational ratio.
- *(24, 8, 1).* The space has 31 members against 28.0 expected, which is normal.
  - The sampler-weighted expectation is 70.7, but only 43 were observed.
  - A direct-membership replay reproduces the 43.
  - Five fresh seeds give 61, 69, 78, 82 and 83.
  - So this is a roughly 3σ low outlier (binomial tail ≈ 2.6·10⁻⁴; codex). It is **unexplained**, and points
    in the negative direction.

**Reading.**
- Where B = 16 permits independent sampling, the model predicts the measurements within Poisson error.
- The best measured costs per relation (s = 2, B = 16) are:
  - 120 at n = 24, against 2.6·10⁴ for brute force;
  - 9.9·10³ at n = 32, against 6.3·10⁶;
  - 1.08·10⁴ at n = 40, against 6.0·10⁶. This row uses V = 65535, so it is not comparable with the others.
- Three bit sizes cannot establish N^{1−o(1)}; that claim rests on the model.

### 5.2 E2: multiplier lattices

21 cells × 100 lattices; D = N^d in all 2100.

**Evaluation.**
- **P6 held**: median ‖b₁‖/GH from 1.01 to 1.15.
- **P7 held**: residual exponents 0.019–0.072 below d/n.
- **P8, first part, held.** There were 0 useful vectors at d = 2. There was one at (32, d1, j2), a Dixon-type
  vector.
- **P8, second part, FAILED in 5 of the 10 eligible cells**, all with fewer roots than the baseline:

  | (n, d, j) | observed | baseline | factor |
  |---|---|---|---|
  | (32,1,2) | 1 | 12.8 | 12.8 |
  | (32,1,4) | 22 | 68.8 | 3.13 |
  | (48,1,4) | 4 | 16.9 | 4.24 |
  | (48,2,10) | 10 | 53.8 | 5.38 |
  | (64,2,10) | 2 | 11.1 | 5.53 |

- *Post-hoc diagnosis.*
  - The registered baseline resampled each coefficient within its own magnitude, which overweights small
    coefficients.
  - Against a uniform-height baseline the rates agree within about 1.5× wherever the expected count is at least
    10. The exception is (32,1,2): 1 vs 5.3, a small count.
  - Agreement on this one statistic does not show that the multipliers are random.
- **The kill signal was not triggered.**

### 5.3 E3: inheritance

Here N is prime and F̃ is a monic factor of P − r mod N (diagnostic harness).

| (k, M, R) | tested | registered P9(e) violations | (H) holds | corrected IH violations |
|---|---|---|---|---|
| (1,4,8) | 23 | 0 | 23 | 0 |
| (1,5,10) | 4 | 0 | 4 | 0 |
| (1,6,11) | 8 | 0 | 8 | 0 |
| (2,4,8) | 61 | **15** | 61 | 0 |
| (2,5,8) | 93 | **10** | 91 | 0 |
| (2,6,6) | 6 | **6** | 6 | 0 |
| (3,5,8) | 21 | 0 | 0 | — |
| (3,6,6) | 30 | 0 | 8 | 0 |

**Evaluation.**
- **P9(a)–(d) held** (246/246). Clause (b) is tautological.
- **P9(e) failed in 31 pairs.** The harness tested the data-dependent form ⌈(n_supp − (k−1))/deg D⌉.
- **Post hoc:**
  - with denominator 2 deg D the data-dependent form has 0 violations;
  - the theorem form of Prop IH(b) and (c) has 0 violations in the 201 pairs that satisfy (H);
  - exhaustive tests under (H) also pass;
  - codex's 12 pairs without (H) are consistent with the corrected statement.
- **P10 held** (35/35).

---

## 6. Verdict and residual open statement (Task 4)

**Settled.**
- **Reformulation.** Evaluated O1 is exactly Π. Non-exact O1 is contained in Π under the scope given in Prop EQ.
- **Why the exact methods fail (Prop BRK).** Lemma R needs differences below 2^{M/2}, while non-exact differences
  have height at least N.
- **Counting.**
  - There is no scarcity obstruction for natural catalogues: Prop AB unconditionally (doubly split), Prop AB′
    under GRH (one-sided).
  - Rigorous small-box counting saves at most n^{O(1)}, which is vacuous in the target window (Prop CS).
  - The determinant-method theorems do not apply directly (Prop BP-NO).
- **Inheritance.** It is computably relation-valued under (H) (Prop IH). It is absent for N-independent fibres
  (Prop NF′). Independence across a family is not proved.
- **Coherence.** Linearly coherent families reduce to O2 (Prop COH). The seed pencils obey the alignment
  dichotomy of Remark PEN.
- **Lattices.** Exhaustive linearised enumeration costs N^{1−o(1)} per relation (Prop LIN, model). Multiplier
  lattices give generic multipliers.

**Residual open statement Π_β.** Fix 0 < β < 1/3 and ℓ = x^{2β}(log x)^{1−2β}. Is there a uniform classical
algorithm with the following properties?
- **Input and cost.** Given only N, it runs in time L[β]^{1+o(1)}, including the cost of building any
  presentation.
- **Output.** L[β]^{1+o(1)} congruences ∏_{a∈A}x_a ≡ v (mod N).
- **Admissibility.**
  - Each x_a is a unit from a factor base of at most L[β] elements: integers up to L[β], or m − a at a computed m
    with |a| ≤ L[β].
  - Multipliers are absent or split over the factor base.
  - Either 0 < |v| ≤ e^{O(ℓ)}, or v is itself split over the factor base.
  - u ≠ v as integers, and no output is an exact integer identity.
- **Rank.** Written over the primes, with known identities between pieces identified and small sides factored,
  the outputs have F₂-rank at least (1 − o(1))·#columns.

Π_β is only a listing problem. A factoring algorithm also needs fruitful square dependencies and the costs of
smoothness testing and linear algebra (heurnfs H; Croot–Granville–Pemantle–Tetali 2012).

**Relation to smoothenum.**
- One-sided Π_β is an algebraic-side lister in the sense of T6.
- Its doubly split form is related to T7.
- T7's infimum 0 is not polynomial time by itself. T8 needs further hypotheses: a polylogarithmic factor base,
  prefix-output, joint smoothness and fruitful extraction.
- The draft's "would give heuristic polynomial time" is **withdrawn**.

**Probability (assessment, subjective).** Below 1% (ours); about 5% (codex, blind).

---

## 7. References

**Lattice points.**
- Minkowski, *Geometrie der Zahlen* (1896); Cassels, *An Introduction to the Geometry of Numbers*.
- Henk, Rend. Circ. Mat. Palermo (2) Suppl. 70 (2002) 377–384, Thm 1.5.
- Barroero–Widmer, IMRN 2014(18) 4932–4957, Thm 1.3.
- Fincke–Pohst, Math. Comp. 44 (1985) 463–471.

**Character sums and smooth numbers in progressions.**
- Burgess, Mathematika 4 (1957) 106–112; Proc. LMS (3) 12 (1962) 179–192; Proc. LMS (3) 13 (1963) 524–536.
- Pólya–Vinogradov.
- Granville–Soundararajan, JAMS 14 (2001) 365–397.
- Montgomery–Vaughan, *Multiplicative Number Theory I*, Ch. 13.
- Iwaniec–Kowalski, *Analytic Number Theory*, Ch. 5.
- Soundararajan, CRM Proc. LN 46 (2008) 115–128, Thm 1.
- Harper, JNT 132 (2012) 182–199.

**Short-interval multiplicative congruences.**
- Cilleruelo–Garaev, GAFA 21 (2011) (arXiv 1007.1526).
- Bourgain–Garaev–Konyagin–Shparlinski, Proc. Steklov Inst. 280 (2013) (arXiv 1203.0017); J. Anal. Math. 124 (2014)
  (arXiv 1210.6429).
- Chang et al., arXiv 1111.1543.

**Determinant method.**
- Bombieri–Pila, Duke 59 (1989) 337–357, Thm 5.
- Heath-Brown, Annals 155 (2002) 553–595.
- Castryck–Cluckers–Dittmann–Nguyen, "The dimension growth conjecture, polynomial in the degree and without
  logarithmic factors", ANT 14 (2020).
- Binyamini–Cluckers–Novikov, arXiv 2302.04209, Thm 2, eq. (27).

**Lattices and factoring.**
- Lenstra–Lenstra–Lovász, Math. Ann. 261 (1982).
- Vallée, Math. Comp. 56 (1991) 823–849.
- Schnorr, EUROCRYPT '91, LNCS 547, 281–293; ePrint 2021/232.
- Adleman, "Factoring and lattice reduction" (unpublished, 1995).
- Ducas 2021.
- Aboumrad–Widdows–Kaushik, arXiv 2308.07804.
- Coppersmith, J. Cryptology 10 (1997) 233–260, Thm 1.
- Boneh, STOC 2000, 265–272; JCSS 64 (2002) 768–784.
- Dixon, Math. Comp. 36 (1981) 255–260.
- Morrison–Brillhart, Math. Comp. 29 (1975).
- Croot–Granville–Pemantle–Tetali, Annals 175 (2012) 1507–1550.

**Generic groups.**
- Wagner, CRYPTO 2002, LNCS 2442, 288–304, Thm 2 (W. Dai) and Cor. 1.
- Nechaev, Math. Notes 55 (1994).
- Shoup, EUROCRYPT '97.

## 8. Reviews

**Codex blind round (`nonexact_1.out`).** Same verdict, probability about 5%. Its contributions are credited
throughout:
- the Smith determinant;
- Henk and Barroero–Widmer;
- conditional abundance;
- the query bound and Wagner;
- the existence half of IH;
- Zariski density;
- the View A degeneracy;
- meet in the middle and Coppersmith;
- the containment (eq. 51);
- useful rows;
- the residual statement.

**Codex adversarial round (`nonexact_2.out`): MAJOR REVISION. All 19 points applied.**
1. EQ multiplier.
2. Π_β admissibility.
3. BRK certificates and the k = 1 wording.
4. AB headline.
5. LAT determinant.
6. VA gcd weights.
7. CS quantifiers; Cilleruelo–Garaev.
8. BP-NO downgraded.
9. IH hypothesis (H).
10. IH existence; independence wording.
11. PEN rank-1 claim withdrawn; dichotomy instead.
12. NF′/COH scope.
13. LIN scoped to exhaustive enumeration.
14. P3 rescored: 4 failures.
15. P8: 5 failures.
16. E1 counts and n = 40 figures.
17. (24, 8, 1) unexplained; model-based p-values.
18. Corrected IH helper and tests.
19. Verdict separates reformulation from exclusion; T8 withdrawn.

Codex also re-scored P1–P10 independently and reproduced the exhaustive counts.

**Opus proof critic: ACCEPT WITH FIXES.** All applied:
- EQ multiplier and containment scope;
- the AB overclaim;
- IH (H), with counterexamples;
- E1/E2 numbers, including the missing P8 cell;
- VA coprimality;
- LAT determinant, upper bound only;
- CS vacuity and exponent;
- BRK details;
- PEN wording;
- LIN wording;
- BP-NO box degree;
- P9(b) tautological;
- stale code;
- labels.

The critic verified NF′, COH, the LIN inequalities (numerically), the PEN algebra, the IH support count and the
2 deg D denominator.

**Opus literature critic: OK WITH FIXES.** All applied:
- CS follows from GRH;
- BCN c(n);
- BGKS and Cilleruelo–Garaev;
- Vallée and Schnorr–Adleman;
- Boneh;
- Burgess and GRH wording;
- Harper;
- Henk and Barroero–Widmer wording;
- Wagner–Dai, Nechaev and Shoup;
- tool references;
- unused references dropped.

It found no overclaimed novelty.

No confirmation round was run on this revision.

## 9. Reproduction

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_nonexact.py
PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment E1   # ~3 min
PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment E2   # ~6 s
PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment E3   # ~3 s
PYTHONPATH=. .venv/bin/python -m src.nonexact.posthoc             # replay: distinct members (~2 min)
PYTHONPATH=. .venv/bin/python -m src.nonexact.posthoc exhaustive  # exhaustive member counts (1.5 GB peak)
# e2_baseline2, e3_reasons, rescore, e3_theorem_check in src/nonexact/posthoc.py (post hoc)
```
