# ptebarrier: a rigidity barrier for route 1's PTE escape

Agent `ptebarrier`, phase 16 (route 1, barrier side), 2026-10-20.

**Files.**
- Code: `src/ptebarrier/`.
  - `rigidity.py`: thresholds of Lemma R. A step-by-step certificate of its proof. The Theorem W invariant. The Prop S small-value set.
  - `enum.py`: exhaustive enumeration of split polynomials.
  - `run_experiment.py`: the runs.
- Tests: `tests/test_ptebarrier.py` (10).
- Pre-registration: `notes/ptebarrier_prereg.md` (commit d6d8073, made before any run).
- Config: `experiments/configs/ptebarrier.json`.
- Raw results: `experiments/raw/ptebarrier_E{1,1b,3}.jsonl` (pre-registered) and `ptebarrier_E1c_posthoc.jsonl` (post hoc, disclosed).
- Codex (GPT-6.1-sol): blind round `codex/ptebarrier_1.{md,out}`; adversarial round `codex/ptebarrier_2.{md,out}` (§8).

**Labels.**
- **PROVED**: proof given here.
- **KNOWN**: published, with a citation.
- **HEURISTIC**: the assumption is named.
- **measured**.
- **open**.
- **assessment**.

Novelty is never claimed; "no prior statement found by us" is the most we say.

**Notation.**
- x = log N; L[β] = exp(Θ(x^β (log x)^{1−β})).
- A multiset A of M integers has the split polynomial P_A = ∏_{a∈A}(X − a), with e_i = e_i(A).
- "Top part" means the coefficients of degree ≥ k, i.e. e₁..e_{M−k}.
- s = M − k.
- ‖·‖ is the max-coefficient height.

---

## 0. Verdict

**Plain language.** The "dense split family" escape asks for *one* N-dependent congruence that makes *many* split (smooth-by-construction) polynomials reduce to small ones. That is how the finite-field algorithms get q³ relations from one presentation. Over the integers this amplification is impossible.

The reason is archimedean. Take two split integer polynomials of large degree M with bounded roots. If they differ by a polynomial that is small and of low degree, they are equal. With a small multiplier allowed, they are equal up to that multiplier.

So every systematic compressed family on one presentation yields **the same relation for all of its members**. In the integer-evaluation view its split sides are not even dense: they use at most deg F distinct pieces. A family of L[α] relations therefore needs L[α] separately N-conditioned presentations. That is ordinary relation search, and the compression contributes nothing.

**Proved.**

- **Lemma R (rigidity) [PROVED, elementary, uniform in every parameter].** Let P, P' be monic of degree M, split over Z with roots in [−B, B]. Let Q, Q' ≠ 0 have degree ≤ j, and let G = QP − Q'P' have degree < K and height ≤ h. Put L = ⌊log₂(2^{K−1}√K·h·B^{K−1})⌋. If M ≥ K + 2L + 4j + 3, then QP = Q'P'.
  - For j = 0 this is an exponential form of ptecons' factorial bound CW, which is sharper there.
  - The multiplier version (j ≥ 1) is what covers BGJT-shaped families.
- **Theorem W (one presentation, one relation) [PROVED].** Let F ∈ Z[X] be arbitrary: any degree, any height, possibly depending on N. Every family of exact identities h_A·P_A = F·Q_A + r_A has the following property, provided
  - the multipliers h_A, Q_A have degree ≤ j,
  - the small sides r_A have degree < k and height ≤ H,
  - the regime of Lemma R holds, with K = k + j and h = 2(j+1)·H·max‖Q_A‖.

  Then the rational function h_A·P_A/r_A is the same for all members, so every member yields *the same* relation.
- **Corollary P (the PTE escape exactly as posed) [PROVED].** Let the members share e₁..e_{M−k}, and let the presentation have degree d ≥ k. Then the quotient is common, and an exact compressed family has **at most one member** once M ≥ k + 2⌊log₂(2^k√k·H·B^{k−1})⌋ + 3. ptecons' CW gives a sharper factorial threshold.
- **Prop S (integer evaluation: no density) [PROVED].** Suppose N | F(m) ≠ 0 and every root lies within Y of m. Every exact identity hP = FQ + r has the split side's roots in S_F ∪ Z(r), where
  - |S_F| ≤ deg F, and S_F is fixed by F;
  - h and Q are arbitrary;
  - r has degree ≤ ρ and height ≤ H;
  - the condition is (ρ+1)·H·(|m|+Y)^ρ·(deg F + 1)·Y^{deg F} < N.

  After cancelling common factors, every relation through F uses only the ≤ deg F pieces m − c, c ∈ S_F.
- **Prop NX (non-exact compressions do not share) [PROVED].** In the PTE setting with d ≥ k, distinct members have distinct N-defects. Every member therefore pays its own congruence mod N; nothing is inherited from the rest of the family.
- **Determinant method, k = 1 [PROVED from KNOWN].** A PTE fiber with one free coefficient and members with ≥ 3 distinct roots has at most C·B^{1/2}(log B)^κ + 2 members. This holds uniformly in M, for *arbitrary* constant differences (including non-exact ones). Here C and κ are absolute, from the uniform Bombieri–Pila bound of Binyamini–Cluckers–Novikov. So one ideal-type fiber cannot fill a factor base of size ≍ B.

**What the counting tools give.**
- *Vinogradov (BDG/Wooley).* It bounds the second moment of fiber sizes, i.e. PTE collisions. Its best pointwise consequence is max fiber ≤ J_{M/2,s}(B) ≪_{M,s,ε} B^{M/2+ε}, whose constants are not uniform in M and s.
  - It is weaker than the elementary injection bound C(2B+k, k) whenever k < M/2, i.e. in the entire compression range.
  - **Vinogradov is the wrong tool**: the barrier lives in the archimedean height of the remainder, not in fiber cardinality.
- *Dimension growth.* Dimension-growth bounds (Salberger; Castryck–Cluckers–Dittmann–Nguyen) have constants depending on the ambient dimension M and polynomial in the degree s!. They are useless here except for k = 1.

**What it rules out (§6).** Take B = L[β] with β < 1/3 and M ≥ (1 − o(1))x/log B. Suppose the small-side degree, multiplier degree and log-heights satisfy (k + j)·log B + log(H·H_Q) = o(x^{1−β}). Every BGJT-type family qualifies: k = 4, j ≤ 3, small heights. Then route 1's open escape *in its systematic form* (exact identities over Z[X] on one presentation) provides **no amplification**: ≤ 1 relation per presentation. It therefore cannot by itself give α < 1/3, nor any improvement over the method used to produce the presentations.

**Not closed (open, stated precisely in §6.3).**
- (O1) *Non-exact or evaluated compressions.* The small side is small only modulo N, or only in value at m. This is generic relation search plus a common-offset structure: ptecons' NC primitive, and smoothenum's doubly-smooth listing.
- (O2) *Linear systems of two or more unrelated presentations* (pencils F₀ + νF₁). The missing lemma is rigidity for small three-term combinations of split polynomials.
- (O3) *Parameter corners.* k·log B ≳ x^{1−β}, which needs a tiny m and β ≥ 1/4. Or large multiplier heights.
- (O4) *Other root sets.* Non-integer roots in non-discrete rings (real quadratic) and multivariate presentations. Gaussian or Eisenstein roots go through with changed constants.

**Experiments (pre-registered).**
- E1/E1b check Lemma R exhaustively: 318 + 147 cells, **0 violations**. E1c (post hoc, larger B) adds 30 cells, also 0 violations.
- *Disclosed:* these are weak tests. Only 2 pre-registered cells (and 8 post-hoc cells) are non-vacuous, i.e. have a nonzero h_R and contain a pair at all. For B ≤ 5 no pairs with an equal top part exist once M − K ≥ 5. P2 holds only vacuously. The lemma rests on its proof. The P5 expectation (structured fibers persisting above √(2M)) **failed**: there were none at all for B ≤ 4.
- E3: the Dickson presentation has exactly one split member in a 401² window.

**Assessment.** Probability < 1% that route 1's PTE escape yields sub-L[1/3]. The surviving open cases O1 and O2 are not compressions in the BGJT sense: O1 is the original relation-finding problem, and O2 has no candidate.

---

## 1. Model and cost targets

### 1.1 Two views of a presentation

A **presentation** is a polynomial F ∈ Z[X] \ {0} together with a ring homomorphism φ: Z[X] → Z/N with φ(F) = 0.
- **View B (residue).** φ(X) = μ ∈ Z/N. The linear factor base is {μ − a : |a| ≤ B} together with small primes. This is the literal analogue of F_{q²}[X] → F_{q^{2n}}.
- **View A (integer evaluation).** φ(X) = m ∈ Z, so N | F(m). The pieces m − a are actual integers, factored over the primes ≤ Y. This is systeq §6's (u, v) formulation and the ptesplit/ptecons model: u = P(m), v = r(m).

**A compressed split family on F** is a set of members (A, h_A, Q_A, r_A) where:
- A is a multiset of M integers. All roots lie in [−B, B], or within Y of m in View A.
- h_A, Q_A ∈ Z[X] \ {0} are multipliers of degree ≤ j.
- r_A ∈ Z[X] \ {0} is the small side, of degree < k and height ≤ H.

The member's relation is φ(h_A)·φ(P_A) = φ(r_A) in Z/N.
- The member is **exact** if h_A·P_A = F·Q_A + r_A holds in Z[X]. This is the "systematic" case: the analogue of the BGJT identity h₁P = (aX + b)(h₁X^q − h₀) + r.
- Otherwise it is **non-exact**: the identity holds only modulo (N, F).

The PTE escape as posed in systeq §6 is the case where the members share e₁..e_{M−k}. With a presentation of degree d ≥ k, the quotient Q_A is then common (Corollary P).

**Lemma EX (systeq §5.5)** says that a congruence modulo (N, F̃) whose reduced difference has coefficients below N is exact. So "exact" is the regime of small reductions, and "non-exact" members have reduced defects of height ≥ N.

### 1.2 Cost targets

These are HEURISTIC, in the random-residual model of ptesplit_1 / ptecons §1.
- *Size of the split side.* Pieces ≤ Y = L[β] and a non-exact integer congruence force M·log Y ≥ (1 − o(1))x.
- *Size of the small side.* For an L[β] algorithm the small side's log-size must be O(x^{2β}(log x)^{1−2β}). So log|v| = o(x^{2/3}) alone is not enough for a fixed β < 1/3.
- *Family size.* The family must supply ≳ L[β] relations at amortised cost L^{o(1)} each.

Theorem W makes the family-size target moot: one presentation supplies one relation, whatever is needed.

## 2. Heuristic count (Step 1) [HEURISTIC + measured]

With q = 2B + 1, there are C(q+M−1, M) ~ q^M/M! multisets (for M² = o(B)). The vector (p₁, …, p_s) has a range of size about ∏_{i≤s} M·B^i = M^s·B^{s(s+1)/2}. Under a local-limit model in which the moment-body density and the singular series are absorbed into a factor 𝒜_{M,s}(c) (codex ptebarrier_1 eq. (1)), the expected fiber is

  λ(c) ≈ (q^M/M!)·B^{−s(s+1)/2}·𝒜_{M,s}(c).

λ drops below 1 at s(s+1)/2 ≈ M − log(M!)/log B, i.e. s ≈ √(2M) when log M ≪ log B. So dense fibers generically need s ≲ √(2M) prescribed functions, i.e. k ≥ M − O(√M) free coefficients: no compression.

A compressed family also needs the k low coefficients (which have size up to C(M,i)·B^{M−i}) to land in a window of height H. That multiplies λ by a factor of order ∏_{i<k} min(1, H/(C(M,i)B^{M−i})), which is astronomically small. Lemma R turns this heuristic into a theorem for the exact case: the true count is 1.

**Measured (E2, B ≤ 4, all M of E1).**
- No nontrivial fiber has s ≥ √(2M) + 2.
- Max fibers sit far below the injection bound C(q+k−1, k). For example, (B, M, k) = (4, 6, 4) has max fiber 10 against a bound of 495.

The pre-registered P5 expectation that structured PTE fibers would persist above √(2M) at small M **failed**: at B ≤ 4 there are none. Known ideal PTE solutions need entries in the hundreds.

## 3. Vinogradov's mean value theorem (Step 2)

Put R(c) = #{a ∈ [1, B]^M : p_i(a) = c_i, i ≤ s} (ordered tuples) and f(α) = Σ_{a≤B} e(α₁a + … + α_s a^s).

- **KNOWN.** BDG 2016 (Ann. Math. 184), with Wooley's efficient congruencing (cubic 2016; all degrees, Proc. LMS 2019):

  J_{t,s}(B) = ∫|f|^{2t} ≪_{t,s,ε} B^{t+ε} + B^{2t − s(s+1)/2 + ε}.
- **What it bounds [PROVED].** Σ_c R(c)² = J_{M,s}(B), the number of ordered PTE collisions. Hence:
  - #{c : R(c) ≥ T} ≤ min(B^M/T, J_{M,s}/T²): *most* fibers are small;
  - #{c : R(c) > 0} ≥ B^{2M}/J_{M,s}.
- **Pointwise [PROVED, codex ptebarrier_1].** R(c) = ∫ f^M e(−c·α) ≤ ∫|f|^M = J_{M/2,s}(B) for even M (Cauchy–Schwarz for odd M). So max fiber ≪_{M,s,ε} B^{M/2+ε} + B^{M − s(s+1)/2+ε}.
  - The suggested "≤ J^{1/2}" (Cauchy–Schwarz on the second moment) is weaker: J_{M,s}^{1/2} ≥ √(M!)·B^{M/2}.
- **Comparison [PROVED].** The injection bound "keep the k largest entries, and Newton recovers the rest" gives max fiber ≤ C(q+k−1, k) ≤ q^k, with no constants. Against B^{M/2}, Vinogradov wins only when k > M/2, i.e. outside the compression range.
- **Uniformity [PROVED + KNOWN].** The implied constants depend on (M, s, ε). The diagonal alone gives J_{M,s} ≥ M!·(B)_M, and log M! ≫ x^β in our range, so constants cannot be suppressed.
  - Explicit forms exist but carry exp(O(s³ log s))-size constants: Wooley 1993 (Vinogradov II); Ford's explicit Theorem 3 (s ≥ 200); Steiner's efficient-boxing bounds (codex ptebarrier_1 §3).
- **Conclusion.** Vinogradov bounds the *number of collisions* and shows that dense fibers are *rare*. The escape needs only *one* dense, N-compatible fiber, so this is the wrong tool. The tool that works (§5) uses the archimedean size of the remainder, which Vinogradov does not see.

## 4. Determinant method / dimension growth (Step 3)

- **Geometry [PROVED].** V_c = {e₁ = c₁, …, e_s = c_s} ⊂ A^M is a complete intersection of dimension k = M − s and degree ≤ s!.
- **Uniform bounds [KNOWN].**
  - Salberger's dimension growth and Castryck–Cluckers–Dittmann–Nguyen (Algebra Number Theory 14 (2020)) give N(X, B) ≤ c(n)·D^{e(n)}·B^{dim X} (projective, D ≥ 5). The constants depend on the ambient dimension n = M.
  - With D = s! ≈ e^{M log M} = L[1 − β]^{Θ(1)}, the degree factor alone exceeds every L[β] budget. **Not useful for k ≥ 2** [PROVED, by comparison].
- **k = 1 [PROVED from KNOWN].** Write Q₀ for the common top part, so the members are Q₀ − t splitting over Z.
  - If Q₀ − t has M ≥ 3 distinct integer roots in [−B, B], every ordered pair (x, y) of distinct roots lies on the plane curve C: (Q₀(x) − Q₀(y))/(x − y) = 0. C has degree M − 1, and distinct t give disjoint pairs.
  - *Linear components.* Q₀(λx + ν) ≡ Q₀(x) forces λ^M = 1, so λ = ±1. λ = 1 forces ν = 0, which does not lie on C. So C has at most one line, y = ν − x, carrying ≤ M pairs per t.
  - *Components not defined over Q* carry ≤ e² rational points each (Bézout with a conjugate).
  - *Components over Q of degree e ≥ 2* carry ≤ c·e²·B^{1/e}(log B)^κ points of height ≤ B. This is the Binyamini–Cluckers–Novikov uniform affine bound; codex cites it as c(n)D²B^{1/D}(log B)^κ with n = 2 fixed. UNVERIFIED by us beyond codex's citation; see the literature critic, §8.
  - Summing: |T|·M(M−2) ≤ c(M−1)²B^{1/2}(log B)^κ + (M−1)². So **|T| ≤ 2c·B^{1/2}(log B)^κ + 2** for M ≥ 4, uniformly in M.
  - This holds for *any* constant differences, exact or not.
  - Members with ≤ 2 distinct roots are of the form (X − a)^p(X − b)^{M−p}. With e₁ and e₂ fixed there are O(M) of them.
- **Consequence.** A single ideal-type fiber supplies ≤ B^{1/2+o(1)} members, below the ≍ B needed to fill a factor base. Without small remainders Lemma R does not apply; this is the k = 1 complement.

## 5. Rigidity (Steps 4–5)

### 5.1 Lemma R [PROVED]

**Statement.**
- Let P, P' ∈ Z[X] be monic of degree M, split over Z, with all roots in [−B, B] (B ≥ 1).
- Let Q, Q' ∈ Z[X] \ {0} have degree ≤ j.
- Put G := QP − Q'P', and suppose deg G ≤ K − 1 and ‖G‖ ≤ h.
- Let Δ := 2^{K−1}·√K·h·B^{K−1} and L := ⌊log₂ Δ⌋.

Then:
- if M ≥ K + 2L + 4j + 3, we have QP = Q'P';
- for j = 0 (Q = Q' = 1) the condition M ≥ K + 2L + 3 suffices.

**Proof.** Suppose G ≠ 0.

1. **Coprime parts.** Write P = gU and P' = gU', where g is the monic gcd (the common sub-multiset) and U, U' are coprime and split, of degree r = M − deg g. Then G = g·G̃ with G̃ := QU − Q'U' ≠ 0, and deg g ≤ deg G ≤ K − 1, so r ≥ M − K + 1.
2. **Degrees.**
   - If deg(QU) ≠ deg(Q'U'), or the leading coefficients differ, then deg G̃ ≥ r. That gives M = r + deg g ≤ deg G ≤ K − 1, which is false.
   - So deg Q = deg Q' =: j' and lc(Q) = lc(Q').
3. **Height.**
   - Mahler measure is multiplicative and M(g) ≥ 1, so M(G̃) ≤ M(G) ≤ ‖G‖₂ ≤ √K·h.
   - Every coefficient of a degree-n polynomial f satisfies |c_i| ≤ C(n, i)·M(f), so ‖G̃‖₁ ≤ 2^{K−1}√K·h.
   - Hence |G̃(a)| ≤ ‖G̃‖₁·B^{K−1} = Δ for every integer |a| ≤ B.
4. **Local smallness.**
   - Call a root a of U *good* if Q'(a) ≠ 0. At a good root, G̃(a) = −Q'(a)U'(a), with Q'(a) a nonzero integer and U'(a) ≠ 0 (coprimality). So 1 ≤ |U'(a)| = ∏_{a'}|a − a'| ≤ Δ.
   - Every factor is a nonzero integer, and each factor ≥ 2 doubles the product. So **at most L roots of U' (with multiplicity) lie at distance ≥ 2 from a**, and supp U' ⊆ {a − 1, a + 1} ∪ E with |E| ≤ L.
   - The same holds with the roles of U and U' exchanged (good roots of U' are those with Q(a') ≠ 0).
5. **Moments.**
   - QU/lc and Q'U'/lc are monic of degree r + j' and differ by a polynomial of degree ≤ K − 1 − deg g. So the elementary symmetric functions e₁..e_J of their complex root multisets agree, with J = r + j' − K + deg g ≥ M − K.
   - By Newton's identities (characteristic 0), the power sums p₁..p_J agree.
   - Let μ = Σ_{roots of QU} δ − Σ_{roots of Q'U'} δ (with multiplicity). Then μ ≠ 0, and Σ_z μ(z)z^i = 0 for 0 ≤ i ≤ J.
   - If |supp μ| ≤ J + 1, the Vandermonde system would force μ = 0. Hence **|supp μ| ≥ J + 2 ≥ M − K + 2.**
6. **Support count.** supp μ ⊆ supp U ∪ supp U' ∪ Z_C(Q) ∪ Z_C(Q'); the last two sets have ≤ j' elements each.
   - *Case 1: U and U' both have good roots.* Then |supp U|, |supp U'| ≤ 2 + L, so |supp μ| ≤ 4 + 2L + 2j.
   - *Case 2: U has no good root and U' does.* Then supp U ⊆ Z(Q'), with ≤ j values.
     - Two good roots a'₁, a'₂ of U' with |a'₁ − a'₂| ≥ 3 would put ≥ r − L roots of U into each of the disjoint sets {a'_i ± 1}. That gives r ≤ 2L, contradicting r ≥ M − K + 1 ≥ 2L + 4.
     - So the good roots of U' take ≤ 3 values, and the bad roots of U' (roots of Q) take ≤ j values.
     - Hence |supp μ| ≤ 4j + 3. The symmetric case is the same.
   - *Case 3: neither has a good root.* Then |supp μ| ≤ 4j.
7. **Contradiction.** The hypothesis gives M − K + 2 ≥ 2L + 4j + 5, which exceeds every case bound. For j = 0 only Case 1 occurs, and M ≥ K + 2L + 3 suffices. ∎

**Remarks.**
- **The reason is archimedean.** In F_q[X] there is no size, every nonzero constant is a unit, and the BGJT family exists. Over Z, two split polynomials cannot be close in their low coefficients unless they are equal. The step-4 argument (a product of nonzero integer gaps is small only if almost all the gaps are 1) uses only discreteness of Z in R. It therefore extends verbatim to roots in Z[i] or Z[ω]: replace {a ± 1} by the ≤ 8 lattice points at distance < 2. It does not extend to non-discrete rings such as Z[√2].
- **Comparison with ptecons CW (§3.4 there; codex ptecons_1).** For j = 0, CW proves ‖P − P'‖ ≥ K_s/(t·2^{t−1}·A^{t−1}) with K_s = ⌊s/2⌋!⌈s/2⌉!, s = ⌊(M−t)/2⌋ + 1.
  - Its proof replaces our steps 4–6 by a positivity argument: each cofactor has ≥ s distinct roots, and distinct gaps multiply to a factorial.
  - CW is **sharper** (factorial against our exponential 2^{(M−K−3)/2}·(2B)^{−(K−1)}).
  - Lemma R's contribution is the multiplier case j ≥ 1, where W = ∏(X − a)² positivity fails on the complex roots of Q and Q' but the Vandermonde count still works.
  - Both thresholds have the same shape in the L-range: M ≳ K·log B.
- **Proof certificate.** `rigidity.lemma_R_certificate` recomputes steps 1–6 on any pair and checks each inequality. `tests/test_ptebarrier.py::test_certificate_random_pairs` runs it on 800 random pairs with j ∈ {0, 1}.

### 5.2 Theorem W: one presentation, one relation [PROVED]

**Statement.** Let F ∈ Z[X] \ {0} be arbitrary. Consider members (A, h_A, Q_A, r_A) with:
- h_A·P_A = F·Q_A + r_A in Z[X];
- deg h_A ≤ j₁ and deg Q_A ≤ j₂, with Q_A ≠ 0 and ‖Q_A‖ ≤ H_Q;
- r_A ≠ 0, deg r_A < k and ‖r_A‖ ≤ H;
- roots in [−B, B].

Put K = k + j₂, j = j₁ + j₂, h = 2·min(k, j₂+1)·H_Q·H. If M ≥ K + 2⌊log₂(2^{K−1}√K·h·B^{K−1})⌋ + 4j + 3, then:
- the rational function Ψ_A := h_A·P_A/r_A ∈ Q(X)^× is the same for all members;
- for every homomorphism φ with φ(F) = 0, all members define the same relation φ(h_A)φ(P_A)φ(r_A)^{−1} = 1. The relation is the same as an exponent vector over the irreducible factors in Z[X] (View B) and over the primes after evaluation at m (View A).

The family's relation lattice has rank ≤ 1.

**Proof.**
1. For two members, Q_{A'}·(h_A·P_A) − Q_A·(h_{A'}·P_{A'}) = Q_{A'}·r_A − Q_A·r_{A'} =: G. Here deg G < k + j₂ and ‖G‖ ≤ h.
2. The multipliers Q_{A'}h_A and Q_A h_{A'} are nonzero, of degree ≤ j.
3. By Lemma R, G = 0, so Q_{A'}r_A = Q_A r_{A'}, and therefore Ψ_A = Ψ_{A'}.
4. Z[X] is a UFD, so equal rational functions have equal exponent vectors. Evaluation at m gives equal rational numbers. ∎

**Notes.**
- **The heights of h_A and of F never enter.** F may have N-sized coefficients.
- **BGJT shape.** The BGJT identity is h₁P = (aX + b)(h₁X^q − h₀) + [(aX + b)h₀ + (cX + d)h₁]: j₁ = 2, j₂ = 1, k = 4. In F_q[X] it gives q³ distinct Ψ's. Over Z, Theorem W allows one.
- **Forced structure.** By the proof, the members satisfy h_A P_A = Q_A·Φ₀/Φ₁ with Φ₀/Φ₁ = F + r/Q fixed. When deg F ≥ k they share the factor Φ₀ of degree ≥ M − j. They are one relation dressed in trivial identities. This is tested in `test_theorem_W_forced_structure_has_one_invariant`.

### 5.3 Corollary P: the PTE escape as posed [PROVED]

**Statement.** Let the members share e₁..e_{M−k}. Let F have degree d ≥ k, and let the members be exact: P_A = F·Q_A + r_A with deg r_A < d and ‖r_A‖ ≤ H. Then:
- Q_A is common to all members;
- P_A − P_{A'} = r_A − r_{A'} has degree < k and height ≤ 2H;
- there is at most one member once M ≥ k + 2⌊log₂(2^k√k·H·B^{k−1})⌋ + 3.

**Proof.**
- F(Q_A − Q_{A'}) = (P_A − P_{A'}) − (r_A − r_{A'}) has degree < d, so Q_A = Q_{A'}.
- Apply Lemma R with j = 0. ∎

**Two remarks.**
- The heights of F and of the common quotient Q never enter.
- The Dickson identity (X² − 4)(X² − 1)² = F + (X + 3), with F = X⁶ − 6X⁴ + 9X² − X − 7, is the k = 2 case at M = 6. That is out of regime, and E3 finds it to be the only split member with |a|, |b| ≤ 200.

### 5.4 Prop S: integer evaluation forbids density [PROVED]

**Statement.** Let m ∈ Z and F ∈ Z[X] of degree ≤ d, with F(m) ≠ 0 and N | F(m). Let Y, H ≥ 1 and ρ ≥ 0. Put B_m = |m| + Y and suppose

  (ρ+1)·H·B_m^ρ·(d+1)·Y^d < N.

Let S_F := {a ∈ Z : |a − m| ≤ Y, |F(a)| ≤ (ρ+1)H·B_m^ρ}. Then:
- (i) |S_F| ≤ d;
- (ii) for every identity hP = FQ + r in Z[X] with P monic and split, all roots within Y of m, h and Q arbitrary, and r ≠ 0 of degree ≤ ρ and height ≤ H: every root of P lies in S_F ∪ Z(r). If gcd(P, r) = 1 (a reduced relation), every root lies in S_F.

**Proof.**
- (i) Suppose a₀..a_d were d + 1 distinct points of S_F. Lagrange gives F(m) = Σ_i F(a_i)∏_{l≠i}(m − a_l)/(a_i − a_l). Each product is at most Y^d in absolute value, so |F(m)| ≤ (d+1)·(ρ+1)H·B_m^ρ·Y^d < N ≤ |F(m)|, a contradiction.
- (ii) At a root a, h(a)P(a) = 0, so F(a)Q(a) = −r(a). If r(a) ≠ 0, then F(a) and Q(a) are nonzero integers, so |F(a)| ≤ |r(a)| ≤ (ρ+1)H·B_m^ρ, and a ∈ S_F.
- Reduction: if (X − a) divides both P and r and F(a) ≠ 0, then (X − a) | Q, and the factor m − a cancels from both sides of the relation. ∎

**Consequences.**
- In View A, every reduced exact relation through F has u = ∏_{c∈S_F}(m − c)^{n_c}: at most d distinct pieces, the same d for every member.
  - A "dense" split side (≍ x/log Y distinct pieces) needs d log Y ≳ x − ρ log B_m − log H, i.e. a presentation of degree comparable to M.
  - At that degree the family on F is a PTE fiber, which Corollary P (or CW) collapses to one member.
- With PTE structure, two members give a signed measure on ≤ d + 2ρ points. So there is ≤ 1 member as soon as M − k ≥ d + 2ρ − 1, for *any* quotient heights. This closes the d < k case of the PTE escape in View A, which Corollary P leaves open.
- Relation to ptecons LEDGER: LEDGER is the norm identity ∏F(a) = Res(F, r). Prop S instead uses |F(m)| ≥ N through Lagrange interpolation. It needs no Mahler-measure bound on F, so it applies to N-dependent presentations of any height.

### 5.5 Non-exact members and single rings [PROVED]

**Prop NX (defects are injective).**
- *Setting.* PTE setting, F monic of degree d ≥ k, members with P_A ≡ r_A mod (N, F) and ‖r_A‖ ≤ H < N/2. Put D_A := red_F(P_A) − r_A ∈ N·Z[X]_{<d} (Lemma EX).
- *Injectivity.* For two members, D_A − D_{A'} = (P_A − P_{A'}) − (r_A − r_{A'}), since P_A − P_{A'} has degree < k ≤ d and is its own reduction. If D_A = D_{A'}, Lemma R (j = 0, height 2H) forces A = A'. So **A ↦ D_A is injective** on any family in the regime.
- *Common part.* The coefficients of D_A in degrees ≥ k are common. They equal those of red_F(A₀) shifted by the unique centred r-coefficient. So the family size is at most the number of attainable defect cells:

  |𝓕| ≤ ∏_{i<k} (⌊2W_i/N⌋ + 1),  with W_i = |red_F(A₀)_i| + C(M,i)·B^{M−i} + H.

- *Independence costs heights* (codex ptebarrier_1 (14)). r linearly independent defects satisfy N^r ≤ r^{r/2}·∏ ‖D_j‖.
- *Reading.* Exact members collapse (Corollary P). Every non-exact member needs its *own* congruence on its low coefficients, which are of size ≍ B^M. The PTE prescription of the top part constrains those residues mod N in no way that we can see.
  - Bounding how many fiber members hit the cells would need equidistribution of e_{s+1}..e_M mod N over a PTE fiber. That is **open**. The cell bound above is weak (≈ (B^M/N)^{O(k)}).

**Same remainder [PROVED, ptesplit_1 codex].** If P ≡ P' mod (N, F) with deg(P − P') < k ≤ d and some prime p | N exceeds 2B, then P = P'. The proof is unique factorisation in F_p[X]. This covers non-exact families with one *common* small side. Prop NX covers distinct small sides.

**ONE-RING (View B).**
- *Setting.* Suppose all relations of a family are exact in one ring A_F = Z[X]/(F), with F squarefree over Q. Their exponent vectors lie in K₀(F) = ker(Z^𝔅 → (Q[X]/F)^×), which depends on F only.
- *Squares.* An even combination gives x ∈ A_F with x² = c², so x/c is a sign vector on the Q-irreducible factors F_i of F. Its image mod N reveals at most gcd(N, F_i(μ)), which is computable from F alone.
- So in View B a single presentation gives *no* factoring information. This is systeq's EX in relation form. Theorem W adds that it gives at most one relation anyway, and Prop S does the same for View A, where single presentations *are* informative through N | F(m).

## 6. Main theorem, scope, and what it rules out (Step 5)

### 6.1 Theorem PB (PTE barrier) [PROVED]

**Hypotheses.** Fix N and a presentation F ∈ Z[X] \ {0} of any degree d and any height, possibly constructed from N, with φ(F) = 0 for φ as in §1.1. Let 𝓕 be a compressed split family on F with parameters (M, B, k, j, H, H_Q).

**Conclusions.**
- **(a) Exact families, any view.** If M ≥ K + 2⌊log₂(2^{K−1}√K·h·B^{K−1})⌋ + 4j + 3, with K = k + j₂ and h = 2·min(k, j₂+1)·H_Q·H, then 𝓕 yields at most one relation.
- **(b) PTE families** (shared e₁..e_{M−k}), exact, d ≥ k. 𝓕 has at most one member if M ≥ k + 2⌊log₂(2^k√k·H·B^{k−1})⌋ + 3. ptecons CW gives a sharper factorial version.
- **(c) View A** (φ(X) = m, roots within Y of m), exact, any j and any multiplier heights. If (k)·H·(|m|+Y)^{k−1}·(d+1)·Y^d < N, every reduced member uses only the ≤ d fixed pieces m − c, c ∈ S_F. With PTE structure and M − k ≥ d + 2k − 3, 𝓕 has at most one member.
- **(d) Non-exact PTE families, d ≥ k.** Members have pairwise distinct N-defects, so there is no shared congruence. Their number is at most the number of defect cells.
- **(e) k = 1 fibers, any constants.** At most C·B^{1/2}(log B)^κ + O(M) members, uniformly in M. This rests on the KNOWN uniform Bombieri–Pila bound.

### 6.2 What it rules out (the L-range)

**Setting.** Pieces ≤ Y = L[β] with β < 1/3, and M ≥ (1 − o(1))x/log Y ≍ x^{1−β}(log x)^{β−1}.

**View B (roots in [−B, B], log B ≍ log Y).** (a) and (b) hold as soon as

  (k + j)·log B + log(H·H_Q) = o(x^{1−β}).  (★)

(★) holds for:
- every BGJT-shaped family: k ≤ 4, j ≤ 3, heights B^{O(1)};
- any small side of degree k = o(x^{1−2β}) and log-height o(x^{1−β}).

Note that x^{1−2β} > x^{1/3} and x^{1−β} > x^{2/3}.

**View A (φ(X) = m, roots within Y of m).** Two cases.
- *d < k, or more generally d·log Y + (k−1)·log(|m|+Y) + log H < x − O(log d).* Then (c) applies, with no condition at all on the multipliers.
- *d ≥ k.* Then (b) applies after centring at m: substitute X ↦ X + m. The roots move to [−Y, Y] and the small-side heights become H' ≤ H·(2 + 2|m|)^{k−1}.
  - Condition (★) becomes k·log Y + (k−1)·log(2 + 2|m|) + log H = o(x^{1−β}).
  - If |m| ≥ 2H + 1, every small side satisfies |r(m)| ≥ |m|^{deg r}/2. The L[β] cost target log|v| = O(x^{2β}(log x)^{1−2β}) then forces (k−1)·log|m| = O(x^{2β}·polylog), and hence k·log Y = O(x^{2β}·polylog) whenever |m| ≥ Y^{Ω(1)}.
  - Both are o(x^{1−β}) because 3β < 1. So in View A, for |m| ≥ max(2H + 1, Y^{Ω(1)}) and log H = o(x^{1−β}), **every parameter compatible with the L[β] target is in the regime.**
- The excluded corner is O3: |m| tiny together with a small side of degree ≍ x^{2β}. It can matter only for β ≥ 1/4.

**Conclusion.** In model 𝓜_sys — exact (systematic) compressions over Z[X], one presentation, small polynomial side — route 1's open escape yields ≤ 1 relation per presentation, and in View A not even a dense one.
- An algorithm reaching L[α] relations must produce L[α] presentations, each N-conditioned separately.
- A presentation together with its split member *is* a relation found by other means.
- So the escape gives no route to α < 1/3, and no improvement of any constant, beyond the underlying relation-finding method. For the PTE family exactly as systeq §6 posed it — one prescribed top part, one presentation — this is unconditional.

### 6.3 What is not closed (open)

- **O1. Non-exact or evaluated compressions.**
  - *Non-exact.* By Prop NX each member pays its own N-congruence; whether many members of one PTE fiber can hit the defect cells is an equidistribution question we cannot settle.
  - *Evaluated* (View A, r small only in value at m, r = (X − m)S + v). This is the integer congruence u ≡ v mod F(m) with a common-offset structure u_A − u_{A'} = D(m). It is ptecons' NC primitive and smoothenum's doubly-smooth listing. It is *not* a BGJT compression; it is the original relation-finding problem.
- **O2. Linear systems of ≥ 2 unrelated presentations** (members P_A = F₀ + ν_A·F₁ + r_A). Each ν gives ≤ 1 member by Corollary P, but the number of ν could be large.
  - *Missing lemma:* rigidity for small integer three-term combinations c₀P₀ + c₁P₁ + c₂P₂ = G of split polynomials (G small, low degree).
  - The obvious pencils degenerate. P_ν = (X + ν)G is split for all ν, but has r = 0. With small r ≠ 0, Lemma R against (X + ν)G forces r = 0.
  - *We have no proof.* In View A, Prop S still makes each member use ≤ d + ρ pieces, so O2 matters only for d ≳ M.
- **O3. Corners.** k·log B ≳ x^{1−β}, which requires |m| tiny and log|v| ≈ x^{2β} with β ≥ 1/4. Or multiplier heights log H_Q ≳ x^{1−β} in View B. (View A needs no multiplier bound, by (c).)
- **O4. Other root sets.** Roots in non-discrete rings, and multivariate presentations (ptesplit_1's A_N with several variables).

## 7. Experiments against the pre-registration

**E1 (exhaustive R0: 318 cells; B = 1..5, M ≤ 30/22/16/13/11, K = 1..4).**
- **P1 holds: 0 violations.**
- *Disclosed: this is a weak test.* Only 51 of the 318 cells contain any pair with an equal top part (always with M − K ≤ 4).
- Only **2** cells are non-vacuous tests, i.e. have h_R ≥ 1 and contain a pair: (B, M, K) = (4, 4, 1) and (5, 4, 1). Both have h_min = 36 against h_R = 1.
- P2 ("h_min ≥ 10·h_R at the largest M of each row") therefore holds only vacuously.

**E1b (exhaustive R with j = 1: 147 cells; Q = c₁X + c₀, |c_i| ≤ 2).**
- **P1b holds: 0 violations.**
- 25 cells contain pairs (all with M ≤ 6), and **none** is a non-vacuous test (h_R = 0 wherever pairs exist).
- The j ≥ 1 case of Lemma R is therefore supported by its proof and by the random-pair certificate tests, not by E1b.

**E1c (POST HOC, disclosed; added because E1 was nearly vacuous).** 30 cells with larger B and small M, where PTE pairs exist: (B, M_max) = (8, 8), (12, 6), (16, 5), (30, 4).
- 0 violations, with 8 non-vacuous cells, all with K = 1.
- For K = 1 the minimal constant difference is 4, 36, 5040, 14400 at M = 3, 4, 5, 6, against h_R ≤ 3. These are the familiar ideal-PTE constants, which are divisible by large factorials, consistent with CW's factorial bound.
- For K = 2 (B = 8) the minima are 1, 2, 36, 432, 14400 at M = 3..7 (h_R = 0 throughout).
- **Bottom line:** the exhaustive runs found no bug, but at feasible sizes they barely reach the lemma's regime. The lemma rests on its proof.

**E2 (fiber statistics; same enumeration).**
- P3 (injection bound): 0 violations.
- P4 (k = 1 members pairwise disjoint): 0 violations.
- **P5 failed as expected-stated:** there are zero nontrivial fibers with s ≥ √(2M) + 2 for B ≤ 4. This is consistent with the heuristic count, not with our expectation of persistent structured solutions.

**E3 (Dickson presentation, (a, b) ∈ [−200, 200]²).** Exactly one split member, (a, b) = (1, 3): the systeq example itself. This agrees with the expectation.

**Resources.** One process, peak RSS 577 MB, 17 s for E1 + E1b + E3, about 40 s for E1c.

## 8. Reviews

- **Codex blind round (`ptebarrier_1.out`, before our analysis).** It reached the same model and the same partial barrier, with no full barrier. It contributed:
  - the uniform injection and packing bounds;
  - the pointwise Vinogradov bound J_{M/2,s};
  - the uniformity counterexamples;
  - the shared-defect example (X + 1)(X + 6)² − (X + 2)²(X + 9) = 8X at M = 3, which is out of regime for Lemma R and kept as a test;
  - the Hadamard payment for independent defects.

  It identified the missing step as "an exceptional-fiber bound incorporating the small-remainder window". Lemma R / Corollary P supply that step for exact members, and Prop NX gives the injectivity for non-exact members.
- **Codex adversarial round (`ptebarrier_2`).** Pending at the time of writing; see the addendum.
- **Opus proof critic and literature critic.** Pending; see the addendum.

## 9. Reproduction

```bash
.venv/bin/python -m src.ptebarrier.run_experiment all   # E1 (+E2 stats), E1b, E3  (~20 s)
.venv/bin/python -m src.ptebarrier.run_experiment E1c   # post-hoc larger-B check (~40 s)
.venv/bin/python -m pytest -q tests/test_ptebarrier.py
```
