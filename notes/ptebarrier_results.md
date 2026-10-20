# ptebarrier: a rigidity barrier for route 1's PTE escape

Agent `ptebarrier`, phase 16 (route 1, barrier side), 2026-10-20. **Revised** after codex ptebarrier_2 (MAJOR REVISION), the Opus proof critic (ACCEPT WITH FIXES) and the Opus literature critic (§8).

**Files.**
- Code: `src/ptebarrier/`.
  - `rigidity.py`: thresholds of Lemma R. A step-by-step certificate of its proof. The Theorem W invariant. The Prop S small-value set.
  - `enum.py`: exhaustive enumeration of split polynomials.
  - `run_experiment.py`: the runs.
- Tests: `tests/test_ptebarrier.py`.
- Pre-registration: `notes/ptebarrier_prereg.md` (commit d6d8073, made before any run).
- Config: `experiments/configs/ptebarrier.json`.
- Raw results:
  - pre-registered: `experiments/raw/ptebarrier_E{1,1b,3}.jsonl`;
  - post hoc, disclosed: `ptebarrier_E1c_posthoc.jsonl`.
- Codex (GPT-6.1-sol): blind round `codex/ptebarrier_1.{md,out}`; adversarial round `codex/ptebarrier_2.{md,out}`.

**Labels.**
- **PROVED**: proof given here.
- **PROVED (sketch)**: the argument is given but not every hypothesis is written out.
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
- ‖·‖ is the max-coefficient height.

---

## 0. Verdict

**Plain language.** The "dense split family" escape asks for *one* N-dependent presentation that makes *many* split (smooth-by-construction) polynomials reduce to small ones. That is how the finite-field algorithms get q³ relations from one presentation.

We prove that over the integers this amplification is impossible **when the reduction is an exact polynomial identity with small multipliers and a small polynomial remainder**, the faithful analogue of the BGJT identity. The reason is archimedean. Take two split integer polynomials of large degree with bounded roots. If they differ by a polynomial that is small and of low degree, they are equal; with small multipliers allowed, they are equal up to those multipliers. So such a family yields the same relation for all its members.

What we do **not** prove: a barrier for reductions that hold only modulo N (non-exact), or only in value at an integer. Those natural versions of systeq §6 remain open (O1). Our results there are partial, and one claimed reading in the draft ("every non-exact member pays its own congruence") was refuted by codex.

**Proved.**
- **Lemma R (rigidity) [PROVED, elementary, uniform in every parameter].**
  - *Statement.* Let P, P' be monic of degree M, split over Z with roots in [−B, B]. Let Q, Q' ≠ 0 have degree ≤ j, and let G = QP − Q'P' have degree < K and height ≤ h. Put L = ⌊log₂(2^{K−1}√K·h·B^{K−1})⌋. If M ≥ K + 2L + 4j + 3, then QP = Q'P'.
  - For j = 0 it is an exponential-strength relative of ptecons' CW and of the classical ideal-PTE divisibility (M − 1)! | C (Kleiman; Rees–Smyth), both sharper there.
  - The multiplier case j ≥ 1 is what covers BGJT-shaped families.
- **Theorem W (one presentation, one relation) [PROVED].**
  - *Setting.* F ∈ Z[X] is arbitrary: any degree, any height, possibly built from N. The members are exact identities h_A·P_A = F·Q_A + r_A, with deg h_A ≤ j₁, deg Q_A ≤ j₂, ‖Q_A‖ ≤ H_Q, r_A ≠ 0 of degree < k and height ≤ H.
  - *Conclusion.* Suppose Lemma R's threshold holds with K = k + j₂, j = j₁ + j₂ and h = 2·min(k, j₂+1)·H_Q·H. Then h_A·P_A/r_A is the same rational function for all members, so every non-degenerate member yields the same relation.
- **Corollary P (exact PTE families) [PROVED].**
  - *Setting.* Exact members P_A = F·Q_A + r_A with h_A = 1, sharing e₁..e_{M−k}, on a presentation of degree d ≥ k.
  - *Conclusion.* The family has at most one member once M ≥ k + 2⌊log₂(2^k√k·H·B^{k−1})⌋ + 3. The heights of F and of the quotient never enter.
- **Prop S (integer evaluation) [PROVED].**
  - *Setting.* N | F(m) ≠ 0, every root lies within Y of m, and (ρ+1)·H·(|m|+Y)^ρ·(d+1)·Y^d < N.
  - *Conclusion.* Every exact identity hP = FQ + r, with h and Q arbitrary and r of degree ≤ ρ and height ≤ H, has all roots of P in S_F ∪ Z(r), where |S_F| ≤ d and S_F is fixed by F.
  - *Consequence.* An exact PTE family through such an F, with *any* multipliers and quotients, has ≤ 1 member when M − k ≥ d + 2ρ − 1.
- **k = 1 fibers, arbitrary constants [PROVED from KNOWN].**
  - *Statement.* Members Q₀ − t (all split over Z, roots in [−B, B]) of a PTE fiber with one free coefficient number at most 2c·B^{1/2}(log B)^κ + 4, for M ≥ 4 and absolute c and κ.
  - *Source.* The squarefree count uses the uniform Bombieri–Pila bound of Binyamini–Cluckers–Novikov. Rolle's theorem allows at most two non-squarefree levels.
  - *Consequence.* This holds for exact and non-exact families alike, so one ideal-type fiber cannot fill a factor base of size ≍ B.

**Partial (non-exact).**
- *Prop NX [PROVED].* In the PTE setting with d ≥ k, distinct members have distinct reduced N-defects.
- *What it does not give.* Codex showed this does **not** imply independence. Two members can inherit compatibility from one identity P − P' ≡ small (mod N).
- *Remark NX1 [PROVED].* For k = 1, such an inherited compatibility is itself a split relation of the same kind: a 2B-smooth integer congruent to a small one.

**What the counting tools give.**
- **Vinogradov (BDG/Wooley).**
  - It bounds the second moment of fiber sizes (PTE collisions).
  - Its pointwise consequence is max fiber ≤ J_{M/2,s}(B) ≪_{M,s,ε} B^{M/2+ε} + B^{M−s(s+1)/2+ε}, with constants not uniform in (M, s).
  - In the compression range k ≪ M this is far weaker than the elementary injection bound C(2B+k, k).
  - It is the wrong tool: the barrier lives in the archimedean height of the remainder.
- **Dimension growth.** Its constants depend on the ambient dimension M and grow polynomially in the degree s! (assessment: useless here for k ≥ 2).

**What it rules out (§6).** Take β < 1/3, pieces ≤ L[β] and M ≍ x/log Y.
- *Scope.* Any exact compressed family with small-side degree, multiplier degree and log-heights satisfying (k + j)·log B + log(H·H_Q) = o(M) = o(x/log B). Every BGJT-shaped family qualifies: k ≤ 4, j ≤ 3, heights B^{O(1)}.
- *Conclusion.* Such a family provides ≤ 1 relation per presentation. The exact (systematic) version of route 1's PTE escape therefore gives no amplification.
- *Exact case closed.* For the exact PTE escape as posed (one prescribed top part, one presentation) at most one member exists. In View A this holds for every presentation degree, under the conditions of §6.2.
- *What remains.* Producing many presentations is, as far as we can see, no easier than producing relations. That is an assessment, not a reduction theorem.

**Not closed (open, §6.3).**
- **O1.** Non-exact or evaluated compressions.
- **O2.** Linear systems of ≥ 2 unrelated presentations.
- **O3.** Exact families with large quotient or multiplier degree j ≳ M/4 in View B, and in View A beyond Prop S. Also parameter corners such as tiny m with a high-degree small side.
- **O4.** Other root sets and multivariate presentations.

**Experiments (pre-registered).**
- E1/E1b check Lemma R exhaustively: 318 + 147 cells, **0 violations**. E1c (post hoc, larger B) adds 30 cells, also 0 violations.
- *Disclosed:* these are weak tests. Only 2 pre-registered cells and 8 post-hoc cells are non-vacuous. The lemma rests on its proof; both external reviewers re-derived it.
- The P5 expectation (structured fibers persisting above √(2M)) **failed**: there were none at B ≤ 4.
- E3: the Dickson presentation has exactly one split member in a 401² window.

**Assessment.** Probability < 1% that route 1's PTE escape yields sub-L[1/3]. Codex: about 2%. The surviving case O1 is where any construction would have to live: congruences mod N that are coherent across a family without being exact over Z. No candidate is known (cf. ptecons §3.5–3.6: NC, HALF).

---

## 1. Model and cost targets

### 1.1 Two views of a presentation

A **presentation** is F ∈ Z[X] \ {0} together with a ring homomorphism φ: Z[X] → Z/N with φ(F) = 0.
- **View B (residue).** φ(X) = μ ∈ Z/N. The linear factor base is {μ − a : |a| ≤ B} together with small primes. This is the literal analogue of F_{q²}[X] → F_{q^{2n}}.
- **View A (integer evaluation).** φ(X) = m ∈ Z, so N | F(m). The pieces m − a are integers, factored over the primes ≤ Y. This is systeq §6's (u, v) model and that of ptesplit/ptecons: u = P(m), v = r(m).

**A compressed split family on F** consists of members (A, h_A, Q_A, r_A):
- A is a multiset of M integers, all in [−B, B] (in View A, within Y of m; after centring, in [−Y, Y]).
- h_A, Q_A ∈ Z[X] \ {0} are multipliers.
- r_A ∈ Z[X] \ {0} is the small side, of degree < k and height ≤ H.

The member's relation is φ(h_A)φ(P_A) = φ(r_A). It is **non-degenerate** if φ(h_A), φ(P_A), φ(r_A) are units (View B), or the integers h_A(m), P_A(m), r_A(m) are nonzero and prime to N (View A). Otherwise a gcd either reveals a factor or the member carries no relation.
- The member is **exact** if h_A·P_A = F·Q_A + r_A holds in Z[X]. This is the "systematic" case. The word is meant in the BGJT sense, not that of the classical F_{2^n} systematic equations of Blake–Fuji-Hara–Mullin–Vanstone; BGJT's own footnote makes the same distinction. It is the analogue of the BGJT identity h₁P = (aX + b)(h₁X^q − h₀) + r.
- It is **non-exact** if the identity holds only modulo (N, F).
- It is **evaluated** (View A) if only P_A(m) ≡ r_A(m) (mod N) holds, with r_A small in value at m but not in height.

systeq §6 poses the escape as "a common reduction to a small polynomial modulo a presentation (N, F)". **That includes non-exact and evaluated members.** Our theorems cover the exact members. By Lemma EX (systeq §5.5), a non-exact member's reduced defect red_F(P_A) − r_A is a nonzero element of N·Z[X], and hence has height ≥ N.

### 1.2 Cost targets [HEURISTIC; random-residual model of ptesplit_1 and ptecons §1]

- *Split side.* Pieces ≤ Y = L[β] and a non-trivial integer congruence force M·log Y ≥ (1 − o(1))x, so M ≍ x^{1−β}(log x)^{β−1}.
- *Small side.* For an L[β] algorithm its log-size must be O(x^{2β}(log x)^{1−2β}).
- *Family.* It must supply ≳ L[β] relations at amortised cost L^{o(1)}.

## 2. Heuristic count (Step 1) [HEURISTIC + measured]

With q = 2B + 1 there are C(q+M−1, M) ~ q^M/M! multisets (for M² = o(B)). The vector (p₁, …, p_s) has a range of size about M^s·B^{s(s+1)/2}. Absorb the moment-body density and the singular series into a factor 𝒜_{M,s}(c), as in codex ptebarrier_1 eq. (1). The expected fiber is then

  λ(c) ≈ (q^M/M!)·B^{−s(s+1)/2}·𝒜_{M,s}(c).

λ drops below 1 at s(s+1)/2 ≈ M − log(M!)/log B, i.e. s ≈ √(2M) when log M ≪ log B. So dense fibers generically need k ≥ M − O(√M) free coefficients: no compression.

A compressed family also needs its k low coefficients to land in a window of height H. These coefficients have size up to C(M,i)·B^{M−i}. That costs a further factor of order ∏_{i<k} min(1, H/(C(M,i)B^{M−i})). Lemma R makes this rigorous for exact families: the count is 1.

**Measured (E2, B ≤ 4).**
- No nontrivial fiber has s ≥ √(2M) + 2.
- Max fibers sit far below the injection bound. For example, (B, M, k) = (4, 6, 4) gives 10 against 495.

The pre-registered expectation P5, that structured fibers would persist above √(2M), **failed**. Known ideal PTE solutions need entries in the hundreds.

## 3. Vinogradov's mean value theorem (Step 2)

Put R(c) = #{a ∈ [1, B]^M : p_i(a) = c_i, i ≤ s}, counting *ordered* tuples, and f(α) = Σ_{a≤B} e(α₁a + … + α_s a^s).

- **KNOWN.**
  - Bourgain–Demeter–Guth (Ann. Math. 184 (2016)) proved J_{t,s}(B) = ∫|f|^{2t} ≪_{t,s,ε} B^{t+ε} + B^{2t − s(s+1)/2 + ε}.
  - Wooley's efficient congruencing gives the cubic case (2016) and all degrees by nested efficient congruencing (Proc. LMS 2019).
- **What it bounds [PROVED].**
  - Σ_c R(c)² = J_{M,s}(B), the number of ordered PTE collisions.
  - Hence #{c : R(c) ≥ T} ≤ min(B^M/T, J_{M,s}/T²): *most* fibers are small.
  - #{c : R(c) > 0} ≥ B^{2M}/J_{M,s}.
- **Pointwise [PROVED, codex ptebarrier_1].**
  - R(c) = ∫ f^M e(−c·α) ≤ ∫|f|^M = J_{M/2,s}(B) for even M, with Cauchy–Schwarz between adjacent moments for odd M.
  - So max R ≪_{M,s,ε} B^{M/2+ε} + B^{M − s(s+1)/2+ε}.
  - Cauchy–Schwarz on the second moment (max R ≤ J_{M,s}^{1/2}) is weaker. Already the diagonal gives J_{M,s}^{1/2} ≥ (M!·(B)_M)^{1/2}.
- **Comparison [PROVED].** For multisets, the injection bound gives max fiber ≤ C(q+k−1, k) ≤ q^k: keep the k largest entries, and Newton recovers the rest. For multisets, Vinogradov's bound becomes R/(orbit size) ≫ B^{M/2}/M! up to constants. It beats q^k only when roughly k log q > (M/2)·log B − log M!, i.e. only far outside the compression range k ≪ M.
- **Uniformity [PROVED + KNOWN].** The implied constants depend on (M, s, ε). Since log M! ≫ x^β in our range, they cannot be suppressed. Explicit versions carry exp(O(s³ log s))-size constants: Wooley, "On Vinogradov's mean value theorem II" (1993); Ford's explicit Theorem 3; Steiner's efficient-boxing estimates.
- **Assessment.** VMVT shows that dense fibers are *rare* (collision counts). The escape needs only *one* dense N-compatible fiber, so VMVT is the wrong tool. This does not conflict with its standard use for Tarry's problem, which is about existence thresholds (Hua; Wooley 1996).

## 4. Determinant method / dimension growth (Step 3)

**Geometry [PROVED].** V_c = {e₁ = c₁, …, e_s = c_s} ⊂ A^M is a complete intersection of dimension k = M − s and degree ≤ s!.

**Uniform bounds [KNOWN; assessment].** Castryck–Cluckers–Dittmann–Nguyen (Algebra Number Theory 14 (2020)) give N(X, B) ≤ c(n)·D^{e(n)}·B^{dim X} for *integral* projective X ⊂ P^n_Q of degree D ≥ 5; V_c may be reducible, so the bound is applied per component of the projective closure (≤ s! components). BCN Thm 3 extends dimension growth to D = 3, 4. The constants depend on the ambient dimension n = M, and D = s! ≈ e^{M log M} = L[1 − β]^{Θ(1)}. By comparison with every L[β] budget these bounds are useless for k ≥ 2. This is an assessment about known bounds, not a theorem.

**k = 1 [PROVED from KNOWN].** The members are Q₀ − t, split over Z with roots in [−B, B], M ≥ 4.

1. **Non-squarefree levels.** All the polynomials Q₀ − t in question are real-rooted. The set of t for which Q₀ − t is real-rooted is an interval [t₋, t₊], from the maximal local minimum to the minimal local maximum (Rolle and interlacing). A repeated root occurs only at a critical value, and inside the interval only at its endpoints. So **at most two members are non-squarefree** (codex ptebarrier_2).
2. **Squarefree members.** Every ordered pair (x, y) of distinct roots of one member lies on the plane curve C: (Q₀(x) − Q₀(y))/(x − y) = 0, of degree M − 1. Distinct t give disjoint pairs.
3. **Rational lines.** A rational line y = λx + ν on C forces Q₀(λx + ν) ≡ Q₀(x), hence λ^M = 1, so λ = ±1. λ = 1 forces ν = 0, and the diagonal is not on C. So C has at most one rational line, y = ν − x, which carries ≤ M pairs per member.
4. **Other components.**
   - Geometric components not defined over Q carry ≤ e² rational points (Bézout with a conjugate).
   - Rational components of degree e ≥ 2 carry ≤ c·e²·B^{1/e}(log B)^κ points of height ≤ B. This is Binyamini–Cluckers–Novikov, "Bounds for rational points on algebraic curves, optimal in the degree, and dimension growth" (arXiv 2302.04209), affine bound Theorem 2 / eq. (27), with c = c(2) and κ absolute.
5. **Sum.** Over components, Σe ≤ M − 1. Hence #squarefree·M(M − 2) ≤ c(M−1)²B^{1/2}(log B)^κ + (M−1)².
6. **Total.** **≤ 2c·B^{1/2}(log B)^κ + 2 + 2 members**, uniformly in M ≥ 4.

This holds for any constants, exact or not. For generic Q₀ the curve C is irreducible (Davenport–Lewis–Schinzel 1961; Fried 1973), which would give B^{1/(M−1)}; the stated B^{1/2} is the uniform worst case. M ≥ 3 is necessary: {X² − a²} has B + 1 members.

## 5. Rigidity (Steps 4–5)

### 5.1 Lemma R [PROVED]

**Statement.**
- Let P, P' ∈ Z[X] be monic of degree M, split over Z, with all roots in [−B, B] (B ≥ 1).
- Let Q, Q' ∈ Z[X] \ {0} have degree ≤ j.
- Put G := QP − Q'P', with K ≥ 1, deg G ≤ K − 1 and ‖G‖ ≤ h (if h < 1 then G = 0).
- Let Δ := 2^{K−1}·√K·h·B^{K−1} and L := ⌊log₂ Δ⌋.

Then:
- if M ≥ K + 2L + 4j + 3, we have QP = Q'P';
- for j = 0 (Q = Q' = 1) the condition M ≥ K + 2L + 3 suffices.

**Proof.** Suppose G ≠ 0.

1. **Coprime parts.** Write P = gU and P' = gU', where g is the monic gcd (the common sub-multiset) and U, U' are coprime and split, of degree r = M − deg g. Then G = g·G̃ with G̃ := QU − Q'U' ≠ 0, and deg g ≤ deg G ≤ K − 1, so r ≥ M − K + 1.
2. **Degrees.** If deg(QU) ≠ deg(Q'U'), or the leading coefficients differ, then deg G̃ ≥ r and deg G ≥ M > K − 1. So deg Q = deg Q' =: j' and lc(Q) = lc(Q').
3. **Height.**
   - Mahler measure is multiplicative and M(g) ≥ 1, so M(G̃) ≤ M(G) ≤ ‖G‖₂ ≤ √K·h.
   - Each coefficient of a degree-n polynomial f satisfies |c_i| ≤ C(n, i)·M(f) (Mahler 1960). So ‖G̃‖₁ ≤ 2^{K−1}√K·h, and |G̃(a)| ≤ Δ for every integer |a| ≤ B.
   - Dividing by g need not preserve coefficient height; the Mahler measure avoids that issue.
4. **Local smallness.**
   - A root a of U is *good* if Q'(a) ≠ 0. There G̃(a) = −Q'(a)U'(a), with Q'(a) ∈ Z \ {0} and U'(a) ≠ 0 (coprimality). So 1 ≤ |U'(a)| = ∏_{a'}|a − a'| ≤ Δ.
   - All factors are nonzero integers, and each factor ≥ 2 doubles the product. So **at most L roots of U' (with multiplicity) lie at distance ≥ 2 from a**: supp U' ⊆ {a − 1, a + 1} ∪ E with |E| ≤ L.
   - The same holds with the roles exchanged (good roots of U' are those with Q(a') ≠ 0).
5. **Moments.**
   - QU/lc and Q'U'/lc are monic of degree r + j' and differ in degree ≤ K − 1 − deg g. So e₁..e_J of their complex root multisets agree, with J = r + j' − K + deg g = M + j' − K ≥ M − K.
   - By Newton (characteristic 0), p₁..p_J agree.
   - Let μ = Σ_{roots(QU)} δ − Σ_{roots(Q'U')} δ. Then μ ≠ 0 (equal root multisets and equal leading coefficients would give G = 0), and its moments of orders 0..J vanish.
   - Vandermonde over C gives **|supp μ| ≥ J + 2 ≥ M − K + 2.** This is the classical fact that a degree-J PTE solution has size ≥ J + 1 (Borwein, *Computational Excursions*, Ch. 11), here with complex multiplier roots allowed.
6. **Support count.** supp μ ⊆ supp U ∪ supp U' ∪ Z_C(Q) ∪ Z_C(Q'); the last two sets have ≤ j' elements each.
   - *Case 1: both U and U' have good roots.* |supp μ| ≤ 4 + 2L + 2j.
   - *Case 2: U has no good root and U' does.* Then supp U ⊆ Z(Q').
     - Two good roots of U' at distance ≥ 3 would put ≥ r − L roots of U into each of two disjoint sets {a'_i ± 1}. That gives r ≤ 2L, contradicting r ≥ M − K + 1 ≥ 2L + 4.
     - So the good roots of U' take ≤ 3 values, and |supp μ| ≤ 4j + 3. The symmetric case is the same.
   - *Case 3: neither has a good root.* |supp μ| ≤ 4j.
   - The proof critic notes that the true Case 2/3 bounds are 2j + 3 and 2j (we double-count Z(Q) and Z(Q')). So M ≥ K + 2L + 2j + 3 already suffices. We keep the conservative constant.
7. **Contradiction.** M − K + 2 ≥ 2L + 4j + 5 exceeds every case bound. For j = 0 only Case 1 occurs. ∎

**Remarks.**
- **The reason is archimedean.** In F_q[X] there is no size, so the BGJT family exists (Joux 2013; Göloğlu–Granger–McGuire–Zumbrägel 2013; BGJT 2014). Compare GKZ, "Indiscreet logarithms", closing open-problems paragraph: for prime fields "fundamentally new ideas will be required". The step-4 argument uses only that Z is discrete in R.
- **Other root sets.** For roots in Z[i] or Z[ω], the analogue should hold with lattice-specific constants: 8 resp. 12 nonzero points of modulus < 2, a revised Case 2, and |a| ≤ √2·B in step 3. This is HEURISTIC (not written out). It fails for non-discrete rings such as Z[√2].
- **j = 0: prior and sibling results.**
  - For the ideal case K = 1, the classical result (M − 1)! | C (Kleiman; Rees–Smyth, Prop. 2.1, as attributed in Caley 2013 §5) is far sharper.
  - For general K, ptecons' CW (codex ptecons_1) gives a factorial bound by a positivity argument. CW is sharper than our exponential bound.
  - Lemma R's role is the multiplier case j ≥ 1. There positivity of ∏(X − a)² fails at the complex roots of Q and Q', but the Vandermonde count still works.
  - Both thresholds have the shape M ≳ K·log B.
- **Reviews.** Codex ptebarrier_2 (CORRECT, with no missing case, after a 42-cell multiplier search with 0 violations) and the Opus proof critic (all steps correct) each re-derived every step.
- **Certificate.** `rigidity.lemma_R_certificate` recomputes steps 1–6 on any pair. It is tested on 800 random pairs with j ∈ {0, 1}.

### 5.2 Theorem W: one presentation, one relation [PROVED]

**Statement.** Let F ∈ Z[X] \ {0} be arbitrary. Consider members (A, h_A, Q_A, r_A) with:
- h_A·P_A = F·Q_A + r_A in Z[X];
- h_A ≠ 0 of degree ≤ j₁;
- Q_A ≠ 0 of degree ≤ j₂ and ‖Q_A‖ ≤ H_Q;
- r_A ≠ 0 of degree < k and ‖r_A‖ ≤ H;
- roots in [−B, B].

Put K = k + j₂, j = j₁ + j₂ and h = 2·min(k, j₂+1)·H_Q·H. If M ≥ K + 2⌊log₂(2^{K−1}√K·h·B^{K−1})⌋ + 4j + 3, then:
- **(i)** the rational function Ψ_A := h_A·P_A/r_A ∈ Q(X)^× is the same for all members;
- **(ii)** all non-degenerate members define the same relation, i.e. the same exponent vector:
  - over ±1, the rational primes and the primitive irreducibles of Z[X] (View B: φ(Ψ) = 1 in (Z/N)^×);
  - over ±1 and the primes after evaluation at m (View A: Ψ(m) ∈ Q^×).

So the family's relation lattice has rank ≤ 1.

**Proof.**
1. For two members, Q_{A'}(h_A P_A) − Q_A(h_{A'}P_{A'}) = Q_{A'}r_A − Q_A r_{A'} =: G.
2. Here deg G < k + j₂. Each product coefficient has ≤ min(k, j₂+1) summands, so ‖G‖ ≤ h.
3. The multipliers Q_{A'}h_A and Q_A h_{A'} are nonzero, of degree ≤ j.
4. Lemma R gives G = 0, so Q_{A'}r_A = Q_A r_{A'}, and therefore Ψ_A = Ψ_{A'}.
5. Q(X)^× = Q^× × (free group on monic irreducibles), and Q^× = ±1 × (free group on the primes). So equal rational functions have equal exponent vectors. Evaluation of equal rational functions at m gives equal rational numbers, whenever numerators and denominators are nonzero there. ∎

**Notes.**
- *What never enters.* The heights of h_A and of F never enter; F may have N-sized coefficients.
- *BGJT shape.* The BGJT identity h₁P = (aX + b)(h₁X^q − h₀) + [(aX + b)h₀ + (cX + d)h₁] has j₁ = 2, j₂ = 1, k = 4. In F_q[X] it gives q³ distinct Ψ. Over Z, Theorem W allows one.
- *Forced structure.* ρ := r_A/Q_A is common to all members, so h_A P_A = Q_A·(F + ρ). If deg F ≥ k, the numerator Φ₀ of F + ρ has degree ≥ deg F ≥ M − j₂. Each P_A contains a sub-factor of Φ₀ of degree ≥ M − j, not necessarily the same one. The members are one relation multiplied by trivial identities.
- *A cautionary example (codex).* Members on one F may share the same Ψ while having different P_A. For example, P'·P = F − 1 = P·P' with F = PP' + 1 has two exact members with multipliers h = P', P. This agrees with (i), but it shows why Corollary P needs h_A = 1 to conclude "one member".

### 5.3 Corollary P: exact PTE families [PROVED]

**Statement.** Let the members share e₁..e_{M−k}, let F have degree d ≥ k, and let the members be exact with h_A = 1: P_A = F·Q_A + r_A, with deg r_A < d and ‖r_A‖ ≤ H. Then:
- Q_A is common to all members;
- P_A − P_{A'} = r_A − r_{A'} has degree < k and height ≤ 2H;
- there is at most one member once M ≥ k + 2⌊log₂(2^k√k·H·B^{k−1})⌋ + 3.

**Proof.**
- F(Q_A − Q_{A'}) = (P_A − P_{A'}) − (r_A − r_{A'}) has degree < d, so Q_A = Q_{A'}. This holds for any nonzero F, monic or not.
- Apply Lemma R with j = 0. ∎

**Remarks.**
- In the threshold one may take k = 1 + max deg(r_A − r_{A'}).
- The Dickson identity (X² − 4)(X² − 1)² = F + (X + 3), with F = X⁶ − 6X⁴ + 9X² − X − 7, has k = 2 and M = 6. That is out of regime, and E3 finds it to be the only split member with |a|, |b| ≤ 200.

### 5.4 Prop S: integer evaluation [PROVED]

**Statement.** Let m ∈ Z and F ∈ Z[X] of degree ≤ d, with F(m) ≠ 0 and N | F(m). Let Y, H ≥ 1 and ρ ≥ 0. Put B_m = |m| + Y and suppose (ρ+1)·H·B_m^ρ·(d+1)·Y^d < N. Let S_F := {a ∈ Z : |a − m| ≤ Y, |F(a)| ≤ (ρ+1)H·B_m^ρ}. Then:
- (i) |S_F| ≤ d;
- (ii) for every identity hP = FQ + r in Z[X], with P monic and split, all roots within Y of m, h and Q arbitrary, and r ≠ 0 of degree ≤ ρ and height ≤ H, every root of P lies in S_F ∪ Z(r).

**Proof.**
- (i) Suppose d + 1 distinct points of S_F existed. Lagrange interpolation at those points gives |F(m)| ≤ (d+1)·(ρ+1)H·B_m^ρ·Y^d < N ≤ |F(m)|, a contradiction. Each Lagrange coefficient has absolute value ≤ Y^d, since the denominators are nonzero integers.
- (ii) At a root a, F(a)Q(a) = −r(a). If r(a) ≠ 0, then F(a) and Q(a) are nonzero integers and |F(a)| ≤ |r(a)| ≤ (ρ+1)H·B_m^ρ. ∎

**Corrections from review.**
- The draft claimed that after cancelling gcd(P, r) every root lies in S_F. That is **false**: cancellation can raise the remainder's height. Codex's counterexample is F = 200(X − 1) − (1 + X + … + X⁴)², P = (X − 1)³, m = 0, N = 201.
- The draft's density consequences ("u uses only d pieces"; "at d ≍ M the family is a PTE fiber") overreached and are withdrawn.

**Valid consequence (PTE).** Take two members of an exact PTE family through such an F, with arbitrary h, Q. Their difference measure is supported on S_F ∪ Z(r_A) ∪ Z(r_{A'}), which has ≤ d + 2ρ points, and moments 0..M−k vanish. So **there is at most one member when M − k ≥ d + 2ρ − 1**. This covers View A presentations with d < k, where Corollary P is silent, and also families with arbitrary multipliers h_A.

### 5.5 Non-exact members and single rings

**Prop NX (defects are injective) [PROVED].**
- *Setting.* PTE setting, F monic of degree d ≥ k, members with P_A ≡ r_A mod (N, F), ‖r_A‖ ≤ H < N/2, and M ≥ k + 2⌊log₂(2^k√k·H·B^{k−1})⌋ + 3. Put D_A := red_F(P_A) − r_A ∈ N·Z[X]_{<d} (Lemma EX).
- *Injectivity.* D_A − D_{A'} = (P_A − P_{A'}) − (r_A − r_{A'}). If D_A = D_{A'}, Lemma R (j = 0, height 2H) gives A = A'. So A ↦ D_A is injective.
- *Cell bound.* Counting the differences D_A − D_{A₀}: their coefficient of degree i < k is bounded by 2C(M,i)B^{M−i} + 2H. Hence |𝓕| ≤ ∏_{i<k}(2⌊(2C(M,i)B^{M−i} + 2H)/N⌋ + 1), and the coefficients of degree ≥ k must agree mod N. This bound is weak.
- *Hadamard (codex ptebarrier_1).* r linearly independent defects satisfy N^r ≤ r^{r/2}·∏‖D_j‖.

**What NX does not give (codex ptebarrier_2, MAJOR, accepted).**
- Injective defects do **not** mean independent congruence conditions. A shared congruence can be inherited from one identity.
- *Example.* P = X²(X² − 49)², P' = (X² − 64)(X² − 25)(X² − 9), with P − P' = 14400, N = 25, F = P + 26, μ = 4, r = r' = −1. Both members are non-exact; they have distinct defects (−25, −14425) and are compatible because 25 | 14400.
- The draft's sentence "every member pays its own congruence" is withdrawn.

**Remark NX1 (inheritance costs a relation, k = 1) [PROVED].**
- *Statement.* If P − P' = C is constant and the members share a non-exact compression, then C ≡ r − r' (mod N) with |r − r'| ≤ 2H.
- *Why that is a relation.* For any root b of P', C = P(b) = ∏_{a∈A}(b − a). So C is a product of M integers of absolute value ≤ 2B (ptesplit_1 codex), and the inheritance is itself a congruence "smooth by construction ≡ small (mod N)" of the same type as the relations sought.
- *k ≥ 2.* The difference D = P − P' is a polynomial, and no such reduction is proved.

**Same remainder [PROVED, ptesplit_1 codex].** If P ≡ P' mod (N, F) with deg(P − P') < k ≤ d and some prime p | N exceeds 2B, then P = P' (unique factorisation in F_p[X]).

**ONE-RING (View B) [PROVED (sketch)].**
- *Hypotheses.* F is monic and squarefree over Q, with Q-irreducible factors F_i. The factor-base elements are units in Q[X]/F. N is prime to the denominators of the idempotents of Q[X]/F, which divide resultants of the F_i.
- *Claim.* Exact relations then lie in K₀(F) = ker(Z^𝔅 → (Q[X]/F)^×), which depends only on F.
- *Squares.* An even combination gives x with x² = 1 in ∏Q[X]/F_i, i.e. a sign vector. Its image mod N reveals at most gcd(N, F_i(μ)) and gcds with those resultants: information computable from F.
- Codex's example F = 15(3X+1)(5X+1) shows why contents and leading coefficients must be handled. That is why we assume F monic.
- So in View B a single presentation carries no factoring information (systeq EX in relation form). Theorem W adds that it carries at most one relation anyway.

## 6. Main theorem, scope, and what it rules out (Step 5)

### 6.1 Theorem PB [PROVED]

Fix N and a presentation F ∈ Z[X] \ {0} of any degree d and any height (possibly built from N), with φ(F) = 0.

- **(a) Exact families, any view.** Under the hypotheses of Theorem W, the family yields ≤ 1 relation.
- **(b) Exact PTE families with h_A = 1 and d ≥ k.** Under the threshold of Corollary P, ≤ 1 member.
- **(c) View A, exact PTE families with any multipliers.** If (ρ+1)H(|m|+Y)^ρ(d+1)Y^d < N and M − k ≥ d + 2ρ − 1, with ρ = max deg r_A: ≤ 1 member.
- **(d) Non-exact PTE families** (F monic, d ≥ k, regime of NX). Reduced defects are pairwise distinct; the family size is at most the number of defect cells. Independence is **not** implied.
- **(e) k = 1 fibers, M ≥ 4, any constants.** At most 2c·B^{1/2}(log B)^κ + 4 members.

### 6.2 The L-range

**Setting.** Pieces ≤ Y = L[β] with β < 1/3, and M ≍ x/log Y = x^{1−β}(log x)^{β−1}. The relevant scale is **M, not x^{1−β}**: the draft's condition o(x^{1−β}) was off by (log x)^{1−β} (proof critic and codex, MAJOR, fixed).

**(a)/(b), View B, roots in [−B, B] with log B ≍ log Y.** They apply when

  (k + j)·log B + log(H·H_Q) = o(x/log B).  (★)

(★) holds for:
- every BGJT-shaped family: k ≤ 4, j ≤ 3, heights B^{O(1)};
- generally, k + j = o(x/log² B) = o(x^{1−2β}(log x)^{2β−2}) and log(H·H_Q) = o(x/log B).

**View A.**
- **(c)** covers exact PTE families whenever d·log Y + ρ·log(|m| + Y) + log H < x − O(log d) and d + 2ρ ≤ M − k + 1. This holds for d up to roughly (1 − o(1))·x/log Y with ρ, log H small. It needs no condition on the multipliers.
- **(b) after centring at m** (X ↦ X + m): the roots move to [−Y, Y], the heights become ≤ H·(2 + 2|m|)^{ρ}, and (★) becomes (ρ+1)·log Y + ρ·log(2 + 2|m|) + log H = o(x/log Y).
  - Choose ρ = max deg r_A (only deg(P_A − P_{A'}) ≤ ρ enters).
  - If |m| ≥ 2H + 1, then |r_A(m)| ≥ |m|^{deg r_A}/2. The cost target log|v| = O(x^{2β}(log x)^{1−2β}) then gives ρ·log|m| = O(x^{2β}(log x)^{1−2β}).
  - If moreover |m| ≥ Y^{c} for a fixed c > 0, then ρ·log Y = O(x^{2β}(log x)^{1−2β}). That is o(x^{1−β}(log x)^{β−1}) because 3β < 1.
  - So, for |m| ≥ max(2H + 1, Y^c) and log H = o(x/log Y), every exact PTE family compatible with the L[β] target is in the regime of (b) or (c).
  - (For d ≥ k with h_A ≠ 1 one needs (c) or Theorem W with centred multiplier heights.)

**Conclusion (what is ruled out).**
- In the **exact, small-coefficient model** above, route 1's PTE escape yields ≤ 1 relation per presentation (≤ 1 member for PTE families with h_A = 1, or in View A).
- The BGJT amplification mechanism — one presentation, many split members, small polynomial remainders — therefore **does not exist over Z** in the L-range β < 1/3 (indeed for any β < 1/2 at bounded k and j).
- An L[α] algorithm built on route 1 would need L[α] separately constructed presentations, or one of the open mechanisms below.
- We do **not** prove that producing presentations is as hard as producing relations. That is an assessment.

### 6.3 What is not closed (open)

- **O1. Non-exact and evaluated compressions.** This is systeq's "modulo (N, F)" and "small value at m" versions.
  - Non-exact members have pairwise distinct defects (NX). But they can share compatibility through one difference identity (codex's example). For k = 1 that identity is itself a split relation (NX1); for k ≥ 2 nothing is proved.
  - Evaluated members (r small in value only) are the common-offset structure u_A − u_{A'} = D(m) of ptecons' NC heuristic.
  - This is where any construction would have to live: congruences mod N that are coherent across a family without being exact over Z. **No barrier and no candidate.**
- **O2. Linear systems of ≥ 2 unrelated presentations** (P_A = F₀ + ν_A F₁ + r_A).
  - Each ν gives ≤ 1 member by Lemma R with j = 0.
  - *Missing lemma:* rigidity for small integer three-term combinations c₀P₀ + c₁P₁ + c₂P₂ = G of split polynomials with G small and of low degree. Only the obvious pencils, such as (X + ν)G, are seen to degenerate.
  - Theorem W gives no lower bound on producing compatible presentations jointly.
- **O3. Exact families outside the regime.**
  - Large multiplier or quotient degree j ≳ M/4 in View B (proof critic #2). In View A, Prop S still bounds the *distinct roots* of each member by d + ρ, but not the number of relations.
  - Large H_Q in View B.
  - The View A corner |m| ≤ max(2H + 1, Y^{o(1)}) with ρ·log Y ≳ M. There the value target does not bound the degree: at m = 0, r = X^t + 1 has value 1 for every t.
- **O4. Other root sets and multivariate presentations.** Roots in non-discrete rings such as Z[√2]. Roots in imaginary quadratic orders (Remark §5.1, heuristic). For prior work on PTE over such rings see Alpers–Tijdeman 2007, Caley 2013 and Croot–Mao–Yip 2026. Multivariate presentations such as ptesplit_1's A_N in several variables.

## 7. Experiments against the pre-registration

**E1 (exhaustive R0: 318 cells; B = 1..5, M ≤ 30/22/16/13/11, K = 1..4).**
- **P1 holds: 0 violations.**
- *Disclosed: this is a weak test.* Only 51 of the 318 cells contain any pair with an equal top part, always with M − K ≤ 4.
- Only **2** cells are non-vacuous: (B, M, K) = (4, 4, 1) and (5, 4, 1), each with h_min = 36 against h_R = 1.
- P2 holds only vacuously.

**E1b (exhaustive R with j = 1: 147 cells; Q = c₁X + c₀, |c_i| ≤ 2).**
- **P1b holds: 0 violations.**
- 25 cells contain pairs (all with M ≤ 6), and none is a non-vacuous test.
- Codex's independent 42-cell multiplier search found 10,985 low-degree differences and 0 violations, also outside the regime.

**E1c (POST HOC, disclosed; added because E1 was nearly vacuous).** 30 cells: (B, M_max) = (8, 8), (12, 6), (16, 5), (30, 4).
- 0 violations, with 8 non-vacuous cells, all K = 1.
- For K = 1 the minimal constant differences are 4, 36, 5040, 14400 at M = 3..6, against h_R ≤ 3. These are the ideal-PTE constants, consistent with (M − 1)! | C.
- For K = 2 (B = 8) the minima are 1, 2, 36, 432, 14400 at M = 3..7.

The exhaustive runs found no bug, but at feasible sizes they barely reach the regime. The lemma rests on its proof.

**E2 (fiber statistics).**
- P3 (injection bound): 0 violations.
- P4 (k = 1 members pairwise disjoint): 0 violations.
- **P5 failed**: there are zero nontrivial fibers with s ≥ √(2M) + 2 at B ≤ 4. This is consistent with the heuristic count, not with our expectation.

**E3 (Dickson presentation, (a, b) ∈ [−200, 200]²).** Exactly one split member, (1, 3), as expected.

**Resources.** One process, peak RSS 577 MB, 17 s for E1 + E1b + E3, about 40 s for E1c.

## 8. Reviews

**Codex blind round (`ptebarrier_1.out`).** Same model, a partial barrier, no full barrier. It contributed:
- the injection and packing bounds;
- the pointwise Vinogradov bound J_{M/2,s};
- the uniformity counterexamples;
- the shared-defect example at M = 3 (out of regime, kept as a test);
- the Hadamard bound.

It named the missing step: "an exceptional-fiber bound incorporating the small-remainder window". Lemma R and Corollary P supply it for exact members.

**Codex adversarial round (`ptebarrier_2.out`): MAJOR REVISION.** Lemma R: CORRECT, no missing case, 0 violations in a 42-cell search. All findings were accepted and applied:
- Theorem W's specialisation needs unit or nonzero images; contents and signs are now included.
- Prop S's cancellation claim is false (counterexample); the density consequences overreached. Both are withdrawn.
- NX injectivity does not imply independence (counterexample); the reading is withdrawn and Remark NX1 added.
- ONE-RING needs monic F and the unit and denominator hypotheses; it is relabelled as a sketch.
- k = 1 needed squarefreeness; fixed by the two-endpoint Rolle argument. M ≥ 3 is necessary. There is one *rational* line.
- PB(b) needs h_A = 1 (counterexample F = PP' + 1).
- (★) had the wrong scale (fixed). The View-A degree was fixed via the maximum actual degree.
- *FATAL to the claimed closure:* 𝓜_sys excludes the non-exact and evaluated versions. The draft's "unconditional for the escape exactly as posed" is retracted; the barrier is now stated for the exact model only, with O1 open.
- Prior art: (M − 1)! | C; Borwein–Lisoněk–Percival.

**Opus proof critic: ACCEPT WITH FIXES, no FATAL.** All steps of Lemma R, Theorem W, Corollary P and Prop S(i)/(ii) were checked. Findings applied:
- the (★) polylog error;
- the missing exact large-j case (now O3);
- the false §5.4 PTE-fiber claim;
- the narrowing of "as posed";
- the k = 1 squarefree gap;
- the threshold degree k = 1 + max deg r;
- the unsupported "β ≥ 1/4" corner claim (removed);
- the NX cell factor 2;
- View-A non-degeneracy;
- the ONE-RING label;
- the §0/§5.2 parameter mismatch;
- the "forced structure" wording;
- the Gaussian-roots label;
- the O2 citation;
- the §3/§4 labels.

**Opus literature critic: OK WITH FIXES, no overclaimed novelty.**
- No published lower bound was found on ‖∏(X−a_i) − ∏(X−b_i)‖ for differences of degree < K with K > 1, with or without multipliers. Searched: Rees–Smyth, Caley's paper and thesis, Borwein–Ingalls, Borwein–Lisoněk–Percival, Filaseta–Markovich, Chen's 2025 survey, Croot–Mao–Yip 2026.
- No published integer analogue of BGJT compression and no published obstruction were found.
- Every citation checked is accurate in substance, including BCN eq. (27) and the Caley §5 attribution.
- Fixes applied:
  - the reference list (§10), which was the MAJOR finding;
  - attribution of the classical ingredients of Lemma R (Mahler's coefficient bound; PTE size ≥ degree + 1);
  - the CCDN integrality and per-component application;
  - a disambiguation of "systematic";
  - BGJT/Joux/GGMZ citations;
  - Davenport–Lewis–Schinzel/Fried for k = 1;
  - prior work for O4 (Alpers–Tijdeman, Caley, Croot–Mao–Yip);
  - the GKZ citation form.

## 9. Reproduction

```bash
.venv/bin/python -m src.ptebarrier.run_experiment all   # E1 (+E2 stats), E1b, E3  (~20 s)
.venv/bin/python -m src.ptebarrier.run_experiment E1c   # post-hoc larger-B check (~40 s)
.venv/bin/python -m pytest -q tests/test_ptebarrier.py
```

## 10. References

**Smooth-by-construction families in finite fields (the mechanism we test).**
- Joux, "A new index calculus algorithm with complexity L(1/4+o(1)) in small characteristic", SAC 2013 (ePrint 2013/095).
- Göloğlu–Granger–McGuire–Zumbrägel, CRYPTO 2013.
- Barbulescu–Gaudry–Joux–Thomé, EUROCRYPT 2014 (arXiv 1306.4244).
- Granger–Kleinjung–Zumbrägel, "Indiscreet logarithms in finite fields of small characteristic", Adv. Math. Commun. 12 (2018) 263–286 (arXiv 1604.03837).
- Kleinjung–Wesolowski, J. AMS 35 (2022).
- Joux–Odlyzko–Pierrot, in *Open Problems in Math. & Comp. Sci.*, Springer 2014, 5–36.

**Prouhet–Tarry–Escott problem.**
- Kleiman, J. reine angew. Math. 278/279 (1975) 48–51.
- Rees–Smyth, "On the constant in the Tarry–Escott problem", LNM 1415 (1990) 196–208, Prop. 2.1.
- Caley, "The Prouhet–Tarry–Escott problem for Gaussian integers", Math. Comp. 82 (2013) 1121–1137, §5 (arXiv numbering).
- Borwein–Ingalls, Enseign. Math. 40 (1994) 3–27.
- Borwein–Lisoněk–Percival, Math. Comp. 72 (2003) 2063–2070.
- Borwein, *Computational Excursions in Analysis and Number Theory*, Springer 2002, Ch. 11.
- Alpers–Tijdeman, J. Number Theory 123 (2007) 403–412.
- Filaseta–Markovich, J. Number Theory 174 (2017) 384–400.
- Croot–Mao–Yip, arXiv 2609.05061 (2026).
- Hua, "On Tarry's problem", Quart. J. Math. 9 (1938) 315–320.
- Wooley, "Some remarks on Vinogradov's mean value theorem and Tarry's problem", Monatsh. Math. 122 (1996) 265–273.

**Vinogradov's mean value theorem.**
- Bourgain–Demeter–Guth, Ann. of Math. 184 (2016) 633–682.
- Wooley, Adv. Math. 294 (2016) 532–561 (cubic case).
- Wooley, "Nested efficient congruencing and relatives of Vinogradov's mean value theorem", Proc. LMS (3) 118 (2019) 942–1016.
- Wooley, "On Vinogradov's mean value theorem, II", Michigan Math. J. (1993); the volume (believed 40) is to be confirmed.
- Ford, Proc. LMS 85 (2002) 565–633, Thm 3 (arXiv 1910.08209).
- Steiner, "Effective Vinogradov's mean value theorem via efficient boxing", J. Number Theory (2019), doi 10.1016/j.jnt.2019.04.010.

**Counting rational points (determinant method, dimension growth).**
- Castryck–Cluckers–Dittmann–Nguyen, "The dimension growth conjecture, polynomial in the degree and without logarithmic factors", Algebra Number Theory 14(8) (2020) 2261–2294.
- Binyamini–Cluckers–Novikov, "Bounds for rational points on algebraic curves, optimal in the degree, and dimension growth", arXiv 2302.04209 (2023), Thm 2, eq. (27).
- Davenport–Lewis–Schinzel, Quart. J. Math. 12 (1961).
- Fried, Illinois J. Math. 17 (1973).

**Polynomial heights.**
- Mahler, Mathematika 7 (1960).

