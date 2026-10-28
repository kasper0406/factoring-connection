# fruit: NFS fruitfulness revisited — hybrid descent, square-root randomisation, involutions (phase 30)

Agent `fruit`, 2026-10-28. **v2**: codex adversarial round (`fruit_2`, MAJOR REVISION) and the Opus proof critic (ACCEPT WITH FIXES, 3 MAJOR) applied; §7 lists every finding and its disposition.

Files:
- this note; pre-registration `notes/fruit_prereg.md` (commit `9f11355`, code frozen there);
- code `src/fruit/`: `linalg.py` (Lemma HY ranks, kernel samplers), `qs_hybrid.py` + `harness_e1.py` (E1),
  `nfs_hybrid.py` + `harness_e2.py` (E2, E3), `theory.py` (exact descent ladders, Prop IV-LV counting), `analyze.py`;
- tests `tests/test_fruit.py`;
- data `experiments/raw/fruit_e{1,2}.jsonl`, `experiments/raw/fruit_summary.json`;
- codex: `codex/fruit_1.{md,out}` (blind round), `codex/fruit_2.{md,out}` (adversarial round).

Labels: **proved** (complete proof here, elementary unless a citation is given); **proved, conditional on X**;
**model** (exact consequence of the heuristic smoothness/size model of heurnfs/lemmaM, leading order in the exponent);
**heuristic**; **open**; **experimentally supported** (toy sizes); **UNVERIFIED** (literature not checked by us).
Notation: L = L_N, D = log N/log log N, ω = ω(N), μ₂ = μ₂(Z/N) ≅ {±1}^ω.

---

## 0. Verdict

**Plain language.** We did not prove fruitfulness, and we have no rigorous L[1/3] factoring algorithm. The new tools do
not reach the missing step. What we can now say precisely is *why*, and what the cheapest sufficient hypothesis is.

1. **Hybrid descent works at the level of linear algebra, and needs only one descended row** (proved: Lemma HY and
   Theorem HYB). Add h Dixon-descended rows to Lee–Venkatesan's (LV) relation matrix. A dependency chosen from public
   data among those that use a descended row is fruitful with probability exactly 1 − 2^{1−ω} ≥ 1/2, for every N with
   ω ≥ 2, whatever the LV relations are. Such a dependency exists in the parity kernel as soon as the descended rows
   are linearly dependent modulo the row space V of the LV rows. For h = 1 this means the row lies in V
   ("absorption"). A *true* square dependency using a descended row needs r ≥ e₀ + 1, where e₀ bounds the residual
   character obstruction; a public square test then finds one (Lemma SEL). The dependency must use a descended row, and
   a public-data rule can force that without breaking Theorem DD.
2. **But one descended row already costs as much as the whole problem** (model, plus a conditional reduction).
   - In the *direct-append* architecture (one LV field; descended rows appended to its matrix), the row must be
     expressed in LV's own factor base and field: degree ≍ D^{1/3}, factor base L[1/3].
   - In the project's descent size model, a k-level descent into such a field costs L[α_k]. Pure descent, with a free
     degree and factor base, costs less at every depth. Both reach 1/3 only as k → ∞.

     | special-q sieving | hybrid (into LV field) | pure descent |
     |---|---|---|
     | linear | (k+2)/(3k+3): 1/2, 4/9, 5/12, 2/5, … | (k+1)/(3k+2): 2/5, 3/8, 4/11, … |
     | degree-t (proof critic) | 2^{k+1}/(3(2^{k+1}−1)): 4/9, 8/21, 16/45, … | 2^k/(3·2^k−1): 2/5, 4/11, 8/23, … |

     Codex, working blind, computed the linear hybrid ladder.
   - Conversely, a rigorous single descent at L[1/3, c] in a common field (Hypothesis SD₁, with the extra algorithmic
     guarantees of §1.5) would already give rigorous L[1/3] factoring by pure descent, with no LV fruitfulness at all.
     The hybrid would improve only the constant, to max((64/9)^{1/3}, c).

   So **in the direct-append architecture, the hybrid reduces NFS fruitfulness to rigorous L[1/3] descent, which is the
   project's other open problem; at the exponent level it makes neither easier.** Multi-field architectures with
   cross-field elimination are not covered by this barrier (codex fruit_2 #12).
3. **Randomising the square-root step is useless** (proved, Prop SR). The algebraic square root of a square in K is
   unique up to sign, so the class of X/Y is a function of (N, f, m, S). Any procedure that outputs a different valid Y
   with probability p is itself a factoring algorithm with success p, because ±φ(v) is computable.
4. **The single-relation involution criterion is vacuous for LV** (proved, Prop IV-LV). An exactly x-preserving
   involution pairs relations into two-term dependencies.
   - For LV's random m, the probability that *any* two non-proportional box pairs have |(a−bm)(a′−b′m)| a perfect
     square is ≤ L[2/3, −1/(2δ) + o(1)].
   - Outside this event, equal square-class vectors force proportional pairs. Those give trivial dependencies. So β is
     constant on every x-fibre, lemmaM's fibre-mixing quantity is η* = 0, and criterion (I) holds for no δ′ > 0.
   - LV's box is not primitive, which is why proportional pairs must be handled (codex).
   - Dixon's involution works because Dixon relations carry a hidden coordinate (the root). NFS relations carry none.
   - This does not exclude involutions on whole batches, on (field, relations) jointly, or modulo a subspace.
5. **Equidistribution in the product sense is false; in the Fourier sense it is Lemma M itself** (proved,
   Cor IV-TV). Since β is a function of x on X, the law of (x, β) sits on a graph. It is at total-variation distance
   ≥ min(p, 1−p) from every product law(x) ⊗ Bernoulli(p); this is non-trivial whenever β is non-degenerate. The only meaningful version is "β is far from every linear function of x",
   which is nfsfruit's Lemma M verbatim. The project's decoupling and BV tools (Prop DC, Lemma 8.5H, design DW) control
   first-order supply: smooth values in progressions, Dirichlet-character twists. Lemma M needs uniformity over the 2^rank
   elements of a binary code, i.e. 2-adic S-unit information. No tool in the repo touches that.
6. **Minimal sufficient hypotheses** (§5):
   - (F) nfsfruit's Conjecture F, which is NFS-intrinsic and untouched here;
   - or SD₁(c) + ABS′: rigorous fibre-blind L[1/3, c] descents of uniform squares into LV's factor base, of which
     (e′ + e₀ + 1)·L^{o(1)} suffice to give a descended true dependency. This is a conditional reduction (sketch), with the
     ingredients listed in §1.5, to L[1/3, max(1.923, c)]. SD₁ with those ingredients already gives rigorous L[1/3] by
     pure descent, with a worse constant.
7. **Experiments** (pre-registered; all predictions met; see §6).
   - E1: a QS base matrix that is fruitless *by construction*: 0/3600 base-only congruences split. Adding descended rows
     and forcing their use gives a fruitful rate of 989/1853 (0.534). Every pre-registered |z| was below its threshold.
   - The surprise was absorption. One descended row lay in the base row space only 6–9% of the time. The cause is
     structural: descended rows carry primes that QS relations can never contain ((N/p) = −1). Given no new column,
     absorption held in 258/265 cases.
   - E2 (cubic NFS, 80 instances, rows descended into the same field): the forced fruitful rate matched
     1 − 2^{1−ω} in all 12 cells (max |z| = 2.0). The NFS absorption rate was 80%. 0 non-square congruences out of 4589.
     The square roots via three inert primes agreed up to sign in 792/792 cases.
   - E3: 0 duplicate square-class vectors among 102 691 NFS relations (base-m and LV-style). Rational-square pairs
     involved 1.4–3.1% of relations, falling with size.

Feasibility of rigorous L[1/3] through fruitfulness: unchanged from lemmaM/nfsfruit, **low**. The hybrid is a clean
constant-level device that becomes useful only once rigorous L[1/3] descent exists.

---

## 1. Angle 1: hybrid / partial descent

### 1.1 Setting

Relations are pairs (X_i, c_i) with X_i ∈ (Z/N)^×, c_i in a certificate group 𝒢 = Q^× × K^× and
ev(c) = q·φ(γ)^{−1}, so X_i² = ev(c_i) (lemmaM §1).
- An **ordinary** NFS relation (a, b) has X_i = 1 and c_i = (a − bm, a − bα); indeed ev(c_i) = 1.
- A **descended** row j comes from r_j uniform on [1, N−1] and R_j = r_j². The procedure sees only R_j and
  independent coins. It returns W_j and c_j with ev(c_j) = R_j W_j², and stores X_j = r_j W_j.
  - `src/fruit/nfs_hybrid.py` gives the explicit form for a monic cubic, with no inversions mod N:
    X_j = r_j v_s C_u φ(P_v) and X_j² ≡ Z_j φ(Π_j). In the ev-convention the certificate is c_j = (Z_j, Π_j^{−1}),
    which is the same modulo squares.
  - The descent must gcd-check Z_j, v_s C_u and φ(P_v), so that W_j is a unit (proof critic). The toy code omits these
    checks; at toy size they would only produce factor exits.

Each relation has an *algorithmic* vector x(c) ∈ F₂^C: rational parities (LV use a constant column rather than a sign
column), ideal parities and LV's character values. Where it matters we write x_exact for the true square class in
(Q^×/Q^{×2}) × (K^×/K^{×2}). For the algorithmic vector, equal vectors give parity-kernel elements (K), not
necessarily true squares (K0). Let A be the M × C matrix of ordinary rows, D the h × C matrix of descended rows, V = rowspace(A),
K the left kernel of [A; D], and π: K → F₂^h the projection to the D-coordinates.

### 1.2 Lemma HY [proved]

rank π = r := h − (rank[A; D] − rank A) = h − dim((V + rowspace D)/V). Consequently:
- (a) a uniform element of K uses a descended row with probability exactly 1 − 2^{−r};
- (b) r ≥ 1 iff the rows of D are linearly dependent modulo V; for h = 1, iff d ∈ V;
- (c) max(0, h − (C_act − rank A)) ≤ r ≤ h, where C_act counts the coordinates charged by A or D.

*Proof.*
- ker π = {(s, 0) : sA = 0} has dimension M − rank A, and dim K = M + h − rank[A; D]. Subtracting gives rank π.
- π is a group homomorphism, so π(S) is uniform on im π when S is uniform on K. This gives (a).
- (b) and (c) are restatements, since dim((V + rowspace D)/V) ≤ min(h, C_act − rank A). ∎

### 1.3 Theorem HYB [proved; it is Theorem DD with an F-measurable cofactor]

**Lemma DD′.** In paper2 Theorem DD, replace X = ∏_{i∈S} r_i by X = G·∏_{i∈S} r_i with G an F-measurable unit. The
conclusion is unchanged: X/Y is uniform on μ₂ given F, and Pr[1 < gcd(X − Y, N) < N | F] = 1 − 2^{1−ω}.
*Proof.* In paper2's proof, X/Y = ∏ε_i · (G∏ρ(R_i)/Y), and the second factor is still F-measurable and lies in μ₂. ∎

**Theorem HYB.** Let N be odd with ω ≥ 2. Let F be generated by:
- the ordinary relations;
- the squares R_j and all descent coins;
- the character and square-root data;
- the selection coins.

Assume **transcript invariance** (T): everything the algorithm computes before forming X_S is a function of the R_j
and of coins independent of the roots, so it is F-measurable. This includes the ordinary relations, splits, children,
columns, characters, square tests, the selection of S, and the algebraic square root. (T) is a hypothesis about the
implementation; Remark T of paper2 checks it line by line for SPLIT, and the same check applies here (codex fruit_2 #2).
Suppose a true dependency S (∏_S c = g² in 𝒢) is chosen F-measurably with S ∩ D ≠ ∅, and Y = ev(g) is a unit. Then
X_S = ∏_S X_i = (∏_{j∈S∩D} r_j)·(∏_{j∈S∩D} W_j), and
**Pr[X_S ≢ ±Y (mod N) | F] = 1 − 2^{1−ω}**, for every odd N with ω ≥ 2 (including non-squarefree N such as p²q),
and for every set of ordinary relations. Prime powers are excluded by ω ≥ 2. A non-unit R_j, W_j or Y is a factor exit.
- *Non-monic f (LV).* With t = f_d, use θ = tα, its monic polynomial g(Z) = t^{d−1}f(Z/t) and g′(θ), whose image is
  t^{d−2}f′(m). Add LV's leading-coefficient/parity column. This is LV/BLP bookkeeping and is not written out here.

*Proof.*
- Every quantity except the r_j is a function of the R_j and of coins independent of the roots. This covers the
  ordinary relations, the split, the special-q children, the columns, the characters, the kernel, the selection of S,
  and the algebraic square root (Prop SR: it is a function of S).
- Apply Lemma DD′ with G = ∏ W_j, which is F-measurable and a unit after the gcd checks. ∎

**Lemma SEL (selecting a true descended dependency) [proved].** Let K0 ⊆ K be the true square dependencies of the
*fixed* augmented batch [A; D]. LV's characters may leave a residual quotient K/K0. Let e₀ be a computable upper bound
on dim(K/K0): e₀ ≤ 1 at fixed degree on paper2's character event (Lemma E: residual group of size ≤ 2), and
e₀ = O(log log N) for LV's exceptional set.
- Suppose r ≥ e₀ + 1. Choose a subspace U ⊆ K of dimension e₀ + 1 on which π is injective; this is public linear
  algebra.
- Then dim(U ∩ K0) ≥ 1. Every nonzero T ∈ U ∩ K0 has π(T) ≠ 0, i.e. uses a descended row.
- The algorithm finds such a T by exact square tests on the 2^{e₀+1} elements of U. These tests are public (they act on
  certificates only), so Theorem HYB applies to the T found.
- Cost: 2^{e₀+1} = (log N)^{O(1)} square-root tests.

Codex fruit_2 #3 notes that LV's own procedure, fresh batches plus pairing, does not literally produce a uniform
element of a fixed K0. Lemma SEL replaces it for the hybrid: keep one fixed augmented batch.
- Uniform sampling, for comparison: a uniform T ∈ K0 uses a descended row with probability 1 − 2^{−r₀}, where
  r₀ = rank(π|K0) ≥ r − dim(K/K0).

So **one descended row (h = 1) suffices** when it is absorbed and r₀ ≥ 1, i.e. when the residual functional is not π
itself. At fixed degree h ≤ 2 absorbed rows suffice. In general h ≥ (#rows needed for r ≥ e₀ + 1). Lemma HY(c) bounds
this deterministically by (C − rank A) + e₀ + 1 with C the *total* column count; the support-based count is circular
when every new row brings new columns. **h ≥ 1 is necessary** for any argument by hidden-root symmetry
(Prop T, §1.6).

### 1.4 Absorption

**Lemma ABS-dom [proved].** Let the ordinary rows be M i.i.d. samples from a law μ on F₂^C. This is nfsfruit Prop R:
LV's kept relations of a fixed field are i.i.d. from the kept law. Let V_M be their span.
- Missing mass: E[μ(V_M^c)] ≤ C/M.
- If a descended row is d = y₁ + … + y_k with each y_i independent of the ordinary rows and of marginal law ≤ K·μ, then
  Pr[d ∉ V_M] ≤ kKC/M.

*Proof.*
- Put q_j = μ(V_j^c), which is nonincreasing. dim V_{j+1} = dim V_j + 1 exactly when x_{j+1} ∉ V_j, an event of
  conditional probability q_j. So Σ_{j<M} E q_j = E dim V_M ≤ C, and M·E q_{M−1} ≤ C.
- Then Pr[y_i ∉ V_M] = E[law(y_i)(V_M^c)] ≤ K·E q_M ≤ KC/M. Take a union bound over i. ∎

**What blocks domination [experimentally supported, plus an elementary reason].** A descended row charges coordinates
that ordinary relations charge rarely or never.
- In QS, the split halves u, v contain small primes with (N/p) = −1, which divide no x² − N.
- In NFS, the children are special-q points, not box points. Their laws are not dominated by μ.

E1 shows the effect is decisive at toy size. One descended row was absorbed in only 6–9% of cases, and a new column
was present in essentially every failure: given no new column, absorption held in 258/265 cases. E2 measures the NFS
version (§6).

**Lemma ABS-red (multi-row absorption from support plus reachable deficit) [proved; corrected in v2].** Work in
coordinates.
- Let 𝒮 = supp(A) be the coordinates charged by the batch.
- Let Γ_sq = ev^{−1}((Z/N)^{×2}) inside the S-unit certificate group Γ (S-units of Q and of K for the factor base).
- Let R_𝒮 := x(Γ_sq) ∩ F₂^𝒮 be the *reachable* space: descended rows always lie in x(Γ_sq).
- Let e′ := dim R_𝒮 − dim(V ∩ R_𝒮) be the batch's deficit inside the reachable space.

(a) Call a descended row *support-compatible* if supp(d) ⊆ 𝒮. This is a public test, so rejecting incompatible rows
keeps Theorem HYB valid. Then **e′ + 1 support-compatible rows always give r ≥ 1**, since they live in R_𝒮, where V
has codimension e′.
- This is a multi-row statement (ABS′). It does *not* give single-row absorption: codex's example V = ⟨(1,1)⟩,
  d = (1,0) has both coordinates charged and d ∉ V.
- v1 used e = |𝒮| − rank A. The proof critic showed that this deficit is structurally large for LV. LV's ≫ log N
  character columns are charged by essentially every row, but rank A ≤ dim x(Γ) ≤ #parity columns + O(d + log N).
  So e ≥ #characters − O(d + log N). Hence v1's (ABS2) was unsatisfiable, and e′ replaces it. E2's median e = 44 is
  this effect.
(b) Here is what full span would mean structurally. Let Γ be the certificate group of factor-base S-units, in Q^× and
  in K^×, and Γ₁ = {c ∈ Γ : ev(c) = 1}. Ordinary relations lie in Γ₁. A descended row
  has ev(c_R) = R W_R². If some g ∈ Γ has ev(g) ∈ r W_R·μ₂, which holds whenever ev(Γ) = (Z/N)^×, then c_R g^{−2} ∈ Γ₁.
  Hence x(c_R) ∈ x(Γ₁), because x kills squares. So **a descended row fails to be absorbed only if the batch does not
  span x(Γ₁) on its support**, i.e. only if box relations fail to generate the mod-2 relation module.
  - This generation statement (G2) is the NFS analogue of the generation lemmas proved, under GRH, for
    *random-exponent* relations in class-group algorithms (Seysen 1987; Hafner–McCurley 1989). Box relations are not
    random-exponent relations, and we know no proof.
  - Condition ev(Γ) = (Z/N)^× asks the factor base to generate (Z/N)^×. Unconditionally this is not known for
    B = L[1/3]: it is a least-character-nonresidue question, as in lemmaM's B1 barrier. It is only sufficient, not
    necessary.
(c) QS shows that (G2) can fail structurally. QS values x² − N are divisible only by primes with (N/p) = 1, so the
  coordinates of non-residue primes are never charged. A descended row whose split halves contain such primes is
  support-incompatible. This is exactly the E1 failure mode.

**Hypothesis ABS′.** For a batch of C·L^{o(1)} LV relations, and under the joint law of successful descents and
batches:
- (ABS1) a descended row is support-compatible with probability ≥ L^{−o(1)}. This is a divisibility statement of the
  kind LV's AP-distribution lemmas give.
- (ABS2′) e′ ≤ (log N)^{O(1)} with probability ≥ 1/2. This is a mod-2 generation statement: box relations span the
  reachable space, cf. (b).

By (a) and Lemma SEL, (ABS1) + (ABS2′) give a descended true dependency from (e′ + e₀ + 1)·L^{o(1)} descents.
Single-row ABS (Pr[d ∈ V] ≥ L^{−o(1)}) is a separate, stronger hypothesis.
- ABS′ concerns x only, never β. But rank and generation are not "first-order" in the analytic sense (codex fruit_2 #4).
  They are global linear-algebra properties of the relation module.
- We do not prove ABS′.
- The triage `codex/r29_triage.out` (lead B) independently named absorption as the first lemma, with the same warning
  that nullity alone gives nothing.
- E2 (cubic NFS) is encouraging: a single descended row was absorbed in 638/798 cases (80%), and introduced only 0.19
  new columns on average, although the naive deficit e had median 44 (character columns).
- **Absorption is not the binding obstacle** (§1.5): even a perfectly absorbed single row must first be produced in
  L[1/3].

### 1.5 What one descended row costs [model; codex fruit_1 agrees blind on the linear ladder]

**Architecture.** In the *direct-append* architecture, one LV field is used and the descended rows are appended to its
matrix. Then the row must live in LV's columns, so its last-level children must be smooth over LV's factor base L[1/3]
in LV's field, of degree ≍ D^{1/3}.
- A second field brings its own ideal columns, which one row cannot cancel. Codex fruit_2 #12 points out that
  multi-field architectures could still cancel them: complete certificates in a second field, eliminated blocks, or
  MNFS-style shared rational columns. No argument here charges cross-field elimination, so the barrier below is stated
  **only for direct append**.
- Two rational sides force a √N-size coefficient (nfsfruit L1 / Prop B).

**Size model** (lemmaM §2 / heurnfs §3.2):
- the split needs 1 − β₁ ≤ α;
- children that are L[β_{i+1}]-smooth cost L[n(β_i) − β_{i+1}], where n(β) is the norm exponent of a special-q step
  from primes L[β].
  - *Linear* special-q sieving (a − bα), as in lemmaM/lv5unif: n(β) = γ + β, where γ is the degree exponent.
  - *Degree-t* sieving in the special-q lattice, the usual individual-log practice (proof critic): t ≍ D^{(β−1/3)/2}
    balances both sides at n(β) = (1+β)/2.
- Recognition is by recursive self-calls; batch smoothness cannot test L[β₁] primes with β₁ > α. So HECM's L[4/9]
  objection does not apply (codex triage). The recursive version needs a Lemma T_rec analogue (paper2 §recursive).

| sieving | hybrid (γ = φ = 1/3) | pure (γ, factor base φ = α free) |
|---|---|---|
| linear | (k+2)/(3k+3): 1/2, 4/9, 5/12, 2/5, 7/18, … | (k+1)/(3k+2): 2/5, 3/8, 4/11, … (= multidescent/pc3) |
| degree-t | 2^{k+1}/(3(2^{k+1}−1)): 4/9, 8/21, 16/45, … | 2^k/(3·2^k−1): 2/5, 4/11, 8/23, … |

- All four ladders are checked exactly in `theory.chain_feasible{,_t}` (tests).
- Codex's blind table reproduces the linear hybrid ladder.
- In both sieving models, the hybrid at depth k is worse than pure descent at depth k, and every finite depth stays
  above 1/3.
- These are model optima, not algorithmic lower bounds. In particular, k → ∞ does not by itself give a uniform
  L[1/3, c] algorithm.

**Hypothesis SD₁(c)** (strengthened after codex fruit_2 #13). There is a field law 𝒜, samplable in time
L[1/3, c_A], with the following property. With probability ≥ L^{−o(1)} over A ∼ 𝒜, a procedure that sees only R,
uniform in (Z/N)^{×2}, and coins independent of its root returns, with probability ≥ L^{−o(1)} and in time L[1/3, c], a
certificate of R over the factor base of A that is:
- sparse: (log N)^{O(1)} nonzero coordinates;
- of height L^{O(1)};
- passed by the characters' admissibility conditions.

The procedure can be repeated with fresh coins.

**Proposition HX [conditional reduction, sketch; model for (iii)].**
- (i) SD₁(c) implies rigorous factoring of every N with ω ≥ 2 in time L[1/3, max(c_A, β + c, 2β)] by pure descent:
  - retain the first field with a success (paper2 SPLIT, Prop select);
  - take C + e₀ + 1 descended rows;
  - characters and square roots for multi-level atoms at growing degree (paper2 LV6 / lv5unif §9 cover only one-level
    atoms, so their heights must be rechecked);
  - sparse linear algebra, Lemma SEL, Theorem HYB.

  No LV fruitfulness is involved. The items listed are proof obligations, not written out here.
- (ii) In a common field, SD₁(c) + ABS′ give the hybrid in L[1/3, max((64/9)^{1/3}, c_A, c) + o(1)]. When 𝒜 is LV's
  search law, c_A ≤ (64/9)^{1/3}.
- (iii) In the direct-append architecture, at every finite depth, the hybrid's exponent exceeds pure descent's (model).

So in the direct-append architecture the hybrid is a constant-level device. Fruitfulness via appended hidden roots is,
at the exponent level, no easier than rigorous L[1/3] by pure descent.
- Heuristically c ≈ 3^{1/3} ≈ 1.44 (Commeine–Semaev's individual-log constant; heurnfs §3.1). Barbulescu's early-abort
  improvements lower it [UNVERIFIED transfer]. Then the maximum is 1.923, and the hybrid would attain LV's constant.

**Verdict on angle 1:**
- *How many rows:* one, if absorbed and r₀ ≥ 1; at fixed degree at most two absorbed rows. In general, enough rows for
  r ≥ e₀ + 1, which is at most C − rank A + e₀ + 1.
- *Must the dependency use them:* yes. Lemma SEL forces it from public data.
- *Barrier:* the cost of producing even one row. At L[2/5] one row already exceeds an L[1/3] budget, as the brief
  anticipated, and the direct-append ladders are worse still.

**Refuted for this model (direct append) as a route to L[1/3] by itself. It is a proved linear-algebra/symmetry
reduction of fruitfulness to (SD₁ + ABS′).**

### 1.6 Proposition T: transcript-preserving symmetry needs coins outside the transcript [proved, formal; scope narrowed in v2]

Let the congruence (X, Y) be computed from coins ξ, and let τ(ξ) be a chosen exposed transcript.
- A *pointwise symmetry argument for τ* is a measure-preserving T on coin space with τ∘T = τ and
  (X/Y)∘T = ε·(X/Y), for some ε ∈ μ₂ ∖ {±1}.
- If X/Y is τ-measurable, no such T exists. This is a tautology.
- Its content is in the choice of τ. For LV with τ = (m, f, pairs, characters, kernel coins), Prop SR makes X/Y a
  function of τ. So every symmetry argument for LV must either use coins that LV does not have (hidden roots), or
  preserve a *coarser* transcript.

Hidden roots can enter in several ways:
- by appended descended rows (this note);
- through a shared hidden root at polynomial selection (nfsfruit Theorem L1; Prop B is a counting obstruction for
  globally uniform roots, not a general barrier);
- through other architectures not classified here.

Arguments preserving a coarser transcript are not excluded: transformations of whole batches, of field parameters
preserving the dependency structure, and coordinate-permuting actions (codex fruit_2 #8). v1's sentence "a symmetry
proof needs … a Dixon or descended row" is withdrawn as an overclaim.

## 2. Angle 2: randomising the square-root step — Prop SR [proved; refuted as a route]

Let f be monic (for LV's non-monic f use θ = f_dα as in §1.3), and let S ∈ K0 with γ_S = f′(α)²∏_S(a_i − b_iα) ∈ K^{×2}.
- (a) The square roots of γ_S in K are exactly ±v, since (w − v)(w + v) = 0 in a field; odd degree is not needed for
  this. And v ∈ Z[α], because f′(α)O_K ⊆ Z[α]. So [u/φ(v)] ∈ μ₂/{±1} is a function
  of (N, f, m, S).
- (b) Every exact method returns ±v, hence ±φ(v):
  - Couveignes (odd degree, sign fixed by the norm);
  - Montgomery's lattice method;
  - Hensel lifting at any inert prime;
  - the choice of archimedean embedding, which only bounds coefficients.
  Choosing signs at different primes inconsistently gives a CRT lift w′ with w′² ≡ γ_S only modulo the auxiliary primes.
  Its image φ(w′) is not a square root of φ(γ_S) mod N except by accident, and then it falls under (d).
- (c) A unit twist ηv squares to η²γ_S ≠ γ_S unless η = ±1. Multiplying the product by a known square η² is a
  per-relation square multiplier, which leaves Φ unchanged (nfsfruit codex §3.1). For even degree, roots of unity
  ζ ∈ μ_K give valid roots only when ζ² = 1. Other ring maps α ↦ m′ (other known roots of f mod N) break the
  rational-side congruence. They give MNFS-type congruences, each again deterministic.
- (d) **Any** procedure that outputs Y with Y² ≡ φ(γ_S) (mod N) and Y ∉ {±φ(v)} with probability p factors N with
  probability p, at the extra cost of computing φ(v): if gcd(Y − φ(v), N) = 1, equality of squares forces
  Y = −φ(v). This is the classical Rabin reduction (MIT/LCS/TR-212, 1979). So randomness in the square-root step can help only if that
  step is already a factoring algorithm.
- E2 Q5 confirms (a) on every checked dependency: the roots via three different inert primes agree up to sign.

Theorem DD's hypothesis "square root computed from the transcript" is about *invariance*, not randomness, and (a)
supplies it for free. **Angle 2: refuted (proved).**

## 3. Angle 3: the involution criterion

**Prop IV [proved].** Let x = x_exact, so that equal vectors give true square dependencies. Suppose ι: X → X preserves
x and flips β on a set F ⊆ X. Then every z ∈ F lies in the two-term dependency {z, ι(z)} ∈ K0, with
Φ({z, ιz}) = β(z) + β(ιz) = 1. So there are |F|/2 fruitful two-term dependencies (under counting measure; a weighted
version holds in general), and each factors N with one square root. With the algorithmic x (finitely many characters)
one only gets {z, ιz} ∈ K.

**Prop IV-LV [proved; v2 handles LV's non-primitive box].**
- *Setting.*
  - E = L[1/3, σ].
  - LV's box (Definition 4.1) is 𝔅 = {(a, b) ∈ Z² : 0 ≤ a < |b|, |b| ∈ [E/2, E]}, with *no* coprimality condition.
  - m is uniform on the integers of LV's interval I = [(N/2)^{1/d}, N^{1/d}]. Put L_m := N^{1/d} = L[2/3, 1/δ + o(1)],
    so |I| ≥ L_m(1 − 2^{−1/d}) − 1 ≥ L_m/(3d).
- *Claim.* Pr_m[∃ non-proportional z, z′ ∈ 𝔅 with |(a − bm)(a′ − b′m)| a nonzero perfect square]
  ≤ 12·E⁴·τ*·√(2E L_m)·d/L_m = L[2/3, −1/(2δ) + o(1)], where τ* = max_{0<|Δ|≤2E²} 2^{ω(Δ)} = E^{o(1)}.
- *Search law.* Let the search try fields until it stops. If T is the number of fields tried, Wald's identity gives
  Pr[output m ∈ 𝒜] ≤ E[T]·Pr_unif[𝒜]. For LV's stochastic deepening E[T] ≤ L[1/3]^{O(1)}, so the bound survives
  (proof critic's fix).

*Proof.*
- For non-proportional pairs, Δ = ab′ − a′b ≠ 0.
- Suppose the absolute product is a nonzero square. Then a − bm and a′ − b′m have squarefree parts equal up to sign.
  Write a − bm = d u² with d squarefree and signed. Then d divides both values, hence divides
  b′(a − bm) − b(a′ − b′m) = Δ.
- u ≤ √((|a| + |b|L_m)/|d|) ≤ √(2E L_m), and each (d, u) determines at most one m.
- So at most 2·2^{ω(Δ)}·√(2E L_m) values of m work for a given pair. Divide by |I| and sum over at most E⁴ pairs. ∎
- `theory.rational_square_bound` implements the per-pair bound; a test checks it against exact counts.

**Proportional pairs [proved].** Let z′ = s·z as pairs, with s ∈ Q^×. Then a′ − b′m = s(a − bm) and
a′ − b′α = s(a − bα). If x_exact(z) = x_exact(z′), then s is a square in K, hence in Q because d is odd: s^{d−1} is a
square and s^d = N_{K/Q}(t)². So s = t² with t ∈ Q. The pair's roots are U = t(a − bm) and V = t(a − bα), with
φ(V) = U. The dependency {z, z′} is trivial, so β(z) = β(z′). The example (9−18m)(16−32m) = 144(1−2m)² from codex shows
such pairs exist in LV's literal box.

**Consequences [proved].** Off the event of the Claim:
- (i) x_exact(z) = x_exact(z′) implies that z and z′ are proportional and β(z) = β(z′). So β is constant on every
  x-fibre, and lemmaM's fibre-mixing quantity η* = Σ_x min(Pr[x, β=0], Pr[x, β=1]) is 0.
- (ii) No exactly x-preserving map flips β anywhere: criterion (I) holds for no δ′ > 0. On the primitive sub-box, x is
  even injective, so the only x-preserving map is the identity.
- (iii) This excludes only the single-relation, exactly-x-preserving criterion. Involutions on batches, on
  (field, relations) jointly, coordinate-permuting actions, or actions preserving x modulo a subspace are not
  excluded; they would need their own compatibility proofs (codex fruit_2 #6).
- (iv) Cor IV-TV, §4.
- η* = 0 does not refute fruitfulness (codex fruit_1 §13): what dies is the involution *route* to Lemma M.

**The candidates, individually** (each also falls under (ii) when it acts on one field's X):
- *N·w multipliers* (paper2, lv5unif) and *design DW's integer multiplier* change the field. They map relations of one
  field to relations of another, so x lives in different coordinates. These are field-level designs for supply moments,
  not maps on X.
- *Galois automorphisms.* An odd-degree K with S_d Galois group has none. A Möbius symmetry fixing m preserves the root
  ratio for even |S|, or changes it by an explicit constant, which by codex's square-twist identity would factor N
  (lemmaM §5).
- *(a, b) ↦ (−a, b).* This moves a − bm to a different value, so x changes entirely. (a, b) ↦ (−a, −b) shifts x by the
  fixed sign/(−1) vector and leaves Φ unchanged on K0, because (−1)^{|T|} with |T| even.
- *The unit twist of QS/CFRAC.* It is the one genuine example. It displaces by a full period, about N^{1/2}/h (lemmaM E2).

**Why Dixon escapes.** Dixon's involution r ↦ εr fixes the *visible* relation, because R = r² is
unchanged, and moves only the hidden root. NFS relations have no hidden coordinate. A visible-data, exactly
x-preserving involution would need non-trivial x-collisions, and Prop IV-LV shows there are none: the only collisions
are proportional pairs, and those are trivial. **Angle 3: refuted for LV's single-relation criterion (proved).**
Adding hidden coordinates (angle 1) is one way around this. Batch-level or field-level symmetries are not excluded.

## 4. Angle 4: equidistribution of (x, β)

**Cor IV-TV [proved].** Suppose β is a function of x on X. By Prop IV-LV(i) this holds for LV off an event of
probability L[2/3, −1/(2δ)+o(1)]. Then the joint law of (x, β) under any μ on X is supported on a graph, and
TV(law(x, β), law(x) ⊗ Bern(p)) = Σ_x μ(x)(1 − p_{β(x)}), where p₁ = p and p₀ = 1 − p. Minimised over p this is exactly
min(μ(β=0), μ(β=1)). So "close to a product" in total variation fails whenever β is non-degenerate. Codex fruit_1 §14
gives the same identity in the form 1/2 − η*.

**The Fourier version is Lemma M [identity].** The weakest sensible version is max_w E(−1)^{β + ⟨w,x⟩} ≤ 1 − 2δ.
The absolute value may be added when x includes the constant feature, which makes the code closed under
complementation.
- By duality, the restrictions to X of the tests (−1)^{⟨w,x⟩} are exactly the ±1 functions on X that are trivial on
  K0(X).
- So the test family is a binary code of dimension rank(x(X)), and "β is far from the code" *is* nfsfruit's Lemma M
  (conn_alg C8). Product-likeness adds nothing.

**What the project's tools control.**
- *Prop DC* (pc3) bounds dead mass in descent trees by relative variances of child counts.
- *Lemma 8.5H* (uprime) is a restricted BV for B-smooth modulus cofactors with B₂-smooth values.
- *Design DW* (c12alg) turns algebraic two-point blocks into lattice problems.

All three are first-order statements about how smooth values are distributed in residue classes. The twists they can
average over are Dirichlet characters to moduli ≤ Q: a structured family of size about Q². Lemma M needs uniformity
over 2^{rank} arbitrary prime-sign twists. On the |X|-point set these include *every* function trivial on K0. And β
itself is the 2^{c+1}-adic Bockstein lift (nfsfruit Prop G), which carries S-unit data modulo 4. Codex Theorem B
(lemmaM §4.3) already shows that fixed-test averaging is compatible with δ = 0 for every N.

**Averaged product laws prove nothing [proved; codex fruit_1 §16].** Take β = ⟨θ, x⟩ with θ ∈ F₂^C hidden and
uniform, and x ≠ 0 almost surely (e.g. constant feature present). Every instance is fruitless (Φ ≡ 0), yet averaged
over θ the law of (x, β) is exactly law(x) ⊗ Bern(½).
So equidistribution statements about *marginal* laws, averaged over N or over field parameters, cannot give
fruitfulness without per-instance control. Averaging a suitable nonlinear kernel-success quantity over fields is not
excluded (codex fruit_2 #10). This is the information-theoretic form of codex Theorem B (lemmaM). Conversely η* = 0 does not refute
fruitfulness (codex fruit_1 §13): Prop IV-LV kills the *involution route*, not Lemma M.

**One place a product law does hold exactly:** descended rows. There x(d) is a function of R, and the root sign ε is
uniform and independent of it given everything (Theorem DD). The "decoupling" that angle 4 asks for is precisely the
hidden-root fibre.

[assessment] The natural language for Φ ≢ 0 is 2-Selmer/Cohen–Lenstra-type statistics of S-unit groups in random
fields. Proven results of that type (Bhargava-style averages of 2-torsion in class groups of odd-degree fields) are
for *fixed* degree [UNVERIFIED details]. LV's degree grows like D^{1/3}, and the needed statement concerns a specific
kernel character, not a class group. **Angle 4: refuted for this model in the TV sense (proved); the Fourier sense is
Lemma M, open.**

## 5. Barrier and minimal hypotheses

**What is proved about the routes** (not an exhaustive classification of proofs; v1's "must take one of two forms" is
withdrawn):
- **Symmetry inside LV's transcript is impossible.** Prop SR shows the square root carries no entropy. Prop IV-LV shows
  there is no exactly x-preserving β-flip and η* = 0. Prop T: no transcript-preserving symmetry exists without coins
  outside LV's transcript.
- **Distributional routes** (Conjecture F, or Lemma M via Prop R) need a property of 2-adic S-unit data that is
  uniform over a binary code of dimension ≍ L[1/3]. Averaged marginal product laws cannot give it (§4).
- **Hidden-root routes in the direct-append architecture** need at least one descended row in LV's factor base. Its
  cost in the size model is ≥ L[α_k] at depth k (§1.5 table), and a rigorous L[1/3] version (SD₁) would already give
  rigorous L[1/3] by pure descent (Prop HX).
- Not covered: multi-field hidden-root architectures, hidden roots at polynomial selection, batch-level or field-level
  symmetries [assessment: we found no candidate in any of these].

**Minimal sufficient hypotheses we can state** (each gives rigorous L[1/3] for every N with ω ≥ 2; Lemma H /
Theorem DD remove any Blum condition):
- **(F)** nfsfruit Conjecture F. With LV Thm 2.1 and Lemma H it gives L[1/3, (64/9)^{1/3} + o(1)], modulo LV's
  non-monic bookkeeping and HECM reliance.
- **(SD₁ + ABS′)**, for the hybrid (§1.4–1.5). This is a conditional reduction (sketch): L[1/3, max((64/9)^{1/3}, c)]
  modulo LV's bookkeeping and the character and square-root adaptation for multi-level descended rows.
  - SD₁ alone, with the same adaptations, gives L[1/3] by pure descent.
  - SD₁ is a multi-level, single-target version of the project's relation-supply theorems. The project has one level
    at L[2/5] (lv5unif, its own field) and, conditionally, two levels at L[3/8] (multidescent/pc3/c12alg).
  - The unbounded depth needed for 1/3 is the open part. Codex's triage notes that precision at depth i must reach
    D^{−iα}, while current estimates stop at the D^{−2α} scale.

## 6. Experiments (pre-registered, `9f11355`; exact arithmetic)

Commands:
```
PYTHONPATH=. .venv/bin/python -m src.fruit.harness_e1 experiments/configs/fruit_e1.json experiments/raw/fruit_e1.jsonl
PYTHONPATH=. .venv/bin/python -m src.fruit.harness_e2 experiments/configs/fruit_e2.json experiments/raw/fruit_e2.jsonl
PYTHONPATH=. .venv/bin/python -m src.fruit.analyze experiments/raw/fruit_e1.jsonl experiments/raw/fruit_e2.jsonl experiments/raw/fruit_summary.json
```

**E1** (QS base with β ≡ 0 by construction, 60 Blum instances, 48/64 bits, 3 base sizes × 20 rounds × h ∈ {1,2,4}).

| prediction | result |
|---|---|
| P1 base-only fruitful = 0 | 0 / 3600 |
| P2 forced = ½, cells (n ≥ 20) \|z\| ≤ 3.5 | 9 cells, z ∈ [−0.79, 2.86]; pooled h=1 147/258 (z 2.24), h=2 274/527 (0.91), h=4 568/1068 (2.08): **met** |
| P3 uniform vs Σ p_i | 518 vs 490.8 (z 1.44): **met** |
| P4 identities | 0 violations |
| P5 absorption (h = 1) increasing in base size | C+5: 74/1200 = 0.062 < 2C: 0.068 < 4C: 0.086: **met (weakly)** |

- The pooled forced rate is 989/1853 = 0.534. All z values are positive in the C+5 cells and near 0 elsewhere. The
  pre-registered thresholds hold. Theorem DD is exact, so we read this as chance: there are 9 correlated cells.
- *Post hoc (not pre-registered).* A descended row adds on average 2.7–3.2 coordinates that the base never charges.
  These are small primes with (N/p) = −1 in u, v. Given *no* new coordinate, a single row was absorbed in 258/265 cases
  (74/81, 81/81, 103/103 by level). The batch deficit on its own support had median 1–4. So in QS the whole absorption
  failure is the support mismatch of Lemma ABS-red(c).

**E2** (monic cubic NFS base, 80 instances: blum/mixed/both1/tri × 40/44 bits × 10 seeds; 10 rounds × h ∈ {1,2,3}).

| prediction | result |
|---|---|
| Q1 forced = 1 − 2^{1−ω}, cells \|z\| ≤ 3.5 | 12 cells, z ∈ [−2.00, 1.63]; pooled h=1: semiprimes 226/480 (z −1.28), tri 117/158 (z −0.28): **met** |
| Q2 non-square ≤ 5% | 0 / 4589 |
| Q3 uniform ≥ lower bound | 1323 vs 818.4 (the base relations are themselves fruitful: 464/798 base-only samples): **met** |
| Q4 identities | 0 violations |
| Q5 square roots via 3 inert primes agree up to sign | 792 / 792 |
| Q6 (exploratory) absorption, h = 1 | 40 bits 335/399 = 0.84; 44 bits 303/399 = 0.76; Pr[r ≥ 1] = 0.80 / 0.96 / 0.99 for h = 1/2/3; 0.19 new coordinates per row; batch deficit on support median 44 |

- Two descents hit an ACCIDENTAL factor: a gcd check on u or v met a 13–14-bit prime factor of a 3-prime N (6397,
  16217). This matches the null rate ≈ 2/p per split; the two rounds were dropped.
- Total run time was 355 s, with peak RSS well under 1 GB.

**E3** (x-injectivity).

| prediction | result |
|---|---|
| R1 identical full parity vectors: 0 pairs | 0 pairs, base-m (56 017 relations) and LV-style (46 674 relations) |
| R2 relations in a rational-square pair ≤ 5% | base-m 3.06% (40 b) → 2.10% (44 b); LV-style 2.21% → 1.42%: **met**, decreasing |

At toy size the box is far from Prop IV-LV's regime: the bound is vacuous when E⁴√(EL)d/L ≫ 1. R2 is reported only
as a trend.

All predictions were met. These are correctness checks of Lemma HY, Theorem HYB and Prop SR at toy size, plus
descriptive absorption data. They carry no complexity claim.

## 7. Reviews

**GPT-6.1-sol blind round (`codex/fruit_1.out`, xhigh, launched before any analysis here).** It reached the same
verdict independently.
- Its §2–4 give the hybrid threshold: "one row can suffice; sharp worst case = deficit + 1", forcing by public data,
  and an allowance for the exceptional set. This matches Lemma HY and Theorem HYB.
- Its §7 gives the hybrid ladder min_a max{1 − a, a/2, γ_k(a)} = (k+2)/(3(k+1)) = 1/2, 4/9, 5/12, 2/5, matching §1.5.
  It adds that a *raw* polynomial-size target has an L[4/9] floor at the first norm, while presmoothed targets do not.
- Its §8: the hybrid improves constants, never the exponent over pure descent.
- Its §9–10: square-root randomness is useless, as in Prop SR.
- Its §11: x(z) = x(z′) forces simultaneous square twists, so an efficiently constructible flipping involution
  factors N. This is a weaker form of Prop IV; it did not find the m-averaged injectivity Prop IV-LV.
- Its §13 and §16 add two points, merged into §4:
  - η* = 0 does not refute fruitfulness. Its example has four vectors with singleton fibres whose sum has label 1.
  - Countermodel: β = ⟨θ, x⟩ with θ hidden and uniform. Every instance is fruitless, yet the θ-averaged joint law is an
    exact product. So averaged product laws cannot imply fruitfulness.
- Its §20 cites literature: Guillevic 2019 (heuristic smoothing), Al Aswad–Pierrot 2023, and Kleinjung–Wesolowski
  (rigorous quasi-polynomial in fixed characteristic, a different geometry). None gives rigorous descent in an LV field.
- Its exact checks covered 44 494 matrices and all unit-square fibres modulo 15, 21, 45, 105 (codex scratch, not
  committed).

**GPT-6.1-sol adversarial round (`codex/fruit_2.out`): MAJOR REVISION.** Each finding and what we did:

| # | finding (severity) | disposition |
|---|---|---|
| 1 | HY, DD′ correct (MINOR) | kept |
| 2 | HYB: state transcript invariance as a hypothesis; non-monic normalisation; gcd exits (MAJOR) | applied (§1.3) |
| 3 | "Choosing S": LV's fresh-batch pairing does not give uniform T ∈ K0; e not public; K = K0 misattributed (MAJOR) | replaced by Lemma SEL (fixed augmented batch, computable e₀, public square tests) |
| 4 | Support + deficit ⇏ single-row absorption; "first-order" misleading (MAJOR) | ABS-red restated as multi-row (ABS′); single-row ABS kept separate; label fixed |
| 5 | IV-LV: LV's box is not primitive, so the literal event has probability 1 (FATAL to the literal application) | proportional pairs handled: they are trivial dependencies, so β is constant on fibres and η* = 0 still; exact m-interval; no sign column |
| 6 | IV-LV consequences: narrower scope (batch, field and subspace actions not excluded) (MAJOR) | applied (§3, §0) |
| 7 | Prop IV needs x_exact; counting measure (MINOR) | applied |
| 8 | Prop T is a tautology; the exhaustive dichotomy is wrong (MAJOR) | narrowed; §5 dichotomy withdrawn |
| 9 | SR: odd degree unnecessary; non-monic; Rabin cost (MINOR) | applied |
| 10 | TV needs β non-degenerate; abs-value Fourier needs constant feature; countermodel needs x ≠ 0 (MINOR) | applied |
| 11 | Ladders correct within the model; not lower bounds (MINOR) | stated |
| 12 | "Row must live in LV's field" false in general (MAJOR) | restricted to direct append |
| 13 | HX/SD₁: residual classes, sparsity, heights, field retention, search cost (MAJOR) | SD₁ strengthened; HX relabelled a conditional reduction (sketch) with its obligations listed |
| 14 | Prior art: Dixon, Pomerance, BLP (Lemma 8.2, obstruction groups), Couveignes, Montgomery; Seysen p. 777 UNVERIFIED; Gordon's 2/5 attribution UNVERIFIED | recorded in §7.1 |
| 15 | Verdict: replace "only way", "must", "exactly as hard" | applied |

Its independent checks confirmed 2324 primitive-pair bounds and both linear thresholds through depth 100 (codex
scratch, not committed).

**Opus proof critic (read-only): ACCEPT WITH FIXES, 3 MAJOR.**
- MAJOR 1: (ABS2) is unsatisfiable for LV, because LV's character columns lie in the support but not in reachable rank.
  Applied: the reachable deficit e′ (§1.4).
- MAJOR 2: HX(i) was under-specified (same-field retention, residual rows, multi-level characters and heights,
  recursive recognition, field-search cost). Applied (§1.5).
- MAJOR 3: the ladders are specific to linear sieving. Degree-t special-q sieving gives 4/9, 8/21, 16/45 (hybrid) and
  2/5, 4/11, 8/23 (pure). Applied: both models in §0 and §1.5, checked in `theory.chain_feasible_t`. The qualitative
  conclusions (finite depth > 1/3; hybrid worse than pure at equal depth) hold in both.
- MINOR, all applied:
  - the units W_j need gcd checks;
  - the ev sign convention, c_j = (Z_j, Π_j^{−1});
  - paper2 Lemma E gives |E| ≤ 1, not K = K0, so h ≤ 2 at fixed degree;
  - the circular C_act bound;
  - x_exact versus the algorithmic x;
  - the Wald form of the search-law extension;
  - renaming L_m;
  - TV non-degeneracy;
  - x ≠ 0 in the countermodel;
  - Prop T, "injecting is descent" and the §5 dichotomy relabelled as assessment.
- The critic hand-verified the nfs_hybrid algebra (X_j, Z_j, Π_j, final X, Y) and found it correct.

### 7.1 Opus literature critic (orchestrator-run, 2026-10-28): OK WITH FIXES

No theorem (HY, DD′, HYB, SEL, SR, IV-LV, IV-TV) and no verdict changes.

**Sources read:**
- LV arXiv 1805.08873v1 source, so the numbering is arXiv v1;
- LNM 1554 (BLP; Couveignes);
- Guillevic arXiv 1505.07553;
- EGT arXiv 0905.2177;
- Stange arXiv 2211.06821;
- the MIT record of Rabin's TR-212.

The Dixon PDF was blocked, so Dixon is checked only via secondary sources.

**Lee–Venkatesan.**
- Thm 2.1 gives a congruence x² ≡ y² (mod n), possibly trivial, in L_n[1/3, (64/9)^{1/3} + o(1)].
- **Thm 2.3 gives fruitfulness only for semiprimes with both factors ≡ 3 (mod 4), conditional on their Conjecture 7.1.** This note's results do not need either restriction.
- Thm 2.5 is relation supply. Thm 2.6 is the character reduction, valid for all but an L(2/3, ·)^{−1} fraction of f.
- Def 4.1 matches §1's box: 0 ≤ a < |b| ∈ [E/2, E], with no coprimality.

**e₀.** LV §6 bounds the kernel of Ψ by (4/3) log log n. So dim(K/K₀) ≤ log₂(1 + (4/3) log log N). That is tighter than e₀ = O(log log N), and Lemma SEL then needs only O(log log N) square tests.

**Prop SR(a) has precedent** in LV §7 and BLP §9.
- LV §7: "since χ_n(−1) = 1, the definition of χ_P is not dependent on which square roots are taken". They add that if Conj 7.1 fails, NFS "can never find a non-trivial congruence of squares".
- BLP §9: the Hensel root is unique up to sign.
- Parts (d) and the barrier framing are this note's own.

**Attributions.**
- *Rabin.* That two square roots not related by ± split N is elementary (Fermat/Legendre/Kraitchik; BLP §2). Rabin's MIT/LCS/TR-212 (1979) is the reduction form.
- *Square root.*
  - Couveignes, LNM 1554, pp. 95–102, is confirmed (odd degree; the sign is fixed by the norm).
  - Add Montgomery, Proc. Sympos. Appl. Math. 48 (1994) 567–571, and Nguyen, ANTS-III, LNCS 1423 (1998) 151–168, for the lattice method.
  - f′(α)O_K ⊆ Z[α] is correct; LV cite Lang, Ch. III, Prop. 2.
- *Dixon.* Read the fruitfulness bound as "Dixon (≥ 1/2); the exact 1 − 2^{1−ω} follows from the same count over μ₂ ≅ {±1}^ω", conditional on a uniform root.
- *BLP.*
  - Lemma 8.2 is confirmed: k + r uniform vectors span a k-dimensional space with probability ≥ 1 − 2^{−r}. The obstruction analysis is §6 and Thm 6.7.
  - BLP propose the random-character version that LV made rigorous.
- *Descent ladders.*
  - Commeine–Semaev's L_p[1/3, 3^{1/3}] is confirmed.
  - Barbulescu's thesis gives booting at L_p[1/3, 1.23] and descent at L_p[1/3, 1.21] (via Guillevic, Asiacrypt 2015 §1); applying these to factoring is a transfer we have not checked.
  - Gordon 1993's individual-log cost is L_p[1/3, 9^{1/3}]. No source was found for a "Gordon 2/5"; it stays UNVERIFIED.
  - The halving recurrence β_{i+1} − 1/3 = (β_i − 1/3)/2 is the special-Q descent of EGT (and Coppersmith 1984).
  - The fixed-depth closed forms (k+2)/(3k+3), (k+1)/(3k+2), 2^{k+1}/(3(2^{k+1}−1)) and 2^k/(3·2^k−1) were found in no source.
  - **Scope:** the ladders assume this note's no-ECM recognition model. The heuristic literature boots in L[1/3] using ECM, and rigorous HECM costs L_p[2/3] per prime, which is where the L[4/9] objection comes from.
- *Novelty.*
  - No prior work appends hidden-root (Dixon-descended) rows to an NFS matrix to get provable non-triviality.
  - Nearest precedents: Dixon/Rabin fibre symmetry, BLP §9, LV §7, and Stange (arXiv 2211.06821, order-finding index calculus with an NFS variant, Hyp. 3.1).
- *Fixed degree.* Suggested, still unchecked: Bhargava–Varma, Duke 164 (2015), for cubic fields; Ho–Shankar–Varma, Duke 167 (2018), for odd degree.


## 8. Candidate cards

```text
ID: FR-1 hybrid partial descent (LV bulk + h appended hidden-root rows; direct-append architecture).
Proved: Lemma HY (rank r of the projection), Lemma SEL (public selection of a true descended dependency when r >= e0+1),
  Theorem HYB (that dependency is fruitful w.p. exactly 1 - 2^{1-omega}, every odd N with omega >= 2, given transcript invariance).
Rows needed: 1 if absorbed and r0 >= 1 (<= 2 at fixed degree); in general enough for r >= e0 + 1 (<= C - rank A + e0 + 1).
Absorption: Lemma ABS-dom (domination => Pr[d notin V] <= kKC/M); ABS-red (multi-row: e'+1 support-compatible rows).
  E1 (QS): 6-9% single-row absorption, failures = support mismatch; E2 (cubic NFS): 80%.
Cost: one descent into LV's field: model L[(k+2)/(3k+3)] (linear) / L[2^{k+1}/(3(2^{k+1}-1))] (degree-t) at depth k,
  above pure descent at equal depth; > 1/3 at every finite depth.
Barrier: rigorous single-target L[1/3] descent (SD_1) is open; SD_1 (+ obligations) alone gives L[1/3] by pure descent.
Status: proved reduction; refuted for this model (direct append) as an independent route; constant-level device (max(1.923, c)).

ID: FR-2 square-root randomisation (embedding, unit twist, Couveignes/Montgomery signs, inert primes).
Status: refuted (proved, Prop SR): roots are +-v; the class of X/Y is a function of (N, f, m, S); any other valid Y factors N by itself (Rabin).

ID: FR-3 exactly x-preserving, beta-flipping involutions on LV's smooth set (single relations).
Status: refuted (proved, Prop IV-LV): for all but L[2/3, -1/(2 delta)+o(1)] of LV's m, equal exact square classes force
  proportional pairs, which are trivial; beta is constant on x-fibres; criterion (I) vacuous; eta* = 0.
  Batch-/field-level symmetries not excluded.

ID: FR-4 product-law equidistribution of (x, beta).
Status: refuted in TV when beta is non-degenerate (proved, Cor IV-TV); averaged product laws cannot imply fruitfulness
  (countermodel); the Fourier version is Lemma M (open); DC / 8.5H / DW control first-order supply only.
```
