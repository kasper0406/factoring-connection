# below16: can the N^{1/6} floor be escaped?

Agent `below16`, 2026-10-18. Pre-registration: `notes/below16_prereg.md` (commit `7592b7c`,
before any pre-registered run). Code: `src/below16/`. Tests: `tests/test_below16.py`.
Raw data: `experiments/raw/below16_e{1,2,4}.jsonl`. Summary:
`experiments/analysis/below16_summary.json`. GPT-6.1-sol: blind round `codex/below16_1.out`,
adversarial review `codex/below16_2.out`. Independent Opus critic and codex review were both
applied (§11).

## 0. Verdict (plain language)

**No escape below N^{1/6} was found. The floor survives the move to the true-factor
(average-case) task, which is the task Theorem D′ actually solves, for every *local* degree-1
design.**
- Escapes 1, 2 and 5 fail inside the local models T1–T3 below. For escape 2 at degree ≥ 2 the
  proof is conditional.
- Escape 3 fails for a separate, partly heuristic reason: there is no true-factor identity in a
  class group.
- Escape 4 (relation methods) is a real escape. It does not help at the 1/6 level, because it
  jumps straight to subexponential time, and it depends on a derandomization theorem nobody has.

**Main theorem (T1; proved; reviewed by an independent Opus critic and GPT-6.1-sol, fixes
applied).**
- Take a *residue-local* degree-1 design in the 𝔐⁺ cost model, with labels charged:
  - each generator's weights (A, B), class modulus M, key modulus K and unit D are fixed;
  - its class set and window may depend on N only through N mod K;
  - the window may be an *arbitrary finite set*, not only an interval;
  - the hypotheses are D/|A| ≤ C₁M², K ≤ C₂M², and M, K ≤ P/2.
- Then the fraction of balanced semiprimes whose true factor is covered is at most
  C·n²·X^{3/2}·N^{−1/4}. Success on a constant fraction needs X ≥ c·N^{1/6}/(log N)^{4/3}.
  - X is the per-input budget, or the budget averaged over uniform n-bit N. Budgets averaged
    only over S_n are not covered.
  - The proof needs no flatness, no height and no equidistribution of primes.

**Further results.**
- **T2 (proved, reviewed)** adds size information: N to precision ≥ P·K on any partition into
  intervals of length ≥ PK, with interval windows.
- **T3 (GPT-6.1-sol; checked here and by the critic)** allows interval-type windows given by
  algebraic formulas in N (for example "near the stationary point"). Its key is N mod M_τ and
  its width caps are fixed.
- **D′.** Theorem D′'s algorithm is dominated, up to a constant factor in cost, by a T1 design
  (§2.4). So among T1/T2/T3 designs the true-factor exponent is exactly 1/6.
  - D′ is within N^{5ε/6}·polylog of the floor, and Corollary D″ is within polylog.
  - This replaces the overclaim "D′ is optimal in its model", which Thm A alone (full-candidate
    covering) does not support.
- **Not covered by T1–T3:**
  - keys finer than C₂M² (e.g. choosing modulus-ℓ₁ targets from N mod ℓ₂);
  - sparse windows combined with size information;
  - keys finer than N mod M in T3;
  - several labels sharing one giant (MB-b).

| escape | broken assumption of Thm A | verdict | label |
|---|---|---|---|
| 1. bases with known relations (2,3,6; β, β^k; succinct powers) | one base; interval windows | reduces to 1-D targets with arbitrary finite windows (Prop. MB). T1 counts window *cardinality*, so sparse sum-of-AP windows do not help | **refuted for this model**: true-factor task, rational bases, residue-local boxes. Open: size-dependent boxes, worst-case E(d) beyond short bases, algebraic bases, giant multiplicity (§3.3) |
| 2. average higher-order flatness | full-candidate covering | average ≤ sup. Degree 1: no flatness effect on average (T1). Degree ≥ 2: N^{1/5} on average (Cor. F) | **refuted for this model**; the degree ≥ 2 part is conditional on the real-fibre hypothesis (*)_e |
| 3. class-group collisions | the group (Z/N)^* | no small-exponent identity in the true factor. The order-2 class (p, 0, q) is the factor itself. A filter cannot raise coverage. A hybrid needs a new identity C_N^{h(p,q)} = D_N | filter: **refuted** (proved, trivial). Class-number search: N^{1/5} (heuristic). Hybrid: open |
| 4. relation methods | sumset collision search | a genuine escape: deterministic QS heuristically costs L[1/2]. It needs H1 (smooth-value concentration on S_n) and H2 (a nontrivial square root). For Blum N, H2 ⇔ the Jacobi character (X_S Y_S/N) is nontrivial on the kernel. E4: signs look like fair independent coins at 40–80 bits | **conjectured** (heuristic L[1/2]); exact open sub-question O-4 |
| 5. only the true factor needs covering | full-candidate covering | T1/T2/T3 cover local adaptive selectors within their stated keys. A *global* selector shrinking the span to N^{1/2−η} would give N^{1/6−η/3}; none is known | **refuted for local designs** (as scoped); global selector open (O-5) |

**Novelty:** the average-case/local-design floor (T1) and the Jacobi form of H2: "no prior
formulation found by us".
- Codex r2 points to a close recent comparison: Ward–Maestro, *Floors for deterministic integer
  factoring* (preprint, Sept 2026). It studies covering floors including N^{1/6}, in a
  full-range covering framework. It was not read by us.
- Brown 2010 (*lower bounds for straight-line factoring*) works in a different model. The Jacobi criterion for Blum integers is classical (Rabin/Blum). A
deterministic QS is classical, and so is its heuristic L[1/2] analysis.

## 1. Setting and notation

- n is the bit length, P = 2^{n/2}, and S_n = {(p,q): p < q ≤ 2p primes, 2^{n−1} ≤ pq < 2^n}.
- Then p ∈ [P/2, P), q ∈ (P/√2, √2P), and |S_n| ≥ c₀P²/n² for n ≥ n₀ (dcltheorem §1).
- "Covered" means the algorithm's *certified* collision (the identity of multimod §2.1–2.3)
  fires for the true pair. Accidental order collisions are treated separately (§7.3).

**Cost model (multimod 𝔐⁺, unchanged).** For an input N, the design has a set of *targets*
(labels) τ. Each target has a finite window S_τ of exponent values, and a shared baby list
(any sumset structure) covers them. Then #targets ≤ X and Σ_τ|S_τ| ≤ X². The generic-group
argument justifies this: every useful collision test involves one target handle and tests one
value, so ≥ Σ|S| tests cost ≥ √(Σ|S|) group operations, and every handle costs one operation.

## 2. The true-factor floor for local designs

### 2.1 Generators

A **T1 generator** ι consists of fixed integers A, B (not both 0), M ≥ 1, K with M | K, D ≥ 1, a
side ∈ {p, q}, and arbitrary maps κ ↦ R(κ) ⊆ Z/M and κ ↦ S(κ) ⊂ Z (finite) on keys
κ ∈ Z/K. Write x for the class variable (x = p on side p, x = q on side q) and y for the
partner. On input N = pq the generator issues |R(κ)| targets, κ = N mod K, and these cover
(p, q) iff

    x mod M ∈ R(κ)   and   ⌊(Bx + Ay)/D⌋ ∈ S(κ).                                  (2.1)

**Why this form is general for degree-1 targets.**
- By multimod Lemma K′ (§3.1(i)), every identity-valid degree-1 two-term target is
  g = (Bx + Ay − c)/D, with D = C_d M², c ∈ [0, D) fixed by the class and N mod D.
- So for true pairs g = ⌊(Bx + Ay)/D⌋, and any window placement is some finite set S.
- One-term targets (A = 0 or B = 0) are included.
- HH-type designs whose weights depend on the class and on N mod M decompose into one
  generator per weight vector.
- Validity of the lock is never used in the count. The bound holds for every rule of shape
  (2.1), valid or not.

**Locality hypothesis (L1).** D/|A| ≤ C₁M² for A ≠ 0, and D/|B| ≤ M for A = 0 (one-term
validity gives D | BM). These are scale-invariant and are implied by Lemma K′ for reduced
representations. Also K ≤ C₂M²: the design looks at N only modulo K. The K ≤ C₂M² hypothesis is needed: GPT-6.1-sol's counterexample (§7.1) shows
that an unrestricted key can encode N itself. Size conditions: M ≤ P/2, K ≤ P/2, K ≤ 2^{n−2}.

**Lemma 1 (fixed-key partner count; proved).** Fix x and a unit key κ. The partners y ≡ κx^{−1}
(mod K) of x that satisfy (2.1) number at most

    φ(κ) := min( |S(κ)|·(⌊D/(|A|K)⌋ + 1),  ⌊ℓ_y/K⌋ + 1 )          (A ≠ 0),

where ℓ_y ≤ x is the length of the partner range.

*Proof.* The y's form an arithmetic progression of step K. For each s ∈ S, ⌊(Bx + Ay)/D⌋ = s
confines y to a half-open interval of length D/|A|, which contains at most ⌊D/(|A|K)⌋ + 1 terms.
The progression has at most ⌊ℓ_y/K⌋ + 1 terms in the partner range. ∎

**Theorem T1 (residue-local designs; proved).** Under (L1), for every T1 generator,

    Count(ι) := #{(p,q) ∈ S_n covered by ι} ≤ 2c_s (P/M)·√(2P(D/|A| + K))·E_N[ |R(N mod K)|·|S(N mod K)|^{1/2} ]
                                           ≤ c₁ P^{3/2} E_N[ |R|·|S|^{1/2} ].

- Here c_s = 1 (side p) or 5/4 (side q), and c₁ = max((5/2)√(2(C₁+C₂)), 5√(C₁+1)). The second
  entry covers A = 0 (critic fix 1).
- E_N is the uniform average over **all** n-bit integers N, not over S_n.
- (L1) holds *after canonical normalisation* of each target (Lemma K′'s representation over
  C_d M²). Unreduced representations can be renormalised first. The budgets are the two
  *moment* budgets E#τ ≤ X and EΣ|S| ≤ X². An expected running time E C_N ≤ X does not imply
  the second.
- Not covered: adaptation driven by many small residue queries, whose CRT combination can
  determine N, and globally chosen weights whose locks happen to be local.
- For A = 0 the same holds with 4c_s(P/M)√(P(D/|B| + M)) in front.

Consequently, for any design (any set of generators) obeying the cost model per input or on
average over uniform n-bit N (E_N[#targets] ≤ X, E_N[Σ|S_τ|] ≤ X²),

    Pr_{(p,q) ∈ S_n}[true pair covered] ≤ (c₁/c₀)·n²·P^{−1/2}·E_N[Σ_τ |S_τ|^{1/2}] ≤ (c₁/c₀)·n²·X^{3/2}·N^{−1/4}.      (2.2)

So success probability σ forces X ≥ (σc₀/c₁)^{2/3}·n^{−4/3}·N^{1/6}.

*Proof.*

1. **Sum by keys.** Both primes exceed K, so they are units mod K. For each x and each unit
   κ there is exactly one partner class a = κx^{−1}. Lemma 1 gives

       Count(ι) ≤ Σ_x Σ_{κ: x mod M ∈ R(κ)} φ(κ) = Σ_κ φ(κ)·#{x : x mod M ∈ R(κ)},

   with x ranging over the integers of the class-variable range. That range has length
   ≤ P/2 (side p) or ≤ 0.71P (side q), so #{x : x mod M ∈ R(κ)} ≤ c_s|R(κ)|·P/M, using M ≤ P/2.
2. **Bound φ.** Use min(u, v) ≤ √(uv), with u ≤ |S|(D/|A| + K)/K and v ≤ ℓ_y/K + 1 ≤ 2P/K
   (since K ≤ P/2). This gives φ(κ) ≤ √(2P|S|(D/|A| + K))/K.
3. **Pass to E_N.** Every residue class mod K contains ≥ ⌊2^{n−1}/K⌋ ≥ 2^{n−2}/K integers
   N ∈ [2^{n−1}, 2^n). Hence Σ_κ f(κ) ≤ 2K·E_N f(N mod K) for f ≥ 0. Combining with steps 1–2
   gives the first inequality.
4. **Closed form.** D/|A| ≤ C₁M² and K ≤ C₂M² give √(D/|A| + K) ≤ M√(C₁+C₂). This is the second
   inequality.
5. **A = 0.** *Correction (codex r2).* The earlier direct bound used P/M without c_s and fails
   on side q: n = 20, M = K = 348, R(65) = {47}, S(65) = {0} give 3 admissible q against
   P/M ≈ 2.94.
   - *Repair: swap sides.* Define R′(κ) = {κr^{−1} mod M : r ∈ R(κ)}. Put the generator on the
     opposite side with (A′, B′) = (B, 0) and the same D, K, S.
   - Coverage is identical, |R′| ≤ |R|, and A′ ≠ 0, so the A ≠ 0 case applies with the same c₁.
   - (L1) for this case reads D/|B| ≤ C₁M².
6. **Sum over the design.** The covered pairs are contained in the union over generators.
   Class-dependent windows are handled by splitting R into one generator per class index.
   Summing step 4 over generators gives Σ_ι E_N[|R_ι||S_ι|^{1/2}] = E_N[Σ_τ |S_τ|^{1/2}],
   since each class in R is one target.
7. **Cost budget.** Cauchy–Schwarz per input, then in expectation:

       E[Σ_τ √|S_τ|] ≤ E[√(#τ·Σ|S_τ|)] ≤ √(E#τ·EΣ|S|) ≤ X^{3/2}.

   Divide by |S_n| ≥ c₀P²/n². ∎

**Remarks.**
- The proof counts *integer* pairs, so primality costs only the n² normalisation. No prime
  equidistribution is used; prime-density thinning saves only logarithms (r16).
- There is no H and no flatness. On average, the class lock's "whole class" bound P/M and the
  per-value bound |S|·D/(AK) combine by a geometric mean. (Windows placed using √N, e.g. at
  the stationary point, are outside T1 and are handled by T2/T3. They also give only √|S|.)
- Arbitrary finite windows are allowed. A sparse window of cardinality s counts as s, never as
  its hull. This is what closes escape 1 (§3).

### 2.2 Coarse size information (T2)

A **T2 generator** has side p and A, B ≥ 1. Its key is (κ, b) = (N mod K, ⌊(N − 2^{n−1})/Δ⌋),
with Δ = 2^j ≥ PK, and its window is an interval [U, U + W) depending on the key. More
generally, the design may use any partition of [2^{n−1}, 2^n) into intervals of length between
PK and 2PK; Lemma 2's proof is unchanged (critic). Require K ≤ P/2. Side q is analogous. That is about n/2 bits of N in total: log K bits of residue and n/2 − log K bits of size.

**Lemma 2 (transfer; proved, any target type).** For a key-box let C⁺ be the set of integers x
in the class set that have *some* partner in the box satisfying the window. Then

    Count(ι) ≤ (3·2^n/P)·E_N[ |C⁺(key(N))| ],   i.e.   Pr ≤ (3/c₀) n² E_N[Σ_τ |C⁺_τ|]/P.

*Proof.*
- Per (x, key), the partners lie in an interval of length ≤ Δ/x ≤ 2Δ/P in one class mod K.
  There are at most 2Δ/(PK) + 1 ≤ 3Δ/(PK) of them.
- Each key-box contains ≥ Δ/(2K) integers N.
- So Σ_{keys} |C⁺| ≤ (2K/Δ)·2^{n−1}E_N|C⁺|. ∎

**Theorem T2 (proved).** For T2 generators,

    |C⁺| ≤ |R|·( 2√(2 I P/A)/M + 2 ),   I = DW + 2AΔ/P.

Hence

    Pr[covered] ≤ (3/c₀) n² ( 4√C₁·X^{3/2}P^{−1/2} + 4√(2C₂)·X P^{−1/2} + 2X/P )

(using D ≤ C₁M² and Δ ≤ 2PK ≤ 2C₂PM²). So success again needs X ≥ c·N^{1/6}/n^{4/3}.

*Proof.*
- If x ∈ C⁺, some N' in the box has F_{N'}(x) = Bx + AN'/x ∈ D[U, U+W).
- For N₀ the box's left end, F_{N₀}(x) lies in an interval I′ of length ≤ DW + 2AΔ/P,
  because x ≥ P/2.
- F_{N₀} is convex with F″ = 2AN₀/x³ ≥ A/P on [P/2, P].
- A convex function with F″ ≥ κ maps a set of measure ≤ 2√(2|I′|/κ) into I′: two intervals,
  and √a − √b ≤ √(a − b).
- Each class mod M meets two intervals of total length L in ≤ L/M + 2 integers.
- Sum with Cauchy–Schwarz as in T1. ∎

### 2.3 Algebraic windows (T3, GPT-6.1-sol) and degree ≥ 2 (Cor. F)

GPT-6.1-sol (below16_1, Theorem 3) proves a third variant. I checked the proof (sketch below).

**Model.** These are the full hypotheses of below16_1 Theorem 3, restored at codex r2's request.
- There is a fixed schedule of ≤ X slots.
- Each slot τ has a fixed modulus M_τ ≤ P, one unit class r_τ(s) per residue s = N mod M_τ,
  and a region R_{τ,s} in the (p, q) box. The region is semialgebraic of bounded complexity
  and may be given by formulas in N, for example "|Bp + Aq − 2√(ABpq)| < DW".
- There is a residue-independent width cap w_τ ≥ 1 with Σ w_τ ≤ X². This is stronger than the
  active-target budgets of T1.
- The semialgebraic complexity is bounded uniformly in n.
- The *uniform real-fibre estimate* holds: length{x : (x, u/x) ∈ R_{τ,s}} ≤ C M_τ(√(P w_τ) + w_τ)
  for every real u ≍ P².

**Conclusion.** The covered fraction is ≤ C n²(X^{3/2}P^{−1/2} + X²/P).

*Proof.* Coarea in (u = pq, x = p) gives area(R_s) ≤ C M(P^{3/2}√w + Pw). A lattice translate
of MZ² meets a bounded-complexity region in ≤ area/M² + C(P/M + 1) points. Sum over the
≤ M residues.

For degree-1 targets the fibre estimate *is* the real form of multimod (1), which is proved.
T3 complements T1 and T2:
- T1: arbitrary sets, residue keys only.
- T2: intervals, n/2 bits of size.
- T3: intervals given by algebraic formulas in N, with fixed width caps.

**Corollary F (average-case degree ≥ 2; proved under the uniform real-fibre hypothesis
(*)_e).** Assume
Write k for the support size, to avoid a clash with the key modulus K. Assume one of:
- (*)_e: fibre length ≤ C M P^{1−e/k} w^{1/k} for degree-e, support-k targets, for every real
  u ≍ P². This is the real-variable form of multimod (2).
- In the T2 model: every window is at least as wide as its *slack*, defined as
  sup_{N in the key-box} sup_x |g_N(x) − g_{N₀}(x)|, where N₀ ≡ κ (mod K) is a point of the box at
  which Lemma K′ applies.

Then targets of degree ≥ 2 alone cover a fraction ≤ C n²(Σ_{e,k} P^{−e/k} X^{1+1/k} + X/P), so
positive success needs X ≥ N^{1/5}/(log N)^{O(1)}.

*Proof.* T3's area-plus-lattice argument (which contributes the X/P term), or Lemma 2 with
multimod (2) applied to C⁺ on a window of width ≤ 3w. Then Hölder:
Σ w^{1/k} ≤ X^{1+1/k}, and e/(k+1) ≥ 2/5 for k ≤ 2e. ∎

**So escape 2 fails, conditionally at degree ≥ 2.** It is the *local counting hypotheses*
that prevent a selector from always choosing the right small set, not "average ≤ sup" by
itself (codex r2). Small per-input coverage is compatible with perfect true-factor coverage if
the selector knows p. Degree 1 stays at 1/6 (proved). Degree ≥ 2 stays at 1/5 under (*)_e or the T2
slack condition.

### 2.4 Theorem D′ lies in the model

DCL-det's targets are (A, B, ℓ, root) with D = ℓ² and class set = the roots of Br² ≡ AN
(mod ℓ) (key N mod ℓ). Its window depends on N mod ℓ² and on √N.
- Replace each window by the n-only superset [⌊F_min/ℓ²⌋ − 1, ⌊F_max/ℓ²⌋ + 1] over all n-bit N.
  This gives a T1 design with K = ℓ.
- *Correction (Opus critic).* An earlier draft said the superset "at most doubles Σ|S|". That
  is false. At n = 80 the ratio is 2.3–3.7 in total and up to ≈ 29 for single weights with
  A/B ≈ 1/√2, whose per-N span is only ≈ 0.025·B√N.
- *Repair.* The superset span is ≤ ≈ 2.12·HP < span̄. So each superset width is ≤ m̄_k + 1, and
  Σ|S| ≤ Σ_k T_k(m̄_k + 1) ≤ (Σ_k(m̄_k + T_k))². That is the square of the list cost D′ already
  pays, so the T1 budget is X = O(cost of D′).
- Hence D′ (cost N^{1/6+5ε/6}·polylog) is within N^{5ε/6}·polylog of the T1 floor
  c N^{1/6}/n^{4/3}, and Corollary D″ (N^{1/6}·polylog) is within polylog.

**Corollary (proved for degree-1 designs).** Among degree-1 residue-local (T1), size-local (T2)
or algebraic-window (T3) designs in the 𝔐⁺ cost model with labels charged, the true-factor
exponent is exactly 1/6, attained by D″. Designs that mix in degree ≥ 2 targets are covered
conditionally on Cor. F's hypothesis.

## 3. Escape 1: several bases with known relations

### 3.1 Reduction (Prop. MB; proved)

**Setup.**
- Let the bases α_1..α_k be positive rationals (any number, any known relations).
- The reduction treats the *exact rational multiplicative relations among the fixed bases*:
  the kernel of the valuation matrix. A false fixed rational equality survives modulo only
  finitely many primes.
- Input-dependent modular relations are excluded from this reduction. Examples: bases 2 and
  N + 2 agree mod N, and 2^{p−1} ≡ 1 has a prime-dependent exponent.
- Every useful target needs a nonconstant valuation projection. Bases (2, 4) with forms
  (2p, −p) give a constant product with all projections constant (codex r2).
- In prime-valuation coordinates (one per prime dividing some α_i), an identity-provable
  collision of ∏α_i^{u_i} with the baby/giant product means E^T u = 0, where E is the k × r
  valuation matrix.

**Statement.** A target with integer exponent forms F_i(p, q), offsets c_i and a box (or any
set) Q ⊂ Z^k covers (p, q) iff E^T(F(p,q) − c) ∈ E^T(Q). For each coordinate j with nonconstant
L_j(F) = Σ_i E_{ij}F_i, the coordinate target covers a superset, with window
S_j = π_j(E^T Q) and |S_j| ≤ |E^T Q| ≤ |B|·|G_τ|.

**Consequence.**
- The multi-base design maps to a 1-D design with the same number of targets and
  Σ|S_j| ≤ X², but with *arbitrary finite* windows (sums of APs, e.g. {u + Kv} for bases
  β, β^K). The offsets go into the window: S_j = π_j(E^T(Q + c)).
- *Missing step (Opus critic; now supplied).* The individual F_i need not be valid on their
  own: multimod E(d)'s worry was that they need only be valid modulo Λ_known. What T1 needs is
  validity of the coordinate form w_j := (E^T F)_j.
  - The identity says ∏α_i^{e_i(N)} ≡ ∏α_i^{F_i(p,q)} (mod p), where e(N) are the handle
    exponents.
  - In valuation coordinates this is (E^T e(N))_j − w_j(p, q) ∈ (p − 1)·Z for every j with the
    virtual base π_j.
  - So w_j is itself a one-base, degree-1, integer-valued target on the lock class, and Lemma K′
    applies to it.
- *Hypotheses this imposes.* The denominator of w_j is ≤ C_d·M′², where M′ is the lcm of the
  coordinate lock moduli, and the key is ≤ C·M′² (lcm(K_i) ≤ C₂^k·M′²). Under these, T1 applies
  to (w_j, S_j). Succinct powers are fine: a huge k only enlarges |A′|, and T1's bound
  decreases in |A|.

**Verdict.** In the true-factor model, multi-base searches with rational bases and arbitrary
known relations, whose boxes are **residue-local** (they depend on N only through N mod K), do
not beat N^{1/6}/polylog. This includes succinct powers β^k with huge k, since T1 uses |S|, not
hull(S).
- This closes multimod Thm E(d) *for residue-local boxes in the true-factor task only*.
- Natural multi-base boxes depend on √N, and they project to sparse windows with size
  information. That combination is the §8 model gap: T2/T3 need intervals.
- E(d) itself (worst case, full candidates) stays open beyond GPT-6.1-sol's Theorem 1 for short
  bases.

- E2 checks the reduction exactly (§9.2).
- GPT-6.1-sol's Theorem 1 (below16_1) is the *worst-case* full-candidate analogue. It holds for
  short bases (b-bit numerators and denominators) and loses √(1 + kb). For huge powers its
  hull argument fails. T1 is the average-case statement that does not.

### 3.2 Why not Smith normal form

The orchestrator's hint suggests Smith normal form (SNF).
- SNF gives an abstract free basis, but intervals become sums of APs, so the interval
  geometry is lost.
- Prime coordinates give the same reduction without a basis change.
- The important point is that T1 does not need intervals.

### 3.3 What remains open (precise)

- **(MB-a) Algebraic bases.** Units or elements of a number field K. For split p the
  reduction is identical. For inert p the group has order p^f − 1, and handles involve an
  unknown Frobenius. Cyclotomic handles are covered by cycloharvey D1; general non-Galois K is
  untested.
- **(MB-b) Giant multiplicity (affects one-base designs too).** The cost model charges labels,
  not stored group elements, and #labels ≤ X is a *model hypothesis* of T1–T3. The
  generic-group remark in §1 does not prove it: an operation creates one element but may carry
  several certified interpretations. One
  giant G could carry R labels: certified identities G ≡ γ^{g_k(p)} on different classes.
  - Then Count grows by R at no extra cost.
  - In the HH-lock family, two labels share a giant only if their integer handle exponents
    coincide *on that input* (heuristic, not a size-separation proof). For distinct band primes this forces ℓ'²A = ℓ²A′, which is impossible. For equal
    ℓ it forces equal A and a congruence on c that holds only for a ~1/ℓ² fraction of N
    (heuristic).
  - So multiplicity is negligible there. A general bound on labels per giant is
    **not proved**.

## 4. Escape 2: average-case flatness

Covered by §2.3 (Cor. F, T1 remark). In one line: average coverage ≤ sup coverage. Locality
gives the uniformity step, and the degree-1 average bound does not even see flatness.
**Refuted for this model.** It is conditional on the real-fibre form (*)_e for degree ≥ 2, or on
windows at least as wide as the slack in T2.

What is left is a selector that places flat neighbourhoods *around the true factor*. That is
escape 5's global selector (O-5).

## 5. Escape 3: class groups

1. **No small true-factor identity in the constructions examined.** This is an observation,
   not a theorem: the earlier heading "proved" overstated it.
   - The unit-group identity works because N ≡ q (mod p − 1) is public arithmetic in the
     exponent.
   - Cl(−4N) (or Cl(−4kN)) has order h(−4kN) ≈ √N·L(1, χ) with no expression in p, q. There is
     no CRT projection of Cl(−4N) to a factor-dependent group with a known exponent.
   - The factor does give an order-2 class [(p, 0, q)] = [I_p], I_p = (p, √−N) (GPT-6.1-sol
     proof: I_p² = (p), and I_p is not principal since a² + Nb² = p is impossible). But writing
     it down *is* knowing p.
2. **Filtering cannot help (proved, trivial).** A second group used to test unit-group targets
   accepts a subset of their covered set. Coverage, and hence T1, is unchanged.
3. **Class-number search (heuristic, for this one search architecture; not a lower bound).**
   - The generic route computes h (or the order of a form) to find ambiguous forms. Shanks's
     BSGS on h ∈ [ĥ − E, ĥ + E] costs x + √E, where ĥ is the Euler product to x.
   - On average over discriminants, the tail Σ_{x<ℓ≤2x} χ(ℓ)/ℓ has variance ≍ 1/(x log x) and is
     asymptotically independent of χ(ℓ) for ℓ ≤ x (orthogonality of quadratic characters).
     So E ≫ √N/√(x log x) on a positive fraction of inputs.
   - So min_x (x + N^{1/4}x^{−1/4}) = N^{1/5}, the GRH exponent of Schoof (1982).
   - A real lower bound would need conditional anti-concentration of the estimation error,
     given everything the algorithm computed, on discriminants −4pq from S_n. Variance and
     pairwise orthogonality do not supply this, and other estimators and non-generic algorithms
     are not addressed (codex r2).
   - Generic collisions in Cl cost h^{1/2} ≈ N^{1/4} (r16).
4. **The exact open hybrid (O-3).** A public relation C_N^{h(p,q)} = D_N in some group, with
   C_N, D_N computable from N, h(p, q) in a searchable set cheaper than N^{1/6}, and one-sided
   extraction. None is known.

## 6. Escape 4: relation methods (the real escape)

### 6.1 Conditional theorem (heuristic in its hypotheses)

The deterministic single-polynomial QS of `src/below16/qs_pilot.py` is uniform, deterministic
and has a fixed stopping rule. It uses Q(x) = (⌈√N⌉ + x)² − N, a factor base of −1, 2 and the
primes ℓ ≤ y with (N/ℓ) = 1, exact trial division, and exact F₂ elimination. With
y = L[1/2, 1/(2√2)] and the interval L[1/2, 3/(2√2) + ε] (GPT-6.1-sol's dense-elimination
schedule):

**Proposition 6.1 (proved conditional on H1, H2; the costs are proved).** Suppose that on a
1 − o(1) fraction of S_n:
- (H1) the interval contains ≥ π(y) + 2 smooth values;
- (H2) some dependency gives a nontrivial square root of 1.

Then the algorithm factors those inputs in L_N[1/2, √2 + ε + o(1)] = N^{o(1)}.

**Exact form of H2 for Blum N (proved).**
- Let p ≡ q ≡ 3 (mod 4). For a kernel vector S, put X_S = ∏(s + x_i) and
  Y_S = (−1)^{e_{−1}(S)/2}∏_{ℓ≥2} ℓ^{e_ℓ(S)/2} (mod N).
- Assume gcd(Y_S·∏(s + x_i), N) = 1; otherwise a factor is already exposed.
- Then gcd(X_S − Y_S, N) is proper ⇔ the Jacobi symbol (X_S Y_S / N) = −1.
  - X_S/Y_S is a square root of 1. The nontrivial roots of 1 have Jacobi symbol −1, the
    trivial ones +1 (because (−1/N) = 1).
- S ↦ (X_S Y_S/N) is a character of the kernel under symmetric difference, since the overlap
  correction is (s + x_i)²/Q(x_i) ≡ 1 (mod N).
- Scope: the Jacobi detection needs Blum N, which are ≈ 1/4 of S_n. For other N, H2 must be
  phrased through the sign pair itself.
- **So H2 ⇔ the computable character ψ_N(S) = ∏_{i∈S}((s+x_i)/N) · ∏_ℓ (ℓ/N)^{e_ℓ(S)/2} is
  nontrivial on ker.**
- H2 is thus a statement about Jacobi symbols of explicit integers, and heuristically fails
  with probability 2^{−dim ker}.

### 6.2 What a proof needs (the exact open sub-question O-4)

1. **H1 on S_n.** For almost all (p, q) ∈ S_n, at least π(y) + 2 of the Q(x) are y-smooth.
   - First moment: smooth numbers in the intervals [2|x|√N, 2|x|√N + O(√N)] over the
     archimedean distribution of N. This is plausible via Hildebrand-type short-interval
     results, but they need y ≥ exp((log N)^{5/6+ε}), larger than L[1/2]. For a mere
     sub-1/6 exponent, polynomial-size y (u fixed) suffices: §6.3.
   - Second moment: joint smoothness of Q(x), Q(x′) as N varies. This is a binary additive
     problem (semiprime + smooth = square) and is open.
2. **H2.** Non-triviality of ψ_N on ker for almost all N. Because ker depends on N, this needs
   independence between the Jacobi symbols ((s + x)/N) and the smoothness pattern. Standard
   tools (Heath-Brown's quadratic large sieve) handle *fixed* integers, not N-dependent kernels.
   **Open.**
3. **Non-uniform shortcut (known; it does not help).** By averaging, a randomized Las Vegas
   algorithm (Lenstra–Pomerance, L[1/2, 1 + o(1)]) yields, for each n, *some* fixed tape that
   works on all but δ of S_n. By Adleman's argument some tape even works on every n-bit N.
   That is a non-uniform deterministic L[1/2] algorithm. The project's question is about
   *uniform* algorithms, and no efficient way to find the tape is known.
4. **Trap (GPT-6.1-sol).** Fixed samples x_j = j < √N give x_j² mod N = j², so every
   dependency is trivial. Deterministic schedules must create real wrap-around.

### 6.3 A polynomial-level version (relevant to "below 1/6", heuristic)

Below 1/6 does not need subexponential time. Take y = N^{c}.
- Q(x) for |x| ≤ N^{c+o(1)} has size N^{1/2+c+o(1)}. It is y-smooth with probability ≈ ρ(u),
  u = (1/2 + c)/c = 1/(2c) + 1, which is a *constant*.
- So π(y) + 2 relations come from an interval of length N^{c+o(1)}.
- Deterministic dense elimination over F₂ by fast matrix multiplication costs N^{ωc+o(1)},
  ω < 2.372 (LSP/LUP decomposition, Ibarra–Moran–Hui 1982, is deterministic). The sparse
  methods, Wiedemann and block Lanczos, are randomized.
- Smoothness detection must be exact (sieve by all prime powers ≤ max|Q|), at cost
  N^{c+o(1)}. Trial division of every value would cost N^{2c}.
- Total cost N^{max(1, ω)c + o(1)}, which is < N^{1/6} for c < 1/(6ω) ≈ 0.070 (u ≈ 8.1).
  Plain cubic elimination needs c < 1/18. Sparse/structured elimination at N^{2c} has no proved
  deterministic bound, so c < 1/12 would be heuristic even in its cost. An earlier draft's
  "c < 1/12, N^{3c}" was wrong (critic, codex).
- A proof at this level needs only H1 at **fixed u** (a constant smooth density for
  N^{1/2+c}-sized values of one explicit polynomial, on average over N ∈ S_n) plus H2. That is a
  cleaner target than L[1/2], but it is still **open**: the second moment is the same binary
  problem.

### 6.4 Experiment E4 (pre-registered)

See §9.3. The sign statistics behave like fair independent coins.

## 7. Escape 5: only the true factor needs covering

1. **Locality is necessary (GPT-6.1-sol; proved).** With an unrestricted key modulus
   (M = 2^{n+1} > N), "depends only on N mod M" means "depends on N". An uncharged selector
   could factor N and place a single width-1 target (A, B, r) = (p, q, p). Hence (L1) and
   K ≤ P/2.
2. **Local adaptivity is useless within T1/T2/T3.** Choosing classes, weights or windows from
   the stated local data does not beat N^{1/6}/polylog in certified coverage. That data is:
   residues mod K ≤ C₂M² with arbitrary windows; or about n/2 bits of size with intervals; or
   algebraic formulas with fixed caps. Combinations of these, and combined residue queries,
   are not covered.
3. **False positives / accidental collisions (explicit heuristic for DCL-type differences).**
   - A non-certified collision mod p needs ord_p(β) | D for a tested difference D. In the
     random-order model the accidental rate is ≤ (#tested pairs)/ord_p ≤ X²·n^{O(1)}/P, which is
     o(1) for X ≤ N^{1/4−δ}. This *assumes* tested differences are uniform modulo ord_p(β). Even
     for random β, Pr[β^d = 1] = gcd(d, p−1)/(p−1), which structured annihilators (p − 1 methods)
     exploit. D′'s data agree: 83 accidental first hits against a null of 164.
   - Two-sided (FULL) collisions mod N give a multiple of ord_N(β), hence a factor by Miller's
     reduction when the 2-adic orders differ. They are rare in D′ (≤ 0.045%), which does not
     bound other designs.
   - Smooth-order annihilation (p ± 1, ECM schedules) is a separate route, outside the
     collision model. A deterministic ECM schedule needs a smooth-order hitting-and-separation
     theorem for almost all p. Open, like O-4.
4. **Global selector (O-5, open).**
   - A cheap public procedure that localises p to an interval of length N^{1/2−η} would give
     N^{1/6−η/3} plus the selector cost (r16; GPT-6.1-sol §5).
   - With η ≥ 1/4 it is already a polynomial-time factorer (Coppersmith).
   - The selector must output one interval (or a controlled number of intervals with
     efficiently supplied endpoints). A "set of small measure" is vacuous, since a finite set
     has measure 0.
   - Any such selector must go beyond T1–T3's local data.
5. **Per-target lattice solving (open; not analysed).** *Correction (codex r2).* The draft
   called this "information-equivalent to known p mod M", which is wrong: a narrow window adds
   archimedean information. For example, with M = 1, a singleton target 2p + q = s is solved
   from 2p² − sp + N = 0.
   - The residue and window restrictions must be analysed jointly, then the selected targets
     costed. M ≈ N^{1/4} is only the standard *sufficient* threshold from a residue alone.
   - No claim either way.

## 8. Precise open sub-questions

- **O-4 (relation route; one concrete conditional escape examined here).** Prove H1 and H2 for one explicit
  deterministic schedule on almost all of S_n. The polynomial-level variant (y = N^{c},
  c < 1/(6ω) with deterministic fast-matrix-multiplication elimination, or c < 1/18 with cubic
  elimination; constant u) suffices for N^{1/6−δ}. The smoothness detection must be exact
  (sieve by all prime powers): the pilot's integer-weight filter can discard smooth prime powers
  (codex r2 example, Q = 2^{30}). So H1 must be stated for the retained relations, or the sieve
  made exact. In H2's Jacobi form: show that
  ψ_N(S) = ∏((s+x_i)/N)∏(ℓ/N)^{e_ℓ(S)/2} is nontrivial on the kernel for almost all Blum N.
- **O-5 (global selector).** Find any public, N^{1/6−δ}-time computable map N ↦ (one interval,
  or N^{o(1)} intervals, of total length N^{1/2−η}) that contains p for most of S_n.
- **MB-a / MB-b.** Algebraic (non-rational) bases; a labels-per-giant bound (also for one base).
- **Smooth-order route.** A deterministic ECM / p ± 1 schedule with a proved hitting theorem.
- **O-3.** A class-group (or other group) identity C_N^{h(p,q)} = D_N with a cheaply searchable
  h.
- **Model gap.** Prove T1 with keys N mod K for K up to P^{1−δ}, instead of K ≤ C₂M², and with
  arbitrary windows plus fine size information simultaneously.

## 9. Experiments (pre-registered, notes/below16_prereg.md)

### 9.1 E1 — exact bound checks

Run: `experiments e1`, 1424 generators over n = 20–30 (exhaustive S_n). The families were
DCL with K = ℓ and K = ℓ²; 900 random T1 generators with ≈ 12k interval, 12k AP, 12k 2-dim GAP
and 12.5k random-set windows; 240 oracle-local; and 32 oracle T2.

- **P1a, P1b and P1c hold: 0 violations** of count ≤ core, core ≤ closed and count ≤ closed.
  No generator failed the size conditions.
- The closed form is loose by about 10³. max count/closed ≤ 5.2·10⁻³ in every family, and it
  falls with n, as the n² prime-density factor predicts.
- The adversarial oracle-local windows (best class and best s values per key) come no closer
  to the bound than DCL does: max count/closed is 3.6·10⁻⁴–5.2·10⁻³.

### 9.2 E2 — multi-base reduction

300 targets over 5 base sets: **0 superset failures** (P2a). |π_j E^T box| ≤ |E^T box| ≤ |box|
always holds (P2b).

### 9.3 E4 — deterministic QS sign statistics

All cells are complete: n = 40–72 with 400 Blum + 200 general instances each, and n = 80 with
200 + 100.
- **P4a:** the Jacobi rule held for 121 921/121 921 basis vectors.
- **P4b:** the per-vector nontrivial rate was 0.481–0.507 in every cell.
- **P4c:** the minimum prefix p-value per cell ranged from 0.004 to 0.38, so every cell met
  the frozen threshold of 0.001. The smallest values came from a mild excess of all-trivial
  prefixes at 40 bits (Blum, j = 3: 70 observed against 50 expected); there was no excess at
  48–80 bits.
- **P4d:** success was 100% in all 12 cells (3 000 instances).
- **P4e:** there were 0 ACCIDENTAL hits.

**Decision (pre-registered rule): sign independence (H2) is experimentally supported at 40–80
bits** for this deterministic schedule. The 40-bit excess is noted as a possible small-size
effect; it is unexplained and was not tested further.

**Caveat (codex r2).**
- The pilot's integer-weight sieve filter can discard smooth prime-power values, so the
  relations it retains are a subset of the smooth values. E4 tests H2 on the *retained*
  relations.
- "2^{−dim ker}" is a random-character heuristic: basis signs are not intrinsically
  independent.
- Blum N are ≈ 1/4 of S_n. The general (non-Blum) cells were checked by gcd only.

## 10. Candidate cards (BRIEF §4)

**Card B16-1: the true-factor floor for local designs (an obstruction).**
- *ID and mechanism.* Any sumset collision search whose targets are chosen locally (from N mod
  K, ≤ n/2 bits of size, or algebraic formulas in N) covers ≤ C n² X^{3/2} N^{−1/4} of S_n.
- *Target / input family.* Balanced semiprimes S_n. This is the task of Theorem D′, the true
  factor counted over pairs.
- *Representation.* (Z/N)^*, one or several rational bases, exponent windows of any shape.
- *Exact operation.* Shared baby list plus giants per target (𝔐⁺ cost model).
- *Lemmas used.* multimod Lemma K′ (denominators), Lemma 1 (fixed-key partner count),
  convexity (T2), coarea plus lattice counting (T3).
- *Status.* T1 and T2 **proved**. T3 proved by GPT-6.1-sol and checked here. Cor. F is proved
  under (*)_e. Experimentally checked: 1424 generators, 0 violations (E1).
- *Cost.* None (it is a lower bound). It is tight up to polylog, since D′ attains it.
- *What would refute it.* A local design beating (2.2). E1's adversarial "oracle-local"
  windows got nowhere near the bound.
- *What escapes it.* Global information about N (O-5), relation methods (O-4), or polynomial
  label multiplicity per giant (MB-b).

**Card B16-2: multi-base searches with known relations (escape 1).**
- *Mechanism.* Bases 2, 3, 6 / β, β^k / γ_ℓ = β^{L/ℓ²} share one baby list in the free group
  they generate.
- *Operation.* Prime-valuation coordinates turn the search into 1-D targets whose windows are
  sums of APs.
- *New lemma.* Prop. MB (**proved**), plus T1's cardinality form.
- *Outcome.* **Refuted for this model**. E2 checked 300 targets with 0 failures.
- *Repair if one wanted to revive it.* Find bases whose relations are known only mod N (not in
  Q^*) or polynomial giant multiplicity. Neither is known.

**Card B16-3: deterministic quadratic sieve on average (escape 4).**
- *Mechanism.* Q(x) = (⌈√N⌉ + x)² − N, a fixed factor base, exact F₂ kernel, and a congruence
  of squares.
- *Target.* Almost all of S_n, deterministic and uniform.
- *Extraction.* gcd(X_S − Y_S, N). For Blum N this is nontrivial ⇔ (X_S Y_S/N) = −1.
- *Cost.* L[1/2, √2 + ε] with dense elimination, heuristic. At the polynomial level,
  y = N^{c} with c < 1/(6ω) ≈ 0.070 gives N^{ωc+o(1)} < N^{1/6}, also heuristic (§6.3).
- *Randomness.* None. The open part is proving success for almost all inputs.
- *Closest prior method.* QS and CFRAC (Pomerance; Morrison–Brillhart). Their analysis is
  heuristic. Lenstra–Pomerance is rigorous but randomized.
- *Missing lemmas.* H1 (smooth-value concentration on S_n) and H2 (the Jacobi character
  nontrivial on the kernel).
- *Experiment (E4).* §9.3.
- *Refuting behaviour.* Sign bias, or excess all-trivial kernels. Not observed.

**Card B16-4: class-group hybrid (escape 3).**
- *Mechanism.* Pair (Z/N)^* targets with Cl(−4kN) collisions.
- *Status.*
  - As a filter: **refuted** (coverage cannot rise).
  - As class-number search: N^{1/5} on average (heuristic variance argument), and generic
    collisions cost N^{1/4}.
  - As a hybrid: it needs a public identity C^{h(p,q)} = D with cheap h. **Open** (O-3); no
    candidate is known.

**Card B16-5: global selector / per-target lattice solving (escape 5).**
- *Mechanism.* Localise p cheaply, then run DCL on the reduced span. Or solve each target by
  Coppersmith instead of BSGS.
- *Status.*
  - Selector: N^{1/6−η/3} if a span-N^{1/2−η} selector exists (**proved** calculation). No
    selector is known (O-5).
  - Lattice: equivalent to known p mod M, threshold N^{1/4}. **Not an escape** (heuristic).


## 11. Reviews

- **GPT-6.1-sol below16_1 (blind).** It reached the same verdict independently and supplied
  T3, Cor. 2 (= Cor. F), Theorem 1 (worst-case short bases), the locality counterexample and
  the x_j = j trap. These are used above with attribution.
- **GPT-6.1-sol below16_2 (adversarial review; `codex/below16_2.out`).** Verdict "major
  revision"; T1's counting and T2 survive, and no counterexample to (2.2) was found in 1 600 T1
  and 240 T2 cases. Adjudication: all 18 findings accepted. Changes:
  - (2) A = 0 step: counterexample at n = 20; repaired by the side swap (§2.1 step 5).
  - (3) (L1) holds after canonical normalisation; the budgets are moment budgets; combined
    residue queries are not covered.
  - (4) T2 box conditions stated.
  - (5) T3's full hypotheses restored; X/P term added to Cor. F.
  - (6) The D′ "factor 2" claim was wrong; repaired (same as the critic).
  - (7) MB restricted to exact rational relations, with a nonconstancy requirement.
  - (8) Giant multiplicity is a model hypothesis, for one base too.
  - (9) Jacobi details fixed.
  - (10) §6.3 arithmetic fixed. I kept c < 1/(6ω) because deterministic O(n^ω) elimination
    exists (Ibarra–Moran–Hui), and added c < 1/18 for cubic elimination.
  - (11) The exact-sieve requirement is added to O-4.
  - (12) The class-number argument is relabelled heuristic and not a lower bound.
  - (13) The accidental-collision uniformity assumption is made explicit; "not an exponent
    source" removed.
  - (14) The lattice "information-equivalence" claim was wrong; the item is now open.
  - (15) "Only live escape" → "one concrete conditional escape"; smooth-order route added;
    "average ≤ sup" no longer used as the explanation.
  - (16) O-5 restated with intervals; the "N mod m carries no information" claim removed.
  - (17) Ward–Maestro 2026 and Brown 2010 cited (not read).
  - (18) The harness's int64 limitation is noted in §12.
- **Independent Opus critic (read-only).** It ran ≈ 600 greedy T1 designs at n = 16–22 and 120
  T2 designs at n = 16–20, with no violation (worst count/closed 0.038 and 0.016). Verdict
  **ACCEPT WITH FIXES**; all 11 fixes applied:
  1. c₁ covers A = 0;
  2. average over uniform n-bit N, not S_n;
  3. hypotheses stated in §0; (L1) put in scale-invariant form;
  4. generator split per class; stationary-point remark moved to T2/T3;
  5. §2.4 counterexample (2.3–3.7×, up to 29×) and the m̄_k repair; D′ is within N^{5ε/6}, D″
     within polylog;
  6. T2 allows any partition with cells ≥ PK;
  7. Cor. F: X/P term, k for support size, slack defined, conclusion conditional;
  8. Prop. MB: missing validity step supplied, verdict restricted to residue-local boxes;
  9. §6.1 overlap/sign/gcd fixes and the Blum fraction;
  10. §6.3 arithmetic;
  11. scope qualifiers in §0 and the table.
- **Status after review.** T1 and T2 are **proved and reviewed** (both reviewers checked them
  line by line and numerically) under the hypotheses now stated. T3 and Cor. F have their
  hypotheses as stated. Prop. MB is proved for exact rational relations and residue-local
  boxes.

## 12. Reproduction

```
.venv/bin/python -m pytest -q tests/test_below16.py
PYTHONPATH=. .venv/bin/python -m src.below16.experiments e1
PYTHONPATH=. .venv/bin/python -m src.below16.experiments e2
PYTHONPATH=. .venv/bin/python -m src.below16.experiments e4
PYTHONPATH=. .venv/bin/python -m src.below16.experiments summarize
```

Harness limitation (codex r2): `avgbound.py` evaluates windows in NumPy int64, so huge
succinct-power coefficients need range guards. The cardinality proof does not depend on the
harness.
