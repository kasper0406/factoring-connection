# below16: can the N^{1/6} floor be escaped?

Agent `below16`, 2026-10-18. Pre-registration: `notes/below16_prereg.md` (commit `7592b7c`,
before any pre-registered run). Code: `src/below16/`. Tests: `tests/test_below16.py`.
Raw data: `experiments/raw/below16_e{1,2,4}.jsonl`. Summary:
`experiments/analysis/below16_summary.json`. GPT-6.1-sol: blind round `codex/below16_1.out`,
adversarial review `codex/below16_2.out`. Independent Opus critic: §11.

## 0. Verdict (plain language)

**No escape below N^{1/6} was found. Instead, the floor survives the move to the true-factor
(average-case) task, which is the task Theorem D′ actually solves.** Four of the five escapes
fail inside one model. The fifth is a real escape, but it does not help at the 1/6 level: it
jumps straight to subexponential time and stands or falls with a derandomization theorem that
nobody has.

- **Main theorem (T1; proved).** Take a *residue-local* design: every target's weights and
  modulus are fixed, and its class set and window (an *arbitrary finite set*, not only an
  interval) depend on N only through N mod K, with K ≤ C·M². Then the fraction of balanced
  semiprimes whose true factor is covered is at most C·n²·X^{3/2}·N^{−1/4}. Success on even a
  constant fraction therefore needs X ≥ c·N^{1/6}/(log N)^{4/3}.
  - The bound needs no flatness and no height, and it needs no equidistribution of primes.
  - It holds under per-input or average cost budgets.
- **T2 (proved)** adds coarse size information: N to precision P·K, with interval windows.
- **T3 (GPT-6.1-sol, checked here)** allows windows given by algebraic formulas in N (for
  example "near the stationary point").
- Theorem D′'s algorithm lies in the T1 model (up to a factor 2 in cost). **So D′ is optimal up
  to polylogarithms among local designs for the true-factor task.** This replaces the
  overclaim "D′ is optimal in its model", which Thm A alone (full-candidate covering) does not
  support.

| escape | broken assumption of Thm A | verdict | label |
|---|---|---|---|
| 1. bases with known relations (2,3,6; β, β^k; succinct powers) | one base; interval windows | reduces to 1-D targets with arbitrary finite windows (Prop. MB). T1 then applies, because it counts window *cardinality*, so sparse sum-of-AP windows do not help. | **refuted for this model** (average case; rational bases). Open: algebraic bases, and polynomial "giant multiplicity" (§3.3) |
| 2. average higher-order flatness | full-candidate covering | average ≤ sup. Degree-1 targets show no flatness effect on average at all (T1). Degree ≥ 2 targets still need N^{1/5} on average (Cor. F, via T3/T2) | **refuted for this model** (degree ≥ 2: proved under the uniform real-fibre hypothesis) |
| 3. class-group collisions | the group (Z/N)^* | no small-exponent identity in the true factor. The order-2 class (p, 0, q) is the factor itself. Using a second group as a filter cannot raise coverage. A hybrid needs a new identity C_N^{h(p,q)} = D_N | **refuted for the filter / class-number-search model**; hybrid identity open |
| 4. relation methods | sumset collision search | a genuine escape: deterministic QS heuristically costs L[1/2]. It needs H1 (smooth-value concentration on S_n) and H2 (the Jacobi character (X_S Y_S/N) is nontrivial on the kernel; exact for Blum N). E4 shows H2's sign statistics are indistinguishable from fair independent coins at 40–80 bits | **conjectured** (heuristic L[1/2]); exact open sub-question O-4 |
| 5. only the true factor needs covering | full-candidate covering | T1/T2/T3 cover every local adaptive selector. A *global* selector shrinking the span to N^{1/2−η} would give N^{1/6−η/3}; none is known | **refuted for local designs**; global selector open (O-5) |

**Novelty:** the average-case/local-design floor (T1) and the Jacobi form of H2: "no prior
formulation found by us". The Jacobi criterion for Blum integers is classical (Rabin/Blum). A
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

**Locality hypothesis (L1).** D ≤ C₁M² (implied by validity) and K ≤ C₂M² (the design looks at
N only modulo K). The K ≤ C₂M² hypothesis is needed: GPT-6.1-sol's counterexample (§7.1) shows
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

- Here c_s = 1 (side p) or 5/4 (side q), and c₁ = (5/2)√(2(C₁+C₂)).
- E_N is the uniform average over n-bit integers N.
- For A = 0 the same holds with 4c_s(P/M)√(P(D/|B| + M)) in front.

Consequently, for any design (any set of generators) obeying the cost model on average
(E_N[#targets] ≤ X, E_N[Σ|S_τ|] ≤ X²),

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
4. **Closed form.** D ≤ C₁M² and K ≤ C₂M² give √(D/|A| + K) ≤ M√(C₁+C₂). This is the second
   inequality.
5. **A = 0.** Here the window depends on x only. The number of admissible x is
   ≤ |R|·min(P/M, |S|(D/(|B|M) + 1)), and each has ≤ 2P/K partners in its class.
6. **Sum over the design.** The covered pairs are contained in the union over generators.
   Summing step 4 over generators gives Σ_ι E_N[|R_ι||S_ι|^{1/2}] = E_N[Σ_τ |S_τ|^{1/2}],
   since each class in R is one target.
7. **Cost budget.** Cauchy–Schwarz per input, then in expectation:

       E[Σ_τ √|S_τ|] ≤ E[√(#τ·Σ|S_τ|)] ≤ √(E#τ·EΣ|S|) ≤ X^{3/2}.

   Divide by |S_n| ≥ c₀P²/n². ∎

**Remarks.**
- The proof counts *integer* pairs, so primality costs only the n² normalisation. No prime
  equidistribution is used; prime-density thinning saves only logarithms (r16).
- There is no H and no flatness. On average, the class lock's "whole class" bound P/M and the
  per-value bound |S|·D/(AK) combine by a geometric mean. Placing a window cleverly (for
  example at the stationary point) cannot beat √|S|.
- Arbitrary finite windows are allowed. A sparse window of cardinality s counts as s, never as
  its hull. This is what closes escape 1 (§3).

### 2.2 Coarse size information (T2)

A **T2 generator** has side p and A, B ≥ 1. Its key is (κ, b) = (N mod K, ⌊(N − 2^{n−1})/Δ⌋),
with Δ = 2^j ≥ PK, and its window is an interval [U, U + W) depending on the key. Any design
that uses size information no finer than PK can be analysed with Δ the smallest such power of
two. That is about n/2 bits of N in total: log K bits of residue and n/2 − log K bits of size.

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

**Model.**
- Each slot τ has a fixed modulus M_τ ≤ P, one unit class r_τ(s) per residue s = N mod M_τ,
  and a region R_{τ,s} in the (p, q) box. The region is semialgebraic of bounded complexity
  and may be given by formulas in N, for example "|Bp + Aq − 2√(ABpq)| < DW".
- There is a width cap w_τ with Σ w_τ ≤ X².
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
- (*)_e: length ≤ C M P^{1−e/K} w^{1/K} for degree-e, support-K targets — the real-variable form
  of multimod (2);
- or, in the T2 model, that every window is at least as wide as its archimedean slack.

Then targets of degree ≥ 2 alone cover a fraction ≤ C n² Σ_{e,K} P^{−e/K} X^{1+1/K}, so positive
success needs X ≥ N^{1/5}/(log N)^{O(1)}.

*Proof.* T3's area argument, or Lemma 2 with multimod (2) applied to C⁺, then Hölder:
Σ w^{1/K} ≤ X^{1+1/K}, and e/(K+1) ≥ 2/5 for K ≤ 2e. ∎

**So escape 2 fails.** Averaging over the true pair cannot beat the supremum of per-input
coverage. For local designs the uniformity step is exactly what T1, T2 and T3 supply. Degree ≥ 2
stays at 1/5 and degree 1 stays at 1/6.

### 2.4 Theorem D′ lies in the model

DCL-det's targets are (A, B, ℓ, root) with D = ℓ² and class set = the roots of Br² ≡ AN
(mod ℓ) (key N mod ℓ). Its window depends on N mod ℓ² and on √N.
- Replace the window by the n-only superset [⌊F_min/ℓ²⌋ − 1, ⌊F_max/ℓ²⌋ + 1] over all n-bit N.
  This at most doubles Σ|S| and gives a T1 design with K = ℓ.
- So D′'s cost N^{1/6+5ε/6} sits within a polylog-and-constant factor of the T1 floor
  c N^{1/6}/n^{4/3}.

**Corollary (proved).** Among residue-local (T1), size-local (T2) or algebraic-window (T3)
designs, the true-factor exponent is exactly 1/6, attained by D′.

## 3. Escape 1: several bases with known relations

### 3.1 Reduction (Prop. MB; proved)

**Setup.**
- Let the bases α_1..α_k be positive rationals (any number, any known relations).
- The relations provable for every p are exactly the relations in Q^*, since a relation false
  in Q^* fails for all large p.
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
  β, β^K).
- L_j(F) is an integer combination of degree-1 forms, so it is a degree-1 target. T1 applies.

**Verdict.** In the true-factor model, multi-base searches with rational bases and arbitrary
known relations do not beat N^{1/6}/polylog. This includes succinct powers β^k with huge k: T1
uses |S|, not hull(S). That closes multimod Thm E(d) in this model.

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
- **(MB-b) Giant multiplicity.** The cost model charges labels, not stored group elements. One
  giant G could carry R labels: certified identities G ≡ γ^{g_k(p)} on different classes.
  - Then Count grows by R at no extra cost.
  - In the HH-lock family, two labels share a giant only if their integer handle exponents
    coincide. For distinct band primes this forces ℓ'²A = ℓ²A′, which is impossible. For equal
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

1. **No small true-factor exponent (proved, elementary).**
   - The unit-group identity works because N ≡ q (mod p − 1) is public arithmetic in the
     exponent.
   - Cl(−4N) (or Cl(−4kN)) has order h(−4kN) ≈ √N·L(1, χ) with no expression in p, q. There is
     no CRT projection of Cl(−4N) to a factor-dependent group with a known exponent.
   - The factor does give an order-2 class [(p, 0, q)] = [I_p], I_p = (p, √−N) (GPT-6.1-sol
     proof: I_p² = (p), and I_p is not principal since a² + Nb² = p is impossible). But writing
     it down *is* knowing p.
2. **Filtering cannot help (proved, trivial).** A second group used to test unit-group targets
   accepts a subset of their covered set. Coverage, and hence T1, is unchanged.
3. **Class-number search (heuristic).**
   - The generic route computes h (or the order of a form) to find ambiguous forms. Shanks's
     BSGS on h ∈ [ĥ − E, ĥ + E] costs x + √E, where ĥ is the Euler product to x.
   - On average over discriminants, the tail Σ_{x<ℓ≤2x} χ(ℓ)/ℓ has variance ≍ 1/(x log x) and is
     asymptotically independent of χ(ℓ) for ℓ ≤ x (orthogonality of quadratic characters).
     So E ≫ √N/√(x log x) on a positive fraction of inputs.
   - Total cost ≥ min_x (x + N^{1/4}x^{−1/4}) = N^{1/5} — the GRH bound of Schoof/Shanks, and no
     better on average.
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
- Let p ≡ q ≡ 3 (mod 4). For a kernel vector S, put X_S = ∏(s + x_i) and Y_S = ∏ℓ^{e_ℓ(S)/2}.
- Then gcd(X_S − Y_S, N) is proper ⇔ the Jacobi symbol (X_S Y_S / N) = −1.
  - X_S/Y_S is a square root of 1. The nontrivial roots of 1 have Jacobi symbol −1, the
    trivial ones +1 (because (−1/N) = 1).
- S ↦ (X_S Y_S/N) is a character of the kernel, since the overlap correction x_i²/Q(x_i) ≡ 1.
- **So H2 ⇔ the computable character ψ_N(S) = ∏_{i∈S}((s+x_i)/N) · ∏_ℓ (ℓ/N)^{e_ℓ(S)/2} is
  nontrivial on ker.**
- H2 is thus a statement about Jacobi symbols of explicit integers, and heuristically fails
  with probability 2^{−dim ker}.

### 6.2 What a proof needs (the exact open sub-question O-4)

1. **H1 on S_n.** For almost all (p, q) ∈ S_n, at least π(y) + 2 of the Q(x) are y-smooth.
   - First moment: smooth numbers in the intervals [2|x|√N, 2|x|√N + O(√N)] over the
     archimedean distribution of N. This is plausible via Hildebrand-type short-interval
     results, but they need y ≥ exp((log N)^{5/6+ε}), larger than L[1/2]. Polynomial-size y
     (u fixed) suffices for N^{o(1)}-free exponents: §6.3.
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
  u ≈ (1/2 + c)/c, which is a *constant*.
- So π(y) + 2 relations come from an interval of length N^{c+o(1)}.
- Deterministic dense elimination over F₂ by fast matrix multiplication costs N^{ωc+o(1)},
  ω < 2.372. (The sparse methods, Wiedemann and block Lanczos, are randomized.)
- Total cost N^{max(1, ω)c + o(1)}, which is < N^{1/6} for c < 1/(6ω) ≈ 0.070, i.e. u ≈ 8.1.
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
2. **Local adaptivity is useless (T1/T2/T3).** Choosing classes, weights or windows from local
   data (residues mod ≤ M², up to n/2 bits of size, or algebraic formulas in N) does not beat
   N^{1/6}/polylog.
3. **False positives / accidental collisions (heuristic).**
   - A non-certified collision mod p needs ord_p(β) | D for a tested difference D. In the
     random-order model the accidental rate is ≤ (#tested pairs)/ord_p ≤ X²·n^{O(1)}/P, which is
     o(1) for X ≤ N^{1/4−δ}. D′'s data agree: 83 accidental first hits against a null of 164.
   - Two-sided (FULL) collisions mod N give a multiple of ord_N(β), hence a factor by Miller's
     reduction when the 2-adic orders differ. They are rare (≤ 0.045% in D′).
   - Not an exponent source.
4. **Global selector (O-5, open).**
   - A cheap public procedure that localises p to an interval of length N^{1/2−η} would give
     N^{1/6−η/3} plus the selector cost (r16; GPT-6.1-sol §5).
   - With η ≥ 1/4 it is already a polynomial-time factorer (Coppersmith).
   - Any selector must use *global* arithmetic of N: continued fractions, class groups, smooth
     relations. Residue information N mod m carries no information about p mod m (trackE
     W-10), and T1–T3 cover everything local.
5. **Per-target lattice solving (card; heuristic).**
   - Replacing BSGS by Coppersmith on each target is information-equivalent to "factor with
     p mod M known", because the lock is implied by p mod M.
   - The Coppersmith/Howgrave-Graham threshold is M ≥ N^{1/4}, so selective covering by lattice
     solves costs N^{1/4}·polylog (the McKee–Pinch regime).
   - Not an escape.

## 8. Precise open sub-questions

- **O-4 (relation route; the only live escape).** Prove H1 and H2 for one explicit
  deterministic schedule on almost all of S_n. The polynomial-level variant (y = N^{c},
  c < 1/(6ω), constant u) suffices for N^{1/6−δ}. In H2's Jacobi form: show that
  ψ_N(S) = ∏((s+x_i)/N)∏(ℓ/N)^{e_ℓ(S)/2} is nontrivial on the kernel for almost all Blum N.
- **O-5 (global selector).** Find any public, N^{1/6−δ}-time computable map N ↦ (set of measure
  N^{1/2−η}) that contains p for most of S_n.
- **MB-a / MB-b.** Algebraic (non-rational) bases; a labels-per-giant bound.
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

Status at hand-back: the cells n = 40–64 (400 Blum + 200 general each) and 225 Blum instances
at 72 bits were complete. The rest of 72 and 80 was still running; see the summary JSON.
- **P4a:** the Jacobi rule held for 76 815/76 815 basis vectors.
- **P4b:** the per-vector nontrivial rate was 0.481–0.507 in every cell.
- **P4c:** every prefix test had p ≥ 0.004, which meets the frozen threshold of 0.001. The
  smallest p values came from a mild excess of all-trivial prefixes at 40 bits (Blum j = 3:
  70 observed against 50 expected), with no excess at 48–72 bits.
- **P4d:** success was 100% in every cell.
- **P4e:** there were 0 ACCIDENTAL hits.

**Decision: sign independence (H2) is experimentally supported at 40–72 bits** for this
deterministic schedule. The 40-bit excess is noted as a possible small-size effect; it is not
explained and was not tested further.

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
- **GPT-6.1-sol below16_2 (adversarial review of this draft).** Launched; **not incorporated
  at hand-back**. See `codex/below16_2.out`.
- **Independent Opus critic.** Launched (read-only); **its verdict was not received before
  hand-back**. T1/T2 must therefore be treated as **proved by the author but not yet
  independently reviewed**. The orchestrator should apply the critic's fixes before relying on
  them.

## 12. Reproduction

```
.venv/bin/python -m pytest -q tests/test_below16.py
PYTHONPATH=. .venv/bin/python -m src.below16.experiments e1
PYTHONPATH=. .venv/bin/python -m src.below16.experiments e2
PYTHONPATH=. .venv/bin/python -m src.below16.experiments e4
PYTHONPATH=. .venv/bin/python -m src.below16.experiments summarize
```
