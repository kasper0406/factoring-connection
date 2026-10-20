# BLIND COPY for auditRA reviewers: the claims of notes/recurse_results.md (Theorem RA, Corollary P), with all review sections, review verdicts and growing-degree material removed. Line content is otherwise verbatim from recurse_results.md (commit 6b4bfbc).

Everything here rests on Theorem A (`notes/thmA_proof.md`). Theorem A is assembled from internal notes that no
outside referee has checked. A result labelled **proved** below is proved *from the proof of Theorem A*: it
inherits every caveat of `notes/thmA_audit.md`, and adds no new external input.


## Claim summary (from the note's verdict section)
**Fixed degree: the 1/√2 floor is gone.**
- Lenstra–Pomerance was the shape test, factoring the halves u₀, v₀ of size ≈ √N. Replacing it by a recursive call of
  the algorithm itself is rigorous, given Theorem A's proof.
- **Theorem RA [proved from Theorem A].** For each fixed (d, k, μ, η₀) admissible for Theorem A (θ < 1), and each
  fixed a > 0, the self-recursive algorithm completely factors every N. Its expected time is L_N[1/2, c(a) + o(1)],
  with c(a) = max(a + θ/(2a), (2 + √2)a).
  - The minimum is at a* = √((√2 − 1)θ/2).
  - The minimum value is **c\* = √((1 + √2)θ)**.
- Constants (infimum over fixed k and η₀; exact in `src/recurse/constants.py`):
  | d | 7 | 9 | 11 | 13 | 19 | 29 | 49 | 99 |
  |---|---|---|---|---|---|---|---|---|
  | c\* | 0.9515 | **0.8510** | 0.7769 | 0.7193 | 0.6018 | 0.4914 | 0.3806 | 0.2691 |
  | Theorem A | 0.9723 | 0.9192 | 0.8839 | 0.8586 | 0.8132 | 0.7778 | 0.7495 | 0.7283 |
  - As d → ∞, c\* → 0.
  - The threshold for c < 1 is unchanged: c\* < 1 ⟺ θ < √2 − 1 ⟺ (1 + θ)/√2 < 1. So the minimal k is the same as
    before (k ≥ 13 at d = 9, where c\* = 0.99521 instead of 0.99720).
- **Corollary P (portfolio) [proved from Theorem A; the argument is GPT-6.1-sol's r19, checked here].** One Las Vegas
  algorithm factors every N in expected time exp(o((log N log log N)^{1/2})) = L_N[1/2, o(1)]. The rate is not
  effective.
  - The same device makes Theorem A's 0.9192 at d = 9 a genuine "+o(1)" (thmA audit caveat 1).
- **What binds.** The shape test was the only place 1/√2 entered.
  - With recursion, the recursive shape test and the special-q search T_SQ bind together.
  - Linear algebra is slack: 2a\*/c\* = 2 − √2 ≈ 0.586.
  - The split probability costs L^{o(1)} and does not bind.
- **(T)** survives. The recursive calls see only u₀ and v₀, which are functions of R, and they use fresh coins

## 1. Fixed degree

### 1.1 The algorithm FACTOR_rec[a]

Fix (d, k, μ, η₀) as in `thmA_proof.md` §1 with θ < 1, and fix a > 0. FACTOR_rec[a] is FACTOR of `thmA_proof.md` §3,
with two changes:
1. B := ⌊L_N[1/2, a]⌋. Theorem A takes a = θ/√2; here a is a free parameter. S² := ⌈L^{θ/(2a)}⌉·⌈log N⌉^{C₃}, as in
   §1, U5 there.
2. In DESCEND_ONCE, the line "factor |u₀| and |v₀| completely with Lenstra–Pomerance" becomes "factor |u₀| and |v₀|
   completely with FACTOR_rec[a], using fresh coins".

FACTOR_rec[a] returns AKS-certified primes, so its output is always correct. Every other line is unchanged, including
the AKS certification of the shape-test output, which is now redundant.

### 1.2 What is used from Theorem A's proof

Write SPLIT^O for SPLIT whose shape test is a procedure O. Here O, given m ≥ 1, always returns the prime factorisation
of m, terminates almost surely, and uses coins independent of everything else. Theorem A's proof, read with a general
a > 0, gives the following facts for every preprocessed N ≥ N₀(d, k, μ, η₀, a). The bounds are finite for every such
N and have the stated asymptotics.

- **(A1)** E[τ(N)] ≤ L_N^{a + o(1)}, where τ is the number of DESCEND_ONCE calls in SPLIT^O. Each call makes at most
  two O-calls, on integers ≤ √(2N).
  - Source: Steps 6, 7 and 10.
  - M_rows = B·L^{o(1)} rows; E[1/q_A | kept] ≤ 1/μ_sel = L^{o(1)}; at most 2 expected LV6 rounds, each with fresh
    rows; field selection L^{o(1)}; absorbing factor exits.
- **(A2)** E[non-O work of SPLIT^O] ≤ L_N^{max(2a, a + θ/(2a)) + o(1)}.
  - Special-q searches: at most 2k·E[τ] = L^{a+o(1)} in expectation. Charge up to 2k searches to every
    DESCEND_ONCE attempt; this simplification is from codex recurse_2, and it makes K(iii)'s waste ratio unnecessary
    here. The absorbing factor-exit bookkeeping of thmA Step 6 integrates q_A away: ∫(q_A/p)(C/q_A) ≤ C/p.
  - Each search costs Õ(B + S²) = L^{max(a, θ/(2a)) + o(1)} (Prop BS).
  - LV6 costs L^{2a+o(1)}.
  - Gauss reduction, gcd checks and certificates are polylog per trial.
- **(A3)** SPLIT^O returns a proper factor almost surely.

**Is a really free?** The orchestrator, relaying r19, flagged this as essential. I checked the source statements
myself:
- FM: `pcest_results.md` l. 304, "Fix odd d ≥ 3, k > (d+1)/4, **a > 0**, a compact σ-range…".
- LD*: `ldbound_results.md` l. 127, the same list.
- Lemma 8.5′: `hypD_results.md` l. 396, B = L[1/2, a] with a generic.
- LV6: `lv6_results.md` l. 88, B = L[1/2, a].
- Theorem R: `pcest_results.md` l. 478.

The value a = θ/√2 enters only at the optimisation step: pcest l. 762, ldbound l. 485, lv6 l. 645, thmA_proof §1 and
Step 10. The one side remark that uses it is Prop BS's "B ≤ T_SQ". That remark concerns cost only, and it is replaced
above by Õ(B + S²) = L^{max(a, θ/(2a))}.
- For every a, the σ-range for S = L^{θ/(4a)} is a fixed point.
- The condition S² ≥ ρ_r^{−2}L^{c} holds for any fixed 0 < c < (h + (d+1)ν)/(2a), because
  log(S²ρ_r²) = ((ξ_F − ξ_r)/(2a) + o(1))Λ with ξ_F − ξ_r = h + (d+1)ν > 0. The margin is a power of L, not of N; the
  first draft's justification was wrong (codex recurse_2), but the conclusion stands. This condition is used only by
  Theorem R, which is off the critical path.
- Lemma Ω needs 2^ω = L^{o(1)}, and ω ≤ log N/log B for every a > 0.

So (A1)–(A3) hold for every fixed a > 0, with N₀ depending on a.

### 1.3 Proof of Theorem RA

**Lemma OS (oracle substitution) [proved].** Let every admissible O return the factorisation in a *canonical form*:
the sorted list of primes with multiplicities. Fix added after the Opus proof critic: otherwise the order of the
loop "for each ℓ | u₀v₀" could depend on O's coins. Under this convention, the joint law of the outer transcript of
SPLIT^O is the same for every admissible O. In particular this holds for τ, the non-O work and the output. The outer transcript is the *logical* one: the
outer coins and every value the mathematical algorithm computes outside O. It excludes timing instrumentation, which
is a wording fix after codex recurse_2: the toy code records elapsed seconds, and its comparison excludes them. The
algorithm has no timeout, race or branch that depends on O's running time or on O's coins.

*Proof.*
- The outer computation is a deterministic function of the outer coins and of the values O returns. Every call on m
  returns the prime factorisation of m, a function of m alone by unique factorisation.
- Induction over the steps of the outer computation shows that every outer value is a function Φ of the outer coins
  only, with Φ independent of O. ∎

Hence (A1)–(A3), proved in Theorem A for O = Lenstra–Pomerance, hold for O = FACTOR_rec[a].

**Lemma W (Wald for adaptive sub-calls) [proved; standard].**
- Setting: a randomised procedure makes sub-calls 1, 2, …. The decision to make call i, and its input m_i, are
  measurable with respect to 𝓕_{i−1}, the history before call i. Call i uses fresh coins and costs C_i, with
  E[C_i | 𝓕_{i−1}, m_i] ≤ F̄ almost surely.
- Conclusion: E[Σ_{i≤τ} C_i] ≤ F̄·E[τ].

*Proof.* E Σ_i C_i 1{τ ≥ i} = Σ_i E[1{τ ≥ i} E[C_i | 𝓕_{i−1}, m_i]] ≤ F̄ Σ_i Pr[τ ≥ i], by monotone convergence. ∎

**The recursion inequality.** Put F(x) := sup_{2≤n≤x} E[cost of FACTOR_rec[a](n)]. F is nondecreasing.
- **Finiteness.** F(x) < ∞ for every x, by strong induction on n.
  - For n < N₀ the cost is O(1).
  - For n ≥ N₀ the inequality below holds with finite P(n) and D(n), and it calls F only at √(2n) < n.
- **Accounting.** FACTOR(N) builds a factor tree with at most 2Ω(N) − 1 ≤ 2 log₂N nodes.
  - Each node n | N costs O(L_n^{a+o(1)}) for AKS, perfect-power detection and trial division to B(n), plus one
    SPLIT^O(n′) on its cofactor.
  - Node values are random, but each node's conditional expected cost given its value is at most the supremum over
    n ≤ N.
  - The O-calls inside SPLIT^O(n′) are FACTOR_rec[a] on inputs ≤ √(2n′) ≤ √(2N), with fresh coins. Their conditional
    expected cost is ≤ F(√(2N)).
- By Lemma W and (A1)–(A2), for x ≥ x₀:
  **F(x) ≤ 2 log₂x·[P(x) + 2D(x)·F(√(2x))],** with P(x) = L_x^{max(2a, a+θ/(2a)) + o(1)} and D(x) = L_x^{a+o(1)}.
  - The suprema over n ≤ x of bounds of the form L_n^{p+ε(n)}, with ε(n) → 0, are again L_x^{p+o(1)}.

**Lemma F (solving the recursion) [proved].** Throughout, write L_x := exp(√(log x·max(1, log log x))). This is the
usual L_x for x ≥ e^e, and it is ≥ 1 everywhere. Fix after the critic: L_x is undefined for x < e.
- Hypotheses: F: [2, ∞) → [0, ∞) is nondecreasing and finite, and for x ≥ x₀,
  F(x) ≤ L_x^{p+ε(x)} + L_x^{δ+ε(x)}F(√(2x)), with ε(x) → 0.
- Conclusion: F(x) ≤ L_x^{c+o(1)}, where c = max(p, δ/(1 − 1/√2)) = max(p, (2 + √2)δ).

*Proof.*
1. **Scaling.** Put y = √(2x). Then log y·log log y ≤ ½(log x + log 2)log log x, so
   L_y^{b} ≤ L_x^{b/√2}·exp(b·O(√(log log x/log x))) ≤ e·L_x^{b/√2} for x ≥ x₁(b).
2. **Margin.** Fix η > 0 and put η′ := η(1 − 1/√2)/4. Choose X ≥ max(x₀, x₁(c+η), 4) such that for x ≥ X:
   ε(x) ≤ η′, L_x^{η−η′} ≥ 2 and L_x^{3η′} ≥ 2e. Put K := max(1, F(X)).
3. **Base.** F(x) ≤ K ≤ K·L_x^{c+η} on [2, X].
4. **Step.** Suppose the bound holds on [2, Y] with Y ≥ X, and let x ∈ (Y, Y²/2]. Then √(2x) ≤ Y and
   F(x) ≤ L_x^{p+η′} + K e L_x^{δ+η′+(c+η)/√2}.
   - Since δ ≤ c(1 − 1/√2), the second exponent is ≤ c + η/√2 + η′ = c + η − 3η′. That term is ≤ (K/2)L_x^{c+η}.
   - The first term is ≤ L_x^{c+η}/2.
   - So the bound extends to [2, Y²/2], and Y²/2 > Y because Y ≥ 4.
5. Hence F(x) ≤ K_η L_x^{c+η} for all x, and since η was arbitrary, F(x) ≤ L_x^{c+o(1)}. ∎

**Theorem RA [proved from Theorem A].** For fixed admissible (d, k, μ, η₀) with θ < 1 and fixed a > 0, FACTOR_rec[a]
is Las Vegas and completely factors every N ≥ 2 in expected time
  L_N[1/2, c(a) + o(1)],  c(a) := max(a + θ/(2a), (2 + √2)a),
where the o(1) depends only on (d, k, μ, η₀, a).

*Proof.*
- Apply Lemma F with p = max(2a, a + θ/(2a)) and δ = a. The 2 log₂x factors are L^{o(1)}.
- Then c = max(2a, a + θ/(2a), (2 + √2)a) = c(a), because (2 + √2)a > 2a.
- Correctness: every output prime is AKS-certified, and termination is almost sure because the expectation is
  finite. ∎

**Optimisation [proved; exact].**
- a + θ/(2a) is decreasing on (0, √(θ/2)], and (2 + √2)a is increasing. They cross where (1 + √2)a = θ/(2a), that is
  at a\* = √((√2 − 1)θ/2) < √(θ/2).
- So c\* = min_a c(a) = (2 + √2)a\*, and c\*² = (2 + √2)²(√2 − 1)θ/2 = (1 + √2)θ.
- Every comparison used here reduces to a rational inequality in θ, checked in `constants.py` and the tests:
  - c\* < 1 ⟺ (1 + θ)² < 2;
  - c\* < (1 + θ)/√2 ⟺ 8θ² < (1 + θ²)² ⟺ θ < √2 − 1 (for 0 < θ < 1);
  - 2a\* < c\* always.
- So the recursion helps exactly when Theorem A already beats Lenstra–Pomerance, and never makes things worse there.
  When θ ≥ √2 − 1, use the min of the two.
- a\* is algebraic, so B = ⌊L^{a\*}⌋ is computable to any precision. Alternatively, a rational a slightly off a\* gives
  c(a) − c\* = O(|a − a\*|).

**Fixed-depth recursion** (Lenstra–Pomerance at the bottom).
- c₀ = 1 and c_{j+1} = g(c_j), where g(c) = min_a max(2a, a + c/√2, a + θ/(2a)) = (θ/c + c)/√2 for c² ≥ θ.
- c₁ = (1 + θ)/√2 is Theorem A itself. c\* is the attracting fixed point, with g′(c\*) = (1 − θ/c\*²)/√2 = √2 − 1.
- At d = 9 (θ → 3/10): 1, 0.9192, 0.8808, 0.8637, 0.8563, 0.8532, … → 0.8510.
- The self-recursive FACTOR_rec[a\*] reaches c\* directly, so no depth limit is needed.

### 1.4 Transcript invariance

**Lemma T_rec [proved].** Theorem DD's hypothesis (T) holds for SPLIT^O with O = FACTOR_rec[a].

*Proof.*
- The shape test's input |u₀|, |v₀| is a function of R.
- Its coins are fresh and independent of every r_i.
- Its entire internal transcript, including its own hidden roots, fields, rows and dependencies, is a function of
  (|u₀|, |v₀|, its coins).
- Hence the joint law of (outer exposed transcript, inner transcripts) is a function of (R_i)_i and coins independent
  of (r_i)_i. Conditional on all of it, each r_i is still uniform on r_i·μ₂, independently.
- DD's proof (lemmaM §1, steps 1–4) applies verbatim.
- The inner calls' own DD arguments concern their own hidden roots. These are *not* independent of the outer squares:
  their law depends on the inputs |u₀|, |v₀|. The first draft said "independent"; corrected after codex recurse_2.
  What holds, and what is needed, is conditional invariance: the inner transcripts are unchanged under r_i ↦ εr_i on
  each outer square-root fibre.
- Nothing reads r_i beyond R_i: the inner gcds are with |u₀|, |v₀| and their divisors, not with N. ∎

Lemma OS gives more: the outer transcript does not depend on which exact oracle is used at all. E1 checks this
bit-for-bit (P2).

### 1.5 What binds, and where 1/√2 entered

The cost ledger of `thmA_proof.md` Step 10, with exponents in units of Λ = (log N log log N)^{1/2}:

| term | Theorem A (LP shape test) | recursive, at a\* |
|---|---|---|
| rows × shape test | a + 1/√2 | a + c\*/√2 = **c\*** |
| rows × T_SQ | a + θ/(2a) | **c\*** |
| rows × batch B | 2a | 2a\* = (2 − √2)c\* |
| linear algebra (LV6) | 2a | (2 − √2)c\* |
| trial division, characters | a | a\* |
| split acceptance 1/q_A, field selection | L^{o(1)} | L^{o(1)} |

- The 1/√2 entered only through T_shape = L_{√N}[1/2, 1]. No other term contains it.
- After recursion the binding pair is {shape test, T_SQ}. Linear algebra and batch smoothness are slack by the factor
  2 − √2. The split probability is L^{o(1)} at fixed k.
- So the next improvement must lower θ (better design floors) or the √N size of the split halves. Better linear
  algebra or a better split law would not help.
- **The halves.** Gauss reduction makes |u₀v₀| ≍ N. Any 2-term split has halves of total size N, so the recursive
  cost L_{√N}[1/2, c] is the minimum for 2-term splits. That is a heuristic remark, not a theorem.

### 1.6 Corollary P (portfolio) and the infimum over d

**Corollary P [proved from Theorem RA; argument from r19 item 1, re-checked].**
- The algorithm: let (A_j)_{j≥1} be the algorithms FACTOR_rec[a\*_j] for a computable sequence of admissible
  parameter sets with c\*_j → 0. One such sequence: d_j = 2j + 7; k_j = ⌈d_j²/4⌉ + j; μ_j the exact optimum; η₀
  small and rational; a_j rational near a\*_j. The critic checked admissibility for j ≤ 7; at j = 7, d = 21,
  k = 118, θ = 0.144 and c\* ≈ 0.590.
- **Restart schedule.** This replaces the first draft's live dovetailing, after the critic: interleaving unboundedly
  many live processes has a model-dependent overhead. In phase s = 1, 2, …, run A_i for i ≤ s, each from scratch with
  fresh coins and step budget 2^{s−i}. Stop at the first certified output. Phase s costs ≤ 2^s steps up to a polylog
  bookkeeping factor.
- **The bound.** Fix j and put s₀ := j + ⌈log₂(2E[T_j(N)])⌉.
  - In each phase s ≥ s₀, A_j's run fails to finish within budget 2^{s−j} ≥ 2^{s−s₀}·2E[T_j] with probability
    ≤ 2^{−(s−s₀)−1}, by Markov. Runs in different phases are independent.
  - So E[T(N)] ≤ polylog·Σ_s 2^s·Pr[no output before phase s]
    ≤ polylog·(2^{s₀+1} + Σ_{m≥1} 2^{s₀+m}·2^{−m(m+1)/2}) = O(2^{j}E[T_j(N)]·polylog).
  - Hence limsup_N log E[T(N)]/Λ_N ≤ c_j for every j, so E[T(N)] = exp(o(Λ_N)). ∎
- Remarks:
  - No uniformity in j is needed, so the rate of o(1) is not effective.
  - The same dovetailing over k at fixed d = 9, without recursion, gives L[1/2, 0.9192 + o(1)]. This answers caveat 1
    of `thmA_audit.md` (the k-diagonalisation) without any uniform estimate.
  - It is **not** an L[α] bound for any α < 1/2.
- The shape tests inside A_j could also call the portfolio itself. That is not needed.

**Infimum within fixed-d theorems.** inf_{k,η₀} c\* = √(3(1 + √2)/(d + 1)) → 0. So each fixed-d theorem gives a
positive constant, and their infimum is 0.

