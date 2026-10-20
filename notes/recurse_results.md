# recurse: the recursive shape test (fixed degree) and the growing-degree question

Agent `recurse`, phase 11, 2026-10-20.

Files:
- Pre-registration: `notes/recurse_prereg.md` (commit `ea5953b`).
- Code: `src/recurse/`.
  - `constants.py`: exact constants.
  - `pipeline_rec.py`: toy recursive pipeline (operational).
  - `harness_e1.py` and `analyze_e1.py`.
- Tests: `tests/test_recurse.py`.
- Data: `experiments/raw/recurse_e1.jsonl` and `_summary.json`.
- Codex: `codex/recurse_1.{md,out}` (blind round; the full report is in `notes/recurse_codex1_report.md`) and
  `codex/recurse_2.{md,out}` (adversarial referee).

Everything here rests on Theorem A (`notes/thmA_proof.md`). Theorem A is assembled from internal notes that no
outside referee has checked. A result labelled **proved** below is proved *from the proof of Theorem A*: it
inherits every caveat of `notes/thmA_audit.md`, and adds no new external input.

## Verdict

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
  (Lemma T_rec). The toy run confirms the stronger exact statement, oracle substitution (Lemma OS).
  - In 22/22 instances the outer transcript was bit-identical under the recursive oracle and under Pollard–Brent.
  - All pre-registered predictions held (§2). The DD rate was 0.484 with z = −2.2, inside the frozen |z| ≤ 3.
- **Reviews.** Two GPT-6.1-sol rounds (blind and adversarial) and two independent Opus critics (proofs; literature)
  found no fatal or major error in RA, P or S″. Every fix is applied (§5).

**Growing degree: no rigorous L[α] with α < 1/2 is obtained.**
- Recursion does remove lemmaM's barrier (c), and every other cost barrier from smoothness or shape testing, at every
  level. A shape test on √N-size halves costs T(√N) = L_N[α, 2^{−α}c]. Every intermediate smoothness test costs
  L_N[α]^{o(1)} (§3.1).
- What binds instead is the descent geometry and the precision of the relation-supply estimates.
  - **Size model, one-level descent [heuristic arithmetic; re-derived independently, agreeing with lemmaM (b) and
    r19 item 4].**
    - If the descendability density satisfies γ ≥ (log N)^{−O(1)} uniformly in d and k, the floor is α ≥ 2/5.
    - If instead log(1/γ) ≍ d^τ log log N, the floor is α ≥ (2 + τ)/(5 + 2τ). For τ = 1, which covers losses
      (log N)^{O(d)} or d^{O(d)}, this gives **3/7**.
  - **Multi-level degree-1 descent in one field.**
    - With a bounded number of levels, the size model reaches every α > 1/3, using ⌈(1 − 2α)/(3α − 1)⌉ levels.
    - With O(log log N) levels whose constants contract geometrically, it reaches **α = 1/3 itself**. My first draft
      missed this by counting exponents only; the GPT-6.1-sol blind round recurse_1 §2.5 found it, and I re-checked
      it (§3.4).
    - This needs a *nested* descendability statement: values whose prime factors lie in a field-dependent good set.
      There is no tool for that (§3.4).
- **Lemma S″ [proved, elementary].** Below L[1/2], trial division reaches only L[α], so 2^ω is no longer L^{o(1)}.
  Lemma S″ removes the 2^ω loss of Theorems S/S′: Pr[accepted k-split] ≥ |S|²/(2φ(N)) once |G| ≥ 4k³2^ω (§3.2). This
  is the one growing-degree obstacle that can be closed at once.
- **Inventory of what must be re-proved (§3.5).**
  - FM and LD* at growing (d, k), with polylogarithmic precision whose exponent is uniform in d and k. This is the
    decisive gap, and it requires an explicit-error rewrite of LV §5 and of Kern.
  - Lemma 8.5′ at Y = L[3/5].
  - Theorem K with growing k. Lemma S″ handles the 2^ω in the split count. The waste ratio is avoided by charging
    2k·T_SQ to every attempt.
  - LV6 at growing d. LV's own §6 is at growing degree. Our Lemma X, however, has a Stark (2d)! threshold that defeats
    a polynomial character-bit budget. A plausible repair costs L[α]^{o(1)}, but it is not written (codex recurse_2).
  - A shrinking slack η₀(N), with every η₀-margin re-checked.
  - LV prove only a first moment for unrestricted boxes, at L^{o(1)} precision (their Thm 5.4 with Lemma 5.21). They do not treat special-q lattices,
    second moments or descent.
- **Answer to the task's questions.**
  - Recursion plus one-level descent gives a rigorous L[2/5] *only if* that precision-uniform FM/LD* is proved. It is
    not proved.
  - L[3/7] is the conservative target, if the constants grow like (log N)^{O(d)}.
  - Multi-level descent with recursion at each level reaches L[1/3] in the size model (with O(log log N) levels). A
    rigorous version needs the nested-descendability lemma, for which no tool exists. The closest are LP92 Thm 6.1
    and de Boer–Pellet-Mary–Wesolowski 2025 (ERH), and neither applies (§3.4).

**Prior art (§4).** We found no rigorous analysis that calls the factoring algorithm recursively on auxiliary
numbers, and no fixed-point L-constant equation. This is a search result, not a novelty claim.
- The closest precedent: Pomerance 1987, pp. 138–139, completely factors the auxiliary p − 1 with his own rigorous
  algorithm REC inside a discrete-log algorithm. That is a single call and not a recursion.
- All rigorous L[1/2] analyses use an L^{o(1)} smoothness test (ECM via Pomerance 1987, or HECM) and recurse only on
  factors already found.
- The recursive variant itself was first computed inside this project, as a remark in hypD (Corollary K, "recursive
  variant", from a critic). GPT-6.1-sol r19 item 1 independently re-derived the fixed point and proposed the
  portfolio corollary.

---

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

## 2. Experiment E1 (pre-registered `ea5953b`; toy; correctness only)

Commands:
- `.venv/bin/python -m src.recurse.harness_e1 experiments/configs/recurse_e1.json experiments/raw/recurse_e1.jsonl`
- `.venv/bin/python -m src.recurse.analyze_e1 experiments/raw/recurse_e1.jsonl experiments/raw/recurse_e1_summary.json`

The run: 22 instances, each run with oracle `pb` (Pollard–Brent, private coins) and oracle `rec` (the pipeline
itself). Every pre-registered prediction holds.

| prediction | result |
|---|---|
| P1: correct complete factorisations | pb 22/22, rec 22/22 |
| P2: outer transcript bit-identical (Lemma OS) | **22/22** |
| P3: recursion exercised in ≥ 80% of semi/tri ≥ 50 bits | 14/14 |
| P4: row-identity violations, all depths and both oracles | 0 |
| P5: DD rate, ω = 2, all `rec` splits | 2390/4936 = 0.484, z = −2.22 (threshold \|z\| ≤ 3: holds) |
| P5: inner splits only | 2315/4792 = 0.483, z = −2.34 (holds) |

Descriptive:
- There were 12 596 inner SPLIT calls, all at depth 1, on 25–28-bit products of two primes in (2^{12}, 2^{14}].
  - 597 ended STRUCTURAL (by Theorem DD).
  - 11 999 ended ACCIDENTAL, i.e. by a gcd hit on an object of a tiny modulus. This is expected at 25–28 bits.
- The recorded maximal FACTOR depth is 2. The depth-2 calls were all resolved by primality or trial division, without
  a SPLIT.
- The `rec`/`pb` time ratio is 1.1–37. Pollard–Brent is far cheaper than any L-type algorithm at toy sizes, so the
  ratio has no asymptotic meaning.
- **The z ≈ −2.2 in P5.** It is inside the frozen threshold, but it is the largest deviation seen so far in the
  project's DD checks (thmA E1: −0.82, +0.29). Nearly all of the tested dependencies come from the tiny inner moduli.
  The 8 dependencies per LV6 round are linearly independent, so they are independent fair tests under DD. We see no
  mechanism for a bias and record it as a fluctuation (two-sided p ≈ 0.026).

E1 checks the code and the exact Lemma OS. It cannot test any asymptotic claim.

## 3. Growing degree

Notation: Λ_α := (log N)^α(log log N)^{1−α}, so L_N[α, c] = exp(cΛ_α). D := log N/log log N.

### 3.1 Recursion at L[α]: the shape-test barrier disappears [proved, conditional on the outer algorithm]

Suppose some design has its non-shape work bounded by L_N[α, p + o(1)] and makes L_N[α, δ + o(1)] expected shape
tests on inputs ≤ √(2N), with all the structure of §1.2. Then Lemma F holds with Λ_α in place of Λ.
- The scaling step becomes L_{√(2N)}[α, b] ≤ e·L_N[α, 2^{−α}b].
- The conclusion is c = max(p, δ/(1 − 2^{−α})).

So the cost of complete factorisation of the halves never forces α. This is lemmaM's barrier (c): HECM at
B1 = L[2/3] costs L[4/9]. The recursive replacement costs L_{√N}[α] = L_N[α, 2^{−α}c].

Every intermediate smoothness test, on a value x = exp((log N)^{γ}) with γ < 1, costs
L_x[α, c] = exp(c(log N)^{γα}…) = L_N[α]^{o(1)}. So multi-level designs need neither HECM nor batch smoothness at the
large intermediate bounds.

**Caveat.** This replaces every cost barrier of testing. It does nothing for relation supply or descent geometry.

### 3.2 Lemma S″: split counting without the 2^ω loss [proved]

At L[α], α < 1/2, preprocessing reaches only B = L[α]. So ω(N) ≤ log N/log B can make
2^ω ≈ exp(c(log N)^{1−α}(log log N)^{α−1}) ≫ L[α]. Theorems S and S′ lose exactly 2^ω.

**Lemma S″.** Let N be odd, h := [(Z/N)^× : (Z/N)^{×2}] = 2^{ω(N)}, and k ≥ 1.
- Let G be a set of primes coprime to N.
- Let S be the set of products of exactly k distinct primes of G that lie in a window [Y/2, Y], with 2Y² < N and
  |S| ≥ C(|G|, k)/(2k).
- Assume |G| ≥ 4k³h.
- Let R be uniform on (Z/N)^{×2}, and let E be the event that the primitive vector of the lattice
  {(u, v) : u ≡ Rv (mod N)} is ±(u, v) with u, v ∈ S coprime.
Then Pr[E] ≥ |S|²/(2φ(N)).

*Proof.*
1. **Classes.** For a class c of (Z/N)^× modulo squares, put S_c := S ∩ c. A pair (u, v) ∈ S² has uv^{−1} a square iff
   u and v lie in the same class. Cauchy–Schwarz over the at most h classes gives Σ_c|S_c|² ≥ |S|²/h.
2. **Shared primes.** For u ∈ S, the number of v ∈ S that share a prime with u is
   ≤ k·C(|G| − 1, k − 1) = (k²/|G|)C(|G|, k) ≤ (2k³/|G|)|S|.
3. **Coprime same-class pairs.** Their number is therefore
   P′ ≥ Σ_c|S_c|² − (2k³/|G|)|S|² ≥ |S|²(1/h − 2k³/|G|) ≥ |S|²/(2h).
4. **Injectivity.** Map a coprime same-class pair to R := uv^{−1} mod N ∈ (Z/N)^{×2}. Suppose two such pairs give the
   same R. Then uv′ ≡ u′v and |uv′ − u′v| < Y² < N, so uv′ = u′v, and since both are coprime and positive they are
   equal.
5. **Primitive vector.** For such R, the shortest lattice vector s satisfies |det(s, (u, v))| ≤ ‖s‖·‖(u, v)‖ ≤ 2Y² < N.
   Every determinant of two lattice vectors is a multiple of N, so this one is 0. Then (u, v) = m·s, and m = ±1 since
   (u, v) is primitive. This is Theorem S's step 5.
6. **Count.** Pr[E] ≥ P′/|(Z/N)^{×2}| = P′h/φ(N) ≥ |S|²/(2φ(N)). ∎

Remarks:
- The 2^ω of Theorem S comes from the Cauchy–Schwarz route, whose Σn(R)² bound does not use the same-class
  restriction. Restricting to coprime pairs makes n(R) ∈ {0, 1}. This is the "Dixon pointwise" route that hypD §2
  mentions but did not write.
- Lemma S″ also removes c_k and (1 + log Y), so it is sharper than S′ even at fixed d.
- The hypothesis |G| ≥ 4k³2^ω is mild *for good fields*, those with g_A above the density lower bound; it is not
  automatic for every retained field.
  - At fixed d it holds trivially: |G| = N^{Θ(1)} and 2^ω = L^{o(1)}.
  - At L[α] with J ≈ N^{1/(2k)} and k ≍ D^α: log|G| ≈ (log N)^{1−α}(log log N)^α/(2κ), while
    ω log 2 ≲ (log N)^{1−α}(log log N)^{α−1}/a. They differ by a factor log log N.
- **Theorem K(iii)'s waste ratio.** It had 2^{2ω}. With Lemma S″ one factor 2^ω remains:
  Pr[shape]/Pr[accept] ≤ 4·2^ω(|T_J|/|S|)² (correction after the critic). Removing it would need an upper bound on
  Pr[shape] that does not lose 2^ω, and that needs the square classes of the J-primes to be equidistributed mod N,
  which is far beyond Burgess.
  - At L[α] one simply charges a special-q batch to every attempt. Cost per accepted row ≤ (T_shape + 2k·T_SQ)/Pr[accept],
    with 1/Pr[accept] = exp(O(k(log log N + log(1/g_A)))). That is L[α, O(1)] under the precision of §3.3.

### 3.3 Size model for one-level descent [heuristic arithmetic, exponents of log N only]

Use design (a‴_k; μ) with d, k growing and B = L[α, a]. The requirements are as follows.
- **(R1) Split acceptance, Theorem K with Lemma S″.**
  - Per-attempt success is ≈ (g_A/(4 log N))^{2k}.
  - So k(log log N + log(1/γ)) = O(Λ_α). With γ ≥ (log N)^{−O(1)} this means k = O(D^α).
  - The primes ℓ ≈ N^{1/(2k)} are then L[1 − α]. This is lemmaM (b)'s β₁ = 1 − α.
- **(R2) Smoothness.** u = θ log N/(aΛ_α). At the scale θ ≍ D^{2α−1} this is u ≍ D^α, so log u ≈ α log log N, and
  log(1/(ρ_rρ_F)) ≈ α·θ·log N·log log N/(aΛ_α). The critic used (1 − α) log log N, which is the value of log u for
  fixed θ. Either way, this is O(Λ_α) iff θ = O(D^{2α−1}).
- **(R3) Design floors** (exact in `src/thmA/params.py`).
  - θ = 2μ + h + (d − 1)ν with h ≥ max(2ν, μ − (d − 2)ν, 1 − dμ + ν).
  - The window and level floors give θ ≥ 3/(d + 1). The Kern floor gives θ ≥ (d + 1)ν = (d + 1)/(4k).
  - So d ≳ D^{1−2α} and k ≳ d/θ ≳ D^{2(1−2α)}.
- **(R4) Box.** S^d = exp(O(dΛ_α)) must be ≤ N^{o(θ)} if the one-level norm-size formula N^{θ} is to stay valid.
  This holds iff α > 1/3, and it does not bind at 2/5. At α = 1/3 the term is merely one more L[2/3]-scale summand,
  so it is no barrier in the multi-level model of §3.4.

Combining (R1) and (R3): 2(1 − 2α) ≤ α, that is **α ≥ 2/5**. Equality holds at d ≍ D^{1/5}, k ≍ D^{2/5}, with ℓ = L[3/5].

**The precision variant.** Suppose the proved moments only give log(1/γ) ≍ d^τ·log log N. Then (R1) becomes
k·d^τ = O(D^α), which gives (2 + τ)(1 − 2α) ≤ α, that is **α ≥ (2 + τ)/(5 + 2τ)**:

| loss in γ | τ | α floor |
|---|---|---|
| (log N)^{O(1)} uniform in d, k | 0 | 2/5 |
| (log N)^{O(d)} or d^{O(d)} | 1 | **3/7** |
| (log N)^{O(d²)} | 2 | 4/9 |

I derived this independently, before reading r19 item 4. r19 reaches the same 2/5 and 3/7.

### 3.4 Multi-level descent [heuristic arithmetic]

Use one field of degree d = D^δ. A level-j prime ℓ_j = L[β_j] has special-q values of log-size
≈ (log N)^{1−δ} + (log N)^{δ+β_j} (degree-1 elements a − bα). These must be L[β_{j+1}]-smooth at cost L[α], so:
- 1 − δ ≤ α + β_{j+1} for all j, including the last level, where β = α. This forces δ ≥ 1 − 2α.
- δ + β_j ≤ α + β_{j+1}, so β_{j+1} ≥ β_j − (α − δ) ≥ β_j − (3α − 1).

Starting from β₁ = 1 − α (forced by the split), the descent reaches α after ⌈(1 − 2α)/(3α − 1)⌉ levels whenever
α > 1/3:
- α = 2/5: 1 level;
- α = 0.38: 2 levels;
- α = 0.35: 6 levels.

**α = 1/3 itself (correction after the codex blind round recurse_1 §2.5, re-checked).** The count above tracks only
exponents of log N. At α = 1/3 the decrement is 0 in the exponent, but the *constants* can contract.
- Setup:
  - E := Λ_{1/3}; Q := t^{2/3}(log t)^{1/3} with t = log N, so Q/E = D^{1/3} =: D′.
  - d = δD′, log M ≍ log‖f‖ ≍ Q, B = e^{bE}, search radius S = e^{sE}.
  - Targets q have z = log q ≤ b₁Q.
- At a node with z ≥ bE, require all children to be ≤ e^{λz} for a fixed λ < 1. Then
  log X_F ≲ κQ + (d/2)z + dsE and log X_r ≲ mQ + z/2 + sE, so
  u_F + u_r ≤ [(κ + m)/(λb) + δ/(2λ) + δs/(λb) + o(1)]·D′
  uniformly in the tree.
- Since D′ log D′ ≍ E/3, each node costs L[1/3, O(1)]. Taking b large makes the s-coefficient smaller than the
  available 2s, so a finite s exists.
- The tree:
  - z falls from b₁Q to bE in O(log(Q/E)) = O(log log N) levels.
  - Each node has at most polynomially many children in t.
  - So the tree has exp(O((log log N)²)) = L[1/3]^{o(1)} nodes.
- The top-level values have size N^{δb₁/2+o(1)}. Recursion on them is not circular provided δb₁/2 < 1, and its cost
  contracts by (δb₁/2)^{1/3}.
- **So in the size model, recursion plus O(log log N)-level degree-1 descent in one field of degree ≍ D^{1/3} attains
  L[1/3, c] for a finite, unoptimised c.**
- This is the heuristic NFS-DL picture: β′ = β/2 + 1/6 recursions for higher-dimensional descents (Enge–Gaudry–Thomé
  §2.2, curves; Joux–Lercier–Naccache–Thomé 2009 for degree t > 1 sieving), and constant-contracting 2-D descents.
  The literature critic could not verify which papers use 2-D versus higher-dimensional descent for the medium primes.
- My first draft said "floor 1/3, not attained". That was wrong *in the model*. The rigorous status below is
  unchanged.

**The rigorous obstacle below 2/5 is not cost.** Recursion makes every test cheap. The obstacle is *endogenous
children*.
- At level j ≥ 2, the primes to be descended are the prime factors of special-q values.
- Theorem K's split counting (Theorem S/S″) works for *any* target set, because R is uniform on a group.
- Polynomial values have no such structure. The needed statement is a lower bound for the number of box points whose
  values factor over a field-dependent subset G^{(j+1)} of relative density ≥ (log N)^{−C}, averaged over the field.
- FM and LD* count smoothness over *all* primes ≤ y. Neither we nor the literature critic found a lower bound for
  polynomial values composed of primes from an arbitrary dense subset, and none at all in the endogenous setting.
- Two near-misses should be cited:
  - Lenstra–Pomerance 1992, Thm 6.1 (p. 499), bounds ψ(x, y; 𝒫) from below for integers in an interval whose primes
    lie in an arbitrary set 𝒫, applied to the ECM-recognisable primes. This is the *integer-interval* analogue.
  - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588, 2025), Thms 1.1/9.5, sample under ERH, for an arbitrary family
    of ideals in a number field of any degree, β ∈ 𝔞 with (β)𝔞^{−1} in that family times B-smooth ideals, with
    probability ≳ δ/3. Thm 17.7 there is a rigorous (ERH) one-step descent at L(1/2) cost.
  - Neither closes the gap. The first is for intervals, not polynomial values in a special-q lattice. The second
    needs ERH, is one-sided (ideal classes, not the two-sided a − bα / a − bm values with a shared rational side), and
    has one level. dBPW is still the closest existing tool for a "field-dependent good set", and a natural place to
    look for an (ERH-conditional) nested lemma.
- r19 item 4 makes the same point ("whole descent trees … Theorem K … does not automatically solve that
  endogenous-child problem").

### 3.5 Inventory: what must be re-proved at growing degree

Status column: "fixed d" means proved internally at fixed d only. "LV" means Lee–Venkatesan prove an analogue at
growing degree (d ≍ (log N/log log N)^{1/3}). LV's coverage was checked by the literature critic (§5).

| ingredient | needed at growing (d, k) | status | LV at growing d? |
|---|---|---|---|
| recursion / shape test | Lemma F at Λ_α | **proved** (§3.1) | n/a (LV use HECM for L[2/3] norms) |
| split counting | Theorem K(ii) with k ≍ D^α, no 2^ω | **proved** (Lemma S″) | no |
| K(iii) waste ratio | charge T_SQ per attempt | **proved** (bookkeeping) | no |
| FM: special-q first moment | precision (log N)^{−O(1)}, exponent uniform in d, k | **open**: proved at fixed d with polylog loss, but exponents and N₀ depend on (d, k) | partly: LV prove a first moment E_{m,f}\|X\| ≥ L[1/3, τ] for the **unrestricted** box at L^{o(1)} precision (LV Thm 5.4, Lemma 5.21 proved in §8; equidistribution steps Lemmas 5.13/5.18/5.20) |
| Lemma 8.5′ (restricted BV, moduli ℓ·s) | Y = L[1 − α] (L[3/5] at 2/5) | **open**: proved only in regimes (A) Y = exp((log N)^{3/4}) and (B) Y = N^{1/(2k)}, k fixed | LV Lemma 8.5 is the model (moduli without the ℓ structure) |
| LD*: second moment | the same precision; Kern's C_d, the "≤ 12 partners", L2_fam's (log N)^{O(C₀)}, all uniform | **open** | no: LV have no second moment |
| Lemma W (Brun–Titchmarsh for w) and the slack η₀ | θ → 0 is incompatible with fixed η₀ > 0, so η₀ = η₀(N) must shrink, and the constant 24/η₀ and every η₀-margin must be re-checked | open (bookkeeping) | — |
| Prop PZ ⇒ D_k | γ ≥ (log N)^{−O(1)} | follows if FM and LD* do | no |
| field rejections (Lemma IR, f′(m), content) | Pr(𝓡) = o(γ) | likely routine; IR's density 1/d must be re-checked at growing d | yes: LV Lemma 5.1 (Turán sieve), Pr[f reducible] ≤ L[2/3]^{−c} |
| LV6 (Lemmas C, H, E, X, SQ) | log\|H\| = O(d log N); character primes; Stark's exceptional zero | **gap as written** (codex recurse_2): Lemma X's zero-free threshold ν₁ ≍ 1/((2d)! log Δ̄) (lv6 §3.3) defeats a polynomial bit budget at d → ∞. Plausible repair: character primes of (2d)!·(log N)^{O(1)} bits, which costs exp(O(d log d)) = L[α]^{o(1)} since δ < α; not written | **yes**: LV §6 is at growing degree, unconditional, via Stark's effective Brauer–Siegel (LV Fact 6.22/Cor 6.23) |
| batch smoothness (Prop BS) | as is | proved | — |
| Theorem DD and (T) | as is | proved | — |
| multi-level (α < 2/5) | nested descendability (endogenous children) | **no tool** | no |

**Conclusion.**
- Recursion plus one-level descent would give a rigorous L[2/5, c] *if* the four "open" rows were proved with
  polylogarithmic precision uniform in (d, k), together with the LV6 and η₀ repairs below. These are sufficient
  inputs, not necessary conditions; the first draft's "iff" was unsupported (codex recurse_2). We expect C_d growth to make 3/7 the realistic first target.
- The decisive missing input is a quantitative, explicit-error version of LV §5, plus Kern with dimension-free
  constants, at special-q lattices. LV's L^{o(1)}-precision first moment is not enough, because γ is raised to the
  power 2k with k ≍ D^α.
- One possible route around the precision requirement is UNVERIFIED. It uses several fields that share the rational
  side G = g₁x − g₀y (same g, independent (w_j, c_j), Coppersmith-style). Then only the rational-side statement (PC_rat,
  Theorem R) needs polylog precision, and the algebraic side needs only L^{−o(1)}. This needs conditional-on-g moments
  that nobody has proved.

## 4. Prior art (Opus literature subagent, primary sources where accessible; 2026-10-20)

- **Lenstra–Pomerance, JAMS 5 (1992)** [read]. Algorithm 10.1 factors the form coefficient a ≈ √n with the ECM-based
  Algorithm 7.2, at L^{o(1)} per number (Thm 7.3). The only recursion is on factors already found: Algorithm 7.2, and
  Algorithm 10.4 / Thm 10.5 on at most log n divisors. The authors note the time is exactly L[1/2, 1 + o(1)] and that
  "one cannot improve this by choosing the parameter y differently" (p. 509).
- **Pomerance 1987** (Discrete Algorithms and Complexity) [read].
  - Thm 2.1: an ECM smoothness procedure that recurses on found factors.
  - pp. 138–139: inside the discrete-log algorithm, one option is "to apply algorithm REC … to p − 1, completely
    factoring it". This is the closest precedent for using one's own rigorous factoring algorithm as a subroutine, but
    it is a single call on a different number, not a recursion that feeds a fixed point.
- Seysen 1987: not accessible. LP cite L[1/2, √(5/4)] under GRH. Lenstra 1987 (2.12) says Seysen's smoothness test is
  Pollard ρ / Pollard–Strassen.
- Lenstra–Pila–Pomerance I/II: HECM, with no recursion. Part III, which LPP cite for details, was apparently never
  published.
- Schnorr–Lenstra 1984: "recursion" there is a group walk. Dixon 1981 and Pomerance 1982: early-abort induction on the
  number of aborts, not on N.
- Vallée 1991: abstract only.
- Lee–Venkatesan 2018 [read]: HECM (Fact 3.18, Cor. 3.19) for L[2/3] norms. No recursion.
- Discrete-log analogues:
  - Enge–Gaudry 2002 [read]: a one-off t_f factorisation of the group order.
  - Kleinjung–Wesolowski: recursion on polynomial degree, not integer factoring.
  - Hafner–McCurley, Diem, Buchmann: not accessed.
- Harsha–Mitropolsky–Rosen (arXiv 2209.10509) §5 calls "an algorithm for factoring that uses oracle calls on smaller
  numbers" an open problem in the *polynomial-reduction* (TFNP) sense. An L-overhead recursion such as ours is not
  covered by that remark.
- **Summary.** No rigorous analysis that we found recurses the main algorithm on √N-size auxiliaries, and no
  self-referential L-constant equation. The reason is clear: when the auxiliary needs only *smoothness*, an L^{o(1)}
  test exists, and full recursive factorisation would cost more. Recursion pays only when the test needs complete
  factorisation into N^{Θ(1)}-size primes, which is Prop H's horn 2. Not accessed: Seysen, Vallée, Dixon and
  Hafner–McCurley full texts; Bach–Shallit; Crandall–Pomerance §6.

## 5. Reviews

All findings were adjudicated, and every valid fix is applied in the text above.

**GPT-6.1-sol blind round `codex/recurse_1` (prompt `recurse_1.md`; the full report is copied to `notes/recurse_codex1_report.md` from
`codex/scratch/recurse_1/recursive_factoring_report.md`, git-ignored).** The first launch hit a "model at capacity"
error with no output; the relaunch completed. Its findings:
- Independently proved the recursive reduction c = √((1 + √2)θ), "including stopping times, recursive inputs,
  restarts, and transcript invariance".
- Found the one-level 2/5 floor.
- Proved a disjoint-products lemma equivalent to S″. Its hypothesis is n ≥ 2·2^ω k² for the full k-subset family; mine
  is |G| ≥ 4k³2^ω, because of the window restriction.
- **Found my error on multi-level descent.** With O(log log N) levels of geometric *constant* contraction, the size
  model attains α = 1/3 itself. My exponent-only count had said "1/3 not attained". This is accepted and corrected in
  §3.4 and the verdict.
- Weakest point, which I agree with: "a quantitative relation-supply theorem for the entire adaptive descent tree in
  one growing-degree field".

**GPT-6.1-sol adversarial round `codex/recurse_2`.**
- Verdicts: RA CORRECT WITH FIXES; P CORRECT WITH FIXES; S″ CORRECT; growing-degree floors CORRECT, with a GAP in the
  inventory's sufficiency claim.
- Fixes applied:
  1. The S² ≥ ρ_r^{−2}L^c margin is a power of L, not of N.
  2. (A2) is simpler by charging 2k searches per attempt.
  3. Lemma OS concerns the *logical* transcript, without timing.
  4. In T_rec, "independent" becomes conditional fibre invariance.
  5. Corollary P needs a real scheduler, or restarts; finite mean does not give finite E[T log T].
  6. **LV6 at growing d has a definite gap.** Stark's factorial threshold ν₁ ≍ 1/((2d)! log Δ̄) defeats a polynomial
     character-bit budget. A plausible repair costs exp(O(d log d)) = L[α]^{o(1)}, but it is not written.
  7. Fixed η₀ is incompatible with θ → 0.
  8. "iff" becomes "if" in the inventory conclusion.
- Its independent table reproduces the c\* values. Its weakest point is inherited: Lemma 8.5′'s absolute-character
  majorant (thmA audit item 7).

**Opus proof critic (read-only, independent).**
- Verdicts: RA ACCEPT WITH FIXES; P ACCEPT WITH FIXES; S″ ACCEPT. The size model is CORRECT WITH FIXES.
- It checked every place a enters the sources: "nothing it uses from Theorem A's proof is secretly tied to a = θ/√2".
  The a-dependent constants (C₁ in FM 2(c), C₀ in L2_fam and DS) are free fixed constants that change only polylog
  exponents and N₀.
- Fixes applied:
  1. Canonical oracle output in Lemma OS.
  2. The domain of L_x in Lemma F.
  3. The restart schedule for Corollary P. This was its weakest point: live dovetailing has a model-dependent
     overhead.
  4. The K(iii) waste ratio keeps one 2^ω after S″, not 2^{2ω}.
- Rejected, with reason: it read log u ≈ (1 − α) log log N in (R2). That holds for fixed θ. At the relevant scale
  θ ≍ D^{2α−1}, u ≍ D^α and log u ≈ α log log N. The conclusion is the same either way, and the text now says both.

**Opus literature critic (read-only, independent).**
- Every attribution it could check against a primary source is CONFIRMED: LP92 Algorithms 7.2/10.1/10.4, Thm 7.3,
  the p. 509 remark; LV degree, Thm 2.1, the first moment and stochastic deepening (no second moment), Fact 3.18 /
  Cor 3.19, Lemma 5.1, §6 via Stark, Lemma 8.5; Harsha–Mitropolsky–Rosen §5.
- No recursive fixed-point precedent found.
- The current rigorous record is still LP92 (Harvey 2020 §1; Pomerance 2008; Galbraith ch. 15). The project's claim
  would beat it, and it rests on unrefereed internals.
- Additions applied: LP Thm 6.1 (ψ(x, y; 𝒫)), de Boer–Pellet-Mary–Wesolowski Thms 1.1/9.5/17.7 as the nearest tool
  for good sets, LV Thm 5.4 / Lemma 5.21, Enge–Gaudry–Thomé and JLNT 2009 for the descent recursions.
- Pomerance 1987 was UNVERIFIABLE for this critic. My prior-art subagent read an OCR scan of it, so the pp. 138–139
  statement rests on one reading.
- lemmaM's attribution of a "one-stage L[2/5]" to Gordon 1993 is UNVERIFIABLE: Guillevic quotes Gordon's individual
  logs at L[1/3, 9^{1/3}]. That attribution is lemmaM's, not ours; it is noted here and not used.

**GPT-6.1-sol r19_transfer item 1/4 (orchestrator's transfer round, read after my own derivations of §§1, 3.3).** It
proposed the a-freedom check, the fixed point, the portfolio corollary, and the 2/5 and 3/7 targets. Corollary P is
credited to it.

## 6. Candidate card

- **Mechanism.** Recursive complete factorisation of the Dixon-split halves inside the algorithm's own shape test.
- **Model.** Theorem A's design (a‴_k; μ), with B = L^a free.
- **Result.** Proved from Theorem A: L[1/2, √((1 + √2)θ)] at fixed d, 0.851 at d = 9, → 0 as d → ∞, and a portfolio
  exp(o(Λ)).
- **Failure scope.**
  - Below L[1/2], the descent geometry gives 2/5 with uniform polylog precision and 3/7 with (log N)^{O(d)} losses.
  - Multi-level descent reaches 1/3 in the size model, using O(log log N) contracting levels, but it is blocked
    rigorously by endogenous children.
- **Next step.** An explicit-error FM/LD* (Kern, Lemma 8.5′ at Y = L[3/5]) uniform in (d, k), aiming at L[3/7]
  first.
