# coverage38: can worst-case coverage beat c = 3/8? New designs, reciprocity, dichotomy, hypotheses

Agent `coverage38`, 2026-10-18 (phase 6).
- Pre-registration: `notes/coverage38_prereg.md` (commit 757af72, before any registered run).
- Code: `src/coverage38/{theory,lehman,oracle,experiments,analyze,explore}.py`. Tests:
  `tests/test_coverage38.py`.
- Raw data: `experiments/raw/coverage38_{e2,e2s,e3}.jsonl` (pre-registered) and
  `coverage38_explore.jsonl` (exploratory). Verdicts:
  `experiments/analysis/coverage38_summary.json`.
- GPT-6.1-sol: blind round `codex/coverage38_1.out`; adversarial review `codex/coverage38_2.out`.
  Opus critic: §9.
- Notation is that of notes/dchworst_results.md (Theorems W, R; barriers B1–B4).

## 0. Verdict (plain language)

**No worst-case exponent below 1/5 is proved, and no log-factor gain either. What is new is a
cleaner and more design-independent form of the barrier. It is a statement about this
algorithm family under explicit modelling assumptions, not a lower bound for factoring.**

- **Theorem U (proved in the cost model).**
  - *Assumptions.* An N-independent list of (weight, modulus) groups, each charged ≥ 1
    operation per input. The window hypothesis WIN, which is verified for full windows and
    for stationary-point (Lehman) windows.
  - *Statement.* For such a design the worst-case cost satisfies

      X ≥ c · N^{1/5} · (λ⁴/J)^{1/5}.

    Here J is the number of distinct moduli and λ is the expected number of covering events
    for a random pair.
  - *Who found it.* GPT-6.1-sol derived the same inequality blind. GPT-6.1-sol's referee
    points to close prior art: Ward–Maestro's Sept-2026 report uses the same
    Cauchy–Schwarz/reciprocal-height architecture for its 1/5 and 1/6 floors. Our addition
    is the multi-modulus Hölder step and the (λ, J) form.
  - *Consequence.* One dyadic band attains the bound up to constants. Every design in the
    family therefore faces the same (λ, J) threshold. A sub-1/5 exponent through Theorem R
    needs coverage of every pair by a design with λ < J^{1/4}. That is Theorem W's c > 3/8,
    and it applies to wide ranges, multiscale bands, composite moduli, prime powers, signed
    weights and stationary windows alike.
  - *Exponent.* Coverage at λ = J^θ gives exponent (1+θ)/(6+θ): 1/6 at θ = 0, 1/5 at
    θ = 1/4.
  - *Log gains.* A log gain also needs λ ≤ J^{1/4}(log N)^{−a}, which buys (log N)^{4a/5}.
  - *What it does not say.* Theorem U does not rank how hard the coverage statements of
    different designs are to prove.
- **Proposition L (proved; an elementary concentration lemma).**
  - *Statement.* Put J points uniformly at random in the complement of the target sets, so
    that they have zero hits. Then with probability ≥ 1/4 all Q prescribed bounded linear
    statistics stay within 2λ + 2√(J log 8Q) of their uniform-model values.
  - *Consequence.* Generic simultaneous deviation bounds at that level cannot certify
    coverage. Nor can the standard phase-discarding Fourier argument. Erdős–Turán with a
    common bound κ needs λ ≳ κ‖Ŝ‖₁/q ≥ √J.
  - *Applied to known bounds.* Every known or conjectured *pointwise* Kloosterman bound is at
    or above that level: Hooley's R*, Karatsuba, Bourgain–Garaev, Korolev. Theorem U needs
    λ < J^{1/4}, so none of them gives a gain through these arguments.
  - *Scope.* This is not a theorem that exponential sums can never help. The real
    configurations x_ℓ = p·ℓ̄ mod q are arithmetic, and a proof can use their structure.
- **Literature status.**
  - For prime-variable Kloosterman sums mod a prime, nothing nontrivial is known at our length
    q^{0.32}. Bourgain, Baker and Korolev all start at q^{1/2+ε}.
  - For integer variables only log savings are known: Y·log log q/log q (Korolev 2016).
- **What would suffice (conjectured, not standard).** A uniform positivity count of band
  primes ℓ with p·ℓ̄ mod q in a Farey-arc set of measure J^{−3/4}.
  - Equivalently, one of H dilates Bp of every residue class mod q is a product
    (band prime) × (integer ≤ Hq/Y), with product size ≈ q^{1.04}.
  - The known uniform sufficient range for small products in residue classes is q^{3/2}.
  - I found no standard named hypothesis with a verified implication. The checks covered
    GRH, Elliott–Halberstam, R*, Montgomery, Cramér-type, Erdős–Odlyzko–Sárközy and
    equidistribution of t_ℓ.
  - This is a literature finding, not a non-implication theorem.
- **Dichotomy (proved in part; no gain in the data).**
  - *Linear families are easy.* Every fixed linear family (p, q) = (ax+b, cx+d), including
    (x, x+1), has a constant δ = Bp − Aq at one weight, so Lehman/Fermat factors such N
    quickly. **Barrier B3 therefore does not constrain "covered or Lehman-easy" proofs.**
  - *B3's survivor.* B3 survives in a sparse non-interval form: (x, ⌊αx⌋) with α badly
    approximable (GPT-6.1-sol).
  - *Data.* The failure sets are not disjoint.
    - Pre-registered with the uniform-window Lehman score: removing Lehman-easy pairs changes
      the exhaustive worst case over S_32 and S_36 by 0 in all 14 cells. At P ≈ Y³
      (n = 60, 72) the uncovered rate of easy pairs is 0.97–1.01 times that of hard pairs.
    - Exploratory, with the corrected Dirichlet-optimal score L*: hard pairs are uncovered
      slightly *more* often at n = 72 (R = 0.93–0.94, CI excludes 1). The dichotomy removes
      the wrong pairs.
- **Public intervals.**
  - Coverage of the true pair is a statement about the hidden integers δ_v.
  - The natural N-public certificate is full-candidate covering: every x in [√(N/2), √N] has
    some Bx² − AN with a band divisor. That is sufficient, not necessary.
  - It is a covering by N-determined quadratic classes, not the "classes 0 mod ℓ" statements
    that Jutila / Friedlander–Lagarias / Harman / Ford handle (assessment).
- **e3 (pre-registered, all of O1–O5 pass).**
  - On real configurations at n = 48 (J = 1612), the pointwise statistics of uncovered and
    covered pairs agree to within 3.3%.
  - They also agree with the zero-hit random surrogates of Prop L to within 7%.
  - The absolute Erdős–Turán error exceeds λ by a factor 113–287.
- **Conditional results.** No conditional result from a standard number-theoretic hypothesis.
  - Under P = BPP, padding plus Lenstra–Pomerance already gives deterministic
    L[1/2, O(1)] factoring. That is unrelated to coverage.
  - Known conditional deterministic results are Schoof 1982 (N^{1/5+ε}, GRH) and Umans–Wang
    (1/6, SPDC, stated for infinitely many sizes).

## 1. Setting

- N = pq, P = √N, (p, q) ∈ S_n (p < q ≤ 2p primes).
- A weight is a primitive v = (A, B), B ≥ 1, A ≠ 0, of height h(v) = max(|A|, B). Write
  δ_v = Bp − Aq.
- **Designs.** A design is an N-independent list 𝒢 of *groups* g = (v, M) (distinct pairs).
  - On input N, group g computes the roots of Br² ≡ AN (mod M), which may be none. Each root
    gives a target with the class lock mod M² (HH 2022 Prop 3.3) and a window W_τ(N) of
    g-values, whose placement may depend on N.
  - This covers D′, Theorem W's designs, Lehman designs, multiscale and composite designs:
    every design in the project.
  - Weights or moduli chosen adaptively from N are *outside* this model, as in below16's
    locality caveats.
  - **Charging convention.** Every group in 𝒢 is charged ≥ 1 operation on every input, its
    root test, even when it has no roots. below16's 𝔐⁺ charges only issued labels, so this
    is an explicit extra convention. Theorem R's construction-cost term |V||𝒫_k| does
    charge the rejected combinations. The convention describes algorithms that enumerate
    their groups. It is not a lower bound for every batching or early-exit strategy.
- **Coverage weight.** π_g is the probability, over uniform (p, q) ∈ S_n, that some target
  of g covers the true pair.
- **Notation.** λ = Σ_g π_g. J is the number of distinct moduli in 𝒢. T = |𝒢|.
  w̄_g = E_N Σ_{τ∈g} |W_τ(N)| is the mean window mass of g over the population.
- **Cost model (𝔐⁺, below16 §1; dchworst Theorem R (iii)).**
  - Every group costs ≥ 1 operation on every input, namely its root test. Theorem R charges
    the rejected (weight, modulus) combinations explicitly (its |V||𝒫_k| term).
  - The window mass costs (Σ_τ|W_τ(N)|)^{1/2} on input N.
  - So the worst-case cost satisfies X := max_N cost(N) ≥ max(T, (max_N Σ|W(N)|)^{1/2})
    ≥ max(T, (Σ_g w̄_g)^{1/2}).
  - Theorem R realises Õ(T + (Σ|W|)^{1/2}) + N^{1/6+o(1)} by choosing the baby length m, so
    the model is tight for this family.
- **Window hypothesis (WIN).** w̄_g ≥ c₁ · π_g² · h(v) · P.
  - *Full windows* (D′, Theorem W).
    - For prime M, the group is active (has roots) for about half of all N, with 2 roots.
      Each root window has ≍ (S-range)/M² g-values.
    - The S-range of Bx + AN/x over [√(N/2), √N] is ≥ (1 − 2^{−1/4})²·h(v)P ≈ 0.0253·h(v)P
      for every positive weight, with the minimum at A/B = 1/√2. Signed weights do better
      (both referees).
    - Hence w̄_g ≍ h(v)P/M². For composite M the activity is ≈ 2^{−ω(M)} with 2^{ω(M)}
      roots, which gives the same estimate.
    - π_g ≍ 1/M needs Brun–Titchmarsh over the prime population, valid for M ≤ P^{1−ε}.
      For M near P only π_g ≲ log P/M holds. These log factors matter for U4.
  - *Lehman / stationary-point windows*: only |δ_v| ≤ Δ is searched.
    - S_v − 2√(ABN) = δ_v²/(S_v + 2√(ABN)), so the window has ≍ Δ²/(h(v)PM²) g-values.
    - The coverage probability is π ≍ Δ/(M h(v) P). This needs an upper-bound sieve and loses
      up to a log log factor when Δ/M = O(1).
    - Both relations give WIN up to these constants. This needs A, B > 0 comparable, the
      stationary point in range, and Δ ≥ M (a window with Δ < M sees no multiple of M).
  - *Arbitrary N-dependent placements are NOT covered.*
    - A "magic" singleton window containing the true g, e.g. at (A, B, M) = (1, 1, 3) with
      the window {(p+q−c)/9}, has mass O(1) and positive coverage probability. It violates
      WIN (GPT-6.1-sol, Opus critic).
    - So WIN for all efficiently computable placements is a computational-hardness
      assumption, not a lemma.
    - **Theorem U is a model theorem for full windows and stationary-centred windows**, and
      for any other window rule for which WIN is proved.

## 2. Theorem U (universal cost bound) — proved in the model

**Theorem U.** Under 𝔐⁺ and WIN, every N-uniform design satisfies

  X ≥ (c₁P/(2√3))^{2/5} · λ^{4/5} · J^{−1/5},  and also  X ≥ (c₁P/(2√3))^{1/3} · λ^{2/3}.

*Proof.*
- Groups are indexed by distinct pairs (v, M). All windows of one weight at one modulus form
  one group, and WIN is assumed for their union, which is a single window set.
- Fix a modulus M, carrying t_M groups (distinct weights), and put u_M = Σ π_g over them.
- Cauchy–Schwarz gives

    Σ_v π_v² h(v) ≥ u_M² / Σ_v h(v)^{−1}.

- **Height count.** At most 3h² primitive weights have height ≤ h. So the k-th smallest
  height is ≥ √(k/3), and Σ_v 1/h(v) ≤ √3 Σ_{k≤t_M} k^{−1/2} ≤ 2√3 √t_M.
  (`test_kth_height_bound` checks this.)
- Hence Σ_g w̄_g ≥ (c₁P/(2√3)) Σ_M u_M² t_M^{−1/2}.
- Hölder with exponents (2, 4, 4):

    λ = Σ_M (u_M² t_M^{−1/2})^{1/2} · t_M^{1/4} · 1 ≤ (Σ_M u_M² t_M^{−1/2})^{1/2} T^{1/4} J^{1/4}.

- So Σ_g w̄_g ≥ (c₁P/(2√3)) λ²/(TJ)^{1/2}.
- With X ≥ T and X ≥ (Σ_g w̄_g)^{1/2} this gives X^{5/4} ≥ (c₁P/(2√3))^{1/2} λ J^{−1/4}.
- The second bound follows from J ≤ T ≤ X. ∎

The moduli sizes cancel out entirely. Nothing about bands, primality or scale enters.
Caveat: `test_theorem_U_lower_bound_random_designs` imposes WIN with equality, so it checks
only the algebra.

**Attainment (proved by direct optimisation; `test_single_band_optima_match_universal_formula`).**
- **One band, full windows.** T = λY, Σ|W| = λHP/Y and H² = λY/J. The optimum is at
  Y = (P²/(λJ))^{1/5}, where X = P^{2/5}λ^{4/5}J^{−1/5}. This is feasible while
  J ≤ (P²/λ)^{1/6}.
- **One band, Lehman windows.** T = JH² and Σ|W| = λ²P/(JH), with Y free in [J, JH²/λ]. The
  optimum is H = (λ√P)^{2/5}J^{−3/5}, with the same X.
- So **one dyadic band attains Theorem U up to constants.**
- The bound is not exactly sharp: Cauchy–Schwarz equality would need π_v ∝ 1/h(v), which
  fails for full windows.
- Uniqueness is not claimed. For stationary windows, (M, Δ) ↦ (tM, tΔ) leaves the modelled π
  and window mass unchanged, so several scales can tie (GPT-6.1-sol).

**Corollaries (proved in the model).**
- **U1 (the required statement).** A design that covers every pair at cost X = N^{1/5−η} has
  λ ≤ C·J^{1/4}·N^{−5η/4}. So any worst-case proof through Theorem R must establish coverage
  for a design whose expected count is below the fourth root of its number of moduli. For
  J ≍ Y that is Theorem W's c > 3/8 (H ≤ Y^{1/8−δ}, P ≍ Y^{25/8}).
  - Every design faces the same (λ, J) threshold.
  - Theorem U does not compare how hard the different designs' coverage statements are to
    prove.
- **U2 (θ-exponent).** Suppose coverage holds for a design with λ = J^θ. (A stronger
  coverage theorem or another design is not excluded.)
  - A single band then gives N^{(1+θ)/(6+θ)}: 1/6, 9/49, 1/5 and 3/13 at θ = 0, 1/8, 1/4 and
    1/2.
  - For θ ≥ 1/4, Theorem U gives X ≥ cN^{1/5}J^{(4θ−1)/5} ≥ cN^{1/5} for every design. The
    unconditional fallback (J = 1, Minkowski) is exactly 1/5.
- **U3 (the 1/6 floor).** For λ ≍ 1, X ≥ cN^{1/6}.
  - This is the average-case floor of multimod Thm A and below16 T1, recovered in this model
    by a two-line argument.
  - It also means no design covering a fraction σ of pairs beats N^{1/6}σ^{2/3}, in this
    *worst-case* budget. It is not an average-running-time statement, since E√W cannot be
    replaced by √(EW).
- **U4 (log factors).** λ = J^{1/4}(log N)^{−a} would buy N^{1/5}(log N)^{−4a/5} in the list
  cost of an attaining design. That is a model scale, not an algorithm.
  - Theorem R's own polylog overhead must be tracked against Harvey–Hittmeir's
    N^{1/5}log^{16/5}N/(log log N)^{3/5}.
  - No such coverage statement is proved. A log gain needs the same J^{1/4} beating of the
    trivial union bound.

**Design families (orchestrator item 1).** Each row is a corollary of Theorem U (under WIN).
The Ford remark is heuristic.

| design | Theorem W analogue | verdict |
|---|---|---|
| one dyadic prime band, full windows | λ⁴ < J (c > 3/8) | extremal |
| Lehman / stationary-point windows, any Δ (span ∝ m² near the stationary point) | λ_Δ⁴ < J, with λ_Δ counting only the \|δ_v\| ≤ Δ events | same threshold. Using the full-window λ for Lehman designs is wrong (GPT-6.1-sol) |
| wide range (Y, Y^{1+κ}] | λ ≍ H² log(1+κ); bound unchanged | no better than one band up to constants |
| multiscale bands with per-scale weights | additive λ, J | no better up to constants |
| composite / smooth moduli | same inequality (root multiplicity and factorisation of M charged) | no better in λ. *Heuristic:* likely worse in Pr[covered]. A single integer has a divisor in (Y, 2Y] with probability ≍ (log Y)^{−0.086}(log log Y)^{−3/2} (Ford 2008) while its expected number of such divisors is log 2, so composite events cluster. Ford's theorem concerns single integers, not the joint law of the H² values δ_v |
| prime powers ℓ^k | J ≤ #primes near Y^{1/k} | worse (fewer moduli for a given scale) |
| signed weights, (1,1), both factors, rectangles | constants in λ; a square box minimises the window for a given count | same threshold |

**A note on the class lock for composite moduli (proved; also GPT-6.1-sol).**
- D′ Lemma 1 uses only invertibility. Write p = r + Mx with Br² ≡ AN (mod M) and
  gcd(M, ABN) = 1. Then Bp + Aq ≡ Br + ANr^{−1} (mod M²).
- Theorem R's Lemmas O, G and X never use primality of the moduli; they use u_M^{M²} = γ.
- So Theorem R extends to squarefree moduli whose prime factors are ≡ 3 (mod 4) and > H.
  - The roots are deterministic by CRT, with one target per CRT root.
  - The extension needs M < √(N/2), or an explicit gcd(M, N) test, and the factorisation of
    M with its root multiplicity 2^{ω(M)} charged.
  - The prime-only accounting "exactly one of ±AN/B is a square" fails for composite M. For
    example N = 989, M = 21: neither sign has roots (GPT-6.1-sol). The per-group charge
    covers this.
- **Consequence (narrowed after review).**
  - *Prime moduli are unnecessary for extraction.* That does not make a prime-detecting
    (Harman) sieve useless for a coverage proof: admissible composites still need their own
    distribution theory.
  - Harman's sieve can exploit Type I/II cancellation before absolute values. It still needs
    integer-moduli counts of the target-set hits at precision o(λ), so the needed precision
    remains the obstacle (assessment).

## 3. Reciprocity and the linear-statistics barrier

### 3.1 Reciprocal form (proved)

- Suppose ℓ | Bp − Aq, say Bp − Aq = ℓk. Then k ≡ Bp ℓ̄ (mod q), and the exact condition is
  ‖Bpℓ̄/q − Bp/(qℓ)‖ ≤ H/ℓ. The centre moves with (B, ℓ); this correction is GPT-6.1-sol's.
- **Sufficient zero-centred test.** Suppose B ≤ H/4, ℓ ∈ (Y, 2Y] and
  |Bpℓ̄ mods q| ≤ qH/(4Y). Then the corresponding A has |A| ≤ 3H/4, so the pair is covered at
  height H.
- So for a fixed pair, coverage follows when some band point x_ℓ = p ℓ̄ mod q lies in the
  **Farey-arc set**

    S = {x : ∃ B ≤ H′, |Bx mods q| ≤ k₀},

  which has measure μ ≍ H′k₀/q and expected hits λ = Jμ.
- `test_reciprocal_hit_is_coverage` checks the implication exactly.
- Equivalently, some dilate Bp (B ≤ H′) is congruent mod q to ℓ·k with |k| ≤ k₀. That is a
  product of a band prime and an integer of size ≤ Hq/Y. At the threshold the product size
  is ≈ q^{1+1/25}.

### 3.2 Proposition L (proved)

**Setting.**
- Let q ≥ 2 and J² ≤ q/4.
- For j = 1..J let S_j ⊂ Z/q have |S_j| ≤ q/2, and put λ = Σ_j |S_j|/q.
- Let 𝓕 be a family of Q statistics F(x) = Σ_j f_{F,j}(x_j) with |f_{F,j}| ≤ 1. These are
  bounded linear statistics. Examples:
  - exponential sums Σ_j a_j e(h x_j/q) with |a_j| ≤ 1;
  - character sums;
  - sums over sub-progressions of moduli.

**Claim.** There are pairwise distinct x_j ∉ S_j such that, for every F ∈ 𝓕,

  |F(x) − E_unif F| ≤ 2λ + 2√(J log 8Q).

*Proof.*
- Draw x_j independently and uniformly from C_j = Z/q ∖ S_j.
- **Mean shift.** E_unif f = (1−μ_j)E_{C_j} f + μ_j E_{S_j} f. So
  |E_{C_j} f − E_unif f| = μ_j |E_{C_j} f − E_{S_j} f| ≤ 2μ_j, and summed over j the shift
  is ≤ 2λ.
- **Deviation.** Apply Hoeffding to the real and imaginary parts, whose terms lie in
  [−1, 1]: P(|F − E F| ≥ s) ≤ 4e^{−s²/(4J)}. A union bound over 𝓕 at s = 2√(J log 8Q) has
  failure probability ≤ 1/2.
- **Distinctness.** A collision has probability ≤ J²/(q/2)·(1/2) ≤ 1/4.
- So all three hold with probability ≥ 1/4. ∎

**Complements (proved).**
- **(a) Parseval.** For J distinct points, Σ_{h≠0}|K(h)|² = J(q − J). So
  max_{h≠0}|K(h)| ≥ √(J(q−J)/(q−1)), and no uniform pointwise bound is below ≈ √J.
- **(b) Erdős–Turán sufficiency.** Count = λ + (1/q)Σ_{h≠0} Ŝ(−h)K(h). So pointwise bounds
  |K(h)| ≤ κ_h certify coverage once λ > (1/q)Σ_{h≠0}|Ŝ(h)|κ_h.
  - ‖Ŝ‖₁/q is *not* ≍ log H′ (corrected). At q = 200003 and k₀ = 100 the measured values are
    3.1, 6.0, 16.1, 27.9 for H′ = 1, 4, 16, 32 (Opus critic). In e3 the value is ≈ 9.5 for
    H′ = 4 and ≈ 23.9 for H′ = 16. In general 1 ≲ ‖Ŝ‖₁/q ≲ H′ log k₀.
  - This makes the Erdős–Turán requirement λ ≳ √J·‖Ŝ‖₁/q *worse* than √J.
  - Applied per B separately the budget is ≍ Hκ (GPT-6.1-sol).
- **(c) What Prop L excludes, precisely.** It excludes deductions whose only premises are
  estimates of the form |F(x) − E_unif F| ≤ κ_F for F ∈ 𝓕, with every
  κ_F ≥ 2λ₀ + 2√(J log 8Q). From such premises no coverage at λ ≤ λ₀ follows, because the
  zero-hit configuration satisfies them.
- **Not excluded:**
  - uniform bounds κ ∈ [√J, 2√(J log 8Q)) used through other routes;
  - frequency-dependent bounds below the random level on a sparse set of frequencies.
    Parseval forbids only *uniform* sub-√J bounds;
  - estimates centred at a different main term;
  - statistics that see the arithmetic structure of x_ℓ = p·ℓ̄. For example
    Σ_j e_q(ℓ_j x_j) = J·e_q(p) on every real configuration but not on generic ones
    (GPT-6.1-sol).
- **Bernstein remark (heuristic resolution calculation, not proved here).** An indicator
  statistic of sets of measure ν fluctuates by ≈ √(Jν log 8Q) + log 8Q. Resolving the shift
  λ needs ν ≲ λ²/(J log 8Q).
  - At λ = J^{1/4} this is ν ≲ J^{−1/2}/log, which is still a power J^{1/4} larger than the
    targets' μ = J^{−3/4}.
  - So random-level estimates on such sets do not suffice. Coarser equidistribution is
    excluded only *via random-level estimates*.

**Consequence (restricted form, after both reviews).**
- Generic simultaneous deviation bounds at the level of Prop L cannot, by themselves,
  certify coverage among unrestricted configurations.
- A Fourier argument that discards phases and uses a common bound κ ≥ √J (Parseval: no
  smaller uniform κ exists) needs λ ≳ κ‖Ŝ‖₁/q ≥ √J.
- Theorem U needs λ < J^{1/4}. So **these arguments give no worst-case gain in any design of
  the family, whatever the quality of the pointwise bound up to square-root cancellation.**
  This makes B1 precise and design-independent.
- It does not exclude arithmetic proofs that use the structure of x_ℓ = p·ℓ̄ or phase
  information.

### 3.3 Literature status of the sums (orchestrator item 2a)

Our sums are K(a) = Σ_{ℓ∈(Y,2Y]} e(aℓ̄/q), q prime ≈ P, at the threshold length Y = q^{8/25}.

| sum | best bound found | applies at q^{0.32}? |
|---|---|---|
| ℓ prime, q prime | Bourgain 2005 (threshold q^{1/2+ε}). Baker 2012 (Acta Arith. 156, Thm 1): Y^{1−δ⁴/2000} for q^{1/2+δ} ≤ Y ≤ q^{3/4+δ}, 0 < δ ≤ 1/24. Garaev; Fouvry–Shparlinski (Acta Arith. 150, 2011, general moduli, q^{3/4} ≤ Y ≤ q^{4/3}): (Y^{15/16} + q^{1/4}Y^{2/3})q^ε, nontrivial for Y > q^{3/4+ε}. Korolev 2019 (Math. Notes 106, 89–97, Thm 1): Y q^{−ε²(1−7ε)} for q^{1/2+ε} ≤ Y ≤ q^{3/4}, 0 < ε ≤ 1/53 | **no**: nothing nontrivial below q^{1/2} |
| ℓ prime, averaged over q | Irving (arXiv 1301.6372): Thm 1.1 (max over a) from Q^{2/3+δ}; Thm 1.3 (fixed a, size-dependent prefactor) from Q^{1/2+δ} | no (one fixed q here) |
| n integer, q prime | Karatsuba's method. Bourgain–Garaev, *Sumsets of reciprocals …*, Izv. Math. 78 (2014), Thm 16: Y(log log q)³log q/(log Y)^{3/2}. Korolev 2016 (Math. Notes 100, 820–827) Cor. 1: Y·log log q/log q for q^{1/4+δ} ≤ Y ≤ q^{1/2−δ}, 0 < δ ≤ 1/8 | yes, but only a **log** saving |
| conjectural (Hooley R*: integer variable, intervals) | √Y·(a,q)^{1/2}q^ε | integer moduli only, no automatic transfer to primes. A pointwise bound: only λ ≳ √J through the arguments of §3.2 |

- Garaev–Shparlinski (Mathematika 69, 2023; arXiv 2304.07953) flag a gap in a different
  Korolev 2016 Sbornik lemma. That is not the Math. Notes corollary used above.
- The "small products in residue classes" literature is the right comparison for §3.1's
  product form.
  - xy ≡ a (mod p) with x ≤ X, y ≤ Y′ is solvable for XY′ ≫ p^{3/2}log²p (Weil;
    Browning–Haynes; Shparlinski's uniform count, error p^{1/2+o(1)}). That is a uniform
    *sufficient* range, not a necessary barrier.
  - Garaev noted that Ford's theorem makes almost all classes *not* products x·y with
    x, y ≤ p^{1/2}(log p)^{0.086−ε}.
  - Our need is product size Hq ≈ q^{1.04}, with a prime factor in a short interval and an
    H-fold union of dilates. This is beyond the known unconditional uniform ranges. R*, for
    integer moduli and via phase-discarding arguments, gives it only at H ≈ √Y for one
    dilate (Minkowski level), or λ ≳ √J for the union.

### 3.4 The exact estimate that would give c > 3/8 (orchestrator item 2c)

**Hypothesis K(θ) (conjectured; not a standard named hypothesis).**
- *Statement.* Fix θ > 0. For all large q: for every prime q ∈ (P/√2, √2P), every
  a ≢ 0 (mod q), and every Farey-arc set S of measure μ with J_Y μ ≥ J_Y^θ:

  #{ℓ ∈ (Y, 2Y] prime : a ℓ̄ mod q ∈ S} ≥ 1.

- K(θ) for some θ < 1/4 gives worst-case N^{(1+θ)/(6+θ)+o(1)} by Theorems R and U2.
- The Poisson-type version replaces J^θ by C·log q and gives N^{1/6+o(1)}.
  - The literal mean-one version is false. GPT-6.1-sol's example: q = 1 000 003, Y = 83,
    J = 15, H′ = 1, k₀ = 33 333 gives J|S|/q = 1.000002, yet p = 500 069 has no hits.
- GPT-6.1-sol's bilinear form: E = Σ_{h≠0} c_h Σ_{B≤B₀} Σ_ℓ e_q(hBpℓ̄) with
  |E| ≪ √λ (log q)^C. This is sufficient for positivity, not equivalent to it. Arbitrary-coefficient square-root bilinear bounds would be false.
  The coefficients must be the prescribed Farey-arc Fourier weights.
- By §3.2(c), K(θ) cannot be deduced from random-level equidistribution estimates. It is
  essentially DCH_pair in reciprocal coordinates.

## 4. Dichotomy designs (orchestrator item 3)

**Proposition D (proved).**
- **(a) Lehman side, operational.** Take weights 1 ≤ A ≤ B ≤ h₀ and scan
  S ∈ [⌈2√(ABN)⌉, ⌈2√(ABN)⌉ + X₀] for S² − 4ABN = D². This finds p whenever some such
  weight has δ_v² ≤ 4X₀√(ABN). The reason is S_v − 2√(ABN) = δ_v²/(S_v + 2√(ABN)).
  - The scan costs ≤ h₀²(X₀+1) square tests. It needs gcd(AB, N) = 1, which holds for
    h₀ < p.
  - The BSGS-shared version needs Theorem R's extraction/order machinery.
  - `lehman_factor` (sees only N) implements it. `test_lehman_score_implies_success` and
    `test_dirichlet_bound_all_pairs_small` (all of S_24) check it.
- **(b) Dirichlet, corrected after the Opus critic.**
  - *The first draft was false.* It claimed the uniform-window score
    L(p,q) = min_v max(h², δ_v²/(4√(ABN))) satisfies L ≤ C·P^{2/5}.
    - Counterexample: q/p ≈ 1 + 1/a with a ≈ P^{0.24}.
    - Measured L/P^{0.4} = 2.4, 5.6, 12.8, 27.5, 48.3 at n = 60…140, with log_P L → 0.48.
      `test_uniform_score_not_dirichlet_but_Lstar_is` checks this.
  - *The correct score* uses height-dependent stationary windows, whose BSGS cost per weight
    is ≍ |δ_v|/√(hP):

      L*(p,q) = min_v max(h², |δ_v|·√(h/P)).

  - *Bound.* Take the convergent with B_k ≤ T < B_{k+1}. Then |δ_k| ≤ q/B_{k+1} ≤ 2P/T, and
    T = P^{1/5} gives L* ≤ 2P^{2/5}. Measured: max L*/P^{2/5} = 0.79, 0.91 and 0.95 at
    n = 32, 60 and 72.
  - *Attribution.* The exponent 1/5 is the BSGS-accelerated stationary-window method
    (Harvey 2021; HH 2022). Lehman's own 1974 method is 1/3. Dirichlet-tight pairs sit
    exactly at the critical scale.
- **(c) Linear families are easy.** Take p = ax + b, q = cx + d, g = gcd(a, c), and the
  weight (A, B) = (a/g, c/g). Then δ_v = (cb − ad)/g is constant.
  - The scan in (a) costs O(h²(1 + δ²/(hP))) with h = max(a, c)/g. That is polylog for
    fixed (or polylog) coefficients.
  - In particular (x, x+1), with δ_{(1,1)} = −1, is Fermat-easy.
- **(d) What remains of B3.**
  - B3 obstructs only a *pure* coverage theorem valid for all coprime pairs. A theorem
    "covered **or** Lehman-easy" is not affected. The two have different quantifiers.
  - But B3 has a non-linear survivor (GPT-6.1-sol). For α ∈ (1,2) badly approximable, the
    pairs (x, ⌊αx⌋) satisfy |Bx − A⌊αx⌋| ≥ c_α x/h − O(h). So they are Dirichlet-tight
    throughout the Lehman range.
  - An all-coprime residual theorem therefore implies band-prime coverage of the sparse set
    {Bx − A⌊αx⌋ : |A|, B ≤ H} for every x with gcd(x, ⌊αx⌋) = 1. That is again ≍ H² values and not an interval.
  - The prime-pair-only version is not obstructed this way.
- **(e) Residual statement.** The needed theorem is "every Lehman-hard balanced pair is
  covered at λ < J^{1/4}". Theorem U is unchanged: dichotomy alters which pairs must be
  covered, not the λ⁴ < J requirement, provided the residual design keeps the accounting and
  WIN with λ defined over the original population. Prop L is an abstract statement and
  supplies no arithmetic obstruction specific to the residual population.
- **(f) Other partners.** I found no link from uncovered-ness to the structure of p − 1,
  q − 1 or orders. Uncoveredness is a condition on (p q̄ mod ℓ)_ℓ, i.e. on p in ≈ |V|
  N-determined classes per band prime. p − 1 BSGS, Fermat and Hittmeir see unrelated data.
  No provable disjointness was found (UNVERIFIED that none exists).

**Experiments: the failure sets are not disjoint (§7.1–7.3).**
- **Pre-registered, uniform-window score L.** h* is attained by Lehman-hard pairs in every
  cell. In three cells the recorded argmax pair itself is Lehman-easy (Opus critic). At
  P ≈ Y³ the uncovered rate does not depend on easiness.
- **Exploratory, corrected score L*.** At n = 72 hard pairs are ≈ 6% *more* often uncovered.
- **Heuristic explanation (not proved).** Weights with 0 < |δ_v| < Y can never be covered.
  Asymptotically there are O(1) of them even for Lehman-easy pairs, which matches dchworst
  e3's +5% "one dead weight" excess for near-rational families.

## 5. Public-interval reformulation (orchestrator item 1, last bullet)

**Proposition P (corrected after GPT-6.1-sol's two counterexamples).**
- **Proved, elementary.** For full windows that use all roots, the true pair is covered iff
  p mod M ∈ R_M(N) for some modulus M of the design. Here
  R_M(N) = {x : M | Bx² − AN for a group (v, M)} is computable from N. This holds because
  Bp² − AN = p·δ_v and gcd(p, M) = 1.
  - For Lehman or partial windows the window condition must be added, which gives a
    two-dimensional (residue, size) condition. Counterexample to the window-free form:
    N = 187, M = 3, A = B = 1, window {4} misses g = 3.
- **The natural public certificate (sufficient, not necessary).**

    𝒞(N): every integer x ∈ [√(N/2), √N] has M | Bx² − AN for some group (v, M).

  - 𝒞(N) implies coverage of the true pair. Coverage does not imply 𝒞(N): N = 91, M = 3
    covers (7, 13), but x = 9 fails.
- **Assessments (not theorems).**
  - 𝒞(N) is a complete covering of an interval by N-dependent quadratic classes
    ±√(AN/B) (mod ℓ). Its heuristic failure count is ≈ P e^{−λ}.
  - The short-interval theorems concern the classes 0 (mod ℓ), i.e. multiplicative
    structure of integers in intervals: Jutila; Friedlander–Lagarias Lemma 2.2; Harman's
    prime-detecting sieves; Baker–Harman; Ford.
  - Sieve methods prove non-emptiness of sifted intervals (Jacobsthal-type, Iwaniec), not
    emptiness at this density.
  - We found no design whose coverage reduces to an unconditional short-interval theorem.
  - The (x, x+1) route of B3 yields U(x, H) only as an *implied* statement (a barrier).

## 6. Conditional results (orchestrator item 4)

| hypothesis | what it gives here | verdict |
|---|---|---|
| Hooley R* (integer variable: Σ_{n∈I} e(an̄/q) ≪ (a,q)^{1/2}\|I\|^{1/2}q^ε; Iwaniec 1982 p. 110) | for integer moduli, pointwise κ = √Y q^ε ⇒ coverage at λ ≳ √J via Erdős–Turán. With one dilate, H ≳ √Y | **this use is insufficient** (Theorem U needs λ < J^{1/4}). Not a proof that no other consequence of R* could help |
| GRH (Dirichlet L) | character sums over primes, bounded linear statistics. The multiplicative route needs product size q^{2}. Korolev's slides: GRH gives nontrivial Kloosterman sums over primes only for X ≥ q^{1+ε} | insufficient |
| Elliott–Halberstam, Bombieri–Vinogradov, Montgomery AP | averaged distribution of primes in progressions. Here one fixed pair's points x_ℓ are needed; the relevant set {δ_v} has only H² elements | parameter mismatch; no deduction found |
| Cramér-type gaps for 𝒮_Y = {n : n has a prime factor in (Y, 2Y]} | applies to contiguous intervals. The relevant sets are rank-2 GAPs (or APs of difference q) | needs a non-standard GAP/AP version, which ≈ DCH_pair |
| Erdős–Odlyzko–Sárközy | products of two primes ≤ m in every class. Our box is ℓ ≍ q^{0.32}, k ≲ q^{0.72} | parameter mismatch; insufficient |
| "equidistribution of t mod ℓ" | sufficient only at Poisson precision on sets of measure ≍ J^{−3/4} (§3.2c) | that is K(θ) |
| P = BPP (complexity) | padding + Lenstra–Pomerance ⇒ deterministic L_N[1/2, O(1)] factoring | true but unrelated to coverage (GPT-6.1-sol) |

- **Literature check.**
  - Schoof 1982 (GRH) gives a deterministic N^{1/5+ε}.
  - Umans–Wang (arXiv 2511.10851, Thm 5.5) give 1/6 under SPDC(1/3, 1/3), and < 1/5 for
    max(α, β) < 2/5. GPT-6.1-sol notes a quantifier issue: their conjectures are stated for
    infinitely many sizes.
  - Carella 1308.2891 claims 1/6. The project's r13 referee judged its proof incomplete.
  - We found no deterministic sub-1/5 factoring theorem under a standard number-theoretic
    hypothesis, and the route here does not produce one. The table records failed
    deductions, not proved non-implications.

## 7. Experiments (pre-registered; notes/coverage38_prereg.md)

All runs used the code frozen at 757af72. The analysis is `src/coverage38/analyze.py`.
- Post-run fix (GPT-6.1-sol): the analyzer now checks grid completeness before declaring a
  pass. The verdicts are unchanged and all grids are complete.

### 7.1 e2 — exhaustive dichotomy worst case (S_32: 6 185 272 pairs; S_36: 78 295 606 pairs)

Band: 440 primes ≡ 3 (mod 4) in (2^13, 2^14]. Signed weights plus (1,1).

| n | J | G | easy frac | h*_all | h*_hard | h_pred(hard) | tail enrichment (easy) |
|---|---|---|---|---|---|---|---|
| 32 | 16 | 4 | 0.234 | 67 | 67 | 65 | 0.98 |
| 32 | 64 | 4 | 0.234 | 38 | 38 | 39 | 1.01 |
| 32 | 128 | 4 | 0.234 | 29 | 29 | 29 | 1.06 |
| 32 | 440 | 4 | 0.234 | 19 | 19 | 18 | 1.02 |
| 32 | 16 / 64 / 128 / 440 | 16 | 0.050 | 67 / 38 / 29 / 19 | 67 / 38 / 29 / 19 | 65 / 40 / 29 / 19 | 0.97 / 0.90 / 0.84 / 1.02 |
| 36 | 16 / 128 / 440 | 4 | 0.256 | 71 / 33 / 19 | 71 / 33 / 19 | 68 / 31 / 19 | 1.00 / 1.00 / 1.00 |
| 36 | 16 / 128 / 440 | 16 | 0.057 | 71 / 33 / 19 | 71 / 33 / 19 | 68 / 31 / 19 | 0.98 / 0.99 / 1.01 |

- **D1 passes**: the easy-fraction ratio is 0.214 (n = 32) and 0.222 (n = 36), against the
  model's 0.18 and the range [0.10, 0.35].
- **D2 passes**, with Δ = 0 in all 14 cells.
- **D3 passes**: |h*_hard − h_pred| ≤ 3.
- **D4 passes**.
- **Decision: the Lehman dichotomy, with the pre-registered uniform-window score, does not
  change the exhaustive worst case (experimentally supported).** h* is attained by a
  Lehman-hard pair in every cell, although in three cells the recorded argmax is Lehman-easy.

### 7.2 e2s — sampled, P ≈ Y³ (Y = 2^{n/6}; 200 000 pairs per n)

| n | J | G | easy frac | H_ref | unc (all) | unc (easy) | unc (hard) | R | z |
|---|---|---|---|---|---|---|---|---|---|
| 60 | 71 | 4 | 0.251 | 7 | 0.0324 | 0.0321 | 0.0325 | 0.987 | −0.45 |
| 60 | 71 | 16 | 0.044 | 7 | 0.0324 | 0.0315 | 0.0325 | 0.968 | −0.55 |
| 72 | 231 | 4 | 0.244 | 8 | 0.0335 | 0.0331 | 0.0336 | 0.986 | −0.50 |
| 72 | 231 | 16 | 0.041 | 8 | 0.0335 | 0.0339 | 0.0335 | 1.013 | +0.22 |

- S1 and S2 pass. The Poisson rate at H_ref is 0.0315 / 0.0330.
- Power: with ≈ 8–9k easy pairs at G = 16, SE(R) ≈ 6%. So S1 excludes only effects of
  ≳ 25–30% at G = 16, and ≳ 10% at G = 4 (SE ≈ 3%).
- **Decision: the failure sets of Lehman and coverage are not disjoint at P ≈ Y³; no
  dichotomy gain (experimentally supported).**

### 7.2′ Exploratory (not pre-registered): the corrected Lehman score L*

`src/coverage38/explore.py`, `experiments/raw/coverage38_explore.jsonl`. "Easy" here means
L* ≤ P^{2/5}/G; max L*/P^{2/5} = 0.79 / 0.91 / 0.95 at n = 32 / 60 / 72.

| set | G | easy frac | result |
|---|---|---|---|
| S_32 exhaustive, J = 16 / 64 / 128 / 440 | 4 | 0.29 | h*_hard = h*_all = 67 / 38 / 29 / 19 |
| S_32 exhaustive, same J | 2 | 0.72 | h*_hard = 67 / 37 / 28 / 17 against h*_all = 67 / 38 / 29 / 19, with Poisson h_pred(hard) = 64 / 38 / 29 / 17 |
| n = 60 sample | 2 / 4 | 0.79 / 0.36 | R = 1.014 [0.987, 1.041] / 1.023 [0.982, 1.063] |
| n = 72 sample | 2 / 4 | 0.79 / 0.37 | R = 0.931 [0.906, 0.955] / 0.941 [0.904, 0.979] |

- The G = 2 reduction at n = 32 matches the smaller Poisson extreme of a population shrunk
  to 28%. It is a sample-size effect, not structure.
- At n = 72, Lehman-easy pairs are ≈ 6% *less* often uncovered. A dichotomy removes the pairs
  that coverage already handles slightly better.
- No dichotomy gain.

### 7.3 e3 — the linear-statistics barrier on real configurations (n = 48, J = 1612)

Three primes q ≈ 1.4–1.6·10⁷, with 1000 primes p per q. The band is all 1612 primes in
(2^14, 2^15]. k₀ is tuned so that λ_q ≈ 3. Each cell also has 100 zero-hit random-complement
surrogates (Prop. L).

| q | H′ | k₀ | λ_q | uncovered (of 1000) | min E_abs/λ | E_abs/E_pred | unc/cov: E_abs, M_low | surrogate/real: E_abs, M_low, M_all | M_all/√J |
|---|---|---|---|---|---|---|---|---|---|
| 13925123 | 4 | 4798 | 2.96 | 57 | 113 | 1.000 | 1.000, 0.998 | 1.001, 1.006, 1.008 | 4.01 |
| 13925123 | 16 | 1315 | 3.00 | 39 | 282 | 1.000 | 1.000, 1.008 | 1.000, 0.997, 1.000 | 4.01 |
| 14925149 | 4 | 5141 | 2.96 | 42 | 114 | 1.000 | 1.000, 1.003 | 1.000, 1.000, 0.972 | 4.15 |
| 14925149 | 16 | 1410 | 3.00 | 43 | 285 | 1.000 | 1.000, 1.002 | 1.000, 1.002, 0.970 | 4.15 |
| 15925153 | 4 | 5482 | 2.96 | 58 | 115 | 1.000 | 1.000, 1.033 | 1.001, 1.011, 0.938 | 4.29 |
| 15925153 | 16 | 1504 | 3.00 | 58 | 287 | 1.000 | 1.000, 1.005 | 1.000, 1.001, 0.937 | 4.29 |

- **O1–O5 all pass.**
  - O1: the absolute Erdős–Turán error is 113–287 times λ, so it certifies no pair.
  - O2: E_abs matches the random-point prediction to 10⁻⁴.
  - O3: uncovered and covered pairs agree to within 3.3%.
  - O4: surrogates match real configurations within 7%, and 100% of surrogates satisfy Prop.
    L's bound.
  - O5: pooled uncovered counts against e^{−λ}: 157/155.3 (H′ = 4) and 140/149.3
    (H′ = 16).
- One small systematic effect: real configurations have a 3–7% *larger* M_all than the
  surrogates for two of the three q.
- Sampling disclosure: p = prevprime(uniform integer), i.e. gap-weighted rather than uniform
  over primes (GPT-6.1-sol).
- M_all = max_{h≠0}|K_p(h)| does not depend on p. Since K_p(h) = K_1(hp), it is a function
  of q and the band only, ≈ 4√J here. So its failure to distinguish covered from uncovered
  pairs is algebraic. Any bound uniform in the numerator a is blind to p.
- **Decision: pointwise exponential-sum information does not distinguish uncovered from
  covered pairs (experimentally supported). The barrier of Prop. L is realised by the actual
  configurations p·ℓ̄ mod q.**

## 8. Candidate card

- **Name.** Universal (λ, J) cost bound and linear-statistics barrier for DCL worst-case
  designs.
- **Proved.**
  - In the cost model, with the per-group charge and WIN for full or stationary windows:
    Theorem U and Corollaries U1–U4. These give the threshold λ⁴ < J for every design,
    attained up to constants by one band, the exponent (1+θ)/(6+θ), and the
    N^{1/6}σ^{2/3} floor. The proof architecture is close to Ward–Maestro (prior art).
  - Prop L, an elementary concentration lemma, with complements (a)–(c) in their restricted
    form.
  - Prop D: the operational Lehman side, linear families easy, B3 circumvented for
    dichotomies, B3's non-linear survivor, and the corrected Dirichlet bound for L*.
  - Prop P in corrected form.
  - The class lock and Theorem R's extraction extend to squarefree moduli.
- **Experimentally supported.**
  - The dichotomy gives no worst-case change (e2) and no rate change at P ≈ Y³ (e2s). The
    L* exploration even shows a slight opposite effect.
  - Pointwise statistics of real configurations are indistinguishable between covered and
    uncovered pairs, and from zero-hit random surrogates (e3).
- **Refuted for this model.**
  - Any design of the family needing a weaker (λ, J) threshold than λ < J^{1/4}.
  - Phase-discarding pointwise exponential-sum arguments, including R* used this way, as a
    source of gain.
  - The uniform-window Lehman score as a Dirichlet-bounded quantity (our own first draft).
- **Conjectured.** K(θ) for some θ < 1/4: true heuristically (Poisson), not implied by any
  standard hypothesis we found.
- **Heuristic.** Composite moduli are worse via Ford clustering; dead weights are O(1).
- **Open (most promising question).** A uniform positivity count of inverses of band primes,
  p ℓ̄ mod q, in a Farey-arc set of measure J^{−3/4}. It must use the arithmetic of
  {ℓ̄ : ℓ prime ∈ (Y, 2Y]} beyond random-level pointwise bounds. Equivalently: one of
  H ≈ q^{1/25} dilates of every residue class mod q is a product of a prime in
  (q^{8/25}, 2q^{8/25}] and an integer ≤ q^{18/25}.

## 9. Reviews

**Opus critic (read-only subagent): ACCEPT WITH FIXES.**
- *Verified.* The Theorem U proof (all steps; the height count even holds with 2h²), WIN for
  full windows including composites (minimum S-range 0.0253hP), both attainment optima,
  U1–U3, Prop L's constants, §3.1, the composite class lock, D(a) and D(d), and the §7.1–7.2
  tables against the raw files.
- *Fixes applied:*
  - **(major) D(b) false**: replaced by the L* score and bound, with a test;
  - **(major)** e3 not to be claimed before its data, now filled;
  - **(major)** the Prop L Consequence restricted to (i) centred estimates with (ii)
    tolerance ≥ 2λ₀ + 2√(J log 8Q);
  - ‖Ŝ‖₁/q corrected;
  - the Bernstein remark relabelled heuristic, and "almost as small as targets" withdrawn;
  - Prop P scoped to full windows, with "the only" changed to "the natural";
  - the composite extension's extra hypotheses;
  - WIN constants (Brun–Titchmarsh range, log log loss);
  - WIN for arbitrary placements declared a hardness assumption (magic window);
  - "extremal" changed to "attains up to constants";
  - U1's "weakest" withdrawn;
  - Ford and dead weights labelled heuristic;
  - the §7.1 argmax sentence;
  - e2s power stated;
  - the tautological test disclosed;
  - D(c) cost made explicit.

**GPT-6.1-sol adversarial review (`codex/coverage38_2.out`): REJECT of the first draft's
universal reading; U and L as conditional/elementary results judged sound.**
- *Fixes applied:*
  - the charging convention made explicit;
  - WIN restricted to specified window rules, with the magic-window obstruction;
  - single-band uniqueness withdrawn, since stationary windows have the (tM, tΔ) symmetry;
  - composite-moduli costs added (factorisation, roots, rejected groups; N = 989, M = 21
    example);
  - "Harman's sieve not needed" narrowed to "prime moduli unnecessary for extraction";
  - the Prop L Consequence replaced by the restricted statement, including the arithmetic
    statistic Σe_q(ℓ_j x_j) = J e_q(p) that generic configurations fail;
  - Prop P's necessity claim deleted, with two counterexamples (N = 187 window; N = 91
    full-candidate);
  - R* scoped to integer variables;
  - citation ranges corrected (Baker, Fouvry–Shparlinski, Korolev 2019, Irving,
    Bourgain–Garaev title);
  - EH/BV wording changed to "parameter mismatch";
  - K(0) mean-one counterexample, and K(θ) restated;
  - the conditional table relabelled as failed deductions;
  - the analyzer completeness check;
  - prevprime sampling disclosed;
  - **prior art: Ward–Maestro (Sept 2026, igent.ai report §§3, 13) has the same
    Cauchy–Schwarz + reciprocal-height architecture and the 1/5 and 1/6 floors.** We have
    not read it ourselves. Our contribution is limited to the multi-modulus Hölder step and
    the (λ, J) formulation: "no prior formulation found by us" of that exact form.
- *Not resolved.* A fully arithmetic WIN for arbitrary efficiently computable windows. This is
  out of reach and is stated as an assumption.

**GPT-6.1-sol blind round (`codex/coverage38_1.out`).** It independently found:
- Theorem U's inequality, (λ⁴ ≲ JTS²/P²);
- the oracle construction;
- linear families being easy;
- the (x, ⌊αx⌋) survivor;
- the zero-centred correction;
- the Korolev 2016 and Bourgain–Garaev bounds;
- the P = BPP remark.

## 10. Reproduction

```
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_coverage38.py
PYTHONPATH=. .venv/bin/python -m src.coverage38.experiments e2    # ~25 min, ~7 GB (n = 36)
PYTHONPATH=. .venv/bin/python -m src.coverage38.experiments e2s   # ~1 min
PYTHONPATH=. .venv/bin/python -m src.coverage38.experiments e3    # ~2 h (FFT of prime length ~1.4e7)
PYTHONPATH=. .venv/bin/python -m src.coverage38.analyze
PYTHONPATH=. .venv/bin/python -m src.coverage38.explore     # exploratory L* rerun, ~5 min
```
