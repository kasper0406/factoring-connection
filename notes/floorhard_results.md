# floorhard: hardness side of the Theorem L_j floor-sum oracle (phase 17)

Agent `floorhard`, 2026-10-20. This is the revised version after the codex adversarial round and two
Opus critics (§9).

> **Erratum (2026-10-27): ETH.** OpenAI, "The Exponential Time Hypothesis is False" (preprint, 17 Oct 2026, rev. 25 Oct 2026), https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/paper.pdf gives a deterministic #k-SAT algorithm in time 2^{O(n log log n/√log n)}. This refutes ETH, randomized ETH, #ETH and Gap-ETH (unrefereed; Theorem 1.1 there). Every ETH-conditional statement below is therefore restated under **ETH⁻**: *for every ε > 0, 3-SAT on n variables with O(n) clauses has no 2^{O(n^{1−ε})}-time algorithm* (randomized ETH⁻ for randomized reductions). On such instances the paper's bound is 2^{O(n/√log n)} = 2^{n^{1−o(1)}}, so it is consistent with ETH⁻. ETH⁻ is stated for O(n)-clause instances because sparsification costs 2^{εn}, which no power-saving hypothesis can absorb. Conclusions here lose an ε in the exponent. Exact degree-2 floor sums have no deterministic 2^{L^{1/2−ε}} algorithm, for any ε > 0, unless ETH⁻ fails; the earlier statement was no 2^{o(√L)} under ETH. The same applies to the Jansen–Klein–Lassota variant and to the uniform-oracle consequence in §4.3. The P = NP and FP = #P consequences are unaffected.


**Files.**
- Pre-registration: `notes/floorhard_prereg.md`, commit `6b3a37e`, made before any code existed.
- Code: `src/floorhard/{core,run_experiments}.py`. Tests: `tests/test_floorhard.py` (10 pass).
- Raw data: `experiments/raw/floorhard_e{1,2,2c,3}.jsonl`. Summary:
  `experiments/analysis/floorhard_summary.json`.
- GPT-6.1-sol: blind round `codex/floorhard_1.{md,out}`; adversarial round `codex/floorhard_2.{md,out}`.

**Labels.** KNOWN (with reference), PROVED (proof here; all elementary; at most "no prior formulation
found by us"), HEURISTIC, CONJECTURED, OPEN, experimentally supported, UNVERIFIED.

**Consistency with `floorsum`.** floorsum (`notes/floorsum_prereg.md`, algorithm side) has two claims
that agree with §4 here:
- C1: the per-block convexity scale sums to N^{1/3} for every j;
- C5: Coppersmith-type lattices detect divisors in j = 2 blocks.

Nothing here depends on floorsum's experimental outcomes.

**Oracle classes.** "𝒪_j, uniform" means one algorithm that is correct for every N and every query in
the class of Theorem L_j. The class consists of:
- Taylor polynomials P_j of N/(x+h) at any integer x with N^{1/(j+2)}/c ≤ x ≤ √N;
- arbitrary rational shifts;
- subintervals of [0, ⌊cx/B⌋).

"𝒪_j^top(κ)" is the same algorithm promised only x ∈ [√N/κ, √N]: a uniform *promise* problem. Costs
T(H) are per query, as a function of the queried interval length H. Polynomial dependence on the bit
length is absorbed into T.

## 0. Verdict (plain language)

**The oracle hypothesis of Theorem L_j, as stated, is not a plausible target, for any j ≥ 2.** It asks
for exact Σ⌊P(h)⌋ in polylog time for all rational P of degree ≤ j, with shifts and subintervals.
- That problem is #P-hard and NP-hard.
- It stays hard when restricted to the ladder's own Taylor polynomials, as long as the algorithm must
  work for every N (or, with randomized conclusions, every prime N).
- A uniform oracle with per-query cost T(H) also factors any m directly, in O(log² m · T(m)).
- For power-law T(H) = H^δ this direct route beats the ladder whenever δ < 2/(j+4).

**But the hardness does not reach the part of the problem that matters.**
- Every hard instance we or codex could build sits in blocks at x ≤ N^{3/8+o(1)}, for specially
  structured N.
- A balanced semiprime needs only the top blocks x ∈ [√N/κ, √N]. There:
  - the j = 2 arcs wrap only O(κ) times, and Coppersmith–May already detects divisors in them;
  - for j ≥ 3 the arcs lie in an intermediate regime, beyond Coppersmith's range for windows yet far
    shorter than any known NP encoding needs.
- All attempts to embed hard problems in top blocks failed: square padding, rational-slope frames,
  and continued-fraction frames (us and codex).
- Exact-value queries at top blocks are within univariate Coppersmith bounds (§4.2 T3′).
- The NP/#P/ETH families themselves are solvable in H^{o(1)}. So they exclude no positive-power cost
  T_j = N^ε, which is all a fixed-j ladder needs.

**Where the hardness lies.** In uniformity over structured N, together with arcs that wrap many times
around a modulus with many known prime factors.
- For the top-block promise, no reduction from any believed-hard problem is known. Its only known
  consequence is Theorem L_j itself: balanced-semiprime factoring in N^{1/(j+2)}·T_j. That is
  "factoring with (1/(j+2))·log₂N known top bits of p", beyond Coppersmith's ¼ for j ≥ 3.

**Plausibility (HEURISTIC, personal judgement).** We put a deterministic top-block oracle with
T_j = N^{o(1)}, at some fixed j ≥ 4, at about 3–5%. Codex estimates about 5% (range 1–20%).
- The NP-hardness is not evidence here.
- The weight comes from the lack of any exact-counting mechanism. Fourier methods need ≍ H terms; hull
  walks only redistribute the global N^{1/3}; one Euclid step turns the arc into a general conic.

## 1. Basic identity

**Lemma W (PROVED).** Let f be an integer polynomial, D ≥ 1, and 0 ≤ u < v ≤ D. Then
#{h ∈ [h₁, h₂) : f(h) mod D ∈ [u, v)} = F(f − u) − F(f − v), where F(g) = Σ_{h₁≤h<h₂}⌊g(h)/D⌋.
- For real P and 0 ≤ α < β ≤ 1: #{h : {P(h)} ∈ [α, β)} = Σ(⌊P − α⌋ − ⌊P − β⌋).
- With shifts, a floor-sum oracle **is** a residue-distribution oracle. Codex (blind (1)–(2)) states the
  same.

## 2. Reductions from fast degree-2 floor sums

### 2.1 Arbitrary intervals

**Theorem SS (PROVED; #SubsetSum → interval root counting, parsimonious).**
- Given w₁, …, w_k ≥ 1 (k ≥ 1) and 0 ≤ T ≤ S = Σw_i, one constructs deterministically in polynomial
  time a, b and ≤ 2 intervals I such that #{x ∈ I : x² ≡ a (mod b)} = #{δ ∈ {0,1}^k : Σδ_iw_i = T}.
- For T ∉ [0, S] the count is 0. **Without this check the construction aliases mod M.** Counterexample
  (critic): w = [1, 2], T = 9 returns 1.
- The bit length is L = O(k log(kS)).

*Construction.*
1. **Moduli.** M = 4S + 4. Take ℓ_i = the i-th odd prime and q_i = ℓ_i^{e_i}, the least power exceeding
   8kM. b = ∏q_i. (This deterministic repair is codex's; the first draft searched for primes > kM,
   which is not a proved poly-time step.)
2. **The θ_i.** u_i = round(q_iw_i/M), plus 1 if ℓ_i | u_i. Then 1 ≤ u_i < q_i/4 and
   |u_i/q_i − w_i/M| ≤ 3/(2q_i). θ_i = (b/q_i)u_i.
3. **The roots.** a = (Σθ_i)² mod b. Units mod odd prime powers have exactly two square roots, so the
   roots of x² ≡ a (mod b) are x_ε = Σε_iθ_i mod b, ε ∈ {±1}^k, all distinct.
4. **Their positions.** x_ε/b ≡ Σε_iw_i/M + η (mod 1), with |η| < 3k/(16kM) < 1/(2M).
5. **Separation.** The values Σε_iw_i = 2Σδ_iw_i − S share a parity and lie in [−S, S] ⊂ (−M/4, M/4).
   Distinct values are therefore ≥ 2/M apart on the circle, without wrap. The window
   [(2T′−1)/(2M), (2T′+1)/(2M)), with T′ = 2T − S, lies in (−¼, ¼) for 0 ≤ T ≤ S. It contains exactly
   the roots with Σδ_iw_i = T, and no root lies on an endpoint.
6. **The floor sums.** The count is 2 or 4 floor sums of x² − a and x² − a − 1 over D = b (Lemma W). ∎

- E2a (pre-registered, prime moduli): 40/40. E2b (CRT core, k = 4–10): 200/200.
- E2c (POST-HOC, deterministic prime-power variant with out-of-range targets): 300/300.

**Corollaries.**
- **(i) #P (PROVED from KNOWN pieces).** #SubsetSum is #P-complete under parsimonious reductions.
  - The CLRS/Sipser 3SAT → SubsetSum reduction, with clause target 4 and slacks {1, 2}, is
    parsimonious: 1, 2 or 3 true literals force slack 3, 2 or 1, each uniquely. There are no carries.
    #3SAT is #P-complete (Valiant 1979).
  - Exact degree-2 floor sums are therefore #P-hard under constant-query (subtractive) Turing
    reductions. Deterministic polynomial time would give FP = #P, hence P = NP.
- **(ii) ETH⁻ (PROVED from KNOWN; v1 used ETH, which is refuted — see the erratum).** ETH implies that Subset Sum with k items and weights 2^{O(k)} has no
  2^{o(k)} algorithm. The standard reduction from 3SAT with O(n) clauses gives Subset Sum with k = O(n) items, with no
  sparsification step, so ETH⁻ ⇒ no 2^{k^{1−ε}} algorithm (Etscheid–Kratsch–Mnich–Röglin, JCSS 2017, Lemma 3; Jansen–Land–Land, SIAM J.
  Discrete Math. 30 (2016)).
  - Here L = O(k²). So exact degree-2 floor sums have no deterministic 2^{L^{1/2−ε}} algorithm (any ε > 0) unless
    ETH⁻ fails. (v1: no 2^{o(√L)} unless ETH fails.)
- **Prior work (KNOWN).** NP-completeness of bounded quadratic congruences is Manders–Adleman (JCSS 16
  (1978) 168–184; Garey–Johnson AN1), including with the factorisation of b given.
  - Jansen–Klein–Lassota (Math. Program. 197 (2023), arXiv 2008.12928) give a squarefree known-modulus
    version with O((n+m)² log(n+m)) bits. That yields ETH⁻ ⇒ no 2^{L^{1/2−ε}} (v1: ETH ⇒ no 2^{o(√(L/log L))}).
  - New here: the parsimonious/#P form, the interval-count form through floor sums, and the dropped log
    factor. Manders–Adleman's own reduction has log b = Θ(n² log n) or more (codex, lit critic). **No
    2^{Ω(L)} bound follows from any known construction.**

**Prop BIS (folklore; the counting analogue of Rabin 1979).** Let m be odd with ω(m) ≥ 2.
- Root counts of x² ≡ 1 (mod m) on subintervals (two calls each) bisect [2, m − 1) down to a
  nontrivial root r. Then gcd(r − 1, m) splits m.
- That is O(log m) calls of length ≤ m for one split. Complete factorisation takes O(log² m) calls,
  plus perfect-power and primality (AKS) tests.
- Codex found this independently.

### 2.2 Full periods: class numbers, residuosity

**Prop FP (PROVED; the prime case is classical: Dirichlet, Berndt 1976, Hudson–Williams 1982).** For
odd squarefree m and gcd(a, m) = 1,

  Σ_{h=0}^{m−1} ⌊ah²/m⌋ = a(m−1)(2m−1)/6 − (m−1)/2 + Σ_{d | m, d ≡ 3 (4)} (a/d)·h*(−d),

with h*(−3) = 1/3 and h*(−d) = h(−d) for d > 3.

*Proof.*
1. #{h : ah² ≡ r} = Σ_{d|m}(ār/d) (CRT).
2. Σ_{t<m} t(t/d) = (m/d)Σ_{s<d} s(s/d). This is 0 for even characters and −m·h*(−d) for odd ones.
3. Subtract from aΣh² and divide by m. ∎

- Checked: E1a, 6 520 pairs with 0 failures. The critics re-checked 1 049 + 27 092 pairs, up to
  ω = 5, also 0 failures.
- Codex's blind (6) is the same identity. Squarefreeness is needed.

**Corollaries.**
- **(FP1).** For prime p ≡ 3 (mod 4), one call gives h(−p).
  - Shifted full periods give the incomplete sums Σ_{t≤X}(t/p) (codex (4)).
  - For p ≡ 1 (mod 4) they give h(−4p) = 2Σ_{t<p/4}(t/p).
  - Best known for h: L[½] under GRH (Hafner–McCurley, J. AMS 2 (1989); Biasse–Erukulangara, Adv.
    Math. Commun. 17 (2023)). Super-polylog hardness of h(−p) is CONJECTURED, not proved.
- **(FP2) Class-number decoder (PROVED).** Let m = pq, p ≠ q odd primes > 3, and sample a uniformly
  mod M subject to the character conditions below. Write C_M(a) for the class part.
  - **m ≡ 3 (mod 4):** take M = m and (a/m) = 1. Then C_m(a) = h(−m) + s·h(−q), where q is the factor
    ≡ 3 (mod 4) and s = ±1.
  - **m ≡ 1 (mod 4):** take M = 3m and (a/3) = (a/M) = 1. Then C_M(a) − 1/3 = h(−M) + s·X, where
    X = h(−p) + h(−q) if p ≡ q ≡ 3 (mod 4), and X = h(−3p) + h(−3q) if p ≡ q ≡ 1 (mod 4). The draft's
    formula was wrong (codex, critic).
  - Since X > 0, exactly two values occur, each sign with probability ½. Their mean is h(−M).
- **Splitting M.** h(−4M) = h(−M)(2 − (−M/2)) holds for M > 3. Splitting M from h(−4M) by ambiguous
  forms is KNOWN (Shanks 1971; Schoof 1982, MC Tract 154) as an algorithm. A rigorous polynomial
  running-time bound needs assumptions on sampling class-group elements (GRH-type). That step is
  **experimental / conditional** here.
  - E1b (pre-registered; a sampled from [2, 1000), as run): h(−M) correct 30/30, median 3 calls (max 6),
    m split 29/30.
- **(FP3).** The window [0, 1) over a full period counts the square roots of c, ∏(1 + (c/p)). This
  decides quadratic residuosity mod a composite (Goldwasser–Micali, JCSS 28 (1984)).
  - Windows [0, Y) give Σ_{d|m}(ā/d)·(incomplete Jacobi sums mod d): selection of local characters.

**Ranges used by our reductions.** All of them use intervals covering a large fraction of the period
of h² mod b:

| Reduction | Interval |
|---|---|
| Thm SS (CRT encoding) | ≈ b/M ≈ b^{1−O(1/k)} |
| BIS | subintervals of [0, m) |
| FP | full period |

Short *shifted* intervals are not automatically trivial (11² = 1 + 8·15 with I = [11, 12)). The
relevant tool there is Coppersmith's univariate theorem for shifted polynomials.

## 3. The ladder's restricted class

**Theorem TE (PROVED).** Fix j ≥ 2, A ≥ 1, b ≥ 2, c′ ≥ 1, and an integer t with t ≥ 2c′/b,
t > 8Ac′³/b and t ≥ (4/c)⁴Ac′⁴/b². Put N = Ab²t³ and x = bt. Then:
- [x, x + c′) lies inside the ladder block [x, x + ⌊cx/B⌋), and x is admissible;
- for 0 ≤ h < c′, P_j(h) = Abt² − Ath + Ah²/b + R(h), with |R(h)| < 1/(4b);
- hence, for every r ∈ [0, b) and I ⊆ [0, c′), #{h ∈ I : Ah² ≡ r (mod b)} =
  #{h ∈ I : {P_j(h) + 1/(4b)} ∈ [r/b, r/b + 1/(2b))}. That is two queries.

*Proof.*
1. P_j(h) = Σ_{i≤j}(−1)^i A h^i t^{2−i}/b^{i−1}.
2. The terms i ≥ 3 total at most (Ac′³/(b²t))/(1 − c′/(bt)) ≤ 2Ac′³/(b²t) < 1/(4b).
3. The block condition: x/B ≥ b^{1/2}t^{1/4}A^{−1/4} for j ≥ 2. ∎

**Variants.**
- **Prime N (PROVED; randomized).** Let N = Ab²t³ + e with 0 ≤ e < E = bt²/(16c′), and t > 16Ac′³/b.
  Queries subtract the constant e/x; the remaining perturbation varies by < 1/(8b).
  - **Correction:** E = N₀^{2/3}/(16A^{2/3}b^{1/3}c′), where N₀ = Ab²t³, is *not* ≍ N^{2/3}. At the
    minimal t ≈ b² it is only ≍ N^{1/2}.
  - Take A = 1, c′ ≤ b, t ≍ b³. Then E ≍ N^{6/11} > N^{0.525}, and Baker–Harman–Pintz (Proc. LMS 83
    (2001)) give ≫ E/log N primes in the range. Random e plus AKS find one in expected O(log N) trials.
    In this case x ≈ N^{4/11}.
  - Prime-only polylog evaluation gives zero-error randomized consequences: NP ⊆ ZPP, and randomized
    exact counting. It does not give deterministic FP = #P. The ETH⁻ conclusion needs randomized ETH⁻.
- **Codex's one-call embedding (PROVED in codex/floorhard_1.out (7)–(11); critic re-checked 852
  cases).** Take x = 16m⁴ (use x = 16m⁴⌈c^{−4}⌉ for small c) and N = (ax³ − bx²)/m + 1. Then
  gcd(N, x) = 1, and one query computes any Σ_{h<m}⌊(ah² + bh + c)/m⌋.
- **Block starts only.** If only the ladder's block starts x_i are allowed, re-centre: the extra error
  is ≤ ρ_i < 1/(8b) once t ≳ c^{12}Ab² (critic; a sketch, not written out).
- **Experiments.** E3a 60/60 (j = 2, 3, 4); E3b 20/20 prime N. These were pre-registered at small
  sizes, where the code's t sufficed. E3c: 40/40 composites ≤ 20 000 factored through embedded windows
  (≤ 32 calls).

**Corollary DOM (PROVED; the slogan of the draft is withdrawn).**
- If 𝒪_j (uniform) costs T(H) per query, deterministically, then every m is completely factored in
  O(log² m · T(m)) (BIS through TE, N ≈ m⁸).
- Theorem L_j gives N^{1/(j+2)}·T(N^{1/2−1/(j+2)}).
- For T = H^δ the direct route is strictly better iff δ < 2/(j+4); at equality the logarithmic factors
  decide. For larger δ the ladder is better: j = 2 and δ = ½ give ⅜ against ½.
- **Against Harvey's N^{1/5}:**
  - the ladder beats it iff δ < 2(j−3)/(5j); the direct route iff δ < 1/5;
  - for j ≤ 6, any uniform power-law oracle that lets the ladder beat N^{1/5} also lets the direct
    route beat it;
  - for j ≥ 7 there is a window 1/5 ≤ δ < 2(j−3)/(5j) where only the ladder does (critic).
- So uniform padding gives a stronger factoring reduction for small query exponents. The top-block
  promise avoids every known padding construction.

**Degree ≥ 3.** The class of degree ≤ j contains degree 2, and TE absorbs the terms i ≥ 3. So every
statement in §§2–3 holds for all j ≥ 2.

## 4. Where the hardness lives

### 4.1 Location of the constructed hard blocks (PROVED for these constructions)

- **TE.** With c′ = b^{1−o(1)} (as in the CRT encodings), t ≥ Ab^{2−o(1)}. So x = (Nb/A)^{1/3} =
  N^{3/8+o(1)} at minimal t, decreasing towards N^{1/3} as t grows. The prime variant gives N^{4/11}.
- **Codex.** x = N^{4/11+o(1)} for a = m^{o(1)}, tending to N^{1/3} for a ≍ m.
- **Scope.** This describes our constructions; it is not a theorem that hardness needs small x.
  Nothing excludes hard encodings with √b ≪ c′ ≪ b^{1−o(1)}.

### 4.2 Top blocks

For N = pq with p ∈ [√N/κ, √N], the ladder scans only x ∈ [√N/κ, √N].

- **(T1) j = 2: bounded wraps (PROVED).** H = cx/N^{1/4}, so N H²/x³ = c²√N/x ≤ c²κ. Thus ⌊P₂⌋ is a
  linear floor plus a count of {ℓ(h)} under O(κ) parabolic thresholds.
- **(T2) j = 2: detection is known (KNOWN + PROVED comparison).**
  - May's theorem: thesis 2003, Thm 7. Building on Coppersmith (J. Cryptology 10 (1997)) and
    Howgrave-Graham (1997), it finds a divisor ≥ N^β within distance c·N^{β²} in poly(log N, c) time.
  - With β = log x/log N ≥ ¼: β² − (β − ¼) = (β − ½)² ≥ 0, so every j = 2 block lies in that range.
    The j = 2 ladder reproduces Coppersmith's N^{1/4+o(1)}.
  - For j = 3 the comparison fails for β ∈ (0.276, 0.724).
  - This is detection, not exact floor sums.
- **(T3) j ≥ 3: no obstruction proved (draft claim REJECTED, now HEURISTIC remark).**
  - The quadratic part wraps ≈ c²N^{1/2−2/(j+2)} ≫ 1 times. For j ≥ 5 the cubic also wraps many times.
    These blocks sit in the intermediate regime √b ≪ H ≪ b^{1−o(1)}.
  - The draft's "CRT encodings need period fraction b^{−O(1/k)}" excludes only the SS positioning. It
    is not a general obstruction.
  - Codex's counterexample is square padding: y² ≡ a (mod b), y < c becomes z² ≡ q²a (mod bq²),
    z < qc. It makes a hard family occupy a short fraction of a larger period.
  - Codex could not push square padding, or rational-slope and continued-fraction frames, into the
    top Taylor class: the unwanted points between multiples of q are not excluded, and a unimodular
    image is a general conic. Neither could we.
- **(T3′) Exact-value queries at top blocks are Coppersmith-easy (PROVED given Coppersmith's univariate
  theorem; codex).**
  - Let D = x^{j+1}/gcd(N, x^{j+1}) ≥ x^{j−1}/κ². Fixing one fractional value of P_j is a monic
    degree-j congruence mod D (divide by N mod D).
  - H ≪ x^{j/(j+2)} ≤ x^{(j−1)/j} ≲ D^{1/j}, so all roots in the block are found in polynomial time,
    up to constant subdivision.
  - The ladder needs windows of positive width ρ ≈ 1/H. Those are the open part.
- **(T4) What the top-block oracle gives (PROVED one way; the draft's "equivalence" is withdrawn).**
  - The ladder (candidate location plus division) turns 𝒪_j^top with cost T into detection of a divisor
    in [x, x + cxN^{−1/(j+2)}): factoring a balanced N in N^{1/(j+2)}·T.
  - The draft equated windows of P_j with hyperbola strips. That holds for N/d, not for its Taylor
    polynomial. Counterexample (codex): N = 101·103, x = 100, h = 1. N/(x+h) = 103, while
    P₂(1) = 103 + 103·10⁻⁶, so the window [0, 1/20000) separates them.
  - Even for exact hyperbola windows the k-range [αd, βd) depends on d.
  - No converse reduction is known.

### 4.3 What the NP/#P/ETH hardness can and cannot exclude (PROVED)

- The sparse-3SAT image of Theorem SS has log q_i = O(k) and k ≍ √L. Its moduli are found by trial
  division in b^{O(1/k)} = b^{o(1)}, and its 2^k roots enumerated in b^{o(1)}.
- So that family is solvable in H^{o(1)}. More generally, with a supplied factorisation and unit
  targets, 2^{ω(b)} = b^{o(1)} (codex).
- This does **not** apply to the BIS family, whose hardness is factoring itself. That is exactly what
  DOM exploits.
- Hence these reductions provide no positive-power runtime exclusion: no "T_j ≤ H^δ ⇒ X believed hard"
  for any δ > 0. A hypothetical quasi-linear-size reduction would still give only 2^{o(L/polylog L)}.

## 5. Degree ≥ 3, cubic sums, cubic residuosity

- **(C1) (PROVED; codex blind (13)).** For squarefree m and gcd(a, m) = 1,
  Σ_{h<m}⌊ah³/m⌋ = am(m−1)²/4 − (m−1)/2, by the odd symmetry h → −h. This fails for non-squarefree m,
  e.g. m = 9. The pure cubic full period carries no class-number information, unlike degree 2.
- **(C2) Shifts (PROVED; corrected).** For integer 0 ≤ c < m,
  F(c) − F(0) = #{h : (ah³ mod m) ∈ [m − c, m)} = Σ_{r∈[m−c, m)} r₃(ā r).
  - r₃ expands into products of cubic characters of the prime factors ≡ 1 (mod 3). So shifted full
    periods give incomplete cubic character sums.
  - The window [0, 1) counts the cube roots of a unit: 0 or 3^s, where s = #{p | m : p ≡ 1 (mod 3)}
    (e.g. mod 55 every unit has one cube root; mod 35, 0 or 3). For m = pq with p ≡ q ≡ 1 (mod 3) this
    decides cubic residuosity.
- **(C3).** The complete Weyl sum of analytic's C3 is a Fourier coefficient of r₃. Floor sums give its
  interval sums. Codex's O(m)-query reconstruction is an upper bound only; no polylog reduction in either
  direction is known (OPEN).
- **(C4).** Q(√−3) has class number 1. We prove no reduction to class numbers of cubic fields or to
  Hecke L-values.
- **Answer to Task 2.** For j ≥ 3, polylog exact degree-j floor sums imply everything degree 2 implies
  (P = NP, FP = #P, factoring). No additional consequence is established here.

## 6. Literature: exact lattice-point counting under curves (hardness results we verified)

| Problem | Best known | Lower bounds / equivalences |
|---|---|---|
| D(x) = Σ_{n≤x}τ(n) | x^{1/3+o(1)}: Polymath / Tao–Croot–Helfgott, Math. Comp. 81 (2012) 1233–1246, §2.1 (Prop. 2.2 is x^{0.49+o(1)}). Sladkey, arXiv 1206.3369: O(x^{1/3}) time, O(log x) space. Li, arXiv 2609.29949 (Sept. 2026, unaudited preprint): N^{9/28+o(1)} expected, isolated full values, via Voronoi | none found. Hull vertices of {xy ≥ N}: Ω(N^{1/3}) and O(N^{1/3} log N) (Alcántara–Blanco–Criado–Santos, arXiv 2501.19193); Θ(N^{1/3} log N) (Balog–Bárány, arXiv 2602.06897). Both bound enumerating methods only |
| truncated Σ_{d∈[u,v)}⌊N/d⌋ | same | interval access factors N by bisection (folklore) |
| Σ_{n≤X} r₂(n) | O(X^{1/3} log X): Har-Peled, CGTA 10 (1998), plus Pick; disc hull Θ(r^{2/3}) (Balog–Bárány, SoCG 1991) | none; sector counts factor N (p ≡ q ≡ 1 mod 4) |
| π(x) | x^{2/3}: LMO, Math. Comp. 44 (1985). x^{1/2+ε}: Lagarias–Odlyzko, J. Algorithms 8 (1987). Õ(√x): Hirsch–Kessler–Mendlovic, arXiv 2212.09857 (journal venue UNVERIFIED) | none. TCH: parity of π on short intervals in x^{1/2−c} |
| Mertens M(x) | O(x^{3/5}(log x)^{8/5+ε}): Helfgott–Thompson, Res. Number Theory 9 (2023); Õ(√x): HKM | none |
| fixed-dimension polytopes | polynomial: Barvinok, Math. Oper. Res. 19 (1994) | — |
| parabola in a box (succinct coefficients) | — | NP-complete: Manders–Adleman 1978; JKL 2023. #P-hard: Thm SS |
| factoring baselines | Hiary, Math. Comp. 85 (2016), N^{1/3}; Harvey, Math. Comp. 90 (2021), N^{1/5}; Harvey–Hittmeir, arXiv 2105.11105, N^{1/5} with log-log savings; Lenstra–Pomerance, J. AMS 5 (1992), rigorous randomized L[½] | — |

**Answer to Task 3.**
- Neither we, nor codex, nor the literature critic found a conditional lower bound or a complexity
  equivalence for exact D(x), r₂ sums, or π(x).
- Li's "1/3 is characteristic of local linear methods" remark is cited from the abstract and
  introduction (UNVERIFIED in the body). It is explicitly not a lower bound.
- The hardness results we verified for lattice counting under algebraic curves are those listed in the
  table. We do not claim the list is exhaustive.

## 7. The sharp barrier statement

For fixed j ≥ 2:

1. **Uniform 𝒪_j (PROVED):**
   - in deterministic poly(log N) time: FP = #P and P = NP;
   - in deterministic 2^{(log N)^{1/2−ε}} time, for some ε > 0: ETH⁻ fails (v1: 2^{o(√log N)} ⇒ ETH fails; erratum);
   - in deterministic per-query time T(H): every m is factored in O(log² m · T(m)), which strictly beats
     the ladder iff T = H^δ with δ < 2/(j+4);
   - correct only for prime N: the randomized analogues (NP ⊆ ZPP; randomized ETH⁻).
2. **The constructed hard instances (PROVED for these constructions):**
   - they lie at x ≤ N^{3/8+o(1)};
   - the NP/#P/ETH families are solvable in H^{o(1)}, so they exclude no T_j = N^δ.
3. **Top-block promise 𝒪_j^top(κ) (OPEN):**
   - no NP/#P-hardness reduction and no factoring reduction stronger than the ladder is known;
   - exact-value queries are Coppersmith-easy (T3′), and j = 2 detection is known (T2);
   - a subpolynomial top-block oracle gives N^{1/(j+2)+o(1)} factoring for balanced semiprimes, not
     automatically for worst-case N.

**Bottom line.** The literal L_j hypothesis is an NP-hard target. The only part the ladder needs is a
uniform promise problem with no known hardness. Its plausibility rests on mechanisms, and we have
none: hull walks, Farey/Euclid, Fourier (barrier F) and lattice detection (T2) all fail beyond
j = 2.

## 8. Pre-registration outcomes

| Prediction | Result |
|---|---|
| P1.1: 0 failures of Prop FP (E1a) | ✓ 6 520 pairs |
| P1.2: decoded h(−M) correct 30/30 | ✓ |
| P1.3: factored ≥ 29/30 | ✓ 29/30 (miss at M = 3m: ambiguous forms returned only 3) |
| P1.4: median calls ≤ 4, max ≤ 12 | ✓ 3 / 6 |
| P2.1 (E2a) | ✓ 40/40 (32 with solutions) |
| P2.2 (E2b) | ✓ 200/200 (129 with solutions; b ≤ 181 bits) |
| P3.1 (E3a) | ✓ 60/60 |
| P3.2 (E3b) | ✓ 20/20 |
| P3.3 (E3c) | ✓ 40/40, ≤ 32 calls |

**Disclosures.**
- **E1b sizes:** m of 15–22 bits, not 16–22. The generator splits bits as ⌊n/2⌋ + ⌈n/2⌉.
- **Peek:** one smoke run before the experiments (m < 200, three subset-sum instances).
- **E3b prime variant:** it used the code's t (four times the TE minimum). That suffices at these sizes,
  but it is below the asymptotic t ≍ b³ required for the BHP argument (correction in §3).
- **Post-hoc changes, after review:**
  - E2c was added;
  - the code gained the T ∈ [0, S] check and the deterministic prime-power construction;
  - the decoder now samples uniformly mod M by default (E1b reproduces with `amax=1000`).
- All oracles are brute force: these runs test identities and reductions, never costs. Toy factoring is
  a correctness check.

Reproduction:
```bash
.venv/bin/python -m pytest -q tests/test_floorhard.py
for e in e1 e2 e3 e2c summary; do PYTHONPATH=. .venv/bin/python -m src.floorhard.run_experiments $e; done
```

## 9. Reviews and fixes applied

**Codex blind (floorhard_1).** Independently found:
- the NP-hardness (Manders–Adleman, two calls) and BIS;
- a one-call Taylor embedding;
- Prop FP (its (6)) and C1;
- "NP-hardness does not exclude H^{o(1)}".

It flagged the narrower oracle as OPEN.

**Codex adversarial (floorhard_2): MAJOR REVISION.** Applied:
- SS: target range and deterministic prime powers; Turing (not parsimonious) floor-sum step;
- FP2: corrected M = 3m formula, uniform sampling, ambiguous-form step conditional;
- TE: c-dependence, e/x, corrected BHP exponent with t ≍ b³, randomized prime-only consequences;
- DOM: slogan withdrawn, strictness, promise terminology;
- T1: shifted intervals;
- T3: rejected as an obstruction; square padding; T3′ added;
- T4: counterexample, equivalence withdrawn;
- §4.3: scope, no 2^{Ω(L)} from Manders–Adleman;
- §5: C2 identity, 3^s counts, no Ω(m) claim;
- §6: TCH §2.1, Balog–Bárány, Harvey–Hittmeir, exhaustiveness wording;
- §7: rewritten along its replacement text;
- its estimate (≈5%, 1–20%) is quoted.

**Opus proof critic: ACCEPT WITH FIXES.** Same main issues (DOM slogan, §0 Coppersmith-regime
sentence, BHP exponent, SS range, FP2 formula). Also applied:
- T1 tightened to c²κ;
- FP1's p/4 formula restricted to p ≡ 1 (mod 4);
- C2's {0, 9} hypothesis;
- §4.1 "need" changed to "our constructions";
- §4.3 restricted to the NP/#P/ETH families;
- re-centring at block starts noted.

**Opus literature critic: OK WITH FIXES.** Applied:
- JKL 2023 cited and the novelty of SS reworded;
- ABCS bound corrected and Balog–Bárány 1991 moved to the disc;
- TCH pinpoint;
- Subset-Sum ETH citations (EKMR, Jansen–Land–Land);
- Hiary 2016, Berndt, Hudson–Williams, Rabin, Goldwasser–Micali, Coppersmith 1997, Howgrave-Graham,
  May, Shanks, Schoof, Helfgott–Thompson venue, Harvey–Hittmeir;
- HKM venue and Li's remark marked UNVERIFIED;
- BIS relabelled folklore.

Not applied: Chamberland–Dilcher (arXiv 2510.04387) is already cited in `analytic_results.md` for the
prime case and is not repeated here.
