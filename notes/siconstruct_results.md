# siconstruct: constructing smooth integers in short intervals near kN (attack on Prop SI)

Agent `siconstruct`, phase 20, 2026-10-21. DRAFT v1 (before reviews).

**Files.**
- Code: `src/siconstruct/` — `core.py` (exact sieve, N-oblivious smooth lists, Hildebrand–Tenenbaum model),
  `constructors.py` (brute force, MITM, 2-D aim), `theory.py` (cost formulas; Prop CV inequalities),
  `run_exp.py` (E1–E3), `diag_e3.py` (post-hoc E3 diagnostic).
- Tests: `tests/test_siconstruct.py` (8 pass).
- Pre-registration: `notes/siconstruct_prereg.md`, commit 9d75362, made before every run.
- Raw data: `experiments/raw/siconstruct_E{1,2,3,3diag}.jsonl`.
- GPT-6.1-sol: blind round `codex/siconstruct_1.{md,out}` (§8).

**Labels.** PROVED (proof here, elementary unless said) · KNOWN (cited) · HEURISTIC (named assumption) ·
MODEL (exact consequence of named heuristics) · SCENARIO (stipulated algorithmic assumptions) · MEASURED ·
POST-HOC · assessment. "No prior statement found by us" is the strongest novelty claim made.

**Notation.** x = ln N, F_α = x^α (ln x)^{1−α}, L[α, c] = exp((c + o(1))F_α). Factor base y = L[1/3, κ], n = π(y);
window h = L[2/3, ν]; I_k = [kN, kN + h]; multipliers k ≤ K = L[1/3, C]; ε := h/(kN), the relative window.
Target: a y-smooth u ∈ I_k at cost L^{o(1)} per output (Prop SI ⇒ L[1/3, 1.2114]).

---

## 0. Verdict

**Plain language.**
- **We found no constructor.** Every approach we analysed costs at least L[2/3] per output (sieving-type methods),
  or exponentially much (meet-in-the-middle, k-list, knapsack algorithms), or far more (lattice CVP in dimension
  L[1/3]). None comes near L^{o(1)}, or even near the weaker L[1/3, 0.757] per output that would still beat MNFS.
- **Meet-in-the-middle does not amortise.** One pair of sorted lists does serve every multiplier k at once, but the
  probability that a product lands in *some* window stays (h+1)/N per pair, whatever the number of k served. So R
  outputs need |A||B| ≈ RN/h pairs and cost ≥ 2√(RN/h). Per output that is ≈ √(N/(Rh)) = N^{1/2−o(1)} for the
  L[1/3] outputs Prop SI needs. The per-output cost becomes L^{o(1)} only if Rh ≥ N^{1−o(1)}, i.e. with exponentially
  many outputs. Measured (E1): over K = 1…256 served multipliers the hit rate per pair stays at 0.92–1.09 of (h+1)/N.
- **Why everything fails: the problem is scale-invariant.** If u = d·j lies in [X, X + h], then d lies in a window
  of the *same relative width* h/X around X/j. Splitting u into parts never relaxes the precision. A factor can be
  aimed cheaply only when its partner is a *tested* (random) cofactor, and then the tested part is ≥ N/(Kh) =
  N^{1−o(1)} (Lemma RC). That is lattice sieving, and it gives L[2/3] per output.
- **The lattice view is exact.** A squarefree Prop SI output is an approximate closest vector with factor
  1 + O(1/ln N) in Adleman's prime lattice of dimension π(y) = L[1/3], and conversely (Prop CV). The instance sits
  at distance ≈ λ₁/√2 from the lattice. So Prop SI is near-exact CVP in dimension L[1/3]; every known CVP algorithm
  at that quality costs 2^{Θ(n)} = exp(L[1/3]).
- **One-sided is not easier than Schnorr's two-sided target.** It is the same lattice problem with a larger residue.
  Prop SI's constant 1.2114 is exactly schnorrcvp's Theorem G constant at the oracle slack τ = 2/(3κ) = 1.10.
- **Choosing k freely buys nothing at L-scale.** K = L[1/3] multipliers give K times more witnesses, which is
  L[2/3, o(1)]. We did not solve, or approach, the one-sided Bohr-set lister.
- **Johnson-radius / Boneh decoding** keeps exponent 2/3: certified outputs have a known divisor d ≥ h, so the
  untested cofactor is ≥ N/(Kh), and the gain is L[1/3, O(1)] on an L[2/3] baseline.

**Main results.**

| Result | Status |
|---|---|
| Prop CV: squarefree Prop SI output ⇔ CVP_{1+O(1/x)} answer (scale C = kN/h); dist ≈ λ₁/√2; slack t costs c(t) | PROVED (inequalities); cost reading MODEL |
| Lemma SR (scale invariance) and Lemma RC (constructed factor + tested cofactor ⇒ cofactor ≥ N/(K(h+1))) | PROVED (SR, RC(i)); RC(ii) HEURISTIC H_eq |
| Lemma MA (two-list amortisation): hits = \|A\|\|B\|(h+1)/N for any K; per output ≥ 2√(N/(R(h+1))) | HEURISTIC H_eq; MEASURED (E1, E2) |
| k-lists / Wagner / representations cannot remove the exponential cost | MODEL (Prop W of smoothenum + representation count) |
| Per-output budget: L[1/3, τ] beats MNFS iff τ < 0.757; τ = κ (linear in n) gives 6^{1/3} = 1.817 | MODEL (exact) |
| Best toy constructors: MITM beats interval sieving at all tested sizes (32–80 bits), model crossover ≈ 310 bits | MEASURED / MODEL |

---

## 1. What is needed, and a softer target

- Prop SI (suffmap §1.4): with a constructor of per-output cost L[1/3, τ] the factoring constant is
  c(τ) = min_κ max(κ + τ + 2/(9κ²), 2κ) (`theory.si_cost`, exact sympy).

| τ | 0 | 0.1 | 0.25 | 0.5 | 0.757 | κ (linear in n = π(y)) |
|---|---|---|---|---|---|---|
| c(τ) | (16/9)^{1/3} = 1.2114 | 1.2819 | 1.4031 | 1.6447 | 1.9019 (= MNFS) | 6^{1/3} = 1.8171 |

- The 0.757 threshold re-optimises κ (at κ = (4/9)^{1/3}). At the Prop SI optimum κ = (2/9)^{1/3} held fixed, the
  threshold is τ < 1.9019 − (16/9)^{1/3} = 0.6905 (codex siconstruct_1).
- So the real target is weaker than L^{o(1)}: anything below L[1/3, 0.757] per output beats MNFS, and even a
  constructor costing O(π(y)) per output would give 1.817 (MODEL, linear algebra at ω = 2).
- None of the approaches below reaches L[1/3, anything] per output.

## 2. Structural results

### 2.1 Prop CV: Prop SI is near-exact CVP in dimension π(y) [PROVED inequalities]

Lattice (Adleman weights, as in schnorrcvp §1): for e ∈ Z^n, b(e) = (√ℓ₁e₁, …, √ℓ_ne_n, C⟨e, ℓ⟩), ℓ_i = ln p_i,
target t_k = (0, …, 0, C ln(kN)), scale **C = kN/h**. Write u = ∏_{e_i>0}p_i^{e_i}, v = ∏_{e_i<0}p_i^{−e_i},
Q = Σe_i²ℓ_i, Δ = ⟨e, ℓ⟩ − ln(kN), so ‖b(e) − t_k‖² = Q + C²Δ².

- **(a) Lower bound.** For every e: ‖b(e) − t_k‖² ≥ ln(kN) − 1/(4C²).
  *Proof.* Q ≥ ln(uv) = ln u + ln v, and ln u = ln(kN) + Δ + ln v. So Q ≥ ln(kN) + Δ + 2 ln v ≥ ln(kN) + Δ, and
  min_Δ (Δ + C²Δ²) = −1/(4C²). ∎
- **(b) Witnesses are near-closest.** A squarefree y-smooth u ∈ I_k has Q = ln u ≤ ln(kN) + h/(kN) and 0 ≤ CΔ ≤ 1.
  So ‖b(e) − t_k‖² ≤ ln(kN) + 1 + h/(kN), and the output is a γ-CVP answer with γ² ≤ 1 + (1 + o(1))/ln(kN).
- **(c) Converse.** Suppose a witness exists (dist² ≤ ln(kN) + 2) and e is a γ-CVP answer with
  γ² ≤ 1 + c/ln(kN). Then ‖b(e) − t_k‖² ≤ ln(kN) + c′ with c′ = c + 2 + o(1), and:
  - 2 ln v ≤ c′ + |Δ|, so v = O(1);
  - |Δ| ≤ √(ln(kN) + c′)/C;
  - w := u − v·kN = v·kN(e^Δ − 1) satisfies |w| ≤ v·h·√(ln(kN) + c′)·(1 + o(1)).
  So u ≡ w (mod N) with u y-smooth (by construction) and |w| ≤ O(√x)·h, at multiplier vk. The factor O(√x) is
  L^{o(1)} and is absorbed in ν. *Proof:* the three displayed inequalities follow from (a)'s chain and Q ≥ ln(uv).
  ∎ (`test_prop_cv_converse_small_v` checks them at 80 digits on random vectors.)
- **(d) Slack.** A CVP oracle with γ² − 1 = 2t·F_{2/3}/x returns v ≤ L[2/3, t] and |w| ≤ h·L[2/3, t]. The residue
  then has log-size (ν + t)F_{2/3}, and the factoring constant is
  c(t) = min_κ max(κ + 2/(9κ²) + t/(3κ), 2κ) (`theory.si_cost_slack`): 1.2114, 1.2481, 1.3937, 1.5702, 1.8982 at
  t = 0, 0.1, 0.5, 1, 2. In dimension terms (ln n = κF_{1/3}): F_{2/3}/x = (ln x/x)^{1/3} = κ ln x/ln n ≈
  3κ ln ln n/ln n. So **any L[1/3] gain needs γ − 1 = O(ln ln n/ln n), and 1.2114 needs γ − 1 = O(1/x) =
  O((ln ln n)²/(ln n)³)** in dimension n = L[1/3, κ].
- **(e) Geometry.** Every nonzero lattice vector has ‖b(e)‖² ≥ 2 ln C − ln ln C − O(1) = (2 − o(1))x.
  *Proof.* If C|Δ| ≤ √(2 ln C), then u ≠ v, |Δ| = |ln(u/v)| ≥ 1/(2M) with M = max(u, v), so M ≥ C/(2√(2 ln C)); and
  Q ≥ ln(uv) ≥ 2 ln M − |Δ|. Otherwise C²Δ² > 2 ln C. ∎ Conversely, two squarefree witnesses u₁ ≠ u₂ in I_k give
  e₁ − e₂ with ‖b(e₁ − e₂)‖² ≤ ln(u₁u₂) + 1 ≤ 2 ln(kN) + 2, so λ₁² = (2 + o(1))x. With (a)–(b): the target lies at distance
  (1 + o(1))λ₁/√2, with L[2/3, ν − 2/(3κ)] witnesses at that distance. The Gaussian heuristic, which predicts
  λ₁ ≈ √(n ln y/(2πe)) = L[1/3, κ/2], is off by an L[1/3] factor (as schnorrcvp §1 warned).
- **Reading.**
  - A Prop SI constructor (squarefree outputs) is a promise-CVP_{1+O(1/x)} solver on this family, and conversely, up
    to an O(√x) factor in h and an O(1) factor in k.
  - Known algorithms for CVP at γ = 1 + o(1) cost 2^{Θ(n)} (exact CVP 2^{n+o(n)}: Aggarwal–Dadush–Stephens-Davidowitz,
    FOCS 2015, arXiv 1504.01995 Cor. 6.5, checked by codex; heuristic sieving 2^{0.29n}). Here that is exp(L[1/3, κ]).
  - BKZ-β with Babai gives approximation factors ≥ 1 + Ω(1) for every β < n (SCENARIO: root-Hermite factor
    δ_β^{Θ(n)}; codex's warning in schnorrcvp that a root-Hermite factor is not a CVP ratio applies). By (d) that
    inflates the residue to N^{Θ(1)}: useless.
  - Even reading the basis costs n = L[1/3, κ] per call, and LLL costs ≥ n² = L[1/3, 2κ] > L[1/3, 0.757].
  - **Smaller dimension does not help.** Restricting to a random 1/r of the primes keeps a witness only if
    r^{u} ≤ L[2/3, ν − 2/(3κ)], u = x/ln y, i.e. r ≤ x^{κ(ν − 2/(3κ))} = poly(x) (HEURISTIC count; uses
    F_{1/3}F_{2/3} = x ln x). So the dimension stays L[1/3, κ].
- **One-sided versus two-sided (quantified).** schnorrcvp's Theorem G: a benign two-sided oracle at slack
  γ² − 1 = 2τF_{2/3}/x costs L[1/3, 2√(τ/3)]. Prop SI's 1.2114 is that formula at τ = 2/(3κ) = 1.1006:
  2√(1.1006/3) = 1.2114. One-sided needs t = 0 to reach the same constant (its residue starts at h = L[2/3, 2/(3κ)]
  instead of B). So in oracle terms one-sided is *harder*, not easier: it is the same oracle with v forced to O(1),
  and the price is an additive 2/(3κ) in the residue exponent. And two-sided exact CVP already gives polynomial time
  in polylog dimension (schnorrcvp Theorem H). The only sense in which one-sided is easier is that its *witnesses*
  are standard smooth-number-in-interval objects, and that the interval sieve solves it at L[2/3] (no CVP algorithm
  does as well).

### 2.2 Lemma SR (scale invariance) and Lemma RC (random cofactor) [PROVED; (ii) HEURISTIC H_eq]

- **Lemma SR.** If u = d·j ∈ [X, X + h], then d ∈ [X/j, (X + h)/j], an interval of relative width h/X containing at
  most ⌊h/j⌋ + 1 integers. ∎ Every factor of u faces the same relative precision ε = h/X. Recursion or splitting
  never lowers ε; it only moves the problem to a smaller scale, where the window has fewer than one integer once the
  partner exceeds h.
- **Lemma RC.** Fix a constructed integer s (smooth, a square, anything) and K. Let
  T_s = {(k, w) : 1 ≤ k ≤ K, 0 ≤ w ≤ h, s | kN + w}.
  - (i) [PROVED] |T_s| ≤ K(⌊h/s⌋ + 1), and every (k, w) ∈ T_s has tested cofactor m = (kN + w)/s ≥ kN/s.
  - (ii) [HEURISTIC H_eq: s chosen independently of N mod s] E|T_s| = K(h + 1)/s·(1 + o(1)). So s yields a pair
    only if s ≲ K(h + 1), and then m ≥ N/(K(h + 1)).
  - (iii) [MODEL] The cofactor has log-size x − νF_{2/3} − O(F_{1/3}) and is y-smooth with probability
    L[2/3, −2/(3κ)]·L[1/3, 2ν/(3κ)]. The gain over a random u ≍ N is L[1/3, 2ν/(3κ)] (derivative of u ln u is
    (2/3 + o(1)) ln x). Per output: **L[2/3, 2/(3κ)]·L[1/3, −2ν/(3κ)]**.
  - (iv) [MODEL] As a factoring method (cofactor and residue both tested, total random size N/K) Theorem SS of
    suffmap gives at best L[1/2, √2].
- **Covers:** continued fractions / Farey / Stern–Brocot approximations of kN/s (a 2-D lattice of determinant s),
  special multipliers k ≡ −wN^{−1} (mod s), the base-m expansion of N, "smooth times arbitrary square" u = s²m
  (s² is the constructed part; QS is the case where the *square* is the whole of u and the residue is √N),
  and the Boneh/Coppersmith certified divisor of §4.4.
- **Measured (E3).** Successful cofactor tests match the H-T density *at the cofactor's size* (pooled 0.95, cells
  0.49–1.25), i.e. the aim gives exactly the smoothness gain of a smaller random cofactor and nothing more.

### 2.3 Lemma MA (amortised two-list merging) [HEURISTIC H_eq; MEASURED]

- *Setting.* Lists A, B of y-smooth integers generated without reference to N (e.g. sorted by ln). One merge serves
  every k: all pairs with a·b ∈ ∪_k [kN, kN + h].
- *H_eq.* For N-oblivious lists, a·b mod N is equidistributed: E[#pairs with ab ≥ N, ab mod N ≤ h] =
  |A||B|(h + 1)/N·(1 + o(1)), **for every product range** [N, KN] (the windows have density (h+1)/N among integers
  in any range of length ≫ N).
- *Cost.* Generating and sorting the lists costs ≥ |A| + |B|. In the Prop SI regime (Kh ≪ √N ≈ min A) each a has at
  most O(1) admissible b-intervals with an integer in them, findable by a 2-D lattice/continued-fraction step, so the
  merge is O(|A| + |B| + hits) independently of K. (At toy sizes h > √N and our implementation pays per (a, k);
  E1 records both counters.)
- *Consequence.* R outputs need |A||B| ≥ RN/(h + 1), so cost ≥ 2√(RN/(h + 1)), i.e. per output
  **≥ 2√(N/(R(h + 1)))**.
  - L^{o(1)} per output iff R(h + 1) ≥ N·L^{−o(1)}.
  - Prop SI uses R = L[1/3, κ + ν/(3κ)] and h = L[2/3, ν], so Rh = N^{o(1)} and the per-output cost is
    N^{1/2 − o(1)}. The total, ≥ √(RN/h), is exponential.
  - The amortisation over k is real (one list, all k) but worthless: what is amortised is the *search*, while the
    *supply* of hits is fixed at (h + 1)/N per pair.
- *Distinct outputs.* A smooth u in a window is hit once per factorisation u = ab with a ∈ A, b ∈ B, so distinct
  outputs ≤ hits (E1: 2453 distinct of 2901 hits; codex's toy run: 264 of 4381). Multiplicity only lowers the
  distinct yield.
- sihard E3 tests the same two-list bound exactly (average over window positions) — see §8 once its note exists.

### 2.4 k lists, Wagner, Schroeppel–Shamir, representations [MODEL]

- Hits among k lists are ≤ ∏|L_i|·(h + 1)/N. Generic k-list algorithms in the low-solution regime (one solution per
  ≈ N/h of the tuple space) find a hit at ≥ (N/h)^{1/2} cost (Schroeppel–Shamir time, BCJ-type exponents are also
  constant fractions of log(N/h)).
- Wagner's k-tree needs lists of size (N/h)^{1/(1+lg k)}, i.e. tuple space (N/h)^{k/(1+lg k)}, which exceeds the
  available ≤ N^{1+o(1)} ordered factorisations for k ≥ 3 (smoothenum Prop W, extended here to products in any range:
  the number of k-tuples of positive integers with product ≤ KN is ≤ KN(ln KN)^{k−1}).
- **Representations** (Howgrave-Graham–Joux, BCJ) exploit many factorisations per solution. The average number of
  ordered k-factorisations of a y-smooth u in the window is k^{u(1+o(1))}, u = x/ln y, i.e. L[2/3, o(1)] for fixed k
  and L[2/3, O(1)] for k ≤ x (ρ_k(u) ≈ ρ(u)k^u). That is N^{o(1)}: it cannot remove an N^{Θ(1)} cost.

## 3. Per-approach cost table

Per-output cost of producing one y-smooth u ∈ I_k (any k ≤ K), y = L[1/3, κ], h = L[2/3, ν], R = L[1/3] outputs.

| # | Approach | Per-output cost (heuristic) | L^{o(1)}? | Total factoring cost via this constructor | Status |
|---|---|---|---|---|---|
| 0 | Interval sieve (baseline) | 1/ρ(x/ln y) = L[2/3, 2/(3κ)] | no | L[1/2, √2] (u and v both tested; Theorem SS, λ ≈ x) | MODEL; MEASURED toy |
| 1 | Log-knapsack / CVP in Adleman's lattice | near-exact CVP_{1+O(1/x)} in dim L[1/3]: 2^{Θ(n)} = exp(L[1/3, κ]); BKZ-β, β < n: residue N^{Θ(1)} | no | worse than trial division | Prop CV PROVED; costs SCENARIO |
| 2 | MITM, sorted logs, one list for all k | 2√(N/(R(h+1))) = N^{1/2−o(1)} | no | ≥ √(RN/h): exponential | Lemma MA; MEASURED E1/E2 |
| 2′ | k-list / k-tree / representations | ≥ (N/h)^{1/2}/R^{O(1)} | no | exponential | MODEL |
| 3 | Free k (joint (k, u) search) | factor K = L[2/3, o(1)] more witnesses; no change in per-pair or per-test rates | no | — | MODEL |
| 4 | CRT / Boneh / Coppersmith / GSS | certified d ≥ h ⇒ Lemma RC: L[2/3, 2/(3κ)]·L[1/3, −Θ(1)] | no | L[1/2, √2] | MODEL (pilist Lemma S≥V; Prop JR does not apply: congruence scaling) |
| 5a | 2-D aim: CF / Farey / Stern–Brocot / special k / base-m / u = s·m | L[2/3, 2/(3κ)]·L[1/3, −2ν/(3κ)] | no | L[1/2, √2] | Lemma RC; MEASURED E3 |
| 5b | Oblivious smooth families: products of consecutive integers, binomials, factorial ratios, norm forms | hit probability (h+1)/N per member ⇒ N/h per output; with MITM as row 2 | no | exponential | HEURISTIC H_eq |
| 5c | Square × smooth (QS-type aim, u = s²m) | aim precision ≈ √N: √N/h per output | no | QS itself is the m = 1, residue √N case | MODEL |
| 5d | Exact interval counting + binary search; global samplers; log-grid DP | counting: O(log h) oracle calls, but no cheap exact count known (needs relative precision h/X); samplers ≈ N/h; DP ≥ N^{1−o(1)} bins | only with a new counting algorithm | — | reduction PROVED; rest MODEL (§4.6) |

## 4. Approaches in detail

### 4.1 Log-knapsack (approach 1)
- u = ∏p_i^{e_i} ∈ I_k ⇔ Σe_i ln p_i ∈ [ln kN, ln kN + ε], ε ≈ h/(kN) = e^{−x + O(F_{2/3})}.
- Exponent vectors: ‖e‖₁ ≤ x/ln 2; a typical witness has ≈ x/ln y = x^{2/3}(ln x)^{−2/3}/κ nonzero entries, mostly 1
  (a constant fraction of witnesses is squarefree, HEURISTIC).
- As a subset-sum instance: n = L[1/3] items, weight ≈ u = x^{2/3−o(1)} (sublinear), x bits of target precision.
  Candidates Ψ(N, y) ≈ N·L[2/3, −2/(3κ)], so the solution count is only L[2/3, ν − 2/(3κ)] = 2^{o(x)}: a
  **density-≈1 instance**, where every generic algorithm (MITM, Schroeppel–Shamir, HGJ/BCJ, Wagner) costs
  2^{Θ(x)} = N^{Θ(1)}. The "high density" n ≫ x is illusory: the bound Σe_i ln p_i ≤ x leaves only N^{1−o(1)}
  candidates.
- Lattice: Prop CV. Exact CVP is needed, in dimension L[1/3].

### 4.2 Meet-in-the-middle (approach 2): the amortisation analysis
- Lemma MA. The per-pair hit probability (h + 1)/N is fixed by the target set, not by how many k one list serves.
  Serving K multipliers enlarges the target set K-fold but also requires products K-fold larger; the density of the
  windows among integers is (h + 1)/N at every scale.
- **Precise amortised statement.** For R outputs the optimum is |A| = |B| = √(RN/(h+1)), cost 2√(RN/(h+1)) + R:
  - per output: 2√(N/(R(h+1))) + 1;
  - L^{o(1)} per output ⇔ R ≥ N/(h·L^{o(1)}), i.e. exponentially many outputs (unless h = N^{1−o(1)}, where the
    residue v is no longer small and the method is Dixon-like, L[1/2, √2]);
  - Prop SI needs only R = L[1/3, κ + ν/(3κ)]: the amortisation saves a factor √R = L[1/3] on N^{1/2}.
- **Factory variant** (lists precomputed once for all N of a size): per-N cost |A| ≥ RN/(h|B|), so per-N L^{o(1)}
  needs a precomputed table of size |B| ≈ N/h = N^{1−o(1)}. Exponential precomputation, comparable to a factor table:
  irrelevant (assessment).
- **Toy scale is misleading in MITM's favour.** On the toy L-path h ≈ N^{0.6–0.7}, so N/h is tiny. MITM beats the
  interval sieve at every tested size (E2), and the o(1)-free model crossover is near 310 bits. Asymptotically the
  per-output slope of MITM tends to 1/2 bit per bit of N, while the sieve's tends to 0.

### 4.3 Free choice of k (approach 3)
- The union of K windows has K·h·L[2/3, −2/(3κ)] smooth elements: K = L[1/3, C] is L[2/3, o(1)], so the existence
  threshold ν > 2/(3κ) is unchanged.
- Per-pair (MITM) and per-test (sieve, 2-D aim) success rates do not depend on K. A multi-target CVP gains at most a
  factor K. So free k is worth L[1/3] at most, on top of costs of size L[2/3] or N^{Θ(1)}.
- **Bohr-set check.** None of our methods lists u with u mod N ≤ L[2α], α < 1/3, or ≤ y. The one-sided Bohr-set
  lister (suffmap E2) at α = 1/3 is just Prop SI with free k, and we do not solve it. The hit law (h + 1)/N per
  oblivious pair shows the cost of every listed object grows as h shrinks, as it must.

### 4.4 CRT / Boneh / Coppersmith (approach 4)
- Boneh (STOC 2000) finds all z ∈ [X, X + H] with gcd(z, P) ≥ S when ln H ≲ (ln S)²/ln P. Certified outputs have a
  known divisor d ≥ H (pilist Lemma S≥V, PROVED there). Here H = h, so d ≥ h, and Lemma RC applies to d: the untested
  cofactor is ≥ N/(Kh) = N^{1−o(1)}.
- **Full smoothness is impossible inside the decoding range.** For d ≈ N one needs ln P ≤ x²/ln h ≈
  x^{4/3}(ln x)^{−1/3}/ν, so P involves only O(x/ln x) primes. The window then contains a T-smooth integer
  (T = those primes) with probability ≈ h·L[2/3, 1/(3κ)]/N: exponentially small.
- **Best partial use** (S ≈ h): about h·D(S)/S ≈ D(S) = L[1/3, 2ν/(3κ)] certified outputs per call (D(S) = number of
  divisors of P near S), each with an untested cofactor ≈ N/h. Per output: L[2/3, 2/(3κ)]·L[1/3, −Θ(1)]. Larger S
  loses L[2/3] in output count for an L[1/3] cofactor gain.
- **Relation to Prop JR.** pilist Prop JR (L^{o(1)} batching gain) is proved in NFS scaling, and explicitly excludes
  this congruence scaling, where one call gains g(d) = L[1/3, Θ(1)] (pilist critic M2). We confirm the L[1/3, Θ(1)]
  gain here, and that it leaves the exponent 2/3 untouched. So Prop SI's regime differs from JR's only in the
  lower-order gain.

### 4.5 Own ideas (approach 5)
- **Continued fractions / Farey / Stern–Brocot.** Convergents p/q of kN/s give s·p − kN·q small, a 2-D lattice of
  determinant s: Lemma RC.
- **Products of consecutive integers.** m(m+1)…(m+r−1) has smoothness probability ρ(u/r)^r = ρ(u)·r^u (better by
  r^u), but every factor is tested and consecutive values are spaced ≈ rN^{1−1/r} ≫ h. For fixed k there is at most
  one candidate, with probability ≈ h/(rN^{1−1/r}): exponential.
- **Binomial coefficients, factorial ratios** C(n, j), n ≤ y: smooth by construction, but an N-oblivious family.
  Products of them are a sub-semigroup of the y-smooth integers, so aiming them is the same knapsack (H_eq: hit
  probability (h + 1)/N per member).
- **Norm forms** (u = N(α), α a product of small primes of Z[i] or Z[√r]): the norm is multiplicative, so again the
  knapsack over split primes. Units do not change the norm.
- **Smooth-twin ladders.** Multiplying by ratios (j + 1)/j (j | u) steps ln u by ≈ 1/j ≥ 1/y. Precision e^{−x}
  needs exponentially many combinations: a knapsack with small items, MITM-bounded.
- **Square × smooth.** u = s²m with s arbitrary is enough for the linear algebra (s² has even exponents), but it is
  Lemma RC with the constructed part s² ≤ Kh, or QS (m = 1, residue √N).
- **NFS reading.** Prop SI is two-polynomial "NFS" with f₁ = X, f₂ = X − N and a skewed box (a ∈ [kN, kN + h],
  b = k), with the f₁ side constructed. It is the degree-1 analogue of smoothenum's algebraic-side lister
  (infimum 0.7631 there).
- **Modulus-oblivious.** A constructor that sees only X = kN is a generic "smooth integer in [X, X + L_X[2/3]]"
  algorithm. It would also speed up prime-field discrete-log relation collection (suffmap §4, row C1; not checked in
  detail). That is further evidence that it is hard (assessment).

### 4.6 Counting, sampling and dynamic programming (codex siconstruct_1, checked by us)
- **Counting ⇒ constructing [PROVED].** An exact oracle for Ψ(b, y) − Ψ(a − 1, y) gives a constructor by binary
  search: O(log h) calls find the least smooth integer in I_k, and the factorisation follows by bisecting on the
  prime cutoff. So a cheap exact interval-count would suffice.
- **No cheap count is known.** Subtracting two global estimates of size Xρ needs relative error below
  h/X = e^{−x+O(F_{2/3})}; saddle-point or Perron-integral estimates give nothing of the kind. Exact counting is
  lattice-point counting in a simplex of dimension π(y) = L[1/3] (Barvinok-type algorithms are polynomial only in
  fixed dimension); Bernstein's exact counting uses a meet-in-the-middle whose cost is governed by the two generated
  sets (codex: Bernstein, "Enumerating and counting smooth integers", §3, Algorithm 3). For an *explicit arbitrary*
  set of primes the existence question is NP-hard by sihard's X3C reduction (see sihard's note), so a fast counter
  would have to use that the factor base is an initial segment of the primes.
- **Sampling.** Uniform or approximate generators of factored smooth numbers (Bach; Kalai; Bach–Sorenson,
  arXiv 2006.07445 Thms 3.1, 5.1–5.3) sample globally. Conditioning on I_k costs ≈ N/h trials, and their
  enumeration order is lexicographic in the prime factors, not by size, so "select the rank near Ψ(kN, y)" does not
  select an integer near kN.
- **Log-grid dynamic programming.** A DP over rounded logarithms with per-item error η and s items needs sη ≲ h/X,
  hence ≥ (ln X)/η = N^{1−o(1)} bins (codex; for that explicit-grid class only).
- **Existence theorems** are far from this scale: the best all-interval results need y ≥ exp(C(ln X)^{2/3}(ln ln X)^{4/3})
  (Jain, arXiv 2502.10530, Thms 1.1–1.2, as reported by codex; not checked by us), while here ln y ≍ (ln X)^{1/3}.

## 5. Total factoring cost and heuristic R

- Every constructor we found has per-output cost ≥ L[2/3]; plugged into Prop SI that gives L[2/3]. Used directly as
  relation generators (u and the residue both tested, total random size N^{1−o(1)}), they are at best Dixon-like,
  L[1/2, √2] (Theorem SS), which is worse than QS.
- MITM: total ≥ √(RN/h), exponential, on top of linear algebra L[1/3, 2κ].
- **Diversity (codex).** Multiplying one relation u₀ ≡ v₀ by a smooth s gives su₀ ≡ sv₀, a new valid output with the
  same exponent-difference vector: worthless. A constructor must produce L[1/3, κ] *independent* relations, and if one
  output per multiplier is used, the multiplier range must have C ≥ κ + ν/(3κ).
- **Fruitfulness (heuristic R).** Not tested here: no constructor reached the regime where it matters. MITM outputs
  are N-oblivious products whose residues are determined by (a·b mod N); sicheck's Prop B (Blum N: fruitfulness is a
  parity condition computable from N) applies to any constructor, and sicheck's `beta0` example shows a constructor
  can be useless while producing valid relations.

## 6. Experiments (pre-registered, commit 9d75362)

Toy L-path (o(1) dropped): y = round(e^{κF_{1/3}}), h = round(e^{νF_{2/3}}), κ = (2/9)^{1/3}, ν = 1.376.
At 32/48/64/80 bits: y = 37/90/184/342, log₂h = 22.8/31.2/38.8/45.9.

### E1. MITM hit law and amortisation over k (48 bits, 3 N, K = 1…256, |A| = |B| = 4096)
- Pooled hits 2901 vs predicted 2888.2: **ratio 1.004** (P1a [0.85, 1.15]: PASS).
- Cells 0.776–1.189 (P1b [0.6, 1.5]: PASS).
- Per-K pooled ratios: K = 1: 0.944, 4: 1.044, 16: 1.085, 64: 0.924, 256: 1.025; max/min 1.174 (P1c ≤ 1.4: PASS).
  Serving 256 times more multipliers with one pair of lists leaves the per-pair hit rate unchanged.
- Distinct u: 2453 of 2901 hits (15% repeats through multiple factorisations).
- Cost counters (seed 0): queries per hit 41.9, 136.5, 498.4, 2492.7, 8621.8 for K = 1…256. At toy size h > √N, so
  each a sees ≈ 3K nonempty b-intervals; in the Prop SI regime (Kh ≪ √N) this K-dependence disappears (§2.3).

### E2. Per-output cost along the toy path (32–80 bits, R₀ = 64)
| bits | log₂h | M = \|A\| = \|B\| | distinct outputs | MITM cost/output | model 5M/R₀ | brute force per output: H-T at N / measured (k = 1..8) |
|---|---|---|---|---|---|---|
| 32 | 22.8 | 157, 166 | 44, 41 | 14, 15 | 12, 13 | 1.5e4 / 3.5e4 |
| 40 | 27.1 | 600, 587 | 52, 49 | 40, 42 | 47, 46 | 1.6e5 / 3.6e5 |
| 48 | 31.2 | 2550, 2312 | 56, 40 | 156, 199 | 199, 181 | 1.3e6 / 3.2e6 |
| 56 | 35.1 | 11101, 9441 | 44, 54 | 866, 600 | 867, 738 | 1.3e7 / 3.6e7 |
| 64 | 38.8 | 44076, 40820 | 52, 45 | 2897, 3100 | 3443, 3189 | 7.0e7 / — |
| 72 | 42.4 | 185928, 207848 | 55, 52 | 11537, 13636 | 14526, 16238 | 4.4e8 / — |
| 80 | 45.9 | 1052431, 879911 | 59, 45 | 61001, 66877 | 82221, 68743 | 3.6e9 / — |

- P2a: pooled hits/pred 822/897 = 0.916, cells 0.734–1.109: PASS.
- P2b: cost per distinct output / model 0.74–1.17 (factor 2 allowed): PASS.
- **P2c: FAIL as registered.** Measured sieve density / H-T density *at N* = 0.36–0.50. Cause: the registered
  comparison evaluated the model at N, but the sieved windows sit at kN, k = 1..8, where the density is up to 4×
  lower at toy sizes. POST-HOC, position-matched (mean of the H-T density at kN, k = 1..8): ratios 0.76–1.08.
  The brute-force column above uses the H-T value at N for the comparison with MITM (whose products lie in [N, 4N]).
- P2d: MITM cheaper than brute force at every size: PASS (80 bits: 6.1e4–6.7e4 vs ≈ 3e9).
- Fitted slopes of log₂(per-output cost) per bit: MITM 0.254 (model 0.261), brute force (H-T) 0.365. At toy sizes
  h grows like N^{0.57+}, so the sieve's cost grows faster per bit. Asymptotically the slopes tend to 1/2 and 0:
  the model crossover is ≈ 310 bits (MITM 94.2 vs sieve 93.0 log₂ per output at 320 bits).

### E3. 2-D aim, constructed s and tested cofactor (40/48/56 bits, K = 16, 64 values of s, ≤ 2048 offsets per (s, k))
- P3b (Lemma RC(i)): every tested cofactor has bit length ≥ that of N/(2S0): PASS (all 12 cells, equality).
- **P3a: FAIL as registered for the small-s cells** (S0 = 2^{⌊lg h⌋/2}): distinct outputs / Σδ_HT(m) = 0.38–0.75.
  PASS for the large-s cells (S0 = h/16): 0.83–1.24.
  - Cause (design flaw): with small s, the 64 values of s scan the same short stretches [kN, kN + 2048s], and one
    smooth u is found through several s (multiplicity 1.6–2.4 at 40–48 bits).
  - POST-HOC diagnostic (`diag_e3.py`, same seeded lists): *successful tests with multiplicity* / Σδ_HT(m) = pooled
    5549/5815 = 0.954, cells 0.49–1.25. So each cofactor is smooth like a random integer of its size.
- P3c (tests per output vs brute force): large-s cells 33–140 tests per output against 6.3e5–4.9e7 for the sieve (H-T model at KN/2). At
  toy sizes s ≈ h/16 leaves a cofactor of only 16–24 bits. That is the toy artefact h = N^{0.6+}. At Prop SI scale
  the cofactor is ≥ N/(Kh) = N^{1−o(1)} (on this path with K = L[1/3, 1]: ≈ 67% of the bits at 1024 bits, 79% at 4096).

**What the experiments show.** The laws behind the negative analysis hold at toy size: (h + 1)/N per oblivious pair
regardless of K (E1, E2), and smoothness of the tested part exactly at its size (E3). They do not and cannot show
asymptotic statements. MITM and 2-D aim are excellent at toy sizes only because toy h is a large power of N.

## 7. Deviations and post-hoc items
- P2c and P3a failed as registered; both causes are design errors in the *comparison*, explained above, with
  post-hoc corrected comparisons labelled POST-HOC.
- `diag_e3.py` is an unregistered diagnostic, run after E3.

## 8. Reviews
- **Codex blind round** (`codex/siconstruct_1.out`, GPT-6.1-sol, xhigh, run before our analysis was written down).
  It reached the same verdict independently: interval sieve L[2/3, 2/(3κ)]; construct-one-factor-and-test has a
  cofactor ≥ N/h (its version of Lemma RC); two explicit sorted lists cost N^{1/2+o(1)} per output for
  L[1/3, O(1)] outputs (its version of Lemma MA, with the exact sorted-list accounting Õ(m + n + Km + H));
  Wagner's tree violates M^r ≤ KN for r ≥ 4; Boneh's and GSS's recovery conditions fail by a huge factor; exact
  CVP costs 2^{m+o(m)} (ADS 2015, Cor. 6.5, checked by codex), m = L[1/3, κ]. Items we adopted:
  - a *binary* embedding e_p = Σ2^j z_{p,j} that enforces e ≥ 0, for which exact CVP solves Prop SI and
    γ² ≤ 1 + 3/(m + 1) suffices (an alternative to Prop CV that handles non-squarefree outputs; dimension
    m = O(π(y) ln x));
  - the fixed-κ threshold τ < 0.6905 (§1);
  - the diversity caveat (§5);
  - the counting-to-search reduction, the precision obstruction for approximate counts, the log-grid DP bound and the
    Jain and Bach–Sorenson references (§4.6).
  Codex labels general non-existence UNVERIFIED, as we do: these are barriers for named method classes.

## 9. Reproduction
```
.venv/bin/python -m pytest -q tests/test_siconstruct.py
.venv/bin/python -m src.siconstruct.run_exp E1   # ~1 min
.venv/bin/python -m src.siconstruct.run_exp E2   # ~10 min, < 1 GB
.venv/bin/python -m src.siconstruct.run_exp E3   # ~5 min
.venv/bin/python -m src.siconstruct.diag_e3      # ~5 min (post-hoc)
```
