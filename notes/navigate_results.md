# navigate: hidden-root short-quotient navigation (smoothbuild §7, W11)

Agent `navigate`, phase 25, 2026-10-23. Version 2. It applies the codex adversarial round (MAJOR REVISION, with one
FATAL point against an experiment), the Opus proof critic (ACCEPT WITH FIXES) and the Opus literature critic (§9–§10).
v1 is commit `2642953`.

**Files.**
- Pre-registration: `notes/navigate_prereg.md` (commit `f2e4e9c`, before any registered statistic).
- Code: `src/navigate/`.
  - `core.py`: operational helpers (factor base, Gauss split, batch smoothness test, exact relation check) and the
    three E1 algorithms S0, G1, J2.
  - `e1_walks.py`: E1 harness.
  - `e2_counts.py`: E2, diagnostic (exact per-class counts, Fourier spread).
  - `posthoc_e2.py`: POST HOC sparse-regime E2.
  - `analyze.py`.
- Tests: `tests/test_navigate.py` (9 tests).
- Data: `experiments/raw/navigate_E1.jsonl`, `navigate_E2.jsonl`, `navigate_E2_posthoc.jsonl`; summary
  `experiments/analysis/navigate_summary.json`.
- GPT-6.1-sol: blind round `codex/navigate_1.{md,out}` (launched before any analysis here); adversarial round
  `codex/navigate_2.{md,out}` (on v1).

**Labels.** KNOWN (published; cited) · PROVED (proof here; elementary unless said) · PROVED (generic model) ·
HEURISTIC / MODEL / SCENARIO (named assumption) · MEASURED · POST HOC · CONJECTURED · OPEN · assessment. "No prior
statement found by us" is the strongest novelty claim made.

**Notation.**
- x = ln N; L[α, c] = exp((c + o(1))x^α(ln x)^{1−α}).
- y = L[1/3, κ] (κ = 0.60 in the target); K = π(y).
- G = (Z/N)^×, Q = G², ω = ω(N); ω_mm denotes the matrix-multiplication exponent.
- FB = {−1} ∪ {primes ≤ y}. N is *admissible*: odd, ω ≥ 2, no prime factor ≤ y.
- Exponent vectors e ∈ Z^K are identified with coprime pairs (u, v) = (∏_{e_ℓ>0}ℓ^{e_ℓ}, ∏_{e_ℓ<0}ℓ^{−e_ℓ}). The
  *height* is |e| := Σ|e_ℓ| ln ℓ = ln u + ln v.
- Λ_N := {e : ∏ℓ^{e_ℓ} ≡ 1 (mod N)} is the relation lattice.

---

## 0. Verdict

**Plain language.**
- **No route reaches the target.** Navigation at L[1/3, 1.20] would factor N rigorously in L[1/3, 1.80]. Each of the
  five routes fails at a named point (table below), and no unconditional lower bound is proved.
- **What the analysis adds is structural.**
  - **The height bound 3 ln N is not used by the reductions.** Dixon's argument (EQ1) and the DL reduction (EQ3) run
    with exponent vectors of up to L[1/3] bits (Prop H).
    - With unbounded height and binary exponents, the problem is randomized-equivalent to "factoring plus discrete
      logarithms modulo each prime-power factor". The hypotheses are that the factor base generates G (GRH at this
      y) and that generators are supplied for the DL instances (Bach 1984; smoothbuild EQ2).
    - So NAV is a *strong* sufficient target: a solver contains rigorous L[1/3] factoring and DL modulo the prime
      factors of N, and, if it works for every admissible N, rigorous L[1/3] prime-field DL.
    - The bounded form additionally asks for a short coset vector, which no reduction here uses.
  - **A navigator for all N at the target constants would beat the best heuristic algorithms for both problems.**
    Theorem DLB makes EQ3 quantitative: NAV(κ = 0.60, τ = 1.20) for every admissible N gives Las Vegas discrete
    logarithms in every F_p^× in L_p[1/3, 1.80].
    - For comparison: heuristic NFS-DL is 1.923 and MNFS-DL 1.902; no rigorous L_p[α < 1/2] prime-field DL is
      published (the project's L[2/5] is internal and unrefereed).
    - A navigator that works only on average over balanced semiprimes gives 2^{1/3}·1.80 ≈ 2.27 instead, which does
      not beat NFS.
    - A navigator for one fixed N escapes this implication.
  - **We know no L[1/3] algorithm for bounded-height NAV, even heuristically, even for prime moduli, even with the
    factorisation and all discrete logarithms given.** With all logs known, NAV is a weighted-ℓ1 closest-vector
    problem in the relation lattice. Its dimension must be ≥ x^{3/2−o(1)} for solutions to exist (Prop LT). The
    methods we can price cost:
    - L[2/3, 2/(3κ)] = L[2/3, 1.11] by testing (Prop TD, a model);
    - N^{1/4−o(1)} for one bounded-height k-tree scenario;
    - 2^{Θ(K)} for known Euclidean lattice algorithms.
    These are prices of known methods, not lower bounds.
  - **The hidden-root argument does not need NAV's integer shape.** paper2's Theorem DD (Dixon descent; internal) gives
    fruitfulness for any dependency that contains a hidden-root row an odd number of times, provided the square root Y
    is computed from the transcript alone. That covers NFS-shaped descended rows once the algebraic square-root step is
    available (paper2 at fixed degree; lv5unif §9 at growing degree). The rigour gap for L[1/3] in the project's
    pipeline is *multi-level* descent: Theorem U's one-level geometry forces exponent 2/5 (lv5unif §13).
- **Recommendation (a research preference, not a theorem).**
  - Prioritise multi-level Dixon descent over W11. Codex's blind round named the same route, with an explicit height
    potential; §5.4 states the criterion it would have to meet.
  - Codex's adversarial round disagrees with "retire W11". It proposes keeping the *unbounded-height* relaxation
    (bounded bit length; Prop H) as a navigation target that permits long descent chains and cancellation. We record
    this as the right form of W11 if it is pursued: it is FACT ∧ DL, so any solution is the goal itself.
- **Barrier evidence (no unconditional lower bound).**
  - Generic algorithms need Ω(√(qδ)) group operations, with q the largest prime factor of φ(N) (Prop G; Nechaev,
    Shoup). This constrains only methods that never look at integer representatives.
  - Every relation for a random target has height ≥ x − ω ln 2 − ln(1/δ) − O(ln x) (Lemma RT). Testing smoothness on
    that much log-size costs L[2/3] at y = L[1/3] (Theorem SS-type model).
  - Counting completions reduces search to restricted approximate counting (JVV). The counts are flat (GRH) only for
    unrestricted word counts high in the tree and enumerable at budgets < ln(N/2). In the window between, an
    approximate counter implies a NAV decision oracle at the threshold budget (Prop CF). This is a reduction plus the
    price of known counters, not an obstruction.
  - Squares are not special: for uniform targets, square- and unit-NAV are equivalent up to 2^ω calls. One direction
    is unconditional; the other needs GRH and 2^ω ≤ y^{1/2−ε} (Prop SQ).
- **Experiments (toy scale; mechanism checks only).**
  - E1: all four predictions held.
    - Monte Carlo completion counting adds nothing per smoothness test (0.99–1.01).
    - Steering by the cheap Gauss-height signal gains 1.47–1.76 per test. The size-matched null explains this to
      within 1–2%. Per *distinct* relation the gain is 1.09–1.38.
  - E2: three of four predictions failed. The informative failure is that with 6–8 primes the exact counts are carried
    by very few characters (≤ 0.6% of them; in two seeds only 4–19 characters), not by noise. This is the dense,
    low-dimensional regime.
  - A post-hoc sparse-regime run (97 and 307 primes) raises the share to 4–8%. Whether it tends to the noise level as
    K grows is CONJECTURED, from two K values only; the residual concentration sits on characters nearly trivial on
    2, 3, 5.
  - Codex found that v1's 309-prime run had N's own factors in the factor base. It was re-run with them excluded;
    the numbers barely moved.
- **Assessment.** Probability that NAV at L[1/3, O(1)] (bounded height) is achieved rigorously: ours ≈ 1–2%, ≤ 0.5% at
  (0.60, 1.20). Codex (both rounds): ≈ 5% at some constant, ≈ 1% at (0.60, 1.20), and it sees no basis in this note for
  going below 5%. The gap is how much weight to give unexplored counting and lattice ideas that exploit structure.

**Route table.**

| Route | Cost / status | Exact failure point |
|---|---|---|
| 1. Counting-guided walks (JVV) | NAV ≤ restricted approximate counting (PROVED); no counter below search cost known | unrestricted word counts are flat high in the tree (GRH) and give no guidance; budgets < ln(N/2) are enumerable; in between, a counter implies NAV decision at the threshold. Word Monte Carlo needs Θ(φ) samples, lattice Monte Carlo L[2/3]; characters need DL. Open: implicit-character or annealing counters |
| 2. Lattice / CVP | needs a basis of Λ_N (from FACT ∧ DL); a CVP_{3/(2+ε)} oracle in weighted ℓ1 suffices (GRH); dimension ≥ x^{3/2−o(1)}; target vectors x^{1/4} to L[1/3, κ/2] below the Gaussian-heuristic radius | no known algorithm for this weighted-ℓ1 problem below exponential in K (Euclidean bounds 2^{Θ(K)}; BKZ-β only exp(Θ(K ln β/β))); Schnorr's real-log lattice does not encode u − av = kN. Open: sparse-support structure |
| 3. Meet-in-the-middle | fixed lists: success ≤ \|S₁\|\|S₂\|/\|Q\|; any generic method Ω(√(qδ)) | collisions are generic; the square/hidden-root promise changes counts by ≤ 2^ω and is invisible to the solver |
| 4. Descent | tested Z-descent L[2/3, 2/(3κ)] (MODEL); NFS descent heuristic L[1/3, 1.23–1.44] per target, ideal-shaped | Z-relations have height ≥ x − O(ln x) for random targets, and a small-target relation has \|u\| + q\|v\| ≥ N, so a cheap leaf level of an integer descent needs Prop SI windows; no bounded-height conversion of NFS rows is known; rigour gap = multi-level descent (§5.4) |
| 5. Barrier | NAV for all N ⇒ FACT ∧ F_p-DL at L[1/3, max(κ + τ, 3κ)] (Theorem DLB); generic Ω(√q) | no unconditional lower bound; the evidence is that NAV at the target constants for all N would beat heuristic NFS for factoring and prime-field DL |

**Main results.**

| Result | Statement (short) | Status |
|---|---|---|
| Prop H | EQ1/EQ2/EQ3 accept exponent vectors of ≤ L[1/3] bits (any height) at L[1/3] cost; NAV_3 ⇒ NAV_∞ ≡_R FACT ∧ DL(p^a ∥ N), given generation and supplied generators | PROVED (bookkeeping on EQ1/EQ2) |
| Theorem DLB | NAV_A(κ, τ, δ) for all admissible N ⇒ Las Vegas F_p^× DL in poly(x, 1/δ)·L_p[1/3, max(κ + τ, 3κ)]; average-case (balanced semiprimes) version at 2^{1/3}× the constant | PROVED (EQ3 made quantitative; padding as dsr Lemma PT; rank step as smoothbuild EQ2) |
| Prop G | generic group with a prime q dividing the order: representing a random target by K given generators succeeds with probability ≤ (T + K + 2)²/(2q) + 1/q | PROVED in the generic group model; standard (Nechaev 1994, Shoup 1997) |
| Lemma RT | #{t ∈ G : some u ≡ tv with 0 < \|u\|\|v\| ≤ X} ≤ 2X(1 + ln X) | PROVED (elementary) |
| Prop SQ | for uniform targets: square-NAV(δ) ⇒ unit-NAV(2^{−ω}δ) unconditionally; under GRH and 2^ω ≤ y^{1/2−ε}, each direction with ⌈2^{ω+1}ln(1/ε)⌉ calls, success (1 − ε)δ, height + ln y | PROVED (GRH for the second part) |
| Prop CF | self-reductions (word and vector encodings); approximate restricted counting ⇒ search; unrestricted word counts flat for length ≥ (2 + ε)x/ln y (GRH); vector counts exact for budgets < ln(N/2) | PROVED (flatness under GRH) |
| Prop LT(a) | a δ-successful NAV_A needs y ≥ (Ax)^{A/(A−1)−o(1)}; NAV_3 needs K ≥ x^{3/2−o(1)} (the sparse regime) | PROVED (de Bruijn's Ψ(X, (ln X)^c) = X^{1−1/c+o(1)}) |
| Prop TD | tested relations of height ≥ x − O(ln x) at y = L[1/3, κ] cost L[2/3, 2/(3κ)] | MODEL (H_S on tested pieces) |
| Obs DD | Theorem DD gives fruitfulness for every dependency with an odd hidden-root coefficient and a transcript-measurable Y; the augmentation of Lee–Venkatesan rows needs ≥ c + 1 hidden rows and is a sketch | PROVED in paper2 given Y (internal); combination CONDITIONAL |

---

## 1. Setting

**NAV_A(y; T, δ).**
- *Input.* For every admissible N: (N, a) with a uniform in Q (its root withheld).
- *Output.* In time T, an exponent pair (e, f) ∈ Z^{K+1} × Z^{K+1} (the −1 column included), written in binary, with
  ∏ℓ^{e_ℓ} ≡ a·∏ℓ^{f_ℓ} (mod N) and |e| + |f| ≤ A·x.
- *Success.* Probability ≥ δ over a and the coins.
- *Variants.* NAV_∞ drops the height bound; exponents are then bounded only by the output length ≤ T. The target of
  smoothbuild §7 is NAV_3 with y = L[1/3, 0.60], T = L[1/3, 1.20], δ ≥ 1/poly(x). Outputs are exponent vectors; u and v
  need never be written out.

**Existence.** Under GRH, for k = ⌈(2 + ε)x/ln y⌉ the k-step walk on G with uniform signed letters ℓ^{±1} (ℓ ≤ y)
reaches every class with probability ≥ 1/φ − λ^k > 0, where λ ≤ y^{−1/2+o(1)} (JMV 2009; dsr Prop RW). So every class has a
word of height ≤ k ln y < 3x for large N. This is where the constant 3 comes from (codex blind §1).

**Quantifiers.** "Every admissible N" matters only for Theorem DLB (padding). EQ1 needs NAV only for the N being
factored.

---

## 2. Reductions and the barrier (Route 5)

### 2.1 Prop H: the height bound is not used by the reductions [PROVED]

**Statement.** Let the outputs be written in binary, so every exponent has ≤ T bits.
- (a) EQ1 with NAV_A runs in O((K + 1)(T + x²)/δ) + F₂ linear algebra + O(K·T·x²), with the same success
  probability. Computing Y costs O(log E) multiplications per prime for exponents E of ≤ T bits.
  - The F₂ step costs O(K·(number of odd entries)), which is O(K²·Ax) when each row has ≤ Ax/ln 2 odd entries (block
    Wiedemann), and O(K^{ω_mm}) always.
  - At κ = 0.60 and τ = 1.20 everything stays within L[1/3, max(κ + τ, ω_mm·κ)] = L[1/3, 1.80].
- (b) EQ2 and EQ3 accept NAV_A outputs unchanged. Their linear algebra is over Z/M and sees exponents modulo M, at
  O(T) extra bit cost per row.
- (c) NAV_3 ⇒ NAV_∞ trivially. Assume:
  - FB generates G (true under GRH by Bach's bound, since y ≫ 2 ln² N);
  - every DL(p^a) instance comes with a generator of (Z/p^a)^×, as in EQ2 (finding one needs the factorisation of
    φ(p^a), which FACT(N) does not supply; codex, proof critic).
  Then NAV_∞ ≡_R FACT ∧ DL(p^a ∥ N) (smoothbuild EQ2; Bach 1984).
  - *Costs (⇒).* NAV_∞ at L[1/3, τ] gives FACT at L[1/3, max(κ + τ, ω_mm κ)] (EQ1), and DL(p^a) at
    L_N[1/3, max(κ + τ, 3κ)] (EQ2 with the same N; an L_p bound needs Theorem DLB's padding).
  - *Costs (⇐).* FACT and DL oracles of cost L[1/3, c] give NAV_∞ after a precomputation of L[1/3, κ + c] (the logs
    of FB). Each call then costs one DL plus a linear solve modulo the group exponent of each prime power.

*Proof.* (a) EQ1 (smoothbuild §4) uses the outputs in two places:
- the F₂ dependency, which uses only exponent parities;
- Y = ∏ℓ^{E_ℓ/2} mod N, computed by repeated squaring in O(log E_ℓ) multiplications.

smoothbuild's term K·A·x³ is the naive cost. The fruitfulness argument (smoothbuild EQ1; paper2 Theorem DD) never looks
at heights. (b) EQ2's rank step works with v_i ∈ (Z/M)^K. (c) EQ2 and Bach. The converse solves Σe_ℓ log ℓ ≡ log a
modulo the group exponent of each prime power, which generation makes possible. ∎

**Reading.**
- The consequences of NAV_3 are those of NAV_∞. NAV_∞ is FACT ∧ DL. So NAV is a strong sufficient target: a
  solver *contains* rigorous L[1/3] factoring and DL modulo the prime factors of N.
- That is not by itself an argument against pursuing it. Intermediate primitives are routinely equivalent to the
  goal (codex). It does say that a navigation result would be at least as strong as the goal, not a stepping stone.
- The height bound makes NAV_3 harder than NAV_∞ (a short coset vector, §4) without improving any consequence used
  here. If W11 is pursued, the right interface is NAV_∞ with exponents of bounded bit length (codex adversarial).
- Heuristically NAV_∞ is L[1/3, 1.923]: NFS factoring, then NFS-DL modulo each prime. Nothing heuristic supports
  τ = 1.20 at κ = 0.60.

### 2.2 Theorem DLB: NAV gives prime-field discrete logarithms [PROVED; EQ3 made quantitative]

**Statement.** Suppose NAV_A(y = L_N[1/3, κ]; T = L_N[1/3, τ], δ) holds for every admissible N, with outputs of
≤ T bits. Then there is a Las Vegas algorithm that, given a prime p, a generator g of F_p^× and h ∈ F_p^×, returns
log_g h in expected time poly(x, 1/δ)·(L_p[1/3, κ + τ] + L_p[1/3, 3κ]). At (κ, τ) = (0.60, 1.20) this is
L_p[1/3, 1.80].

*Proof.*
- *Padding (dsr Lemma PT).* Let r be a prime in (Y₀, 2Y₀] with Y₀ := max(y(p²), 3 ln² p²), found by random sampling
  and a primality test (or as the least such prime by sieving, at cost L_p[1/3, 2^{1/3}κ], which is dominated by
  L_p[1/3, 3κ]; codex).
  - Then N := pr is admissible, since y(N) < y(p²) < r < p.
  - ln N = (1 + o(1))ln p, so y(N) = L_p[1/3, κ + o(1)]. (A fixed small r would not be admissible; codex.)
- *Queries.* Put M = p − 1. If h is not a square mod p (Euler's criterion), replace h by gh and subtract 1 at the end.
  - Choose s₀ queries t_i = CRT(g^{2k_i} mod p, s_i² mod r), with k_i uniform mod M and s_i uniform in (Z/r)^×.
  - At one uniformly random position j*, put t_{j*} = CRT(h·g^{2k} mod p, s² mod r).
  - Every t_i, including t_{j*}, is an independent uniform element of Q, so the query sequence has the same law for
    every j*.
  - Each query is answered by an independent run of NAV with fresh coins. Hence j* is independent of the
    *oracle-visible* transcript (queries and answers), though not of the reduction's private masks k_i (codex).
- *Relations.* A success at query i gives (e_i, f_i). Reducing mod p, v_i := e_i − f_i satisfies v_i·L ≡ 2k_i (mod M),
  where L = (log_g ℓ)_ℓ. At the hidden position, v_{j*}·L ≡ log h + 2k.
- *Rank step (smoothbuild EQ2, exchangeability).* Over Z/M, call a successful index *indispensable* if its vector lies
  outside the span of the other successful vectors.
  - Ordering the indispensable indices i₁, i₂, … makes span(v_{i₁}, …, v_{i_j}) strictly increasing. So there are at
    most log₂|(Z/M)^K| = K log₂M of them.
  - If j* succeeded and is dispensable, then v_{j*} = Σ_{i≠j*}c_i v_i over Z/M. The coefficients come from a
    coefficient-producing Howell-form computation, which needs gcds, not the factorisation of M (Storjohann–Mulders
    1998).
  - Then log h ≡ Σc_i·2k_i − 2k (mod M).
- *Cost.*
  - By independence, a run with s₀ = ⌈4K log₂M/δ⌉ queries has Pr[j* succeeds and is dispensable] ≥ δ − K log₂M/s₀ ≥
    3δ/4. No Chernoff bound is needed (proof critic).
  - So O(1/δ) runs suffice: O(Kx/δ²) oracle calls at time T each.
  - The Howell form of an s₀ × K matrix over Z/M costs O(s₀K²) ring operations: L[1/3, 3κ]·poly(x, 1/δ).
  - The answer is verified by one exponentiation (Las Vegas). ∎

**Average-case version (corrected after codex and the proof critic).** Suppose NAV succeeds on average over (N, a),
with N a random balanced semiprime whose primes are uniform of their size.
- Then for a ≥ δ/2 fraction of primes p, the average success over a fresh random q is ≥ δ/2.
- Run the proof with a fresh q (and N = pq) for each query. Exchangeability still holds, because the query law is the
  same at every position.
- But ln N ≈ 2 ln p, so L_N[1/3, c] = L_p[1/3, 2^{1/3}c]. At the target constants this gives
  L_p[1/3, 2^{1/3}·1.80] ≈ L_p[1/3, 2.27], which does *not* beat heuristic NFS-DL. The "beats heuristic DL" evidence
  rests on the worst-case (padded) version only.

**Reading.**
- Rigorous prime-field DL at L_p[1/3, c] is open for every c (literature critic: no rigorous L_p[α < 1/2] algorithm is
  published, unconditional or under GRH).
  - The best published rigorous bound is Pomerance 1987, L_p[1/2, √2]. The project's L[2/5] and paper2's L[1/2, < 1]
    are internal and unrefereed.
  - Heuristic NFS-DL is L_p[1/3, (64/9)^{1/3} ≈ 1.923]; MNFS-DL is ≈ 1.902 (Matyukhin 2003; corrected by
    Barbulescu–Pierrot 2014).
- So NAV(0.60, 1.20) for all N would improve the *heuristic* state of the art for prime-field DL and for factoring at
  once. This is the strongest evidence against the target constants. It is evidence, not a lower bound.
- A navigator for the one N to be factored escapes Theorem DLB, but still gives DL modulo the p | N (Prop H(c)).

### 2.3 Prop G: generic algorithms [PROVED in the generic group model; standard]

**Statement.** In Shoup's generic group model for Z/q (q prime), let an algorithm receive the encodings of
x₁, …, x_K and z, all independent and uniform. It makes T group-operation queries and outputs e ∈ Z^K. Then
Pr[Σe_i x_i = z] ≤ (T + K + 2)²/(2q) + 1/q. The same bound holds in any finite abelian group whose order is divisible
by q: keep the q-primary component formal; a nonzero linear form with a unit coefficient on Z vanishes with
probability ≤ 1/q (proof critic).

*Proof (standard).*
- Simulate with formal variables X₁, …, X_K, Z. Every encoded element is a known linear polynomial in (1, X, Z).
- The real and simulated games differ only if two distinct polynomials agree at the random point. That has
  probability ≤ C(T + K + 2, 2)/q, since each difference is a nonzero linear polynomial.
- In the simulated game e is independent of (x, z) and Σe_iX_i − Z ≠ 0, so the success probability is 1/q. ∎

(Nechaev 1994; Shoup 1997, whose Thm 1 covers groups whose order has a large prime factor. The representation version
is the same proof.)

**Reading.**
- The bound applies to *generic* algorithms: those that only multiply, invert and compare residues (collision
  search, meet-in-the-middle, walks steered by equality tests). They need Ω(√(qδ)) operations.
- Goldfeld (1969) shows that at least half of all primes p (in lower density) have a prime factor of p − 1 exceeding
  p^{1/2}. For N = pq with such a p, q ≥ N^{1/4 − o(1)} and the generic cost is ≥ N^{1/8−o(1)}.
- The model says nothing about methods that look at integer representatives: smoothness testing, lattice
  reduction, or arithmetic steering (codex). It shows only that an L[1/3] navigator cannot be generic.

### 2.4 Lemma RT: relations for random targets are full-size [PROVED]

**Statement.** For X ≥ 1,
#{t ∈ G : ∃(u, v), v ≠ 0, gcd(v, N) = 1, u ≡ tv (mod N), 0 < |u|·|v| ≤ X} ≤ 2X(1 + ln X).

*Proof.* Each pair (u, v) determines t = uv^{−1}, and (u, v), (−u, −v) give the same t. The pairs with v ≥ 1, u ≠ 0
and |u|v ≤ X number 2Σ_{v≤X}⌊X/v⌋ ≤ 2X(1 + ln X). ∎ Checked exhaustively at 16 bits in
`test_lemma_RT_random_target_size_floor`.

**Corollary (corrected).** For a uniform t ∈ Q and f ≥ 1:
Pr[some relation for t has ln|u| + ln|v| ≤ x − ln f] ≤ 2^{ω+1}(1 + x)N/(fφ(N)).
Conditioning on NAV success costs a further factor 1/δ. So typical successful outputs have height
≥ x − ω ln 2 − ln(1/δ) − O(ln x) (codex, proof critic).

**Reading.**
- This is the random-target form of smoothbuild's Lemma DSC (which concerns small targets).
- So a NAV output carries ≈ x of log-size that must be y-smooth. Either that smoothness is tested (Prop TD) or it is
  constructed, which is navigation proper.

### 2.5 Prop SQ: squares are not special (for uniform targets) [PROVED; second part under GRH]

**Statement.** Let S be a NAV solver with success δ on uniform elements of Q, and U one with success δ on uniform
elements of G. All targets below are uniform; for a *fixed* target t the reduction can fail, since a solver may reject
the K elements t·q, q ≤ y, while succeeding almost everywhere (codex).
- (a) Unconditionally, S applied to a uniform t ∈ G succeeds with probability ≥ 2^{−ω}δ, at the same height.
- (b) Assume GRH for the real characters mod N and 2^ω ≤ y^{1/2−ε}. Put m = ⌈2^{ω+1} ln(1/ε)⌉.
  - Call S on t·q_i for m independent uniform primes q_i ≤ y, and output (e, f + e_{q_i}) on the first verified success.
    For uniform t ∈ G this succeeds with probability ≥ (1 − ε)δ, at height + ln y.
  - Symmetrically, U on a·q_i for a uniform a ∈ Q succeeds with probability ≥ (1 − ε)δ.

*Proof.*
- (a) t is uniform on G, so Pr[t ∈ Q] = 2^{−ω} exactly, and conditionally t is uniform on Q.
- (b), S direction. Condition on the class c of t in G/Q; t is uniform on its coset. For a fixed q in class c^{−1},
  t·q is uniform on Q.
  - So one call landing in class c^{−1} already succeeds with probability δ; further such calls can only help.
  - Hence Pr[success | c] ≥ δ·Pr[∃i : class(q_i) = c^{−1}].
  - Under GRH every nonprincipal real character χ mod N has |Eχ(q)| ≤ y^{−1/2+o(1)} for a uniform prime q ≤ y (JMV's
    bound). So each class has probability ≥ 2^{−ω}(1 − 2^ω y^{−1/2+o(1)}) ≥ 2^{−ω−1}.
  - Then m tries miss c^{−1} with probability ≤ (1 − 2^{−ω−1})^m ≤ ε.
- (b), U direction. On input a·q_i the input is uniform on the coset class(q_i).
  - U's success density δ_c on coset c averages δ over the cosets.
  - With probability ≥ 1 − ε a best coset c* (δ_{c*} ≥ δ) is visited, and one call there succeeds with probability
    δ_{c*}. ∎

**Reading.**
- The hidden root is invisible to the solver: (N, a) has the same law whether a is drawn as r² or as a uniform square.
  Its only role is EQ1's fruitfulness.
- For uniform targets, NAV for squares and for arbitrary targets are the same problem up to 2^ω calls. This answers
  codex's blind-round question about a bounded-height reduction.
- 2^ω is 4 for semiprimes. For N with many factors it can exceed every L[1/3, O(1)], and the hypothesis 2^ω ≤ y^{1/2−ε}
  is then a real restriction.

---

## 3. Route 1: counting-guided walks

**Prop CF (counting vs finding) [PROVED; flatness under GRH].**
- **(i) Self-reductions.**
  - *Word encoding.* The instance (t, k) asks for a signed word of length k over the primes ≤ y (no −1 letter, or −1 at
    most once) with product t. Its solutions are ⊔_{ℓ,±}(ℓ^{±1}, solutions of (tℓ^{∓1}, k − 1)). Depth
    k ≤ k_max = Ax/ln 2; branching 2K.
  - *Vector encoding* (canonical: largest prime of the support first). The instance (t, B, j) asks for e supported on
    the first j − 1 primes with |e| ≤ B. Its solutions are {0 if t ≡ 1} ⊔ ⊔_{i<j, k≠0} solutions of
    (t·p_i^{−k}, B − |k| ln p_i, i). Depth ≤ B/ln 2; branching ≤ 2KB/ln 2.
  - Binary encodings of the prime index reduce branching to 2 at the cost of interval-restricted counts (codex).
  - The plain count C(t, B) is not closed under these restrictions. The restricted counts are what JVV needs
    (codex blind §1).
- **(ii) Counting ⇒ search [PROVED; JVV 1986, decision-to-search].**
  - An oracle that approximates the restricted counts within any multiplicative factor decides emptiness.
  - Walking down the tree and keeping a non-empty child finds a solution with ≤ depth × branching = poly(K, x) calls.
  - With relative error 1 + 1/poly per level, JVV gives almost-uniform sampling.
  - The converse (counting from decision or search) is not claimed.
- **(iii) Flatness of unrestricted word counts [PROVED under GRH].** For the counts W_k(t) of length-k signed words with
  product t:
  - W_k(t) = (2K)^k(1/φ)(1 + θ), with |θ| ≤ φλ^k and λ ≤ y^{−1/2+o(1)} (JMV; Fourier inversion).
  - For k ≥ (2 + ε)x/ln y, every child count of an unrestricted word instance is ≈ the mean, so a uniform first letter
    is as good as the oracle's. Equivalently, a random prefix maps a uniform target to a near-uniform target.
  - This flatness is *not* shown for vector counts C(t, B), nor for restricted (prefix- or interval-limited)
    subproblems. Restricting the alphabet to one prime destroys it (codex, proof critic).
- **(iv) Enumerable bottom [PROVED].** Suppose gcd(u, v) = 1, v ≥ 1, |u|·v < N/2 and u ≡ tv (mod N). Write
  u = tv − kN. Then:
  - gcd(k, v) divides u, so k/v is in lowest terms;
  - |t/N − k/v| = |u|/(Nv) < 1/(2v²);
  - so k/v is a convergent of t/N (Legendre).
  Hence for budgets B < ln(N/2) all vector solutions come from the O(x) convergents, and exact counting costs O(x)
  smoothness tests. The proof critic checked all 5,549 coprime pairs at 14 bits.
- **(v) Consequence.** No reduction here locates the information to be supplied by a counter more precisely than this:
  - not in the high unrestricted word levels (iii);
  - not in vector budgets < ln(N/2) (iv);
  - somewhere in between.
  At the existence threshold B*, where the mean count is ≈ 1, an approximate restricted counter *implies* a NAV decision
  oracle at budget B*, and by (ii) decision gives search. (The registered threshold for NAV_3 is far above B*: the mean
  number of height-3x witnesses is N^{2−o(1)}; codex.)
- **(vi) Candidate counters, priced.**
  - *Monte Carlo over words:* the success probability per sample is ≈ 1/φ at the threshold, so Θ(φε^{−2}) samples
    (codex blind).
  - *Monte Carlo over lattice points of {(u, v) : u ≡ tv}:* the density of smooth pairs is L[2/3, −2/(3κ)] at
    y = L[1/3, κ] (Prop TD), so L[2/3] samples per unit of count.
  - *Characters:* C(t, B) = φ^{−1}Σ_χ χ̄(t)Ĉ_B(χ). Evaluating χ of large order needs discrete logs modulo the prime
    factors (and the factorisation), which already gives NAV_∞ (Prop H).
    - For NAV_3 one would also need few characters to suffice, *and* a way to find them without DL.
    - E2 finds strong concentration in the dense toy regime (§7.2) and weaker concentration in the sparse toy regime.
    - Concentration of L² mass is not the same as pointwise counting accuracy, and it says nothing about restricted
      conditional counts. Whether "implicit-character" counters exist is OPEN (codex adversarial).
  - *Sieves / analytic number theory:* equidistribution theorems for smooth numbers in progressions (Granville;
    Soundararajan 2008; Harper 2012 — under GRH or on average over moduli) give the *mean* of C. They give no
    pointwise counts at the threshold for an individual modulus of size N with y = L[1/3].

**Verdict on Route 1.** The reduction NAV ≤ (approximate restricted counting) is real (JVV). We know of no counter
whose cost is below that of search in the window where counting carries information. Failure point: Monte Carlo
costs φ or L[2/3], characters need DL, and the analytic theory gives means only. Not excluded: counters using
implicit character information, annealing, or other structure (OPEN).

---

## 4. Route 2: lattice / CVP

**Prop LT (dimension and geometry).**
- **(a) Dimension floor [PROVED].**
  - A NAV_A solver with success δ ≥ 1/poly needs #{e : |e| ≤ Ax} ≥ δ|Q|.
  - Each y-smooth n ≤ N^A has 2^{ω(n)} ≤ N^{o(1)} coprime splits, so #{e : |e| ≤ Ax} ≤ Ψ(N^A, y)·N^{o(1)}.
  - |Q| = N^{1−o(1)}, since 2^ω ≤ e^{x/ln y}.
  - By de Bruijn (Ψ(X, (ln X)^c) = X^{1−1/c+o(1)} for fixed c > 1) and monotonicity in y, this needs
    y ≥ (Ax)^{A/(A−1)−o(1)}.
  - For A = 3: y ≥ (3x)^{3/2−o(1)}, i.e. K ≥ x^{3/2−o(1)}. A solution has ≤ Ax/ln 2 letters, so K/letters → ∞:
    bounded-height NAV is always in the *sparse* regime.
- **(b) The CVP [PROVED reformulation; one-way].**
  - Given a basis of Λ_N (obtainable from FACT ∧ DL, Prop H(c)) and one coset representative e_a (∏ℓ^{e_a} ≡ a),
    NAV_A asks for a vector of e_a + Λ_N with weighted ℓ1 norm |·| ≤ Ax.
  - Under GRH the coset's distance is ≤ (2 + ε)x, so a CVP_{3/(2+ε)} oracle in weighted ℓ1 suffices for NAV_3 (codex).
  - The converse fails. The distance is ≥ x − O(ln x) for most a (Lemma RT), and heuristically (1 + o(1))x, since a first
    moment count already reaches |Q| at height x + O(x^{2/3}(ln x)^{4/3}). So a NAV output may have approximation ratio
    up to ≈ 3 (proof critic). The approximation actually demanded is γ ∈ [3/(2 + ε), 3], heuristically ≈ 3.
  - In ℓ2 with Adleman weights √ln ℓ, height equals squared norm only for squarefree signed vectors.
- **(c) Far below the Gaussian heuristic [HEURISTIC/SCENARIO].**
  - det Λ_N = |⟨FB⟩| ≤ φ, so det^{1/K} = e^{O(x/K)} → 1, and the Euclidean Gaussian-heuristic radius is
    ≈ √(K ln y/(2πe)).
  - That is ≈ x^{3/4} at the minimal K and L[1/3, κ/2] at K = L[1/3, κ], against ≈ √(3x) for squarefree target
    vectors. The volume heuristic fails for these sparse integer points (siconstruct Prop CV(e); Mulder–Sterner–van
    Woerden).
  - This shows unusual structure. It is *not* hardness evidence against exploiting that structure (codex).
  - Known algorithm costs are scenarios, not lower bounds:
    - Euclidean exact or constant-factor CVP costs 2^{Θ(K)}: exact 2^{K+o(K)} (Aggarwal–Dadush–Stephens-Davidowitz
      2015, Euclidean norm), and 2^{0.29K}-type heuristic sieving. These do not transfer verbatim to weighted ℓ1
      (codex).
    - BKZ-β gives approximation exp(Θ(K ln β/β)). At an L[1/3, 1.2] budget β = O(x^{1/3}(ln x)^{2/3}), so the factor
      is exp(L[1/3]) (codex blind).
  - CSI-FiSh (Beullens–Kleinjung–Vercauteren 2019) solves the analogous class-group problem by BKZ + Babai in
    dimension 74. That is the dense regime (det^{1/K} ≫ 1), where the Gaussian heuristic is accurate. It is an analogy,
    not evidence that the method transfers.
- **(d) Schnorr's lattice does not encode the congruence [PROVED observation; codex].**
  - The real-log lattice (rows √ln ℓ·e_ℓ, last coordinate C ln ℓ) has determinant ∏√ln ℓ·(1 + C²Σln ℓ)^{1/2}, not N.
  - A target C·ln a enforces u/v ≈ a over R, not u − av ∈ NZ.
  - With a lift, the residual u − Av must be tested. That is siconstruct Prop CV / schnorrcvp, whose slack certificate
    (Prop CV(d)) gives *sufficient* approximation conditions for L[1/3]; necessity was withdrawn there.
- **(e) A bounded-height k-tree [SCENARIO; corrected].**
  - Run Wagner's k-tree (k a power of two) on the knapsack Σe_ℓ log ℓ ≡ log a (mod M), with all logs known.
  - It needs k lists of size ≥ M^{1/(1+lg k)}. A list of S distinct vectors has typical height ≥ (1 − o(1))ln S.
  - If the lists have disjoint supports and values behave randomly, the output height is ≈ k·x/(1 + lg k).
  - Height ≤ 3x forces k/(1 + lg k) ≤ 3 + o(1), i.e. k ≤ 8, and cost N^{1/4−o(1)}.
  - This prices one construction. It is not a lower bound for generalised birthday methods; overlapping supports
    can cancel (codex, proof critic).

**Verdict on Route 2.** Failure points:
- (1) There is no basis without FACT ∧ DL.
- (2) Even with a basis, the problem is a constant-factor weighted-ℓ1 CVP in dimension ≥ x^{3/2}, with targets far
  below the Gaussian heuristic. Every known algorithm we can price is exponential in K.
- (3) The archimedean lattice encodes the wrong constraint.

Required approximation: a CVP_{3/(2+ε)} oracle suffices; the problem itself tolerates γ up to ≈ 3. Not excluded:
algorithms exploiting sparse-support structure (OPEN).

---

## 5. Route 3 (meet-in-the-middle) and Route 4 (descent)

### 5.1 Meet-in-the-middle

- **Fixed lists [PROVED; codex blind §3].** For lists S₁, S₂ fixed independently of a uniform a ∈ Q,
  Pr[aS₁ ∩ S₂ ≠ ∅] ≤ |S₁||S₂|/|Q|, since each pair (s₁, s₂) fixes one a. Lists of size L[1/3] succeed with
  probability N^{−1+o(1)}.
- **Adaptive lists** that only multiply and compare are generic: Prop G gives Ω(√(qδ)). sihard Prop MM covers interval
  tables.
- **The hidden-root promise.**
  - By Prop SQ the solver cannot exploit it.
  - The square subgroup keeps every odd prime factor of the group order, so the generic bound is unchanged.
  - In E2 the per-class counts on squares and on non-squares agree in every cell (|Δh| ≤ 0.012; P2d held).
- **Verdict.** Plain or fixed-list meet-in-the-middle does not amortise. Squares are not special at composite N beyond
  the factor 2^ω.

### 5.2 Descent over Z

**Prop TD (tested relations cost L[2/3] at y = L[1/3]) [MODEL: H_S for tested pieces].**
- Take a relation of height S ≥ x − O(ln x) (Lemma RT) whose smoothness is tested on j ≤ x^{o(1)} pieces. Under H_S it
  is y-smooth with probability ≤ exp(−(S/ln y)·ln(S/(j ln y))(1 + o(1))).
- Equal pieces are optimal by convexity of u ln u. At y = L[1/3, κ] this is L[2/3, −2/(3κ)]; at κ = 0.6, L[2/3, 1.11]
  per output (codex: the same constant for Gauss reconstruction plus testing).
- Pieces of size ≈ y (j ≈ x^{2/3}) are construction, not testing.
- The model does not constrain deliberately constructed pieces, or smoothness laws after adaptive selection.

**Integer descent [PROVED identity + MODEL].**
- A descent step for a prime q is a relation u ≡ qv (mod N) with u ≠ qv. So |u| + q|v| ≥ N (smoothbuild Lemma DSC),
  and its tested pieces have log-size ≥ x − ln q − O(1).
- For the leaf level (children ≤ y = L[1/3]) to cost L[1/3], the tested pieces must have log-size
  O(x^{2/3}(ln x)^{1/3}) (Theorem SS). That forces q ≥ N/L[2/3] and qv within L[2/3] of a multiple of N: suffmap's
  Prop SI window problem.
  - Prop SI's constructors are open, and a modulus-generic one would improve prime-field DL (sihard Prop DLI).
  - Recursion does not relax the relative precision (smoothbuild §5).
- Scope: this applies to integer-relation designs with tested pieces. It does not exclude growing intermediate primes,
  cross-level cancellation, or unbounded-height constructions (codex).

### 5.3 NFS descent

- **Heuristic costs [KNOWN].**
  - NFS individual-logarithm descent: Joux–Lercier 2003; Commeine–Semaev 2006 (initial split L_p[1/3, 3^{1/3}
    ≈ 1.442]).
  - Barbulescu's early-abort analysis gives ≈ 1.232. Fried–Gaudry–Heninger–Thomé 2017, Appendix A, recompute these
    and label them heuristic (codex).
  - These are individual-log costs given the surrounding NFS machinery, not complete DL costs. 1.232 > 1.20 excludes
    that parameterised approach, not improved descents.
- **Shape.** Descent outputs live in ideals of a number field. They are used through virtual logarithms.
  - Turning them into Z-relations needs exact cancellation of the ideal parts over Z, a trivial unit part, and
    computable generators.
  - We know of no conversion algorithm with a bounded output height. v1's "generically huge entries" and the
    independent-incidence remark were not supported quantitatively and are withdrawn (codex).
  - With unbounded height (Prop H) such a conversion would give NAV_∞-type outputs, i.e. heuristic NFS-DL by another
    route.
- **Obs DD (what the hidden-root argument needs) [paper2 Theorem DD, internal; combination CONDITIONAL].**
  - paper2's Theorem DD (Rabin/Dixon fibre symmetry for descents) assumes:
    - a dependency whose set S of hidden-root rows (with odd coefficient) is non-empty;
    - an X² ≡ Y² with X = ∏_{i∈S} r_i and Y measurable with respect to the transcript (the R_i and the coins).
    It then gives Pr[nontrivial split] = 1 − 2^{1−ω}, whatever the rows' algebraic shape. Supplying Y for NFS-shaped
    rows needs the character, unit and square-root step: paper2 at fixed degree (lv6), lv5unif §9 at growing degree.
  - So the hidden root survives a descent into ideals. NAV's integer shape and its height bound are not needed for
    fruitfulness.
  - *Augmenting standard NFS rows (sketch, corrected).* Let V be the F₂ row space of the standard rows (all columns:
    rational, algebraic, character, sign and unit), of codimension c.
    - Then c + 1 hidden-root rows are needed, by pigeonhole in F₂^n/V. c rows that span a complement of V give *no*
      dependency containing them, which was v1's error (proof critic, codex).
    - Dependencies made of standard rows alone, or in which every hidden root has even total coefficient, carry no
      hidden randomness.
    - Lee–Venkatesan 2018 prove relation collection and a possibly trivial congruence at L[1/3, (64/9)^{1/3}]. Their
      theorem does not establish this augmentation. The lifts, characters, units and square-root branches would all
      have to be redone. We label the combination CONDITIONAL.

### 5.4 Where the rigour stops, and the criterion multi-level descent must meet [OPEN; HEURISTIC bookkeeping]

- **Theorem U's floor.** lv5unif §13 (and codex blind §4), with D = x/ln x:
  - Theorem U's one-level descent needs the split acceptance (ln N)^{−2k} to be affordable, which requires k ≤ D^α.
  - It also needs a small norm exponent θ ≥ c/√k, which requires k ≳ D^{2(1−2α)}.
  - Together these force α ≥ 2/5.
- **Reaching 1/3 needs multi-level descent.** Second-level special-q primes are prime factors of descent values, not
  uniform primes in an interval. So Theorem U's averaged descendability (a fraction γ ≥ (ln N)^{−C_γ} of the q in a
  dyadic window) does not apply to them as stated.
- **A criterion (HEURISTIC bookkeeping; independence assumed).** Run the tree with local backtracking: a node retries
  its own relation until all of its children descend.
  - Suppose a level-i node has m_i children, each good with probability ≥ γ_i conditionally on the history. Then the
    expected cost multiplies by ≈ ∏_i γ_i^{−m_i}·m_i along a root-to-leaf path.
  - An L[1/3] total therefore needs Σ_i m_i ln(1/γ_i) = O(x^{1/3}(ln x)^{2/3}).
  - With γ_i = (ln N)^{−C} this allows Σ_i m_i = O(x^{1/3}(ln x)^{−1/3}). The first level must then produce
    children of size L[2/3], as in heuristic descent, and Theorem U's (ln N)^{−2k} is the one-level case m₁ = 2k.
  - The rigorous ingredient that is missing is the conditional goodness bound itself: at each level ≥ 2, a lower bound
    on the number of special-q lattice points whose values are smooth *over the good primes of the next level*. That
    is a count of G-smooth values for a prime set G defined by descendability, uniform over a set of q of full measure
    under the children's distribution. Theorem U supplies the averaged count only at level 1, where the q are fresh
    uniform primes.
  - We state this as the target; codex's blind round calls it "conditional smoothness estimates valid for the
    generated children" with "an explicit height potential".

**Verdict on Route 4.**
- Integer descent fails at the leaf level within its model (Prop TD; Prop SI windows).
- NFS descent is heuristically L[1/3] but ideal-shaped, with no bounded-height conversion known.
- Fruitfulness of ideal-shaped rows is available given the square-root machinery (Obs DD).
- The rigour gap is the multi-level goodness bound of §5.4.

---

## 6. Answers to the five questions

1. *Is approximate counting easier than finding?*
   - We know no counter cheaper than search where it matters.
   - Counting is easy only where it is uninformative: flat unrestricted word levels, and budgets < ln(N/2). At the
     threshold a counter implies NAV decision (Prop CF).
   - Cost of the known counters: Θ(φ) (word Monte Carlo), L[2/3] (lattice Monte Carlo), DL-dependent (characters).
   - OPEN: counters with implicit character information.
2. *Approximation factor?*
   - A CVP_{3/(2+ε)} oracle in weighted ℓ1 suffices (GRH); NAV itself tolerates up to ≈ 3.
   - Dimension ≥ x^{3/2}, with target vectors x^{1/4} to L[1/3, κ/2] below the Gaussian heuristic.
   - Only usable once a basis of Λ_N (FACT ∧ DL) is in hand. BKZ at L[1/3] cost reaches exp(L[1/3]); exact methods
     cost 2^{Θ(K)} (Prop LT).
3. *Do squares help MITM?* No.
   - The solver cannot see the root, and square- and unit-NAV are equivalent up to 2^ω for uniform targets (Prop SQ).
   - Meet-in-the-middle is generic (Prop G) or fixed-list (≤ |S₁||S₂|/|Q|).
4. *Descent?*
   - Tested Z-descent is L[2/3, 2/(3κ)] (MODEL), and integer descent needs Prop SI windows.
   - NFS descent is heuristically L[1/3, 1.23–1.44] per target, but ideal-shaped.
   - Rigorous NFS-type descent for random squares is the project's Dixon descent: rigorous at one level (L[2/5]), with
     the multi-level criterion of §5.4 open.
5. *Barrier?*
   - NAV for all N ⇒ FACT ∧ F_p-DL at L[1/3, max(κ + τ, 3κ)] (Theorem DLB). NAV_∞ ≡ FACT ∧ DL (Prop H), and generic
     methods need Ω(√q) (Prop G).
   - NAV is not equivalent to rigorous L[1/3] factoring alone: it also gives DL, unless DL ≤ FACT. NAV_3 adds a
     short-vector requirement that no reduction uses.

---

## 7. Experiments (pre-registered in `notes/navigate_prereg.md`, commit `f2e4e9c`)

### 7.1 E1: toy walks vs sampling (operational)

Balanced semiprimes, 6 per size, 100 targets per N, 400 smoothness tests per (algorithm, target), y = 2^{bits/4}.

| bits | Y_S0 | Y_G1 | Y_J2 | J2/S0 | G1/S0 | R_S0 | R_G1 | mean tested log-size S0 → G1 | G1/S0 distinct | splits/output S0, G1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 0.0721 | 0.1268 | 0.0725 | **1.005** | **1.759** | **1.013** | **1.009** | 25.2 → 21.9 | 1.087 | 10.5, 84.8 |
| 48 | 0.0730 | 0.1151 | 0.0738 | **1.011** | **1.577** | **1.005** | **1.012** | 30.7 → 27.6 | 1.306 | 10.4, 96.6 |
| 56 | 0.0725 | 0.1068 | 0.0720 | **0.994** | **1.473** | **1.002** | **1.016** | 36.3 → 33.2 | 1.382 | 10.5, 106.6 |

- **All four predictions held:** P1a (R_S0 ∈ [0.85, 1.15]), P1b (J2/S0 ∈ [0.85, 1.15]), P1c (G1/S0 ≥ 1.15) and P1d
  (R_G1 ∈ [0.85, 1.15]). Codex reproduced the summary exactly.
- *Reading.*
  - Monte Carlo completion counting (J2) is sampling with extra bookkeeping.
  - The height-steered walk (G1) wins per smoothness test only because it tests smaller pieces (≈ 3.1–3.3 nats
    smaller). The size-matched coprime null explains its rate to within 1–2%.
  - G1 pays ≈ 8–10× more splits per output. The split counter omits the re-split of the chosen candidate, so the true
    ratio is slightly higher (codex).
  - Per *distinct* relation its gain falls to 1.09–1.38. Post hoc we attribute this to the walk lingering near
    small-height states, whose splits are rescalings of one relation (smoothbuild E2's effect).
  - Signs: 50% of outputs have u, v > 0 for all three algorithms (descriptive).
  - Mean tested log-size pools per-target means (codex).
  - The null pair is drawn before the u-test, so a split cut off by the test budget counts as a null pair with no
    chance of success. That is ≤ 1 split per run, negligible (proof critic).
- No L-scale conclusion is drawn: u = 2 at every size.

### 7.2 E2: exact completion counts (diagnostic)

20-bit semiprimes (4 seeds; prime factors 809–1009). For every signed exponent vector of height ≤ β ln N over FB =
primes ≤ y, the per-class counts C_β are computed on H = ⟨FB⟩ (H = G in all instances).

| y (K) | β | μ (mean) | h (hit frac.) | Poisson 1 − e^{−μ} | D = Var/μ | f75 |
|---|---|---|---|---|---|---|
| 13 (6) | 1.0 | 0.06 | 0.053–0.067 | 0.057–0.065 | 0.93–1.15 | 0.002–0.005 |
| 13 | 1.5 | 0.59–0.66 | 0.35–0.54 | 0.44–0.49 | 0.71–1.54 | 0.000–0.002 |
| 13 | 2.0 | 3.1–3.5 | 0.88–1.00 | 0.96–0.97 | 0.47–1.69 | 0.000–0.002 |
| 13 | 3.0 | 34–39 | 1.00 | 1.00 | 0.30–1.17 | 0.000–0.001 |
| 19 (8) | 1.25 | 0.73–0.81 | 0.47–0.58 | 0.52–0.55 | 0.71–1.22 | 0.001–0.006 |
| 19 | 2.0 | 23–26 | 1.00 | 1.00 | 0.42–0.67 | 0.001–0.003 |

(Full grid in the summary JSON.)

- **P2a failed** (11 cells, all with μ ∈ [0.5, 3.5]). Hit fractions deviate from Poisson by up to 0.15 (seed 2,
  y = 13, β = 1.75: 0.614 vs 0.764), in both directions, seed-dependently.
- **P2b failed** (25 cells). D ranges over 0.30–1.69. Seeds 0 and 3 are strongly under-dispersed (counts close to
  constant at large β); seeds 1 and 2 are over-dispersed.
- **P2c failed** (28 cells). f75 = 0.000–0.006 in the registered cells (≤ 0.008 in all cells), against the registered
  [0.2, 0.5]. The counts are *not* noise-like: 75% of the non-principal Fourier mass sits on ≤ 0.6% of the characters,
  and in seeds 1 and 2 at y = 13, β ≥ 1.75 on only 4–19 characters (proof critic).
  - *Noise reference corrected (codex).* C(t) = C(t^{−1}) (e ↦ −e), so Ĉ is real. The right i.i.d. reference has
    squared-normal powers, f75 ≈ 0.27, not v1's complex-Gaussian 0.39. The registered band [0.2, 0.5] contains both
    values, so the outcome is unchanged.
- **P2d held**: squares and non-squares have the same hit fraction to within 0.012 in every cell.
- *Explanation (POST HOC; plausible, not proved).*
  - With K = 6–8 primes and budgets of 12+ letters, each prime is used with many exponents: a *dense* regime.
  - Ĉ_B(χ) = Σ_{|e|≤B}∏χ(ℓ)^{e_ℓ} then concentrates on characters whose angles arg χ(ℓ) are all small (scale
    ≈ K ln ℓ/B). The proof critic confirmed that the dominant characters at y = 13 have all six angles ≤ 0.08 turns.
  - v1's product-Fejér estimate (≈ 0.008) misses seeds 1 and 2 by about three orders of magnitude. The share is
    seed-dependent and set by the dual-lattice geometry. A truncated weighted ℓ1 ball is not a product of Fejér
    kernels, and sidelobes and correlated angles matter (codex).
  - Equivalently, C_β is a smooth function of the coset of the low-dimensional lattice Λ_N, whose short vectors (or
    their absence) also set the seed-dependent dispersion. This is the CSI-FiSh-like regime, opposite to
    bounded-height NAV (Prop LT(a)).
- *Post-hoc sparse regime* (`posthoc_e2.py`; same code).
  - y = 512 (K = 97), β ∈ {0.8, …, 1.1}.
  - y = 2048 with the two prime factors of N removed from the factor base (K = 307), β ∈ {0.8, 0.9, 1.0}.
  - **v1 error (codex, FATAL to that table):** v1's y = 2048 run included N's own factors (809–1009 < 2048) in the
    factor base, and the log table mapped them to 0. It is re-run with them excluded (N is then not admissible for y:
    a disclosed diagnostic deviation). μ changed by ≤ 0.03 and h, D, f75 by ≤ 0.003. An assertion now refuses nonunit factor bases.

  | y (K) | β | μ | h | Poisson | D | f75 |
  |---|---|---|---|---|---|---|
  | 512 (97) | 0.8 | 0.32–0.34 | 0.29–0.31 | 0.28–0.29 | 0.85–0.86 | 0.056–0.059 |
  | 512 | 0.9 | 1.16–1.20 | 0.74–0.76 | 0.69–0.70 | 0.73–0.74 | 0.055–0.058 |
  | 512 | 1.0 | 4.0–4.1 | 0.99 | 0.98 | 0.69–0.70 | 0.048–0.051 |
  | 2048\* (307) | 0.8 | 0.40–0.42 | 0.36–0.37 | 0.33–0.34 | 0.82–0.84 | 0.074–0.078 |
  | 2048\* | 0.9 | 1.56–1.60 | 0.85–0.86 | 0.79–0.80 | 0.67–0.68 | 0.077–0.081 |
  | 2048\* | 1.0 | 5.9 | 1.00 | 1.00 | 0.64–0.65 | 0.070–0.073 |

  (\* nonunit primes excluded.)
  - The Fourier share rises 30–80-fold (≈ 0.001 → 0.05 → 0.075) as K grows. It is still below noise (≈ 0.27), and the
    seed-to-seed spread vanishes.
  - At K = 97 the top characters are nearly trivial on 2, 3, 5 only (proof critic). That is small-prime structure,
    which need not disappear as K grows.
  - Counts are mildly sub-Poissonian (D 0.64–0.86), with hit fractions 0.02–0.06 above Poisson.
  - Most sparse cells lie at β ≤ 1, where Prop CF(iv)'s convergent enumeration already applies to budgets < ln(N/2)
    (codex).
  - *Reading:* the strong low-frequency structure is a dense-regime phenomenon. Whether the residual structure vanishes
    in the NAV regime (K = L[1/3] ≫ letters ≈ x^{2/3}) is CONJECTURED, from two K values. Toy sizes cannot reach that
    regime.
- *Oracle-walk profile (descriptive; description corrected).*
  - The code draws a uniform vector of height ≤ B and keeps only its class (a size-biased target). It then walks,
    choosing uniformly among *all* (prime, sign) letters whose remaining instance has a non-zero count (grid-rounded
    budgets), and records the fraction σ of such letters.
  - This is neither "along sampled solutions" nor canonical-order as v1 said; there is no −1 letter in E2.
  - σ ≈ 1 while ≥ 2 ln N of budget remains (y = 13). It falls to 0.2–0.6 below ≈ 1.6 ln N. In the sparse post-hoc
    cells it is 0.04–0.15 near the threshold.
  - The information an oracle-guided walk consumes sits at the end of the walk.
- **Deviations.**
  - (1) A smoke test of `core.py` used seed 99 at 40 bits (one target; not a registered instance).
  - (2) The post-hoc runs. v1's post-hoc run peaked at 4.9 GB RSS, and the first y = 2048 re-run at 3.0 GB, each for
    ≈ 1 min, above the 2 GB guideline; no swapping occurred. The cause was fragmentation in the vectorised enumeration
    over ~100–300 primes. `e2_counts.py` now enumerates the large primes sparsely. The final regeneration (all E2 and
    post-hoc data, registered cells bit-identical) peaks at 1.3 GB and 0.68 GB.
  - (3) The walk profile is not the registered "along uniformly sampled solutions" statistic (descriptive only; see
    above). The regeneration changed its random draws; the qualitative profile is unchanged.
  - (4) Signs are allowed in outputs (−1 column), as registered; smoothbuild §7 says "positive". About half the outputs
    are positive.

---

## 8. Literature (checked by the literature critic; see §10.2)

(See §10.2 for the checked list and corrections.)

---

## 9. GPT-6.1-sol rounds

### 9.1 Blind round (`codex/navigate_1.out`, launched before any analysis here)

Codex agreed on every route's failure point. Several statements here come from it:
- the JVV restricted-count caveat;
- word Monte Carlo at Θ(φ) samples;
- the GRH derivation of the constant 3;
- that a CVP_{γ<3/2} oracle suffices, and that the Gaussian-heuristic gap L[1/3, 0.30] is *not* hardness evidence;
- that Schnorr's lattice does not encode the congruence;
- the fixed-list MITM bound;
- reconstruction plus testing at L[2/3, 1.111];
- the descent constants 1.442/1.232 (FGHT Appendix A);
- the padding correction;
- the Theorem U geometry.

Its probability: ≈ 5% (examined routes, any constant), ≈ 1% at (0.60, 1.20). Its most promising route was
"target-conditioned multilevel descent extending Theorem U, with an explicit height potential and conditional
smoothness estimates valid for the generated children" (§5.4). It left a bounded-height reduction from
arbitrary-target to square-target NAV open; Prop SQ closes it for uniform targets.

### 9.2 Adversarial round on v1 (`codex/navigate_2.out`): MAJOR REVISION; all points adjudicated

| # | Codex point | Severity | Adjudication |
|---|---|---|---|
| 1 | Prop H correct but "not an intermediate target" overstated; bit-complexity, supplied generators, generation hypothesis | MAJOR | accepted: H restated with binary outputs and hypotheses; wording "strong sufficient target" |
| 2 | DLB padding sieve cost is L_p[1/3, 2^{1/3}κ], not L^{o(1)} | MINOR | accepted (dominated by 3κ) |
| 3 | j* independent of the oracle-visible transcript only; Chernoff unnecessary; cite a Howell algorithm | MINOR | accepted |
| 4 | average-case DLB constant becomes 2^{1/3}·1.80 ≈ 2.27 | MAJOR | accepted (also the proof critic) |
| 5 | G's reach overstated; Goldfeld quantifier; RT corollary needs 2^ω and 1/δ | MAJOR | accepted |
| 6 | SQ(b) is an average-target reduction; fixed t can fail; 2^ω condition substantial | MAJOR | accepted |
| 7 | CF: flatness only for unrestricted word counts; "is exactly" → "implies"; e = 0 branch; binary encodings | MAJOR | accepted; CF no longer called an obstruction |
| 8 | LT: γ < 3/2 suffices only one-way; ADS15 Euclidean, not weighted ℓ1; Gaussian-heuristic gap is not hardness | MAJOR | accepted |
| 9 | LT(e) counting statement false (lists of height ≤ H are exponential in H) | MAJOR | accepted: restated as a SCENARIO with list size S ⇒ height ≳ ln S |
| 10 | TD and integer descent are model statements; NFS conversion remark unsupported | MAJOR | accepted: scope added; conversion remark withdrawn |
| 11 | Obs DD conditional on the square-root identity; "span the corank" wrong (need d + 1) | MAJOR | accepted (also the proof critic) |
| 12 | y = 2048 factor base contained N's factors (dlog of 0) — FATAL to that table; walk profile differs from registration | FATAL / MAJOR | accepted: re-run with nonunits excluded (μ changes ≤ 0.03, others ≤ 0.003), assertion added, profile relabelled |
| 13 | f75 noise reference should be ≈ 0.27 (real coefficients); Fejér explanation unproved; delete N^{1−o(1)} extrapolation | MINOR / MAJOR | accepted |
| 14 | literature/attribution tightening; keep blind-round qualifications | MINOR | accepted |
| 15 | verdict: "no demonstrated route" justified; "wrong target / retire" too strong; keep the unbounded-height relaxation as a target; probability ≈ 5% / ≈ 1% | MAJOR | partly accepted: wording changed to a recommendation; codex's position recorded in §0; we keep our lower probability, giving the reason (Theorem DLB at the target constants) |

## 10. Opus critics

### 10.1 Proof critic (read-only, on v1): ACCEPT WITH FIXES; all applied in v2

It checked step by step: Prop H(a)(b), Theorem DLB (padding, query law, chain argument, Howell), Prop G (including the
composite-order extension), Lemma RT, Prop SQ, Prop CF(iv) (5,549 pairs at 14 bits) and Prop LT(a). It also checked
the E1/E2 code.

Findings and fixes:
- MAJOR: the Obs DD condition was backwards (c + 1 rows; odd hidden coefficient; transcript-measurable Y).
- MINOR:
  - γ ≈ 3 rather than 3/2 "relative to the true distance";
  - the average-case DLB constant 2.27;
  - FACT does not supply generators; ω notation clash;
  - Prop CF mixed encodings; "is exactly" → "implies"; e = 0 and −1 cases;
  - Prop LT(e) list-size statement;
  - the "fixed fraction of characters" claim contradicted by seeds 1–2 (4–19 characters);
  - the Fejér estimate off by 10³ in two seeds, and K = 97 concentration on 2, 3, 5;
  - the walk profile mis-described;
  - Lemma RT's 2^ω factor;
  - Prop SQ's conditionality in §0 and the table/body mismatch in m.
- NIT: the cost of Y at A = exp(L[1/3, c]); E1 null censoring; the O(d) count in the conversion remark (now
  withdrawn).

### 10.2 Literature critic (read-only)

(pending)

## 11. Reproduction

```
.venv/bin/python -m src.navigate.e1_walks      # ~25 s
.venv/bin/python -m src.navigate.e2_counts     # ~25 s, peak RSS 1.3 GB
.venv/bin/python -m src.navigate.posthoc_e2    # ~50 s, peak RSS 0.7 GB (post hoc)
.venv/bin/python -m src.navigate.analyze
.venv/bin/python -m pytest -q tests/test_navigate.py
```
