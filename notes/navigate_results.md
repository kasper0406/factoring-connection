# navigate: hidden-root short-quotient navigation (smoothbuild §7, W11)

Agent `navigate`, phase 25, 2026-10-23. Version 1 (before the codex adversarial round and the Opus critics; §9–§10
will record them).

**Files.**
- Pre-registration: `notes/navigate_prereg.md` (commit `f2e4e9c`, before any registered statistic).
- Code: `src/navigate/`.
  - `core.py`: operational helpers (factor base, Gauss split, batch smoothness test, exact relation check) and the
    three E1 algorithms S0, G1, J2.
  - `e1_walks.py` (E1 harness), `e2_counts.py` (E2, diagnostic: exact per-class counts, Fourier spread), `posthoc_e2.py`
    (POST HOC sparse-regime E2), `analyze.py`.
- Tests: `tests/test_navigate.py` (7 tests).
- Data: `experiments/raw/navigate_E1.jsonl`, `navigate_E2.jsonl`, `navigate_E2_posthoc.jsonl`; summary
  `experiments/analysis/navigate_summary.json`.
- GPT-6.1-sol: blind round `codex/navigate_1.{md,out}` (launched before any analysis here); adversarial round
  `codex/navigate_2.{md,out}`.

**Labels.** KNOWN (published; cited) · PROVED (proof here; elementary unless said) · PROVED (generic model) ·
HEURISTIC / MODEL (named assumption) · MEASURED · POST HOC · OPEN · assessment. "No prior statement found by us" is the
strongest novelty claim made.

**Notation.** x = ln N; L[α, c] = exp((c + o(1))x^α(ln x)^{1−α}); y = L[1/3, κ] (κ = 0.60 in the target), K = π(y);
G = (Z/N)^×, Q = G², ω = ω(N); FB = {−1} ∪ {primes ≤ y}. N is *admissible*: odd, ω ≥ 2, no prime factor ≤ y. Exponent
vectors e ∈ Z^K are identified with coprime pairs (u, v) = (∏_{e_ℓ>0}ℓ^{e_ℓ}, ∏_{e_ℓ<0}ℓ^{−e_ℓ}); the *height* is
|e| := Σ|e_ℓ| ln ℓ = ln u + ln v. Λ_N := {e : ∏ℓ^{e_ℓ} ≡ 1 (mod N)} is the relation lattice.

---

## 0. Verdict

**Plain language.**
- **No route reaches the target, and the target turns out to be the wrong one.** Navigation at L[1/3, 1.20] would
  factor N rigorously in L[1/3, 1.80]. Each of the five routes fails at a named point (table below). What the
  analysis adds is structural:
  - **The height bound 3 ln N does no work.** Dixon's argument (EQ1) and the DL reduction (EQ3) run unchanged with
    heights up to exp(L[1/3]) (Prop H). With unbounded heights the problem is equivalent to "factoring plus discrete
    logarithms modulo each prime factor" (Bach 1984; smoothbuild EQ2). So NAV is not an intermediate step toward
    rigorous L[1/3] factoring. It contains rigorous L[1/3] factoring *and* rigorous L[1/3] prime-field DL, and its
    bounded form adds a closest-vector problem that buys nothing.
  - **A navigator at the target constants would beat the best *heuristic* algorithms for both problems.**
    Theorem DLB makes EQ3 quantitative: NAV(κ = 0.60, τ = 1.20) for all admissible N gives Las Vegas
    discrete logarithms in every F_p^× in L_p[1/3, 1.80]. Heuristic NFS/MNFS-DL is 1.923/1.902; no rigorous
    L_p[α < 1/2] DL algorithm for prime fields is published (the project's L[2/5] is internal and unrefereed).
  - **Bounded-height NAV has no known L[1/3] algorithm even heuristically, even for prime moduli, even when the
    factorisation and all discrete logarithms are given.** With all discrete logs known, NAV is a weighted-ℓ1 closest-vector problem
    in the relation lattice. Its dimension must be ≥ x^{3/2−o(1)} just for solutions to exist (Prop LT), and the
    solutions lie far below the Gaussian-heuristic radius. The best methods we know cost L[2/3, 2/(3κ)] = L[2/3, 1.11] (testing,
    Prop TD), N^{1/4−o(1)} (generalised birthday at bounded height, Prop LT(e)) or 2^{Θ(K)} (lattice algorithms).
  - **Rigorous L[1/3] factoring does not need NAV.** The project's Theorem DD (Dixon descent, paper2) already makes
    NFS-shaped descended rows fruitful, because the hidden root survives a descent into ideals. What is missing is a rigorous
    *multi-level* descent: Theorem U's one-level geometry forces exponent 2/5 (lv5unif §13). Codex's blind round
    independently named the same place as the most promising route.
- **Barrier evidence (no unconditional lower bound).**
  - Generic algorithms need Ω(√(qδ)) group operations, with q the largest prime factor of φ(N) (Prop G; Shoup/Nechaev).
  - Every relation for a random target has height ≥ x − O(ln x) (Lemma RT). Obtaining that much smoothness by testing
    costs L[2/3] at y = L[1/3] (Theorem SS model).
  - Counting completions is at least as hard as finding one at the threshold budget, and trivial elsewhere (Prop CF).
  - Squares are not special: square- and unit-targets are equivalent up to 2^{ω} calls (Prop SQ).
- **Experiments (toy scale; mechanism checks only).**
  - E1: all four predictions held. Monte Carlo completion counting adds nothing per smoothness test (0.99–1.01).
    Steering by the cheap Gauss-height signal gains 1.47–1.76 per test, entirely explained by smaller tested pieces
    (size-matched null 1.01–1.02); per *distinct* relation the gain is 1.09–1.38.
  - E2: three of four predictions failed. The informative failure: with 6–8 primes, the exact counts are carried by a
    very few characters (0.0–0.8% of them), not by noise. That is the dense, low-dimensional regime of CSI-FiSh-type
    lattices. A post-hoc sparse-regime run (97 and 309 primes) raises the share to 4–8%, trending toward noise.
    NAV at bounded height is necessarily in the sparse regime.
- **Assessment.** Probability that NAV at L[1/3, O(1)] is achieved rigorously: ≤ 1% (ours). Codex (blind): ≈ 5%
  at any constant, ≈ 1% at (0.60, 1.20). Probability that a rigorous L[1/3] factoring algorithm comes *through the NAV
  interface* rather than through multi-level Dixon descent: lower still. We recommend retiring W11 as a target in favour
  of multi-level Dixon descent.

**Route table.**

| Route | Cost / status | Exact failure point |
|---|---|---|
| 1. Counting-guided walks (JVV) | reduction NAV ≤ restricted counting PROVED; no counter below search cost known | the counter is trivial above (2 + ε)x (GRH flatness) and below ln(N/2) (convergents); in between, at the existence threshold, approximate counting *is* NAV decision. Monte Carlo: 1/φ per word sample, L[2/3] per lattice sample; characters need DL and (toy E2) a constant fraction of all characters |
| 2. Lattice / CVP | needs a basis of Λ_N (≡ FACT ∧ DL) and CVP_γ, γ ≈ 3/2 (weighted ℓ1), in dimension ≥ x^{3/2−o(1)}, far below the Gaussian-heuristic radius | no known algorithm below 2^{Θ(K)}; BKZ-β reaches only exp(Θ(K ln β/β)); Schnorr's real-log lattice does not encode u − av = kN |
| 3. Meet-in-the-middle | fixed lists: success ≤ \|S₁\|\|S₂\|/\|Q\|; any generic method Ω(√(qδ)) | collisions are generic; the square/hidden-root promise changes counts by ≤ 2^ω and is invisible to the solver |
| 4. Descent | tested Z-descent L[2/3, 2/(3κ)]; NFS descent heuristic L[1/3, 1.23–1.44] per target but ideal-shaped | Z-relations have height ≥ x − O(ln x) (Lemma RT/DSC), so a cheap leaf level needs Prop SI windows; NFS rows convert to Z only with unbounded heights (fine by Prop H, but then NAV_∞ ≡ FACT ∧ DL); rigour gap = multi-level descent |
| 5. Barrier | NAV ⇒ FACT ∧ F_p-DL at L[1/3, max(κ + τ, 3κ)] (Theorem DLB); generic Ω(√q) | no unconditional lower bound; the evidence is that NAV at the target constants would beat heuristic NFS for factoring and for prime-field DL |

**Main results.**

| Result | Statement (short) | Status |
|---|---|---|
| Prop H | EQ1 and EQ3 accept heights up to exp(L[1/3]) at the same L[1/3] cost; NAV_3 ⇒ NAV_∞ ≡_R FACT ∧ DL(p^a ∥ N) | PROVED (bookkeeping on EQ1/EQ2) |
| Theorem DLB | NAV_A(κ, τ, δ) for all admissible N ⇒ Las Vegas F_p^× DL in poly(1/δ)·L_p[1/3, max(κ + τ, 3κ)]; average-case version over random semiprimes | PROVED (EQ3 made quantitative; padding as dsr Lemma PT, rank step as smoothbuild EQ2) |
| Prop G | generic group of order with prime factor q: representing a random target by K given generators needs (T + K + 2)² ≥ (δ − 1/q)·2q queries | PROVED in the generic group model; standard (Nechaev 1994, Shoup 1997) |
| Lemma RT | #{t ∈ G : some u ≡ tv with 0 < \|u\|\|v\| ≤ X} ≤ 2X(1 + ln X) | PROVED (elementary) |
| Prop SQ | square-target NAV(δ) ⇒ unit-target NAV(2^{−ω}δ) unconditionally; under GRH, equivalence up to O(2^ω ω log(1/ε)) calls and height + ln y when 2^ω ≤ y^{1/2−ε} | PROVED (GRH for the second part) |
| Prop CF | NAV is self-reducible (depth ≤ k_max, branching 2K); approximate restricted counting ⇒ search; counts flat above (2 + ε)x (GRH), exact below ln(N/2) | PROVED (flatness under GRH, for word counts) |
| Prop LT | a δ-successful NAV_A needs π(y) ≥ (Ax)^{A/(A−1)−o(1)}; NAV_3 needs K ≥ x^{3/2−o(1)}, i.e. the sparse regime | PROVED (from de Bruijn's Ψ(X, (ln X)^c) = X^{1−1/c+o(1)}) |
| Prop TD | tested Z-relations of height ≥ x − O(ln x) at y = L[1/3, κ] cost L[2/3, 2/(3κ)] | MODEL (H_S on tested pieces) |
| Obs DD | Theorem DD (paper2) gives fruitfulness for ideal-shaped descended rows; NAV's Z-shape is not needed for the hidden-root argument | PROVED in paper2 (internal); observation here |

---

## 1. Setting

**NAV_A(y; T, δ).** A randomized algorithm that, for every admissible N, on input (N, a) with a a uniformly random
element of Q (its root withheld), outputs in time T an exponent pair (e, f) ∈ Z^K × Z^K with
∏ℓ^{e_ℓ} ≡ a·∏ℓ^{f_ℓ} (mod N) and |e| + |f| ≤ A·x, with probability ≥ δ over a and its coins. NAV_∞ drops the height
bound (exponents of L^{O(1)} bits). The target of smoothbuild §7 is NAV_3 with y = L[1/3, 0.60], T = L[1/3, 1.20],
δ ≥ 1/poly(x). Outputs are given as exponent vectors; u and v need never be written out.
- *Existence.* Under GRH, for k = ⌈(2 + ε)x/ln y⌉ the k-step walk on G with uniform signed letters ℓ^{±1} (ℓ ≤ y)
  reaches every class with probability ≥ 1/φ − λ^k > 0, λ ≤ y^{−1/2+o(1)} (JMV 2009; dsr Prop RW). So every class has
  a word of height ≤ k ln y < 3x for large N. This is where the constant 3 comes from (codex blind §1).
- *Variants.* The "worst-case N" quantifier matters only for Theorem DLB (padding); EQ1 needs NAV only for the N being
  factored.

---

## 2. Reductions and the barrier (Route 5)

### 2.1 Prop H: the height bound is not load-bearing [PROVED]

**Statement.** Let A ≤ exp(L[1/3, c]) for some c, possibly much larger than 3.
- (a) EQ1 with NAV_A costs O((K + 1)(T + x²)/δ) + (F₂ linear algebra) + O(K·log(KAx)·x²), with the same success
  probability. The F₂ step is O(K²·Ax) by block Wiedemann when rows have ≤ Ax/ln 2 odd entries, and O(K^ω) always.
- (b) EQ2 and EQ3 (discrete logarithms) accept NAV_A outputs unchanged; the linear algebra is over Z/M and sees
  exponents only modulo M.
- (c) Hence NAV_3 ⇒ NAV_∞, and NAV_∞ ≡_R FACT ∧ DL(p^a) for every p^a ∥ N (smoothbuild EQ2; Bach 1984), given that FB
  generates G. The equivalence costs: NAV_∞ at L[1/3, τ] ⇒ FACT at L[1/3, max(κ + τ, ωκ)] and DL(p^a) at
  L[1/3, max(κ + τ, 3κ)]. Conversely, FACT and DL oracles of cost L[1/3, c] give NAV_∞ at L[1/3, κ + c] precomputation
  (logs of FB) plus one DL and a K-term linear solve per call.

*Proof.* (a) EQ1 (smoothbuild §4) uses the outputs in two places. The F₂-dependency uses only exponent parities.
Computing Y = ∏ℓ^{E_ℓ/2} mod N with E_ℓ ≤ (K + 1)Ax/ln 2 costs O(log E_ℓ) multiplications per ℓ by repeated squaring.
smoothbuild's term K·A·x³ is the naive cost and is not needed. The fruitfulness argument (Theorem DD of paper2;
smoothbuild EQ1) does not look at heights. (b) EQ2's rank step works with the vectors v_i ∈ (Z/M)^K. (c) EQ2 and the
standard Bach-type argument; the converse direction solves Σ e_ℓ log ℓ ≡ log a modulo the group exponent per prime
power, which is possible by generation. ∎

**Reading.**
- The payoff of NAV (rigorous factoring) is the payoff of NAV_∞, and NAV_∞ *is* FACT ∧ DL. So NAV is not a weaker
  intermediate target: an L[1/3] NAV algorithm contains a rigorous L[1/3] factoring algorithm and rigorous L[1/3]
  discrete logs modulo the prime factors.
- The height bound makes NAV_3 *harder* than NAV_∞ (a short coset vector, §4) without improving any consequence.
- Heuristically NAV_∞ is L[1/3, 1.923] (NFS factoring, then NFS-DL modulo each prime; or NFS relations, exact integer
  elimination of ideals and units, then Z-relations of huge height, §5.3). Nothing heuristic supports τ = 1.20.

### 2.2 Theorem DLB: NAV gives prime-field discrete logarithms [PROVED; EQ3 made quantitative]

**Statement.** Suppose NAV_A(y = L_N[1/3, κ]; T = L_N[1/3, τ], δ) holds for every admissible N, with A ≤ exp(L[1/3]).
Then there is a Las Vegas algorithm that, given a prime p, a generator g of F_p^× and h ∈ F_p^×, returns log_g h in
expected time poly(x, 1/δ)·(L_p[1/3, κ + τ] + L_p[1/3, 3κ]). At (κ, τ) = (0.60, 1.20) this is L_p[1/3, 1.80].

*Proof.*
- *Padding (dsr Lemma PT).* Let r be the least prime > max(y(p²), 3 ln² p²), found by sieving in time L^{o(1)}.
  N := pr is admissible, ln N = (1 + o(1))ln p, and y(N) = L_p[1/3, κ + o(1)]. (A fixed small r would not be
  admissible; codex.)
- *Queries.* Put M = p − 1. If h is not a square mod p (Euler's criterion), replace h by gh and subtract 1 at the end.
  Choose s₀ queries t_i = CRT(g^{2k_i} mod p, s_i² mod r), with k_i uniform mod M and s_i uniform in (Z/r)^×; at one
  uniformly random position j* put t_{j*} = CRT(h·g^{2k} mod p, s² mod r). Every t_i, including t_{j*}, is a uniform
  element of Q, and the joint law of the query multiset does not depend on j*.
- *Relations.* A success at query i gives (e_i, f_i); reducing mod p, v_i := e_i − f_i satisfies v_i·L ≡ 2k_i (mod M),
  where L = (log_g ℓ)_ℓ, and v_{j*}·L ≡ log h + 2k.
- *Rank step (smoothbuild EQ2, exchangeability).* Over Z/M, the vectors that lie outside the span of the other
  successful vectors form a strict chain of submodules, so at most K·log₂M of them exist. The oracle's view does not
  depend on j*, so Pr[j* is one of them | s successes] ≤ K log₂M/s. Otherwise v_{j*} = Σ_{i≠j*} c_i v_i over Z/M
  (Howell form; M need not be factored), and log h ≡ Σ c_i·2k_i − 2k (mod M).
- *Cost.* j* is independent of the whole transcript (queries and oracle answers), because the query sequence has
  the same law for every j*. A run with s₀ = ⌈4K log₂M/δ⌉ queries has ≈ δs₀ successes (Chernoff), and the target query
  is a successful, dispensable one with probability ≥ δ − K log₂M/s₀ − o(1) ≥ δ/2. So O(1/δ) independent runs suffice:
  O(K x/δ²) oracle calls at time T each, i.e. poly(x, 1/δ)·L[1/3, κ + τ]. The Howell form of an s₀ × K matrix over
  Z/M costs O(s₀K²) ring operations: L[1/3, 3κ]·poly(x, 1/δ). The answer is verified by one exponentiation (Las
  Vegas). ∎

**Average-case version.** If NAV succeeds on average over (N, a) with N a random balanced semiprime (each prime
uniform of its size), then for a ≥ δ/2 fraction of primes p of that size it succeeds on a ≥ δ/2 fraction of N = pq;
running the proof with N = pq for fresh random primes q gives DL for those p.

**Reading.**
- Rigorous prime-field DL at L_p[1/3, c] is open for every c. The best published rigorous bound is Pomerance 1987,
  L_p[1/2, √2] (the project's L[2/5] and paper2's L[1/2, < 1] are internal and unrefereed). Heuristic NFS-DL is
  L_p[1/3, 1.923], MNFS-DL 1.902 (Matyukhin 2003; Barbulescu–Pierrot 2014).
- So NAV(0.60, 1.20) for all N would improve the *heuristic* state of the art for prime-field DL and for factoring at
  once. This is the strongest evidence against the target constants. It is evidence, not a lower bound.
- A navigator that works only for the one N to be factored escapes Theorem DLB, since padding needs N = p·r. It still
  gives DL modulo p | N by EQ2 (Prop H(c)).

### 2.3 Prop G: generic algorithms [PROVED in the generic group model; standard]

**Statement.** In Shoup's generic group model for Z/q (q prime), let an algorithm receive the encodings of
x₁, …, x_K and z, with all of them independent and uniform, make T group-operation queries, and output e ∈ Z^K. Then
Pr[Σe_i x_i = z] ≤ (T + K + 2)²/(2q) + 1/q.

*Proof (standard).* Simulate with formal variables X₁, …, X_K, Z. Every encoded element is a known linear polynomial
in (1, X, Z). The real and simulated games differ only if two distinct polynomials agree at the random point, which
has probability ≤ C(T + K + 2, 2)/q (each difference is a nonzero linear polynomial). In the simulated game e is
independent of (x, z), and Σe_iX_i − Z ≠ 0, so the success probability is 1/q. ∎ (Nechaev 1994; Shoup 1997, whose
Thm 1 covers cyclic groups of order n with largest prime factor q; the representation version is the same proof.)

**Reading.**
- Any NAV method that only multiplies and compares residues (collision search, meet-in-the-middle, walks steered by
  equality tests) needs Ω(√(qδ)) operations. For at least half of all primes p, the largest prime factor of p − 1
  exceeds p^{1/2} (Goldfeld 1969; from memory), so for most balanced semiprimes q ≥ N^{1/4} and the generic cost is
  ≥ N^{1/8−o(1)}.
- So an L[1/3] navigator must use the *integer* structure of residues, in practice smoothness of representatives:
  index calculus. The model treats factor-base logarithms as random. It says nothing about non-generic methods.

### 2.4 Lemma RT: relations for random targets are full-size [PROVED]

**Statement.** For X ≥ 1, #{t ∈ G : ∃(u, v), v ≠ 0, gcd(v, N) = 1, u ≡ tv (mod N), 0 < |u|·|v| ≤ X} ≤ 2X(1 + ln X).
Hence for a uniform t ∈ G (or Q) and f ≥ 1, Pr[some relation for t has ln|u| + ln|v| ≤ x − ln f] ≤
2(1 + x)·N/(f·φ(N)) (times [G : Q] for t ∈ Q).

*Proof.* Each pair (u, v) determines t = uv^{−1}, and (u, v), (−u, −v) give the same t. The pairs with v ≥ 1, u ≠ 0,
|u|v ≤ X number 2Σ_{v≤X}⌊X/v⌋ ≤ 2X(1 + ln X). ∎ (Checked exhaustively at 16 bits in `test_lemma_RT_random_target_size_floor`.)

**Reading.** With f = x², every relation for a random target has height ≥ x − 2 ln x with probability 1 − O(1/x). This
is the random-target form of smoothbuild's Lemma DSC (which needs a small target). So a NAV output always carries
≈ x of log-size that must be y-smooth. Either it is tested, which is Prop TD, or it is constructed, which is
navigation proper.

### 2.5 Prop SQ: squares are not special [PROVED; second part under GRH]

**Statement.** Let S be a NAV solver with success δ on uniform elements of Q, and U one with success δ on uniform
elements of G.
- (a) Unconditionally, S applied to a uniform t ∈ G succeeds with probability ≥ 2^{−ω}δ, at the same height.
- (b) Assume GRH for the real characters mod N and 2^ω ≤ y^{1/2−ε}. Given t ∈ G, call S on t·q_i for
  m = ⌈2^{ω+1} ln(1/ε)⌉ independent uniform primes q_i ≤ y and output (e, f + e_{q_i}) on the first verified success.
  This succeeds with probability ≥ (1 − ε)δ, with height + ln y. Symmetrically, U on a·q_i (a ∈ Q, same m) succeeds
  with probability ≥ (1 − ε)δ.

*Proof.* (a) t is uniform on G, so Pr[t ∈ Q] = 2^{−ω} exactly, and conditionally t is uniform on Q. (b) Condition on
the class c of t in G/Q; t is uniform on its coset. For a fixed q in class c^{−1}, t·q is uniform on Q, so the inputs
t·q_i with class(q_i) = c^{−1} are each uniform on Q. The union of the translates S_good·q_i^{−1} is a subset of the
coset of density ≥ δ. Hence Pr[success | c] ≥ δ·Pr[∃i : class(q_i) = c^{−1}]. Under GRH every nonprincipal real
character χ mod N has |E χ(q)| ≤ y^{−1/2+o(1)} for a uniform prime q ≤ y (JMV's bound), so each class has
probability ≥ 2^{−ω}(1 − 2^ω y^{−1/2+o(1)}) ≥ 2^{−ω−1}, and m tries miss c^{−1} with probability
≤ (1 − 2^{−ω−1})^m ≤ ε. For U: on input a·q_i the input is uniform on the coset class(q_i); U's success density δ_c on
coset c averages δ over cosets, and with probability ≥ 1 − ε a best coset c* (δ_{c*} ≥ δ) is visited; the union of the
translates over the visits to c* has density ≥ δ_{c*}. ∎

**Reading.** The hidden root is invisible to the solver: (N, a) has the same law whether a is drawn as r² or as a
uniform square. Its only role is EQ1's fruitfulness. Up to 2^ω calls, NAV for squares and for arbitrary targets are
the same problem. This answers codex's open question (blind §3) about a bounded-height reduction. 2^ω is not
L^{o(1)} for every N (codex), but it is 4 for semiprimes.

---

## 3. Route 1: counting-guided walks

**Prop CF (counting vs finding).**
- **(i) Self-reducibility [PROVED].**
  - *Word encoding:* the instance (t, k) asks for a signed word of length k over FB with product t; its solutions are
    ⊔_{ℓ,±}(ℓ^{±1}, solutions of (tℓ^{∓1}, k − 1)). Depth k ≤ k_max = Ax/ln 2, branching 2K.
  - *Vector encoding* (canonical: largest prime of the support first): (t, B, j) asks for e supported on the first
    j − 1 primes with |e| ≤ B. Its solutions are ⊔_{i<j, k≠0} solutions of (t·p_i^{−k}, B − |k| ln p_i, i).
    Depth ≤ B/ln 2, branching ≤ 2KB/ln 2. The plain count C(t, B) is not closed under these restrictions; the restricted
    counts are what JVV needs (codex blind §1).
- **(ii) Counting ⇒ search [PROVED; JVV 1986, decision-to-search].** An oracle that approximates the restricted counts
  within any factor decides emptiness. Walking down the self-reduction tree and keeping any non-empty child finds a
  solution with ≤ depth × branching = poly(K, x) oracle calls. With relative error 1 + 1/poly per level, JVV gives almost-uniform
  sampling.
- **(iii) The top levels carry no information [PROVED under GRH].** For word counts W_k(t) (length-k signed words
  with product t), W_k(t) = (2K)^k(1/φ)(1 + θ), |θ| ≤ φλ^k, λ ≤ y^{−1/2+o(1)} (JMV; Fourier inversion). For
  k ≥ (2 + ε)x/ln y the oracle answers "≈ mean" for every child, so a uniform choice is as good as the oracle's.
  Mixing also means that a random prefix maps a uniform target to a near-uniform target, which it already was.
- **(iv) The bottom levels are trivial [PROVED].** If |u|·v < N/2 and u ≡ tv (mod N), then u = tv − kN and
  |t/N − k/v| = |u|/(Nv) < 1/(2v²), so k/v is a convergent of t/N (Legendre). For budgets B < ln(N/2), all solutions
  come from the O(x) convergents, and exact counting costs O(x) smoothness tests.
- **(v) Consequence.** All the information sits in the window B ∈ [ln(N/2), (2 + ε)x]. Near the existence threshold
  B*, where the mean count is ≈ 1, an approximate counter is exactly a NAV *decision* oracle at budget B*, and by
  (ii) decision gives search. So the counting route asks for something at least as strong as NAV at the threshold
  budget, never weaker.
- **(vi) Candidate counters, priced.**
  - *Monte Carlo over words:* the success probability per sample is ≈ 1/φ at the threshold, so Θ(φε^{−2}) samples
    (codex).
  - *Monte Carlo over lattice points of {(u, v) : u ≡ tv}:* the density of smooth pairs is ρ-type,
    L[2/3, −2/(3κ)] at y = L[1/3, κ] (Prop TD), so L[2/3] samples per unit of count.
  - *Characters:* C(t, B) = φ^{−1}Σ_χ χ̄(t)Ĉ_B(χ). Evaluating χ of large order needs discrete logs modulo the prime
    factors (and the factorisation), i.e. FACT ∧ DL, which already gives NAV_∞ (Prop H). For NAV_3 one would also
    need few characters to suffice. E2 found few characters suffice in the dense toy regime (0.0–0.8% of all), but the
    share rises with the factor base (4–8% at 97–309 primes, POST HOC), and bounded-height NAV is always sparse
    (Prop LT). In every measured cell the number of needed characters is a fixed fraction of φ, i.e. N^{1−o(1)}.
  - *Sieves / analytic number theory:* equidistribution theorems for smooth numbers in progressions (Granville,
    Soundararajan, Harper; GRH or averaged over moduli) give the *mean* of C. They give no pointwise counts at the
    threshold for an individual large modulus with y = L[1/3].

**Verdict on Route 1.** The reduction NAV ≤ (approximate restricted counting) is real (JVV). The counting oracle is
trivial exactly where it is useless (flat top, enumerable bottom) and equivalent to NAV decision where it is needed.
Failure point: no counter below search cost exists for the threshold window; Monte Carlo is φ or L[2/3], and characters
need DL.

---

## 4. Route 2: lattice / CVP

**Prop LT (dimension and geometry) [PROVED for (a); (b)–(c) HEURISTIC/SCENARIO as labelled].**
- **(a) Dimension floor [PROVED].** A NAV_A solver with success δ ≥ 1/poly must have
  #{e : |e| ≤ Ax} ≥ δ·|Q|, and #{e : |e| ≤ Ax} ≤ Ψ(N^A, y)·N^{o(1)} (each smooth n ≤ N^A gives 2^{ω(n)} = N^{o(1)}
  coprime splits). By de Bruijn (Ψ(X, (ln X)^c) = X^{1−1/c+o(1)} for fixed c > 1), this needs y ≥ (Ax)^{A/(A−1)−o(1)}.
  For A = 3: y ≥ (3x)^{3/2−o(1)}, so K ≥ x^{3/2−o(1)}. The number of letters in a solution is ≤ Ax/ln 2 and typically
  ≈ Ax/ln y, so K/letters → ∞: bounded-height NAV is always in the *sparse* regime.
- **(b) The CVP [PROVED reformulation].** Given a basis of Λ_N (obtainable from FACT ∧ DL, Prop H(c)) and one coset
  representative e_a (∏ℓ^{e_a} ≡ a), NAV_A is: find a vector of e_a + Λ_N with weighted ℓ1 norm |·| ≤ Ax. Under GRH
  the distance is ≤ (2 + ε)x, so NAV_3 is CVP_γ with γ = 3/(2 + ε) ≈ 3/2 in weighted ℓ1 (codex). Lemma RT gives
  distance ≥ x − 2 ln x for most a. In ℓ2 with Adleman weights √ln ℓ, squarefree vectors keep their norm², so the
  target radius is ≈ √(3x).
- **(c) Below the Gaussian heuristic [HEURISTIC/SCENARIO].** det Λ_N = |⟨FB⟩| ≤ φ, so det^{1/K} = e^{O(x/K)} → 1
  and GH ≈ √(K ln y/(2πe)). That is ≈ x^{3/4} at the minimal K and L[1/3, κ/2] at K = L[1/3, κ], against a target
  radius √(3x). The sought vectors are therefore x^{1/4} to L[1/3, κ/2] times shorter than the Gaussian-heuristic
  radius; the volume heuristic is simply wrong for these sparse integer points (siconstruct Prop CV(e);
  Mulder–Sterner–van Woerden).
  - Every known constant-γ CVP algorithm costs 2^{Θ(K)}: exact 2^{K+o(K)} (Aggarwal–Dadush–Stephens-Davidowitz
    2015), sieving 2^{0.29K}-type heuristically. With K ≥ x^{3/2}, that is ≥ 2^{x^{3/2}} ≫ N.
  - BKZ-β gives approximation exp(Θ(K ln β/β)); at an L[1/3, 1.2] budget β = O(x^{1/3}(ln x)^{2/3}), so the factor is
    exp(L[1/3]) (codex).
  - CSI-FiSh (Beullens–Kleinjung–Vercauteren 2019) solves the analogous class-group problem by BKZ + Babai in dimension
    74. That is the dense regime (det^{1/K} ≫ 1), where the Gaussian heuristic is accurate.
- **(d) Schnorr's lattice does not encode the congruence [PROVED observation; codex].** The real-log lattice (rows
  √ln ℓ·e_ℓ, last coordinate C ln ℓ) has determinant ∏√ln ℓ·(1 + C²Σln ℓ)^{1/2}, not N. A target C·ln a enforces
  u/v ≈ a over R, not u − av ∈ NZ. With a lift, the residual u − Av must be tested, which is Prop CV / schnorrcvp
  (γ = 1 + O(1/x) needed, refuted quantitatively there).
- **(e) Even with all logarithms, generalised birthday fails at bounded height [HEURISTIC, standard k-tree count].**
  Wagner's k-tree on the knapsack Σe_ℓ log ℓ ≡ log a (mod M) needs k lists of size ≥ M^{1/(1+lg k)}. A list of
  vectors of height ≤ H has ≤ H^{1+o(1)} elements, so the output height is ≥ k·x/(1 + lg k)·(1 − o(1)). Height ≤ 3x forces
  k/(1 + lg k) ≤ 3 + o(1), i.e. k ≤ 8, and cost ≥ N^{1/4−o(1)}.

**Verdict on Route 2.** Failure points: (1) no basis without FACT ∧ DL; (2) with a basis, a constant-factor CVP far
below the Gaussian heuristic in dimension ≥ x^{3/2}, which every known algorithm solves only in 2^{Θ(K)}; (3) the
archimedean lattice encodes the wrong constraint. Required approximation factor: γ ≈ 3/2 in weighted ℓ1, or in
Gaussian-heuristic units a factor x^{−1/4} to L[1/3, −κ/2]. BKZ and sieving at L[1/3] cost are exp(L[1/3]) away.

---

## 5. Route 3 (meet-in-the-middle) and Route 4 (descent)

### 5.1 Meet-in-the-middle

- **Fixed lists [PROVED; codex blind §3].** For lists S₁, S₂ fixed independently of a uniform a ∈ Q,
  Pr[aS₁ ∩ S₂ ≠ ∅] ≤ |S₁||S₂|/|Q|, since each pair (s₁, s₂) fixes one a. Lists of size L[1/3] succeed with
  probability N^{−1+o(1)}.
- **Adaptive lists** that only multiply and compare are generic: Prop G gives Ω(√(qδ)). sihard Prop MM covers interval
  tables.
- **The hidden-root promise.** By Prop SQ the solver cannot exploit it. The square subgroup has the same generic
  hardness, since q is odd, and the counts per class are the same for squares and non-squares in every E2 cell
  (|Δh| ≤ 0.012, P2d held).
- **Verdict.** Plain or structured meet-in-the-middle does not amortise; squares are not special at composite N
  beyond the factor 2^ω.

### 5.2 Descent over Z

**Prop TD (tested relations cost L[2/3] at y = L[1/3]) [MODEL: H_S for tested pieces].** A relation of height
S ≥ x − O(ln x) (Lemma RT) whose smoothness is tested on j ≤ x^{o(1)} pieces is y-smooth with probability
≤ exp(−(S/ln y)·ln(S/(j ln y))(1 + o(1))) = L[2/3, −2/(3κ)] at y = L[1/3, κ]. At κ = 0.6 that is L[2/3, 1.11] per
output (codex: the same constant for Gauss reconstruction plus testing). Pieces of size ≈ y (j ≈ x^{2/3}) are
construction, not testing.

**Z-descent [PROVED identity + MODEL].** A descent step for a prime q is a relation u ≡ qv (mod N) with u ≠ qv, so
|u| + q|v| ≥ N (smoothbuild Lemma DSC), and its tested pieces have log-size ≥ x − ln q − O(1). For the leaf level
(children ≤ y = L[1/3]) to cost L[1/3], the pieces must have log-size O(x^{2/3}(ln x)^{1/3}) (Theorem SS), which
forces q ≥ N/L[2/3] and qv within L[2/3] of a multiple of N. That is suffmap's Prop SI window problem. Its
constructors are open, and a modulus-generic one would improve prime-field DL (sihard Prop DLI). Recursion does not
relax the relative precision (smoothbuild §5).

### 5.3 NFS descent

- **Heuristic costs [KNOWN].** NFS individual-logarithm descent: Joux–Lercier 2003; Commeine–Semaev 2006
  (L_p[1/3, 3^{1/3} ≈ 1.442] for the initial split); Barbulescu's early-abort analysis ≈ 1.232; Fried–Gaudry–Heninger–
  Thomé 2017, Appendix A, recompute these and label them heuristic (codex). 1.232 already exceeds NAV's 1.20 budget.
- **Shape.** Descent outputs live in ideals of a number field: t ≡ (rational smooth) × φ(product of factor-base
  ideals' generators), usable through virtual logarithms. Turning them into Z-relations needs an exact integer kernel
  vector of the ideal-exponent matrix, with the unit part made trivial. Such vectors exist (rank), but their entries
  are generically huge (Hadamard-size, L[1/3] digits). So NFS gives NAV_∞-type outputs heuristically, which is
  NFS-DL by another name (Prop H(c)), and not NAV_3. *Remark [MODEL]:* a bounded-height conversion needs O(d)
  relations whose algebraic factorisations cancel exactly. Under independent incidence over L[1/3] ideals, each relation
  carrying ≈ x^{1/3} ideals, the expected number of such combinations → 0.
- **Obs DD (what rigorous L[1/3] factoring actually needs).** paper2's Theorem DD (Rabin/Dixon fibre symmetry
  applied to a descent; internal, unrefereed, checked by AI referees) shows: if each row is a descent of a uniform
  square r_i² into an NFS factor base, then every dependency is fruitful with probability 1 − 2^{1−ω}, whatever the
  rows' algebraic shape. The hidden root survives a descent into ideals. (The combination with Lee–Venkatesan's rows
  below is a sketch: we have not checked it against their square-root and character framework.) So the hidden-root ingredient that NAV
  packages is already available for NFS-shaped rows. The Z-shape and the height bound are unnecessary. With
  standard NFS rows (Lee–Venkatesan 2018 prove relation collection and a possibly trivial congruence rigorously at
  L[1/3, 1.923]), it suffices that the hidden-root rows span the corank of the standard rows. In the worst case that
  is K + 1 descents, still L[1/3] if each descent is.
- **Where the rigour stops.** lv5unif §13 (and codex blind §4): Theorem U's one-level descent needs a split acceptance
  (ln N)^{−2k} affordable (k ≤ D^α) and a norm exponent θ ≥ c/√k small (k ≳ D^{2(1−2α)}), with D = x/ln x, which
  forces α ≥ 2/5. Exponent 1/3 needs multi-level descent with control of the distribution of the target-dependent
  children (second-level special-q primes are not uniform in an interval). That is a different, narrower open problem
  than NAV, and it is the one we recommend.

**Verdict on Route 4.** Z-descent fails at the leaf level (Prop TD; Prop SI windows). NFS descent is heuristically
L[1/3] but ideal-shaped. NAV's Z-shape is not needed (Obs DD), and the rigour gap is multi-level descent geometry.

---

## 6. Answers to the five questions

1. *Is approximate counting easier than finding?* No method is known. Counting is trivial where uninformative and
   equivalent to NAV decision at the threshold (Prop CF). Characters need DL, and at toy size they need a constant
   fraction of all characters. Cost of the counting oracle: Θ(φ) (word Monte Carlo) or L[2/3] (lattice Monte Carlo).
2. *Approximation factor?* γ ≈ 3/2 in weighted ℓ1 relative to the true distance, in dimension ≥ x^{3/2}, which is
   x^{1/4} to L[1/3, κ/2] below the Gaussian heuristic, and only once a basis of Λ_N (FACT ∧ DL) is in hand. Outside
   BKZ's reach; sieving costs 2^{Θ(K)} (Prop LT).
3. *Do squares help MITM?* No: the solver cannot see the root, and square- and unit-NAV are equivalent up to 2^ω
   (Prop SQ). MITM is generic (Prop G).
4. *Descent?* Tested Z-descent is L[2/3, 2/(3κ)]; integer descent needs Prop SI windows. NFS descent is heuristically
   L[1/3, 1.23–1.44] per target but ideal-shaped. Rigorous NFS-type descent for random squares is exactly the project's
   Dixon descent, rigorous only at one level (L[2/5]).
5. *Barrier?* NAV ⇒ FACT ∧ F_p-DL at L[1/3, max(κ + τ, 3κ)] (Theorem DLB); NAV_∞ ≡ FACT ∧ DL; generic Ω(√q). It is not
   equivalent to rigorous L[1/3] factoring alone (it also gives DL, unless DL ≤ FACT), and NAV_3 adds a CVP layer.

---

## 7. Experiments (pre-registered in `notes/navigate_prereg.md`, commit `f2e4e9c`)

### 7.1 E1: toy walks vs sampling (operational)

Balanced semiprimes, 6 per size, 100 targets per N, 400 smoothness tests per (algorithm, target), y = 2^{bits/4}.

| bits | Y_S0 | Y_G1 | Y_J2 | J2/S0 | G1/S0 | R_S0 | R_G1 | mean tested log-size S0 → G1 | G1/S0 distinct | splits/output S0, G1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 0.0721 | 0.1268 | 0.0725 | **1.005** | **1.759** | **1.013** | **1.009** | 25.2 → 21.9 | 1.087 | 10.5, 84.8 |
| 48 | 0.0730 | 0.1151 | 0.0738 | **1.011** | **1.577** | **1.005** | **1.012** | 30.7 → 27.6 | 1.306 | 10.4, 96.6 |
| 56 | 0.0725 | 0.1068 | 0.0720 | **0.994** | **1.473** | **1.002** | **1.016** | 36.3 → 33.2 | 1.382 | 10.5, 106.6 |

- P1a (R_S0 ∈ [0.85, 1.15]) **held**; P1b (J2/S0 ∈ [0.85, 1.15]) **held**; P1c (G1/S0 ≥ 1.15) **held**; P1d
  (R_G1 ∈ [0.85, 1.15]) **held**.
- *Reading.*
  - Monte Carlo completion counting (J2) is sampling with extra bookkeeping.
  - The height-steered walk (G1) wins per smoothness test only because it tests smaller pieces (≈ 3.1–3.3 nats
    smaller): the size-matched coprime null explains its rate to within 1–2%. It pays ≈ 8–10× more splits per
    output.
  - Per *distinct* relation its gain falls to 1.09–1.38. Post hoc we attribute this to the steered walk lingering
    near small-height states, whose splits are rescalings of one relation (smoothbuild E2's effect).
  - Signs: 50% of outputs have u, v > 0 for all three algorithms (descriptive).
- No L-scale conclusion is drawn (u = 2 at every size).

### 7.2 E2: exact completion counts (diagnostic)

20-bit semiprimes (4 seeds); every signed exponent vector of height ≤ β ln N over FB = primes ≤ y; per-class counts
C_β on H = ⟨FB⟩ (H = G in all 8 instances).

| y (K) | β | μ (mean) | h (hit frac.) | Poisson 1 − e^{−μ} | D = Var/μ | f75 |
|---|---|---|---|---|---|---|
| 13 (6) | 1.0 | 0.06 | 0.053–0.067 | 0.057–0.065 | 0.93–1.15 | 0.002–0.005 |
| 13 | 1.5 | 0.59–0.66 | 0.35–0.54 | 0.44–0.49 | 0.71–1.54 | 0.000–0.002 |
| 13 | 2.0 | 3.1–3.5 | 0.88–1.00 | 0.96–0.97 | 0.47–1.69 | 0.000–0.002 |
| 13 | 3.0 | 34–39 | 1.00 | 1.00 | 0.30–1.17 | 0.000–0.001 |
| 19 (8) | 1.25 | 0.73–0.81 | 0.47–0.58 | 0.52–0.55 | 0.71–1.22 | 0.001–0.006 |
| 19 | 2.0 | 23–26 | 1.00 | 1.00 | 0.42–0.67 | 0.001–0.003 |

(full grid in the summary JSON.)

- **P2a failed** (11 cells, all with μ ∈ [0.5, 3.5]): hit fractions deviate from Poisson by up to 0.15 (seed 2,
  y = 13, β = 1.75: 0.614 vs 0.764), in both directions, seed-dependently.
- **P2b failed** (25 cells): D ranges 0.30–1.69. Seeds 0 and 3 are strongly under-dispersed (counts close to
  constant at large β); seeds 1 and 2 are over-dispersed.
- **P2c failed** (28 cells): f75 = 0.000–0.006 in the registered cells (≤ 0.008 in all cells), against the registered
  [0.2, 0.5] (noise ≈ 0.39). The counts are *not* noise-like: 75% of the non-principal Fourier mass sits on ≤ 0.6% of
  the characters.
- **P2d held**: squares and non-squares have the same hit fraction to within 0.012 in every cell.
- *Explanation (POST HOC; computation below).* With K = 6–8 primes and budgets of 12+ letters, each prime is used
  with many exponents: a *dense* regime. Ĉ_B(χ) = Σ_{|e|≤B}∏χ(ℓ)^{e_ℓ} is then a K-fold Fejér-type kernel in the angles
  arg χ(ℓ), of width ≈ K ln ℓ/B per prime (typical |e_ℓ| ≈ B/(K ln ℓ)). Only characters with all K angles small
  carry mass, a fraction ≈ ∏_ℓ(K ln ℓ/B) ≈ 0.008 at K = 6, B = 21 (an order-of-magnitude heuristic; observed
  0.000–0.002 there). Equivalently, C_β is a smooth function of the coset of the
  low-dimensional lattice Λ_N, whose geometry (short relation vectors or not) also sets the seed-dependent
  dispersion. This is the CSI-FiSh regime, the opposite of bounded-height NAV (Prop LT(a)).
- *Post-hoc sparse regime* (`posthoc_e2.py`, same code; y = 512 (K = 97), β ∈ {0.8, …, 1.1}; y = 2048 (K = 309),
  β ∈ {0.8, 0.9, 1.0}):

  | y (K) | β | μ | h | Poisson | D | f75 |
  |---|---|---|---|---|---|---|
  | 512 (97) | 0.8 | 0.32–0.34 | 0.29–0.31 | 0.28–0.29 | 0.85–0.86 | 0.056–0.059 |
  | 512 | 0.9 | 1.16–1.20 | 0.74–0.76 | 0.69–0.70 | 0.73–0.74 | 0.055–0.058 |
  | 512 | 1.0 | 4.0–4.1 | 0.99 | 0.98 | 0.69–0.70 | 0.048–0.051 |
  | 2048 (309) | 0.8 | 0.40–0.42 | 0.36–0.37 | 0.33–0.34 | 0.82–0.84 | 0.074–0.078 |
  | 2048 | 0.9 | 1.57–1.60 | 0.85–0.86 | 0.79–0.80 | 0.67–0.68 | 0.077–0.081 |
  | 2048 | 1.0 | 5.9–6.0 | 1.00 | 1.00 | 0.64–0.65 | 0.070–0.073 |

  - The Fourier share rises 30–80-fold (0.001 → 0.05 → 0.075) as K grows. It is still below noise (0.39), and the
    seed-to-seed spread vanishes.
  - Counts are mildly sub-Poissonian (D 0.64–0.86), with hit fractions 0.02–0.06 above Poisson.
  - *Reading:* the low-frequency structure is a dense-regime artefact. Toy sizes cannot reach the regime of NAV at
    L[1/3] (K = L[1/3] ≫ letters ≈ x^{2/3}), so the extrapolation is CONJECTURED, not shown.
  - Resource deviation: this post-hoc run peaked at 4.9 GB RSS for ≈ 1 min (head enumeration over 95 primes), above
    the 2 GB guideline; no swapping occurred (memory 84% free before). Not re-run.
- *Oracle-walk profile (descriptive).* Along uniformly sampled solutions, the fraction σ of next letters that keep a
  non-zero completion count is ≈ 1 while ≥ 2 ln N of budget remains (y = 13). It falls to 0.2–0.5 below ≈ 1.6 ln N.
  In the sparse post-hoc cells it is 0.03–0.1 near the threshold. The information a navigator must supply sits at the
  end of the walk, as Prop CF predicts.

**Deviations.**
- (1) A smoke test of `core.py` used seed 99 at 40 bits (one target; not a registered instance).
- (2) The post-hoc run above, and its memory overrun.
- (3) Signs are allowed in outputs (−1 column), as registered; the target in smoothbuild §7 says "positive". About
  half the outputs are positive.

---

## 8. Literature notes (from memory unless marked; to be checked by the literature critic)

- **JVV / counting.** Jerrum–Valiant–Vazirani, TCS 43 (1986): approximate counting ⇔ almost-uniform generation for
  self-reducible relations; Sinclair–Jerrum (1989) bootstrapping.
- **Generic lower bounds.** Nechaev 1994; Shoup, EUROCRYPT '97.
- **Schnorr '91** (EUROCRYPT '91): factoring/DL via Diophantine approximation in the prime-number lattice. Relation
  shape u ≡ v·(stuff); refuted quantitatively in schnorrcvp.
- **JMV 2009**: GRH expanders on (Z/N)^× with small-prime generators.
- **Rigorous index calculus.** Pomerance 1987 (L[1/2, √2] factoring and DL in F_p, F_{2^n}); Lenstra–Pomerance 1992
  (L[1/2, 1] factoring); Enge–Gaudry 2002 (rigorous L[1/2, √2] DL in hyperelliptic Jacobians of large genus); the
  L[1/3] curve results of Enge–Gaudry 2007 / Enge–Gaudry–Thomé 2011 are heuristic (from memory).
- **NFS descent.** Joux–Lercier 2003; Commeine–Semaev 2006; Barbulescu (thesis 2013) early abort; Fried–Gaudry–
  Heninger–Thomé (EUROCRYPT 2017, arXiv 1610.02874) Appendix A.
- **Rigorous NFS.** Lee–Venkatesan, J. Number Theory 2018 (as summarised in lv5unif §14).
- **Relation lattices modulo composites.** Stange, arXiv:2211.06821 (2022): factoring from multiplicative relations
  modulo n (our Λ_N), under her Hypothesis 3.1; paper2 cites her as the closest precedent for hidden random targets.
- **Class-group analogue.** Beullens–Kleinjung–Vercauteren, CSI-FiSh, ASIACRYPT 2019: short representations in a
  class-group relation lattice of dimension 74 by BKZ + Babai.
- **Generalised birthday.** Wagner, CRYPTO 2002.
- **Smooth numbers.** de Bruijn 1966 (Ψ(x, (log x)^c)); Granville; Soundararajan 2008 and Harper 2012 (smooth numbers in
  progressions, GRH).
- **Bach 1984**: DL modulo N ≡ FACT ∧ DL modulo primes.

---

## 9. GPT-6.1-sol rounds

### 9.1 Blind round (`codex/navigate_1.out`, launched before any analysis here)

Codex agreed on every route's failure point, and several statements here come from it:
- the JVV restricted-count caveat;
- word Monte Carlo at Θ(φ) samples;
- the GRH derivation of the constant 3;
- γ ≈ 3/2 and the Gaussian-heuristic discrepancy L[1/3, 0.30];
- Schnorr's lattice not encoding the congruence;
- the fixed-list MITM bound;
- reconstruction + testing at L[2/3, 1.111];
- the descent constants 1.442/1.232 (FGHT Appendix A), and that 1.232 > 1.20;
- the padding correction (r > y; a fixed small prime is inadmissible) and the conservative L[1/3, 1.80] for EQ3;
- the Theorem U geometry k ≤ D^α vs k ≳ D^{2(1−2α)}.

Its probability: ≈ 5% (any constant), ≈ 1% at (0.60, 1.20). Its most promising route, "target-conditioned multilevel
descent extending Theorem U", is our Obs DD recommendation. It left "a bounded-height reduction from arbitrary-target
NAV to square-target NAV" open; Prop SQ closes it. It did not note that EQ1 tolerates unbounded heights (Prop H).

### 9.2 Adversarial round

(pending)

## 10. Opus critics

(pending)

## 11. Reproduction

```
.venv/bin/python -m src.navigate.e1_walks      # ~25 s
.venv/bin/python -m src.navigate.e2_counts     # ~20 s, peak RSS 1.3 GB
.venv/bin/python -m src.navigate.posthoc_e2    # ~60 s, peak RSS 4.9 GB (post hoc)
.venv/bin/python -m src.navigate.analyze
.venv/bin/python -m pytest -q tests/test_navigate.py
```
