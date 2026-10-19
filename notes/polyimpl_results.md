# polyimpl: what a fast factoring algorithm would imply, and what that says about blockers and routes

Agent `polyimpl`, phase 7, 2026-10-19.
- Pre-registration: `notes/polyimpl_prereg.md` (commit `1a9c3de`, made before E1 ran).
- Code: `src/polyimpl/{residual,primefind}.py`. Tests: `tests/test_polyimpl.py`.
- Raw data: `experiments/raw/polyimpl_e1.jsonl`. Summary:
  `experiments/analysis/polyimpl_e1_summary.json`.
- GPT-6.1-sol: blind round `codex/polyimpl_1.{md,out}`; adversarial review
  `codex/polyimpl_2.{md,out}` (§9).
- Independent Opus critic: §9.
- Status labels: **KNOWN** (with reference), **PROVED** (proof here; routine unless said
  otherwise), **CONDITIONAL**, **SPECULATIVE**, **FALSE**, **experimentally supported**.
- "Codex §x" refers to sections of `codex/polyimpl_1.out`. Its proofs were re-checked here.

## 0. Verdict (plain language)

**Fast factoring would change almost nothing in mathematics, and that is the main finding.**

What each hypothesis would and would not give:
- Polynomial-time factoring, deterministic or randomized, would give:
  - fast versions of every factorization-derived function: μ, φ, σ, squarefree parts, rings of
    integers, order finding, point counts mod N, τ(N);
  - the death of RSA, Rabin and BBS.
- It would not give:
  - class groups, prime-field discrete logs, or deterministic square roots mod p;
  - any new theorem about primes;
  - a contradiction with any unconditional theorem or conjecture we checked, including
    Sarnak's conjecture and the τ-conjecture.
- **No free lower bound for factoring comes from mathematics.**

**No blocker from prime construction.** The natural candidate fails on arithmetic (Lemma PF,
§3).
- Every known "large prime factor in a short interval" theorem is strong enough that the large
  prime factor can be stripped out *without factoring anything*.
- So a polynomial-time, or N^{1/6−δ}, factoring algorithm improves no deterministic
  prime-finding bound. The barrier there is purely analytic.

**One conditional upper bound that changes the framing.**
- Theorem DR: under the standard derandomization hypothesis (E needs exponential-size circuits,
  Impagliazzo–Wigderson 1997), Lenstra–Pomerance's rigorous randomized algorithm becomes a
  **deterministic L_N[1/2, O(1)] factoring algorithm for every N**.
- So the project's "deterministic" frontier (N^{1/5} worst case, N^{1/6} average) is the
  frontier of *unconditional derandomization*, not of factoring.
- Conversely, an unconditional N^{Ω(1)} lower bound for deterministic factoring would prove
  that E has subexponential circuits infinitely often, which is believed false.

**Escape signature.**
- Every model floor in the repository and in Ward–Maestro has a known algorithm that escapes
  it (§5, table).
- What has no witness is the *conjunction*: classical + deterministic + rigorous + almost all
  (or worst case) + below N^{1/6}.
- Two new combined barriers pin down the shape such an algorithm must have:
  - **Theorem CB** (ours; conditional on a Fouvry-type sieve bound). Local collision search
    plus *any* smooth-exponent annihilation in (Z/N)^* fails on a positive proportion of
    balanced semiprimes below N^{1/6}/polylog. Smoothness bounds up to N^{1/3} and polylog
    many bases are allowed.
  - **Theorem CS** (codex, re-checked). Local search plus a fixed list of tests with R prime
    signatures needs N^{1/6}R^{−2/3}/polylog in the worst case.
- So a sub-1/6 deterministic algorithm must do one of:
  - take order information from groups other than (Z/N)^*, such as tori, curves or class
    groups, *with a hitting theorem*;
  - use relations with a deterministic smooth supply and fruitfulness;
  - use a non-local selector.
- Known families supply each ingredient separately, but none supplies the conjunction
  unconditionally.

**Experiment E1 (pre-registered; verdict "mixed" under the frozen rule).**
- Cyclotomic exposure tests d ∈ {1, 2, 3, 4, 6} at the N^{1/6} budget leave a residual of
  primes:
  - 8.3% at 20 bits, 10.2% at 26 bits, 10.9% at 32 bits;
  - pairs: 0.69%, 1.08%, 1.13%.
- The residual *grows* toward the independence-model value 14.1%.
- Three of the five frozen ranges were missed on the low side at 20 bits, and the stability
  ratio was 1.31, just above its band, because of growth, not decay.
- Read post hoc: no decay, so a cyclotomic exposure portfolio cannot supply the polynomially
  small residual that codex's Rank-4 route needs (§4.4, §7).

**Ranked directions (§7).**
1. Fixed-u deterministic quadratic sieve on almost all N. Below 1/6 needs only *constant*
   smoothness density, plus fruitfulness and deterministic elimination.
2. Group-order hitting theorems outside (Z/N)^*. Theorem CB says this is forced; Prop LW says
   the worst-case version is a least-witness problem.
3. Non-local selectors with a key-size certificate (codex Rank 1).
4. Promise-preserving partial-information reductions (codex Rank 2).

None is cheap *and* decisive beyond E1.

---

## 1. Hypotheses and how consequences transfer

The hypotheses:
- (H1) deterministic polynomial time;
- (H2) randomized polynomial time (classical);
- (H3) L[1/3 − δ] or L[1/4], heuristic or rigorous;
- (H4) either of:
  - (H4w) deterministic worst case N^{1/5−δ};
  - (H4a) deterministic average or almost-all N^{1/6−δ} on balanced semiprimes;
- (H5) partial information only.

**PROVED (transfer rule; codex §2).** A reduction that makes poly(n) factoring calls on integers
≤ N^c turns an N^α factorizer into N^{cα+o(1)}, and an L[a, C] factorizer into L[a, Cc^a].
- Average-case hypotheses (H4a) do **not** transfer to reductions that query integers outside
  the ensemble (shifted values, p − 1, short intervals). An exceptional set of density
  N^{−Ω(1)} can contain any prescribed short interval.

**PROVED (no collapse).** Factoring is a total search problem with a unique, efficiently
verifiable answer.
- It lies in NP ∩ coNP (Pratt; AKS). It lies in UP ∩ coUP (Fellows–Koblitz 1992; citation from memory, UNVERIFIED). It lies in
  BQP (Shor 1997).
- So H1 collapses no complexity class and needs no NP-hardness to fail. H1 does not imply
  P = BPP.

**Theorem DR (PROVED; a routine composition of KNOWN results; no explicit statement found by us).**
- *Hypothesis.* Some language in E = DTIME(2^{O(m)}) has circuit complexity ≥ 2^{εm} for all
  large m (Impagliazzo–Wigderson, STOC 1997, main theorem).
- *Conclusion.* There is a deterministic algorithm that outputs the complete factorization of
  every N in time L_N[1/2, c] for an absolute constant c.
- *Corollary.* Under the same hypothesis, H2 ⇒ H1.

*Proof.*
- Lenstra–Pomerance (JAMS 5 (1992), main theorem, §1) is a Las Vegas algorithm with expected
  time T = L_N[1/2, 1+o(1)], with no GRH.
- Truncate it at 2T. By Markov it outputs a full factorization with probability ≥ 1/2, using at
  most 2T random bits.
- The output is checked (multiply, then AKS), so the predicate "seed r gives a verified
  factorization" is computed by a circuit D_N of size S = T^{O(1)}, with N hardwired.
- The IW97 hardness–randomness trade-off gives a generator G : {0,1}^{O(log S)} → {0,1}^{S}
  that fools every size-S circuit with error < 1/4. It is computable in time poly(S), because
  the hard function at input length O(log S) is evaluable in time 2^{O(log S)}.
- Some seed s therefore has D_N(G(s)) = 1. Enumerate all S^{O(1)} seeds.
- Total time S^{O(1)} = L_N[1/2, O(1)]. The H2 ⇒ H1 corollary is the same argument with T
  polynomial. ∎

*Contrapositive.* Suppose deterministic factoring needs time > L_N[1/2, C] for every C on
infinitely many n.
- Then for every ε > 0, every language in E has circuits of size < 2^{εm} for infinitely many
  m: an "i.o. subexponential circuits for E" upper bound, believed false.
- So a superpolynomial-in-L lower bound for deterministic factoring is out of reach for the
  same reason circuit upper bounds for E are.

Relation to what was known:
- below16 §6.2(3) and Ward–Maestro §35 already note the *non-uniform* deterministic L[1/2]
  (Adleman's argument). Theorem DR is the uniform, conditional version.
- The prime-construction analogue (P = BPP ⇒ primes in polynomial time) is noted in Polymath4.
- Factoring has a unique output, so "pseudodeterministic" factoring is just Lenstra–Pomerance;
  compare Chen–Lu–Oliveira–Ren–Santhanam 2023 for primes.

---

## 2. Consequence catalog by hypothesis

### 2.1 (H1) deterministic polynomial time

| Area | Consequence | Label / reference |
|---|---|---|
| Arithmetic functions | μ, λ (Liouville), ω, Ω, d, φ, σ_k, λ (Carmichael), rad, squarefree part and decomposition are all in FP | PROVED (formulas; codex §3.1) |
| Reverse reductions | σ ⇒ factoring, randomized | KNOWN: Bach–Miller–Shallit, SIAM J. Comput. 15 (1986), Thm 2 |
| | φ, λ ⇒ factoring, randomized; deterministic under ERH (Miller 1976); deterministic for N = pq (s = N + 1 − φ) | KNOWN; PROVED for N = pq |
| | O(log N) φ-queries ⇒ deterministic full factorization in L[1/3, 1+o(1)] | KNOWN: Żrałek, Math. Comp. 79 (2010) 513–533 |
| | d(N) ⇒ μ with one call, d(N^{q−1}) | KNOWN: Shallit–Shamir, IPL 20 (1985), Thm 1 |
| | d ⇔ exponent multiset | KNOWN: Shallit–Shamir, Thm 2 |
| | factoring ≤ μ, ≤ squarefreeness, ≤ ω, ≤ d: **no reduction known** | KNOWN-open |
| Rings of integers | maximal order ⇔ largest squarefree divisor of the discriminant, so H1 ⇒ O_K in deterministic polynomial time | KNOWN: Buchmann–Lenstra, JTNB 6 (1994), Thms 1.3, 1.4, 6.12; Chistov 1989 |
| Class groups, units | **not** implied: relation finding, regulators and principal-ideal tests remain L[1/2] (Hafner–McCurley 1989, GRH). H1 does give the 2-rank (genus theory: ω(D) − 1) | KNOWN / SPECULATIVE for "not implied" |
| QR, square roots | QR decision mod N is in P | PROVED (Euler + CRT) |
| | square-root *construction* is in randomized polynomial time only; deterministic needs a non-residue (GRH: Ankeny–Bach) or a fixed radicand (Schoof 1985). H1 does not derandomize Tonelli–Shanks | PROVED / KNOWN |
| | Rabin: square roots ⇒ factoring, randomized | KNOWN |
| Order finding | ord_N(a) in FP (factor N, then each p − 1) | PROVED (codex §3.2) |
| | order ⇒ factoring, randomized (Miller/Shor) | KNOWN |
| Discrete logarithms | DLP mod N ⇒ factoring; DLP mod N ⇐ factoring + DLP mod p. H1 leaves prime-field DLP at L[1/3] | KNOWN: Bach, UCB/CSD-84-186 (1984) |
| r₂(N), two squares | r₂ in FP; one representation constructed deterministically (Schoof √−1 mod p + Cornacchia) | PROVED + KNOWN |
| Point counts | #E(Z/N) in FP (CRT + Schoof) | PROVED |
| | for N = pq, #E and #E′ (twist) ⇒ deterministic factoring | KNOWN: Dieulefait–Urroz, arXiv 1911.11004, Thm 1 |
| | Kunihiro–Koyama 1998: probabilistic, with a gap noted by Dieulefait–Urroz | KNOWN |
| Modular-form coefficients | τ(N) in FP (Edixhoven–Couveignes + multiplicativity) | KNOWN / PROVED |
| | computing τ(N) for RSA moduli is factoring-hard | UNVERIFIED citation: Bach–Charles 2007, "The hardness of computing an eigenform" |
| Explicit constructions | primitive roots mod p become *testable* in P, but finding one is still p^{1/4+ε} deterministic (Burgess / Shparlinski 1996), or polylog under GRH (Shoup) | KNOWN |
| | no change to irreducible-polynomial construction | SPECULATIVE |
| Prime construction | **no new bound** (§3) | PROVED for all known short-interval inputs |
| Complexity | see §4 | |
| Analytic NT | **nothing**. Faster evaluation of μ(n), Λ(n) or [P⁺(n) ≤ y] gives no cancellation in sums. Mertens sums are already sub-x^{2/3} without factoring (Deléglise–Rivat x^{2/3}; Helfgott–Thompson 2023 ≈ x^{3/5}, exponent UNVERIFIED by us). Only finite verification gets cheaper: odd-perfect bounds, aliquot sequences, squarefreeness of the Fermat numbers F_m for small m, Cunningham tables | PROVED (codex §6) + SPECULATIVE (none found) |
| Crypto | RSA, Rabin, BBS, Paillier, Goldwasser–Micali, RSA-group VDFs and accumulators die. Class-group, lattice, hash and code-based schemes survive (H1 does not compute class groups) | KNOWN |

### 2.2 (H2) randomized polynomial time

- Everything in §2.1 holds with "randomized".
- Under the IW97 hypothesis, H2 ⇒ H1 (Theorem DR). Without it, nothing about derandomization
  follows.
- Natural proofs:
  - Naor–Reingold's factoring-based TC⁰ pseudorandom function (JACM 51 (2004)) loses its
    hypothesis.
  - The TC⁰ natural-proofs barrier survives through DDH or lattice PRFs (Banerjee–Peikert–Rosen
    2012).
- Learning hardness: Kearns–Valiant 1994 and Kharitonov 1993 lose their factoring instances.
  Lattice-based versions (Klivans–Sherstov 2009) survive. **KNOWN.**
- Quantum evidence: factoring stops being evidence that BQP ⊄ BPP. Simon, Forrelation and
  sampling-type evidence remain. KNOWN (folklore).

### 2.3 (H3) L[1/3 − δ] or L[1/4]

- There are no mathematical consequences.
- Complexity: assumptions of the form "factoring is hard for circuits of size 2^{n^ε}"
  (Kharitonov; Naor–Reingold natural proofs) already need ε < 1/3 because of NFS. H3 would
  force ε < 1/3 − δ, or ε < 1/4.
- The only historical analogue of an L[1/3] → L[1/4] jump is discrete logs in small
  characteristic (Joux 2013; Barbulescu–Gaudry–Joux–Thomé 2014). It used the Frobenius identity
  x^q − x = ∏(x − a) of a subfield, which has no analogue over Z.
  - SPECULATIVE: an L[1/4] factoring algorithm would need a "Frobenius substitute" for Z.
  - Rigorous L[1/3] would be new for factoring *and* for prime-field DLP. The repository's
    barriers are nfsfruit (Lemma M), lemmaM (c) (HECM splitting forces L[4/9] in that
    architecture) and hypD (Prop H).

### 2.4 (H4) the project's own frontiers

**(H4w) Deterministic worst case N^{1/5−δ}.**
- It gives deterministic N^{1/5−δ} versions of:
  - φ, λ, μ, σ, the squarefree part;
  - order finding mod N (p − 1 adds only N^{1/10});
  - strong-root extraction in (Z/N)^*.
- These already beat the generic-group root-extraction floor (Damgård–Koprowski 2002) at
  exponent 1/5, through Harvey–Hittmeir: a witness that representation-dependent arithmetic
  escapes generic floors.
- It gives no prime-construction gain (§3) and no other mathematical consequence.
- Inside the project, it must violate the hypotheses of:
  - Ward–Maestro Thm 3.4 and Thm 9.4;
  - cycloharvey LB;
  - multimod Thm A′/B (it needs N^{Ω(1)} moduli, or must leave the model);
  - coverage38 Theorem U (it needs coverage of *every* pair with λ < J^{1/4}).
- Through Theorem R, it would prove a covering statement for which no proof method is known
  (dchworst B1–B3).

**(H4a) Average or almost-all N^{1/6−δ} deterministic.**
- It must leave:
  - T1/T2/T3 (below16);
  - Theorem CB and Theorem CS (§4.3);
  - multimod Thm A (any full-candidate cover).
- It implies nothing outside factoring, because transfer fails off-ensemble (§1).

### 2.5 (H5) partial information

| Oracle | What it gives on S_n (balanced semiprimes) | Label |
|---|---|---|
| ω, Ω, d, μ, μ², rad on N | constant (2, 2, 4, 1, 1, N): **0 bits** | PROVED |
| the same on coprime masks Nk, Nk² | determined by k | PROVED (codex §9.3) |
| "N is a sum of two squares", r₂(N) | for N ≡ 1 (mod 4): ε = (−1/p), 1 bit. For N ≡ 3 (mod 4): constant | PROVED (codex §3.3) |
| | equivalent to φ(N) mod 8 | PROVED (codex §9.2) |
| | no reduction from factoring is known (Seurin 2-Φ/4-Hiding; codex r6) | KNOWN-open |
| φ(N) mod m, all prime m ≤ Cn | φ(N) by CRT, so deterministic factoring for pq | PROVED (codex §9.2) |
| | general N: randomized (BMS) or deterministic L[1/3] (Żrałek) | KNOWN |
| φ(N) mod 2^k, k ≥ n/4 + v₂(q − p) | Coppersmith, 2 decodes | PROVED (Track E P-6) |
| φ(N) mod one fixed m | s mod m only. Information floor: per-modulus symmetric data beyond N mod m needs consistent labels | PROVED (Track E F1/F2) |
| Jacobi-type bits from N mod m | 0 bits (Track E F1). The eps battery found no predictor of ε | PROVED / experimentally supported |
| μ in TC⁰ | refutes Kalai's conjecture that μ is asymptotically orthogonal to TC⁰ (blog post, 15 Dec 2024). So factoring in TC⁰ would refute it. The conjecture implies factoring ∉ TC⁰ | PROVED (μ from a factorization is in TC⁰: equality tests and counting) |

The lesson of H5: one-bit and multiplicative-function oracles on S_n are either constant or
are ε. The only small-output oracles that factor are **uniform families** of residues of a
symmetric function (φ mod m for many m), which together encode s.

---

## 3. Deterministic prime construction (task 1(b))

**KNOWN.** Polymath (Tao–Croot–Helfgott), Math. Comp. 81 (2012):
- the best deterministic construction of a prime ≥ z costs z^{1/2+o(1)} (Lagarias–Odlyzko);
- Theorem 1.2 computes the *parity* of π on intervals of length z^{1/2+c} in z^{1/2−c+o(1)};
- Conjecture 1.1 is the existence version.
- The paper says nothing about a factoring oracle.
- The Polymath4 research thread II (9 Aug 2009), "Strategy 3", assumes an O(1) factoring
  oracle. It asks for intervals [n, n + a] guaranteed to contain a non-S-smooth number, and
  notes that S > a² would beat the RH bound.
- **No bound obtained with the oracle is stated anywhere.**

**Lemma PF (PROVED; code `src/polyimpl/primefind.py`; codex §4.2 for the first route).**
- *Input theorem T(β, θ).* Every interval [X, X + X^β] contains an n with P⁺(n) > X^θ.
- *Factoring cost.* Deterministic factoring of n ≤ 2X costs X^{α+o(1)}.
- *Conclusion.* A prime p > z = X^θ (p ≤ 2X) is found in time z^{e+o(1)} with e = min(e_fac, e_strip):
  - e_fac = (β + α)/θ: factor every integer of the interval;
  - e_strip = max(β, 1 − θ)/θ, available when θ > 1/2. Remove all prime factors ≤ 2X^{1−θ}
    from all integers of the interval with one product/remainder tree, at cost
    (X^β + X^{1−θ})^{1+o(1)}. When P⁺(n) > X^θ > 2X^{1−θ}, the cofactor n/P⁺(n) < 2X^{1−θ} is
    removed completely and the remainder is the prime P⁺(n). This uses **no factoring algorithm**.

*Proof.* Count the operations; the prime is verified by AKS. ∎

**Corollary PF.**
- (i) A recipe of this kind beats 1/2 for some α ≥ 0 only if θ > 2β, i.e. intervals of length h
  contain a prime factor > h^{2+δ}.
- (ii) The factoring hypothesis matters only when also e_strip ≥ 1/2. That means θ ≤ 1/2, or
  max(β, 1 − θ) ≥ θ/2. With (i) this forces β < 1/3 and 2β < θ ≤ 2/3.
- (iii) For H4w at α = 1/5 − δ, one also needs θ > 2β + 2/5 − 2δ, so β < 2/15 + δ.

*Proof.* Each claim is algebra on the exponents above. ∎

**Known inputs** (literature agent and codex §4.3, read from the sources; exact arithmetic in
`primefind.py`):

| Theorem (β, θ) | best e (any α) | does factoring matter? |
|---|---|---|
| Baker–Harman–Pintz 2001: a prime in X^{0.525} (21/40, 1) | 0.525 (primality test only) | no |
| Merikoski 2021, MPCPS 170, Thm 1: x^{1/2}log^{1.39}x, P > x^{18/19} (1/2, 18/19) | 19/36 = 0.528 | no (strip) |
| Jia–Liu 2000, Acta Arith. 95: x^{1/2+ε}, P > x^{25/26} | 13/25 = 0.52 | no |
| R. Li 2026 (arXiv 2601.00910, preprint): x^{1/2+ε}, P > x^{35/36} | 18/35 = 0.514 | no |
| Sander 1994 (k^{1+ε} ≤ x short ranges): θ = β + cβ³ | ≥ 1 | no |
| Sylvester–Schur / Laishram–Shorey: P(Δ) > 1.95k | ≈ 1 | no |

**Conclusion (PROVED for these inputs).**
- No known short-interval theorem has θ > 2β (all have β ≥ 1/2 or θ ≈ β). So H1 and H4 give
  **no** new prime-construction bound by this route.
- The route would need intervals of length h ≤ x^{1/3} that contain a prime factor in
  (h^{2+δ}, x^{2/3}]. Below β = 1/2 only Sander-type bounds θ = β + O(β³) are known.
- This is not a blocker, and it confirms "factoring hardness is not tied to prime-finding".
- Other routes: factoring ⌊z⌋! + 1, or 2^k − 1, costs poly(z) in bits (codex §4.4).
- H1 + a smooth-run conjecture (every [x, x + (log x)^A] contains a non-x^c-smooth integer)
  ⇒ polynomial-time primes. The conjecture is far beyond current methods. **CONDITIONAL.**

---

## 4. Complexity (task 1(c))

### 4.1 Straight-line programs and the τ-conjecture (both directions)

- **⇒ factoring:**
  - uniform polylog SLPs for k! ⇒ factoring ∈ FP by gcd binary search (Lipton 2009; PROVED,
    codex §5.1);
  - Lipton 1994 (in the form of Koiran–Portier–Rojas, Thm 1.14(II)): polynomials with
    R > exp(τ^ε) integer roots ⇒ non-uniform, i.o., positive-fraction factoring circuits;
  - Shub–Smale / BCSS Thm 3: τ-conjecture ⇒ P_C ≠ NP_C;
  - Bürgisser 2009, Main Thm 1.2: τ(Per_n) polynomial ⇒ τ(k!) polylog.
- **⇐ factoring: nothing.** H1 produces no polynomial with many integer roots and no short SLP
  for k!, so the τ-conjecture is consistent with H1 (codex §5.1). **PROVED (scope).**
- Blocker for *specific* mechanisms. An algorithm that itself evaluates interval products
  Π(x + i) (a Pollard–Strassen or BGS type) by a uniform SLP of exponent γ is bounded by the
  τ-conjecture: exponent c forces γ ≥ 1/c (Ward–Maestro §7.1). With c = 2, Strassen's 1/2 is
  optimal. So such an algorithm would refute τ with c ≤ 3 if it got below N^{1/6}.
  - WM (i): an exponent below 2/5 for the rising factorial would give below 1/5.

### 4.2 Möbius, AC⁰ and Sarnak

- **KNOWN.**
  - Allender–Saks–Shparlinski, JCSS 62 (2001): MAJORITY AC⁰-reduces to primality,
    squarefreeness and gcd, so none of them is in AC⁰[p].
  - Green, CPC 21 (2012), Thm 1: μ ⊥ AC⁰(d) for d ≤ c log n/log log n.
  - Bourgain, Israel J. Math. 197 (2013): the Walsh spectrum of μ is spread.
- Consequence: a factoring algorithm in AC⁰[p] is impossible, but every real algorithm uses
  multiplication or gcd anyway, so this barrier is escaped trivially.
- **Kalai's TC⁰ conjecture** (15 Dec 2024) implies factoring ∉ TC⁰ (§2.5). Modular powering
  with n-bit exponents and gcd are not known to be in TC⁰. Every known algorithm uses them.
- **Sarnak.** H1 gives μ ∈ P. This does not contradict Möbius disjointness:
  - P-computable sequences can have positive entropy (codex §5.2 gives an explicit example);
  - Sarnak's own slides ("Möbius randomness and dynamics", ICTS) say: "Perhaps ξ ∈ P ⇒ μ is
    orthogonal to ξ? I don't believe so since I believe factoring and μ itself is in P." He
    then asks for a bounded ξ ∈ P with (1/N)Σμ(n)ξ(n) → α ≠ 0.
  - Under H1, ξ = μ answers that question trivially.

### 4.3 Generic models, and the combined barriers (task 2)

- **KNOWN.**
  - Aggarwal–Maurer (Eurocrypt 2009) and Jager–Schwenk (J. Cryptol. 2013) prove
    *equivalences*: generic-ring RSA and generic Jacobi symbols reduce to factoring.
  - Jacobi symbols are easy by reciprocity, which is the canonical witness that generic-ring
    restrictions omit representation-dependent operations.
  - Damgård–Koprowski 2002 is a genuine generic-group floor, escaped by Harvey–Hittmeir, NFS
    and anything that reads N.
- No unconditional theorem is contradicted by H1–H5. **PROVED (scope; codex §5.3).**

**Lemma KE (key enlargement; codex §8.1; re-checked).**
- In T1 (below16 §2.1), replace the key bound K ≤ C₂M² by κ := max(1, K/M²).
- Then the covered fraction of any balanced prime box is ≤ C n² √κ X^{3/2} N^{−1/4}.
- So beating the exponent by δ through larger keys alone needs κ ≥ N^{3δ}/polylog.

**Theorem CS (combined signature barrier; codex §8.2; re-checked).**
- *Model.* Model X is T1 with key parameter κ. Model Y is a fixed, N-independent list of
  certificates (V_i). It splits pq only when p and q have different signatures
  σ(p) = (1[p | V_i])_i, or richer labels, as in Ward–Maestro Prop 32.2. R is the number of
  distinct signatures on a prime interval I with Π ≍ P/n primes.
- *Conclusion.* On uniform pairs from I, any X + Y algorithm failing on a fraction ε satisfies
  C n²√κ X^{3/2}N^{−1/4} ≥ Coll − ε − O(1/Π), where Coll = Σ_j (a_j/Π)² ≥ 1/R is the
  signature collision probability.
- Worst case: X ≫ N^{1/6}κ^{−1/3}R^{−2/3}n^{−4/3}.

*Proof.* Equal-signature pairs are never split by Y. They form a fraction ≥ Coll. Model X must
cover them, so apply KE. ∎

- *Use.* A worst-case sub-1/6 algorithm of this type needs κR² ≥ N^{3δ}/polylog: either keys or
  label richness of polynomial size.
- Ward–Maestro Conjecture 33.1 (an explicit ECM list with distinct labels on all primes) has
  R = Π. It is *consistent* with this barrier, and gives N^{1/(2u)}.

**Theorem CB (local collision search + (Z/N)^* smooth-order exposure; ours; CONDITIONAL on F_iv).**
- **Hypothesis F_iv** (Fouvry-type, dyadic form). For fixed 0 < a < b there is c(a,b) > 0 with
  #{p ∈ [aX, bX] prime : P⁺(p − 1) > p^{2/3}} ≥ c(a,b) X/log X for X large.
  - Fouvry, Invent. Math. 79 (1985) proves #{p ≤ x : P⁺(p−1) > x^{0.6687}} ≫ x/log x.
  - The fixed-interval form is the natural output of the same sieve. Ward–Maestro §35 uses it in
    a nearly balanced range. We did not check the interval form in the source: **UNVERIFIED**.
  - E1 measures the density as 0.33–0.36 at 20–32 bits.
- **Model.** The algorithm may do two things.
  - (X) Run T1-local certified collision searches with budget X. Classes and windows depend on N
    only through N mod K, K ≤ C₂M², and labels are charged (below16 §2.1). Group elements
    a^E mod N may serve as T1 handles, as in D′'s γ = β^{∏M²}, but not as selectors.
  - (Y₁) Compute gcd(F, N) for any F = ∏_j (a_j^{E_j} − 1), and also ± variants
    (a^E + 1 | a^{2E} − 1), where:
    - every a_j lies in a fixed base set A with Σ_{a∈A} log|a| ≤ P^{1/3}/n²;
    - every E_j is B₂-smooth, with B₂ ≤ (P/4)^{2/3} (so up to ≈ N^{1/3});
    - the choices may be adaptive, and there is no cost bound on Y₁.
  - "Success" means a certified T1 hit or a proper gcd from Y₁. Accidental collisions are
    excluded, as in T1.
- **Conclusion.** There are c₀ > 0 and C such that every such algorithm fails on at least
  (c₀ − C n² X^{3/2} N^{−1/4})·|S_n| inputs of S_n.
- In particular, success on a (1 − o(1)) fraction (indeed on a fraction > 1 − c₀/2) forces
  X ≥ c N^{1/6}/n^{4/3}.

*Proof.*
1. *A good box.* Fix the sub-box p ∈ I_p = [0.6P, 0.65P], q ∈ I_q = [0.84P, 1.19P]. Every such
   pair has p < q ≤ 2p and P²/2 ≤ pq < P², so it lies in S_n. By F_iv, the set G₀ of pairs with
   p, q ∈ F := {ℓ prime : P⁺(ℓ−1) > ℓ^{2/3}} has |G₀| ≥ c₁P²/n². Since |S_n| ≤ C P²/n², this is
   ≥ c₂|S_n|.
2. *Remove small orders.* Fix a ∈ A with |a| ≥ 2. A prime ℓ ≥ P/2, ℓ ∤ a, with
   ord_ℓ(a) ≤ (2P)^{1/3} divides a^k − 1 for some k ≤ (2P)^{1/3}. Such a^k − 1 has
   ≤ k log|a|/log(P/2) prime factors ≥ P/2, so there are ≤ (2P)^{2/3}log|a|/log(P/2) such
   primes. Summing over A gives ≤ 2^{2/3}P/(n² log(P/2)). That is o(P/log P), so the pairs with
   p or q among these primes are o(P²/n²). (If ℓ | a, then a^E ∓ 1 ≡ ∓1 is a unit mod ℓ, and
   nothing needs removing.) Call the rest G; then |G| ≥ c₀|S_n| with c₀ = c₂/2.
3. *Y₁ never fires on G.* Take (p, q) ∈ G and write ℓ_p = P⁺(p−1) > p^{2/3}.
   - Then ℓ_p² > p, so ℓ_p ∥ p − 1, and (p−1)/ℓ_p < p^{1/3} ≤ (2P)^{1/3}.
   - Since ord_p(a) > (2P)^{1/3}, ord_p(a) ∤ (p−1)/ℓ_p, so ℓ_p | ord_p(a).
   - ℓ_p > p^{2/3} ≥ (0.6P)^{2/3} > (P/4)^{2/3} ≥ B₂, so ℓ_p ∤ E for every B₂-smooth E. Hence
     a^E ≢ 1 and, since 2E is also B₂-smooth, a^E ≢ −1 (mod p). The same holds mod q.
   - So every factor a_j^{E_j} ∓ 1 is a unit mod N, and every gcd(F, N) = 1.
4. *Only T1 remains on G.* Every Y₁ answer on G is the constant 1. The adaptive algorithm
   therefore runs, on G, as a fixed T1 design. T1 (below16 §2.1, with the budget conventions
   there) covers at most C n² X^{3/2} N^{−1/4}|S_n| pairs of S_n, hence of G. ∎

Remarks.
- (a) Theorem CB is the almost-all counterpart of Theorem CS for the specific group (Z/N)^*.
  Its Coll is bounded below by the *both-unexposed* class, which has positive density even when
  R is huge.
- (b) It explains rigorously why Pollard p − 1 succeeds on a positive proportion below N^{1/6}
  (Ward–Maestro Prop 35.3) but cannot be part of an almost-all sub-1/6 algorithm built on
  (Z/N)^* alone.
- (c) **Extension CB⁺ (CONDITIONAL / SPECULATIVE).** Allow all cyclotomic tori T_d,
  d ∈ {1, 2, 3, 4, 6}, i.e. Lucas sequences and their analogues. The same proof needs a positive
  proportion of p with P⁺(Φ_d(p)) > p^{2/3} *simultaneously* for d = 1, 2. That is a
  two-shift large-prime-factor statement: p ≡ 1 (mod ℓ₁) and p ≡ −1 (mod ℓ₂) with
  ℓ₁ℓ₂ > p^{4/3}, which is beyond level-1 sieve methods. E1 measures the joint residual at
  12% (32 bits).
- (d) What escapes Theorem CB, and so what a sub-1/6 almost-all algorithm must contain:
  - order information from other groups (tori with CB⁺ open, elliptic curves, class groups),
    with an almost-all hitting theorem (ecmavg HIT; Ward–Maestro Conj 33.1);
  - relations;
  - selectors that read N beyond N mod M² (Lemma KE quantifies how far);
  - use of a^E values as selectors.

**Prop LW (Ward–Maestro Conj 33.1 is a least-witness statement; PROVED, by pigeonhole).**
- *Assume* Conj 33.1 (fixed u, C; B₁ = x^{1/u}, B₂ = B₁²): the first C_{u,C} log x Suyama
  parameters give distinct label vectors (WM Def. 32.1) on all primes in [x, Cx].
- *Conclusion.* For all but at most one prime p ∈ [x, Cx], some σ among the first
  C_{u,C} log x has λ_σ(p) ≠ 0. That is, p is degenerate on that curve, or the point order is
  ≤ B₂. So the least exposing Suyama index w(p) is O(log p) for all p except one per range.

*Proof.* Label 0 means "unexposed". Distinctness allows at most one all-zero vector. ∎

- *Comparison.* For least quadratic non-residues, a bound O(log p) for all p is false under GRH
  (Montgomery: ≫ log p log log p infinitely often). GRH gives 2 log² p (Bach).
- For smooth-order witnesses no GRH-type bound of any polylog strength is known. Conj 33.1 is
  heuristically consistent: the events are not multiplicative, and independence gives
  max w ≈ log π(x)/|log(1−θ)|. But proving it is a least-witness theorem at log scale.
  **Blocker** (for that route to a worst-case result).
- The almost-all version needs only Pr_p[w(p) > L] → 0 with a rate. That is ecmavg's HIT
  hypothesis, open.

---

## 5. Escape-signature table (task 2)

Columns:
- **must break**: the assumption that any polynomial-time, or deterministic sub-1/6, algorithm
  of the barrier's type must violate;
- **witnesses**: known algorithms that already escape it.

Abbreviations:
- PS: Pollard–Strassen / BGS;
- LP: Lenstra–Pomerance class groups (and Schnorr–Lenstra, Seysen);
- Copp: Coppersmith;
- p−1⁺: Pollard p − 1 on a positive proportion (Ward–Maestro Prop 35.3).

Full model statements are in codex §7.1 (re-checked against the notes).

| Barrier (floor) | Model assumptions | Must break | Witnesses |
|---|---|---|---|
| WM Cor 3.3 (N^{1/3} candidates), Thm 3.4 (N^{1/5}, materialised giants), Thm 9.4 (fixed degree), Thm 2.8 (hypothesis thinning) | Oblivious coverage of *every* real p in J_C by Lehman cells aq + bp; materialised giants; one common modulus | non-oblivious exclusion; free giants; non-cell tests | D′/AV (many moduli, true factor only); PS; p−1⁺; ECM; QS/NFS; LP; Shor; Copp |
| WM Thm 3.6 / multimod Thm A (N^{1/6}, free giants, any moduli) | sumset collision search; full-candidate coverage; Σ W ≤ X² | full coverage → true-factor only; or leave sumset search | p−1⁺, ECM, QS/NFS, LP, Shor (D′ attains; does not break) |
| WM Thm 4.4 / Prop 4.5 (explicit difference covers; PS optimal for divisibility tests) | explicit B − G representation; start-statistic envelope | implicit product evaluation (WM 7.2) | none below 1/5 known |
| WM Thm 10.1 (piecewise evaluation ≥ N^{1/4}/√H) | pieces priced ≥ √length separately | shared computation across pieces | none known |
| multimod A′/B (deg ≥ 2 ⇒ N^{1/5}; J moduli ⇒ N^{1/5}J^{−1/5}) | as Thm A | degree-1 targets with N^{Ω(1)} moduli | D′ |
| below16 T1/T2/T3 (true factor, local, N^{1/6}/polylog) | key N mod K, K ≤ CM²; ≤ n/2 size bits; charged labels; degree 1 | locality (Lemma KE: κ ≥ N^{3δ}); or leave sumset | p−1⁺, specialfam/Copp families, LP, QS (heuristic) |
| **Theorem CB** (new; T1 + (Z/N)^* smooth exponents) | as T1, plus B₂-smooth annihilation with B₂ ≤ N^{1/3} | orders from other groups with a hitting theorem; relations; non-local selection | ECM (randomized, heuristic L_p[1/2]); LP; QS/NFS; Shor |
| **Theorem CS** (codex; T1 + fixed list with R signatures) | fixed-list certificates | R ≥ N^{3δ/2} label richness, or κ large | WM Conj 33.1 lists (conjectural); randomized ECM (Prop 33.2) |
| coverage38 Thm U (worst case for DCL designs) | N-independent (weight, modulus) groups; WIN | coverage at λ < J^{1/4} for every pair; or WIN-violating windows | none |
| coverage38 Prop L; dchworst B1 | phase-discarding linear statistics; error ≳ √J | arithmetic structure x_ℓ = pℓ̄ mod q | none (proof-method barrier) |
| dchworst B2 | single-modulus pigeonhole | union over moduli | D′ (on average) |
| dchworst B3 | proof valid for all coprime pairs | "covered or easy" dichotomy | coverage38 dichotomy (partial) |
| sieve2d Prop 5 | only ≤ J-prime statistics | structure beyond moments | none |
| unbalanced Prop O′ (min(√p, q^{1/3})) | 1-D Lehman-linear full coverage | average-case / non-𝔏 | PS (first branch), ECM, NFS, p−1⁺ |
| Track C O1–O4 | finite-algebra powers, power maps, subgroup zero sets for unit predicates | non-faithful or nonlinear readouts; different task | PS/BGS (√ cost); Tal (needs p); every factoring algorithm avoids the task |
| Track A A1 | rank-D linear algebra observables on F_N | not reading F_N; sparse circuits | every algorithm |
| Track E information floor | per-modulus data | non-per-modulus computation | every algorithm; Copp needs n/4 coset-shaped bits |
| eps battery | finite menu | untested statistics | CFRAC/QS relations (L[1/2]) |
| nfsfruit Lemma M; lemmaM B1, (c) | NFS fruitfulness; even-|S| two-form sieves; one-level HECM descent (L[4/9]) | fibre symmetry; multilevel descent | Dixon / Theorem DD (randomized); LP |
| hypD Props G, H | slicing; bounded-k shape test | new averaged estimates (FM, PC) | none rigorous below L[1/2, 1] |
| generic group (Shoup; Damgård–Koprowski) | representation-oblivious | read the integer N | HH, NFS, everything |
| generic ring (AM09, JS13) | equivalence, not a floor | — | Jacobi via reciprocity |
| AC⁰[p] (ASS, Green, Bourgain) | bounded depth | multiplication or gcd | every algorithm |
| τ-conjecture (as a floor for SLP product evaluation) | uniform SLP for Π(x+i) | not a single SLP; branching on representation | ECM, QS (branch on lifts) |
| Theorem DR contrapositive | — | — | conditional deterministic L[1/2] escapes *all* model floors at once |

**Escape signature of a deterministic N^{1/6−δ} algorithm (almost all or worst case).**
Combining rows 2, 6–8 and Theorem DR:
- (S1) It is not a full-candidate sumset cover (Thm A / WM 3.6).
- (S2) If it uses sumset collisions, either its selector is non-local with κ ≥ N^{3δ} (KE), or
  collisions carry only part of the success.
- (S3) The rest of the success comes from:
  - (S3a) group orders outside (Z/N)^* (CB), with an almost-all hitting theorem, and for the
    worst case label richness R ≥ N^{3δ/2} (CS) and a least-witness theorem (LW); or
  - (S3b) relations with deterministic smooth supply, deterministic fruitfulness and
    deterministic linear algebra; or
  - (S3c) something outside all listed models.
- (S4) It must use representation-dependent operations (not generic) of super-constant depth
  (not AC⁰[p]; not TC⁰ if Kalai's conjecture holds).

**Witnesses for each ingredient alone.**
- (S1)–(S2) with sub-1/6 on a positive proportion: p−1⁺, specialfam.
- (S3a) randomized: ECM (heuristic), LP (rigorous, class groups).
- (S3b) randomized: Dixon, LP (rigorous); deterministic heuristic: det-QS (below16 §6).
- Conditional deterministic: Theorem DR (E-hardness), WM Conj 33.1 (N^{1/(2u)}), det-QS
  (H1 + H2 of below16).
- **No family has the whole signature unconditionally.**

---

## 6. Blockers (task 3(i))

1. **Least witness (Prop LW).** Worst-case deterministic ECM-type routes, including
   Ward–Maestro Conj 33.1, imply O(log p) least-witness bounds for smooth-order exposure. No
   analogous bound is known even under GRH, and for QNRs it is false.
2. **Two-shift sieve (CB⁺).** Even the barrier for cyclotomic portfolios needs joint large prime
   factors of p − 1 and p + 1, beyond level 1. Conversely, an algorithm that proves its
   cyclotomic portfolio exposes almost all p would prove a joint *smoothness* statement of the
   same depth. E1 says the joint residual does not decay (12% at 32 bits), so no such
   statement is true at the needed strength.
3. **Coverage theorems (H4w via Theorem R).** dchworst B1–B3 and coverage38 U/L, unchanged: a
   sparse covering statement U(x, Y^{1/8}) beyond known methods.
4. **τ-type.** A product-evaluation (PS/BGS-type) algorithm below N^{1/6} refutes the
   τ-conjecture with exponent c ≤ 3 (§4.1).
5. **Lower bounds are blocked (Theorem DR).** Any proof that deterministic factoring needs
   N^{Ω(1)} would put E in i.o.-SIZE(2^{εm}) for every ε. So the repository's model floors
   can never be "upgraded" to general deterministic lower bounds.
6. Not blockers:
   - prime construction (§3);
   - Sarnak, τ, μ ⊥ AC⁰ (§4);
   - analytic NT (§2.1).

## 7. New directions (task 3(iii)), ranked

**D1. Fixed-u deterministic QS on almost all of S_n.** Feasibility low–medium; payoff high.
- Below 1/6 needs a factor base y = N^c with c < 1/(6ω), so the smoothness density ρ(u) is a
  *constant* (below16 §6.3). The fixed-u regime is far more amenable than L[1/2]: constant
  densities, and first and second moments over N.
- *First lemma.* The first moment over N ∈ S_n of #{t ≤ T : (⌈√N⌉+t)² − N is y-smooth}. For
  fixed t, the values over N with ⌈√N⌉ = m fill [2mt + t², 2m(t+1) + t²). Over all m the count
  is a long-range smooth count, times a Bombieri–Vinogradov input for semiprimes in
  {x² − v}.
- *Second step.* The second moment over **all** integers N. For fixed (t, t′) the map
  (m, m² − N) ↦ (v, w) is linear with determinant 2(t − t′). So the pair count is smooth-pair
  counting on a sublattice of index 2(t′ − t), i.e. smooth numbers in APs mod 2k ≤ 2T ≈ 2y.
  That is within Harper's y^{4√e−δ} individual-modulus range and the BDH-type large sieve.
- *Transfer (PROVED, trivial).* An exceptional set of density ε_n among all n-bit integers is a
  fraction ≤ Cε_n n² of S_n. So almost-all-integers with ε_n = o(n^{−2}) suffices, and the
  binary "semiprime + smooth = square" problem of below16 is bypassed.
- *Kill criteria.*
  - The error terms available for smooth numbers in short multiplicative windows of relative
    width 1/T² (needed for k = t′ − t small) give only (log)^{−c} savings with c ≤ 2. Then the
    transfer fails.
  - The common-mode α-variation over N mod Q (ecmavg-type) cannot be stratified away.
- Even if D1 succeeds, fruitfulness (below16 H2) remains. Codex §9.4 turns it into one rank test.

**D2. Group-order hitting outside (Z/N)^*** (forced by CB). Feasibility low.
- The cheapest instance: one CM-curve family, or a small set of tori, together with a
  *two-shift* large-prime-factor theorem.
- *First lemma.* Whether Fouvry's method gives a positive proportion of p with P⁺(p−1) and
  P⁺(p+1) both > p^{1/2+η}. That is the weakest version CB⁺ needs at budget N^{1/8}.
- *Kill criterion.* Level-of-distribution arithmetic: ℓ₁ℓ₂ > p^{1+2η} exceeds BV, so it is
  likely dead unless one shift is handled by a bilinear (Type II) structure.
- The worst-case version is a least-witness problem (LW): kill.

**D3. Non-local selector audit (codex Rank 1).** Feasibility high for an audit, low for a gain.
- KE turns any proposed selector into one inequality, κ ≥ N^{3δ}.
- *Kill criterion.* The key compresses back to M²N^{o(1)}, or window placement needs factor
  information.

**D4. Promise-preserving partial-information reductions (codex Rank 2).** Feasibility low.
- ε ⇔ φ mod 8. Target: a reduction from Jacobi-(+1) QR to ε on constructed semiprimes.
- *Kill criterion.* Every transformation is a coprime mask (codex §9.3) or needs fresh factors.

**D5. Unconditional targeted derandomization of Lenstra–Pomerance or Dixon at the N^{1/6}
scale.** Feasibility very low.
- Theorem DR shows that only derandomization separates us from L[1/2]. At the N^{1/6} scale a
  generator needs seed ≤ n/6 bits against one *specific* test.
- For Dixon that test is "r² mod N is y-smooth". The structured seed set r = ⌈√(kN)⌉ + t *is*
  det-QS, so D5 reduces to D1. For LP it is "h(−kN) is smooth" for explicit multipliers k, an
  open short-family smoothness problem.
- *Kill criterion.* Reduce to D1, or to a short-family Frobenius-equidistribution problem
  (families of size p^{o(1)}; known methods need ≥ p^{1/2}).

**D6. Residual-tail portfolios (codex Rank 4).** **Refuted for cyclotomic portfolios
(heuristic + E1).**
- The mean bound N^a + ηN^{1/5} needs η ≤ N^{−1/30−δ}. Cyclotomic residues stay ≈ 1–1.1% of
  pairs (E1).
- Only list-length-growing ECM portfolios remain: ecmavg, WM Part 5.

No direction has a test that is both cheap and decisive beyond E1. D3 and D4 are audits of
specific proposals, and D1/D2 need proofs, not data.

---

## 8. Experiment E1 (pre-registered, `notes/polyimpl_prereg.md`, commit `1a9c3de`)

**Setup.**
- 20 000 uniform primes per size k (bits).
- B₂ = ⌊2^{2k/3}⌋, i.e. X² with X = N^{1/6}.
- Generous exposure: p is exposed by d iff P⁺(Φ_d(p)) ≤ B₂.
- 10 000 balanced pairs per size.

| k | B₂ | r₁⁽¹⁾ (P⁺(p−1) > B₂) | r₁⁽¹,²⁾ | r₁ (all d ∈ {1,2,3,4,6}) | r₂ (pairs) | unexposed by d = 3, 4, 6 |
|---|---|---|---|---|---|---|
| 20 | 10 321 | 0.330 | 0.095 | 0.083 ± 0.002 | 0.0069 ± 0.0008 | 0.961, 0.958, 0.965 |
| 26 | 165 140 | 0.350 | 0.115 | 0.102 ± 0.002 | 0.0108 ± 0.0010 | 0.964, 0.959, 0.963 |
| 32 | 2 642 245 | 0.358 | 0.123 | 0.109 ± 0.002 | 0.0113 ± 0.0011 | 0.958, 0.955, 0.959 |
| model | | 0.405 | 0.164 | 0.141 | 0.020 | 0.951 |

**Frozen checks.**
- P1 holds.
- P2, P3, P4 fail at k = 20 only, on the low side: 0.095 < 0.12, 0.083 < 0.09,
  0.0069 < 0.008. They hold at k = 26 and 32.
- P5 fails: r₁(32)/r₁(20) = 1.31 > 1.3, by **growth**.
- **Pre-registered verdict: "mixed".** Neither branch of the decision rule applies: the residual
  did not shrink, but the stability band was missed.

**Post hoc reading (labelled).**
- The residual increases with k toward the independence model, as finite-size Dickman
  corrections predict. Pr[P⁺(n) > n^{2/3}] approaches ln(3/2) from below.
- p − 1 and p + 1 are nearly independent: r₁⁽¹,²⁾ = 0.1226 against r₁⁽¹⁾r₁⁽²⁾ = 0.1269 at
  k = 32.
- The failure direction is the opposite of the one that would undercut CB⁺ or rescue D6. So
  CB⁺ remains experimentally supported only informally (not by the frozen rule), and D6 is
  refuted for cyclotomic portfolios on heuristic grounds.

Reproduction: `PYTHONPATH=. .venv/bin/python -m src.polyimpl.residual run` (≈ 2 s);
`PYTHONPATH=. .venv/bin/python -m src.polyimpl.primefind`; `.venv/bin/python -m pytest -q tests/test_polyimpl.py`.

## 9. Reviews

See §9.1 (Opus critic) and §9.2 (GPT-6.1-sol adversarial review), filled in after the reviews.
