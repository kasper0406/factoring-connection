# polyimpl: what a fast factoring algorithm would imply, and what that says about blockers and routes

Agent `polyimpl`, phase 7, 2026-10-19.
- Pre-registration: `notes/polyimpl_prereg.md` (commit `1a9c3de`, made before E1 ran).
- Code: `src/polyimpl/{residual,primefind}.py`. Tests: `tests/test_polyimpl.py`.
- Raw data: `experiments/raw/polyimpl_e1.jsonl`. Summary:
  `experiments/analysis/polyimpl_e1_summary.json`.
- GPT-6.1-sol: blind round `codex/polyimpl_1.{md,out}` (the `.out` is the full report);
  adversarial referee `codex/polyimpl_2.{md,out}`.
- Independent read-only Opus critic: findings in §9, all applied.
- Status labels:
  - **KNOWN**, with a reference;
  - **PROVED**, with the proof here, routine unless said otherwise;
  - **CONDITIONAL**, **SPECULATIVE**, **FALSE**, **experimentally supported**.
- "Codex §x" means a section of `codex/polyimpl_1.out`. Its proofs were re-checked here.

## 0. Verdict (plain language)

**We found no consequence of fast factoring outside factorization-derived computation.** This is
an editorial assessment of the literature we checked, not a theorem.

What polynomial-time factoring (H1, H2) would give:
- fast versions of the factorization-derived functions μ, φ, σ, the squarefree parts, maximal
  orders, order finding, #E(Z/N) and τ(N);
- the end of RSA- and Rabin-type cryptography.

What it would not give, as far as we can tell:
- No implication is established for class groups, prime-field discrete logs, deterministic
  square roots mod p, or any theorem about primes.
- No contradiction appears with any unconditional theorem or named conjecture we checked:
  Sarnak, τ, μ ⊥ AC⁰, Kalai's TC⁰ conjecture (unless the algorithm is in TC⁰).
- So mathematics supplies **no free lower bound** for factoring.

**No blocker from prime construction (Lemma PF, §3).**
- The square-root-length "large prime factor in a short interval" theorems we checked need no
  factoring at all: the large prime is stripped out with one remainder tree.
- The shorter-interval theorems known (Sylvester–Schur and Sander type) are too weak.
- So none of these inputs turns an H1 or H4 algorithm into a better deterministic prime
  construction than z^{1/2+o(1)}.

**A folklore conditional upper bound, worth keeping in view (Theorem DR, §1).**
- Assume almost-everywhere circuit hardness for E (Impagliazzo–Wigderson 1997).
- Then derandomizing Lenstra–Pomerance gives deterministic L_N[1/2, O(1)] factoring of every N.
- So the project's unconditional deterministic exponents, 1/5 worst case and 1/6 average, sit
  far above what the standard derandomization hypothesis predicts.
- An unrestricted deterministic N^{Ω(1)} lower bound would refute that hypothesis.

**Escape signature (§§4.3, 5).**
- Every model floor we list has some known algorithm outside its model, except the proof-method
  rows (coverage38, dchworst B1, sieve2d Prop 5), which have none.
- But no known classical algorithm is unconditionally deterministic, rigorous and below N^{1/6}
  on almost all (or all) balanced semiprimes.
- Two combined model barriers constrain such an algorithm. Both are routine combinations of
  known ingredients, "no prior formulation found by us":
  - **Theorem CB** (CONDITIONAL on the Fouvry-type interval bound F_iv).
    - The algorithm runs a T1-local collision search with a per-input budget X.
    - It may also run the (Z/N)^* tests gcd(∏(a^E ∓ 1), N), with B₂-smooth E (B₂ ≤ N^{1/3}),
      up to ≈ N^{1/6}/n² bases, adaptively. Order collisions a^E ≡ a^{E′} are excluded.
    - Under a stated closure condition, it fails on a positive proportion of S_n unless
      X ≥ N^{1/6}/polylog.
  - **Theorem CS** (codex, re-checked).
    - The algorithm combines a T1-local search with a fixed certificate list having R prime
      signatures, and the collision design does not adapt to the list's answers.
    - Then worst-case success needs X ≫ N^{1/6}R^{−2/3}/polylog.
- Inside these models a sub-1/6 almost-all algorithm needs at least one of:
  - order information from other groups, with a hitting theorem;
  - relations with a deterministic smooth supply and fruitfulness;
  - a non-local selector;
  - something outside every listed model.

**Experiment E1 (pre-registered; frozen verdict "mixed").**
- It measures whole-order smoothness of the five cyclotomic tori at the N^{1/6} budget.
- The residual of unexposed primes is 8.3% at 20 bits, 10.2% at 26 bits and 10.9% at 32 bits.
- For pairs it is 0.69%, 1.08% and 1.13%.
- Post hoc: no decreasing trend at these sizes. That is evidence against a fixed cyclotomic
  portfolio having the polynomially small residual codex's Rank-4 route needs. It is not an
  asymptotic statement.

**Directions (§7).** Two ranked lists:
- **Ours:**
  1. D1, fixed-u deterministic QS on almost all N, via all-integer averaging plus the transfer
     lemma;
  2. D3, an audit of non-local selectors;
  3. D2, hitting outside (Z/N)^*;
  4. D5, targeted derandomization;
  5. D4, promise-preserving partial-information reductions.
- **Codex's referee:** D3, D1, D5, D2, D4, D6.
- No direction has a test that is both cheap and decisive beyond E1.

---

## 1. Hypotheses, transfer, and Theorem DR

The hypotheses:
- (H1) deterministic polynomial time;
- (H2) randomized polynomial time (classical);
- (H3) L[1/3 − δ] or L[1/4], heuristic or rigorous;
- (H4w) deterministic worst case N^{1/5−δ};
- (H4a) deterministic average, or almost-all, N^{1/6−δ} on S_n. "Average" and "almost all" are
  different hypotheses: an unspecified exceptional density does not bound the mean with a slow
  fallback;
- (H5) partial-information oracles.

**PROVED (transfer rule; codex §2).** A reduction making poly(n) factoring calls on integers
≤ N^c turns N^α into N^{cα+o(1)}, and L[a, C] into L[a, Cc^a].
- H4a does **not** transfer to reductions that query integers outside its ensemble.
- An exceptional set of density N^{−Ω(1)} can contain any prescribed short interval.

**PROVED (no known collapse).** The threshold language "N has a factor ≤ k" and its complement
have unique certificates: the sorted factorization, verified with AKS.
- So that language is in UP ∩ coUP ⊆ NP ∩ coNP. Factoring is also in BQP (Shor 1997).
- No known class collapse follows from H1, and none of the reductions here gives P = BPP.

**Theorem DR (PROVED; a routine corollary of two KNOWN theorems; folklore).**

*Hypothesis (IW).* There are H ∈ E = DTIME(2^{O(m)}), η > 0 and m₀ such that, for **all**
m ≥ m₀, H restricted to length m needs circuits of size ≥ 2^{ηm}. This is the hypothesis of
Impagliazzo–Wigderson, STOC 1997, pp. 220–229.

*Conclusion.* A uniform deterministic algorithm outputs the complete factorization of every N
in time L_N[1/2, c], for a constant c that depends on H. Under the same hypothesis, H2 ⇒ H1.

*Proof.*
1. Lenstra–Pomerance (JAMS 5 (1992), unnumbered main theorem of §1) is Las Vegas, with no GRH.
   Fix an explicit clock U(n) = exp(C₀√(n log n)) that bounds its expected time on every n-bit
   input.
2. Truncate the algorithm after 2U(n) steps, counting random bits as steps. By Markov it
   succeeds with probability ≥ 1/2.
3. Add verification of the factorization: multiply out and run AKS. With N hardwired, this is a
   circuit D_N of size S = U^{O(1)} on a random tape padded to length S.
4. IW gives G : {0,1}^{O(log S)} → {0,1}^S, computable in time S^{O(1)}, which fools every
   size-S circuit with error < 1/4.
5. So some seed s has D_N(G(s)) = 1. Enumerate all S^{O(1)} seeds.
6. The total cost is S^{O(1)} = L_N[1/2, O(1)]. The leading constant 1 is not preserved. ∎

*Contrapositive.* Suppose **no** deterministic algorithm factors every N in time L_N[1/2, C],
for any constant C. Then for every H ∈ E, every η > 0 and every m₀, some m ≥ m₀ has
circuit complexity of H at length m below 2^{ηm}.
- That is, E ⊆ ⋂_η i.o.-SIZE(2^{ηm}), language by language and exponent by exponent.
- So an unrestricted deterministic lower bound would refute the IW hypothesis.

*Relation to known statements.*
- below16 §6.2(3) and Ward–Maestro §35 note the non-uniform deterministic L[1/2] (Adleman).
- Polymath4 notes the prime-construction analogue.
- Factoring has a unique output, so "pseudodeterministic factoring" is simply Lenstra–Pomerance.
  Compare Chen–Lu–Oliveira–Ren–Santhanam 2023 for primes.

---

## 2. Consequence catalog by hypothesis

Throughout, "no implication established" means none is known to us, not a non-implication
theorem.

### 2.1 (H1) deterministic polynomial time

| Area | Consequence | Label / reference |
|---|---|---|
| Arithmetic functions | μ, Liouville λ, ω, Ω, d, φ, σ_k (fixed k), Carmichael λ, rad(N), squarefree decomposition N = ab² — all in FP | PROVED (codex §3.1) |
| Reverse reductions | σ(N) ⇒ factoring, randomized | KNOWN: Bach–Miller–Shallit, SIAM J. Comput. 15 (1986), Thm 2 |
| | φ or λ ⇒ factoring, randomized (deterministic under ERH: Miller 1976). For N = pq, φ gives s = N + 1 − φ, deterministically | KNOWN; PROVED for pq |
| | Żrałek, Math. Comp. 79 (2010), Cor. 8.2: O(log N) queries along the *totient chain* N, φ(N), φ(φ(N)), … ⇒ deterministic full factorization in L[1/3, 1+o(1)]. It is not one supplied value φ(N); Thm 8.1 assumes φ(N) is already factored | KNOWN |
| | Shallit–Shamir, IPL 20 (1985): one call d(N^{q−1}), q prime > 1 + log₂N, decides squarefreeness and gives μ. d ⇔ the exponent multiset, for an all-integer oracle | KNOWN (codex §9.3 proof) |
| | no known reduction of factoring to μ, squarefreeness, ω or d (general binary inputs) | KNOWN-open, scoped |
| Maximal orders | Buchmann–Lenstra, JTNB 6 (1994): Thm 1.3, maximal-order *construction* ⇔ rad(d) for arbitrary d; Thm 6.12, *recognition* ⇔ squarefreeness. So H1 ⇒ deterministic polynomial-time maximal orders. A squarefreeness bit alone does not suffice | KNOWN |
| Class groups, units | no reduction from class-group or unit computation to factoring is known. Imaginary quadratic class groups stay at L[1/2] under GRH (Hafner–McCurley 1989). For negative fundamental D, H1 gives the 2-rank ω(D) − 1 | KNOWN (scoped) |
| Residuosity, roots | QR decision mod N is in P | PROVED |
| | root *construction* given the factors is randomized polynomial time. Tonelli–Shanks still needs a non-residue (GRH: Ankeny–Bach; Lamzouri–Li–Soundararajan sharpen the constant). Fixed radicands are deterministic (Schoof 1985) | PROVED / KNOWN |
| | root search ⇒ factoring, randomized (Rabin 1979) | KNOWN |
| Order finding | ord_N(a) in FP, using further H1 calls on p − 1 | PROVED (codex §3.2) |
| | order ⇒ factoring, randomized | KNOWN |
| DLP | Bach 1984 (randomized reductions): composite DLP ⇒ factoring; composite DLP ≤ factoring + prime DLP. No prime-field DLP algorithm follows | KNOWN |
| Two squares | r₂(N) in FP. One representation is deterministic (Schoof √−1 + Cornacchia) | PROVED + KNOWN |
| Point counting | #E(Z/N), good reduction at all p ∣ N: #E(Z/p^e) = p^{e−1}#E(F_p), then CRT + Schoof | PROVED (codex §3.4) |
| | Dieulefait–Urroz (arXiv 1911.11004), Thm 4: #E and #E^d mod N = pq, with (d/p) = −1 at a factor, ⇒ deterministic factoring. Thm 5 finds d under ERH. Kunihiro–Koyama 1998 is probabilistic and has a gap they point out | KNOWN |
| τ(N) | τ(N) in FP (Edixhoven–Couveignes at primes, recurrence, multiplicativity) | KNOWN / PROVED |
| | Bach–Charles (arXiv 0708.1192), Thm 2.2: a uniform algorithm for a(n) of a fixed non-CM eigenform, queried at n and n², factors a density-one set of RSA moduli. Not single-value worst-case hardness | KNOWN (checked by codex) |
| Constructions | primitive roots become testable. Finding one: p^{1/4+ε} deterministic (Burgess; Shparlinski 1996), polylog under GRH (Shoup). No reduction for irreducible polynomials found | KNOWN / SPECULATIVE |
| Prime construction | no new bound from the checked inputs (§3) | PROVED for those inputs |
| Analytic NT | faster single evaluations of μ, Λ or [P⁺(n) ≤ y] give no cancellation in sums. Mertens sums are already O(x^{3/5}(log x)^{8/5+ε}) without factoring (Helfgott–Thompson 2023, verified by codex). Only finite verification gets cheaper (odd-perfect bounds, aliquot sequences, Cunningham tables). We found no other consequence | PROVED (codex §6) + assessment |
| Crypto | RSA, Rabin, Goldwasser–Micali, Paillier and RSA-group trapdoors are broken on the covered keys. BBS loses its security assumption; an exposed final state is reversible (codex §5.4). Class-group, lattice, hash and code assumptions are not refuted by this reduction | KNOWN / PROVED |

### 2.2 (H2) randomized polynomial time

- §2.1 holds with "randomized". Under the IW hypothesis, H2 ⇒ H1 (Theorem DR).
- Naor–Reingold, JACM 51 (2004), has both factoring-based and DDH-based PRFs in TC⁰. Only the
  factoring-based instance of the natural-proofs barrier loses its premise. The other
  instances survive only under their own assumptions, at the required (subexponential)
  security.
- Learning hardness: only the factoring-based instances lose their premise (Kearns–Valiant
  1994; Kharitonov 1993).
- Quantum evidence: factoring stops being evidence for BQP ⊄ BPP. The other evidence is
  oracle-relative or conjectural. KNOWN.

### 2.3 (H3) L[1/3 − δ] or L[1/4]

- No mathematical consequence is established. A rigorous H3 would be a new proved factoring
  bound; it would not by itself give new prime-field DLP algorithms.
- An L[a, C] attack rules out "factoring is hard for circuits of size 2^{n^ε}" when ε > a.
  Heuristic NFS (a = 1/3) is not a proof, so it contradicts nothing unconditionally.
- Analogy only (SPECULATIVE): the L[1/4] and quasi-polynomial jumps for small-characteristic
  DLP (Joux 2013; BGJT 2014) used the subfield Frobenius x^q − x. Z has none.
- In-repo barriers for a *rigorous* L[1/3]: nfsfruit Lemma M (an open sufficient condition),
  the lemmaM (c) HECM accounting (L[4/9] for that architecture), and hypD Prop H.

### 2.4 (H4) the project's own frontiers

**(H4w) Deterministic worst case N^{1/5−δ}.**
- It gives deterministic N^{1/5−δ} versions of the §2.1 functions on the same inputs, including
  order finding. For balanced N, the extra calls on p − 1 cost N^{1/10}; on general inputs,
  charge the H4w calls on p − 1.
- No other consequence is established, and in particular no prime-construction gain (§3).
- Obligations apply only to algorithms that use the respective coverage certificates and cost
  models:
  - Ward–Maestro Thm 3.4 and Thm 9.4, and cycloharvey LB;
  - multimod A′/B: N^{Ω(1)} moduli;
  - coverage38 U: coverage at λ < J^{1/4} for every pair;
  - dchworst B1–B3, via Theorem R.

**(H4a) Average or almost-all N^{1/6−δ}.**
- In-model obligations:
  - T1/T2/T3 (each locality class separately);
  - Theorem CB (per-input budgets, closure condition);
  - Theorem CS, only when its failure fraction is o(c_eq);
  - multimod A for full-candidate covers.
- Expected time on S_n does not by itself satisfy T1's all-integer moment budgets.

### 2.5 (H5) partial information

| Oracle | On S_n | Label |
|---|---|---|
| ω, Ω, d, μ, μ², rad on N | constant (2, 2, 4, 1, 1, N): 0 bits | PROVED |
| the same on coprime masks Nk, Nk² | determined by k | PROVED (codex §9.3) |
| two-square decision, r₂(N) | for N ≡ 1 (mod 4): ε = (−1/p), equivalently φ(N) mod 8 ∈ {0, 4}. For N ≡ 3 (mod 4): constant "no", and (−1/p) depends on which factor is called p | PROVED (codex §§3.3, 9.2) |
| | no reduction from factoring known (Seurin's 2-Φ/4-Hiding; codex r6) | KNOWN-open |
| φ(N) mod m for all primes m ≤ Cn | CRT recovers φ(N), so pq is factored deterministically. General N needs randomization (BMS) or totient-chain access (Żrałek) | PROVED / KNOWN |
| φ(N) mod 2^k, k ≥ n/4 + v₂(q − p) | Coppersmith, 2 decodes | PROVED (Track E P-6) |
| φ(N) mod one fixed m; Jacobi-type bits from N mod m | s mod m only; 0 bits beyond N mod m (information floor) | PROVED (Track E F1/F2) |
| factoring in polynomial-size TC⁰ (all inputs) | gives μ in TC⁰, refuting Kalai's conjecture (15 Dec 2024 post, item 3). Take f = μ where μ ≠ 0 and f = 1 elsewhere; then the mean of μf is 6/π². Nothing follows for sub-1/6-time algorithms or for promise algorithms on S_n | PROVED (codex) |

Among this oracle menu, the oracles that factor S_n are uniform families of residues of the
symmetric function s (φ mod m for many m). One-bit and multiplicative-function oracles on N are
constant, or equal ε.

---

## 3. Deterministic prime construction (task 1(b))

**KNOWN.** Polymath, "Deterministic methods to find primes", Math. Comp. 81 (2012).
- It cites the deterministic constructions z^{0.525+o(1)} (BHP gaps) and z^{1/2+o(1)} (RH, or
  Lagarias–Odlyzko).
- Its main result is Theorem 1.2, a parity computation; Conjecture 1.1 is the existence version.
- Nothing on factoring oracles is stated there.
- Polymath4 research thread II (9 Aug 2009), "Strategy 3": with an O(1) factoring oracle, find
  intervals [n, n + a] guaranteed to contain a non-S-smooth number. "S > a²" would beat the RH
  bound. No oracle bound was obtained.

**Lemma PF (PROVED; `src/polyimpl/primefind.py`; the first route is codex §4.2).**
- *Input theorem T(β, θ).* **Every** sufficiently large interval [X, X + X^β] contains n with
  P⁺(n) > X^θ. Almost-all-interval results have the wrong quantifier.
- *Factoring cost.* Deterministic full factorization of n ≤ 2X costs X^{α+o(1)}.
- *Conclusion.* A prime p with z := X^θ < p ≤ 2z^{1/θ} is found in time z^{e+o(1)},
  e = min(e_fac, e_strip):
  - e_fac = (β + α)/θ, by factoring every integer of the interval;
  - e_strip = max(β, 1 − θ)/θ, for θ > 1/2 and **with no factoring**. Let H = ⌈2X^{1−θ}⌉ and
    Q = (∏_{ℓ≤H} ℓ)^{⌈log₂ 2X⌉}. A remainder tree gives Q mod n for all n in the interval;
    replace n by n/gcd(n, Q^∞). When P⁺(n) > X^θ > H, the remainder is the prime P⁺(n). The
    cost is (X^β + X^{1−θ})^{1+o(1)}. ∎

**Corollary PF (necessary conditions only; not a classification).**
- (i) Beating 1/2 for some α ≥ 0 needs θ > 2β.
- (ii) A factoring algorithm is *needed* to beat 1/2 only if, in addition, the strip route
  fails. With (i) this means 2β < θ ≤ 2/3, hence β < 1/3. (When e_strip < 1/2, factoring
  can still improve on it if α < 1 − θ − β.)
- (iii) H4w at α = 1/5 − δ needs θ > 2β + 2/5 − 2δ, hence β < 2/15 + δ.

**Checked inputs** (the literature agent read Merikoski's survey table and the Shorey–Tijdeman
survey; codex checked the primary sources):

| Theorem (β, θ) | best e | factoring needed? |
|---|---|---|
| Baker–Harman–Pintz 2001: a prime in [X − X^{0.525}, X] | 0.525 (primality test) | no |
| Merikoski 2021, MPCPS 170, Thm 1: (1/2, 18/19), up to log factors | 19/36 = 0.528 | no (strip) |
| Jia–Liu 2000, Acta Arith. 95: (1/2+ε, 25/26−ε) | 13/25 + O(ε) | no |
| R. Li, arXiv 2601.00910 (preprint): (1/2+ε, 35/36−ε) | 18/35 + O(ε) = 0.514 | no |
| Sander 1994, "On numbers with a large prime power factor" (j = 1 case): θ = β(1 + cβ²), for β ≥ (log x)^{−1/2+o(1)}-type ranges | 1/(1 + cβ²) at α = 0 | **yes** (θ < 1/2) |
| Sylvester–Schur / Laishram–Shorey: P(Δ) > 1.95k | ≈ 1 | — |

**Conclusion.**
- None of the checked inputs gives e < 1/2. For the square-root-length inputs this holds
  whatever α is.
- For Sander-type inputs, H1 would give e = 1/(1 + cβ²). That is < 1/2 only if cβ² > 1 in the
  admissible range. We did not determine Sander's constant c (**UNVERIFIED**). For small β we
  expect cβ² ≪ 1, but we have not checked.
- This is not a blocker for H1 or H4 on the checked inputs. It is a literature check, not a
  proof that no route exists.
- Other routes: factoring ⌊z⌋! + 1 or 2^k − 1 costs poly(z) in bits (codex §4.4).
  Polynomial and multiplicative families need their own count, e = (b + dα)/λ (codex). An
  infinitude theorem gives no search bound.
- **CONDITIONAL:** H1 plus "every [x, x + (log x)^A] contains a non-x^c-smooth integer" ⇒
  polynomial-time primes.

---

## 4. Complexity (task 1(c))

### 4.1 Straight-line programs and τ (both directions)

**⇒ factoring.**
- Uniform polylog SLPs for k! ⇒ factoring ∈ FP, by gcd binary search (Lipton 2009; codex §5.1,
  PROVED).
- Lipton 1994, in the form of Koiran–Portier–Rojas Thm 1.14(II): polynomials with R > exp(τ^ε)
  integer roots of bounded height ⇒ non-uniform, i.o., positive-fraction factoring circuits.
- The τ-conjecture ⇒ P_C ≠ NP_C (BCSS Thm 3).
- Bürgisser 2009, Main Thm 1.2: polynomial τ(Per_n) ⇒ polylog τ(k!).

**⇐ factoring.** No implication is established. H1 supplies no polynomial with many integer
roots and no SLP for k!. We exhibit no contradiction between H1 and τ (codex §5.1).

**Mechanism-specific obligation (Ward–Maestro §7.1).**
- A uniform division-free SLP evaluating the rising factorial P_k with exponent γ factors in
  N^{γ/2}.
- The τ-conjecture with exponent c forces γ ≥ 1/c.
- So a *P_k-evaluation* algorithm below N^{1/6} (γ < 1/3) would refute τ for every c ≤ 3.
- This applies only to algorithms that evaluate that polynomial.

### 4.2 Möbius, bounded depth, Sarnak

**KNOWN.**
- Allender–Saks–Shparlinski, JCSS 62 (2001), Thm 1: MAJORITY AC⁰-reduces to primality,
  squarefreeness and gcd. So polynomial-size AC⁰[p] cannot compute them exactly.
- Green, CPC 21 (2012), Thm 1: μ ⊥ AC⁰(d) for d ≤ c log n/log log n.
- Bourgain, Israel J. Math. 197 (2013): the Walsh spectrum of μ is spread.

These concern exact or correlated *polynomial-size* circuit families on all inputs. They do not
constrain sub-1/6-time algorithms, which have circuits of exponential size, or promise
algorithms on S_n, which has density Θ(n^{−2}).

**Kalai's TC⁰ conjecture** (μ orthogonal to TC⁰) implies that factoring is not computable by
polynomial-size TC⁰ circuits on all inputs (§2.5).

**Sarnak.**
- H1 gives μ ∈ P. No contradiction with zero-entropy Möbius disjointness is exhibited:
  P-computable sequences can have positive entropy (codex §5.2).
- Sarnak's Mahler Lectures 2011 slides ("Möbius randomness and dynamics") say: "Perhaps ξ ∈ P ⇒
  μ is orthogonal to ξ? I don't believe so since I believe factoring and μ itself is in P." He
  then asks for a bounded ξ ∈ P with (1/N)Σμ(n)ξ(n) → α ≠ 0.
- Under H1, ξ = μ·1[μ ≠ 0] + 1[μ = 0] answers this.

### 4.3 Generic models, and the combined barriers (task 2)

**KNOWN.**
- Aggarwal–Maurer (Eurocrypt 2009) and Jager–Schwenk (J. Cryptol. 2013) give *equivalences*:
  generic-ring RSA and generic Jacobi symbols each reduce to factoring.
- The Jacobi symbol is easy by reciprocity, the canonical example of a representation-dependent
  operation outside the generic-ring model.
- Damgård–Koprowski 2002 is a generic-group floor for hidden-order root extraction. An explicit
  modulus is outside its model, so it gives no factoring floor.

**Lemma KE (key enlargement; codex §8.1; re-checked).**
- *Setting.* Keep T1's fixed generators, M | K, M, K ≤ P/2, the normalised D/|A| ≤ C₁M², and
  T1's moment budgets over all n-bit N. Replace K ≤ C₂M² by κ := max_ι max(1, K_ι/M_ι²).
- *Conclusion.* The covered pairs number O(P^{3/2}√κ X^{3/2}).
- So, in any balanced prime box of fixed relative size with ≫ P²/n² pairs, the covered fraction
  is ≤ C n²√κ X^{3/2}N^{−1/4}.
- *Use.* An N^{1/6−δ} design of this kind covering a constant fraction needs κ ≥ N^{3δ}/n⁴.

**Theorem CS (combined signature barrier; codex §8.2; re-checked, premises restored).**

*Model.*
- X is the enlarged-key T1 model of KE.
- Y is a fixed, N-independent list of certificates with prime-local meaning: no N-dependent
  exponents or curves. Y splits pq only when σ(p) ≠ σ(q). For richer labels, equality of labels
  must rule out every Y-split.
- I is a prime interval with Π ≍ P/n primes, and I × I lies in one n-bit balanced range.
- R is the number of distinct signatures on I, and a_j are the class sizes.
- **Aggregate closure.** All certified X-successes, after any Y-transcript, lie in **one**
  KE-admissible design with the stated global budgets.
- There are no cross-stage correlations and no useful accidental collisions.

*Conclusion.* Let ε be the failure fraction on distinct unordered pairs from I, and
c_eq = Σ_j a_j(a_j − 1)/(Π(Π − 1)) = (Π·Coll − 1)/(Π − 1), where Coll = Σ(a_j/Π)² ≥ 1/R. Then

  C n²√κ X^{3/2} N^{−1/4} ≥ c_eq − ε.

*Proof.* Equal-signature pairs are never split by Y. Closure puts their coverage inside one KE
design. ∎

- *Non-vacuous range.* R = o(Π) and ε = o(1/R). Then X ≫ N^{1/6}κ^{−1/3}R^{−2/3}n^{−4/3}, and
  X = N^{1/6−δ} forces κR² ≫ N^{3δ}/n⁴.
- *Closure matters* (Opus critic; codex). Equal-signature pairs from different classes can still
  get different Y-transcripts: a certificate made of exactly one class's primes returns N on that
  class and 1 elsewhere. An X that adapts to the transcript can run R designs. Union bounds then
  give only X ≫ N^{1/6}κ^{−1/3}R^{−4/3}/polylog.
- Ward–Maestro Conj 33.1 (distinct labels on all primes) has c_eq = 0. It is consistent with CS.

**Theorem CB (local collision search + (Z/N)^* smooth-exponent tests).**
- **CONDITIONAL on F_iv.**
- It is a routine combination of T1 with the small-order-removal mechanism of Ward–Maestro
  Prop 35.3(ii). We extend that mechanism to many bases, ± tests and adaptive smooth exponents.
  No prior formulation of the combination was found by us.

**Hypothesis F_iv.**
- *Statement.* For fixed 0 < a < b there is c(a, b) > 0 with
  #{p ∈ [aX, bX] prime : P⁺(p − 1) > p^{2/3}} ≥ c(a, b)X/log X for large X.
- *Source.* Fouvry, Invent. Math. 79 (1985) 383–407, gives the cumulative count
  #{p ≤ x : P⁺(p − 1) > x^{0.6687}} ≫ x/log x. Its localisation to every fixed interval is
  **UNVERIFIED** in the source.
- Ward–Maestro Prop 35.3 uses the range (x^{1−ε}, x] and does not supply F_iv.
- *Unconditional variant (sketch, unreviewed).* Pigeonholing Fouvry's count over dyadic, then
  2^{1/4}-adic, blocks gives a good block [z, 2^{1/4}z] at a scale z ≍ x. Put both p and q in
  it. Then CB holds for a set of n with bounded gaps.
- E1 measures the density of P⁺(p − 1) > B₂ as 0.33–0.36 at 20–32 bits.

**Model.**
- **(X)** T1-local certified collision search, with the **per-input** budget #τ ≤ X and
  Σ|S_τ| ≤ X² on every input.
  - Classes, windows and activations depend only on N mod K (K ≤ C₂M²) and on the **answer
    bits** of Y₁.
  - They do not depend on Y₁'s queries, on intermediate residues, on timing, or on group values.
  - Values a^E mod N may serve as T1 handles, as in D′'s γ = β^{∏M²}.
- **(Y₁)** Tests gcd(F, N) with F = ∏_j (a_j^{E_j} − 1) or ∏_j (a_j^{E_j} + 1), where:
  - every a_j lies in a fixed base set A with |a| ≥ 2 and Σ_{a∈A} log|a| ≤ P^{1/3}/n²;
  - the E_j are positive and B₂-smooth, with 2 ≤ B₂ ≤ (P/4)^{2/3};
  - the choices may be adaptive, and Y₁ has no cost bound.
- **Success** means a certified T1 hit, or a proper gcd from a Y₁ test of exactly this form.
  - *Excluded:* order collisions gcd(a^E − a^{E′}, N), which fire when ord_p(a) | E − E′, and
    other accidental hits, as in T1.
  - With a cost-free Y₁, the exclusion is essential: otherwise Y₁ would contain order finding.
    So CB is **not** about every use of (Z/N)^*.

**Conclusion.** There are c₀ > 0 and C such that every such algorithm fails on
≥ (c₀ − C n² X^{3/2} N^{−1/4})|S_n| inputs of S_n. Success on a fraction > 1 − c₀/2 forces
X ≥ cN^{1/6}/n^{4/3}.

*Proof.*
1. **Box.** Let I_p = [0.6P, 0.65P] and I_q = [0.84P, 1.19P].
   - Then pq ∈ [0.504P², 0.7735P²] and q/p ≤ 119/60 < 2, so every such pair lies in S_n.
   - F_iv on these two disjoint intervals gives ≫ P²/n² pairs with p, q ∈ F, where
     F = {ℓ : P⁺(ℓ − 1) > ℓ^{2/3}}. That is ≥ c₂|S_n|.
2. **Small orders.**
   - Fix a ∈ A and put U = (2P)^{1/3}. The primes r ≥ P/2 with r ∤ a and ord_r(a) ≤ U all
     divide ∏_{k≤U}(a^k − 1). So there are ≤ Σ_{k≤U} k log|a|/log(P/2) ≪ U² log|a|/log P of
     them.
   - Summed over A, this is O(P/(n² log P)) = o(P/log P).
   - The primes dividing some a ∈ A number ≤ Σ log|a|/log(P/2); these are negligible too.
   - Remove all of them, and let G be the remaining pairs, with |G| ≥ c₀|S_n|.
3. **Y₁ is silent on G.** Take p ∈ G and ℓ = P⁺(p − 1) > p^{2/3}.
   - ℓ² > p, so ℓ ∥ p − 1, and (p − 1)/ℓ < p^{1/3} ≤ U.
   - ord_p(a) > U, so ℓ | ord_p(a).
   - ℓ > (0.6P)^{2/3} > B₂, so ℓ ∤ E and ℓ ∤ 2E. Hence a^E ≢ ±1 (mod p).
   - The same holds for q, so every F is a unit mod N and every answer is 1.
4. **Closure.**
   - Selectors depend only on N mod K and the answer bits. So all inputs with key κ and
     all-one answers run one design D*(κ). Set D*(κ) = ∅ if no input has key κ and all-one
     answers.
   - D*(κ) is actually executed on such an input, so it obeys the per-input budget. Hence D* is
     a T1 design with T1's budgets over all n-bit N.
   - On G, the algorithm's certified coverage equals that of D*. T1 counts coverage through
     (2.1) alone, so handles are irrelevant.
   - T1 bounds D*'s total covered pairs in S_n by C n² X^{3/2} N^{−1/4}|S_n|. Intersect with G.
     ∎

**Remarks.**
- (a) *Averaged budgets.* Under T1's averaged-budget convention the argument fails (Opus
  critic). An algorithm could stop cheaply on inputs that Y₁ splits and spend ≫ X on the rest.
  Closing this needs a rough-numbers-in-progressions bound. So CB does **not** cover "average
  time" H4a.
- (b) *Scope.* CB shows that, within this model, Pollard-type exposure in (Z/N)^* cannot be
  completed to an almost-all sub-1/6 algorithm. Ward–Maestro Prop 35.3 already gives
  positive-proportion success below N^{1/6}.
  - Escaping CB while staying in (Z/N)^* is possible by enlarging B₂ or the base portfolio, by
    admitting order collisions or other identities, or by non-local selectors.
  - "Other groups" is therefore one escape, not a forced one.
- (c) *CB⁺ (SPECULATIVE).* Extending CB to the tori T_d, d ∈ {2, 3, 4, 6}, needs three things:
  - a positive-density joint residual for all allowed groups;
  - control of element orders inside each torus (the 41·43 example in §8 shows why);
  - the same closure.
  The direct progression-counting route for p ≡ 1 (mod ℓ₁), p ≡ −1 (mod ℓ₂) with
  ℓ₁ℓ₂ > p^{4/3} is beyond level 1. That diagnoses one strategy; it is not an impossibility.

**Prop LW (WM Conj 33.1 implies a least-exposure bound; PROVED, by pigeonhole).**
- *Assume* WM Conj 33.1 for fixed u and C: B₁ = x^{1/u}, B₂ = B₁², and the first C_{u,C} log x
  Suyama parameters give distinct label vectors on the primes of [x, Cx].
- *Conclusion.* All but at most one prime p in each such range have some σ among the first
  C_{u,C} log x with λ_σ(p) ≠ 0 (WM Def. 32.1). That is, a construction or singularity label,
  or ord_{E(F_p)}([L]P) ≤ B₂ for the stage-1 residual point [L]P, L = lcm{ℓ^e ≤ B₁}.
- So the least exposing index is O(log p), except for one prime per range.

*Proof.* Label 0 means unexposed, and distinctness allows at most one all-zero vector. ∎

- *Strength, not refutation.*
  - LW is necessary for Conj 33.1, not equivalent to it.
  - It measures the uniform hitting theorem the conjecture needs: a log-scale least-witness
    bound for a non-character event. No GRH-type theorem of that kind is known for an explicit
    curve family.
  - The QNR comparison is only an analogy. Least non-residues are O(log p) for no fixed
    constant (Graham–Ringrose, unconditional, i.o.; Montgomery under GRH). But the cause is
    multiplicativity, and an i.o. result does not refute "all but one per range". Under
    independence, a maximum waiting time ≍ log p is plausible.
  - Polynomial-length lists (any fixed-list route needs only w(p) ≤ N^{1/6−δ−1/(2u)}) are in a
    Burgess-type regime that LW does not touch.

---

## 5. Escape-signature table (task 2)

Column guide:
- **Outside-model examples**: known methods that are not in the model. They are not necessarily
  fast for the model's task; each keeps its own rigor and success population.
- "Promise" means the method needs auxiliary information or a special family.
- Full model statements, with the codex row audit, are in codex/polyimpl_1.out §7.1 and
  codex/polyimpl_2.out item 6.

| Barrier | Model assumptions | What a faster in-family method must change | Outside-model examples |
|---|---|---|---|
| WM Cor 3.3 / Thm 3.4 / Thm 9.4 / Thm 2.8 | oblivious coverage of every real p by Lehman cells; materialised giants (3.4); fixed-degree forms (9.4); one enumerated hypothesis modulus (2.8) | free giants (3.4 only); several moduli plus true-factor-only coverage; non-cell tests | D′/AV (distributional, ≥ 1/6); PS (1/4); p−1⁺; ECM, QS/NFS (heuristic); LP (randomized); Shor; Copp (promise) |
| WM Thm 3.6 / multimod Thm A | bounded degree; identity-certified **full-candidate** coverage; #τ ≤ X, ΣW ≤ X² | true-factor-only coverage, or leave sumset search | as above; D′ attains, it does not escape |
| WM Thm 4.4 | explicit difference covers; start-statistic envelope D_max ≤ r^η | different family or envelope; implicit product evaluation (WM 7.2) | none below 1/5 in this architecture |
| WM Prop 4.5 | divisibility tests by polynomially bounded integers | a different task representation | PS attains 1/4 |
| WM Thm 10.1 | pieces priced ≥ √length separately | shared computation | none supplied |
| multimod A′ / B | A′: insufficient degree-1 coverage. B: J moduli, height conventions | N^{Ω(1)} moduli | D′ (the many-moduli ingredient only) |
| below16 T1 / T2 / T3 (three separate locality classes) | T1: key N mod K ≤ CM², arbitrary windows. T2: plus ≤ PK size boxes. T3: algebraic windows keyed by N mod M. Moment budgets over all n-bit N | KE: κ ≥ N^{3δ}/n⁴ for T1; T2 and T3 by their own terms | p−1⁺; specialfam / Copp (promise); LP (randomized) |
| Theorem CB (CONDITIONAL on F_iv) | T1 (per-input) + Y₁ smooth-exponent tests, as stated | larger B₂ or base set; order collisions; other groups; relations; non-local selection | ECM (heuristic); LP; QS/NFS; Shor |
| Theorem CS | T1 (KE) + fixed certificate list; aggregate closure | κR² ≥ N^{3δ}/n⁴; adaptivity | WM Conj 33.1 (conjectural); WM Prop 33.2 (randomized lists under exposure-probability hypotheses) |
| coverage38 Thm U | N-independent (weight, modulus) groups, all charged; WIN where established | coverage with λ < J^{1/4} for every pair (an achievement inside the model); adaptive groups or other charging | none |
| coverage38 Prop L; dchworst B1 | phase-discarding / error-budget arguments | keep phases or exploit arithmetic correlations | none (proof-method rows) |
| dchworst B2 | single-modulus pigeonhole | a worst-case union-over-moduli theorem | D′ (average only) |
| dchworst B3 | proof valid for all coprime pairs | use primality; prove U(x, Y^{1/8}); or a covered-or-easy dichotomy | none |
| sieve2d Prop 5 | fixed J; only ≤ J-prime event statistics | growing-order or structural information | none |
| unbalanced Prop O′ | 1-D Lehman-linear full / constant-fraction candidate coverage | true-factor coverage, or leave 𝔏 | PS *attains* the √p branch; ECM, NFS, p−1⁺ change the task |
| Track C O1–O4 | specific unit-predicate representations (finite-algebra powers, power maps) | nonlinear or non-faithful readouts, or a different task | PS/BGS at √ cost; Tal needs p (different input) |
| Track A A1 | faithful linear observables on dense F_N coordinates | avoid that representation | none supplied (no new fast witness) |
| Track E information floor | per-modulus information | non-per-modulus computation | Copp needs promised coset-shaped n/4 bits |
| eps battery | finite negative experiment, not a theorem | — | CFRAC/QS relations reveal ε at L[1/2] |
| nfsfruit Lemma M | an **open sufficient** fruitfulness condition, not a floor | — | Dixon / Theorem DD use the fibre symmetry (randomized) |
| lemmaM B1 | even-|S| two-form degree-1 sieves | other sieve shapes | QS; hidden-root multipliers |
| lemmaM (c) | one-level HECM-priced descent: L[4/9] | multilevel descent; other smoothness tests | LP (L[1/2], randomized) |
| hypD G / H | counting template; fixed-k shape test | new averaged estimates (FM, PC) | none rigorous below L[1/2, 1] |
| generic ring (AM09, JS13) | equivalences, not floors | — | Jacobi via reciprocity |
| generic group (Shoup; DK02) | different tasks and representations; no factoring floor | — | — |
| AC⁰[p] (ASS; Green; Bourgain) | exact or correlated **polynomial-size** circuits on all inputs | irrelevant to exponential-time or promise algorithms | — |
| τ-conjecture (WM 7.1) | uniform SLP evaluating P_k | do not evaluate P_k | ECM, QS |
| Theorem DR | a conditional upper bound, not a floor | — | conditional deterministic L[1/2] |

**A defensible escape signature for an unconditional deterministic N^{1/6−δ} algorithm**
(almost all, or worst case):
- **S1.** It is not a charged, bounded-degree, identity-certified full-candidate cover (multimod
  A / WM 3.6).
- **S2.** If an enlarged-key T1 component alone covers a constant fraction, then
  κ ≫ N^{3δ}/n⁴ (KE). General non-local selectors are outside KE.
- **S3.** Within repaired CB/CS, the complementary mechanism must escape one of their explicit
  assumptions. Examples:
  - other groups with a hitting theorem (ecmavg HIT);
  - relations with a deterministic supply and fruitfulness;
  - larger B₂ or order collisions with their own analysis.
  Signature richness (CS) and least witnesses (LW) matter only for fixed-list, worst-case
  splitting.
- The earlier draft's S4 (super-constant depth) is **withdrawn** as a necessity. It concerns
  polynomial-size all-input circuits only.

**Witnesses.**
- No known classical algorithm is unconditional, uniform, deterministic and proved below N^{1/6}
  on almost all of S_n, or on every input. This is a literature finding, not a completeness
  theorem.
- Ingredients exist separately:
  - p−1⁺ and specialfam (positive proportion or promise);
  - LP and Dixon (rigorous, randomized);
  - ECM and det-QS (heuristic);
  - Theorem DR (conditional on IW);
  - WM Conj 33.1 (conditional; N^{1/(2u)} for bounded ratio).

---

## 6. Blockers (task 3(i))

1. **Least witness (Prop LW).** Conj 33.1's *log-length* lists need an O(log p) least-exposure
   bound for an explicit curve family. No theorem of that kind is known even under GRH. This
   blocks that route; it does not block polynomial-length lists.
2. **Joint residuals (CB⁺).** A barrier for torus portfolios needs joint large prime factors.
   Conversely, a proof that a torus portfolio exposes almost all p would need a joint
   smoothness statement of comparable depth. E1 (post hoc) shows no decay at 20–32 bits;
   small-size data cannot settle the asymptotics.
3. **Coverage theorems (H4w via Theorem R).** dchworst B1–B3 and coverage38 U/L, unchanged.
4. **τ-type.** A P_k-evaluation algorithm below N^{1/6} refutes τ for every c ≤ 3 (§4.1).
5. **Lower bounds.** By Theorem DR, an unrestricted deterministic N^{Ω(1)} lower bound would
   refute the IW hypothesis. The repository's model floors cannot become such lower bounds
   unless that hypothesis fails.
6. **Not blockers** on the evidence checked: prime construction (§3), Sarnak, τ, μ ⊥ AC⁰
   (§4.2), analytic number theory (§2.1).

## 7. New directions (task 3(iii))

**D1. Fixed-u deterministic QS on almost all of S_n.** Feasibility low–medium; payoff high.
- Below 1/6 needs y = N^c with c < 1/(6ω), so the smoothness densities are constants
  (below16 §6.3). The remaining obligations are fruitfulness and deterministic elimination.
- **Transfer lemma (PROVED, trivial; codex confirms).** If at most ε_n·2^{n−1} n-bit integers
  are bad, the bad fraction of S_n is O(ε_n n²). So ε_n = o(n^{−2}) suffices for an almost-all
  statement. For a mean-time statement with an N^{1/5} fallback one needs instead
  η ≤ N^{−1/30−δ}.
- **Coordinates.**
  - Put m = ⌈√N⌉, r = m² − N, v = r + 2mt + t² and w = r + 2mt′ + t′².
  - The map (m, r) ↦ (v, w) is affine with determinant 2(t − t′). With k = t′ − t:
    w − v ≡ t′² − t² (mod 2k), and m = (w − v − (t′² − t²))/(2k).
  - The domain is an affine lattice coset intersected with a **thin coupled region**: for fixed
    v, w lies in a window of relative width O(k/T²).
  - So the pair sum is a two-variable smoothness correlation, not a one-variable count of smooth
    numbers in an AP. (The earlier draft's appeal to Harper was also wrong: Harper's theorem
    needs log x/log q → ∞, which fails at q ≍ T ≍ y.)
- **First lemma (target).** An all-integer lower-tail bound. Find a count Z_N with mean ≍ T
  and a fixed gap above the relation count, with Var(Z_N) = o(T²/n²) after justified
  stratification by N mod Q, or a direct lower-tail bound of strength o(n^{−2}).
  - Near-diagonal pairs |t − t′| ≤ K contribute only O(K/T) relatively, so they can be
    discarded.
- **Kill criterion.** The thin-region correlation for k in a positive-density range of pairs
  cannot be bounded with o(n^{−2}) relative error by any available smooth-number technique:
  short-interval, AP or bilinear.
- Fruitfulness then reduces to one rank test (codex §9.4), but its *frequency* still has to be
  proved.

**D3. Non-local selector audit (codex Rank 1).** High feasibility for audits.
- For enlarged-key T1, KE turns a proposal into κ ≥ N^{3δ}/n⁴.
- Kill the proposal if the key compresses back to M²N^{o(1)}, if placement needs factor
  information, or if construction cost cancels the gain.

**D2. Hitting outside (Z/N)^*.** Feasibility low.
- Separate the constructive question (a hitting or separation theorem for an explicit torus or
  curve family) from the obstructive one (CB⁺).
- First lemma: positive density of p with P⁺(p ± 1) both > p^{1/2+η}, the CB⁺ input at budget
  N^{1/8}.
- Kill criterion: only for strategies whose moduli product exceeds the level of distribution.
  That is a limitation of the strategy, not an impossibility.

**D5. Targeted derandomization (Theorem DR, unconditionally, for one test).** Feasibility very
low.
- To get below 1/6 one needs seed length (1/6 − δ)n, with N^{o(1)} cost per seed and a proved
  hitting property.
- Dixon's test with structured seeds r = ⌈√(kN)⌉ + t is det-QS (D1).
- Lenstra–Pomerance's success predicate is more complex than "h(−kN) is smooth"; it is not
  reducible to D1.

**D4. Promise-preserving partial-information reductions (codex Rank 2).**
- Target: Jacobi-(+1) QR ≤ ε on constructed semiprimes, with the promise and query lengths
  proved.
- Coprime masks are one failure mode (codex §9.3); they are not the only one.

**D6. Residual-tail portfolios (codex Rank 4).**
- Keep this as a quantitative hypothesis, η ≤ N^{−1/30−δ}.
- E1's fixed five-group proxy leaves about 1% of pairs at the sampled scales. That is evidence
  against a fixed cyclotomic portfolio at these scales, not an asymptotic obstruction, and it
  says nothing about point-order or growing portfolios.

No direction offers a cheap decisive experiment beyond E1. D1, D2 and D5 need proofs, and D3
and D4 audit specific proposals.

---

## 8. Experiment E1 (pre-registered, `notes/polyimpl_prereg.md`, commit `1a9c3de`)

**Event.**
- "Exposed by d" means P⁺(Φ_d(p)) ≤ B₂, i.e. **whole-order smoothness** of the torus T_d(F_p).
- The proxy is generous relative to whole-order stage-1/2 criteria, because it ignores stage-1
  cost.
- It is **not** an upper bound on all exposures: an element may have smooth order inside a
  non-smooth group.
  - Codex's example: at B₂ = 3, all five Φ_d values of 41 and of 43 are non-smooth, yet
    gcd(9^4 − 1, 41·43) = 41.
  - CB removes small orders for fixed integer bases. E1 has no such step for the other tori.

**Sampling.** 20 000 uniform primes per dyadic block k, without replacement. B₂ = ⌊2^{2k/3}⌋.
10 000 disjoint matched pairs per block (balanced, but not uniform on one S_n).

| k | B₂ | r₁⁽¹⁾ | r₁⁽¹,²⁾ | r₁ (all five d) | r₂ (pairs) | unexposed by d = 3, 4, 6 |
|---|---|---|---|---|---|---|
| 20 | 10 321 | 0.3304 | 0.0951 | 0.0829 | 0.0069 | 0.961, 0.958, 0.965 |
| 26 | 165 140 | 0.3497 | 0.1148 | 0.1016 | 0.0108 | 0.964, 0.959, 0.963 |
| 32 | 2 642 245 | 0.3575 | 0.1226 | 0.1086 | 0.0113 | 0.958, 0.955, 0.959 |
| Dickman model | | 0.405 | 0.164 | 0.141 | 0.020 | 0.951 |

The descriptive binomial errors are ≈ 0.002 for r₁ and ≈ 0.001 for r₂. They ignore the
finite-population correction and the matching.

**Frozen checks.**
- P1 holds at every k.
- **P2 fails at k = 20 and k = 26** (0.0951, 0.1148 < 0.12) and holds at 32.
- P3 fails at 20 (0.0829 < 0.09) and holds at 26 and 32.
- P4 fails at 20 (0.0069 < 0.008) and holds at 26 and 32.
- P5 fails: r₁(32)/r₁(20) = 1.310 > 1.3, by **growth**.
- **Pre-registered verdict: "mixed"** (neither branch of the rule applies).

**Post hoc (labelled).**
- The residuals increase with k toward the Dickman-model values, which are heuristics for these
  prime-polynomial values, not theorems.
- p − 1 and p + 1 are close to independent: 0.1226 against 0.3575 × 0.3549 = 0.1269 at k = 32.
- The data support only "no decreasing trend at these sizes". They give no support to CB⁺ as a
  theorem, and no refutation of D6 beyond the fixed five-group proxy.

**Reproduction.**
- `PYTHONPATH=. .venv/bin/python -m src.polyimpl.residual run` (≈ 2 s)
- `PYTHONPATH=. .venv/bin/python -m src.polyimpl.primefind`
- `.venv/bin/python -m pytest -q tests/test_polyimpl.py`

## 9. Reviews

**Independent Opus critic** (read-only; all findings applied).
- DR: correct with fixes. Folklore; wording of the contrapositive.
- CB: gap plus overclaim.
  - The arithmetic, box and small-order count were all verified.
  - Fixes: per-input budget only (the averaged convention fails); selectors depend only on N mod
    K and answer bits; the exact Y₁ form, with order collisions excluded; the false WM §35
    attribution for F_iv removed; downgraded to a routine combination with WM 35.3(ii) cited.
  - The critic suggested the unconditional bounded-gap variant.
- CS: gap. Codex's premises had been dropped. Restored: aggregate closure, the exact
  distinct-pair fraction, the ε = o(1/R) range, and the adaptive R^{−4/3} caveat.
- KE: correct with fixes (retained conditions).
- PF: correct with fixes. "Needed to beat 1/2"; the literature-completeness claim scoped.
- LW: correct; blocker overclaimed. Polynomial lists are unaffected; Graham–Ringrose; the i.o.
  versus all-but-one mismatch.
- Catalog: Kalai item correct. Overclaims fixed (S3c, base count, CONDITIONAL row label,
  balanced-only N^{1/10}).
- E1: §8 honest; §§0, 6, 7 had overclaimed. Fixed.

**GPT-6.1-sol adversarial referee** (`codex/polyimpl_2.out`; all valid fixes applied).
- DR: correct with fixes. The almost-everywhere hypothesis made explicit, the explicit clock,
  and the quantifiers of the contrapositive.
- CB: correct with fixes under an explicit closure hypothesis, now in the model. Bases 0, ±1
  excluded; primes dividing bases removed; F_iv kept conditional.
- KE/CS: correct with fixes. Fixed-relative-size box; c_eq; aggregate closure; prime-local
  certificates.
- PF: correct with fixes.
  - The Sander row was wrong: e = 1/(1 + cβ²), not ≥ 1, and factoring *is* needed there.
  - The quantifier is "every interval".
  - The output range is z < p ≤ 2z^{1/θ}.
  - Merikoski and Li confirmed.
- LW: correct with fixes. The label is the order of the stage-1 residual point; one exception
  per range; the QNR comparison is an analogy only.
- Escape table: overclaims fixed row by row:
  - O′: PS attains, not escapes;
  - AC⁰ row scope;
  - generic-group row;
  - nfsfruit row split;
  - eps row a finite experiment;
  - S4 withdrawn.
- Catalog:
  - Żrałek: totient chain;
  - Dieulefait–Urroz: Thm 4/5 and the twist condition;
  - Bach–Charles: verified, density-one, a(n) and a(n²);
  - class groups: scope;
  - H5 classification limited to the menu;
  - Sarnak source: Mahler Lectures 2011;
  - Helfgott–Thompson verified;
  - NR has a factoring-based PRF;
  - H3 quantifiers;
  - "no known collapse".
- E1: P2 fails at 26 bits as well (corrected). Event renamed whole-order smoothness. CB⁺ and D6
  claims withdrawn to "no decreasing trend at these sizes".
- D1: thin-region correlation; Harper's missing hypothesis; kill criterion replaced by the
  lower-tail target. Codex's own ranking: D3, D1, D5, D2, D4, D6.
