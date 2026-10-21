# factorbits: prime-exponent addition, hardcore bits of factors, and the information floor (phase 21)

Agent `factorbits`, 2026-10-21.
- Pre-registration: `notes/factorbits_prereg.md` (commit `3c9fe3f`, before E1 ran).
- Code: `src/factorbits/{equiv,harness,pure,gl,run_experiments}.py`. Tests: `tests/test_factorbits.py`.
- Raw data: `experiments/raw/factorbits_e1.jsonl`.
- GPT-6.1-sol: blind round `codex/factorbits_1.{md,out}`; adversarial round `codex/factorbits_2.{md,out}` (§8).
- Two read-only Opus critics (§9).
- Labels:
  - **KNOWN** (with reference);
  - **PROVED** (proof here; routine unless stated; "no prior formulation found by us" is the strongest novelty claim);
  - **HEURISTIC**;
  - **CONDITIONAL**;
  - **OPEN**;
  - **experimentally supported**;
  - **UNVERIFIED**.

Notation:
- n = bit length of N.
- N = pq is a balanced semiprime with p < q.
- A *vector* is the prime-exponent representation {prime: exponent} of a positive integer.
- FACT is complete factorisation of arbitrary integers. "Splitting" means finding a nontrivial divisor.

---

## 0. Verdict (plain language)

**1. The equivalence is real and cheap.**
- Fast addition in the prime base is the same problem as factoring.
  - An addition oracle factors N in at most popcount(N) − 1 calls. A successor oracle x ↦ x + 1 is enough.
  - Heuristically it needs only one call. Write N = u + v with v prime (primality is easy) and u small enough to
    trial-divide.
  - Conversely, one factoring call adds.
- The Zech-logarithm analogue in F_q is exact: discrete logs need ≤ ⌊log₂p⌋ + 1 Zech calls, and Zech needs one
  discrete-log call.
- Partial forms:
  - Every output from which the support of a+b can be read is as hard as factoring. This includes exponents mod m,
    up to a quasi-polynomial overhead.
  - Sums of many smooth terms (one call: the binary expansion) are as hard.
  - Easy forms: the smooth part, a smooth-promised result, residues of a+b, linear functionals that are Dirichlet
    characters or log, and approximation to relative error 1/poly.
  - Open forms: ω, Ω mod 2, μ, squarefreeness, the numerical radical, and two smooth summands.
  - The equivalence needs a size convention. With exponents written in binary, "addition" would factor Fermat
    numbers from O(k)-bit inputs.

**2. No natural single bit of the factors is known to be as hard as factoring.**
- Every known hardcore bit of factoring is a bit of something else, or uses auxiliary randomness:
  - the Rabin/RSA message bits of ACGS and Håstad–Näslund;
  - the Goldreich–Levin parity ⟨p, r⟩ with a public random r.
- Two artificial single-bit functions of the integer alone are factoring-hard:
  - Prop GLM (new here, average case, amplifies weak predictors);
  - codex's B(M) (worst case).
  - Both work only because the reduction can smuggle fresh randomness or a bit index into the cofactor.
- For every natural fixed quantity (bit i of p, p mod 4, ε = (−1/p), bit 1 of the smaller factor, ω, μ,
  squarefreeness, φ(N) mod ℓ, fixed-base QR), the status is OPEN. No reduction is known in either direction.
- Things that are as hard as factoring:
  - ⌊p/2^k⌋ for k ≤ n/4 + O(log n);
  - p mod M for log₂M ≥ n/4 + O(log n);
  - φ(N) mod 2^k for k ≥ n/4 + v₂(q−p);
  - "N has a factor ≤ x" with x as input.
  - All of these go through Coppersmith or binary search.

**3. The information floor is a theorem for N-oblivious algebraic algorithms.**
- Theorem PT (PROVED; elementary; no prior formulation found by us) covers *pure* generic-ring decision trees: integer
  straight-line programs whose constants do not depend on N, with zero tests mod N, of any size, adaptive. It covers
  Pollard rho, p−1, p+1, ECM, Pollard–Strassen and Lipton-type SLPs.
  - Such a tree either splits N or ends in a leaf that is a function of a *partition class* of p, and that class is
    also the class of q.
  - Consequently, for iid factors, any advantage δ in predicting *anything* about the factors forces the same tree,
    with gcds, to split N with probability ≥ δ/2.
- With side information N mod m, the only non-split leak is the local fiber {p mod m, q mod m}. Even that is paid for:
  a pure tree that decides ε on N ≡ 1 (mod 4) with advantage δ splits N ≡ 1 or N ≡ 3 (mod 4) with probability
  ≥ δ²/2.
- Floor-type statements also hold for:
  - functions of N mod m (trackE F1);
  - abelian reciprocity symbols (nonsmooth §5.4);
  - multiplicative automatic sequences (Konieczny–Lemańczyk–Müllner).
- OPEN: finite automata and AC⁰ on N's bits, and every non-pure, N-dependent class (Fermat/MR exponents, continued
  fractions, forms of discriminant −4N, sieving, lattices). All known sub-exponential and sub-N^{1/4} methods live
  there.
- E1 confirmed the lemma and both bounds on real executions: 0 mismatches in 681 990 instrumented runs, 0 violations.

**4. Leverage: a weak predictor for a fixed bit of p is not known to amplify, and the reason is structural.**
- ACGS and Blum–Micali amplify because the secret lives in a group that acts on instances: x ↦ ax for RSA, x ↦ ax + b
  for DL. For factoring the secret is the *support* of the exponent vector, and the computable instance maps only add
  known vectors.
- Lemma PSR (PROVED): no polynomial map N ↦ h(N) carries the factor p to any g(p) other than c·p^j.
  - Divisibility sequences (2^N − 1 ∋ 2^p − 1) do move p, but at exponential size.
- So the minimal known sufficient leakage is:
  - ≈ n/4 + O(log n) bits of p in a coset shape (Coppersmith; trackE measured n/4 + n/(4d) − 2);
  - or εn arbitrary yes/no queries (Maurer 1995, conditional).
- O(log n) bits are worthless.
- A single fixed bit would need a self-reduction that nobody has; we found none and rule out the polynomial ones.

---

## 1. The equivalence (question 1)

### 1.1 Size convention (PROVED; codex agrees)
The vector of a has size s = Σ(log p_i + log e_i), but a itself has Σ e_i log p_i bits. For example (2: 2^k) has O(k)
bits and denotes a (2^k + 1)-bit integer. An addition oracle on succinct inputs would factor the Fermat number
2^{2^k} + 1 from an O(k)-bit query. So **ADD_n** is defined on inputs with log₂a, log₂b ≤ n, and every complexity is
measured in n. Equivalently, bound the weighted exponent sum Σ e_i log p_i.

### 1.2 Prop E1 (ADD ⇒ FACT). PROVED, with heuristic and conditional strengthenings

| Oracle | Calls to fully factor any N | Status |
|---|---|---|
| successor S(x) = vector of x+1 | exactly popcount(N) − 1 ≤ ⌊log₂N⌋ (Horner: doubling adds 1 to the exponent of 2 for free; queries only even x ≤ N) | PROVED |
| ADD, base 2^w windows | number of nonzero base-2^w digits after the leading one ≤ ⌈n/w⌉ − 1; with w = Θ(log n) that is O(n/log n), digits trial-divided | PROVED |
| ADD, short-interval prime trick | O(log log N) calls, randomized. Write X = u + v with v a random prime in [X − X^θ, X], θ > 7/12, found in expected O(log X) primality tests (Huxley's short-interval PNT). Then u ≤ X^θ is factored recursively, down to poly(n) size | PROVED (randomized, Las Vegas) |
| ADD, one call | Take v = the largest prime below N with u = N − v trial-factored. The prime gap is O(log²N) under Cramér's model, and ≤ N^{0.525} by Baker–Harman–Pintz | HEURISTIC (1 call); tested 10/10 at 64 bits |
| S, one call | find t with tN − 1 prime, call S on it, divide out t | HEURISTIC. CONDITIONAL on GRH: random t ≤ N·polylog(N) works with probability ≈ 1/ln, by the GRH form of PNT in progressions. Tested 10/10 at 64 bits |
| many-summand ADD of 2-smooth terms | 1 call: N = Σ_{b_i=1} 2^i | PROVED |

**Why the integers need fewer calls than F_q.** Primes are free in the prime base: they are dense, their vector is
themselves, and recognising them is polynomial (AKS). F_q has no comparable dense set of elements whose discrete log
is free and recognisable.

### 1.3 Prop E2 (FACT ⇒ ADD). PROVED
- With a complete-factorisation oracle, addition takes 1 call.
- With a splitting oracle and primality tests, recursive splitting takes exactly Ω(a+b) − 1 ≤ log₂(a+b) − 1 calls.
- ADD_n is therefore polynomial-time equivalent to FACT (complete factorisation of n-bit integers). It is not known to
  be equivalent to factoring *semiprimes*: the ladder queries arbitrary prefixes ⌊N/2^j⌋.

### 1.4 Prop Z (Zech ⇔ DL in F_q). PROVED (classical in substance)
Z(k) = log_g(1 + g^k), with a separate zero symbol.
- **DL ⇒ Zech:** 1 DL call after computing 1 + g^k.
- **Zech ⇒ DL in F_p:** log 2 = Z(0) costs 1 call. Horner on the bits of h then costs popcount(h) − 1 calls: doubling
  adds log 2, and adding 1 is one Zech call. The total is popcount(h) ≤ ⌊log₂p⌋ + 1 for h > 1. No partial value
  vanishes, because all prefixes are < p as integers.
- **Windowing** costs 2^w − 2 calls of digit-log precomputation and O(log p/w) Zech calls, so O(log p/log log p) overall.
- **F_{p^k}** in a basis 1, g, …, g^{k−1} (codex): δ + Σ_{c_j≠0}(wt₂(c_j) − 1) + (s − 1) ≤ k⌊log₂p⌋ calls, where s
  is the number of nonzero coordinates and δ = 1 if some coordinate is > 1. Partial sums of independent basis terms
  never vanish.

The dictionary:

| integers | F_q^* |
|---|---|
| free abelian group ⊕_p Z, the vector | Z/(q−1), log_g |
| FACT | DL |
| ADD | Zech |
| N mod m (positional floor) | Pohlig–Hellman projections (power-residue symbols) |
| secret = the support of the vector | secret = an element of a group |
| no random self-reduction known | x ↦ ax + b random self-reduction (Blum–Micali; Håstad–Näslund) |

The last row is the substance of §4.

### 1.5 Partial forms (PROVED unless marked; the oracle takes full vectors and returns partial output)

| Partial ADD | Status | Why |
|---|---|---|
| support of a+b (the list of primes) | **≡ FACT**, popcount − 1 calls | exponents by repeated division (`factor_via_support`) |
| exponents mod m, m > log₂(a+b) | **≡ FACT** | no exponent wraps (codex) |
| exponents mod a fixed m, zero residues dropped | **≡ FACT under quasi-polynomial reductions**: ≤ (2n)^{log_m n + 1} calls; O(n) calls when every prefix's m-th-power part has a trial-divisible root (typical) | the hidden part is an m-th power r^m with r ≤ x^{1/m}; recurse on r (`factor_via_expmod`, tested on adversarial square parts). A polynomial bound is OPEN |
| largest prime factor P⁺(a+b) only | for semiprime targets: 1 call HEURISTIC (the prime trick returns q). In general: OPEN as a polynomial reduction | the recursion on the unfactored cofactor x/P⁺^e is quasi-polynomial only under Dickman-type heuristics |
| numerical radical rad(a+b) | OPEN | squarefree a+b gives no information; computing rad is itself open (Adleman 1994 list, below) |
| ω, Ω mod 2 (Liouville), μ, squarefreeness, d(a+b) | OPEN in general; constant on semiprimes, hence useless for RSA moduli | Shallit–Shamir 1985 relate these to d(·) but give no factoring reduction |
| d(a+b) mod 2 | easy | it is 1 iff a+b is a square |
| B-smooth part, B ≤ poly(n) | easy | trial division of the ordinary sum |
| B-smooth part, larger B | HEURISTIC L_B[1/2, √2] by ECM | — |
| result promised y-smooth | easy for y ≤ poly; ECM in general | as above |
| linear functional ⟨v, w⟩, w = log | easy | log(a+b) |
| ⟨v, w⟩, w finitely supported | easy | repeated division |
| ⟨v, w⟩, w a Dirichlet-character index (w_p = ind χ(p)) | easy | it equals χ(a+b) |
| ⟨v, w⟩, w ≡ 1 mod 2 (Liouville λ) | OPEN | — |
| ⟨v, w⟩, w = an arbitrary set of primes, e.g. bit i of p | OPEN; on semiprimes it is bit_i(p) ⊕ bit_i(q) | — |
| two y-smooth summands, y ≤ poly(n) | OPEN | It contains Cunningham factoring (2^k ± 1). Using it needs a representation kN = a + b with a, b smooth. That is relation-finding: ≈ π(y) such representations factor N *without* the oracle by Dixon combination (HEURISTIC: fruitfulness can fail for structured producers, sicheck Prop B) |
| ordinary residue (a+b) mod m | easy | — |
| vector of the least residue (a+b) mod m | easy for m ≤ poly; ≡ FACT for m > N | — |
| approximate: any factored c with \|c − (a+b)\| ≤ E | E ≥ (a+b)/poly: easy, c = 2^j⌊(a+b)/2^j⌋ with a polynomial leading coefficient (codex). In general no harder than factoring one (n − log₂E)-bit number, by the same formula. E < 1 is exact addition. Not a function of the factors; no reduction *from* factoring known. With E ≈ N^θ it is "construct a factored integer in a short interval" (cf. sihard, Prop SI) | PROVED / OPEN |

**Is there a classification of the "easy" linear functionals?**
- These are the cheaply computable completely multiplicative functions.
- KNOWN for automatic sequences: Konieczny–Lemańczyk–Müllner, *Multiplicative automatic sequences*, Math. Z. 300
  (2022) 1297–1318. Away from finitely many primes these are a Dirichlet character or vanish, extending Klurman–Kurlberg
  (C. R. Math. 357 (2019)) and Schlage-Puchta (Integers 11 (2011) A31). So for finite automata, "addition-free
  homomorphisms" are exactly characters, i.e. floor information.
- For AC⁰, Green (CPC 21 (2012)) proves that μ is orthogonal to AC⁰.
- For polynomial time no classification is known. The natural conjecture is "poly-time completely multiplicative
  functions are χ(n)·n^{it} up to finitely many primes". It would imply that λ and μ are not in P, so it is CONJECTURED
  at most.

---

## 2. Hardcore bits and partial information about the factors (question 2)

### 2.1 What the classical theorems are about (KNOWN)
- **Alexi–Chor–Goldreich–Schnorr**, *RSA and Rabin functions: certain parts are as hard as the whole*, SIAM J. Comput.
  17(2) (1988) 194–209.
  - The LSB of the RSA/Rabin *message* x is hardcore. For Rabin this holds under factoring.
- **Fischlin–Schnorr**, J. Cryptology 13 (2000) 221–244: tighter reductions.
- **Håstad–Näslund**, *The security of all RSA and discrete log bits*, JACM 51(2) (2004) 187–230.
  - Every individual bit of the RSA message (Thm 8.1), and of an injective Rabin variant (Thm 9.1), is hardcore.
  - DL bits: Thm 10.2, excluding the easy 2-primary low bits.
- **Goldreich–Levin**, STOC 1989. For any distribution of canonical pairs on which factoring is average-case hard,
  B(N, r) = ⟨p, r⟩ mod 2 with r uniform and public is hardcore for (p, q, r) ↦ (pq, r).
  - The proof gives noticeable advantage ⇒ list-decode p on a noticeable fraction of N ⇒ verify by division.
  - It needs the varying public mask r.
- **Blum–Blum–Shub / Rabin LSB**: the standard factoring-based hardcore bit. Again it is a bit of a square root, not
  of p.

None of these is a statement about a fixed bit of p as a function of N.

### 2.2 Two single-bit functions of the integer alone that ARE factoring-hard (PROVED; routine)
**Prop GLM** (this note; the GL randomness is embedded in a cofactor; no prior formulation found by us). Define
f(M) = parity(P⁺(M) AND M), the GL inner product of the largest prime factor with the low bits of M itself.
- Setup: N balanced, L = bitlen(q), r ∈ {0,1}^L, k = r·N^{−1} mod 2^L, M = Nk.
- If P⁺(k) < q, which fails for only O(1/L) of r, then P⁺(M) = q and f(M) = ⟨q, r⟩ mod 2.
- Hence a predictor with advantage ε for f on the distribution of M = N·k(r) gives:
  - a GL oracle for q, on a ≥ ε/2 fraction of N, with advantage ≥ ε/2 − O(1/L);
  - Rackoff's version of GL list-decodes q with L·2^l queries, 2^l = O(L/ε²);
  - each candidate is verified by division.
- Toy check: `src/factorbits/gl.py`. Predictor accuracy 0.55 (ε = 0.05) at n = 32 recovers the factor 5/5 times with
  16·2^12 queries. Tests pass.
- *What this shows.* "Single bit ⇒ factoring with amplification" is possible for a function of M alone, but only
  because the reduction controls a free public part of the instance (k) that the bit mixes with the secret.

**Codex's B(M)** (codex_1 §4; worst case, no amplification). B(M) = bit_{v₂(M)}(spf(M/2^{v₂(M)})). Queries
M = 2^iN read p bit by bit.

### 2.3 Status table (question 2; RSA modulus N = pq, p < q, n bits)

| Quantity (from N alone unless stated) | Status | Reduction / reference |
|---|---|---|
| bit i of p, fixed nontrivial i (incl. bit 1 = p mod 4) | **OPEN** | nothing found (codex search UNVERIFIED as a negative) |
| bit i of p with i an input | ≡ FACT | query all i |
| ⌊p/2^k⌋, k ≤ n/4 + O(log n) | ≡ FACT | Coppersmith J. Cryptology 10 (1997) Thms 4–5 (≈ n/4 top bits) plus enumeration of O(log n) bits; Rivest–Shamir 1985 needed n/3 |
| ⌊p/2^k⌋, n/4 + ω(log n) ≤ k < n/2 − O(1) | **OPEN** (sufficient ⇐ ≈ n/4 bits; no converse) | topblock: lattice comparison theorem |
| top O(1) bits of p | **OPEN** for exact computation; trivially predictable with 1/poly success by guessing; partly fixed by p ≤ √N | — |
| p mod M, log₂M ≥ n/4 + O(log n) | ≡ FACT | Coppersmith/May with known modulus (trackE E1) |
| p mod ℓ (ℓ small) | **OPEN**; = s mod ℓ (symmetric) + a 1-bit label | trackE F2 |
| s mod ℓ = (N + 1 − φ(N)) mod ℓ | **OPEN**; decides "ℓ* is a square mod N" when (N/ℓ) = 1 | nonsmooth Prop RANK |
| φ(N) mod 2^k, k ≥ n/4 + v₂(q−p) | ≡ FACT (2 decodes) | trackE E2 |
| φ(N) mod 8 (N ≡ 1 mod 4) | ⇔ ε | trackE E2′ |
| ε = (−1/p), N ≡ 1 (mod 4) | **OPEN**; = "−1 is a square mod N" = [r₂(N) = 16] | fixed-base QR; Seurin's 2-Φ/4-hiding (PKC 2014; eprint 2013/256) |
| bit 1 of the smaller factor, N ≡ 3 (mod 4) | **OPEN** (an orientation bit: which of p, q is ≡ 3 mod 4) | — |
| [ℓ \| φ(N)] (Φ-hiding) | **OPEN**; conjectured hard on suitable distributions | Cachin–Micali–Stadler, EUROCRYPT 1999 |
| QR mod N, random base with Jacobi +1 | **OPEN** (the classical QR-vs-factoring problem) | — |
| Ω(N) mod 2, ω(N), d(N), μ(N), squarefreeness, rad(N) | easy on the RSA promise (0, 2, 4, 1, yes, N); **OPEN** unrestricted | listed open in Adleman, FOCS 1994 (C7–C9: squarefree part, squarefreeness, ω), the companion of Adleman–McCurley ANTS-I 1994 (exact ANTS numbering UNVERIFIED); Shallit–Shamir IPL 20 (1985); Landau, Inf. Comput. 78 (1988); Woll, Inf. Comput. 72 (1987) (details UNVERIFIED) |
| d(N) mod 2 | easy (square test) | — |
| φ(N), λ(N), σ(N), ord_N(a) | ≡ FACT (randomized) | Miller 1976; Bach–Miller–Shallit SIAM J. Comput. 15 (1986) |
| "N has a factor ≤ x", x an input | ≡ FACT (binary search, ≤ n calls) | — |
| smallest / largest prime factor | ≡ FACT | strip and recurse |
| #E(Z/N), chosen curves | ≡ FACT (randomized; ≤ 6 calls expected via twists, codex) | Kunihiro–Koyama EUROCRYPT 1998 |
| #E(Z/N), one given curve | OPEN | — |
| r₂(N) | = 16·[ε = 1] on N ≡ 1 (4): **OPEN**; 0 on N ≡ 3 (4) | — |
| one sum-of-two-squares representation | OPEN; two inequivalent ones factor N | — |
| ⟨p, r⟩, r public uniform | **hardcore** (average case) | Goldreich–Levin |
| f(M) = parity(P⁺(M) ∧ M) on M = N·k | **hardcore**, amplifies | Prop GLM |
| message bits of RSA/Rabin | hardcore | ACGS; Håstad–Näslund |

**Answer to "is there any 'one fixed bit of p is as hard as factoring' theorem?"**
- No. Codex's search found none; that negative is UNVERIFIED.
- All positive single-bit results take an extra input (r, x, y, i) or embed one in the instance (GLM, B(M)).
- For the natural fixed bits, §3 and §4 explain where a proof would have to come from.

---

## 3. The middle ground: representations that leak (question 3)

Everything computed from N is a function of N. So "leak" is a computational notion: does a poly-time R(N) predict a
function of p better than the *floor*, the best predictor from N mod m? (trackE F1.)

### 3.1 Classes with a proved floor
- **(A) R factors through N mod m.** PROVED (trackE F1; codex). This holds for any distribution of (p, q), and
  conditioned on C = N mod m, I((p mod m, q mod m); R | C) = 0. With m ≤ (log X)^A, Siegel–Walfisz makes the labelled
  pair asymptotically uniform on the fiber F_c = {(u, c/u)}, so the floor itself is ≈ 0 bits about ε, bits of p, etc.
  Caveat (codex): many small moduli together are one large modulus, so the bound must apply to lcm(m_i).
- **(B) Abelian reciprocity aggregates.** Characters and power-residue symbols χ(N) of fixed conductor, and Artin
  symbols of (N) in abelian fields, are functions of N mod the conductor. This reduces to (A) (nonsmooth §5.4,
  Kronecker–Weber). It covers products ∏χ(r)^{e_r}, not component selections like χ(p) + χ(q).
- **(C) Multiplicative automatic sequences.** These are Dirichlet characters up to finitely many primes (KLM 2022), so
  they reduce to (A).
- **(D) Pure generic-ring decision trees. Theorem PT, below.**

### 3.2 Theorem PT (pure trees: information costs splitting)

**Model.**
- A *pure* tree T is a decision tree, possibly adaptive and of any size, possibly randomized with coins independent
  of N.
- Node v carries an integer c_v given by a straight-line program over Z (+, −, ×, integer constants). The program and
  constants may depend on the path so far but not otherwise on N.
- On input N the node computes c_v mod N. It reports a split if 1 < gcd(c_v mod N, N) < N, and otherwise branches on
  [c_v ≡ 0 (mod N)].
- The output is a function of the leaf.
- Division a/b is two nodes: invertibility of b, and then the test on a. Equality tests are tests of differences.
- *Included:* Pollard rho, p−1, p+1, ECM with N-independent random curves, Pollard–Strassen/BGS products, Lipton's SLP
  factoring.
- *Excluded:* anything whose program depends on more of N (Fermat/Miller–Rabin exponents N−1, √N, x² − N, forms of
  discriminant −4N, CF of √N, lattices containing N), and anything that reads the integer representative of a residue.
- A side-information variant (PT-m) lets the root read c = N mod m and pick a pure tree T_c.

**Partition lemma (PROVED).** For a prime x, the *diagonal path* of x answers [x | c_v] at each node; let its leaf be
class(x). For distinct primes x, y, the instrumented run of T on N = xy splits N iff class(x) ≠ class(y). Otherwise it
ends in the leaf class(x) = class(y).

*Proof.* gcd(c_v, xy) ∈ {1, x, y, xy}, and it is x or y iff [x | c_v] ≠ [y | c_v]. Otherwise the answer
[xy | c_v] equals [x | c_v] = [y | c_v]. By induction the run follows both diagonal paths while they agree, and at the
first node where they disagree it splits. ∎

So the non-split behaviour of every pure tree is a partition of the primes into classes S_ℓ.

**Theorem PT (PROVED).** Let X, Y be iid with law μ on primes. Let F(X, Y) be any target, e.g. bit i of min(X, Y),
or X mod 4. Let π_ℓ = μ(S_ℓ), η = 1 − max_ℓ π_ℓ and P_split = 1 − Σ_ℓ π_ℓ² (on X = Y, define the run to end in
class(X)). Then P_split ≥ η, and

  Pr[T outputs F(X, Y) or splits N] ≤ max_b Pr[F = b] + P_split + η² ≤ max_b Pr[F = b] + 2·P_split.

*Proof.*
- Σπ² ≤ max π, so P_split ≥ η.
- Let ℓ* be the heaviest class. Pr[no split and class ≠ ℓ*] = Σ_{ℓ≠ℓ*} π_ℓ² ≤ (max_{ℓ≠ℓ*} π_ℓ)·η ≤ η².
- On {no split, class ℓ*} the output is the fixed value o(ℓ*), correct with probability ≤ max_b Pr[F = b].
- Sum the three cases. ∎

**Corollary PT1.** Advantage δ over the trivial predictor (counting splits as successes) forces P_split ≥ δ/2. So
within this model *every* predicate of the factorisation is hardcore with loss 2. Weak predictors are already
factoring algorithms.

**Corollary PT2 (ε; the fiber leak is paid on the cross population). PROVED.**
- Population A (N ≡ 1 mod 4): ε = ±1 with probability 1/2, with X, Y iid μ₁ (primes ≡ 1 mod 4) or iid μ₃.
- Population B (N ≡ 3 mod 4): X ~ μ₁, Y ~ μ₃.
- Let α_ℓ = μ₁(S_ℓ) and β_ℓ = μ₃(S_ℓ) for the tree T₁ that PT-m uses on N ≡ 1 (mod 4).

*Proof.*
- On non-split runs the output is g(ℓ), so Pr_A[correct or split] ≤ ½Σ_ℓ max(α_ℓ², β_ℓ²) + P_split(A).
- ½Σ max(α², β²) = ½(1 − P_split(A)) + ¼Σ|α² − β²|.
- Σ|α² − β²| ≤ ‖α − β‖₂‖α + β‖₂ ≤ 2‖α − β‖₂.
- ‖α − β‖₂² = 2P_split(B) − 2P_split(A), where P_split(B) = 1 − Σα_ℓβ_ℓ is the split rate of *the same pure tree*
  T₁ run on the mixed N ≡ 3 (mod 4). Hence

  Pr_A[correct or split] ≤ ½ + ½P_split(A) + ½√(2P_split(B)).

- Advantage δ therefore forces P_split(A) ≥ δ or P_split(B) ≥ δ²/2. ∎

The general form for any two fiber points u, u′: ‖π^{(u)} − π^{(u′)}‖₂² ≤ 2·s(μ_u × μ_{u′}). So PT-m leaks fiber
information only as fast as T_c splits the cross population, whose N lie in residue class u·u′ mod m. Example: the
pure tree c = ∏_{x ≤ (p−1)/2}(x² + 1) decides ε exactly and splits *every* mixed N. Its cost is exponential, which PT
ignores.

**Low-height strengthening (codex_1 §3, PROVED there).**
- Suppose there are t tests with integer heights ≤ H, and primes lie in [X, 2X].
- Then the transcript equals the formal integer-zero transcript except with probability O(tH/X), or
  O(tHφ(m)/X) with side information. For polylog t, H there is no leak at all.
- PT is the complement for unbounded height, where repeated squaring makes H exponential: there, information equals
  splitting.

**Relation to the literature.**
- The core observation, "a component-separating zero is a gcd factor", is the standard engine of generic-ring
  reductions: Aggarwal–Maurer EUROCRYPT 2009; Jager–Schwenk J. Cryptology 26 (2013).
- The quantitative predicate version for N-oblivious programs (PT, PT2) is no prior formulation found by us.
- Caveat (Jager–Schwenk): generic-ring hardness is not ordinary hardness. Computing the Jacobi symbol of a given
  element is generically as hard as factoring yet easy by reciprocity. PT says nothing about non-pure algorithms, and
  that is where ε could be decided cheaply if it can be at all.

**E1 (pre-registered; all predictions met).**
- Prime sets: all b-bit primes with the top two bits set, b = 14, 16, 18 (n = 28, 32, 36).
- Trees: nested p−1; adaptive p−1/p+1/ECM; three ECM curves; and at b ≤ 16 a Lipton interval product and full and
  partial ∏(x² + 1).
- Partition lemma against real instrumented execution mod N: 0 mismatches in 681 990 runs (all 92 665 pairs × 6 trees
  at n = 28; 20 000 or 2 000 sampled pairs per tree at n = 32, 36).
- PT violations 0/75 and ε-cross violations 0/15.
- Descriptive: for the natural trees the non-split information is small.
  - At n = 36 the largest non-split excess over the "no information" rate P_split + (1 − P_split)·baseline is 0.003
    (adaptive), 0.006 (ECM) and 0.025 (nested p−1), against P_split = 0.42, 0.73 and 0.43.
  - The p−1 excess comes from the smooth-order class correlating with p mod 4 (p ≡ 1 (4) makes p−1 smoother). It is
    within η² = 0.075.
  - `qr_sniff_full` decides ε exactly with P_split(A) = 0, P_split(B) = 1, as designed.
  - `qr_sniff_partial` (n = 32): ε advantage 0.036 with P_split(A) = 0.035, P_split(B) = 0.036.

### 3.3 What is open
- **(E) Finite automata and AC⁰ on the bits of N**, for factor bits of RSA moduli. OPEN.
  - Green's AC⁰ theorem (μ), Bernasconi–Shparlinski (STACS 1999; squarefreeness) and Allender–Saks–Shparlinski
    (JCSS 62 (2001); primality ∉ AC⁰[r]) concern other targets.
  - By Linial–Mansour–Nisan, an AC⁰ floor for ε would follow from bounds on correlations between ε(p) and parities of
    polylog many digits of pq. That is a type-II bilinear digit sum and looks within reach of Mauduit–Rivat/Bourgain
    methods. We did not attempt it.
  - No factor-bit lower bound was located in Shparlinski's 2003 book (UNVERIFIED negative).
- **(F) Non-pure algebra**: Fermat/Miller–Rabin/Lucas tests with exponent N ± 1, Jacobi sums, integer
  representatives of residues. OPEN as a class.
  - The eps battery found no predictor among 66 such statistics at |ρ| ≥ 0.022 (notes/eps_results.md).
- **(G) Euclid-type and lattice representations of N**: CF of √N, forms of discriminant −4N, Coppersmith lattices,
  sieving. OPEN.
  - Leakage observed here so far comes only through relations, i.e. smooth numbers. Examples: the eps-battery n = 32
    CF hits, and nonsmooth's (U3) sums, which are divisor detection.

So the sharpest statement is a trichotomy:
- **floor proved** for local, abelian, automatic and low-height pure classes;
- **"leak ⇒ split"** (Theorem PT) for all pure algebraic classes of any height;
- **open** for everything that uses N's integer structure beyond N mod m. That is exactly where every algorithm faster
  than N^{1/4} lives (Harvey/Hittmeir's N^{1/5}, CFRAC, QS, NFS).

---

## 4. Leverage and self-reduction (question 4)

### 4.1 Minimal leakage that is known to factor
- **Bits of p in a coset shape** (top bits, low bits, p mod M): n/4 + O(log n) bits; trackE measured n/4 + n/(4d) − 2
  in dimension d. n/4 − t bits cost a 2^t enumeration, so ω(log n) missing bits are not polynomially recoverable by
  this route.
- **Random positions of p and q**: density > 1/2 in the low n/4 positions, about n/4 bits total (trackE E6;
  Heninger–Shacham CRYPTO 2009, heuristic).
- **Arbitrary yes/no queries**: εn queries for every fixed ε > 0, in polynomial time, CONDITIONAL on an ECM-type
  conjecture (Maurer, Computational Complexity 5 (1995) 237–247; ETH lists 1996). Rivest–Shamir (EUROCRYPT '85) need
  n/3.
- **Noisy bits** with independent errors at every position of (p, q): error < 0.084 (Henecka–May–Meurer CRYPTO 2010,
  heuristic). The capacity bound is 0.110, i.e. 1 − H(δ) ≥ 1/2 (Paterson–Polychroniadou–Sibborn ASIACRYPT 2012).
  - A per-bit advantage ε gives ≈ (2/ln 2)ε² bits per position, so ≈ n/4 bits of information need ε ≳ 0.3 per bit.
    Inverse-polynomial advantages are far below every known decoder.
- **Lower side**: O(log n) bits of anything are worthless (enumerate). A general minimal threshold is OPEN.

### 4.2 Why ACGS-type amplification has no analogue for fixed bits of p
- **Where amplification comes from.**
  - RSA: c ↦ a^e c gives instances with secret ax.
  - Rabin: multiply by a².
  - DL: (g^x)^a g^b gives secret ax + b.
  - So a predictor for one fixed bit position, queried on many transformed instances, yields many *different*
    functionals of one secret. That is the input ACGS, Håstad–Näslund and Blum–Micali decode (codex_1 §4 agrees).
- **What is missing for factoring.**
  - The secret is the support {p, q} of the exponent vector.
  - The instance-computable maps are N ↦ N·k (adds a known vector) and N ↦ N^j (scales). Both leave the unknown
    support invariant, so a predictor for bit i of "the" factor returns the same bit every time.
  - Amplification would need a computable action that *moves* p.
- **Lemma PSR (PROVED).** Let g, h ∈ Z[x], g nonconstant, h ≠ 0. Suppose that for infinitely many primes q there are
  infinitely many primes p with g(p) | h(pq). Then g = c·x^j.
  - *Proof.* Fix such a q. Write D·h(qx) = A(x)g(x) + R(x) with A, R ∈ Z[x] and deg R < deg g. Then g(p) | R(p).
    Since |R(p)| < |g(p)| for large p, R(p) = 0 infinitely often, so R = 0 and g(x) | h(qx) in Q[x].
  - A root α ≠ 0 of g would make h(qα) = 0 for infinitely many q, contradicting h ≠ 0. ∎
  - Codex proved the "all large pairs" version and also obtained x^j | h.
  - So polynomial instance maps can only expose p itself (or p^j). They cannot expose p + 1, 2p + 1, or any re-masked
    p. Randomized polynomial maps are covered coefficient-wise. General poly-time maps are NOT covered.
  - Divisibility sequences (a^N − 1 ∋ a^p − 1, Lucas u_N ∋ u_p, elliptic divisibility sequences) do carry p to a
    different divisor, but at exponential size. Reducing them modulo anything loses the divisibility.
- **What does work: varying a public mask.**
  - Goldreich–Levin (⟨p, r⟩) and Prop GLM (r embedded as N·k) amplify from advantage 1/poly, because the
    reduction varies a public string that the bit mixes with p.
  - In the pure generic model (Theorem PT) there is nothing to amplify: any advantage already is a splitting
    probability.
- **Answer to question 4.**
  - A weak per-bit predictor for a FIXED natural bit position of the factor is not known to amplify. Not even a
    perfect one is known to give anything beyond that single bit per N.
  - The obstruction is the lack of a factor-moving self-reduction. We proved there is none among polynomial maps,
    and that pure algebraic predictors are already splitters.
  - The question "computing bit i of p for random N-related inputs ⇒ factoring" reduces exactly to finding a
    poly-time map family N ↦ N′ whose factors are known, varied functions of p (not p itself).
  - Status: OPEN. We know no candidate. Lemma PSR excludes the polynomial ones, and divisibility sequences are too
    large.

---

## 5. Candidate cards (rule 9)

| Mechanism | Model | Outcome |
|---|---|---|
| ADD/successor/partial-output oracles | oracle reductions | all equivalent to FACT except the OPEN rows of §1.5; no mechanism for *computing* ADD faster |
| pure trees (p−1, p+1, ECM, rho, products) as ε / factor-bit predictors | Theorem PT | **refuted for this model**: information ⇒ split, quantitatively |
| polynomial self-reduction for bits of p | Lemma PSR | **refuted for this model** |
| GL / embedded-cofactor bits | Prop GLM | works, but needs a free public component; not a natural bit |
| automaton / AC⁰ readers of N | — | not yet tested (OPEN; LMN route sketched) |

## 6. Most decisive open question
**Is there a polynomial-time map family N ↦ N′ (any algorithm, not polynomial in N) whose prime factors include known
and varied functions of p other than p itself?** A positive answer would give an ACGS-style amplification for fixed
factor bits. Divisibility sequences are the only natural examples, and they are exponentially large. The negative
side is equally open. In the pure algebraic world the question is moot (Theorem PT). Outside it, the first concrete
target is unchanged from trackE §7: a non-pure poly-time observable with 1/poly advantage on ε.

## 7. Reproduction
```
.venv/bin/python -m pytest -q tests/test_factorbits.py                 # 48 tests
.venv/bin/python -m src.factorbits.run_experiments e1 14 16 18          # E1 -> experiments/raw/factorbits_e1.jsonl (~3 min)
```
