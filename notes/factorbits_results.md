# factorbits: prime-exponent addition, hardcore bits of factors, and the information floor (phase 21)

Agent `factorbits`, 2026-10-21.
- Pre-registration: `notes/factorbits_prereg.md` (commit `3c9fe3f`, before E1 ran).
- Code: `src/factorbits/{equiv,harness,pure,gl,run_experiments}.py`. Tests: `tests/test_factorbits.py` (51 pass).
- Raw data: `experiments/raw/factorbits_e1.jsonl`, plus `factorbits_e1_rerun.jsonl` (§3.2, deviation D1).
- GPT-6.1-sol: blind round `codex/factorbits_1.{md,out}`; adversarial round `codex/factorbits_2.{md,out}`
  (MAJOR REVISION; applied, §8).
- Two read-only Opus critics (§9).
- Labels:
  - **KNOWN** (with reference);
  - **PROVED** (proof here; routine unless stated; "no exact prior formulation located" is the strongest novelty claim);
  - **HEURISTIC**;
  - **CONDITIONAL**;
  - **OPEN**;
  - **experimentally supported**;
  - **UNVERIFIED**.
- "OPEN" means that no polynomial algorithm and no reduction *from* factoring is known to us. Factoring always
  computes these quantities. Negative literature conclusions come from our and codex's searches and are UNVERIFIED as
  negatives.

Notation:
- n = bit length of N.
- N = pq is a balanced semiprime with p < q and bitlen(p) = bitlen(q) = n/2 (the convention of `balanced_semiprime`).
- A *vector* is the prime-exponent representation {prime: exponent} of a positive integer.
- FACT is complete factorisation of arbitrary integers.

---

## 0. Verdict (plain language)

**1. The equivalence is real and cheap.**
- Fast addition in the prime base is the same problem as factoring arbitrary integers.
  - A successor oracle x ↦ x + 1 factors N in popcount(N) − 1 ≤ log₂N calls (Horner; doubling is free).
  - A full addition oracle needs O(log log N) calls, randomized and unconditional (write X = u + v with v a random
    prime in a short interval).
  - Heuristically it needs one call.
  - Conversely, one factoring call adds.
- The Zech-logarithm analogue is exact: a discrete log in F_p costs ≤ ⌊log₂p⌋ Zech calls, and a Zech value costs one
  discrete-log call.
- The structural difference: primes are free, dense and recognisable in the prime base. F_q has no such set of
  elements whose logs are free.
- Partial forms:
  - As hard as factoring: the support of a+b, its exponents mod any fixed m (randomized polynomial reduction, from
    codex), the largest prime factor (one call on semiprimes heuristically; quasi-polynomial in general), and sums of
    many smooth terms.
  - Easy: the smooth part; a smooth-promised result; the residue of a+b; functionals ⟨v, w⟩ that are logs, or
    character indices *taken mod the character order*; approximation to relative error 1/poly.
  - Open: ω, Ω mod 2, μ, squarefreeness, the numerical radical, and two smooth summands.
  - Everything needs a size convention. With exponents written in binary, "addition" would factor Fermat numbers
    from O(k)-bit inputs.

**2. No natural single bit of the factors is known to be as hard as factoring.**
- Every known factoring-hardcore bit is a bit of something else, or uses auxiliary randomness:
  - the RSA/Rabin message bits of ACGS and Håstad–Näslund;
  - the Goldreich–Levin parity ⟨p, r⟩ with a public random r.
- Two single-bit functions of the integer alone are factoring-hard, but both are artificial:
  - Prop GLM (here): f(M) = parity(P⁺(M) AND M). It is average case, and any 1/poly advantage amplifies to factoring
    with no loss.
  - codex's B(M): worst case.
  - Both work because the reduction can hide fresh randomness, or a bit index, in a cofactor.
- For every natural fixed quantity the status is OPEN: bit i of p, p mod 4, ε = (−1/p), bit 1 of the smaller factor,
  ω, μ, squarefreeness, φ(N) mod ℓ, fixed-base QR.
- Things that are as hard as factoring:
  - ⌊p/2^k⌋ for k ≤ n/4 + O(log n);
  - p mod M for log₂M ≥ n/4 + O(log n);
  - φ(N) mod 2^k for k ≥ n/4 + v₂(q−p);
  - "N has a factor ≤ x" with x as input.
  - All of these go through Coppersmith or binary search.

**3. The information floor is a theorem for N-oblivious algebraic algorithms.**
- Theorem PT (PROVED; elementary; no exact prior formulation located) covers pure generic-ring decision trees: finite
  trees of integer straight-line programs whose constants do not depend on N, with zero tests mod N, of any size and
  height, adaptive and randomized. With fixed schedules this covers Pollard rho, p−1, p+1, ECM and Pollard–Strassen.
  - Such a tree either splits N or ends in a leaf determined by a *partition class* of p, which is also the class
    of q.
  - Consequently, for iid factors, an advantage δ in predicting *anything* about the factors forces the same tree,
    with gcds, to split N with probability ≥ δ/2.
- With side information N mod m (PT-m), the tree learns, without splitting, at most what the local fiber
  {p mod m, q mod m} gives, up to 2·P_split.
- For ε the fiber information has a price: the tree T₁ that a pure algorithm uses on N ≡ 1 (mod 4) decides ε with
  advantage δ only if T₁ splits N ≡ 1 (mod 4) with probability ≥ δ, or splits N ≡ 3 (mod 4) with probability
  ≥ δ²/2 when run on them.
- PT is **not a time barrier**: ECM is pure and heuristically sub-exponential. It says pure methods are
  all-or-nothing.
- Floor statements also hold for:
  - functions of N mod m (trackE F1);
  - abelian reciprocity symbols (nonsmooth §5.4);
  - multiplicative automatic sequences (Konieczny–Lemańczyk–Müllner).
- OPEN: finite automata and AC⁰ on N's bits, and every non-pure, N-dependent class (Fermat/MR exponents, CF of √N,
  forms of discriminant −4N, sieving, lattices).
- E1: the partition lemma matched 681 990 real instrumented runs mod N with 0 mismatches. The two inequalities held
  in all 90 exact numerical checks.

**4. Leverage: a weak predictor for a fixed bit of p is not known to amplify, and the reason is structural.**
- ACGS and Blum–Micali amplify because the secret lives in a group that acts on instances: x ↦ ax for RSA, x ↦ ax + b
  for DL. For factoring the secret is the *support* of the exponent vector, and the computable instance maps (N·k,
  N^j) only add known vectors. A predictor for bit i of "the" factor therefore returns the same bit every time.
- Lemma PSR (PROVED, also in a balanced-modulus version): no fixed polynomial map N ↦ h(N) carries p to any g(p) other
  than c·p^j.
- Divisibility sequences (2^N − 1 ∋ 2^p − 1) do move p, but at exponential size.
- Known sufficient leakage:
  - ≈ n/4 + O(log n) bits of p in coset shape (Coppersmith; trackE measured n/4 + n/(4d) − 2);
  - εn arbitrary yes/no queries (Maurer 1995, conditional).
- O(log n) bits are worthless.
- A single fixed bit would need a factor-moving self-reduction. We found none, and we rule out the polynomial ones.

---

## 1. The equivalence (question 1)

### 1.1 Size convention (PROVED; codex agrees)
The vector of a has encoding size Σ(log p_i + log e_i), but a has Σ e_i log p_i bits. For example (2: 2^k) has O(k)
bits and denotes a (2^k + 1)-bit integer. On such succinct inputs an addition oracle would factor 2^{2^k} + 1 from an
O(k)-bit query. So **ADD_n** takes inputs with log₂a, log₂b ≤ n, and complexity is measured in n. Equivalently, bound
the weighted exponent sum.

### 1.2 Prop E1 (ADD ⇒ FACT)

| Oracle | Calls to factor any N completely | Status |
|---|---|---|
| successor S(x) = vector of x+1 | exactly popcount(N) − 1 ≤ ⌊log₂N⌋ for this algorithm. Horner: doubling adds 1 to the exponent of 2 for free. Queries only even x < N | PROVED (`factor_via_successor`) |
| ADD, base 2^w windows | number of nonzero base-2^w digits after the leading one ≤ ⌈n/w⌉ − 1. With w = Θ(log n) this is O(n/log n); digits are trial-divided | PROVED |
| ADD, short-interval prime trick | O(log log N) calls. Write X = u + v with v uniform in [X − X^θ, X), fixed 7/12 < θ < 1, until v is prime: expected O(log X) certified (AKS) primality tests, by Huxley's asymptotic π(X) − π(X − X^θ) ~ X^θ/log X (Invent. Math. 15 (1972)). Factor u ≤ X^θ recursively down to poly(n) size: depth log_{1/θ}(n/log n) | PROVED, Las Vegas (`factor_via_addition_shortint`; ≤ 6 calls at 96–128 bits in tests) |
| ADD, one call | v = the largest prime < N with u = N − v trial-factored. Prime gaps are O(log²N) under Cramér's model | HEURISTIC; 10/10 at 64 bits |
| S, one call | t ≤ poly(n) with tN − 1 prime; S on it; read p, q among the primes of tN | HEURISTIC; 10/10 at 64 bits |
| S, one call under GRH | GRH gives π(x; N, −1) = Li(x)/φ(N) + O(√x log(xN)). With t uniform ≤ T = N(log N)^C, C > 4, tN − 1 is prime with probability (1 + o(1))·N/(φ(N) log(NT)). Queries have ≈ 2n bits. t need not be factored: p and q are read off the returned factorisation of tN by testing which primes divide N | CONDITIONAL (GRH) (codex) |
| many-summand ADD of 2-smooth terms | 1 call: N = Σ_{b_i=1} 2^i | PROVED |

### 1.3 Prop E2 (FACT ⇒ ADD). PROVED
- With a complete-factorisation oracle, addition costs 1 call.
- With a splitting oracle and primality tests, recursive splitting costs exactly Ω(a+b) − 1 ≤ log₂(a+b) − 1 calls.
- So ADD_n ≡ FACT under polynomial-time reductions.
- FACT ⇒ ADD needs factoring of arbitrary integers a + b, so this is not an equivalence with factoring *semiprimes*.
  The converse direction does not help either: the ladder queries arbitrary integers. An ADD promised only on
  semiprime outputs gives one-call factoring of N only together with a factored decomposition N = u + v. That is
  available heuristically by the prime trick.

### 1.4 Prop Z (Zech ⇔ DL in F_q). PROVED (classical in substance)
Z(k) = log_g(1 + g^k), with a separate zero symbol. Addition in log coordinates is log(g^u + g^v) = u + Z(v − u).
- **DL ⇒ Zech:** 1 DL call.
- **Zech ⇒ DL in F_p:**
  - log 2 = Z(0) costs 1 call.
  - Horner on the bits of h then costs popcount(h) − 1 calls: doubling adds log 2, and "+1" is one Zech call.
  - The total is popcount(h) for h > 1. Since 1 ≤ h ≤ p − 1 forces popcount(h) ≤ ⌊log₂p⌋, the total is
    ≤ ⌊log₂p⌋.
  - No partial value vanishes, because all prefixes are integers in [1, p).
- **Windowing:** 2^w − 2 calls of digit-log precomputation (these are charged, unlike integer digit factorisations)
  plus one call per nonzero digit, O(log p/log log p) in total.
- **F_{p^k} in the basis 1, g, …, g^{k−1}** (codex): δ + Σ_{c_j≠0}(wt₂(c_j) − 1) + (s − 1) = Σ wt₂(c_j) − 1 + δ
  ≤ k⌊log₂p⌋ calls. Here s is the number of nonzero coordinates, δ = [some c_j > 1], and the bound uses
  wt₂(c) ≤ ⌊log₂p⌋. Partial sums of independent basis terms never vanish.

| integers | F_q^* |
|---|---|
| free abelian group ⊕_p Z (vector) | Z/(q−1) (log_g) |
| FACT | DL |
| ADD | Zech |
| primes: free, dense, recognisable vectors | no analogue: an element's log is free only if we built it as g^k |
| N mod m, characters (positional floor) | Pohlig–Hellman projections, power-residue symbols |
| secret = the support of the vector | secret = an element of a group |
| no random self-reduction known | x ↦ ax + b random self-reduction (Blum–Micali; Håstad–Näslund) |

The primes row explains why integer ADD needs O(log log N) calls rigorously and one heuristically, while Zech uses
Θ(log p/log log p) by the known method. We know no lower bound for either. The last row is §4.

### 1.5 Partial forms
In this subsection the oracle takes full vectors of a and b and returns partial output about a + b. Status is PROVED
unless marked.

| Partial ADD | Status | Why |
|---|---|---|
| support of a+b (list of primes) | ≡ FACT, same call counts as full ADD | exponents by repeated division (`factor_via_support`) |
| exponents mod m, m > log₂(a+b) | ≡ FACT | nothing wraps |
| exponents mod a fixed m ≥ 2, zero residues dropped | ≡ FACT under randomized polynomial reductions (codex). Query (u, v) as in the short-interval trick and strip the returned primes. The rest is r^m with r ≤ X^{1/m}; recurse on r. C(b) ≤ C(θb) + C(b/m) + 1 = O(b^{α}) with θ^α + m^{−α} = 1, α < 2 for θ = 3/5 | `factor_via_expmod_shortint`, tested on hidden square/cube parts. The deterministic prefix ladder gives only n^{⌈log_m n⌉} |
| largest prime factor P⁺(a+b) only | semiprime N: 1 call HEURISTIC (prime trick; it returns q). In general: randomized quasi-polynomial b^{O(log b)} (codex: strip P⁺ repeatedly, preparing each query by the short-interval recursion). A polynomial reduction is OPEN | — |
| numerical radical rad(a+b) | OPEN | squarefree a+b gives no information |
| ω, Ω mod 2 (Liouville), μ, squarefreeness, d of a+b | OPEN | constant on semiprime outputs, but an unrestricted oracle could still be queried on other sums. Shallit–Shamir (IPL 20 (1985)) relate these to d(·) without recovering primes |
| d(a+b) mod 2 | easy | it is 1 iff a+b is a square |
| B-smooth part, B ≤ poly(n); result promised poly-smooth | easy | trial division |
| B-smooth part, larger B | HEURISTIC L_B[1/2, √2] (ECM) | — |
| ⟨v, w⟩ ∈ Z for w = log p | easy | log(a+b) |
| ⟨v, w⟩ ∈ Z for w finitely supported | easy | repeated division |
| ⟨v, w⟩ **mod ord χ** for w_p = ind χ(p) | easy | it equals χ(a+b) |
| the same sum over **Z** | OPEN | for χ₄ with w_p ∈ {0, 1}, on N ≡ 1 (4) it equals 2·[ε = −1] |
| ⟨v, 1⟩ mod 2 (Liouville λ) | OPEN | — |
| ⟨v, 1_S⟩ mod 2 for a set S of primes | OPEN in general. On semiprimes, with S = {primes whose bit i is set}, it is bit_i(p) ⊕ bit_i(q) | — |
| two y-smooth summands, y ≤ poly(n) | OPEN | It contains Cunningham factoring (2^k ± 1). To use it on N we need kN = a + b with a, b smooth, i.e. a multiplicative relation mod N. About π(y) such relations factor N *without* the oracle by Dixon combination (HEURISTIC; fruitfulness can fail for structured producers, sicheck Prop B) |
| residue (a+b) mod m | easy | — |
| vector of the least residue | easy for m ≤ poly; ≡ FACT for m > N | — |
| approximate: some factored c with \|c − (a+b)\| ≤ E, 1 ≤ E < a+b | Easy for E ≥ (a+b)/poly: c = 2^j⌊(a+b)/2^j⌋ with j = ⌊log₂E⌋ (codex). In general no harder than factoring one (n − log₂E + O(1))-bit number by the same formula, an upper reduction only. E < 1 is exact addition. It is not a function of the factors and no reduction from factoring is known. With E ≈ N^θ it is "construct a factored integer in a short interval" (cf. sihard, Prop SI) | — |

**Classifying the "easy" homomorphisms.**
- The linear functionals with values in a group A are exactly the completely multiplicative maps to A.
- KNOWN for automata: multiplicative automatic sequences are, away from finitely many primes, a Dirichlet character
  or vanish (Konieczny–Lemańczyk–Müllner, Math. Z. 300 (2022) 1297–1318). This extends Klurman–Kurlberg (C. R. Math.
  357 (2019) 752–755) and Schlage-Puchta (Integers 11 (2011) A31). For finite automata the addition-free homomorphisms
  are exactly characters, i.e. floor information.
- Green (CPC 21 (2012) 942–951) shows that μ is asymptotically orthogonal to AC⁰.
- For polynomial time no classification is known. The natural conjecture, "poly-time computable completely
  multiplicative functions are χ(n)n^{it} away from finitely many primes", would imply that λ is not in P. It is
  CONJECTURED at most.

---

## 2. Hardcore bits and partial information about the factors (question 2)

### 2.1 What the classical theorems are about (KNOWN)
- **Alexi–Chor–Goldreich–Schnorr**, SIAM J. Comput. 17(2) (1988) 194–209. The LSB of the RSA/Rabin *message* is
  hardcore; for Rabin this rests on factoring.
- **Fischlin–Schnorr**, J. Cryptology 13 (2000) 221–244: tighter reductions.
- **Håstad–Näslund**, JACM 51(2) (2004) 187–230. Every individual bit of the RSA message is hardcore (Thm 8.1), as is
  every bit of an injective Rabin variant (Thm 9.1). DL bits: Thm 10.2.
- **Goldreich–Levin**, STOC 1989 25–32. On any distribution of canonical pairs where factoring is average-case hard,
  ⟨p, r⟩ mod 2 with r public and uniform is hardcore for (p, q, r) ↦ (pq, r). It needs the varying public mask r.
- **BBS / Rabin LSB**: a bit of a square root, not of p.

None of these concerns a fixed bit of p as a function of N.

### 2.2 Single-bit functions of the integer alone that are factoring-hard (PROVED; routine)

**Prop GLM** (here; no exact prior formulation located, and the decoding is standard GL). The function is
f(M) = parity(P⁺(M) AND M).
- Setup: L = bitlen(q), so q ≥ 2^{L−1}. Take r ∈ {0,1}^{L−1}, k = r·N^{−1} mod 2^{L−1}, M = Nk, with k = 0 handled
  arbitrarily.
- Then k < 2^{L−1} ≤ q, so **P⁺(M) = q always**, and M ≡ r (mod 2^{L−1}). Hence
  f(M) = ⟨q mod 2^{L−1}, r⟩ ⊕ bit_{L−1}(M).
- Take a predictor with advantage ε for f on the distribution of M = N·k(r) (N from the balanced distribution,
  r uniform). It yields the GL oracle B(r) = pred(M) ⊕ bit_{L−1}(M) for x = q mod 2^{L−1}, with exactly the same
  advantage, since the embedding is exact.
- So on ≥ ε/2 of the N, B has advantage ≥ ε/2. Rackoff's version of GL list-decodes x with (L−1)·2^l queries,
  2^l = O(L/ε²), and q = 2^{L−1} + x is verified by division.
- Any 1/poly advantage therefore amplifies to factoring a 1/poly fraction of N.
- L = n/2 under the equal-length convention; otherwise try both candidate lengths.
- An earlier draft used k mod 2^L, which has an O(1/L) embedding error from prime k > q. Both reviewers flagged it,
  and the (L−1)-bit version removes it.
- Toy check (`gl.py`): a random-noise planted predictor of accuracy 0.55 recovers the factor 5/5 at n = 32 with
  15·2^12 queries. The code includes the empty subset J = ∅, which is harmless.
- *Meaning.* "One bit ⇒ factoring, with amplification" is possible for a function of M alone, because the reduction
  controls a free public part k of the instance that the bit mixes with the secret.

**codex's B(M)** = bit_{v₂(M)}(spf(M/2^{v₂(M)})). Queries M = 2^iN read p bit by bit. This is worst case, with no
amplification.

### 2.3 Status table (question 2)

| Quantity (from N alone unless stated) | Status | Reduction / reference |
|---|---|---|
| bit i of p, fixed nontrivial i (incl. bit 1 = p mod 4) | **OPEN** | none found |
| bit i of p, i an input | ≡ FACT | query all i |
| ⌊p/2^k⌋, k ≤ n/4 + O(log n) | ≡ FACT | Coppersmith J. Cryptology 10 (1997) Thms 4–5, plus 2^{O(log n)} enumeration; Rivest–Shamir 1985 needed n/3 |
| ⌊p/2^k⌋, n/4 + ω(log n) ≤ k < n/2 − O(1) | **OPEN** | topblock: lattice comparison theorem only |
| p mod M, log₂M ≥ n/4 + O(log n) | ≡ FACT | Coppersmith/May, known modulus (trackE E1) |
| p mod ℓ, ℓ small odd prime | **OPEN** | = (s mod ℓ) + a 1-bit label (trackE F2) |
| s = p+q mod ℓ, equivalently φ(N) mod ℓ | **OPEN**; decides "ℓ* is a square mod N" when (N/ℓ) = 1 | nonsmooth Prop RANK |
| φ(N) mod 2^k, k ≥ n/4 + v₂(q−p) | ≡ FACT (2 decodes) | trackE E2 |
| φ(N) mod 8, N ≡ 1 (4) | ⇔ ε | trackE E2′ |
| trivial residues: p mod 2, φ(N) mod 4, [2 \| φ(N)] | easy | — |
| ε = (−1/p), N ≡ 1 (4) | **OPEN**; = [−1 is a square mod N] = [r₂(N) = 16] | Seurin, PKC 2014 (2-Φ/4-hiding) |
| bit 1 of the smaller factor, N ≡ 3 (4) | **OPEN** (orientation bit) | — |
| [ℓ \| φ(N)], ℓ odd (Φ-hiding) | **OPEN**; conjectured hard on suitable distributions | Cachin–Micali–Stadler, EUROCRYPT 1999 |
| QR mod N, random base with Jacobi +1 | **OPEN** (classical) | — |
| fixed-base QR | **OPEN** in nontrivial cases | — |
| Ω(N) mod 2, ω, d, μ, squarefreeness, rad | trivial on the RSA promise (0, 2, 4, 1, yes, N); **OPEN** for general N | Adleman, FOCS 1994 (C7–C9: squarefree part, squarefreeness, ω), companion of Adleman–McCurley ANTS-I (LNCS 877) 1994, exact ANTS numbering UNVERIFIED; Shallit–Shamir 1985; Landau Inf. Comput. 78 (1988) and Woll Inf. Comput. 72 (1987) reduce *to* totient-type oracles, not from these (direction per codex; Du–Volkovich FSTTCS 2021 §1.1) |
| d(N) mod 2 | easy (square test) | — |
| φ(N), λ(N), σ(N); ord_N(a) with a selectable base | ≡ FACT (randomized) | Miller 1976; Bach–Miller–Shallit SIAM J. Comput. 15 (1986) |
| "N has a factor ≤ x", x an input | ≡ FACT | binary search, ≤ n calls |
| smallest / largest prime factor (of an arbitrary input) | ≡ FACT | strip and recurse |
| #E(Z/N), chosen curves | ≡ FACT (randomized, ≤ 6 calls expected via twists, codex) | Kunihiro–Koyama EUROCRYPT 1998 (with a distributional assumption) |
| #E(Z/N), one given curve | OPEN | — |
| r₂(N) | = 16·[ε = 1] on N ≡ 1 (4), **OPEN**; 0 on N ≡ 3 (4) | — |
| one sum-of-two-squares representation | OPEN; two inequivalent ones factor N | — |
| ⟨p, r⟩ with r public and uniform | hardcore (average case) | Goldreich–Levin |
| f(M) = parity(P⁺(M) ∧ M) on M = N·k | hardcore, loss-free embedding | Prop GLM |
| RSA/Rabin message bits | hardcore | ACGS; Håstad–Näslund |

**"Is there any theorem that one fixed bit of p is as hard as factoring?"**
- None found.
- All positive single-bit results take an extra input (r, x, y, i) or embed one in the instance (GLM, B(M)).

---

## 3. The middle ground: representations that leak (question 3)

Everything computed from N is a function of N. So "leak" is computational: does a poly-time R(N) predict a function
of p better than the *floor*, the best predictor from N mod m?

### 3.1 Classes with a proved floor
- **(A) R factors through N mod m.** PROVED (trackE F1; codex). This holds for any distribution of (p, q):
  I((p mod m, q mod m); R | N mod m) = 0. For m ≤ (log X)^A, Siegel–Walfisz makes the labelled pair asymptotically
  uniform on the fiber {(u, c/u)}, so the floor itself carries ≈ 0 bits about ε or about p's bits. The bound applies
  to lcm(m_i), not to each m_i separately.
- **(B) Abelian reciprocity aggregates** χ(N) and Artin symbols of (N) of fixed conductor. They reduce to (A)
  (nonsmooth §5.4). This covers products ∏χ(r)^{e_r}, not selections like χ(p) + χ(q).
- **(C) Multiplicative automatic sequences.** They reduce to (A) by KLM 2022.
- **(D) Low-height pure programs** (codex_1 §3). Take t tests, nonzero test integers of ≤ H bits, and a factor law
  that does not concentrate on a few primes (e.g. uniform on the primes of [X, 2X]). The transcript then equals the
  formal integer-zero transcript except with probability O(tH/X): a nonzero H-bit integer has ≤ H/log₂X prime
  divisors there. No leak.
- **(E) Pure trees of any height: Theorem PT.**

### 3.2 Theorem PT: in pure trees, information costs splitting

**Model (pure tree).**
- A finite decision tree T. Its structure and every program may depend on public parameters that are constant on the
  support of the instance law (e.g. n under the equal-length convention) and on coins drawn from an N-independent
  distribution, but on nothing else about N.
- Each node v carries an integer c_v, given by an integer straight-line program over Z (+, −, ×, integer constants,
  reuse of earlier values allowed).
- On input N the node computes c_v mod N. If 1 < gcd(c_v mod N, N) < N it outputs a *split*. Otherwise it branches
  on [c_v ≡ 0 (mod N)].
- Outputs depend only on the leaf.
- Inversion is formalised by numerator–denominator pairs: before using 1/b, test b (a split, or a failure branch if
  b ≡ 0 in both components). Equalities a/b = c/d become tests of ad − bc.
- Equivalently, classes are complete transcripts: leaves of the unrolled tree, not output labels. (The v1 code keyed
  classes by labels; codex found the counterexample c = 3 with both leaves labelled alike. Fixed, deviation D1.)
- **Covered**, when schedules (bounds B1/B2, step counts, curve and seed ranges) are fixed by the size class and
  coins come from an N-independent range: Pollard rho with a polynomial iterate, p−1, p+1, projective ECM,
  Pollard–Strassen/BGS products, and any SLP-based product test.
- **Not covered**:
  - sampling uniformly from Z/N;
  - reading integer representatives of residues;
  - programs depending on more of N: Fermat/MR exponents N − 1, √N, x² − N, forms of discriminant −4N, CF, lattices.
- **PT-m**: the root reads c = N mod m and picks a pure tree T_c.

**Partition lemma (PROVED).** For a prime x, its *diagonal transcript* answers [x | c_v] along the tree; class(x) is
that transcript. For distinct primes x, y, the instrumented run on N = xy splits iff class(x) ≠ class(y), and
otherwise ends with transcript class(x). On N = x² a non-split run also follows the diagonal transcript of x.

*Proof.* gcd(c_v, xy) ∈ {1, x, y, xy}, and it equals x or y iff [x | c_v] ≠ [y | c_v]. Otherwise the answer
[xy | c_v] equals [x | c_v] = [y | c_v]. By induction the run follows both diagonal transcripts while they agree,
and it splits at the first disagreement. For N = x²: gcd(c_v, x²) = x is a split; otherwise [x² | c_v] = [x | c_v].
∎ (Unit-tested, including codex's labelled-leaf counterexample and square inputs.)

**Theorem PT (PROVED).**
- Setting: T is deterministic (fix the coins); X, Y are iid with law μ on primes; F(X, Y) is any target. Let
  π_ℓ = μ(class ℓ), η = 1 − max π_ℓ and s₀ = 1 − Σπ_ℓ², and let s ≥ s₀ be the actual split probability of the
  instrumented run (on X = Y it can also split).
- Then s₀ ≥ η, and

  Pr[T outputs F(X, Y) or splits] ≤ max_b Pr[F = b] + s + η² ≤ max_b Pr[F = b] + 2s.

- For randomized T, η² is replaced by E_R[η_R²] ≤ E_R[s_R], and the last form holds with E[s].

*Proof.*
- Σπ² ≤ max π gives s₀ ≥ η.
- Every non-split run ends in the common class. So Pr[no split and class ≠ ℓ*] ≤ Σ_{ℓ≠ℓ*}π_ℓ² ≤ η·η, where ℓ* is the
  heaviest class.
- On {no split, class ℓ*} the output is a fixed value o(ℓ*), correct with probability ≤ max_b Pr[F = b]. ∎

**Corollary PT1.** Advantage δ (counting splits as successes) forces s ≥ δ/2. In this model, on iid factor laws, every
predicate of the factorisation is hardcore with loss 2 (relative to the trivial baseline). Weak predictors already
are factoring algorithms.

**Theorem PT-m (fiber floor; PROVED; the lemma is from the Opus critic).**
- Setting: X ~ μ_u and Y ~ μ_{u′} independent, with class laws α, β, and ℓ* = argmax α.
- Then Pr[no split, class ≠ ℓ*] = Σ_{ℓ≠ℓ*}α_ℓβ_ℓ ≤ 1 − α_max ≤ 1 − Σα_ℓβ_ℓ = s₀(u, u′).
- So, conditional on N ≡ c (mod m) (mixture over the fiber (u, c/u)):

  Pr[correct or split | c] ≤ E_{u | c}[ max_b Pr[F = b | u] + 2·s(u, c/u) ].

- Without splitting, a pure tree with side information N mod m learns at most the *fiber floor*: the best prediction
  from the unobserved fiber point u = p mod m.
- Example: ε is a function of u for m = 4, so this does not exclude ε.

**Corollary PT2 (what deciding ε costs; PROVED).**
- Population A (N ≡ 1 mod 4): ε = ±1 with probability ½, X, Y iid μ₁ (primes ≡ 1 mod 4) resp. iid μ₃.
- Population B: X ~ μ₁, Y ~ μ₃.
- Let α, β be the class laws of the tree T₁ that the algorithm uses on N ≡ 1 (mod 4), and let s_B = 1 − Σα_ℓβ_ℓ be
  the split rate of **T₁ forcibly run on population B**. This is not necessarily what the PT-m algorithm does on
  N ≡ 3 (mod 4), where it may run some other T₃ (codex's counterexample: T₃ = stop).
- Statement:

  Pr_A[correct or split] ≤ ½ + ½s_A + ½√(2 s_B),

  so advantage δ forces s_A ≥ δ or s_B > δ²/2.

*Proof (virtual-diagonal convention for X = Y).*
- Pr_A[…] ≤ ½Σ_ℓ max(α_ℓ², β_ℓ²) + s_A.
- ½Σmax(α², β²) = ½(1 − s_A) + ¼Σ|α² − β²|.
- Σ|α² − β²| ≤ ‖α − β‖₂‖α + β‖₂ ≤ 2‖α − β‖₂.
- ‖α − β‖₂² = 2s_B − 2s_A. ∎

**Meaning.** A pure tree can extract ε (or any fiber function) without splitting its own inputs only if the same tree
splits the cross population. In the ε case that population is the mixed moduli. The extreme example is
c = ∏_{x ≤ B}(x² + 1) with B ≥ (max supp μ − 1)/2. It decides ε exactly (s_A = 0) and splits every mixed N (s_B = 1).
Its cost is exponential, which PT ignores. The general two-law inequality ‖π^{(1)} − π^{(2)}‖₂² ≤ 2s(μ₁ × μ₂) holds
for any two factor laws.

**What PT is not.**
- It is not a running-time barrier. ECM is pure and heuristically finds p in L_p[1/2, √2]. Pure methods can be fast;
  they just learn nothing about the factors except by splitting.
- It is not a statement about non-pure algorithms. Jager–Schwenk's warning applies: computing the Jacobi symbol of a
  given element is generically as hard as factoring, yet it is easy by reciprocity.

**Literature.**
- The "component-separating zero ⇒ gcd factor" engine is standard in the generic-ring model.
- The closest prior statement (codex_2) is Jager–Schwenk, J. Cryptology 26 (2013) 225–245, Lemmas 3–4 and Thm 1. In
  a fixed-modulus, challenge-input setting, generic subset-membership advantage gives factoring with probability
  ≥ Adv/[2(m² + 4m + 3)].
- Also related: Aggarwal–Maurer EUROCRYPT 2009 §3.2.3.
- PT's input is only N, with the modulus varying over the instance law. The partition and fiber statements (PT, PT-m,
  PT2) are no exact prior formulation located by us.

**E1 (pre-registered; all predictions met).**
- Prime sets: all b-bit primes with the top two bits set, b = 14, 16, 18.
- Trees: nested p−1; adaptive p−1/p+1/ECM; three ECM curves; and at b ≤ 16 a Lipton interval product and full and
  partial ∏(x² + 1).
- **P1:** the partition lemma against real instrumented execution mod N gave 0 mismatches in 681 990 runs. These were
  all 92 665 pairs × 6 trees at n = 28, and 20 000 or 2 000 sampled pairs per tree at n = 32 and n = 36.
- **P2/P3:** the PT bound held 75/75 and the ε cross bound 15/15. These are exact numerical evaluations from the
  class partition (virtual diagonal, Bayes-optimal leaf labels), i.e. checks of proved inequalities, not of real runs.
- **Deviation D1.** After codex_2 the class key was changed from leaf labels to transcripts, and E1 was rerun
  (`factorbits_e1_rerun.jsonl`). All 15 records are identical, because every E1 tree has distinct leaf labels.
- **Descriptive.**
  - The theorem's slack, success − baseline − s, is negative in all 75 cells, so the bound holds with room.
  - The non-split information gain Pr[correct ∧ no split] − (1 − s)·baseline at n = 36 is:
    - 0.0027 for the adaptive p−1/p+1/ECM tree (s = 0.415);
    - 0.0057 for ECM (s = 0.733);
    - 0.0247 for nested p−1 (s = 0.434).
  - The p−1 gain comes from the smooth-order class correlating with p mod 4: p ≡ 1 (4) makes p − 1 smoother.
  - `qr_sniff_partial` at n = 32 has ε advantage 0.0356 = 0.0349 from splits + 0.0007 non-split.

### 3.3 What is open
- **(F) Finite automata and AC⁰ on the bits of N**, as predictors of factor bits of RSA moduli: OPEN.
  - Green (μ ⊥ AC⁰), Bernasconi–Shparlinski (STACS 1999, squarefreeness) and Allender–Saks–Shparlinski (JCSS 62
    (2001), primality ∉ AC⁰[r]) treat other targets.
  - By Linial–Mansour–Nisan, an AC⁰ floor for ε would follow from bounds on correlations of ε(p) with parities of
    polylog many binary digits of pq. That is a type-II bilinear digit sum, plausibly within reach of
    Mauduit–Rivat/Bourgain methods. Not attempted.
- **(G) Non-pure algebra**: Fermat/MR/Lucas tests with exponent N ± 1, Jacobi sums, integer representatives. OPEN.
  The eps battery found no predictor among 66 such statistics at |ρ| ≥ 0.022.
- **(H) Euclid-type and lattice structure of N**: CF of √N, forms of discriminant −4N, Coppersmith lattices, sieving.
  OPEN. All leakage seen in the project so far is relation leakage: the eps n = 32 CF hits, and nonsmooth's (U3)
  sums.

The sharpest statement is a trichotomy:
- the floor is **proved** for the local, abelian, automatic and low-height-pure classes;
- **leak ⇒ split** (PT, PT-m, PT2) holds for pure algebraic trees of any height;
- **open** for everything that uses N's integer structure beyond N mod m, as all known methods below L[1/2] do.

---

## 4. Leverage and self-reduction (question 4)

### 4.1 Minimal leakage known to factor
- **Coset-shaped bits of p** (top, low, p mod M): n/4 + O(log n). trackE measured n/4 + n/(4d) − 2 in dimension d.
  n/4 − t bits cost 2^t, so ω(log n) missing bits are not recoverable by this route.
- **Random positions of p and q**: density > 1/2 in the low n/4 positions, about n/4 bits in total (trackE E6;
  Heninger–Shacham CRYPTO 2009, heuristic).
- **Arbitrary yes/no queries**: εn for every fixed ε > 0, CONDITIONAL on an ECM-type conjecture (Maurer, Computational
  Complexity 5 (1995) 237–247). Rivest–Shamir (EUROCRYPT '85) needed n/3.
- **Noisy bits with independent errors at every position of p and q**: error < 0.084 (Henecka–May–Meurer CRYPTO 2010,
  heuristic). The capacity limit is 0.110 (Paterson–Polychroniadou–Sibborn ASIACRYPT 2012).
  - A per-bit advantage ε carries ≈ (2/ln 2)ε² bits, so ≈ n/4 bits of information need ε ≳ 0.3 per bit, under
    independence.
  - 1/poly advantages are far below every known decoder.
- **Lower side**: O(log n) bits of anything are worthless. A general threshold is OPEN.

### 4.2 Why ACGS-type amplification has no analogue for fixed bits of p
- **Where amplification comes from.**
  - RSA: c ↦ a^e c gives instances with secret ax.
  - Rabin: multiply by a².
  - DL: (g^x)^a g^b gives secret ax + b.
  - Many transformed instances give many *different* functionals of one secret.
- **What factoring lacks.**
  - The secret is the support {p, q}.
  - The instance maps N ↦ Nk and N ↦ N^j add or scale known vectors and leave the unknown support fixed.
  - Amplification would need a computable action that *moves* p.

**Lemma PSR (PROVED).** Let g, h ∈ Z[x] be fixed, with coefficients and degrees independent of N, g nonconstant and
h ≠ 0. Then g = c·x^j in each of the following cases:
- **(a)** For infinitely many primes q there are infinitely many primes p with g(p) | h(pq).
- **(b)** (balanced version, added after the Opus critic noted that (a) cannot hold for balanced reductions) For all
  sufficiently large primes p and all primes q ∈ (p, 2p), g(p) | h(pq). A fixed positive proportion of such q also
  suffices.

*Proof of (a).*
- Fix q. Write D·h(qx) = A(x)g(x) + R(x) with A, R ∈ Z[x] and deg R < deg g. Then g(p) | R(p).
- |R(p)| < |g(p)| for large p, so R has infinitely many roots. Hence R = 0 and g(x) | h(qx) in Q[x].
- A root α ≠ 0 of g gives infinitely many roots qα of h. ∎

*Proof of (b).*
- Write g = x^jG with G(0) ≠ 0. Let d = deg h, H = max|h_i| and C₀ = max(H, 8d + 2).
- Let p > |G(0)| be a large prime, and let ℓ > C₀ be a prime dividing G(p); then ℓ ∤ p.
- The polynomial y ↦ h(py) mod ℓ is nonzero: some h_i ≠ 0 has |h_i| < ℓ, and ℓ ∤ p. Its degree is ≤ d, so it has
  ≤ d roots. Yet it vanishes at every prime q ∈ (p, 2p). So those ≥ p/(2 ln p) primes lie in ≤ d residue classes
  mod ℓ. That is impossible for large p:
  - ℓ ≥ p: the primes are pairwise distinct mod ℓ.
  - √p < ℓ < p: each class holds ≤ p/ℓ + 1 integers of (p, 2p).
  - C₀ < ℓ ≤ √p: by Brun–Titchmarsh each class holds ≤ 4p/(φ(ℓ) ln p) primes, and φ(ℓ) > 8d.
- So for every large prime p, G(p) is C₀-smooth.
- But the values G(p) for primes p ≤ Y include ≥ π(Y)/deg G − O(1) distinct integers of size ≤ Y^{O(1)}, while there
  are only O((log Y)^{π(C₀)}) C₀-smooth integers of that size. So G is constant. ∎

**Scope.**
- Fixed polynomials only. Families that vary with the input length, randomized coefficient distributions that grow
  with N, and general poly-time maps are NOT covered.
- Codex's example h_a(t) = t + a with a = q has g = x + 1 dividing h_a(pq) for every p. It shows that randomizing
  over unbounded families escapes the lemma's quantifiers, though not that this yields an efficient reduction.
- The v1 claim "randomized polynomial maps are covered coefficient-wise" was wrong and is withdrawn.
- Divisibility sequences (a^N − 1 ∋ a^p − 1, Lucas u_N ∋ u_p, elliptic divisibility sequences) do move p, but at
  exponential size, and reducing them modulo anything loses the divisibility.

**What does amplify.**
- A *varying public mask*: Goldreich–Levin (⟨p, r⟩), Prop GLM (r embedded in N·k, loss-free).
- In the pure model there is nothing to amplify: any advantage is already a splitting probability (PT1).

**Answer to question 4.**
- A weak per-bit predictor for a FIXED natural bit position of the factor is not known to amplify. Even a perfect one
  is not known to give more than that bit.
- "Bit i of p on random N-related inputs ⇒ factoring" needs a poly-time family N ↦ N′ whose factors are known, varied
  functions of p other than p itself.
- Lemma PSR excludes fixed polynomial maps, in balanced and unbalanced form. The only natural examples
  (divisibility sequences) are exponentially large.
- OPEN. We have no candidate.

---

## 5. Candidate cards (rule 9)

| Mechanism | Model | Outcome |
|---|---|---|
| ADD/successor/partial-output oracles | oracle reductions | equivalent to FACT except the OPEN rows of §1.5; no mechanism for *computing* ADD faster |
| pure trees (p−1, p+1, ECM, rho, products) as ε / factor-bit predictors | Theorems PT, PT-m, PT2 | **refuted for this model**: information ⇒ split, quantitatively |
| polynomial self-reduction for bits of p | Lemma PSR (a), (b) | **refuted for this model** (fixed polynomials) |
| GL / cofactor-embedded bits | Prop GLM | works, but needs a free public component; not a natural bit |
| automaton / AC⁰ readers of N | — | not yet tested (OPEN; LMN route sketched) |

## 6. Most decisive open question
**Is there a polynomial-time family N ↦ N′ whose prime factors include known, varied functions of p other than p
itself?**
- A positive answer would give ACGS-style amplification for fixed factor bits.
- Divisibility sequences are the only natural examples, and they are exponentially large.
- In the pure algebraic world the question is moot (PT).
- Outside it, the first concrete target is still trackE §7: a non-pure poly-time observable with 1/poly advantage
  on ε.

## 7. Reproduction
```
.venv/bin/python -m pytest -q tests/test_factorbits.py                 # 51 tests, ~1 min
.venv/bin/python -m src.factorbits.run_experiments e1 14 16 18          # E1 -> experiments/raw/factorbits_e1.jsonl (~3 min)
FACTORBITS_RAW=experiments/raw/factorbits_e1_rerun.jsonl .venv/bin/python -m src.factorbits.run_experiments e1 14 16 18
```

## 8. Codex adversarial round (codex_2: MAJOR REVISION), applied

| Item | Response |
|---|---|
| model formalisation | done: transcripts not labels, finite trees, inversions, public parameters, squares via the actual split rate, randomized E[η²] |
| PT2 measures T₁ forced onto B | stated |
| GLM loss O(1/L) | removed by the (L−1)-bit embedding (also the Opus critic's fix) |
| fixed-m exponents | ≡ FACT by a randomized polynomial reduction (adopted, implemented, tested) |
| GRH successor row | corrected range T = N(log N)^C, C > 4; no factoring of t needed |
| randomized PSR claim | withdrawn; balanced version (b) added |
| N^{1/4} trichotomy | wrong (ECM is pure, heuristically L[1/2]); replaced |
| reduction directions, trivial residues, ord with a selectable base | fixed in §2.3 |
| low-height floor needs a non-concentrated law | stated in §3.1(D) |

## 9. Opus critics
- **Proof critic: ACCEPT WITH FIXES.**
  - It checked the partition lemma, PT, PT1, PT2, PSR(a), E1–E2, Z and the GL decoder.
  - It recomputed the E1 numbers and reran GL at ε = 0.05 (5/5).
  - All 16 items applied:
    - GLM embedding;
    - the character row's value group;
    - PSR for balanced moduli;
    - the randomized-PSR overclaim;
    - squares;
    - randomized trees;
    - algorithm-fit conditions;
    - the PT-m lemma;
    - E1 wording;
    - the slack quantity;
    - the ⌊log₂p⌋ bound;
    - the product example's B;
    - Las Vegas needs AKS;
    - the §1.3 clarification;
    - the fixed-m wording;
    - the gl.py notes.
- **Literature critic:** pending.
