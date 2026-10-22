# selfred: self-reductions of factoring, the open partial-information bits, and Theorem PT's open classes (phase 22)

Agent `selfred`, 2026-10-22.
- Pre-registration: `notes/selfred_prereg.md` (commit `86d5f7c`, before E1/E2 ran).
- Code: `src/selfred/{rigidity,reductions}.py` (OPERATIONAL), `src/selfred/{harness,census,run_experiments}.py`
  (HARNESS). Tests: `tests/test_selfred.py` (28 pass).
- Raw data: `experiments/raw/selfred_e1.jsonl`, `selfred_e1_posthoc.jsonl`, `selfred_e2.jsonl`.
- GPT-6.1-sol: blind round `codex/selfred_1.md`, full report `codex/selfred_1_report.md` (+ `selfred_1_check_claims.py`,
  output alongside); adversarial round `codex/selfred_2.{md,out}` (§8).
- Two read-only Opus critics (§9).
- Labels: **KNOWN** (with reference), **PROVED** (proof here; routine unless stated; "no exact prior formulation located"
  is the strongest novelty claim), **CONDITIONAL**, **HEURISTIC**, **OPEN**, **experimentally supported**,
  **REFUTED (prediction)**, **UNVERIFIED**.
- "OPEN" = we know no polynomial algorithm and no reduction *from* factoring. Negative literature findings come from our
  and codex's searches and are UNVERIFIED as negatives.
- "codex_1 §x" = `codex/selfred_1_report.md` §x; its proofs were re-checked here (§9 records the critics).

Notation: N = pq, p < q primes of equal bit length, n = log₂N, X ≍ √N the factor scale. ε = (−1/p) on N ≡ 1 (mod 4).
cont(h) is the content of h ∈ Z[x]. A *moved divisor* is g_r(p) | N′ with g_r a public, non-monomial function.

---

## 0. Verdict (plain language)

**No genuine self-reduction of factoring was found. Every algebraic way to build one either fails or already factors N.**
That is the structural reason why ACGS-style amplification has no analogue for fixed bits of p. It is a barrier for
algebraic maps, not a lower bound.

1. **N-oblivious polynomial maps cannot move p (Theorem SR1, PROVED).**
   - Take any h ∈ Z[x] of any height, with degree ≤ D. If G | h(pq) for a δ-fraction of the primes q ∈ (p, 2p), then
     G is p^a times a divisor of cont(h) times a smooth number S:
     - S has all prime factors < 1 + 4D/δ (Brun–Titchmarsh, here);
     - S has ≤ D(D+1)/2·log₂(O(D log p/δ)) bits (Vandermonde, codex).
   - A poly-time map has poly(n) output bits, so D ≤ poly(n) (Lemma DEG).
   - So randomized, growing-degree, arbitrary-height families never produce a large prime that moves with p.
   - E2 (pre-registered, exact): in all 3200 (h, p) cases the gcd minus its p-part equals the p-independent "unit fixed
     divisor" U(h), the same value for all 80 primes p.
2. **Rigidity (Theorem SR2, PROVED).** Suppose a poly-time map, N-dependent or not, outputs N′ with an exact public
   identity N′ = A(p)·B(q) for polynomials A, B of poly(n) degree and height.
   - Either A(x)B(N/x) is constant (A = a·x^j, B = b·x^j, the multiplier class), or root-finding over Z on one
     univariate polynomial factors N deterministically.
   - Torus orders (p+1)(q+1), (p²+p+1)(q²+q+1), shifted products (p+r)(q+r) and point counts with polynomial formulas
     are all instances.
   - The contrast with RSA: in RSA the secret moves inside Z/N, where root-finding is hard. A factor moves inside Z,
     where a second integer equation pins it down.
   - Escapes: non-algebraic moves (#E(F_p) = p+1−a_p, modular moves (a·p mod M) with unknown carries, bit operations)
     and moved divisors hidden in large unknown cofactors. #E(Z/N) is the one natural non-algebraic example, and for
     chosen twists it is factoring-equivalent (codex_1 §2.6).
3. **Meta-barrier (SR3, PROVED, trivial).** An unconditional proof that no such map exists, or that factoring does not
   reduce to ε or to squarefreeness, would prove factoring ∉ FP. So only restricted-model barriers are possible.
4. **Open bits.** ω, μ, λ, squarefreeness, rad, p mod 4, ε, φ(N) mod ℓ (fixed ℓ) and fixed-base QR all stay OPEN. We
   found no reduction from factoring and no separation. The status table (§2.1) is now precise:
   - rad ≡ squarefree decomposition ≡ parity kernel, by deterministic Turing reductions with ≤ log₂N calls (PROVED).
   - rad is also equivalent to maximal-order construction (Buchmann–Lenstra Thm 1.3).
   - With one d(·)-oracle call one gets μ, and with ≤ log₂N calls the whole signature (Shallit–Shamir). These
     reductions go *to* d, not from factoring.
   - Adleman FOCS 1994 C7–C9/O7–O9 are the squarefree part, squarefreeness and ω problems, still open.
   - φ(N) mod ℓ with ℓ an *input* is ≡ FACT: about 0.35n/ln n calls and CRT recover s = p+q (PROVED). For fixed ℓ it
     is OPEN.
   - Signature oracles queried on multiplicative constructions from N learn nothing. Unequal valuations of p and q
     split N by gcds (Lemma SIG). Arbitrary coprime queries remain open.
5. **Theorem PT's open classes.**
   - **Positional floor (Theorem L, codex; PROVED under GRH):** every function of N mod m, of any complexity, has
     advantage o(1) for ε (and for any function of p mod m) for all m ≤ N^{1/2−η}, under the independent-prime law.
   - **Unconditionally** the same holds for the low j bits of N with 2^j ≤ N^c (Bourgain's prescribed binary digits of
     primes).
   - So there is **no secrecy threshold at N^{1/4}**: the residue p mod m factors N above N^{1/4}, yet N mod m reveals
     nothing about it up to N^{1/2}.
   - E1 (exact census, n = 28, 32, 36): the Bayes-optimal advantage of N mod 2^j on ε obeys the proved bound in 93/93
     cells. Once the post-hoc correction (D2) is applied, it matches a random-label null.
   - Finite automata and AC⁰ remain OPEN. Codex gives a restricted extension for carry/Fourier-property sequences via
     Müllner's type-II bound (UNVERIFIED by us) and a precise LMN sufficient criterion.
6. **No positive leak** beyond what already factors N. The only honest advantage found is exponentially tiny: one Fermat
   step factors an Ω(N^{−1/4}) fraction of instances (codex_1 §5). That is a split, not a leak.

**Most decisive open question (refined).** Is there a poly-time map whose moved divisor is *non-algebraic* (like
#E(F_p)) or *hidden in an unknown cofactor*? These are the only two escapes from SR1–SR2. Concretely: is computing
#E(Z/N) for one fixed curve E (no twists) as hard as factoring, or easier?

---

## 1. Self-reductions (question 1)

### 1.1 What is being asked
- RSA blinding c ↦ c·r^e moves the secret x to x·r mod N by a public, efficiently invertible relation. Many such moves
  give many different functionals of x, which is what ACGS / Håstad–Näslund amplify.
- For factoring, the analogue is a poly-time (randomized, auxiliary-output) map N ↦ (N′, aux) with g_r(p) | N′ for a
  public non-monomial g_r. A *self-reduction to the same distribution* additionally needs N′ to be a balanced semiprime
  g_r(p)·h_r(q).
- Four requirements must be kept apart (codex_1 §2.1):
  - factor movement;
  - inversion-preserving self-reduction;
  - random self-reducibility over moduli;
  - information-theoretic instance hiding.
- All barriers below concern moved *divisors*, so they also apply to the stronger prime-factor requirement.

### 1.2 Candidate table

| Candidate | Status | Reason |
|---|---|---|
| kN, kN^j (multipliers) | **PROVED: monomial class only** | factors are p, q and those of k; GLM (factorbits §2.2) shows this suffices only for artificial predicates |
| h(N), h ∈ Z[x] N-independent (any height, poly degree, random coefficients) | **PROVED barrier (SR1)** | moved divisors are p^a · (divisor of cont h) · (O(D/δ)-smooth) |
| Φ_m(N) | **PROVED: SR1 applies**; p ∤ Φ_m(N) since Φ_m(0) = ±1 | bit length Θ(φ(m) n), so m ≤ poly(n) (codex_1 §2.2) |
| Norm_{K/Q}(α), α with coordinates polynomial in N, fixed K | **PROVED: SR1 applies** | the norm is a determinant, hence a polynomial in N |
| Norms with N-dependent algebraic data (e.g. NFS: f(m) ≡ 0 mod N) | **OPEN in general; observation** | p enters via N | Norm(θ − m) (multiplier) or via accidental gcd splits; other prime ideals are unrelated to p |
| Resultants Res(f, g) with a common root mod N | **PROVED: multiplier** | N \| Res; the cofactor is not a function of p |
| Discriminants with p \| disc f, q ∤ disc f | **PROVED: already splits (Lemma EUC)** | Euclid on (f, f′) over Z/N meets a zero divisor because deg gcd differs mod p and mod q (tested) |
| (p+1)(q+1), (p²+p+1)(q²+q+1), (p+r)(q+r), polynomial point counts | **PROVED: ≡ FACT (SR2)**; generalises nonsmooth §5.3 | exact algebraic identity ⇒ root-finding over Z |
| #E(Z/N) for chosen curves/twists | **KNOWN/PROVED: ≡ FACT (randomized)** | 4 twist counts sum to 4(p+1)(q+1) with prob ≥ 3/8 (codex_1 §2.6; Kunihiro–Koyama 1998 and later) |
| #E(Z/N), one fixed curve | **OPEN** | #E(F_p) = p+1−a_p is non-algebraic; SR2's approximate form fails (error a_p·q ≈ N^{3/4} ≫ s) |
| Lucas / torus recurrences mod N | **PROVED: residues ≠ orders** | evaluating u_k mod N does not give the integer group order |
| a^N − 1 ⊇ a^p − 1; Lucas u_N ⊇ u_p; EDS | **PROVED: exponential size**; reduction mod M destroys divisibility | the only natural moved divisors known |
| N′ = g(p)h(q) a balanced semiprime with algebraic g, h | **PROVED: ≡ FACT (SR2)** | exactly the "self-reduction to the same distribution" |

### 1.3 Theorem SR1 (oblivious polynomial maps cannot move p). PROVED

**Setting.** h ∈ Z[x] \ {0}, D = deg h ≥ 1, c = cont(h) > 0, p prime, and I = (x₀, x₀ + y] an interval with y ≥ 2. Let
Q ⊆ I be a set of primes with h(pq) ≠ 0 for q ∈ Q, and let G ≥ 1 satisfy G | h(pq) for all q ∈ Q.

**(a) Support bound (here).** Let ℓ be a prime with ℓ ≠ p and v_ℓ(G) > v_ℓ(c) (in particular any ℓ | G with ℓ ∤ c·p).
Then Q lies in at most D residue classes mod ℓ. Consequently:
- if ℓ ≥ y, then |Q| ≤ D;
- if √y ≤ ℓ < y, then |Q| ≤ D(√y + 1);
- if ℓ < √y, then |Q| < 4Dy/((ℓ − 1) ln y).

So if |Q| > D(√y + 1), every such ℓ satisfies ℓ < 1 + 4Dy/(|Q| ln y). If moreover |Q| ≥ δ·y/ln y, then ℓ < 1 + 4D/δ.

*Proof.*
- Write h = c·h₁ with h₁ primitive. Since v_ℓ(G) > v_ℓ(c), ℓ | h₁(pq) for every q ∈ Q.
- The polynomial f(Y) = h₁(pY) has coefficients h₁,ᵢpⁱ. Because ℓ ≠ p and h₁ is primitive, f mod ℓ is a nonzero
  polynomial of degree ≤ D. It vanishes at every q ∈ Q, so Q lies in ≤ D classes mod ℓ.
- Case ℓ ≥ y: the elements of I are pairwise incongruent mod ℓ, so |Q| ≤ D.
- Case √y ≤ ℓ < y: each class meets I in ≤ y/ℓ + 1 ≤ √y + 1 integers.
- Case ℓ < √y: use Brun–Titchmarsh in Montgomery–Vaughan's form: π(x₀+y; ℓ, a) − π(x₀; ℓ, a) < 2y/(φ(ℓ) log(y/ℓ))
  for ℓ < y (Mathematika 20 (1973), Thm 2). Here log(y/ℓ) > ½ ln y. ∎

**(b) Size bound (codex_1 §2.3, Theorem H; re-checked).**
- Let s = |Q| > D, and let R be the part of G/gcd(G, c) coprime to p.
- There are q₀ < … < q_D in Q inside a window of length H ≤ ⌊D·y/(s − D)⌋ with R | Δ = ∏_{i<j}(q_j − q_i). (Q may
  even consist of arbitrary integers here.)
- Hence R is H-smooth and log₂R ≤ D(D+1)/2·log₂H.
- *Proof:* the adjugate identity adj(V)·V·a = det(V)·a for the Vandermonde matrix V = (q_iʲ) and the coefficient
  vector a of h₁(pY), which has a unit coordinate mod every ℓ ≠ p. The window comes from pigeonholing the s − D blocks
  of D + 1 consecutive samples.
- With s ≥ δπ-density, H ≪ D log p/δ, which carries a log p factor. (a) avoids it for the support, and (b) adds the
  bit-length bound.

**Lemma DEG (output size forces degree; PROVED).**
- Let h ∈ Z[x] have degree D ≥ 1, and let |h(N)| < 2^T on D + 1 consecutive integers. Then D! < 2^{D+T}, i.e.
  D·log₂(D/(2e)) < T.
- *Proof:* Δ^D h ≡ D!·lead(h) is a nonzero multiple of D!. Also |Δ^D h(x)| = |Σᵢ(−1)^{D−i}C(D,i)h(x+i)| < 2^D·2^T. ∎
- So any polynomial map computed in poly time on all n-bit inputs, with outputs of T = poly(n) bits, has D ≤ poly(n),
  whatever circuit computes it. Exponential-degree circuits must have exponential output.

**Corollary SR1′ (self-reduction form).**
- Let (h_r) be N-independent polynomials (they may depend on n and on coins r) computed in poly time.
- Suppose public g_r satisfy g_r(p) | h_r(pq) for a δ-fraction of balanced pairs.
- By averaging, for a ≥ δ/2-fraction of p, the set Q_p of good partners has density ≥ δ/2 in the balanced range. That
  range is an interval of length y ≍ p containing ≍ p/ln p primes. So SR1(a) applies with δ′ ≍ δ.
- Hence every prime factor of g_r(p) is p, a prime of cont(h_r) (independent of p), or < 1 + O(D/δ) = poly(n).
- In particular no such family yields balanced semiprimes N′ = g_r(p)·h′(q) with g_r(p) a moved large prime.
- Fixing a seed is legitimate: the theorem is about one polynomial. Averaging only controls |Q_p|.

**Refinement: exact p-independence (PROVED for p ≥ p₀(h); E2).**
- Let u_ℓ(h) = min_{z: ℓ∤z} v_ℓ(h(z)) and U(h) = ∏_ℓ ℓ^{u_ℓ(h)}. Only primes ℓ ≤ D+1 or ℓ | c contribute: otherwise a
  nonzero polynomial of degree ≤ D < ℓ − 1 misses some unit.
- If Q = all primes in (p, 2p) and p ≥ p₀(h), then G_p := gcd_{q∈Q} h(pq) = p^{v_p(G_p)} · U(h)/p^{v_p(U(h))}.
  - "≥" always holds.
  - "≤": the minimum valuation over units is attained on a class mod ℓ^{u_ℓ+1}, and the set {pq mod ℓ^{u_ℓ+1}} meets
    every unit class once Siegel–Walfisz applies (ℓ^{u_ℓ+1} ≤ (log p)^A).
  - For ℓ | c with ℓ − 1 > D, |Q| > D(y/ℓ + 1) suffices.
- The smooth part of G is then literally the same integer for every p: the factors never move.
- Codex's caveat stands: for a *subset* Q (density δ < 1), the small smooth part can carry residue information about p.
  For example, h = x − 1 with Q = {q ≡ 1 mod 4} detects p mod 4. This is ≤ O(log(1/δ)·log(D/δ)) bits of small-prime
  data, not a moved large factor.

**E2 (pre-registered, exact).**
- 40 polynomials: D₀ ∈ {1, 2, 3, 5, 8}; coefficients up to 2^64; content up to 3 primes < 2^20; half planted with
  x·∏_{i≤D₀}(x − i), so total degree ≤ 17.
- 80 primes p (40 after 10⁴, 40 after 10⁵), with Q = all primes in (p, 2p). That gives 3200 cases.
- **P1 (proved bound):** held in 3200/3200 cases. The hypothesis |Q| > D(√p + 1) held in 2880 of them; in the rest the
  bound is vacuous. The largest prime ℓ ∤ c dividing G_p was 7, against a minimum bound L = 5.14 among hypothesis cases.
  The worst ratio ℓ/L was 0.39.
- **P2:** G_p/p^{v_p} = U(h) in 3200/3200 cases. Each h produced exactly one value across all 80 primes p.
- Runtime 11 s, RSS 80 MB.

### 1.4 Theorem SR2 (rigidity: exact algebraic self-reductions already factor). PROVED

**Statement.**
- Let N = pq with p ≠ q primes. Let A, B ∈ Z[x] be nonzero, with degrees a, b and coefficients of ≤ H bits; they may
  depend on N arbitrarily.
- Let N′ ∈ Z with N′ = A(p)·B(q) (or the same with p, q swapped).
- Put P(x) = A(x)·Σ_j B_j N^j x^{b−j} − N′x^b ∈ Z[x]. Then P(p) = 0, deg P ≤ a + b, and the coefficients of P have
  O(H + b·n + log|N′|) bits.
- Either (i) P ≡ 0, which happens iff A = αx^j and B = βx^j for the same j (and then N′ = αβN^j), or (ii) the integer
  roots of P include p. They can be found deterministically in time poly(a + b, H, n, log|N′|), which factors N.

*Proof.*
- P ≡ 0 iff the Laurent polynomial A(x)B(N/x) is the constant N′.
- If A had a root α ≠ 0, then A(α)B(N/α) = 0, forcing N′ = 0 and A(x)B(N/x) ≡ 0. That is impossible for nonzero A, B.
- So A = αx^i. Symmetrically, a nonzero root β of B gives a zero at x = N/β, so B = βy^k. The product αβN^k x^{i−k} is
  constant iff i = k.
- In case (ii), integer roots are found by exact real-root isolation or by LLL factorisation over Q (LLL 1982). Our
  implementation is `rigidity.integer_roots`: recursive sign-change bracketing of the squarefree part, using only
  integer evaluations. Each root r is tested for r | N. ∎

**Extensions and exceptions** (codex_1 §2.4; re-checked):
- **Rational functions.** The degenerate case is "A(x)B(N/x) ≡ const after cancellation", which is broader than the
  monomial case. Example: A = x + 1, B(y) = y/(N + y). The general statement applies to the whole rational expression
  R_N(X) = U/V. If R_N is nonconstant, F = U − N′V is a nonzero polynomial with root p.
- **Cofactors.** An unknown cofactor K with N′ = K·A(p)B(q) is handled if K takes ≤ poly(n) values (enumerate them).
  The same holds for an additive error with ≤ poly(n) candidate values.
- **Approximation.** Relative approximation alone is NOT covered: N′ = N approximates (p+1)(q+1) to relative error
  O(N^{−1/2}).
- **Escapes.** A modular move u ≡ ap (mod M) has an unknown carry, exponentially many possibilities. Floors, bit
  operations, non-algebraic g (point counts with a_p) and large unknown cofactors also escape. No barrier is claimed
  for them.
- **Tests.** `test_rigidity_factors` covers 12 seeds at 64–128 bits, degrees ≤ 4, 40-bit coefficients and torus maps,
  with both orientations. `test_rigidity_degenerate_multiplier_class` covers the monomial case.

**Why RSA escapes (PROVED; codex_1 §2.5 agrees).**
- The RSA move's relation (rx − kN)^e ≡ c·r^e (mod N) is a consequence of x^e ≡ c. It is not a new integer equation,
  and lifting it to Z introduces an unknown carry k.
- Root-finding modulo a composite of unknown factorisation is the hard step. For e = 2 it is factoring itself (Rabin).
- In factoring, the instance is the integer equation pq = N. A second, independent, exact algebraic integer equation
  reduces to univariate root-finding over Z, which is easy. Multipliers give the dependent equation kN = k·pq.
- This explains the algebraic asymmetry. It does not exclude a group action on some *other* encoding of the factoring
  secret.

### 1.5 SR3 (meta-barrier). PROVED, trivial (codex agrees)
- If factoring ∈ FP, then N ↦ (p + r)(q + r) is a poly-time moved-divisor map, and every oracle problem trivially
  reduces from factoring with zero queries.
- So an unconditional proof that (i) no self-reduction exists, or (ii) factoring ≰ ε / squarefreeness / μ / …, implies
  factoring ∉ FP (resp. ∉ BPP for randomized reductions), hence P ≠ NP.
- Barriers can only be restricted-model (SR1, SR2, PT) or relativized.

### 1.6 Random self-reducibility and instance hiding: what is known (KNOWN; details codex_1 §2.7)
- **Abadi–Feigenbaum–Kilian**, JCSS 39 (1989) 21–50. §6 *poses* as Question 1 whether factoring is encryptable, i.e.
  random-self-reducible. Their negative results concern NP-hard functions.
- **Feigenbaum–Fortnow**, SICOMP 22(5) (1993) 994–1005.
  - Thm 3.1: an NP language with a nonadaptive random self-reduction lies in AM_poly.
  - Cor. 3.3: a PH collapse for NP-complete sets.
  - Factoring's decision version ("N has a prime factor ≤ b") is in NP ∩ coNP (Pratt certificates). The NP-completeness
    hypothesis therefore fails, and FF says nothing about factoring.
- **Beaver–Feigenbaum**, STACS 1990: multi-oracle instance hiding exists for every function. It is not a map to
  independent factoring instances.
- **Bogdanov–Trevisan**, FOCS 2003 / SICOMP 36 (2006): an NP-complete barrier for nonadaptive worst-to-average
  reductions. Nothing about factoring.
- **Per-modulus RSR (Rabin 1979; Goldwasser–Micali 1984).** For a fixed N, square roots and QR randomize over residues
  y ↦ y·r². An α-fraction square-root oracle factors that N with probability α/2. This is the only "partial random
  self-reducibility" of factoring we found: over residues, not over moduli.
- **Not located:** any worst-to-average reduction that changes the modulus, including for Blum integers. Also not
  located: any result giving random self-reducibility "for smooth parts". Randomising known smooth cofactors is the
  multiplier class.
- **Accurate statement:** factoring is not known (to us) to be random self-reducible over the modulus. AFK list this as
  open.

---

## 2. The open partial-information bits (question 2)

### 2.1 Status table (oracles on their full stated domain; "≤" = poly-time Turing reduction)

| Quantity | Value on N = pq | Status | Known reductions / references |
|---|---|---|---|
| ω(M) | 2 | **OPEN** | ω ≤ d (signature interpolation, ≤ log₂M calls; Shallit–Shamir IPL 20 (1985) 151–153). Adleman FOCS 1994 C9/O9 |
| μ(M) | 1 | **OPEN** | μ ≤ d with **one** call: D = d(M^{r−1}), r > log₂M + 1 prime, is a pure r-power iff M is squarefree (Shallit–Shamir; proof codex_1 §3.3). μ = λ·SF |
| Liouville λ(M) | 1 | **OPEN** | λ ≤ d (signature); not the factoring-equivalent Carmichael λ |
| squarefreeness SF(M) | yes | **OPEN** (Adleman FOCS 1994 C8/O8) | SF ≤ φ deterministic (Woll, Inf. Comput. 72 (1987)); SF ≡ maximality testing of orders (Buchmann–Lenstra JTNB 6 (1994) Thm 6.12). Not in AC⁰[ℓ] for any prime ℓ (Allender–Saks–Shparlinski JCSS 62 (2001) Cor. 2); superpolynomial AC⁰ size (Bernasconi–Shparlinski STACS 1999). Booker–Hiary–Keating (Duke 164 (2015)): under GRH, Prop. 2.3 certifies squarefreeness or finds a square factor in exp((log M)^{θ*+o(1)}) with θ* = 2/3 HEURISTIC; Thm 1.1 proves RSA-210/220/230/232 *not squarefull* (not squarefree); rigorous elementary M^{1/6+o(1)} |
| rad(M) | N | **OPEN** | rad ≡ squarefree decomposition ≡ parity kernel s(M) (M = u²s), deterministic, ≤ ⌊log₂M⌋ calls (Prop RAD, PROVED; `squarefree_decomposition_via_rad`, tested). rad ≡ ring of integers (Buchmann–Lenstra Thm 1.3; Chistov 1989). rad ≤ φ, ≤ Carmichael λ (Landau Inf. Comput. 78 (1988)). Adleman FOCS 1994 C7/O7a (algorithm?), O7b (randomized reduction from factoring?) |
| p²q special case | rad = pq, s = q | these *values* factor p²q; computing them OPEN | squarefreeness on {p²q} ∪ {pqr} OPEN; nothing polynomial known for r = 2 |
| p mod 4 (smaller factor) | — | **OPEN** | = ε on N ≡ 1 (4); orientation bit on N ≡ 3 (4) |
| ε = (−1/p), N ≡ 1 (4) | ±1 | **OPEN** | = QR_{−1}(N) = [φ(N) ≡ 0 mod 8] = [N is a sum of two squares]; Seurin PKC 2014 (2-Φ/4-hiding). h(−4N) mod 8 does **not** determine ε (Rédei; 4 counterexamples, codex_1 §3.7) |
| φ(N) mod ℓ, ℓ fixed odd prime | s mod ℓ | **OPEN**; decides ℓ* QR mod N when (N/ℓ) = 1 (nonsmooth RANK) | Φ-hiding [ℓ \| φ(N)]: Cachin–Micali–Stadler 1999, an assumption |
| φ(N) mod ℓ, ℓ an **input** prime | — | **≡ FACT (PROVED)** | CRT over odd primes with ∏ℓ > 2^{⌈n/2⌉+2} gives s exactly; ≈ 0.35n/ln n calls (`factor_via_phi_mod_ell`, tested) |
| φ(N) mod 2, 4 | 0, 0 | trivial | — |
| fixed-base QR_a(N) | — | **OPEN** for nontrivial a | a function of p mod 4\|a\| on the fiber; random-base QR = Goldwasser–Micali assumption |
| d(M) | 4 | OPEN for general M | Shallit–Shamir: d is equivalent to the exponent signature |

**Two corrections to earlier project notes.**
- Allender–Saks–Shparlinski do cover squarefreeness: TC⁰ ⊆ AC⁰[Square-Free] (Thm 1) and SF ∉ AC⁰[ℓ] (Cor. 2). The
  factorbits wording "primality ∉ AC⁰[p]" undersold it.
- The Boneh–Durfee–Frankel residue result is **Corollary 1** (p. 28) of the ASIACRYPT 1998 paper, not "Cor. 2.2"
  (codex_1 §4.4). The factorbits row should be read accordingly.

### 2.2 Lemma SIG (signature oracles and multiplicative queries). PROVED
- For any nonzero X computed from N, strip copies of N while gcd(X, N) = N.
  - If at some step 1 < gcd < N, N is split without the oracle.
  - Otherwise X = N^a·Y with gcd(Y, N) = 1, so sig(X) = sig(Y) ⊎ {a, a}.
- `valuation_split`, tested.
- Consequence (codex_1 §3.4): queries built by products and powers from N and fully factored public integers either
  split N directly or have a fully known signature. Signature oracles (ω, μ, λ, SF, d) then give nothing.
- A useful reduction must query integers Y coprime to N whose *factorisation pattern* depends on p. That is a
  self-reduction question again (§1). For Y = h(N) with N-independent h, SR1 forbids moved large primes. It does
  **not** forbid p-dependence of the signature of Y in general, so this is a reduction of the question, not a barrier.

### 2.3 Lemma EUC (discriminants). PROVED, standard
- Monic Euclid in (Z/N)[x] either meets a non-unit leading coefficient, which gives a gcd split, or computes h with
  h mod p = gcd mod p and h mod q = gcd mod q, so the degrees agree.
- So any constructed f with p | disc f and q ∤ disc f (f monic) splits N.
- `poly_gcd_mod_N`, tested on f = (x − r)(x − r − p).

### 2.4 Separations
- No relativized/oracle separation between factoring and any of these natural functions was found by us or codex.
- By SR3, unrestricted non-reducibility is out of reach.
- Restricted-model "separations" are: PT (generic ring: ε-advantage ⇒ splitting), Lemma SIG (multiplicative queries)
  and SR1 (oblivious polynomial queries cannot move p). None is an oracle separation of the real integer problems.

---

## 3. Theorem PT's open classes (question 3)

### 3.1 Positional floor: functions of N mod m

**Theorem L (codex_1 §4.2; PROVED under GRH; proof re-checked).**
- Setting: P, Q iid uniform primes in [X, 2X], P ≠ Q; for ε, condition on PQ ≡ 1 (mod 4); L = lcm(4, m) < 4X.
- Under GRH mod L: TV(law of P mod L, uniform on units) ≪ √(φ(L)/X)·log X·log(XL).
- Consequently every predictor F(N mod m), of any complexity and possibly randomized, has advantage on ε at most
  O(√(m/X) log²X) + O(log X/X).
- The same bound holds for any target f(p mod m) relative to its best constant; ordering is handled by fixed-length
  sub-interval bins.
- Proof:
  - GRH gives |μ̂(χ)| ≪ log X log(XL)/√X for each nonprincipal χ (Montgomery–Vaughan, Multiplicative Number Theory I,
    Thm 13.7 / (13.21)).
  - Parseval gives Σ_a |μ(a) − 1/φ|² = φ^{−1}Σ_{χ≠χ₀}|μ̂(χ)|². Cauchy–Schwarz turns this into TV.
  - In the uniform model the pair (U, V) given UV = c is uniform on the fiber, and χ₄(U) is a fair coin independent of
    c (trackE F1).
- **Range:** m ≤ X^{1−η} = N^{1/2−η}.
- This supersedes our first, pointwise argument, which reached only m ≤ N^{1/4}/poly(n). That argument produced a
  spurious "meeting" with the Coppersmith threshold.

**Unconditional low bits (codex_1 §4.3; PROVED from a KNOWN theorem).**
- Bourgain (Israel J. Math. 206 (2015) 165–182, §1 Theorem, (1.1)–(1.2)): primes with r < γb prescribed binary digits
  have the asymptotic density 2^{−r}.
- Hence the law of p mod 2^j is o(1)-close to uniform for j ≤ γ′·n/2.
- So any function of the low j bits of N has advantage o(1) on ε for 2^j ≤ N^c, with an unspecified absolute c > 0.
- Barban–Linnik–Tshudakov (Acta Arith. 9 (1964), range x^{3/8}) and Gallagher (Invent. Math. 16 (1972), x^{2/5}) are
  stated for fixed **odd** prime powers in the checked secondary source (Baier–Pujahari, Acta Arith. 204 (2022)). We do
  not use them for 2^j. (Our pre-registration had cited BLT for 2^j; corrected.)

**No secrecy threshold at N^{1/4} (PROVED distinction).**
- p mod M with log₂M ≥ n/4 + O(log n) factors N together with N (BDF Cor. 1; trackE E1).
- N mod m reveals essentially nothing about p mod m up to m ≈ N^{1/2−η} (GRH).
- These are different experiments, because the fiber point is not given. The N^{1/4} threshold belongs to the hint,
  not to the residue of N.

**E1 (pre-registered, exact census).**
- Population: ordered iid pairs from S_b = all primes in [3·2^{b−2}, 2^b), b = 14, 16, 18 (n = 28, 32, 36), with
  |S_b| = 431, 1491, 5296.
- For every j ≤ n we compute the exact Bayes-optimal advantage A*_j of any function of N mod 2^j on ε (restricted to
  N ≡ 1 mod 4), the exact TV₁ and TV₂, and the proved bound 3TV₂/(1 − 2TV₂).

| b | j | A*_j | null (dedup, D2) | TV₁ |
|---|---|---|---|---|
| 16 | 4 | 6.7e−4 | 1.1e−3 | 0.012 |
| 16 | 8 | 1.8e−3 | 4.3e−3 | 0.077 |
| 16 | 12 | 1.1e−2 | 1.7e−2 | 0.417 |
| 16 | 16 | 6.5e−2 | 6.8e−2 | 0.954 |
| 16 | 20 | 0.255 | 0.256 | 0.997 |
| 16 | 24 | 0.470 | 0.470 | 1.000 |
| 18 | 8 | 1.9e−3 | — | 0.038 |
| 18 | 12 | 3.3e−3 | — | 0.198 |
| 18 | 19 | 5.3e−2 | — | 0.980 |

- **P1 (proved bound): 93/93 cells hold.** The bound is non-vacuous only while TV₂ < 1/2, i.e. 2^j ≲ |S|/8. There the
  empirical TV₁ scales like √(2^{j}/|S|), which is the random-model shape of Theorem L's √(φ/X).
- **P2: REFUTED as registered (prediction error, not a leak).**
  - (a) "A* ≤ 0.02 for j ≤ n − 2log₂b − 2" ignored finite-population overfitting. Once classes hold few pairs, even
    random labels give a large Bayes advantage. The registered null itself exceeds 0.2 there.
  - (b) "ratio to null ∈ [0.5, 2]" failed for b = 14, j = 2–4 (ratios 5.3, 3.8, 2.7). There A* is the population's
    global ε imbalance (#primes ≡ 1 vs 3 mod 4 in S₁₄; a single class at j = 2), whose size is set by TV₁ and is
    covered by P1.
  - The registered null also double-counted (X, Y) and (Y, X), which share N and ε.
- **Post-hoc (D2, labelled):** with one coin per unordered pair, A*/null lies in [0.42, 1.0] for every j ≥ 6 at b = 14
  and every j ≥ 3 at b = 16, and equals 1.00 (to two decimals) for j ≥ 19. The excess at the smallest j (b = 14:
  3.8, 2.7, 1.9, 1.3 for j = 2–5) is the global ε imbalance of the finite prime set, bounded by P1. Beyond those j, A* sits at or below random-label fluctuation. Below the null is
  plausible, because ε(X) = ε(Y) structure makes classes more balanced.
- **Descriptive:** j_{0.05} = 13, 16, 19 for n = 28, 32, 36. That is ≈ n/2 − 1, the onset of finite-population
  fitting, not of information.
- Nothing in E1 is evidence of ε leakage from low bits. It is consistent with Theorem L extending, empirically, to all j
  where classes are populated.

### 3.2 Automata, AC⁰, integer-structure algorithms

| Class | Status |
|---|---|
| multiplicative automatic sequences | floor (KLM 2022 Thm 1: on large semiprimes a(N) = f₂(N), a function of N mod a fixed period) — **PROVED**, codex_1 §4.5 |
| sequences with Mauduit–Rivat carry + Fourier properties, linear Fourier decay | **codex claim**: Müllner's type-II bound (Duke 166 (2017), arXiv v3 Prop. 5.5, (29)) with a_u = 1_prime·χ₄, b_v = 1_prime gives o(1) correlation with ε on balanced semiprimes. **UNVERIFIED by us** (hypotheses not re-checked); not all automata satisfy them |
| all finite automata | **OPEN** |
| AC⁰ | **OPEN**. Sufficient criterion (codex_1 §4.6, PROVED): \|E_S[εF]\| ≤ √(Σ_{j≤t}C(b,j))·β_t + 2√(η/ρ), where β_t = max over parities of ≤ t bits of the ε-twisted semiprime correlation, η = LMN tail, ρ ≍ 1/log² the slice density. Needs β_t ≤ b^{−ω(t)} uniformly, a bilinear digit problem Σχ₄(p)e(b·pq/2^k) with denominators up to the full range. Not proved |
| AC⁰[ℓ] for SF itself | KNOWN lower bound (ASS) — but SF is constant on semiprimes; no implication for ε |
| CF of √N, forms of disc −4N, lattices, Fermat/MR with exponent N ± 1 | **OPEN**; all observed leakage is relation/splitting leakage (eps battery; §4 below) |

---

## 4. Positive leaks (question 4)

- **Fermat one-step (codex_1 §5.1, PROVED).** With a = ⌈√N⌉, a² − N is a square whenever q − p ≤ √X. That happens on
  an Ω(N^{−1/4}) fraction of independently sampled same-class pairs (Cauchy–Schwarz over √X-bins; no prime-pair
  conjecture).
  - On a fair mixture of the ε = ±1 populations, "factor if Fermat succeeds, else guess" has advantage Ω(N^{−1/4}) on ε.
  - It is a split, it is exponentially small, and it does not amplify.
  - Artefact check (codex): independent sampling, hidden label, verified factors.
- **E1/E2** are exact censuses that read the factors by design (HARNESS). No operational code sees p.
  - `rigidity.py` and `reductions.py` import neither sympy, flint, harness nor census (AST test).
  - Their oracles are passed in from `harness.py`.
- **No inverse-polynomial non-splitting predictor** of any open bit was found.

---

## 5. Candidate cards (rule 9)

| Mechanism | Model | Outcome |
|---|---|---|
| polynomial / norm / cyclotomic / resultant self-reduction | N-oblivious h ∈ Z[x], any height, poly degree | **refuted for this model** (SR1 + DEG): no moved large prime; smooth part p-independent for full Q (E2) |
| algebraic self-reduction to semiprimes, torus/point-count orders | exact public identity N′ = A(p)B(q) | **refuted as an amplifier**: it is already a factoring algorithm (SR2) |
| discriminant-based moves | p \| disc f, q ∤ disc f | **refuted**: Euclid splits (EUC) |
| signature oracles on multiplicative queries | products/powers of N and factored integers | **refuted for this model** (SIG) |
| low bits / residues of N as ε predictors | any function of N mod m | **refuted** for m ≤ N^{1/2−η} (GRH) and 2^j ≤ N^c (unconditional) |
| #E(Z/N) one curve, non-algebraic / hidden-cofactor moves | — | **not yet tested / OPEN** |

## 6. Most decisive open question

The question is no longer "is there any poly-time N ↦ N′ moving p". SR1 and SR2 show that the algebraic versions are
either impossible (oblivious) or already factoring algorithms (exact). The remaining room is:
- (a) **non-algebraic moved divisors.** The prototype is N ↦ #E(Z/N) for one fixed curve, whose moved factor
  #E(F_p) = p+1−a_p varies with p but satisfies no polynomial identity. Is this one-curve count factoring-equivalent? An
  algorithm for it would be the first factor-moving self-reduction.
- (b) **moved divisors inside large unknown cofactors** (divisibility-sequence-like) of polynomial size.

A barrier for (a) would need a model of "non-algebraic" maps. We have none.

## 7. Reproduction
```
.venv/bin/python -m pytest -q tests/test_selfred.py                       # 28 tests, ~1 s
.venv/bin/python -m src.selfred.run_experiments e1 14 16 18                # E1 -> experiments/raw/selfred_e1.jsonl (~90 s; b=18 peak RSS 3.5 GB, see D3)
.venv/bin/python -m src.selfred.run_experiments e1posthoc 14 16            # D2 -> selfred_e1_posthoc.jsonl
.venv/bin/python -m src.selfred.run_experiments e2                         # E2 -> experiments/raw/selfred_e2.jsonl (11 s)
```

**Deviations (disclosed).**
- **D1.** The E1 null uses de Moivre's closed form exactly for class sizes k ≤ 4000 and log-gamma double precision
  above. It is a reference only; P1 is exact.
- **D2.** A post-hoc unordered-pair null (`e1_posthoc_dedup`), because the registered null double-counts orderings.
- **D3.** E1 at b = 18 peaked at 3.5 GB RSS, over the 2 GB guideline. It ran alone, with memory 80% free beforehand,
  and finished in 85 s. A chunked rewrite would avoid this; not redone.
- **D4.** The pre-registration cited BLT for 2-power moduli. BLT and Gallagher are stated for odd prime powers, so the
  unconditional statement uses Bourgain instead.

## 8. Codex adversarial round
PENDING.

## 9. Opus critics
PENDING.
