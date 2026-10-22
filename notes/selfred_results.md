# selfred: self-reductions of factoring, the open partial-information bits, and Theorem PT's open classes (phase 22)

Agent `selfred`, 2026-10-22.
- Pre-registration: `notes/selfred_prereg.md` (commit `86d5f7c`, before E1/E2 ran).
- Code: `src/selfred/{rigidity,reductions}.py` (OPERATIONAL), `src/selfred/{harness,census,run_experiments}.py`
  (HARNESS). Tests: `tests/test_selfred.py` (33 pass).
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

**We found no self-reduction of factoring that changes the modulus, and no new polynomial-time predictor of any open
bit that works without splitting N.** Two precisely specified algebraic mechanisms now come with rigorous
restrictions:
- a polynomial shared across many partners q (SR1);
- an exact public algebraic equation (SR2).

This explains why the *algebraic* routes to ACGS-style amplification for fixed bits of p fail. It is not a lower bound
and not a classification of all maps. (v2: narrowed after codex_2 MAJOR REVISION and the Opus proof critic; §8–9.)

1. **Shared-polynomial maps cannot move a large prime (Theorem SR1, PROVED).**
   - Take h ∈ Z[x] of any height and degree D, independent of N (it may depend on n and on coins). Suppose G | h(pq)
     for a set Q of primes q ∈ (p, 2p) with |Q| ≥ δ·κ_p·p/ln p and |Q| > D(√p + 1).
   - Then G = p^a · gcd(G, cont h)_{p′} · S, with S := (G/gcd(G, c))_{p′}, and:
     - every prime factor of S is < 1 + 4D/(κ_pδ) (Brun–Titchmarsh; here);
     - log₂S ≤ D(D+1)/2 · log₂(D p/(|Q| − D)) (Vandermonde; codex).
   - Here κ_p = (π(2p) − π(p)) ln p/p → 1.
   - Poly-time maps are covered after replacing h by its remainder modulo the vanishing polynomial of the input range,
     which has degree poly(n) (Lemma DEG).
   - E2 (pre-registered; exact; rerun with the fixed code): in 3200/3200 (h, p) cases, G/p^{v_p(G)} equals the
     p-independent unit fixed divisor U(h), the same integer for all 80 primes p.
2. **Rigidity (Theorem SR2, PROVED).** Suppose a map (N-dependent or not) outputs N′ and public polynomials A, B of
   poly(n) degree and height with N′ = A(p)·B(q) exactly.
   - Then either A = a·x^j and B = b·x^j (the multiplier class), or exact root-finding over Z on one univariate
     polynomial factors N deterministically.
   - Instances: torus orders, (p+r)(q+r), and polynomial point counts.
   - In RSA the secret moves inside Z/N, where root-finding is hard. Here a second exact *integer* equation pins p down.
   - Not covered: non-algebraic moves (#E(F_p) = p+1−a_p), modular moves with unknown carries, implicit or
     exponential-size representations, N-dependent divisibility constructions, and divisors hidden in large unknown
     cofactors. No claim is made that these are the only escapes.
3. **Meta-point (SR3, PROVED, trivial).** An unconditional proof that no moved-divisor map exists, or that factoring does
   not reduce to ε or to squarefreeness, would prove factoring ∉ FP (resp. BPP), hence P ≠ NP. So the barriers within
   current reach are restricted-model ones.
4. **Open bits: all still OPEN.** This covers ω, μ, Liouville λ, squarefreeness, rad, p mod 4, ε, φ(N) mod ℓ (fixed ℓ)
   and fixed-base QR. We found no reduction from factoring and no separation (§2.1).
   - rad ≡ squarefree decomposition ≡ parity kernel (Adleman's C7) by deterministic Turing reductions with ≤ log₂N calls.
   - rad ≡ maximal orders (Buchmann–Lenstra).
   - μ needs one call to d(·); the whole prime signature needs ≤ log₂N calls (Shallit–Shamir). These go *to* d.
   - φ(N) mod ℓ with ℓ an *input* is ≡ FACT, by CRT recovery of s = p + q (≈ 0.35n/ln n calls asymptotically). For a
     fixed ℓ it is OPEN.
5. **Theorem PT's open classes.**
   - **Positional floor (Theorem L, codex; PROVED under GRH).** Setting: P, Q independent uniform primes in [X, 2X], a
     fixed public modulus m, and only N mod m observed. Then every function of N mod m has o(1) advantage on ε (and on
     any function of p mod m) for m ≤ X^{1−η}.
   - **Unconditionally** the same holds for the low j bits of N with 2^j ≤ N^c (Bourgain's prescribed binary digits).
   - So the residue of N stays uninformative far past N^{1/4}, the point where the residue *of p*, given alongside N,
     already factors N (BDF Cor. 1). These are different experiments; there is no common threshold.
   - E1 (exact census, n = 28, 32, 36): the proved bound holds in all 21 informative cells (72 cells are vacuous).
     With an exact unordered-pair random-label null (post-hoc), A*/null ≤ 1.003 for j ≥ 9 (b = 14, 16) and j ≥ 10
     (b = 18). Below that, A* equals the prime set's global ε imbalance, which the proved bound covers.
   - Finite automata and AC⁰ remain OPEN. There is a restricted extension via Müllner's type-II bound (codex; hypotheses
     not re-verified by us) and an LMN sufficient criterion.
6. **No positive leak that does not split N.** One Fermat step factors an Ω(N^{−1/4}) fraction of instances (codex).
   That gives an exponentially small advantage *by splitting*.

**Most decisive open question (refined).** Is there a poly-time moved-divisor map outside SR1/SR2?
- Prototype: N ↦ #E(Z/N) for one fixed **ordinary (non-CM) curve E**, on a promise domain of good reduction. Its moved
  factor #E(F_p) varies with p but obeys no polynomial identity.
- CM curves can collapse. For example y² = x³ − x on p, q ≡ 3 (mod 4) has a_p = 0, so the count is (p+1)(q+1), which
  factors N (codex_2).
- The recovery step is reportedly known. Dieulefait–Urroz (arXiv:1911.11004 §3; preprint, not re-checked by us) factor
  N deterministically from the *factorisation* of N′ = #E(Z/N) for a single curve, via Coppersmith on
  #E(F_p) = p+1−a_p with |a_p| ≤ 2√p.
- So a poly-time one-curve count would reduce factoring N to factoring one N′ ≈ N whose prime factors move with p: a
  malleability-type reduction (Paillier–Villar 2006). The open part is computing that count without the factors.

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
| h(N), h ∈ Z[x] N-independent and shared over a δ-dense set of partners q (any height, random coefficients; degree reduced via DEG) | **PROVED barrier (SR1)** | divisors are p^a · (divisor of cont h) · (O(D/δ)-smooth) |
| Φ_m(N) | **PROVED: SR1 applies**; p ∤ Φ_m(N) since Φ_m(0) = 1 (m > 1) | bit length Θ(φ(m) n), so m ≤ poly(n) (codex_1 §2.2) |
| Norm_{K/Q}(α), α with coordinates polynomial in N, fixed K | **PROVED: SR1 applies** | the norm is a determinant, hence a polynomial in N |
| Norms with N-dependent algebraic data (e.g. NFS: f(m) ≡ 0 mod N) | **OPEN in general; observation** | p enters via N \| Norm(θ − m) (multiplier) or via accidental gcd splits; no moved divisor known |
| Resultants Res(f, g) with a common root mod N | **PROVED: N \| Res** (multiplier); cofactor: OBSERVATION only | no moved divisor is exhibited; that the cofactor is unrelated to p is not proved |
| Discriminants with p \| disc f, q ∤ disc f | **PROVED: already splits (Lemma EUC)** | Euclid on (f, f′) over Z/N meets a zero divisor because deg gcd differs mod p and mod q (tested) |
| (p+1)(q+1), (p²+p+1)(q²+q+1), (p+r)(q+r), polynomial point counts | **PROVED: ≡ FACT (SR2)**; generalises nonsmooth §5.3 | exact algebraic identity ⇒ root-finding over Z. KNOWN special cases: (p+1)(q+1) = σ(N) and φ(N) give p + q. In general σ ⇒ factoring (Bach–Miller–Shallit, SIAM J. Comput. 15(4) (1986) 1143–1154) and φ ⇒ factoring (Miller 1976 under ERH; randomized unconditionally). SR2 generalises this folklore |
| #E(Z/N) for chosen curves/twists | **KNOWN/PROVED: ≡ FACT (randomized)** | 4 twist counts sum to 4(p+1)(q+1), with all four sign patterns hit with probability 3/8 per attempt (codex_1 §2.6; = eq. (1) of Dieulefait–Urroz, arXiv:1911.11004, preprint). Kunihiro–Koyama (EUROCRYPT 1998, LNCS 1403, 47–58) use many curves and, per Dieulefait–Urroz Remark 2, an inaccurate distributional assumption. Martín–Morillo–Villar (Appl. Math. Lett. 14 (2001) 341–346): orders of points mod N ⇒ factoring. Dieulefait–Urroz claim a deterministic reduction from one curve plus one twist (not re-checked) |
| #E(Z/N), one fixed ordinary (non-CM) curve, good-reduction promise | **OPEN** | #E(F_p) = p+1−a_p is non-algebraic; SR2's approximate form fails (error a_p·q ≈ N^{3/4} ≫ s). CM curves can collapse: y² = x³ − x with p, q ≡ 3 (4) gives (p+1)(q+1) ⇒ FACT (codex_2) |
| Lucas / torus recurrences mod N | **PROVED: residues ≠ orders** | evaluating u_k mod N does not give the integer group order |
| a^N − 1 ⊇ a^p − 1; Lucas u_N ⊇ u_p; EDS | **PROVED: exponential size**; reduction mod M destroys divisibility | natural moved divisors, besides group orders, that we know of |
| N′ = g(p)h(q) a balanced semiprime with algebraic g, h | **PROVED: ≡ FACT (SR2)** | exactly the "self-reduction to the same distribution" |

### 1.3 Theorem SR1 (shared-polynomial maps cannot move a large prime). PROVED

**Setting.** h ∈ Z[x] \ {0}, D = deg h ≥ 1, c = cont(h) > 0, p prime, and I = (x₀, x₀ + y] an interval with y ≥ 2 and
x₀ ≥ √y (true for I = (p, 2p]). Let
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
  for reduced classes a and ℓ < y (Mathematika 20 (1973), Thm 2; any shifted interval). Here log(y/ℓ) > ½ ln y. The
  class 0 mod ℓ contains no prime of I, since x₀ ≥ √y > ℓ; without x₀ ≥ √y, add 1 to the count. ∎
- **Constant.** "A δ-fraction of the primes in (p, 2p)" means |Q| ≥ δκ_p·p/ln p with κ_p := (π(2p) − π(p)) ln p/p < 1
  (≈ 0.95 at p = 10⁴, → 1). So the bound reads ℓ < 1 + 4D/(κ_pδ) = 1 + (4 + o(1))D/δ. All statements also need
  h(pq) ≠ 0 on Q (drop the ≤ D exceptional q).

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
- **Application (corrected after both reviewers).** DEG does *not* bound the degree of whatever polynomial a circuit
  computes. For example h = N + ∏_{n-bit M}(N − M) equals N on every n-bit input. The correct use is this.
  - Let V(x) = ∏_{M ∈ [2^{n−1}, 2^n)}(x − M), which is monic. Replace h by r = h mod V ∈ Z[x]; division by a monic
    polynomial keeps integrality, and cont(h) | cont(r).
  - r agrees with h on every n-bit input. If the outputs have ≤ T bits, DEG (with D + 1 ≤ 2^{n−1} consecutive inputs)
    gives deg r ≤ poly(T).
  - SR1 applies to r, with Q restricted so that pq stays in the n-bit range.
  - Content primes of r may differ from those of h, but they remain independent of p.

**Corollary SR1′ (self-reduction form).**
- Let (h_r) be N-independent polynomials (they may depend on n and on coins r) computed in poly time.
- Suppose public g_r satisfy g_r(p) | h_r(pq) for a δ-fraction of balanced pairs.
- Averaging (per seed r, or over (p, r)): if the success rate a_p of p has E[a_p] ≥ δ under the actual marginal law of p,
  then Pr[a_p ≥ δ/2] ≥ δ/(2 − δ).
- For equal-bit-length pairs, discard the p whose partner range is shorter than ε·p; their pair-mass is O(ε²), so take
  ε ≍ √δ. The remaining good p have Q_p of relative density ≥ δ/2 (up to constants) among the primes of their partner
  interval, of length y ≥ ε·p. SR1(a), applied with that y, gives ℓ < 1 + O(D/δ). The trimming only guarantees
  |Q_p| > D(√y + 1) for δ ≥ 1/poly(n) and large n.
- Hence every prime factor of g_r(p) is p, a prime of the content (from a fixed, p-independent *catalogue*), or
  poly(n)-small. Codex_2 adds that an arbitrary public g_r may still *select* among content primes depending on p, so
  the content part is a fixed catalogue, not a fixed value.
- In particular no such family yields balanced semiprimes N′ = g_r(p)·h′(q) whose large prime g_r(p) moves with p
  outside the content catalogue. The catalogue has ≤ log₂ cont(h_r) = poly(n) primes.
- Fixing a seed is legitimate: the theorem is about one polynomial. Averaging only controls |Q_p|.

**Refinement: exact p-independence (PROVED for p ≥ p₀(h); E2).**
- Let u_ℓ(h) = min_{z: ℓ∤z} v_ℓ(h(z)) and U(h) = ∏_ℓ ℓ^{u_ℓ(h)}. Only primes ℓ ≤ D+1 or ℓ | c contribute: otherwise a
  nonzero polynomial of degree ≤ D < ℓ − 1 misses some unit.
- If Q = all primes in (p, 2p) and p ≥ p₀(h), then G_p := gcd_{q∈Q} h(pq) = p^{v_p(G_p)} · U(h)/p^{v_p(U(h))}.
  - "≥" holds once pq is a unit at every candidate prime ℓ, e.g. p > max(D+1, primes of c up to 2p). Counterexample
    for tiny p (codex_2): h = x² − 1, U = 24, p = 2, Q = {3} gives 35.
  - "≤": first use SR1(a) to confine the possible excess primes to a finite set. For each, the minimum valuation over
    units is attained on a class mod ℓ^{u_ℓ(h₁)+1}, and {pq mod ℓ^{u+1}} meets every unit class once Siegel–Walfisz
    applies to these finitely many fixed moduli. That is a p₀(h); there is no uniformity for growing h_n.
  - For ℓ | c with ℓ − 1 > D, |Q| > D(y/ℓ + 1) suffices.
- The smooth part of G is then literally the same integer for every p: the factors never move.
- KNOWN context: U(h) is Pólya's fixed divisor restricted to units; for primitive h of degree D the fixed divisor
  divides D!. Lemma DEG's Δ^D h = D!·lead(h) is classical. Neither is new here.
- Codex's caveat stands: for a *subset* Q (density δ < 1), the small smooth part can carry residue information about p.
  For example, h = x − 1 with Q = {q ≡ 1 mod 4} detects p mod 4. The only proved bound on this data is SR1(b)'s bit
  length; it is small-prime data, not a moved large factor. (v1's "O(log(1/δ)·log(D/δ)) bits" had no proof and is
  withdrawn.)

**E2 (pre-registered, exact).**
- 40 polynomials: D₀ ∈ {1, 2, 3, 5, 8}; coefficients up to 2^64; content up to 3 primes < 2^20. The 16 polynomials
  with D₀ ≥ 2 and odd k are planted with x·∏_{i≤D₀}(x − i), so total degree ≤ 17. (v1 said "half"; D₀ = 1 is never
  planted.)
- 80 primes p (40 after 10⁴, 40 after 10⁵), with Q = all primes in (p, 2p). That gives 3200 cases.
- **P1 (proved bound):** held in 3200/3200 cases; 2880 satisfy the hypothesis |Q| > D(√p + 2) and the rest are vacuous.
  The worst ratio ℓ/L was 0.39.
- **P2:** G_p/p^{v_p} = U(h) in 3200/3200 cases. Each h produced exactly one value across all 80 primes p.
- **Rerun (D5).** `selfred_e2_rerun.jsonl` uses an exact cap-free U(h) (gcd of h(a + ℓt), 1 ≤ a < ℓ, 0 ≤ t ≤ D;
  codex_2). P1 is now checked on the support of S = (G/gcd(G, c))_{p′}, so excess content-prime valuations are covered,
  with a certified integer log bound. Results are identical: P1 3200/3200 (certified), P2 3200/3200, every U and G_p
  unchanged, largest prime of S = 7.
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
- In case (ii), integer roots are found by exact real-root isolation (Sturm) or by LLL factorisation over Q (LLL 1982).
  Each root r is tested for r | N. ∎
- **Implementation.** v2: `rigidity.integer_roots` uses Sturm-sequence bisection on the squarefree part, with exact
  integer/rational arithmetic; cost poly(deg) · log(hi − lo) evaluations.
  - v1's derivative-bracketing routine was **wrong**. It missed roots when P′ had two roots in one unit interval:
    P = 2(x+1)(7x² + 3x − 1) (Opus critic) and x(x+1)(2x+1) (codex_2). `solve_rigid` failed on valid promises (2 of 40
    random 82-bit cases; codex_2's N = 35).
  - Both counterexamples, a 3000-polynomial brute-force fuzz and a 200-case adversarial SR2 stress test (0 failures)
    are now in the tests. The theorem was never affected.

**Extensions and exceptions** (codex_1 §2.4; re-checked):
- **Rational functions.** The degenerate case is "A(x)B(N/x) ≡ const after cancellation", which is broader than the
  monomial case. Example: A = x + 1, B(y) = y/(N + y). The general statement applies to the whole rational expression
  R_N(X) = U/V. If R_N is nonconstant, F = U − N′V is a nonzero polynomial with root p.
- **Cofactors.** An unknown cofactor K ≠ 0 with N′ = K·A(p)B(q) is handled if K ranges over a *publicly enumerable*
  set of ≤ poly(n) values. The same holds for an additive error. Denominators must be nonzero at the actual factors.
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
- So barriers currently within reach are restricted-model (SR1, SR2, PT) or relativized. SR3 does not show that
  unrestricted lower bounds are unprovable, only that they would be major results.

### 1.6 Random self-reducibility and instance hiding: what is known (KNOWN; details codex_1 §2.7)
- **Abadi–Feigenbaum–Kilian**, JCSS 39 (1989) 21–50. §6 Question 1 asks, separately, whether factoring is encryptable
  (leaking at most the input length) and whether it is random-self-reducible. Their negative results concern NP-hard
  functions.
- **Feigenbaum–Fortnow**, SICOMP 22(5) (1993) 994–1005.
  - Thm 3.1: if S ∈ NP has a nonadaptive random self-reduction, then the *complement* S̄ ∈ AM_poly. (codex_2 caught the
    complement bar, which the extracted PDF text loses; codex_1 had dropped it too.)
  - Cor. 3.3: a third-level PH collapse for NP-complete sets.
  - Factoring's decision version ("N has a prime factor ≤ b") is in NP ∩ coNP (Pratt certificates), so it is not known
    to be NP-complete, and FF gives no obstruction for factoring.
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
- **Closest known notion: malleability of factoring.**
  - Paillier–Villar (ASIACRYPT 2006, LNCS 4284, 252–266) ask whether an oracle factoring moduli N′ ≠ N (coprime to N)
    helps to factor N. They conjecture that typical key generators are non-malleable.
  - Dieulefait–Jiménez Urroz (*Small primitive roots and malleability of RSA moduli*, J. Comb. Number Theory 2(2) (2010)
    171–179; arXiv:0801.0030) refute a form of the conjecture. Their N′ is explicit but exponential-size, built from
    m^N − 1 in base m (cf. our a^N − 1 row).
  - This is a modulus-changing reduction *to a factoring oracle*, not a worst-to-average self-reduction. The "not
    located" statement above stands only in that narrower sense.
- **Accurate statement:** factoring is not known (to us) to be random self-reducible over the modulus. AFK list this as
  open.

---

## 2. The open partial-information bits (question 2)

### 2.1 Status table (oracles on their full stated domain; "≤" = poly-time Turing reduction)

| Quantity | Value on N = pq | Status | Known reductions / references |
|---|---|---|---|
| ω(M) | 2 | **OPEN** | ω ≤ d (signature interpolation, ≤ log₂M calls; Shallit–Shamir IPL 20 (1985) 151–153). C9/O9 (FOCS OCR; UNVERIFIED in ANTS) |
| μ(M) | 1 | **OPEN** | μ ≤ d with **one** call: D = d(M^{r−1}), r > log₂M + 1 prime, is a pure r-power iff M is squarefree (Shallit–Shamir; proof codex_1 §3.3). μ = λ·SF |
| Liouville λ(M) | 1 | **OPEN** | λ ≤ d (signature); not the factoring-equivalent Carmichael λ |
| squarefreeness SF(M) | yes | **OPEN** (C8/O8: FOCS OCR; UNVERIFIED in ANTS) | SF ≤ φ deterministic (Woll, Inf. Comput. 72 (1987), per abstract); SF ≡ maximality testing of orders (Buchmann–Lenstra JTNB 6 (1994) Thm 6.12). Not in AC⁰[ℓ] for any prime ℓ (Allender–Saks–Shparlinski JCSS 62 (2001) Cor. 2); SF ∉ AC⁰ (Bernasconi–Damm–Shparlinski 1998; Bernasconi–Shparlinski, STACS 1999, LNCS 1563, 47–56). SF ≤ the count of automorphic representations of one weight on Γ₀(N), and factoring ≤ newform counts, probabilistically (Gu–Martin, Canad. Math. Bull. 62(1) (2019) 81–97). Booker–Hiary–Keating (Duke 164 (2015)), under GRH for quadratic Dirichlet L-functions: Prop. 2.3 returns a square factor or a conditional squarefreeness certificate, with running time exp((log M)^{θ*+o(1)}) stated **for squarefree inputs** (θ* = 2/3 is HEURISTIC). Thm 1.1 proves RSA-210/220/230/232 *not squarefull* (not squarefree). Rigorous elementary bound M^{1/6+o(1)} |
| rad(M) | N | **OPEN** | rad ≡ squarefree decomposition ≡ parity kernel s(M) (M = u²s), deterministic, ≤ ⌊log₂M⌋ calls (Prop RAD, PROVED; `squarefree_decomposition_via_rad`, tested). rad ≡ ring of integers (Buchmann–Lenstra Thm 1.3; Chistov 1989). rad ≤ φ, ≤ Carmichael λ (Landau Inf. Comput. 78 (1988)). Adleman–McCurley, ANTS-I, LNCS 877 (1994) 291–322, problem C7 (M = u²s, the parity kernel; equivalent to rad by Prop RAD, not literally rad) with O7a (algorithm?) and O7b (randomized reduction from factoring?), corroborated via the CRYPTREC ESIGN report §3.1. The same list appears in Adleman FOCS 1994 (read via OCR by codex only) |
| p²q special case | rad = pq, s = q | these *values* factor p²q; computing them OPEN | squarefreeness on {p²q} ∪ {pqr} OPEN; nothing polynomial known for r = 2 |
| p mod 4 (smaller factor) | — | **OPEN** | = ε on N ≡ 1 (4); orientation bit on N ≡ 3 (4) |
| ε = (−1/p), N ≡ 1 (4) | ±1 | **OPEN** | ε = +1 ⇔ −1 is a QR mod N ⇔ φ(N) ≡ 0 mod 8 ⇔ N is a sum of two squares; Seurin PKC 2014 (2-Φ/4-hiding). h(−4N) mod 8 does **not** determine ε (Rédei; 4 counterexamples, codex_1 §3.7) |
| φ(N) mod ℓ, ℓ fixed odd prime | determines s mod ℓ | **OPEN**; decides ℓ* QR mod N when (N/ℓ) = 1 (nonsmooth RANK) | Φ-hiding [ℓ \| φ(N)]: Cachin–Micali–Stadler 1999, an assumption |
| φ(N) mod ℓ, ℓ an **input** prime | — | **≡ FACT (PROVED)** | balanced promise: CRT over odd primes with ∏ℓ > 2^{⌈n/2⌉+2} gives s exactly (general N: bound N + 1). Calls = π(B) − 1 with θ(B) ≈ (n/2) ln 2, ~ 0.347n/ln n asymptotically; actual 12, 26, 75, 131 at n = 80, 256, 1024, 2048 (`factor_via_phi_mod_ell`, tested; deterministic trial-division primality for the small ℓ) |
| φ(N) mod 2, 4 | 0, 0 | trivial | — |
| fixed-base QR_a(N) | — | **OPEN** for nontrivial a | a function of p mod 4\|a\| on the fiber; random-base QR = Goldwasser–Micali assumption |
| d(M) | 4 | OPEN for general M | Shallit–Shamir: d is equivalent to the exponent signature |

**Two corrections to earlier project notes.**
- Allender–Saks–Shparlinski do cover squarefreeness: TC⁰ ⊆ AC⁰[Square-Free] (Thm 1) and SF ∉ AC⁰[ℓ] (Cor. 2). The
  factorbits wording "primality ∉ AC⁰[p]" undersold it.
- The Boneh–Durfee–Frankel residue result is **Corollary 1** (p. 28) of the ASIACRYPT 1998 paper, not "Cor. 2.2"
  (codex_1 §4.4; checked by codex in a proceedings scan, not independently re-verified). The factorbits row should be
  read accordingly.
- Do not cite Bach–Shallit "Chapter 14": Volume I has Chapters 1–9, and Chapter 14 belongs to the unpublished Volume
  II (literature critic; codex_1 §3.5 was wrong). The note does not rely on it.

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
  and SR1 (shared-polynomial queries cannot move a large prime). None is an oracle separation of the real integer
  problems.
- KNOWN restricted-model results:
  - In the generic ring model, breaking RSA ≡ factoring (Aggarwal–Maurer, EUROCRYPT 2009, LNCS 5479).
  - Computing the Jacobi symbol and deciding QR are each ≡ factoring generically (Jager–Schwenk, ASIACRYPT 2009;
    J. Cryptology 26(2) (2013) 225–245), although the Jacobi symbol is easy in the standard model. So generic-ring
    hardness (cf. PT) is weak evidence.
  - Boneh–Venkatesan (EUROCRYPT 1998, LNCS 1403, 59–71) separate low-exponent RSA from factoring for algebraic
    (straight-line) reductions: the template for a restricted separation.

---

## 3. Theorem PT's open classes (question 3)

### 3.1 Positional floor: functions of N mod m

**Theorem L (codex_1 §4.2; PROVED under GRH; proof re-checked).**
- Setting: P, Q iid uniform primes in [X, 2X], P ≠ Q; for ε, condition on PQ ≡ 1 (mod 4); L = lcm(4, m) < 4X.
- Under GRH mod L: TV(law of P mod L, uniform on units) ≪ √(φ(L)/X)·log X·log(XL).
- Consequently every predictor F(N mod m), of any complexity and possibly randomized, has advantage on ε at most
  O(√(m/X) log²X) + O(log X/X).
- For any target f(p mod m), relative to its best constant, the advantage is **o(1)** (not the same rate). Ordering is
  handled by K sub-interval bins at cost O(K√(m/X) log²X + 1/K). For ε no bins are needed, since χ₄(P) = χ₄(Q) on E.
- Proof:
  - GRH gives |μ̂(χ)| ≪ log X log(XL)/√X for each nonprincipal χ (Montgomery–Vaughan, Multiplicative Number Theory I,
    Thm 13.7 / (13.21)).
  - Parseval gives Σ_a |μ(a) − 1/φ|² = φ^{−1}Σ_{χ≠χ₀}|μ̂(χ)|². Cauchy–Schwarz turns this into TV.
  - In the uniform model the pair (U, V) given UV = c is uniform on the fiber, and χ₄(U) is a fair coin independent of
    c (trackE F1).
- **Range:** m ≤ X^{1−η} ≍ N^{(1−η)/2}, i.e. any fixed power below N^{1/2}. Restrictions: the independent-prime law,
  a modulus fixed (or chosen independently of the instance), observation of N mod m only, and a single modulus (for
  several moduli, apply it to their lcm).
- This supersedes our first, pointwise argument, which reached only m ≤ N^{1/4}/poly(n). That argument produced a
  spurious "meeting" with the Coppersmith threshold.

**Unconditional low bits (codex_1 §4.3; PROVED from a KNOWN theorem).**
- Bourgain (Israel J. Math. 206 (2015) 165–182, §1 Theorem, (1.1)–(1.2)): primes with r < γb prescribed binary digits
  have the asymptotic density 2^{−r}.
- Prescribing the top bit and the j − 1 nonconstant low bits uses r = j positions; codex_2 checked that position b − 1
  is allowed. Bourgain's asymptotic is uniform over all assignments, which the TV sum over 2^{j−1} classes requires.
  Hence the law of p mod 2^j on a dyadic interval is o(1)-close to uniform for j ≤ γ′·n/2.
- So any function of the low j bits of N has advantage o(1) on ε for 2^j ≤ N^c, with an unspecified absolute c > 0.
- Barban–Linnik–Tshudakov (Acta Arith. 9 (1964), range x^{3/8}) and Gallagher (Invent. Math. 16 (1972), x^{2/5}) are
  stated for fixed **odd** prime powers in the checked secondary source (Baier–Pujahari, Acta Arith. 204 (2022)). We do
  not use them for 2^j. (Our blind codex prompt `codex/selfred_1.md` had suggested BLT/Gallagher for 2^j; the
  pre-registration did not cite them.)

**No secrecy threshold at N^{1/4} (PROVED distinction).**
- p mod M with log₂M ≥ n/4 + O(log n) factors N together with N (BDF Cor. 1; trackE E1).
- N mod m (alone, fixed modulus, independent-prime law) reveals essentially nothing about p mod m up to any fixed
  power below N^{1/2} (GRH).
- These are different experiments, because the fiber point is not given. The N^{1/4} threshold belongs to the hint,
  not to the residue of N.

**E1 (pre-registered, exact census).**
- Population: ordered iid pairs from S_b = all primes in [3·2^{b−2}, 2^b), b = 14, 16, 18 (n = 28, 32, 36), with
  |S_b| = 431, 1491, 5296.
- For every j ≤ n we compute the exact Bayes-optimal advantage A*_j of any function of N mod 2^j on ε (restricted to
  N ≡ 1 mod 4), the exact TV₁ and TV₂, and the proved bound 3TV₂/(1 − 2TV₂).

| b | j | A*_j | exact unordered null (D2′) | A*/null | TV₁ | proved bound |
|---|---|---|---|---|---|---|
| 16 | 4 | 6.7e−4 | 1.07e−3 | 0.63 | 0.012 | 0.054 |
| 16 | 8 | 1.8e−3 | 4.28e−3 | 0.42 | 0.077 | 0.415 |
| 16 | 12 | 1.13e−2 | 1.71e−2 | 0.66 | 0.417 | vacuous |
| 16 | 16 | 6.5e−2 | 6.8e−2 | 0.96 | 0.954 | vacuous |
| 16 | 20 | 0.255 | 0.256 | 1.00 | 0.997 | vacuous |
| 16 | 24 | 0.470 | 0.470 | 1.00 | 1.000 | vacuous |
| 18 | 4 | 1.9e−3 | 3.0e−4 | 6.27 | 0.004 | 0.016 |
| 18 | 8 | 1.9e−3 | 1.21e−3 | 1.57 | 0.038 | 0.185 |
| 18 | 12 | 3.3e−3 | 4.82e−3 | 0.67 | 0.198 | vacuous (1.89) |
| 18 | 19 | 5.3e−2 | 5.4e−2 | 0.98 | 0.980 | vacuous |

- **P1 (proved bound): holds in all 21 informative cells** (bound < 1/2: 5, 7, 9 cells for b = 14, 16, 18).
  - The remaining 72 cells are vacuous: 63 have TV₂ ≥ 1/2, so no bound is defined, and 9 have a bound ≥ 1/2. The code
    scores these P1_ok = True trivially. v1 said "93/93"; the opus critic and codex_2 flagged it.
  - Proof of the bound: in the uniform product model Q(E) = 1/2 and Q(F = ε ∧ E) = 1/4. Hence P(F = ε ∧ E) ≤ 1/4 + d
    and P(E) ≥ 1/2 − d with d = TV₂, so A* ≤ 3d/(1 − 2d). TV₂ ≤ 2TV₁ by the triangle inequality; the TV₁ form needs
    TV₁ < 1/4.
  - Where it is informative, the empirical TV₁ scales like √(2^j/|S|): the random-model shape of Theorem L's √(φ/X).
- **P2: REFUTED as registered (a prediction error, not a leak).**
  - (a) "A* ≤ 0.02 for j ≤ n − 2log₂b − 2" ignored finite-population fitting. When classes hold few pairs, even random
    labels give a large Bayes advantage; the registered null itself exceeds 0.2 there.
  - (b) "Ratio to null ∈ [0.5, 2]" failed for b = 14, j = 2–4 (5.3, 3.8, 2.7).
  - At small j, A* is the population's global ε imbalance (#primes ≡ 1 vs 3 mod 4 in S_b), constant in j and within
    P1. At b = 18 the same effect gives registered ratios 17.7 … 2.2 for j = 2–8 (A* = 0.00189 throughout, bound
    0.0028 at j = 2). These were outside P2's scope only because A_null < 10⁻³.
  - The registered null also double-counted (X, Y) and (Y, X), which share N and ε.
- **Post-hoc D2′ (labelled; replaces v1's approximate D2).** One fair coin per unordered pair, with weights 2
  (off-diagonal) and 1 (diagonal) per class, exact: ½E|2S_o + S_d| (Fraction for o + d ≤ 3000, FFT convolution above,
  cross-checked). File: `selfred_e1_posthoc_v2.jsonl`; b = 18 included.
  - A*/null ≤ 1.003 for all j ≥ 9 (b = 14, 16) and j ≥ 10 (b = 18). The maximum is 1.0034 at b = 14, j = 19.
  - Below those j the excess is the global imbalance described above (b = 14: 3.8, 2.7, 1.9, 1.3 at j = 2–5; b = 18:
    12.5 … 1.1 at j = 2–9).
  - Matching a random-label expectation is descriptive, not a calibrated test.
- **Descriptive:** j_{0.05} = 13, 16, 19 for n = 28, 32, 36 (≈ n/2 − 1).
  - This is where classes become small enough for *finite-population lookup*. N carries full information
    (A*_n = 1/2), so the growth reflects the information becoming usable by table lookup over a finite set, not a
    computationally useful predictor.
  - Nothing in E1 is evidence of ε leakage from low bits beyond the imbalance of the finite prime set.

### 3.2 Automata, AC⁰, integer-structure algorithms

| Class | Status |
|---|---|
| multiplicative automatic sequences | floor (KLM 2022 Thm 1: on large semiprimes a(N) = f₂(N), a function of N mod a fixed period) — **PROVED**, codex_1 §4.5 |
| sequences with the carry and Fourier properties (Müllner Defs 4.1–4.2, generalising Mauduit–Rivat, *Prime numbers along Rudin–Shapiro sequences*, J. Eur. Math. Soc. 17(10) (2015) 2595–2642), linear Fourier decay | **codex claim**: Müllner's type-II bound (Duke 166 (2017), arXiv v3 Prop. 5.5, (29)) with a_u = 1_prime·χ₄, b_v = 1_prime gives o(1) correlation with ε on balanced semiprimes. **UNVERIFIED by us** (hypotheses not re-checked); not all automata satisfy them |
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
  - `rigidity.py` and `reductions.py` import neither sympy, flint, harness nor census (AST test). codex_2 also
    inspected them by hand and found no factor access; the AST test alone would not establish this.
  - Their oracles are passed in from `harness.py`.
- **No inverse-polynomial non-splitting predictor** of any open bit was found.

---

## 5. Candidate cards (rule 9)

| Mechanism | Model | Outcome |
|---|---|---|
| polynomial / norm / cyclotomic self-reduction | one N-independent h ∈ Z[x] shared over a δ-dense set of partners q; any height; values reduced mod the input-range vanishing polynomial (DEG) | **refuted for this model** (SR1): no moved large prime outside a fixed content catalogue; the smooth part is p-independent for full Q (E2) |
| algebraic self-reduction to semiprimes, torus/point-count orders | exact public polynomial identity N′ = A(p)B(q), poly(n) degree/height | **refuted as an amplifier**: it is already a factoring algorithm (SR2) |
| discriminant-based moves | monic f with p \| disc f, q ∤ disc f | **refuted**: Euclid splits (EUC) |
| signature oracles on multiplicative queries | products/powers of N and factored integers | **refuted for this model** (SIG) |
| residues / low bits of N as ε predictors | any function of N mod m, fixed m, independent-prime law | **refuted** for m ≤ X^{1−η} (GRH) and 2^j ≤ N^c (unconditional) |
| #E(Z/N) for one ordinary curve; non-algebraic or hidden-cofactor moves; N-dependent divisibility constructions | — | **not yet tested / OPEN** |

## 6. Most decisive open question

The algebraic versions of factorbits' question are settled in two precise models:
- shared-polynomial divisibility is impossible beyond small primes (SR1);
- exact public algebraic identities are already factoring algorithms (SR2).

What remains open includes, but is not claimed to be limited to:
- (a) **Non-algebraic moved divisors.** The prototype is N ↦ #E(Z/N) for one fixed ordinary (non-CM) curve E on a
  good-reduction promise. Its moved factor #E(F_p) = p+1−a_p varies with p but satisfies no polynomial identity.
  - Chosen twists make it ≡ FACT. CM curves can collapse to (p+1)(q+1).
  - Is the one-curve count factoring-equivalent? Recovery is reportedly known: Dieulefait–Urroz (arXiv:1911.11004 §3,
    preprint, not re-checked) factor N from the factorisation of #E(Z/N) by Coppersmith. So an efficient count would
    give a malleability-type reduction of factoring N to factoring one N′ ≈ N. The open part is the count itself.
  - Whether such an N′ supports bit-*amplification* (a weak predictor of a fixed bit of the smaller factor of N′
    transferring to a varied function of p) is a further question.
- (b) **Moved divisors inside large unknown cofactors** of polynomial size (divisibility-sequence-like).
- (c) **N-dependent divisibility constructions and adaptive queries**, which neither SR1 (one shared polynomial) nor
  SR2 (exact identity) covers.

A barrier for (a)–(c) would need a model of "non-algebraic" maps. We have none.

## 7. Reproduction
```
.venv/bin/python -m pytest -q tests/test_selfred.py                  # 33 tests, ~2 s
.venv/bin/python -m src.selfred.run_experiments e1 14 16 18           # E1 -> experiments/raw/selfred_e1.jsonl (~90 s; peak RSS 3.5 GB, D3)
.venv/bin/python -m src.selfred.run_experiments e1posthoc 14 16 18    # D2' -> selfred_e1_posthoc_v2.jsonl (~130 s; peak RSS 4.1 GB, D3)
.venv/bin/python -m src.selfred.run_experiments e2                    # E2 (v2 code) -> selfred_e2_rerun.jsonl (11 s)
```
- `selfred_e2.jsonl` (v1 code, original pre-registered run) and `selfred_e1_posthoc.jsonl` (v1 approximate D2, b = 14,
  16) are kept for the record.
- The current code reproduces the E2 file's values exactly (D5), but not the v1 posthoc file, whose approximate null
  was replaced.

**Deviations (disclosed).**
- **D1.** The registered E1 null uses de Moivre's closed form exactly for class sizes k ≤ 4000 and log-gamma double
  precision above. It is a reference only; P1 is exact.
- **D2 / D2′.** A post-hoc unordered-pair null, because the registered null double-counts orderings.
  - v1 (D2) scaled an unweighted null by a global mean weight, which codex_2 showed is not the weighted null; its ratios
    reached 1.0038.
  - v2 (D2′) is the exact per-class weighted null, with b = 18 added.
- **D3.** Memory guideline exceeded twice, each time running alone with ≥ 74% of memory free beforehand: E1 at b = 18
  (3.5 GB, 85 s) and D2′ at b = 18 (4.1 GB, 133 s). Chunking would avoid this; not redone.
- **D4.** BLT/Gallagher are stated for odd prime powers, so the unconditional 2-power statement uses Bourgain. The BLT
  suggestion came from our codex prompt, not from the pre-registration (v1 misattributed it).
- **D5.** After review, the E2 code changed:
  - exact cap-free U(h);
  - the support check on S = (G/gcd(G, c))_{p′};
  - a certified integer log bound.
  E2 was rerun (`selfred_e2_rerun.jsonl`) with identical values and verdicts.
- **D6.** The OPERATIONAL root finder was replaced (§1.4). The v1 test suite had not exercised the failing case.

## 8. Codex adversarial round (codex/selfred_2.{md,out}: MAJOR REVISION), applied in v2

| Item | Response |
|---|---|
| `integer_roots` wrong (x(x+1)(2x+1); solve_rigid(35, 6, …) = None) | replaced by Sturm bisection; counterexamples, fuzz and stress tests added (D6) |
| DEG inference false (degree-2^{n−1} polynomial equal to 1 on inputs) | corrected: apply DEG and SR1 to the remainder mod the input-range vanishing polynomial |
| SR1 headline dropped |Q| > D(√p+1) (counterexample h = L + ∏(x − pq)) | hypotheses restored in §0 |
| SR1′: averaging under the actual law of p, boundary trimming, per-seed quantifiers; content primes are a catalogue | rewritten |
| p-independence: finite-support step first; "≥" needs units; information bound unproved | fixed; bound withdrawn |
| zero class in BT; x₀ ≥ 0 | x₀ ≥ √y assumed (zero class then prime-free) |
| SR2: enumerable cofactor set, K ≠ 0, nonzero denominators | added |
| SR3: "only restricted barriers possible" does not follow | reworded to "within current reach" |
| AFK "i.e."; FF Thm 3.1 is about the complement S̄ ∈ AM_poly; "NP-completeness not known" | fixed (codex_1 had the same slip) |
| BHK runtime only for squarefree inputs; Adleman C7 is the parity kernel | fixed |
| φ mod ℓ: balanced promise, 3027 example, deterministic small-prime test, asymptotic count only | fixed in code and table |
| Theorem L: bins give o(1), not the same rate; X^{1−η} ≍ N^{(1−η)/2}; state restrictions | fixed |
| E1: "93/93" vacuous cells; "non-vacuous iff TV₂ < 1/2" wrong; TV₁ form needs TV₁ < 1/4 | 21 informative / 72 vacuous |
| D2 not the weighted null; "at or below" false (1.0038) | exact D2′; ratios reported (max 1.0034 beyond the imbalance region) |
| census: U(h) caps (2^20 + 2^20 x raised), content-prime support not checked, float log, int64 guard | all fixed; E2 rerun identical (D5) |
| "onset of fitting, not information" | reworded (finite-population lookup) |
| §0/§6 overclaims: "every algebraic way", "only two escapes", one-curve OPEN without CM exceptions, map ≠ self-reduction | all narrowed |
| CORRECT per codex_2 | SR1(a) cases incl. MV-BT on shifted intervals; Theorem H; SR2 mathematics and rational exception; RSA contrast; SIG; EUC; RAD; μ one-call; Rabin per-modulus; Theorem L's Parseval step; Bourgain corollary; KLM; Müllner row honestly labelled; Fermat leak; E1 inequality; E1 census and TV₂ code; P2 REFUTED verdict |

## 9. Opus critics

**Proof critic: ACCEPT WITH FIXES.**
- It independently recomputed E1 b = 14 cells exactly (A*, TV₁, TV₂ match as Fractions), 12 E2 cases, the E1
  inequality and the de Moivre null. It found the same root-finder bug independently (P = 2(x+1)(7x²+3x−1); 2/40 SR2
  failures).
- All 12 items applied:
  - root finder;
  - DEG corollary (remainder mod ∏(x − M));
  - SR1 hypotheses (zero class, κ_p, large-p condition, the S decomposition);
  - SR1′ averaging;
  - refinement caveats;
  - Theorem L rate;
  - E1 21/72 counts;
  - the post-hoc overclaim and b = 18 small-j ratios;
  - the D4 attribution;
  - the φ mod ℓ call counts;
  - the E2 description (16 planted, content primes, ratio only);
  - the resultant row relabelled OBSERVATION.

**Literature critic: OK WITH FIXES.**
- It found no wrong theorem number or venue among the references the note cites.
- Verified from primary text: FF Thm 3.1 (complement) and Cor 3.3; AFK §6 Q1 (verbatim); ASS Thm 1 / Cor 2;
  Buchmann–Lenstra Thms 1.3 and 6.12; Shallit–Shamir; BHK Thm 1.1 / Prop 2.3 (runtime only for squarefree inputs; 2/3
  conjectural); MV Thm 13.7 / (13.21) / Cor 13.8; Bourgain §1 (position b − 1 allowed); Baier–Pujahari; KLM Thm 1;
  Müllner Prop 5.5 / (29); Adleman–McCurley C7 (via CRYPTREC).
- Unverifiable: the FOCS C8/C9 numbering (paywall/OCR) and BDF Cor. 1 (host down; codex's scan only).
- All 10 fixes applied:
  - malleability literature (Paillier–Villar; Dieulefait–Jiménez Urroz) in §1.6;
  - Dieulefait–Urroz one-curve recovery in §0/§6, which changes the open question to "compute the count";
  - Kunihiro–Koyama qualified, plus Martín–Morillo–Villar;
  - Bach–Miller–Shallit and Miller as SR2's known special cases;
  - generic-ring results (Aggarwal–Maurer, Jager–Schwenk) and Boneh–Venkatesan in §2.4;
  - Gu–Martin and Bernasconi–Damm–Shparlinski in the SF row;
  - the Müllner framework attributed to Mauduit–Rivat JEMS 2015, not Annals 2010;
  - C7 via ANTS/CRYPTREC, with C8/C9 marked FOCS-OCR only;
  - fixed-divisor context;
  - BHK "quadratic Dirichlet L-functions"; BDF marked as codex-checked.
- It found no known modulus-changing RSR, no reduction from factoring to SF/μ/ω/ε, no oracle separation, no
  semiprime-specific automatic-sequence result, and no prior statement of Theorem L (routine; credited to codex).
