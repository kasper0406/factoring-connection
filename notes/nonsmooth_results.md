# nonsmooth: factoring paradigms not based on smooth numbers (phase 16, route 3)

Agent `nonsmooth`, 2026-10-20.
- Pre-registration: `notes/nonsmooth_prereg.md` (commit `3b1b19e`, before E1/E2 ran).
- Code: `src/nonsmooth/{equiv,forms,cheng,census,run_experiments}.py`. Tests:
  `tests/test_nonsmooth.py` (67 pass; 46 at pre-registration, more added since, see §7).
- Raw data: `experiments/raw/nonsmooth_e{1,2}.jsonl`. Summary:
  `experiments/analysis/nonsmooth_summary.json` (commit `4ea0295`).
- GPT-6.1-sol: blind round `codex/nonsmooth_1.{md,out}`; adversarial referee
  `codex/nonsmooth_2.{md,out}`. Its findings are applied throughout, listed in §8.
- Two read-only Opus critics: correctness/cost and literature. Both returned ACCEPT WITH FIXES;
  all fixes are applied (§9).
- Labels: **KNOWN** (with reference); **PROVED** (proof here; all routine or folklore, "no prior
  formulation found by us" at most); **CONDITIONAL**; **HEURISTIC**; **experimentally supported**;
  **SPECULATIVE**; **REFUTED for this model**; **UNVERIFIED**.

## 0. Verdict (plain language)

**Among the paradigms examined, no smoothness-free route to general factoring survives below
L[1/3], or even below exponential time.** This is a verdict about known methods and the mechanisms
we could construct. It is not a lower bound.

The organising question was: which factoring-equivalent invariant of N might have an
unexpectedly cheap computation that does not factor N first? We mapped about 25 such invariants.
They fall into three classes:
- **(a) Easy at primes, products of local factors at composite N.** φ, λ, σ, ψ = #P¹(Z/N), r₄,
  Pillai's Σ gcd, #E(Z/N), tori, τ(N), a_f(N), modular-form dimensions at level N, and the Hurwitz
  aggregate C(N).
  - The fast prime methods we examined (Schoof, Edixhoven–Couveignes, p-adic methods) all use the
    Frobenius at p as an operator on an étale algebra. This is a taxonomy, not a theorem.
  - Prop FROB (§5.2, proved) shows that producing that operator over Z/N factors N. It covers
    Schoof's own algebra (Z/N)[x]/ψ_ℓ whenever (N/ℓ) = −1.
  - Every global formula we found for these invariants is a sum of length ≥ N^c whose evaluation
    amounts to divisor detection.
- **(b) Hard even at primes.** h(−p), regulators, ord_p(a). A route through these must solve a
  problem not known to be easier than factoring.
- **(c) Computable by reciprocity.** Jacobi symbols, power-residue symbols at (N), and Artin
  symbols of (N) in abelian fields. These depend only on N mod the conductor (§5.4), so they carry
  0 bits beyond it (Track E F1).

**The non-smooth methods that exist are exponential, special-purpose or quantum.**
- Exponential: Fermat/Lehman, Pollard–Strassen/BGS, SQUFOF, Rubinstein and Hiary (N^{1/3}),
  Coppersmith (N^{1/4}), Hittmeir (N^{2/9}), Harvey (N^{1/5}, Math. Comp. 2021) and
  Harvey–Hittmeir's log-log refinement (2022).
- Special-purpose: CM/anomalous curves (Cheng's 4p−1 method).
- Quantum: Shor.
- The subexponential infrastructure variants (SQUFOF2; Murru–Salvatori's regulator route) are
  smoothness-based.

**The most promising candidate was codex's survivor. It is reduced, not refuted (§4.1).**
- The idea: approximate the Hurwitz class-number aggregate C(N) = Σ_t H(4N − t²) = 2N + 2q to
  within N^{1/4} cheaply, then apply Coppersmith.
- Since C(N) = 2N + 2q exactly, such an approximation *is* a high-bits hint for q. The candidate
  restates the hint problem; it adds no mechanism.
- For one natural analytic expansion this is exact. In the Ramanujan-sum series for σ(N), every
  term with modulus c < p vanishes identically (proved, tested). Truncating returns N and misses
  p + q + 1 entirely.
- That the other representations fail the same way (Rademacher/Kloosterman, mock-modular, Shintani)
  is HEURISTIC.

**Order-annihilation methods (p−1, p+1, ECM, tori, abelian varieties).**
- Prop OA (§5.1, proved, for scalars m = N^k·S) says either the part of #G(F_p) coprime to pq is
  smooth, or G(F_p) is a {p, q}-group. The second case is an anomalous variety, or a "q | #G(F_p)"
  coincidence.
- The anomalous branch is real. Our operational Cheng-type factorer splits 180/180 special N up to
  512 bits and 0/180 random controls (E1).
- But it is thin. Coverage decays like 0.25·√(Y/p) when the CM discriminant is restricted to
  D ≤ Y. The pre-registered slope test measured −0.127 against a prediction of −0.125 (E2).

**Minimal primitive (§5.5).** For the *linear* divisor formulas (φ, σ, ψ, r₄), knowing the value
modulo one small odd prime ℓ already decides a fixed-base quadratic-residuosity bit: "is ℓ* a
square mod N", e.g. "is N = x² + 3y²". That is the Track E ε-type bit (Seurin's 2-Φ/4-hiding for
base −1). No algorithm for it faster than factoring is known, and no reduction from factoring to it
is known.

## 1. Equivalence map

Column key:
- **→F**: what the value, or an oracle for it, gives towards factoring. R = randomized,
  D = deterministic.
- **F→**: computing the invariant from the factorisation.
- **Prime-easy**: computable at a prime argument or level.
- **Cls**: the class (a), (b) or (c) from §0.

| Invariant | →F | F→ | Prime-easy | Cls | Reference / where checked |
|---|---|---|---|---|---|
| φ(N) | R all N; D under ERH (Miller 1976, JCSS 13); D for pq (s = N+1−φ); D unconditional for many N via LLL (Morain–Renault–Smith, arXiv 1802.08444) | D | p−1 | a | `reduce_phi` |
| λ(N), or any multiple | R; D under ERH | D | yes | a | Miller 1976; Morain–Renault–Smith |
| σ(N) | R all N, no ERH | D | yes | a | Bach–Miller–Shallit, SIAM J. Comput. 15 (1986) |
| ψ(N) = #P¹(Z/N) = [SL₂Z : Γ₀(N)] | D for pq | D | yes | a | `reduce_psi` |
| r₄(N) = 8σ(N), N odd | as σ | D | yes | a | Jacobi; brute-force test |
| Pillai Σ_{k≤N} gcd(k, N) = (2p−1)(2q−1) | D for pq | D | yes | a | `reduce_pillai`. N terms, each polylog by Euclid |
| genus g(X₀(N)) = dim S₂(Γ₀(N)) | D for pq (≤ 4 candidates) | D | yes | a | `reduce_genus_X0` |
| dim S₁₆(Γ₀(N)) − dim S₄(Γ₀(N)) = ψ(N) | D for pq | D | yes | a | codex nonsmooth_1 (Stein, Prop 6.1) |
| N/φ(N) = ζ_{Spec Z/N}(1) | as φ | D | yes | a | codex nonsmooth_1 |
| #E(Z/N), good reduction | D from a curve and one twist (Dieulefait–Urroz, arXiv 1911.11004, Thm 1); any single curve under ERH (their Thm 5); R from four twists, no ERH (codex). Kunihiro–Koyama (Eurocrypt 1998) is randomized and rests on a distribution assumption (D–U Remark 2) | D (Schoof + CRT) | Schoof 1985 | a | polyimpl §2.1; also order of a point on E(Z/N) (Martín–Morillo–Villar 2001) |
| #T(Z/N), tori. Conic x²+y² = 1: (p − ε_p)(q − ε_q), ε = (−1/·) | D for pq (two-case analysis; quasi-polynomial extension of §5.3) | D | yes | a | `reduce_conic`, brute-force test |
| τ(N), τ(N²) | D; fails only if τ(p) = τ(q) = 0. Density one unconditionally; all RSA N under Lehmer's conjecture (Bach–Charles Thm 2.2, Cor 2.8) | D (E–C, Ann. Math. Stud. 176, Thm 15.1.1, plus the factorisation) | E–C, poly(log p) | a | Bach–Charles arXiv 0708.1192; `reduce_tau_pair` (denominator form, tested incl. N = 15, 35) |
| a_f(N), a_f(N²), fixed non-CM eigenform | same | D given the factored index (Bruin 2011, GRH, Las Vegas) | yes | a | Bach–Charles; Bruin, Publ. Math. Besançon 2011 (via codex) |
| C(N) = Σ_{t²<4N} H(4N−t²) = Σ_{ab=N} max(a, b) = 2N + 2q | D for pq | D | C(p) = 2p | a | Kronecker–Hurwitz; Popa–Zagier arXiv 1604.02822 Thm 1; `reduce_hurwitz_sum`, brute force for 7 N |
| Frobenius operator Φ_B on étale Z/N-algebras | D, or R plus CRT for arbitrary B (Prop FROB) | D | x ↦ x^p | a | §5.2; `reduce_quadratic_frobenius`, `reduce_frobenius_fixed_space` |
| h(−4N) with class-group arithmetic | R (ambiguous forms: Shanks, Schnorr–Lenstra); D under GRH (Schoof 1982 §6, via codex) | not known | **no**: h(−p) is L[1/2] under GRH (Hafner–McCurley 1989) | b | `reduce_class_number` |
| square roots in Cl(−4N) | R (two roots differ by an ambiguous class) | R | given the factors | a | Rabin-type folklore |
| square roots mod N; discrete log mod N | R | R (+ prime DLP) | yes / no | a / b | Rabin 1979; Bach 1984; Shmuely 1985 |
| regulator of Q(√N) | R (ambiguous form at half period) | not known | **no** | b | Shanks; Murru–Salvatori arXiv 2409.03486 (smooth relations) |
| ord_N(a), random a | R; D under ERH | D, plus factoring p−1 | **no** | b | Miller 1976; Morain–Renault–Smith (Shor's oracle) |
| Σ_{a<N/2} (a/N) | gives h(−N) | — | **no**: at a prime it is h(−p) | b | `class_number_half_sum` (tested against enumeration) |
| Jacobi (d/N); power-residue symbols at (N) | **0 bits** beyond N mod the conductor | — | yes | c | Jager–Schwenk, J. Cryptol. 2013 (generically factoring-hard, easy by reciprocity); de Boer–Pagano ISSAC 2017 |
| d(N), μ(N), ω(N) | not known (constant on RSA N) | D | — | — | Shallit–Shamir 1985; polyimpl §2.5 |

## 2. Organising principle (taxonomy, not a theorem)

What can be computed from N in polynomial time? Three tools:
- ring operations in Z/N-algebras, plus gcds;
- the integer representation through Euclid-type recursions (reciprocity, floor sums, continued
  fractions);
- lattice reduction.

GPT-6.1-sol's diagnosis (nonsmooth_1 §4), adopted as a taxonomy: these compute **aggregates that
forget the CRT components**. Recovering a component-sensitive quantity needs an "unpaid" operation:
- **(U1)** select a local component;
- **(U2)** construct local Frobenius or order data;
- **(U3)** evaluate a large arithmetic sum without expanding it.

How each is blocked:
- **(U2), Frobenius operators.** Factoring by Prop FROB.
- **(U2), orders.** Prop OA reduces scalar annihilation to smoothness or to an anomalous
  coincidence.
- **Reciprocity.** The one known non-generic shortcut (Jager–Schwenk) evaluates only abelian
  aggregates, which depend on N mod m (§5.4).
- **(U3).** The only door left open. Every concrete (U3) sum we found is divisor detection in
  disguise:
  - Pillai, r₄, C(N);
  - codex's Jacobi autocorrelation A_N(H) = H − p;
  - the Heilbronn–Merel sums for T_N;
  - the Ramanujan-sum series for σ.

  The best (U3) evaluators known are N^{1/3} (Hiary; Rubinstein). Lattice and collision methods
  dominate them (N^{1/4}, then N^{1/5}).

## 3. Paradigm verdicts

| Paradigm | Verdict | Precise reason |
|---|---|---|
| Edixhoven–Couveignes / Bruin at composite N | **no route**; the operator step is **REFUTED for this model** | a(N) = a(p)a(q) is a product of traces at two Frobenius classes, not one element's trace. The only p-dependent step of E–C is identifying Frob_p in Gal(K_ℓ/Q). Producing the Frobenius operator on an algebra with (disc/N) = −1 (e.g. ψ_ℓ with (N/ℓ) = −1) factors N (FROB(b, d)). Whether E–C's *output* (a conjugacy class) needs the full operator is not proved. Bruin needs the index factored (KNOWN). Saxena–Venkatesh (arXiv 2511.05272) is again prime-only |
| a(N) mod ℓ, composite N | no route | Exceptional ℓ: σ_k(N) mod ℓ is a degree-k polynomial in s mod ℓ (≤ k candidates). Non-exceptional ℓ: Chebotarev data of a pair of classes. Enough ℓ reconstruct a(N), which is Bach–Charles |
| Eichler–Shimura over Z/N | **REFUTED for this model** | T_N = T_pT_q. T_p ≡ F + V mod p uses Frob_p; T_q is a correspondence of degree q+1 ≈ √N |
| Dokchitser resolvents | no route (scoped) | They need x ↦ x^p, and x ↦ x^N does not permute the roots (codex). FROB covers producing Φ, not resolvent values directly |
| Schoof/SEA "pretending N is prime" | **REFUTED for this model** | The Frobenius operator on (Z/N)[x]/ψ_ℓ factors N for every ℓ with (N/ℓ) = −1 (FROB(b, d)). Running it with x^N leaks only at the primality-test rate (W-13/W-14) |
| Kedlaya / Lauder / Satoh | no route | A Frobenius lift mod N^k reduces mod N to Φ (FROB). The deformation base point needs Gauss sums over F_p; Gross–Koblitz needs Γ_p, whose "N-adic" analogue is the factorial unit predicate (Track C). Big-Witt Frobenius exists over any ring but is not the hidden operator (codex) |
| Harvey average polynomial time | **REFUTED for this model** | The amortisation is over all primes ≤ X with a shared tree; reaching an unknown p needs X ≥ p. A composite N gives a product of length N |
| Trace formulas | reduced to divisor detection (§4.1) | Geometric side: √N class numbers. Spectral side: Hecke eigenvalues at N plus divisor sums. Charles (Ramanujan J. 11, 2006) computes τ(N) factor-free in N^{1/2+ε} (via codex; GRH status UNVERIFIED) |
| Dequantising Shor | no route | Query bounds: Ω(r^{1/3}/√m) (Cleve, Inf. Comput. 192, 2004), and Ω̃(√r) for polynomially bounded jumps (Chakraborty–García-Soriano–Matsliah, ICS 2011). These are query-model bounds, not bounds for explicit modular multiplication (codex). Tang-style dequantisation needs low rank. The exploited non-black-box structures are smoothness and index calculus. Cheap-order elements (±1, public roots of unity) have equal local orders; a nontrivial √1 already is a factor |
| Lattices / Coppersmith | no invariant supplies the hint | It needs ≈ n/4 coset-shaped bits (Track E E1); per-modulus data give 0 bits (F1). Class-(a) values mod small ℓ give at most s mod ℓ (§5.5); s mod ℓ for all ℓ ≤ 0.35n gives s by CRT (polyimpl §2.5), but no single residue is computable |
| Varieties with global formulas | **REFUTED** for (quasi-)polynomial counts (§5.3) | Either #V(Z/N) = c²N^r, or it is factoring-equivalent. Other counts: no factor-free method known |
| Generic rings | equivalence, not a floor | Aggarwal–Maurer 2009; Jager–Schwenk 2013; Boneh–Venkatesan 1998 |
| Class groups / infrastructure | smoothness-based or exponential | Hafner–McCurley, Schnorr–Lenstra, SQUFOF2 (Bradford–Wagstaff), Murru–Salvatori: smooth relations. SQUFOF: N^{1/4}. Shanks/Schoof under GRH: N^{1/5} |
| Analytic / lattice-point counting | exponential | Hiary (Math. Comp. 2016) and Rubinstein (Integers 2013): N^{1/3}. Taylor of degree j plus exact degree-j floor sums would give N^{1/(2j+2)} *if* those sums were polylog. But complete quadratic floor sums mod p compute h(−p) (class b), and Coppersmith/HH beat this architecture. Umans–Wang (arXiv 2511.10851) give a conjectural deterministic N^{1/6} |
| CM / anomalous curves | **real but thin** (§4.2) | Prop OA; E1/E2 |
| Reciprocity amplification (codex) | **REFUTED** | A_N(H) = H − p: an interval divisor selector (codex nonsmooth_1, proved) |
| Cartier compression of Δ mod ℓ | **REFUTED** in the uniform form | Mod 2, depth-k sections have rank ≥ 2^{k−3} (codex, proved) |

## 4. Candidates developed

### 4.1 Hurwitz aggregate (codex's sole survivor): reduced to the Coppersmith hint

**Sketch.** Write C(N) = Σ_t H(4N − t²) through L-values or a Rademacher-type series. Compress the
t-sum (theta, contour or Shintani methods) with certified error E. Since C(N) = 2N + 2q, an error
E ≤ 2N^{1/4} gives q to within N^{1/4}, and Coppersmith (Track E E1) finishes.

**Analysis.**
1. **(PROVED) Restatement.** C(N) is an exactly linear proxy for q: q = (C − 2N)/2, as codex notes.
   So "approximate C(N) to N^{1/4}" *is* "find the top half of q". The candidate supplies a name
   for the hint, not a mechanism.
   - It is better conditioned than approximating σ(N). Going from s back to p is ill-conditioned
     near p = q (codex), so the two hints are not literally identical.
2. **(PROVED) One expansion fails exactly.** Ramanujan's σ(N) = Nζ(2) Σ_c c_c(N)/c² can be
   rewritten as σ(N) = N + Nζ(2) Σ_c (c_c(N) − μ(c))/c². If gcd(c, N) = 1 then c_c(N) = μ(c), so
   every term with c < p vanishes. Truncating at any X < p returns N exactly, missing p + q + 1;
   reaching c = p costs ≍ √N (codex nonsmooth_2; tested). With C(N) = 2σ(N) − λ₁(N), the same
   holds for C.
3. **(HEURISTIC) Euler-factor version.**
   - For odd ℓ ∤ N, the t-average of the imprimitive local factor at ℓ is
     1 − (1 + (N/ℓ))/(ℓ(ℓ² − 1)) (codex; checked by critic 1). It depends on N only through N mod ℓ,
     which is the information floor.
   - At ℓ = p the complete character sum is p − 1, against −1 at a typical prime: the mean local
     factor is ≈ 1 + 1/p.
   - That bias is where p-information enters. It is not a derivation of the exact 2q, because
     Hurwitz numbers include imprimitive forms, and the λ₁ correction cancels a p-contribution
     (critic 1, codex).
   - This matches Gekeler's (IMRN 2003) local-density description of Σ_t H(4p − t²) at primes.
4. **(HEURISTIC) Other representations.**
   - S₂(SL₂Z) = 0 removes only the holomorphic cuspidal term of the trace formula.
   - Hurwitz numbers have genuine Rademacher / mock-modular representations (Duncan–Mertens–Ono
     arXiv 1702.03516, Prop 5.2; arXiv 2505.05574).
   - Their Kloosterman sums S(·, N; c) see N only mod c; we expect the 2q to sit at c ≡ 0 (mod p),
     as in item 2, but have not proved it. No complexity bound for any of them is known.

**Verdict.** Not a mechanism. Any fast certified approximation of C(N) is a fast factoring
algorithm, by the identity. Nothing in the representations checked suggests one exists, and the
one fully analysed expansion provably needs moduli ≥ p.

### 4.2 Anomalous CM (the only smoothness-free order-annihilation branch)

**Algorithm (KNOWN).** References: Cheng 2002, ePrint 2002/109; Shirase 2017, ePrint 2017/403;
Aikawa–Nuida–Shirase, IEICE E102-A (2019); Sedlacek et al., SECRYPT 2019.
- If p | N and 4p − 1 = Ds², then t_p = ±1. So E_D mod p, or its twist, has exactly p points.
- [N]P on a random x-only point vanishes mod p with probability 1/2 (1/6 for D = 3).
- For h(−D) > 1 the method works over (Z/N)[j]/H_D(j).

`src/nonsmooth/cheng.py` implements the h = 1 cases (D ∈ {3, 11, 19, 43, 67, 163}) operationally.
D = 7 is impossible, since 4p − 1 ≡ 3 (mod 8) forces D ≡ 3 (mod 8). That is proved, and holds for
all 80 000 census primes.

**E1 (pre-registered; PASS).**

| arm | n | success | structural | tries at the right D | mul / try / bit |
|---|---|---|---|---|---|
| special | 128, 256, 512 | 60/60 each | 180/180 | 2.03–2.50 (D > 3), 6.27 (D = 3) | 23.6–23.9 |
| control | 128, 256, 512 | 0/60 each | — | — | — |

Predictions P1–P3 held as written. At 512 bits a success takes 0.9 s in Python.

**E2 (pre-registered; PASS).** 20 000 random primes per size; D = squarefree part of 4p − 1.

| bits | median log D / log p | F(0.25) | F(0.5) | F(0.75) | F(0.9) |
|---|---|---|---|---|---|
| 24 | 1.084 | 4e−4 (8) | 4.5e−3 | 3.14e−2 | 0.122 |
| 32 | 1.063 | 5e−5 (1) | 1.2e−3 | 1.68e−2 | 0.081 |
| 40 | 1.050 | 0 | 5e−5 | 8.9e−3 | 0.064 |
| 48 | 1.042 | 0 | 5e−5 | 3.7e−3 | 0.049 |

- Slopes of log₂F per bit:
  - θ = 0.75: −0.127, against a predicted −0.125 (window [−0.19, −0.06]).
  - θ = 0.9: −0.053, against −0.05.
  - θ = 0.5: −0.30, against −0.25; reported only.
- Prefactor c₀.₇₅ = 0.23–0.27.
- Q1–Q4 held as written.
- So F(θ) ≈ 0.25·p^{−(1−θ)/2} is **experimentally supported**, consistent with the elementary
  upper bound O(√(XY)) (codex).

**Cost.**
- Covering a fraction F of primes needs D up to Y ≈ 16F²p. The pre-registration wrote Y ≈ F²p; the
  constant 16 comes from the measured prefactor 0.25 (disclosed refinement).
- A portfolio that processes every class polynomial does arithmetic in algebras of dimension
  h(−D). That costs ≥ Σ_{D≤Y} h(−D) ≍ Y^{3/2} (the classical Gauss–Mertens average), i.e.
  ≳ 64F³p^{3/2}, with H_D supplied free.
- Computing the H_D themselves pushes the cost towards Y².
- Even at O(1) per D the cost is ≥ #D ≍ Y ≍ F²p.
- These bound portfolios of this shape only, not every compressed CM algorithm (codex).

**Verdict.** A correct special-purpose family, of density zero, exponential for constant
coverage.

The non-CM variant ("anomalous ECM", a_p = 1) has coverage ≍ p^{−1/2} per random curve, under a
random-curve model (HEURISTIC; Lang–Trotter constants can vanish, e.g. curves with rational
2-torsion have a_p even). That is cost ≍ N^{1/4}, no better than ρ.

### 4.3 Smaller candidates, killed on the spot
- **Hecke action on modular symbols.** The Heilbronn–Merel sum over x mod N has N terms. It is a
  (U3) sum that would compute τ(N); we found no structure that shortens it.
- **Jacobi interval sums.** These compute h(−N), so they are class (b).
- **Ideals strictly between (N) and O_K.** Naming one needs a root of f mod p (Rabin-type).
  Power-residue symbols at (N) are abelian (§5.4).

## 5. Barriers (statements and proofs)

### 5.1 Prop OA (order annihilation). PROVED

*Setting.*
- G is a smooth commutative group scheme of finite type over Z[1/M₀] with connected fibres, and
  p ∤ M₀ is odd.
- The method computes [m]x for x ∈ G(Z/N) and a fixed integer m, chosen independently of x mod p.
  It tests [m]x ≡ 0 (mod p) via a gcd.
- x mod p is uniform on G(F_p). For x-only mixing of twists, treat each twist separately.

(a) #G(F_p) = #T(F_p)·p^{dim U}·#A(F_p). Here #T(F_p) is prime to p and ≥ (p − 1)^{dim T}.

*Proof.*
- Chevalley gives 0 → L → G → A → 0, with L connected linear and A an abelian variety.
- Over the perfect field F_p, L ≅ T × U with U split unipotent, so #U(F_p) = p^{dim U}.
- Lang gives H¹(F_p, L) = 0, so #G(F_p) = #L(F_p)·#A(F_p).
- #T(F_p) = |det(p − F | X*(T))|, where F is the finite-order Galois action, with root-of-unity
  eigenvalues ζ. Hence ∏|p − ζ| ≥ (p − 1)^{dim T}, and det(p − F) ≡ det(−F) = ±1 (mod p). ∎

(b) Write G(F_p) ≅ ∏ Z/d_i. Then Pr_x[[m]x = 0] = ∏ gcd(m, d_i)/d_i (codex). If a prime ℓ | #G(F_p)
does not divide m, this is ≤ ℓ^{−v_ℓ(#G(F_p))} ≤ 1/ℓ. ∎

(c) Corollary for m = N^k·S, S composed of primes ≤ B.
- Success probability > 1/B forces every prime ℓ > B dividing #G(F_p) to lie in {p, q}.
- So #G(F_p) = p^a q^b × (B-smooth).
- Necessary, not sufficient: exponents of small primes and separation at q also matter.
- With no smooth part (#G(F_p) a {p, q}-group, p ≥ 3): dim T = 0, and #A(F_p) is a {p, q}-number.
  - Either A(F_p) is a p-group (anomalous),
  - or q | #A(F_p). Codex's example: N = 35, #E(F₅) = 7 = q.
- For dim A = 1 and balanced N, q | #E(F_p) forces q ≤ p + 1 + 2√p, the Fermat regime.
- The unipotent part is killed at both p and q, so it never splits.

*Scope (codex).* Not covered: arbitrary computable m (heuristically, a "random" integer of size
≈ p^{dim} must divide m, which favours highly composite m); structured exponents; stage two;
subgroup selection; collision and order searches. This is not a factoring lower bound. The torus
branch is Bach–Shallit's cyclotomic framework (Math. Comp. 52, 1989), and the abelian branch is
Lenstra's ECM (Ann. Math. 126, 1987). The packaging with the anomalous branch: "no prior
formulation found by us".

(d) Coverage.
- CM with D ≤ Y: §4.2.
- Fixed non-CM E: #{p ≤ X : a_p = 1} ~ C_{E,1}√X/log X with C_{E,1} ≥ 0 (Lang–Trotter,
  conjectural).
- Dimension g ≥ 2: no universal exponent. A = E^g is anomalous exactly when E is (codex), so the
  answer depends on the family; UNKNOWN.

### 5.2 Prop FROB (Frobenius operators over Z/N). PROVED

*Setting.*
- N = pq is odd and squarefree.
- B is a finite étale (Z/N)-algebra with an explicit basis.
- Φ_B is the componentwise Frobenius, i.e. Φ_B ⊗ F_p = Frob_p. At prime N it is x ↦ x^N.

(a) B = (Z/N)[y]/(y² − d), (d/N) = −1. Then Φ(y) = εy, where ε = CRT((d/p), (d/q)) is a
nontrivial √1. So gcd(ε − 1, N) is a factor (the Rabin/Miller split).
- Finding d: randomized, expected O(1) trials; deterministic polynomial time under ERH.
- `reduce_quadratic_frobenius`, tested.

(b) B = (Z/N)[x]/(f), f monic with disc f a unit mod N, and (disc f / N) = −1. Then a suitable
O(n³) elimination on Φ_B − 1 over Z/N extracts a factor.

*Proof.*
- B ⊗ F_p = ∏_{i=1}^{r_p} F_{p^{d_i}}, so dim ker(Frob − 1) = r_p. (This is the Berlekamp Q − I
  matrix.)
- Frob permutes the roots with sign (−1)^{deg f − r_p}. By Pellet–Stickelberger (Swan, Pacific
  J. Math. 12, 1962) that sign is (disc f / p).
- Since (disc/p) ≠ (disc/q), we get r_p ≠ r_q, so ker(Φ_B − 1) ≅ F_p^{r_p} × F_q^{r_q} is not
  free.
- Elimination with unit pivots alone would give a free kernel. So some entry met is a nonzero
  nonunit. ∎
- `reduce_frobenius_fixed_space`, tested on Schoof's algebra (Z/N)[x]/ψ₅.

(c) B ⊇ (Z/N)[ζ_ℓ] with ζ_ℓ explicit, ℓ ∤ N, (N/ℓ) = −1, so p ≢ q (mod ℓ). For a ≡ p (mod ℓ),
Φ(ζ) − ζ^a vanishes mod p only, and the gcds of its coordinates with N factor it. This takes
poly(ℓ, log N) (cf. Lenstra 1985; Bach–Shallit 1989).

(d) **Schoof's algebra.** Assume N odd, gcd(ℓΔ_E, N) = 1, p, q > 3, and ψ_ℓ made monic.
- The sign character of GL₂(F_ℓ) on (F_ℓ² − 0)/±1 is Legendre(det).
  - For ℓ ≥ 5, SL₂(F_ℓ) is perfect, so the character factors through det: it is trivial or
    Legendre(det).
  - Evaluating diag(g, 1) for a generator g decides which. The line y = 0 gives one cycle of length
    (ℓ−1)/2. The (ℓ−1)/2 classes with y ≠ 0 each give an (ℓ−1)-cycle. So the sign is −1 = (g/ℓ).
  - Checked by brute force for ℓ = 5, 7, 11, 13 (tests) and every odd ℓ ≤ 31 (critic 1).
- Hence the sign of Frob_p on the roots of ψ_ℓ is (p/ℓ), and (disc ψ_ℓ / N) = (N/ℓ).
  - Tests: ψ₅ and ψ₇ on two curves each, 100 primes each, 0 mismatches. A further 400-prime
    scratch check and critic 1's check (ℓ = 3, 5, 7, three curves) also showed none.
  - Brunyate–Clark (Ramanujan J. 37, 2015) is a Zolotarev-type analogue; Atkin's parity rule in SEA
    (Schoof, JTNB 7, 1995) is the modular-polynomial analogue.
- So for every ℓ with (N/ℓ) = −1, the Frobenius operator Schoof needs on (Z/N)[x]/ψ_ℓ factors N
  via (b).

*Scope (codex).*
- This is hardness of producing the *operator*, as a matrix or as images of generators. B = Z/N
  has Φ = id, and arbitrary B need a separation argument; (d) supplies one for Schoof's algebra.
- Outputs that are only numbers or conjugacy classes (a(N) mod ℓ, Dokchitser resolvent values,
  p-adic matrices) are not covered.
- Lenstra 1985 and Grantham 2001 are the nearest primality-testing literature. FROB is folklore.

### 5.3 (Quasi-)polynomial point-count dichotomy. PROVED

**Polynomial case (codex nonsmooth_1 §4).** Suppose #V(F_p) = P(p) = Σ_{i=r}^{d} c_i p^i with
c_r c_d ≠ 0.
- r = d: #V(Z/pq) = c_r²N^r, which carries no information.
- r < d: #V(Z/pq) is a polynomial in s of degree d − r with leading coefficient c_r c_d N^r. Its
  integer roots give s.

**Quasi-polynomial case.** Suppose P depends on p mod m through a character, e.g. the conic
(p − ε_p)(q − ε_q) with ε_pε_q = (−1/N) known. Enumerate the finitely many sign patterns consistent
with the known character values; each gives s or q − p. This is `reduce_conic`. Routine (cf.
Bach–Miller–Shallit; Bach–Shallit).

### 5.4 Abelian aggregate theorem. PROVED (codex; Kronecker–Weber)

For K/Q abelian of conductor m and (N, m) = 1, Art_{K/Q}((N)) depends only on N mod m. With
Track E F1, character values and power-residue symbols at (N) carry 0 bits beyond N mod m.

*Scope.* This is a redundancy statement for these observables. It does not classify everything
computable by reciprocity. Also tr ρ(Art(N)) ≠ a(N): for example σ_k(pq) = 1 + N^k + p^k + q^k,
whereas the abelian aggregate is 1 + N^k.

### 5.5 Prop RANK (minimal primitive, linear divisor formulas). PROVED (elementary)

**Statement.** Let ℓ be an odd prime with ℓ ∤ N. From φ(N), σ(N) or ψ(N) mod ℓ (or r₄(N) mod ℓ)
one gets s mod ℓ, and hence {p, q} mod ℓ as the roots of x² − sx + N. When (N/ℓ) = 1 this decides
the common value (p/ℓ) = (ℓ*/p), i.e. whether ℓ* is a square mod N. For ℓ = 3 this is whether N is
of the form x² + 3y². The 2-adic analogue (s mod 8) is ε (Track E E2′).

**Not universal (codex counterexamples).**
- A curve with rational ℓ-torsion has #E(Z/N) ≡ 0 (mod ℓ) always.
- For G_m², φ(N)² mod 7 confuses {2, 4} with {6, 6}.
- Degree-k Eisenstein coefficients (τ mod 691 = σ₁₁ mod 691) give ≤ k candidates.

**Status of fixed-base QR.** Factoring decides it trivially. No reduction from factoring to it is
known (for base −1 this is Seurin's 2-Φ/4-hiding assumption, PKC 2014 / ePrint 2013/256), and no
algorithm for it faster than factoring is known.

### 5.6 Prop EIS. Summary of §4.1

| Claim | Status |
|---|---|
| C(N) = 2σ(N) − λ₁(N) = 2N + 2q | PROVED (brute-force checked) |
| An approximation of C(N) is a linear q-hint | PROVED |
| The Ramanujan-sum expansion has zero terms below p | PROVED, tested |
| The Euler-factor picture | HEURISTIC |
| Other representations behave the same way | HEURISTIC |

## 6. What would change the verdict

- A (U3) evaluator that is sub-N^{1/5} for a divisor-detecting sum and does not localise at p.
  Examples:
  - an exact polylog algorithm for degree-2 floor sums Σ_{h<H}⌊P(h)⌋ with H ≪ the denominator;
  - a Rademacher/Kloosterman truncation for C(N) whose tail can be bounded without moduli
    c ≡ 0 (mod p).
- A sub-factoring algorithm for fixed-base QR mod N (§5.5).
- A factor-free method that outputs a(N) mod ℓ without producing the Frobenius operator (FROB
  scope).

Nothing in the literature we checked points to any of these.

## 7. Reproduction and deviations

```
.venv/bin/python -m pytest -q tests/test_nonsmooth.py
.venv/bin/python -m src.nonsmooth.run_experiments e1   # ~4 min
.venv/bin/python -m src.nonsmooth.run_experiments e2   # ~1 min
```

Disclosed after the pre-registration:
- Tests grew from 46 to 67:
  - FROB(b) on ψ₅; the disc ψ_ℓ checks; the GL₂ sign;
  - the Ramanujan-sum test; conic counts; the extra class-number cases.
- Code changes after E1/E2, none touching cheng.py or census.py:
  - `reduce_class_number` now also tries a ∓ b/2 for ambiguous forms (a, b, a). This was a bug
    found by critic 1: about 10% of order-2 hits were missed.
  - `reduce_tau_pair` now uses the Bach–Charles denominator, so it no longer fails on N = 35.
- E2 pricing uses Y ≈ 16F²p, not the pre-registered F²p. This is a refinement from the measured
  prefactor (§4.2).
- θ = 0.25 counts were added to the table.

## 8. GPT-6.1-sol adversarial round (codex/nonsmooth_2.out): adjudication

| Codex finding | Adjudication |
|---|---|
| OA(a) true | kept |
| OA(b) needs the sampling hypothesis | added |
| **OA(c) false** (N = 35, #E(F₅) = 7 = q) | restated with {p, q}-parts and as a necessary condition |
| "Only smoothness-free case" overclaimed | scope paragraph added |
| Lang–Trotter constant may vanish | fixed |
| Dimension g ≥ 2 "rarer" false (E^g) | now UNKNOWN |
| FROB(a) "deterministic" needs ERH | fixed |
| FROB(b) true; attribute Pellet–Stickelberger | done |
| FROB(c) needs ℓ ∤ N and explicit ζ | added |
| FROB(d) needs per-algebra separation | given for Schoof via the GL₂ sign; other outputs scoped out |
| EIS identity true | kept |
| **EIS(i) false** (imprimitive vs Hurwitz; average depends on N) | Euler-factor argument relabelled HEURISTIC, with the correct average |
| **EIS(ii) false** conclusion | removed |
| **EIS(iii) false** as a kill | rephrased as a reduction; Ramanujan-sum exact limitation added |
| **RANK false** universally | restricted to linear divisor formulas, counterexamples listed |
| Global exclusion is a literature verdict only | §0 says so |
| Missed SQUFOF2 and Murru–Salvatori | added; both smoothness-based |

## 9. Opus critics: adjudication

**Critic 1 (correctness/cost; ACCEPT WITH FIXES).**
- Applied:
  - EIS relabel (1, 2);
  - OA(c) restatement (3);
  - FROB hypotheses and the diag(g, 1) cycle proof (4);
  - RANK restriction (5);
  - quasi-polynomial conic (6);
  - Gauss–Mertens attribution and the "free class polynomial" line (7);
  - prereg disclosures (8);
  - class-number bug (9);
  - Harvey 2021 attribution (11).
- Item 11's doubt about the ICS 2011 citation: we verified the title and authors on the ICS 2011
  site, and critic 2 confirmed them. Kept.
- E1/E2 numbers were recomputed by the critic from raw data and agree.

**Critic 2 (literature).**
- Applied:
  - Kunihiro–Koyama is randomized; Dieulefait–Urroz Thm 1 is deterministic;
  - Bach–Charles fails only at τ(p) = τ(q) = 0, and the code was upgraded;
  - Popa–Zagier replaces the unverified Gross–Zagier (4.4) citation;
  - Harvey 2021, Hittmeir 2018;
  - Cleve's jump parameter;
  - Morain–Renault–Smith; discrete log and square roots mod N;
  - Swan, Lenstra, Bach–Shallit, Brunyate–Clark, Atkin parity, Gekeler, Seurin;
  - "no reduction in either direction" corrected;
  - §0 overgeneralisations softened;
  - Rubinstein, Shanks–Schoof N^{1/5}, Umans–Wang, Saxena–Venkatesh.
- Left as UNVERIFIED: whether Charles 2006's N^{1/2+ε} bound is GRH-conditional.
- Critic 2 found no missed non-smooth subexponential method.
