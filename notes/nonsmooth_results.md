# nonsmooth: factoring paradigms not based on smooth numbers (phase 16, route 3)

Agent `nonsmooth`, 2026-10-20.
- Pre-registration: `notes/nonsmooth_prereg.md` (commit `3b1b19e`, before E1/E2 ran).
- Code: `src/nonsmooth/{equiv,forms,cheng,census,run_experiments}.py`. Tests:
  `tests/test_nonsmooth.py` (52 pass).
- Raw data: `experiments/raw/nonsmooth_e{1,2}.jsonl`. Summary:
  `experiments/analysis/nonsmooth_summary.json` (commit `4ea0295`).
- GPT-6.1-sol: blind round `codex/nonsmooth_1.{md,out}`; adversarial referee
  `codex/nonsmooth_2.{md,out}` (§8).
- Opus critics: §9.
- Labels: **KNOWN** (with reference), **PROVED** (proof here; all routine, "no prior formulation
  found by us" at most), **CONDITIONAL**, **HEURISTIC**, **experimentally supported**,
  **SPECULATIVE**, **REFUTED for this model**, **UNVERIFIED**.

## 0. Verdict (plain language)

**No smoothness-free route to general factoring below L[1/3], or even subexponential, survives.**
The organising question was: which factoring-equivalent invariant of N might have an unexpectedly
cheap computation that does not factor N first? We mapped 20 such invariants and tested every
listed paradigm against them. The invariants fall into three classes:
- **(a) easy at primes, products of local factors at composite N.** φ, λ, σ, ψ = #P¹(Z/N), r₄,
  Pillai's Σgcd, #E(Z/N), tori, τ(N), a_f(N), dim S_k(Γ₀(N)), the Hurwitz aggregate C(N).
  Every fast method at a prime p consumes the **Frobenius at p**, as an operator on an étale
  algebra. Prop FROB (§5.2, proved) shows that this operator over Z/N factors N. That covers
  Schoof's own algebra (Z/N)[x]/ψ_ℓ whenever (N/ℓ) = −1. Every known global formula for these
  invariants is a sum of length ≥ N^c whose evaluation is divisor detection.
- **(b) hard even at primes.** h(−p), the regulator, ord_p(a). A route through them must solve a
  problem not known to be easier than factoring.
- **(c) computable by reciprocity.** Jacobi and power-residue symbols at (N), and more generally the
  Artin symbol of (N) in an abelian field. These depend only on N mod the conductor (abelian
  aggregate theorem, §5.4), so they carry zero bits beyond it (Track E F1).

**The genuinely non-smooth methods that do exist are exponential, special-purpose or quantum.**
- Exponential: Fermat/Lehman, Pollard–Strassen/BGS, SQUFOF, Hiary's N^{1/3}, Coppersmith N^{1/4},
  Harvey–Hittmeir N^{1/5}.
- Special-purpose: CM/anomalous curves (Cheng's 4p−1 method).
- Quantum: Shor.

**The most promising candidate is codex's survivor, and it dies (§4.1).** The idea: approximate the
Hurwitz class-number aggregate C(N) = Σ_t H(4N − t²) = 2N + 2q to within N^{1/4} cheaply, then run
Coppersmith. But the 2q term is exactly the Euler-factor bias at ℓ = p. Weight 2, level 1 has no
cusp forms, so C(N) has no representation other than "class numbers = divisor sum". The candidate
is therefore the Coppersmith hint (s = p + q to within N^{1/4}) in disguise.

**Order-annihilation methods** (p−1, p+1, ECM, tori, abelian varieties) obey a proved dichotomy
(Prop OA, §5.1): either the prime-to-p part of #G(F_p) is smooth, or G(F_p) is a p-group, an
anomalous abelian variety. The second branch is real: our operational Cheng-type factorer splits
180/180 special N up to 512 bits and 0/180 random controls (E1). But it is thin. Coverage decays
like √(Y/p) when the CM discriminant is restricted to D ≤ Y. The pre-registered slope test
measured −0.127 against a predicted −0.125 (E2). Constant coverage would cost at least p^{3/2}.

**The single minimal primitive** that every class-(a) route needs, even modulo one small prime ℓ,
is a fixed-base quadratic-residuosity bit: "is ℓ* a square mod N", or "is N = x² + 3y²" (§5.5).
That is the Track E ε-type bit, an open problem with no known algorithm faster than factoring.

## 1. Equivalence map

Columns:
- "→F": what the invariant's value (or oracle) gives towards factoring.
- "F→": the reverse direction.
- "Prime-easy": a known method computes it at a prime argument or prime level.
- "Class": (a), (b) or (c) as in §0.

| Invariant | →F | F→ | Prime-easy | Class | Reference / where checked |
|---|---|---|---|---|---|
| φ(N) | R all N; D under ERH; D for pq (s = N+1−φ) | D | φ(p)=p−1 | a | Miller 1976 JCSS 13; `equiv.reduce_phi` |
| λ(N), or any multiple of it | R; D under ERH | D | yes | a | Miller 1976 |
| σ(N) | R all N, no ERH | D | yes | a | Bach–Miller–Shallit, SIAM J. Comput. 15 (1986) |
| ψ(N) = #P¹(Z/N) = [SL₂Z : Γ₀(N)] | D for pq | D | yes | a | `reduce_psi` |
| r₄(N) = 8σ(N) (odd N) | as σ | D | yes | a | Jacobi; brute-force test |
| Pillai Σ_{k≤N} gcd(k,N) = (2p−1)(2q−1) | D for pq | D | yes | a | `reduce_pillai`. Each term is polylog by Euclid; there are N terms |
| genus g(X₀(N)) = dim S₂(Γ₀(N)) | D for pq (≤ 4 candidates) | D | yes | a | `reduce_genus_X0` |
| dim S₁₆(Γ₀(N)) − dim S₄(Γ₀(N)) = ψ(N) | D for pq | D | yes | a | codex nonsmooth_1 (Stein, Prop 6.1). The elliptic and cusp terms cancel |
| N/φ(N) = ζ_{Spec Z/N}(1) | as φ | D | yes | a | codex nonsmooth_1 |
| #E(Z/N), good reduction | R from counts of E and three twists, with no ERH (codex four-twist argument); deterministic versions in Kunihiro–Koyama (Eurocrypt 1998, LNCS 1403) and Dieulefait–Urroz (arXiv 1911.11004; codex flags a gap in its projective proof) | D (Schoof + CRT) | Schoof 1985 | a | polyimpl §2.1 |
| #T(Z/N), algebraic tori (e.g. the conic x²+y² = 1: Π(p − (−1/p))) | D for pq by the polynomial-count dichotomy §5.3 | D | yes | a | §5.3 |
| τ(N), together with τ(N²) | D, on a density-one set of RSA moduli | D (E–C Thm 15.1.1 + factorisation) | E–C, poly(log p) | a | Bach–Charles arXiv 0708.1192, Thm 2.2; `reduce_tau_pair`. Fails only if both p, q are non-ordinary (p ∣ τ(p)), e.g. N = 35 |
| a_f(N), a_f(N²), fixed non-CM eigenform | same | D given the factored index (Bruin 2011, GRH, Las Vegas) | yes | a | Bach–Charles; Bruin, Publ. Math. Besançon 2011 (via codex) |
| C(N) = Σ_{t²<4N} H(4N−t²) = Σ_{ab=N} max(a,b) = 2N + 2q | D for pq | D | C(p) = 2p | a | Kronecker–Hurwitz; Gross–Zagier 1985 (4.4) per codex; `reduce_hurwitz_sum`, brute-force verified for 7 N |
| Frobenius operator Φ_B on étale Z/N-algebras | D (Prop FROB) | D | x ↦ x^p | a (the operator behind every prime method) | §5.2; `reduce_quadratic_frobenius`, `reduce_frobenius_fixed_space` |
| h(−4N) and class-group arithmetic | R (ambiguous forms; Shanks, Schnorr–Lenstra); D under GRH (Schoof 1982 §6) | not known | **no**: h(−p) is not known in poly time (L[1/2] under GRH, Hafner–McCurley) | b | `reduce_class_number` |
| square roots in Cl(−4N) | R (two roots differ by an ambiguous class) | R, given the factors | yes, given the factorisation of the discriminant | a/b | folklore, the Rabin analogue |
| regulator of Q(√N) (infrastructure) | R (ambiguous form at the half period; Shanks) | not known | **no** | b | Shanks / SQUFOF |
| ord_N(a), random a | R; D under ERH | D, plus factoring p−1 | **no**: ord_p(a) needs the factorisation of p−1 | b | Miller 1976; Cleve 2004 (query model) |
| interval sums Σ_{a<N/2} (a/N) | gives h(−N), so as row "h(−4N)" | — | **no**: the same sum at prime p is h(−p) | b | `class_number_half_sum` (tested against enumeration) |
| Jacobi (d/N), power-residue symbols at (N) | **0 bits** beyond N mod the conductor | — | yes | c | Jager–Schwenk J. Cryptol. 2013: generically as hard as factoring, yet easy by reciprocity; de Boer–Pagano ISSAC 2017 |
| d(N), μ(N), ω(N) | not known (constant on RSA inputs) | D | — | — | Shallit–Shamir 1985; polyimpl §2.5 |

Notes:
- (i) Class (a) is closed under one rule. If #V(F_p) = P(p) is a polynomial with more than one
  monomial, then #V(Z/pq) ⇔ factoring (§5.3).
- (ii) Every class-(a) entry is a product of local factors at the unknown primes; every class-(b)
  entry is a global invariant of a field or group built from N; every class-(c) entry is abelian.

## 2. Organising principle

Everything polynomial-time computable from N is computed by some combination of three things:
- ring operations in Z/N-algebras, plus gcds;
- the integer representation through Euclid-type recursions (reciprocity, floor sums, continued
  fractions);
- lattice reduction.

GPT-6.1-sol's diagnosis (nonsmooth_1 §4, SPECULATIVE there; adopted here as a taxonomy, not a
theorem): these operations compute **aggregates that forget the CRT components**. Recovering a
component-sensitive quantity needs one of three "unpaid" operations:
- **(U1)** select a local component;
- **(U2)** construct local Frobenius or order data;
- **(U3)** evaluate a large arithmetic sum without expanding it.

How the paradigms fall under them:
- Prop FROB (§5.2) shows (U2) is factoring for Frobenius operators.
- Prop OA (§5.1) shows (U2) for orders reduces to smoothness or to the anomalous coincidence.
- The abelian aggregate theorem (§5.4) shows that reciprocity, the one known non-generic shortcut
  (Jager–Schwenk), computes only N mod m.
- (U3) is the only door left open. Every concrete (U3) sum we found is divisor detection in
  disguise: Pillai, r₄, C(N), codex's Jacobi autocorrelation A_N(H) = H − p, and the Heilbronn/
  Merel sums for T_N. The best (U3) evaluator known is N^{1/3} (Hiary, convex hull or theta sums),
  and lattice and collision methods (N^{1/4}, N^{1/5}) dominate it.

## 3. Paradigm verdicts

| Paradigm | Verdict | Precise reason |
|---|---|---|
| Edixhoven–Couveignes / Bruin at composite N | **REFUTED for this model** (it needs Φ) | a(N) = a(p)a(q) is a product of traces at two Frobenius classes, not the trace of one element. The only p-dependent step of E–C, identifying Frob_p in Gal(K_ℓ/Q), consumes Φ on an étale algebra, which factors N (Prop FROB(b, c): K_ℓ ⊃ Q(ζ_ℓ), and the sign of Frob on ψ_ℓ-roots is (p/ℓ)). Bruin requires the index in factored form (KNOWN) |
| a(N) mod ℓ, composite N | no route | For an exceptional ℓ (691, 2^11, …) it is σ_k(N) mod ℓ, a polynomial in s mod ℓ, so it contains the QR bit of §5.5. For a non-exceptional ℓ it is Chebotarev data of a pair of classes. Enough ℓ reconstruct a(N), which is Bach–Charles |
| Eichler–Shimura over Z/N | **REFUTED for this model** | T_N = T_pT_q. Mod p, T_p ≡ F + V uses Frob_p (Prop FROB); T_q is a correspondence of degree q+1 ≈ √N |
| Dokchitser Frobenius resolvents | **REFUTED for this model** | They need x ↦ x^p. x ↦ x^N is not a ring endomorphism and does not permute the roots (codex). Prop FROB(b) |
| Schoof/SEA "pretending N is prime" | **REFUTED for this model** | The Frobenius it needs on (Z/N)[x]/ψ_ℓ factors N for every ℓ with (N/ℓ) = −1 (Prop FROB(b) + (disc ψ_ℓ / p) = (p/ℓ), verified). Running it with x^N leaks only at the primality-test rate (W-13/W-14) |
| Kedlaya / Lauder / Satoh (p-adic) | **REFUTED for this model** | A Frobenius lift mod N^k reduces to Φ mod N (Prop FROB). The deformation base point needs Gauss sums over F_p (characters of F_p^*), and Gross–Koblitz needs Γ_p, whose "N-adic" analogue is the factorial unit predicate (Track C) |
| Harvey average-polynomial-time batches | **REFUTED for this model** | The amortisation is over all primes ≤ X with a shared product tree. Reaching an unknown p needs X ≥ p, i.e. trial-division cost. One composite N gives a product of length N (Track C) |
| Trace formulas (Eichler–Selberg, Kronecker–Hurwitz) | **REFUTED as a mechanism** (§4.1) | The geometric side is √N class numbers; the spectral side is Hecke eigenvalues at N plus divisor sums. In weight 2, level 1 there is only the divisor sum, and C(N) = 2N + 2q. Charles 2006 computes τ(N) from the trace formula in N^{1/2+ε} randomized (codex), which is exponential |
| Dequantising Shor | **no route** | The classical order-finding query bound is Ω(r^{1/3}) (Cleve, Inf. Comput. 192, 2004), improved to Ω̃(r^{1/2}) for bounded jumps (Chakraborty–García-Soriano–Matsliah, ICS 2011). Tang-style dequantisation needs low rank; the periodic state has none. The only non-black-box structures of a ↦ a^x known to be exploited are smoothness (p−1) and the representation (index calculus). Elements with cheap orders (±1, roots of unity in public algebras) have the same order on every component (codex), and a nontrivial √1 already is a factor |
| Lattices / Coppersmith | **no invariant supplies the hint** | The need is ≈ n/4 coset-shaped bits (Track E E1). Per-modulus data carry 0 bits (F1). The class-(a) invariants mod a small ℓ give at best s mod ℓ (§5.5); s mod ℓ for all ℓ ≤ 0.35n determines s exactly by CRT (polyimpl §2.5), but no single such residue is computable |
| Varieties with global count formulas | **REFUTED** (proved dichotomy §5.3) | Polynomial-count V: either #V(Z/N) = c²N^r (no information) or #V(Z/N) ⇔ factoring. Non-polynomial counts (curves) need Frobenius (FROB) |
| Generic ring / black-box ring | equivalence, not floor | Aggarwal–Maurer 2009; Jager–Schwenk 2013. Reciprocity is the canonical non-generic shortcut and is abelian (§5.4) |
| Class groups / infrastructure | smoothness-based or exponential | Hafner–McCurley, Schnorr–Lenstra: smooth ideals or orders, L[1/2]. SQUFOF: N^{1/4}. h(−p) is class (b) |
| Analytic / lattice-point counting | exponential | Hiary 2016 N^{1/3} (Math. Comp., arXiv 1408.2608): interval divisibility by exponential sums. Polynomial-Taylor plus floor/Weyl sums of degree j would give N^{1/(2j+2)} *if* exact degree-j floor sums were polylog. Complete quadratic floor sums mod p compute h(−p) (class b), and Coppersmith/HH already do better in the same architecture |
| CM / anomalous curves | **real but thin** (§4.2) | Prop OA; E1/E2 |
| Reciprocity amplification (codex: incomplete Jacobi autocorrelation) | **REFUTED** | A_N(H) = H − p: an interval divisor selector (codex nonsmooth_1 §3, proved) |
| Cartier/automatic compression of Δ mod ℓ | **REFUTED in the uniform form** | Mod 2 the depth-k sections have rank ≥ 2^{k−3} (codex, proved) |

## 4. Candidates developed

### 4.1 The Hurwitz aggregate (codex's sole survivor): killed

**Sketch.** Write C(N) = (1/π) Σ_t √(4N − t²) · L(1, Δ_t), with Δ_t = t² − 4N. Compress the t-sum
geometrically (theta, contour or Shintani methods) with certified error E. Since C(N) = 2N + 2q,
an error E ≤ 2N^{1/4} gives q to within N^{1/4}, and Coppersmith (Track E E1) finishes. The claimed
payoff would be L[1/4] or polynomial time, whatever the compression costs.

**Cost analysis and kill (Prop EIS, §5.6; proved for the representation, HEURISTIC beyond it).**
Expand each L(1, Δ_t) as an Euler product.
- For a prime ℓ ∤ 4N, the average of the local factor over t mod ℓ is an N-independent constant:
  the complete sum Σ_{t mod ℓ} ((t² − 4N)/ℓ) = −1.
- For ℓ = p, Δ_t ≡ t² (mod p) is a nonzero square for every t ≢ 0, so the local factor is
  (1 − 1/p)^{−1} for a 1 − 1/p fraction of t: the complete sum is p − 1. This bias, applied to the
  main term ≈ 2N, is the 2q.

So every approximation built from local factors at primes ℓ ≤ X < p (the singular series)
equals 2N·const and has error ≍ 2q ≍ √N. Capturing the 2q needs the local factor at the unknown p,
or a sum over all ℓ up to p, i.e. trial-division cost. Codex's own kill test asked whether the
certificate retains the 2q or absorbs it into an O(√N) error. The answer: it absorbs it, unless
the factor at p is known.

**No third representation.** M₂(SL₂Z) has no cusp forms. The Kronecker–Hurwitz identity is
therefore the full Eichler–Selberg trace formula at weight 2: its spectral side is the divisor sum
alone, and C(N) is the N-th coefficient of a (mock/quasi-modular) Eisenstein object.

Rademacher-type exact series need Kloosterman sums S(m, N; c). For prime c these are computable,
but their N-dependence enters only through N mod c, which is the information floor (F1). The
2q term again comes only from c divisible by p. This is HEURISTIC; codex nonsmooth_2 was asked to
break it (§8).

**Verdict.** The candidate is "approximate σ(N) to within N^{1/4}", which is the Coppersmith hint
itself. It secretly needs the factorisation: specifically, the Euler factor at p.

### 4.2 Anomalous CM (the only smoothness-free order-annihilation branch)

**Algorithm (KNOWN: Cheng 2002, ePrint 2002/109; Shirase 2017; Aikawa–Nuida–Shirase,
IEICE E102-A (2019); Sedlacek et al., SECRYPT 2019).** Let −D be a CM discriminant and E_D a CM
curve.
- If p | N with 4p − 1 = D s², then t_p = ±1, so E_D mod p or its twist has exactly p points.
- [N]P on a random x-only point vanishes mod p with probability 1/2 (1/6 for D = 3).
- For h(−D) > 1, Cheng works over (Z/N)[j]/H_D(j), at cost poly(h(−D), log N).

`src/nonsmooth/cheng.py` implements h = 1 (D ∈ {3, 11, 19, 43, 67, 163}), operationally. D = 7
is impossible: 4p − 1 ≡ 3 (mod 8) forces D ≡ 3 (mod 8), proved and confirmed in all 80 000 census
primes.

**E1 (pre-registered; PASS).**

| arm | n | success | structural | mean tries at the right D | mul / try / bit |
|---|---|---|---|---|---|
| special | 128, 256, 512 | 60/60 each | 60/60 | 2.0–2.5 (D > 3), 6.3 (D = 3) | 23.6–23.9 |
| control | 128, 256, 512 | 0/60 each | — | — | — |

- Predictions P1–P3 all held: 180/180 special; 0/180 control; 24 multiplications per ladder bit.
- At 512 bits a success costs 0.9 s in Python.

**E2 (pre-registered; PASS).** 20 000 random primes per size; D = squarefree part of 4p − 1.

| bits | median log D / log p | F(0.5) | F(0.75) | F(0.9) |
|---|---|---|---|---|
| 24 | 1.084 | 4.5e−3 | 3.14e−2 | 0.122 |
| 32 | 1.063 | 1.2e−3 | 1.68e−2 | 0.081 |
| 40 | 1.050 | 5e−5 | 8.9e−3 | 0.064 |
| 48 | 1.042 | 5e−5 | 3.7e−3 | 0.049 |

Checks against the pre-registered predictions:
- Fitted slopes of log₂F per bit: θ = 0.75: −0.127 (predicted −0.125, window [−0.19, −0.06]);
  θ = 0.9: −0.053 (predicted −0.05); θ = 0.5: −0.30 (predicted −0.25; reported only, as counts
  are 1–89).
- Prefactor c₀.₇₅ = 0.23–0.27, inside [0.02, 1].
- Q1–Q4 all hold.

So F(θ) ≈ 0.25·p^{−(1−θ)/2} is **experimentally supported**. It matches the proved upper bound
O(√(BX)) (codex nonsmooth_1).

**Cost.**
- Covering a fraction F of primes needs D up to Y ≈ 16F²p.
- Trying all D ≤ Y costs ≥ Σ_{D≤Y} h(−D) ≍ Y^{3/2} (Siegel), i.e. ≳ F³p^{3/2}.
- Balanced N gives ≳ F³N^{3/4}.
- Even a free class polynomial leaves cost Y ≍ F²p.

**Verdict: a correct special-purpose family, density zero, exponential for constant coverage.**
The non-CM variant ("anomalous ECM": random curves, a_p = 1) covers ≈ K/(√p log p) of primes with
K curves under Lang–Trotter (HEURISTIC), i.e. cost ≍ N^{1/4}: no better than ρ.

### 4.3 Smaller candidates, killed on the spot
- **Hecke action on modular symbols.** T_N acts by Heilbronn–Merel matrices, whose sum
  Σ_{x mod N} (continued-fraction product of x/N) has N terms. This is a (U3) sum whose fast
  evaluation would compute τ(N). Sums of Euclid statistics over x mod N encode divisor data
  (Heilbronn's formula; Pillai). No sub-N evaluator is known, and we found no structure.
- **Jacobi interval sums.** These compute h(−N), so they are class (b): at a prime they would give
  h(−p).
- **Ideals between (N) and O_K.** Naming one that is a nontrivial proper product of primes above
  N needs a root of f mod p (Rabin-type). Power-residue symbols at (N) are abelian (§5.4).

## 5. Barriers (statements and proofs)

### 5.1 Prop OA (order-annihilation dichotomy). PROVED (routine)

*Setting.*
- G is a smooth commutative group scheme of finite type over Z[1/M₀] with connected fibres.
- A method computes [m]x for x ∈ G(Z/N) and an integer m it can write down, and detects
  [m]x ≡ 0 (mod p) by a gcd.
- x mod p is uniform on G(F_p), or on a union of twists, each treated separately.
- Write #G(F_p) = p^e·m_p with p ∤ m_p.

(a) #G(F_p) = #T(F_p)·p^{dim U}·#A(F_p), and #T(F_p) is prime to p with #T(F_p) ≥ (p−1)^{dim T}.

*Proof.* Chevalley: 0 → L → G → A → 0 with L connected linear and A an abelian variety. Over the
perfect field F_p a commutative connected L is T × U, and U is split unipotent, so
#U(F_p) = p^{dim U}. Lang: H¹(F_p, L) = 0, so #G(F_p) = #L(F_p)#A(F_p). Next,
#T(F_p) = |det(p − F | X*(T))| with F of finite order, so its eigenvalues ζ are roots of unity.
Then ∏|p − ζ| ≥ (p−1)^{dim T}, and the integer det(p − F) ≡ det(−F) = ±1 (mod p). ∎

(b) If a prime ℓ ≠ p divides m_p and ℓ ∤ m, then Pr_x[[m]x ≡ 0 (mod p)] ≤ 1/ℓ.

*Proof.* The ℓ-Sylow S ≠ 1 of G(F_p) is a direct factor. [m] is an automorphism of S, so x must
have trivial S-component, which happens with probability 1/#S. ∎

(c) Corollary. Take m = N^k·S with S composed of primes ≤ B, and suppose success has probability
> 1/B. Then m_p | S·q^k. Excluding the coincidence q | #G(F_p) (balanced N, dim G = 1: impossible
unless q − p ≤ 2√p + 1, the Fermat regime), this says m_p is B-smooth: the Pollard p−1, Williams,
Bach–Shallit, ECM regime. The smoothness-free case m_p = 1 forces:
- dim T = 0 (p ≥ 3);
- #A(F_p) a power of p, i.e. A anomalous;
- for a split, the same must fail at q.

The unipotent part is killed at both p and q by N^k, so it never splits.

*Scope.* Statement (c) covers only scalars of the shape N^k·S. For an arbitrary computable m,
success needs m_p | m. That means a "random-looking" integer of size ≈ p^{dim} must divide m,
which needs d(m) large, i.e. a smooth m (HEURISTIC). This is not a lower bound for non-group
methods.

(d) Coverage. Anomalous elliptic curves with CM discriminant D ≤ Y cover ≍ √(Y/p) of primes
(proved upper bound O(√(BX)), codex; measured in E2). Fixed non-CM curves have a_p = 1 for
≍ √X/log X primes (Lang–Trotter, conjectural). Higher dimension is rarer: P_A(1) must hit p^j
exactly (HEURISTIC; exponent left open).

### 5.2 Prop FROB (Frobenius operators over Z/N). PROVED

Let B be an étale (Z/N)-algebra given by structure constants. Let Φ_B be its componentwise
Frobenius: the unique automorphism with Φ_B ⊗ F_p = Frob_p for each p | N. At prime N this is
x ↦ x^N, which is cheap.

(a) B = (Z/N)[y]/(y² − d) with (d/N) = −1. Then Φ(y) = εy with ε = CRT((d/p), (d/q)), so
gcd(ε − 1, N) ∈ {p, q}. The least such d is found deterministically.
`reduce_quadratic_frobenius`, tested.

(b) B = (Z/N)[x]/(f), f monic and separable mod N, with (disc f / N) = −1. Then a deterministic
O(n³) elimination on Φ_B − 1 over Z/N meets a nonzero nonunit.

*Proof.* B ⊗ F_p = ∏_{i=1}^{r_p} F_{p^{d_i}}, so dim_{F_p} ker(Frob − 1) = r_p. Frob permutes the
roots with cycle type (d_i), so its sign is (−1)^{deg f − r_p}. By Stickelberger that sign is
(disc f / p). Since (disc/p) ≠ (disc/q), we get r_p ≠ r_q, and ker(Φ_B − 1) ≅ F_p^{r_p} × F_q^{r_q}
is not free. If every pivot were a unit, the matrix would be equivalent to diag(1, …, 1, 0, …, 0),
with a free kernel. So some entry met is a nonzero nonunit. ∎
`reduce_frobenius_fixed_space`, tested on Schoof's algebra (Z/N)[x]/ψ₅.

(c) B ⊇ (Z/N)[ζ_ℓ] and (N/ℓ) = −1, so p ≢ q (mod ℓ). Then Φ(ζ) − ζ^a, at a ≡ p (mod ℓ), vanishes
mod p only. Trying a = 1, …, ℓ−1 factors N.

(d) Schoof's algebra. The sign character of GL₂(F_ℓ) on (F_ℓ² − 0)/±1 is Legendre(det): checked
for ℓ = 5, 7, 11, 13, and for ℓ ≥ 5 it is forced to factor through det because SL₂ is perfect.
Hence the sign of Frob_p on the roots of ψ_ℓ is (p/ℓ), and (disc ψ_ℓ / N) = (N/ℓ). This was
checked against sympy discriminants of ψ₅ and ψ₇ for 400 primes, with 0 mismatches. So for every
ℓ with (N/ℓ) = −1 (half of all ℓ), the Frobenius operator Schoof needs on (Z/N)[x]/ψ_ℓ factors N
via (b).

*Scope.* FROB is about producing the *operator* (as matrix or images). It says nothing about
algorithms that output only a number such as a(N) mod ℓ; for those, use Bach–Charles. Lenstra 1985
("Galois theory and primality testing") and Grantham's Frobenius pseudoprimes are the nearest
literature. FROB is probably folklore.

### 5.3 Polynomial point-count dichotomy. PROVED (codex nonsmooth_1 §4; re-checked)

Let #V(F_p) = P(p) = Σ_{i=r}^{d} c_i p^i with c_r c_d ≠ 0.
- If r = d, then #V(Z/pq) = c_r²N^r, which is no information.
- If r < d, then #V(Z/pq) is a polynomial of degree d − r in s with leading coefficient
  c_r c_d N^r. Its roots give s, so the count is deterministically equivalent to factoring.

The proof writes p^j + q^j = U_j(s, N) by the Lucas recurrence. This also answers the question of
"varieties with a global formula": only the monomial counts have one.

### 5.4 Abelian aggregate theorem. PROVED (codex; Kronecker–Weber)

For K/Q abelian of conductor m and (N, m) = 1, Art_{K/Q}((N)) depends only on N mod m. With
Track E F1, characters and power-residue symbols evaluated at (N) carry 0 bits beyond N mod m.

*Scope* (codex): this is a redundancy statement, not a lower bound. tr ρ(Art(N)) ≠ a(N), even for
reducible ρ: σ_k(pq) = 1 + N^k + p^k + q^k, whereas 1 + N^k is the aggregate. "Abelian origin" does
not make a Hecke eigenvalue easy.

### 5.5 Minimal primitive (Prop RANK). PROVED (elementary)

For the reducible (Eisenstein) members of class (a), the value mod a small prime ℓ is a polynomial
in s mod ℓ given N. Examples: σ(N) mod ℓ; τ(N) mod 691 = σ₁₁(N) mod 691; φ(N) mod ℓ. Such a
polynomial has at most d candidate roots. Knowing s mod ℓ fixes {p, q} mod ℓ as the roots of
x² − sx + N. When (N/ℓ) = 1 it therefore fixes the common value (p/ℓ) = (ℓ*/p), i.e. whether ℓ* is
a square mod N. Examples: ℓ = 3 asks whether N = x² + 3y², and ℓ = 4 (s mod 8) gives ε, the
two-squares bit (Track E E2′).

Whether "fixed-base QR mod N" is easier than factoring is open; no reduction is known in either
direction (codex r6; Seurin's 2-Φ/4-hiding). So the lowest rung of every class-(a) route is an open
problem that is not known to be easier than factoring itself. Non-Eisenstein ℓ do not avoid this:
enough residues reconstruct the integer value, which is Bach–Charles.

### 5.6 Prop EIS (Hurwitz aggregate). PROVED for the stated representation

- Identity: C(N) = 2σ(N) − Σ_{d|N} min(d, N/d) = 2N + 2q for N = pq. Kronecker–Hurwitz,
  brute-force checked for N ∈ {35, 51, 77, 143, 377, …}.
- The Euler-factor bias computation of §4.1: complete sums −1 for ℓ ∤ 4N and p − 1 for ℓ = p.
- The singular-series truncation error is exactly 2q when X < p.
- "No other representation evades it" is HEURISTIC; see §8.

## 6. What would change the verdict

A (U3) evaluator that is genuinely sub-N^{1/5} for a divisor-detecting sum and does not localise at
p would change it. Examples: an exact polylog algorithm for degree-2 floor sums Σ_{h<H} ⌊P(h)⌋ with
H ≪ denominator; or a Rademacher/Kloosterman truncation for C(N) whose tail can be bounded without
c ≡ 0 (mod p). Either would be a new analytic algorithm, and nothing in the literature we checked
points to one. Also open: whether the fixed-base QR bit (§5.5) is easier than factoring.

## 7. Reproduction

```
.venv/bin/python -m pytest -q tests/test_nonsmooth.py
.venv/bin/python -m src.nonsmooth.run_experiments e1   # ~4 min
.venv/bin/python -m src.nonsmooth.run_experiments e2   # ~1 min
```

## 8. GPT-6.1-sol adversarial round (nonsmooth_2)

(pending)

## 9. Opus critics

(pending)
