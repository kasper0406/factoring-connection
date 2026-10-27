# conn_cmp: complexity, operator algebras, quantum, crypto, ML, OT, numerics, sparse recovery, games, arithmetic statistics ↔ factoring (phase 28)

Agent `conn_cmp`, 2026-10-27. Brief: [`AGENT_PHASE28_SURVEY.md`](AGENT_PHASE28_SURVEY.md). This is a survey. I ran two tiny
sanity scripts (`scratchpad/conn_cmp/{phi2_rabin,aff_bc}.py`, each under a second). The GPT-6.1-sol blind round is in
`codex/conn_cmp_1.{md,out}`; its candidates are merged into the tables and marked **[cx]**. The Opus read-only critic's
review is in §12.

**Verdict labels** (from the brief):
- **RESTATEMENT**: needs a factor to write down.
- **AGGREGATE**: depends only on N mod m.
- **FACT-EQ**: factoring-equivalent (or DL-equivalent).
- **→M / →A / →Q / →⟨barrier⟩**: reduces to one of the three leak families of [`leakmap`](leakmap_results.md) §0 (M = group
  order / smoothness, A = archimedean, Q = quadratic genus), or to a named project barrier.
- **KNOWN**: already used in factoring.
- **OPEN**: worth a look.

★ marks something surprising or new to the project. Claim labels: PROVED (routine, proof here), KNOWN (with a
reference), HEURISTIC, UNVERIFIED.

---

## 0. Verdict (plain language)

**No field gives a new mechanism, and nothing here survives as OPEN at the asymptotic level** (below L[1/3], below
rigorous L[2/5], or below deterministic N^{1/5}). In 60 candidates over ten fields (7 of them contributed only by codex), every
asymptotic route fails at one of six points:
1. **The object needs the factors.** Examples: Rédei symbols, ambiguous forms, zero-encodings in multilinear maps, the
   Galois-averaged Bost–Connes state values, and the Plancherel decomposition of Aff(Z/N).
2. **The object is reachable only through ring arithmetic mod N**, so Theorem PT applies. Examples: PIT, crossed
   products, planted isogeny walks, generic delay functions.
3. **The object is a group order or a relation (M).** Examples: Tal's sub-√ factorial (it needs a divisor of p − 1),
   the tensor-network cost of Shor, which depends on the order r and its factors, and the classical side of Regev's
   lattice, which is a doubly-smooth congruence.
4. **The object is archimedean or an exact analytic value (A, barrier F).** Examples: θ-series, L(1, χ), CM values, and
   sparse-FFT recovery of 1[gcd(x, N) > 1].
5. **It is a distributional statement with no per-N algorithm.** Examples: Cohen–Lenstra, Malle, Haddad–Koukoulopoulos
   couplings, Wasserstein bounds, and ML batteries.
6. **It is a constant-factor or polylog lever only.** Examples: FFT and product trees, learned polynomial selection,
   cost models, and sub-quadratic GF(2) linear algebra.

What is left is **constant-level and already known to the project**:
- **Sub-quadratic sparse GF(2) dependency finding.** This changes only the MNFS constant (GNFS stays at 1.9230; MNFS
  falls to 1.8595 for ω ≲ 1.667); see §7 and the attribution note below.
- **Codex's three OPEN items.** On inspection each is a known project problem: Lee–Venkatesan fruitfulness, LA
  structure, and output-sensitive smooth-pair enumeration (§11).

**Small items new to the project** (none changes a bound):
- ★ **Isogeny neighbours over Z/N (§4.1).** This is KNOWN in part: Altuğ–Chen (ePrint 2018/926, Thm 3.3) show that
  finding *all* ℓ-neighbours factors N, and so does finding *one* if the solver returns a random root. The small new
  step here removes the random-output assumption for ℓ with genus-0 X_0(ℓ) (ℓ = 2, 3, 5, 7, 13). Planting the root
  through the Fricke parametrisation makes *any* root-finder for Φ_ℓ(j, X) mod N a randomized factoring algorithm.
  For ℓ = 2 the toy check gives 352/795 = 0.443 splits against a predicted 4/9. So a CSIDH/CRS-style walk over Z/N
  can be navigated only along planted (ring-arithmetic) parametrisations, which are PT.
- ★ **Bost–Connes (§2.1).** Individual KMS_β values Li_β(e(a/N))/ζ(β) are analytic functions of the real number a/N.
  Their Galois sum is Σ_{d|N} μ(N/d)d^{1−β} = (p^{1−β} − 1)(q^{1−β} − 1) (checked numerically), which is FACT-EQ. The
  phase transition sits exactly on the line between "computable" and "factor-revealing".
- ★ **Aff(Z/N) for squarefree N (§2.5).** The group has exactly N conjugacy classes (checked for N ≤ 35), so
  rank K_0(C*(Aff(Z/N))) = N carries nothing. Its irreducible dimensions {1, p−1, q−1, φ(N)} are FACT-EQ.
- ★ **Regev's lattice (§3.5).** Its classical side is a relation-lattice target of the E3(ii) type in
  [`suffmap`](suffmap_results.md), with the norm relaxed to 2^{O(√n)} in ≈ √n dimensions. Lattice vectors are
  congruences u ≡ v (mod N) between products of tiny primes. No classical sub-L[1/3] method for this relaxed target is
  known.
- ★ **Learned smoothness filters (§5.3)** cannot move the MNFS constant at ω = 2, because the test term has slack
  ([`heurnfs`](heurnfs_results.md)). Fixed-depth transformers lie in TC⁰, where [`leakmap`](leakmap_results.md)'s AC⁰
  theorem does not reach. TC⁰ contains factoring-based PRFs (Naor–Reingold), so "learn factor bits" is a FACT-EQ
  target.
- ★ **Attribution note for [`heurnfs`](heurnfs_results.md) row 20** (found by codex; checked by the Opus critic, who
  re-ran heurnfs's model). heurnfs's MNFS-vs-LA-exponent curve is **KNOWN for ω ∈ [1.78, 2]**:
  Lenstra–Shamir–Tomlinson–Tromer 2002 (§2.7) analyse MNFS with matrix cost L(β)^{2ε} and get 1.9019, 1.8870, 1.8719
  and 1.8689 at 2ε = (3+√17)/4 ≈ 1.7808. LSTT call 1.8689 "optimal relation collection". It is a point on their
  constrained curve (sieve = ECM cost), not the global minimum. In the same model the cost keeps falling below 1.78:
  1.8647 at 1.75, 1.8604 at 1.7, and flat at **1.85954 for ω ≲ 1.667**. So heurnfs is right, and its part below 1.78 is
  a correction to LSTT's wording. Its novelty claim should be narrowed to that part. I did not edit heurnfs.

Subjective: P(any item here leads to an exponent below 1/3, or a rigorous one below 2/5) < 0.5%.

---

## 1. Algebraic complexity

| # | Connection (object from N; feature or primitive) | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 1.1 | Short SLPs for t! or the Pochhammer polynomial ∏(x−i): gcd(⌊√N⌋! mod N, N) = p. τ-conjecture / VP⁰ vs VNP⁰ | suffmap tier E; slpcirc Thm SP | KNOWN (tier E) | If the permanent has poly-size constant-free circuits, then τ(n!) = polylog and factoring ∈ P/poly. A VNP *collapse* would help; GCT and lower bounds cannot. Believed false | Lipton, ANTS-I (1994); Shub–Smale, Duke 81 (1995); Bürgisser, Comput. Complex. 18 (2009) |
| 1.2 | Tal's sub-√ modular factorial (ECCC TR26-211), Õ(q^c + √(p/q)) for q \| p−1 (and cyclotomic variants); quantum version arXiv 2607.29453 | codex_r1_tal_bgs (full read) | RESTATEMENT + →M | The method needs the prime p, a finite-field character of order q, and a divisor q of p−1 (or of Φ_r(p)), i.e. group-order structure at the unknown prime. Even a componentwise saving of t^α beats N^{1/5} only if α < 2/5 **[cx]** | Tal 2026; BGS, SICOMP 36 (2007) |
| 1.3 | Pollard–Strassen / BGS / Harvey and Harvey–Hittmeir via polynomial multipoint evaluation | many notes | KNOWN | These give the deterministic records (N^{1/4}, then N^{1/5}) | Harvey 2021; HH 2022 |
| 1.4 | PIT / ring morphisms: polynomials over Z/N whose identity behaviour differs mod p and mod q. Finding a nontrivial ring automorphism or isomorphism | factorbits PT | FACT-EQ / →PT | Kayal–Saxena: finding a nontrivial automorphism of a ring given in basis form is randomized-equivalent to factoring. PIT only *tests* a supplied polynomial, and over Z/N it uses N as a modulus only, so Theorem PT applies **[cx]** | Kayal–Saxena, Comput. Complex. 15 (2006); Agrawal–Biswas, JACM 50 (2003) |
| 1.5 | Hitting sets / PIT derandomisation for deterministic factoring | zroute DSR; dsr Thm SB | →DSR barrier | Every black-box derandomiser fails on DSR (dsr Thm SB). Deterministic PIT would not supply DSR's constructive lemma | dsr_results |
| 1.6 | Tensor rank / bilinear complexity of multiplication (in Z, Z/N, F_{p^k}); fast integer multiplication; ω for dense matrices | suffmap Thm SS(c) | KNOWN (polylog) | Multiplication Z/N × Z/N → Z/N has bilinear rank 1 over Z/N **[cx]**. Splitting it into CRT components needs the idempotents. Gains are polylog and absorbed into L^{o(1)}; single-family sampling ignores ω ≤ 2 | Harvey–van der Hoeven, Ann. Math. 193 (2021) |
| 1.7 | **[cx]** Witness-counting polynomials (VNP form) Σ_{x,y∈{0,1}^n}[xy = N]·(prefix of x): fast evaluation would recover the factor by prefix queries | — | FACT-EQ | This is #P-style counting of factor witnesses. GCT gives a lower-bound framework, not uniform constructions | Mulmuley–Sohoni, SICOMP 31 (2001) |
| 1.8 | **[cx]** ★ "Ultimate complexity" of a multiple m_t·t! (Cheng): short SLPs for some nonzero multiple of t! | — | →multiplier contamination | m_t may contain the other prime, so the gcd becomes N. Weaker compression ≠ a modular-factorial algorithm | Cheng, TCS 2004 (UNVERIFIED details) |
| 1.9 | Generic-ring lower bounds (algebraic complexity of factoring itself) | factorbits PT; leakmap | KNOWN | Breaking RSA, or speeding up repeated squaring generically, is equivalent to factoring | Aggarwal–Maurer, EUROCRYPT 2009; Rotem–Segev, CRYPTO 2020 |

## 2. Operator algebras / noncommutative geometry

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 2.1 | ★ Bost–Connes KMS_β states (β > 1) evaluated on e(a/N) | newfamily row 48 (polylogs only) | individual values: no factor structure; Galois sum: FACT-EQ | φ_β(e(a/N)) = Li_β(e(a/N))/ζ(β) is an analytic function of the real a/N, computable to any precision but blind to the arithmetic of N. The Q-rational (Galois-summed) value is Σ_{d\|N}μ(N/d)d^{1−β} = (p^{1−β}−1)(q^{1−β}−1) (PROVED via Ramanujan sums; checked at N = 33, 35 to 10⁻¹²). It gives p^{1−β} + q^{1−β}, but it is a φ(N)-term sum, or needs the factors. For β ≤ 1 the unique KMS state takes exactly the Galois-averaged value N^{−β}∏_{p\|N}(1−p^{β−1})/(1−p^{−1}) on e(a/N) (critic), which is FACT-EQ. So the phase transition separates "computable but blind" from "factor-revealing" | Bost–Connes, Selecta 1 (1995); Connes–Marcolli book (2008) |
| 2.2 | **[cx]** Bost–Connes Gibbs sampling: draw k with Pr ∝ k^{−β}, test gcd(k, N) | — | →A (trial division) | Pr[p \| k] + Pr[q \| k] ≈ p^{−β} + q^{−β}: a rare event. The symmetry breaking at β = 1 gives no conditioning on the unknown primes | Connes–Consani |
| 2.3 | Connes' adelic trace formula / Weil positivity / RH programme; primes as orbits of length log p | geometry row 50 | RESTATEMENT / →barrier F | The objects are N-independent. The explicit formula reproduces sieving, and isolating one unknown orbit is a separate inversion problem. RH would give no algorithm | Connes, Selecta 5 (1999) |
| 2.4 | Connes–Marcolli–Ramachandran imaginary-quadratic system for K = Q(√−N): ground-state values generate K^ab | genusleak §4 (CM) | →A (exponential) | The values needed are singular moduli of discriminant −4N, with a class polynomial of degree h(−4N) ≈ √N. CM constructions cost \|D\|^{1/2+o(1)} | CMR, Selecta 11 (2005) |
| 2.5 | ★ C*(Aff(Z/N)), Aff = Z/N ⋊ (Z/N)^*; Cuntz's ax+b-semigroup algebra; Cuntz–Li ring C*-algebras | — | K-theory: trivial; dimensions: FACT-EQ; operations: →PT | For squarefree N, Aff(Z/N) has exactly N conjugacy classes (multiplicative; Aff(F_p) has p), so rank K_0 = N (PROVED; checked for N ≤ 35). The irreducible dimensions {1, p−1, q−1, φ(N)} factor N, but computing central idempotents needs the CRT idempotents of Z/N. **[cx]** C(Z/N) ⋊ Z/N ≅ M_N(C) has trivial centre: the noncommutative hull forgets the CRT structure | Cuntz 2008; Cuntz–Li, JNCG (2010); Li arXiv 0905.4861 |
| 2.6 | Subfactors: intermediate subfactors of R ⊂ R ⋊ Z/N | — | RESTATEMENT | The intermediate lattice is the subgroup lattice of Z/N, i.e. the divisor lattice of N | Jones 1983 |
| 2.7 | Connes–Consani arithmetic site / F_1 point counts | geometry row 50 | RESTATEMENT | Counting functions reproduce the explicit formula, and there is no N-specific object | — |

## 3. Quantum information (classical side only)

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 3.1 | ★ MPS / tensor-network simulation of Shor: bond dimension of Σ_x\|x⟩\|a^x mod N⟩ | periodfind §5 (not TN) | →M | The Schmidt rank across the exponent/output cut is the order r, with near-flat weights **[cx]**. Dang et al. find that the space cost depends on r *and its factors*. Our gloss (not their wording): structured r is cheaper, which is the p−1 (M) family resurfacing inside the simulation | Dang–Hill–Hollenberg, Quantum 3 (2019), arXiv 1712.07311; Vidal, PRL 91 (2003) |
| 3.2 | Stabiliser / Clifford+T / matchgate simulation | periodfind §5 | KNOWN, →exponential | The T-count of modular exponentiation is poly(n), and Bravyi–Gosset is exponential in T. Small matchgate extensions are universal **[cx]** | Bravyi–Gosset, PRL 116 (2016); Jozsa–Miyake (2008) |
| 3.3 | Dequantisation (Tang; CGLLTW; Gharibian–Le Gall) | periodfind §5 | →periodfind (Prop BASIS) | The exponentiation matrix has stable rank r. Querying eigenvector entries is a discrete log. Fourier sampling is the missing primitive | Tang, STOC 2019; CGLLTW, JACM 2022 |
| 3.4 | Semiclassical / approximate QFT is classically simulable | periodfind §5 | KNOWN, bottleneck elsewhere | The QFT is easy on product or low-entanglement inputs. All the difficulty is in the modular-exponentiation state (Theorem T) **[cx]** | Griffiths–Niu, PRL 76 (1996); Yoran–Short 2007; ALM 2006 |
| 3.5 | ★ Classical side of Regev's multidimensional algorithm: L = {e : ∏ b_i^{e_i} ≡ 1 (mod N)} with b_i small primes | schnorrcvp (Regev only as a lattice author) | →M (suffmap tier E) | A short e ∈ L is a congruence u ≡ v (mod N) between two products of tiny primes, i.e. a doubly-smooth relation. This is an E3(ii)-type relation-lattice target (suffmap), **with the norm relaxed to 2^{O(√n)} in ≈ √n dimensions**: u and v are ≈ exp(2^{√n}), far bigger than N, and Regev also needs vectors outside the trivial sublattice (the analogue of fruitfulness). No classical sub-L[1/3] method is known. Pilatte proved Regev's number-theoretic conjecture | Regev, JACM (2025) / arXiv 2308.06572; Pilatte, Forum Math. Pi 14, e5 (2026); Ekerå–Gärtner, IACR CiC 1 (2024) |
| 3.6 | Tensor-network Schnorr sieving (Tesoro et al., PRA 113, 032418 (2026); "polynomial scaling" seen numerically up to 100 bits) | schnorrcvp | →schnorrcvp | A TN solver for CVP in Schnorr's lattice inherits that lattice's yield, which schnorrcvp refutes quantitatively (γ → 1 needed). Numerics at ≤ 100 bits cannot separate poly from L[·]. Exact contraction of an array-multiplier network has treewidth Θ(n) | Markov–Shi, SICOMP 38 (2008) |
| 3.7 | Classical simulation of noisy Shor circuits (Pauli paths) | — | →periodfind (analogy) | AGLLV concerns *random* circuits, not Shor. By analogy, at constant noise the period peaks decay, and simulability would come exactly from losing the order information | Aharonov–Gao–Landau–Liu–Vazirani, STOC 2023 |

## 4. Cryptographic constructions as sources of structure

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 4.1 | ★ ℓ-isogeny graphs over Z/N (CSIDH/CRS style): vertices j ∈ Z/N, i.e. the tensor product G_p × G_q. Navigation = roots of Φ_ℓ(j, X) mod N | genusleak §4 (CM); tracedata | FACT-EQ (KNOWN in part) / →PT | Altuğ–Chen Thm 3.3: all neighbours factor N, and one neighbour does if the solver's output is random. PROVED here (routine; principle Rabin 1979), for genus-0 X_0(ℓ): plant a root with the Fricke parametrisation, for ℓ = 2 j_0 = (t+16)³/t, j_1 = (t+256)³/t². Given j_0, the root j_1 is uniform among the F_p-rational roots at each prime, independently by CRT. So *any* root-finder returns a CRT-mixed root with constant probability, and gcd(r − j_1, N) splits N. Under the planted law Φ₂(j_0, X) has 1 or 3 roots mod p with probability ≈ 1/2 each, so Pr[oracle matches mod p] = 2/3 and Pr[split] = 2·(2/3)(1/3) = 4/9; toy 352/795 = 0.443. For ℓ with X_0(ℓ) of positive genus a planting method is still needed. Walks along planted parametrisations are ring arithmetic (PT). Altuğ–Chen conjecture variants of neighbour search hard **[cx]**; Galbraith–Perrin–Voloch's level structure reduces to CSIDH, not to factoring | Altuğ–Chen, ePrint 2018/926 (ASIACRYPT 2019); GPV, ePrint 2023/1726 |
| 4.2 | SIDH / SQIsign structures over Z/N: E_{1728} is supersingular at both primes of a Blum N; quaternion norm equations; four-square representations of N | geometry (sphere ladder) | →M / RESTATEMENT | The group orders are p+1 and q+1 (Williams p+1, M). The quaternion algebras B_{p,∞} and B_{q,∞} differ, but no computable invariant reads them. Four-square representations of N are easy (Rabin–Shallit 1986) and carry no factor without a second, "different" one (r_4 = 8σ(N)) **[cx]**: small-CM examples return to order-based methods | Rabin–Shallit, CPAM 39 (1986) |
| 4.3 | ★ CLT13 multilinear maps: zeroizing attacks factor x_0 = ∏p_i from zero-test values (eigenvalues = CRT components) | — | RESTATEMENT + →A | The attack needs encodings of zero, which only the key holder (who knows the p_i) can make. The mechanism is smallness of CRT residues (approximate GCD), i.e. A | Cheon–Han–Lee–Ryu–Stehlé, EUROCRYPT 2015 (ePrint 2014/906) |
| 4.4 | iO / obfuscated programs hiding factors | — | RESTATEMENT | Any factor-bearing program is constructed with the witness. Attacks exploit planted algebra | — |
| 4.5 | Hidden-order groups: RSA groups (delay functions, strong RSA), class groups of unknown order | genusleak Thm GB; leakmap | KNOWN, FACT-EQ / →M | A generic speed-up of repeated squaring is equivalent to factoring. A multiple of λ(N) factors. Class-group order reduces to relations. Codex warns: RSA inversion, strong RSA and class-group DL are not all equivalent **[cx]** | Rotem–Segev, CRYPTO 2020; Katz–Loss–Xu, TCC 2020; Damgård–Koprowski 2002 |

## 5. Machine learning

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 5.1 | ★ Learning factor bits / ε / QR from N | eps (66 non-CF stats null at \|ρ\| ≥ 0.022); newfamily (84); classE (130); leakmap AC⁰ theorem | AGGREGATE (empirical) / FACT-EQ (as a target) | Trained transformers on arithmetic learn divisibility by small primes and by divisors of the base: Charton's GCD models output the largest such divisor, an AGGREGATE rule. Log-precision transformers lie in uniform TC⁰. leakmap's AC⁰ theorem does not extend there, but TC⁰ contains Naor–Reingold PRFs (secure if factoring is hard), so learning TC⁰ is factoring-hard (Kearns–Valiant, Kharitonov). A random-parity predictor would factor by Goldreich–Levin **[cx]**; one fixed bit has no such reduction (factorbits) | Charton, ICLR 2024; Merrill–Sabharwal, TACL 2023; Naor–Reingold, FOCS 1997; Kearns–Valiant, JACM 41 (1994) |
| 5.2 | Learned polynomial selection / sieve parameters (Murphy E, root properties) | heurnfs | KNOWN, constant only | A fixed multiplicative gain in yield is L^{o(1)}. Murphy α captures only local densities **[cx]** | Murphy thesis (1999) |
| 5.3 | ★ Learned smoothness filter before cofactorisation | heurnfs §4.4 | →heurnfs (null at ω = 2) | It changes only the test term, which has slack (1.777 vs 1.902). It matters only if it raises yield, and then it is a sampler (suffmap tier C, worth 1.6386 if achieved; smoothenum found no mechanism that achieves it) or only after ω < 2 | heurnfs; suffmap |
| 5.4 | Generative models for smooth numbers near kN | sihard Thm J; siconstruct | →sihard | Every oblivious constructor costs L[2/3] per output. A neural sampler is such a constructor | sihard |
| 5.5 | **[cx]** Length-generalising arithmetic transformers (scratchpads) | — | →uniformity barrier | Accuracy at bounded size gives no growing-n success bound and no total-cost bound | ICLR 2025 arithmetic transformers |
| 5.6 | ML-driven algorithm search (FunSearch, AlphaEvolve) | — | meta, no mechanism | A search method; whatever it finds faces the same barriers | Romera-Paredes et al., Nature 2024 |

## 6. Optimal transport

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 6.1 | Wasserstein control of relation supply | multidescent/pc3/c12* (TV and covariance chain) | →relation-supply chain | 1_{y-smooth}(n) is not Lipschitz in any translation-invariant metric (n and n+1 are independent events). Smoothness has probability L[1/3, −c], and an O(1) transport error gives no relative error **[cx]**. The chain needs TV/Fourier (BV, large sieve) control, which is where AB is stuck | — |
| 6.2 | **[cx]** ★ Coupling of prime-factor profiles with Poisson–Dirichlet (expected log-discrepancy O(1)) | — | distributional only (FACT-EQ if made conditional on N) | It couples laws, not the given integer. A per-N sampler for the factor profile would be factoring | Haddad–Koukoulopoulos, J. Éc. polytech. 12 (2025) |
| 6.3 | Positive coupling of rational and algebraic NFS norms (raise joint smoothness) | heurnfs lever (a); smoothenum | KNOWN constant / →smoothenum | Beyond shared local densities (Murphy α, Coppersmith's rational-side reuse = MNFS), the events are heuristically independent. A *constructive* L[1/3, η] gain is heurnfs lever (a) / smoothenum's open enumeration problem **[cx] OPEN item 3** | Coppersmith, J. Cryptology 6 (1993) |
| 6.4 | Couplings in proofs: fruitfulness by fibre symmetry; Chen–Stein Poisson approximation for child boxes | README key idea 1; pc3 Prop DC | KNOWN / →AB | The project's fruitfulness argument already *is* a measure-preserving involution coupling. Chen–Stein's error term is the pair-dependence term, i.e. exactly C₁/C₂ (proved rationally by c12rat; algebraic block AB open). No new input | Arratia–Goldstein–Gordon, Ann. Probab. 17 (1989) |

## 7. Numerical analysis

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 7.1 | Sub-n² sparse GF(2) dependency finding | heurnfs §4.4(b), suffmap C5, codex r21_linalg_scope | OPEN (constant only; known open problem) | GNFS stays at (64/9)^{1/3} for every ω' ≤ 2 (Thm SS(c)). MNFS falls: 1.9019 at ω' = 2, 1.8870 at 1.9, 1.8719 at 1.8, 1.8689 at (3+√17)/4 (LSTT §2.7; codex re-derived this), then 1.8604 at 1.7 and flat at 1.8595 for ω' ≲ 1.667 (heurnfs; confirmed by the critic, §0). The state of the art is ≈ n² (Wiedemann/Lanczos); full-nullspace methods (Eberly et al. 2007, Casacuberta–Kyng 2022) are superquadratic. heurnfs found that block-angular MNFS structure gives nothing at ω = 2 | LSTT 2002, cs-people.bu.edu/tromer/papers/meshc; Wiedemann 1986; Coppersmith 1994 |
| 7.2 | Hierarchical matrices / nested dissection on the relation matrix | r21 | →no structure | These need small separators or low-rank off-diagonal blocks. The relation hypergraph is random-like (Bender–Canfield), so separators are Θ(n) and GF(2) blocks have full rank. Alon–Yuster planar solvers do not apply | Alon–Yuster, FOCS 2010; Bender–Canfield 1999 |
| 7.3 | FFT / FMM for sieving; product/remainder trees for batch smoothness | many notes | KNOWN (polylog) | Sieving is already O(log log B) per candidate. Divisibility kernels have no smooth far field to compress **[cx]**. Bernstein's trees cost polylog per input and do not create joint smoothness | Bernstein 2004 (smoothparts); Greengard–Rokhlin 1987 |
| 7.4 | Real or rational sparse solvers on lifted exponent matrices (Q-kernel ⇒ GF(2) kernel) | r21 | →exactness | Peng–Vempala and Nie give n^{2.27–2.33} > n² and need conditioning. An approximate real null vector encodes no mod-2 dependency **[cx]** | Peng–Vempala, SODA 2021; Nie, STOC 2022 |
| 7.5 | High-precision numerics of analytic invariants (θ-series r_4 = 8σ(N), L(1, χ_{−4N}) ⇒ h) | analytic barrier F; genusleak | →A / barrier F | Exact values are needed. Euler-product or analytic routes cost \|D\|^{1/5} (genusleak) | Shanks; Srinivasan 1998 |
| 7.6 | Cost-model changes (mesh or circuit linear algebra) | heurnfs | KNOWN | These change the AT-cost constant, not operation counts | Bernstein 2001; LSTT 2002 |

## 8. Compressed sensing / sparse recovery / inverse problems

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 8.1 | Group testing for smoothness | — | KNOWN | This is Bernstein batch smoothness (product-tree gcds), and the sieve itself pools by arithmetic progressions. Unit-cost tests hide the product's bit length **[cx]** | Bernstein 2004 |
| 8.2 | Sparse exponent-vector recovery from log y = Σ e_ℓ log ℓ (one real measurement) | schnorrcvp | →schnorrcvp | This is a knapsack/CVP in Schnorr's lattice, quantitatively refuted for factoring | Schnorr 2021 |
| 8.3 | ★ Sparse FFT on hidden-period signals; inverse problem "recover p from f(x) = 1[gcd(x, N) > 1]" | periodfind (sft.py) | →A (trial division) / →periodfind | f = 1[p\|x] + 1[q\|x] − 1[x = 0] is (p+q−1)-sparse in frequency (every multiple of 1/p and 1/q, flat weights; critic), so sparse FFT already needs Ω(k) ≈ √N samples. Its support has density ≈ 2N^{−1/2}, and any nonzero sample already factors. Prony/Berlekamp–Massey would recover the linear complexity p+q−1 (hence p+q) but needs ≈ 2(p+q) samples: again √N. Observables of a^x need a large-order character (Theorem T) **[cx]** | Hassanieh–Indyk–Katabi–Price, STOC 2012 |
| 8.4 | Sparse (low-weight) dependencies | — | irrelevant / NP-hard in general | NFS dependencies need not be sparse. Minimum-weight codeword is NP-hard **[cx]** | Berlekamp–McEliece–van Tilborg 1978 |

## 9. Game theory / interactive proofs

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 9.1 | Decision factoring ∈ NP ∩ coNP (and UP ∩ coUP) | — | KNOWN | It is not NP-hard unless NP = coNP. A prover's certificate is the factor | Pratt 1975 |
| 9.2 | ★ TFNP: factoring ≤ PPA and ≤ WeakPigeon (PPP) by *randomized* reductions; in PPA ∩ PPP under GRH (derandomised) | — | KNOWN | Membership in PPAD is open, so equilibrium solvers do not factor. Downward self-reducibility would place it in UEOPL (selfred: none in algebraic classes) **[cx]** | Jeřábek, JCSS 82 (2016); Buresh-Oppenheim 2006; Harsha–Mitropolsky–Rosen, ITCS 2023 |
| 9.3 | ★ ZK / IP certificates of N's structure (Blum, product of two quasi-safe primes) | factorbits (ZK mention) | KNOWN (a no-leak theorem in itself) | The ZK simulator computes the whole transcript from N alone, so the certificate carries no factor information. Extraction needs a successful prover **[cx]** | van de Graaf–Peralta, CRYPTO 1987; Gennaro–Micciancio–Rabin, CCS 1998; Camenisch–Michels, EUROCRYPT 1999 |
| 9.4 | **[cx]** Sum-check / PCP for factor-witness counts | — | FACT-EQ | A cheap verifier does not bound prover cost. Replacing the prover is the hard computation | LFKN, JACM 39 (1992) |
| 9.5 | Delegation (GKR) and mechanism design for distributed NFS | — | constant | These handle verification and incentives only | — |

## 10. Arithmetic statistics

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| 10.1 | Cohen–Lenstra / genus theory for Cl(−4N): 2-rank = ω(N) − 1; 4-rank via the Rédei matrix; 8-rank via Rédei symbols / governing fields | genusleak §2.5; newfamily row 18; codex_r4 | →Q (4-rank) / RESTATEMENT (Rédei) / →M | The 2-rank is an aggregate count. The 4-rank is built from (p/q)-type symbols. The 8-rank uses non-abelian (Rédei) Frobenius data, but evaluating it needs the factors, or the class group (relations) | Fouvry–Klüners, Invent. 167 (2007); Stevenhagen, MPCPS 172 (2022); Smith, arXiv 1702.02325 |
| 10.2 | Negative Pell x² − Ny² = −1 (sign of N(ε)), governed by 4- and 8-ranks | genusleak Prop HR | →M / →A | The continued-fraction period parity is exponential, and infrastructure (Buchmann) needs relations. When N(ε) = +1, the period midpoint already factors N (72%, k = 1) | Chan–Koymans–Milovic–Pagano, "The 8-rank of the narrow class group and the negative Pell equation", Forum Math. Sigma 10, e46 (2022) |
| 10.3 | Bhargava parametrisations: 2×2×2 cubes (Gauss composition); binary cubic forms ↔ 3-torsion of Cl(−4N) (Davenport–Heilbronn); pairs of ternary forms (quartic rings) | — | RESTATEMENT / →GB | 2-torsion classes are the factorisations of N, so constraining a cube to output one presupposes them. 3-torsion gives neither an ambiguous class nor a CRT separator. Random classes are already easy to sample, and the quartic-resolvent link concerns *cubic* class groups **[cx]** | Bhargava, Ann. Math. 159 (2004) I, II; Bhargava–Shankar–Tsimerman (2013) |
| 10.4 | Smoothness of h(−kN) (Cohen–Lenstra-type heuristics) | KNOWN in the project | KNOWN, M | SPAR (Schnorr–Lenstra); Lenstra–Pomerance class-group factoring | Schnorr–Lenstra, Math. Comp. 43 (1984); LP 1992 |
| 10.5 | 2-Selmer of the congruent-number curve y² = x³ − N²x (Monsky matrices; Smith) | — | →Q / RESTATEMENT | The data are Legendre symbols among the prime factors, and 2-descent needs the primes dividing N | Heath-Brown, Invent. 118 (1994); Smith 2017 |
| 10.6 | Malle counts / prescribed ramification **[cx]**; smooth values of polynomials (Dickman for polynomial values) | heurnfs | AGGREGATE / KNOWN | Asymptotics over bounded discriminants do not isolate discriminant N. Smooth-value statistics are the NFS heuristic itself (rigorous only for linear forms) | Alberts–Lemke Oliver–Wang–Wood (Crelle 2025) |
| 10.7 | **[cx]** Statistical nontriviality of NFS congruences: for Blum N, ν(v) = (X_v Y_v^{−1} \| N) is a computable character on the relation kernel K, nontrivial ⇔ some congruence splits | [`nfsfruit`](nfsfruit_results.md) (Lemma H: fruitfulness Φ̄ is a homomorphism for *every* odd composite N; Conjecture F; Bockstein/Jacobi row NF-3); lemmaM (Thm DD); nested §10; detGRH | KNOWN open (project's central fruitfulness problem) | The Jacobi-symbol detector for Blum N is correct and elementary, but Lemma H makes the Blum restriction unnecessary. What is needed is Conjecture F, the batch form of LV Conj 7.1 (nfsfruit prereg: Conj 7.1 ⇔ rank Φ ≥ 1), which asks only for probability L[1/3, −o(1)]. The project sidesteps it with Dixon descent at the price of 2/5. Not a new connection | Lee–Venkatesan, JNT 187 (2018) |

---

## 11. Ranked shortlist

**No candidate survives as OPEN-WORTH-A-LOOK at the asymptotic level.** Codex marked three items OPEN. On checking each
against the repository, all three are known project problems:

1. **NFS congruence nontriviality** (10.7, codex #1). This is the strongest *rigorous* lever, but it is Lee–Venkatesan's
   Conj 7.1, already mapped in [`nested`](nested_results.md) §10 and [`detGRH`](detGRH_results.md).
   - Precise first question (codex, sharpened by the critic): for any odd semiprime N (Lemma H removes the Blum
     restriction) and randomized polynomial selection, prove nfsfruit's Conjecture F: the kernel character Φ̄ is
     nontrivial with probability L[1/3, −o(1)].
   - nfsfruit shows that Conjecture F plus LV Thm 2.1 gives L[1/3, (64/9)^{1/3}] for every semiprime. That holds modulo
     LV's non-monic bookkeeping and their reliance on HECM (LPP Part III is unpublished; lv6 §6.5).
   - Not new.
2. **Sub-quadratic exact linear algebra for the MNFS matrix** (7.1, codex #2). It is constant-only: the floor is
   1.8595 for ω ≲ 1.667 (heurnfs; LSTT's 1.8689 is a constrained point).
   - First question: can the projected block-Krylov sequence U^T M^j V (j < 2n/k) be produced in O(n^{2−ε}) bit
     operations for a matrix with one shared rational block and separate field blocks?
   - heurnfs already found that block-angular structure gives nothing at ω = 2, and r21 estimates < 1% for ω ≤ 1.67.
3. **Constructive joint-smoothness transport** (6.3, codex #3). This is heurnfs lever (a), output-sensitive smooth-pair
   enumeration, worth 1.6386 if achieved. [`smoothenum`](smoothenum_results.md) found no mechanism that achieves it. OT contributes
   language only: a Wasserstein-small coupling does not construct the sampler.

If one item from this survey deserves a short follow-up, it is **bookkeeping, not research**:
- Narrow heurnfs row 20's novelty claim. The curve for ω ∈ [1.78, 2] is LSTT 2002 §2.7; the part below 1.78, falling to
  1.8595, corrects LSTT's "optimal relation collection" wording.
- Record the isogeny-root-finding equivalence (4.1, extending Altuğ–Chen Thm 3.3) in APPROACHES as a closed door.

## 12. Reviews

- **GPT-6.1-sol blind round** (`codex/conn_cmp_1.out`, xhigh). Same overall verdict: no fourth family, nothing below 1/3.
  Its 12 candidates not already in our list are merged above and marked [cx]. It independently re-derived the
  MNFS-vs-ω' curve and located LSTT §2.7. Its three OPEN items are re-labelled as known project problems in §11.
- **Opus read-only critic** (verdict **OK WITH FIXES**; no mislabelled OPEN, and nothing killed that should be OPEN).
  All of its findings are applied:
  - M1: the LSTT/heurnfs floor reversed. heurnfs's 1.8595 is right in the shared model, and LSTT's 1.8689 is a
    constrained point.
  - M2: §11.1 now uses the nfsfruit consequence, not nested Theorem NR.
  - M3: 10.7 coverage now cites nfsfruit Lemma H, Conjecture F and lemmaM.
  - M4: 3.5 is an E3(ii)-type target with relaxed norm, not E3 itself.
  - M5: 4.1 is KNOWN in part (Altuğ–Chen Thm 3.3) and restricted to genus-0 X_0(ℓ).
  - Minor: the 4/9 derivation; the 8.3 sparsity (p+q−1) and Prony; the β ≤ 1 KMS value; randomized Jeřábek reductions
    and UEOPL; AGLLV as an analogy only; the Dang et al. gloss; Pilatte 2026; the CKMP title; the smoothenum wording.
  - The critic verified Kayal–Saxena, Bürgisser 2009, Haddad–Koukoulopoulos, Cheon et al. and the LSTT §2.7 numbers
    against sources. It checked §2.1 and §2.5 mathematically.
  - Missed connections it suggested, none asymptotic:
    - Prony/Berlekamp–Massey (added to 8.3).
    - The distribution of ord_N(a) and λ(N) (Erdős–Pomerance–Schmutz), which is the M statistic behind 3.1.
    - A qudit of dimension N, where x ↦ ax is Clifford but controlled exponentiation is not (3.2).

---

### Reproduction (scripts inlined because scratchpad/ is git-ignored)
Outputs: `phi2_rabin.py` → trials=795 splits=352 rate=0.443 (predicted 4/9); `aff_bc.py` → Aff(Z/N) has N classes for N ∈ {6,10,14,15,21,35}; Galois sums 0.861391809069 (N=35), 0.785414868519 (N=33) equal (p^{1−β}−1)(q^{1−β}−1) at β = 2.5.

```python
# scratchpad/conn_cmp/phi2_rabin.py
# Toy check: a root-finder for Phi_2(j0, X) mod N (returning a uniformly random root, CRT-combined)
# splits N when queried at j0 planted via the genus-0 X_0(2) parametrisation j0=(t+16)^3/t, j1=(t+256)^3/t^2.
import random, math
from sympy import randprime
from sympy.polys.galoistools import gf_from_int_poly
def roots_mod(coeffs, p):  # brute-force roots of integer poly (highest first) mod small p
    return [x for x in range(p) if sum(c*pow(x, len(coeffs)-1-i, p) for i, c in enumerate(coeffs)) % p == 0]
def phi2_coeffs(j):
    # classical modular polynomial Phi_2(X,Y) as polynomial in Y with X=j
    X = j
    return [1, -X**2 + 1488*X - 162000, 1488*X**2 + 40773375*X + 8748000000,
            X**3 - 162000*X**2 + 8748000000*X - 157464000000000]
random.seed(1); splits = trials = 0
for _ in range(40):
    p, q = randprime(300, 600), randprime(600, 900); N = p*q
    for _ in range(20):
        t = random.randrange(1, N)
        if math.gcd(t, N) != 1: continue
        ti = pow(t, -1, N)
        j0 = (t+16)**3 * ti % N; j1 = (t+256)**3 * ti * ti % N
        c = phi2_coeffs(j0)
        rp, rq = roots_mod(c, p), roots_mod(c, q)
        assert j1 % p in rp and j1 % q in rq
        a, b = random.choice(rp), random.choice(rq)   # oracle output, independent of t given j0
        r = (a*q*pow(q, -1, p) + b*p*pow(p, -1, q)) % N
        trials += 1
        g = math.gcd((r - j1) % N, N)
        if 1 < g < N: splits += 1
print(f"trials={trials} splits={splits} rate={splits/trials:.3f}")
```

```python
# scratchpad/conn_cmp/aff_bc.py
# (1) number of conjugacy classes of Aff(Z/N) for small squarefree N; (2) Galois sum of Li_beta(e(a/N)) vs prod(p^{1-beta}-1)
import math, cmath
def classes(N):
    U=[u for u in range(1,N) if math.gcd(u,N)==1]; G=[(b,u) for u in U for b in range(N)]
    seen=set(); k=0
    for g in G:
        if g in seen: continue
        k+=1; b,u=g
        for (c,v) in G:  # conj by x->vx+c : (c,v)(b,u)(c,v)^-1
            vi=pow(v,-1,N); # h^{-1}: x-> vi(x-c)
            # h g h^{-1}(x) = v(u(vi(x-c))+b)+c = u x + (v b + c - u c)
            seen.add(((v*b+c-u*c)%N,u))
    return k
print([(N,classes(N)) for N in (6,10,14,15,21,35)])
import mpmath as mp
beta=2.5
for (p,q) in [(5,7),(3,11)]:
    N=p*q
    S=sum(mp.polylog(beta, mp.exp(2j*mp.pi*a/N)) for a in range(1,N) if math.gcd(a,N)==1)/mp.zeta(beta)
    print(N, mp.nstr(S,12), mp.nstr((p**(1-beta)-1)*(q**(1-beta)-1),12))
```
