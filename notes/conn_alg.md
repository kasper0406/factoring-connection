# conn_alg: algebra, arithmetic geometry and dynamics ↔ factoring (phase 28 survey)

Agent `conn_alg`, 2026-10-27. Survey under `notes/AGENT_PHASE28_SURVEY.md`; no pre-registered experiments.
- GPT-6.1-sol blind round: `codex/conn_alg_1.{md,out}` (candidates merged below, marked **[cx]**).
- Opus read-only critic on the OPEN items and KNOWN attributions: §12.
- Tiny exact check (≤ 1 s, no background process): `scratchpad/conn_alg/hurwitz.py` (git-ignored), re-run command in §13.
- Verdict labels as in the brief: RESTATEMENT · AGGREGATE · FACTORING-EQUIVALENT (FACT) · REDUCES-TO-⟨M / A / Q / PT /
  named barrier⟩ · KNOWN · OPEN-WORTH-A-LOOK. Leak families as in `leakmap` §0: **M** group order/smoothness,
  **A** archimedean approximation, **Q** quadratic genus/residue symbol. **PT** = Theorem PT / PT-X (`factorbits` §3.2,
  `leakmap` §0): methods using N only as a modulus learn about the factors only by splitting. "Size N" = the only
  known construction of the object has size ≥ N^{c}, c > 0.
- Claims are labelled PROVED (here, elementary), KNOWN (reference), HEURISTIC, assessment, or UNVERIFIED.

N = pq, p < q odd primes, unless stated.

---

## 0. Verdict

**Plain language: no connection in these ten fields gives, or credibly points to, a faster factoring algorithm.**
We examined 78 concrete connections (64 ours, 14 added from GPT-6.1-sol's blind round; its other 17 candidates
coincided with or were merged into ours). Every one fails at an identifiable step, and the steps are always the same few:

1. **The invariant is the factorisation in disguise (FACT / RESTATEMENT).** Dimension formulas, group orders,
   character counts, K-group orders, motive ranks, Hecke-orbit sizes, trace-formula terms and class numbers of
   level-N objects are functions of p + q (ψ(N), φ(N), (p^j − 1)(q^j − 1), (p² − 1)(q² − 1), …). Every known way to
   evaluate them either uses the formula (needs p) or builds an object of size ≍ N^c.
2. **The computable part only sees N mod m or quadratic symbols (AGGREGATE / Q).** Traces of the Weil
   representation, twists, Artin symbols, root numbers, Fermat-surface Picard numbers, Selmer 2-ranks, K₂ 2-ranks.
3. **The computable part is a pure computation mod N (PT).** N-adic analysis, Witt vectors, heights,
   expander walks, sum-product growth: they reveal p exactly when a gcd would.
4. **The algorithmic content is a known method (KNOWN / M).** K-theory of finite fields = the cyclotomic tori of
   p ± 1 / Bach–Shallit; Deligne–Lusztig tori = Williams' p + 1; Picard–Fuchs recurrences = Bostan–Gaudry–Schost;
   abelian varieties = ECM; Diophantine approximation of √N = CFRAC/SQUFOF; determinant method = Coppersmith;
   zero divisors in quaternion algebras = Rónyai (factoring-hard).

Verdict counts (primary label per row; several rows carry a second label): FACT or RESTATEMENT (incl. factoring-hard) 28;
REDUCES-TO a named family or project barrier 25 (M 7, PT 5, Q 5, A 2, FROB 2, a specific project note 4); KNOWN 16;
AGGREGATE 6;
OPEN-WORTH-A-LOOK 2 (C6, C8) plus 1 OPEN boundary test (L7). All three are *existing* project openings; none is new.

**Shortlist (§11):** (1) nfsfruit's fruitfulness Lemma M (would give L[1/3, (64/9)^{1/3}] via Lee–Venkatesan, modulo
nfsfruit's stated caveats, beating the published rigorous L[1/2, 1] and the project's internal L[2/5]; it does not
beat heuristic NFS); (2) coverage38's positivity count, an
additive-combinatorics statement that would push the deterministic worst case below N^{1/5}; (3) one fixed
form's coefficient a_N (weight 2, or the weight-1 S₃ form η(z)η(23z)), a boundary test with no speed consequence.
None is new to the project; feasibility of each with current tools is low (assessment).

**Surprising or new to the project** (all elementary consequences of classical formulas; no prior factoring-use
statement found for the ★ items, UNVERIFIED as a negative):
- ★ L2: the Kronecker–Hurwitz relation gives **Σ_{t²<4N} H(4N − t²) = 2(N + q)** exactly, and the weight-4 trace
  formula gives Σ(t² − N)H(4N − t²) = −2(1 + p³) (both checked exactly on 120 / 60 semiprimes). The trace formula
  moves the factor from the divisor term into a √N-long class-number moment; it does not hide it.
- R1 / K6 (codex): one scalar k(GL₂(Z/N)) = (N + 1)² − (p + q)², or (p^j − 1)(q^j − 1) = #K_{2j−1}(Z/N), factors N
  (for K6, j = 1 is the classical φ(N); not surprising, but it shows how short a FACT invariant can be).
- ★ A1: canonical heights are PT objects: their non-archimedean part is computable without factoring
  (Müller–Stoll), and log p appears only on events a gcd already detects.
- ★ R2: the irreducible pieces of the Weil representation of SL₂(Z/N) are cut out by the central element with
  CRT components (−1, +1), i.e. a non-trivial square root of 1 mod N.
- ★ H2: the Picard number of the degree-N Fermat surface is 3(N − 1)(N − 2) + 1 for gcd(N, 6) = 1 (aggregate).
- ★ L10: σ(N) mod ℓ = (1 + p)(1 + q) mod ℓ (an Eisenstein coefficient mod ℓ) is an O(1)-bit fiber function giving
  p + q mod ℓ; computable for all small ℓ it would CRT to p + q. No factor-free computation is known. (G6, the
  analogous Galois-cohomology count, carries only the ℓ-support of p − 1 and q − 1: group-order data, M.)

Subjective probability that anything in these fields yields sub-L[1/3] heuristic factoring: < 1% (ours and
codex's verdicts agree).

---

## 1. Langlands program / automorphic forms

| # | Connection (object from N; carrier of factor information) | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| L1 **[+cx L1, L2]** | Hecke eigenvalue a_N(f) of one fixed eigenform (τ(N), a_N(E)); a_N = a_p a_q | tracedata §0 (OPEN for one form), ecount | KNOWN for (a_N, a_{N²}); single a_N: REDUCES-TO tracedata open item | for a fixed *non-CM* eigenform, a_N and a_{N²} (which encodes a_N(Sym²)) together factor most N (density 1, via Serre's bound on a_p = 0) because a_{N²} carries det Frob_p = p (Bach–Charles). Explicitly (cx): with A = a_N, B = a_{N²}, S = A² + N − B = q·a_p² + p·a_q² (weight 2), and the roots of T² − ST + NA² are q·a_p², p·a_q²; gcd with N splits. One coefficient has no determinant (tracedata Prop DEC); codex rates the single-coefficient question OPEN-WORTH-A-LOOK, we keep it as tracedata's existing open item. Edixhoven–Couveignes computes τ(p) for prime p only: its p-dependent step is identifying Frob_p, which over Z/N factors N (nonsmooth Prop FROB). | Bach–Charles, arXiv 0708.1192 (Contemp. Math. 463, 2008); Edixhoven–Couveignes, Ann. Math. Stud. 176 (2011) |
| L2 | **Eichler–Selberg trace formula for T_N at level 1.** For weights k with S_k(1) = 0 (k = 4, 6, 8, 10, 14) it reads Σ_{t²<4N} P_k(t, N)·H(4N − t²) = −Σ_{dd′=N} min(d, d′)^{k−1}. At k = 2 (Kronecker–Hurwitz): **Σ_t H(4N − t²) = 2σ(N) − Σ_{d\|N} min(d, N/d) = 2(N + q)** | not in repo (geometry has the genus/area version) | FACT, KNOWN identity; size √N · N^{1/4} | The elliptic (class-number) side is a sum of ~4√N Hurwitz class numbers of unrelated discriminants 4N − t²; the hyperbolic side is a divisor sum. So the trace formula moves the factors from the divisor term into a long class-number moment, and nothing computes that moment faster than term by term. Checked exactly for 120 semiprimes (k = 2) and 60 (k = 4: Σ(t² − N)H = −2(1 + p³)). | Kronecker 1860, Hurwitz 1885; Zagier's appendix to Lang, *Introduction to Modular Forms*; Cohen 1975 |
| L3 **[=cx L3]** | dim S_k(Γ₀(N)), genus of X₀(N), trace of T_ℓ on S_k(Γ₀(N)) for small ℓ | geometry §0 | dim/genus: FACT, size N. tr T_ℓ (ℓ non-square): REDUCES-TO Q | The main term (k−1)ψ(N)/12 is p + q in disguise; the elliptic terms at fixed ℓ count roots of x² − tx + ℓ mod N, i.e. Legendre symbols of t² − 4ℓ (Eichler–Shimura turns this into #X₀(N)(F_ℓ), the geometry QR tier). Modular-symbol computations of these spaces have dimension ≍ N. | geometry_results §0; Cohen–Strömberg, *Modular Forms* (2017) ch. 12 |
| L4 | Jacquet–Langlands / Brandt matrices: definite quaternion algebra, Eichler order of level N | geometry (Kirschmer–Voight) | FACT, size N | Class number of the Eichler order = (ψ(N)(ℓ−1))/12 + small terms; Brandt matrices have size ≍ N/12. Kirschmer–Voight already reduce factoring to these class numbers. | Kirschmer–Voight, SIAM J. Comput. 39 (2010) Prop 5.8 |
| L5 | Twisted L-values L(f ⊗ χ_N, 1) (Waldspurger/Tunnell: squares of half-integral-weight coefficients c(N), i.e. ternary representation numbers) | geometry (r₃ = 12h), genusleak | REDUCES-TO class groups (KNOWN route: M, L[1/2]) | Root numbers are aggregates (w(f⊗χ_N) = w(f)χ_N(−M)). The value needs conductor ≍ N² terms in the approximate functional equation, or a class-number-type ternary count (r_Q(N) ∝ h(−cN) via Siegel–Weil). Knowing a class number factors via ambiguous forms, which is the class-group method. | Waldspurger 1981; Tunnell 1983; Schnorr–Lenstra 1984; Lenstra–Pomerance, JAMS 1992 |
| L6 | Functoriality transfers that might "carry p" from N: Sym^m, base change to Q(√N), Asai, Rankin–Selberg | tracedata §2.2 | single-index data: RESTATEMENT or REDUCES-TO L1 | Any coefficient at an index divisible by N² (Sym², a_{N³}, …) carries det = p^{k−1} and closes the algebra (BC). Base change to K = Q(√N): N ramifies, the eigenvalue at 𝔭 = (p, √N) equals a_p, but writing 𝔭 down needs p. L(BC f, s) = L(f, s)L(f⊗χ_N, s) adds only L5. | Bach–Charles; Langlands *Base change for GL(2)* (1980) |
| L7 **[+cx S5]** | Weight-1 Artin forms / non-abelian Chebotarev: a_N(ρ) = tr ρ(Frob_p)·tr ρ(Frob_q) for a fixed S₃ or dihedral ρ (e.g. ρ from x³ − x − 1, i.e. η(z)η(23z)); the root count r_N(f) = (1 + a_p)(1 + a_q) is a related, not equivalent, fiber function | tracedata Prop NA / Prop IL (OPEN) | OPEN (tracedata item; boundary test, not a speed route; codex: OPEN-WORTH-A-LOOK) | Abelian ρ is aggregate (Theorem RA). For S₃, a_N ∈ {4, 1, 0, −2} is a non-abelian fiber function: by the argument of Prop NA (stated there for r_N(f) and M mod ℓ) it is not determined by N mod m, but it is O(1) bits, so not factoring-hard under N-local reductions unless factoring is easy (Prop IL). Even an oracle for many ρ gives class data (Chebotarev sets), not CRT-decodable residues. | tracedata §1, §4 |
| L8 | Function-field analogue: over F_q[t] factoring is easy (Berlekamp; Drinfeld modules) | nonsmooth (FROB) | KNOWN analogy; RESTATEMENT of FROB | The function-field algorithms use the global Frobenius x ↦ x^q, which exists because the base field is known. Over Z the analogue x ↦ x^N is Frobenius only when N is prime (AKS). The gap is exactly FROB. | Berlekamp 1970; Cantor–Zassenhaus 1981 |
| L9 **[cx L4]** | Fixed twists f⊗χ as local separators | — | AGGREGATE | a_N(f⊗χ) = χ(N)a_N(f): the twist contributes χ(N) only, a function of N mod cond(χ). | standard twisting |
| L10 | Eisenstein coefficients mod ℓ: σ(N) ≡ (1 + p)(1 + q) (mod ℓ), the N-th coefficient of E₂ mod ℓ (reducible Galois rep 1 ⊕ χ_ℓ) | tracedata fixed-ℓ counts | FACT in aggregate; per ℓ REDUCES-TO PT-m fiber floor | Each ℓ gives p + q mod ℓ, i.e. the hidden fiber point up to swap (O(1) bits); all ℓ ≤ n would CRT to p + q. No computation of one value from N alone is known; this is the reducible case of tracedata's fixed-ℓ counts. | tracedata §4 |

## 2. Representation theory

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| R1 **[+cx R1]** | Character degrees of GL₂(Z/N) = GL₂(F_p) × GL₂(F_q): principal series (p+1)(q+1) = ψ(N), cuspidal (p−1)(q−1) = φ(N), Steinberg pq = N | not in repo as such | FACT, RESTATEMENT | The degree set contains p, p ± 1, q, q ± 1 (e.g. St_p ⊗ a character at q) and ψ(N), φ(N); so is the number of irreducibles k(G) = (p² − 1)(q² − 1) = (N + 1)² − (p + q)² (cx, ★ as a scalar). Producing a representation of the right dimension requires a CRT idempotent (a factor), and the group has order ≍ N⁴: Schreier–Sims or orbit computations cost ≥ N. | Fulton–Harris §5.2; Piatetski-Shapiro, *Complex Representations of GL(2, K)* (1983) |
| R2 **[+cx R2]** | Weil representation of SL₂(Z/N) on C[Z/N]: it is ω_p ⊗ ω_q; irreducible constituents of dimension (p±1)(q±1)/4 | newfamily row 32 (traces: aggregate) | traces: AGGREGATE; constituents: FACT; projectors: RESTATEMENT | Traces are Gauss sums (Gauss's sign formula, computable). The isotypic projectors need the central element acting as (−1, +1) on the CRT components, i.e. diag(e, e) with e² = 1, e ≠ ±1: a non-trivial square root of unity, i.e. a factor. | Gérardin 1977; Cliff–McNeilly–Szechtman 2000 (Weil representations over finite rings); CRT |
| R3 | Weil representation = Clifford group on one qudit of dimension N; Shor uses it | none | KNOWN; REDUCES-TO PT | Clifford circuits mod composite d are classically simulable (Gottesman–Knill for arbitrary d). Shor's speedup needs modular exponentiation, which is not Clifford. So the representation-theoretic part is easy and carries nothing. | Hostens–Dehaene–De Moor, PRA 71 (2005) 042315; Shor 1994 |
| R4 | Wedderburn decomposition: Q[C_N] ≅ ⊕_{d\|N} Q(ζ_d); number of simple components τ(N). Decomposing a semisimple Q-algebra given by structure constants | none | FACTORING-HARD, KNOWN | Rónyai: finding zero divisors in (explicit) quaternion algebras over Q, hence decomposing Q-algebras, is at least as hard as factoring squarefree integers (randomized reduction). Q[C_N] itself has dimension N, so it is not even an instance of polynomial size. | Rónyai, STOC 1987 / J. Algorithms 1988; Ivanyos–Rónyai–Schicho 2012 |
| R5 | Modular representation theory of C_N or GL₂ over F_ℓ: decomposition of F_ℓ[C_N] ↔ factorisation of x^N − 1 over F_ℓ; orbit lengths of ×ℓ on Z/N = ord_d(ℓ) | geometry (NECK/CAT) | REDUCES-TO M | The number of simple F_ℓ[C_N]-modules is Σ_{d\|N} φ(d)/ord_d(ℓ): a necklace count (FACT, geometry Prop NECK) and ord_p(ℓ) is a group order. | geometry_results §5 |
| R6 | Deligne–Lusztig: class numbers and unipotent degrees of GL_n(Z/N) are products of polynomials at p and q, not the polynomial at N (k(GL₂(F_x)) = x² − 1, so k(GL₂(Z/N)) − (N² − 1) = −p² − q² + 2) | none | FACT, size N | The "generic" character table evaluated at N differs from the true one by a symmetric function of p, q; computing the true count is Burnside counting on a group of order ≍ N⁴ (geometry Prop CAT). | Deligne–Lusztig 1976; Green 1955 |
| R7 | Quaternion/LPS picture: two random representations of N as sums of four squares (Rabin–Shallit) and their right gcd in the Hurwitz order | none | REDUCES-TO PT (birthday ⇒ Pollard rho N^{1/4}) | The left ideal Oα corresponds to a pair of lines in P¹(F_p) × P¹(F_q); a gcd is non-trivial only when two lines coincide at one prime (probability ~1/p). Lines are hidden CRT coordinates, so a collision search is Pollard rho. | Rabin–Shallit, CPAM 39 (1986); Lubotzky–Phillips–Sarnak 1988 |
| R8 **[cx R3]** | Deligne–Lusztig tori T_d = ker(Norm: ((Z/N)[√d])^× → (Z/N)^×), local order r − (d/r) | — | KNOWN: M | This is Williams' p+1 group; representation-theoretic repackaging does not change the smoothness distribution. | Williams 1982 |

## 3. Motives, Hodge theory, periods

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| H1 | Gamma values Γ(a/N) (periods of Fermat curves, Chowla–Selberg): distribution relations exist for every divisor m \| N, so integer relations among {Γ(a/N)} "see" p | none | RESTATEMENT (relation length ≥ p) | The multiplication formula for m = p links p values Γ(a/N + jq/N), j < p, with step q/N: writing it needs q, and it has p terms. Rohrlich–Lang predicts no shorter relations, so LLL on high-precision Γ values would need dimension ≥ p. | Kubert 1979; Deligne (Hodge cycles on Fermat), LNM 900; Rohrlich–Lang conjecture |
| H2 | Picard number / Hodge cycles of the Fermat surface x^N + y^N + z^N + w^N = 0 | none | AGGREGATE | ρ = 3(N−1)(N−2) + 1 whenever gcd(N, 6) = 1 (Shioda; Aoki–Shioda: lines span NS ⊗ Q; Degtyarev proved integral generation, answering Shioda): a fixed polynomial in N on that residue class, so it carries no divisor information. The divisor-dependent extra Hodge classes (Shioda's β(m)) occur only through the small divisors 2, 3. | Shioda 1981; Aoki–Shioda 1983; Degtyarev, J. Number Theory 147 (2015), arXiv 1305.3073 |
| H3 | Motive of X₀(N): old/new decomposition h¹ = ⊕_{d\|N} new(d)^{σ₀(N/d)}; Atkin–Lehner involutions W_p | geometry | RESTATEMENT | W_N (z ↦ −1/Nz) is computable; W_p needs an integer matrix of determinant p normalising Γ₀(N), i.e. p. The old/new projectors are built from the degeneracy maps at d \| N. | Atkin–Lehner 1970 |
| H4 | Mahler measures of polynomials built from N, e.g. m(x + 1/x + y + 1/y − N) = r_N·L′(E_N, 0) with r_N rational (CONJECTURED, Boyd; proved for some N, e.g. Rogers–Zudilin) | none | REDUCES-TO A | m(P_N) is computable to high precision (hypergeometric series), but turning it into arithmetic data needs L(E_N, 2) to relative precision ~cond(E_N)^{−1} ≍ N^{−O(1)} (Euler products converge like 1/T), and the conductor sees rad(N² − 16)·N — a different factoring problem, not p. | Boyd, Exp. Math. 7 (1998); Rodriguez-Villegas 1999; Rogers–Zudilin 2014 |
| H5 | Arakelov/period invariants of X₀(N): Faltings height of J₀(N), self-intersection ω² of the dualizing sheaf; exact formulas for squarefree N contain Σ_{p\|N} (p+1)/(p−1)·log p | geometry (area/genus) | FACT, size N | Same as the area of X₀(N): every evaluation goes through the formula (needs factors) or through a model of genus ≍ N/12 (Green's functions, spectral data). | Abbes–Ullmo, Crelle 484 (1997); Ullmo 2000 |
| H6 | Kontsevich–Zagier periods as a complexity class: is p "a short period expression in N"? log p, π, Γ(1/N) … | none | RESTATEMENT | Every known period formula for the factors (Dirichlet class number formula, areas) has an integrand of size ≍ N or a sum over residues mod N; KZ periods are closed under ring operations, so the question is the same as an archimedean (A) formula for p. | Kontsevich–Zagier, *Periods* (2001) |
| H7 **[cx H1]** | Rank of the Artin motive / Betti realisation of Q(ζ_N) | — | FACT | rank = φ(N) gives p + q; explicit bases have ≍ N elements. | Washington, *Cyclotomic Fields* |
| H8 **[cx H2]** | CM periods (Chowla–Selberg) of discriminant −4N to locate ambiguous classes | genusleak | REDUCES-TO M | Period identities do not locate the order-two class; finding it needs class-group relations or an annihilator. | Chowla–Selberg 1967; Shanks 1971 |
| H9 **[cx H3]** | Picard–Fuchs / Cartier–Manin recurrences mod N; scalar case m! mod N, m = ⌊√N⌋ | detGRH, cycloharvey | KNOWN | Bostan–Gaudry–Schost's fast evaluation of polynomial-coefficient recurrences is the Pollard–Strassen N^{1/4} engine; cohomological origin gives no faster evaluation. | Bostan–Gaudry–Schost, SIAM J. Comput. 36 (2007) |

## 4. Algebraic K-theory

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| K1 | Quillen: K_{2i−1}(F_p) ≅ Z/(p^i − 1), so K_*(Z/N) = K_*(F_p) ⊕ K_*(F_q) has orders built from Φ_k(p)Φ_k(q) | detGRH, cycloharvey (cyclotomic orders) | KNOWN: M | The K-groups of Z/N are exactly the algebraic tori whose orders p^i − 1 the p−1 / p+1 / cyclotomic methods exploit. Higher K-theory adds no new group orders. | Quillen 1972; Bach–Shallit, Math. Comp. 52 (1989) |
| K2 | K-theory of Z/N^k via TC (Antieau–Krause–Nikolaus): #K_{2i−1}(Z/p^k)/#K_{2i−2}(Z/p^k) = p^{i(k−1)}(p^i − 1) (this is #K_{2i−1} itself above their even-vanishing range) | none | KNOWN shape: M | Same tori plus p-power parts; the computation works one prime at a time (needs p). | Antieau–Krause–Nikolaus, arXiv 2204.03420 (announcement; Cor. extending Angeltveit) |
| K3 **[+cx K2]** | K₁(Z[1/N]) = units ±p^a q^b (rank ω(N)); K₂(Z[1/N]) ≅ K₂(Z) ⊕ ⊕_{p\|N} F_p^* (Tate), order exactly 2φ(N); tame symbols ∂_p | none | RESTATEMENT / FACT | A basis of the S-units, a tame symbol at p, or #K₂(Z[1/N]) all need p (#K₂ gives φ(N)). | Milnor, *Introduction to Algebraic K-Theory* (1971) §11 |
| K4 **[+cx K3]** | Tame kernel K₂(O_F), F = Q(√N): Tate/Browkin–Schinzel 2-rank formula r₁ + g(2) − 1 + (2-rank of a class group of O_F[1/2]); Birch–Tate #K₂O_F = w₂ζ_F(−1) with ζ_F(−1) = ζ(−1)L(−1, χ_D), D = disc Q(√N) ∈ {N, 4N} | none | REDUCES-TO Q (2-rank) and A (order) | The 2-rank is a genus count (ω(N)); 4-ranks are Rédei-type Legendre matrices (as for Cl). L(−1, χ_D) is a rational of size N^{3/2}; exact evaluation by functional equation + Euler product at s = 2 needs N^{3/2} primes, or ≍ N^{1/2} terms. And unlike Cl there is no efficient group law on K₂O_F, so even a known order would not give an ambiguous-form analogue. | Browkin–Schinzel, Crelle 331 (1982); Tate 1976; Belabas–Gangl, K-Theory 31 (2004) |
| K5 | Whitehead group of the cyclic group: rank Wh(C_N) = ⌊N/2⌋ + 1 − τ(N) | none | AGGREGATE (for semiprimes) | τ(N) = 4 for every semiprime, so the rank formula carries nothing; torsion of K-groups of Z[C_N] needs the cyclotomic components Z[ζ_d], d \| N. | Bass 1965; Oliver, *Whitehead Groups of Finite Groups* (1988) |
| K6 **[cx K1]** | The scalar #K_{2j−1}(Z/N) = (p^j − 1)(q^j − 1) for any fixed j | — | FACT (j = 1 is φ(N), classical; larger j are the Bach–Shallit torus orders) | p^j + q^j = N^j + 1 − C_j, then a quadratic and an integer root. Quillen's formula is evaluated after the local decomposition; no factor-free evaluation exists. Combined with K1: the *orders* are FACT, the *groups* are the M-family tori. | Quillen 1972 |

## 5. Arithmetic geometry (heights, Arakelov, abc, Chabauty, rational points)

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| A1 | Néron–Tate heights on curves with bad reduction at p and q: ĥ(P) = λ_∞ + Σ_v λ_v, and λ_p at an I_n fibre carries rational multiples of log p | none | REDUCES-TO PT | The non-archimedean part is computable without factoring (Müller–Stoll: gcd/denominator arithmetic, quasi-linear). The p-specific correction appears only when P meets the singular point mod p but not mod q, which is a gcd split of N already. So heights are PT objects: they reveal log p exactly when a gcd would. | Silverman, Math. Comp. 66 (1997); Müller–Stoll, LMS JCM 19 (2016), arXiv 1509.08748 |
| A2 **[=cx G3]** | abc / Szpiro / S-unit counting near kN | abcsunit (full study) | REDUCES-TO abcsunit (no barrier, no target) | abc never excludes primitive relations of height ≤ N^{1+ε}; S-unit listing = Λ_N short vectors. | abcsunit_results §0 |
| A3 **[=cx G2]** | Integral points on norm-form / Thue curves F(x, y) = N: two non-associated solutions α, β of N_{K/Q} = N give gcd ideal of norm p (Euler's two-representations method) | geometry (x² + y² = N) | KNOWN; RESTATEMENT for finding points | Two essentially different points factor N (Euler; generalised to any norm form). Finding one point needs a root of the form mod N (Cornacchia), already as hard as a square root of −D mod N; Thue–Mahler solvers factor the right-hand side first. | Euler; McKee, "Turning Euler's factoring method into a factoring algorithm", BLMS 1996 |
| A4 | Rational points and cusps of X₀(N): rational cusps ↔ divisors d \| N (2^{ω(N)} of them) | geometry | RESTATEMENT | The cusp 1/d needs d. Non-cuspidal rational points exist only for finitely many N (Mazur, Kenku; none for N > 163), so for large N nothing else is rational. | Ogg 1974; Mazur 1977; Kenku 1981 |
| A5 | 2-Selmer groups of congruent-number twists E_N: y² = x³ − N²x (and quadratic twists by N in general) | geometry (linking number, open) | REDUCES-TO Q | Monsky's matrix: the 2-Selmer rank is the corank of a matrix of Legendre symbols among −1, 2, p, q. The descent itself lives in Q(S, 2) with S ∋ p, q (needs factors); parity/root numbers are aggregates. The non-aggregate content is (p/q), the "linking number" left open in geometry. | Monsky (appendix to Heath-Brown, Invent. 1994); Rédei 1934 |
| A6 | Chabauty–Coleman / Kim on curves built from N (y² = x⁵ + N, superelliptic twists) | none | RESTATEMENT | The p-adic integration is at an auxiliary good prime ℓ ∤ N chosen by the solver; the Mordell–Weil input (2-descent) needs the factorisation of the discriminant, which contains N. Chabauty consumes factorisations, it does not produce them. | Coleman 1985; Balakrishnan et al. 2019 |
| A7 | Arakelov class group of Q(√N) (class group + regulator) | genusleak, classE | KNOWN: class-group methods (M, Q) | Computing the Arakelov class group is the subexponential class-group machinery; ambiguous classes give the factor. | Schoof, "Computing Arakelov class groups" (2008); Buchmann 1990 |
| A8 **[cx G1]** | Abelian varieties / Jacobians over Z/N with smooth multipliers | — | KNOWN: M (ECM, HECM) | Success is smoothness of one local group order; higher dimension gives more orders but costlier groups and no better order distribution. | Lenstra 1987; Cosset 2010 |

## 6. p-adic analysis, perfectoid spaces, prismatic cohomology

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| P1 | N-adic analysis: Z_N = lim Z/N^k ≅ Z_p × Z_q; N-adic logs, Hensel lifts, Fermat quotients (a^{N−1} − 1)/N | factorbits, leakmap NE1 | REDUCES-TO PT (and NE1 for exponent N−1) | Ring operations in Z/N^k with zero tests mod N^k are a pure tree over a fixed Z/N-algebra (PT-X): component-blind. Fermat quotients are the Fermat test, a group-order (M) event. | factorbits §3.2; leakmap §0 |
| P2 **[=cx P3]** | Frobenius lifts / δ-ring structures "at N" (φ ≡ x^p mod p and x^q mod q); Buium's arithmetic derivative at N | nonsmooth Prop FROB | REDUCES-TO FROB (FACT) | Such a φ is the hidden componentwise Frobenius; producing it on an étale Z/N-algebra factors N. The big-Witt Frobenius exists over any ring but is not this operator. | nonsmooth §5.2; Buium, *Arithmetic Differential Equations* (2005) |
| P3 | Inseparability: x^N − 1 ≡ (x^q − 1)^p (mod p) and (x^p − 1)^q (mod q); its derivative N x^{N−1} ≡ 0 mod N, so over Z/N it "is" a p-th power at p and a q-th power at q | none | RESTATEMENT | Extracting the root, i.e. writing x^N − 1 = G(x^p), needs the exponent p. Perfectoid tilting exploits exactly the bijectivity of x ↦ x^p in characteristic p; with two unknown characteristics there is nothing to tilt. | Scholze 2012 (perfectoid spaces) for the tilting mechanism |
| P4 | Kedlaya/Dwork/Lauder point counting, Gross–Koblitz (Γ_p), Morita Γ at composite N | nonsmooth (Kedlaya/Lauder/Satoh row), track C | REDUCES-TO FROB | A Frobenius lift mod N^k reduces mod N to Φ (FROB); Gross–Koblitz needs Γ_p, whose N-adic analogue is the factorial unit predicate. | nonsmooth §4 |
| P5 **[=cx P1]** | Prismatic/syntomic cohomology and TC of Z/N or of Z/N-algebras | none | REDUCES-TO K1/K2 (M) | Spec Z/N is a disjoint union of Spec F_p and Spec F_q; prismatic cohomology is defined relative to one prime at a time, so every invariant is a product of per-prime invariants and is computed only after choosing p. Global invariants (orders of TC/K groups) are the cyclotomic tori of K1. | Bhatt–Scholze, Ann. Math. 196 (2022); AKN (K2) |
| P6 **[cx P2]** | Bounded truncations of big Witt vectors W_S(Z/N), ghost map | nonsmooth | REDUCES-TO PT | Universal polynomial operations; inverting the ghost map fails exactly when an index integer is a non-unit, which is a gcd split; enumerating indices up to p is a search. | Borger 2011 |

## 7. Model theory and o-minimality

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| M1 **[+cx O1]** | Pila–Wilkie / Bombieri–Pila counting of rational points on sets built from N, and the determinant method for integer points near xy = N | floorsum, ptesplit, nonexact (Bombieri–Pila), topblock, coverage38 | KNOWN (exponential regime); Pila–Wilkie itself: RESTATEMENT | The divisors are integer points on the *algebraic* curve xy = N, which Pila–Wilkie removes by definition (the algebraic part). The algebraic counting/determinant method is the auxiliary-polynomial machinery of Coppersmith (known-bits factoring), i.e. the exponential N^{1/4}-type regime and the topblock barrier. | Bombieri–Pila, Duke 59 (1989); Pila–Wilkie, Duke 133 (2006); Coppersmith, J. Cryptology 10 (1997); topblock_results |
| M2 | Point counts of definable sets over Z/N: by Chatzidakis–van den Dries–Macintyre, \|X(F_p)\| = μp^d + O(p^{d−1/2}), so \|X(Z/N)\| = \|X(F_p)\|·\|X(F_q)\| (e.g. squares: (p+1)(q+1)/4) | tracedata (EXT), geometry | FACT | Counting a definable set over Z/N is a point count of the type #E(Z/N), which factors N (tracedata EXT; for squares directly). No definable counting procedure avoids the factors. | Chatzidakis–van den Dries–Macintyre, Crelle 427 (1992) |
| M3 | Ax–Kochen / ultraproduct transfer from F_p((t)) to Q_p; decidability of the theory of all finite fields | none | RESTATEMENT (no algorithmic transfer) | Transfer moves first-order truths for almost all p, not running times; "N is composite" is already decidable (AKS). The function-field analogue of factoring is easy because of Frobenius (L8). | Ax–Kochen 1965; Ax 1968 |
| M4 **[+cx O2]** | Effective/sharp Pila–Wilkie (Binyamini–Novikov, polylog in complexity) for Pfaffian sets encoding N | none | RESTATEMENT | Any definable set whose rational points are the divisors contains them in its algebraic part (xy = N); transcendental sets built from N (y = N^x, sin(πN/x) on bounded domains) have only trivial rational points. | Binyamini–Novikov, Ann. Math. 2019 |
| M5 **[cx O3]** | Constant-size formula ∃e (e² = e, e ≠ 0, 1) over Z/N | — | FACT | Deciding the formula is trivial on the promise; any witness gives gcd(e, N). Low logical complexity says nothing about witness extraction. | — |
| M6 **[cx O4]** | Ax–Kochen-type local formulas in N at fixed known primes ℓ | — | AGGREGATE | Bounded-precision answers factor through N mod ∏ℓ^{k_ℓ}; using the unknown primes requires naming them. | Ax–Kochen 1965 |

## 8. Additive combinatorics

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| C1 **[+cx C1]** | Inverse sum-product in Z/N: sets with small \|A + A\| + \|AA\| that are not concentrated on residue classes of a large divisor expand (Bourgain for arbitrary q; Bourgain–Chang for few prime factors), so non-expanding sets sit on the structure of a divisor | smoothbuild §3 (sum-product gives distribution, not navigation) | RESTATEMENT | The exceptional \"approximate subrings\" are pZ/N and qZ/N; any set that triggers the inverse theorem already lives on a divisor's residue structure, so writing it down needs p. Sets built from N by ring operations are pure trees (codex: REDUCES-TO PT): growth fails only at a non-unit or an unequal local collision, which is a gcd split. | Bourgain, J. Anal. Math. 106 (2008); Bourgain–Chang, GAFA 16 (2006) (exact composite-modulus reference UNVERIFIED) |
| C2 | Expansion of SL₂(Z/N) Cayley graphs, squarefree N (Bourgain–Gamburd–Sarnak); random walks mod N | prim, leakmap | REDUCES-TO PT (and M for orders) | Walks are pure trees: uniform expansion hides the CRT components. Escape to I mod p but not mod q has probability ≍ p^{−3}; element orders mod p divide p(p²−1) (p±1 methods). | Bourgain–Gamburd–Sarnak, Invent. 179 (2010); Williams 1982 |
| C3 | Polynomial method (Stepanov, Croot–Lev–Pach, Ellenberg–Gijswijt, slice rank) | topblock, floorsum | KNOWN where algorithmic (Coppersmith = auxiliary polynomial); otherwise RESTATEMENT | The auxiliary-polynomial idea is already Coppersmith's method (exponential regime). CLP/EG bound sizes of progression-free sets; no search step and no N-dependent object. | Coppersmith, J. Cryptology 10 (1997); Ellenberg–Gijswijt, Ann. Math. 185 (2017) |
| C4 | Higher-order Fourier analysis: Gowers norms / nilsequences of functions g(a^x) or of residues mod N | periodfind Theorem T | REDUCES-TO M | Period-concentrated (linear) Fourier mass exists only if g predicts a large-order character (Theorem T); a nilsequence correlation of period ord_p(a) is again an order (M). The Jacobi symbol's Fourier transform vanishes exactly on non-units, but evaluating it there already needs gcd(ξ, N) > 1. | Green–Tao 2012; periodfind_results |
| C5 **[=cx C2]** | Freiman/Bohr-set structure of smooth numbers mod N (smooth u with u mod N small) | suffmap, sihard, siconstruct | REDUCES-TO sihard (L[2/3] per output) | This is the Bohr-set lister / Prop SI; oblivious constructors cost L[2/3]. | suffmap Prop SI; sihard Theorem J |
| C6 | Products of band primes in short residue sets mod q (p·ℓ̄ mod q in Farey arcs): Bourgain–Garaev, BGK multiplicative energy | coverage38 | **OPEN-WORTH-A-LOOK** only as coverage38's existing conjecture (no literature input known) | This is the deterministic worst-case covering statement (λ < J^{1/4}) of the project's internal Theorems R/U (published record: N^{1/5}, Harvey–Hittmeir). coverage38 checked Bourgain–Garaev, Korolev: all at q^{1/2+ε} or log savings. Pointwise inverse-prime Kloosterman bounds are uniform in the dilate p, so they are exactly the generic statistics coverage38's Prop L rules out; a proof must use the arithmetic of x_ℓ = p·ℓ̄. | coverage38_results §0 |
| C7 **[cx C3]** | CLP / slice rank to force square products among smooth exponent vectors in (Z/4)^r | heurnfs, suffmap | REDUCES-TO size–entropy | r + 1 parity vectors already force a dependency; progression thresholds add nothing and still need the smooth inputs. | Croot–Lev–Pach 2017; Balister–Bollobás–Morris, Ann. Math. 188 (2018) |
| C8 **[cx C4]** | **Fruitfulness of rigorous NFS (nfsfruit Lemma M), in Fourier form.** nfsfruit reduces rigorous NFS fruitfulness to δ = min_w Pr[β(z) ≠ ⟨w, x(z)⟩] ≥ L[1/3, −o(1)] on a set of (m, f) of measure L[1/3, −o(1)] under the search-output law (Lemma M), where x(z) ∈ F₂^C are the square-class features of a kept pair (including the constant feature when \|S\| parity is imposed) and β is the explicit Bockstein function. Since E[(−1)^{β+⟨w,x⟩}] = 1 − 2δ_w identically, this reads max_w E[(−1)^{β(z)+⟨w,x(z)⟩}] ≤ 1 − 2L[1/3, −o(1)] (no absolute value needed when the constant feature is present). | nfsfruit (OPEN; \"feasibility low\") | **OPEN-WORTH-A-LOOK** (existing opening; the Fourier form is a restatement, not a new handle) | A proof would give L[1/3, (64/9)^{1/3}] for every semiprime via Lee–Venkatesan + nfsfruit Lemma H (modulo LV's non-monic bookkeeping and codex's restart lemma; needs \|X\| ≫ C/δ), beating the published rigorous L[1/2, 1] and the project's internal L[2/5] (Theorem U, not human-refereed) and conditional L[3/8]; it does not beat heuristic NFS. The failing step: LV's equidistribution is first-order (mod-2 characters); β is a mod-4 / second-order quantity, and uniformity over all 2^C linear phases (C ≍ L[1/3]) under two-sided smoothness conditioning is out of reach of moment bounds. nfsfruit §4 concludes a structural (involution-type) argument is needed; that is the reason to look, since it is the largest payoff on this list. | nfsfruit_results §3–4; Lee–Venkatesan, arXiv 1805.08873 |

## 9. Ergodic theory and homogeneous dynamics

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| E1 **[=cx D1]** | Equidistribution of Hecke orbits T_N·z (Clozel–Oh–Ullmo, Eskin–Oh) | geometry | FACT, size N | The orbit has ψ(N) points (p + q in disguise); its construction is X₀(N) → X(1), degree ψ(N). | Clozel–Oh–Ullmo, Invent. 144 (2001) |
| E2 **[=cx D2]** | Duke: Heegner points / closed geodesics of discriminant −4N, 4N; genus cosets | classE (Heegner "floor" withdrawn), genusleak | REDUCES-TO Q / genusleak Theorem GB | Equidistribution is a statement about the whole class group; the factor lives in ambiguous classes (2-torsion), reachable only through group orders or relations. | Duke, Invent. 92 (1988); Harcos–Michel 2006 |
| E3 | Furstenberg ×a×b: orbits of ⟨a, b⟩ ⊂ (Z/N)^* acting on N^{−1}Z/Z (effective Rudolph–Johnson) | none | REDUCES-TO M | Effective equidistribution (BLMV; Bourgain for general N) fails exactly when ⟨a, b⟩ mod p is small, i.e. when ord_p is small: the p−1 method event. | Bourgain–Lindenstrauss–Michel–Venkatesh, ETDS 29 (2009) |
| E4 | Horocycle / Euclid points (a + i)/N on Γ\H; equidistribution of a/N continued fractions | classE (Prop VIS, Theorem F-SQ) | REDUCES-TO classE | Reducing (a + i)/N is the Euclid path of a/N; a convergent sees p only if a ≡ uq ± 1. | classE_results |
| E5 | Cat-map / toral automorphism periods mod N; mixing of x ↦ Ax on (Z/N)² | none | KNOWN: M (p+1/Lucas) | Period = order of A in GL₂(Z/N) = lcm of orders mod p, q; dividing p − 1 or p + 1: Williams' p+1 method. | Williams, Math. Comp. 39 (1982) |
| E6 | Littlewood / Diophantine approximation of √N | — | KNOWN (CFRAC, SQUFOF) | The continued fraction of √N is the classical source of small quadratic residues. | Morrison–Brillhart 1975; Shanks (SQUFOF) |
| E7 | Spectral gap / Ratner for congruence covers Γ(N)\H | none | AGGREGATE (uniform) | Selberg/Kim–Sarnak gaps are uniform in N; equidistribution rates carry no N-specific information beyond the cover's degree (FACT, size N). | Kim–Sarnak 2003 |
| E8 **[cx D3]** | Littlewood/EKL measure rigidity for lattices from (a/N, b/N): short vectors with a proper gcd | classE | REDUCES-TO A | Rational points with public denominator N; a useful return needs an extra divisibility condition that approximation does not enforce, and rigidity gives no constructive hitting time. | Einsiedler–Katok–Lindenstrauss, Ann. Math. 164 (2006) |
| E9 | Arithmetic dynamics: iterate x ↦ x² + c mod N and detect a cycle mod p (Pollard rho; Floyd/Brent) | throughout (pure trees) | KNOWN (exponential, ≈ p^{1/2} = N^{1/4}); REDUCES-TO PT | The most direct dynamics-to-factoring link: the orbit mod p is a random-mapping orbit whose cycle appears after ≈ √p steps and is detected by gcd. Its rigorous analysis is itself open (Bach 1991; Kim–Montenegro–Peres–Tetali for the discrete-log rho), but the exponent is fixed by the birthday bound. | Pollard 1975; Bach, Inf. Comput. 90 (1991) |

## 10. Galois cohomology and class field theory

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| G1 **[+cx S3]** | Zero divisors in quaternion algebras over Q / rational points on conics ax² + by² = Nz² | none | FACTORING-HARD, KNOWN | Rónyai: factoring squarefree integers reduces (randomized) to finding zero divisors in quaternion algebras over Q; the explicit isomorphism with M₂(Q) is at least as hard as factoring. Known algorithms factor the discriminant first (Simon; Ivanyos–Rónyai–Schicho). | Rónyai, STOC 1987 and J. Algorithms 9 (1988) 494–506 (venue from memory); Kutas, arXiv 1606.01053 (secondary summary) |
| G2 | Genus field / Hilbert class field of Q(√−N); ray class field of conductor N (Q(ζ_N), degree φ(N)) | genusleak | RESTATEMENT / FACT | The genus field is Q(√p*, √q*): its defining polynomial gives p. [Q(ζ_N):Q] = φ(N). CM construction of the Hilbert class field costs h ≍ √N. | Cox, *Primes of the form x² + ny²* §6 |
| G3 **[+cx S4]** | Class groups and Cohen–Lenstra: 2-rank = ω − 1, 4-rank = Rédei matrix of (p/q) | genusleak, geometry | KNOWN (class-group methods, M + Q); 4-rank: REDUCES-TO Q (open (p/q)) | Class-group factoring (Schnorr–Lenstra, Lenstra–Pomerance rigorous L[1/2]) is known; Cohen–Lenstra is a distributional heuristic with no search procedure. | Schnorr–Lenstra 1984; Lenstra–Pomerance 1992; Rédei 1934 |
| G4 **[=cx S1]** | Reciprocity: Artin symbols of the rational ideal (N) | newfamily Theorem RA | AGGREGATE | Depends only on N mod conductor. | newfamily §0 |
| G5 | Higher residuosity of fixed small integers mod N (is 2 a cube mod N?); cubic/quartic symbols | leakmap E2 (cubic mod 7 null), crypto literature | REDUCES-TO Q-type fiber function (believed hard; a standard assumption) | The symbol of the rational N is trivial by RA; the residuosity itself is a hidden fiber function like ε. It is an assumption in Benaloh / Naccache–Stern style cryptosystems, not a speed route. | Benaloh 1994; Naccache–Stern, CCS 1998 |
| G6 | Galois cohomology of Z[1/N]: for odd primes ℓ ∤ N, dim H¹(Z[1/N], Z/ℓ) = #{p \| N : p ≡ 1 mod ℓ} = [ℓ \| p − 1] + [ℓ \| q − 1] (cyclic degree-ℓ fields unramified outside N, by Kronecker–Weber); with Z/ℓ^k coefficients, the unordered pair {v_ℓ(p − 1), v_ℓ(q − 1)} | none | REDUCES-TO M (FACT in aggregate) | This is group-order data: the ℓ-support of p − 1 and q − 1, about (log ℓ)/ℓ bits per small ℓ, so small ℓ alone can never pin p. Over all prime powers ℓ^k (including the largest prime of φ(N)) the data give φ(N): FACT. Per ℓ the only known computation builds the extensions from p. | Neukirch–Schmidt–Wingberg ch. 8 |
| G7 | Selmer groups of twists by N (see A5); BSD/Tunnell | geometry | REDUCES-TO Q / class groups | See A5 and L5. | — |
| G8 **[cx S2]** | Explicit 2-descent matrices for E_N: y² = x(x − N)(x + N) | — | RESTATEMENT (for the algorithm); scalar rank: REDUCES-TO Q | Explicit descent builds square-class groups supported at the bad primes, i.e. uses p, q. The scalar Selmer rank is Monsky's Legendre-symbol corank (A5). | Poonen–Schaefer 1997; Monsky 1994 |

## 11. Ranked shortlist (OPEN-WORTH-A-LOOK)

No new connection survives. The items below are existing project openings that these fields touch; they are
ranked by what they would buy, not by how likely they are. Only #1 and #2 are OPEN-WORTH-A-LOOK.

**1. C8 — fruitfulness of rigorous NFS (nfsfruit Lemma M; codex's #1).** Would give L[1/3, (64/9)^{1/3}] for every
semiprime via Lee–Venkatesan + nfsfruit Lemma H (modulo LV's non-monic bookkeeping and codex's restart lemma),
beating the published rigorous L[1/2, 1] (Lenstra–Pomerance) and the project's internal, not human-refereed
L[2/5]; no effect on heuristic NFS.
- *Restatement.* Under the kept-pair law μ on the smooth set X_{n,m,f}, with square-class features x(z) ∈ F₂^C and
  the Bockstein function β (nfsfruit Prop G), show
  max_{w ∈ F₂^C} E_μ[(−1)^{β(z) + ⟨w, x(z)⟩}] ≤ 1 − 2L[1/3, −o(1)]
  for (m, f) in a set of search-output measure L[1/3, −o(1)]. This is an identity-level restatement of Lemma M
  (E[(−1)^{β+⟨w,x⟩}] = 1 − 2δ_w), not a new handle; nfsfruit §4 argues moment/Fourier bounds over 2^C phases cannot
  reach it and that a structural, involution-type argument is needed.
- *First lemma (one week).* Find an involution z ↦ z′ on a subset of X of μ-measure ≥ L[1/3, −o(1)] that preserves
  x(z) and flips β (nfsfruit notes Dixon's ε-multiplication is such an involution for Dixon; none is known for
  NFS pairs). Codex's variant: for an irreducible cubic field, express β modulo the row space through ideal
  valuations mod 4 and at most two unit characters, *without* an S-unit basis, or isolate the class-group
  obstruction. Feasibility: low (nfsfruit; assessment).

**2. C6 — coverage38's positivity count.** Would push the deterministic worst case below N^{1/5} (Theorem R).
- *Statement* (coverage38 §0): a uniform positivity count of band primes ℓ with p·ℓ̄ mod q in a Farey-arc set of
  measure J^{−3/4}, i.e. λ < J^{1/4}.
- *First question.* Do Bourgain–Garaev-type multiplicative-energy bounds for {ℓ̄ mod q : ℓ prime in a band}
  extend from length q^{1/2+ε} to q^{0.32} for *positivity in a union of arcs* (not asymptotics)? coverage38's
  Prop L shows generic deviation bounds cannot suffice, so the input must use the arithmetic of x_ℓ = p·ℓ̄.
  Feasibility: low; no literature handle known.

**3 (boundary test, not OPEN-WORTH-A-LOOK). L1 / L7 — one fixed form's coefficient at N** (tracedata's open item;
codex's #2 and #3).
- *Question (codex).* For η(z)η(23z) = ½(Θ_{x²+xy+6y²} − Θ_{2x²+xy+3y²}), write a complete factor-free
  specification for b_N mod 2 and audit the fastest quadratic-form route: which modular square roots or
  representation searches does it call, and can they be removed? For weight 2 (e.g. 37a1): given A = a_pa_q ≠ 0
  and its factorisation, can one extract p on an explicit positive-density family?
- *Why it is not a speed route.* b_N is O(1) bits, so it is not factoring-hard under N-local reductions unless
  factoring is easy (tracedata Prop IL). A weight-2 a_N has about n/2 bits, so Prop IL does not apply; but the
  question there is hardness (does one a_N factor N?), and no algorithm computing a_N without the factors is known,
  so neither answer would speed up factoring.

## 12. Critic

One read-only Opus critic reviewed the OPEN items, every KNOWN attribution and the ★ claims (2026-10-27).
**Verdict: ACCEPT WITH FIXES; all fixes applied.** It re-verified L2 numerically (55 semiprimes; k = 4, 6, 8, 10,
14), the Fermat Hodge count at m = 5, 7, 25, 35, and the theta identity for η(z)η(23z) through q³⁰⁰. It checked
Bach–Charles, AKN, Müller–Stoll, Degtyarev, Kutas, Kirschmer–Voight and Neuhauser against the sources; the remaining
citations were checked from its own knowledge.
- MAJOR, G6: the cohomology count gives the ℓ-support of p − 1 and q − 1, not p mod ℓ, so it cannot CRT to p.
  Relabelled M (FACT only in aggregate, over all ℓ^k); the ★ now covers L10 only.
- MAJOR, C8/§11.1: the payoff was under-qualified. Fixed: it beats the published L[1/2, 1] and the internal,
  unrefereed L[2/5]; nfsfruit's caveats are restated, and so is the measure quantifier.
- MAJOR, C8: the Fourier form is an identity-level restatement, and nfsfruit §4 says moment bounds will not reach
  it. Fixed: C8 is no longer called a "new handle", and the reason to look is now the involution criterion.
- MINOR, all fixed:
  - L7 relabelled as a boundary test, with the exact value set; r_N is "related", not equivalent; NA now reads
    "by the same argument".
  - §11.3 weight-2 sentence corrected.
  - K2: the AKN statement is a ratio.
  - H2: attributions fixed (Shioda / Aoki–Shioda over Q; Degtyarev integral), and "mod 6" replaced by "polynomial
    in N".
  - R1 degree set corrected.
  - R2: the Neuhauser reference is replaced by Cliff–McNeilly–Szechtman.
  - K6: the ★ is dropped (j = 1 is φ(N)).
  - A4: Kenku added, with N > 163.
  - K3: #K₂ = 2φ(N) exactly.
  - K4: χ_D, D = disc Q(√N).
  - H4: rationality marked conjectural.
  - L1: non-CM hypothesis added.
  - C1: reference marked UNVERIFIED.
  - G1: primary Rónyai reference added.
  - E9: Pollard rho row added.
  - C6: internal vs published deterministic record noted.
- No row was upgraded. The critic agreed that C6 and C8 are the only OPEN-WORTH-A-LOOK items.

## 13. Reproduce

```bash
.venv/bin/python -I scratchpad/conn_alg/hurwitz.py   # Kronecker–Hurwitz (120 semiprimes) and k = 4 trace identity (60)
```
The script is in the git-ignored scratchpad; its content is reproduced here for permanence:
H(D) is the Hurwitz class number (reduced forms of discriminant −D, weights 1/2 and 1/3 for the forms
(a, 0, a) and (a, a, a)); the checks are Σ_{|t|<2√N} H(4N − t²) = 2σ(N) − Σ_{d|N} min(d, N/d) = 2(N + q) and
Σ (t² − N)H(4N − t²) = −2(1 + p³), for N = pq with 3 ≤ p < q < 60.
