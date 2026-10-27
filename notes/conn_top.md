# conn_top: knots, Floer, homotopy, categories, combinatorics, expanders, tropical, tilings, logic, automata (phase 28)

Agent `conn_top`, 2026-10-27. Survey brief: `notes/AGENT_PHASE28_SURVEY.md`. **Version v1 (codex blind round merged; Opus
critic ACCEPT WITH FIXES, all applied).**
- GPT-6.1-sol blind round: `codex/conn_top_1.{md,out}` (launched before any analysis here; merged in §11, credited
  "cx" in the tables).
- Scratch: `scratchpad/conn_top/checks.py` (tiny exact sanity checks, quoted where used; all passed).
- Labels: **KNOWN** (reference), **PROVED** (routine proof here; "no prior formulation found by us" is the strongest
  novelty claim), **HEURISTIC**, **OPEN**, **UNVERIFIED**. Negative literature statements ("not found") are UNVERIFIED
  as negatives.

Notation: N = pq, p < q odd primes, n = log₂N. Verdict codes as in the brief: **RESTATE**, **AGG**, **FACT**
(randomized-polynomial-time equivalent to factoring; "FACT-hard" if only one direction is known), **→M / →A / →Q**
(reduces to group order/smoothness, archimedean approximation, quadratic genus), **→⟨barrier⟩** (a named project
barrier), **KNOWN**, **OPEN-WORTH-A-LOOK**. "EXP" marks a feature whose only known evaluation has size N^{Ω(1)}.

---

## 0. Verdict

**No connection in these ten fields gives, or plausibly leads to, factoring below L[1/3, 1.902], below the rigorous
L[2/5], or below the deterministic N^{1/5}.** We examined 81 candidates: 69 of ours, 8 new ones from GPT-6.1-sol's
43-candidate blind round (its other 35 overlap ours) and 4 added by the Opus critic. Every one fails at an identifiable step. Most land in tiers the project already knows:
- **RESTATE:** writing the object down needs a factor. Examples: tiles of Z_N (|A| ∈ {p, q}), torus knots T(p, q),
  matroid characteristic sets, Müger decomposition of C(Z/N).
- **AGG:** fixed-level TQFTs, by Ng–Schauenburg congruence, and fixed tiles.
- **FACT, and computable only at size N:** the unitary Cayley graph spectrum, level-N TMF ranks, K-theory of Z/N,
  the E₈ theta coefficient σ₃(N), subobject counts (p + 3)(q + 3).
- **Families M / A / Q:** image-of-J denominators, Wendt/dimer resultants, sampling-based spectral estimates, lens
  homotopy type.
- **Known barriers:** factorial/Pochhammer for tiling counts, Theorem PT for Floer and matroids over Z/N, suffmap
  Tier C for spectral or matroidal linear algebra.

Logic gives only one-way transfers. Proof-system properties or provabilities that would yield factoring algorithms
(Krajíček–Pudlák, Bonet–Pitassi–Raz, Liu–Pass) are believed false for exactly that reason. The October 2026 ETH
refutation changes nothing for factoring: the SAT encoding of a known multiplier circuit has m = Ω(n log n) variables,
so 2^{m/√log m} ≥ 2^{Ω(n√log n)} is worse than trial division (L5).

**New to the project, or surprising.** None of these is a route to faster factoring.
1. **Lens-space homotopy type is quadratic residuosity (K1).** L(N, a) ≃ L(N, b) iff ab ≡ ±n² (mod N) (Reidemeister
   1935, Whitehead 1941). Deciding *oriented* homotopy equivalence of poly-size triangulated lens spaces is therefore
   exactly QRP under a randomized reduction (the unoriented version is QRP when p or q ≡ 1 mod 4, and Jacobi-trivial
   for Blum N), and building an explicit equivalence gives square roots, which is FACT. A topological home for family Q.
2. **Tilings (D1–D5, D12; the user's question).**
   - Tilings of Z_N restate factoring.
   - Sparse probe sets (fewer than p elements) are provably blind to the cyclotomic divisors of x^N − 1 (Lam–Leung
     weights).
   - Fixed tiles give aggregates.
   - Tiling counts of compressed regions fall into a dichotomy (P3). Transfer-matrix counts are N-exponent gates
     (family M; codex's periodic-dimer/Wendt candidate tested: all 18 of 18 asymmetric splits are small-order
     events). Product-formula counts (MacMahon, Kasteleyn) are factorial ratios, so a fast mod-N counter would solve
     Strassen's factorial problem. For example, the (1, x, x) hexagon has C(2x, x) lozenge tilings, which splits N
     for x ∈ [p/2, min(p, q/2)).
3. **d-invariant collisions of L(N, a) are Fermat/Lehman (codex's candidate; K10; verified).** Every
   non-conjugate collision d_i = d_j splits N, and about 9–20·√N of them exist at n ≤ 16. But an exact collision is
   x² − y² = kN with k bounded by the lower recursion level: at a = 2 these are exactly Fermat's pairs, and at
   a ≈ N^{1/6} the bound reproduces Lehman's N^{1/3} (critic). So topology rediscovers difference-of-squares
   methods.
4. **Fixed-level quantum invariants (P1).** These are aggregates for a structural reason: the SL₂(Z) representation
   of a modular tensor category has a congruence kernel. This sharpens newfamily row 31.
5. **Combinatorial restatements.** Unitary Cayley graph = K_p ⊗ K_q, with clique number p; maximum sum-free density
   of Z/N gives the least prime factor ≡ 2 (mod 3); Bruck–Ryser for planes of order N is the project's open bit ε.

**Shortlist (§13).** No item survives as OPEN-WORTH-A-LOOK. The only item aimed at an actual cost is codex's
near-regular-matroid elimination for NFS matrices (X11): a constant-only (suffmap Tier C) target that we expect to
die on a cheap test, stated in §13.

Subjective probability that any row here leads to a record: < 0.5%.

---

## 1. Four cross-field principles (what kills most candidates)

These recur across the ten fields; the tables cite them as P1–P4.

- **P1 (fixed finite target ⇒ aggregate).** If N enters a construction only as an exponent, twist count or surgery
  coefficient, and the invariant is read through a fixed finite object (finite group, F_ℓ-homology, a TQFT at fixed
  level, a fixed tile's transfer matrix), the invariant depends on N mod m for a fixed m (geometry Prop AGG-H, Prop
  CYC). New instance here: for a modular tensor category the SL₂(Z) representation has congruence kernel of level
  ord(T) (Ng–Schauenburg, CMP 300 (2010)), so fixed-level Reshetikhin–Turaev invariants of lens spaces and torus
  bundles built from N are AGG up to a framing phase e(c·Φ(M)/24), where Φ is the Rademacher function of the CF of
  a/N (class E). KNOWN ingredients; consistent with newfamily row 31.
- **P2 (sampling and local estimation are archimedean).** Densities over Z/N that see the factors differ from the
  generic value by ≍ 1/p (φ(N)/N = 1 − (p+q−1)/N, eigenvalue masses (p−1)/N of the unitary Cayley graph, character
  sums of size √N among N terms). Monte Carlo, random-walk moments and sublinear spectrum estimation
  (Cohen-Steiner–Kong–Sohler–Valiant, KDD 2018: Wasserstein error ε from (1/ε)^{O(1/ε)} samples) resolve such
  terms only with N^{Ω(1)} samples. Family (A), EXP.
- **P3 (tiling and dimer counts of compressed regions).** A region described in poly(n) bits but with area ≍ N^{Ω(1)}
  has a tiling count computed either (a) by a fixed-width transfer matrix, which is a linear recurrence: then
  count(N) mod p is an N-exponent gate (leakmap class G / NE1, Lehmer–Pierce Res(x^N − 1, P), family M; Bach–Shallit,
  Math. Comp. 52 (1989)); or (b) by a product formula (MacMahon, Aztec, Kasteleyn on rectangles), which is a ratio of
  factorials/Pochhammer products, so evaluating it mod N is the factorial problem: Pollard–Strassen N^{1/4}, and a
  poly-size SLP would be suffmap's Tier-E Pochhammer target. Example (PROVED, Kummer): the lozenge count of the
  (1, x, x) hexagon is C(2x, x), and gcd(C(2x, x) mod N, N) = p for x ∈ [p/2, min(p, q/2)) (one base-p carry, no
  base-q carry). For x ≥ p a carry is not guaranteed (C(2p, p) ≡ 2 mod p; critic: p = 11, q = 37 fails for 6 of 13
  x). Check: p = 1009, q = 1201 split for x ∈ [505, 600], `checks.py`.
- **P4 (logic transfers hardness one way).** Every proof-complexity or bounded-arithmetic link we found has the form
  "property X of a proof system / provability of Y in a weak theory ⇒ factoring in P/poly (or P)", where X or Y is believed false
  precisely because of this implication (Krajíček–Pudlák; Bonet–Pitassi–Raz; Liu–Pass). Such links give no algorithm
  unless one proves X, which is harder than factoring.

---

## 2. Knot theory and low-dimensional topology

Prior coverage: geometry §5 and rows 22–32 (lens spaces, cyclic covers, Prop AGG-H/CYC/KUM, linking numbers OPEN-2),
newfamily rows 13, 31–35 (d-invariants, WRT, Kashaev, Habiro).

| # | Connection (object from N; feature) | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| K1 | **Homotopy type of lens spaces.** L(N, a) ≃ L(N, b) iff ab ≡ ±n² (mod N) (orientation-preserving: ab ≡ n²). Deciding homotopy equivalence of two triangulated lens spaces with H₁ = Z/N is therefore the quadratic residuosity problem mod N (oriented version: exactly QRP; unoriented: c ∈ ±□, which is QRP when p or q ≡ 1 mod 4 but Jacobi-trivial for Blum N); producing an explicit equivalence (its π₁-degree n) is a square root mod N | new to project (not in geometry) | **→Q** (≡_R QRP for decisions; construction FACT by Rabin) | Reduction QRP ≤ HOMOTOPY: given c with (c/N) = 1, pick a with small CF partial-quotient sum and b = c·a⁻¹r² (rerandomise r until b's CF is short), so both have poly(n)-size layered triangulations; L(N,a) ≃⁺ L(N,b) iff c is a square. The converse direction reads the linking form λ(g,g) from any triangulation (Smith form). A topological home for QRP, not a leak | Reidemeister 1935; Whitehead 1941; Yao–Knuth 1975 (mean partial-quotient sum ≍ (6/π²)ln²N, so rerandomisation succeeds by Markov); Jaco–Rubinstein (layered triangulations, size = partial-quotient sum) |
| K2 | Homeomorphism of lens spaces (q ≡ ±q′^{±1}); Reidemeister torsion, d-invariants, Casson–Walker of L(N, a) | geometry #23, newfamily #13 | AGG / class E | Euclid recursions along a/N; Casson–Walker is a Dedekind sum and its local part is FACT (newfamily Prop DSF) | Ozsváth–Szabó 2003; Brody 1960 |
| K3 | Fixed-level quantum invariants (WRT/RT, Turaev–Viro at fixed r) of L(N, a), of surgeries with coefficient N, of 2-bridge knots b(N, a) (Jones at a fixed root of unity) | newfamily #31; P1 | AGG (+ class-E framing phase) | P1: the mapping-class action of the CF word factors through SL₂(Z/ord T) (Ng–Schauenburg); the leftover phase is Rademacher Φ. (cx: Turaev–Viro on *arbitrary* triangulations is #P-hard already at r = 4, so "fixed level" is only cheap for these fixed-shape families) | Ng–Schauenburg 2010; Jeffrey 1992; Burton–Maria–Spreer 2018 |
| K4 | Quantum invariants at ζ_N (level growing with N): WRT of Seifert manifolds at e(1/N), Kashaev/colored Jones at e(1/N), Habiro values in Z[ζ_N] | newfamily #33–35 | class E or EXP | Seifert: Gauss + Dedekind sums (AGG + E); hyperbolic: O(N) recurrences; exact value has φ(N) coordinates and its trace over Q(ζ_N) is a divisor sum | Lawrence–Zagier 1999; Habiro 2004 |
| K5 | Representation counts #Hom(π₁(M), G_N) for a fixed 3-manifold M and N-dependent target: G = Z/N gives ∏ gcd(tors_i, N)·N^{b₁}; G = SL₂(Z/N) gives the product of F_p- and F_q-point counts of the character variety | geometry Prop AGG-H (fixed G only) | AGG (abelian) / FACT + EXP (SL₂) | abelian: gcds with fixed torsion numbers. SL₂(Z/N): CRT-multiplicative point count of a fixed variety over Z/N, i.e. geometry Prop MF (FACT), computable only by enumeration | — |
| K6 | Arithmetic topology: linking number of the "knots" p, q (Legendre (p/q)), Milnor/Rédei invariants, arithmetic Chern–Simons of Q(√N) with finite gauge group | geometry #29–30, OPEN-2 | →Q / CLASS | 4-rank of Cl⁺(Q(√N)) (Rédei–Reichardt); arithmetic CS needs the class field theory of a field of discriminant ≍ N (L[1/2]) | Morishita 2012; Kim 2015; Chung–Kim–Kim–Park–Yoo |
| K7 | Lens-space surgery realisation: is L(N, a) integer surgery on a knot in S³ (Greene's changemaker classification, Berge families) | not covered | →Q / RESTATE | Berge-type conditions are a ≡ ±k^{±2} (mod N) with k ranging over divisors of N ± 1 (Berge types I/II: N = ik ± 1, a factoring problem for a *different* number) or k² + k + 1 ≡ 0 (mod N): small-k cases are Coppersmith-findable, the generic ones are QR decisions whose two distinct witnesses would factor N; no witness is produced by the topology | Greene, Ann. Math. 177 (2013) |
| K8 | Cyclic branched covers Σ_N(K) of a fixed knot; Fox colourings of b(N, a) | geometry #25, #28 | AGG / RESTATE / EXP | Prop CYC; \|H₁\| = ∏_{d\|N}\|Res(Δ, Φ_d)\| has ≍ N digits; m-colourable iff gcd(m, N) > 1 | Fox; Riley 1990 |
| K9 | Polynomial-time 3-manifold algorithms (normal surfaces, S³/lens recognition in NP, FPT Turaev–Viro) applied to manifolds built from N | not covered | AGG / RESTATE | these decide topology of a poly(n)-size triangulation; when N is in the description, the topology is π₁ = Z/N-type data (K1–K3), and covers of degree d exist iff d \| N (geometry Prop COV: RESTATE) | Burton–Maria–Spreer 2018; Lackenby–Schleimer 2022 |
| K10 | **d-invariant collisions (cx ★).** d(L(N,a), i) is cheap for each label i (Ozsváth–Szabó recursion); a collision d_i = d_j with j ≢ a − 1 − i (mod N) gives N \| 2(i − j)(i + j − a + 1) and hence a split | cx (Jabuka–Robins–Wang, IMRN 2013); not in repo | **RESTATE: Fermat/Lehman** | Verified (`dinv.py`, a = ⌈√N⌉, 18 semiprimes, n = 12–16): ≈ 9–20·√N nontrivial collisions, *every* one splits N. Critic: d_i ≡ ā·x_i²/(4N) − 1/4 (mod 1), ā = a⁻¹ mod N, x_i = 2i + 1 − N − a, so a collision is x² ≡ y² (mod N); an exact collision is x² − y² = kN with k bounded by the range of the lower recursion term (≲ a²). At a = 2 the collisions are exactly Fermat's X² − Y² = ±N (checked at N = 143, 217, 221); a ≈ N^{1/6} reproduces Lehman's N^{1/3}. A topological rediscovery of difference-of-squares methods | Ozsváth–Szabó 2003; Jabuka–Robins–Wang 2013; Lehman 1974 |

## 3. Symplectic geometry and Floer theory

No prior coverage beyond geometry's lens-space Heegaard Floer rank (= N, AGG) and Prop CAT (periodic points).

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| S1 | HMS for the elliptic curve: Lagrangian of slope (a, N) on T² ↔ vector bundle of rank N, degree a; HF ranks = \|det\| intersection numbers; Massey products = theta functions with characteristics in (1/N)Z | none | RESTATE / AGG | Floer ranks are \|det\| of slope pairs (= N or gcd-type data); Atiyah: indecomposability ↔ gcd(N, a) = 1; theta values at a/N are archimedean functions of a/N with no divisor structure | Polishchuk–Zaslow 1998; Atiyah 1957 |
| S2 | Floer homology of iterates HF(φ^N) of a fixed symplectomorphism (= Fix(φ^N) for pseudo-Anosov / toral maps); Smith-type inequalities for Z/p-equivariant Floer | geometry #37–38 (Prop CAT) | AGG mod m / FACT + EXP | ranks ≍ λ^N are exponential integers; mod a fixed ℓ periodic (P1); primitive part is a Möbius sum (FACT, one query); Smith inequality needs the Z/p sub-action, i.e. φ^q (RESTATE) | Seidel–Smith 2010; Hendricks 2012 |
| S3 | ECH capacities / ellipsoid embeddings E(1, N) → B(λ), McDuff–Schlenk staircases; lattice-point counts in triangles of slope a/N | none | →A / class E | above the Fibonacci staircase (a ≥ (17/6)²) the capacity is √a (volume); ECH capacities are lattice counts in triangles: Ehrhart/Dedekind-sum data along CF (E) | McDuff–Schlenk, Ann. Math. 175 (2012) |
| S4 | Lagrangian Floer / Fukaya categories with coefficients in Z/N (Novikov over Z/N) | none | →PT | any finite chain-level computation over Z/N is a generic-ring algorithm: factorbits Theorem PT (no info without splitting) | factorbits §3.2 |
| S5 | Symplectic Khovanov / Heegaard Floer of torus knots T(p, q), Brieskorn Σ(p, q, r) | geometry #27 | RESTATE | the object needs p, q to be written down | — |

## 4. Homotopy theory

No prior coverage (no hits for TMF, Morava, chromatic, image of J).

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| H1 | **Image of J at stem 2N − 3.** \|im J\| in π^s_{4k−1} is the denominator of B_{2k}/4k; with 2k = N − 1, p divides it iff (p − 1) \| (N − 1) iff (p − 1) \| (q − 1) | none | →M (rare Fermat-liar event) | von Staudt–Clausen: the primes are those ℓ with (ℓ−1) \| (N−1); for p \| N this is the Korselt-type condition (p−1) \| (q−1), a group-order event of density o(1). Computing the denominator also needs the divisors of N − 1 (factoring N − 1) | Adams, J(X) IV (1966); Quillen; von Staudt–Clausen |
| H2 | Agoh–Giuga sums N·B_{N−1} mod N = Σ_{a<N} a^{N−1} mod N | none | →M, EXP | for p \| N the sum is q·Σ_{a<p} a^{N−1} ≡ −q or 0 (mod p) according to (p−1) \| (q−1): same Korselt event; no fast evaluator (Faulhaber needs B_k up to N) | Agoh 1995; Giuga 1950 |
| H3 | TMF with level structure TMF_0(N), TMF_1(N) (defined after inverting N); rational homotopy = modular forms of level N | geometry #19a (dimensions) | FACT + EXP | ranks are dim M_k(Γ_1(N)), whose differences give ψ(N) or J₂(N) (FACT); building the étale cover of the moduli stack has degree ψ(N) | Hill–Lawson, Invent. Math. 203 (2016) |
| H4 | Chromatic height: formal group of E/(Z/N) has height 2 at a supersingular prime and 1 at an ordinary one; [N]-series leading degree p^h | none | →M / EXP | detecting height = Hasse invariant mod p (exponent p−1 coefficient, cf. newfamily #54 C(2N,N)); the [N]-series' informative term has degree ≥ p. Over Z/N the formal-group computation is ECM's group law in disguise | Lubin–Tate; Silverman AEC IV |
| H5 | Morava K-theory / HKR characters of BZ/N; K-theory of Z/N (Quillen: K_{2i−1}(Z/N) = Z/(p^i−1) ⊕ Z/(q^i−1)) | conn_alg (K-theory) | RESTATE / FACT | K(n) at the prime ℓ sees only the ℓ-part of N (must pick ℓ = p); \|K₁(Z/N)\| = φ(N), \|K₃\| = (p²−1)(q²−1): FACT | Hopkins–Kuhn–Ravenel 2000; Quillen 1972; cx (\|K_{2j−1}\| = (p^j−1)(q^j−1) gives p^j + q^j) |
| H6 | Equivariant stable homotopy of Z/N: Burnside ring A(Z/N) (rank d(N)), tom Dieck splitting indexed by subgroups, R(Z/N) ⊗ Q = ∏_{d\|N} Q(ζ_d) | none | ω-tier / RESTATE | ranks count divisors (4 for every semiprime); idempotents of R(Z/N) ⊗ Q are the cyclotomic factors of x^N − 1 (degree-N object, and the factor Φ_p needs p) | tom Dieck; Dress 1969 |
| H7 | Adams operations ψ^N = ψ^p ψ^q; e-invariants; vector fields on spheres (Radon–Hurwitz) | none | AGG | ψ^k on K(S^{2m}) is multiplication by k^m (an integer power of N); vector-field numbers depend on v₂ only | Adams 1962 |

## 5. Category theory and higher categories

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| C1 | Krull–Schmidt: Z/N ≅ Z/p ⊕ Z/q; idempotents of End(Z/N); K₀ / G₀ of Z/N-modules ≅ Z² | geometry #40 | RESTATE / ω | a decomposition is an idempotent is a factor; ranks give ω only | — |
| C2 | **Prime decomposition of modular tensor categories.** The pointed MTC C(Z/N, x²/N) factors (Müger) as C(Z/p) ⊠ C(Z/q) | none | RESTATE | finding a modular subcategory = finding a subgroup of Z/N | Müger, Proc. LMS 87 (2003) |
| C3 | Gauss sums / higher central charges ξ_n of C(Z/N, x²/N); Witt class in W(Q/Z) = ⊕_ℓ W(F_ℓ) | none | AGG / →Q | for gcd(n, N) = 1 these are Jacobi-symbol Gauss sums (AGG, Milgram); the ℓ-components of the Witt class are ⟨q⟩ ∈ W(F_p), ⟨p⟩ ∈ W(F_q), i.e. the Legendre symbols (q/p), (p/q) (= geometry OPEN-2); n with gcd(n, N) = p needs p | Ng–Rowell–Wang–Zhang 2022; Milgram |
| C4 | Decategorifications equal to σ(N), τ(N), d(N): index-N sublattices (σ), Hecke T_N traces (τ(N) = τ(p)τ(q)), Burnside counts | geometry #10, ecount (Bach–Charles) | FACT / OPEN-known | σ(N) FACT; τ(N) together with τ(N²) factors N (Bach–Charles, arXiv 0708.1192), τ(N) alone OPEN (tracedata's a_N question); Edixhoven–Couveignes needs n factored | Bach–Charles 2007; Edixhoven–Couveignes 2011 |
| C5 | Magnitude / Euler characteristic of categories and graphs built from N (divisor poset, unitary Cayley graph) | none | FACT + EXP / →A | magnitude of the unitary Cayley graph is N/(1 + t·φ(N) + t²(N − 1 − φ(N))) (degree φ(N)): FACT, and estimating φ(N) by sampling is P2 | Leinster 2008, 2013 |
| C6 | Categorified "Spec Z/N": topoi, F₁-geometry, Connes–Consani arithmetic site evaluated at N | none | RESTATE | no constructible finite object; points/components are the primes | Connes–Consani 2014 |
| C7 | Subobject and weighted counts (cx): #{submodules of (Z/N)²} = (p + 3)(q + 3) = N + 3(p + q) + 9; Σ over subgroups of C_N of their order = σ(N); splitting the unit of Perf(Z/N) | cx | FACT | compact definitions of factoring-strength counts; evaluation enumerates or uses CRT. An idempotent splitting of the tensor unit is e² = e, i.e. gcd(e, N) | cx |

## 6. Extremal combinatorics, Ramsey theory, matroids

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| X1 | **Matroid characteristic sets.** Von Staudt / Mnëv constructions give, from the binary expansion of N, a poly(n)-size matroid M_N whose realisation space is Spec Z/N: M_N is representable over a field F iff char F \| N | rankdegen (unrelated matroid use) | RESTATE (FACT) | the combinatorial matroid is free data; its characteristic set {p, q} is exactly the factorisation, and its Tutte polynomial etc. are N-generic combinatorics | Mnëv 1988; Kahn 1982; Sturmfels |
| X2 | **Max sum-free density of Z/N** = 1/3 + 1/(3p′) with p′ the least prime factor ≡ 2 (mod 3); 1/3 − 1/(3N) when all prime factors ≡ 1 (mod 3) | none | RESTATE + EXP | an extremal number of a size-N optimisation; reveals p only if p ≡ 2 (mod 3), and computing it is the size-N problem (check Z/15: 6 = 15(1/3 + 1/15)) | Diananda–Yap 1969; Green–Ruzsa 2005 |
| X3 | Kneser / Cauchy–Davenport in Z/N: \|A + B\| < \|A\| + \|B\| − 1 forces a nontrivial period subgroup | none | RESTATE + EXP | any deficient pair is a union of cosets of a subgroup of order ≥ p, so witnesses have size ≥ p; below p, Z/N behaves like Z/(prime) | Kneser 1953 |
| X4 | White-box Ramsey: find a clique/independent set of size n/2 in a succinct graph on 2^n vertices | suffmap (PWPP) | KNOWN (FACT-hard direction) | hard if collision-resistant hashing exists (factoring-based CRHF via claw-free permutations), so a Ramsey solver would factor; no algorithmic content | Komargodski–Naor–Yogev, FOCS 2017; Goldwasser–Micali–Rivest 1988 |
| X5 | Chevalley–Warning / Combinatorial Nullstellensatz / EGZ over Z/N; PPA_q | suffmap, detGRH (Jeřábek) | KNOWN | existence arguments put FACT in PPA/PWPP (Jeřábek 2016); Chevalley–Warning is PPA_p-complete; no algorithm. Counting zeros over Z/N is CRT-multiplicative (geometry Prop MF) | Göös–Kamath–Sotiraki–Zampetakis, CCC 2020; Jeřábek, JCSS 82 (2016) |
| X6 | Designs: projective planes of order N (Bruck–Ryser), MacNeish/orthogonal Latin squares, Paley-type constructions on Z/N | none | →Q | Bruck–Ryser excludes N ≡ 1, 2 (mod 4) not a sum of two squares: for semiprimes this is ε = (−1/p) (r₂(N) > 0), the project's open Q bit; constructions need prime-power orders | Bruck–Ryser 1949 |
| X7 | Lattice tilings by crosses ↔ splittings of Z/N by {±1, …, ±k} ↔ subgroups of (Z/N)^* | none | →M / FACT | a splitting set is a union of cosets of a subgroup of (Z/N)^* of controlled index; existence depends on ord_p of small integers (M) and construction needs the CRT structure | Stein–Szabó, *Algebra and Tiling* (1994) |
| X8 | Jacobi-coloured Ramsey (cx): colour {x, y} by ((x − y)/N); seek homogeneous sets with a non-unit difference | cx | →Q / RESTATE | homogeneity forces no non-unit difference; a difference divisible by p would already split | cx; KNY 2017 |
| X9 | Binomial gcds (cx): gcd(C(N, k), N) is N for 1 ≤ k < p and q at k = p | cx; newfamily #54 | RESTATE / →⟨factorial barrier⟩ | the cancellation occurs exactly at the factor scale (Kummer); fast evaluation is the factorial problem (P3) | Kummer 1852 |
| X10 | Matroids over Z/N (Fink–Moci) from a public matrix: local matroids over F_p, F_q differ | cx | →PT | a minor vanishing mod exactly one prime is a gcd split; finding such a minor in a succinct huge matrix is the whole problem (Theorem PT for small matrices) | Fink–Moci, JEMS 18 (2016) |
| X11 | **Near-regular binary matroids in NFS relation matrices (cx ★).** Delete few columns so the column matroid is regular, then eliminate via Seymour's decomposition (graphic/cographic pieces) instead of Block Wiedemann | cx; suffmap Tier C (sparse F₂ linear algebra) | →⟨suffmap Tier C⟩, implausible | target is the sparse linear-algebra exponent ω′ (ω′ = 2 is today's cost, MNFS 1.9019; the floor with free linear algebra is disputed, LSTT 1.8689 vs heurnfs 1.8595; GNFS stays 1.923): constant only. Obstacles: filtered NFS matrices are random-like sparse binary matrices, and random binary matroids of growing rank contain Fano minors (we expect many disjoint ones; UNVERIFIED for the NFS law); regular-matroid recognition/decomposition (Truemper) is cubic, above Block Wiedemann's B^{2+o(1)}; the graphic special case (weight ≤ 2 columns, large-prime cliques) is already eliminated by CADO's merge | Seymour, JCTB 28 (1980); Truemper 1990 |
| X12 | Sum-product in Z/N (Bourgain): sets with small \|A + A\| + \|A·A\| must be near a subring, and the proper subrings of Z/N are the factor ideals (critic) | suffmap (sum-product unrelated) | RESTATE + EXP | the obstruction set is ≍ a coset of pZ/N, so it has size ≥ q; finding one is finding a factor | Bourgain, Katz–Tao; Bourgain 2008 (Z/N) |

## 7. Spectral graph theory and expanders

Prior coverage: geometry row 20 and §7 (coset-graph walks on P¹(Z/N), Weyl law), suffmap/navigate (Cayley-graph
navigation is sufficient and FACT ∧ DL-hard), genusleak Lemma GIRTH.

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| G1 | **Unitary Cayley graph** X_N = Cay(Z/N, (Z/N)^*) = K_p ⊗ K_q: edges poly-time from N; spectrum {φ(N)¹, −(q−1)^{p−1}, −(p−1)^{q−1}, 1^{φ(N)}}; clique and chromatic number = p | none in repo | RESTATE + →A | every spectral or extremal quantity reveals p, but the graph has N vertices; closed-walk/return statistics differ from uniform by O(1/p) relative terms (P2); clique number is an NP-type optimum whose value is the answer (checked N = 15, 35, 77) | Klotz–Sander, EJC 14 (2007) |
| G2 | Jacobi-symbol Cayley graph Cay(Z/N, {x : (x/N) = 1}); small-subgraph densities | newfamily #47 | →A / EXP | eigenvalues at gcd(k, N) = 1 are AGG Gauss sums; subgraph counts are character sums such as Σ_{x mod N}((x³+ax+b)/N) = a_p·a_q (verified on 2 semiprimes, `checks.py`), i.e. tracedata's open a_N, but as an N-term sum whose estimation needs ≍ N² samples (P2) | Jacobsthal; tracedata |
| G3 | Random walks on (Z/N)^*, SL₂(Z/N) with fixed generators: returns, mixing, girth | genusleak GIRTH, geometry row 20 | →M / KNOWN | a word ≡ I mod p but not mod q is a relation in SL₂(F_p); short words have entries < p (girth ≍ log p), long ones are birthday/rho searches on group orders: Pollard rho, Williams p+1 | Pollard 1975; Williams 1982 |
| G4 | LPS / Morgenstern-type Cayley graphs on PSL₂(Z/N) (Ramanujan only for prime modulus; expanders over Z/N by strong approximation) (built from N using a root of x²+y²+1 ≡ 0, Pollard–Schnorr); path finding | navigate, smoothbuild | FACT-hard direction | for prime modulus path finding is poly-time (Petit–Lauter–Quisquater 2008) because it solves norm equations with square roots mod q; at composite N those square roots are Rabin-hard, so navigation is at least FACT-hard (UNVERIFIED in print); spectrum multiplicities (dim² of ρ_p ⊗ ρ_q) need the N³-vertex graph | LPS 1988; Petit–Lauter–Quisquater, SCN 2008 |
| G5 | Brandt / Pizer graphs of level N (supersingular isogeny graphs with Γ_0(N) structure): closed-walk counts tr B(ℓ^k) via Eichler's trace formula | geometry #16–17, OPEN-1; Kirschmer–Voight | →Q (odd k) / FACT + EXP (even k) | tr B(ℓ^k) contains the main term (r−1)ψ(N)/12 only when ℓ^k is a square; otherwise class numbers × ∏_{r\|N}(1 + (D/r)): the X_0(N)(F_ℓ) vs F_{ℓ²} dichotomy in spectral form; the graph has ≍ ψ(N) vertices | Eichler; Pizer 1980; Kirschmer–Voight 2010 |
| G6 | Uniform expansion of SL₂(Z/N) for squarefree N (Bourgain–Gamburd–Sarnak affine sieve) | none | KNOWN-irrelevant | gives *existence* of almost-prime values in thin orbits (the opposite of smooth values); suffmap: existence results do not move the cost | BGS, Invent. Math. 179 (2010) |
| G7 | Spectral replacement for the linear-algebra step (Laplacian solvers, spectral sparsification) | suffmap Tier C | KNOWN constant-only | the NFS matrix is over F₂, not SDD over R; even ω < 2 linear algebra changes only the MNFS constant | Spielman–Teng 2014; suffmap §2 |
| G8 | Sandpile (critical) groups of circulants / wheels built from N; spanning-tree counts | geometry #46 | →M (class G) / EXP | closed forms are Lucas-type sequences at index N (e.g. wheel: L_{2N} − 2), so mod p they are N-exponent gates; Smith forms of N×N Laplacians are EXP | Biggs 1999 |

## 8. Tropical geometry

Prior coverage: smoothbuild W3 (kill: smoothness is a valuation condition, the tie to N is archimedean; tropicalised
problem = log-lattice CVP, schnorrcvp), geometry #45.

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| T1 | Tropicalisation of xy = N w.r.t. v_p (the line X + Y = v_p(N)) or the archimedean amoeba (X + Y = log N) | smoothbuild W3 | RESTATE | the p-adic tropicalisation needs p; the archimedean one carries only log N; integral points on the amoeba are the divisors | — |
| T2 | Berkovich spectrum M(Z): the function x ↦ log\|N\|_x is nonzero exactly on the branches of p \| N | none | RESTATE | evaluating on a branch needs its prime | Berkovich 1990 |
| T3 | Newton polytopes in Coppersmith/lattice methods (Blömer–May, Jochemsz–May strategies) | coverage38, resbarrier | KNOWN | tropical/polytope optimisation of lattice shapes is standard and stays at N^{1/4}-type exponents | Blömer–May, EUROCRYPT 2005; Jochemsz–May, ASIACRYPT 2006 |
| T4 | Tropical (min-plus) algebra as a sieve or linear-algebra primitive; Newton-polygon root properties in NFS polynomial selection | heurnfs, smoothbuild W3 | KNOWN / →⟨Theorem SS⟩ | sieving is already additive accumulation of log p over residue classes; v_p(F(a, b)) depends on residues, not on a tropical polynomial in (a, b); Murphy's α is polynomial selection (constant only) | Murphy 1999 |

## 9. Discrete geometry, sphere packing, tilings (treated carefully)

Prior coverage: geometry rows 48–50 (Sturmian/cut-and-project = E/AGG; tiling decision problems = RESTATE;
Deninger), geometry #8–9 (polytopes), newfamily #49–50 (Ehrhart/Barvinok). The user asked about tilings, so this
field has the most rows.

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| D1 | **Tilings of Z_N.** A ⊕ B = Z_N | geometry #49 (generic) | RESTATE (PROVED, trivial) | \|A\|·\|B\| = N, so every tiling of Z_pq by a tile with 1 < \|A\| < N has \|A\| ∈ {p, q}: writing down any nontrivial tile *is* a factor. Hence "count/decide tilings of Z_N by a given tile" cannot even be posed without one | — |
| D2 | Cyclotomic tiling criteria: A ⊕ B = Z_N iff \|A\|·\|B\| = N and every Φ_d, 1 < d \| N, divides A(x) or B(x) (elementary); Coven–Meyerowitz (T1)/(T2) settle tiling of Z when \|A\| has ≤ 2 prime factors; Łaba–Londner proved (T2) for three odd primes | none | RESTATE + PROVED blindness of probe sets | as tiles have \|A\| ∈ {1, p, q, N}, the content is about arbitrary *probe* sets: for a 0/1 polynomial A with fewer than p terms, gcd(A(x), x^N − 1) = 1 up to Φ_1, because a vanishing sum of d-th roots of unity with nonnegative coefficients has weight in pℕ + qℕ (Lam–Leung, J. Algebra 224 (2000)). Sparse (poly(n)-term) sets never see the cyclotomic divisors of x^N − 1; AP-union sets see them only through gcd(k_i, N) | Coven–Meyerowitz, J. Algebra 212 (1999); Łaba–Londner, Invent. Math. 232 (2023); Lam–Leung 2000 |
| D3 | Hajós/Sands: which cyclic groups are "good" (every factorisation has a periodic factor); Fuglede/spectral sets in Z_N | none | RESTATE (ω-pattern) | goodness depends only on the exponent pattern of N (Z_N good iff N \| p^a q, p²q², pqr, p²qr or pqrs); every Z_{pq} is good; Fuglede is known in Z_{pq} | Hajós 1941; Sands 1962; Szabó–Sands, *Factoring Groups into Subsets* (2009) |
| D4 | Fixed tile, tilings of Z_N or of the N-periodic strip/cylinder: counts = tr(T^N) for a fixed transfer matrix | geometry Prop AGG-H | AGG / →M | P1 and P3(a): an exponent gate on a fixed matrix | Newman 1977 (periodicity) |
| D5 | **Kasteleyn / dimer and lozenge counts of compressed regions** (area ≍ N^{Ω(1)}, poly(n)-bit description) | none | →M (transfer) or →⟨Pochhammer SLP barrier⟩ (product formulas) | P3: transfer matrices give exponent gates; product formulas are factorial ratios, so a fast mod-N evaluator would be Strassen's factorial problem (deterministic N^{1/4} via Pollard–Strassen; poly-size SLPs = suffmap Tier E). Lozenge count of the (1, x, x) hexagon = C(2x, x) splits N for x ∈ [p/2, q/2) | Kasteleyn 1961; MacMahon; Strassen 1976; suffmap §2 |
| D6 | Substitution tilings / automatic sequences with N-dependent inflation (constant length N; base-N odometer Z_N ≅ Z_p × Z_q) | none | RESTATE / AGG / EXP | the N-adic odometer splits by CRT (idempotents = factors); Dekking height is defined by gcds with N; Cobham: base-N and base-p automaticity are incompatible except for periodic sequences; learning the p-automatic component (e.g. C(2n, n) mod p) needs ≍ p states | Dekking 1978; Cobham 1969; Christol 1979 |
| D7 | Viazovska-type modular-form constructions at level N (magic functions, Fourier interpolation coefficients, theta series of lattices of determinant N) | geometry #19a | FACT + EXP / →Q | a basis of M_k(Γ_0(N)) has dimension ≍ ψ(N)/12 (dimension itself FACT); theta series of x² + Ny² give genus-theoretic r_Q(m) (Q); Radchenko–Viazovska coefficients are weakly holomorphic (Rademacher series, ≍ √N digits, archimedean) | Viazovska 2017; Radchenko–Viazovska 2019 |
| D8 | Squared rectangles / tilings of an a × b rectangle by squares or by 1 × k bars (de Bruijn–Klarner: k \| a or k \| b) | geometry #49 | RESTATE / class E | bar tilings test divisibility of chosen k; square tilings follow the Euclidean algorithm of a/b (CF, class E) | de Bruijn 1969; Kenyon 1996 |
| D9 | Greenfeld–Tao aperiodic / undecidable translational tilings (encode p-adically structured functions) | none | unrelated | the encoding goes arithmetic → tiling (they force p-adic structure for chosen primes); no tiling problem with N as data is shown to be easier than factoring | Greenfeld–Tao, Ann. Math. 200 (2024) |
| D10 | Lattice/packing invariants of 2-D lattices from N (x ≡ ay mod N; forms of discriminant −4N): minimal vectors, kissing, densest member | geometry #11 | class E / CLASS | Gauss reduction = CF of a/N (E); reduced-form counts are h(−4N) (CLASS, L[1/2]) | — |
| D11 | E₈ theta coefficient (cx): r_{E₈}(2N) = 240σ₃(N), and σ₃(N) − N³ − 1 = p³ + q³ | cx | FACT | Viazovska's inequalities certify density; they give no fast isolated coefficient (same as r₄, analytic) | Viazovska 2017; cx |
| D12 | **Periodic dimers / generalised Wendt resultants (cx ★).** D_N(c) = Res(X^N − 1, (X + c)^N − 1) is a twisted periodic honeycomb determinant, and D_pq(c) ≡ D_q(c) (mod p) | cx | →M (tested) + EXP | p \| D_q(c) iff c ∈ μ_q − μ_q inside F_{p^k}, k = ord_q(p): for large k this has probability ≈ q²/p^k. Test (`wendt.py`, all 153 pairs 11 ≤ p < q < 80, c = 1, 2, 3): 18 asymmetric splits, all with min(ord_q p, ord_p q) ≤ 3, while only 21 of the 153 pairs have such a small order (critic), i.e. an order event (NE1 class G). Evaluation mod N is a degree-N resultant | Helou, Math. Comp. 66 (1997); cx |
| D13 | Tijdeman's dilation theorem: if A tiles Z_N (or Z), so does tA for gcd(t, \|A\|) = 1 (critic) | none | AGG | unit dilations act on tiles without changing the tiling structure, so the tiling question is invariant under (Z/N)^*-type rescaling; no factor information | Tijdeman 1995 |
| D14 | Kolountzakis–Matolcsi algorithms for translational tiling (cx, critic) | none | RESTATE / EXP | polynomial in the diameter of the tile, not in its bit length; for tiles of Z_N the diameter is ≍ N | Kolountzakis–Matolcsi, arXiv 0810.4338 |

**Tilings bottom line.** Tilings of Z_N are factorisations of Z_N, so they restate factoring (D1); sparse probe sets are
provably blind to the cyclotomic structure that the Coven–Meyerowitz conditions use (D2, Lam–Leung); fixed tiles give
aggregates (D4); and counts for large, compressed regions are either exponent gates or factorial ratios (D5, P3). The
last point is the only one with algorithmic content, and it lands on the known factorial/Pochhammer barrier: a fast
exact tiling counter for hexagons of side ≍ √N modulo N would be Strassen's problem solved.

## 10. Logic and proof complexity

Prior coverage: suffmap §2 (Jeřábek PPA/PWPP; Krajíček–Pudlák and Bonet–Pitassi–Raz feasible interpolation),
detGRH §8, sihard/floorhard ETH errata.

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| L1 | Krajíček–Pudlák: an NP pair provably disjoint in S¹₂ is P/poly-inseparable if RSA is secure; EF has no feasible interpolation unless RSA is broken; EF plus "primes are not composite" tautologies has no feasible interpolation unless factoring and DL are in P/poly | suffmap §2 | KNOWN (P4) | the positive properties would be (nonuniform) factoring algorithms, which is why they are believed false | Krajíček–Pudlák, LCC'94 LNCS 960 / Inf. Comput. 140 (1998) |
| L2 | Feasible interpolation for TC⁰-Frege fails if factoring Blum integers is hard for poly-size circuits (via Diffie–Hellman modulo a Blum integer) | suffmap §2 | KNOWN (P4) | a positive interpolation theorem would put Blum-integer factoring in P/poly; believed false | Bonet–Pitassi–Raz, SIAM J. Comput. 29 (2000) |
| L3 | WPHP and TFNP: FACT ∈ PPA ∩ PWPP (randomized; deterministic under GRH_q) | suffmap, detGRH | KNOWN | membership, not an algorithm; PPA/PWPP are believed hard | Jeřábek, JCSS 82 (2016) |
| L4 | Proof-complexity lower bounds (cx adds: short refutations, even if they exist, need not be findable; the composite instance is satisfiable) for "N is prime" multiplication-circuit CNFs (Resolution, Cutting Planes) | none located | KNOWN-irrelevant (literature UNVERIFIED) | lower bounds for refuting prime instances say nothing about finding factors of composite ones; CDCL factoring is empirically exponential | (search found no specific theorem) |
| L5 | **ETH refutation (OpenAI preprint, Oct 2026: #k-SAT in 2^{O(m log log m/√log m)})** applied to factoring circuits | sihard, floorhard errata | no change (PROVED, routine) | a factoring CNF has m ≥ n variables, m = Θ(n log n) to Θ(n²) by multiplier; 2^{m/√log m} ≥ 2^{n√log n·(1+o(1))} for m ≍ n log n, worse than trial division 2^{n/2}. Beating NFS via SAT would need 2^{m^{1/3+o(1)}} on these instances, which ETH⁻ does not address either way | errata in notes/sihard, floorhard |
| L6 | Meta-complexity: average-case time-bounded Kolmogorov complexity ⇔ OWFs, so an average-case K^t algorithm inverts multiplication on average | none | KNOWN (P4) | the direction is K^t-easy ⇒ factoring-easy; no evidence K^t is easier | Liu–Pass, FOCS 2020 |
| L7 | Decision procedures: Skolem arithmetic (ℕ, ×), existential Presburger with divisibility (∃x: x \| N ∧ 1 < x < N) | none | RESTATE | Skolem arithmetic works in the prime-exponent representation (factorbits: addition there ≡ FACT); the divisibility sentence's witness is a factor and its truth value is AKS | Lechner–Ouaknine–Worrell, LICS 2015 |

## 11. Cellular automata and computation models

Prior coverage: prior_art_map (BSS), slpcirc/suffmap (factorial SLPs), factorbits (prime-exponent representation).

| # | Connection | Coverage | Verdict | Why / where it fails | Refs |
|---|---|---|---|---|---|
| A1 | Unit-cost arithmetic RAM: factoring in O(log N) operations; real RAM with floor solves PSPACE | prior_art_map, suffmap Tier E | KNOWN | the honest restriction is bit cost (or constant-free SLP size), where it becomes the Pochhammer/τ-conjecture question | Shamir, IPL 8 (1979); Schönhage 1979; Hartmanis–Simon 1974 |
| A2 | **1-D CA preimage (Garden-of-Eden) search through a real-time multiplier array** | none | EXP for the de Bruijn method (light-cone argument) | predecessors for one step of a 1-D CA are found in poly time via the de Bruijn graph, but multiplication needs T ≥ n/2 steps (output bit n depends on all input bits), and T steps form a radius-T CA whose de Bruijn graph has \|Σ\|^{2T} = 2^{Θ(n)} states: no better than 2^{cn}. This bounds only the de Bruijn method, not preimage search in general; Atrubin's array has serial real-time I/O, so the model is "a 1-D CA run T steps from an initial configuration encoding x, y". In 2-D the preimage problem is NP-complete (Durand, TCS 134 (1994)) | Atrubin, IEEE TEC 14 (1965); Sutner 1995 |
| A3 | Reversible CA / FRACTRAN implementing x ↦ ax mod N or prime-register machines | factorbits | AGG / FACT | FRACTRAN reads only the exponents of the finitely many primes in its fractions (AGG on the large primes); periods of mod-N dynamics are orders (FACT, Miller) | Conway 1987 |
| A4 | DNA / tile self-assembly factoring | none | KNOWN, exponential | nondeterministic assembly; success probability 2^{−Θ(n)} per assembly, i.e. exponentially many molecules | Brun, TCS 395 (2008) |
| A5 | Analog / memcomputing / dynamical-system "solvers" | none | unsupported | claims of polynomial-time factoring hide exponential precision or time in the ODE simulation | Traversa–Di Ventra, Chaos 27 (2017) |
| A6 | Rule 110 / universal CA | none | RESTATE | universality makes CA a computation model with the same complexity (up to polynomial factors) as Turing machines | Cook 2004 |
| A7 | Special-purpose NFS hardware cost models: Bernstein's circuits (area-time), the LSTT 2002 mesh for linear algebra, TWIRL sieving (critic) | heurnfs, suffmap | KNOWN, constant only | changing the cost model (area × time) changes the L[1/3] constant, not the exponent; the most relevant computation-model item, already priced in the project's constant tables | Bernstein 2001; Lenstra–Shamir–Tomlinson–Tromer, ASIACRYPT 2002; Shamir–Tromer, CRYPTO 2003 |

---

## 12. GPT-6.1-sol blind round (cx), merged

`codex/conn_top_1.out`, launched before our analysis (effort xhigh). Its overall verdict matches ours ("no candidate
justifies beating the stated records"); of its 43 candidates it marks three OPEN-WORTH-A-LOOK, weakly (adjudicated below). Its additions, credited in the tables:
- K10 (d-invariant collisions);
- the Turaev–Viro #P-hardness caveat on K3;
- C7 (subobject counts; splitting the unit of Perf);
- X8–X11 (Jacobi-coloured Ramsey, binomial gcds, matroids over rings, near-regular NFS matroids);
- D11 (E₈ theta σ₃);
- D12 (Wendt/periodic dimers);
- K-theory p^j + q^j;
- the proof-search point in L4.

Disagreements, adjudicated:
- **D12 (codex: OPEN, weakly).** We tested it, and it is an order event. On 153 prime pairs × 3 shifts, all 18
  asymmetric splits have min(ord_q p, ord_p q) ≤ 3, as predicted by the μ_q − μ_q argument. Verdict →M.
- **K10 (codex: OPEN).** Confirmed numerically (every collision splits). v0 read it as a class-E remainder; the critic
  showed an exact collision is x² − y² = kN with bounded k (Fermat at a = 2, Lehman for a ≈ N^{1/6}). Verdict:
  RESTATE (difference-of-squares methods); dropped from the shortlist.
- **X11 (codex: best prospect).** Kept as a constant-only Tier-C target with stated obstacles (Fano minors in
  random binary matroids; cubic decomposition).
- **D-tilings.** Codex makes the same |A| ∈ {p, q} point and adds that the Coven–Meyerowitz tests of
  Kolountzakis–Matolcsi are polynomial in the diameter, not in the bit length. We add the Lam–Leung blindness lemma
  (D2) and the P3 dichotomy.

## 13. Ranked shortlist

**No candidate survives as OPEN-WORTH-A-LOOK.** v0 listed two "least-dead" items. One of them, K10, is
Fermat/Lehman (critic), which leaves one item aimed at an actual cost:
1. **Near-regular structure of NFS relation matrices (X11, cx).** At best a constant (suffmap Tier C), and expected
   to die. *Kill test:* take CADO's post-merge residual matrix from a real run (cadobench data). Merge has already
   eliminated the graphic, weight ≤ 2 part. Search it for F₇ / F₇* minors after deleting the heaviest B^{1−δ}
   columns. If they are abundant (expected), the idea is dead. Otherwise ask for a decomposition in B^{2−ε} bit
   operations; Truemper's recognition is cubic.

The OPEN items this survey touches are the project's existing ones, and it adds no evidence for or against them:
- geometry OPEN-1: Brandt/Pizer odd-length walk counts, i.e. small-discriminant QR (G5);
- geometry OPEN-2: the linking number (p/q) (K6/C3);
- tracedata's a_N (G2, C4).

## 14. Critic

Independent Opus critic (read-only; checks in its own scratch area), run on v0. **Verdict: ACCEPT WITH FIXES. No OPEN
item is mislabelled alive.** All fixes are applied in v1.

MAJOR:
1. **P3 hexagon claim.** "gcd = p for all x ∈ [p/2, q/2)" was false once x ≥ p (C(2p, p) ≡ 2 mod p). Corrected to
   [p/2, min(p, q/2)).
2. **K1.** The unoriented problem is Jacobi-trivial for Blum N. Now stated: oriented = exactly QRP; unoriented = QRP
   only when p or q ≡ 1 mod 4. The rerandomisation is sound; Yao–Knuth and Jaco–Rubinstein are now cited.
   Przytycki–Yasuhara removed (it is a homeomorphism proof).
3. **K10 is Fermat/Lehman.** The fractional-part formula is corrected. Verdict changed to RESTATE and the item
   removed from the shortlist.

MINOR:
- X11 payoff wording (ω′, LSTT vs heurnfs floor) and kill test (post-merge matrix).
- §0/L5: the ETH argument now needs m = Ω(n log n).
- L1: Krajíček–Pudlák restated from the abstract.
- L2: Bonet–Pitassi–Raz assumption named (Blum-integer factoring, P/poly); P4 says P/poly.
- H3: Hill–Lawson volume is 203.
- D2: the criterion is the elementary one; Lam–Leung applies to probe sets, not tiles.
- A2: scope limited to the de Bruijn method; Durand 1994 cited.
- K7: Berge types use divisors of N ± 1.
- G4: Ramanujan only for prime modulus.
- D12: base rate 21/153 added.
- Missing rows added: A7 (NFS hardware cost models), X12 (sum-product in Z/N), D13 (Tijdeman), D14
  (Kolountzakis–Matolcsi).

Checked correct by the critic:
- H1, H2, H5;
- the D12 resultant algebra;
- the Reidemeister/Whitehead classification;
- Lam–Leung; Coven–Meyerowitz; Łaba–Londner (Invent. Math. 232 (2023)); Sands' list; Diananda–Yap / Green–Ruzsa;
  Klotz–Sander;
- Jeřábek 2016; KNY 2017; Liu–Pass; Ng–Schauenburg;
- Atrubin, Shamir, Brun;
- PLQ 2008, Greene 2013, BGS 2010;
- Seymour, Truemper (cubic), Fink–Moci, McDuff–Schlenk.


## 15. Reproduction

```bash
.venv/bin/python notes/conn_top_scripts/checks.py   # unitary Cayley spectra, Jacobi sum = a_p a_q, C(2x,x) splits, sum-free
.venv/bin/python notes/conn_top_scripts/dinv.py     # lens-space d-invariant collisions (n = 12–16), structure for N = 101·307
.venv/bin/python notes/conn_top_scripts/wendt.py    # Wendt/periodic-dimer resultant asymmetry vs multiplicative order
```
(Run from the project root; scratch copies lived in the git-ignored scratchpad/conn_top/.)
