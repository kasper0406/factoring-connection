# nativeinv: native computation cost of factoring-equivalent invariants (phase 29)

Agent `nativeinv`, 2026-10-27. **Version v2.1 (final: codex blind round merged; Opus proof critic ×2, Opus literature critic and codex adversarial round (MAJOR REVISION) applied; §5 lists the changes).**
- Pre-registration: `notes/nativeinv_prereg.md` (commit `9b8c88e`, before any code ran).
- Code: OPERATIONAL `src/nativeinv/trichotomy.py` (sees N, an invariant value and the public description of the local
  data; exact arithmetic, Sturm root isolation, no factorisation calls; an AST test checks this). HARNESS
  `src/nativeinv/{harness,run_experiments}.py`. Tests: `tests/test_nativeinv.py` (9 pass, < 1 s).
- Raw data: `experiments/raw/nativeinv_e2.jsonl`. Command: `.venv/bin/python -m src.nativeinv.run_experiments`.
- GPT-6.1-sol: blind round `codex/nativeinv_1.{md,out}`; adversarial round `codex/nativeinv_2.{md,out}`.
- Labels: **KNOWN** (with reference), **PROVED** (proof here; routine unless stated; "no prior formulation found by us"
  is the strongest novelty claim), **HEURISTIC**, **experimentally supported**, **OPEN**, **UNVERIFIED**.

Notation: N = pq, p < q odd primes, n = log₂N, s = p + q. "FACT" = randomized-polynomial-time equivalent to factoring
N. "Native algorithm" for an invariant X(N), operationally: an algorithm whose specification never calls a
factoring or divisor-search step; we record per algorithm whether p, q nevertheless appear in its registers (S below
holds ⌊N/p⌋ = q; class-group algorithms output factors as a by-product). This is a description of implementations,
not a complexity class (codex: "never factors N" has no stable formal definition).

---

## 0. Verdict (plain language)

**No factoring-equivalent invariant in the repository has a known native algorithm that beats factoring, and the
"new views" collapse into very few distinct questions.** The inventory has 55 rows (many rows merge several notes'
items). Sorted by what a native algorithm would actually have to compute:

1. **21 of the 55 rows are aliases of one target: p + q.** For every invariant of the form f(N) = F(p)F(q) with polynomial (or
   character-twisted polynomial) local data F — φ, σ, ψ, Pillai, the K-group orders #K_{2j−1}(Z/N), the class counts of
   GL₂, SL₂, the Heisenberg group, the lens-space count, the necklace count mod 2^{n+1}, areas, genera, point counts —
   Prop TRI (proved) shows the value factors N outright unless F is a single monomial. So K-theory, representation
   theory, topology and Burnside counting here are p + q in different clothes. Their best native algorithm is a
   summatory difference: Õ(N^{1/2}) by the hyperbola method (codex's algorithm S), Õ(N^{1/3}) by codex's weighted
   extension of Sladkey's lattice-hull method (SKETCH, §1.3) — the level of Hiary/Lehman, far above factoring.
2. **Four are divisor sums that need a p-th power** (cat-map primitive orbits, necklaces with fixed colours,
   cyclic-cover homology). Only N-th powers are available natively, and those leak only group orders (leakmap).
3. **Fifteen are search problems whose answer is a factor** (idempotents, square roots, tiles, zero divisors).
4. **Two are quadratic-residue bits, not factoring** (lens-space homotopy type, ε): Prop TRI's fiber tier.
5. **The genuine remainder is small:** Frobenius traces at composite index (a_N; for `nativemod`), class numbers
   (for `nativeana`; the only cluster with a native subexponential algorithm, L[1/2], which is itself an index
   calculus), and root-type exponential sums (Salié).

**Why every known native algorithm (except class-group index calculus, L[1/2]) costs N^{Ω(1)}.** This is a
description of the implementations we know, not a lower bound; none is proved. Every definition is a sum over a set
of size N^{Θ(1)}; the compressions we found (CRT; N-exponent gates; ring algebra mod N; reciprocity along a/N;
smoothness relations; weighted lattice-point summation; quotient recursions; compressed linear algebra — the list is
not exhaustive) either need the factors, are covered by a project barrier, are index calculus, or stop at N^{1/3}.
In the point-query model "evaluate f from oracle local data and N as an integer", Prop LOC (proved) shows that
success exposes a query at p or q. A superpolynomial lower bound against nonuniform constant-free SLP splitters
would imply constant-free permanent lower bounds (Bürgisser; codex_r2), so lower bounds of that kind are not in
reach. Prop AUT (known in substance: Schlage-Puchta 2003, Coons 2010): no digit automaton
computes φ, σ, Pillai, J₂ or k(SL₂) on all integers.

**Shortlist (§4), none a route:** S1 exact summatory functions with Möbius weights (Φ, the twisted totient giving
the ring class number h(−7N²)); the unweighted/polynomial-weight case reaches N^{1/3}; S2 digit automata correct on
semiprimes only, for a fiber function such as φ(N) mod 3; S3 compressed integer linear algebra (Smith normal form of
small explicit matrices) with torsion φ(N), where codex proved a fixed-size obstruction. Runner-up: Euclid-computable
homotopy invariants of lens spaces (QR through class E), after a counterexample forced a corrected statement.
Deprioritised by the brief's NP-hardness rule: tilings, matroid realisability, presentation-level group homology,
unitary-Cayley chromatic number, witness counting.

---

## 1. The inventory

### 1.1 Column conventions

- **Native cost (definition).** The cost of the most direct evaluation that does not use p, q.
- **Best native.** An audited upper bound for a method that does not call factoring (not a literature-optimality
  claim). "none" means nothing better than the definition is known to us. Every cluster-I row also has S (Õ(N^{1/2}))
  and its hull refinement (Õ(N^{1/3}), sketch) through s. The fastest method overall is always "factor N (NFS,
  L[1/3, 1.923]) and use the local formula".
- **Gen.** Complexity of the natural generalisation.
  - **≤F** means: the same definition for every N lies in FP^FACT (with a general factoring oracle, so also
    ord_p(a), which needs p − 1 factored). Then NP-hardness under deterministic Turing reductions would put
    NP ⊆ P^FACT ⊆ NP ∩ coNP, i.e. NP = coNP (sihard Remark FO); under randomized reductions, NP ⊆ coNP/poly and
    PH collapses. So ≤F rows are *not* NP-hard unless PH collapses. This needs an actual upper reduction; "summatory
    function of a multiplicative function" is not ≤F (factoring x does not give Φ(x)).
  - Otherwise the class of the natural wider problem (input variety, input group, input tile) is named.
- **Barrier** (what covers the native algorithms we know):
  - **PT**: Theorem PT (factorbits §3.2): pure mod-N ring algebra plus gcd learns nothing without splitting.
  - **PRIM**: N-dependent exponent gates (leakmap NE1, prim): x^{E(N)} leaks only family M. The rows marked PRIM
    need a p-exponent (x^p, ψ^p, A^p, Frob_p) where only N-exponents are available.
  - **SE**: size–entropy law (suffmap Theorem SS): approximate versions are direct targets, i.e. Coppersmith hints.
  - **TRI**: Prop TRI (§3.1): the invariant is s (or a fiber bit) in disguise, so a native algorithm *is* a factoring
    algorithm.
  - **FROB**: nonsmooth Prop FROB: producing the componentwise Frobenius operator on an étale Z/N-algebra factors N.
  - **P1**: sampling floor (conn_ana): the summand deviates from its generic value with probability ≍ 1/p.
  - **PW**: precision wall (nativeana's topic): an individual value is a jump of a summatory or analytic quantity of
    size N^{a}, so it needs absolute error < 1/2, i.e. relative precision N^{−a} (a = 2 for H₁(N)).
  - **EXP**: the exact integer, or the specified explicit model (triangulation, Brandt matrix, Hecke space), has size
    N^{Θ(1)}; this bounds that representation, not every representation.
  - **none**: no project barrier covers the known native algorithms.

### 1.2 The table (55 rows; deduplicated across conn_alg/ana/top/cmp, geometry, newfamily, ecount, tracedata,
factorbits, analytic, nonsmooth, selfred, genusleak, navigate, smoothbuild)

**Cluster I. Quasi-polynomial local data: Prop TRI's FACT tier (one native question: s).**

| # | Invariant | Where | Reduction | Native cost (definition) | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| I1 | φ(N) = #(Z/N)^* = #K₁(Z/N); Res_{s=1}L(s,χ₀) = φ/N; [Q(ζ_N):Q]; degree / top eigenvalue of the unitary Cayley graph; rank of the Artin motive of Q(ζ_N); Knauf spin-chain degeneracy | factorbits §2; conn_alg K6, H7; conn_ana R3, R6, M5; conn_top C5 | s = N + 1 − φ (deterministic; general N: Miller 1976, folklore) | N (gcd per residue) | **S** (Õ(N^{1/2}), §1.4); Φ(N) − Φ(N−1) with the summatory totient Φ in Õ(N^{2/3}) (Pawlewicz–Pătraşcu 2009); Smith: det G_N/det G_{N−1} = φ(N) for G_m = (gcd(i,j)) | ≤F | PT, TRI, P1 |
| I2 | σ(N) = #index-N sublattices of Z²; r₄(N) = 8σ(N); σ₃(N) = r₈(N)/16 = [q^N]E₄/240; ψ(N) = #P¹(Z/N) = [SL₂(Z):Γ₀(N)] = 3·area(X₀(N))/π = deg Φ_N = #T_N-orbit; dim S_{k+12} − dim S_k(Γ₀(N)) (= ψ for even k ≥ 4) | geometry #10, #12, #19a; newfamily 28; conn_alg L3, E1; conn_ana S1, S7, E1; conn_top D7, D11 | s = σ − N − 1 (general N: Bach–Miller–Shallit 1986, σ_k and r₄, r₈) | N (divisor sum; modular symbols; any triangulation has ≥ ψ/6 − 7 edges) | **S**: σ_k(N) = H_k(N) − H_k(N−1), hyperbola method, Õ(N^{1/2}) (codex, §1.4); r₄, r₈ by theta-series multiplication Õ(N) | ≤F | TRI, PW, EXP |
| I3 | genus g(X₀(N)); dim S_k(Γ₀(N)); TMF₀(N) ranks; Farey-symbol size | geometry #13; conn_top H3 | four candidates for ψ | Θ(g) cells | S (via ψ) | ≤F | TRI, EXP |
| I4 | J₂(N) = (p²−1)(q²−1) = k(GL₂(Z/N)) = #K₃(Z/N) = \|SL₂(Z/N)\|/N; irreducible-degree set of GL₂(Z/N) (contains ψ, φ, N); g(X(N)), g(X₁(N)) | conn_alg R1, R6, K6; conn_top H5; geometry #19a | p² + q² = N² + 1 − J₂ | k(GL₂): Σ_{t∈Z/N, d∈(Z/N)^*} 2^{ω(gcd(t²−4d,N))} (Nφ(N) terms; E1) or N⁸ commuting pairs | S (via s) | ≤F; k(G) for black-box matrix groups: no hardness result known to us (UNVERIFIED) | TRI, P1 |
| I5 | #K_{2j−1}(Z/N) = (p^j−1)(q^j−1), j ≥ 1; #K₂(Z[1/N]) = 2φ(N) | conn_alg K3, K6; conn_top H5 | p^j + q^j, then s | the definition gives no finite algorithm by itself; group homology of GL_m(Z/N) (\|G\| ≈ N^{m²}) | S (via s); no K-theoretic shortcut (§2.1) | ≤F | TRI, PRIM (ψ^p vs ψ^N), FROB |
| I6 | Pillai P(N) = Σ_g gcd(g,N) = #{xy ≡ 0 mod N} = (2p−1)(2q−1); star-polygon census; cycles of the shear (x, y) ↦ (x + y, y) on (Z/N)² (codex row 22); necklace count at k = 1 + R mod R² (Prop NECK, one query) | geometry #34, #36 | s | N | S | ≤F | TRI, P1 |
| I7 | Two-colour necklace count Neck(N, 2) mod 2^{n+1}: N·Neck ≡ 2φ(N) (mod 2^{p}) | codex nativeinv_1 (blind); proof critic | φ exactly when p > n + 1 (deterministic); with only 2^{n/4+O(log n)}: Coppersmith/BDF needs k ≥ n/4 + v₂(q−p) | N terms; the integer has ≈ N bits | S (via φ) | ≤F | TRI, SE |
| I8 | k(SL₂(Z/N)) = (p+4)(q+4) = Σ_{t∈Z/N} 3^{ω(gcd(t²−4,N))} (N odd squarefree) | KNOWN field count + CRT; gcd-sum form here; E1 checks N ≤ 35 | s | N | S | ≤F | TRI, P1 |
| I9 | Aff(Z/N): irreducible-degree multiset {1^{φ}, (p−1)^{q−1}, (q−1)^{p−1}, φ^{1}} (k(Aff) = N carries nothing) | conn_cmp 2.5; E1 | max degree = φ | Clifford theory: orbits of (Z/N)^* on the dual group = divisors of N | S (via s) | ≤F | TRI (CRT idempotent) |
| I10 | k(Heis(Z/N)) = (p²+p−1)(q²+q−1); k(U_n(Z/N)) wherever k(U_n(F_r)) is a non-monomial polynomial in r (Higman's conjecture; verified n ≤ 13 by Vera-López–Arregi 2003, n ≤ 16 by Pak–Soffer, arXiv:1507.00411) | KNOWN field count + CRT; E1 check | s | N³ enumeration | S (via s) | ≤F for fixed n; pattern groups: class counts can be non-polynomial (Halasi–Pálfy 2011); no complexity result known to us | TRI |
| I11 | #subgroups of (Z/N)² = (p+3)(q+3); Σ over subgroups of C_N of their orders = σ(N) | conn_top C7; E1 (N = 15) | s | HNF enumeration over divisors | S (via s) | ≤F | TRI |
| I12 | #V(Z/N) for V with quasi-polynomial counts: G_m, A¹∖{0,1}, SL₂, Grassmannians, split tori (Prop MF); the conic x²+y² = 1 (Prop MF′); Markoff X*(Z/N) = ∏(r² ± 3r); #P²(Z/N) = (p²+p+1)(q²+q+1) (definable sets in general are only μr^d + O(r^{d−1/2}), Chatzidakis–van den Dries–Macintyre: not in this row) | geometry §4.1; conn_alg M2; conn_ana D3; nonsmooth §5.3; selfred SR2 | Prop TRI | N^{dim V} | S (via s) | fixed quasi-polynomial-count V: ≤F; V in the input: #P-hard (counting solutions of polynomial systems over F_2) | TRI, PT |
| I13 | #X₀(N)(F_{ℓ²}) = (ℓ−1)ψ/12 + bounded QR terms; class number of an Eichler order of level N; Brandt-matrix size; tr B(ℓ^{2k}) | geometry Prop X0-3; Kirschmer–Voight 2010 Prop 5.8; conn_top G5; conn_ana R7 | ψ | Brandt matrices of size ≍ N/12 | S (via s) | ≤F | TRI, EXP |
| I14 | #X₁(N)(F_ℓ) = φ(N)/2 (only cusps, once every r \| N exceeds ℓ+1) | geometry Prop X1 | φ | point count of a curve of genus ≍ N²/24 | S (via s) | ≤F | TRI |
| I15 | Bianchi covolume of Γ₀(N·O_K): ∏_{r\|N}(local factor) | geometry #21 | s | Swan/Riley domains ∝ index | S (via s) | ≤F | TRI, EXP |
| I16 | Bost–Connes Galois sum (p^{1−β}−1)(q^{1−β}−1), real β > 1 (zero at β = 1; precision/recovery theorem not supplied: UNVERIFIED); main terms of moments of L(½,χ), χ mod N; φ*(N) = (p−2)(q−2) primitive characters | conn_cmp 2.1; conn_ana R1, R2 | s (archimedean decode or exact) | N L-functions | S (via s) | ≤F | TRI, PW |
| I17 | Number of lens spaces with π₁ = Z/N up to homeomorphism = (φ(N) + 4 + 4[−1 ∈ (Z/N)^{*2}])/4 (Burnside for a ↦ ±a^{±1}) | this note, §2.3 (PROVED, elementary) | φ (two candidates) | φ(N) residues | S (via s) | ≤F | TRI |
| I18 | Faltings height of J₀(N), ω² (Abbes–Ullmo); heat-trace / Weyl / Selberg-zeta leading terms on X₀(N) | conn_alg H5; conn_ana S6, S7 | area terms: s; Faltings height / ω²: exact reduction not supplied (UNVERIFIED) | genus-≍N/12 models | S (via s) | ≤F (real precision) | TRI, EXP, PW |
| I19 | Ring class number h(−7N²) = #{x² + 7y² = 1 over Z/N} (needs h(−7) = 1, unit index 1, 7 ∤ N) = (p − χ(p))(q − χ(q)), χ = (−7/·); generally #{x² − Dy² = 1}(Z/N), fixed D | codex nativeinv_1 rows 30–31 (blind; Cox Thm 7.24) | one scalar gives s or \|p − q\| (Prop TRI(a); the sign pattern is a fiber point) | reduced forms of disc −7N²: ≈ N^{1+o(1)} | twisted summatory recursion F_χ(N) − F_χ(N−1), Õ(N^{2/3}) (codex) | ≤F | TRI, PW |
| I20 | Smith's gcd determinant ratio det G_N / det G_{N−1} = φ(N), G_m = (gcd(i,j))_{i,j≤m}; G_N = E·diag(φ(1..N))·E^T | codex row 36; Smith 1875 | φ | N×N matrices | as I1 | ≤F (explicit determinants are in P; succinct dimension N is the obstacle) | TRI |
| I21 | Counts that are φ or J_k in topological dress: faithful 2-dim irreps of D_{2N} (φ/2); faithful rank-1 local systems on L(N,1) (φ); surjective labelled Fox N-colourings of T(2,N) (Nφ); generators of the torsion of H₁ of the mapping torus of [[1,N],[0,1]] (φ); epimorphisms H₁(S_N) → C_N for the Euler-number-N circle bundle over T² (J₃) | codex rows 16–20 | φ / J₃ | N^{O(1)} | S | ≤F | TRI (homology is C_N, easy; the primitive/surjective filter is the hard part) |

**Cluster II. Möbius / divisor sums over exponent data (the p-exponent pattern).**

| # | Invariant | Where | Reduction | Native cost | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| II1 | Primitive period-N points of a toral automorphism, Prim_t(N) = Σ_{d\|N}μ(N/d)(V_d(t)−2); trace t = 2 + R (one query, Prop CAT); fixed t = 3 (figure-eight monodromy, multi-query small moduli) | geometry #37, #26 | J₂ (t = 2+R); s mod r per small r (fixed t) | needs V_p(t), V_q(t); integer has ≍ N bits | Fix(A^N) = V_N(t) − 2 is poly-time (all periods dividing N), but carries nothing | ≤F | PRIM, EXP |
| II2 | (fixed-value questions, fiber tier; FACT only with selectable moduli/colours) Lyndon/necklace counts with fixed colours mod small M; cyclotomic cosets Σ_{d\|N}φ(d)/ord_d(2) (cyclic/BCH/duadic codes); simple F_ℓ[C_N]-modules | newfamily 43, 45; conn_alg R5; conn_ana I6 | s mod (period) per modulus; orders | N | none | ≤F | PRIM, M |
| II3 | Primitive parts of cyclic-cover homology Σ μ(N/d)\|H₁(X_d)\|; Floer HF(φ^N) primitive part | geometry #26; conn_top S2 | Möbius over divisors | needs \|H₁(X_p)\| | each single \|H₁(X_N)\| mod M is AGG | ≤F | PRIM, EXP |
| II4 | λ(N); ord_N(a) with selectable a (FACT); period of one fixed cat map or of x² + c mod N (no reduction proved); #periodic points of x ↦ x² on (Z/N)^* = oddpart(φ(N)) (codex row 25) | factorbits; newfamily 36, 37 | Miller 1976 | BSGS O(N^{1/2}); rho | none native below factoring | ≤F | M, PRIM |

**Cluster III. Non-abelian Frobenius local data (deep dive: `nativemod`).**

| # | Invariant | Where | Reduction | Native cost | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| III1 | #E(Z/N) for chosen curves / twists | factorbits; selfred; ecount | twist counts sum to 4(p+1)(q+1) (Dieulefait–Urroz eq. (1); Kunihiro–Koyama) | Schoof over Z/N needs Frobenius; group order of E(Z/N) | BSGS on E(Z/N) O(N^{1/2}) (a factoring algorithm); ECM factors in L_p[1/2] | ≤F | FROB, PT |
| III2 | #E(Z/N) for one fixed curve | tracedata Theorem EXT | ≡_R FACT on every good N | as III1 | none | ≤F | FROB |
| III3 | (a_N, a_{N²}) (Bach–Charles; fails iff a_p = a_q = 0); (a_N, a_{N³}) (works iff a_N ≠ 0); τ(N) with τ(N²); a_N alone OPEN | tracedata §2.2; conn_top C4 | Bach–Charles 2007 | Edixhoven–Couveignes needs N factored | none | ≤F | FROB (operator only) |
| III4 | Class-number moments from the Eichler–Selberg formula: Σ_t P_k(t,N)H(4N−t²) = −Σ_{dd′=N}min(d,d′)^{k−1} when S_k(1) = 0 (the trace itself is then 0); Kronecker–Hurwitz Σ_t H(4N−t²) = 2(N+q); weight 4: −2(1+p³) | conn_alg L2 | hyperbolic (divisor) term | ≈ 4√N class numbers | none | ≤F through the divisor formula | PW (`nativeana`, `nativemod`) |
| III5 | Roots of Φ_ℓ(j, X) mod N (genus-0 X₀(ℓ)); ℓ-isogeny neighbours over Z/N | conn_cmp 4.1; Altuğ–Chen Thm 3.3 | planted Fricke root (randomized) | — (search) | none | ≤F | PT (search = splitting) |
| III6 | #Hom(π₁M, SL₂(Z/N)) = character-variety point counts over F_p, F_q | conn_top K5 | TRI when the count is a non-monomial quasi-polynomial; otherwise III-type | N^{3g} | none | fixed M: ≤F only for quasi-polynomial counts (poly(log r) point counting for general fixed varieties is not known) | TRI / FROB |

**Cluster IV. Global, non-multiplicative invariants (class groups; deep dive: `nativeana`).**

| # | Invariant | Where | Reduction | Native cost | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| IV1 | h(−4N), r₃(N) (FACT-hard only on N ≢ 7 mod 8) (= 12h(−4N) for N ≡ 1, 2 mod 4; 24h(−N) for N ≡ 3 mod 8; 0 for N ≡ 7 mod 8), the class group of disc −4N (ambiguous forms), Fricke fixed points on X₀(N), reduced lattices of disc −4N | geometry #5, #11, #19; newfamily 27 | FACT-hard direction: h ⇒ elements of order 2 ⇒ ambiguous forms (randomized) | O(√N) reduced forms | **native and subexponential**: Hafner–McCurley L[1/2] (GRH); Shanks O(N^{1/4}), Lenstra–Schoof N^{1/5} (ERH) | not known ≤F (hard at prime discriminants too) | none (index calculus escapes TRI and PT) |
| IV2 | Partial sums of Jacobi symbols: for N ≡ 3 (mod 4) squarefree, Σ_{a<N/2}(a/N) = (2 − χ(2))h(−N) (class-number formula; zero for N ≡ 1 mod 4) | this note §3.4 | ⇒ h ⇒ FACT-hard | N terms, each poly-time by reciprocity | via the class group, L[1/2] | not known ≤F | PW |
| IV3 | Genus field / Hilbert class field of Q(√−N) (explicit field structure; the Rédei 4-rank / (p/q) are fiber-type bits, not FACT) | conn_alg G2; geometry #30 | genus field contains √p* | class group | L[1/2] | not known ≤F | none |

**Cluster V. Root-type sums and reciprocity local parts.**

| # | Invariant | Where | Reduction | Native cost | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| V1 | Salié S(a,b;N); Kloosterman K(a,b;N²) | newfamily 21, 22 (Prop SAL) | FACT at precision N^{−4} (PROVED); at δ = O(1/log N) (HEURISTIC) | N terms | none; at a prime modulus Salié has a closed form through √(ab) mod p (Tonelli–Shanks) | ≤F | none of PT/PRIM (archimedean shadow of √ mod N) |
| V2 | Cubic Weyl / Gauss / Jacobi sums mod N | newfamily 24 | FACT-hard on p ≡ q ≡ 1 (mod 3) with selectable queries | N | none | ≤F | — |
| V3 | Dedekind local part s(a,p) + s(a,q); unit part U(a,N) (Prop DSF needs gcd(2a(a²+1), N) = 1; e.g. N = 65, a = 8 gives 0) | newfamily 7 (Prop DSF) | FACT (one value, deterministic) | reciprocity at p needs p | s(a,N) itself is poly-time (class E) but not FACT | ≤F | RA / E |
| V4 | binom(2N, N) mod N | newfamily 54 | FACT-hard on q < 3p/2 | integer of 2N bits; mod-N evaluation meets factorials mod N | none native (Pollard–Strassen/BGS-type N^{1/4} would already factor) | ≤F | Pochhammer/τ barrier |

**Cluster VI. Quadratic (fiber) tier, listed because a deep-dive item lives here; NOT known FACT.**

| # | Invariant | Where | Status | Native cost | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| VI1 | Lens spaces: oriented homotopy L(N,a) ≃⁺ L(N,b) ⇔ a/b ∈ (Z/N)^{*2} (QRP; 4 classes); unoriented ⇔ ±squares (2 + 2[−1 ∈ □] classes); the map on H₁ of an equivalence (multiplication by a square root n; the degree is ±1) | conn_top K1; §2.3 | QRP; a suitably quantified square-root search oracle factors (Rabin) | linking form from a triangulation: poly-time; its isometry class needs Legendre symbols | a normalised Gauss sum of the linking form gives the Jacobi symbol (AGG); the full Witt class retains local data | ≤F, in NP ∩ coNP | PT (generic-ring QR ≡ FACT, Jager–Schwenk), RA, TRI-FIBER |
| VI2 | ε = (−1/p), r₂(N), ν₂, ν₃, #X₀(N)(F_ℓ), φ(N) mod ℓ, [ℓ \| φ(N)] (Φ-hiding) | geometry X0-2; factorbits; nonsmooth RANK | fiber bits (O(1) bits each) | — | none | ≤F | PT-m fiber floor |

**Cluster VII. Witness, search and structure invariants (RESTATEMENT; a native algorithm is a factoring algorithm).**

| # | Invariant | Where | Reduction | Native cost | Best native | Gen. | Barrier |
|---|---|---|---|---|---|---|---|
| VII1 | Finding a nontrivial idempotent of Z/N (deciding ∃e is trivially true on the promise); splitting the unit of Perf(Z/N); nontrivial automorphism of a ring (Z/N)[x]/(f), f a random monic cubic, given in basis form (Kayal–Saxena Thm 8.1: finding ≡_ZPP FACT; Z/N itself has none) | conn_alg M5; conn_cmp 1.4; conn_top C7 | witness = factor | — | — | ring isomorphism: GI-hard, in NP ∩ coAM; deciding existence of an automorphism is in P; finding one ≡_ZPP FACT (Kayal–Saxena, Comput. Complexity 15 (2006)) | PT |
| VII2 | Square roots mod N (Rabin); Weil-representation isotypic projectors of SL₂(Z/N); H₁-maps of lens-space homotopy equivalences, as a search oracle on random instances | conn_alg R2; conn_top K1 | Rabin 1979 | — | — | ≤F | PT |
| VII3 | Two inequivalent representations N = x² + y² (or by a fixed norm form, under a promise that two exist); Thue points F(x,y) = N | conn_alg A3 | Euler; McKee 1996 | Cornacchia needs √−D mod N | — | ≤F for binary quadratic forms; general Thue solving: not justified | Q / SQRT |
| VII4 | Zero divisors in quaternion algebras over Q; Wedderburn decomposition of Q-algebras; points on ax² + by² = Nz² | conn_alg R4, G1 | Rónyai 1987/88 (randomized) | — | poly-time given factoring (Ivanyos–Rónyai–Schicho 2012) | ≤F | — |
| VII5 | Frobenius operator on étale Z/N-algebras satisfying FROB's separation hypotheses (on B = Z/N, Φ = id) (Schoof's (Z/N)[x]/ψ_ℓ with (N/ℓ) = −1); Frobenius lifts / δ-structures at N | nonsmooth Prop FROB; conn_alg P2 | FROB(a–d) | — | — | ≤F | FROB |
| VII6 | Matroid M_N with characteristic set {p, q} (every finite set of primes is a characteristic set: Kahn 1982) | conn_top X1 | restatement | — | — | realisability over R is ∃R-complete (Mnëv 1988; Shor 1991; Schaefer 2009); over C / finite characteristic it is Nullstellensatz-type (NP-hard, in PSPACE; AM under GRH, Koiran 1996) | DEPRIORITISE |
| VII7 | Producing a tile of a nontrivial tiling of Z_N (its size is a divisor); counting tilings; Coven–Meyerowitz cyclotomic data | conn_top (tilings) | restatement (producing); counting and deciding are different tasks | counts ≥ q^{p−1}, i.e. ≥ Θ(p log q) bits | — | tiling a finite region with given tiles: NP-complete (Garey–Johnson; Moore–Robson 2001); whether one tile tiles a given periodic subset of Z^d is undecidable (Greenfeld–Tao, JEMS 2025) | DEPRIORITISE |
| VII8 | Addition / successor in the prime-exponent representation | factorbits Props E1, E2 | ≡ FACT | — | — | ≤F | — |
| VII9 | Partial factor information: p mod M with log M ≥ n/4; φ mod 2^k; top n/4 bits of p; "N has a factor ≤ x"; smallest prime factor | factorbits table | Coppersmith; Boneh–Durfee–Frankel | — | — | ≤F | SE |
| VII10 | Hidden-order groups: a multiple of λ(N); a generic speed-up of repeated squaring (Rotem–Segev 2020) | conn_cmp 4.5 | Miller | — | — | ≤F | PT (generic) |
| VII11 | Witness-counting polynomials Σ_{x,y∈{0,1}^{n}}[xy = N]·(prefix of x); sum-check provers for them | conn_cmp 1.7, 9.4 | prefix queries | 2^n | — | #CircuitSAT: #P-complete | DEPRIORITISE |
| VII12 | Short-quotient navigation, CLASS_∞, hidden-root samplers | navigate; smoothbuild | FACT ∧ DL | — | — | ≤F+DL | SE |
| VII13 | Integer relations among Γ(a/N); Ramanujan expansion of σ(N) (terms with modulus < p are blind) | conn_alg H1; nonsmooth §4.1 | relation length ≥ p | — | — | — | PW |
| VII14 | Order of SL₃(Z/N), given by an O(log N)-size presentation computable from N alone: N³(p²−1)(p³−1)(q²−1)(q³−1) | this note §2.1 (KNOWN ingredients) | s | Todd–Coxeter on a vector stabiliser enumerates ≈ J₃(N) ≍ N³ cosets | S (via s) | the promised family is ≤F; general finitely presented groups: undecidable (Adian–Rabin), which does not obstruct this family | presentation-level: no barrier |
| VII15 | Chromatic (= clique) number of the unitary Cayley graph = least prime factor (Klotz–Sander 2007) | codex row 24 | restatement | — | — | graph colouring: NP-hard | DEPRIORITISE |

Total: 21 + 4 + 6 + 3 + 4 + 2 + 15 = 55 rows (several rows merge 3–6 invariants from different notes). **Rows are
not independent targets** (codex): I20 repeats I1, I4's K₃ repeats I5, I19 overlaps I12, I21's J₃ is I5's K₅. Row
percentages are percentages of interpretations, not of computational problems.

### 1.3 The one native algorithm that every cluster-I row has (codex blind round)

**Algorithm S** (codex nativeinv_1, PROVED there and re-checked by us at the formula level). With
H_k(x) = Σ_{m≤x}σ_k(m), P_k(t) = Σ_{j≤t}j^k, h = ⌊√x⌋, the hyperbola identity
H_k(x) = Σ_{d≤h} d^k⌊x/d⌋ + Σ_{j≤h} P_k(⌊x/j⌋) − hP_k(h) gives σ_k(N) = H_k(N) − H_k(N−1) in Õ(N^{1/2}) without any
divisibility test. Hence s = σ(N) − N − 1, and by Prop TRI every cluster-I invariant, natively in Õ(N^{1/2}).
- This is native only in a syntactic sense: the registers contain ⌊N/d⌋ for every d ≤ √N, including ⌊N/p⌋ = q. Its
  cost is trial division's. It is the honest "definition-level" baseline for cluster I, and it shows that the
  question for cluster I is the *summatory* question: can H_k(x) (or Φ(x), or a twisted F_χ(x)) be evaluated
  exactly much faster than x^{1/2}? Known: Φ(x) in Õ(x^{2/3}) (Pawlewicz–Pătraşcu 2009); the divisor summatory
  function in Õ(x^{1/3}) by lattice-hull methods (Sladkey, arXiv:1206.3369, unrefereed; D. Li, arXiv:2609.29949,
  x^{9/28}, randomized, preprint; conn_ana). No published Õ(x^{1/3}) algorithm for Σ_{n≤x}σ(n) was found, but **codex (adversarial round) gives
  the extension**: run Sladkey's unimodular decomposition carrying the count and the first-coordinate moment; under
  its affine coordinate changes the weight d stays affine, and on a triangle T_t = {u, v ≥ 0, u + v ≤ t} the moments
  are closed forms (Σu = t(t+1)(t+2)/6), so each piece still costs O(1) operations. Hence H₁(x) = Σ_{ab≤x} a, and
  σ(N) = H₁(N) − H₁(N−1), in Õ(N^{1/3}) natively. **SKETCH (codex), judged credible by the round-2 proof critic; not
  implemented or checked by us.** Both linear moments (Σa and Σb) and the weighted versions of every piece
  (rectangles, shared boundary lines, the brute-force strip) are needed; all are closed forms. Weighted hull walks
  for Σ f(i)⌊n/i⌋ are probably folklore in competitive programming; no novelty is claimed. That is the Hiary/Lehman level, far slower than factoring. Möbius-weighted sums (Φ, F_χ) do not
  get this extension automatically: shortlist S1.
- The coefficient sequence of E₄ (σ₃) is not D-finite (folklore; codex's argument) (E₄ blows up at every rational cusp,
  so its q-series has infinitely many singular roots of unity). So no fixed-order polynomial-coefficient
  recurrence gives σ₃(N); this closes that representation (and with it the generic Chudnovsky/BGS Õ(N^{1/2})
  evaluation of P-recursive sequences), not digit algorithms or nonlinear recursions.

### 1.4 What the table says, in one paragraph

Of the 55 rows, 21 (cluster I) are Prop TRI's FACT tier: the invariant equals F(p)F(q) for a quasi-polynomial F, so it
*is* s in disguise (as a target; a native algorithm for one interpretation need not respect another's restriction) and every rewriting of it (homology, representation theory, point counting, Burnside) is a
rewriting of s. Four more (cluster II) are divisor sums whose terms need a p-th power. Fifteen (cluster VII) are
search problems whose witness is a factor. Two (cluster VI) are QR fiber bits, not FACT. What remains is small:
non-abelian Frobenius data (cluster III, `nativemod`), class groups (cluster IV, `nativeana`) and root-type sums
(cluster V). Exactly one cluster has known native algorithms that are subexponential and do not need the factors:
class groups (IV1–IV3), L[1/2], historically one source of subexponential factoring.

---

## 2. Deep dives (combinatorial and algebraic items)

Each item asks: is there a poly(log N)-size representation with a native algorithm that avoids enumerating size-N
objects? Tools tried: transfer matrices, automata over the digits of N, closed-form generating functions,
Burnside/character identities without divisors, homology of small complexes. For each item: **the obstruction
encountered by the constructions examined** (v1 said "forcing step"; codex rightly objects that no necessity is
proved outside the restricted theorems named), or OPEN.

### 2.1 K-groups of Z/N

- **Facts.** Z/N ≅ F_p × F_q, so K_*(Z/N) = K_*(F_p) ⊕ K_*(F_q) and #K_{2j−1}(Z/N) = (p^j − 1)(q^j − 1) (Quillen 1972).
  j = 1 is φ(N). By Prop TRI each j is the FACT tier, and all j are deterministic-polynomially inter-reducible
  (each gives s).
- **Compression attempts.**
  1. *Adams operations.* Quillen: K(F_r) is the homotopy fibre of ψ^r − 1 on BU. The only operation available from
     N alone is ψ^N, and hofib(ψ^N − 1) has π_{2j−1} = Z/(N^j − 1): the "generic" answer, an aggregate. The true
     answer needs ψ^p and ψ^q separately. This is the same split as the generic class number (Green 1955; Steinberg)
     (k(GL₂(F_x)) = x² − 1 evaluated at x = N is wrong by p² + q² − 2) and as V_N versus V_p + V_q in Prop CAT.
     **Obstruction: a p-exponent (PRIM).**
  2. *Small complexes.* SL₃(Z/N) has a presentation of size O(log N) computable from N (VII14): by the congruence
     subgroup property for SL₃(Z) (Mennicke 1965; Bass–Milnor–Serre 1967: SL₃(Z, NZ) = E(3, Z, NZ), the relative
     elementary group, i.e. the normal closure in E₃(Z) = SL₃(Z) of the e_ij(N); all e_ij(N) are conjugate to
     x₁₂^{±N} by signed permutation matrices), Γ(N) is the normal closure of x₁₂(N), so
     SL₃(Z/N) = SL₃(Z)/⟨⟨x₁₂^N⟩⟩. SL₃(Z) has the Steinberg presentation plus (x₁₂x₂₁^{−1}x₁₂)⁴ = 1 (Milnor,
     *Introduction to Algebraic K-theory*, Cor. 10.3, number UNVERIFIED), and x₁₂^N costs O(log N) auxiliary generators y_{k+1} = y_k². So a
     poly(log N)-size 2-complex X_N with π₁ = SL₃(Z/N) is computable from N alone (KNOWN ingredients; the assembly
     is routine). From X_N, H₁ (Smith normal form) is free to compute, and it is trivial. The order of the group,
     and the group homology H₃(G; Z) (which carries K₃ after a stability comparison), are not. (X_N is
     2-dimensional, so its own H₃ is 0; a finite group has no finite-dimensional K(G,1), so higher resolution data
     are required.)
     - **Obstruction encountered (HEURISTIC).** Extending X_N needs the module π₂(X_N) of identities among
       relations. The identities we know for SL_n(F_r) are built from norm elements Σ_{c∈C} c of cyclic subgroups of
       orders dividing r, r ± 1, r² + r + 1. Norm sums do compress when the order is known (S_{2r}(g) = S_r(g) +
       g^r S_r(g), codex; e.g. the global order N of an elementary unipotent), so the obstruction is the unknown
       *local* orders in the constructions examined, not a theorem. Todd–Coxeter on the stabiliser of a vector
       materialises ≈ J₃(N) cosets.
     - **Small homology in general (codex, adversarial round; PROVED).** On the semiprime promise, deterministic
       poly-time factoring ⇔ a deterministic poly-time construction of an explicit poly-size integer chain complex
       whose specified torsion order is φ(N) (⇐: Smith normal form; ⇒: boundary map [φ(N)]). So no general
       obstruction can be proved without separating factoring from FP. (The torsion of H_i is the torsion of
       coker ∂_{i+1}, so one matrix suffices.) **Restricted obstruction (codex; round-2 critic fixes):** a fixed-size
       matrix A(N) with entries in Z[N] cannot have |tors coker A(N)| = φ(N) even on all but finitely many semiprimes.
       Let ρ be the rank of A over Q(N). Outside the finitely many roots of one nonzero ρ×ρ minor, rank A(N) = ρ and
       the torsion order is the gcd d_ρ of the ρ×ρ minors. Write those minors as h(N)q_i(N) with h their gcd in Z[N]
       (Gauss's lemma: q_i ∈ Z[N], coprime over Q[N]); Bézout gives Σ b_i q_i = C ≠ 0 with b_i ∈ Z[N], so
       g(N) = gcd_i q_i(N) depends only on N mod C, and the order is |h(N)|·g(N mod C): a quasi-polynomial in N. By
       Prop TRI(c) (with the finite-exception remark there) F(r) = r − 1 is not of that form. Growing matrices, digit-
       or floor-dependent entries and Euclidean recursions are not covered (shortlist S3).
     - The natural generalisation (order or homology of a finitely presented group) is undecidable (Adian–Rabin;
       Gordon 1995 for H₂ = 0, LMS Lecture Note Ser. 204), so as a *presentation* problem this is deprioritised. The black-box
       matrix-group version is ≤F (+DL): Babai–Beals–Seress (STOC 2009) compute orders of matrix groups over F_q of odd
       characteristic in randomized polynomial time given factoring and DL oracles (over Z/N one first splits by CRT).
  3. *Automata / generating functions.* §3.3 (Prop AUT): no weighted automaton over base-k digits computes φ (or any
     cluster-I invariant) for all n.
- **Verdict.** The constructions examined stop at a p-exponent (Adams) or at local orders (identities among
  relations). As a target it is s. Small-complex homology for φ(N) is open beyond the fixed-size obstruction (S3).

### 2.2 Conjugacy-class counts and irreducible degrees of GL₂(Z/N), SL₂(Z/N), Aff(Z/N)

- **Facts (KNOWN field counts + CRT; E1 brute-force checks).** k(SL₂(Z/N)) = (p+4)(q+4) (N ∈ {15, 21, 33, 35}); k(GL₂(Z/N)) = (p²−1)(q²−1)
  (N ∈ {15, 21}); k(Aff(Z/N)) = N; k(Heis(Z/N)) = (p²+p−1)(q²+q−1). All checks held, as pre-registered.
- **Native formula (PROVED, elementary; N odd and squarefree).** Classes of SL₂ over a field split by the characteristic polynomial
  x² − tx + 1: one class when t² ≠ 4, three when t = ±2 (central plus two unipotent). Over Z/N, by CRT,

  k(SL₂(Z/N)) = Σ_{t∈Z/N} 3^{ω(gcd(t²−4, N))} = Σ_{d\|N} 2^{ω(d)}·(N/d)·#{t mod d: t² ≡ 4} = N∏_{r\|N}(1 + 4/r).

  Likewise k(GL₂(Z/N)) = Σ_{t, d∈(Z/N)^*} 2^{ω(gcd(t²−4d, N))} (checked for N = 15, 21; both formulas fail for non-squarefree N, e.g. N = 9: k(SL₂) = 25 vs 21, k(GL₂) = 78 vs 72 — proof critic). So the native object is a
  **gcd-sum** Σ_{x∈(Z/N)^k} w(gcd(P(x), N)) with fixed P, w. Pillai (P = x, w = id), φ (w = [· = 1]), necklaces
  (w = k^·), Burnside counts and point counts all have this shape.
- **Compression attempts.**
  1. *Burnside / character identities.* k(G) = #{commuting pairs}/|G| = Σ_χ 1 replaces one size-N^{Θ(1)} sum by
     another. For **radial** summands w(gcd(x, N)) (Pillai, φ, necklaces), the functions on Z/N invariant under
     (Z/N)^* form a τ(N)-dimensional algebra (the rational Schur ring of Z/N; KNOWN), whose orbit basis and
     Ramanujan-sum basis are both indexed by the divisors of N, and Burnside/Fourier rewritings stay inside it. The
     class-count summands w(gcd(P(t), N)) with nonlinear P are **not radial** (codex: at N = 15, t = 1, 2 are units
     but gcd(t² − 4, 15) = 3, 15); they are pushforwards of counting measure under P, whose distribution would have
     to be computed. A divisor-indexed basis does not force an algorithm to construct that basis (here it has
     dimension 4). Obstruction encountered: the known evaluations use the divisor-indexed coefficients (CRT).
  2. *Clifford theory for Aff(Z/N).* The degrees come from the orbits of (Z/N)^* on the dual of Z/N, and these
     orbits are the divisors. k(Aff(Z/N)) = Σ_{d\|N}φ(d) = N collapses (the Gauss identity), which is exactly why it
     carries nothing; the degree multiset does not collapse (Prop TRI FACT tier at the maximal degree φ(N)).
  3. *Generic class polynomials.* For GL_n and SL_n, k(G(F_r)) is a polynomial (quasi-polynomial for SL_n) in r
     (Green; Deligne–Lusztig; Lehrer). Evaluating it at N gives the AGG value; the truth is P(p)P(q). Prop TRI.
  4. *Sampling.* The summand 3^{ω(gcd(t²−4,N))} differs from 1 with probability ≈ 2/p + 2/q: P1 floor.
- **Generalisation.** k(G) for a permutation group given by generators: neither a polynomial-time count nor a
  hardness result was found (literature critic; listing class representatives is not a count, since k(G) can be
  exponential, Garonzi–Maróti 2015). Pattern groups U_P(F_q): class counts can be non-polynomial in q
  (Halasi–Pálfy, J. Group Theory 14 (2011)); no complexity result is known to us. (v0 attributed a universality
  theorem to Pak–Soffer; the critic checked: they prove none. Removed.)
- **Verdict.** The evaluations examined use CRT idempotents / divisor-indexed coefficients. As a target it is s.

### 2.3 Lens-space invariants

- **Facts.**
  - H₁(L(N,a)) = Z/N and the linking form λ(g, g) = ±a^{±1}/N mod 1 (up to convention) are computable from any triangulation in poly time
    (Smith normal form plus one intersection computation).
  - Oriented homotopy type: L(N,a) ≃⁺ L(N,b) ⇔ ab^{−1} ∈ (Z/N)^{*2}; unoriented: ⇔ ab^{−1} ∈ ±(Z/N)^{*2} (Whitehead;
    Olum). So the complete homotopy invariant is the isometry class of a rank-1 linking form, i.e. the pair of
    Legendre symbols ((a/p), (a/q)) up to the known product: **QRP**, not FACT.
  - Homeomorphism: b ≡ ±a^{±1} (Reidemeister, Brody): four congruences, trivially poly-time. So the *finer*
    classification is easy and the *coarser* one is QR-hard. The reason: homotopy asks ∃n (ab ≡ n²), an existential
    over a square root. This is for **oriented** homotopy. Unoriented (±□): if p ≡ q ≡ 3 (mod 4), −1 is a
    Jacobi-(+1) nonsquare, ±□ = {(a/N) = 1} and the decision is poly-time; if p ≡ 1, q ≡ 3 (mod 4) it reveals (a/p).
  - **Counts (PROVED, Burnside; likely folklore).** Up to (unoriented) homeomorphism there are (φ(N) + 4 + 4[−1 ∈ □])/4 lens spaces with π₁ = Z/N
    (the Klein four-group {a ↦ ±a^{±1}} acting on (Z/N)^*; fixed points: a² ≡ 1 (4 of them), a² ≡ −1 (0 or 4),
    a ≡ −a (none, N odd)). This is FACT (row I17). Up to unoriented homotopy there are |(Z/N)^*/±□| = 2 + 2[−1 ∈ □]
    classes: a fiber bit (row VI1); up to oriented homotopy, 4. So the natural counts land in different tiers of
    Prop TRI.
- **Compression attempts.**
  1. *Small complexes.* The standard CW structure has one cell per dimension; its homology with any fixed
     coefficients is an aggregate (geometry Prop AGG-H). Cup products with Z/r coefficients need r \| N.
  2. *Quantum invariants and torsion.* WRT/Turaev–Viro at fixed level are aggregates. Reidemeister torsion,
     Casson–Walker and Dedekind sums are computed by Euclid recursions (newfamily 13, 31: AGG/E) and are finer than
     homotopy invariants, so they are not constant on QR classes: s(1, 15) = 91/90 but s(4, 15) = 19/90 although
     4 is a square (codex). A single normalised quadratic Gauss sum of the linking form gives the Jacobi symbol (AGG);
     the full Witt class retains local (primary) information.
  3. *Reciprocity.* Theorem RA: reciprocity invariants of (N) are aggregates. The local invariants of linking forms
     (Kawauchi–Kojima) are Legendre symbols at each r \| N.
- **Obstruction encountered.** Isometry of rank-1 forms over Z/N needs the local symbols at p and q individually;
  the computable *homotopy* invariants we know (normalised Gauss sums) see only their product, and the computable
  finer invariants (torsion, Dedekind sums) are not homotopy invariants. **OPEN = QRP.** The lens-space view adds
  one question (runner-up in §4): a square-invariant Euclid-computable function that separates Jacobi fibres
  uniformly would decide QR.

### 2.4 Necklace and Burnside counts (Pillai Σ gcd)

- **Facts.** N·Neck(N,k) = Σ_{g∈Z/N} k^{gcd(g,N)} = Σ_{d\|N}φ(d)k^{N/d} = k^N + (p−1)k^q + (q−1)k^p + φ(N)k.
  - Prop NECK: k = 1 + R, modulus R², R > 4N gives Pillai (2p−1)(2q−1). One query, FACT.
  - **codex (blind round):** fixed k = 2 already suffices. Modulo 2^{p} only the d = N term survives, so
    N·Neck(N,2) ≡ 2φ(N) (mod 2^{min(p,q)}). Simplest use (proof critic): reduce mod 2^{n+1}; if p > n + 1 (else trial
    division), N is invertible mod 2^{n+1} and φ(N) mod 2^n = φ(N) exactly. So the *standard* two-colour necklace
    count, reduced mod 2^{n+1}, factors N deterministically (row I7; identity checked at N = 143, 221, 323 by codex).
    (With only 2^{n/4+O(log n)} one needs k ≥ n/4 + v₂(q − p) for the Coppersmith route; the worst case fails.)
  - Fixed k and a fixed small modulus M: (p−1)k^q + (q−1)k^p mod M is a function of (p, q) mod lcm(M, ord_M k): fiber tier.
- **Compression attempts.**
  1. *Transfer matrix.* Words of length N are closed walks of length N in the complete graph; tr(T^N) = k^N is a
     poly-time N-exponent gate. Necklaces need the orbit structure under rotation, i.e. tr(T^{N/r}) for r \| N:
     **p-exponent (PRIM).**
  2. *Generating function.* ∏_n(1 − x^n)^{−M(n,k)} = 1/(1 − kx) (cyclotomic identity) and Σ_n Neck(n,k)x^n =
     −Σ_d (φ(d)/d) log(1 − kx^d). Extracting the coefficient of x^N costs N (power series to degree N). A D-finite
     closed form would give the generic Õ(N^{1/2}) upper bound (Chudnovsky–Chudnovsky / BGS); some D-finite
     sequences have much faster isolated evaluation, so this is not a lower bound.
  3. *Lattice points.* Pillai = #{(x,y) ∈ [0,N)²: N \| xy}: the lattice points on the N hyperbolas xy = mN. Fast
     lattice-point counting (Dirichlet divisor problem, Õ(N^{1/3})) counts regions, not this arithmetic subset.
  4. *Ramanujan expansion.* Σ_g w(gcd(g,N)) expanded in Ramanujan sums c_c(N): terms with c < p see gcd(c,N) = 1
     and are blind (nonsmooth §4.1 for σ). **Obstruction encountered: Möbius inversion over the divisors of N,
     equivalently a p-exponent.**
- **Verdict.** As a target: s (FACT tier) or fiber bits (fixed small modulus; for gcd(k, M) > 1 use eventual
  periodicity with a transient). Not OPEN beyond s.

### 2.5 Cat-map orbit counts

- **Facts.** Prop CAT: Prim_t(N) = V_N(t) − V_p(t) − V_q(t) + V_1(t) (for N = pq); at t = 2 + R mod R² it gives
  J₂(N). For fixed t it gives V_p(t) + V_q(t) mod M for each modulus: a fiber function of p mod the period of V
  mod M. The *total* periodic count Fix(A^N) = V_N(t) − 2 is poly-time (Lucas powering) and is an aggregate.
- A different invariant with the same name: the number of orbits of A acting on the finite torus (Z/N)². For
  A irreducible modulo **both** p and q it is 1 + (p²−1)/o_p + (q²−1)/o_q + (p²−1)(q²−1)gcd(o_p, o_q)/(o_p o_q), with o_r = ord(A mod r) |
  r ∓ 1, **provided p, q ∤ t² − 4** (PROVED, elementary: orbits of a product action; without that hypothesis A mod r
  is ±unipotent and the formula fails, e.g. t = 3, N = 35: 59 orbits vs 37 — proof critic). This is family M data; we did not settle whether one
  such count factors N (not claimed FACT).
- **Compression attempts.** Transfer matrix = A itself: A^N is poly-time, A^p is the target. Zeta function
  ζ_A(z) = exp(Σ Fix(A^m) z^m/m) is rational (Artin–Mazur/Manning), so the total periodic data are trivial; it does
  not give efficient primitive extraction, which is the Möbius transform over divisors.
- **Obstruction encountered.** p-exponent (PRIM; leakmap NE1 classifies N-exponent gates: only family M).

### 2.6 Tiles of Z_N

- A nontrivial tiling A ⊕ B = Z_N has |A|·|B| = N, so writing a tile down exhibits a divisor: RESTATEMENT.
  Deciding whether Z_N has a nontrivial tiling is compositeness (AKS). Producing, counting and deciding are
  different tasks; counting tilings gives integers ≥ q^{p−1} (complete residue systems), i.e. ≥ Θ(p log q) bits —
  exponential in n, though not shown to be Θ(N) bits. Coven–Meyerowitz's conditions (T1, T2) are stated through the
  cyclotomic divisors Φ_d, d \| N, of the tile's mask polynomial: they need divisors (conn_top: sparse probe sets are
  blind to cyclotomic divisors, Lam–Leung).
- Generalisation: tiling a finite region by given tiles is NP-complete; deciding whether one finite tile tiles a given
  periodic subset of Z^d (d in the input) is undecidable (Greenfeld–Tao, JEMS 2025; monotiling of Z^d itself is open). **Deprioritised** per the brief.

---

## 3. Meta-question: a common reason, a model, and the escapees

### 3.1 Prop TRI (quasi-polynomial local data: AGG / FIBER / FACT). PROVED (routine extension of nonsmooth §5.3 and of Bach–Miller–Shallit, SIAM J. Comput. 15 (1986), which prove σ_k (k ≠ 0) and r₄ factoring-equivalent for all N)

*Setting.* m ≥ 1. For each unit u mod m a polynomial c_u(x) = Σ_i c_{u,i} x^i ∈ Q[x] (clear denominators when needed). F(r) := c_{r mod m}(r) for
primes r ∤ m, and f(pq) := F(p)F(q). For a nonzero c_u let d(u) ≤ D(u) be its lowest and highest degree.

*Statement.* Let u = p mod m and v = q mod m.
- (a) **FACT pair.** If c_u, c_v ≠ 0 and not (d(u) = D(u) = d(v) = D(v)), then {p, q} is computable in deterministic
  polynomial time from (N, f(N)), polynomial in n, the maximal degree, the coefficient sizes and the size of the
  explicit residue table (m fixed, or the full table of φ(m) polynomials counted in the input; a succinct
  description with binary-encoded variable m is a different problem). Integer roots of the degree-≤ 2D polynomial
  are found by exact Sturm sequences plus bisection, polynomial in its bit size for bounded degree (any standard
  polynomial-time real-root isolation, e.g. Collins–Akritas, also works).
- (b) **FIBER pair.** If d(u) = D(u) = d(v) = D(v) = k, then f(N) = c_{u,k}c_{v,k}N^k: a function of N and of the
  fiber point (u, v) over N mod m.
- (c) **AGG.** The following are equivalent:
  (i) there are Q_w ∈ Q[x], w ∈ (Z/m)^*, with f(pq) = Q_{pq mod m}(pq) for all primes p ≠ q coprime to m (f is a
  quasi-polynomial in N, identically on all prime pairs);
  (ii) for each w, either every pair (u, v) with uv ≡ w is ZERO, or every such pair is FIBER with one common k and one
  common value c_{u,k}c_{v,k};
  (iii) F ≡ 0 on the units, or F(r) = c·χ(r)·r^k with c ≠ 0, one k, and χ a real (±1-valued) character mod m.
  (Note: "f(N) is a function of N" is vacuous, since N determines {p, q}; (i) is the meaningful notion.)

*Proof.* (a) Put q = N/x. For each of the ≤ φ(m) pairs (u′, v′) with u′v′ ≡ N (mod m), form
G(x) = x^{D(v′)}(Σ_{i,j} c_{u′,i}c_{v′,j}N^j x^{i−j} − f(N)) ∈ Z[x]. For the true pair the coefficient of
x^{D(u)−d(v)+D(v)} is c_{u,D(u)}c_{v,d(v)}N^{d(v)} ≠ 0 (it is the only pair (i, j) with i − j = D(u) − d(v)), and
likewise for x^{d(u)−D(v)+D(v)}; if not all four degrees coincide, one of D(u) − d(v), d(u) − D(v) is nonzero, so G is
not a constant multiple of x^{D(v)}, hence G ≠ 0 with p among its ≤ deg G integer roots. Isolate integer roots
exactly (Sturm) and test divisibility. (b) Immediate. (c) (i)⇒(ii)/(iii): each residue class contains infinitely
many primes (Dirichlet), so the prime pairs in classes (u, v) are Zariski-dense in the plane and c_u(x)c_v(y) =
Q_{uv}(xy) in Q[x, y]. Only diagonal monomials x^i y^i can occur, so supp c_u = supp c_v = {k}, with one k for all
units. Zero is all-or-nothing: if a_u := c_{u,k} = 0 then Q_{uv} = 0 for every v, so a_v² = 0 for all v. Otherwise
g(u) = a_u/a_1 satisfies g(u)g(v) = g(uv) (compare (u, v) with (uv, 1)), so g is a homomorphism (Z/m)^* → Q^*; its
image is a finite subgroup of Q^*, i.e. ⊆ {±1}: a real character. (iii)⇒(i): c·χ(r)r^k gives c²χ(N)N^k. (ii) is (iii)
restated pairwise. The proof of (i)⇒(iii) survives removing finitely many prime pairs (if P vanishes on S×T∖E,
multiply P by a polynomial vanishing on the finite set E), so (i) may be weakened to "all but finitely many". ∎

*Experiment E2 (pre-registered, all predictions held).* 15 local-data families × 120 seeded semiprimes (n = 32, 64,
128): every FACT-tier family split 120/120 (φ, ψ = σ, Pillai, k(SL₂), J₂ = k(GL₂) = #K₃, k(Heis), subgroups of Z_r²,
unit circle r − χ₄(r), #K₅, (r − 1)²); the four FIBER families (1 + χ₄, 1 + (−3/r), r(1 + χ₄), ε-type 1/2) gave 0
splits; the AGG family χ₄(r)r² matched χ₄(N)N² on 120/120. Total 4.8 s, 29 MB.

*Reading.* Prop TRI is why cluster I (21 of 55 rows, ≈ 40%) "collapses" to s. Geometry's three tiers (FACT, QR,
aggregate) are instances of (a), (b) and (c); (b) is "fiber bits over N mod m", which is QR-type only when the
coefficients are quadratic characters (in general it is any function of (p mod m, q mod m)). AGG is a global
property: AGG data are FIBER pairs whose products depend only on uv. The FIBER tier is where ε, r₂, ν₂, ν₃,
#X₀(N)(F_ℓ), φ(N) mod ℓ and the homotopy count of lens spaces live; by Theorem PT-m a pure tree learns there at
most the fiber floor. Prop TRI says nothing about routes other than evaluating F(p)F(q) as a polynomial.

### 3.2 Prop LOC (generic local data force a factor query). PROVED (elementary)

*Model.* The algorithm gets N = pq and oracle access to local data F: {primes} → Z/M (M prime), drawn uniformly
at random; it may do arbitrary computation on N. Target: f(N) = F(p)F(q).

*Statement.* With h = Pr[the algorithm queries p or q] (event Q), Pr[output = F(p)F(q)] ≤ h + (1 − h)(2M − 1)/M²
(codex's sharpening of v1's h + (2M − 1)/M²). Hence h ≥ (Pr[success] − β)/(1 − β), β = (2M − 1)/M².

*Proof.* On ¬Q the run, and the event ¬Q itself, are functions of the coins and of F away from {p, q}, so they are
independent of (F(p), F(q)). If neither p nor q is queried, (F(p), F(q)) is uniform on (Z/M)² and independent of the view; the most likely
product is 0, with probability (2M − 1)/M². ∎

*Reading.* (M must be prime: for M = 4, Pr[ab ≡ 0] = 1/2 > 7/16.) Equivalently the algorithm queries p or q
with probability ≥ Pr[success] − (2M−1)/M², which is meaningful for M large or success bounded above 2/M. In the
point-query oracle model, successful evaluation of a generic multiplicative function exposes a factor query. It
says nothing about structured functions evaluated through identities, generating functions or summatory
recursions. The model cannot give an
N^c lower bound unless factoring itself has one: it allows every factoring algorithm. Any escape must exploit
*structure* in F. Prop TRI says which polynomial structures help (only c·χ·r^k), and the escapees (§3.4) are the
structures outside its scope.

### 3.3 Prop AUT (no digit automaton for cluster I on all integers). A corollary of KNOWN results (Schlage-Puchta, J. Integer Seq. 6 (2003) Cor. 3; Konieczny, arXiv:1905.11981 Thm 1.1: multiplicative automatic sequences agree off finitely many primes with a Dirichlet character or 0; with Allouche–Shallit, TCS 98 (1992) Cor. 2.4); our Wirsing/Berstel proof is an alternative

*Statement.* Let f be multiplicative, integer-valued, with f(r) = F(r) at all but finitely many primes for a polynomial F ∈ Z[x], and
suppose some prime ℓ has 0 < #{u ∈ (Z/ℓ)^*: F(u) ≡ 0} < ℓ − 1. Then f is not k-regular for any k ≥ 2 (no Z-weighted
automaton over the base-k digits computes f(n) for all n). Examples: φ (ℓ = 3), σ (ℓ = 3), Pillai (ℓ = 3),
J₂ (ℓ = 5), k(SL₂) (ℓ = 3).

*Proof.* If f were k-regular, f mod ℓ would be k-automatic (Allouche–Shallit 1992, Cor. 2.4), so A = {n: ℓ ∤ f(n)} would be
k-automatic. Its indicator is multiplicative and equals 1 on a set of primes of relative density 1 − c with
c = #roots/(ℓ − 1) ∈ (0, 1) (Dirichlet). By Wirsing's theorem #A ∩ [1, x] ~ C x(log x)^{−c} with C > 0 (the indicator is 0/1 multiplicative, mean 1 − c > 0 on primes), so the number of
members with exactly n digits is ~ C′ k^n n^{−c}. For a regular language the number of words of length n is, on each residue class of
n modulo some period, either eventually 0 or ~ C n^d λ^n with d ∈ Z_{≥0} (rational generating functions with
nonnegative coefficients; Berstel's theorem on dominant poles). λ = k with exponent n^{−c}, c ∉ Z, is impossible;
λ < k is too small. ∎

*Scope.* This rules out automata that are correct on *all* n. Automata correct only on semiprimes are not
excluded, nor are non-multiplicative interpretations of cluster-I targets; that question is newfamily row 53 /
leakmap's automaton class, OPEN (S2). Prior art: φ is not k-regular (Coons, JTNB 22 (2010) Thm 3.3); σ_m, φ, μ mod v are not automatic (Yazdani, JTNB
13 (2001)); the general multiplicative criterion is Schlage-Puchta 2003. Regular-language asymptotics: Berstel–
Reutenauer, *Noncommutative Rational Series with Applications* (2011); Flajolet–Sedgewick Thm V.3 (numbers UNVERIFIED). J₂ needs ℓ = 5 (at ℓ = 3 both units are roots); k(SL₂) = r + 4 fails only at r = 2, hence "all but finitely many primes".

### 3.4 The common reason, stated carefully

What the inventory shows (a non-exhaustive taxonomy of implementations, not a theorem; no general native-time lower
bound is proved anywhere in this note):
1. **Every definition-level algorithm is a sum or a product over a set of size N^{Θ(1)}** (residues, divisors,
   group elements, cells, Hecke spaces).
2. **The compressions we found, and what limits each:**
   - CRT / multiplicativity: needs the factors (Prop LOC for generic local data; Prop TRI: for quasi-polynomial data the value itself factors N);
   - N-exponent gates (x^N, ψ^N, A^N, Frob at N): leak only family M (leakmap NE1, prim);
   - pure mod-N ring algebra: Theorem PT;
   - reciprocity / Euclid along a/N: aggregates or class E (Theorem RA; Prop DSF shows the FACT part of a Dedekind
     sum is its local part, which reciprocity does not reach);
   - smoothness relations in a global group (class groups): subexponential, L[1/2] = N^{o(1)}, so *below* every N^c;
     it is an index calculus, i.e. a factoring algorithm of the familiar kind, slower than NFS;
   - weighted lattice-point summation and quotient recursions (S, its hull refinement, Φ, F_χ): no barrier; they stop
     at N^{1/3} (unweighted; polynomial weights by the SKETCH of §1.3) or N^{2/3} (Möbius weights) as far as known (S1);
   - compressed integer linear algebra (Smith's determinant, small chain complexes): no barrier beyond the fixed-size
     obstruction of §2.1 (S3).
3. **A dichotomy that organises the escapes (HEURISTIC principle).** Fourier coefficients at index N of Hecke
   objects (σ, r₄, a_N, dimensions, point counts) are multiplicative: FACT tier or fiber tier. Values at the cusp
   a/N (Dedekind sums, modular symbols, quantum invariants of lens spaces, Gauss sums) have SL₂(Z)-reciprocity:
   computable, but aggregates or class E. The bridge between them is the circle method (Rademacher, Petersson):
   coefficients are sums over cusps a/c of Kloosterman sums of N mod c, and the information about p enters only at
   c ≥ p (nonsmooth §4.1, proved for the Ramanujan expansion of σ). This is the precision wall in another form.
4. **Lower bounds of the PT kind are out of reach.** Precisely (codex_r2_generic_ring §3, [proved] there): a
   superpolynomial lower bound, in arithmetic operations, against **all nonuniform, modulus-independent,
   constant-free, division-free SLPs followed by one gcd with N** with constant success on balanced prime distributions would imply that τ(m!) is
   not polylogarithmic (a short SLP for a factorial cutoff near the median prime would split with probability
   ≈ 1/2), and hence, by Bürgisser (ECCC 2006 Thm 1.1; Comput. Complexity 18 (2009) Main Thm 1.2), superpolynomial
   constant-free, division-free permanent complexity. Uniform-only, bit-cost or worst-case bounds do not
   automatically give this implication. The only model
   theorems available are the restricted ones above: Prop LOC (generic local data), Prop AUT (all-n automata), the
   Ramanujan-expansion blindness (nonsmooth), and Theorem PT itself (split-or-nothing, not a time bound).

**Escapees** (native algorithms not of the CRT / N-exponent / PT / reciprocity type):
- **E1. Non-abelian Frobenius local data** (cluster III): a_N(E), τ(N), #E(Z/N) for one curve. F(r) = a_r is not
  quasi-polynomial, so Prop TRI does not apply; Prop FROB covers only the operator. `nativemod`.
- **E2. Global non-multiplicative invariants** (cluster IV): h(−4N), r₃(N), Jacobi prefix sums, Hurwitz sums. Their
  native algorithm is an index calculus (L[1/2]) — the one escape that historically produced factoring algorithms.
  IV2 is the cleanest form: each summand (a/N) is poly-time by reciprocity, the prefix sum Σ_{a<N/2}(a/N) is a class
  number. `nativeana`.
- **E3. Root-type sums** (cluster V1): Salié S(a,b;N) = (Gauss factor)·Σ_{x²≡ab (N)} e(2x/N) is the archimedean
  shadow of the four square roots; at prime modulus it is computed by Tonelli–Shanks, at composite N nothing
  better than N terms is known. Not covered by PT or PRIM (family A).
- **E4. Reciprocity-class functions** (class E): computable, but no FACT member is known (classE: no leak).
- **E5. Weighted lattice-point summation and compressed linear algebra** (algorithm S and its hull refinement;
  Smith's determinant): outside PT and outside CRT; the best known is N^{1/3} (SKETCH for σ).

---

## 4. Shortlist

Screen: (i) the native question is genuinely open, (ii) no NP-hardness result is known for the natural
generalisation (for S1 this is "no hardness known", not a proved ≤F screen: summatory functions are not in FP^FACT
as far as known), (iii) nothing known rules out a polylog (or sub-L[1/3]) native algorithm. None of the three is a
new route; each is a precise sub-question that the views in this note make sharper.

**S1. Exact summatory functions with Möbius weights (summatory totient; twisted totient giving h(−7N²)).**
- Why it passes: h(−7N²) = F_χ(N) − F_χ(N−1) with F_χ(x) = Σ_{m≤x} m∏_{ℓ|m}(1 − χ(ℓ)/ℓ), χ = (−7/·) of fixed
  conductor (codex); φ(N) = Φ(N) − Φ(N−1). These are native recursions on the quotient states ⌊x/j⌋, not CRT
  evaluations, so neither Prop LOC nor Theorem PT applies. Best known Õ(x^{2/3}). For the unweighted and
  polynomial-weight analogues (D(x), H₁(x)) the lattice-hull method already gives Õ(x^{1/3}) (Sladkey; codex's
  weighted extension, §1.3, SKETCH), so native σ(N) should cost Õ(N^{1/3}).
- First question: do Φ(x) and F_χ(x) admit a hull/Euclid-type evaluation in Õ(x^{1/3}) (closing the gap to σ), and
  which quotient states survive the difference F(N) − F(N−1)? Any exponent below 1/4 would beat Pollard–Strassen;
  below 1/5, Harvey–Hittmeir.
- Strongest reason to fail: the difference at N is a single-integer jump, and every known summatory speed-up below
  x^{1/3} is a smoothing whose correction needs the arithmetic of integers near N (conn_ana's de-smoothing window;
  Li 2026). Möbius weights already cost x^{2/3}. Deep dive: `nativeana`.

**S2. Digit automata on semiprimes for a fixed-modulus fiber function.**
- Why it passes: for fixed ℓ, φ(N) mod ℓ (or σ₃(N) mod ℓ) is a Prop TRI fiber function (O(1) bits; not FACT
  alone), but over ℓ ≤ O(n) the residues CRT to s. A uniformly constructible automaton on the base-k digits of N,
  correct on semiprimes, with poly(ℓ) states, would factor N. Prop AUT rules out automata correct on *all* integers;
  leakmap's AC⁰ theorem does not cover automata (parity is regular).
- First question: is there a k-automatic g with g(pq) = (p − 1)(q − 1) mod 3 for all primes 3 < p < q? On the
  promise this is [p ≡ q ≡ 2 (mod 3)] — for N ≡ 1 (mod 3), the bit "N is not of the form p·q with p ≡ q ≡ 1 mod 3",
  i.e. the (−3/p) fiber bit. **Partial negatives (codex, PROVED; checked by the round-2 critic):** no eventually
  periodic g works (take large prime pairs ≡ 1 and ≡ −1 mod lcm(3, M)); no multiplicative automatic g works
  (Konieczny's classification, applied to g viewed in {0, ±1} ⊂ C; checked against the abstract only). Konieczny
  is genuinely needed: g(n) = [every prime factor of n is ≡ 2 mod 3] is multiplicative and matches on all semiprimes. The arbitrary
  extension is CONJECTURED impossible. Route: Möbius disjointness for automatic sequences (Müllner, Duke 166 (2017))
  and prime number theorems for automatic sequences (Mauduit–Rivat; Adamczewski–Drmota–Müllner), extended to a
  bilinear (semiprime) statement — not available off the shelf.
- Strongest reason to fail: it is the fiber-bit problem (Φ-hiding-type) restricted to a weak model; a NO would be a
  new barrier, not an algorithm; a YES would be astonishing. One fixed modulus would not by itself factor.

**S3. Compressed integer linear algebra for φ(N) (small chain complexes, Smith's gcd determinant).**
- Why it passes: explicit linear algebra is in P (Kannan–Bachem), so the generalisation is not hard; φ(N) is the last
  Schur complement of Smith's matrix G_N = E·diag(φ)·E^T and the torsion order of a chain complex is a Smith-form
  output. Codex's equivalence (§2.1) shows that the unrestricted question *is* deterministic factoring, so the
  research content is in restricted classes. Its fixed-size obstruction (entries in Z[N]) is PROVED.
- First question: extend the obstruction to matrices of size poly(n) with entries in Z[N] of bounded degree, and
  then to entries that may also use ⌊N/2^i⌋ (digits). Does a poly(n)-size Z[N]-matrix family exist whose torsion
  order is φ(N) on all semiprimes? (The fixed-size argument used finitely many minors and one Bézout identity.)
- Strongest reason to fail: with growing size, Smith forms can simulate gcd computations, and the target is just
  factoring; an obstruction for poly(n)-size Z[N]-matrices would be a new restricted lower bound, not an algorithm.

Runner-up (corrected from v1's S3). **Euclid-computable homotopy invariants of lens spaces.** v1 conjectured that
every square-invariant class-E function is a function of the Jacobi symbol and N mod m. Codex refuted this as
stated: H(N, a) = (a/5) if 5 | N, else 0, is square-invariant and class E (a bounded coboundary reads the terminal
CF matrix mod 5), and separates a = 1, 3 at N = 35. Corrected question: with gcd(N, m₀) = 1 for the level primes
m₀ and finitely many exceptional N excluded, is there a precisely defined class-E function (terminal terms and
postprocessing specified) that separates a Jacobi fibre for **every** admissible N? (Then QR follows by comparing
a, 1 and ab, b with (b/N) = −1.) Known: raw Dedekind sums are not square-invariant (s(a,N) ≤ s(1,N) with equality
only at a ≡ 1 mod N, by the rearrangement inequality; so not square-invariant once N ∤ 24); the standard
square-invariant phase 12N·s(a,N) ≡ N + 1 − 2(a/N) (mod 8) (N odd, gcd(a, N) = 1; Rademacher–Grosswald; verified
for all odd N < 80) is Jacobi.

Cross-reference (not mine): the single most natural open native question in the whole inventory is a_N for one
fixed non-CM curve (cluster III), which `nativemod` owns.

---

## 5. Review record

- **codex blind round** (`codex/nativeinv_1.out`): independent list of 36 invariants; contributed algorithm S, the
  binary-necklace residue (I7), the ring class number h(−7N²) and its twisted summatory recursion (I19), Smith's gcd
  determinant (I20), topological φ/J₃ counts (I21), the shear (I6) and squaring (II4) counts, the unitary-Cayley
  chromatic number (VII15), the E₄ non-D-finiteness argument, and its own shortlist (σ₃ digit kernel, h(−7N²)
  recursion, elliptic twist counts mod ℓ). Same verdict: many formulas, few targets, no native algorithm beating
  factoring, no lower-bound theorem in an unrestricted model.
- **Opus proof critic** (ACCEPT WITH FIXES): Prop TRI(c) was vacuous as written (restated as quasi-polynomial
  identity on prime pairs, with the missing direction proved; characters are real); the cat-map orbit formula
  needed p, q ∤ t² − 4 (counterexample N = 35, t = 3); gcd-sum class formulas need N squarefree (N = 9 fails); the
  necklace reduction needed 2^{n+1} (or v₂(q − p) in the Coppersmith form); oriented vs unoriented homotopy;
  randomized-reduction caveat in the ≤F screen; several row fixes. All applied.
- **Opus literature critic** (OK WITH FIXES): Prop AUT is known in substance (Schlage-Puchta 2003; Coons 2010;
  Yazdani 2001); the "Pak–Soffer universality" claim was a misattribution (removed; Halasi–Pálfy for non-polynomial
  pattern groups); Greenfeld–Tao statement corrected; Kayal–Saxena example corrected; Kahn 1982 for characteristic
  sets; Bach–Miller–Shallit cited for σ_k, r₄. All applied.
- **codex adversarial round** (`codex/nativeinv_2.out`, MAJOR REVISION): lower-bound language removed (§0, §3.4: a
  description of implementations, not a theorem; L[1/2] class-group algorithms are below every N^c); "forcing step"
  replaced by "obstruction encountered"; weighted hull extension gives native σ(N) in Õ(N^{1/3}) (S1 rewritten);
  Schur-ring claim restricted to radial summands (counterexample t = 1, 2 at N = 15); lens-space invariant claims
  corrected (Witt class is finer than the Jacobi symbol; torsion/Dedekind sums are not homotopy invariants,
  s(1,15) = 91/90 ≠ s(4,15) = 19/90); v1's S3 conjecture refuted by H(N,a) = (a/5)[5 | N] and moved to a corrected
  runner-up; small-homology equivalence and fixed-size obstruction added (new S3); S2 partial negatives added;
  inventory fixes (III4, IV1, IV2, IV3, V3, VII1, VII3, VII5, I16, I18, III3, II2, II4); permanent-barrier
  quantifiers restored; aliases flagged (row counts are not target counts). We re-checked numerically: the
  mod-8 Dedekind formula (8 values of N), s(1,15), s(4,15), and Σ_{a<N/2}(a/N) = (2 − (2/N))h(−N) for 9 values of
  N ≡ 3 (mod 4).
- **Opus proof critic, round 2** on the codex-supplied claims (no FATAL): weighted hull relabelled SKETCH; small-
  homology obstruction repaired (ρ×ρ minors for generic rank ρ, finitely many exceptions, h the Z[N]-gcd via Gauss's
  lemma, TRI(c) finite-exception remark added); S2 negatives and the runner-up counterexample confirmed; mod-8 and
  rearrangement facts verified for odd N < 80; LOC independence step written out; permanent-barrier model
  description aligned with codex_r2 ("division-free, followed by one gcd"). All applied.
