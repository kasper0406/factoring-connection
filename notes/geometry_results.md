# geometry: geometric constructions from N, and where the factors hide (phase 26)

Agent `geometry`, 2026-10-24. **Version v1.1 (codex blind round integrated; before the codex adversarial round and the Opus critics).**
- Pre-registration: `notes/geometry_prereg.md` (commit `5ae824e`, before any code ran; amendment A1 for E4 committed
  before E4 code was written).
- Code: `src/geometry/{modcurve,curves,harness,decoders,run_experiments}.py`. Tests: `tests/test_geometry.py`.
- Raw data: `experiments/raw/geometry_e{1,2,3,4}.jsonl`.
- GPT-6.1-sol: blind round `codex/geometry_1.{md,out}`; adversarial round `codex/geometry_2.{md,out}`.
- Labels: **KNOWN** (with reference), **PROVED** (proof here; routine unless stated; "no prior formulation found by
  us" is the strongest novelty claim), **HEURISTIC**, **experimentally supported**, **OPEN**, **UNVERIFIED**.
  Negative literature statements are from our and codex's searches and are UNVERIFIED as negatives.

Notation: N = pq with p < q odd primes unless stated; s = p + q; ψ(N) = N∏_{r|N}(1 + 1/r); ω(N) = number of
prime factors; ε = (−1/p) on N ≡ 1 (mod 4). "FACT" = randomized-polynomial-time equivalent to factoring N.
"QR-type" = a function of N and of Legendre symbols (D/p), (D/q) at small fixed D; deciding these is quadratic
residuosity modulo N for small arguments, which is not known to be easy and not known to be as hard as factoring.
"Aggregate" = a function of N mod m for a fixed m (`notes/factorbits_results.md` §3.1(A)).

---

## 0. Verdict (plain language)

**Yes, there are many genuine connections between geometry and factorisation, and several are beautiful. One of
them is a real algorithm, already known: Lenstra's elliptic-curve method (ECM, 1987) factors by doing geometry
(the group law of an elliptic curve) modulo N, and elliptic curves also prove primality (Goldwasser–Kilian,
Atkin–Morain). The number field sieve is geometry of numbers. Beyond these known methods, no geometric
construction we examined gives a faster way to factor, and we can say exactly why in each case.**

The key point first. Primality itself is already easy: it can be decided in polynomial time (AKS 2002). So a
geometric picture of primality only helps with *factoring* if three things hold at once:
1. the geometric object can be built from N alone, without knowing the factors;
2. some feature of it encodes the factors (or even one bit about them);
3. that feature can be computed geometrically (by counting, linear algebra on a triangulation, volume, fixed-point
   formulas, …) faster than factoring.

Every construction we examined fails one of the three, in one of four recurring ways:
- **Trivial restatement.** "N is prime iff the only integer rectangle of area N is 1 × N", "iff every star polygon
  {N/k} is connected", "iff the lens space L(N, 1) has no connected covers except the 3-sphere". These are
  definitions in disguise; the geometry adds no algorithm.
- **The feature encodes the factors, so computing it is as hard as factoring.** The cleanest example: the
  **modular curve X_0(N)**. It is defined from N alone. Its hyperbolic area is (π/3)·ψ(N), and for N = pq,
  ψ(N) = (p + 1)(q + 1) = N + (p + q) + 1. So the area gives p + q, and p + q with pq = N gives p and q. The genus
  of X_0(N) does the same (up to four candidates). But every way to compute the area or genus either uses the
  formula (which needs the factors) or builds the surface explicitly. An explicit surface has at least about
  N/3 triangles, because a surface of genus g needs ≥ 4g − 4 triangles and g ≈ N/12. That is exponential in the
  number of digits of N.
- **The feature is computable, but it only sees "local" information about N** — N modulo a fixed number, or
  quadratic-residue symbols. Lens-space invariants, homology of cyclic covers with fixed coefficients, counts of
  covers of bounded degree, quantum invariants at fixed level and polytope counts with fixed shape are of this kind.
- **The feature sees a little more, but only quadratic-residue information.** Again X_0(N) gives the sharpest
  example. Its number of points over a small finite field F_ℓ is a weighted count of square roots of the numbers
  t² − 4ℓ modulo N, so it reveals Legendre symbols (D/p) for |D| ≤ 4ℓ and nothing else. Its *cone points* do the same:
  X_0(N) has 4 cone points of angle π if −1 is a square mod N and none otherwise.

**The sharpest observation (PROVED from classical formulas; verified 800/800; the F_4 case was found first and
independently by GPT-6.1-sol from the moduli interpretation).** For X_0(N) the line between
"quadratic symbols only" and "as hard as factoring" falls exactly between F_ℓ and F_{ℓ²}:
- #X_0(N)(F_ℓ) is QR-type: it is a function of 2^{ω(N)} and the symbols ((t² − 4ℓ)/p), ((t² − 4ℓ)/q).
- #X_0(N)(F_{ℓ²}) is factoring-equivalent: it equals (ℓ − 1)ψ(N)/12 plus bounded QR-type terms. A single number,
  the count of points of X_0(N) over the 4-element field, factors N.
- The geometric reason: over F_{ℓ²} the supersingular elliptic curves have Frobenius equal to the scalar ±ℓ. Every
  one of the ψ(N) cyclic subgroups of order N is then rational, so the whole fibre of X_0(N) → X(1) above each
  supersingular point is counted. Over F_ℓ, Frobenius is never scalar (4ℓ is not a square), and only the subgroups
  fixed by a non-scalar Frobenius count, i.e. square roots of t² − 4ℓ mod N.
- The same split holds for X_1(N): it has exactly φ(N)/2 points over F_ℓ once N > (√ℓ + 1)², all of them cusps, so
  #X_1(N)(F_2) alone factors N.

**The ladder of lattice points on spheres** shows the same three tiers (cross-referencing `analytic`):
- the circle x² + y² = N: r₂(N) ∈ {0, 16} is one QR bit; two essentially different points factor N; finding one
  point is equivalent to finding √−1 mod N (OPEN);
- the 2-sphere x² + y² + z² = N: r₃(N) = 12h(−4N) (Gauss), a class number. A known class number usually leads to
  an ambiguous form and a factor, though one value is not known to factor every N;
- the 3-sphere: r₄(N) = 8σ(N) = 8(N + p + q + 1) (Jacobi) factors N at once.

**General principle (a proved statement in each model, not a theorem about all geometry).** Invariants of a compressed
object built from N fall into exactly the tiers the project already knows:
- invariants computed through a fixed finite target (homomorphisms to a fixed group, homology with F_ℓ
  coefficients, covers of bounded degree) depend only on N mod m (Prop AGG-H below);
- local fixed-point counts give quadratic symbols (Lefschetz-type data → QR);
- global counts of a size-N set (area, Euler characteristic, number of points of P¹(Z/N), Burnside counts,
  primitive periodic orbits) are multiplicative or Möbius sums over the divisors of N and are factoring-equivalent
  (Props MF, NECK, CAT);
- cycle structures and periods are orders, hence factoring-equivalent by Miller/Shor-type reductions.

**What remains OPEN** (§8): (i) whether the counts #X_0(N)(F_ℓ), ℓ ≤ L, for a growing L ≈ (log N)^{O(1)} — that is,
the Legendre symbols of p at all small discriminants — let one find p in subexponential time. No method is known;
this is the project's QR/genus question in geometric dress. (ii) the "linking number" of the two primes of N in
the arithmetic-topology dictionary, which is a Rédei symbol (4-rank of the class group of Q(√N)); no algorithm from
N alone is known.

**Bottom line for the user.** Geometry is an excellent *language* for factorisation: divisors are rectangles,
sublattices, covers and fixed points; primes are the "knots" of arithmetic topology; the number field sieve is
geometry of numbers. But in every construction we know, the geometric feature that sees the factors is either
defined through them or costs as much as the arithmetic it re-expresses. We found no route to faster factoring, and
none of the candidates below is a live lead.

---

## 1. Sorting principle

A construction is a map N ↦ X_N (a polygon, surface, 3-manifold, variety, group, dynamical system) with a
description computable from N in poly(log N) time, plus a feature F(X_N). If F(X_N) is computable from the
description in time T, then T bounds the cost of F as a function of N. So "a geometric feature that encodes the
factors and is fast to compute" is *literally* "a fast factoring algorithm"; the geometry contributes only if it
suggests one (codex calls this the composition principle; it rules out nothing unconditionally). Conversely,
artificial positive examples are trivial: run AKS and output S¹ × S² if N is prime and S³ otherwise. The useful
questions are therefore (a) which tier F lies in, and (b) where exactly a geometric
computation of F pays size N^{Ω(1)} or needs the factors.

Tiers (the project's taxonomy: `notes/factorbits_results.md` §3, `notes/leakmap_results.md` §0,
`notes/newfamily_results.md` §0):
- **R** trivial restatement (no algorithmic content);
- **AGG** aggregate (N mod m);
- **QR** quadratic-residue type;
- **E** class E (continued-fraction / Dedekind-sum statistics of a/N; computable, non-aggregate, OPEN);
- **ω** depends on ω(N) or d(N) only (OPEN in general, trivial for semiprimes under the promise);
- **FACT** factoring-equivalent;
- **EXP** computable geometrically only with size N^{Ω(1)} objects.

---

## 2. Candidate table

(i) built from N alone? (ii) feature and how it encodes the factors; (iii) best computation without the factors;
(iv) verdict. "§x" points to the proof below; cross-references to earlier notes are not redone.

### 2.1 Lattice points, polygons, polytopes

| # | Construction | (i) | (ii) feature → factors | (iii) cost without factors | (iv) |
|---|---|---|---|---|---|
| 1 | Rectangles of area N; hyperbola xy = N | yes | lattice points = divisors | trial division | R |
| 2 | Lattice points under xy = N in a strip, D_N(X) = Σ_{d≤X}⌊N/d⌋ | yes | D_N(X) − D_{N−1}(X) = #{d ≤ X : d \| N}; bisection | Õ(N^{1/3}) (convex-hull walk, Θ(N^{1/3}) hull vertices) | FACT, EXP (`analytic` §2; `newfamily` row 50) |
| 3 | Full Dirichlet count D(N) = Σ_{m≤N} d(m) | yes | first difference d(N) | Õ(N^{1/3}) | ω (useless for semiprimes) |
| 4 | Circle x² + y² = N | yes | r₂(pq) ∈ {0, 16} decides ε when N ≡ 1 (4); two inequivalent points factor N (Euler) | one point ⇔ √−1 mod N (Cornacchia) | QR (count); OPEN (one point); FACT (two) |
| 5 | Sphere x² + y² + z² = N | yes | r₃(N) = 12h(−4N) (N ≡ 1, 2 mod 4, squarefree; Gauss) | h(−4N): N^{1/5+o(1)} under ERH | class number: a known h(−4N) yields an ambiguous form with constant probability on typical N, but the reachable order-2 class can be the non-splitting (2, 2, ·) class for some 2-Sylow shapes, so one value is not known to factor every N (codex). FACT-likely, HEURISTIC |
| 6 | 3-sphere in Z⁴ | yes | r₄(N) = 8σ(N) = 8(N + s + 1), N odd (Jacobi) | ≍ N samples via theta extraction | FACT (`analytic` table) |
| 7 | Conic x² + y² ≡ 1 (mod N) (points over Z/N) | yes | # = (p − χ₄(p))(q − χ₄(q)) = N − χ₄(p)s + 1 or N ± (q − p) − 1 | count over Z/N | FACT (§4.1, Prop MF′) |
| 8 | Lattice polygons with poly(log N)-bit vertices (Pick); Barvinok counts in fixed dimension | yes | Pick: A = I + B/2 − 1, B = Σ gcd(edge) | poly time (Barvinok 1994) | AGG / E (`newfamily` rows 9, 49) |
| 9 | Polytopes with N-dependent normals, Θ(N^{1/3}) vertices | yes | can detect divisors | size of the hull | EXP (`newfamily` row 50) |
| 10 | Sublattices of Z² of index N (N-sheeted covers of the torus) | yes | number = σ(N) | — | FACT (§4.1) |
| 11 | Reduced lattices / flat tori of discriminant −4N (CM points on X(1)); cusps of the Bianchi orbifold of Q(√−N) | yes | number = h(−4N); 2^{ω−1} genera | O(√N) enumeration; N^{1/5} (ERH) | as #5 |

### 2.2 Modular curves and hyperbolic surfaces (§3)

| # | Feature of X_0(N) / X_1(N) | (ii) | (iii) | (iv) |
|---|---|---|---|---|
| 12 | Hyperbolic area (π/3)ψ(N) = index [SL₂(Z):Γ_0(N)] = \|P¹(Z/N)\| = degree of X_0(N) → X(1) = degree of the modular polynomial Φ_N | ψ(pq) = N + s + 1 | any explicit cell structure has ≥ ψ/6 − O(1) faces; Φ_N has degree ψ(N) | FACT (Prop X0-1), EXP |
| 13 | Genus g(X_0(N)) = 1 + ψ/12 − ν₂/4 − ν₃/3 − ν∞/2 | 4 candidates for ψ | triangulation ≥ 4g − 4 triangles; Farey symbol (Kulkarni) Θ(g + cusps); modular-symbol space dimension 2g + ν∞ − 1 | FACT, EXP |
| 14 | Cusps ν∞ = Σ_{d\|N} φ(gcd(d, N/d)) = 2^{ω(N)} (squarefree) | ω | — | ω |
| 15 | Cone points ν₂ = ∏(1 + (−1/r)), ν₃ = ∏(1 + (−3/r)) | ν₂ ∈ {0, 4} decides ε on N ≡ 1 (4) | — | QR |
| 16 | Points over F_ℓ: #X_0(N)(F_ℓ) = 2^{ω} + ½Σ_{t²<4ℓ}Σ_f h_w((t²−4ℓ)/f²)∏_{r\|N}(1 + ((t²−4ℓ)/r)) | Legendre symbols at \|D\| ≤ 4ℓ | E3: the counts for ℓ ≤ 13 determine ε and all ≈ 7 unknown symbols, 4000/4000 | QR (Prop X0-2) |
| 17 | Points over F_{ℓ²}: (ℓ − 1)ψ/12 + bounded QR terms | ψ | E4: 800/800 split from one count | FACT (Prop X0-3) |
| 18 | Points of X_1(N) over F_ℓ, N > (√ℓ+1)²: φ(N)/2 rational cusps, no other points | φ | — | FACT (Prop X1) |
| 19 | Lefschetz number of the Fricke involution w_N: #fix = class numbers h(−4N) (+ h(−N) if N ≡ 3 mod 4) (Ogg 1974; Fricke) | h | as #5 | as #5 |
| 19a | Genus of X_1(N); dimensions dim S_k(Γ_0(N)) | g₁ = 1 + J₂(N)/24 − ¼Σ_{d\|N}φ(d)φ(N/d) (N ≥ 5); dim S_{k+12} − dim S_k = ψ(N) | Θ(N²) cosets for Γ_1 | FACT (codex; Gu–Martin 2019 for dimension counts) |
| 20 | Spectrum of the Laplacian / Selberg zeta / random walks on Γ_0(N)\H or its coset graph | Weyl law: area; return counts of one word W: fixed points of W on P¹(Z/N) = ∏(1 + ((tr²W − 4)/r)) | discretisation of size ψ(N) | Weyl term FACT+EXP; return counts QR; cycle structure of W = orders, FACT (Shor-type) |
| 21 | Bianchi orbifolds Γ_0(N·O_K)\H³ | covolume × ∏_{r\|N}(r² + 1 or (r + 1)²) | Swan/Riley domains have size ∝ index | FACT, EXP |

### 2.3 Topology and arithmetic topology (§5)

| # | Construction | (ii) | (iii) | (iv) |
|---|---|---|---|---|
| 22 | Lens space L(N, q) (Heegaard genus 1, O(log N) description) | connected d-fold covers exist iff d \| N; N prime iff the only nontrivial connected cover is S³ | — | R |
| 23 | Invariants of L(N, q): H₁ = Z/N, Reidemeister torsion, Casson–Walker, d-invariants, WRT, spectrum | Dedekind sums, Gauss sums, lattice counts | Euclid recursions | AGG / E (`newfamily` rows 13, 31; Prop DSF: the local part of s(a, N) is FACT) |
| 24 | Groups/manifolds with N in exponents or surgery coefficients (L(N,q), M(N,1) on a fixed knot, mapping tori of φ^N): #Hom(π₁, G) for fixed finite G; covers of degree ≤ k; Fox colourings; Dijkgraaf–Witten with fixed gauge group | divisors appear only for target size ≥ p | poly(log N) | AGG (Prop AGG-H) |
| 25 | N-fold cyclic (branched) covers Σ_N(K) of a fixed knot K | \|H₁\| = ∏_{d\|N, d>1}\|Res(Δ_K, Φ_d)\|; dim H₁(Σ_N; F_ℓ) | the integer has ≍ N·m(Δ_K) digits; mod-M reductions | AGG with fixed coefficients (Prop CYC); EXP as integers |
| 26 | Primitive part Σ_{d\|N} μ(N/d)\|H₁(X_d)_tors\| for fibred K (= primitive periodic points of the monodromy) | Möbius over divisors | — | FACT (Prop CAT) |
| 27 | Torus knot T(p, q), Brieskorn spheres Σ(p, q, r) | Δ_{T(p,q)} = Φ_{pq} when p, q prime; genus φ(N)/2 | need p, q to build | (i) no |
| 28 | Torus knot T(2, N), 2-bridge knots b(N, q) | Δ_{T(2,N)} = (t^N + 1)/(t + 1) factors over Q as ∏_{d\|N, d>1}Φ_{2d}; b(N,q) is Fox m-colourable iff gcd(m, N) > 1 | polynomial factorisation of degree N − 1; gcd | R / EXP |
| 29 | Mazur–Morishita: Spec Z[1/N] as the complement of the "link" of N | components ↔ primes; H¹_et(Spec Z[1/N], μ₂) = (Z/2)^{1+ω}; basis of Z[1/N]^× = the primes | Kummer theory | dimension ω; basis FACT (Prop KUM) |
| 30 | Linking numbers of the two "knots" of N = pq (Legendre (p/q)); Milnor triple invariants (Rédei symbols) | 4-rank of Cl⁺(Q(√N)) (Rédei) | needs the class group | OPEN (§8) |
| 31 | Hyperbolic Dehn fillings M(N, 1), mapping tori of φ^N: volume | vol(M) − c/N² + … analytic in 1/N (Neumann–Zagier 1985); N·vol(M_φ) | — | AGG (no arithmetic) |
| 32 | Kashaev / colored Jones at e(1/N) | — | O(N) recurrences | EXP (`newfamily` row 34) |

### 2.4 Symmetry, combinatorics, dynamics (§4)

| # | Construction | (ii) | (iii) | (iv) |
|---|---|---|---|---|
| 33 | Regular N-gon, star polygons {N/k} | {N/k} has gcd(k, N) components | Euclid per k; finding a disconnected one = finding a factor | R |
| 34 | Total number of components of all {N/k}, k = 1..N (Pillai) | Σ_k gcd(k, N) = (2p − 1)(2q − 1) | — | FACT |
| 35 | Rotations of order exactly N (generators of D_N) | φ(N) | — | FACT (Miller 1976) |
| 36 | Necklaces of length N in k colours (Burnside) | N·Neck = Σ_{g∈Z/N} k^{gcd(g,N)}; at k = 1 + R mod R²: N + R·Pillai(N) | — | FACT, one query (Prop NECK) |
| 37 | Primitive periodic points of a toral automorphism (cat map) at period N | Möbius sum; trace 2 + R mod R² gives R·J₂(N) = R(p² − 1)(q² − 1) | — | FACT, one query (Prop CAT); fixed trace 3 (figure-eight monodromy): small-modulus decoder 40/40 (exploratory) |
| 38 | Periodic points of period N (all, not primitive) | Fix(A^N) = V_N(t) − 2 | powering, poly(log N) | AGG mod M (class G exponent gate, `leakmap` NE1) |
| 39 | Automorphism groups / subgroup lattices of size-N objects | divisor lattice of N | — | R / ω / FACT |

### 2.5 Algebraic and arithmetic geometry (§6)

| # | Construction | (ii) | (iii) | (iv) |
|---|---|---|---|---|
| 40 | Spec Z/N: a finite set of ω(N) points | connected components ↔ idempotents ↔ factors | — | R (definition of factoring) |
| 41 | Point counts #X(Z/N) of a scheme with polynomial count F(r) | F(p)F(q) is a nonconstant polynomial in s unless F is a monomial | count over Z/N | FACT (Prop MF) |
| 42 | #E(Z/N) for one fixed curve | — | — | FACT (`ecount`, `tracedata`) |
| 43 | Frobenius operators on étale Z/N-algebras | — | — | FACT (`nonsmooth` Prop FROB) |
| 44 | Normalisation / ramification locus of Spec Z[x]/(x² − d) | primes dividing the discriminant | ⇔ squarefree kernel (Buchmann–Lenstra 1994; Chistov 1989) | FACT-like (rad) |
| 45 | NFS: geometry of numbers in number fields; toric supports; tropicalisation | — | L[1/3] | covered (`resbarrier`, `heurnfs` §5, `smoothbuild` W3) |
| 46 | Lattices / Coppersmith / Lehman | — | N^{1/4}, N^{1/3} | covered (`coverage38`, `resbarrier`) |
| 47 | Class groups of discriminant −4N (geometry of forms) | ambiguous forms | — | covered (`genusleak`) |

### 2.6 Tilings, quasicrystals, other (§7)

| # | Construction | (iv) |
|---|---|---|
| 48 | Sturmian words / cut-and-project sets of slope a/N or √N | E or AGG (continued fractions of the slope) |
| 49 | Tiling decision problems encoding divisibility (rectangles, Wang tiles) | R; the reduction runs factoring ≤ tiling, giving no algorithm |
| 50 | Deninger/Connes "primes as closed orbits" | no constructible object; counting functions reproduce the explicit formula (sieving) |

(Codex's additional candidates are listed in §7.)

---

## 3. Modular curves: the worked example

### 3.1 The orbifold signature of Γ_0(N)\H (KNOWN formulas; the reading is ours)

For any N (Shimura, *Introduction to the Arithmetic Theory of Automorphic Functions* (1971) §1.6; Diamond–Shurman,
*A First Course in Modular Forms* (2005) §3.9):
- μ = ψ(N) = |P¹(Z/N)|, the index, the degree of X_0(N) → X(1), and 3/π times the hyperbolic area;
- ν₂ = ∏_{r|N}(1 + (−1/r)) (0 if 4 | N), ν₃ = ∏_{r|N}(1 + (−3/r)) (0 if 9 | N);
- ν∞ = Σ_{d|N} φ(gcd(d, N/d));
- g = 1 + μ/12 − ν₂/4 − ν₃/3 − ν∞/2.

For N = pq (odd), each part of the signature lands in a different tier:

| invariant | value for N = pq | tier |
|---|---|---|
| area, index, degree | (p + 1)(q + 1) = N + s + 1 | FACT |
| genus | ((p + 1)(q + 1) − 3ν₂ − 4ν₃)/12 − 1 | FACT |
| cone points of order 2 | 4 if p ≡ q ≡ 1 (mod 4), else 0 | QR (ε) |
| cone points of order 3 | 4 if p ≡ q ≡ 1 (mod 3), else 0 | QR |
| cusps | 4 | ω |

**Prop X0-1 (PROVED, elementary).** From N = pq and g = g(X_0(N)) one recovers p, q in polynomial time.
*Proof.* ν∞ = 4 and ν₂, ν₃ ∈ {0, 4}, so ψ = 12(g + 1) + 3ν₂ + 4ν₃ takes one of four values. For each, set
s = ψ − N − 1 and test whether s² − 4N is a square with root of the right parity; the correct ψ passes and gives
p, q = (s ∓ √(s² − 4N))/2. A wrong candidate either fails or gives the same s. ∎ (E2: 400/400.)

*General N (codex, PROVED from the same formulas).* For N coprime to 11·23, multiplying the level by 11 or 23 kills
both elliptic terms ((−1/11) = (−3/11) = −1, likewise for 23) and doubles the cusps, so
g(23N) − g(11N) = ψ(N); a known ψ(N) splits N in randomized polynomial time because r + 1 | ψ(N) for every r | N
(Lucas-type exponent splitting in (Z/N)[X]/(X² − D)). For X_1: J₂(N) = 2g₁(7N) − 3g₁(5N) + 1 (gcd(N, 35) = 1).
Prior art on factoring from modular-form dimension counts: M. Gu and G. Martin, *Factorization tests and algorithms
arising from counting modular forms and automorphic representations*, Canad. Math. Bull. 62 (2019) 81–97 (codex:
their Thm 10 uses two automorphic-representation counts and a newform count; it is not a one-value theorem).

**Where a geometric computation pays.** Each geometric route to g or μ needs either the formula (the factors) or an
object of size ≥ cN:
- *Triangulations / cell structures.* For a closed triangulated surface, 3F = 2E gives F = 2(V − χ) ≥ −2χ = 4g − 4
  triangles. Since g ≥ ψ/12 − 10/3 (as ν₂, ν₃, ν∞ ≤ 4 for N = pq), any triangulation of X_0(pq) has
  ≥ ψ/3 − 18 > N/3 − 18 faces (PROVED, two lines). Any CW structure has V − E + F = 2 − 2g with V, F ≥ 1, hence
  ≥ 2g ≥ ψ/6 − 7 edges.
- *Presentations of π₁.* Rank H₁ = 2g forces ≥ 2g generators in any presentation of π₁(X_0(N)) (PROVED).
  Compressed descriptions (the implicit permutation action of S, T on P¹(Z/N), which has size O(log N)) do exist.
  But reading g off them by Riemann–Hurwitz needs the cycle counts c(S), c(ST), c(T) *and* μ = |P¹(Z/N)|, and
  counting P¹(Z/N) is computing ψ(N).
- *Farey symbols* (Kulkarni, Amer. J. Math. 113 (1991); Lascurain Orive 2002; Kurth–Long 2007). The image of
  Γ_0(N) in PSL₂(Z) is a free product of 2g + ν∞ − 1 copies of Z, ν₂ of Z/2 and ν₃ of Z/3 (Kulkarni), and a Farey
  symbol yields one independent generator per side pairing, so it has ≥ 2g + ν∞ − 1 paired sides: size Θ(N). The
  algorithms enumerate cosets, so their cost is Θ(μ) = Θ(N) or more.
- *Modular symbols.* Manin symbols are indexed by P¹(Z/N); the space M_2(Γ_0(N)) has dimension 2g + ν∞ − 1 (E1b
  checks this for 188 levels). Linear algebra on it costs poly(ψ(N)).
- *The modular polynomial Φ_N(X, Y)*, a plane model of X_0(N), has degree ψ(N) in each variable.

So the genus is a genuine geometric feature that encodes the factors, and it cannot be computed fast without them
by any known geometric route: the explicit routes are exponential and the compressed ones reduce to counting
P¹(Z/N).

### 3.2 Points over finite fields: the F_ℓ / F_{ℓ²} dichotomy

Let N be squarefree, ℓ a prime not dividing N, and every prime r | N larger than 2√(ℓ^k) when working over F_{ℓ^k}.

**The trace formula (KNOWN).** Eichler–Selberg, in Cohen's form (H. Cohen, *Trace des opérateurs de Hecke sur
Γ_0(N)*, Sém. Théorie des Nombres Bordeaux 1976–77; Cohen–Strömberg, *Modular Forms: A Classical Approach*, AMS GSM
179 (2017), Thm 12.4.11). For k = 2, trivial character and gcd(n, N) = 1:
Tr(T_n | S_2(Γ_0(N))) = A1 + A2 + A3 + A4, where
- A1 = ψ(N)/12 if n is a square, else 0;
- A2 = −½ Σ_{t² < 4n} Σ_f h_w((t² − 4n)/f²) μ(t, f, n), with μ = ψ(N)/ψ(N/N_f) · #{x mod N : x² − tx + n ≡ 0
  (mod N·N_f)}, N_f = gcd(N, f);
- A3 = −½ Σ_{d|n} min(d, n/d) · Σ_{τ|N} φ(gcd(τ, N/τ));
- A4 = σ(n).

E1/E4 verified this implementation against: level 1 (trace 0, 50 values); nine Cremona curves at levels 11–33
including old forms (399 values of a_ℓ, plus a_{ℓ²} and brute-force #E(F_{ℓ²})); and an independent Manin-symbol
computation of T_ℓ and T_{ℓ²} for every squarefree N with ψ(N) ≤ 400 (1501 + 387 traces). 0 mismatches.

**Prop X0-2 (PROVED from the trace formula; formula for squarefree N, r > 2√ℓ for r | N).**

  #X_0(N)(F_ℓ) = 2^{ω(N)} + ½ Σ_{t² < 4ℓ} Σ_f h_w((t² − 4ℓ)/f²) ∏_{r|N} (1 + ((t² − 4ℓ)/r)).

*Proof.* Eichler–Shimura gives #X_0(N)(F_ℓ) = ℓ + 1 − Tr T_ℓ. Here A1 = 0 (ℓ is not a square), A3 = −2^{ω}
(Σ_{d|ℓ} min = 2, and gcd(τ, N/τ) = 1), A4 = ℓ + 1. Since every r | N exceeds 2√ℓ ≥ f, N_f = 1 and
μ = #{x mod N : x² − tx + ℓ ≡ 0} = ∏_r (1 + ((t² − 4ℓ)/r)) by CRT. ∎
- *Moduli reading:* 2^ω rational cusps, plus, for each isogeny class over F_ℓ with trace t, the number of
  Frobenius-stable cyclic subgroups of order N, i.e. the number of roots of x² − tx + ℓ mod N.
- *For N = pq:* ∏ = 1 + χ(N) + χ(p)(1 + χ(N)) with χ = ((t² − 4ℓ)/·). Given N it is 0 when χ(N) = −1 and
  2(1 + χ(p)) otherwise. So the count is an explicit integer combination of the indicators [t² − 4ℓ is a square mod
  N]. Computing it is no harder than quadratic residuosity mod N at these small arguments. Conversely, E3 decoded
  ε, and every unknown symbol, from the counts at ℓ ≤ 13 in 4000/4000 instances (experimentally supported; no general
  uniqueness proof).
- *Verdict:* QR. The point counts of X_0(N) over prime fields see Legendre symbols of p at discriminants of size
  ≤ 4ℓ and nothing else.

**Prop X0-3 (PROVED from the trace formula; E4 800/800).** For ℓ prime, N squarefree and odd with all r | N > 2ℓ:

  #X_0(N)(F_{ℓ²}) = ℓ² + 1 − Tr T_{ℓ²} + ℓ·g(N) = (ℓ − 1)ψ(N)/12 + (bounded QR-type terms).

Precisely, for N = pq: #X_0(N)(F_{ℓ²}) = (ℓ − 1)ψ/12 + ℓ² + 1 − σ(ℓ²) + 2(ℓ + 2) + ½ Σ_{t² < 4ℓ²} Σ_f h_w(·)
∏_r(1 + χ_t(r)) − ℓ(1 + ν₂/4 + ν₃/3). Hence one count over F_{ℓ²} determines ψ(N) up to the 2^{≤ 6} sign patterns
of the symbols involved (at most 61 patterns in E4), and so factors N.

*Proof.* Let α, ᾱ be the Frobenius eigenvalues on H¹, paired with a Hecke eigenbasis of S_2(Γ_0(N)) (T_ℓ is
semisimple there for ℓ ∤ N). Then Σ α² = Σ_f (a_f(ℓ)² − 2ℓ) = Tr(T_ℓ²) − 2ℓg, and T_ℓ² = T_{ℓ²} + ℓ for ℓ ∤ N.
So #X_0(N)(F_{ℓ²}) = ℓ² + 1 − Tr T_{ℓ²} + ℓg. In Tr T_{ℓ²}: A1 = ψ/12; A3 = −½(ℓ + 2)2^ω; A4 = σ(ℓ²);
A2 is QR-type as in X0-2 (the term t = ±2ℓ, i.e. D = 0, is excluded from A2 and is exactly A1). Substitute
g = ψ/12 + 1 − ν₂/4 − ν₃/3 − ν∞/2. The coefficient of ψ is ℓ/12 − 1/12 = (ℓ − 1)/12 > 0. ∎

*Geometric reading.* The term A1 = ψ/12 is the supersingular locus. Over F_{ℓ²} (with the model in which every
supersingular curve has Frobenius ±ℓ) Frobenius is a scalar on E[N], all ψ(N) cyclic subgroups are rational, and the
whole fibre over each supersingular j counts. By the Eichler mass formula, Σ 1/|Aut E| = (ℓ − 1)/24 over
supersingular E, so these fibres contribute ≈ (ℓ − 1)ψ/12 points. Over F_ℓ no Frobenius is scalar.

*Odd extensions (PROVED).* Over F_{ℓ^k} with k odd, Σα^k is an integer combination of Tr T_{ℓ^j} for odd j
(Newton's identities plus T_ℓ T_{ℓ^j} = T_{ℓ^{j+1}} + ℓT_{ℓ^{j−1}}). None of these n = ℓ^j is a square, so A1 never
appears, the g term never appears, and #X_0(N)(F_{ℓ^k}) is QR-type. For even k the ψ term appears. The general
even-k coefficient is (ℓ − 1)/12 by the same supersingular count (sketch; only k = 2 is checked).

**Prop X1 (PROVED modulo the standard description of the cusps; not machine-checked).** If N is squarefree,
ℓ ∤ N, every prime factor of N exceeds ℓ + 1, and N > (√ℓ + 1)², then #X_1(N)(F_ℓ) = φ(N)/2.
*Proof.* For N ≥ 4, X_1(N) is a fine moduli space over Z[1/N]. A non-cuspidal F_ℓ-point is a pair (E/F_ℓ, P) with
P ∈ E(F_ℓ) of order N, impossible since #E(F_ℓ) ≤ (√ℓ + 1)² < N (Hasse). Cusps: on the Tate curve G_m/q^Z the
points of order N are ζ^a q^{b/N}. Those with gcd(b, N) = 1 give φ(N)/2 cusps (up to ±) defined over Q. Those with
gcd(b, N) = d > 1 are defined over a field containing Q(ζ_m)^+ with m ∈ {p, q, N}. For ℓ < m − 1, ℓ ≢ ±1 mod m,
so they are not F_ℓ-rational. (Diamond–Shurman §9.3 and Ogg's description of the rational cusps give the fields of
definition; which half is rational depends on the model, but the count φ(N)/2 does not.) ∎ So #X_1(N)(F_2) = φ(N)/2 for
N = pq with p, q ≥ 5, and Miller's φ-reduction factors N. The same moduli argument is classical in torsion bounds
("Y_1(N)(F_ℓ) = ∅ for N > (√ℓ + 1)²").

**Fricke involution (KNOWN, Ogg 1974, Fricke).** The fixed points of w_N on X_0(N) (N > 4) number
h(−4N) + h(−N)·[N ≡ 3 mod 4] up to the classical corrections. This is a Lefschetz number; it is a class number,
hence FACT (randomized, as row 5). Exact form UNVERIFIED in detail.

### 3.3 Summary of where X_0(N) needs the factors

- Area / genus / degree: either the formula or an object of size Θ(N).
- Cone points / F_ℓ points / return counts of random walks: Legendre symbols of small discriminants (QR).
- F_{ℓ²} points / X_1 points / Fricke fixed points: factoring-equivalent, so no geometric shortcut can exist unless
  factoring itself is easy.
- Cycle structure of a hyperbolic element acting on P¹(Z/N) (the periods): orders of eigenvalues mod p, i.e. order
  finding, FACT (Shor-type).

---

## 4. Counting, symmetry and dynamics

### 4.1 Point counts over Z/N

**Prop MF (PROVED, elementary; the φ and σ cases are KNOWN: Miller 1976, Bach–Miller–Shallit 1986).** Let X be
a scheme of finite type over Z with #X(F_r) = F(r) for all primes r ∤ c, F ∈ Z[x]. Write F = Σ_{i∈I} a_i x^i
with a_i ≠ 0. If |I| ≥ 2, then from N = pq (p, q ∤ c) and #X(Z/N) one recovers p, q in time poly(log N, deg F, log
max|a_i|).
*Proof.* X(Z/N) = X(F_p) × X(F_q) (CRT on points of a finite-type scheme), so the count is
F(p)F(q) = Σ_i a_i² N^i + Σ_{i<j} a_i a_j N^i (p^{j−i} + q^{j−i}). Each power sum p^k + q^k is a monic polynomial of
degree k in s with coefficients in Z[N]. So G(s) := F(p)F(q) is a polynomial in s of degree max I − min I ≥ 1
with leading coefficient a_{max I}a_{min I}N^{min I} ≠ 0. Find the integer roots of G(s) − #X(Z/N) (polynomial time) and test
each as in Prop X0-1. ∎
- Examples: P¹ (ψ), G_m (φ), A¹ ∖ {0, 1} ((p − 2)(q − 2)), SL₂ (N(p² − 1)(q² − 1)), tori, Grassmannians.
- **Prop MF′ (PROVED).** The conic x² + y² = 1 has (r − χ₄(r)) points over F_r, so
  #(Z/N) = N − χ₄(q)p − χ₄(p)q + χ₄(N). If χ₄(p) = χ₄(q) this is N ∓ s + 1; otherwise it is N ± (q − p) − 1.
  Either way one of four candidates for s (or for q − p, with q + p = √((q − p)² + 4N)) works.
- **Contrast.** The *integer* points on the circle x² + y² = N give r₂(N), one QR bit. Points over Z/N on the
  *unit* circle give the factors.

### 4.2 Burnside counts and periodic orbits

**Prop NECK (PROVED; one query).** Let O(k, M) = Σ_{g∈Z/N} k^{gcd(g,N)} mod M. This is N times the number of
k-colour necklaces of length N (Burnside), and it equals Σ_{d|N} φ(d)k^{N/d}. Take R > 4N, k = 1 + R, M = R².
Then O ≡ N + R·Σ_g gcd(g, N) (mod R²), so O gives Pillai's function P(N) = Σ_g gcd(g, N) = (2p − 1)(2q − 1) exactly,
hence s and the factors.
*Proof.* (1 + R)^m ≡ 1 + mR (mod R²). ∎ Geometrically, P(N) is the total number of cycles of all rotations of the
N-gon, i.e. the total number of components of all star polygons {N/g}. So the Burnside count "linearises" to the
star-polygon census.
- Small moduli also suffice (E2, oracle K): from O(k, r) at primes r ≤ 60 and k ∈ {2, 3, 5, 7}, the decoder
  recovers s mod r(r − 1) and CRTs. 400/400 splits, ≤ 24 queries (experimentally supported; the query bound is
  heuristic).

**Prop CAT (PROVED; one query).** Let A ∈ SL₂(Z) be hyperbolic with trace t ≥ 3, acting on the torus. The number of
points of exact period N is Prim_t(N) = Σ_{d|N} μ(N/d)|det(A^d − I)| = Σ_{d|N} μ(N/d)(V_d(t) − 2), where
V_d(t) = λ^d + λ^{−d}. Take t = 2 + R with R > N² and M = R². Since V_d(2 + R) ≡ 2 + d²R (mod R²),
Prim_{2+R}(N) ≡ R·J₂(N) (mod R²), with J₂(pq) = (p² − 1)(q² − 1). This gives p² + q², hence s.
*Proof.* V_d(t) = 2T_d(t/2) with T_d the Chebyshev polynomial, and T_d'(1) = d². So V_d(2 + c) ∈ Z[c] has
constant term 2 and linear coefficient d², i.e. V_d(2 + c) ≡ 2 + d²c (mod c²) in Z[c]. Put c = R, sum with μ(N/d),
and read J₂(N) = Prim/R mod R exactly, as R > J₂(N). ∎
- For the *fixed* figure-eight monodromy (t = 3) the small-modulus decoder split 40/40 at n = 64 (exploratory,
  not pre-registered). With t ∈ {3, 4, 5, 6} it split 400/400 (E2, oracle C).
- Topological reading: for a fibred knot with monodromy A on H₁(fibre) = Z², |H₁(X_d)_tors| = |det(A^d − I)|,
  where X_d is the d-fold cyclic cover of the complement. So the Möbius primitive part of the cover homology
  orders is factoring-equivalent, while each single order |H₁(X_N)_tors| mod M is an aggregate (row 38).

**Moral of §4.** Global counts over Z/N or over the symmetry group Z/N (orbits, Burnside, primitive orbits,
cycles) are divisor sums. Divisor sums of non-monomial weights are factoring-equivalent. The only cheap counts are
those that do not sum over divisors: the full fixed-point count Fix(A^N), a single gcd, a Jacobi symbol.

---

## 5. Topology and arithmetic topology

**Prop COV (PROVED, trivial).** π₁(L(N, q)) = Z/N. Connected d-fold covers correspond to index-d subgroups, which
exist iff d | N; the cover is L(N/d, q). Hence N is prime iff L(N, 1) has no connected cover other than itself and
S³, and the least degree of a nontrivial non-universal cover is the least prime factor of N. This is a restatement
(row 22).

**Prop AGG-H (PROVED).** Let Π_N be a group given by a presentation of fixed shape in which N enters only as an
exponent (relators x^N, or words whose exponents are SLP-computed from N). Examples: L(N, q) (π₁ = Z/N), Dehn
fillings M(N, 1) of a fixed knot (relator μ^N λ), and mapping tori of φ^N. For every fixed finite group G,
#Hom(Π_N, G), and hence the number of connected covers of degree ≤ k, Fox colourings, and Dijkgraaf–Witten
invariants with gauge group G, depend only on N mod exp(G) (mod lcm(1, …, k) for covers).
*Proof.* g^N = g^{N mod exp G} in G, so the set of solutions of the relators in G^{gens} depends only on N mod exp
G. ∎ To see a prime p | N one needs targets of size ≥ p, which the description does not supply. Targets over Z/N
(e.g. SL₂(Z/N)) put the computation in the generic-ring model of Theorem PT (`factorbits` §3.2), where every bit
learned without splitting is paid for by a deviation (`leakmap` Theorem PT-X).

**Prop CYC (PROVED, standard).** For a fixed knot K and a fixed prime ℓ, dim H₁(Σ_N(K); F_ℓ) depends only on N mod
ℓ·m, where m is the order of t on the finite module A_K ⊗ F_ℓ (A_K the Alexander module).
*Proof.* H₁(Σ_N) ≅ A_K/ν_N A_K with ν_N = 1 + t + … + t^{N−1} (Fox; for knots). Tensoring with F_ℓ (right-exact;
H₁(Σ_N; F_ℓ) = H₁(Σ_N) ⊗ F_ℓ by universal coefficients), H₁(Σ_N; F_ℓ) ≅ (A_K ⊗ F_ℓ)/ν_N. On the finite module
A_K ⊗ F_ℓ, t has finite order m. Then ν_{N+ℓm} = ν_N + t^N ν_{ℓm} and ν_{ℓm} = ν_m·(1 + t^m + … + t^{(ℓ−1)m}) = ℓν_m = 0
there. ∎ The integer |H₁(Σ_N)| = ∏_{d|N, d>1}|Res(Δ_K, Φ_d)| does "see" the
divisors of N, but it has ≍ N·log M(Δ_K) digits (M the Mahler measure; Silver–Williams 2002), and its
divisor-indexed factors are exactly the Möbius parts of Prop CAT.

**Prop KUM (PROVED, Kummer theory).** For N odd and squarefree, H¹_et(Spec Z[1/N], μ₂) ≅ Z[1/N]^×/(Z[1/N]^×)² ⊕
Pic(Z[1/N])[2] = (Z/2)^{1+ω(N)}, with basis {−1} ∪ {primes r | N}. In the Mazur–Morishita dictionary (Morishita,
*Knots and Primes*, Universitext 2012, Ch. 2–4) this is H¹ of the complement of the "link" of N. Its dimension is
ω(N) (ω tier) and a basis *is* the factorisation.

**Linking numbers and Rédei symbols (KNOWN; Morishita Ch. 4, 8).** For the two "knots" p, q of N the linking
number is the Legendre symbol (p/q) (as lk mod 2), and the triple Milnor invariants are Rédei symbols. For
p ≡ q ≡ 1 (mod 4), (p/q) = (q/p) controls the 4-rank of the narrow class group of Q(√pq) (Rédei 1934). So "the
linking number of the primes of N" is a class-group 4-rank datum: OPEN from N alone (§8).

**Hyperbolic volumes.** Volumes of manifolds built from N by Dehn filling are real-analytic in the filling slope
(Neumann–Zagier, Topology 24 (1985)). Volumes of mapping tori of φ^N are N·vol(M_φ). Neither carries arithmetic.
Volumes that *do* carry the factors are covolumes of congruence subgroups (rows 12, 21), and those are indices.

---

## 6. Algebraic geometry over Z/N

- Spec Z/N is ω(N) points. Its connected components are the idempotents of Z/N; one nontrivial idempotent e gives
  gcd(e, N). Here geometry *is* the definition of factoring (row 40).
- Point counts: Prop MF (polynomial counts), ecount/tracedata (elliptic curves), FROB (Frobenius operators). All
  FACT.
- Normalisation and ramification: the ramified primes of Q(√d) are the primes of the discriminant; computing a
  maximal order is equivalent to computing the squarefree kernel of the discriminant (Chistov 1989; Buchmann–Lenstra,
  JTNB 6 (1994), Thm 1.3, already cited in `polyimpl`).
- NFS / toric / tropical: covered by `resbarrier`, `heurnfs` §5, `smoothbuild` W3; not redone.

---

## 7. Wild ideas (ours and GPT-6.1-sol's)

Ours (each with a quick kill or keep):
- **Persistent homology of the modular hyperbola {(x, x⁻¹ mod N)}.** The point cloud has N points; any Rips complex
  is EXP. Kill.
- **Weyl law on X_0(N).** The leading eigenvalue-counting term is area/4π, i.e. ψ(N). It is FACT, and spectra need
  discretisations of size ψ(N). Kill.
- **Random walks on the coset graph of P¹(Z/N).** The t-step return count is Σ_{|W|=t} #fix(W). For a hyperbolic
  word W the fixed points are the roots of c x² + (d − a)x − b mod N, ∏_r(1 + ((tr²W − 4)/r)): QR-type. A walk
  of poly(log N) steps sees only these. Kill (QR tier).
- **Deninger/Connes dynamical systems with closed orbits of length log p.** No constructible object; the counting
  functions reproduce the explicit formula, i.e. sieving. Kill.
- **Tilings with divisibility constraints.** They encode factoring (so they are NP-type), but the reduction goes the
  wrong way. Kill.

GPT-6.1-sol's blind list (`codex/geometry_1.out`, 30 candidates) agrees with ours on every shared item. Items it
added or sharpened, now in §2:
- the explicit F_4 count #X_0(N)(F_4) = 2^r + (ψ(N) + 3R_{−4} + 8R_{−3})/12 + R_{−7} + 2R_{−15}, gcd(N, 210) = 1,
  R_D = ∏_{r|N}(1 + (D/r)), derived from the moduli problem (the supersingular curve over F_4 with Frobenius
  [−2], Aut/±1 = A_4, Burnside) and found blind. It agrees with our trace-formula computation on 50 semiprimes and
  4 composite levels (0 mismatches). Codex flagged it as "the most interesting candidate" but did not note that the
  F_ℓ counts are only QR-type;
- the general-N genus reductions g(23N) − g(11N) = ψ(N), J₂(N) = 2g₁(7N) − 3g₁(5N) + 1, and the Gu–Martin
  reference;
- Gaussian Bianchi covers: index I(pq) = ((p+1)(q+1))² (both primes split in Z[i]) or (N − 1)² + s² (both inert);
- Brieskorn spheres S³_{1/N}(trefoil) = Σ(2, 3, 6N − 1) with Casson invariant N; lens-space Heegaard Floer rank N
  and d-invariants: AGG;
- Spec(Z/N): b₀ = ω(N), Pic(P¹_{Z/N}) = Z^{ω(N)}; codex's C2 asks for a compressed chain complex with a Betti
  number equal to ω(N) (OPEN; quick kill: every known construction builds the idempotents);
- a caution that a class number alone is not known to factor every N (applied to rows 5, 11, 19);
- ECM and ECPP as the genuine geometric algorithms (now in §0).

---

## 8. OPEN items, stated precisely

**OPEN-1 (X_0 point counts = small-discriminant QR).** Let L = L(n) → ∞. Given N = pq and the counts
#X_0(N)(F_ℓ) for all primes ℓ ≤ L, find p in time subexponential in n. By Prop X0-2 the input is equivalent to the
Legendre symbols (D/p) for the fundamental D with |D| ≤ 4L and (D/N) = 1 (decoding is unique in all 4000 E3
instances at L = 13). Information-theoretically, L ≈ n² log n should pin p down. No algorithm is known that uses
quadratic-residue patterns of a hidden divisor (Lenstra's divisors-in-residue-classes needs a modulus > N^{1/4} with
a single class, not a half-density pattern). This is the project's QR/genus question (`genusleak` Prop GEN) in
geometric dress, not a new route.

**OPEN-2 (linking numbers of N).** Given N = pq with p ≡ q ≡ 1 (mod 4), decide (p/q) (the mod-2 linking number of
the two "primes" in the arithmetic-topology dictionary, equivalently the 4-rank of Cl⁺(Q(√N))) in polynomial time,
or reduce factoring or QR_N to it. Neither direction is known to us.

Neither OPEN item is a route to faster factoring: even a positive answer to OPEN-2 is one bit, and OPEN-1 asks for a
new algorithm for a problem no geometric construction makes easier.

---

## 9. Experiments (pre-registered; all predictions met)

| Exp | What | Result |
|---|---|---|
| E1a | trace formula vs level 1 (46 l), six genus-1 curves (267 a_ℓ), old-form levels 22, 26, 33 (132) | 0 mismatches |
| E1b | trace formula vs independent Manin-symbol traces, 1501 (N, ℓ) pairs, ψ(N) ≤ 400; dim = 2g + 2^ω − 1 | 0 mismatches (31 s) |
| E1c | orbifold signature from explicit S, ST, T on P¹(Z/N), 1448 levels (all N ≤ 600, squarefree ≤ 2000) | 0 mismatches |
| E2 | decoders G, A, P, L, K, C, 200 instances each at n = 32, 64 | 200/200 for all; K ≤ 24 queries, C ≤ 52 |
| E3 | ε and the symbol pattern from #X_0(N)(F_ℓ), ℓ ≤ 13, N ≡ 1 (4), 2000 instances each at n = 32, 64 | ε determined and correct 2000/2000 at both sizes; ≈ 7 unknown symbols |
| E4a (A1) | T_{ℓ²} trace formula: level 1, curves (incl. brute-force #E(F_{ℓ²})), Manin T_{ℓ²}; 491 checks | 0 mismatches |
| E4b (A1) | factor from #X_0(N)(F_4) alone, and from #X_0(N)(F_9) alone; 200 instances each at n = 32, 64 | 200/200 in all four cells; ≤ 61 sign patterns |

Deviations: smoke tests on seeds 0–5 before the frozen runs (disclosed in prereg A1); the t = 3-only cat-map
decoder (40/40) and the one-query linearised NECK/CAT checks (20/20) are exploratory, not pre-registered.

These are correctness checks. All oracles are computed from p and q; nothing here is a factoring algorithm.

---

## 10. Reproduction (project root)

```
.venv/bin/python -m src.geometry.run_experiments e1 e2 e3 e4
.venv/bin/python -m pytest -q tests/test_geometry.py
```
