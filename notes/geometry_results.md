# geometry: geometric constructions from N, and where the factors hide (phase 26)

Agent `geometry`, 2026-10-24. **Version v4 (final: codex blind and adversarial rounds, Opus proof critic and Opus literature critic all applied).**
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
Atkin–Morain). The number field sieve works in number fields and uses lattices in places. Beyond these known methods, no geometric
construction we examined gave us a faster way to factor. For each one we can point to the specific place where it
needs the factors or size-N work, in the representations we considered. That is not a proof that no geometric
algorithm exists.**

The key point first. Primality itself is already easy: it can be decided in polynomial time (AKS 2002). So a
geometric picture of primality only helps with *factoring* if three things hold at once:
1. the geometric object can be built from N alone, without knowing the factors;
2. some feature of it encodes the factors (or even one bit about them);
3. that feature can be computed geometrically (by counting, linear algebra on a triangulation, volume, fixed-point
   formulas, …) faster than factoring.

Almost every construction we examined fails one of the three, in one of the following recurring ways (class numbers, row 5, form a separate tier):
- **Trivial restatement.** "N is prime iff the only integer rectangle of area N is 1 × N", "iff every star polygon
  {N/k} is connected", "iff the lens space L(N, 1) has no connected covers except the 3-sphere". These are
  definitions in disguise; the geometry adds no algorithm.
- **The feature encodes the factors, so computing it is as hard as factoring.** The cleanest example: the
  **modular curve X_0(N)**. It is defined from N alone. Its hyperbolic area is (π/3)·ψ(N), and for N = pq,
  ψ(N) = (p + 1)(q + 1) = N + (p + q) + 1. So the area gives p + q, and p + q with pq = N gives p and q. The genus
  of X_0(N) does the same (up to four candidates). But every known way to compute the area or genus either uses the
  formula (which needs the factors) or builds the surface explicitly. An explicit surface has at least about
  N/3 triangles, because a surface of genus g needs ≥ 4g − 4 triangles and g ≈ N/12. That is exponential in the
  number of digits of N.
- **The feature is computable, but it only sees "local" information about N** — N modulo a fixed number, or
  quadratic-residue symbols. Lens-space invariants, homology of cyclic covers with fixed coefficients, counts of
  covers of bounded degree, quantum invariants at fixed level and polytope counts with fixed shape are of this kind.
- **The feature sees a little more, but only quadratic-residue information.** Again X_0(N) gives the sharpest
  example. Its number of points over a small finite field F_ℓ is a weighted count of square roots of the numbers
  t² − 4ℓ modulo N, so it is computed from Legendre symbols (D/p) of the discriminants t² − 4ℓ (|t| < 2√ℓ) and nothing else (for p, q > 4ℓ). Its *cone points* do the same:
  X_0(N) has 4 cone points of angle π if −1 is a square mod N and none otherwise.

**The sharpest observation (the ingredients are KNOWN; the framing as a factoring statement is ours, and no prior
formulation of that framing was found by us; verified 800/800; the F_4 case was found first and independently by
GPT-6.1-sol from the moduli interpretation).** The closest prior art: Kirschmer–Voight (SIAM J. Comput. 39 (2010),
Prop 5.8) reduce factoring to computing class numbers of definite quaternion orders. The class number of an Eichler
order of level N is the number of supersingular points of X_0(N) mod ℓ (Deuring; Ribet 1990; Gross 1987), and its
main term is ψ(N)(ℓ − 1)/12. The supersingular contribution to #X_0(N)(F_{ℓ²}) is classical (Ihara 1981;
Tsfasman–Vlăduţ–Zink 1982). For X_0(N) the line between
"quadratic symbols only" and "as hard as factoring" falls exactly between F_ℓ and F_{ℓ²}:
- #X_0(N)(F_ℓ) is QR-type: it is a function of 2^{ω(N)} and the symbols ((t² − 4ℓ)/p), ((t² − 4ℓ)/q).
- #X_0(N)(F_{ℓ²}) is factoring-equivalent: it equals (ℓ − 1)ψ(N)/12 plus bounded QR-type terms. A single number,
  the count of points of X_0(N) over the 4-element field, factors N.
- The geometric reason: over F_{ℓ²} the supersingular elliptic curves have Frobenius equal to the scalar ±ℓ. Every
  one of the ψ(N) cyclic subgroups of order N is then rational, so the whole fibre of X_0(N) → X(1) above each
  supersingular point is counted. Over F_ℓ, Frobenius is never scalar (4ℓ is not a square), and only the subgroups
  fixed by a non-scalar Frobenius count, i.e. square roots of t² − 4ℓ mod N.
- For X_1(N) there is no QR tier at all: once N > (√ℓ + 1)² and every prime factor of N exceeds ℓ + 1, X_1(N) has
  exactly φ(N)/2 points over F_ℓ, all of them cusps. So #X_1(N)(F_2) alone factors N.

**The ladder of lattice points on spheres** shows the same three tiers (cross-referencing `analytic`):
- the circle x² + y² = N: r₂(N) ∈ {0, 16} is one QR bit; two essentially different points factor N; finding one
  point is equivalent to finding √−1 mod N (OPEN);
- the 2-sphere x² + y² + z² = N: for N ≡ 1, 2 (mod 4), r₃(N) = 12h(−4N) (Gauss), a class number. A known class number usually leads to
  an ambiguous form and a factor, though one value is not known to factor every N;
- the 3-sphere: r₄(N) = 8σ(N) = 8(N + p + q + 1) (Jacobi) factors N at once.

**Pattern (HEURISTIC; only the first and third bullets are proved, and only in their stated models).** Invariants of a
compressed object built from N fell into the tiers the project already knows:
- invariants computed through a fixed finite target (homomorphisms to a fixed group, homology with F_ℓ
  coefficients, covers of bounded degree) depend only on N mod m (Prop AGG-H, Prop CYC; PROVED);
- fixed-point counts of a single *non-scalar* transformation (Frobenius over a prime field on X_0(N), one hyperbolic
  word acting on P¹(Z/N)) give quadratic symbols. Fixed-point counts of scalar or near-identity transformations
  (Frobenius over F_{ℓ²} on X_0, Frobenius on X_1, the identity word) count a whole size-N set and are FACT;
- global counts of a size-N set (area, Euler characteristic, number of points of P¹(Z/N), Burnside counts,
  primitive periodic orbits) are multiplicative or Möbius sums over the divisors of N. The ones in Props MF, NECK and
  CAT are factoring-equivalent (PROVED). Not every divisor sum is: Σ_{d|N} φ(d) = N;
- cycle structures and periods are orders of elements, which we expect to be factoring-hard (Shor-type), though
  for a single fixed matrix this is not Miller's random-element reduction and is not proved here.

**What remains OPEN** (§8): (i) whether the counts #X_0(N)(F_ℓ), ℓ ≤ L, for a growing L ≈ (log N)^{O(1)}, which are
computed from Legendre symbols of p at small discriminants, let one find p in polynomial time, or faster than
existing factoring algorithms. No method is known; this is the project's QR/genus question in geometric dress.
(ii) the "linking number" of the two primes of N in the arithmetic-topology dictionary, which is the Legendre symbol
(p/q) and governs the 4-rank of the class group of Q(√N); no polynomial-time algorithm is known.

**Bottom line for the user.** Geometry is an excellent *language* for factorisation: divisors are rectangles,
sublattices, covers and fixed points; primes are the "knots" of arithmetic topology; ECM is geometry of elliptic
curves. But in every construction we know, the geometric feature that sees the factors is either
defined through them or, in every representation we know, costs as much as the arithmetic it re-expresses. We
obtained no faster factoring algorithm. The two OPEN items are genuine unresolved questions, but neither is a
promising lead.

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
- **CLASS** a class number of discriminant ≍ N (computable in L[1/2] under GRH; factor extraction from one value heuristic);
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
| 5 | Sphere x² + y² + z² = N | yes | r₃(N) = 12h(−4N) for N ≡ 1, 2 (mod 4), 24h(−N) for N ≡ 3 (mod 8), 0 for N ≡ 7 (mod 8) (N squarefree > 3; Gauss) | h(−4N): L[1/2] under GRH (Hafner–McCurley 1989) | class number: a known h(−4N) yields an ambiguous form with constant probability on typical N, but the reachable order-2 class can be the non-splitting (2, 2, ·) class for some 2-Sylow shapes, so one value is not known to factor every N (codex). CLASS: a separate tier. Extraction of a factor from one value is HEURISTIC, and the reverse direction (h from the factors) is not known to be polynomial either, so this is not an equivalence |
| 6 | 3-sphere in Z⁴ | yes | r₄(N) = 8σ(N) = 8(N + s + 1), N odd (Jacobi) | ≍ N samples via theta extraction | FACT (`analytic` table) |
| 7 | Conic x² + y² ≡ 1 (mod N) (points over Z/N) | yes | # = (p − χ₄(p))(q − χ₄(q)) = N − χ₄(p)s + 1 or N ± (q − p) − 1 | count over Z/N | FACT (§4.1, Prop MF′) |
| 8 | Lattice polygons with poly(log N)-bit vertices (Pick); Barvinok counts in fixed dimension | yes | Pick: A = I + B/2 − 1, B = Σ gcd(edge) | poly time (Barvinok 1994) | AGG / E (`newfamily` rows 9, 49) |
| 9 | Polytopes with N-dependent normals, Θ(N^{1/3}) vertices | yes | can detect divisors | size of the hull | EXP (`newfamily` row 50) |
| 10 | Sublattices of Z² of index N (N-sheeted covers of the torus) | yes | number = σ(N) | — | FACT (§4.1) |
| 11 | Reduced lattices / flat tori of discriminant −4N (CM points on X(1)); cusps of the Bianchi orbifold of Q(√−N) | yes | number = h(−4N) (order class number); the Bianchi cusp count is the maximal-order h(d_K) instead (h(−35) = 2 vs h(−140) = 6); 2^{t−1} genera, t = number of prime discriminants dividing the discriminant | O(√N) enumeration; N^{1/5} deterministic (ERH; Lenstra 1982, Schoof); L[1/2] randomized (Hafner–McCurley 1989; Biasse–Erukulangara 2021) | CLASS (see #5) |

### 2.2 Modular curves and hyperbolic surfaces (§3)

| # | Feature of X_0(N) / X_1(N) | (ii) | (iii) | (iv) |
|---|---|---|---|---|
| 12 | Hyperbolic area (π/3)ψ(N) = index [SL₂(Z):Γ_0(N)] = \|P¹(Z/N)\| = degree of X_0(N) → X(1) = degree of the modular polynomial Φ_N | ψ(pq) = N + s + 1 | any explicit CW structure has ≥ ψ/6 − 7 edges; Φ_N has degree ψ(N) | FACT (Prop X0-1), EXP |
| 13 | Genus g(X_0(N)) = 1 + ψ/12 − ν₂/4 − ν₃/3 − ν∞/2 | 4 candidates for ψ | triangulation ≥ 4g − 4 triangles; Farey symbol (Kulkarni) Θ(g + cusps); modular-symbol space dimension 2g + ν∞ − 1 | FACT, EXP |
| 14 | Cusps ν∞ = Σ_{d\|N} φ(gcd(d, N/d)) = 2^{ω(N)} (squarefree) | ω | — | ω |
| 15 | Cone points ν₂ = ∏(1 + (−1/r)), ν₃ = ∏(1 + (−3/r)) | ν₂ ∈ {0, 4} decides ε on N ≡ 1 (4) | — | QR |
| 16 | Points over F_ℓ: #X_0(N)(F_ℓ) = 2^{ω} + ½Σ_{t²<4ℓ}Σ_f h_w((t²−4ℓ)/f²)∏_{r\|N}(1 + ((t²−4ℓ)/r)) | Legendre symbols at \|D\| ≤ 4ℓ | E3: the counts for ℓ ≤ 13 determine ε in 4000/4000; the full ≈ 7-symbol pattern in 3461/4000 (otherwise 2 patterns) | QR (Prop X0-2) |
| 17 | Points over F_{ℓ²}: (ℓ − 1)ψ/12 + bounded QR terms | ψ | E4: 800/800 split from one count | FACT (Prop X0-3; closest prior art Kirschmer–Voight 2010 Prop 5.8 via supersingular counts) |
| 18 | Points of X_1(N) over F_ℓ, N squarefree, N > (√ℓ+1)², all r \| N > ℓ + 1: φ(N)/2 rational cusps, no other points | φ | — | FACT (Prop X1) |
| 19 | Lefschetz number of the Fricke involution w_N: #fix = class numbers h(−4N) (+ h(−N) if N ≡ 3 mod 4) (Ogg 1974; Fricke) | h | as #5 | as #5 |
| 19a | Genus of X_1(N); dimensions dim S_k(Γ_0(N)) | g₁ = 1 + J₂(N)/24 − ¼Σ_{d\|N}φ(d)φ(N/d) (N ≥ 5); dim S_{k+12} − dim S_k = ψ(N) for even k ≥ 4 (ψ − 1 at k = 2) | Θ(N²) cosets for Γ_1 | FACT (codex; the identity is elementary from the dimension formula, cf. G. Martin, J. Number Theory 112 (2005), UNVERIFIED). Related prior art: Gu–Martin 2019 (newform / automorphic-representation counts) |
| 20 | Spectrum of the Laplacian / Selberg zeta / random walks on Γ_0(N)\H or its coset graph | Weyl law: area; return counts of one word W: fixed points of W on P¹(Z/N) = ∏(1 + ((tr²W − 4)/r)) | discretisation of size ψ(N) | Weyl term FACT+EXP; fixed points of one word W ≢ ±I mod every r \| N: QR; total return counts include words ≡ ±I mod p (or trivial) and, by mixing, tend to 1/ψ(N): FACT+EXP (proof critic); cycle structure of W = orders (expected FACT, not proved) |
| 21 | Bianchi orbifolds Γ_0(N·O_K)\H³ | covolume × ∏_{r\|N}(r² + 1 inert, (r + 1)² split, r(r + 1) ramified) | Swan/Riley domains have size ∝ index | FACT, EXP |

### 2.3 Topology and arithmetic topology (§5)

| # | Construction | (ii) | (iii) | (iv) |
|---|---|---|---|---|
| 22 | Lens space L(N, q) (Heegaard genus 1, O(log N) description) | connected d-fold covers exist iff d \| N; N prime iff the only nontrivial connected cover is S³ | — | R |
| 23 | Invariants of L(N, q): H₁ = Z/N, Reidemeister torsion, Casson–Walker, d-invariants, WRT, spectrum | Dedekind sums, Gauss sums, lattice counts | Euclid recursions | AGG / E (`newfamily` rows 13, 31; Prop DSF: the local part of s(a, N) is FACT) |
| 24 | Groups/manifolds with N in exponents or surgery coefficients (L(N,q), M(N,1) on a fixed knot, mapping tori of φ^N): #Hom(π₁, G) for fixed finite G; covers of degree ≤ k; Fox colourings; untwisted Dijkgraaf–Witten with fixed gauge group | divisors appear only for target size ≥ p | poly(log N) | AGG (Prop AGG-H; twisted DW needs a separate periodicity argument, codex: the C₃ theory differs at L(15,1), L(21,1)) |
| 25 | N-fold cyclic (branched) covers Σ_N(K) of a fixed knot K | when finite (Res ≠ 0), \|H₁\| = ∏_{d\|N, d>1}\|Res(Δ_K, Φ_d)\|; dim H₁(Σ_N; F_ℓ) | the integer has ≍ N·log M(Δ_K) digits if M(Δ_K) > 1 (periodic if M = 1, e.g. trefoil); mod-M reductions | AGG with fixed coefficients (Prop CYC); EXP as integers |
| 26 | Primitive part Σ_{d\|N} μ(N/d)\|H₁(X_d)_tors\| for fibred K (= primitive periodic points of the monodromy) | Möbius over divisors | — | for the figure-eight knot (t = 3): small-modulus decoder 40/40 (exploratory, multi-query; experimentally supported). Prop CAT's one-query trace 2 + R is not a knot monodromy (Δ(1) = 2 − t ≠ ±1) |
| 27 | Torus knot T(p, q), Brieskorn spheres Σ(p, q, r) | Δ_{T(p,q)} = Φ_{pq} when p, q prime; genus φ(N)/2 | need p, q to build | (i) no |
| 28 | Torus knot T(2, N), 2-bridge knots b(N, q) | Δ_{T(2,N)} = (t^N + 1)/(t + 1) factors over Q as ∏_{d\|N, d>1}Φ_{2d}; b(N,q) is Fox m-colourable iff gcd(m, N) > 1 | polynomial factorisation of degree N − 1; gcd | R / EXP |
| 29 | Mazur–Morishita: Spec Z[1/N] as the complement of the "link" of N | components ↔ primes; H¹_et(Spec Z[1/(2N)], μ₂) = (Z/2)^{2+ω}; basis of Z[1/N]^× = the primes | Kummer theory | dimension 2 + ω; basis FACT (Prop KUM) |
| 30 | Linking numbers of the two "knots" of N = pq (Legendre (p/q)) | 4-rank of Cl⁺(Q(√N)) (Rédei–Reichardt); triple Milnor invariants (Rédei triple symbols) belong to three primes / 8-ranks | no polynomial algorithm known | OPEN (§8) |
| 31 | Hyperbolic Dehn fillings M(N, 1), mapping tori of φ^N: volume | vol(M) − c/N² + … analytic in 1/N (Neumann–Zagier 1985); N·vol(M_φ) | — | no arithmetic: smooth (fillings) or linear (N·vol) in N; not AGG in the N mod m sense, but carries no factor information |
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
| 44 | Normalisation / ramification of Spec Z[x]/(x² − d) and general orders | listing the ramified primes is factoring the discriminant; for squarefree d the maximal order is explicit (d mod 4) and does not factor d | general maximal orders ⇔ squarefree kernel (Buchmann–Lenstra 1994, Thm 1.3; Chistov 1989) | R (listing) / rad-equivalent (general maximal orders) |
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
*Proof (p ≥ 5; for p = 3, ν₃ ∈ {0, 2} and two more candidates are needed).* ν∞ = 4 and ν₂, ν₃ ∈ {0, 4}, so
ψ = 12(g + 1) + 3ν₂ + 4ν₃ takes one of four values. For each, set
s = ψ − N − 1 and test whether s² − 4N is a square with root of the right parity; the correct ψ passes and gives
p, q = (s ∓ √(s² − 4N))/2; a wrong candidate fails the test. ∎ (E2: 400/400.)

*General N (codex, PROVED from the same formulas).* For N coprime to 11·23, multiplying the level by 11 or 23 kills
both elliptic terms ((−1/11) = (−3/11) = −1, likewise for 23) and doubles the cusps, so
g(23N) − g(11N) = ψ(N) (proof critic: 0 failures for admissible N < 600). For squarefree N, ψ(N) = σ(N), and
Bach–Miller–Shallit (SIAM J. Comput. 15 (1986)) factor N from σ(N) in randomized polynomial time. For non-squarefree N
we have not written the ψ → factors reduction out (the naive "r + 1 | ψ" Lucas argument fails when also r − 1 | ψ for
every r, e.g. N = 15, ψ = 24). For X_1: J₂(N) = 2g₁(7N) − 3g₁(5N) + 1 (gcd(N, 35) = 1).
Related prior art: M. Gu and G. Martin, *Factorization tests and algorithms arising from counting modular forms and
automorphic representations*, Canad. Math. Bull. 62 (2019) 81–97. Their counts are A(k, N) (automorphic
representations, i.e. newforms of level dividing N) and newform counts, not dim S_k(Γ_0(N)). By Gekeler's theorem
A(k, N) is computable without factoring for squarefree N, so the newform count is what factors there (literature
critic). Closest prior art for "a class-number-type count with main term ψ(N) factors N": Kirschmer–Voight, SIAM J.
Comput. 39 (2010) 1714–1747, Prop 5.8, which also kills the elliptic terms with an auxiliary prime (13), as the
11N/23N trick does here.

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
- *Farey symbols* (Kulkarni, Amer. J. Math. 113 (1991) 1053–1133; Lascurain Orive, Conform. Geom. Dyn. 6 (2002) 33–60; Kurth–Long, Advances in Algebra and Combinatorics, World Sci. 2008, 225–242, arXiv 0710.1835). The image of
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
179 (2017), Thm 12.4.11 (theorem number UNVERIFIED); restated in K. Popa, arXiv:1706.02691, Thm 3). For k = 2, trivial character and gcd(n, N) = 1:
Tr(T_n | S_2(Γ_0(N))) = A1 + A2 + A3 + A4, where
- A1 = ψ(N)/12 if n is a square, else 0;
- A2 = −½ Σ_{t² < 4n} Σ_f h_w((t² − 4n)/f²) μ(t, f, n), with μ = ψ(N)/ψ(N/N_f) · #{x mod N : x² − tx + n ≡ 0
  (mod N·N_f)}, N_f = gcd(N, f);
- A3 = −½ Σ_{d|n} min(d, n/d) · Σ_{τ|N, gcd(τ, N/τ) | (n/d − d)} φ(gcd(τ, N/τ)) (the divisibility condition is vacuous for squarefree N, the only case used here);
- A4 = σ(n).

E1/E4 verified this implementation against: level 1 (trace 0, 50 values); nine Cremona curves at levels 11–33
including old forms (399 values of a_ℓ, plus a_{ℓ²} and brute-force #E(F_{ℓ²})); and an independent Manin-symbol (Merel's Heilbronn matrices: L. Merel, *Universal Fourier expansions of modular forms*, LNM 1585 (1994) 59–94)
computation of T_ℓ and T_{ℓ²} for every squarefree N with ψ(N) ≤ 400 (1501 + 451 traces). 0 mismatches. These runs
used the brute-force root-count branch; the product branch used by the E3/E4 oracles was compared with it by the
proof critic on 6524 (N, n) cases and in `tests/test_geometry.py` (0 mismatches).

**Prop X0-2 (PROVED from the trace formula; formula for squarefree N, r > 2√ℓ for r | N).**

  #X_0(N)(F_ℓ) = 2^{ω(N)} + ½ Σ_{t² < 4ℓ} Σ_f h_w((t² − 4ℓ)/f²) ∏_{r|N} (1 + ((t² − 4ℓ)/r)).

*Proof.* Eichler–Shimura gives #X_0(N)(F_ℓ) = ℓ + 1 − Tr T_ℓ. Here A1 = 0 (ℓ is not a square), A3 = −2^{ω}
(Σ_{d|ℓ} min = 2, and gcd(τ, N/τ) = 1), A4 = ℓ + 1. Since every r | N exceeds 2√ℓ ≥ f, N_f = 1 and
μ = #{x mod N : x² − tx + ℓ ≡ 0} = ∏_r (1 + ((t² − 4ℓ)/r)) by CRT. ∎
- *Moduli reading:* 2^ω rational cusps, plus, for each isogeny class over F_ℓ with trace t, the number of
  Frobenius-stable cyclic subgroups of order N, i.e. the number of roots of x² − tx + ℓ mod N.
- *For N = pq with p, q > 4ℓ* (so that p, q ∤ t² − 4ℓ; e.g. N = 17·101, ℓ = 13, D = −51 violates this):
  ∏ = 1 + χ(N) + χ(p)(1 + χ(N)) with χ = ((t² − 4ℓ)/·). Given N it is 0 when χ(N) = −1 and 2(1 + χ(p)) otherwise.
  (If χ(N) = 0 then gcd(t² − 4ℓ, N) is a factor; the decoders now check this.) So the count is an explicit integer combination of the indicators [t² − 4ℓ is a square mod
  N]. Computing it is no harder than quadratic residuosity mod N at these small arguments. Conversely, E3 decoded
  ε from the counts at ℓ ≤ 13 in 4000/4000 instances. The full pattern of the ≈ 7 unknown symbols was unique in
  1745/2000 (n = 32) and 1716/2000 (n = 64) instances, with exactly 2 consistent patterns otherwise (recorded
  after the proof critic found that v1's "every unknown symbol" was not supported; experimentally supported, no
  general uniqueness proof).
- *Verdict:* QR. The point counts of X_0(N) over prime fields see Legendre symbols of p at discriminants of size
  ≤ 4ℓ and nothing else.

**Prop X0-3 (PROVED from the trace formula; E4 800/800).** For ℓ prime, N squarefree and odd with all r | N > 4ℓ²:

  #X_0(N)(F_{ℓ²}) = ℓ² + 1 − Tr T_{ℓ²} + ℓ·g(N) = (ℓ − 1)ψ(N)/12 + (bounded QR-type terms).

Precisely, for N = pq: #X_0(N)(F_{ℓ²}) = (ℓ − 1)ψ/12 + ℓ² + 1 − σ(ℓ²) + 2(ℓ + 2) + ½ Σ_{t² < 4ℓ²} Σ_f h_w(·)
∏_r(1 + χ_t(r)) − ℓ(1 + ν₂/4 + ν₃/3). Hence one count over F_{ℓ²} determines ψ(N) up to 2^{O(ℓ²)} sign patterns of the symbols involved (≤ 2⁶ for
ℓ ≤ 3; at most 61 patterns tried in E4), and so factors N.

*Proof.* Let α, ᾱ be the Frobenius eigenvalues on H¹, paired with a Hecke eigenbasis of S_2(Γ_0(N)) (T_ℓ is
semisimple there for ℓ ∤ N). Then Σ α² = Σ_f (a_f(ℓ)² − 2ℓ) = Tr(T_ℓ²) − 2ℓg, and T_ℓ² = T_{ℓ²} + ℓ for ℓ ∤ N.
So #X_0(N)(F_{ℓ²}) = ℓ² + 1 − Tr T_{ℓ²} + ℓg. In Tr T_{ℓ²}: A1 = ψ/12; A3 = −½(ℓ + 2)2^ω; A4 = σ(ℓ²);
A2 is QR-type as in X0-2 (the term t = ±2ℓ, i.e. D = 0, is excluded from A2 and is exactly A1). Substitute
g = ψ/12 + 1 − ν₂/4 − ν₃/3 − ν∞/2. The coefficient of ψ is ℓ/12 − 1/12 = (ℓ − 1)/12 > 0. ∎

*Geometric reading (KNOWN mechanism: Ogg, Bull. SMF 102 (1974), §3 Thm 3 uses exactly y² + y = x³ over F_4 with
Frobenius −2, all cyclic subgroups rational, and the ψ(N)/12 contribution; codex).* A1 = ψ/12 is the identity term
of the trace formula, and its ψ-dependence comes from the supersingular points. It is not literally the supersingular
count: at N = 221, ℓ = 2, A1 = 21 while the supersingular contribution over F_4 is 22, the difference being
automorphism-orbit corrections (codex). Over F_{ℓ²} (with the model in which every supersingular curve has
Frobenius ±ℓ) Frobenius is a scalar on E[N], all ψ(N) cyclic subgroups are rational, and the whole fibre over each
supersingular j counts. By the Eichler mass formula (Deuring/Eichler; Gross 1987), Σ 1/|Aut E| = (ℓ − 1)/24 over
supersingular E, so these fibres contribute (ℓ − 1)ψ/12 + O(1) points. Over F_ℓ no Frobenius is scalar.

*Odd extensions (PROVED).* Over F_{ℓ^k} with k odd, Σα^k is an integer combination of Tr T_{ℓ^j} for odd j
(Newton's identities plus T_ℓ T_{ℓ^j} = T_{ℓ^{j+1}} + ℓT_{ℓ^{j−1}}). None of these n = ℓ^j is a square, so A1 never
appears, the g term never appears, and #X_0(N)(F_{ℓ^k}) is QR-type. *Even extensions (PROVED; proof critic).*
#X_0(N)(F_{ℓ^k}) = ℓ^k + 1 − Tr T_{ℓ^k} + ℓ·Tr T_{ℓ^{k−2}}; both traces have A1 = ψ/12, so the ψ coefficient is
(ℓ − 1)/12 for every even k. The critic checked k = 3 (Manin, 160 levels) and k = 4 (#E(F_16)).

**Prop X1 (PROVED modulo the standard description of the cusps; checked by the proof critic with Γ_1(N) Manin
symbols at (N, ℓ) = (15,2), (21,2), (21,5), (35,2), (35,3), (39,2), (39,5), (55,3), (65,3)).** If N is squarefree,
ℓ ∤ N, every prime factor of N exceeds ℓ + 1, and N > (√ℓ + 1)², then #X_1(N)(F_ℓ) = φ(N)/2.
*Proof.* For N ≥ 4, X_1(N) is a fine moduli space over Z[1/N]. A non-cuspidal F_ℓ-point is a pair (E/F_ℓ, P) with
P ∈ E(F_ℓ) of order N, impossible since #E(F_ℓ) ≤ (√ℓ + 1)² < N (Hasse). Cusps: on the Tate curve G_m/q^Z the
points of order N are ζ^a q^{b/N}. Those with gcd(b, N) = 1 give φ(N)/2 cusps (up to ±) defined over Q. Those with
gcd(b, N) = d > 1 are defined over Q(ζ_r) when d = r is a prime factor of N, and over Q(ζ_N)^+ when d = N. Since
every r | N exceeds ℓ + 1, ℓ ≢ 1 mod r and ℓ ≢ ±1 mod N, so they are not F_ℓ-rational. (The hypothesis is needed:
#X_1(15)(F_7) = 8 ≠ 4 and #X_1(33)(F_7) = 20 ≠ 10 (proof critic). For ω = 2 it can be weakened to ℓ ≢ 1 mod each
r | N and ℓ ≢ ±1 mod N.) (Diamond–Shurman §9.3 and Ogg's description of the rational cusps give the fields of
definition; which half is rational depends on the model, but the count φ(N)/2 does not.) ∎ So #X_1(N)(F_2) = φ(N)/2 for
N = pq with p, q ≥ 5 (p = 3 is covered by the weakened ω = 2 condition, 2 ≢ 1 mod 3, if one accepts it), and φ gives
s = N + 1 − φ directly. Ingredients KNOWN: Y_1(N)(F_ℓ) = ∅ by Hasse is standard in torsion bounds (Kamienny, Invent.
Math. 109 (1992); Merel, Invent. Math. 124 (1996)); rational cusps: Ogg, Proc. Symp. Pure Math. 24 (1973). The
count as a factoring statement: no prior formulation found by us. The same moduli argument is classical in torsion bounds
("Y_1(N)(F_ℓ) = ∅ for N > (√ℓ + 1)²").

**Fricke involution (KNOWN, Ogg 1974, Fricke).** The fixed points of w_N on X_0(N) (N > 4) number
h(−4N) + h(−N)·[N ≡ 3 mod 4] (A. Ogg, Hyperelliptic modular curves, Bull. SMF 102 (1974) 449–462, Prop 3). This is a Lefschetz number and a class number:
CLASS tier (as row 5; extraction HEURISTIC). Exact form checked by the literature critic via secondary sources.

### 3.3 Summary of where X_0(N) needs the factors

- Area / genus / degree: either the formula or an object of size Θ(N).
- Cone points / F_ℓ points / fixed points of one non-scalar word: Legendre symbols of small discriminants (QR).
- F_{ℓ²} points / X_1 points: factoring-equivalent, so no geometric shortcut can exist unless factoring itself is
  easy. Fricke fixed points: a class number (CLASS tier). Total random-walk return counts: FACT (they include the
  ψ(N) fixed points of words that are ±I mod p).
- Cycle structure of a hyperbolic element acting on P¹(Z/N) (the periods): orders of eigenvalues mod p, i.e. order
  finding (expected FACT, Shor-type; not proved for a fixed element).

---

## 4. Counting, symmetry and dynamics

### 4.1 Point counts over Z/N

**Prop MF (PROVED, elementary; the φ and σ cases are KNOWN: Miller, JCSS 13 (1976) (deterministic under ERH; the randomized unconditional version is folklore), Bach–Miller–Shallit 1986; cf. Woll, Inform. Comput. 72 (1987)).** Let X be
a scheme of finite type over Z with #X(F_r) = F(r) for all primes r ∤ c, F ∈ Z[x]. Write F = Σ_{i∈I} a_i x^i
with a_i ≠ 0. If |I| ≥ 2, then from N = pq (p, q ∤ c) and #X(Z/N) one recovers p, q in time poly(log N, deg F, log
max|a_i|).
*Proof.* X(Z/N) = X(F_p) × X(F_q) (CRT on points of a finite-type scheme), so the count is
F(p)F(q) = Σ_i a_i² N^i + Σ_{i<j} a_i a_j N^i (p^{j−i} + q^{j−i}). Each power sum p^k + q^k is a monic polynomial of
degree k in s with coefficients in Z[N]. So G(s) := F(p)F(q) is a polynomial in s of degree max I − min I ≥ 1
with leading coefficient a_{max I}a_{min I}N^{min I} ≠ 0. Find the integer roots of G(s) − #X(Z/N) (polynomial time) and test
each as in Prop X0-1. ∎
- Examples: P¹ (ψ), G_m (φ), A¹ ∖ {0, 1} ((p − 2)(q − 2)), SL₂ (N(p² − 1)(q² − 1)), split tori, Grassmannians.
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
- Topological reading (genus-one fibred knots only have t = 3, the figure-eight; trace 2 + R is not a knot monodromy,
  so for knots only the multi-query small-modulus result applies): |H₁(X_d)_tors| = |det(A^d − I)|,
  where X_d is the d-fold cyclic cover of the complement. So the Möbius primitive part of the cover homology
  orders is factoring-equivalent, while each single order |H₁(X_N)_tors| mod M is an aggregate (row 38).

**Moral of §4 (a summary, not a theorem).** Global counts over Z/N or over the symmetry group Z/N (orbits,
Burnside, primitive orbits, cycles) are divisor sums. The specific ones in Props MF, NECK and CAT are
factoring-equivalent. Divisor sums in general need not be: Σ_{d|N} φ(d) = N (codex). The cheap counts we met are
the ones that collapse like that, or that do not sum over divisors at all: the full fixed-point count Fix(A^N), a
single gcd, a Jacobi symbol.

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
invariants with gauge group G, depend only on N mod exp(G) (mod lcm(1, …, k) for covers) when N enters as an
exponent of group elements. For mapping tori of φ^N, the dependence is on N modulo the order of φ* acting on the
finite set Hom(π₁(F), G) (e.g. period 3 for the cat map and G = Z/2): still an aggregate.
*Proof.* g^N = g^{N mod exp G} in G, so the set of solutions of the relators in G^{gens} depends only on N mod exp
G; for mapping tori, Hom(Π_N, G) is determined by the fixed points of (φ*)^N twisted by G, which is periodic in N. ∎ To see a prime p | N one needs targets of size ≥ p, which the description does not supply. Targets over Z/N
(e.g. SL₂(Z/N)) put the computation in the generic-ring model of Theorem PT (`factorbits` §3.2), where every bit
learned without splitting is paid for by a deviation (`leakmap` Theorem PT-X).

**Prop CYC (PROVED, standard).** For a fixed knot K and a fixed prime ℓ, dim H₁(Σ_N(K); F_ℓ) depends only on N mod
ℓ·m, where m is the order of t on the finite module A_K ⊗ F_ℓ (A_K the Alexander module).
*Proof.* H₁(Σ_N) ≅ A_K/ν_N A_K with ν_N = 1 + t + … + t^{N−1} (Fox; for knots). Tensoring with F_ℓ (right-exact;
H₁(Σ_N; F_ℓ) = H₁(Σ_N) ⊗ F_ℓ by universal coefficients), H₁(Σ_N; F_ℓ) ≅ (A_K ⊗ F_ℓ)/ν_N. On the finite module
A_K ⊗ F_ℓ, t has finite order m. Then ν_{N+ℓm} = ν_N + t^N ν_{ℓm} and ν_{ℓm} = ν_m·(1 + t^m + … + t^{(ℓ−1)m}) = ℓν_m = 0
there. ∎ The integer |H₁(Σ_N)| = ∏_{d|N, d>1}|Res(Δ_K, Φ_d)| does "see" the
divisors of N, but it has ≍ N·log M(Δ_K) digits (M the Mahler measure, when M(Δ_K) > 1; for M = 1, e.g. torus knots, the orders are periodic (Gordon); Riley 1990; González-Acuña–Short 1991; Silver–Williams, Topology 41 (2002) 979–991), and its
divisor-indexed factors Res(Δ_K, Φ_d) come from *multiplicative* Möbius inversion of the orders. Prop CAT's primitive
count is the *additive* Möbius sum. These are different numbers: for the figure-eight at N = 15 the additive count is
1,860,360 and the multiplicative cyclotomic factor is 961 (codex).

**Prop KUM (PROVED, Kummer theory).** For N odd and squarefree, put X = Spec Z[1/(2N)] (2 must be inverted for μ₂
to be étale; codex). The Kummer sequence gives H¹_et(X, μ₂) ≅ O(X)^×/(O(X)^×)² ⊕ Pic(X)[2] = (Z/2)^{2+ω(N)}, with
basis {−1, 2} ∪ {primes r | N}. In the Mazur–Morishita dictionary (Morishita,
*Knots and Primes*, Universitext 2012, Ch. 2–4) this is H¹ of the complement of the "link" of N. Its dimension is
2 + ω(N) (ω tier) and a basis *is* the factorisation.

**Linking numbers and Rédei symbols (KNOWN; Morishita Ch. 4, 8).** For the two "knots" p, q of N the linking
number is the Legendre symbol (p/q) (as lk mod 2); triple Milnor invariants are Rédei symbols, relevant for three
primes or the 8-rank. For p ≡ q ≡ 1 (mod 4), (p/q) = (q/p) controls the 4-rank of the narrow class group of Q(√pq)
(Rédei, J. Reine Angew. Math. 171 (1934); Rédei–Reichardt, ibid. 170 (1934); UNVERIFIED). Chapter numbers refer to
Morishita's 1st edition. So "the
linking number of the primes of N" is a class-group 4-rank datum: OPEN from N alone (§8).

**Hyperbolic volumes.** Volumes of long hyperbolic Dehn fillings built from N have an asymptotic expansion in the
filling slope near the complete structure (Neumann–Zagier, Topology 24 (1985). Volumes of mapping tori of φ^N are N·vol(M_φ). Neither carries arithmetic.
Volumes that *do* carry the factors are covolumes of congruence subgroups (rows 12, 21), and those are indices.

---

## 6. Algebraic geometry over Z/N

- Spec Z/N is ω(N) points. Its connected components are the idempotents of Z/N; one nontrivial idempotent e gives
  gcd(e, N). Here geometry *is* the definition of factoring (row 40).
- Point counts: Prop MF (polynomial counts), ecount/tracedata (elliptic curves), FROB (Frobenius operators). All
  FACT.
- Normalisation and ramification: *listing* the ramified primes of Q(√d) is factoring the discriminant. For
  squarefree d the maximal order is explicit and does not factor d. For general orders, computing the maximal order
  is equivalent to computing the largest squarefree divisor of the discriminant (Chistov 1989; Buchmann–Lenstra,
  JTNB 6 (1994), Thm 1.3, already cited in `polyimpl`).
- NFS / toric / tropical: covered by `resbarrier`, `heurnfs` §5, `smoothbuild` W3; not redone.

---

## 7. Wild ideas (ours and GPT-6.1-sol's)

Ours (each with a quick kill or keep):
- **Persistent homology of the modular hyperbola {(x, x⁻¹ mod N)}.** The point cloud has N points; any Rips complex
  is EXP. Kill.
- **Weyl law on X_0(N).** The leading eigenvalue-counting term is area/4π, i.e. ψ(N). It is FACT, and spectra need
  discretisations of size ψ(N). Kill.
- **Random walks on the coset graph of P¹(Z/N).** For one word W ≢ ±I mod every r | N the fixed points are the roots
  of c x² + (d − a)x − b mod N, ∏_r(1 + ((tr²W − 4)/r)): QR-type. But the total t-step return count Σ_{|W|=t}
  #fix(W) includes trivial words (each contributes ψ(N)) and words ≡ ±I mod p only (length O(log p) exists), and by
  expander mixing the return probability tends to 1/ψ(N): FACT, and computing it needs the ψ(N)-vertex graph.
  (v1 labelled this QR; corrected after the proof critic.) Kill (FACT + EXP).
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

**OPEN-1 (X_0 point counts = small-discriminant QR).** Let L = L(n) → ∞ with L ≤ n^{O(1)}. Given N = pq and the counts
#X_0(N)(F_ℓ) for all primes ℓ ≤ L, find p in polynomial time, or in time below the best known factoring bounds
(heuristic L[1/3, 1.923]; rigorous L[1/2, 1 + o(1)], Lenstra–Pomerance 1992). "Subexponential" alone would be
trivial (codex). Proved: by Prop X0-2 the counts are computable from N, 2^ω and the symbols (D/p) for the fundamental
discriminants D of the orders of discriminant t² − 4ℓ (t² < 4ℓ, ℓ ≤ L) with (D/N) = 1. This is not every |D| ≤ 4L:
at L = 13 the set is {−3, −4, −7, −8, −11, −19, −20, −24, −35, −40, −43, −51, −52}, which omits −23, −31, −47. At
L = 13 the counts determined ε in all 4000 E3 instances and the full pattern in 87%. HEURISTIC: L ≈ n² log n should
pin p down information-theoretically. No algorithm is known that uses
quadratic-residue patterns of a hidden divisor (Lenstra's divisors-in-residue-classes, Math. Comp. 42 (1984), needs a modulus > N^{1/3}, and Coppersmith–Howgrave-Graham–Nagaraj, Math. Comp. 77 (2008), > N^{1/4}, with
a single class, not a half-density pattern). This is the project's QR/genus question (`genusleak` Prop GEN) in
geometric dress, not a new route.

**OPEN-2 (linking numbers of N).** Given N = pq with p ≡ q ≡ 1 (mod 4), decide (p/q) (the mod-2 linking number of
the two "primes" in the arithmetic-topology dictionary, equivalently whether the 4-rank of Cl⁺(Q(√N)) is 1) in
polynomial time. Or give a reduction from factoring or from QR_N to an oracle answering this for *arbitrary* such
inputs N′. A reduction that only uses the bit at the same N is trivial: try both values. No polynomial-time algorithm
and no such reduction is known to us.

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
| E3 | ε and the symbol pattern from #X_0(N)(F_ℓ), ℓ ≤ 13, N ≡ 1 (4), 2000 instances each at n = 32, 64 | ε determined and correct 2000/2000 at both sizes; ≈ 7 unknown symbols; pattern unique 1745 / 1716, else 2 patterns (re-run with the pattern count recorded) |
| E4a (A1) | T_{ℓ²} trace formula: level 1, curves (incl. brute-force #E(F_{ℓ²})), Manin T_{ℓ²}; 491 checks | 0 mismatches |
| E4b (A1) | factor from #X_0(N)(F_4) alone, and from #X_0(N)(F_9) alone; 200 instances each at n = 32, 64 | 200/200 in all four cells; ≤ 61 sign patterns |

Deviations (disclosed):
- smoke tests on seeds 0–5 before the frozen runs (prereg A1);
- oracle C was registered with trace 3 only but ran with t ∈ {3, 4, 5, 6} per small prime; oracle K was registered with
  k a primitive root mod r but ran with k ∈ {2, 3, 5, 7};
- E3 was re-run unchanged (same seeds, same results) after v1 to record the number of consistent symbol patterns,
  which the frozen code discarded;
- exploratory, not pre-registered, archived as `run_experiments explore` → `experiments/raw/geometry_explore.jsonl`:
  the t = 3-only cat-map decoder (40/40, ≤ 13 queries, n = 64) and the one-query linearised NECK/CAT checks (20/20 each).

These are correctness checks. All oracles are computed from p and q; nothing here is a factoring algorithm.

---

## 10. Reviews

**Opus proof critic (read-only): ACCEPT WITH FIXES.** It independently re-checked the genus identities (N < 600), the
product vs brute branches of the trace formula (6524 cases), codex's F_4 formula (883 levels), odd and even extensions
(k = 3 by Manin symbols at 160 levels, k = 4 against #E(F_16)) and Prop X1 (Γ_1 Manin symbols at 9 levels).
All 20 findings are applied above. The three MAJOR ones:
1. v1's "the counts determine every unknown symbol" was false: 13–14% of instances have two consistent patterns. Only
   ε is always determined. E3 was re-run with the pattern count recorded.
2. Total random-walk return counts are FACT, not QR.
3. The §0 "general principle" was overclaimed. It is now a HEURISTIC pattern, with the scalar-Frobenius exceptions
   stated.

MINOR fixes: the X_1 hypotheses and cusp fields, edges vs faces, the knot-monodromy caveat for Prop CAT, the
AGG-H modulus for mapping tori, the p, q > 4ℓ condition, p = 3 in X0-1, Fricke consistency, k ≥ 4 for dimension
differences, provable even extensions, the KUM dimension, the check counts, the disclosure of deviations, the r₃
cases and the L[1/2] cost, "every known way", the ramified Bianchi factor, split tori, and σ = ψ with
Bach–Miller–Shallit.

**Opus literature critic (read-only, web access): OK WITH FIXES.** All 16 findings are applied. The main one is
prior art for the F_{ℓ²} statement. Kirschmer–Voight (2010) Prop 5.8 reduces factoring to class numbers of definite
quaternion orders, which by Deuring are supersingular counts on X_0(N) mod ℓ, with main term ψ(N)(ℓ − 1)/12. The
supersingular contribution to #X_0(N)(F_{ℓ²}) is classical (Ihara; Tsfasman–Vlăduţ–Zink). Prop X0-3 is therefore
labelled as a derivation of known ingredients, and only the F_ℓ / F_{ℓ²} factoring framing is claimed, as "no prior
formulation found by us". Other fixes: the A3 divisibility condition; Gu–Martin describes newform and representation
counts, not dim S_k; Merel cited; X_1 ingredients KNOWN; Lenstra's N^{1/3}; Ogg Prop 3; "NFS is geometry of numbers"
softened; r₃ cases; Mahler measure M = 1; Rédei details; full Farey references. The critic found no prior statement
that Pillai or necklace counts are FACT, and no reduction for √−1 mod N (an UNVERIFIED negative).
Shallit–Shamir (IPL 20 (1985)) covers only d, μ, Ω and exponent-multiset functions.

**GPT-6.1-sol adversarial round (`codex/geometry_2.out`, reviewing v1.1): MAJOR REVISION.** It confirmed the central
formulas independently:
- X0-2;
- the F_{ℓ^k} identity #X_0(N)(F_{ℓ^k}) = ℓ^k + 1 − Tr T_{ℓ^k} + ℓ·Tr T_{ℓ^{k−2}} (46 extension-field curve counts for
  k ≤ 5, plus Manin checks);
- X1, counting cusps by their Frobenius action in 426 admissible cases with N < 300, all φ(N)/2;
- NECK, CAT, MF, MF′, COV and CYC.

Of its 21 issues, about half had already been fixed after the proof critic: E3 uniqueness, random walks and scalar
words, X_1 hypotheses, CW faces, KUM wording, p = 3, the r₃ cases, and the general-principle overclaim. Applied in
v4:
- OPEN-1 as first stated ("subexponential") was trivially solved by existing algorithms. It is now restated with a
  polynomial-time or beat-the-record target and the exact discriminant set.
- AGG-H's modulus is fixed for mapping tori, and twisted DW is excluded (C₃ counterexample).
- Prop CAT applies only to the variable-trace family, not to a fixed knot.
- Additive and multiplicative Möbius inversion are now separated.
- H₁ finiteness and Mahler-measure hypotheses are added.
- Class numbers are a separate CLASS tier, not FACT, because neither direction is a proven reduction.
- A1 is the identity term, not literally the supersingular count. Ogg 1974 §3 is cited as prior art for the
  F_4 mechanism.
- Row 11 (order vs maximal-order class number; genus count), the cost column (L[1/2] randomized), and the
  normalisation row are corrected.
- N·vol is "no arithmetic", not AGG.
- OPEN-2 now states its oracle model.
- §0 is softened: "we point to where each construction needs the factors in the representations considered", not
  "we can say exactly why"; no "no live lead" dismissal.

Codex's closing remark is adopted as the reading guide. This note distinguishes four kinds of statement: (a) known
algorithms (ECM, ECPP); (b) oracle reductions (X0-1, X0-3, X1, MF, NECK, CAT); (c) output-size obstructions for
specified representations (§3.1); (d) unresolved problems (OPEN-1, OPEN-2). None of (c) is a lower bound for
factoring.

---

## 11. Reproduction (project root)

```
.venv/bin/python -m src.geometry.run_experiments e1 e2 e3 e4 explore
.venv/bin/python -m pytest -q tests/test_geometry.py
```
