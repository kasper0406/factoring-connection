# The public ellipse count modulo 16

Research note and reproducible benchmarks, 26 November 2026. Base revision `9e70736b4abaacff226c93e4fdddaf5e6fcad5f5` of [factoring-connection](https://github.com/kasper0406/factoring-connection/tree/9e70736b4abaacff226c93e4fdddaf5e6fcad5f5). The original definitions, eligibility conditions and conditional coverage statement are in the [canonical-index note](https://github.com/kasper0406/factoring-connection/blob/9e70736b4abaacff226c93e4fdddaf5e6fcad5f5/notes/shared_collision_nov26.md).

The radial evaluator remains open. A small extension of it has a stronger payoff: counting the same points on a prefix of one public row coordinate suffices to factor the input from a single eligible branch. The finite Fourier transform preserves the mixed lines, including their coefficients modulo 16. On the independent lifting track, inverse Gauss duplication and enumeration of the 2-torsion recovered all the expected roots in the exact tests.

## Normalization and assumptions

Use the repository's principal form `Q(x,y) = x²+bxy+dy²`, with `b = Δ mod 2`, `d=(Δ+b²)/4`, and discriminant `−Δ < −4`. The branch gives odd distinct primes `p<q≤2p`, `N=pq`, a public `ℓ,c>0` such that `4N=c²+Δℓ²`, and `gcd(N,ℓΔ)=1`. The factors are hidden from the algorithms for the public count. Assume `B<min(p,q)` and the repository's largest-lift convention when using its bound `Δ<4√N`.

Write `R(B)` for the count of all signed integer vectors with `0<Q≤NB` and `N|Q`. The repository's `U_N` counts up to simultaneous negation. The raw query is `R mod 16`; `U_N = R/2 mod 8`. If `P(B)` and `V(B)` count nonzero vectors up to sign through B in the principal and mixed scaled forms, then exactly

    R(B) = 4 P(B) + 4 V(B).

All four CRT lines are disjoint on the counted nonzero points. An intersection at one prime would make the positive norm a multiple of at least `pN`. Both principal lines can be counted exactly. The two mixed scaled forms are inverse classes and have the same theta series. At first entry, their minimum is denoted by `a`. It is the minimum of the **squared** class.

## A one-coordinate prefix oracle factors N

Write `m=(c−bℓ)/2`. The vector `α=(m,ℓ)` is primitive and has norm N. Define its integer determinant against a vector by

    L(x,y) = m y − ℓ x.

Consider a query returning the count modulo 8 with the single extra condition `0<L≤t`. Denote the response by `A(B,t)`. Both cutoffs are public. An algorithm for A gives a proper factor using O(log N) queries on a single eligible branch, provided the first doubling endpoint stays below p, as it does for the sufficiently large balanced inputs in the repository. This direction requires neither a factorization of Δ nor a grid hit.

1. Set `T(B)=floor(2N√(B/Δ))`. The Gram determinant bound gives `|L|≤T(B)`. All vectors with `L=0` are integral multiples of α. Negation pairs the two strict half-planes, giving `R(B) = 2 A(B,T(B)) + 2 floor(√B) mod 16`.
2. Recover the first mixed norm a by doubling and binary search, after subtracting `4P(B)`. At a doubling endpoint below `2a` there are at most three mixed pairs in each scaled lattice. In any binary lattice, two different antipodal pairs of norm below twice the minimum cannot occupy the same parity class: their half-sum and half-difference would have total norm below twice the minimum. Only three nonzero parity classes are available. Thus the first nonempty prefix is visible modulo 4 after normalization.
3. The reduced-form inequality gives `a≤√(Δ/3)<Δ/4` for all eligible Δ. Principal vectors at this cutoff are exactly the signed multiples of α and of its conjugate. The first direction has L=0; on the conjugate direction the step in |L| is `cℓ`. At the first mixed shell the exact positive mixed count through t is therefore the least nonnegative residue of

       A(a,t) − min(floor(√a), floor(t/(cℓ)))   modulo 8.

   There are at most six positive mixed vectors across both mixed lattices and at least one. Hence this residue cannot wrap. It is a monotone integer count at this shell.
4. Binary-search the first t with positive mixed count. At a mixed vector, its CRT line agrees with the line of α at exactly one prime. Therefore `gcd(t,N)` is that prime. No mixed root or point coordinates are recovered by the reduction.

The clipped oracle in the benchmark is simulated by slow public enumeration. Only its residues reach the factor-extraction logic. The 194 checks all passed, using at most 25 queries per branch. For N=11009 the first positive mixed coordinate was 3924 with gcd 109. For N=28829 it was 4767 with gcd 227. This identifies precisely what a floor-sum reciprocity that also supports row prefixes would deliver. The radial oracle alone does not supply those prefixes.

## Two exact floor-sum formulations

**Swapping the norm for its quotient.** The unsplit point indicator is the difference of the integer floors of `Q/N` and `(Q−1)/N`. Summing by the positive quotient `j=Q/N` gives

    (2x+b y)² = 4Nj − Δy²,      1≤j≤B.

For each j the implementation visits `0≤y≤floor(2√(Nj/Δ))`. An integer square and a parity check give the signed multiplicity, with both zero-coordinate ties handled. No modular square root is used. The exact number of tests is the sum of those row counts, of order `B + B^(3/2)√(N/Δ)`. Even at B=1 this is Ω(N^(1/4)) under the repository's Δ bound. At `B≈√Δ` it is of order `√N Δ^(1/4)`. The modulus has been replaced with an integer-square constraint.

**Swapping into the public norm-N basis.** Extend α to an integral unimodular basis. In its coordinates the form is `Ns²+hsy+ky²`, with `h²−4Nk=−Δ` and h a unit modulo N. Here y is L. The entire congruence becomes `N | y(hs+ky)`. Put `M=N/gcd(y,N)`. On row y the admissible s occupy the single progression `s=−h⁻¹ky mod M`. Intersect it with the exact ellipse interval by two floors. This gives a correct row formula and mutates the form without using the mixed root. It visits `floor(2N√(B/Δ))` positive row indices; row gcds can already reveal a factor.

More generally, if the first basis vector is v, the ellipse has real half-height `2√(Q(v)NB/Δ)` in the other coordinate. The minimum possible integral Q(v) is 1. A change of unimodular basis alone cannot make a full row scan polylogarithmic. These bounds concern the stated scans; a batched recurrence would need its own analysis.

## Fourier transform and the 2-adic test

For s=p,q, let λ_s be the sum of the indicators of the two isotropic lines over the prime field, so λ_s is 2 at the origin. Define W=λ_p λ_q, the sum of the four CRT-line indicators. In the allowed positive cutoff W is the ordinary congruence indicator. Globally it is a different periodic weight. It can be evaluated pointwise from public data: it is zero outside the congruence, and on it has values 1, 2, 4 according as `gcd(x,y,N)` is 1, a proper divisor, or N. Subtract the indicators of the two public lines to obtain the mixed weight M.

Use the unnormalized finite Fourier transform with kernel `exp(2πi det(u,v)/N)`. On a primitive cyclic line H of order N, character orthogonality gives

    F(1_H)(u) = N 1_H(u).

Indeed the determinant annihilator of H is H itself. Hence `F(M)=N M` and `F(W)=N W`. The identity is integral in the cyclotomic ring and persists on reduction modulo 16. On the primitive points of the two mixed lines, the coefficient is the unit N modulo 16. There are exactly `2φ(N)` such frequency indices in a full period. Listing the complete nonzero spectrum costs Ω(N). Locating even one primitive mixed index permits a proper gcd using L. This does not bound the cost of aggregate evaluation.

For the **literal** global cone indicator, inclusion-exclusion instead gives `F(1_cone)=(pλ_p−1)(qλ_q−1)`. It equals 1 outside both dual cones. The extra intersection terms are irrelevant to the positive sharp cutoff below min(p,q), but cannot simply be dropped from an unrestricted Gaussian sum.

Let `D₀=Δ/4` be the determinant of the Gram matrix of Q. For the completed weight M, including its value at the origin, set `Θ(t)=Σ_v M(v) exp(−πtQ(v)/N)`. The Gaussian Fourier transform, the preceding finite identity, and the 90-degree rotation taking the inverse binary Gram matrix to the original matrix divided by D₀ give

    Θ(t) = 1/(t√D₀) · Θ(1/(D₀ t)).

Thus the smooth scale B=1/t exchanges with `D₀/B`. At `B=√Δ/2` the inversion has a fixed point. For a sharp disk the transformed radial kernel involves `J₁(s)/s`, with `s=2π√(Bj/D₀)` and `j=Q(u)/N` on the mixed support. One may smooth first, or take a cutoff between consecutive norm levels. The sharp formula does not have a finite Gaussian tail one can discard modulo 16.

In the unscaled integer frequency metric `√Q(u)`, the first nonzero mixed term is at **exactly** `√(Na)`: each mixed lattice has scaled minimum a. The zero mode alone cannot recover the residue. For `Δ=307, ℓ=13, B=1`, the public branches `(N,c)=(35621,301)` and `(14453,77)` give raw counts 4 and 8. Their squared forms are `[7,1,11]` and `[1,1,77]`. The smooth zero mode depends on Δ and the scale and is the same for both. This rules out the direct low-height truncation, while allowing other ways to sum the tail. Standard background: [lattice Poisson summation](https://people.math.harvard.edu/~elkies/M272.19/sep16.pdf).

## Lifting with a supplied factorization of Δ

Given the decoded proper form F (either orientation), find a rational point on `F(x,y)=z²` by passing its doubled Gram matrix and the supplied factorization of its absolute determinant `2Δ` to PARI. The choice of orientation preserves the candidate minima. The lifting routines have no N, p, q, ℓ or c input.

Normalize the point and extend its primitive (x,y) to an oriented integral basis. In this basis F is `[z²,b,d]`. If `gcd(z,b)=1`, Gauss duplication shows that `[z,b,zd]` squares to F. Every proposed seed and every resulting root is verified by composition.

**Repair at conductor and ramified primes.** An arbitrary conic solution can have `gcd(z,b)>1`. The implementation moves it along one rational secant. For a slope S/T in these coordinates, another point is

    X = dT²−S²,    Y=2zST−bT²,    Z=z(dT²+S²)−bST.

Choose S,T by CRT over the primes of Δ so that Z is coprime to Δ after division by the common content. At a good prime p with `p∤z`, take `(S,T)=(1,0) mod p`. At a bad `p^e || Δ`, p divides z and b, and primitivity makes d a unit. Pick w and a unit s with `s²=F(w,1) mod p^(e+1)`, in the new basis. For odd p take w=0 and Hensel lift the square unit d. For p=2 try w modulo 8 to make F(w,1) equal 1 modulo 8, then lift. These local choices exist for a square class: a representative of a half-class can be chosen with leading coefficient prime to 2Δ by choosing a primitive first column locally and lifting it. Duplicating it gives F a primitive square representation coprime to 2Δ. In the bad coordinates that representation has y a unit at p, which yields the asserted local square condition.

Take `S=s−zw, T=1 mod p^(e+1)`. Then `X=wY, Z=sY` at that precision. The product of the expressions for Y using either sign of an exact local square root is `b²−4z²d=−Δ`, so `v_p(Y)≤e`. Removing the gcd leaves a unit Z. This handles nonfundamental orders without assuming the first conic solution is a suitable duplication point.

All remaining roots are the seed composed with the primitive reduced ambiguous classes. These have b=0, b=a, or a=d. Enumerate the factorizations `uv=Δ` from the supplied discriminant factors: the three cases correspond to `(u,v)=(a,4d),(a,4d−a),(2a−b,2a+b)`, respectively. At odd primes a shared factor of u and v would make the form imprimitive, so assign each complete odd prime power to one side. At 2 try all exponent splits. Keep exactly the primitive reduced candidates. The construction tests at most `3(v₂(Δ)+1)2^ωodd(Δ)` candidates; the number of proper lifts returned is `|Cl(−Δ)[2]|`. The square test then costs one trial per distinct minimum per public grid weight. The output size and the cost of obtaining the discriminant factorizations across a band still matter. Background: [Shanks on the square-root coset](https://doi.org/10.1090/S0025-5718-1971-0297737-4), [PARI on supplying the determinant factorization](https://pari.math.u-bordeaux.fr/archives/pari-users-2312/msg00067.html).

## Measurements and reproduction

The repository documents that the historical 192-row raw sample was unavailable at consolidation. Only the two explicit historical collision rows are identified here as reproduced. A separate 192-row replacement sample uses seed 20261126, distinct balanced odd prime pairs in [41,241], and one eligible branch chosen from primes 3 through 19. Every branch and cutoff is in `results.json`. Its grid uses all 359 coprime ordered positive weights with u,v≤24; this is a toy extraction experiment, not the asymptotic canonical grid or a measure of band coverage.

| Exact comparison | Preserved collision pair | Replacement sample |
|---|---:|---:|
| Distinct semiprimes | 2 | 192 |
| Radial cutoffs checked | 40 | 6,377 |
| Norm-quotient square checks to generate all prefixes | 7,593 | 503,876 |
| Positive rows in the public norm-N basis, at the maximum cutoff | 79,676 | 3,151,327 |
| Squared forms decoded and independently checked | 2 | 192 |
| Nonprincipal squared forms | 0 | 89 |
| Factors recovered by the simulated clipped oracle | 2 | 192 |
| Branches with a successful integer square test after lifting | 2 | 157 |
| Integer square tests across distinct minima and grid weights | 1,436 | 113,085 |

Both historical inputs have identical radial counts for all 20 tested cutoffs, including R(5)=24. Both have squared form [1,0,5] and candidate minima 1 and 2. The integer tests succeed uniquely at minimum 1 with (u,v)=(2,3) for N=11009, and at minimum 2 with (7,9) for N=28829. In every successful replacement branch only the true minimum passed. The remaining 35 branches still contained their true minimum among the lifts. In the replacement sample the numbers of branches with 1, 2 and 4 lifts were 73, 101 and 18.

The finite mixed Fourier identity was checked by exact polynomial reduction modulo the cyclotomic polynomial at every frequency for N=15,35,55 (4,475 checks). These three synthetic Fourier controls use Δ=4N−1, for which the public lines have slopes 0 and −1. The transform identity requires splitting and unramifiedness, independently of the largest-lift convention. Against exhaustive class enumeration for all admissible 4<Δ≤1200, all 3,217 square classes gave exactly their full root sets, comprising 5,878 preimages. The initial rational point needed local repair in 502 cases. A separate 33-case supplied-factorization stress test covered 24- to 128-bit discriminants. The 45-bit order Δ=29,682,952,539,240 yielded 2,048 roots; the post-conic enumeration and verification took about 0.02 s in the recorded run. Its output size illustrates why a uniform polylogarithmic per-branch lift budget needs a separate argument.

For the 194 semiprime cases, the recorded batch times were about 0.003 s for the small public discriminant factorizations, 0.044 s for the conics including one GP startup, 0.003 s for inversion, torsion enumeration and composition checks, and 0.028 s for the grid square tests. These exclude enumeration of the public counts and the separate exhaustive controls. Timings are descriptive; the integer operation counts and the exact agreements are the reproducible comparisons.

Unpack the archive at the repository root. From that checkout run:

    python experiments/ellipse_global/benchmark.py --output results-rerun.json

The run used Python 3.12.14, NumPy 2.3.5, SymPy 1.14.0 and PARI/GP 2.19.1. Put `gp` on PATH or set `GP` to the executable. The script imports the base repository's form composition and small enumeration controls. Result fields recording wall time will vary.

A radial reciprocity remains worth pursuing if it batches across rows, across shells or across Fourier support. The one-coordinate prefix reduction sets a concrete test of how strong such a recurrence is. Establishing the conditional almost-all exponent also requires the independent band-coverage audit and an aggregate budget for discriminant factoring and lifts.
