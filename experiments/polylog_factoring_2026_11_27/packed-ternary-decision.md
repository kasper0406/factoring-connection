# Factoring from yes/no questions on one sphere

Given a product $N$ of two distinct odd primes, we previously constructed a pair of vectors with three prescribed inner products such that 48 of its 72 representations reveal a factor. Here the information requested from the evaluator is a single bit. Is an integer $M$ a sum of three squares with the three coordinates in one prescribed residue class? All questions for an input have the same $M$. They can be submitted simultaneously. Each is solvable over the reals and over every ring of $p$-adic integers. The factor comes from which classes have an *integer* point. [Original 72 representations](mixed-ternary-embeddings.md) · [Implementation](packed_ternary_decision.py) · [Independent checks](check_packed_ternary_decision.py)

For a public mixed parameter $x$, the bit sizes are $O(\log N)$ for $0\le x<N$, and the work surrounding the queries is polynomial in $x+1$ and $\log N$. A short $x$ therefore suffices. Given distinct primes in $[R,2R]$, for each fixed $D>1$ all but at most $R^{2/D+o(1)}$ of their products pass a public scan through $(\log R)^D$ or give a factor during the scan. This is unconditional. A bound for every such product follows under GRH. The small-input reference evaluator below enumerates points. Its running time grows as a power of $N$. [Short-parameter proof and distribution](short-jacobi-grams.md)

## Putting two vector conditions in one norm

Set $A=x^2+1$ and assume the Jacobi symbol $(-A/N)=-1$. Let $d$ be the unique prime factor of $N$ at which $-A$ is a nonsquare, and put $r=N/d$. The marked frames satisfy

$$
u\cdot u=A,\qquad u\cdot v=-Nx,\qquad v\cdot v=N^2 .
$$

The 24 coordinate frames are signed permutations of $((x,1,0),(-N,0,0))$. All 48 others have $\gcd(v)=d$ and form one orbit under the 48 signed permutations. Indeed, an isometry fixing two independent vectors is either the identity or reflection in their plane. A signed permutation reflection has primitive normal of length $1$ or $\sqrt 2$, whereas the primitive normal of a noncoordinate frame has length $r$. [Frame count](mixed-ternary-embeddings.md) · [Orbit and trace identities](ternary-frame-moments.md)

Choose $Q$ as the least power of 7 at least $8(x+1)(N+1)$ and let

$$
M=(x+QN)^2+1 .
$$

For a fixed $u$ of squared norm $A$, the map $W=u-Qv$ takes its partners bijectively onto the integral points of squared norm $M$ with $W\equiv u\bmod Q$. To prove the reverse direction, set $v=(u-W)/Q$. The triangle inequality gives $\|v\|\le(\sqrt M+\sqrt A)/Q<N+1$. Expanding the norm yields

$$
Q(\|v\|^2-N^2)=2(u\cdot v+xN).
$$

The absolute value on the right is at most $2(\sqrt A\,(N+1)+xN)<Q$. Both sides vanish. This isolates the required inner product and the second norm exactly.

Enumerate one nonnegative, sorted representative of each signed-permutation orbit of possible first vectors $u$. Two nested loops through coordinates at most $x+1$ and an exact square test suffice. For each test prime $\ell\nmid 210AMN$ and each such $u$, enumerate all triples $z\bmod\ell$ with

$$
z\cdot z=N^2,\qquad u\cdot z=-Nx\pmod\ell .
$$

At least one entry of $u$ is a unit modulo $\ell$. Looping over two entries of $z$ determines the third by the linear condition. For every remaining triple ask just

$$
\exists W\in\mathbb Z^3:\quad W\cdot W=M,\qquad
W\equiv u-Qz\pmod {Q\ell}\,?
$$

The packing identity identifies a yes precisely with a partner $v\equiv z\bmod\ell$. There are at most 72 points in any queried class. Taking one $u$ per orbit preserves every signed-permutation invariant used below.

## Extracting a factor from the bits

For a yes triple set $b_i=z_i^2$ in $\mathbb F_\ell$ and form the two elementary symmetric polynomials

$$
E_2=b_1b_2+b_1b_3+b_2b_3,\qquad E_3=b_1b_2b_3 .
$$

The coordinate frames give $(0,0)$. All noncoordinate frames give the same pair, the reduction of two integers $(E_2^*,E_3^*)$. Keep the unique nonzero pair seen for this prime, or use $(0,0)$ if none appears. In the latter case the noncoordinate frames still exist, so their invariants really reduce to zero. This handles collisions with the coordinate frames.

Choose the consecutive eligible primes through the first one making their product exceed $N^6$. Both integer invariants are nonnegative and smaller than $N^6$, so ordinary CRT recovers them exactly. They satisfy

$$
\gcd(N,E_2^*,E_3^*)=d .
$$

Indeed, writing $v=dw$ makes both invariants divisible by $d$. If $r$ also divided both, all three elementary symmetric polynomials of $w_1^2,w_2^2,w_3^2$ would vanish in $\mathbb F_r$, since the first is $\|w\|^2=r^2$. The monic cubic with those roots would be $T^3$, forcing $r$ to divide every coordinate of the primitive $w$.

The skipped primes divide an integer of size $N^{O(1)}$ for $x<N$. By the Chebyshev bound on the product of the primes up to $t$, all required test primes are $O(\log N)$. Thus there are polynomially many queries for short $x$. At $N=77,x=1$ one noncoordinate $v$ has absolute entries $14,63,42$. The residues of the corresponding two symmetric polynomials already determine 7. The evaluator receives only the norm, the modulus and the residue triple.

The converse also holds on the promised queries. From $(M,m,c)$ recover $Q$ as the 7-part of $m$ and divide $\sqrt{M-1}$ by $Q$ to recover the quotient $N$ and the remainder $x$. Given the factorization of $N$, the Jacobi signs single out $r$. Find a square root of $-A\bmod r$ in expected polynomial bit time. The rank-three unimodular completion at $r$, followed by scaling the second vector by $d$, gives one noncoordinate frame. Its signed permutations and the public coordinate frames give every possible answer, including multiplicity. Thus testing these particular classes and factoring the original inputs give reductions in both directions in the short-$x$ regime. [Explicit completion and bit bound](prescribed-ternary-embedding.md)

## Even parity suffices

The evaluator may instead return the parity of the number of points in the class. In fact it need only return no on empty classes and yes on singletons, with either answer permitted on classes containing more points. To see this, let $S$ be the set of distinct nonzero squared coordinates of one noncoordinate $v$ and $k=|S|\le3$. Every member lies strictly between $0$ and $N^2$. Its positive square root is $d$ times an integer strictly between $0$ and $r$. An axis for $v$ would put the whole frame in a coordinate plane.

For each prime collect all squared coordinates of the yes triples, removing $0$ and $N^2\bmod\ell$. For primes coprime to

$$
\Delta=\prod_{a\in S}a(N^2-a)\prod_{a<b\,\in S}(b-a),
$$

the nonzero signed coordinate values of the one orbit remain distinguishable, including from zero and the axes. Thus with $u$ fixed, its triple $z$ identifies each coordinate of $v$ uniquely. Each such class has zero or one point, and the collected set has size $k$. At any prime dividing $\Delta$ the collected set has fewer than $k$ entries, even if every class with multiple points is accepted. This allows us to recognize the full sets by their maximum observed cardinality.

Use a longer prime list with product exceeding $N^{24}$. Since $|\Delta|<N^{18}$, primes attaining the maximum have product exceeding $N^6$. Reduce the product of the members of each full set modulo its prime and combine by CRT. The recovered positive integer equals $\prod_{a\in S} a<N^6$, and its gcd with $N$ is $d$. Queries remain nonadaptive. In particular the same conclusion holds for a single coefficient modulo two of the product of three series $\sum_{n\equiv c_j\bmod Q\ell} t^{n^2}$, one series for each requested coordinate $c_j$. This parity asks about a prescribed class at a growing modulus.

## Local tests and known bounds

Every submitted residue already has squared norm $M\bmod Q\ell$. At 7, $M\equiv x^2+1\ne0$, and at $\ell$ the choice of prime makes $M$ a unit. In both cases at least one coordinate has invertible derivative, so Hensel lifting extends the prescribed residue to a solution at every higher power. At other primes the congruence imposes no restriction, and $(x+QN,1,0)$ represents $M$. This also gives a real solution. Any test that checks only local solubility therefore accepts *every* query, including the globally empty ones.

For each fixed $u,\ell$ there are exactly $\ell-\left(\frac{-A}{\ell}\right)$ local candidate triples. Translation by $(-Nx/A)u$ reduces the two constraints to a binary quadratic form of nonzero value $N^2/A$ on $u^\perp$. Its determinant has squareclass $A$. A split plane has $\ell-1$ points of a specified nonzero norm, and an anisotropic plane has $\ell+1$ by the finite-field norm map. At most 72 of these candidates lift globally for this $u$. A rule that always answers no can therefore have high accuracy on large test primes while discarding every factor-revealing answer.

For short $x$ and $m=Q\ell$ one has $m/M^{1/4}\asymp\ell\sqrt{x+1}$. Ellenberg, Michel and Venkatesh's Theorem 1.4 gives lifting for almost all residue classes in a range $m\le M^{1/4-\nu}$, under conditions including squarefreeness of $M$. Here we need a specified class just above that scale and $M$ may have a square divisor. Sardari proves general congruence-constrained sphere decision NP-complete under randomized reductions. His polynomial-time construction for special four-square congruences assumes an integer-factoring algorithm and a further conjecture. Neither result supplies an algorithm or a complexity lower bound specific to these ternary instances. [Ellenberg, Michel and Venkatesh, Theorem 1.4](https://arxiv.org/pdf/1001.0897#page=6) · [Sardari, *Complexity of Strong Approximation on the Sphere*](https://doi.org/10.1093/imrn/rnz233)

The earlier predictor work already decoded supplied congruences by CRT. The ternary frame count, rational rotations and quartic traces concern the same underlying representations. The additional statement here is the reduction to ordinary scalar congruence decisions with a common near-square norm and bounded multiplicity, including parity. A new evaluation method must distinguish individual classes after all the local tests have passed. [Earlier CRT ledger](README.md) · [Rational rotations](ternary-frame-rotations.md) · [Quartic traces](ternary-frame-moments.md)

## Reproducing the decisions

The native reference takes only $(M,m,c)$ and enumerates the first two coordinates of $W$ in their prescribed progressions. For each pair it tests the remaining integer for being a square and checks both signs. It uses $O((1+\sqrt M/m)^2)$ square tests per query. Here $\sqrt M/m\asymp N/\ell$. The same enumeration counts points for the parity experiment.

| $N$ | Selected $x$ | Factor | Exact-existence queries | Parity queries |
|---:|---:|---:|---:|---:|
| 15 | 0 | 3 | 116 | 832 |
| 35 | 0 | 7 | 180 | 1,248 |
| 65 | 4 | 5 | 472 | 3,516 |
| 77 | 1 | 7 | 248 | 1,906 |
| 91 | 0 | 7 | 264 | 2,004 |
| 183 | 0 | 3 | 264 | 2,352 |

The six inputs use 1,544 exact-existence queries and 11,858 parity queries. Every point count on the longer menus was cross-checked against an independent enumeration of frames. At $N=91,x=0,\ell=17$ one queried residue contains exactly two points, so existence and parity give different answers. Both runs recovered 7. [Run the checks](check_packed_ternary_decision.py)

A shorter decoder reads a single signed coordinate when $x\ge2$ is coprime to $N$. Its existence queries need test primes of product only $>2N$. A direct evaluator for parity also saves work by testing for a coordinate symmetry and then searching one coordinate on a conic. Finally, every residue at a fixed test modulus lies in the same proper spinor genus. The exact differences among the individual counts are coefficients of cusp forms orthogonal to unary theta series. [Signed-coordinate decoder](packed-ternary-linear.md) · [Native parity evaluator](native-ternary-parity.md) · [Spinor-genus calculation and explicit occupied and empty pair](ternary-coset-spinor-genus.md)
