# A fixed diagonal disk and a native floor-sum evaluator

For an odd integer $N\ge5$, write $M=(N-1)/2$. The relevant scalar is

$$P_N(T)=\#\{1\le x\le T:x^2\equiv1\pmod N\}\pmod2,\qquad 0\le T\le M.$$

It has a public floor-sum evaluator using $O((1+T/N^{1/3})\operatorname{polylog}N)$ bit operations. At the central prefix this still costs a power of $N$. A stationary linear recurrence over the entire period has exponential size in the bit length even after subtracting the two known roots.

## Two fixed diagonal forms encode every prefix

Set $A=N+1$, $d=AN\lceil M^2/N\rceil-1$, $d'=d+AN^2$. Both forms $Q=Ax^2+dy^2$ and $Q'=Ax^2+d'y^2$ are primitive and reduce to $x^2-y^2$ modulo $N$, giving public slopes $\pm1$ and public discriminant roots $\pm2$. For each threshold take

$$B_T=\left\lceil\frac{d+A T^2}{N}\right\rceil.$$

At this cutoff the $y=0$ rows cancel, including for nonsquarefree $N$. The second form has no other row since $d'>NB_T$. In the first form only $y=\pm1$ remain. Indeed $NB_T<d+AM^2+N<4d$, and $0\le NB_T-d-AT^2<N<A$ gives the exact range $|x|\le T$ on both rows. The difference of the signed disk counts is **four times the integer prefix count**, so two queries modulo 8 recover its parity. Both forms stay fixed as $T$ varies. Their cutoffs are of order $N^2$. For the usual squared-minima bound write $t=\sqrt{d/A}\in[M,M+1]$. Then $(t-M)^2+1\le t$ for $M\ge2$, giving $t^2+M^2+1\le Nt$. Hence $B_T<\sqrt{Ad}<\sqrt{|\operatorname{disc}(Q)|/3}$.

For $N$ with at least two distinct prime divisors, opposite roots give $2^{\omega(N)-1}$ roots in the lower half. Removing 1 leaves odd parity. Bisect $[2,M]$, keeping an interval of odd parity. The terminal $x$ gives the proper divisor $\gcd(x-1,N)$. The sibling [all-input construction](general_scope.py) uses residues modulo 4 and changes the form along the path. Neither construction supplies a polynomial-bit evaluator.

## One fixed form suffices for complete factorization

At each odd input $N\ge5$, use only $Q$. First query lower cutoffs $C_s=\lceil As^2/N\rceil$, $1\le s\le M$. Since $M$ is coprime to $N$, $d-AM^2\ge A-1=N$. Therefore $NC_s<d$ and these queries contain only the axis $y=0$, up to the exact radius $s$. On that axis $N\mid x^2$ iff $S\mid x$, where $S=\prod_{p^e\parallel N}p^{\lceil e/2\rceil}$. The signed nonzero count is $2\lfloor s/S\rfloor$. This formula justifies the search without requiring $S$ as an input.

Start with $s=1<S$ and repeatedly double $s$, capping it at $M$. After a previous radius below $S$, the next radius is strictly below $2S$. Its signed residue can only be 0 or 2. The first 2 therefore brackets $S$, and bisection with the same queries returns $S$ exactly. For nonsquarefree odd $N$, $N/S$ has an odd prime factor, so $1<S\le N/3\le M$. This step splits every such input, including odd prime powers.

If no radius produces 2, $N$ is squarefree. For each of the subsequent $B_T$, $NB_T<AN^2$ follows from $d/A<(M+1)^2$ and $(M+1)^2+M^2+1<N^2$. Thus even the largest disk has no nonzero axis point, as an axis solution now needs $N\mid x$. Its signed count alone modulo 8 gives $4P_N(T)$. At $T=M$ the residue is 4 for a prime and 0 for a composite, by the exact count $2^{\omega(N)-1}$. In the latter case bisect as above. Remove powers of 2 directly, use 3 as a prime leaf, and recursively split both proper divisors. Each node takes $O(\log N)$ queries. There are $O(\log N)$ nodes since every prime factor with multiplicity contributes at least one input bit.

All raw coefficients satisfy $A=N+1$ and $d<A(M+1)^2$, the absolute discriminant satisfies $4Ad<4A^2(M+1)^2$, and the cutoffs are at most $B_M<\sqrt{Ad}$. All query descriptions therefore have $O(\log N)$ bits. This is a one-form query reduction. The subenumerative parity algorithm below addresses the prefix component. Using only that algorithm to recover a squarefree split does not require any general disk evaluator.

## A public computation using Euclidean floor sums

For consecutive $x=a+s$, $0\le s<L$, any integral point with $Ny=x^2-1$ is in the parallelogram

$$0\le Ny-(a^2-1+2as)\le (L-1)^2.$$

Take a power of two $L=\Theta(N^{1/3})$ small enough that $2(L-1)^3<N$. Its area is less than $1/2$. Three noncollinear integer lattice points span a triangle of area at least $1/2$, so all candidate points in the parallelogram lie on one line. Its vertical width is less than 1, allowing at most one integer $y$ per $x$.

With $v=2a,u=a^2-1,D=(L-1)^2$, the exact number of candidates in the first $j$ columns equals

$$\sum_{s=0}^{j-1}\left(\left\lfloor\frac{vs+u+D}{N}\right\rfloor-
\left\lfloor\frac{vs+u-1}{N}\right\rfloor\right).$$

Each of the two sums uses the Euclidean floor-sum recurrence. Binary search for the first two occupied columns when present. One point can be tested directly. Two determine a line whose intersection with $Ny=x^2-1$ uses a single quadratic discriminant test and endpoint checks. This works even if the line contains many candidates. Each block takes polynomially many bit operations in $\log N$. Visiting all blocks costs $O((1+T/N^{1/3})\operatorname{polylog}N)$. It searches locally for intersections using only public data. This exact geometrical method offers no complexity advantage from accumulating the result modulo 2.

## What binary Fourier cancellation buys

Let $h(x)=1_{x^2=1\bmod N}$, viewed as an $N$-periodic binary sequence. Write $m_i$ for the pairwise coprime maximal prime powers in $N$ and $e_i$ for their CRT idempotents. In the proof only, its group-algebra polynomial factors as

$$H(X)=\prod_i (X^{e_i}+X^{-e_i})\pmod{X^N-1}.$$

At $X=\zeta^j$, where $\zeta$ has order $N$ over an extension of $\mathbb F_2$, the $i$-th factor vanishes exactly when $m_i\mid j$. Thus exactly $N\prod_i(1-1/m_i)$ additive frequencies survive. Since $N$ is odd, $X^N-1$ has simple roots. Each surviving frequency must be a root of the period sequence's annihilator. Equivalently the minimum stationary linear recurrence order over $\mathbb F_2$ is $N-\deg\gcd(H,X^N-1)$, which equals this support size.

For a semiprime, subtracting the public roots leaves $X^u+X^{-u}$, where $u$ is a unit. Its support comprises all $N-1$ nontrivial frequencies. Differencing consecutive prefix parities recovers that mixed indicator, so a stationary recurrence valid for its periodically extended prefixes also needs order at least $N-1$. These exact bounds apply to a full-period linear representation. They do not bound other algorithms for supplying individual prefixes, such as the local floor sums above.

## Reproduction

Run from the repository root.

    python experiments/ellipse_disk_mod8/native_dyadic.py 11009 5504
    python experiments/ellipse_disk_mod8/native_dyadic_check.py

The central prefix at $N=11009$ has parity 0 using 344 blocks of length at most 16. Independent checks verify 2,079 signed diagonal disk differences, 2,463 native prefixes, and every block on four larger moduli (10,440 blocks). A polynomial gcd check verifies the binary spectral ranks for the 255 odd moduli from 3 through 511.

The independent signed disk oracle also exercises the complete single-form split on these same 255 odd inputs, extracting a proper divisor on all 159 composites and recognizing all 96 primes in 2,921 total queries. Each query checks the squared-minima bound and that its coefficients equal the fixed pair for that input.

## A larger exact block via the standard univariate Coppersmith construction

For the adjacent-offset quadratic floor difference, only the roots in the
block matter. A monic small-root lattice gives a direct public evaluator
on blocks of length proportional to $\sqrt N$. Here are explicit
parameters to avoid a fixed epsilon loss in that exponent. Use
$m=\lceil\log_2 N\rceil$, $H=\lfloor\sqrt N/16\rfloor$, and scan
directly when $H=0$. At a block starting at $a$ set $f(X)=(a+X)^2-1$.
Use coefficient rows of the following polynomials after substituting
$X=HZ$, with $0\le j<m$ and $i=0,1$:
$$N^{m-j}X^i f(X)^j.$$
At any desired offset every row polynomial before scaling is divisible
by $N^m$ in value. Ordered by degree the scaled matrix is triangular,
with determinant $N^{m(m+1)}H^{m(2m-1)}$. Standard LLL gives an
integer combination whose unscaled polynomial $g$ satisfies, for
$|x|\le H$,
$$|g(x)|\le\sqrt{2m}\,2^{(2m-1)/4}N^{(m+1)/2}H^{(2m-1)/2}<N^m.$$
For the strict inequality, divide by $N^m$, insert the upper bound on
$H$, use $N^{1/4}\le2^{m/4}$ and $\sqrt{2m}\le2^m$. The basis has
full rank, so $g$ is nonzero. Every desired offset is an integer root of
$g$ and can be found by rational polynomial factorization and checked
in the block. Dimension, degree and coefficient bit lengths are
polynomial in the input bit length for $0\le a<N$.

Summing disjoint blocks evaluates an interval of length $T$ in
$O((1+T/\sqrt N)\operatorname{polylog}N)$ bit operations. This improves
the direct tangent-block bound above at long cutoffs using a classical
small-root theorem. A central query still uses order $\sqrt N$ blocks.
[The certified implementation](small_root_window.py) centers each block,
tries increasing lattice dimensions, and accepts a row only after its
coefficient 1-norm passes the strict bound. Exhaustive scanning agreed on
126 windows over ten odd moduli from 3 to 124,568,113, including primes,
semiprimes and repeated factors. The windows included known-root
neighborhoods as well as pseudorandom starts. This checks correctness
on those inputs; it is not a runtime comparison with the tangent method.

The same block works on a modular progression $x=a+sj$ for a unit
step $s$. Put $u=s^{-1}\bmod N$ and solve the monic congruence
$(j+au)^2-u^2=0\pmod N$ in the index interval. The determinant and
radius bound are unchanged. The implementation matched scans on 102
progressions across six odd moduli, including wraparound and negative
steps. Thus a public list of polylogarithmically many progressions of
length $O(\sqrt N)$ known to hit a mixed root suffices to split $N$.
For a uniformly random start and fixed unit step on an odd semiprime
with distinct factors, a length-$L$ progression hits either mixed root
with probability at most $2L/N$ by the union bound. Choosing the step
alone cannot improve that random-start guarantee.

**Nonlinear changes of coordinate.** The existing [two-evaluation
identity](../polylog_factoring_2026_11_27/single-quadratic-count.md#two-evaluations-classify-polynomial-maps-on-the-roots-of-one)
also classifies short windows applied after a public polynomial map
$F$ with unrestricted output. Put $A=(F(1)-F(-1))/2$ and
$B=(F(1)+F(-1))/2$ modulo odd $N$. Since $r^2=1$, on every root
$F(r)=Ar+B$. The same formula holds for a rational circuit when each
denominator is a unit at both public signs. If $A$ is zero the window
takes the same answer on every root. A proper gcd of $A$ and $N$ splits
the modulus. Otherwise the preimage on the root set of a cyclic window
$[l,l+w)$ is exactly its intersection with the progression starting at
$(l-B)/A$ with step $1/A$ and length $w$. Both progression parameters
are obtained with two public evaluations and modular inversion.

A check on 600 sampled degree-at-most-six polynomials over six odd
moduli gave 357 unit slopes (all interval preimages agreed with explicit
progressions), 7 constant restrictions and 236 proper gcds. These small
moduli were chosen for checking the cases, not estimating factoring
success rates. This reuses the earlier identity and closes the nonlinear
single-coordinate window variant without a new search algorithm.
