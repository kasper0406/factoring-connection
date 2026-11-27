# Scope and cost audit of the general-input reductions

This audit separates the reductions that cover arbitrary odd inputs from the
balanced-semiprime fixed-cap evaluator. It also checks whether the
$N^{2/5+o(1)}$ tangent evaluator can be inserted into a general reduction.
The conclusion is negative for that insertion. The general reductions remain
valid oracle reductions, and no polynomial-bit evaluator is established.

## Deterministic fixed-window route

For an odd current modulus with at least two distinct prime divisors, the
fixed-window construction needs $O(\log N)$ adaptive parity queries. In its
two-form version, each parity query is the difference of two signed counts
modulo four. In the common-discriminant version, one midpoint count supplies
the invariant axis baseline and each later parity uses one signed count.
The bisection does not assume monotonicity. It maintains two starts with
different parities. At adjacent starts, exactly the removed and inserted
endpoints can change the answer, and direct square and gcd tests identify a
proper divisor.

Prime and prime-power leaves have no mixed square root of one. Polynomial-time
primality testing and exact integer-root tests distinguish them. Removing
powers of two and recursing therefore gives complete factorization. A binary
factorization tree has at most $O(n)$ nodes for
$n=\lceil\log_2 N\rceil$. Charging $O(\log M)$ queries at a node of
value $M$ gives the conservative complete bound $O(n^2)$ signed-count calls.
All coefficients, cutoffs and class exponents have $O(n)$ bits.

For the common-discriminant form, $L=\Theta(N^2)$,
$D=\Theta(N^6)$, the numerical cutoff is $\Theta(N^5)$ and the
absolute discriminant is $\Theta(N^{10})$. These are powers of the value of
$N$, but their binary descriptions still have length $O(n)$. All query forms
are computable powers in a public cyclic subgroup of order $L$. Computing the
power or composing the class is polynomial-bit work. It does not count the
short vectors in that class.

The reference programs reproduced 6,478 shifted windows on all odd inputs from
5 through 319. Both variants split all 84 inputs having at least two distinct
prime factors and treated the other 74 inputs as prime or prime-power leaves.
The cyclic encoding independently reproduced all 6,478 query classes.

## Nonadaptive and weak-output routes

Affine isolation covers any number of distinct odd prime divisors and repeated
prime powers. The strengthened pair bound uses a constant number of trials at
each of $O(n)$ public dyadic widths, hence $O(n)$ **weighted** queries for
constant success probability per split. This supersedes the older $O(n^2)$
count in affine-isolation.md. Amplification costs
$O(\log(1/\varepsilon))$. Reconstructing an isolated coordinate
nonadaptively from unweighted residues costs $O(n^2)$ counts per trial and
$O(n^3)$ counts per split. A conservative complete-factorization bound is
$O(n^4)$ unweighted counts. These are query counts, not evaluator running
times.

The quarter-window decoder and its weak-prediction Fourier reduction apply
when the current modulus has exactly two distinct prime divisors, allowing
prime powers. They do not supply a uniform decoder when three or more distinct
prime factors remain. Under their stated promise, an inverse-polynomial
prediction advantage in $n$ is enough, and the final gcd verifies every
candidate.

The one-sided quadratic gap comparison is the broadest relaxed target. It
inherits the affine input isolation and therefore covers repeated factors and
arbitrarily many distinct primes. Each oracle input has ten bounded integer
coordinates, one rank-two quadratic equation, and $O(n)$-bit coefficients.
The comparison needs additive resolution on the order of $N$ around an ambient
count as large as order $wN^3$. Fixed dimension alone therefore does not make
ordinary volume approximation sufficient. Exact Barvinok counting for a
linear rational polytope also does not evaluate the added quadratic equality.
The repository reduction proves that a polynomial-time evaluator for these
comparisons would yield randomized polynomial-time factoring. It gives no
such evaluator.

There is also a sharp reason to avoid strengthening this target back to exact
affine-window counting. A parsimonious reduction from Subset Sum encodes every
choice as a CRT sign pattern and uses a public unit scale to place its signed
sum on a circle. One affine half-period window gives the cumulative subset
count. Thus arbitrary public affine root-window counts are #P-hard, and their
parities are $\oplus$P-hard, even when the modulus factorization is supplied.
This statement leaves the promised one-sided gap comparison open. It also
leaves the unscaled deterministic windows and the doubling-only quarter family
open because the hardness construction chooses the CRT scale independently at
every prime-power component. [Proof and controls](affine-root-hardness.md)

## Why the tangent evaluator does not transfer

The $N^{2/5+o(1)}$ tangent-progressions algorithm uses the special balanced
branch

$$
w^2+\Delta y^2=4Nj,\qquad
\Delta=O(\sqrt N),\quad
j\le A=O(N^{1/4}).
$$

It factors the public discriminant and sums a bounded set of norm levels. The
fixed-window forms instead deliberately have only the rows $y=0,\pm1$ inside
the disk. On the nonaxis rows their congruence is the clipped problem
$x^2=1\pmod N$. Substituting their parameters into the shell method gives up
to $\Theta(N^5)$ norm levels and a public diagonal term
$D=\Theta(N^6)$ whose factorization is not supplied. Exploiting the three-row
geometry returns to the clipped-root problem. The certified univariate
Coppersmith cover currently costs $\Theta(\sqrt N)$ blocks for a central
window.

The balanced scalar evaluator also requires a product of two distinct primes
with comparable size for its shape lemma and decoder. Multiplying a general
input by a known prime does not create that promise. Finally,
$N^{2/5+o(1)}$ is asymptotically slower than the repository's cited rigorous
general $N^{1/5+o(1)}$ algorithm. Taking a minimum with the tangent evaluator
therefore gives no improved general bound.

## Priority

For a uniform polynomial-bit algorithm, the highest-value relaxed target is
the **one-sided affine-isolated quadratic gap comparison**. It covers all odd
inputs after the explicit prime-power branch, tolerates a gap and returns a
heavy frequency whose gcd is verified. A proposed evaluator must exploit the
inverse-graph structure and integer interval boundaries together. Generic
fixed-dimensional volume estimates, separate local CRT counts, or the
balanced tangent routine do not meet the required additive resolution. The
new exact-count hardness result makes the promised gap and weak prediction
essential parts of this target, rather than conveniences to discard.

The fixed-window parity remains the smallest exact falsification target. It
should be used to reject proposed exact count mechanisms before lifting them
to the ten-coordinate comparison.
