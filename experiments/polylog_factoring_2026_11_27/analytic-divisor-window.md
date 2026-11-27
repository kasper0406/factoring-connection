# A Poisson kernel for divisor windows

For an interval I of W positive integers x≤N, set
w_x=exp((−1+2πiN)/x), with r_x=|w_x|. The normalized real part of
(1+w_x)/(1−w_x), multiplied by (1−r_x)/(1+r_x), is the positive weight

    K(x) = sinh²(1/(2x)) / (sinh²(1/(2x)) + sin²(πN/x)).

It equals one on divisors. At any other x the nearest nonzero residue
to N mod x has absolute size at least one, so |sin(πN/x)|≥2/x.
Also sinh(1/(2x))≤1/x for x≥1. Hence K(x)≤1/5. With k the smallest
positive integer satisfying 5^k≥4W, the sum of K(x)^k differs from the
exact divisor count by at most 1/4. An absolute error smaller than 1/4
in this sum suffices for rounding. Pointwise evaluation to O(1/W)
absolute precision costs polynomially many bits and operations in
log N: the denominator is at least 1/(4x²), phases can be reduced
modulo x, and k=O(log(2W)). Direct evaluation still visits W positions.

On the intervals [2, floor(sqrt(N))] at N=77,91,10403,2989601243, the
computed sums were respectively 1.000312171334, 1.000350365958,
1.00002006019 and 1.000000000175, with exact count one in each case.
The last interval contained 54,676 integer positions.

For a possible summation transform, the reciprocal phase N/x may be
replaced by N/x+x at integer x. Near sqrt(N) the linear change cancels
and x+N/x−2sqrt(N)=(x−sqrt(N))²/x. This is also the geometry exploited
by Fermat's method. A new algorithm would need a bound on the work of
summing the sharp powered kernel *and* an absolute error below 1/4.
The [prior analytic audit](../../notes/conn_ana.md) addresses the
difficulty of recovering an exact divisor jump from smooth estimates.
This calculation supplies a pointwise gap; it does not supply the
global fast transform or assert priority for the kernel.

## Moving the sum into the complex plane

For positive integer N, extend the kernel meromorphically away from zero,

$$
K_N(z)=\frac{\cosh(1/z)-1}{\cosh(1/z)-\cos(2\pi N/z)}.
$$

The denominator vanishes when the two cosine arguments agree up to sign
and a multiple of $2\pi$, using $\cosh u=\cos(iu)$. Thus in the right
half-plane the poles are exactly
$z_j^{\pm}=(N\pm i/(2\pi))/j$ for positive integers $j$.
They are simple because the two solution families are disjoint and the
sine at each solution is nonzero. The numerator is nonzero there since
$1/z$ has positive real part. Taking a positive integer power $k$
raises their order to $k$.

A contour integral with sampling factor $\pi\cot(\pi z)$ recovers the
integer samples by residues in a sufficiently thin rectangle. Moving its
upper edge above the upper poles introduces one residue per $j$ with
$N/U<j<N/L$ for vertical sides at $L,U$ (choose endpoints off poles
and integers). For $L=\sqrt{N/2}$ and $U=\sqrt N$, that is
$(\sqrt2-1)\sqrt N+O(1)$ distinct poles. Their real parts $N/j$
parameterize complementary factor candidates. For $N=2,989,601,243$
the integer indices run from 54,678 through 77,325, giving 22,648
upper poles. Each has height $1/(2\pi j)$.

This rules out a polylogarithmic cost for *listing* all crossed poles
in this window. It leaves open a collective evaluation of their
contributions. A candidate identity would need to bound total bit cost
and absolute error uniformly, including endpoint contributions.
