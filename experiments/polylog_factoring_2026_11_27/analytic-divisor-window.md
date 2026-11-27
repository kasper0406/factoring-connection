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
