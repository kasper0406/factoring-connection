# Mixed cubic labels from short rational seeds

There is an input-only family where the predetermined seed $\lfloor(8N)^{1/3}\rfloor$ produces a mixed cubic label and splits the input. The seed height is of order $N^{1/3}$, the first possible scale for this mechanism. Two elementary integer square roots recover the same factors on this family. More generally, at the first feasible scale, reconstruction of a homogeneous quadratic gives a factor even if the seed has a denominator. These statements describe the arithmetic promise behind this example. [Witness fusion](shanks-witness-fusion.md) · [Implementation and controls](shanks_mixed_slope.py)

## A two-variable family detected from the input

Take coprime positive $X,Y$ with $3X/2\le Y\le2X$. Define
$$
 A=8XY+1,\qquad P=AY-X,\qquad Q=AX+X+Y,\qquad N=PQ.
$$
For example take $X$ even, $Y$ odd to make both factors odd. The identity
$$
 8N-A^3=8A(Y^2-X^2)-A-8(X^2+XY)
$$
lies strictly between zero and $3A^2$: the lower bound follows from $Y^2-X^2\ge5X^2/4$, and the upper bound from $Y^2-X^2\le3XY/2$. Consequently the seed $A=\lfloor(8N)^{1/3}\rfloor$ is public. The original rational direction $X/Y$ and the seed $A$ agree at $P$, whereas $A$ agrees with the branch $-(X+Y)/X$ at $Q$ after clearing denominators.

For the root-direction Shanks parameter, let $F=X^3-3XY^2-Y^3$, $D=XY(X+Y)$ and $D_0=X^2+XY+Y^2$. The homogeneous Shanks cubic evaluated at $A$ with parameter $(F,D)$ factors into the three branch differences, two of which are $P$ and $Q$. Thus $(F,D)$ is in the public parameter lift lattice whenever the seed denominators are invertible. The norm of this primitive parameter vector is $D_0^3$. Since $D_0\le(7/2)XY$ and $N>A^3/8>64(XY)^3$, we have $D_0^3<N$. The lattice has covolume $3\sqrt3 N/2$ in the norm plane, so $(F,D)$ is its unique shortest direction. The cubic is rationally split, and its two local labels yield proper gcds. A preliminary nonunit at $A,A+1,A^2+A+1$, or at the short denominator, is already a proper factor when it occurs on these growing inputs. One can alternatively factor the rationally split cubic and compare each primitive rational root directly with $A$ by gcd. [Shortest-lift construction](shanks-conductor-relations.md#rational-lifts-reduce-the-uniform-norm-bound-from-quadratic-to-linear)

There is a simpler direct recognition. From $A$ recover $U=(A-1)/8=XY$ and
$$
 E=\left\lceil\frac{8N-A^3}{8A}\right\rceil=Y^2-X^2.
$$
The ceiling is exact because the difference between the true integer and this fraction is $1/8+(X^2+XY)/A$, strictly between zero and one. Then $\sqrt{E^2+4U^2}=X^2+Y^2$ gives $X,Y$ by integer square roots. Verify the resulting factors by division. This two-parameter family has at most $O(T^{1/3})$ values with $N\le T$ since $X,Y=O(T^{1/6})$ and comparable. The construction gives no density claim for prime $P,Q$. Four selected semiprimes of 56, 56, 86 and 126 bits split via the cubic, the direct identity and the general quadratic reader below. A further 218 parameter pairs passed the two input-only family checks. [Executable controls](shanks_mixed_slope.py)

## Which rational seeds can have mixed labels?

Write the modular seed as $A/B$ with $\gcd(A,B)=1$ and $B$ invertible modulo $N=pq$. Let $X/Y$ be a primitive rational direction whose ordinary Shanks parameter has a globally split short lift. Use $H=\max(|X|,|Y|)$ and $S=\max(|A|,|B|)$. The three branch comparisons are
$$
 L_0=AY-BX,\qquad L_1=AX+B(X+Y),\qquad L_2=A(X+Y)+BY=L_0+L_1 .
$$
A mixed seed means precisely that $p,q$ divide **different** $L_i$. At a good parameter the three local branches are distinct. A zero integer $L_i$ would be a globally shared label, so each of the two mixed comparisons is nonzero. Since $|L_i|\le3SH$, both input primes are at most $3SH$. The lift bound $D_0^3\le3N$ implies $H=O(N^{1/6})$. Thus for balanced inputs no representation of a mixed modular seed can have height $o(N^{1/3})$, even if a particular reduction algorithm chooses another representative. The preceding family attains this order with $B=1$.

There is a separate coverage limit for fixed menus. Let both primes lie in $[R,2R]$, allow $H=O(R^{1/3})$, and allow any set of $M$ rational directions with polynomial height $R^{O(1)}$ across this band. A fixed direction and a fixed $(X,Y)$ have at most a constant number of successful semiprimes: each nonzero $L_i$ is $R^{O(1)}$ and has only $O(1)$ prime divisors at least $R$. There are $O(H^2)$ target directions. The whole set therefore covers $O(MH^2)$ pairs out of $\asymp R^2/(\log R)^2$. This bound permits an arbitrary choice among the menu after inspecting the input. In particular, even the full integer seed range up to $N^{2/3-\epsilon}$ covers a vanishing fraction of comparable prime pairs through this particular mixed-label mechanism. This count concerns the union of all allowed directions over the band. A short list computed afresh from each $N$ need not have a short union. Ordinary rational reconstruction at an unknown modular seed gives height $O(\sqrt N)$ and allows many more numerator-denominator pairs.

## Native completion from either side at the first feasible scale

**Given the seed $A/B$.** If the labels are mixed, two primitive linear forms in $(A,B)$ have values $kp,\ell q$, for some nonzero signed $k,\ell$ with $|k|,|\ell|\le K=3SH/\min(p,q)$. Under the balanced promise use the public lower bound $\sqrt{N/2}$ in setting $K$. Their product has the form
$$
 c_2 A^2+c_1 AB+c_0 B^2=mN,\qquad 0<|m|\le K^2,
$$
with $|c_2|,|c_0|\le4H^2$ and $|c_1|\le8H^2$. Swap coordinates if necessary to make $|A|=S$. At each signed $m$, first enumerate the bounded $c_0=m N/B^2\bmod |A|$. After subtraction and division by $A$, enumerate $c_1=((mN-c_0 B^2)/A)/B\bmod |A|$. Determine $c_2$ by exact division. At the genuine coefficients the homogeneous quadratic factors into the two branch forms. A square-discriminant test supplies their primitive linear factors; their evaluations recover a divisor by gcd. If a multiplier itself contains an input prime it may be checked during enumeration; at the limiting scale all multipliers are bounded. Successful branch forms have nonzero leading and trailing coefficients at good ordinary roots. The work is
$$
 O\bigl((1+K^2)(1+H^2/S)^2\operatorname{poly}(\log N,\log S,\log H)\bigr).
$$
At $S\asymp N^{1/3}$ and $H=O(N^{1/6})$ both enumeration factors are bounded. This is homogeneous digit-polynomial reconstruction in a rational base and applies to any two short linear factors at that scale. The cubic relation between their coefficients becomes relevant for larger $K$, where the Shanks shortest-lift test still checks its narrower promise in polynomial time. The checker also exercised 267 coprime branch-pair cases with nonintegral seeds. [Previous carry-polynomial reconstruction and prior art](carry-polynomial.md)

**Given the short root direction $X/Y$, without the seed.** Rotate the global root by its order-three substitution, and if necessary interchange $p,q$, to put the two labels in positions $L_0=kp,L_1=\ell q$. This rotation changes the height by at most a factor of two. Solving the two linear equations gives
$$
 A=\frac{(X+Y)kp+X\ell q}{D_0},\qquad B=\frac{Y\ell q-Xkp}{D_0}.
$$
Since $X,Y$ are units modulo $D_0=X^2+XY+Y^2$, integrality of $B$ implies integrality of $A$. If $D_0$ is coprime to $N$, the condition is equivalent to
$$
 X k p^2=Y\ell N\pmod {D_0}.
$$
Conversely any such congruence produces the displayed integer pair, before primitivity and valid seed checks. Thus for a prime $D_0>K$ there are at most two classes for $p$ per signed pair $k,\ell$. To include either factor as the smaller prime, add $N/r\bmod D_0$ for each candidate residue $r$. The three rotated root directions suffice for the unknown branch locations. For factored composite $D_0$, charge root extraction and the actual list size, dividing a common $\gcd(k,D_0)$ first when it also divides $\ell$. A proper gcd with $N$ already splits the input. When $D_0\gg N^{1/4}$ and there are polynomially many classes, the balanced known-residue algorithm finishes in polynomial bit time, with the usual randomized square-root extraction or supplied square-root data at large prime moduli. For $N=48{,}652{,}964{,}904{,}198{,}163$, $X/Y=246/371$ gave four classes, including reciprocal classes, at the prime $D_0=289423$ with $k=\ell=1$. The native reader recovered $179612351$. The label was supplied for this last test. [Native residue-list theorem](hh-residue-lists.md) · [Checks](shanks_mixed_slope.py)
