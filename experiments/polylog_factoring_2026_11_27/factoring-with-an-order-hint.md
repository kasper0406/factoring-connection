# Factoring with a partially shared group order

Given only a balanced semiprime $N=pq$ and an integer $e$ of size about $\sqrt N$, one can deterministically find a factor in polynomial bit time, for sufficiently large $N$, if $e$ shares a divisor larger than $N^{1/(2\sqrt2)+\eta}$ with either prime's predecessor or successor, for fixed $\eta>0$. The algorithm requires neither the shared divisor nor the prime factorization of $e$. For example, $N^{0.36}$ suffices asymptotically at this exponent size. The result converts a partial order hint into a factor. Producing a hint with this property on every balanced input remains open. [Weighted calculation](weighted_order_hint.py) · [Exact checks](check_weighted_order_hint.py)

The repository already handled a supplied exponent missing at most a polynomial-in-bit-length cofactor of $p-1$ by padding the exponent. At the threshold above, that cofactor can grow as a positive power of $N$. The lattice bounds used below come from classical small-root recovery and weighted Chinese remaindering. [Previous scalar-order reduction](bounded-quadratic-roots.md#completing-an-exponent-that-covers-one-hidden-prime) · [Previous CRT decoding audit](../../notes/crtdecode_results.md) · [Relation to powered interval sampling](power-window-orders.md)

## A single approximate multiple

Write $r$ for either hidden prime and $s=N/r$ for the other. Suppose $g=\gcd(e,r-1)$. Then
$$
 gs\mid eN,\qquad gs\mid N-s=s(r-1).
$$
For $g=\gcd(e,r+1)$ replace $N-s$ by $N+s$. The public center is $N$, the modulus is $eN$, and the unknown offset is always within $\sqrt{2N}$. The affine small-root algorithm lists every qualifying gcd at that radius. The gcd of the qualifying value with $N$ gives $s$, since $\gcd(N,N\pm s)=s$. The offset zero can give the whole original modulus and is discarded. The same identity for a fixed linear expression is $gs\mid aN\pm bs$ when $g\mid e$ and $g\mid ar\pm b$. [Affine order-hint decoder](affine_order_hint.py) · [Checks](check_affine_order_hint.py)

For the size calculation write $e=N^{\alpha+o(1)}$, $g=N^{\gamma+o(1)}$, with fixed $\alpha>0$. An affine polynomial modulo a divisor of size $D$ of modulus $M$ admits polynomial-time small-root recovery at radius strictly below $\exp((\log D)^2/\log M)$ by a fixed power of $M$. Here $D=gs=N^{1/2+\gamma+o(1)}$ and $M=eN=N^{1+\alpha+o(1)}$. Sufficient coverage of the partner prime follows from
$$
  \frac{(1/2+\gamma)^2}{1+\alpha}>\frac12 .
$$
In particular if the missing cofactor $(r\pm1)/g=N^{o(1)}$, then $\gamma=1/2$. Every $e\le N^{1-\epsilon}$ works under that promise for fixed $\epsilon>0$ and sufficiently large $N$. The interval-concentration argument forces precisely this small-cofactor property at some prime when powered centers have inverse-polynomial split probability. [Affine bound and dyadic parameter selection](affine-divisor-scales.md) · [Coset concentration](power-window-orders.md#what-a-frequently-successful-window-implies)

## Giving the two divisibility conditions different weights

We can also seek $r$ directly. Check $\gcd(N,e)$ first and use the following construction for coprime moduli. Fix a sign $u\in\{-1,1\}$. Choose $C\equiv0\pmod N$, $C\equiv u\pmod e$ by CRT. At the desired root $x=r$, the value $f(x)=x-C$ is divisible both by $r$ and by $g=\gcd(e,r-u)$. These divisors are coprime. This is partial agreement at two known composite moduli. The classical soft-decision CRT decoder treats agreement at supplied coordinates. Here divisibility by unknown factors at those coordinates is enough, as the following determinant calculation shows. [Guruswami, Sahai and Sudan](https://people.csail.mit.edu/madhu/papers/2000/gss-conf.pdf) · [Related small-height lattice methods](https://cr.yp.to/ntheory/smallheight-20080502.pdf)

Use $D$ rows indexed by $0\le i<D$ with coefficients from
$$
 f(x)^i\,N^{\max(m-i,0)}e^{\max(n-i,0)}.
$$
Scale column $j$ by $X^j$, for $X=\lfloor\sqrt{2N}\rfloor$. Evaluations at $r$ are multiples of $r^m g^n$. The determinant is
$$
 N^{m(m+1)/2}e^{n(n+1)/2}X^{D(D-1)/2}.
$$
Using public lower bounds $L=\lceil\sqrt{N/2}\rceil$ and $G\le g$, a reduced row with $D\,\|v\|^2<L^{2m}G^{2n}$ gives an integer polynomial vanishing at $r$. Isolate all integer roots and verify the divisors. The same lattice allows a public center $s_0$ by replacing $x$ by $s_0+t$ and bounding $|t|$.

In the LLL determinant bound choose $m/D$ near $\log L/\log N$ and $n/D$ near $\log G/\log e$. The leading terms, after letting $D$ be sufficiently large for a fixed exponent margin, give the sufficient condition
$$
 \log X < \frac{(\log L)^2}{\log N}+\frac{(\log G)^2}{\log e}
$$
with a gap of at least a fixed positive multiple of $\log N$. Indeed the two normalized multiplicity contributions are $2z\log L-z^2\log N$ and $2w\log G-w^2\log e$. They attain their maxima at the displayed choices. Integer rounding and the linear terms in the determinant exponents cost $O((\log N+\log e)/D)$ in this comparison, and the LLL approximation contributes $O(1)$ for fixed $D$. For a balanced input this condition reads $1/2<1/4+\gamma^2/\alpha$, or $\gamma>\sqrt\alpha/2$. For fixed $\alpha$ it is at least as strong as giving both coordinates the same multiplicity, by Cauchy–Schwarz. All finite certificates are checked with integers. The threshold is sufficient, and a failure to reach it leaves open what other algorithms might recover.

The difference between the two determinant tests can be observed in a planted 383-bit example. Here $e=3\cdot 2^{113}$ and $g=2^{113}$ at one prime. With the full balanced interval, the weighted basis finds that prime in dimension seven with multiplicities three and six. The test verifies two integer inequalities showing that *no dimension* passes the sufficient determinant criterion when the multiplicities are required to be equal, for the same public lower bounds and search radius. This statement concerns the guarantee from the determinant. An individual reduced row can outperform its guaranteed length. A public value of $G$ can also be selected by minimizing the integer bound that passes the test at a chosen dimension. [Weighted implementation and controls](check_weighted_order_hint.py)

With a public center $s_0$, a radius $N^{\beta+o(1)}$ suffices for any candidate it contains when $\beta<1/4+\gamma^2/\alpha$. Thus even $e>N$ can help on inputs where a smaller interval is available. Centering at $\lfloor\sqrt N\rfloor$ gave a finite example with a qualifying prime farther than $N^{1/4}$ from the center. It adds an interval assumption to the global result.

### Combining several weak congruences

Suppose the public hints give pairwise-coprime moduli $e_i$, each coprime to $N$, residues $a_i$, and individual lower bounds $G_i\le g_i=\gcd(e_i,r-a_i)$ at the same hidden divisor $r$. Set $e_0=N$, $G_0=L$ and $a_0=0$. Choose the common CRT representative $C$ of the $a_i$. Row $j$ now has polynomial $(x-C)^j\prod_i e_i^{\max(m_i-j,0)}$. Its value at $r$ is divisible by $r^{m_0}\prod_{i\ge1}g_i^{m_i}$. Choosing each multiplicity independently gives the sufficient bound
$$
 \log X < \frac{(\log L)^2}{\log N}+\sum_{i\ge1}\frac{(\log G_i)^2}{\log e_i}
$$
with fixed exponent slack, for a fixed number of hints of polynomial size in $N$. This uses the multiplicities of classical weighted CRT decoding. The implementation verifies partial divisibility inside each composite coordinate. [Guruswami, Sahai and Sudan](https://people.csail.mit.edu/madhu/papers/2000/gss-conf.pdf) · [Multi-hint implementation](weighted_order_hint.py)

If only a lower bound $H$ on $\prod g_i$ is given, combining the hint residues by CRT gives one modulus $E=\prod e_i$ and one contribution $(\log H)^2/\log E$. This assigns all hint coordinates a common multiplicity. A preselected unequal multiplicity vector can also be run without individual bounds. For a reduced row $v$ the candidate $r$ is guaranteed to appear whenever $(r^{m_0}\prod g_i^{m_i})^2>D\|v\|^2$. The row norm is public, while the individual gcds in this last condition depend on $r$. For a fixed number $k$ of hints, trying all $D^{k+1}$ multiplicity vectors exhausts the search at dimension $D$.

A planted 257-bit balanced semiprime has one hint at $r-1$ and another at $r+1$. From the public bounds $2^{61}$ and $2^{127}$ the decoder recovers a prime in the full balanced interval at $D=7$ with multiplicities $(3,6,1)$. Either singleton fails the sufficient determinant criterion in every dimension. Combining both hints with their **exact** gcd product still fails that criterion in every dimension with a common hint multiplicity. The check proves these exclusions by integer power comparisons, certifies primality of the planted factors, and repeats recovery with only the chosen weights. The moduli were constructed using the target prime. [Reproduce the three comparisons and the recovery](check_multi_order_hints.py)

There is a list-size obstruction to choosing a different optimal weight vector for every candidate integer. Take $t^2$ distinct primes in $[t^4,2t^4]$ for sufficiently large $t$. There are at least $t^t$ products $d$ of $t$ distinct primes from the list. Choose a prime $B$ between $t^{t/4}$ and $2t^{t/4}$. Some nonzero class $c\bmod B$ contains at least $t^t/B$ of the products. Take $A\equiv c^{-1}\bmod B$ coprime to every listed prime, with $1\le A\le(t^2+1)B$. Among $t^2+1$ consecutive representatives of that class each listed prime rules out at most one. Put $X=2A(2t^4)^t$. For every product in the class, $x=2Ad-1<X$ equals $-1$ modulo the primes in $d$ and modulo $A$, and equals $+1$ modulo $B$. Omit the coordinate $A$ if it equals one. These are pairwise coprime coordinates. Full matches alone have weight at least $\log(AdB)\ge\log X+\log B-(t+1)\log 2>1.05\log X$ for all large $t$. Thus a single list of signs and moduli of polynomial bit length in $t$ admits superpolynomially many such integers. Listing every integer by its individually optimized score would already take superpolynomial time to print. The candidates in this construction are general integers. Requiring them to divide a supplied semiprime is an additional restriction. [Earlier CRT list analysis](../../notes/crtdecode_results.md) · [Availability of supported signs at actual prime divisors](sparse-order-support.md)

## Recovering the rational order relation

A second decoder recovers the small ratio hidden in the order hint. For $g=\gcd(e,r-1)$, set $h=(r-1)/g$, $k=e/g$. Then $\gcd(h,k)=1$ and $eh+k=kr$. For $r+1$, use $eh-k=kr$. The integer $kr$ itself divides $eN$ and is at least $e/2$. At public bounds $h\le H, k\le K$, homogenizing the affine basis to degree $D-1$ scales column $j$ by $K^j H^{D-1-j}$. The usual affine determinant for modulus $M$ now has radius $HK$. Dehomogenizing one certified short row gives a polynomial with rational root $k/h$. The implementation isolates its roots and verifies each candidate by exact integer division. Castagnos, Joux, Laguillaumie and Nguyen give the homogeneous unknown-divisor theorem for equal coordinate bounds. The separate column scales here implement unequal bounds. [Published homogeneous theorem, Section 3](https://www.lirmm.fr/~laguillaum/NiceFactor.pdf) · [Projective decoder](projective_order_window.py) · [Checks](check_projective_order_window.py)

With $M=N$ and divisor $r$, the quarter-root bound covers $HK\le\lfloor N^{1/4}/8\rfloor$. With $M=eN$ and divisor $kr$, it covers larger exponents. When $h=N^{o(1)}$ and $\alpha\ge1/2$, the unknown product has exponent $\alpha-1/2+o(1)$ in base $N$. The second modulus gives capacity $\alpha^2/(1+\alpha)+o(1)$. The difference is $(1-\alpha)/(2(1+\alpha))$, positive for $\alpha<1$. Public boxes can be selected using $k\le eH/(L-1)$ and the exact dyadic divisor-size bounds. In a planted 1,043-bit control the recovered order cofactor was $2^{52}-1$. The exponents in these controls were constructed using one prime to verify the theorem's promise. The decoders themselves received only $N,e$ and the sign.

## Extension to other polynomial orders

Let $P\in\mathbb Z[T]$ be fixed of positive degree $d$, with $P(0)=\pm1$. For $N=rs$ define the monic polynomial
$$
 F(X)=P(0)X^d P(N/X),\qquad M=N^d e .
$$
If $g=\gcd(e,P(r))$, then $B=s^d g$ divides both $M$ and $F(s)$. This holds with repeated prime factors in $e$ and without knowing them. For balanced primes, $|P(r)|=N^{d/2+o(1)}$. If $|P(r)|/g=N^{o(1)}$, then $B=N^{d-o(1)}$. With $e=N^{\alpha+o(1)}$, the standard degree-$d$ unknown-divisor Coppersmith threshold $M^{(\log_M B)^2/d}$ gives radius $N^{d/(d+\alpha)-o(1)}$. This exceeds $\sqrt{2N}$ by a fixed power for $\alpha\le d-\epsilon$. Public dyadic guesses of the lower bound on $B$ suffice, and every candidate is checked against $N$. This is a direct application of the published small-root theorem to the partner prime. [Cohn and Heninger, Theorem 1.1 and historical references](https://arxiv.org/pdf/1008.1284#page=2)

The order $P(r)=r^2+r+1$ of the degree-three finite-field norm-one group illustrates how the hint can be larger than $N$. In a 69-bit planted control with $N<e<N^2$, a degree-two reciprocal polynomial and a ten-dimensional lattice returned the partner prime 22,333,829,989. The missing cofactor was 13, so earlier group-order padding also applies to this particular input. Eighteen exhaustive windows tested the reciprocal implementation, including four with several roots. [Runnable reciprocal-polynomial check](check_reciprocal_order_polynomial.py) · [Earlier extension-field orders](translated-batches.md)

## Can the exponent come from the size of the input?

Test $e=\lfloor N^\alpha\rfloor$ for fixed $1/2\le\alpha<1$. Let distinct primes be chosen uniformly in $[R,2R]$. Holding $r$ fixed, $s\mapsto(rs)^\alpha$ has derivative of order $R^{2\alpha-1}$, bounded away from zero. Its rounded values occupy an integer interval of width $O(R^{2\alpha})$, and any value has $O(1)$ preimages. For either sign, $\gcd(e,r\pm1)\ge N^\gamma$ requires that $e$ be divisible by a divisor of $r\pm1$ of size at least $R^{2\gamma}$. Each such divisor allows at most $O(1+R^{2\alpha-2\gamma})$ partner integers. There are $R^{o(1)}$ divisors and $\Theta(R/\log R)$ prime partners. Thus the fraction meeting the order criterion at either prime is at most
$$
 R^{-1+o(1)}+R^{2\alpha-2\gamma-1+o(1)}.
$$
For $\gamma>\sqrt\alpha/2$ this tends to zero as a fixed negative power of $R$, since $2\alpha-1<\sqrt\alpha$ for $1/2\le\alpha<1$. This is an average bound for the explicit exponent rule. A different function of the input digits need not have this distribution.

A seeded test used 256 distinct prime pairs in each of four cohorts, at 16, 24, 32 and 40 bits per prime. For $\alpha=1/2,2/3,3/4$ the respective tested $\gamma$ were $9/25,5/12,7/16$. The square-root rule met its threshold twice across the 1,024 inputs. The other rules never met theirs. The factors were withheld from the exponent calculation and used to score the promises. In a separate check, 546 exhaustive small cases verified the count of admissible partners. [Unplanted benchmark and finite bounds](check_root_scale_order_hints.py)

To test end-to-end recovery, use a four-dimensional weighted lattice on a partition of the balanced interval into about one piece per input bit and try both signs. This choice was made after inspecting the two initial square-root hits. The public rule recovered both factors on the original 1,024 inputs. We then froze it and ran it on 1,024 fresh prime pairs with no overlap. The table counts verified splits from $N$ alone.

| Bits per prime | Original 256 pairs per row | Held-out 256 pairs per row |
| ---: | ---: | ---: |
| 16 | 1 | 1 |
| 24 | 1 | 0 |
| 32 | 0 | 0 |
| 40 | 0 | 0 |

The held-out split was $2{,}144{,}751{,}821=39{,}829\cdot53{,}849$. Its public exponent was $46{,}311=\lfloor\sqrt N\rfloor$. It shared 1,077 with $53{,}849+1$, above the public bound 922 for that input. In total the public rule made 113,775 weighted-window calls on the original set and 113,856 on the held-out set. The full runs grade the inputs using their factors only after the solver returns. [Replay both cohorts](check_root_scale_order_hints.py)

The determinant gives a bound for an actual decoding rule as well. Suppose all its search radii are $N^{1/2+o(1)}$, as when the balanced interval is partitioned into a number of equally sized pieces polynomial in the bit length. A passing certificate necessarily has
$$
 \log X < (\log L)^2/\log N+(\log G)^2/\log e .
$$
To see this, drop the positive LLL penalties from the certificate. For each multiplicity $m$ between zero and $D-1$, use $D(D-1)z^2-2Dmz+m(m+1)\ge0$ with $z=\log L/\log N$, and apply it again for $e,G$. The inequality follows by completing the square, or from its nonpositive discriminant. Therefore $G\ge N^{\sqrt\alpha/2-o(1)}$. The weighted decoder returns only divisors passing its chosen $G$. For $e=\lfloor\sqrt N\rfloor$, the counting argument bounds the mean success of this rule on balanced prime pairs by $N^{-1/(2\sqrt2)+o(1)}$. Its preliminary direct gcd costs no extra exponent in that bound. Conditional on one prime $r$, it requires a multiple of $r$ among $O(R)$ consecutive values of $e$, each attained by $O(1)$ partners. An input-dependent exponent with a different distribution remains a separate question.
