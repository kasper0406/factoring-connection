# Fourier dual and autocorrelation of the quadratic-scale Salié sum

Averaging the product of two normalized Salié values over all scales returns the original signal at the shift, multiplied by a scalar near one on balanced semiprimes. Completing the quadratic phase removes one summation variable exactly. The resulting formula retains the original Salié sum. [Normalized analytic target](salie-quarter.md)

Let $T(c)=S_N(1,c^2/4)/G_N$ for odd $N$, including nonunit scales $c$. Use the normalized additive Fourier transform. Expanding the Salié sum, its coefficient at $k$ is
$$\widehat T(k)=\frac1{NG_N}\sum_{m\in U_N}\chi(m)e_N(m)
\sum_{c\bmod N}e_N(c^2/(4m)-kc).$$
Completing the square in $c$ makes the inner sum $G_N\chi(m)e_N(-mk^2)$. The two Jacobi factors cancel, giving
$$\widehat T(k)=\frac{c_N(1-k^2)}N,$$
where $c_N$ is the ordinary Ramanujan sum.

For a squarefree semiprime, most dual coefficients equal $1/N$. Any frequency with a proper gcd of $1-k^2$ already exposes a factor. If $k^2=1\pmod N$, then either $k$ is public $\pm1$ or a gcd of $k-1$ gives a factor. At the public frequencies the coefficient is $\varphi(N)/N$, close to one on balanced inputs. Thus generic pointwise access to this dual table is cheap. Locating the hidden peaks is the remaining search.

## Exact autocorrelation

Define $A(h)=\mathbb E_{c\bmod N}T(c)\overline{T(c+h)}$. For a unit shift $h$, expand this as a double sum with unit indices $m,l$ and put $d=l-m$. The quadratic coefficient in $c$ is $d/(4ml)$ and the linear coefficient is $-h/(2l)$. When $d$ is a nonunit, translation by $N/r$ for a prime divisor $r$ of $\gcd(d,N)$ leaves the quadratic term invariant and changes the linear phase by a nontrivial root of unity. Those pairs contribute zero.

For unit $d$, completing the square gives residual phase $-d-h^2/(4d)$ and combined character $\chi(d)$. Both are independent of $m$. There are
$$Q_N=N\prod_{r\mid N}(1-2/r)$$
admissible $m$ for every unit $d$, since $m$ and $m+d$ must avoid two distinct residues at each prime. The remaining $d$ sum is the complex conjugate of the original Salié sum at scale $h$. Since $|G_N|^2=N$ and $T$ is real, the result is
$$A(h)=\frac{Q_N}N T(h).$$
On balanced semiprimes, $Q_N/N=1-O(N^{-1/2})$ and $|T(h)|\le4$ for unit $h$. Thus a constant-precision evaluator for this averaged correlation, with sufficient slack in its error, would also meet the earlier factoring criterion. This identity does not supply such an evaluator.

At $h=0$ on a squarefree semiprime, Parseval instead gives
$$A(0)=(2-3/p)(2-3/q)=4-\frac{6(p+q)-9}{N}.$$
Indeed at one prime $r$, the squared Ramanujan coefficient gives two terms $(r-1)^2$ and $r-2$ terms of value 1 before division by $r^2$. Therefore an exact energy, or an approximation with error less than $3/N$, determines $p+q$ by rounding. Constant accuracy for this zero-shift energy alone gives no such rounding.

## Checks

The diagnostic generated the Salié table by direct finite Fourier transforms on five moduli, including two with repeated prime factors. It verified every dual coefficient and every unit-shift autocorrelation. The three squarefree controls also passed the zero-shift energy formula. Run `python check_salie_autocorrelation.py` with NumPy and SymPy.

## Correcting exceptional shifts by a three-query identity

For distinct odd primes N=pq put C(x)=cos(2 pi ux/N), where u
is either mixed root. For any target h and mask r one has exactly

    C(h) = 2 C(r) C(h-r) - C(2r-h).

This follows from the cosine product identity before averaging.
With r uniform modulo N, all three query arguments are individually
uniform, since 2 is invertible. No independence among the three
arguments is required. On nonunits use the public cosine directly,
since ux equals x or -x there.

For example suppose a normalized Salié estimate on a uniform unit
is within 0.038 with probability at least 0.9 (including fresh
internal randomness). Subtract the public cosine pair, halve, and
clip to [-1,1] to estimate C. Each masked cosine query now has
error <=0.019 except with probability <=0.1. With probability
>=0.7 all three are good. The resulting estimate of C(h) has
error <=5*0.019=0.095: the product of two clipped estimates has
error <=2*0.019, then double it and add the third error. This is
inside the direct phase decoder's 0.1 cosine tolerance.

Use independent masks and fresh oracle coins across trials. A median
of O(log(Q/delta)) trials succeeds with failure at most delta/Q
by Hoeffding, where Q=O(log N) is a bound on the number of decoder
queries. Fresh trials after each adaptive query allow a union bound.
One can therefore use the deterministic phase backtracking with
overall success at least 1-delta on this weaker coverage promise.
This is a local correction procedure given the average-case
estimator; it does not construct that estimator. The separate
heavy-Fourier recovery in salie-quarter.md tolerates a looser RMS
promise when only factoring, rather than answering a chosen phase
query, is required.
