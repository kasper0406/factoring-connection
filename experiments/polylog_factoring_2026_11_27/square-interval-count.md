# Factoring from an unweighted inverse phase interval

For an odd input with exactly two distinct prime divisors, an approximate unweighted interval count modulo $N^2$ suffices for randomized polynomial-time factoring. The interval is fixed, the lower digit is fixed, and only the public upper digit varies. The proof uses exact counting in arithmetic progressions and the sampled Fourier decoder. It needs no character weights or prime-modulus exponential-sum bounds. [Square-modulus phase identity](square-kloosterman.md) · [Sampled decoder](quarter-fourier.md)

Let $M=N^2$, let $I$ be the centered interval $[-\lfloor M/4\rfloor,\lfloor M/4\rfloor]$ modulo $M$, and put $c=1+Nj$ for $j\bmod N$. Ask for the count
$$C(j)=\#\{x\in(\mathbb Z/M)^\times:x+c^2/x\in I\pmod M\}.$$

## Separate the blocks by their exact periods

Write $x=x_0+Ny$ for public residue classes $x_0\bmod N$ and $y\bmod N$. The phase is affine in $y$ modulo $M$:
$$x_0+x_0^{-1}+2Nj,x_0^{-1}+Ny(1-x_0^{-2}).$$
Let $d=\gcd(x_0^2-1,N)$. As $y$ varies, the phase runs over one residue class modulo $Nd$, each point repeated $d$ times. This class depends on $j$ only modulo $d$. Its contribution to the interval count therefore has period dividing $d$. Between two shifts its count can change by at most $d$, because the number of integers from a contiguous circular interval in any such class differs by at most one.

When $d<N$, the block has zero Fourier coefficient at all unit frequencies in $j$. This follows from its exact period. When $d=N$, $x_0$ is a root of 1. Let $U$ be its unique lift with $U^2=1\pmod{N^2}$. Its block contributes exactly $N I(2cU)$, evaluated at that lift. The public roots $\pm1$ contribute $2NI(2c)$ and can be subtracted directly.

## The mixed pair supplies a large coefficient

Let $u=U\bmod N$ be a mixed root. The remaining two stationary blocks have the same membership bit by reflection. After normalization by $N$, they contribute $2I(2U+2Nju)$. Within the residue class of $2U\bmod N$, the fixed interval selects a circular interval in the upper digit with length $N/2+O(1)$. Its first normalized Fourier coefficient has magnitude $1/\pi+O(1/N)$ by a geometric sum. Thus the count after subtracting the public blocks, divided by $N$, has a coefficient of magnitude $2/\pi+O(1/N)$ at the unit frequency $2u$. No other blocks contribute at that frequency.

Subtract the corrected count at one reference shift before normalization. This does not change the nonzero coefficient, and its output is bounded polynomially in the bit length. In fact the variation divided by $N$ is at most
$$\frac1N\sum_{x_0\in U_N}\gcd(x_0^2-1,N)
=\frac{\varphi(N)}N\prod_{r^e\parallel N}(2e+1).$$
The identity follows from $d=\sum_{a\mid d}\varphi(a)$. At a prime power $r^e$, the number of residues satisfying $x_0^2=1\pmod{r^k}$ is $2r^{e-k}$ for $1\le k\le e$. For two distinct primes the bound is $O((\log N)^2)$, and for a squarefree semiprime it is at most 9.

If a public evaluator estimates the corrected count difference with RMS error at most $N/8$ over uniform $j$, its normalized signal still has a constant-size coefficient at $2u$. Clip it to the public polynomial bound and divide by that bound. The sampled Fourier recovery finds a coefficient of inverse-polynomial size. Test the recovered frequencies minus 2 by gcd to obtain a verified split. The same guarantee follows from sufficiently accurate uniform pointwise estimates of the raw count and one shared reference count.

A direct scan can compute $C(j)$ with $N$ block queries by clipping each arithmetic progression with floors. The grouping explains the cancellation, but the scan still has a power-of-$N$ cost. A polynomial-bit-time evaluator of the count difference remains open.

## Checks

Direct enumeration modulo $N^2$ verified all 418 upper digits across four inputs including $N=225$. It checked that removing the mixed blocks leaves zero at every unit frequency, and that the hidden coefficient gives a proper gcd. Run `python check_square_interval_count.py`.
