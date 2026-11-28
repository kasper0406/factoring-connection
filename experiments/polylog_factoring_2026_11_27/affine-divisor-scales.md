# One random center at several divisor sizes

The balanced two-prime window can also be used as an ordinary splitting trial without a factor-size promise. For every composite $N\ge4$ the implementation has success probability at least $N^{-1/4}/256$ on a uniform center. Every center costs polynomial bit time. This applies to prime powers, products with several primes, and very unbalanced inputs. It is a randomized use of the classical approximate-common-divisor bound. Repetition has a worse exponent than the best general deterministic factoring algorithms. [Code](affine_divisor_scales.py) · [Small-root bound and prior work](coupled-rank-windows.md) · [Harvey](https://arxiv.org/abs/2010.05450)

Let $B=\lceil\log_2N\rceil$ and $L=\log_2 N$. For a possible proper divisor $b\ge\sqrt N$, set $k=\lfloor\log_2 b\rfloor$ and $R=2^k$. The code tries every relevant $k$ on the same center $x$. Take $X=2^h$ if $h=\lfloor k(k-1)/(B-1)\rfloor-2$ is nonnegative, and $X=0$ otherwise. For each band, the reader returns every $\gcd(x+y,N)\ge R$ with $|y|\le X$. Any proper gcd is verified and returned. If the window already contains a multiple of $N$, every other offset has gcd at most $2X<R$ and this band can be skipped.

For $X>0$ choose affine lattice dimension $D=B$ and multiplicity $m=k$. Relative to a divisor lower bound $R$, the base-two logarithm of the LLL and evaluation bound is at most
$$\tfrac12\log_2 B+\tfrac{B-1}{4}+\tfrac{k(k+1)}2+\tfrac{B-1}{2}\log_2 X.$$
It is strictly less than $k^2=\log_2 R^m$ by the chosen two-bit margin. The search for the smallest sufficient dimension therefore terminates. For fixed dimension the bound is best at $m=\lfloor\log_N R^D\rfloor$ clamped to the allowed integer range, as seen by comparing adjacent multiplicities. The solver's exact integer certificate covers every offset, including when $b$ is composite or a prime power.

Write $u=\log_2 b$ and $\beta=u/L$. Since $k\le u<k+1$, $B-1<L\le B$, and $b\le N$, we have $u^2/L-k^2/B\le3$. Also $k^2/B-k(k-1)/(B-1)=k(B-k)/(B(B-1))\le1/2$. Accounting for the floor and the two-bit margin gives
$$W=2X+1\ge N^{\beta^2}/128.$$
If $h<0$, the same inequality follows from $k(k-1)/(B-1)<2$.

To bound success for any composite, let $b=N/p_{\min}$, so $\sqrt N\le b\le N/2$. The $W$ offsets give distinct residues modulo $b$. With probability $W/b$ exactly one of them produces a multiple of $b$. Conditional on that event it is a multiple of $N$ with probability $b/N$, by uniformity on the $N/b$ lifts of each residue. Every other outcome on this event yields a proper divisor. Thus one band alone succeeds with probability at least
$$\frac Wb(1-b/N)\ge\frac{1}{256}N^{-\beta(1-\beta)}\ge\frac{1}{256}N^{-1/4}.$$
The algorithm does not require $b$ or its band as input. The probabilities at different bands need no independence assumption. The tests enumerate every center on 14 small moduli, validate the sufficient LLL inequalities on larger bands and recover verified divisors on unbalanced and repeated-factor examples. [Reproducible checks](check_affine_divisor_scales.py)
