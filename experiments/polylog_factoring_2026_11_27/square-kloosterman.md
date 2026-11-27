# Factoring from square modulus Kloosterman approximation

For an odd input with exactly two distinct prime divisors, constant-precision approximation of an ordinary Kloosterman sum modulo $N^2$ suffices for randomized polynomial-time factoring. Only its upper base-$N$ digit varies. This gives an alternative to the Salié target with no Jacobi twist. [Salié approximation](salie-quarter.md)

Let
$$K_M(1,b)=\sum_{x\in(\mathbb Z/M)^\times} e_M(x+b/x),\qquad e_M(t)=e^{2\pi it/M}.$$
For $M=N^2$ and unit $c$, write $x=x_0+Ny$. Modulo $N^2$,
$$x+c^2/x=x_0+c^2/x_0+Ny(1-c^2/x_0^2).$$
Summing over $y\bmod N$ gives zero unless $x_0^2=c^2\pmod N$, in which case it gives $N$ times the first phase. Each stationary $x_0$ has a unique lift $r$ with $r^2=c^2\pmod{N^2}$, by the unit derivative $2x_0$. At the lifted value the phase is $2r$. Consequently
$$\frac{K_{N^2}(1,c^2)}N=\sum_{r^2=c^2\bmod N^2}e_{N^2}(2r).$$

## Vary just one digit

Let $u$ be a mixed root modulo $N$ and $U$ its lift modulo $N^2$. Choose $c=1+Nj$ with $j\bmod N$. The public root pair $\pm c$ contributes $2\cos(4\pi c/N^2)$. The remaining pair contributes
$$2\cos(4\pi U/N^2+4\pi uj/N).$$
Thus the unknown component is a sinusoid in $j$ with frequency $2u\bmod N$, even though all the scales have the same lower digit.

Suppose a public algorithm approximates $\operatorname{Re}K_{N^2}(1,(1+Nj)^2)/N$ to RMS error at most $1/2$ over uniform $j$. Subtract the public term, divide by two and clip to $[-1,1]$. The resulting signal has RMS error at most $1/4$ from a unit-amplitude cosine. Its correlation with that cosine is at least $1/4$ by Cauchy–Schwarz and the exact mean square $1/2$ of the cosine. This forces a Fourier coefficient at $2u$ of magnitude at least $1/4$ despite the unknown phase. The public sampled Fourier search finds the frequency. Testing its difference from 2 by gcd gives a proper divisor. The sum identity and lifting also apply when $N$ has repeated factors.

## What the block cancellation saves

The lift sum reduces $N^2$ terms to a scan over at most $N$ residue blocks, with a polynomial-bit-time stationarity test per block. Scanning only for the unknown stationary residues is exactly a search for the mixed roots modulo $N$. Once the two public stationary residues have been removed, a uniform block hits either unknown block with probability $2/N$. Independent uniform block sampling therefore does not provide the needed precision efficiently. An evaluator that aggregates or locates stationary blocks without this scan remains open.

## Checks

Direct Kloosterman enumeration verified the identity and the two upper-digit frequencies on all 179 upper digits across five small moduli, including two repeated-factor inputs. The diagnostic also verified a proper gcd from both recovered frequencies. Run `python check_square_kloosterman.py`.
