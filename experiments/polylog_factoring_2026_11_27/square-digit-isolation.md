# Unweighted digit filtered counts for arbitrary squarefree moduli

The square-modulus counting reduction extends to every odd squarefree input. It uses two unweighted counts distinguished by a public bit filter on the lower phase digit. An evaluator of the signed count difference with additive RMS error $N/8$ gives randomized polynomial-time factoring. The number of distinct primes need not be known. This is an oracle reduction. A polynomial-time evaluator for the counts remains open. [Two-prime interval](square-interval-count.md) · [Affine isolation bound](affine-pair-bound.md)

## Two digit filters

Choose a uniform unit $c_0$ and uniform offset $b\bmod N$. For upper input digit $j$, write $c=c_0+Nj$ and decompose the phase of a unit $x\bmod N^2$ as
$$x+c^2/x=a+Nz\pmod{N^2},\qquad 0\le a,z<N.$$
Retain the phase if $(a+b)\bmod N<w$ and $z\le(N-1)/2$. Independently choose pairwise independent signs $\epsilon(a)$, for example by a random binary inner product with the bits of $a$, plus an independent constant bit. Count separately the retained phases with sign $+1$ and sign $-1$. Subtract the counts, remove the contributions of the public stationary roots $\pm c$, and subtract the same corrected quantity at reference $j=0$. Normalize by $N$ to obtain $F(j)$.

For a fixed block $x=x_0+Ny$, put $d=\gcd(x_0^2-c_0^2,N)$. Its lower phase digit is constant, and its upper digit visits one class modulo $d$, each point repeated $d$ times. As $j$ varies the class depends only on $j\bmod d$. Its normalized contribution has variation at most $d/N$, since an interval intersects two such classes in counts that differ by at most one before the multiplicity is applied. If $d<N$ the block has no unit Fourier frequencies in $j$.

For a stationary root $u^2=1\pmod N$, the block has lower phase $2c_0u\bmod N$, and upper phase $z_0+2uj\bmod N$ for a constant $z_0$. If the lower window contains exactly one mixed root, the normalized hidden block contributes one signed interval indicator in the upper phase. Its coefficient at $2u\bmod N$ has magnitude $1/\pi+O(1/N)$. The sign and the reference shift do not change that magnitude. The previous affine isolation proof applies to the lower phases $2c_0u\bmod N$, since $2c_0$ is uniform over units.

## Variance after random signs

Let $Q=2^{\omega(N)}$ and group blocks by lower phase $a$. Define
$$B_a=\frac1N\sum_{x_0:\ x_0+c_0^2/x_0=a}\gcd(x_0^2-c_0^2,N).$$
Within one group the change of its normalized count between two upper digits is at most $B_a$. Removing the public stationary terms preserves this bound, since it removes entire blocks. The pairwise independent signs eliminate cross terms when squaring. Thus, for fixed $c_0$, averaging over the signs and the translation gives
$$\mathbb E_{b,\epsilon}\|F\|_2^2\le\frac wN\sum_a B_a^2.$$

The sum factors by CRT. At a prime $p$, the two endpoint phases $\pm2c_0$ have one preimage each of local weight one. Every other phase with preimages has two, each of local weight $1/p$. The number of these other phases is $(p-3)/2$, by pairing the remaining units under $x_0\mapsto c_0^2/x_0$. Hence the local squared sum is $2+2(p-3)/p^2\le2(1+1/p)$. It follows that
$$\sum_a B_a^2\le Q\prod_{p\mid N}(1+1/p)\le Q(1+\omega(N)).$$
The last bound follows by ordering the primes and telescoping with $p_i\ge i+1$.

Use the dyadic width from the uniform isolation proof. When $N/(8192Q)\ge1$, this width satisfies $w/N\le1/(8192Q)$ and has singleton probability at least $1/(8\cdot8192)$. Markov's inequality shows that with probability at least half that lower bound, a trial is a singleton and $\|F\|_2^2\le16(n+1)$, where $n$ is the bit length. If $N/(8192Q)<1$, use width one. Its singleton probability is already constant from below in this case, and $B_a\le1$ by the local group formula. Trying all dyadic widths is public and incurs only $O(n)$ width choices.

## Decode a good trial

Clip the normalized estimate at $64(n+1)$. On a good trial the clipping changes any Fourier coefficient by at most $\|F\|_2^2/(4\cdot64(n+1))\le1/16$. An additional normalized RMS error $1/8$ still leaves a constant-size coefficient at $2u$ for sufficiently large inputs. Divide by the clipping bound and apply the sampled Fourier search. A gcd of frequency minus 2 with $N$ supplies a verified divisor. Finite small cases can use trial division. There is no need to identify which random trial was good. Recurse on the divisors under the squarefree promise.

The RMS premise is on the signed count contrast for each fixed trial over uniform upper digits. Alternatively, a uniform pointwise additive error of at most $N/32$ in each unweighted count suffices, since the contrast uses the two counts at the chosen digit and at the reference. The proof bounds query construction, clipping and decoding by polynomial bit time. It leaves the count evaluation open and uses squarefreeness in the local variance calculation.

## Checks

The diagnostic verified the local squared norm formula and the unit-frequency isolation on products of two through five primes. The largest input was 15,015. Direct summation modulo $N^2$ cross-checked four upper digits on each of the two smallest fixtures. Run `python check_square_digit_isolation.py` with NumPy.
