# One bit from a fixed dimensional count suffices for factoring

The modular interval counting reduction can use a one-bit gap comparison instead of a numeric approximation. For each public trial and upper digit, let $D(j)$ be the corrected signed count difference before division by $N$. The oracle returns zero when $|D(j)|\le N/4$ and one when $|D(j)|\ge3N/4$, with either answer permitted between the two. This suffices for randomized polynomial-time factoring, including repeated prime factors. [Affine interval count](affine-interval-variance.md) · [Sampled Fourier decoder](quarter-fourier.md)

## Make the nonstationary variation small

Use the lower input isolation window and the affine interval signs. The contrast subtracts the two public stationary blocks and the reference at upper digit zero. In a singleton trial write its normalized value as
$$D(j)/N=H(j)+R(j),$$
where the one remaining stationary block contributes $H(j)=\epsilon(I(z_0+2uj)-I(z_0))\in\{-1,0,1\}$. The term $R$ collects only the nonstationary blocks.

Before sampling, try trial division up to $P=2^{18}C(n+1)$ with $C=8192$ and $n=N.\mathrm{bit\_length}()$. Any factor found suffices for a split. In the remaining case every prime exceeds $P$. Let $Q=2^{\omega(N)}$. In the affine-sign variance proof, the global pair kernel is at most $Q\prod_{p\mid N}(1+5/p)$. Removing the stationary inputs removes at least the $Q$ diagonal terms of weight one. The kernel is nonnegative term by term. Therefore at the isolation width, where $wQ/N\le1$,
$$\mathbb E\|R\|_2^2
\le \prod_{p\mid N}(1+5/p)-1
\le 10\omega(N)/P.$$
This estimate concerns the nonstationary remainder even on trials with several stationary blocks. The last inequality uses $\exp(x)-1\le2x$ for $0\le x\le1$.

The dyadic width has singleton probability at least $1/(8C)$. Markov's inequality makes $\|R\|_2^2>1/512$ occur with probability below $1/(16C)$, so singleton isolation and the small remainder occur together with constant probability. The public reduction tries all the allowed widths, since $Q$ is unknown.

## Decode the comparison bit

On a good trial, at most $1/32$ of upper digits have $|R(j)|>1/4$. On the other digits the gap comparison returns the indicator of $H(j)\ne0$. This indicator is the stationary half-interval bit or its complement, depending on the reference value. Its Fourier coefficient at $2u$ has magnitude $1/\pi+O(1/N)$. An arbitrary error on at most $1/32$ of the digits changes the coefficient by at most $1/32$. Thus the oracle bit has a heavy coefficient of absolute constant size. The sampled Fourier algorithm finds it without clipping or magnitude estimation. Testing the recovered frequency minus 2 by gcd returns a proper divisor. Repeat trials until a verified split, and use primality and perfect-power tests for the leaves.

This is a decision target on a contrast of ordinary fixed-dimensional quadratic counts. The oracle need not supply the counts themselves. No polynomial-time comparison algorithm is established.

## Checks

The diagnostic obtained the comparator from exact block counts on five moduli, three with repeated factors. It verified the hidden coefficient and a proper gcd on each. The nonstationary RMS on these fixtures ranged from zero to 0.055. These controls use smaller factors than the conservative asymptotic preprocessing threshold. Run `python check_count_comparison.py`.
