# A weak quarter-window predictor is enough to factor

For an odd $N$ with exactly two distinct prime divisors, a predictor for the scaled quarter-window bit with advantage $\delta$ on uniformly random unit scales gives a randomized factoring algorithm with bit cost polynomial in $\log N$, $1/\delta$, and the cost of a prediction. The predictor may have input-dependent errors. It can also be randomized, with fresh randomness on every call. Repeated prime factors are allowed. [Exact quarter-window reduction](nonadaptive-quarter.md)

## Correlation exposes a hidden frequency

Let $u$ be a mixed root of 1. For a unit scale $c$, the quarter-window bit for the hidden pair $\pm cu$ is the sign
$$g(c)=g_0(cu\bmod N),\qquad g_0(x)=2[\min(x,N-x)\le\lfloor N/4\rfloor]-1.$$
A prediction of the signed disk bit gives a prediction of $g$ by removing the axis and the known pair $\pm c$. The premise applies after these public corrections, or directly to the hidden-pair bit.

Write $F(c)$ for the expected predicted sign at a unit $c$, and extend it by zero elsewhere. The success premise gives
$$\mathbb E_{c\bmod N} F(c)g(c)\ge 2\delta\,\varphi(N)/N\ge\delta.$$
The extension is evaluable by gcd. A proper gcd already splits $N$. The value at zero can be set to zero.

Use normalized discrete Fourier coefficients. For $h\ne0$, the geometric sum on the centered quarter interval gives $|\widehat g_0(h)|\le1/|h|_N$, where $|h|_N$ is the smaller absolute representative. The constant coefficient has magnitude $1/N$. Put $K=\lceil32/\delta^2\rceil$. If $N\le2K+1$, trial division is polynomial in $1/\delta$. Otherwise the Fourier tail beyond $K$ has squared norm at most $2/K$. By Parseval and Cauchy–Schwarz its contribution to the correlation is at most $\delta/4$. The constant contribution is also at most $\delta/4$. Since $\widehat g(hu)=\widehat g_0(h)$, some nonzero frequency $k=hu\bmod N$ with $|h|\le K$ has
$$|\widehat F(k)|\ge \frac{\delta}{4(1+\log K)}.$$

Recover all frequencies above this threshold by the interval search below. For each frequency $k$ try every nonzero $h$ in the prescribed range. If $h$ is a nonunit, its gcd can already give a factor. For the correct pair with unit $h$, $\gcd(k-h,N)=\gcd(u-1,N)$ is a proper divisor. All candidate divisors are verified.

## Find heavy frequencies with public interval estimates

For an integer frequency interval $I$, put $D_I(t)=\sum_{k\in I}\exp(2\pi i kt/N)$. Its energy is
$$\sum_{k\in I}|\widehat F(k)|^2
=\mathbb E_{x,t\bmod N} F(x)\overline{F(x+t)}D_I(t).$$
To estimate this, choose a distance shell for $t$ uniformly among $\{0\}$ and the nonempty shells $2^j\le |t|_N<2^{j+1}$, then sample uniformly in that shell. For its probability $q(t)$, multiply each sample by $D_I(t)/(Nq(t))$. This has magnitude at most the number of shells. For $t\ne0$, use $|D_I(t)|\le N/(2|t|_N)$ and at most $2^{j+1}$ points in the shell. When the predictor is randomized, take independent calls for the two factors, including when $t=0$. Compute the kernel by the geometric sum formula to additive precision below the required error in the bounded sample weight. Modular reduction of the phase arguments and polynomial-bit trigonometric approximation suffice.

For threshold $\tau$, estimate each interval energy to error $\tau^2/4$ with high confidence by bounded sampling. Start with the whole frequency range, bisect, and retain an interval if its estimated energy is at least $\tau^2/2$. Every interval containing a frequency of magnitude at least $\tau$ survives. On accurate estimates, each retained interval has energy at least $\tau^2/4$, so Parseval bounds their number at each level by $4/\tau^2$. Cap the work at that level bound, and choose confidence per estimate by a union bound over the resulting polynomial number of tests. This gives bounded-error recovery in time polynomial in $\log N$ and $1/\tau$.

## Checks and remaining target

The diagnostic planted noisy quarter-window predictors on four semiprimes. The interval search using exact spectral energies recovered factoring frequencies on all four. The shell estimator separately agreed with direct energies on 12 tested intervals. These tests use fixture roots to build the planted predictors and exact Fourier transforms as controls. They do not provide a predictor on unfactored input. Run `python check_quarter_fourier.py`.

A candidate evaluator can now be tested for advantage on independently sampled unit scales at a fixed modulus. An inverse-polynomial advantage in bit length would suffice, without requiring correct answers at all scales or finding a reliable prefix of the hidden root directly.

## Sampled end-to-end recovery

The executable search in `quarter_fourier_recovery.py` builds interval energies only from sampled predictor calls. It includes a conservative sample budget calculated from the advantage and a separate calibration with a stronger planted signal. In the calibration, flip each planted answer independently with probability 12%, disable the direct gcd probes on sample inputs, and search with threshold $1/4$ and 8,192 samples per interval. The search recovered a fundamental hidden frequency and a verified factor on all five inputs, ranging from 8 to 36 bits and including one repeated-prime input. It used about 421,000 to 2.31 million predictor calls per input. These calibration parameters demonstrate recovery for the planted signal. The general advantage guarantee uses the larger bound from the proof. Neither run supplies a predictor on unfactored input.

## Reuse predictor calls across the interval tree

The latest implementation draws one pool of pairs $(x,x+t)$ from the shell distribution. For each pair it stores the sampled prediction product times its importance weight, aggregated by $t$. Each interval applies a different Dirichlet kernel to those same weighted products. No interval evaluation needs a fresh predictor call.

For a fixed interval the reusable estimate has the same concentration bound as before. To handle the adaptive choice of intervals, take a union bound over the entire predetermined binary partition tree, which has $2N-1$ nodes. The logarithm of this node count is only $O(\log N)$. For threshold $\tau$ and failure probability $\epsilon$, a single pool of $O((\log N)^2\tau^{-4}\log(N/\epsilon))$ pairs therefore suffices for all intervals at once. The interval search still computes kernels separately for each visited interval, so this improves predictor calls and adds storage for the pool.

On the five planted fixtures above, the shared pool of 16,384 pairs recovered a fundamental hidden frequency and a verified divisor on every input, with 25,507 to 32,768 predictor calls. The 8,192-pair shared pool missed the 36-bit fixture in the seeded calibration. The doubled pool compares with about 421,000 to 2.31 million predictor calls for the earlier per-interval estimator. These are calibration settings. The conservative sample bound uses the full tree size.
