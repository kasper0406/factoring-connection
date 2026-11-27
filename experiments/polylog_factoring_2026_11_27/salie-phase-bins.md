# Factoring from one signed phase interval

For balanced squarefree semiprimes, one fixed phase interval per scale suffices for the Salié factoring reduction. The query is a signed integer sum with no oscillatory weights. RMS error of at most $\sqrt N/8$ on uniform unit scales suffices for randomized polynomial-time factoring. A separate constant-bin approximation recovers the original complex Salié target in mean square. Both statements concern the evaluation target. No polynomial-time interval evaluator is known here. [Phase coefficients](salie-phase-weights.md) · [Normalized Salié target](salie-quarter.md)

## The interval and its hidden frequency

For a unit $c$, let
$$W_c(t)=\sum_{m\in U_N:\ m+c^2/(4m)=t}(m/N).$$
Let $I$ be the indicator of the following public interval.

* For $N\equiv1\pmod4$, use the centered interval $[-\lfloor N/4\rfloor,\lfloor N/4\rfloor]$ modulo $N$.
* For $N\equiv3\pmod4$, use $[0,(N-3)/2]$.

The query is the integer $A(c)=\sum_t I(t)W_c(t)$. The second interval avoids cancellation when the Jacobi character is odd.

To analyze the Fourier transform in $c$, extend the function from units by the four-translate formula $W_c(t)=\sum_{r\in\{1,-1,u,-u\}}\chi(t+cr)$, where $u$ is a mixed root. This is a proof device on nonunit scales. For the normalized extension $A(c)/\sqrt N$, the coefficient at frequency $u$ is
$$
\frac{G_N}{\sqrt N}\left(
\hat I(1)+\chi(-1)\hat I(-1)
+\chi(u)\hat I(u)+\chi(-u)\hat I(-u)
\right).
$$
Indeed the Jacobi transform at $k$ is $G_N\chi(-k)/N$, and multiplication by any $r$ permutes the frequencies. The first two terms have magnitude $2/\pi+O(1/N)$. The centered interval gives the real part when the character is even. The interval beginning at zero gives the imaginary part when it is odd. The other two terms are $O(N^{-1/2})$: the least absolute representative of a mixed root has square at least $N+1$, and the interval Fourier coefficients are $O(1/|k|_N)$.

The normalized extension is uniformly $O(\log N)$ by the absolute Fourier sum of an interval and the Jacobi Gauss bound. Replacing its nonunit values by zero changes the coefficient by $O((\log N)/\sqrt N)$ on balanced semiprimes. An RMS error of $1/8$ on the normalized integer query still leaves a constant-size coefficient at $u$. Clip the approximated function at the public $O(\log N)$ bound and divide by that bound for the sampled Fourier recovery. Its heavy coefficient then has inverse-logarithmic size. Recover the frequency and take its difference from 1 by gcd. [Sampled recovery](quarter-fourier.md)

## Constant-bin approximation to the full phase

Alternatively replace $e_N(t)$ by $w_L(t)=\exp(2\pi i\lfloor Lt/N\rfloor/L)$ using $L$ bins. If $E=w_L-e_N$, then $\|\widehat E\|_2=\|E\|_2\le2\pi/L$. For unit $c$, the normalized Fourier transform of $W_c$ in its phase variable at unit frequency $k$ is $\chi(k)\sum_{r\in\{\pm1,\pm u\}}e_N(krc)$, in the convention with positive additive phase. At nonunit frequencies it is zero, since the primitive local character sums to zero at a prime where the phase vanishes. Expanding the bin error using these transforms, each fixed $r$ permutes the frequencies and contributes at most $\|\widehat E\|_2$ to the RMS over all residue scales, using the displayed formula to extend the error. Thus the RMS error over unit scales is at most $8\pi\sqrt{N/\varphi(N)}/L$. A sufficiently large absolute constant $L$ suffices for the normalized approximation target. This estimate uses cancellation across scales. A uniform-in-scale bound instead follows from the Fourier absolute sum and costs $O(\log N)$ bins.

## Checks

The diagnostic verified the single-interval formula against grouped public terms at 96 unit scales on four semiprimes. On all four inputs the full extension's coefficient at $u$ had magnitude between 0.629 and 0.640. It also checked finite phase bin approximations at 16, 64 and 256 bins. Run `python check_salie_phase_bins.py` with NumPy and SymPy.

## Constant clipping and sampled recovery

The normalized four-translate extension $B(c)=A(c)/\sqrt N$ has $L^2$ norm at most 2. Each translate separately has squared norm at most $\|I-\mathbb E I\|_2^2\le1/4$: its Fourier coefficients are those of the interval indicator, multiplied by unit-magnitude Gauss factors at unit frequencies and zero elsewhere. Then apply the triangle inequality across the four translates. Zeroing the nonunit inputs preserves the bound.

Clip at absolute value 4. The change in any Fourier coefficient is at most $\mathbb E(|B|-4)_+\le\mathbb E B^2/16\le1/4$, using $(x-a)_+\le x^2/(4a)$ for $x\ge0$. Clipping is Lipschitz, so the promised normalized RMS oracle error costs at most a further $1/8$. After division by 4, a constant-size coefficient at the hidden frequency survives on sufficiently large balanced semiprimes. This removes the logarithmic normalization from the first decoder.

The noisy-interval diagnostic uses the factor-supplied array as an oracle, adds independent bounded integer errors, clips at 4, and runs the sampled frequency search with the direct gcd probes disabled. Run `python check_salie_interval_recovery.py`.
