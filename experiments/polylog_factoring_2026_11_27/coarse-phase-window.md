# Native averages of cubic and quartic magnitudes

Let $N=PQ$ have exactly two distinct odd prime-power components, put $M=N^2$, and take the fixed cubic $3x-x^3$ or quartic $x^4-2x^2$. For the cubic assume $3\nmid N$. Average the normalized complete-sum magnitude, or its square, over the unit scales in an integer interval $[a,a+L)$. For $L\ge N/\operatorname{poly}(\log N)$ and inverse-polynomial additive error, the window average has a native polynomial-bit algorithm from $N,L,a$. It may discover a factor when a hidden frequency is unusually low. On other inputs the returned approximation is a public diagonal contribution. Single-scale magnitudes remain open.

## Reduce low frequencies to ordinary quadratics

At unit scales, stationary phase gives a product of two local expressions. For the squared magnitudes the local Laurent coefficients are
$$\begin{aligned}
W_3(z)&=2+z^4+z^{-4},\\
W_4(z)&=5+2z+2z^{-1}.
\end{aligned}$$
Evaluate at $z=e_M(e_Pt)$ and $z=e_M(e_Qt)$, respectively, for the two CRT idempotents. Use this product to extend the stationary expression to all integer scales. For any finite local list of labels with weights $w_j$, the extension has total absolute Fourier mass at most $W=(\sum_j |w_j|)^2$. Each line has a frequency $f$ which is a local label $j$ modulo $P^2$ and a local label $k$ modulo $Q^2$. If $j=k$, the line is public with coefficient $w_j^2$.

For $j\ne k$ the centered integer frequency satisfies
$$ (f-j)(f-k)=uM.$$
Trial-divide up to a bound that exceeds all label differences and is large enough for the desired unit-density error. If this finds a prime, strip its full power and use both components. Otherwise the bases of $P,Q$ exceed every label difference. Then every mixed line gives a proper component by $\gcd(f-j,N)$. To find all possible mixed lines of absolute frequency at most $B$, try each pair of labels and $0\le u\le (B+\max(|j|,|k|))^2/M$. Test whether $(j-k)^2+4uM$ is a square, recover the two integer candidates, and verify a proper gcd. Negative $u$ cannot occur after the trial bound since the quadratic minimum is $-(j-k)^2/4>-M$. This enumeration includes every true line and verifies every returned divisor. It uses only integer square roots and gcds.

If a component is found, compute the full list of frequencies by CRT. Average each by the geometric sum formula, including and excluding multiples of the two prime bases to keep exactly the unit scales. Exact integer root tests and primality suffice to obtain those bases from the verified components. A geometric sum on any arithmetic progression has a closed sine ratio. When the denominator is nonzero it has magnitude at least $2/M$, so polynomially many working bits suffice for the requested numerical accuracy. The diagnostic uses high-precision arithmetic for these kernels.

## Error and total cost

Write $r=\lceil1/\varepsilon\rceil$. A line of centered nonzero frequency $f$ has an interval average of magnitude at most $M/(2L|f|)$. Choose $B=\min(\lfloor M/2\rfloor,\lceil WrM/L\rceil)$. If the search finds no component, it certifies that all mixed lines lie outside the cutoff. Their omitted contribution is at most $1/(2r)$. If the interval is shorter than the trial-division bound, use $B=\lfloor M/2\rfloor$ and recover the components from any mixed pair. This short-window branch is charged in the bound below.

In the unfactored branch both primes and $L$ exceed $32Wr$ up to endpoint rounding. The proportion of nonunit scales is at most $1/p+1/q+2/L$, where $p,q$ are the two prime bases. As the extended Laurent statistic lies in $[-W,W]$, the difference between all-scale and unit-only averaging is at most $1/(4r)$. The public diagonal average then has error below $1/r$ with numerical error below $1/(4r)$.

For squared magnitude there are three local labels. The quotient bound costs $O(r^2(1+N/L)^2)$ iterations, with polynomial bit overhead. The forced full cutoff at shorter lengths still falls under this bound because $L=O(r)$ there. All trial divisions, label-pair construction, and geometric sums are charged.

For the unsquared cubic magnitude use the local series $2|\cos(4\pi x)|$. Its constant coefficient is $4/\pi$ and its coefficients at $\pm4h$ are $(4/\pi)(-1)^{h+1}/(4h^2-1)$. Its total absolute weight is below three and its discarded weight beyond $H$ is $(4/\pi)/(2H+1)$. Taking $H=12r$ changes the product of two factors by less than $1/(2r)$. For the quartic magnitude use the local formula $2|b(z)|^2$, where $b(z)=(1+z/2)^{1/2}$. The coefficients of $b$ have absolute value at most $2^{-j}$. Truncate at $H=\lceil\log_2(512r)\rceil$. The resulting product of two local expressions has uniform error below $1/(2r)$ and total absolute Fourier weight at most 64. Apply the finite-list algorithm at precision $2r$ to either truncation. Including all pairs of labels, the total work is bounded by $O(r^4(1+N/L)^2)$ times polynomial bit overhead. Polylogarithmic requested accuracy and a window of length $N/\operatorname{poly}(\log N)$ give polynomial bit complexity.

This evaluates a coarse average of the same magnitude seen by the decoder. The error proof removes the hidden high-frequency terms. On the remaining inputs, repeated coarse queries return explicit public diagonal sums and do not expose those phases. Any use of a smaller error after deconvolution must charge the resulting precision and window length.

The evaluator also supports a supplied Fourier center at the same polynomial-bit scale. Translated quadratic lattice blocks replace the quotient scan in general bands. A hit yields a factor, while independently sampled centers and small rational centers have small coverage. [Modulated averages and coverage](modulated-phase-window.md)

## Checks

The diagnostic compares 36 averages with direct summation over unit scales, including wraparound and a repeated prime. It also tests a low mixed line recovered without trial division. Ten further long-window cases compare the unfactored approximation with factor-supplied Fourier sums. [Evaluator and controls](coarse_phase_window.py)

```sh
uv run --no-project --with sympy --with python-flint==0.9.0 --with gmpy2 python experiments/polylog_factoring_2026_11_27/coarse_phase_window.py
```
