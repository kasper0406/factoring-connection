# Predicting a dense bit of a shifted count

For $N=PQ$ with two complete odd prime-power components, let $S=P+Q$ and $P_N(A)=(P+A)(Q+A)$. Take a prime $\ell\ge N^{1/4}$ with $\gcd(\ell,N)=1$. A predictor with inverse-polynomial advantage for the sign of $\cos(2\pi P_N(A)/\ell)$ at a uniform $A\bmod\ell$ is sufficient to split $N$ in randomized polynomial bit time, including the progression search. A predictor that uses only $N\bmod\ell$ has exponentially small advantage on the balanced prime distribution below. The same obstruction covers short summaries of $N$ computed before seeing $A$.

## Recovery from weak prediction

Write $m=N\bmod\ell$, $e_\ell(t)=\exp(2\pi it/\ell)$, and $w(t)=\operatorname{sgn}\cos(2\pi t/\ell)$. The desired bit is $f_S(A)=w(A^2+SA+m)$. Supply the public value $w(m)$ at $A=0$ and use the ordinary shifted count at all other residues. Let $F(A)\in[-1,1]$ be the expected predicted sign. If the success probability is at least $(1+\rho)/2$, then $\langle F,f_S\rangle\ge\rho$. The errors can depend on the shift. Independent predictor randomness on each call suffices.

Let $c_h$ be the normalized Fourier coefficients of $w$. A geometric sum on the centered quarter interval gives $|c_0|=1/\ell$ and $|c_h|\le1/|h|_\ell$ otherwise, using centered representatives. Truncate at $H=\lceil16/\rho^2\rceil$. If $\ell\le2H+1$, enumeration of the sum residues and verification by the lattice cost polynomial bit time in $1/\rho$ and $\log N$. Otherwise the squared tail norm is at most $2/H$ on the field. The map $A\mapsto A^2+SA+m$ has at most two preimages per value, so the tail changes the correlation by at most $2/\sqrt H\le\rho/2$. The constant term contributes at most $\rho/4$.

For each $h$ set $F_h(A)=F(A)e_\ell(-h(A^2+m))$. The retained correlation is $\sum_{|h|\le H}\overline{c_h}\widehat F_h(hS)$. Some $1\le h\le H$ has
$$|\widehat F_h(hS)|\ge\frac{\rho}{8(1+\log H)},$$
by the $\ell^1$ bound and conjugation symmetry. For this $h$, at least one of the real and imaginary parts of $F_h$ has a Fourier coefficient of at least half that magnitude. Both expected functions can be sampled by multiplying the predicted sign by the appropriate public trigonometric value.

Use the shared-sample Fourier interval search at threshold $\tau=\rho/(16(1+\log H))$. For each of the $2H$ real functions, a sampled pool of $O(\tau^{-4}\log^2\ell\log(H\ell/\epsilon))$ pairs suffices for failure probability at most $\epsilon/(2H)$. The energy bound limits each search to $O(\tau^{-2}\log\ell)$ intervals. Invert $h$ on the returned frequencies and solve $X^2-SX+m$ in the field for each candidate $S$. The roots for the correct sum are the component classes. Use modular square roots with a random nonresidue and the native lattice progression search at step $\ell$. Verify each recovered component by division. This charges all $2H$ searches and all candidate verification. The resulting cost is polynomial in $\log(N\ell)$, $1/\rho$, $\log(1/\epsilon)$, and the prediction cost. [Fourier interval search](quarter-fourier.md) · [Native progression search](large_mod_residue_evaluator.py)

## Information in the public congruence

In the field model choose the two component residues uniformly and independently among units. Condition on their product $m$. The first is still a uniform unit $u$, and the second is $m/u$. For each nonzero shift $A$ and harmonic $h\ne0$, the conditional average of $e_\ell(h(A^2+A(u+m/u)+m))$ is a unit phase times a normalized Kloosterman sum. Its magnitude is at most $2\sqrt\ell/(\ell-1)$ by Weil's bound. The Fourier $\ell^1$ norm of $w$ is $O(\log\ell)$, so
$$|\mathbb E[f_S(A)\mid m]|=O(\ell^{-1/2}\log\ell)\qquad(A\ne0).$$
This bound is uniform in the public $m$ and nonzero $A$ in this model. [Weil's Kloosterman bound](https://doi.org/10.1073/pnas.34.5.204)

For uniform independent $A$, suppose the predictor also receives an arbitrary $b$-bit summary of the full input, chosen before seeing $A$. Fix $m,\ell$ and any one resulting bounded response function $G(A)$, and put $R(S)=\langle G,f_S\rangle$. For each nonzero $h$, Parseval and the bijection $S\mapsto hS$ give
$$\sum_{S\bmod\ell}\left|\widehat{G(A)e_\ell(-h(A^2+m))}(hS)\right|^2\le1.$$
Minkowski across harmonics, including the $O(\ell^{-1/2})$ constant contribution to the vector norm, yields $\sum_S|R(S)|^2=O(\log^2\ell)$. Each $S$ has probability at most $2/(\ell-1)$ given $m$, since $u+m/u=S$ is a quadratic equation. Bound the maximum over all $2^b$ response functions by the square root of the sum of their squared correlations. It follows that even the best summary for each input gives mean absolute correlation $O(2^{b/2}\ell^{-1/2}\log\ell)$ under the field model. The response functions and the summary can be arbitrarily expensive to compute. The summary must be independent of the sampled shift.

## Lift to balanced semiprimes

Let $X$ be dyadic, $I=[X,5X/4]$, and choose distinct primes $p,q$ uniformly in $I$. This fixes the bit length of $N=pq$. Independently sample a uniform prime $\ell\in[2\sqrt X,4\sqrt X]$. These auxiliary primes are large enough for the lattice and smaller than all factors for large $X$. Write $\Pi$ for the number of primes in $I$ and $\nu_\ell$ for their residue distribution. Character orthogonality and Cauchy–Schwarz give
$$d_{\rm TV}(\nu_\ell,U_\ell)^2\le\frac{1}{4\Pi^2}\sum_{\chi\ne1\bmod\ell}\left|\sum_{p\in I}\chi(p)\right|^2,$$
where $U_\ell$ is uniform on the units. Every nonprincipal character at a prime is primitive. The multiplicative large sieve bounds the sum of the right-hand numerators across the chosen $\ell$ by $O(X\Pi)$. By the prime number theorem, $\Pi\asymp X/\log X$ and the number of such $\ell$ is $\asymp\sqrt X/\log X$. The mean total variation is thus $O(X^{-1/4}\log X)$. The product distributions double the bound, and conditioning on distinct primes costs at most another $1/\Pi$. [Multiplicative large sieve](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

The field-model correlation bound transfers with that total-variation error. For the summary result, transfer the maximum over all response functions, a bounded function of the two residues and $\ell$. This allows the actual summary to depend on all digits of $N$. Over this balanced distribution the mean prediction advantage is
$$O\!\left((2^{b/2}+1)X^{-1/4}\log X\right).$$
An inverse-polynomial mean advantage therefore needs $b\ge\tfrac14\log_2 N-O(\log\log N)$ in this model. A fixed polynomial number of possible summary values supplies at most $O(\log\log N)$ bits. The bound leaves open algorithms that use the full input during each query.

An additive large-sieve argument on semiprime products also bounds long CRT transcripts from independently sampled auxiliary primes. The combined modulus can reach $N^{1-\eta}$ for fixed $\eta>0$, while conveying negligible additional information about fixed shifts. [Several public congruences](congruence-summary.md)

## Diagnostic

The sampled decoder split all four fixture moduli using their planted exact signs and again using planted fifth-harmonic signs. The second predictor had measured agreement from 56.7% to 61.4%, and the diagnostic selected the fifth harmonic explicitly. The proof searches all harmonics up to the bound when the predictor's errors are unknown. The fixtures reach 60 bits and include a repeated prime. A separate exhaustive control computes the maximum field-model bias at four small auxiliary primes. A native predictor that sees the full input remains open. [Sampled decoder and controls](dense_count_prediction.py)

```sh
uv run --no-project --with sympy --with python-flint==0.9.0 --with gmpy2 python experiments/polylog_factoring_2026_11_27/dense_count_prediction.py
```
