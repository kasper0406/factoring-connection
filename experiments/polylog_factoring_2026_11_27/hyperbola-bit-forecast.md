# How much do the size of the product and a small modulus predict?

Draw $p,q$ independently and uniformly from a common dyadic interval $[R,2R)$ and condition on $p\ne q$. Here $R$ is a power of two and $N=pq$. Fix a positive odd $h$ and a dyadic period $D$ dividing $R$. Write $\lambda=hR/D$, $H(u)=1_{\{u\}\ge1/2}$, and $Y=H(hp/D)H(hq/D)$. At an eligible level this is the Boolean binomial reading. The input-only score based on $\sqrt N$ has a genuine correlation with $Y$. Here we compute a stronger benchmark and determine how much information a broad class of scores can extract in this regime. [Square-root phase](upper-bit-square-root-phase.md) · [Binomial bit identity](two-fixed-and-streams.md) · [Earlier tests that certify feasible bits](cofactor-bit-completions.md)

**The forecast.** Take two continuous variables on $[R,2R]$ of density proportional to $1/\log x$. Condition on their product being exactly $N$. Let $\mu(N)$ be the resulting probability of both upper-half bits. This function is computable from the public input in $O(1+\lambda)$ intervals, with polynomial bit cost per interval at polynomial precision. In the common dyadic band model, $R=2^{\lfloor(\log_2N)/2\rfloor}$ is itself determined by $N$. The broader promise $p<q<2p$ also allows the factors to straddle a dyadic boundary.

For the true prime sample, for each fixed $A,C,B>0$, suppose $\lambda,m\le(\log R)^A$ and $m$ is a predetermined modulus. Allow an arbitrary function $G_c(N/R^2)$ for each residue $c=N\bmod m$, bounded by one and of total variation at most $(\log R)^C$ on $[1,4]$. Then uniformly over these scores,
$$
\left|\mathbb E[(Y-\mu(N))G_{N\bmod m}(N/R^2)]\right|=o((\log R)^{-B}).
$$
Thus the gain from the square-root score is explained by the real product constraint in this model. For every eligible forecast $G\in[0,1]$, its mean squared loss minus that of $\mu$ is $\mathbb E[(G-\mu)^2]+o((\log R)^{-B})$. This covers nonlinear forecasts and arbitrary combinations of tests modulo a common small modulus within the variation bound. Full integer algorithms with fine-scale dependence on $N$ are outside the result.

**Amount of information.** When $\lambda\to\infty$ within the stated range, the gain in mean squared loss of $\mu$ over the best constant forecast is
$$
\Theta\!\left(\frac{\log\lambda}{\lambda}\right).
$$
The positive and negative bounds concern the same prediction target and input distribution. The gains are summable as $D$ decreases dyadically. A better scoring rule confined to these inputs for the score cannot turn the successive upper levels into a long run of reliable binomial readings. This is a bound on marginal squared-error gains. It makes no assertion about an arbitrary joint-transcript decoder.

## Computing the conditional forecast

In the continuous model the first coordinate $x$, conditional on the product $N$, has density proportional to $1/[x\log x\log(N/x)]$ on $[\max(R,N/(2R)),\min(2R,N/R)]$. An antiderivative, up to the common factor $1/\log N$, is
$$
F_N(x)=\log\!\left(\frac{\log x}{\log(N/x)}\right).
$$
Intersect the intervals in which $hx/D$ lies in the upper half with the reciprocal images of the same intervals for $hN/(Dx)$. Both lists are sorted and have $O(1+\lambda)$ members. Summing $F_N(b)-F_N(a)$ on their intersections and dividing by the total gives $\mu(N)$. The reference program compares rational endpoints exactly and uses adjustable decimal precision on the integrals. [Implementation](hyperbola_bit_forecast.py)

At the upper corner, $N>4R^2-RD/h$ certifies $Y=1$ because both primes exceed $2R-D/(2h)$. Predicting YES on that corner and NO elsewhere therefore improves classification accuracy over always answering NO by $(1+o(1))/(8\lambda^2)$. In coordinates $p/R=2-u$, $q/R=2-v$, the corner has area $(1+o(1))/(8\lambda^2)$ and the normalized prime densities tend uniformly to one. The corresponding lower corner has $\mu=0$. These are hard decisions on small regions. An explicit integer rule using additional short intervals attains the optimal asymptotic gain in hard-decision accuracy for the stated model. [Hard decisions and their exact constant](hyperbola-hard-decisions.md)

## Transfer to prime inputs and the role of the small modulus

Use the log coordinate $t=\log(x/R)$. On $[0,\log2]$ let $f(t)$ be the normalized density proportional to $e^t/\log(Re^t)$ and let $g(t)=f(t)H(\lambda e^t)$. Extend both by zero. For $s=\log(N/R^2)$ the forecast is $(g*g)(s)/(f*f)(s)$. The numerator and denominator are Lipschitz with constants $O(1+\lambda)$ and $O(1)$, respectively. Indeed translation in $L^1$ of a bounded-variation function costs at most its variation times the shift. The denominator is bounded below by a constant times $\ell(s)=\min(s,2\log2-s)$, the length of the overlap of the log intervals. In one-sided neighborhoods of size $c/\lambda$ of the two endpoints the ratio is constant, because $\lambda$ is integral and both factors lie in their respective first or last half cells. Extend by those limits at the endpoints. Elsewhere division using $0\le\mu\le1$ gives $\operatorname{Lip}(\mu)=O((1+\lambda)^2)$, also in the product coordinate on $[1,4]$.

For each pair of unit residue classes modulo $m$, partial summation replaces the prime sums in each variable by the integrals on $[R,2R]$ against $dx/(\varphi(m)\log x)$. The Siegel–Walfisz theorem makes the normalized discrepancies smaller than any fixed inverse power of $\log R$, uniformly for the stated $m$. For fixed values of the other coordinate the variations of the integrands are bounded by a fixed power of $\log R$, by the half-cell count and the preceding Lipschitz bound. Sum over the at most $m^2$ residue pairs. The smooth model uses the same independent density for every pair of classes, and its conditional residual has mean zero at each product. Removing equal primes costs $O(\log R/R)$. Applying the same estimate with $G-\mu$ proves the squared-loss identity. [Siegel–Walfisz, in Koukoulopoulos, *The Distribution of Prime Numbers*, chapter 12](https://dms.umontreal.ca/~koukoulo/documents/publications/primes.pdf)

## Size of the best possible gain in this model

Put $z=N/R^2$, use dimensionless $x\in[a,b]=[\max(1,z/2),\min(2,z)]$, and put $w(x,z)=(\log R)^2/[x\log(Rx)\log(Rz/x)]$. This density is bounded above and below and has uniformly bounded first derivatives on the window. Its integral is comparable to $\ell(\log z)$. The unnormalized residual is
$$
J(z)=\int_a^b w(x,z)\left[H(\lambda x)H(\lambda z/x)-1/4\right]dx.
$$
Uniformly in the product, $|J(z)|\ll\lambda^{-1/2}$. To see this, expand $H$ in its constant coefficient $1/2$ and coefficients $i/(\pi m)$ for odd $m\ne0$. This may first be done with Fejér polynomials, then passed to the limit by bounded convergence. Nonconstant terms have phase $\lambda(mx+nz/x)$. If one frequency is zero, their signs are opposite, or $|m/n|$ is outside $[1/4,4]$, the derivative is bounded below in magnitude by $c\lambda\max(|m|,|n|)$ since $z/x^2\in[1/2,2]$. Integration by parts and summation of the absolutely weighted bounds cost $O(1/\lambda)$. For the remaining terms, the second derivative test gives $O((\lambda |n|)^{-1/2})$ uniformly in the endpoints. Summing $1/|m|$ over comparable frequencies costs $O(1)$ per $n$, and $\sum_{n\ge1}n^{-3/2}$ converges.

Consequently $|\mu-1/4|\ll\min(1,1/(\sqrt\lambda\,\ell))$. The continuous product density in log coordinates is $O(\ell)$. Integrating the square with a split at $\ell=1/\sqrt\lambda$ gives $O((1+\log\lambda)/\lambda)$. The same order bounds the optimal gain over the best constant on the prime sample by the transfer above.

For a matching lower bound let $\delta=\lambda^{-1/16}$, $L(z)=1/(z-1+\delta)$, and $S(z)=-\cos(4\pi\lambda\sqrt z+\pi/4)$. Test the residual $Y-1/4$ against $W=S L$. The two-dimensional stationary-phase calculation for the square-root score applies with the additional smooth amplitude $L(xy)$. Its first two derivatives cost $O(\delta^{-3})$. The stationary harmonic $(1,1)$ contributes a positive constant times
$$
\frac1{\sqrt\lambda}\int_1^{\sqrt2}
     \frac{u^2 w_R(u^2)^2}{u^4-1+\delta}\,du
\asymp\frac{\log\lambda}{\sqrt\lambda},
$$
where $w_R$ is the normalized single-coordinate prime density on $[1,2]$. With a Fejér cutoff $K=\lceil\lambda^2\rceil$ the smoothing, nonstationary, and stationary remainders are $O(\delta^{-1}\log K/K+\delta^{-2}\log K/\lambda+\delta^{-3}/\lambda)$. These are $o(\lambda^{-1/2})$. Also $\mathbb E W^2\ll\log\lambda$ by the product density near the corner. Conditional Cauchy–Schwarz therefore gives $\mathbb E(\mu-1/4)^2\gg(\log\lambda)/\lambda$. The actual mean of $Y$ differs from $1/4$ by $O(1/\lambda)$, so centering at that mean preserves the bound. All weights and their variations cost only powers of $\log R$ in the prime comparison. [Stationary-phase calculation and conventions](upper-bit-square-root-phase.md)

## Controls

An independent implementation split the integration at every change in either bit and matched the decimal implementation on 654 product/period settings. A separate seeded sample used 20,000 pairs of distinct primes in $[2^{19},2^{20})$. The table reports squared-loss gains over the best constant fitted to the sample. Both forecasts use the public product. Factor values supply the labels used for scoring.

| $h$ | $R/D$ | Square-root score | Hyperbola forecast |
|---:|---:|---:|---:|
| 1 | 4 | 0.009616 | 0.024991 |
| 1 | 16 | 0.002930 | 0.006453 |
| 1 | 64 | 0.000535 | 0.001244 |
| 3 | 4 | 0.003361 | 0.007500 |
| 3 | 16 | 0.000959 | 0.001873 |
| 3 | 64 | 0.000186 | 0.000559 |

The measured correlations of the centered residual $Y-\mu$ with the simple cosine were within 1.7 Monte Carlo standard errors of zero in each setting. The asymptotic scope follows from the proof, independently of this calibration. [Reproduce the checks](check_hyperbola_bit_forecast.py)
