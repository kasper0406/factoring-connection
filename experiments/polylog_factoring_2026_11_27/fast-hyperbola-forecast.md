# A fast evaluator for the continuous upper-bit forecast

The geometric amplitude $J(z)$ in the [forecast proof](hyperbola-bit-forecast.md) has a short approximation at its natural scale even when the bit period is far below the factor size. For every $R\ge4$, product ratio $z=N/R^2\in[1,4]$, $\lambda=hR/D\ge1$ and odd-frequency cutoff $K\ge1$, the finite sum below approximates $\sqrt\lambda J(z)$ with **uniform absolute error**
$$
O(K^{-1/2}+\lambda^{-1/2}).
$$
It involves $O(K^2)$ terms independently of $\lambda$. When $\lambda$ is small one can integrate the $O(1+\lambda)$ half-cell intersections directly. Combining the two methods gives an approximation to $\sqrt\lambda J(z)$ of error $O(\epsilon)$ in time polynomial in $1/\epsilon$ and the input bit length. Equivalently it resolves the *unnormalized* centered geometric probability to $O(\epsilon/\sqrt\lambda)$. The conditional forecast divides this by the known mass of the product window. That denominator shrinks near the two edges. The object computed here is an integral over real pairs with the continuous prime-density weight. A single actual binomial bit depends on the prime divisors of the integer $N$. [Implementation](fast_hyperbola_forecast.py) · [Direct real integral](hyperbola_bit_forecast.py)

## Why only boundedly many oscillations need evaluating

Write $a=\max(1,z/2)$, $b=\min(2,z)$ and $w(x,z)=(\log R)^2/[x\log(Rx)\log(Rz/x)]$. Expand the bit functions in $J$ into Fourier modes. As in the proof of the forecast bound, the combined contribution of the zero-frequency axes, the opposite-sign pairs and the noncomparable pairs is $O(1/\lambda)$. For comparable positive odd $m,n$ with $1/4\le r=m/n\le4$, set
$$
u(x)=\sqrt{rx}-\sqrt{z/x},\qquad x_0=\sqrt{z/r},\qquad
G(u)=w(x,z)\frac{2x}{\sqrt{rx}+\sqrt{z/x}} .
$$
Then $mx+nz/x=2\sqrt{mnz}+nu(x)^2$. On the interval between $x_0$ and the endpoints, $x$ stays between $1/2$ and $4$ and the reciprocal $z/x$ stays in the same bounded range. The derivative of $u$ is bounded away from zero. For $R\ge4$ the logarithms in $w$ are bounded away from zero relative to $\log R$. Consequently $G$ and its first two derivatives are bounded absolutely, uniformly in $z,r,R$. These facts also apply when the stationary point $x_0$ lies just outside the window.

Replace $G(u)$ by $G(0)$ under the quadratic oscillation. The difference divided by $u$ extends differentiably through zero and has bounded variation. One integration by parts using $(e(\lambda n u^2))'=4\pi i\lambda n u e(\lambda n u^2)$ bounds the remainder by $O(1/(\lambda n))$. Here $e(t)=\exp(2\pi i t)$. Put $F(s)=\int_0^s e(v^2)dv$. After multiplying by $\sqrt\lambda$ and combining conjugate terms, the summand is
$$
-\frac{2G(0)}{\pi^2 mn\sqrt n}\operatorname{Re}\!\left[
 e(2\lambda\sqrt{mnz})
 \left(F(\sqrt{\lambda n}\,u(b))-F(\sqrt{\lambda n}\,u(a))\right)\right].
$$
Keeping the endpoint difference handles stationary points crossing either boundary. The weighted integration errors sum to $O(\lambda^{-1/2})$ in the scaled expression because $\sum n^{-2}\sum_{m\asymp n}1/m$ converges. The contribution of comparable frequencies with $\max(m,n)>K$ is $O(K^{-1/2})$ by the boundedness of $F$ and $\sum_{n>K/4}n^{-3/2}$. Fejér expansion followed by bounded convergence justifies the collective expression for discontinuous bit indicators.

To obtain error $O(\epsilon)$, use direct half-cell integration for $\lambda\le\epsilon^{-2}$ and otherwise use $K=\lceil\epsilon^{-2}\rceil$. Trigonometric phases only need polynomially many bits in the input length and $\log(1/\epsilon)$. The Fresnel integral has limits $\pm(1+i)/4$ and an $O(1/|s|)$ tail by integration by parts. It therefore suffices to use those limits for large arguments. On a bounded interval of length polynomial in $1/\epsilon$ ordinary quadrature at polynomially many points achieves the required accuracy by a first-derivative bound. The implementation uses multiprecision Fresnel functions and a separate tail guard. The approximation theorem bounds analytic truncation error, while the decimal precision setting controls numerical evaluation. [Fresnel functions and asymptotic expansions](https://dlmf.nist.gov/7.12)

## How far this informs the prime problem

For the same dyadic prime model, the ordinary prime number theorem with classical error $O(R\exp(-c\sqrt{\log R}))$ allows the preceding distributional forecasts to be transferred even if $\log\lambda=o(\sqrt{\log R})$. Indeed partial summation only multiplies the normalized error by a fixed polynomial in $\log R,\lambda$ for these scores and their variations. The same observation applies with a predetermined modulus of polylogarithmic size using the Siegel–Walfisz estimate. For the corresponding comparison of arbitrary scores one must bound their variation as well, for example by a fixed power of $\lambda+\log R$. In this range the error is still $o(\lambda^{-B}(\log R)^{-B})$ for every fixed $B$. Only the polylogarithmic $\lambda$ regime guarantees an inverse-polynomial squared-loss signal in the input length. The fast integral algorithm itself also works for much larger $\lambda$. A comparable claim about prime-pair prediction at those deeper levels would require a new prime distribution argument.

The direct and accelerated routines provide an input-only control for numerical tests. The comparison checks interior and boundary stationary points at several finite scales and evaluates a case with $\lambda=2^{80}$ at a fixed frequency cutoff. The latter is a complexity and precision control for the formula. [Reproduce the numerical checks](check_fast_hyperbola_forecast.py)
