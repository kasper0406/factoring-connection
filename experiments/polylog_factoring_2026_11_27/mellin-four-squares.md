# What Mellin smoothing reveals about the positive scalar count

The positive prefix counts sum $8f_N(t)$, where $f_N(t)=\sigma(\gcd(t,N)^2)$. For odd $N$ and $\Re s>1$, its normalized Dirichlet series is the finite polynomial

$$D_N(s)=\frac{1}{\zeta(s)}\sum_{t\ge1}\frac{f_N(t)}{t^s}
=\sum_{d\mid N}\eta(d)d^{2-s},\qquad \eta(d)=\prod_{r\mid d}(1+1/r).$$

Indeed $f_N(t)=\sum_{d\mid\gcd(t,N)}d^2\eta(d)$. Absolute convergence allows summation over multiples of each $d$, which cancels $\zeta(s)$. Both lines studied here lie entirely in this domain. The unweighted Fourier transform of logarithms of divisors has established uses in divisor distribution. The earlier repo kernel tested each candidate integer for divisibility and then required summing a long window. Here the question is the numerical accuracy available for one aggregate on different vertical lines. [Positive-prefix identity](four-square-prefix.md) · [Mangerel, §5.3](https://doi.org/10.1007/s40687-023-00376-0) · [Earlier Poisson kernel](analytic-divisor-window.md)

| Quantity | Computation and scope |
| --- | --- |
| $\lvert D_N(3+i\pi/\log a)\rvert$, any odd $N$, integer $a\ge2$ | **Proved and implemented.** A deterministic algorithm using $N,a,\epsilon$ returns rational bounds whose ratio is at most $1+\epsilon$, in time polynomial in the input lengths and $1/\epsilon$. |
| $D_N(\sigma)$, real $\sigma>1$, $N=pq$ with distinct odd primes $p<q\le2p$ | **Proved.** A formula in $N$ and $\sigma$ has relative error $O(1/\log^2 N)$ uniformly in $\sigma$. |
| $\lvert D_N(2+i\pi/\log a)\rvert$ on the same semiprimes | **Oracle reduction proved and implemented.** A factor-of-two estimate at each of $O(\log N)$ adaptively chosen $a$ recovers $p$. The required frequencies have size $O(1/\log N)$. Computing these estimates at the small magnitudes remains open. |

## A public evaluator at real part three

Put $z_r=r^{-1-i\pi/\log a}$. For $r^e\parallel N$ the local factor is $1+(1+1/r)\sum_{j=1}^e z_r^j$. Its distance from 1 is at most $(1+1/r)/(r-1)$. In particular every local magnitude is at least $1/3$ for odd $r$.

Trial-divide through $T=\lceil16b/\epsilon\rceil$, with $b$ the bit length of $N$ and $0<\epsilon\le1$. Use interval arithmetic to evaluate the known factors and retain the unfactored cofactor. The sum $u$ of the distance bounds over its primes is at most $2b/T\le\epsilon/8$. The unknown product of magnitudes lies in $[1-u,1/(1-u)]$ by the product inequality $\prod(1-x_i)\ge1-\sum x_i$. Multiply these bounds by an enclosure of the known product with upper-to-lower ratio at most $1+\epsilon/4$. The resulting ratio is at most $(1+\epsilon/4)/(1-\epsilon/8)^2\le1+\epsilon$.

The implementation uses outward-rounded elementary functions and Horner evaluation with $|z_r|\le1/3$. There are at most $b$ local powers to process, the phases are $O(b)$, and each local magnitude stays above $1/3$. Polynomially many precision bits suffice for the enclosure, including the argument reduction and accumulated product error. The trial budget scales with inverse tolerance. [Public evaluator](mellin_right_line.py)

## What positive real smoothing computes cheaply

For the balanced semiprime promise, set $R=\sqrt N$, $v=(2-\sigma)\log R$ and $c=\log(q/p)/\log N$. Dropping $1+1/p,1+1/q$ from the two local factors changes the product by a ratio between 1 and $(1+\sqrt{2/N})^2$. The remaining product divided by the public $(1+e^v)^2$ equals

$$1+\frac{\cosh(cv)-1}{1+\cosh v}\le1+c^2.$$

The inequality follows termwise from the series for $\cosh$ and $0\le c\le1$. Thus the public expression is a lower estimate for all real $\sigma>1$ with a common upper-to-lower ratio of $(1+\sqrt{2/N})^2(1+(\log 2/\log N)^2)$. Varying the real exponent preserves this guarantee.

## Locating a prime from the magnitude at real part two

For integer $a$ between $\lceil\sqrt{N/2}\rceil$ and $\lfloor R\rfloor$, set $t=\pi/\log a$. The two Euler factors have magnitudes $m_r=|1+(1+1/r)e^{-it\log r}|$ for $r=p,q$. Each satisfies

$$m_r^2=1/r^2+4(1+1/r)\sin^2\!\left(\frac{\pi\log(r/a)}{2\log a}\right).$$

At $a=r$ its magnitude is $1/r$. Nearby, the sine term increases with the distance to $r$. For $N\ge2^{16}$ put $W_r=b+10|r-a|$. Combining both terms gives $W_r/2\le br m_r\le3W_r/2$. To see the constants, write $(brm_r)^2=b^2+u^2(r-a)^2$. For $r\ne a$, $u$ is the product of $\pi b/\log a$, $\sqrt{1+1/r}$, $x\log(x)/(x-1)$ at $x=r/a$, and $\sin y/y$ at $y=\pi\log(x)/(2\log a)$. Here $1/\sqrt2\le x\le2$, $|y|\le\pi/15$ and $2/\log2\le b/\log a\le34/(15\log2)$. These bounds give $7\le u\le15$. The asserted comparison follows by Cauchy–Schwarz and the triangle inequality.

Scale the measured magnitude to $V=Nb^2 m_p m_q$, so $W_p W_q/4\le V\le9W_p W_q/4$. Start at $a=\lceil\sqrt{N/2}\rceil\le p$ and always approach $p$ from below. The public scale $S=b+10(\lceil R\rceil-a)$ satisfies $S\le W_q\le(5/2)S$, since $\lceil R\rceil\le q\le N/a$ and $N/a-a\le(1+\sqrt2)(R-a)$. If the oracle supplies $V/K\le Y\le K V$, then

$$\frac{b+10(p-a)}{4K}\ \le\ Z:=Y/S\ \le\ 6K(b+10(p-a)).$$

If $Z\le12Kb$, there are only $O(K^2 b)$ integers left to check by gcd. Otherwise advance $a$ by $\lceil(Z-6Kb)/(60K)\rceil$. This step never passes $p$ and covers at least $(p-a)/(48K^2)$. Consequently $O(K^2\log N)$ oracle calls and gcd checks suffice, with small $N$ handled directly. Even $K$ polynomial in the bit length is enough for polynomial-time recovery when the oracle runs in polynomial time and gives polynomial-length rational answers. [Exact decoder](mellin_magnitude.py)

**Accuracy at the first cancellation.** Writing $d=\frac12\log(q/p)$, the public approximation $(1+R^{-it})^2$ to $D_N(2+it)$ has absolute error at most $t^2 d^2+2(1/p+1/q)+1/N$. This follows by expanding both Euler factors and using $2(1-\cos(td))\le t^2d^2$. At $t=\pi/\log R$ the approximation is zero. An $O(1/\log^2 N)$ absolute error bound holds there regardless of the factors. By contrast, at $a=p$ the two factors have magnitudes $1/p$ and at most $1/q+O(\log(q/p)/\log p)$. Resolving the latter product to a multiplicative tolerance can require exponentially smaller absolute error in the bit length. The right-line evaluator's omitted factors can undergo precisely this cancellation if used at real part two.

**Connection to the earlier interval decoder (proof).** Suppose a critical-line measurement has an additional absolute error allowance $E$ in the *scaled* value, so $V/2-E\le Y\le2V+E$. While $Y>2E$, it is within a factor of four of $V$ and the same walk works. If $Y\le2E$, then $V\le6E$. Since $a\le p\le q$, the inequality $V\ge W_pW_q/4$ puts $p$ in the public interval $[a,a+\sqrt{24E}/10]$. For fixed $0<\gamma<1/4$, an additive allowance of $N^{-1/2-2\gamma}/b^2$ in the *unscaled* magnitude makes this interval $O(N^{1/4-\gamma})$ wide. The existing Coppersmith affine-hint result then recovers $p$ in polynomial time. This corollary prices the accuracy needed by the combination. The numerical checks below exercise the multiplicative-error walk. [Affine-hint decoder and its hypotheses](../../notes/trackE_candidates.md)

**Checks and research decision.** Run `python experiments/polylog_factoring_2026_11_27/check_mellin_magnitude.py`. The test uses the factors to enclose the oracle answers independently, including adverse multiplicative errors on 60- and 90-bit inputs. It also checks the public right-line enclosures against independently bounded values, including inputs with large untreated cofactors. The first unmet operation is a factor-blind estimate of the magnitude near cancellation. More real exponents or low-frequency measurements with only the above absolute guarantee repeat the normalization problem. [Verification script](check_mellin_magnitude.py)
