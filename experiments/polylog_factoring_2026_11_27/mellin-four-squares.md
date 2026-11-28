# What Mellin smoothing reveals about the positive scalar count

The positive prefix counts sum $8f_N(t)$, where $f_N(t)=\sigma(\gcd(t,N)^2)$. For odd $N$ and $\Re s>1$, its normalized Dirichlet series is the finite polynomial

$$D_N(s)=\frac{1}{\zeta(s)}\sum_{t\ge1}\frac{f_N(t)}{t^s}
=\sum_{d\mid N}\eta(d)d^{2-s},\qquad \eta(d)=\prod_{r\mid d}(1+1/r).$$

Indeed $f_N(t)=\sum_{d\mid\gcd(t,N)}d^2\eta(d)$. Absolute convergence allows summation over multiples of each $d$, which cancels $\zeta(s)$. Both lines studied here lie entirely in this domain. The unweighted Fourier transform of logarithms of divisors has established uses in divisor distribution. The earlier repo kernel tested each candidate integer for divisibility and then required summing a long window. Here the question is the numerical accuracy available for one aggregate on different vertical lines. [Positive-prefix identity](four-square-prefix.md) · [Mangerel, §5.3](https://doi.org/10.1007/s40687-023-00376-0) · [Earlier Poisson kernel](analytic-divisor-window.md)

| Quantity | Computation and scope |
| --- | --- |
| $\lvert D_N(3+i\pi/\log a)\rvert$, any odd $N$, integer $a\ge2$ | **Proved and implemented.** A deterministic algorithm using $N,a,\epsilon$ returns rational bounds whose ratio is at most $1+\epsilon$, in time polynomial in the input lengths and $1/\epsilon$. |
| $D_N(\sigma)$, real $\sigma>1$, $N=pq$ with distinct odd primes $p<q\le2p$ | **Proved.** A formula in $N$ and $\sigma$ has relative error $O(1/\log^2 N)$ uniformly in $\sigma$. |
| $\lvert D_N(2+i\pi/\log a)\rvert$ on the same semiprimes | **Oracle reduction proved and implemented.** A factor-of-two estimate at each of $O(\log N)$ adaptively chosen $a$ recovers $p$. The required frequencies have size $O(1/\log N)$. Computing these estimates at the small magnitudes in polynomial bit time remains open. |

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

**Checks and research decision.** Run `python experiments/polylog_factoring_2026_11_27/check_mellin_magnitude.py`. The test uses the factors to enclose the oracle answers independently, including adverse multiplicative errors on 60- and 90-bit inputs. It also checks the public right-line enclosures against independently bounded values, including inputs with large untreated cofactors. The first unmet operation is a polynomial-bit-time estimate of the magnitude near cancellation from the unfactored input. More real exponents or low-frequency measurements with only the above absolute guarantee repeat the normalization problem. [Verification script](check_mellin_magnitude.py)

## What continuation costs

The earlier Poisson-kernel experiment put a divisor window into the complex plane and counted the poles crossed by moving its contour. The Mellin decoder above located a balanced prime from an assumed magnitude at real part two. Here we test whether actual measurements at real part three can supply that missing information. If the measurements have inverse-polynomial accuracy in $\log N$, they can all equal 1 on balanced semiprimes. With sufficiently many precision bits proportional to $\log N$, measurements on one fixed frequency interval suffice for general factoring. [Earlier contour calculation](analytic-divisor-window.md#moving-the-sum-into-the-complex-plane) · [Earlier coefficient-extraction bound for a different observation model](../../notes/nativeana_results.md#45-prop-hw-hankel-rank-uncertainty-for-coefficient-extraction-proved)

**What the cheap readings reveal.** For $N=pq$ with $p<q\le2p$, put $R=\sqrt N$ and $u_r=(1+1/r)/r$. At every real frequency $t$,
$$
|D_N(3+it)-1|\le u_p+u_q+u_pu_q\le 4/R\qquad(R\ge8).
$$
Indeed $r\ge R/\sqrt2$ for both primes. Squaring the magnitude gives $\big||D_N(3+it)|^2-1\big|\le10/R$ for $R\ge8$. A source permitted that absolute error can return 1 to every query, including adaptive queries at arbitrary frequencies. Any algorithm guaranteed to factor using those answers must still succeed when supplied with that constant. This statement concerns information guaranteed by the measurement interface. The input integer remains available for other arithmetic.

**The semiprime case already works at zero frequency if enough precision is supplied.** On a squarefree input each local factor of $D_N(3)$ is $(r^2+r+1)/r^2$. Thus $N^2D_N(3)=\sigma(N^2)$, exactly the zero-scalar count from the earlier four-square calculation divided by eight. On $N=pq$ with distinct primes and $S=p+q$, this integer equals $N^2+NS+S^2-N+S+1$. Consecutive candidate sums change $D_N(3)$ by more than $1/N$. One real reading to absolute error $1/(4N)$ recovers $S$ by binary search, and the roots of $X^2-SX+N$ give the primes. This uses the zero-frequency option of the extended evaluator. On semiprimes the obstacle was already visible in that single real number. [Earlier four-square identity](four-square-prefix.md)

### Two real readings per split

Write $H(m)=m^2D_m(3)=\sum_{d\mid m}(m/d)^2\psi(d)$ with $\psi(d)=d\eta(d)$. This is an integer. In the earlier four-square construction with three annihilating scalars, the number of tuples is precisely $8mH(m)$, by summing the counts at each content. That construction already gives an efficient uniform sampler. Here the required operation is to compute its total count to the last integer. [Earlier sampler and content counts](four-square-incidence.md#which-points-carry-a-factor)

The function $H$ is multiplicative. Its value at $r^e$ is
$$
H(r^e)=r^{e-1}\left[r^{e+1}+(r+1)\sum_{j=0}^{e-1}r^j\right].
$$
When $e\equiv1\pmod 3$, the bracket is divisible by $\Phi_3(r)=r^2+r+1$. Modulo this integer, $r^3\equiv1$, and every block $1+r+r^2$ vanishes. Thus the sum is congruent to 1 and $r^{e+1}$ is congruent to $r^2$.

For a noncube $N$ some $r^e\parallel N$ has $3\nmid e$. Either $e$ or $2e$ is 1 modulo 3, so $\Phi_3(r)$ divides $H(N)H(N^2)$. In fact the product covers every prime with $3\nmid e$. For each $m=N,N^2$, a real reading within $1/(4m^2)$ of $D_m(3)$ recovers $H(m)$ by rounding. The product has polynomial bit length because $H(m)\le m^3$. A single degree-three norm-one group suffices in the existing randomized splitter whenever at least two distinct primes remain. Strip powers of two, recognize primes and perfect powers by the usual polynomial-time algorithms, and repeat on composite pieces. This yields general expected-polynomial-time factoring conditional on a supplier whose cost is polynomial in the number of precision bits. The splitting algorithm comes from earlier work. The arithmetic observation here is the divisibility supplied by the pair of inputs. [Earlier proof for a supplied exponent](translated-batches.md#factoring-with-an-exponent-for-a-norm-one-group) · [Bach, Miller and Shallit](https://doi.org/10.1137/0215083)

Passing to powers also preserves the accuracy obstruction. For a balanced $N=pq$, at every $t$ and every integer $j\ge1$ each local factor of $D_{N^j}(3+it)$ differs from 1 by at most $(1+1/r)/(r-1)\le2/r$. With $R=\sqrt N\ge8$, the entire product differs from 1 by at most $4\sqrt2/R+8/R^2\le7/R$. Thus adding powers or frequencies with the same coarse tolerance cannot force a reading distinguishable from 1.

The check uses 52 pairs of real readings, including errors toward opposite ends of the promised intervals. For $N=77{,}175=3^2 5^2 7^3$, none of the three cubic orders divides $H(N)$. The second reading supplies the orders at 3 and 5 and the existing coefficient gcd returns 225. For the exact cube $3375=15^3$, the product misses both orders, and taking the integer cube root first gives a valid split of 15. [Two-reading routine](mellin_continuation.py) · [Controls](check_mellin_continuation.py)

### What the complex readings add

**Continuation from magnitude alone.** Fix $t$ and form the entire function
$$F_t(x)=D_N(3+i(t+x))D_N(3-i(t+x)).$$
For real $x$ it is the measured squared magnitude. At $x=ih$ it equals $D_N(3-h+it)D_N(3+h-it)$. For fixed $h=1$ or $3$, use the $m$ real nodes $j/m$, $0\le j<m$, and evaluate their Lagrange interpolant at $ih$. Each term of $F_t$ has frequency $\log(e/d)$ for divisors $d,e$ of $N$. The absolute frequencies are at most $\log N$, the sum of the absolute coefficients is at most $N^2$ since $\eta(d)\le d$, and the growth through height $h$ is at most $N^h$. With $b$ the bit length, the interpolation error for exact samples is therefore at most
$$
N^{2+h}\frac{((h+1)b)^m}{m!}.
$$
This follows from the integral formula for complex divided differences on the convex hull of the nodes and $ih$. The sum of the absolute Lagrange weights is at most $(2\sqrt{1+h^2}\,m)^{m-1}/(m-1)!$, bounded by $10^m$ at height one and $24^m$ at height three. To reach absolute error $\tau$, take $m=O(b+\log(1/\tau))$ to make the remainder at most $\tau/2$, then bound each reading error by $\tau/(2\cdot 24^m)$. Both the number of readings and their precision bits are linear in this bound. The rational weights at these two heights have polynomial bit length. This is classical interpolation applied to this finite divisor sum. [de Boor, formula (52)](https://arxiv.org/abs/math/0502036) · [Demanet and Townsend on extrapolation from perturbed samples](https://arxiv.org/abs/1605.09601)

**Uniform relative accuracy at real part two.** For any odd $N$ and real $t$, $|D_N(2+it)|\ge1/N$. To see this, write the factor at $r^k\parallel N$ as $L=1+c\sum_{j=1}^k z^j$ with $c=1+1/r$, $z=e^{i\theta}$ and $|\theta|\le\pi$. For $k=1$ the reverse triangle inequality gives $|L|\ge1/r$, attained when $z=-1$. For $k\ge2$ and $|\theta|\le\pi/(k+1)$ the sum of the cosines is nonnegative by the geometric-sum formula, giving $\Re L\ge1$. Outside that arc, $|1-z|\ge2\sin(\pi/(2(k+1)))\ge3/(k+1)$. The identity $(1-z)L=1+z/r-cz^{k+1}$ gives
$$
|L|\ge\frac{c-|1+z/r|}{|1-z|}
=\frac{|1-z|}{r(c+|1+z/r|)}
\ge\frac{3}{2(r+1)(k+1)}\ge r^{-k}.
$$
For the last inequality, $3r^k/[2(r+1)(k+1)]$ increases in both arguments and exceeds one at $r=3,k=2$. Multiply over the primes. At real part four the local distances from 1 sum to at most $\sum_{j\ge3}1/[j(j-1)]=1/2$, so $|D_N(4-it)|\ge1/2$. At height one the continued product thus has magnitude at least $1/(2N)$. Continuing to absolute error $\epsilon/(50N)$ and dividing its magnitude enclosure by a relative enclosure on the fourth line gives a relative enclosure on the second line. The implementation accepts a rational offset from $\pi/\log a$, or a rational frequency itself. [Public readings](mellin_right_line.py) · [Certified continuation](mellin_continuation.py)

**Why precise readings on a fixed interval suffice for all odd inputs.** Set $t=0$ and $h=3$. Now $F_0(3i)=D_N(0)D_N(6)$ is positive and real. For each $r^k\parallel N$, the local factor of $D_N(0)$ is $1+r+\cdots+r^{2k}$, so $D_N(0)=\sigma(N^2)$. In $D_N(6)$ the denominator of the summand $\eta(d)/d^4$ divides $d^4\operatorname{rad}(d)$, which divides $N^5$. Hence
$$
E=N^5 F_0(3i)\quad\text{is a positive integer multiple of }\sigma(N^2).
$$
Continuing to absolute error at most $1/(4N^5)$ recovers $E$ by rounding. This uses $O(b)$ squared-magnitude readings on $[0,1)$, each with $O(b)$ bits after the binary point. The integer $E$ has $O(b)$ bits and is divisible by $(r^{2k+1}-1)/(r-1)$ for every $r^k\parallel N$. The existing norm-one group algorithm can try each odd degree through $2b+1$ to split a composite with two distinct prime divisors. Use primality and perfect-power tests and recurse on the pieces to get full factorization in expected polynomial time given the readings. Conversely, the prime powers allow each reading at the specified accuracy by finite geometric sums in polynomial bit time. The extraction of a factor from the exponent uses the earlier Bach–Miller–Shallit method. [Existing splitter and its arbitrary-input proof](translated-batches.md#factoring-with-an-exponent-for-a-norm-one-group) · [Bach, Miller and Shallit, Algorithm B](https://www.cs.cmu.edu/~glmiller/Publications/Papers/BMS86.pdf#page=6)

**Research decision.** The present evaluator pays for inverse tolerance by trial-dividing. At the precision used in the fixed-interval reduction its cutoff exceeds the smaller prime on balanced inputs, so the numerical experiment already finds that prime while generating its readings. A new evaluator with polynomial cost in the requested precision bits would change the result. Another way of decoding the same readings leaves that cost in place. [Reproducible checks](check_mellin_continuation.py)

The checks recovered the exact exponent on 23 inputs and confirmed the real-axis identity on 44 squarefree inputs. They also compared six critical-line enclosures against independently bounded divisor sums. For $N=1{,}242{,}079=1009\cdot1231$ at $a=1009$, the continuation used 164 readings and 609 working precision bits for a one-eighth relative enclosure. The native source had already reached 1009 by trial division. Run `python experiments/polylog_factoring_2026_11_27/check_mellin_continuation.py` to reproduce these controls.
