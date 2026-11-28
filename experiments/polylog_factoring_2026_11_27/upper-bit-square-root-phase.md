# An input-only square-root phase predicts an upper binomial bit weakly

Take distinct primes $p,q$ uniformly in $[R,2R]$, write $N=pq$, and fix $h=1$ or $3$. At dyadic period $D$, let $Y$ indicate that both $hp$ and $hq$ fall in the upper half modulo $D$. In the range below this is exactly the Boolean reading $N\mid\binom{2t}{t}$ at $t=\lfloor hN/D\rfloor$. Put $\lambda=hR/D$. The following score uses only $N,h,D$.

$$
 S(N)=-\cos\!\left(\frac{4\pi h\sqrt N}{D}+\frac{\pi}{4}\right).
$$

As $R,\lambda\to\infty$ with $\log\lambda=o(\sqrt{\log R})$,

$$
 \mathbb E[(Y-\tfrac14)S(N)]
 =\frac{4(2\sqrt2-1)}{3\pi^2\sqrt2}\frac{1+o(1)}{\sqrt\lambda}.
 \tag{1}
$$

The constant is $0.174663\ldots$. The public score detects one harmonic of the hidden sum through the identity $p+q-2\sqrt{pq}=(\sqrt p-\sqrt q)^2$. The coherent contribution comes from near the diagonal in the sample of prime pairs. The correlation has inverse-polylogarithmic size when $\lambda\le(\log R)^A$ for a fixed $A$, covering upper $O(\log\log R)$ levels. The asymptotic extends farther with a smaller signal. The existing full-range decoder asks for many more bits and a much lower error rate. [Binomial reading and decoding hypotheses](two-fixed-and-streams.md#exact-search-and-the-binomial-observations)

**Probability forecast.** With $\widehat y=\frac14+\frac14\sqrt{D/h}\,N^{-1/4}S(N)$, the reduction in mean squared error relative to the best constant forecast $\mathbb E Y$ is

$$
 \frac{1}{\lambda}\left(\frac{1}{\sqrt2\pi^2}
                  -\frac{(\sqrt2-1)^2}{8}+o(1)\right).
 \tag{2}
$$

The leading coefficient is positive. For large $\lambda$ this forecast lies between zero and one. Computing a bounded-precision score costs polynomially many bit operations in $\log N$ for the accuracy required here. For example, to retain $w$ phase bits, compute $\lfloor\sqrt{(2h\,2^w)^2N}\rfloor$, divide by $D$ with integer flooring, and retain the low $w$ bits. The cosine on a single period can be computed to the same accuracy.

## Derivation

Write $e(z)=\exp(2\pi i z)$ and rescale $x=p/R$, $y=q/R$. The prime number theorem with its classical error replaces the two prime averages by integrals on $[1,2]^2$ with the normalized product density $w_R(x)w_R(y)$, where $w_R(x)$ is proportional to $1/\log(Rx)$. In each variable the integrands have total variation $O(1+\lambda)$, counting all half-interval jumps and the derivative of the exponential. Partial summation in each coordinate therefore costs $O((1+\lambda)\log R\,\exp(-c\sqrt{\log R}))$. Excluding the diagonal changes bounded expectations by $O(\log R/R)$. These errors are smaller than every inverse power of $\lambda$ in the stated range. Also $w_R=1+O(1/\log R)$ and its first two derivatives are uniformly bounded. [Prime number theorem with error term](https://www.math.leidenuniv.nl/~evertsejh/ant18-1.pdf)

Let $H$ be the periodic indicator of $[1/2,1)$. We evaluate the complex average of $(H(\lambda x)H(\lambda y)-1/4)e(-2\lambda\sqrt{xy})$. The Fourier coefficients of $H$ are $a_0=1/2$, $a_m=i/(\pi m)$ for odd $m$, and zero for the other frequencies. Use the Fejér polynomial of degree $K=\lceil\lambda^2\rceil$ for each indicator. Its one-period $L^1$ error is $O(\log K/K)$ because the $L^1$ norm of a translation difference of $H$ is $O(|u|)$ and the Fejér kernel is $O(\min(K,1/(K u^2)))$ on $|u|\le1/2$. The same error bound holds after composition with $\lambda x$ on $[1,2]$. Both approximants are bounded by one.

The phase of the $(m,n)$ term is $\lambda(mx+ny-2\sqrt{xy})$. For $m\ne1$ its $x$ derivative has magnitude at least $c\lambda(1+|m|)$ because $\sqrt{y/x}\in[1/\sqrt2,\sqrt2]$. Integration by parts bounds this term, before Fourier weights, by $O(1/(\lambda(1+|m|)))$. If $m=1$ and $n\ne1$, integrate in $y$. Summing the weighted bounds gives $O(\log K/\lambda)$. The constant subtraction also obeys the nonstationary bound. These estimates allow any fixed smooth amplitude with uniformly bounded first two derivatives.

The only stationary harmonic is $(1,1)$, with coefficient $a_1^2=-1/\pi^2$ up to $O(1/K)$ from smoothing. For smooth amplitude $W(x,y)$, put $u=\sqrt x,v=\sqrt y$ and integrate along $z=u-v$. The resulting integral is $\int g(z)e(\lambda z^2)\,dz$, supported on $|z|\le\sqrt2-1$, with

$$
 g(0)=4\int_1^{\sqrt2}u^2W(u^2,u^2)\,du .
$$

The amplitude $g$ is piecewise twice continuously differentiable with uniform bounds on each side of zero. Subtract $g(0)$ on its support. On each side the residual divided by $z$ has bounded variation, so one integration by parts gives $O(1/\lambda)$. The constant piece contributes $g(0)e^{i\pi/4}/\sqrt{2\lambda}+O(1/\lambda)$ by the Fresnel integral and integration by parts in the tails. Taking $W=w_R(x)w_R(y)$ gives $g(0)=4(2\sqrt2-1)/3+O(1/\log R)$. Multiplying the complex average by $-e^{-i\pi/4}$ and taking its real part proves (1).

For (2), let $a(x,y)=(xy)^{-1/4}$. The same calculation with amplitude $a\,w_R(x)w_R(y)$ replaces $4(2\sqrt2-1)/3$ by $4\int_1^{\sqrt2}u\,du+o(1)=2+o(1)$. Also $\mathbb E[a^2S^2]=2(\sqrt2-1)^2+o(1)$, by the identity for $\cos^2$ and a nonstationary integral for its oscillating term. Finally $\mathbb E Y-1/4=O(1/\lambda)$ from averaging the half-intervals. Expanding the squared error with $\widehat y=1/4+aS/(4\sqrt\lambda)$ gives (2). Among linear adjustments using $S$ alone, the best possible squared-error gain on this sample is $O(1/\lambda)$ since $\operatorname{Var}(S)\to1/2$. These gains are geometrically summable as the dyadic periods decrease within the proved range.

## Scores centered at other factor ratios

For a fixed positive odd integer $v$, use $S_v=-\cos(4\pi h\sqrt{vN}/D+\pi/4)$. The same argument gives $\mathbb E[(Y-1/4)S_v]=(C_v+o(1))/\sqrt\lambda$, with

$$
 C_v=\frac{4}{3\pi^2\sqrt2\,v^2}
 \sum_{\substack{mn=v,\ m,n>0\\1/2<m/n<2}}
 \big((2\min(m,n))^{3/2}-(\max(m,n))^{3/2}\big).
$$

Indeed the stationary harmonics are precisely these ordered divisor pairs. Their phase is $(\sqrt{mx}-\sqrt{ny})^2$, and all have the same Fourier coefficient $-1/(\pi^2v)$. After changing variables to $\sqrt{mx},\sqrt{ny}$, their diagonal overlaps from $\sqrt{\max(m,n)}$ to $\sqrt{2\min(m,n)}$ and the Jacobian on that diagonal is $4u^2/v$. Terms with one sufficiently large frequency are bounded as before. Among the remaining finitely many nonstationary terms the gradient has positive minimum, so integrating by parts along the gradient gives $O_v(1/\lambda)$ including the square's boundary terms. For instance $v=15$ detects the ratios $q/p=5/3$ and $3/5$ and $C_{15}>0$. The fundamental signal corresponds to $C_1$. The constants satisfy $C_v\ll\tau(v)v^{-5/4}$ since each contributing $m,n$ is comparable to $\sqrt v$. Thus even their unsigned leading coefficients sum to a finite number over all odd $v$. The asymptotic is for each fixed $v$.

The checker took all 2,596,932 ordered distinct prime pairs in $[2^{14},2^{15}]$ and all 9,177,870 in $[2^{15},2^{16}]$. It used both multipliers and $R/D\in\{16,32,64\}$, without fitting any coefficient. Across the twelve cases $\sqrt\lambda\,\mathbb E[(Y-1/4)S]$ ranged from 0.1475 to 0.1769. The forecast improved on the best constant on each cohort. These finite-sample numbers use the factors only to score the correct Boolean answers. [Reproduce the twelve comparisons](check_upper_bit_square_root_phase.py)
