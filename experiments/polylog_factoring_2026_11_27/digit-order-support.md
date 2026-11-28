# Can the binary digits of $N$ supply primes from a factor's order?

Choose distinct primes $p,q$ uniformly in $[R,2R)$, put $X=R^2$, and write $N=pq$. At every binary cut let $H=2^k$ and form the two public integers $A=\lfloor N/H\rfloor$ and $B=N\bmod H$. Use $1\le k<\operatorname{bitlength}(N)$, so both outputs are positive. Even if the prime divisors of a chosen output are available for free, they cover a part of $p^2-1$ of size $N^t$ on at most $X^{-t+o(1)}$ of the inputs, for $0<t\le1/4$. This holds for the best single output over all cuts and either hidden prime. Pooling prime support from multiple cuts gives weaker bounds below. [Earlier polynomial-support theorem and mixed-sign reader](sparse-order-support.md) · [Earlier root-scale exponent test](factoring-with-an-order-hint.md#can-the-exponent-come-from-the-size-of-the-input) · [Executable checks](check_digit_order_support.py)

For a set $S$ of primes, the covered part is $[p^2-1]_S=\prod_{\ell\in S}\ell^{v_\ell(p^2-1)}$. An output divisible by $\ell$ only once still receives credit for the full valuation. The theorem counts availability under this generous scoring rule. It assumes no algorithm for factoring the binary outputs.

## Uniform spacing at every cut

The two end pieces admit the spacing argument below. A word cut from the interior has a separate collision bound, which permits pooling the prime divisors of *all* proper words at all positions. [Interior-block proof and direct gcd count](interior-order-support.md)

Fix $p$, $H\ge2$ a power of two, and $m\ge2$ dividing $p^2-1$. For each of the two outputs, at most $1+4R/m$ **prime** partners $q\in[R,2R)$ give a positive multiple of $m$. This includes $H$ comparable to $\sqrt N$ and cuts far from the middle.

For the quotient, write $A=mu>0$. The fractions $u/q$ lie in an interval of width $1/(mR)$ ending at $p/(Hm)$. Each is reduced. Indeed, $u=zq$ with $z\ge1$ would give $B=(p-Hmz)q>0$ because $p$ is odd and $H$ even. It would also give $H<p/m<R$, whereas $q\le B<H$. Distinct prime denominators therefore give distinct fractions spaced at least $1/(4R^2)$ apart. This proves the count.

For the remainder, $B$ is odd, so only odd $m$ can occur. The slopes $A/q$ lie in an interval of width $1/R$ ending at $p/H$. The congruence $HA\equiv pq\bmod m$ implies $m\mid Aq'-A'q$ for every pair of these slopes, since $H$ is invertible modulo $m$. Nonintegral slopes are reduced and distinct. An integral slope $A=zq$ forces $m=q$. If $z=0$, use $m\mid pq$ and $\gcd(p,m)=1$. If $z\ge1$, then $0<(p-Hz)q<H\le p<2R$, forcing $p-Hz=1$ and $B=q$. Thus there is at most one integral slope, and all distinct slopes have spacing at least $m/(4R^2)$. No information on prime counts in arithmetic progressions is needed.

## Support from a single digit output

Let $T=X^t$ and fix the target $p$. It suffices to count coverage of $T$, since $N\ge X$. If a supporting prime $\ell>T$ occurs, fix $p$ first. There are $O(\log R)$ primes dividing $p^2-1$, all at most $p+1\le2R+1$. The spacing estimate gives $O(R/\ell)$ possible prime partners for each cut, including the constant term at this scale. There are $O(\log R)$ cuts. Dividing by the $\Theta(R/\log R)$ available prime partners costs only logarithmic factors.

Otherwise multiply the supported prime factors, with repetitions, until the first divisor $d$ reaches $T$. Then $T\le d\le T^2\le R$, $p^2=1\bmod d$, and $m=\operatorname{rad}(d)$ divides the chosen digit output. There are at most $2^{1+\omega(d)}$ roots of unity modulo $d$, including the power of two. They allow $O(2^{\omega(d)}R/d)$ first primes by counting integers. The spacing estimate then allows $O(R/m)$ partner primes for each side of each cut. The root-count factor is $X^{o(1)}$ uniformly for $d\le R$.

Normalize by $\Theta(R^2/\log^2R)$ ordered distinct prime pairs. For $0<\sigma<1$ the remaining sum is at most
$$
 \sum_{d\ge T}\frac{1}{d\,\operatorname{rad}(d)}
 \le T^{-\sigma}\prod_\ell
 \left(1+\frac{\ell^{-2+\sigma}}{1-\ell^{-1+\sigma}}\right).
$$
This product converges. Letting $\sigma$ approach one gives $X^{-t+o(1)}$ even after choosing the more successful factor, cut and side for each input. The previous polynomial-value proof uses roots of $F(pq)$ at each supporting prime. Here division and discarding binary digits change the residue constraints, and the spacing argument supplies the needed partner count. [Polynomial support comparison](sparse-order-support.md#can-we-obtain-the-bases-from-the-input)

At $t=1/4$ the exponent is attained already by the fixed cut $H=2$, whose quotient is $(N-1)/2$. The earlier Bombieri–Vinogradov argument constructs at least $R^{3/2-o(1)}$ ordered prime pairs having some $d\mid p-1,q-1$ with $d\gg\sqrt R(\log R)^c$ for a fixed $c>0$. Then $d/\gcd(d,2)$ divides the quotient and $p^2-1$, and exceeds $N^{1/4}$ for large $R$. This matches the upper bound for the best-single-output family. The lower bound uses one end cut and provides a logarithmic margin at the quarter threshold. [Pairs with a common divisor at that scale](sparse-order-support.md#can-we-obtain-the-bases-from-the-input)

## Pooling every prefix and remainder

Give one pool the prime divisors of all $O(\log R)$ positive outputs. Suppose some $s$ of these primes account for a covered part at least $N^t$, where $1\le s=O(\log\log X)$. The probability is at most
$$
 \exp\!\left(-\big(t-o(1)\big)\frac{\log X}{s}\right).
$$
For $s=1$ the single-output result applies. For $s\ge2$, handle a supported prime above $T$ as before. A first-crossing divisor $d\le R$ has $\omega(d)\le s$. Its largest prime factor is at least $\operatorname{rad}(d)^{1/s}$ and must divide at least one output. Test this one prime over all cuts with the spacing bound. After counting the possible $p$ as before, the normalized sum is at most a polylogarithmic factor times
$$
 2^s \sum_{d\ge T}\frac1{d\,\operatorname{rad}(d)^{1/s}}.
$$
For $L=\log X$ take $\sigma=(1-1/\log L)/s$ in Rankin's inequality. The Euler product after pulling out $T^{-\sigma}$ has log at most $O(\sum_\ell\ell^{-1-1/(s\log L)})=O(\log(2s\log L))$, since $\sigma\le1/2$. This proves the bound uniformly in the stated range of $s$. It allows the scorer to choose the best subset after seeing both factors.

Even with **unrestricted pooled support**, the probability is at most
$$
 \exp\!\left(-\big(\sqrt{t/2}-o(1)\big)
                     \sqrt{\log X\,\log\log X}\right).
$$
Put $y=\exp(a\sqrt{L\log L})$ for fixed $a>0$. At fixed $p$, each potential supporting $\ell>y$ divides one of $O(\log R)$ outputs for at most $O(R/\ell)$ prime partners per output. Their total probability is at most $\exp(-(a-o(1))\sqrt{L\log L})$. If all supporting primes are at most $y$, a first-crossing $d$ is $y$-smooth and satisfies $T\le d<Ty<R$. Disregard the digit supports and count all such divisors of $p^2-1$ by the root bound. Apart from logarithmic factors, the resulting probability is bounded by
$$
 \sum_{\substack{d\ge T\\P^+(d)\le y}}\frac{2^{\omega(d)}}d
 \le T^{-\sigma}\prod_{\ell\le y}
       \left(1+\frac{2\ell^{-1+\sigma}}{1-\ell^{-1+\sigma}}\right).
$$
Choose $\sigma=((\log L)/2-2\log\log L)/\log y$. For large inputs this is positive and tends to zero. The product has log $O(y^\sigma\log\log y)=o(\sqrt{L\log L})$, by the prime harmonic sum. The exponent from $T^{-\sigma}$ is $-(t/(2a)-o(1))\sqrt{L\log L}$. Taking $a=\sqrt{t/2}$ balances the two bounds.

## Finite input checks

Enumeration verified both conditional counts at 18,622 combinations of $p,m,H$, including the reduced fractions and the divisibility of every remainder-slope determinant. There were three integral remainder slopes across these checks. A comparison with $R=2^{20}$, $p=1{,}863{,}671$, $H=4096$ and squarefree $m=25{,}878\mid p^2-1$ found 1,791 **composite** partners of the form $q=455j$ with $m\mid A$. The prime-partner bound here is below 164. An independent powered-gcd calculation verified the complete pooled part on all 724 ordered distinct prime pairs for $R=16,32,64,128$.

Each table row uses 128 seeded distinct products of primes with the indicated bit length. The threshold is $N^{1/4}$ at either factor. The middle cut uses $H=2^{\lfloor\operatorname{bitlength}(N)/2\rfloor}$. The first two columns choose one positive output. The last two pool the supports of all cuts and permit the best two bases or every base respectively.

| Bits per prime | Either half at middle cut | Any one half | Best two pooled bases | Full pool |
| ---: | ---: | ---: | ---: | ---: |
| 16 | 24/128 | 106/128 | 110/128 | 124/128 |
| 32 | 0/128 | 14/128 | 5/128 | 81/128 |
| 48 | 0/128 | 1/128 | 0/128 | 31/128 |
| 64 | 0/128 | 0/128 | 0/128 | 18/128 |

The positive counts measure coverage using the known synthetic factors. The unrestricted pool was still successful on 18 of the 128 largest test inputs. These finite results give neither a method to identify the relevant bases nor a success rate for a public factoring algorithm. [Check and sampling code](check_digit_order_support.py)
