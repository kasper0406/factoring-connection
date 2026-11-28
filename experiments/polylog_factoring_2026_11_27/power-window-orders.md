# What a power exponent tells the affine interval reader

Let $N=pq$ with distinct odd primes $p<q\le 2p$. Fix a positive integer $e$ and a unit $a$, allowed to depend on $N$. Draw a fresh uniform unit $y$ and apply the earlier affine interval reader to $a y^e\bmod N$. A sufficiently concentrated power image can help the reader. For $e\le N^{1-\epsilon}$, any inverse-polynomial success on this input also lets us recover a divisor directly from $N,e$, for sufficiently large inputs and fixed $\epsilon>0$. The multiplier and the successful seed are unnecessary for that recovery. [Affine reader](coupled-rank-windows.md) · [Factoring from the exponent](factoring-with-an-order-hint.md)

The earlier power-sampling proof was for the decoder that checks whether a square residue has a short integer square root. Its forest count made one direct gcd competitive for that decoder. The affine reader checks proximity to a multiple of a hidden prime. Its set of successful residue pairs differs. [Earlier power-sampling proof](power-sampling.md)

## What can happen on one input

Write $G_r=(\mathbb F_r^*)^e$, $t_r=(r-1)/\gcd(e,r-1)$, and $X=\lfloor N^{1/4}/8\rfloor$. The local output is uniform on the coset $aG_r$. Let $m_r$ count the elements whose centered representative is in $[-X,X]$, and put $u_r=m_r/t_r$. The CRT makes the two draws independent. The product returned by the reader is a proper divisor with probability
$$
u_p(1-u_q)+u_q(1-u_p).
$$
The direct test $\gcd(y^e-1,N)$ has the same formula with local probabilities $1/t_p,1/t_q$.

Use the same public exponent $e=\operatorname{lcm}(1,\ldots,10)=2520$ and multiplier one in both examples. The sweep tests $\gcd((y^e)^j-1,N)$ for every $1\le j\le\operatorname{bitlength}(N)$.

| Selected prime pair | Radius | Window split probability | One gcd | Order sweep |
| --- | ---: | ---: | ---: | ---: |
| $587,593$ | 3 | $796/21682$ | $365/21682$ | $657/21682$ |
| $8191,8209$ | 11 | $466/1482$ | $125/1482$ | $1481/1482$ |

At the first pair, $y=32$ gives powered residues $-2\bmod 587$ and $152\bmod 593$. The reader returns 587. Both local orders exceed 19, so the entire sweep misses that seed. At the second pair the first image has order 13 and contains the four small elements $1,2,4,8$. The sweep splits on every image pair except $(1,1)$. These are deliberately selected inputs. The exact counts compare methods on each input. [Reproducible checks](check_power_window_cosets.py)

## What a frequently successful window implies

Fix one hidden prime $r$ and let $A$ be the signed integer representatives in the windowed coset. There are $m_r^2$ ordered products. Their residues lie in $a^2G_r$, a set of size $t_r$. Because $2X^2<r$, distinct integer products remain distinct modulo $r$. Each nonzero product $v$ occurs at most $2\tau(|v|)$ times, where $\tau$ counts positive divisors. Consequently
$$
u_r\le\sqrt{\frac{2\max_{1\le v\le X^2}\tau(v)}{t_r}}.
$$
For $X=0$ the local window is empty. The classical divisor bound makes the numerator $r^{o(1)}$, uniformly in the exponent and coset. Thus inverse-polynomial split probability forces $t_r=N^{o(1)}$ for at least one prime. In particular on any family where both image orders exceed a fixed power of $N$, every choice of unit multiplier has exponentially small success in the input bit length. Bourgain, Konyagin and Shparlinski study this same intersection using multiplicative independence and smooth numbers. [Cosets of small subgroups, Theorem 1](https://arxiv.org/pdf/1103.0567) · [Hardy and Wright, Theorem 315](https://blngcc.wordpress.com/wp-content/uploads/2008/11/hardy-wright-theory_of_numbers.pdf)

Constant local hit probability gives a tighter implication. If $u_r\ge\delta>0$, take the most numerous sign in $A$ and let $B$ be the absolute values of those representatives. Then $|B|\ge\delta t_r/2$ and $|BB|\le t_r$. Represent elements by their vectors of prime valuations, and let $s$ be the affine dimension of this finite set. Freiman's dimension lemma gives $|BB|\ge (s+1)|B|-s(s+1)/2\ge(s+1)|B|/2$. Projection onto a set of $s$ coordinates is injective on the affine span, and each coordinate is at most $\log_2 X$. Hence
$$
 t_r\le (2/\delta)(1+\log_2 X)^{4/\delta}.
$$
For fixed $\delta$ the image order is already polynomial in the input bit length, for any magnitude of $e$. Given an upper bound, the earlier result pads $e$ with the least common multiple of the integers up to the order bound and uses successive squaring to factor with constant probability. This handles the case where both images collapse as well. For an inverse-polynomial $\delta$, the dimension estimate alone gives a weaker conclusion. [Ruzsa on Freiman's dimension lemma, Theorem 1.1](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM2006.3/ICM2006.3.ocr.pdf) · [Earlier order-padding reduction](bounded-quadratic-roots.md#completing-an-exponent-that-covers-one-hidden-prime)

## Choosing an exponent

For a fixed nonzero polynomial exponent $e=|F(N)|$, the repository's earlier divisor count also explains why the promise rarely holds. Indeed $F(pq)\equiv F(q)\bmod(p-1)$. Holding $q\in[R,2R]$ fixed and counting $p$ in the same range, the sum of $\gcd(F(q),p-1)$ is $O(R\tau(|F(q)|))$, by summing over divisors of $F(q)$. For fixed $F$ this yields a mean kernel size of $R^{o(1)}$ as $R\to\infty$. Any local interval of width $O(\sqrt R)$ therefore has mean hit probability at most $R^{-1/2+o(1)}$, even allowing an input-dependent translation or unit scale fixed before the fresh draw. The remaining constructive question is how to select a useful exponent from the arithmetic of a general input. [Earlier polynomial-exponent bound](bounded-quadratic-roots.md#fixed-polynomial-exponents-on-separated-balanced-primes)

The accompanying work gives an explicit threshold for recovery from an order hint and tests three proposals constructed from the digits of an unfactored input. [Proof and constructor tests](factoring-with-an-order-hint.md)
