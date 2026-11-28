# How predictable are binary multiplication carries?

For two random distinct $b$-bit primes, predicting the parity of every outgoing carry in the first $b-2$ free columns succeeds with probability at most $N^{-0.2075187+o(1)}$ if each prediction uses the product bits through its own column. The predictor can even inspect all earlier bits of both factors. Access to the whole low half of $N$ before making the first prediction gives a weaker upper bound of $N^{-0.1135460+o(1)}$. For every fixed $\delta<1/4$, both bounds remain exponential in $b$ when up to $\delta(b-2)$ wrong predictions are allowed. [Carry-polynomial reader](carry_polynomial.py) · [Exact finite checks](check_carry_prediction_bounds.py)

These estimates apply directly to the sparse changes used by the reader. For first or second differences, at least $(1/4-o(1))(b-2)$ changes in the free low columns are odd with probability tending to one. A uniformly short list of carry shapes cannot generally supply the information that the reader consumes.

## One column with past information

Start with independent uniform odd $b$-bit integers $P,Q$, put $R=2^{b-1}$ and $N=PQ$. Their leading and trailing bits are one. The pairs $(p_i,q_i)$ at the $m=b-2$ free columns $1\le i\le m$ are independent and uniform among four choices. Write $C_i$ for the outgoing carry at column $i$, with $C_{-1}=0$, and $n_i$ for the $i$-th product bit. The raw convolution at a free column is $w_i+p_i+q_i$, where $w_i=\sum_{j=1}^{i-1}p_jq_{i-j}$ is determined by earlier factor bits. Thus
$$
 n_i+2C_i=C_{i-1}+w_i+p_i+q_i.
$$
After seeing the earlier factor bits, $n_i$ tells the predictor whether $p_i=q_i$. Opposite bits give sum one. Equal bits leave $00$ and $11$ equiprobable, and the two possible outgoing carries differ by one. Before seeing $n_i$, the conditional error probability of any rule that will see $n_i$ is therefore at least $1/4$. This conditioning is on earlier factor bits and any independent coins, without future product bits. On the equal-bit case alone the error probability is $1/2$.

Let $E$ count the parity errors. In this model every past error is already measurable before the next bit pair is drawn. The successive errors therefore dominate independent Bernoulli trials of parameter $1/4$: at each history a uniform auxiliary coin couples a Bernoulli $1/4$ below the actual conditional error. In particular
$$
 \Pr(E\le h)\le\Pr(\operatorname{Bin}(m,1/4)\le h).
$$
The bound is exact for the stronger model with earlier factor bits supplied. Always predict the carry for sum one if the bits differ, or for $00$ if equal. Its errors are precisely the independent events $p_i=q_i=1$. At $h=0$ the bound is $(3/4)^m$. At $h=\lfloor\delta m\rfloor$, $0\le\delta<1/4$, the usual Chernoff bound gives $2^{-mD_2(\delta\Vert 1/4)}$, using binary relative entropy.

For any fixed integer polynomial $F$ with odd constant coefficient, predict $C_i\bmod 2$ by setting $(FC)_i=0\bmod2$ and using all earlier true carries. An error is exactly an odd coefficient of $FC$ in that column. Taking $F=(1-X)$ or $F=(1-X)^2$ proves the stated result for changes of carries or of slopes. It counts actual nonzero changes regardless of their magnitude. The same tail bound holds for a portfolio of these rules after multiplication by the number of distinct candidates. In particular, arbitrary lists of $2^{o(b)}$ fixed templates, or of causal predictors, preserve exponential decay for every fixed $\delta<1/4$.

## Seeing all low product bits at once

Now grant the entire word $u=N\bmod 2^{b-1}$ and, as each prediction is made, the earlier true factor bits. Given $u$, the free bits of $P$ remain independent and uniform. Indeed its odd residue modulo $2^{b-1}$ is uniform and for every such residue there is one compatible residue $Q=u/P\bmod 2^{b-1}$. At depth $i$, $q_i$ equals $p_i$ or its complement, with the choice determined before $p_i$ is drawn: the displayed carry equation and $n_i$ fix $p_i+q_i\bmod 2$. On an equal-bit occasion, toggling $p_i$ therefore toggles $C_i\bmod2$. Even after the full low word and the past are fixed, the parity error at such an occasion is a fair coin.

Let $K$ count equal-bit occasions. Unconditionally $K\sim\operatorname{Bin}(m,1/2)$, since the original bit pairs are independent. For any $1\le k\le m$,
$$
 \Pr(E\le h)\le
 \Pr(\operatorname{Bin}(m,1/2)<k)+
 \Pr(\operatorname{Bin}(k,1/2)\le h).
$$
To see the second term, fix $u$ and run the predictor. At an equal-bit occasion its next error is fair independently of all history through the preceding column. Extend the process after column $m$ with independent fair errors until there are at least $k$ occasions. The first $k$ marked errors are then independent fair bits, by iterated conditioning at these predictable occasions. On $K\ge k$, all of them are genuine errors or correct predictions within the original $m$ columns. The first term uses the unconditional law of $K$; its conditional law given $u$ can differ.

For $h=0$, take $k$ near $\alpha m$, where the solution of $\alpha=1-H_2(\alpha)$ in $(0,1/2)$ is $0.2270921952\ldots$, with $H_2$ the binary entropy. This bounds the success probability by $2^{-\alpha m+o(m)}$. For $0\le\delta<1/4$, choose $2\delta<\alpha<1/2$ instead. At $h=\lfloor\delta m\rfloor$, the two terms decay with respective base-two exponents per column at least $1-H_2(\alpha)$ and $\alpha[1-H_2(\delta/\alpha)]$. The one-step **unconditional** mean error is still at least $1/4$ at each free column for the odd-integer model. The probability of the equal-bit occasion is $1/2$ and its conditional error is $1/2$. Future low bits can correlate the locations of errors with past success, so the causal binomial domination above does not follow from this mean.

To pass either tail bound to distinct $b$-bit primes, condition on $P,Q$ both prime and unequal. This costs $O(b^2)$ by the prime number theorem on $[R,2R)$. Since $\log_2 N=2b+O(1)$, the exact causal and low-word bounds yield the exponents in the introduction. The tail estimates imply $\liminf \mathbb E[E]/m\ge1/4$ in the prime experiment. They do not assert a lower bound of $1/4$ at each specified column under prime conditioning. If a public rule receives $z$ additional bits of advice, union over its $2^z$ possible values. Thus $z=o(b)$ preserves these exponential bounds for fixed $\delta<1/4$.

## A high product word

A predictor may instead see the whole high word $h=\lfloor N/W\rfloor$ for a predetermined $W\asymp R$. Within the first $k-1$ free columns, with $k=\lfloor(\log_2R)/4\rfloor$, allow the high word as well as the usual causal information. The probability of perfect parity prediction there is at most $N^{-\log_2(4/3)/8+o(1)}=N^{-0.0518796+o(1)}$ on prime pairs. A fraction $\delta<1/4$ of errors **in this shorter block** still gives exponential decay.

Here is the counting step controlling dependence on the high word. Let $M=2^k\le R^{1/4}$. Uniformly in $h$ and residue classes $a,c\bmod M$, there are $O(R/M^2)$ integer pairs with $P,Q\in[R,2R)$, $P=a\bmod M$, $Q=c\bmod M$, and $hW\le PQ<(h+1)W$. Empty product ranges need no count. Otherwise $hW\asymp R^2$. Write $P=a+Mj$, $Q=c+Mv$ and $K_0=R/M$. The allowed integers $v$ lie in
$$
 \left[f(j), f(j)+O(1/M)\right],\qquad
 f(j)=\frac{hW}{M(a+Mj)}-\frac cM.
$$
The $j$ form an interval of $O(K_0)$ integers, and $f''(j)\asymp1/K_0$. For $M$ larger than a constant, there is at most one $v$ for each $j$. A Fejér kernel of degree $J-1$, with $J\asymp M$ chosen sufficiently small relative to $M$, is $\gg J$ whenever the distance from $f(j)$ to an integer is $O(1/M)$ within the displayed interval width. Its Fourier expansion bounds the number of valid $j$ by
$$
 O\!\left(\frac{K_0}{J}+\frac1J\sum_{s=1}^{J}
       \left|\sum_j \exp(2\pi i s f(j))\right|\right).
$$
The second derivative estimate of van der Corput gives $O(\sqrt{sK_0}+\sqrt{K_0/s})$ for the inner sum. In this setting it follows by splitting according to the distance of $(sf)'$ to an integer. With $\lambda=s/K_0$, the derivative is monotone with slope comparable to $\lambda$. There are $O(K_0\lambda+1)$ crossings or adjacent relevant integers. Portions within $\sqrt\lambda$ of these integers contain $O((K_0\lambda+1)(1+1/\sqrt\lambda))$ points, and on each complementary interval the first derivative bound gives $O(1/\sqrt\lambda)$. This yields the asserted estimate. The resulting occupancy is $O(K_0/M+\sqrt{K_0M})=O(R/M^2+\sqrt R)=O(R/M^2)$. Bounded $M$ is covered by the elementary $O(R)$ count in the same strip. [Graham and Kolesnik, Theorem 2.2](https://assets.cambridge.org/97805213/39278/excerpt/9780521339278_excerpt.pdf)

For each fixed $h$, the causal tree accepts at most $(M/2)^2\Pr(\operatorname{Bin}(k-1,1/4)\le d)$ ordered odd low-residue pairs with at most $d$ errors. Multiply by the occupancy bound, sum over $O(R)$ high words, and divide by the number of ordered prime pairs. This costs only $O(b^2)$ beyond the binomial probability.

## Scope and checks

An unrestricted table indexed by the **whole** product can return the carries of its prime factorization with zero error. A fixed-length program can compute the same table entries by exhaustive factorization if its runtime is unrestricted. The bounds concern limited product information or predeclared candidate families. They place no lower bound on the running time of an arbitrary full-input predictor.

For factor lengths three through nine, the check enumerates every ordered pair of odd integers and reproduces the exact causal binomial histogram, each column's $1/4$ error rate, and the first- and second-difference lower tails. For the low-word experiment, dynamic programming computes the actual optimum Hamming-budget success counts over every possible word and every allowed past-dependent policy. At nine bits the optimal causal rule made 2,187 perfect predictions out of 16,384 pairs. Full low-word lookahead raised the optimum to 2,634 on the same pairs. All finite binomial inequalities passed. Twelve exhaustive high-word occupancy tables had largest observed cell size at most $9R/(4M^2)$. [Reproduce the checks](check_carry_prediction_bounds.py)
