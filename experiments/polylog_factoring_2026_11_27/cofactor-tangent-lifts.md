# What survives the exact floor-quotient check?

Suppose $N=pq$ with distinct odd primes $p<q<2p$. Choose one or two binary positions above the units bit, write $M$ for their mask and $P=2^{1+\lfloor\log_2 M\rfloor}$ for their period. The dyadic binomial question asks whether both $p$ and $q$ have a one in $M$. Consider searching for the answer using odd, balanced integer candidates $u,v$ that obey
$$
v=\lfloor N/u\rfloor,\qquad uv\equiv N\pmod P.
$$
For every fixed $\delta>0$ and all sufficiently large $N$, these constraints admit **both answers** for every locally ambiguous mask with $P\le N^{1/4-\delta}$. Here locally ambiguous means that both answers occur on odd residue pairs with the specified product modulo $P$. Only a few singletons at the first four positions have a forced answer. The candidates we construct also have $\sqrt{N/2}<u<\sqrt N$ and $u<v<2u$. They need not be primes. An exact divisor requires the remainder $N-uv$ to vanish. This result bounds the deductions available from the displayed constraints. [Binomial identity](cofactor-bit-binomials.md) · [Local classification and a lift for every prescribed residue pair up to a sixth power](cofactor-bit-completions.md) · [Short-interval proof of the uniform exponent](cofactor-tangent-short-intervals.md)

## An integer point near a tangent line

For an odd integer $L$, consider $3x+2y=L$ together with $xy=N$ modulo $P$. There is a unique odd residue pair. Indeed $x\mapsto3x+2N/x$ permutes the odd residues, since its difference at $x,z$ is $(x-z)(3-2N/(xz))$ with an odd second factor. One can compute $x$ by Newton lifting the root $1$ modulo 2 of $3x^2-Lx+2N=0$, whose derivative $6x-L$ is odd, then set $y=N/x\bmod P$.

Take $L$ in an interval of length $R$ starting at $\lceil\sqrt{24N}\rceil$. On the real line $3u+2v=L$, the lower root of $uv=N$ is $\alpha=(L-\sqrt{L^2-24N})/6$. Among $u\equiv x\bmod P$, exactly one residue class modulo $2P$ makes $v=(L-3u)/2$ congruent to $y\bmod P$. Choose the largest such $u\le\alpha$. For $e=\alpha-u$ we have $0\le e<2P$ and
$$
N-uv=\tfrac12 e\sqrt{L^2-24N}+\tfrac32 e^2
       \ll P N^{1/4}\sqrt R+P^2
$$
when $1\le R\ll\sqrt N$. The remainder is divisible by $P$. If $P^2 R\le c\sqrt N$ for a sufficiently small absolute $c$ and $R=o(\sqrt N)$, then $u=(\sqrt{2/3}+o(1))\sqrt N$, $v/u\to 3/2$, and the remainder is less than $u$. This realizes the residue answer on two balanced integers with the exact floor quotient. Every line costs polynomially many bit operations in $\log N$. [Construction](cofactor_tangent_lifts.py)

## Complete sums give an elementary uniform range

Write $e_P(z)=\exp(2\pi i z/P)$. The indicator of either half of $\mathbf Z/P$ has constant Fourier coefficient $1/2$, all its remaining frequencies odd, and sum of absolute normalized coefficients $O(\log(2P))$, by the geometric sum. For a positive mask answer, it suffices that $x$ and $y$ both occupy the upper half at the largest selected bit. Count the pattern as zero on even $L$ and subtract $(1/4)\mathbf1_{L\ {\rm odd}}$. At each frequency $t$ the residual normalized Fourier coefficient is bounded by $O(P^{-1/2}\log^2(2P))$. To see this, expand the two halves. Any term outside the constant pair has an odd $h$ or $k$ (each is odd or zero) and leads, up to the sign convention for $t$, to
$$
\frac{1}{P}\sum_{x\bmod P,\ x\ {\rm odd}}
  e_P\big((h+3t)x+(k+2t)N/x\big).
$$
If $k$ is odd the reciprocal coefficient is odd. Otherwise $h$ is odd. Even $t$ then gives an odd first coefficient, while odd $t$ gives a reciprocal coefficient divisible by exactly 2. Thus at most one common factor of 2 needs to be removed. At the remaining modulus $2^s$, write $x=a+2^{\lceil s/2\rceil}w$ with $a$ odd. Summing in $w$ requires $A a^2\equiv B\bmod 2^{\lfloor s/2\rfloor}$ for the primitive coefficient pair $(A,B)$. This has at most four odd solutions modulo $2^{\lfloor s/2\rfloor}$, by the odd-square root count, giving $O(2^{\lceil s/2\rceil})$ for the unnormalized sum. Bounded moduli are absorbed in the constant. [General background on modular hyperbolas](https://arxiv.org/abs/1103.2879)

For a negative singleton answer it suffices that $x$ lies in its lower half, with baseline $(1/2)\mathbf1_{\rm odd}$. For two bits, require both corresponding bits of $x$ zero. If $B$ is the lower-half indicator at period $P$ and $C$ at the smaller period $P_0$, write $BC=(B-1/2)C+C/2$. All frequencies from $C$ lifted to period $P$ are even, so the first term has only odd frequencies and Fourier norm $O(\log^2(2P))$. The other term is a singleton at period $P_0$ and has exactly its baseline on a full cycle because the line map permutes odd residues. This gives the same bound for each component at its own period.

On any interval of $R\le P$ integer lines the unnormalized interval transform has sum of absolute coefficients $O(P\log(2P))$. The count thus differs from its odd-line baseline by $O(\sqrt P\log^3(2P))$. For a smaller-period component discard its full cycles first. A sufficiently large constant times $\sqrt P\log^3(2P)$ lines guarantees both answers for large $P$. Taking $P\le c' N^{1/5}/(\log N)^{6/5}$ suffices for the geometric lift. For finitely many bounded periods, a full cycle suffices for each locally feasible answer. The stronger short-interval result replaces the required length by $P^\varepsilon$ for each fixed $\varepsilon>0$, which gives the uniform exponent in the opening paragraph. Its constants depend on $\varepsilon$. [Short-interval argument](cofactor-tangent-short-intervals.md)

## Polylogarithmic search on almost all balanced prime pairs

For fixed odd $N\bmod P$, average the squared error of the pattern count around its odd-line baseline over the $P$ possible starting lines. For $R\le P$ the result is $O(R\log^4(2P))$. Indeed Parseval multiplies the squared uniform Fourier bound above by $\sum_t|\sum_{j=1}^R e_P(tj)|^2=PR$. A smaller-period component repeats evenly in this average. Its complete cycles have zero error, so the same argument applies to its remainder. Consequently at most an $O(\log^4(2P)/R)$ fraction of starting residues can miss either outcome when $R$ is sufficiently large.

Let $X=U^2$ and sample two distinct primes uniformly from $[U,2U]$. We need the *specific* start $\lceil\sqrt{24N}\rceil$. For integers $N\in[X,4X]$, each pair of residues
$$
(N\bmod P,\ \lceil\sqrt{24N}\rceil\bmod P)
$$
is attained $O((1+\sqrt X/P)^2)$ times. The second coordinate fixes $O(1+\sqrt X/P)$ possible integer ceilings. Each restricts $N$ to an interval of length $O(\sqrt X)$, in which the first coordinate allows $O(1+\sqrt X/P)$ integers. For $R\le P\le cX^{1/4}/\sqrt R$ this bounds the number of exceptional integers per mask by $O(X\log^4 X/R)$. There are $\asymp X/\log^2 X$ distinct prime products in this family by the prime number theorem and unique factorization. A union bound over $O(\log^2 X)$ masks gives exceptional proportion $O(\log^8 X/R)$. Periods below $R$ are handled by one whole cycle, provided each outcome is locally feasible.

Taking $R=\lceil(\log X)^{10}\rceil$ shows that outside an $O(1/\log^2 X)$ fraction of pairs **all** locally ambiguous singleton and pair questions with $P\le c X^{1/4}/(\log X)^5$ have both witnesses in the specified initial interval. Finding these witnesses then takes polynomial bit time per question. They explain why the floor relaxation cannot distinguish the answers even after enforcing the low product bits. They do not evaluate the binomial.

## Reproducible checks

The checker verified the inverse permutation in 9,936 direct comparisons and the window second moment against a separate autocorrelation calculation. At 80 seeded balanced semiprimes with 32 to 128 bits per prime, it tried a singleton and a pair per input. With highest bit $\lfloor\operatorname{bitlength}(N)/4\rfloor-5$, both answers appeared in all 160 trials within 22 odd lines. At the quarter-bit index itself, both appeared in 119 of 160 trials within a budget of 512 odd lines. No reported point had zero remainder. A separate run imposing $u\equiv1$, $v\equiv N\bmod 105$ found both answers in 146 of 160 trials at the lower index specified in the checker. That run is a numerical control for an extra sieve. The theorems use the power-of-two period alone. [Code and parameter choices](check_cofactor_tangent_lifts.py)
