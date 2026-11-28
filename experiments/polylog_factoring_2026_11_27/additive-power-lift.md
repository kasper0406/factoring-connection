# Reading successive digits of the additive power defect

Let $N=pq$ for distinct odd primes $p<q$ and put $F_N(x)=(x+1)^N-x^N-1$. Compute $F_N(x)\bmod N^j$ by modular powering. At each digit take the gcd with $N$. A proper divisor splits $N$. A gcd of $N$ lets us divide the residue by $N$ and read the next digit, up to the chosen precision. A gcd of 1 ends the test. For a rational seed with unit denominator the same calculation uses modular inversion. A proper denominator gcd already gives a factor.

The first digit is the scalar additive Frobenius test. The October record already identifies its exponent period, cites the single-prime Cauchy–Mirimanoff polynomial used below, and studies divided derivatives modulo $N^2$. The November scalar tests cover uniform seeds, the prescribed square-root seed and six rational symmetries. Here we divide the **value** by $N$ when the first gcd equals $N$. This value is publicly computable by modular powering. The canonical cyclotomic digit used by the earlier predictor has a different definition and still lacks a public evaluator. [October propositions B4 and D1](../../notes/trackBD_analysis.md) · [November scalar test](https://github.com/kasper0406/factoring-connection/commit/d61476b) · [Canonical cyclotomic digit](canonical-cyclotomic-lifts.md)

## Uniform seeds at a fixed modulus

Choose $x$ uniformly modulo $N$. If **any digit after the first** provides a split, then $N\mid F_N(x)$ and $r^2\mid F_N(x)$ for at least one prime $r\mid N$. The second condition admits a bound independent of the gap between the primes.

Write $s=N/r$ and $b=x^s\bmod r$. Modulo $r^2$, taking an $r$-th power depends only on the base modulo $r$. Thus $x^N=b^r\bmod r^2$. If the first digit vanishes at $r$, then $(x+1)^s=b+1\bmod r$, and $(x+1)^N=(b+1)^r\bmod r^2$. A prime-square hit therefore requires

$$ (b+1)^r-b^r-1=0\pmod{r^2}. $$

For $L_r(T)=\sum_{i=1}^{r-1}T^i/i$ over $\mathbb F_r$, binomial expansion gives $((b+1)^r-b^r-1)/r=-L_r(-b)\bmod r$. Mit’kin and Heath-Brown showed that each value of $L_r$ has $O(r^{2/3})$ preimages in $\mathbb F_r\setminus\{0,1\}$. Including those two points preserves the bound. [Mit’kin’s original paper](https://www.mathnet.ru/eng/mzm4630) · [Heath-Brown, Lemma 4](https://ora.ox.ac.uk/objects/uuid%3A455db511-6d3c-433c-b391-7b7dc22b723c/files/m73afbaf3523093bee3f76580d0c32d24)

* If $p\nmid q-1$, the map $x\mapsto x^s$ permutes $\mathbb F_r$ for both $r=p$ and $r=q$. A union bound gives probability $O(p^{-1/3}+q^{-1/3})$ for an additional split at any precision. This applies to every $p<q<2p$.
* If $p\mid q-1$, write $q=kp+1$. The necessary first-digit test at $q$ is $(x+1)^p-x^p-1=0\bmod q$, a polynomial of degree $p-1$ with leading coefficient $p\ne0\bmod q$. Its pass probability is at most $(p-1)/q$, already below $1/p$ when $k\ge p-1$. For $k<p-1$, reduction modulo $p-1$ gives $q\equiv k+1$. The corresponding first-digit test at $p$ is $(x+1)^{k+1}-x^{k+1}-1=0\bmod p$, whose degree is $k$ since $k+1<p$. The CRT components of the uniform seed are independent. Consequently the probability of passing both first-digit tests is at most $(k/p)((p-1)/q)<1/p$. This bounds all additional splits.

On every distinct odd semiprime the probability of a split after digit one is therefore $O(p^{-1/3})$, uniformly in the chosen precision. A success rate of at least $\delta$ for additional splits implies $p=O(\delta^{-3})$, within reach of trial division in time polynomial in $\log N$ and $1/\delta$. The probability is measured before checking the first digit. A deterministic choice of canonical seed depending on $N$ remains open.

The same argument applies to $x=P(t)/Q(t)$, for a fresh uniform $t\bmod N$, when the rational map has degree at most $D$ and remains nonconstant at each prime. Its coefficients may depend on $N$ and earlier trials. On regular inputs each local value has at most $D$ preimages, by counting roots of $P-yQ$. The nonexceptional bound gains a factor $D$ and the exceptional bound at most $D^2$. These are probabilities per draw, counting poles as missed draws. Checking a denominator gcd separately adds at most $D/p+D/q$ when $Q$ is a nonzero polynomial of degree at most $D$ at both primes.

## Deliberately lifting a seed that already passed twice

Suppose an integer $x$ satisfies $N^2\mid F_N(x)$. Put $U=F_N(x)/N^2$ and $G=(x+1)^{N-1}-x^{N-1}$. For $P(T)=F_N(x+NT)/N^2$, binomial expansion gives $P\in\mathbb Z[T]$ and $P(T)\equiv U+GT\bmod N$ coefficient by coefficient. The linear term follows from $F'_N(x)=NG$. The quadratic coefficient is divisible by $N$ because $N$ is odd and $N\mid\binom N2$. The higher coefficients contain a factor $N$ from the substitution.

A proper gcd of $G$ and $N$ already splits the input. If $G$ is a unit, a common root $T_*\bmod N^k$ of $P$ can be computed in time polynomial in $k\log N$. Start at $-U/G\bmod N$. At each stage correct the next digit with the inverse of $G\bmod N$, evaluating $P$ by modular powering before the exact division by $N^2$. For every $t$, the polynomial divided difference between $t$ and $T_*$ is $G\bmod N$ and hence a unit. At either prime, the valuations of $P(t)$ and $t-T_*$ consequently agree when capped at $k$. Searching among the translated seeds $x+Nt$ at this precision gives exactly the same factors at the same extra digits as applying the digit gcds to the public integer $t-T_*$. The case $G=0\bmod N$ lies outside this unit-slope calculation.

## Exhaustive small control

The control tried all 10,888 canonical residues across the 24 semiprimes formed from odd primes below 38 with $p<q<2p$. Each split is counted at the first digit that exposed it. [Reproducible check](check_additive_power_lift.py)

| First splitting digit | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| Seeds | 2,829 | 98 | 54 | 2 | 0 |

For every split at digits 2 through 5 in this cohort, at least one of $x$, $x+1$ and $x^2+x+1$ also gave a proper gcd after stripping common powers of $N$ from the nonzero integer.

A separate example was selected using known factors. For $N=8{,}293{,}697$ and $x=113{,}370$, the public computation gives $F_N(x)=7{,}765{,}881N\bmod N^2$ and $\gcd(7{,}765{,}881,N)=3881$. The other factor is $2137$. Here the three elementary gcds all return 1. This example shows the later digit can contribute additional information on a particular seed. Its selection gives no estimate for finding such a seed on a fresh input.

## How often the first divisibility test is reached

Choose two distinct primes uniformly without replacement from $[R,2R]$. At any precision, extracting a factor by successive gcds of the digits requires at least one prime to divide the original defect. This gives average bounds even for the first digit, complementary to the per-input bound above.

Fix one prime $r$ and a residue $u$ different from $0,-1$. Write $v=1+u$ and $L=\mathrm{lcm}(\mathrm{ord}_r u,\mathrm{ord}_r v)$. The partner prime $s$ passes at $r$ precisely when its class belongs to
$$S=\{e\bmod L:v^e-u^e=1\}.$$
For a nonzero difference $h\bmod L$, there is at most one ordered pair $e,e+h$ in $S$. Indeed, the two equations imply $(v^h-u^h)u^e=1-v^h$. A vanishing coefficient would force $u^h=v^h=1$. Otherwise both $u^e$ and $v^e=1+u^e$ are determined, which fixes $e\bmod L$. Hence $|S|(|S|-1)\le L-1$. Each class contains at most $1+R/L$ integers of the interval and there are $\asymp R/\log R$ eligible partner primes. Since $L\le 2R$, the conditional hit probability is $O(\log R/\sqrt L)$. This uses the same local orders identified in [Track B, Proposition B4](../../notes/trackBD_analysis.md#proposition-b4--where-the-frobenius-coincidence-sits-proved).

**Uniform seeds.** When $x$ is uniform modulo $N$, its residue at each prime is uniform regardless of the partner. The endpoints contribute $2/r$. For the other residues, $L\ge\mathrm{ord}_r u$, and
$$\sum_{u\in\mathbb F_r^\times}\frac{1}{\sqrt{\mathrm{ord}_r u}}
=\sum_{d\mid r-1}\frac{\varphi(d)}{\sqrt d}\le \tau(r-1)\sqrt r.$$
Here $\tau(m)$ counts positive divisors and is $m^{o(1)}$. Averaging over both primes and bounding the chance of a hit on either side gives $R^{-1/2+o(1)}$ per trial. Polynomially many uniform draws in the input bit length preserve a negligible mean success probability.

**Predetermined rational seeds.** Fix $x=a/b$ in lowest terms and let $H=\max(2,|a|,|b|,|a+b|)$. Endpoints give identically zero defect. Away from endpoints the homogeneous form $(a+b)^N-a^N-b^N$ is the negative sum of the odd powers of $a,b,-a-b$. We may permute the triple and reverse its sign to arrange $a,b>0$. These operations preserve its prime divisibility. The primes dividing the original product $ab(a+b)$ account for probability $O(\log(2H)/R)$, including any denominator gcd.

For $a\ne b$, the positive rationals $G=a/b$ and $J=(a+b)/b$ are multiplicatively independent. Indeed $a,b,a+b$ are pairwise coprime, and prime factors of at least two of them force both exponents in any relation $G^iJ^j=1$ to vanish. If the group they generate modulo $r$ has size $L\le Y$, two of the values $G^iJ^j$ with $0\le i,j\le\lceil\sqrt Y\rceil$ coincide modulo $r$. Thus $r$ divides the nonzero numerator of $G^hJ^k-1$ for one of $O(Y)$ nonzero pairs with $|h|,|k|\le\lceil\sqrt Y\rceil$. Each numerator has $O(\sqrt Y\log(2H))$ bits. Counting its prime divisors of size at least $R$ and dividing by the number of eligible primes gives
$$\Pr(L\le Y,\ r\nmid ab(a+b))\ll \min\{1,Y^{3/2}\log(2H)/R\}.$$
Here $L$ is indeed the group size since $\mathbb F_r^\times$ is cyclic. Summing the bound on dyadic ranges below $(R/\log(2H))^{2/3}$ and using $L^{-1/2}$ directly above that threshold gives mean $L^{-1/2}\ll (\log(2H)/R)^{1/3}$ over the regular primes, with the trivial bound when the threshold is below one. This is the elementary small-reduction argument for multiplicatively independent generators. [Murty and Petersen, Proposition 3.2](https://www.math.fsu.edu/~aluffi/archive/paper371.pdf#page=6)

When $a=b$ the reduced positive pair is $(1,1)$, representing the orbit of $1,-2,-1/2$. At odd $r$, a hit requires $s\equiv 1\bmod\mathrm{ord}_r 2$. Primes of order at most $Y$ divide $\prod_{d\le Y}(2^d-1)$, whose bit length is $O(Y^2)$. The same dyadic estimate with threshold $\sqrt R$ yields mean reciprocal order $O(R^{-1/2})$, giving hit probability $O(\log R/\sqrt R)$.

For any polynomial-size menu of rational seeds with polynomial bit length, fixed before drawing the primes (the menu can depend on $R$), a union bound gives at most $R^{-1/3+o(1)}$ mean success at any precision. Choosing adaptively *among those same seeds* satisfies the bound. A rule computing a fresh seed from the digits of the particular $N$ can range over many more values and remains outside this argument.

The exhaustive check also enumerates the complete exponent sets at every eligible residue for primes through 53. It checked 312 ordered nonzero differences across 349 sets and 4,886 prime-partner comparisons.
