# Small first vectors for the 72-frame construction

For $N=pq$ with distinct odd primes $p,q\in[R,2R]$, scan $x=0,1,\ldots$, accepting the first $x$ with $(-(x^2+1)/N)=-1$, and return any proper gcd encountered first. **For every fixed $D>1$, at most $R^{2/D+o(1)}$ products require searching beyond $(\log R)^D$.** For every fixed $p$, at most $R^{1/D+o(1)}$ partner primes $q$ require it. The mean number of probes over ordinary balanced pairs remains bounded even without a cutoff. These unconditional estimates make the entire list of first vectors in the mixed ternary Gram accessible in polynomial bit time on all but an $R^{-2+2/D+o(1)}$ fraction of pairs. One can use the public cutoff $\lceil(\log N)^D\rceil$. The earlier uniform-residue sampler finds a valid Jacobi sign quickly but can give a large first norm. [Original 72-frame result](mixed-ternary-embeddings.md) · [Public selector and first-vector enumeration](short_jacobi_grams.py)

### The character tests forced by a failed scan

Write $s_x=-(x^2+1)$. Work modulo rational squares and suppose $Y\ge3$. The product $s_1s_2s_3=-100$ gives the squareclass of $-1$. Together with $s_1$ this also gives $2$. Let $\ell\le2Y+1$ be prime with $\ell\equiv1\bmod4$. A root $1\le x\le(\ell-1)/2$ of $x^2=-1\bmod\ell$ gives $x^2+1=\ell m$ with $m<\ell$. Every odd prime of $m$ is $1\bmod4$ and smaller than $\ell$. Induction shows that all these $\ell$ belong to the squareclass span of the tested values. For the polylogarithmic cutoff $Y$, both hidden primes exceed $Y^2+1$ for large $R$, so all tests are units. A failed scan forces $(\ell/p)=(\ell/q)$ on this whole prime set.

Set $Y=\lfloor(\log R)^D\rfloor$. Choose $k=\lfloor\log R/\log Y\rfloor$, and let $\mathcal S$ consist of all products of $k$ distinct primes congruent to 1 modulo four from $[Y/2,Y]$. These labels are odd, squarefree and at most $R$. By the prime number theorem in the fixed progression, there are $\asymp Y/\log Y$ available primes. Taking the binomial coefficient gives $|\mathcal S|=R^{1-1/D+o(1)}$.

Apply Heath-Brown's quadratic large sieve with $a_d=(d/p)$ for $d\in\mathcal S$ and zero elsewhere. With stars restricting both indices to odd squarefree integers, Theorem 1 gives
$$
\sum_{q\le 2R}^{*}\left|\sum_{d\le R}^{*} a_d(d/q)\right|^2
\ll_\epsilon R^{1+2\epsilon}|\mathcal S|.
$$
For an exceptional prime $q$, each summand inside the absolute value is $+1$. Dividing by $|\mathcal S|^2$ and allowing arbitrary fixed $\epsilon>0$ proves the stated count uniformly in $p$. This proof uses the positive-integer labels only for analysis. The public scan computes ordinary Jacobi symbols directly from $N$. [Heath-Brown, Theorem 1](https://doi.org/10.4064/aa-72-3-235-275)

For the count over all products, double the allowed logarithmic label size. Take $k=\lfloor2\log R/\log Y\rfloor$, so $d\le R^2$ and $|\mathcal S|=R^{2-2/D+o(1)}$. Use $a_d=1$ for $d\in\mathcal S$, with outer odd squarefree index $m\le4R^2$. Each failing $m=pq$ contributes $|\mathcal S|^2$ since $(d/N)=1$ on the entire label set. The sieve upper bound $O_\epsilon(R^{2+4\epsilon}|\mathcal S|)$ gives at most $R^{2/D+o(1)}$ such distinct products. Counting ordered prime pairs changes this by at most two.

### Smaller cutoffs and worst cases

There is also an exponential tail for very small $Y$. At least $(1/2-o(1))Y$ of the values $x^2+1,\ 1\le x\le Y$, have a prime factor greater than $2Y$. To see this, the logarithm of their product is $2Y\log Y+O(Y)$. An odd prime $\ell\le2Y$ that occurs is $1\bmod4$. At most $2Y/\ell^e+2$ terms are divisible by $\ell^e$, by the two lifted roots. Summing over powers $\ell^e\le Y^2+1$ and such primes, their contribution is at most
$$
2Y\sum_{\ell\le2Y,\;\ell\equiv1(4)}\frac{\log\ell}{\ell-1}+O(\pi(2Y)\log Y)
= Y\log Y+o(Y\log Y).
$$
The prime 2 contributes $O(Y)$. Prime factors above $2Y$ occur to exponent one and in at most one row, since $\ell\mid(x-z)(x+z)$ would force $x=z$. Each row contributes at most $\log(Y^2+1)$ to the remaining logarithmic mass. The rows with such a prime are therefore linearly independent squareclasses.

Everest and Harman's Theorem 1.4 shows that for all sufficiently large $Y$, more than $0.5324Y$ rows have a primitive prime divisor that has not occurred earlier. For $x>1$ their Proposition 1.3 shows this prime exceeds $2x$, hence has valuation one. The row at $x=1$ has the new prime 2 with valuation one. Triangular independence gives the rank bound $r_Y>0.5324Y$, even if a prime divides later rows. [Everest–Harman](https://arxiv.org/pdf/math/0701234)

Here $r_Y$ is the squareclass rank of all tested values. The associated quadratic characters have image of size $2^{r_Y}$ on units of the common modulus $M=8\prod_{\ell\text{ odd},\,\ell\mid\prod_{x\le Y}(x^2+1)}\ell$. This follows from the primitive quadratic character of every nontrivial squareclass. For fixed $p$, a failed scan allows one such sign pattern for $q$, occupying $\varphi(M)/2^{r_Y}$ residue classes. For fixed $Y$, the limiting fraction of partners is exactly $2^{-r_Y}$ by the prime number theorem in progressions. Uniformly when $Y\log Y\le(1/2-\delta)\log R$, we have $M\le R^{1-\delta}$ for large $R$. Brun–Titchmarsh in each class bounds the fraction by $O_\delta(2^{-r_Y})$. A sufficiently large multiple of $\log\log R$ gives any prescribed inverse-logarithmic failure rate. [Montgomery–Vaughan](https://doi.org/10.1112/S0025579300004708)

No fixed cutoff works universally. Let $P$ be the number of primes in $[R,2R]$ and take $Y=\lfloor\log_2(P-1)\rfloor$. For large $R$ all factors exceed $Y^2+1$. Among the $P$ primes two agree on all $Y$ signs $s_1,\ldots,s_Y$, hence also on $s_0=-1$. Their product misses the entire scan with $Y=(1-o(1))\log_2 R$.

An unconditional upper bound for either an accepted sign or a proper gcd is $x\le(p+q-2)/2$. Otherwise the periodic local words $x\mapsto(s_x/p)$ and $x\mapsto(s_x/q)$ coincide throughout the symmetric interval from $-(p+q-2)/2$ to $(p+q-2)/2$. This includes zeros because a zero at only one prime yields a proper gcd. The common word has periods $p,q$ and length $p+q-1$. Fine–Wilf forces period one. Each local word is nonconstant, since its period sum is $-(-1/p)$ or $-(-1/q)$, obtained by counting $(y-x)(y+x)=1$. For $p\le q\le2p$, this gives a bound of $3\sqrt N/(2\sqrt2)-1$. Under GRH, effective Chebotarev in $\mathbb Q(i,\sqrt N)$ instead supplies an unramified $\ell\equiv1\bmod4$ with $(N/\ell)=-1$ of size $O(\log^2 N)$. The discriminant has logarithm $O(\log N)$. Reciprocity and the squareclass induction then give a polylogarithmic scan for every pair under that hypothesis. [Fine–Wilf](https://doi.org/10.1090/S0002-9939-1965-0174934-9) · [Lagarias–Montgomery–Odlyzko](https://doi.org/10.1007/BF01390234)

For uniformly chosen distinct primes in $[R,2R]$, the limiting mean number of probes, including the successful one, is
$$
1+\sum_{Y\ge0}2^{-r_Y}.
$$
The finite sum through $Y=120$ is 2.15750158371793 to the precision shown, with $r_{120}=86$. For each fixed $Y$ the summand is the tail probability above. To justify summing the limits, fix $D>2$. Up to $K\log\log R$ the tails are bounded by a summable geometric sequence. Between there and $(\log R)^D$ their sum is $o(1)$ by taking $K$ large. Beyond this cutoff the contribution is $O(R)\,R^{-2+2/D+o(1)}=o(1)$ by the all-products count and the periodicity bound. The same argument, choosing $D>2/(2-\alpha)$, gives convergence of every fixed moment of order $0<\alpha<2$ of the stopping index.

### What remains after the first vector is known

For the selected $x$, put $A=x^2+1$. Enumerate the nonnegative triples $i\le j\le k$ of squared norm $A$ by looping over $i\le\sqrt{A/3}$ and $j\le\sqrt{(A-i^2)/2}$ and testing whether $A-i^2-j^2$ is a square. This takes $O(A)$ integer-square tests. Signed permutations give every possible first vector $u$.

For an uncapped algorithm, a faster theoretical enumeration gives polynomial expected bit cost averaged over both the balanced input distribution and the algorithm's random bits. For each of the $O(1+x)$ first coordinates $i$, factor the positive auxiliary integer $m=A-i^2$ in $m^{1/5+o(1)}$ deterministic bit operations. For each prime $\ell\equiv1\bmod4$ of $m$, a uniformly random nonzero residue $b$ is a nonresidue with probability one half. In this case $b^{(\ell-1)/4}$ squares to $-1$ modulo $\ell$. A Gaussian Euclidean gcd with $\ell$ finds a Gaussian prime above $\ell$. Unique factorization in $\mathbb Z[i]$ now lists all signed two-square representations, by distributing split prime powers between conjugates, using the ramified factor $1+i$, and requiring inert exponents even. There are at most $4\tau(m)$ choices. Filtering for ordered coordinates yields the required triples. Thus the expected enumeration cost at first index $x$ is $(1+x)^{7/5+o(1)}$, apart from polynomial factors in bit length. The bounded stopping moments apply. The supplied code uses the shorter nested square test for the cutoff guarantee. [Harvey–Hittmeir, deterministic integer factorisation](https://doi.org/10.1090/mcom/3708)

Each remaining search can be put in one explicit binary norm form. Put $u_0=(x,1,0)$. For $u\ne u_0$, set $w=u-u_0,\ d=A-u\cdot u_0$ and $T=dI-ww^\top$. Then $1\le d\le2A,\ T^2=d^2I,\ Tu=du_0$. For $u=u_0$, take $d=1,T=I$. When $2A<\min(p,q)$, $d$ is a unit modulo $N$, and $Tv$ preserves the common gcd of the coordinates with $N$. This size condition holds in the polylogarithmic regime above. The transformed coordinates have the form
$$
Tv=(-dN+t,-xt,z),\qquad z^2=t(2dN-At).
$$
If $z=0$, the vector is divisible by $N$, since $A$ is a unit modulo $N$. Otherwise write $z/t=h/s$ in lowest terms. Then $F=h^2+A s^2\mid2dN$, since $t=2dNs^2/F$ and $\gcd(F,s)=1$. The inert prime for $-A$ cannot divide $F$. For a mixed frame the split prime must divide $F$, since $F\mid 2d$ would make every transformed coordinate divisible by $N$. Thus a mixed frame gives a primitive norm value $F$ equal to the split prime times a divisor of $2d$. Conversely any such norm value immediately reveals the prime through $\gcd(F,N)$. This calculation identifies a public finite menu of multipliers. It leaves the task of finding the binary coordinates producing the norm. At $N=35,x=4$, the noncanonical first vector $(2,-3,-2)$ with partner $(-10,30,15)$ gives $d=12,\ F=(-2)^2+17\cdot1^2=21$, and the gcd is 7.

### Checks on ordinary inputs

Four held-out seeds generated 256 balanced semiprimes in each of the 48-, 64-, 96- and 128-bit tiers. The 512 primes in each tier are distinct. No gcd shortcut occurred. The table counts the first accepted $x$ across all 1,024 moduli.

| $x$ | 0 | 1 | 2 | 4 | 5 | 6 | 9 | 10 | 12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Moduli | 511 | 263 | 131 | 56 | 22 | 22 | 12 | 4 | 3 |

The construction also found a balanced pair $73{,}477\cdot80{,}557=5{,}919{,}086{,}689$ whose first accepted $x$ is 27. Other controls check 1,744 squareclass inclusions, 16 equal-size sign classes on 256 units, 1,218 selector comparisons, the periodicity bound on every distinct pair of odd primes below 120, and 13 complete signed first-vector orbits. [Reproducible checks](check_short_jacobi_grams.py)
