# What the quartic tells us when the discriminant grows

At discriminant $-4027$, two or three quartic root-parity queries determine the hidden squared ideal class up to inversion. At discriminant $-1579$, two explicit semiprimes have identical answers to **all** unramified cubic-character queries and different hidden squared minima. The difference is the structure of the class group. Efficient public parity evaluation remains open in both examples.

The [discriminant −59 experiment](cubic59.md) used a cubic depending on the public branch coordinate $c$. Here we construct quartics from the discriminant's cubic fields. Their coefficients are independent of $c$ and $N$. Branch generation, finding and labelling fields, evaluating parity on the composite input, and extracting a divisor all have separate costs.

## A quartic for any supplied unramified cubic character

Let $K=\mathbb Q(\sqrt{-D})$ have fundamental discriminant $-D$, and let $4N=c^2+\ell^2 D$ be a primitive public branch, with $\gcd(c\ell D,N)=1$ and $(c+\ell\sqrt{-D})/2$ integral. Assume $N=pq$ for distinct odd primes. The public element has norm $N$. Its two prime ideal factors over $p,q$ therefore have inverse classes $C,C^{-1}$. An unramified cyclic cubic extension of $K$ gives a character to $\mathbb F_3$ whose values on these prime ideals are opposites. Equivalently its $S_3$ cubic over $\mathbb Q$ splits completely at both primes when the character vanishes, and is irreducible at both otherwise. For the examples below, an irreducible monic cubic of discriminant $-D$ with $D>3$ prime yields just such an extension: inertia at $D$ is a tame transposition and intersects $A_3$ trivially. [Class-field correspondence between cubic fields and index-three subgroups](https://msp.org/tunis/2020/2-2/tunis-v2-n2-p02-p.pdf)

Suppose a defining cubic is $F(Y)=Y^3+aY+b$ with $ab\ne0$. A general monic cubic can be depressed by $Y=3X+s$ where $s$ is its quadratic coefficient. Work away from $6ab\,\mathrm{disc}(F)$. For $\theta$ the image of $Y$, set $u=-\theta^3/b=1+(a/b)\theta$. Its norm is one and its three geometric components are distinct. Its trace is $3$ and $\mathrm{Tr}(u^{-1})=3+a^3/b^2$ by the elementary symmetric functions. The four norm-one square roots of $u$ over the algebraic closure have distinct traces, exactly the roots of

$$H_F(X)=b^2(X^4-6X^2-8X-3)-4a^3.$$

Indeed, $X=\mathrm{Tr}(v)$, $Y=\mathrm{Tr}(v^{-1})$ for $v^2=u$ give $X^2-2Y=\mathrm{Tr}(u)$ and $Y^2-2X=\mathrm{Tr}(u^{-1})$. Two allowed sign choices giving the same trace would force two components of $u$ to coincide. When $b=1$ one can instead take $u=-\theta$ and use the simpler $H_F=X^4-8X-4a$. In particular, $X^4-8X-8$ gives the same target bit on the discriminant $-59$ promise as the earlier quartic, subject to the good-prime restriction.

At an irreducible prime $r$, the norm-one group has odd order $r^2+r+1$. The quartic has exactly one root there. At a split prime the three coordinates of $u$ either admit four norm-one square roots in total or none. The bijection by trace commutes with Frobenius, so it counts rational roots in both cases. Consequently the composite root count is 1, 0 or 16 and its parity is precisely $[\chi(C)\ne0]$. A change of cubic polynomial or norm-one generator may change 0 into 16 without changing parity.

## Exactly which class information survives?

For $G=\mathrm{Cl}(K)$ every such bit factors through $V=G/G^3$. All of them together determine exactly $\{v,-v\}$ in $V$, where $v$ is the image of $C$. If two nonzero vectors span distinct lines, a linear functional vanishes on one and not the other. Conversely a linear functional vanishes on both $v$ and $-v$ or on neither. The same unordered image is obtained for the mixed class $C^2$, as multiplication by 2 is negation in $\mathbb F_3$.

For $r=\dim_{\mathbb F_3}V\ge1$, at most $2r-1$ adaptive queries suffice **given labelled fields representing the required functionals**. Test the $r$ basis coordinates for zero. If one is nonzero, choose its sign as $+1$. For each other nonzero coordinate, test whether its sum with this reference coordinate is zero. This fixes all relative signs. Rank zero gives no queries and no information. This query count does not bound the cost of finding or representing the fields.

At $D=4027$ enumeration gives nine classes and cubing every class gives the identity. The table uses four irreducible polynomials of discriminant $-4027$ and small split primes to label the four lines. The prime in each row splits completely for the corresponding polynomial and is inert for the other three. Thus the four kernels are distinct. [Independent arithmetic reference for the $X^3+10X+1$ field](https://bibliotekanauki.pl/articles/1393414.pdf)

| Cubic | Calibration prime | Nonprincipal reduced class in its kernel, up to inversion |
|---|---:|---|
| $X^3-44X+113$ | 13 | $[13,-9,79]$ |
| $X^3-8X+15$ | 17 | $[17,-11,61]$ |
| $X^3+X^2-7X+12$ | 19 | $[19,-1,53]$ |
| $X^3+10X+1$ | 29 | $[29,-27,41]$ |

The decoder tested the first two kernels and used the third when both characters were nonzero. It returned the correct hidden class up to inversion on **288 of 288** balanced inputs, covering all five possible outputs. It took two queries on 159 inputs and three on 129. The expected classes came independently from reducing $[p,c,q]$. All four quartic counts were evaluated on the **control side from the factors** to simulate the parity oracle. A fixed five-way class label alone gives no asymptotic factoring speedup under instance-local reductions, since all five guesses can be run and any returned divisor verified. [Prior instance-local bound](../../notes/tracedata_results.md)

At $D=1579$ the group is cyclic of order nine. There is just one kernel of index three. The polynomial $X^3+4X+7$ defines its cubic. The following two inputs therefore agree for every unramified cubic character at this discriminant.

| $N$ | $c$ | Quartic parity | Hidden squared class up to inversion | Minimum |
|---:|---:|---:|---|---:|
| 789,827 | 1,777 | 1 | $[5,-1,79]$ | 5 |
| 759,907 | 1,743 | 1 | $[17,-11,25]$ | 17 |

This rules out a fixed-discriminant decoder of the squared minimum using just the cubic-character answers. An algorithm using $N$ in other ways or obtaining observations at several discriminants has access to additional information.

## Coverage at growing discriminant

Suppose $P\le p<q<2P$ and $0\le2\sqrt N-c\le K$, $1\le K\le P$. Consider the direct catalogue of discriminants $D=4a^3+27b^2$ with $a\ge0$ and arbitrary integral $b$. Permit every positive conductor multiplier $\ell$ in $4N=c^2+\ell^2D$. The gap is at most $8PK$. For each $c,\ell$, the two coefficients have respectively $O((PK/\ell^2)^{1/3})$ and $O((PK/\ell^2)^{1/2})$ choices when the gap can contain a positive $D$. There are $O(P)$ possible $c$. Since $\sum_{\ell\ge1}\ell^{-5/3}$ converges and each tuple fixes $N$, the fraction covered among balanced pairs is

$$O\!\left(P^{-1/6}K^{5/6}(\log P)^2\right).$$

We used the prime number theorem only for the total number of pairs in the denominator. For $K$ polylogarithmic in $P$, even **all** such nonnegative-coefficient cubics with **all** conductors cover a vanishing fraction. To keep positive coverage in this model the allowed distance must be at least of order $P^{1/5}/(\log P)^{12/5}$. There is no positive coverage guarantee at that distance.

The same coefficient box applies to signed $a$ under a height restriction $4|a|^3+27b^2\le M\,i^2D$ where the polynomial discriminant is $-i^2D$. Summing over $i\le I$ costs a factor $O(M^{5/6}I^{8/3})$ in the bound. In particular polylogarithmic cancellation ratios and indices do not repair coverage for this method. Large cancellation, large indices, general binary cubic forms and branches far from $2\sqrt N$ are outside this estimate.

As a finite check, choose the largest odd $c\le2\sqrt N$ and all positive $\ell$, prime $D$, and irreducible positive-coefficient cubics of exact discriminant $-D$. Exhaustive balanced dyadic intervals gave:

| $P$ | Prime pairs | Covered | Distinct $D$ used | Bit 0 / bit 1 |
|---:|---:|---:|---:|---:|
| 256 | 903 | 18 | 4 | 7 / 11 |
| 1,024 | 9,316 | 145 | 18 | 53 / 92 |
| 4,096 | 107,416 | 1,303 | 44 | 470 / 833 |

The two local splitting types agreed on every hit. The quartic parity matched the prime cubic splitting test on all 2,352 distinct prime/polynomial checks across the three disjoint intervals. This measures one catalogue, independently of the repository's original conditional canonical-grid coverage claim. [Earlier canonical index and coverage hypothesis](../../notes/shared_collision_nov26.md)

**Implication for general factoring.** The cubic observables recover the full class in some fields and provably lose the required minimum in others. A route using them across fields needs an affordable way to supply labelled cubic fields on enough branches, a decoder connecting the accessible quotients to a divisor, and an unfactored-input parity evaluator. The three-query experiment isolates the last step on a rank-two control. None of these results gives a uniform factoring bound.

From the repository root run:

    python experiments/ellipse_disk_mod8/growing59.py --output rerun.json

[Implementation](growing59.py) · [Exact results](growing59_results.json)
