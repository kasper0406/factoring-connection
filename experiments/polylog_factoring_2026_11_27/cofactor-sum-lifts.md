# How many factor-sum bits does a dyadic state buy?

For an odd semiprime $N=pq$, one can search directly for a square along power-of-two progressions constructed from $N$. The method splits whenever the odd part of either $q-p$ or $q+p$ is at most $\sqrt{2B\sqrt N}$, with $B$ square tests per progression and $O(\log N)$ progressions. It needs no binomial answers or lattice. The hypothesis is restrictive. For two random primes in $[R,2R]$, even the full success set of this search has density $O(\sqrt{B/R}\log^2 R)$ for $1\le B\le R$. A second result below explains why a supplied factor residue and a short interval of possible sums do not extend the existing affine completion range. [Search and state congruence](cofactor_sum_lifts.py) · [Checks](check_cofactor_sum_lifts.py)

## Two square tests from the input alone

Write $A=(p+q)/2$, $C=(q-p)/2$ and $U=\sqrt N$, so $A^2-C^2=N$. For each modulus $M=4^j$, lift the at most four residues of $x^2=N\pmod M$ and test the first $B$ centers $x\ge U$ in each progression for whether $x^2-N$ is a square. Also lift $y^2=-N\pmod M$ and test its first $B$ nonnegative centers for whether $y^2+N$ is square. Include $M=1$ and increase $j$ through the bit length of $N$. A proper difference of the two square roots is a verified divisor.

The exact one-indexed positions of the true factorization in its fastest positive and negative progressions are
$$
1+\left\lfloor\frac{A-U}{2^{2v_2(C)}}\right\rfloor,
\qquad
1+\left\lfloor\frac{C}{2^{2v_2(A)}}\right\rfloor .
$$
Indeed the largest valid $M$ at the true positive center divides $C^2$, and at the negative center it divides $A^2$. Smaller valid moduli can only increase those ranks. For distinct prime factors the product has just one proper representation, so the smaller rank gives the exact budget required. If $d$ is the odd part of $q-p$, the first quotient equals $d^2/(A+U)$. If $d$ is the odd part of $q+p$, the second equals $d^2 C/(C^2+N)$. Both bounds are strictly below $d^2/(2U)$ for distinct primes. In particular $d^2\le2BU$ suffices in either sign.

For the density bound, fix $p$ in $[R,2R]$. A positive hit requires $q-p$ to lie in $2^j\mathbf Z$ and have magnitude $O(2^j\sqrt{BR})$ for some $j$, since $A+U=O(R)$. The number of possible $q$ per $j$ is $O(1+\min(R/2^j,\sqrt{BR}))$. A negative hit requires $q+p$ to lie in $2^j\mathbf Z$ and $q-p=O(B4^j)$, giving $O(1+\min(R/2^j,B2^j))$ candidates per $j$. The $O(\log R)$ levels give $O(\sqrt{BR}\log R)$ candidate integers in total. Divide by the number $\asymp R/\log R$ of prime choices, with an absolute constant for distinctness and ordering. At polynomial $B$ in $\log N$ the fraction tends to zero.

## Exact sum precision in a supplied state

Suppose one of the factors is known to equal the odd residue $a\bmod m$, $m=2^k$. On that state set $F(x)=x+N/x$ over the 2-adic units. Put $t=v_2(a^2-N)$, taking $v_2(0)=\infty$. The largest power of two modulo which **every** local value of $F(x)$ agrees with $F(a)$ is

| Case | Common modulus |
|---|---:|
| $t<k$ | $m2^t$ |
| $t=k$ | $2m^2$ |
| $t>k$ | $m^2$ |

For $x=a+mz$, the difference is $mz(a^2-N+amz)/(ax)$. The denominator is odd. In the first and last cases $z=1$ attains the asserted minimum valuation. In the middle case $z((a^2-N)/m+az)$ is always even, and $z=2$ gives valuation exactly one. All congruences use only $N,a,k$. The larger modulus on a diagonal state does not double the number of independently known factor bits. The reverse conversion from sum to factor residues and the loss near repeated roots were worked out earlier. [Prior inverse calculation, §4.1](../../notes/trackE_analysis.md) · [Earlier 2-adic shift calculation, §6.2](../../notes/signsplit_results.md)

There is also a ceiling on using other features of the exact local sum image. For $k\ge2$, that image contains a full 2-adic residue class of step at most $8m^2$. Choose $x_0$ from $a,a+m$ such that $t_0=v_2(x_0^2-N)\le k+1$. Such a choice exists since the two squares differ by a number of valuation $k+1$. On $x_0+4m\mathbf Z_2$ the derivative $1-N/x^2$ has constant valuation $t_0$. After subtracting $F(x_0)$ and dividing by $4m2^{t_0}$, the map on the free coordinate is bijective modulo every power of two. In fact its difference at two coordinates $z_1,z_2$ equals $z_1-z_2$ times an odd 2-adic unit, by the same difference identity. Thus it takes every residue at each precision. This proves the image contains the stated class.

If an interval of width $W$ contains at most $T$ locally admissible integer sums, then $W\le8m^2(T+1)$. For $x_1<x_2\le U$, the corresponding sum difference is $(x_2-x_1)(N/(x_1x_2)-1)\ge(x_2-x_1)^2/U$. The smaller-factor interval therefore has width at most $\sqrt{UW}$. For a balanced semiprime and a known unordered state, try both endpoints as residues of the smaller factor. Within each residue there are $O(1+\sqrt T)$ affine windows of radius $cN^{1/4}$ in the divided coordinate, for a fixed certified $c>0$. Thus the same data already permit polynomial affine factoring when $T$ is polynomial. In particular if a public interval containing the true sum excludes a hypothesized child modulo $2m$ through this full local sum test, $W\le32m^2$ and constant-many affine windows suffice from the parent state. [Certified affine reader](coupled-rank-windows.md)

## Checks and earlier coverage

The direct search passed 2,213 exact-budget comparisons on small semiprimes. In six planted cases the true signed progression required one check, with prime lengths 48, 96 and 160. The inputs were chosen so the prime gap was about $2^{3k/5}$ and the known-to-the-generator shared or opposite suffix had only $\lfloor k/3\rfloor$ bits for $k$-bit primes. At 160 bits the ordinary Fermat progression would require 906,921,140 and 881,814,417 checks respectively. The signed search received only $N$ and $B=1$. On 100 unplanted deterministic controls per size it split 17, 1 and 0 at prime lengths 16, 32 and 64. A further 280 direct 2-adic state enumerations agreed with the maximal common modulus in the table.

The earlier shared-suffix lattice screen already exploited $p\equiv q\pmod {2^k}$ once the common suffix by itself reached its certified threshold. The construction here trades a shorter suffix against an ordinary size bound or a small odd component and evaluates progressions by square testing. These are forms of Fermat sieving and partial-factor completion. The results do not supply a generic source of the first unknown binomial answer. [Shared-suffix screen](check_shared_suffix_split.py) · [Earlier signed-gap screen](check_sparse_suffix.py) · [Adaptive binomial questions](adaptive-binomial-lift.md)
