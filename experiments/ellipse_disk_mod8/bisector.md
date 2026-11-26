# One norm shell, all slopes

At $B=1$ the 8,192-slope groups from the [centered recurrence](centered.md) can be evaluated together. The method recovers **every** positive mixed point of the parent ellipse, even outside the selected chart. It then answers slope-prefix counts by filtering this short output. On the largest chosen positive fixture it examined 48,688 polynomial values in the collision search; a direct square test used 8,388,610 rows. Setup, including factoring the public discriminant, is charged separately below.

For an eligible public branch with balanced $N=pq$, the deterministic bit bound is
$$\widetilde O\!\left(\Delta^{1/4}+\tau(\Delta)\bigl(1+(N/\Delta)^{1/4}\bigr)\right).$$
Here $\tau$ counts positive divisors and the tilde suppresses factors polynomial in the bit length. For $\Delta=N^{1/2+o(1)}$ this gives $N^{1/8+o(1)}$ for the norm-$N$ disk. On a positive instance the existing final gcd factors $N$ within the same bound. The separate coverage argument below shows why fast norm-$N$ evaluation alone cannot factor almost all balanced inputs on the $\ell=1$ near-square branch.

## The equal-norm identity

Write the public point as $(c,\ell)$ and another positive point as $(w,y)$; both have norm $4N$ in the metric $w^2+\Delta y^2$. Their sum is parallel to the angular bisector. Write the reduced ratio $(w+c)/(y+\ell)=s/t$, put $u=\gcd(s,\Delta)$, $v=\Delta/u$, $X=s/u$, $Y=t$ and
$$D=uX^2+vY^2,\qquad A=uX^2-vY^2.$$
Then $\gcd(uX,Y)=\gcd(X,vY)=1$. Reflecting the public point in the bisector, or expanding the equal-norm equation, gives
$$w=\frac{cA+2\Delta\ell XY}{D},\qquad y=\frac{2cXY-\ell A}{D}.$$

Two identities put short bounds on the search.

1. $D(\ell w+cy)=8NXY$. Since $\gcd(D,XY)=1$, we have $D\mid8N$. For each prime $h\mid N$, the same identity gives $h\mid D$ exactly when $h\nmid \ell w+cy$. A mixed point agrees with opposite public signs at the two primes, and its local vector is nonzero at both since its norm is $4N$. Hence precisely one of $p,q$ divides $D$. In particular $D\le 8\sqrt{2N}$. When $\Delta\equiv3\bmod8$, primitiveness of $(s,t)$ gives $v_2(s^2+\Delta t^2)\in\{0,2\}$ and $u$ is odd, so $4\sqrt{2N}$ suffices.
2. $(c-w,y-\ell)=2(cY-\ell uX)(vY,X)/D$. The two entries of $(vY,X)$ are coprime. Therefore $D\mid 2(cY-\ell uX)$. The odd hidden prime dividing $D$ makes $\ell uX-cY$ a nonunit modulo $N$.

The usual parametrization of a conic by a known rational point underlies the reflection. The two divisibilities are what shorten the search on this shell. [Rational parametrization background](https://www.math.gordon.edu/ntic/ntic2021/section-rational-conics.html)

## How the entire group is computed

For each public divisor $u$ of $\Delta$, let $K$ be the bound on $D$. Search the rectangle $1\le X\le\lfloor\sqrt{K/u}\rfloor$, $1\le Y\le\lfloor\sqrt{K/v}\rfloor$. Its area is at most $K/\sqrt{\Delta}$. The linear coefficients $\ell u,c$ are units modulo $N$. Normalize the longer coordinate so that the test is whether $i-\rho j$ shares a factor with $N$. Write $L,S$ for the long and short side lengths.

* Choose a power-of-two block width $k\asymp\sqrt{LS}$. Build $P_k(z)=\prod_{i=1}^k(z-i)$ over $\mathbb Z/N\mathbb Z$. Evaluate it at all $\rho j-b$ where $b$ is the start of a long-side block. This takes $O(S\lceil L/k\rceil)$ values. A unit value certifies that the block contributes no mixed bisector.
* Split every nonunit block in two and perform the next batched evaluation with $P_{k/2}$. A zero product can arise from different positions at $p$ and $q$; both children are tested. For sufficiently large balanced $N$, even the padded long side is $O(N^{1/4})<p$. Each short-side index has at most two nonunit positions, one per prime. Every level after the first has $O(S)$ active blocks, including padding. The implementation trims children past the rectangle and checks actual positions at the leaves.
* At each surviving pair perform the reflection with exact division. Check primitiveness, positivity, parity, the norm, and both public lines, and retain the resulting mixed point. The gcds in the search are used only as unit tests. No divisor of $N$ is used to construct a mixed CRT root. For the positive group, the downstream factoring gcd is available after verification.

Monic product and remainder trees work over this composite ring. With fast polynomial arithmetic, all stages for one rectangle take $\widetilde O(1+\sqrt{LS})$ bit operations, including building polynomials and testing individual leaves. The polynomials of a given width are reused across divisors. To obtain all $u$, the implementation deterministically factors $\Delta$ by a separate block product over possible divisors through the square root, recursively splitting the public integer. Width of order its fourth root bounds both setup and the final trial block. The factor tree has $O(\log\Delta)$ nodes, giving $\widetilde O(\Delta^{1/4})$ plus divisor enumeration. This gives the stated full bound. [Strassen product construction](https://web.maths.unsw.edu.au/~davidharvey/talks/factoring.pdf) · [FLINT multipoint evaluation](https://flintlib.org/doc/fmpz_mod_poly/_fmpz_mod_poly_evaluate_fmpz_vec.html)

Each positive $(w,y)$ gives four signed vectors on the parent shell. In the odd-$y$ chart use $x=(\pm w-(\Delta\bmod2)y)/2$ and reduce $x/y$ only modulo the public stride $2^k$. Each chart hit with positive $y$ represents its two simultaneous signs. Add the public point $(c,\ell)$ separately. Thus all residues, or both halves of any dyadic group, are available without traversing the fine slopes. This evaluates child groups supplied through their public parent branch and projective chart. An arbitrary form supplied without its parent need not have this data. The complete parent signed count modulo 16 is $4(1+\#\text{positive mixed points})\bmod16$.

## Opposite outcomes, actual work

| $N$ | Selected group mod 8 | Mixed slopes mod 65,536 in group | Group signed count mod 16 | All slopes signed count mod 16 | Polynomial values, public $\Delta$ / search | Direct $y$ rows, all slopes |
|---|---:|---|---:|---:|---:|---:|
| 33,119,086,993 | 2 | 12,122 | 4 | 8 | 21 / 1,461 | 6,391 |
| 38,942,974,393 | 3 | none | 2 | 4 | 59 / 362 | 636 |

For the first input the two children mod 16 separate cleanly. Residue 2 contains the public pair, and residue 10 contains the mixed pair. Its positive point has $w=53105$ and $y=6323$. The gcd of $(363973y-w)/2$ with $N$ is $206993$; the complementary factor is $160001$. The second input has no mixed point on this shell in any chart. The first fixture was deliberately chosen for a hit, so neither row estimates success frequency.

The [controls](bisector_check.py) check exact point sets on all 194 base branches, ten deliberate tie branches, 560 exhaustive small balanced N-only branches and 598 eligible even-$\ell$ branches. They also compare signed counts for every residue modulo 8 and its two lifts on these inputs by testing the original linear congruences on signed representations. There were 368 independent modular-window controls, including clipped final blocks and multiple nonunit locations, and 1,024 checks of public-discriminant divisor generation. Exactly 105 of the 194 base branches had a hit at $B=1$. The tie fixtures have their first mixed points at larger cutoffs.

Four further constructed positive instances, 47 to 91 bits, use primes $p,q=p+X+1$ and $\Delta=4p-X^2$ with $X$ odd, chosen so $\Delta$ is comparable to $\sqrt N$. The known and mixed representations come from the products of the two norm elements with real coordinates $X,X+2$ and imaginary coordinate 1, choosing the same or opposite signs. Only $N$ goes into the timed evaluator. The independent direct check evaluates the same full shell. On these small measured cases the faster implementation varied with the input. Counts of polynomial evaluations and square tests are different units; the raw data include polynomial construction and single-run wall times. [Results](bisector_results.json)

## How frequently can small cutoffs work?

Take distinct primes $p,q\in[P,2P]$ and $\ell=1$. Allow **all** positive $c$ with $4N=c^2+\Delta$ and $0<\Delta\le P^{1+o(1)}$. This includes the N-only choice and polynomially many in $\log N$ neighboring values of $c$. Even after pooling these branches, the proportion with a mixed point at $B=1$ is at most $P^{-1/2+o(1)}$, or $N^{-1/4+o(1)}$ on this scale.

Indeed, a positive mixed representation has $y>1$ and
$$(c-w)(c+w)=\Delta(y^2-1).$$
For given $\Delta,y$ the right side is positive and $O(P^2)$. It has $P^{o(1)}$ divisors, each specifying at most one pair $(c,w)$ and hence one $N$. The number of choices of $y$ is $O(P/\sqrt{\Delta})$. Summing over $\Delta$ bounds the number of inputs admitting a hit by $P^{3/2+o(1)}$. The denominator has order $P^2/\log^2 P$ by the [prime number theorem](https://dlmf.nist.gov/27.12). This counts actual mixed representations; a unit-test implementation that chooses to return an incidental factor could succeed for a different reason.

**Extension to early nonsquare shells.** For $1\le j\le B$, $B\le P^{O(1)}$, the equation becomes
$$w^2-jc^2=-\Delta(y^2-j).$$
The right side $M$ is nonzero at a mixed point: zero forces $j=k^2$, $y=k$, $w=kc$, the public multiple. For fixed $\Delta,j,y$, $|M|\le P^{O(1)}$. If $j$ is square, the same divisor-pair argument gives $P^{o(1)}$ choices of $(w,c)$, allowing either sign of $M$. Otherwise $w+c\sqrt j$ is an algebraic integer in the real quadratic field $\mathbb Q(\sqrt j)$ with norm $M$. In the maximal order there are at most $\tau(|M|)$ integral ideals of norm $|M|$: each prime exponent $e$ allows at most $e+1$ ideal choices by prime ideal factorization. Positive generators of the same ideal with the specified norm differ by a positive norm-one unit. Such a unit larger than 1 has integer trace at least 3 and hence size at least $(3+\sqrt5)/2$. All first embeddings here lie between 1 and $O(P\sqrt j)$. There are therefore $O(\log(Pj))$ possible generators per ideal, uniformly in $j$, also when $j$ has square factors. [Background on ideal factorization](https://www.math.iitb.ac.in/~dprasad/ant.pdf)

Summing $O(P\sqrt j/\sqrt\Delta)$ possible positive $y$ over $\Delta,j$ bounds the proportion of pairs with a mixed point through $B$ by
$$\min\{1,\ P^{-1/2+o(1)} B^{3/2}\}.$$
In particular, for fixed $\epsilon>0$ and $B\le N^{1/6-\epsilon}$ this proportion tends to zero, even allowing every near-square $c$ in the stated $\ell=1$ range. This is a bound on early witnesses in this family. It does not require an algorithm to visit every shell, impose a runtime lower bound, or settle coverage for other $\ell$ in the canonical grid.

As a diagnostic, a seeded sample of 64 distinct balanced pairs in each prime interval $[2^k,2^{k+1})$, $k=8,12,16,20$, yielded 42, 15, 8 and 3 first-branch hits at $B=1$ on the unadjusted floor choice (skip ineligible branches). Allowing that value of $c$ and its two predecessors gave 50, 20, 10 and 4. These finite samples illustrate the coverage issue; the asymptotic statement follows from the count above.

## Reproduce and continue

    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/bisector.py 33119086993
    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/bisector.py 38942974393
    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/bisector_check.py --output results-rerun.json

For a nonsquare shell $Nj$, the known public norm-$N$ point has the wrong length for reflection. Squares allow a rescaling, but the coverage bound still applies to any bounded search through early shells. The open target is a public aggregate of the varying equations $w^2-jc^2=-\Delta(y^2-j)$ at the common $c$, with all work for $j\le B$ charged. The two opposite-outcome slope fixtures now serve as solved baselines. A competitive end-to-end claim needs a mechanism that survives growing nonsquare cutoffs.
