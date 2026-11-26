# Two public mechanisms for the fixed-cap first moment

The [single-moment decoder](weighted_recurrence_single.py) needs the sum of the positive long coordinates $w$ modulo $N$ over $F(w,y)=w^2+\Delta y^2\le4NA$ with $4N\mid F(w,y)$. Two proposed ways to avoid scanning this disk are to construct a vanishing polynomial on its points or to telescope the weighted sum along rows. Both allow whole-disk batching. The following tests identify the particular information each proposal loses.

## Small-root certificates on the whole disk

Put $M=4N$. The ordinary integer shifts for multiplicity $m\ge1$ generate $I^m=(M,F)^m$ in $\mathbb Z[y,w]$. Every member has value divisible by $M^m$ on the required congruence. The proposed algorithm would find a nonzero member with absolute value strictly less than $M^m$ on a box containing the positive disk, so that every actual point would satisfy one integer equation. Standard lattice reduction proves the absolute bound by bounding scaled coefficients. Take arbitrary integer translations of both coordinates and keep the translated $F$ monic of degree two in $w$.

**Leading coefficient bound.** If $H\in I^m$ has degree $i$ in $w$, every coefficient in its leading $w$-coefficient is divisible by $M^{\max(0,m-\lfloor i/2\rfloor)}$. Repeated Euclidean division by monic $F$ gives the unique expansion $H=\sum_k (a_k(y)+w b_k(y))F^k$. Expanding each shift $M^{m-j}F^j S$ for $0\le j\le m$ in this basis shows that both coefficients at $k$ are multiples of $M^{\max(0,m-k)}$. The $w$-degrees of the basis monomials are distinct, so the highest $w$-coefficient has the claimed divisor.

For halfwidths $X\ge\sqrt M$ and $Y\ge1$, at least one scaled monomial coefficient has absolute value at least

$$M^{\max(0,m-\lfloor i/2\rfloor)} X^i\ \ge\ M^m.$$

This prevents the coefficient-norm certificate at every degree, even with exact shortest-vector search. The row $y=1$ extends to $\lfloor\sqrt{MA-\Delta}\rfloor$. A cover using only rectangles with $X<\sqrt M$ therefore needs $\Omega(\sqrt A)$ rectangles when $A$ grows on the near-square branch. The analogous univariate bound applies when a row is first specialized and the resulting polynomial must be nonzero. The repository’s [discriminant size estimate](../ellipse_disk_mod8/conductor_band.md#why-this-does-not-remove-a-power-of-n) gives $A=N^{1/4-o(1)}$ in its almost-all dyadic sense, so even these width subdivisions grow as $N^{1/8-o(1)}$. This is a cost for this cover. A method with prior row exclusions or additional congruences would require its own analysis.

**Integer grids.** Bounding the true maximum on integer points changes the constant in the width barrier for polynomials of polylogarithmic degree. Let $i>0$ and suppose a rectangular grid contains more distinct $y$ values than the degree in $y$ of the leading coefficient of $H$. At some grid row that coefficient is a nonzero integer divisible by $C=M^{\max(0,m-\lfloor i/2\rfloor)}$. Let $L$ be the difference between the final and first integer $w$ in the row. If $L\ge2i$, take the $i$-th finite difference with step $h=\lfloor L/i\rfloor$. Its absolute value is at least $C i! h^i$, and its coefficients have total absolute value $2^i$. The largest sampled absolute value is therefore at least $C(L/(4e))^i$, using $i!\ge(i/e)^i$. Width $L\ge4e\sqrt M$ precludes the strict bound here as well. For $i=0$, the nonzero leading value already has absolute value at least $M^m$. This elementary result concerns integer-coefficient shifts and ordinary rectangles. Compare the broader univariate capacity results and the repo’s partial-divisor analysis. [Chinburg et al.](https://arxiv.org/abs/1605.08065) · [Top-block analysis](../../notes/topblock_results.md)

## A direct test of first-moment telescoping

The full sum of the roots of $w^2\equiv c^2 y^2\pmod N$ over one period is zero modulo odd $N$, by $w\leftrightarrow N-w$. Could the weighted, sharply clipped partial sums have an antiderivative that uses only modular arithmetic? A function calculated from its endpoint by ring operations modulo $N$ and division by units takes equal values modulo $p\mid N$ at arguments congruent modulo $p$. Its row differences have the same property, including if one advances two steps to retain the prescribed parity. This holds at every degree and for coefficients depending on the two public inputs $N,y$.

On the existing balanced control $N=38{,}942{,}974{,}393$, the public parameters are $c=394679$, $\Delta=384531$, $A=358$. At $y=8499$ the cap allows $w\le5{,}290{,}605$. Compare two odd coordinates on this same row.

| $w$ | $F(w,y)\bmod4N$ | Mixed first-moment increment modulo $153247$ |
|---:|---:|---:|
| 2,038,973 | 0 | 46,762 |
| 1,732,479 | 90,247,771,288 | 0 |

The coordinates differ by twice the prime $153247\mid N$. Neither is on a public line. The first has norm level 205, and the second fails the common modulus test. The weighted increments differ modulo the prime although the inputs agree there. Consequently such an antiderivative cannot realize all these row prefixes, even after removing the principal contribution. This includes ordinary polynomial and unit-denominator rational telescoping. Integer representatives, floors, nonunit detection or an algorithm for only the fixed endpoint can operate outside this model. The full period's cancellation alone gives no partial row formula.

The [checks](native_cap_check.py) verify the row from the public form, the finite differences, and 96 leading-coefficient divisibilities, including cases with forced cancellation of high-degree terms. The factor in the row witness certifies a contradiction to the proposed formula. No efficient public evaluator follows from it.

```sh
python experiments/ellipse_aggregate/native_cap_check.py
```
