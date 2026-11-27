# Exact floor-square reciprocity at the fixed cap

The fixed-cap absolute moment has an exact formulation as a difference of two-dimensional floor-square sums. The formulation tests the proposed floor-sum and integer-carry route directly. Its literal size-reducing swap does **not** reduce the public parameters: it exchanges the two axes of the same ellipse, while the unit finite difference retains the exact square test. Direct evaluation takes

$$
\Theta\!\left(\sqrt{\frac N\Delta}\,A^{3/2}\right)
=\Theta(\sqrt N\,\Delta^{1/4})
$$

terms when $A=\lfloor\sqrt{\Delta/3}\rfloor$ grows. This is a cost of this exact expansion, not a lower bound for every algorithm. The tangent-progression evaluator remains strictly better. [Scalar decoder](../ellipse_aggregate/weighted_recurrence_single.md) · [Tangent evaluator](tangent-progressions.md)

## A unit finite difference selects the boundary

For an integer $t$, put

$$
G(t)=\binom{\lfloor\sqrt{\max(t,0)}\rfloor+1}{2},
$$

with $G(t)=0$ when $t<0$. At every integer $m\ge0$,

$$
G(m)-G(m-1)=
\begin{cases}
\sqrt m,&m\text{ is a positive square},\\
0,&\text{otherwise}.
\end{cases}
$$

Consequently the exact input-only residue before reduction modulo $N$ is

$$
T_N=\sum_{j=1}^{A}\ \sum_{1\le y<\sqrt{4Nj/\Delta}}
\bigl(G(4Nj-\Delta y^2)-G(4Nj-\Delta y^2-1)\bigr).
$$

The strict upper limit only avoids the zero coordinate. The difference is nonzero precisely when $4Nj-\Delta y^2=w^2$ for a positive integer $w$, and its value is then the required weight $w$. This formula preserves the integer carries exactly. It does not split the congruence modulo $N$ into its hidden linear branches.

Define the cumulative weighted ellipse count

$$
C(H)=\sum_{y\ge1}G(H-\Delta y^2).
$$

Expanding the triangular number and exchanging the two finite sums gives the reciprocal form

$$
C(H)=\sum_{w=1}^{\lfloor\sqrt{H-\Delta}\rfloor}
w\left\lfloor\sqrt{\frac{H-w^2}{\Delta}}\right\rfloor.
$$

Thus

$$
T_N=\sum_{j=1}^{A}\bigl(C(4Nj)-C(4Nj-1)\bigr).
$$

This is the exact analogue of swapping the axes in an ordinary two-dimensional floor sum. The new summation has the longer $w$ axis. Swapping again returns the first formula. No modulus or cutoff decreases, because the boundary is quadratic rather than a rational line. Evaluating the two large cumulative values separately also discards the useful cancellation: their difference is only the sparse boundary coefficient.

## Carry coordinates do not remove the square boundary

Let $c^2+\Delta=4N$ and set

$$L=(cy-w)/2.$$

The prescribed parity makes $L$ integral, and

$$
N\mid L(cy-L),\qquad
j=y^2-\frac{L(cy-L)}N.
$$

For fixed $y$, grouping $L$ by the quotient of the concave quadratic does describe the integer carries. Requiring zero remainder at a quotient boundary is exactly the condition that

$$w^2=4Nj-\Delta y^2$$

be a square. A group for each possible quotient therefore reproduces the $(j,y)$ terms above. The number of such pairs is

$$
R=\sum_{j\le A}\left\lfloor\sqrt{\frac{4Nj-1}{\Delta}}\right\rfloor
=\Theta\!\left(\sqrt{\frac N\Delta}A^{3/2}\right).
$$

For a generic near-square gap $\Delta=\Theta(\sqrt N)$, this is $\Theta(N^{5/8})$. A useful carry method must batch the **zero-remainder boundaries** across many quotient cells. Summing the cell populations or their endpoints by ordinary affine Dedekind reciprocity leaves that condition unresolved.

## Checks and decision

The implementation verified the cumulative swap and every unit boundary difference on 11,280 small $(H,\Delta)$ pairs. It matched an independent positive-coordinate disk enumeration on 253 balanced branches. On

$$N=38{,}942{,}974{,}393$$

it visited 2,879,887 floor-square terms, found the same 19 boundary points, and reproduced $T_N=69{,}529{,}082$. The tangent/progression implementation uses far fewer arithmetic tests on this fixture and has the better asymptotic bound.

This closes the literal two-dimensional floor-sum swap as a route to polynomial bit complexity. Keep a future proposal only if it supplies a new operation that contracts exact square boundaries across many $(j,y)$ cells. Integer-carry coupling remains admissible. Repackaging the boundary as the difference of two cumulative ellipse counts is insufficient.

Run from the repository root:

```sh
python experiments/polylog_factoring_2026_11_27/fixed_cap_carries.py
```
