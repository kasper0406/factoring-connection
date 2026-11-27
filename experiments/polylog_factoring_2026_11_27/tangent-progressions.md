# Tangent batching inside discriminant progressions

Combining modular-root progressions with rational-tangent blocks inside each shell gives an expected $N^{2/5+o(1)}$ bit bound for the public fixed-cap scalar on balanced semiprimes. The previous crossover applied the two methods separately and gave $N^{4/9+o(1)}$. The new evaluator is implemented, including the crossover decision. Its exponent concerns the scalar and remains larger than those of established general factoring algorithms. [Discriminant progression analysis](discriminant-progressions.md)

## A translated lattice in each shell

Write the shell as $x^2+D y^2=H$ as in the progression evaluator. Each root class restricts $x$ to $r+m\mathbb Z$. Split the ellipse into the halves $x^2\le H/2$ and $D y^2<H/2$. On each graph, write the independent and dependent physical coordinates as
$$X=x_0+h k,\qquad Z=z_0+\ell f(k),\qquad aX^2+cZ^2=H,$$
where $aX^2\le H/2$. Both the graph index $k$ and the target ordinate $f(k)$ must be integers. The steps $(h,\ell)$ are $(m,1)$ on the first graph and $(1,m)$ on the second.

On the graph,
$$|f''(k)|\le\frac{3a h^2}{\ell\sqrt{cH}}.$$
Use blocks of length $K=\max(1,\lfloor(\ell\lfloor\sqrt{cH}\rfloor/(4ah^2))^{1/3}\rfloor)$. For $K>1$, choose a continued-fraction slope $u/v$ with $v\le K$ and $|v f'(k_0)-u|\le 9/(8K)$ at the block start. The implementation computes the slope against a rational approximation obtained by scaling the square root before taking its integer floor. With the chosen scaling, the error in the exact derivative is at most $1/(8K^2)$, and the last convergent with denominator at most $K$ contributes at most $1/K$ to the multiplied slope error.

For any $k$ in the block, the change in $v f(k)-u k$ from its value at $k_0$ is bounded by $9/8+3/8<2$. Thus every integral point is on one of the five lines with intercepts from $\lfloor v f(k_0)-u k_0\rfloor-2$ through $\lfloor v f(k_0)-u k_0\rfloor+2$. Intersections solve one integer quadratic each. Check the block limits and the dependent coordinate progression. When $K=1$, check individual graph positions. The implementation uses integer arithmetic for the block construction and intersection tests. [Original tangent method](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_row_batches/arc_batches.py)

## Work at the public cap

For one class, divide each half's coordinate range by its block length. Both halves give $O(1+H^{1/3}D^{-1/6}m^{-1/3})$ blocks or direct positions. The additive term charges short and empty progressions. For example the first half has coordinate length $O(\sqrt H/m)$ and block length on the scale $(\sqrt{DH}/m^2)^{1/3}$. When the latter is bounded, the direct count is also bounded by the stated expression. The same calculation for the other half uses length $O(\sqrt{H/D})$ and block length on the scale $(m\sqrt H/D)^{1/3}$.

At shell $j$, the progression step is $m=D/s_j$ and there are at most $4^{\omega(D)}$ classes. The earlier bound $\sum_{j\le A}s_j\le A D^{o(1)}$ also bounds $\sum_{j\le A}s_j^{1/3}$. Thus all tangent work, root queries and trial division of $D$ together cost
$$N^{o(1)}\left(\sqrt\Delta+N^{1/3}\Delta^{1/6}\right).$$
Take tangent progressions if $\Delta^5<N^2$. Otherwise use direct progressions, with cost $N^{o(1)}(\sqrt\Delta+\sqrt N\,\Delta^{-1/4})$. Both main terms are $N^{2/5}$ at the crossover, and $\sqrt\Delta=O(N^{1/4})$. Randomness is confined to modular-root setup. This is a bound for computing the exact scalar on the public input.

## Checks

The translated graph code matched independent enumeration on 1,600 seeded graph queries. The exact scalar and divisor matched the direct progression evaluator on 560 small balanced semiprimes. The two methods also gave identical point lists on both historical auxiliary branches and on the 36-bit and 47-bit controls. A separate fixture exercised the actual tangent selection at $N=4{,}268{,}393$, $\Delta=148$, returning divisor 2,003 and trace 24,570.

The switch follows asymptotic bounds. On the 36-bit control the forced tangent variant took 122,270 line tests plus 4,230 direct positions. The hybrid selected direct progressions there. Run `python tangent_progressions.py 4268393` or `python tangent_progressions.py 38942974393` with SymPy available.
