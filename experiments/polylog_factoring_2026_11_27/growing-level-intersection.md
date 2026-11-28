# Computing one intersection on a modular curve of growing level

Given $N\ge1$ and two integer matrices $g,h\in\Gamma_0(N)$, `level_intersection(N,g,h)` computes the signed algebraic intersection of their loops on the compact modular curve $X_0(N)$. The cost is polynomial in the bit lengths of the three inputs, without factoring $N$ or constructing a homology basis. For $N>1$, `unit_chord_reduction(N,g)` either returns a proper divisor of $N$ found among path denominators, or expresses the loop as a short sum of chords indexed by public units. The latter case identifies exactly which modular inverse comparisons the scalar measures. [Implementation](growing_level_intersection.py) · [Independent checks](check_growing_level_intersection.py)

Use row-major matrices $(a,b,c,d)$ of determinant one with $N\mid c$. A loop $[g]$ goes from $z$ to $gz$ with endpoints identified on the curve. Its class in absolute homology is independent of the basepoint, including a cusp. The quadrilateral joining two choices has opposite identified sides. For example, at level 11, $(7,-1,-55,8)$ and $(5,-1,-44,9)$ give intersection $+1$ in the complex orientation.

## Counting signed crossings without listing the edges of the cover

Set $S(z)=-1/z$, $T(z)=z+1$ and orient the base Farey edge $e$ from 0 to $\infty$. Express each matrix as a word in $S$ and powers of $T$ by nearest-integer Euclid in its first column. Removing $T^aS$ sends the absolute lower-left entry to at most half its previous value. Thus each word has $O(1+\log\|g\|)$ or $O(1+\log\|h\|)$ pieces, with the maximum entry norm used here. The upper triangular remainder is a single $T$ power in $\mathrm{PSL}_2$.

Realize $[g]$ as a path of Farey edges from $\infty$ to $g\infty$. A $T$ power fixes $\infty$. Each $S$ step with prefix $G$ contributes the single edge $-G e$. For the second loop choose a point inside the triangle with vertices $0,1,\infty$, and follow the successive images of this triangle along the word of $h$. Each adjacent pair permits a connecting arc crossing the common edge once. For a prefix $H$ their crossings and their signs against $H T^j e$ are

| Next piece | Crossed edge $H T^j e$ | Sign against that oriented edge |
|---|---|---:|
| $S$ | $j=0$ | $+1$ |
| $T^a$, $a>0$ | $1\le j\le a$ | $-1$ |
| $T^a$, $a<0$ | $a+1\le j\le0$ | $+1$ |

For instance, a positive translation crosses an upward oriented vertical edge from left to right. The intersection sign in the order (upward edge, rightward path) is negative. This accounts for the orientation of the dual walk. The primal edges contribute the further minus sign above.

The stabilizer of an oriented Farey edge in $\mathrm{PSL}_2(\mathbb Z)$ is trivial. Two matrices represent the same oriented edge upstairs modulo $\Gamma_0(N)$ precisely when their primitive bottom rows are projectively equal modulo $N$. For rows $(x,y)$ of $G$ and $(c,d)$ of $H$ the test for $G e$ to equal $H T^j e$ is therefore

$$xcj\equiv yc-xd\pmod N.$$

This is just the lower-left congruence for a matrix quotient. For the reverse orientation use the bottom row $(y,-x)$ of $GS$. Each congruence has zero solutions or one progression modulo $N/\gcd(xc,N)$, so an interval count costs an extended gcd and two divisions. Summing forward and reverse counts with the indicated signs gives the pairing.

The cusp at $H\infty$ has width $w=N/\gcd(c^2,N)$ because the lower-left entry of $HT^jH^{-1}$ is $-jc^2$. A full turn of width $w$ is parabolic and contracts after filling that cusp. Each $T$ interval can be shortened to its residual length modulo $w$. Its oriented edges then have distinct labels, so each fixed primal edge contributes at most one in absolute value on that piece. An $S$ step has the same bound. In particular

$$|[g]\mathbin{\cdot}[h]|\ll (1+\log\|g\|)(1+\log\|h\|).$$

The computed value is an intersection pairing of specified cycles. Determining the full genus or the trace of a growing Hecke operator still requires additional information. [Earlier global genus and triangulation analysis](../../notes/geometry_results.md) · [Full Heilbronn–Merel coefficient calculation](../../notes/cuspcoef_results.md)

**Cusps and elliptic points.** Choose the dual crossings at generic interior points of edges, away from elliptic fixed points, and smooth cusp endpoints of the primal cycle in neighborhoods disjoint from the compact dual walk. If a group element identifies a Farey edge with its reverse, the projected edge traverses the same arc in both directions. Both orientation tests hold and cancel. Order-three stabilizers occur inside triangles, where the connecting arcs can avoid them. Parabolic and finite-order elements give zero in the absolute homology of the compact surface. This argument uses closed matrix loops. Relative chains with uncancelled cusp boundary can pair nontrivially with a full cusp turn.

## What is left in the absence of a denominator factor

For $N>1$, test the gcd with $N$ of every cusp denominator on the primal path. If none gives a proper factor, every vertex lies at cusp 0 (unit denominator) or at cusp $\infty$ (denominator divisible by $N$) on the quotient. A primitive $a/c$ with $\gcd(c,N)=1$ is equivalent to 0 by Bézout, while $N\mid c$ identifies it with $\infty$. The cusp $\infty$ has just one incident unoriented edge orbit, since its neighbors upstairs are the integers and $T\in\Gamma_0(N)$ permutes them transitively. This edge joins $\infty$ to 0. Its two orientations have equal total multiplicities in the closed path and cancel.

Each remaining $-Ge$ has both entries $c,d$ of its bottom row invertible modulo $N$. Its row is equivalent to $(1,k)$, where $k=d/c\bmod N$. The matrix with rows $(0,-1),(1,k)$ sends $e$ to the edge from $-1/k$ to 0. Thus $-Ge$ gives the loop from 0 to $-1/k$. A matrix representing that loop is $g_k=(u,-1,1-uk,k)$ for $uk\equiv1\bmod N$. Applying this procedure to both inputs expresses the general answer as a sum of pairings between their returned indices. A denominator factor concerns the chosen path. Earlier work computed its frequency for one regular continued-fraction construction. [Denominator visibility, Proposition VIS](../../notes/classE_results.md)

At the cusp 0, use $S$ to lift a small horocycle to a horizontal line, with coordinate considered modulo the cusp width $N$. The two incident rays of the unit edge indexed by $a$ are at $a$ and $a^*=-a^{-1}\bmod N$. To see the second coordinate, carry $g_a0$ back by $g_a^{-1}$ before applying $S$. Distinct Farey edges have disjoint interiors in the tessellation. Their algebraic intersection can therefore be computed by smoothing the two loops near cusp 0 and counting a single signed chord crossing. For the open counterclockwise arc $A=(a^*,a)$ on this circle,

$$[g_a]\mathbin{\cdot}[g_b]=\mathbf1_A(b^*)-\mathbf1_A(b)$$

with $A$ empty for $a=a^*$. When endpoints coincide, the two edges coincide up to reversal or a chord is folded. The displayed open-arc formula gives zero in all these cases. For fixed $a$ and $b$ this needs modular inverses and cyclic order. The prime-level formula is Merel's chord lemma. The same cusp calculation gives this unit-edge version at composite level. General intersection formulas for modular symbols and Farey-symbol algorithms precede this implementation. The sparse input here is a pair of $\Gamma_0(N)$ matrices, with no preprocessing at their level. [Merel, Lemma 4](https://doi.org/10.1007/s002220050059) · [Farey-symbol pairing](https://doi.org/10.5802/jtnb.1143) · [Magma intersection algorithm](https://docs.magma-maths.org/ModularArithmeticGeometry/ModularSymbols/intersection_pairing.html)

### Conditioning a random unit edge on a hidden component

Let $N=pq$ with distinct primes. Fix $t$ public unit chords, independently of a uniform unit $b\bmod N$. Conditional on any $b\equiv\beta\bmod p$, the **joint law of all $t$ signed pairings** differs in total variation from a common public law by

$$O\!\left(t^2 q^{-1/2}\log^2(2q)\right).$$

For the reference law take two independent uniform points on the circle and use them as the ordered endpoints of the test chord. For one probe, if its arc occupies fraction $v$ of the circle, each sign $+1,-1$ has probability $v(1-v)$ and the rest gives zero. For a public matrix loop with a unit chord reduction, all its pairings against $g_b$ are obtained by summing the individual signed values. The same bound therefore applies with $t$ equal to the number of chords in all fixed reductions combined.

Indeed write $\gamma=b\bmod q$. On the unit torus, CRT expresses $(b/N,b^*/N)$ as a translation depending on $\beta$ of $(\bar p\gamma/q,-\bar p\gamma^{-1}/q)$, with $\bar p=p^{-1}\bmod q$. Each rectangle has discrepancy $O(q^{-1/2}\log^2(2q))$ from its area uniformly in the translation. To verify this, complete the interval indicators on the $q$-grid, whose normalized Fourier coefficients have $\ell^1$ norm $O(\log(2q))$. When both frequencies are nonzero, the complete prime Kloosterman sum has absolute value at most $2\sqrt q$. With just one nonzero frequency its magnitude is 1. The grid means differ from the arc lengths by $O(1/q)$. The $2t$ fixed endpoints cut the circle into $O(t)$ arcs, so all signs are constant on each of $O(t^2)$ rectangles off their boundaries. The probability of hitting a fixed endpoint is $O(1/q)$. Summing rectangle discrepancies proves the joint bound. It also bounds the difference between any two hidden rows up to a factor of two. This statement uses fixed public probes and a random unit test edge. It gives no corresponding bound for input-dependent choices of both paths. [Prime modular hyperbola estimates](https://arxiv.org/abs/1103.2879) · [Earlier unconditional inverse-endpoint bounds](../../notes/classE_results.md)

For the signed expectation alone the row error is smaller. If $[g]$ has $t$ unit chords, the scalar has the form $F(b^*)-F(b)$ with $F$ a sum of $t$ open-arc indicators. Its unconditional mean is exactly zero since $b\mapsto b^*$ is an involution. On a fixed hidden row each endpoint separately is uniform on a translated $q$-grid with one point removed. The counts of two such grids in any open arc differ by at most two, including when a boundary is a grid point. Therefore $|\mathbb E([g]\mathbin{\cdot}[g_b]\mid b\equiv\beta\bmod p)|\le 2t/(q-1)$.

## Checks

Run `python experiments/polylog_factoring_2026_11_27/check_growing_level_intersection.py`. The program recovers the absolute pairing 1 for Stein's integral basis of $X_0(11)$ and all 16 ordered entries in Magma's genus-two example at level 37. It checks the independent chord formula on 15,957 unit pairs across levels 2 through 47, bilinearity and skew symmetry in 216 general matrix cases, explicit elliptic and parabolic cycles, 4,097-bit cusp turns, and 25 other crossings at a 648-bit composite modulus. Of the general path tests, 165 returned a verified chord reduction and 51 returned an actual proper divisor at a denominator. Those small-level counts describe the fixtures. The size and hidden-component estimates above are uniform theorems. [Stein, Example 3.9](https://www.wstein.org/books/modform/modform/weight_two.html) · [Magma level-37 basis and pairings](https://docs.magma-maths.org/ModularArithmeticGeometry/ModularSymbols/intersection_pairing.html) · [Runnable checks](check_growing_level_intersection.py)
