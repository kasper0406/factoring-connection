# Which integer does the circle factor?

The branch already used the circle $x^2+y^2=a^2+1$ and its public point $(a,1)$. When the input to factor is $a^2+1=pq$ with distinct odd primes, finding a point outside the symmetry orbit of $(a,1)$ is equivalent to factoring that input. Asking whether such a point lies below a given height recovers it by binary search. This was the earlier predictor. [Circle-prefix equivalence](https://github.com/kasper0406/factoring-connection/commit/bf0ebd0) · [Earlier two-squares audit](../../notes/codex_r4_sos_classgroup.md) · [Integer and modular circles in the geometry survey](../../notes/geometry_results.md)

The variation reviewed here sets $a=N$ and asks for a factor of **$N$**. Even complete information about the representations of $N^2+1$ would need a way to transfer information to this different integer. For example, $493^2+1=139^2+473^2$, and $\gcd(139+473,493)=17$. How often can any short formula in the coordinates of a single point do this?

## Coverage on balanced semiprimes

Choose an unordered pair of distinct primes $p,q$ uniformly in $[R,2R]$ and put $N=pq$. Let $H=R^{o(1)}$, for example any fixed power of $\log R$. Supply **all** integer representations $x^2+y^2=A N^2+k$ for every $1\le A\le H$ and $0<|k|\le H$. Let $d\ge1$ be fixed. Allow all integer polynomials $F(x,y)$ of total degree at most $d$ and coefficients bounded in absolute value by $H$, even if the choice depends on the supplied points.

The fraction of inputs for which some point and some such $F$ give $1<\gcd(F(x,y),N)<N$ is at most

$$R^{-1/2+o(1)}.$$

This holds for affine tests such as $x+y$ and for any fixed degree. The bound already grants the points for free. It concerns a polynomial in the coordinates of one point. The cost of acquiring points and methods that combine information from several points are separate questions.

### Proof

Fix one hidden prime $r$ and count possible partners $s$ with $r\mid F(x,y)$ but $s\nmid F(x,y)$. Fix $A,k,F$. All estimates below are uniform for these choices. There are $R^{o(1)}$ choices in total since $d$ is fixed.

**Residues at the hidden prime.** Modulo $r$, a candidate satisfies $x^2+y^2=k$ and $F(x,y)=0$. The projective closure of the conic is smooth for large $R$ since $r$ is odd and $r\nmid k$. Unless $F$ is divisible by the conic modulo $r$, Bézout gives at most $2d$ residue pairs. In the divisible case, divide $F$ over the integers by $x^2+y^2-k$, monic in $y$. The remainder has degree at most one in $y$ and coefficients of size $R^{o(1)}$ because $d$ is fixed. Divisibility modulo $r$ then makes this remainder zero over the integers. On the supplied circle, $A N^2$ divides $F(x,y)$, which yields gcd $N$.

**Lifting a residue pair.** Fix representatives $(x_0,y_0)$ in $[0,r)$ and write $(x,y)=(x_0,y_0)+r(X,Y)$. Reduction of the norm equation modulo $r^2$ gives

$$2x_0 X+2y_0Y\equiv (k-x_0^2-y_0^2)/r\pmod r.$$

This congruence is independent of $s$. Its solutions form a translate of a lattice $\Lambda$ of index $r$, since $(x_0,y_0)\ne(0,0)\bmod r$. The coordinates $X,Y$ are $O(\sqrt A R+1)$. Let $w=(w_1,w_2)$ be a shortest nonzero vector of $\Lambda$. Pigeonholing a square grid of side $\lceil\sqrt r\rceil+1$ among the $r$ cosets gives $\|w\|=O(\sqrt r)$.

All values of $\det(w,(X,Y))$ belong to one residue class modulo $r$, because the determinant of any two vectors of $\Lambda$ is divisible by its index. Their absolute values are at most $R^{3/2+o(1)}$. There are therefore at most $R^{1/2+o(1)}$ lines to count. Each line fixes $V=w_1 y-w_2 x$ exactly. Set $D=w_1^2+w_2^2$ and $U=w_1x+w_2y$. On that line,

$$U^2-DA\,N^2=Dk-V^2.$$

**Nonzero right-hand side.** Write $K=Dk-V^2$. For $K\ne0$, the ternary quadratic form $U^2-DA\,N^2-KT^2$ and its binary part at $T=0$ are both nonsingular. With $T=1$, the uniform binary conic bound gives $O_\epsilon(B^\epsilon)$ integer pairs $(U,N)$ in a box of side $O(B)$, independent of the coefficients. Here $B=R^4$ suffices for all large $R$. This gives $R^{o(1)}$ partners per line. [Browning and Gorodnik, Lemma 4.1 and its application with $T=1$, p. 33](https://arxiv.org/pdf/1606.06342#page=33)

**Zero right-hand side.** There are at most two exceptional values of $V$. An integral solution requires $k>0$ and $DA,Dk$ both squares. Write $D=\delta e^2$, $A=\delta j^2$ and $k=\delta l^2$ with $\delta$ squarefree. Solving for the coordinates gives, for some signs $\sigma,\tau$,

$$e(x,y)=\sigma j N w+\tau l(-w_2,w_1).$$

Here $e=O(\sqrt R)$, so it is coprime to $N$. Let $P_0=\tau l(-w_2,w_1)/e$. This rational point has norm squared $k$. The integer $C=e^d F(P_0)$ satisfies $e^d F(x,y)\equiv C\pmod N$ and $|C|\le e^d R^{o(1)}$. If $C=0$, the gcd is $N$. Otherwise a proper hit at $r$ requires $r\mid C$.

Count the possible $C$ across **all** choices of $r$. The fixed $A$ determines its squarefree part $\delta$. For each $e=O(\sqrt R)$, the number of choices for $w$ is bounded by the number of representations of $\delta e^2$ as a sum of two squares. The two-squares formula bounds this by four times the number of positive divisors of $\delta e^2$, hence by $R^{o(1)}$. The vector $w$ and the sign of $P_0$ determine $C$ for fixed $A,k,F$. Each nonzero $C$ has only $O_d(1)$ prime divisors in $[R,2R]$ by its size. Thus at most $R^{1/2+o(1)}$ primes $r$ occur in a proper hit on any exceptional line, each with at most $O(R)$ partners.

The nonzero slices have at most $R^{1/2+o(1)}$ partners per prime. The zero slices have at most $R^{3/2+o(1)}$ ordered pairs altogether. Summing over the $R^{o(1)}$ choices of $A,k,F$ preserves the exponent. The number of unordered pairs of distinct primes in the interval is $R^{2+o(1)}$ by the prime number theorem, giving the stated fraction.

The nonzero-shift hypothesis matters. At $k=0,A=1$, if $s=u^2+v^2$ is prime with $s\equiv1\pmod4$, then $(r(u^2-v^2),2ruv)$ has norm squared $N^2$ and its first coordinate reveals $r$. This example concerns coverage with all points supplied.

### Relation to published work

The estimate per nonsingular conic above is published. The deduction here uses the particular congruence modulo $r^2$ and accounts for exceptional lines that could yield a proper gcd. Cao and Huang use $x^2+1=y_1^2+y_2^2$ itself as an example in Remark 1.8, concerning existence of representations. Browning and Verzobio study integral points on affine surfaces subject to polynomial congruences. Sica studies reconstructing factors using supplied factorizations of many neighboring integers. The branch records the coverage result as a restriction on a factoring proposal, with no claim of literature priority. [Cao and Huang, Remark 1.8](https://archive.mpim-bonn.mpg.de/id/eprint/4585/1/mpim-preprint_2021-27.pdf#page=8) · [Browning and Verzobio](https://arxiv.org/abs/2408.11453) · [Sica and the repository's prior-art entry](../../notes/prior_art_map.md#38-sica-factoring-with-hints-can-require-exponentially-many-hints)
