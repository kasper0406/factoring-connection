# Linear query count for affine mixed root isolation

The affine isolation reduction needs only $O(\log N)$ nonadaptive weighted queries for constant success probability. The earlier bound was $O((\log N)^2)$. Averaging over the CRT sign patterns gives an absolute constant number of trials per dyadic width. The explicit constant below is conservative, so the implementation uses the earlier small-input guarantee when it requires fewer trials. [Weighted disk reduction](affine-isolation.md)

Let $k=\omega(N)$, $Q=2^k$, and $R=Q-2\ge Q/2$ be the mixed-root count. For a random unit multiplier and a uniform translation, let $X$ count the mixed roots in a window of width $w$. Put $f=w/N$, so $\mathbb E X=Rf$.

## Average over sign differences

An ordered pair of distinct roots, including the two public roots in this upper bound, differs at some nonempty subset of the prime-power components of $N$. Let $m$ be the product of those components. Each subset occurs for exactly $Q$ ordered pairs because there is one partner for every first root. The difference after uniform unit multiplication has the form $(N/m)v$, with $v$ uniform over the units modulo $m$.

The two images can both land in the window only if the difference has circular distance less than $w$. If $m<1/f$, this is impossible. Otherwise the number of possible units $v$ at either end of the period is at most
$$2\left(f\varphi(m)+2^{\omega(m)}\right).$$
Indeed, the count of units up to $fm$ equals $f\varphi(m)$ plus an error at most $2^{\omega(m)}$ by inclusion-exclusion. Independently of the multiplier, a translation puts the first image in the window with probability $f$. Hence
$$\mathbb E[X(X-1)]\le 2Qf\left(Qf+\sum_{m\ge1/f}\frac{2^{\omega(m)}}{\varphi(m)}\right),$$
where the sum ranges over products of complete prime-power components of $N$.

## A uniform bound on the tail

Multiplying each summand by $(fm)^{1/2}$ bounds the tail by
$$\sqrt f\prod_{p^e\parallel N}\left(1+\frac{2\sqrt{p^e}}{\varphi(p^e)}\right)\le 9\sqrt{Qf}.$$
For the last inequality, order the odd primes so that $p_i\ge2i+1$. The ratio of the product to $\sqrt Q$ is bounded by the product of $(1+\sqrt{2i+1}/i)/\sqrt2$. Each factor decreases in $i$ and is below one from $i=13$ onward. The product of the first twelve is less than 9, as can be certified by rational upper bounds on their square roots.

Take $C=8192$. If $N/(CQ)\ge1$, choose the largest dyadic $w\le N/(CQ)$, so $1/(2C)<Qf\le1/C$. Then
$$\frac{\mathbb E[X(X-1)]}{\mathbb E X}
\le4\left(Qf+9\sqrt{Qf}\right)<\frac12.$$
Thus $\Pr(X=1)\ge\mathbb E X-\mathbb E[X(X-1)]>1/(8C)$. If $N/(CQ)<1$, width one gives singleton probability $R/N>1/(2C)$ directly. Both cases use a width already in the public dyadic list. Taking $8C$ independent trials at each width gives success probability at least $1-e^{-1}$ even without knowing $Q$. Amplification requires only a further $O(\log(1/\varepsilon))$ factor in the query count. The scalar and the disk cutoff are unchanged. Computing the individual queries in polynomial bit time remains open.

The accompanying check certifies the numerical constant and verifies the ordered-pair identity by exact enumeration for 27 widths over five moduli, including repeated factors. Run `python check_affine_pair_bound.py`.
