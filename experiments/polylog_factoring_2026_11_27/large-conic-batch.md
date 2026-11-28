# A native batch for a large first coordinate

For a balanced odd semiprime $N=pq$, $p<q<2p$, fix the public first vector $(x,1,0)$. Its noncoordinate partners on a mixed Jacobi input exist exactly when $r$ or $2r$ has a representation $a^2+(x^2+1)b^2$, $b\ge1$, for the splitting prime $r\mid N$. Given a single $x$ of size comparable to $N^{1/4}$ where such a representation exists, one can **factor deterministically in $N^{1/8+o(1)}$ bit operations**. The representation and the factors are both hidden. The method works with ordinary integer polynomials and gcds. A density bound below limits how often *any* choice of large $x$ can work, even when it depends on the entire input. [Fixed-vector representation criterion](short-first-vector.md) · [Rotation proof](ternary-frame-rotations.md) · [Implementation](large_conic_batch.py)

## Public computation

Put $H=\lfloor\sqrt{8N}\rfloor$, $A=x^2+1$, and $B=\lfloor\sqrt{H/A}\rfloor$. For each $1\le b\le B$, search $0\le a\le\lfloor\sqrt{H-Ab^2}\rfloor$ for a nontrivial gcd of $a^2+Ab^2$ and $N$. Both $2p$ and $2q$ are in range. Also $H<4p$ and $H<3q$, so every actual hit has $a^2+Ab^2=er$ for $e\in\{1,2,3\}$. The only extra possibility beyond the original representation is $3p$. These statements use strict balance and positive terms. No Jacobi test is needed to run the batch.

For a row of $K$ consecutive values of $a$ starting at zero use $s=\lfloor\sqrt K\rfloor$. Form
$$
P(T)=\prod_{0\le j<s}\big((T+j)^2+Ab^2\big)\pmod N
$$
with a balanced product tree and evaluate at $0,s,2s,\ldots$ for the complete blocks with a monic remainder tree. Inspect block gcds in order. In the first nonunit block, at most $s$ individual gcds recover the witness. This works even if the block gcd is $N$: every term lies in $[1,H]$ with $H<N$, so every nonunit term gives a proper divisor. Check the final incomplete block directly. Hence the routine returns the same first witness as an exhaustive scan in $(b,a)$ order.

There are $O(N^{1/4})$ possible $a$ per row. Fast polynomial multiplication and monic division give time $B N^{1/8+o(1)}$ when $B\ge1$ and $N^{1/8+o(1)}$ bits of working memory. When $B=0$ the test ends after computing the limits. The implementation uses the existing exact polynomial backend with a dyadic FFT fallback. Block evaluations are polynomial operations and cost more than scalar gcds. For comparison, scanning the original parameter $t$ uses $\lfloor2N/(x^2+1)\rfloor$ square tests. [Polynomial backend](sqrt_prime_marks.py) · [Earlier interval-product baseline](../../src/trackC/baselines.py) · [Bostan, Gaudry and Schost on fast polynomial products](https://mathexp.eu/bostan/publications/BoGaSc07.pdf)

## How many large settings can work?

Let $p,q\in[R,2R)$ and $x\in[X,2X)$ with $X\ge1$. Any hit has $1\le b\le2\sqrt R/X$. Fix a prime $r$ in this range. For each $b,e$, successful $x$ are a subset of the solutions
$$
a^2+c^2=er-b^2,\qquad c=bx>0 .
$$
The right side is positive for a solution. Jacobi's two-square formula bounds its number of ordered signed representations by four times its number of divisors. Thus each $r$ has at most $(2\sqrt R/X)R^{o(1)}$ successful integer settings. At $X\asymp\sqrt R$ a uniformly sampled integer $x$ hits on at most $R^{-1/2+o(1)}$ of trials on **each input**, before any conditioning on the Jacobi sign. For a fixed $x$, counting $(a,b)$ gives at most $O(R/X)$ possible primes. A menu of $J$ fixed settings in the band, visited in any order, therefore reaches at most $O(J\log R/X)$ of balanced pairs. This last bound concerns menus listed independently of the input. [Jacobi's sum-of-two-squares theorem](https://mathworld.wolfram.com/SumofSquaresFunction.html)

A separate sieve bounds **all input-dependent settings** simultaneously. If $1\le b\le L$, the fraction of prime pairs admitting even one hit with such a $b$ is
$$
\ll \min\{1,L/\sqrt{\log R}\}.
$$
In particular all $x\ge c\sqrt R$ together fail on a density-one subset for each fixed $c>0$. More generally this holds for $x\ge2\sqrt R/L$ when $L=o(\sqrt{\log R})$. The estimate counts the extra $e=3$ hits and therefore applies to the actual gcd algorithm. It also allows for both prime sides. It bounds the coverage of this particular scan.

To prove it, fix $e\in\{1,2,3\}$ and $b\ge1$. If $er-b^2=a^2+c^2>0$, its valuation at a prime $\ell\equiv3\bmod4$ is even. Sieve all integers $r\in[R,2R)$, excluding multiples of $\ell$ at every $\ell\le z$. For $\ell\equiv3\bmod4$ not dividing $2eb$, also exclude $v_\ell(er-b^2)=1$. Among the $\ell^2$ classes modulo $\ell^2$, the first restriction excludes $\ell$ classes. At these good inert primes the second excludes $\ell-1$ further unit classes. Thus the surviving density at a prime is $(1-1/\ell)^2$ at good inert primes and $1-1/\ell$ elsewhere.

For squarefree $m$, CRT gives the intersection of excluded sets main density $g(m)=\prod_{\ell\mid m}g(\ell)$ with interval-count error at most $m2^{\omega(m)}$, where $g(\ell)$ is the excluded density. Use the fundamental lemma of sieve theory in its upper-bound form with level $D=R^{1/3}$ and $z=D^{1/s}$ for a fixed sufficiently large $s$. Its dimension hypothesis holds uniformly with dimension two since $1-g(\ell)\ge(1-1/\ell)^2$. The summed errors are $O(D^2\log D)$. All primes $r$ under consideration exceed $z$. The main term is bounded by
$$
R\prod_{\ell\le z}(1-g(\ell))
\ll {R\over(\log R)^{3/2}}\,{eb\over\varphi(eb)} .
$$
Here the first product over all primes contributes $1/\log z$ and the second over primes $3\bmod4$ contributes $1/\sqrt{\log z}$ by Mertens' theorem in fixed progressions. Removing primes dividing $2eb$ costs at most a constant times $eb/\varphi(eb)$. The error is smaller than the displayed bound. Finally $\sum_{b\le L} b/\varphi(b)\le L\prod_\ell(1+1/(\ell(\ell-1)))=O(L)$. Sum over the three $e$, normalize by the number of primes in the interval and allow either side of the pair. [Tao, Notes 4, Corollary 19, upper-bound sieve for arbitrary forbidden sets](https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/) · [Tao and Teräväinen, Lemma 5.1](https://www.utupub.fi/server/api/core/bitstreams/8b6b144a-c064-4a7b-8b57-ed10ecf691d9/content#page=15)

## Exact checks

A planted 54-bit example is $N=10412720415486989$, $x=8192$. We selected prime $p=4026^2+(8192^2+1)=83317541$ and a balanced prime $q=124976329$ with negative local sign. Receiving only $N,x$, the batch found $p$ using 122 block evaluations of a degree-244 polynomial and one individual gcd, following at most 122 block gcds. The original $t$ scan has 310,323,246 positive positions. At $N=91,x=4$, the terms 21 and 26 instead saturate a single block, and its first individual gcd returns 7 from the additional multiplier 3. [Run the checks](check_large_conic_batch.py)

For each row below we used **all** primes in $[R,2R)$, all distinct unordered pairs, and all integers $\sqrt R\le x<2\sqrt R$. The genuine columns count solutions with $e=1,2$, without imposing a mixed Jacobi sign. The operational columns count any proper gcd the batch will return, including $e=3$. A pair is covered if at least one of its two prime sides has a successful setting. The checker runs the actual batch and independent exhaustive scan on every one of the 624 small settings at $R=64$. At larger $R$ it exhaustively counts the equivalent equations $er=a^2+(x^2+1)b^2$ using the known cohort primes to grade coverage.

| $R$ | Primes | Primes with genuine setting | Genuine $(r,x)$ | Pairs | Pairs with genuine setting | Pairs with any operational hit | Operational $(N,x)$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 13 | 8 | 14 | 78 | 68 | 74 | 198 |
| 256 | 43 | 25 | 49 | 903 | 750 | 815 | 2,298 |
| 1,024 | 137 | 83 | 154 | 9,316 | 7,885 | 8,268 | 24,379 |
| 4,096 | 464 | 237 | 465 | 107,416 | 81,765 | 88,235 | 260,327 |

At $R=4096$ the 260,327 operational successes occupy 3.79% of 6,874,624 tested settings, while 82.1% of pairs have at least one. The density-zero estimate is asymptotic and permits substantial coverage at these scales.
