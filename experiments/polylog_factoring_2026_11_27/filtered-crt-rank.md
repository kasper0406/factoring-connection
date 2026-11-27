# Exact CRT rank after imposing the ellipse congruence

For the public near-square branch of a balanced semiprime, the rank of the mixed weighted disk table equals the number of mixed positive points. This gives an exponential rank example even after imposing the congruence. Its source is repeated multiples of one point, and those multiples admit a closed-form sum once the base point is available.

## Exact rank

Let $N=pq$ with distinct odd primes $p<q\le2p$. Take the public $c$ coprime to $N$, $\Delta=4N-c^2$ and $A=\lfloor\sqrt{\Delta/3}\rfloor$ with $3\le A<p$. Use positive points with $w^2+\Delta y^2\le4NA$ and $4N\mid w^2+\Delta y^2$. Remove the public ray $w=cy$. In a matrix over $\mathbb F_p$, index the rows by $(w,y)\bmod p$ and the columns by $(w,y)\bmod q$. Put weight $w\bmod p$ at each remaining point and zero elsewhere.

The cap gives $w<N$ from $4A<N$, and $y^2\le4N/(3A)<p^2$. Therefore $y$ is a unit modulo $N$ and the congruence implies $w$ is a unit as well. The congruence modulo 4 forces $w\equiv cy\pmod2$. If the signs of $w/(cy)$ agree at the two primes, then $w\mp cy$ is divisible by $2N$ for the corresponding sign. Cauchy–Schwarz gives
$$|w\mp cy|^2\le(w^2+\Delta y^2)(1+c^2/\Delta)\le16N^2 A/\Delta<4N^2.$$
Here $\Delta\ge3A^2>4A$. Agreement of signs thus gives only the public ray among positive points.

At a fixed $y$, the two possible mixed residues of $w$ are negatives modulo $N$. They cannot both occur for positive $w<N$ of the same parity, since their sum would have to be the odd integer $N$. Thus each remaining point has a distinct $y$. As $y<p<q$, their row indices are pairwise distinct and their column indices are pairwise distinct. Every weight is nonzero, so the compressed table is a permutation matrix with nonzero weights. Its rank equals the number of mixed positive points. [Fixed-cap setup](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_aggregate/weighted_recurrence_single.md)

## An infinite family with large rank

Write $p=m-d$ and $q=m+d$. For any fixed $d$ and sufficiently large $m$, the public branch has $c=2m-1$ and $\Delta=4m-1-4d^2$. A mixed norm-one point is
$$w_0=2m-4d^2,\qquad y_0=2d.$$
Direct expansion gives $w_0^2+\Delta y_0^2=4N$. Also $w_0-cy_0=2(1-2d)q$ and $w_0+cy_0=2(1+2d)p$, so the signs are mixed for sufficiently large $p,q$. Its positive multiples through $k=\lfloor\sqrt A\rfloor$ all survive.

In fact these are all the mixed points for these inputs. Reflect the ordinate of each point as needed so that any point and $(w_0,y_0)$ lie on the same mixed slope modulo $N$. Their determinant is divisible by $N$, and coordinate parity makes it even. By the norm bounds its square is at most $16N^2 A/\Delta<4N^2$, so the determinant vanishes. If a collinear point scales the norm-one vector by a positive rational $a/b$ in lowest terms, its norm index is $a^2/b^2$, a positive integer, which forces $b=1$. The rank is therefore exactly $\lfloor\sqrt A\rfloor=\Theta(N^{1/8})$. The bounded prime gaps theorem implies that some fixed $d$ occurs infinitely often. [Maynard, Small gaps between primes](https://arxiv.org/abs/1311.4600)

## Checks and consequence

The independent shell scan and compressed-rank check passed on 545 small branches satisfying the hypotheses. It also recovered ranks 3, 6, 10 and 1 at $N=11{,}009$, $1{,}022{,}117$, $100{,}160{,}063$ and $38{,}942{,}974{,}393$, respectively. Run `python check_filtered_rank.py` in this directory.

An explicit sum of rank-one matrices representing the full mixed table can be exponentially large even after filtering. The prime-gap examples themselves are easy for Fermat factoring. They also explain how an aggregate can evade the table bound. The repeated points contribute $w_0\lfloor\sqrt A\rfloor(\lfloor\sqrt A\rfloor+1)/2$ in closed form. The fixed-cutoff question is how to obtain the needed direction contributions from public input without enumerating candidate directions or finding the factors first.
