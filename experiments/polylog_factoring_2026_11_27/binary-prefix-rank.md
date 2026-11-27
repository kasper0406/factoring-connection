# Binary matrix width for the general root prefix

The root-prefix parity in the general-input roadmap needs large width in a binary matrix product that reads the endpoint bits once in order. On products of the first $k$ odd primes, the needed width is $2^{\Omega(n/\log n)}$ for input bit length $n$. The result allows the matrices to depend arbitrarily on the modulus. It addresses the general-input oracle. A semiprime has only two roots in the lower half and does admit a short representation once its mixed root is supplied. [Root-prefix reduction](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_disk_mod8/native_dyadic.md)

## Local spacing by distance between CRT signs

Let $N$ be an odd squarefree integer. A root of $x^2=1\pmod N$ corresponds to a vector of signs at the primes. Give the coordinate for a prime $p$ weight $\log p$. For distinct roots $x,y$, the product of primes where their signs agree divides $x-y$. The weighted Hamming distance is therefore at least $\log N-\log(|x-y|)$.

The same argument applies to distinct roots in one residue class modulo a power of two $m$. In that case the product of agreeing primes divides $(x-y)/m$, since it is odd. View that progression in its integer quotient coordinate. In either setting, an interval of diameter at most $\sqrt N/2$ contains $O(\log N)$ roots. Indeed, for $R$ roots in the interval the pairwise weighted distance is at least $\frac12\log N+\log2$. The sum of distances is at most $\frac14 R^2\log N$ because each coordinate contributes at most $R^2/4$ differing pairs. Combining the two bounds gives $(R-1)\log2\le\frac12\log N$.

## Rank across the middle cut

Take $m$ to be the greatest power of two at most $\sqrt N$, and write $x=mj+i$ with $0\le i<m$. Keep only the complete columns with $0\le x\le (N-1)/2$. A fixed column is an ordinary interval of length at most $m$. A fixed row is a power-of-two progression whose quotient range has length $O(\sqrt N)$. Split each range into a constant number of intervals of diameter at most $\sqrt N/2$. By the spacing bound, the bipartite support graph of the resulting binary root table has maximum degree $d=O(\log N)$. There are $2^{\omega(N)-1}$ roots in the lower half. Omitting the final incomplete column removes at most $O(\log N)$ roots by the same interval bound.

A graph of maximum degree $d$ with $E$ edges has an induced matching of size at least $E/(2d^2)$. Greedily select an edge and delete all edges incident to neighbors of either endpoint. At most $2d^2$ edges are deleted per selection, and each later edge has no cross-edge with a selected edge. The induced rows and columns of the root table form a permutation submatrix. Thus its binary rank is at least $(2^{\omega(N)-1}-O(\log N))/O(\log^2 N)$.

For the prefix parity table, take adjacent differences down the rows in each column. This gives the root table except for a correction supported in the first row, namely the prefix at the last position of the preceding column. Therefore the rank of the root table is at most the prefix rank plus one. A binary matrix product over $\mathbb F_2$ reading the endpoint bits from either end factors the prefix table through its internal width at this cut, regardless of the matrices' dependence on $N$ and the layer. Its width is at least the prefix rank. On odd primorials $n=\Theta(k\log k)$ by the prime number theorem, yielding the stated bound.

The proof leaves adaptive algorithms, different variable orders, nonlinear representations and evaluators for only selected endpoints open.

## Checks

The check script constructs the support by CRT on 11 odd primorials, finds induced matchings and explicitly compares the binary ranks on the five moduli below one million. The 48-bit case had 2,048 matched roots. Run `python check_binary_prefix_rank.py`. Factors are supplied to this diagnostic only.
