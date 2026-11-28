# Slicing and projecting an even-dimensional sublattice

Let $R=\mathbb Z/N$, let $V\subset R^{2k}$ be a uniform free rank-$k$ direct summand, and put $H=R^{2k-1}\times0$, $K=0\times R$. Intersecting with $H$ gives $A=V\cap H$ and dropping the last coordinate gives $B=\operatorname{pr}_H(V)$. Identify $H$ with $R^{2k-1}$, so $A\subseteq B$. This experiment fixes the free rank before sampling. [Public construction](sublattice_sampler.py)

For a prime $p\mid N$, reduction gives a uniform $k$-subspace $\bar V$ of $\mathbb F_p^{2k}$. The events $\bar V\subseteq\bar H$ and $\bar K\subseteq\bar V$ are disjoint and both have probability $1/(p^k+1)$. For the second event, symmetry gives $(p^k-1)/(p^{2k}-1)$ by counting the nonzero vectors in $\bar V$. Annihilator duality gives the first. Local choices are independent, because a uniform change of basis over $R$ has independent uniform prime-power components. The probability of avoiding both events at every prime is therefore

$$
\prod_{p\mid N}\left(1-\frac{2}{p^k+1}\right).
$$

On this good event the last coordinate map $V\to R$ is onto and $V\to B$ is an isomorphism with split image. Indeed, modulo each prime the first map is nonzero and the second has rank $k$. A unit minor at each local component lifts the corresponding splitting over every prime power. It follows that $A$ has free rank $k-1$, $B$ has free rank $k$, and $B/A\simeq R$. Conditioned on this event the flag $A\subset B$ is uniform over such split flags. For each flag, $V$ is the graph of a surjection $B\to R$ with kernel $A$, and the $\varphi(N)$ isomorphisms $B/A\to R$ give all preimages.

## Reading the interval from an integer basis

For the inverse-image lattice of $V$, take its column Hermite normal form and drop the last row. The first $2k-1$ remaining columns span the inverse-image lattice $\widetilde A$ in $H$, and adjoining the last truncated column $b$ gives $\widetilde B$. For the first assertion one can choose a lift with last coordinate zero since the lattice contains $N\mathbb Z^{2k}$. Upper triangularity then forces the last basis coefficient to vanish. Thus $B/A$ is cyclic on every input, and its order $D=[\widetilde B:\widetilde A]$ is the ratio of the two basis determinants and divides $N$.

The intermediate groups have inverse-image lattices $\widetilde A+\mathbb Z\,t b$, one for each $t\mid D$, with index $t$ in $\widetilde B$. On a good sample $D=N$, so producing a strict intermediate group is equivalent to producing a proper divisor of $N$. Either endpoint gives a common rank across the prime components. On a squarefree semiprime $D$ equals the product of the primes whose reductions are good. A mixed good/bad sample already yields a proper factor through this index. If $g_p=1-2/(p^k+1)$, the split probability from $D$ alone is $g_p(1-g_q)+(1-g_p)g_q$. Increasing $k$ reduces this probability for large balanced primes.

The exact diagnostic found 104 good subspaces among the 130 planes in $\mathbb F_3^4$, and 1,085 among 1,395 three-spaces in $\mathbb F_2^6$. Over $\mathbb Z/9$ the twelve free lines gave index 1 twice, 3 four times, and 9 six times. At $N=15$ a good rank-two example in dimension four gave the four intermediate indices $1,3,5,15$. Run `python experiments/polylog_factoring_2026_11_27/check_sublattice_slicing.py` from the repository root. [Exact checks](check_sublattice_slicing.py)
