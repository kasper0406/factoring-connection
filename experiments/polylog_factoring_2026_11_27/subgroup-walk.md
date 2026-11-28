# An add-and-intersect walk on all subgroups

Let $G=(\mathbb Z/N)^d$. A public Markov step can preserve the uniform law on every subgroup of $G$, including subgroups at repeated prime powers. For $N=pq$ with distinct primes and $d=2k+1$, start at any subgroup with the same rank modulo $p$ and $q$. The expected time to the first **accepted** subgroup with unequal ranks is

$$\Theta\!\left(\frac1{p^{-(k+1)}+q^{-(k+1)}}\right).$$

The constants are absolute. In dimension three and for balanced primes this takes order $N$ steps. The subgroup index yields a verified factor at the first such state. [Public walk](subgroup_walk.py) · [Index decoder](sublattice-sampling.md)

## Transitions and stationary measure

Choose between two proposals by a fair coin. From $X$, add a uniform vector $v\in G$ to form $Y=X+\langle v\rangle$, or draw a uniform character $\chi$ and set $Y=X\cap\ker\chi$. A uniform row in $(\mathbb Z/N)^d$ represents the character by its dot product. For a genuine upward edge $X\subset Y$, the quotient $Y/X$ is cyclic, say of order $h$. Conditional on choosing the upward direction, its proposal probability is $|X|\varphi(h)/|G|$, since each generator of $Y/X$ has $|X|$ lifts. The reverse downward proposal has probability $\varphi(h)/|Y|$. There are $\varphi(h)$ faithful characters of $Y/X$, and each character of $Y$ has $|G|/|Y|$ extensions to $G$ by finite abelian duality.

Accept an upward proposal with probability $\min(1,|G|/(|X||Y|))$ and a downward proposal with the reciprocal ratio clipped at one. The totient cancels. Each oriented edge then has the same transition probability as its reverse, giving detailed balance for the uniform law. Every subgroup can be reached from zero by adding generators, all such edges have positive reverse probability, and drawing the zero vector gives a holding probability. Thus the stationary law is unique. [Hastings, 1970](https://doi.org/10.1093/biomet/57.1.97)

For a column Hermite normal form $H$ of the inverse-image lattice, $|X|=N^d/\det H$. Adjoining the vector as a column gives the upward HNF. For a character row $a$, take the column HNF of the block matrix

$$\begin{pmatrix} H&0\\aH&N\end{pmatrix}.$$

Its first $d$ columns span the section at last coordinate zero. Their upper $d$ entries give the intersection lattice. Both updates and the exact rational acceptance test take polynomial bit time.

## First unequal rank

Write $\epsilon=p^{-(k+1)}+q^{-(k+1)}$. At an aligned state of rank $(r,r)$, the probability, conditional on an upward proposal, of staying put at $p$ and increasing rank at $q$ is $p^{r-d}(1-q^{r-d})$. Its acceptance probability is

$$\min(1,p^{d-2r}q^{d-2r-1}).$$

For $r\le k$ the proposal probability is already at most $p^{-(k+1)}$. For $r\ge k+1$ the accepted probability is bounded by $p^{-r}q^{d-2r-1}\le p^{-(k+1)}$. Interchanging the primes accounts for the other mixed move. Taking annihilators sends rank $r$ to $d-r$ and exchanges the two directions. Thus the conditional chance of leaving the aligned states in each complete step is at most $\epsilon$. For the first exit time $T$, this gives $\mathbb E T\ge1/\epsilon$.

For the reverse bound, track the distance from $r$ to $\{k,k+1\}$ until a split occurs, with distance zero on absorption. Below the central pair, proposing up and increasing both ranks has probability at least $1/8$ and is accepted. An outward move decreasing both ranks has probability at most $N^{-4}/2$, by the downward acceptance ratio $N^{2r-1-d}$. Annihilator duality gives the same bounds above the pair. The expected decrease in distance outside the center is at least $1/9$. Summing this drift up to the first central visit or absorption gives an expected time of at most nine times the initial distance.

At rank $k$ an upward proposal increases just one local rank with probability $a(1-b)+b(1-a)$, where $a=p^{-(k+1)}$ and $b=q^{-(k+1)}$. Each such proposal is accepted. Accounting for the direction coin gives at least $\epsilon/4$ per visit. At rank $k+1$ use the downward proposal. A nonsplitting departure from the central pair has distance one, so one central attempt and a return or absorption cost at most ten steps in expectation. The expected number of attempts is at most $4/\epsilon$. Altogether,

$$\frac1\epsilon\le\mathbb E T\le 9d+\frac{40}\epsilon.$$

Since $d\epsilon$ is absolutely bounded for distinct primes, this proves the stated order. Under the uniform law the local ranks are independent. At either prime the weights of ranks $r$ and $d-r$ agree, and odd $d$ pairs every rank with a different one. Each atom is therefore at most $1/2$. The stationary probability of unequal ranks is at least $1/2$. From an aligned start, total variation distance at step $t$ is at least $1/2-t\epsilon$.

**Reading proposals before acceptance.** At equal rank $r$, either raw direction produces an unequal pair with probability at most $1/p+1/q$. From rank one, a raw intersection attains this order. At $N=15$, $d=3$, the raw step from $(1,1)$ to $(1,0)$ has probability $2/15$, including the direction coin. The walk accepts it with probability $1/75$. An implementation inspecting every raw index can therefore find a factor before its first accepted rank split.

Fulman studies the related fixed-rank down-up walk over one finite field. The earlier scalar-annihilator experiment studied a different Gibbs update. Here the CRT rank calculation uses the common proposal direction and the cardinality-dependent acceptance ratio. [Fulman, §8.2](https://arxiv.org/pdf/0712.1375#page=31) · [Earlier conditional sampler](four-square-incidence.md)

The diagnostic enumerated 6,582 raw proposals on 81 subgroup states, including $\mathbb Z/4$ and $\mathbb Z/9$ components, and verified detailed balance exactly. It solved 16 first-hit systems over the rationals. For $N=101\cdot103$ and $d=3$ the expected first split from rank $(1,1)$ is approximately $10{,}401.010$ steps. [Reproducible check](check_subgroup_walk.py)
