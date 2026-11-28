# Uniform intermediate lattices in even and odd dimension

Given $N\ge2$, consider **all** full-rank sublattices $L$ satisfying $N\mathbb Z^d\subseteq L\subseteq\mathbb Z^d$, including both endpoints. Each lattice gets one vote across all possible indices. Its image $H=L/N\mathbb Z^d$ is a subgroup of $(\mathbb Z/N)^d$. For even $d$, we can sample this distribution exactly in expected time polynomial in $d$ and $\log N$. For odd $d$, an exact uniform sampler on arbitrary inputs is equivalent under randomized polynomial-time reductions to factoring. Returning an integer Hermite basis makes the lattice index and the inclusion conditions publicly checkable. [Sampler](sublattice_sampler.py) · [Checks](check_sublattice_sampler.py)

Earlier repository work counted sublattices of a prescribed index and computed individual Hecke neighbors by specifying a projective line. Those restrict the index or the local type. Here even prime powers admit every subgroup, including subgroups that have no direct complement. The local counting formula dates to Birkhoff. Aivazidis gives sharp bounds for the total counts at even and odd rank, including the role of the central terms. We use this established counting theory to construct a sampler when the prime powers of the input have not been supplied. [Earlier fixed-index count](../../notes/geometry_results.md) · [Birkhoff, 1935](https://doi.org/10.1112/plms/s2-38.1.385) · [Aivazidis, Proposition 4.2](https://arxiv.org/html/2606.03320v1#S4)

## The even-dimensional construction

Fix one prime-power component $p^e$ for the analysis. The conjugate partition of the subgroup type consists of integers $d\ge r_1\ge\cdots\ge r_e\ge0$. Here $r_j$ counts cyclic factors of the subgroup of order at least $p^j$. Put $r_0=d,r_{e+1}=0$ and $C_t=\prod_{i=1}^t(1-p^{-i})$. Birkhoff's count for this type has the equivalent forms
$$
W(\mathbf r)=\prod_{j=1}^e p^{r_{j+1}(d-r_j)}
{d-r_{j+1}\brack r_j-r_{j+1}}_p
=p^{\sum_{j=1}^e r_j(d-r_j)}
\frac{C_d}{\prod_{j=0}^{e} C_{r_j-r_{j+1}}}.
$$
The second expression follows by cancelling adjacent factors in the Gaussian binomials. The first expression supplies an integer dynamic program with $O(ed^2)$ transitions. It samples a type in proportion to its count, and can exclude one specified type directly by subtracting its completion weight. [Local subgroup counts, equation 11](https://doi.org/10.1016/j.disc.2013.03.021)

Let $d=2k$. The central type $(k,\ldots,k)$ is a free rank-$k$ direct summand. There is a bound on the probability $\delta_p$ of all other types that is uniform in both dimension and exponent. To see it, move any noncentral coordinate one step toward $k$, choosing the end of a constant run that preserves the nonincreasing order. This raises $\sum r_j(d-r_j)$ by at least one. If the adjacent gaps that shrink and grow were $g\ge1$ and $h\ge0$, the quotient of the remaining factors is $(1-p^{-g})/(1-p^{-(h+1)})\ge1-1/p$. After $T=\sum_j|r_j-k|$ such moves the type has reached the center. Its original weight was at most $(p-1)^{-T}$ times the central weight.

At a given $T\ge1$, the positive deviations at the left and the negative deviations at the right determine the entire type. Removing zeros gives an ordered list of signed positive integers of sum of absolute values $T$. There are at most $2\cdot3^{T-1}\le4^T$ such lists. Summing the weight ratios and normalizing gives
$$\delta_p\le4/(p-1)\qquad(p>5).$$

The sampler treats $2,3,5,7$ by trial division and ordinary local draws. On the rest of $N$, take ten independent uniform residues $x$. For each, trial-factor $g=\gcd(x,N)$. Every remaining prime divisor is discovered independently with probability $h_p=1-(1-1/p)^{10}\ge10/(p+10)\ge4/(p-1)$, since $p\ge11$. For a discovered prime, extract $p^e$ from $N$, compute $\delta_p$ exactly, and draw a noncentral type with probability $\delta_p/h_p$. Use the central type at every other component. This gives independent local types with exactly the target probabilities.

The cost of trial-factoring these random gcds is included. A trial divisor used before the remaining largest prime is returned is at most the largest prime factor of the original $g$. That largest prime is bounded by $\sum_{p\mid N}p\,\mathbf1_{p\mid g}$, whose expectation is $\omega(N)\le\log_2 N$. Removing powers adds $O(\log N)$ divisions. In particular the algorithm does occasionally discover a large prime, at its correct probability.

The central blocks at every unexceptional component are the same $k$ coordinate directions. At the known exceptional components, replace the diagonal generators by the powers of $p$ prescribed by the sampled types. Smith normal form over $\mathbb Z/p^e$ shows that the local general linear group acts transitively on subgroups of a given type. One uniform invertible matrix modulo $N$ therefore produces independent uniform subgroups in the selected local orbits. The column Hermite form of their lifted generators together with $NI$ is the desired integer lattice. Drawing the basis by determinant rejection and computing the integer normal form both take expected polynomial bit time. [Uniform basis routine](idempotent_sampler.py)

## What its index reveals

For an odd-dimensional space over a prime field, orthogonal complementation pairs subspaces of rank $r$ with those of rank $d-r$. Their parities are opposite, and each local parity is fair. On a squarefree semiprime the two local choices are independent. With probability at least $1/2$ their dimensions differ. The returned index $D=[\mathbb Z^d:L]$ then has different valuations at the two primes and immediately gives a proper divisor.

The same index decoder works over arbitrary nonsquare $N$ with at least two prime divisors. Write $e_p=v_p(N)$, $T_p=v_p(D)$. The finite group has a perfect pairing, so taking the annihilator of a subgroup interchanges $T_p$ and $de_p-T_p$. Some $e_p$ is odd. For that prime and odd $d$, the two values are distinct, hence each atom of $T_p$ has probability at most $1/2$. Fixing the local outcome at any other prime fixes at most one integer $j$ for which $D=N^j$ could still hold. The probability of an index that is a power of $N$ is consequently at most $1/2$.

Divide $D$ by $N$ as long as possible. If the remaining integer is greater than one, its gcd with $N$ is a proper divisor, since the original $D$ divides $N^d$. Exact squares and prime powers are recognized directly. Repeating this operation and recursing on the divisors yields full factorization given an odd-dimensional sampler. In the reverse direction, the known factorization and Birkhoff recursion sample the entire law for any dimension.

For even dimension on $N=pq$, the index is $N^{d/2}$ whenever both local types are central. A factor from the **index alone** thus occurs with probability $O(1/p+1/q)$, uniformly in even $d$. For example, at $d=2$ the local fiber has $p+1$ lines and two endpoint subspaces. At $N=19\cdot29$, exactly 102 of its 704 intermediate lattices yield a factor from the index. Increasing the even dimension cannot remove the asymptotic bound on balanced inputs. [Index readout and reproduction](check_sublattice_sampler.py)

## Trying to pass from even to odd

For a free middle-dimensional summand $V$ in even dimension, both intersection with a coordinate hyperplane and projection onto it are public lattice operations. In the generic case they give nested odd-dimensional submodules $A\subset B$ with $B/A\simeq\mathbb Z/N$. The interval between them has exactly one intermediate subgroup per divisor of $N$. Choosing a strict intermediate already supplies a proper divisor by its relative index. This explains why the two easy endpoints do not implement the uniform odd-dimensional draw. [Proof and explicit slicing checks](sublattice-slicing.md)
