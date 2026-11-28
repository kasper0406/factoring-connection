# Uniform matrix projectors over an unfactored integer

Let $R=\mathbb Z/N$. A projector is a matrix $E$ satisfying $E^2=E$. For even $d>0$ there is an exact uniform sampler for the entire set of projectors in $M_d(R)$, with expected bit cost polynomial in $d$ and $\log N$. For odd $d$, a uniform draw splits any odd semiprime with probability exactly $1/2$ using just $\det(I-2E)$. In particular increasing the even dimension by one changes the computational problem. The earlier scalar-idempotent experiment counted solutions, and the earlier four-square incidence experiment had a different local counting law. [Scalar count](idempotent_count.py) · [Four-square incidence](four-square-incidence.md) · [Binary matrix sampler](idempotent_sampler.py)

## More than two labeled summands

An ordered decomposition consists of $s$ matrices $E_i$ with $E_i^2=E_i$, $E_iE_j=0$ for $i\ne j$, and $\sum E_i=I$. Zero summands are allowed. A single projector specifies the pair $(E,I-E)$. The exact sampling classification for positive $s,d$ is

| Dimensions | Uniform sampling from $N$ |
|---|---|
| $s\mid d$ | An exact sampler takes expected time polynomial in $s,d,\log N$. |
| $s\nmid d$ | A sampler for arbitrary $N$ is equivalent under randomized polynomial-time reductions to factoring. On a modulus with two distinct prime divisors, a sample gives a verified proper divisor with probability at least $1-1/s$. |
| Any $s,d$, restricted to free summands | Exact uniform sampling under this restriction takes expected polynomial time. The restriction means each summand has the same rank at every prime dividing $N$. |

The field orbit counts and the classification of projectors over finite local rings are classical. The algorithms below use those counts while taking the composite as their only number-theoretic input. [Ellerman, direct-sum counts over finite fields](https://arxiv.org/abs/1603.07619) · [Dolžan, idempotents in finite rings](https://doi.org/10.1080/00927872.2021.1929277)

### The local law

Over $\mathbb Z/p^e$, every projective summand is free. A decomposition has a profile $\mathbf r=(r_1,\ldots,r_s)$ of nonnegative ranks summing to $d$. All decompositions of that profile are conjugates of the standard block projectors. Their stabilizer is the product of the general linear groups on the blocks. Write
$$G_j=p^{e j^2}C_j(p),\qquad C_j(p)=\prod_{h=1}^j(1-p^{-h}),\qquad G_0=1.$$
The number with profile $\mathbf r$ is $w(\mathbf r)=G_d/\prod_iG_{r_i}$. A uniform decomposition over $R$ therefore has independent prime-power profiles with these weights. Given the factorization, the integer coefficient recursion for $(\sum_{r=0}^d (G_d/G_r)x^r)^s$ samples the local profiles in polynomial bit time. One uniform element of $\mathrm{GL}_d(R)$ then conjugates the combined block projectors. Integer matrix rejection uses the determinant gcd. Its acceptance probability is at least $(\varphi(N)/N)\prod_{j\ge2}\zeta(j)^{-1}$, since the local probabilities multiply.

### How to sample when $s\mid d$

The central profile is $(m,\ldots,m)$ with $m=d/s$. Write $a_i=r_i-m$. Since $\log C_j(p)$ is a convex sequence, the weight ratio to the central profile is at most $p^{-e\sum a_i^2}$. Set $L=\sum_{a_i>0}a_i$. For a nonnegative integer $t$, $p^{e(t^2-t)}\ge t!$. Summing over positive and negative parts separately and applying the multinomial formula bounds the total relative weight at distance $L$ by $(s^2/p^{2e})^L/(L!)^2$. The central probability is consequently at least
$$\exp(-s^2/p^{2e})\ \ge\ (1-p^{-2})^{s^2}.$$
This supplies a way to find **every exceptional local profile required for one exact draw**, while leaving the other prime divisors unidentified.

1. Draw $s^2$ independent pairs $(x,y)$ uniformly modulo $N$. Trial-factor each $g=\gcd(x,y,N)$. Every prime divisor $p$ is discovered independently with probability $h_p=1-(1-p^{-2})^{s^2}$, by CRT.
2. For a discovered prime, extract its full power in $N$. Compute its exact probability $\delta_p$ of a noncentral profile from the local weights. With probability $\delta_p/h_p$, use the coefficient recursion to draw a noncentral profile, excluding the central tuple directly. All remaining prime powers get the central profile.
3. Install the exceptional blocks with CRT and conjugate the projectors by one uniform invertible matrix. Uniformity follows from independence across primes and transitivity within each orbit.

The trial factorization has a polynomial **expected** cost. Indeed $\mathbb E\sqrt g\le\sum_{b\mid N}\sqrt b\,\Pr[b\mid g]=\sum_{b\mid N}b^{-3/2}\le\zeta(3/2)$. This accounts for rare large gcds as well. The binary implementation also gives an exact sampler for a sparse random set of integers with independent inclusion probabilities $1/j^2$ using the telescoping product of $1-1/j^2$. Two passes suffice in the binary case, as shown below. Kalai's descending-chain construction is a related way of generating sparse independent prime events. [General implementation](projector_decompositions.py) · [Binary implementation](idempotent_sampler.py) · [Kalai's construction](https://www.cs.cmu.edu/afs/cs/user/akalai/genfactor/web/factor/factor.html)

For $s=2$, changing the local rank of $E$ by $j$ from $d/2$ has relative weight at most $p^{-2ej^2}$. The sum of all noncentral weights relative to the central one is at most $2/(p^{2e}-1)$. Thus $\delta_p\le2/(p^{2e}+1)\le 1-(1-p^{-2})^2$. For balanced semiprimes, even granting a proper divisor whenever either rank is exceptional gives only $O(1/p^2+1/q^2)$ success. On $N=35$ a public sample with seed 8 does give $\det(I-2E)=6\pmod{35}$, and $\gcd(6-1,35)=5$. [Exact binary checks](check_idempotent_sampler.py)

### How a draw identifies a prime component when $s\nmid d$

Every profile then has at least two distinct entries. Permuting summand labels preserves weight and gives at least $s$ profiles in its orbit: the positions of any rank of multiplicity $b$ yield at least $\binom{s}{b}\ge s$ choices for $1\le b<s$. No local profile has probability exceeding $1/s$. For any two distinct prime components their independent profiles disagree with probability at least $1-1/s$.

A rank disagreement of one $E_i$ gives a factor using its characteristic polynomial. Over a prime-power component where its rank is $r$, this polynomial is $T^{d-r}(T-1)^r$. Scan coefficients from constant term up. At the first degree whose coefficient is nonzero modulo $N$, that coefficient is $\pm1$ at every component of largest rank and zero at each component of smaller rank. A gcd separates the two groups, including their full prime powers. This works even for ranks larger than a residue characteristic. Repeating on the resulting composite factors and handling prime powers gives a randomized complete factorization if uniform sampling is available on arbitrary inputs.

For binary projectors with odd $d$ and odd $N$, $E\mapsto I-E$ pairs ranks $r,d-r$ of opposite parity at every prime power. The signs of $\det(I-2E)=(-1)^r$ are independent fair bits at the different prime components. On $N=pq$ the gcd with $\det(I-2E)-1$ splits exactly half the time. [Rank decoder and exact checks](check_projector_decompositions.py)

### Sampling exactly up to the first factor-revealing distinction

For any $s,d$, a decomposition over $R$ has all summands free precisely when its local profiles agree. There is a public sampler for the exact conditional uniform law on this whole event. This includes rank profiles arbitrarily far from the most likely ones. [Sampler](free_projector_sampler.py)

Write $a=\lfloor d/s\rfloor$ and $b=d\bmod s$. A balanced profile has $b$ entries $a+1$ and $s-b$ entries $a$. There are $K=\binom{s}{b}$ such profiles, of equal weight. For a common profile the global orbit count is $|\mathrm{GL}_d(R)|/\prod_i|\mathrm{GL}_{r_i}(R)|$. Moving one rank unit from a larger entry to a smaller entry differing by at least two raises this count by at least $N^2$, by the local order formula. Every unbalanced profile can be reached from a balanced one by reversing these transfers. There are at most $s(s-1)$ choices per step. A geometric series bounds the probability of an unbalanced profile under the conditional law by $h=\min(1,s(s-1)/N^2)$. For $h<1$ the ratio of the total unbalanced weight to the balanced weight is at most $h/(1-h)$, giving this bound after normalization.

Start with a uniformly permuted balanced profile. With probability $h$, trial-factor $N$, compute the true unbalanced probability $\delta$ from the global orbit counts, and replace the profile with probability $\delta/h$. The replacement is an exact draw conditioned on being unbalanced. The coefficient recursion subtracts the balanced suffix mass at each step. Conjugate by one uniform global basis. The expected cost of trial division is $O(h\sqrt N)$ divisions and is polynomial in $s,\log N$. The sampler requires no supplied factors and can explicitly factor the input on its rare correction branch.

Let $U$ denote the unrestricted uniform law and $F$ the common-profile event. The conditional sampler has total variation distance exactly $1-U(F)$ from $U$. This is optimal among distributions supported on $F$. More generally any sampler within distance $\epsilon$ of $U$ reveals a rank mismatch with probability at least $1-U(F)-\epsilon$. At a balanced semiprime and fixed $s$, $U(F)=1/K+O_s(p^{-2}+q^{-2})$. Indeed at each prime the total mass outside the balanced profiles is $O_s(p^{-2})$ by the same transfer argument. Two balanced profiles agree with conditional probability $1/K$. Two exceptional profiles can agree only on an event of probability at most the product of the two exceptional probabilities.

At $N=551=19\cdot29$, the exact values of $U(F)$ are about $0.992489$ for $(d,s)=(2,2)$, $0.499996$ for $(3,2)$, and $0.332072$ for $(4,3)$. The tests check 408 exceptional profiles at both endpoints of their exact sampling intervals and 192 whole-modulus draws. [Reproduce the exact fractions](check_free_projectors.py)

An explicit attempt to use the easy even sampler to obtain the hard odd sampler is to remove one coordinate. Two Schur complements of the same even projector give two adjacent odd-dimensional ranks. Interpolating between the two at different hidden primes requires a scalar idempotent, from which a gcd immediately splits $N$. This identifies the missing operation for that construction. [Compression and checked examples](projector-compression.md)
