# Sampling four squares with scalar annihilators

For an odd positive integer $N$, consider the finite set of tuples $(x,t_1,\ldots,t_k)$ where $x$ is an ordered, signed integer quadruple of squared norm $N^2$, each $t_i$ is a residue modulo $N$, and $t_i x_j=0\pmod N$ for every coordinate $j$. Put the uniform distribution on the **tuples**. A point admitting more choices of the scalars occurs more often.

| Number of scalars | Exact uniform sampling from the public input |
| --- | --- |
| 0, 1, or any $k\ge3$ | Expected polynomial time in $\log(2N)$ and the output length. |
| 2 | Equivalent to randomized polynomial-time factoring. For $N=pq$ with distinct odd primes, one uniform tuple reveals a proper factor with probability $1/2-1/[2(2p+1)(2q+1)]$. |

The earlier quaternion experiment compared lines obtained from **two independent points of norm $N$**. The question here concerns the distribution of **one point of norm $N^2$** when equipped with scalar annihilators. The fixed-norm ideal sampler has published precedents. The deduction below concerns how to sample across different coordinate gcds, including for nonsquarefree inputs. [Earlier experiment, row R7](../../notes/conn_alg.md) · [SQIsign2D–West, §3.3](https://www.math.u-bordeaux.fr/~damienrobert/pro/publications/articles/sqisign2d.pdf#page=14) · [Pollack and Treviño, §4](https://math.colgate.edu/~integers/sjs15/sjs15.pdf)

## Which points carry a factor?

Write $g=\gcd(x_1,\ldots,x_4)>0$ for the content of $x$. Then $g\mid N$. Write $m=N/g$ and $\eta(m)=\prod_{r\mid m}(1+1/r)$ over primes, with $\eta(1)=1$. Jacobi's formula gives $8\sigma(s)$ representations of an odd $s$ by four squares. Removing common square divisors by Möbius inversion gives $8s\eta(s)$ primitive representations. At $r^a$ with $a\ge2$ the local difference is $\sigma(r^a)-\sigma(r^{a-2})=r^a+r^{a-1}$, and the same formula holds for $a=1$. [Jacobi's four-square theorem](https://doi.org/10.1017/S144678870002440X)

There are therefore $8m^2\eta(m)$ points of content $g$. Exactly $g$ scalars annihilate such a point, namely the multiples of $m$ modulo $N$. The exact probability of content $g$ is proportional to

$$g^k\,8m^2\eta(m)\ \propto\ g^{k-2}\eta(N/g),\qquad g\mid N.$$

For $N=pq$ the two decisions $p\mid g$ and $q\mid g$ are independent. At $r\in\{p,q\}$ the inclusion probability is $r^k/(r^2+r+r^k)$. Thus on balanced pairs with $p,q\in[R,2R]$ the nearby cases $k=1$ and $k=3$ each produce a proper content only with probability $O(1/R)$. At $k=1$ most samples have $g=1$. At $k=3$ most have $g=N$. At $k=2$ the two primes each have probability close to one half of dividing $g$.

## Sampling when $k\ne2$

Set $h=|k-2|\ge1$. Take $v$ to be the gcd of $N$ and $h$ independent uniform residues modulo $N$. Set $g=v$ for $k<2$, and $g=N/v$ for $k>2$. The proposal probability for a particular $v\mid N$ is

$$v^{-h}\prod_{r\mid N/v}(1-r^{-h}).$$

Two rejection steps correct this weight. For a divisor $s\mid N$ let $T(s)=N/\gcd(N,s^b)$, where $b$ is the bit length of $N$. This selects the full prime powers of $N$ at primes absent from $s$ and can be computed by modular powering and a gcd.

1. Draw $h$ uniform residues modulo $T(N/v)$. Accept if their joint gcd with that modulus is one. This inserts all missing factors $1-r^{-h}$.
2. Modulo $T(m)$, draw a uniform pair $(u,w)$ conditioned on $\gcd(u,w,T(m))=1$. Accept if $u$ is a unit. Locally the acceptance probability is $(1-1/r)/(1-1/r^2)=1/(1+1/r)$. The product is $\eta(m)/\eta(N)$.

The accepted mass is $v^{-h}\eta(m)\prod_{r\mid N}(1-r^{-h})/\eta(N)$. Its dependence on $g$ is exactly $g^{k-2}\eta(m)$ on either side of two. Restart on rejection. The $v=1$ term alone gives an acceptance probability at least $\prod_{r\mid N}(1-1/r)/\eta(N)\ge (\varphi(N)/N)^2$. Thus the expected number of proposals is $O((1+\log\log(3N))^2)$ by the standard totient bound. The primitive pair in step 2 takes at most $\zeta(2)$ draws on average. Modulus one is permitted throughout.

After choosing $g$, sample a uniform primitive point of squared norm $m^2$, multiply it by $g$ and draw each scalar as $m$ times an independent uniform residue modulo $g$. Every tuple in the fibre now has equal probability. [Factor-free implementation](four_square_incidence.py)

**Primitive point subroutine.** For an odd modulus $s$, find $a^2+b^2=-1\pmod s$ by the auxiliary-prime method in Pollack–Treviño. The code draws $u,v$ and puts $t=-(u^2+v^2)\bmod s$. When $t$ is a unit and we find a verified representation $t=z^2+w^2$, dividing $u+iv$ by $z+iw$ gives the required pair. The case $t=1$ also works. The verified Cornacchia calculation succeeds frequently enough on primes $t\equiv1\bmod4$. The paper's Lemma 3 gives $s\prod_{r\mid s}(1-(\frac{-1}{r})/r)$ preimages for each unit $t$. The prime number theorem modulo four therefore gives an expected polynomial number of trials. No primality assertion about an accepted $t$ is needed.

This pair identifies the Hurwitz quaternion order modulo $s$ with $M_2(\mathbb Z/s)$. A uniform pair of jointly primitive residues spans a uniform projective line. Take the preimage of the right ideal of matrices with columns on that line, using a left Euclidean gcd with $s$. Its index in the order is $s^2$, so its generator has norm $s$. Locally the columns of the generator span a direct-summand line, which makes the generator primitive. Conversely a primitive quaternion of norm $s$ has local Smith invariants $(1,0)$ modulo each exact prime power dividing $s$, giving back its line. Each ideal has 24 generators, obtained by right multiplication by the Hurwitz units. Eight have integer coordinates. Indeed every quaternion is associated to an integral one, and for an integral quaternion of odd norm the eight integral units preserve integrality while the sixteen half-integral units give half-integral coordinates. Choosing among the eight yields the required uniformity. Euclidean division shrinks the norm by a constant factor using the nearer of the integer and half-integer grids. [Primitive sampler](four_square_sampler.py) · [Pollack and Treviño, especially Lemmas 2–4](https://math.colgate.edu/~integers/sjs15/sjs15.pdf)

## Why two scalars recover factoring

At $k=2$ the power of $g$ cancels. We can convert a draw into an **exact uniform divisor** of $N$. Draw a uniform jointly primitive residue pair modulo $m=N/g$. Retain $g$ if the first coordinate is a unit. The acceptance probability $1/\eta(m)$ cancels the remaining weight and costs $O(1+\log\log(3N))$ trials in expectation. At least a third of the divisors of every composite integer are proper and nontrivial. Recursion with primality testing gives complete randomized polynomial-time factoring, handling even factors separately. For distinct odd semiprimes the direct probability in the table follows from the two local inclusion probabilities, even without rejection.

Conversely, given the factorization, choose the valuations of $g$ independently. At $r^a\parallel N$ give each exponent below $a$ weight $r+1$ and exponent $a$ weight $r$. Generate the point and the two scalars from the chosen content as above. This gives an exact sampler with the factors supplied.

**Trying to obtain the joint sample from efficient conditionals.** Fix the first scalar $t$ and let $D=\gcd(t,N)$. Write $x=(N/D)y$. Conditional on $t$, draw $(y,z)$ from the one-scalar sampler at $D$, then lift $z$ uniformly to the second scalar modulo $N$. This conditional update is public and exact. Alternating the choice of fixed scalar while resampling the point and the other scalar preserves the desired two-scalar law. From content one both scalars are zero. Each update then has the one-scalar law at $N$, with probability only $O(1/R)$ of leaving content one for balanced semiprimes. After $o(R)$ steps the probability of still having content one tends to one, while its stationary probability tends to $1/4$. This rules out a polylogarithmic mixing guarantee for that update from this starting class. A public way of generating the joint law at $k=2$ remains open. [Conditional implementation](four_square_incidence.py)

**Finite checks.** Enumerating all 170,488 signed ordered points for odd $N\le45$ matched all 115 exact content distributions for $k=0,1,2,3,4$, including the factor-supplied control at two. The script also checked all 529 fixed-scalar conditional laws in that range. At $N=493=17\cdot29$ the exact proper-content probabilities for $k=1,2,3$ are $48/589$, $1032/2065$, and $23808/267397$, respectively. [Reproducible check](check_four_square_incidence.py)
