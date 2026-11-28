# A dyadic factorial ratio reads bits of the hidden cofactor

Let $N=pq$ for distinct odd primes $p<q<2p$. For every public integer $c\ge2$, put $t=\lfloor N/c\rfloor$, $u=\lfloor t/2\rfloor$, and
$$
 F_c(N)=\frac{t!}{(u!)^2}=
 \begin{cases}
 \binom{2u}{u}&t=2u,\\
 t\binom{2u}{u}&t=2u+1.
 \end{cases}
$$
Then the prime valuations are exactly
$$
v_p(F_c(N))=\lfloor q/c\rfloor\bmod2,\qquad
v_q(F_c(N))=\lfloor p/c\rfloor\bmod2 .
$$
Thus at $c=2^j$ each prime sees the $j$th binary digit of the *other* prime, with the units digit numbered zero. For $N=77$ the tests at $c=4,8$ give gcds 11 and 7 respectively. This calculation arose from asking whether signed periodic factorial products retain information that a long, unweighted product loses. [Previous full-product test](collective-norm-product.md)

For any integer $N$ and any prime divisor $s$, Legendre's factorial valuation and nested integer division give
$$
 v_s(F_c(N))=\sum_{a\ge1}
 \left(\left\lfloor\frac{N}{c s^a}\right\rfloor-
 2\left\lfloor\frac{N}{2c s^a}\right\rfloor\right).
$$
Every summand is zero or one. In the balanced case $N/(c p^2)<1$, so only $a=1$ can contribute. Taking $c=2^j$ in general, each term is bit $j$ of $\lfloor N/s^a\rfloor$. This is a specialization of the classical prime valuations of factorials and binomial coefficients. [Granville on Legendre, Kummer and Lucas](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/binomial.html)

## The same menu separates every pair of odd prime divisors

The publicly specified menu of $F_{2^j}(N)$, $1\le j<\operatorname{bitlength}(N)$, separates **all distinct prime divisors of any odd $N$**. Unequal prime sizes and repeated prime factors are allowed. For a prime divisor $s$, its pattern of positive valuations is the bitwise OR of the integers $\lfloor N/s^a\rfloor$ for $a\ge1$, with the unused units bit dropped.

To prove separation, take any two prime divisors $p<q$ and set $t=N/p$, $u=N/q$. Let $j$ be the most significant bit in which $t$ and $u$ differ. These are odd integers, so $j\ge1$. Bit $j$ of $t$ is one, while bit $j$ of $u$ is zero. Moreover
$$
 2^j>\frac{t-u}{2}=\frac{N(q-p)}{2pq}>\frac{N}{q^2},
$$
where $q-p\ge2$. Every higher-power term at $q$ therefore vanishes in this test. At $p$ any additional terms only increase the valuation. Consequently $p\mid F_{2^j}(N)$ and $q\nmid F_{2^j}(N)$. Above this index $q$ has no higher-power contributions, so any change to equal higher cofactor bits can only favor $p$.

Given the modular residues of this nonadaptive menu, there is deterministic polynomial-time postprocessing. For each residue $f$, compute $g=\gcd(N,f^{\operatorname{bitlength}(N)})$ by modular powering. It contains full prime powers for exactly the prime divisors of $f$. Refine the current list of coprime divisors using gcds with $g$ and exact division. Pairwise separation leaves one prime power per block. Exact integer roots at the possible exponents then give its prime and multiplicity. An input with a single prime divisor is already one such block. Powers of two in an even input can be removed first. The code constructs the binomials directly for small checks. Polynomial evaluation of these residues on unbounded input remains open. [Reproducible checker and decoder](check_cofactor_bit_binomials.py)

Combining all tests with their natural dyadic weights loses the separation:
$$
 \prod_{j\ge1} F_{2^j}(N)^{2^{j-1}}=\lfloor N/2\rfloor! .
$$
The identity telescopes since the denominator factorial in one ratio is the next numerator. This helps explain why a computation of the combined product alone does not give the individual divisibility bits.

## Relation to earlier binomial tests

The repository already showed that a central binomial at a square-root-scale public index can split balanced semiprimes, and that locating a factorial or binomial divisibility threshold suffices for factoring. An adjacent Gaussian CM coefficient also distinguishes mixed prime classes modulo four. Here one fixed geometric menu gives the separate binary digits on balanced inputs and pairwise prime-power separation on arbitrary odd inputs. This strengthens the description of the observable. A fast general evaluator for even one useful binomial was already the open step. [Square-root-index coefficient](../../notes/direct_split_poly.md) · [Doubling and threshold experiments](../doubling/ledger.md) · [Adjacent Gaussian CM coefficient](cartier-paths.md#the-coefficient-beside-the-gaussian-cm-middle-term)

On balanced inputs, if none of the queries for $1\le j<k$ gives a proper gcd, then $p,q$ have the same $k$ low bits, including their common odd bit. Thus $p^2\equiv N\bmod 2^k$. For $k\ge3$ there are four possible odd roots when a solution exists. They can be enumerated by successively lifting and testing the two children of each root. Each candidate gives one affine progression for the unknown divisor. The earlier certified small-root method applies once the known modulus exceeds the quarter-power threshold with its stated margin. This observation reduces how many differing bits are needed from a hypothetical evaluator. It leaves the cost of its largest binomial queries. [Earlier affine divisor recovery](list-order-window.md)

Tests checked 720 balanced valuations at dyadic and other cutoffs and recovered the complete factorization of 604 odd inputs, including unbalanced cases and repeated prime powers. They also verified 544 ordered prime-pair distinctions and the telescoping equality through the tested small range.
