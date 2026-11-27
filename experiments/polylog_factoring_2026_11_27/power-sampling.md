# A direct gcd matches power sampling for the square residue decoder

Let $N=pq$ for distinct odd primes $p<q$. Suppose a trial takes independent uniform units $u_i$ and forms $v=\prod_i u_i^{e_i}\bmod N$ for public exponents, then applies the square-residue decoder to $z=cv$ for a public unit $c$. The direct test $\gcd(v-1,N)$ succeeds with at least half the probability of the square-residue test whenever the power image is nontrivial. The bound places no restriction on the degrees of the power maps. If the power image is trivial, the candidate is always the public value $c$.

## Successful square residues form a forest

Put $B=\lfloor\sqrt N\rfloor<q$. Consider a bipartite graph with vertices the unit residues modulo $p$ on the left and modulo $q$ on the right. Each unit $z$ which decodes to a proper factor is an edge, indexed by its two reductions. Such an edge has mixed signs relative to some ordinary integer $1\le t\le B$ with $\gcd(t,N)=1$. At a right vertex $j\in[1,q-1]$, the only possible $t$ are $j$ and $q-j$. The respective left endpoints are $-j$ and $q-j$ modulo $p$.

Any right vertex of degree two thus has $q-B\le j\le B$. This interval contains at most $p-1$ integers, since $p+q-2B\ge2$. The residues $-j\bmod p$ of those vertices are distinct. Contract each degree-two right vertex to an edge from $-j$ to $-j+q$ on the left. These edges are a proper subset of the length-$p$ cycle obtained by repeatedly adding $q$ modulo $p$. Degree-zero and degree-one right vertices cannot create cycles. Therefore the successful-edge graph is a forest.

## Compare probabilities

The product of powers is a group homomorphism on the input units. It has uniform output on its image. By CRT and independent local inputs, this image is the product of two local subgroups, of sizes $a$ and $b$. Multiplication by $c$ gives local cosets of the same sizes. The subgraph on their $a+b$ vertices is also a forest, so at most $a+b-1$ of the $ab$ candidates decode.

The direct gcd succeeds precisely when one local component of $v$ equals 1 and the other differs from 1. There are $a+b-2$ such pairs, and their probability is $(a+b-2)/(ab)$. Thus for $a+b>2$ the ratio of square-decoder success probability to direct-gcd success probability is at most
$$1+\frac1{a+b-2}\le2.$$
If $a=b=1$, then $v=1$ and the square decoder can be applied directly to $c$. Neither the local sizes nor the factorization are needed to run both the gcd and the check on $c$. A nonunit in the original samples or scale can be checked by gcd before using the unit comparison.

This reduction also applies to any sampler that is uniform on a Cartesian product of local sets, if it supplies a public anchor in that product. Replace $v-1$ by the difference between the sampled residue and the anchor. The proof leaves general high-degree maps and CRT-dependent sampling open. A successful exponent choice for the power sampler would already give a comparably successful direct order test.

## Checks

The check script exhausts exponents over one full period and all unit scales on six semiprimes. All 81,108 nonconstant comparisons obeyed the bound, and the largest observed ratio was 1. Run `python check_power_sampling.py`.
