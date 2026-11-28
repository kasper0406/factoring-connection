# Recovering a prime from yes/no questions about dyadic binomials

For every odd $N>1$, at most $\operatorname{bitlength}(N)-2$ adaptive yes/no answers recover its **largest prime divisor**. The only permitted question here is whether the whole $N$ divides a product of selected dyadic factorial ratios raised to a public power. The decoder receives one bit per call. For a product of two distinct primes, queries using at most two ratios suffice, without the power. On balanced semiprimes the two-ratio method takes at most $2\operatorname{bitlength}(N)$ calls in two rounds. Computing the answers in polynomial bit time from $N$ remains open. [Decoder](cofactor_zero_queries.py) · [Reference evaluators and checks](check_cofactor_zero_queries.py) · [Prime valuation identity](cofactor-bit-binomials.md)

## What a single prime would see

Write $b=\operatorname{bitlength}(N)$, $j=1,\ldots,b-2$, $t_j=\lfloor N/2^j\rfloor$, and $F_j=t_j!/\lfloor t_j/2\rfloor!^2$. Each ratio is a central binomial coefficient times either 1 or $t_j$. For a prime $s\mid N$, form the integer
$$
B_s=\bigvee_{a\ge1}\lfloor N/s^a\rfloor,
$$
where $\vee$ denotes binary bitwise OR. Legendre's valuation formula gives $s\mid F_j$ exactly when bit $j$ of $B_s$ is one. Bit zero of $B_s$ is always one because $N/s$ is odd. The higher quotients may add low bits to $N/s$. Every added bit is below the next power of two above $\lfloor N/s^2\rfloor$. Consequently
$$
0\le B_s-N/s<2N/s^2,\qquad 0\le s-N/B_s<2 .
$$
If $\lfloor N/s^2\rfloor=0$, there are no added bits. Thus $s$ is the unique odd integer among $\lceil N/B_s\rceil$ and its successor. Distinct odd primes are at least two apart, so the masks decrease strictly as $s$ increases. This rounding step handles contamination by all higher prime powers.

## Reading the smallest mask with Boolean queries

On a subset $J$ of the bit positions, ask whether
$$
N\ \bigm|\ \left(\prod_{j\in J}F_j\right)^b .
$$
The answer is yes exactly when $J$ hits the support of *every* prime of $N$. The power saturates its multiplicities. Begin with $J$ empty and process the positions from highest to lowest. Keep a tentative position in $J$ exactly when the answer stays no. A no-answer means at least one prime has zeros throughout $J$. Feasible $J$ are therefore the union of subsets of the complements of individual prime masks. Greedy maximization gives the numerically largest feasible $J$, the complement of the smallest individual mask. Restore its known units bit, divide $N$ by this decoded mask, and use the rounding rule above. The returned prime is the largest one. Remove its powers by exact division and repeat on the remaining modulus to get the entire factorization in $O(b^2)$ Boolean calls. Powers of two are removed directly.

## Restricting each question to two ratios

For $N=pq$, $p<q$ distinct, query every singleton. Together with bit zero, the positive answers give the intersection $C=B_p\mathbin{\&}B_q$. Here $B_q=p$ exactly, since $q^2>N$. If one mask contains the other, their intersection is $p$. Test division by $C$ first. In the remaining case, both masks have at least one exclusive position. Outside their intersection, a pair of positions gives a positive answer precisely when one belongs exclusively to each mask. The graph of positive pairs is complete bipartite between the two exclusive parts. An incident vertex and all its neighbors recover one entire mask by adjoining its neighbors to $C$. The reciprocal rounding rule recovers a prime. The complete nonadaptive list has $O(b^2)$ singleton and pair queries. Squarefreeness makes the saturating power unnecessary.

If $q<2p$, then $B_p=q$ as well, because the only possible extra bit is bit zero. Write $p=C+u$, $q=C+v$, where the nonshared bits $u,v$ are disjoint. If $C$ did not already divide $N$, both $u,v$ are nonzero with different 2-adic valuations. The least set bit of $N-C^2=C(u+v)+uv$ is exactly the least set bit of $u+v$, since $C$ is odd. This gives an exclusive pivot from the singleton answers alone. Only the pairs incident to it have to be queried.

For example, at $N=77$ the singleton answers give $C=7\mathbin{\&}11=3$. The least set bit of $77-3^2=68$ is $4$. The positive cross-pair joins bit weights $4$ and $8$. The resulting mask $3+8=11$ gives $77/11=7$. Each query was just a whole-modulus zero test.

## An exact evaluator for dense subsets

There is a native algorithm for a specified query with only $k$ unselected bit positions. If the answer is no, at least one prime mask is supported entirely on those $k$ positions plus the known units bit. Enumerate their $2^k$ submasks, adjoin one, and use the reciprocal rounding rule to get a single possible odd prime per submask. Check primality, divisibility and the actual mask by successive integer division. If every candidate fails the answer is yes. With a polynomial-time primality test the bound is $2^k\operatorname{poly}(\log N)$. The Python default uses a deterministic test through 64 bits and Lucas certificates above that range. It is exact but has no general polynomial-time certificate-generation bound. [Implementation](cofactor_zero_queries.py)

On a balanced 256-bit fixture, one selected group omitted only bit positions 121, 125 and 127. The input-only evaluator certified nondivisibility in at most eight candidate tests, and certified divisibility on selecting every position. These particular positions were chosen using the test factors. Across 3,413 small queries, the complement evaluator agreed with direct modular products on the entire mask domain. Early greedy queries, and singleton or pair queries at large $N$, leave many positions free, so the dense-query bound leaves their evaluation cost open.

## Previous scope and verification

The first singleton also reads the earlier sum-of-two-squares bit. For balanced $N=pq\equiv1\bmod4$, $N\mid F_1$ precisely when both $p,q$ are $3\bmod4$. The opposite answer means both are $1\bmod4$ and $N$ is a sum of two squares. Producing enough higher Boolean answers from the input is the unsolved part. [Earlier sum-of-two-squares bit](../../notes/eps_results.md)

Classical factorial-divisibility thresholds already give factoring reductions. For squarefree $N$, the least $k$ with $N\mid k!$ is its largest prime factor. On balanced semiprimes a variable-cutoff binomial threshold was already in the repository. The buildbit study also has a binary zero-test decoder with adaptive numerators under an additional residue promise. The distinction here is the fixed list of dyadic indices, products of at most two of them on semiprimes, and extraction of an arbitrary largest prime from a one-bit group test even with repeated factors. This analysis does not reduce the known native cost to evaluate these factorial ratios. [Earlier thresholds](../doubling/ledger.md) · [Earlier zero-test decoder](../buildbit/ledger.md) · [Classical binomial valuations](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/binomial.html)

Exact direct products recovered every factorization from 1 to 2,500. All 491 odd distinct-semiprime cases in the range were also split using only singleton/pair Boolean answers. A separate 311-bit control supplied the expected answers using known factors and verified reconstruction in 688 calls. The larger control exercises the decoder only. [Reproduce the checks](check_cofactor_zero_queries.py)
