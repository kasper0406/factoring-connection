# Generating the smaller prime from two fixed binomial readings per level

Let $N=pq$ with distinct odd primes $p<q<2p$. Choose the indices from $N$ alone. With $h=1,3$ put $t_{h,j}=\lfloor hN/2^{j+1}\rfloor$ and request the Boolean value $[N\mid {2t_{h,j}\choose t_{h,j}}]$. On a density-one set of balanced prime pairs, these readings permit deterministic polynomial-time recovery of the smaller prime. This includes generating all surviving prefixes and verifying the divisor by integer division. The runtime for producing the readings is additional. One may also use $h=1,5$ and request the second multiplier on only three of every five bit positions. [Implementation](fixed_binomial_streams.py) · [Exact carry proof](two-fixed-and-streams.md) · [Sparse schedule and its limitations](sparse-fixed-multipliers.md)

## Choosing the depth and finishing by division

Set $L=\lceil\sqrt{N/2}\rceil$, $U=\lfloor\sqrt N\rfloor$ and $k=\lfloor\log_2L\rfloor$. Start the first multiplier at $j=1$, and the second when $2^{j+1}\ge4h$. Stop both at $j=k-1$. Every query satisfies $2t_{h,j}<p^2$ and $2^j<p$. Thus the bit is exactly bit $j$ of $(hp)\mathbin{\&}(hq)$ by Kummer's theorem and the rounding calculation in the carry proof. All queries are specified before the first answer. The sparse multiplier-5 schedule reads levels with $j\bmod5\in\{1,3,4\}$. [Binomial valuations](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/binomial.html)

An admissible odd $x\bmod 2^j$ has exactly two lifts as a first coordinate. At the next modulus $D=2^{j+1}$ its partner is $y=Nx^{-1}\bmod D$. The two candidate pairs are $(x,y)$ and $(x+2^j,y+2^j)$ modulo $D$. Keep a child precisely when its predicted bits agree with the readings. The true ordered pair survives independently of all complexity estimates. At the last depth each surviving class has at most one representative in the factor interval, because
$$
U-L\le(\sqrt2-1)L<2^k.
$$
That representative is $L+((x-L)\bmod 2^k)$ if it is at most $U$. A division tests it. This procedure works on every input satisfying the promise if one allows the exact search to finish.

## What “almost all” means at full factor depth

Choose two primes uniformly from $[V,2V]$, conditioned on being distinct, and sort them. There are $\asymp V^2/(\log V)^2$ ordered distinct pairs before sorting. The input-dependent $k$ has only a bounded number of possible values, each with $2^k\asymp V$. For each such depth a residue pair modulo $2^k$ has only $O(1)$ integer representatives in the dyadic square. Therefore the probability of any exceptional prefix-tree event at that depth in the prime sample is at most $O((\log V)^2)$ times its probability for independent uniform odd endpoints. Apply this separately to each possible depth before summing. Sorting changes an upper bound by at most a constant.

For the full second stream the uniform-model list bound is $\Pr(F_j>S)\ll S^{-\eta}$ with $\eta=\log(10/9)/\log4$. A union over at most $k$ levels and the preceding prime-count factor give failure $O(k^{3-42\eta})=o(1)$ at width $k^{42}$. This exponent is a conservative certificate. A larger fixed exponent works for the sparse schedule, using its per-level tail. These are distributional statements about balanced prime inputs. The carry proof exhibits exponentially large fibers for particular low-product observations even using every member of a fixed finite multiplier menu.

## Checks and comparison

There are 199 balanced pairs among distinct odd primes below 137. For both settings the tests used literal integer binomials, 2,160 readings across 398 runs. They also checked the generated width at every level against independent exhaustive grouping. For larger inputs, the table uses 16 seeded prime pairs per size. These readings were computed with the known test factors by independent factorial valuations. They measure only the decoder. [Reproduce the checks](check_fixed_binomial_streams.py)

| Bits per prime | Full multiplier 3, median / maximum peak width | Sparse multiplier 5, median / maximum peak width |
|---:|---:|---:|
| 24 | 9 / 30 | 14 / 38 |
| 48 | 21 / 328 | 44 / 3,414 |
| 96 | 35 / 358 | 59 / 590 |
| 128 | 41 / 216 | 125 / 466 |

Earlier work here recovered factors on all balanced inputs with adaptively chosen binomials. One fixed multiplier left exponentially many local candidates on typical inputs. Adding a second fixed multiplier permits a nonadaptive reduction and an elementary finisher. Bitwise lifting with supplied factor or key bits also has a substantial prior literature. Here the observations are deterministic joint bits at prescribed integer multiples, encoded in specified binomial tests. The unresolved algorithmic target is to produce these bits, or a sufficiently accurate version of them, from $N$. [Adaptive reduction](adaptive-binomial-lift.md) · [Single multiplier](shared-and-fiber-entropy.md) · [Reconstruction from random key bits, Heninger and Shacham](https://hovav.net/ucsd/papers/hs09.html) · [Noisy joint readings](noisy-fixed-streams.md)
