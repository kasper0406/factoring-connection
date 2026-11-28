# One binary window of square roots of unity encodes a subset count

Given positive binary weights $w_1,\ldots,w_k$, there is a randomized construction of an odd **squarefree** $N$, its full factorization and a power of two $D$. For every $1\le T\le S=\sum_i w_i$, simultaneously,

$$
\#\{x\in[TD,(T+1)D):x^2\equiv1\bmod N\}
=\#\{A\subseteq\{1,\ldots,k\}:\sum_{i\in A}w_i=T\}.
$$

The intervals are aligned on the ordinary integer grid, lie in the lower half and exclude $\pm1$. The construction is always correct and has expected polynomial bit cost. Thus exact counting is #P-hard, and parity counting is ⊕P-hard, under randomized expected-polynomial-time count-preserving reductions. The earlier construction chose the right-hand side of the quadratic congruence. Here it is fixed to 1 for every input. The positive-weight version of #SubsetSum already admits the usual parsimonious reduction from #3SAT. Targets outside $[0,S]$ have count zero, and $T=0$ has count one. [Earlier variable-target construction and complexity reduction](../../notes/floorhard_results.md#21-arbitrary-intervals)

## Programming the CRT coefficients

Take $L$ to be the least power of two at least $16(S+1)$. Choose $k$ distinct odd primes $q_i>32kL$ and put $Q=\prod_iq_i$. The next section specifies their sizes and a polynomial-time sampling procedure. Set $b_i=\lceil q_iw_i/L\rceil$. These are units modulo $q_i$ with

$$0\le b_i/q_i-w_i/L<1/q_i.$$

Choose a prime $P$ in the unit residue class given by

$$P\equiv \big[b_i(Q/q_i)\big]^{-1}\pmod{q_i}\qquad(1\le i\le k).$$

Choose its size so that $H<N=PQ\le H(1+1/(4S))$ for a power of two $H\ge L$. Then $D=2H/L$ is an integer power of two. The CRT idempotent for $q_i$ has representative $e_i=(N/q_i)b_i$. Indeed it is 1 modulo $q_i$ and 0 at every other prime of $N$.

Start with the root $-1$. Keep the sign at $P$ negative and flip the signs at $q_i$ for $i\in A$. For a nonempty subset of weight $t$, the resulting integer is

$$z_A=-1+2\sum_{i\in A}e_i,\qquad
2Nt/L-1\ \le z_A < 2Nt/L-1+N/(16L).$$

The lower bound exceeds $tD-1$, so $z_A\ge tD$ by integrality. For the upper bound use $2Nt/L\le (t+1/4)D$ and $N/(16L)\le(1+1/(4S))D/32$. Hence $z_A<(t+1)D$. All these integers are canonical roots since $(S+1)D\le H/8<N/8$. Distinct subsets give distinct sign vectors. The empty subset gives $N-1$, and changing the sign at every prime accounts for all roots with positive sign at $P$. Those reflected roots lie above $7N/8$ except for 1. This proves the count identity.

There is also a two-query consequence for the translated half-window evaluator. Let $m=(N-1)/2$ and $W(j)$ count roots in $[j,j+m)$. For $u=TD$, $v=(T+1)D$, the difference $W(u)-W(v)$ counts the roots in $[u,v)$ minus those in $[u+m,v+m)$. The latter interval lies in the empty middle of this construction. The identity holds modulo two and for mixed roots alone. [Fixed-length factoring reduction](fixed-window-split.md)

## Finding the auxiliary prime in expected polynomial time

Use two independent choices of core primes from **disjoint** sets. For example, take one prime from each of $2k$ disjoint intervals $[2^j B,2^{j+1}B)$ for $B=32kL$, $1\le j\le2k$. Sample uniform integers and test primality. The prime number theorem gives an expected $O(\log(2^jB))$ trials per interval. All candidates and the two core products $Q_1,Q_2$ have polynomial bit length. Compute a CRT class $a_i$ for each core as above. Either core gives the same subset counts.

Put $R=\max(Q_1,Q_2)$, $h=\lceil\log_2(2R)\rceil^3$ and $H=2^h$. We need a prime in just one of the two progressions, in the interval

$$H/Q_i < P\le (1+1/(4S))H/Q_i.$$

The endpoints are at most $H$. We use the prime number theorem for progressions with a possible exceptional real character. In the quantitative form of Drappeau and Fiorilli, Lemma 2.2, the error, uniformly for all cutoffs $y\le H$ and reduced residue classes for $Q_i\le R$, is at most

$$E=O\!\left(H e^{-\delta\sqrt{\log H}}+R\sqrt H(\log H)^C\right)$$

after accounting for the possible exceptional term, with absolute fixed constants. Its conductor exceeds one and divides at most one of the coprime $Q_i$. On the other core the main term for the von Mangoldt sum is $y/\varphi(Q_i)$ at **both** endpoints. Subtraction gives main term $H/(4S Q_i\varphi(Q_i))$ with error $O(E)$. Our choice of $h$ ensures $E=o(H/(S R^2))$, since $S<R$ and $\sqrt{\log H}$ grows faster than $\log R$. Contributions from proper prime powers at most $H$ are $O(\sqrt H\log^2 H)$, also negligible. This controls the short interval even when $1/S$ is small in terms of the input length. [Drappeau and Fiorilli, Lemma 2.2 and equations (1.4)–(1.6)](https://arxiv.org/pdf/2003.02201#page=4)

For this core the interval contains $\gg H/(S Q_i\varphi(Q_i)\log H)$ primes in the chosen class, among $O(H/(S Q_i^2))$ candidates. Select the core with a fair coin, sample a progression member uniformly from its interval and test primality. The success probability is $\Omega(1/h)$ for large inputs. Primality can be decided deterministically in polynomial bit time. For a uniform procedure, after $h^2$ unsuccessful trials double $h$ and repeat. The failure probability per sufficiently large stage is at most $\exp(-\Omega(h))$, so the doubling rule preserves expected polynomial cost and handles the finite small exceptions. Every accepted prime satisfies the exact residue and band conditions. [Agrawal, Kayal and Saxena, *PRIMES is in P*](https://annals.math.princeton.edu/2004/160-2/p12)

## Verification and implications for factoring

The diagnostic uses deterministic prime searches on six small fixtures and checks all 92 CRT sign vectors and all 24 positive target blocks. For four unit weights, a single modulus gives counts $4,6,4,1$ in consecutive aligned blocks of width $2^{111}$. The weights $(2,2,3)$ also exercise empty blocks. Exhaustive root enumeration and direct subset enumeration agree. The finite search illustrates the encoding. The expected-time guarantee uses the two-core sampler above. [Reproduce the construction](check_fixed_one_subset_encoding.py)

A polynomial-time algorithm for parity on **every** such interval would compute the parity of the number of satisfying assignments of an arbitrary Boolean formula in randomized expected polynomial time. This is a complexity implication, and an unconditional time lower bound remains open. The reduction uses a growing number of distinct primes. With two distinct odd prime factors supplied there are four roots to inspect on a semiprime.

The starting interval matters even for a search-only guide. Suppose it must run in polynomial time from any supplied aligned interval and return a root whenever the interval is occupied. Run it on $[TD,(T+1)D)$ and verify the root. This decides positive SubsetSum, since an empty interval cannot produce a valid witness. A child selector with that guarantee for arbitrary starting intervals gives the same algorithm by descent to a singleton.

Factoring needs a weaker operation. Its descent starts at the **entire lower half** and may follow whichever mixed root is easiest to obtain. With the factors supplied one can choose a mixed sign vector by CRT and immediately guide that descent, regardless of how many other roots any child contains. On an unfactored input the unresolved constructor is a guide that selects an occupied child along its own path down to the certified local search radius. The arbitrary-interval hardness proof gives no hardness result for that promise. [Path formulation and local search](README.md#working-ledger-for-the-quarter-parity-evaluator)
