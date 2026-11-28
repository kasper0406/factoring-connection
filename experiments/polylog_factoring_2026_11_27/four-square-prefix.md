# Including zero changes the complexity of approximate prefix counts

Let $N$ be odd and $0\le B<N$. Count the pairs $(x,t)$ with $x\in\mathbb Z^4$, $\sum_j x_j^2=N^2$ and $t x_j=0\pmod N$ for all $j$. Order and signs of the coordinates matter. Let $C_+(B)$ count $0\le t\le B$ and $C_-(B)$ count $1\le t\le B$.

| Status | Query | Result |
| --- | --- | --- |
| Proved and implemented | Approximate $C_+(B)$ from $N,B,\epsilon$ | A deterministic lower estimate $L$ with $L\le C_+(B)\le(1+\epsilon)L$, in time polynomial in $\log(2N)$ and $1/\epsilon$. |
| Proved and implemented | Approximate $C_-(B)$ within 25% for requested cutoffs | An oracle with that guarantee finds the least prime factor of every odd composite in $O(\log N)$ calls. Given the factorization, a deterministic fully polynomial approximation scheme supplies the oracle. |

The earlier result samples the unrestricted one-scalar space uniformly. The prefix problem asks for the size of a specified subset. Earlier repo work already identified the first nonunit index and the role of restricted approximate counting in guided search. Here the same integer count has both an input-only approximation algorithm and a factoring-equivalent variant, depending on one endpoint. [Unrestricted sampler](four-square-incidence.md) · [First surviving index](../../BRIEF.md) · [Earlier guided-search audit](../../notes/navigate_results.md) · [Implementation](four_square_prefix.py)

## Counting each scalar

For a fixed $t$, put $D=\gcd(t,N)$. Then $x=(N/D)y$ with $\sum y_j^2=D^2$. Jacobi's formula gives exactly $8\sigma(D^2)$ choices, where $\sigma$ sums positive divisors. At $t=0$ this gives $8\sigma(N^2)$ points at once.

Define $\eta(d)=\prod_{r\mid d}(1+1/r)$ over primes and $w(d)=d^2\eta(d)$, including $w(1)=1$. We have $\sigma(D^2)=\sum_{d\mid D} w(d)$, by checking prime powers. Exchanging the two sums gives

$$C_-(B)/8=\sum_{d\mid N}w(d)\lfloor B/d\rfloor,\qquad C_+(B)/8=\sum_{d\mid N}w(d)(1+\lfloor B/d\rfloor).$$

Gcd-sum expansions of this kind are classical. The algorithms below use the size of the boundary contribution in this particular sum. [Jacobi's four-square formula](https://doi.org/10.1017/S144678870002440X) · [Tóth, weighted gcd sums](https://cs.uwaterloo.ca/journals/JIS/VOL14/Toth/toth9.pdf)

## A deterministic public approximation when zero is included

Set $M(N)=\sum_{d\mid N}d\eta(d)$. An exact decomposition into nonnegative terms is

$$C_+(B)/8=B M(N)+\sum_{d\mid N}d\eta(d)\bigl(d-(B\bmod d)\bigr).$$

Use bit length $b$, trial-divide through $L=\lceil12b/\epsilon\rceil$, and enumerate the ordinary integers $e\le K=\lceil4/\epsilon\rceil$. These two short ranges suffice.

1. The factor of $M(N)$ at $r^a\parallel N$ is $r^a+2(r^a-1)/(r-1)$. Include this factor at each known prime. For the unprocessed cofactor, its leading factor is the cofactor itself. Omitting the remaining correction factors changes $M(N)$ by a ratio at most $\exp(3b/L)$, since each is at most $1+2/(r-1)$ and at most $b$ distinct primes remain.
2. In the second sum retain $d=N/e$ for every $e\le K$ dividing $N$. Compute $d\eta(d)$ using the known primes of $d$ and omit the other factors $1+1/r$. The same ratio bound applies to each retained term. All arithmetic in the implementation is integral.
3. Before rounding the prime factors, the omitted terms weigh at most $N^2\eta(N)\sum_{e>K} e^{-2}\le N^2\eta(N)/K$. In the original zero-inclusive count, $d=N$ alone contributes $N^2\eta(N)$ to $C_+(B)/8$ for every cutoff. The missing terms therefore weigh at most a $1/K$ fraction of the target.

The returned count is at least $(1-1/K)\exp(-3b/L)C_+(B)\ge C_+(B)/(1+\epsilon)$ for $0<\epsilon\le1$. The last inequality uses $(1-\epsilon/4)\exp(-\epsilon/4)\ge1-\epsilon/2\ge1/(1+\epsilon)$. Operands have polynomial bit length. The work includes $O(b/\epsilon)$ trial candidates and $O(1/\epsilon)$ potential complementary divisors, without enumerating the divisors of the unfactored cofactor.

## Why removing zero recovers a factor

Let $p$ be the smallest prime factor of the odd composite $N$. For $1\le B<p$ every scalar is a unit, so $C_-(B)=8B$. For $B\ge p$ each multiple of $p$ contributes at least $8(p^2+p+1)$ instead of the baseline eight, giving

$$\frac{C_-(B)}{8B}\ \ge\ 1+\frac{(p^2+p)\lfloor B/p\rfloor}{B}\ \ge\ 1+(p+1)/2\ \ge\ 3.$$

Thus an estimate within 25% of the true count is below $12B$ before the first factor and above $12B$ from then on. Binary search through $[1,\lfloor\sqrt N\rfloor]$ gives $p$, also for repeated or unbalanced factors. Verify the divisor and recurse, stripping even factors and testing primality as usual. For bounded-error estimates, amplify confidence across the polynomial number of calls.

**Reverse direction with the factors supplied.** Process each $r^a\parallel N$ with the recurrence $F_i(u)=\sum_{j=0}^a w(r^j)F_{i-1}(\lfloor u/r^j\rfloor)$, starting from $F_0(u)=u$ for the positive prefix. Compress each nondecreasing function into a step table. Keep its initial value, then binary-search for the first argument where the value exceeds the previous retained one by a factor $1+\delta$. The first positive value is kept even when there is an interval of zeros. All thresholds act on values, so floor discontinuities in the argument are preserved to the stated relative accuracy. The use of lower tables and nonnegative weights accumulates a factor at most $(1+\delta)^h$ over $h$ primes. Taking $\delta=\epsilon/(4\max(1,h))$ suffices. There are $O(1+\log(2N)/\delta)$ retained values per stage since all counts have $O(\log(2N))$ bits. The time is polynomial in the input bit length and $1/\epsilon$. This uses the monotone-table technique also used for counting knapsack solutions. [Halman's approximation scheme](https://pluto.huji.ac.il/~halman/Papers/CountingTCS.pdf)

## Why subtracting the zero count loses the guarantee

At $N=1009\cdot1013=1{,}022{,}117$, the positive prefix is 8,064 at $B=1008$ and 8,160,792 at $B=1009$. Including zero adds 8,374,343,605,224 to **each** count. The first-factor jump is exactly $1/1{,}027{,}183$ of the zero contribution. At 25% tolerance the public algorithm returns the same value, 8,357,785,293,512, for both cutoffs and at $B=0$. Subtracting its two estimates gives zero across the jump.

More generally, for growing balanced semiprimes, both primes eventually exceed every fixed polynomial in the bit length. At any fixed inverse-polynomial tolerance, the algorithm then sees just $e=1$ in its correction and uses $M(N)\approx N$. The two terms simplify to $BN+N(N-B)=N^2$ for every cutoff. A way to control the *relative error of the positive part* across cutoffs would give a factoring algorithm. [Reproducible checks](check_four_square_prefix.py)

**Verification.** Run `python experiments/polylog_factoring_2026_11_27/check_four_square_prefix.py` from the repository root. It enumerates 170,488 signed ordered points for odd $N\le45$, checks 3,174 estimates at 529 cutoffs, and tests larger examples where trial division leaves an unfactored cofactor. The separate oracle test recovers the least factor on 70 odd composites under three error patterns, including one that directs every 25% error toward the decision threshold. Its test oracle uses known factorizations.
