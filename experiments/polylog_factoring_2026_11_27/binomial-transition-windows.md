# Choosing a multiplier or shifting a binomial query

The adaptive lift allows several odd numerators $h$ for the same next bit of a factor. Some states force **every aligned numerator** to give a binomial index of order $N$. Nearby shifts have an exact description in terms of one arithmetic progression. For the exceptional opposite states this description gives a factor directly, using a product over a shorter interval. Its length can still grow as a power of $N$. [Adaptive single-binomial lift](adaptive-binomial-lift.md) · [Earlier use of the consecutive-coefficient recurrence, §R4](../buildbit/ledger.md)

## A state in which all aligned indices are large

Take a power $d\ge16$ of two and the candidate lifts $\{1,-3\}$ and $\{1+d/2,-3+d/2\}$ modulo $d$. For each odd $0<h<d/4$, $h$ lies in the lower half while $-3h\bmod d=d-3h$. Thus alignment holds **exactly** for $d/6<h<d/4$. Adding $d/2$ to both endpoints interchanges their halves and preserves agreement. Every normalized aligned choice has $t=\lfloor hN/d\rfloor>N/6-1$. For each prescribed $d$, either lift occurs on arbitrarily large balanced semiprimes by taking distinct primes in its two odd progressions in $[X,2X]$ for large $X$. This uses the prime number theorem for fixed progressions and asserts no uniform bound on the least such $X$. The result bounds the index that an optimizer can attain within this family. A large index by itself gives no lower bound on indirect evaluation.

There is also no close dyadic index from a preceding level. Unequal reduced dyadic fractions, one at denominator $d$ and the other at a smaller power of two, differ by at least $1/d$. Their integer indices therefore differ by at least $\lfloor N/d\rfloor$. Any two different normalized odd numerators at denominator $d$ are separated by at least $\lfloor 2N/d\rfloor$ index steps. This spacing concerns the unshifted centers.

The *smallest* aligned numerator for any state $(a,b)$ can itself be found in polynomial bit time by fixed-dimensional integer programming. For a bound $H$, let $h=2u+1$ with $0<h\le H$ and $h<d/4$. For each of the two same-half choices $[A,B]=[0,d/2-1]$ or $[d/2,d-1]$, impose
$$
A\le ah-dz\le B,\qquad A\le bh-dw\le B .
$$
There are only three integer variables $u,z,w$, so Lenstra's algorithm decides feasibility in polynomial time. Binary search recovers the minimum. A proposed interval $[V,W]$ of inexpensive binomial indices can be imposed by two more linear inequalities $dV\le hN<d(W+1)$. A polynomial-sized list of such intervals is handled one interval at a time. In the family above even this exact selector cannot lower the index below $N/6-1$. [Lenstra, *Integer Programming with a Fixed Number of Variables*](https://doi.org/10.1287/moor.8.4.538)

## Exact singularities and a uniform interval of unit ratios

Write $B_u=\binom{2u}{u}$, $m=d/2$, and $e=hN\bmod d$. For $t=\lfloor hN/d\rfloor$ the exact recurrence is
$$
(t+j)B_{t+j}=2(2(t+j)-1)B_{t+j-1}.
$$
For odd $N$, powers of two are units. Since $dt=hN-e$, the two gcds on the right and left are respectively
$$
\gcd(N,m(2j-1)-e),\qquad \gcd(N,2mj-e).
$$
Successive steps sweep the two interleaved halves of **one progression** $mi-e$. To test all $K$ steps in the increasing direction, form one polynomial $\prod_{i=1}^{b}(X+di)$, $b=\lfloor\sqrt K\rfloor$. Evaluate it at the block starts of both progressions and inspect the partial end blocks. The standard product and remainder trees need $\widetilde O(\sqrt K\,\mathrm{poly}(\log N))$ bit operations for $1\le K\le N$. The code computes both product gcds with $N$. [Implementation](binomial_transition_windows.py)

Put $L=\lceil\sqrt{N/2}\rceil$ for balanced $N=pq$. For $J\ge1$ with $dJ+m\le L$, every step between $t-J$ and $t+J$ has both coefficients invertible modulo $N$. Indeed its progression indices run from $1-2J$ to $2J$, so each $mi-e$ has absolute value at most $dJ+m-1$. Each is odd and hence nonzero. As both primes are at least $L$, there are no singularities. The Boolean zero answer and its whole gcd with $N$ are constant throughout this interval. This is uniform over every odd choice of $h$. In particular a shift polynomial in $\log N$ cannot change the zero answer at any center through the balanced quarter-root stopping level for sufficiently large $N$. It could still make some other formula for that same answer available.

## A factor on both opposite lifts

For a true opposite state in the balanced algorithm take $d\ge16$, $d/2<p$, and the folded $h$ from the shifted-lift proof. The two scaled cofactor residues are either the center pair $\{d/2-1,d/2+1\}$ or the edge pair $\{1,d-1\}$. Put $K=\lfloor(\lfloor\sqrt{2N}\rfloor+1)/d\rfloor+2$. For each prime $s\mid N$, its initial units digit is $\lfloor sr/d\rfloor$ for the corresponding cofactor residue $r$. We have $s/d+1<K\le2s/d+2$ and $3s/d+2<s/2$. The last inequality uses $s>d/2$ and $d\ge16$. Inspect the next $K$ increments.

* In the center pair, the lower digit crosses $(s+1)/2$ exactly once. The upper digit stays above the midpoint and below $s$. Exactly one prime appears in the numerator progression and neither appears in the denominator.
* In the edge pair, the high digit crosses $s$ exactly once and ends below the next midpoint. The low digit stays below the midpoint throughout. Exactly one prime appears in the denominator progression and neither appears in the numerator.

The two product gcds are therefore $(s,1)$ or $(1,s)$ for a proper factor $s$ of $N$. Both outputs have arithmetic certificates determined by $N$ and the proposed state. If a false state yields a proper factor it is still usable. Any output outside the two required patterns rules out the proposed state under the same balance and size assumptions. This evaluates the shifted decision on a true opposite state while recovering a factor. The length $K=O(\sqrt N/d+1)$ is the same order as searching for a factor in a prescribed residue class in the balanced interval. The computation still does not decide a general aligned Boolean query in polynomial bit time. [Shift bounds](adaptive-binomial-lift.md#a-shifted-single-test-removes-branching-on-balanced-inputs)

The checker matched 4,867 batched pairs of gcds against literal recurrence products, verified 4,019 integer binomial shifts in the unit intervals, and checked the exact factor patterns on 40 center and 22 edge states. It also exercised 180 false state representatives and checked the full large-index multiplier set for thirteen powers of two through $65,536$. [Reproduce](check_binomial_transition_windows.py)
