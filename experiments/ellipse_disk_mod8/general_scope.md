# General integer factoring from centered counts modulo four

A deterministic algorithm can factor every positive integer given a polynomial-time evaluator for the centered counts on the explicit family below. The reduction uses primitive forms, a public unit square root of each discriminant modulo the current odd input, and scalar residues modulo four. It also works with residues modulo eight or sixteen. Its reference evaluator scans linearly many positions at the first child. An unconditional polynomial-time evaluator for even this family remains open.

The family widens the [centered decoder](centered.md). Its initial cutoff is approximately $N/2$, and its discriminants grow up to $O(N^4)$. It has no promise of a public representation of $N$. Thus the result does not transfer an evaluator restricted to the earlier near-square branch to arbitrary factoring. The [earlier floor-sum hardness work](../../notes/floorhard_results.md) uses exact interval counts of square roots of unity. Here centered, unmasked disks encode the parity of selected roots, including for moduli with repeated prime factors.

## A disk with three possible rows

For odd $N\ge3$, put $B=(N-1)/2$, $d=N(B+1)-1$, and set
$$
F(x,y)=x^2-Nxy+dy^2,\qquad \Delta=4d-N^2=N^2+2N-4.
$$
The discriminant is $-\Delta\equiv4\pmod N$. On the row $y=1$,
$$
F(x,1)-NB=(x-1)(x-(N-1)),\qquad F(x,1)\equiv x^2-1\pmod N.
$$
The eligible points on that row are exactly the roots $1\le x\le N-1$ of $x^2\equiv1\pmod N$. Their negatives give the row $y=-1$. For $|y|\ge2$, completing the square gives $F(x,y)\ge\Delta y^2/4\ge\Delta>NB$. The axis $y=0$ can contribute when $N$ has repeated prime factors.

For a power of two $m\ge2$ and an even residue $0\le r<m$, substitute $(x,y)=(mu+rv,v)$. The resulting form has coefficients
$$
[m^2,\;m(2r-N),\;r^2-Nr+d].
$$
Make a second form by replacing $d$ with $d+N^2$. This raises every nonaxis value by at least $N^2>NB$ and leaves all axis values unchanged. Normalize each form separately. If its content is $g$, then $g\mid m$ because $2r-N$ is odd. Divide all coefficients by $g$ and take cutoff $\lfloor B/g\rfloor$. This preserves the count because $g$ is invertible modulo $N$ and a counted value divided by $g$ must be an integer multiple of $N$. The resulting discriminant is $-(m/g)^2\Delta$ in the main form and $-(m/g)^2(\Delta+4N^2)$ in the comparator. In both cases $2m/g\bmod N$ is a public unit square root. The leading coefficient is a unit modulo $N$, so a root in the chart's slope coordinate is also public by the quadratic formula.

Let $R$ and $R_0$ be the corresponding scalar counts modulo four. Both count all signed nonzero integer vectors with normalized form value at most $N\lfloor B/g\rfloor$ and divisible by $N$. Then
$$
P_N(m,r)=\frac{(R-R_0-2\mathbf1_{r\equiv N-1\bmod m})\bmod4}{2}
$$
is the parity of the number of *nontrivial even* roots $x\in[1,N-1]$ with $x\equiv r\bmod m$. The indicator subtracts the public even root $N-1$. Each remaining root contributes a pair of opposite signed vectors.

Since $m<2N$ throughout the search, all normalized coefficients have $O(\log N)$ bits and the discriminants are $O(N^4)$. All cutoffs also satisfy $\lfloor B/g\rfloor\le\sqrt{|\operatorname{disc}|/3}$, by $3B^2\le\Delta$. The numerical cutoff and discriminant themselves can grow as powers of $N$. A bound stated in terms of discriminant magnitude would need to account for that growth. For nonsquare $N\ge5$, even the starting form $F$ has no representation of $N$ because its nonaxis values are at least $\Delta/4>N$ and its axis values are squares.

## Detecting prime powers and splitting all other inputs

For any odd prime power $p^e$, either $p^e\mid(x-1)$ or $p^e\mid(x+1)$ whenever $x^2\equiv1\pmod {p^e}$, since $\gcd(x-1,x+1)\mid2$. The Chinese remainder theorem gives $2^{\omega(N)}$ roots, where $\omega(N)$ is the number of distinct prime factors. Pairing $x$ with $N-x$ selects one even representative from each pair. The number of hidden even representatives is $2^{\omega(N)-1}-1$. Hence $P_N(2,0)=0$ exactly when $N$ is a prime power.

If that parity is one, maintain an even residue $r\bmod m$ with $P_N(m,r)=1$, starting at $m=2,r=0$. Its two children partition its roots. Query the left child's parity $P_N(2m,r)$ and choose it if odd. Otherwise choose $r+m\bmod 2m$. When the stride reaches $N$ there is at most one integer in the residue class between 1 and $N-1$. Its parity is odd, so it is a nontrivial root. Taking $\gcd(x-1,N)$ recovers a proper divisor and its cofactor contains the complementary prime powers. This takes $2(\lceil\log_2 N\rceil-1)$ scalar calls after the initial parity test.

For any positive integer, first remove factors of two. At each odd node with $P_N(2,0)=1$, split and recurse. At an odd prime-power leaf, find the largest $e\le\lfloor\log_2 N\rfloor$ for which $N=a^e$ by exact integer root testing. The base is its prime factor, including $e=1$. The recursion has $O(\log N)$ nodes and $O((\log N)^2)$ scalar queries. All remaining bit arithmetic is polynomial in $\log N$. This gives a full conditional factorization algorithm without an independent primality oracle.

## What the evaluator must recover

On an odd semiprime with distinct factors, there is just one hidden even root. Its dyadic residue class is exactly what the corrected parity tests. Conversely, given the factors on this three-row family, all four roots come from two CRT signs at each prime, so each query can be answered in polynomial time in $\log N$ by counting representatives in intervals. This is a factoring equivalence for the semiprime three-row family. The complete even-domain parity is always one on these inputs, regardless of which child contains the root.

For example, $39=3\cdot13$ and $879=3\cdot293$ agree modulo 840 and their respective prime factors agree modulo eight. Their hidden even roots are 14 and 292, in opposite classes modulo four. These local data therefore cannot decide the first child. The comparison addresses that particular shortcut. It supplies no lower bound for other ways of computing the parity.

For repeated prime powers the comparator matters. On $147=3\cdot7^2$, the first child $m=4,r=0$ gives $R=2$, $R_0=2$, with no public even root in that child. Its hidden parity is zero. Treating the main count alone as a parity would send the search to the wrong child.

## A silent regime on the original near-square branch

There is a concrete reason the earlier capped decoder does not directly extend to small least factors. Write $N=pq$ with $p<q$ odd primes. On an eligible original branch let $Q$ have discriminant $-\Delta_0$ and let $Q(\alpha)=N$ be the public representation. Each mixed CRT lattice shares its line modulo $q$ with either $\alpha$ or its conjugate $\bar\alpha$. At $p$ the vector vanishes after multiplication by $p$. Hence the lattice contains $v=p\alpha$ or $p\bar\alpha$, of norm $Np^2$.

A mixed CRT lattice has index $N$. If $w$ is independent of $v$ and belongs to the same lattice, then $|\det(v,w)|\ge N$. The Gram determinant gives
$$
Q(v)Q(w)\ge \frac{\Delta_0}{4}\det(v,w)^2,
\qquad Q(w)\ge \frac{N\Delta_0}{4p^2}.
$$
Consequently, for $4p^2B<\Delta_0$ every mixed-lattice point of norm at most $NB$ lies on the indicated public ray. Because $N$ is squarefree, $\alpha$ is primitive in integer coordinates; every integer point on its line is an integer multiple of it. The same determinant argument applied to $\alpha$ shows that all points in its public lattice up to this cutoff lie on its ray, and likewise for its conjugate. The entire divisible-point set in the disk therefore consists of public multiples. Every centered sublattice intersection of this disk inherits that equality, even if counted over the integers.

Two checks compute the full physical point sets independently. At $N=10{,}713=3\cdot3571$ the actual near-square branch has $c=206$, $\Delta_0=416$ and cap 11. All 12 nonzero signed points through the cap are public multiples. At $N=652{,}055=5\cdot130{,}411$ the branch has $c=1614$, $\Delta_0=3224$ and cap 32. All 20 points are public. This proves failure of mixed-point extraction on those disks and gives a sufficient condition for it. It does not exclude changing the branch or using other information about $N$.

## Reproducible checks

The [factorer](general_scope.py) uses only the Python standard library. The [checks](general_scope_check.py) use SymPy for independent prime factors and seeded prime generation. The factor-backed checker solves the normalized quadratic modulo $N$ on each of its two outer rows and counts the axis interval using the prime-power divisibility condition. Every scalar residue for the exhaustive split cohort and the historical sample was compared with direct enumeration using only the public form and cutoff. The separate complete-factorization cohort checks the final prime factorization. The [recorded results](general_scope_results.json) include the axis witness, the first-bit collision and the exact cohort sizes.

| Population | Validation | Result |
|---|---|---|
| Every integer from 1 to 2,000 | Complete recursive public-input factorer | 2,000 complete factorizations, 33,294 scalar calls |
| Every eligible odd modulus through 2,000 | Public enumeration checked against CRT per query | 676 splits, 12,372 matches, 168 nonsquarefree moduli, 149 with at least three different primes |
| Existing 194 semiprime inputs, using the new forms | Public enumeration checked against CRT per query | 194 splits, 5,260 matches |
| Two historical collision inputs | Same as above | $11{,}009=101\cdot109$, root 2,726. $28{,}829=127\cdot227$, root 7,492 |
| $30{,}021=3\cdot10{,}007$ | Public enumeration | 28 scalar queries and 72,468 visited positions |
| 72 seeded unbalanced inputs up to 322 bits | Factor-backed oracle, decoder validation only | 72 splits, 27,990 scalar calls, including 54 nonsquarefree inputs |

From the repository root run

    python experiments/ellipse_disk_mod8/general_scope.py 147
    python experiments/ellipse_disk_mod8/general_scope.py --complete 132300
    python experiments/ellipse_disk_mod8/general_scope_check.py

The row count is small, but the first dyadic child contains $\Theta(N)$ candidate positions on its outer rows. A polynomial-time solution still requires a factor-free method to compute their congruence parity in $\operatorname{poly}(\log N)$ time. The recursion and the axis comparator do not supply it.
