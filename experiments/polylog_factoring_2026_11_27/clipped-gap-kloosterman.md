# Deterministic equivalence of the clipped gap and factoring on semiprimes

For $N=pq$ with distinct odd primes, factoring and the repository's prescribed one-sided clipped-count comparison reduce to each other in deterministic polynomial bit time. For the forward reduction, one can fix the lower scale $c=1$, use constant positive weights and find a mixed root by binary search with $O(\log N)$ threshold queries. For the converse, the supplied factors suffice to answer **each** promised comparison with arbitrary allowed affine sign, lower window, scale, upper digit and polarity.

The analytic step is uniform cancellation of the nonstationary lower digits. In particular, on balanced large inputs the actual clipped contrast differs from its stationary contribution by $O(N^{3/4}\log^4 N)$ on every query. Computing the mixed stationary contribution from unfactored input is still open. The older reduction handles arbitrary odd inputs by random isolation and a mean-square estimate. The older deterministic half-window bisection uses an exact parity evaluator. Here the one-sided gap itself gives the interval membership needed for bisection. [Original signed count and thresholds](single-quadratic-count.md) · [Earlier randomized gap reduction](count-comparison.md) · [Known half-window root parity](single-quadratic-count.md#deterministic-starting-interval-with-scale-one)

## Exact contribution of a proper gcd

Let $n$ be the bit length of $N$, $H=(N-1)/2$, $I$ the cyclic lower-digit window and

$$
A(a)=1_I(a)\big(2[(\alpha a+\beta\bmod N)\le H]-1\big).
$$

The change of upper scale from $c$ to $c+Nj$ produces the raw signed count difference at $j$ relative to zero. Subtract the two known contributions from $a=\pm c\bmod N$ to obtain $D(j)$. At polarity $\sigma$ the comparison must say zero if $\sigma D(j)\le N/4$ and one if $\sigma D(j)\ge3N/4$. Assume $\gcd(c,N)=1$ as in the reduction.

For a unit lower digit $a=x\bmod N$, set $b=(a+c^2/a\bmod N)$ in $[0,N)$ and $g=\gcd(a^2-c^2,N)$. As $x$ runs over its $N$ lifts, the phase $x+(c+Nj)^2/x\bmod N^2$ has lower digit $b$. Its upper digits occupy one residue class modulo $g$, each with multiplicity $g$, by expansion to first order in $N$. The clipped multiplicity changes by $g$ times the change in the indicator of the first $(g+1)/2$ upper residues modulo $g$. Nonunits have no solutions and $g=1$ contributes no change.

Consider $g=d\in\{p,q\}$ and put $m=N/d$. The stratum splits by $a\equiv s c\bmod d$ with $s=\pm1$. Apart from $a\bmod m=0,\pm c$, each such digit has $g=d$. Let $b_0$ be the least residue of $2sc$ modulo $d$ and write $b=b_0+dk$, $0\le k<m$. Since $x\equiv s(c+Nj)\bmod d$, the identity

$$
x+(c+Nj)^2/x-2s(c+Nj)=(x-s(c+Nj))^2/x
$$

fixes the phase modulo $d^2$. Consequently its upper digit $z$ satisfies

$$
z\equiv 2sj+C-m^{-1}k\pmod d,\qquad C=m^{-1}(2sc-b_0)/d\pmod d,
$$

where in $C$ the integer division by $d$ precedes multiplication by $m^{-1}\bmod d$. Write $B(k)$ for the difference of the two $(d+1)/2$-interval indicators at these residues, using $j$ and $0$. With $v=a\bmod m$, we have

$$
k(v)=d^{-1}(v+c^2/v-b_0)\pmod m,
$$

using the least representative modulo $m$. Thus the entire sign stratum contributes $d\sum_{v\ne0,\pm c}A(a(v))B(k(v))$, where $a(v)$ is given by CRT.

## A uniform bound for every window and upper digit

Use normalized Fourier coefficients on each finite cyclic group, and put $L=n+2$. An interval indicator on $\mathbb Z/t$ has Fourier coefficient sum in absolute value at most $L$ for $t\le N$, by the geometric-series bound and the harmonic sum. The same sum for the length-$m$ word $\exp(2\pi i\theta k)$ is at most $4L$ for every real $\theta$. Indeed the coefficient at $h$ has magnitude at most $\min(1,1/(2\|m\theta-h\|_m))$, where $\|\cdot\|_m$ is circular distance to zero. Sum by unit distance shells.

On the progression $a=a_*+di$, $0\le i<m$, the lower window is a cyclic interval in $i$. The affine sign asks whether $\alpha i+\gamma\bmod m$ lies in a prefix interval, since the residue of $\alpha a+\beta$ modulo $d$ is fixed. Changing coordinates from $i$ to $v$ is an affine permutation. The Fourier coefficient sum for $A(a(v))$ is at most $3L^2$, by multiplying the window and sign. This also covers $\alpha\equiv0\bmod m$. On $\mathbb Z/d$ the corresponding norm of $B$ is at most $2L$. Restricting its Fourier expansion to integer representatives $0\le k<m$ and using the $4L$ bound for each fractional frequency gives norm at most $8L^2$ on $\mathbb Z/m$.

For each nonzero Fourier frequency $h$ of $B$ and every frequency $\ell$ of $A$, the sum over $v\ne0$ of the associated product is a phase times

$$
\sum_{v\in\mathbb F_m^\times} e_m\big((\ell+h/d)v+(hc^2/d)/v\big),
\qquad e_m(t)=\exp(2\pi i t/m).
$$

Its reciprocal coefficient is nonzero. The prime-field Kloosterman bound gives magnitude at most $2\sqrt m$, also when the linear coefficient vanishes. The nonzero frequencies of $B$ therefore contribute at most $48L^4\sqrt m$ before multiplying by $d$. [Weil, *On Some Exponential Sums*](https://doi.org/10.1073/pnas.34.5.204)

The zero frequency of $B$ has magnitude at most $8/m$. Full periods of length $d$ in its defining sum cancel. For the remaining $t=m\bmod d$ terms, let $v_0=m^{-1}\bmod d$. Then $v_0/d=h_0/t+1/(dt)$ for an integer $h_0$ coprime to $t$. The progression of $t$ points is within $1/d<1/t$ of a permuted complete $t$-grid. Each half-interval count differs by at most four from $t$ times its length, allowing for the two endpoints and grid rounding. Both intervals have equal length. Finally we delete $v=\pm c$, costing at most two. The absolute value of each sign stratum is bounded by

$$
d\big(48L^4\sqrt m+10\big).
$$

Summing the four bounds proves the stated balanced estimate. The bound used all upper digits and lower signs simultaneously. It applies to an individual fixed query and costs no rerandomization of that query.

## Evaluate the gap given the factors

Set $T=2^{22}L^8$. In each stratum with $m<T$, enumerate its $m$ possible residues using the exact formula above. When $m\ge T$, omit the stratum. Its error is below $N/32$ since $\sqrt m\ge 2048L^4$. Compute the two mixed $g=N$ stationary changes exactly by CRT and one modular inverse per root. The known $a=\pm c$ terms are already removed from $D$. Thus the integer estimate $\widetilde D$ has $|\widetilde D-D|<N/8$. Return one when $\sigma\widetilde D\ge N/2$. This answers the prescribed promise deterministically on every distinct odd semiprime, including unbalanced products. The procedure uses $O(T)$ arithmetic operations on polynomial-bit integers. It concerns the gap between thresholds and supplies no algorithm for the exact packed count with the factors given.

## Find the factor using the gap alone

Trial-divide through $T$ first, stopping early if a divisor is found. On the remaining semiprimes both cofactors of proper strata exceed $T$. Uniformly on every allowed window, the normalized corrected contrast is within $1/8$ of the mixed stationary contribution divided by $N$.

Set $c=1$. Use positive weights by taking $\alpha=\beta=0$, or $\alpha=1,\beta=0$ since every chosen interval is contained in $[2,H]$. The interval $I=[2,H]$ contains exactly one mixed root $a$ of $a^2=1\bmod N$, since the two mixed roots form an opposite pair. It contains neither known root. Its stationary contrast is $N$ times

$$
1_{[0,H]}(z+sj\bmod N)-1_{[0,H]}(z),\qquad s=2/a\bmod N ,
$$

for an unknown upper digit $z$. The unit $s$ guarantees a sign change at one of the public shifts $j=\pm 2^i$, $0\le i\le n$. Indeed, let $m_0$ be the smaller distance of $s$ to zero on the residue circle. Double it until the first $\delta=2^i m_0\ge N/4$. Then $\delta\le H$. In either of the two consecutive runs of the original bit, length at most $H+1$, the three residues $z,z-\delta,z+\delta$ cannot all lie in that run. In coordinates within the run both differences would have to be the ordinary differences $-\delta,+\delta$, since $N-\delta>H$. This would require $2\delta\le H$.

Try both polarities at the public shifts until a one is returned. No promise ambiguity occurs on these queries, since the stationary contribution is always $0$ or $\pm N$. Fix this shift and polarity. On every subinterval of $I$, the response is now precisely the predicate that it contains $a$. Bisect to recover $a$, then return $\gcd(a-1,N)$. If only dyadic widths are admitted, start with the two width-$2^{\lfloor\log_2|I|\rfloor}$ intervals at the ends of $I$. They cover $I$ and are contained in it. Try both in the shift search, then bisect the one that returned one. The whole reduction needs $O(n)$ queries and deterministic polynomial bit work. The preliminary cutoff is an effective asymptotic bound with a large constant.

The checker matched 72 original modular-hyperbola contrasts, 432 proper gcd strata and 4,500 grid sums. Forty narrow-window checks exercised actual omission of two or four strata at 131 and 234 bits. Its forward procedure receives only $N$ and a comparison callable. With the factor-backed comparator standing in for that oracle, it split four fixtures in 15, 18, 24 and 249 queries. The three small bisection fixtures deliberately bypassed preliminary trial division to exercise the search. [Reproduce the checks](check_clipped_gap_kloosterman.py)
