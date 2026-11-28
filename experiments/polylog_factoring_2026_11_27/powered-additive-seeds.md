# Additive tests after repeated powering

Write $F_e(x)=(1+x)^e-x^e-1$. A public input-dependent choice that escapes the earlier fixed menu is $x_j=a^{N^j}\bmod N$. Each step takes one modular exponentiation. There are two ways to use it.

* At the new seed, compute $\gcd(F_N(x_j),N)$. This requires forming $1+x_j$ before the final power.
* Alternatively, propagate the original pair $a,1+a$ side by side and compute $\gcd(F_{N^j}(a),N)$. This requires no addition between power steps.

The tests agree at their starting comparison $F_N(a)$. Thereafter they can give different factors. For $a=2$, $N=29{,}792{,}417$, the first powered seed is $25{,}590{,}566$, and its defect modulo $N$ is $13{,}500{,}524$. The gcd gives $4111$, with cofactor $7247$. The original seed, the five propagated comparisons and the gcd controls below all miss this input. The example was selected from the exhaustive scan. [Previous fixed-seed and uniform-seed bounds](additive-power-lift.md)

## A complete counterexample to steering by powering alone

For $N=31{,}439=149\cdot211$, put $H(x)=x(x+1)(x^2+x+1)$. Exhausting the residues verifies the following identities for **every** seed, including nonunits, with $y=x^N\bmod N$.
$$\gcd(F_N(x),N)=\gcd(H(x),N)=\gcd(H(y),N).$$
Applying the equality repeatedly shows that the first-digit gcd is constant along every powered orbit, at every depth. A seed that already gives a proper gcd exposes the same factor by evaluating the degree-four polynomial $H$ once. The four gcd classes $1,149,211,N$ contain $30{,}429$, $414$, $588$ and $8$ residues respectively.

At this input 2 is a primitive root at both primes, with multiplicative orders 148 and 210. Its full public orbit is
$$2,\ 12865,\ 12425,\ 17282,\ 30047,\ 240,\ 2.$$
The product of its six additive defects is $24709\bmod N$, a unit, so even arbitrary precision cannot produce a split from any of these six seeds by successive digit gcds. The complete period-six orbits from 3 and 5 likewise have unit products $5996$ and $8386$. Starting a fresh seed changes which gcd class one is in. Taking further powers of the same seed preserves the class.

The local explanation is that the complete zero set of $F_N$ modulo 149 consists of $0,-1$. Modulo 211 the only additional zeros are the roots $14,196$ of $x^2+x+1$. Powering by $N$ permutes each local field and preserves multiplicative order. The factorization explains the certificate. To verify its conclusion directly, loop over $0\le x<N$ and check the two displayed gcd equalities by modular powering. [Factor-free certificate and orbit checks](check_powered_additive_seeds.py)

## Measurements

For every unordered prime pair in $[R,2R)$, the scan uses $a=2$. The original seed is $x_0$. The other four fresh additive tests use $x_1$ through $x_4$. Controls check the gcds of $x_j$, $x_j+1$, $x_j^2+x_j+1$, $x_j^N-x_j$ and $(x_j+1)^N-(x_j+1)$ for $0\le j\le4$. They also check all pairs of the six orbit values $x_0,\ldots,x_5$ for equality under each of the six rational symmetries of the additive defect, with unit denominator. Each count below is a count of moduli. [Reproduce the complete scan](check_powered_additive_seeds.py)

| $R$ | Prime pairs | Split at $x_0$ | Split at $x_1$ | Extra at $x_1$ after controls and $x_0$ | Extra at any $x_1,\ldots,x_4$ after controls and $x_0$ |
|---:|---:|---:|---:|---:|---:|
| 64 | 78 | 3 | 3 | 0 | 0 |
| 256 | 903 | 10 | 12 | 2 | 7 |
| 1024 | 9,316 | 20 | 21 | 13 | 40 |
| 4096 | 107,416 | 54 | 67 | 34 | 147 |

All 147 in the last cell also escaped the five propagated comparisons. The propagated comparisons at depths 2 through 5 found 3,468 additional moduli beyond depth 1 at the same scale. For only 102 of those moduli was there a successful step at which $2^{N^j}$, reduced modulo the returned prime, lay outside the six transforms of 2. On a separate 1,000 draws each at $R=2^{20},2^{32},2^{48}$ with $p\in[R,2R)$ and $q\in[3R,4R)$, neither sequence produced a split. Those larger draws used separated prime ranges. The counts for both designs describe the tested inputs.

At $R=4096$ each prime was also labelled by the index $(r-1)/\mathrm{ord}_r(2)$, using its factorization data only for this retrospective count. A split is attributed to the prime returned by the proper gcd. Run the diagnostic with `--order-bins` to reproduce these columns.

| Index | Prime sides tested | Returned at $x_1$ | Returned at $x_1$ after all controls and $x_0$ |
|---|---:|---:|---:|
| 1 | 81,025 | 11 | 10 |
| 2 | 60,190 | 18 | 15 |
| 3–8 | 49,541 | 10 | 6 |
| Greater than 8 | 24,076 | 28 | 3 |

This locates 25 of the 34 surviving hits on the index-1 and index-2 sides. Those two bins also account for 141,215 of the 214,832 sides tested.

## Propagating both bases to polynomial depth

Draw distinct odd primes $r,s$ uniformly without replacement from $[R,2R]$, and set $N=rs$. Use a fixed positive rational seed $a$, with numerator and denominator of polynomial bit length in $\log R$. Fix $A\ge0$. Across all integers $1\le j\le(\log R)^A$, the probability that **any** propagated test $F_{N^j}(a)$ splits the modulus is at most
$$R^{-1/(3(A+1))+o(1)}.$$
A menu of polynomially many such seeds, fixed at the input size scale before drawing the primes, gives the same exponent. The bound allows the depth to depend on the observed input. It also bounds extracting unequal prime valuations at higher precision, since that first requires a hit at one prime.

To see this, fix $r$ coprime to the numerator and denominator of both $a$ and $1+a$. Put $L=\mathrm{lcm}(\mathrm{ord}_r a,\mathrm{ord}_r(1+a))$. The earlier unique-difference argument gives at most $1+\sqrt L$ exponent classes $S\bmod L$ with $(1+a)^e-a^e=1$. For $a\ne1$, its small-order estimate gives mean $L^{-1/2}\le R^{-1/3+o(1)}$ over regular primes, counting exceptions to regularity separately at cost $R^{-1+o(1)}$. For $a=1$, the hit exponents form the single class $1\bmod\mathrm{ord}_r 2$, with mean reciprocal order $R^{-1/2+o(1)}$. [Proofs of both order estimates](additive-power-lift.md#how-often-the-first-divisibility-test-is-reached)

At $r$, $N^j\equiv s^j\bmod L$. Also $s$ is coprime to $r-1$, since $r-1<2s$ and $r-1=s$ is impossible by parity. The map $s\mapsto s^j$ on units modulo $L$ has fibers of size at most $2j^{\omega(L)}$, where $\omega(m)$ counts distinct prime divisors. This follows by decomposing into prime powers. The unit group is cyclic at odd prime powers, and the 2-power contributes at most one extra factor 2. A set of $b$ allowable classes mod $L$ meets at most $b(1+R/L)$ integers in the interval. Dividing by the number of partner primes thus costs only $O(\log R)$ times its density $b/L$.

For any constant $c>0$, the fraction of primes $r$ in the interval for which $\omega(r-1)>c\log R/\log\log R$ is at most $R^{-c+o(1)}$. Indeed, take $k$ to be the ceiling of that threshold. Count multiples of products of $k$ distinct primes among all positive integers through $2R$. Their count, with multiplicity, is at most $2R(\sum_{\ell\le2R,\ \ell\text{ prime}}1/\ell)^k/k!$. The prime reciprocal sum is $O(\log\log R)$ and Stirling’s estimate gives the claim after division by $\asymp R/\log R$. On the remaining primes every fiber has size at most $R^{Ac+o(1)}$. The mean hit probability on one side is bounded by $R^{-c+o(1)}+R^{-1/3+Ac+o(1)}$. Choose $c=1/[3(A+1)]$, count both sides, and take a union bound over all depths and seeds.

**Polynomial exponents.** The same argument covers $F_{P(N)}(a)$ for a menu of polynomially many integer polynomials $P$ fixed before drawing the primes, where $P(N)>0$, $1\le\deg P\le(\log R)^A$ and the nonconstant coefficients have gcd 1. The resulting bound is $R^{-1/[3(A+3)]+o(1)}$. Indeed, $P(N)\equiv P(s)\bmod L$. For each prime $\ell\mid L$, the reduction of $P$ is nonconstant, so a prescribed value has at most $\deg P$ roots modulo $\ell$. Lifting all possibilities to modulus $L$ costs at most $L/\mathrm{rad}(L)\le(r-1)/\mathrm{rad}(r-1)$, where the radical is the product of the distinct primes. If this last quantity exceeds $R^{2c}$, then $r-1$ is divisible by $d^2$ for some $d>R^c$, by taking half of each prime exponent rounded down. At most $\sum_{d>R^c}\lfloor2R/d^2\rfloor=O(R^{1-c})$ positive integers through $2R$ qualify. Outside both exceptional sets the fiber cost is at most $R^{(A+2)c+o(1)}$. Taking $c=1/[3(A+3)]$ proves the bound. The earlier exponent study looked for powers covering a multiplicative group order. Here the target is the additive equality at the fixed pair of bases. [Earlier polynomial exponent test](bounded-quadratic-roots.md#fixed-polynomial-exponents-on-separated-balanced-primes)

## Seeds obtained by raising a fixed base to a power of the input

Fix a rational $a=A/B\notin\{0,1,-1\}$ in lowest terms and a fixed integer $j\ge1$. Draw distinct odd primes $r,s$ uniformly without replacement from $[R,2R]$, put $N=rs$, and compute $x=a^{N^j}$ modulo $N$. For each fixed $0<\epsilon<1$, the joint probability that $r\mid F_N(x)$ and $\operatorname{ord}_r(a)\le R^{1-\epsilon}$ is at most $R^{-\min(1/6,\epsilon)+o(1)}$ on primes where the order is defined. In particular the same rate, up to a factor two, bounds the probability of a split at any digit **with both local orders at most $R^{1-\epsilon}$**. A split at any depth requires a first-digit hit at one of the primes. The estimate controls a contribution to the joint probability over prime pairs. Conditioning on a potentially rare order range would require its probability as a denominator.

Write $M=r-1$, $T=\operatorname{ord}_r(a)$ and $H=\langle a\rangle$. We may exclude primes dividing $AB(A^2-B^2)$, a finite set for this fixed base. Thus $T>2$. The partner $s$ is coprime to $M$ since $M<2s$ and $M=s$ is impossible by parity. Modulo $r$ the seed is $a^{s^j}$. Set $c=s\bmod T$. This fixes $x=a^{c^j}$ of order $T$ and also fixes $x^s=x^c$. A hit now requires $(1+x)^s=1+x^c$. Both sides are units because $x$ and $x^c$ have order $T>2$. Writing $U=\operatorname{ord}_r(1+x)$, this equation permits at most one residue of $s$ modulo $U$. Together with $s=c\bmod T$, it permits at most $M/d$ residues modulo $M$, where $d=\operatorname{lcm}(T,U)$.

Let $\mu_d$ denote the subgroup of order $d$ in $\mathbb F_r^\times$. For each $d$ above, $x\in H$ and $1+x\in\mu_d$. A resulting $x$ can come from at most $\kappa_j(T)$ values of $c$, the kernel size of $c\mapsto c^j$ on $(\mathbb Z/T)^\times$. The decomposition of units modulo prime powers gives $\kappa_j(T)\le j^{\omega(T)+1}=R^{o(1)}$ for fixed $j$. Hence the number of eligible exponent residues divided by $M$ is at most
$$
\kappa_j(T)\sum_{T\mid d\mid M}
\frac{|H\cap(\mu_d-1)|}{d}.
$$

For $d<M/(M^{1/4}+1)$, the Garcia–Voloch shifted-subgroup estimate and $H\subseteq\mu_d$ bound the summand by $4d^{-1/3}\le4T^{-1/3}$. For $d\ge M/(M^{1/4}+1)\ge M^{3/4}/2$, expand the two subgroup indicators in multiplicative characters. Extending characters by zero at zero, the two principal characters contribute $(r-2)Td/M^2\le Td/r$. Each other Jacobi sum has absolute value at most $\sqrt r$, and the absolute coefficient weights sum to one. Thus the summand is at most $T/r+\sqrt r/d=T/r+O(r^{-1/4})$. [Garcia and Voloch, *Fermat curves over finite fields*](https://doi.org/10.1016/0022-314X%2888%2990007-8) · [Exact subgroup form and cutoff, equation (1)](https://arxiv.org/pdf/1102.1172#page=1)

Each residue modulo $M$ has $O(1)$ integer representatives in $[R,2R]$. Dividing their count by the number of available prime partners costs $O(\log R)$ relative to the fraction of exponent residues. Since $\tau(M)=R^{o(1)}$, the local hit probability at fixed $r$ is at most
$$
R^{o(1)}\bigl(T^{-1/3}+T/r+r^{-1/4}\bigr).
$$
For $T\le Y$, $r$ divides the nonzero integer $\prod_{k=1}^{\lfloor Y\rfloor}(A^k-B^k)$, whose log absolute size is $O_a(Y^2)$. Counting its prime divisors at least $R$ gives $\Pr(T\le Y)\ll_a\min(1,Y^2/R)$ at regular primes. Summing $T^{-1/3}$ over dyadic ranges below $\sqrt R$ and using $T^{-1/3}\le R^{-1/6}$ above gives $\mathbb E[T^{-1/3}]\ll_a R^{-1/6}$. On the order range of the claim, $T/r\le R^{-\epsilon}$. This proves the stated bound.

When $T=r^{1-o(1)}$, including when $a$ is a primitive root, the term $T/r$ can preclude a power saving. This calculation gives no negligible unconditional hit rate across all balanced prime pairs. For all $1\le j\le(\log R)^A$ simultaneously, the bound on $\omega(r-1)$ in the previous section gives a joint probability at most $R^{-b/(A+1)+o(1)}$ on the same restricted order range, where $b=\min(1/6,\epsilon)$. Indeed, outside an exceptional fraction $R^{-c+o(1)}$, all kernel sizes are at most $R^{Ac+o(1)}$. Apply the average estimate to the parenthesized expression above, take $c=b/(A+1)$ and sum over the polynomially many depths. This extension still leaves the near-full-order case open.

## Normalizing a failed pair

One proposed update after a miss is $U(a)=a^N/((1+a)^N-a^N)$. Write $g(t)=-t/(1+t)$. Since $N$ is odd, $g$ is an involution and $U(a)=g(g(a)^N)$. Moreover $F_N(g(t))=-F_N(t)/(1+t)^N$. Wherever the denominators are units, all successive normalized tests therefore have exactly the same divisibilities as the powered seeds starting at $g(a)$.

The denominator itself sometimes exposes a prime. At depth $j$ the homogeneous denominator is $D_j=(1+a)^{N^j}-a^{N^j}$. For a uniform $a\bmod N$, $N=pq$ with distinct odd primes $p<q$, the probability of a proper denominator gcd is exactly
$$\frac{p^{\min(j,v_p(q-1))}-1}{q}.$$
Indeed, at a prime $r\mid N$, a pole requires $a\ne0$ and $((1+a)/a)^{N^j}=1$. As $a$ ranges over nonzero residues the ratio ranges over the field minus 1. This gives $\gcd(N^j,r-1)-1$ solutions, zero at $p$. The pole sets grow with $j$. For example $N=1379$ and $a=17$ give a unit defect and a unit denominator. At $U(a)=325$ the defect is again a unit, while the denominator gives $\gcd(788,1379)=197$. Here $197=4\cdot 7^2+1$.

If this denominator scan succeeds with probability at least $\delta>0$, then $q=kp^h+1$ for $1\le h\le\log_2 N$ and $k<1/\delta$. For each such pair of integers $k,h$, binary search for an integer solution of $kt^{h+1}+t=N$ and verify it divides $N$. The left side is strictly increasing for positive $t$. This gives a direct deterministic split in time polynomial in $\log N$ and $1/\delta$ on every such input. For $p<q<2p$ the denominator success probability is zero.
