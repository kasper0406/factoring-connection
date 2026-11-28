# The next cyclotomic digit at a canonical base

There are two different inputs to distinguish even when their low digits agree. A public power map generates bases at which $\Phi_N$ is explicitly evaluable modulo $N^j$ for every desired precision $j$. Reducing its output to the ordinary representative $1\le x<N$ can change the next digit. For a product of distinct odd primes, predicting that digit on even an inverse-polynomial fraction of uniform canonical units suffices for randomized polynomial-time factoring. On balanced inputs the decoder usually needs just one correct guess. It verifies a candidate factorization directly.

The earlier result used the derivative at a specified base or a chosen lift such as $1+N$. Here the evaluator receives only an ordinary integer in $[1,N)$, with a uniform marginal distribution. [Earlier value and derivative result](bounded-quadratic-roots.md#value-and-derivative-modulo-the-index) · [Earlier audit of Fermat quotients](../../notes/user_ideas/2026-10-18_paradigms_4_6.md)

## A public family with exact values to any precision

Given $N>1$, an integer $a$ and $j\ge1$, take $K=\max(j,\operatorname{bitlength}(N))$ and $b=a^{N^K}\bmod N^j$. If $N$ has at least two distinct prime divisors, then

$$\Phi_N(b)\equiv1\pmod{N^j}.$$

If $N=p^e$, the answer is $p\bmod N^j$ when $a\equiv1\bmod p$, and 1 otherwise. Recognition of prime powers uses integer root and primality tests. Generating the base and the residue takes time polynomial in the bit lengths of $a$, $N$ and the precision $j\log N$.

For the proof, fix $p^e\Vert N$ and write $m=N/p^e$. If $p\mid a$, then $b=0\bmod p^{ej}$ and $\Phi_N(0)=1$. Otherwise $b^{p-1}=1\bmod p^{ej}$ by Euler's theorem, since $p^{ej-1}\mid N^K$. Also the order of $b$ modulo $p$ is coprime to $N$. Raising $a$ to $N^K$ removes every prime power in its order supported on primes dividing $N$, since the original order divides $p-1<N$. If $m>1$, the denominator in the identity

$$\Phi_N(b)\Phi_m\!\left(b^{p^{e-1}}\right)=\Phi_m\!\left(b^{p^e}\right)$$

is a unit at $p$. In characteristic $p$, the roots of $\Phi_m$ have exact order $m$ since $p\nmid m$. Its two arguments agree modulo $p^{ej}$ because $b^p=b$ there. Division and CRT prove the claim. For $m=1$, the prime-power cyclotomic polynomial is a sum of $p$ powers of $b^{p^{e-1}}$. This sum is 1 if $b-1$ is a unit at $p$, and $p$ if $b=1\bmod p^{ej}$. On unit inputs $b\bmod p=a\bmod p$ in this prime-power case.

## What reducing the base discards

Write $N=pq$ with distinct odd primes $p<q$, and put $s=p+q$. For a unit residue $x$, let $\omega(x)\bmod N^2$ be its lift whose components satisfy $t^{p-1}=1\bmod p^2$ and $t^{q-1}=1\bmod q^2$. Such a lift exists uniquely and is multiplicative. Locally the lift of $b\bmod p$ is $b^p\bmod p^2$, independent of the representative. Knowledge of the primes is only used to describe this map.

Call the case $\gcd(x^N-1,N)=1$ generic. Here $\Phi_N(x)=1\bmod N$ and $\Phi_N(\omega(x))=1\bmod N^2$. The latter follows locally from $\Phi_N(t)\Phi_q(t)=\Phi_q(t^p)$ at $p$, since $t^p=t\bmod p^2$ and $t^q-1$ is a unit. Exchange $p,q$ for the other component. Write $E(x)=(\Phi_N(x)-1)/N\bmod N$. Logarithmic differentiation of the semiprime formula gives

$$x\Phi'_N(x)=A(x)-sB(x)\pmod N,\qquad
A(x)=\frac{x}{x-1},\quad B(x)=\frac{x^N}{x^N-1}.$$

Indeed the terms with coefficients $N,p$ vanish modulo $p$, and $x^N=x^q\bmod p$. Exchange the primes and use CRT. Taylor expansion about $\omega(x)$ now gives

$$E(x)=\frac{x/\omega(x)-1}{N}\,[A(x)-sB(x)]\pmod N.$$

The division by $N$ takes a residue known to be $1\bmod N$, subtracts 1 and reads its next digit. Every such quotient in this note is well defined modulo $N$.

## One canonical query on the coprime-order case

Draw a uniform unit $z\bmod N$ and compute $w=z^N\bmod N^2$, $x=w\bmod N$ and $h=\lfloor w/N\rfloor/x\bmod N$. Then $w=\omega(x)$, because raising to $N$ kills both local principal-unit components. All of $w,x,h$ are public, and

$$E(x)=-h[A(x)-sB(x)]\pmod N.$$

A predicted value of $E(x)$ gives $s$ by one linear equation when $h$ is a unit. Check that its canonical representative has square discriminant $s^2-4N$ and produces two proper factors with product $N$. Any proper gcd encountered also gives a factor.

For example, $N=77$ and $z=3$ give $w=3155$ and $x=75$. The exact residues at the two bases are $\Phi_{77}(3155)=1$ and $\Phi_{77}(75)=771$ modulo $77^2$. Thus $E(75)=10$. With $h=57$, $A=52$ and $B=5$, the decoder obtains $(10+57\cdot52)/(57\cdot5)=18\bmod77$. The discriminant is $18^2-4\cdot77=16$, giving 7 and 11.

If $p\nmid q-1$, powering by $N$ permutes the units modulo $N$, so $x$ is uniform. This includes all $p<q\le2p$. The loss from $h=0\bmod N$ is at most $1/(p-1)$ by the distinct-carry argument below. The generic checks also give a split unless $x=1$ whenever $\gcd(x^N-1,N)>1$. In the case $x^N=1\bmod N$, the order of $x\bmod p$ divides $\gcd(N,p-1)=1$, so $\gcd(x-1,N)$ detects the split except at $x=1$. Thus a candidate routine correct with probability $\delta$ on uniform canonical units gives a verified split with probability at least $\delta-1/(p-1)-1/\varphi(N)$ per trial in this case.

## Three canonical queries for every distinct odd semiprime

An evaluator may be wrong on the entire image of $z\mapsto z^N$ when $p\mid q-1$. Sample a uniform unit $u\bmod N^2$, draw three independent $z_i$ as above and use $w_i=u z_i^N\bmod N^2$. Put $x_i=w_i\bmod N$ and $h_i=\lfloor w_i/N\rfloor/x_i\bmod N$. Each $x_i$ is marginally uniform on the canonical units. All three share the unknown coordinate $\theta=(u/\omega(u\bmod N)-1)/N\bmod N$. Multiplicativity gives $(x_i/\omega(x_i)-1)/N=\theta-h_i\bmod N$.

For a generic query set $r_i=A(x_i)/B(x_i)$ and use the returned digit to compute $y_i=E(x_i)/B(x_i)+h_i r_i$. Correct answers obey a common linear equation

$$ y_i=\theta r_i-\theta s+h_i s\pmod N.$$

Check all differences $r_i-r_j$ for proper gcds. If they all vanish modulo $N$, subtract rows to recover $s$ from a unit difference of the $h_i$, checking any proper gcd. Otherwise at least one $r$ difference is a unit. Subtract the first row from the other two and use the two-by-two determinant to eliminate $\theta$. A unit determinant gives $s$, and a proper gcd of the determinant gives a factor. Verify every result. With the factors supplied, generic values are efficiently available from the semiprime geometric sums. Even on exceptional inputs one can evaluate the polynomial quotient at sufficient precision: the identity $\Phi_{pq}(x)=S_q(x^p)/S_q(x)$ with $S_q(x)=1+\cdots+x^{q-1}$ always has unit denominator at $p$ for $p<q$, and at $q$ except when $x=1\bmod q$. At that last component, $\Phi_{pq}(x)=S_p(x^q)/S_p(x)$ has unit denominator at $q$. Both sums use binary powering and CRT.

There are two probability issues.

1. **Correctness on chosen cosets.** Given $u$, the three bases are independent on one coset of the power image. If the conditional success rate of the evaluator is $d_u$, then $\mathbb E[d_u]=\delta$ by the uniform marginals, and $\mathbb E[d_u^3]\ge\delta^3$. This allows errors to concentrate on any input set. A randomized evaluator uses fresh randomness on each call.
2. **Distinct carries.** Condition on $u$ and all three residues $x_i\bmod q$. Each $w_i\bmod q^2$ is then fixed, since $z_i^N$ is the unique lift of its residue at $q$ with order dividing $q-1$. Meanwhile $z_i^N\bmod p$ is uniform on all nonzero residues since $\gcd(N,p-1)=1$. For a fixed $x_i\bmod q$ the canonical candidates have the form $b+jq$, $0\le j<p$, omitting the one divisible by $p$. At $q$,

   $$ h_i\equiv (w_i-x_i)/(q p x_i)\pmod q,$$

   where the numerator is first divided by $q$ using its value modulo $q^2$. Changing $j$ gives $p-1$ distinct equally likely values because $p<q$ and $p x_i$ is a unit modulo $q$. The three choices are conditionally independent. On residues generic at $q$, the $r_i\bmod q$ are already fixed. If they are not all equal, the determinant has a nonzero coefficient for at least one $h_i\bmod q$ and vanishes there with probability at most $1/(p-1)$. If all are equal at $q$, a difference nonzero at $p$ splits $N$, or the slope in $h$ suffices unless all $h_i$ agree at $q$. The same probability bound suffices. This estimate is made **before** conditioning on honest oracle answers.

In total the chance of a verified split is at least $\delta^3-1/(p-1)-3/\varphi(N)$. Trial division to a constant times $\delta^{-3}$ removes the small-prime case, leaving an expected $O(\delta^{-3})$ batches. The resulting bit cost is polynomial in $\log N$, $1/\delta$ and the cost of a call. At no point does the decoder supply the high lift, its carry or its private randomizer to the evaluator. The result assumes inverse-polynomial success for the modulus being factored. It supplies no predictor with that accuracy.

**A single prescribed base.** At $x=2$ on the generic branch, let $F\bmod r$ be the Fermat quotient $(2^{r-1}-1)/r\bmod r$ for each prime $r\mid N$. The lift at $r$ is $2^r\bmod r^2$, so the normalized difference $(2/\omega(2)-1)/N$ equals $-F/s\bmod N$ by CRT. Hence $E(2)=F(B(2)-2/s)\bmod N$. At $N=2099\cdot3511$ this digit vanishes. The Fermat quotient vanishes at 3511, while $s B(2)-2$ vanishes at 2099. The probabilistic decoder above uses a newly sampled canonical base with a known carry; factoring from one fixed value $\Phi_N(2)\bmod N^2$ remains open here.

## Which public formulas miss the digit?

Consider a proposed digit whose value modulo $p$ depends only on $x\bmod p$, and whose value modulo $q$ depends only on $x\bmod q$. This includes arithmetic modulo $N$ using fixed additions, multiplications, powers and unit divisions with input-independent coefficients, even when the coefficients depend on $N$. On uniform generic canonical units, any such formula agrees with $E(x)$ with probability $O(1/p)$. A menu of $J$ such formulas, even with an input-dependent choice of index, gives at most $O(J/p)$. A formula that reads an integer carry can fall outside this class.

To see this, fix $r\in\{p,q\}$, let $\ell=N/r$, and fix the residue $b\bmod r$ of a generic input. The canonical lifts are $b+jr$ with $0\le j<\ell$, with nonunits excluded. Taylor expansion modulo $r^2$ makes $E(b+jr)\bmod r$ affine in $j$ of slope $\Phi'_N(b)/\ell\bmod r$. If the slope is nonzero, a fixed local answer can agree on at most $\lceil\ell/r\rceil$ of the lifts. Summing the contribution from nonzero slopes at either prime gives probability $O(1/p)$ among all units.

It remains to bound the residues where both slopes vanish. At $r$ a zero slope requires

$$ (\ell-1)b^\ell-\ell b^{\ell-1}+1=0\pmod r. $$

When $\ell\not\equiv1\bmod r$, each ratio $v\ne1$ occurs for at most one ordered pair of nonzero roots $b,vb$. Subtracting their equations and dividing by $b^{\ell-1}$ gives $(\ell-1)(v^\ell-1)b=\ell(v^{\ell-1}-1)$. Both coefficients could vanish only for $v=1$. A set of $k$ roots therefore satisfies $k(k-1)\le r-2$. Applying this at both primes and using CRT gives $O(1/\sqrt N)$ for the simultaneous-zero event. If $q\equiv1\bmod p$, the equation at $p$ instead has at most $\gcd(q-1,p-1)\le\min(p-1,(q-1)/p)$ nonzero roots. The ratio bound still applies at $q$. The combined fraction is $O(\min(1,q/p^2)/\sqrt q)=O(1/p)$. Finally, the generic domain has at least a quarter of all units. At $p$ only 1 is excluded, and at $q$ the excluded subgroup has size $\gcd(N,q-1)\le(q-1)/2$. Conditioning on this domain preserves the bound.

## Reproducible controls

* The factor-free entry point `evaluate_power_lift(n, a, precision)` constructs the base and returns its exact residue. Direct modular Horner evaluation of explicit polynomials verified all 3,129 full-depth cases at $2\le N\le150$, precisions 1 through 3 and seven signed bases. Another 3,129 checks used only the precision as the power depth, which left 20 exceptional local residues on indices with several distinct prime factors. The local Taylor formula passed 18,112 checks on 9,056 generic inputs across 24 balanced semiprimes. [Public evaluator and exact checks](check_teichmuller_cyclotomic.py)
* On $N=1009\cdot1013$, 128 single-query trials with the factor-backed exact oracle produced 126 verified factor sums and two splits at preliminary gcd checks. With persistent errors on about a quarter of queried inputs, the same draws produced 91 verified sums, two gcd splits and 35 rejected candidates. The canonical decoder receives no factors. The reference oracle uses them to supply test digits. [Canonical-base decoder and controls](check_canonical_cyclotomic_lift.py)
* On $N=11\cdot67$, an adversarial test oracle answered incorrectly on the entire image of $z\mapsto z^N$. The one-query control had zero verified sums. Over 128 random-coset batches the three-query decoder recovered 13 verified sums, split the input by gcd in 102 more, rejected 12 and encountered one singular system. One successful batch queried 542, 329 and 381 and recovered $p+q=78$. These finite rates illustrate the role of the coset choice. The probability guarantee above covers arbitrary locations of the correct digits.
