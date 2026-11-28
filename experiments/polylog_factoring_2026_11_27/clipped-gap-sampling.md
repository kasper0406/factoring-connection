# Sampling the clipped gap given the prime powers

With the complete factorization $N=\prod p^e$ supplied, the original one-sided comparison has a bounded-error polynomial-bit evaluator whenever

$$T=\prod_{p^e\parallel N}\big(1+2e(1-1/p)\big)=(\log N)^{O(1)}.$$

This includes every fixed number of distinct primes with arbitrary exponents. It also covers squarefree inputs with $O(\log\log N)$ prime factors. It works on each prescribed affine sign, window, upper digit and polarity, including highly atypical queries. [Sampler](clipped_gap_sampler.py) · [Checker](check_clipped_gap_sampler.py) · [Deterministic result on distinct semiprimes](clipped-gap-kloosterman.md)

Write $c$ for the known lower scale and $g(a)=\gcd(a^2-c^2,N)$. For one lower residue $a$, the normalized contrast summand is zero on nonunits and the corrected public roots. Otherwise it is the signed, windowed difference $(C_j(a)-C_0(a))/N$ defined in the original comparison. The numerator is $0$ or $\pm g(a)$. Sampling a residue uniformly would rarely hit a mixed stationary root. [Clipped block formula](square-input-isolation.md)

Use instead the probability $P(a)=g(a)/(NT)$ on **all** residues. To generate it locally modulo $p^e$, use one ticket class of weight $p$ and $e$ further classes of weight $2(p-1)$ each. In the first class draw a uniform residue modulo $p^e$. In class $i\in[1,e]$ choose a fair sign and a uniform lift of $\pm c\bmod p^i$. The classes have the required total mass because

$$\gcd(a^2-c^2,p^e)=1+\sum_{i=1}^e\varphi(p^i)\,1_{a=\pm c\bmod p^i}.$$

The two signs are disjoint at every positive depth. Local choices are independent and CRT gives the claimed global law. With this proposal, return $T$ times the signed windowed difference divided by $g(a)$, or zero at excluded residues. Its expectation is exactly $D(j)/N$ and its absolute value is at most $T$. Hoeffding's inequality yields additive error $\varepsilon$ except with probability $\delta$ from $O(T^2\varepsilon^{-2}\log(2/\delta))$ independent samples. Fixed $\varepsilon<1/4$ suffices at the midpoint of the original thresholds. Each local draw and contrast evaluation takes expected polynomial bit time. [Earlier full unit-mass identity](square-interval-count.md)

This proposal also identifies how much factor information the sampler uses. At a prime $p\mid N$,

$$s_p=\Pr[p\mid a-c]=\frac{1/p+e(1-1/p)}{1+2e(1-1/p)}\in[1/3,1/2].$$

For two distinct prime divisors $p,q$, independence gives $s_p(1-s_q)+s_q(1-s_p)\ge4/9$. On that event $\gcd(a-c,N)$ is a proper divisor. Consequently, a public sampler approximating this proposal to total variation distance $\eta<4/9$ splits with probability at least $4/9-\eta$ per draw on every odd input with two or more distinct primes. Constructing the biased sampler from the integer alone already supplies a factoring algorithm when that probability is inverse-polynomial.

Enumeration checked the exact proposal law on $N=15,25,27,45$. Forty-four contrasts matched literal inversion modulo $N^2$, with public corrections and repeated-prime cases. A separate 306-bit example checked the arithmetic of sample aggregation. Having many high prime powers may make $T$ large. The unrestricted worst-case gap can then express Subset Sum even when the prime factors are supplied. [Hard query construction](clipped-gap-subset-hardness.md)
