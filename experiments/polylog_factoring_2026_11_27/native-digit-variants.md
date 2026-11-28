# Can a small multiplier give the digit polynomial useful factors?

Let $p,q$ be two distinct primes sampled uniformly from the $b$-bit primes, let $N=pq$, and put $R=2^{b-1}$. Try any integer multiplier $1\le |a|\le H$, write $aN$ in base $2$ or $-2$, factor its digit polynomial over $\mathbb Z$, and evaluate a chosen polynomial divisor at the same base. For positive base take the digits of the absolute value and keep the sign separately. Put $\gamma=\log_2(4/3)/4=0.103759\ldots$. For a fixed $0\le\delta<\gamma$ and $H\le R^{2\delta}$, the probability **any** of those evaluations gives a proper gcd with $N$ is at most

$$
N^{-\gamma+\delta+o(1)} .
$$

Here powers of $N$ describe the common input scale $\log_2N=2b+O(1)$. The bound counts the existence of a successful multiplier, so the multiplier can depend on the whole input. In particular all multipliers of size polynomial in $b$ together still have the exponent $-\gamma+o(1)$. Multiplication by a power of two only inserts zero digits at the low end. The corresponding polynomial factor $X$ has gcd one with odd $N$. In negative base its removal can flip the sign, which is already allowed.

The unmodified method appears explicitly in Breuillard and Varjú. Choosing a multiple before selecting a polynomial also occurs in the special number field sieve, where the evaluation point and coefficients can be chosen differently. The bound here concerns the canonical zero-one digits at the two indicated evaluation points. [Breuillard and Varjú, §1.1](https://arxiv.org/pdf/1810.13360#page=2) · [Elkenbracht-Huizing et al., §§1–2](https://homes.cerias.purdue.edu/~ssw/cnta5.pdf) · [Unmodified binary case and prior art](native-digit-polynomial-bound.md)

## Which factorizations can work

For a monic binary polynomial $P$ of degree $d$ and constant term one, write $P^*(X)=X^dP(1/X)$. Choose its interior coefficients uniformly and independently. Kogler's Theorem 1.3 and Lemma 3.1 imply that, except on a fraction at most $8(3/4)^{\lfloor(d-1)/4\rfloor}$ of words, every decomposition $P=FG$ over the integers has $F=\pm F^*$ or $G=\pm G^*$. Indeed, reversing just one nonreciprocal factor of a pair of nonreciprocal factors produces a distinct binary polynomial after normalizing the sign, with the same reciprocal product. [Kogler, Theorem 1.3 and Lemma 3.1](https://arxiv.org/html/2609.26771v2)

Suppose $P(\beta)$ is a multiple of $N$, where $\beta=2$ or $-2$, and $\gcd(F(\beta),N)=p$. Then $p\mid F(\beta)$ and $q\mid G(\beta)$. If $F$ is reciprocal up to sign, $p\mid P^*(\beta)$. If $G$ is reciprocal up to sign, $q\mid P^*(\beta)$. Hence a successful split outside the exception requires a $b$-bit prime dividing both $P(\beta)$ and $P^*(\beta)$.

## Counting words with a common large prime

Use a central interval of $w\le b-1$ digit positions invariant under reversal $i\mapsto d-i$, choosing the maximum such width with the correct parity. For the degrees used below $d\ge2b-O(1)$, and this interval is interior. Fix all digits outside it, its center digit if present, and the sum of the bits in each of its $\lfloor w/2\rfloor$ mirror pairs. This determines
$$
L=P(\beta)+P^*(\beta).
$$
It is nonzero. At base two both values are positive. At negative two both have sign $(-1)^d$, because the leading term dominates the sum of the absolute values of all lower powers. Also $|L|\le 2^{d+2}$. For $d=O(b)$ there are therefore at most $O(1)$ candidate $b$-bit prime divisors of $L$.

Let $J$ count pairs whose sum is one. Conditional on the fixed data, all $2^J$ orientations are equally likely. For each candidate odd prime $r$, at most one orientation makes $r\mid P(\beta)$. The difference between two evaluations is $\beta^s$ times a nonzero integer of absolute value at most $2^w-1<R$, where $s$ is the low end of the interval. Nonvanishing follows from uniqueness of the base expansion, and $\beta$ is a unit modulo $r$. Before conditioning, $J$ is binomial with parameters $\lfloor w/2\rfloor,1/2$. Averaging gives
$$
\Pr(\text{a common }b\text{-bit prime})\le O(1)\,\mathbb E[2^{-J}]
=O((3/4)^{b/2}).
$$

## Counting all multiples once

After stripping powers of two it suffices to consider odd signed integers $z$ with $R^2\le |z|\le4HR^2$. Each contributes at most one distinct balanced semiprime $N$ dividing $z$, since $4HR^2<R^3$ for all sufficiently large $b$ in the stated range. For $\beta=-2$, a monic word has $2^d/3\le |P(-2)|\le 4\cdot2^d/3$, by summing the geometric series for positions of opposite and matching parity. Thus for either radix all degrees satisfy $2b-O(1)\le d\le 2b+\log_2H+O(1)$.

For each degree, multiply the fraction from the twin and common-prime bounds by the $2^{d-1}$ possible words. Summing the geometric progression through the largest degree counts $O(HR^2(3/4)^{b/2})$ exceptional words. Signs and bases cost a constant factor. Divide by the $\asymp R^2/b^2$ distinct balanced prime products to get $O(b^2H(3/4)^{b/2})$, as claimed. [Finite checks for both bases](check_native_digit_variants.py)

## A larger multiplier range from the reciprocal factor degree

Kogler's Lemma 3.5 bounds the probability that a uniform binary word has a monic reciprocal factor of degree at least $L$ by $8\cdot2^{-L/2}$. For either base, a factor $F$ of degree $l$ obeys $|F(\beta)|\le3^l\mathcal M(F)$, where $\mathcal M$ is Mahler measure. This follows from $|\beta-z|\le3\max(1,|z|)$ for each root. By multiplicativity and monicity, $\mathcal M(F)\le\mathcal M(P)\le\sqrt{d+1}$, the last bound by Jensen's inequality and Parseval on the unit circle. A factor whose value is divisible by a $b$-bit prime therefore has degree at least $(b-O(\log b))/\log_2 3$ for $d=O(b)$. The evaluated factors are nonzero since $P(\beta)=aN\ne0$. Both monic factors have constant term one, so reversal up to sign here means exact equality. [Kogler, Lemmas 3.1 and 3.5](https://arxiv.org/html/2609.26771v2)

For $H\le R^{2\delta}$ the reciprocal-factor contribution after counting all words is at most $N^{\delta-1/(4\log_2 3)+o(1)}$. For the twin contribution, keep Kogler's degree-dependent bound inside the sum. Summing $O(2^{(1-\gamma)d})$ words per degree up to $2b+\log_2H+O(1)$ and dividing by the number of balanced semiprimes gives $N^{-\gamma+\delta(1-\gamma)+o(1)}$. The sum of the two bounds decays as a power of $N$ for every fixed $\delta<\gamma/(1-\gamma)=0.1157\ldots$. At $H=1$ or $H=b^{O(1)}$ the larger contribution is still $N^{-\gamma+o(1)}$.

The finite check runs both bases, reconstructs 12,290 signed or unsigned canonical expansions, and checks 5,616 families of mirror orientations. At $N=551=19\cdot29$, using $a=5$ gives polynomial factors with values $5,19,29$ at two, while the reversed word evaluates to $3125$, coprime to $N$. This example exercises the twin alternative in the proof. [Check](check_native_digit_variants.py)
