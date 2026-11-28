# Factoring the binary digit polynomial of a balanced semiprime

Choose two distinct $b$-bit primes $p,q$ uniformly and put $N=pq$. Write $U_N(X)=\sum_j n_j X^j$ for the binary digits of $N$. The probability that $U_N$ is reducible over $\mathbb Q$, or equivalently that factoring $U_N$ over the integers and evaluating its factors at $2$ splits $N$, is at most
$$
 b^{O(1)}(3/4)^{b/2}
 = N^{-\log_2(4/3)/4+o(1)} = N^{-0.103759\ldots+o(1)}.
$$
The exponent comes from two events with separate proofs: a nontrivial polynomial with the same autocorrelations, and a common large prime divisor of an integer and its binary reversal. This result concerns the unmodified digit polynomial over $\mathbb Z$. [Carry-corrected polynomials](carry-polynomial.md) · [Public multipliers and other variants](native-digit-variants.md)

## Why every nontrivial native factorization splits the integer

Every zero $z$ of a monic binary polynomial satisfies $\Re z<3/2$. In the root bound in the proof of Brillhart–Filaseta–Odlyzko, Theorem 3, the lower coefficients have magnitude at most one. The resulting bound is $\max\{(1+\sqrt5)/(2\sqrt2),r_2\}$, where $r_2$ is the positive root of $x^3-x^2-1$. Both numbers are less than $3/2$. Degree one is immediate. Since $U_N(1)>0$, every nonconstant monic integer factor $G$ has
$$
 |G(2)|>|G(1)|\ge1
$$
by multiplying $|2-z|>|1-z|$ over its roots. Also $G(2)>0$, since $U_N$ has no positive real root and $G$ is monic. Thus if $U_N=AB$ with both factors nonconstant and monic, their integer values at two are precisely $p,q$ in some order. Gauss's lemma gives such a factorization whenever $U_N$ is reducible over $\mathbb Q$. The same inequality gives the prime-input criterion in their Corollary 2. [Brillhart, Filaseta and Odlyzko, pp. 1055–1058](https://doi.org/10.4153/CJM-1981-080-0)

## A polynomial twin or a shared reversal prime

First take a uniformly random odd word $M$ of degree $d$, so the endpoints of $P=U_M$ are one and the $d-1$ interior bits are independent fair bits. Write $P^*(X)=X^dP(X^{-1})$ and $\operatorname{Rev}(M)=P^*(2)$. Kogler's Theorem 1.3 bounds the probability of a different binary word $Q\notin\{P,P^*\}$ with $QQ^*=PP^*\bmod 4$ by $8(3/4)^{\lfloor(d-1)/4\rfloor}$. [Kogler, Theorem 1.3](https://arxiv.org/html/2609.26771v2)

For a monic factorization $P=AB$, the two constant coefficients are each $+1$. They multiply to one, and a monic factor with constant coefficient $-1$ would have a positive real root. As in Kogler's Lemma 3.1, take $Q=AB^*$. It is monic with constant coefficient one and satisfies $QQ^*=PP^*$ over the integers. Its coefficient sum equals $P(1)$. Its sum of squared coefficients also equals $P(1)$, by the middle coefficient of $PP^*$ and the binary coefficients of $P$. Consequently $\sum_j q_j(q_j-1)=0$, proving that every integer $q_j$ is zero or one. Cancellation gives $Q=P$ exactly when $B=B^*$, and $Q=P^*$ exactly when $A=A^*$. Thus, absent a nontrivial twin, every such split has a reciprocal factor. On $M=N$, that factor divides $P^*$ as well as $P$, so its evaluation gives a $b$-bit prime common to $N$ and $\operatorname{Rev}(N)$.

## Counting the shared prime

Fix $d=2b-2$ or $2b-1$, the possible degrees of $U_{pq}$. In a uniform odd word of length $d+1$, choose the central contiguous block of the largest width $w\le b-1$ with $w\equiv d+1\bmod 2$. It contains $l=\lfloor w/2\rfloor=b/2+O(1)$ mirrored pairs. Retain every bit outside the block, the sum of each mirror pair inside it, and any central singleton. These data determine $H=M+\operatorname{Rev}(M)$, since swapping the two members of any mirror pair leaves its contribution to $H$ fixed.

There are $O(1)$ primes $r\ge2^{b-1}$ dividing $H<2^{d+2}$. Conditional on the retained data, the $2^J$ orientations of the $J$ pairs with sum one are equiprobable. For each eligible odd $r$, every orientation gives a different residue $M\bmod r$: the difference between two words is $2^s h$, where $s$ is the offset of the block and $0<|h|<2^w\le r$. Hence the conditional probability that any eligible $r$ divides $M$ is at most $O(2^{-J})$. Before conditioning, $J$ is binomial with parameters $l$ and $1/2$. Averaging gives
$$
 \Pr\{\exists\ b\text{-bit prime } r\mid M,\ r\mid\operatorname{Rev}(M)\}
 \ll \mathbb E[2^{-J}]=(3/4)^l.
$$

For the transfer, count the exceptional words in both degree bands. Together they contain $O(2^{2b})$ words. By the prime number theorem there are $\Theta(2^{2b}/b^2)$ distinct products of two different $b$-bit primes, each with the same sampling probability. The exceptional density therefore increases by at most $O(b^2)$. The twin and reversal bounds both give $O((3/4)^{b/2})$ before this transfer. Palindromes are already covered by the shared-prime event.

## Transfer from irreducibility of the noncyclotomic part

Kogler's Theorem 1.1 gives the unconditional reducibility probability $\sqrt{2/(\pi d)}+O(d^{-1})$ for the uniform odd binary word. That estimate alone loses its force under the $O(b^2)$ transfer. Theorem 1.2 gives an alternative proof with a weaker final rate: outside a set of probability at most $C\exp(-cd/(\log d)^4)$, $P=\Phi I$ with $\Phi$ a product of cyclotomic polynomials and $I$ irreducible and noncyclotomic. [Kogler, Theorems 1.1–1.2, September 2026 preprint, revised October 2026](https://arxiv.org/html/2609.26771v2)

If $\Phi_m\mid U_N$, then $m\ge2$ because $U_N(1)>0$. Every primitive $m$-th root $\zeta$ satisfies $|2-\zeta|>1$, so $\Phi_m(2)>1$ and some $r\in\{p,q\}$ divides $2^m-1$. Moreover $\varphi(m)\le d$ and $\varphi(m)\ge\sqrt{m/2}$, by the prime-power formula for $\varphi$. The product of the integers $2^m-1$ over $2\le m\le2d^2$ has $O(d^4)$ bits and contains at most $O(d^4/b)$ distinct $b$-bit primes. The chance of either hidden prime appearing in this list is $O(b^{O(1)}2^{-b})$. Theorem 1.2 and the degree-band count therefore already imply a reducibility probability $\exp(-\Omega(b/(\log b)^4))$.

Breuillard and Varjú obtained the uniform-word irreducibility asymptotic under the Riemann hypothesis for the relevant Dedekind zeta functions. Their Theorem 2 gives the single irreducible noncyclotomic part with failure $\exp(-\Omega(\sqrt d/\log d))$, and their Corollary 4 gives the leading reducibility term. The present exponential bound for balanced semiprimes uses Kogler's unconditional twin estimate and the additional arithmetic condition $M=pq$. [Breuillard and Varjú, *Acta Mathematica* 223 (2019), §1](https://doi.org/10.4310/ACTA.2019.v223.n2.a1)

## Finite checks

Run `python experiments/polylog_factoring_2026_11_27/check_native_digit_bound.py` from the repository root. It factors the digit polynomial for every distinct prime pair of equal bit length from three through eight, checks the signed factor identities and the twin-or-reversal implication, and verifies injectivity of word residues in every conditional orientation family in both degree bands for each bit length. Among 364 semiprimes, 20 had reducible digit polynomials. Eighteen of those shared a prime with their reversal, and the remaining two produced nontrivial twins. The occupancy check passed for 1,848 eligible-prime cases across 14,444 conditional word families. [Exact check](check_native_digit_bound.py)
