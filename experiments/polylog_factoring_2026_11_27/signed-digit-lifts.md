# Signed factors of an ordinary binary digit polynomial

Write $U_n(X)$ for the polynomial with the binary digits of $n$ as coefficients, so $U_n(2)=n$. Factoring $U_n$ over $\mathbb Z$ and evaluating the factors at 2 is a known algorithm. Breuillard and Varjú describe it in §1.1 of their 2018 preprint, published in 2019. The example below distinguishes this integer calculation from searching for the *binary digit polynomials of the designated factors* as divisors over $\mathbb F_2$. [Breuillard–Varjú](https://arxiv.org/pdf/1810.13360v1#page=2) · [Native digit polynomial bound](native-digit-polynomial-bound.md) · [Carry parity](carry-parity.md)

For $n=2{,}006{,}801$, direct integer polynomial factorization gives $U_n=fg$ with
$$
\begin{aligned}
f&=X^{10}+X^7-X^6+X^5+X^2-X+1,\\
g&=X^{10}+X^9+X^8-X^3+X+1.
\end{aligned}
$$
Both factors are irreducible over $\mathbb Z$, and their values at 2 are the primes 1123 and 1787. Over $\mathbb F_2$, dividing $U_n$ by $U_{1123}$ leaves $X^9+X^7+X^5+X^4+X^3+X^2+1$. Dividing by $U_{1787}$ leaves $X^5+X^4+X^3+X^2+X$. Thus neither designated binary factor can be obtained by grouping factors of $U_n\bmod 2$. The signed polynomials $f,g$ themselves still reduce to divisors there.

### An infinite family

For $t$ a power of two, let $a=f(2^t),\ b=g(2^t)$, and $N_t=ab$. Binary place values give
$$
U_{N_t}(X)=U_n(X^t)=f(X^t)g(X^t).
$$
In fact the displayed factors are the two irreducible factors over $\mathbb Q$ for every power of two $t$. Exact factorization over $\mathbb Z$ verifies the finite base cases $f(X^4)$ and $g(X^4)$, both irreducible of degree 40. For a root $\alpha$ of either original polynomial, the field-degree tower gives irreducibility of $Y^4-\alpha$ over $\mathbb Q(\alpha)$. Capelli's binomial criterion says this excludes $\alpha$ being a square or $-4$ times a fourth power in that field. The same two exclusions imply $Y^{2^j}-\alpha$ is irreducible there for every $j\ge1$. Applying the tower again gives the claim. Integer polynomial factorization from $N_t$ alone therefore returns precisely the two designated values $a,b$ on evaluation at 2. Those integers can themselves be composite. [Capelli's criteria, Jones and Phillips §2](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/9B22D4A35CB8C35FCE9EBB8BEEF4A38F/S0004972717000648a.pdf/an-infinite-family-of-ninth-degree-dihedral-polynomials.pdf)

For every such $t\ge4$, neither $U_a$ nor $U_b$ divides $U_{N_t}$ in $\mathbb F_2[X]$. Indeed $U_{N_t}=U_n(X)^t\bmod2$, so the product of its distinct irreducible factors has degree at most 20. Put $T=X^t$ and $S=1+X+\cdots+X^{t-1}=(1+T)/(1+X)$. Expanding the binary blocks in the two integer factors gives
$$
U_a=A+(T+T^6)S,\quad A=1+T^5+T^{10},\qquad
U_b=S P,\quad P=(1+X)(1+T^9)+T^3+\cdots+T^7.
$$

Since $T'=0$, the derivative of $(1+X)U_a$ equals $A$. This is coprime to $(1+T)(T+T^6)$, because $T+T^6=T(1+T^5)$ and $A=1$ at $T=0$ or $T^5=1$. Hence $(1+X)U_a$, and thus $U_a$, is squarefree. Its degree $10t$ exceeds 20. Likewise $P'=1+T^9$ is coprime to $T^3+\cdots+T^7=T^3(1+T+\cdots+T^4)$. A common nonzero root would have order dividing both 9 and 5, and the sum of five terms at 1 equals 1. Thus the squarefree divisor $P\mid U_b$ has degree $9t+1>20$. Both exclusions follow.

The check also verifies the two nondivisibilities directly at $t=1,2$, and for every power of two through 128. It runs integer polynomial factorization using only $N_t$ through $t=16$, giving verified splits at 21, 41, 81, 161 and 321 bits. At the original semiprime a zero carry-parity template returns no factor. At the dilations the theorem excludes the two designated binary divisors and leaves other integer divisors open. [Runnable checks](check_signed_digit_lifts.py)
