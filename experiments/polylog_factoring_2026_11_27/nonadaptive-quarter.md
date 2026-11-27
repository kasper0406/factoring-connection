# Nonadaptive fixed cutoff queries for two prime divisors

An odd input with exactly two distinct prime divisors can be split by $O(\log N)$ nonadaptive signed disk queries modulo 8 at one common cutoff. The two factors can be prime powers. Unlike shift bisection, all forms can be generated before any residues are available. The forms share a leading coefficient and are diagonal, with the other coefficient varying between queries.

## Read a mixed root by doubling

Let $u$ be the lower of the two mixed roots of 1 modulo $N$. For $c_i=2^i\bmod N$, the roots of $x^2=c_i^2\pmod N$ are exactly $\pm c_i,\pm c_i u$. Put $K=\lfloor N/4\rfloor$. Suppose a query gives parity of the number of root pairs having a representative in $[1,K]$. Subtract the known pair's contribution, leaving a bit $q_i$ for the mixed pair.

Write the binary expansion of $u/N$ as $.b_1b_2\dots$, with $b_1=0$. The fraction has no dyadic ambiguity because $N$ is odd and $u$ is a unit. The fractional part of $2^i u/N$ starts with $b_{i+1},b_{i+2}$. A representative of its signed pair is within the quarter window exactly when those two bits are equal. Hence
$$b_{i+2}=b_{i+1}\mathbin{\oplus}1\mathbin{\oplus}q_i.$$
With $n=N.\mathrm{bit\_length}()$, the $n-1$ queries with $0\le i<n-1$ recover the first $n$ bits. Their dyadic interval has width $2^{-n}<1/N$ and therefore contains a unique fraction with denominator $N$. Recover $u$ exactly by rounding the lower endpoint times $N$ up to an integer. The gcd of $u-1$ with $N$ is a proper prime-power divisor.

## Diagonal encoding

Take $A=(N+1)^3$, $C_0=A(K+1)^2$, $R=N(N+1)$ and
$$B_0=\left\lceil (C_0+AK^2+R)/N\right\rceil.$$
For each $c_i$, choose $d_i\in[C_0,C_0+R)$ with $d_i\equiv-c_i^2\pmod N$ and $d_i\equiv1\pmod{N+1}$. The diagonal form $Ax^2+d_i y^2$ is primitive and reduces to $x^2-c_i^2y^2$ modulo $N$. At the common cutoff $H=NB_0$, the only non-axis rows are $y=\pm1$, each with $|x|\le K$. Indeed $H<2C_0$ since $A(2K+1)>R+N$, and already at $|x|=K+1$ on either row the form is at least $2C_0$. All $|y|\ge2$ have norm at least $4C_0$. The cutoff also satisfies $B_0<2C_0/N<A(K+1)<\sqrt{4Ad_i/3}$.

Use one reference form with the same $A$, cutoff $B_0$, and last coefficient greater than $H$ and congruent to 1 modulo $N+1$. This contains only the axis, which is identical for all queries. Subtraction cancels it, including when $N$ is nonsquarefree. Each signed root pair contributes four points on the remaining rows. Therefore the difference modulo 8 yields exactly the requested parity. All coefficients and the cutoff have $O(\log N)$ bits.

## Checks and scope

Direct enumeration of the signed disks split 215 inputs using 2,158 queries, including five inputs with a repeated prime factor. Run `python nonadaptive_quarter.py`. Enumeration is a correctness control and still costs a power of $N$. With more than two distinct primes, these bits give the XOR of the adjacent-bit differences across the mixed pairs, so this decoder does not recover an individual root in that setting.
