# Which residues of the real divisor sum can we compute from the product?

The two-reading reduction uses $H(m)=m^2D_m(3)$. We tested a concrete supplier on $N=pq$ with distinct primes greater than three. Given a modulus $M$ coprime to $N$, let the hidden factors vary over **all units of product $N$ modulo $M$**. If all completions give the same $H(N)$, this residue is available from $N$ alone. For example, $N=3239$ gives $H(N)\equiv75\pmod{1152}$. The maximal such modulus is small enough to compute from $N\bmod48$. [Earlier divisor-sum reduction](mellin-four-squares.md#two-real-readings-per-split) · [Implementation](mellin_local_residues.py)

| Class of $N$ | Largest forced power of two |
|---|---:|
| $1\bmod4$ | 8 |
| $3\bmod8$ | 32 |
| $15\bmod16$ | 64 |
| $7\bmod16$ | 128 |

Multiply the table entry by 9 when $N\equiv2\pmod3$. The resulting $M$ is the largest possible coprime modulus in this model, and the residue is $3(N^2+N+1)\bmod M$.

## Why the table ends there

Put $f(x)=x^2+x+1$, let $x,y$ be local units with $xy=N$, and set $S=x+y$. The discrepancy from the reference completion $(1,N)$ is
$$
\Delta=f(x)f(y)-3f(N)=-(x-1)(y-1)(S+2N+2).
$$
For odd $x,y$ each factor on the right is even. If $N\equiv3\bmod4$, take $x\equiv1\bmod4$ by interchanging the factors. Then $4\mid x-1$ and $4\mid S$. For $N\equiv7\bmod8$, the identity $xS=x^2+N$ gives $8\mid S$. When also $N\equiv7\bmod16$, either $x\equiv1\bmod8$, giving another factor 2 in $x-1$, or $x\equiv5\bmod8$, giving $16\mid S$. These prove the four lower bounds. In table order, the choices $x=-1,5,5,9$ give exactly the stated valuations of $\Delta$ with $y=N/x$ in the 2-adic units.

If $N\equiv1\bmod3$, taking $x=1$ and $-1$ gives different answers modulo 3. Otherwise one unit is 1 modulo 3 and its $f$-value is 3 modulo 9. The other has $f$-value 1 modulo 3. At $x=4$, $16\Delta=9(4-N)(8+3N)$ has exact 3-adic valuation two. For a prime $\ell\ge5$, clear the denominator to get
$$
x^2\Delta=(x-1)(x-N)\big[x^2+2(N+1)x+N\big].
$$
This monic quartic cannot vanish on all the units when $\ell\ge7$. Modulo 5 it would have to equal $x^4-1$, forcing both $N=-1$ and $N^2=-1$ modulo 5. The Chinese remainder theorem now gives the maximum. For any fixed congruence modulus the prime number theorem in arithmetic progressions realizes the contrasting local completions with distinct primes in common intervals $[X,2X]$ for large $X$. They give different integers $N$ with the same observed congruence. Arithmetic on the full input can use more information.

## What raising the input to a power buys

There is also a polynomial-time procedure in $j$ and $\log N$ for the **largest** coprime modulus forced for $H(N^j)$, for every integer $j\ge1$. With $s=j+1$, put
$$
B_j(X)=X^{j+1}+(X+1)\sum_{a=0}^{j-1}X^a.
$$
The desired value is $N^{j-1}B_j(p)B_j(q)$. After removing the unit $N^{j-1}$, subtract $B_j(1)B_j(N)$ and multiply by $X^s$ with $q=N/X$. The result is a monic integer polynomial of degree $d=2s$, since $B_j$ is monic with constant term one.

For each prime $\ell$, take the first $d+1$ positive integers $x_i$ coprime to $\ell$. The units cycle in blocks of length $\varphi(\ell^a)$ modulo every $\ell^a$. Thus every Lagrange denominator at these nodes has valuation at most $v=\sum_{a\ge1}\lfloor d/\varphi(\ell^a)\rfloor$. If all the values of the monic polynomial vanished modulo $\ell^{v+1}$, the interpolation formula for its leading coefficient would make 1 divisible by $\ell$. Only primes $\ell\le d+1$ require testing.

The same nodes give the exact valuation, up to this cap. In each residue class modulo $\ell^a$ the node counts differ by at most one. At any unit argument $u$, the numerator of the $i$-th Lagrange weight contains at least as many factors divisible by $\ell^a$ as its denominator. If $u\equiv x_i$, these counts agree. Otherwise the denominator excludes $x_i$ from its class, while the numerator retains every node congruent to $u$. All weights are therefore $\ell$-adically integral. Taking the gcd of the node values with $\ell^v$ captures exactly the forced prime power. This is fixed-divisor interpolation on the units. [Bhargava on fixed divisors and generalized factorials](https://doi.org/10.1006/jnth.1998.2220)

Since $\varphi(\ell^a)\ge\ell^a/2$, the combined modulus divides $(2d)!$ by Legendre's formula. Its bit length is $O(j\log(2j))$. For example, at $N=3239$ and $j=10$ the maximal residue is 21 modulo 10,080. With $j$ polynomial in $\log N$, these individual moduli are too small for direct CRT reconstruction from the size bound $N^{2j}\le H(N^j)\le N^{3j}$. This comparison concerns reconstructing one powered value from forced congruences. It leaves open other computations on $N$ or combinations of different powers.

The checks compare 1,965 balanced prime pairs, exhaust all unit completions for 212 local product classes at $j=1$, and test the general interpolation method against 270 exhaustive powered-modulus comparisons. [Reproducible checks](check_mellin_local_residues.py)
