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

### How much do local unit restrictions buy in a lattice determinant?

In a triangular basis of degrees $0,\ldots,D-1$, enforcing the unit restriction at every prime up to $D$ increases the radius allowed by the coefficient determinant by a factor tending to **3.410875**. The comparison is with polynomials already integer-valued on all integers. For an affine order hint, a hidden divisor $r\mid N$ is automatically a unit at each auxiliary prime $\ell\nmid N$. Writing $r=x+y$ merely translates the excluded residue class of $y$. The unknown-divisor vanishing conditions of the order-hint lattice still have to hold. [Factoring with an order hint](factoring-with-an-order-hint.md)

For a set $S$ of auxiliary primes, write $A_D(S)$ for the logarithm of the extra leading-coefficient denominators. The unit ordering above attains the local valuation at each degree $j$, while Legendre's formula gives the unrestricted valuation of $j!$. Consequently
$$
A_D(S)=\sum_{\ell\in S}\log\ell\sum_{j=0}^{D-1}\sum_{k\ge1}
\left(\left\lfloor\frac j{\varphi(\ell^k)}\right\rfloor-
      \left\lfloor\frac j{\ell^k}\right\rfloor\right).
$$
Only $\ell\le D$ contribute to this finite sum. For $C_S=\sum_{\ell\in S}\log\ell/(\ell-1)^2$, the error for every $D\ge2$ satisfies
$$
\left|A_D(S)-\tfrac12 D^2 C_S\right|\le 7D\log(2D)\qquad\text{uniformly in }S.
$$

To see the floor error, put $F(t)=\sum_{v\ge1}\max(D-vt,0)$. At integer $t$ this equals $\sum_{j<D}\lfloor j/t\rfloor$. For $a=\varphi(\ell^k)$ and $b=\ell^k$, subtract $D^2(1/a-1/b)/2$ from $F(a)-F(b)$. The absolute remainder is at most $b/8$, since for $s=D\bmod t$ one has $F(t)=D^2/(2t)-D/2+s(t-s)/(2t)$. It is also at most $(D/2)\log(b/a)\le D/[2(\ell-1)]$. This follows by integrating the slope of $F$ and using $|\lfloor z\rfloor(\lfloor z\rfloor+1)-z^2|\le z$.

For $\ell\le D$, let $K$ be the last $k$ with $a<D$. Summing $b/8$ for $k<K$ costs at most $D/[4(\ell-1)]$ by the geometric series. Use the second bound at $K$. For $k>K$ the floors vanish and the continuous tail is $D^2/[2(\ell-1)^2\ell^K]\le D/[2(\ell-1)]$. The total local error is at most $5D/[4(\ell-1)]$ before multiplying by $\log\ell$. Summing uses $\sum_{\ell\le D}\log\ell/(\ell-1)\le4\log D$, obtained from $\log(D!)\ge\sum_{\ell\le D}\lfloor D/\ell\rfloor\log\ell$. For $\ell>D$ the continuous terms alone contribute at most $D(\log(2D)+1)$ by comparison with the integral of $\log(2t)/t^2$ from $D-1$ to infinity. These bounds give the stated constant.

Scaling the variable by a radius $X$ multiplies the determinant by $X^{D(D-1)/2}$. The corresponding radius multiplier is therefore $\exp(2A_D(S)/[D(D-1)])=\exp(C_S+O(\log D/D))$. For all primes, $C_S=1.2269688\ldots$ is the known constant in Díaz's asymptotic for the prime Bhargava factorial. The unit denominator at degree $j$ is $(j+1)!_{\mathbb P}$ in that convention. The finite multiplier at $D=8192$ is 3.40966908. This bounded improvement leaves the quarter-power exponent unchanged. [Bhargava's generalized factorials](https://doi.org/10.2307/2695734) · [Díaz, Theorem 2](https://arxiv.org/pdf/2010.13645#page=2) · [Exact finite-degree checks](check_unit_factorials.py)
