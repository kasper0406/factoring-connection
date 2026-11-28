# Certifying unequal hidden cubic labels in an Eisenstein torus

For $N=pq\equiv-1\pmod3$ with distinct odd primes, a cubic Jacobi symbol computed from $N$ alone decides whether a norm-one point has **different cubic endpoint labels at the two primes**. At $N\equiv2,5\pmod9$ the public point $\omega$ already has unequal labels. At $N\equiv8\pmod9$ a uniform norm-one point succeeds with probability $2/3$. This supplies a witness whose local powers disagree. Computing their combined value modulo $N$ remains a separate operation. [Native symbol and probe](eisenstein_torus.py) · [Exact checks](check_eisenstein_torus.py)

## The local-to-global formula

Let $\omega^2+\omega+1=0$ and $T_r=\{z\in(\mathbf F_r[\omega])^\times:z\bar z=1\}$ at each $r\mid N$. For $r\ne3$, put $\epsilon_r=1$ when $r=1\bmod3$ and $-1$ otherwise. The group $T_r$ is cyclic of order $r-\epsilon_r$. For any $z\in T_N=T_p\times T_q$ define $L_r\in\{0,1,2\}$ by
$$
z_r^{(r-\epsilon_r)/3}=\omega^{L_r}.
$$
Write $J(z)\in\mathbf F_3$ for the exponent of the cubic Jacobi symbol $(z/(N))_3$. Then
$$
J(z)=-\epsilon_p L_p-\epsilon_q L_q\pmod3 .
$$

At a split prime choose $w$ as the image of $\omega$ in one residue field. The two coordinates of $z$ are $s,s^{-1}$. If $s^{(r-1)/3}=w^{L_r}$ then $(s^{-1})^{(r-1)/3}=(w^{-1})^{L_r}$, so **both** prime ideals give the same exponent. Their sum $2L_r$ is $-L_r\bmod3$. At an inert prime the residue field has size $r^2$. Its Euler exponent $(r^2-1)/3$ is $(r-1)(r+1)/3$, and $r-1=1\bmod3$, so its contribution is $L_r$.

| Public residue | Meaning of $J(z)=0$ | Unequal labels conditional on $J$ for uniform $z$ |
|---|---|---|
| $N=2\bmod3$ | $L_p=L_q$ | Impossible at zero and certain at either nonzero value |
| $N=1\bmod3$ | $L_p=-L_q$ | Probability $2/3$ at **each** of the three values |

Indeed each label is independently uniform on $\mathbf F_3$ under the uniform torus distribution. In the second row each prescribed sum has three solutions, exactly one with equal coordinates. A cubic Jacobi exponent only tests this one relation. For example the first row detects a noncube modulo at least one prime whenever $J\ne0$, without naming that prime.

## Public sampling and rotations

Choose $w=a+b\omega$ uniformly modulo $N$ and test the gcd of $N$ and $a^2-ab+b^2$. A proper gcd is a factor. If the norm is a unit, use $z=w^2/\operatorname{Norm}(w)=w/\bar w$. The image is uniform on $T_N$ because this is a surjective group homomorphism with equal-size fibers in each residue algebra. Surjectivity in the split case follows by dividing the two coordinates. In the inert case the map is $w\mapsto w^{1-r}$ and its image has order $r+1$.

The supplement $J(\omega)=(N^2-1)/3\bmod3$ gives the public rotation rule $J(z\omega^i)=J(z)+iJ(\omega)$. For mixed prime signs the table for $N=2,5,8\bmod9$ gives $J(\omega)=1,2,0$ respectively. In the first two cases every rotation triple has two unequal-label witnesses, including the direct choice $z=\omega$ when starting from 1. In the third case rotations leave the original decision fixed. Sampling conditional on a unit then needs $3/2$ torus points in expectation. The element $z=\omega$ itself has the **same** literal component at either prime. Its hidden powers with exponents $(r-\epsilon_r)/3$ give the certified difference.

The implementation computes the symbol by Euclidean division in $\mathbf Z[\omega]$ without factoring $N$. For a primary denominator $c+d\omega$ ($3\mid d$, $3\nmid c$), remove powers of $1-\omega$ from the numerator remainder and multiply by a power of $\omega$ to make that remainder primary. The two supplementary exponents for $1-\omega$ and $\omega$ are $-(c^2-1)/3$ and $(c^2-cd-1)/3$. Cubic reciprocity exchanges the primary pair. Rounding the two quotient coordinates gives a remainder of norm at most three quarters of the previous denominator's norm, giving polynomial bit complexity. [Bach and Sandlund, equations (2)–(6)](https://arxiv.org/pdf/1807.07719)

The checks matched 5,103 arbitrary numerator and denominator cases against prime-field Euler symbols, including repeated prime factors of the denominator and nonunits. Enumerating 6,696 torus points over 28 products of two primes verified every entry of the table. Three products of 60-bit primes exercised the input-only probe at each possible mixed residue modulo 9. At the $8\bmod9$ example, the public symbol exponent was 2, the independent prime-factor check gave labels $(0,1)$, and the sampled point's cube differed from 1.

## Which shifts give a second reading?

Conjugating the Eisenstein numerator at a rational denominator negates its exponent, so $J(w/\bar w)=2J(w)=-J(w)$. In particular for $c=\pm\omega^j$ and unit $1+cz$, its quotient by its conjugate equals $cz$, hence
$$
J(1+cz)=-J(c)-J(z).
$$
There were 17,526 independent unit-shift checks. More generally let $f^*(X)=X^d\bar f(1/X)=\lambda f(X)$ for a unit $\lambda$ and polynomial $f$ of degree $d$. At a norm-one $z$ with $f(z)$ a unit, $f(z)/\overline{f(z)}=\lambda^{-1}z^d$, and therefore $J(f(z))=J(\lambda)-dJ(z)$. These shifted symbols contain the same one relation. Nonunits can instead be checked for a proper norm gcd. Other rational functions and integer representatives fall outside this identity.

The archive's earlier fixed-numerator cubic symbols reduce to rational residue classes of $N$. An earlier Cardano numerator had cubic Jacobi exponent zero on its specified input family. Here the numerator is an actual sampled norm-one point, and the criterion concerns its own two labels. The existing Lucas analysis also gives exponent identities for norm-one groups and a route to a split when a suitable local order multiple is supplied. Those results and this symbol do not provide the combined local endpoint as an element of the public ring. [Fixed numerators](../../notes/newfamily_results.md) · [The Cardano example](../ellipse_disk_mod8/cubic59.md) · [Earlier norm-one exponent identities](../../notes/primality_analysis.md)
