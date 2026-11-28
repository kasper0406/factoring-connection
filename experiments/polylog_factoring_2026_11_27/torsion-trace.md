# What the trace pairing sees in elliptic torsion

Let $N=pq$ for distinct odd primes, and let $E$ have good reduction at both. The factors of a division polynomial can have different degrees at $p$ and $q$ even when its trace pairing has full rank and square determinant at both. The earlier Schoof investigation identified the missing local Frobenius operator. Here we check what taking its permutation sign can retain. Schmidt gives exact formulas for multiplication-polynomial resultants and discriminants. The argument below uses only the discriminant square class. [Earlier Schoof investigation](../../notes/user_ideas/2026-10-18_schoof_semigroups_and_cross_difference.md) · [Schmidt, *Journal of Number Theory* (2015)](https://doi.org/10.1016/j.jnt.2014.10.007)

## Odd torsion and the determinant character

Fix an odd prime $\ell\nmid N$. The roots of the monic division polynomial $f=\psi_\ell/\ell$ are the $x$-coordinates of the nonzero $\ell$-torsion points, indexed by $(\mathbb F_\ell^2\setminus\{0\})/\{\pm1\}$. Over either hidden field $\mathbb F_r$, the algebra $A=\mathbb F_r[x]/f$ is étale. Its trace pairing $\operatorname{Tr}_A(uv)$ has rank $(\ell^2-1)/2$. Frobenius acts on a Vandermonde square root of the discriminant by the sign of its permutation of the roots.

More generally, suppose a torsion configuration is obtained from a fixed finite $\mathrm{GL}_2(\mathbb Z/m)$-set by the Galois action on $E[m]$, with $m$ odd and coprime to $N$. Its permutation sign factors through the determinant. Indeed elementary unipotents have odd order and generate $\mathrm{SL}_2(\mathbb Z/m)$. Hence the local discriminant sign is $\chi(r)$ for a character $\chi:(\mathbb Z/m)^\times\to\{\pm1\}$ determined by that set. Here the Weil pairing gives $\det(\mathrm{Frob}_r)=r\bmod m$. This applies to full configurations such as tuples and flags. Choosing a curve-specific Galois-stable subset requires a separate analysis. [Schoof, *Mathematics of Computation* (1985)](https://pages.cs.wisc.edu/~cs812-1/Schoof85.pdf)

For the $x$-coordinate set, write a torsion vector as $(u,v)$ up to sign and take $g$ primitive modulo $\ell$. On $v\ne0$, $\operatorname{diag}(g,1)$ has one even-length cycle in each of $(\ell-1)/2$ rows. The row $v=0$ is a single cycle of length $(\ell-1)/2$. The total sign is $-1$, so $\chi(r)=(r/\ell)$. By quadratic reciprocity the determinant has the square class of $(-1)^{(\ell-1)/2}\ell$ in $\mathbb F_r$. Rank and determinant classify nonsingular symmetric forms over this field, so the local isometry class is independent of $E$.

## An explicit pair of reductions

Take $N=703=19\cdot37$ and $E:y^2=x^3-3x$. Its monic 3-division polynomial is $f=x^4-6x^2-3$, with discriminant $-2^{12}3^3$. Modulo 19 it factors as $(x^2+6x-4)(x^2-6x-4)$, with both quadratics irreducible. Modulo 37 its four roots are $\pm11,\pm12$. Both primes are 1 modulo 3, so both trace forms have square determinant despite the different numbers of rational roots.

## Twisting the pairing gives a point-order test

Multiplication by $x-t$ changes the determinant of the trace pairing by $\operatorname{Norm}_A(x-t)=(-1)^{\deg f}f(t)$. For $E:y^2=x^3+ax+b$, put $c=t^3+at+b$. If $c$ is a unit modulo $N$, the public twist $Y^2=X^3+c^2aX+c^3b$ contains $P'=(ct,c^2)$. Over a field containing $\sqrt c$ this is the image of $(t,\sqrt c)$ on $E$. Thus the twisted pairing is singular at $r$ exactly when $[\ell]P'=O$ there. A proper $\gcd(c,N)$ already splits the input. If $c=0$ at both primes, the point has order two at both and the odd division polynomial is a unit there. [Lenstra, *Factoring integers with elliptic curves* (1987)](https://pages.cs.wisc.edu/~cs812-1/Lenstra1987.pdf)

In the example, choosing the known root $t=11$ gives $\gcd(703,f(11))=37$, and $c=1298$ is a unit. This verifies the equivalence on a successful witness. Finding such a witness from $N$ is the remaining point-order problem.
