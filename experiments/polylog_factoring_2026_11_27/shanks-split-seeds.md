# Choosing a cyclic cubic field from a public modular root

For an odd distinct semiprime $N=pq$ with $p,q>3$, a public seed produces a cyclic cubic in which both primes split. Shanks's family is
$$
 f_t(X)=X^3-tX^2-(t+3)X-1,\qquad L_t=t^2+3t+9.
$$
For integral $t$ this polynomial is irreducible over $\mathbf Q$: its only possible rational roots are $\pm1$, and $f_t(-1)=1$, $f_t(1)=-2t-3$. Its discriminant is $L_t^2$. The splitting field is therefore cyclic of degree three. The repository had already used this same cubic and its order-three substitution on a special near-square branch where deciding whether it split was the unknown bit. Here we choose $t$ from a root. [Shanks, *The simplest cubic fields*](https://doi.org/10.1090/S0025-5718-1974-0352049-8) · [Earlier discriminant −59 branch](../ellipse_disk_mod8/cubic59.md)

Choose $x\bmod N$. Check that $x,x+1,x^2+x+1$ are units, returning any proper gcd encountered. Set
$$
 t=\frac{x^3-3x-1}{x(x+1)}\pmod N
$$
and take its centered integer lift. The orbit under $x\mapsto-(x+1)/x$ consists of three distinct roots modulo each prime. The chosen units certify this directly: the fixed points obey $x^2+x+1=0$, and $L_t x^2(x+1)^2=(x^2+x+1)^3\bmod N$. Thus both hidden primes are unramified and split completely. In the Eisenstein algebra with $\omega^2+\omega+1=0$, the same formula gives $(x-\omega)^3=x(x+1)(t-3\omega)$. Its quotient by the conjugate gives a cube in the norm-one group. This distinguishes it from a bare assertion that a rational integer is a cube: every integer modulo a prime $2\bmod3$ has a cube root, while complete splitting in a cyclic cubic cuts out one third of the unit residue classes at a prime conductor.

## When a residue list is explicit

If $L_t$ is a rational prime $\ell$, the field discriminant is $\ell^2$. Indeed its only other possibility by the power-basis index formula is 1, excluded by Minkowski's bound. The conductor-discriminant formula for a cyclic cubic then gives conductor $\ell$. The field is the unique cubic subfield of $\mathbf Q(\zeta_\ell)$. For a prime $s\ne\ell$ its splitting test is $s^{(\ell-1)/3}=1\bmod\ell$. This gives an explicit list of $(\ell-1)/3$ admissible unit residues. In particular, the size of $\ell$ matters even when its primality is certified. For composite $L_t$ the actual field conductor can be smaller, and the root construction alone does not give a prime-conductor list.

At $N=713$, the public seeds $7$ and $19$ give $t=-7$ and $-14$, respectively. Their $L_t$ are the primes $37$ and $163$. The restrictions allow 12 and 54 unit residues, compared with 36 and 162 before the tests. CRT therefore gives $648$ classes modulo $6031$ instead of $5832$ unit classes. Both primes $23,31$ pass. The calculation constructing the two fields uses only $N$ and the seeds. The checker compares roots with the residue test on 1,310 prime fixtures. [Executable example](check_shanks_split_seeds.py)

## Why targeting small coefficients costs work

For every unramified $t$ that splits at both primes, exactly nine residues $x\bmod N$ map to it. There are three roots locally at each prime. No root can equal $0$ or $-1$, and an unramified root cannot be fixed by the order-three action. Consequently a uniform seed finds any of $B$ specified centered parameters with probability at most $9B/N$ per draw. Conditional on passing all the unit checks, it is uniform on the attainable parameters. At $s=p,q$, there are $s-4$ eligible seeds for $s=1\bmod3$ and $s-2$ for $s=2\bmod3$, by removing $0,-1$ and the primitive cube roots of 1. The checker exhausts all seeds on 36 composites.

There is likewise a height obstruction even for deliberately small rational seeds. Write $x=a/b\bmod N$ with $b$ a unit, $|a|,|b|\le H$, and suppose that the integral lift of the parameter satisfies $|t|\le T$. Then $N$ divides the nonzero integer $a^3-ta^2b-(t+3)ab^2-b^3$. Its nonvanishing follows from irreducibility and its absolute value is at most $(2T+5)H^3$. Both heights therefore cannot be polynomial in $\log N$ as $N$ grows. These bounds leave open a targeted way to compute a large root for a small field or recognize a small conductor at a large parameter. Simply repeating independently generated large-conductor fields need not produce a short residue list: the simultaneous kernel for independent prime conductors $\ell_i$ has $\prod_i((\ell_i-1)/3)$ classes.

## Scope of the rational map

Suppose a single-variable rational map $T=R(X)$ over $\mathbf Q$ generates an extension $\mathbf Q(X)/\mathbf Q(T)$ which is itself cyclic Galois. The deck group lies in $\mathrm{PGL}_2(\mathbf Q)$ because every rational automorphism of $\mathbf Q(X)$ is fractional linear. If a representing matrix $M$ has finite projective order $e>1$, its eigenvalue ratio is a primitive $e$th root $\zeta$ and $\zeta+\zeta^{-1}=\operatorname{tr}(M)^2/\det(M)-2$ is rational. This rational algebraic integer belongs to $\{-2,-1,0,1,2\}$, so $e$ belongs to $\{2,3,4,6\}$. The extension degree equals $e$. Hence this particular way to generate a cyclic splitting certificate by choosing one unrestricted rational parameter has bounded degree. Composita, covers with nonrational source and specialization-specific Galois groups are outside the argument. [Beauville, *Finite subgroups of* $\mathrm{PGL}_2(K)$, §1](https://math.univ-cotedazur.fr/~beauvill/pubs/PGL%282%29.pdf)

Using a cyclotomic base of dimension $\ell-1$ gives a different construction with arbitrary prime degree, as well as explicit bounds on random candidate rejection. [Higher-degree certificate](cyclic-kummer-certificates.md)
