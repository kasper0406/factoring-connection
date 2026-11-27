# Excluding public roots preserves the quadratic shift ideal

Let $N$ be odd with at least two distinct prime divisors, and let $F(X)=(X-b)^2-a^2$ for a unit $a\bmod N$. A mixed root is a root of $F$ modulo $N$ other than the two public roots $b\pm a$. For an integer $m\ge1$ and a polynomial $G\in\mathbb Z[X]$ of degree less than the smallest prime factor of $N$, the following are equivalent.

1. $N^m\mid G(x)$ for every integer lift $x$ of every mixed root.
2. $G\in(N,F)^m$.

Thus excluding the public roots does not enlarge the low-degree lattice of universal shift certificates. A certificate that uses information on the integer representative or only needs to vanish at roots in the queried interval is outside the hypothesis. [Existing small-root decoder](bounded-quadratic-roots.md)

## Proof

Fix a prime power $p^e\parallel N$. At this component both signs occur among the mixed roots, since one can choose an opposite sign at another component. Fix one sign, with integer representative $r=b+a$ or $b-a$. The CRT lets $x$ range over every value $r+p^e t\bmod p^{em}$ while fixing the signs at other components so it remains a mixed root modulo $N$. Under the first condition the polynomial $G(r+p^e t)$ therefore vanishes modulo $p^{em}$ for all integer $t$.

Evaluate at $t=0,1,\ldots,\deg G$. The Vandermonde determinant is a unit modulo $p^{em}$, since every nonzero difference is smaller than $p$. Hence every coefficient of $G(r+p^e T)$ is divisible by $p^{em}$. Expanding at $r$, this is exactly membership in $(p^e,X-r)^m$. The two sign ideals are comaximal because their linear generators differ by $2a$, a unit modulo $p^e$. Their intersection is $(p^e,F)$ by reduction modulo $p^e$, so the intersection of their $m$-th powers is $(p^e,F)^m$.

Now expand $G=\sum_i (u_i+v_iX)F^i$ by repeated division by the monic quadratic. Membership in $(p^e,F)^m$ is equivalent to divisibility of both coefficients at index $i$ by $p^{e\max(0,m-i)}$. To see the forward direction, expand arbitrary polynomial multiples of the generators $p^{e(m-j)}F^j$ in the same basis. The reverse direction follows term by term. Combining the divisibilities for all prime powers gives $N^{\max(0,m-i)}$ in place of $p^{e\max(0,m-i)}$, which is exactly membership in $(N,F)^m$. Conversely, each generator of $(N,F)^m$ evaluates to a multiple of $N^m$ at any root modulo $N$.

## Consequence for interval certificates

If a nonzero $G$ of degree $d$ has the universal mixed-root property, its leading coefficient is divisible by $N^{\max(0,m-\lfloor d/2\rfloor)}$ by the same monic expansion. This holds after any integer translation of $X$. Consequently at scaling radius $X_0\ge\sqrt N$, its scaled leading coefficient alone has absolute value at least $N^m$. The usual sufficient certificate by the sum of absolute scaled coefficients cannot prove the strict bound $|G(x)|<N^m$ on such an interval. Omitting the two public roots therefore leaves that barrier intact. [Related bivariate coefficient test](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_aggregate/native_cap.md)

There is also a direct check for a polynomial promised to vanish modulo $N$ at one mixed root. Its evaluations at the two public roots are each divisible by at least one prime-power component where its sign agrees. If either evaluation is nonzero modulo $N$, its gcd with $N$ therefore splits the input. If both evaluations are zero modulo $N$, monic division gives $F\mid G$ in $(\mathbb Z/N\mathbb Z)[X]$, because the degree-one remainder has two roots whose difference is a unit.

## Checks

The exact diagnostic checked 1,800 polynomials across five moduli, including a repeated-prime input. It compared the monic expansion criterion against evaluation on all mixed-root lifts modulo $N^m$ for $m=1,2$ and three affine maps. Run `python check_mixed_root_certificates.py`.
