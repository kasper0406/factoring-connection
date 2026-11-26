# Five residues on one fixed disk suffice to factor

For every balanced product of distinct odd primes, the minimum search in the [weighted-disk reduction](README.md) can be replaced by five polynomial moments at a single public cutoff. The decoder takes polynomial time in the input bit length. The forward evaluator in [weighted_recurrence.py](weighted_recurrence.py) uses the existing residue join, which still scans up to $2\sqrt{NA}$ long-coordinate values. Producing the five residues efficiently from $N$ alone remains open.

## Observable

Use the [public construction](README.md#public-branch-construction) with $\ell=1$, $\Delta=4N-c^2>4$ and $A=\lfloor\sqrt{\Delta/3}\rfloor<p$. A nonunit $c$ already yields a divisor. Sum over every positive pair $(w,y)$ in the one disk $w^2+\Delta y^2\le 4NA$ satisfying $4N\mid w^2+\Delta y^2$. Write $j=(w^2+\Delta y^2)/(4N)$. The five raw residues modulo $N$ have weights

$$1,\quad j,\quad j^2,\quad y^2,\quad jy^2.$$

These are exact integer quotients on the summation domain. To express the task using only polynomials in $w,y$ with integer coefficients, a weight containing $j^h$ can equivalently use the numerator $(w^2+\Delta y^2)^h$ modulo $(4N)^h N$, followed by exact division by $(4N)^h$. This costs $O(\log N)$ bits per residue. The observable is stronger than the original unweighted count modulo 16.

Subtract the public pairs $(kc,k\ell)$ through $K=\lfloor\sqrt A\rfloor$. Their contributions are $(K,P_2,P_4,\ell^2P_2,\ell^2P_4)$, where $P_i=\sum_{k=1}^K k^i$ have fixed-degree closed forms. Call the resulting residues $(M,J_1,J_2,Y_0,Y_1)$. Under the [mixed-orbit correspondence](../ellipse_row_batches/README.md#one-global-second-moment-suffices-to-factor), they count unsigned vectors in one mixed lattice, whose norm divided by $N$ is a primitive integral positive form $F$ of discriminant $-\Delta$. Each vector retains its physical coordinate $y$ up to sign. Conjugation switches the two mixed lines, so positive pairs count precisely one copy of each antipodal pair in that lattice.

## Shape at the cap

For independent $u,v$ of norms at most $A$, the Gram identity gives

$$\Delta\det(u,v)^2=4F(u)F(v)-[F(u+v)-F(u)-F(v)]^2\le4 A^2\le4\Delta/3.$$

Here the determinant uses a basis of the scaled lattice. It follows that $u,v$ form a basis. Any third primitive direction must be $u+v$ or $u-v$ up to sign, because both of its determinants against the basis vectors have absolute value one. Let the three norms be $s\le t\le z$. The same identity, now symmetric in all three norms, gives

$$\Delta=2(st+tz+zs)-s^2-t^2-z^2\le3z^2\le3A^2\le\Delta.$$

For the first inequality, fixing $t,z$ and raising $s$ to $t$ strictly increases the expression unless $s=t$. With $s=t$ the expression is $4tz-z^2\le3z^2$. Thus $s=t=z=A$. In the basis $u,v$ the form would be $[A,\pm A,A]$. Primitivity forces $A=1$, contradicted by $\Delta>4$. Primitivity here also follows locally. At a prime dividing $N$ the discriminant is a unit. At any other prime the index-$N$ lattice change and the scale $N$ preserve the primitive content of the public form.

If two independent directions occur, a nontrivial multiple of either cannot fit. Indeed the primitive vector on that direction would have norm $\le A/4$, and pairing it with a vector of norm $\le A$ in the Gram identity would give $\Delta\le A^2$. The disk therefore has one or two points up to sign, or it has exactly the multiples $u,2u,\ldots,Mu$ of a single primitive vector. In every case $1\le M\le\max(2,\sqrt A)<N$, so the count is recovered from its residue.

## Decode the first norm and its moment

Write $a$ for the minimum norm and $S$ for the sum of the physical $y^2$ on it. The [existing decoder](../ellipse_row_batches/README.md#one-global-second-moment-suffices-to-factor) proves $0<S<N$ and extracts a verified divisor from $(a,S)$, including ties. The fixed-cap arithmetic is as follows.

* If $M=1$, take $(a,S)=(J_1,Y_0)$.
* If $M=2$, let the two norms be $a\le b$. Their sum is the exact integer $J_1=a+b<N$, since $A<p$. Their gap is $g=\sqrt{(2J_2-J_1^2)\bmod N}$. Indeed the radicand is congruent to $(b-a)^2$ and $0\le (b-a)^2<A^2<p^2<N$. Thus $a=(J_1-g)/2$. If $g=0$, then $S=Y_0$. Otherwise $g<p$ is a unit and $S=(bY_0-Y_1)/g\bmod N$.
* If $M\ge3$, only one ray occurs. Set $U=M(M+1)(2M+1)/6$. Its three numerator factors are all at most $2M+1\le M^2\le A<p$, so $U$ is invertible modulo $N$. Then $a=J_1/U\bmod N$ and $S=Y_0/U\bmod N$.

The divisions, square roots and gcds have polynomial bit cost. The norm minimum and individual points are never input to the decoder.

## Converse with supplied factors

For promised distinct odd primes $p<q\le2p$ given in binary, the [factor-supplied evaluator](weighted_recurrence_control.py) computes all five residues in polynomial bit time. Its input includes the two primes. Primality testing is unnecessary under this promise. Construct one mixed root by CRT and its index-$N$ lattice in the integral $(x,y)$ coordinates of the public form. The CRT basis has $O(\log N)$-bit entries. Track the physical second coordinate of each basis vector while Gauss-reducing its divided norm to $[a,b,d]$ with $|b|\le a\le d$. For a reduced-coordinate vector $(u,v)$, completing the square bounds $v^2\le 4aA/\Delta\le4/3$ on the disk. The $v=0$ line gives one closed-form ray. For $v\ne0$ choose the simultaneous sign giving $v=1$. The least norm on this row is $d$, attained at $u=0$. If $d\le A$, then $4aA\ge4ad\ge\Delta\ge3A^2$ gives $a\ge3A/4$. At $|u|\ge2$ the norm is at least $a(u^2-|u|)+d\ge3a>A$. It suffices to test $u=-1,0,1$. The change of basis gives the physical $y$ for each vector. Add the principal terms. Two-dimensional Gauss reduction takes a polynomial number of bit operations on these inputs. The remaining bounds, fixed-degree sums, arithmetic and tests use $O(\log N)$-bit integers. This converse uses the factors exactly at the construction of the mixed line. It proves a polynomial-time equivalence for the promised family; it supplies no input-only method for that step.

## Boundary of a linear-recurrence shortcut

Let $s(B)$ be the principal-subtracted $y^2$ sum modulo $N$ at each integer cutoff. Suppose $s(B)=u^T T^B v$ on $0\le B\le A$ for a fixed $d$-dimensional matrix over $\mathbb Z/N\mathbb Z$, allowed to depend on the input. By Cayley–Hamilton, the first $d$ zero terms force all subsequent terms of this representation to vanish. On the other hand $s(0)=\cdots=s(a-1)=0$ and $s(a)\ne0$. Therefore $d>a$. This also covers scalar homogeneous constant-coefficient recurrences of order $d$. The repository's [early-witness counting bound](../ellipse_disk_mod8/bisector.md) gives $a>N^{1/6-\epsilon}$ for almost all pairs in the specified dyadic balanced-prime family, for any fixed $0<\epsilon<1/6$. The obstruction applies to a common linear state producing *all* integer prefixes after public subtraction. It leaves room for algorithms that compute only the fixed cap above, change state by cutoff, or use nonlinear operations.

## Verification

From the repository root run

```sh
uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_aggregate/weighted_recurrence_checks.py
```

The [checks](weighted_recurrence_checks.py) compare all five residues against an independent integer-square scan of each shell, validate $(a,S,M)$, and verify each returned divisor. The 560 exhaustive N-only cases with $q\le241$ included 157 single-point, 194 unequal-two-norm, 10 tied, and 199 three-or-more-point cases. All 194 historical public branches and 10 additional tied controls matched. The historical $\ell=11$ collision cases $N=11009,28829$ gave $(a,S)=(1,1849),(1,1444)$ respectively. On $N=38{,}942{,}974{,}393$ the single full disk at $A=358$ gave $a=205$, $S=72{,}233{,}001$, and factors $153{,}247$ and $254{,}119$. It used 7,467,686 long-coordinate probes and 12,042 short-coordinate insertions for 19 hits. These costs belong to the evaluator; the constant query count only bounds the decoder's calls.

The factor-supplied method matched the independent shell enumerator and public join for all 764 small branches and controls. We also generated a separate 241-bit balanced semiprime from two primes verified using FLINT. The factor-supplied evaluator took 33 Gauss steps at $A=884415271641753766$. The five residues decoded the first norm $134202694750453207$ and recovered a supplied prime by the verified gcd. Only the reverse direction was exercised at this size. [Full check output](weighted_recurrence_results.json)
