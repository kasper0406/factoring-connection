# Computing parity on one packed sphere class

Two elementary computations answer a single packed parity query directly from its stated inputs $(M,m,c)$. A sign or coordinate symmetry gives the answer in polynomial bit time whenever a residue coordinate is zero or two coordinates agree up to sign. In every remaining class, the near-square packing reduces the reference search from two nested coordinate loops to one. On six earlier semiprimes the symmetry answered 740 of 11,858 queries, and the two methods together used 18,666 candidate square tests against 116,555 in the original scalar enumeration. Each answer was compared with that enumeration. [Packed reduction](packed-ternary-decision.md) · [Native evaluator](native_ternary_parity.py) · [Exact comparison](check_native_ternary_parity.py)

Here $m=Q\ell$, $Q$ is a power of seven and $7\nmid\ell$. The input $M=(QN+x)^2+1$ recovers $N,x$ by integer square root and division by $Q$. The centered residues of $c\bmod Q$ give $u$, and $z=(u-c)/Q\bmod\ell$ gives the residue of the partner $v$. The evaluator checks $\|u\|^2=A=x^2+1$, $\gcd(A,N)=1$ and the packing size bound. Each solution $W$ to the scalar query corresponds to exactly one $v$ with $\|v\|^2=N^2,\ u\cdot v=-Nx,\ v=z\bmod\ell$. The implementation falls back to enumeration for a query outside its recognized input family.

## Classes fixed by a mirror

For the unary series $\theta_j(T)=\sum_{k=j\bmod m}T^{k^2}$ at odd $m$, one has $\theta_{-j}=\theta_j$, $\theta_0=1\bmod2$, and $\theta_j(T)^2=\theta_j(T^2)\bmod2$. On the packed three-square query the corresponding sign flip or signed transposition preserves $c\bmod m$. It fixes $u$ over the integers, since $Q>2\max_i|u_i|$.

Consequently the mirror permutes the partners $v$ in the queried class. An integer point fixed by the mirror would have both $u$ and $v$ in its reflection plane. The primitive normal of such a plane has length $1$ or $\sqrt2$. The latter is impossible because $\|u\times v\|=N$. In the former case the gcd of the minors is $N$, and the unit condition on $A$ implies $N\mid v_i$ for all $i$. Thus every fixed point comes from an axis partner $v=\pm N e_i$. Conversely the set of possible axis partners is invariant under the mirror. Its parity, computed by testing these six vectors, gives the full answer.

This also identifies exactly which parity queries escape these symmetries. If all $c_i$ are nonzero and distinct up to sign modulo $m$, the signed-permutation stabilizer of $c$ is trivial. Every coordinate frame has a zero residue in the unused coordinate. Thus a generic class of odd parity contains a noncoordinate frame. For the mixed semiprime promise each such frame supplies a proper gcd. Parity alone does not provide its coordinates.

## One-dimensional count for a generic class

Choose $j$ minimizing $D=A-u_j^2>0$ and write $k,l$ for the remaining indices. At a candidate $v_j=t$, set $h=-Nx-u_j t$ and $w=u_kv_l-u_lv_k$. Then
$$
 w^2=D(N^2-t^2)-h^2,\qquad
 (At+Nx u_j)^2+A w^2=N^2D.
$$
The second identity restricts $At+Nx u_j$ to the interval $[-\lfloor N\sqrt D\rfloor,\lfloor N\sqrt D\rfloor]$. Enumerate only $t=z_j\bmod\ell$ there. When the first right-hand side is a square, test each distinct sign of $w$ by
$$
 v_k=(u_kh-u_lw)/D,\qquad v_l=(u_lh+u_kw)/D.
$$
Integer divisibility and the two remaining congruences determine exactly which partners count. This uses $O(1+N\sqrt D/(A\ell))$ candidates, with polynomial bit work per candidate. For $u=(x,1,0)$ and $x\ge1$ the bound is $O(1+N/(A\ell))$. For fixed $x$ and $\ell$ of logarithmic size the worst-case bound still grows as a power of $N$. The binary CM formula available at $x=1$ evaluates a sum over all partners; an individual angular residue generally distinguishes the summands. [Fixed-first-vector binary norm](short-first-vector.md) · [Quartic traces and the CM coefficient](ternary-frame-moments.md)

## Exact checks

| $N$ | $x$ | Mirror queries | Odd mirror | Generic queries | Odd generic | Native square tests | Scalar square tests |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 15 | 0 | 104 | 72 | 728 | 144 | 468 | 536 |
| 35 | 0 | 136 | 92 | 1,112 | 184 | 1,360 | 2,660 |
| 65 | 4 | 74 | 25 | 3,442 | 50 | 682 | 14,171 |
| 77 | 1 | 78 | 52 | 1,828 | 104 | 1,896 | 9,316 |
| 91 | 0 | 164 | 108 | 1,840 | 208 | 4,336 | 15,628 |
| 183 | 0 | 184 | 124 | 2,168 | 232 | 9,924 | 74,244 |

The reference counted all signed ordered sphere points in each prescribed class. All six end-to-end parity runs recovered a proper factor. At $N=15,x=1,\ell=1,u=(0,1,1)$ an additional mirror check counts six points and returns even, illustrating the difference from an existence test. The comparison also covers 256 ordinary sphere queries outside the packing. Run `python experiments/polylog_factoring_2026_11_27/check_native_ternary_parity.py` from the repository root.
