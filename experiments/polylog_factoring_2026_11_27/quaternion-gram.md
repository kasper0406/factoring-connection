# A binary Gram matrix embedded in four squares

Let $L=\left(\begin{smallmatrix} A&B\\B&C\end{smallmatrix}\right)$ be
positive definite with $AC-B^2=N^2$ for a positive odd $N$. There is an
unconditional expected-polynomial-time algorithm that returns a primitive
ordered pair $U,V\in\mathbf Z^4$ with Gram matrix $L$. When
$\gcd(A,B,C)=1$, it samples **all** such primitive ordered pairs exactly
uniformly. This includes cases where both diagonal entries share factors
with $N$. The algorithm uses a primitive point on the sphere of radius
$N$ and one Euclidean quaternion gcd. [Implementation](quaternion_gram.py) ·
[Primitive sphere sampler](ternary-sphere.md)

## Construction and uniformity

At least one diagonal entry is odd, since $N^2+B^2$ has 2-adic valuation
at most one. Swap the columns if needed so $A$ is odd. Draw a uniform
primitive integer point $(x,y,z)$ on $x^2+y^2+z^2=N^2$ and put
$$
T=B+xi+yj+zk,\qquad \operatorname{Norm}(T)=AC .
$$
In the Hurwitz order $\mathcal H$, compute a generator of the right ideal
$I=T\mathcal H+A\mathcal H$ by a greatest common left divisor. Its norm
is $A$. Indeed at an odd $r\mid A$ the local order is
$M_2(\mathbf Z_r)$. Primitivity of $(x,y,z)$ makes $T$ nonzero modulo
$r$, so its local Smith invariants have valuations $0$ and $v_r(AC)$.
Adding $A\mathbf Z_r^2$ to the column lattice of $T$ gives index
$r^{v_r(A)}$. There are two columns, hence the local ideal has index
$r^{2v_r(A)}$. Elsewhere $A$ is a unit. The full index $A^2$ equals the
square of the norm of any generator.

Choose an integral generator $\alpha$ uniformly among its eight integral
right-unit associates. An odd-norm Hurwitz generator always has eight
such associates, as in the four-square sampler. Put
$$
U=\bar\alpha,\qquad V=\bar\alpha T/A .
$$
Since $T\in\alpha\mathcal H$, $V$ is Hurwitz integral. Its coordinate
denominators also divide the odd integer $A$, so they are one. Now
$\bar U V=T$, which gives $U\cdot V=B$ as well as the prescribed norms.
The imaginary coordinates of $\bar U V$ are signed sums of the $2$-by-$2$
minors of $(U,V)$. Any common divisor of all the minors would divide
$x,y,z$, so their gcd is one. [Integral associates](four-square-incidence.md)

The eight choices exhaust the fiber at every $T$. In fact, for any other
pair with $\bar U V=T$, the right ideal $\bar U\mathcal H$ contains $T$
and $A=\bar UU$. Its index is again $A^2$, giving equality with $I$.
Thus the attainable subset consists exactly of the pairs with primitive
$\operatorname{Im}(\bar U V)$. It has $8$ times the number of primitive
sphere points, and all have equal mass. Swapping the columns negates the
imaginary part and preserves that law.

To see why this covers *every* primitive pair for primitive $L$, suppose
$r$ divides all three imaginary coordinates. Their squares sum to $N^2$,
so $r\mid N$ is odd. At least one of $A,C$ is a unit modulo $r$, since
$r\mid A,C$ would force $r\mid B$. When $A$ is a unit, the congruence
$\bar U V\equiv B$ implies $V\equiv (B/A)U\bmod r$. With $C$ a unit
one gets the reverse dependence by conjugation. Either makes every minor
divisible by $r$, which contradicts primitivity. For imprimitive $L$
the attainable subset can be proper. The matrix $3I_2$ has 384 primitive
embeddings, exactly 192 of which have primitive imaginary part.

## Earlier results and verification

Mordell studied quaternion constructions for binary forms as four linear
squares in 1937. Hunsucker’s Theorem 1 gives the primitive existence
criterion. For $L=eL_0$ with primitive Gram $L_0$ of determinant $d$, it
is $e^2d\equiv3\bmod8$, or $d\equiv1,2\bmod4$ with $4\nmid e$.
An odd square determinant satisfies the latter condition. Hunsucker’s
Theorems 2 and 3 count primitive embeddings using genera and local
congruences. The proof invokes Pall’s completion algorithm. These papers
establish the prior representation theory. The sampler above specifies
how to draw its sphere witness without knowing the factors and proves a
bit-complexity and uniformity guarantee. [Mordell, *An Application of Quaternions to the Representation of a Binary Quadratic Form as a Sum of Four Linear Squares*](https://doi.org/10.1093/qmath/os-8.1.58) ·
[Pall, §3](https://doi.org/10.4153/CJM-1949-032-9) ·
[Hunsucker, pp. 321–325](https://doi.org/10.4064/aa-19-4-321-325)

Independent enumeration at nine Gram matrices recovered exactly 2,448
supported frames out of 3,024 primitive frames across the fixtures. They
include $\mathrm{diag}(9,25)$ at $N=15$, an even first diagonal entry,
negative inner product, imprimitive Grams, and $N=1$. The public sampler
also passed 88 seeded calls, including determinant $(1009\cdot1013)^2$.
[Reproduce the checks](check_quaternion_gram.py)
