# Rotating one point on the packed sphere

Let $N=pq$ be a product of distinct odd primes and take $A=x^2+1$ with Jacobi symbol $(-A/N)=-1$. Use the packing scale $Q$ from the scalar decision reduction. Write $y=x+QN$, $M=y^2+1$, and $W_0=(y,1,0)$. Suppose $W=u-Qv$ is one of the 48 factor-revealing packed points, where the pair $u,v$ has Gram $(A,-Nx,N^2)$. Write $r=N/\gcd(v)$, the prime where $-A$ is a square.

Exactly one rational proper rotation from $W_0$ to $W$ has denominator below $M/(2r)$. Its reduced denominator is $r$. Other rotations to this **same single point** can have denominator coprime to $N$. One always exists with denominator at most $M/r+r$. For a supplied target, every possible reduced denominator follows from a single reduced binary quadratic form of discriminant $-4M$. [Original frame rotations](ternary-frame-rotations.md) · [Packed sphere](packed-ternary-decision.md) · [Implementation](packed_point_rotations.py)

## The binary form giving every denominator

Put $u_0=(x,1,0)$, $v_0=(-N,0,0)$, so $W_0=u_0-Qv_0$. There is a unique frame rotation $R_f$ taking both $u_0$ to $u$ and $v_0$ to $v$. Its columns are $-v/N$, $u+xv/N$, $(u\times v)/N$. Let $h_0$ be its primitive integral Euler quaternion. Its norm is $e r$, where $e\in\{1,2,4\}$. For any primitive integral Euler quaternion the reduced rotation denominator is the odd part of its norm. At an odd prime dividing the norm, the quaternion is a nonzero rank-one matrix and its conjugation numerator cannot vanish on all traceless matrices. At two, primitive quaternion norms have valuation at most two, and the corresponding power of two cancels from the Euler matrix. [Quaternion denominator proof](ternary-frame-rotations.md) · [Classical cubic coincidence rotations](https://journals.iucr.org/paper?a10612)

Identify three-vectors with pure imaginary quaternions. Since $h_0 W_0=W h_0$, the kernel line of $h_0\bmod r$ is $W_0$-stable. The action of $W_0$ on the quotient line is scalar. Hence there is exactly one residue $B\bmod r$ for which $h_0(B+W_0)\equiv0\pmod r$. Its square is $-M\bmod r$. Use the representative with $|B|<r/2$ and set

$$
h_1=h_0(B+W_0)/r,\qquad C=(B^2+M)/r.
$$

The rational quaternions intertwining $W_0$ and $W$ are the plane $h_0\mathbb Q(W_0)$, since the centralizer of $W_0$ is $\mathbb Q(W_0)$. Thus representatives of all its rotations are $h=c h_0+b h_1$ with coprime integers $c,b$, unique up to a common sign. Multiplicativity of the norm gives

$$
\operatorname{Norm}(h)=e F(c,b),\qquad
F(c,b)=r c^2+2Bcb+C b^2=\frac{(rc+Bb)^2+Mb^2}{r}.
$$

No odd prime divides every coordinate of $h$. Indeed $\bar h_0 h=e(rc+Bb+bW_0)$, whose $j$-coordinate is $eb$. A common odd prime must divide $b$, and then $c h_0$ modulo that prime forces it to divide $c$ as well. Primitive normalization can therefore divide only by a power of two. **The reduced denominator is exactly the odd part of $F(c,b)$.** The form $[r,2B,C]$ is primitive and reduced, since $r\nmid B$, $2|B|<r$ and $C\ge r$. The argument covers all three possible values of $e$.

The pair $(rc+Bb,b)$ is primitive at two. As $M=y^2+1$ equals 1 or 2 modulo 4, we have $4\nmid F(c,b)$. If $b=0$, the rotation is $R_f$ and the denominator is $r$. If $b\ne0$, the denominator is at least $M b^2/(2r)$. This explains why a single-point rotation of small denominator that stays in an informative packing class must already be the frame rotation.

For the upper bound, take $b=1$. Choose $c=0$ unless $r\mid C$, in which case use $c=-\operatorname{sgn}(B)$. The two candidates $B$ and $B-r\operatorname{sgn}(B)$ differ by $r$ and $2B$ is a unit modulo $r$, so the chosen numerator $(rc+B)^2+M$ is divisible by $r$ exactly once. At the other prime of $N$, $-M\equiv-A$ is a nonsquare and the numerator is a unit. Both candidate representatives have absolute value less than $r$. This proves the coprime denominator bound and constructs the quaternion using the **supplied target point**.

## Public small-prime neighbors

Given $W_0$ and a small odd prime $\ell$, enumerate primitive quaternions of norms $\ell,2\ell,4\ell$ up to sign. Keep each integral point $hW_0\bar h/\operatorname{Norm}(h)$. These are all integral images under rotations of reduced denominator $\ell$. This is a deterministic operation with work polynomial in $\ell$ and $\log M$. From $W_0$, any $\ell<M/(2r)$ distinct from $r$ yields zero informative packed points. This is a bound on individual prime steps and their total reduced denominator. A sequence may travel through points outside the packed classes and return with large denominator. The quaternionic class group and Hecke trajectories on integer spheres are classical in Linnik's method. [Ellenberg, Michel and Venkatesh, §§2 and 5](https://www.epfl.ch/labs/tan/wp-content/uploads/2018/10/pointsonspheres.pdf) · [Public neighbor routine](packed_point_rotations.py)

For $N=77,x=1,Q=2401$ one has $W_0=(184878,1,0)$, $M=34179874885$. The prime 19 gives 48 distinct integral neighbors, all outside the informative packed classes. The hidden good prime 11 gives exactly the 48 informative points. For the particular target $W=(33615,151264,-100842)$ the form is $[11,6,3107261354]$. Its $(0,1)$ rotation has prime denominator $1553630677=\lceil M/22\rceil$, attaining the lower bound after rounding. This verifies concretely that a direct transporter need not carry a factor of $N$ in its denominator.

The checks recovered the quadratic form at all 336 informative points on seven mixed inputs and compared 5,712 rotations with exact Euler matrices. The original frame-quaternion norm multipliers $1,2,4$ occurred 56, 168, 112 times respectively. Ninety-six coprime-denominator constructions required the second representative of $B$. [Reproduce the tests](check_packed_point_rotations.py)
