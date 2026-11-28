# Gaussian residues for perpendicular shapes on the integer sphere

Let $N>1$ be odd. A uniformly chosen primitive signed, ordered point $v\in\mathbb Z^3$ with $v\cdot v=N^2$ gives a **uniform oriented perpendicular form in one genus** of discriminant $-4N^2$. Explicit Gaussian and quaternion Euclidean algorithms map a residue of norm $-1\bmod N$ to such a point and recover its residue up to a global sign. They work with repeated prime factors and use only the input radius. [Implementation](ternary_embedding.py) · [Uniform source of points](ternary-sphere.md) · [Prescribed-class inversion and factoring](prescribed-ternary-embedding.md)

The squaring of the Picard action on perpendicular lattices is classical. Ellenberg, Michel and Venkatesh describe the square action and its coset of shapes. Khayutin's Remark 8.15 gives the map $\lambda\mapsto\lambda/\bar\lambda$, including the quotient at the prime above two when two ramifies. The specialization below identifies the actual residue returned by the existing Hopf sampler, gives an inverse from a marked point, and accounts for the integral rotations at this Gaussian conductor. [Ellenberg, Michel and Venkatesh, *Linnik's ergodic method and the distribution of integer points on spheres*](https://www.epfl.ch/labs/tan/wp-content/uploads/2018/10/pointsonspheres.pdf#page=24) · [Khayutin, Remark 8.15, pp. 3001–3002](https://ems.press/content/serial-article-files/32576?nt=1#page=54)

## Computing the coordinate

Write $O=\mathbb Z[i]$ and $h=a+bj$ in Hamilton's quaternions, with $a,b\in O$ and $|a|^2+|b|^2=N$. Then $v=h i\bar h$ has longitudinal coordinate $|a|^2-|b|^2$ and transverse complex coordinate $-2iab$, in the oriented basis $i,j,ij$. The existing function $\mathrm{hopf\_point}$ uses the Gaussian pair $a,i\bar b$, or integer tuple $(\Re a,\Im a,\Im b,\Re b)$. To put an arbitrary sphere point in this convention, cyclically rotate its unique odd coordinate into first position.

For a primitive point, $\gcd_O(a,\bar b)$ is a Gaussian unit. Indeed a Gaussian prime common to both has underlying rational prime dividing both norms and the two coordinates of $ab$. Conversely, if an odd prime divides the three coordinates of $v$, it divides $N$, both norms, and both coordinates of $ab$. Locally at an inert prime this forces both $a,b$ to vanish. At a split prime, examining the two conjugate components shows that $a,\bar b$ vanish together in at least one. The prime above two cannot divide both when $N$ is odd.

A Gaussian Bézout identity therefore gives the unique $t\bmod NO$ satisfying
$$
 ta=b,\qquad t\bar b=-\bar a .
$$
The compatibility is $|a|^2+|b|^2=0\bmod N$. Conjugating the equations gives $(t\bar t+1)a=(t\bar t+1)\bar b=0$, so $t\bar t=-1\bmod N$. From a point $v=(x,y,z)$ in the chosen orientation, write $y/2+iz/2=U\bar V$ with $U=a,\ V=i\bar b$. A Gaussian gcd of $(N+x)/2=|U|^2$ and $U\bar V$ recovers $U$ up to a unit, and division recovers $V$. The four common unit choices for $(U,V)$ change $t$ only by a global sign.

Conversely, given $t\bar t=-1$, the pairs $b=ta\bmod N$ define a left quaternion ideal modulo $N$, generated there by $1+tj$. It has index $N^2$: $a$ is arbitrary and $b$ is determined in the quotient. The Hurwitz and integral-coordinate lattices have the same local order at the odd conductor. Lift the left ideal by adjoining $N$, conjugate, and use the Hurwitz Euclidean gcd. The generator has norm $N$. Of its 24 unit associates exactly eight have integral coordinates. For such a generator, $r\mid a,b$ at a rational prime $r\mid N$ is impossible: after division by $r$ the congruence still holds modulo $N/r$, so the new norm would be divisible by $N/r$, whereas it equals $N/r^2$. Since $t$ is a unit, any Gaussian prime dividing $a,\bar b$ would imply this rational-prime divisibility. Thus all eight Hopf outputs are primitive. [Hurwitz gcd and integral associates](four-square-incidence.md)

## Identifying the oriented form

For a Gaussian residue $t$, put $J_t=\mathbb Z t+N O$, using any lift. If $t$ is a unit modulo $N$, $J_t$ has index $N$ in $O$. The oriented isometry $F(\alpha)=h(\alpha j)\bar h/N$ from the Gaussian plane to $v^\perp$ sends $1,i$ to the rotations of $j,ij$, whose cross product points along $v$. The numerator is
$$
 2i\,\Im(a\bar b\alpha)+(a^2\alpha+b^2\bar\alpha)j .
$$
For $\alpha=t$ both components vanish modulo $N$ by the two residue equations. Thus $F(J_t)\subseteq v^\perp\cap\mathbb Z^3$. Both lattices have covolume $N$, the latter because $v$ is primitive and $\|v\|=N$. Equality gives the norm form on $J_t$ exactly with this orientation.

Let $R=\mathbb Z/N$, $P=(O/NO)^\times/R^\times$, and $G=P/\langle[i]\rangle=\operatorname{Pic}(\mathbb Z+NO)$, using the Gaussian maximal order's class number one. The class of $J_t$ is given by the projective residue of $t$. In particular, two such lattices of the same index are properly isometric exactly when their generating residues differ by an invertible real scalar and a global Gaussian unit. Indeed a rational-plane orientation-preserving isometry is multiplication by a norm-one Gaussian rational. Extending each $J_t$ to $O$ gives $O$, so that multiplier must preserve $O$ and is a Gaussian unit. For unoriented isometry one may also conjugate.

Fix any public $t_0$ of norm $-1$. The possible residues are $t_0 w$ with $w\bar w=1$. The map $u\mapsto u/\bar u$ identifies $P$ with the norm-one residues. Locally at each odd prime dividing $N$, one can solve $u=w\bar u$ by $u=\lambda+w\bar\lambda$ with $u$ a unit: modulo the prime this is a nonzero rank-one linear map onto a line of units apart from zero, and a unit solution lifts with the same formula. CRT combines the choices, and the kernel consists of the real unit scalars. Projectively $[u/\bar u]=[u]^2$. The occupied oriented classes are therefore exactly $[t_0]G^2$.

## Uniformity, rotations and checks

Each norm-minus-one residue has eight equally weighted quaternion generators. Each marked sphere point with its first coordinate odd has four Hopf lifts: any other differs by a right multiplier $c\in\mathbb Q(i)$ of norm one, and integrality with $\gcd(a,\bar b)=1$ forces $c\in O^\times$. Cyclic rotation gives three choices of odd axis. The existing sampler consequently gives the uniform law on the coset, with six times the number of residue representatives per oriented form.

With $\omega$ distinct primes dividing $N$, each oriented class has $2^{\omega+1}$ residues of norm $-1$. Indeed their ratios are $i^j\lambda$, with $j=0,1$ and $\lambda^2=1\bmod N$. There are $2^\omega$ independent CRT sign choices for odd $N$. An oriented class therefore has $6\cdot 2^{\omega+1}$ signed ordered points. The 24 proper integral rotations act freely for $N>1$: a nontrivial cube rotation fixes only a coordinate, edge or vertex axis, whose primitive squared length is $1,2$ or $3$. This leaves $2^{\omega-1}$ rotation orbits per occupied oriented class. For a semiprime with distinct odd factors there are 48 points and two rotation orbits per class. Independent unconditioned draws have probability $1/|G^2|$ of sharing an oriented class. The prescribed-class problem asks for a point after the class has already been specified.

From the repository root run

~~~sh
python experiments/polylog_factoring_2026_11_27/check_ternary_embedding.py
~~~

The check compared the positive-determinant HNF norm form of $J_t$ to the oriented perpendicular form for 329 residues and all 1,974 primitive points on 14 odd radii through 65. It checked 497 matched scalar roots, testing both outcomes against ordinary form reduction, and ran 20 seeded quaternion checks including $1009\cdot1013$ and $3^{23}5^{17}$. At radius 15 the two points $(5,-14,2)$ and $(11,-2,-10)$ give the same oriented form. For $z=(2,5)$ the matcher returns $1$ and $11\bmod15$; the gcd of their difference with 15 is 5. [Reproduce the checks](check_ternary_embedding.py)
