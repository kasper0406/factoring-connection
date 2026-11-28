# Prescribing the perpendicular lattice on a sphere

Fix an odd $N>1$. An integer point on the sphere of radius $N$ is primitive if its three coordinates have gcd one. We can sample such points exactly and uniformly using $N$ alone. Prescribing the shape of the full lattice of integer vectors perpendicular to the point changes the computational problem. An algorithm returning a matching point on a fraction $\epsilon$ of the **uniformly sampled attainable oriented shapes at this same $N$** gives a modular square-root algorithm with success at least $\epsilon$ on every unit square modulo $N$. On a semiprime with distinct odd factors, one call then gives a factor with probability at least $\epsilon/2$. [Public sphere sampler](ternary-sphere.md) · [Shape correspondence](ternary-shapes.md) · [Reduction](ternary_form_reduction.py)

## What a primitive embedding determines

Use Gram entries $(A,B,C)$ for the binary form $AX^2+2BXY+CY^2$, where $A>0$, $AC-B^2=D=N^2$ and $\gcd(A,N)=1$. A marked embedding gives $u,v\in\mathbb Z^3$ with inner products $A,B,C$. Require $w=u\times v$ primitive. In that case $u,v$ span the full perpendicular lattice of $w$, which has squared length $D$.

Choose $l\in\mathbb Z^3$ with $l\cdot w=1$ by Bézout. Since $(w\times u)\times w=Du$, reduction modulo $D$ gives
$$
w\times u\equiv \lambda w\pmod D,\qquad \lambda=l\cdot(w\times u).
$$
Indeed any vector whose cross product with a unimodular vector vanishes is a scalar multiple of it, as follows by contracting the vanishing minors with $l$. Taking the cross product with $u$ once more gives $-A w\equiv\lambda^2w\pmod D$. Thus $\lambda^2\equiv-A\pmod{N^2}$. The computation uses the supplied vectors only.

Conversely, a root $s^2\equiv-A\pmod N$ has a unique lift modulo $D$ by the unit derivative $2s$. Start with an abstract copy of the input lattice on $e_1,e_2$ and an orthogonal vector $f$ of norm $D$, where norm means squared length. Adjoin
$$
(-Be_1+Ae_2+s f)/D .
$$
Its norm $(A+s^2)/D$ and all pairings are integral. Its exact order modulo the original lattice is $D$, so the enlarged rank-three lattice has determinant one. Its intersection with the input plane is exactly the original binary lattice, by the unit $s$. Integral Gram LLL with parameter $15/16$ gives an orthonormal basis in polynomial bit time. In ranks at most four, the first squared length is at most $(16/11)^{(d-1)/2}<2$, hence equals one, and size reduction makes the other vectors perpendicular to it. Apply the same argument on the remainder. In an oriented basis $f=u\times v$, and the root read above is $-s$. [Implementation](gram_embedding.py) · [Prescribed-sign and 512-bit checks](check_gram_embedding.py)

In particular **testing existence** of a primitive embedding here is the test whether $-A$ is a square modulo $N$. To encode an arbitrary unit $a$, find a Gaussian integer $z$ with $|z|^2\equiv-a\pmod N$ by the public Gaussian norm algorithm. Set $h=\gcd(\Re z,\Im z)$ and $e_1=z/h$ as an integral plane vector. Complete it to an oriented unimodular pair $e_1,t$, and take the Gram of $e_1,Nt$. It has determinant $N^2$, and $h$ times a root for $-A$ gives a root of $a$. This supplies both directions for decision and for search. A decision procedure alone has only been reduced to quadratic residuosity here. The factoring reduction uses a witness. [Norm algorithm](four_square_sampler.py) · [Construction and checks](check_ternary_form_reduction.py)

## Only the shape goes to the oracle

Write $O=\mathbb Z[i]$ and $J_z=\mathbb Z z+NO$. Its oriented shape can be supplied as a Gauss-reduced binary quadratic form. Represent a unit square $a$ by a uniform Gaussian residue $z$ of norm $-a$ modulo $N$. To make an initial solution uniform, multiply by $w/\bar w$ for a uniform Gaussian unit $w$ modulo $N$. By the local norm-one parametrization, this ranges uniformly over the entire fiber. Multiplying a residue by a unit scalar modulo $N$ preserves $J_z$. Since $a$ is a square, this shape distribution equals the one for norm $-1$. It is exactly uniform on the attained oriented classes. The query contains only $N$ and the reduced form. [Proof of the distribution and public maps](ternary-shapes.md)

From a returned point compute its Gaussian residue $t$ of norm $-1$ by the inverse Hopf map. The perpendicular lattice has the shape of $J_t$. Isometric Gaussian line lattices differ by a global Gaussian unit, allowing conjugation when orientation is forgotten. Indeed an oriented rational plane isometry is multiplication by a norm-one Gaussian rational, and both lattices generate all of $O$ over $O$. Thus, after trying the four global units and conjugation, $z=c t\pmod{N}$ with $c$ a scalar unit. Taking norms gives $c^2=a\pmod N$. Finding the root requires no preselected point or matching basis. [Matcher](ternary_embedding.py)

To factor, take a uniform secret unit $s$ and request a root of $a=s^2$. Conditioned on $a$, the returned root is independent of which of the $2^{\omega(N)}$ roots was chosen as $s$. Comparing it with $s$ by gcd gives a proper divisor with probability $1-2^{1-\omega(N)}$ on a successful call. An oracle may always return the same representative of each shape. Its chance of returning some correct point is the only oracle assumption.

## Prior work and checks

Pall's codimension-one completion uses a row satisfying a congruence modulo the Gram determinant. His equations (19)–(21) specialize to the gluing above. Khayutin describes the squared class-group action on perpendicular shapes. Bosma and Stevenhagen give algorithms for generating the principal genus and taking class-group square roots with knowledge of the discriminant factorization. Here the executable calculation crosses from a prescribed form at square determinant to ordinary scalar roots, including queries distributed uniformly among attainable shapes. [Pall, §3 (1949)](https://doi.org/10.4153/CJM-1949-032-9) · [Khayutin, §8](https://ems.press/content/serial-article-files/32576) · [Bosma and Stevenhagen](https://jtnb.centre-mersenne.org/articles/10.5802/jtnb.170/)

The reference oracle used one representative from an exhaustive small sphere table. On every challenge over $N=15,35,45$ the exact split fraction, counting every secret root, was $1/2$. At $N=105$ it was $4608/6144=3/4$. In each case the challenge forms were exactly uniform among the attained oriented classes. The separate marked test passed 256 public input-and-root checks. The table oracle uses enumeration. Directly sampling a four-dimensional embedding and deleting a coordinate gives a second, entirely public test with a much lower success rate. [Reproduce the counts](check_ternary_form_reduction.py) · [Four-dimensional construction and exact rate](quaternary-embedding.md)
