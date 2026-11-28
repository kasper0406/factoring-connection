# A ternary Gram with 48 factor-revealing representations

For every product $N=pq$ of distinct odd primes, there is a public randomized construction of a binary Gram matrix with exactly 72 marked representations in $\mathbb Z^3$. In 48 of them, the gcd of the second vector is the same proper prime factor of $N$. The remaining 24 are signed coordinate permutations of a representation available directly from $N$. Consequently, **any representation outside a coordinate plane returns a factor**, and one uniform representation succeeds with probability $2/3$. The matrix construction takes expected polynomial bit time. The open step is obtaining one of those 48 representations at comparable cost. [Executable construction and checks](check_mixed_ternary_embeddings.py)

## Constructing the matrix from the input

Draw $x$ uniformly modulo $N$ until the Jacobi symbol $\left(\frac{-(x^2+1)}N\right)=-1$. Put

$$
A=x^2+1,\qquad B=-Nx,\qquad C=N^2.
$$

The Gram determinant $AC-B^2$ is $N^2$. The Jacobi test ensures $A$ is a unit modulo $N$ and exactly one of $p,q$ admits a square root of $-A$. The routine uses the Jacobi symbol directly. A zero symbol also permits an immediate gcd test.

The cost bound needs no assumption about primes in progressions. For an odd prime $s$ let $\chi_s$ be its quadratic character and $h_s=\chi_s(-1)$. The polynomial $x^2+1$ has $1+h_s$ zeros, and $\sum_{x\bmod s}\chi_s(x^2+1)=-1$. Indeed $(y-x)(y+x)=1$ has $s-1$ solutions. By CRT the number of accepted residues is exactly

$$
\frac{(p-1-h_p)(q-1-h_q)-h_p h_q}{2}.
$$

The acceptance probability exceeds $1/8$ for every pair in scope and tends to $1/2$ as both primes grow. Each local nonzero fraction is at least $3/5$ and $pq\ge15$. This works for both primes $1\bmod4$, both $3\bmod4$, and opposite residue classes. For $N\equiv3\bmod4$ the deterministic choice $x=0$ already passes.

## Why there are 72 representations

A marked representation consists of ordered signed vectors $u,v$ with $u\cdot u=A$, $u\cdot v=B$ and $v\cdot v=C$. Let $d$ be the positive gcd of their $2\times2$ minors. Their squared minors sum to $N^2$, so $d\mid N$. For each coordinate $j$,

$$
A v_j-Bu_j=\sum_i u_i(u_i v_j-u_j v_i).
$$

Since $d\mid B$ and $A$ is a unit modulo $d$, this gives $d\mid v_j$. Conversely $\gcd(v)$ divides every minor. Hence **$d=\gcd(v)$**. These identities hold in any ambient dimension.

Set $m=N/d$. Dividing $v$ by $d$ gives a primitive marked representation of the smaller Gram $(A,-mx,m^2)$, where primitive means the minors have gcd one. Multiplying the second vector of any primitive representation of that matrix by $d$ gives the inverse. In three dimensions the primitive representations of the smaller Gram are counted by

$$
24\,\#\{s\bmod m:s^2\equiv-A\pmod m\}.
$$

For completeness, the oriented normal of a primitive pair has squared length $m^2$. Adjoin that formal normal to the input plane. The plane's discriminant group is cyclic of order $m^2$, with generator $(-B_m e_1+A e_2)/m^2$ of squared norm $A/m^2$, where $B_m=-mx$. The normal line's dual quotient has generator of squared norm $1/m^2$. The primitive orthogonal inclusions in a unimodular ternary lattice glue these two groups by an isomorphism. The possible gluings are exactly the scalars $t$ with $t^2\equiv-A\pmod{m^2}$. Each unit root modulo odd $m$ lifts uniquely, and each resulting positive unimodular ternary lattice has an oriented orthonormal basis. An orientation-preserving isometry of marked frames fixes both the input generators and the oriented normal, so two different scalars give different orbits. Each orbit has 24 members because a proper signed permutation fixing two independent vectors is the identity. At $m=1$ there is one gluing and the same count of 24. This specializes Pall's representation count and his codimension-one construction. [Prescribed-root construction](prescribed-ternary-embedding.md) · [Pall, §3 and Theorem 4, p. 359](https://doi.org/10.4153/CJM-1949-032-9)

Write $r$ for the unique prime factor of $N$ where $-A$ is a square. The only contributing divisors are $m=1$, with one root, and $m=r$, with two. The former gives 24 representations with $d=N$, the latter 48 with $d=N/r$. A coordinate plane contains both vectors exactly when $d=N$, since the primitive normal then has length one. For any sampler within total variation $\eta$ of uniform, the success probability is at least $2/3-\eta$. Sampling uniformly among the 24 public coordinate frames attains distance exactly $2/3$ and yields zero proper factors. The public coordinate example is $u=(x,1,0)$ and $v=(-N,0,0)$. On $N=77$, $x=1$ passes although both primes are $3\bmod4$. The representation $u=(1,1,0)$, $v=(-14,-63,42)$ has inner products $(2,-77,5929)$ and immediately yields $\gcd(v)=7$.

## Relation to earlier circle and sampling tests

The earlier circle-prefix calculation used inputs $a^2+1=pq$ and identified the same hidden orbit as the original predictor. The shifted-circle calculation bounded gcd tests on one supplied point of squared norm $KN^2+k$ with small nonzero shift. Here the semiprime is arbitrary, and the public Jacobi condition forces opposite local signs for a *pair* of vectors with three prescribed inner products. At $x=0$ the noncoordinate case reduces to a non-axis sum of two squares of $N^2$, recovering the familiar Pythagorean-triple special case. Historical two-squares factoring also uses inequivalent norm representations. The count of representations of a binary form by three squares goes back to Pall. The input construction and the distinction between a freely chosen representation and a factor-revealing one are the focus here. [Earlier predictor and circle ledger](README.md) · [Shifted-circle scope](shifted-circle-gcd.md) · [Euler's factoring method and its generalizations, Riesel, Chapter 5](https://doi.org/10.1007/978-1-4612-0251-6_5)

A uniform sampler of **all** four-dimensional representations is available using public data. On this same Gram its probability of returning a proper gcd is $(\lambda(p)+\lambda(q))/((1+\lambda(p))(1+\lambda(q)))$, where $\lambda(s)=s-\chi_s(-1)$. This is of order $1/p+1/q$. At $N=77$ it is $20/117$. The ternary count reaches $2/3$; knowing the number 72 alone gives no indication which prime divides the vector. Scalar three-square algorithms prescribe one norm, whereas these 72 ordered pairs must satisfy both lengths and their inner product. [Four-dimensional construction](quaternary-embedding.md) · [Pollack and Schorn on the scalar problem](https://doi.org/10.1090/mcom/3349)

The check enumerates all signed triples $v$ of squared length $N^2$, then all possible $u$ in the coordinate intervals implied by $\|u+(x/N)v\|=1$. Across 39 challenges and controls on eight semiprimes it verified the content of all 3,576 frames, including 23 mixed challenges and eight controls for each equal-sign case. This enumeration scales with $N$. A separate 234-bit control used a **supplied prime factor** to construct 16 noncoordinate frames and checked the public gcd on each. [Reproduce both checks](check_mixed_ternary_embeddings.py)
