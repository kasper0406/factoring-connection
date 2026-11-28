# Sampling all primitive four-dimensional embeddings of a binary form

Let $L$ have positive integral Gram matrix $\left(\begin{smallmatrix}A&B\\B&C\end{smallmatrix}\right)$ with $AC-B^2=N^2$, $N>1$ odd, and $A$ coprime to $N$. There is an exact, factor-free sampler for **all ordered primitive** embeddings of $L$ in $\mathbb Z^4$ in expected polynomial bit time. Primitive means that the six two-by-two minors of the output vectors have gcd one. This applies even if the input form has no primitive embedding in three dimensions. The construction below also prescribes the complementary plane by a public Gaussian residue. A quaternion construction samples directly and handles primitive Grams with nonunit leading coefficient. [Gaussian sampler](ternary_form_reduction.py) · [Quaternion construction](quaternion-gram.md)

Existence, counts and the genus of the complementary binary form are classical. In the present range Hunsucker's Theorem 3 gives $96\,|\operatorname{Pic}(\mathbb Z+N\mathbb Z[i])|$ ordered primitive embeddings. We use this count to verify completeness of an explicit sampling procedure. [Hunsucker, Theorems 1–3 (1971)](https://doi.org/10.4064/aa-19-4-321-325) · [Pall's matrix completion (1949)](https://doi.org/10.4153/CJM-1949-032-9)

## Constructing the embedding

Put $D=N^2$, $O=\mathbb Z[i]$. Sample $z\bmod N$ uniformly with $|z|^2\equiv-A\pmod N$. One verified solution comes from the existing modular Gaussian norm routine. Multiplication by $w/\bar w$ for a uniform unit $w\bmod NO$ gives uniformity. Multiply an integer lift of $z$ by $1+kN$ to get $|z|^2\equiv-A\pmod D$, using the unit $2|z|^2$ to determine $k\bmod N$. The choice of lift leaves the final orbit unchanged. [Modular norm routine](four_square_sampler.py) · [Local norm fibers](ternary-shapes.md)

In the rational span of $L$ with basis $e_1,e_2$ take $d=(-B e_1+A e_2)/D$. This generates $L^*/L$, has norm $A/D$, and pairs with $e_1,e_2$ as $0,1$. In a perpendicular Euclidean plane put $J_z=\mathbb Z z+NO$. This lattice has index $N$ in $O$. Enlarge $L\oplus J_z$ by the vector
$$
(d, iz/N).
$$
The second coordinate pairs integrally with $z$ and $NO$, and the whole vector has integral norm $(A+|z|^2)/D$. Each coordinate has exact order $D$ modulo its own lattice. For the second coordinate, write $z=h e$ for an integer $h$ coprime to $N$ and a primitive Gaussian vector $e$. Take $t$ with $\det(e,t)=1$. In the basis $e,Nt$ of $J_z$, its coordinates are $h(-B',A')/D$ for $A'=|e|^2$, $B'=N(e\cdot t)$. Since $hA'$ is a unit modulo $D$, the claim follows. The enlarged lattice is integral unimodular, and both summands are primitive in it. Gram LLL yields an orthonormal basis in fixed dimension, so the first two vectors give an exact integral embedding. The code randomizes by a uniform signed permutation at the end. [Orthogonalization](gram_embedding.py)

## Why the sampler is uniform

The number of norm-$-A$ residues modulo $N$ is
$$
\lambda(N)=N\prod_{r\mid N}\left(1-\frac{\chi_{-4}(r)}r\right).
$$
It follows by CRT from the norm fibers in the split and inert quadratic algebras, including the $r$ lifts per added prime power. If $z'=z+Nw$ is another lift of the required norm modulo $D$, then $z\cdot w=0\pmod N$. The two vectors $z,iz$ form a basis modulo $N$. Thus $iw\in J_z$, so the glued lattice is identical. Multiplication of the seed by a global Gaussian unit gives an isometric marked construction.

Conversely, an oriented isometry fixing the ordered input vectors takes $J_z$ to $J_{z'}$ by a plane rotation. Both lattices generate $O$ as Gaussian ideals, since their generators include $NO$ and an element of unit norm modulo $N$. The rotation must be multiplication by a global Gaussian unit $\varepsilon$. The first component of the glue has order $D$, so its coefficient in the second gluing must be $1\pmod D$. Comparing the plane components gives $i(\varepsilon z-z')/N\in J_{z'}$, and in particular $z'\equiv\varepsilon z\pmod N$. For a fixed seed the only such stabilizing unit is $1$, since $N>1$ is odd and $z$ is a unit.

There are consequently $\lambda(N)/4$ distinct proper rotation orbits, each of size $192=2^3\,4!$. Their total $48\lambda(N)$ agrees with Hunsucker's count because $|\operatorname{Pic}(\mathbb Z+NO)|=\lambda(N)/2$. All seeds and all orbits have equal mass. Using all signed permutations at the final step also preserves the resulting uniform distribution.

## Dropping a coordinate

A coordinate can be deleted from both vectors exactly when the complementary lattice contains an integral unit vector. In the standard ambient lattice these are just the signed coordinate axes. In $J_z$ they are the Gaussian units. This is equivalent to the line of $z$ being the real or imaginary axis modulo $N$. Writing $\rho_N(-A)$ for the number of scalar square roots of $-A\bmod N$, the exact probability is
$$
\frac{2\rho_N(-A)}{\lambda(N)} .
$$
For distinct primes $p,q$ of equal bit length and a three-dimensional solution, this is $8/[(p-\chi_{-4}(p))(q-\chi_{-4}(q))]$, of order $1/N$. The native four-dimensional solver includes every primitive embedding and nevertheless seldom attains one of the three-dimensional solutions by this operation. The same count gives zero when the ternary instance is impossible. For example $\operatorname{diag}(1,9)$ embeds primitively in four dimensions with vectors $(1,0,0,0)$ and $(0,-2,2,-1)$, while $-1$ has no square root modulo $3$.

Allowing longer rational directions in the complementary plane has a precise arithmetic requirement. If $A'X^2+2B'XY+C'Y^2=T^2$ with $\gcd(T,N)=1$, then $(A'X+B'Y)^2\equiv A'T^2\pmod N$. As $h^2A'\equiv-A\pmod N$, this supplies a root of $-A$. An algorithm finding such a direction on general soluble instances would supply exactly the information needed for the ternary construction. [Ternary equivalence](prescribed-ternary-embedding.md)

The finite test compared all 1,344 ordered embeddings for four selected Grams with exhaustive four-dimensional shells and recovered each with the same multiplicity. It checked 200 individual embeddings for the coordinate-deletion criterion and ran the native sampler at odd radii of up to 216 bits. These are tests of the distribution and identities. [Reproducible checks](check_ternary_form_reduction.py)

## Sampling every embedding of the Jacobi-selected form

The three-dimensional factoring test uses ordinary embeddings, including imprimitive ones, of $A=x^2+1$, $B=-Nx$, $C=N^2$, where $\gcd(A,N)=1$. We can sample this **same set of constraints in four dimensions**, including every content, in expected polynomial bit time. Let $g$ be the gcd of the six minors. Cauchy–Binet gives $g\mid N$. The identity $Av_j-Bu_j=\sum_i u_i(u_i v_j-u_j v_i)$ shows that $g\mid v_j$ for all $j$. Thus $g=\gcd(v_1,v_2,v_3,v_4)$. Writing $m=N/g$, division of $v$ by $g$ identifies its fiber with the primitive embeddings of $(x^2+1,-mx,m^2)$. Their number is $48\lambda(m)$, including $48$ at $m=1$. [Four-dimensional count and bijection](quaternion-gram.md)

The existing one-scalar four-square content sampler proposes $g$ with weight $m\prod_{r\mid m}(1+1/r)$. Accept if a primitive four-square sample of squared norm $m$ maps to a primitive point under the Hopf map. Its acceptance probability is $\prod_{r\mid m}(1-\chi_{-4}(r)/r)/(1+1/r)$, at least $(\varphi(N)/N)^2$. This gives exactly the desired $\lambda(m)$ law without factoring. After acceptance sample the primitive frame independently and multiply its second vector by $g$. [Implemented sampler](ternary_form_reduction.py) · [Proposal](four-square-incidence.md) · [Hopf acceptance proof](ternary-sphere.md)

For distinct odd primes $N=pq$ the probability of a proper gcd is exactly $(\lambda(p)+\lambda(q))/[(1+\lambda(p))(1+\lambda(q))]$. For comparable primes this is of order $1/\sqrt N$. In three dimensions, choosing the Jacobi symbol negative makes the corresponding exact probability $2/3$ for the same $N,x$. At $N=77$, $x=1$, the four-dimensional probability is $20/117$. The unresolved step is obtaining a factor-bearing three-dimensional embedding from the public input. [Three-dimensional count](mixed-ternary-embeddings.md)
