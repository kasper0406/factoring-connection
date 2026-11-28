# Counting a specified three-square residue with a factored conductor

Let $A=x^2+1$, let $T=pq$ have two known distinct odd prime factors prime to $A$, and let $\gcd(m,T)=1$. We can count the points of squared norm $AT^2$ in **one specified coordinate class modulo $m$**, including when $m$ is even or shares primes with $A$. The algorithm can count all points or just primitive points. After enumerating the short base shell of norm $A$, it uses $B(p+1)$ forward images and $q+1$ backward images for one admissible query, where $B$ is the number of base points used. The forward table serves every query at the same modulus. For $p,q$ of size about the original factoring input $N$ this still costs a power of $N$. Unrestricted primitive sampling is already possible with an unfactored $T$. [Individual counter](factored_ternary_hecke.py) · [Whole-shell sampler](fixed_field_walk.py) · [Fixed-field reduction](fixed-field-packing.md)

## The one-prime correspondence

Identify triples with pure Hamilton quaternions. For an odd prime $s$, choose one integral quaternion $h$ of norm $s$ for each right-associate class under the eight integral units. There are $s+1$ classes, one for each possible image line in the matrix algebra modulo $s$. The implementation enumerates them deterministically by tabulating pairs of integers by squared norm. Among four coordinates of total norm $s$, one has parity different from the other three. Put that coordinate first and fix the overall sign. This chooses exactly one quaternion per class and uses $O(s)$ integer operations up to logarithmic factors. Quaternion rotations and the Hecke recursion for unrestricted three-square counts are classical. [Quaternion image-line correspondence](four-square-incidence.md) · [Ellenberg, Michel and Venkatesh](https://www.epfl.ch/labs/tan/wp-content/uploads/2018/10/pointsonspheres.pdf) · [Hirschhorn and Sellers](https://doi.org/10.1016/S0012-365X(98)00288-X)

Fix a source norm $n$ prime to $s$. Write $C_n(c)$ for its count in the class $c\bmod m$, either ordinary or primitive throughout. Put $\chi_s=(-n/s)$ and set $k_s=\chi_s$ for ordinary counts and $k_s=1+\chi_s$ for primitive counts. Then
$$
 C_{ns^2}(c)=\sum_{[h]} C_n(s^{-2}\bar h c h)-k_s C_n(s^{-1}c).
$$
All inverses in the queries are taken modulo $m$. To prove the formula, push each source vector $v$ forward by $v\mapsto h v\bar h$. For a target $w$ primitive at $s$, its trace-zero matrix modulo $s$ is a nonzero nilpotent, with a unique possible image line for $h$. To see that this line suffices, write $h=U\operatorname{diag}(s,1)V$ over $\mathbb Z_s$. In the resulting basis the matrix of $w$ has entries $(a,b;c,-a)$ with $s\mid a,b$ and $c$ a unit. Its determinant $ns^2$ implies $s^2\mid b$. Thus the entries of the inverse image are integral, being conjugate to $(a/s,b/s^2;c,-a/s)$. Only $s$ can occur in a denominator.

For an imprimitive target $w=s t$, the inverse image is $h^{-1}th$, which is integral exactly when the image line of $h$ is invariant under $t\bmod s$. The characteristic polynomial is $Y^2+n$ and has $1+\chi_s$ invariant lines. The raw sum therefore counts each such target $1+\chi_s$ times. Away from $s$ the map preserves coordinate primitivity. This proves both corrections, including the inert case where ordinary counting adds the missing $s t$.

For $p$ the algorithm pushes the base shell forward and accumulates residues with the correction above. For $q$ it reads the corrected table at the inverse images of the requested residue. Either prime can go first. The stored rotations and table occupy $O(p+q+B(p+1))$ integer words during construction. Summing over all residues gives the ordinary count $r_3(A)\prod_{s\mid T}(s+1-\chi_s)$ and the primitive count $r_3^*(A)\prod_{s\mid T}(s-\chi_s)$, now with $\chi_s=(-A/s)$. The symbols can differ at the two primes, even if both split in the Gaussian integers.

## The common genus average at primitive residues

For a residue admitting a local lift of norm $M=AT^2$ at every prime of $m$, assume its coordinates generate the unit ideal modulo $m$. Such residues form one proper genus orbit. At an odd prime $\ell^a\Vert m$, primitive trace-zero matrices of fixed determinant are cyclic modulo $\ell$ and hence conjugate to the same companion matrix over $\mathbb Z/\ell^a$. Integral lifts of the conjugating matrices act by proper isometries. The number of residues in this orbit is
$$
 \ell^{2a-2}\begin{cases}
 \ell^2+\ell(-A/\ell),&\ell\nmid A,\\
 \ell^2-1,&\ell\mid A .
 \end{cases}
$$
Indeed each primitive solution modulo $\ell$ has $\ell^2$ lifts per higher level. At level one, a nonzero norm has the displayed ternary count, and the punctured zero cone has $\ell^2-1$ points.

For $2^a\Vert m$ the orbit size is $3\cdot4^{a-1}$. Here $M$ is $1,2$ or $5\bmod8$. Proper signed permutations connect the eligible patterns modulo four, except that at $M=1\bmod8$ the two even coordinates can both be zero or both be two modulo four. Conjugation by the quaternion $(1,1,1,0)$ of norm three takes an axial residue in the first pattern to the second. For two exact $2$-adic vectors $v,w$ agreeing modulo four, let $z=(v+w)/2$. Its norm is a unit or twice a unit. The proper integral isometry $-s_z$, with $s_z$ reflection in $z$, takes $v$ to $w$. There are three eligible parity patterns and twelve residues modulo four. For $a\ge2$, a class modulo $2^a$ lifts to exact norm $M$ precisely when its norm agrees modulo $2^{a+1}$. Of its eight children four satisfy the next norm bit, since an odd coordinate is present. This proves the count by lifting.

Let $S$ be the product of the local orbit sizes. Genus weights are inverse proper automorphism counts. Because the ambient lattice has one proper class, orbit–stabilizer makes the genus average the ordinary uniform average over these $S$ residues. Its primitive numerator is the full primitive shell count. For the ordinary numerator replace $r_3(A)$ in the ordinary formula by the number of base points with coordinate gcd coprime to $m$. Indeed, both rotation and scaling by $s\mid T$ preserve local primitivity at every prime of $m$. The implementation evaluates these two averages given the factorization of $m$.

For a fixed packed first vector and a test prime $\ell$, the candidate classes also share their proper spinor genus. They agree at every conductor prime except $\ell$. There the norm is a unit, and the integral stabilizer's spinor norms contain every local unit, so the change of class at $\ell$ has trivial spinor obstruction. This argument concerns classes at a fixed first vector. [Coset genus weights and the unit stabilizer calculation](ternary-coset-spinor-genus.md)

## Exact controls

All 47,526 small primitive and ordinary residue counts agreed with direct three-square enumeration. These tests also verify the orbit sizes and genus averages at even, odd, ramified and mixed residue moduli. At $N=77,x=1,T=257\cdot397$, a single table used 3,096 forward images and had 3,072 occupied keys for primitive counts. The 14 candidate residues modulo $1240\cdot13$ used 398 backward rotations each and contained six points in all. The full list of forward quaternion pairs would use 1,232,208 images. Fine classes may also admit a cheaper direct coordinate scan. At $N=15,T=29\cdot101$, the fixed-field reduction recovered 3 using 140 individual primitive decisions at six moduli. [Reproduce all checks](check_factored_ternary_hecke.py)
