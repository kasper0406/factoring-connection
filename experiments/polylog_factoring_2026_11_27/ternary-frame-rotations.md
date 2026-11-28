# Rational rotations of the 72 ternary frames

Let $N=pq$ have distinct odd prime factors. Choose $x$ with $\left(\frac{-A}{N}\right)=-1$, where $A=x^2+1$. Write $r$ for the prime where $-A$ is a square. The Gram $(A,-Nx,N^2)$ has 24 coordinate frames and 48 other frames. The coordinate frame $u_0=(x,1,0)$, $v_0=(-N,0,0)$ is public. Exactly **two reduced denominators** occur among the rational rotations to the 72 frames, namely 1 and $r$. A rotation keeping $u_0$ has a more restrictive binary norm test, which can exclude all 48 at once. [Frame count and public input construction](mixed-ternary-embeddings.md) · [Exact rotation checks](check_ternary_frame_rotations.py)

## Denominator of the transition

For a target $u,v$ the unique proper rotation taking $(u_0,v_0)$ to $(u,v)$ has columns

$$
-v/N,\qquad u+xv/N,\qquad (u\mathbin\times v)/N.
$$

The first two columns follow from the two prescribed images and are orthonormal by the Gram identities. Their cross product gives the third. The first column already has reduced common denominator $N/\gcd(v)$, and this clears every entry in the remaining columns. The frame count gives 24 rotations with denominator 1 and 48 with denominator $r$. Starting from any other coordinate frame just composes with an integral signed permutation. A walk that stays on integral frames therefore has to use a step of reduced denominator $r$ on its first exit from the coordinate frames. Even allowing rational intermediate frames, a product of rotations whose denominators are coprime to $N$ cannot make that exit.

## Keeping the first vector

Rational rotations have the form $z\mapsto h z\bar h/\operatorname{Norm}(h)$ on pure imaginary quaternions, with $h$ a primitive integral quaternion unique up to sign. The reduced denominator has the same odd part as $\operatorname{Norm}(h)$. Indeed, modulo an odd prime $\ell$ dividing the norm, $h$ is a nonzero rank-one matrix in the split quaternion algebra. Its adjugate is nonzero with image $\ker h$. Some matrix moves that image outside $\ker h$. The identity contributes zero to $h I\bar h$, so some trace-zero matrix has nonzero image under $T\mapsto hT\bar h$. Thus $\ell$ cannot cancel from every entry of the rotation numerator. This denominator is also the classical coincidence-site lattice index for cubic lattices. [Grimmer, Bollmann and Warrington (1974)](https://journals.iucr.org/paper?a10612) · [Grimmer on coincidence rotations (1974)](https://journals.iucr.org/paper?a10740)

Identify $u_0$ with $U=xi+j$. A rotation fixes $U$ precisely when $h$ commutes with $U$. Thus $h=(a,bx,b,0)$ with coprime integers $a,b$. For $b\ne0$ the norm is $a^2+A b^2$. This primitive sum of three squares is never divisible by four. When the norm is even, every entry in the usual quaternion rotation numerator is even. It follows that a noncoordinate frame with $u=u_0$ exists **if and only if**

$$
a^2+A b^2\in\{r,2r\},\qquad b\ne0,\qquad \gcd(a,b)=1.
$$

For the reverse implication the rotation fixes $U$, and its denominator $r\mid N$ ensures its image of $v_0$ is integral. Denominator $r>1$ then certifies that the result is noncoordinate. In particular $A\le 2r$ is necessary.

The test applies to the entire signed-permutation orbit of $u_0$, since its zero coordinate permits an orientation adjustment. In fact the orbit contains the first columns of **all 48 noncoordinate frames or of none**. For $x\ge2$ it has 24 vectors. One noncoordinate rotation fixing $u_0$ supplies a second via its inverse. They are distinct because a half-turn would have $a=0$, whereas the norm equation and the unit $A\bmod r$ rule this out. The signed permutations thus give at least $2\cdot24$ solutions. For $x=0,1$ every vector of norm $A=1,2$ belongs to the public orbit already. The checks also count the frames fixing $u_0$ exactly via the binary norm equations, identifying $h$ and $-h$.

For uniform $x\in[0,N)$ conditioned on the accepted Jacobi sign, the probability of passing $A\le2r$ is at most $8(1+\sqrt{2\max(p,q)})/N$, since the acceptance probability exceeds $1/8$. On balanced inputs this is $O(N^{-3/4})$. This estimate concerns the uniform challenge distribution. A deliberately short $x$ can behave differently. For example at $N=77,x=1$, the good prime is $r=11$, and $11=3^2+2\cdot1^2$. All 48 frames use the public first-column orbit. At $N=65,x=4$ the good prime is $13$ and $26=3^2+17\cdot1^2$, whereas at $N=65,x=15$ none of the 48 uses the public orbit. Finding a qualifying binary norm without knowing $r$ remains a search problem. [Public search for small Jacobi-negative parameters](short_jacobi_grams.py)

Sampling a fresh primitive first column uniformly from the sphere of squared radius $A$ has a separate limitation. At most 48 such columns can extend to one of the noncoordinate frames, even with a perfect completion test. A first column of a valid frame is primitive because $\gcd(A,Nx)=1$. For $A>1$ in this family Gauss's formula gives $r_3^{\rm prim}(A)=12h(-4A)$. The ring-class formula and Siegel's lower bound give $h(-4A)\gg_\epsilon A^{1/2-\epsilon}$, including nonfundamental discriminants. Indeed for $-4A=f^2D_0$ the conductor contributes $\gg_\epsilon f^{1-\epsilon}$, and the fundamental class number contributes $\gg_\epsilon |D_0|^{1/2-\epsilon}$. Summing the resulting upper bound over $x<N$, and allowing the single $x=0$ separately, shows that this uniform-first-column method succeeds with probability $O_\epsilon(N^{-1+\epsilon})$ averaged over the uniform accepted challenges. The constant from Siegel's theorem is ineffective. [Gauss's primitive formula in Voight, Theorem 30.1.3](https://jvoight.github.io/quat-book.pdf) · [Ring-class formula in Cox, Theorem 7.24](https://doi.org/10.1002/9781118400722)

## What a public norm-$N$ quaternion detects

One can draw a uniformly random primitive integral quaternion $h$ of norm $N$ without knowing either prime. The numerator $T=h U\bar h$ gives the native trial $\gcd(N,T_i,T_j,T_k)$. At a prime $s\mid N$, $T$ vanishes exactly when the kernel line of the rank-one matrix $h\bmod s$ is invariant under $U$. The kernel line is uniform on $\mathbb P^1(\mathbb F_s)$ by the quaternion ideal sampler, or by applying conjugation to its uniform image line. The characteristic polynomial of $U$ is $X^2+A$. Thus the good prime has two invariant lines among $r+1$, and the other prime has none. The gcd gives $r$ with probability **exactly $2/(r+1)$** and otherwise equals 1. The norm-$N$ rotation itself never takes the coordinate frame to an integral frame. [Uniform quaternion sampler and proof](mellin-four-squares.md) · [Native gcd trial](check_ternary_frame_rotations.py)

The same local condition has a simpler scalar implementation. In companion coordinates the eigenlines are the two roots of $k^2+A$ at $r$, so a uniform affine $k\bmod N$ gives a proper gcd with probability $2/r$. The quaternion transition supplies no probability improvement here. This comparison leaves open input-dependent selection of a line or of a different first vector.

Across all 119 Jacobi-negative challenges on six small semiprimes, direct enumeration checked 8,568 rotations and both binary norm equations. Ten targets kept the public first-column orbit for all 48 noncoordinate solutions, and 109 kept it for none. All 82,624 norm-$N$ quaternion numerator tests matched the exact predicted frequencies. A separate run of 150 native random trials at $N=77,x=1$ returned the prime 11 on 27 draws. [Reproduce the checks](check_ternary_frame_rotations.py)
