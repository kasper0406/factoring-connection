# Which plane would contain a factor-revealing point?

Let $N$ be an odd distinct semiprime and $A=x^2+1$ with $(-A/N)=-1$. Fix an integral $u$ of squared length $A$. Packed points have the form $W=u-Qv$ with $u\cdot v=-Nx$ and $v^2=N^2$. There are 48 noncoordinate frames, each with $\gcd(v)=d$, where $d$ is the prime at which $-A$ is a nonsquare. Write $r=N/d$. The packing yields $W^2=M=(x+QN)^2+1$. The results below address the idea of generating a point by representing $M$ with a binary quadratic form in some rational plane. [Packing identity](packed-ternary-decision.md) · [Frame classification](mixed-ternary-embeddings.md)

## Two gcds for a homogeneous plane

If any noncoordinate packed witness belongs to $w\cdot W=0$, for primitive integral $w$, then **either $\gcd(w\cdot u,N)$ or $\gcd(w^2,N)$ is a proper factor**. Both are computable before searching for a point in the plane. This works at every height of $w$. First test $\gcd(Q,N)$ to handle $7\mid N$.

Indeed $d\mid v$ gives $d\mid w\cdot u$. If $r$ also divides the pairing, $w$ is orthogonal to $u$ modulo both primes. At $r$ the plane equation now implies $w\perp v$. The two vectors $u,v$ are independent there, because $u^2=A$ is a unit and $v$ is a nonzero null vector with $u\cdot v=0$. Their common orthogonal line is spanned by $v$, so $w^2=0\bmod r$. At $d$ the binary quadratic form on $u^\perp$ is anisotropic: its determinant is $A$ up to squares, and $-A$ is a nonsquare. The primitive $w$ is nonzero there, giving $d\nmid w^2$. In particular a coordinate plane cannot help, since $W_j=0$ implies $u_j=v_j=0$ from $|u_j|<Q$.

Allowing explicit polynomial dependence on the input does not make a short prelisted menu of such planes cover ordinary balanced products. Suppose its integer normal is $w(N)$, where the three entries are polynomials. Divide out the largest shared formal power of $N$ and put $c=w(0)\ne0$. At a prime $s\mid N$, primitive normalization of the evaluated normal preserves the line of $c\bmod s$ unless $s$ divides every entry of $c$. The latter primes already divide $c^2$. If $c\cdot u\ne0$, a successful plane forces $d\mid c\cdot u$ or $d\mid c^2$. If $c\cdot u=0$, either a side divides $c^2$ through common content or the two-gcd proof gives $r\mid c^2$. Thus, per menu item $(u,w)$, only $O(H/\log R)$ possible sides in $[R,2R)$ matter, where $H$ bounds the bit lengths of the two nonzero integers being counted. Across $J$ items, the prime number theorem gives $O(JH/R)$ coverage, including adaptive selection among them. Normals computed by unrestricted input-dependent procedures are outside this bound.

## An exact solution for shifted planes

For a rational affine plane clear denominators to write $w\cdot W=h$. If $e=(w\cdot u-h)/Q$ is nonintegral the section is empty. Otherwise put $B=-Nx$, $C=N^2$, $L=u\times w$, $D=L^2$ and $c=u\cdot w$. A partner has $w\cdot v=e$. If $D>0$, its Gram determinant gives
$$
(L\cdot v)^2=DC-B^2 w^2 +2Bec-Ae^2.
$$
Check this integer for being a nonnegative square $t^2$. The only candidates are
$$
Dv=(Bw^2-ec)u+(Ae-Bc)w\ \pm t L .
$$
Checking integrality recovers every partner in polynomial bit time, using no factors. If $D=0$, the extra plane is inconsistent or is exactly the public conic plane $u\cdot W=A+QNx$. Thus a proposed affine plane gives a direct two-candidate test in every transverse case. A rational curve on the sphere of explicit degree $k$ admits the analogous intersection: substitute its parameter in the public conic plane equation, find the rational roots and include parameter infinity. This takes polynomial time in $k$ and the coefficient bit lengths for a supplied parametrization; if the equation is identically zero its image lies on the original conic.

## What a publicly factored square radius contributes

At $x=0,1$ the fixed-vector conic is exactly $(N-At)^2+Az^2=N^2$. At $A=2$ the first coordinate is automatically odd. Write $\mathcal O=\mathbb Z[\sqrt{-A}]$, a norm-Euclidean ring in these two cases, and $S_L=\{z\in\mathcal O:\operatorname{Norm}(z)=L^2\}/\mathcal O^\times$. For every $T$ coprime to $N$, multiplication gives an exact bijection $S_N\times S_T\longrightarrow S_{NT}$. Given any padded point $z$, compute $z_T=\gcd_{\mathcal O}(z,T^2)$ and $z_N=z/z_T$. Unique factorization shows this picks out the full prime-ideal contribution over $T$, including ramified and inert primes. It works without the factorization of either integer. Conversely, given the rational factorization of $T$, uniformly sampling $S_T$ takes expected polynomial time by splitting its known split primes via Cornacchia and drawing their conjugate valuations uniformly. Ramified and inert prime contributions have a single orientation. [Cornacchia's algorithm](https://mathworld.wolfram.com/CornacchiasAlgorithm.html)

On a mixed input the three orientations in $S_N$ are $N$, $d\pi^2$ and $d\bar\pi^2$ up to units, where $r=\pi\bar\pi$. Padding leaves the factor-revealing fraction at $2/3$. More strongly, the content gcd with $N$ stays exactly the same on each padded point since the pad has invertible norm modulo $N$. At any auxiliary modulus coprime to $NT$, the distributions of the normalized residues $z/L$ in the norm-one group modulo units obey $\mu_{NT}=\mu_N*\mu_T$. Each Fourier coefficient is multiplied by the corresponding known padding coefficient. Existence information can disappear completely. For $N=15$, $A=1$, modulo 7 the norm-one group modulo Gaussian units has order two. With $T=13$ the element $(3+2i)^2/13$ reduces to $2+2i$, whose square is $i$. Both classes now have lifts coming from the public coordinate branch alone.

The checks compare 1,152 exact affine sections with exhaustive frame lists. On 1,536 homogeneous plane incidences the first gcd split 1,116 and the squared-normal gcd split the other 420. They also exhaust 592 ordered binary norm product pairs, verifying coprime padding and factor preservation in both rings. [Plane routines](packed_plane_sections.py) · [Run the checks](check_packed_plane_sections.py)
