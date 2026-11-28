# Norm relations between auxiliary primes near $\sqrt N$

Let $N=pq$ for distinct odd primes $p<q$, with $N>16$. Choose any set of odd prime powers $A$ coprime to $N$ in $[\sqrt N,2\sqrt N]$. Join $A$ and $B$ when

$$
AB=x^2+Ny^2,\qquad x\ge0,\quad y>0.
$$

**Each vertex has at most one neighbor, allowing a loop to count as its own neighbor.** Given a square root $t$ of $-N$ modulo $A$, choose the even representative $x$ of the pair $t,-t$ in $[0,A)$. The only candidate partner is $B=(x^2+N)/A$. Check whether it belongs to the chosen vertex set. An unsuccessful membership test exhausts that vertex.

Indeed, $AB<4N$ since $N$ is nonsquare, so $y=1$. Also
$$
x^2=AB-N\le 2A\sqrt N-N=A^2-(A-\sqrt N)^2<A^2 .
$$
The same argument gives $x<B$. The oddness of $A,B,N$ forces $x$ even. At a coprime odd prime power the two roots, whenever they exist, are complementary and have opposite parity in this interval. These facts give uniqueness at both endpoints, including when the vertices are chosen after observing earlier edges. The degenerate representations with $y=0$ just give the ordinary square $AB=x^2$.

### Exactly when a product of edges splits

Pick any multiset of edges for which the product of their endpoint labels is a square $s^2$. Then $z=\prod_e x_e/s\bmod N$ satisfies $z^2=1$, so $\gcd(z-1,N)$ tests for a split. A repeated edge contributes sign $+1$ in pairs. Each underlying prime supports at most one vertex, since consecutive positive powers differ by a factor at least three. On a nonloop edge $A=\ell^i,B=r^j$, either odd exponent forces the edge multiplicity to be even. The edges that can occur with odd multiplicity are exactly the loops and the nonloops with both exponents even. In each case $m=\sqrt{AB}$ is an integer and

$$
(m-x)(m+x)=N,\qquad \sqrt N\le m\le2\sqrt N.
$$

The factor pair $1,N$ would give $m=(N+1)/2>2\sqrt N$. Thus $m=(p+q)/2$, and this single edge gives both primes as $m-x,m+x$. Every such edge contributes the same mixed sign up to a global sign. With positive coordinates it contributes exactly $(q-p)/(p+q)\bmod N$. Consequently a multiset splits exactly when it uses an odd total multiplicity of these edges. For prime vertices the only possible splitting edge is the loop at the midpoint itself. A midpoint in the vertex set and band always has this loop, via $x=(q-p)/2$. Balanced primes $q\le 2p$ always put their midpoint in the band.

Even prime-power exponents matter. For $N=2{,}823{,}013$, the distinct vertices $41^2,43^2$ share the edge $x=534$. From either endpoint and its root one recovers $41\cdot43-534=1229$ and $41\cdot43+534=2297$. Supplying an arbitrary nearby prime power does not ensure that its product with the computed partner is a square.

### A single seed captures cycle signs outside the band

Take any odd prime $\ell\nmid N$ with a supplied root $t^2=-N\bmod\ell$. The remaining vertices may now be arbitrary odd primes coprime to $N$. Edges may use any $y>0$ and multiple representations are allowed. For every connected component reached from this seed, one two-dimensional lattice suffices to either split $N$ or certify that every cycle in the component has global sign $\pm1$.

Lift $t$ to a root $t'$ modulo $\ell^2$ by Hensel's formula. Gauss-reduce the lattice $x=t'y\bmod\ell^2$ with the positive norm $x^2+Ny^2$. Every norm is a multiple of $\ell^2$. If the minimum equals $\ell^2$, take a shortest vector and compute $z=x/\ell\bmod N$. A mixed $z$ splits $N$. A larger minimum or a global $z$ certifies every cycle through the component has global sign. Given $t$, this requires a polynomial number of bit operations in $\log(N\ell)$. A mixed seed supplies its own self-loop even if the proposed graph omits that edge.

To prove the certificate, work in $\mathcal O=\mathbb Z[\sqrt{-N}]$. The seed ideal $\mathfrak l=(\ell,\sqrt{-N}+t)$ has square represented by the lifted lattice above. Odd ideals coprime to $N$ in this construction are prime to the conductor and have unique prime-ideal factorization. A vector of norm $\ell^2$ generates $\mathfrak l^2$ by equality of indices. The only units are $\pm1$.

For a closed walk based at $\ell$, each oriented edge gives an element $\alpha=x+y\sqrt{-N}$ with norm the product of its two endpoints. Even for a loop its ideal is a product of two primes above the endpoints. A loop with $y>0$ cannot use both conjugate ideals, since then $\alpha/\ell$ would be a unit. Choose $\alpha$ or its conjugate successively so the two ideal factors at each intermediate vertex are conjugates. Dividing their product by the product of the visited primes, one per step, leaves the unit ideal, $\mathfrak l^2/(\ell)$, or its inverse. The normalized element therefore equals a global unit, or differs from the normalized seed generator or its inverse by a global unit. Reduction modulo $\sqrt{-N}$ sends conjugate edges to the same $x\bmod N$. Every cycle sign thus lies in $\{\pm1,\pm z\}$, with the second pair available only when the seed square is principal. A path to a different cycle base followed by the reverse path contributes sign $+1$, which extends the statement throughout the connected component. Square products for prime vertices decompose into closed walks by parity of the degrees.

The repository previously identified ambiguous reduced forms that produce a factor and characterized relations among split prime ideals by smooth norms. The band matching and the one-seed certificate specify which relations can contribute a new sign in this two-endpoint model. The argument does not cover relations with three or more independent prime labels. [Ambiguous forms and smooth norm relations](../../notes/genusleak_results.md)

The exhaustive band check uses 421 balanced semiprimes with $7\le p<q<200$. Among 9,290 prime-power vertex occurrences it found 558 edges and 57 loops. Restricting to primes gave 435 edges and 44 loops. Outside the band, direct norm enumeration agreed with the seed lattice at all 132 split primes tested over eight moduli. There were 103 nonprincipal seed squares, 11 global seed signs and 18 mixed seed signs. All 495 graph closure checks agreed with the seed classification, including 38 mixed closures. [Runnable checks](check_near_square_prime_norms.py)
