# Local ranks of the matrix power map

Given $N=pq$ with distinct primes and a public $d$-by-$d$ matrix $A$ modulo $N$, the full Jacobian of $A\mapsto A^N$ is computable by binary powering. Its column for the matrix unit $E$ is the upper-right block of $\left(\begin{smallmatrix}A&E\\0&A\end{smallmatrix}\right)^N$, read in row-major order. This takes polynomial time in $d$ and $\log N$ and uses no characteristic-specific exponent. At $N=91$, the fixed matrix $A=\left(\begin{smallmatrix}0&-1\\1&3\end{smallmatrix}\right)$ has $A^{91}=\left(\begin{smallmatrix}25&13\\78&77\end{smallmatrix}\right)\bmod91$. The Jacobian has ranks 2 modulo 7 and 0 modulo 13. The gcd of its entries with 91 returns **13**. The same factor is already visible in the off-diagonal entries and diagonal difference of $A^{91}$. [Public computation](check_matrix_power_jets.py) · [Classical block derivative formula, Mathias](https://doi.org/10.1137/S0895479895283409)

## Which perturbations survive

Fix $r\in\{p,q\}$ and $K=\mathbb F_r$. Let $R$ be a finite-dimensional commutative local $K$-algebra with residue field $K$ and maximal ideal $\mathfrak m$ satisfying $\mathfrak m^s=0$. Let $A'$ be **any** matrix over $R$ reducing to $A$ over $K$, including perturbations that do not commute with $A$. If $b$ is the largest Jordan block size of $A$ over an algebraic closure and $bs\le r$, then

$$\{(A')^N:A'\equiv A\bmod\mathfrak m\}=\{G A^N G^{-1}:G\equiv I\bmod\mathfrak m\}.$$

In particular $ds\le r$ suffices uniformly over $A$. To prove the claim, put $f_0=\operatorname{minpoly}(A)$. Every entry of $f_0(A')$ lies in $\mathfrak m$, so $f_0(A')^s=0$. After scalar extension, primary decomposition for $f_0(X)^s$ puts $A'$ into blocks $\lambda I+U$ with $U^{bs}=0$. Then $(\lambda I+U)^r=\lambda^r I$. Thus $C=(A')^N$ is annihilated by the separable minimal polynomial $f$ of $B=A^N$, whose roots are the distinct $\lambda^N$.

The conjugator exists over $R$ itself. Write $(f(X)-f(Y))/(X-Y)=\sum h_{ij}X^iY^j$ and $V=\sum h_{ij}C^iB^j$. Telescoping gives $CV=VB$, and reduction modulo $\mathfrak m$ gives $V\equiv f'(B)$. Since $f$ is separable, $f'(B)$ is invertible. Thus $G=V f'(B)^{-1}$ is an invertible intertwiner reducing to the identity. The converse uses $A'=GAG^{-1}$.

For $s=2$, use $R=K[\epsilon]/\epsilon^2$. Conjugating $B$ by $I+\epsilon X$ gives $B+\epsilon(XB-BX)$. Consequently the Jacobian image is exactly the space of commutators with $B$ if $2b\le r$. Its rank is

$$d^2-\sum_{\xi}m_\xi^2,\qquad m_\xi=\sum_{\lambda^N=\xi}\operatorname{mult}_A(\lambda),$$

because a matrix commuting with semisimple $B$ preserves each eigenspace. Entries in a fixed coordinate system can still change throughout this orbit.

For $2d\le\min(p,q)$, a difference in the two ranks requires different multiplicity profiles of the characteristic polynomial or an eigenvalue collision on one side. The multiplicities can be compared by squarefree decomposition using polynomial gcds modulo $N$. Either the same degree profile results locally or a nonunit leading coefficient supplies a factor. In characteristic $p$, two distinct eigenvalues with the same $N$th power are nonzero with ratio of order $q$. At least one eigenvalue then has multiplicative order divisible by $q$ in an extension of degree $\le d$. Hence $\operatorname{ord}_q(p)\le d$. The roles interchange at $q$. This is the existing Vandermonde obstruction applied to the entire noncommutative differential. [Earlier power-substitution rank test](cartier-paths.md#minimum-degree-for-a-local-rank-drop-in-power-substitution)

For uniformly chosen distinct primes in $[L,2L]$, the fraction with either exceptional order is $O(d^2\log L/L)$. For fixed $q$, at most $\sum_{e\le d}\varphi(e)\le d(d+1)/2$ residues modulo $q$ have order at most $d$, and each has at most two representatives in the interval. The prime number theorem gives the estimate after counting both directions. It applies even when $A$ depends on $N$. It concerns additional rank differences after the multiplicity check, and leaves individual matrix entries and other functions of $A$ open.

## The first precision that guarantees a split

In $(\mathbb Z/N)[\epsilon]/\epsilon^s$ let $U$ be the size-$d$ cyclic shift whose wraparound coefficient is $\epsilon$. Then $U^d=\epsilon I$, and on all $D=ds$ scalar coordinates it is a single nilpotent shift. With $d\le p<q$ and $p<D\le q$, $(I+U)^N-I$ has local ranks $D-p$ and 0 on those coordinates. Indeed modulo $r$ the first term is $(N/r)U^r$, followed by higher powers of $U^r$. Its cofactor is invertible. The unperturbed block at $\epsilon=0$ powers to $I$ on both sides. This proves sharpness of the uniform $ds$ threshold, including examples where a size-three Jordan block has nonzero differential in characteristic five although its powered matrix is scalar.

One public choice is $D=\lfloor\sqrt N\rfloor+1$, $d=1$, $s=D$. It always satisfies $p<D\le q$. Expanding the positive coefficients of $(1+z)^N-1\bmod z^D$ and taking their gcd with $N$ returns $q$. The check implements this directly with truncated binary powering. Changing the layout to $d$ matrix rows and $s$ coefficient layers stores the same $D$ basis elements. This gives an explicit square-root-sized construction. A method of compressing the needed rank or coefficient remains open. The scalar prefix already has the known identity $N/\gcd(N,\operatorname{lcm}(1,\ldots,D-1))$. [Earlier binomial prefix proof](cartier-paths.md#formal-multiplication-on-the-gaussian-cm-curve)

Run `python experiments/polylog_factoring_2026_11_27/check_matrix_power_jets.py` from the repository root. It verifies 108 local Jacobian images against actual commutator spaces, four explicit higher-order conjugators (including nonsplit eigenvalues), four cases across the sharp precision boundary, and the public shift computation on all 91 products of distinct odd primes through 47.
