# Schur compression of a free projector

Let $B=\mathbb Z/N$ with $N\ge2$. Choose $E$ uniformly among the conjugates of $\operatorname{diag}(I_k,0)$ under $\mathrm{GL}_d(B)$, where $0<k<d$. Deleting one coordinate gives two natural projectors in dimension $d-1$. When $d=2k$, their ranks are the two closest integers to $(d-1)/2$. The construction and its complete conditional law use a single sample. [Exact diagnostic](check_projector_compression.py)

## The two endpoints

Write $E=\left(\begin{smallmatrix}A&u\\v&e\end{smallmatrix}\right)$ and suppose $e(1-e)$ is a unit. The block equations in $E^2=E$ give $A^2+uv=A$, $Au=(1-e)u$, $vA=(1-e)v$, and $vu=e(1-e)$. Define

$$
P_-=A-\frac{uv}{e},\qquad
P_+=A+\frac{uv}{1-e},\qquad
R=P_+-P_-=\frac{uv}{e(1-e)}.
$$

Here $R=u\ell$ with $\ell=v/(e(1-e))$ and $\ell u=1$. Thus $R$ is an idempotent with free rank-one image and trace one. The block equations give $P_-^2=P_-$ and $RP_-=P_-R=0$. Consequently $P_+=P_-+R$ is also idempotent. Eliminating the last row and column of $E$ with pivot $e$ gives a matrix equivalent to $\operatorname{diag}(P_-,e)$. Reduction at each prime then gives rank $k-1$ for $P_-$ and rank $k$ for $P_+$. Their images are free over $B$ because their ranks are constant across its prime-power components.

Conditional on the good pivot, $(P_-,R,I-P_+)$ is exactly uniform among ordered decompositions of $B^{d-1}$ of free ranks $(k-1,1,d-k-1)$. This holds even conditional on any particular good value of $e$. To count its fiber, choose any generator $u$ of the given rank-one image and write $R=u\ell$ for the unique such row $\ell$. Set $A=P_-+(1-e)R$ and $v=e(1-e)\ell$. These choices reconstruct an idempotent in the prescribed rank-$k$ conjugacy class. They give all preimages, exactly $\varphi(N)$ for every pair and every good $e$. Conjugacy also shows that either endpoint alone is uniform in its fixed free-rank orbit.

## How often the pivot works

Over $\mathbb F_p$, write $E=G\operatorname{diag}(I_k,0)G^{-1}$, $x=G^{-1}e_d$, and $\lambda=e_d^\top G$. The vector $x$ is uniform nonzero, and conditional on $x$ the functional $\lambda$ is uniform subject to $\lambda(x)=1$. Write $x=x_1+x_0$ along the two coordinate blocks. The pivot is $e=\lambda(x_1)$. If $x_1=0$ or $x_0=0$ it equals zero or one, respectively. Otherwise $x_1,x_0$ are linearly independent and $e$ is uniform on $\mathbb F_p$. Hence

$$
\Pr[p\nmid e(1-e)]
=\frac{(p^k-1)(p^{d-k}-1)}{p^d-1}\,\frac{p-2}{p}.
$$

The factor depends only on the reduction modulo $p$. Reducing a uniform basis at a prime power gives a uniform basis over $\mathbb F_p$, and the reductions at distinct primes are independent. Multiply these factors over the primes dividing $N$ for the global probability. At $p=2$ the good event is empty. The local failure probability is at most $4/p$ since each event $x_1=0$ or $x_0=0$ has probability at most $1/p$. A bad pivot may yield a proper divisor via either $\gcd(e,N)$ or $\gcd(1-e,N)$. Their individual gcds can split $N$ even if $\gcd(e(1-e),N)=N$.

## Interpolating a mixed rank

For $P(t)=P_-+tR$, orthogonality gives $P(t)^2-P(t)=(t^2-t)R$. Taking traces shows this vanishes exactly when $t$ is itself a scalar idempotent modulo $N$. Moreover $t=\operatorname{tr}(P(t)-P_-)\bmod N$. The global choices $0,1$ give free endpoint ranks. A mixed choice has $t=0$ on some prime-power components and $t=1$ on the rest, so $\gcd(t,N)$ already splits the input.

At $N=15$, a rank-two input in dimension four compresses to $P(t)=\operatorname{diag}(1,t,0)$. The choices $t=6,10$ give ranks $(1,2)$ and $(2,1)$ modulo $(3,5)$, respectively. The public choices $t=0,1$ give ranks $(1,1)$ and $(2,2)$. Exhausting all 10,530 rank-two projectors in dimension four over $\mathbb F_3$ gives 2,808 good pivots and the predicted fiber sizes. The diagnostic also covers modulus 9 and a squarefree composite. Run `python experiments/polylog_factoring_2026_11_27/check_projector_compression.py` from the repository root. [Reproduce the counts](check_projector_compression.py)
