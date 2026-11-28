# Two neighboring line ideals suffice at discriminant 11

Let $\mathcal O=\mathbb Z[\omega]+\mathbb Z[\omega]j$ be the discriminant-11 order, with $j^2=-1$ and $(1+j)\overline{(1+j)}=2$. Fix $n$ coprime to 22 and a matrix splitting modulo $n$. Write $I_n(L)$ for the right ideal whose matrices have both columns in the projective line $L$. Let $C_n(L)$ indicate that it has a generator of norm $n$. Given **any** line, one of $L$ and $(1+j)L$ supplies a generator. `repair_line` tests the first ideal and tests the second only if needed. Both shell tests take polynomial bit time, including for prime powers. Returning the first generator found makes this a deterministic construction after the splitting has been supplied. The earlier sampler gives the uniform distribution when that is required. [Repair](level11_neighbors.py) · [Splitting and shell test](level11-sampling.md)

## Why the second test always works

For $g=1+j$, multiplication on the left gives $g I_n(L)\subset I_n(gL)$. The latter ideal has index $n^2$ in $\mathcal O$, and $g I_n(L)$ has index $4n^2$. The inclusion is a norm-two neighbor relation, with source class that of $I_n(gL)$ and target class that of $I_n(L)$. The Brandt matrix at 2, with the four-unit principal class first and the six-unit class second, is

$$B_2=\begin{pmatrix}1&2\\3&0\end{pmatrix}.$$

Indeed $\mathcal O$ has four elements of norm 2, giving one principal subideal of the principal class. Both rows sum to three and detailed balance uses stationary masses $3/5,2/5$. The zero entry excludes the case where both ideals are nonprincipal. This uses **one fixed** transformation of the supplied line. The three outgoing Brandt edges count all norm-two subideals of an ideal class. Also $g^2=2j$, so $g^2$ acts on projective lines as the global unit $j$ and $C_n(g^2L)=C_n(L)$. The entire orbit of class bits has period dividing two. [Level-11 Brandt matrices](../../notes/nativemod_results.md#4-q4--brandt-matrices-and-quaternary-theta-series)

For a balanced or unbalanced $n=pq$ of two distinct primes away from 2 and 11, put $\pi=3/5$, $\lambda_q=a_q/(q+1)$ for the level-11 weight-two coefficient, and fix the component $L_p=l$ of a uniform line. Define $A=\pi+\lambda_q(C_p(l)-\pi)$ and $B=\pi+\lambda_q(C_p(gl)-\pi)$. These are the two individual conditional means, since multiplication by $g$ permutes the $q$-lines. The support condition above gives the entire conditional law.

| $(C_n(L),C_n(gL))$ | $(0,0)$ | $(0,1)$ | $(1,0)$ | $(1,1)$ |
|---|---:|---:|---:|---:|
| Probability given $L_p=l$ | $0$ | $1-A$ | $1-B$ | $A+B-1$ |

When $a_q=0$, these probabilities are exactly $0,2/5,2/5,1/5$ on every row. For $551=19\cdot29$ this gives counts $0,12,12,6$ in every row of length 30. On this family even an arbitrary Boolean function of the two bits has the same conditional law on all those rows. [Original marginal identity](level11-sampling.md#can-the-ideal-test-distinguish-the-hidden-primes)

## What several quaternion translates can reveal about one hidden component

Take primitive integral quaternions $g_i$ of positive norms $m_i$ coprime to $11pq$. Include the identity if the original bit is used. Their action must come from these rational quaternions in the chosen splitting. The assertion below does not cover arbitrary matrices chosen modulo $pq$. Write $K=\widehat{\mathcal O}^{\times}$, and let $h_i$ be the finite idele equal to $g_i^{-1}$ at primes dividing $m_i$ and 1 elsewhere. Form the compact open subgroup and its index

$$K'=K\cap\bigcap_i h_i K h_i^{-1},\qquad I=[K:K'].$$

For any Boolean function of $(C_{pq}(g_iL))_i$, the absolute difference between its probability conditional on a particular $p$-line and its unconditional probability is at most

$$\sqrt{5I/2}\left(\frac{\sqrt q}{q+1}+\frac{2\sqrt{pq}}{(p+1)(q+1)}\right).$$

Here $I\le\prod_i\psi(m_i)$, where $\psi(m)=m\prod_{r\mid m}(1+1/r)$. If all $m_i$ divide a common $M$ coprime to $11pq$, then $I\le |\operatorname{PGL}_2(\mathbb Z/M)|\le M^3$, independently of the number of translates. Both estimates follow locally from the Smith form of a primitive matrix of determinant $m_i$. At $r^e\parallel m_i$ the intersection of its two maximal compact groups stabilizes a vertex at distance $e$ from the base vertex in the lattice tree. Its index is $r^e+r^{e-1}$. The subgroup of units scalar modulo $M$ stabilizes every such lattice. Thus for balanced factors of scale $R$ and $M=R^{o(1)}$ the change in any row is at most $R^{-1/2+o(1)}$. Adaptive choices from the same fixed collection are included if the output uses only the observed bits.

**Reduction to a Hecke step.** Let $B$ now denote the rational quaternion algebra, and $X=B^\times(\mathbb Q)\backslash B^\times(\mathbb A_f)/K'$. Each map $[x]\mapsto[xh_i]_K$ to the two ordinary ideal classes is well-defined. Represent $I_{pq}(L)$ by local matrices $x$ supported at $p,q$. Multiplying $I_{pq}(g_iL)$ on the left by the rational $g_i^{-1}$ recovers $I_{pq}(L)$ locally at $p,q$ and gives the lattice $g_i^{-1}\mathcal O_r$ at $r\mid m_i$. Thus the lifted state $[x]$ determines every requested bit. Let $F$ be their Boolean function on $X$. For fixed $L_p=l$, the conditional mean is $(P_q F)([x_l])$, where $P_q$ averages the $q+1$ Hecke neighbors. Globally the mean is $(P_pP_q F)([1])$. These are exact identities because $K'$ is maximal at $p,q$.

**Determinant classes and mixing.** Give $X$ its mass measure $\mu(x)$ proportional to $1/|\Gamma_x|$, where $\Gamma_x=B^\times(\mathbb Q)\cap xK'x^{-1}$. Every stabilizer has size at most 6. The unnormalized mass is $I(1/4+1/6)$, by lifting the two ordinary classes to level $K'$. Thus $\mu(x)\ge2/(5I)$. Let $E$ average over each reduced-norm class in $\mathbb Q_{>0}^{\times}\backslash\mathbb A_f^\times/\operatorname{nrd}(K')$. All $r$-neighbors of a vertex lie in the same next norm class. The weight-two Ramanujan bound on the orthogonal complement of the functions of this class gives

$$\|P_r(F-EF)\|_2\le\frac{2\sqrt r}{r+1}\|F-EF\|_2\qquad(r=p,q).$$

Under the ideal-class description of supersingular curves with level structure this is the bound on the complement of the Weil-pairing classes. The level subgroup here contains all scalars and a principal congruence subgroup. The weighted adjacency operator is normal in general, and it preserves this orthogonal complement. These facts also cover levels whose determinant subgroup is proper. [Codogni and Lido, §§2.2–2.3, especially Theorem 2.3.8](https://doi.org/10.1016/j.jnt.2026.02.006)

The $EF$ terms in $(P_qF)([x_l])$ and $(P_pP_qF)([1])$ are identical, since both end in the same norm class. The Boolean condition gives $\|F-EF\|_2\le1/2$. Bound evaluation at one vertex by $\mu(x)^{-1/2}$ times the $L^2$ norm, and apply the spectral estimate once or twice, respectively. This gives the stated bound without treating distinct translates as independent.

**Computing additional bits from one generator.** Let $\alpha$ generate $I_n(L_0)$, and let $h\in\mathcal O$ have norm $m$ coprime to $n$. Then $C_n(hL_0)=1$ exactly when there is a $\zeta\in\mathcal O$ of norm $m$ with $h\alpha\bar\zeta\in m\mathcal O$. A generator $\beta$ of the target gives $h\alpha=\beta\zeta$. Conversely $\beta=h\alpha\bar\zeta/m$ is integral of norm $n$, and its line at $n$ is $hL_0$ because $\zeta$ is invertible modulo $n$. The congruence needs only $\alpha\bmod m$ and the precomputable norm-$m$ shell. If the repair selected $L_0=g_0L$ with $g_0=1$ or $1+j$, use $h=g_i\bar g_0$ to recover the bit at $g_iL$. In a collection with all $m_i\mid M$, the one generator modulo $2M$ and the choice of $g_0$ suffice for all these tests. Enumerating the small shells costs polynomial time in $M$ and hence polynomial bit time when $M$ is polylogarithmic in $n$.

The finite diagnostic compares the repair with exhaustive norm shells and tests the two conditional rows in both CRT directions. The new output is a guaranteed generator after at most two tests and a bound for these correlated probes. A method for splitting arbitrary balanced semiprimes remains open. [Runnable checks](check_level11_neighbors.py)
