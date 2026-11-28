# Equal quaternion frames on two independent lines

Let $N=pq$ for distinct odd primes other than 11. Draw two independent uniform projective lines modulo $N$, and keep the pair when the lines are distinct and have the same full discriminant-11 frame at $M=2^k$. The gcd of their determinant with $N$ gives a proper factor exactly when they share one CRT component. For $R\le p,q\le 2R$ and every $\epsilon>0$, the probability of a proper factor **conditional on keeping the pair** is

$$O_\epsilon(R^{-1+\epsilon}),$$

uniformly in $k$. For $M>2\sqrt{\max(p,q)}$, the number of successful pairs is zero, including when the conditioning event is empty. Thus even direct uniform sampling from the distinct equal-frame pairs, with no cost to find a match, gives at most $O_\epsilon(R^{-1+\epsilon})$ success per pair. This samples a frame with weight proportional to $m_z(m_z-1)$, where $m_z$ is its occupancy. A procedure that deliberately chooses a particular frame or a correlated pair of lines requires a separate analysis. [Public frame classifier](level11-marking.md) · [Earlier unmarked collision law](level11-sampling.md#can-the-ideal-test-distinguish-the-hidden-primes)

For $h$ independent uniform lines, even inspecting every pair gives a baseline success bound of $\binom h2(p+q)/((p+1)(q+1))$ by the union bound. Here $(p+q)/((p+1)(q+1))$ is the exact two-line probability of sharing exactly one component. A polylogarithmic sample therefore has vanishing success before any frame tests. The conditional theorem also bounds the success of $h$ uniform distinct equal-frame pairs by $h\,O_\epsilon(R^{-1+\epsilon})$.

The use of quaternion norm equations and congruences to control short paths goes back to Lubotzky, Phillips and Sarnak. Here the graph has discriminant 11, the step prime is one of the unknown factors, and the count needed for a gcd is the proportion of matches sharing one of the two prime steps. [Lubotzky, Phillips and Sarnak, §3](https://courses.grainger.illinois.edu/cs598sgt/sp2015/LPS.pdf) · [Supersingular graphs with level structure](https://doi.org/10.5802/jtnb.1283)

## Repeated targets in one prime step

Fix a marked source at level $M\ge4$ in the supersingular graph over $\overline{\mathbb F}_{11}$. Its $r+1$ outgoing edges correspond to the order-$r$ cyclic subgroups of an elliptic curve $E$, where $r$ is any prime away from $2$ and $11$. Write $T_r(E,M)$ for the number of **ordered pairs of distinct edges** with the same marked target. A uniform bound over all sources and levels is

$$T_r(E,M)\ll_\epsilon \frac{r^{2+\epsilon}}{M^3},\qquad T_r(E,M)=0\quad\text{if }M>2\sqrt r.$$

Indeed, identify the two targets by an isomorphism preserving their transported frames, and call the resulting quotient isogenies $f_1,f_2:E\to E'$. They agree on $E[M]$. Their composite $\alpha=f_1^\vee f_2$ belongs to $\operatorname{End}(E)$, has norm $r^2$, and agrees with $[r]$ on $E[M]$. Since $[M]$ is the quotient by $E[M]$, we can write $\alpha=r+M\beta$ with $\beta\in\operatorname{End}(E)$. Expanding the norm and using $\gcd(M,r)=1$ gives an integer $u$ satisfying

$$\operatorname{nrd}(\beta)=ru,\quad \operatorname{trd}(\beta)=-Mu,\quad 0\le u\le 4r/M^2.$$

The last inequality follows from $\operatorname{trd}(\beta)^2\le4\operatorname{nrd}(\beta)$ in a definite quaternion algebra. The case $u=0$ gives $\alpha=r$, hence $f_1=f_2$ by composition with $f_1$ and cancellation of $[r]$. This proves the injectivity cutoff.

Each remaining $\alpha$ determines both edges. If $\alpha$ killed $E[r]$, then $\alpha=r v$ for a norm-one unit $v$. The congruence $v\equiv1\bmod M\operatorname{End}(E)$ forces $v=1$ because $\operatorname{nrd}(v-1)\le4<M^2$. Consequently $\alpha$ has rank one on $E[r]$, where it has determinant zero. Its kernel is $\ker f_2$ and its image is $\ker f_1$. The latter follows from $f_1\alpha=[r]f_2$. Distinct ordered pairs therefore give distinct $\alpha$, even if the unmarked target has automorphisms.

For fixed $u>0$, the trace-zero quaternion $w=2\beta+Mu$ lies in the integral trace-zero lattice with

$$\operatorname{nrd}(w)=t=u(4r-M^2u)\le 4r^2/M^2.$$

Here $t>0$, since $M^2$ cannot divide $4r$. The lattice and its ternary norm form $Q$ come from one of only two underlying endomorphism orders at discriminant 11. For any fixed positive definite integral ternary form,

$$r_Q(t)=\#\{w:Q(w)=t\}\ll_{Q,\delta}t^{1/2+\delta}\qquad(t\ge1).$$

To see this, diagonalize over $\mathbb Q$ and clear denominators. Integral solutions inject into those of $AX^2+BY^2+CZ^2=ht$ for fixed positive integers $A,B,C,h$. There are $O(\sqrt t)$ choices of $Z$. For $T=ht-CZ^2>0$, each binary solution produces a distinct algebraic integer $AX+\sqrt{-AB}\,Y$ of norm $AT$ in the fixed imaginary quadratic field. The number of ideals of norm $AT$ is at most the divisor count $\tau(AT)$ by prime-ideal factorization. There are at most six generators per principal ideal. The divisor bound $\tau(AT)\ll_\delta t^\delta$ and the single binary solution for $T=0$ give the claim.

For each $u$, the value of $w$ recovers $\beta$. Summing over at most $4r/M^2$ choices of $u$ bounds the number of $\alpha$ by

$$\ll_\delta \frac r{M^2}\left(\frac{4r^2}{M^2}\right)^{1/2+\delta}\ll_\epsilon \frac{r^{2+\epsilon}}{M^3}.$$

In particular repeated unmarked automorphisms, embeddings of nonmaximal quadratic orders, and prime-dependent choices of source require no extra multiplicity factors in this count.

## Conditioning on an equal frame

There are $D=(p+1)(q+1)$ lines. Their labels lie in a common determinant sector of size $S=5M^3/16$. Fixing the $p$-component leaves a $q$-step from a fixed marked source, since the modifications at $p$ and $q$ commute. Conditional on that component agreeing for the two lines, their ordered $q$-directions are independent uniform choices. The row bound gives

$$\Pr(\text{equal frames, distinct lines, shared component})
\ll_\epsilon \frac1{M^3}\left(\frac{q^\epsilon}{p+1}+\frac{p^\epsilon}{q+1}\right).$$

The second term interchanges the roles of $p$ and $q$. If $m_z$ counts lines with frame $z$, then $\sum_z(m_z/D)^2\ge1/S$ by Cauchy–Schwarz. The probability that the lines are identical is $1/D$. The denominator of the requested conditional probability is therefore at least $1/S-1/D$. In the balanced range, whenever $M\le 2\sqrt{\max(p,q)}$ we have $S/D=O(R^{-1/2})$. Dividing the last two bounds proves the theorem for all sufficiently large $R$. Above the cutoff, both row directions are injective. The finitely many smaller $R$ can be absorbed into the constant whenever the conditioning event has positive probability.

At $M=1,2$ the unmarked or mod-two labels have two or five possible values. The same denominator argument has $S\le5$, and the unconditioned probability of sharing a component is at most $1/(p+1)+1/(q+1)$. This covers these levels without assumptions on marked stabilizers.

An independent, weaker estimate follows from the one-step Ramanujan bound. Within the output determinant sector the squared $\ell^2$ norm of a normalized row law is at most $1/S+(1-1/S)4r/(r+1)^2$. Removing identical edges subtracts $1/(r+1)$. For $S<D$ the conditional factor probability is consequently at most

$$\frac{1/(p+1)+1/(q+1)+6S/D}{1-S/D}.$$

This alone, with the injectivity cutoff, yields $O(R^{-1/2})$ at all levels on balanced inputs. The trace-zero norm count strengthens the exponent by controlling off-diagonal coincidences directly. [Hecke operator bound, Codogni–Lido, Theorem 2.3.8](https://doi.org/10.1016/j.jnt.2026.02.006)

At $N=19\cdot29$ the native classifier gives the following ordered counts among the 600 lines. Sharing a 19-line fixes the 19-component, so the differing edges have degree 29. [Exact fixture and independent frame projection](check_level11_depth.py)

| $M$ | Distinct equal-frame pairs | Share 19-line | Share 29-line | Proper gcd given a distinct match |
|---:|---:|---:|---:|---:|
| 4 | 17,640 | 600 | 600 | $10/147$ |
| 8 | 2,072 | 0 | 96 | $12/259$ |
| 16 | 228 | 0 | 0 | 0 |
| 32 | 44 | 0 | 0 | 0 |
