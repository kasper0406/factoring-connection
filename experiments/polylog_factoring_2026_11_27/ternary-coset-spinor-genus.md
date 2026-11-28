# One proper spinor genus for unit-norm three-square residues

Let $Q(u)=u_1^2+u_2^2+u_3^2$, $L=\mathbb Z^3$, and $m$ be odd. For a unit $b\bmod m$, put
$$
 S_m(b)=\{c\in L/mL:Q(c)\equiv b\pmod m\},\qquad C_c=c+mL .
$$
**All $C_c$, for $c\in S_m(b)$, belong to the same proper genus, and this genus contains exactly one proper spinor genus.** In particular, if $M=y^2+1$ is prime to $m$, the proper spinor genus of every queried coset represents $M$ primitively. Its common genus contains the coset of $(y,1,0)$. The statement includes $M$ twice a square.

## Local equivalence of residues

Write $p^e\Vert m$. Fix an integer $B\equiv b\bmod m$. For an integral lift $a$ of a queried residue, $Q(a)$ is a $p$-adic unit. Multiply $a$ by the square root of $B/Q(a)$ congruent to $1\bmod p^e$, which exists by Hensel's lemma at odd $p$. The result $v$ has exact norm $B$ and still belongs to the same local residue class. Take two such vectors $v,w$. Since $Q(v-w)+Q(v+w)=4B$ is a unit, at least one of these two norms is a unit. The reflection in a vector $z$ of unit norm is integral, with formula $s_z(t)=t-2(t\cdot z)z/Q(z)$.

If $Q(v+w)$ is a unit, $-s_{v+w}$ is a proper integral isometry taking $v$ to $w$. Otherwise $s_{v-w}$ takes $v$ to $w$. Composing with an integral reflection in a unit-norm vector of $w^\perp$ fixes $w$ and gives determinant one. Such a vector exists because $w^\perp\cap L_p$ is a unimodular binary plane. This proves proper local equivalence at each conductor prime. At other finite places the cosets both equal $L_p$.

## Unit spinor norms at every finite place

Write $K_p=\operatorname{SO}(C_{c,p})$, and identify spinor norms with subgroups containing squares. At odd $p\mid m$, choose any integral representative $a$ of $c$. Since $Q(a)$ is a unit, $L_p=\mathbb Z_pa\perp W$ with $W$ unimodular of rank two. Its reduction represents every nonzero element of $\mathbb F_p$. In the split case the binary form is a nonzero multiple of $xy$, and in the nonsplit case it is a multiple of the surjective quadratic-field norm. Each such representation is nonsingular and lifts to any prescribed $p$-adic unit with that residue. For $\epsilon\in\mathbb Z_p^\times$, choose vectors in $W$ of norms $1$ and $\epsilon$. Their two reflections fix $a$, and their product lies in $K_p$ with spinor norm $\epsilon$. At odd $p\nmid m$, $C_{c,p}=L_p$ and the same argument applies with $a=(1,0,0)$.

At $2$, oddness of $m$ gives $C_{c,2}=L_2$. Integral reflections in $(1,0,0),(1,1,1),(1,2,0)$ have spinor norms $1,3,5$. Products of pairs give proper rotations of spinor norms $3$ and $5$, which generate $\mathbb Z_2^\times/\mathbb Z_2^{\times2}$.

Kane–Kim's coset form of the spinor-genus index is
$$
\#\{\text{proper spinor genera in }\operatorname{gen}^+(C_c)\}
=[\mathbb I_{\mathbb Q}:\mathbb Q^\times\textstyle\prod_v\theta(K_v)].
$$
All finite local unit groups are in the denominator. For any idèle, divide by a rational with the prescribed finite valuations and the same real sign. Its remaining finite components are units and its real component is positive, hence belongs to the real spinor norm group modulo squares. The index is one. The condition that $b$ be a unit also gives full conductor $m$, as used in the cited coset convention. [Kane–Kim, Proposition 2.1](https://arxiv.org/pdf/2211.03987#page=5)

## The exact common representation counts

The underlying positive definite unimodular lattice has one proper class in its genus. Indeed, Minkowski's convex body theorem in rank three gives a nonzero vector of squared norm less than $2$. It has norm $1$ and splits off an orthogonal summand. Its unimodular rank-two complement has a norm-one vector by the same argument, giving an orthonormal integral basis. A basis sign change adjusts orientation. For any coset in the genus, its ambient lattice is its difference lattice divided by $m$, so it may be identified with $L$. In $L$ its shift is integral and has norm $b\bmod m$, by local equivalence. Thus the proper classes in the genus are exactly the $\operatorname{SO}(L)$-orbits on $S_m(b)$.

Let $r_3(n)$ count all ordered signed representations and $r_3^*(n)$ those with coprime coordinates, for $n>0$. Weight each proper class by the reciprocal of its proper automorphism-group order, and divide by the sum of the weights. By orbit–stabilizer, this is the uniform average over the residues $c\in S_m(b)$. Each integral representation of $n\equiv b\bmod m$ contributes in precisely its own residue. Hence
$$
 r_{\operatorname{spn}^+(C_c)}^{(*)}(n)=r_{\operatorname{gen}^+(C_c)}^{(*)}(n)
 =\begin{cases}r_3^{(*)}(n)/|S_m(b)|,& n\equiv b\bmod m,\\
 0,& n\not\equiv b\bmod m .
 \end{cases}
$$
The same formula holds with or without the primitive restriction. Its denominator is
$$
 |S_m(b)|=m^2\prod_{p\mid m}\left(1+\frac{(\frac{-b}{p})}{p}\right).
$$
For an odd prime $p$, the number of $(u,v)$ with $u^2+v^2=t$ is $p-(\frac{-1}{p})$ at $t\ne0$ and $p+(p-1)(\frac{-1}{p})$ at $t=0$. Summing over $t=b-w^2$ gives $p^2+p(\frac{-b}{p})$. Each solution has nonzero gradient and has $p^2$ lifts per higher power of $p$. The Chinese remainder theorem gives the displayed product.

For a fixed $m$ and $M$, even exact primitive genus and spinor-genus counts give identical answers across the queried residues. Information in the common number $r_3^*(M)$ and differences between representation counts of individual proper classes are separate questions.

## Where the difference occurs in the theta series

The repository already investigated subtracting theta series. The October level-11 computation expresses a fixed elliptic newform as one quarter of the difference of two quaternary lattice counts. Phase 30 then studied what a coefficient of a fixed non-CM form at $N=pq$ reveals. Here the form, the modulus and the index vary with the input. The statement specific to these ternary cosets is that the Eisenstein parts agree for every queried residue at the same modulus and every unary component vanishes. The October single-coefficient question remains open. [Earlier level-11 calculation, §4](../../notes/nativemod_results.md) · [October 28 fixed-coefficient investigation](../../notes/cuspcoef_results.md)

Kane and Kim's Theorem 1.2 decomposes the theta series of an individual ternary lattice coset into an Eisenstein series, a linear combination of unary theta series, and a cusp form orthogonal to the unary theta series. The first is the genus average. The second is the spinor-genus average minus the genus average and vanishes here. For two residues at the same unit norm and odd modulus, the difference of their full theta series therefore belongs entirely to the third space. This conclusion concerns the ordinary series counting all integral representations. It gives no running-time estimate for computing a coefficient in this space. [Kane and Kim, Theorem 1.2](https://arxiv.org/pdf/2211.03987#page=2)

For example, take the packed sphere for $N=77,x=1$, with squared norm $M=34\,179\,874\,885$ and modulus $m=40\,817=7^4\cdot17$. The class $(21\,609,28\,813,33\,615)$ contains exactly one representation, from a noncoordinate frame. The class $(31\,213,2\,402,19\,209)$ contains none. The difference of their $M$-th coefficients is exactly one. Both give the same Eisenstein series and have zero unary part.

The check constructs proper rotations of square spinor norm at both prime powers taking each submitted residue to the public class $(\sqrt{M-1},1,0)$. It tests the determinant, the orthogonality of the columns and the image of the residue. Across the six exact-existence fixtures it did this for all 1,544 queries, including 456 occupied classes and 1,088 empty ones, counted independently by an integer square test. [Construct the rotations and check the integer counts](check_ternary_coset_spinor.py) · [Factoring reduction and full fixture table](packed-ternary-decision.md)
