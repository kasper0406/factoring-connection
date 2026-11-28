# Short tangent-line intervals for cofactor bit tests

For every fixed $0<\delta<1/4$ and all sufficiently large odd $N$, the balanced floor-quotient test below admits every locally possible answer to every singleton or pair of bit positions of period $P\le N^{1/4-\delta}$. The statement holds for each such $N$, including every balanced odd semiprime. It concerns the answers to individual bit queries. The witnesses can change with the query and the answer.

Number bits from zero and select one or two distinct positions, all at least 1. Let $M$ be their bit mask, let $k$ be one plus the highest position, and write $P=2^k$. A **positive** answer means each of the two integers has at least one selected bit set. A **locally possible** answer occurs for some odd $x,y\bmod P$ with $xy\equiv N\bmod P$. For each such answer the theorem gives odd integers $u,v$ attaining it with

$$
 u<v<2u,\qquad v=\lfloor N/u\rfloor,\qquad uv\equiv N\pmod P.
$$

The actual hidden divisors on a balanced semiprime satisfy all three conditions. The singleton negative outcome is forced locally precisely at $N\equiv3\bmod4$ for $P=4$, $N\equiv5,7\bmod8$ for $P=8$, $N\equiv11,13\bmod16$ for $P=16$, and $N\equiv19,29\bmod32$ for $P=32$. Every other singleton and every pair admit both answers locally. The preceding result also gives a lift for **each prescribed pair** of product residues through $P\ll N^{1/6}$. [Local classification and prescribed-residue lift](cofactor-bit-completions.md)

## A uniform short-interval estimate

Fix $0<\epsilon<1$, an odd integer $n$, and an interval $I$ of $R$ consecutive integers, $P^\epsilon\le R\le P$. For every odd $L$, let $x=x(L),y=y(L)$ be the unique odd residues modulo $P$ with

$$
xy=n,\qquad 3x+2y=L\pmod P.
$$

Indeed $f(x)=3x+2n/x$ permutes the odd residues, since $f(x)-f(z)=(x-z)(3-2n/(xz))$ and the second factor is odd. Write $e_P(s)=\exp(2\pi i s/P)$. If each of $h,j$ is either zero or odd, and at least one is odd, then for some $\sigma=\sigma(\epsilon)>0$,

$$
 \left|\sum_{\substack{L\in I\\ L\ {\rm odd}}}
     e_P\bigl(hx(L)+jy(L)\bigr)\right|
 \ll_\epsilon R P^{-\sigma}.
 \tag{1}
$$

The bound is uniform in $n,h,j$ and the location of $I$. It suffices to prove it for all sufficiently large $k$ in terms of $\epsilon$.

**Local power series.** Work in the 2-adic integers. The equation $Ly=3n+2y^2$ has the odd solution

$$
 y(L)=\frac{3n}{L}\sum_{\nu\ge0} C_\nu
             \left(\frac{6n}{L^2}\right)^\nu ,
 \qquad x(L)=\frac{L-2y(L)}3,
$$

where $C_\nu$ are the integral Catalan numbers. Their generating series satisfies $C(z)=1+zC(z)^2$, and the argument $6n/L^2$ is even, so substitution converges. Reduction gives the unique pair at every finite precision.

At any odd $c$, expand $y(c+s)=\sum_{m\ge0}a_m s^m$ for even $s$. **Every $a_m$ is odd.** The term with $\nu=0$ contributes $3n(-1)^m c^{-m-1}$. For $\nu\ge1$ the expansion of $(c+s)^{-2\nu-1}$ has 2-adically integral coefficients and the multiplier contains $2^\nu$. This also justifies summing the expansions coefficientwise. In the phase

$$
 G(L)=h x(L)+j y(L)=hL/3+(j-2h/3)y(L)
$$

each coefficient of degree $m\ge2$ therefore has valuation $v=v_2(3j-2h)$. Here $v=0$ for odd $j$, and $v=1$ for $j=0$ and odd $h$.

**Polynomial reduction.** Take the fixed integer $D=\lceil 2/\epsilon\rceil$ and $a=\lceil k/(D+1)\rceil$. Partition the odd integers in $I$ into odd residue classes $c\bmod 2^a$. On writing $L=c+2^a t$, degrees greater than $D$ in $G$ disappear modulo $P$. Choose ordinary integer representatives for the remaining 2-adic coefficients modulo $P$. This gives an ordinary real polynomial phase $G(c+2^a t)/P$ modulo one whose degree-$D$ coefficient has exact reduced denominator

$$
 Q=2^{k-Da-v}\asymp_D P^{1/(D+1)}
$$

for sufficiently large $k$. Each resulting $t$ interval has length $H\asymp R/2^a$. Translating it preserves the leading coefficient.

For fixed $\eta>0$, the classical Weyl inequality for this rational leading coefficient bounds each sum by

$$
 \ll_{D,\eta} H^{1+\eta}
       (Q^{-1}+H^{-1}+QH^{-D})^{\,2^{1-D}}.
$$

The three terms are respectively $O_D(P^{-1/(D+1)})$, $O_D(P^{1/(D+1)-\epsilon})$, and $O_D(P^{1-D\epsilon})$. All exponents are negative. Taking $\eta$ sufficiently small in terms of $\epsilon$ gives a fixed power saving relative to $H$. Summing over classes proves (1). [Weyl's original paper](https://doi.org/10.1007/BF01475864) · [Trevor Wooley, *Translation invariance, exponential sums, and Waring’s problem*, Theorem 7.1 with the displayed bound](https://www.math.purdue.edu/~twooley/publ/2014icm.pdf#page=15)

## Detecting both answers

Write $A_P(z)$ for the indicator that the representative $0\le z<P$ is in the upper half. Its mean is $1/2$. In the normalized additive Fourier expansion, $A_P-1/2$ has only odd frequencies by translation by $P/2$. The geometric progression formula bounds the absolute value at $1\le h<P$ by $O(1/\min(h,P-h))$, so the coefficients have absolute sum $O(\log P)$. The lower-half indicator $B_P=1-A_P$ has the same properties.

Let $T$ count odd $L\in I$. A sufficient positive witness has $A_P(x)A_P(y)=1$, so its highest selected bit is set in both coordinates. Expand the product in Fourier series. The constant term gives $T/4$, and (1) bounds the sum over all remaining terms by $O_\epsilon(RP^{-\sigma}\log^2 P)$. Thus these witnesses occur in every such interval for sufficiently large $P$.

A sufficient negative witness has every selected bit of $x$ clear. For a singleton its count is $T/2+O_\epsilon(RP^{-\sigma}\log P)$. For a pair let $P_0<P$ be the smaller bit's period and expand the sufficient witness indicator as

$$
 B_P(x)B_{P_0}(x)
   =(B_P(x)-1/2)B_{P_0}(x)+\tfrac12 B_{P_0}(x).
$$

On modulus $P$, the first product has only odd frequencies since all frequencies lifted from $P_0$ are multiples of the even integer $P/P_0$. Convolution bounds their absolute sum by $O(\log^2 P)$. Equation (1) with $j=0$ gives $O_\epsilon(RP^{-\sigma}\log^2 P)$ for its contribution.

For the last term, $x(L)\bmod P_0$ depends only on $L\bmod P_0$. In every complete $P_0$ cycle the count of odd $L$ with $B_{P_0}(x)=1$ is exactly half the number of odd $L$, because $P_0\ge4$ and the inverse line map permutes the odd residues. If $P_0\le\sqrt R$, the incomplete cycle costs $O(\sqrt R)$. For $P_0>\sqrt R$, an incomplete cycle of length less than $P_0^{\epsilon/2}$ also costs $O(\sqrt R)$ since $P_0^\epsilon\le P^\epsilon\le R$. For a longer incomplete cycle apply (1) and the half-interval Fourier expansion at modulus $P_0$ with exponent $\epsilon/2$. The resulting error is $O_\epsilon(RP_0^{-\sigma(\epsilon/2)}\log P_0)=o_\epsilon(R)$ uniformly as $P\to\infty$. This gives $T/4+o_\epsilon(R)$ for the sufficient two-bit negative witnesses, uniformly in the lower bit.

## From residues to floor quotients

Take $n=N$ and $I$ to start at $\lceil\sqrt{24N}\rceil$. Use $R=\lceil P^\delta\rceil$ and $\epsilon=\delta$. For large $P$ the preceding result supplies each answer. For $P$ bounded in terms of $\delta$, use an interval of $P$ consecutive integers. It contains a line for every locally possible pair since $f$ permutes the odd residues. All estimates below also hold in this bounded case with $R=P$.

For a line $L$ giving the desired answer, set $d=\sqrt{L^2-24N}$ and $t=(L-d)/6$. There is one class $u\bmod 2P$ with $u\equiv x\bmod P$ and $3u\equiv L-2y\bmod 2P$. Indeed the latter congruence has one solution modulo $2P$ and reduces to $u\equiv x\bmod P$. Choose the largest $u\le t$ in that class and put $s=t-u$, $v=(L-3u)/2$. Then $0\le s<2P$, $v\equiv y\bmod P$, and expansion about the lower root gives the exact identity

$$
 N-uv=\tfrac12(ds+3s^2).
$$

Since $L=\sqrt{24N}+O(R)$ and $R=o(\sqrt N)$, we have $d=O(N^{1/4}\sqrt R)$ and

$$
 0\le N-uv\ll P N^{1/4}\sqrt R+P^2=o(\sqrt N).
$$

For the last step, $P\sqrt R\ll N^{(1/4-\delta)(1+\delta/2)}=o(N^{1/4})$ uniformly over growing $P\le N^{1/4-\delta}$. Also $u\sim\sqrt{2N/3}$ and $v/u\to3/2$ uniformly. Consequently $0\le N-uv<u$ and $u<v<2u$ for all sufficiently large $N$. The congruence $uv\equiv N\bmod P$ and the selected bits follow from the chosen residues.

The result leaves open periods at the quarter-power endpoint, windows of polylogarithmic length, and algorithms that use the exact divisibility condition $u\mid N$. For fixed $\epsilon$ the proof permits scanning $O(P^\epsilon)$ lines, with a saving exponent in (1) that depends on $\epsilon$. [Integer witness construction](cofactor_tangent_lifts.py)

As a finite check on the interval calculation, enumerate all odd $x\bmod P$, put $y=n/x$ and order the pairs cyclically by $3x+2y\bmod P$. Over $n\in\{1,3,P/2-1,P-1\}$, the longest run of successive odd lines missing the sufficient positive witness was 46, 25, and 43 for $k=8,12,16$ respectively. For the sufficient negative witness, maximizing also over singleton and pair masks with top position $k-1$, the corresponding run lengths were 20, 38, and 41.
