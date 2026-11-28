# Additive Fourier mass of shifted cubic Jacobi readings

Let $N=pq$ for distinct primes greater than three, and put $\omega^2+\omega+1=0$. At each $x\bmod N$, read the cubic Jacobi symbols of

$$
 z(x+a_i)=\frac{x+a_i-\omega}{x+a_i-\bar\omega},\qquad 1\le i\le J,
$$

at fixed public shifts. These symbols can be evaluated by Eisenstein reciprocity. When all norm denominators are units, write the label vector in $G=(\mathbb Z/3)^J$ using $\omega$ as generator. Given any complex function $F$ on $G$, form a scalar state by applying $F$ to the vector. Give the state amplitude zero at inputs with any singular denominator. A proper gcd found at such an input already splits $N$. [Public implementation](eisenstein_torus.py) · [Bach and Sandlund on cubic Jacobi evaluation](https://arxiv.org/abs/1807.07719)

On balanced inputs, if $3^J J^2\le c\sqrt N$ for a sufficiently small absolute $c$, the fraction of the state's Fourier mass at nonzero frequencies divisible by $p$ or $q$ is $O(J/\sqrt p+J/\sqrt q)$. This holds for **every** $F$ producing a nonzero state, including indicators of rare label patterns. In the other direction, $(2+\epsilon)\log_3N$ independent uniform shifts give an injective signature on all $x\bmod N$ with probability tending to one. At that length the existence of a suitable table and the cost of constructing or evaluating it are separate issues.

The Fourier transform here is on the additive parameter $x\bmod N$. Earlier work treated multiplicative orbits, shifted quadratic symbols on root labels, and a pair of norm-Jacobi probes for another decision problem. A state supported on observable nonunits has a different normalization and already yields factors when sampled in the input basis. [Multiplicative-orbit audit](../../notes/periodfind_results.md) · [Earlier norm-Jacobi probe](../ellipse_disk_mod8/cubic59.md) · [Nonunit-state calibration](nonunit-fourier-sampling.md)

## What survives averaging at a prime?

For $r=p,q$ let $T_r$ be the norm-one torus over $\mathbb F_r$. Its cyclic group of rational points has order $r-1$ in the split case and $r+1$ in the inert case. The cubic Jacobi symbol at the ideal $(r)$ restricts to a nontrivial cubic character on this group. At a split prime the two conjugate prime ideals give equal cubic symbols on norm-one elements, whose product is their square. At an inert prime the exponent $(r^2-1)/3$ restricts to $(r-1)(r+1)/3$ on $T_r$, where $r-1\equiv1\bmod 3$. The global label is therefore the sum of the two local labels.

For one shift the normalized sum over finite $x\bmod r$ of this cubic character, with zero at exceptional $x$, is exactly $-1/r$. The projective parameter gives every torus element once, the parameter infinity gives the identity, and the full nontrivial character sum is zero. Write $t_r$ for the fraction of valid parameters. For the single complex-valued symbol, the normalized mass at frequencies divisible by $q$ but not $p$ is exactly

$$
 \frac{(t_p-p^{-2})q^{-2}}{t_pt_q}.
$$

Indeed this frequency set measures the variance over $x_p$ of the average over $x_q$. The complex cubic symbol is a pure product of local functions. Interchanging $p,q$ gives the other part of the mass.

For several shifts, remove identical shifts modulo $N$. If two remaining shifts coincide modulo one prime, a gcd of their difference already recovers a factor. Assume the shifts are therefore distinct modulo each prime and $J<\min(p,q)$. For a nonzero vector $e\in G$ consider the rational torus function $f_e=\prod_i z(X+a_i)^{e_i}$, interpreting exponents as integers. Over an algebraic closure its divisor modulo three consists of edges between $\bar\omega-a_i$ and $\omega-a_i$, with common displacement $\Delta=\omega-\bar\omega$. A cube would require the function assigning weight $e_i$ at the starting point of each edge, and zero elsewhere, to be invariant under translation by $\Delta$. Every nonempty orbit of that translation has $r$ points. At most $J<r$ points have nonzero weight, so invariance is impossible. This handles coincidences of a numerator root with a denominator root, even when $(a_i-a_j)^2+3$ vanishes at a prime.

Let $\pi_r$ be the law of the local label vector conditional on a valid input and set

$$
 b_r=\max_{e\ne0}\left|\mathbb E_{\pi_r}\omega^{e\cdot U}\right|.
$$

We have $t_r\ge1-2J/r$ and $b_r=O(J/\sqrt r)$ for $J\le r/4$. Here is a proof including the inert torus. For each of three representatives $v$ of $T_r/T_r^3$, take the smooth projective model of $y^3=v f_e(X)$ in the torus. It is an $\mathbb F_r$-curve, geometrically connected by the divisor argument. In geometric multiplicative coordinates its degree-three map has at most $2J$ ramified points, all totally ramified. Thus Riemann–Hurwitz gives genus at most $2J-2$. Hasse–Weil bounds the point count of each twist by $r+1+O(J\sqrt r)$. Removing points above infinity or any of the excluded zeros and poles costs $O(J)$. At every remaining $X$ there are three rational points on exactly the twist for its cubic class, since $T_r$ has three rational cube roots of unity. Subtracting these class counts with weights $1,\omega,\omega^2$ cancels the main term and proves the bias bound. Estimates of this type for norm groups have established prior literature. [Wen-Ching Winnie Li, *Character sums over norm groups*](https://doi.org/10.1016/j.ffa.2004.06.006)

## Every nonlinear table in the small-menu range

Set $K=3^J$ and $\eta=(K-1)b_pb_q$. For $\eta<1$ the fraction of nonzero Fourier mass on the $q$-multiple axis, divided by the **total** Fourier mass including zero, is at most

$$
 (1-t_p)+\frac{(1+(K-1)b_p)b_q^2}{1-\eta}. \tag{1}
$$

For the proof use the uniform probability measure on $G$. Fourier inversion gives the density of $\pi_p$ at most $1+(K-1)b_p$. Independence of the CRT inputs makes the joint label law on valid parameters $\pi_p*\pi_q$. Its density is at least $1-\eta$ since each of its nonconstant Fourier coefficients is a product of two biases. In particular every nonzero $F$ has mean-square on this law at least $(1-\eta)\|F\|_2^2$.

Let $B(u)=\mathbb E_{\pi_q}F(u+V)$. Under the uniform measure its mean is the uniform mean $\bar F$ of $F$. Parseval gives $\|B-\bar F\|_2^2\le b_q^2\|F\|_2^2$. The stated upper density of $\pi_p$ therefore bounds $\operatorname{Var}_{\pi_p}B$ by $(1+(K-1)b_p)b_q^2\|F\|_2^2$. On averaging the actual state over $x_q$, the function of $x_p$ obtained is zero at invalid parameters and $t_q B(U)$ otherwise. Its variance has contributions $t_pt_q^2\operatorname{Var}_{\pi_p}B$ and $t_p(1-t_p)t_q^2|\mathbb E_{\pi_p}B|^2$. Divide by the state norm $t_pt_q\mathbb E_{\pi_p*\pi_q}|F|^2$. Jensen's inequality bounds the second ratio by $1-t_p$, proving (1). The $p$-multiple frequencies give the same expression with the indices exchanged.

When $\eta\le1/2$, the terms $(K-1)b_pb_q^2$ and $(K-1)b_qb_p^2$ are bounded by $b_q/2$ and $b_p/2$. The claimed small-menu bound follows. It includes arbitrary nonlinear postprocessing of the valid cubic symbols. Assigning nonzero amplitude to a singularity flag permits a state supported on nonunits and is outside (1).

## Can the cubic readings predict the scaled quarter-window bit?

The repository's quarter-window reduction factors a modulus with two distinct prime divisors given a predictor for $g_0(ux)$ with inverse-polynomial advantage on uniform **unit** scales $x$. Here $u$ is a mixed square root of 1 and $g_0(y)=2[\min(y,N-y)\le\lfloor N/4\rfloor]-1$ for $0\le y<N$. Its first harmonic has frequency $u$, which is a unit. Control of the factor-divisible frequencies in (1) alone says nothing about this harmonic. [Quarter-window reduction](quarter-fourier.md)

**Bound.** Keep the distinct-shift hypotheses above. Fix any unit $u\bmod N$ and any readout $F$ of the $J$ shifted cubic symbols with $|F|\le1$. The shifts and table may depend on $N$ and $u$, and are held fixed across $x$. A bounded readout on singular inputs can be assigned arbitrarily. Then

$$
\left|\mathbb E_{x\in(\mathbb Z/N)^\times}F(\operatorname{labels}(x))g_0(ux)\right|
\ll \frac{J^2\sqrt{3^J}\log N}{\sqrt N}+
J\left(\frac1p+\frac1q\right).
$$

For binary predictions the success advantage over $1/2$ is at most half this bound. The same holds for randomized predictions by taking their conditional expected sign. For balanced primes and $3^J\le N^{1-\epsilon}$, any fixed $\epsilon>0$, the bound decays faster than every inverse power of $\log N$. Shifts chosen adaptively as a function of $x$ lie outside this estimate.

**Local additive twists.** Let $L_r(x)$ be the vector of local labels at $r=p,q$, and let $W_{e,r}(x)=\omega^{e\cdot L_r(x)}$ on valid inputs and zero elsewhere. For every nonzero $e\in G$ and every frequency $t\bmod r$,

$$
\left|r^{-1}\sum_{x\bmod r}W_{e,r}(x)\exp(2\pi i tx/r)\right|\ll J/\sqrt r .
$$

At $t=0$ this follows from the three cubic class counts above. For $t\ne0$, weight the rational points on each of the same three twists by $\exp(2\pi i tX/r)$. On its smooth projective model $X$ has three simple geometric poles. Indeed, they lie over infinity, where $f_e(\infty)=1$ and the cubic cover is unramified. A function with a simple pole cannot equal $H^r-H$ plus a constant over the algebraic closure. Bombieri's additive-character estimate on a curve of genus $g$ bounds the sum away from the poles by $(2g-2+3+3)\sqrt r$. Removing fibers above any of the at most $2J$ finite excluded parameters costs at most $6J$ points per twist. On valid parameters exactly one twist contributes three points. Taking the class-weighted sum and dividing by three therefore gives the claimed estimate for both split and inert primes. [Bombieri, *On exponential sums in finite fields*, Theorem 5](https://doi.org/10.2307/2373048)

**From the twists to the predictor.** First average over all $x\bmod N$, giving value zero to the predictor on singular inputs. For each $e\ne0$, its global signature word on valid inputs, extended by zero, factors as $W_{e,p}(x_p)W_{e,q}(x_q)$. Every additive Fourier coefficient of this word is $O(J^2/\sqrt N)$ by CRT and the local estimate, including when one or both local frequencies vanish. The geometric-series bound for the quarter interval gives $\sum_h|\widehat g_0(h)|=O(\log N)$. Multiplication of the argument by $u$ permutes these coefficients. The correlation of each nonconstant word with $g_0(ux)$ is consequently $O(J^2\log N/\sqrt N)$.

Expand $F$ on $G$ using normalized coefficients. Parseval and Cauchy–Schwarz give $\sum_{e\ne0}|\widehat F(e)|\le\sqrt{3^J-1}$. Its constant coefficient, of magnitude at most one, multiplies the validity indicator. On the complete residue space the mean of $g_0(ux)$ is $\pm1/N$, since $u$ permutes the residues. Removing the invalid fraction changes this mean by at most $2J/p+2J/q$. Allowing arbitrary bounded predictions on those same inputs costs at most their fraction once more. Finally the nonunit fraction is at most $1/p+1/q$ and $N/\varphi(N)$ is bounded, giving the result. This accounts for the majority baseline and for singularities without multiplying their costs by the number of signatures.

If all singular readings return a single common value, the final $J(1/p+1/q)$ term can be sharpened to $O((J+1)^2/N)$, also on unbalanced pairs satisfying the distinct-shift hypothesis. Include the local unit condition $x\ne0\bmod r$ directly in each word. Deleting that one point changes each local twist by at most $1/r$. The constant word now discards at most $2J+1$ residue classes at each prime. On each such class the points $ux$ form an equally spaced grid in the complete residue interval, so its **unnormalized** quarter-sign sum is $O(1)$. Each intersection of an excluded class modulo $p$ with one modulo $q$ is a single point. Inclusion–exclusion bounds the constant word's normalized correlation by $O((J+1)^2/N)$. The common sentinel multiplies a difference of the unit and valid-unit indicators, both covered by this calculation.

The checker evaluates all local additive twists, verifies their CRT products and computes the best possible table by a majority vote on every signature at each fixture. The 36 unit-scale checks cover two public scales and one mixed root per menu. On $101\cdot103$ with shifts $0,1,3$, allowing a perfect prediction on every singular unit scale gives the optimum mixed-root correlation $0.084510$ (advantage $0.042255$). The bound using the computed Fourier magnitudes is $0.893270$ on that fixture. [Reproduce the finite checks](check_cubic_jacobi_fourier_axes.py)

## Adaptive shifts chosen from earlier cubic answers

Suppose the next shift and each bounded leaf output depend only on $N$, private coins independent of $x$, and previous cubic answers. Fix the coins. The three regular answers give a ternary tree with $D$ terminal leaves (including empty ones), $I=(D-1)/2$ internal nodes, and maximum depth $J<\min(p,q)$. On the first singular answer one may return any bounded prediction, even the correct sign. Remove repeated shifts along a path by reusing the previous answer. Assume every remaining pair of shifts on the same path is distinct at both primes. A pair coinciding at just one prime gives a proper gcd of its difference with $N$, a separate direct factor witness.

For every unit $u$, the absolute correlation with $g_0(ux)$ on uniform units is

$$
\ll \frac{D J^2\log N}{\sqrt N}+
      \frac{(J+1)^2}{N}+ I\left(\frac1p+\frac1q\right).
$$

To see this, fix a clean leaf of depth $j$ prescribing the successive label vector $v$. Its path determines one fixed shift list even though different leaves use different shifts. On valid inputs the leaf indicator is the point mass $1_{L=v}$ on $(\mathbb Z/3)^j$. Its $3^j$ normalized Fourier coefficients each have magnitude $3^{-j}$. Their sum of absolute values is **one**. Imposing validity at every queried shift and $x\ne0$ at both primes, the nonconstant words contribute $O(j^2\log N/\sqrt N)$ per leaf by the twisted-sum bound above. Summing their contributions costs at most $D$ times the maximum. For the constant word the valid-unit mask contributes $O((j+1)^2/N)$ multiplied by $3^{-j}$. The ternary leaf weights sum to one, including in an unbalanced tree. The bounds remain of the same order upon dividing by $\varphi(N)/N$.

Each internal node names a fixed shift. On unit scales the chance that either prime divides its norm denominator is at most $2/(p-1)+2/(q-1)$. A union bound over all $I$ nodes pays for every path that ever hits a singularity, without assuming independence from earlier answers. Averaging over the private coins preserves the estimate whenever the size and depth bounds hold for every realized tree.

On balanced semiprimes this excludes every tree of size polynomial in $\log N$, as well as every tree of depth $J\le(1/2-\epsilon)\log_3N$ for fixed $\epsilon>0$, since $D\le3^J$. The stronger fixed-menu estimate above applies when all leaves read from one list of $J$ shifts. If a query uses the numerical value of $x$, a computed ring element, or another $x$-dependent input to choose a shift, its leaf need no longer prescribe a fixed list and this tree argument does not apply.

## Why a restriction on the menu or the readout is necessary

Now use a single sentinel for singular readings. Fix distinct $x,y\bmod N$ and a uniformly random shift $a$. Fill any singular **local** label by zero temporarily. At a prime where $x\ne y$, the mean cube root of unity of their local label difference has magnitude $O(1/\sqrt r)$ by the two-shift estimate, with the filled exceptions changing the mean by $O(1/r)$. At the other prime the magnitude is at most one. Chinese remaindering multiplies the means. Orthogonality on $\mathbb Z/3$ and the probability $O(1/p+1/q)$ of encountering a singularity yield

$$
 \Pr[\operatorname{label}(x+a)=\operatorname{label}(y+a)]
 \le\tfrac13+O(1/\sqrt{\min(p,q)}).
$$

For $J$ independent shifts, take the $J$th power and sum over fewer than $N^2$ pairs. With $J=\lceil(2+\epsilon)\log_3N\rceil$ this tends to zero on balanced inputs. Consequently an unrestricted table on the labels can represent any prescribed function of $x$, including a single pure factor-frequency phase, when its signature is injective. The argument does not give a fast formula for that table from $N$. Even if one discards every row containing a sentinel, injection persists on the retained subset; when both primes are inert every row is retained. Two finite inert examples already have complete separation using only consecutive shifts. [Reproducible Fourier and injection checks](check_cubic_jacobi_fourier_axes.py)
