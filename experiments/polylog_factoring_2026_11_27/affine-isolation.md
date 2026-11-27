# Nonadaptive weighted disk reduction by mixed root isolation

For any odd input with at least two distinct prime divisors, a randomized batch of $O(n^2)$ weighted disk queries gives a split with constant probability, where $n=N.\mathrm{bit\_length}()$. All queries use one public numerical cutoff. The weight is the fixed polynomial $xy$, and the output is taken modulo $N$. Amplification gives failure probability at most $\varepsilon$ using $O(n^2\log(1/\varepsilon))$ queries. This extends the nonadaptive approach beyond two prime divisors by using a stronger weighted oracle. [Two-prime-divisor parity reduction](nonadaptive-quarter.md)

## Isolation in an affine window

Let $S$ be the set of mixed roots of 1 modulo $N$, and let $R=|S|$. Independently choose a uniform unit $a$ and a uniform offset $b\bmod N$. A root $r$ maps to the representative $z=ar+b\bmod N$ in $[0,N)$. Keep it if $0\le z<w$. A uniform multiplier can be sampled by gcd, with any proper nonunit already supplying a split.

For distinct roots $r,r'$, put $d=\gcd(r-r',N)$ and $m=N/d$. Their image difference is uniform on $d(\mathbb Z/m\mathbb Z)^\times$. If $w\le d$, two images cannot both fall in the interval. Otherwise at most $2w/d$ possible differences can fit. Since $m/\varphi(m)\le\omega(m)+1\le n+1$ by sorting its prime factors and telescoping, the conditional probability of a difference that fits is at most $2(n+1)w/N$. Independently, for each fixed multiplier the offset puts the first image in the interval with probability $w/N$. Thus for the number $X$ of retained mixed roots,
$$\mathbb E X=\lambda=Rw/N,\qquad
\mathbb E[X(X-1)]\le 2(n+1)\lambda^2.$$

If $N/(4(n+1)R)\ge1$, choose the largest dyadic $w$ below that ratio, allowing equality. Then $1/(8(n+1))<\lambda\le1/(4(n+1))$ and $\Pr(X=1)\ge\mathbb E X-\mathbb E[X(X-1)]>1/(16(n+1))$. Otherwise $w=1$ gives singleton probability $R/N>1/(4(n+1))$ directly. The suitable width always appears among dyadic widths at most $N/4$. Since $R$ is unknown, prepare independent trials at every such width. A singleton is sufficient, and every returned candidate is verified by gcd.

## One weighted scalar reads the singleton

For a width $w$, set $K=\lfloor w/2\rfloor$ and $s=K-b$. Choose a primitive quadratic form $Q(x,y)=Ax^2+Bxy+Cy^2$ which reduces to $(x+sy)^2-a^2y^2$ modulo $N$. Set $A=(N+1)^4$, $J=\lfloor N/8\rfloor$, $P=N(N+1)$ and $C_*=A(J+1)^2$. At every width use the same budget
$$H=N\left\lceil(2C_*+PJ+N)/N\right\rceil.$$
Take $B\in[0,P)$ with $B\equiv2s\pmod N$, $B\equiv1\pmod{N+1}$. Choose $C$ in $[C_K,C_K+N)$ with $C\equiv s^2-a^2\pmod N$, where $C_K=H-AK^2-PK-N$. Then $C_K\ge C_*$. At $y=1$ all $|x|\le K$ fit. At $|x|\ge K+1$, the excess over $H$ is bounded below by $(A-P)(2K+1)-N>0$. Completing the square excludes every $|y|\ge2$ because $H<3C_*$ and $4C_*-P^2/A>3C_*$. The form is primitive since $\gcd(A,B)=1$. Its cutoff $H/N$ is below $\sqrt{(4AC-B^2)/3}$, using $H<3C_*$ and $3(J+1)<N$.

On the $y=1$ row the congruence is exactly $x+K=ar+b\bmod N$ for a root $r$. The opposite point $(-x,-y)$ has the same weight $xy$. Hence half the signed weighted residue modulo $N$ is the sum of $x=z-K$ over all roots with $0\le z\le2K$. The axis has weight zero, even for repeated factors. For even $w$, test and remove the extra endpoint at $z=w$. Remove the predictable contributions of $r=\pm1$. If exactly one mixed root remains, adding $K$ recovers $z$ modulo $N$, and then $(z-b)/a$ recovers the root. All queries, coefficients and weights have polynomial bit descriptions. A polynomial-bit-time evaluator for the weighted disk remains open.

## Checks

The independent signed-disk scan agreed with the affine window moment on 12,800 queries across ten inputs, including repeated factors and inputs with four distinct prime divisors. It decoded every one of the 2,122 singleton windows. Run `python affine_isolation.py`.
