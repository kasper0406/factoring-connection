# When does the public product justify a hard YES?

Take two distinct uniformly sampled primes in a common dyadic band $[R,2R)$, an odd fixed multiplier $h$, and a dyadic period $D\mid R$. Write $N=pq$, $z=N/R^2$, $\lambda=hR/D$, $H(t)=1_{\{t\}\ge1/2}$, and $Y=H(hp/D)H(hq/D)$. The constant NO classifier is the best constant at large $\lambda$. For the continuous prime-density model conditioned on the product, its optimal accuracy improvement is

$$
 \left(\frac1{24}+o(1)\right)\frac{\log\lambda}{\lambda^2}.
 \tag{1}
$$

Here and below $\log$ is natural and $\lambda,R\to\infty$. The same leading term is attained on actual prime pairs by a public binary rule when $\lambda\le(\log R)^A$ for fixed $A$. It is optimal among product-and-residue classifiers of polylogarithmically bounded variation, with a variation allowance large enough to include this rule. The earlier conditional forecast gains $\Theta((\log\lambda)/\lambda)$ in squared loss. Turning a modest change in probability into a justified YES is much rarer. [Sampling model, conditional forecast and prime transfer](hyperbola-bit-forecast.md)

**A rule with integer comparisons.** On $z\ge2$ set $K=\lambda(4-z)/2$. Predict YES if $K\le1/2$, or if there is an integer $j\ge1$ with

$$
 j^3\le\lambda,\qquad |K-j-\tfrac12|\le\frac1{4j}.
 \tag{2}
$$

Predict NO otherwise. The only possible $j$ is $\lfloor K\rfloor$. The interval $K\le1/2$ is a sure YES in the upper corner, up to irrelevant boundary points. The additional windows account for the logarithmic factor. All quantities in (2) can be compared exactly using $N,h,D$ and the common dyadic band, which is inferable from $N$ on this model. [Public decision and finite checks](check_hyperbola_hard_decisions.py)

## Hard-positive contours concentrate in small windows

Let $w_R(x)$ be the normalized density proportional to $1/\log(Rx)$ on $[1,2]$. For fixed $z$, put $a=\max(1,z/2)$, $b=\min(2,z)$, $L=b-a$, so that $ab=z$. On $x\in[a,b]$ set $W_z(x)=w_R(x)w_R(z/x)/x$. This contour weight is bounded above and below by positive constants and has uniformly bounded first derivative. Write $d(z)=\int W_z(x)\,dx$ and $u(z)=\int W_z(x)H(\lambda x)H(\lambda z/x)\,dx$ on that interval. The unconditional gain of the optimal hard rule over NO is exactly $\int_1^4(2u-d)_+\,dz$.

For $1\le z\le2$ this integrand vanishes. Indeed $a=1$, $\lambda$ is an integer, and every prefix integral of $H(t)-1/2$ starting there in the scaled variable is nonpositive. The weight decreases in $x$ since its logarithmic derivative is $[-1-1/\log(Rx)+1/\log(Rz/x)]/x<0$ for large $R$. Integration by parts gives a one-bit marginal of at most $d/2$, which also bounds the joint hit $u$.

For the upper corner put $K=\lambda L$ as in (2). At any contour with $2u>d$, the uniform estimate $|u-d/4|\ll\lambda^{-1/2}$ from the conditional forecast and $d\asymp L$ imply $K=O(\sqrt\lambda)$. Each one-bit marginal differs from $d/2$ by $O(1/\lambda)$, by integrating its mean-zero periodic function against the smooth weight. The two marginal integrals are equal by $x\mapsto z/x$. The identity

$$
 2BC-1=(B-\tfrac12)+(C-\tfrac12)-|B-C|,\quad B,C\in\{0,1\},
$$

therefore gives $2u-d\ll1/\lambda$ and, on a hard-positive contour,

$$
 \int_{\lambda a}^{\lambda b}
   |H(t)-H(\phi(t)-t)|\,dt=O(1),\qquad
 \phi(t)=t+\lambda^2z/t.
$$

To see where this can hold, freeze $\phi$ on each full period of $t$. For constant $v$ the length of the mismatch set in a full period is $2\|v-1/2\|$, where $\|\cdot\|$ is distance to the nearest integer. Changing the phase by at most $\epsilon$ can affect only the points within $\epsilon$ of its two jumps. Here $\sup|\phi'|=O(K/\lambda)$, so summing the freezing errors costs $O(1+K^2/\lambda)=O(1)$ including the incomplete ends. The distance function itself is Lipschitz. This yields $\int\|\phi(t)-1/2\|\,dt=O(1)$.

On the rightmost quarter of the contour, $\phi'>0$ is comparable to $K/\lambda$. Indeed the stationary point $x=\sqrt{ab}$ is at or below $(a+b)/2$. The phase covers an interval of length $\Delta\asymp K^2/\lambda$ there. The distance to a half-integer integrated over any interval of length $\Delta$ is at least a constant times $\min(\Delta,\Delta^2)$. A change of variable then gives a lower bound $\gg\min(K,K^3/\lambda)$ for the preceding distance integral. Hence $K=O(\lambda^{1/3})$.

At the endpoint the phase is $\phi_e=\lambda(a+b)=4\lambda-K$, and throughout the contour

$$
 0\le\phi_e-\phi(\lambda x)
   =\lambda\frac{(b-x)(x-a)}{x}\ll K^2/\lambda.
$$

For $K\ge1$, averaging the distance integral gives $\|\phi_e-1/2\|\ll1/K+K^2/\lambda\ll1/K$. Thus positivity requires $K$ within $O(1/j)$ of $j+1/2$ for $1\le j=O(\lambda^{1/3})$, apart from a bounded initial range of $K$. The width in $z$ is $O(1/(\lambda j))$ and the signed height is at most $O(1/\lambda)$. This gives $O((\log\lambda)/\lambda^2)$ even before computing the constant.

## The leading constant

It suffices to evaluate the bands $\log\lambda\le j\le\lambda^{1/3}/\log\lambda$. The omitted positive part has mass $O((1+\log\log\lambda)/\lambda^2)$ by the window bound. In these bands write $K=j+1/2+\delta$ with $|\delta|\le C/j$ for a fixed sufficiently large $C$. The unweighted integral of $2(H(t)-1/2)$ across the contour is $1/2-|\delta|$, since its upper endpoint $2\lambda$ is an integer. With $\phi=\phi_e$, each of the $j$ complete periods gives mismatch $2|\delta|$ and the last partial period has mismatch $O(|\delta|)$. Replacing $\phi_e$ by the actual phase changes mismatch on at most $O((K+1)K^2/\lambda)=o(1)$ total length. Thus the signed unweighted height in the $t$ coordinate is

$$
 1/2-2j|\delta|+o(1).
$$

For this calculation one can use twice the first single marginal minus the mismatch. The weighted marginal equality is exact. Replacing the weight by $W_z(2)$ introduces at most $O(K^2/\lambda)=o(1)$ error in scaled units, by its bounded derivative. Also $W_z(2)=1/2+o(1)$ uniformly here. These facts yield

$$
 (2u-d)_+=\lambda^{-1}
     \left(\tfrac12(1/2-2j|\delta|)_++o(1)\right)
$$

uniformly in the indicated windows. The area of the positive triangle in $\delta$ is $1/(8j)$ before the factor $1/2$. Since $|dz/dK|=2/\lambda$, the total contribution is $\lambda^{-2}\sum_j(1/(8j)+o(1/j))$. The harmonic sum over these bands is $(1/3+o(1))\log\lambda$, proving (1).

The rule (2) takes exactly the triangular support on the middle bands. In its additional windows $j^3\le\lambda$, the magnitude of its signed height in scaled units stays $O(1)$. This follows from the same frozen-phase comparison, now with error $O(1)$ in place of $o(1)$. Their signed contribution is $O((1+\log\log\lambda)/\lambda^2)$. This proves (1) for the explicit rule as well.

## What transfers to primes

For every allowed binary classifier $G_c(z)$ indexed by $c=N\bmod m$ for a fixed $m\le(\log R)^A$, the preceding prime transfer applies to the gain $\mathbb E[(2Y-1)G_c]$. No such classifier with polylogarithmic variation exceeds the continuous optimum except by $o((\log R)^{-B})$ for every fixed $B$ in the setting of the conditional forecast. One can use the Lipschitz positive part of $2\mu-1$ for this upper bound, where $\mu=u/d$. The explicit rule has variation $O(\lambda^{1/3})$ and therefore transfers with a sufficient fixed power of $\log R$. If an independently imposed variation budget $(\log R)^C$ is smaller, restricting $j$ in (2) to at most a sufficiently small constant times $(\log R)^C$ still attains order $\Theta_{A,C}((\log\lambda)/\lambda^2)$ for every fixed $A,C>0$. The full leading constant requires enough accepted windows. These statements concern the indicated product-scale score class.

The independent checker compared 105 continuum contour heights against the triangular approximation with explicit finite tolerances, and matched six small prime-pair censuses by direct enumeration. Larger censuses used binary search to count all ordered pairs in the predicted windows. At $R=2^{21}$, $h=1$, $\lambda=128$, the rule selected 3,753,088 pairs, including 1,979,116 true positives. The sure upper corner contained 142,146 pairs. At $h=3$, $\lambda=384$, the corresponding counts were 559,098, 287,418, and 14,046. These finite gains can differ considerably from the leading asymptotic constant. [Reproduce the comparisons](check_hyperbola_hard_decisions.py)
