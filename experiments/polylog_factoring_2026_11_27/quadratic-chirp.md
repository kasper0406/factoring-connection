# Complete quadratic weighting of the fixed phases

Write $M=N^2$, $e_m(x)=\exp(2\pi i x/m)$, and $C_P(t)=N^{-1}\sum_{x\bmod M}e_M(tP(x))$ for the literal cubic or quartic. At unit scales the CRT stationary formula gives an extension $W(t)=\sum_f w_f e_M(ft)$ defined on all scales. Its values at nonunits can differ from $C_P(t)$. The formulas below specify the domain and normalization of each average.

## A native evaluator for projected literal sums

Take the outer quadratic coefficient $Na$ with $a$ a unit modulo $N$. For the literal raw sum define
$$J=\frac1{N\sqrt N}\sum_{t\bmod M} C_P(t)e_M(Nat^2+Bt).$$
Writing $t=u+Nv$ projects onto $N\mid P(x)+B$. The surviving kernel is the complete quadratic Gauss sum modulo $N$ with linear coefficient $(P(x)+B)/N$. Put $d=(4a)^{-1}\bmod N$ and $\eta=\epsilon_N(a\mid N)$, where $\epsilon_N$ is one for $N\equiv1\pmod4$ and $i$ otherwise. On writing $x=x_0+Nw$ the exact identity is
$$J=\eta\sum_{P(x_0)+B\equiv0\ (N)}\frac1N\sum_{w\bmod N}e_N\big(-d(\beta_0+wP'(x_0))^2\big),\qquad \beta_0=\frac{P(x_0)+B}{N}.$$
This expression is used for the proof and the factor-supplied control.

For squarefree $N=pq$, a term with derivative nonzero at either prime has magnitude at most $1/\sqrt{\min(p,q)}$. There are at most $\deg(P)^2$ terms. Terms with zero derivative at both primes are exactly the stationary extension under the same quadratic weight. Test the gcds of $B+f$ with $N$ for each public local label $f$. A proper gcd yields a divisor. Otherwise the projection retains a sole public diagonal when $B\equiv-f\pmod N$, or gives zero. For the diagonal compute its weight and the one standard Gauss phase via the Jacobi symbol. Trial division through $(2\deg(P)^2/\varepsilon)^2$ either gives a factor or bounds the discarded terms by $\varepsilon/2$. Including the numerical phase, this gives polynomial bit complexity in $\log N$ and $1/\varepsilon$ with no supplied factors. [Classical Gauss evaluations](https://mathweb.ucsd.edu/~revans/contents.html) · [Stationary phase](https://doi.org/10.4064/aa-80-1-1-48)

If all $B+f$ are units, each quotient root is simple at every component. Then the inner change of variable is bijective and the two Gauss factors cancel their phases. The exact value is $J=N^{-1/2}\#\{x_0:P(x_0)+B=0\pmod N\}$, also for repeated prime powers. This identifies the additional root-count information at precision $O(N^{-1/2})$. For instance, if the offset is constructed with a public factorization of the shifted polynomial into distinct linear factors modulo $N$, and all differences are units, the count is exactly $\deg(P)^2$ on two-component inputs. A general critical offset with repeated factors needs a separate calculation, since its discarded terms need not be small.

## What changes for unit outer coefficients

For $A$ a unit and $D=(4A)^{-1}\bmod M$, the stationary extension obeys
$$\frac1N\sum_{t\bmod M} W(t)e_M(At^2+Bt)=e_M(-DB^2)\sum_f w_f e_M(-Df^2-2DBf).$$
If the sum uses only $t$ coprime to $N$, insert the exact mask $\gcd(B+f,N)=1$. The block $t=u+Nv$ survives at the unique stationary $u\equiv-(B+f)/(2A)\pmod N$, which proves the mask. Test the gcds for each local label. A proper gcd splits the input, and otherwise a global collision deletes that label at every component. In particular, the unit-scale raw cubic at $B\equiv\pm2$ and raw quartic at $B\equiv1$ leave a single public line.

For a spectrum satisfying $f^2+u f+v=0\pmod M$, the *all-scale extension* transform is a point evaluation at $(u-2B)D$ times the public phase $e_M((v-B^2)D)$. Raw cubic and quartic have $(u,v)=(0,-4)$ and $(1,0)$ respectively. A generic unit-coefficient weighted evaluator of the extension therefore retains the original magnitude problem.

For the unit-scale squared magnitude, local labels are $0,\pm h$, with $h=4$ for cubic and $h=1$ for quartic. Set $B=-t_0/(2D)$ for a desired unit scale $t_0$. When the mask retains all labels, the corrected transform as a function of $D$ is a sequence of at most $R=2^k$ unit-circle modes for $k$ prime-power components. Indeed $f^2$ has labels $0$ or $h^2$ componentwise. At $D=0$ the sequence gives the squared magnitude. For fixed $k$, choose $2R$ samples $D=jD_0$ at one of $D_0=1,2R+1,2R+2$. The three sets are disjoint. For large $N$, at most two values of $D$ have a global label collision, a proper gcd yields a split, and small $N=O(R^2)$ can be trial-divided in time polynomial in $R$.

Extrapolation is stable without separating the modes. A unit-mode sequence $y_j$ of rank at most $R$ has a backward characteristic recurrence $y_j+\sum_{i=1}^R b_i y_{j+i}=0$ with coefficient sum at most $2^R$, padding the roots if needed. With samples $\tilde y_1,\ldots,\tilde y_{2R}$ at absolute error $\delta$, fit $b$ on the $R$ observed rows using least squares with penalty $\delta^2\|b\|_2^2$. The true recurrence gives an objective at most $(R+1)4^R\delta^2$. For the fitted $b$, the true residual sequence is itself of rank at most $R$. Its first $R$ residuals have size $O(R2^R\delta)$, so its unobserved residual at zero has size $O(R4^R\delta)$ by the true recurrence. The estimate $-\sum_i b_i\tilde y_i$ has the same error bound. Exact Gaussian-rational linear algebra charges finite precision and avoids unstable root fitting. Eight samples suffice for two components. A fixed inverse-polynomial point accuracy requires inverse-polynomial chirp accuracy in that case.

## Literal all-scale quartic square

The actual $|C_P(t)|^2$ includes nonunit contributions. For the quartic and a unit $A$, at offset $B=0$ its exact transform is
$$T(A)=\frac1N\sum_{t\bmod M}|C_P(t)|^2 e_M(At^2)=\prod_{q=p^e\mid\mid N}\big(C_q+4e_M(-D e_q)\big),\qquad C_{p^e}=p^{e-1}\big((5e-2)p-5e-\chi_p(-1)\big),$$
where $e_q$ is the CRT idempotent for the $q^2$ component. Here $C_q$ counts pairs with equal quartic value modulo $q$.

Locally at $p^{2e}$, the sector $p^e\mid t$ gives $C_q$ by summing over the lower-modulus scales and using additive orthogonality. At valuation $h<e$ the complete quartic reduces to modulus $p^{2e-h}$. Stationary phase gives squared magnitude $p^h$ times a constant five plus two opposite harmonics with weights $2$ for even $h$ and $2\chi_p(-2)$ for odd $h$. The constants average to zero with the restricted unit quadratic weight. For $h>0$ each harmonic has unit linear coefficient and a quadratic coefficient divisible by $p$, so cancels on the upper digit. At $h=0$ the two harmonics give the equal Gauss value $4e_M(-D e_q)$. This proves the transform over repeated factors as well.

For the count, the equality factors as $(x-y)(x+y)(x^2+y^2-2)=0$. Modulo $p$ the two lines and conic give $3p-5-\chi_p(-1)$ points, with five transverse crossings and no triple crossing. Each smooth point lifts to $p^{e-1}$ points modulo $p^e$. At a crossing the two functions are local coordinates. Count $p^e\mid uv$ with $p\mid u,v$. The zero $u$ contributes $p^{e-1}$ choices of $v$. Each valuation $1\le v_p(u)<e$ contributes $(p-1)p^{e-1}$ pairs. Combining the five crossing counts with the smooth points yields the displayed $C_q$.

With two components $P,Q$, five values at $D=jt_0$, $1\le j\le5$, recover the original quartic point magnitude. Put $x=e_M(-t_0 e_P)$, $y=e_M(-t_0 e_Q)$, and $z=xy=e_M(-t_0)$, which is public. Subtract $16z^j$ and take consecutive differences, scaled by $\sqrt N$. The four differences $h_j$ satisfy $h_{j+2}-(x+y)h_{j+1}+zh_j=0$. From samples at absolute error $\delta\sqrt N$, estimate $x+y$ by the ratio using the larger measured denominator $h_2,h_3$ when it exceeds $\delta^{2/3}$. The slope error is $O(\delta^{1/3})$. If both are below threshold, the recurrence bounds $|h_0|$ by $5|h_2|+2|h_3|$. Since $-\Re h_0=4[C_Q(1-\Re x)+C_P(1-\Re y)]/\sqrt N$ and $\max(C_P,C_Q)\ge\sqrt N$, one of $x,y$ is within $O(\delta^{1/3})$ of one. The other is then close to $z$. In both branches the magnitude estimate has error $O(\delta^{1/3})$, including imbalanced and repeated components. The large unknown constant $C_P C_Q$ cancels in the differences.

At the same zero offset the *unit-only* squared quartic transform is public given the component count: the mask retains only labels $\pm1$ at every component, so the value is $4^k e_M(-D)$. These observations leave the native unit-coefficient literal evaluator and the original point magnitude open. They separate projections where the useful phase has been removed from weighted averages that retain it.

The diagnostic matched 24 literal projected chirps, nine literal squared chirps, 12 unfactored coarse projections, six noisy eight-sample extrapolations, and nine noisy five-chirp point evaluations. [Evaluator and checks](quadratic_chirp.py)

```sh
python experiments/polylog_factoring_2026_11_27/quadratic_chirp.py
```
