# One absolute-trace residue on one disk is equivalent to balanced factoring

Let $N=pq$ with distinct odd primes $p<q\le 2p$, in binary. Use the [public branch](README.md#public-branch-construction) $\ell=1$, $\Delta=4N-c^2>4$ and fixed public $A=\lfloor\sqrt{\Delta/3}\rfloor<p$. A nonunit $c$ directly gives a divisor. The following scalar, at this one cutoff, suffices to factor every remaining input:

$$T_N=\sum_{\substack{w,y>0\\w^2+\Delta y^2\le4NA\\4N\mid w^2+\Delta y^2}}w\pmod N.$$

Equivalently this is one quarter of the full signed disk sum of $|2x+b y|$, with $b=\Delta\bmod 2$, because both axes are empty on the punctured disk. The quadrant and absolute value matter: a signed linear moment over the full disk cancels. This changes both the weight and the output modulus of the original unweighted mod-16 problem.

## Decoder

Write $r=c/\ell\bmod N$ for the public root. The argument also works on the repository's eligible auxiliary branches with invertible $c,\ell$ and the same cap and orbit hypotheses. Set $K=\lfloor\sqrt A\rfloor$, subtract the public ray and normalize:

$$z=r^{-1}\bigl(T_N-cK(K+1)/2\bigr)\pmod N.$$

Let $h$ be either mixed root; its signs relative to $r$ are opposite at $p$ and $q$. It is used only in the proof. Every remaining positive pair has $w\equiv\epsilon h y\pmod N$ for a sign $\epsilon$. Put $D=\sum\epsilon y$. Hence $z\equiv (h/r)D\pmod N$ and $z^2\equiv D^2\pmod N$.

We claim $0<|D|<\sqrt N$. Every individual $y^2\le4NA/\Delta<N$ since $A<\Delta/4$. The [fixed-cap shape lemma](weighted_recurrence.md#shape-at-the-cap) ensures at least one point. On one primitive ray all $m$ points are its consecutive positive multiples, so their relative signs agree and $D\ne0$. For $m\ge2$ the first norm index $a$ satisfies $m^2a\le A$, and $y_0\le2\sqrt{Na/\Delta}$ gives

$$|D|=\frac{m(m+1)}{2}y_0\le (1+1/m)\sqrt{\frac{N}{3a}}<\sqrt N.$$

There remain two independent points. Cancellation would require two opposite-slope points with the same positive $y$, so $w_1+w_2\equiv0\pmod N$. Since $4A<N$ (use $A<p$ and $q\ge5$), both $w_i<N$. Their sum would have to be $N$. But $w_1,w_2$ both have parity $by$, which makes that sum even. Thus $D\ne0$.

Give the two vectors physical coordinates $(w_i,\epsilon_i y_i)$ in the same mixed lattice. They form a basis by the determinant bound. Write $s,t\le A$ for their norm indices and $b_1$ for the integral Gram cross coefficient, with $b_1^2=4st-\Delta\le A^2$. Their sum has index $s+t+b_1$. Since $w_1+w_2>0$,

$$D^2=\frac{4N(s+t+b_1)-(w_1+w_2)^2}{\Delta}<\frac{4N(s+t+|b_1|)}{\Delta}\le\frac{12NA}{\Delta}.$$

For $A\ge4$, $\Delta\ge3A^2\ge12A$. For $A=1$, independence would require $\Delta\le4st\le4$, excluded. For $A=2$, $\Delta\ge12$ forces $s=t=2$. Primitivity and $|b_1|\le2$ give $|b_1|=1$ and $\Delta=15$. But $4N=c^2+15\ell^2$ and $\gcd(N,15)=1$ imply $N\equiv1\pmod3$. The putative shell $w^2+15y^2=8N$ then fails modulo 3. For $A=3$, $\Delta\ge27$ similarly forces $s=t=3$. Primitivity gives $|b_1|=1$ or $2$, so $4(s+t+|b_1|)\le\Delta$ in both cases. This proves the strict bound on $D^2$ throughout.

If $1<\gcd(z,N)<N$ the gcd is already a divisor. The bound $0<|D|<\sqrt N$ and the fact that $h/r$ is a unit exclude $\gcd(z,N)=N$. Otherwise $s=z^2\bmod N$ in $[0,N)$ is the ordinary integer square $D^2$. Take its exact integer square root $t=|D|$. Since $z/t\bmod N$ has opposite signs at $p,q$, $\gcd(z-t,N)$ is a proper divisor. This takes one ordinary square root and a constant number of modular operations on $O(\log N)$-bit inputs.

## Evaluation in both directions

Given the factors, the [reduced-basis control](weighted_recurrence_control.py) constructs one mixed lattice and reduces it in polynomial bit time. Its [cap bound](weighted_recurrence.md#converse-with-supplied-factors) leaves one multiples ray and at most three other candidate vectors to check. At any such vector of norm index $j$ with physical ordinate $y$, $w=\sqrt{4Nj-\Delta y^2}$ is an exact nonnegative integer. The ray contributes $w_0 m(m+1)/2$ in closed form. Add the known public ray and reduce modulo $N$. Hence computing $T_N$ on the branch constructed from $N$ and factoring the promised family are polynomial-time equivalent. The input-only evaluator still uses the public residue join. Computing $T_N$ faster is the unresolved step.

All 764 exhaustive small inputs, public fixtures and additional tie controls agreed between the input-only join, an independent exact-shell scan and the factor-supplied control. The scalar decoder returned a factor on the first square test in each. On the 36-bit input $N=38942974393$, the input-only join computed $T_N=69529082$ with 7,467,686 probes and 12,042 insertions. Public subtraction and normalization gave $z=22336814480$. Its square modulo $N$ was $72233001=8499^2$, and $\gcd(z-8499,N)=254119$ gave the complementary prime $153247$ by division. On a separate 241-bit input the factor-supplied evaluator took 33 Gauss steps and its scalar decoded a supplied prime. The large calibration uses the factors in the evaluator; the 36-bit run starts with $N$ alone. [Checks](weighted_recurrence_checks.py) and [recorded output](weighted_recurrence_results.json)

Run the input-only scalar factorer from the repository root:

```sh
uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_aggregate/weighted_recurrence_single.py 38942974393
```
