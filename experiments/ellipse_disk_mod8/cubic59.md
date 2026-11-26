# The discriminant −59 bit on an unfactored semiprime

On the canonical branch $4N=c^2+59$, take $N=pq$ with distinct balanced primes coprime to $6\cdot59$, and $c=\lfloor 2\sqrt N\rfloor$. The cubic $f(Y)=Y^3+2Y+1$ either splits completely at both primes (bit 0) or is irreducible at both (bit 1). The previous experiment derives this promise from inverse ideal classes in the Hilbert class field of discriminant $-59$. [Class-character reduction](class_trace.md)

**Result.** The public branch gives a Cardano numerator that contains the exact bit at either prime. Its cubic Jacobi symbol at $N$ is identically 1. A second formulation gives an explicit quartic with exactly one root at bit 1 and either zero or sixteen at bit 0. We tested two ways to extract the bit. A pair of shifted norm-Jacobi symbols has provably vanishing distinguishing advantage for uniform samples. Substituting $N$ for a prime in the quartic root calculation missed all 221 unique roots in the benchmark. There is still no polylogarithmic decision algorithm here.

## What the repository already knew

| Prior work | Consequence for this test |
|---|---|
| [Coefficient experiment at −23](../../notes/coef_results.md) and [kernel experiment](../../notes/kern_results.md) | A fixed class-number-three form can encode the bit. Its known representation count scans to order $\sqrt N$. In a cubic algebra, certifying a norm-one nonsquare would expose the split case. |
| [Hidden class of a prime, target T2](../../notes/hidclass_results.md) | For a cyclic cubic with both local types equal, a uniform norm-one element is a square with probability 1 in the inert case and $1/16$ in the split case. Its earlier signed norm probe used a different numerator and has a scalar-rescaling symmetry. |
| [Rational-ideal reciprocity, Theorem RA](../../notes/newfamily_results.md) and [prime-symbol implementation](../../src/newfamily/diag.py) | For a **fixed** numerator the global residue symbol factors through a fixed conductor. Our numerator varies with $c$, so it requires a separate calculation. We reuse the local symbol code only for the controls. |
| [Instance-local oracle bound, Prop IL](../../notes/tracedata_results.md) | One bit at one $N$ can be guessed both ways and candidate factors verified. A fast factoring reduction using that bit needs more information or different queries. Its fixed-field Chebotarev argument on unrestricted primes does not prove a density statement on this near-square family. |

## Exact cancellation in cubic reciprocity

Let $\omega^2+\omega+1=0$, $\delta=1+2\omega=\sqrt{-3}$, and work in the Eisenstein integers. Since $3\nmid N$, the branch identity forces $3\mid c$ and $N\equiv-1\bmod9$. The primes have opposite residues modulo 3. Put

$$
\beta=(c+3\delta)/2,\qquad \beta\bar\beta=(c^2+27)/4=N-8 .
$$

At $r=p,q$, the usual Cardano radicand for $f$ is $b=(-9+c\delta)/18$, with $\delta^3b=\beta$. Its two conjugates have product $(-2/3)^3$. Since $c^2=-59\bmod r$, the split-discriminant cubic criterion says that $f$ splits exactly when $b$ is a cube in the field with $\omega$ adjoined. At a split Eisenstein prime the two embeddings have the same cube status. Thus $\beta$ is a cube in the full algebra at $r$ exactly for bit 0. [Dickson's criterion, Theorem 1](https://eprint.iacr.org/2009/457.pdf)

Nevertheless $(\beta/N)_3=1$. Remove the factors of $1-\omega$ and multiply by a unit to make a primary associate $\gamma$ of the remaining part of $\beta$. Every removed factor and unit has cubic symbol 1 at rational $N\equiv-1\bmod9$ by the supplementary laws. Cubic reciprocity now exchanges $\gamma$ and $N$. Since $\gamma\mid\beta\bar\beta=N-8$, the reverse numerator is $8=2^3\bmod\gamma$. This proves cancellation for every promised input, including either local bit. In the benchmark, the two rational-prime symbol exponents were always $(0,0)$, $(1,2)$ or $(2,1)$. [Cubic reciprocity and supplements, equations (3)–(6)](https://arxiv.org/pdf/1807.07719)

## A cyclic cubic with a public automorphism

Set $t=(c-3)/2$, so $t^2+3t+17=N$. The Shanks polynomial

$$S_t(Z)=Z^3-tZ^2-(t+3)Z-1$$

has discriminant $(t^2+3t+9)^2=(N-8)^2$, a unit modulo $N$. Its root transformation $\sigma(Z)=-(Z+1)/Z$ has order three. To check that its local splitting is the desired bit, put $U=\beta/\bar\beta$. A cube root of $\beta$ in the quadratic algebra gives a norm-one cube root $v$ of $U$ by dividing the root by its conjugate. Conversely $v^3=U$ gives $(-\beta/(2v))^3=\beta$ since $\beta\bar\beta=-8\bmod N$. The substitution $v=(Z-\omega^2)/(Z-\omega)$ identifies $v^3=U$ with $S_t(Z)=0$, by clearing denominators. This is the usual rational parametrization of the norm-one conic. The omitted point at infinity gives $v=1$, while $U\ne1$ at either prime. The possible poles are disjoint from the roots since $(t-3\omega)(t-3\omega^2)=N-8$ is a unit. Hence $S_t$ has three or zero roots at each prime with the same bit as $f$. [Shanks cubic family](https://www.math.toronto.edu/graduate/PhDthesis_PeterCho.pdf)

Write $A=(\mathbb Z/N)[Z]/S_t$. A uniform unit $a$ gives a uniform norm-one element $u=\sigma(a)/a$ by the finite-field or split-algebra form of Hilbert 90 at each prime. Locally the norm-one group has odd order $r^2+r+1$ in the inert case. In the split case it is $(\mathbb F_r^\times)^2$, so exactly one quarter of its elements are squares. Squareness inside the torus agrees with squareness in $A^\times$ for norm-one elements: a square root of norm $-1$ locally can be replaced by its negative. This transfers the repository's $1$ versus $1/16$ square-probability test to the public branch algebra.

## A nonlinear public probe and its limit

We tried the two ordinary Jacobi symbols $J_m=(\operatorname{Nm}(1+u^m)/N)$ for $m=1,2$ **on the same draw**. The shift removes the immediate constant $\operatorname{Nm}(u)=1$. Even their joint distribution approaches two independent fair signs for either bit. For balanced primes the total variation error is $O(N^{-1/2})$ per draw, including zero signs. The following proof concerns only the reported signs or gcd events from uniform independent draws. It makes no assertion about algorithms that also use the sampled coefficients or choose correlated inputs.

Over $\mathbb F_r$, write the three geometric coordinates of $a=b+s$ as $b_i+s$, with $\sigma$ cycling them and $s\in\mathbb F_r$. Averaging over uniform $b$ and an independent uniform $s$ keeps $a$ uniform. The quadratic character for $m=1$ is the character of the product of the three denominator and three numerator linear factors. Its six roots in $s$ are $-b_i$ and $-(b_i+b_{i+1})/2$. For $m=2$ the denominators are squares. The roots of the three numerator quadratics are $(\xi b_{i+1}-b_i)/(1-\xi)$ for $\xi^2=-1$. The two sets give twelve distinct linear forms in the three independent geometric coordinates in characteristic greater than three. Forms supported on the same pair have different coefficients.

A collision of two forms imposes a nonzero linear equation on $b$ over the algebraic closure. The rational points of the algebra span its scalar extension, so each equation holds for at most $r^2$ of the $r^3$ values of $b$, also in the inert case. Off at most $\binom{12}{2}$ such hyperplanes, each of the two polynomials and their product is squarefree of degree at most 12. The quadratic character sum over $s$ then has magnitude at most $11\sqrt r$. Removing nonunits and zeros changes a local average by $O(1/r)$, by the degree bound on each line. Thus each local sign and their product have mean $O(r^{-1/2})$. Independence across $p,q$ makes the three global means $O(N^{-1/2})$. Those three moments determine the sign-pair law up to the zero probability $O(1/p+1/q)$. Polynomially many independent repetitions in $\log N$ still give vanishing advantage from this statistic. [Multiplicative Weil bound, Theorem 5](https://page.mi.fu-berlin.de/szabo/PDF/SubgraphsPNG.pdf)

## A deterministic quartic for the next test

The Shanks generator $z$ has norm 1. For a norm-one square root $v$ of $z$, put $X=\operatorname{Tr}v$, $Y=\operatorname{Tr}(v^{-1})$. The elementary symmetric functions give $X^2-2Y=t$ and $Y^2-2X=-(t+3)$. Eliminating $Y$ gives

$$
H_t(X)=(X^2-t)^2-8X+4(t+3)
\equiv X^4-2tX^2-8X+t-5\pmod N .
$$

Over an algebraic closure there are four norm-one square roots, obtained by the four even sign changes. Their traces are distinct because equality would force two components of $z$ to coincide. The trace map therefore identifies these four roots equivariantly with the roots of $H_t$. At an inert prime the norm-one group has odd order and the quartic has exactly one root. At a split prime it has four roots if $z$ is a square and zero otherwise. By CRT the quartic consequently has **one** root modulo $N$ at bit 1 and **zero or sixteen** at bit 0. Mere existence can coincide across the two cases.

At an inert prime $r$, $v=z^{(r^2+r+2)/2}$ gives the unique root by taking its trace. We verified this at both prime factors of every bit-1 control. The public substitution $r=N$ gave no quartic root and no gcd factor on any of the 320 inputs.

There is also a limit on closed-form root guesses. Let $\tau$ satisfy $\tau^2+3\tau+17=0$ in $K=\mathbb Q(\sqrt{-59})$. The monic quartic $H_\tau$ has no root in $K$: modulo the prime $(17,\tau)$ its reduction is $X^4-8X-5$, which has no root in $\mathbb F_{17}$ by testing the 17 residues. A root in $K$ would be integral there. For a fixed rational expression $R(\tau)$ defined in $K$, clear denominators in $H_\tau(R(\tau))$. The result is a fixed nonzero algebraic integer. Wherever the public expression is defined, $H_t(R(t))=0\bmod N$ would force $N$ to divide its fixed nonzero integer norm. Thus any **finite fixed list** of rational formulas in $t$ can produce roots for only finitely many such moduli. This leaves variable-length recurrences using $N$ open.

## Verification and decision

| Modulus bits | Inputs | Bit 0 / bit 1 | Nonsquare torus draws at bit 0 | Squares at bit 1 | Quartic root counts $0/1/16$ |
|---|---:|---:|---:|---:|---:|
| 15–31 | 128 | 34 / 94 | 512 of 543 | 1,504 of 1,504 | 32 / 94 / 2 |
| 37–38 | 64 | 22 / 42 | 337 of 352 | 672 of 672 | 21 / 42 / 1 |
| 61 | 64 | 20 / 44 | 304 of 320 | 704 of 704 | 18 / 44 / 2 |
| 81 | 64 | 23 / 41 | 339 of 368 | 656 of 656 | 23 / 41 / 0 |

Each tier takes consecutive eligible examples from a deterministic prime search with a different starting point. There were 16 reproducible draws per modulus, with one nonunit draw omitted from the square statistics. Two small inputs yielded incidental gcd factors during sampling. All larger tiers yielded zero. Local cubic symbols, polynomial roots, and square tests use fixture factors only on the control side. The public arithmetic uses $N$ alone and fixed-degree rings. Independently, exhaustive enumeration at 5 and 17 gave respectively 31 squares among 31 norm-one units and 64 among 256.

Three small regression inputs cover all outcomes. The quartic has zero roots for $N=252{,}521$ ($t=501$), one for $N=31{,}877$ ($t=177$), and sixteen for $N=29{,}014{,}397$ ($t=5{,}385$).

**Decision.** Keep $H_t$ as a small benchmark for a genuinely different public root-parity or unique-root method. Charge all setup and require it to handle the sixteen-root controls. The immediate hurdle is a variable-length computation exploiting $t^2+3t+17=N$ which distinguishes the unique-root case. The two shifted Jacobi symbols, fixed rational guesses and prime-exponent substitution above have specific obstructions. Even solving the fixed family would leave growing discriminants and factoring coverage open: $4N=c^2+59$ describes only $O(\sqrt X)$ integers $N\le X$.

Run from the repository root with python-flint, SymPy, NumPy and gmpy2 installed:

    python experiments/ellipse_disk_mod8/cubic59.py 252521 --draws 16
    python experiments/ellipse_disk_mod8/cubic59_check.py --output rerun.json

[Public probes](cubic59.py) · [Control generator and independent checks](cubic59_check.py) · [Results including all sign-pair histograms](cubic59_results.json)
