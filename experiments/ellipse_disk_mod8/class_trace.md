# A class-character test for the primitive nonsquare band

**Decision.** Park class and level enumeration as a route to classical polylogarithmic factoring. We can sum the entire primitive nonsquare band with a finite-field trace formula. On the two prime-discriminant controls, obtaining the inputs to that formula still visits every inverse pair of classes. A pair of public inputs also rules out replacing the full-order character information by all the lower-order traces on these controls. The ellipse reductions and the power-of-$N$ benchmarks remain useful independently. [Previous experiment and definition of the band](conductor_band.md)

## What a supplied class lets us compute

Assume $4N=c^2+\ell^2\Delta$, with $N=pq$ a product of distinct odd primes and $\gcd(c\ell,N)=1$. Here $\Delta>3$ is a prime congruent to 3 modulo 8, and the class group of discriminant $-\Delta$ is cyclic of odd order $h>1$. The command-line experiment chooses $\ell=1$ from $N$ alone. The small control set also uses supplied public branches. Prime discriminant gives odd class number. Cyclicity is a separate restriction.

For a reduced representative $F=[a,b,d]$, let $b_B(F)$ be the parity of primitive vectors up to sign with nonsquare norm $B/4<t\le B$, where $B\le\sqrt{\Delta/3}$. In this family the conductor is one, so a vector is primitive exactly when its two basis coordinates are coprime. Only four norms need testing:

$$a,\qquad d,\qquad a+b+d,\qquad a-b+d.$$

Indeed, completing the square gives $F(x,y)\ge\Delta y^2/(4a)$, and $a\le\sqrt{\Delta/3}$ excludes $|y|\ge2$. At $y=0$ primitivity gives $|x|=1$. At $y=1$, $|x|\ge2$ gives $F(x,1)\ge2a+d\ge2a+\Delta/(4a)>\sqrt{\Delta/3}$. Repeated norms from different directions are counted with their multiplicity. In particular $b_B$ is invariant under inversion of the class. The principal form has zero label throughout this band.

## One band sum in characteristic two

Write $G=\langle g\rangle$, $b(e)=b_B(g^e)$, and take a primitive $h$th root of unity $z$ over $\mathbb F_2$. Group nonzero indices $k\bmod h$ by $k\mapsto\pm 2^j k$. For one representative per orbit put

$$
d_k=\min\{d>0:2^dk=\pm k\bmod h\},\qquad
u_k(e)=z^{ke}+z^{-ke},\qquad
w_k=b(0)+\sum_{j=1}^{(h-1)/2} b(j)u_k(j).
$$

The ordinary Fourier inverse on the cyclic group becomes

$$
b(e)=b(0)+\sum_{k\text{ representatives}}
 \operatorname{Tr}_{\mathbb F_{2^{d_k}}/\mathbb F_2}\big(w_k u_k(e)\big).
$$

Indeed, $h=1$ as a field scalar, and inverse pairs cancel in the zero-frequency coefficient. Pairing opposite frequencies in the inverse gives $w_k u_k(e)$. Moving around the remaining orbit squares that product. Both factors belong to the displayed trace field. Since the orbit contains $2d_k$ indices, $\sum_k d_k=(h-1)/2$. This is exactly the number of nonprincipal inverse pairs. The identity removes the loop over norm levels in the contraction, but the generic list of field coordinates can still be as long as the class table.

## The coefficient that depends on the unfactored input

The public element $\alpha=(c+\ell\sqrt{-\Delta})/2$ has norm $N$. Its two prime ideal factors over $p,q$ have inverse classes $C,C^{-1}$. The four integral ideals of norm $N$ therefore have classes $1,1,C^2,C^{-2}$ with multiplicity. If $C^2=g^{e_*}$, summing the class character $\chi_k(g)=z^k$ over those four ideals gives exactly $u_k(e_*)$ in characteristic two. This is the norm-$N$ coefficient of the ideal-character theta series, using ideals or representations up to units. Evaluating the character on $(\alpha)$ alone always gives 1. The usual fast modular-form coefficient theorem takes the index in factored form. Its general factoring example uses an Eisenstein series and does not settle the complexity of this particular one-bit coefficient. [Theta series from class characters](https://www.math.mcgill.ca/darmon/pub/Articles/Research/79.DHRV/michigan.pdf) · [Bruin, Theorem 1.1 and Remark 1.3](https://pub.math.leidenuniv.nl/~bruinpj/coefficients.pdf)

Our public evaluator counts norm-$N$ representations up to sign for each reduced class up to inversion. For $[a,b,d]$ it solves

$$
w^2+\Delta v^2=4aN,\qquad w\equiv bv\pmod {2a}.
$$

It uses at most two progressions for $w$ from a square root modulo the public prime $\Delta$. Each class needs one root test and $O(1+\sqrt{aN}/\Delta)$ integer square tests. Choosing $v>0$ or $v=0,w>0$ gives the unit normalization. The transform of those counts yields all the $u_k(e_*)$ without factors or a mixed root modulo $N$. But this class-by-class count already identifies the hidden inverse pair. Building the public form catalogue takes exactly $\lfloor\sqrt{\Delta/3}\rfloor$ leading-coefficient probes, before generator selection, field construction and cutoff calibration. At $\Delta\asymp\sqrt N$, the probes alone cost $\Theta(N^{1/4})$. A reusable table is charged when $\Delta$ changes.

## Two transitions and two narrower obstructions

| Public $N$ | $\Delta$ | $h$ | Cutoffs | Band parity | Trace degrees | Leading-coefficient probes | Norm-$N$ square tests |
|---:|---:|---:|---:|:---:|:---|---:|---:|
| 138,468,065,327 | 1,410,683 | 311 | 578, 579 | 0, 1 | 155 | 685 | 5,476 |
| 39,866,375,420,117 | 1,639,243 | 119 | 280, 281 | 0, 1 | 24, 24, 4, 3, 4 | 739 | 31,908 |

One readout serves both cutoffs of a given $N$. For context, the earlier inert-prime sieve followed by targeted roots made 1,694 and 4,299 final square tests respectively at the larger cutoff. It also factors candidate levels. These are counts of one kind of operation, not timings of equivalent work. [Baseline measurements](conductor_band_results.json)

1. **A single proper quotient loses information.** The jump at $579=3\cdot193$ is supported on four oriented classes. The jump at the prime 281 is supported on two. A subgroup stabilizing either support acts freely on it, so its order divides both the support size and the odd class number. Its order must be one. If both cutoffs were functions of a common proper quotient of $G$ for every possible input class, the quotient kernel would stabilize their difference.
2. **Even the collection of all proper-order unordered traces can collide.** At $h=119=7\cdot17$, $N=39{,}866{,}375{,}420{,}117$ and $N'=1{,}120{,}192{,}977{,}083$ give the same $\Delta$ and the same unordered character values for every character of order 7 or 17. The orders recover local signs only up to independent reversal. At cutoffs 280 and 281 the answers are $(0,1)$ and $(0,0)$ respectively. Both trace vectors were recomputed from the unfactored inputs. An evaluator that also uses other information about $N$ lies outside this counterexample.

For an independent test we factored the *small candidate norms* on the control side, chose primitive orientations of their split prime powers, and accumulated them by class composition. All 860 class-and-cutoff predictions agreed. We also enumerated physical disk points at all 298 cutoffs on all 28 eligible branches in the repository's 194-branch sample. The sample includes the two preserved collision inputs and the replacement cohort of 192. Nineteen more branches have trivial class group and 147 lie outside the prime-discriminant family. In particular both historical collision branches have $\Delta=20$ and are covered by the earlier general-band tests. [Small-sample provenance](../ellipse_row_batches/README.md) · [Earlier general-band checks](conductor_band.md#results-and-controls)

## A one-bit test that avoids growing class fields

Fix $\Delta=59$, where $h=3$. The irreducible cubic $f=X^3+2X+1$ has discriminant $-59$. Its splitting field has group $S_3$ over $\mathbb Q$, and the quadratic subfield is $\mathbb Q(\sqrt{-59})$. At 59 inertia is a transposition (tame ramification and cubic discriminant exponent one), so the cyclic cubic extension of the quadratic subfield is unramified. The reduced forms give class number three, making it the Hilbert class field. At either prime factor of an eligible $N$, the cubic thus splits completely for the principal class and is irreducible for the other two classes. The order-three trace bit is 0 or 1 respectively.

We tried the public computation $\gcd(f,X^N-X)$ over $\mathbb Z/N\mathbb Z$ with modular exponentiation and monic Euclidean division, recording any nonunit pivot as a factor. On the first 128 eligible balanced controls from the deterministic generator, it gave degree zero every time and no divisor. The correct bits were 0 on 34 and 1 on 94. For instance $N=252{,}521$ has bit 0 and $N=31{,}877$ has bit 1. We checked the bits using both the public norm profile and the two prime Frobenius computations on the control side. Modulo $p$, even a root $r\in\mathbb F_p$ of the split cubic obeys $r^{pq}=r^q$; it need not be fixed by the exponent $N$. This experiment rejects that substitution on its own.

**When to reopen.** A concrete new test is a direct algorithm for this one-bit $\Delta=59$ problem on an unfactored promised input, with a polynomial-in-$\log N$ bit-cost argument. This would warrant testing a growing discriminant. The full reduction would still need affordable cutoff weights and trace fields, additional disk contributions, and enough accessible queries to decode a factor. An algorithm that obtains the disk bit directly from public data could bypass this representation altogether.

The [follow-up on the discriminant −59 bit](cubic59.md) gives a quartic with one versus zero or sixteen roots and checks public reciprocity and square tests on 320 inputs through 81 bits.

From the repository root with python-flint, SymPy and NumPy installed:

    python experiments/ellipse_disk_mod8/class_trace.py 138468065327 578 579
    python experiments/ellipse_disk_mod8/class_trace.py 39866375420117 280 281
    python experiments/ellipse_disk_mod8/class_trace_check.py --output rerun.json

[Evaluator](class_trace.py) · [Checks and control generator](class_trace_check.py) · [Results](class_trace_results.json)
