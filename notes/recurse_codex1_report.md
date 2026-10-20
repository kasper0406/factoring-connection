# Recursive shape tests in Dixon descent: a blind assessment

20 October 2026

The recursive substitution is mathematically sound as a reduction. Assuming the fixed-degree relation-supply and LV6 estimates needed by the proposed architecture, its optimal constant is

\[
\boxed{c_{\rm rec}(\theta)=\sqrt{(1+\sqrt2)\theta}},\qquad
\boxed{a_{\rm rec}=\sqrt{\frac{\sqrt2-1}{2}\theta}}.
\]

Random subproblem inputs, unbounded Las Vegas running times, and the ordinary factorization recursion do not invalidate this reduction. With growing degree, the corresponding size/smoothness model has minimum exponent **2/5 for one level** and **1/3 for several levels**. The latter requires a descent tree, including algebraic prime ideals, in a common field. Neither exponent is established as an unconditional factoring theorem by the notes or by the recursion argument below.

I read only the six requested notes, disregarding their assessment sections and proof-status labels. I did not inspect the excluded research files. Published sources were checked separately. In particular, I do not promote an internal theorem to an established result merely because it is cited by another internal theorem.

**Labels.** **PROVED** refers to a proof given here, or to an identified published result. A conditional implication is explicitly labelled as such. **HEURISTIC** identifies a model assumption about smooth values or descent supply. **OPEN** identifies the mathematical input still needed for an unconditional conclusion. The distinction matters: proving the recursive reduction does not independently prove all the analytic estimates underlying Theorem A.

## 1. Fixed degree

### 1.1 The precise input needed from the architecture

Write

\[
F(x)=\sqrt{\log x\log\log x},\qquad L_x=e^{F(x)}.
\]

Fix admissible \(d,k,\mu,\eta_0\), and put

\[
\nu=\frac1{4k},\quad
h=\max\{2\nu,\mu-(d-2)\nu,1-d\mu+\nu\}+\eta_0,
\quad \theta=2\mu+h+(d-1)\nu<1.
\]

Allow a fixed positive smoothness constant \(a\), instead of fixing it to the original roadmap's value. The following is the interface needed for the recursion theorem:

1. Preprocessing, field sampling, and all arithmetic outside shape tests cost the amounts stated below, uniformly for every sufficiently large input integer.
2. A field trial produces a row or a proper factor with probability \(L_x^{-o(1)}\). Retaining the first successful row-field permits collection of \(L_x^{a+o(1)}\) rows at an amortized cost of \(L_x^{o(1)}\) descent attempts per row.
3. A descent attempt factors at most two integers of size at most \(\sqrt{2x}\), followed by at most \(2k\) searches. Each search costs
   \[
   L_x^{s+o(1)},\qquad s=\max\{a,\theta/(2a)\}.
   \]
   Here the \(a\) accounts for the batch smoothness setup.
4. LV6 produces the required public square root, with choices independent of the stored hidden roots, in expected time \(L_x^{2a+o(1)}\).

Every \(o(1)\) must be uniform over input integers, with the fixed parameters held fixed. For the optimum below, \(a^2<\theta/2\), so \(s=\theta/(2a)\).

These hypotheses are weaker than accepting the whole Theorem A statement, but stronger than knowing only its final running-time bound. A black-box bound for the old algorithm at one particular \(a\) does not prove this interface at a new \(a\).

For clarity, the proposed analytic route to item 2 is as follows. With \(m=S^2\rho_r\rho_F\), one needs

\[
\mathbb EN_\ell\ge mL_x^{-o(1)},\qquad
\mathbb E[N_\ell(N_\ell-1)]\le m^2L_x^{o(1)},
\qquad m=L_x^{o(1)}.
\]

Paley–Zygmund then gives \(\gamma=\Pr(N_\ell>0)\ge L_x^{-o(1)}\). Fixed \(k\) split counting converts this to the row probability in item 2. Field rejections must have probability \(o(\gamma)\); proper-factor exits are successes. FM, Lemma 8.5′, LD*, and the LV6 adaptation are substantive inputs to this route, not consequences of recursion.

**Conclusion — PROVED, conditional reduction.** Under this interface, the recursive algorithm has the bound derived below. **OPEN here:** an independent unconditional certification of all the internal analytic inputs. This report establishes the reduction, not a new unconditional replacement for the published general factoring bound.

### 1.2 Random inputs and stopping times

Let \(T(n)\) be the expected time for complete factorization of the particular integer \(n\), including every auxiliary recursive call, and let

\[
\overline T(x)=\max_{2\le n\le x}T(n).
\]

Suppose a parent computation, after a history \(\mathcal H\), chooses an integer \(U\le\sqrt{2n}\). Give the child an independent random tape. Then

\[
\mathbb E[\text{child cost}\mid\mathcal H,U]=T(U)
\le\overline T(\sqrt{2n}). \tag{1}
\]

The distribution of \(U\) can be arbitrarily biased. No average-case assumption on integers being factored is used. For two adaptively chosen children the same argument applies successively.

A useful replacement for an informal invocation of Wald's identity is this elementary lemma. Let independent trials have pairs \((C_i,I_i)\), where \(C_i\ge0\) is the whole cost and \(I_i\) indicates success. The cost and success of one trial may be correlated. If \(\tau\) is the first success and \(p=\Pr(I_i=1)>0\), then Tonelli's theorem gives

\[
\begin{aligned}
\mathbb E\sum_{i=1}^{\tau}C_i
&=\sum_{i\ge1}\mathbb E[1_{\{\tau\ge i\}}C_i]\\
&=\sum_{i\ge1}\Pr(\tau\ge i)\,\mathbb EC_1
=\frac{\mathbb EC_1}{p}. \tag{2}
\end{aligned}
\]

Only independence of trial \(i\) from the event that the previous trials failed is used. A conditional version bounds the result by \(D/p\) when each trial's conditional mean cost is at most \(D\) and its conditional success probability is at least \(p\).

In a fixed field \(A\), let the row probability be \(q_A\), and let \(b=\mathbb E_Aq_A\). The first retained row-field has law

\[
d\mathbb P_{\rm keep}(A)=\frac{q_A}{b}\,d\mathbb P(A),
\quad
\mathbb E_{\rm keep}(1/q_A)=\frac{\Pr(q_A>0)}b\le\frac1b. \tag{3}
\]

Equations (1)–(3) apply even when shape-test costs have long tails. Acceptance depends on the factors, which are uniquely determined, and not on whether the factorization happened quickly. The calls run to completion; no timeout-based acceptance rule is needed.

Proper-factor exits need a small refinement. If a raw field trial exits with a factor with probability \(\delta\), delivers a row with probability \(b\), and \(p=\delta+b\), then the branch that retains a row-field is reached with probability \(b/p\). Its expected later cost is at most \(CD/b\) for \(C\) rows, so its unconditional contribution is at most \(CD/p\). Field selection costs at most \(D/p\). This avoids conditioning away successful factor discoveries.

Finally, after an observed trivial gcd, discard those rows and obtain fresh hidden roots. The conditional success probability of each completed round is at least \(1/2\), as proved in §1.4. Thus the expected number of rounds is at most two. Independence of a round's duration from its gcd outcome is unnecessary; apply the conditional version of (2).

**Conclusion — PROVED.** Expected-time bounds on every smaller input suffice. Neither Wald's identity with unjustified independence nor a tail estimate for the recursive shape call is required.

### 1.3 The recurrence, ordinary split recursion, and the \(o(1)\)

The expected local work of one SPLIT, counting auxiliary shape calls but excluding subsequent factorization of its output divisors, is bounded by

\[
D(n)\le L_n^{o(1)}
\left[L_n^a\overline T(\sqrt{2n})
+L_n^{a+\theta/(2a)}+L_n^{2a}\right]. \tag{4}
\]

The constant number of shape calls and the fixed \(2k\) searches are absorbed into the uniform subexponential error. For this upper bound one may charge every attempt for all its searches; a sharp shape/success ratio is unnecessary.

An ordinary factorization tree on \(n\), excluding auxiliary shape-test trees, has at most \(\lfloor\log_2n\rfloor\) prime leaves counted with multiplicity and fewer internal split nodes. Its random divisors need not be balanced. Conditional expectation and a uniform bound for each divisor show that the total local work is at most a polynomial in \(\log n\) times the maximum local bound for integers at most \(n\). Hence (4) yields

\[
T(n)\le A(n)
\left[e^{aF(n)}\overline T(\sqrt{2n})
+e^{(a+\theta/(2a))F(n)}+e^{2aF(n)}\right],
\quad \log A(n)=o(F(n)). \tag{5}
\]

For bounded small inputs, enlarge the constant in (5). The ordinary tree bound holds for prime powers and repeated prime factors as well, provided multiplicities are recorded and perfect powers are handled correctly.

Here is a genuine induction, rather than a substitution of an assumed asymptotic. Suppose

\[
c_0\ge 2a,\qquad c_0\ge a+\frac\theta{2a},\qquad
c_0\ge a+\frac{c_0}{\sqrt2}. \tag{6}
\]

Fix \(\varepsilon>0\), set \(C=c_0+\varepsilon\), and aim to prove

\[
T(n)\le K_\varepsilon e^{CF(n)}. \tag{7}
\]

For \(n>2\), \(\sqrt{2n}<n\), so strong induction is legitimate. Uniformly over the smaller inputs in (5),

\[
F(\sqrt{2n})=\left(\frac1{\sqrt2}+o(1)\right)F(n).
\]

The recursive term has exponent at most

\[
a+C/\sqrt2+o(1)
\le c_0+\varepsilon/\sqrt2+o(1)
=C-\varepsilon(1-1/\sqrt2)+o(1).
\]

The other terms have exponents at most \(c_0+o(1)=C-\varepsilon+o(1)\). Their strict margins absorb \(A(n)\), the finite sum, and every logarithmic tree factor for all sufficiently large \(n\). Choose \(K_\varepsilon\) to cover the remaining finite set. This proves (7). Since it holds for every \(\varepsilon>0\), for the same fixed algorithm,

\[
T(n)\le L_n[1/2,c_0+o(1)]. \tag{8}
\]

Finite expectation and almost-sure termination follow simultaneously by induction: each child input is smaller, and each parent renewal has positive success probability. A deterministic small-input cutoff starts the induction.

Shape-call chains have \(O(\log\log n)\) successive square-root reductions. Ordinary divisor chains can have \(O(\log n)\) nodes. The entire auxiliary recursion tree need not have only polynomially many nodes: its cost is already counted in (5). Charging a separate exponential “recursion-depth penalty” would double count that work.

**Conclusion — PROVED, under the interface of §1.1.** The complete factorization recursion, including arbitrarily unbalanced proper factors \(g,n/g\), has the bound (8). Polynomial overhead and recursion depth add no leading constant.

### 1.4 Transcript invariance and fruitfulness

Let \(H_N=\{\epsilon\in(\mathbb Z/N)^\times:\epsilon^2=1\}\). For odd \(N\), CRT gives \(|H_N|=2^{\omega(N)}\), including when \(N\) has repeated prime factors. Conditional on \(R=r^2\), a uniform unit \(r\) is uniform on its fiber \(rH_N\).

Fix all random tapes of every recursive invocation. Replacing a parent root \(r\) by \(\epsilon r\) leaves \(R\), the Gauss split, and both shape-call inputs unchanged. The recursive calls therefore execute identically: their own fields, roots, running times, gcds on their own moduli, and returned prime factorizations are unchanged. Their entire transcripts can be exposed without revealing the parent's fiber element.

This argument is pointwise in the tapes. It covers caching, input-dependent runtimes, and nested recursion. Use independent tapes, or disjoint independent portions of a master tape. The parent must not pass its stored root as an extra hint or use it to select relations. A canonical shortest-vector tie rule should depend only on \(R,N\).

Condition on the public transcript before the final gcd. Selection of a nonempty dependency \(S\) and its computed root \(Y\) uses only that transcript. Fix all hidden fiber elements but one belonging to \(S\). Then

\[
X_S/Y=\epsilon Z,
\]

with \(\epsilon\) uniform on \(H_N\), so \(X_S/Y\) is uniform on \(H_N\). Exactly two of its elements are \(\pm1\). Therefore

\[
\Pr(1<\gcd(X_S-Y,N)<N\mid\mathcal T)
=1-2^{1-\omega(N)}\ge\frac12. \tag{9}
\]

The public transcript here excludes the stored parent roots and the yet-unobserved final gcd. After observing a trivial gcd, refresh the rows before making another use of (9).

**Conclusion — PROVED.** Theorem DD's hypothesis (T) survives recursive shape tests using their own coins. This proof does not use a heuristic about independence between factoring difficulty and acceptance.

### 1.5 Exact optimization of the cost bound

For a fixed \(a\), (6) is equivalent to

\[
c_0\ge\max\left\{2a,\ a+\frac\theta{2a},\ 
\frac{a}{1-1/\sqrt2}\right\}.
\]

Put \(q=1-1/\sqrt2\). Since \(1/q=2+\sqrt2>2\), the linear-algebra term is redundant. Thus minimize

\[
\max\{(2+\sqrt2)a,\ a+\theta/(2a)\}.
\]

The two terms meet at

\[
a^2=\frac\theta{2(1+\sqrt2)}
=\frac{\sqrt2-1}{2}\theta.
\]

This point lies below \(\sqrt{\theta/2}\), where the second term is decreasing. To its right the first term is increasing; to its left the second is larger. Consequently the minimum is exactly

\[
\boxed{c_0=\sqrt{(1+\sqrt2)\theta}},\qquad
\boxed{a=qc_0}.
\tag{10}
\]

Equivalently, feasibility requires \(a\le qc_0\) and \(\theta\le2a(c_0-a)\). Maximizing the latter expression on \(0<a\le qc_0\) gives \(\theta\le(\sqrt2-1)c_0^2\), with equality in (10).

This is an optimum of the certified cost envelope, not a lower bound on the actual running time of every implementation.

For comparison,

\[
c_0<1\quad\Longleftrightarrow\quad\theta<\sqrt2-1.
\]

Thus the parameter threshold for beating the published constant 1 is unchanged. The recursive constant improves the old envelope in that regime. For larger \(\theta\), the pure recursive envelope can exceed 1. A separate robust variant can interleave the recursive algorithm and the published complete factoring algorithm, returning the first verified answer. That variant has an unconditional \(L_n[1/2,1+o(1)]\) fallback and, conditional on §1.1, constant \(\min\{1,c_0\}\). This does not prove the improved bound without the supply hypotheses. The published fallback is [Lenstra–Pomerance, Theorem 10.5](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf).

### 1.6 All parameter constraints

The change in \(a\) does not change \(\mu,h,\nu,\theta\). The relevant checks are:

| Requirement | Check at the recursive optimum |
|---|---|
| Odd fixed \(d\ge3\), fixed integer \(k>(d+1)/4\) | Unchanged. |
| Rational \(\mu\), \(\nu+\eta_0\le\mu<1\), \(\eta_0>0\) | Unchanged; rational approximations suffice when taking infima. |
| Kernel floor \(h\ge2\nu+\eta_0\) | \(H/V^2=N^{h-2\nu}L^{O(1)}\to\infty\). |
| FM level | \((\mu+h+d\nu)/2-(\mu+\nu)\ge\eta_0/2\). |
| Base-expansion/window floor | \(h\ge1-d\mu+\nu+\eta_0\), with the same polynomial margin. |
| Rational fibers and Möbius error | \(M/V=N^{\mu-\nu}L^{O(1)}\), \(\mu-\nu\ge\eta_0\). |
| Batch cost \(B\le S^2\) | \(a^2=(\sqrt2-1)\theta/2<\theta/2\). |
| \(S^2\ge\rho_r^{-2}L^{c'}\) for some \(c'>0\) | Exponent difference is \([h+(d+1)\nu]/(2a)>0\). |
| \(S<\sqrt\ell\), \(B^{C_0}<Y\), \(2S^2<Y\) | Fixed powers of \(L\) are smaller than \(N^{2\nu}\). Choose the fixed \(C_0\) needed by DS/L2_fam. |
| \(H\gg d\omega_0V\) for the FM truncation \(\omega_0=L^{C_1}\) | The positive fixed power \(h-\nu\) dominates any such fixed \(L\)-power. |
| Compact positive range for \(\log S/\log L\) | \(\sigma=\theta/(4a)\) is a fixed positive constant. |
| Prime multiplier below \(N\) | \(\mu+\nu+\eta_0/2\le\theta/2<1/2\). |
| Leading coefficient and polynomial height | \(\mu+h<\theta<1\); fixed constants are harmless asymptotically. |

The inequalities involving \(S\) hold with its computable choice \(S^2=L^{\theta/(2a)+o(1)}\). Its exact smoothness mass can be smaller than every inverse logarithmic power; only \(m=L^{o(1)}\) is used at fixed \(k\).

**Conclusion — PROVED arithmetic.** No geometric or analytic size restriction displayed in §1 of the roadmap forces \(a=\theta/\sqrt2\). That equality was the optimization imposed by the old shape tester. The constraint actually relevant to its batch-cost simplification is \(a^2\le\theta/2\), satisfied strictly by (10). The fixed-degree analytic estimates must still be available with their constants depending on this new fixed \(a\).

### 1.7 Minimizing \(\theta\) and numerical constants

At zero slack write

\[
f(\mu)=\max\{2\mu+(d+1)\nu,\ 3\mu+\nu,\
1-(d-2)\mu+d\nu\}.
\]

The first two affine functions increase and the third decreases. Their intersections give

\[
\begin{array}{ll}
\mu_*=(1-\nu)/d,&
\theta_*=\displaystyle\frac2d+\frac{d^2+d-2}{d}\nu,
\quad k\le(d^2+1)/4,\\[4pt]
\mu_*=[1+(d-1)\nu]/(d+1),&
\theta_*=\displaystyle\frac3{d+1}+\frac{4d-2}{d+1}\nu,
\quad k\ge(d^2+1)/4.
\end{array} \tag{11}
\]

For example, at the first intersection the increasing branch \(2\mu+(d+1)\nu\) dominates precisely when \(\mu_*\le d\nu\), equivalent to \(\nu\ge1/(d^2+1)\). This proves the division into cases; the formulas agree at the boundary.

For prescribed positive \(\eta_0\), the lower bound on \(\mu\) must also be imposed. The constrained minimizer is

\[
\widehat\mu=\max\{\nu+\eta_0,\mu_*\},\qquad
\theta_{\min}=f(\widehat\mu)+\eta_0,
\tag{12}
\]

provided \(\widehat\mu<1\) and \(\theta_{\min}<1\). Equation (12), or its infimum over rational \(\mu\), covers finite slack instead of silently optimizing outside the admissible set.

Consequently

\[
\inf_{k,\eta_0,\mu}\theta=\frac3{d+1},\qquad
\inf c_0=\sqrt{\frac{3(1+\sqrt2)}{d+1}}.
\tag{13}
\]

| Fixed \(d\) | Infimum \(\theta\) | Infimum recursive \(c_0\) | Corresponding infimum-limit \(a\) |
|---:|---:|---:|---:|
| 7 | 3/8 | 0.951488352998 | 0.278684486373 |
| 9 | 3/10 | 0.851037054841 | 0.249262982322 |
| 11 | 1/4 | 0.776886987015 | 0.227544930281 |
| 19 | 3/20 | 0.601774072519 | 0.176255545099 |

These are infima: each theorem uses finite fixed \(k\) and positive slack. In particular, writing one of the infimum constants with an \(o(1)\) for a single fixed parameter tuple would be inaccurate.

The infimum over fixed-degree theorems as \(d\to\infty\) is **0**. This does not imply any particular exponent smaller than \(1/2\). Under valid fixed-degree inputs, a computable portfolio of fixed-parameter algorithms can even yield \(\exp(o(F(N)))\): give algorithm \(j\) a fixed positive fraction of the work, then apply each fixed theorem before letting \(j\) increase. That diagonal statement supplies no quantitative power saving in \(\log N\).

**Conclusion — PROVED optimization, conditional complexity interpretation.** Equations (10)–(13) give the exact requested constants in the stated envelope.

### 1.8 What binds

At (10),

\[
T_{\rm shape}=L_N^{c_0/\sqrt2+o(1)},\qquad
T_{\rm SQ}=L_N^{\theta/(2a)+o(1)}
=L_N^{c_0/\sqrt2+o(1)}.
\]

Multiplying by \(L_N^a\) rows gives \(L_N^{c_0+o(1)}\). Linear algebra has exponent

\[
2a=(2-\sqrt2)c_0<c_0.
\]

At fixed parameters the split loss, good-field selection, and expected gcd-round count contribute no leading constant. The binding terms are the special-q searches and the recursive shape calls, balanced against the number of rows.

The original explicit \(1/\sqrt2\) shape cost came from applying the published constant-1 algorithm at size \(\sqrt N\). Recursion retains the geometric contraction \(F(\sqrt N)/F(N)\to1/\sqrt2\); it replaces the coefficient 1 by \(c_0\). Other square roots in optimized formulas come from smoothness optimization. They should not be treated as an additional compulsory shape-test cost.

**Conclusion — PROVED in the cost envelope.** There is no residual universal floor \(c_0\ge1/\sqrt2\).

## 2. Growing degree

### 2.1 A precise size/smoothness model

Put \(t=\log N\), \(l=\log t\), and

\[
E_\alpha=t^\alpha l^{1-\alpha},\qquad
B=\exp(bE_\alpha).
\]

In this section a smoothness model means that relevant, nonuniform polynomial values have the same leading smoothness exponent as random integers of their size, that the required special-q searches have adequate supply in one field, and that the associated losses are sufficiently uniform. These are **HEURISTIC assumptions** for the proposed descent family.

The numerical rule is elementary: if \(X=L_N[\xi,A]\) and \(y=L_N[\beta,b]\), \(\xi>\beta\), then

\[
\log\rho(\log X/\log y)^{-1}
=\left(\frac{A}{b}(\xi-\beta)+o(1)\right)
t^{\xi-\beta}l^{1-\xi+\beta}. \tag{14}
\]

Thus a smoothness search has exponent \(\xi-\beta\). This follows by inserting \(u=(A/b)(t/l)^{\xi-\beta}\) into \(\log\rho(u)^{-1}\sim u\log u\). Applying (14) to the actual polynomial-value distribution is the heuristic step.

Recursive testing has a different composition rule. If an algorithm has bound \(L_x[\alpha,c]\) for every integer \(x\), then completely factoring a candidate of size \(L_N[\xi,A]\) costs

\[
L_N[\alpha\xi,\ cA^\alpha\xi^{1-\alpha}+o(1)]. \tag{15}
\]

For \(\xi<1\), this is lower order than \(L_N[\alpha]\). For candidates as large as \(N^\rho\), the cost is instead \(L_N[\alpha,c\rho^\alpha+o(1)]\); a genuine contraction \(\rho<1\) is then essential.

The shape inputs \(\sqrt{2N}\) cost \(L_N[\alpha,c2^{-\alpha}+o(1)]\). Their recursion can be closed by a strict-margin induction exactly as in §1.3, whenever the other work leaves a finite leading constant.

### 2.2 One level: exponent 2/5

Let the initial prime bound be \(B_1=L_N[\beta,b_1]\), and take

\[
d\asymp(t/l)^\delta.
\]

Factoring the two Gauss halves and requiring their prime factors to be at most \(B_1\) has a split loss of exponent \(1-\beta\). Therefore a budget \(L_N[\alpha]\) requires

\[
\beta\ge1-\alpha. \tag{16}
\]

A regular two-dimensional special-q lattice has coordinate scale \(q^{1/2}\) times its search inflation. The logarithm of an algebraic value contains the two size terms

\[
\frac{t}{d}\quad\hbox{and}\quad d\log q.
\]

For \(q\le B_1\), their exponents are \(1-\delta\) and \(\beta+\delta\). Coefficient inflation for the level condition changes constants in the balanced regime. Minimizing their maximum gives

\[
\delta=\frac{1-\beta}{2},\qquad
\xi=\frac{1+\beta}{2}. \tag{17}
\]

Requiring the resulting values to be \(B\)-smooth with cost exponent at most \(\alpha\) gives

\[
\xi-\alpha\le\alpha,
\quad\text{or}\quad \beta\le4\alpha-1. \tag{18}
\]

Combining (16) and (18),

\[
\boxed{\alpha\ge2/5}.
\]

Equality has \(\beta=3/5\), \(\delta=1/5\), and norm-size exponent \(\xi=4/5\). Recursive factoring of these norms costs exponent \(\alpha\xi=8/25<2/5\). Hence replacing the testers cannot improve this one-level threshold: the supply of sufficiently smooth special-q values binds.

In particular, at a proposed \(\alpha=1/3\), the initial split forces \(\beta\ge2/3\); the best one-level norm exponent is then \(5/6\), whose \(B\)-smoothness cost has exponent \(1/2\).

**Conclusion — PROVED optimization within the model; HEURISTIC factoring prediction.** One level has minimum exponent \(2/5\), even with recursive shape and smoothness tests. A rigorous growing-degree factoring theorem remains **OPEN**.

### 2.3 The HECM barrier

The old test-composition obstruction was

\[
L_{B_1}[2/3]\rightsquigarrow L_N[2\beta/3],
\]

which at \(\beta=2/3\) gives exponent \(4/9\). A recursive call on a Gauss half has exponent \(\alpha\), with constant multiplied by \(2^{-\alpha}<1\), so it does not impose the inequality \(2\beta/3\le\alpha\).

**Conclusion — PROVED about the recurrence.** This particular tester obstruction disappears. The one-level size obstruction in §2.2 survives. Removing the former is not an unconditional construction of an \(L[1/3]\) factoring algorithm.

### 2.4 Several levels: lower bound 1/3 in the model

For a degree-one special-q ideal \(\mathfrak q=(q,\alpha-s)\), use the determinant-\(q\) lattice \(a\equiv bs\pmod q\). A relation has cofactor norm approximately

\[
\log X\lesssim\log\|f\|+
(d/2-1)\log q+d\log S+O(d\log d).
\tag{19}
\]

There is also a rational value; its logarithm has terms \(\log M+\tfrac12\log q+\log S\). Dividing a norm by one \(q\) does not remove the leading \(d\log q/2\) term when \(d\to\infty\).

At the bottom of a descent, a target just above \(B\) must be eliminated over primes at most \(B\). The generic size model has exponent at least

\[
\xi_{\rm bottom}\ge\max\{1-\delta,\alpha+\delta\}.
\]

Its smoothness cost can have exponent at most \(\alpha\) only if

\[
1-\delta\le2\alpha,\qquad \alpha+\delta\le2\alpha.
\]

Thus \(\delta\ge1-2\alpha\) and \(\delta\le\alpha\), implying

\[
\boxed{\alpha\ge1/3}. \tag{20}
\]

This is a limitation of the generic two-dimensional degree-one descent size model. It is not a lower bound for all factoring algorithms or all possible lattice constructions.

For \(\alpha>1/3\), choose \(1-2\alpha<\delta<\alpha\). A prime-size exponent \(\beta_j\) can then be reduced according to

\[
\beta_{j+1}\ge\max\{1-\delta,\beta_j+\delta\}-\alpha.
\tag{21}
\]

Above the bottom range the decrease can be \(\alpha-\delta>0\), so finitely many levels suffice for each fixed \(\alpha>1/3\).

At equality \(\alpha=\delta=1/3\), (21) alone does not decrease exponents. It would be a mistake to conclude that this excludes \(L[1/3]\): constants in the logarithmic prime bounds can decrease through \(O(\log\log N)\) levels.

### 2.5 Attaining 1/3 by constant contraction in the model

Set

\[
E=t^{1/3}l^{2/3},\quad D=(t/l)^{1/3},\quad
Q=DE=t^{2/3}l^{1/3},\quad d=\delta D.
\]

Use one field with

\[
\log M=mQ,\qquad \log\|f\|\le\kappa Q,
\qquad B=e^{bE}.
\]

The coefficient/window restrictions require suitable fixed constants, for example \(m\delta\gtrsim1\); one cannot choose the height independently of the degree and \(N\). Choose initial targets with \(\log q\le b_1Q\), fixed \(0<\lambda<1\), and search inflation \(S=e^{sE}\).

At a target with \(z=\log q\ge bE\), require its replacement factors to have logarithms at most \(\lambda z\). From (19), the two relevant logarithmic sizes have envelopes

\[
\log X_F\le\kappa Q+\frac d2z+dsE+o(dz+Q),
\quad
\log X_r\le mQ+z/2+sE+o(Q+z).
\]

Consequently, uniformly throughout the tree,

\[
u_F+u_r
\le\left[\frac{\kappa+m}{\lambda b}
+\frac\delta{2\lambda}+\frac{\delta s}{\lambda b}+o(1)\right]D.
\tag{22}
\]

Since \(D\log D=(1/3+o(1))E\), (14) assigns a bounded \(L_N[1/3]\) search constant at every level. Search volume and the \(dsE\) inflation must be solved consistently: one can choose \(b\) sufficiently large that the coefficient of \(s\) on the right side of the resulting search-budget inequality is smaller than the available volume coefficient \(2s\), then choose a finite \(s\). This avoids assuming that enlarging the box leaves the norms unchanged.

Replacing \(z\) by \(\lambda z\) reaches \(bE\) after

\[
O\left(\log\frac{b_1Q}{bE}\right)=O(\log\log N)
\]

levels. A node has at most polynomially many children in \(t\): the number of large prime factors is at most its logarithmic norm divided by \(\log B\). Therefore the entire tree has at most

\[
\exp(O((\log\log N)^2))=L_N[1/3,o(1)]
\]

nodes. Multiplying a bounded leading search cost by this tree size leaves the exponent \(1/3\).

There is an additional requirement if every candidate is completely factored recursively. At the top, \(\log X_F\) can be of order \((\delta b_1/2)t\). Choose constants with

\[
\delta b_1/2<\rho<1.
\]

The \(Q\)-scale height terms are lower order than \(t\); hence all candidate integers are eventually at most \(N^\rho\), or are explicitly rejected by that bound. Such a rejection must be included in the supply model. The leading recursive cost then contracts by at most

\[
r=\max\{2^{-1/3},\rho^{1/3}\}<1.
\]

If all nonrecursive work has a finite constant \(C\), an inequality of the form \(c\ge C+rc\) has a finite solution. Near the bottom, norm factoring has exponent \(2/9\) by (15). Without the strict upper bound \(N^\rho\), invoking recursion on polynomial-sized norms would be circular and could call the algorithm on larger inputs.

**Conclusion — PROVED size/budget construction under the model's supply assumptions; HEURISTIC complexity prediction.** Several levels can attain \(L_N[1/3,c]\) with a finite, unoptimized constant in this model. The minimum model exponent is exactly \(1/3\), not merely every exponent greater than \(1/3\). Its unconditional implementation is **OPEN**.

### 2.6 Which objects must be descended, and in which field

The first special-q steps concern rational primes dividing the Gauss halves. A lattice condition on the rational side forces that prime to divide \(G(a,b)\). Subsequent algebraic norms introduce prime ideals of the selected field.

Away from the leading coefficient, discriminant, and index, a prime ideal dividing a primitive linear element \(a-b\alpha\) has residue degree one: when \(p\nmid b\), its residue satisfies \(\alpha=a/b\), so its residue field is \(\mathbb F_p\); when \(p\mid b\), primitivity excludes divisibility. This proves the degree-one claim for the good primes actually arising from these linear elements. It does not say that every rational prime has a degree-one prime above it.

An algebraic target must therefore retain its root label \((q,s)\). Descending only the rational integer \(q\), while forgetting which ideal above \(q\) occurs, does not eliminate the corresponding ideal column. Rational large primes and algebraic degree-one ideals both need compatible relations. Bad primes require explicit treatment or a rigorously bounded residual square-class group.

The same field permits all these identities to be combined and then processed by a common square-root computation. Choosing a fresh unrelated field for each algebraic target does not provide identities in the original field. A multi-field construction would need a separate algebraic argument and a new LV6 analysis.

The good-target set could be defined as the primes whose *whole capped descent tree* succeeds. Initial split counting could then use that set. Proving that this set has sufficient density is much stronger than proving a marginal first moment for an independently chosen prime. Intermediate targets depend on the field and on previous successful norm factorizations; unconditional averages cannot be reapplied as if those targets were independent.

**Conclusion — PROVED algebraic requirements; OPEN supply theorem.** A common-field degree-one ideal descent is the relevant several-level model. Merely recursively factoring the integer norms does not supply the needed ideal relations.

### 2.7 The \(2^{\omega(N)}\) loss, and an elementary repair

Trial division through \(B=e^{bE_\alpha}\) gives only

\[
\omega(N)\le\frac{t}{\log B},\qquad
\log2^{\omega(N)}=O(t^{1-\alpha}l^{\alpha-1}). \tag{23}
\]

For \(\alpha<1/2\), this is larger than \(E_\alpha\). Therefore the original Theorem S/S′ bound with this loss does not certify the desired worst-case runtime. Equation (23) is an obstruction to that proof, not a lower bound on the algorithm's actual split probability.

Here the special shape of the accepted set permits a stronger elementary lemma.

**Disjoint-products lemma — PROVED.** Let \(G\) contain \(n\) primes coprime to odd \(N\). Let \(\mathcal S\) be all products of \(k\) distinct primes of \(G\), each at most \(X\), with \(2X^2<N\). Put \(h_N=2^{\omega(N)}\). If

\[
n\ge2h_Nk^2,
\]

then a uniform square \(R\) satisfies \(R\equiv u/v\pmod N\) for an accepted primitive pair with probability at least

\[
\boxed{\frac{\binom nk^2}{2\varphi(N)}}. \tag{24}
\]

*Proof.* Partition \(\mathcal S\) into the \(h_N\) cosets of the square subgroup. The number of ordered same-coset pairs is at least \(|\mathcal S|^2/h_N\). Two independently chosen \(k\)-subsets of an \(n\)-set intersect with probability at most \(k^2/n\), by a union bound. Thus at least

\[
|\mathcal S|^2/h_N-|\mathcal S|^2k^2/n
\ge|\mathcal S|^2/(2h_N)
\]

same-coset pairs have \(\gcd(u,v)=1\). Their ratios modulo \(N\) are distinct: equality implies \(uv'=u'v\), since the determinant has absolute value below \(N\); equality of positive reduced fractions implies equality of the pairs. Divide their count by the number \(\varphi(N)/h_N\) of squares. The primitive pair is also the shortest lattice vector up to sign: a shorter independent vector would have determinant of absolute value at most \(2X^2<N\), impossible in this lattice. Hence the actual Gauss membership test accepts it. ∎

The strict \(2X^2<N\) can always be ensured by shrinking the dyadic bound by a fixed factor.

For the growing-degree split sizes, \(\log q\asymp t^{1-\alpha}l^\alpha\), whereas (23) has one less factor of \(l\). If the good-prime density is at least a fixed inverse logarithmic power, then

\[
\log|G|\asymp t^{1-\alpha}l^\alpha\gg\log h_N+\log k.
\]

Thus the hypothesis of (24) holds for every sufficiently large preprocessed \(N\), including integers with many prime factors. The harmful \(2^{\omega(N)}\) loss is removed without trial division to \(L[1-\alpha]\).

**Conclusion — PROVED.** Keeping the old split bound creates a barrier below \(1/2\), but that particular barrier can be repaired for these large prime-product sets. It need not be the ultimate obstruction.

### 2.8 Precision of \(\gamma\)

Split counting pays a term comparable to

\[
k\log(1/g_A),
\]

in addition to combinatorial and prime-density losses of order \(k\log\log N\). At the boundary split scale,

\[
k\asymp(t/l)^\alpha,\qquad E_\alpha/k\asymp l.
\]

Consequently, sufficient quantitative requirements are

\[
\log(1/\gamma)=O(\log\log N)
\]

for a finite additional leading constant, and

\[
k\log(1/\gamma)=o(E_\alpha)
\]

if that loss is to disappear into the \(o(1)\). Thus a bound of the form \(\gamma=L_N[\alpha,-o(1)]\) is insufficient by itself. Inverse logarithmic density with an exponent growing uncontrollably with \(d,k\) is also insufficient.

This also affects the computable box size. If the moment bounds have a common loss \(P_N\), their useful form is

\[
\gamma\ge\frac1{P_N^3+P_N/m},\qquad m=S^2\rho_r\rho_F.
\]

A leading-term choice of \(S\) giving only \(m=\exp(o(E_\alpha))\) does not yield the required precision. One possible repair is to inflate the search volume by \(\exp(\varepsilon E_\alpha)\), with fixed positive \(\varepsilon\), so that a suitable lower smoothness bound certifies large \(m\). This costs a finite extra leading constant and avoids computing \(\Psi\) exactly. It works only if the moment estimates remain uniform for the inflated box and \(P_N\) itself has sufficiently small growth. It does not repair uncontrolled degree dependence in \(P_N\).

At \(\alpha=2/5\), \(k\asymp(t/l)^{2/5}\), \(d\asymp(t/l)^{1/5}\). At \(\alpha=1/3\), both can be of order \((t/l)^{1/3}\). Constants of size \(\exp(O(d\log d))\) inside the density estimate can become fatal after raising their reciprocal to the \(k\)-th power. They cannot be hidden as constants depending on the degree.

**Conclusion — PROVED precision requirement.** A sufficiently uniform density theorem is **OPEN**. This is a stronger demand than the fixed-parameter \(\gamma\ge L^{-o(1)}\) statement.

### 2.9 What must be proved uniformly

The following is a checklist of mathematical statements needed to upgrade the model. It is not supplied by substituting \(d=d(N)\) in fixed-degree notation.

| Internal ingredient | Required growing-degree statement |
|---|---|
| GB/GB_w and size ledger | Explicit degree dependence of coefficient bounds, discriminants, and bit costs; admissible shrinking slack; candidate sizes below the parent bound required by recursion. |
| Candidate/lattice lemmas 1.2, 1.3 | Uniform enumeration, regularity and recognition for rational targets and for labelled degree-one ideal targets at every level. |
| Lemma Kern | Explicit dimension dependence of fiber counts, lattice cells, wedge bounds and the height floor; adequate moment constants. |
| DS and smoothness estimates | Uniformity as interval lengths, \(u\), moduli and powers \(C_0\) change; losses small enough for §2.8. |
| Lemma 8.5′ | A restricted-modulus BV estimate for a large prime times a smooth part, with conductors and relative errors in the new scales. |
| FM | First moments for the actual special-q distribution, primitive/prime-coordinate filters, moving residues and all recursive levels. |
| W and L2_fam | Uniform multiplier spread, divisor switching and large-sieve bounds in the new target ranges. |
| LD* | Sufficiently precise second moments and descendability density, including the degree-dependent geometric constants. |
| Theorem K and field selection | Growing-\(k\) product counting with explicit constants, valid good-field selection, and whole-tree/adaptive target supply. The lemma in §2.7 supplies one useful part. |
| IR, content and derivative checks | Rejection probability negligible compared with the required density, with every proper factor exit counted as success. |
| LV6 | Uniform column count, exceptional square-class analysis, character sampling, sparse linear algebra, exact square root, and certificate heights/sparsity for entire trees. |
| DD and recursive testing | Independent tapes and root-invariant selection; strict subproblem size contraction. Their probabilistic argument is degree independent, but the supplied certificates must satisfy it. |

There are concrete range changes. In Lemma 8.5′'s fixed-\(k\) proof, a character nonprincipal on \(\ell\) has conductor at least a fixed power of the smooth-number argument, so Harper's large-conductor estimate applies. With growing \(k\), that exponent can tend to zero. New intermediate conductor ranges appear. In parts of a several-level tree, the condition \(Q\ge x^\varepsilon\) for fixed \(\varepsilon>0\) can also fail. An error estimate uniform only for fixed \(\varepsilon\) cannot simply be reused.

Likewise, L2_fam uses \(B^{C_0}<Y\) and a unique isolated large-prime divisor in its modulus partition. Near the bottom, the special prime is comparable with an \(L[\alpha]\) bound; those separations require new checking. At algebraic targets, rational-side coprimality arguments from the initial special-q search must be replaced by the appropriate ideal-side arguments.

For LV6, writing \((\log N)^{O_d(1)}\) is not a uniform algorithmic estimate. The fixed-degree proof's exceptional-zero window contains a factor \((2d)!\); its character-prime sampling must be rebalanced or analyzed quantitatively. Number-field square-root and matrix costs must include tree-expanded sparsity and coefficient heights. These are mathematical tasks, although none supplies a size-model exponent barrier beyond (20).

### 2.10 What Lee–Venkatesan already supplies

The relevant paper is Jonathan D. Lee and Ramarathnam Venkatesan, *Rigorous analysis of a randomised number field sieve*, **Journal of Number Theory 187 (2018), 92–159**, DOI [10.1016/j.jnt.2017.10.019](https://doi.org/10.1016/j.jnt.2017.10.019), [arXiv:1805.08873](https://arxiv.org/pdf/1805.08873).

Their degree grows as \(\delta(\log N/\log\log N)^{1/3}\). Theorem 2.5 proves relation supply for their ordinary randomized NFS; Theorem 2.6 and §6 provide growing-degree congruence extraction. Lemmas 5.13, 5.18 and 5.20 analyze coefficient randomization; Lemma 5.21 and §8 address good moduli. Lemma 8.5 is a quantitative estimate for smooth moduli. Lemma 6.10 handles exceptional square classes. These are genuine growing-degree precedents.

They do not prove the internal FM, LD*, K, or 8.5′ for this special-q family, nor an adaptive several-level descent theorem. Theorem 2.1 produces a congruence of squares; fruitful factorization is addressed conditionally in Theorem 2.3. These distinctions are explicit in [their statements](https://arxiv.org/pdf/1805.08873).

**Conclusion — PROVED bibliographic scope; OPEN transfer.** Their existing growing-degree results are valuable inputs, but do not fill the missing rows of §2.9 for this architecture.

## 3. Prior art

Three different meanings of “recursive factoring inside smoothness testing” should be separated:

1. Reapply a splitting method to the two factors it finds, until a smooth integer is completely factored.
2. Use a different factoring method as the smoothness-testing subroutine of an outer relation algorithm.
3. Call the entire outer relation algorithm on its smaller auxiliary integers, and close its own expected-time recurrence to improve a complexity constant.

The first two have clear rigorous precedents. I did not identify the third construction in the cited classical analyses. A negative literature search is not a proof of historical priority.

### 3.1 Lenstra–Pomerance 1992

H. W. Lenstra, Jr. and Carl Pomerance, *A rigorous time bound for factoring integers*, **Journal of the American Mathematical Society 5(3) (1992), 483–516**, DOI [10.1090/S0894-0347-1992-1137100-0](https://doi.org/10.1090/S0894-0347-1992-1137100-0); [paper](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf).

Algorithm 7.2 recursively applies ECM to discovered factors in its smoothness tester; Theorem 7.3 guarantees recognition of the specified recognizable smooth numbers. Algorithm 10.1 uses that tester in the class-group relation method. Algorithm 10.4 recursively factors the output divisors; Theorem 10.5 proves complete factorization. These are explicit instances of meanings 1 and 2. The inner smoothness tester is not Algorithm 10.1 called on its own auxiliary coefficient. No recurrence of the form (5) is used there.

**Conclusion — PROVED:** recursive factor extraction within a rigorous smoothness analysis predates this proposal. **OPEN historical claim:** priority for the stronger self-bootstrap.

### 3.2 Pomerance 1987

Carl Pomerance, *Fast, rigorous factorization and discrete logarithm algorithms*, in D. S. Johnson, T. Nishizeki, A. Nozaki and H. S. Wilf, eds., **Discrete Algorithms and Complexity**, Academic Press, Orlando, 1987, **119–143**; [author's scan](https://math.dartmouth.edu/~carlp/disclog.pdf), [publisher record](https://www.sciencedirect.com/science/chapter/edited-volume/abs/pii/B9780123868701500149).

This is the classical fast rigorous ECM smoothness-recognition input, including the treatment of suitable prime factors, used to improve random-squares factoring and index calculus. LP's discussion of its Theorem 7.3 specifically points to Pomerance's Theorem 2.1. The publisher states the rigorous random-squares factoring constant \(\sqrt2\). This establishes a factoring-subroutine precedent; I found no outer-algorithm self-bootstrap in this source.

**Conclusion — PROVED:** relevant rigorous smoothness-testing prior art. Identification with the stronger proposal is unsupported.

### 3.3 Seysen 1987 and A. K. Lenstra 1988

Martin Seysen, *A probabilistic factorization algorithm with quadratic forms of negative discriminant*, **Mathematics of Computation 48(178) (1987), 757–780**, DOI [10.1090/S0025-5718-1987-0878705-X](https://doi.org/10.1090/S0025-5718-1987-0878705-X); [paper](https://scispace.com/pdf/a-probabilistic-factorization-algorithm-with-quadratic-forms-4sxy4xmvw9.pdf).

Seysen's class-group algorithm is analyzed under ERH. Algorithm 4.2 uses trial division of reduced leading coefficients. §§6–7 introduce early aborts and Pollard–Strassen small-factor extraction. Theorem 6.2 separates relation-recognition cost from the outer algorithm; §6 explicitly suggests ECM as an improvement. This is not the outer algorithm factoring its own coefficients recursively.

Arjen K. Lenstra, *Fast and rigorous factorization under the generalized Riemann hypothesis*, **Nederl. Akad. Wetensch. Proc. Ser. A 91 = Indagationes Mathematicae 50(4) (1988), 443–454**, DOI [10.1016/S1385-7258(88)80022-2](https://doi.org/10.1016/S1385-7258(88)80022-2); [publisher record](https://www.sciencedirect.com/science/article/pii/S1385725888800222).

This combines Seysen's method with Pomerance's ECM smoothness test to obtain constant 1 under GRH. The outer and inner algorithms are different. The two volume numbers above are parallel journal-series numbering, not two different papers.

**Conclusion — PROVED:** rigorous conditional analyses with improved inner factoring tests; no identified instance of the stronger bootstrap.

### 3.4 Schnorr–Lenstra

C.-P. Schnorr and H. W. Lenstra, Jr., *A Monte Carlo factoring algorithm with linear storage*, **Mathematics of Computation 43(167) (1984), 289–311**, DOI [10.1090/S0025-5718-1984-0744939-5](https://doi.org/10.1090/S0025-5718-1984-0744939-5); [author-hosted paper](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1984c/artc.pdf).

Its random-class-group method exploits smooth class numbers and group operations to obtain ambiguous forms. Its proposed fast general runtime rests on a class-number smoothness assumption. This is different from recognizing smooth leading coefficients by a recursively invoked outer relation algorithm. It is not an unconditional precedent for the claimed recursive bound.

**Conclusion — PROVED scope of the precedent.** The superficial presence of class groups and smoothness does not establish the proposed self-recurrence.

### 3.5 Vallée

Brigitte Vallée, *Generation of elements with small modular squares and provably fast integer factoring algorithms*, **Mathematics of Computation 56(194) (1991), 823–849**, DOI [10.1090/S0025-5718-1991-1068808-2](https://doi.org/10.1090/S0025-5718-1991-1068808-2); [author-uploaded text](https://www.researchgate.net/publication/253163661_Generation_of_Elements_with_Small_Modular_Squares_and_Provably_Fast_Integer_Factoring_Algorithms).

The improvement comes from a rigorous near-uniform sampler of residues with squares of size \(N^{2/3+o(1)}\), combined with smoothness recognition. It yields the classical constant \(\sqrt{4/3}\). The smaller auxiliary-number size is relevant to the present cost discussion, but does not supply the recursive fixed point (10).

**Conclusion — PROVED:** earlier rigorous use of smaller auxiliary values. I found no self-bootstrap complexity argument in this work.

### 3.6 Lenstra–Pila–Pomerance

H. W. Lenstra, Jr., Jonathan Pila and Carl Pomerance:

- *A hyperelliptic smoothness test. I*, **Philosophical Transactions of the Royal Society of London A 345 (1993), 397–408**, DOI [10.1098/rsta.1993.0138](https://doi.org/10.1098/rsta.1993.0138); [author's scan](https://math.dartmouth.edu/~carlp/hyperI.pdf).
- *A hyperelliptic smoothness test. II*, **Proceedings of the London Mathematical Society (3) 84(1) (2002), 105–146**, DOI [10.1112/plms/84.1.105](https://doi.org/10.1112/plms/84.1.105); [paper](https://math.dartmouth.edu/~carlp/PDF/paper129.pdf).

Part I supplies short-interval smooth-number estimates; Part II supplies quantitative genus-two Jacobian-order results, particularly Theorem 1.1. They are the analytic basis for hyperelliptic small-factor methods. Part II points to a forthcoming third paper for further algorithmic details. Consequently a claim that these two papers alone contain every detail of an all-input smoothness recognizer should be checked separately. They do not establish the outer-algorithm self-bootstrap considered here.

**Conclusion — PROVED:** relevant inner-test prior art, with the scope of the published statements kept explicit.

### 3.7 Other relevant precedents

Daniel J. Bernstein, *How to find small factors of integers*, version **23 September 2002**, and *How to find smooth parts of integers*, draft **10 May 2004**, concern batch extraction and recognition. They are listed under “Recognizing and factoring smooth numbers” on [the author's factorization page](https://cr.yp.to/factorization.html). These are a separate way of reducing inner-test cost; they are relevant to Prop BS, rather than evidence for recurrence (5).

For several-level descent, Daniel M. Gordon, *Discrete logarithms in GF(p) using the number field sieve*, **SIAM Journal on Discrete Mathematics 6(1) (1993), 124–138**, DOI [10.1137/0406010](https://doi.org/10.1137/0406010), and Antoine Joux and Reynald Lercier, *Improvements to the general number field sieve for discrete logarithms in prime fields. A comparison with the Gaussian integer method*, **Mathematics of Computation 72(242) (2003), 953–967**, [paper](https://cr.yp.to/bib/2003/joux.pdf), provide the relevant NFS individual-logarithm context. The latter explicitly develops individual logarithms without a separate very large linear system. These are descent precedents, with heuristic smoothness/supply analyses; they do not give a rigorous all-integer factoring theorem for the present construction.

**Conclusion — OPEN historical priority for the specific stronger proposal.** Recursive factor extraction and factoring-based smoothness tests have established prior art. The narrower claim that the *entire outer factoring algorithm* was previously used as its own auxiliary shape test, with a rigorous improved self-consistent constant, was not verified in this search. It would be inaccurate either to call all recursion novel or to cite the classical papers as proving this particular recurrence.

## 4. Claims established, and the weakest point

The claims I consider **PROVED** are:

- The conditional recursive reduction, including random inputs, stopping times, field retention, proper-factor exits, the uniform induction, and complete-factorization recursion.
- Preservation of transcript invariance and the exact fruitfulness probability \(1-2^{1-\omega(N)}\).
- The optimal envelope \(c_0=\sqrt{(1+\sqrt2)\theta}\), its parameter feasibility, and the fixed-degree infima in §1.7.
- The one-level threshold \(2/5\) and several-level threshold \(1/3\) **within the explicitly stated size/smoothness model**, including the constant-contraction construction at the boundary.
- The disjoint-products split lemma removing the \(2^{\omega(N)}\) loss when its explicit cardinality condition holds.
- The need for much stronger quantitative density control when \(k\) grows, and the distinction between published growing-degree NFS results and the missing special-q descent results.

An unconditional improved all-integer factoring theorem is **OPEN on the evidence certified here**. At fixed degree the remaining qualification is the independently unverified supply/LV6 interface; at growing degree the additional uniformity and adaptive-target requirements are substantial.

**The single weakest point is a quantitative relation-supply theorem for the entire adaptive descent tree in one growing-degree field, with a descendability density strong enough to survive the growing initial target list.** Marginal smoothness estimates and recursive factorization of the tested integers do not establish that statement.
