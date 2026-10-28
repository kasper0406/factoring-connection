# audit2L: GPT-6.1-sol blind referee report (verbatim)

Prompt: `codex/audit2L_1.md`. Run: `codex/run.sh audit2L_1 xhigh`, 10:37–11:08 CET, 2026-10-28. The raw output is
`codex/audit2L_1.out`; it is reproduced below unchanged. Its finite counterexample script is
`codex/scratch/audit2L_1/referee_checks.py` (git-ignored scratch). Codex reports that a heading-based filter showed it one
inline legacy sentence in escapeR2 §7.3, and that it did not use it.

---

(i) \(3/8<\alpha<2/5\): **MAJOR REVISION**.

(ii) \(\alpha=3/8\), including \(a>2c_2\): **MAJOR REVISION**.

Locations below refer to the blind excerpt. I found a precision-critical gap in the parent rational block. I did **not** establish a counterexample to the final complexity claim.

| Severity | Location | Issue | Effect on exponent/constant | Required fix |
|---|---|---|---|---|
| **MAJOR** | B §8.3, [lines 1597–1607](/Volumes/git/factoring-connection/notes/audit2L_blind_claims.md:1597) | The proof obtains \(\mathrm{polylog}\,(\delta_e+\mathcal V_0)\), then asserts \(\lambda^C(\delta_e+\mathcal V_0)\). These are different precision classes. | The claimed exponent is unproved; this cannot be absorbed into \(c(\alpha)\). | Prove parent normalization and weighted straddling estimates with **only powers of \(\lambda\)** lost. Retain the \(X_F/q\) density throughout. |
| **MINOR** | B Lemma EXC(c), [lines 1532–1537](/Volumes/git/factoring-connection/notes/audit2L_blind_claims.md:1532) | \(s_i\) contains mid-prime factors, but the stated covering argument treats only a large \(P_0\)-smooth part. | No exponent change after adding the missing Euler-product bound. | Include the mid-prime part in the Rankin argument explicitly. |
| **MINOR** | B Lemma PRIM, [lines 1555–1557](/Volumes/git/factoring-connection/notes/audit2L_blind_claims.md:1555) | For parents, “\(f'f''=0\) unless \(p\le y\)” is false: a common content prime can belong to \(J\). It need not divide \(K\). | Apparently negligible, but the parent proof is not verbatim. | Bound the missing \(J\)-prime content contribution. |
| **MINOR** | A §13, [line 816](/Volumes/git/factoring-connection/notes/audit2L_blind_claims.md:816) | W-P(P3) is an aggregate lower bound, not the pointwise lower bound cited in the weight bridge. | No apparent exponent change. | Cite and apply the pointwise smoothness lower bound from OPV-A step 6, then perform the acceptability subtraction. |

**The single weakest link: RAT-2’s parent normalization.**

The decisive displayed line is

> \[
> \sum_{P,P'}c_Pc_{P'}r_P
> \le (m^{\mathrm{par}})^2\cdot\mathrm{polylog}\cdot(\delta_e+\mathcal V_0).
> \]

The excerpt explicitly defines “polylog” as \(e^{O(\lambda)}=(\log N)^{O(1)}\), **not** \(\lambda^{O(1)}\).

Put \(v=3-7\alpha\). Even granting the proposed child bound \(A\ll D^{-2\alpha}\lambda^C\), the displayed parent argument only gives, for an unspecified fixed \(K\),

\[
B_{\rm rat}\ll D^K\lambda^C
 \left(D^{-2\alpha}+D^{-2v}+\delta_e\right).
\]

With \(\eta\asymp D^{-v}\), the available bound for the mixed decoupling term is consequently

\[
\frac{\sqrt{AB}}{\eta}
 \ll D^{K/2-\alpha}\lambda^C.
\]

Already **one extra factor of \(\log N\)**, \(K=1\), leaves \(D^{1/2-\alpha}\lambda^C\). At the endpoint this is \(D^{1/8}\lambda^C\), rather than a quantity tending to zero. Enlarging \(A_0\) can suppress the \(\delta_e\) term; it cannot remove the loss multiplying \(\mathcal V_0\).

There is a concrete normalization issue behind this bookkeeping. Parent smoothness is root-conditioned, so its scale is \(X_F/q\). HT’s local ratio gives

\[
\frac{\widehat\rho_{B_2}(X_F/q)}
     {\widehat\rho_{B_2}(X_F)}
 =q^{1-\beta_2}\bigl(1+O(1/u_2)\bigr).
\]

For \(q\asymp B_2\), this factor is of order

\[
B_2^{1-\beta_2}\asymp u_2\log u_2,
\]

which contains a power of \(D\). Thus one cannot simply relabel the “polylog” in SPR as \(\lambda^C\) while keeping its stated \(S_{\rm typ}\), based on \(\rho_2(X_F)\). A §9 correctly introduces \(\rho_a(q)=\widehat\rho_{B_2}(X_F/q)\); B’s parent argument must do the same before normalizing.

A sufficient repair would prove, uniformly in labels, environment cells and relevant pair types, bounds of the form

\[
c_P\le \lambda^C C_q,\qquad
m^{\rm par}_{q,e}\ge \lambda^{-C}n_pC_q,
\qquad
C_q=R_{\rm typ}\widehat\rho_{B_2}(X_F/q),
\]

together with a corresponding weighted straddling bound. These look approachable from the stated inputs. They are nevertheless stronger than the bounds actually supplied. W-P’s \(D^{\pm C}\) bounds suffice for exponentially small character errors, but do not close RAT-2.

**1. Cell covariance and one-point moments**

The central covariance device survives the basic algebraic checks.

For an interior algebraic cell, CF(ii) gives

\[
\#\{\text{cell points with coordinate 1 smooth}\}
 =n(R,\delta)d_1(1+\mathrm{tiny}).
\]

For rational cells, TYPE’s diagonal invariance gives constant residue-fibre degrees, and FIB supplies the same identity with its inclusion–exclusion error. Consequently, weighting an interior cell by its actual lattice-point count gives the correct marginal identity. There is no missing fibre-length factor **provided the stated count weights and the individual \(Z_y\) denominators are retained**.

With the atom \(*\), the covariance identity is also correct. Its contribution to the centered second moment is \(\widehat\nu(*)c_i^2\). Boundary terms can be bounded using moments over the enlarged convex region. Replacing those moments by maxima would introduce an unacceptable spread factor; the displayed Hölder argument avoids that.

For OPV-A, the important ingredients are consistent:

- CL handles full prime-power classes correctly.
- The local factor has exactly \(E_{\pi_0}\tau=1\).
- Archimedean variation and the FO error enter the second moment quadratically.
- The stated squared \(\beta\)-drift errors fit the budgets below \(2/5\).
- \(h^{1/2}=D^{-5}\) and \(e^{-\lambda}=1/(D\lambda)\) are sufficiently small.

OPV-A step 6’s revised lower-bound argument is valid in structure:

\[
E_{\widehat\nu}(\Psi_1/c_1)
 \ge 1-\widehat\nu(*)-\mathcal V^{1/2}.
\]

It uses the exact interior single-count formula and does not require PU-A or a boundary smoothness estimate. I found no surviving first-order term in this algebraic covariance calculation.

There is, however, a missing case in rational EXC(c). For example, take

\[
P_0=5,\quad S_0=7,\quad y=17,\quad q=19,
\]
\[
P_1=(1,1),\quad P_2=(1,210),\quad g=(1,23).
\]

Then \(K=209=11\cdot19\), and the values are \(22,-187\), both \(17\)-smooth. Their \(P_0\)-parts are \(2,1\), but their type factors are \(s_1=22,s_2=11>S_0\). A large small-prime part does not cover these types.

This omission is repairable: for \(\sigma=1/(2\log P_0)\), the additional mid-prime Euler product is bounded by

\[
\prod_{\substack{p\mid K\\P_0<p\le y}}
 (1-p^{-1+\gamma+\sigma})^{-1}
 \le \exp\!\left(O(\omega(K)P_0^{-1+\gamma+\sigma})\right)
 =1+o(1).
\]

But that argument should appear in the proof.

**2. Character errors**

Drappeau’s Lemme 5 has the claimed \(1/\varphi(k)\) weighting, conductor restriction, and effective constant for \(A<1\). The L¹-to-L² reduction used in EN/EN-A preserves that normalization. [Drappeau, Lemme 5](https://arxiv.org/pdf/1307.7554).

The required conductor ranges fit:

| Block | Order of \(\log k\) | Order of \(\log x\) |
|---|---:|---:|
| Rational children | \(D^{4\alpha-1}\lambda\) | \(D^{2\alpha}\lambda\) |
| Rational parents | \(D^{1-\alpha}\lambda\) | \(D^{2\alpha}\lambda\) |
| Algebraic children | \(D^\alpha\lambda\) | \(D^{2\alpha}\lambda\) |
| Algebraic parents | \(D^{1-\alpha}\lambda\) | \(D^{2-3\alpha}\lambda\) |

Thus \(k=x^{o(1)}\), including after the stated endpoint, divisor and primorial reductions. The lower condition \((\log x)^{c_D}\le y\) also has ample margin.

WA itself is the correct elementary inequality. Its content is precisely the domination quantity \(\sum_j\max_k\varphi(k)H_j(k)\).

I found no additional exponent obstruction in the family accounting:

- Rational moduli retain \(q\), allowing the divisor counting over possible window primes.
- The comparison \(Q_2\max_q\widehat W_q\ll\mathrm{polylog}\sum_q\widehat W_q\) is legitimate for the explicit dominating weights on a dyadic window.
- Algebraic moduli are \(\Delta'/\delta\); they need not retain \(q\). FAM-A’s different summation argument is appropriate.
- DOM’s pointwise upper bounds only need polynomial losses here.
- Separate parent thresholds are necessary and give the stated saving.

Indeed,

\[
\theta=2\alpha-\tfrac12<\alpha,\qquad
\theta_p=1-2\alpha<3-7\alpha
\]

throughout the claimed interval. Exponential character savings absorb the polynomial losses in these arguments. They cannot absorb the polynomial loss in RAT-2’s **main covariance term**.

**3. Thin pairs and skewed labels**

The multiplicity argument avoids a divisor-function loss. The child thin-pair fraction \(O(s/\vartheta_0)\) follows from counting partners with a fixed determinant.

For parents, fixing each permitted prime \(b'\) gives the required fraction

\[
O(1/\vartheta_0^p)+O(q/V_1);
\]

one should use this count, rather than compare all lattice points in a strip against a prime-window denominator.

The losses in PUB/PU-A are affordable for these exponentially rare sets. Applications of BTS to intervals crossing zero need the usual split, or an explicit additional \(\Psi(W,y)\) term; this does not change the loss exponent.

AB §10’s skew-label repair is sound in principle. Its second moment and lower mean are both expressed through the **same sum of per-point weights**. It does not require a comparison between the shortest and longest points of a skewed box. Combining that bound with LS and DOM gives the claimed negligible skew contribution once \(C_T\) is sufficiently large.

**4. Parents**

The hard filter requires the separate treatment given in RAT-2. Inactive parents can have zero density, so the child relative formulation cannot be reused unchanged. Centering at \(\pi_Pc_P\), and charging straddling parents, is the right approach. The unresolved step is the normalization identified above.

W-P(P1)–(P3) provide enough precision for FAM-A’s parent character errors. They do not supply the stronger normalization needed in RAT-2.

The root-conditioned acceptability subtraction has the correct critical scale:

\[
\frac{D^C}{B\,\rho_2(X_F)^2}
 =
 \exp\!\left[-aD^\alpha\lambda
       +2c_2D^{3-7\alpha}\lambda(1+o(1))\right].
\]

Thus \(a>2c_2\) is a sufficient endpoint constraint for this crude subtraction. The weight bridge only divides by one smoothness density, so its endpoint requirement is weaker.

The parent version of PRIM also needs an additional case. Its quoted vanishing assertion fails, for example, with

\[
y=7,\quad J=(20,40],\quad q=13,
\]
\[
P_1=(4,5),\quad P_2=(3,7),\quad g=(23,23).
\]

Both points have the same label modulo \(13\), \(K=13\), and their rational values are \(-23,-92\). Both belong to \(23\mathcal S_7\), although the common content prime \(23>y\) does not divide \(K\). Hence EXC(b) does not already remove this case.

The missing contribution should be exponentially small: the content probabilities contribute \(\sum_{\ell\in J}\ell^{-2}\), while the density enhancement costs powers of \(\log Y\). That estimate must be added explicitly.

**5. Assembly, endpoint and cost**

The elementary probabilistic reductions check out:

- UB is an incidence subtraction followed by Paley–Zygmund.
- DC is Cauchy–Schwarz under the root-conditioned law.
- DC\(_e\) follows by conditioning and summing.
- VD is total variance.

None requires independence between different children.

Lemma \(\eta\)’s ratio argument avoids the fatal normalization loss of the older proof. On the long interval pieces, the relevant prime sum is

\[
\sum_{B<q\le B_2}\frac1q
 \frac{\rho(u-\log q/\log B_2)}{\rho(u)}
 \le (1+o(1))u+O(\lambda).
\]

The omitted near-zero mass is small enough even after the deterministic child cap. The short-interval range is appropriate. The vague log-concavity attribution can be replaced by a precise published reference: Balister–Bollobás–Morris, Theorem 3.3, attributes it to Hildebrand’s proof of Lemma 1. [Published statement](https://annals.math.princeton.edu/wp-content/uploads/annals-v188-n1-p02-p.pdf).

For U′ under DW, I found no independent exponent failure in the range restatement. At every relevant scale,

\[
u_2\asymp D^{3-7\alpha}\gg\lambda^4.
\]

The cited U1₂, DS-U′ and CS′ arguments use local HT estimates, fixed-factor shifts and Rankin bounds; they do not intrinsically need the particular threshold \(D^{1/4}\). The class-(B)/(C) margins and Kern-U’s height condition also remain adequate.

At \(3/8\), the new DW choice gives

\[
\eta_0\log N\asymp D^{5/8}\lambda,
\]

so its TV terms are negligible against the \(L[3/8]\)-scale densities. The restrictions arising from the older design with \(\eta_0=D^{-1/2}\) should not be imported into this design.

The exponent arithmetic is consistent:

\[
(3\alpha-\tfrac12)-(4\alpha-1)=\tfrac12-\alpha>0,
\]
\[
(3-7\alpha)-(1-2\alpha)=2-5\alpha>0.
\]

The two-level feasibility window closes below \(3/8\), while the parent margins and degree branch become borderline at \(2/5\).

If the claimed second-order bounds were established, PB would have ample margin even at the endpoint:

\[
A/\eta\ll D^{3-9\alpha}\lambda^C,\qquad
\sqrt{AB}/\eta\ll D^{-\alpha}\lambda^C.
\]

The later computational accounting also appears consistent, conditional on the missing relation supply:

- field rejection is negligible compared with a polylogarithmic supply;
- two-level row cancellation preserves the required identity;
- LV6’s \(\exp(O(d\log d))\) overhead remains below \(L[\alpha]\);
- canonical factorization outputs preserve transcript invariance;
- factoring parent values with Theorem U costs \(L[\alpha]^{o(1)}\) when \(\alpha>4/11\);
- the square-root recursion gives the stated denominator \(1-2^{-\alpha}\).

Accordingly, the conditional cost formula is consistent:

\[
c(\alpha)=\max\!\left(
2a,\ a+s+\max(a,\sigma_1+\sigma_2),\
\frac{a+s}{1-2^{-\alpha}}
\right).
\]

The stated \(\sigma_1,\sigma_2\), including the additional \(c_2\) in \(\sigma_1(3/8)\), agree with the leading value sizes. But the present proof of RAT-2 does not establish the supply needed to invoke this accounting.

**6. Sanity**

The historical benchmark is not a contradiction: the assumed internal Theorem U already goes beyond the published Lenstra–Pomerance bound. [Lenstra–Pomerance’s theorem](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf).

These checks reveal no fundamental obstruction to two-level descent at \(L[3/8]\). They reveal an incomplete precision argument. The cell device can square first-order errors; it does not permit an unspecified \((\log N)^C\) loss afterward. A rewritten parent RAT-2 proof with uniform, root-conditioned normalization is required before accepting either theorem.

The finite examples and independent exponent checks are reproducible in [referee_checks.py](/Volumes/git/factoring-connection/codex/scratch/audit2L_1/referee_checks.py).

*Blindness limitation: a heading-based exclusion filter inadvertently displayed one inline legacy verdict/confidence sentence in escapeR2 §7.3. It was not used in these findings.*