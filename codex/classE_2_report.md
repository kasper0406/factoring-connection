# Referee report on `notes/classE_results.md`, v1

**Recommendation: MAJOR REVISION.**

The elementary identities and the main battery's reported z-scores mostly survive scrutiny. The conclusions do not. The draft asserts an invalid universal Heegner floor, promotes a truncated tail theorem to a larger class of path readouts, and calls a collection of classical Dedekind sums a test of genuinely three-dimensional content. Its original validity rule also fails in two tb cells. These problems require changes to the mathematical statements and experimental interpretation, rather than a bibliography update.

I treated the repository as read-only. Audit programs and their outputs are in this scratch directory. This report targets the v1 text initially read; the notes were being edited during the audit. Line references below refer to that initial text, and missing-file statements refer to the audited snapshot. “FATAL” means the stated claim must be withdrawn; it does not mean that the entire project is unsalvageable.

## Numbered findings

### 1. [MAJOR] Prop SPLIT proves a character identity, not an exhaustive classification (§0; §1.3, especially its “Meaning”; §4.1)

For **proper** indices `0 ≤ k < r`, the two asserted formulas are correct:

\[
(-1)^k s_k r_{k+1}\equiv a s_k^2\pmod N,
\qquad
\frac{N}{a_{k+1}+2}\le s_k r_{k+1}<N.
\]

The first yields the stated Legendre-symbol identity when the factors are units. The second follows from the continuant and remainder recurrences. Neither implication uses, or establishes, an exhaustion theorem for computable statistics.

In particular, “a split form carries no character constraint except one bit” is too broad. For **every** integer `d ≥ 2`, the same public path data satisfy

\[
W_{k,d}:=(-1)^k s_k^{d-1}r_{k+1}
\equiv a s_k^d\pmod N.
\]

Thus, for any local character `χ` of order dividing `d`, `χ(W_{k,d})=χ(a)` at unit steps. If `a=x^d mod N` with a known root, then `W_{k,d}=(xs_k)^d mod N` has a publicly computable power certificate. For example, `fermres=2^(N−1)` has a known **fourth** root throughout this population, because `N≡1 mod 4`. These identities need not provide an efficient new factor leak; they do refute the literal exclusion of higher-character constraints. Evaluating an unknown prime's higher character is a separate computational question.

There are also argument-specific quadratic identities. With `a=r=floor(sqrt(N))`, `D=N−r²`,

\[
r_{k+1}^2+D s_k^2
=N\big((s_k-rh_k)^2+D h_k^2\big),
\qquad (-D|p)=(-D|q)=1.
\]

This is already in cx1 §1.2. It is not a new independent genus bit, but the proposed classification must explain its scope rather than assert that no other constraint exists.

**Exact fix:** retain Prop SPLIT as an elementary identity. Replace “only two exact ways,” “no genus theory beyond,” “strictly weaker CFRAC,” and “covers their character content” with a restricted statement: *the displayed quadratic product relation is the usual congruence-of-squares mechanism; no exhaustive classification or algorithmic domination theorem is proved.* A square discriminant places the form outside ordinary nonsquare quadratic-order genus theory; that fact is not a computational lower bound.

### 2. [MAJOR] The size and CFRAC interpretation needs hypotheses; a principal-cycle sign is wrong (§1.3)

A lower bound depending on the next digit is not a uniform estimate `~N`. For `a=floor(sqrt(N))`, the `k=0` product is exactly `a~sqrt(N)` and its next digit is `~sqrt(N)`, deterministically. Large quotients are built into several structured arguments in the menu. A random-numerator large-quotient tail cannot be applied to them without proof.

Nor does the magnitude of an integer determine its smoothness probability: square factors, products of two separately accessible integers, correlations between steps, and relation yield all matter. No theorem comparing factoring complexity with CFRAC follows from this inequality. Even CFRAC's “size about `2sqrt(kN)`” should be an upper scale for its complete-quotient norms, not an assertion that all norms have that size.

The comparison also drops a sign. The principal norm identity gives

\[
((-1)^iQ_i|p)=1,
\quad\text{hence}\quad (Q_i|p)=(-1|p)^i,
\]

not `(Q_i|p)=1` individually. For the admissible product `N=103·107=11021`, `floor(sqrt(N))=104`, `Q_1=205`, and `(205|103)=−1`, whereas `(−205|103)=1`.

**Exact fix:** restore the sign everywhere. Qualify the `N^{1-o(1)}` relation-size discussion by an explicit bound on the relevant digits, and label the smoothness/runtime comparison heuristic. Do not infer “only smoothness can reveal a character” from the congruence.

### 3. [MAJOR] F-SQ controls the stated tail data, not arbitrary costs of the advertised full prefix (§2.2, lines 249–255, 279–284)

Let `A=r+t`, and let `u_j/v_j` be convergents to `s′/r`. The corresponding full convergents are

\[
\frac{v_j}{A v_j+u_j},
\qquad s_{j+1}=A v_j+u_j.
\]

The theorem's data include **only** `r mod m₀`, `t`, and the stopped tail word. They do not include `A` itself. Consequently they cannot evaluate an arbitrary class-E first-digit cost `F(A,γ)`. Only its residue-dependent part is covered. Likewise, exact remainders require additional labels such as `N mod m₀`; they do not follow merely from convergent-matrix residues.

For a concrete collision within one admissible window `I=[192,256)`, take `N=193·229=44197` and `N=197·233=45901`. With `m₀=1,Q=2`, both have `t=0` and stopped tail word `[2]`, while their first digits are respectively 210 and 214. Thus an arbitrary first-digit cost is not measurable from the theorem’s data.

The theorem must use a fixed function of its stated data, without unrestricted additional dependence on `N`. Allowing the latter would immediately destroy a target floor.

**Correct replacement:** put `X=P²` and choose integer parameters `m₀,Q` deterministically from the ensemble scale. For any fixed bounded function of

\[
(r\bmod m_0,\ t,\ \mathcal P_Q(s'/r)),
\]

where the overshooting digit is excluded, the two epsilon-class expectations differ by

\[
O(m_0Q^2X^{-1/6}\log X),\qquad m_0\le X^{1/6}.
\]

This is nontrivial when `m₀Q²X^(−1/6)log X=o(1)`. For fixed `m₀`, choosing `Q≤X^(1/12−δ)` controls the **tail word and the matrix states obtained from `A mod m₀`**, at full denominators of order `X^(7/12−δ)`. It does **not** control arbitrary costs of the large initial digit, all remainder labels, or arbitrary postprocessing using the public input. State the exact cutoff with `A v_j+u_j`, rather than identify it with `AQ`.

With the draft's weaker tb exponent, the corresponding scale is only `Q≤X^(1/24−δ)` and full denominators `X^(13/24−δ)`, not `X^(7/12−δ)`.

Finally, this is not the entire forward half of the tail fraction. That half extends to tail denominators about `sqrt(r)=X^(1/4)`, well beyond the useful `X^(1/12)` cutoff. The inverse-tail problem is open, but it is not the only uncovered part after the proved cutoff.

### 4. [MINOR] The F-SQ large-sieve calculation survives, but its proof and references need an actual statement (§2.2, lines 258–273)

The reduction uses a **comparison** step function, not an exact parametrization by `y`. On each of the `m₀` residue branches, `δ={m₀y}`; the map `{2δ}` traverses the unit interval twice. There are `O(m₀Q²)` partition boundaries because the stopped prefix is constant on a refinement by Farey fractions of order `2Q`. The displayed bound `m₀(|F_{2Q}|+3)` omits an inessential factor for the two traversals.

The displacement `δ²/r≤1/r` produces exceptional intervals of width `O(1/(m₀P))` in `y`, including wraparound and displaced `t` thresholds. The comparison `t₀=floor(2δ)` has its interior jump at `δ=1/2`; the exact threshold is `sqrt(r²+r)−r`. There is also a possible `t=2` layer near `δ=1`. On the unit semiprime population that endpoint case cannot occur, but it must be handled or explicitly excluded. Cylinder length alone does not guarantee agreement for every point: distance from the boundary is what matters.

For bounded product weights, write `M=π(I)`, `x_p=sqrt(p)−c₀`, `y_q=(h/m₀)sqrt(q)`. Absorbing `e(c₀y_q)` into the second weight is valid. Both close-pair counts are `O(M(1+m₀/h))`; these are counts of **pairs**, not assertions of prime-density estimates in short intervals. The one-dimensional double large sieve gives

\[
|S_h|^2\ll (1+hP/m_0)\,M^2(1+m_0/h)^2.
\]

In the stated range `hP/m₀≫1`. Dividing by `M²` **after taking the square root** yields exactly the draft's normalized estimate

\[
|E e(h\sqrt N/m_0)|
\ll \log P\,P^{-1/2}
\left(\sqrt{h/m_0}+\sqrt{m_0/h}\right).
\]

Erdős–Turán therefore gives

\[
D_c\ll H^{-1}+\log P\,P^{-1/2}
\left(\sqrt{H/m_0}+\sqrt{m_0}\right)+O(\log P/P).
\]

Set `H=P^(2/3)`, `m₀≤P^(1/3)`. This is `O(P^(−1/3)log P)`. Multiplying by the partition complexity and adding its boundary-layer mass proves the corrected epsilon statement in finding 3.

For class `u mod 4`, the denominator is `M_u(M_u−1)`, not literally `M²`. The fixed-modulus prime number theorem gives `M_u~M/2`. Removing the diagonal costs `O(1/M_u)=O(log P/P)` in bounded expectations. Sorting the two primes does not affect these symmetric data. Thus the two conditionings really are harmless here, but the justification should be given.

**Exact reference:** Bombieri–Iwaniec, *On the order of ζ(1/2+it)*, Ann. Scuola Norm. Sup. Pisa (4) **13** (1986), 449–472, **Lemma 2.4, pp. 452–453**. This is the required multidimensional inequality specialized to dimension one. [Original paper](https://numdam.org/item/ASNSP_1986_4_13_3_449_0.pdf).

**Exact fix:** quote that inequality and delete “Iwaniec–Kowalski §7.?”; I did not verify a matching numbered IK lemma and will not supply a guessed one. The original reference suffices.

### 5. [MAJOR] The tb approximation is repairable, but an unspecified trigonometric approximation is not a proof (§2.2, lines 274–276)

For unsorted primes, the correct symmetric region is

\[
|\log p_1-\log p_2|\le\log(8/7).
\]

Its indicator cannot be approximated uniformly by continuous trigonometric polynomials. A Lebesgue `L¹` error `O(1/T)` does not automatically give the same error for this discrete prime-pair distribution.

A repair is available. Use a specified majorant/minorant or a smoothed cutoff with controlled tails and a period containing the whole log-ratio range. A width-`η` strip at a boundary has at most `O(M(Pη+1))` prime pairs, by counting possible integers for the second prime. Its normalized mass is `O(log P(η+1/P))`; that bound also controls a suitably decaying smoothing error. Each Fourier term factors into prime weights of modulus one. Even the crude coefficient sum `O(T)` gives

\[
D_{tb}\ll T\,P^{-1/3}(\log P)^{O(1)}
+T^{-1}(\log P)^{O(1)}+P^{-1}(\log P)^{O(1)}.
\]

Taking `T=P^(1/6)` recovers the claimed weak rate `X^(−1/12+o(1))`, after normalizing by each target class's positive limiting probability. Better coefficient bounds may improve this, but are unnecessary for that statement.

**Exact fix:** insert the approximation lemma, discrete error estimate, coefficient norm, treatment of `N≡1 mod 4`, and class normalization. Keep the weaker coverage in finding 3 unless a stronger estimate is actually written down.

### 6. [FATAL] The random-class “every bounded shape statistic” floor is false (§4.2, especially lines 382–385; §0)

Weak equidistribution does not imply total-variation indistinguishability or agreement for every bounded statistic. For a fixed discriminant, different genus cosets are disjoint sets of exact reduced forms. Their exact shape coordinates recover `a,b`, and then `c=(b²−D)/(4a)`. A bounded function recognizing a genus therefore distinguishes the cosets perfectly. No Duke-type weak-convergence theorem contradicts this.

The coefficient residues are an additional, concrete problem. When `−4` is a genus discriminant factor, on forms with odd leading coefficient,

\[
\chi_{-4}([a,b,c])=(-4|a)=(-1|a),
\]

so a residue modulo 4 can itself carry the prescribed genus sign. There is no general theorem making the joint laws of `(a mod m₀,b mod m₀)` independent of every genus coset.

**Exact fix:** withdraw the stated consequence. Replace it by equidistribution against **fixed continuous or specified regular geometric tests**, with constants depending on the test and with any finite local components and admissible genus conditions retained. Do not conclude a computational floor for arbitrary coefficient statistics or the FW battery from it.

### 7. [MAJOR] The Heegner bibliography mixes three distinct results; nonfundamental orders and level are not automatic (§4.2)

The relevant checked statements are:

- **Duke (1988), Theorem 1(i), p. 75:** fundamental negative discriminants; fixed convex regions with piecewise smooth boundary; ineffective quantitative constants. This supports the fundamental-discriminant geometric statement. [Duke's paper](https://www.math.ucla.edu/~wdduke/preprints/hyperbolic.pdf).
- **Michel (2004), Theorem 3 and Theorem 10, pp. 189 and 224:** equidistribution after supersingular reduction / on **definite** quaternionic Gross-point sets. This is not itself the asserted archimedean modular-surface theorem. The genus factorization appears explicitly in §6, equation (6.2) and its ensuing discussion. [Michel, Annals 160, 185–236](https://annals.math.princeton.edu/wp-content/uploads/annals-v160-n1-p05.pdf).
- **Harcos–Michel (2006), Theorem 2, equation (1.4):** maximal-order Heegner points on `X₀(1)`, for subgroup orbits of index `o(|D|^(1/23042))`, against fixed continuous tests. Taking the square subgroup gives the needed genus-coset geometric result when its index is sufficiently small; for fundamental discriminants it is `|D|^{o(1)}`. Their **Theorem 6, Remark 6.1, and §6.4** give the specified Shimura-curve / squarefree modular-level extensions, under the Heegner condition. These do not say that arbitrary coefficient congruences can simply be appended. [Harcos–Michel, Invent. Math. 163, 581–655](https://www.renyi.hu/~gharcos/rankinselberg.pdf).

Write `D=f²D_K` and specify primitive forms/proper ideal classes before discussing a nonmaximal order. In the actual FW battery `D=−4N` **is fundamental**, since `N` is squarefree and `N≡1 mod 4`. For squarefree fixed `k`, `k=1,5` also give fundamental `−4kN`; `k=3` instead gives conductor `f=2` over discriminant `−3N`. It is incorrect either to apply the maximal-order theorem verbatim to all `−4kN`, or to imply that the actual `k=1` battery needs a nonfundamental extension.

Nonmaximal-order equidistribution involves finite-place conductor factors; a useful precise account is Einsiedler–Lindenstrauss–Michel–Venkatesh, **§2.3, equations (3)–(4), pp. 8–9**, which explicitly separates the field twist estimate from Clozel–Ullmo's local estimate for nonmaximal orders. This explains the additional input; it does not supply the draft's full joint genus-and-residue assertion. [ELMV](https://people.math.ethz.ch/~einsiedl/ELMV3.pdf).

A level lift must specify its CM embedding, local orientation and admissible orbit. `Γ₀(m)` and `Γ(m)` are different covers; the quoted squarefree-level statement is not a theorem for every principal level, especially at ramified primes such as 2. For a nonfundamental genus-coset theorem with the precise arbitrary level/residue readout in the draft, the supplied citations do **not** constitute a checked reference.

**Exact fix:** use the fundamental, level-one, regular-test statement that applies here. State any order/level extension separately with its hypotheses and an actual supporting theorem or proof. Replace the informal product-of-L-values equality by the correctly normalized **squared Weyl-sum** formula, including local factors and the Eisenstein contribution. Ordinary Gross–Zagier concerns a derivative/height formula and should not be used as a substitute for that period formula.

### 8. [MAJOR] “Q-leaks are cusp excursions” confuses discriminant signs and arithmetic events (§4.2, lines 385–388)

For a positive-definite form, `Im τ=sqrt(|D|)/(2a)` is correct. For a threshold `a≤A` corresponding to a cusp height `Y≥1`, normalized hyperbolic mass is `3/(πY)`, hence proportional to `A/sqrt(|D|)`. This is a statement about a **region**, not the mass of a single coefficient.

The eps/leakmap principal `sqrt(kN)` paths have **positive** discriminant and correspond to closed geodesics, not those negative-discriminant Heegner points. Small coefficients can have a related cusp interpretation for geodesic arcs, but it requires the appropriate correspondence and sampling measure. Moreover, squareness and smoothness are arithmetic conditions, not equivalent to being high in a cusp. Fixed-region equidistribution does not give relative estimates in the shrinking regions relevant to very small coefficients.

**Exact fix:** distinguish definite Heegner points from indefinite geodesics, and smallness from smoothness/squareness. Retain the geometric analogy only with those qualifications. Short initial power walks `f^j`, `j≤32`, are not large subgroup orbits covered by Harcos–Michel.

### 9. [MAJOR] The entire advertised shifted d=3 toy family collapses to classical Dedekind sums (§3.5; `toy3.py:48`)

The draft itself gives

\[
F_3(n;a,b,c)-F_3(n+c;a,b,c)=F_2(n;a,b).
\]

Every toy tuple uses `c=1`, `n=1`. Therefore

\[
\begin{aligned}
F_3(1;a,b,1;N)
&=F_3(0;a,b,1;N)-F_2(0;a,b;N)\\
&=-\frac{N-1}{8N}
+\frac12\{s(ab^{-1},N)-s(a,N)-s(b,N)\}.
\end{aligned}
\]

All 16 toy statistics are postprocessing of three ordinary Dedekind sums and public data. They are not a test of genuinely three-dimensional shifted content. I independently checked this and the following identity in **356 exact cases**, including even moduli and multiples of 3:

\[
F_2(1;a,1;N)=-\frac{N-1}{4N}-s(a,N).
\]

Thus the eight pairs of `dr1_*` main-menu features also add no new shifted-sum primitive. The FD2 pairs with arbitrary `a,b` do not all collapse by this particular one-step argument.

**Exact fix:** relabel the toy as a classical Dedekind-combination check. Withdraw the inference about genuinely shifted three-dimensional content. A new experiment must choose shifts/arguments avoiding short telescoping to unshifted lower-dimensional sums, with that algebraic audit performed before labels are examined.

### 10. [MAJOR] The claim that only the spatial local part factors is false (§3.4, lines 346–349)

The two decompositions are different, as already stated explicitly in cx1 §3.1. But failure of the **same modular decoder** to apply to the weighted split does not imply failure to factor from that split.

For a direct counterexample, let

\[
L_{\rm freq}=s(1,p)/q+s(1,q)/p,\quad H=p+q.
\]

Using `12M s(1,M)=(M−1)(M−2)`,

\[
X:=12NL_{\rm freq}=H^2-2N-3H+4,
\qquad
H=\frac{3+\sqrt{8N-7+4X}}2.
\]

So this actual frequency local part factors `N` immediately. I checked the identity at `p=103,q=107`, where `X=21432` and the decoder returns `H=210`.

**Exact fix:** say *the unweighted Rademacher decoder does not automatically decode the weighted local oracle*. Do not claim only one decomposition factors. Credit cx1 §3.1 for the Dedekind spatial/frequency distinction itself, not merely a related observation.

### 11. [MAJOR] The toy result is incomplete and contains unreported factor exposures (§0, §3.5; prereg §6)

At this audit, only `classE_toy3_n20.jsonl` and its summary exist. The `n=22` raw and analysis files advertised in §0 are absent, and §3.5 literally contains `TOY3_PLACEHOLDER`.

The `n=20` data contain all **324** products in the registered population; I independently enumerated and matched them. The missing `n=22` population has **1089** products. For `n=20`, recomputed maxima are `1.6405445` for epsilon and `1.8023630` for tb. These support the limited reported n=20 non-detection.

The toy routine silently zeroes nonunit tuples. There are two such products:

- `N=674413=769·877`, with `floor(N^(2/3))=7690=10·769`;
- `N=727393=773·941`, with `hash2=286783=371·773`.

These are direct public gcd exposures. The raw toy records do not log them as auxiliary splits. In particular, a general “size-impossible for the other structured arguments” rationale is false for the two-thirds argument.

**Exact fix:** report only completed populations, disclose those zero-coded exclusions/exposures, and remove the n=22 claim until results exist. Complete the section. An exhaustive population has no Monte Carlo sampling uncertainty; its permutation z is a descriptive randomization score whose null requires justification. Prime reuse does not by itself prove the error is always anti-conservative. None of this repairs the algebraic collapse in finding 9.

### 12. [MAJOR] D1 is a justified first-order correction, but it cannot be called a pass under the original registered rule (§5; §9)

For the planted construction `x=y` with probability `(1+ρ)/2`, write `S=Σ_h S_yy,h`. Conditional on the observed labels,

\[
E[U]=\rho S,
\quad
E[V]=\frac{1-\rho^2}{4}S
+\rho^2\sum_h\frac{S_{yy,h}^2}{n_h-1}.
\]

Consequently `2ρsqrt(S)` is the correct **leading approximation** to the expected z. It is not an exact expectation of the random ratio `U/sqrt(V)`. At `ρ=.05` the omitted terms are too small to explain the original discrepancy; the variance adjustment is substantively right.

My raw recomputation gives:

| n | planted tb z | corrected first-order expectation |
|---:|---:|---:|
| 48 | 16.100760 | 18.179780 |
| 64 | 16.751711 | 18.209328 |
| 128 | 18.255111 | 18.198817 |

The registered `ρsqrt(200000)=22.36068`, with tolerance 3, fails at 48 and 64. The amended expectation passes. This is a defensible mathematical bug correction, not evidence of arbitrary threshold tuning, but it remains a **post hoc amendment** to a decision rule.

**Exact fix:** report original validity status and amended status separately. Replace “Conclusion (registered rule)” for those two cells by a qualified conclusion under disclosed amendment D1. A prospective fresh run can restore a fully prospective amended validity claim. There was no registered obligation to replicate a nonexistent hit; the issue is the changed validity rule.

### 13. [MAJOR] The tb MDE confuses planted flip strength with residual correlation (§5, lines 417–423)

For each feature, set `Sxx=Σ_h S_xx,h`, `Syy=Σ_h S_yy,h` and define the pooled residual correlation `R=U/sqrt(Sxx Syy)`. Exactly,

\[
z=R\sqrt{\frac{Sxx\,Syy}{V}},
\qquad n_{\rm eff,feature}=\frac{Sxx\,Syy}{V}.
\]

The quantity `4Syy` describes the attenuation of the particular planted **binary flip** construction. It is not a universal effective sample size for every feature's residual correlation. In the tb raw data, for example, `isqrt2_ms_mod3` has feature-effective sample size about **198397** at n=48, rather than `4Syy≈132202`. Its nominal normal-approximation MDE for `R` is about `.01536`, not `.0188`. Other features can have different heteroskedastic weighting or sparse support.

The `.0188` number may be conservative for many features, but it is not the advertised common correlation calibration. Also, 80% detection power is not a confidence statement excluding every effect at that scale. A null result at a detection threshold does not logically “refute” an effect.

**Exact fix:** define the effect parameter. Report feature-specific power/intervals for residual correlation, or label the amended number as an MDE for the specified planted flip-strength model. Use “no registered hit” unless an actual equivalence or simultaneous confidence analysis is supplied. Keep the existing warning that numeric residue codes are only tested linearly; that warning is necessary.

### 14. [MINOR] Top-12 size strata are a reproducible operational floor, not a theorem removing archimedean information (§5)

The registered floor is legitimate as the definition of the test. It does not remove within-bin dependence on exact size, nor establish a floor conditional on arbitrary public residues jointly with size. Top-16 and joint `(top-8, N mod 16)` checks are useful sensitivity analyses, but the latter has **coarser size resolution** than top-12 and is not a refinement of it. Finer strata also lose observations/information through small or pure-label cells.

The draft's sensitivity results are reassuring at the tested resolution. My independent test of linear normalized size within top-12 strata gives z `−1.1235`, `.6899`, `−.6663` at n=48,64,128. This supplies no evidence of that particular residual size confound; it is not a proof that every nonlinear size proxy has been removed.

**Exact fix:** state the conclusion conditional on the specified coarse floor. If claiming a broader archimedean floor, use a prespecified geometry-aware smooth baseline and joint residue adjustment, with calibration and power reported for that analysis.

### 15. [MAJOR] Several nulls are forced or nearly forced by the feature construction (§5; `battery.py`)

Beyond the classical-sum collapse in finding 9:

- **`fw_minbits` mostly measures the chosen starting prime.** It equals the bit length of the least odd split prime `ell` in 199998/200000 n=48 records, and in **every** n=64 and n=128 record. A walk beginning at a very small coefficient makes its running minimum uninformative about later coefficient smallness. The two n=48 exceptions reach coefficient 1.
- **`fw_sq` includes unreduced squares by construction.** At n=128, all **99860** records with starting prime `ell=3` have `fw_sq=16`: even powers among the first 32 give square leading coefficients before reduction becomes necessary. This is not a generic random-class square-frequency test.
- **`div3` is a bounded-residue/size family.** If `a=floor(N/3)−1`, then `N=3a+d`, `d∈{4,5}` for the present large prime products. Its remaining Euclidean computation has bounded modulus after one large quotient. Several nonconstant statistics merely encode fixed public residues, despite being numerically nondegenerate. In these data `div3_len`, `div3_n_a1`, `div3_rem1mod4`, `div3_phimod5`, and `div3_dmod7` are all constant within `N mod 4200` classes.
- The n=20 toy `isqrt_invsqrt_Z3sign` is `−1` on **all 324** records. Its zero z is a degeneracy, not evidence against association of a substantive feature.

These are useful controls if labeled correctly. They weaken the claim that the menu broadly probes independent mechanisms.

**Exact fix:** list algebraic redundancies and conditional degeneracies separately from informative tests. For future FW tests, use endpoint/window features after the forced small initial portion, or condition on the starting prime and separately record unreduced powers. Do not relabel old tests as independent confirmation after inspecting their outcomes.

### 16. [MAJOR] The verdict still asserts more than either the proofs or the experiment establish (§0; §6)

A defensible verdict is:

> No main-menu scalar association crossed the registered threshold in the supplied datasets. The elementary VIS/SPLIT identities and a precisely scoped F-SQ tail theorem are supported. A full structured-path floor, an exhaustive character-channel classification, genuinely three-dimensional shifted-sum evidence, and an unrestricted reduced-coefficient floor are not established. Two tb validity cells require amendment D1.

Specific further fixes:

- The random-prefix law is independent of `N` **at fixed cutoff `T`**. Substituting `T=T(N)` changes the reference law with `N`; include that public scale in the baseline or choose a common dyadic cutoff. The residue-labeled extension also retains its public `N mod m` baseline. Therefore “beyond nothing at all” needs qualification. This was already explicit in cx1 §2.5.
- An arbitrary coefficient readout that retains `N` is computationally universal only in the unrestricted sense that it can run a factoring algorithm. Do not imply a proved polynomial-time completeness reduction.
- The Coppersmith comparison needs a nontriviality condition: recovering a small root of `a−c−x≡0 mod q` factors only when the resulting gcd with `N` is proper. For `a=1,c=0`, the root `x=1` is already a global integer root and the gcd is `N`; mere closeness modulo `q` does not always factor.
- The `N^(2/3)` extension is expressly an unwritten sketch. It should not inherit the PROVED label of the square-root theorem.
- FD3-0 is already **cx1 equation (37)**; the spatial/frequency Dedekind distinction is already **cx1 §3.1**. Preserve those credits. “No prior formulation found by us” is not a novelty proof for these elementary consequences or the general weighted Barvinok method.
- The shifted d=2 integer is not automatically the ordinary **Barkan–Hickerson–Knuth** functional for arbitrary shifts. Polynomial-time floor-sum evaluation does not prove membership in a fixed-state additive grammar. State the endpoint/shift state actually retained; use BHK directly only in the classical cases, including the collapsed `dr1` case above.

### 17. [MINOR] Correct the formal and reproducibility details

- `L_a={(s,r):as≡±r mod N}` is a **union of two lattices**, generally not a lattice. Use `L_a={(s,u):as≡u mod N}` with signed remainder coordinate `u=(-1)^k r_{k+1}`. Modulo a factor this is one line. The unsigned union is two lines.
- The actual 130 roles are **117 `menu`, 5 `menuG`, 5 NC, 2 epsilon PCs, 1 tb PC**. The “113 menu ... plus 4 degenerate NC” accounting is wrong. Among the 117 menu features, the numbers numerically nondegenerate under the registered floors are 105/103 for epsilon/tb at n=48, 103/103 at n=64, and 100/100 at n=128.
- The named aggregate-congruence test does **not** test all even moduli and multiples of 3: it explicitly keeps odd `N` coprime to 6. The cx1 proof covers all moduli, and my independent 240 cleared-congruence checks cover these excluded cases, but the test attribution/count should be truthful. The SPLIT test currently uses `fermres` and a random argument, not the two exponent arguments claimed in the results text.
- The AST isolation test checks direct imports only. `battery.py` imports `modsym`, which imports `sympy`. This does not demonstrate factor-label leakage, but the test does not prove transitive isolation from everything in its banned list. Say “direct-import guard and N-only entry point.”
- The raw n=48 run contains **199998 distinct products** among 200000 records. This is consistent with the sampler, which does not deduplicate. Prime reuse is substantial at 48, slight at 64, absent in the supplied 128 sample. Independent draws of pairs remain independent under the stated sampling law; reuse alone does not invalidate their population inference. A prime-disjoint design would answer a different generalization question.
- Registered P4's heuristic `8·40·2/sqrt(N)` is neither the exact VIS probability nor negligible after multiplying by 200000 at n=48 (it predicts order 8 events). Use the actual uniform-unit VIS probability for the randomized reference and separate arguments for structured inputs.
- Remove `TOY3_PLACEHOLDER`, `CODEX_PLACEHOLDER`, and `CRITIC_PLACEHOLDER` before presenting this as a completed results document.

## What I verified correct

1. **VIS and the elementary CF structure:** the exact cylinder qualification, absence of a proper q-divisible convergent denominator in the balanced range, the unit square-root gcd obstruction, the product-size inequality, signed character identity, and inverse/reversal conventions are consistent with the supplied proofs. Do not extend the square-root obstruction to the two-thirds or hash arguments.
2. **The epsilon F-SQ analytic core:** the centered double large sieve, its spacing counts, prime-pair normalization, Erdős–Turán exponent, Farey partition complexity and boundary-layer method work after the statement/proof repairs above.
3. **FD3-0 and the cleared congruences:** the signs and constants in the displayed unshifted formula, `12MF₂` and `24MF₃` congruences, and cx1 Theorem 4's proof are correct. My additional checks found **0 failures in 240 cases**, including even moduli and multiples of 3.
4. **The unweighted DSF decoder:** `D_d N(F_d(p)+F_d(q))≡C_d(p+q) mod N` is correct. If `C_d` is a unit, the least residue recovers `p+q`; a proper gcd of `C_d` with `N` already factors; `C_d≡0` gives no conclusion. The weighted oracle is different, and can sometimes factor too. I checked 24 small DSF congruences exactly.
5. **Exact computation:** the floor-sum evaluator and the fixed-dimension weighted-lattice reduction are sound at the stated scope. The O(N) toy implementation is exponential in input bit length; the polynomial Barvinok route remains unimplemented. The supplied test file passes: **11 tests, 0 failures**, run with bytecode and pytest cache writes disabled.
6. **Main raw-data arithmetic:** each supplied main file has 200000 records. Independent grouped covariance/variance recomputation reproduces the reported maxima:

| n | epsilon max `|z|` | feature and signed z | tb max `|z|` | feature and signed z |
|---:|---:|---|---:|---|
| 48 | 3.146240 | `cf3_sumbits`, −3.146240 | 3.027202 | `isqrt2_ms_mod3`, −3.027202 |
| 64 | 2.273956 | `cbrt2_phimod5`, −2.273956 | 2.782637 | `cbrt2_n_a1`, −2.782637 |
| 128 | 2.616632 | `isqrt_rem1mod4`, +2.616632 | 2.745399 | `div3_midbits`, +2.745399 |

I also recomputed the complete operational feature vector from `N` for the first, middle and last row of each main file: **9/9 exact matches**. The auxiliary records show no main nonunit arguments, no rational-path exposures, and exactly three n=48 form-walk split records. These checks support “no observed main-menu threshold hit”; they do not support the broader mathematical exclusions.

Audit artifacts: `audit_classE.py/.json`, `algebra_checks.py/.json`, `shift_collapse_checks.py/.json`, and `degeneracy_checks.py/.json`, all confined to this scratch directory.
