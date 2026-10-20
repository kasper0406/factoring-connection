REFEREE REPORT: nested_results.md
20 October 2026

Recommendation: MAJOR REVISION. The Haar-average calculation in AF is sound. The numerical backtracking calculation and the proposed recursion constant are sound within their stated models. The large-children lower bound is wrong, several advertised barriers exceed what is proved, and NR is not yet a complete conditional theorem as written. Its main algebraic mechanism is repairable without FM, LD*, or Theorem K. The revision received during this review addresses irregular target lattices by enlarging the boxes; this changes the correct descent exponent from δb1/2 to δb1. The remaining main NR gaps are the admissible-field probability and the cost of its specified Hensel-prime selection.

Scope. I read nested_results.md in full, thmA_proof.md, recurse_results.md §§1 and 3, lv6_results.md §§2–3 and the relevant field-check discussion, nested_prereg.md, the summary JSON, and the E1 analysis and raw data. The note changed during the review. Issues 1–13 describe the 577-line snapshot saved as nested_results.reviewed.txt, SHA-256 74f2a9674356ccf381e2c62a07238d6539ffa0e0f611ae7631bdef74d47ad623. I then checked the full diff against the 609-line revision saved as nested_results.updated.txt, SHA-256 72020164d603f2de5de7e258c0bd341c0b4406bcdb3b06843874fafa59a9b5eb. Issue 14 explicitly credits and assesses those fixes and supersedes the affected earlier findings. The verdict register below reflects that later revision. References use section and lemma names to avoid moving line numbers. I made no repository edits.

I fetched the original dBPW v2 PDF through the web PDF reader and checked Theorems 7.1 and 9.5 and the chain (34)–(50). A separate attempt to save the PDF with curl failed because the shell could not resolve arxiv.org; there is no local PDF binary. Source: [dBPW v2](https://arxiv.org/pdf/2512.01588v2). I also checked the GRH estimate against [Lagarias–Odlyzko, Theorem 1.1, retypeset original](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf), and the character-cost remark against [Lee–Venkatesan, Remark 6.7](https://arxiv.org/pdf/1805.08873).

Verdict register, corresponding to the seven requested items:

1. AF's displayed inequalities and auxiliary Gamma(r1+r2−1) integral: CORRECT. The literal cofactor-norm Gamma law and tail gloss: WRONG. Corollary AF and the ray-condition remark: CORRECT WITH FIXES. The unrestricted sampler/impossibility conclusions: WRONG.
2. SI: CORRECT WITH FIXES in its application. FB's local lemma: CORRECT WITH FIXES; the earlier universal conclusion was WRONG and is now qualified. R: CORRECT WITH FIXES. R2: GAP as a claimed general floor, including within the revised counting model unless its missing assumptions are added. OD's fixed-point linear dependence: CORRECT. The children-count claims: WRONG. The broader method exclusions: GAP.
3. BT: CORRECT WITH FIXES as a restart model; it is not a general descent lower bound.
4. X_ERH: CORRECT WITH FIXES. The revised text correctly calls its character primes poly(t)-bit primes and repairs the 1/(5d) margin. The correction to recurse §3.5: CORRECT at the L-scale claimed.
5. NR: GAP as stated. RW, TS, and T_multi: CORRECT WITH FIXES. Split acceptance and revised descent sizes: CORRECT WITH FIXES / CORRECT, respectively. The recursion inequality and c_NR are CORRECT with the revised τ_F=δb1, once the remaining implementation and field-law gaps are repaired. Revised AD is a precise, strong hypothesis; the reduced-basis box addresses the earlier geometric defect, but its joint random-integer model remains an assumption.
6. E1: CORRECT WITH FIXES. P3 is honestly reported as failed; the numerical table and preregistered bootstrap intervals reproduce exactly. The independence interpretation and one rate label need correction.
7. Verdict and Candidate card: WRONG in their universal barrier language. Retain the conclusion that this work establishes no ERH-only improvement, with the qualifications below.

I use FATAL for a defect destroying an advertised central theorem with no identified repair, MAJOR for a substantive invalid inference, missing hypothesis, or algorithmic gap, and MINOR for a local correction. I find no demonstrated unrepairable FATAL obstruction to the conditional NR mechanism. This does not certify NR as presently proved.

1. MAJOR — AF's displayed comparison is correct; its literal Gamma cofactor-norm law and tail gloss are wrong.

Write r=r1+r2−1≥1, X=ρ^n, w=|μ_K|, and A0=2^{r1}(2π)^{r2}X/√|Δ_K|. The box-count estimate gives e^{−1/4}A0≤M(a)≤e^{1/4}A0 uniformly in a. For a fixed integral ideal c, let U_c(a) be its generator count in the box. Then p_c(a)=U_c(a)/M(a). Uniform bounds on M, rather than an unjustified ratio-of-expectations argument, give

  e^{−1/4} E[U_c]/A0 ≤ E[p_c] ≤ e^{1/4} E[U_c]/A0.

Use weighted logarithmic coordinates Log(u)=(n_ν log|u_ν|)_ν, n_ν=1 or 2. The relevant simplex has volume √(r+1)(log(X/Nc))_+^r/r!, while the log-unit torus has volume √(r+1)R_K. Thus the two square-root factors cancel, and

  E[U_c] = w (log(X/Nc))_+^r / (h_K R_K r!).

The analytic class-number formula cancels the factors w,h_K,R_K,√|Δ_K| and 2^{r1}(2π)^{r2}, leaving C(Nc)/(ρ_K X), with C(v)=(log(X/v))_+^r/r!. In particular, (2π)^{r2} is correct in the paper's metric, whose lattice covolume is √|Δ_K|. Mixing it with the unscaled complex-coordinate covolume 2^{−r2}√|Δ_K| would be incorrect.

For A_S(u)=#{c∈S:Nc≤u}, summation by parts gives

  Σ_{c∈S} C(Nc)/(ρ_K X)
   = [1/(ρ_K X(r−1)!)] ∫_1^X A_S(u)(log(X/u))^{r−1} du/u
   = ∫_0^{log X} δ̃_S(Xe^{−v}) v^{r−1}e^{−v} dv/(r−1)!.

The boundary at X vanishes; the norm-one ideal, when present, is included by the Stieltjes jump at 1. This proves AF. The Gamma shape is r, not n, r+1, or r1+2r2−1. The radius hypothesis is precisely what supplies the uniform denominator estimate; Haar distribution of the class alone does not supply it. Γ_K≥1 gives X≥(8n²)^n√|Δ_K|. Real quadratic fields, added in the reviewed snapshot, are allowed. Treat the zero sample as failure; conditioning it away changes the normalization slightly.

Concrete fixes: call AF a comparison up to e^{±1/4}, with an exact integral expression for the comparison measure. The actual ideal law weights each ideal by C(Nc). The auxiliary Gamma variable is not itself the actual log-norm deficit: δ̃ uses a cumulative count up to a norm bound, not a uniform shell. In particular, the remark assigning a Gamma(r) upper tail to cofactor norms below Xe^{−V} is not a consequence of AF and is incorrect as a norm-law statement.

One can see the discrepancy even in the uniform-ideal-count model A_all(u)=ρ_K u. Put T=Xe^{−V}. Then δ̃_{S_T}(Xe^{−v}) equals e^{v−V} for v<V and 1 for v≥V. Ignoring the norm-one truncation, whose contribution disappears as X grows, AF's comparison measure assigns

  Pr[Nc≤T] ≈ e^{−V} V^r/r! + Pr[Gamma(r,1)≥V]
             = e^{−V} Σ_{j=0}^{r} V^j/j!
             = Pr[Gamma(r+1,1)≥V].

Equivalently, integrating the single-ideal weight C(u) against ρ_K du gives log(X/Nc) density e^{−v}v^r/r!, the Gamma(r+1,1) density. Thus the displayed expectation really has shape r, while the actual norm deficit in the smooth ideal-count model has shape r+1. The extra exponential component comes from choosing a norm below the cumulative cutoff. For r=1 the difference is already e^{−V} versus e^{−V}(1+V). This does not change the leading L-exponent argument in Corollary AF, but it matters for the claimed “cofactor law” and barrier.af_weight interpretation.

The description “exactly a uniformly random ideal of norm ≥√|Δ|” is too literal as well. The mathematical comparison above is already available in [dBPW, equations (34)–(48)](https://arxiv.org/pdf/2512.01588v2), as the revised note acknowledges.

2. MAJOR — Corollary AF is a reasonable conditional heuristic, not a lower-bound theorem for all usable samplers or fields.

Suppose log X=at+o(t), a≥ζ/2, and log y=b t^α(log t)^{1−α}. Then

  u=(a/b)t^{1−α}(log t)^{α−1}(1+o(1)),
  u log u=(a(1−α)/b+o(1))t^{1−α}(log t)^α.

Under a uniform Dickman-type smooth-ideal model, moving the logarithmic norm by v increases smooth density by approximately exp(v log u/log y). Since log u/log y=o(1), that tilt does not compensate the Gamma factor e^{−v}. With r≤n=t^{1/3+o(1)}, a shift O(n log n) is o(t). The model therefore preserves an N^{a+o(1)} norm scale and gives at least the displayed ζ(1−α)/(2b) smoothness exponent at the smallest allowed scale. Balancing a factor-base exponent α with this search exponent yields α≥1/2. The exponent arithmetic is right.

Three restrictions must stay visible. First, a lower bound log|Δ_K|≥ζt is an assumption, not a consequence of the height upper bound log||f||=O(Q). A field discriminant can be much smaller than a polynomial discriminant because of the index. Second, dBPW's proved lower bounds for smooth-ideal density do not prove the smooth-ideal upper bound required for a computational barrier. Third, “near-uniform input class” must also retain uniform box sampling and the radius hypothesis; class marginal uniformity does not restrict an arbitrary alternative output rule.

For an explicit warning against universal field claims, x^d−2 is Eisenstein, and if N divides m^d−2 then α↦m defines the required residue map. The field has discriminant dividing that polynomial's discriminant, of size at most d^d2^{d−1}, independently of m. This is not a general factoring algorithm, but it disproves “every field we can use for factoring has log discriminant comparable to log N.” The note itself admits that the needed lower bound was not proved.

The ray-condition construction is algebraically valid. For a surjective ψ:R→Z/N, its kernel n has index N. Take a unit residue R0=r² and a trivial input ideal, and impose β≡R0 mod n. Then ψ(β)=r² is a legitimate one-sided row. The general modulus statement imposes the appropriate coprimality/class hypotheses, and for an order one must include its conductor as in the order variant. The certified scales are at least N√|Δ| for the Haar box hypothesis and N|Δ|^{3/2}, times the explicit degree factors, for Algorithm 2. The latter output also contains its known walk factors. These are sampling scales and upper norm envelopes, not pointwise lower bounds on every output norm. Keep the claim that this direct implementation has the Dixon-scale obstruction under the smoothness model; delete an unconditional “L[1/2] at best” claim about all possible designs.

The ±ε transfer concerns the exposed cofactor c after removing the recorded walk-prime product. For an arbitrary family S, the larger event that the total output belongs to S·S_B need not be equivalent to c∈S.

3. MAJOR — SI is sound, but FB does not prohibit making algebraic-child labels exogenous by resampling a constrained polynomial.

SI(i) is the elementary support argument claimed. Deleting at most W values removes at most W from every cumulative count; since u≥x/e^n, the stated local-density error e^nW/(ρ_K x) follows. SI(ii) follows by integrating a pointwise minorization against an A-independent family. Neither statement says that the actual recursively defined good family is adversarial, or that every structured field-dependent supply statement is impossible.

FB's quotient-algebra proof is correct provided γ is integral at q, so its residue class exists. If 1 and γ̄ are linearly independent over F_q, ā−b̄γ̄=0 forces ā=b̄=0. Thus a primitive degree-one element cannot lie in the whole extended ideal under those hypotheses. Forgetting an ideal and retaining its rational norm likewise loses information about conjugate factors. This does not rule out other elements, higher-dimensional searches, explicit transfer certificates, or changes of the field in which retained certificates are interpreted.

Here is a concrete omitted design. Retain a parent point P0=(a0,b0), its integer values, and G. Put L0(x,y)=b0x−a0y and sample a new integral binary form H of degree d−2. Define

  F_H = F + G L0 H.

Then F_H(P0)=F(P0) and F_H(g0,g1)=F(g0,g1)≡0 mod N. Thus the parent's residue identity and norm factorization can be reinterpreted in the new field. For every good child prime p of F(P0), the label (p,a0/b0 mod p) is fixed before H is drawn and remains a root label of F_H mod p. Restrict H, if desired, to preserve the leading coefficient too. Reject reducible or locally bad new polynomials. The remaining coefficients supply fresh randomness conditional on the parent certificate. For finitely many points, use the product of their vanishing linear forms instead of L0. A bounded number of factors of height e^{O(Q)} retains the height scale e^{O(Q)}.

This makes algebraic-child labels exogenous relative to unused randomness in a constrained new field law. It does not transport the old abstract ideal through a field homomorphism, which is what FB discusses. It also does not complete an L[1/3] design: retaining all nodes of a wide tree consumes degrees and coefficient height, and a supply estimate for the constrained law is still missing. It is nevertheless a valid route that the blanket conclusion “algebraic children cannot be made exogenous” misses.

R's multiple-field algebra is correct, with a separate algebraic block, parity bookkeeping, and square root for each field. Independence of a new field makes a rational target fixed conditional on history. It does not by itself prove SI(ii)'s value-distribution minorization or a simultaneous supply theorem for all those targets. State these as outstanding hypotheses, not automatic consequences of freshness. Also account for the number of field blocks before using the construction at scale.

4. MAJOR — R2 has the right AM–GM and exponent calculation, but its lower bounds are model assumptions.

Let l=log ℓ=o(t). The exact minimization relevant to the display is

  t/(d+1)+(d/2)l
   = t/u+(u/2)l−l/2 ≥ √(2tl)−l/2,  u=d+1.

Hence the stated (√2−o(1))√(tl) is right. If l has L-exponent 1−α, this norm logarithm has scale t^{1−α/2}(log t)^{α/2}. Dividing by log B and multiplying by log u_smooth gives a search exponent of scale

  t^{1−3α/2}(log t)^{3α/2}.

Requiring this to fit t^α(log t)^{1−α} gives 1−3α/2≤α, equivalently α≥2/5. At equality the logarithmic powers agree as well. There is no AM–GM or exponent-arithmetic mistake.

However, a polynomial having a root modulo N does not imply ||f||≥N^{1/(d+1)}. Already f=x³−2 has height 2 and the root 5 modulo N=123, while 123^{1/4}>3.3. More generally take divisors of m^d−2. Nor is |F(P)|≥||f||V^d a pointwise inequality: the familiar elementary inequality has the opposite direction, up to d+1. Root neighborhoods, coefficient shape, skew, and target-conditioned boxes matter. The claimed uniform sublevel-set estimate is not proved and would need growing-degree constants; it cannot be inserted as a theorem by calling it an area argument.

Concrete fix: formulate a conditional proposition with a balanced coefficient/root-height model, a typical norm-size assumption, a random-integer smoothness model, and the split scale as explicit premises. Then “2/5 floor in this model” is correct. “Every rational-only or fresh-field design has this floor” is not established. Designs retaining large algebraic columns for recombination, absorbing allowed bad primes in H, changing search dimension, or selecting exceptional small values are outside that formulation. No faster rigorous algorithm follows merely from listing them.

5. MAJOR — OD's children lower bound is wrong; the resulting mandatory-overdetermination claims do not follow.

For fixed g,w and fixed points, the displayed dependence in c is correct: an integer relation among s>d rational vectors v(P_i) gives the required identity. A cleaner full-law observation is that a binary form of degree d has d+1 coefficients. At fixed P_i, its evaluation vector belongs to the image of an s×(d+1) matrix, regardless of how g,w,c encode those coefficients. Thus for s>d+1 there is already a fixed linear dependence. This avoids treating a discrete digit-construction map as an algebraic parametrization. For adaptive points, these statements must be applied after specifying what is fixed or conditioned on.

The inference from log|F(P)|≥(d/2)z to at least d/(2λ) algebraic children is invalid. It counts all norm mass as if it came from distinct primes above B. Factor-base primes may carry all that mass, and repeated large prime factors create one target, not many independent target values. For an algebraic target one must also remove its prescribed q factor before counting cofactor mass.

Precisely, if C is the algebraic cofactor and C_large is its part above B, then the cap p≤e^{λz} implies only

  Ω(C_large) ≥ log(C_large)/(λz),

where Ω counts prime factors with multiplicity. The note has no lower bound for log(C_large), and this inequality is not a lower bound for the number of distinct children. Memoization can further identify children shared between different branches.

Consequently, “one node below 2/5 already has more children than parameters,” “every depth-3 configuration,” and the corresponding uses of barrier.min_children are not proved; literally stated, they are wrong. A thin or terminating tree is not excluded. Under an extra typical-factorization model with order-d distinct children, the numerical comparisons m≈d/(2λ), 2m>d for λ<1, and m²≫d are fine. Label them as such and remove “every.”

Affine dependence prevents a proof that simply multiplies independent-coordinate probabilities. It does not prevent averaged correlation theorems, conditioned moments, or other field-randomization arguments. After enough values determine c, further events are deterministic given c; one can still study the measure of the set of c satisfying those events. The subsection's explicit “method barrier, not impossibility” sentence is appropriate. The Verdict and Candidate card do not preserve that limitation.

6. MAJOR — BT's L-scale arithmetic is right; neither its necessity claim nor the union-bound exclusion is proved.

Accept the model recurrence T_j≈m p_{j+1}^{−m}T_{j+1}: each restart tests an independent batch, all m children must lie in the next good set, and retained successes do not change the stipulated law. For h levels,

  log overhead = h log m + Σ_j m log(1/p_j).

With m≈D′, h≤A log t, and p_j≥t^{−C}, this is O(D′(log t)²)=O(E log t). Thus the guaranteed model bound is exp(O(E log log N))=L[1/3]^{O(log log N)}, with the displayed (log log N)^{5/3} power. This is generally not a fixed-constant L[1/3,c] bound. If h≈log t and p_j=t^{−C} throughout, the order is attained within this model. A lower bound p_j≥t^{−C} alone gives O, not asymptotic equality: p_j=1 is an immediate counterexample to the printed “≈” conclusion.

The exact fixed-constant condition in this model is Σ_j m_j log(1/p_j)=O(E), plus the branching terms. With comparable m_j and h≈log t, it means bounded average log(1/p_j), i.e. a geometric-mean density bounded below. It does not require every level's density to have a separate positive constant bound. Marginal densities without the independence/restart assumptions do not determine a tree cost: fully correlated children can cost p_j^{−1} instead of p_j^{−m}, and structural selection or reuse can change it further.

The union-bound paragraph also overreaches. There are at most (d+1)e^{b1 Q} candidate targets, so an exponentially small per-target dead probability would suffice for a direct union bound. But Θ(1) relative variance of a raw count does not imply a Θ(1) dead probability, or rule out useful conditional moments. Take Z_A=M c_A with c_A equal to 1 or 2: its relative variance is bounded away from zero, yet Pr[Z_A=0]=0. Field-dependent local factors can be conditioned on or normalized. E1's finite-size rate variation is not an asymptotic variance theorem for the number of good box points. Replace “averaging cannot certify AD by moments” with the narrower statement that the particular unconditioned Chebyshev estimate described here is insufficient. The asserted necessity of a specific e^{−Ω(kE)} relative error also needs an actual moment argument.

7. MAJOR — X_ERH works after constant and size-language repairs; the unconditional factorial correction is right.

Let B0=log Δ̄. The quoted GRH Chebotarev error is the one in [Lagarias–Odlyzko, Theorem 1.1](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf). For L_h/K, its group has order 2, so the interval main term is (1/4+o(1))x/log x, with error O(√x(B0+2d log x)+B0). It is enough that √x dominate (B0+2d log x)log x. The proposed polynomial analytic threshold is more than sufficient for large parameters, with a sufficiently large absolute C0, and must also be applied at x/2.

Prime ideals of residue degree at least 2 number at most dπ(√x)≤d√x. Chebotarev already omits primes ramified in L_h/K; subtracting them again is harmless. Rejecting rational primes dividing D can remove up to d·ω(D) prime ideals, not just ω(D) ideals. All these losses are negligible here. For a rigorous 1/(5d), retain, for example, at least 0.23x/log x first-degree inert primes and divide by at most 1.1d x/log x. The written proof first weakens the numerator to x/(5 log x), then subtracts losses and divides by 1.1d; those displayed bounds only give approximately 1/(5.5d). The desired constant follows from the unused margin in the true main term, not from that chain of inequalities.

The sampler which draws integers p in the interval, tests primality, finds roots, accepts with probability (#roots)/d and picks a uniform root is uniform over surviving degree-one primes. The lower count gives O(d log x) expected integer trials. Uniform finite-field algorithms have polynomial cost in d and log x. Under GRH no exceptional near-one zero needs a class exception. For each h≠1 the escape probability ≥1/(5d) gives survival at most exp(−tχ/(5d)); a union bound over H proves the claimed kernel probability. State GRH for ζ_K as well if it is invoked; when H is nontrivial it also follows from GRH for a quadratic L_h through the factorization into ζ_K and an entire quadratic Hecke L-function. The trivial-H case needs no characters.

The atom-avoidance argument is correct for accepted tree atoms. Every prime of F(P) outside D is either its target prime or an allowed cofactor prime, hence at most e^{b1 Q}; factors from the leading coefficient are in D. Since x/2>e^{b1 Q} and p∤D, a character prime divides no atom norm. Failed scan points are not row atoms.

However, this forces x≥e^{b1 Q+1}: these primes have L[2/3]-size, not polylogarithmic size. Their bit length is log x=O(Q+log t)=poly(t), which is all the complexity proof needs. Correct the Verdict and the sentence after X_ERH accordingly. The analytic Chebotarev threshold alone is polylogarithmic. Truly small primes would require a different treatment of nonunit atomic evaluations, or additional rejection accounting.

Finally, Stirling gives, for d=δD′,

  log((2d)!) = (2δ/3+o(1))E.

Thus a bit length (2d)!·poly(t), including any extra polynomial logarithmic factors needed to suppress the zero term, is L[1/3,2δ/3+o(1)]. Polynomial arithmetic in that bit length incurs L[1/3,Θ(δ)] overhead, not e^{o(E)}. The Θ(δ) statement correctly leaves the arithmetic exponent unspecified. At the modeled degrees for α>1/3, d log d=o(Λ_α), so that overhead is negligible. The reviewed snapshot appropriately acknowledges the related observation in [Lee–Venkatesan, Remark 6.7](https://arxiv.org/pdf/1805.08873); its particular 4δ/3 constant is not an automatically transferred constant for this implementation.

8. MAJOR — NR omits quantitative sampling and admissible-field hypotheses.

“Any Las Vegas-samplable law” is not a running-time bound. A sampler may have enormous expected cost even when its output description is short. Require a single uniform sampler for the family (F_N)_N, with expected e^{o(E)} sampling/certification cost, or include an explicit sampling exponent in the theorem. Fixed constants and threshold computations should be computable. LV6′ as currently worded does not specify the cost of this arbitrary field law.

AD is measured under the stated field law before NEST's field gcd rejections. If a sampled field has gcd(disc f,N)=N, the algorithm rejects it. AD's statement that its targets have good points does not say it passes that check. An irreducibility certificate does not fix this. The claimed probability ≥1/2 of a usable good field therefore does not follow as written. This is especially relevant because the law is allowed to be arbitrary, rather than the particular product law for which rejection estimates might be proved.

Concrete fix: require the law to output admissible fields with gcd(g1 f_d,N)=1 and ψ(g′(β)) a unit, together with any additional full-discriminant check the algorithm chooses; or allow a factor exit and impose a positive uniform lower bound on “factor exit or admissible AD-good field.” An alternative is to prove suitable rejection estimates for one specified raw field law and formulate AD for its accepted law. Requiring full disc(f) coprimality is stronger than the local derivative condition LV6 actually needs.

With that correction, the restart argument works. On each new field trial there is probability at least 1/2 of a good admissible field. It never fails in DESC. Each bad field either causes a restart or produces valid rows. Uniform split acceptance and the worst-case box/tree caps bound work before a restart as well as work in a good field. Fresh-root rounds have success probability at least 1/2, so a field that continues producing rows does not trap the algorithm indefinitely. Summing conditional per-trial cost bounds gives at most a factor 2 for field trials; it is not necessary to recognize good fields.

This is a gap in the present “proved conditional” statement, rather than a need for FM or LD*. State ERH + AD + the explicit uniform LV6 and field-generation assumptions in the plain-language Verdict too.

9. MAJOR — AD is logically definable, but the “standard heuristic” explanation misses irregular lattices.

The quantifier pattern is mostly good: fix the parameters and one uniform field-law family; for every sufficiently large preprocessed N, with probability at least 1/2 over A, every target in T(A) has at least one specified good point. It is a finite, decidable property of each field, even though testing it globally is expensive. It is not circular, since “good” is defined by an explicit factor-size condition, not by algorithmic eventual success.

State the cofactors by target type and take absolute values. The actual allowed bound is y(q)=max(B,q^λ), with the added exclusion p∤D_A for p>B. Calling both cofactors merely q^λ-smooth is inaccurate near the bottom. For a good algebraic cofactor prime p>B, the exclusion p∤D_A gives a simple root; primitivity and p∤f_d ensure p∤b, so (p,a/b) is defined. Large excluded primes cannot be children, and small bad primes are deliberately retained in H. These are sufficient closure conditions for the defined descent.

There is a separate geometric issue. Take a rational target q|g0, with q∤D_A and q>e^{2sE}. Then g1 is invertible mod q, and its lattice is a≡0 mod q. Since V=⌊√q e^{sE}⌋<q, every box point has a=0. Primitivity and b>0 force the unique point (0,1). It is good only if g0/q and f(0) satisfy the actual cofactor restrictions. There are not e^{2sE} independent trials. The analogous algebraic target (q,0), for q|f(0) and q∤D_A, has exactly the same one-point box. If f(0) has two comparable allowed-range primes q1,q2 above B and above e^{2sE}, the cofactor at at least one of these targets violates the q^λ cap; for comparable primes and fixed λ<1, both do.

More generally, a lattice vector P0 sufficiently shorter than q/V forces every box point to be collinear with it: det(P0,P) is a multiple of q and has absolute value <q. There is then at most one primitive point in the half-box. The Poisson estimate using e^{2sE} trials per target is therefore not a uniform heuristic for AD as stated. These irregular targets are naturally endogenous factors of small-point values; they cannot be dismissed by modeling each target lattice as an independent generic lattice.

This does not prove that AD is false for every possible law. It does invalidate the asserted justification that the usual independent-box-point model makes failure doubly exponentially small for every target. Call AD a stronger all-target hypothesis and address irregular targets explicitly. If targets are restricted to regular lattices, prove that descent children remain admissible or provide a priced escape routine. Excluding q|g0 repairs one example, not the general short-vector problem. A weaker visited-target hypothesis could be more appropriate; SI does not rule it out merely because it is field-dependent.

10. MINOR — RW, TS, and T_multi work, but record the following details and correct the exponent notation.

Use T=f_d to avoid overloading t=log N. At a rational target q with current rational exponent e, write G(P)=σqn. The substitution

  q^e ≡ [σ(g1/T)ψ(z_P)/n]^e

introduces z_P^e and n^{−e}. At an algebraic target qfrak already occurring in Z with exponent a, take the next element with exponent e=−a and multiply the rational side by [T G(P)/g1]^{−e}. This preserves the residue identity while cancelling that ideal. By Lemma C, q∤D and primitivity give a unique degree-one prime above q in z_P. Goodness implies q∤F(P)/q, since q>B and q>q^λ, hence v_qfrak(z_P)=v_q(F(P))=1. The cancellation is exact, including signs and all integer multiplicities.

After all expansions the scalar is

  q_row = σ_row (g1/T)^{E_row} ∏_{p≤B}p^{a_p},  E_row=Σ_v e_v.

The reviewed snapshot now states this useful formula. Sign, small-prime, and E_row parity columns make that scalar a rational square on a dependency. Remaining algebraic valuations are supported at the small good primes and the bad primes included in H. No principal-ideal generator, class-group assumption, FM, LD*, or Theorem K is hidden in this induction.

The printed exponent bound |e_v|≤t^{depth} is not literally correct: the lemma has just redefined t=f_d, and root multiplicities can already exceed 1 at depth zero. A correct uniform bound is |e_v|≤(C log N)^{depth+1}, since every relevant integer valuation is O(log N). With depth O(log t) this is exp(O((log t)²)). Negative exponents cause no integrality problem: if Z=A/B, the integral representative A B equals Z B². On a dependency use the accumulated denominator B as an explicit square correction when constructing Y; replacing Z by AB without that correction would change the residue row.

TS's depth and branching upper bounds are valid, independently of the false children lower bound in §5. A target log-norm contracts by λ until it reaches log B; depth is O(log t). Each node has O(t/E) children by counting prime factors of its two values with multiplicity. Hence the entire tree has exp(O((log t)²)) nodes, and the same cap applies before failure in a bad field.

There is a short uniform proof of the claimed primitive box count. Order its M primitive points by angle around the origin. No two lie on the same positive ray. Each consecutive determinant is a nonzero multiple of q, so each intervening triangle has area at least q/2. These triangles have disjoint interiors inside the half-square of area 2V². Thus (M−1)q/2≤2V² and M≤1+4V²/q≤5e^{2sE}. Primitivity is essential; the bound is false for all lattice points in a very skinny box. Specify efficient enumeration as well: Gauss reduction detects the one-line case, where only primitive generators need be inspected; otherwise the second minimum is O(V), and reduced-coordinate enumeration has O(1+V²/q) points up to absolute constants. A naive loop through all b≤V would not meet the cost bound.

T_multi is correct for the current batch before its final gcd is revealed. Conditional on its squares and all root-invariant data, the roots remain independent uniform elements of their μ2 fibres. The selected dependency is nonempty and root-invariant, so its product divided by Y is uniform on μ2(Z/N). For odd N with ω≥2, exactly two of its 2^ω elements give trivial factors, yielding 1−2^{1−ω}. Previous failed gcds may be conditioned on because each new batch has fresh roots. The current X or current gcd outcome must not be included in the conditioning transcript. Canonical recursive factorization outputs and fresh subcall coins preserve this invariance; timing-dependent oracle output order must not be used.

11. MAJOR — Repair the Hensel prime; then the split, linear algebra, size bounds, and c_NR close without circularity.

The inherited SQ algorithm chooses the least unramified prime p∈(B,2B] and assumes Γ is a p-adic unit. That assumption followed from B-smooth atoms in the original LV6 setting. It no longer follows here. Clearing negative exponents replaces Z by AB and introduces square factors from large child primes. Those factors can lie in (B,2B]. Even valuations do not make Γ a unit, and the inverse-square-root Newton iteration is then undefined.

Concrete fix: reuse any character prime from X_ERH as the Hensel prime. It avoids D and every accepted atom norm. This costs only poly(t) bits and gives a unit Γ at an unramified prime. This option exists in LV6's discussion, but it must be selected explicitly in NEST's multi-level adaptation; the default prime is not justified.

LV6′ should also specify the actual uniform cost required. SQ's sign enumeration has a factor 2^{d−1}, not a polynomial in d. This is harmless here: d=O(D′), so 2^d=e^{o(E)}. The coefficient bounds for the monic polynomial of β have logarithm O(dQ)=O(t), and the residual-group bound gives log|H|=O(dt), hence tχ=O(d²t). The sum of absolute atomic multiplicities over a dependency is e^{(b+o(1))E}; its logarithmic square-root coefficient bound is that size times a polynomial in d,t. Product formation and quasi-linear Hensel lifting, multiplied by 2^{O(d)}, therefore cost e^{(b+o(1))E}. An unspecified “polynomial-time square root” in the exponentially long input would not suffice to preserve this constant; require the quasi-linear dependence just described.

For completeness, the split estimate does remove the 2^ω loss. Write h=2^ω, n=|P|, U=the k-element prime products. Same-square-class pairs number at least |U|²/h. Pairs sharing a prime number at most (k²/n)|U|²; n≥4k³h is more than sufficient to leave at least |U|²/(2h) coprime pairs. Distinct coprime pairs have distinct ratios mod N because their cross-product difference has magnitude <N. Gauss reduction returns their primitive vector because its determinant with a shortest vector is a multiple of N of magnitude <2X²<N. Dividing by the square subgroup size φ(N)/h gives probability at least |U|²/(2φ(N)). The lower endpoint of S″'s original window is not needed in this version; state that slight generalization.

Here log|D_A|=O(t). Removing its primes from (Y,2Y] removes O(t/Q) primes. Also log h=O(t/E)=O(Q/log t), whereas log n=b1 Q+o(Q), so the required n≫k³h holds. With k=D′/(2b1)+O(1) and log Y=t/(2k)−log 2+O(1/k),

  log Pr[accept] ≥ −2k log(k log Y)−O(k)
                  = −(1/b1+o(1))E.

The note's intermediate numerical constant 4 does not follow from its coarser n≥Y/(3 log Y) together with the factor 2 in Y; use 6+o(1), or use the stronger valid n≥Y/(2 log Y). This changes nothing at L-scale. This argument is elementary S″ counting, not a hidden appeal to Theorem K.

The factorization inputs in a box satisfy

  log|F(P)| ≤ κQ+d(z/2+sE)+O(log d)
            ≤ (δb1/2)t+O(Q),
  log|G(P)| = O(Q).

Thus every such recursive input is at most N^{τ_F+o(1)} with τ_F=δb1/2<1. Split halves are at most √(2N). This applies to unsuccessful scan points too, not just good ones.

There are e^{(b+o(1))E} columns and rows, each row has e^{o(E)} nonzero entries, and the character-column count is polynomial. Verified Wiedemann/KS linear algebra costs e^{(2b+o(1))E}; square roots and character evaluations fit the smaller e^{(b+o(1))E} budget after the preceding fixes. With admissible-field sampling as in issue 8, the cost recurrence is therefore

  C(N) ≤ e^{(2b+o(1))E}
       + e^{(b+1/b1+o(1))E} C(√(2N))
       + e^{(b+2s+o(1))E} C(N^{τ_F+o(1)}).

Expected costs of adaptive subcalls are bounded by conditioning on their inputs and using fresh coins, then summing by monotone convergence. Dependence of the decision to call on earlier canonical outputs causes no problem. Complete factorization after a proper split adds only a polynomial factor when charged through its divisor tree. The trial-division race supplies termination on the finite small-input range.

For any fixed c larger than

  max{2b, (b+1/b1)/(1−2^{−1/3}), (b+2s)/(1−τ_F^{1/3})},

all recursive terms have strict exponent slack after substituting C(n)≤exp(cΛ_{1/3}(n)). Strong induction absorbs the uniform o(1) terms and polynomial factors. Thus c_NR is correct. There is no factoring circularity: all auxiliary inputs are smaller for sufficiently large N, and the finite remainder is covered by the race. Under the repaired hypotheses, neither Theorem A's smoothness moments nor FM, LD*, K, or an unproved endogenous density estimate reappears. The numerical example gives 24.2366105093. Its δ=1 is outside the printed open interval; δ just below 1, with the other displayed parameters, gives the same maximum from the split term.

12. MINOR — E1's preregistered outcomes are reported honestly; correct the uniform-rate label and interpretation.

I recomputed the summary directly from all 12 JSONL instances, including the 2,000 parent-bootstrap resamples with seed 0. P1's n=14333, both rates, difference, and both confidence endpoints match exactly. P2 and P4's n=1383 and their intervals also match exactly. The child histogram, 18 degenerate children, and 144.6381979-second total match. The M and parent-count changes are disclosed in Addendum A; the report uses its S=32 and does not silently substitute the original parameters.

P3 has 211 all-good parents among 4571 eligible parents, observed rate 0.04616057755 versus pooled expectation 0.02697910155, ratio 1.71097534368. It fails the frozen [0.75,1.33] criterion, and both the JSON and note say so. There is no suppression or relabeling of that failure.

Using each field's matched-uniform type rates gives ratio 1.05127478127, agreeing with the exploratory 1.05. However, the matched-uniform algebraic rates range from 0.1098570354 to 0.5873786408. The printed 0.108–0.612 endpoints are the endogenous rates, 0.1075996990 to 0.6116504854. Fix the label/numbers.

“Within a field, children behave independently” is not established by one aggregate ratio near 1. Independence is a joint-distribution claim; child norms and multiplicities, repeated targets, parent selection, and estimation of the reference rates matter. State that the exploratory adjustment makes the aggregate all-good statistic consistent with the product model. The positive P2/P4 differences and small norm mismatch are descriptive evidence for a matching effect, not a proof that no endogeneity effect exists.

The parent bootstrap follows the preregistration and is correctly reported. It does not provide an uncertainty estimate over the field law: the experiment has 12 field clusters, with many dependent observations within each. Any population-level conclusion about fields needs field-level or hierarchical uncertainty, explicitly exploratory if added now. In particular, E1 cannot prove the asymptotic moment obstruction claimed in §6.

13. MAJOR — Rewrite the Verdict to match the actual scope and acknowledge the omitted designs and geometry.

The defensible result is: this note does not prove an ERH-only L[α] factoring bound for any fixed α<1/2; the direct dBPW substitution has a discriminant-scale cost under the stated additional norm/density model; the independent-coefficient proof strategy lacks the joint supply estimates it would need. It is also useful that the corrected NR mechanism reduces the problem to an explicit, very strong relation-supply hypothesis plus uniform algebraic and field-generation conditions.

The following phrases should be removed or restricted:

  • “every field we can use” has N-power discriminant scale — no such lower bound is proved, and small-field special inputs contradict the universal wording;
  • “below 2/5 a descent must create algebraic children” — only a particular size model is analyzed, with unjustified lower bounds;
  • “algebraic children can only be descended in the same field” — FB excludes one direct transfer, not constrained polynomial resampling or all other certificate designs;
  • “every depth-3 configuration” is overdetermined — the children lower bound is false;
  • “all field-moment methods beyond depth 1” fail — dependence does not imply impossibility of correlation estimates;
  • “equivalently, AD” — AD for every target is a sufficient strengthening, not an equivalent or necessary formulation of visited-target relation supply;
  • “one clean hypothesis” makes everything rigorous — the text also assumes LV6′, lacks a sampling/admissibility condition, and must repair the Hensel implementation;
  • “ERH does nothing for endogenous children” — it does not presently furnish the needed two-sided supply theorem; the adaptive one-sided lemma already shows that a literal “nothing” is false.

The missing short-vector case in issue 9 is a substantive additional design barrier. Conversely, conditional polynomial resampling in issue 3, higher-dimensional searches, conditional/local-factor moments, and global relation recombination are not excluded by the arguments given. None is established here as a working faster factoring algorithm. Distinguish absence of a known applicable theorem from a proof that no theorem or design can exist. The literature-search conclusion should retain “no known general rigorous or ERH-conditional bound found,” rather than becoming a mathematical impossibility statement.

Concrete acceptance conditions: correct the false children inference; narrow R2/FB/OD and the Verdict; make the field law and LV6 cost assumptions explicit; select an atom-avoiding Hensel prime; state the integral-representative denominator correction; and replace the all-target Poisson justification of AD with an account of irregular targets. The experiment can remain, with the small numerical label and inference corrections above.

14. MAJOR — Assessment of the 609-line revision received during review; this supersedes the affected earlier findings.

The following changes are valid fixes, and I withdraw the corresponding objections to the revised wording: the sampling law now has an expected e^{o(E)} cost bound; the Verdict acknowledges LV6′; the character-prime bit length and the Chebotarev constant margin are corrected; FB's concluding field-bound claim is qualified; BT is presented as a naive strategy bound; and SQ's factor 2^d is now explicitly allowed as e^{o(E)}. The AF cumulative-density wording is improved. The literal norm-deficit Gamma-tail sentence remains and still needs issue 1's correction.

The new reduced-basis box is a substantive change, not a cosmetic repair of AD. For a Gauss-reduced integral basis u1,u2 of a determinant-q lattice, |det(u1,u2)|=q and the angle condition gives ||u1||||u2||≤2q/√3. A nonzero integer vector has length at least 1; with ||u1||≤||u2|| this gives ||u2||≤2q/√3 and ||u1||+||u2||≤3q. Hence a coefficient box |i|,|j|≤e^{sE} contains at most (2e^{sE}+1)² points and has ambient coordinates at most 3q e^{sE}. Enumeration by (i,j) has the advertised cost. Filtering primitivity, the zero vector, and G=0 is exact and preserves the cap. The one-point sup-norm examples of issue 9 no longer apply to this new box: e.g. basis (0,1),(q,0) supplies many points (qj,i).

The corresponding worst-case bounds are now

  log|F(P)| ≤ κQ+d(z+sE+O(1)) ≤ δb1 t+O(Q),
  log|G(P)| = O(Q).

Thus τ_F=δb1<1 is the correct condition. The revision makes this change both in the parameters and in the proof, so the same recursion argument works with

  c_NR=max{2b, (b+1/b1)/(1−2^{−1/3}), (b+2s)/(1−(δb1)^{1/3})}.

The former δ=b1=1 example is unavailable with these worst-case boxes. The new δ=1/2,b1=1,b=4,s=3/4 example has τ_F=1/2 and c_NR=5.5/(1−2^{−1/3})=26.6602715602. Its revised heuristic coefficient φ=3+1+0.1875=4.1875 gives 2s=1.5>φ/3=1.3958333. Those arithmetic checks are correct. This is a different and stronger-in-search-region AD statement; E1 still tests its original sup-norm boxes and is not experimental validation of the new reduced-basis AD. The assertion about a particular number of irregular targets per field needs a model or counting argument, not just the one-point example.

The revised Hensel rule now says “the least odd prime p>e^{b1 Q} with p∤D.” This repairs the unit condition, but it does not provide an algorithm with the claimed cost for finding that least prime. Scanning successive integers has no proved poly(t) bound here; even the usual GRH prime-gap bound would allow an L[2/3]-scale scan. It is unnecessary to find or certify the least prime. Reuse a sampled character prime, or sample uniformly in an interval of comparable size with primality and gcd checks. The former already exists in the transcript and gives the required expected polynomial bit cost. This remains a MAJOR implementation gap until that policy is stated.

The accepted-field probability gap in issue 8 also remains: the new sampler-cost bound does not imply that an AD-good raw field survives a full-discriminant gcd rejection. Require the law to output admissible fields, or give the positive “factor exit or admissible AD-good field” mass explicitly. The corrected recurrence is sound after these two fixes, without FM/LD*/K. The global exponent-notation correction in RW, explicit denominator square correction, and E1 interpretation fixes remain as well.

The revised R2 correctly admits that the old coefficient-height lower bound was wrong and substitutes a counting heuristic. This is a useful narrowing, but its new infeasible-window argument still is not a proved all-design floor. If Λ=log B and U≈Λ/log t, its endpoint arithmetic is consistent: A=Λ²U²/t and H=t/U give A/H≈Λ²U³/t², with t-exponent 5α−2; A/Λ has t-exponent 3α−1. However, several substantive steps need additional premises. A rational cofactor has log-size approximately log H_g−z/2 plus box inflation, so the omission of the −z/2 term requires an explicit uniform range argument. The dichotomy 2zU>t/2 by itself does not show that an arbitrary feasible node above the window has log H_g≳t/4; one must combine actual algebraic feasibility and coefficient-size inequalities, retaining constants and the degree restriction. Most directly, u log u≤Λ bounds the norm-to-cutoff ratio u; it does not say that a smooth cofactor has at most U prime factors. The claim “splits into at most U pieces, so has a large child” is the same multiplicity/size confusion identified in issue 5. Highly smooth values can have many small factors. To prove a barrier for a probabilistic model, bound the cost of the exceptionally smooth jumps instead of asserting a deterministic child-count or infinite-branch dichotomy. The new height-counting premise must also remain explicitly heuristic; rare or structured polynomial selection is not excluded by candidate-counting alone.

Finally, adding “ignoring the B-smooth part” to OD does not turn its lower bound into a consequence of the actual descent definitions. It introduces a further model assumption, and multiplicities and distinct-child reuse still matter. Conditional polynomial resampling in issue 3 remains an omitted exogeneity design; it is compatible with the newly qualified FB discussion but not with the Candidate card's broad failure scope. The recommendation remains MAJOR REVISION; NR remains GAP as implemented, with a concrete repair and a correct conditional cost calculation rather than a demonstrated fatal obstruction.

Reproduction artifacts in this scratch directory: check_review.py and review_checks.json. They contain independent arithmetic and raw-data checks, not a numerical validation of the asymptotic theorems. Both reviewed note snapshots and their complete diff are retained to distinguish superseded findings from remaining ones.
