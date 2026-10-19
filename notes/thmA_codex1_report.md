Independent analytic referee report on Theorem A
19 October 2026

My recommendation from the assigned analytic half is ACCEPT WITH FIXES. I find no substantive gap in the FM -> LD* -> Paley–Zygmund chain for the final design. There are incorrect ancillary assertions and several statements that need their quantifiers or hypotheses corrected. The strongest incorrect assertion is the inverse-polylogarithmic bound for gamma in thmA_proof.md §6 U4: it fails for the explicit leading-term choice of S allowed in that document. The weaker bound gamma >= L^{-o(1)}, which is sufficient for the stated running-time exponent, survives.

This is an analytic-half recommendation, not an independent certification of all the number-field, character-sampling, square-root, and linear-algebra machinery in LV6. I read the requested underlying sections for their interfaces, but those algebraic results remain inputs to the whole-theorem recommendation. I did not use the excluded verdicts, reviews, candidate cards, or prior audit reports.

Notation below: v = nu = 1/(4k), Lambda = sqrt(log N log log N), L = exp(Lambda), y = B, x = X*, Xr = X_r, and rho(z) = Psi(z,B)/z. All constants can depend on the fixed design parameters. Equations involving V use bounded dyadic envelopes where necessary. N_l denotes the filtered count, and averages in the moment bounds are over BOTH l and the raw field law.

1. Final-design sizes and reduced residues — CORRECT.

Location: ldbound_results.md §1, Lemma GB_w, “Sizes”, Lemma 1.2_w and (1.1); thmA_proof.md §4 Step 1.

The base expansion gives Fhat(g0,g1)=Nw. Its large constant coefficient is accounted for: dividing its contribution by the random window scale gives

  (Nw/M^d)/(MH) = N^{1-d mu+v+eta0/2-h+o(1)}
                  <= N^{-eta0/2+o(1)}.

Thus the window floor h >= 1-d mu+v+eta0 is precisely what is needed. One cannot use the earlier bound Fhat(P)=O(MV^d), but the final argument does not require it.

For a counted point, t=|G(P)|=l|n| with n B-smooth, b prime >B, b<l, and l not dividing N. Hence (b,t)=1. Since (g0,g1)=1, any prime dividing both g1 and t would divide b, which is impossible. Preprocessing gives (N,t)=1, and w>t gives (w,t)=1. Homogeneity therefore yields

  Fhat_{g,w}(P) = N b^d g1^{-d} w (mod t).

The multiplier c_P=N b^d g1^{-d} is indeed a unit independent of w. The balanced representatives used in constructing Fhat may vary nonlinearly with w; this does not affect this congruence. The rational smoothness restriction is essential here.

The two smooth-size exponents are

  xi_r = mu-v > 0,      xi_F = mu+h+d v,
  theta = xi_r+xi_F = 2mu+h+(d-1)v.

2. Published smooth-number inputs and their actual ranges — CORRECT WITH FIXES.

Location: pcest_results.md §2.1, especially (H), (LV8.5), (T), (Hi), (C8.16), and (LV5).

I checked the PDF numbering of Lee–Venkatesan: Lemma 8.5 is on printed p. 41, Fact 8.11 on p. 43. The HTML rendering renumbers them, so its displayed numbering should not be used to identify the citations. Lemma 8.5 states a maximum-residue estimate; its proof, equation (8.2) and the subsequent estimates, bounds the stronger absolute-character majorant W needed here. The proof's absorption of logarithms is valid in the present regime. [Lee–Venkatesan, §8](https://arxiv.org/pdf/1805.08873).

Harper's Proposition 2 does have the two sums over primitive conductors and inducing moduli used in pcest. Its hypotheses include 0<eta<=1/80, y<=x^{9/10}, and x^eta<=Q<=sqrt(x). No restriction that the inducing modulus itself be smooth is imposed. These hypotheses hold here with fixed power margins. [Harper, Proposition 2, pp. 19–20](https://arxiv.org/pdf/1208.5992).

There is a genuine transcription error in (T). Tenenbaum's condition is

  omega(t) <= y^{c/log(1+u)},

not y^c/log(1+u). This is a much smaller allowed range in general. It still contains every modulus used here: log omega(t)=O(log log N), whereas log y/log(1+u) is of order sqrt(log N/log log N). The relative error in the displayed main term is O((log log N)^2/sqrt(log N log log N))=o(1). For t=l s, one can apply the smooth-modulus theorem to s directly; excluding the prime l changes the totient factor by 1+O(1/l). [Tenenbaum, main theorem, condition (Omega), p. 2](https://tenenb.perso.math.cnrs.fr/PPP/Psi_q(x,y).pdf).

The Hildebrand short-interval range is ample for intervals of fixed relative length, and the bounded-factor and B^{O(1)} changes of scale used below cost at most powers of log N. The latter statement can also be recovered by applying the saddle-point ratio estimate a fixed number of times. It is not permissible to replace an arbitrary N^{o(1)} change of scale by a polylogarithmic loss. [Hildebrand–Tenenbaum, Corollary 1.7 and Theorem 5.1](https://www.numdam.org/item/JTNB_1993__5_2_411_0.pdf).

3. Restricted BV, Lemma 8.5' — CORRECT in the final regime.

Location: pcest_results.md §2.2; ldbound_results.md §2, FM Step 3.

Put q0=mu+v and xi=mu+h+d v. Uniformly over the whole grid,

  Q=N^{q0+o(1)},    x=N^{xi+o(1)},
  Q/l=N^{mu-v+o(1)},
  xi/2-q0 = [h-mu+(d-2)v]/2 >= eta0/2.

Consequently Q/l >= x^epsilon for some fixed epsilon>0, and Q <= sqrt(Psi(x,B)) with a fixed polynomial saving. The same statements hold at the lowest grid point x=X* omega^{-2}, because omega=L^{C1} changes log x only by O(Lambda).

The decomposition into t=l s is unique: s is B-smooth and l>B. Therefore

  #T_Q = sum_l [Psi(Q/l)-Psi(Q/(2l))]
        asymp Q rho(Q/Y)/log Y.

For characters principal at l, smooth summation variables are automatically coprime to l. This gives exactly the W(Q/l)/(l-1) contribution in the note. Summing the main terms produces rho(Q/Y)/log Y; the large-sieve error is even smaller by a factor of order Y.

For characters nonprincipal at l, the primitive conductor contains l. Choose a fixed eta smaller than both 1/80 and v/xi. Then l>=x^eta and Q>=x^eta. Harper applies by positivity to the enlarged family. Dividing its error by Psi(x,B) delta_T gives terms bounded by

  N^{q0-xi/2+o(1)}  and  N^{-xi eta+o(1)},

because rho(x) and delta_T are N^{-o(1)}. Both are polynomially small. Here the use of N^{o(1)} is legitimate: it is absorbed by a fixed negative power of N, not silently charged as L^{o(1)}.

The small-conductor/Page-character part inherited from LV needs logarithmic factors to be absorbable into exp(-c u/log^2 u). In this application

  u/log^2 u asymp sqrt(log N)/(log log N)^{5/2} >> log log N.

Thus the potential endpoint issue in the broader published range is absent. I found no failed external hypothesis at this interface.

4. FM fibres, induced law, and coprimality — CORRECT.

Location: pcest_results.md §2.3 Step 1 and Step 2(a), transferred in ldbound_results.md §2.

For fixed primitive P=(a,b), a fibre G_g(P)=l n is an integer line with primitive direction (a,b). Its intersection with the g-square contains O(M/V) points. The continuous convolution supplies a one-sign interval of n-values of length asymp MV/l, bounded away from zero by a fixed multiple of MV/l, on which the unprimitive fibre length is at least c M/V.

For a divisor e of l n, the condition e|gcd(g0,g1) selects exactly one residue class of the line parameter modulo e. This follows by choosing a Bezout combination of a,b equal to 1; the determinant condition ensures compatibility. Möbius inversion therefore gives

  #G_P(n) >= c(M/V) phi(l n)/|l n| - O(2^{omega(l n)}).

Since mu-v>=eta0, M/V=N^{mu-v}L^{O(1)} dominates the uniform divisor error N^{O(1/log log N)}. The remaining lower bound is c M/(V log log N).

This is the correct induced law. There is no assumption that t is uniform among the special moduli. Its actual fibre weights are controlled above and below. Hildebrand applies to the one-sign interval after reflection if needed. The total rational mass per P is at least rho_r/(l (log N)^{O(1)}).

The O(1) primes l in J that divide N can be omitted uniformly: their proportion is N^{-2v+o(1)}, and even the cap O(S^2) makes their contribution negligible compared with S^2 rho_r rho_F.

5. FM kernel, flattening, signs, grid, and bad-modulus weights — CORRECT WITH FIXES.

Location: pcest_results.md §2.3 Steps 2(b)–4; ldbound_results.md §2, FM re-check.

The coefficient-to-value smoothing can be reconstructed directly. For l_i=|a|^{d-1-i}b^i, add an independent uniform remainder of length l_i using a lift supported in a cube of radius O(V). Then l_i times a uniform coefficient plus this remainder is uniform on a consecutive interval of length (2H+1)l_i. Combining the d lifts perturbs the cube law in total variation by O_d(V/H). Swapping a,b when necessary and changing coefficient signs are harmless.

The resulting convolution has maximum mass O(1/(H V'^{d-1})) and centre displacement O_d(V^{d-1}). Flattening inside Y_omega=(2H+1)V'^{d-1}/omega changes its total mass by O(1/omega); it does not discard that central region. Since H/(omega V) tends to infinity polynomially, the two tails are monotone away from zero and admit the claimed positive mixture of one-sided uniform intervals.

For every such uniform interval,

  |Fhat(P)|/(t Y_omega) << omega N^{-eta0/2}.

The sign-crossing correction is adequate. One may keep a same-sign portion covering absolute values between |Fhat(P)| and an endpoint of size t y-|Fhat(P)|. The loss is O(|Fhat(P)|/t+1) terms. Choosing a grid point in (t y/4,t y/2] and decreasing the harmless numerical constant if necessary gives smooth density at least rho_F/(log N)^C. Both the translation error and V/H are polynomially smaller than rho_F; choose C1>xi_F/(2a) with a fixed margin to absorb 1/omega as well.

The grid is independent of P,g,w,t and the residue. A modulus is declared good for all reduced residues and all grid points, so conditioning on its moving residue costs nothing. Lemma 8.5' bounds the number of bad moduli by

  O((log omega) epsilon_N) #T_Q.

Multiplying by the maximal fibre weight O(M/V), and then comparing with the good fibre lower bound M/(V log log N), loses only a log log N factor. This still tends to zero. The nonuniform P-dependence is harmless because candidate counts across J differ by bounded factors and the fibre bound is uniform in P.

Two formulations should be corrected. First, the final FM conclusion is averaged over l; the notation E_A N_l in thmA_proof.md §4 Step 3 must not suggest a lower bound for each individual l. The next step needs only the average. Second, the dyadic cover of [MV/32,2MV] is O(1) after using V_- and V_max. If one instead literally enlarges it by an unspecified L^{o(1)} factor, an O(1) cover is not justified. Using the actual envelopes fixes the wording.

The result consumed by PZ is therefore valid:

  E_{l,A} N_l >= S^2 rho_r rho_F/(log N)^C.

6. Lemma W and Corollary W — CORRECT.

Location: ldbound_results.md §3; their use in §5 Step 4.

The interval Brun–Titchmarsh bound has denominator log(W/t), and requires t<W. It gives

  Pr_w[w=a (mod t)] <= 4 log W/[phi(t) log(W/t)]
                      <= 24/[eta0 phi(t)].

Indeed log(W/(2MV_max))=(eta0/2)log N-O(Lambda), so the required lower bound (eta0/3)log N holds uniformly. In the theorem's range theta<1, log W<log N for large N; the displayed looser bound 2log N is sufficient. Nonreduced classes have probability zero because every sampled prime w exceeds t. The case t=1 is trivial.

I could not retrieve the original 1973 scan for a page-by-page check. The precise interval formula and condition t<W are also recorded in Motohashi's research paper, equation (5.4), and match the invocation here. This is a source-access limitation, not an unidentified extra hypothesis. [Motohashi, equation (5.4)](https://www.mathnet.ru/php/getFT.phtml?jrnid=spm&paperid=35&what=fullt).

Since c_P is a fixed unit, multiplication by it permutes the reduced classes. Corollary W follows for every fixed nonnegative Phi. Its use would be invalid if Phi still depended on the offset Fhat_{g,w}; Step 3 of LD* removes that dependence before applying it.

7. Lemma DS — CORRECT in the ranges used.

Location: pcest_results.md §3, Lemma DS; ldbound_results.md §4(A) and §5 Step 5.

For e=(r,q), remove e, and put q'=q/e. The remaining progression is reduced modulo q'. Set D=min(A/q,X/(e B^{C0})). Every relevant smooth integer n/e has a smooth divisor D/B<d<=D. Since d is coprime to q', its complementary factor occupies one progression, with at most 2A/(dq) possibilities. Summing over possible d is an upper bound even though a number can have several such divisors. Dyadic summation gives

  count <= C(A/q) log(2B) rho(D/B).

As X/(D/B)<=B^{2C0+1}, density comparison costs only (log X)^{O(C0)}. This proves the claimed upper bound, including nonreduced residues. In the applications D/B=N^{mu-v+o(1)} or N^{xi_F+o(1)}, so issues with a very small divisor cutoff do not arise.

8. Lemma Kern — CORRECT.

Location: pcest_results.md §4.1, parts (i)–(iii); ldbound_results.md §5 Step 1.

A unimodular change sends primitive P to (1,0), and sends P' to (u,Delta) with (u,Delta)=1. The gcd of the Veronese 2-by-2 minors is then |Delta|. This proves both the image index and the asserted congruence description; primitivity of P' also proves that kappa is a unit modulo Delta.

The kernel consists exactly of integral multiples of (b x-a y)(b' x-a' y). Gauss's lemma supplies integrality, and its monomial translates give a basis with coefficients O(V^2). The fundamental-cell argument is valid because H/V^2=N^{h-2v}L^{O(1)} tends to infinity polynomially. Its covolume is ||v(P) wedge v(P')||/|Delta|. The adjacent bottom minor is exactly Delta (bb')^{d-2}; hence

  sup_z Pr_c[Phi(c)=z] << 1/(H^2 V^{2d-4}).

No unspecified compactness constant is needed for the actual candidate set.

9. Additive Parseval and the partition in L2_fam — CORRECT.

Location: ldbound_results.md §4, “Additive Parseval”, rearrangement, and partition (A)/(B)/(C).

Orthogonality gives Q(m)=m^{-1} sum_{h mod m}|S(h/m)|^2. Reducing h/m gives Q(m)=m^{-1}sum_{r|m}f(r), with f(1)=#S*^2. This includes signed smooth integers and the diagonal exactly.

The common X* is essential: it makes f(r) independent of l. Since B and 2S^2 are both smaller than Y, a contributing r either contains no J-prime, or contains exactly the sampled l to exponent 1 or 2. The three classes are disjoint and exhaustive for the rearranged sum. A denominator with two J-primes has zero multiplicity.

The needed large-sieve estimate is the additive one for a sequence supported on an interval of length O(X*), with zero coefficients at nonsmooth integers. It is valid for negative indices by translation; it does not require a multiplicative sequence. Even an absolute-constant version of (K+Q^2) suffices. [Vaughan, large-sieve notes, Theorem 3 and proof](https://personal.science.psu.edu/rcv4/LargeSieve.pdf).

10. L2_fam, class (A) — CORRECT.

Location: ldbound_results.md §4(A).

DS bounds each residue class modulo r<=R0 in the outer interval. The central interval contributes at most O(X*/(R0 r)+1), which is absorbed because R0^{-1}<=rho_F and R0<<Psi(X*). Consequently f(r)<<Psi(X*)^2(log N)^C.

Writing r2=(r,q), r1=r/r2, the divisibility condition r|nq is equivalent to r1|n. If r1 is not B-smooth there are no n. Otherwise its multiplicity is at most

  C Psi(X'_+)/r1 times a power of log N.

The assignment r -> (r1,r2) is injective and r2|q. Harmonic summation therefore contributes tau(q)(1+log R0), not an uncontrolled divisor factor depending on n.

11. L2_fam, classes (B), (C), and every margin — CORRECT.

Location: ldbound_results.md §4(B), §4(C), and deduction of (4.1).

Write X=X'_+, Z=X*, r_r=rho(X), r_F=rho(Z). In (B), a dyadic block has multiplicity at most C Xq/R. Applying the large sieve and summing geometric series gives

  (B) << Xq Z Psi(Z)/R0 + X^2 q^2 Psi(Z).

Relative to Psi(X)Psi(Z)^2, these are

  O(q/(R0 r_r r_F)),
  O(X q^2/(Z r_r r_F))
      = N^{-h-(d+1)v} L^{O(1)}.

The first is absorbed by the prescribed R0; the second has a fixed negative power.

In (C), only one l contributes for any fixed denominator. Thus the two useful multiplicity bounds are C Psi(X)/pi_J and C XqY^2/(pi_J R). Applying them separately to the X* and R^2 terms of the same large-sieve bound is legitimate: if nu<=A and nu<=B, then nu(K+R^2)<=AK+BR^2. The resulting terms are

  O((log N) Psi(X) Z Psi(Z)/pi_J),
  O(X^2 q^2 Y^4 Psi(Z)/pi_J).

Their relative sizes are respectively

  N^{-2v} L^{O(1)},
  N^{(5-d)v-h} L^{O(1)}.

The last margin is strictly positive in the required direction:

  h-(5-d)v >= (d-3)v+eta0 > 0.

This includes d=3. Averaging over l supplies the indispensable factor 1/pi_J. Nothing here establishes the corresponding assertion for a fixed l, and no such assertion is needed.

For filtered n, m=l^2|n|q=l q t and phi(t)>>t/log log N, while t>=MV_-/32. Hence

  m phi(t) >> Y q (MV_-)^2/log log N.

This proves exactly the denominator in (4.1), including its extra factor Y. There is no lost power of l.

12. LD*, Steps 0–4 — CORRECT.

Location: ldbound_results.md §5, Steps 0–4.

Step 0 uses the raw product law; pairs with l not dividing Delta cannot occur because primitive g is nonzero modulo l. Distinct primitive candidates with b,b'>0 have nonzero determinant.

Step 1 applies Kern with the correct wedge lower bound and counts possible (z,z') in its image lattice. Step 2 is the pointwise inequality

  sum_j A_j(P) A_{kappa^{d-1}j}(P')
       <= (sum_j A_j(P)^2+sum_j A_j(P')^2)/2,

because the multiplier permutes the residue classes. No statistical independence is asserted.

For Step 3, F=Fhat+Gz is injective and z=z' (mod |Delta|) is equivalent to F=F' (mod t|Delta|), regardless of the sign of G. Enlarging to all signed smooth values in the common interval removes every dependence of Q_0 on w except its argument r0. The scale ratio X*/(tT) is bounded under the filter.

Step 4 may therefore average over w using Corollary W. Summing over reduced classes counts a subset of all pairs congruent modulo t|Delta|, so the result is at most C_W Q(t|Delta|)/phi(t). Symmetry deals with the P' term. This is why a joint equidistribution law for the two residues under their shared w is unnecessary.

13. LD*, Steps 5–7, including the power count — CORRECT.

Location: ldbound_results.md §5, Steps 5–7.

After averaging w, the algebraic factor depends on g only through signed n=G(P)/l. On a fixed fibre g=g*+jP, the partner's n' advances in steps of |Delta|/l=q through an interval of length O(q M/V). The filter keeps |n'| between fixed multiples of Xr. DS is applicable after enlarging a short interval if necessary, because

  Xr/(R0 q) / (M/V) asymp S^2/(R0 q) <= 1.

It supplies at most (M/V_-)rho_r(log N)^C partners on a fibre. Dropping the primitivity restriction on g only increases this upper bound.

For fixed P,q there are at most 12 geometric partners P': each of the two determinant equations defines an integer line with spacing |P|>V/2 inside a segment of length at most 2 sqrt(2)V. This removes the geometric partner sum without a factor V or l.

Substituting (4.1), Psi(X'_+)<< (MV/Y)rho_r times polylogarithms, and Psi(X*)<<MHV^d rho_F times polylogarithms, the non-density factors are

  [1/(H^2 V^{2d-4})] [1/M^2] [V^2] [M/V]
       [MV/Y] [(MHV^d)^2] [1/(Y(MV)^2)]
    = V^4/Y^2 asymp S^4.

All M and H powers cancel. The remaining arithmetic sum is sum_{q<=2S^2}tau(q)/q=O((1+log S)^2), which is a power of log N. Thus

  E[N_l(N_l-1)] << (log N)^C S^4 rho_r^2 rho_F^2.

14. Computable S and the extension of LD* — CORRECT WITH FIXES.

Location: ldbound_results.md §5, theorem hypotheses and last “PC” line; thmA_proof.md §1, §4 Step 5, and §6 U5.

LD* is formally stated for the implicit choice S^2=(log N)^{C3}/(rho_r rho_F), whereas the algorithm permits an explicit computable leading-term choice. Its proof does not use that equality in Steps 0–7. It needs only S=L^{sigma+o(1)} in a fixed compact positive range and the consequent fixed choice of C0. The equality is used at the end only to assert E N_l>=1. Amend the theorem statement to make the uniform family of S explicit, and remove that final assertion for general S.

Put m=S^2 rho_r rho_F and z=E N_l. FM and LD* give, for a fixed polylogarithmic P,

  z>=m/P,     E[N_l(N_l-1)]<=P m^2.

Paley–Zygmund then gives

  gamma >= z^2/(P m^2+z) >= 1/(P^3+P/m).

Since the proposed computable S makes m=L^{o(1)}, this proves gamma>=L^{-o(1)} uniformly without assuming z>=1, without computing Psi, and without concealing a fixed L-power. Also

  S^2/rho_r^{-2}=L^{[h+(d+1)v]/(2a)+o(1)},

so the extra lower bound on S mentioned in §1 has room, although Theorem R is unnecessary.

15. The inverse-polylogarithmic assertion in U4 — WRONG as an ancillary assertion; repair does not change Theorem A.

Location: thmA_proof.md §6 U4, final clause “because ... gamma >= 1/(log N)^{O(1)}”, when used together with the explicit example for S in §1 and U5.

Here is a quantitative reason, rather than merely the observation that FM is a lower bound. For the explicit choice

  S^2=ceil(L^{theta/(2a)}) ceil(log N)^{C3},

write T=log N and lambda=log T. The smoothness expansion

  log rho(X_i)=-u_i(log u_i+log log u_i-1+O(log log u_i/log u_i)),
  u_i=(xi_i/a)sqrt(T/lambda)+O(1),

gives

  log m= -[theta/(2a)]sqrt(T/lambda) log lambda
          +O(sqrt(T/lambda)).

The leading Lambda terms cancel, but the displayed negative secondary term dominates every multiple of log log N. Thus m is smaller than every inverse power of log N while still equal to L^{o(1)}.

A matching upper bound sufficient to apply this observation follows from ingredients already available here. The smoothed one-point coefficient law has maximum mass O(1/(HV^{d-1})), with total-variation error O(V/H). For fixed g,P satisfying the rational filter, enlarge the value window to a fixed multiple of X*. Lemma W then bounds its averaged smooth mass by

  [C/(HV^{d-1})] Psi(CX*)/phi(t)+O(V/H)
       << rho_F (log N)^C+O(V/H).

The rational first moment is at most C S^2 rho_r times a polylogarithm, by the fibre upper bound and Psi(Xr). Consequently

  gamma <= E N_l <= m(log N)^C+N^{-c}L^{O(1)}

for some c>0. This rules out an inverse-polylogarithmic lower bound for gamma for this explicit S. The required repair is simply to use gamma^{-O(k)}=L^{o(1)}, proved in item 14. Fixed k is essential. The stated exponent does not require the stronger assertion.

16. Filtered PZ, capped search, and rejection interface — CORRECT WITH FIXES.

Location: hypD_results.md §4 Prop PZ; pcest_results.md §1 Lemma 1.3; lv6_results.md §6.2 Prop BS; thmA_proof.md §4 Steps 2, 5–6.

PZ is the exact Cauchy–Schwarz inequality for any nonnegative integer count. It must be applied to the filtered count used by the algorithm. The regularity lemma makes the complete candidate enumeration O(S^2): an irregular lattice has no primitive candidate with b>V/2, and a regular lattice has O(S^2) box points. Batch smoothness recognizes every qualifying value exactly. Thus the analytic existence event is the event consumed by the search, apart from successful factor exits, which only help.

Subtracting rejection probabilities at the bounded success-indicator level is the correct order. The bounds exp(-log H/(13d log log H))+O(1/H+1/W) are L^{-omega(1)}, hence negligible against L^{-o(1)}. They must not be subtracted after multiplication by the maximum relation count.

One wording issue remains: removing proper-factor exits does not automatically leave the conditional law of nonterminating fields with gamma/2. The rigorous argument is a two-case one. Either factor exits carry at least a fixed fraction of the raw success mass, and the algorithm succeeds quickly through those exits, or the remaining accepted fields retain that fraction. This repairs the assertion without adding a hypothesis or an L-power.

17. Theorem R is not on the critical path — CORRECT.

Location: ldbound_results.md §5 remarks; thmA_proof.md §4 Step 5.

The needed dependency chain is FM + Kern + W + DS + L2_fam -> LD* -> PZ. The partner's rational smoothness is estimated directly along a fibre in LD* Step 5. Neither the statement nor the proof of rational pair-correlation Theorem R is needed. Lemma Q's determinant observation and Lemma DS happen to be located in the same section as R; that does not make R an input. I do not use R to justify extending the allowed S.

18. Exponent optimization and k=13 — CORRECT.

Location: thmA_proof.md §1 and §4 Step 10; src/thmA/params.py.

Eliminating h at zero slack gives the upper envelope of three affine functions:

  theta_0(mu)=max(2mu+(d+1)v, 3mu+v, 1+(2-d)mu+d v).

The first two increase and the third decreases. Their relevant intersections give

  k<(d^2+1)/4:
    mu=(1-v)/d,
    theta_0=2/d+[(d^2+d-2)/d]v;

  k>=(d^2+1)/4:
    mu=[1+(d-1)v]/(d+1),
    theta_0=3/(d+1)+[(4d-2)/(d+1)]v.

Positive slack adds eta0 to theta, provided the chosen mu also satisfies mu>=v+eta0. These formulas are zero-slack minima, not a claim that the same mu remains feasible for every larger eta0.

Independently minimizing this affine envelope with exact rational arithmetic agreed with params.py for every admissible k<=500 and odd 3<=d<=51. This check did not call the module's optimization routine to construct the independent candidates.

At d=9:

  k=12: mu=47/432, h_*=1/24, theta_0=23/54,
        c=1.0082818917... >1;
  k=13: mu=17/156, h_*=1/26, theta_0=16/39,
        c=0.9972018709... <1.

The zero-slack minimum decreases with k, so no smaller admissible k improves on 1. For k=13 the extra condition is

  0<eta0<sqrt(2)-1-16/39=0.0039571521...,

which is compatible with mu>=v+eta0.

The total cost exponent is max(2a,a+1/sqrt(2),a+theta/(2a)). For theta<1 its minimum occurs at a=theta/sqrt(2): the two relation terms meet there, and 2a is strictly smaller. This yields c=(1+theta)/sqrt(2). The shape-test benchmark is supported by the complete-factorization theorem of Lenstra–Pomerance; replacing an input of size N by O(sqrt N) multiplies its L-constant by 1/sqrt(2). [Lenstra–Pomerance, main theorem, p. 483](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf).

For d=9 the infimum of theta is 3/10, and the infimum of the RUNNING-TIME CONSTANT is

  (13/10)/sqrt(2)=0.9192388155425... .

It is an infimum over fixed parameters. It is not attained by any fixed k,eta0, and no growing-k diagonalization has been established here.

19. Uniformity and scope of the recommendation — CORRECT WITH FIXES.

Location: thmA_proof.md §2 and §6; ldbound_results.md §§2–5.

Every N^{o(1)} divisor or density factor that appears in my reconstruction is absorbed by a fixed polynomial margin: mu-v, h-v, eta0/2, 2v, h+(d+1)v, or h-(5-d)v. The remaining losses are fixed polylogarithms or L^{o(1)}. In particular, moving the value window by B^{O(1)}, summing the grid, and summing tau(q)/q do not introduce a fixed L-power. The explicit S supplies the needed uniform relation-supply bound through item 14, with item 15's weaker gamma estimate used throughout the cost ledger.

The analytic arguments use no squarefreeness of N, no assumptions on its residue classes, and no exclusion of an unspecified exceptional set of N. The only coprimality hypotheses come from the stated preprocessing, the rational smoothness event, and l not dividing N. Fixed parameters give all the required positive margins.

I have not independently certified the external Stark/Lagarias–Odlyzko statements or all the algorithms and cost claims in LV6: those belong to the other half of this review. I therefore make no claim that this report alone certifies the entire factoring theorem or its historical priority. Within the assigned analytic chain, I found corrections rather than a major mathematical gap.

WHOLE-THEOREM VERDICT: ACCEPT WITH FIXES
WEAKEST LINK: pcest_results.md §2.2 Part 1, importing the absolute-character majorant W from the proof of LV Lemma 8.5. It requires more than the lemma's stated maximum-residue conclusion, including uniform treatment of induced and exceptional characters. The checked fixed-power regime supports the import, but this is the step most deserving a fully rewritten proof before claiming a new rigorous factoring bound.
