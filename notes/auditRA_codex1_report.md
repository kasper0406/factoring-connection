**Independent referee report on RA and Corollary P**

I find the recursive implication mathematically sound, conditional on the supplied proof of Theorem A at its original parameter being correct. Changing the smoothness constant introduces no new range obstruction in the estimates examined. The expected-time recurrence and its solution also check out. Corollary P follows after making its parameter sequence and uniform implementation explicit.

This is conditional acceptance with fixes, not independent certification of Theorem A. In particular, accepting the base proof is doing substantial work: the restricted smooth-modulus estimate, its transfer into FM, the second moment, and the number-field post-processing remain inherited mathematical inputs. The published literature alone does not establish this project's relation supply.

I respected the exclusions in the prompt. I did not open the prohibited reports, the original recursive-results file, or other material under `codex/` outside this authorized scratch directory. Source verdict/review sections were excluded from the reading copies; assessment labels and candidate cards played no part in the assessment. Locations below are original source line numbers, not line numbers in those copies. Numerical checks were independently calculated in [audit_calculations.py](/Volumes/git/factoring-connection/codex/scratch/auditRA_1/audit_calculations.py).

The classifications distinguish a mathematical defect from a repair to the formulation or implementation. A status of CORRECT for a source estimate below concerns its extension to fixed general `a` and its use by RA, under the permitted base-proof assumption.

**1. The freedom to choose a — CORRECT WITH FIXES.**

Locations: [blind claims:39](/Volumes/git/factoring-connection/notes/auditRA_blind_claims.md:39), lines 49–88; [Theorem A:25](/Volumes/git/factoring-connection/notes/thmA_proof.md:25), lines 25–70 and 202–253; [ldbound:127](/Volumes/git/factoring-connection/notes/ldbound_results.md:127).

Write `t = log N`, `l = log log N`, `Lambda = sqrt(t l)`, and `L = exp(Lambda)`. Set

\[
 \nu=\frac1{4k},\qquad
 \xi_r=\mu-\nu>0,\qquad
 \xi_F=\mu+h+d\nu,\qquad
 \theta=\xi_r+\xi_F.
\]

For each fixed `a > 0`, independently of whether it is the old optimum,

\[
 B=L^{a+o(1)},\quad S=L^{\sigma+o(1)},\quad
 \sigma=\frac{\theta}{4a}>0,\quad
 \rho_r=L^{-\xi_r/(2a)+o(1)},\quad
 \rho_F=L^{-\xi_F/(2a)+o(1)}.
\]

The factors `S` in the arguments of the smooth-number densities change their logarithms by `O(Lambda)`, and hence do not change these leading constants. In particular

\[
 m:=S^2\rho_r\rho_F=L^{o(1)}.
\]

There is no requirement that `a` be uniformly bounded away from zero as the parameters range over different algorithms. There is a requirement that it be fixed while taking the asymptotic limit for one algorithm. The latter is precisely what RA and the portfolio proof use.

The potentially dangerous constants can be chosen as follows:

* Choose the FM truncation constant `C1 > xi_F/(2a)`, with a fixed positive margin. Then `1/omega`, where `omega = L^C1`, is negligible against `rho_F` times any fixed inverse polylogarithm. The required `H > d omega V` holds because `h - nu > 0` has a polynomial margin.
* In L2_fam choose an integer
  \[
    C_0>\frac1a+\frac{\theta}{a^2}.
  \]
  This gives `B^C0 >= 64 S^2 L/(rho_r rho_F)` eventually. `C0` can become large when `a` is small; it is still a fixed constant. Also `B^C0 = N^{o(1)} < Y`, and is smaller than the polynomial lengths needed in divisor switching.
* For Theorem R alone, the weaker choice `C0 > theta/(2a^2)` suffices to cover `q <= 2 S^2`.
* In the separate proof of Lemma 8.5′, `K = 6 + ceil(2 xi_max/a^2)` is fixed. Thus
  \[
    B^K=N^{o(1)}<Y=N^{2\nu+o(1)},
    \qquad B^K\le z^{\eta_D}
  \]
  eventually for every argument `z >= x^(1/2)` used there.

Here are the geometric and distribution margins explicitly:

\[
\begin{aligned}
 M/V &=N^{\mu-\nu+o(1)},& H/V^2&=N^{h-2\nu+o(1)},\\
 W/(2MV_{\max})&=N^{\eta_0/2+o(1)},&
 X_F^{1/2}/(MV)&=N^{(h-\mu+(d-2)\nu)/2+o(1)}.
\end{aligned}
\]

The respective margins are at least `eta0`, `eta0`, `eta0/2`, and `eta0/2` in the exponent where the h-floor is used. Consequently `B < V/2`, `V < ell`, `2 S^2 < Y`, and `S < sqrt(ell)` all hold eventually. These remain true for both very small and very large fixed `a`.

For `x = N^(xi+o(1))` in the fixed compact range of the FM grid,

\[
 u=\frac{\log x}{\log B}
   =\left(\frac{\xi}{a}+o(1)\right)\sqrt{\frac t l},
 \qquad
 \frac{u/\log^2 u}{\log\log x}\longrightarrow\infty.
\]

Thus `B > (log x)^C` for every fixed `C`, `B < x^(1/log log x)`, and the Hildebrand short-interval range holds. No endpoint failure appears when `a` changes. The large-conductor error, divided by the required restricted-family main term, has polynomial savings: the two relevant ratios are `N^(q0-xi_F/2+o(1))` and `N^(-eta+o(1))`, where `q0 = mu + nu < xi_F/2` and `eta > 0` is fixed. The densities lost in the denominator are only `N^{o(1)}`. These are also the ranges in [Harper, Proposition 2 and §3.3](https://arxiv.org/pdf/1208.5992), rather than the much narrower range of some older smooth-number results.

The individual estimates have the following status.

| Estimate and exact location | Classification for RA | Reason |
|---|---|---|
| FM: `pcest_results.md:304–406`; transfer `ldbound_results.md:214–240` | CORRECT | Its smoothing, kernel, moving-residue, and level conditions have the margins above. `C1` may depend on `a`. The needed conclusion is averaged over both `ell` and the field. |
| Lemma DS: `pcest_results.md:448–466` | CORRECT | `C0` is allowed to be any fixed constant. The interval after removing `B^O(1)` factors is still polynomial in `N`; the divisor-switching argument covers non-reduced residues as well. |
| Theorem R: `pcest_results.md:478–520` | CORRECT | Its extra requirement holds: `log(S^2 rho_r^2) = ((xi_F-xi_r)/(2a)+o(1)) Lambda`, and `xi_F-xi_r = h+(d+1)nu > 0`. Any smaller fixed positive exponent gives the advertised margin. R is not needed in the LD* route. |
| Kern: `pcest_results.md:553–592` | CORRECT | `H/V^2 = N^(h-2nu+o(1))` tends to infinity. This also absorbs the fixed factors 2 or 4 in the more general norm convention. The explicit divided minor `(bb')^(d-2) Delta` supplies the wedge bound. |
| Spread Lemma W: `ldbound_results.md:263–279` | CORRECT | `log(W/t_mod) >= (eta0/3) log N` eventually, even after the `S` factor. Thus the Brun–Titchmarsh constant `24/eta0` is unaffected by freeing `a`. |
| L2_fam: `ldbound_results.md:303–350` | CORRECT | The larger `C0` above is permissible. Its other errors have exponents `-h-(d+1)nu`, `-2nu`, and `(5-d)nu-h`, all strictly negative. Fixed `B`-powers and `S`-powers are subpolynomial. |
| LD*: `ldbound_results.md:396–458` | CORRECT WITH FIXES | Its stated hypothesis uses the implicit density-defined `S`. Its actual proof needs only the displayed size/range inequalities and therefore extends to the explicit `S`. State that extension; do not retain the assertion `E N >= 1` at line 459 for the explicit choice. |
| Theorem 4.3 / Prop PZ: `pcest_results.md:644–660`; `hypD_results.md:340–344` | CORRECT WITH FIXES | Use PZ's inequality directly, without its sufficient condition `E N >= 1`. If `E N >= m/P0` and `E N(N-1) <= P0 m^2`, then `gamma >= 1/(P0^3+P0/m) = L^{-o(1)}`. Here `P0` is a fixed polylogarithm. |
| S / S′: `hypD_results.md:161–217` | CORRECT | These are finite set-counting statements. They have no smoothness-constant restriction. The product sets and their windows still meet the size hypotheses. Preprocessed `N` is odd, so `2 X^2 < N` for `X = floor(sqrt(N/2))`. |
| K and field selection: `hypD_results.md:221–283` | CORRECT | `gamma = L^{-o(1)}` implies that fields with `g_A >= gamma/2` have polynomially many good primes, eventually exceeding `4k^2`. All powers of `gamma^{-1}` have fixed exponents because `k` is fixed. |
| Lemma Omega: `hypD_results.md:153–159`; stronger preprocessing `thmA_proof.md:188–193` | CORRECT | After trial division to `B`, `omega(N) <= log N/log B`, so `log(2^omega)/Lambda <= (log 2)/(a log log N)+o(1)`, which tends to zero. |
| Lemma 8.5′: `pcest_results.md:257–299`; separate proof `lem85_results.md:82–112, 329–412` | CORRECT | The polynomial gaps give `x^lambda Y <= Q <= x^(1/2-lambda)`. The larger fixed `K` still gives `B^K < Y` and `B^K <= z^eta`. The small-conductor polylogarithms are absorbed because `u/log^2 u >> log log x`. |
| LV6: `lv6_results.md:114–147, 213–383`; transfer `ldbound_results.md:519–537` | CORRECT | Fixed degree and `N^O(1)` field heights are sufficient. For every fixed `a`, `B/log B >> log N`, so the Hensel prime exists; character primes of logarithm `O_d((log N)^2)` exceed `B`. The costs are `L^(2a+o(1))` and `L^(a+o(1))` as stated. |
| Prop BS: `lv6_results.md:551–567` | CORRECT WITH FIXES | Its general cost is `soft-O(B+S^2 log N)`. The shortcut `B <= S^2` is false for `a^2 > theta/2`. RA correctly restores the `B` term. A historical time cap of order `S^2` must not truncate this batch operation. The pseudocode has a point-count cap, which remains valid. |

Two source imprecisions do not create a new general-`a` gap. The per-`ell` FM conclusion displayed at `ldbound_results.md:238` is stronger than the bad-modulus argument establishes; Theorem A's roadmap correctly uses only the joint average. Also the Tenenbaum condition is the corrected one in `thmA_proof.md:213–215`, not the typo in `pcest_results.md:243`. The corrected condition holds for all these fixed `a`, since `omega(t_mod) = O(log N)` whereas `B^(c/log(1+u))` eventually exceeds any fixed power of `log N`.

The computable choice of `S` need not make `m` bounded below by a polylogarithm. De Bruijn's next term can make `m` smaller than every inverse polylogarithm. The direct PZ calculation above handles that; replacing `L^{-o(1)}` by a polylogarithmic lower bound would be unjustified.

The corrected cost ledger is therefore

\[
 D(N):=\mathbb E[\#\mathrm{DESCEND\_ONCE}]=L^{a+o(1)},
 \qquad
 T_{SQ}=L^{\max(a,\theta/(2a))+o(1)},
\]

and

\[
 P(N):=\mathbb E[\text{work outside the shape oracle}]
 \le L^{\max(2a,\ a+\theta/(2a))+o(1)}.
\]

There are at most `2k` special-q searches per descent attempt, whether or not that attempt supplies a row. This gives the displayed bound directly and avoids dependence on K(iii)'s more delicate waste-ratio statement. Sampling, certification, gcds, and Gauss reduction add only the already allowed factors. With an LP oracle, the unrestricted ledger would instead be `max(2a, a+1/sqrt(2), a+theta/(2a))`.

**2. Recursion, stopping times, and finiteness — CORRECT WITH FIXES.**

Locations: [blind claims:108](/Volumes/git/factoring-connection/notes/auditRA_blind_claims.md:108), lines 108–159; [Theorem A:111](/Volumes/git/factoring-connection/notes/thmA_proof.md:111), lines 111–134 and 147–156.

The new recursive inputs satisfy `|u0|, |v0| <= sqrt(2n) < n` for `n > 2`. They need not have the intended product-of-primes shape: the oracle completely factors them first. It consequently must work on arbitrary smaller integers, including primes, prime powers, and highly composite integers. Strong induction supplies exactly that guarantee.

Define `FACTOR(1) = []` with bounded cost, and include that cost in the envelope, for example

\[
 F(x)=\max\left(1,\sup_{1\le n\le\lfloor x\rfloor}
        \mathbb E[T(n)]\right).
\]

The omission of 1 from the displayed definition in the notes is harmless after this adjustment. Empty products must be permitted by the “at most k primes” rule. A prime power returned by an inner call is valid factorisation output; it usually causes the outer distinct-prime shape test to fail, as it should. Perfect powers are processed before SPLIT. For unit `R` and `n > 2`, a zero coordinate would force a lattice vector of length at least `n`, contradicting the Gauss bound. Thus the oracle is not asked to factor 0.

For one input `N`, the ordinary factor tree has at most `2 Omega(N)-1 <= 2 log_2 N-1` nodes. This tree counts recursive splitting of divisors of `N`. Auxiliary shape calls have their own complete factor trees inside `F(sqrt(2N))`; counting them as divisor nodes of the original `N` would be wrong. The notes keep these two accountings separate.

The adaptive Wald bound at lines 108–114 is correct. The event that call `i` is made and its input are known before that call's fresh tape is sampled. Tonelli and conditional expectation give

\[
 \mathbb E\sum_i {\bf1}_{\{i\text{ is made}\}}C_i
 \le F(y)\sum_i\Pr(i\text{ is made}).
\]

No independence between call count, input, and cost is required. Inner field rejections, unsuccessful rows, failed dependencies, and retries are all included in `T(m)`. There is no missing multiplier for those failures.

Condition on the value at each ordinary factor-tree node, take the supremum over its possible inputs, and apply this bound to the at most two shape calls per descent. One obtains

\[
 F(x)\le 2\log_2 x\left[P(x)+2D(x)F(\sqrt{2x})\right],
\]

after absorbing bounded base costs into `P`. This is an upper bound, not an asserted tight recurrence. Suprema over smaller inputs preserve the displayed `L` exponents: split the supremum into a fixed finite initial range and the range where the error in the exponent is uniformly small.

The finiteness induction is legitimate once the small-input implementation described in item 3 is supplied. First assume all smaller integer inputs have finite expected cost and canonical, correct outputs. They then form admissible shape oracles for the current input; the finite outer expectation bounds apply. This proves finiteness at the current input without presupposing termination of the same-size call. Almost-sure termination follows from finite expectation.

Theorem A's attempt count is independent of which exact oracle is used. The oracle's outputs determine precisely the same shape acceptance events. Retaining a field on its first successful attempt gives the `q_A`-biased law and the bound on `E[1/q_A]`; it does not retain the fastest field. Factor exits only shorten the relevant collection process, with the absorbing accounting at `thmA_proof.md:261–268`. The conditional fruitfulness bound controls the number of fresh-row rounds; independence of their costs is unnecessary.

**Lemma F — CORRECT.** The proof at `blind claims:131–148` has the required margin, including its constant `K`. For `b = c+eta` and sufficiently large `x`,

\[
 L_{\sqrt{2x}}^b\le e L_x^{b/\sqrt2}.
\]

With `eta' = eta(1-1/sqrt(2))/4`, the recursive term in the induction is at most

\[
 KeL_x^{\delta+\eta'+(c+\eta)/\sqrt2}
 \le KeL_x^{c+\eta-3\eta'}
 \le (K/2)L_x^{c+\eta}.
\]

The additive term is at most `(1/2)L_x^(c+eta)`, and `K >= 1`, so their sum is at most `K L_x^(c+eta)`. Starting at `X >= 4`, the successive endpoints `Y -> Y^2/2` tend to infinity. This establishes the induction over the entire half-line, rather than just a sparse sequence. The final constant `K_eta` can be absorbed into `L^{o(1)}`. Polynomial logarithmic factors in the original recurrence contribute to the vanishing exponent and do not change the answer.

**The product-size aside — WRONG as a universal statement, with no effect on RA.** At `blind claims:219–220`, Gauss reduction does not imply `|u0 v0| asymp N`, nor must every two-term split have total logarithmic size `log N`. For `R = 1`, `(u0,v0) = (1,1)` is a counterexample. The cost proof needs only the upper bounds on each coordinate. Delete the universal wording; a statement about typical sizes would need a separate probabilistic qualification.

**3. Oracle substitution and transcript invariance — CORRECT WITH FIXES.**

Locations: [blind claims:92](/Volumes/git/factoring-connection/notes/auditRA_blind_claims.md:92), lines 92–106 and 182–196; [DD:99](/Volumes/git/factoring-connection/notes/lemmaM_results.md:99), lines 99–119.

OS is correct for the uncapped logical algorithm specified in the notes. Couple the outer random tapes. Each oracle call returns the same sorted prime list with multiplicities. Induction over successive outer instructions then gives identical outer values and branches, even if the calls have different, random running times. A countable union of probability-zero nontermination events remains null. Random time is irrelevant because no outer instruction uses it.

Canonical ordering matters if prime-processing order can change which factor exit occurs first. It should be an explicit output convention of every recursive call, not merely a convention in the proof. Certification and sorting cost only polynomial logarithmic overhead in the output size here.

T_rec also holds if the recursive calls' full internal transcripts, including their running times, are exposed. Conditional on the outer squares, each inner computation is a function of its input and its independent tape. Replacing an outer root by `epsilon r` changes neither auxiliary input. Coupling the inner tapes therefore leaves the entire inner computation unchanged. Its timing can reveal information about the square `R`, which is already public; it cannot reveal which root on that square's fibre was drawn. Unconditional independence from the outer roots is neither true nor required.

The original algorithm draws uniformly from `[1,N-1]`, with a factor exit on a nonunit (`thmA_proof.md:148`). Thus its successful unit draws have the exact fibre law required by DD. The older `[2,N-2]` sampling caveat in `lemmaM_results.md:92–93` is not inherited by this pseudocode. Fresh rows after a trivial final gcd are necessary and are present. Inner DD arguments use their own roots conditionally on their own input; no independence of the inner input from the outer squares is assumed.

The actual pseudocode has no timeout or race around the shape oracle. Historical caps on smoothness recognition do not authorize one there. The complete-box BS search must cover all counted candidates, including its `B` work. Timing instrumentation in toy code is outside the logical transcript asserted by OS.

**The N0 branch needs a uniform implementation repair.** `thmA_proof.md:117` and U3 at lines 342–343 use an analytic threshold. RA's own hypotheses at `blind claims:53` correctly add `a` to its dependencies, whereas “every other line is unchanged” does not explicitly update that branch. One must use `N0(d,k,mu,eta0,a)`. More significantly for the portfolio, separate existence of a finite threshold for each parameter set is not enough to compile a computable sequence of programs containing those thresholds.

An explicit repair avoids this question entirely. At each invocation, run the recursive factoring core and deterministic complete trial division in equal shares, returning the first verified complete factorisation in canonical order. Give the core arbitrary safe failure behaviour for small inputs where its sampling sets or prerequisites are unavailable. This wrapper is one uniform program and makes every finite input terminate. For sufficiently large inputs, its cost is at most twice the cost of the core using the same wrapped smaller-input oracles. The recurrence merely acquires fixed factors, which Lemma F absorbs.

This wrapper introduces a race at the wrapper level. Do not claim OS's identical outer transcript for that race. Apply OS to the uncapped SPLIT core, and dominate the wrapper's time by that core's time. A wrapper inside a child oracle still returns the same canonical factorisation and is fully compatible with OS in its caller. Thus the repair does not require a timing-invariance assertion that would be false.

There is a genuine distinction here: `lem85_results.md:438–460` itself does not give an entirely unqualified proof of effectivity of every analytic constant. That is an inherited implementation issue, not an obstruction to the conditional asymptotic argument. The wrapper resolves it without using an uncomputed threshold. Also restrict algorithm parameters to computable values, or formulate arbitrary real `a` as an exponent bound realised by a computable choice. Literal exact floors involving arbitrary noncomputable real parameters are not an algorithm. Computable rounding with the same `L` exponent suffices; no test of exact equality with a transcendental threshold is needed.

**4. Corollary P — CORRECT WITH FIXES.**

Locations: [blind claims:224](/Volumes/git/factoring-connection/notes/auditRA_blind_claims.md:224), lines 224–244.

The example at lines 225–228 is incomplete. Holding a small positive `eta0` fixed gives `theta_j >= eta0`, preventing `cstar_j -> 0`. Merely checking finitely many parameter sets would not prove an infinite sequence admissible. “A rational near a-star” also needs an accuracy convention that ensures the resulting constants tend to zero.

A fully explicit replacement is

\[
\begin{aligned}
 d_j&=2j+7,&
 k_j&=\lceil d_j^2/4\rceil+j,&
 \nu_j&=1/(4k_j),\\
 \mu_j&=\frac{1+(d_j-1)\nu_j}{d_j+1},&
 \eta_{0,j}&=\frac1{100d_j^2},&
 h_j&=h_*(\mu_j)+\eta_{0,j}.
\end{aligned}
\]

Here `j >= 1`, so `d >= 9`, and `4k = d^2+2d-11`. Since `k >= (d^2+1)/4`, the level and window floors meet at the displayed `mu`. Direct substitution gives

\[
 \theta_j=\frac3{d_j+1}
       +\frac{4d_j-2}{d_j+1}\nu_j
       +\frac1{100d_j^2}=O(1/d_j).
\]

For every such `d`, `nu < 1/d^2` and

\[
 \mu-\nu=\frac{1-2\nu}{d+1}>\eta_0,
 \qquad
 \theta<\frac3{10}+\frac4{81}+\frac1{8100}<1.
\]

The remaining strict h-floor inequalities follow from adding `eta0`; `mu+h < theta < 1` and the `w < N` condition follow as in Theorem A. This proves admissibility for every `j`, rather than extrapolating the finite checks. The independent rational-arithmetic checks through `j = 1000` are supplementary checks of these formulas.

Use the computable algebraic optimum `a-star_j`, with harmless computable rounding of `B`, or rational approximations with relative error tending to zero. An even simpler rational choice sufficient for P is `a_j = 1/ceil(sqrt(d_j))`. Then `a_j -> 0` and `theta_j/a_j -> 0`, so `c(a_j) -> 0`. The optimum is unnecessary for the existence conclusion. All parameter construction can be carried out with overhead polynomial in `j`; all setup inside a run is charged against that run's budget. Use the uniform wrapper of item 3 for every `A_j`.

**Restart calculation — CORRECT.** Put `M_j = max(1,E T_j(N))` and

\[
 s_0=j+\left\lceil\log_2(2M_j)\right\rceil.
\]

In phase `s0+r`, the budget for `A_j` is at least `2^(r+1) M_j`, so Markov gives failure probability at most `2^(-r-1)`. Fresh runs in distinct phases are independent. Reaching phase `s0+m` therefore has probability at most

\[
 \prod_{r=0}^{m-1}2^{-r-1}=2^{-m(m+1)/2}.
\]

The phase budgets sum to `2^s-1`. Including ordinary polynomial bookkeeping or universal-simulation overhead, the tail series is bounded by a fixed polynomial in `s0` times

\[
 2^{s_0}\left(1+\sum_{m\ge1}2^{m-m(m+1)/2}(1+m)^C\right).
\]

Consequently, for each fixed `j`,

\[
 E T(N)\le O\left(2^jM_j\,
       \operatorname{poly}(j+\log(M_j+2)+\log N)\right).
\]

The polynomial factor is `L^{o(1)}` when `j` is fixed. This more explicit expression is what the notes' unspecified “polylog” should mean. Expected runtime is finite for every `N`, even below every analytic asymptotic threshold, because any one wrapped `A_j` has finite mean there. Every returned answer is verified.

For every fixed `j` and sufficiently large `N`, the RA estimate now gives

\[
 \limsup_{N\to\infty}\frac{\log E T(N)}{\Lambda_N}
 \le c(a_j).
\]

Letting `j` tend to infinity proves the stated `exp(o(Lambda_N))` upper bound. No simultaneous estimate for increasing `j` as a function of `N`, and no uniform asymptotic threshold, is used.

**Finite means and the extra-moment issue — CORRECT for this schedule.** Finite `E T` alone does not imply finite `E[T log T]`: a distribution assigning probability proportional to `1/(n^2 2^n)` to `T = 2^n` is a counterexample. Thus charging an unexamined `log T` overhead to a live dovetail can be invalid. The restart schedule bounds the total elapsed-time tail afresh, using Markov at increasing budgets, and never makes that invalid inference.

**Effectivity wording — CORRECT WITH FIXES.** The portfolio is computable; its proof does not require computing the crossover thresholds. The argument as stated supplies an asymptotic limit without an extracted modulus of convergence. It does not prove that every modulus is noncomputable. If effective error bounds and thresholds are supplied uniformly by the source proofs, an effective, extremely poor rate could in principle be extracted. Replace “the rate is not effective” by “no effective rate is supplied by this argument” unless a stronger claim is separately proved. Nor does `L[1/2,o(1)]` imply an `L[alpha]` bound for any fixed `alpha < 1/2`.

**5. Literature, barriers, and the heuristic comparison — CORRECT; no cited impossibility theorem refutes RA.**

The literature supports a barrier for particular random-squares models, not a universal lower bound for all algorithms that happen to factor a square-root-size auxiliary integer.

For independent random integers of size `x = N^xi`, the elementary optimistic smoothness ledger is

\[
 \max\left(a+\frac{\xi}{2a},2a\right).
\]

Its minimum is `sqrt(2xi)`, attained at `a = sqrt(xi/2)`. Thus the familiar constant 1 occurs at `xi = 1/2`. This calculation presupposes obtaining the final square dependence from the sampled integers themselves. Faster recognition of smoothness alone cannot remove its relation-count cost. Pomerance explicitly ties the common `L[1/2,1]` behaviour of the methods in his discussion to their square-root-size smoothness auxiliaries, and then explains how NFS changes the auxiliary sizes. [Pomerance, *The role of smooth numbers in number theoretic algorithms*, §5](https://math.dartmouth.edu/~carlp/PDF/paper104.pdf).

[Croot–Granville–Pemantle–Tetali, *On sharp transitions in making squares*](https://annals.math.princeton.edu/2012/175-3/p10) study independent uniform draws from `[1,x]` and the first subset having an integer-square product. Their bounds sharpen Pomerance's random-squares threshold. [Balister–Bollobas–Morris, *The sharp threshold for making squares*](https://annals.math.princeton.edu/2018/188-1/p02) subsequently identify its sharp asymptotic location. Neither result is a lower bound for number-field descent with externally supplied relations eliminating the large primes.

The remark on [Lenstra–Pomerance, p. 509](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf) states that Algorithm 10.1 has expected time equal to `L_n[1/2,1+o(1)]` and that choosing its parameter `y` differently does not improve it. It does not state that all rigorous factoring algorithms, or every algorithm factoring square-root-size auxiliary numbers, require this time. Reading that remark as a general impossibility theorem is incorrect.

RA does not demand that the Gauss halves be `B`-smooth. It accepts products of a fixed number of much larger primes, factors arbitrary rejected halves as well, and eliminates accepted primes using special-q number-field relations. The smooth cofactors in that descent have combined exponent `theta`, not `1/2`. Therefore neither the displayed model ledger nor the random-squares subset theorem applies to its complete relation construction. This is a legitimate distinction if the claimed first and second moments are true; calling an operation “descent” does not by itself prove those moments.

The published [Lee–Venkatesan paper](https://arxiv.org/pdf/1805.08873), Theorem 2.1 and §7, establishes congruences of squares for its randomized NFS and separately discusses the difficulty of proving nontrivial factor extraction. It does not give the particular fruitful fixed-degree descent supply required here. Thus the novelty and mathematical burden lie upstream of the elementary recurrence. I cannot explain historical non-discovery by a theorem forbidding the recurrence. Accepting a new, sufficiently strong rigorous supply theorem would indeed have major consequences.

For a fair heuristic check, take the balanced nonmonic degree-`d` NFS construction, with rational and algebraic values each of polynomial exponent `1/(d+1)` in a subpolynomial sieve box. The nonmonic parameter improvement from `N^(2/d)` to `N^(2/(d+1))` is discussed in [Buhler–Lenstra–Pomerance, §12.9](https://www.researchgate.net/publication/225539588_Factoring_integers_with_the_number_field_sieve). The following is my fixed-degree heuristic derivation, not a claimed published rigorous bound.

The joint smoothness penalty is `L^(1/((d+1)a)+o(1))`. To collect `B` relations, the search costs

\[
 L^{a+1/((d+1)a)+o(1)},
\]

while sparse linear algebra costs `L^(2a+o(1))`. The minimum is

\[
 c_{\mathrm{NFS},d}=\frac2{\sqrt{d+1}},
 \qquad a=\frac1{\sqrt{d+1}}.
\]

If one imposes the older monic base-`m` convention, the analogous constant is `2/sqrt(d)`. The convention matters when `d` is fixed. With the balanced convention, RA's infimum satisfies

\[
 \frac{c_{RA,d}}{c_{\mathrm{NFS},d}}
 =\frac{\sqrt{3(1+\sqrt2)}}2>1.
\]

RA therefore remains slower than this heuristic at every displayed degree. Its limit as the fixed degree increases is not a contradiction with heuristic NFS.

**Prior art — limited finding, not a novelty certificate.** Auxiliary cofactorisation and recursive descent are established ideas; routine recursive splitting of the factors of `N` should not be confused with this self-application on unrelated auxiliaries. An [INRIA 2009 activity report, §5.4 on TIFA](https://citeseerx.ist.psu.edu/document?doi=835d3a51b88211d31afb9a34dcb6a4a2e8ecba31&repid=rep1&type=pdf) explicitly considers recursive NFS for its smaller byproducts and rejects it on practical efficiency grounds. Bernstein's [smooth-parts paper, history following Algorithm 2.1](https://cr.yp.to/factorization/smoothparts-20040510.pdf) supplies a separate precedent for recursive batch auxiliary processing, with a different recurrence. My searches did not locate an earlier theorem giving RA's exact self-recursive exponent transformation. That search result does not establish originality.

**6. Constants and fixed-depth iteration — CORRECT.**

Locations: [blind claims:150](/Volumes/git/factoring-connection/notes/auditRA_blind_claims.md:150), lines 150–178; the two supplied constants modules. The following calculations are independent of the implementation of those modules.

Apply Lemma F with

\[
 p=\max(2a,a+\theta/(2a)),\qquad \delta=a.
\]

Since

\[
 \frac1{1-1/\sqrt2}=2+\sqrt2>2,
\]

the resulting upper-bound constant is exactly

\[
 c(a)=\max\left(a+\frac{\theta}{2a},(2+\sqrt2)a\right).
\]

The first branch decreases up to `sqrt(theta/2)` and the linear branch increases. Their crossing occurs before that minimum, at

\[
 a_*^2=\frac{\theta}{2(1+\sqrt2)}
        =\frac{(\sqrt2-1)\theta}{2}.
\]

At this point

\[
 c_*=(2+\sqrt2)a_*,\qquad c_*^2=(1+\sqrt2)\theta,
 \qquad \frac{2a_*}{c_*}=2-\sqrt2.
\]

These identities are exact. They optimize the proved upper-bound ledger; they are not lower bounds for the algorithm's actual running time. Also

\[
 c_*<1\iff \theta<\sqrt2-1\iff(1+\theta)^2<2.
\]

For `0 < theta < 1`, the comparison with Theorem A is

\[
 c_*<\frac{1+\theta}{\sqrt2}
 \iff 2\sqrt2\theta<1+\theta^2
 \iff \theta<\sqrt2-1.
\]

Thus the minimal admissible `k` for a constant below 1 is unchanged. Above that threshold, self-recursion at its own optimum can have a worse bound than the original LP-oracle construction; taking the better algorithm is legitimate.

If the bottom-level oracle costs `L_m[1/2,c+o(1)]`, one more fixed level has optimized constant

\[
 g(c)=\min_{a>0}\max\left(2a,a+c/\sqrt2,a+\theta/(2a)\right)
 =\begin{cases}
   (c+\theta/c)/\sqrt2,&c^2\ge\theta,\\
   \sqrt{2\theta},&c^2<\theta.
  \end{cases}
\]

In the first case the minimizing choice is `a = theta/(sqrt(2)c)`. Starting from `c0 = 1` with `theta < 1`, every iterate satisfies the first-case condition: `g(c) >= sqrt(2theta) > sqrt(theta)`. Hence the claimed iteration is valid at every depth. Its positive fixed point is `c-star`, with derivative `sqrt(2)-1`. It decreases from 1 to that point when `theta < sqrt(2)-1`; it increases when `theta > sqrt(2)-1`.

At `theta = 3/10` the first values are

\[
 1,\ 0.919238816,\ 0.880769231,\ 0.863646494,\
 0.856313977,\ 0.853232352,\ldots\to0.851037055.
\]

For the requested fixed-degree infima:

| `d` | Infimum of `theta` | Exact `c-star` | Decimal `c-star` | Balanced heuristic NFS |
|---:|---:|---|---:|---:|
| 7 | `3/8` | `sqrt(3(1+sqrt(2))/8)` | 0.951488353 | 0.707106781 |
| 9 | `3/10` | `sqrt(3(1+sqrt(2))/10)` | 0.851037055 | 0.632455532 |
| 11 | `1/4` | `sqrt(1+sqrt(2))/2` | 0.776886987 | 0.577350269 |
| 13 | `3/14` | `sqrt(3(1+sqrt(2))/14)` | 0.719257588 | 0.534522484 |
| 19 | `3/20` | `sqrt(3(1+sqrt(2))/20)` | 0.601774073 | 0.447213595 |

These are infima as fixed `k` increases and positive slack tends to zero, not values attained by one finite parameter set at the displayed infimum.

At `d = 9, k = 13`, exact zero-slack values are

\[
 \nu=1/52,\quad \mu=17/156,\quad h_*=1/26,\quad
 \theta_0=16/39.
\]

With the required positive slack, `theta = 16/39 + eta0`. At zero slack,

\[
 a_* = \sqrt{\frac{8(\sqrt2-1)}{39}}
       =0.291490797\ldots,\qquad
 c_* =4\sqrt{\frac{1+\sqrt2}{39}}
       =0.995211832\ldots.
\]

Theorem A's corresponding value is `55/(39 sqrt(2)) = 0.997201871...`. Both constants remain below 1 precisely when

\[
 0<\eta_0<\sqrt2-1-16/39
            =0.00395715211668\ldots.
\]

The rational certificate `(55/39)^2 = 3025/1521 < 2`, with gap `17/1521`, confirms the strict zero-slack threshold without floating-point comparison.

The repairs required for these conditional conclusions are explicit and local: allow the input 1, canonicalize oracle outputs, include `a` in thresholds or remove those thresholds with the uniform verified fallback, state the general-`S` PZ argument, and give a computable portfolio sequence with vanishing slack. None changes `c(a)` or the portfolio limit. The strongest remaining mathematical vulnerability is inherited from the relation-supply proof, not hidden in the recursion algebra.

VERDICT RA: ACCEPT WITH FIXES (conditional on Theorem A at a = theta/sqrt2 being correct)
VERDICT P: ACCEPT WITH FIXES (same condition)
WEAKEST LINK: The inherited FM transfer from restricted smooth-modulus distribution to the weighted moving-residue family; a wrong relative-error or induced-weight estimate there would destroy the relation supply, whereas the recursive and restart arguments do not repair it.
