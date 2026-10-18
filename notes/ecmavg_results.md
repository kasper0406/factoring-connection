# ecmavg: an explicit ECM family as a prefilter for deterministic average-case factoring

Agent `ecmavg`, phase 6, 2026-10-18. Starting point: codex/r18_ideas.out §5. Inputs: Theorem AV
(notes/avgtime_results.md), the n⁴/H⁴ coverage bound (notes/sieve2d_results.md, Thm 2), Theorem R
(notes/dchworst_results.md §3) and T1/T2 (notes/below16_results.md).
- Pre-registration: `notes/ecmavg_prereg.md`, commit `4013d9e`, made before any run.
- Code: `src/ecmavg/{exponents,curves,labels,experiments,explore}.py`.
- Tests: `tests/test_ecmavg.py` (9 tests).
- Raw data: `experiments/raw/ecmavg_labels.jsonl`.
- Analyses: `experiments/analysis/ecmavg_{summary,explore,tail}.json`.
- Reviews: GPT-6.1-sol `codex/ecmavg_1.out` (blind + adversarial) and an Opus critic (§8).

Notation: S_n is the set of balanced pairs, uniform. N^x means N^{x+o(1)}. B = N^b is the stage-1
cutoff, J = N^θ the number of curves, and δ the density of a subset A ⊆ S_n.

## 0. Verdict (plain language)

**The route does not beat 1/6 unconditionally, and I could not prove any restricted version of
the moment hypothesis it needs. Everything deterministic in it is now proved, and it is sharper
than codex stated. The probabilistic hypothesis (EC2), however, is the wrong one: as stated it is
heuristically false for every polynomial number of curves. The hypothesis the route really needs
is a power-saving "hitting" theorem for a fixed, small ECM family. Its strength is close to a
deterministic average-case ECM theorem, and known techniques are far from it.**

- **Fallback, Theorem F (proved, given Theorem R, sieve2d Thm 2 and HH).** Use Theorem R's DCL-wc
  ladder, not AV's D′ ladder. Start it at H₀ = N^{θ/4}, cap it at N^{1/50}, then run HH. For
  *every* A ⊆ S_n of density δ = N^{−θ}:

  E[1_A·T_fallback] ≤ N^{o(1)}·min(δN^{1/5}, N^{1/6}δ^{7/12}).

  The worst case stays at HH. Codex's additive N^{2/15} term is an artefact of two choices: AV's
  ladder (its E_c and C(T) sets pay the worst case) and stopping at N^{1/60}. GPT-6.1-sol reached
  the same conclusion independently.
- **Re-optimisation (proved, conditional on the hypothesis).** 23/152 is the optimum of one
  particular cost model: B fixed at N^{1/8}, and every input pays for all J curves. It is not the
  optimum of the route.
  - The prefilter can instead run its curves in order and stop at the first divisor. Suppose the
    failure fraction after j curves is ≤ j^{−κ} for all j ≤ N^Θ. Then the exponent is
    max(b + (1−κ)₊Θ, F(κΘ)), where F(θ) = min(1/5 − θ, 1/6 − 7θ/12).
  - Under EC2 (κ = 1) up to J = N^{1/38}: 23/152 for every b ≤ 23/152.
  - Up to N^{1/14}: 1/8.
  - Up to N^{1/5−b}: **b, for any b > 0, so arbitrarily small exponents.**
  - Batch model with b → 0: 1/10, not 23/152.
  - Weaker moments, Var ≤ J^{2−η} or failure ≤ J^{−κ}: min(2(1−κ)/(12−5κ), (1−κ)/5).
  - Stage 2 changes constants only.
- **EC2 is heuristically false.**
  - **Proved:** Var Z ≥ Var(E[Z | shared data]).
  - **Heuristic:** the family-average exposure rate fluctuates with p (p mod small primes, the
    position of p in its shell) by a relative ≍ 1/log B. So Var Z ≳ J²/log²B, which exceeds
    J·N^{o(1)} once J ≫ log²B.
  - **Data (exploratory):** the mean cross-curve covariance is positive at z = 3–29 in 22 of 24
    cells. c̄/s²·log₂²B stays at 0.3–1.6 from b = 12 to 22, i.e. no decay beyond 1/log²B.
    Conditioning on p mod 840 and the size of p removes only part of it.
  - (Heuristic.) Chebyshev with the true variance therefore gives only the upper bound
    failure ≲ 1/log²B, which is no power saving and no exponent gain.
- **The correct hypothesis is HIT.** Prefix failure ≤ J^{−κ}. With staged gcds plus an
  all-but-one exponent step, a pair fails only if both primes are unexposed, or if both have
  *equal* point orders (Lemma S′, proved). At toy sizes the pair failure is ≈ P(W=0)².
- **What is provable (literature, §6).**
  - The transfer theorems go from all curves mod p to explicit boxes. They need family size
    AB > p^{1+ε} (Banks–Shparlinski) or A, B ≥ √m (David–Smith). That is J ≥ N^{1/2}, so running
    the family costs more than trial division.
  - Fixed-curve smoothness is open for non-CM curves. For CM curves it is known only under EH
    (Barbulescu–Jouve), and at each p the CM curves over Q take ≤ 32 values of #E(F_p), so they
    cannot supply a growing number of independent group orders.
  - GRH alone gives fixed-level divisibility, not smoothness.
  - **No restricted covariance theorem obtained.**
- **Experiments (pre-registered, toy sizes b ≤ 22).**
  - K1 is CONSISTENT with EC2: R_u ≤ 1.61 in every K1 cell.
  - K2 does not fire: the staged failure drops ≥ 10⁴× from J = 8 to J = 64.
  - Q2, Q3 (where defined) and Q5 hold.
  - Q1 failed (I predicted INCONCLUSIVE), Q4 failed in 1 of 10 defined cells (0.897 < 0.90),
    and Q6 failed (the excess decays by a factor 3.4, not < 3).
  - **No toy-size result bears on the asymptotic claim.**

| statement | label |
|---|---|
| Theorem F (fallback on any subset), Lemmas S, S′, Corollary D (failure decomposition) | proved (given cited theorems) |
| Theorem C (exponent under HIT), all optimisation values, the 2/15-floor diagnosis | proved (conditional implication; exact arithmetic) |
| EC2 ⇒ HIT with κ = 1 | proved (Chebyshev) |
| Var Z ≥ Var E[Z∣𝒞] | proved (trivial) |
| relative common-mode fluctuation ≍ 1/log B, so EC2 fails for J ≫ log²B | heuristic (standard local model) + exploratory data |
| unconditional or GRH proof of HIT for any explicit family with J ≤ N^{1/5} | not obtained; out of reach of the cited methods (§6) |
| r18's 23/152 as the route's exponent | refuted for this model (artefact of fixed b + batch cost) |

## 1. The algorithm ECMAV(b, Θ) (n-only parameters)

1. **Prefilter** (`src/ecmavg/curves.py`, operational: it sees only N).
   - Fixed explicit family: Suyama σ = 6, 7, …, Montgomery x-only.
   - For j = 0, 1, …, J_max − 1 with J_max = ⌈N^Θ⌉:
     - (step 0) gcd(D_j, N), where D_j is the degeneracy integer;
     - then the stage-1 prime powers l^{e_l} ≤ X (l ≤ B, in increasing order), each followed by
       gcd(Z, N).
   - Output any proper gcd and stop.
   - If gcd = N occurs at step i of curve j, run the **all-but-one refinement** (Lemma S′); if that
     gives nothing, go on to curve j + 1.
2. If no divisor was found: **fallback FB_{κΘ}** (Theorem F).
3. Finish as in AV: validate the divisor, run AKS on both cofactors, and run HH on anything else.
   - On S_n the finish never calls HH.

## 2. Theorem F (fallback on an arbitrary subset) — proved

**FB_θ** (θ ∈ [0, 1/5] n-only).
- If θ > 2/25, run HH.
- Otherwise run DCL-wc (Theorem R) with D′'s design at heights H_i = 2^i·H₀, where
  H₀ = ⌈n·2^{θn/4}⌉, for as long as H_i ≤ H_max = 2^{n/50}. Then run HH.

**Theorem F.** There are absolute n₀ and C such that for n ≥ n₀ the following hold.
- **(i)** FB_θ is complete.
- **(ii)** Its worst case is HH(N) + O(N^{1/5}polylog) = N^{1/5+o(1)}. Capping at
  H_max = N^{1/50}/n^A instead gives HH(N)·(1 + o(1)) and changes no exponent (Opus critic).
- **(iii)** For every A ⊆ S_n with |A| ≤ N^{−θ}|S_n|,

  (1/|S_n|)·Σ_{(p,q)∈A} cost(pq) ≤ N^{o(1)}·min(N^{1/5−θ}, N^{1/6−7θ/12}).

**Proof.**
- (i) is immediate: every step validates its output, and HH is complete.
- (ii): each rung costs ≤ C n^c N^{1/6}H_i^{5/3} (Theorem R(iii) with D′'s list cost, plus
  base selection N^{1/6+o(1)}). H_i ≤ N^{1/50} gives N^{1/6+1/30} = N^{1/5}, and there are O(n)
  rungs.
- (iii), case θ > 2/25: HH runs only on A, giving δN^{1/5}. Here 1/5 − θ < 1/6 − 7θ/12, so the
  min is this branch.
- (iii), case θ ≤ 2/25:
  - **Reach sets.** R_0 = A. For i ≥ 1, R_i ⊆ A ∩ FAIL_wc(H_{i−1}) ⊆ A ∩ U(H_{i−1}). Here U(H)
    is the set of pairs uncovered by every block, and FAIL_wc(H) ⊆ U(H) is Theorem R(ii).
  - sieve2d Thm 2, valid for 2 ≤ H ≤ N^{1/15}/C₀, gives
    |U(H)|/|S_n| ≤ C n⁴H^{−4} + C n⁶H^{1/3}N^{−1/6}.
  - No nesting and no independence are used.
  - **Rung 0** contributes δ·n^c N^{1/6}H₀^{5/3} = N^{1/6−θ+5θ/12+o(1)} = N^{1/6−7θ/12+o(1)}.
  - **Rung i ≥ 1** contributes ≤ n^c N^{1/6}H_i^{5/3}(Cn⁴H_{i−1}^{−4} + Cn⁶H_i^{1/3}N^{−1/6}).
    - The first part is ≤ C′n^{c+4}N^{1/6}H₀^{−7/3}2^{−7(i−1)/3}. It is geometric and sums to
      N^{1/6−7θ/12+o(1)}.
    - The second part is ≤ n^{c+6}H_max² = N^{1/25+o(1)}.
  - **HH** contributes N^{1/5}·min(δ, |U(H_max)|/|S_n|) ≤ N^{1/5−max(θ,2/25)+o(1)}
    = N^{3/25+o(1)}. This is ≤ N^{1/6−7θ/12} because θ ≤ 2/25.
  - **Base selection** contributes δ·N^{1/6+o(1)} per rung.
  - Summing these terms gives (iii). ∎

**Remarks.**
- (a) **Why not AV's ladder.** D′'s failures also include E_c(H,T) and the detector coincidences
  C(T), and these are not controlled by coverage. In AV's proof they pay W = N^{1/5} at fraction
  N^{−1/18}.
  - Optimising that accounting, with HH after the ladder, gives the floor
    z ≥ 3/40 + 7h/2 with h ≥ (1/5 − z)/4, i.e. **z = 2/15 at h = 1/60, t = 14/45**
    (`exponents.fallback_AV_floor`). This is exactly r18's second term.
  - Theorem R removes E_c and C(T), and the floor disappears.
- (b) Starting the ladder lower than H₀ only adds a geometric sum. So FB need not know δ, except
  for the HH branch (GPT-6.1-sol, item 1).
- (c) The closed form F(θ) = max(0, min(1/5 − θ, 1/6 − 7θ/12)) equals a brute-force minimum over
  all (h₀, h_max) ladders on a rational grid (`fallback_R_search`, tested at 40 random θ).

## 3. Separation: Lemmas S, S′ and the failure decomposition — proved

Labels (`src/ecmavg/labels.py`, diagnostic) are defined as follows.
- k_p(j) = 0 if p | D_j.
- k_p(j) = i if the index of P⁺(ord(P_j mod p)) among the primes ≤ B is i and the order is
  B-smooth. Then Z vanishes mod p first after step i.
- k_p(j) = INF otherwise.

**Lemma S.** Curve j yields a proper divisor of N = pq in the staged prefilter iff
k_p(j) ≠ k_q(j). The divisor is the prime with the smaller label.

*Proof.*
- Z_i ≡ 0 (mod p) iff [M_i]P_j = O mod p, where M_i = ∏_{i′≤i} l_{i′}^{e}. This is monotone in i.
  It holds iff ord | M_i, since the exponents e satisfy l^e ≤ X < l^{e+1} and ord ≤ X, i.e. iff
  every prime of the order is ≤ l_i.
- The x-only Montgomery ladder is correct whenever the difference point is neither O nor the
  2-torsion point (0,0) (Costello–Smith, "Montgomery curves and their arithmetic", 2017; here it is also checked against an
  affine group-law brute force on 12 curves × all 100 primes in [101, 700) × 2 families).
- During step 1 the difference is P_j, with x₀ ≢ 0 because p ∤ D_j. From step 2 on, it is a point
  of odd order ≠ O.
- Once Z ≡ 0, it stays 0.
- So gcd(Z_i, N) = p exactly for k_p(j) ≤ i < k_q(j). ∎ (`test_prefilter_matches_labels_on_pairs`:
  150 pairs.)

**Lemma S′ (all-but-one refinement).** Suppose curve j gives gcd = N at step i. Compute
y_l = [M_i/l^{e_l}]P_j for every l ≤ l_i. The all-but-one exponent product tree (as in dchworst
Lemma G) does this in Õ(B) ring operations. To keep x-only arithmetic valid (no xADD whose
difference is the 2-torsion point (0,0)), first compute y₂ = [M_i/2^{e₂}]P_j by one ladder from P_j,
and for odd l pre-multiply by 2^{e₂} with xDBL only before the odd-part tree (Opus critic fix). Then take gcd(Z([l^t]y_l), N) for t = 0, …, e_l.
This yields a proper divisor unless ord_p(P_j) = ord_q(P_j).

*Proof.* y_l has order l^{v_l(ord_p P_j)} mod p and l^{v_l(ord_q P_j)} mod q. The first t at
which Z vanishes is the l-adic valuation. If the valuations differ for some l, the gcd at the
smaller one is a single prime. Equal valuations for every l means equal orders. ∎

(Lemma S′ is analysed but not implemented in `curves.py`; the experiments use Lemma S.)

**Corollary D (failure decomposition).** Let I be the primes in [2^{(n−2)/2}, 2^{(n+1)/2}]. Let
η(J) be the fraction of p ∈ I that none of the first J curves exposes. Let
ν = max_{j<J, d} #{p ∈ I : ord(P_j mod p) = d}/|I|. Include in ν the class "p | D_j" (O(log|D_j|) primes, so
O(n log J/P) of I). Then the fraction of S_n failing prefilter + S′ is ≤ C(η(J)² + ν).

*Proof.*
- |{(p,q) ∈ S_n : p ∈ U, q ∈ V}| ≤ |U||V|, |S_n| ≥ c₀P²/n² and |I| ≤ CP/n.
- A failing pair has, on every curve, either both primes unexposed or equal orders.
- If p is exposed by some curve, let j* be the first such curve; then q has the same order on j*.
  Sum over (j*, d). ∎

**Coincidence bounds (proved, partial).**
- For d ≤ D: #{p > P/2 : ord_p(P_j) = d} ≤ (2ĥ(P_j)d² + C_j)/log(P/2). The primes divide the
  denominator of x([d]P_j).
- For d > P^{1/2+τ}: ≤ 2(4P/d + 1)·9√P. Each value m = #E(F_p) is shared by ≤ 9√P primes (the
  Hasse window), and E and its twist give the factor 2. So ν ≤ Cn·P^{−τ}.
- The first bound already gives ν ≲ n ĥ d²/P ≤ n P^{−τ} for d ≤ P^{(1−τ)/2}. So the uncontrolled
  window is P^{(1−τ)/2} ≲ d ≲ P^{1/2+τ}, orders near √P, and it is **not controlled**. (GPT-6.1-sol: "hitting alone is
  insufficient; equal orders give identical annihilation profiles".)
- So a route needs both a hitting statement and an equal-order coincidence bound.

## 4. Theorem C (exponents under a hitting hypothesis) — proved conditional implication

**HIT(b, κ, Θ).** Let B = 2^{⌊bn⌋}. For all n ≥ n₀ and all 1 ≤ J ≤ N^Θ, the fraction of S_n on
which ECMAV's prefilter with its first J curves outputs no divisor is ≤ N^{o(1)}·J^{−κ}. The
N^{o(1)} is uniform in J.

**Theorem C.** Under HIT(b, κ, Θ) with Θ + b ≤ 1/5, ECMAV(b, Θ) followed by FB_{κΘ} is complete
and deterministic. Its worst case is N^{1/5+o(1)} (O(HH) with the n^{−A} caps), and its average cost over S_n is N^{e+o(1)} with

  e(b, κ, Θ) = max(b + (1−κ)₊Θ, F(κΘ)), where F(θ) = min(1/5 − θ, 1/6 − 7θ/12).

*Proof.*
- **Prefilter.** Curve j costs B·N^{o(1)}, including S′. It is reached only by pairs failing the
  first j curves, so E[prefilter] ≤ B·N^{o(1)}Σ_{j≤J_max} min(1, j^{−κ}). This is
  B·N^{o(1)}·J_max^{(1−κ)₊}.
- **Fallback.** The failure set has density ≤ N^{−κΘ+o(1)}. Apply Theorem F. No independence
  between ECM failures and coverage failures is needed.
- **Worst case.** J_max·B ≤ N^{1/5}. ∎

**EC2 ⇒ HIT(b, 1, Θ)**, provided EC2 holds for every prefix (or a doubling grid of prefixes):
- Zs ≥ Zu, so P(Zs_J = 0) ≤ P(Zu_J = 0) ≤ Var Zu/(E Zu)² ≤ N^{o(1)}/(c²J).
- EC2 at a single terminal J gives no stopping-time gain (GPT-6.1-sol).

**Values** (`src/ecmavg/exponents.py`, exact Fractions; `tests/test_ecmavg.py`):

| model | exponent |
|---|---|
| r18: batch cost J·B, B = N^{1/8}, AV fallback | 23/152 = 0.15132 (optimum of that model, θ = 1/38) |
| batch, b → 0, AV fallback (codex's (*) taken literally) | 2/15 = 0.1333. With the reach factor min(δ, N^{−1/15}) on E_c ∪ C(T) and the HH-only branch restored, the AV ladder also gives 1/10 (Opus critic) |
| batch, b → 0, R fallback | **1/10** (θ = 1/10, HH branch); GPT-6.1-sol's 2/19 uses the ladder branch only |
| stopping, κ = 1, Θ = 1/38 | 23/152 for every b ≤ 23/152 |
| stopping, κ = 1, Θ = 1/14, b = 1/8 | 1/8 |
| stopping, κ = 1, Θ = 1/5 − b | **b** (any fixed b > 0) |
| stopping, 0 < κ < 1, Θ free, b → 0 | min(2(1−κ)/(12−5κ), (1−κ)/5); the HH branch wins iff κ > 2/5. κ = 1/4: 6/43; κ = 1/2: 1/10; κ = 3/4: 1/20 |
| Var ≤ J^{2−η} | κ = η in the previous row |
| κ = 0 (e.g. Chebyshev with the common-mode variance, §5) | 1/6: no gain |

- **Is 23/152 optimal under EC2?** No. It is optimal only for fixed b = 1/8 and batch cost.
- **What EC2 really implies.** With κ = 1 over all prefixes up to N^{1/5−b}, the exponent is b for
  every fixed b > 0, i.e. N^ε average for each ε.
  - Letting c = c(b) enter EC2 explicitly with B = L_N[1/2] would even give L[1/2] average.
  - So EC2 is as strong as a deterministic average-case ECM theorem. That alone marks it as far
    beyond present technique.
  - One N^{o(1)} algorithm would also need uniformity as b → 0 (GPT-6.1-sol).
- **Stage 2** (BSGS/FFT continuation to B₂ ≍ B²) costs Õ(B + √B₂) = B^{1+o(1)}. It enlarges the
  accepted orders to "B-smooth times one prime ≤ B₂", which changes c(b), not exponents.
  (Proved accounting; the stage-2 correlations would need their own hypothesis.)
- **Escape from T1/T2.** The prefilter is not a collision search with charged labels. It succeeds
  by smooth-order annihilation. T1/T2 concern degree-1 key/window designs and say nothing here,
  consistent with below16's open-question list, where the smooth-order (ECM) hitting theorem is named as open.

## 5. Why EC2 is the wrong hypothesis

**Prop V (proved, trivial).** For any partition 𝒞 of S_n,
Var Z ≥ Σ_C P(C)(E[Z∣C] − E Z)².
- Suppose two classes C₁, C₂ of mass ≥ μ each have |E[Z∣C₁] − E[Z∣C₂]| ≥ εJ.
- Then Var Z ≥ μ[(m₁−m)² + (m₂−m)²] ≥ μ(m₁−m₂)²/2 ≥ με²J²/2.
- If ε, μ ≥ (log N)^{−K}, EC2 fails for every J = N^θ with θ > 0 (GPT-6.1-sol item 3; μ condition added by the Opus critic).

**Heuristic V′ (standard local model: Lenstra / Chebotarev at small primes, Dickman).**
- Small primes. Conditional on p mod l, the distribution of v_l(#E(F_p)) over a family changes.
  - For random curves, Pr[l | #E] = l/(l²−1) or (l²−2)/((l−1)(l²−1)), according as
    p ≢ 1 or p ≡ 1 (mod l).
  - That shifts the B-smooth probability of the remaining cofactor by a relative
    (−ρ′/ρ)(u)·ΔE[v_l]·log l/log B. For l = 2, 3 this is ≍ 1/log B.
- Size. The position of p in its dyadic shell shifts u by 1/log₂B, another relative
  ≍ (−ρ′/ρ)/log₂B.
  - For a pair, this cancels to first order in the symmetric separation probability (GPT-6.1-sol),
    but the local effects do not cancel.
- Families with forced torsion (Suyama: 12 | #E) couple the curves further through p mod 8 and
  similar conditions.
- So Var Zu ≍ J·s′ + J²c̄ with c̄/s² ≍ 1/log²B.
- Chebyshev's upper bound on the failure is then only ≍ c̄/s² ≍ 1/log²B. That is a polylog, not a
  power, so this route's exponent stays at 1/6. (This is no lower bound on the failure.)
- **Second-moment information caps κ at 1 anyway** (Cantelli: a distribution with given mean and
  variance can have P(Z=0) = Var/(Var + (EZ)²)). With stopping, κ = 1 suffices, so the real issue
  is the common mode, not the order of the moment.

**Data** (exploratory, `explore.py`, 100 permutations per null).
- The cross-curve covariance c̄ of the exposure indicators is positive and significant.
  - Its z-score against a row-permutation null is 3.1–29 in 22 of 24 cells (S2/S3, β ∈ {1/3, 1/4}).
  - The exceptions are `mont` at b = 18 and b = 22 with β = 1/4, at z = 1.6 and 2.5.
- The normalised c̄/s²·log₂²B lies in [0.34, 1.56] from b = 12 to 22, with no trend. That is
  consistent with ≍ 1/log²B.
- Conditioning on p mod 120, p mod 840, 16 size buckets, or (p mod 120) × 4 size buckets removes
  part of it. Typically 30–70% remains for `suyama`, at within-class z up to 22. At b = 22,
  β = 1/3, the remainder is z = 0.5.
- Extrapolation (heuristic): take B = N^{1/8} at n = 1024 (log₂B = 128) and s ≈ ρ(4) ≈ 0.005.
  - The J² term overtakes the J term at J ≈ log₂²B/(0.7s) ≈ 2^{22}.
  - That is already below r18's J = N^{1/38} ≈ 2^{27}.
- **The lower tail behaves much better than Chebyshev.** P(W_J = 0) is within a factor 1.04–4.33 of
  the independent-curves model ∏(1 − s_j) at J = 64 (`ecmavg_tail.json`). Pair failure ≈ P(W=0)².
  - So at toy sizes the common mode spoils the variance, not the hitting.

**The correct hypotheses.**
- (a) HIT itself.
- (b) EC2_𝒞: a conditional second moment E[(Z − E[Z∣𝒞])²] ≤ J·N^{o(1)}, together with
  min_C E[Z∣C] ≥ cJ, for a σ-algebra 𝒞 of "local data". Then Chebyshev holds class by class, and
  HIT follows with κ = 1.
  - In the model, the residual common mode after conditioning on Frobenius at l ≤ L is
    ≈ Σ_{l>L} l^{−5}(log l/log B)² ≍ L^{−4}/log²B.
  - So 𝒞 must include all l ≲ J^{1/4}. That is a growing level, beyond any Chebotarev range.
- (c) The single-prime version: η(J) ≤ J^{−κ/2}, plus the equal-order coincidence bound of
  Corollary D.

## 6. The analytic attack: what is and is not provable

Pairwise covariance across curves j ≠ j′ over primes p is the crux. Tools, with exact ranges
(own search plus GPT-6.1-sol §4. The theorem numbers and ranges are as GPT-6.1-sol reported them;
I checked only the arXiv abstracts myself, so the ranges carry that caveat):

| input | what it gives | why it does not give HIT for J ≤ N^{1/5} = p^{2/5} |
|---|---|---|
| Chebotarev in Q(E_j[m], E_{j′}[m]) (unconditional for fixed m; GRH for m up to a small power) | covariance of *fixed-level divisibility* m ∣ #E_j, m ∣ #E_{j′} | smoothness is an all-levels event; truncating to fixed levels "leaves the central problem" (r18) |
| Banks–Shparlinski 2009, Thm 12 / Cor. 15 (per prime) and Thms 16–18 (horizontal) | transfer from all curves mod p to boxes \|a\|≤A, \|b\|≤B of any isomorphism-invariant property | needs min(A,B) > p^{1/4+ε} and AB > p^{1+ε} (horizontal: AB > x^{1+ε}), i.e. J ≥ N^{1/2+ε}: running the family costs more than trial division; no prefix structure |
| David–Smith 2013 (+ corrigendum), Martin–Pollack–Smith 2014 | box averages of #{p : #E(F_p) = m} | A, B ≥ √m (log m)^{…}; asymptotics under BDH; fixed m, not smoothness |
| David–Koukoulopoulos–Smith 2017, Thm 1.3 | vertical Sato–Tate in short intervals (full population 𝒞_p) | full population mod p, trace windows, not small boxes or smoothness weights |
| Matomäki–Radziwiłł 2016; MR II (arXiv 2007.04290) Cor. 1.1 | for all but O(X^{1−c(u)}) prime centres, the Hasse interval has a positive density of X^{1/u}-smooth integers (GPT-6.1-sol's deduction; UNVERIFIED by me) | with Lenstra's Prop. 2.8 this controls *uniformly random* curves mod p, not a fixed family of J ≪ p curves; the box transfer above is the missing step |
| Barbulescu–Jouve (Acta Arith. 2024) Thm 1.2 | ψ_E(x,y) ∼ ρ(u)x/log x for fixed CM E, under EH for quadratic fields; non-CM needs a stronger EDS hypothesis | conditional result is the asymptotic. Over Q there are 13 CM j-invariants, and twists give ≤ 6/4/2 orders, so at each p all CM curves over Q take ≤ 11·2 + 4 + 6 = 32 values of #E(F_p) (standard facts). For *group-order* smoothness a CM family therefore offers boundedly many chances per p, hence no power saving in J. For *point-order* smoothness this is only heuristic: a point order can be smooth when the group order is not |
| Lenstra 1987 | rigorous reduction (Prop. 2.8); optimal L[1/2,√2] only conjectural (Conj. 2.10) | randomized; random curves mod p |

- **Conclusion (labelled: not obtained; out of reach of these methods).** Every unconditional
  transfer cited above, from random curves mod p to an explicit family, needs J ≥ p^{1+ε}
  (David–Smith needs AB ≥ m^{3/2}, which is larger still).
- HIT asks for families of size p^{2θ} ≤ p^{2/5} with prefix control. That is below every
  threshold above by a power of p, both per prime (Burgess-type p^{1/4} per side) and averaged
  over p.
- A rigorous "almost-all with any power saving in J" for an explicit family in this range would
  be a new theorem about reductions of a sparse set of rational curves. Combined with the Corollary D
  coincidence bound, such a theorem with κ ≥ 1 over all prefixes up to N^{1/5−b} would give
  deterministic N^ε-average factoring (Theorem C). With κ < 1 it gives
  min(2(1−κ)/(12−5κ), (1−κ)/5). I found no route to it under GRH
  either: GRH controls Chebotarev at fixed or small levels, not smoothness.
- The only restricted statements I can prove are fixed-level divisibility covariances
  (standard Chebotarev) and the coincidence bounds of §3. Neither moves an exponent.

## 7. Experiments (pre-registered; run exactly as frozen)

- Commands: `.venv/bin/python -m src.ecmavg.experiments labels` (92 s), then
  `… analyze` (21 s).
- Families: `suyama` (σ = 6..69) and `mont` (A = 3..66, x₀ = 2), with J ≤ 64.
- B = ⌊2^{bβ}⌋ with b the prime bit length.
- Statistics over all unordered pairs in each sample: S2 is all primes at b = 12, 14, 16 (up to
  4.6M pairs); S3 is 2048 primes at b = 18, 20, 22.

**K1 (literal EC2).** R_u = [Var Zu_64/64]/[Var Zu_8/8] for `suyama`:

| cell | β = 1/3 | β = 1/4 |
|---|---|---|
| S2, b = 14 | 1.20 | 1.61 |
| S2, b = 16 | 1.21 | 1.45 |
| S3, b = 22 | 1.05 | 1.13 |

All are ≤ 2, so **CONSISTENT** with EC2 at toy sizes. The variance excess is real but small at
J ≤ 64. Pure J² growth would give 8.

**K2 (route at toy sizes).** fail_s(64)/fail_s(8) = 1.9·10⁻⁵/0.19 = 10⁻⁴ (S2, b = 16, β = 1/4).
**Does not fire.**

**Predictions.**
- **Q1 FAILED.** I predicted INCONCLUSIVE; the verdict is CONSISTENT. R_u ∈ [1.1, 4) held in 5 of
  6 cells.
- **Q2 HOLDS at J = 64** (the J used in the summary). c̄ > 0 in every S2/S3 cell, for both families
  and every β. The prereg did not fix J (disclosed): at J = 16, c̄ = −2.3·10⁻⁴ < 0 in S3 `suyama`,
  b = 22, β = 1/2. Q2 holds at J = 8, 32 and 64.
- **Q3 HOLDS where defined.** In 3 of the 24 S2/S3 cells with β ∈ {1/2, 1/3},
  fail_s(8) = fail_s(64) = 0, so the ratio is undefined. The "where defined" handling was not
  frozen (disclosed). Every defined ratio is ≤ 2·10⁻⁴.
- **Q4 FAILED in 1 of 10 defined cells.** S2 `suyama` b = 16, β = 1/4: 0.897 of staged failures
  are double-unexposed. The others are 0.909–1.000.
- **Q5 HOLDS.** `mont`'s R_u lies within ×2 of `suyama`'s in every K1 cell (0.94–1.09 vs
  1.05–1.61).
- **Q6 FAILED.** R_u − 1 for `suyama`, β = 1/3, fell from 0.17 (b = 12) to 0.05 (b = 22), a factor
  3.4. That exceeds my "< 3". The endpoint ratio matches the 1/log²B ratio (4/7.33)⁻² ≈ 3.4, but
  the series is non-monotone (0.165, 0.197, 0.214, 0.191, 0.124, 0.049), so no decay law is
  inferred.

**Other frozen statistics** (`ecmavg_summary.json`).
- **Mean.** E Zu/J is flat in J, as it must be (it is the average of per-curve rates). For `suyama`,
  β = 1/4, the separation rate E Zu/J goes 0.376 → 0.183 → 0.092 at b = 12 → 16 → 22. The
  exposure rates s go 0.251 → 0.102 → 0.048, falling toward the asymptotic ρ.
- **S1 (codex's exact experiment, 128 primes).** Var Zu/J rises with J for `suyama`, e.g. 0.249
  → 0.491 at b = 12, β = 1/4. For `mont` it is flat (0.10–0.29).
- **Identical signatures.** Pairs with identical unstaged 64-curve signatures occur at
  β = 1/4 for b ≥ 14 (≤ 2.7·10⁻²), and also in S3 `mont` b = 22, β = 1/3 (2.9·10⁻⁶). They are almost
  entirely all-zero signatures. Their fraction equals fail_u by definition, and ≈ P(W=0)².

**Not inferable.** Every quantity is pre-asymptotic: s is far above ρ(u) at these sizes. Neither
EC2 nor HIT is tested asymptotically, and finite data cannot "kill" the route (GPT-6.1-sol).
What the data do show is that the common mode exists, is not purely local, and is of the size
the heuristic predicts.

## 8. Review

- **GPT-6.1-sol `codex/ecmavg_1.out`** (blind + adversarial, xhigh). Independently it found:
  - Theorem F, with R rather than AV, and the min(δN^{1/5}, N^{1/6}δ^{7/12}) form;
  - the 2/15 term as an artefact;
  - the stopping formula and 1/8 at θ = 1/14;
  - that N^ε follows if EC2 holds up to 2/7 (that range breaks the worst-case condition Θ + b ≤ 1/5;
    superseded here by Θ = 1/5 − b);
  - the weak-κ formula (ladder branch only; the HH branch (1−κ)/5 for κ > 2/5 is added here);
  - the common-mode objection to EC2 (Prop V), and that hitting alone is insufficient because of
    equal orders (adopted in Corollary D);
  - the literature table.

  Its corrections are adopted: the sieve2d additive error term is kept, the ladder is capped at
  N^{1/50}, the HH branch is added, and the finite-data "kill" is disallowed.
- **Opus critic:** see the addendum below.

## 9. Reproduction

```
.venv/bin/python -m pytest -q tests/test_ecmavg.py
.venv/bin/python -m src.ecmavg.experiments labels && .venv/bin/python -m src.ecmavg.experiments analyze
.venv/bin/python -m src.ecmavg.explore
.venv/bin/python -c "from src.ecmavg.explore import lower_tail; lower_tail()"
.venv/bin/python -c "from src.ecmavg.exponents import *; print(codex_point(), optimise(0,1,False,F(1,5),'R'))"
```

## 10. Candidate card

- **Mechanism.** Explicit ECM family (Suyama), staged gcds and the all-but-one refinement, early
  stopping, then a DCL-wc ladder and HH.
- **Proved.**
  - Theorem F: fallback on any subset, N^{1/6}δ^{7/12} ∧ N^{1/5}δ.
  - Lemmas S and S′ (with the 2-adic fix), and Corollary D.
  - Theorem C: average N^{max(b+(1−κ)₊Θ, F(κΘ))} under HIT.
- **Refuted for this model.** 23/152 as the route's exponent. EC2 as a plausible hypothesis
  (heuristic + data).
- **Obstruction.**
  - The needed HIT is a power-saving hitting theorem for J ≤ p^{2/5} explicit curves.
  - All known transfer theorems start at J ≥ p^{1+ε}.
  - CM/EH routes have boundedly many group-order values per p.
  - Separation also needs an equal-order coincidence bound in the window P^{(1−τ)/2} ≲ d ≲ P^{1/2+τ}.
- **Open question.** For an explicit family (Suyama σ ≤ J), prove for all but J^{−κ} of primes
  p ~ X that some σ ≤ J has #E_σ(F_p) X^{1/u}-smooth, for any κ > 0 and any J ≤ X^{2/5}. Even a
  GRH-conditional version would be new to us.

## Addendum: Opus critic (independent, read-only)

**Verdict: ACCEPT WITH FIXES.** The critic re-ran the 9 tests, recomputed every §4 value exactly,
and checked every frozen verdict against the JSON files. No main exponent or verdict changed.

**Confirmed correct.**
- Theorem F: reach sets, sieve2d range, every exponent, and R's per-rung cost terms.
- The 2/15 derivation.
- Lemma S.
- Corollary D.
- Theorem C and the whole table, apart from one row.
- EC2 ⇒ HIT.
- Prop V and the Cantelli remark.
- K1, K2 and Q1–Q6 as evaluated.

**Fixes applied.**
1. Worst case: N^{1/5+o(1)}, or HH·(1+o(1)) with n^{−A} caps; "O(HH)" was wrong.
2. The batch/AV table row is now labelled "(*) taken literally". With reach factors restored, it
   is 1/10.
3. Lemma S′ needed a 2-adic ordering to avoid xADD with difference (0,0).
4. Corollary D: the degeneracy class is now included in ν.
5. The uncontrolled coincidence window is P^{(1−τ)/2} ≲ d ≲ P^{1/2+τ}. It was misstated as starting
   at P^{1/3}.
6. Chebyshev wording: the bound is an *upper* bound ≍ c̄/s². §0 is now tagged heuristic.
7. Prop V needs μ ≥ (log N)^{−K}.
8. CM: "≤ 32 values of #E(F_p) per p", not 32 functions. Point-order smoothness is heuristic.
9. "Any power saving in J ⇒ N^ε" is now restricted to κ ≥ 1 together with the coincidence bound.
10. "Every transfer cited" replaces "every unconditional transfer".
11. Barbulescu–Jouve's result is the asymptotic.
12. Q2: J = 64 is disclosed, with the J = 16 exception.
13. Q3: 3 of 24 cells are undefined.
14. Q6: the series is non-monotone.
15. Corrected ranges: z 3.1–29, within-class z up to 22, tail factor 1.04–4.33.
16. Separation rates vs exposure rates.
17. The identical-signature cell list.
18. GPT's 2/7 range is now noted as superseded.

**Code.** `src/ecmavg/exponents.fallback_AV` still encodes (*) literally, to reproduce r18. This is
documented here and not changed.
