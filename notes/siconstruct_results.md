# siconstruct: constructing smooth integers in short intervals near kN (attack on Prop SI)

Agent `siconstruct`, phase 20, 2026-10-21. Version 4: after the Opus proof/cost critic, the codex adversarial round and the Opus literature critic (§8).

**Files.**
- Code: `src/siconstruct/` — `core.py` (exact sieve, N-oblivious smooth lists, Hildebrand–Tenenbaum model),
  `constructors.py` (brute force, MITM, 2-D aim), `theory.py` (cost formulas; Prop CV inequalities),
  `run_exp.py` (E1–E3), `diag_e3.py` (post-hoc E3 diagnostic).
- Tests: `tests/test_siconstruct.py` (8 pass).
- Pre-registration: `notes/siconstruct_prereg.md`, commit 9d75362, made before every run.
- Raw data: `experiments/raw/siconstruct_E{1,2,3,3diag}.jsonl`.
- GPT-6.1-sol: blind round `codex/siconstruct_1.{md,out}` (§8).

**Labels.** PROVED (proof here, elementary unless said) · KNOWN (cited) · HEURISTIC (named assumption) ·
MODEL (exact consequence of named heuristics) · SCENARIO (stipulated algorithmic assumptions) · MEASURED ·
POST-HOC · assessment. "No prior statement found by us" is the strongest novelty claim made.

**Notation.** x = ln N, F_α = x^α (ln x)^{1−α}, L[α, c] = exp((c + o(1))F_α). Factor base y = L[1/3, κ], n = π(y);
window h = L[2/3, ν]; I_k = [kN, kN + h]; multipliers k ≤ K = L[1/3, C]; ε := h/(kN), the relative window.
Target: a y-smooth u ∈ I_k at cost L^{o(1)} per output (Prop SI ⇒ L[1/3, 1.2114]).

---

## 0. Verdict

**Plain language.**
- **We found no constructor.** Every approach we analysed costs at least L[2/3] per output (sieving-type methods),
  or exponentially much (meet-in-the-middle, k-list, knapsack algorithms), or far more (generic lattice CVP in
  dimension L[1/3]). None comes near L^{o(1)}, or even near the weaker L[1/3, 0.757] per output that would still
  beat MNFS. The barriers below are first-moment statements in named models (flat products, independently placed
  factors, Johnson-type certificates). They are not lower bounds against arbitrary adaptive algorithms; constructors
  that select factors *from the target residues* more cheaply than sieving remain OPEN.
- **Meet-in-the-middle does not amortise.** One pair of sorted lists does serve every multiplier k at once, but
  (assuming the products are flat modulo N) the probability that a product lands in *some* window stays (h+1)/N per
  pair, whatever the number of k served. So R
  outputs need |A||B| ≈ RN/h pairs and cost ≥ 2√(RN/h). Per output that is ≈ √(N/(Rh)) = N^{1/2−o(1)} for the
  L[1/3] outputs Prop SI needs. The per-output cost becomes L^{o(1)} only if Rh ≥ N^{1−o(1)}, i.e. with exponentially
  many outputs. Measured (E1): over K = 1…256 served multipliers the hit rate per pair stays at 0.92–1.09 of (h+1)/N.
- **Why everything fails: the problem is scale-invariant.** If u = d·j lies in [X, X + h], then d lies in a window
  of the *same relative width* h/X around X/j. Splitting u into parts never relaxes the precision. An independently
  chosen factor s can be placed cheaply only when s ≲ Kh and its partner is a *tested* cofactor ≥ N/(Kh) =
  N^{1−o(1)}; larger s pay a placement penalty s/(Kh) that outgrows the cofactor gain (Lemma RC, a cost statement in
  expectation). That is lattice sieving, and it gives L[2/3] per output. Factors selected from the target (partial
  sieving) escape the lemma but, correctly charged, still cost L[2/3].
- **The lattice view is exact** (its geometry is Micciancio's 2001 locally dense prime lattice). With the target at
  the centre of I_k, squarefree Prop SI is equivalent (up to
  constants in h) to promise CVP with factor 1 + O(1/ln N) in Adleman's prime lattice of dimension π(y) = L[1/3]
  (Prop CV(c*), repaired after codex found that the draft's left-end target returns u < kN). The instance sits at
  distance ≈ λ₁/√2 from the lattice when two witnesses exist. So Prop SI is near-exact CVP in dimension L[1/3]; every known *generic* CVP
  algorithm at that quality costs 2^{Θ(n)} = exp(L[1/3]). This is a reformulation, not hardness evidence: by Prop CV
  the interval sieve is itself an L[2/3] promise-CVP solver for this family.
- **One-sided is not easier than Schnorr's two-sided target in the benign-oracle cost model.** Its residue starts
  at the smooth-gap floor h, so it needs a stronger certified approximation to reach a given constant. Prop SI's
  1.2114 is schnorrcvp's Theorem G constant at slack τ = 2/(3κ) = 1.10 (both are suffmap Theorem SS(b) at residue
  log-size 1.1006·F_{2/3}). This compares cost models; it is not a reduction ordering.
- **Choosing k freely buys nothing at L-scale.** K = L[1/3] multipliers give K times more witnesses, which is
  L[2/3, o(1)]. We did not solve, or approach, the one-sided Bohr-set lister.
- **Johnson-radius / Boneh decoding** keeps exponent 2/3: inside the decoding range a call certifies O(1)·L^{o(1)}
  points, each with an untested cofactor ≈ N/h; sihard's Theorem J proves the averaged per-call bound
  L[2/3, −2/(3κ) + o(1)] smooth outputs for the whole Johnson class.

**Main results.**

| Result | Status |
|---|---|
| Prop CV: inequalities (a)–(e); centred-target exact identity (c*): squarefree Prop SI ⇔ promise CVP_{1+c/x} up to constants in h; dist ≈ λ₁/√2 given two witnesses; slack t certificate gives c(t) | (a), (b), (e) KNOWN in substance (Micciancio 2001 §5; Schnorr 2013); (c*) PROVED (codex), no prior statement found by us; cost reading MODEL |
| Lemma SR (scale invariance); Lemma RC: independently placed factor s, tested cofactor ⇒ per output ≳ max(1, s/(K(h+1)))/Pr[m smooth], optimum s ≈ Kh, L[2/3] | PROVED (SR, RC(i)); RC(ii)–(iii) HEURISTIC H_eq (placement flatness) |
| Lemma MA (two-list amortisation): flat products give \|A\|\|B\|(h+1)/N hits for any K; per output ≥ 2√(N/(R(h+1))) | HEURISTIC H_flat (sihard Prop MM proves an averaged-over-N form); MEASURED (E1, E2) |
| k-lists / Wagner / representations: no applicable fast construction established | MODEL (positive-product obstruction M^k ≤ KN; representation counts under a typical-shape assumption) |
| Per-output budget: L[1/3, τ] beats MNFS iff τ < 0.757; τ = κ (linear in n) gives 6^{1/3} = 1.817 | MODEL (exact) |
| Best toy constructors: MITM beats interval sieving at all tested sizes (32–80 bits), model crossover ≈ 300–310 bits (R₀ = 64) | MEASURED / MODEL |

---

## 1. What is needed, and a softer target

- Prop SI (suffmap §1.4): with a constructor of per-output cost L[1/3, τ] the factoring constant is
  c(τ) = min_κ max(κ + τ + 2/(9κ²), 2κ) (`theory.si_cost`, exact sympy).

| τ | 0 | 0.1 | 0.25 | 0.5 | 0.757 | κ (linear in n = π(y)) |
|---|---|---|---|---|---|---|
| c(τ) | (16/9)^{1/3} = 1.2114 | 1.2819 | 1.4031 | 1.6447 | 1.9017 (≈ MNFS; exact threshold 0.75719) | 6^{1/3} = 1.8171 |

- The 0.757 threshold re-optimises κ (at κ = (4/9)^{1/3}). At the Prop SI optimum κ = (2/9)^{1/3} held fixed, the
  threshold is τ < 1.9019 − (16/9)^{1/3} = 0.6905 (codex siconstruct_1).
- So the real target is weaker than L^{o(1)}: anything below L[1/3, 0.757] per output beats MNFS, and even a
  constructor costing O(π(y)) per output would give 1.817 (MODEL, linear algebra at ω = 2).
- None of the approaches below reaches L[1/3, anything] per output.

## 2. Structural results

### 2.1 Prop CV: Prop SI is near-exact CVP in dimension π(y) [(a), (b), (e) KNOWN in substance; (c), (c*), (d) PROVED here]

**Prior work (literature critic).** The geometry is Micciancio's *locally dense* prime lattice: Micciancio, "The
shortest vector problem is NP-hard to approximate to within some constant", FOCS 1998 / SICOMP 30(6) 2001, §5.
His Lemma 2 (every nonzero vector has ‖·‖_p^p > 2 ln α) is our (e) lower bound; his Lemma 3 (a squarefree smooth
g ∈ [β, β(1 + 1/α)] lies within (ln β + 2)^{1/p} of the target (0, …, α ln β)) is our (b); with α = β^{1−ε} (our
C = X/h) he notes many lattice points within ≈ λ/2^{1/p} of the target *provided* [β, β + β^ε] contains many
squarefree smooth numbers, and his Conjecture 1 (§6.2) is a squarefree-smooth-in-short-intervals hypothesis. The
λ₁ bound with the prime 2 is Schnorr, "Factoring integers by CVP algorithms", LNCS 8260 (2013), Lemma 2 (restated in
Oberwolfach Report 35/2014, pp. 1981–82). Liu–Lyubashevsky–Micciancio (APPROX-RANDOM 2006) use these lattices to
show BDD_α NP-hard for α > 1/√2 (checked only second-hand), so Prop SI sits exactly at that threshold; their
targets are planted (products of random prime subsets), so nothing transfers to targets ln(kN). **In Micciancio's
terms Prop SI asks to *find* a point of his locally dense configuration for the unplanted target ln(kN); his
reductions only need such points to exist.** Mulder–Sterner–van Woerden (arXiv 2509.17699, ANTS XVII; their
Def. 5.1 is the same lattice) observe that shortest vectors often correspond to smooth twins (the converse half of
(e)) and that the Gaussian-heuristic ratio is O(1/√log n) (their Prop. 6.4). New here, as far as we found: the
centred identity (c*), the multiplier/k reading and the cost bookkeeping (no prior statement found by us).


Lattice (Adleman weights, as in schnorrcvp §1): for e ∈ Z^n, b(e) = (√ℓ₁e₁, …, √ℓ_ne_n, C⟨e, ℓ⟩), ℓ_i = ln p_i,
target t_k = (0, …, 0, C ln(kN)), scale **C = kN/h**. Write u = ∏_{e_i>0}p_i^{e_i}, v = ∏_{e_i<0}p_i^{−e_i},
Q = Σe_i²ℓ_i, Δ = ⟨e, ℓ⟩ − ln(kN), so ‖b(e) − t_k‖² = Q + C²Δ².

- **(a) Lower bound.** For every e: ‖b(e) − t_k‖² ≥ ln(kN) − 1/(4C²).
  *Proof.* Q ≥ ln(uv) = ln u + ln v, and ln u = ln(kN) + Δ + ln v. So Q ≥ ln(kN) + Δ + 2 ln v ≥ ln(kN) + Δ, and
  min_Δ (Δ + C²Δ²) = −1/(4C²). ∎
- **(b) Witnesses are near-closest.** A squarefree y-smooth u ∈ I_k has Q = ln u ≤ ln(kN) + h/(kN) and 0 ≤ CΔ ≤ 1.
  So ‖b(e) − t_k‖² ≤ ln(kN) + 1 + h/(kN), and the output is a γ-CVP answer with γ² ≤ 1 + (1 + o(1))/ln(kN).
- **(c) Converse.** Suppose a witness exists (dist² ≤ ln(kN) + 2) and e is a γ-CVP answer with
  γ² ≤ 1 + c/ln(kN). Then ‖b(e) − t_k‖² ≤ ln(kN) + c′ with c′ = c + 2 + o(1), and:
  - 2 ln v ≤ c′ + |Δ|, so v = O(1);
  - |Δ| ≤ (√c′ + o(1))/C (from (a)'s chain: Δ + C²Δ² ≤ c′);
  - w := u − v·kN = v·kN(e^Δ − 1) satisfies |w| ≤ v·h·√c′·(1 + o(1)).
  So u ≡ w (mod N) with u y-smooth (by construction) and |w| = O(h), at multiplier vk; w may be negative (a
  two-sided window) and u need not be squarefree. *Proof:* the displayed inequalities follow from (a)'s chain
  Q + C²Δ² ≥ ln(kN) + Δ + 2 ln v + C²Δ² and Q ≥ ln(uv). ∎ (`test_prop_cv_converse_small_v` checks the first
  inequality and the identity for w at 80 digits on random {−1, 0, 1} vectors, which are far from closest; it is a
  sanity check, not a test of the near-closest regime.)
  - **This converse does not return a literal Prop SI output** (codex siconstruct_2, FATAL to the draft's "⇔"): the
    target sits at the *left end* of I_k, so u just below kN is as close as u just above. Codex's example: N = 143,
    y = 7, h = 67, k = 1; the witness 210 ∈ [143, 210] has D² = 6.0197, but exact CVP returns 105 < 143 (D² = 5.0886).
- **(c\*) Repaired converse with a centred target [PROVED; codex siconstruct_2, checked by us].** Put X = kN,
  m = (ln X + ln(X + h))/2, C = X/h, target C·ln T with ln T = m + 1/(2C²), and E_rep = Σ(e_i² − |e_i|)ℓ_i ≥ 0. Then
  for every e, exactly,
  D_T² = m + 1/(4C²) + 2 ln v + E_rep + C²(ln(u/v) − m)².
  *Proof:* Q = (ln u − ln v) + 2 ln v + E_rep; expand C²(δ − m − 1/(2C²))² with δ = ln(u/v). ∎
  - Any vector with v ≥ 2 or a repeated prime has excess ≥ 2 ln 2 = 1.386. A squarefree u ∈ I_k has excess
    ≤ (Cε/2)² ≤ 1/4 (ε = ln((X + h)/X) ≤ h/X).
  - So **exact CVP returns a squarefree y-smooth u ∈ I_k whenever one exists** (the squarefree smooth u closest to
    the centre in log scale). A γ-CVP answer with γ²(m + 1/2) ≤ m + 2 ln 2 has v = 1, is squarefree, and lies within
    O(h) of the centre of I_k (residue w = u − kN ∈ [−O(h), O(h)], which serves Prop SI equally well).
  - Hence, up to constants in h: **squarefree Prop SI (two-sided window) ⇔ promise CVP_{1+c/x}, c < 2 ln 2 − 1/2,
    on the centred instance.** Checked by exhaustive exact CVP on 6 toy instances, including codex's
    (`test_prop_cv_centred_exact_cvp_returns_interval_witness`).
- **(d) Slack (a sufficient certificate, not a necessity).** A CVP oracle with γ² − 1 = 2t·F_{2/3}/x returns
  v ≤ L[2/3, t] and |w| ≤ h·L[2/3, t]. The residue
  then has log-size (ν + t)F_{2/3}, and the factoring constant is
  c(t) = min_κ max(κ + 2/(9κ²) + t/(3κ), 2κ) (`theory.si_cost_slack`): 1.2114, 1.2481, 1.3937, 1.5702, 1.8982 at
  t = 0, 0.1, 0.5, 1, 2. In dimension terms (ln n = κF_{1/3}): F_{2/3}/x = (ln x/x)^{1/3} = κ ln x/ln n ≈
  3κ ln ln n/ln n. So in dimension n = L[1/3, κ] **this certificate gives an L[1/3] method when
  γ − 1 = O(ln ln n/ln n), and 1.2114 + ε when γ − 1 ≤ δ(ε)·ln ln n/ln n**, i.e. (γ² − 1)x = o(F_{2/3}). These are
  sufficient conditions from a residue *upper* bound. Whether a worse oracle must fail is the tightness question that
  schnorrcvp left open (its withdrawn Prop A); we do not claim necessity. Three requirements differ: bounded v
  (γ − 1 = O(1/x)), a literal one-sided output (c\*), and preserving the leading constant (t → 0).
- **(e) Geometry.** Every nonzero lattice vector has ‖b(e)‖² ≥ 2 ln C − ln ln C − O(1) = (2 − o(1))x.
  *Proof.* For a lattice vector write δ = ⟨e, ℓ⟩ = ln(u/v), so ‖b(e)‖² = Q + C²δ². If C|δ| ≤ √(2 ln C), then
  u ≠ v (unique factorisation), |δ| ≥ 1/(2M) with M = max(u, v), so M ≥ C/(2√(2 ln C)); and
  Q ≥ ln(uv) ≥ 2 ln M − |δ|. Otherwise C²δ² > 2 ln C. ∎ Conversely, *if* two squarefree witnesses u₁ ≠ u₂ lie in
  I_k, then e₁ − e₂ has ‖b(e₁ − e₂)‖² ≤ ln(u₁u₂) + 1 ≤ 2 ln(kN) + 2, so λ₁² = (2 + o(1))x. With (a)–(b): the target
  then lies at distance (1 + o(1))λ₁/√2. (The expected L[2/3, ν − 2/(3κ)] witnesses per interval is the local
  counting HEURISTIC, not a consequence of the inequalities.) The Gaussian heuristic, which predicts
  λ₁ ≈ √(n ln y/(2πe)) = L[1/3, κ/2], is off by an L[1/3] factor (as schnorrcvp §1 warned; Mulder–Sterner–van
  Woerden Prop. 6.4 make the same point).
- **Reading.**
  - With the centred target (c\*), squarefree Prop SI and promise CVP_{1+c/x} on this family are equivalent up to
    constants in h. With the left-end target (b)–(c), only a two-sided bounded-denominator proximity result holds.
  - These are **reformulations, not hardness evidence.** The cited bounds are *upper* bounds for generic algorithms
    on arbitrary lattices: exact CVP in 2^{n+o(n)} (Aggarwal–Dadush–Stephens-Davidowitz, FOCS 2015, arXiv
    1504.01995 Cor. 6.5, checked by codex), heuristic sieving 2^{0.29n} (Becker–Ducas–Gama–Laarhoven 2016 for SVP; Laarhoven 2016, Doulgerakis–Laarhoven–de Weger 2019 for CVP(P)); here exp(L[1/3, κ]). They are not lower
    bounds for this structured family, which the interval sieve already solves in L[2/3].
  - BKZ-β with Babai gives approximation factors ≥ 1 + Ω(1) for every β < n (SCENARIO: root-Hermite factor
    δ_β^{Θ(n)}; codex's warning in schnorrcvp that a root-Hermite factor is not a CVP ratio applies). By (d) that
    inflates the residue to N^{Θ(1)}: useless.
  - Per-call setup is not the obstacle: in a CVP-with-preprocessing formulation one lattice with C = N/h serves all
    k ≤ K (witnesses keep CΔ ≤ 1/k), at the price of an L[1/3] factor in the residue. The obstacle is the CVPP
    *quality* required (critic m3).
  - **Random thinning of the factor base does not help** (adaptive choice of a small prime support is not covered).
    Restricting to a random 1/r of the primes keeps a witness only if
    r^{u} ≤ L[2/3, ν − 2/(3κ)], u = x/ln y, i.e. r ≤ x^{κ(ν − 2/(3κ))} = poly(x) (HEURISTIC count; uses
    F_{1/3}F_{2/3} = x ln x). So the dimension stays L[1/3, κ].
- **One-sided versus two-sided (a comparison of two benign-oracle cost models, not a reduction ordering).**
  schnorrcvp's Theorem G: a benign two-sided oracle at slack
  γ² − 1 = 2τF_{2/3}/x costs L[1/3, 2√(τ/3)]. Prop SI's 1.2114 is that formula at τ = 2/(3κ) = 1.1006:
  2√(1.1006/3) = 1.2114 (by construction: both are suffmap Theorem SS(b) at residue log-size 1.1006·F_{2/3}).
  One-sided needs t → 0 to reach the same constant: its residue starts at the smooth-gap floor h = L[2/3, 2/(3κ)]
  instead of B. So in this comparison the one-sided model needs a stronger certified approximation quality to match
  a given constant. It is not "the same oracle": the settings differ in witness promises, scales C, denominator
  allowances and possibly dimension; a two-sided proximity oracle may keep returning points below the interval, and
  a doubly-smooth oracle also imposes smoothness of the residue (codex siconstruct_2 #11). Two-sided exact CVP heuristically gives polynomial time in polylog
  dimension (schnorrcvp Theorem H, HEURISTIC). Caveats: "harder" is a statement about the oracle quality required on
  lattices with different scales C, not a computational comparison; and the one-sided instance has L[2/3] solutions,
  which the interval sieve exploits at L[2/3] (no generic CVP algorithm does as well).

### 2.2 Lemma SR (scale invariance) and Lemma RC (placed factor, tested cofactor) [SR, RC(i) PROVED; RC(ii)–(iv) HEURISTIC H_eq / MODEL]

- **Lemma SR.** If u = d·j ∈ [X, X + h], then d ∈ [X/j, (X + h)/j], an interval of relative width h/X containing at
  most ⌊h/j⌋ + 1 integers (at most one once j > h, since its length is then below one). ∎ Every factor of u faces
  the same relative precision ε = h/X. SR is an identity: it explains why splitting does not relax the precision; it
  is not a lower bound against recursive algorithms.
- **Lemma RC** is the 1-D special-q bound: a placed factor s ≲ Kh buys only the smoothness of a cofactor smaller by
  s (special-q in QS: Davis–Holdridge; MPQS: Silverman 1987; Pollard 1993, "The lattice sieve"). Hittmeir's
  subset-sum and Smooth Subsum Search formulations (arXiv 2205.10074; arXiv 2301.10529, IJFCS 35(8) 2024) construct
  divisors in the same class. Fix a constructed integer s (smooth, a square, anything) and K. Let
  T_s = {(k, w) : 1 ≤ k ≤ K, 0 ≤ w ≤ h, s | kN + w}.
  - (i) [PROVED] |T_s| ≤ K(⌊h/s⌋ + 1), and every (k, w) ∈ T_s has tested cofactor m = (kN + w)/s ≥ kN/s.
  - (ii) [HEURISTIC H_eq: s generated independently of N mod s] E|T_s| = K(h + 1)/s·(1 + o(1)). This is a statement
    in expectation: rare placements with s > K(h + 1) happen (codex: N = 143, K = h = 1, s = 16, 144 = 16·9), but they
    are paid for. In a dyadic range of independently generated s the per-output cost is
    ≳ max{1, s/(K(h + 1))}/Pr[m is y-smooth], m ≈ kN/s (codex siconstruct_2 #5).
  - (iii) [MODEL] At s = Kh·L[2/3, σ] the placement penalty is L[2/3, σ] while the cofactor's smoothness improves only
    by L[1/3, O(1)], so σ = 0 is optimal. There the cofactor has log-size x − νF_{2/3} − O(F_{1/3}), and the ratio of
    its smoothness probability to that of a random integer ≍ N is ρ((x − νF_{2/3})/ln y)/ρ(x/ln y) = L[1/3, 2ν/(3κ)]
    (derivative of u ln u is (2/3 + o(1)) ln x). Per output: L[2/3, 2/(3κ)] up to that L[1/3] ratio. (The hidden
    error in a standalone L[2/3] estimate can exceed the whole F_{1/3} term; only the ratio is meaningful.)
  - (iv) [MODEL] As a direct relation generator (cofactor and residue both tested, total random size N^{1−o(1)}),
    Theorem SS of suffmap, re-optimised as a sampling program, gives at best L[1/2, √2].
- **Scope (narrowed after codex siconstruct_2 #10).**

| Route | Status |
|---|---|
| Continued fractions / Farey / Stern–Brocot approximations of kN/s, special multipliers k ≡ −wN^{−1} (mod s), base-m expansion, u = s²m with s² the placed part | reduce to RC with a larger multiplier range (CF choices try to violate H_eq, but any placement is still a lattice point of determinant s) |
| Several independently constructed factors; unequal splits | combine them into s; charge placement and generation; MA's AM–GM already allows unequal list sizes |
| **Adaptive partial sieving** (extract the smooth part of kN + w over half the factor base, keep those with a divisor d ≍ √N) | **escapes RC** (d depends on kN + w and the cofactor is ≍ √N), but correctly charged it costs L[2/3, 1/(3κ)] to find d and another L[2/3, 1/(3κ)] for the cofactor: total L[2/3, 2/(3κ)], the sieve |
| Factors chosen jointly from target residues more cheaply than sieving | OPEN (this is the problem) |
| Boneh / Coppersmith certified divisors | not covered by RC(ii) (the decoder selects d from the residues); see §4.4 and sihard Theorem J |
| Bivariate / whole-box Coppersmith | not covered; would need its own certificate and output analysis |
| Target-dependent norm forms or ideal selection | not covered; no cheap constructor known; ordinary NFS relations do not meet SI's height bounds |

- **Measured (E3, post-hoc diagnostic).** Successful cofactor tests roughly match the H-T density *at the cofactor's
  size* (pooled 0.95, records 0.49–1.25; the 56-bit small-s cell is low, pooled 0.64). This is consistent with the
  placement giving the smoothness gain of a smaller random cofactor and nothing more.

### 2.3 Lemma MA (amortised two-list merging) [HEURISTIC H_flat; MEASURED]

- *Setting.* Lists A, B of y-smooth integers; one merge serves every k: all pairs with a·b ∈ ∪_k [kN, kN + h].
- *H_flat (product flatness).* For the joint distribution of (A, B, N) considered,
  p := Pr[ab ∈ ∪_{k≤K} [kN, kN + h]] ≈ (h + 1)/N per pair, for any product range [N, KN]. Independence from N
  alone does not imply this: a structured multiplicative distribution need not be uniform mod N, and the density of
  the windows among integers is not their density among products of arbitrary lists (codex siconstruct_2 #7). sihard
  Prop MM proves the averaged-over-N (and over window position) form for fixed lists.
- *Consequence (lower bound).* R distinct outputs need |A||B|·p ≥ R, so the cost, at least |A| + |B|, satisfies
  (|A| + |B|)/R ≥ 2/√(Rp) = **2√(N/(R(h + 1)))** per output. No assumption on the matching time is needed, and
  repeated factorisations of one u (distinct outputs < hits) only make it worse.
  - L^{o(1)} per output iff R(h + 1) ≥ N·L^{−o(1)}.
  - Prop SI uses R = L[1/3, κ + ν/(3κ)] and h = L[2/3, ν], so Rh = N^{o(1)}: per output N^{1/2 − o(1)}; total
    ≥ √(RN/h), exponential.
  - The amortisation over k is real (one list, all k) but worthless: what is amortised is the *search*, while the
    *supply* of hits is fixed at (h + 1)/N per pair.
- *Matching cost (upper bound).* Our implementation costs Õ(|A| + |B| + K|A| + hits) (bisection per (a, k)). The
  draft claimed K-independence when Kh ≪ a; that is withdrawn: Kh < a only makes the candidates (k, b) for a given a
  collinear (a(k₁b₂ − k₂b₁) = k₁w₂ − k₂w₁ = 0), and a resonant ray can carry many of them (codex: N = ab₀ − 1 gives
  b = jb₀ for all j ≤ min(K, h)). The lower bound above does not depend on this.
- *Distinct outputs.* E1: 2453 distinct of 2901 hits; codex's toy run: 264 of 4381.
- *Prior art.* The two-list walk itself is Bernstein's exact counting/enumeration method (1995 "epsi" draft,
  Algorithm 3 / Lemma 3.2: a prime-power partition T ⊔ U, cost < 3(Ψ(x, T) + Ψ(x, U)), with unique factorisation and
  no multiplicity). The flatness lower bound for windows near kN (and sihard Prop MM) is the new part. The known
  N-dependent escapes from flat products, Schroeppel's linear sieve ((H + a)(H + b) − N small, H = ⌈√N⌉) and the
  Coppersmith–Odlyzko–Schroeppel cubic sieve, *test* their factors and leave residues N^{1/2} or N^{1/3}: they fall
  under RC / Theorem SS, not MA. They are analogues, not this note in disguise (suffmap §7.2).
- *Precomputation.* For an explicitly stored list B scanned against fresh A, per-N cost |A| ≥ RN/(h|B|), so per-N
  L^{o(1)} needs |B| ≈ N/h = N^{1−o(1)} — exponential, far larger even than a √N prime table. This is not a lower
  bound for every compressed data structure; a universal preprocessing/query lower bound is unproved.

### 2.4 k lists, Wagner, Schroeppel–Shamir, representations [MODEL]

- Hits among k lists number ∏|L_i|·p. For low-solution knapsacks there is **no generic Q^{1/2} lower bound**:
  representation algorithms (HGJ, BCJ 2^{0.291n}; Bonnetain et al., ASIACRYPT 2020, arXiv 2002.05276, 2^{0.283n}) beat it on their random
  subset-sum distributions. Their transfer to this sparse, correlated prime-log instance is unproved; applied to a
  binary encoding of all π(y) exponents they stay exponential in L[1/3]. We claim only that no applicable fast
  construction is established.
- **Positive-product obstruction for Wagner** (codex siconstruct_1, smoothenum Prop W): if every combination of k
  leaf factors must stay below KN and each leaf list has M distinct positive integers, then M^k ≤ KN, while
  Wagner's k-tree with k = 2^d needs M ≈ (N/h)^{1/(d+1)}; for k ≥ 4 these conflict. Oversized intermediate products
  need their own filtering analysis (sihard's bit-budget bound: cost ≥ √X/h in Wagner's filter semantics).
- **Representations.** The number of ordered k-factorisations of z is d_k(z) = ∏_{p^e∥z} C(e + k − 1, k − 1). For a
  *typical-shape* y-smooth z in the window (squarefree part dominant, Ω(z) ≈ x/ln y, HEURISTIC) this is
  k^{(1+o(1))x/ln y} = L[2/3, o(1)] for fixed k: N^{o(1)}, too few to remove an N^{Θ(1)} cost. It is not uniform in
  k: for z = 2^{⌊x/ln 2⌋} and k ≈ x, d_k(z) = N^{1.65+o(1)} (codex), but such z are not in typical windows.

## 3. Per-approach cost table

Per-output cost of producing one y-smooth u ∈ I_k (any k ≤ K), y = L[1/3, κ], h = L[2/3, ν], R = L[1/3] outputs.
Costs are for the named strategy in its stated model; the sieve is always available as a fallback. "Total as a
direct relation generator" re-optimises the sampling program under Theorem SS (it does not keep Prop SI's (κ, ν)).

| # | Approach | Per-output cost (heuristic) | L^{o(1)}? | Total cost (via Prop SI / as a direct relation generator) | Status |
|---|---|---|---|---|---|
| 0 | Interval sieve (baseline) | 1/ρ(x/ln y) = L[2/3, 2/(3κ)] | no | L[2/3] / L[1/2, √2] (u and v both tested; Theorem SS, λ ≈ x) | MODEL; MEASURED toy |
| 1 | Log-knapsack / CVP in Adleman's lattice | near-exact promise CVP in dim L[1/3]: generic algorithms 2^{Θ(n)} = exp(L[1/3, κ]); BKZ-β, β < n: residue N^{Θ(1)} | no | worse than trial division (generic algorithms) | Prop CV PROVED; costs SCENARIO (upper bounds of named algorithms) |
| 2 | MITM, sorted logs, one list for all k | ≥ 2√(N/(R(h+1))) = N^{1/2−o(1)} | no | ≥ √(RN/h): exponential | Lemma MA (H_flat); MEASURED E1/E2 |
| 2′ | k-list / k-tree / representations | no applicable construction below N^{Θ(1)} known; Wagner blocked by M^k ≤ KN | no | exponential (known methods) | MODEL |
| 3 | Free k (joint (k, u) search) | factor K = L[2/3, o(1)] more witnesses; no change in per-pair or per-test rates | no | — | MODEL |
| 4 | CRT / Boneh / Coppersmith / GSS | O(1)·L^{o(1)} certified points per call, cofactor ≈ N/h: L[2/3, 2/(3κ)]·L[1/3, −Θ(1)] | no | L[2/3] / L[1/2, √2] | MODEL; sihard Theorem J (PROVED, averaged) |
| 5a | 2-D aim: CF / Farey / Stern–Brocot / special k / base-m / u = s·m; adaptive partial sieving | L[2/3, 2/(3κ)] (up to an L[1/3] ratio) | no | L[2/3] / L[1/2, √2] | Lemma RC (expectation, H_eq); partial sieving charged directly; MEASURED E3 |
| 5b | Oblivious smooth families: products of consecutive integers, binomials, factorial ratios, norm forms of products of small primes | hit probability (h+1)/N per member ⇒ N/h per output; with MITM as row 2 | no | exponential | HEURISTIC H_flat |
| 5b′ | Target-dependent norm forms / ideal selection; bivariate Coppersmith; structured k (Pell-type) | not analysed beyond: standard approximations give residues ≈ √N | unknown | — | OPEN |
| 5c | Square × smooth (QS-type aim, u = s²m) | aim precision ≈ √N: √N/h per output | no | QS itself is the m = 1, residue √N case | MODEL |
| 5d | Exact interval counting + binary search; global samplers; log-grid DP | counting: O(log h) oracle calls (+ polynomially many for the factorisation), but no cheap exact (or zero-distinguishing multiplicative) interval count known; subtracting global estimates needs relative precision h/X; samplers ≈ N/h; explicit log-grid DP ≥ N^{1−o(1)} bins | only with a new counting algorithm | — | reduction PROVED; rest MODEL (§4.6) |

## 4. Approaches in detail

### 4.1 Log-knapsack (approach 1)
- u = ∏p_i^{e_i} ∈ I_k ⇔ Σe_i ln p_i ∈ [ln kN, ln kN + ε], ε ≈ h/(kN) = e^{−x + O(F_{2/3})}.
- Exponent vectors: ‖e‖₁ ≤ x/ln 2; a typical witness has ≈ x/ln y = x^{2/3}(ln x)^{−2/3}/κ nonzero entries, mostly 1
  (a constant fraction of witnesses is squarefree, HEURISTIC).
- As a subset-sum instance: n = L[1/3] items, weight ≈ u = x^{2/3−o(1)} (sublinear), x bits of target precision.
  Candidates Ψ(N, y) ≈ N·L[2/3, −2/(3κ)], so the solution count is only L[2/3, ν − 2/(3κ)] = 2^{o(x)}: a
  **density-≈1 instance**, where every generic algorithm (MITM, Schroeppel–Shamir, HGJ/BCJ, Wagner) costs
  2^{Θ(x)} = N^{Θ(1)}. The "high density" n ≫ x is illusory: the bound Σe_i ln p_i ≤ x leaves only N^{1−o(1)}
  candidates.
- Lattice: Prop CV. Exact CVP is needed, in dimension L[1/3].

### 4.2 Meet-in-the-middle (approach 2): the amortisation analysis
- Lemma MA. The per-pair hit probability (h + 1)/N is fixed by the target set, not by how many k one list serves.
  Serving K multipliers enlarges the target set K-fold but also requires products K-fold larger; the density of the
  windows among integers is (h + 1)/N at every scale.
- **Precise amortised statement.** For R outputs the optimum is |A| = |B| = √(RN/(h+1)), cost 2√(RN/(h+1)) + R:
  - per output: 2√(N/(R(h+1))) + 1;
  - L^{o(1)} per output ⇔ R ≥ N/(h·L^{o(1)}), i.e. exponentially many outputs (unless h = N^{1−o(1)}, where the
    residue v is no longer small and the method is Dixon-like, L[1/2, √2]);
  - Prop SI needs only R = L[1/3, κ + ν/(3κ)]: the amortisation saves a factor √R = L[1/3] on N^{1/2}.
- **Factory variant** (lists precomputed once for all N of a size): in the scheme "precompute B, generate A per N",
  the per-N cost is |A| ≥ RN/(h|B|). In general any explicitly stored precomputed set (or sorted product table) needs ≥ RN/h
  elements for per-N L^{o(1)}: N^{1−o(1)} precomputation, far larger even than a √N prime table. Irrelevant
  (assessment).
- **Toy scale is misleading in MITM's favour.** On the toy L-path log₂h/bits runs from 0.71 (32 bits) to 0.57
  (80 bits), so N/h is tiny. MITM beats the
  interval sieve at every tested size (E2), and the o(1)-free model crossover with R₀ = 64 (a model extrapolation)
  is near 300–310 bits
  (310 with the cost 2√(N/(R₀h)), 300 with E2's measured cost model 5M/R₀). Asymptotically the
  per-output slope of MITM tends to 1/2 bit per bit of N, while the sieve's tends to 0.

### 4.3 Free choice of k (approach 3)
- The union of K windows has K·h·L[2/3, −2/(3κ)] smooth elements: K = L[1/3, C] is L[2/3, o(1)], so the existence
  threshold ν > 2/(3κ) is unchanged.
- Per-pair (MITM) and per-test (sieve, 2-D aim) success rates do not depend on K. A multi-target CVP gains at most a
  factor K. So free k is worth L[1/3] at most, on top of costs of size L[2/3] or N^{Θ(1)}.
- **Bohr-set check.** None of our methods lists u with u mod N ≤ L[2α], α < 1/3, or ≤ y. The one-sided Bohr-set
  lister (suffmap E2) at α = 1/3 is just Prop SI with free k, and we do not solve it. The hit law (h + 1)/N per
  oblivious pair shows the cost of every listed object grows as h shrinks, as it must.

### 4.4 CRT / Boneh / Coppersmith (approach 4)
- Boneh (JCSS 64 (2002), Thm 3.1: |I| < ¼·T^{log T/log S − 2.5/d}) finds all z ∈ [X, X + H] with gcd(z, P) ≥ S when
  ln H ≲ (ln S)²/ln P. His own §4 remarks that "the current bounds on generalized CRT list decoding are insufficient
  for improving the basic quadratic sieve" and that factoring needs a fully smooth integer, not a large smooth
  factor; the bullets below are a quantitative version of that remark at the Prop SI scale. Certified outputs have a
  known divisor d ≥ H (pilist Lemma S≥V, PROVED there). Here H = h, so d ≥ h. That bounds the untested cofactor from
  *above* only, and RC(ii) does not apply (the decoder selects d from the residues). A lower bound comes instead from
  the decoding range: d | P and ln P ≤ (ln S)²/ln h, so P, hence d, is small when S is (below), and the cofactor is
  ≥ N/P = N^{1−o(1)}. For the full factor-base primorial (ln P ≈ y) the decoding condition already fails badly:
  ln h ≲ x²/y = o(1) (codex siconstruct_1). Restricted moduli are analysed below; an independently charged count
  over the whole class (empty calls, modulus construction, target-dependent selection) is sihard's Theorem J.
- **Full smoothness is impossible inside the decoding range.** For d ≈ N one needs ln P ≤ x²/ln h ≈
  x^{4/3}(ln x)^{−1/3}/ν. If P's primes are near y, P involves O(x/ln x) primes and the window contains a T-smooth
  integer (T = those primes) with probability ≈ h·L[2/3, 1/(3κ)]/N; with small primes, O(x^{4/3}(ln x)^{−4/3})
  primes and probability ≈ h·N^{−3/4+o(1)} (critic). Exponentially small either way.
- **Best partial use** (S ≈ h^{1+δ}): the decoding range forces P ≲ h^{(1+δ)²}, and the expected number of certified
  points per window is ≈ h^{−δ}·#{m | P : m ≤ h^δ} = O(1)·L^{o(1)} (critic M2; our draft wrongly gave L[1/3] here).
  Each has an untested cofactor ≈ N/h, smooth with the RC(iii) gain L[1/3, 2ν/(3κ)]. Per output:
  L[2/3, 2/(3κ)]·L[1/3, −Θ(1)].
- **sihard Theorem J** (PROVED, averaged over the window position, any weights, batched over k) bounds the whole
  Johnson class by L[2/3, −2/(3κ) + o(1)] smooth certified outputs per call, attained by the trivial divisor trick.
- **Relation to Prop JR.** pilist Prop JR (L^{o(1)} batching gain) is proved in NFS scaling and explicitly excludes
  this congruence scaling, where one call gains g(d) = L[1/3, Θ(1)] (pilist critic M2). Our analysis is consistent
  with an L[1/3, Θ(1)] gain that leaves the exponent 2/3 untouched; sihard's Theorem J leaves the exact L[1/3] factor
  open.

### 4.5 Own ideas (approach 5)
- **Continued fractions / Farey / Stern–Brocot.** Convergents p/q of kN/s give s·p − kN·q small, a 2-D lattice of
  determinant s: Lemma RC.
- **Products of consecutive integers.** m(m+1)…(m+r−1) has smoothness probability ρ(u/r)^r = ρ(u)·r^u (better by
  r^u), but every factor is tested and consecutive values are spaced ≈ rN^{1−1/r} ≫ h. For fixed k there is at most
  one candidate, with probability ≈ h/(rN^{1−1/r}): exponential.
- **Binomial coefficients, factorial ratios** C(n, j), n ≤ y: smooth by construction, but an N-oblivious family.
  Products of them are a sub-semigroup of the y-smooth integers, so aiming them is the same knapsack (H_flat: hit
  probability (h + 1)/N per member).
- **Norm forms** (u = N(α), α a product of small primes of Z[i] or Z[√r]): the norm is multiplicative, so again the
  knapsack over split primes. Units do not change the norm.
- **Smooth-twin ladders.** Multiplying by ratios (j + 1)/j (j | u, j and j + 1 smooth) steps ln u by ≈ 1/j. For
  j ≤ y the step is ≥ 1/y, and precision e^{−x} needs exponentially many combinations (a knapsack with small items,
  MITM-bounded). For large j the step is small, but finding smooth (j, j + 1) is again the target problem; it is
  exactly the smooth-twin problem of Mulder–Sterner–van Woerden (arXiv 2509.17699), who find twins through short
  vectors of the same lattice.
- **Square × smooth.** u = s²m with s arbitrary is enough for the linear algebra (s² has even exponents), but it is
  Lemma RC with the constructed part s² ≤ Kh, or QS (m = 1, residue √N).
- **NFS reading.** Prop SI is two-polynomial "NFS" with f₁ = X, f₂ = X − N and a skewed box (a ∈ [kN, kN + h],
  b = k), with the f₁ side constructed. It is the degree-1 analogue of smoothenum's algebraic-side lister
  (infimum 0.7631 there).
- **Modulus-oblivious.** A constructor that sees only X = kN is a generic "smooth integer in [X, X + L_X[2/3]]"
  algorithm (sihard Lemma EQ). It would also give prime-field discrete logarithms in L_p[1/3, ≤ 1.442]
  (sihard Prop DLI, HEURISTIC). That is further evidence that it is hard (assessment).

### 4.6 Counting, sampling and dynamic programming (codex siconstruct_1, checked by us)
- **Prior positioning.** Granville, "Smooth numbers: computational number theory and beyond", MSRI 44 (2008), §5.7
  gives the interval-sieve cost ≪ z log log y + uy, conjectures that any algorithm needs ≫ π(y), and asks whether z
  can be removed from the running time "by a more constructive approach" (citing Boneh). That is the question of this
  note. His §5.8 calls exact Ψ(x, y) "likely impractical" and surveys approximations (Hunter–Sorenson 1997, factor
  1 + O(1/u + log y/y); Bernstein 2002, "Arbitrarily tight bounds on the distribution of smooth integers";
  Parsell–Sorenson 2006).
- **Counting ⇒ constructing** (standard self-reducibility). An exact oracle for Ψ(b, y) − Ψ(a − 1, y) gives a
  constructor by binary search: O(log h) calls find the least smooth integer in I_k, and the factorisation follows by
  bisecting on the prime cutoff (polynomially many further calls). So a cheap exact, or zero-distinguishing
  multiplicative, interval count would suffice.
- **No cheap count is known.** Subtracting two global estimates of size Xρ needs relative error below
  h/X = e^{−x+O(F_{2/3})}; saddle-point or Perron-integral estimates give nothing of the kind. Exact counting is
  lattice-point counting in a simplex of dimension π(y) = L[1/3] (Barvinok-type algorithms are polynomial only in
  fixed dimension); Bernstein's exact counting uses a meet-in-the-middle whose cost is governed by the two generated
  sets (codex: Bernstein, "Enumerating and counting smooth integers", §3, Algorithm 3). For an *explicit arbitrary*
  set of primes the existence question is NP-hard by sihard's X3C reduction (see sihard's note), so a fast counter
  would have to use that the factor base is an initial segment of the primes.
- **Sampling.** Uniform or approximate generators of factored smooth numbers (Bach; Kalai; Bach–Sorenson,
  arXiv 2006.07445 Thms 3.1, 5.1–5.3) sample globally. Conditioning on I_k costs ≈ N/h trials, and their
  enumeration order is lexicographic in the prime factors, not by size, so "select the rank near Ψ(kN, y)" does not
  select an integer near kN.
- **Log-grid dynamic programming.** This is Bernstein's 2002 rounding method (Granville §5.8). A DP over rounded
  logarithms with per-item error η and s items needs sη ≲ h/X, hence ≥ (ln X)/η = N^{1−o(1)} bins (codex; for that
  explicit-grid class only): its bounds sandwich Ψ between Ψ(x) and Ψ(x^{1+1/(M ln 2)}) for M bins, so resolving I_k
  needs M ≳ x·kN/h.
- **Existence theorems** are far from this scale. All-interval results (Harman 1991; Lenstra–Pila–Pomerance 1993;
  Jain, arXiv 2502.10530, Thm 1.2) need h ≥ x^{1/2+o(1)} and y ≥ exp((log x)^{2/3+o(1)}). Jain's almost-all Thm 1.1
  allows h = L[2/3] (h ≥ exp((1 + ε)((11/8)u log u + 4 log log X)), i.e. L[2/3, 11/(12κ)] at y = L[1/3, κ], the same
  shape as Prop SI's window) but needs y ≥ exp(C(log X)^{2/3}(log log X)^{4/3}), versus ln y ≍ (ln X)^{1/3} here
  (checked in arXiv 2502.10530v1 by the literature critic).

## 5. Total factoring cost and heuristic R

- Every constructor we found has per-output cost ≥ L[2/3]; plugged into Prop SI that gives L[2/3]. Used directly as
  relation generators (u and the residue both tested, total random size N^{1−o(1)}), they are at best Dixon-like,
  L[1/2, √2] (Theorem SS), which is worse than QS.
- MITM: total ≥ √(RN/h), exponential, on top of linear algebra L[1/3, 2κ].
- **Diversity (codex).** Multiplying one relation u₀ ≡ v₀ by a smooth s gives su₀ ≡ sv₀, a new valid output with the
  same exponent-difference vector: worthless. A constructor must produce L[1/3, κ] *independent* relations, and if one
  output per multiplier is used, the multiplier range must have C ≥ κ + ν/(3κ).
- **Fruitfulness (heuristic R).** Not tested here: no constructor reached the regime where it matters. MITM outputs
  are N-oblivious products whose residues are determined by (a·b mod N); sicheck's Prop B (Blum N: fruitfulness is a
  parity condition computable from N) and Cor. B (a constructor seeing only kN and N can make every dependency
  trivial) apply to any constructor. sicheck's Lemma I (gcd(u_k, v_k) = gcd(u_k, k)) explains the inherited duplicate
  relations behind our diversity caveat, and its measured fresh-relation rate (10–69% of naive H_S at 40–80 bits,
  recovering asymptotically) applies to any constructor whose outputs look like interval-sieve outputs. sicheck also
  notes that k ≍ N^θ changes the constant to (16(1+θ)/9)^{1/3}, so the multiplier range must stay N^{o(1)}.

## 6. Experiments (pre-registered, commit 9d75362)

Toy L-path (o(1) dropped): y = round(e^{κF_{1/3}}), h = round(e^{νF_{2/3}}), κ = (2/9)^{1/3}, ν = 1.376.
At 32/48/64/80 bits: y = 37/90/184/342, log₂h = 22.8/31.2/38.8/45.9.

### E1. MITM hit law and amortisation over k (48 bits, 3 N, K = 1…256, |A| = |B| = 4096)
- Pooled hits 2901 vs predicted 2888.2: **ratio 1.004** (P1a [0.85, 1.15]: PASS).
- Cells 0.776–1.189 (P1b [0.6, 1.5]: PASS).
- Per-K pooled ratios: K = 1: 0.944, 4: 1.044, 16: 1.085, 64: 0.924, 256: 1.025; max/min 1.174 (P1c ≤ 1.4: PASS).
  Serving 256 times more multipliers with one pair of lists leaves the per-pair hit rate unchanged.
- Distinct u: 2453 of 2901 hits (15% repeats through multiple factorisations).
- Here K places the products in [KN, 4KN] (about 3K multipliers served); it is not a literal cutoff k ≤ K.
- Cost counters (seed 0): queries per hit 41.9, 136.5, 498.4, 2492.7, 8621.8 for K = 1…256. At toy size h > √N, so
  each a sees ≈ 3K nonempty b-intervals. Whether the matching can be made K-independent in the Prop SI regime is
  open (§2.3); the lower bound does not need it.

### E2. Per-output cost along the toy path (32–80 bits, R₀ = 64)
| bits | log₂h | M = \|A\| = \|B\| | distinct outputs | MITM cost/output | model 5M/R₀ | brute force per output, seed 0: H-T at N / measured (k = 1..8) |
|---|---|---|---|---|---|---|
| 32 | 22.8 | 157, 166 | 44, 41 | 14, 15 | 12, 13 | 1.5e4 / 3.5e4 |
| 40 | 27.1 | 600, 587 | 52, 49 | 40, 42 | 47, 46 | 1.6e5 / 3.6e5 |
| 48 | 31.2 | 2550, 2312 | 56, 40 | 156, 199 | 199, 181 | 1.3e6 / 3.2e6 |
| 56 | 35.1 | 11101, 9441 | 44, 54 | 866, 600 | 867, 738 | 1.3e7 / 3.6e7 |
| 64 | 38.8 | 44076, 40820 | 52, 45 | 2897, 3100 | 3443, 3189 | 7.0e7 / — |
| 72 | 42.4 | 185928, 207848 | 55, 52 | 11537, 13636 | 14526, 16238 | 4.4e8 / — |
| 80 | 45.9 | 1052431, 879911 | 59, 45 | 61001, 66877 | 82221, 68743 | 3.6e9 / — |

- P2a: pooled hits/pred 822/897 = 0.916, cells 0.734–1.109: PASS.
- P2b: cost per distinct output / model 0.74–1.17 (factor 2 allowed): PASS.
- **P2c: FAIL as registered.** Measured sieve density / H-T density *at N* = 0.36–0.50. Cause: the registered
  comparison evaluated the model at N, but the sieved windows sit at kN, k = 1..8, where the density is up to 4×
  lower at toy sizes. POST-HOC, position-matched (mean of the H-T density at kN, k = 1..8; `analyze.py`): ratios
  0.76–1.08.
  The brute-force column above uses the H-T value at N for the comparison with MITM (whose products lie in [N, 4N]).
- P2d: MITM cheaper than brute force at every size: PASS (80 bits: 6.1e4–6.7e4 vs ≈ 3e9).
- Fitted slopes of log₂(per-output cost) per bit: MITM 0.254 (model 0.261), brute force (H-T) 0.365. At toy sizes
  h grows like N^{0.57+}, so the sieve's cost grows faster per bit. Asymptotically the slopes tend to 1/2 and 0:
  the model crossover (R₀ = 64) is ≈ 310 bits with the cost 2√(N/(R₀h)) (94.2 vs 93.0 log₂ per output at 320
  bits), ≈ 300 bits with the measured cost model 5M/R₀.

### E3. 2-D aim, constructed s and tested cofactor (40/48/56 bits, K = 16, 64 values of s, ≤ 2048 offsets per (s, k))
- P3b (Lemma RC(i)): every tested cofactor has bit length ≥ that of N/(2S0): PASS (all 12 cells, equality).
- **P3a: FAIL as registered for the small-s cells** (S0 = 2^{⌊lg h⌋/2}): distinct outputs / Σδ_HT(m), pooled over
  seeds per cell, = 0.405 (40 bits), 0.544 (48), 0.553 (56). PASS for the large-s cells (S0 = h/16): 0.879, 1.146,
  1.083.
  - Cause, partly (design flaw): with small s, the 64 values of s scan the same short stretches [kN, kN + 2048s],
    and one smooth u is found through several s (multiplicity 1.6–2.4 at 40–48 bits). At 56 bits the multiplicity is
    only 1.09–1.28, so it does not explain that cell.
  - POST-HOC diagnostic (`diag_e3.py`, same seeded lists): *successful tests with multiplicity* / Σδ_HT(m) = pooled
    5549/5815 = 0.954, records 0.49–1.25; the 56-bit small-s cell stays low (0.82 and 0.49, pooled 0.64). Tests are
    clustered (all s scan the same stretches near kN), so the effective sample size is much smaller than the test
    count. Reading: roughly consistent with cofactors being smooth like random integers of their size; not a sharp
    confirmation.
- **P3c: FAIL as registered.** Registered: tests per output ≈ 1/δ_HT(N/S0)·(1 ± 0.6). Measured
  tests per output × δ_HT(N/S0) = 1.64–2.49 (large s) and 2.95–5.90 (small s): all 12 records outside [0.4, 1.6].
  Cause: the same design error as P2c (the model is evaluated at N/S0, while the cofactors are ≈ kN/s with
  k ≤ 16, i.e. larger and less often smooth), plus the multiplicity above. The position-matched comparison is the
  diagnostic of P3a. The registered second clause (gain = δ(N/S0)/δ(N)) was not evaluated separately.
- Context (not a registered statistic): large-s cells need 33–140 tests per output, against 6.3e5–4.9e7 for the
  sieve (H-T model at KN/2). At toy sizes s ≈ h/16 leaves a cofactor of only 16–24 bits. That is the toy artefact
  h = N^{0.6+}. At Prop SI scale the cofactor is ≥ N/(Kh) = N^{1−o(1)} (on this path with K = L[1/3, 1]: ≈ 67% of
  the bits at 1024 bits, 79% at 4096).

**What the experiments show.** The laws behind the negative analysis hold at toy size: (h + 1)/N per oblivious pair
regardless of K (E1, E2, sharp), and smoothness of the tested part roughly at its size (E3, post-hoc, noisy). They do not and cannot show
asymptotic statements. MITM and 2-D aim are excellent at toy sizes only because toy h is a large power of N.

## 7. Deviations and post-hoc items
- P2c, P3a (small-s cells) and P3c failed as registered. P2c and P3c are design errors in the model *comparison*
  (wrong evaluation point); P3a is partly explained by multiplicity. Corrected comparisons are labelled POST-HOC.
- `diag_e3.py` and `analyze.py` (position-matched P2c, P3c numbers, pooled P3a) are unregistered, written after the
  runs.

## 8. Reviews
- **Codex blind round** (`codex/siconstruct_1.out`, GPT-6.1-sol, xhigh, run before our analysis was written down).
  It reached the same verdict independently: interval sieve L[2/3, 2/(3κ)]; construct-one-factor-and-test has a
  cofactor ≥ N/h (its version of Lemma RC); two explicit sorted lists cost N^{1/2+o(1)} per output for
  L[1/3, O(1)] outputs (its version of Lemma MA, with the exact sorted-list accounting Õ(m + n + Km + H));
  Wagner's tree violates M^r ≤ KN for r ≥ 4; Boneh's and GSS's recovery conditions fail by a huge factor; exact
  CVP costs 2^{m+o(m)} (ADS 2015, Cor. 6.5, checked by codex), m = L[1/3, κ]. Items we adopted:
  - a *binary* embedding e_p = Σ2^j z_{p,j} that enforces e ≥ 0, for which exact CVP solves Prop SI and
    γ² ≤ 1 + 3/(m + 1) suffices (an alternative to Prop CV that handles non-squarefree outputs; dimension
    m = O(π(y) ln x));
  - the fixed-κ threshold τ < 0.6905 (§1);
  - the diversity caveat (§5);
  - the counting-to-search reduction, the precision obstruction for approximate counts, the log-grid DP bound and the
    Jain and Bach–Sorenson references (§4.6).
  Codex labels general non-existence UNVERIFIED, as we do: these are barriers for named method classes.
- **Opus proof/cost critic** (read-only; recomputed all constants, proofs and data). Verdict: minor revision, no
  FATAL. MAJOR items, all fixed:
  - M1: P3c had no verdict and fails; now reported (§6, §7).
  - M2: §4.4 claimed L[1/3, 2ν/(3κ)] certified outputs per call (that is the cofactor gain, not a count; the correct
    count is O(1)·L^{o(1)}), and inferred a cofactor *lower* bound from d ≥ h; both repaired, sihard Theorem J cited.
  - M3: "1.2114 needs γ − 1 = O(1/x)" overstated the precision; corrected to o(ln ln n/ln n).
  Minor items adopted: converse of Prop CV holds with |w| = O(h) (the O(√x) was an artefact); quantifiers of Prop CV;
  the lattice view is a reformulation, not hardness evidence (the sieve is an L[2/3] promise-CVP solver); CVPP
  remark; caveats on "one-sided is harder"; RC(ii) in expectation with the σ trade-off; Lemma MA K-independence
  untested and Õ; sihard Prop MM; factory remark; k-list exponent (N/h)^{Ω(1)} and tuple count X(1 + ln X)^{k−1};
  twin ladders; split total-cost column; pooled P3a and the 56-bit cell; seed-0 brute column; crossover 300–310
  bits; the weak test of the converse is now described as a sanity check.
- **Codex adversarial round** (`codex/siconstruct_2.{md,out}`, GPT-6.1-sol, xhigh, on draft v1). Recommendation:
  major revision; no general constructor found. All items addressed in this version:
  - FATAL #1: the draft's "⇔" in Prop CV returned u < kN (explicit counterexample). Fixed by codex's centred-target
    identity, Prop CV(c*), checked by us and by an exhaustive test.
  - FATAL #5: Lemma RC turned an expectation into a pointwise claim. Restated as a per-output cost in expectation,
    with the placement penalty.
  - FATAL #6: the Boneh deduction had the inequality backwards and violated RC's independence hypothesis (also
    critic M2). Repaired through the decoding range; sihard Theorem J cited for the class-wide count.
  - MAJOR: §2.1(e) notation (δ vs Δ) and the two-witness promise; generic CVP runtimes are upper bounds, not
    barriers; slack conditions are sufficient, not necessary; H_flat instead of "N-oblivious"; the K-independent
    matching claim withdrawn (resonant rays); no generic Q^{1/2} bound for knapsacks, exact d_k(z) formula; scope
    table for RC (adaptive partial sieving escapes RC but costs L[2/3]); "one-sided is harder" reworded as a
    cost-model comparison; P3c disclosed (as the critic also found).
  - MINOR: threshold rounding; toy h range 0.57–0.71; crossover labelled extrapolation; E1's K semantics; the weak
    converse test now complemented by the exact centred test; literature additions (Schnorr's 2014 report on the
    prime lattice and λ₁, Granville 2008 §5.7, Bernstein's Algorithm 3), applied in v4 after the literature critic.
- **Opus literature critic** (read-only; checked Micciancio 2001, Schnorr's Oberwolfach report, Granville MSRI 44,
  Boneh JCSS, Jain v1, Bernstein epsi, ADS, Mulder–Sterner–van Woerden). Verdict: no correctness change, no published
  constructor, no prior statement of the 1.2114 reduction found. Applied:
  - M1: Prop CV (a), (b), (e) and "target at λ₁/√2" relabelled KNOWN in substance (Micciancio 2001 Lemmas 2–3,
    Schnorr 2013 Lemma 2, LLM 2006); Mulder–Sterner–van Woerden cited; new part = (c*), the k reading, bookkeeping.
  - M2: the literature additions promised after codex round 2 are now in the body.
  - M3: §4.6 positioned against Granville §5.7–5.8, Hunter–Sorenson, Bernstein 2002; counting ⇒ search called
    standard self-reducibility.
  - Minor: Jain's h-conditions; Boneh JCSS Thm 3.1 and his §4 remark; Bernstein's Algorithm 3 as the MA upper-bound
    mechanism; RC as special-q/MPQS accounting (Davis–Holdridge, Silverman, Pollard); Schroeppel and COS sieves as
    analogues; Hittmeir; lattice-sieving and knapsack citations.

## 9. Reproduction
```
.venv/bin/python -m pytest -q tests/test_siconstruct.py
.venv/bin/python -m src.siconstruct.run_exp E1   # ~1 min
.venv/bin/python -m src.siconstruct.run_exp E2   # ~10 min, < 1 GB
.venv/bin/python -m src.siconstruct.run_exp E3   # ~5 min
.venv/bin/python -m src.siconstruct.diag_e3      # ~5 min (post-hoc)
.venv/bin/python -m src.siconstruct.analyze      # registered + post-hoc statistics
```
