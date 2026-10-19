# detqs: deterministic quadratic sieve on almost all N (polyimpl D1)

Agent `detqs`, phase 8, 2026-10-19.
- Pre-registration: `notes/detqs_prereg.md`, commit `8a6b789`, made before any experiment run.
- Code: `src/detqs/` (operational: `sieve.py`, `linalg.py`, `qs.py`; model: `localmodel.py`;
  harness: `harness.py`, `analysis.py`). Tests: `tests/test_detqs.py`.
- Raw data: `experiments/raw/detqs_e{1,2,3}_n*.{jsonl,npz}`. Summary:
  `experiments/analysis/detqs_summary.json`.
- GPT-6.1-sol: blind round `codex/detqs_1.{md,out}` (scratch note
  `codex/scratch/detqs_1/fixed_u_moment_note.md`); adversarial round `codex/detqs_2.{md,out}`.
- Reviews: §10.
- Labels: **PROVED**; **PROVED IN OUTLINE** (every step is given, and the routine verifications still
  owed are listed); **CONDITIONAL**; **OPEN**; **experimentally supported**; **FALSE**. "Codex" means
  `codex/detqs_1.out`.

## 0. Verdict (plain language)

**No. We did not prove that a deterministic QS factors almost all balanced semiprimes below N^{1/6},
let alone subexponentially.** Of the two obligations, relation supply (H1) is now essentially in
hand at the polynomial level. Fruitfulness (H2) is untouched by every averaging argument we have and
remains the whole obstruction.

**H1 (relation supply) at fixed u: reduced to one open estimate, Proposition V.** The codex referee
rates V a GAP (§10), so H1-S is CONDITIONAL. What is proved in outline is H1 for almost all
*integers*, with an O(2^n/n²) exceptional set. The transfer to S_n misses by a constant.
- **Theorem H1-S (§4; CONDITIONAL on Proposition V).**
  - Setting: fix c ∈ (0, 1/2) and K, and put y = T = ⌊2^{cn}⌋.
  - Claim: all but o(2^n/n²) integers N ∈ [2^{n−1}, 2^n) have at least π(y) + K y-smooth values among
    (⌈√N⌉ + t)² − N, 1 ≤ t ≤ T.
  - By the transfer lemma this holds for all but o(|S_n|) balanced semiprimes.
- **Inputs:**
  - Bombieri–Vinogradov for smooth numbers at fixed u: Fouvry–Tenenbaum 1996 Théorème 1, or Harper's
    preprint arXiv 1208.5992 Thm 1 (ineffective any-A version);
  - the uniform Dickman asymptotic Ψ(x, y) = xρ(u)(1 + O(1/log y)) on compact u-ranges.
- **Four ideas carry the proof:**
  1. *The region is not thin.* After summing over t at fixed gap k = t′ − t, the pair sum runs over a
     *fat* 2-D region with lattice index 2k (Lemma C). It is a Barban–Davenport–Halberstam/BV-type sum
     at level ≈ T, far below √Ψ. polyimpl §7's "thin coupled region" holds only for a fixed pair
     (t, t′); it does not describe the aggregate. This settles the D1 kill criterion: **not
     triggered.**
  2. *The unstratified variance has exact order n^{−2} (Prop L).* So unstratified Chebyshev gives
     exceptional sets of size Θ(2^n/n²), which miss the transfer by a constant. This is a genuine
     local (Knuth–Schroeppel-type) effect, not a proof artefact.
  3. *Stratify* by N mod W (W a product of prime powers) and by the position of N in its dyadic block
     (resolution 1/J).
  4. *Average the local covariance over k (proof strategy for V; GAP).* Within a stratum, the
     main-term covariance at the primes of W should cancel. The primes outside W and the archimedean slot leave r(W, J)/log² y with
     r → 0. An iterated limit then gives o(n^{−2}).
- **Codex reached ideas 1–3 independently.** It stated the stratified estimate (its (*)) as
  UNVERIFIED. Idea 4 is our proof of a sufficient form of (*).

**H2 (fruitfulness): OPEN. We have only small results and a precise target.**
- **Lemma P (PROVED, elementary).** Any dependency of size 2 is fruitful: if Q_tQ_{t′} = Y² with
  t ≠ t′ and T = o(√N), then N < X − Y < X + Y < 2N. Sizes 1 and 2 are the only sizes where a size
  argument can work (§6.2).
- **Prop S (PROVED).** For c < 1/2, a pair dependency exists for at most a B^{c−1/2+o(1)} fraction of
  N. So the size mechanism is dead for almost all N.
- **Barrier.** lemmaM Theorem B blocks every fixed-test averaging argument, because the kernel depends
  on N. Theorem A balances Jacobi signs on a *fixed* set, but the smooth set depends on N. The
  coupling is exactly what nobody can control.
- **Missing statement H2\*(c):** #{N ∈ S_n : H1 holds and ψ_N is trivial on ker A_N} = o(|S_n|).

**Conditional consequence (Theorem DQ, §7).**
- H1-S plus H2\*(c) give a uniform deterministic algorithm. It costs N^{ωc+o(1)} on every input and
  factors all but o(|S_n|) of S_n.
- With c < 1/(6ω) ≈ 0.0703 this is below the project's 1/6. With H2\*(c) for every c it is N^ε for
  every ε > 0.
- It is an *almost-all* statement, not an average-time one: the failure set is only o(|S_n|), far
  above the N^{−1/30} needed to beat AV in mean.

**Subexponential (y = L[1/2]).**
- Fixed-u H1-S does not extend. Any residue stratification admissible for the BV step leaves a
  residual variance ≫ (log n)²/n² against the n^{−2} needed (Prop B, PROVED IN OUTLINE).
- So the residue-stratified second-moment method fails there by a factor (log n)². Either higher
  moments (multi-point smooth correlations on a quadratic family) or a multiplicative local model
  would be needed.

**Experiments (pre-registered; frozen verdicts "mixed" for both H1 structure and H2).**
- The local-factor theory is confirmed strongly:
  - the model L(N) explains 97–99% of Var(Z_N) at 40–72 bits (corr 0.985–0.996);
  - the pair-correlation ratios by k mod ℓ match the model (e.g. n = 72: 1.115/1.108, 1.059/1.058,
    1.069/1.062 for ℓ = 2, 3, 5);
  - no large-k anomaly.
- The frozen failures are all small-size local extremes. One example is a 40-bit semiprime that is a
  non-residue modulo every prime 3–19 and has Z = 0. H1 fails for 97, 24, 1 and 2 of 1000 semiprimes
  at 40, 48, 56 and 64 bits, and for 0 at 72.
- H2:
  - the Jacobi rule held in 100% of cases;
  - 203/203 size-≤2 dependencies were fruitful (Lemma P);
  - post hoc, the per-vector trivial rate is 0.470–0.514 at 56–80 bits, but 0.566 at 40 bits (a
    small-size excess, also seen by below16 E4).

**Literature (§9).**
- We found no deterministic almost-all (or average-case) factoring algorithm below N^{1/5} other than
  the project's own D′/AV. This is a search result, not a novelty certificate.
- We also found no rigorous almost-all H1 theorem for the QS polynomial.

---

## 1. Setting

- B = 2^n, I = [B/2, B) ∩ Z, P = 2^{n/2}.
- S_n = {pq : p < q ≤ 2p primes} ∩ I, with |S_n| ≍ B/n².
- f = 1_{P⁺(·) ≤ y}.
- m = ⌈√N⌉, Q_t = (m+t)² − N, and Z_N = Σ_{t≤T} f(Q_t).
- Parameters: y = T = ⌊B^c⌋, common to the whole dyadic block (codex: if T varied with N, its
  variation would enter the variance). u = log(√B·T)/log y → 1 + 1/(2c).
- **Transfer lemma (PROVED, trivial).** If |{bad N ∈ I}| = ε_n|I|, then the bad fraction of S_n is
  O(ε_n n²).

## 2. Exact coordinates: the fat region

**Lemma C (PROVED; codex has the same identities).**
- Put x = m + t and v = x² − N.
- The map (N, t) ↦ (x, v) is a bijection from I × [1, h] onto
  D_h = {(x, v) : x² − v ∈ I, 2x − 1 ≤ v ≤ 2x(h+1) − (h+1)² − 1}.
- For fixed t, v runs over [2xt − t², 2x(t+1) − (t+1)² − 1]. These intervals adjoin as t grows.
- For a gap k ≥ 1, put w = Q_{t+k} = v + 2kx + k². Then

  Σ_{N∈I} Z_N² = Σ_N Z_N + 2 Σ_{k<T} Σ_{(x,v)∈D_{T−k}} f(v) f(v + 2kx + k²).

- In the (v, w)-plane the k-th sum runs over R_k ∩ Λ_k, where:
  - Λ_k = {w − v ≡ k² (mod 2k)} has index 2k;
  - x = (w − v − k²)/(2k);
  - R_k is a region of dimensions ≍ √B T (along v) by ≍ √B k (along w − v), cut out by the inequalities
    of D_{T−k} and x² − v ∈ I.

**Consequences.**
- *Fixed (t, t′) is thin; the sum over t is not.* For fixed (t, t′) the image of (m, r) is a sliver of
  relative width ≍ k/t² (polyimpl §7). Summing over t at fixed k tiles these slivers into R_k.
- For k ≍ T, R_k is fat in both directions. For K₀ < k ≪ T its aspect ratio is k/T ≥ (log B)^{−3}.
  These intervals have relative length (log B)^{−O(1)}, and arbitrary log savings absorb them.
- *Level.* The congruence moduli are q = 2kW ≤ 2TW = B^{c+o(1)}. The relevant smooth counts have size
  Ψ ≍ √B T ρ. So q ≤ Ψ^{1/2}B^{−(1/4 − c/2)+o(1)}, a polynomial margin below Harper's level √Ψ for
  every c < 1/2.
- *Correction to polyimpl §7 / D1.* The kill criterion ("the thin-region correlation cannot be
  bounded…") is **not triggered**. The pair sum is a BV/BDH-type sum, not a binary additive problem.
  below16 §6.2 called the second moment "a binary additive problem (semiprime + smooth = square)". That
  was a statement about averaging over *semiprimes*. Over all integers, plus the transfer, it becomes
  linear-congruence sums at low level.

**First moment (PROVED; codex §2).** E Z_N = T(ρ(u) + O_c(1/log B)).
- Exchange x and v in Lemma C.
- The x-multiplicity is a piecewise-smooth weight, up to O(1) per v; edge effects are O(T/√B).
- Then apply partial summation of Ψ(x, y) = xρ(u)(1 + O(1/log y)).

## 3. The unstratified variance has exact order n^{−2}

**Prop L (PROVED IN OUTLINE; codex §1 and §4, re-derived here).**
- *Local mean.* Fix an odd prime ℓ. Then

  (E[Z | N ≡ a (ℓ)] − E Z)/T = −ρ′(u) log ℓ (a/ℓ) / ((ℓ − 1) log y) + o_ℓ(1/log y).

  - The number of x mod ℓ with ℓ | x² − N is 1 + (N/ℓ).
  - Smooth v ≡ 0 (ℓ) number Ψ(z/ℓ, y) ≈ ℓ^{−1}Ψ(z, y)(1 − ρ′(u) log ℓ/(ρ(u) log y)).
- *Variance.* By the law of total variance,

  liminf (log y)² Var(Z)/T² ≥ ρ′(u)²(log ℓ)²/(ℓ(ℓ − 1)) > 0.

- **So "Var Z/(E Z)² = o(n^{−2})" is FALSE.**
- *Upper bound.* Codex also proves Var Z = O_c(T²/log² y). Its route is the §2 identity, Harper's
  Theorems 1–2 and first-order expansions; the linear local terms cancel. This gives:
  - **Theorem H1-int (PROVED IN OUTLINE; codex):** #{N ∈ I : Z_N < π(y) + K} = O_c(B/(log B)²).
  - That is H1 for almost all *integers*. But the exceptional set is too large for the transfer to S_n
    by exactly a constant.
- *Source of the n^{−2}.*
  - (i) The Legendre symbols (N/ℓ) at small ℓ change how often ℓ divides Q_t. This is the effect
    Knuth–Schroeppel multipliers exploit. For p-adic valuations ≥ 2 the relevant data is N mod ℓ^a,
    not N mod ℓ (codex).
  - (ii) The *archimedean* place: Q_t ≈ 2t√N, so log Q_t moves by up to (1/2) log 2 across the
    dyadic block. That moves ρ(u) by a relative ≍ ξ(u)/log y, with ξ = −ρ′/ρ.
  - Both have relative size ≍ 1/log y per place, so squares ≍ n^{−2}. Neither is a fluctuation of
    smoothness itself.

## 4. Theorem H1-S: stratification and averaging over k

**Strata.**
- Let W = ∏_{p ≤ z} p^{a} with z, a ≥ 1 fixed, and J ≥ 1 fixed.
- A stratum σ = (s, j) is {N ∈ I : N ≡ s (mod W), N ∈ I_j}, where I_j = [B/2·(1 + j/J),
  B/2·(1 + (j+1)/J)).
- Put M_σ = E[Z | σ].

**Proposition V (OPEN: a GAP per the codex referee; the strategy below is not certified).** For fixed c ∈ (0, 1/2), W and J,

  limsup_{n→∞} (log B)² · E_N[(Z_N − M_{σ(N)})²]/T² ≤ C_c · r(W, J),

  r(W, J) := Σ_{p > z} (1 + log p)²/p² + Σ_{p ≤ z} Σ_{e ≥ a} (1 + e log p)² p^{−e/2} + J^{−2}.

  The p^{−e/2} (not p^{−e}) is needed because when p^a | N, x² ≡ N (mod p^e) only needs p^{⌈e/2⌉} | x.

**Proposition M (PROVED IN OUTLINE).** For fixed W and J, max_σ |M_σ/T − ρ(u)| = O_{W,J}(1/log B).
- Proof: first moments in progressions of the fixed modulus W, on fixed sub-intervals. Harper's
  Theorem 1 at the single modulus W, or Fouvry–Tenenbaum, gives them. Every local factor at a fixed
  place is 1 + O(1/log y).

**On the range of c.** The BV step needs only c < 1/2. Codex states its theorem for c < 1/6, which is
the range relevant to beating 1/6, not a constraint of the proof. The < 1/6 algorithm (§7) needs
c < 1/(6ω).

**Theorem H1-S (CONDITIONAL on V; the deduction from V and M is PROVED).** Fix c ∈ (0, 1/2), K ≥ 0 and ε > 0.
1. Choose W, J with 32 C_c r(W, J) < ερ(u)². This is possible because r → 0 as z, a, J → ∞.
2. For n ≥ n₀, Prop M gives M_σ ≥ ρT/2 for every σ. Also π(y) + K ≤ ρT/4, since T = y.
3. A bad N therefore has |Z − M_σ| ≥ ρT/4.
4. Chebyshev and V give #{bad}/|I| ≤ 16 E(Z − M_σ)²/(ρT)² ≤ ε/(log B)² for n ≥ n₁.
5. Since ε is arbitrary, #{bad} = o(B/n²). The transfer lemma then gives almost all of S_n. ∎
- The constants in Harper's arbitrary-A savings are ineffective. So n₀(ε) is ineffective. The
  algorithm itself is uniform and explicit.

**Proof of Proposition V (outline, eight steps).** Fix σ = (s, j). Write E_σ and Cov_σ for averages over
N ∈ σ, and put f_t = f(Q_t(N)).
1. **Expansion.** Var_σ(Z) = Σ_t Var_σ(f_t) + 2 Σ_{k=1}^{T−1} Σ_{t≤T−k} Cov_σ(f_t, f_{t+k}). The first
   sum is ≤ T = o(T²/log² B).
2. **Near-diagonal.** For k ≤ K₀ := T/(log B)³, use |Cov| ≤ 1. The total is ≤ K₀T = o(T²/log² B).
3. **Coordinates and boxes.**
   - For k > K₀, the coupled sum Σ_t E_σ[f_t f_{t+k}] is |I_σ|^{−1} Σ f(v) f(w) over (v, w) ∈ R_k ∩ Λ_k
     with x² − v ∈ σ (Lemma C).
   - The stratum adds the congruence x² − v ≡ s (mod W). Here x mod W is fixed by (w − v) mod 2kW.
   - So the conditions are (v, w) ≡ (a, b) (mod q), q = 2kW, for (a, b) in a set 𝒜_k(s).
   - Each a has at most ν(W) partners b, and each b at most ν(W) partners a. Here ν(W) bounds the
     number of square roots mod W.
   - Cover R_k by product boxes J × J′ of side √B k/(log B)^{A₁}. The boundary boxes carry a fraction
     O((log B)^{−A₁}) of the lattice points; with f ≤ 1 and A₁ = 3 they cost o(T²/log² B) in total.
   - The archimedean slot j is an interval condition on x² − v and is absorbed the same way.
4. **BV replacement.**
   - On a box, Σ_{(a,b)∈𝒜} A_J(a) A_{J′}(b), where A_J(a) = #{v ∈ J smooth, v ≡ a (q)}.
   - Replace A by Ā_J(a) := Ψ_{q/g}(J/g, y)/φ(q/g), g = gcd(a, q), which is the mean over the classes
     with the same g.
   - For g ≤ G := (log B)^R, the errors E = A − Ā are reduced-class errors at modulus q/g and scale
     J/g. The bounded multiplicities ν(W) give:
     - |ΣĀE′| ≤ ν(W) max Ā Σ_b|E′_b|;
     - |ΣEE′| ≤ ν(W)‖E‖‖E′‖.
   - Summed over k (distinct k give distinct q ≤ 2TW), Harper's Theorem 1 at the scales J/g bounds
     everything by ≪_A Ψ(J)Ψ(J′)(log B)^{−A+O(R)}. Theorem 2 is an alternative for the ‖E‖ terms.
   - Classes with g > G: the v ∈ J with gcd(v, q) > G number ≪ Ψ(J) Σ_{g|q, g>G} g^{−α}, by a crude
     Hildebrand–Tenenbaum bound. Averaged over k this is ≪_W G^{−α+o(1)}. It is negligible for
     R > 3/α.
   - The product term Σ_t E_σ f_t E_σ f_{t+k} gets the same replacement. Each E_σ f_t is a
     long-interval average: for fixed t, v ranges over ≍ t√B values as N runs over σ. So Harper at the
     single modulus W suffices, with error (log B)^{−A} uniformly in t ≥ K₀.
5. **Main term = a local covariance.**
   - Write δ_J(a) := q Ā_J(a)/Ψ(J), a function of g = gcd(a, q) only.
   - On each box, (main coupled) − (main product) = (Ψ(J)Ψ(J′)|𝒜|/q²)·Cov_μ(δ_J(a), δ_{J′}(b)). Here
     μ is the uniform law on 𝒜_k(s), and its two marginals equal the product-term marginals.
   - The archimedean dependence of δ is constant up to 1 + O(1/(J log y)) on a slot. That is the J^{−2}
     term.
6. **Splitting δ.**
   - Let q = q_W q′, with q_W the W-primes part and q′ coprime to W. Then q′ | 2k and g = g_W g′.
   - Since w − v ≡ 0 (mod k), g′ = gcd(v, q′) = gcd(w, q′): v and w carry the same g′.
   - Write δ = D_W(g_W)·κ(g_W, g′), where:
     - D_W(g_W) = q_W Ψ_{q_W/g_W}(J/g_W)/(φ(q_W/g_W)Ψ(J)) is the W-only density;
     - κ is the q′-correction.
   - From Ψ(X/d, y) = (X/d)ρ(u − log d/log y)(1 + O(1/log y)), with inclusion–exclusion over the
     divisors of q_W (boundedly many prime factors) and of q′, we get:
     - (LD-a) |D_W − 1| ≤ C_W/log y;
     - (LD-b) |κ(g_W, g′)/κ(g_W, 1) − 1| ≤ C Π(q′)(1 + log g′)/log y for g′ ≤ G, with
       Π(q′) := ∏_{p|q′}(1 + C/p);
     - (LD-c) |κ(g_W, 1)/κ(1, 1) − 1| ≤ C Π(q′) σ(q′)(1 + log g_W)/log y, with σ(q′) := Σ_{p|q′} 1/p.
       The q′-dependence enters only through divisors d > 1 of q′, of weight 1/d.
7. **Averaging over k.**
   - Let g′ have its law under μ: it is independent of the W-data by CRT, with Pr(g′) ≈ 1/g′. Put
     κ̄ = E_{g′}κ and Cov_μ(δ, δ′) = T₁ + T₂ + T₃:
     - T₁ = Cov over the W-data of D_W κ̄ and D′_W κ̄′;
     - T₂ = E[D_W D′_W (κ − κ̄)(κ′ − κ̄′)], with the same g′ in both factors;
     - T₃ = interaction terms.
   - **T₂ (large primes).** By (LD-b), |T₂| ≲ Σ_{g′|q′, g′>1} Π²(1 + log g′)²/(g′ log² y) plus the
     square of the corresponding first moment. Averaged over k, Pr(g | k) = 1/g. Codex's convergent
     sums Σ_g(1 + log g)²/g² and Σ_{d,e}(d,e)(1 + log de)²/(d²e²), restricted to g, d, e > 1 coprime
     to W, give ≲ r(W)/log² y.
   - **T₁ (W-primes).**
     - κ̄ = λ_k(1 + μ_k(g_W)), with λ_k = 1 + O(σ(q′)) and |μ_k| ≲ Π σ(q′)(1 + log g_W)/log y by
       (LD-c).
     - The μ_k parts are products of two first-order quantities: ‖μ_k‖₂ · ‖D′_W − E D′_W‖₂ ≲
       σ(q′)/log² y. Averaged over k, σ(q′) has mean ≲ Σ_{p>z} 1/p², which is part of r(W).
     - The main part is λ_kλ′_k · Cov(D_W(v_W), D′_W(w_W)), with (v_W, w_W) = (x² − N, (x + k)² − N) at
       the W-primes and x uniform. Truncate D_W at valuation a: the error is ≲ Σ_{p|W} Σ_{e≥a}
       (1 + e log p)²p^{−e/2}/log² y, part of r. This uses |δ(e) − 1| ≲ e log p/log y and the
       convergence of Σ_e p^{−e/2}p^{2e(1−α)}, which holds since α → 1. The truncated covariance then depends on k only through
       k mod W.
     - λ_k depends only on the primes ∤ W of k, so it is asymptotically independent of k mod W.
     - **Averaging k over complete residue systems mod W makes (x, x + k) a uniform independent pair.**
       The averaged covariance is then Cov_N(E_x D_W(x² − N), E_{x′}D′_W(x′² − N)) over N in the
       stratum. This is the variance of local data finer than N mod W, again ≲ Σ_{p|W} Σ_{e≥a} p^{−e/2}(1 + e log p)²/log² y. After truncation it is in fact exactly 0: the
       truncated D_W depends on v only mod W, and v ≡ x² − s (mod W) is fixed by x and s.
     - The archimedean slot contributes ≲ J^{−2}/log² y in the same way.
   - **T₃** collects the products of (κ − κ̄) with D_W-fluctuations. Each is a product of two
     first-order quantities, one of which carries σ(q′) or g′ > 1. It is ≲ r(W)/log² y on average.
8. **Summation.** Sum over the boxes of R_k and over k ∈ (K₀, T). The weights Ψ(J)Ψ(J′)|𝒜|/q² add up to
   ≍ (ρT)²·|I_σ|, so Σ_k Σ_t Cov_σ ≤ (C r(W, J) + o(1))T²/log² y. ∎

**Known gaps in the proof of V (codex detqs_2 #4–#11). They are substantive, not routine.**
- (#6) Step 6 cannot use Ψ = xρ(u)(1 + O(1/log y)). On boxes of relative length (log B)^{−A}, endpoint
  subtraction loses the main term. A differentiable de Bruijn–Saias (Λ) approximation is required.
- (#7) LD-a is false uniformly. The W-part of q carries the valuations of k, so |D_W − 1| ≍
  (1 + log g_W)/log y. Weighted moment bounds are needed instead.
- (#5) Step 5's "marginals coincide" is false as written: the product term averages over independent
  N, N′ with t-dependent weights. Both weighted main terms must be derived and the archimedean
  covariance kept.
- (#8, #9) LD-b/LD-c need the cancellation shown before taking absolute values. The "asymptotic
  independence" of λ_k from k mod W needs a divisor expansion with summable tails and CRT.
- (#10) The valuation truncation must use conditional projections modulo p^a, uniformly over singular
  strata. The parity condition w − v ≡ k² (mod 2k) and 2 | W must be carried throughout.
- (#2, #4, #11) Discard k close to T and small t; use endpoint grids for Harper; account for repeated
  moduli q/g; give a weighted pair bound for g > G; compute the archimedean J^{−2} term explicitly.

**Remarks.**
- Codex's (*) asks more than V: second-order coefficients and an explicit centering. V uses exact
  conditional means (stratification), so no second-order coefficient is ever evaluated. Only
  first-order bounds (LD-a)–(LD-c) and the exact cancellation of step 7 are needed.
- The Hildebrand short-interval theorem and Matomäki–Radziwiłł are not needed after aggregation (codex
  agrees).
- The algorithm is unchanged: deterministic, no multiplier, a single polynomial. Only the analysis is
  stratified.

## 5. The subexponential regime: a barrier for the method

**Prop B (PROVED IN OUTLINE; a barrier for the method, not for H1).**
- Take y = L[1/2, a], T = L[1/2, b], and any stratification by N mod W, with W ≤ B^{1/4−ε} so that
  step 4's level condition holds.
- *Variance from unstratified primes (repaired after codex #13).*
  - Use only primes p ≍ D log B, with D a large fixed constant. Such p are < y, so the local mean formula
    of Prop L applies uniformly.
  - A modulus with log W = O(log B) contains only O(log B/log log B) of these primes, so most of the
    band stays unstratified.
  - Remove from each Legendre character (N/p) its conditional projection onto N mod W. Bessel's
    inequality for the resulting nearly orthogonal family then gives

    E Var(Z | N mod W)/(E Z)² ≫ ξ(u)²/log² y · Σ_{p≍D log B, p∤W} (log p)²/p² ≍ ξ(u)²/log² y · log log B/log B.

- *Size.* With ξ ≍ log u ≍ (1/2) log n and log² y ≍ a² n log n, this is ≍ (log n)²/n². That exceeds
  n^{−2} by the factor (log n)². The constant is UNVERIFIED.
- So Chebyshev at this variance cannot give o(n^{−2}). This is a barrier for the variance-based
  certification, not for H1.
- *Other obstacles at L[1/2]* (codex): close-pair deletion loses 1/ρ(u), and the first-order local
  perturbation squared is ≍ (log n)/n.
- *Routes not ruled out:*
  - a 4th or 6th moment, which needs multi-point smooth correlations along v_i = v + k_i(2x + k_i), a
    quadratic family no BV theorem covers;
  - a multiplicative local model Z ≈ ρT∏_{p≤y}λ_p(N)(1 + η), whose lower tail is a large-sieve
    estimate for character sums. This needs main terms at moduli ∏_{p≤y} p, which is beyond every
    distribution theorem.
- Both are **OPEN**.

## 6. H2: fruitfulness

### 6.1 Proved facts

**Lemma P (PROVED; codex's proof).**
- *Statement.* Let √N < x < x′ ≤ √N + T with T = o(√N), and (x² − N)(x′² − N) = Y². Put X = xx′. Then
  N < X − Y < X + Y < 2N. Hence gcd(X − Y, N) is a proper divisor of N. In particular, no such pair
  exists when N is prime.
- *Proof.*
  - (X − N)² − Y² = N(x′ − x)² > 0, so X − Y > N.
  - X ≤ (⌈√N⌉ + T)² and Y ≤ (Q + Q′)/2 ≤ (⌈√N⌉ + T)² − N, so X + Y < 2N once
    2(⌈√N⌉ + T)² < 3N. This holds eventually since T = o(√N).
  - N divides (X − Y)(X + Y) = N(N + Q + Q′) but divides neither factor. So the gcd is neither 1
    (that would force N | X + Y) nor N. ∎
- *Size 1.* Q_t = s² gives N = (m + t − s)(m + t + s) (Fermat). The factor m + t − s exceeds 1 because
  m + t + s < N.
- Tests: `test_lemma_P_*`. E2 observed 203/203 size-≤2 dependencies fruitful.
- Lemma P is a special case of Lehman's nontriviality argument with multiplier (x′ − x)² (literature
  critic). We found no explicit QS-phrased statement.

**Prop S (PROVED).** For fixed c < 1/2,

  #{(N, t, t′) : N ∈ I, t < t′ ≤ T, Q_tQ_{t′} = □} ≤ B^{1/2+c+o(1)}.

- So a pair dependency exists for at most a B^{c−1/2+o(1)} = o(1) fraction of N.
- *Proof.*
  - In the coordinates of Lemma C, a triple is a pair (v, w) with v, w ≤ 3√B T and vw = □, plus a
    factorisation w − v = k(2x + k).
  - With V := 3√B T, the pairs with vw = □ number Σ_{s sqfree ≤ V} ⌊√(V/s)⌋² ≤ Σ_{s≤V} V/s ≪
    √B T log B.
  - Each such (v, w) has ≤ τ(w − v) = B^{o(1)} factorisations. ∎
- *Data.* The observed counts fall with n: 138, 43, 16, 5, 1, 0 small dependencies at 40–80 bits.

**Prop SZ: WITHDRAWN (FALSE; codex detqs_2 #15).**
- The draft claimed that no conjugation pattern gives a size certificate for |S| ≥ 3. A size certificate
  does not need A < N. It needs only that A − Y and A + Y lie strictly between the same pair of
  consecutive multiples of N.
- *Counterexample (codex, re-verified).* N = 1572001 = 1097·1433, x = (1255, 1257, 1258), so
  Q = (3024, 8048, 10563) with product Y², Y = 507024.
  - The conjugation (+, +, −) gives A = 13249776 ≡ X (mod N).
  - Then 8N < A − Y < A + Y < 9N, so gcd(A − Y, N) = 1097.
- So for |S| ≥ 3 the size mechanism is not ruled out in principle.
- What remains true is weaker. For |S| ≥ 3, |A| is typically ≥ N^{1+c}, so the window between
  consecutive multiples of N has relative width ≍ N^{−c}. Nothing forces A ± Y into the same window.
  We state no theorem about it.

**Jacobi form (KNOWN; below16 §6.1).**
- For Blum N, H2 ⇔ ψ_N(S) = (X_S Y_S/N) is a nontrivial character on ker A_N.
- By reciprocity, for odd x_i: ψ_N(S) = ∏_{i∈S}(−Q_i/x_i)·(2/N)^{e_2(S)/2}.
- Expanding (ℓ/x_i) = (x_i/ℓ)(±1), the quartic symbol of the chosen root s_ℓ cancels because e_ℓ(S)
  is even. What remains is the parity of how many i ∈ S use the root −s_ℓ (mod ℓ), over the
  factor-base primes ℓ ≡ 3 (mod 4), plus explicit mod-8 signs.
- So H2 fails exactly when this computable "root-label" vector lies in the column space of A_N
  (corrected by the 2-adic term). This is a reformulation, not progress. It shows that H2 is a
  statement about the ideal-level refinement of the exponent matrix in Z[√N].

### 6.2 Barriers and the missing statement

- **Fixed-test averaging is blocked.** lemmaM Theorem B (codex, proved): fixed tests can cancel for
  almost all n while an N-dependent predictor annihilates every n. Our kernel ker A_N depends on N
  through the smoothness pattern. Theorem A (Jacobi balance on a *fixed* smooth set) therefore does not
  apply. Writing Pr[H2 fails] = E_N[|ker_N|^{−1}Σ_{S∈ker_N}ψ_N(S)] makes the obstacle explicit: a fixed
  S gives a fixed-test average, but the weight |ker_N|^{−1} and the support ker_N are N-dependent.
- **Rank alone does not help.** Codex: a kernel of growing dimension need not contain a vector of
  weight 2. Prop S shows that weight-2 vectors are absent for almost all N.
- **B1 is not an obstruction here.** lemmaM B1 concerns even-|S| two-form sieves without a known root.
  For QS it fails: N = 33, 7² − 33 = 16. The least-nonresidue phenomenon does appear, but only in H1's
  local factors: the Z = 0 semiprime of §8 is a QNR modulo 3–19.
- **Multipliers and MPQS** (deterministic choice of k or of polynomials) multiply the relation families.
  They supply no independence for ψ: the same coupling problem recurs for each family.
- **H2\*(c) (OPEN).** For fixed c ∈ (0, 1/2):

  #{N ∈ S_n : Z_N ≥ |FB| + 2 and every kernel vector S has 1 < gcd(X_S − Y_S, N) fail} = o(|S_n|).

  - For Blum N this is: ψ_N is trivial on ker A_N for o(|S_n|) of N.
  - It suffices to test a basis.
  - The fair-coin heuristic predicts a failure rate 2^{−dim ker} ≤ 2^{−(T ρ − π(y))}, which is
    super-polynomially small.

## 7. Conditional algorithm, cost, and comparison with AV

**Theorem DQ (CONDITIONAL on H2\*(c); the cost is PROVED).** Algorithm A_c, given N:
1. Set n = bit length and y = T = ⌊2^{cn}⌋.
2. Check perfect squares and trial-divide by the primes ≤ y.
3. Sieve Q_t exactly for t ≤ T. This costs T·(log log y)·n^{O(1)}.
4. Compute a basis of ker A_N by deterministic elimination: O(y^ω) with LSP/LUP (Ibarra–Moran–Hui), or
   O(y³) by Gauss.
5. Test every basis vector by gcd.

- *Cost:* N^{ωc+o(1)} on every input.
- *Correctness:* A_c factors every N ∈ S_n outside (H1-S exceptions) ∪ (H2\*(c) exceptions), i.e. all
  but o(|S_n|).
- With c < 1/(6ω) ≈ 0.0703 (cubic elimination: c < 1/18), the exponent is < 1/6. If H2\*(c) holds
  for every fixed c, then for every ε > 0 there is a deterministic algorithm factoring almost all of
  S_n in time N^ε.
- **Not** average time: mean cost with an HH fallback needs failure density ≤ N^{−1/30−δ}. H1-S gives
  only o(n^{−2}), and AV's N^{1/6} mean is not improved.
- **Comparison.** AV is unconditional, mean N^{1/6}polylog, worst case N^{1/5}. DQ is conditional on
  H2\*, almost-all, at any fixed exponent. Theorem DR (IW hypothesis) gives conditional deterministic
  L[1/2] for *every* N; DQ's hypothesis is unrelated, an equidistribution statement about one sieve.
- Unconditionally, nothing below 1/6 follows. We do not improve on AV.

## 8. Experiments (pre-registered; `notes/detqs_prereg.md`, commit 8a6b789)

**Disclosed deviations and errors.**
- The pooled E2 binomial p-value overflowed in float (M = 4600). After the data existed, the function
  was rewritten in log space; the test is unchanged.
- **Model error (found by codex detqs_2 #20).** `localmodel.u_blocks` used m ≈ 1.2247·2^{n/2}. The
  typical value is √(3/4)·2^{n/2} ≈ 0.866·2^{n/2}, so the frozen model overstated m by √2. That shifts
  u_t by about (1/2) log 2/log y.
  - The frozen L(N), V_L and MODEL_PAIR_RATIO values contain this error. We did not re-fit them; any
    corrected fit would be post hoc.
  - The checks that depend on it (P1b, P1c, P1d, P3b) passed anyway.
- "Explained variance" below means 1 − V_res/V for the multiplicative residual r = Z/(cL). It is not a
  regression R². corr² (0.97–0.99) gives a similar figure.

**E1 (variance and lower tail of Z_N).** y = 2^⌊n/5⌋, T = 2^{16..20}.

| n | Z̄ | V = Var(Z/Z̄) | Poisson 1/Z̄ | V/V_L | corr(Z, L) | V_res/V | V_strat/V | min r (int / sp) | H1 fails (sp) |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 116 | 0.713 | 0.0086 | 0.88 | 0.985 | 0.031 | 0.65 | 0.14 / 0.00 | 97/1000 |
| 48 | 283 | 0.548 | 0.0035 | 0.95 | 0.990 | 0.019 | 0.69 | 0.34 / 0.47 | 24/1000 |
| 56 | 867 | 0.267 | 0.0012 | 1.00 | 0.994 | 0.012 | 0.63 | 0.71 / 0.68 | 1/1000 |
| 64 | 1463 | 0.203 | 0.0007 | 1.02 | 0.995 | 0.010 | 0.61 | 0.80 / 0.62 | 2/1000 |
| 72 | 5248 | 0.127 | 0.0002 | 1.02 | 0.996 | 0.008 | 0.54 | 0.87 / 0.90 | 0/500 |

- Frozen checks:
  - P1a–P1d hold at every n;
  - P1h holds at every n (mean r over semiprimes 0.985–0.998);
  - **P1e fails at 40–64** (V_strat/V = 0.61–0.69 > 0.6);
  - **P1f fails at 40 and 48**;
  - **P1g fails at 40–64**.
- **Frozen verdict: "mixed".**
- Post hoc (labelled):
  - The local model explains 97–99% of the variance. The residual is 2.6–5.5× Poisson.
  - Every lower-tail failure is a local extreme. For the Z = 0 semiprime at 40 bits, brute force confirms
    Z = 0: N is a QNR modulo every prime 3–19, and the factor base is {2, 23, 29, …}.
  - Such N have density ≍ 2^{−π(z)} for fixed z. Asymptotically their effect on M_σ is 1 + O(1/log y)
    (Prop M). At 40–64 bits, log y is 5.5–8.3 and the effect is a factor 3–10.
  - Restricting the strata to (N mod 8, 3, 5) explains less than the full model (P1e). Most of the
    local variance sits at primes 7 ≤ ℓ ≤ y. This is consistent with the T₂/r(W) structure of §4.

**E3 (pair correlation ratio by k mod ℓ, k ∈ [T/8, T/2); measured / model).**

| n | ℓ = 2 | ℓ = 3 | ℓ = 5 |
|---|---|---|---|
| 40 | 1.477 / 1.512 | 1.277 / 1.261 | 1.258 / 1.293 |
| 48 | 1.350 / 1.396 | 1.215 / 1.206 | 1.231 / 1.229 |
| 56 | 1.213 / 1.202 | 1.119 / 1.108 | 1.123 / 1.118 |
| 64 | 1.177 / 1.168 | 1.100 / 1.090 | 1.100 / 1.097 |
| 72 | 1.115 / 1.108 | 1.059 / 1.058 | 1.069 / 1.062 |

- P3a and P3b hold at every (n, ℓ).
- The pair structure is local and model-accurate, with no thin-region or large-k anomaly. The
  excess decays as predicted by the local model's bias²/ℓ, which tends to 0 like 1/log² y.

**E2 (H2 rank event on S_n).**
- P2a held: the Blum Jacobi rule was right for 1098 runs with at least one tested vector, and for every
  vector.
- P2f held: 203/203 size-≤2 dependencies were fruitful.
- **P2b failed** (23/30 per-(n, j) tests; pooled failed), **P2c failed** (107 > 2) and **P2d failed** (H1
  fails in 104, 30 and 2 runs at 40, 48 and 56 bits).
- **Frozen verdict: "mixed".**
- Post hoc (labelled):
  - Of the 107 runs with no fruitful vector, 105 had too few relations (H1 failures), so they are
    censored runs.
  - The other **two are genuine finite H2 failures** (both at 40 bits, non-Blum):
    - N = 560915181371 (Z = 25, |FB| = 23);
    - N = 738702714041 (Z = 26, |FB| = 23).
    Each exhausted its relations with a 4-dimensional kernel, every basis vector trivial. The expected
    count of such runs under the fair-coin model is ≈ 1000·2^{−4}·Pr(dim = 4), so two is not
    surprising.
  - Per tested kernel vector, the trivial rate is:

    | n | 40 | 48 | 56 | 64 | 72 | 80 |
    |---|---|---|---|---|---|---|
    | rate | 0.566 | 0.522 | 0.497 | 0.514 | 0.510 | 0.470 |
    | z | +6.1 | +1.9 | −0.3 | +1.2 | +0.6 | −1.2 |

  - The 40-bit excess of trivial vectors matches below16 E4's 40-bit excess. It is unexplained and
    vanishes by 56 bits.
  - The data are consistent with H2\* and with the fair-coin heuristic at 56–80 bits. They say nothing
    about the asymptotic statement.

**Reproduction.**
- `PYTHONPATH=. .venv/bin/python -m src.detqs.harness e1`: about 5 min.
- `… e2`: about 1.5 min.
- `… summary`.
- `.venv/bin/python -m pytest -q tests/test_detqs.py`.

## 9. Literature (own search, codex §4, and the Opus literature critic, §10)

- **Rigorous relation methods are randomized.**
  - Dixon 1981 (random squares);
  - Croot–Granville–Pemantle–Tetali, Annals 175 (2012), give the sharp threshold for Pomerance's
    random-squares problem;
  - Vallée, Math. Comp. 56 (1991) 823–849, rigorous randomized L[1/2, √(4/3)];
  - Seysen 1987 and A. K. Lenstra 1988 (class-group relations, rigorous under GRH);
  - Lenstra–Pomerance, JAMS 5 (1992) 483–516, rigorous randomized L[1/2, 1 + o(1)], no GRH.
- **Schnorr–Lenstra 1984** (Math. Comp. 43) is a randomized group-order method. Its conjectured L[1/2]
  bound was withdrawn (Lenstra–Pomerance 1992), so it is not a rigorous relation method.
- **Lee–Venkatesan**, J. Number Theory 187 (2018), is rigorous for every n because it randomizes the
  NFS polynomial. It is the closest precedent for "make a sieve rigorous by averaging over a family".
  Our averaging is over N, with a deterministic algorithm.
- **QS analyses are heuristic.**
  - Pomerance 1982 (Math. Centre Tracts 154) and 1985 (LNCS 209), and "Smooth numbers and the quadratic
    sieve" (MSRI 44, 2008), which works out the *conjectured* complexity.
  - Granville's MSRI 44 survey (§4.3) states that a rigorous basic QS needs ≈ y^{1+o(1)} smooth values of
    one given quadratic, and calls this open.
  - The smooth-values-of-polynomials literature (Hmyrova, Balog–Ruzsa, Hildebrand 1989, Dartyge 1996,
    Dartyge–Martin–Tenenbaum 2001, Martin 2002, Bober–Fretwell–Martin–Wooley arXiv 1710.01970) treats a
    fixed polynomial with the variable far beyond the height. It never averages over the constant term N
    with t ≤ N^c.
- **Deterministic rigorous factoring.**
  - Exponent 1/5: Harvey, Math. Comp. 90 (2021) 2937–2950, building on Hittmeir's 2/9 (Math. Comp.
    2018).
  - Harvey–Hittmeir, Math. Comp. 91 (2022), remove a log-log factor.
  - Lehman 1974 gives N^{1/3}; Hart 2012 is heuristic for general N.
- **Closest "deterministic almost-all" formulation.** Dąbrowski–Pomykała–Shparlinski, "On oracle
  factoring of integers" (arXiv 1912.00345). It is deterministic polynomial time and finds one nontrivial
  factor of almost all n, but only *with an elliptic-curve point-count oracle*. It differs from our
  setting in three ways: oracle, one factor, and all integers rather than S_n.
- **Smooth-number tools used.**
  - Harper, "Bombieri–Vinogradov and Barban–Davenport–Halberstam type theorems for smooth numbers",
    arXiv 1208.5992, an **unpublished preprint** (v1 2012), Theorem 1. Statement checked by the critic:
    log^K x ≤ y ≤ x, Q ≤ √Ψ, max over reduced classes, main term Ψ_q/φ(q). At fixed u only the
    **ineffective** any-A version saves anything, so n₀ in H1-S is necessarily ineffective.
  - Primary published alternative at fixed u: **Fouvry–Tenenbaum, Proc. LMS (3) 72 (1996) 481–514,
    Théorème 1**. Its √x Q term is B^{1/4+3c/2+o(1)} ≪ Ψ/log^A exactly when c < 1/2, the same
    condition. Granville 1993 (Phil. Trans. A 345) is another alternative.
  - Theorem 2 is not needed: Theorem 1 with a max bound handles the ‖E‖ terms.
  - Dickman asymptotic Ψ = xρ(u)(1 + O(log(u+1)/log y)): Hildebrand, J. Number Theory 22 (1986)
    289–307; Hildebrand–Tenenbaum survey (JTNB 5, 1993) Theorem 1.1. The survey's Theorem 1.8 is
    Saias's Λ-approximation.
  - Ψ(X/g, y) ≪ g^{−α}Ψ(X, y): Hildebrand–Tenenbaum 1986 (survey Corollary 1.7).
- **Lemma P** is a special case of Lehman's nontriviality argument with a square multiplier d². Indeed
  (X − N)² − Y² = N(x′ − x)² is a Lehman representation A² − B² = kN with k = (x′ − x)². Size 1 is
  Fermat; the CFRAC analogue is SQUFOF's square form. We found no explicit QS-phrased statement.
- **Search result.** We found no unconditional, oracle-free, deterministic almost-all (or average-case)
  factoring algorithm below N^{1/5} other than the project's D′/AV. We found no almost-all H1 theorem for
  (⌈√N⌉ + t)² − N. "No prior formulation found by us" at most.
- **Caveat.** Finding *one* factor of almost all integers is easy (trial division). *Complete*
  factorization of almost all integers is not trivial, and S_n is the case where the n^{−2} transfer
  bites.

## 10. Reviews

**GPT-6.1-sol adversarial referee** (`codex/detqs_2.out`). Valid findings applied:
- V/H1-S downgraded from "proved in outline" to GAP/CONDITIONAL; gaps listed in §4.
- Prop SZ withdrawn (counterexample re-verified).
- Prop B proof repaired (primes ≍ log B, Bessel); its rate (log n)²/n² is confirmed.
- Lemma P's quantitative condition stated.
- §8:
  - model error disclosed;
  - E3 scope narrowed;
  - the two genuine 40-bit H2 failures reported;
  - rate range corrected.
- Its labels: C, M, transfer, P, S PROVED; L/H1-int correct outline needing a stronger stated input;
  V/H1-S GAP; B correct barrier; SZ FALSE; DQ conditional on H1-S and H2\*; H2\* OPEN; experiments
  "mixed".
- Not yet applied: in DQ, deterministic root preparation should be charged y^{2+o(1)}, and the theorem
  must be distinguished from the implementation, which uses incremental Gauss with a 30-vector cap.

**Opus literature critic** (applied in §9):
- Harper is an unpublished preprint; only its ineffective version saves at fixed u.
- Fouvry–Tenenbaum 1996 is the published alternative.
- The Dickman citation is HT survey Thm 1.1 / Hildebrand 1986.
- Schnorr–Lenstra is not rigorous; the exponent 1/5 is due to Harvey 2021.
- Added Dąbrowski–Pomykała–Shparlinski (oracle) and Lee–Venkatesan.
- Lemma P is a special case of Lehman's argument.
- "Trial division" applies to one factor only.

**Opus proof critic:** launched but had not reported when this note was frozen. Its findings are not
included. No theorem here is claimed as passing two critics.

## 11. Candidate card (BRIEF §4)

**Card DQ-1: deterministic single-polynomial QS on almost all S_n.**
- *Mechanism.* Exact sieve of (⌈√N⌉ + t)² − N for t ≤ T = y = N^{c}, deterministic elimination, gcd
  tests.
- *Input family.* S_n, almost all.
- *Proved.*
  - H1 for all but o(|S_n|) (Theorem H1-S, in outline, conditional on Harper Thm 1 and the uniform
    Dickman asymptotic);
  - the size mechanism for H2 (Lemma P), available only for o(1) of N (Prop S, Prop SZ);
  - cost N^{ωc+o(1)}.
- *Open.* H2\*(c). The barrier is lemmaM Theorem B: N-dependent kernels.
- *Experiments.* The local-factor structure is confirmed. The H2 fair-coin rate holds at 56–80 bits, with a
  40-bit excess.
- *Kill criterion of polyimpl D1.* Not triggered for H1 (the region is fat). For H2, D1's own
  criterion does not apply. The obstruction is the N-dependent kernel, for which no averaging
  technique is available.
- *Most promising open question.* H2\*(c) for a single fixed c < 1/(6ω), even for a positive
  proportion of S_n. A positive-proportion version would already give a deterministic
  positive-proportion N^{1/6−δ} result, which the project does not have below 1/6.
