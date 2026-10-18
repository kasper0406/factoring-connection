# cycloharvey: can cyclotomic Frobenius handles beat Harvey's N^{1/5}?

Agent `cycloharvey`, 2026-10-18. Pre-registration: `notes/cycloharvey_prereg.md`
(commit `1d92a48`, made before the search). Code: `src/cycloharvey/`. Tests:
`tests/test_cycloharvey.py`.

## 0. Verdict

**No.** This is refuted for the model and family we froze, and in that model it is
proved.

1. **The cyclotomic handles add no new search targets** (proved, §3.1). Every
   exponent F(p,q) that σ_r makes reachable is already reachable in plain (Z/N)^*,
   where a^{F(p,q)} ≡ a^{F(1,N)} (mod p) for every integer polynomial F. Any division
   trick that comes from the guessed residue r ≡ p (mod m) is also available in
   (Z/N)^* by guessing the same residue. The premise that "higher-degree Taylor
   cancellation is not possible in F_p^*" is therefore false. It was always available;
   the calibrated algorithms just do not use it. The ring adds only cost:
   D = φ(m)-dimensional arithmetic, and φ(m) residue guesses where (Z/N)^* is free to
   skip them.
2. **Higher-degree Taylor cancellation does not help** in the
   Lehman/Hittmeir/Harvey model with any bounded degree.
   - **Proved lower bound** (§3.4). Every covering family of bounded degree, with no
     guess or with residue guesses modulo a single modulus M, costs ≥ N^{1/5−o(1)}.
     An independent Opus reviewer found it correct after the scope fixes, which are
     applied here.
   - **Scaling heuristic** (§3.3). It reproduces the known exponents and gives exponent
     S/(K(K+3)) for a family with K free coefficients on exponents of total size S.
     Over all 63 degree-≤3 supports this is minimised uniquely by Lehman's {q, p},
     at 1/5. The next best is 3/14 for {q², q, p, p²}.
   - **Why the user's hoped-for δ = 5/2 is not reachable.** Cancelling the
     second-order term needs coefficients on p² or q². Those are quantised at scale
     N rather than √N, which inflates the window. The effective δ falls below 3/2;
     for example K=4 gives 4/3.
   - **Guessing residues is exactly exponent-neutral** (heuristic, with a baby list
     shared across classes). The integer-valued lattice gains index M^{K(K+1)/2},
     which cancels the φ(M) ≈ M guesses. The heuristic reproduces Harvey–Hittmeir's
     (log log N)^{3/5} gain. The proved bound only limits the gain to M^{o(1)}.
3. **Experiments agree.**
   - The exact constructive search over all degree-≤3 exponents (LLL plus exact
     verification) finds nothing better than Lehman's family. See §4.
   - An operational Harvey baseline matches the model's list counts exactly: 112/112
     instances, 32–56 bits.
   - The cyclotomic version factors too, but at ≈ φ(m)·D times the work.

## 1. Calibration gate (passed)

Model (§2) cost = list elements. All window widths are exact integers
(`src/cycloharvey/calibration.py`). Exponents come from two sources:
- exact rational LP (`exponent_lp`): **proved** inside the model;
- least-squares fits of log2(cost) on log2 N for exact model counts at
  N = 2^n + 2^{n/2} + 1, with parameters optimised on a 2^{1/4} grid.

| algorithm | model mechanism | LP exponent | fitted slope (upper half) | sizes n |
|---|---|---|---|---|
| BSGS on p+q (Pollard 1978; Hittmeir 2018 without sieving) | 1 target, W ≈ (c^{−1/2}+c^{1/2}−2)√N | **1/4** | 0.2500 | 40–160 |
| Lehman 1974 | targets k ≤ r, windows √N/(4r√k), square test (m = 1) | **1/3** | 0.3334 | 30–66 |
| Hittmeir 2021 | per-(a,b) BSGS for ab ≤ ξ; one giant step each for ξ < ab ≤ η | **2/9** (ξ = N^{1/9}, η = N^{2/9}) | 0.2342 | 40–96 |
| Harvey 2021 | all ab ≤ r, ⌈W_ab/m⌉ giant steps, shared baby list m | **1/5** (r = m = N^{1/5}) | 0.2099 | 40–120 |
| Harvey–Hittmeir 2022 | φ(M)·#σ₀-buckets single giant steps, baby list 2λ+1, λ = ⌈4√N/(M m₀)^{3/2}⌉ | **1/5** (T = M m₀ = N^{1/5}) | 0.1997 | 40–128 |

- All fits lie within 0.012 of the LP exponent, which passes the pre-registered gate of
  ±0.03. The residual excess for Hittmeir and Harvey is the log-factor drift: Harvey
  carries lg r from Σ_{ab≤r} 1, and Hittmeir is tiered.
- The HH slope sits slightly below Harvey's. The φ(M)/M ~ 1/log log N gain is
  visible as a downward drift.
- Sources were read from the arXiv PDFs: Harvey 2010.05450 §§3–4, Hittmeir
  2006.16729 §5, Harvey–Hittmeir 2105.11105 §§3–4. Lehman is read via Harvey's
  Lemma 3.3 and Remark 3.4. We did not read the 1974 paper itself, so the
  constant-level details of the Lehman row are **UNVERIFIED** against the original.
- Raw data: `experiments/raw/cycloharvey_calibration.jsonl`. Summary:
  `experiments/analysis/cycloharvey_calibration.json`.

### 1a. The user's T = MK heuristic against the real proofs

The user's derivation is in `notes/user_ideas/2026-10-18_mod_M2_lattice_T_MK.md`:
B ~ (MK)^{1/2}, R ~ M^{1/2}K^{−3/2}N^{1/2}, S = R/M², and T = MK = N^{1/5}.

**Where it agrees.** It is exactly HH's Prop. 3.3 bookkeeping, with K = m₀ and
M = m (checked line by line):
- The lattice of (a,b) with b ≡ (N/σ²) a (mod m) has index m. Gauss reduction gives
  |aN/σ₀| ≤ 2N^{1/2}(m m₀)^{1/2}, which is the coefficient inflation B ~ (MK)^{1/2}.
- The second-order term t₀²·δ·aN/σ₀ with t₀ < 1/m₀ is ≤ 2N^{1/2}m^{1/2}m₀^{−3/2}.
  This is R, with constant 2.
- The baby list is λ = 4N^{1/2}/(m m₀)^{3/2} = R/M² (step 9, units of m²).
- The cost is φ(m)·m₀·lg N + λ, i.e. T + S. The balance gives T = N^{1/5}.

**Where the heuristic and the proofs differ** (all checked):
1. **First order is balanced, not cancelled.** HH do not kill the first-order term.
   The 2-D Minkowski box makes |−aN/σ₀ + bσ₀|·t₀ the same size as the second-order
   term, also 2N^{1/2}m^{1/2}m₀^{−3/2}, so the bound in (3.3) is 4N^{1/2}m^{1/2}m₀^{−3/2}.
   The exponents match the heuristic; the mechanism is a two-term balance.
2. **Geometric buckets.** The "K Taylor centres of width √N/K" are geometric buckets
   σ₀ → ⌈(1+1/m₀)σ₀⌉ in HH, so there are O(m₀ lg N) of them, not K. For balanced p
   only O(m₀) are needed.
3. **T is not MK.** T = φ(M)·m₀ (one branch per unit class), not M·m₀. The ratio
   φ(M)/M ≍ 1/log log N is the entire HH improvement.
   - The volume heuristic with κ = K(K+1)/2 and a shared baby list gives exactly a
     (φ(M)/M)^{3/5} gain, which is HH's (log log N)^{3/5} (§3.3). This is
     heuristic.
   - The proved lower bound (§3.4) is weaker. It allows a gain of up to M^{o(1)},
     which includes O(1)^{ω(M)}τ(M) factors that can be much larger than
     (log log N)^{3/5} for primorial M. Optimality of HH's log-log exponent is
     therefore **not** proved here (the reviewer caught this overclaim in an earlier
     draft).
4. **Harvey 2021 uses no residue classes.** Its branches are all pairs ab ≤ r with
   windows √N/(4r√(ab)) (Lemma 3.3, Lehman/Dirichlet). Per-branch steps vary with ab;
   the exponent is the same.
5. **The heuristic omits two costs that are binding exactly at 1/5.**
   - Deterministic large-order element (Prop. 2.7, √D with D ≥ N^{2/5}): N^{1/5}.
     Harvey says "good enough … but only just".
   - Removal of small factors p ≤ √(N/r) by Strassen in Harvey: (N/r)^{1/4} = N^{1/5}.
   Any sub-1/5 search would also need these improved. Hittmeir's Remark 6.4 suggests
   D = N^{1/3+o(1)} suffices, giving N^{1/6}. A randomised base removes the first cost.
6. **Polylog factors.** lg^{16/5} N from the product tree/Bluestein and the per-branch
   lattice reduction. These do not affect exponents.

**Generalised yardstick.** Suppose per-branch steps are S ~ N^{a}/T^{δ}. Then T =
N^{a/(1+δ)}. Here a = 1/2 is not fixed: a support with K coefficients and exponent
mass S_I has a = S_I/(2K) and δ = (K+1)/2 (§3.3). Folding the inflation of a into
δ gives δ_eff = K(K+3)/(2S_I) − 1, and the exponent is 1/(2(1+δ_eff)).

## 2. The abstract cost model (step 1 components)

- **(a) Window.** A target is a function g(p) = F(p, N/p) with a window
  {U, …, U+W−1}. With t = p/√N, g = Σ_k b_k t^k with b_k ∈ (√N)^{|k|}Z (**proved**:
  p^i q^j = N^{min(i,j)} p^{i−j}). Coverage is required for every integer p in the
  balanced range [√(N/c), √N]. This is the "identity model": the correctness proof
  uses only the polynomial identity, as Lehman, Hittmeir, Harvey and HH all do.
- **(b) Strides / residue locks.** Guessing p ≡ r (mod M) (φ(M) runs) admits the
  integer-valued lattice L_{r,M} ∋ (F−c)/M^e. HH Prop. 3.3 is the K=2 case, the
  index-M sublattice fixed mod M². The index of L_{r,M} over Λ_I is M^κ, computed
  exactly by HNF (`lattice.class_lattice_index`).
- **(c) Batching.** cost = m + Σ_τ ⌈W_τ/m⌉, with one shared baby list. This is Harvey's
  "single BSGS over the whole space". It is optimal among sumset strategies, since
  any baby set B and giant set G_τ covering W_τ candidates need |B|·|G_τ| ≥ W_τ.
  Hittmeir's two-tier scheme and Lehman's square test are restrictions of it.
- **(d) Cyclotomic overhead.** φ(m) guesses, plus every ring operation on
  D = φ(m) coefficients (≥ D base operations; schoolbook D², FFT D log D). The handle
  T_{r,F} costs O(deg F · log N) ring multiplications. Constructing elements whose
  order divides Φ_d(p) costs a further product of σ-conjugates.
- **(e) One-sidedness.** No cost in the model (**proved**). A two-sided collision
  (≡ mod N) is caught by sort-and-match (Harvey Step 3, HH Step 12) and gives the
  integer u = g(p), from which p follows by solving g(p) = u. A one-sided collision is
  caught by Alg. 4.1. Experiment E3 (§5) confirms codex's observation, and also finds
  that p² ≡ q² (mod m) suffices for F = p²+q² to be two-sided. Neither changes the
  cost.

## 3. Theory

### 3.1 D1: cyclotomic handles are dominated by (Z/N)^* handles (proved)

**Claim.** Let r ≡ p (mod m), α ∈ R^*, and T = ∏ σ_r^{i−j}(α)^{c_ij N^j}. Suppose
the identity α^{u(p)} ≡ T (mod p) holds for every α and every p ≡ r. Then u(p) ≡
F(p,q) (mod p^f−1), where f = ord_m(p) and R/p ≅ ∏ F_{p^f}. Any identity-model target
built this way is therefore a function u ∈ F + (p^f−1)·Z[p,q], or (F−c)/M^e if α
is an M^e-th power. Every such u lies in Z[p, N/p], or in the integer-valued lattice
L_{r,M}.

**Matching target in (Z/N)^*.** For the same u, take a = b^{M^e}. Then
a^{u} = b^{F(p,q)−c} ≡ b^{F(1,N)−c} (mod p), because p ≡ 1 and q ≡ N (mod p−1). This
handle is computable from N and the guess, and the baby list is a^i. Windows, branch
set and correctness proof are identical. Therefore:

  e_cyc(family) ≥ e_{F_p^*}(same family) + log_N D,

and (Z/N)^* may also drop the guesses when they are not used.

Two caveats from the independent review:
- **Base order.** "Same windows" needs ord_p(b^{M^e}) to exceed the baby-list length.
  In (Z/N)^* that order can be as small as (p−1)/gcd(M^e, p−1), whereas F_{p^f} has
  more room. The lower bound below does not use the matched handle; it only uses the
  inclusion of the target sets (cyclotomic targets ⊆ Λ ∪ L_{r,M}). So the cost
  comparison is stated for parameter ranges where the (Z/N)^* order is not binding.
  This holds for M^e ≪ p/W, which covers all calibrated algorithms.
- **Multi-base searches.** Handles that combine several σ_s with s ∉ ⟨r⟩ as
  *independent* search directions are multi-base vector-exponent searches. They lie
  outside the 1-D cost model, as do multi-base searches in (Z/N)^*. See §3.5 for
  their separable part.

**The (Z/N)^* handle collapse is a feature, not a loss.** Z[p,q]/(pq−N, p−1) ≅ Z:
one handle a^{c} serves every F with F(1,N) = c simultaneously.

**Genuinely independent cyclotomic handles.** In the subgroup of order dividing
Φ_d(p) (d | f), the exponent ring is Z[p,q]/(pq−N, Φ_d(p)) ≅ Z[x]/Φ_d(x), because
q = N x^{−1} and x^{−1} is a polynomial mod Φ_d. It has rank φ(d). For example, for
d = 4: p² ≡ −1, p³ ≡ −p, q ≡ −Np, q² ≡ −N². So only σ^0, …, σ^{φ(d)−1}, weighted by
powers of N, are independent handles. This sharper structure means fewer targets
share a handle, not more.

Cost of making elements of order dividing Φ_d(p): β^{(p^d−1)/Φ_d(p)} =
β^{∏_{e|d,e<d}Φ_e(p)}, a σ-product of degree < d (for d=4, codex's σ²(β)/β). It is
useful only when the residue degree f is a multiple of d, and that is not guaranteed by
the guess. All of this is pure overhead in the model, because of D1.

E2 checks the identity T_{r,F} ≡ α^{F(p,q)} (mod p) for 1440 cases, including
negative coefficients via the norm inverse α^{−1} = ∏_{s≠1}σ_s(α)/Nm(α): 0 failures.

### 3.2 Lemmas for the lower bound

**Lemma W (Wronskian, proved; exhaustively checked for all 63 supports).**
det[(i)_k]_{k=1..K, i∈I} = ∏_{i∈I} i · ∏_{i<i'}(i'−i) ≠ 0 for 0 ∉ I.

*Proof.* (i)_k is a monic degree-k polynomial in i with zero constant term. Row
reduction gives det[i^k]_{k=1..K}, which is ∏i times a Vandermonde.

Consequence: for t in a compact subinterval of (0, ∞),
(G'(t), …, G^{(K)}(t)) = M(t)·b with M(t) invertible, so Σ_{k≤K}|G^{(k)}(t)| ≥ c_I‖b‖_∞.

**Lemma C (coverage, proved).** Let G = Σ_{k∈I} b_k t^k, |I| = K, and t ∈ J =
[c^{−1/2}, 1]. For any measurable E ⊂ ℝ,

  |{t ∈ J : G(t) ∈ E}| ≤ C_I·(|E|/‖b‖_∞)^{1/K}.

*Proof.*
- Split J into O_I(1) intervals on which G is monotone with G' monotone (the zeros of
  the Laurent polynomials G' and G'' are bounded in number).
- On each piece the preimage density 1/|G'| is monotone in the value. So for fixed
  |E| the preimage measure is maximised by an interval E* at the flat end.
- Further split into O_I(1) pieces on which some |G^{(k)}| ≥ c_I‖b‖/K (Lemma W).
  Apply the van der Corput sublevel lemma |{|G − y| ≤ ε}| ≤ C_k(ε/λ_k)^{1/k}.
- For |E| ≤ ‖b‖, (|E|/‖b‖)^{1/k} ≤ (|E|/‖b‖)^{1/K}. If |E| > ‖b‖ the bound is
  trivial with C_I ≥ |J|.

Each window {U..U+W−1} is a real set of measure W. Integers p ≡ r (mod M) in a
p-interval of length √N·μ number ≤ √N μ/M + 1.

**Lemma K′ (denominators of the guessed lattice, proved).** If v = Σ α_k p^k (with
α_k for k<0 the coefficient of q^{|k|}) lies in L_{r,M}, then C_d·M^K·α ∈ Z^K, where C_d
depends only on d = max|k|.

*Proof, at a prime ℓ^a ‖ M.* The function f(x) = v(r + ℓ^a x) = Σ_n T_n ℓ^{an} x^n is
ℓ-integral-valued on Z_ℓ, where T_n = Σ_k α_k N^{max(0,−k)} binom(k,n) r^{k−n}
(N^{|k|} is an ℓ-unit). So v(T_n) ≥ v* := min v(α_k). Since (T_1..T_K) = W(r)α with v(det W) = c_W (= 0 if ℓ > 2d), some j ≤ K
has v(T_j) ≤ v* + c_W.

Take j₀ minimising v(T_j) for j ≤ K and use Δ^{j₀}f(0) = j₀! Σ_{n≥j₀} S(n,j₀)T_nℓ^{an}
∈ Z_ℓ. For a > c_W the n = j₀ term dominates, so v(T_{j₀}) ≥ −a j₀ − v(j₀!). Hence
v* ≥ −aK − v(K!) − c_W. If a ≤ c_W, pass to the subclass mod ℓ^{c_W+1}; the extra is
bounded. ∎

The HNF table (`experiments/analysis/cycloharvey_kappa.json`) shows the denominator
exponent is exactly K(K+1)/2 per prime ℓ > 2d. It is larger by a bounded amount for
ℓ ∈ {2, 3}, and for ℓ = 5 when d = 3.

### 3.3 Scaling heuristic, δ yardstick, and guess-neutrality (heuristic)

The heuristic matches BSGS (1/4), Lehman (1/3), Harvey/HH (1/5) and HH's
(log log)^{3/5}. Hittmeir's 2/9 comes from his two-tier restriction of the cost model
(§1), not from this formula.

**Volume heuristic.** A target covers a t-interval of length ℓ around t₀ iff its
Taylor coefficients satisfy |c_j(t₀)| ≲ Wℓ^{−j} for j = 1..K. Such a lattice vector
exists for typical t₀ iff the box volume W^K ℓ^{−K(K+1)/2} ≳ covol = (√N)^{S}, with
S = Σ_{k∈I}|k|. This gives

  ℓ = (W^K N^{−S/2})^{2/(K(K+1))},  #targets = 1/ℓ,  cost ≍ 1/ℓ + √(W/ℓ),

and so

  **e(I) = S/(K(K+3))**, with a = S/(2K), δ = (K+1)/2, δ_eff = K(K+3)/(2S) − 1.

Per-branch steps are W ≍ N^{S/(2K)}·T^{−(K+1)/2}. Reproduced cases:
- {1} gives 1/4 (BSGS on p, or on p+q).
- {−1,1} gives 1/5 (Harvey/HH). For Lehman's typical a ≈ √r pairs the heuristic ℓ
  equals Harvey's window bookkeeping.
- The same formula with the Lehman m = 1 restriction gives 1/3.

All 63 supports of degree ≤ 3 (`search.py analytic`, exact rationals):

| support I (k: p^k, k<0: q^{|k|}) | K | S | e = S/(K(K+3)) | raw δ=(K+1)/2 | a=S/(2K) | δ_eff |
|---|---|---|---|---|---|---|
| {−1, 1} (Lehman: aq+bp) | 2 | 2 | **1/5** | 3/2 | 1/2 | **3/2** |
| {−2,−1,1,2} | 4 | 6 | 3/14 ≈ 0.214 | 5/2 | 3/4 | 4/3 |
| {−2,−1,1}, {−1,1,2} | 3 | 4 | 2/9 ≈ 0.222 | 2 | 2/3 | 5/4 |
| {−3,−2,−1,1,2,3} (all of degree ≤ 3) | 6 | 12 | 2/9 ≈ 0.222 | 7/2 | 1 | 5/4 |
| {±1, ±2, 3} etc. (K=5) | 5 | 9 | 9/40 = 0.225 | 3 | 9/10 | 11/9 |
| {1} or {−1} | 1 | 1 | 1/4 | 1 | 1/2 | 1 |
| {−3,−1,1,3} (codex's s³−12Ks) | 4 | 8 | 2/7 ≈ 0.286 | 5/2 | 1 | 3/4 |

**Effective δ for the cyclotomic families** (the orchestrator's yardstick, after
all costs; m = N^μ, D = φ(m)):

  δ_eff,cyc(I) = 1/(2(e(I) + μ_D + μ_guess)) − 1,

where:
- μ_D = log_N D comes from the D-dimensional arithmetic (D1);
- μ_guess = 0 with a baby list shared across the φ(m) guesses, or
  μ·(2/(K+3)) with per-guess lists;
- coefficient inflation is already in e(I) through a = S/(2K);
- one-sidedness costs 0 (§2(e)).

So δ_eff,cyc ≤ δ_eff(I) ≤ 3/2 for every family. Equality requires I ⊇ {−1,1} used as
Lehman vectors, with bounded m.

**Every support other than {−1,1} is strictly worse** (test
`test_min_heuristic_exponent_is_lehman_support`).

**The user's δ = 5/2.** This is the raw δ of K = 4: three Taylor orders cancelled
and the fourth as remainder. It would need a = 1/2, i.e. S = K. But four distinct
nonzero exponents force S ≥ 1+1+2+2 = 6, so a = 3/4. The p² and q² coefficients are
quantised at scale N, not √N, and that blow-up exactly undoes the gain: δ_eff = 4/3 <
3/2. This is the general form of codex's s³−12Ks failure.

**Guess-neutrality.** With M = N^μ and lattice index M^κ, the heuristic gives

  e(I, μ) = S/(K(K+3)) + μ·(1 − 2(K+κ)/(K(K+3))).

With κ = K(K+1)/2 (HNF-measured, §3.2), the μ-coefficient is exactly 0 for every K.
Guesses are exponent-neutral.

- This requires the baby list to be **shared across the φ(M) classes**, as in HH,
  where the baby steps β^{m²i} do not depend on σ.
- With a separate baby list per run, the μ-coefficient is 2/(K+3) > 0, and guessing
  hurts.
- κ = K(K+1)/2 is measured by HNF, not proved. Lemma K′ proves only that the
  per-coordinate denominators are ≤ C_d M^K. The lower bound in §3.4 needs only
  Lemma K′. Keeping φ(M) separate from M gives a gain factor
(φ(M)/M)^{power}, which for K=2 is (φ(M)/M)^{3/5}. The cyclotomic version adds +μ
(D ≈ m).

### 3.4 Theorem LB: lower bound in the identity model (proved within the model's scope; independently reviewed)

**Statement.** Fix d ≥ 1 and c > 1. Consider any covering family for the balanced
range [√(N/c), √N]:
- targets are Laurent polynomials of support ⊆ [−d,d]∖{0}, from Λ or from L_{r,M}
  for guessed classes modulo **one** modulus M (any subset of classes);
- windows are arbitrary integer intervals;
- total cost is X, either per-run (Σ_runs (m + Σ_τ⌈W_τ/m⌉)) or with one baby list
  shared across runs. The proof uses only "#distinct target vectors ≤ X and
  Σ W_τ ≤ X²", which holds in both.

Then X ≥ c_{d,c}·N^{1/5}·M^{−o(1)}.

**Proof.** Write P = √N. Group targets with the same (class, vector) and use Lemma C
for the union of their windows. The number of distinct vectors is ≤ X and
Σ W_τ ≤ X·max m ≤ X². Covering needs Σ_τ (P μ_τ/M + C_I) ≥ c′P, since the preimage of a window is a union
of ≤ C_I intervals; if X < c′P/2 this means
Σ_τ μ_τ ≥ c″ M. Here μ_τ ≤ C(W_τ/‖b_τ‖)^{1/K_τ}. Let e = max|k| over the support of
b_τ, so K_τ ≤ 2e.

- **e ≥ 2.** By Lemma K′, ‖b‖ ≥ P^{e}/(C_d M^{K}), so μ_τ ≤ C M P^{−e/K} W_τ^{1/K}.
  Hölder over ≤ X vectors with Σ W ≤ X² gives Σ μ ≤ C M P^{−e/K} X^{(K+1)/K}. M
  cancels, and e/(K+1) ≥ e/(2e+1) ≥ 2/5. So this class needs X ≥ c·P^{2/5}, or it
  covers only a small fraction.
- **e = 1** (K ≤ 2, vectors (Aq+Bp)/(C M²), ‖b‖ ≍ P·max(|A|,|B|)/M²).
  - For fixed (A,B), the number of classes r (mod M) for which it is integer-valued is
    ≤ O(1)^{ω(M)}·gcd(A,B,M), because each class must satisfy B r² ≡ A N modulo a
    divisor of M. The reviewer brute-forced this for M ≤ 12; it held. So #{(class, vector): ‖b‖ ≤ V} ≤ M^{o(1)}(VM²/P)².
  - Hence the X smallest have Σ 1/‖b‖ ≤ M^{o(1)} X^{1/2} M²/P.
  - Cauchy–Schwarz: Σ μ ≤ C(ΣW)^{1/2}(Σ 1/‖b‖)^{1/2} ≤ C M^{1+o(1)} P^{−1/2} X^{5/4}.
  - Requiring this ≥ c″M gives X ≥ c P^{2/5} M^{−o(1)}.
- **Combining.** Summing the finitely many (e, K) classes: if X ≤ δ P^{2/5}M^{−o(1)},
  each class covers ≤ C δ^{const}·M, which is too little for small δ. Hence
  X ≥ c N^{1/5} M^{−o(1)}. ∎

**Scope (rule 5).** This is a lower bound *in the model*:
- identity-proved coverage of all p in the range;
- bounded-degree polynomial exponents;
- a single guess modulus;
- sumset-type BSGS with one handle per target and one base.

With several moduli M_j the e = 1 count gains a factor #moduli, and the bound
degrades to about N^{1/6}. That case is **not covered**.

It does not cover proofs that use arithmetic of the specific N, non-polynomial
exponent maps, or other collision mechanisms. It is not a factoring lower bound.

The e = 1 case with M = 1 is the statement that Harvey's 1/5 is optimal for Lehman's
family under full sharing, consistent with Harvey's remark that N^{1/6} would need a
"fully square-root speedup" for r ≍ N^{1/3}.

**Review.** An independent read-only Opus 5.5 reviewer checked D1, W, C, K′, LB and
the heuristic.
- Verdicts: CORRECT for W, C and K′; D1 CORRECT within the model, with the two scope
  caveats now added; LB CORRECT for a single modulus.
- Required fixes, all applied above: restrict LB to one modulus; reconcile the
  shared and per-run wording; withdraw the (φ/M)^{3/5}-optimality claim; add the
  constant C_I per window; add the N^{|k|} unit in K′.
- Status: **proved, single reviewer**.

### 3.5 Separability (user's cross-difference law) — enforced by the model (proved)

Every test in the model has the separable form E(τ, j; i) = [g_τ(p) − jm] − i.
- The shared baby list is α^i with the single exponent h(i) = i.
- Each giant step is T_τ·α^{−jm}.
- This is exactly the user's condition E(j,u) = f(j) − h(u).

**Can Frobenius handles enter separably?**
- **Giant (handle) side: yes, in full.** σ^k(α) appears only inside T_τ, which is the
  whole point of the handle. It changes nothing beyond the target set, and by D1
  that set is already available in (Z/N)^*.
- **Baby side: only by relabelling.**
  - Using σ^k(α)^i = α^{p^k i} as baby steps makes the collision i·p^k ≡ g_τ(p)
    (mod ord). That is i ≡ p^{f−k}g_τ(p) (mod p^f − 1), a separable search for the
    target p^{f−k}g_τ, still in Λ (or L_{r,M}) but of higher degree, hence worse
    (§3.3).
  - Mixing several σ^k(α) powers in one baby list (i₀ + i₁p + …) makes h(i) depend on
    the unknown p. This is non-separable (cross-difference ≠ 0) and is not
    multipoint-testable with one product tree. Splitting it into separable pieces
    gives a multi-base search outside the 1-D model, priced additively.
- **Decoupling.** The user's decoupling argument (a_i ≡ a₀ (mod M), with M-inflation
  cancelling the M² gain) is the K = 2 instance of guess-neutrality: index M^{3} from
  M² plus first order, against M runs.

### 3.6 Elkies-prime analogue (qualitative; not part of the frozen search)

At an Elkies prime ℓ for E mod p, guessing the eigenvalue λ makes the substitution
y(X) ↦ y(x([λ]X)) computable modulo the kernel polynomial h_ℓ. It agrees with the
p-power Frobenius on A = (Z/N)[X]/h_ℓ (mod p), and is a ring endomorphism of
dimension D = (ℓ−1)/2.
- **For exponent targets F(p,q) it changes nothing.** The argument of D1 applies
  verbatim. A/p ≅ ∏F_{p^f}, so the identity-model targets are again Z[p,N/p] (or
  L_{r,ℓ}), all realisable in (Z/N)^*. The costs are about ℓ guesses (λ, and whether
  ℓ is Elkies for p) and D-dimensional arithmetic, i.e. dominated, as for σ_r.
- **What it adds.** The guess carries information about the trace: λ + p/λ ≡ t_p
  (mod ℓ). That matters only for targets in the *curve* group, with exponents
  involving #E(F_p) = p+1−t_p, a new unknown of size N^{1/4} that is unrelated to q's
  curve.
- **Status.** Such a "Z[p, q, t_p]" target lattice is outside our frozen family and
  is **not yet tested**. A naive guess is that the new monomial t_p carries
  quantisation P^{1/2} and the same volume obstruction applies, but we have not
  computed it.

### 3.7 Relation to the Coppersmith exponent bridge (user note; codex r8)

The user's bridge takes H in L_{m,d} with k = H(x₀)/p^m = Σ c_ij x₀^i q^{m−j} and
a^{H(1−p₀)} ≡ a^k (mod p). It is exactly a target of our lattice.
- k is a Laurent polynomial in p with q^{j}-coefficients divisible by N^{j}, i.e. an
  element of Λ_I with I ⊆ [−m, d−m]. The bucket centre p₀ is our Taylor centre, and
  the handle is our F(1,N)-type handle.
- Therefore:
  - m = 1 is the e = 1 (Lehman) class, and m ≥ 2 brings q² and q³ terms (e ≥ 2).
  - Our LB (§3.4) covers every fixed (m, d) family in the identity model. Its e ≥ 2
    branch gives cost ≥ c·P^{e/(K+1)} ≥ c·N^{1/5}, with equality of the *crude*
    exponent exactly at {±1,±2}. The heuristic gives 3/14 there.
- Codex r8 (`notes/codex_r8_copp_bridge.md`) finds that the L_{m,d} determinant
  bounds cannot certify < 1/5: their minimum over m ≥ 2 at β = 3/10 is exactly 1/5,
  attained at (m,d) = (2,4) and (2,5). The (2,4) case is the same equality case as
  ours.
- Codex stresses that determinant bounds are upper bounds and do not by themselves
  give a no-go. LB supplies the missing lower bound for bounded-degree families in
  this model, subject to its scope (single modulus, identity coverage, 1-D separable
  BSGS).
- The user's claim 2 (the phase-transition slope) is consistent with this but is not
  re-derived by us (**UNVERIFIED** here). Claim 3 ("exponent vs base ring") is refuted
  by codex r8; we agree, since Harvey's Alg. 4.1 is itself a base-ring polynomial
  evaluation.
- **Prior art for the bridge** (codex r8): Gao–Feng–Hu–Pan, arXiv 2512.19076,
  §§1.1, 4, 6. They convert a Coppersmith polynomial's hidden quotient into an
  exponent collision mod p.

## 4. Pre-registered search results

**Analytic objective (A).** We evaluated S/(K(K+3)) exactly over all 63 supports
(`search.py analytic`).
- The minimum is 1/5, attained only by {−1, 1}. The next best is 3/14.
- With guesses at the HNF-measured κ = K(K+1)/2 (primes M > 2d), the μ-coefficient is 0
  in (Z/N)^*. This assumes a baby list shared across classes, as in HH.
- The cyclotomic μ-coefficient is +1, from the D-arithmetic.
- No (I, μ) reaches 1/5 − 0.01. **Refuted for this model and space.**

**Constructive objective (B)** (exact covers; `experiments/raw/cycloharvey_constructive.jsonl`,
summary `experiments/analysis/cycloharvey_search.json`). The table gives best cost
min_W min_m (m + n_T⌈W/m⌉), with the optimal W and target count in brackets.

| support | n=40 | n=48 | n=56 | n=64 | n=72* | n=80* | slope 56–80 | slope all |
|---|---|---|---|---|---|---|---|---|
| {−1,1} (Lehman) | 168 | 517 | 1555 | 4651 | 14496 | 42656 (W=2^14, 26272 tg) | **0.1997** | 0.1998 |
| {−1,1,2}, {−2,−1,1}, {−1,1,2,3}, {−3,−1,1,3} | identical to {−1,1} | | | | | | 0.1997 | 0.1998 |
| {−2,−1,1,2} | 165 | 518 | 1549 | 4653 | 14480 | 42656 | 0.1998 | 0.2003 |
| **all degree ≤ 3** {−3..3} | 165 | 518 | 1549 | 4653 | 14480 | 42656 | **0.1998** | 0.2003 |
| {1}, {1,2} (BSGS on p) | 1116 | 4466 | 17864 | 71314 | 285888 | 1143424 | 0.2500 | 0.2500 |

\* sampled sizes (pre-registered: 4 sub-intervals, each of relative length 2^{−6}, extrapolated).

- **Every support containing {−1,1} gives the same slope 0.1997–0.2003.**
- LLL, run on the full 6-dimensional degree-≤3 lattice at every location, almost always
  returns a Lehman vector (a q + b p). Its rare higher-degree picks change the cost by
  < 2% and only at n ≤ 72.
- Pre-registered tests:
  - s_full − s_(−1,1) = +0.0002, not < −0.02;
  - no s_I < 0.19;
  - analytic and constructive agree within 0.001 for {−1,1} and {1}.
  Nothing is "promising".
- For the larger supports the constructive slope (0.200) lies *below* the analytic
  3/14 or 2/9. That is expected rather than a discrepancy: the family contains the
  Lehman sub-lattice, and the heuristic e(I) is for vectors with full support I. The
  family exponent is min over sub-supports, which is 1/5.
- **Deviation from the pre-registration.** 8 of 711 tasks were stopped unfinished.
  All of them have support {1} or {1,2} at the smallest W (n ∈ {64, 72, 80}); each
  needs ≥ 10^7 targets. Their missing W values cannot be the minimum, since cost ≥ n_T
  is far above the reported optimum, so the slopes are unaffected.

**Guessed classes** (M ∈ {5, 12, 60}; class r = 1; cost = φ(M) × per-class cost,
with a separate baby list per class).

| support | M | n=48 total (per class) | n=56 total (per class) | unguessed n=48 / n=56 |
|---|---|---|---|---|
| {−1,1} | 5 | 776 (194) | 2320 (580) | 517 / 1555 |
| {−1,1} | 12 | 368 (92) | 992 (248) | 517 / 1555 |
| {−1,1} | 60 | 704 (44) | 1792 (112) | 517 / 1555 |
| all degree ≤ 3 | 12 | 296 (74) | 824 (206) | 518 / 1549 |
| all degree ≤ 3 | 60 | 512 (32) | 1600 (100) | 518 / 1549 |

- Guessing changes constants by factors 0.6–1.5. Two sizes cannot give a slope, and M
  is bounded, so none of this is an exponent effect. The M = 12 gain is in line with
  φ(12)/12 = 1/3, partly offset by per-class baby lists.
- The pre-registered W grid floor 2^{⌈0.12n⌉} was hit for the guessed runs: the
  optimum W sat at the floor. These costs are therefore slight overestimates.
- The cyclotomic version would multiply every entry by at least D = φ(m) (D1). For
  m = 12, D = 4, so it would be ≥ 4× the guessed (Z/N)^* cost. That is 2.4–3.6× the
  unguessed cost at these sizes.

## 5. Step 4 experiments (exact, 32–56 bits; calibration only)

**E1. Operational Harvey baseline** (`operational.harvey_search`, randomised base,
parameters r, m from the model optimum for N alone; 16 balanced semiprimes per size).

| bits | success | STRUCTURAL / ACCIDENTAL | impl. giant list == model count | mean list cost m+giant | mean mults (incl. powering) | mean gcds | null accidental/instance |
|---|---|---|---|---|---|---|---|
| 32 | 16/16 | 16/0 | 16/16 | 240 | 286 | 117 | 0.44 |
| 40 | 16/16 | 16/0 | 16/16 | 967 | 1024 | 480 | 0.45 |
| 48 | 16/16 | 12/4 | 16/16 | 3482 | 108463 | 1640 | 0.35 |
| 56 | 16/16 | 14/2 | 16/16 | 11903 | 447076 | 5433 | 0.26 |

(36, 44 and 52 bits are in `experiments/analysis/cycloharvey_experiments.json`.)

- The implementation's giant-list length equals the model's exact count in **112/112**
  instances. This is an exact calibration of model (c) against an implementation.
- The list-cost slope is 0.225 per bit over 40–56 bits. That is above 1/5 at these
  tiny sizes because the optimum is r = 1 (BSGS on p+q) up to 40 bits, and because of
  log factors. The model count itself reaches 0.20 only above about 80 bits (§1).
- The multiplication count is dominated by the r·lg r handle exponentiations, a
  polylog factor per handle.
- STRUCTURAL means the colliding giant step's target equation i + jm + shift =
  a(N/g) + bg holds for the factor g found. ACCIDENTAL hits (6/112) are below the
  null of 0.26–0.45 per instance from m·n/√N random pair collisions. A hit is
  classified as ACCIDENTAL only when it is the first one found.

**E2. Cyclotomic-handle version** (`operational.cyclo_harvey_search`; m ∈ {5, 8, 12},
all D = 4; all φ(m) = 4 guesses; same (a, b) targets; collisions detected on a random
Z/N-linear projection, with α^i ≡ v (mod p) ⇒ λ(α^i) ≡ λ(v) (mod p); 8 instances per
cell, 32–48 bits).
- It factors 120/120.
- Giant list = φ(m) × the (Z/N)^* list: ratio 2.0–4.0, below 4 only when an exact match
  ends a run early.
- R-multiplications run 2.5–6.5× the (Z/N)^* multiplications. The ratio exceeds φ(m)
  because each cyclotomic handle needs three powerings.
- Base operations (lower bound D per R-mult) run **10–26×** the (Z/N)^* operations.
- ACCIDENTAL hits are 22–62% of successes. The null m·n/√N is 1.0–2.1 per instance,
  because projections collide at rate 1/p and there are φ(m) times more giant steps.
  So the cyclotomic version finds more *chance* collisions but no additional
  structural ones.
- **Experimentally supported:** the cyclotomic version is a constant-factor
  (≥ φ(m)·D) slowdown of the same search, as D1 predicts.

**E3. Handle identity and one-sidedness** (360 instance×m cases at 32/40/48 bits, m ∈ {5, 8, 12}, 4 exponents).

| F | p ≡ q (mod m) | identity mod p | also holds mod q (two-sided) |
|---|---|---|---|
| 3q+2p | no / yes | 278/278, 82/82 | 0 / 0 |
| q³+2p (cubic, asymmetric) | no / yes | 278/278, 82/82 | 0 / 0 |
| p+q | no / yes | 278/278, 82/82 | 0/278, **82/82** |
| p²+q² | no / yes | 278/278, 82/82 | **216/278**, 82/82 |

- The identity T_{r,F} ≡ α^{F(p,q)} (mod p) has 0 failures in 1440 cases.
- Symmetric F at p ≡ q (mod m) is always two-sided, confirming codex.
- **New:** for p² + q² it suffices that p² ≡ q² (mod m). This is automatic for m ∈ {8, 12}
  (every unit squares to 1) and holds when p ≡ ±q for m = 5, giving 29/91.
- In general a symmetric F is two-sided whenever the subgroup generated by the
  exponents k of the σ^k it uses cannot tell p from q.
- None of this affects cost (§2(e)).

## 6. Candidate card (BRIEF §4 template, abbreviated)

- **Mechanism.** One-sided collision search α^u ≡ T_{r,F} (mod p) with cyclotomic
  Frobenius handles, u in a short window of F(p,q). Higher-degree Taylor cancellation
  is supposed to shorten the window.
- **Representation.** R = (Z/N)[X]/Φ_m, D = φ(m), guess r ≡ p (mod m).
- **Factor feature.** u* = F(p,q) in a window. **Extraction.** Product
  tree/Bluestein gcd, or sort-and-match.
- **Outcome.**
  - Refuted for this model: dominated by (Z/N)^* (D1).
  - Exponent ≥ 1/5 − o(1) for all bounded-degree F, with no guess or a
    single-modulus guess (LB, proved, reviewed once).
  - Heuristic exponent table: all higher supports are worse.
  - Constructive search and experiments agree (§§4–5).
- **Scoped obstruction.** In the identity model, adding monomials p^k or q^k with
  |k| ≥ 2 to the target lattice costs a covolume factor (√N)^{|k|}, which outweighs the
  extra Taylor order. Residue guesses gain lattice index M^{K(K+1)/2}, exactly offset
  by the φ(M) runs.
- **What would escape.**
  - Correctness proofs that are not uniform over p, using arithmetic of N itself.
  - Exponents outside Z[p, N/p], such as non-polynomial maps of p.
  - A collision mechanism with sub-sumset cost, which is not possible for
    equality-only detection mod an unknown prime.
  - At the margin: improving the deterministic order-finding step (N^{1/5} now).
    This matters only if the search itself went below 1/5.

## 7. Reproduction

```
.venv/bin/python -m pytest -q tests/test_cycloharvey.py
.venv/bin/python -m src.cycloharvey.run_calibration            # gate (HH sizes 136+ cut; see note)
.venv/bin/python -m src.cycloharvey.run_calibration summarize
.venv/bin/python -m src.cycloharvey.kappa_table
.venv/bin/python -m src.cycloharvey.search analytic
.venv/bin/python -m src.cycloharvey.search constructive --procs 26
.venv/bin/python -m src.cycloharvey.search guessed --procs 26
.venv/bin/python -m src.cycloharvey.search summarize
.venv/bin/python -m src.cycloharvey.experiments harvey
.venv/bin/python -m src.cycloharvey.experiments cyclo
.venv/bin/python -m src.cycloharvey.experiments onesided
.venv/bin/python -m src.cycloharvey.analyze
```
