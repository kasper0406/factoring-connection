# lv5unif: the relation-supply chain at growing degree, with explicit errors — rigorous L[2/5]

Agent `lv5unif`, phase 12, 2026-10-20.

Files:
- Pre-registration: `notes/lv5unif_prereg.md` (E1, committed before the run).
- Code: `src/lv5unif/` — `kern.py` (exact fibre counts, Lemma Kern-U), `e1_kern.py` (E1), `params.py` (exact
  asymptotic bookkeeping of every margin, the τ table as data, and the cost constant with 2^{2/5} bracketed by
  rationals). Tests: `tests/test_lv5unif.py`.
- Data: `experiments/raw/lv5unif_e1.jsonl`, `experiments/raw/lv5unif_e1_summary.json`.
- Codex: `codex/lv5unif_1.{md,out}` (blind round on items 1–2), `codex/lv5unif_2.{md,out}` (adversarial referee
  of this note; §15).

Base texts. The fixed-degree chain is cited in its cleanest audited form, the preprint `paper2/main.tex`
("[pp]", lemma labels as there), with `notes/lem85_results.md` ("[lem85]") for Lemma 8.5′ and
`notes/recurse_results.md` ("[rec]") for Lemma S″, Lemma F, Lemma OS and Lemma T_rec. Everything here is
**internal and unrefereed**; "proved" means proved from those internal results plus the published theorems listed
in §14.

---

## Verdict

**Yes, conditionally on the project's internal chain: a rigorous Las Vegas L_N[2/5, c] factoring algorithm for every N.**

- **Theorem U [proved, conditional on the internal unrefereed chain; see §10 for the exact dependency list].**
  Fix rationals a, κ, d₀ > 0. The algorithm FACTOR_U[a, κ, d₀] of §1 is Las Vegas and outputs the complete
  prime factorisation of every integer N ≥ 2. Its expected running time is
  L_N[2/5, c + o(1)], where c = c(a, κ, d₀) < ∞ is given by (10.3) in terms of one exponent C_γ = C_γ(a, κ, d₀)
  (the polylogarithmic precision loss of the relation supply, §7).
  - Parameters: d ≍ D^{1/5} (odd), k ≍ D^{2/5}, B = exp(aΛ), where D = log N/log log N and
    Λ = (log N)^{2/5}(log log N)^{3/5}. Special-q primes are L[3/5], norms are N^{θ}, θ ≍ D^{−1/5}.
  - For illustration (exact rational grid search, `params.optimise`): c ≤ 3.52, 4.05, 5.04 and 6.50 for
    C_γ = 0, 1, 5 and 20. We do not compute C_γ numerically: it involves the unspecified absolute constants of
    Hildebrand–Tenenbaum's saddle-point estimates.
  - Hence also L[α, o(1)] for every α > 2/5, and in particular L[3/7].
- **The precision exponent is τ = 0 for every ingredient** (table below). This is the decisive question of
  recurse §3.5, and the answer is that the feared d^{O(d)}, c^d and (log N)^{O(d)} losses are all artifacts:
  - of two wide windows (b ∈ (V/2, V] and the dyadic special-q range, each raised to the power ≈ 2d in LD*);
  - of an exponential constant 2^{2d+5} kept in the norm envelope X_F although it is multiplied by N^{−η₀/2};
  - of a cube-section volume bound in Lemma Kern, which costs ((d+1)H)^{d−2}/H^{d−2}.
  The repairs are a b-window of relative width 1/d, one common box radius V for all special-q primes, the
  envelope X_F := (2d+1)MHV^d, and a dimension-free triangular fibre count (Lemma Kern-U, §3). They cost only
  polynomial factors in d, which are (log N)^{O(1)}.
  - The size model of recurse §3.3 then gives α = (2 + τ)/(5 + 2τ) = **2/5**, attained (not just approached),
    with a finite constant.
- **GPT-6.1-sol's blind round (`codex/lv5unif_1`) reached the same verdict independently**: "the apparent τ = 1
  barrier is removable … the repaired one-level framework reaches 2/5; 3/7 is not forced by precision". It found
  the same three repairs (with a narrow ℓ-window instead of a common V). It also supplied a counterexample showing
  that Lemma Kern at the bare height H = V² really has C_d ≥ (3/2)^d/9, so the height margin H ≥ dV² is needed,
  and a failure mode of Lemma W when η₀ shrinks too fast (η₀ = (log N)^{−2}). Both are respected here.
- **Below 2/5: not reached.** The floor is the one-level descent geometry (split acceptance (log N)^{−2k} forces
  k = O(D^α); the norm exponent θ ≥ c/√k forces k ≳ D^{2(1−2α)}), not precision (§13). Multi-level descent is
  the open problem, as recurse §3.4 says.
- **The shared-rational-side multi-field trick works (Proposition MF, §11) [moment part proved; multi-field
  post-processing proved in outline].** Conditional on the rational data, the fields are independent, and a
  pointwise first-moment bound at one good rational candidate suffices. So the algebraic side needs no second moment
  at all and only exp(o(Λ))·ρ_F precision, at the price of ≈ 1/ρ_F = L[2/5]^{O(1)} fields and one extra height floor
  log H ≥ 13(σ + 1)dΛ log log N (needed so that rejected fields cannot absorb the pointwise mass; codex lv5unif_2 #3).
  It gives a second route to L[2/5] that avoids LD*, Lemma Kern, L2_fam and Lemma W. It does not go below 2/5,
  because it does not change the geometry.
- **LV6 at growing degree [proved, §9].** Stark's factorial is harmless: character primes of
  s = ⌈(2d)!·(log N + log Δ̄ + d)⁴⌉ bits satisfy every requirement of Lemmas E and X, and cost exp(O(d log d)) =
  L[2/5]^{o(1)}, because d log d ≍ D^{1/5}log log N = o(Λ).
- **η₀ bookkeeping [proved, §12].** η₀ := 1/⌈D^{1/2}⌉. Then every N^{η₀}-margin is L_N[1/2, 1/3]-large, which beats
  every L[2/5]^{O(1)}; Lemma W's constant is 24/η₀ = O((log N)^{1/2}); and η₀ = o(θ), so the slack does not move c.
- **Literature [§14].** We found no rigorous L[α, c] with α < 1/2 for factoring arbitrary N. The record we know is
  Lenstra–Pomerance 1992, L[1/2, 1 + o(1)]. Lee–Venkatesan 2018 prove L[1/3] bounds for finding a *congruence of
  squares*, and their factoring theorem assumes their Conjecture 7.1. This is a search result, not a novelty claim.

**The weakest points** (my assessment, shared by codex):
1. Everything rests on the internal fixed-degree chain (hypD, lemmaM, pcest, lv6, ldbound, lem85, recurse), which
   no outside referee has read. Theorem U does not add a new kind of input; it re-proves that chain with explicit
   constants and three design changes.
2. The one external analytic input of substance is still Drappeau's Lemme 5 (Compositio 2015), whose published
   proof refers to §3.3 of Harper's unpublished preprint 1208.5992 ([lem85] §2).
3. Bookkeeping: LD*'s matching of geometric envelopes against smoothness densities is where an error would most
   likely hide. §6 does this matching explicitly.

### τ table (precision exponent at growing (d, k); loss log(1/γ) ≍ d^τ log log N)

| ingredient | where d, k enter as written in the fixed-degree sources | exponential? | repair | τ after repair |
|---|---|---|---|---|
| FM (first moment) | #𝒫 ≍ V²/(d log V); one-point smoothing: TV ≤ 3d²V/H, \|M₀\| ≤ dV^{d−1}; ϱ-argument shifts 2^{O(d)}ϖ² (cost exp(O(d log log N/Λ)) = 1 + o(1)); C_d = 2^{2d+5} only against N^{−η₀/2} | no | b-window and common V for matching with LD*; X_F = (2d+1)MHV^d | **0** |
| Lemma 8.5′ (restricted BV) | none, once restated in a regime (R′) that does not mention L[1/2] | no | regime (R′), §5 | **0** |
| Lemma Kern | cube-section volume ((d+1)H)^{d−2}; (V²/bb′)^{d−2} ≤ 4^{d−2}; codex: C_d ≥ (3/2)^d/9 at H = V² | **yes**, d^{O(d)} | triangular fibre count + H ≥ dV² + b-window: Pr ≤ e/((2H+1)²(bb′)^{d−2}) | **0** |
| L2_fam | (log N)^{O(C₀)} with C₀ ≍ log(1/(ρ_rρ_F))/log B; the old cushion L_N[1/2,1] would force C₀ ≳ D^{1/10} | only via the cushion | cushion B instead of L[1/2] (codex) | **0** |
| LD* (second moment) | Kern's C_d; (V_max/V₋)^{2d} ≤ 2^d; (2^{2d+5})² in Ψ(X_F)²; "≤ 12 partners" (absolute); C_W | **yes**, 2^{O(d)} | common V; b-window; X_F = (2d+1)MHV^d; Kern-U | **0** |
| Lemma W | C_W = 24/η₀ | no, if η₀ log N ≫ Λ + log S | η₀ = D^{−1/2} | **0** |
| Theorem K / S″ | the design's split cost p ≥ (γ/4)(γ/(16 log N))^{2k}/(32k²); no C(2k, k), no 2^ω | k enters only as the (by design) L[2/5]-size split cost | Lemma S″ replaces S′ | **0** |
| Lemma CS (computable S) | the old explicit S lost de Bruijn's second-order term (thmA U5), and S^d must be inside X_F | — | Rankin value at an approximate saddle point, evaluated at the actual X_F(S) | **0** |
| LV6 | Stark ν₁ ≍ 1/((2d)! log Δ̄); SQ's 2^{d−1} sign vectors; t_χ = poly(d, log N) | cost only, exp(O(d log d)) = L[2/5]^{o(1)} | factorial-bit character primes | not a γ-loss |
| Lemma IR, rejections | exp(−log H/(13 d log log H)) = exp(−(log N)^{2/5−o(1)}) | no | — | not a γ-loss |

So α = (2 + 0)/(5 + 0) = 2/5.

---

## 1. Parameters and the algorithm FACTOR_U

**Notation.** For the integer N handed to SPLIT (after preprocessing): λ := log log N, D := log N/λ,
Λ := D^{2/5}λ = (log N)^{2/5}(log log N)^{3/5}, and L_N[2/5, c] = e^{cΛ}. "L^{o(1)}" now means exp(o(Λ)),
uniformly in N. All "C" below depend only on (a, κ, d₀) and on the absolute constants of the cited theorems.
**They never depend on d, k or N.**

**Fixed parameters.** Rationals a, κ, d₀ > 0.

**N-dependent parameters** (all computable in time (log N)^{O(1)}, except S, which costs B(log N)^{O(1)}; §2):
- d := the least odd integer ≥ max(3, d₀D^{1/5}); k := max(⌈κD^{2/5}⌉, d + 1); ν := 1/(4k).
- η₀ := 1/⌈D^{1/2}⌉.
- μ := μ_opt(d, k) of [pp] Lemma `lem:optmu` (exact rational); h := h_*(μ) + η₀ with
  h_*(μ) = max(2ν, μ − (d−2)ν, 1 − dμ + ν); θ := 2μ + h + (d−1)ν.
  Then θ·D^{1/5} → θ₀(d₀, κ) := 3/d₀ if κ ≥ d₀²/4, and 2/d₀ + d₀/(4κ) otherwise (`params.theta0`). (At fixed
  d the k-term of the level regime is lower order: 4ν = 1/k ≍ D^{−2/5}.)
- B := ⌊e^{aΛ}⌋. X := ⌊√(N/2)⌋; j maximal with 2^{kj} ≤ X; Y := 2^{j−1}; J := primes in (Y, 2Y]. So
  log Y = log N/(2k) + O(1) = Λ_{3/5}/(2κ)(1 + o(1)), i.e. ℓ = L_N[3/5, 1/(2κ)].
- M := ⌊N^μ⌋, H := ⌊N^h⌋, W := ⌈N^{μ+ν+η₀/2}⌉ (integer roots of N^p for the rational exponents).
- S: the first S = 2^j in the scan j = 0, 1, 2, … with S²·R(X_r(S))·R(X_F(S)) ≥ 1 (Lemma CS), where
  **V(S) := ⌊√(2Y)·S⌋, common to all ℓ ∈ J**, X_r(S) := 2MV/Y and **X_F(S) := (2d + 1)·MHV^d**.
  Put V := V(S), X_r, X_F at this S, ρ_r := ϱ(X_r), ρ_F := ϱ(X_F) (ϱ(z) = Ψ(z, B)/z).
- Candidates (one set for all ℓ): **𝒫 := {(a, b) ∈ Z² : b prime, (1 − 1/d)V < b ≤ V, |a| ≤ V, gcd(a, b) = 1}.**
- Filter: |G(P)| ≥ MV/32. Counted points, N_ℓ(A), descendable set 𝒟_A, δ_A and γ as in [pp] §3.3, with 𝒫 and V.

The analysis-only parameters (never used by the algorithm): ϖ := B²S² (truncation, FM), the grid
𝒳 := {X_F ϖ^{−2}(3/2)^i} ∩ [X_Fϖ^{−2}, 16dX_F], and R₀ := B^{C₀} with C₀ the least integer such that
B^{C₀} ≥ 64S²B/(ρ_rρ_F) (L2_fam; note the cushion B, not L_N[1/2, 1], following codex lv5unif_1 §1(d)).

**The algorithm FACTOR_U[a, κ, d₀]** is [pp] Algorithm `alg:factor` (FACTOR, with SPLIT and trial division run in
alternation, so there is no N₀ branch) with SPLIT, SAMPLE_FIELD, DESCEND_ONCE, SQ_SEARCH and POSTPROCESS of [pp]
§3.4 and §8, changed as follows:
1. the parameters above, with the common V and the b-window in 𝒫;
2. DESCEND_ONCE's shape test calls **FACTOR_U itself** on |u₀| and |v₀| with fresh coins (recursion, [rec] §1.1);
3. POSTPROCESS uses character primes of s := ⌈(2d)!·(log N + log Δ̄ + d)⁴⌉ bits (§9);
4. the row target M_rows := C_col + 12 as in [pp] (eq:Mrows), with t_χ of §9.

The shape test accepts iff |u₀| and |v₀| are each products of at most k distinct primes, all in J, as in [pp].

*Common V in the special-q search.* In SQ_SEARCH and in [pp] Lemma `lem:reg`, V_ℓ is replaced by V and S by
S_ℓ := V/√ℓ ∈ [S/√2·(1 − o(1)), √2·S]. Lemma `lem:reg` then reads: if λ₁ < √ℓ/(2S_ℓ) = ℓ/(2V), the lattice has no
candidate (the line spacing ℓ/|w₀| exceeds 2V ≥ √2V, and |w₀| < ℓ/(2V) < (1 − 1/d)V < b because ℓ ≤ 2Y < V²/2);
otherwise it has ≤ (8V² + 32λ₂²)/ℓ ≤ 179V²/ℓ ≤ 400S² points in [−V, V]², since λ₂² ≤ (4/3)ℓ²/λ₁² ≤ (16/3)V² and
V² ≤ 2YS². The pair bound of [pp] Lemma `lem:pairs` becomes q = |Δ|/ℓ ≤ 2V²/ℓ ≤ 4S².

**Lemma P (the parameters are consistent) [proved; leading-order asymptotic arithmetic, encoded in
`params.margins`].** For N ≥ N₀(a, κ, d₀), all of the following hold, each with a ratio that is
exp(c(log N)^p(log log N)^q) with c > 0 and (p, q) > (0, 0); the right-hand "C" is any fixed constant. `params.margins`
encodes these comparisons and the further margins named in §§4–6 and 11 (22 checks after codex lv5unif_2 #6); it
compares leading terms only, and is a check of the notes' arithmetic, not an independent proof:
(P1) Y ≥ B^C; (P2) N^{η₀/4} ≥ B^C, ≥ S^C, ≥ 2^{Cd}; (P3) u_F := log X_F/log B ≥ (log log N)^4; (P4) d log d ≤ Λ/C;
(P5) d log S = o(θ log N) and η₀ = o(θ); (P6) log Y ≥ ω(N) log 2 + C log log N, since ω(N) ≤ log N/log B;
(P7) log H/(13d log log H) ≥ C log log N; (P8) μ − ν ≥ η₀ and h − 2ν ≥ η₀ and h − μ + (d−2)ν ≥ η₀;
(P9) 2S² < Y, B < V/2, V < ℓ; (P10) H ≥ dV²·N^{η₀/2}, H ≥ 3dϖ·V, M/V ≥ N^{η₀/2}.
*Proof.* (P8) is the definition of h and the size μ ≈ 1/d ≫ ν + η₀. The others are the listed comparisons of
leading terms: for example, (P1) compares log Y ≍ (log N)^{3/5}(log log N)^{2/5} with log B ≍ (log N)^{2/5}(log log N)^{3/5};
(P2) compares η₀ log N = (log N log log N)^{1/2} with Λ; (P4) compares d log d ≍ (log N)^{1/5}(log log N)^{4/5} with Λ;
(P6) compares log Y with log N/log B ≍ (log N)^{3/5}(log log N)^{−3/5}. ∎

---

## 2. Smooth-number inputs in the new range

Throughout, z is **in range** if B^{D^{1/3}} ≤ z ≤ N. Every argument of ϱ used below is in range: the smallest is
X_r B^{−2C₀−1}, with log X_r/log B ≍ D^{2/5}.

**Lemma U1 (local behaviour of ϱ) [proved].** For N ≥ N₀ and z₁ ≤ z₂ ≤ z₁B^{C} with z₁ in range:
(i) ϱ(z₂) ≤ 2ϱ(z₁); (ii) ϱ(z₁) ≤ (log N)^{C_U(C+1)}ϱ(z₂), with C_U absolute;
(iii) for fixed Z ≥ 1, Ψ(z₁(1 + 1/Z)) − Ψ(z₁) ≥ Ψ(z₁)/(3Z), and Ψ(z₁/2) ≥ Ψ(z₁)/3;
(iv) log(1/ϱ(z)) = u_z log u_z (1 + o(1)), u_z = log z/log B.

*Proof.* As [pp] Lemma `lem:local`, whose proof uses only Hildebrand–Tenenbaum's Theorems 2(i) and 3 and u → ∞:
(1 − α(z, B))log B = log log z − log(1 + log z/B) + O(log log B) ≤ log log N + C log log B, so
B^{1 − α} ≤ (log N)^{1+C}; and 1/u_z ≤ D^{−1/3} → 0 uniformly in range. (iv) is Canfield–Erdős–Pomerance, valid
since u_z ≤ log z/(2 log log z) for z in range (here u_z ≍ D^{2/5} while log z/log log z ≫ D^{4/5} for z ≥ X_r^{1/2}). ∎

**Lemma CS (a computable S with m := S²ρ_rρ_F polylogarithmic) [proved].** For z in range let α̂(z) be the root of
Σ_{p≤B} log p/(p^α − 1) = log z, computed by bisection to absolute accuracy (log N)^{−3}, and
R(z) := z^{α̂−1}∏_{p≤B}(1 − p^{−α̂})^{−1}, computed to relative accuracy 2. Then
  R*(z)/2 ≤ R(z) ≤ 2R*(z), where R*(z) := z^{α̂−1}∏_{p≤B}(1 − p^{−α̂})^{−1} exactly, and
  ϱ(z) ≤ R*(z) ≤ (log N)^{C_R}ϱ(z)   (C_R absolute).
Moreover (existence, proved first, without using the size of S): uniformly for 1 ≤ S² ≤ e^{2σΛ},
  log(R(X_r(S))R(X_F(S))) = −σΛ(1 + o(1)),  σ := 2θ₀/(5a).
Hence the scan of §1 (S = 2^j for j = 0, 1, 2, …, stopping at the first success; there is no binary search and no
monotonicity assumption) stops at some S with S² ≤ e^{(σ+o(1))Λ}; at that S, (log N)^{−2C_R}/4 ≤ m ≤ 64. Computing S
costs O(Λ)·π(B)(log N)^{O(1)} = B^{1+o(1)}.

*Proof.*
- **Upper bound (Rankin).** For every α ∈ (0, 1], Ψ(z, B) ≤ Σ_{P⁺(n)≤B}(z/n)^α = z^α∏_{p≤B}(1 − p^{−α})^{−1}.
- **Lower bound.** Hildebrand–Tenenbaum 1986, Theorem 1: Ψ(z, B) = z^{α}ζ(α, B)/(α√(2πφ₂(α, B)))·(1 + O(1/u + log B/B))
  uniformly for z ≥ B ≥ 2, at the saddle point α = α(z, B), with φ₂(α, B) = Σ_{p≤B} p^α(log p)²/(p^α − 1)²
  ≍ log z·log B. So the Rankin value at α exceeds ϱ(z) by the factor α√(2πφ₂)(1 + o(1)) ≤ log N. The function
  g(α′) := log(z^{α′}ζ(α′, B)) is convex with g′(α) = 0 and g″ = φ₂. For |α′ − α| ≤ 1/log B every term of φ₂ changes by
  a factor ≤ e³, so if |α̂ − α| ≤ (log N)^{−3}, then g(α̂) − g(α) ≤ e³φ₂(α)(log N)^{−6}/2 ≤ 1. Bisection is valid
  because the defining sum is decreasing in α. Hence R*(z) ≤ e·log N·ϱ(z)(1 + o(1)); take C_R := 2 if HT's implied
  constant is absorbed, and in general C_R absolute.
- **Existence (codex lv5unif_2 #2).** For 1 ≤ S² ≤ e^{2σΛ}, d log S ≤ dσΛ = o(θ log N) directly from the growth
  rates d ≍ D^{1/5}, Λ = D^{2/5}λ, θ log N ≍ D^{4/5}λ; so log X_r + log X_F = θ log N(1 + o(1)) uniformly on this range,
  and Lemma U1(iv) with the bounds ϱ ≤ R* ≤ (log N)^{C_R}ϱ gives the displayed asymptotic. At S² = e^{(σ+ε)Λ} the
  product S²R_rR_F is therefore ≥ e^{(ε+o(1))Λ} ≥ 1, so the scan stops at or before that S. This argument uses
  nothing about the S where the scan stops, which removes the circularity of the first draft.
- **m.** At the chosen S, S²R(X_r)R(X_F) ≥ 1, so m = S²ρ_rρ_F ≥ S²R*_rR*_F(log N)^{−2C_R} ≥ S²R_rR_F/(4(log N)^{2C_R})
  ≥ (log N)^{−2C_R}/4. If S ≥ 2, the inequality failed at S/2, so
  S² < 4/(R_r(S/2)R_F(S/2)) ≤ 16/(R*_r(S/2)R*_F(S/2)) ≤ 16/(ρ_r(S/2)ρ_F(S/2)) ≤ 64/(ρ_r(S)ρ_F(S)),
  using Lemma U1(i) (X(S/2) ≤ X(S) ≤ 2^{d+1}X(S/2) ≤ BX(S/2)). So m ≤ 64. (If S = 1, m ≤ 1 trivially.) The value
  of σ: u log u summed over the two sides is (θ₀/a)D^{2/5}·(2/5)λ(1 + o(1)) = σΛ(1 + o(1)).
- **Why S^d must be inside X_F.** S^d = B^{Θ(d)}, and a shift of the argument of ϱ by B^{Θ(d)} changes ϱ by
  (log N)^{Θ(d)} (Lemma U1(ii) with C ≍ d). Evaluating ρ_F at the box size without S would therefore manufacture
  a (log N)^{Θ(d)} = τ = 1 discrepancy (codex lv5unif_1 §1 (2) makes the same point). The scan over S = 2^j evaluates
  R at the actual X_F(S) each time; the map S ↦ S²R(X_r(S))R(X_F(S)) is not assumed monotone. ∎

*Remark.* thmA's explicit S² = ⌈L^{θ/(2a)}⌉ gave only m = L^{o(1)} (thmA U5); at L[1/2] that sufficed. At L[2/5]
the split raises γ to the power 2k ≍ D^{2/5}, so m must be polylogarithmic, and Lemma CS is needed.

---

## 3. Lemma Kern-U (the dimension-free kernel count)

**Lemma Kern-U [proved; elementary].** Let d ≥ 3, and let P = (a, b), P′ = (a′, b′) be primitive with b, b′ ≥ 1 and
Δ := ab′ − a′b ≠ 0. Let Φ: Z^d → Z², c ↦ (⟨c, v(P)⟩, ⟨c, v(P′)⟩), with v(P) = (a^{d−1−i}b^i)_{0≤i<d}, and
C_H = [−H, H]^d ∩ Z^d. Then for every z ∈ Z²,
  #{c ∈ C_H : Φ(c) = z} ≤ (⌊2H/(bb′)⌋ + 1)^{d−2}.
If moreover bb′ ≤ (2H + 1)/d, then Pr_{c∈C_H}[Φ(c) = z] ≤ e/((2H + 1)²(bb′)^{d−2}).

*Proof.*
1. **The kernel.** Identify c with the form R_c = Σ c_i x^{d−1−i}y^i, so ⟨c, v(P)⟩ = R_c(P). [pp] Lemma `lem:kern`(ii)
   (Gauss's lemma) gives ker Φ = Q·Z[x, y]_{d−3} with Q = (bx − ay)(b′x − a′y), and Q is primitive, so
   R ↦ QR is a bijection from Z[x, y]_{d−3} ≅ Z^{d−2} onto ker Φ. Its monomial images e_m := Q·x^{d−3−m}y^m
   (0 ≤ m ≤ d−3) are therefore a Z-basis of ker Φ. In coordinates, e_m has entry bb′ at index m, −(ab′ + a′b) at
   m + 1 and aa′ at m + 2, and zero elsewhere.
2. **Triangular parametrisation.** A nonempty fibre is {c₀ + Σ_m λ_m e_m : λ ∈ Z^{d−2}} ∩ C_H. For i ≤ d − 3,
   coordinate i of c₀ + Σλ_m e_m equals c₀,i + bb′λ_i + (ab′ + a′b)(−λ_{i−1}) + aa′λ_{i−2} (terms with negative
   index absent). Given λ₀, …, λ_{i−1}, the condition |c_i| ≤ H confines λ_i to a real interval of length
   2H/(bb′), which contains at most ⌊2H/(bb′)⌋ + 1 integers. Induction over i = 0, …, d−3 gives the count.
   (The last two coordinates only remove further points.)
3. **Probability.** ⌊2H/(bb′)⌋ + 1 ≤ (2H + bb′)/(bb′), so
   Pr ≤ (2H + bb′)^{d−2}/((bb′)^{d−2}(2H + 1)^d) = (1 + (bb′ − 1)/(2H + 1))^{d−2}/((bb′)^{d−2}(2H + 1)²),
   and (1 + 1/d)^{d−2} ≤ e. ∎

*Remarks.*
- Compared with [pp] Lemma `lem:kern`(iii): no determinant ‖v ∧ v′‖, no cube-section volume, and the constant
  is e for all d. The old proof enlarged the cube to [−(d+1)H, (d+1)H]^d, which costs a factor of order
  ((d+1)/1)^{d−2}·(section constant); that is a d^{O(d)} loss.
- The bound (bb′)^{d−2} ≥ V^{2d−4}(1 − 1/d)^{2d−4} ≥ e^{−8/3}V^{2d−4} needs the b-window (1 − 1/d)V < b. With
  b ∈ (V/2, V] one only gets 4^{−(d−2)}V^{2d−4}: that is the second exponential loss of the old chain.
- **The height margin is necessary (codex lv5unif_1 §1(c)).** At H = b², P = (1, b), P′ = (1, b′) with b′/b → 1, the
  3^{d−2} forms Q·(±1, 0 coefficients) all lie in ker Φ ∩ C_H, so any bound C_d/(H²(bb′)^{d−2}) needs
  C_d ≥ (3/2)^d/9. Lemma Kern-U's hypothesis bb′ ≤ (2H + 1)/d excludes this; in the design H/V² ≥ dN^{η₀/2} (P10).
- E1 (§16) checks the count exactly: 0 violations in 228 pairs and 2 052 fibres, 47 pairs cross-checked by brute
  force; the constant (2H + 1)²(bb′)^{d−2}·max Pr reached at most 1.18 (the lemma allows e ≈ 2.72 with
  bb′ ≤ (2H+1)/d).

---

## 4. Theorem FM-U (the first moment, explicit and uniform)

**Lemma SM1 (one-point smoothing) [proved in [pp] Lemma `lem:smooth1`; uniformity checked].** For P = (a, b) with
gcd(a, b) = 1, ab ≠ 0, V′ := max(|a|, |b|), H ≥ 3dV′: there is a probability measure ϑ on Z with
TV(Law⟨c, v(P)⟩, ϑ) ≤ 3d²V′/(2H + 1), ϑ = Law(Σ Y_i) with Y_i independent uniform on integer intervals of lengths
(2H + 1)l_i (l_i = |a|^{d−1−i}|b|^i), ϑ symmetric unimodal about M₀ with |M₀| ≤ dV′^{d−1}, and
max ϑ ≤ 1/((2H + 1)V′^{d−1}).
*Uniformity check.* Every constant in that proof is explicit and polynomial in d: the kernel basis u_j = be_j − ae_{j+1}
has ℓ^∞-cells of radius |a| + |b|; the d smoothing shifts have ℓ^∞-size ≤ 3V′; the total-variation cost is
d·3dV′/(2H + 1); log-concavity is dimension-free. The exponentially large *cardinality* of the smoothing sets is
neither a probability loss nor a cost (they are analysis devices). This replaces LV Lemmas 5.13, 5.18 and 5.20
(LV §5 is the growing-degree model; we need none of its L^{o(1)}-precision statements). ∎

**Lemma GM-U (algebraic smoothness at a good modulus) [proved].** Let ℓ ∈ J, P = (a, b) ∈ 𝒫, g ∈ 𝒢 and w be such that
t := |G(P)| ∈ [MV/32, 2MV], t is good ([pp] Definition `def:good`, with the grid 𝒳 of §1), and
r₀ := F̂_{g,w}(P) mod t is a unit. Then Pr_c[F_A(P) ∈ 𝒮_B] ≥ ρ_F/(log N)^{C_g}, with C_g = C(a, κ, d₀).

*Proof.* [pp] Lemma `lem:goodmod`, with these explicit inputs.
- V′ ∈ ((1 − 1/d)V, V], so V′^{d−1} ≥ V^{d−1}/e; Y_ϖ := (2H + 1)V′^{d−1}/ϖ; |M₀| ≤ dV′^{d−1} ≤ Y_ϖ/2 by (P10).
- Small layers have ϑ-mass ≤ 3/ϖ.
- **Norm envelope.** By [pp] Lemma `lem:valuesizes`, |F̂(P)| ≤ 2^{2d+5}MHV^dN^{−η₀/2} ≤ MHV^dN^{−η₀/4} by (P2).
  This is the only place the exponential constant of the base expansion occurs, and it is absorbed here.
  Every value F̂ + Gz with |z| ≤ dHV^{d−1} has absolute value ≤ (2d + 1)MHV^d = X_F.
- One large layer: ty′ ≥ (MV/32)(Y_ϖ/2 − 1) ≥ MHV^d/(64eϖ) ≥ 4|F̂| by (P2); ty′/4 ≥ X_Fϖ^{−2} because
  ϖ ≥ 256e(2d + 1); ty′/2 ≤ 16dX_F. So some x₁ ∈ 𝒳 lies in (ty′/4, ty′/2].
- Good t and Tenenbaum's theorem (ω(t) ≤ log N; [pp] Fact `fact:Ten`) give
  Ψ(x₁; t, r′) ≥ Ψ_t(x₁)/(2φ(t)) ≥ Ψ(x₁)/(3t) ≥ (y′/12)ϱ(x₁). Here X_F/B^{C₁} ≤ x₁ ≤ 16dX_F ≤ BX_F, with
  C₁ := ⌈σ/a⌉ + 3 fixed (since X_F/x₁ ≤ 256e(2d + 1)ϖ). If x₁ ≤ X_F, Lemma U1(i) gives ϱ(x₁) ≥ ϱ(X_F)/2; if
  x₁ > X_F, Lemma U1(ii) gives ϱ(x₁) ≥ ϱ(X_F)(log N)^{−2C_U}. Either way ϱ(x₁) ≥ ρ_F/(log N)^{2C_U}.
- The error terms: (|F̂|/t + 1)/y ≤ 33eϖN^{−η₀/4} and 3d²V/(2H + 1) ≤ N^{−η₀}, both ≤ ρ_F·B^{−1} by (P2), (P10).
  So the proportion of smooth values in each large layer is ≥ ρ_F/(24(log N)^{2C_U}) − 33eϖN^{−η₀/4} (each half
  [0, y′], [−y′, 0] contributes y′ϱ(x₁)/12 of the 2y + 1 points), and
  Pr_c ≥ (1 − 3/ϖ)(ρ_F/(24(log N)^{2C_U}) − 33eϖN^{−η₀/4}) − N^{−η₀} ≥ ρ_F/(log N)^{C_g}, C_g := 2C_U + 5. ∎

**Lemma FIB (fibres) [proved in [pp] Lemma `lem:fibres`, with V_ℓ := V].** |𝒢_P(n)| ≤ 2M/V for every n; there is an
interval I_P of length ≥ MV/(16ℓ) − 1 on which ℓn has constant sign, MV/32 ≤ ℓ|n| ≤ 2MV and
|𝒢_P(n)| ≥ M/(CV log log N). The Möbius error is 2^{ω(ℓ|n|)} ≤ exp(O(log N·μ/log log N)) ≤ (M/V)^{o(1)} (codex's
corrected form), which is negligible against M/V ≥ N^{η₀/2}. ∎

**Theorem FM-U [proved].** For N ≥ N₀(a, κ, d₀),
  E_{ℓ∈J}E_{A∼𝒜}N_ℓ(A) ≥ S²ρ_rρ_F/(d·(log N)^{C_FM}),  C_FM = C(a, κ, d₀).

*Proof.* [pp] Theorem `thm:FM`, with:
- The primes of J dividing N are dropped; there are at most log N/log Y ≤ 2k of them, a fraction ≤ 2k/π_J of J
  (Opus proof critic #6; [pp] says O(1) at fixed k).
- #𝒫 ≥ V·(V/(2d log V)): primes in ((1 − 1/d)V, V] number ≥ V/(2d log V) by the prime number theorem with
  de la Vallée Poussin's error term, since 1/d ≥ exp(−c√(log V)) (log V ≍ (log N)^{3/5}(log log N)^{2/5}); each b
  has ≥ V coprime a with |a| ≤ V.
- Good mass: for each P, Σ_{n∈I_P∩𝒮_B}|𝒢_P(n)| ≥ cM²ρ_r/(ℓ log log N), by Lemma U1(iii) and FIB. With
  ℓ ≤ 2Y and V² ≥ YS² this gives a contribution ≥ cS²ρ_rρ_F/(d log V log log N (log N)^{C_g}).
- Bad mass: the moduli ℓ|n| lie in 𝒯_Q for the 7 dyadic Q covering [MV/32, 2MV] (one range for all ℓ, by
  the common V). Lemma 8.5′(R′) (§5) applies to each Q and each x ∈ 𝒳 (by (P8): q₀ − 2ν = μ − ν ≥ η₀ and
  ξ/2 − q₀ ≥ η₀/2, ξ := μ + h + dν), and [pp] Corollary `cor:bad` gives
  π_J^{−1}Σ_ℓ 𝔟_ℓ ≤ C|𝒳|ε(N)(MV/Y)ρ_r(log N)^{C}, with |𝒳| = O(log ϖ) = O(Λ) and
  ε(N) = exp(−c u/log²u), u ≍ D^{2/5}. This is ≤ (log N)^{−C−5} times the good mass. ∎

*Losses (all polynomial in d or polylogarithmic):* 1/d and 1/log V from 𝒫, log log N from φ(m)/m, and
(log N)^{C_g} from Lemma GM-U. Compare [pp] Remark after `thm:FM`.

---

## 5. Lemma 8.5′ in regime (R′)

The hypothesis "y = ⌊L_N[1/2, a]⌋, log x/log N in a compact set" of [lem85] §1 is used in its proof only through
the following consequences. We take them as the hypotheses.

**Regime (R′)** (quantitative gaps; no fixed λ — repair after codex lv5unif_2 #1). Fix K₀ ≥ 1, and let U₀ ≥ 8 be an
absolute constant such that the error term of Hildebrand–Tenenbaum's Theorem 3 is ≤ 1/100 for u ≥ U₀ − 1. Let y ≥ 2
and x, Y, Q satisfy:
- (R′1) log y ≥ 4(log log x)², (log x)^{c_D} ≤ y;
- (R′2) u := log x/log y ≥ max((log log x)⁴, 4(2K₀ + 6)/η_D);
- (R′3) u log(4u log y) ≤ K₀ log y;
- (R′4) Y ≥ y^{K}, K := 2K₀ + 6, and the **level gap** Γ := min(log(Q/Y), log(√x/Q)) satisfies
  Γ ≥ (2K₀ + U₀)log y + 4 log log x.

**Lemma 8.5′(R′) [proved; re-check of [lem85]].** Under (R′), for x ≥ x₀(K₀),
  Σ_{t∈𝒯_Q} max_{(r,t)=1}|Ψ(x, y; t, r) − Ψ_t(x, y)/φ(t)| ≤ C·Ψ(x, y)·δ_𝒯·ε,
  ε := exp(−c u/log²(u + 1)) + y^{−1/2},
and Qδ_𝒯/10 ≤ #𝒯_Q ≤ 3Qδ_𝒯 (δ_𝒯 = ϱ(Q/Y, y)/log Y), where C, c, x₀ depend only on K₀ and on absolute constants
(Drappeau's η_D, c_D, δ_D, C_D; HT's; C_dBT, C_α).

*Proof (every use of regime (R) and of λ in [lem85], replaced).*
- (α1)/(α2) ([lem85] §2): 1 − β(z, y) ≤ (log(u log y) + C_α log log y)/log y ≤ 1/4 for 2y ≤ z ≤ x, by (R′1).
- E7: ϱ(z, y) ≥ y^{−1}exp(−u_z log(4u_z log y)) ≥ y^{−1−K₀} for 2y ≤ z ≤ x, by (R′3) (u_z ≤ u). This replaces
  "u log(4u log y) = (ξ/a²)(1 + o(1))log y".
- Lemma SM (a): Drappeau's Lemme 5 at z = x/d ≥ x^{1/2} (d ≤ Q′ ≤ √x) needs (log z)^{c_D} ≤ y ≤ z and
  R₁ = y^K ≤ z^{η_D}; the latter is K log y ≤ (η_D/2)log x, i.e. u ≥ 2K/η_D, which is (R′2). The transfer factor is
  R₁^{1−β′} ≤ (log x)^{K(1+C_α)} by (α2). The saving exp(−δ_D(u/2)/log²(u/2 + 1)) beats (log x)^{K(1+C_α)+1}
  because u/log²u ≥ (log log x)⁴/(16(log log log x)²) by (R′2).
- Lemma SM (b) and Part 2: Corollary LS′ gives ≤ 1.5·10³π(y)(log x)²(y^{3/2}x/R + y^{1/2}x^{1/2}Q₀) with R = y^K and
  Q₀ = Q′ (resp. R = Y and Q₀ = Q). The target is Ψ(x, y)ϱ(Q′, y) (resp. Ψ(x, y)δ_𝒯), and by E7 both are
  ≥ x·y^{−2−2K₀}/log Y. Relative to the target, the first term is ≤ C(log x)³y^{2K₀+9/2}/R ≤ y^{−1/2} because
  K = 2K₀ + 6 and Y ≥ y^K, and the second is ≤ C(log x)³y^{2K₀+7/2}·Q₀/√x ≤ C(log x)³y^{2K₀+7/2}e^{−Γ} ≤ y^{−1/2}
  by (R′4) (this is the uniform saving codex lv5unif_2 #1 asks for).
- Part 1: Lemma SM applies to Q′ = Q/ℓ, since 2y ≤ e^{Γ}/2 ≤ Q/(2Y) ≤ Q′ ≤ Q ≤ √x·e^{−Γ}.
- (i), the count of 𝒯_Q: HT Theorem 3 at z = Q/(2ℓ) ≥ e^{Γ}/4 ≥ y^{U₀−1}, so u_z ≥ U₀ − 1, and α(z, y) ≥ 0.99 by (α2)
  and (R′1); the constants 0.8 and 4.1 of [lem85] follow.
- E2 (lift averaging), E3, E4, E5, E6, E8, Lemma LS, Perron: regime-free. ∎

**In the design.** y = B, x ∈ 𝒳 (so log x = θ log N(1 + o(1)) and u ≍ D^{2/5}), Y as in §1, Q one of the dyadic
values covering [MV/32, 2MV]. (R′1): log B ≍ Λ ≫ (log log N)². (R′2): u ≍ D^{2/5} ≫ (log log N)⁴.
(R′3): u log(4u log B)/log B → (4/5)θ_F₀/a², so K₀ := 1 + ⌈θ₀/a²⌉ works. (R′4): log Y ≍ Λ_{3/5} ≫ K log B; and, for
every x ∈ 𝒳 (log x ≥ (μ + h + dν)log N − 2 log ϖ) and every such Q,
  log(Q/Y) ≥ (μ − ν)log N − O(Λ) ≥ η₀ log N − O(Λ),
  log(√x/Q) ≥ ((h − μ + (d−2)ν)/2)log N − O(Λ) ≥ (η₀/2)log N − O(Λ)
by the floors (P8). So Γ ≥ (η₀/2)log N − O(Λ), which exceeds (2K₀ + U₀)log B + 4 log log N by (P2). No step uses a
fixed λ; the gap grows like (log N log log N)^{1/2}, and every constant depends only on K₀ = K₀(a, κ, d₀).
Codex lv5unif_1 §1(b) notes that B = L_x[1/2, Θ(1)]; codex lv5unif_2 #1 notes that the earlier fixed-λ formulation did
not justify the application with λ_N = η₀/(4ξ) → 0, which this version removes. ∎

---

## 6. Theorem LD*-U (the second moment) and its inputs

**Lemma W-U [proved; [pp] Lemma `lem:W`].** For every 1 ≤ t ≤ 2MV and r, Pr_w[w ≡ r (mod t)] ≤ C_W/φ(t) with
C_W := 4 log W/log(W/t) ≤ 24/η₀ = O(D^{1/2}).
*Proof.* Brun–Titchmarsh (Montgomery–Vaughan) as in [pp]. log(W/t) ≥ (η₀/2)log N − log(4S) ≥ (η₀/3)log N by
(P2); log W ≤ 2 log N. Codex's counterexample (η₀ = (log N)^{−2}, where t can exceed W/S^{1/2}) does not apply,
because η₀ log N ≫ log S. ∎

**Lemma DS-U [proved; [pp] Lemma `lem:DS`].** For fixed C₀, X in range, 1 ≤ q ≤ B^{C₀}, and an interval
I ⊂ [XB^{−C₀}, X] of length A ≥ XB^{−C₀}: #{n ∈ I : n ≡ r (q), P⁺(n) ≤ B} ≤ (A/q)ϱ(X)(log N)^{C_U(2C₀+2)+2}.
*Proof.* [pp]; the only non-elementary step is Lemma U1(ii) with the shift B^{2C₀+1}. ∎

**Theorem L2-U (= L2_fam) [proved; [pp] Theorem `thm:L2`].** With 𝒮* := 𝒮_B ∩ [−X_F, X_F] ∖ {0},
X′₊ := 2MV/Y, 𝒩_ℓ := {n ∈ 𝒮_B : MV/32 ≤ ℓ|n| ≤ 2MV}, and R₀ = B^{C₀} ≥ 64S²B/(ρ_rρ_F), uniformly for
1 ≤ q ≤ 4S²:
  E_{ℓ∈J}Σ_{n∈𝒩_ℓ} 𝔔(ℓ²|n|q)/φ(ℓ|n|) ≤ τ(q)(log N)^{C_L}·Ψ(X′₊)Ψ(X_F)²/(Yq(MV)²).
*Proof and uniformity.* [pp]'s proof, line by line:
- X_F, hence f(m′), is common to all ℓ (common V).
- Class (A): DS-U, (log N)^{O(C₀)}, with C₀ ≤ 2σ/a + 3 fixed (because S² and 1/(ρ_rρ_F) are B^{O(1)}); Σ1/m₁ ≤ τ(q)(1 + log R₀).
- Class (B): ratios q/(R₀ρ_rρ_F) ≤ 1/(16B) and X′₊q²/(X_Fρ_rρ_F) = N^{−h−(d+1)ν}B^{O(1)}d^{O(1)} = o(1).
- Class (C): ratios log N/(π_Jρ_F) ≤ N^{−2ν}B^{O(1)} and X′₊q²Y⁴/(π_JX_Fρ_rρ_F) ≤ N^{(5−d)ν−h}B^{O(1)}; both o(1)
  by (P1), (P8).
- (4.1): φ(ℓ|n|) ≥ ℓ|n|/(C log log N).
- **The cushion.** [pp] fixes R₀ ≥ 64S²L/(ρ_rρ_F) with L = L_N[1/2, 1]. At L[2/5] that would force
  C₀ ≥ log L/log B ≍ D^{1/10}, and DS's (log N)^{O(C₀)} would become exp(Θ(D^{1/10}log log N)): a spurious τ > 0
  (codex lv5unif_1 §1(d)). Class (B) needs only q/(R₀ρ_rρ_F) ≤ 1, so the cushion B suffices. ∎

**Theorem LD*-U [proved].** For N ≥ N₀(a, κ, d₀),
  E_{ℓ∈J}E_{A∼𝒜}[N_ℓ(A)(N_ℓ(A) − 1)] ≤ d·η₀^{−1}(log N)^{C_LD}·S⁴ρ_r²ρ_F²,  C_LD = C(a, κ, d₀).

*Proof.* [pp] Theorem `thm:LD`, Steps 0–7, with these changes and explicit constants.
- **Pairs.** For P ≠ P′ in 𝒫 ∩ Λ_ℓ(g): ℓ | Δ, q := |Δ|/ℓ ≤ 2V²/ℓ ≤ 4S²; for fixed P and q at most 10 partners
  P′ (two lines parallel to P, spacing |P| ≥ b > 2V/3, segments of length ≤ 2√2V). [pp] Lemma `lem:pairs`.
- **Step 1 (kernel).** Lemma Kern-U with bb′ ≤ V² ≤ (2H + 1)/d (P10) and (bb′)^{d−2} ≥ e^{−8/3}V^{2d−4} (b-window):
  Pr_c[Φ(c) = z] ≤ e^{11/3}/((2H + 1)²V^{2d−4}).
- **Steps 2–3.** Unchanged: Z ≤ ½(𝔮(P) + 𝔮(P′)); 𝔮(P) ≤ 𝔔₀(t, |Δ|; r₀); all values in [−X_F, X_F] (Lemma GM-U's
  envelope).
- **Step 4.** Lemma W-U: E_w𝔔₀ ≤ (C_W/φ(t))𝔔(tℓq).
- **Step 5.** For fixed P, P′, n: #{g : G_g(P) = ℓn, Y_{P′}(g) = 1} ≤ U := (3M/V)ρ_r(log N)^{C}, by DS-U with
  modulus q ≤ 4S² ≤ B^{C₀} and X_rB^{−C₀}/q ≤ M/V (because X_r/(M/V) = 2V²/Y ≤ 4S²).
- **Step 7 (assembly; the exact powers).** Since the n-sum does not depend on P,
  E[N(N−1)] ≤ (#𝒫/|𝒢|)·10·(e^{11/3}C_W/((2H+1)²V^{2d−4}))·U·Σ_{q≤4S²}τ(q)(log N)^{C_L}Ψ(X′₊)Ψ(X_F)²/(Yq(MV)²).
  Insert #𝒫 ≤ 4V²/(d log V) (Brun–Titchmarsh on the b-window), |𝒢| ≥ M²/8, Ψ(X′₊) = X_rρ_r (X′₊ = X_r),
  Ψ(X_F) = (2d + 1)MHV^dρ_F, and Σ_q τ(q)/q ≤ (1 + log 4S²)². (FM-U's ρ_r is ϱ at the I_P scale ∈ [X_r/128, X_r],
  within a factor 2 of ϱ(X_r) by Lemma U1(i).)
  The powers of M, H and V cancel exactly:
    V²·M^{−2}·H^{−2}V^{−(2d−4)}·(M/V)·(MV/Y)·(MHV^d)²/(Y(MV)²) = V⁴/Y² ≤ 4S⁴,
  and the d-dependent constants are (1/d)·(2d + 1)² ≤ 9d. **There is no ratio V_max/V₋** (common V) and **no
  2^{O(d)}** (Kern-U, the envelope, the b-window). ∎

*Where the old chain lost exponentially, explicitly.* Old Step 1: C_d of the cube section (d^{O(d)}) and
(V/b)^{2d−4} ≤ 4^{d−2}. Old Step 7: (V_max/V₋)^{2d−4+…} ≤ 2^{d+O(1)} and (C′_d)² = (2^{2d+5} + 2d)² in Ψ(X_F)². None of
these can be compensated by the smoothness densities, because the corresponding argument shifts change ϱ only by
1 + o(1) (codex lv5unif_1 §1(e) makes the same observation).

---

## 7. Relation supply: γ ≥ (log N)^{−C_γ}

**Theorem D-U [proved].** γ = Pr_{ℓ,A}[N_ℓ(A) ≥ 1] ≥ (log N)^{−C_γ} with C_γ = C(a, κ, d₀).
*Proof.* Paley–Zygmund as in [pp] Theorem `thm:Dk`: with m = S²ρ_rρ_F, P₁ := d(log N)^{C_FM} and
P₂ := dη₀^{−1}(log N)^{C_LD}, γ ≥ (E N)²/(E[N(N−1)] + E N) ≥ 1/(P₁²P₂ + P₁/m). By Lemma CS, m ≥ (log N)^{−2C_R}; and
d ≤ log N, η₀^{−1} ≤ (log N)^{1/2}. So P₁²P₂ ≤ (log N)^{2C_FM+C_LD+4}, P₁/m ≤ 4(log N)^{C_FM+2C_R+1}, and
C_γ := 2C_FM + C_LD + 2C_R + 5 works. ∎

**Field rejections [proved; [pp] Lemma `lem:IR`, `lem:reject`].** Lemma IR's proof is uniform in d with the
absolute prime threshold r₀ = 16 (the irreducible count ((r^d − 1)/(r − 1) − 2r^{d/2})/d ≥ r^{d−1}/(2d) holds for
all d ≥ 3, r ≥ 16), and the near-uniformity factor (1 + Q_IR/(2H + 1))^d ≤ e needs Q_IR ≤ (2H+1)/d, true since
Q_IR = H^{1/2+o(1)}. So Pr(𝓡) ≤ exp(−log H/(13d log log H)) + 3/(2H + 1) ≤ γ/2 by (P7). ∎

---

## 8. Split counting and field selection at growing k

**Theorem K-U [proved].** Let A be a field returned by SAMPLE_FIELD with δ_A ≥ γ/4. Then one call of
DESCEND_ONCE(A) returns a row or makes a factor exit with probability ≥ (δ_A/(4 log N))^{2k}/(32k²).
*Proof.* |𝒟_A| ≥ γπ_J/4 ≥ 4k³2^{ω(N)} by (P6) and Theorem D-U. The products of k distinct primes of 𝒟_A lie in
(Y^k, (2Y)^k] ⊂ (X2^{−2k}, X], so some window [Y′/2, Y′] of the 2k dyadic windows holds a set 𝒮 of
≥ C(|𝒟_A|, k)/(2k) of them, with 2Y′² < N. Lemma S″ ([rec] §3.2) gives Pr_R[E] ≥ |𝒮|²/(2φ(N)) for R uniform on
Q_N, with E: the primitive vector of Λ_R is ±(u, v), u, v ∈ 𝒮 coprime. So Pr_r[E] ≥ |𝒮|²/(2N). With
|𝒮| ≥ (|𝒟_A|/k)^k/(2k), |𝒟_A| ≥ δ_AY/(2 log Y), Y^{2k} ≥ N/(4·16^k) and 2k log Y ≤ log N, this is
≥ (δ_A/(4 log N))^{2k}/(32k²). On E the shape test passes (exactly k distinct primes of 𝒟_A ⊂ J on each side) and
every special-q search returns an atom ([pp] Prop `prop:BS`). ∎

*Remarks.* No 2^{−ω(N)} (Theorem S′'s loss; 2^ω can be L_N[3/5]^{Θ(1/log log N)}, far above L[2/5]) and no
C(2k, k) (S′'s c_k). The shape-test waste ratio of Theorem K(iii) is not used: every attempt is charged a full
shape test and 2k special-q searches ([rec] §3.2, codex recurse_2).

**Proposition select-U [proved; sharper than [pp] Prop `prop:select`(ii)].** Let p be the probability that one trial
(SAMPLE_FIELD + DESCEND_ONCE) ends in a factor exit or a row, C_trial a uniform bound on the expected cost of a
trial or of one DESCEND_ONCE, and T_post that of POSTPROCESS. Then
(i) p ≥ (γ/4)(γ/(16 log N))^{2k}/(32k²); (ii) E[cost of SPLIT] ≤ (2M_rows + 1)C_trial/p + 2T_post.
*Proof.* (i): [pp] Prop `prop:select`(i) with Theorem K-U and the set 𝒢* = {A ∉ 𝓡 : δ_A ≥ γ/4} of measure ≥ γ/4.
(ii): selection costs ≤ C_trial/p (Wald). The kept field A* has law ς_A𝒜(dA)/b, and given A* = A the expected cost
is ≤ 2(M_rowsC_trial/ς_A + T_post) (≤ 2 expected rounds, Theorem DD). Averaging,
(b/p)·2(M_rowsC_trial E[1/ς_{A*}] + T_post) ≤ 2M_rowsC_trial/p + 2T_post, since E[1/ς_{A*}] ≤ 1/b and b ≤ p.
([pp] divides T_post by p too; at L[1/2] that was harmless, here it would add s to the 2a term.) ∎

---

## 9. LV6 at growing degree

[pp] §8 (= lv6 Prop LV6, transferred by ldbound §6.1) is stated for fixed d. Its constants, re-checked:

- **Lemma C, rows, columns:** degree-free. Computing the roots of f mod r for all r ≤ B costs π(B)·poly(d, log B).
- **Lemma H:** log₂|H| ≤ d log₂|D| + d + ½log₂|disc g| + (d−1)log₂(d + log₂|disc g|). With ‖f‖_∞ ≤ 3MH
  (Lemma GB_w and (P2)): log|D| = log|f_d disc f| = O(d(μ + h)log N + d log d) = O(log N), since d(μ + h) = O(dθ) = O(1);
  and disc g = f_d^{(d−1)(d−2)}disc f gives log|disc g| = O(d²(μ + h)log N + d log d) = O(d log N) (codex lv5unif_2 #5).
  So log|H| = O(d log N) and t_χ := ⌈5d(ln|H|_bd + 2 ln log N)⌉ = O(d² log N). Polynomial.
- **Lemma E (≤ 1 exceptional class per field):** valid for every odd d. Stark's Lemma 8 with c(L) = 4(2d)! for
  L_h = K(√h) (non-normal over Q); Landau with conductors ≤ Δ̄^{1/d}.
- **The factorial threshold and its repair [proved].** ν₁ = min(1/(4(2d)! log Δ̄), c_L/(2 log Δ̄)). Lemma X needs the
  exceptional term 2x^{1−ν₁} to be ≤ x/(40 log x), i.e. ν₁ log x ≥ log(80 log x). The old choice
  log x ≈ (log N)² fails as soon as (2d)! ≫ log N. Take
    s := log x := ⌈(2d)!·(log N + log Δ̄ + d)⁴⌉ (codex's choice; ours was 8(2d)!(log Δ̄)² + (log N)²).
  Then ν₁s ≥ (log N + log Δ̄)³/(4 log Δ̄) ≫ log(80s) = O(d log d + log log N), and LO's conditions
  s ≥ 20d(log Δ̄)² (log x ≥ 10n_L(log d_L)², n_L = 2d) and c₁x·exp(−c₂s^{1/2}/√(2d)) ≤ x/(40 log x) hold. The base-field
  sampler ([pp] "Sampler") needs ζ_K to have no zero in [1 − 1/(4·d!·log|Δ_K|), 1) beyond what s absorbs; the same s
  works. So Lemma X holds verbatim: Pr_𝔭[χ_𝔭(h) = −1] ≥ 1/(5d) for every h ∉ E ∪ {0}.
- **Cost of the repair.** Each character prime has s = exp(O(d log d))·(log N)^{O(1)} bits. Sampling (AKS,
  Cantor–Zassenhaus root finding) and each Legendre symbol cost s^{O(1)}. In total
  t_χ·(expected O(d s) draws + 2kM_rows evaluations)·s^{O(1)} = M_rows·exp(O(d log d))(log N)^{O(1)} = L^{a+o(1)}, by (P4).
- **Lemma SQ:** 2^{τ−1} ≤ 2^{d−1} sign vectors, each a Newton lift with O(d²) multiplications of log A-bit numbers,
  log A = O(k M_rows d log N + d² log(dMH)). Cost exp(O(d))·L^{a+o(1)}.
- **Linear algebra:** Wiedemann/Kaltofen–Saunders, M_rows²·poly(d, log N); verification exact; three tests
  s₁, s₂, s₁ + s₂ (|H₀| ≤ 2).
- **(T):** unchanged ([pp] Remark `rem:T`); the recursive shape test keeps it ([rec] Lemma T_rec).

**Proposition LV6-U [proved].** For every field returned by SAMPLE_FIELD, POSTPROCESS is Las Vegas with expected cost
T_post ≤ M_rows²(log N)^{O(1)}·exp(O(d log d)) + B^{1+o(1)} = L^{2a+o(1)}, and its output satisfies (T). ∎

*Polynomial-bit alternative (not used).* Codex lv5unif_1 §3: Thorner–Zaman, Algebra & Number Theory 13 (2019),
Thm 1.1/1.4, gives a lower bound relative to Li(x) − θ₁Li(x^{β₁}); with polynomial s the detection probability is
≫ min(1, ν₁s)/d, so factorially many characters are needed, again exp(O(d log d)). We have not checked those
statements against the paper. **UNVERIFIED**, and not needed.

---

## 10. Recursion at Λ_{2/5} and Theorem U

**Lemma F_α [proved; [rec] Lemma F with √2 replaced by 2^α].** Let 0 < α < 1 and L_x := exp(Λ_α(x)) with
Λ_α(x) := (log x)^α(max(1, log log x))^{1−α}. If F: [2, ∞) → [0, ∞) is nondecreasing and finite and, for x ≥ x₀,
F(x) ≤ L_x^{p+ε(x)} + L_x^{δ+ε(x)}F(√(2x)) with ε(x) → 0, then F(x) ≤ L_x^{c+o(1)}, c := max(p, δ/(1 − 2^{−α})).
*Proof.* Scaling: for y = √(2x), (log y)^α ≤ 2^{−α}(log x)^α + α2^{−α}(log 2)(log x)^{α−1} and log log y ≤ log log x, so
Λ_α(y) ≤ 2^{−α}Λ_α(x) + o(1) and L_y^b ≤ e·L_x^{2^{−α}b} for x ≥ x₁(b). Then [rec]'s steps 2–5 verbatim with
η′ := η(1 − 2^{−α})/4: in step 4 the second exponent is δ + η′ + (c + η)2^{−α} ≤ c + η2^{−α} + η′ = c + η − 3η′. ∎

**The recursion inequality.** Let F(x) := sup_{2≤n≤x} E[cost of FACTOR_U(n)]. Finiteness holds by the alternation
with trial division ([rec] §1.1 item 3). For x ≥ x₀, by [rec]'s accounting (factor tree with ≤ 2log₂x nodes,
Lemma W of [rec] for adaptive sub-calls, Lemma OS for oracle substitution) and Proposition select-U:
  F(x) ≤ 2log₂x·[P(x) + D(x)·F(√(2x))],
  P(x) = L_x^{max(2a, a + s + max(a, σ)) + o(1)},  D(x) = 2(2M_rows + 1)/p = L_x^{a + s + o(1)},   (10.1)
with s := 2κ(1 + C_γ), because:
- M_rows ≤ (d + 1)π(B) + O(d² log N) = L^{a+o(1)};
- 1/p ≤ 4γ^{−1}(16 log N/γ)^{2k}·32k² = exp((2k + 1)(log 16 + (1 + C_γ)log log N) + O(log k)) = L^{s+o(1)},
  since 2k log log N = 2κΛ(1 + o(1));
- C_trial = 2F(√(2N)) [shape test] + 2k·T_SQ + (log N)^{O(1)}poly(d), with T_SQ = Õ(B + S²) = L^{max(a,σ)+o(1)}
  ([pp] Prop `prop:BS`; the batch has ≤ 200S²·2 numbers of ≤ log X_F = O(log N) bits);
- T_post = L^{2a+o(1)} (Prop LV6-U); preprocessing L^{a+o(1)}; S costs B^{1+o(1)} (Lemma CS).

**Theorem U [proved, conditional as stated in the Verdict].** For fixed rationals a, κ, d₀ > 0, FACTOR_U[a, κ, d₀] is
Las Vegas and completely factors every N ≥ 2 in expected time L_N[2/5, c + o(1)] with
  c = max(2a, a + s + max(a, σ), (a + s)/(1 − 2^{−2/5})),  s = 2κ(1 + C_γ),  σ = 2θ₀(d₀, κ)/(5a).   (10.3)
*Proof.* Lemma F_α with α = 2/5 applied to (10.1); the 2log₂x factors are L^{o(1)}. Correctness: AKS-certified primes;
almost-sure termination from the finite expectation. ∎

**Dependencies of Theorem U** (all internal results are unrefereed):
- internal: Theorem DD (lemmaM §1), Lemma S″, Lemmas OS, W, F, T_rec (recurse), Lemma 8.5′'s proof (lem85), [pp]
  §§3, 5–9 (= hypD, pcest, ldbound, lv6, thmA, audited at fixed d), and this note's Lemmas P, U1, CS, Kern-U, GM-U,
  W-U, DS-U, L2-U, Theorems FM-U, LD*-U, D-U, K-U, Prop select-U, Prop LV6-U, Lemma F_α;
- published: Hildebrand–Tenenbaum 1986 (Thms 1, 2(i), 3), Canfield–Erdős–Pomerance 1983, Tenenbaum's Ψ_q theorem,
  Drappeau 2015 Lemme 1(i) and Lemme 5, Iwaniec–Kowalski Thm 7.13, Montgomery–Vaughan 1973 (large sieve,
  Brun–Titchmarsh), the prime number theorem with error term, Lagarias–Odlyzko 1977, Stark 1974 Lemma 8,
  Landau's theorem, AKS, Wiedemann/Kaltofen–Saunders, Cantor–Zassenhaus, Bernstein's batch smooth parts.

**Exact cost arithmetic (`params.py`, rationals only).** 1/(1 − 2^{−2/5}) ∈ [4.12981, 4.12983] (bracketed by
checking p⁵ against 4q⁵ exactly). The optimum of (10.3) over a rational grid, by C_γ: 3.517 (C_γ = 0;
a = 3/5, κ = 1/8, d₀ = 1), 4.047 (C_γ = 1), 5.038 (C_γ = 5), 6.498 (C_γ = 20). At each optimum the special-q term
a + s + max(a, σ) and the recursion term (a + s)/(1 − 2^{−2/5}) nearly balance (for C_γ = 0: 3.5167 against 3.5103, so
the special-q term binds; for C_γ = 5 the recursion term binds); linear algebra (2a) never binds. (Corrected after
codex lv5unif_2 #7.) These numbers illustrate the shape of c; C_γ itself is not
computed (it contains HT's implied constants), and C_γ also depends on (a, κ, d₀) through C₀ and C₁.

---

## 11. The shared-rational-side multi-field trick

**Design MF.** One rational side g (as before) and J_f independent algebraic sides (w_j, c_j), j ≤ J_f, giving
F_j := F̂_{g,w_j} + G·R_{c_j}, all with F_j(g₀, g₁) = Nw_j and the common root m = g₀/g₁ mod N. (w_j = 1, i.e. design
(a″_k), also works: Lemma W is not needed.) A prime ℓ is MF-descendable if some P ∈ 𝒫 ∩ Λ_ℓ(g) has Y_P = 1
(rational, filtered) and F_j(P) ∈ 𝒮_B for some j.

**Design MF, with an extra height floor.** As above, and additionally H ≥ exp(13(σ + 1)·d·Λ·log log N) (that is,
h ≥ h_*(μ) + η₀ + 13(σ + 1)dΛ log log N/log N; the extra term is ≍ (log N)^{−2/5}(log log N)^{7/5} = o(θ), so it does
not move the exponent). This is needed because the pointwise bound below concerns raw fields, and post-processing needs
accepted ones (codex lv5unif_2 #3: in the kernel-binding regime Lemma IR's rejection bound exp(−log H/(13d log log H))
= exp(−Θ(D^{2/5})) is far *larger* than ρ_F = exp(−Θ(D^{2/5}log log N)), so without the extra floor all the
pointwise smoothness mass could sit in rejected fields).

**Proposition MF [moment part proved; multi-field post-processing proved in outline].** Let
γ_r := Pr_{g,ℓ}[∃P ∈ 𝒫 ∩ Λ_ℓ(g) with Y_P = 1 and |G(P)| good]. If J_f ≥ 4(log N)^{C_g}/ρ_F, then
  γ_MF := Pr[ℓ is MF-descendable through an *accepted* field] ≥ γ_r/2,  and  γ_r ≥ (log N)^{−C}.
*Proof.*
- **Decoupling.** Condition on (g, ℓ) and on the event that some P has Y_P = 1 with t = |G(P)| good; let P₀ be the
  first such P in a fixed order (a function of (g, ℓ)). By Lemma GM-U, for every w_j,
  Pr_{c_j}[F_j(P₀) ∈ 𝒮_B] ≥ ρ_F/(log N)^{C_g} =: p₁ (its residue c_{P₀}w_j is a unit, [pp] Lemma `lem:coprime`). The
  rejection event of field j depends on (g, w_j, c_j) only, and by Lemma IR (stated for fixed g, w, with c uniform) and
  [pp] Lemma `lem:reject`, Pr_{c_j}[rejected | g, w_j] ≤ exp(−log H/(13d log log H)) + 3/(2H + 1) + 1/W ≤ p₁/2 by
  the extra floor. So Pr_{c_j}[F_j(P₀) ∈ 𝒮_B and field j accepted | g, ℓ, w_j] ≥ p₁/2. (A field that makes a factor
  exit ends the algorithm successfully.) The pairs (w_j, c_j) are independent of each other and of (g, ℓ), so
  Pr[no accepted j with F_j(P₀) smooth | g, ℓ] ≤ (1 − p₁/2)^{J_f} ≤ e^{−2} < 1/2.
- **The rational moment.** N_r^{good} := #{P : Y_P = 1, t good}. E N_r^{good} ≥ E N_r/2 ≥ S²ρ_r/(log N)^C (the
  rational part of FM-U, with the bad-mass bound of §4), and E[N_r^{good}(N_r^{good} − 1)] ≤ E[N_r(N_r − 1)] ≤
  (log N)^C S⁴ρ_r² ([pp] Prop `prop:R` = Theorem R; it is rational-only and uniform: Lemma `lem:reg` gives ≤ 200S²
  partners, and its divisor switching is DS-U). With S² ≥ (log N)^{−C}/ρ_r (the rational-only analogue of Lemma CS),
  Paley–Zygmund gives γ_r ≥ (log N)^{−C}. ∎
- **Why no algebraic second moment is needed.** The one-point bound is a statement about one fixed candidate, and
  independence across fields turns it into a probability. Codex lv5unif_1 §2 gives a different valid argument
  (conditional PZ on the full count, needing an averaged algebraic second moment with exp(o(Λ)) loss, and only
  exp(o(Λ)) fields). Its counterexample — fields that all succeed on the same deterministic target set — is
  excluded here by the pointwise lower bound p₁ > 0 at every (g, ℓ) with a good rational candidate.
- **Post-processing with J_f fields (outline, repaired after codex lv5unif_2 #4).**
  - Columns: the shared rational ones (sign, p ≤ B), and for each field j its parity, (r, s) and character columns.
    Per field, t_χ^{(j)} := ⌈5d(ln|H_j|_bd + 2 ln log N + ln(4J_f))⌉ characters, so that, by Lemma X and a union bound,
    with probability ≥ 3/4 every field has |H_{j,0}| ≤ 2 simultaneously. Rows: M_rows ≥ C_col^{MF} + J_f + 12, so the
    kernel has dimension ≥ J_f + 12.
  - Linearity: φ_j is F₂-linear on the kernel K of all columns (Lemma H's φ, restricted), with values in H_{j,0}. On
    the event |H_{j,0}| ≤ 2, H_{j,0} ⊂ F₂, so b_j(s) := [Γ_{j,s} is not a square] = φ_j(s) is an F₂-linear functional
    on K.
  - Take J_f + 2 independent uniform kernel vectors s_i (independent with probability ≥ 1 − 2^{−9}); compute the
    J_f × (J_f + 2) bit matrix b_j(s_i) with Algorithm SQ in each field; take a nonzero x in its right kernel (it
    exists, J_f equations in J_f + 2 unknowns) and S := Σx_is_i ≠ 0. Then φ_j(S) = 0 for every j, i.e. every field's
    product is a square. **Verify exactly** with SQ in every field; on failure (characters bad for some field),
    resample characters and kernel vectors. Expected O(1) repetitions; Las Vegas.
  - Congruence: Y := √q_S·∏_j t_j^{−E_j(S)/2}ψ(η_j)(ψ(γ⁻_{j,S})ψ(g_j′(β_j)))^{−1}, with one shared ψ: x ↦ m. (T) holds:
    every quantity is a function of certificates and coins.
  - Cost: J_f(J_f + 2) SQ calls and (J_f dB)² linear algebra; all exp(O(Λ)).
- **Cost and exponent.** J_f = L^{σ_F + o(1)}, σ_F := log(1/ρ_F)/Λ. Per special-q: ≤ 200S_r² rational values with
  S_r² ≈ 1/ρ_r, then J_f algebraic values for each of ≈ S_r²ρ_r rationally smooth candidates. All exponents are
  finite, so Theorem U's conclusion α = 2/5 also follows from Proposition MF, with a different constant, once the
  outlined post-processing is written in full. **The trick does not lower α**: (R1) and (R3) of §13 are unchanged.

So the answer to the task's item 2: yes, conditioning on g decouples the algebraic sides; the conditional moment
needed is only the one-candidate first moment Lemma GM-U (provable, proved), and the algebraic precision can be as
weak as ρ_F·exp(−o(Λ)) (the loss only increases J_f). It is an alternative route to 2/5, not a route below it.

---

## 12. η₀(N) → 0 and Lemma W

- η₀ := 1/⌈D^{1/2}⌉. Every margin that the fixed-degree chain wrote as N^{cη₀} becomes exp(cD^{−1/2}log N) =
  exp(c(log N log log N)^{1/2}) = L_N[1/2, c]: it beats every L[2/5]^{O(1)}, every 2^{O(d)} and every S^{O(1)}
  (Lemma P (P2)).
- Lemma W-U's constant C_W ≤ 24/η₀ = O((log N)^{1/2}) is polylogarithmic; it enters P₂ once.
- θ = θ_* + η₀ + O(dν·0) with θ_* ≍ D^{−1/5} ≫ η₀ = D^{−1/2}, so the slack changes σ by o(1).
- Failure mode avoided (codex): if η₀ log N ≤ log S, then W ≤ t for some moduli t ≤ 2MV and Lemma W fails; (P2)
  excludes this. Codex's alternative η₀ ≍ 1/d also works (then C_W = O(d)).
- No other η₀-dependence: the level and window floors use η₀ only as an additive slack, and Lemma 8.5′(R′) uses it
  through the level gap Γ ≥ (η₀/2)log N − O(Λ) (no fixed λ; §5).

---

## 13. Why not below 2/5 in this design

These are statements about this algorithm family (one-level special-q descent of Dixon-split primes), not lower
bounds for factoring.
- (R1) [heuristic upper bound on the split] The shape test accepts only if |u₀|, |v₀| ≤ √(2N) are products of ≤ k
  primes of size N^{1/(2k)}; the density of such integers is about (C/log N)^k each, so the expected number of
  attempts per row is ≥ (log N/C)^{2k}, and L[α] cost forces k = O(D^α).
- (R3) [proved for the design] θ = 2μ + h + (d−1)ν with h ≥ 2ν and the window/level floors gives θ ≥ max(3/(d+1),
  (d+1)ν) − o(1), so θ² ≥ 3/(4k)(1 − o(1)); and the smoothness cost of norms N^{θ} at B = L[α] is
  exp((1 − α)θ log N·log log N/log B·(1+o(1))) by CEP, which is L[α]^{O(1)} only if θ = O(D^{2α−1}).
- Together k ≳ D^{2(1−2α)} and k = O(D^α), so α ≥ 2/5. Theorem U attains it.
- Below 2/5 one needs multi-level descent and a lower bound for smooth values over a field-dependent good set
  (endogenous children; recurse §3.4). Agent `nested` (phase 12) works on that under ERH. Nothing here helps.

---

## 14. Literature (has anyone proved a rigorous L[α < 1/2] for factoring?)

[Filled in after the Opus literature critic; see §15.]

---

## 15. Reviews

All findings were adjudicated; every valid fix is applied in the text above.

**GPT-6.1-sol blind round `codex/lv5unif_1` (xhigh; items 1–2 of the task, before my draft existed).** Independent
verdict: τ = 0 after repairs; α = 2/5 reached; 3/7 not forced by precision; the shared-side trick works but needs more
than conditional independence; factorial-bit characters are a valid LV6 repair. Its specific contributions, adopted:
the Kern counterexample C_d ≥ (3/2)^d/9 at H = V² (so the height margin is needed); the Lemma W failure mode for
η₀ = (log N)^{−2}; the L2_fam cushion (B instead of L_N[1/2,1]); S^d inside X_F in both moments; the corrected Möbius
error 2^{ω(t)} ≤ exp(O(log t/log log t)); the multi-field post-processing repair (bit matrix over J_f + 2 kernel
vectors). Its mixed-moment version of the MF amplification is recorded in §11 as an alternative. It used a narrow
ℓ-window where we use a common V; both work.

**GPT-6.1-sol adversarial referee `codex/lv5unif_2` (xhigh): ACCEPT WITH FIXES for Theorem U** ("I found no surviving
exponential loss in the displayed relative moment calculation"; it confirmed Kern-U, the GM-U/FM-U envelope, the exact
cancellation V⁴/Y² ≤ 4S⁴ in LD*-U Step 7, L2-U/W-U/DS-U, K-U/select-U, LV6-U and Lemma F_α). Findings:

| # | severity | finding | adjudication |
|---|---|---|---|
| 1 | MAJOR | Lemma 8.5′(R′) fixed λ, applied with λ_N → 0; "x^{−λ}y^{2K₀+7}(log x)² ≤ x^{−λ/2}" did not follow | accepted: §5 rewritten with the quantitative level gap Γ, no λ; constants depend on K₀ only |
| 2 | MAJOR | Lemma CS: size of S proved circularly through (P5); "binary search" with no monotonicity | accepted: existence proved first on S² ≤ e^{2σΛ}; linear scan over powers of 2 |
| 3 | MAJOR | Prop MF: pointwise bound is for raw fields; Lemma IR's rejection exp(−Θ(D^{2/5})) can exceed ρ_F | accepted: extra height floor log H ≥ 13(σ+1)dΛ log log N (o(θ log N)); pointwise bound now for accepted fields |
| 4 | MAJOR | MF post-processing: character union bound over J_f fields; kernel dimension; linearity of the square indicator | accepted: per-field ln(4J_f) characters, rows ≥ C_col + J_f + 12, linearity on |H_{j,0}| ≤ 2, final exact SQ verification; still labelled outline |
| 5 | MINOR | log|disc g| = O(d log N), not O(log N) | accepted (§9); conclusions unchanged |
| 6 | MINOR | `params.margins` not a complete certificate | accepted: 8 further checks added (22 in all), wording narrowed |
| 7 | MINOR | "recursion term binds in every case" false | accepted: the two terms balance at the optimum (§10) |

**Opus proof critic (read-only, independent; reviewed the pre-fix snapshot): ACCEPT WITH FIXES, no FATAL.** It
re-derived every constant of Kern-U, CS, U1, GM-U, FM-U, LD*-U Step 7 (power count M⁰H⁰V⁴Y^{−2}, d-constant ≤ 9d),
L2-U, W-U, D-U, K-U, select-U, (10.1)/(10.3), Lemma F_α and LV6-U, and ran its own brute-force Kern-U check (300 cases,
0 violations, maximum 1.18). Findings: MF rejection across fields (MAJOR; same as codex #3; fixed as there); Lemma 8.5′(R′)
λ → 0, too-weak (R′4), zero slack with λ = η₀/(2ξ) (same as codex #1, and my own pre-review fix to η₀/(4ξ); now
superseded by the gap formulation); binding term (codex #7); disc g (codex #5); common V in Lemma `lem:reg` (fixed in §1:
≤ 400S² points, q ≤ 4S²); primes of J dividing N number ≤ 2k (stated in FM-U); stale strings in `params.TAU_TABLE` and
GM-U's per-layer constant 1/24 (fixed). Its weakest point: §5, now re-proved in gap form.

**Opus literature critic.** See §14.

---

## 16. Experiment E1 (pre-registered, `notes/lv5unif_prereg.md`, commit 6c8b1d8)

Command: `.venv/bin/python -m src.lv5unif.e1_kern experiments/raw/lv5unif_e1.jsonl experiments/raw/lv5unif_e1_summary.json`
(210 s).

| prediction | result |
|---|---|
| P1: 0 violations of the triangular bound | **0/228 pairs (2 052 fibres)**; 47 brute-force pairs (d = 3, 4) also cross-check the lattice enumeration exactly (0 failures) |
| P2: median central/bound ≥ 0.05 on window pairs with H = dV² | **holds** where testable: d = 3: 1.0, d = 4: 0.85, d = 5: 0.86, d = 7: 0.62; d = 6 and 9 had no such pair (disclosed) |
| P3: max (2H+1)²(bb′)^{d−2}·max Pr ≤ 6.6 over H ≥ dV² | **holds**: maximum 1.18 (d = 5) |

Disclosures: 12 (d, V, H) cells were skipped because every drawn pair had a bound above the 3·10⁶ enumeration cap
(2 441 redraws in total); at d = 9 only H < dV² cells ran. As the pre-registration says, in lattice mode P1 is a
check of the parametrisation, and it is an independent test of the inequality only in brute-force mode. P2 shows the
bound is tight to within a factor ≈ 1.6 in the design's regime. E1 tests an exact combinatorial lemma only.

## 17. Candidate card

- **Mechanism.** Theorem A's relation supply (special-q first and second moments, Dixon descent, LV6) at growing
  degree, with a recursive shape test.
- **Model.** Design (a‴_k; μ) with d ≍ D^{1/5}, k ≍ D^{2/5}, a b-window of relative width 1/d, one box radius for all
  special-q primes, B = L[2/5, a].
- **Result.** τ = 0 for every ingredient; rigorous L[2/5, c] (conditional on the internal chain); also via the
  multi-field route.
- **Failure scope.** Below 2/5 the one-level geometry fails (§13), whatever the precision.
- **Next step.** Nested descendability for multi-level descent (L[1/3] in the size model).
