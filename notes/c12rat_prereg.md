# c12rat pre-registration (phase 27): toy checks for the rational two-point blocks of C₁/C₂

Agent `c12rat`, 2026-10-24. Code frozen at the commit of this file: `src/c12rat/{geom,bts,e3,experiments}.py`
(analysis script `src/c12rat/analyze.py` may be written after the runs but may only compute the statistics named
here). Toy scale only: nothing here is asymptotic, and no prediction tests the analytic theorems themselves. The
deterministic predictions (marked **[det]**) are exact checks of proved finite lemmas; the others are heuristic.

**Disclosure of smoke runs before freezing.** Two smoke runs of `e3.label_run`, at settings NOT used below:
(a) M = 300, S = 6, label (701, 123), primitive g; (b) M = 1200, S = 8, labels (701, 123) and (997, 5), primitive g
only; (c) M = 500, S = 5, label (503, 77), both laws. Run (b) showed a systematic negative within-cell residual
R ≈ −0.08…−0.12 under the primitive law. I identified it (before run (c)) as the primitivity correlation at the
untyped primes 7 < p ≤ y, a second-order local term that the theory puts into the types for p ≤ P₀; the model
`e3.local_ratios` was written for it, and run (c) gave R = −0.053 against the model −0.051, and R = +0.006 under the
all-g law. The E3 predictions below were written after these smoke runs. Geometry sanity runs (E1-like, 13 labels)
were also made while writing `geom.py`.

## E1 — label shapes and thin pairs (exact; Lemmas LS, MULT, TF of the results note)

Setup: q ∈ {10007, 100003}, ALL r mod q (shape); q = 10007, 1000 uniform r, child box S = 10 (thin pairs),
K ∈ {2, 8, 32}. σ := λ₁²/q, t² := λ₁/λ₂ for the Gauss-reduced basis.
- **P1a [det]** F(σ) := Pr_r[λ₁² ≤ σq] ≤ π(√(σq) + 1)²/(2q) for every q and σ ∈ {1/256, 1/64, 1/16, 1/4}.
- **P1b** F(σ)/σ ∈ [0.85, 1.05] for σ ∈ {1/64, 1/16, 1/4} at q = 100003 (heuristic 3/π ≈ 0.955: one short
  primitive vector per skewed label).
- **P1c [det]** For every label, max_Δ #{unordered pairs with coordinate determinant Δ} ≤ 5S² + S, and the thin
  fraction f(K) := Pr_pairs[κ > K] (κ = |P′||P″|/|det|) is ≤ the Lemma TF bound
  min(1, ⌈2√3S²/(t²K)⌉(10S² + 2S)/(n(n − 1))).
- **P1d** Label-averaged thin fractions decay like 1/K: E_r f(8)/E_r f(2) and E_r f(32)/E_r f(8) both in
  [0.15, 0.45].
- **P1e** Skewed labels (σ < 1/16) carry a share ∈ [0.3, 0.9] of E_r f(32) (thin pairs concentrate in skewed labels,
  but regular labels also contribute).

## E2 — Lemma BT-S (exact inequality)

Setup: y = 50, smooth table to 5·10⁷, 2000 random instances: D uniform in [10, 10⁴], a a random unit mod D,
Z/D log-uniform in [10², 10⁴], X₀ uniform in [10⁶, X − Z), W = Z/D.
- **P2a [det]** #{n ∈ (X₀, X₀ + Z] : n ≡ a (D), P⁺(n) ≤ y} ≤ Σ_{m₁ ∈ (W/y, W] smooth, (m₁, D) = 1}(⌊Z/(m₁D)⌋ + 1)
  in all instances.
- **P2b** (equidistribution sanity) the median of count/((Z/φ(D))·local smooth density) is in [0.8, 1.25].
- Reported only: the median loss rhs/count.

## E3 — cell decomposition of child-box pair covariances (rational block)

Setup: y = 150, g = ALL integer points ("all" law) or all primitive points ("prim" law) of [600, 1200]²
(exact enumeration), child box S = 8, labels q prime in (500, 1000]: 12 with r uniform, 4 with σ < 1/16
(rejection-sampled), seed 20261027. f(n) = 1[|n| is 150-smooth]. Cells of a value: log-interval of ratio 1.5 ×
(v₂ ≤ 3, v₃ ≤ 2, v₅ ≤ 1, v₇ ≤ 1, capped) × (min(v_p, v_p(D)) for the primes 7 < p ≤ y dividing D). For each pair
Pair − S′S″ = Cov_ν(ψ′, ψ″) + R exactly, with ψ′ the one-point cell mean. Per label:
ϕ_ex = Σ_{P′≠P″}(Pair − S′S″)/m², C_ν, B₁ = Σ√(V′V″)/m², R_lab = Σ R/m², and (prim law) the model
R_prim = Σ S′S″(∏_{7<p≤y, p∤D} ratio_p − 1)/m² of the primitivity correlation at untyped primes.
- **P3a [det]** |C_ν| ≤ B₁ for every label and law (Cauchy–Schwarz: the theorem's one-point bound).
- **P3b** all law: |R_lab| ≤ 0.015 for every label (the within-cell residual cancels over the box), skewed labels
  included.
- **P3c** prim law: |R_lab − R_prim| ≤ 0.02 for every label (the residual is the computable second-order
  primitivity term of the untyped primes, which the theory types for p ≤ P₀).
- **P3d** (cancellation) all law: the median over labels of sd_pairs(R/(S′S″))/max(|R_lab|, 10⁻⁴) is ≥ 3.
- **P3e** (heuristic, direction only) all law: mean ϕ_ex(skew) − mean ϕ_ex(uniform) ≥ −0.02 (skewed boxes are
  effectively one-dimensional in g, so their archimedean size correlation is not smaller).

Commands:
```
.venv/bin/python -m src.c12rat.experiments E1 experiments/raw/c12rat_e1.jsonl
.venv/bin/python -m src.c12rat.experiments E2 experiments/raw/c12rat_e2.jsonl
.venv/bin/python -m src.c12rat.experiments E3 experiments/raw/c12rat_e3.jsonl
```
