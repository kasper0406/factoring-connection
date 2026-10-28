# cov38b pre-registration (frozen before any registered run)

Agent `cov38b`, phase 30, 2026-10-28. Code: `src/cov38b/{design,fourier,experiments,analyze}.py`.
Raw output: `experiments/raw/cov38b_e{1,2,3}.jsonl`. Analysis: `experiments/analysis/cov38b_summary.json`.
All arithmetic that decides coverage is exact integer arithmetic; FFTs (E2) are floating point and are used only for
the Fourier statistics they report.

Notation (coverage38 / dchworst): band 𝒫(q) = primes ℓ ≡ 3 (mod 4) in (Y, 2Y] with **Y = ⌊q^{8/25}⌋** (the c = 3/8
scaling of Theorem W, P ≈ Y^{25/8}); J = |𝒫|. Signed primitive weights (A, B), B ≥ 1, A ≠ 0, height max(|A|, B), (1,1)
included. μ_ℓ(t) = least height of a weight with A ≡ tB (mod ℓ). For a ∈ [1, q−1]: h(a) = min_{ℓ∈𝒫} μ_ℓ(a q̄ mod ℓ)
(μ_ℓ(0) = ∞). For a prime p ∈ (q/2, q) this is exactly "the pair (p, q) is covered at height h(p)".

## E1 — worst-case coverage height at the c = 3/8 scaling (exhaustive in p)

- q: for Q0 ∈ {10^5, 10^6, 10^7, 10^8, 10^9}, the three primes q = nextprime(⌊Q0⌋), nextprime(⌊1.3Q0⌋),
  nextprime(⌊1.7Q0⌋) (15 values).
- For every prime p ∈ (q/2, q): h(p). Record h*(q) = max_p h(p), the histogram of h(p), n_p, J, Y.
- "All residues" variant for the 9 values with Q0 ≤ 10^7: every a ∈ [1, q−1]; h*_all(q).
- Model (independent residues): s(H) = ∏_ℓ (1 − ν_ℓ(H)/(ℓ−1)), ν_ℓ(H) = #{t ≠ 0 : μ_ℓ(t) ≤ H} (exact).
  H_α(n) = least H with exp(−n·s(H)) ≥ α. Empirical exponent c_emp = log(√Y/h*)/log J.
- **P1a (proved sanity):** h*(q) ≤ ⌊√ℓ_min⌋ for every q (Thue).
- **P1b:** h*(q) ∈ [H_{0.05}(n_p) − 1, H_{0.95}(n_p) + 1] for at least 12 of the 15 q.
- **P1c:** h*(q) < ⌊√ℓ_min⌋ (strictly below the Thue height) for every q with Q0 ≥ 10^7.
- **P1d:** 0 ≤ h*_all(q) − h*(q) ≤ 3 for all 9 q of the all-residues variant.

## E2 — Erdős–Turán certificates and the Wiener norm of Farey-arc sets

- q = nextprime(10^5), nextprime(10^6), nextprime(10^7); band as above.
- K(h) = Σ_{ℓ∈𝒫} e(h ℓ̄/q) for all h (FFT). κ = max_{h≠0} |K(h)|.
- Farey-arc set S(H′, k₀) = {B̄k mod q : 1 ≤ B ≤ H′, |k| ≤ k₀}. Two nonnegative functions supported on S:
  g_sh = 1_S and g_tri(x) = max_B max(0, 1 − |k|/(k₀+1)) over representations x = B̄k.
  W(g) = (1/q) Σ_{h≠0} |ĝ(h)|; LB(g) = (Σg²)^{3/2}/E(g)^{1/2}, E(g) = (1/q)Σ_h |ĝ(h)|⁴.
- Wiener scan: H′ ∈ {1, 2, 4, 8, 16, 32} with k₀ = ⌊q/(50·H′⁴)⌋ (so 4H′³k₀ < q).
- Certificate scan at the design: for H = 4, 8, 12, …, ≤ 2Y: H′ = ⌊H/4⌋ (≥ 1), k₀ = ⌊qH/(4Y)⌋, and the two shapes
  "farey" (g_tri on S(H′, k₀)) and "arc" (Fejér triangle on the single arc B = 1, |k| ≤ k₀).
  cert(H) := J·Σg/q − κ·W(g) > 0. H_ET = least certified H for each shape.
- **P2a:** κ/√(J ln q) ∈ [0.5, 2] for each q.
- **P2b:** for H′ ∈ {2, …, 32}: W(g_sh)/H′ ∈ [0.2, 5] and W(g_tri)/H′ ∈ [0.2, 5]; the least-squares slope of
  log W(g_tri) on log H′ over H′ ∈ {4, 8, 16, 32} lies in [0.7, 1.15] for each q.
- **P2c:** LB(g) ≤ W(g) + 1 (proved sanity) and LB(g_tri)/H′ ∈ [0.1, 2] for H′ ∈ {2, …, 32}.
- **P2d:** H_ET ≥ ⌊√Y⌋ for both shapes and every q (no phase-discarding certificate below the Thue height),
  while h*(q) from E1 (same q) is < ⌊√Y⌋ for q = nextprime(10^7).

## E3 — the Lehman prong at sub-1/5 budget

- n ∈ {48, 64, 80}; 200 000 pairs of S_n each (`src.dchworst.experiments.sample_pairs`, seed0 = 3838 + n).
- L*(p, q) = min_v max(B², |Bp − Aq|·√(B/P)) (`src.coverage38.explore.lstar_vec`), P = √(pq).
- frac(G) = fraction of pairs with L* ≤ P^{2/5}/G for G ∈ {2, 4, 8, 16, 32, 64}.
- Model (Theorem U with J = 1): frac(G) ≍ G^{−5/4}.
- **P3a:** least-squares slope of log frac(G) on log G over G ∈ {8, 16, 32, 64} lies in [−1.5, −1.0] for each n.
- **P3b:** frac(64) ≤ 0.02 for each n.

## Decision rules

- E1: P1a–P1d all pass ⇒ "at the c = 3/8 scaling the true worst case follows the independent-residue extreme and
  lies below the Thue height (heuristic c → 1/2); no structured obstruction". P1b failing on the *high* side for
  ≥ 4 q ⇒ report a candidate structured obstruction.
- E2: P2a–P2d pass ⇒ "the Wiener norm of Farey-arc minorants grows like H′ (Lemma F), and phase-discarding
  Erdős–Turán certificates do not reach the Thue height at these sizes".
- E3: P3a–P3b pass ⇒ "the Lehman prong at budget P^{2/5}/G covers a fraction ≍ G^{−5/4}, as Theorem U predicts".
- Deviations will be disclosed as such; no prediction will be edited after a run.
