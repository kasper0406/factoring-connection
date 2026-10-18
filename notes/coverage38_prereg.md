# coverage38 pre-registration (frozen before any registered run)

Agent `coverage38`, 2026-10-18. Code: `src/coverage38/{theory,lehman,oracle,experiments}.py`
(frozen at the commit that adds this file). Raw output: `experiments/raw/coverage38_{e2,e2s,e3}.jsonl`.
Run: `PYTHONPATH=. .venv/bin/python -m src.coverage38.experiments e2|e2s|e3`.

**Disclosed smoke tests before freezing** (different parameters, output to the session
scratchpad, not to experiments/raw): e2 at n = 24 (degenerate: band primes exceed p),
e2s at n = 48 with 20 000 pairs, e3 at n = 40 with 1 q, 60 p, 10 surrogates. They were used
only to debug the code (an index bug in e2s; adding the weight (1,1)).

## Theory being tested (stated in notes/coverage38_results.md)

- **Theorem U** (cost model): any design with J distinct moduli and expected covering count
  λ costs ≥ c·N^{1/5}(λ⁴/J)^{1/5}. So a worst-case gain needs coverage at λ < J^{1/4}.
- **Proposition L** (linear-statistics barrier): configurations drawn uniformly from the
  complement of the target set (zero hits) have all pointwise exponential sums within
  2λ + 2√(J log 8Q) of the uniform model. Prediction: the real configurations
  x_ℓ = p·ℓ̄ mod q are pointwise indistinguishable from these surrogates, and the
  absolute-value Erdős–Turán bound is ≫ λ.
- **Dichotomy** (Lehman catches well-approximable q/p). Asymptotically only O(1) weights
  are "dead" (0 < |δ_v| < Y) even for Lehman-easy pairs, so the dichotomy should not change
  the uncovered rate at P ≈ Y³ (e2s). At the exhaustive sizes P ≈ 8Y–32Y (e2) dead weights
  are not negligible, so a toy-size effect is possible there.

## e2 — exhaustive dichotomy worst case (S_32, S_36)

Band: primes ≡ 3 (mod 4) in (2^13, 2^14] (440), nested prefixes J ∈ {16, 64, 128, 440} (n = 32)
and {16, 128, 440} (n = 36). Signed weights plus (1,1), Hcap = 110. Lehman score
L(p,q) = min_B max(B², δ²/(4√(ABN))); Lehman-easy at gain G ∈ {4, 16}: L ≤ ⌊P^{2/5}/G⌋.

- **D1.** frac_easy(G=16)/frac_easy(G=4) ∈ [0.10, 0.35] at both n (model: 4^{−5/4} ≈ 0.18).
- **D2.** Δ := h*_all − h*_hard ≤ 3 in every (n, J, G) cell (the dichotomy does not lower the
  exhaustive worst case by more than the Poisson extreme-value fluctuation).
- **D3.** Poisson check on the hard population: |h*_hard − h_pred(|hard|)| ≤ 4 in every cell.
- **D4.** Tail enrichment of Lehman-easy pairs (G = 4) among the top 0.1% of h_J lies in
  [0.5, 2.5] in every cell (wide: toy regime).
- Decision: D2 ∧ D3 ⇒ "Lehman dichotomy does not change the exhaustive worst case
  (experimentally supported)". D2 failing ⇒ report the size and the structure of the
  removed worst pairs, and whether the effect is the dead-weight effect (check
  0 < |δ_v| < Y for the removed worst pairs).

## e2s — sampled dichotomy at the D′-type scaling (n = 60, 72; 200 000 pairs each)

Band primes ≡ 3 (mod 4) in (Y, 2Y], Y = 2^{n/6}; signed weights plus (1,1); H_ref = least H
with uncovered fraction ≤ 5%.

- **S1.** R = Pr[uncovered | easy]/Pr[uncovered | hard] ∈ [0.8, 1.25] for G = 4, and
  ∈ [0.7, 1.4] for G = 16, at both n.
- **S2.** frac_easy(16)/frac_easy(4) ∈ [0.10, 0.35].
- Decision: S1 ⇒ "failure sets of Lehman and coverage are not disjoint; no dichotomy gain at
  P ≈ Y³ (experimentally supported)". S1 failing with R < 0.8 ⇒ report a possible dichotomy
  effect (would need theory).

## e3 — the linear-statistics barrier on real configurations (n = 48)

Three primes q (first ≥ 0.83·2^24, then spaced by ~10^6); per q, 1000 random primes p with
pq ∈ S_48; band = all primes in (2^14, 2^15] (J = 1612). Target S = Farey-arc set
{x : ∃B ≤ H, |Bx mods q| ≤ k0}, H ∈ {4, 16}, k0 tuned so λ_q = J|S|/q ≈ 3.
Statistics: E_abs = (1/q)Σ_{h≠0}|Ŝ(h)||K_p(h)| (absolute Erdős–Turán error), M_low =
max_{1≤h≤⌈q/k0⌉}|K_p(h)|, M_all = max_{h≠0}|K_p(h)|; 100 random-complement surrogates per
cell.

- **O1.** E_abs/λ_q ≥ 5 for every real pair (absolute E–T certifies no pair).
- **O2.** median(E_abs, real)/E_pred ∈ [0.8, 1.25], E_pred = (‖Ŝ‖₁ − |S|)/q·(√π/2)√(J(1−J/q)).
- **O3.** Uncovered vs covered real pairs (hits = 0 vs > 0): ratio of medians of E_abs and of
  M_low ∈ [0.9, 1.1] in every cell.
- **O4.** Surrogates vs real: ratio of medians of E_abs, M_low and M_all ∈ [0.8, 1.25] in every
  cell; every surrogate has 0 hits (by construction) and max|K| ≤ 2λ + 2√(J log 8q) in at
  least 50% of surrogates (Prop. L).
- **O5.** Uncovered fraction among real pairs ∈ [0.5, 2]·e^{−λ_q} pooled over the cells of each H.
- Decision: O1–O4 ⇒ "pointwise exponential-sum information does not distinguish uncovered
  from covered pairs at these sizes (experimentally supported); the barrier is realised by
  actual arithmetic configurations."
