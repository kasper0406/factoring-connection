# abcsunit pre-registration (phase 19)

Agent `abcsunit`, 2026-10-21. Frozen before any registered run. Code: `src/abcsunit/{core,model,run_experiment}.py`;
config `experiments/configs/abcsunit.json`; tests `tests/test_abcsunit.py` (11 pass at freeze).

## Design-stage disclosure

Before freezing, the following (cell, seed) pairs were run during code development. They were used to design the
counting model M1 (the local-factor / shift singular series) and the E2 relation selection (primitive
spanning-forest relations; random combinations of the dependency basis):
(n, κ, K, seed) = (24, 1.25, 9/8, 0), (24, 1.25, 5/4, 0), (24, 1.5, 9/8, 0), (24, 1.5, 5/4, 0).
Observed there: pairs/M1 = 0.70, 0.84, 0.97, 1.02; pairs/M0 = 0.20–0.60; one early E2 design (all chain edges incl.
non-primitive ones, 23 relations for 11 primes) gave rank 8 < 9 and 0/15 fruitful dependencies, which led to the
odd-index observation (Lemma FR in the results note). These four runs are part of the registered grid; they are
scored but flagged "design" in the results, and every prediction is also reported with them excluded.
The last design record is kept in `experiments/raw/abcsunit_design_smoke.jsonl`.

## Instances and grid

- N = `balanced_semiprime(n, seed)` from `src/common/instances.py`, seeds 0, 1, 2.
- y = round((n ln 2)^κ) (depends only on the public bit length); S = primes ≤ y; X = ⌊N^K⌋ (exact integer root).
- Cells (n, κ, K): (24: 1.25 × {9/8, 5/4, 3/2}; 1.5 × {9/8, 5/4, 3/2}; 2.0 × {9/8, 5/4}),
  (28: 1.25 × {9/8, 5/4, 3/2}; 1.5 × {9/8, 5/4}), (32: 1.25 × {9/8, 5/4}; 1.5 × 9/8), (36: 1.25 × 9/8).
  17 cells × 3 seeds = 51 runs. Ψ cap 3·10⁷ (a capped run is reported as skipped, never silently dropped).
- K = 1 is excluded by design: smooth numbers ≤ N are pairwise incongruent mod N (they are distinct residues).
- Memory: one job at a time; `memory_pressure` checked before the run; ≤ 2 GB RSS.

## Quantities

- E1: `pairs` = #{u < v ≤ X, both y-smooth, u ≡ v (mod N)} (exact). Models (floating point, model evaluation only):
  - M0 = (1/|H|)·Σ_k Σ_{u ≤ X−kN} d(u + kN) (equidistribution in H = ⟨S⟩ ⊂ (Z/N)^*, no local factors);
  - **M1** = M0 with each shift kN weighted by the singular series 𝔖_y(k) = Π_{ℓ≤y} σ_ℓ(v_ℓ(k)) · (ℓ > y factors),
    σ_ℓ from Haar-averaged local weights w_ℓ(j) = ℓ^{j(1−α)}(1 − ℓ^{−α})/(1 − 1/ℓ), α = α(X, y) the saddle point
    (definition in `model.py`), u subsampled to 2·10⁴ (seeded);
  - quadratic-character correction factor 1 + Σ_χ (T_χ/Ψ)² over the quadratic characters χ ∈ {χ_p, χ_q, χ_pχ_q}
    nontrivial on H (R1b = pairs/(M1·corr)).
  |H| and the characters are DIAGNOSTIC (they use p, q).
- E2: relations = primitive pairs (gcd(u, v) = 1) forming a spanning forest of each residue class, sorted by v; the
  first R = 3π(y) + 30 are used. Rows = parities of e(u) + e(v). rank over F₂; dependency basis; 40 seeded random
  combinations of the basis; `fruitful` = #combinations with 1 < gcd(W − V, N) < N (STRUCTURAL); z-image dimension
  (DIAGNOSTIC) = F₂-rank of the (z mod p = ±1, z mod q = ±1) sign vectors of the basis dependencies.
  `char_span` = dim span{c_p, c_q}, c_p = indicator of the primes ℓ ≤ y with (ℓ/p) = −1.
- E3: `triply_all` = congruent pairs whose multiplier k = (v − u)/N is also y-smooth; `triply` = the coprime ones
  (solutions of the S′-unit equation u + kN = v, S′ = S ∪ {p, q}). Model: M1 restricted to smooth k (M1_smoothk).
  abc quality q = log v / log rad(u·v·kN) for coprime pairs (exact radicals); `abcv_violations` = coprime pairs with
  q > 1.1 that satisfy NEITHER (i) k·gcd(rad k, N)/rad(k) > N^{0.1}/2 NOR (ii) u/rad(u) > N/2 (Prop ABCV's
  necessary condition, taken with K_ε = 1).
- E4 (N-free): primitive a + b = c, a ≤ b, a, b, c S-units, S = first n primes (n = 1..5), c ≤ H = 2^{16,24,32,40}.

## Predictions (all registered now)

- **P1** (E1, M1). R1 = pairs/M1 ∈ [0.6, 1.5] in every run with M1 ≥ 100; median R1 over those runs ∈ [0.8, 1.2].
- **P1′** (local factors matter). R0 = pairs/M0 < 0.85 in every run with K = 9/8.
- **P1b** (quadratic characters negligible). |R1b − R1| ≤ 0.05 in every run.
- **P2.** |H| = φ(N) in every run.
- **P3** (Prop AB with |H|; identity). Σ_r c_r² ≥ Ψ²/|H| in every run.
- **P4** (E2, rank). rank ≤ π(y) − char_span in every run (theorem); rank = π(y) − char_span in every run where
  ≥ R primitive forest relations exist.
- **P5** (E2, fruitfulness). Over runs reaching the predicted rank: pooled fruitful fraction ∈ [0.40, 0.60]; z-image
  dimension = 2 in every such run; at least one fruitful combination in every such run.
- **P6** (E3). R3 = triply_all/M1_smoothk ∈ [0.6, 1.5] in every run with M1_smoothk ≥ 100 and ⌊X/N⌋ > y (the runs
  where k-smoothness is not automatic).
- **P7** (E3, abc). q_gt1eps/coprime_pairs ≤ 10⁻³ in every run; abcv_violations = 0 in every run (a nonzero value
  means a bug or a wrong proof).
- **P8** (E3, Evertse). triply ≤ 3·7^{2(π(y)+2)+3} (theorem; consistency) and triply ≥ π(y) in every run with
  M1_smoothk ≥ 100 (a fixed factor base supplies more triply-smooth relations than columns already at toy size).
- **P9** (E4). Counts are nondecreasing in H, never exceed the von Känel–Matschke totals (1, 4, 17, 63, 190), and
  equal them at H = 2^40.

## Commands

```
memory_pressure | tail -1
.venv/bin/python -m src.abcsunit.run_experiment E123   # -> experiments/raw/abcsunit_E123.jsonl
.venv/bin/python -m src.abcsunit.run_experiment E4     # -> experiments/raw/abcsunit_E4.jsonl
```

Every failed prediction will be reported as failed; post-hoc analyses will be labelled post hoc.
