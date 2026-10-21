# siconstruct: pre-registration (phase 20, 2026-10-21)

Agent `siconstruct`. Written and committed before any registered run. Code: `src/siconstruct/` (core, constructors,
theory, run_exp). Raw output: `experiments/raw/siconstruct_E{1,2,3}.jsonl`. Commands:
`.venv/bin/python -m src.siconstruct.run_exp E1` (then E2, E3), one job at a time, memory ≤ 2 GB.

## Setting

Toy Prop SI path with the o(1) terms dropped (`core.lpath_params`): y = round(exp(κ F_{1/3})),
h = round(exp(ν F_{2/3})), κ = (2/9)^{1/3}, ν = 1.25·2/(3κ) = 1.376. N = `balanced_semiprime(bits, seed)`.

**Caveat, stated in advance.** At toy sizes this path gives h ≈ N^{0.57–0.71}, so N/h is small (2^{9} at 32 bits,
2^{34} at 80 bits). Every method that pays (N/h)^{1/2} looks cheap here. The experiments test *laws* (hit counts,
cost scaling, independence of cofactors), not asymptotic competitiveness.

## E1: the MITM hit law and amortisation over k

- bits = 48; seeds 0, 1, 2; K ∈ {1, 4, 16, 64, 256}.
- A: 4096 distinct N-oblivious y-smooth numbers in [2^20, 2^21]; B: 4096 in [B0, 2B0], B0 = ⌈KN/2^20⌉; products lie
  in [KN, 4KN], so about 3K multipliers k are served by one pair of lists.
- `mitm_construct` finds every pair with a·b ∈ [kN, kN + h] for some k (exact; checked against brute force in tests).
- Model (heuristic H_eq: a·b mod N uniform for N-oblivious lists): hits = |A||B|(h+1)/N ≈ 150 per cell,
  **independent of K**.
- **P1a.** Pooled Σhits/Σpred over the 15 cells lies in [0.85, 1.15].
- **P1b.** Every cell ratio lies in [0.6, 1.5].
- **P1c.** No amortisation gain: the per-K pooled ratios (over seeds) have max/min ≤ 1.4.

## E2: per-output cost along the toy path, MITM versus brute force

- bits ∈ {32, 40, 48, 56, 64, 72, 80}; seeds 0, 1; target R0 = 64 outputs.
- MITM: |A| = |B| = M = ⌊√(R0 N/(h+1))⌋ + 1, A ⊂ [A0, 2A0], A0 = 2^{⌈bits/2⌉}, B ⊂ [B0, 2B0], B0 = ⌈N/A0⌉.
  Cost = |A| + |B| + (a, k) queries + hits.
- Brute force: sieve [kN, kN + W) for k = 1..8, W = min(h+1, 2^26), measured for bits ≤ 56. The Hildebrand–Tenenbaum
  saddle-point local density (`ht_local_density`) is the model for all bits.
- **P2a.** Pooled hits/pred over the 14 cells in [0.85, 1.15]; every cell in [0.6, 1.5].
- **P2b.** MITM measured cost per distinct output lies within a factor 2 of (2M + 3M)/R0 ≈ 2.5·2√(N/(R0(h+1)))
  (3 queries per a, since products lie in [N, 4N]).
- **P2c.** For bits ≤ 56 the measured brute-force density is within a factor 1.5 of the H-T model.
- **P2d (direction).** MITM is cheaper per output than brute force at every tested size (model prediction). The
  model crossover with R0 = 64 is near 310 bits (computed before the run: log2 per-output brute 93.0 vs MITM 94.2 at
  320 bits). Asymptotically MITM costs (N/(hR))^{1/2} = exp(Θ(x)) per output, brute force L[2/3, 2/(3κ)].

## E3: the 2-D aim (one constructed factor, tested cofactor; Lemma RC)

- bits ∈ {40, 48, 56}; seeds 0, 1; K = 16; S: 64 N-oblivious y-smooth s ∈ [S0, 2S0], S0 = 2^{e}, e ∈ {⌊lg h⌋/2,
  ⌊lg h⌋ − 4} (lg h = bit length of h minus 1). For each (s, k), at most 2048 offsets w ≡ −kN (mod s), w ≤ h;
  the cofactor m = (kN + w)/s is tested.
- **P3a (model).** Distinct outputs / Σ_tests δ_HT(m) lies in [0.6, 1.6] pooled over seeds per (bits, e) cell, i.e. the
  cofactor is smooth like a random integer of its size (no gain from the aim beyond shrinking the tested part).
- **P3b (proved, checked).** Every tested cofactor is ≥ kN/(2S0), so its size is never below N/(K h) when useful
  (s ≤ Kh).
- **P3c.** Tests per output ≈ 1/δ_HT(N/S0)·(1 ± 0.6); the gain over brute force equals δ(N/S0)/δ(N).

No other experiment is registered. Any extra run will be reported as exploratory.
