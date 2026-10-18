# below16 pre-registration (agent `below16`, 2026-10-18)

Committed before any pre-registered run. Code: `src/below16/{avgbound,multibase,qs_pilot,experiments}.py`.
Tests: `tests/test_below16.py` (11 pass). Config: `experiments/configs/below16.json`
(= `CONFIG` in `experiments.py`). Raw output: `experiments/raw/below16_e{1,2,4}.jsonl`;
summary `experiments/analysis/below16_summary.json`.

**Disclosure of what was run before this file.**
- The unit tests (n ≤ 20 bound checks, 40/48-bit QS determinism and the Jacobi rule) — all pass.
- A timing-only job for E1 components at n = 24, 26 (outputs not inspected beyond run time).
- Two QS smoke jobs (10 instances at 40–64 bits, seeds unrelated to E4's). All 10 factored.
  The per-instance fractions of nontrivial kernel basis vectors I saw were
  72/185, 51/105, 32/56 (40 bits, an earlier version that sieved the whole interval),
  339/658, 196/402, 301/559 (48 bits, same version), and 17/32, 20/40, 35/66, 27/50
  (56/64 bits, with the frozen stopping rule). The E4 band below was set knowing these.

## Questions

1. **Theory (no experiment can prove it; E1 checks the proof's inequalities exactly).** An
   average-case / true-factor analogue of multimod Thm A for *local* designs. Each target's
   class, weights and window depend on N only through N mod K (T1, arbitrary finite windows),
   or through (N mod K, ⌊N/Δ⌋) with Δ ≥ P·K (T2, interval windows), K ≤ C·M². Claimed
   bound: the covered fraction of S_n is ≤ C n² P^{−1/2} E_N[Σ_τ |S_τ|^{1/2}], hence
   ≤ C n² X^{3/2} P^{−1/2} under the sumset cost model.
2. **Multi-base (escape 1).** With rational bases, identity-provable collisions live in
   prime-exponent coordinates. Every coordinate gives a 1-D target with a finite (sum-of-APs)
   window, whose coverage contains the multi-base target's coverage.
3. **Relation methods (escape 4).** A deterministic, uniform relation method (single-polynomial
   QS, exact F₂ algebra) needs a *sign-independence* hypothesis for a derandomized proof. Does
   the data look like independent fair signs?

## E1 — exact bound checks (exhaustive S_n)

n ∈ {20, 22, 24, 26, 28, 30}. Generators:
- (a) DCL (H = 3, band (Y, 2Y], both roots, n-only whole-range windows, K = ℓ and K = ℓ²);
- (b) 150 random T1 generators per n (random A, B, M, K | M², D ≤ M², side; random class sets;
  random windows: interval / AP / 2-dim GAP / random set);
- (c) 40 "oracle-local" T1 generators per n (per key, the class and the s ∈ {1, 4, 16, 64}
  most frequent window values: an adversarial local design built from the population);
- (d) 8 oracle T2 generators per n ∈ {20, 22, 24, 26} (per (key, box), the best class and
  best interval of width W ∈ {1, 4, 16, 64}).

Predictions (all are consequences of the proofs, so any violation means a proof or code error):
- **P1a.** count ≤ core bound (the exact Lemma-1/Lemma-3 counting bound) for 100% of generators.
- **P1b.** core ≤ closed form for 100% of generators that satisfy the stated conditions
  (M ≤ P/2, K ≤ P/2, K ≤ 2^{n−2}).
- **P1c.** count ≤ closed form for 100% of those generators.
- Descriptive only: the tightness ratios count/closed by family.

## E2 — multi-base reduction (n = 22, exhaustive)

Base sets {2,3,6}, {2,4,8}, {2,3,12,18}, {6,10,15}, {2,9,12}; 60 random targets each
(forms B_i p + A_i q, |A_i|, |B_i| ≤ 3, boxes with sides ≤ 6, centred on a random true pair).
- **P2a.** For every prime coordinate j, coverage(multi-base) ⊆ coverage(coordinate j):
  0 failures.
- **P2b.** |π_j(E^T box)| ≤ |E^T box| ≤ |box| always.

## E4 — deterministic QS statistics

n ∈ {40, 48, 56, 64, 72, 80}. Per n: 400 Blum pairs (p ≡ q ≡ 3 mod 4) and 200 general pairs,
uniform from S_n (rejection sampling, seed `20261020|n`). At n = 80: 200 + 100.
Parameters `qs_pilot.default_config(n)`, a function of n only: y = 2^9 … 2^14, X0 = 2^13 … 2^20,
24 extra relations. The algorithm sees N only.
- **P4a.** The Jacobi rule [gcd proper ⇔ (XY/N) = −1] holds for 100% of the kernel basis vectors
  of Blum instances. (This is a proved classical fact, so the prediction is a code check.)
- **P4b.** The pooled per-basis-vector nontrivial rate lies in [0.45, 0.55] in every (n, family)
  cell.
- **P4c.** For each prefix length j ∈ {1, 2, 3, 4, 5, 6, 8, 10}, the number of instances whose
  first j basis vectors are all trivial is consistent with Binomial(m, 2^{−j}). Criterion: an
  exact two-sided binomial p-value ≥ 0.001 in every cell.
- **P4d.** Success (a proper factor) in ≥ (cell size − 1) instances of every cell.
- **P4e.** No ACCIDENTAL (factor-base trial-division) hits: both factors exceed y.

**Decision rule.** If P4b and P4c hold, the sign-independence hypothesis needed by the open
sub-question of escape 4 is *experimentally supported* (40–80 bits, this deterministic
sequence). A failure of P4b or P4c is reported as a candidate structural obstruction, with no
post-hoc rescue. P4a/P1/P2 failures are bugs; they will be fixed and disclosed, and the run
repeated.
