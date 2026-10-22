# zroute pre-registration (E1–E3)

Agent `zroute`, phase 22, 2026-10-22. Written and committed before any run of the registered experiments.

**Disclosed pilots.** Before freezing, I ran parameter pilots on N = balanced_semiprime(36/40/44 bits, seed 999)
for E1, and on balanced_semiprime(40, 0..2) for E2/E3, only to set box sizes, row caps and budgets. Pilot facts that
shaped the design: at 40 bits the seedless base-m cubic reached full rational rank after ≈ n + C_alg + 6.3n rows
(seed 999), with output rank equal to the number of FB primes that occur on the rational side; orbit (d = 1) rows
reached rank n − 4 after 12n rows, the 4 missing primes never occurring; random-exponent rows reached full rank at
exactly n rows. Seed 999 is not in any registered cell; E2/E3 seeds 0–2 are (disclosed).

Code: `src/zroute/core.py` (operational, deterministic, sees only N), `src/zroute/harness.py` (harness; the random
module is used only by the registered *controls*). Config: `experiments/configs/zroute.json` (frozen with this
file). Outputs: `experiments/raw/zroute_e{1,2,3}.jsonl`. Factor base: FB = {−1} ∪ {primes ≤ fb_bound(N)},
fb_bound = least multiple of 64 that is ≥ 2 ln² N (detGRH's convention), so "PIN" (every prime < 2 ln² N in the
Q-span of the relation lattice) coincides with full rank up to the ≤ 64-wide rounding band.

Everything is toy-sized. These are correctness and mechanism checks, not complexity evidence.

## E1: seedless deterministic NFS Z route

- Instances: balanced semiprimes 36 bits (seeds 0–7), 40 bits (0–7); three-prime products 42 bits (0–3).
- Field: deterministic monic base-m cubic (m = ⌊N^{1/3}⌋). Rows: all coprime (a, b), b = 1..300, a = −1500..1500 in
  this order, with a − bm and F(a, b) both B-smooth (B = fb_bound); pairs with a prime of F(a, b) dividing disc f
  skipped (toy restriction for bad primes, disclosed).
- Stages: prefix M = smallest with M ≥ n + C_alg(prefix) + k·n, k ∈ (2, 4, 7); at each stage XK-Z (exact left
  kernel of the ideal columns; units removed by per-vector LLL with the Blanksby–Montgomery certificate), rational
  relations (verified exactly, sign fixed by e_{−1}), HNF, PIN, D = product of HNF pivots, deterministic Miller with
  bases 2..⌈2 ln² N⌉. Stop at the first stage with PIN.
- Predictions:
  - **P1a** every output relation verifies (ev = 1 after the sign fix): 0 failures; 0 NONUNIT kernel vectors;
    0 undecided LLL certifications.
  - **P1b** output rank = rank of the raw rational parts of the prefix rows, in ≥ 90% of stages.
  - **P1c** PIN is reached within the schedule in ≥ 16/20 instances.
  - **P1d** whenever PIN holds, Miller with D returns a proper factor with a base < 2 ln² N: 100% (a theorem under
    ERH; Prop Z/PIN).
  - **P1e** the number of log-independent pivots found equals r₁ + r₂ − 1 of the field in every instance where
    dim Λ₁ ≥ 3.

## E2: rank growth for d = 1 rows, deterministic orbit targets versus random-exponent seeds

- Instances: balanced semiprimes 40 bits, seeds 0–9.
- Orbit (deterministic): λ_j = 2^j mod N, j = bitlen(N) + 1, …; Gauss split λ = u/v; row j·e₂ − e(u) + e(v) when
  both halves are B-smooth; cap 30n rows or 400 000 targets.
- Random control (Lemma HR setting): x uniform in [0, N²)^n, λ = ev(x); same split test; cap 2n rows.
- Predictions:
  - **P2a** random: rows_to_full (mod 2^61 − 1 rank) = n exactly in ≥ 9/10 instances.
  - **P2b** orbit: rows_to_full > 3n (or never reached within the cap) in ≥ 8/10 instances.
  - **P2c** orbit: at the cap, n − rank = #{FB primes q ≠ 2 occurring in no half} in ≥ 9/10 instances (pinning is
    by occurrence only).
  - **P2d** whenever PIN holds, Miller splits N: 100%.

## E3: pointwise supply uniformity across N

- Instances: balanced semiprimes 40 bits, seeds 0–999. Budget J = 2000 targets per N for each design: orbit
  λ_j = 2^j (j from bitlen(N) + 1) and uniformly random λ (control). Count targets whose Gauss split has both
  halves B-smooth.
- Predictions:
  - **P3a** mean(orbit hits)/mean(random hits) ∈ [0.9, 1.1].
  - **P3b** min/median ≥ 0.5 for each design.
  - **P3c** variance/mean ≤ 3 for each design.
  - **P3d** (descriptive) number of N whose orbit repeats within J targets.

Analysis script: `src/zroute/analyze.py` (to be written after the runs; it only aggregates the JSONL fields named
above). Commands:
```
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e3 experiments/configs/zroute.json experiments/raw/zroute_e3.jsonl
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e2 experiments/configs/zroute.json experiments/raw/zroute_e2.jsonl
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e1 experiments/configs/zroute.json experiments/raw/zroute_e1.jsonl
```

## Addendum A1 (2026-10-22, before any E4 run; E1 still running, E2/E3 finished but not yet analysed)

Motivation: codex zroute_1 §5 proposed diagonal dominance to force rank deterministically (for integer smooth
representatives of q^H). I extend it to d = 1 split rows (Lemma DD1 in the results note). Pilot: one run on
balanced_semiprime(36, 999) (not a registered seed): all 207 primes got a row, the block was diagonally dominant.

## E4: deterministic rank forcing by diagonal dominance (d = 1)

- Instances: balanced semiprimes 40 bits (seeds 0–9), 48 bits (0–9), 56 bits (0–4).
- For each FB prime q (FB = primes ≤ fb_bound(N)), targets a = q^H·w mod N for w = 1, 2, … ≤ 20 000, B-smooth w
  only, in increasing order; H := K + bitlen(2N) + 2, K := max Ω(w) over the budget. Keep the first w whose Gauss
  split has both halves B-smooth; row H e_q + e(w) − e(u) + e(v). Then HNF, D, deterministic Miller.
- Predictions:
  - **P4a** whenever every prime has a row, the prime block is strictly diagonally dominant and the rank is n: 100%
    (a theorem: Ω(w) + Ω(u) + Ω(v) < H).
  - **P4b** whenever every prime has a row, Miller splits N: 100% (ERH theorem; toy check).
  - **P4c** every prime gets a row within the budget in ≥ 9/10 instances at 40 and at 48 bits and ≥ 4/5 at 56 bits.
  - **P4d** (descriptive with a window) at 40 bits the mean number of tries per prime lies in [0.5/δ, 2/δ], where δ
    is E3's random-target hit rate (pooled mean rand_hits/J).
- Command: `PYTHONPATH=. .venv/bin/python -m src.zroute.harness e4 experiments/configs/zroute.json experiments/raw/zroute_e4.jsonl`
