# crtdecode pre-registration (phase 23)

Agent `crtdecode`, 2026-10-22. Written and committed before any experiment code exists or runs.

## Question
Can an ADAPTIVE decoder of the CRT code behind Prop SI (positions j in [0, h], received word -X mod q^e,
agreement ln gcd(X + j, P)) find a y-smooth X + j more cheaply than the oblivious divisor trick? Toy experiments
cannot see L-asymptotics. They test (i) whether the most natural adaptive aggregation rule (pointer heavy-hitters:
positions hit by several smooth divisors get tested first) beats the divisor trick at matched cost, and (ii) the
information/computation split behind Lemma IC, made concrete as a decision tree of decoder calls.

## Common setup
- Configs. C1: X uniform in [2^43, 2^44), y = 2^9, h = 2^14. C2: X uniform in [2^51, 2^52), y = 2^10, h = 2^16.
  (h < X^{1/3} in both, so the toy is not in the h >= sqrt(X) regime.)
- Windows: 300 per config for E1, the first 100 of the same windows for E2. Seeds: `random.Random("crtdecode-E1-C1")`
  and `"crtdecode-E1-C2"` for X; repertoire orders from `random.Random("crtdecode-rep-C1")` / `"crtdecode-rep-C2"`,
  fixed before and independent of X.
- Ground truth: exact sieve of [X, X + h] (all integers, all prime powers q^e <= X + h, q <= y). Exact arithmetic only.
- Cost units (E1): 1 per pointer (computing the multiples of one d in the window), 1 per smoothness test of one
  position (exact trial division over primes <= y). A position is never tested twice by one decoder.
- Budget cap: B = 3000 (C1), 6000 (C2) cost units per window. Metric: cost to the first smooth output; windows not
  solved within B are scored at B ("capped mean"). Windows with no smooth integer are reported and excluded from the
  capped mean (no decoder can succeed there). Success rate within B is also reported.

## E1: adaptive pointer heavy-hitters vs the divisor trick
Decoders (all see only X, h, y and public config; no ground truth):
- `DT`  (oblivious divisor trick): repertoire R1 = all y-smooth d in [h/2, h], fixed random order. For each d: pay 1,
  then test each not-yet-tested multiple of d in the window (pay 1 each).
- `DT8` (control): as DT, with repertoire R2 = all y-smooth d in [h/2, 8h] (fixed random order).
- `PH-1.25`, `PH-1.5` (adaptive heavy-hitter): walk R2 in the same order, paying 1 per pointer and recording for each
  hit position j the known divisor L[j] = lcm of the d that hit it. A position is tested (pay 1) as soon as
  ln L[j] >= theta * ln h (theta = 1.25 or 1.5). If R2 is exhausted before the budget, the untested hit positions are
  tested in decreasing order of L[j].
- References (not adaptive decoders): `RT` = test positions in a fixed random order; `SV` = sieve cost to first, i.e.
  the number of (prime power, position) incidences of the window prefix up to and including the first smooth
  position.
Predictions:
- **P1 (main).** For both PH variants and both configs, capped mean(PH) >= 0.8 * capped mean(DT). (Adaptive
  heavy-hitter aggregation does not beat the divisor trick by more than 20%.)
- **P2 (enhancement).** Among positions tested by PH (both thresholds pooled) and by DT8, the fraction that is
  smooth is at least 2x higher for positions hit by >= 2 pointers than for positions hit by exactly 1 pointer, in both
  configs (counted only if each group has >= 30 tested positions).
- **P3 (ordering of references).** capped mean(DT) < capped mean(RT) < mean SV cost in both configs.

## E2: the information decoder (decision tree of idealised Johnson calls)
Idealised call (Theorem J's class, unweighted): modulus M (y-smooth), region R; certified set
{z in R : (ln gcd(z, M))^2 >= ln|R| * ln M}; the answer is the list of (z, gcd(z, M)).
Tree (adaptive; X enters only through answers):
1. Group the prime powers q^{a_q}, a_q = floor(log_q(h + 1)), into products M <= h + 1 (greedy, descending); for each
   group query M with region [X, X + M - 1] (a sub-window), which certifies exactly the unique multiple of M and so
   reveals X mod M.
2. For each q and e = a_q + 1, a_q + 2, ... <= e_q (e_q = floor(log_q(X + h))): query M = q^e with region = the
   window; stop at the first e with no position of valuation >= e in the window.
3. Reconstruct all valuations in the window from the answers (the "leaf computation": a sieve on known residues,
   counted in incidences) and output a smooth position z*, certified by one last call with M = z*.
Predictions:
- **P4.** The reconstructed smooth set equals the sieve's smooth set in 100% of the 100 windows per config.
- **P5.** Depth (number of calls) <= 3 * pi(y) + 1 in every window (C1: <= 292; C2: <= 517).
- Reported, not predicted: leaf computation cost / h.

## Analysis and reporting
- Raw JSONL: `experiments/raw/crtdecode_E1.jsonl`, `experiments/raw/crtdecode_E2.jsonl`.
- Commands: `.venv/bin/python -m src.crtdecode.run_exp E1|E2`.
- Any change after this commit is a disclosed deviation; post-hoc analyses are labelled POST-HOC.
- Interpretation fixed in advance: P1 holding is consistent with (not evidence for) the asymptotic claim that this
  adaptive class gains at most the divisor trick's L[1/3] factor; P1 failing at toy size would be reported as a toy
  effect to be explained before any asymptotic reading. E2 is a demonstration of Prop DT (depth ~ pi(y), leaf = sieve),
  not evidence about computation.
