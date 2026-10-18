# dchworst: pre-registration (frozen before the pre-registered runs)

Agent `dchworst`, 2026-10-18. Committed together with `src/dchworst/*` and
`tests/test_dchworst.py` before any run of `src/dchworst/experiments.py e1..e5`.

## 0. Disclosure of runs before this file

- Unit tests (n = 20, 22; tiny bands). A timing run of `coverage.sweep` at n = 32 with the first
  5 primes of the band (2^13, 2^14], signed weights: only the wall time was printed (0.23 s),
  no h values were looked at.
- Generator smoke tests: 3 pairs per structured family at n = 72 (membership in S_n only), and
  100 random pairs at n = 40 (count only). No coverage computed.
- **Design-only predictions** (no pairs involved, only the per-prime tables ν_ℓ(H) and the
  population sizes |S_n|) were computed and are frozen below.

## 1. Question

Theorem D′ covers all but N^{−ε} of balanced pairs. A worst-case N^{1/6+o(1)} needs every pair
(p, q) to be covered: some weight (A, B) of small height and some band prime ℓ with ℓ | Bp − Aq.
(1) What do uncovered pairs look like, can they be constructed, and is DCH false?
(2) Is there a design with *provable* full coverage, and what is the best provable worst-case
exponent? (3) Experiments: exhaustive worst case at small size, adversarial search.

## 2. Theory program (to be proved or refuted in notes/dchworst_results.md)

- T1 (criterion). Coverage ⇔ some Λ_ℓ = {(A,B) : ℓ | Bp − Aq} = Z(p,q) + ℓZ² has a primitive
  vector of height ≤ H in the weight set.
- T2 (CRT degrees of freedom). How many band constraints can be forced by congruences.
- T3 (necessity of λ → ∞). Union bound + Brun–Titchmarsh: if λ = |V| Σ 1/ℓ is small, a positive
  proportion of S_n is uncovered; hence DCH fails for H ≤ c√log N.
- T4 (Theorem W, cost criterion). A coverage guarantee H = √Y·J^{−c} for every pair gives the
  worst-case exponent min(1/5, (1−c)/(7/2−c)); below 1/5 iff c > 3/8. Minkowski/McKee–Pinch
  gives c = 0 (signed weights).
- T5 (barriers). (a) CRT-local families; (b) "main term − uniform error" proofs with error
  ≥ J^{1/2} give c ≤ 1/4, hence no gain over 1/5.
- T6 (designs). Composite, two-band-prime, smooth, growing-weight, Lehman-hybrid, Umans–Wang
  prime-cover options: cost and provability.

## 3. Experiments (frozen)

Notation: μ_ℓ(t) = least height max(|A|, B) of a weight with A B^{−1} ≡ t = p q^{−1} (mod ℓ).
h_L(p,q) = min_{ℓ∈L} μ_ℓ. Weight sets: *signed* {1 ≤ |A|, B ≤ H, primitive, ≠ (1,1)} and
*positive* (Theorem D′'s V). Poisson model: Pr[h_L > H] = ∏_{ℓ∈L}(1 − ν_ℓ(H)/(ℓ−1)), with the
exact class counts ν_ℓ(H). h_pred = least H with |pop|·Pr[h > H] < 1.

### e1: exhaustive worst case (all pairs of S_n)
- n = 32 (6 185 272 pairs) and n = 36 (78 295 606 pairs); bands L = primes ≡ 3 (mod 4) in
  (Y, 2Y], Y = 2^9, 2^11, 2^13 (37, 128, 440 primes), nested prefixes in increasing order;
  Hcap = ⌊√(2Y)⌋ + 1; signed and positive weights.
- Frozen design-only predictions h_pred(J) (signed): n=32, Y=2^13: J=16→65, 32→52, 64→40,
  128→29, 256→23, 440→19. n=36, Y=2^13: 68, 55, 43, 32, 24, 19. Positive, Y=2^13, n=32:
  93, 73, 56, 42, 31, 26.
- **Predictions.**
  - P1 (proved; sanity): signed h*_1 ≤ ⌊√ℓ₁⌋ (Minkowski).
  - P2: for every (n, Y, sign) and every recorded J with h_pred(J) ≤ Hcap,
    h_pred(J) − 4 ≤ h*_J ≤ h_pred(J) + 3.
  - P3: tail fit: for every recorded (J, H) with expected count E ≥ 20, the observed count is
    within [0.5E, 2E]; ≥ 95% of such cells within [0.7E, 1.4E].
  - P4: fitted local exponent ĉ = −d log h*_J/d log J over J ∈ [16, |L|] within ±0.08 of the
    same fit on h_pred(J) (Y = 2^11, 2^13).
- **Decision rule.** If P2–P4 hold, "the exhaustive worst case is the random-model extreme; no
  structured obstruction at these sizes". If P2 fails with h* > h_pred + 3 anywhere, report the
  worst pairs as a structured family and analyse them.

### e2: sampled tail at the D′ scaling
- n = 60, 72, 84 with 200 000 / 100 000 / 50 000 uniform pairs of S_n; Y = ⌊2^{n/6}⌋, the four
  blocks (2^kY, 2^{k+1}Y], ℓ ≡ 3 (mod 4); per block and all four; Hcap = 64; both weight sets.
- **Prediction.** Same P2/P3 rules (sample size in place of |S_n|).

### e3: structured families (n = 72, block 0 of Y = 4096, 231 primes)
- Families (≈ 20 000 pairs each): twin_gap (q = nextprime(p)), near_2p, near_3_2, near_5_4,
  cong_p_M30 (q ≡ p mod ∏_{r≤29} r), cong_mp_M30 (q ≡ −p), crt_worst2 (t_ℓ forced to the worst
  class at the 2 smallest band primes); baseline 100 000 random pairs.
- Reference heights (design-only): positive H_ref = 8 (Poisson uncovered 0.19), signed
  H_ref = 6 (0.169).
- **Prediction.** For every family, uncovered fraction at H_ref within [0.5, 2]× the random
  baseline, unless the family is a near-rational one (near_*, twin_gap), where a single small
  |Bp − Aq| removes one weight's chances (expected effect ≤ 1/|V| relative). crt_worst2 within
  [0.8, 1.25]× baseline (forcing 2 of 231 primes changes λ by < 2%).
- **Decision.** A family with ratio > 2 and a 3σ binomial excess is reported as a structured
  adversary.

### e4: exhaustive D′-design coverage at small n
- n = 24, 28, 32, H = 2..10; positive weights; 4 blocks; Y = ⌈2^{n/6}/⌊H^{1/3}⌋⌉. Records
  uncovered (no block covers) and uncovered_clean (no clean block covers; clean = no band prime of
  the block divides q − 1); the Minkowski-regime flag H² ≥ Y.
- Frozen Poisson expected uncovered counts: n=32: H=5 → 18 646, H=6 → 4 121, H=7 → 29.4,
  H=8 → 0.13; n=28: H=6 → 271, H=7 → 1.9; n=24: H=5 → 12.1, H=6 → 2.0, H=7 → 0.01.
- **Prediction.** For cells with expected ≥ 20, observed within [0.5, 2]× expected. The minimal
  full-coverage H* is within ±1 of the least H with expected < 0.5.

### e5: random search for pairs uncovered in all four blocks (n = 64)
- 10⁶ uniform pairs; D′ design, positive weights, H = 7 (Y = 1625) and H = 9 (Y = 813).
- Frozen Poisson expectations: H = 7: all-4 uncovered ≈ 4 702 (per block ≈ 2.4·10⁵, 2.3·10⁵,
  2.8·10⁵, 3.0·10⁵); H = 9: all-4 ≈ 75 (per block ≈ 5.7·10⁴, 1.0·10⁵, 1.0·10⁵, 1.3·10⁵).
- **Prediction.** Observed all-4 count within [0.5, 2]× expected (block independence), per block
  within [0.8, 1.25]×.

## 4. Review

Any claimed worst-case theorem goes to an independent read-only Opus critic. GPT-6.1-sol
(codex/dchworst_1) makes a blind attempt in parallel; its output is read after this file is
frozen.

## Addendum A (2026-10-18, restart; frozen before any e6 run)

The e1–e5 runs above were started on the new machine exactly as frozen (code unchanged since
commit bde6c64). This addendum adds one experiment, e6, for a new theorem (Theorem R,
notes/dchworst_results.md §3): **DCL-wc** (`src/dchworst/wcalg.py`) outputs a proper divisor
for every balanced semiprime N = pq whose pair is covered by at least one block. It needs no
clean band, no order hypothesis mod q, and no FULL exclusion.

**Runs before this addendum (disclosed).**
- Unit tests `tests/test_dchworst_wc.py`, including an exhaustive check at n = 20 (both
  designs, both modes): 0 failures on covered pairs.
- Smoke runs at n = 32–48: 40 random pairs per n; 300 covered random pairs per n at
  n = 32, 36, 40 against D′; and pairs with q | 2^k − 1 at n = 36, 40, 44. There were 0
  failures on covered pairs, and D′ also failed on none of them.
- A constructed equal-order family at n ≈ 28–40 was used to exercise the equal-order split.

**e6a (exhaustive).**
- Population: every pair of S_n for n = 26 and 28.
- Design: H = 3, Y = max(3, 2^{⌊n/6⌋}/2), and 4 blocks of primes ≡ 3 (mod 4) in
  (2^kY, 2^{k+1}Y].
- Weight designs: "pos" (D′'s V) and "sgn" (signed weights plus (1,1)).
- Modes: "full" (Lemma O with T = 2^{⌈n/3⌉}, B0 = 30) and "force2" (base 2, no Lemma O).
- D′ is run on (pos, full) for comparison.

**e6b (adversarial).**
- Sizes: n = 40, 48, 56.
- Families:
  - small2_q: q | 2^k − 1 for some k ≤ 240;
  - small2_both: p, q both such primes;
  - eqord: p − 1 = j·ℓ₀·o and q − 1 = k·ℓ₀·o with ord_p(2) = ord_q(2) = ℓ₀o, where ℓ₀ is a
    band prime and o ≤ 60.
- Same designs and modes as e6a.

**Predictions.**
- F1 (proved, Theorem R): in every cell, covered_fail = 0.
- F2: improper outputs = 0 in every cell.
- F3 (proved, Lemma X): max_frontier_over_Ngiants ≤ 1 in every cell.
- F4 (expectation, not proved): in e6b eqord under force2, equal_order_runs > 0 at every n.

**Decision rule.** Any covered_fail > 0 or improper > 0 means Theorem R or its
implementation is wrong. I will then report the counterexample and withdraw or repair the
theorem before reporting it.
