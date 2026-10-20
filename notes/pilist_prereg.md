# pilist pre-registration (phase 18)

Agent `pilist`, 2026-10-20. Written and committed **before** any registered run of E1–E3. Code: `src/pilist/`
(`pencil.py`, `decoder.py`, `tails.py`, `run_experiment.py`), config `experiments/configs/pilist.json`, tests
`tests/test_pilist.py` (15, all on toy inputs). Raw output: `experiments/raw/pilist_E{1,2,3}.jsonl`.

**Disclosed pilots (before this file).**
- Functional checks of the line finder on (M, B) = (3, 3) and (4, 3). These found the M = 3 identity
  X(X+3)² − 2(X+1)³ + (X+2)(X−1)² = 0 (now a unit test). No registered cell was run.
- One toy row of the c60 rational polynomial (b = 37, X = 2^12, z = 60) to check the HG decoder. It returned
  certified thresholds 21.5 → 19.95 nats for m = 4..12 against the Boneh threshold 20.15. It found 0 outputs,
  with 0.024 expected under the model. A 200-row empirical count gave 0.015 per row. This row/cell is not in
  the registered grid.
- Timing of `simulate` at n = 500 and n = 1000 (k = 8, j = 2, one seed each).
- The model predictions below are computed from the model only, with no experimental data.

Memory: every job streams or holds ≤ 10⁴ coefficient vectors; expected RSS < 500 MB. `memory_pressure` is checked
before each run, and runs are sequential.

---

## E1 — pencils and three-term identities of split polynomials (Task 4, PTE O2)

**Setup.** For each (M, B) ∈ {(3,8), (4,6), (5,5), (6,4), (7,4), (8,3), (9,3), (10,3)}:
- enumerate all monic split P = ∏(X − a), a ∈ [−B, B] (multisets);
- find every line in coefficient space with ≥ 3 members:
  - kproj = 0: exact pencils F₀ + tF₁;
  - kproj = 1: lines after dropping the constant coefficient, i.e. O2 systems F₀ + νF₁ + const.
- For every line, record:
  - the fixed part (multiset intersection) and moving degree e;
  - Claim A (pairwise disjoint moving roots) and the number of non-squarefree members;
  - every member triple whose primitive relation c ∝ (t₁ − t₂, t₂ − t₀, t₀ − t₁) has max|c_i| ≤ 2. For each such
    triple: G = Σc_iP_i (0, or a constant), its Lemma TL class A/B/C, and whether it is "in regime". In regime
    means |G| ≥ 1 and M ≥ 1 + 2⌊log₂|G|⌋ + 3 (Lemma R's j = 0, K = 1 threshold), or G = 0 and e ≥ 6.

**Predictions.**
- **P1 (PROVED, Claim A of Lemma PEN).** Every kproj = 0 line has pairwise disjoint moving roots, and ≤ 2e − 2
  members with a repeated moving root. Prediction: 0 violations in all cells.
- **P2 (sharpness of e ≥ 3).** Every cell has kproj = 0 lines with e ≤ 2 and ≥ B + 1 members.
- **P3 (the naive three-term rigidity is false).** For M = 3 and M = 4 there is a kproj = 0 line with e = M
  (no fixed part), ≥ 3 members and a |c| ≤ 2 triple.
- **P4 (exploratory; the in-regime question).** Zero class-C, in-regime, small-c triples in every cell with M ≥ 6.
  Any such triple would be a counterexample to the Lemma R-shaped three-term rigidity, and is reported as a
  positive finding.
- **P5 (exploratory expectation, not theory-backed).** For M ≥ 6, the maximum number of members on a
  kproj = 0 line with e ≥ 3 is ≤ 6.

**Kill signal.** A violation of P1 means a bug or a false lemma. Stop and debug before reporting anything.

## E2 — Johnson-radius decoders on the c60 rational side (Tasks 2–3, Prop JR)

**Setup.**
- Polynomial: `experiments/raw/cadobench_testbed_c60`, G(a, b) = Y₁a + Y₀b with Y₁ = 3961719737 = 53·233·541·593
  and Y₀ = −166188820558869.
- Rows: the first R values of b ≥ 1 with gcd(b, Y₁) = 1. (X, R) ∈ {(2^10, 200), (2^13, 200), (2^16, 64)}; the row
  is |a| ≤ X, V = 2X + 1.
- Smoothness: B = 78682 (lim0). Truth is computed by sieve, then confirmed by trial division, with gcd(a, b) = 1.
- Decoder modulus: P = ∏_{q ≤ z, q ∤ Y₁} q, for z ∈ {13, 20, 30, 60, 100, 160}.
- Threshold: S at Boneh's bound, log S = √(log V · log P). A cell is skipped if log S > log P.
- Outputs per call: the truth set {a : gcd(G, P) ≥ S}.
- Decoder check: on the first 16 rows of each cell, run HG with m ∈ (4, 6, 8, 12, 16) and t = 2m, until the
  certified threshold is within 0.5 nats of log S. Every a with gcd > S_cert must be returned.
- Γ̂ = (smooth outputs per call) / (box smooth rate per coprime point).

**Model predictions** (random-divisibility model; expected outputs per call, from `count_model_fraction`; the
JR Rankin bound uses λ = log(u log u)/log B with u = 3.52, λ = 0.132):

| X \ z | 13 | 20 | 30 | 60 | 100 | 160 |
|---|---|---|---|---|---|---|
| 2^10 | 0.273 (Γ ≤ 8.9) | 0.267 (11.3) | 0.161 (9.8) | 0.0262 (4.1) | 0.0026 (1.07) | 1.3e−4 (0.17) |
| 2^13 | 0.546 (7.9) | 0.324 (15.5) | 0.169 (14.3) | 0.0229 (5.3) | 0.0016 (1.07) | 5e−5 (0.12) |
| 2^16 | skipped | 0.270 (19.2) | 0.212 (20.8) | 0.0221 (7.2) | 0.0011 (1.15) | 2e−5 (0.09) |

**Predictions.**
- **P6 (PROVED; HG lemma).** 0 decoder misses: every a with gcd(G(a, b), P) > S_cert is returned.
- **P7 (count model).** In every cell with pred_total ≥ 10: observed/predicted total outputs ∈ [0.6, 1.6].
- **P8 (Prop JR at the threshold).** Mean outputs per call < 1 in every non-skipped cell.
- **P9 (yield).** Γ̂ ≤ the JR bound of the table in every cell where ≥ 5 smooth outputs are expected
  (pred_total × box rate × JR bound ≥ 5); and Γ̂ ≤ 25 in every cell.

## E3 — adaptive closure of Θ(u) tails (Task 3, Prop TAIL)

**Setup.**
- n uniform band columns (p = 1/n), n ∈ {1000, 2000, 4000, 8000}.
- (k, j) ∈ {(3,0), (5,0), (5,1), (8,0), (8,2), (12,2), (12,4), (16,6)}, 20 seeds each.
- Strategies:
  - `nonadaptive`: k random entries, j = 0;
  - `greedy` (legal): the j forced entries are chosen before the random part, among the columns seen exactly once;
  - `oracle` (ILLEGAL in the model): the forced entries are chosen after the random part, to cancel its residual.
- R* is the number of relations at the first F₂ dependency on the band columns.
- R_lb is the first-moment threshold with Pr[dependency by R_lb] ≤ 0.05 (`r_threshold(k, j, 1/n, 0.05)`), frozen
  values:

| n | 3,0 | 5,0 | 5,1 | 8,0 | 8,2 | 12,2 | 12,4 | 16,6 |
|---|---|---|---|---|---|---|---|---|
| 1000 | 711 | 521 | 177 | 337 | 99 | 92 | 50 | 32 |
| 2000 | 1425 | 1044 | 355 | 675 | 198 | 184 | 100 | 65 |
| 4000 | 2855 | 2089 | 711 | 1351 | 397 | 368 | 200 | 130 |
| 8000 | 5714 | 4180 | 1423 | 2704 | 794 | 736 | 401 | 261 |

**Predictions.**
- **P10 (Prop TAIL, PROVED in the model).** For legal strategies with k ≥ 2j + 2, at most 3 of 20 runs per cell
  have R* ≤ R_lb (expected ≤ 1).
- **P11.** Legal strategies: median R*/n ≥ 0.5 in every cell.
- **P12 (exploratory).** Report oracle against greedy medians. No directional prediction.
