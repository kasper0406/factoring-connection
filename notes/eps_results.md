# Results: pre-registered correlation battery for ε = (−1/p) and QR_2 on D_n (agent `epsbattery`, 2026-10-18)

Pre-registration: `notes/eps_prereg.md`, committed in `779a650` **before** any battery instance was generated. Code: `src/eps/`.
Tables for every statistic: `experiments/analysis/eps_battery_tables.md`. Machine-readable summaries:
`experiments/analysis/eps_*.json`. Raw per-instance data: `experiments/raw/eps_{main,fresh}_n*.jsonl.gz`,
`experiments/raw/eps_pccurve_n*.jsonl.gz`.

## 0. Verdict

| size | target | validated (V1–V3) | confirmed menu hits | pre-registered conclusion |
|---|---|---|---|---|
| 32 | ε | yes (PC1, PC2, PC3) | 9, all CF prefix statistics | claim **not refuted**: CF prefix statistics predict ε |
| 32 | QR_2 | yes (PCQ2, PCQ3) | 5, all CF prefix statistics | claim **not refuted** |
| 48 | ε | yes (PC2, PC3; PC1 underpowered) | 0 | **refuted for this menu**: no effect `|ρ| ≥ 0.0216` |
| 48 | QR_2 | yes (PCQ2, PCQ3) | 0 | **refuted for this menu**: no effect `|ρ| ≥ 0.0306` |
| 64 | ε | yes (PC2, PC3; PC1 never fires) | 0 | **refuted for this menu**: no effect `|ρ| ≥ 0.0216` |
| 64 | QR_2 | yes (PCQ2, PCQ3) | 0 | **refuted for this menu**: no effect `|ρ| ≥ 0.0305` |

* All 7 one-sided positive controls (PC1–PC5, PCQ2, PCQ3) have **0 violations** in 3·10^5 battery instances,
  10^5 fresh instances and 1.6·10^5 supplementary instances (**experimentally supported**; each one is **proved** in the prereg §3).
* Every n = 32 hit is a CF-prefix statistic, and every one decays to `|z| < 3` by n = 48. The registered prediction for the
  mechanism (period closure) was **wrong** as stated: the hits survive on instances whose period does not close within 1024
  steps. A post-hoc analysis (§5) shows that they belong to the same family as PC2/PC3. That family is the restriction
  `((−1)^i Q_i | p) = 1` on the small values `Q_i` met in the cycle, which forbids, for example, square `Q_i` at odd `i` when
  ε = −1. Its strength per step decays like a negative power of `N` (`≈ N^{−1/4}` for squares).
* Interpretation: no Track A section, Track B sketch/rank/residue, BM length or small-form count (66 non-CF menu statistics)
  predicts ε or QR_2 beyond the floor at any size, at `|ρ| ≥ 0.022` (ε) / `0.031` (QR_2). The CF prefix exposes ε and QR_2 at
  n = 32 only through the classical relation mechanism, whose cost to reach a constant advantage grows exponentially or
  subexponentially (§6). Nothing here is a polynomial-time ε predictor.

## 1. What was run

* `D_n`: `balanced_semiprime(n, seed)` with `N ≡ 1 (mod 4)`, the first `10^5` qualifying seeds from `9.2·10^12` (main),
  `9.3·10^12` (fresh re-run, n = 32), `9.4·10^12 + n·10^9` (supplementary curve, `2·10^4` per size), all disjoint from earlier
  runs. QR_2 sample = `N ≡ 1 (mod 8)` part: 49 896 / 49 946 / 50 243 instances.
* 94 statistics (`MENU` in `src/eps/operational.py`; definitions in prereg §2), computed by `compute_stats(N)` (OPERATIONAL);
  labels by `src/eps/harness.py` (HARNESS).
* Test: stratified permutation z with strata `N mod 2^8`; hit iff `|z| ≥ 6.0` (S = 100 Bonferroni, family-wise error
  `≤ 2·10^{−7}` per size×target, `≤ 1.1·10^{−6}` over all 564 tests); replication on fresh seeds iff same sign and `|z| ≥ 4`.

## 2. Validity

**Positive controls (one-sided; `fires / z_obs / z_pred`; "powered" = `z_pred ≥ 8`).** In every row, `|z_obs − z_pred| ≤ 0.4`.

| n | PC1 T=1024 | PC2 T=1024 | PC3 T=1024 | PC2 T=64 | PC4 (MR) | PC5 (SS) | PCQ2 T=1024 | PCQ3 T=1024 |
|---|---|---|---|---|---|---|---|---|
| 32 | 1 426 / 38.0 / 38.0 | 49 216 / 311 / 311 | 49 760 / 315 / 315 | 12 949 / 122 / 122 | 11 / 3.3 | 9 / 3.0 | 24 584 / 212 / 212 | 25 712 / 222 / 222 |
| 48 | 4 / 2.0 / 2.0 (unpowered) | 12 951 / 122 / 122 | 42 746 / 273 / 273 | 898 / 30.0 / 30.0 | 0 | 0 | 4 699 / 72.5 / 72.6 | 22 538 / 204 / 204 |
| 64 | 0 (unpowered) | 941 / 30.8 / 30.9 | 394 / 19.9 / 19.9 | 63 / 8.0 / 8.0 | 0 | 0 | 339 / 18.4 / 18.4 | 213 / 14.5 / 14.6 |

* **PC1 (the control frozen in NEXT_STEPS)** is detected only at n = 32 (z = 38). At n = 48 it fires 4 times in `10^5`
  (`z_pred = 2.0`), and at n = 64 never. **PC1 cannot validate n = 48 or n = 64.** By the NEXT_STEPS rule alone those sizes would be
  invalid. They are validated by PC2/PC3 (and PCQ2/PCQ3), as registered in prereg §3 (V3) before the run.
* PC4/PC5 (MR square root of −1; SS `b = −1`; from GPT-6.1-sol r6, added before the freeze) fired 11 and 9 times at n = 32,
  all with ε = +1, and never at n ≥ 48. Same instance count: MR also produced a nontrivial square root of 1 (a factor) in
  9 instances at n = 32 (`mr_split`; STRUCTURAL: Fermat-pass events, rate ∝ `gcd(N−1,p−1)gcd(N−1,q−1)/φ(N)`).
* Planted control (`ρ = 0.05`): z = 17.2 / 16.0 / 15.5 (ε; expected 15.8) and 11.3 / 11.3 / 10.1 (QR_2; expected 11.2). All pass.
* Within-stratum shuffle: mean `z^2` = 0.69 / 0.98 / 0.93 (ε), 0.96 / 0.90 / 1.03 (QR_2).
* Negative controls: `nc_N_mod8`, `nc_jac2_N`, `nc_N_mod256` are DEGENERATE (exact zero) at every size, which shows that the
  floor removal is exact. The other NC statistics (bits 8–11, odd-modulus and size controls) have `|z| ≤ 2.1` everywhere.

## 3. Floor baseline

* ε: 2-fold cross-validated accuracy of the best `N mod 2^k` predictor, `k = 2..12`: 0.4985–0.5044 (n = 32),
  0.4963–0.5012 (48), 0.4991–0.5054 (64); SE 0.0016. χ² of ε against `N mod 2^8`: `z_WH` = 0.09, −0.80, 0.31. The floor is
  **experimentally null** at all sizes, as predicted by trackE Lemma F1.
* QR_2: accuracy equals the majority rate. At n = 32 that rate is 0.5181. This is a finite-set effect: 16-bit primes with the
  top two bits set number 377/363/368/383 in the classes 1/3/5/7 mod 8, which gives `P(QR_2 = 1 | N ≡ 1 mod 8) = 0.5194`
  (exact count). It is a constant offset, not information carried by `N`. χ² `z_WH` = 0.66, 1.33, −0.63.

## 4. Statistics: hits, nulls, minimum detectable effects

`ρ_MDE` (80% power at `|z| = 6`) = 0.0216 (ε, `neff = 99 936`), 0.0306 / 0.0306 / 0.0305 (QR_2). For a balanced binary
statistic this is a prediction accuracy of 51.1% (ε) and 51.5% (QR_2). For a rare binary statistic firing at rate `f`, it is a
shift `Δ = (ρ_MDE/2)√((1−f)/f)` in `P(y=1 | x=1)`. Full per-statistic z (stratified and raw): `eps_battery_tables.md`.

**Hits (all at n = 32, all CF group; z main → fresh-seed z; ρ at n = 32; z at n = 48 and n = 64):**

| target | statistic | z main | z fresh | ρ (n=32) | z n=48 | z n=64 | two-way cluster-robust z (post-hoc) |
|---|---|---|---|---|---|---|---|
| ε | cf64_maxbits | +7.4 | +8.1 | 0.024 | +0.3 | +0.9 | 6.3 |
| ε | cf64_n_Qeven | +6.2 | +5.1 | 0.020 | −2.0 | +1.0 | 5.4 |
| ε | cf64_n_Q1mod4_odd | +7.2 | +7.4 | 0.023 | +2.9 | +1.0 | 6.2 |
| ε | cf1024_n_a1 | +18.4 | +17.9 | 0.058 | +1.0 | −0.9 | 11.9 |
| ε | cf1024_n_a2 | −13.1 | −13.8 | −0.041 | −1.2 | +2.5 | −11.0 |
| ε | cf1024_sumbits | +7.7 | +6.6 | 0.024 | +0.0 | −0.9 | 6.7 |
| ε | cf1024_maxbits | +29.7 | +30.0 | 0.094 | +2.8 | +0.4 | 14.4 |
| ε | cf1024_n_Qeven | +18.8 | +16.6 | 0.059 | +1.2 | +0.9 | 12.4 |
| ε | cf1024_n_Q1mod4_odd | +24.7 | +23.5 | 0.078 | +2.7 | −0.5 | 16.0 |
| QR_2 | cf1024_n_a1 | +18.5 | +14.5 | 0.083 | −0.1 | +1.2 | 11.2 |
| QR_2 | cf1024_sumbits | −8.4 | −4.1 | −0.038 | +2.0 | +0.2 | −6.4 |
| QR_2 | cf1024_maxbits | +7.0 | +6.5 | 0.032 | +1.6 | +1.9 | 4.5 |
| QR_2 | cf1024_n_Qeven | −14.1 | −12.7 | −0.063 | −1.9 | +0.7 | −9.2 |
| QR_2 | cf1024_n_Q1mod4_odd | +13.4 | +11.8 | 0.060 | +1.2 | −0.4 | 8.7 |
| QR_2 | cf1024_n_sq_even | +6.4 | +2.1 (**not confirmed**) | 0.029 | +1.6 | +2.0 | 4.9 |

Fresh-seed run also had `|z| ≥ 6` for `cf64_n_a1`, `cf1024_n_2sq_odd` (ε) and `cf64_n_Q1mod4_odd` (QR_2); these were
below threshold in the main run (4.6, 5.3, 4.0), so they are not hits under the rule.

**Nulls.** Every one of the 66 non-CF menu and NC statistics is a null at every size and for both targets: max `|z|` 3.5
(n = 32), 2.7 (48), 2.1 (64); mean `z^2` over non-degenerate non-CF statistics 1.50 / 0.91 / 0.80 (ε) and
1.08 / 1.27 / 0.92 (QR_2). The value 1.50 at n = 32 is explained by prime reuse (§5.3; cluster-robust 1.10). At n = 48 and
64 the CF statistics are null too (max `|z|` 2.9). DEGENERATE everywhere: the three 2-adic NC statistics. Degenerate at
n ≥ 48 (generic full rank / full BM length): `Bx3m2_rank`, `bm_Fc`, `bm_coef64`, `bm_euler` (QR_2 n = 48 and both n = 64; ε n = 48 has only
`Bx3m2_rank`). These detect nothing, and the minimum detectable effect does not apply to them.

## 5. Leakage audit of the hits (prereg §6) and post-hoc analysis

### 5.1 Pre-registered steps (`experiments/analysis/eps_audit_n32.json`)
1. Static: all hit statistics come from `src/eps/cf.py::cf_features(N)`; the AST isolation test passes.
2. Recompute from `N` alone in a fresh process: 0 mismatches (2 000 instances per statistic). Post-hoc recomputation of all
   `10^5` with independent code (`src/eps/posthoc.py`): 0 mismatches.
3. Floor sensitivity: strata `N mod 2^12` and `(N mod 2^8, N mod 105)` change z by `≤ 0.6` for every hit.
4. Generator artefacts: hits correlate with `⌊16q/p⌋` and `nc_lead16` at `|z| ≤ 9` (CF statistics depend on the size of `√N`),
   but those quantities are null against both targets (`nc_lead16`: `|z| ≤ 1.8`). Seed parity: `|z| ≤ 2.3`.
5. Fresh seeds: 9/9 ε hits and 5/6 QR_2 hits confirmed (table above). Caveat (§5.3): at n = 32 the fresh instances reuse
   the same 1 491 primes per side.
6. Mechanism, registered diagnostic: on the 97 120 instances whose period does **not** close within 1024 steps, the hits
   remain (e.g. `cf1024_maxbits` 18.0, `cf1024_n_Q1mod4_odd` 27.5, `cf64_n_Q1mod4_odd` 7.1). **The registered period-closure
   explanation is refuted**, so by the registered rule these count as genuine menu hits at n = 32.

### 5.2 Post-hoc mechanism analysis (NOT pre-registered; `src/eps/posthoc.py`, `eps_posthoc_n32.json`)
For odd `i`, `A_{i−1}^2 ≡ −Q_i (mod N)`, so with `Q_i = k^2 m` (`gcd(k, N) = 1`), `(−m | p) = 1`, i.e. `(m | p) = ε`; for even `i`,
`(m | p) = 1`. In particular
`Q_i = k^2` at odd `i` is impossible when ε = −1 (this is PC2). More generally, which small `Q` values can occur at each index
parity depends on ε and on the small characters `(ℓ|p)`. Hypothesis H-SQ: the hits are carried by steps with non-squarefree
`Q_i`. Each count was recomputed over the squarefree steps and over the non-squarefree steps separately:

| ε statistic (T=1024) | all | squarefree steps | non-squarefree steps |
|---|---|---|---|
| n_Q1mod4_odd | +24.7 | −2.0 | +22.4 |
| n_Qeven | +18.8 | (empty: every even `Q_i` was divisible by 4) | +18.8 |
| n_a1 | +18.4 | +8.8 | +2.4 |
| n_a2 | −13.1 | −10.4 | −2.0 |
| maxbits | +29.7 | +14.3 | +27.5 |
| sumbits | +7.7 | −7.0 | +13.3 |

H-SQ explains `n_Q1mod4_odd` and `n_Qeven` (and `n_2sq_odd`). It does **not** explain the partial-quotient statistics
(`n_a1`, `n_a2`, `maxbits` keep `|z|` 9–14 on squarefree steps). Those depend on the frequency of small `Q_i` (large `a_i`).
A plausible but **unproved** mechanism is the parity-dependent admissibility of small squarefree values `m`. For ε = +1, a
small `m` with `(m|p) = +1` is admissible at both parities. For ε = −1, each admissible `m` is admissible at exactly one
parity. That changes how small norms are arranged along the cycle. Status: **conjectured**. For QR_2 the split is reversed
for `n_a1` and `n_Q1mod4_odd` (carried by squarefree steps); the mechanism is not identified (**not yet tested**).

### 5.3 Post-hoc: dependence between instances at n = 32
At n = 32 the primes have 16 bits; only 1 491 qualify per side, so each prime occurs in about 130 of the `10^5` instances.
The permutation variance assumes independent instances. A two-way cluster-robust variance (clusters by `p` and by `q`) gives
the right-hand column of the hit table. It shrinks z by a factor of 1.2–2. Still, 7/9 ε hits and 4/6 QR_2 hits keep `|z| ≥ 6`,
and the non-CF null mean `z^2` falls from 1.50 to 1.10 (max `|z|` 2.4). So the n = 32 conclusion stands, but the n = 32
z values in the main table are anti-conservative. At n ≥ 48 there are ≥ 2.5·10^5 primes per side and reuse is negligible.
This is a limitation of `D_32` itself, not of the test; it was not anticipated in the prereg.

## 6. How much of ε the CF prefix exposes cheaply (supplementary curve, `eps_pccurve_summary.json`)

Firing frequency, `2·10^4` instances per size; `P(ε = +1) ≈ 0.50`. 0 violations in every cell.

| n | PC1 T=1024 | PC1 T=16384 | PC2 T=1024 | PC2 T=16384 | PC3 T=1024 | PC3 T=16384 |
|---|---|---|---|---|---|---|
| 24 | 0.158 | 0.292 | 0.450 | 0.450 | 0.452 | 0.452 |
| 32 | 0.016 | 0.180 | 0.491 | 0.496 | 0.498 | 0.498 |
| 40 | 0.0010 | 0.013 | 0.355 | 0.501 | 0.500 | 0.502 |
| 48 | 0.0002 | 0.00075 | 0.133 | 0.501 | 0.432 | 0.505 |
| 56 | 0 | 0.0001 | 0.037 | 0.342 | 0.095 | 0.494 |
| 64 | 0 | 0 | 0.0096 | 0.133 | 0.0034 | 0.460 |
| 80 | 0 | 0 | 0.0008 | 0.0091 | 0.0001 | 0.0079 |
| 96 | 0 | 0 | 0 | 0.00045 | 0 | 0 |

* PC1 (negative Pell via an odd period ≤ T) falls ~5–16× per 8 bits at T = 1024, roughly consistent with `T/√N` (the period is `≈ √N` in typical
  cases). It exposes almost none of ε beyond n ≈ 40 (**experimentally supported**).
* PC2 (an odd-index square `Q_i`) falls ~4× per 8 bits at fixed T: rate `≈ T·N^{−1/4}`, the SQUFOF scaling.
* PC3 (CFRAC linear algebra over primes ≤ 256) is the strongest: at T = 16 384 it certifies ε = +1 for ~92% of ε = +1
  instances at n = 64, and 0 at n = 96. It needs about as many smooth relations as CFRAC needs to factor `N`, so its cost is
  `L_N[1/2]`-type. It is not cheaper than factoring by CFRAC.
* All three are one-sided (they certify ε = +1). Their decay sets the scale of the n = 32 menu hits, which are weak,
  two-sided versions of the same relations.

## 7. Nonunits (rule 3)

Values whose Jacobi symbol was taken, BM discrepancies and pivots are gcd-checked. Proper divisors found (ACCIDENTAL-null
expectation `Σ checks·(1/p + 1/q)`): n = 32: 2 925 vs 898; n = 48: 18 vs 3.5; n = 64: 0 vs 0.014. By source at n = 32:
`bm_euler` 849, `A16c5` 526, `bm_coef64` 476, `bm_Fc` 415, `Bphi12` 175, `Bphi5` 153, `A8c3` 110, `Bfib` 83, others ≤ 21. Single
values (`A1c*`, `A2c3`) match the null. Larger sections, cyclotomic sketches and BM runs exceed it by 3–20×. That fits the
Frobenius-coincidence rates of Track B (B1: roots appear at rates `∝ (small order)/p`). Classification: STRUCTURAL-candidate
(rate Θ(1/p) with larger constants; mechanism not analysed here). All rates vanish with `p`. One Φ_5 sketch coefficient is
exactly 0 mod N at a position fixed by `N mod 5`: 99 979 / 100 105 / 100 192 non-proper events, as noted in the prereg. MR found
factors in 9 instances at n = 32 (§2). The CF code had 0 gcd rejections, and the factor base was coprime to N in every instance.

## 8. Cost

Main battery: 547 / 555 / 1 064 CPU-s (n = 32/48/64), mean 5.5 / 5.6 / 10.6 ms per instance, ~2 min wall on 30 workers.
Fresh re-run 445 CPU-s; supplementary curve ≈ 53 CPU-min. Total ≈ 1.5 CPU-h including analysis. No statistic was dropped.

## 9. Deviations from the pre-registration

* Added after the freeze: `src/eps/tables.py` (formatting), `src/eps/posthoc.py` (§5.2, §5.3; labelled post-hoc). The frozen
  modules were not changed.
* The registered mechanism prediction for CF hits (§2a of the prereg) was wrong; reported as such (§5.1 step 6).
* The fresh re-run covers only n = 32, the only size with hits. It is not prime-disjoint from the main run (§5.3).

## 10. Interpretation

* **What is excluded (experimentally supported).** For this menu of 94 statistics, on `D_48` and `D_64`, no statistic
  correlates with ε at `|ρ| ≥ 0.0216`, or with QR_2 on `N ≡ 1 (mod 8)` at `|ρ| ≥ 0.031`, beyond the `N mod 2^8` floor
  (80% power; validated by PC2/PC3/PCQ2/PCQ3 and the planted control). This excludes only effects of that size. It says
  nothing about smaller correlations, other statistics, or larger `n`. The Track A section values, Track B sketch
  coefficients, ranks and residues, BM lengths and `a ≤ 50` form counts are null even at n = 32 (`|ρ| < 0.02`).
* **What is found.** At n = 32, CF-prefix statistics of `√N` predict ε (`|ρ|` up to 0.094) and QR_2 (up to 0.083). This holds
  through fresh seeds, finer floors and cluster-robust variance. The effect is the classical relation
  `A_{i−1}^2 ≡ (−1)^i Q_i (mod N)` acting on small `Q_i`. It is gone by n = 48, with per-step strength falling like a negative power
  of `N`. It is not a polynomial-time signal. Its one-sided extremes (PC2, PC3) need `T ≈ N^{1/4}` (SQUFOF) or `L_N[1/2]`
  (CFRAC) steps for a constant advantage. The same relations then factor `N` anyway.
* **Negative Pell (PC1)** is detectable only at n ≤ 32 with `T ≤ 1024`. At n ≥ 40 it exposes a fraction `< 10^{−3}` of the
  instances.
* **Context.** Deciding ε for balanced semiprimes is Seurin's conjectured-hard "2-Φ/4-Hiding" problem (PKC 2014, eprint
  2013/256; see `notes/codex_r6_eps_reduction.md`). His prime windows and mixture weights differ from `D_n`. The nulls here are
  consistent with that conjecture, and nothing here bears on the reduction question (NEXT_STEPS item 3).
* No novelty is claimed. The relation mechanism is the textbook basis of CFRAC and SQUFOF.

## 11. Reproduction (project root)

```
.venv/bin/python -m pytest -q tests/test_eps.py                                       # 20 tests
.venv/bin/python -m src.eps.harness --n 32 48 64 --range main --workers 30            # ~2 min wall
.venv/bin/python -m src.eps.harness --n 32 --range fresh --workers 16
.venv/bin/python -m src.eps.harness --n 24 32 40 48 56 64 80 96 --range supp --curve --count 20000 --workers 30
.venv/bin/python -m src.eps.analysis experiments/raw/eps_main_n32.jsonl.gz experiments/raw/eps_main_n48.jsonl.gz \
    experiments/raw/eps_main_n64.jsonl.gz --out experiments/analysis/eps_battery_summary.json
.venv/bin/python -m src.eps.analysis experiments/raw/eps_fresh_n32.jsonl.gz --out experiments/analysis/eps_fresh_n32_summary.json
.venv/bin/python -m src.eps.analysis experiments/raw/eps_main_n32.jsonl.gz --out experiments/analysis/eps_audit_n32.json \
    --audit cf64_maxbits:eps cf64_n_Qeven:eps cf64_n_Q1mod4_odd:eps cf1024_n_a1:eps cf1024_n_a2:eps cf1024_sumbits:eps \
    cf1024_maxbits:eps cf1024_n_Qeven:eps cf1024_n_Q1mod4_odd:eps cf1024_n_a1:qr2 cf1024_sumbits:qr2 cf1024_maxbits:qr2 \
    cf1024_n_Qeven:qr2 cf1024_n_sq_even:qr2 cf1024_n_Q1mod4_odd:qr2
.venv/bin/python -m src.eps.analysis --curve experiments/raw/eps_pccurve_n{24,32,40,48,56,64,80,96}.jsonl.gz \
    --out experiments/analysis/eps_pccurve_summary.json
.venv/bin/python -m src.eps.posthoc experiments/raw/eps_main_n32.jsonl.gz                      # post-hoc §5.2
.venv/bin/python -c "from src.eps.posthoc import cluster_main; cluster_main('experiments/raw/eps_main_n32.jsonl.gz', None, 'experiments/analysis/eps_cluster_n32.json')"   # post-hoc §5.3
.venv/bin/python -m src.eps.tables
```
