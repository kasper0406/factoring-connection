# fruit: pre-registration (phase 30)

Agent `fruit`, 2026-10-28. Committed before any E1/E2/E3 run counted below. Code is frozen at this commit:
`src/fruit/{linalg,qs_hybrid,nfs_hybrid,harness_e1,harness_e2,theory}.py`; configs
`experiments/configs/fruit_e{1,2}.json`. The analysis script (`src/fruit/analyze.py`) is written after this commit and
implements exactly the tests below.

**Disclosed pilots (timing/debugging only, not analysed):** `scratchpad/fruit/pilot_e{1,2}.py`, one 40-bit and one
44-bit NFS instance, one 64-bit QS instance, two rounds each. They showed (i) the code runs (≈0.4 s per QS instance,
≈10 s per 44-bit NFS instance), and (ii) that at base size M = C + 5 a single descended row is often *not* in the
rowspace of the base matrix (r = 0). (ii) motivated adding the base-size levels 2C and 4C to E1 before freezing.

## What is tested

Lemma HY (proved in the notes): with base rows A (M rows) and h hidden-root rows D, a uniform element of the
left kernel K of [A; D] uses a hidden-root row with probability exactly 1 − 2^{−r}, r = h − (rank[A;D] − rank A).
Theorem HYB (= Theorem DD applied to the hybrid): a kernel element chosen from public data among those that use a
hidden-root row is fruitful with probability exactly 1 − 2^{1−ω(N)}, whatever the base rows are.

Independence: each *round* draws fresh hidden roots (fresh descents) on the same base matrix; one forced and one
uniform kernel sample are drawn per (round, h). By Theorem DD the forced outcomes of different rounds are independent
given all public data. Samples for different h in one round share descents and are analysed per h.

## E1: fruitless-by-construction QS base + h Dixon-descended rows (`python -m src.fruit.harness_e1`)

- N: Blum semiprimes (nfsfruit generator), classes blum1 (N ≡ 1 mod 8) and blum5 (N ≡ 5 mod 8), 48 and 64 bits,
  the first 15 seeds in each class satisfying the residue condition (60 instances).
- Base: QS relations x² − N, B-smooth, with nfsfruit Cor J1's β(x) = 0 (computed from N by Jacobi symbols). By Cor J1
  every base-only dependency is trivial.
- Base sizes M ∈ {Ct + 5, 2Ct, 4Ct}, Ct = 1 + #{p ≤ B : p = 2 or (N/p) = 1}.
- 20 rounds per (instance, level); each round 4 fresh descents (lemmaM QS special-q descent, r uniform on [1, N−1]);
  h ∈ {1, 2, 4}.
- Parameters: 48 bits B = 1000, B1 = 10⁶, J = 3000; 64 bits B = 4000, B1 = 1.6·10⁷, J = 6000.

Predictions:
- **P1** (h = 0): 0 fruitful base-only samples (exact).
- **P2** (forced): for each (level, h) cell with n ≥ 20, the fruitful rate is 1/2 with |z| ≤ 3.5; pooled over all cells
  of each h, |z| ≤ 3.
- **P3** (uniform): pooled fruitful count matches Σ_i p_i, p_i = ½(2^{dK} − 2^{dK−r})/(2^{dK} − 1), |z| ≤ 3.5.
- **P4** (identities): max(0, h − deficit) ≤ r ≤ h in every record; X² ≡ Y² in every congruence (asserted in code).
- **P5** (directional): Pr[r ≥ 1 | h = 1] increases with the base size: C+5 < 2C ≤ 4C (pooled over bits and classes).

## E2: cubic NFS base + h rows descended into the same field (`python -m src.fruit.harness_e2`)

- N: classes blum, mixed, both1 (semiprimes) and tri (3 primes), 40 and 44 bits, seeds 0–9 (80 instances).
- Base: monic base-m cubic NFS exactly as nfsfruit nfs3_base (40 bits: B = 1500, A = 1000, Bb = 400; 44 bits: B = 2500,
  A = 2600, Bb = 1300; 40 quadratic characters; excess 60).
- Descent: r uniform on [1, N−1], R = r², Gauss split u/v, ≤ 2 large primes per side in (B, B1] (B1 = 12000 / 20000),
  special-q in the same cubic field (lattice box ±48), rows with columns/characters shared with the base.
- 10 rounds per instance, 3 fresh descents per round, h ∈ {1, 2, 3}.

Predictions:
- **Q1** (forced): fruitful rate 1 − 2^{1−ω} (½ semiprimes, ¾ tri): |z| ≤ 3.5 in each (bits, ω-group, h) cell with
  n ≥ 20; pooled per ω-group, h = 1: |z| ≤ 3.
- **Q2**: X² ≡ Y² always (asserted); "nonsquare" (character failure) ≤ 5% of attempted congruences.
- **Q3** (uniform, lower bound): pooled fruitful count ≥ Σ (1 − 2^{−r})(1 − 2^{1−ω}) − 3.5σ.
- **Q4**: identities as P4.
- **Q5** (Prop SR): in every checked dependency the square roots computed via the three inert primes agree up to sign.
- **Q6** (exploratory, no prediction): Pr[r ≥ 1 | h = 1], deficits, new columns introduced by descents.

## E3: two-term dependencies (x-injectivity) in NFS relation sets (run inside E2, round 0)

On each E2 base set and on an nfsfruit nfs3_rand set (LV-style randomized cubic, nfsfruit parameters) for the same N:
- **R1**: no two relations have identical full parity vectors (rational + algebraic + characters): 0 pairs in every
  instance.
- **R2**: relations in some rational-square pair ((a−bm)(a′−b′m) a perfect square): ≤ 5% of relations, pooled.
  Trend 40 → 44 bits reported without a directional prediction (toy sizes are far from the asymptotic regime of
  Prop IV-LV).

Commands (after this commit):
```
PYTHONPATH=. .venv/bin/python -m src.fruit.harness_e1 experiments/configs/fruit_e1.json experiments/raw/fruit_e1.jsonl
PYTHONPATH=. .venv/bin/python -m src.fruit.harness_e2 experiments/configs/fruit_e2.json experiments/raw/fruit_e2.jsonl
PYTHONPATH=. .venv/bin/python -m src.fruit.analyze experiments/raw/fruit_e1.jsonl experiments/raw/fruit_e2.jsonl
```
