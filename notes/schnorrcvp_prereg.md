# schnorrcvp pre-registration (phase 19)

Agent `schnorrcvp`, 2026-10-21. Written and committed before E1 and E2 were run; only quick dry runs on
scratch outputs (`--quick`, ≤ 2 instances) were used to debug the code. Code: `src/schnorrcvp/`.

## Disclosure: model numerics computed before this file

During development I computed the H_S model total count E_tot(N, n) (all fac-relations of any size, see
below) at five points:

| bits | n | E_tot |
|---|---|---|
| 40 | 47 | 1.9e14 |
| 200 | 47 | 6.6e−11 |
| 400 | 47 | 1.3e−52 |
| 400 | 95 | 3.0e−31 |
| 800 | 95 | 1.5e−115 |

E0 below re-runs these deterministically and extends them. It is a model computation, not an experiment;
it is listed here so that the sensitivity checks are fixed in advance.

I also ran two toy E1 boxes while debugging (24-bit n=10 and 32-bit n=15, seed 0): observed primitive 7 and
4, against H_S′ predictions 3.7 and 1.5. These two instances are part of the E1 grid and will be re-run.

## Definitions

- Factor base P = first n primes. A fac-relation is (u, v, w) with u = vN + w, u, v, |w| all P-smooth, w ≠ 0.
  It is *primitive* if gcd(v, w) = 1, which makes u, v, w pairwise coprime. Lattice vectors always have
  gcd(u, v) = 1.
- **H_S (naive).** u = vN + w is P-smooth with probability d(vN), the local density of P-smooth integers.
- **H_S′ (refined, primitive).** For coprime smooth (v, w), u is coprime to Q = supp(vw), and residues
  modulo the other primes are uniform. So Pr = dens{smooth m near vN with gcd(m, Q) = 1} / ∏_{q∈Q}(1 − 1/q),
  evaluated by inclusion–exclusion (`psi.expected_primitive_in_box`).
- **H_res (naive)** for a lattice output (u, v, w): Pr[w smooth] = d(|w|).
- **H_res′.** w is coprime to Q = supp(uv) (proved: gcd(uv, w) = 1). Pr = d(|w|)·∏_{q∈Q}(1 − f_q)/(1 − 1/q),
  where f_q = d(|w|/q)/(q·d(|w|)) (independence-product form; `run_e2.Res.refined`).
- Densities come from exact integer bracket tables (`psi.PsiBracket`) using midpoints. Floating point is
  used only for model predictions.

## E0: model numerics (deterministic)

- E_tot(bits, n) via `psi.expected_relations` (H_S, all v with ln v ≤ 1.5 ln N + 50, |w| < vN, both
  signs), for bits ∈ {64, 128, 200, 256, 400, 512, 800, 1024}.
- n_min(bits): the least n on a grid (step 5 for n < 200, step 25 above) with E_tot ≥ n + 1.
- Sensitivity check (pre-registered): the resolution change K = 50 → 100 moves log10 E_tot by ≤ 1 at
  (400, 47) and (800, 95), and lnv_max = 1.5 ln N + 50 versus 3 ln N changes it by ≤ 1%. If either check
  fails, I report the bracket range instead of a point value.
- Expectation (not a test): E_tot(400, 47) ≪ 1 and E_tot(800, 95) ≪ 1, so heuristically *no*
  fac-relations exist at Schnorr's claimed parameters.

## E1: exact relation counts in boxes vs H_S / H_S′

- Grid (bits, n, V, W):
  - (24,10,1e3,1e5), (28,10,1e3,1e5), (32,10,1e3,1e5)
  - (24,15,1e3,1e5), (32,15,1e3,1e5), (40,15,1e3,1e5)
  - (24,25,1e3,2e4), (32,25,1e3,2e4), (40,25,1e3,2e4)
- Instances: balanced semiprimes, seeds 0–9.
- Exact enumeration of every relation with v ≤ V, 1 ≤ |w| ≤ W (`enum_rel.relations_in_box`).
- Command: `.venv/bin/python -m src.schnorrcvp.run_e1` → `experiments/raw/schnorrcvp_E1.jsonl`.
- Predictions:
  - **E1-P1.** Pooled over all cells, Σ obs_primitive / Σ pred_primitive ∈ [1/2, 2]. Per cell (10 seeds
    pooled), the ratio is in [1/3, 3] whenever the pooled prediction is ≥ 10.
  - **E1-P2.** Pooled per cell, pred_naive ≥ obs_primitive in every cell. That is, naive H_S over-predicts
    primitive relations, which supports using it as a heuristic *upper* estimate in E0.

## E2: CVP-output relation yield vs H_res / H_res′

- Grid:
  - bits ∈ {20, 24, 28, 32, 36, 40}, 3 balanced semiprimes each (seeds 0–2);
  - n ∈ {12, 20, 30, 40};
  - five lattice variants:

    | Variant | Diagonal | Scale |
    |---|---|---|
    | S21perm | f(i), a permutation | C = N^{1/(n+1)} (Schnorr 2021) |
    | S21sqrt | √f(i) | C = N^{1/(n+1)} |
    | A0.5 | √(ln p_σ(i)), randomly permuted | C = N^{0.5} |
    | A1.0 | same | C = N^{1.0} |
    | A1.25 | same | C = N^{1.25} |

  - 10 random permutation seeds per cell.
- Outputs per lattice:
  - Babai after LLL (approximate oracle);
  - Babai after BKZ-10;
  - the up-to-20 closest vectors by exact enumeration (mpfr GSO, exactly re-sorted). list0 is the exact
    CVP answer.
- Command: `.venv/bin/python -m src.schnorrcvp.run_e2` → `experiments/raw/schnorrcvp_E2.jsonl` and
  `schnorrcvp_E2_factor.jsonl`.
- Predictions:
  - **E2-P1.** Pooled over all outputs where p_ref is defined, Σ rel / Σ p_ref ∈ [1/2, 2]. Per variant
    (pooled), the ratio is in [1/3, 3] where Σ p_ref ≥ 10.
  - **E2-P2.** The refined model fits better than the naive one: |log(Σrel/Σp_ref)| < |log(Σrel/Σp_naive)|,
    pooled.
  - **E2-P3.** For each (n, variant), the relation yield per output at 40 bits is below that at 20 bits, in
    every cell with ≥ 1 relation at 20 bits.
  - **E2-P4.** Schnorr-2021 scaling fails to target small residuals. For S21perm and S21sqrt at
    bits ≥ 32, the median log2|w|/bits over list outputs is ≥ 0.5.
  - **E2-P5 (descriptive).** Duplication: the number of distinct exact-CVP (list0) relations across the 10
    permutation seeds, versus the number of list0 outputs that are relations, per cell.
  - **E2-P6 (correctness).** For each N with > |P_40| + 1 distinct pooled relations, the GF(2) congruence of
    squares returns a factor of N. Every reported factor divides N. Factors are STRUCTURAL: N's prime
    factors exceed p_40 = 173, so trial-division hits are impossible.

## Resources

- Single process, RSS ≤ 2 GB (expected < 300 MB). memory_pressure is checked first. E1 and E2 run
  sequentially.
- Any deviation will be disclosed in the results note.
