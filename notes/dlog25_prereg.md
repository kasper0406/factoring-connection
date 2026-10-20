# dlog25 pre-registration: E1 (Lemma S_p″) and E2 (Lemma XK4-U at d = 3, 5, 7)

Agent `dlog25`, phase 13, 2026-10-20. Frozen at the commit that adds it, together with `src/dlog25/{params,toy_sp,toy_xk}.py`
and `tests/test_dlog25.py`. Calibration (before freezing, not part of E1/E2): toy_sp on p = 10007, 40009; toy_xk seeds
900–903 (d = 3, 5, 7, 9). The calibration showed that at d = 7, 9 the toy pair boxes give too few kernel vectors, so d = 7
is run as a descriptive block and d = 9 is dropped.

## What is tested (correctness of two exact lemmas only; no asymptotic claim)

- E1 tests the split-counting identity behind Lemma S_p″ (notes/dlog25_results.md §3): for T running over all of F_p^*,
  the Gauss-reduced vector of {u ≡ Tv (mod p)} is ±(u, v) with u, v ∈ S coprime for exactly P′ := #{(u, v) ∈ S²: gcd = 1}
  values of T. Exhaustive over T.
- E2 tests Lemma XK4-U (unit elimination with the EXPLICIT precision κ and relation bound R of §4, not an ad hoc
  precision) at d = 3, 5 (and 7 if the toy produces kernels). Toy deviations: units come from NFS-type pair elements
  a − bα of monic small-coefficient fields with good-prime-only smooth norms, not from descended rows; the kernel basis is
  LLL-reduced before use; Miller–Rabin/sympy isprime for auxiliary primes.

## Blocks

- E1: the 8 cases `toy_sp.CASES` (p ∈ {200003, 1000003, 2000003, 4000037}, k ∈ {1, 2, 3}).
- E2: `toy_xk.BLOCKS`: d = 3 seeds 0–7; d = 5 seeds 10–17; d = 7 seeds 20–23 (descriptive).

## Predictions (frozen)

- **P1 (E1, identity).** identity_ok in 8/8 cases (hits = P′ exactly).
- **P2 (E1, bound).** In every case with |G| ≥ 4k³: P′ ≥ |S|²(1 − 2k³/|G|) and hits/(p − 1) ≥ |S|²/(2(p − 1)).
- **P3 (E2, certification).** In every d ∈ {3, 5} run with status ok: the LLL first vector has n ≠ 0 and is certified
  by the Blanksby–Montgomery half-gap at the explicit κ (certified = True). Any failure refutes the explicit κ.
- **P4 (E2, torsion).** Every certified x* has ψ_q(γ_{x*}) ∈ {±1} at all 5 auxiliary 40-bit primes with a consistent sign
  (μ_K = {±1}, d odd): 0 violations.
- **P5 (E2, units / rank).** 0 non-unit kernel vectors; unit_rank_implied ≤ r_U in every ok run; 0 LLL rows that are
  certified non-relations pass the auxiliary torsion test.
- **P6 (E2, coverage).** status ok in ≥ 6/8 runs at d = 3 and ≥ 6/8 at d = 5. d = 7: descriptive, no prediction.
- **P7 (descriptive).** log₂‖n‖ against log₂R and κ against d.

## Commands

```
.venv/bin/python -m src.dlog25.toy_sp experiments/raw/dlog25_e1.jsonl
.venv/bin/python -m src.dlog25.toy_xk experiments/raw/dlog25_e2.jsonl
.venv/bin/python -m src.dlog25.params > experiments/raw/dlog25_params.txt
.venv/bin/python -m pytest -q tests/test_dlog25.py
```
