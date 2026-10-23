# prim: second-stage predictions (post hoc classes, fixed before the b = 16 post-hoc run)

Agent `prim`, 2026-10-23. Written after the registered E1 runs (b = 12, 14, 16) and after developing the post-hoc
classes on b = 12 and b = 14 only (`src/prim/e1_posthoc.py`), and BEFORE running `e1_posthoc --b 16`.
The registered predictions in `notes/prim_prereg.md` are unaffected and are reported as registered.

Post-hoc classes (developed on b ≤ 14): COSET (registered), FROZEN (all gate outputs at r equal those at an integer
|t| ≤ 64), F1 (gates with a pure base frozen at such t), CLOSURE (a non-constant gate input equals ±(product of the
node's integer constants with total exponent ≤ 2)), PIN (registered), OTHER. The comparison model is the registered
fresh-random model with the shift exclusion u ≠ c (`shift_model`).

Motivation for CLOSURE: exact identity a^r = −(a + 1), r odd ⇒ (a^r + 1)^r = a + 1; control with random bases a
(scratch, b = 14): unit-r solutions of (a^r + 1)^r = w are 2.03× the model for w = a + 1 and 0.87–1.13× for
w ∈ {a + 2, random, (a + 1)², (a + 1)^{−1}}.

Predictions for b = 16 (all 1491 primes; q-world = all unordered pairs via Theorem U):
- Q1: for each of the 8 nested affine nodes (D2a_1, D2a_2, D2a_3, D2b, D2c, D2d, D2e, D3) the full-r ratio
  OTHER / model-on-OTHER-residues lies in [0.8, 1.25].
- Q2: for each nested node whose q-world OTHER model mass is ≥ 20, q-world OTHER / model ∈ [0.6, 1.6].
- Q3: two-sided (non-split) zeros with both sides OTHER, summed over the nested nodes, ≤ 3·(model) + 5.

## Addendum (v2, after codex adversarial round prim_2): independent confirmation on an untouched range

Codex (prim_2 #21) correctly notes that b = 16 was partially exposed (the registered E1 run, including FROZEN-residual
summaries, preceded the second-stage registration). For a clean confirmation we use b = 15 (787 primes in
[3·2^13, 2^15)), on which no E1 or post-hoc computation of any kind has been run. Code is fixed at the v2 state of
`src/prim/e1_posthoc.py` (only change since v1: the q-world model mask now also excludes PIN residues, matching the
numerator; codex prim_2 #20). Predictions for b = 15, identical in form to Q1–Q3:
- Q1': full-r residual ratio ∈ [0.8, 1.25] for each of the 8 nested affine nodes.
- Q2': q-world residual ratio ∈ [0.6, 1.6] for each nested node whose q-world residual model mass is ≥ 20.
- Q3': two-sided OTHER/OTHER zeros summed over nested nodes ≤ 3·model + 5.
