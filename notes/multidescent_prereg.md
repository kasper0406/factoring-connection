# multidescent — pre-registration (E1–E3)

Agent `multidescent`, phase 25, 2026-10-23. Committed before any main run. Code frozen at this commit:
`src/multidescent/{toy,experiments,analyze,model}.py`. All experiments are DIAGNOSTIC toys (48-bit N, d = 3, 5):
they test the shape of the children's law and of conditional goodness at toy size and support nothing rigorous.

**Disclosed pilots (calibration only, not used as data).** Seed 9001 (instance) and seeds 901 (E1–E3 harness) were
used to choose box sizes, factor bases and the E2 pool size so that rates are neither 0 nor 1, and to choose the bands
below. Pilot values: E1 algebraic sub-bucket ratios vs the 1/q label law 0.93–1.11; E2 diff −0.004 to +0.003, zero
fractions 0.02–0.22; E3 (d = 3) alg/null 1.03 in the within-budget stratum, line excess −0.007, fiber excess 0.24;
(d = 5) alg/null 0.99, line excess 0.05, fiber excess 0.05.

**Instances.** `balanced_semiprime(48, s)` for s = 1..6; two independent main fields per instance (design (a‴) toy
fields, `toy.sample_field`), seeds as in `experiments.py`. Commands:
```
.venv/bin/python -m src.multidescent.experiments E1 E1   experiments/raw/multidescent_e1.jsonl 6
.venv/bin/python -m src.multidescent.experiments E2 E2   experiments/raw/multidescent_e2.jsonl 6
.venv/bin/python -m src.multidescent.experiments E3 E3d3 experiments/raw/multidescent_e3d3.jsonl 6
.venv/bin/python -m src.multidescent.experiments E3 E3d5 experiments/raw/multidescent_e3d5.jsonl 6
.venv/bin/python -m src.multidescent.analyze E1|E2|E3 <file>
```
Configurations are the dict `CONFIGS` in `experiments.py` (frozen).

## E1 — the children's law (level-1 → level-2)

Level-1 rational targets ℓ uniform among primes in (2^17, 2^18] (300 per field); ALL box points (S₁ = 16) whose
cofactors G/ℓ and F are B₂-smooth (B = 2^11, B₂ = 2^16) are accepted; children = rational primes and simple
degree-one prime ideals of norm in (B, B₂] (exponent 1). Windows (2^j, 2^{j+1}], j = 11..15, each split into 4
equal-log-width sub-buckets. References: "uniform on labels" (each available label — prime, or root of f mod q —
weight 1) and "1/q labels" (label weight 1/q), both pooled over fields with weight = the field's number of nodes.

- **P1a.** Algebraic children: every sub-bucket with ≥ 150 children has (empirical share)/(1/q-label share) in
  [0.85, 1.15].
- **P1b (domination).** Algebraic children: max over sub-buckets (≥ 150 children) of (empirical share)/(uniform-label
  share) ≤ 1.6, and min ≥ 0.6.
- **P1c.** Rational children: every sub-bucket with ≥ 100 children has ratio vs the 1/q law in [0.7, 1.3].

## E2 — rational children in a pool of fresh fields vs uniform primes

M = 2^22, ℓ in (2^16, 2^17] (400 per field), first accepted point per ℓ; for each rational child q a log-matched
uniform prime q′ ∈ [q/1.05, 1.05q] (q′ ≠ q); pool of 24 fresh fields (same law, independent of the main field);
p_q = fraction of pool fields in which q has a box point (S₂ = 36) with both cofactors B-smooth. Cluster bootstrap
over parents (2000 resamples).

- **P2a.** 95% CI of mean(p_child − p_unif) ⊂ [−0.03, 0.03].
- **P2b.** 95% CI of (zero-fraction child − zero-fraction uniform) contains 0 and ⊂ [−0.06, 0.06].

## E3 — algebraic children under the conditional field law (constrained resampling) and the rational block

For a parent P₀ (first accepted point with ≥ 1 algebraic child; 60 parents per field for d = 3, 40 for d = 5), the
fiber family F_H = F + G·L₀·H (H a binary form of degree d − 2, coefficients in [−Hh, Hh]; all 81 forms for
d = 3, Hh = 4; 81 sampled of 625 for d = 5, Hh = 2) keeps F(P₀), F(g₀, g₁) = Nw and every child label. Fiber fields
in which some label becomes degenerate (q | f_d·disc) are skipped for that parent (disclosed). For each fiber field:
per-child goodness (box S₂, both cofactors B-smooth) and the counts N_𝔮. Null: rational targets with norms matched
(×1.05) to the children, valid in every fiber field (g unchanged): same field quality, no shared congruence with P₀.
J = Pr_H[all children good], Π = ∏ Pr_H[child good]; R_joint = ΣJ/ΣΠ per stratum (m ≤ d − 1 = within the fiber
rank; m > d − 1 = overdetermined); alg_over_null = R_joint(children)/R_joint(null). E3r (rational block): for each
child box, R(j) = #{P′ in the box : G_j(P′) B-smooth} along the line g_j = g + j·P₀, |j| ≤ 40 (G_j(P₀) = G(P₀)).
Excess relative variance = Var/μ² − 1/μ (0 for Poisson), over boxes with μ ≥ 0.5.

- **P3a.** Fiber algebra: 0 violations (F_H(P₀) = F(P₀), F_H(g) = Nw, labels are roots of f_H mod q).
- **P3b.** Within-budget stratum (both d): alg_over_null 95% CI ⊂ [0.85, 1.2].
- **P3c.** Overdetermined stratum (d = 3, m ≥ 3; d = 5 has too few such parents and is not tested): alg_over_null
  point estimate in [0.67, 1.5] (no detectable overdetermination effect beyond field quality).
- **P3d.** Rational block along lines: median excess relative variance in [−0.1, 0.15] (both d).
- **P3e.** Fiber (algebraic) block: median excess relative variance in [−0.1, 0.4] (d = 3) and [−0.1, 0.2] (d = 5);
  median Pr_H[N_𝔮 = 0]/e^{−μ} in [0.8, 1.25] over children with e^{−μ} > 0.02 (both d).

**Interpretation fixed in advance.** P1 supports the "weighted uniform" description of the children's law (the label
weight ∝ 1/q, i.e. domination by uniform up to a factor ≈ 2 per dyadic window). P2 supports Prop RC's mechanism
(rational children are as good as matched uniform primes in fresh fields, including the lower tail). P3b–c test the
heuristic that conditional (fiber) goodness of the children factorises up to field quality, also beyond the fiber rank
— a toy check of the premise behind both Hyp LS_C and Hyp PC₃, not evidence for either asymptotically. P3d–e measure
how close child-box counts are to Poisson under the conditional laws (the quantity Prop UB needs to be small); toy
local factors are expected to make the algebraic excess positive. A failed prediction is reported as failed.

## E0 — exact model arithmetic (not a prediction)
`model.py`: α_j = (j+1)/(3j+2) and its brute-force grid check, the 3/8 parameters, and the union-bound precision
margins; covered by `tests/test_multidescent.py`.
