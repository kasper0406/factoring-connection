# nested: pre-registration of E1 (toy endogeneity test)

Agent `nested`, phase 12, 2026-10-20. Frozen before the E1 run. E1 is a toy-size *diagnostic*. It cannot test any
asymptotic claim and it proves nothing. It asks one empirical question: in one fixed field, are the *endogenous*
children of a special-q relation (the prime factors of the chosen box point's two values) one-step descendable at
the same rate as *uniformly drawn* targets of the same type and size? A clear penalty would be evidence that the
"good set" of a fixed field is adversarial to its own children. No penalty is what the standard descent heuristic
predicts, and it would not prove anything.

**Code (frozen with this file):**
- `src/nested/toy_e1.py`: field construction (imports `src/thmA/numfield.py` unchanged: `make_f`, `hom_eval`,
  `roots_mod`, `irreducible_mod`, `disc`), the one-step test D1, the parent search, the matched baseline.
  It never factors N and never uses any factor of N: N enters only through `make_f`.
- `src/nested/harness_e1.py`, `src/nested/analyze_e1.py`; config `experiments/configs/nested_e1.json`.
- Output `experiments/raw/nested_e1.jsonl` and `experiments/raw/nested_e1_summary.json`.

**Definitions.**
- Field A = (g₀, g₁, w, c) as in design (a‴_k; μ) at toy size: degree d = 3, M = 2^{⌈n/4⌉} (n = bit length of N),
  (g₀, g₁) uniform in [M/2, M]² coprime, w a uniform prime in (2⁵, 2⁶], c uniform in [−H, H]³ with H = 4;
  f = `make_f(N, w, 3, g₀, g₁, c)`. Reject and resample if the content of f is ≠ 1 or no prime r ≤ 200 certifies
  irreducibility. G(a, b) = g₁a − g₀b, F(a, b) = Σ f_i a^i b^{3−i}.
- Targets: *rational* primes q, lattice Λ = {(a, b) : q | G(a, b)}; *algebraic* degree-1 prime ideals 𝔮 = (q, r),
  r a simple root of f mod q with q ∤ f₃·disc(f), lattice Λ = {(a, b) : a ≡ rb (mod q)}.
- Box of a target of norm q: the points P = (a, b) ∈ Λ with |a| ≤ V, 0 < b ≤ V, gcd(a, b) = 1, G(P) ≠ 0, F(P) ≠ 0,
  where V = ⌊√q·S⌋.
- **D1(τ)** (one-step descendable): some box point has both cofactors B-smooth, where the cofactor on the target's
  side is the value divided by q once, and the other side is the full value.
- **Parent node.** For a rational prime ℓ (the analogue of a level-1 target), scan its box (radius ⌊√ℓ·S₁⌋) in the
  order (b, a) and take the first P whose cofactors G(P)/ℓ and F(P) are both Y-smooth and which has at least one
  child. The children are the rational primes p ∈ (B, Y] dividing G(P)/ℓ, and the degree-1 prime ideals
  (p, a·b^{−1} mod p) for the primes p ∈ (B, Y] dividing F(P) with p ∤ b·f₃·disc(f). (Children with p | b·f₃·disc f
  are counted separately as degenerate and excluded.)
- **Matched uniform baseline.** For each child of norm p, draw one target of the same type uniformly from the targets
  with norm in [p/1.25, 1.25p] ∩ (B, Y] (algebraic: uniform over the admissible pairs (p′, r′)). Draws use a private
  stream, independent of the parent search.
- **Independent-field control.** A second field A′ is drawn independently (same law). Each rational child is also
  tested for D1 in A′, together with its own matched uniform draw in A′.

**Parameters.** n = 44 bits, d = 3, B = 2^{10}, Y = 2^{14}, S as fixed by calibration (default 16), S₁ = S,
parents ℓ = the primes in [2^{22}, 2^{22} + 2^{16}) in increasing order, up to 120 parent nodes per instance.
12 balanced semiprime instances `balanced_semiprime(44, seed)`, seeds 1–12.

**Calibration (allowed, disclosed).** Before the main run I may run seeds 900–903 only to choose S ∈ {8, 12, 16, 24,
32} so that the uniform algebraic D1 rate lies in [0.2, 0.8], and to bound the run time. Calibration seeds are
excluded. The chosen S and the calibration numbers will be reported. Nothing else may be changed after calibration.

**Predictions (pooled over the 12 instances; 95% CIs by bootstrap over parent nodes, 2000 resamples, seed 0).**
- **P1 (algebraic children, same field).** The 95% CI of rate_endo − rate_unif lies inside [−0.08, +0.08].
- **P2 (rational children, same field).** The same, inside [−0.08, +0.08].
- **P3 (within-node independence).** Over parent nodes with ≥ 2 non-degenerate children, the observed fraction of
  nodes all of whose children are D1 divided by the mean over nodes of ∏_children p_type (p_type = the pooled matched
  uniform rate for that child's type) lies in [0.75, 1.33].
- **P4 (rational children, independent field A′).** The 95% CI of rate_endo(A′) − rate_unif(A′) lies inside
  [−0.08, +0.08].

**Descriptive only (no prediction):** the distribution of the number of children per parent node, and the fraction of
nodes with more than 6 children (6 = the number of field parameters (g₀, g₁, w, c) at d = 3); the dead fraction
1 − D1 rate of uniform targets as a function of S ∈ {S/2, S, 2S} on the first 3 instances; the rate of degenerate
children; run time.

**Interpretation fixed in advance.** If P1–P4 hold, E1 is consistent with "no endogeneity penalty at toy size", which
is what the standard heuristic predicts and what Hypothesis AD (notes/nested_results.md) would need. It does not
support any rigorous statement. If P1 or P3 fails with a *negative* difference, that is evidence of an endogeneity
penalty in a fixed field, and it will be reported as such.

## Addendum A (2026-10-20, after calibration, before the main run)

Disclosed deviations, fixed before any main-run seed was executed:
1. **M = 2^{20}** (config `M_bits = 20`) instead of M = 2^{⌈n/4⌉} = 2^{11}. With M = 2^{11} the rational cofactor
   G(P)/ℓ ≈ M·S/√ℓ ≈ 2^{4} lies below B, so no parent had a rational child (calibration at S = 16: 0 rational
   children in 4 instances) and P2/P4 would be vacuous. With M = 2^{20} the rational cofactor is ≈ 2^{13–14}, inside
   (B, Y]. This is the design-(a‴) freedom μ (base M = N^μ with μ ≈ 0.45 here).
2. **max_parents = 400** instead of 120, so that the rational-children sample is large enough for the ±0.08 window
   (calibration: about 0.3 rational children per parent).
3. **S = 32**, chosen by the pre-registered rule. Calibration (seeds 900–903, M = 2^{20}, uniform algebraic D1 rate):
   S = 4: 0.05–0.09; S = 8: 0.08–0.14; S = 16: 0.15–0.29; S = 24: 0.19–0.37; S = 32: 0.22–0.48 (mean 0.36). S = 32 is
   the only listed value whose four calibration rates all lie in [0.2, 0.8]. Run time ≈ 3.5 s per instance at
   120 parents.
Calibration file: `experiments/raw/nested_e1_calib.jsonl` (kept; excluded from E1). Predictions P1–P4, the analysis
and all other definitions are unchanged.
