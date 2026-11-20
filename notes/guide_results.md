# guide: is there any signal for smoothness beyond size and local divisibility? (2026-11-21)

**Motivation.** The user asked about guided growth: plant seeds, track densities with another system, and steer. The
α ≥ 1/3 model allows arbitrary adaptive guidance *under* H. Local-congruence guidance is the divisor trick (Theorem J,
one-for-one), and size is lattice geometry. So guidance can only help through a **global** smoothness correlation.
This experiment tests for one directly. A positive result would be the first evidence against H.

**Process.**
- Codex built and ran the pipeline (session 01a1fb1f-8343-70b0-8ee6-a79474139df5).
- Gemini independently wrote the feature and leakage specification (agy bcb713a3): `experiments/guide/guide_spec.md`.
- Round 2 implemented Gemini's full spec.
- Briefs: `codex/guide_{common,codex,gemini,2}.md`. Outputs: `codex/guide_*.out`.
- Code, pre-registrations, reports and summary tables: `experiments/guide/round{1,2}/`. The 3 GB of raw data and
  models were not committed.

## Verdict: null. No evidence of any global smoothness signal and none of gain growing with u.

**Round 1 (codex).**
- Data: 3,610,064 labelled NFS candidates on 104 balanced semiprimes, with 48 moduli and separate spatial regions held
  out.
- Baseline: y-rough cofactor sizes plus small-prime valuations. Test features: geometry and search history.
- Result at B = 262,139: 11,536 vs 11,536 relations at a 20% confirmation budget. Calibration put zero weight on every
  added predictor.
- Planted control: a 2.01 odds ratio was detected at 3.9–5.4 cluster standard errors.
- Shuffle controls behaved as expected.

**Round 2 (Gemini's spec, fresh held-out moduli).**
- Data: 1,077,330 candidates.
- Feature families G1–G8:
  - G1: coordinate S-units.
  - G2: continued fractions and root approximation.
  - G3: Minkowski conjugate balance, real sign chambers, covariants.
  - G4: sanitized shifts and Bezout/Farey neighbours.
  - G5: multiplicative and bilinear relatives.
  - G6: descent parent-seed inheritance.
  - G7: history, cosets and other candidates' digits.
  - G8: raw bits, carries and an MLP.
- Targets: the continuous valuation share W_s and the cross-side residual covariance.
- Leakage checks: L1 rejected 211 of 432 polynomial/shift combinations.
- **None of 108 comparisons passed the adjusted discovery rule** (smallest Holm p = 0.416). Refitting against a
  comparator given every local valuation and projective prime-power column gave the same verdict.
- **G6, descent inheritance ("do smooth parents seed smoother children?").** 320 matched special-q lattice pairs on 20
  fresh moduli. Adjusted child smoothness difference:
  - y = 2000: −0.017 pp, 99% CI [−0.10, +0.06];
  - y = 32768: +0.029 pp, 99% CI [−0.03, +0.11].

  A planted 4× parent effect was detected at both cutoffs. A 2× effect was detected on algebraic-side smoothness but
  missed on joint smoothness.

**Power caveat.** In the uncertain-cofactor stratum there were only 758 / 40 joint hits, so joint-target power is
weak, especially for G3–G5, whose planted-plant power was 0%. Algebraic-side power was good for G1, G2, G6, G7 and
G8.

**Arithmetic positive control.** On 24 extra fresh moduli, forced norm splitting raised algebraic smoothness from 1.36%
to 5.14%. G8 detected it (p = 0.000676); G2 and G3 missed it.

**Bottom line.** At toy sizes the smoothness landscape behaves as H predicts: a golf course beyond size and congruences.
Moderate joint effects for G3–G5 are not excluded.
