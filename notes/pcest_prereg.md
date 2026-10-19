# pcest pre-registration (phase 8, 2026-10-19)

Frozen code: `src/pcest/{geometry,params,toy,e1}.py`, tests `tests/test_pcest.py`.
Config: `experiments/configs/pcest_e1.json`. Output: `experiments/raw/pcest_e1.jsonl`.
This file and the code are committed before E1 is run.

## E1: toy pair correlations in designs 1D (g = x − m) and 2D (g = g₁x − g₀)

**Purpose.** A plausibility and correctness check, not evidence for any asymptotic claim. The
claims of `notes/pcest_results.md` are proved or conjectured statements about N → ∞ at
smoothness L[1/2]; a toy at 56–64 bits can test only that (i) the exact lattice identities hold in
running code, (ii) the 2D-randomized rational side does not change yields or correlations in
practice, and (iii) nothing pathological appears in the pair correlations.

**Design.**
- Instances: balanced semiprimes ("semi") and 3-prime products ("tri") at 56 and 64 bits, seeds 0
  and 1: 8 instances. Same ℓ-list (30 primes in (2¹¹, 2¹²]) for both designs of an instance.
- d = 3, M = ⌊N^{1/4}⌋, R coefficients uniform in [−1024, 1024].
- 1D: m uniform in [M/2, M]; f = (centred base-m expansion) + (x − m)R.
- 2D: (g₀, g₁) uniform coprime in [M/2, M]²; f = gbe(N, 3, g₀, g₁) + (g₁x − g₀)R (Lemma GB).
- Special-q box: hypD's reduced-basis box with S = 16 (≤ 561 points), one sign per ± pair,
  primitive, b ≠ 0. Rational smoothness bound 256 on G/ℓ, algebraic bound 16384 on F. (Two
  bounds, as in LV's (β, β′), so that the rational event is not near-certain at toy sizes.)
- 40 fields per (instance, design); each field × ℓ is one cell. Statistics over cells:
  PC_rat = E[N_r(N_r − 1)]/(E N_r)²; PC = E[N(N − 1)]/(E N)²;
  PC_alg = (E[N(N−1)]/E[N_r(N_r−1)])/(E N/E N_r)² (a conditional proxy); γ = Pr[N > 0].

**Predictions** (registered; "holds" means in at least 7 of 8 instances unless stated).
- P1: PC_rat(2D) ∈ [0.8, 1.5].
- P2: PC(2D) ∈ [0.8, 2.0].
- P3: |PC_rat(2D) − PC_rat(1D)| ≤ 0.25.
- P4: E N(2D)/E N(1D) ∈ [0.7, 1.4] (Lemma GB: sizes unchanged up to constants).
- P5 (exact, all instances): 0 violations of ℓ | G(P) and F(P)g₁^d ≡ N b^d (mod G(P)) on every
  tested rational-smooth point; 0 violations in the Lemma Q checks (20 per 2D record: index
  q = |Δ|/ℓ, graph form, and membership of the actual (G/ℓ, G′/ℓ)).
- P6 (reported, no threshold): γ(2D) versus γ(1D), and PC_alg(2D) versus PC_alg(1D).

**Commands.**
```
.venv/bin/python -m src.pcest.e1 experiments/configs/pcest_e1.json experiments/raw/pcest_e1.jsonl 8
.venv/bin/python -m src.pcest.e1 --analyze experiments/raw/pcest_e1.jsonl
```
Any deviation will be disclosed in the results note.
