# resbarrier pre-registration (frozen before any experiment is run)

Agent `resbarrier`, phase 16, 2026-10-20. Code: `src/resbarrier/`. Raw output: `experiments/raw/resbarrier_E*.jsonl`.
The experiment is exact (integer arithmetic). Floating point appears only in printed summaries.
It is a correctness/sharpness check of a counting statement in notes/resbarrier_results.md. It cannot
show a speed-up, and a pass is not evidence for any algorithm.

## E1. Rank-2 count lemma (Lemma C2) on real polynomial pairs

**Statement under test.** For coprime integral binary forms F₁, F₂ (degrees d₁ + d₂ = D ≥ 3, nonzero
discriminants) with N | Res(F₁, F₂), and every E ≥ 1, X ≥ 1:

  #{(a, b) ∈ Z², |a| ≤ E, 0 ≤ b ≤ E : 0 < |F₁F₂(a, b)| ≤ X}
      ≤ 3B(1/3,1/3) · N^{−2/(D(D−1))} · X^{2/D} + 2(D+1)(2E+1).

The first term uses the Bean area bound A_F ≤ 3B(1/3,1/3)|disc F|^{−1/(D(D−1))} (cited; to be verified by the
literature critic) and |disc(F₁F₂)| ≥ Res² ≥ N². The comparison is done exactly: with c = 16 > 3B(1/3,1/3) = 15.8997…
and S = count − 2(D+1)(2E+1), check S ≤ 0 or (S/c)^{D(D−1)} · N² ≤ X^{2(D−1)}. (As a separate diagnostic we also
report the sharper exact test with |disc(F₁F₂)| in place of N².)

**Instances (seed 20261020).** N = product of two random primes of n/2 bits, n ∈ {40, 56, 72}.
- GNFS: base-m pairs, d ∈ {2, 3, 4}, m = ⌊N^{1/(d+1)}⌋ (+ small shift until the base-m digits give a pair with
  nonzero discriminants), g = x − m. 3 values of N per (n, d).
- SNFS: N = 2^e − 1 times nothing (no factor structure is used), f = x^d − 2^{d·⌈e/d⌉ − e}, g = x − 2^{⌈e/d⌉},
  e ∈ {40, 56, 72}, d ∈ {2, 3, 4}.
- Box E = 256. X ranges over 2^k for k from ⌈log₂ min nonzero |F|⌉ up to log₂ max |F| on the box, step 4.

**Predictions.**
- P1 (theorem check): 0 violations of the displayed bound over all instances and all X. Any violation means
  Lemma C2 or the cited area bound is wrong as stated; that would be reported as a refutation of C2.
- P2 (saturation): define the saturation s = (count − 2(D+1)(2E+1))₊ / (main term). Over the X grid take
  s_max per instance. Prediction: at n = 72, d = 3, median s_max(SNFS) ≥ 8 · median s_max(GNFS). Rationale:
  GNFS discriminants exceed N² by about N^{(2d−2)/(d+1)}, which lowers its count by a factor
  ≈ N^{−2(d−1)/((d+1)·D(D−1))} ≈ 2^{−6} at n = 72, d = 3, whereas SNFS discriminants are N²·2^{O(e/d·…)}.
  A failure of P2 would mean the box/shape terms dominate at these sizes, not that C2 is wrong.

## E2 (dropped before running)
An E2 testing "G for almost all N" by exhaustive pair enumeration was drafted and dropped before any run: at
N₀ ≤ 2^24 the resultant constraint H_f·H_g^d ≳ N₀ already forces |P| ≳ N₀ for every cell we could enumerate, so the
experiment cannot separate G from Lemma RES at feasible sizes. Theorem GAE is stated with its (elementary) proof only.

## Rules
- Commands and raw JSONL are recorded in the results note. Deviations from this file are disclosed there.
- No factoring routine is used on N in E1 (the N of the GNFS cells are built from their primes only to fix the
  instance; the counting code sees F₁, F₂ and N).
