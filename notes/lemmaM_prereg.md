# lemmaM pre-registration (frozen before any E1/E2 run)

Agent `lemmaM`, 2026-10-18. Code frozen with this file: `src/lemmaM/{descent,harness_e1,unit_twist,run_e2,analyze}.py`,
configs `experiments/configs/lemmaM_e{1,2}.json`, tests `tests/test_lemmaM.py`.
Before this commit only the tests and two smoke runs (36 bits, two seeds; 64 bits, one seed with seed 99, which is not in the
experiment range 0..19) were executed. No experiment instance has been run.

## E1. Dixon descent (toy: QS special-q descent)

**What it tests.** Theorem DD (notes/lemmaM_results.md): a relation r² ≡ T (mod N), obtained by descending the uniformly random
square R = r² mod N, where every computation except the final factor r depends on R alone, makes *every* dependency that contains
it fruitful with probability exactly 1 − 2^{1−ω}. There is no Blum condition, no Conjecture F, and no assumption about the
other relations. The toy descent splits R ≡ u/v with a Gauss-reduced (u, v), requires u and v to be B1-smooth (B1 = B²), and
descends each large prime ℓ ∈ (B, B1] with (N/ℓ) = 1 through a QS special-q relation x² − N = ℓ·c with c B-smooth. This is
a correctness check of the theorem and its bookkeeping. It makes no complexity claim, because the special-q values have QS size √N.

**Instances.** `src/nfsfruit/instances.composite`, classes blum / mixed / both1 / tri, bits 48 / 56 / 64 / 72, seeds 0..19.
**Modes.**
- pure: only descent relations, #cols + 10 of them.
- mixed: #cols − 10 plain QS relations plus 20 descent relations.
- plain: Farey-reconstruction Dixon with B1 = B (no descent), 40 relations, used for cost counting only.
Per instance, 60 uniform random nonzero kernel elements. In mixed mode, 60 more come from the kernel of the QS rows alone.
Parameters as in the config.

**Predictions.**
- P1: every congruence satisfies X² ≡ Y² (mod N). 0 non-squares.
- P2 (pure mode): fruitful fraction 1/2 (semiprimes, pooled over blum/mixed/both1) and 3/4 (tri). Cluster-robust |z| ≤ 3 in each
  of the two groups. No Blum/non-Blum difference; reported per class.
- P3 (mixed mode, vectors containing ≥ 1 descent relation): the same as P2. Vectors with no descent relation: reported only.
  Heuristically 1/2 and 3/4; no theorem applies.
- P4: flip check. For each descent-containing vector, recompute the descent for (εr_j)² with a nontrivial ε (harness). The
  output (W, T) must be identical and the class of X/Y must flip. Predicted pass rate 100%.
- P5 (cost bookkeeping, heuristic): r-trials per accepted relation, plain / descent, is ≥ 5 at 64 and 72 bits. Special-q
  success fraction is reported.
Failure of P1, P2, P3 or P4 would indicate a bug or a false theorem. P5 is a heuristic expectation.

**Commands.**
```
.venv/bin/python -m src.lemmaM.harness_e1 experiments/configs/lemmaM_e1.json experiments/raw/lemmaM_e1.jsonl 8
.venv/bin/python -m src.lemmaM.analyze e1 experiments/raw/lemmaM_e1.jsonl
```

## E2. Unit-twist census (candidate β-flipping involution for QS / CFRAC)

**What it tests.** In the real-quadratic picture of QS/CFRAC, the relation set of convergents of √N carries the map k ↦ k + P
(P = period), multiplication by the fundamental unit ε = a + b√N. It keeps Q_k (so the square-class vector) and multiplies the
root A_k by a (mod N), a² = (−1)^P. When a² = 1 and a ≢ ±1, this is exactly a Dixon ε-involution, i.e. x-preserving and β-flipping. It moves
every relation by a full period.

**Instances.** classes blum / mixed / both1, bits 28 / 32 / 36 / 40, seeds 0..19, period cap 2·10^7.
**Predictions.**
- P6: the identities Q_{k+P} = Q_k, A_{k+P} ≡ a·A_k and a² ≡ (−1)^P hold in 100% of instances.
- P7: the fitted slope of log₂ P against bits lies in [0.40, 0.55], i.e. the involution's displacement is ≈ √N (exponential).
- P8 (reported, weak expectation): the fraction of instances with a nontrivial a (a² = 1, a ≢ ±1) is high for blum
  (expected ≥ 0.8, after Legendre's classical treatment of x² − pqy² = 1). Norm −1 units occur only in the both1 class.
  Mixed is reported without a prediction.

**Commands.**
```
.venv/bin/python -m src.lemmaM.run_e2 experiments/configs/lemmaM_e2.json experiments/raw/lemmaM_e2.jsonl 8
.venv/bin/python -m src.lemmaM.analyze e2 experiments/raw/lemmaM_e2.jsonl
```
