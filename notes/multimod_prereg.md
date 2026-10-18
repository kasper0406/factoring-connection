# multimod: pre-registration (frozen before the pre-registered runs)

Agent `multimod`, 2026-10-18. This file is committed together with the code it
references (`src/multimod/{dcl,operational,experiments}.py`,
`tests/test_multimod.py`) before any run of `src/multimod/experiments.py`.

## 0. Disclosure of exploratory runs made before this file

To debug the code and choose the λ grid, I ran these smoke tests (scratch scripts, not
kept):
- DCL model counts on N = 2^n + 2^{n/2} + 1 for n ∈ {40, 48, 56, 64, 80, 96}, λ ∈ {3, ln √N + 2},
  Hmax ∈ {12, 24, 40}. Seen: R1 cost 316 / 783 / 2079 / 5429 / 35841 / 239092 against HH
  267 / 797 / 2297 / 6891 / 60731 (n = 40…80). That is a slope ≈ 0.17 for R1 and ≈ 0.197 for R2,
  against HH ≈ 0.196.
- Exhaustive coverage of all integer candidates for the same special N at n = 36–48. R1 left
  0.1–4% uncovered (below the Poisson prediction e^{−λ}); R2 left 0 uncovered.
- Operational DCL at λ = 3 on `balanced_semiprime(bits, seed)` for bits ∈ {32, 40, 48} and
  seeds 0–5: 18/18 factored. **These seeds overlap E-O below**, so E-O seeds 0–5 at 40 and 48
  bits are not blind.

The E-M counts are deterministic functions of N, so they are not blind either. They are
pre-registered for the record and for the decision rule. The blind parts are E-C on random
semiprimes, E-S (2000/2000/1000/300 random instances), and E-O at 56 and 64 bits (and seeds
6–15 at 40 and 48 bits).

## 1. Question

Cycloharvey's lower bound N^{1/5−o(1)} covers one guess modulus. In the enlarged model,
which has arbitrarily many moduli, one base, a shared or separate baby list, and identity
coverage, does a several-moduli scheme reach a model exponent below 1/5 at matched success?

## 2. Scheme under test (frozen): divisor-cover Lehman (DCL)

- **Targets.** A target is τ = (v = (A,B), M, r), where:
  - gcd(A,B) = 1, max(A,B) ≤ Hmax = 12;
  - M is a prime with M ∤ 2ABN;
  - r is a root of B r² ≡ AN (mod M). No residue is guessed; the roots are computed from N.
- **Lock.** If p ≡ r (mod M), then Bp + Aq ≡ c := Br + A·N·r^{−1} (mod M²) (HH Prop. 3.3).
  So g = (Bp + Aq − c)/M² lies in the whole-range window [g_lo, g_lo + W). The window is
  computed exactly from F = Bp + AN/p over p ∈ [√(N/2), √N].
- **Coverage.** Target τ covers exactly the integers p ≡ r (mod M) in the range.
- **Design.** Choose (v, M) pairs greedily by efficiency ρ/(M·⌈W/m⌉), where ρ is the number
  of roots, until λ := Σ ρ/M reaches the regime value. The baby length m is optimised on a
  2^{1/4} grid (`dcl.best_design`). Everything is computed from N and public parameters.
- **Regimes.**
  - **R1:** λ = 3. Poisson heuristic: success ≈ 1 − e^{−3} ≈ 0.95.
  - **R2:** λ = ln √N + 2. Heuristically full coverage, the "deterministic" regime.
- **Model cost.** m + Σ_τ ⌈W_τ/m⌉, the same counting as cycloharvey's calibrated model.
- **Operational.** One shared baby list γ^i with γ = β^L and L = ∏ M² over the moduli used.
  The handle of τ is (β^{L/M²})^{B+AN−c}·γ^{−g_lo}. All β^{L/M²} come from an exponent
  product tree. Collision detection is cycloharvey's Alg. 4.1 (imported unchanged), and β is
  random from the seed.

## 3. Experiments (frozen; `src/multimod/experiments.py`, config dumped to `experiments/configs/multimod_experiments.json`)

- **E-M (model counts).** N = 2^n + 2^{n/2} + 1, n ∈ {40, 48, 56, 64, 72, 80, 88, 96}.
  - DCL cost for R1 and R2, plus calibrated HH (`best_hh`) and Harvey (`best_harvey`) counts.
  - Slopes are least-squares fits of log2 cost against n over n ∈ {56, …, 96}.
- **E-C (exhaustive coverage).** n ∈ {32, 36, 40, 44, 48}, seeds 0–3 of `balanced_semiprime`,
  plus the special N.
  - Count every integer candidate p ∈ [√(N/2), √N] covered by no target, for R1 and R2.
  - Exact window verification (`verify_windows`).
  - Diagnostic: λ* = the smallest grid λ (step 0.5) that gives zero uncovered candidates.
- **E-S (sampled success, harness).** Random balanced semiprimes at n = 56, 64, 80, 96 (2000,
  2000, 1000, 300 instances).
  - Records whether the true p lies in some target class, for R1 and R2.
  - m is fixed per (n, regime) from the E-M special N.
- **E-O (operational).** n ∈ {40, 48, 56, 64}, seeds 0–15.
  - Runs DCL R1 and R2 with full lists (early_exit = False), and calibrated Harvey
    (`harvey_search`, parameters from `best_harvey`) on the same N.
  - Records success, class (STRUCTURAL / PM1 / ACCIDENTAL), implementation list size
    against model count, multiplications and wall time.
  - PM1 = a factor found from a gcd on the order of γ, i.e. Pollard p−1 by L. It is reported
    separately and is not the DCL mechanism.

## 4. Decision rules (frozen)

- **SUPPORTED: "sub-1/5 in the calibrated model at fixed success."** All of the following
  must hold:
  - (a) the E-M R1 slope is ≤ 0.185;
  - (b) the HH slope minus the R1 slope is ≥ 0.015;
  - (c) E-S R1 success is ≥ 0.93 at every size;
  - (d) E-O R1 implementation list == model count in every run that builds a list, and
    STRUCTURAL+PM1 successes are ≥ 14/16 at every size.
- **REFUTED for this model at these sizes:** the R1 slope is ≥ 0.19, or E-S R1 success is
  < 0.90 at any size.
- **Anything else is INCONCLUSIVE.**
- **Deterministic regime R2.** Prediction: 0 uncovered candidates in every E-C instance, and
  E-S R2 success = 1.
  - Any uncovered candidate or failure is reported, with its count.
  - The R2 slope is reported. Prediction: 0.18–0.20 at these sizes, because λ grows like
    log N and the weight heights grow with it. No sub-1/5 claim is attached to R2 at these
    sizes; its asymptotic exponent is a conjecture (DCH, §5).
- **Discrepancy rule.** An implementation list ≠ model count, or a window violation, is a
  bug. It is reported, not silently fixed.
- **Fairness note (stated in advance).** HH and Harvey are deterministic with 100% success.
  - R1 is matched to a 95% success target, not to 100%. At matched 100% (R2) I predict no
    finite-size advantage at ≤ 96 bits.
  - The comparison is in list elements, as in cycloharvey. Multiplications and wall time
    are reported but are not decision criteria.

## 5. Theory claims to be audited (not assumed)

1. **Thm A.** In the enlarged model, cost ≥ c·N^{1/6}. Linear (e = 1) targets give
   per-target coverage ≤ C(PW)^{1/2}; Cauchy–Schwarz then applies.
2. **Thm A′.** Targets of degree e ≥ 2 cannot go below N^{1/5}, whatever the moduli.
3. **Thm B.** With J distinct moduli, cost ≥ c·N^{1/5}·(J·M^{o(1)})^{−1/5}.
4. **Thm C.** CRT composes the class lattices multiplicatively: index(M₁M₂) =
   index(M₁)·index(M₂), and φ is multiplicative too. All-classes schemes stay exactly
   guess-neutral; only selective covering across moduli escapes.
5. **Thm D.** For each ε, a fixed DCL design of model cost N^{1/6+O(ε)} covers all but an
   O(log N·N^{−ε}) fraction of balanced semiprimes. Proof: the arithmetic large sieve in p
   for each fixed q. This needs separately a base γ of adequate order.
6. **DCH (conjecture).** Every candidate p is covered at cost N^{1/6+o(1)}. Under DCH the
   enlarged-model exponent is exactly 1/6.
7. **Multi-base.** Independent bases reduce to the 1-D bounds coordinate-wise.
