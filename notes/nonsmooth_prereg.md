# nonsmooth: pre-registration (phase 16, route 3)

Agent `nonsmooth`, 2026-10-20. Written and committed BEFORE E1/E2 were run. Code at the commit of
this file: `src/nonsmooth/{equiv,forms,cheng,census,run_experiments}.py`, tests
`tests/test_nonsmooth.py` (46 pass). The tests ran small instances of the same code as smoke tests
only (96-bit Cheng instances, 50-prime census sample); no E1/E2 data were produced.

## Why these experiments

The route-3 analysis (to be written in `notes/nonsmooth_results.md`) reduces every
"order-annihilation" method (compute [m]x in G(Z/N) for an algebraic group G and a computable m)
to a dichotomy: the prime-to-p part of exp G(F_p) must divide m, so either it is smooth (Pollard
p−1, Williams p+1, Bach–Shallit, ECM) or it is trivial, which for commutative G with nontrivial
torus/abelian part forces an *anomalous* abelian variety (#A(F_p) a power of p). The one known
algorithmic instance of the second branch is Cheng's 4p−1 method (CM curve with trace 1). E1 checks
that our operational implementation of that branch works and does nothing on random inputs. E2
measures how thin the branch is: its coverage of random primes, which by the heuristic below
decays like p^{−(1−θ)/2} when the CM discriminant is restricted to D ≤ p^θ.

## E1: operational CM/anomalous factorer (`cheng.cheng_factor`, sees N only)

- Arms. *special*: N = p·q, p = (D s² + 1)/4 prime with D ∈ {3, 11, 19, 43, 67, 163}
  (class number one, D ≡ 3 mod 8), q a random prime; n ∈ {128, 256, 512}; 10 seeds per (n, D)
  = 180 instances. *control*: 60 random balanced semiprimes per n (180 instances).
- Config: tries_per_D = 24 (64 for D = 3), all six D tried in order; the algorithm is not told D.
- Classification (rule 3): STRUCTURAL iff the factor is the special prime, found via the final
  Z-coordinate gcd at the correct D; otherwise ACCIDENTAL. Control successes are STRUCTURAL only if
  p or q happens to be special.
- **Predictions.**
  - P1: special arm success 180/180, all STRUCTURAL. Per-try success ≈ 1/2 for D > 3 (the
    quadratic twist is chosen by the Legendre symbol of x³+ax+b mod p) and ≈ 1/6 for D = 3
    (sextic twist class 1/3 × twist 1/2). So for each D ≠ 3 the number of tries at the right D is
    geometric with mean ≈ 2. Failure probability per instance ≤ (5/6)^64 + 2^{−24} ≈ 10^{−5}.
  - P2: control arm 0/180 successes. Null for an accidental hit: ≤ 184 tries × 2 × 2^{−n/2+1},
    i.e. < 10^{−15} at n = 128; null for a structural hit (a random 64-bit prime being special
    for one of these D): ≈ 6 · 2^{−32}.
  - P3: cost per try is O(log N) ring operations: mean `mul` per try ∈ [20 n, 30 n] (the
    Montgomery ladder does 24 multiplications per bit).
- **Decision rule.** E1 "passes" iff P1 ≥ 178/180 and P2 = 0/180. Anything else is reported as a
  bug or as a refutation of the claimed mechanism.

## E2: coverage census of the anomalous-CM branch (`census`, diagnostic)

- For b ∈ {24, 32, 40, 48}: 20 000 random b-bit primes p (seeded); D = squarefree part of 4p − 1
  (sympy factorint; diagnostic). Record D, log D / log p, and F_b(θ) = #{D ≤ p^θ}/20 000 for
  θ ∈ {0.25, 0.5, 0.75, 0.9} (exact integer comparison D^den ≤ p^num).
- **Heuristic prediction (stated now).** #{p ∈ [X, 2X] prime : D ≤ Y} ≈ Σ_{D ≤ Y} c_D √(X/D)/log X
  ≍ √(XY)/log X, so F(θ) ≈ c_θ · p^{−(1−θ)/2}.
  - Q1: every D ≡ 3 (mod 8) (proved; a deterministic check).
  - Q2: median of log D / log p ≥ 0.85 at every b (D is typically ≍ p).
  - Q3: slope of log₂ F(θ) against b: θ = 0.75 in [−0.19, −0.06] (predicted −0.125); θ = 0.9 in
    [−0.10, 0.0] (predicted −0.05). θ = 0.5 (predicted −0.25) and θ = 0.25 are reported only, as
    counts will be small.
  - Q4: the prefactor c_{0.75} = F(0.75) · p^{(1−0.75)/2} lies in [0.02, 1.0] at every b.
- **Decision rule.** If Q1–Q4 hold, the coverage scaling √(Y/p) is "experimentally supported" and
  is used in the results note to price the branch: covering a fraction F of primes needs D up to
  Y ≈ F² p, at total cost ≥ Σ_{D ≤ Y} h(−D) ≍ Y^{3/2} ring operations (class polynomials of degree
  h(−D) ≍ D^{1/2±o(1)}, Siegel), i.e. ≥ F³ p^{3/2}. If Q3 fails, the note reports the measured
  slope instead and does not use the heuristic.

## Commands

```
.venv/bin/python -m src.nonsmooth.run_experiments e1
.venv/bin/python -m src.nonsmooth.run_experiments e2
```

Memory: both runs are < 200 MB RSS; run one at a time after `memory_pressure`.
