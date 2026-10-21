# Pre-registration: factorbits E1 (Theorem PT on concrete pure generic-ring trees)

Agent `factorbits`, phase 21, 2026-10-21. Committed together with `src/factorbits/{pure,run_experiments}.py` and
`tests/test_factorbits.py` **before** E1 is run at the registered sizes. Runs before freezing: the unit tests (12-bit
primes, exhaustive) and one smoke run of `run_experiments e1 10 12` (10- and 12-bit primes, written to the scratchpad,
not to `experiments/raw/`). Both showed 0 mismatches and 0 bound violations. Nothing else was looked at.

## 1. What is tested

**Model.** A *pure* generic-ring decision tree computes integers c_v by straight-line programs whose constants and
shape do not depend on N. It branches on [c_v ≡ 0 (mod N)] and, being instrumented, reports a split when
1 < gcd(c_v mod N, N) < N. The tree may be adaptive and of any size.

**Partition lemma (claimed proved; proof in notes/factorbits_results.md §3).** For distinct primes x, y the tree run on
N = xy splits N iff the diagonal paths of x and y differ, where the diagonal path of x answers [x | c_v] at each node.
Otherwise it ends in the common leaf.

**Theorem PT (claimed proved).** Let X, Y be iid primes and F(X, Y) any target. Let π_l be the class probabilities,
η = 1 − max π_l and P_split = 1 − Σ π_l². Then P_split ≥ η, and
Pr[tree outputs F, or splits] ≤ max_b Pr[F = b] + P_split + η².

**ε cross bound (claimed proved).**
- Population A: ε = ±1 with probability 1/2 each, X, Y iid uniform on the primes ≡ 1 (resp. 3) mod 4 of the range.
- Population B: X ≡ 1, Y ≡ 3 (mod 4).
- Then Pr_A[tree outputs ε, or splits] ≤ 1/2 + P_split(A)/2 + √(2 P_split(B))/2.

E1 checks the lemma against real instrumented execution mod N, and checks both inequalities exactly. These are
implementation and proof checks, not discoveries.

## 2. Instances and trees

- Prime sets: all primes with b bits and the top two bits set (`pure.primes_top2(b)`), for b ∈ {14, 16, 18}, i.e.
  n = 28, 32, 36. Laws are uniform on these sets.
- Trees are as in `run_experiments.trees_for(b)`. M(B) = ∏_{ℓ≤B} ℓ^{⌊log_ℓ 2^b⌋}.
  - `pm1_nested`: tests 2^{M(16)} − 1, 2^{M(64)} − 1, 2^{M(256)} − 1.
  - `mixed_adaptive`: p−1, p+1 (Lucas, P = 3), ECM (a24 = 7, x0 = 2) and p−1 base 3, branching adaptively.
  - `ecm3`: three N-independent Montgomery curves, M(128).
  - b ≤ 16 only:
    - `lipton_interval`: ∏_{x∈[lo, lo+2^{b−3}]} x;
    - `qr_sniff_full`: ∏_{x ≤ 2^{b−1}} (x²+1);
    - `qr_sniff_partial`: ∏_{x ≤ 2^{b/2+2}} (x²+1).
- Real-execution check:
  - b = 14: all unordered pairs;
  - b = 16, 18: 20 000 seeded random pairs per tree (2 000 for the product trees), with seed string
    `factorbits|e1|{b}|{tree}`.
- Targets for PT: bit 1, bit 2, bit b−3 and bit b−4 of min(X, Y), and X mod 4.

## 3. Predictions and decision rule

- P1: partition mismatches = 0 for every (b, tree).
- P2: PT violations = 0.
- P3: ε-cross violations = 0.
- Any nonzero count means the lemma, the theorem or the code is wrong. The result note must then report it and withdraw
  the affected claim until it is resolved.
- Descriptive only, with no prediction: the leak (success − baseline) against P_split per tree and target.
  `qr_sniff_full` is a designed extreme case. It should decide ε exactly with P_split(A) = 0 and P_split(B) = 1.

## 4. Command

`.venv/bin/python -m src.factorbits.run_experiments e1 14 16 18` writes `experiments/raw/factorbits_e1.jsonl`.
Run it alone after checking `memory_pressure`. Expected RSS is below 1 GB.
