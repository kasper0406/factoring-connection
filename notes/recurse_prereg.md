# recurse: pre-registration of E1 (toy recursive shape test)

Agent `recurse`, phase 11, 2026-10-20. Frozen before the E1 run. It is a correctness smoke test only. Toy sizes cannot
test any asymptotic claim.

**Code (frozen with this file):**
- `src/recurse/pipeline_rec.py` (operational): the thmA toy pipeline (`src/thmA/pipeline.py`, imported unchanged), with
  the shape test delegated to an oracle. The oracle's coins come from a private stream that is consumed once per shape
  test. Oracle `pb` is Pollard–Brent, the thmA stand-in for Lenstra–Pomerance. Oracle `rec` is the pipeline itself,
  called recursively on the non-B-smooth cofactor of each split half.
- `src/recurse/harness_e1.py`: the harness. It knows the labels. It uses sympy only to label ω of inner split inputs.
- `src/recurse/analyze_e1.py`: the analysis.
- `experiments/configs/recurse_e1.json`: the configuration.

**Instances (22):**
- semi: 46 bits ×4, 50 bits ×6, 54 bits ×2;
- tri: 52 bits ×4, 56 bits ×2;
- small_x_semi: 50 bits ×2;
- p2q: 48 bits ×2.

Each instance is run with both oracles, using the same seed.

**Calibration disclosure.** I ran seeds 900–904 before freezing: semi 46, 50 and 56 bits, tri 56 bits. They are
excluded from E1. All factorisations were correct, and the outer transcripts were identical. Recursion was exercised
at 50+ bits, to depth 2. A 56-bit semiprime took 821 s with `rec` against 11 s with `pb`, so 56-bit semiprimes are
omitted.

**Predictions:**
- **P1.** Every complete factorisation is correct, for both oracles (22/22 each).
- **P2 (oracle substitution, Lemma OS).** For every instance, the depth-0 transcript is identical for `pb` and `rec`.
  The transcript comprises each outer split's input, outcome, factor, accept/reject counters and dependency outcomes.
  Predicted 22/22.
- **P3.** The recursion is exercised: at least one inner SPLIT call (depth ≥ 1) occurs in at least 80% of the semi/tri
  instances with ≥ 50 bits. Reported but not predicted: the maximum depth.
- **P4.** There are 0 row-identity violations, over all splits at all depths and both oracles.
- **P5 (DD at all depths).** In the `rec` runs, the fruitful rate over all tested dependencies, grouped by the ω of the
  split input, is within |z| ≤ 3 of 1 − 2^{1−ω}. The same check is made separately for inner splits only, when they
  include at least 20 tested dependencies.

Descriptive only: split kinds by depth (STRUCTURAL / ACCIDENTAL), and the rec/pb time ratio. The ratio is expected to
be > 1 at toy sizes, because Pollard–Brent is cheap there. It carries no asymptotic meaning.
