# smoothenum pre-registration (agent `smoothenum`, phase 15, 2026-10-20)

Committed before any pre-registered run. Code: `src/smoothenum/{nfsdata,mech,run_proto}.py`; tests
`tests/test_smoothenum.py` (14 pass). All runs are single-process, stream JSONL, keep one sieve row or one
lattice in memory (target peak RSS < 1 GB; recorded with `/usr/bin/time -l`), and start only when
`memory_pressure` reports normal.

**Disclosure (engineering runs before this file).** One smoke run (not part of the registered runs):
rational and two-sided line sieve of rows b = 1..8, |a| < 2^14, on the instance below. It found 6855
rational relations in 262144 points, and 216 two-sided ones. The sieve recalled 2498 of the 2571 exactly
tested rational relations in rows 1–2 (97%). One q-lattice was enumerated: 2759 points, against 2758 expected.
The analytic L[1/3] programs (`src/smoothenum/target.py`, `barrier.py`) and the codex blind round
`codex/smoothenum_1.out` were done before this file. They contain no experimental data.

## Instance (real CADO data)

- N = `balanced_semiprime(200, 1501)` (61 digits). CADO-NFS (build of agent cadobench, used read-only)
  `polyselect -degree 4 -P 420 -admax 10000 -incr 60 -nq 64 -keep 10 -t 2`, best size-optimized
  polynomial, then `polyselect_ropt -ropteffort 0.1 -I 10 -Bf 111342 -Bg 78682`.
  Result: `experiments/raw/smoothenum/c61.poly` (degree 4, skew 9594, MurphyE 3.40e-7).
- Parameters = CADO `params.c60`: lim0 = 78682, lim1 = 111342, lpb0/1 = 18/19, mfb0/1 = 17/38.
- Relation tests (exact): side i passes iff |norm| = (lim_i-smooth part)·c, with c ≤ 2^mfb_i and every
  prime of c below 2^lpb_i.
- Box: |a| < 2^21, 1 ≤ b ≤ 2^9, so U = 2^31 points. Seed 20261020.

## Runs and frozen predictions

**calib (prediction input, not a test).** P_s(k) is the fraction of uniform random integers in
[2^{k−½}, 2^{k+½}) passing the rational test. It is measured at k = 14, 16, …, 60 with 20000 samples each. The
algebraic test is measured at k = 20, 24, …, 112 with 60000 samples each. Predictions interpolate log P_s
linearly in k, with P_s = 1 below 2^14.

**E1. Box line-sieve baseline.** The run covers 32 rows b = 1 + 16k with random windows of 2^16 a-values, on
the rational side and on both sides. A recall check exactly tests every coprime point of the first 2 rows.
- **P1:** observed rational relations / Σ_coprime points P_s(|G|) ∈ [0.80, 1.25].
- **P1′:** sieve recall ≥ 0.95.

**E2. Smooth-modulus lattice sweep, rational side.** This tests mechanisms (b) and (c) and Proposition LQ.
- *Moduli.* M is squarefree, built from random distinct rational factor-base primes ≥ 2^8, with
  M ∈ [2^t, 2^{t+1}). The buckets t and lattice counts K_t are
  {8:60, 12:60, 16:60, 20:60, 24:60, 28:60, 31:2000, 34:4000, 37:20000, 40:40000}.
- *Enumeration.* The box points of Λ_M are counted exactly. All points are tested if there are ≤ 2000;
  otherwise 2000 uniform samples are tested.
- *Cost.* Each lattice costs n_M + 1, one unit for the setup.
- *Measured quantities.*
  - Y_t = Σ estimated relations / Σ cost.
  - Y_box = rational relations per box point (from E1).
- **P2 (geometry):** Σ n_M / Σ U/M ∈ [0.95, 1.05] for t ≤ 28. For t ≥ 31 the ratio lies in
  [0.90, 1.10] (Poisson-limited).
- **P3 (cofactor heuristic H_LQ):** in every bucket with ≥ 30 observed relations,
  observed / Σ_{tested coprime points} P_s(|G|/M) ∈ [0.75, 1.33].
- **P4 (the barrier shape).** Ŷ_t is computed by `geo_prediction()` from geometry, P_s and the E1 |G|
  sample only (no E2 data).
  - (i) argmax_t Y_t ∈ {28, 31, 34}, i.e. the peak is within 3 bits of log₂U = 31.
  - (ii) max_t Y_t / Y_box ≤ 1.5 · max_t Ŷ_t/Ŷ_box.
  - (iii) Y_37 < Y_34 and Y_40 < Y_37. The observed log₂-slope of Y_t over t ∈ {34, 37, 40} lies within ±0.3
    per bit of the slope of Ŷ_t.

  This is the finite-size form of Proposition LQ: beyond M ≈ area the yield per unit cost falls like U/M,
  and the peak gain is the cofactor-size gain P_s(y/U)/P_s(y). That gain is u^{O(1)}, which is L^{o(1)}
  asymptotically, even when it is large at this size.

**E3. Bottom-up generation (toy).**
- *Setup.* N = `balanced_semiprime(64, 1503)`, base-m with d = 3 (m ≈ 2^16), box |a| ≤ 2^9,
  1 ≤ b ≤ 2^8.
- *Generation.* 200000 random products of primes ≤ 2^10 with s ∈ [2^20, 2^23).
- *Hit rule.* A hit is a box pair with G = ±s.
- **P5:** hit rate / (2A/m) ∈ [0.90, 1.10]. This is a check of the exact injectivity and density statement;
  the c61 value 2A/|Y0| is reported analytically.

**E4. Two-sided special-(q₁, q₂) vs las at equal output.**
- *las run.* `makefb -poly c61.poly -lim 111342 -maxbits 10 -out c61.roots1.gz -t 2`, then
  `las -poly c61.poly -fb1 c61.roots1.gz -lim0 78682 -lim1 111342 -lpb0 18 -lpb1 19 -mfb0 17 -mfb1 38
  -ncurves0 2 -ncurves1 2 -I 10 -q0 61961 -q1 63961 -t 2 -dup -out c61_las.rels`
  under `/usr/bin/time -l`. The special-q is on side 1, as in CADO c60.
- *K.* K = the number of unique las relations.
- *Our mechanism.* `run_proto e4ours K 3600`. It picks a rational q₁ ≥ 2^13 and an algebraic ideal
  q₂ ≥ 2^13, both random factor-base elements, and uses the lattice of index q₁q₂ in our box. Every box
  point is tested on both sides. The run stops at K distinct relations or after 1 h.
- **P6a:** ≥ 99.5% of the unique las relations pass our exact two-sided test.
- **P6b (wall time):** las is ≥ 20× faster per relation than our Python mechanism (2 threads vs 1).
  This is recorded as an implementation fact, not as evidence about asymptotics.
- **P6c (unit cost, descriptive with one prediction).** Our relations per unit cost exceed E1's two-sided
  box yield per point (Y_box_both) by a factor > 1 and < 10^4. This is a constant-size u^{O(1)} gain. We
  register no tighter band, because the algebraic-side P_s carries an unknown α-correction.

**Decision rule.** A failed P1–P5 is reported as a failure of the corresponding heuristic at this size.
For P3/P4, the barrier's finite-size quantitative form is then **refuted at this size**. No prediction
bears on the L[1/3] constant directly. The asymptotic statements are model theorems (results note §3), and
the experiments test only the per-point heuristic they rest on.
