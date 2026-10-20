# cadobench pre-registration

Agent `cadobench`, phase 15, 2026-10-20. Written and committed before any experiment below is run.
The baseline factorizations (c60 smoke, c70, c80, and c90 in progress) were started before this file.
They are plain measurements with no prediction attached, and their numbers are not used to set
any threshold below. (Their parameters, e.g. lpb, were read in order to design the tests.)

## Scope change (orchestrator, 2026-10-20) — task 4 WITHDRAWN, UNRUN

The project wants asymptotic levers, not constant factors. Every constant-factor experiment in
the original task 4 is **withdrawn and was never run**:
- (a) cofactorization chains and batch smoothness as a speed tweak;
- (b) lpb/mfb/lambda tuning;
- (c) NEON/SME/Metal offload;
- (d) bwc tuning.

None was pre-registered before the withdrawal. The batch-smoothness *test* still appears below, but
only as an exact instrument for counting smooth points, not as a speed-up claim.

## Instruments

- `src/cadobench/run_cado.py`: full factorizations; per-phase wall time from CADO's timestamped log.
- `src/cadobench/memwatch.py`: peak RSS and a kill-switch.
  - Kills a job at tree RSS > 8 GB, available memory < 2 GB, or memory pressure ≥ warn for 5 s.
  - Runs one heavy job at a time.
- `src/cadobench/lasbench.py`: las on a fixed special-q range; funnel, `-T` breakdown, yields.
- `src/cadobench/smoothcount.py` (to be written): exact smoothness of the norms at sampled or
  exhaustive (i, j) points in real special-q lattices.
  - Method: Bernstein product/remainder tree over the primes ≤ B, then x | (P mod x)^{2^e}.
  - Exact integer arithmetic (gmpy2). Floating point is used only for ρ and for the fits.
- Test numbers: `balanced_semiprime(bits, 0)` for c60, c70, c80, c90, c100 (c110 only if c100's
  peak RSS < 6 GB). Polynomials and parameters are CADO's own, from the baseline runs.

Definitions used throughout:
- For a pair (a, b) in the q-lattice region, the norms are
  - n0 = |G(a, b)| (rational side);
  - n1 = |F(a, b)|/q (algebraic side; special-q on side 1, as CADO does here).
- "Point" means a lattice point (i, j) of the sieve region with gcd(a, b) = 1.
- **S_B** (pure smoothness): both norms are 2^{lpb_i}-smooth.
- **S_CADO** (CADO's definition): S_B, and also, on each side, the part of the norm above lim_i is
  ≤ 2^{mfb_i}.
- ρ is Dickman's function.
- u_i = log n_i / log 2^{lpb_i}, using the *actual* norm at each sampled point (not a box maximum).
  The model prediction for a point is ρ(u0)·ρ(u1), averaged over the points.
- α-corrected: u_i' = (log n_i + α_i)/log B_i, with Murphy's α of each polynomial
  (α = 0 for a random integer).

## Experiments and predictions

### V1. Single-side smoothness vs Dickman ρ (all sizes c60–c100)

- Sample 2·10^5 uniformly random coprime points per size, spread over 20 special-q in CADO's
  q-range. The seeds are fixed.
- Measure the fraction of points with n_i 2^{lpb_i}-smooth, separately for i = 0, 1.
- Compare with mean ρ(u_i), and with the α-corrected mean ρ(u_i').
- **Prediction V1a.** For both sides and every size, measured / α-corrected-ρ lies in **[0.8, 3.0]**.
  - Reason: Ψ(x, y)/x exceeds ρ(u) by a second-order factor of order 1 + (1−γ)·log(u+1)/log y
    at these sizes.
- **Prediction V1b.** Uncorrected ρ is further from the measurement than α-corrected ρ on the
  algebraic side at ≥ 4 of the 5 sizes. (CADO's polynomials are root-optimised, so α1 < 0.)

### V2. Independence of the two sides (heuristic H: "joint events multiply")

- c60, c70, c80: exhaustively enumerate windows of the sieve region of the first special-q of 5
  fixed special-q.
  - The window is the full region for c60; a centred 2^20-point window otherwise.
- Count S_B points exactly, and compare with Σ_points 1[n0 smooth]·Pr-estimate, i.e.
  joint / (product of the marginal fractions in the same window).
- **Prediction V2.** joint / (marginal product) lies in **[0.67, 1.5]** at each of the three sizes.
  - This holds once the counts are large enough (≥ 30 joint hits). Below that we report the interval
    only, and no verdict.

### V3. las recall against exact truth

- On the same windows as V2, compare the S_CADO set with the relations las outputs for that special-q
  (production parameters, restricted to the window).
- **Prediction V3.** las finds **≥ 85%** of S_CADO points at c60–c80.
  - Its sieve threshold (lambda) and finite ECM effort lose some; we predict they lose < 15%.

### V4. Sieve work vs box area vs output (c90; c100 if time allows)

- Fix 200 consecutive special-q starting at the CADO qmin of the size.
- Run las single-threaded-equivalent timing (`-t 10`, CPU seconds from las) at I = I_CADO − 1,
  I_CADO, I_CADO + 1. Three repeats each.
- Fit CPU = c_area·(area) + c_out·(survivors entering cofactorization) + c_q·(#q) by least squares.
- Report the per-point share c_area·area / CPU at I_CADO.
- **Prediction V4a.** las CPU time per special-q scales with area with exponent **in [0.85, 1.25]**.
- **Prediction V4b.** At I_CADO the per-point share (norm init + sieving + bucket fill, from `-T`) is
  **≥ 60%** of las CPU, and the cofactorization share (`-T` categories 5, 6, 261) is **≤ 30%**.
  - If both hold, the box term E², not the output term, dominates at real sizes. That is the regime
    where an output-sensitive enumerator would pay off.
- **Prediction V4c.** Relations per special-q grow more slowly than area: yield exponent **< 0.9**.
  - Reason: larger regions have larger norms.

### V5. Relations, matrix and linear-algebra share (from the baseline runs c70–c100, c110 if run)

- **Prediction V5a.** The linear-algebra share (bwc + characters) of total wall time is **< 25%** at
  every size c80–c100.
  - At these sizes CADO is far from the asymptotic "sieve = LA" balance of heurnfs §1.
- **Prediction V5b.** The final matrix dimension / unique relations lies in **[0.05, 0.5]**.
- **Prediction V5c.** LA wall time grows faster than sieve wall time from c80 to c100. The ratio of
  the two growth factors is > 1.
  - Reason: the share is moving toward the asymptotic balance.

### V6. Large-prime incidence (assumption T of heurnfs §3.1)

- Take the relations after singleton removal (purge) at c80–c100. Over the primes > lim_i on each
  side, compute the histogram of the number of relations each prime occurs in.
- **Prediction V6.** The fraction of large primes that occur in exactly 2 relations is **< 50%**,
  and the mean number of occurrences is **≥ 3**.
  - So assumption T ("each large prime in exactly two relations") is far from what a real special-q
    collection looks like. heurnfs says that sharing in special-q collections is not covered by T. This
    measures how far the real data sit from it.

### V7. Empirical L[1/3] constant

- Fit log(total wall seconds) = log A + c·x^{1/3}(log x)^{2/3}, with x = ln N, over c70–c100
  (c110 if run). Do the same for sieve-only wall time.
  - This is an ordinary least-squares fit with 2 parameters. With 4–5 points it is descriptive only.
- **Prediction V7.** The fitted c is **in [1.0, 1.9]**, below 1.923 because of o(1) terms. It is
  reported with the caveat that the fit absorbs every o(1).

## Coordination with smoothenum (task 4 of the new scope)

- Export under `experiments/raw/cadobench_testbed_c{D}/`; the format is documented in
  notes/cadobench_results.md:
  - poly and parameters;
  - special-q lattice bases for a fixed q-range;
  - las relations;
  - the exact S_B / S_CADO sets on the V2 windows.
- smoothenum's prototypes are compared with las on identical inputs (same (f, g), same q-lattices
  or the same (a, b) box, same B and mfb). Metrics: the set of pairs found, recall against exact
  truth, and CPU seconds per pair.
- No prediction is pre-registered for that comparison until smoothenum hands over a prototype.
  An addendum will be committed first.

## Pass/fail rule

Each prediction is reported as **confirmed** or **refuted**, with the measured value. Deviations from
the protocol (sizes skipped, sample sizes changed) are disclosed in the results note.
