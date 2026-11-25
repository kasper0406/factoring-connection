# dig B: adaptive factor-interval refinement down to width N^{1/4}

Working hypothesis (from notes/inp_results.md): a poly-time factorer could be a cheap rule R(N, I) that, given a
candidate interval I for the smaller factor p (balanced N = pq), scores the two halves, and keeps the true factor
among poly(n) surviving intervals all the way down to width N^{1/4}, where Coppersmith finishes. Dig in. Repo (do NOT
edit; scratch in your cwd): notes/inp_results.md and experiments/inp/dual_interval_pilot.py (the dual-Coppersmith
sibling score: perfect at radius N^{0.12}, chance at N^{0.28}, N^{0.34}), notes/topblock_results.md (Prop T: floor
ladder = Coppersmith known-bits), notes/floorsum_results.md, notes/floorhard_results.md (general floor oracle
NP/#P-hard; top blocks open), notes/analytic_results.md (floor-sum ladder L_j, N^{1/(j+2)} given an oracle),
notes/classE_results.md (CF convergents of a/N), notes/leakmap_results.md (archimedean family A; AC0-ε),
notes/below16_results.md, notes/unbalanced_results.md, notes/avgtime_results.md, notes/codex_r8_copp_bridge.md, and the
"bucket" discussion: deciding "does N have a factor in [a,b]" is the decision version of factoring.
Questions (≤ 2000 words):
I1. Information/error analysis: with branching factor 2 and per-step correctness 1/2 + δ, how large must δ be (and how
    must errors correlate) for a beam of width poly(n) to keep the true path over ~n/4 halvings? What does this imply
    about the advantage any score must have at EACH scale (from width √N down to N^{1/4})? Contrast with the
    exponentially decaying leaks found (N^{−1/4}).
I2. Candidate scores at an interval [a, a+X]: exact or approximate counts of lattice points under xy = N with x in the
    interval (Σ_{d∈I} (⌊N/d⌋ − ⌊(N−1)/d⌋) — exact count is the decision problem; what about smoothed / weighted /
    approximate versions computable fast via Dirichlet-hyperbola / Voronoi-type summation (N^{1/3} algorithms;
    truncated Voronoi with error terms)? at what interval widths X does an approximation of error < 1/2 become cheap?),
    Coppersmith-type lattices with various shifts (primal/dual, multiple centers), continued fractions of N/a², Lehman
    multipliers restricted to the interval, quadratic forms / Fermat-type a² − N near the interval, analytic tests
    (Σ_{d∈I} e(N/d)-type exponential sums). For each: what it computes, its cost as a function of X, its signal-to-
    noise at the critical widths, and the closest repo barrier.
I3. The key quantitative question: is there ANY scale X (between N^{1/4} and N^{1/2}) at which a poly-time-computable
    statistic distinguishes "I contains p" from "I does not" with advantage ≥ 1/poly? (The dual pilot says no at
    N^{0.28}; Coppersmith says yes at ≤ N^{1/4}; Lehman/Hart give structure at widths related to p/q.) Prove a barrier
    for a natural class (e.g. all scores that are functions of the Coppersmith lattice Gram data, or of truncated
    Voronoi sums of given length), or find a statistic.
I4. Run a cheap pre-registered experiment testing your best statistic at widths N^{0.26}, N^{0.30}, N^{0.40} on
    48–96-bit N (parent interval supplied, child hidden), with exact verification and chance baselines.
I5. Verdict and probability.
