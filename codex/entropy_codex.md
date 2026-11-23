(Read /Volumes/git/factoring-connection/codex/entropy_common.md.) Role: BUILD AND RUN.
1. Test grammar (start broad; Gemini will extend): compositions of cheap N-only operations — x^N, x^{(N−1)/2}, x^{N±1}
   mod N for small x; Lucas V_N(P,Q) mod N; depth-2 nested gates (x^N + c)^N mod N; floor/remainder gates ⌊N/k⌋,
   N mod k, isqrt(N), frac(√N), frac(N·α) for fixed irrationals; first CF partial quotients of √N and of N/2^k; Jacobi
   symbols of all the above outputs; bits/low bits/parities of all outputs; small sums/products/XORs of these. Enumerate
   ≥ 10^4 distinct tests (bit-valued) by grammar depth.
2. Data: for each size n ∈ {24, 32, 40, 48, 56, 64, 80, 96}, ≥ 2·10^5 fresh random balanced semiprimes (more at small n
   if cheap); split discovery / validation / confirmation sets.
3. Statistic: advantage = Pr[test = target] − ½ (and mutual information) per (test, target, n), with exact CIs;
   multiple-testing control over (tests × targets) on the discovery set; only survivors are re-measured on validation
   and confirmation sets at every n.
4. Scaling: for survivors, fit advantage(n) against c·n^{−k} vs c·2^{−βn}; a "lead" = advantage on confirmation sets
   significantly > 0 at the two largest sizes AND consistent with polynomial decay (pre-register exact criteria).
5. Controls: (a) planted weak leak (a synthetic test agreeing with the target with advantage 1/n² and 1/n) to show
   detection power; (b) low-bit sanity targets (must be null); (c) shuffle controls; (d) accidental-split flagging
   (report with/without).
Pre-register before running. ≤ 3 hours compute. Files in cwd; final summary ≤ 1500 words.
