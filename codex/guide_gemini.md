(Read /Volumes/git/factoring-connection/codex/guide_common.md.) Role: DESIGNER + ADVERSARY (no heavy compute).
Codex is building the pipeline in parallel. Write an independent specification: (1) the exact target, baseline and
the full list of candidate global features worth testing (be creative: what structure COULD correlate with smoothness
beyond size and local divisibility — e.g. arithmetic of a,b, continued-fraction structure of a/b vs roots, Galois/
conjugate relations, multiplicative relations among candidates, "smoothness inheritance" through products/descents,
history-dependent effects); (2) a rigorous leakage checklist (every way a feature can smuggle divisibility by primes in
(y,B] or size information; how to residualize; correct train/test splitting by N and by region; multiple-testing
control); (3) the planted-signal positive control and what effect size would matter for an exponent change (relate to
u log u: a gain that stays O(1) in log-probability is a constant factor; an exponent change needs a log-likelihood
ratio growing like Ω(u log u)·(fraction)); (4) the theoretical prior: known results on correlations of smoothness
(e.g. smooth values of polynomials at consecutive/related arguments, Hildebrand, Granville–Soundararajan,
Elliott–Halberstam-type independence, Balog–Wooley smooth consecutive integers) and which features they predict to be
null vs non-null. Read-only commands plus python3 in your cwd. Final answer ≤ 1600 words, plain text. Write the spec
also to guide_spec.md in your cwd.
