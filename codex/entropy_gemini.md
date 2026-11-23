(Read /Volumes/git/factoring-connection/codex/entropy_common.md.) Role: DESIGNER + ADVERSARY (no heavy compute).
Codex is building the pipeline in parallel. Write an independent spec (also to entropy_spec.md in your cwd):
(1) the richest sensible test grammar beyond codex's starting list — what N-only computations COULD carry information
    on bits of p without splitting N (think: the three known leak families group-order/smoothness, archimedean,
    genus; PRIM nested gates (notes/prim_results.md), class-group / quadratic-form computations of disc −N or 4N at
    small cost, continued-fraction / infrastructure walks, modular-form / cusp coefficient proxies, lattice
    reductions on small N-dependent lattices, elliptic-curve point arithmetic mod N with small multipliers);
(2) targets, especially size-free versions of top bits of p and GL parities;
(3) the statistics protocol: multiple testing over a huge program space, how to separate genuine polynomially-decaying
    advantage from small-n artifacts, sample sizes needed to detect advantage 1/n² at n = 96, how to fit and test
    polynomial vs exponential decay;
(4) a leakage / artifact checklist (accidental splits, size leaks, CRT-correlated targets, prime-generation biases,
    pseudo-random-generator artifacts, sorting p<q);
(5) the theoretical prior per test family (which are provably null: Theorem PT, AC0-ε, PT-J Jacobi; which are open:
    PRIM depth ≥ 2, CF statistics H-GEN).
Read-only commands plus python3 in cwd. Final answer ≤ 1600 words, plain text.
