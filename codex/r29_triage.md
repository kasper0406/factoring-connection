You are an expert in analytic number theory and computational number theory. Work read-only on /Volumes/git/factoring-connection; your cwd is a writable scratch dir. Use web search where useful.

Context: this repo is an AI-assisted research project on the asymptotic complexity of integer factorization. Read README.md first, especially "Top findings", "Key ideas" and "Open problems, ranked". Then read BARRIERS.md, notes/native_results.md and notes/connections_results.md.

The candidate leads are:
(A) Write out Theorem AB (notes/c12alg_results.md), the algebraic two-point blocks under design DW. This would give rigorous L[α] for every α > 3/8. AB₂ would give L[3/8].
(B) NFS fruitfulness: Lemma M / Conjecture F (notes/lemmaM_results.md, notes/nfsfruit_results.md), the batch form of Lee–Venkatesan Conj 7.1. This would give rigorous L[1/3, 1.923].
(C) A single coefficient a_f(N) of a fixed non-CM cusp form at composite N (notes/nativemod_results.md). Does it factor N? Does it have a native algorithm?
(D) Deterministic worst case below N^{1/5}, via the coverage38 positivity count and Bourgain–Garaev-type bounds (notes/coverage38_results.md, notes/conn_alg.md).
(E) Three-level descent at 4/11 (notes/pc3_results.md §7.5).
(F) Constant-level items: sub-quadratic sparse GF(2) linear algebra, matroid elimination, percolation-guided relation selection.
(G) Anything else you consider more promising after reading the repo.

Answer concisely:
1. Rank the leads by (probability of a real result within a focused effort) × (significance). Give a probability and a one-line justification for each.
2. For the top 3, give a concrete attack plan: the first lemma to prove, the key published tools with exact references, and the most likely failure point.
3. For (B) specifically: is there a new angle the project's earlier attempts missed, given the later tools in the repo (Dixon descent, Theorem DD, decoupling Prop DC, hybrid BV Lemma 8.5H, design DW)? For example, can a hybrid design get fruitfulness at L[1/3] by descending only part of each relation, or by randomizing the square-root step?
4. For (C): is there a known reduction from a single a_f(N) to factoring, or a known obstruction?
