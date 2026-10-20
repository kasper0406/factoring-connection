You are an expert in algebraic complexity, sparse linear algebra over finite fields, and the number field sieve. Work read-only on /Volumes/git/factoring-connection; your cwd is a writable scratch dir. Use web search.

Context: notes/heurnfs_results.md (the heuristic cost model). With linear algebra cost L^{ω·β}, an effective exponent ω ≤ ~1.67 instead of 2 would lower MNFS from L[1/3, 1.902] to ~1.86. We want to know whether there is NOVEL, worthwhile research here, before committing a long agent thread.

Answer concisely and concretely:

1. State of the art. Give the best known algorithms, with exact bounds and references, for:
   - finding a nonzero kernel vector (one dependency) of a sparse n×(n+k) matrix over GF(2) with n·polylog nonzeros;
   - solving a general sparse system or computing rank.
   Cover Wiedemann, block Wiedemann/Lanczos, Eberly et al., Peng–Vempala 2021, Casacuberta–Kyng 2022, and any 2022–2026 work. Is anything below n² known for any natural sparse class over finite fields?

2. Known lower bounds and barriers. Is "below n² for sparse GF(2) kernel finding" a recognised open problem? What would it imply? Are there known reasons (black-box Krylov lower bounds, etc.) that Krylov methods need ~n matrix-vector products?

3. NFS-specific structure. NFS relation matrices are divisibility-incidence matrices: row = lattice point (a,b), column = prime ideal dividing its norm. They have:
   - heavily skewed column weights;
   - many large-prime columns of weight 2–3;
   - an arithmetic "row generation" structure.
   Has anyone exploited this beyond standard filtering, structured Gaussian elimination, and large-prime cycle finding? Propose up to 3 concrete NFS-specific ideas that might find ONE dependency in sub-n² time. Examples: a dependency found during sieving via cycle structure; using the lattice structure of rows; exploiting that we can choose WHICH rows to generate (adaptive relation collection aimed at closing a dependency); trading more sieving for a structurally easier matrix. For each, give a quick feasibility check and say whether it is known.

4. Verdict. Is there a plausible novel result that experts would find interesting, at least as a well-defined NFS-specific question with a first lemma? Give a probability and a recommendation: run a thread, or don't.
