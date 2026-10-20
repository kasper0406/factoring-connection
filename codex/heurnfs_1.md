# heurnfs_1 — blind round: the HEURISTIC L[1/3] frontier for general-N factoring

You are a research number theorist. Read-only access to the repo at /Volumes/git/factoring-connection
(you may read any file; scratch work goes in your cwd). Web search is enabled; use it for literature.

## Task
Find any HEURISTIC improvement to factoring a general integer N (no special form), i.e. either
(a) an L_N[1/3, c] algorithm with c < 1.902 (Coppersmith MNFS, (92+26√13)^{1/3}/3), or
(b) a heuristic exponent below 1/3.
Honesty matters more than a positive verdict. Label every statement KNOWN (with citation),
PROVED (give the proof), or HEURISTIC (state the exact heuristic assumption: e.g. "norms behave like
uniform integers of the same size for smoothness", Canfield–Erdős–Pomerance u^{-u(1+o(1))}).

## Required output sections
1. **State of the art, precisely.** GNFS (64/9)^{1/3}; Coppersmith MNFS; factorization factory
   (Coppersmith; amortized 1.639, precomputation 2.007); Bernstein–Lange batch NFS (1.704 per N
   amortized? check); exTNFS / MexTNFS / Kim–Barbulescu / Sarkar–Singh constants (DL only: say why they
   do not apply to factoring); Lee–Venkatesan 2018 rigorous NFS; any 2018–2026 work claiming a
   factoring constant below 1.902 (search arXiv/ePrint; report honestly if none).
   Write the classical GNFS optimisation explicitly: degree d = δ (log N/loglog N)^{1/3},
   coefficient bound, sieve box E = L[1/3,e], factor base B = L[1/3,β], norms
   rational ≈ E·N^{1/d}, algebraic ≈ E^d·N^{1/d}, smoothness u^{-u}, linear algebra B^{2}; the
   resulting system and its optimum; same for MNFS (V algebraic sides with factor base B^{κ}).
2. **Devices to evaluate** (the project's — read the notes for definitions; summarise each in 3 lines
   then compute its effect on the L[1/3] constant in the classical heuristic model):
   - Dixon descent of a random square (notes/lemmaM_results.md): relations from descending random targets
     rather than sieving a box. Does "descent-style relation generation" change the cost balance?
   - "Split the rational side into k primes then descend" (Theorem K; grep notes/ for "Theorem K").
   - Recursion (notes/recurse_results.md): recursive NFS for cofactorization of large primes;
     confirm or refute that large-prime variations only give 1+o(1) in the constant.
   - The N·w multiplier (notes/ldbound_results.md): use kN or many multipliers — any constant gain?
   - 2D-random rational side (notes/pcest_results.md) and multiple rational sides.
   - Batch smoothness (notes/lv6_results.md) — Bernstein batch ECM replacement; any constant gain?
   - Norm-height vs search-entropy separation (codex/r20_wild.out §5) and small-characteristic
     elimination with Witt-vector carries (§6).
   For each: the modified optimisation system, its optimum constant (numerically, to 4 digits),
   and if no gain, the precise reason (which term in the balance is unaffected).
3. **Structurally new attempts at exponent < 1/3.** For each attempt give the precise heuristic
   barrier. Discuss: why 1/3 arises (degree × coefficient-size balance: norm ≥ N^{2/(d+1)}-ish from
   the resultant/Minkowski bound, and the smoothness exponent), the best possible norm product
   for a pair of polynomials with common root mod N (lattice argument: |Res(f,g)| ≥ N, so
   ||f||∞·||g||∞ heights are constrained), lattice/Coppersmith-style constraints, descent-based
   (as in quasi-polynomial small-char DL — why the Frobenius/field-of-definition structure is missing
   over Z), and function-field analogues.
4. **Best candidate**: if you find any variant with a constant < 1.902 or exponent < 1/3, give the full
   derivation and then attack it yourself adversarially (where does the heuristic fail? is it already
   in the literature?). If none survives, say so plainly and rank the near misses.

Be numerically careful: write the optimisation as an explicit function and minimise it (you may run
python in your scratch dir). Keep it under ~4000 words.
