# Blind round: unexpected connections between ten mathematical fields and integer factorization

You are a research mathematician and complexity theorist. The repository (read-only for you) is a research project on
asymptotically faster CLASSICAL integer factorization. Context you may read: README.md, BARRIERS.md,
notes/leakmap_results.md §0, notes/suffmap_results.md §0, notes/factorbits_results.md §3, notes/newfamily_results.md §0,
notes/geometry_results.md §0. Briefly: every known way to make the primes p, q of N = pq "behave differently" is
(M) a group order / smoothness, (A) an archimedean approximation (exponential cost), or (Q) a quadratic genus symbol.
Records to beat: heuristic L[1/3, 1.902] (MNFS), rigorous L[2/5] (this project, internal), deterministic N^{1/5}.

Task: for EACH of the fields below, propose 3–6 CONCRETE, unexpected candidate connections to factoring N. Go beyond
the seed ideas. For each candidate:
1. What object is built from N alone (no factors), and what poly(log N)-computable (or subexponential) feature
   would carry factor information, OR what algorithmic primitive would replace sieving / sparse linear algebra.
2. Prior literature (cite precisely; use web search) and whether it is already used in factoring.
3. Your verdict: RESTATEMENT (needs a factor to write down) / AGGREGATE (depends only on N mod m) /
   FACTORING-EQUIVALENT (or DL-equivalent) / REDUCES-TO-M, -A, -Q or a named barrier / KNOWN (used in factoring) /
   OPEN-WORTH-A-LOOK, with the EXACT step where it fails (needs factors, costs size N, lands in M/A/Q).
4. Flag anything you think is genuinely new or surprising.

Fields (seed ideas in parentheses; do not stop at them):
1. Knot theory / low-dim topology (arithmetic topology Mazur–Morishita, Heegaard Floer d-invariants of lens spaces
   L(N,q), quantum knot invariants at N-th roots of unity, poly-time-computable 3-manifold invariants, e.g. Turaev–Viro
   at fixed level, Alexander polynomials, torus knots T(p,q) vs T(1,N)).
2. Symplectic geometry / Floer theory (Fukaya categories of tori with N-dependent data, HMS for elliptic curves,
   Floer homology ranks depending on divisors, Lagrangian intersection counts).
3. Homotopy theory (stable stems, chromatic picture, Morava K-theory, TMF and level-N modular forms, image of J and
   Bernoulli denominators, equivariant homotopy for Z/N).
4. Category theory / higher categories (categorification of divisibility, Grothendieck groups, categories whose
   decategorification encodes σ(N) or τ(N); fusion categories / modular tensor categories of rank depending on N).
5. Extremal combinatorics / Ramsey / matroids (matroids representable over Z/N, Ramsey-type guaranteed structure in
   sets built from N, sum-free sets, the polynomial method, Paley-type graphs over Z/N).
6. Spectral graph theory / expanders (Cayley graphs of (Z/N)^*, SL_2(Z/N), LPS Ramanujan graphs, eigenvalue
   multiplicities revealing factors, random walks, mixing times, local algorithms that estimate spectra).
7. Tropical geometry (tropicalisation of xy = N, valuations, tropical curves over Q_p, Newton polygons).
8. Discrete geometry / sphere packing / tilings (Viazovska-type modular-form methods, tilings of Z_N:
   Coven–Meyerowitz, Hajós / de Bruijn / Sands factorizations of cyclic groups, Kasteleyn/dimer counts for regions
   built from N, substitution tilings with N-dependent inflation, aperiodic monotiles). TREAT TILINGS CAREFULLY:
   e.g. is "does Z_N admit a non-periodic / non-Hajós tiling with a given tile" or counting tilings of Z_N by a tile
   computable without factoring, and does it reveal factors?
9. Logic / proof complexity (bounded arithmetic: factoring and the weak pigeonhole principle, Krajíček–Pudlák,
   Jeřábek; proof-complexity lower bounds for factoring tautologies; TFNP classes PPA/PPP/PWPP and factoring;
   what a 2026 refutation of ETH for 3-SAT would or would not change).
10. Cellular automata and computation models (Rule 110, reversible CA computing mod-N dynamics, self-assembly,
    analog / BSS models, unit-cost arithmetic RAMs (Shamir 1979: factoring in O(log N) unit-cost ops) and what
    restriction makes such models honest).

Finish with your ranked top 3 OPEN-WORTH-A-LOOK items (if any survive), each with a precise first lemma or question.
Be honest: most candidates should die, and say where. Keep the whole answer under ~5000 words, as a table per field
then the shortlist.
