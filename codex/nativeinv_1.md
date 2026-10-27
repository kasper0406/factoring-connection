# Blind round: factoring-equivalent invariants of N and their NATIVE computation cost

You are a research mathematician and complexity theorist. Context: a research project on asymptotically faster
classical integer factorization (targets to beat: heuristic NFS L_N[1/3, 1.902]; rigorous L_N[2/5];
deterministic N^{1/5}, Harvey–Hittmeir). Rigour and honesty matter more than a positive result.

## The idea
Many invariants X(N) of a composite N = pq (or general N) satisfy "computing X(N) lets you factor N":
e.g. φ(N), σ(N), #E(Z/N) for a fixed curve, class-number-type data, the number of conjugacy classes of a
group built from Z/N. Such a reduction shows X is at least as hard as factoring. It does NOT show factoring is
hard. Each X is a NEW VIEW of the problem; the question is whether X has a **native algorithm** — one that never
factors N (never computes p, q) — that runs in time polynomial in log N, or at least sub-exponentially
faster than the obvious "factor, then use the CRT/multiplicativity formula" or "enumerate Z/N".
Caveat: if the natural generalisation of X (e.g. X for arbitrary inputs of the same shape: arbitrary groups,
arbitrary graphs, arbitrary complexes) is NP-hard or #P-hard in a relevant way, it should be deprioritised.

## Project barriers you must test each candidate against
- **Theorem PT:** algorithms that use N only as a modulus (ring operations in Z/N, plus gcd) learn nothing
  about the factors except by splitting N (generic-ring / black-box CRT argument).
- **PRIM / N-dependent exponent gates:** every known leak between p and q is (M) a group order / smoothness
  (p−1, #E(F_p), class numbers), (A) an archimedean approximation (costs exponential precision), or (Q) a
  quadratic character (learns only N mod small m). Tests of the form x^{E(N)} mod N with E an integer function
  of N are classified (Monier-type liar counts); nested gates open at depth ≥ 2.
- **Size–entropy law:** any "direct target" whose solution would give exponent < 1/3 is poly-time factoring in
  disguise.

## Task (work independently; do not look at the repo's notes for this part — this is a blind round)
1. Give your OWN list of 25–40 factoring-equivalent (or factoring-hard via randomized reduction) invariants
   X(N), across: arithmetic functions; algebraic K-theory and group cohomology of Z/N; finite group theory and
   representation theory of groups built from Z/N (GL_2, SL_2, Aff, Heisenberg, dihedral...); topology
   (lens spaces, Seifert fibred spaces, torus knots, mapping tori); combinatorics (necklaces, Burnside counts,
   tilings, graph invariants of Cayley/Paley-type graphs on Z/N, circulant matrices); dynamics (cat maps,
   x→x^2 mod N, orbit counts); modular forms and L-functions (coefficients at composite index, class numbers);
   algebraic geometry over Z/N (point counts).
   For each: (a) the precise reduction to factoring (cite or prove in 2 lines; say if randomized / needs
   N = pq / needs multiplicity); (b) the definition-level native cost; (c) the BEST KNOWN native algorithm
   that does not factor N, with cost and a reference (or "none better than definition"); (d) complexity of the
   natural generalisation (NP-hard? #P-hard? in P? in NP∩coNP? TFNP?), with reference; (e) which barrier above
   covers its native algorithms, if any.
2. Single out invariants whose native algorithm is NOT obviously a "sum over Z/N" or "multiplicative function
   evaluated by CRT": e.g. invariants computed by homology of a small complex, by a transfer matrix / automaton
   over the digits of N, by a closed-form generating function, by a determinant, by a recursion in N (like
   continued fractions / Euclid / Dedekind-sum reciprocity). Reciprocity-type computations (Dedekind sums are
   poly-time via reciprocity, Jacobi symbols likewise) are exactly the pattern of a native polylog algorithm
   for an invariant that *looks* like a sum over Z/N. Which factoring-equivalent invariants could plausibly
   have such a reciprocity/recursion? Which provably cannot (e.g. because the invariant is multiplicative and
   a reciprocity algorithm would be a PT/ring algorithm)?
3. Propose a precise **model** in which "every native algorithm for a multiplicative arithmetic function
   f(N) costs ≥ N^c" is a theorem (or explain why no such model is natural). Which invariants escape the model?
4. Rank your top 3 invariants where (i) X is factoring-equivalent, (ii) the generalisation is not NP-hard,
   (iii) nothing known rules out a polylog or sub-L[1/3] native algorithm. For each give a precise first
   research question and the strongest reason it might fail.

Be concrete and technical. Use web search for references (authors, year, arXiv/DOI). Label every claim
proved / known (with ref) / conjectured / heuristic. You may write scratch code in your working directory
(Python 3 with sympy is likely available) to check small cases; keep any computation tiny.
Output: a single markdown report (tables welcome).
