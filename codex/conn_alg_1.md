# Blind survey round: unexpected connections between ten fields of pure mathematics and integer factorization

You are a research mathematician and complexity theorist. Goal: find **unexpected connections** between the
fields below and classical integer factorization that could, even speculatively, lead to an asymptotically
faster algorithm than: heuristic NFS L_N[1/3, 1.902] (MNFS), rigorous L_N[2/5] (randomised), deterministic
N^{1/5} (Harvey–Hittmeir). Breadth and creativity are wanted, but verdicts must be honest and technical.

The repository is at the parent of your working directory (read-only; README.md and BARRIERS.md summarise
what is already known; notes/*.md hold details; use `grep -rl <keyword> ../../../notes/` or similar to check
prior coverage). Key known barriers you should test every idea against:
- **Leak-family map:** every known way to make p and q behave differently is (M) a group order / smoothness
  (p−1, #E(F_p), class numbers ...), (A) an archimedean approximation (costs exponential precision), or
  (Q) a quadratic genus/Jacobi symbol (only learns N mod small m type data).
- **Theorem PT:** methods that use N only as a modulus (ring operations in Z/N plus gcd) learn nothing
  about the factors except by splitting N.
- **Size–entropy law:** any "direct target" that would give exponent < 1/3 is poly-time factoring in disguise.
- **Theorem RA:** reciprocity invariants depend only on N mod m or are factoring-equivalent. Dedekind sums,
  Salié sums mod N, #X_0(N)(F_{ℓ²}), #E(Z/N) for a fixed curve, genus/area of X_0(N) are all factoring-equivalent
  (computing them is as hard as factoring). Single-curve a_N (trace of Frobenius at N for a fixed E/Q) is open.
- Constructing smooth numbers near kN costs L[2/3] per output for every constructor class tried.

## Fields (10)
1. Langlands program / automorphic forms (functoriality, base change, Hecke eigenvalues at composite level N, trace formulas).
2. Representation theory (characters of GL₂(Z/N), SL₂(Z/N), Weil representation, Deligne–Lusztig, modular rep theory, dimension formulas).
3. Motives / Hodge theory / periods (periods of varieties built from N, Kontsevich–Zagier periods).
4. Algebraic K-theory (K-groups of Z/N, Z[1/N], number rings; regulators; Lichtenbaum).
5. Arithmetic geometry (heights, Arakelov theory, abc-type statements, Chabauty, rational points).
6. p-adic analysis / perfectoid spaces / prismatic cohomology (Witt vectors, tilting, Frobenius lifts mod N, δ-rings).
7. Model theory / o-minimality (Pila–Wilkie counting, definability, Ax–Kochen; counting rational points on transcendental sets built from N).
8. Additive combinatorics (sum-product, Freiman, polynomial method, Croot–Lev–Pach, higher-order Fourier analysis).
9. Ergodic theory / homogeneous dynamics (Ratner, Littlewood, equidistribution of Hecke orbits, Duke's theorem).
10. Galois cohomology / class field theory (Selmer groups, Cohen–Lenstra, reciprocity, Brauer groups).

## Task
For EACH field give 3–6 concrete candidate connections. Go beyond the obvious. For each:
1. **Connection:** what object is built from N alone, and what computable feature would carry factor
   information — or what algorithmic primitive would replace sieving / linear algebra / smoothness.
2. **Prior art:** literature references (search the web; give authors, year, arXiv/DOI where possible) and
   whether it is already used in factoring.
3. **Verdict:** exactly one of RESTATEMENT (needs a factor to write down); AGGREGATE (depends only on N mod m);
   FACTORING-EQUIVALENT (or DL-equivalent); REDUCES-TO-<M/A/Q or named barrier>; KNOWN (already used in
   factoring, with reference); OPEN-WORTH-A-LOOK.
4. **Where it fails:** 2–5 sentences pinpointing the exact step that needs the factors, costs size N, or lands
   in a known family.

Then rank your **top 3 OPEN-WORTH-A-LOOK** items (if any survive) and give for each a precise first question
or lemma that could be attacked in a week. If none survive honestly, say so. Flag anything you believe is
surprising or not in the literature. Do not pad: a crisp negative verdict with the failing step is more useful
than a vague positive one. Output as markdown, one table per field, then the ranked shortlist.
