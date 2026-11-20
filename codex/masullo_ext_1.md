# masullo_ext: can Masullo's orbital / adelic-θ counting be EXTENDED to the distinct-orbit count LP_exc^rel needs?

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch only in your cwd).
Goal: factoring below L[1/3]. Read: notes/masullo_results.md (and codex/masullo_1_{codex,gemini}.out), notes/
lowerbound_results.md (hatch 2), notes/subthird_results.md (Theorem KH ideal-relative form; the LP_exc^rel statement;
the Q(∛3) norms-1-and-8 example), notes/grdeg_results.md, notes/rank3_results.md. Paper text:
/Volumes/git/factoring-connection/codex/scratch/masullo/music.txt (pypdf; "=====PAGE=====" breaks): §1.2–1.4
(orbital Schauder bases, one-prime-at-a-time orbital decomposition, p-orbital traces vanishing for large p), §2.2–2.4
(global and S-local adelic θ-functions; Thm 2.4.2.1 lattice point counting), App. A.2 (orbital structures on adeles and
ideles), A.4 (S-unit counts, Thm A.4.1.1), Ch. 8 (periodic orbits), Ch. 11 (shell structure, orbital kernels).

The target (LP_exc^rel, simplified): fixed base field k of degree e (unit rank ≥ 1), K, K′ of relative degree n ~ √x
over k (x = log N), a linked O_k-module M = O_k w_1 + … + O_k w_m of relative rank m ~ x^{1/4} inside O_K × O_K′
(linked: the two evaluations to Z/N agree on M), candidates a ∈ {0,±1}^{em} (or a box). Need: an UPPER bound, uniform
in n and m, on the number of O_k^×-INEQUIVALENT candidates with 0 < |N(a)N′(a)| ≤ Y, log Y ~ √(x log x), whose parity
rows are independent modulo proper k-subspaces and lower-norm products — strong enough to force α ≥ 1/3 — OR a reason
it is false (a construction).

Orchestrator's initial view (test it, do not just confirm): Masullo's counts are full-rank (sum over F^× / ideals of a
ray class), fixed-field (constants depend on F), and count representatives within orbits; our set is a LOW-RANK
sublattice (linear section of the torus R_{K/k}G_m restricted to an O_k-span of m generators) with degree growing.
Candidate bridges to evaluate concretely:
 B1. Work over k: the norm form N_{K/k}(Σ a_i w_i) in m variables over k; adelic harmonic analysis on the torus
     T = R_{K/k}G_m (Batyrev–Tschinkel height zeta functions on toric varieties; Poisson summation on T(A_k)/T(k)) —
     the natural generalization of Masullo's GL1 machinery — applied to the linear section. Does any known theory
     (e.g. counting on linear sections of tori, norm-form equations via harmonic analysis, Peyre/Tschinkel) give
     uniform-in-n bounds?
 B2. Masullo's one-prime-at-a-time orbital chopping: does the per-prime decomposition give error terms with controlled
     dependence on degree (as opposed to C^n), which is exactly what we lack?
 B3. Sandwich: bound low-rank orbits by full-rank ray-class orbit counts (Hecke/Landau, with Masullo's GRH giving
     square-root errors) over the intermediate lattices Theorem KH uses (degree-T monomial spans that fill K×0 with
     index divisible by N). Compute what main term / error such a sandwich gives at the L[1/4] balance; is it enough?
 B4. Anything else in the paper (Ch. 8 periodic orbits, shell structure) that counts DISTINCT orbits.
Answer ≤ 1800 words, plain language, label proved / heuristic / speculative; for each bridge: the precise statement it
would give, the ledger at the L[1/4] balance (d ~ √(x/log x), r ~ x^{1/4}(log x)^{3/4}, log Y ~ √(x log x)), the
missing lemma, and a verdict. Finish with the single most promising line and a concrete next step (a lemma to prove or
a toy computation; pre-register + run if cheap).
