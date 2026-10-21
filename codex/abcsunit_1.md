# abcsunit round 1 (BLIND): factoring relations as S-unit equations — abc / S-unit theory as barrier and as target

You are a research mathematician (Diophantine equations, S-unit equations, the abc conjecture, smooth numbers,
analytic number theory, lattice algorithms) collaborating on a speculative project about asymptotically faster
classical integer factorization. Rigour matters more than speed and honesty more than a positive result. Label every
claim PROVED (give the proof), KNOWN (cite precisely: authors, journal, year, theorem number — check with web search
when unsure; say UNVERIFIED if you could not check), HEURISTIC (name the assumption), model-PROVED (exact consequence
of a named model), or OPEN. Do not read other files of this repository; work from this prompt, your knowledge and web
search. Scratch computation (Python, exact integer arithmetic) goes in the current directory. Final answer: a
structured report (≤ 6000 words) with numbered sections matching the tasks.

## Notation and background

x = log N, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}). N = pq a balanced semiprime (also consider general N).
"y-smooth" = all prime factors ≤ y. S = a set of s primes (typically all primes ≤ y), Γ_S = {±∏_{p∈S} p^{e_p}} the
S-units of Q.

**Relations.** A Dixon/QS-type relation is a congruence u ≡ v (mod N), i.e. u − v = kN, where either both u and v are
y-smooth ("doubly smooth"), or one side is a known square/product and the other side v is small and y-smooth
("one-sided": Dixon x² mod N; QS (x+⌊√N⌋)² − N). F₂-dependencies among the exponent vectors give X² ≡ Y² (mod N),
which split N with probability ≈ 1/2 heuristically. NFS uses pairs (a, b) with a − bm (rational side, m a root of f
mod N with f(m) = N in base-m selection) and the norm F(a, b) (algebraic side) both smooth.

**Smoothness-size principle (heuristic).** An L[α] method needs its by-chance-smooth candidates to have log-size
≍ x^{2α}(log x)^{1−2α}; QS/Dixon (α = 1/2) test candidates of size N^{Θ(1)}; NFS (α = 1/3) has norms of log-size
x^{2/3}(log x)^{1/3}.

**The project's "product-listing primitive" Π.** Cheaply LIST (at amortised cost far below the box sieve) many
products of factor-base pieces congruent mod N to something small or to another split product. A doubly-split
version: output L[β]^{1+o(1)} pairs of y-smooth integers u ≠ v with u ≡ v (mod N), y = L[β], spanning F₂-rank
(1 − o(1))·π(y), total time L[β]^{1+o(1)}. A model computation (heuristic NFS cost model, fruitful dependencies,
early stopping allowed) says that an output-sensitive lister of doubly-smooth pairs for two linear forms of height
≈ √N with a POLYLOGARITHMIC factor base y = (log N)^κ would give heuristic polynomial-time factoring (a quantified
form of the Schnorr 1991 / Adleman 1995 lattice observation); listing ALL such pairs would be exponential, so a
prefix/early-stopping version is needed.

## Tasks

1. **Exact translation.** State precisely how a factoring relation maps to an S-unit or abc-type equation a + b = c:
   the factor base S, the role of N (is N or its prime factors part of S? what if k is or is not smooth?), the
   gcd normalisation, the size constraints. Which relation types fit (QS, Dixon, doubly-smooth congruences, the NFS
   rational side, Schnorr-type u − vN small), and which do not (the two-sided NFS relation with an algebraic norm)?
   Is there a uniform formulation (e.g. a kernel lattice of Γ_S → (Z/N)^*, or S-unit congruences in number fields
   modulo an ideal of norm N)?

2. **What abc and S-unit theory RULE OUT.**
   - (a) Does the abc conjecture (or effective/explicit forms: Stewart–Tijdeman, Stewart–Yu, Baker-type bounds) imply
     any lower bound on relation sizes, or any upper bound on relation counts, that constrains factoring methods?
     Try to prove a statement of the form "under abc, no relation family with candidates of log-size < f(N) exists",
     or prove that abc is vacuous here (and explain precisely why). Compare with the smoothness-size principle.
   - (b) Do S-unit counting bounds (Evertse 1984: ≤ 3·7^{d+2s} solutions of x + y = 1; Beukers–Schlickewei 1996;
     size-independent) limit how many relations a factor base of size s can supply? NFS/QS needs about s relations.
     Does any such bound bite, for any regime of s (L[β], polylog, O(log N), o(log N))?
   - (c) Lower bounds for the number of S-unit solutions (Erdős–Stewart–Tijdeman 1988; Konyagin–Soundararajan 2007):
     relevance?
   - (d) Is any barrier you find known in substance (e.g. in Lagarias–Soundararajan)? Be precise.

3. **What algorithmic S-unit solving would GIVE.** (a) If one could efficiently list all solutions of a + b = c in
   S-units with c = kN (k unrestricted, or k also smooth), or all smooth a, b with a ≡ b (mod N), what would factoring
   cost? (b) Compare with Π: are these literally the same target? Make the equivalence/containment precise. (c) Do
   the known S-unit algorithms (de Weger 1987/1989 Baker + LLL reduction; Tzanakis–de Weger; Smart; von
   Känel–Matschke; Alvarado–Koutsianas–Malmskog–Rasmussen–Vincent–West Sage solver) have ANY step that could be made
   N-dependent? Identify exactly which step fails for the N-condition (e.g. the p-adic logarithm reduction for an
   unknown or composite "prime"), and whether that failing step is itself equivalent to factoring. Relation to
   Schnorr's CVP approach?

4. **Lagarias–Soundararajan and smooth abc solutions (most important).** Lagarias–Soundararajan ("Smooth solutions to
   the abc equation: the xyz conjecture", JTNB 2011; and "Counting smooth solutions to the equation A + B = C",
   Proc. LMS 2012) show, conditionally on GRH, many solutions to a + b = c with a, b, c all y-smooth for
   y = (log c)^{8+ε}. That is a triply-smooth relation with tiny smoothness.
   - (a) State their results precisely (hypotheses, the exponent 8, the conjectured threshold, what abc implies, the
     counting asymptotic and its range, primitive vs non-primitive).
   - (b) Can their counting be made N-dependent: smooth a, b, c with c = kN, or smooth a, b with a ≡ b (mod N)? What
     does the count become, what inputs (GRH for which characters) are needed, and is any of it easier or harder than
     their theorem?
   - (c) Is anything constructive? If their counting extended to a congruence mod N AND were constructive, factoring
     would be (heuristically) polynomial time. Find EXACTLY where this fails: non-constructive counting? no
     N-dependence? the congruence condition? Is there any known efficient construction of polylog-smooth abc triples
     even WITHOUT N? (Polynomial identities, Pell equations/Lehmer, Prouhet–Tarry–Escott, ...)
   - (d) Quantify: what does the birthday/pigeonhole construction give, and what would a constructive version of the
     existence statement need to beat?
   Be careful and honest. This is the central question.

5. **Experiments.** Propose 2–4 small exact experiments (falsifiable predictions) at toy scale (N ≤ 2^40, ≤ 2 GB RSS)
   that test the counting predictions for smooth a + b = c with a congruence mod N, and the triply-smooth / Evertse /
   abc side. You may run small scratch checks yourself.

6. **Verdict.** Is "solve S-unit equations / smooth abc with an N-congruence efficiently" a SUFFICIENT high-profile
   target for fast factoring? Is it the same as, weaker than, or stronger than a doubly-smooth relation lister? Does
   abc/S-unit theory give any barrier at all? Subjective probability that this route yields an L[α<1/3] algorithm.
