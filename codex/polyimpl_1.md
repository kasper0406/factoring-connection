# polyimpl_1 — blind research round: consequences of hypothetical fast factoring, as a lens for blockers and ideas

You are a research mathematician (computational number theory + complexity). The repo at
/Volumes/git/factoring-connection is read-only for you (scratch: cwd). Context, read as needed:
README.md, APPROACHES.md (catalog of mechanisms/barriers), and the barrier statements in
notes/multimod_results.md (Thm A), notes/below16_results.md (T1/T2), notes/coverage38_results.md (Thm U,
Prop L), notes/dchworst_results.md §5 (B1–B4), notes/sieve2d_results.md (Prop 5),
notes/unbalanced_results.md (Prop O′), notes/trackC_analysis.md (O1–O4), notes/trackA_analysis.md (Thm A1),
notes/trackE_analysis.md (information floor), notes/eps_results.md, notes/nfsfruit_results.md,
notes/lemmaM_results.md, notes/hypD_results.md. Also Ward–Maestro, "Floors for deterministic integer
factoring", https://igent.ai/research/floors-not-bounds/factoring_report.pdf (Sept 2026).
Project frontier: deterministic worst case N^{1/5} (Harvey 2021, Harvey–Hittmeir 2022); our average-case
deterministic N^{1/6}polylog on balanced semiprimes (Theorem AV); rigorous randomized L[1/2,1] (Lenstra–
Pomerance 1992); heuristic L[1/3, (64/9)^{1/3}] (NFS).

IDEA. Suppose we had a fast factoring algorithm. Which theorems, conjectures, lower bounds or open problems
would be affected? Use the answer to find (i) new blockers (a hypothetical algorithm would imply something
implausible or something that is itself a hard open problem, so proving the algorithm would prove it) and
(ii) new ideas (consequences that look harmless or already true might be provable directly and leveraged).

Hypotheses (consider each separately):
 (H1) deterministic polynomial time; (H2) randomized polynomial time (classical);
 (H3) heuristic or rigorous L[1/3 − δ] or L[1/4];
 (H4) deterministic worst case N^{1/5−δ}; or average / almost-all N^{1/6−δ};
 (H5) algorithms computing only partial information: ω(N), μ(N)/squarefreeness, "is N a sum of two
      squares" (for N = pq: ε = (−1/p) with p ≡ q mod 4 issues), φ(N) mod small m, Jacobi-type bits.

Please answer, labelling EVERY claim KNOWN (exact reference: paper, theorem number) / PROVED (give the
proof) / SPECULATIVE / FALSE:

Q1 (consequence catalog). For each hypothesis, concrete consequences with precise statements:
 (a) computational NT: μ, φ, σ, squarefree part, ring of integers (Buchmann–Lenstra 1994: maximal order ⇔
     squarefree decomposition of disc), class groups/units, quadratic residuosity, square roots mod N,
     order finding / φ(N) (Miller; Bach–Miller–Shallit "Sums of divisors, perfect numbers and factoring"),
     discrete logs mod composites (Bach 1984), r_2(N), point counting, explicit constructions.
 (b) Deterministic prime finding (Polymath4; Tao–Croot–Helfgott "Deterministic methods to find primes",
     Math. Comp. 2012, which gets x^{1/2+o(1)}): EXACTLY what does a factoring oracle give, and at what
     speed? Known: with a factoring oracle for numbers ≤ x, one can find a prime in [x,2x] in polylog time
     IF certain conjectures hold (e.g. Polymath4 discussion: factoring + Cramér-type gaps; or the
     "factoring consecutive integers" / "find an integer with a large prime factor" routes). State what
     is unconditional. Compute precisely: would (H1) imply a new unconditional deterministic prime
     construction bound? Would (H4) (an N^{1/6−δ} or N^{1/5−δ} factoring algorithm, deterministic) imply a
     prime-finding bound better than x^{1/2+o(1)}? E.g. factor all integers in a short interval [x, x+y]
     and use a known theorem guaranteeing an integer there with a prime factor > x^{θ} (Baker–Harman–Pintz;
     Haugland; Merikoski; Jia–Liu largest prime factor in short intervals) to output a prime of size
     x^{θ}; what range of primes / what bound results? If (H4) would give something better than known, that
     is a BLOCKER (proving H4 proves it). Do the arithmetic carefully.
 (c) Complexity: Lipton's "straight-line complexity and integer factorization" (ANTS 1994) and the
     Shub–Smale τ(n!) / τ-conjecture connection (Blum–Cucker–Shub–Smale; Cheng 2003 "On the
     straight-line complexity of factorials"; Bürgisser "On defining integers and proving arithmetic
     circuit lower bounds" 2009): exact statements in BOTH directions (fast SLP for n! ⇒ factoring; does
     fast factoring ⇒ anything about τ(n!)?). Circuit lower bounds for μ / squarefreeness / primality
     (Allender–Saks–Shparlinski 2001; Green 2012 "On (not) computing the Möbius function using bounded
     depth circuits"; Bourgain on AC^0-type Möbius disjointness; Shparlinski's books). Sarnak's Möbius
     disjointness: μ ∈ P under (H1): any contradiction? Generic-ring results (Aggarwal–Maurer 2009
     "breaking RSA generically is equivalent to factoring"; Boneh–Venkatesan 1998; Jager–Schwenk 2013 on
     Jacobi symbols in generic rings; Damgård–Koprowski). Factoring-based PRFs/PRGs (BBS, Naor–Reingold)
     and natural proofs (Razborov–Rudich): which unconditional theorems would lose hypotheses; is ANY
     unconditional theorem contradicted; is any known lower bound for a RESTRICTED model of factoring
     violated (this tells which models any fast algorithm must escape)?
 (d) crypto: one paragraph only.
 (e) analytic NT: does fast deterministic factoring imply any new theorem about primes or smooth numbers
     (π(x) computation, prime gaps, explicit Linnik constants, verification of conjectures)? We expect
     "almost none"; check carefully and identify any exception.

Q2 (escape-signature table). For each barrier/model in the repo (the files above) and in Ward–Maestro's
 floors: state the model assumptions exactly; which assumption any polynomial-time OR deterministic
 sub-N^{1/6} algorithm must break; and which KNOWN algorithms break it (NFS/QS, ECM, Lenstra–Pomerance
 class groups (Schnorr–Lenstra / Seysen / LP 1992), Shor, Pollard–Strassen, Coppersmith, our D′/AV).
 Then identify the intersection of required escapes for a hypothetical DETERMINISTIC N^{1/6−δ} worst-case
 (or almost-all) algorithm, and whether any known family has that whole signature. Try to state a precise
 "combined barrier" theorem of the form: any deterministic algorithm that is local in sense X AND uses only
 group-order information in sense Y needs ≥ N^{c}. Prove it if you can.

Q3 (reverse direction). (i) Blockers: consequences of H1–H5 that look implausible / are open problems.
 (ii) Harmless consequences that might be provable directly and leveraged (decision-to-search,
 partial-information-to-factor reductions; e.g. does an oracle for ε, or for μ(N) on numbers of a special
 shape, or for φ(N) mod m, factor? what is known: Landau/Bach–Miller–Shallit σ ⇒ factoring; Shallit–
 Shamir "number-theoretic functions which are equivalent to number of divisors"; Woll 1987 on
 relativized complexity of μ, squarefreeness; Adleman–McCurley open problems list 1994).
 (iii) 3–6 concrete new research directions, each with a feasibility estimate, a first lemma to prove,
 and a kill criterion. Rank them. Prefer directions that are cheap to test or decisive.

Be concrete and quantitative; do the arithmetic for (b) explicitly with exponents. Output a single
markdown report.
