You are a research mathematician / theoretical cryptographer. Blind round for agent `selfred` (phase 22) of a
research project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you;
scratch files go in your current directory). Rigour and honesty matter more than a positive result. Use web search
to check every literature claim; give exact citations (authors, title, venue, year, theorem/problem number where
possible) and mark anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED
(proof given here) / HEURISTIC / CONJECTURED / OPEN / FALSE. Negative literature findings ("no result found") must
be stated as such, not as facts.

Background you may read (optional): notes/factorbits_results.md (esp. §2.3 status table, §3 Theorem PT for "pure"
generic-ring decision trees, §4.2 Lemma PSR: for FIXED g, h in Z[x], g(p) | h(pq) for all large p and all primes
q in (p, 2p) forces g = c·x^j), notes/trackE_analysis.md §5 (information floor: a function of N mod m carries no
information about (p mod m, q mod m) beyond the fiber {(u, c/u)}), notes/nonsmooth_results.md §5.3–5.5.

Notation: N = pq balanced semiprime, p < q, n = log2 N. epsilon = (−1/p) on N ≡ 1 (mod 4).

QUESTIONS

1. Self-reductions of factoring. RSA has c ↦ c·r^e (secret moves to x·r mod N), which is why ACGS / Håstad–Näslund
   amplification works. Is there a polynomial-time computable (possibly randomized, with auxiliary output) map
   N ↦ N′ whose prime factors are KNOWN, VARIED functions g_r(p) of p other than c·p^j?
   (a) Go through the candidates and classify each (works / covered by a barrier / exponential size /
       factoring-equivalent to compute): multipliers kN; N′ = f(N) for polynomials f (including families whose
       degree and coefficients grow with n, randomized coefficients independent of N); norms Norm(α) from number
       fields with α built from N; resultants and discriminants of polynomials with a known root mod N; values
       Φ_m(N) and Φ_m(x) at x related to N; group orders #E(Z/N), torus orders (p+1)(q+1), p²+p+1 via
       Lucas-type sequences; divisibility sequences a^N − 1.
   (b) Prove barriers that are as broad as you can. Two specific things we want checked or improved:
       (i) A height-free / poly-degree version of PSR: h in Z[x] of degree ≤ D with ARBITRARY coefficients
           (independent of N), p prime, Q a set of primes in (p, 2p) of relative density ≥ δ. If G | h(pq) for all
           q in Q, what can G be? (We expect: G divides cont(h)·p^∞ times an L-smooth number with L = poly(D, n, 1/δ),
           with a bound on the information about p it carries.)
       (ii) A rigidity statement for N-DEPENDENT maps: if any poly-time map outputs N′ together with a description
           of polynomials / rational functions A, B, C of degree ≤ poly(n) and poly-size coefficients such that
           N′ = A(p)·B(q)·C exactly (or up to an unknown cofactor of size ≤ 2^{O(log n)}, or approximately), then
           either A(x)B(N/x) is constant in x (the monomial/multiplier case) or N can be factored in polynomial time
           from (N, N′) directly. Is this right? Exact hypotheses? What escapes it (modular moves like (a·p mod M),
           non-algebraic g, factors hidden in large unknown cofactors)? Explain precisely why the RSA move
           x ↦ x·r mod N escapes the analogous argument and the integer equation pq = N does not.
       (iii) Meta-barrier: is it correct that an unconditional proof that NO polynomial-time self-reduction of the
           desired kind exists would imply factoring ∉ P (e.g. because with a factoring algorithm,
           N ↦ (p+1)(q+1) is one)? Same question for "no reduction from factoring to epsilon / squarefreeness".
   (c) Literature on random self-reducibility of factoring and instance hiding: Abadi–Feigenbaum–Kilian (JCSS 1989),
       Feigenbaum–Fortnow (SIAM J. Comput. 1993), Beaver–Feigenbaum, Bogdanov–Trevisan, and any result on
       random self-reducibility of factoring restricted to smooth parts, Blum integers, or per-modulus (Rabin
       square roots over random residues). What exactly is known? Is "factoring is not known to be random
       self-reducible (over the modulus)" the accurate statement? Does Feigenbaum–Fortnow say anything about
       factoring (its decision version is in NP ∩ coNP)?

2. Open partial-information bits. For each of: ω(N), μ(N), Liouville λ(N), squarefreeness, rad(N) (and squarefree
   part / squarefree decomposition), p mod 4, epsilon, φ(N) mod ℓ, fixed-base quadratic residuosity:
   give the precise current status: known reductions in either direction (with references), any oracle/relativized
   separation or restricted-model barrier, fast algorithms in special cases (e.g. N = p²q, Booker–Hiary–Keating
   "Detecting squarefree numbers" Duke 2015 — what exactly does it achieve and at what cost?). Give the exact
   statements and numbering in Adleman–McCurley "Open problems in number theoretic complexity, II" (ANTS-I, LNCS
   877, 1994) and Adleman FOCS 1994 for squarefreeness / squarefree part / ω; Woll "Reductions among number
   theoretic problems" (Inf. Comput. 72, 1987); Landau "Some remarks on computing the square parts of integers"
   (Inf. Comput. 78, 1988); Shallit–Shamir (IPL 1985); Chistov / Buchmann–Lenstra (ring of integers ⇔ largest
   squarefree divisor); Bach–Shallit's book. Also: given an oracle for rad on arbitrary inputs, is the full
   squarefree decomposition ≤ log N calls? Given an oracle for a function of the prime signature only (ω, μ, λ,
   squarefreeness, d), what can a reduction do with queries X computed from N (note gcd(X, N^∞) = p^a q^b splits N
   when a ≠ b)? Is anything known about class-number congruences (e.g. h(−4N) mod 8 via Rédei) encoding epsilon?

3. Theorem PT's open classes. Are there automata, AC⁰, or "reads N's integer structure" algorithms that leak factor
   information on balanced semiprimes? Use known results: Allender–Saks–Shparlinski (JCSS 2001, what exactly do
   they prove about squarefreeness and AC⁰[p]?), Bernasconi–Shparlinski, Green (μ ⊥ AC⁰), Konieczny–Lemańczyk–
   Müllner, Müllner's automatic-sequences-Sarnak theorem, Mauduit–Rivat, Bourgain's binary digits of primes.
   Specific question: is it provable that ANY function of N mod m (or of the low j bits of N) has advantage o(1)
   for predicting epsilon (or any function of p mod m) as long as m ≤ N^{1/4 − o(1)} under GRH, and unconditionally
   for m = 2^j ≤ N^{c} using equidistribution of primes in progressions to prime-power moduli (Gallagher 1972;
   Barban–Linnik–Tshudakov 1964)? State the exact theorems used. Note p mod m with log2 m ≥ n/4 + O(log n) factors
   N (Coppersmith / Boneh–Durfee–Frankel), so the two thresholds would meet at N^{1/4}. Can an AC⁰ floor for epsilon
   be proved via Linial–Mansour–Nisan plus bilinear (type II) bounds for Σ_{p,q} χ4(p) e(b·pq/2^k)? How far can
   you get?

4. Any positive leak, however tiny: give it with the exact mechanism, and check it is not a leakage artefact.

Deliver a structured report: (1) verdict; (2) self-reduction candidates table + barriers with full proofs;
(3) status table of the open bits with precise references; (4) PT classes; (5) what you could not verify.
