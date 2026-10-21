You are a research mathematician / theoretical cryptographer. Blind round for agent `factorbits` of a research
project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you; scratch
files go in your current directory). Rigour and honesty matter more than a positive result. Use web search to
check every literature claim; give exact citations (authors, title, venue, year, theorem number where possible)
and mark anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (proof
given here) / HEURISTIC / CONJECTURED / OPEN / FALSE.

Background (optional to read): notes/trackE_analysis.md §5 (information floor: any function of N mod m carries 0
bits about (p mod m, q mod m) beyond the fiber {(u, c/u)}), notes/nonsmooth_results.md §§1, 5.4, 5.5 (equivalence
map; abelian/reciprocity aggregates depend only on N mod conductor; Prop RANK: φ(N) mod ℓ decides "ℓ* is a square
mod N"), notes/eps_results.md (no cheap predictor of ε = (−1/p) found), notes/periodfind_results.md §0.

Theme. Factoring is a change of representation. In the prime-exponent representation (an integer is its
exponent vector over the primes) factoring is trivial and multiplication is vector addition, but addition is hard.

QUESTIONS

1. Equivalence. State and prove cleanly, with exact oracle-call counts and the input-size convention that makes
   it true (careful: exponents written in binary make the prime-exponent encoding exponentially more succinct
   than binary for numbers like 2^(2^k); what must be bounded?):
   (a) an addition oracle in the prime-exponent representation factors N in ≤ log2 N calls (and a successor
       oracle x ↦ x+1 already suffices?); fewer calls with windowing?
   (b) conversely, addition costs one complete-factorization call (how many calls of a mere splitting oracle?);
   (c) the F_q analogue: Zech logarithms vs discrete logarithms, exact call counts in F_p and F_{p^k}.
   Then classify PARTIAL forms as "≡ factoring", "easy (poly time)", or "open / not known either way":
   output only the support / exponents mod m / largest prime factor / B-smooth part / ω / Ω mod 2 / μ /
   squarefreeness of a+b; linear functionals ⟨v, w⟩ of the exponent vector of a+b (which w are easy? is there a
   theorem that the cheaply computable completely multiplicative functions are essentially Hecke characters
   χ(n)n^{it}, e.g. for automatic sequences / AC⁰ / finite automata — Klurman–Kurlberg, Konieczny–Lemańczyk–
   Müllner, Schlage-Puchta, Green's AC⁰ Möbius result?); addition only for smooth summands (two summands vs many);
   addition with smooth result; addition modulo m; approximate addition (a factored integer within E of a+b).

2. Hardcore bits / partial information of factors. Pin down what is KNOWN. For each quantity computed from N
   alone, status "equivalent to factoring (with reduction)", "known easy", or "open (no reduction known either
   way)": individual bits of p (fixed position i; p = smaller factor), p mod small m, (−1/p) for N ≡ 1 mod 4,
   bit 1 of the smaller factor for N ≡ 3 mod 4, ⌊p/2^k⌋ for various k, parity of Ω(N), ω(N), d(N), μ(N),
   squarefreeness, φ(N) mod m, [ℓ | φ(N)] (Φ-hiding), quadratic residuosity (fixed base and random base),
   "N has a factor ≤ x", largest/smallest prime factor, #E(Z/N), r_2(N), a sum-of-two-squares representation.
   Literature to check and cite precisely: Alexi–Chor–Goldreich–Schnorr (SIAM J. Comput. 1988); Fischlin–Schnorr
   (J. Cryptology 2000); Håstad–Näslund (all RSA/Rabin bits hardcore, JACM 2004); Goldreich–Levin (STOC 1989) —
   does it give a "bit of the factors" hardcore predicate ⟨p, r⟩ for the function (p,q) ↦ pq, and exactly under
   which distributional hypothesis?; Rivest–Shamir 1985 (n/3 bits); Coppersmith 1996/97 (n/4 bits); Boneh–Durfee–
   Frankel 1998; Maurer 1995 "On the oracle complexity of factoring integers" (exact statement: εn arbitrary
   yes/no queries? conditional on what?); Woll 1987 "Reductions among number theoretic problems"; Shallit–Shamir
   1985 (number of divisors); Landau 1988 (square parts); Adleman–McCurley "Open problems in number theoretic
   complexity II" (ANTS 1994): which of these are listed open?; Bach–Shallit, Algorithmic Number Theory vol. 1;
   Bach–Miller–Shallit 1986 (σ); Kunihiro–Koyama 1998 (#E(Z/N)); Seurin 2014 (2-Φ/4-hiding); Cachin–Micali–Stadler
   (Φ-hiding); Jager–Schwenk (generic ring model, Jacobi symbol); Aggarwal–Maurer 2009 (generic RSA ⇔ factoring);
   Bernasconi–Shparlinski and Allender–Saks–Shparlinski (circuit lower bounds for squarefreeness / primality);
   Shparlinski's book "Cryptographic applications of analytic number theory" (any result on Boolean complexity of
   bits of factors of RSA moduli?); Henecka–May–Meurer 2010 and Paterson–Polychroniadou–Sibborn 2012 (noisy key
   bits); Heninger–Shacham 2009. Is there ANY theorem of the form "computing one fixed bit of p (or p mod 4, or any
   single non-trivial bit of the factorization that is a function of N alone, not using auxiliary randomness r) is
   as hard as factoring"? Or an explicit statement in the literature that this is open?

3. Middle ground. Representations R = poly-time function of N that leak factor information. Formalise an
   "information floor" theorem for broad natural classes and find the sharpest statement. Candidates I want you
   to attack independently: (a) R factoring through N mod m (with m ≤ (log N)^A, prime equidistribution); (b)
   "pure" generic-ring decision trees over Z/N (straight-line programs with integer constants independent of N,
   zero tests mod N, adaptive), possibly with side information N mod m; what is the exact leak? (c) finite automata
   / AC⁰ on the bits of N; (d) reciprocity-computable symbols. Which classes provably leak nothing, which leak
   only the local fiber (a function of p mod m), and which are open?

4. Leverage / self-reduction. Coppersmith needs ≈ n/4 bits of p. What is the minimal leakage that factors in
   polynomial time (bits of p in fixed positions; arbitrary yes/no queries)? Does a WEAK per-bit predictor for a
   FIXED bit position of the factor amplify to factoring, as ACGS amplify for RSA? Find or rule out (in a stated
   model) a random self-reduction: a randomized poly-time map N ↦ N' whose factor p' is a known function of p
   (e.g. polynomial maps h ∈ Z[x] with g(p) | h(pq) for all primes p, q). Explain precisely why ACGS / Blum–Micali
   amplification works for RSA/Rabin/DL and what is missing for bits of p. What does Goldreich–Levin give, and
   what does a noisy-bit decoder (Henecka–May–Meurer, Paterson et al.) need in terms of per-bit accuracy and
   independence?

Deliverable: a structured report (≤ 5000 words) with sections 1–4, a final status table (quantity | status |
reduction / reference), and a list of the 5 claims you are least sure of.
