# slpcirc_1: blind round. Factoring via short straight-line programs (algebraic complexity)

You are a research mathematician and algebraic-complexity theorist with web search. The repository at
/Volumes/git/factoring-connection is read-only context (an AI-assisted factoring research project).
Background files, if you want them: README.md ("Barriers"), notes/trackC_analysis.md (§3: obstructions
O1-O4 for interval/prefix zero sets mod p), notes/polyimpl_results.md §4.1, notes/systeq_results.md
(Lemma CF), notes/nonsmooth_results.md §5.1 (Prop OA, order annihilation). Do NOT re-derive what those
files contain; cite them. Only CLASSICAL algorithms count. Be bold generating, ruthless killing, and
value honesty over a positive result.

## The paradigm

If prod_{i<k}(x+i) or n! mod N had an SLP of size polylog(N), factoring would be polynomial
(gcd(floor(sqrt N)! mod N, N); Lipton 1994/2009, Shub-Smale 1995, Strassen 1976, Bostan-Gaudry-Schost
2007). The tau-conjecture predicts no polylog SLP for n! over Z. Factoring needs much less: a short SLP
only MODULO N, for the specific N; or ANY polynomial F with a short SLP and a computable point x0 with
gcd(F(x0) mod N, N) nontrivial.

## Tasks

1. **Weakest sufficient statement.** Formulate precisely the weakest property of a family of
   short-SLP polynomials (or SLP families over (Z/N)[x], possibly with N-dependent constants, divisions
   with gcd checks, several variables) that yields (a) polynomial and (b) L[alpha], alpha < 1/3,
   factoring. Natural candidate: "for every prime p in the relevant range, the zero set of F mod p has
   density >= delta, while F is not identically 0 mod p" (then random x0 succeeds w.p. ~delta, and one
   needs also non-vanishing mod q). Is there an equivalence: classical poly-time factoring <=> existence
   of a uniformly constructible family of polylog SLPs with root density 1/poly(log N) mod every prime
   p <= sqrt N (and not identically zero)? Prove both directions, or say which direction fails and why.
   (Hint to check: a poly-time factoring algorithm is a Boolean circuit, not an arithmetic one over Z/N;
   can it be converted into a short SLP over (Z/N)[x] with the root-density property? Consider SLPs that
   use the integer N's bits as constants, x -> x^N, x^{N-1}, etc.)
2. **Known results**, with exact references (author, year, venue, theorem number): tau-conjecture and
   variants (Shub-Smale; Buergisser 2009 "tau => VP != VNP"-type implications; Koiran's real
   tau-conjecture and Wronskian work; adelic tau-conjecture (Phillipson-Rojas?); Cheng's WL-conjecture);
   degree/root lower bounds for SLPs (Strassen 1974 degree bound, Schnorr, Heintz-Sieveking 1980,
   Heintz-Schnorr; Borodin-Cook; Grigoriev; Baur-Strassen); "polynomials with many roots mod p" results
   (e.g. lower bounds for SLPs computing polynomials vanishing on a set mod p, Kayal's or Koiran's or
   Cheng's work; Cheng "Primality proving via one round in ECPP and one iteration in AKS"?; Kayal's
   "factoring via SLP"?; Cheng-Gao-Rojas-Wan "sparse univariate polynomials with many roots over finite
   fields" (2017); Bi-Cheng-Rojas "sub-linear root detection"; Kelley "roots of sparse polynomials over
   finite fields"; Kelley-Owen). Which lower bounds actually bind for the weak, mod-N, N-dependent
   version? (I expect: essentially none beyond log2(#roots) <= size, since tau-conjecture-type bounds
   are over Z/characteristic 0.)
3. **Constructions table.** Short-SLP polynomials with many roots mod p: x^n - 1, Dickson/Chebyshev
   D_n(x,a), division polynomials, Lucas/Lehmer, Lattes maps, iterated compositions (x^2+c)^{o k},
   x^p - x analogues, products/gcds/resultants of these, norms from extensions, sparse polynomials
   with structured exponents (Cheng-Gao-Rojas-Wan families), Artin-Schreier-type constructions, Deligne/
   Weil-bounded families (exponential sums give equidistribution, hence density ~1/2 for "is a square"-type
   predicates). For each: root count mod p, which group order governs it (p-1, p+1, #E(F_p), #A(F_p),
   Phi_k(p), ...), and whether it reduces to a known group-order method. Is there ANY short-SLP family with
   many roots mod p whose root count is NOT governed by a group order (or, a fortiori, density not tied to
   smoothness of a group order)? If you think every such family is "secretly" a group-order method, give
   the strongest argument: e.g. Chebotarev/Galois theory (if F splits into linear factors mod many p, its
   splitting field has many totally split primes...), dynamics (periodic points of an SLP map of degree d
   are bounded by ...), or Weil bounds (an absolutely irreducible curve y = F(x) has ~p points, so
   "F(x)=0 has >= p/L solutions" forces F to have >= p/L linear factors over F_p, and a degree-2^s
   polynomial computed by an SLP of size s ...). Note F(x) mod p has at most deg F roots but deg F can be
   2^s with s = SLP size; so root density needs deg F >= p * density, i.e. s >= log p - log(1/density).
4. **Barrier.** Strassen/degree bounds give only #roots <= 2^s, i.e. s >= log2 p for a density-1/polylog
   zero set: no barrier at all. What DOES rule out subexponential methods in a precise, restricted model?
   Candidates: (i) "fixed-p-independent" SLPs over Z (uniform constants) with zero set density delta mod
   every p in a range: a counting argument over the finitely many SLPs of size s vs the primes; (ii) a
   Galois/Chebotarev argument: for an SLP F over Z with no N-dependence, the density of p for which F has
   >= D roots mod p is governed by the Galois group of F; (iii) an "every short SLP is a composition of
   group-law steps" structural statement (false in general?). Try to prove the strongest honest barrier,
   with its exact scope, and say what it leaves open. In particular think about: F depending on N (e.g.
   constants N, x^N), divisions with gcd, and the multivariate case.
5. **Small search proposal.** Suggest an exhaustive tiny-SLP experiment (size <= 6-8, constants from a
   small set, ops +,-,*) that would detect polynomials whose zero-set sizes mod many primes p are NOT
   explained by group orders, and the exact statistic to test.

Output: a structured report. Label each claim PROVED / KNOWN (with reference) / CONJECTURED / UNVERIFIED.
Final section: your single best guess at the "true answer" for this paradigm, and the single most
promising open question.
