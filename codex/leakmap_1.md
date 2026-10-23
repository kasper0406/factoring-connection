You are a research mathematician / theoretical computer scientist. Blind round for agent `leakmap` of a research
project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you; scratch files
go in your current directory). Rigour and honesty matter more than a positive result. Use web search to check every
literature claim; give exact citations (authors, title, venue, year, theorem/lemma number where possible) and mark
anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (full proof given here)
/ HEURISTIC / CONJECTURED / OPEN / FALSE.

Background you should read first (in the repo):
- notes/factorbits_results.md §3 (Theorem PT: "pure trees" = finite adaptive decision trees whose nodes are integer
  straight-line programs with N-independent constants, evaluated mod N, with gcd/zero tests. Partition lemma: such a
  tree on N = pq either splits N or ends in the common "class" of p and q, where class(x) is the diagonal transcript
  [x | c_v] of a single prime. Consequence for iid factors: Pr[correct or split] <= max_b Pr[F=b] + 2s. PT-m, PT2.)
  and §3.3 (open classes F, G, H).
- notes/slpcirc_results.md §4.2 and §4.4 (Theorem SP sparse alignment; Cor SP-N on exponents polynomial in N).
- notes/primality_analysis.md §§0-1 and §5b (CRT criterion, order-pattern lemma, Lemma M2 "N-derived exponents are
  q-driven").
- notes/eps_results.md (only CF-prefix statistics of sqrt(N) leak eps = (-1/p), via A_{i-1}^2 = (-1)^i Q_i mod N).
- notes/selfred_results.md §3 (positional floor; AC0 sufficient criterion), notes/nonsmooth_results.md §§2, 5.

Setting. N = pq balanced semiprime, p, q iid uniform primes in (X, 2X]. A "leak" is a poly-time computation from N
whose output predicts some function F(p, q) of the factors (e.g. eps = (-1/p) on N = 1 mod 4) better than the
fiber floor, in an event where it does NOT split N. The project's map says known leaks come in three families:
(M) multiplicative / group order (p-1, p+1, ECM, class-group orders, relation methods; smoothness),
(A) archimedean / Diophantine approximation (Fermat, Lehman, Coppersmith known bits; always exponential so far),
(Q) quadratic / genus (ambiguous forms, SQUFOF infrastructure, Jacobi-type reciprocity).
The goal is to locate exactly where PT-type "information costs splitting" statements stop being provable.

QUESTIONS (attack each independently; we want theorems with proofs, or precise counterexamples)

1. N as an exponent (class G). Extend PT to trees whose SLPs may raise values to N-dependent exponents E(N)
   (a^{N-1}, a^{(N-1)/2^j}, Lucas V_{N+1}, x^{N} in Z/N[x]/(f) as in Frobenius/AKS tests).
   (a) Prove or refute: every non-split leak of such a tree comes from group-order coincidences of gcd(E(N), p-1)
       type (rate like gcd(N-1,p-1)gcd(N-1,q-1)/phi(N) for Fermat/MR), and find a theorem of the form
       Pr[correct or split] <= baseline + 2s + (explicit coincidence terms). Be careful about: adaptive trees of
       arbitrary depth; pure sub-computations of arbitrary depth feeding the bases; several exponent gates per node;
       extension algebras (unit group exponent p^d - 1 can exceed q); nested exponentiation ((a^N + 1)^N).
   (b) Is there any regime where such trees produce an ARCHIMEDEAN (family A) leak? (Hint to check, not to trust:
       N = q mod (p-1).)
2. Integer representatives. Allow the tree to read the integer representative in [0, N) of residues it computes
   (comparisons, floor division, bits). What does this add? Is the resulting model complete (simulates every
   poly-time algorithm on N)? If so, what is the strongest sub-model where something provable remains, and does
   every leak in it go through a short-vector / Diophantine event (family A) as for Lehman and Coppersmith? Relate
   to the Coppersmith n/4 known-bits barrier and to covering arguments (Lehman/Harvey N^{1/5}-type).
3. Continued fraction of sqrt(N), forms of discriminant 4N or -4N, infrastructure (class H). Do all leaks reduce to
   (Q) ambiguous/square forms or (M) relations? Formulate the strongest statement you can prove (e.g. about what the
   principal cycle reveals through genus characters ((-1)^i Q_i | p) = 1), and the strongest you can only
   conjecture. Is there a natural random model ("principal cycle = random walk subject only to genus constraints")
   and is it known/plausible? Literature on statistics of CF partial quotients of sqrt(N) vs factorization.
4. AC0 / finite automata on the bits of N (class F, secondary). By Linial-Mansour-Nisan an AC0 floor for eps would
   follow from bounds on E[eps(p) * (-1)^{sum_{i in S} bit_i(pq)}] for |S| <= polylog. Is this a type-II bilinear
   Walsh sum with ARBITRARY bounded coefficients on balanced ranges? Do Bourgain ("Moebius-Walsh correlation bounds
   and an estimate of Mauduit and Rivat", J. Anal. Math. 2013), Green ("On (not) computing the Moebius function
   using bounded depth circuits", CPC 2012), Mauduit-Rivat (Gelfond problem, Annals 2010; Rudin-Shapiro JEMS 2015),
   Muellner (Duke 2017) already contain the needed type-II estimate with the needed uniformity (saving b^{-omega(t)}
   for |S| <= t = polylog b)? Give exact statements and either derive the AC0 (or Thue-Morse / automata) floor or
   pinpoint the missing estimate.
5. Synthesis. Either a theorem "in model X every non-split leak is (M), (A) or (Q)" with X as large as you can make
   it, or an explicit computable invariant in classes F/G/H that leaks beyond the fiber floor without splitting via
   any of the three. Name the strongest natural model where the classification is open.
   Also: Jacobi symbols of PURE values, (c|N) = (c|p)(c|q) for an N-independent integer c given by an SLP — what
   does a pure tree augmented with such Jacobi-symbol tests leak about eps without splitting? (Note (-1|N) is known.)

Format: numbered sections matching the questions; each claim labelled; proofs in full; a final table "model ->
status -> leak family -> reference". Keep it under ~6000 words.
