# newfamily_1: BLIND round — a fourth family of factor leakage

You are a research mathematician (analytic/algebraic number theory, complexity). The repo
/Volumes/git/factoring-connection is READ-ONLY context; write scratch files only under the current directory.
Do NOT read notes/newfamily_* (you are blind to our own list). You may read notes/factorbits_results.md §2–3
(Theorem PT and the "fiber floor") and notes/nonsmooth_results.md if useful.

## Setting
N = pq, balanced RSA semiprime, p, q unknown. Every known way to get *any* information about p beyond what
N mod m (m small) already determines falls into one of three families:
- (M) group order / smoothness (p−1, ECM, class-group/relations, NFS);
- (A) archimedean/Diophantine (Fermat, Lehman, Coppersmith lattices);
- (Q) quadratic/genus (ambiguous forms of disc −4N, SQUFOF, infrastructure).
"Fiber floor": a quantity R(N) that factors through N mod m, or is a product of local symbols
χ(p)χ(q) for a character χ of fixed conductor (the Jacobi symbol (a/N) = (a/p)(a/q) being the archetype), carries
no information beyond the floor.

## The question
Produce your OWN list (aim for 15–30 entries, be wild, cross-domain) of candidate invariants I(N) that are
**computable in polynomial (or subexponential) time from N alone** and that might carry information about the factors
beyond the floor (target: e.g. ε = (−1/p) for N ≡ 1 mod 4, or p mod 8, or any bit of p), WITHOUT going through a smooth
group order. Prioritise objects computable by reciprocity / Euclid-type recursion (the way the Jacobi symbol escapes
generic-ring hardness): Dedekind sums, generalised Dedekind/Rademacher sums, Gauss sums G(a;N) via Landsberg–Schaar,
power residue symbols (a/N)_n (de Boer–Pagano), modular symbols, eta multipliers, theta transformations,
Kloosterman/Salié sums, half-integral weight coefficients, quantum invariants at e^{2πi/N}, dynamical periods, etc.

For EACH entry give:
1. Definition, and a concrete poly-time algorithm from N (with reference) or "no known poly-time algorithm".
2. A verdict, with a 2–6 line argument:
   - AGGREGATE (reduces to the fiber floor: show the formula as a function of N mod m and/or products of local
     symbols at p and q, e.g. via the multiplicativity/twisted-multiplicativity formula);
   - FACTORING-EQUIVALENT (give the reduction and reference);
   - NOT COMPUTABLE (best known algorithm needs the factorization);
   - GENUINELY OPEN (poly-time computable AND not visibly an aggregate). These are the gold. Be precise and skeptical.
3. For the open ones: a concrete numeric test we can run on 10^5 random balanced semiprimes of 32–64 bits to
   detect correlation with ε = (−1/p) (or with p mod 8) beyond the N mod 2^8 floor.

Also attempt a THEOREM of the form: "every function f(a, N) computable by a Euclid-type recursion with a cocycle /
reciprocity law of type R is a function of (N mod m, local symbols)". Find the right class R and prove it, or exhibit a
counterexample class (a reciprocity-computable function that is not a local aggregate). In particular decide:
- Is the Dedekind sum s(a, N) (computable by reciprocity s(a,c)+s(c,a) = −1/4 + (a/c + c/a + 1/(ac))/12) an aggregate?
  Note s(a,N) is NOT multiplicative in N; does it nonetheless reduce to the floor? (Hint: 12N·s(a,N) mod N, the
  Rademacher φ function, the Dedekind eta multiplier formula and its relation to the Jacobi symbol.)
- Do the Gauss sums G(a;N), Kloosterman K(a,b;N), Salié S(a,b;N) have known poly-time algorithms (exactly, or to
  1/poly precision)? Salié has a closed form via square roots mod N — what does that imply?

Be honest: "I believe X is an aggregate because ..." is fine; flag every claim as proved / heuristic / unverified.
Output a single markdown document: the candidate table first (one line each: name | computable? | verdict |
reason), then details, then the theorem attempt, then your top-3 open candidates ranked.
