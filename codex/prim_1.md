You are a research mathematician / theoretical computer scientist. Blind round for agent `prim` (phase 25) of a
research project on classical integer factorization (repo: /Volumes/git/factoring-connection, readable; your current
directory is a scratch dir you may write to; you may run python with /Volumes/git/factoring-connection/.venv/bin/python,
which has sympy, gmpy2, python-flint, numpy). Rigour and honesty matter more than a positive result. Use web search to
check literature claims; give exact citations (authors, title, venue, year, theorem/lemma number) and mark anything you
could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (full proof here) / HEURISTIC /
CONJECTURED / OPEN / FALSE.

READ FIRST (in the repo): notes/leakmap_results.md in full (especially §2 Theorem PT-X, §3.2 Lemma NE-C and Theorem
NE1, §3.4, §3.6, §7 Theorem PT-J, §8 "PRIM"), notes/factorbits_results.md §3 (Theorem PT, PT-m, PT2),
notes/slpcirc_results.md §1.2 (Prop R) and §4.1-4.2 (Lemma S1, Theorem SP, Kelley 2016).

SETTING. N = pq, p < q iid uniform primes in I = (X, 2X]. A "leak" is a poly-time computation from N predicting some
F(p, q) (e.g. eps = (-1|p) on N = 1 mod 4, where eps(p) = eps(q)) better than the fiber floor, on runs that do NOT
split N. Known leaks fall in three families: (M) group order / smoothness / gcd(E(N), p-1)-alignment; (A) archimedean
pins (small-height relations between p and q, factor sizes); (Q) quadratic genus symbols. Phase 24 (leakmap) proved
PT-X: in any tree with zero/gcd tests mod N, every non-split bit is charged to a "deviation": one prime's diagonal
transcript tau_p (answers [p | c_v(N)]) differing from a reference kappa(p) depending on p alone; Pr[correct or split]
<= B + 2s + 2 Pr[D and no split]. NE1 (one-layer gates y^{E(N)}, pure base y, node = low-degree polynomial in one gate
output) bounds deviations by root counting: deviation at p needs E(N) = t mod ord_p(y).

MODEL PRIM (the open frontier). Uniform programs of N-independent shape over a fixed free Z-algebra A (e.g. Z, or
Z[x]/(x^2 - bx + c)), computing in A/NA. Primitives: ring operations with N-independent constants; exponent gates
x -> x^{E(N)} with E from a fixed menu (alpha N + beta; floor((N - c)/2^j)); gate outputs may be RETAINED, NESTED (fed
into later gates, e.g. (a^N + 1)^N) and SUMMED (e.g. a^N + b^{N+1} + c); zero tests / gcd with N of computed values;
Jacobi symbols (c|N) of computed values c in Z/N. (Gates must be primitives: square-and-multiply with N-dependent
control is an N-dependent SLP, which is vacuous.)

Useful observations to check and use (we believe them; verify): at a prime p with A tensor F_p = F_p (A = Z), since
p = 1 mod (p - 1), x^{alpha N + beta} = x^{alpha q + beta} in F_p; for balanced p < q < 2p - 2, q mod (p - 1) = 1 + Delta
with Delta = q - p, and at q, p mod (q - 1) = p = 1 - Delta mod (q - 1). So PRIM at p is a program over F_p with the
hidden integer exponent 1 + Delta, and at q one with 1 - Delta. Each residue r mod (p - 1) is hit by <= 2 primes q in
I, so Pr_q[r_q in Z] <= 2|Z|/pi(I) ~ 4 ln X |Z|/(p - 1) for any set Z: deviations at p are controlled by the
"uniform-r density" rho_v(p) = |{r in Z/(p-1): node v deviates at exponent r}|/(p - 1), a pure function of p (this
fails in extension components F_{p^d}, d >= 2, where q is a short interval in Z/(p^d - 1)).

QUESTIONS
1. Two-layer gates. y = x^{E1(N)}, z = (y + c)^{E2(N)} with pure x, c (and variants: products with retained y, three
   layers). (a) With x a FIXED pure base, the zero set {r in Z/(p-1): (g^{a1 r + b1} + c)^{a2 r + b2} = w} mixes a Zech
   logarithm with the exponent: equivalently log(u) * log(u + c) = t-type "Zech-hyperbola" counting. Is there ANY
   known or provable nontrivial bound |Z| <= p^{1 - delta} (cf. Bourgain-Konyagin-Shparlinski 2008 on fixed points of
   the discrete log; Glebsky-Shparlinski 2010; Balog-Broughan-Shparlinski 2011; Cilleruelo-Garaev 2016 on x^x = c;
   Kurlberg-Luca-Shparlinski; anything on Zech logarithms)? Or can nesting manufacture a large zero set in r that is
   not a union of cosets r0 + mZ (m | p - 1) and not a short-interval pin? Try hard to construct one (e.g. with x of
   small order, c chosen so y + c lands in a subgroup, sums of nested gates). (b) With x a FRESH uniform random base
   (as in leakmap's NE1-R), prove the both-pass rate equals the (M) prediction up to an explicit Weil error (we expect
   ~ sqrt(p)/|H1| where H1 = image of the first gate), and find the first depth where Weil-type methods give only
   trivial bounds and why.
2. Summed gates. Node = t-sparse combination sum_m gamma_m prod_j y_j^{E_j(N)} with pure y_j and menu exponents. At p
   it is f(g^r) for a t-sparse Laurent f in Z/(p-1)-exponents (exponents depend on p via discrete logs). Using Kelley
   (LMS J. Comput. Math. 19A (2016), Thms 2.2-2.3) and Canetti-Friedlander-Konyagin-Larsen-Lieman-Shparlinski (Israel
   J. Math. 120 (2000)), prove a STRUCTURE statement: zeros = (union of cosets of subgroups of order >= L on which f
   vanishes identically; each forces Kelley partner congruences, i.e. order alignment (M)) union (a residual of size
   <= explicit bound). Check: is something like |residual| <= 2 n^{1 - c/(t log t)} provable by recursing on Kelley
   (vanishing coset of size C forces exponents to pair up mod C, so the sparsity halves on cosets)? Give the best
   explicit statement and the threshold on t (number of terms) beyond which it says nothing. What about nodes that
   apply a pure circuit of exponential degree to gate outputs, e.g. (g^N + 1)^M - 1 with M = lcm(1..B)?
3. Jacobi of compressed / computed values. (c|N) = chi_p(c) chi_q(c) never splits. PT-J: a single pure c leaks eps iff
   beta(c) beta(-c) != 0, beta(c) = E_x (c|x) over primes x in I. (a) Is there ANY c computable by a poly-size SLP
   (doubly exponential integer, e.g. Fermat numbers, a^M - 1, Lucas values, prod (x^2 + 1) over a range) with
   |beta(c) beta(-c)| non-negligible? Relate to known results on character sums / Legendre symbols of SLP-defined or
   sparse integers, and to the reverse direction (Jager-Schwenk generic Jacobi hardness). (b) For c computed with gates,
   e.g. c = P(y^{E(N)}) with P a pure low-degree polynomial: chi_p(P(y^{E})) - use Weil over subgroups to show only
   small order (M) or P = const * square (reducing to pure Q) can bias it, under a fresh-random-base assumption or in
   the deterministic case.
4. Depth hierarchy. If full PRIM is out of reach: prove the classification for depth-d nesting with explicit loss,
   identify the first d where the proof breaks and why (we suspect: deterministic bases break already at d = 2 via
   coupled exponents; fresh random bases break at d = 3 because two stacked square-root cancellations give a trivial
   bound). Confirm or refute these suspicions with arguments.
5. Literature. Has any generic-ring / hidden-order-group / generic-RSA result classified something like PRIM:
   Damgard-Koprowski (EUROCRYPT 2002), Rotem-Segev (CRYPTO 2020), Aggarwal-Maurer (EUROCRYPT 2009 / IEEE T-IT 2016),
   Jager-Schwenk (J. Cryptology 2013), Dachman-Soled-Loss-O'Neill (ITC 2024?), Leander-Rupp (ASIACRYPT 2006),
   Maurer-Wolf, Boneh-Venkatesan? In particular: is "generic ring algorithms that may raise to the power N" analysed
   anywhere?

DELIVERABLE: a structured answer, one section per question: statements, full proofs of whatever you prove, explicit
constants, counterexamples, and a verdict line per question. If you find an explicit PRIM program that leaks a factor
bit beyond the fiber floor without splitting, give it with a numerical check (small exhaustive prime pairs) - that would
be the most important outcome, but do not force it.
