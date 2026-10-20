# analytic_1: blind round. ANALYTIC (non-relational) factoring algorithms

You are a research mathematician (analytic number theory, computational complexity) with web search.
The repository at /Volumes/git/factoring-connection is read-only context (an AI-assisted factoring
research project). Background you may cite instead of re-deriving: README.md ("Barriers"),
notes/nonsmooth_results.md (equivalence map of factoring-equivalent invariants; "class (a)/(b)/(c)";
Prop FROB; the Ramanujan-sum σ(N) expansion whose terms with modulus < p all vanish),
notes/trackE_analysis.md (information floor F1: anything that is a function of N mod m carries no
information; Coppersmith needs ≈ (log2 N)/4 coset-shaped bits of p).

Only CLASSICAL (deterministic or randomized) algorithms matter. Quantum ideas count only if they yield
a classical algorithm.

## Task

The paradigm: compute a factor-revealing quantity directly from an analytic formula (exponential
sums, theta functions, approximate functional equations, Riemann–Siegel / Odlyzko–Schönhage type
fast multi-evaluation, van der Corput / Weyl differencing), with no smooth numbers.

1. **Hiary's factoring algorithm, dissected.** Locate the exact reference(s): Hiary's deterministic
   O(N^{1/3+ε}) factoring algorithm and his fast algorithms for quadratic and cubic exponential sums
   (and Rubinstein's related N^{1/3} method, Integers 2013, if relevant). State precisely:
   (a) which sum reveals a factor, and how the factor is read off (bisection? exact value?);
   (b) how fast the sum is evaluated, with which subroutine (quadratic sums in polylog time via
   theta-type transformations; cubic sums in ~ K^{4/5} or whatever the correct exponent is);
   (c) WHY the exponent is 1/3: derive the cost balance exactly (Taylor-expansion block length vs
   number of blocks vs degree of the local polynomial vs error control). Which constraint binds:
   sum length, phase degree, precision, or the cost of degree-3 sums? What exponent would you get if
   degree-j exponential/floor sums of length K were computable in polylog time for every fixed j?
   If cubic sums of length K cost K^{β}, what exponent results?
2. **Better factor-revealing analytic quantities.** For each candidate give: the exact identity
   linking it to a factor; the minimal sum length/precision needed (information argument where
   possible); the best known evaluation cost; whether any known fast-summation technique beats
   N^{1/3}, N^{1/5}, or N^ε. Candidates (add others):
   - D(x) = Σ_{n≤x} ⌊N/n⌋ or Σ_{n≤x} ψ(N/n) (sawtooth), Σ_n e(hN/n); jump detection at divisors;
   - the Dirichlet hyperbola method / lattice points under xy = N, Voronoi-type summation;
   - r_2(N), r_4(N), class numbers via theta series near cusps / modular transformations;
   - L(1, χ_{−4N}) or L(1/2, ·) via approximate functional equation (cost ~ √conductor), and
     whether partial L-information (e.g. a few Euler factors, or low-precision values) leaks factors;
   - Σ_{d|N} d via Ramanujan sums / Kloosterman sums (avoid the known failure: all terms with
     modulus < p vanish).
3. **Dequantisation.** Shor's algorithm uses the Fourier transform of x ↦ a^x mod N. Is there any
   classical fast evaluation of sums like Σ_{x<X} e(h·(a^x mod N)/N), or Σ_x e(k x / r)-type period
   detectors, over structured ranges (geometric progressions, subgroups) that extracts period or
   factor information? Identify precisely where classical evaluation fails (e.g. the phase is not a
   low-degree polynomial / not "smooth", so van der Corput, Weyl differencing, Poisson summation and
   theta transformations give no compression; the sum has no small-height algebraic structure). Is
   there a clean theorem-level statement of the obstruction (query lower bounds, Cleve 2004; random
   self-reducibility; the sum's value being itself factoring-hard)?
4. **Subexponential analytic ideas.** Is there ANY analytic sum whose fast evaluation would give
   L[1/2] or better? E.g. sums whose phase becomes a low-degree polynomial after a substitution
   encoding N's factors; multi-dimensional theta sums (quadratic forms in many variables evaluated
   in polylog time via lattice reduction / Siegel modular transformations); hypergeometric or
   creative-telescoping sums. For each give a feasibility calculation and kill it if you can.
5. **Barrier map.** Known lower bounds/obstacles for analytic factoring: uncertainty-principle /
   large-sieve type limits on how short a sum can detect a divisor of size √N; precision Ω(log N)
   bits per term; information-theoretic lower bounds on sum length; degree growth under Taylor
   expansion; the fact that complete quadratic sums mod p encode class numbers. State the strongest
   statements you can PROVE (or cite exactly), with scope.
6. **New mechanisms.** Propose up to four NEW analytic candidates not in the literature as far as you
   can find. For each: mechanism, missing lemma, resulting cost, first kill test, your adversarial
   verdict. Kill as many as you can. Rank survivors.

Label every claim KNOWN (with exact reference), PROVED (give the proof), HEURISTIC, CONJECTURED, or
UNVERIFIED. Do not inflate. Exact citations matter: if you are not sure of a reference, say so.
Write your answer as a self-contained report in Markdown.
