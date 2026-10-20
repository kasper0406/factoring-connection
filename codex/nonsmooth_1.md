# nonsmooth_1: blind round. Factoring paradigms NOT based on smooth numbers

You are a research mathematician and complexity theorist with web search. The repository at
/Volumes/git/factoring-connection is read-only context (an AI-assisted factoring research project).
Relevant background files, if you want them: README.md, APPROACHES.md, notes/polyimpl_results.md
(sections 0, 2, 5), notes/trackE_analysis.md (information floor, Coppersmith price list). Do NOT
re-derive what those files already contain; cite them.

## Task

Every known classical subexponential factoring algorithm (QS, NFS, ECM, Dixon, Lenstra–Pomerance,
class-group methods) is smoothness-based. The aim here is to look systematically for factoring
paradigms that are NOT, and that could plausibly beat L[1/3] or even reach polynomial time. Be bold
when generating, ruthless when killing, and honest about known equivalences.

Organising principle: many invariants of N are polynomial-time equivalent to factoring (φ(N), λ(N),
σ(N), #E(Z/N), τ(N), coefficients a_f(N) of eigenforms (Bach–Charles, arXiv 0708.1192), class
numbers / class groups via genus theory, orders of elements, dimensions of spaces of modular forms of
level N, r_4(N), ...). Factoring gets faster exactly when one of them acquires an unexpectedly
efficient computation that does not secretly factor N first.

Please deliver:

1. **Equivalence map.** A table of factoring-equivalent (or factoring-implying) invariants. For each:
   exact reference (author, year, venue, theorem number where you can), the direction of each
   reduction (randomized/deterministic, under ERH or not, for all N or a density-one set), and
   whether a method computes it WITHOUT factoring for special inputs (prime arguments, special N, or
   via a global formula). Mark invariants that are hard even at prime arguments (e.g. h(−p)), since
   a route through them must solve a problem harder than factoring.
2. **Paradigm verdicts.** For each paradigm below: does the method implicitly need the factorization?
   Where exactly does it break for composite N? Is there a variant that survives? Give the precise
   failure reason in one or two sentences, with references.
   - Galois representations: Edixhoven–Couveignes compute τ(p) for prime p in poly(log p). For
     composite N, Frobenius at N is undefined. Any substitute (global trace of T_N, Bruin's work,
     a(N) mod ℓ for composite N, Eichler–Shimura over Z/N, Dokchitser-style Frobenius resolvents
     over (Z/N)[x]/(f))?
   - Point counting: Schoof/SEA over Z/N "pretending N is prime"; p-adic methods (Kedlaya, Lauder
     deformation, Satoh canonical lifts); Harvey's average-polynomial-time batch a_p for all p ≤ X
     (accumulating remainder trees) — is there any batch analogue that reveals factors of one N?
   - Trace formulas: Eichler–Selberg / Kronecker–Hurwitz class-number relations at composite N; the
     cost of computing τ(N) for composite N by any route avoiding factoring.
   - Dequantising Shor: classical query lower bounds for order finding (e.g. Cleve 2004), what
     non-black-box structure of a ↦ a^x mod N is known to be exploitable, elements of special
     structure whose orders are cheap.
   - Lattices/Coppersmith: minimum information needed; can any invariant supply it cheaply?
   - Algebraic geometry: #V(Z/N) for varieties with a global formula; zeta functions over Z/N.
   - Group/ring structure: generic ring results (Aggarwal–Maurer 2009, Jager–Schwenk 2013),
     reciprocity laws (Jacobi, power residue symbols: de Boer–Pagano) as the canonical non-generic
     structure; any classically exploitable structure beyond them.
   - Analytic / lattice-point methods (e.g. Hiary's deterministic N^{1/3} factoring, divisor
     summatory functions, exponential sums), CM / anomalous-curve special-purpose methods (Cheng's
     4p−1 method), and anything else you know that is genuinely non-smooth.
3. **New mechanisms.** Propose up to five NEW candidate mechanisms that are not smoothness-based, not
   in the literature as far as you can find, and not obviously equivalent to a known method. For each:
   the mechanism, the exact missing lemma that would make it work, the cost it would give, the first
   kill test, and your own adversarial verdict. Then kill as many as you can. Rank survivors by
   (probability of working) × (significance).
4. **Unifying barrier.** Is there a single principle explaining why the computable invariants carry no
   factor information (e.g. "everything computable from N by reciprocity is abelian, hence a function
   of N mod a conductor")? State the strongest version you can actually prove, and its scope.

Rules: label every claim KNOWN (with reference), PROVED (with proof), CONDITIONAL, SPECULATIVE,
FALSE, or UNVERIFIED. Verify references with web search; give arXiv/DOI links. Do not claim novelty;
say "not found by us" at most. Be concise: tables and short paragraphs. If you write scratch code,
put it under codex/scratch/nonsmooth_1/.
