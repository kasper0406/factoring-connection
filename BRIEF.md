# Research agent brief: find a simple classical connection for integer factorization

**Prepared:** 17 October 2026  
**Mode:** Active mathematical research, with exact computational experiments  
**Primary target:** A classical polynomial-time factorization algorithm, or a concrete new mechanism that materially advances toward one

## 1. Mission

Investigate whether a simple connection between explicit integer arithmetic and compressed algebraic representations can produce an efficient classical factorization algorithm.

The working hypothesis is that the eventual mechanism may be conceptually simple and may not resemble sieving or classical simulation of a quantum algorithm. Treat this as a productive hypothesis, not an established fact. Do not dismiss a direction because factoring has been studied for a long time. Equally, do not mistake a short identity, a short program describing an enormous object, or successful toy examples for an efficient algorithm.

**Do actual research.** Read relevant primary sources, formulate candidate mechanisms, derive consequences, implement the most testable constructions, find counterexamples, and revise. A literature survey or a statement that the problem is unresolved is not the requested deliverable.

Begin with the polynomial lead below, but do not become attached to it. You may replace it with a better representation if you explain what new operation becomes possible.

The central question is:

> Which factor-revealing feature has a compact exact representation, and what operation would let us extract that feature without expanding the representation?

### Problem and computational model

Input is an integer `N > 1`, in binary. Write `n = N.bit_length()` and `R_N = Z/NZ`.

The initial hard family is `N = p*q`, with independently generated, distinct primes of similar bit lengths. The factors are not given to the candidate algorithm. Also test prime powers, unbalanced factors, and composites with more than two prime factors.

The first algorithmic goal is to return any `d` with `1 < d < N` and `N % d == 0`. Recursive splitting and primality testing turn a sufficiently uniform splitting algorithm into full factorization; account for those costs when making a final claim.

Use ordinary classical bit complexity. Deterministic polynomial time or randomized expected/bounded-error polynomial time would satisfy the main objective; state which one a proposal targets. Count preprocessing, repetitions, memory, bit precision, coefficient lengths, and representation conversion. A polynomial bound in `N`, `p`, a polynomial degree near `sqrt(N)`, or a group order is not a polynomial bound in `n`.

A restricted-family improvement, a new constructive lemma, a useful reduction with a smaller missing primitive, or a carefully delimited obstruction is worthwhile intermediate progress. Label it accurately and use it to choose the next experiment.

## 2. Mathematical starting point

The elementary claims in this section are supplied with proofs. They are research scaffolding, **not claims of novelty**. Check them independently and search for prior formulations before attributing originality.

### 2.1 The polynomial and its first surviving coefficient

Define

\[
F_N(X)=(1+X)^N-1-X^N\in R_N[X].
\]

The prime/composite polynomial identity is the starting point of the AKS primality test; see [AKS], Lemma 2.1 and Section 2. AKS uses suitable low-degree quotient rings to obtain a primality test. It does not supply the factor-extraction operation sought here.

Let `N` be composite and let `p` be its smallest prime divisor. For `1 <= k < p`,

\[
\binom Nk=\frac{N(N-1)\cdots(N-k+1)}{k!}\equiv0\pmod N,
\]

because `k!` is invertible modulo `N`.

At `k = p`, use the exact integer identity

\[
\binom Np=\frac Np\binom{N-1}{p-1}.
\]

All denominators in the product for the second binomial coefficient are units modulo `N`, so

\[
\binom{N-1}{p-1}
=\prod_{j=1}^{p-1}\frac{N-j}{j}
\equiv(-1)^{p-1}\pmod N.
\]

For odd `p`, the sign is positive. For `p = 2`, `-N/2` and `N/2` are the same residue modulo `N`. Therefore the least nonnegative residue is exactly

\[
\boxed{[X^p]F_N=N/p\pmod N.}
\]

Consequently,

\[
\boxed{\min\{k:[X^k]F_N\ne0\}=p.}
\]

This holds for every composite `N`, including prime powers, not just squarefree semiprimes. Since `p < N`, the removed endpoint `X^N` causes no issue. For prime `N`, `F_N` is the zero polynomial instead.

Examples:

| N | First surviving degree p | Coefficient residue N/p |
|---|---:|---:|
| 9 | 3 | 3 |
| 15 | 3 | 5 |
| 35 | 5 | 7 |
| 45 | 3 | 15 |
| 49 | 7 | 7 |
| 143 | 11 | 13 |

The expression for `F_N` has an arithmetic circuit with `O(n)` gates using repeated squaring. This is a statement about circuit description size, not expanded polynomial size or the cost of coefficient extraction.

### 2.2 A weaker extraction task already suffices: any surviving index

Write `a_k = binom(N,k) mod N` for an interior coefficient. The integer identity

\[
k\binom Nk=N\binom{N-1}{k-1}
\]

gives

\[
k a_k=0\quad\text{in }R_N.
\]

Thus, for `0 < k < N`,

\[
a_k\ne0\Longrightarrow1<\gcd(k,N)<N.
\]

**We do not need the first coefficient, or even its value. Finding any index of a nonzero coefficient gives a factor immediately.**

Target the weaker problem whenever it permits a simpler construction:

> Given the compact circuit for `F_N`, find one exponent in its nonzero support without enumerating the support.

### 2.3 Why individual coefficient queries do not provide the missing information

Suppose an algorithm asks for `a_k` at an explicitly chosen interior index `k`.

If `gcd(k,N) = 1`, its answer is already known to be zero. Otherwise, computing that GCD has already found a nontrivial factor.

Therefore an individual-coefficient oracle can be simulated up to factor discovery: compute the GCD first; either return zero or terminate with a factor. This does not prove factoring hard. It shows that **merely making coefficient evaluation fast is insufficient for this polynomial**. The useful operation must act on many indices at once, reveal a support index, or exploit different arithmetic information.

This is a particularly important filter for proposals based on random coefficient sampling, a binomial-coefficient library, or a fast high-order derivative at one guessed index.

### 2.4 Prefix testing and its factorial connection

For `K >= 1`, define

\[
Z(N,K):=[F_N(X)\equiv0\pmod{N,X^K}].
\]

For composite `N`,

\[
Z(N,K)\Longleftrightarrow K\le p
\Longleftrightarrow\gcd\bigl(N,(K-1)!\bigr)=1.
\]

For finding `p`, it suffices to use `2 <= K <= floor(sqrt(N))+1`. A uniform polynomial-time prefix predicate would recover `p` by `O(n)` binary-search queries. A randomized predicate must have sufficiently amplified error for the adaptive search.

This is an exact connection to bounded divisibility and modular factorials. It is not, by itself, a speedup. Investigate whether the polynomial representation makes the unit/nonunit predicate cheaper than evaluating the factorial, and whether factorial or recurrence techniques suggest a better compression. Read [BGS] and the related classical constructions [Harvey], [HH].

Do not require full factorial evaluation unnecessarily: the predicate asks only whether it is a unit modulo `N`. Conversely, do not treat an efficient unit predicate as already available.

### 2.5 Ordinary differentiation loses the signal

\[
F_N'(X)=N(1+X)^{N-1}-NX^{N-1}=0\quad\text{in }R_N[X].
\]

Every positive-order ordinary derivative is therefore zero. A characteristic-zero method based on the logarithmic derivative or derivative ratios cannot simply be transferred to this ring.

Hasse derivatives retain coefficient information:

\[
D^{[k]}F_N(0)=[X^k]F_N.
\]

But evaluating them at individually selected orders has the limitation in Section 2.3. A potentially useful derivative construction must aggregate orders, preserve additional lifted information, or access a different invariant.

### 2.6 Scalar evaluations and sparse recovery have specific traps

A nonzero polynomial over `R_N` can define the zero function on all scalar residues. For example, `F_561(a) = 0 mod 561` for every scalar `a`: modulo each of `3, 11, 17`, exponentiation to `561` acts as the identity. Yet `F_561` is not the zero polynomial, by Section 2.1.

Therefore scalar random evaluation in `R_N` is not automatically a sound polynomial identity test. Distinguish formal polynomial identity from equality as functions on a finite ring.

Nor is the expanded support polynomially small in general. For `N = p*q`, `p < q`, reduction modulo `q` gives

\[
F_N(X)\equiv(1+X^q)^p-1-X^{pq}\pmod q.
\]

For every `1 <= j < p`, `binom(p,j)` is nonzero modulo `q`, because its prime divisors are at most `p`. Thus `F_N` has at least `p-1` nonzero coefficients. On balanced semiprimes this is exponentially many in `n`.

The relevant opportunity is structured or partial support recovery, not an unsupported assumption of ordinary sparsity.

## 3. Priority research tracks

These are starting directions, not a compulsory architecture. First investigate A and B, connect them to C and D, and keep one genuinely independent alternative alive. The proposed modifications below are hypotheses, not results established by the cited papers.

### Track A — Binary coefficient sections and compressed witness extraction

For any polynomial `P`, define

\[
E(P)(Y)=\sum_j[X^{2j}]P\,Y^j,
\qquad
O(P)(Y)=\sum_j[X^{2j+1}]P\,Y^j.
\]

Then

\[
P(X)=E(P)(X^2)+X\,O(P)(X^2).
\]

If `P` is nonzero, at least one section is nonzero. Repeatedly selecting a nonzero section reconstructs the bits of a support index in `O(log degree(P))` levels.

There are elementary gate rules. Writing `f=f_e(X^2)+Xf_o(X^2)` and similarly for `g`,

\[
E(fg)=f_e g_e+Yf_o g_o,
\qquad
O(fg)=f_e g_o+f_o g_e.
\]

Addition is componentwise. A single section transformation has constant-factor circuit overhead with appropriate sharing. Repeating it `O(n)` times need not: a constant-factor increase at every level can still become exponential.

**Precise conditional target:** Construct a representation for all sections reached along one search path that remains polynomial in `n`, together with a polynomial-time formal zero test for those representations over `R_N`. This would turn the support-index reduction into a factorizer.

Tasks:

1. Derive and implement section transformations on arithmetic circuits, preserving shared subexpressions. Measure growth in unique nodes, not printed expression length alone.
2. Exploit the special repeated-power structure of `(1+X)^N`. Investigate matrix recurrences, rational generating functions, and reusable formulas for sections; do not default to a generic circuit expansion.
3. Determine what information a zero test actually needs. Test whether a weaker certificate or a small family of quotient evaluations can select a nonzero child without reconstructing either child.
4. Study coefficientwise products with interval or bit-pattern masks. The mask may have a short description, but prove or measure the cost of applying it.
5. Look for a class closed under both sectioning and zero testing. It need only contain the required descendants, not arbitrary arithmetic circuits.

Record separately the cost of constructing a section, representing it, evaluating it, and testing its formal identity. Making only one of these cheap is partial progress, not the full algorithm.

**First experiment:** Follow both diagnostic paths chosen with known support information and operational paths chosen without factors. Diagnostic paths measure representation growth only; they are not factoring results.

### Track B — Quotient sketches that retain useful exponent information

Compute compact projections such as

\[
F_N(X)\bmod(X^r-1),\qquad
F_N(X)\bmod(X^r-c),
\]

with polynomially bounded `r` and cheaply chosen `c`. Repeated squaring evaluates each projection without expanding degree `N`.

Modulo `X^r-1`, a coefficient at position `j` is the sum of original coefficients at `j, j+r, j+2r, ...`. A unit-valued sum can mix contributions divisible by different hidden prime factors. A nonzero projected coefficient is therefore not automatically a factor-revealing coefficient of the original polynomial.

Modulo `X^r-c`, folded terms acquire weights `c^t` according to the quotient of their exponent by `r`. Investigate whether several weights or several quotient rings expose a useful exponent, a nonunit, or enough partial information for another method.

Concrete questions:

- Can paired sketches isolate one hidden characteristic's contribution without knowing that characteristic?
- Can small-degree resultants, subresultants, annihilators, ranks, or failed inversions reveal a nontrivial GCD?
- Can a small family of weighted projections identify one support exponent modulo enough small moduli without enumerating exponentially many matching choices?
- Can the relation `F_N' = 0` constrain the information discarded by reduction? Differentiation before and after reduction need not commute: retain the quotient terms in the derivation.
- Can a sketch answer a prefix or binary-section question rather than reconstruct the full polynomial?

Read [AKS] for the distinction between proving compositeness and extracting a witness, and [Sparse] for interpolation mechanisms and their assumptions. Characteristic, number of terms, roots of unity, discrete-log steps, and field-versus-ring behavior must all be audited before transfer.

Do not use generic sketch collisions to claim an impossibility theorem for this family: the decoder receives `N` and can exploit its explicit arithmetic. A negative result must identify the restricted sketch/decoder model it actually addresses.

**First experiment:** Exhaustively inspect several small families of weighted projections, their coefficient GCDs, and basic invariants on tiny semiprimes. Derive the reason for any successful pattern before scaling it.

### Track C — Factorial unit tests and compressed recurrence evaluation

The prefix predicate is precisely a factorial-unit predicate. This links the polynomial lead to algorithms for products of consecutive terms and linear recurrences with polynomial coefficients [BGS].

Investigate whether unit/nonunit information admits a simpler recurrence than the factorial value itself. Candidate tools include product trees, block products, multipoint evaluation, cancellation-aware recurrences, and transformations that isolate an interval containing a divisor.

Every proposed recurrence should state whether its dimension, polynomial degree, or number of blocks grows with `K` rather than `log K`. If a recurrence uses division modulo `N`, compute the GCD of each denominator before inversion: a nonunit can be useful output, not merely an implementation failure.

A particularly relevant new lead is Yann Tal's September 2026 preprint [Tal]. Its abstract describes modular-factorial improvements through Jacobi sums, reconstruction, and recursive reductions. **Its setting supplies a known prime modulus, or specified prime-power settings, and appropriate divisor conditions.** It is not a general unknown-factor-modulus algorithm. Inspect the full construction and ask which identities survive without knowing the prime or its group structure.

A useful adaptation must avoid assuming a factorization of `N`, its totient, auxiliary quantities of comparable unknown difficulty, or a known prime underlying `R_N`. Character sums are useful here only if the proposed way to construct or evaluate them is explicit and affordable.

**First experiment:** Implement a transparent block-product baseline, then one recurrence or invariant that attempts to answer only the unit predicate. Compare actual arithmetic counts and largest representation sizes.

### Track D — Hidden characteristic, lifted information, and annihilators

For `N = p*q`, the polynomial has different Frobenius structure in its two hidden components:

\[
F_N\bmod p=(1+X^p)^q-1-X^{pq},
\]

\[
F_N\bmod q=(1+X^q)^p-1-X^{pq}.
\]

The goal is not to perform these reductions using known factors. It is to find an operation over `R_N` that separates the two structures without knowing them.

Investigate Hasse derivatives, divided-power constructions, formal-series sections, finite differences, recurrence singularities, and annihilating ideals. Consider modest lifts such as modulus `N^t` where `t` is polynomially bounded, but derive what survives the lift and how it is extracted. Ordinary differentiation becomes zero after reduction modulo `N`; differentiation before reduction or in a lift retains other data, whose usefulness must be demonstrated rather than assumed.

A promising failure mode of a field algorithm is encountering a nonunit pivot. Instrument such algorithms to return `gcd(pivot,N)`. However, selecting a pivot guaranteed to become a useful nonunit is part of the research problem.

Avoid circular use of Lucas/Kummer digits in base `p`, a `p`-Frobenius map, or Cartier operators parameterized by an unknown factor. They can guide identities in the proof, but an operational algorithm must construct its own parameters from `N`.

**First experiment:** Build small lifted or sectioned representations and classify the denominators, pivots, and annihilators that arise. Separate structural factor discovery from accidental divisibility by a small tested integer.

### Track E — Transfer partial information into lattices or hidden relations

A sketch does not have to reveal a complete factor. Determine whether it yields a rigorously specified interval, congruence, or correlated set of bits that existing small-root methods can use. Coppersmith's work [Coppersmith] is the relevant reference, but extract the precise theorem and parameter requirements before applying it.

For hidden-relation approaches, use [Regev], [Pilatte], and [EG] to identify the information that classical postprocessing consumes. A short-relation existence theorem is not a procedure for finding those relations; a decoder is not a sampler.

One exact interface is multiplicative characters. With generators `a_i`, define

\[
\mathcal L=\{z\in\mathbb Z^d:\prod_i a_i^{z_i}=1\pmod N\}.
\]

If a computable character satisfies `chi(a_i)=exp(2*pi*i*w_i)`, then

\[
\langle w,z\rangle\in\mathbb Z\quad\text{for every }z\in\mathcal L.
\]

The question is whether explicit arithmetic supplies useful character phases or partial phase constraints efficiently. For squared generators `a_i=b_i^2`, ambient characters taking values only in `{+1,-1}` become trivial. A proposed source of phases must avoid that degeneracy and avoid constructing the unknown group decomposition first.

Also inspect [Gao] for the neighboring classical idea of deliberately avoiding uninformative shortest vectors. Transfer a concrete objective or construction, not merely the instruction to use lattice reduction.

Treat [Cleve] as a model-specific warning about generic oracle order finding, not as a lower bound for explicit integer arithmetic.

**First experiment:** Specify an exact partial-information format and demonstrate its downstream usefulness with a diagnostic oracle. Then separately attempt to produce it from `N` alone. Never combine those two stages into an unsupported end-to-end speed claim.

### Track F — One independent connection

Reserve some effort for a route not suggested by the polynomial lead. Examples include a different low-description-length generating function, a nonunit-producing algebraic identity, or an arithmetic invariant of a low-dimensional object.

The proposal must specify an executable operation and why it sees the hidden factors. Merely restating `pq=N` as optimization, SAT, a spectral problem, or a dynamical system is not enough. Show where the representation gives a cheaper operation or a smaller search problem.

## 4. Research procedure

### Phase 1: establish the mathematical and literature baseline

Reproduce the elementary lemmas and tests in this brief. Read the core sources needed for the two most promising tracks; do not read every listed paper before doing mathematics.

Create a compact map of the nearest prior art. For each relevant result record its exact problem, input representation, theorem assumptions, complexity variables, and the primitive it provides. Use primary sources. Verify title, authors, date, and theorem statement; earlier conversational summaries are leads, not evidence.

Keep discovery work separate from verification. In the first generation pass, propose at least six mechanisms, including one outside the supplied tracks. Then select two or three with the clearest operation and smallest missing lemma.

### Phase 2: write a candidate card before substantial coding

Use this template for each serious candidate:

```text
ID and one-sentence mechanism:
Target problem / supported input family:
Representation and how it is built from N:
Exact operation that extracts information:
Why that information exposes a factor:
Known lemmas being used:
New lemma required:
Correctness argument or conjecture:
Cost in n, including construction and extraction:
Randomness and success probability:
Closest prior method and the actual difference:
Smallest informative experiment:
Counterexample or scaling behavior that would refute this version:
Possible repair if that happens:
```

Reject or repair proposals that name the missing operation without implementing or bounding it. Do not abandon a promising representation merely because the first naive implementation expands it.

### Phase 3: alternate construction and criticism

For each selected mechanism:

1. Derive the strongest statement you can justify.
2. Implement the smallest exact experiment that distinguishes it from the obvious baseline.
3. Search for boundary cases and counterexamples.
4. Audit representation growth and total work.
5. Produce a revised mechanism or a narrow obstruction identifying what a revision must change.

Use a separate mathematical critic and complexity auditor if real subagents are available. Give the critic the candidate and code, not a persuasive narrative about expected success. Without subagents, perform separately documented review passes. Do not claim independent review unless it actually occurred.

The critic must state which assumption fails and give an example, derivation, or missing bound. “Factoring is believed hard” is not a technical critique.

Do not stop after a survey or after one failed proposal. Complete at least one repair cycle on the best candidate, unless a complete proof and independently checked implementation are already obtained.

### Phase 4: keep the result reproducible

Preserve proofs, candidate versions, failed constructions, seeds, raw measurements, and exact commands. A negative result about one representation is useful when its scope is explicit and it points to a different operation to try.

Allocate most effort to mathematical construction and discriminating experiments, not ever-expanding bibliographies or GPU optimization. Use exact integer/ring arithmetic initially. Floating-point exploration is allowed only with a later precision analysis.

## 5. Experimental design and leakage controls

### Input families

Begin with exhaustive small composites for algebraic correctness. Then use reproducibly generated balanced semiprimes at increasing bit lengths, for example `16, 24, 32, 48, 64`, extending only as resources permit. Successful factoring of small inputs is a correctness check, not evidence of polynomial scaling.

Include independent-prime balanced semiprimes, prime powers, composites with several distinct factors, highly unbalanced semiprimes, Carmichael numbers, and prime inputs. Near-equal primes, special algebraic forms, and deliberately smooth `p-1` or `q-1` belong in separately labeled diagnostic families.

### Separate algorithm inputs from diagnostic truth

The candidate receives only `N`, public configuration, and a random seed. A separate test harness may retain factors to compute labels and understand failures.

Never leak factors through filenames, seeds deterministically encoding factors, candidate-specific parameter choices, precomputed group orders, smoothness-selected instances, library factorization calls, or a simulator initialized with the factors.

Use actual random generation of the intended family and record the sampling procedure. Freeze candidate choices before the holdout evaluation. Special test families are useful but must not silently replace the intended distribution.

### Record meaningful costs

For each run record at least:

- Correctness, result, seed, wall time, timeout, and peak memory.
- Ring operations and GCD calls, with operand bit lengths.
- Largest polynomial degree actually materialized and number of stored coefficients.
- Circuit nodes, recurrence dimension, matrix size, section depth, and rank where relevant.
- Number of sketches, candidate combinations, restarts, and rejected trials.
- Numerical precision and error bounds for any approximate component.

Separate preprocessing, information acquisition, decoding, and verification. Count all work, including unsuccessful attempts. Distinguish parallel depth from total classical work.

Compare against transparent trial division and Pollard-rho implementations on small tests, a relevant product/block baseline for Track C, and the nearest prior construction for the actual proposed operation. Use controlled same-language arithmetic comparisons where possible. Do not claim superiority to GNFS from a small-input benchmark or an improvement over trial division.

Empirical fits are diagnostics, not complexity proofs. Show raw counters versus both `n` and the actual structural parameters. In particular, a method polynomial in a hidden factor can look innocuous on tiny bit lengths.

## 6. Mandatory complexity and correctness audit

For every proposed breakthrough, answer the following before calling it an algorithmic improvement.

**Representation.** Is the object genuinely stored in polynomial space, or merely described by a short program that later expands? Are all required operations efficient on that representation?

**Coefficient access.** Does the method require sequential coefficients up to `p`, a dense truncation at `K`, or exponentially many terms after sectioning? Section 2.3 rules out individual coefficient evaluation as the missing advantage here.

**Field assumptions.** Does a theorem require a known prime, an algebraically closed field, a characteristic bound, or invertible denominators? Prove the needed ring version or report factor-revealing nonunits.

**Hidden inputs.** Is any step using `p`, `q`, `phi(N)`, `lambda(N)`, an unknown order, a hidden lattice basis, or informative samples supplied by an oracle? State exactly how it is obtained.

**Numerical cost.** An error tolerance `2^(-poly(n))` uses polynomially many bits. An error tolerance `2^(-p)` needs order-`p` bits. A tiny error written with a short exponent is not automatically cheap to achieve. Conversely, do not call every exponentially small error an exponential-precision requirement.

**State-space cost.** A short diffusion rule, tensor description, or matrix-power expression is not enough. Bound the represented state, ranks, contractions, and the cost of the final observable.

**Degree and root bounds.** Check polynomial identity-test assumptions when formal degree is exponential. Over a finite ring, nonzero polynomials can vanish at every scalar point. Auxiliary prime moduli also change coefficient divisibility: reducing modulo an unrelated prime does not preserve vanishing modulo `N` for free.

**Randomness.** Bound the success probability and cost per attempt, including constructing the distribution and rejecting bad outputs. Distinguish worst-case, average-case, and restricted-family claims.

**Certificates.** A candidate factor is cheaply verifiable. That does not bound search time or certify an average-case/polynomial-time theorem. A decoder's correctness with ideal samples does not validate the source of those samples.

**Novelty.** Compare with the closest specific method. A new notation, a different implementation of a known identity, or a known reduction rediscovered here is not an asymptotic breakthrough. It can still be useful scaffolding for the next modification.

Use these checks to improve proposals, not to replace proposal generation.

## 7. Workspace and required deliverables

Suggested layout:

```text
factoring-connection/
  README.md
  STATUS.md
  references.md
  notes/
    seed_lemmas.md
    prior_art_map.md
    candidate_001.md
    candidate_002.md
    review_and_repairs.md
  src/
    seed_checks.py
    ring_polynomials.py
    circuit_sections.py
    quotient_sketches.py
    candidates/
  tests/
  experiments/
    configs/
    raw/
    analysis/
  REPORT.md
  NEXT_STEPS.md
```

Use Python with exact arithmetic for the first prototypes; use existing computer-algebra packages where appropriate, but inspect any operation that may internally factor a modulus. Optimize only after identifying an informative mechanism. Do not buy cloud resources, publish results, or contact authors without authorization.

A candidate interface should accept `N` and public configuration only, return a factor or explicit failure, and emit an auditable metrics record. Reference factorization functions belong in the harness, never inside the purported new factorizer.

The final report must contain:

1. **Best surviving mechanism:** its precise new operation and the strongest supported claim.
2. **Mathematics:** proved lemmas, conjectures, counterexamples, and any remaining proof obligation, clearly separated.
3. **Experiments:** reproducible commands, raw results, input distributions, failures, and scaling diagnostics.
4. **Comparison:** closest prior art and exactly what changed.
5. **Continuation:** the next most decisive mathematical question or experiment, with a concrete proposed attack.

Use status labels such as `proved`, `experimentally supported`, `conjectured`, `refuted for this model`, and `not yet tested`. A useful report may contain a novel mechanism plus a sharply identified unresolved lemma. Do not inflate that into a completed factorizer, and do not stop at “no breakthrough found” when the work supports a more specific next step.

## 8. Primary-source reading list

These are selected starting points, not a claim of an exhaustive literature search. Bibliographic/abstract pages were checked when preparing this brief; the AKS seed identity was also checked in the paper. Except for the elementary derivations above, do not assume that every cited theorem or implementation has been independently verified here. Open the relevant full papers, record the version consulted, and verify exact conditions before relying on them.

**[AKS]** Manindra Agrawal, Neeraj Kayal, Nitin Saxena. *PRIMES is in P*. Annals of Mathematics 160 (2004), 781–793. [Paper](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/primality_journal.pdf). Focus on Lemma 2.1 and the quotient-ring reduction in Section 2; distinguish primality detection from factor extraction.

**[BGS]** Alin Bostan, Pierrick Gaudry, Éric Schost. *Linear Recurrences with Polynomial Coefficients and Application to Integer Factorization and Cartier–Manin Operator*. SIAM Journal on Computing 36(6) (2007), 1777–1806. [Author-hosted paper](https://mathexp.eu/bostan/publications/BoGaSc07.pdf). This directly connects distant recurrence terms, product computations, and classical factorization. Audit the dependence on the recurrence index.

**[Harvey]** David Harvey. *An exponent one-fifth algorithm for deterministic integer factorisation*. arXiv:2010.05450 (2020). [Record](https://arxiv.org/abs/2010.05450). A rigorous `N^(1/5+o(1))` classical factoring construction; a useful structured-search baseline, not polynomial in input bit length.

**[HH]** David Harvey, Markus Hittmeir. *A log-log speedup for exponent one-fifth deterministic integer factorisation*. arXiv:2105.11105 (2021). [Record](https://arxiv.org/abs/2105.11105). Study the arithmetic structure and the actual batching operation, not just its asymptotic headline.

**[Sparse]** Joris van der Hoeven, Grégoire Lecerf. *Sparse polynomial interpolation: Faster strategies over finite fields*. Author-hosted final version dated 2 May 2024. [Text](https://www.texmacs.org/joris/ffsparse/ffsparse.html). Relevant to projections, support reconstruction, and circuit access. Its field, sparsity, characteristic, and root-of-unity assumptions require explicit checking for our setting.

**[Tal]** Yann Tal. *Computing Modular Factorials Below the Square-Root Barrier*. ECCC TR26-211, September 2026 preprint. [Record and paper access](https://eccc.weizmann.ac.il/report/2026/211/). A recent lead involving known prime/prime-power moduli, Jacobi sums, and recursive factorial computations. Treat the claims as preprint results to inspect, and preserve the known-modulus/divisor conditions. No general unknown-factor-modulus speedup is asserted by this brief.

**[Coppersmith]** Don Coppersmith. *Small Solutions to Polynomial Equations, and Low Exponent RSA Vulnerabilities*. Journal of Cryptology 10 (1997), 233–260. [Publisher record](https://link.springer.com/article/10.1007/s001459900030). Use to specify exactly what partial factor information would suffice, rather than casually assuming that arbitrary correlated bits are usable.

**[Regev]** Oded Regev. *An Efficient Quantum Factoring Algorithm*. arXiv:2308.06572 (2023; consult current version). [Record](https://arxiv.org/abs/2308.06572). A source for the hidden-relation formulation and the separation between quantum information acquisition and classical postprocessing.

**[Pilatte]** Cédric Pilatte. *Unconditional correctness of recent quantum algorithms for factoring and computing discrete logarithms*. arXiv:2404.16450; Forum of Mathematics, Pi 14 (2026), e5. [Record](https://arxiv.org/abs/2404.16450). Structural existence/correctness results for suitable modified constructions; do not equate existence of short relations with constructive classical access to them.

**[EG]** Martin Ekerå, Joel Gärtner. *Extending Regev's factoring algorithm to compute discrete logarithms*. arXiv:2311.05545 (2023; consult current version). [Record](https://arxiv.org/abs/2311.05545). Relevant to what approximate lattice information suffices for classical recovery. Extract the exact noise, distribution, and precision requirements.

**[Gao]** Yiming Gao, Yansong Feng, Honggang Hu, Yanbin Pan. *On Factoring and Power Divisor Problems via Rank-3 Lattices and the Second Vector*. arXiv:2512.19076 (December 2025 preprint). [Record](https://arxiv.org/abs/2512.19076). A concrete neighboring classical lead emphasizing useful rather than merely shortest vectors. Reproduce and inspect the construction before generalizing it.

**[Cleve]** Richard Cleve. *The query complexity of order-finding*. arXiv:quant-ph/9911124 (1999). [Record](https://arxiv.org/abs/quant-ph/9911124). Read the oracle model carefully. Its lower-bound conclusions do not automatically apply to algorithms exploiting the explicit representation of modular integers.

For additional searches, combine terms such as “arithmetic circuit coefficient sections,” “order of vanishing straight-line programs,” “polynomial identity testing composite characteristic,” “Hadamard product compressed polynomials,” “factorial unit test composite modulus,” “holonomic recurrence distant terms,” and “constructive AKS factor extraction.” Treat these as search directions, not as claims that a particular theorem already exists.

## 9. Small exact smoke test

Save this as `src/seed_checks.py`. It is reference/testing code using trial division and explicit coefficients. It is deliberately not a fast factoring algorithm and must not appear as the hidden implementation of a candidate.

```python
"""Exact tests for the elementary seed lemmas; small instances only."""
from math import comb, gcd, isqrt


def least_divisor_reference(value: int) -> int:
    if value < 2:
        raise ValueError("value must be at least 2")
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return divisor
    return value


def check_seed_lemmas(limit: int = 200) -> int:
    if not 4 <= limit <= 2000:
        raise ValueError("use a small test limit between 4 and 2000")
    composites = 0
    for modulus in range(2, limit + 1):
        first = least_divisor_reference(modulus)
        # Only interior coefficients: the endpoints of F_N are removed.
        coeffs = [comb(modulus, k) % modulus
                  for k in range(1, modulus)]
        if first == modulus:
            assert not any(coeffs), (modulus, "prime identity")
            continue

        composites += 1
        assert not any(coeffs[:first - 1])
        assert coeffs[first - 1] == modulus // first
        for k, coefficient in enumerate(coeffs, start=1):
            assert k * coefficient % modulus == 0
            if coefficient:
                assert 1 < gcd(k, modulus) < modulus

        # product = (K-1)! mod modulus at each iteration.
        product = 1
        for cutoff in range(1, isqrt(modulus) + 2):
            prefix_zero = not any(coeffs[:cutoff - 1])
            assert prefix_zero == (cutoff <= first)
            assert prefix_zero == (gcd(product, modulus) == 1)
            product = product * cutoff % modulus

    for point in range(561):
        assert (pow(1 + point, 561, 561)
                - 1 - pow(point, 561, 561)) % 561 == 0
    return composites


if __name__ == "__main__":
    count = check_seed_lemmas()
    print(f"Seed checks passed for {count} small composites.")
```

During preparation, the first-index/value and coefficient-annihilator identities were also checked directly for all 831 composites up to 1,000, and the scalar-zero property of `F_561` was checked at all 561 residues. These checks supplement the proofs; they are not evidence of a fast extraction procedure.

## 10. Start here

Read Sections 1–2, run the smoke tests, and create `STATUS.md`. Then investigate binary sectioning and weighted quotient sketches in parallel or as separate passes. Produce explicit candidate cards, implement the most discriminating small experiment for each, and perform a repair cycle after criticism.

The aim is to find the connection, not merely describe why existing approaches are expensive. Keep asking what information a representation exposes cheaply that was not cheaply accessible before, and make that operation concrete.
