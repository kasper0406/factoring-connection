# Blind round: unexpected connections to integer factorization (fields: complexity, NCG, quantum, crypto, ML, OT, numerics, sparse recovery, games, arithmetic statistics)

You are a research mathematician / theoretical computer scientist. Goal: find **unexpected connections** between the
fields below and **classical integer factorization**, that could plausibly lead to an asymptotically faster
algorithm than the heuristic NFS L_N[1/3, 1.902] (MNFS), the rigorous L_N[1/2] (Lenstra–Pomerance; this project
claims an internal rigorous L[2/5]), or the deterministic N^{1/5} (Harvey 2021, Harvey–Hittmeir 2022).

Work blind: do not read the repo (other than files you create). Use web search for literature checks.

Known barriers to respect (from earlier work on this project; treat as given):
- Theorem PT: methods that use N only as a modulus (ring operations mod N, polynomials over Z/N) learn nothing about
  p, q except by splitting N, i.e. through a gcd with something whose order/structure differs between p and q.
- Every known "leak" that treats p and q differently is one of three families: (M) a group order or smoothness
  (p−1, #E(F_p), class numbers), (A) an archimedean approximation (Coppersmith, Fermat; exponential), (Q) a
  quadratic genus character (Jacobi symbols). A fourth family would be big news.
- Quantities depending only on N mod m for fixed m are AGGREGATE (useless). Quantities like φ(N), σ(N), #X_0(N)(F_ℓ²),
  Dedekind sums s(a,N), Salié sums mod N are FACTORING-EQUIVALENT (no gain unless they have an independent fast
  algorithm).
- Size–entropy law: every "direct target" that would give exponent < 1/3 is poly-time factoring in disguise; building
  smooth numbers near multiples of N costs L[2/3] per output for every constructor class tried.
- Classical Fourier side of Shor: an observable g gives period-concentrated Fourier mass of g(a^x) only if g
  predicts a large-order multiplicative character; additive, ring-op and computable-character observables are flat.

Fields (aim for 3–6 concrete candidate connections per field; go beyond the seeds):
1. Algebraic complexity: VP vs VNP, GCT, complexity of n! mod N (Strassen; Bostan–Gaudry–Schost; any 2025–26
   "sub-square-root modular factorial" results, e.g. ECCC TR26-211 if it exists), polynomial identity testing,
   tensor rank of matrix multiplication / of multiplication mod N.
2. Operator algebras / noncommutative geometry: Bost–Connes system (partition function ζ, KMS states), Connes' RH
   programme, C*-algebras of ax+b groups over Z/N, crossed products.
3. Quantum information, classical side: classical simulation of parts of Shor (tensor networks, stabiliser/Clifford,
   matchgates, low-entanglement period finding), Tang-style dequantisation of period finding / QFT.
4. Cryptographic constructions: isogeny graphs and class-group actions over Z/N (CSIDH-like over composite moduli),
   multilinear maps / iO and their attacks, hidden-order groups (RSA groups, class groups) and known
   reductions to factoring.
5. Machine learning: learning factor bits from N; learned sieving/polynomial selection (constant factor only?);
   transformers on arithmetic. Be honest about asymptotics.
6. Optimal transport: couplings between smooth-number distributions, Wasserstein arguments for relation supply.
7. Numerical analysis: FMM, FFT, hierarchical matrices for sieving, batch smoothness (Bernstein), sparse GF(2)
   linear algebra (Wiedemann / block Lanczos are ~n^2 in the matrix dimension; is sub-n^2 dependency finding possible,
   and what would it buy in the NFS constant?).
8. Compressed sensing / sparse recovery / inverse problems: recovering sparse exponent vectors; group testing for
   smoothness; sparse FFT.
9. Game theory / interactive proofs: verifier–prover structure for factoring-related statements (IP, PCP, zero
   knowledge of N's structure, e.g. proofs that N is a product of two primes / Blum integer).
10. Arithmetic statistics: Cohen–Lenstra, Malle, Bhargava parametrisations (2-torsion / 3-torsion in class groups of
   Q(√−N), binary cubic/quartic forms), counting number fields; does any parametrisation *construct* factor-revealing
   objects (ambiguous forms) more cheaply than relations?

For EACH candidate:
- What object is built from N alone; what computable feature would carry factor information, or what primitive would
  replace sieving / linear algebra.
- Literature status (cite).
- Verdict: RESTATEMENT (needs a factor to write down) / AGGREGATE (depends only on N mod m) / FACTORING-EQUIVALENT
  (or DL-equivalent) / REDUCES-TO-<M, A, Q or a named barrier> / KNOWN (already used in factoring, with reference)
  / OPEN-WORTH-A-LOOK. 2–4 sentences on exactly WHERE it fails.
- Flag anything you think is surprising.

Then: rank your top 3 OPEN-WORTH-A-LOOK items (if any) with a precise first question or lemma each. Also give, for
the numerical-analysis / linear-algebra item, a concrete quantitative answer: if the GF(2) kernel step cost
n^{ω'} with ω' < 2 (n = matrix dimension), what does the NFS constant become (relation/linear-algebra balance), and
is that known in the literature? Be quantitative, short, and honest; no padding. Output markdown.
