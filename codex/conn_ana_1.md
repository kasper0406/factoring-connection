# Blind brainstorm: unexpected connections between analysis/physics/probability fields and integer factorization

You are GPT-6.1-sol, asked for an independent, creative but honest survey. Do NOT read files under notes/ or the
project's results in the repository (this is a blind round). You may use web search for literature.

Goal: find unexpected connections from the fields below to integer factorization that could lead to an
asymptotically faster classical algorithm: beating heuristic L[1/3, 1.902] (MNFS), the rigorous L[1/2] /
L[2/5] records, or the deterministic N^{1/5} (Harvey–Hittmeir) record.

Background barriers you should respect (established by the project; treat as known):
- Every known "leak" that makes p and q behave differently is one of three families: (M) a group order /
  smoothness (p−1, ECM, class numbers), (A) an archimedean approximation (Fermat, Lehman, Coppersmith-type known
  bits; exponential cost), (Q) a quadratic genus symbol (Jacobi symbols, which reveal nothing beyond N mod m).
- Theorem PT: algorithms using N only as a modulus (ring operations mod N, Jacobi symbols) learn nothing about
  the factors except by accidentally splitting N (gcd).
- Anything depending only on N mod m for fixed m is an aggregate and useless. Many invariants (Dedekind sums
  s(a,N), Salié sums mod N, #X_0(N)(F_{ℓ²}), genus of X_0(N), cat-map orbit counts mod N, #E(Z/N)) are
  factoring-equivalent: computing them is as hard as factoring.
- Size–entropy law: any "direct target" that would give α < 1/3 is polynomial-time factoring in disguise.
- Constructing smooth numbers near multiples of N costs L[2/3] per output for every constructor class tried.

Fields (go beyond the seed ideas):
1. Harmonic analysis: restriction, decoupling (Bourgain–Demeter–Guth; Vinogradov mean value), Kakeya (incl. a
   claimed Oct 2026 resolution by OpenAI, unverified), exponential sums over smooth numbers (Harper 2016).
2. Random matrix theory: L-function zeros, GUE; spectral statistics of operators built from N.
3. Spectral theory, trace formulas, quantum chaos: Selberg trace formula on Γ_0(N); quantum cat maps mod N;
   arithmetic QUE.
4. Probability: random walks on (Z/N)* and SL₂(Z/N) (mixing, cutoff), percolation, coupling from the past.
5. Statistical mechanics / integrable systems: Bethe ansatz, Yang–Baxter, transfer matrices with period N,
   partition functions whose zeros/periodicities depend on the factorization; Lee–Yang.
6. QFT and mirror symmetry: Chern–Simons at level N, topological recursion, modularity in physics.
7. Convex geometry and optimisation: SOS certificates for factoring, SDP/LP hierarchies, interior-point on
   lattice problems, polytope relaxations of x·y = N.
8. Information and coding theory: list decoding (CRT codes, Guruswami–Sudan), compressed descriptions, entropy.
9. Dynamical systems and fractals: arithmetic dynamics mod N (orbits, periods; Pollard rho), Collatz-type maps,
   Julia sets over Z/N.
10. PDE / numerical analysis analogies: heat-kernel or theta-function methods for counting divisors; multigrid
   or FFT accelerations of sieving.

For EACH field give 3–6 concrete candidate connections. For each:
 (a) the object built from N alone and the computable feature that would carry factor information, or the
     algorithmic primitive that would replace sieving / linear algebra;
 (b) literature (cite real papers; say UNVERIFIED if unsure);
 (c) a verdict: RESTATEMENT (needs a factor to write down) / AGGREGATE (depends only on N mod m) /
     FACTORING-EQUIVALENT / REDUCES-TO-M, -A or -Q / KNOWN (already used in factoring) / OPEN-WORTH-A-LOOK;
 (d) 2–4 sentences on exactly WHERE it fails (which step needs the factors, costs size N, or lands in a family).
Prioritise surprising connections. Be honest: most will fail; say so precisely.

Finish with your ranked top 3 OPEN-WORTH-A-LOOK items (if any survive), each with a precise first lemma or
question that a researcher could attack in a week. Also give a probability estimate that each leads to any
asymptotic improvement.

Write the final answer as Markdown, ≤ 3500 words.
