You are a research mathematician / computational number theorist. Blind round for agent `nativemod` (phase 29) of a
research project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you;
scratch files go in your current directory; python with sympy/gmpy2/python-flint/cypari2 (PARI) is at
/Volumes/git/factoring-connection/.venv/bin/python). Rigour and honesty matter more than a positive result. Use web
search to check every literature claim; give exact citations (authors, title, venue, year, theorem/section number)
and mark anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (proof given
here) / HEURISTIC / CONJECTURED / OPEN / FALSE. Negative literature findings must be stated as such.

Optional background in the repo: notes/nonsmooth_results.md §5.2 (Prop FROB: producing the Frobenius of an étale
Z/N-algebra factors N), notes/tracedata_results.md (a_N = a_p a_q for a fixed elliptic curve; whether a_N alone
factors N is open there), notes/prim_results.md and notes/leakmap_results.md (the "PRIM" model: straight-line
programs over Z/N with N-dependent exponent gates x -> x^{f(N)}), notes/conn_alg.md (Kronecker–Hurwitz identities
Σ_t H(4N − t²) = 2(N + q) etc. for N = pq).

Setting. N = pq (p < q primes, balanced, unknown). The question: a reduction "computing X factors N" is a new view
of factoring, not a hardness proof. Is there a NATIVE algorithm (one that never factors N) for X when X is a
Fourier coefficient or Hecke trace of a FIXED-level modular form at composite index N? Cost is measured in
n = log2 N; "fast" = poly(n), "interesting" = anything below L_N[1/3] or below the best known factoring cost.

QUESTIONS
1. Edixhoven–Couveignes (and Bruin, Zeng, Peng Tian) compute a(p) for prime p in poly(log p) via mod-ℓ Galois
   representations: compute the ℓ-torsion field / projective polynomial P_ℓ over Q (fixed level, weight, ℓ), reduce
   mod p, compute the Frobenius x -> x^p on F_p[x]/(P_ℓ), read off tr ρ_ℓ(Frob_p) mod ℓ, CRT over ℓ ≤ O(log p).
   At composite N: (a) What does each step become over Z/N: the algebra A = (Z/N)[x]/(P_ℓ), the map x -> x^N on A,
   and its "trace"? (b) x -> x^N on A is a ring endomorphism of A ⊗ F_p and A ⊗ F_q separately equal to
   Frob_p^{q}·... — compute precisely what operator it is on each CRT component (Frob_p^{?} on A⊗F_p ...) and what
   any natural trace/character of it equals in terms of ρ_ℓ(Frob_p), ρ_ℓ(Frob_q). Is it related to a(N) = a_p a_q
   (for a newform with a_N multiplicative), to tr ρ(Frob_p^q) ⊗ ..., or to something else? (c) Can it be computed
   in poly(n) from N alone, and does it factor N? (d) Is this a "PRIM program" (straight-line over Z/N with
   N-dependent exponent gates), or an instance of "produce Frobenius of an étale Z/N-algebra" (factoring-hard), or
   something else?
2. Eisenstein series: r_4(N) = 8σ(N) for odd N (coefficient of θ^4, level 4, weight 2). Is there ANY known native
   algorithm for σ(N), or for the N-th coefficient of a fixed-level Eisenstein series, that does not factor N?
   Consider: modular symbols at fixed level, Hecke operators on a fixed finite-dimensional space, q-expansion
   recurrences/modular equations, theta-function evaluation, lattice-point counting in a 4-ball of radius √N.
   State costs. (Note: σ(N) for N = pq gives p + q, so it factors N; the question is native cost.)
3. Fixed-level Hecke action. T_N acts on the FIXED finite-dimensional S_k(Γ_0(M)). Computing its matrix by the
   definition or via Heilbronn/Merel matrices costs ~σ(N) or ~N log N. Hecke multiplicativity
   T_m T_n = Σ_{d|(m,n)} d^{k−1} T_{mn/d²} makes T_N a polynomial in T_p, T_q. Is there any route to T_N (or to
   its trace) that treats N "as an integer" (binary expansion, repeated doubling, Manin's continued-fraction trick,
   Merel's set for T_N of size ~σ(N)...) without factoring? Make precise why every known route uses the
   factorisation. Compare with Bach–Charles (arXiv:0708.1192, "The hardness of computing an eigenform"): exactly
   which ingredients does their reduction use (a(N) and a(N²)? the Ramanujan bound? a specific form like Δ?), and
   what is OPEN (a(N) alone; a(N) for a weight-2 form whose a_N = a_p a_q).
   Also: Eichler–Selberg Tr T_N on S_k(SL_2(Z)) contains −½ Σ_{dd'=N} min(d,d')^{k−1}; does the trace of T_N on a
   fixed space factor N (for k ≥ 4 where the space is nonzero, and for dimension-0 spaces where the trace formula
   gives an identity among class numbers)? Make the class-number identities explicit for small k (k = 2, 4, 6, 8,
   10, 14 where S_k(1) = 0) and state what they give.
4. Brandt matrices / quaternary forms: counting x in a fixed maximal order of a definite quaternion algebra (or a
   fixed positive quaternary lattice) with nrd(x) = N. Is there any lattice-point counting algorithm for a fixed
   4-dimensional ellipsoid of radius √N faster than ~N (or ~N^{c} for some c < 1) that does not factor N?
   (Gauss-circle-type algorithms: for 2-dim circles the count up to X is computable in ~X^{1/3}; for the shell
   count at exactly N this is a difference of two such counts — what is the cost in dimension 4?) Barvinok is for
   polytopes. Do any of these give a sub-exponential-in-n native algorithm?
5. Verdict: is "compute the N-th coefficient of a fixed modular form in poly(log N) time" (a) ruled out by some
   model barrier, (b) equivalent to a known open problem, (c) genuinely open? State the single sharpest first
   question. Is the natural generalisation (level and weight part of the input) NP-hard or #P-hard? (Expected not;
   check.)

Deliverable: a self-contained answer (≤ 4000 words) with formulas and proofs where you claim PROVED. Include any
small PARI/python checks you ran (commands and outputs) in an appendix.
