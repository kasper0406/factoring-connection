# hidclass: deciding the hidden factor's class in a fixed abelian field (Briefs 1 + 3, merged)

Repo is READ-ONLY for you; work and scratch files go in your cwd. Relevant repo notes (read as needed):
notes/coef_results.md, notes/kern_results.md (u_2 ⇒ generalized QR in a norm-one torus; trace forms lose it),
notes/genusleak_results.md, notes/newfamily_results.md (reciprocity invariants are aggregates or ≡FACT),
notes/leakmap_results.md (leak families M/A/Q), notes/factorbits_results.md (Theorem PT),
notes/entropy_results.md (CF(√7N) leaks a character of p, but decaying like N^{-1/4}), BARRIERS.md.

## Background (established in the orchestrator session, 2026-11-26; spot-check, don't redo)
- For a weight-1 dihedral EIGENFORM of 2-power level, a(N²) ≡ 1 mod 𝔓 for odd N (ray class group is a 2-group, so
  ψ ≡ 1 mod 𝔓; Tate: ρ̄ trivial). Exact: a_ψ(N²) = ψ((N)) ∏_{p|N split} S_{e_p}(u_p), with S_e(u) = Σ_{|j|≤e} u^j and
  u_p = ψ(𝔭/𝔭̄). Non-eigen integral combinations (theta series of lattice cosets) need eigenforms beyond 𝔓, and
  their mod-2 coefficient at N² depends on non-congruence data such as ∏_{p|N}(2/p)_4. Checked: for N with all
  p ≡ 1 mod 8, ∏(2/p)_4 = 1 ⟺ N is primitively x²+64y² (520 N's).
- Character PRODUCTS χ(N) = ∏ χ(p)^{e_p} are free; the class of a single hidden prime (or prime-class counts mod 2)
  is not determined by N mod anything.

## The two targets (the briefs' formulas were lost in transit; these are reconstructions)
- **T1 (Brief 1).** N = rq with r ≡ 1 mod 8 and q ≡ 3 mod 8. Bit b(N) = [r ≡ 9 mod 16]. Verified numerically
  (600 semiprimes): with A(z) = Σ_{n≥0} z^{n²}, b(N) = [z^{N²}] A(z)(A(z^{128}) + A(z^{256})) mod 2. The brief also
  had "U(N) counts integer triples with [a ternary form] = [N² or similar], satisfies an exact linear identity in r,
  and U(N) mod 2^j ⟺ r mod 2^{j+2}". Reconstruct the most natural U (e.g. r_3(N²) = 6∏(σ(p^e) − (−1/p)σ(p^{e−1}))
  gives 6(N+2r) for r ≡ 1, q ≡ 3 mod 4; find the variant matching "mod 2^j ⟺ r mod 2^{j+2}") and state it.
- **T2 (Brief 3).** N = pq ≡ 1 mod 7. Decide (p,q) ≡ ±(1,1) mod 7 (both split completely in the cubic field K of
  conductor 7) versus {p,q} ≡ ±{2,4} (both inert). Cubic: F(u,v) = u³+u²v−2uv²−v³ (min. poly of 2cos(2π/7);
  automorphism σ: x ↦ x²−2). Split ⇒ 9 roots mod N and a nonzero solution of F ≡ 0 mod N with max(|u|,|v|) ≲ √N;
  inert ⇒ N | u and N | v. Equivalently: is N a norm from K / is (Z/N)[x]/F totally split?
Both are instances of **generalized Blum-integer recognition**: given N = pq, decide the (unordered) Frobenius class
pair of p, q in a fixed abelian extension, given their product class.

## Questions
H1. **Literature and status.** Blum-integer recognition (p ≡ q ≡ 3 vs ≡ 1 mod 4 for N ≡ 1 mod 4; ⟺ "is N a sum of
    two squares" for semiprimes), and its generalisations to p mod m and to norm forms of abelian fields. Known
    algorithms, interactive/ZK proofs (van de Graaf–Peralta, Gennaro–Micciancio–Rabin, Camenisch–Michels, ...),
    complexity-class placement, any reduction to/from QR, DCR, factoring. Search the web. Is anything known to make
    T1 or T2 easy?
H2. **Relations between targets (self-reduction / lifting).** Model an oracle O_m(N) returning the unordered class
    pair {p mod m, q mod m} (allow it on all N, including 3-prime N·s with known s if useful). Prove or refute:
    (a) O_16 ⇒ O_32 (Brief 1 "breakthrough" lift, j → j+1); (b) O_4 ⇒ O_8 or O_16; (c) O_7-type (T2) ⇔ O_4-type;
    (d) O_{2^j} for 2^j ≥ N^{1/4+ε} ⇒ factoring (Coppersmith — confirm constants). Note that multiplying N by known
    primes only adds known classes, so per-N information looks bounded; make this into a statement if true.
H3. **Reciprocity / descent obstruction (Brief 1 direction 1–2).** Formalise a "reciprocity-descent" model
    (Euclid-like index reductions whose state updates come from theta multipliers, Gauss sums, Jacobi symbols,
    Hecke/Atkin–Lehner at 2-power level) and prove it computes only character products of N, hence never b(N);
    or find the escape. Specifically: can SL₂(Z) transformations of A(z)A(z^{128}) track the single coefficient at
    N² mod 2 with a growing state (as in the Jacobi algorithm) without expanding Θ(N) coefficients? Where exactly does
    it break? Also check Brief 1 direction 3: does the semiprime promise give a functional identity for U(N) mod 2^j
    that fails for general composites, and is it exploitable?
H4. **T2 root-free decision.** (a) A = (Z/N)[x]/F has a computable ring automorphism σ; in the inert case
    Frobenius equals σ^{±1} on each component, in the split case it is trivial. Is there a Frobenius-free distinguisher
    (power maps, norms, σ-twisted equations a^k = σ(a)^m, unit-group exponent tests, sampling) — or prove each such
    test needs a group-order leak (family M) or a root? Compare kern (quadratic analogue reduced to QR in a torus).
    (b) Show a root mod N factors N with constant probability given a random-root oracle (σ acts diagonally, mixed
    roots give gcds); then is the DECISION "∃ solution with max(|u|,|v|) ≤ N^{0.9}" any easier than finding one?
    (c) Theta/Epstein-sum decision over the solution set: is any such sum computable without the lattice?
H5. **Experiments (pre-register first: write prereg.md in your cwd with hypotheses, statistics, sample sizes and
    the pass threshold BEFORE running anything).** 40–80-bit balanced semiprimes, ground truth from the factors, flag
    any gcd split. Test the 2–4 most promising Frobenius-free statistics for T1 and T2 (e.g. sampling statistics in
    (Z/N)[x]/F using σ; Jacobi-type symbols in Z[i], Z[ζ_8], Z[ζ_16] modulo N; CF(√(kN)) features — expect the
    N^{-1/4} decay seen in entropy, so report the decay exponent). Report the advantage with confidence intervals and
    how it scales with bit size.
H6. **Verdict.** What is proved / known / refuted / measured; the sharpest remaining open statement for T1 and T2;
    probability that either has a poly(log N) algorithm; and whether anything bears on factoring below L[1/3].

≤ 2500 words. Label every claim proved / known (with citation) / heuristic / measured / speculative.
