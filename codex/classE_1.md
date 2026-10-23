# classE_1: BLIND round — continued-fraction functionals of a(N)/N at structured arguments

You are a research mathematician (analytic number theory, continued fractions / Euclidean algorithms, lattice-point
geometry, complexity). The repo /Volumes/git/factoring-connection is READ-ONLY context; write scratch files and any
code only under the current directory (python3 with sympy/gmpy2 is available via
/Volumes/git/factoring-connection/.venv/bin/python). Do NOT read notes/classE_* (you are blind to our own work).
You MAY read notes/newfamily_results.md (§2.3 Prop DS/DSF, §2.5 class E), notes/leakmap_results.md (§4 Lemma GAR,
§5 class H) and notes/genusleak_results.md (§6.2 "open class").

## Setting
N = pq balanced (p < q < 2p, both with the top two bits set), factors unknown. "Floor": anything that is a function
of N mod m (m fixed) plus the real size of N carries no factor information beyond the fiber floor. Targets:
ε = (−1/p) on N ≡ 1 (mod 4), and the archimedean bit [p/q ≥ 7/8].

**Class E.** For a = a(N) computable from N with gcd(a, N) = 1, let a/N = [0; a_1, …, a_r] with convergent matrices
γ_k ∈ SL_2^±(Z). Class E = additive functionals S = Σ_k F(a_k, γ_k mod m₀) with fixed m₀ and bounded-complexity F
(polynomial-size output). It contains Dedekind sums (12 s(a,N) = Σ(−1)^{i+1}a_i + (a+ā)/N − 2 + (−1)^r), the
Rademacher Φ, level-M modular symbols via Manin's trick, Euclid-remainder residue counts, quantum-modular cocycles.
Known: Prop DSF — the "local part" s(a,p)+s(a,q) of s(a,N) = s(a,p)+s(a,q)+U(a,N) factors N. A battery of 84
statistics at a ∈ {⌊√N⌋, ⌊αN⌋ (golden), ⌊N^{2/3}⌋, (⌊√N⌋+1)²−N, 2^{(N−1)/2} mod N, hash} found no signal.

## Questions (answer each with proofs where possible; label proved / heuristic / conjectured / unverified)

Q1 (structure). Write the CF / Euclid path of a/N in terms of public data and (where relevant) p, q for:
 (a) a = ⌊√N⌋; (b) a = ⌊N^{2/3}⌋; (c) a = ⌊N/k⌋ ± c, small k, c; (d) a = ⌊ωN⌋ for fixed real ω;
 (e) a = 2^{(N−1)/2} mod N and a = 2^{N−1} mod N (CRT combinations).
 In particular: for a = ⌊√N⌋ is the CF of a/N related exactly to the CF of √N (or √(4N)) and/or of p/q? Derive it.
 Is there any a(N) computable from N whose CF path contains a convergent/remainder carrying factor information
 (e.g. a convergent denominator divisible by p, or ≈ p, or a convergent ≈ p/q)? Characterise exactly when a convergent
 denominator of a/N is divisible by p (we suspect: iff a ≡ ±r/s mod q with rs ≲ q/p — prove or refute) and compare
 with Coppersmith's N^{1/4} small-root bound.
 Lemma GAR: rep(x)/N ≡ (x_p − x_q)θ + x_q/N (mod 1), θ = (q^{-1} mod p)/p; θ to (3/4)n bits factors N. When can a
 class-E functional of rep(x)/N "see" θ? Is there a general "two-halves" decomposition: the CF of a/N = prefix
 determined by a/N, suffix = reversed prefix of ±ā/N, plus an O(1) junction? State it exactly.

Q2 (floor theorems). (i) For a uniform random unit a mod N, prove that the law of any class-E functional (or of a
 natural truncation, e.g. the parts of the path with convergent denominators ≤ N^{1/2−δ} at either end) depends on the
 factorisation of N only by O(N^{−c}), using fixed-modulus Euclid statistics (Porter, Ustinov, Bykovskii–Frolenkov,
 Aistleitner–Borda–Hauke, David–Kim–Mor–Shapira 2025), or Kloosterman/Weil bounds for the joint law of (a/N, ā/N)
 (note S(h,k;pq) has factor-dependent size when p | gcd(h,k)). Give exact statements and exponents.
 (ii) For structured a(N) (2^{(N−1)/2} mod N etc.), what conditional floor theorem can be proved? Formulate the weakest
 natural hypothesis (equidistribution of a(N)/N and ā(N)/N over semiprimes in a target class) and prove the implication.

Q3 (lattice-cone sums, d = 2, 3). Fourier–Dedekind sums s_n(a_1..a_d; N) = (1/N) Σ_{ζ^N=1, ζ≠1} ζ^n/∏(1−ζ^{a_j})
 (Beck–Diaz–Robins), Dedekind–Rademacher sums, Zagier's higher-dimensional cotangent sums. For N = pq:
 (a) Do they decompose like DSF (local parts at p, q + unit part; does the local part factor N; are there
 Rademacher-type congruences making N^k·s mod N an aggregate)? (b) Is the remaining (non-aggregate) content a
 CF-type functional for d = 2 and a "2-dimensional continued fraction / Barvinok decomposition" functional for d = 3,
 or is there something genuinely new? (c) Give an exact polynomial-time algorithm for d = 3 if you can (e.g. via
 weighted lattice-point counts / multidimensional floor sums Σ_{x,y<N} xy⌊(γ+αx+βy)/N⌋), with pseudocode.

Q4 (reduced-form coefficients). For binary quadratic forms of discriminant −4kN or 4kN, consider non-multiplicative
 statistics of the integer coefficients along reduction / composition / infrastructure paths (beyond the quadratic
 symbols ((−1)^i Q_i | p) = 1 of the principal cycle). Can you prove these reduce to class E (or class H, CF of
 quadratic irrationals) plus genus characters (Q), or find a mechanism that leaks something non-quadratic?

Q5. Propose (do not run at scale) concrete statistics for a battery on 2·10^5 balanced semiprimes (n = 48–128),
 with one-sided positive controls, targeting ε and [p/q ≥ 7/8] beyond the N mod 2^8 floor (plus top-bits floor for
 the archimedean target).

Be skeptical and precise; a scoped obstruction is as valuable as a leak. Output one markdown document: summary table
first, then per-question results with proofs, then a ranked list of the most promising leak mechanisms (if any).
