# analytic: candidate cards (BRIEF §4 template)

Agent `analytic`, phase 17. The verdicts point to `notes/analytic_results.md`.

---
**A1. Floor-sum ladder (reduction L_j).**
- **Mechanism.** Find divisors of N block by block. A degree-j Taylor polynomial P of N/(x+h) plus an
  exact floor-sum oracle counts and locates the few h with {P(h)} within ρ of an integer.
- **Target.** All N.
- **Representation.** Rational degree-j polynomials with denominator x^{j+1}, built from N and the block
  centre x.
- **Extraction.** Differences of two floor sums count candidates; bisection locates them; one division
  tests each.
- **Why it exposes a factor.** Every divisor in the block is a candidate (sign of the remainder).
- **Known lemmas.** Taylor/geometric series; Shiu's Brun–Titchmarsh bound for τ in short intervals.
- **New lemma required.** Exact degree-j floor sums in time T_j = H^{o(1)}: lattice points under
  near-parabolic arcs, below the hull complexity H^{1/3}.
- **Correctness.** PROVED conditional on the oracle (§3).
- **Cost.** N^{1/(j+2)} log² N oracle calls. j = 1 is operational (Euclid): N^{1/3}.
- **Randomness.** Deterministic.
- **Closest prior method.** Hiary (j = 1 analogue); D(N) convex-hull algorithms. The actual difference
  is the explicit degree ladder and candidate localisation.
- **Smallest experiment.** E3, done: slopes 0.336, 0.252, 0.201; O(1) candidates per block.
- **Refutation criterion.** ≍ H candidates per block, or a missed divisor. Neither was seen.
- **Repair if refuted.** Not needed.
- **Verdict.** A survivor as a precise open primitive. No evaluator is known for j ≥ 2.

---
**A2. Complete cubic Weyl sums (reduction C3).**
- **Mechanism.** Nine values S(a^i b^j, N) form a rank-one array 4√N u vᵀ. The cubes of Σu_kω^k and
  Σv_lω^l give π_pπ_q and π_pπ̄_q; their gcd in Z[ω] gives p.
- **Target.** N = pq with p ≡ q ≡ 1 (mod 3).
- **Extraction.** Exact Z[ω] arithmetic on oracle values.
- **Known lemmas.** CRT for exponential sums; cubic Gauss and Jacobi sums (Ireland–Rosen Ch. 9).
- **New lemma.** None. The barrier it gives: approximating S to O(1) is factoring-hard for these N.
- **Cost.** About 15 oracle calls expected, plus polylog. Each oracle call is O(N) classically.
- **Randomness.** 48/81 success per attempt.
- **Closest prior method.** van Dam–Seroussi's Gauss-sum estimation versus discrete log (multiplicative
  characters, quantum).
- **Experiment.** E2, done: 100/100, typical δ* ≈ 3 (rounding decoder).
- **Precision.** Codex: a conic small-root decoder (Coron 2004) tolerates N^{1/6−ε}. So constant precision
  is sufficient but not necessary.
- **Refutation criterion.** Failure on exact values. None was seen.
- **Verdict.** A barrier, not an algorithm. The degree-k extension is REFUTED (k = 4, N = 77).

---
**A3. Hull walk.**
- **Mechanism.** Enumerate the vertices of the integer hull of {xy ≥ N} (Alcántara et al.). Divisors
  are vertices.
- **Cost.** Θ(N^{1/3} log N) vertices (Balog–Bárány) and O(N^{1/3} log² N) operations (Alcántara et al.).
- **Verdict.** KNOWN (Alcántara et al.). Hiary has the same order of cost but is not a hull walk. Barrier H
  bounds enumerating methods only.

---
**A4. Band-limited sawtooth / exponential-sum divisor detector.**
- **Mechanism.** Σ_{d∈I} F(N/d) with F a trigonometric polynomial of degree K, evaluated frequency by
  frequency with Hiary-type theta blocking.
- **Kill.** Barrier F gives K ≥ 3x/8 − 1. Cost ≳ x > |I|.
- **Verdict.** REFUTED for this model.

---
**A5. Theta-series coefficient extraction (r₄ = 8σ, r₂).**
- **Kill.**
  - Sampling aliasing needs ≍ N samples (codex analytic_1 §2.4).
  - Exact singular series hits barrier R: the terms with modulus < p are N-independent.
  - r₂(pq) ∈ {0, 16} carries one bit.
- **Verdict.** REFUTED for these models.

---
**A6. L(1, χ_{−4N}) and the class number.**
- **Cost.** AFE costs √Q = N^{1/2}. Shanks' Euler-product + BSGS hybrid costs Q^{1/5} under GRH.
- **Verdict.** No route below N^{1/5}. h(−p) is class (b).

---
**A7. Gaussian product isolation M₁/Z** (codex analytic_1 §6B).
- **Kill.** The phase is quartic. Fourier linearisation returns a coefficient-extraction integral.
- **Verdict.** Killed.

---
**A8. Classical evaluation of Shor sums Σ_x e(h a^x/N).**
- **Kill.**
  - The phase class is closed under differencing.
  - Bit-cut rank ≍ min(2^s, 2^{m−s}, r).
  - The QFT part is already classical.
- **Verdict.** No route (§5).

---
**A9. Square-sieve / Fermat exponential sums.**
- **Kill.** A second-moment argument needs ≳ T characters for T candidates, with T ≍ √N.
- **Verdict.** Exponential.
