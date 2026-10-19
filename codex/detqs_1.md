# detqs_1 — blind round: deterministic quadratic sieve on almost all N

You are a research mathematician (analytic number theory + computational number theory). Be rigorous and
honest; label every claim PROVED / CONDITIONAL / HEURISTIC / FALSE / UNVERIFIED. The repo is read-only for you;
scratch files go in your current directory. Context files (read them): ../../notes/polyimpl_results.md §7 (D1),
../../notes/below16_results.md §6 (relation methods, H1/H2), ../../notes/lemmaM_results.md §§3–5 (Theorem A/B,
barrier B1), ../../notes/hypD_results.md §4.1′, ../../notes/avgtime_results.md §0.

## Setting
n-bit N, P = 2^{n/2}, S_n = {pq : p < q ≤ 2p primes, 2^{n−1} ≤ pq < 2^n} (density ≍ n^{−2} among n-bit integers).
Deterministic single-polynomial QS: m = ⌈√N⌉, Q_t(N) = (m+t)² − N for 1 ≤ t ≤ T; factor base = primes ≤ y
(with (N/ℓ) ≠ −1); exact sieve; exact deterministic F₂ elimination. Z_N = #{t ≤ T : Q_t(N) is y-smooth}.
- (H1) Z_N ≥ π(y) + 2 (or + K).
- (H2) some kernel vector S gives gcd(X_S − Y_S, N) proper, X_S = ∏(m+t_i), Y_S = √∏Q_{t_i}.
Fixed-u regime: y = N^c, T = N^{c+o(1)}, cost N^{ωc+o(1)} (deterministic fast elimination), < N^{1/6} for
c < 1/(6ω). Subexponential regime: y = L[1/2, a]. Project's best deterministic average/almost-all exponent is
1/6 (Theorem AV); worst case N^{1/5} (Harvey–Hittmeir).
Transfer: a bad set of ε_n·2^n n-bit integers is an O(ε_n n²) fraction of S_n, so all-integer statements need
ε_n = o(n^{−2}).

## Tasks (answer each; proofs where you can)
1. **H1 all-integer second moment.** Compute E_N Z_N and E_N Z_N² over uniform n-bit N. Identify the exact
   change of variables for the pair sum Σ_N f(Q_t)f(Q_{t'}) (f = 1_{y-smooth}). Is the coupled region thin or
   fat after summing over t at fixed k = t′−t? Which smooth-number theorems (Bombieri–Vinogradov or
   Barban–Davenport–Halberstam for smooth numbers — Harper 2012; Hildebrand short intervals;
   Matomäki–Radziwiłł; Drappeau; Fouvry–Tenenbaum; Granville; Soundararajan/Harper APs) give it, at which
   precision? Is Var(Z_N)/(E Z_N)² = o(n^{−2}) true unconditionally, or is there a local (small-prime / N mod ℓ)
   contribution of exact order n^{−2}? If the latter, what stratification fixes it, and does the fix survive
   in the L[1/2] regime (where log y ≍ √(n log n))?
2. **Write the strongest H1 theorem you can prove** (fixed u), with all inputs cited precisely (theorem
   numbers, ranges of y, level Q, error term).
3. **H2 for almost all N.** Any route? Examples to examine: (i) size arguments (are dependencies of size 2
   always fruitful?); (ii) Jacobi form ψ_N(S) for Blum N and reciprocity; (iii) averaging over N with the
   N-dependent kernel (contrast with lemmaM Theorem B); (iv) multipliers kN / MPQS families; (v) barrier B1.
   Either prove H2 for almost all N, or give a precise barrier/obstruction for natural proof methods, and the
   exact missing statement.
4. **Literature.** Is any deterministic almost-all (or average-case) factoring algorithm below N^{1/5}, or
   subexponential, known? Search: Pomerance, Schnorr, Lenstra, Vallée, "deterministic quadratic sieve",
   Hittmeir, "almost all integers" factoring, Rabin/Lehman/Hart one-line. Give exact references.
5. **Verdict**: what is provable now (c, or L[1/2, c]), what is conditional, and the single most promising
   missing lemma.

Keep the final answer ≤ 2500 words, structured by task number.
