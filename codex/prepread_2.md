# prepread_2: adversarial referee round

You are a hostile but fair referee (analytic number theory, Galois theory, arithmetic dynamics, algebraic
complexity, cryptography) with web search. Repository /Volumes/git/factoring-connection is read-only. Read
notes/prepread_results.md (§0–§5; §6–§7 are experiments still being filled in; ignore them) and, as needed,
notes/periodfind_results.md (Theorem T, §1.1; PF §6), notes/slpcirc_results.md (model U, W_U §1.3, Prop BB §1.5,
Theorem SP §4), notes/nonsmooth_results.md §5.1 (OA). Your own blind report for this phase is
codex/prepread_1.out; do not simply repeat it.

For EVERY labelled claim in §1–§5 give a verdict: CORRECT / CORRECT WITH FIXES / WRONG / OVERCLAIMED / KNOWN
(give the exact reference) and the specific fix. Be concrete: counterexamples, missing hypotheses, wrong
constants, wrong references or theorem numbers. Priorities:

1. **Theorem O1-DL and Corollary O1-N (§2.4).** Is the reduction correct (frequency finding via SFT over
   Z_{p−1}, the estimator, the doubling/backward unwrapping, error bounds, the gcd(j, p−1) cofactor, the
   no-generator remark)? In O1-N: the CRT factorisation Ψ̂(χ_pχ_q) = A_p·A_q with the stated multipliers, the
   order bound ord χ_p ≥ (p−1)/(κi), the averaging over pairs, the simulation of the multiplier c = h(pq′)q′^{−1}
   at a single prime, randomised h. Is the conclusion "O1 ⇒ prime-field DL for ≥ μ/2 of primes" right, and is it
   known (AGS 2003, Akavia, Boneh–Venkatesan, HNP literature)?
2. **Prop PH (§2.3)**: correctness of the bit-peeling program, size, the correlation bound, N-independence; is
   "order-only barrier is false even for N-independent programs" fairly stated?
3. **§2.1 Weil width**: the Laurent bound (a + b)√p for nontrivial χ, the cyclic-width definition, the threshold.
4. **Prop GAL-0 (§3.3)**, cases (a)–(c), in particular (c) (existence of ℓ, CRT choice of θ, primitivity, roots
   mod p distinct) and (b) (Eisenstein + types ⇒ S_D). And the consequence "permanent-hard whenever 𝒢 contains
   the trivial group".
5. **Prop CS (§3.2)**: the GRH prime-ideal-theorem error term and the "silent by a factor √X" conclusion.
6. **Prop AB (§3.4)**: Z_p = Σ n_χ χ(p); the interpretation; the GL₂ remark.
7. **Theorem MIN (§4.4)** — the main new claim. Check every step: invariant disc for monic integer polynomials of
   degree ≥ 2 (constants), Northcott count, prefix transfer inequality, the choice of d and m, Burnside on distinct
   roots of rad(P), Chebotarev in dyadic intervals, bad primes, compositions with < m layers, κ. Is the statement
   true? Try hard to break it: e.g. compositions whose prefixes have many rational/low-degree factors (x², x² − 1,
   x² − 2 layers; Chebyshev/power families), repeated roots, PCF maps, words where the prefix is fully split.
   Is it already in the literature (Odoni, Jones, Ferraguti, Benedetto–Juul, Hindes, Looper, "density of prime
   divisors of orbits", Northcott-type arguments in arithmetic dynamics)? Can the depth restriction
   k ≤ log₂X + Δ be relaxed with the same method?
8. **Lemma DEF, Props OUT/PFX, §4.5, Conjecture QC, §5 table.**
9. **§0 verdict**: is anything overclaimed, in particular "the two remainders separate exactly at the order of
   the character", "no per-prime obstruction", "O1 strictly stronger target than factoring"?

Output: a structured report (≤ 3500 words): overall verdict (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION /
REJECT), then the per-claim table, then details. Scratch files only under your working directory.
