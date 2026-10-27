# Phase 28 synthesis: unexpected connections to factoring across 40 fields

Orchestrator summary, 2026-10-27, of four survey notes:
- [`conn_alg`](conn_alg.md): algebra and number theory;
- [`conn_ana`](conn_ana.md): analysis, physics, probability, optimisation;
- [`conn_top`](conn_top.md): topology, combinatorics, tilings, logic;
- [`conn_cmp`](conn_cmp.md): complexity, operator algebras, quantum, crypto, ML, numerics, arithmetic statistics.

The brief is [`AGENT_PHASE28_SURVEY.md`](AGENT_PHASE28_SURVEY.md). Each survey merged a GPT-6.1-sol blind round and passed one Opus critic (all ACCEPT/OK WITH FIXES).

## Verdict

About **300 concrete connections** in 40 fields: 78 + 79 + 81 + 60, including about 38 contributed by GPT-6.1-sol
and the critics. **None gives, or credibly points to, factoring below MNFS L[1/3, 1.902], below the rigorous
L[2/5], or below the deterministic N^{1/5}.** Each fails at a named step. The surveys agree, independently, on why.

## Why they fail: seven recurring mechanisms

1. **Restatement / FACT.** The invariant is p + q, φ(N) or σ(N) in disguise, so computing it is factoring. Examples:
   - areas and genera of modular curves;
   - class numbers of GL₂(Z/N);
   - #K_{2j−1}(Z/N);
   - irreducible dimensions of Aff(Z/N) and the Weil representation;
   - Kronecker–Hurwitz sums;
   - tiles of Z_N (|A| ∈ {p, q});
   - Galois-summed Bost–Connes states.

   Either the object has size ≥ N^c, or you need a factor to write it down.
2. **Aggregate.** The invariant depends only on N mod m: reciprocity, fixed-level quantum invariants (Ng–Schauenburg), Fermat-surface Picard numbers, fixed-ℓ counts.
3. **Ring arithmetic mod N (Theorem PT).** PIT, crossed products, planted isogeny walks and canonical heights (Müller–Stoll) see the factors only by splitting.
4. **Group order / smoothness (family M).** Dynamics (Ritt/Schur: only power, Chebyshev and Lattès maps give p ± 1 / ECM; every other map is rho-like), K-theory of finite fields (Bach–Shallit tori), Deligne–Lusztig tori (Williams' p + 1), Tal's factorial, and fixed-width tiling recurrences.
5. **Sampling precision.** Probabilistic, spectral and physical statistics see the factors only through 1/p corrections, so they need about N samples. Smoothing (heat kernels, Li's divisor-sum algorithm) hides the jump at N.
6. **CRT tensor structure.** Operators built from N split along p ⊗ q. Finding the split is the CRT map, and splitting spectra over Q is factoring x^N − 1.
7. **Distributional statements with no algorithm for a given N.** Cohen–Lenstra, Malle, optimal transport and ML.

## Small facts new to the project (none changes a bound)

- **Lens spaces.** Oriented homotopy equivalence of L(N, a) and L(N, b) is exactly quadratic residuosity mod N, under a randomized reduction (Reidemeister–Whitehead). d-invariant collisions are Fermat/Lehman in disguise. [conn_top]
- **Modular polynomials.** For ℓ with X_0(ℓ) of genus 0, finding a root of Φ_ℓ(j, X) mod N is randomized-factoring-equivalent. This extends Altuğ–Chen Thm 3.3; the toy split rate was 0.443 against a predicted 4/9. [conn_cmp]
- **Class number sums.** Σ_t H(4N − t²) = 2(N + q) (Kronecker–Hurwitz), and Σ_t (t² − N)H(4N − t²) = −2(1 + p³) (weight-4 trace formula). Both are checked exactly; both are FACT. [conn_alg]
- **Lozenge tilings.** The (1, x, x) hexagon has C(2x, x) tilings, which splits N only for x already within a factor of 2 of p. [conn_top]
- **Extension complexity.** The bit-multiplication polytope has extension complexity ≥ 2^{Ω(√n)}, because one face projects onto a cut polytope (Prop XC, a critic's sketch). [conn_ana]
- **Attribution fix.** heurnfs row 20: the MNFS-vs-ω curve for ω ≥ 1.78 is already in LSTT 2002 §2.7. [conn_cmp]

## Cross-field shortlist

Two surveys independently pointed back at the same existing openings. Nothing new survived.

1. **NFS fruitfulness: lemmaM's Lemma M / nfsfruit's Conjecture F (Lee–Venkatesan Conj 7.1, batch form).** Flagged by
   conn_alg and conn_cmp. A proof would give **rigorous L[1/3, (64/9)^{1/3}]** through Lee–Venkatesan's randomized
   NFS. That is far below the published rigorous L[1/2, 1] and this project's L[2/5], though not below heuristic NFS.
   - First lemma: an involution that preserves the square-class features and flips β, on a set of measure
     L[1/3, −o(1)].
   - lemmaM (phase 6) found no proof and refuted averaging over N as a route.
   - Feasibility with current tools was assessed as low.
2. **Deterministic worst case below N^{1/5}: coverage38's positivity count.** Can Bourgain–Garaev-type bounds for
   inverses of band primes give positivity in unions of Farey arcs at length q^{0.32}? [conn_alg]
3. **Constant-level only.** None changes the exponent; the best case is about 1.86.
   - Sub-quadratic sparse GF(2) linear algebra: MNFS ≥ 1.8595 at ω ≲ 1.667.
   - Near-regular-matroid elimination. Kill test: look for Fano or dual-Fano minors in a post-merge CADO matrix.
   - Percolation-guided relation selection.

Subjective probability, per the surveys, that anything here improves a project benchmark: below 1–2%.
