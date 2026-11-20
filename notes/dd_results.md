# dd: the divided-difference test for LP_exc^rel (2026-11-20)

**Question.** In the slow-product-growth regime (divided-difference modules over k = Q(√2) with pole polynomials of
total degree p = O(m)), chase Lemma 2 showed a failure of LP_exc^rel *requires* conditioning loss
log(M_i/t_i) ≥ b_k·m log m at a positive fraction of the D embeddings. Does the shared denominator allow that?

**Runs.**
- Codex: prover and measurement. Session 01a1fabf-adf7-7392-9cb2-c797b5ccd41c; `codex/dd_codex.out`;
  `experiments/dd/codex/` (proof.txt, prereg, 9,096-embedding CSV).
- Gemini: adversary. Conversation aca5d478; `codex/dd_gemini.out`; `experiments/dd/gemini/` (proofs.txt, prereg,
  results).
- Briefs: `codex/dd_{common,codex,gemini}.md`.

## Verdict (both consultants independently)

**Necessary is not sufficient.** The local conditioning loss *can* occur at a positive fraction of embeddings, so the
height-uniform version of the threshold inequality is **refuted**. But no construction produces fresh independent
small-norm candidates. **No sub-L[1/3] yield; the ideal-relative threshold inequality remains open.**

## Proved

- **Cancellation identity (C).** After clearing U = ∏R_j, a ternary row's numerator is G_a = V_a + √2 W_a, with V_a and
  W_a of degree ≤ p+1, and the common denominator cancels exactly in t_i/M_i.
- **Why the resultant / Liouville route stops (C).** For a fixed row, the product over the 2n embeddings is the integer
  |N Res(P−√2, G_a)|. But if one row loses λ at u embeddings, then u(λ + log Cm) ≤ H + D log(Cm). At H ≍ Dm log m this
  allows u ≍ D: other conjugates compensate.
- **Loss at half the embeddings (C), with zero yield.**
  - Construction: q = m−2 quadratic poles T² − r_j². Interpolation plants a column of an Alon–Vũ ill-conditioned
    inverse near zero, and adding a large multiple of T^{n/2}(T−1)U preserves the residues and the link.
  - Result: by Rouché, log(M_i/t_i) ≥ ½ m log m − O(m) at D/2 embeddings, with both sides irreducible and unit
    generator ideals.
  - Yet for any prescribed cap Y < N², the parameters can push every nonscalar row above the cap. Exactly six scalar
    units survive (plus ±√2 if 2^{2n} ≤ Y), so the yield is zero.
- **Loss at a (1 − (m−2)/n) fraction (G).** Linear poles T − j, with residues interpolating an Alon–Vũ / O_k-staircase
  inverse column. Certified losses ≥ 0.75 m log m appeared on 92%, 83%, 75% and 67% of the D = 48 embeddings for
  m = 3..6 at n = 12.
  - Obstruction: the alignment forces root–pole separation ρ ≥ exp(Ω(m log m)). Through the exact identity
    w_{2+j}(α_i) = ∏_{l≠i}(σ(r_j) − α_l) this inflates Q ≥ Ω((1−m/n)·D m² log m) ≫ x/T.

## Measured

- **Codex (pre-registered).** One split semiprime each at 60, 72, 84, 96 and 112 bits; m = 3..8; p = m..2m; 195
  pairs, 9,096 embeddings; certified interval widths ≤ 3·10⁻⁹. The normalized loss r = ℓ/(m log m) had medians 0.49,
  0.38, 0.34, 0.25, 0.17 and 0.14 for m = 3..8, so it falls with m. The upper tail persists: at m = 8, 33 of 45 pairs
  had some embedding with r ≥ ½. Real embeddings account for 212 of the 214 values with r ≥ 1.
- **Gemini (pre-registered).** Six families plus cyclotomic/Pell pole supplements, 60–112 bits, ball arithmetic plus
  exact PARI ideal parities.
  - The low-H family C1b drives Q negative (−10.3 at 112 bits), yet every full-support candidate pays
    log|F| ≥ 1.24 log N (c = 1).
  - **N_ind = 0 across all 61 adversarial boxes at 84–112 bits.**

## Where this leaves Lemma 2

The threshold inequality cannot be proved height-uniformly, because local losses are realizable at Θ(D) embeddings.
The remaining target is the necessary range H = O(x/T): bound the contribution to Σ_i log t_i from the conjugates that
compensate for the local losses. Equivalently, show that loss-aligned modules pay in Q (as in G's ρ-inflation
identity) or in nonscalar norms (as in C's construction). Both constructions failed for exactly these two reasons.
