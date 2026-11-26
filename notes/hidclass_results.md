# hidclass: the hidden factor's class in a fixed abelian field (outside Briefs 1 + 3) (2026-11-26)

**Process.** GPT-6.2-astra via codex (session 01a2176f-6e29-7222-940b-6f50a9634bd6). Brief: `codex/hidclass.md`. Output:
`codex/hidclass_codex.out`. Prereg: `notes/hidclass_prereg.md`. Code and results: `experiments/hidclass/`.
The outside briefs arrived with all formulas stripped; the targets were reconstructed in the orchestrator session.

**Verdict.** Nothing here gets below L[1/3]. Both targets are generalized Blum-integer recognition. The measurements
were null. Codex puts about 5% on a poly-time decider for either target.

## Targets
- **T1.** N = rq, r ≡ 1 mod 8, q ≡ 3 mod 8. The bit is [r ≡ 9 mod 16]. With A = Σ_{n≥0} z^{n²}, it equals
  [z^{N²}] A(z)(A(z^128)+A(z^256)) mod 2 (orchestrator, 600 N). Proved (codex): this is the XOR of r = x²+64y² and
  r = x²+32y², which is Kaplansky's criterion.
- **T2.** N = pq ≡ 1 mod 7. Decide whether p and q both split in the cubic field of conductor 7
  (F = x³+x²−2x−1) or are both inert.

## Proved / known
- **Ternary count (Brief 1's U).** U(N) = #{0<x<y<z : x²+y²+z²=N²}. It satisfies **8U = N + 2r − 21**, which comes
  from Hurwitz's formula r₃(N²) = 6(N+2r). Spot-checked on 5 N by the orchestrator. So U mod 2^j ⟺ r mod 2^{j+2}.
- **Totient forms.** T1 ⟺ v₂(φ(N)) = 4. T2 ⟺ φ(N) mod 7 ∈ {0,4}.
- **Blum recognition is one fixed-base QR query.** It asks whether −1 is a square mod N. Seurin calls the hardness
  assumption 2-Φ/4-hiding. QR on the bases −1, 2 and −2 gives the full mod-8 class.
- **T1 reduces to a 16th-power symbol.** For a uniform unit a, (a⁸|N)₁₆ is always +1 when b = 1, and −1 half the
  time otherwise (via Hittmeir 2016, which covers sums of two squares and v₂(p−1)).
- **T2 reduces to square membership in the norm-one torus.** Take t = σ(a)/a for a uniform unit a.
  - If p, q are inert, T has odd order, so t and t² are identically distributed.
  - If they split, the two distributions have statistical distance 15/16.
  - This matches kern's QR-in-a-torus reduction.
- **The N^{0.9}-box decision for F ≡ 0 mod N is exactly T2.** A root-sampling oracle factors N with probability 2/3.
- **Complexity.** Both targets lie in UP∩coUP and in SZK. No classical decider is known, no reduction from
  factoring is known, and neither target is known to be equivalent to QR or DCR.
- **Oracle lifting (known-multiple model).** H_m(N^e s) = eH_m(N) + H_m(s). So one answer simulates every such
  query, and a fixed-m oracle saves only a constant factor for factoring. Explicit CRT pairs show O₁₆ does not give
  O₃₂, and O₄ does not give O₈. The unrestricted reductions remain open. O_{2^j} with 2^j ≥ N^{1/4} factors N
  (Coppersmith).
- **Descent model.** Ray characters, fixed Gauss sums and fixed-level multiplier tables, with arbitrary
  postprocessing, factor through N mod M, so they exclude T1. Growing Euclidean state (Dedekind-type) escapes this
  model. Theta S/T transformations move the expansion to another cusp and do not shrink the index N².
  Also A = (1+θ)/2, so reducing mod 2 before dividing by 2 loses information. No Ω(N) lower bound is known.

## Measured (null)
24,576 balanced semiprimes, 40–80 bits. Statistics: CF(√7N) quotient-sum parity, CF(√7N) denominator parity, and
the cubic sign from Jacobi symbols of Norm(a−σ(a)). Pooled advantages were +0.59, −0.09 and +0.21 percentage
points, and every 95% CI contains 0. Proved: the cubic sign is exactly fair, since multiplying by a scalar unit of
Jacobi sign −1 flips it.

## Open
For T1, distinguishing 8th powers from 16th powers among units mod N. For T2, distinguishing square classes in the
norm-one torus. Also: an unrestricted lift O₁₆ ⇒ O₃₂, and approximating the Gaussian theta sum over solutions of
F ≡ 0 mod N to additive accuracy e^{−2π}/3.
