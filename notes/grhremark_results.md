# grhremark: what a proof of RH/GRH would change for this project (2026-11-20)

Trigger: a 4-page note at https://alpo.ge/remark.pdf (dated 2026-11-20) announcing a proof of RH and GRH, with ERH
"expected with more time", formalized in Lean (alpo.ge/lean.tar.gz.part01..14, **not checked by us**). It has no
factoring content. Its only technical content is a sketch: a Stepanov-type argument with multiplicity, realised over Z
by a chain of ~√X iterated theta lifts between definite unitary groups; Guinand–Weil positivity with a variance test
function on a "path quadric"; GGP/Bessel periods by seesaw; an Iwasawa tower (Behrend); a Bost slope inequality; a
weakened Vandiver. It does not say which L-functions "GRH" covers.

Process: orchestrator (Claude) + two consultants with repo read access, two rounds (independent answers, then
cross-review). The consultants were GPT-6.2-astra via codex (`codex/grhremark_{1,2}.codex.out`, session
01a1f9b2-7261-7c03-8404-96b2aef52b71) and Gemini 4 Carbon (High) via agy (`codex/grhremark_{1,2}.agy*.out`). Briefs:
`codex/grhremark_{1,2}.md`. Status of this note: assessment; no new theorem.

## Verdict (plain language)

**Even granting RH, GRH and ERH in every common sense, no factoring exponent improves, and nothing points below
L[1/3].** All three parties agree. Every sub-L[1/3] bottleneck in the repo is a construction or counting problem
(resultant/height floors, the size–entropy law Thm SS, constructor costs, G_r independence). None is a statement
about L-function zeros. The sketched *method* gives no identity that inserts smoothness indicators, constructs
smooth numbers near kN, or controls principal-character drift.

## What becomes unconditional (if the claim holds)

- **RH for ζ only.**
  - genus2rh Thm R (L_p[1/2, c]) loses its only zero hypothesis. The checked inputs give c = 3/√2. c = 2 still
    rests on Younis arXiv:2409.05761v1 (unrefereed), a verification issue, not a conjecture.
  - snfsrig (O2), cofactor recognition, is discharged via genus2rh. (O1) and the bookkeeping remain.
- **Dirichlet GRH** (the repo's "ERH" in detGRH/zroute/dsr, Miller–Bach style; imprimitive characters mod N are
  used).
  - The reductions Prop Z, EQ, PIN and Thm DZ, together with Prop GX (existence), become unconditional reductions.
  - The DSR-side results (randomized walks, Thm FN(ii)) lose their GRH qualifier. So does Prop IWD: IW circuit
    hardness then gives deterministic L[1/2].
  - Prop ZA (L[1/2] / L[2/5]) loses GRH but keeps the internal chain.
  - Also unconditional: the selfred Thm L / leakmap PT-J1(b) model barriers, and the factorbits successor
    reduction.
- **Dedekind / Hecke GRH** (the repo's "ERH" in nested/escapeR2).
  - nested Lemma X_ERH (it needs ζ_K and every ζ_{K(√h)}).
  - Splitting-field prime selection (detGRH row 9).
  - Bach's class-group generators.
  - The dBPW sampler, which needs Hecke characters with infinity types.
  - The almostall13 C_χδ constant.
  - Simpler square certification in U-MF (§8.4).
- **Automorphic GL_n GRH.** Nothing further for factoring.

## What does not change

- **Deterministic factoring.** Harvey–Hittmeir N^{1/5+o(1)} worst case; Thm AV N^{1/6} average.
  - Deterministic subexponential factoring still needs **DSR_α**, a constructive lemma blocked by Thm SB, Thm FN and
    Prop BOX-∞.
  - The cov38b phase-discarding cap (4/25 < 7/25) is untouched.
- **Rigorous randomized.**
  - U-MF L[2/5] and 2L-MF (3/8, 2/5) are already GRH-free.
  - Below 3/8, GRH supplies none of the conditional covariances needed.
  - Rigorous L[1/3] NFS still needs fruitfulness (Conj F) or Hyp AD + LV6′.
- **Below L[1/3].** GRH does not touch any open door:
  - grdeg relative-rank ≥ 3 type-D (L[1/4] typical-size model);
  - Prop SI (constant only; RH gives h ≥ X^{1/2+ε}, existence not construction);
  - single a_f(N);
  - PRIM.

## Proposals adjudicated (round 2)

| Proposal | Verdict | P | Missing piece |
|---|---|---|---|
| Pool Lee–Venkatesan fields + XK-Z + PIN ⇒ rigorous Las Vegas L[1/3] (Gemini r1) | Valid conditional reduction; **GRH/remark irrelevant**. Random-base Miller splitting needs no GRH. The gap is exact Q-rank of NFS-derived rational relations (zroute §3.2 "seedless PIN": one row per field can cover every target with all W_j = 0) | 1–5% | Hyp PIN-LV: Pr[Q₀ ⊂ Σ_j W_j + Q·e_{−1}] ≥ 1/2 over J = L[1/3] i.i.d. LV fields |
| Bost/Arakelov slope inequality to close grdeg door (i) (Gemini r1) | Barrier direction only. **As stated it is vacuous**: at the L[1/4] balance the main term ((m−1)/(2n−1))·log N ≍ x^{3/4}, but the error O(dn) ≍ x (codex; still true of Gemini's revised "Relative R4"). The suggested Q(i), r = 6 witness is a power basis, so it does not test the non-power case | 3% (codex) – 35% (Gemini) | Relative R1/R4 over O_k with error o(x^{3/4}), plus control of independent low-norm box points under unit action |
| Remark's variance positivity for the 4/11 three-level target (codex r1) | The remark's method gives a ≥ 0 lower bound, but Prop DC needs an upper bound φ₂ ≤ η₂/4. Non-principal errors already decay, so the D^{−2α} wall is principal-character drift (gcd / β-drift), which zeros do not see. LF-A/OPV-A are upper bounds in the two-level design, not a proven floor at depth 3 | <1% via the remark; 3–15% via classical DC_e conditioning | Coarse gcd/archimedean cells giving E[φ_{2,e}], E[σ²_{2,e}] ≤ η₂/4 at α = 4/11 |

## Recommended next threads (merged; none depends on the remark)

1. Consolidate rigorous L[3/8]: Hyp U′ under DW, plus merging 2L-MF and SPR-P′ into one referee-ready note.
   60–70%, 2–4 weeks.
2. L[4/11] endpoint variance estimate by classical conditioning (Prop DC_e). 3–15%, 4–8 weeks.
3. Type-D relative rank ≥ 3 with non-power generators, in both directions: barrier 3–35%, construction < 1%,
   1–3 months.

## Hygiene

- Keep every (RH)/(GRH)/(ERH) label until the Lean development is independently checked and a specialist has read it.
  Do not swap GRH shortcuts into the GRH-free U-MF/2L-MF/AV arguments.
- A meaningful check for us:
  - Pin the archive hashes and toolchain, and build in a networkless sandbox (Lake scripts run code).
  - Grep for sorry/axiom/opaque/extern/implemented_by/native_decide.
  - Run `#print axioms` to confirm only {propext, Classical.choice, Quot.sound}.
  - Diff the vendored Mathlib definitions (riemannZeta, DirichletCharacter.LFunction).
  - Replay the build through an independent kernel (the leanA tooling: lean4lean / Nanoda / comparator).
- Check the **exact** statements. We need:
  - Mathlib `RiemannHypothesis`;
  - Dirichlet GRH for all moduli and all χ, imprimitive included or with the Euler-factor bridge;
  - for nested/dBPW, Dedekind ζ_K for all number fields and Hecke characters with infinity types.
- A clean axiom report can still certify a conditional implication. The vacuous-CEP lesson in leanA applies.
