# masullo: "The Music of the Primes" (A. M. Masullo, draft, 2026-11-20): anything useful for us?

**Source.** Announced by Ravi Vakil on the Math 216 blog ("Hurried Announcement: Alessandro Masullo GRH"). Brian
Conrad has been checking it for months; the post says it is "all human, no AI". The draft is 940 pages:
https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf (sha256 e039d3f6…97bf5).

**Process.** Reviewed by GPT-6.2-astra xhigh via codex (session 01a1fa70-d434-7841-9a5d-a7f723d4dc09;
`codex/masullo_1_codex.out`) and Gemini 4 Carbon High via agy (conversation 29e4b377; `codex/masullo_1_gemini.out`)
from the extracted text. Brief: `codex/masullo_1.md`. This is not a verification.

## Verdict

- **If correct, it discharges every abelian (RH)/(GRH)/(ERH) hypothesis in the repo**, including the Dedekind/Hecke
  tier the alpo.ge Lean does not reach. **No factoring exponent changes**, and nothing touches α < 1/3.
- **Independent of GRH, it gives no usable tool for our open walls** (H, LP_exc^rel, R/reuse).
- **Keep all GRH qualifiers in the repo.** The deduction of the zeros is not yet written (see the stubs below).

## What is claimed (claimed-in-paper)

- **Scope.** GRH for L(s,χ) for every number field F and every unitary Hecke character χ of C_F, including infinity
  types, Dedekind ζ_F and Dirichlet L-functions. Imprimitive characters follow, since the extra Euler factors have
  zeros on Re s = 0.
- **Function fields** are explicitly *not* covered (§14.2). Residue norms are powers of q, so the polar-multiplicity
  theorem fails.
- **GL_n** (§14.3) is a 10-line stub ("ongoing work with Brian Conrad").
- **Multiplicity bounds.** §14.1 has m ≤ [F:Q] away from 1/2, but the October 2026 letter to Sarnak downgrades this
  to an expectation, and the letter says the a = 0 trace formula is impossible.

**What is written and what is a stub.**
- Chs. 2–11 and Apps. A–E, about 840 pages, are written in detail.
- The Main Theorem 12.5.1.1 (pp. 519–520) is written: a truncated explicit formula whose spectral side runs only over
  zeros with |Re ρ − 1/2| < a, for every a > 0.
- §12.3 (global-to-S-local assembly) is an empty heading marked "Include section".
- Ch. 13, the actual deduction of zero locations, is a 2-page stub: Theorem 13.3.1.2 says "Insert statement". The
  deduction exists only as a sketch in the unfinished introduction (§§1.4, 1.8).
- There are about 15 `??` references and about 20 check-notes.

## Method (claimed-in-paper)

This is a Tate / Connes / Meyer / Deninger-style spectral approach on spaces of adelic θ-functions.
1. Take Hardy–Sobolev completions with strip width a. Only zeros within a of the critical line survive continuously
   in the spectral trace.
2. The geometric kernel is lifted to GL₁(A_F)² and cut **one prime at a time** by "orbital Schauder bases". Using
   Q-linear independence of the log p, Poisson summation on R_{>0}/q_v^Z gives p-orbital traces that vanish for large
   p by compact support. This is the trick that avoids Meyer's divergent Euler product.
3. The S-local remainder is a Connes–Meyer computation. The 1/2 comes from |x|^{1/2}|y|^{1/2} on the diagonal
   converting multiplicative to additive Haar measure.
4. The paper explicitly uses **no Weil positivity**: it needs φ(C_F¹) = 0, which kills ψ∗ψ̃ tests. The geometric side
   is independent of a, and letting a → 0 puts the zeros on the line.

## Pressure points (both reviewers, independently)

1. **The sketched deduction (pp. 28–30 / 53–55) has a concrete problem.** Approximating 1_{(1,∞)}(t)t^{−z} by tests
   vanishing near 1 requires a cutoff-independent bound on ∫ φ(t) t^{−1/2}/(t−1) dt. Fatou shows this fails: the
   integrand behaves like 1/(t−1) near 1. Gemini proposes repairs: test against (log t)t^{−z}, or subtract Weil's
   explicit formula B.4.4.7 from Thm 12.5.1.1. Codex says only that the estimate needs replacing.
2. **The real crux is Chs. 10–12** (§11.3, the missing §12.3, and the unfinished (λ′, λ″) cutoff in Thm 12.4.1.1). An
   off-line zero enlarges the θ-closure and changes the projector P_θ. Everything hinges on whether the orbital
   chopping of this non-trace-class kernel leaves zero defect.

## For this project

**Discharged if true.** Everything in grhremark's inventory, plus the Hecke tier:
- RH for ζ: genus2rh Thm R (c = 3/√2; c = 2 still needs Younis); snfsrig O2.
- Dirichlet GRH: Props Z, EQ, PIN, DZ and GX, the DSR walks, Thm FN(ii), Prop IWD (its circuit premise stays), Prop ZA
  (its internal chain stays), selfred Thm L / leakmap PT-J1(b), and the factorbits successor reduction.
- Dedekind/Hecke GRH: nested Lemma X_ERH, splitting-field prime selection, Bach's generators, the dBPW sampler (it
  needs infinity types), the almostall13 C_χδ constant, and U-MF §8.4.

**Exponents are unchanged** (both reviewers).

**Possible tools, both NO.**
- §2.4 and App. A.4 count S-units ε with εg in the unit ball, O(T^{|S|−1}). These are classical Davenport/Schmidt
  volume counts for a *fixed* field, with ineffective constants. They count representatives of **one** unit orbit,
  whereas LP_exc^rel needs the number of **distinct** orbits contributing independent rows, uniformly as the degree
  grows.
- B.5–B.7 (recurrence for fixed finite Euler products, the classical zero-free region, Mertens on the 1-line) are
  standard.
- The geometric side's periods are log N𝔭 at known places. Nothing there can be evaluated mod a composite N.

**One minor use (codex).** Use App. A.2's log-unit fundamental domain to deduplicate candidates by a common base-field
unit in the multi-cone LP_exc^rel experiment. That is an experimental control, not an exponent gain.
