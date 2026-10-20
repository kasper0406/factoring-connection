# preprint: review adjudication for paper2/main.tex

Agent `preprint`, phase 11, 2026-10-20.
- Paper: `paper2/main.tex`, `paper2/refs.bib`, compiled `paper2/main.pdf`. Built with tectonic 0.17: no errors, no undefined
  references, no overfull boxes.
- Commits, in order:
  - `b82face`: first draft.
  - `22e9a2f`: Appendix A replaced by lem85's proof.
  - Opus fixes.
  - `f808719`: recursion section.
  - `c9e8f49`: codex fixes.
- Reviews:
  - GPT-6.1-sol `codex/preprint_1.{md,out}` (xhigh): whole paper, as of `b82face`.
  - GPT-6.1-sol `codex/preprint_2.{md,out}` (xhigh): the parts revised in `22e9a2f`.
  - One independent read-only Opus critic: whole paper, plus the revised Appendix A after a mid-run update.

## Verdicts

| reviewer | scope | verdict | weakest link named |
|---|---|---|---|
| GPT-6.1-sol preprint_1 | whole paper (old Appendix A) | **ACCEPT WITH MINOR FIXES**. No FATAL or MAJOR. Lemma 8.1 (one-point smoothing, new) "survives the detailed audit". Table 1 recomputed exactly; all entries correct. | old Appendix A's import of LV's majorant (now removed) |
| GPT-6.1-sol preprint_2 | revised Appendix A, citation updates, alternation | **CORRECT WITH FIXES** | Drappeau Lemme 5, whose published proof refers to Harper's unpublished §3.3 |
| Opus critic | whole paper, plus revised Appendix A | **ACCEPT WITH MINOR FIXES**. No FATAL or MAJOR. Every inequality chain traced. Table recomputed. | Fact D (Drappeau Lemme 5); next, Theorem L2 + LD* Step 5 (novel, deserves a human check) |

No reviewer found anything that breaks Theorem 1.1. The new Section 15 ("The recursive shape test": Theorem RA,
Corollary P) was added after these reviews, at the orchestrator's request. It is transcribed from `notes/recurse_results.md`,
whose own reviews were two GPT-6.1-sol rounds and two Opus critics. It is marked in the paper as **under independent audit**
(agent `auditRA`) and is not in the headline.

## Issue table and adjudication

C1 = codex preprint_1, C2 = codex preprint_2, O = Opus critic. Status: FIXED (applied in paper2), REJECTED (with reason),
MOOT (superseded).

| # | location | issue | sev. | found by | status |
|---|---|---|---|---|---|
| 1 | Thm FM, bad mass | β_ℓ counted all \|n\| ≤ 2MV/ℓ, but the dyadic cover only covers [MV/32, 2MV] | MINOR (C2 said MAJOR, but it is locally repairable) | C1 #1, C2 #1, O #1 | FIXED: β_ℓ → 𝔟_ℓ restricted to MV_ℓ/32 ≤ ℓ\|n\| ≤ 2MV_ℓ (only n ∈ I_P enter) |
| 2 | Thm FM, #𝒯_Q bound | cites Lemma local (i) only; needs (ii) and a polylog | MINOR | O #2 | FIXED |
| 3 | Thm FM, good mass | used local (iii) with Z ∈ (0, 32]; right end ≤ 32(λ+1) | MINOR | O #3 | FIXED: two-case argument with Z = 1 and Z = 33 |
| 4 | Lemma DS | D/B can fall slightly below N^{ξ₁} | MINOR | O #4 | FIXED (range with ξ₁/2) |
| 5 | 𝓡 definition, Lemma reject | N \| f′(m) is undefined when gcd(g₁, N) > 1 | MINOR | C1 #2, O #5 | FIXED: 𝓡 defined on gcd(g₁, N) = 1; nonunit g₁ is a factor exit |
| 6 | App. A, Lemma LS bilinear display | claimed false as written | — | O #6 | REJECTED. With M₁ = z^{1/2} and y ≥ 1, (z/M₁ + 4R²) ≤ (yM₁ + 4R²), so the product is ≤ (yz^{1/2} + 4R²)(yz)^{1/2} = y^{3/2}z + 4y^{1/2}z^{1/2}R² exactly. A one-line justification was added to the text. |
| 7 | effectivity of N₀ (§1.3(c), §13, §16) | Drappeau's η_D, c_D, δ_D are not stated effective | MINOR | O #7, C2 #2 | FIXED: "effective provided Drappeau's constants are; the algorithm uses no analytic constant, so Thm 1.1 is unaffected" |
| 8 | title box "all inputs published" | accurate only with the Harper provenance caveat | MINOR | O #8 | FIXED (the box and Remark 1.3(c) state the provenance) |
| 9 | §1.1 Seysen | listed as reaching exponent 1 under GRH | MINOR | O #9 | FIXED: Seysen rigorous under GRH; A. K. Lenstra brought the constant to 1. No constant is given for Seysen (not verified by us). |
| 10 | §1.1 LV bullet | says LV's techniques "are used" | PRES. | O #10 | FIXED: only the randomisation idea and the quoted Stark/Landau statements are used |
| 11 | HECM remark | needs a precise citation | MINOR | O #11 | FIXED: LPP I (1993) and II (2002) cited; Part III has not appeared |
| 12 | Lemma smooth1(c) | centre is (l_i − 1)/2 for either sign; \|M₀\| ≤ dV′^{d−1}/2 | PRES. | O #12 | FIXED |
| 13 | Prop R | Λ_ℓ((a, b)) undefined; partner count rests on Lemma reg | PRES. | O #13 | FIXED |
| 14 | notation clashes | ν(m′), β_ℓ, Γ, … | PRES. | O #14 | FIXED for ν(m′) → mult(m′), β_ℓ → 𝔟_ℓ, Γ (kern) → Γ_Φ. The other overloads are left, since each is local and unambiguous in context. |
| 15 | §13 recursion uniformity | "non-increasing envelope, polylog factor" is unjustified | PRES./MINOR | O #15, C1 #13 | FIXED: ε-argument |
| 16 | Lemma reject (iii) | w > H not guaranteed | — | O #16 | No change: the bound ⌈(2H+1)/w⌉/(2H+1) ≤ 1/w + 1/(2H+1) holds regardless |
| 17 | Lemma 2.4 | α ≤ 1 not shown | MINOR | C1 #3 | FIXED (saddle-point equation, monotonicity) |
| 18 | Lemma GB | tie when g₁ even | MINOR | C1 #4 | FIXED: −g₁/2 < a_i ≤ g₁/2 |
| 19 | §3 example | algebraic exponent printed as 59/156 | MINOR (typo) | C1 #5 | FIXED: 50/156 (sum 64/156 = 16/39) |
| 20 | Lemma goodmod | V/H ≤ N^{−ν−η₀} drops the L-factor | MINOR | C1 #6 | FIXED: N^{−ν−η₀/2} |
| 21 | Prop R covering | w₀ ≥ X_r/(2S²) has the wrong constant | MINOR | C1 #7 | FIXED: X_r/(4S²) |
| 22 | Lemma H | ‖f‖_∞ ≤ M(H+1) is false (the interior coefficient is g₁c_i − g₀c_{i−1}) | MINOR | C1 #8 | FIXED: ≤ 3M(H+1) < N for large N |
| 23 | t_χ | needs a computable \|H\| bound | MINOR | C1 #9 | FIXED: \|H\|_bd = 2^d\|D\|^d X₀(1 + log X₀)^{d−1} |
| 24 | KS operation count | logs omitted | MINOR | C1 #10 | FIXED (Õ) |
| 25 | Prop select (ii) | b = 0 not separated; M_rows depends on the field | MINOR | C1 #11 | FIXED |
| 26 | SPLIT for small N | structural guards (M < 4, J = ∅, …) | MINOR | C1 #12 | FIXED: step 0 "Guards" |
| 27 | old App. A Fact LV85 range | stated over LV's full y-range | MINOR | C1 #14 | MOOT (old appendix replaced) |
| 28 | character χ_𝔭 definition | γ𝔭^{−ord} is not an element | PRES. | C1 #15 | FIXED (uniformiser) |
| 29 | Lemma SQ numerator bound | justification omitted | PRES. | C1 #16 | FIXED (Mahler measure) |
| 30 | intro wording | "expected time at most"; "best single algorithm"; rational η₀ | PRES. | C1 #17 | FIXED |
| 31 | Table 1 η₀ column | rounded thresholds are not exact cutoffs | PRES. | C1 #18 | FIXED (≲, plus the exact expression √2 − 1 − θ in the caption) |
| 32 | novelty wording | "conditional or unconditional" is too broad | PRES. | C1 #19 | FIXED: "unconditional, or GRH-conditional, for factoring every integer". The 0.851 limitation item was replaced by §15 (labelled as an infimum). |
| 33 | bibliography | Joux–Lercier subtitle; Bernstein draft dates | PRES. | C1 #20 | FIXED |
| 34 | App. A paragraph after Lemma SM | "exactly / precisely when" overstated; the R₁^{1−β′} display is not an identity | MINOR | C2 #3 | FIXED (sufficiency only; inequality via (α2)) |
| 35 | App. A | P⁻(1) undefined | PRES. | C2 #4 | FIXED (P⁻(1) = +∞) |
| 36 | App. A constants lemma | the Euler-product wording conflates an upper bound with an identity | PRES. | C2 #5 | FIXED |
| 37 | Fact HT | α approximation is Thm 2(i), not Thm 3 | PRES. | C2 #6 | FIXED |

Not adopted: none, apart from #6 (reviewer error), #16 (no change needed) and the remaining notation overloads (#14).

## Other changes during this phase
- **Appendix A (orchestrator, after agent lem85).** Replaced the TODO appendix by lem85's self-contained proof:
  - lift-averaging identity;
  - Lemma SM (Drappeau Lemme 5 for small conductors, de la Bretèche–Tenenbaum via Drappeau Lemme 1(i));
  - Lemma LS / Corollary LS′ (Iwaniec–Kowalski Thm 7.13 plus a proved Perron bound);
  - the general Lemma A.1 and the deduction of Lemma 8.4.

  The paper now states honestly that Drappeau's Lemme 5 is published but that its proof refers to Harper's unpublished
  arXiv:1208.5992 §3.3. All citations of LV Lemma 8.5 and Harper's Prop. 2 as inputs were removed; the remaining mentions
  are historical. FACTOR's alternation of trial division with SPLIT is now credited with removing the N′ < N₀ branch
  (Remark `rem:noN0` in Appendix A), at a factor-2 cost.
- **Recursion (orchestrator, after agent recurse).** Added:
  - Section 15: Lemma (main proof at general a), Lemma OS, Wald for sub-calls, Lemma T_rec, Lemma F, Theorem RA,
    Remark (constants, what binds, fixed depth, prior art), Corollary P with the restart schedule;
  - a framed audit-status box;
  - an intro subsection with a footnote;
  - an abstract paragraph;
  - the recursion toy run;
  - §16 "Limitations and outlook" with the growing-degree targets 2/5, 3/7 and 1/3.

  The headline remains the audited statement (every c > 1/√2). Theorem RA and Corollary P are explicitly marked as not
  yet independently audited (agent `auditRA` is running). **The recursion section itself has not been reviewed in its
  paper form by these referees**; it was transcribed from the reviewed note.
- **dBPW citation.** The title was checked against arXiv (2512.01588: "Rigorous methods for computational number theory").

## Open items for the orchestrator
1. Report the `auditRA` verdict. If it is ACCEPT, the audit-status box, the intro footnote and the abstract qualifier can
   be relaxed. If it finds a gap, Section 15 needs revision.
2. A paper-level referee pass on Section 15 (as typeset) has not been run.
3. Unverified bibliographic details are listed in the header of `paper2/refs.bib`.

## Round 2 (orchestrator follow-ups): auditRA verdict, RA/P fixes, discrete-log section, codex preprint_3

**Changes.**
- The auditRA verdict (`notes/auditRA_audit.md`) was ACCEPT WITH FIXES for both RA and P, from three blind referees with no FATAL or MAJOR finding.
  - The §15 audit box is replaced by a short status remark: conditional on Theorem A, internally audited by three blind referees, not peer reviewed.
  - The intro footnote is updated. RA and P now appear in the abstract as further consequences. The headline sentence is unchanged.
- The five auditRA MINOR fixes are applied:
  - core-vs-trial-division race, with OS applied to the core inside the strong induction;
  - η₀_j = 1/d_j² and a computable parameter sequence;
  - "no effective rate supplied";
  - 0.851 stated as an infimum (k = 13 gives 0.9952);
  - the "|u₀v₀| ≍ N" claim was never in the paper.
- auditRA's sanity remark was added: RA's constant is 1.346× the heuristic 2/√(d+1), and the reasons the random-squares folklore does not bind.
- New §16, "Application: discrete logarithms in F_p^*", from `notes/dlog_results.md`. It contains:
  - Theorem DL and Corollary DL′;
  - the literature paragraph (Pomerance 1987 √2; Vollmer retraction; Biasse–Erukulangara; Hafner–McCurley; Enge–Gaudry; KW; Bender–Pomerance; dBPW);
  - Algorithm DLOG, the transfer lemma, Lemmas XK0/XK1–5 and HX, and Proposition XK with proofs;
  - the recursion-interaction remark (c_DL = c\* with RA; 3√θ/2 with the portfolio, where the 3a term takes over);
  - a prior-art and non-transfer remark, and the toy check.

**Codex preprint_3** (`codex/preprint_3.{md,out}`, §§15–16 at c8e1ee5). Verdicts: §15 CORRECT WITH FIXES; §16 CORRECT WITH FIXES. No FATAL or MAJOR. All
nine items were FIXED:

| # | issue | sev. | fix |
|---|---|---|---|
| 1 | \|a_j − a\*_j\| → 0 does not give c(a_j) → 0; the d = 9 portfolio needs k_j → ∞ and η₀_j → 0 | MINOR | a\*_j ≤ a_j ≤ (1+1/j)a\*_j; both limits stated |
| 2 | admissibility for all j asserted by "monotone, checked j ≤ 60" | PRES. | closed-form proof (D = d² + 2d − 11; μ − ν ≥ 1/d²; θ < 3/(d+1) + 5/d²) |
| 3 | DLOG may not terminate for tiny p (empty IR interval); no race | MINOR | race with exhaustive enumeration; p = 2 handled |
| 4 | abstract's DL constant lacks θ ≤ 1/2; "every optimised tuple" false (d=3, k=3, μ=7/24, η₀=1/100 gives θ = 581/600, c_DL ≈ 1.476); 3√3/4 is an infimum | MINOR | fixed |
| 5 | "not attained at fixed d" contradicts Corollary P | PRES. | distinguished one parameter tuple from one algorithm scheduling tuples at fixed degree |
| 6 | DL row cost omitted the 2a batch term | MINOR | max(2a, …), dominated by 3a; also in Remark DLrec |
| 7 | T_j "all coins" includes the label coins; step 4 is label-dependent | PRES. | T_j defined explicitly |
| 8 | "every gcd-checked quantity < p" is false for f′(m) | MINOR | stated precisely: no proper-factor exit, p \| f′(m) is a rejection |
| 9 | Biasse co-author; attribution of the Vollmer retraction | PRES. | Biasse–Erukulangara; retraction attributed to their published introduction |

Paper: 54 pages; compiles with no errors, no undefined references and no overfull boxes.
