# ptesplit round 2 (ADVERSARIAL): attack the synthesis on route 1's PTE escape

You are an adversarial referee. Read (read-only) these repository files:
- notes/ptesplit_results.md (the synthesis under review) and notes/ptesplit_formulation.md;
- notes/ptebarrier_results.md (Lemma R, Theorem W, Corollary P, Prop S, Prop NX — still a draft);
- notes/ptecons_results.md (LEDGER, CUN, Kleiman, LIN-R);
- notes/resbarrier_results.md §2.6 and §4.2 (rank-r linked modules, heuristic G_r);
- your own earlier blind answer codex/ptesplit_1.out.

Attack, with concrete counterexamples or computations where possible (scratch dir is cwd):
1. **Lemma NC** (norm compression): correctness of ∏|Ĝ(a_j,b_j)| = |N(r(θ))| ≤ ‖r‖₁^δ‖F̃‖₂^{d−1}, the reduction to
   an irreducible factor G with N | G(m), the unit counts (Thue/Bombieri–Schmidt, Pell), and the claim that members with
   the same non-unit part differ only by uninformative exact unit relations. Is the "element-factor-base NFS" reading right?
2. **Prop NF**: the semiprime counting bound and its proof; is the fraction e^{−x+o(x)} claim correct?
3. **Window lemma**: is k − 1 ≤ x^{3α−1} correctly derived, and is "generic Δe" a reasonable reading? Does Kleiman/the
   B^{1/2} bound plus "amplifier cannot beat its seed source" really close mechanism B, or can seeds be cheap for special Q?
4. **§2.4 coincidence claim**: is route 1's View-B residue really an instance of resbarrier §2.6 (rank-k linked module
   between two nonlinear orders), with View A excluded via KL? Is the converse correctly denied? Is anything in the PTE
   escape NOT covered by either (e.g. multivariate presentations, non-monic F, pieces not of the form bX − a)?
5. **Synthesis faithfulness**: does the note overstate ptebarrier (a draft) or ptecons? Any verdict line labelled
   stronger than its support?
6. **Experiments**: are the E1/E2 prediction evaluations (P1–P9) honest and correctly scored given the raw files
   experiments/raw/ptesplit_E{1,2}.jsonl and the frozen notes/ptesplit_prereg.md?

Output: numbered issues with severity FATAL / MAJOR / MINOR and a concrete fix each; then an overall verdict
(ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT). Under ~2500 words.
