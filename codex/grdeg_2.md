# Adversarial review of grdeg's findings on the degeneracy-driven G_r failure

Review the note `notes/grdeg_results.md` (in the repo you have read-only) and its code `src/grdeg/`. This is a study
of whether the rank-≥3, no-rational-side ("G_r") case of integer factoring by NFS-type relation methods can beat
L[1/3], focusing on frame-less degenerate (subfield / Weil-restriction) modules. Your own earlier blind round is
`codex/grdeg_1.out`; its relative-resultant argument was adopted as Proposition RR.

Check these claims hard, and say for each CORRECT / CORRECT-WITH-FIXES / WRONG with specifics:

1. **Ledger (deliverable #1, §1, src/grdeg/ledger.py).** The achievable factoring exponent is α = 1/(4−c), where the
   N-saving (generator-height) exponent κ_gen has x-exponent 1 + c·ρ − δ: c=0 (κ_gen=1/d, the aggregate L[1/4]) vs
   c=1 (κ_gen=(r−1)/(D−1), giving L[1/3]). Is the exponent algebra right? Is the note honest that the endpoints are
   heuristic (outside rank3's box scope H5_r because the {0,±1} box has lattice error ≫ B)?

2. **Proposition RR (§3, src/grdeg/relres.py).** For a frame-less subfield degeneracy over k (deg e) with relative
   minimal polys f,g sharing the linking root mod N and relative resultant R=Res_k(f,g)≠0: N | Norm_{k/Q}(R), hence
   generator heights h₁+h₂ ≥ e·x/d − O(d), which for a rank-2e module restores κ_gen ~ r·x/(2D) and α≥1/3; and if R=0
   the sides descend to the same form giving trivial square congruences (useless). Is RR1 (ψ(R)=0 ⇒ N|Norm(R)) fully
   rigorous? Is the height consequence RR2/RR3 correct, and does it really restore 1/3? Does RR genuinely subsume the
   (withdrawn) factoring-hardness obstruction via "small generators ⇒ R small ⇒ R=0 or N≤Norm(R)"?

3. **Scope honesty.** The note claims RR closes the SUBFIELD frame-less case at 1/3 (via height, not construction
   hardness), and leaves open: (a) non-subfield frame-less modules (your `(1,s,t,t²)/(1,s+23,t,t²)` finite-volume
   example; decoupled fibre products), and (b) N-dependent subfield / tower-NFS (root free, height N^{1/e}), and (c)
   large boxes beyond B^{1/(r−1)} (LP_exc). Is this the right open list? Is there a subfield frame-less family that
   escapes RR (e.g. R=0 but still useful, or a relative resultant that does not charge N)? Try to construct one.

4. **Verdict.** Can the G_r case beat L[1/3]? If you still think there is a crack, give the most concrete version.
   Flag any remaining overclaim in the note's labels.

Be concrete, use exact small examples, cite real literature. Report as plain text.
