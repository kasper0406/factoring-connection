# cogentic: how do Google Research's Cogentic papers relate to our repo?

Directory: https://github.com/google-research/google-research/tree/master/cogentic_research_papers (posted 2026-11-25;
produced by "Cogentic", arXiv:2609.40324, an agentic framework with Gemini 4 as base model; "no human interaction other
than the initial prompt"). Local copies (READMEs + pypdf text, formulas may be mangled) in
/Volumes/git/factoring-connection/codex/scratch/cogentic/:
  integer_factorization{_README.md,.txt}     — unconditional Las Vegas complete factorization in L_N[1/3, 105],
                                              workspace L_N[1/3, 40]; Lean 4 formalization "in a companion submission"
  deterministic_factoring_l34{_README.md,.txt} — deterministic factoring O(exp(3(ln N)^{3/4})) = L_N[3/4]; algSplit
                                              given a multiple of p−1 in L[1/2]-type time; deterministic QNR in L_p[1/2]
  discrete_logarithm{_README.md,.txt}        — fixed-characteristic F_{p^n} DL in exp((2+o(1))(ln n)²/lnlnln n)
  README.md                                  — full list (also game theory / learning theory papers, not relevant)
Our repo: /Volumes/git/factoring-connection (do NOT edit). Read README.md, BARRIERS.md, STATUS.md (tail), and the
relevant notes: rigorous chain notes/thmA_proof.md, notes/mfU_theorem.md (U-MF L[2/5]), notes/thmAB_results.md,
notes/audit2L_audit.md (L[3/8]), notes/fruit_results.md, notes/lv6* / notes/ldbound*, notes/hypD*, notes/pcest*
(relation supply, fruitfulness, Lee–Venkatesan); deterministic line notes/detGRH_results.md, notes/zroute_results.md,
notes/dsr_results.md (DSR, Prop PIN, Miller), notes/avgtime_results.md (Thm AV), notes/equivtarget_results.md;
the model barrier notes/lowerbound_results.md; GRH notes notes/grhremark_results.md, notes/masullo_results.md.

Answer (≤ 1800 words, plain language; cite paper sections/pages from the .txt files and repo file:line):
C1. What exactly each of the two factoring papers claims, and the core technique (how does the L[1/3] paper get
    rigorous relation supply AND fruitfulness — the two things our rigorous chain could only get at 2/5, 3/8 and
    rigorous L[1/3] needed Conj F / Hyp AD + LV6′? how does the deterministic paper get L[3/4] without ERH — does it
    solve our DSR, or bypass it? how does it get QNRs deterministically in L_p[1/2]?).
C2. Precisely which repo results are superseded, matched, or contradicted (Thm A, U-MF, 2L-final, Thm DL/dlog25,
    detGRH/zroute/DSR, Thm AV, the lowerbound model theorem — note L[1/3] rigorous is consistent with α ≥ 1/3).
    Does anything in the papers contradict a proved repo barrier (e.g. dsr Thm SB / Thm FN, BOX-∞, the fruitfulness
    barriers in fruit, rkrank, nested R2)? If so, which side is wrong?
C3. Any ingredient that could be repurposed toward α < 1/3 (the repo's only goal), e.g. in light of our three escape
    hatches (H, LP_exc^rel, R/reuse) and the characteristic obstruction (notes/equivtarget_results.md) for the DL paper.
C4. Credibility read: internal consistency, suspicious steps, whether constants (105, 40) and the Lean claim look
    plausible, and the 2–3 places a referee should check first.
Use only read-only shell commands plus python3 in your cwd. Final answer ≤ 1800 words, plain text.
