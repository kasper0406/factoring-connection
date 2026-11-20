# grhremark round 1: if RH/GRH is now a theorem, what changes for this project, and is there a path below L[1/3]?

You are a research consultant on the repository at /Volumes/git/factoring-connection (read-only for you; do not
edit files there). Goal of the project: the asymptotically fastest integer factoring (algorithms, not constants).

## Context you must read first (in this order)
1. README.md, BARRIERS.md, then the last ~25 lines of STATUS.md.
2. GRH/ERH-relevant notes: notes/detGRH_results.md, notes/zroute_results.md, notes/dsr_results.md,
   notes/nested_results.md, notes/genus2rh_results.md (skim verdicts and the exact hypothesis used: which
   L-functions, which zero-free statement, explicit or not).
3. Below-1/3 frontier: notes/grdeg_results.md §0, notes/midk_results.md §0, notes/completeness_audit_results.md,
   notes/suffmap_results.md (Thm SS), notes/sicheck_results.md / notes/siconstruct_results.md (Prop SI).
4. Rigorous chain: notes/mfU_theorem.md (L[2/5]), notes/thmAB_results.md (L[3/8] window). Grep `GRH|ERH` in notes/
   for anything else.

## New external input (untrusted, unverified)
A 4-page note at https://alpo.ge/remark.pdf (dated 2026-11-20) claims a proof of RH and GRH, "expect with more time
ERH too", formalized in Lean (tarball at alpo.ge/lean.tar.gz.part01..14, not yet checked by anyone we know). Its
only technical content is a sketch: a Stepanov-type argument "with multiplicity" realised over Z by chaining ~√X
iterated theta lifts between totally definite unitary groups (U(X_0)→…→U(X_m), π ↦ π ⊞ μ_j[t]); positivity from
Guinand–Weil explicit formula with a variance test function and a sum over a "path quadric"
x_1x_2+…+x_{k−1}x_k; GGP/Bessel periods via seesaw (Rallis); Iwasawa tower (Behrend); Bost slope inequality;
a weakened Vandiver for leftover H^1 classes. It does NOT say which L-functions "GRH" covers (Dirichlet? Dedekind?
Hecke? automorphic GL_n?) versus "ERH".

Treat the claim as a hypothetical: "suppose RH + GRH (and, separately, ERH in each common meaning) are theorems".

## Questions (answer each; cite repo files/sections for every claim about the repo)
Q1. Inventory: which project results or reductions become unconditional (or less conditional) under (a) GRH for
    Dirichlet L-functions, (b) GRH for Dedekind zeta / Hecke L-functions, (c) GRH for automorphic L-functions? For
    each: statement, file, what remains conditional afterwards (e.g. DSR, G_r, C₁/C₂, LD*, internal unrefereed chain).
Q2. Deterministic factoring: Prop PIN / Miller splitting / zroute Theorem DZ. With GRH a theorem, exactly what is the
    best deterministic complete-factoring algorithm the repo can now claim unconditionally? Is it still N^{1/5} worst
    case / N^{1/6} average (Thm AV), or does anything subexponential become deterministic?
Q3. Rigorous-randomised side: does GRH (any version) let us lower the rigorous exponent below 2/5 or 3/8, or make
    rigorous L[1/3] NFS possible? Be concrete about which hypothesis in the chain GRH would or would not discharge
    (nested Prop R2 says "no ERH route below 2/5" — re-examine it).
Q4. Below L[1/3]: none of the open doors (grdeg relative-rank ≥3 type-D, G_r favourable failure → L[1/4] model,
    Prop SI constant improvement, a_f(N), PRIM) is an L-function-zero statement. Does GRH-as-theorem touch any of
    them? If not, could the *method* sketched (Stepanov multiplicity via iterated theta lifts, Guinand–Weil
    positivity, √X chain length) plausibly yield tools that do — e.g. smooth numbers in short intervals / near kN,
    pair-correlation or BDH-type estimates for smooth numbers in progressions, point-count or character-sum bounds
    beating Burgess, counting hypotheses like G_r? Be specific and skeptical; name the single most promising
    conversion and what statement it would need to deliver.
Q5. Give your verdict: is there a path post L[1/3]? List at most three concrete next threads for the project, each
    with (target statement, why the remark is relevant or not, probability of success, cost), ranked.
Q6. Credibility hygiene: what should the project do or not do with this claim (e.g. never mark results
    unconditional until the Lean is independently checked; what minimal checks of the Lean would be meaningful for
    us — which exact Lean statement would we need to see).

Be concise: ≤ 1500 words, bullets, plain language. Mark each claim proved/KNOWN/heuristic/speculative.
