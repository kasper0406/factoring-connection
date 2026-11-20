# masullo_ext round 2: check Gemini's Theorem Rel-BT and reconcile with your cover result

Gemini's full answer: /Volumes/git/factoring-connection/codex/masullo_ext_1_gemini2.out (section 1 "Theorem Rel-BT"
statement + proof; section 3 toy checks in /Volumes/git/factoring-connection/codex/scratch/mext_g/). Your own round-1
answer is in this session (cover by embedding thresholds t_i; the 72-bit divided-difference witness).
≤ 1100 words:
R1. Is Theorem Rel-BT correct? Check (i) the divisibility v_{p_k}(Δ_rel) ≥ P(m,2n) = m(m−1)(2n−2)!/(2n−m)! (the
    row-subtraction argument, ramification e_K, e_K′, Galois invariance/integrality of Δ_rel, and whether it really
    answers grdeg's open "relative R1"); (ii) Prop Rel-R4; (iii) the unit-cone count via Masullo A.2 over the FIXED k
    (poly(x) cones; R_unit bound) and the volume bound V_rel·Y^{m/(2n)} from Bean–Thunder / Thunder over k; (iv) the
    boundary error r·D·(2E+1)^{r−1}. Give counterexamples where it fails; say exactly what it proves and in what scope
    (it is stated to close grdeg door (i) only inside the box scope H5_r — confirm that this excludes the L[1/4]
    ternary-box regime).
R2. Is Rel-BT new relative to notes/rank3_results.md (R1/R4), notes/rankdegen_results.md, notes/grdeg_results.md
    (Prop RR, "relative R1 blocked by infinite unit volume"), and Theorem KH? Is Masullo's A.2 actually needed, or is
    this just the standard Dirichlet unit fundamental domain over k?
R3. Your 72-bit witness vs Gemini's P3 observation (all small hits in proper k-subspaces; union spans full rank via
    relative-unit cancellation): what is the single cleanest lemma to attempt next (Gemini proposes a rank-3
    "two-hyperplane intersection" relative Liouville lemma; you propose Σ log t_i ≥ c_k x/T − O_k(Dm))? Pick one,
    state it precisely, and say whether it would close the L[1/4] ternary-box regime if proved.
