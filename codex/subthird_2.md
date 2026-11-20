# subthird round 2: cross-kill (goal unchanged: factoring general N in L[α], α < 1/3)

Common rules: codex/subthird_common.md. Round-1 outputs:
- Consultant C (codex): codex/subthird_1_codex.out; full proof + experiments: codex/scratch/subthird_c/consultant_c.txt
  (and yield.csv, the *.py and prereg files there).
- Consultant G (Gemini): codex/subthird_1_gemini.out (§0–§1 up to "Surviving Adaptive Mechanism") and
  codex/subthird_1_gemini2.out (rest: M2 CG-2D-CRT, §2 kills incl. Lemma FT, M1, M3, literature, ranking).
  G's scripts: codex/scratch/subthird_g/{prereg.md,run_checks.py,results.json}.
(All paths relative to /Volumes/git/factoring-connection.)

You are one of the two; your job now is to try hard to BREAK the other's work and to rescue anything salvageable.
≤ 1200 words:
X1. The other's strongest positive claim(s): is each correct? For C: the linear-Kneser generator-height bound
    H ≥ x/(2⌈2n/s⌉+1) − O_k(m log n) for reduced type-D pairs (check the reduction steps: removing coincident sides and
    repeated subfield norms — is that WLOG? does Kneser apply in the needed generality? does it really cover the
    non-power, relative rank m ≥ 3 case at the L[1/4] balance?). For G: Lemma FT (is it new or Thue-inequality
    counting, e.g. Thunder/Bombieri–Schmidt?), and the kills K2b–K2d.
X2. Does C's height bound kill G's M1? If yes, what exactly survives ("sparse selection of exceptionally small norms
    inside a high-height module") and give the precise counting statement whose failure would still give α < 1/3.
X3. G's M2 (class-group 2D-CRT capacity constructor) and M3 (non-exact split pencils): is the assumed primitive
    (beyond-Johnson quadratic CRT list decoder; N-coherent non-exact pencil generator) simply a known open problem
    in the repo (pilist H_sel, crtdecode, nonexact O1/O2, sihard Thm J) or factoring/DL-equivalent? Does restricting
    to definite forms of disc −4N really block the navigate DLB / suffmap DLT reduction?
X4. Spot-check G's literature items (do the arXiv numbers exist and say that? use web search) — flag any you cannot
    find.
X5. Final: the surviving list (possibly empty) with P(α<1/3) and the single best experiment to run next.
