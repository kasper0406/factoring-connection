# hatch cross-review (rules: codex/hatch_common.md)

Round-1 outputs (relative to /Volumes/git/factoring-connection/codex):
- C-H1 (codex, hatch 1): hatch1_codex.out + scratch/hatch_h1c/report.md (Theorem J with per-depth q^k weights;
  folding obstruction at split primes; quadratic reconstruction measurements).
- C-H2 (codex, hatch 2): hatch2_codex.out + scratch/hatch_h2c/report.txt (determinant + KH obstruction: small norm ⇒
  one of ≤ D proper subspaces or cancellation at Ω(√x/log x) embeddings; Evertse/Evertse–Ferretti assessment;
  S-unit steering experiments 0 hits; exhaustive audit).
- G-H1 (Gemini, hatch 1): hatch1_gemini.out + scratch/hatch_h1g/hatch1_buium_crt_report.md (Theorems 1a–1e: Buium
  p-derivations, Int(Z) capacity wall; Theorem 2a beyond-Johnson Galois-folded O_K CRT decoding for orbit-respecting
  errors; Theorem 2b split-prime kill).
- G-H3 (Gemini, hatch 3): hatch3_gemini.out + hatch3_gemini2.out + scratch/hatch_h3g/hatch3_full_notes.md (Lemma
  PIN-1, Prop THIN-LEDGER, Theorem DAG-POOL removing Assumption T / H_cusp in the shared regime, Theorem
  SINGLE-N-REUSE, Theorem TIED-PARITY Reed–Muller compression, candidate C-RM / ML-RM).

YOUR TASK: break the OTHER consultant's results (codex reviews G-H1 and G-H3; Gemini reviews C-H1 and C-H2), ≤ 1500
words. For each claimed theorem: correct / fixable / wrong, with the specific step and a counterexample where possible.
Prioritise: (1) Theorem DAG-POOL (it would repair a gap in notes/lowerbound_results.md: check disc scaling
Q^{(D−1)(D−2)}, the D^{ω(Q)} sublattice multiplicity count, the pooled profile, and the cusp-budget condition);
(2) Theorem TIED-PARITY and whether C-RM is a genuine new route or reduces to LP_exc^rel; (3) C-H2's determinant
obstruction (statement, the "≤ D proper subspaces" claim, and the √x/log x embedding count); (4) G-H1 Theorem 2a novelty
(is beyond-Johnson decoding of Galois-folded number-field CRT codes known? cf. Guruswami 2003, BHKS) and Theorems 1c–1d;
(5) C-H1's weighted Theorem J. End with: the list of results you would accept into the repo as proved, and the single
best next attack on any hatch.
