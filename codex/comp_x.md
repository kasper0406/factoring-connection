# comp cross-check (rules: codex/comp_common.md)
- C (codex): comp_codex.out + scratch/dd_c/compensation/{split_poles.txt,proof.txt,results.txt} — split-pole
  compensation theorem Q ≥ (m−6)x/n − 2nm log m − pH/(n−d_*) − O(Dm) ⇒ Q ≥ (2−o(1))x/T for integer-split poles with
  integer remainders; resultant/endpoint identity for general poles; common-rank sufficient hypothesis.
- G (Gemini): comp_gemini.out + scratch/dd_g/comp_proofs.txt — Submodule-KH Compensation Lemma (shared ternary basis
  on D−o(D) embeddings ⇒ Q ≥ (1−o(1))x/T − O_k(Dm)); two-place covolume obstruction; support-width lemma; general-pole
  divided-difference identity; families M1–M5 with N_ind = 0 in 55/55 boxes; non-scalar norm-1 units in M2_mix.
(Paths relative to /Volumes/git/factoring-connection/codex.) Check the OTHER consultant's theorems hard (codex checks
G; Gemini checks C), ≤ 900 words: correct / fixable / wrong with the specific step; whether the two remaining open
regimes ("overlapping small pole factors at one embedding" vs "distinct unshared bases at many embeddings") are the
same; and whether combining C's split-pole theorem with G's Submodule-KH lemma already closes a larger class. End
with the precise single remaining statement.
