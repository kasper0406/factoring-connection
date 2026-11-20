# chase cross-review (rules: codex/chase_common.md)

Outputs (relative to /Volumes/git/factoring-connection/codex):
- C-L2 (codex, lemma 2 prover): chase_L2_codex.out + scratch/ch_L2c/proof.txt (threshold inequality Q ≥ (x − D log
  C_kD)/g_h − (e/2) D m log m − O_k(Dm) for every h; "g_h = o(T) ⇒ Q ≥ x/T"; fresh-rank bound B^{o(1)} for candidates
  with o(m/log x) deep embeddings; 72-bit witness explained via its ternary unit representative; sweep 60–112 bits).
- C-L3 (codex, lemma 3): chase_L3_codex.out + scratch/ch_L3c/lemma3.txt (rank formula r = rank M − rank M_L; refutation
  of the rank-to-fresh-residual lemma by a coset-table DAG with zero fresh probes, cost N^{1+o(1)}).
- G-L1 (Gemini, lemma 1): chase_L1_gemini.out (truncated) + scratch/ch_L1g/full_proofs.md (Theorem 1 fiber identity /
  equivalence of branch selection and M2 witness search; Theorem 1(c) random-branch bounds; Theorem 2 two-choice
  weighted Theorem J for M2; proxy experiments).
- G-L2 (Gemini, lemma 2 adversary): chase_L2_gemini.out + scratch/ch_L2g/details.txt (constructions C1–C4; Rank-3
  Two-Hyperplane Relative Liouville Lemma; 72-bit witness covered).

Break the OTHER consultant's results (codex reviews G-L1 and G-L2; Gemini reviews C-L2 and C-L3), ≤ 1400 words. For
each claimed theorem: correct / fixable / wrong (specific step, counterexample if possible), and novelty. Priorities:
(1) C-L2's threshold inequality and fresh-rank bound — hypotheses, the Cramer/product-formula loss (e/2)m log m, the
"deep embedding" counting with exp(o(m)) choices, the KH application to each group; (2) G-L2's two-hyperplane lemma
(the det(B) = ‖v‖ construction, Hadamard bounds) and whether C2's 2^m − 2 orbits really have parity rank ≤ m;
(3) C-L3: is the coset-table DAG a genuine refutation of the model's R/credit rules, or does the model already charge
the table (N^{1+o(1)} work)? What minimal model amendment (e.g. charging modular-equality certificates) restores the
lemma, and is it natural? (4) G-L1 Theorems 1–2. End with: results to accept as proved, and the single best next
attack (or "the hatches are as closed as current tools allow").
