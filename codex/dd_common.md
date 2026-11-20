# dd: the divided-difference test for LP_exc^rel (common preamble)

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch only in your cwd).
Goal: factoring below L[1/3]. Read notes/chase_results.md (Lemma 2 section), experiments/chase/ch_L2c/proof.txt
(§1 eqs (1)–(2) and the "necessary combination" paragraph ~lines 42–51; §4 "Where the divided differences sit"
~lines 127–144), notes/subthird_results.md (Theorem KH, LP_exc^rel), notes/hatch_results.md, experiments/chase/ch_L2c/
sweep.py and experiments/masullo_ext/mext_c/cover.py (existing divided-difference fixtures over k = Q(√2), N = A²−2).

The test. Divided-difference modules: on each side P(α) = √2 (resp. Q(α′) = √2) with Q − P divisible by (T−b)·∏R_j,
pairwise coprime monic pole polynomials R_j of total degree p = O(m); generators 1, T and the quotients
(√2 − (P mod R_j)(α))/R_j(α). After clearing the common denominator ∏R_j(α), products satisfy s_h ≤ h(p+1)+1: exactly
the slow-growth regime where codex's threshold inequality (Q ≥ (x − D log C_kD)/g_h − (e/2)Dm log m − O_k(Dm)) does
not close. A failure of LP_exc^rel here REQUIRES log(M_i/t_i) ≥ b_k m log m at a positive fraction of the D embeddings
(near-maximal Cramer/conditioning loss), where M_i = max_j |L_i(w_j)| and t_i = least threshold at which ternary
vectors with |L_i(a)| ≤ t_i span k^m.

The question: does the shared denominator permit (or forbid) conditioning loss ≍ m log m at a positive fraction of
embeddings? Equivalently: for divided-difference modules, prove Σ_i log t_i − log(N𝔠N𝔠′) ≥ c_k x/T − C_kDm (confining
every capped candidate to D subspaces), or exhibit a family where it fails and count the resulting independent
small-norm candidates (the would-be sub-1/3 yield).

Output ≤ 1800 words, plain language, labels proved / heuristic / measured / speculative, file:line cites; long material
to a file in your cwd; pre-register experiments; exact or certified arithmetic. Final verdict: proved / refuted /
partial + exponent ledger + next step.
