You are an expert referee (computational number theory, deterministic integer factorisation: Harvey 2021, Harvey–Hittmeir 2022/2026, McKee–Pinch, Lehman). Read-only access to /Volumes/git/factoring-connection; write scratch scripts in your current directory.

Read notes/dchworst_results.md §§1–6 (draft), the code src/dchworst/wcalg.py, and for background notes/dcltheorem_results.md §§2–4 (Theorem D′, lemmas 1–10) — Theorem R reuses D′'s targets and class lock.

Your job: try hard to BREAK Theorem R ("DCL-wc factors every balanced N = pq whose pair is covered by at least one block, in cost Õ(Σ(m+T') + B0 T^{1/2} + √N/T + √N/Ψ) with no clean-band / order-mod-q / FULL hypotheses"). Check in particular:
1. Lemma O (one-sided order dichotomy, smooth-subgroup case Λ ≥ Ψ(p−1,B0)); whether it is genuinely new relative to Harvey–Hittmeir arXiv 2601.11131 (cite their lemma numbers) and whether the cost claims hold.
2. Lemma G and the equal-order split (factored multiple E = L·o, valuations argument, D > T because b was selected by Lemma O).
3. Lemma X: uniqueness of p-roots given ord_p(γ) ≥ m, the exploration induction for the true giant, the exact-leaf identification S = Bp+Aq, the batched cost bound Õ(m + #N-giants) per level.
4. Signed weights (A<0) and the diagonal weight (1,1): class lock, window bounds, giants.
5. Corollaries R1 (Minkowski ⇒ N^{1/5} worst case; is this just Harvey's/McKee–Pinch's?), R2 (DCH_pair ⇒ N^{1/6+o(1)} worst case), R3 (D′ exceptional set reduced to uncovered pairs).
6. Theorem N (Brun–Titchmarsh union bound) and T2; the barrier statements B1–B3 in §5 (are they correctly stated and labelled? any known result that gives worst-case coverage with c>0 that we missed?).
Also run the code on any adversarial N you can think of (e.g. p,q with equal orders of all small bases, p−1 and q−1 sharing band primes, tiny ord_p(γ)); report any N where a covered pair fails.

Output Markdown: verdict first (ACCEPT / ACCEPT WITH FIXES / REJECT for Theorem R), then numbered findings with severity, then concrete fixes. Label claims PROVED / UNVERIFIED / REFUTED.
