You are a research mathematician (analytic number theory: smooth numbers, Bombieri–Vinogradov / large sieve for smooth numbers, Drappeau/Harper; rigorous index calculus). Work READ-ONLY on /Volumes/git/factoring-connection (your cwd is a writable scratch dir). Use web search where useful. This is a BLIND round: work the problem yourself; do not look for other agents' answers to this exact task (there are none yet; files notes/uprime_* and notes/pc3_* are off limits even if they appear).

CONTEXT (read these):
- notes/multidescent_results.md, especially §3.2 (design 2L, Prop UB) and §3.3 (Hyp U′ and its four re-proof items), and the review tables (§8).
- notes/lv5unif_results.md in full (Theorem U: Lemma CS, GM-U, FM-U, Lemma 8.5′(R′) §5 with (R′1)–(R′4), W-U, DS-U, L2-U, LD*-U, D-U, Lemma IR, §10 recursion).
- paper2/main.tex: def:good, lem:goodmod, thm:FM, lem:DS, thm:L2, thm:LD, and Appendix A (proof of Lemma 8.5′: Lemmas E1–E8, LS, cor:LS, SM).
- notes/auditU_audit.md; notes/recurse_results.md §1.3 (Lemmas OS, W, F); notes/escapeR2_results.md §7.3 and §8 step 8.

TASK: Hyp U′ of multidescent. Two-level design at α = 3/8: d ≍ D^{1/4}, k ≍ D^{3/8}, B = L[3/8, a] (factor base), B₂ = exp(b₂ D^{1/2} log log N) (child bound), D = log N/log log N. Acc_ℓ(P): rational cofactor G(P)/ℓ B-smooth; algebraic value F(P) = s·q₁⋯q_m with s B-smooth and q_i DISTINCT primes in (B, B₂] not dividing b·D_A. Children 𝔮_i = (q_i, a/b mod q_i), child box = points of Λ_𝔮 = {P′: det(P,P′) ≡ 0 mod q} in a box of radius S₂ (to be specified), Good(P′): q ∥ F(P′), G(P′) and F(P′)/q B-smooth.

Decide, with proofs or exact failing estimates:
1. Hybrid restricted BV. Lemma 8.5′(R′) uses one y. With y = B₂ the modulus family is too large (error e^{−cD^{3/8}/log²D} vs ϱ(z,B₂)/ϱ(z,B) = e^{Θ(D^{3/8} log D)}); with y = B the algebraic u ≍ D^{1/2} violates (R′3). State and prove (or give a precise proof plan and the exact gap of) a hybrid restricted BV: moduli t = ℓ·s with s B-smooth (y₁ = B), values B₂-smooth (y₂ = B₂), error normalised by the density of the ACTUAL modulus family. Go through Appendix A of paper2 step by step (E3, E8, Fact D = Drappeau Lemme 5, dBT, LS, cor:LS, SM, Steps D/E) and say exactly which step changes and what new condition appears. Then: do GM-U, FM-U (good/bad mass), L2-U, DS-U, LD*-U, D-U go through with two smoothness parameters? Where exactly?
2. Mixed (B, B₂) scan replacing Lemma CS: existence, the size of S₁ (exponent σ₁ in terms of a, κ, d₀, b₂), m₁ = S₁²ρ_r ρ_F bounds; and an explicit child radius S₂ (possibly depending on the window of q) with the child main term m₂ ≥ (log N)^{C₂}.
3. Squarefree large part / q ∤ D_A: union bound Pr_c[p² | F(P)] and Pr_c[q | F(P) ∧ q | D_A]: constants, and the resulting constraint on (a, b₂, κ, d₀). Or the AP mod t·p² alternative.
4. Recursive value factorisations up to B₂ (cost; Lemmas OS/W; the bounded-depth search recurrence) and Lemma IR at γ₂.
5. The bullet "E[N_𝔮 | E] large enough for PC₃": what exactly can be proved at polylog loss for the CONDITIONAL child mean (given the parent's acceptance and child), what is a genuinely new 2-point lower bound, and whether that belongs to U′ or to PC₃.
6. If U′ fails at 3/8, which α ∈ (3/8, 2/5) can it deliver?

OUTPUT (≤ 4000 words): verdict first (U′ PROVED / PROVED EXCEPT X / FAILS AT Y), then for each item a precise statement, proof or proof sketch with the exact conditions, and the gaps. Flag every place where you are unsure.
