# pcest_2 — adversarial referee: notes/pcest_results.md

You are a hostile but fair referee (analytic number theory + rigorous factoring). The repo at /Volumes/git/factoring-connection
is read-only; cwd is scratch. Read notes/pcest_results.md IN FULL, plus the context it builds on: notes/hypD_results.md (§1–§4,
Theorem K, Prop PZ, Lemma 8.5′, PC), notes/lemmaM_results.md §1–§3, and LV arXiv 1805.08873 §5 and §8 (Lemmas 5.13, 5.18, 5.20,
8.5, Facts 3.13–3.17, 8.8–8.12, Claims 8.14, 8.16) and Harper arXiv 1208.5992 Theorems 1–2 / Propositions 1–2.

For EVERY claim labelled proved (Lemma GB, Lemma 1.2, Lemma 1.3, Lemma 8.5′, Theorem FM Steps 1–4, Lemma Q, Lemma DS, Theorem R,
Lemma Kern (i)–(iii) incl. det K = ‖v∧v′‖/|Δ| and the compactness bound, Prop 4.2, Theorem 4.3, Prop B5, Prop Z, the cost
arithmetic of Theorem A and A′, the minimal-k numbers) give a verdict CORRECT / CORRECT WITH FIXES / WRONG with the exact reason.

Check in particular:
1. Uniformity in N of every constant; whether any step hides an L[1/2]-size loss (e.g. N^{o(1)} vs L^{o(1)}, divisor functions,
   2^{ω}, number of grid points, log ω with ω = L^C, the choice S² = L^ε/(ρ_r ρ_F) vs Theorem R's S² ≥ ρ_r^{−2}).
2. The exact ranges of every cited theorem (LV Lemma 8.5: x^ε ≤ Q ≤ √Ψ, ω → ∞, ω = y^{O(1)}; Fact 8.11: 0 < η < 1/80, x^η ≤ Q ≤ √x;
   Hildebrand short intervals z ≤ y^{5/12}; Tenenbaum ω(q) ≤ y^c/log(1+u); Claim 8.16 needs u → ∞), and whether LV's proof of
   Lemma 8.5 really bounds the majorant W (8.2) for an arbitrary sub-window, as Part 1 of Lemma 8.5′ uses.
3. Theorem FM Step 2: does the LV 5.13/5.18/5.20 transplant hold with general signs, |a| > b, the truncation ω = L^C, and the
   residue-free grid 𝒳? Is the lower bound for the uniform pieces correct (the |F̂| boundary terms, the sign split)?
4. Theorem R: Lemma Q's index and graph form; the covering of the thin parallelogram Π by squares of side M/V (is the width
   bound right, does the angle restriction from b, b′ ∈ (V/2, V] hold); the near-axis band; Lemma DS's choice of D and the
   non-reduced residues; the pair count via Lemma 1.3; is anything lost when the candidate set is restricted to prime b
   (does PZ still apply, does Theorem K still apply).
5. Lemma Kern (iii) and Prop 4.2: the fibre count, H ≥ V², and whether V²/|Δ| (the bounding-box factor) is correctly averaged.
6. Theorem 4.3 and Theorem A: is LD (as stated) really sufficient together with the proved pieces and hypD's Theorem K, DD,
   field selection, LP shape test? Is anything else silently assumed (LV6, HECM, irreducibility transfer, the g₁ column)?
7. Is Hypothesis LD perhaps provable from known results after all (or false)? Is Prop B5 a fair barrier statement?
8. Literature: is the divisor-switching upper bound (Lemma DS) classical (give the reference), is the 2D-randomized rational
   side in a rigorous analysis new, and has anyone beaten c = 1 rigorously 1992–2026?

Be concrete. Do not pad. ≤ 3500 words.
