# pcest_1 — blind research round: the pair correlation PC for bounded-k Dixon descent

You are a research analytic number theorist. The repo at /Volumes/git/factoring-connection is read-only; use cwd for scratch.
Read first: notes/hypD_results.md (all of it; especially §3 Theorem K and Corollary K, §4 Prop PZ, §4.1′ Lemma 8.5′ and FM,
§4.2 PC_rat / PC_alg / Prop G / RBDH′, §7.4 Prop P) and notes/lemmaM_results.md §1–§3.

GOAL. Prove PC (the second-moment input of Prop PZ) at precision L^{o(1)}, L = L_N[1/2,1], uniformly in N, for the fixed-k
design (a′_k), or find a modified design (same or better constant c = (1+θ)/√2) in which it is provable from published theorems.
Only an UPPER bound E[N(N−1)] ≤ L^{o(1)}·(E N)² is needed, averaged over the field A, the prime ℓ ∈ J ≍ N^{1/(2k)} and all
pairs of box points. The precision needed is extremely weak: L^{o(1)} = exp(o((log N log log N)^{1/2})).

Label every claim PROVED / KNOWN (exact citation, exact ranges of x, y, moduli, residues) / CONJECTURED / FALSE.

Q1 (PC_rat, the binary problem). With the rational side x − m, m uniform, the rational part of E[N(N−1)] is an average over
 b, b′ ≤ √ℓ S, h = Δ/ℓ with |h| ≤ 2S², of #{j : s₀ − bj and s₀′ − b′j both B-smooth}, B = L[1/2, a], values N^{Θ(1)}.
 (a) Is any upper bound within L^{o(1)} of the product of densities known or provable, at least on average over (b, b′, h)
     with the actual weights? Consider Rankin's trick, upper-bound sieves, Shiu / Nair–Tenenbaum bounds, Hildebrand's work on
     smooth pairs, Matomäki–Radziwiłł–Tao-type averaged correlations, de la Bretèche–Drappeau / Pascadi. Say exactly why each
     works or fails (density L^{−c} vs savings).
 (b) Randomize the rational side: g = g₁x − g₀ with (g₀, g₁) random (Kleinjung-style, Res(f, g) = N, f = f̂_g + g·R). Check
     Kleinjung's residue constraints on g₀ mod g₁ (or a generalized base-(g₀/g₁) expansion), the coefficient and norm sizes, and
     whether θ and c change. Does the 2D rational randomness turn PC_rat into a provable statement? If so, from which published
     theorem exactly (state its hypotheses and check them at our parameters, including the moduli range and whether intervals
     are full-length or short)? Watch for hidden losses: lattice index ℓ | Δ, irregular lattices, small q, non-coprime classes,
     short-interval effects.
Q2 (PC_alg). Given the rational data, PC_alg = Pr_c[F, F′ both B-smooth]/(Pr_c[F]Pr_c[F′]) with F = F̂ + G⟨c, v⟩,
 F′ = F̂′ + G′⟨c, v′⟩, G, G′ the rational values. The joint lattice has index |GG′Δ|. Prove an upper bound within L^{o(1)}
 on average over the pair family under the rational-smoothness weights, or reduce it to the weakest clean statement you can,
 and say whether that statement follows from Harper's BV/BDH (arXiv 1208.5992), LV Lemma 8.5 (arXiv 1805.08873) or any
 other published result. Note codex's earlier warning that a 2D LV kernel lemma needs coefficient height H ≫ V²; cost it.
Q3 (alternatives). Two randomized algebraic sides; restricting which pairs are counted; weighted counts; any design that
 makes PZ unnecessary. For each, the exact constant c and the single missing estimate.
Q4 (literature). Has anyone proved a rigorous factoring algorithm for all N with L[1/2, c], c < 1, 1992–2026? Any rigorous
 (not heuristic) second-moment / pair-correlation results for smooth values in NFS-like settings?

Be concrete; give formulas and parameter checks. Do not pad. Final answer ≤ 4000 words.
