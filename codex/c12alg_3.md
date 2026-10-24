You are an adversarial referee (analytic number theory: smooth numbers, BDH/BV for smooth numbers, short intervals, Hildebrand–Tenenbaum; lattice counting). Work READ-ONLY on /Volumes/git/factoring-connection (cwd is a writable scratch dir). Use web search. Rank issues FATAL / MAJOR / MINOR with location and fix.

READ notes/c12alg_results.md v2.1. Your previous adversarial report on v1 is codex/c12alg_2.out; check that its 13 issues are correctly repaired in v2 (§10 table) — in particular Lemma η's new pointwise proof via Hildebrand at B₂ (§1.2), the PP short-fibre truncation and the parent truncation at B₂^{1/2} (§3.3), Lemma PU's repairs, the same-ℓ parent pairs (§3.4(ii)), Prop CLB without angular discard (§4), the exact identity w/z = Pr(R)/#J (§5).

THEN ATTACK the new §3.6 (Theorem AB₂, the cell-covariance device, adapted from c12rat's Prop PD in notes/c12rat_results.md §8.1, an unreviewed working-tree note) and Theorem 2L-DW₂ (§5: conditional L[3/8]):
1. Step 1: is the in-cell factorisation exact as claimed (principal term N(c)d₁d₂ for the class δ; singles exact with no character error; δ common to both coordinates since κ is a unit)?
2. Step 2: is A₁₂ − A₁A₂ = Cov_ν̂(d₁,d₂) + errors correct, including the independence of δ from the rectangle, the boundary cells (Cauchy–Schwarz against fourth moments), the mixture over y (Cov_y term), the TV and the w-mixture?
3. Step 3: is Var(d_i)/(E d_i)² ≤ C(D^{−2α}λ^C + λ^C/u²) right — i.e. are all first-order quantities entering only through variances? Is there any first-order term entering linearly (e.g. normalisation by m_𝔮², the 1/m term, the root conditioning, environment cells, the q | Δ′ case, the parent acceptability conditions, the SQ subtraction at 3/8 where it is borderline)?
4. Step 4 and the parameters K, K₀, T at second order; Hildebrand's range with K = K₀⁴D^{8α}(log N)^C.
5. Is the claim that, given RAT₂ and AB₂, the budget closes on [3/8, 2/5) (A ≤ D^{−2α}λ^C, √(AB) ≪ η) correct, and are all other inputs available at α = 3/8 (uprime's U′ at 3/8 with DW; (SQ) constraint; η₀; margins `params.margins2`)?
6. Overclaims in §0, §5.
OUTPUT (≤ 2500 words): recommendation, ranked issues, probabilities that (a) Theorem AB₂ holds as sketched, (b) Theorem 2L-DW₂ is a correct conditional theorem.
