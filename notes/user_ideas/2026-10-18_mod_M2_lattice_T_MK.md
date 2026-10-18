# User-supplied derivation (2026-10-18): why M (mod-M² lock) and K (Taylor centers) collapse to T = MK = N^{1/5}

- Weighted sums a·p + b·q are fixed mod M² exactly for (a,b) in the index-M sublattice Λ_{r,s} = {a r − b s ≡ 0 (mod M)} ⇔ a r² ≡ b N (mod M).
- Restricting the weights to Λ inflates the coefficient size to B ~ (MK)^{1/2}. With K Taylor centers of width ε = N^{1/2}/K, the second-order remainder is R ~ Bε²/N^{1/2} = M^{1/2} K^{−3/2} N^{1/2}. Steps per branch: S = R/M² = N^{1/2}/(MK)^{3/2}.
- Total: L_baby + MK + MK·S/L_baby. Balancing gives T^{5/2} = N^{1/2}, i.e. T = MK = N^{1/5} (Harvey's exponent).
- "Crack" attempt: decouple r from (a,b). The user's own conclusion is that it fails: each window needs M pairs. Restricting to the sublattice gives an M^{1/2} net gain, the same as K^{1/2}.

Orchestrator generalization: if a scheme gives per-branch steps S ~ N^{1/2}/T^{δ}, the balance gives T = N^{1/(2(1+δ))}. Here δ = 3/2 → 1/5. Cancelling the second-order term with no coefficient blow-up (δ = 5/2) would give N^{1/7}. This is exactly the test for cyclotomic Frobenius handles; codex r7 found that the naive cubic widens the window instead.
