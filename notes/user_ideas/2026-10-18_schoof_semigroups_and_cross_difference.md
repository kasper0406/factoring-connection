# User-supplied results from a parallel session (2026-10-18; scripts not available here)

1. **Two-semigroup obstruction (Schoof / ψ_ℓ):**
   - Group-law powers [N]_E mod ψ_ℓ collapse to [N mod ℓ]_E.
   - Coordinate powers X^N = (X^p)^q apply the q-th power in the coordinate ring. That ring's unit group has order p^{d_i} − 1, so it is not an endomorphism composition.
   - Symmetric functions (traces, resultants) erase the p-Frobenius permutation: Sym(h(x_i)^N) ≡ Sym(h(x_i)^q) (mod p).
   - Orchestrator note: this matches Track B/D B1 (F_N ≡ F_q^p, q mod L_α dependence) and the primality "exponent menu". "Single semigroup" is a model assumption: Schoof itself mixes, and at Elkies primes a guessed eigenvalue λ makes X^p ≡ x([λ]P) mod the kernel polynomial computable. That is an elliptic analogue of the cyclotomic σ_r handle.
2. **Linearized 2D Lehman and the cross-difference law:**
   - Exact identity, checked by hand by the orchestrator: (a0+u)p + (b0+v)q − (2√(a0b0N) + √(b0N/a0)u + √(a0N/b0)v) = (√(a0p) − √(b0q))·((a0+u)√p/√a0 − (b0+v)√q/√b0).
   - Linearization trades Δ_fine² for Δ_coarse·Δ_fine.
   - The Hessian of 2√(abN) is rank 1, and its flat direction is radial, so it leaves a/b unchanged.
   - Multipoint/BSGS collision tests need separable exponents, E(j,u) = f(j) − h(u), i.e. a vanishing cross-difference. Decoupling the residue r from weight index i forces a_i ≡ a0 and b_i ≡ b0 (mod M), and the resulting M-inflation cancels the M² gain. This supports N^{1/5} robustness for separable schemes.
