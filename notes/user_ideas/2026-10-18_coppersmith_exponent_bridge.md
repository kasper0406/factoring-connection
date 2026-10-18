# User-supplied results from a parallel session (2026-10-18): Coppersmith lattice ↔ Harvey/Hittmeir exponents

1. **The bridge.** For H in Coppersmith's shift lattice L_{m,d} spanned by g_{i,j}(x) = x^i (x+p0)^j N^{m−j}, with x0 = p − p0:
   - k = H(x0)/p^m = Σ c_ij x0^i q^{m−j} is an integer.
   - Mod p−1 we have H(1−p0) ≡ k, hence a^{H(1−p0)} ≡ a^k (mod p).
   - At m = 1 with H(x) = (x+p0)² − S0(x+p0) + N: k = p+q−S0 and H(1−p0) − k = (p−1)(q−1). This is exactly the Fermat/Lehman search a^{N+1} = a^{p+q}.
   - Orchestrator checked the integrality, the mod-(p−1) reduction and the m = 1 example by hand.
2. **Phase transition (claimed):** the Minkowski bound gives log|k|/log N = (β − 1/4)·m + (1−2β)/4 at d = 2m, with X = N^β. For β > 1/4 a higher multiplicity m inflates |k|, so m = 1 (order-2 cancellation) is optimal. NOT yet independently verified.
3. **Exponent vs base ring (claimed):** p² − S0 p + N = k·p. Mod p the factor p is a zero divisor and k is annihilated, so Strassen-type base-ring evaluation is stuck at N^{1/4}. Mod p∓1 the factor p is a unit, which exposes the compressed k ~ X²/√N and gives N^{1/5}. A heuristic explanation.
