# User-supplied results (parallel session, 2026-10-18; scripts not available here): Gray-code Harvey implementation and audit

- CRT fiber power tables: Δ_i(a,b) = k_i·M_i² with |k_i| < m_i²/2, so step multipliers come from global tables, with no per-center exponentiations. Zero-carry theorem: carries are absorbed into the baby-step index. Base exponent compressed by q_W = ⌊W/M²⌋.
- 60-bit benchmark, M = 8·3·5·7·11 = 9240:
  - lattice reductions 1920 → 30;
  - base exponentiations 19200 → 300;
  - exponentiation multiplications ×293 fewer;
  - total multiplications ×18.5 fewer;
  - wall clock ×48.7.
- Their own audit: Stage 2 (setup: lattices, square roots, exponentiations) gains 2^{t+2} and Θ(log N). Stage 3 (the subproduct tree over T_pts = φ(M)K leaves) is unchanged except for a constant 2–4× from a mod-64 identity: r0 ↦ r0+4 shifts C0 by ±16 (mod 64), giving 2 buckets. For odd m_i, bucketing or affine-sublattice constraints cost more than they save.
- Orchestrator assessment:
  - Asymptotically, Stage 3 at Θ(T log² N) already dominates Stage 2 at Θ(T log N), so the total gain is at most a log-type factor, not exp(Ω(log N/log log N)).
  - The ×293/×48 figures need a fair baseline. An optimized Harvey implementation walks lattice points incrementally, at one multiplication per giant point, rather than doing T full exponentiations.
  - The odd-prime bucketing failure matches cycloharvey's guess-neutrality.
  - An asymptotic gain requires fewer distinct giant LEAVES, which is exactly the several-moduli gap (agent multimod).
