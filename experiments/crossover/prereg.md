# crossover: hidden bits of a count versus its b-adic precision (pre-registration, 2026-11-27)

Question: for counts that factor N at full precision, how many bits about the factors does C(N) mod b^k carry
beyond public data, as k grows? Is there a count whose middle digits are hidden but natively cheap?

Families (N = pq balanced, p,q in [2^13, 2^14)): r2 (x²+y²), rQ2 (x²+2y²), R23 (principal form x²+xy+6y²),
r4 = 8σ(N), r8 = 16σ₃(N), KH = Σ_t H(4N−t²) = 2(N+q), τ(N) = τ(p)τ(q), and S(N), the floor8 short-square count
for x²+xy+6y², computed from N alone by modrec's hull evaluator. Ground truth for the others uses the factors.

Measure: collision (Rényi-2) conditional entropy Y_b(k) = −log₂ P[C ≡ C' mod b^k | same public cell], where the
public cell is N mod b^{k+1} plus the Jacobi symbols (−1/N), (2/N), (−3/N), (5/N), (−23/N). Bases b = 2, 3, 5.
Control: a public function of N (N² + 3N mod b^k), which must give Y = 0.

Hypotheses.
- H1: every family has Y = 0 up to a small public depth k0 (k0 ≤ 3 for b = 2). Its first hidden bit is a
  quadratic-character (Blum-type) bit.
- H2: Eisenstein families (r4, r8, KH) then gain about log₂ b bits per digit (digit-type). Cusp/class families
  (τ, R23, S) saturate or grow more slowly (character-type).
- H3 (the hope): some family has Y > 0 at a depth where a native evaluator is poly(log N). Expected false.
  Kill: none found.
