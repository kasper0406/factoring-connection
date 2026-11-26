# buildbit: constructive (builder-role) attack on the hidden-class bit (2026-11-26)

**Process.** GPT-6.2-astra via codex in the builder role (session 01a21773-2898-7433-9eac-dd535aa3988d). Brief:
`codex/buildbit.md`. Output: `codex/buildbit_codex.out`. Prereg: `notes/buildbit_prereg.md`. Ledger, code and
counts: `experiments/buildbit/` (`ledger.md`). The user is ~90% sure the bit has a poly(log N) algorithm.

**Verdict.** Nothing is sub-L[1/3] and no public statistic passed. The best construction is an exact
reformulation. A poly-time binomial zero test would factor N. Evaluating that test is the factorials-mod-N wall.

- **R4, binomial carry test (proved; the orchestrator checked T0 on 183 N and T1 on 168 N with 0 errors).** Put
  m = ⌊hN/d⌋.
  - T0: N | C(2m, m) with d = 4, h = 1 ⟺ p ≡ q ≡ 3 mod 4.
  - T1: d = 16, with h = 1 or 3 by N mod 16. Then N | C(2m, m) ⟺ r ≡ 9 mod 16.
  - T2: d = 7 gives a unit, zero or a proper gcd.
  - The mechanism is Kummer: a carry in the low base-p digit.
  - **Lift (proved).** A uniform zero-test oracle over d = 2^j, with h found by bisection on powers of N/a² mod d,
    recovers r in O(log N) queries. A factor-backed oracle did so on 768/768 inputs (40–80 bits).
  - Folding the coefficient mod X^L − 1 (L = 8, 16, 32) is not a public substitute: 2,304 of 2,304 folds gave units.
  - Evaluating C(2m, m) mod N at m ~ N is a factorial-mod-N problem. The best known cost is about N^{1/2}
    (Bostan–Gaudry–Schost). The open next step codex names, a dyadic doubling rule for the zero image of
    2nU_n = (2n−1)U_{n−1}, would be a poly-time factorial-type evaluation.
- **R1, pooled CF(√kN) over 12 multipliers (measured).** This gives public square-root-of-−1 certificates via
  F₂ relations. Over 40/52/64/80 bits the advantage was 4.7/0.32/0.037/0.012%, decaying like N^{−0.215}
  (99% CI [0.12, 0.52]). Pooling added 144/1/0/0 certificates beyond the union of single multipliers. The limit is
  the supply of odd-sign square relations (families M/Q).
- **R2, explicit ambiguous class (proved; 125 small class groups checked).** The form f = (N, N, (5N−1)/4) has
  D = −N(4N−1), and [f] is a square in Cl(D) ⟺ T0 is in the 1,1 case. What remains is square membership for one
  public class, i.e. class-group structure, currently L[1/2].
- **R3, scalar quartic for torus squares (proved).** b = a³/Norm(a) and H(X) = X⁴ − 2TX² − 8X + T² − 4S. A root
  exists locally ⟺ b is a torus square, which is always true inert and true 1/16 of the time split. The resultant
  Res(H, X^N − X) was a unit on all 3,072 inputs. It tests an M-family order condition, and its Jacobi sign was null.
- **Gap.** Every route ends at deciding one public object without the factorization: a binomial zero mod N, a class
  square, or a root of a quartic mod N.
