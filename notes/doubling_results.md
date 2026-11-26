# doubling: fast evaluation of the central-binomial zero test mod N (builder role) (2026-11-26)

**Process.** GPT-6.2-astra via codex (session 01a21793-c085-7b80-b28a-547d48d7fdf5). Brief: `codex/doubling.md`.
Output: `codex/doubling_codex.out`. Prereg: `notes/doubling_prereg.md`. Code, checks and ledger:
`experiments/doubling/`.

**Verdict.** No poly(log N) evaluator was found, and nothing here is sub-L[1/3]. The best new algorithm is a
factor-free T0 decider in Õ(N^{1/4}). That is no faster than factoring and reading off the class (Pollard–Strassen
is N^{1/4}; randomized factoring is subexponential).

- **D1 CM truncation (proved).** T_K = Σ_{j<K} C(4j,2j)C(2j,j)·72^{−j} mod N, with K = ⌊√N/2⌋ + 1. This is the
  truncation of ₂F₁(1/4, 3/4; 1; 8/9).
  - Mod each s | N, 3^{(s−1)/2}·T_K is the Hasse coefficient of y² = x(x+1)(x+2). That is the curve x³ − x
    shifted by 1, so for balanced N: T_K ≡ 0 ⟺ p ≡ q ≡ 3 mod 4, and T_K is a unit otherwise.
  - The Lucas cross term from the earlier B(N) formula drops out.
  - Orchestrator check: 4,456/4,456 balanced semiprimes (q < 2p) with 0 errors. The test fails on unbalanced inputs
    when the cutoff passes p/2; codex's ratio form then catches p through a denominator gcd.
  - BGS block products give Õ(√K) = Õ(N^{1/4}).
- **B(N) Jacobi formula (proved).** With d = (q−p)/4: (B/N) = ((q−p−1)!!/N), checked on 384 N. The
  Gaussian coordinates cancel, but d determines the factors.
- **Eisenstein candidate (refuted).** The E₄-coefficient of E_{N−1} fails at N = 397·2,755,601, where it is
  zero although B is a unit.
- **D2 Gaussian binomial (proved).** gcd([2m choose m]_X, (X^N−1)/(X−1)) is 1 in class (1,1) and Φ_p·Φ_q in class
  (3,3), so evaluation at X = 1 gives 1 or N. The degree gives p + q. Fast evaluation is open.
  - Obstruction (proved): binary sectioning of the carry-support series needs Ω(2^j) gates after j sections.
  - Reflecting the recurrence erases the midpoint distinction.
- **Calibration (checked).** An arbitrary-index nonunit oracle factors N. So does a zero-mod-N oracle for equal-size
  factors, with the first joint zero at (q+1)/2.
- **D3 amortisation (measured).** One product/remainder tree over j ≡ 3 mod 4 in [L, 2L) classified 32,768 of
  32,768 32-bit inputs. That is polylog per member only for batches of about √N targets. A single target pays
  Θ(L log L).
- **D4 (measured, null).** 768 exact residues; 0 of 48 N-only predicates passed; no hidden-coordinate signal
  (family-wise p = 0.98, 0.27, 0.14).
