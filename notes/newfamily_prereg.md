# Pre-registration: newfamily battery (agent `newfamily`, phase 24, 2026-10-23)

This note, `src/newfamily/{recip,modsym,battery,harness,diag}.py` and `tests/test_newfamily.py` are committed
together **before** any instance of the `main` or `fresh` seed ranges below is generated. Before the freeze the only
runs were: unit tests (small primes; `pilot` range for the Salié test), a cost benchmark, and a **label-free** pilot
on 3 000 + 3 000 instances of the `pilot` range (n = 32, 64) that inspected only which statistics are constant. That
pilot caused two menu edits: `dmod3` → `dmod7` (6N·s(a,N) ≡ 0 mod 3 identically when 3 ∤ N, Rademacher) and
`ms_mod5` → `ms_mod3` (2[a/N]⁺ ≡ 2 mod 5 identically: the 5-Eisenstein congruence of 11a). No label was looked at.

## 1. Question, distribution, targets
* Same instance law and targets as the eps battery (`notes/eps_prereg.md` §1): `D_n`: N = pq from
  `balanced_semiprime(n, seed)`, conditioned on N ≡ 1 (mod 4). Target 1: ε = (−1/p) on all of `D_n`.
  Target 2: QR₂ = [2 is a square mod N] on the N ≡ 1 (mod 8) part.
* Question: do the surviving candidates of the newfamily theory pass (reciprocity / Euclid-computable invariants at
  N-dependent arguments; see `notes/newfamily_results.md` §1–2) predict a target beyond the N mod 2⁸ floor?
* Sizes n ∈ {32, 48, 64, 128}; 2·10⁵ instances of `D_n` per size (QR₂ sample ≈ 10⁵).

## 2. Menu (98 statistics; order = `MENU` in `src/newfamily/battery.py`, frozen)
All exact integers from `compute_stats(N)` (OPERATIONAL; N only; AST isolation test).
* **NC (10):** `nc_N_mod8`, `nc_N_mod256` (exactly degenerate under the stratified test), `nc_bit8`, `nc_bit9`,
  `nc_jacN_5` = (N/5), `nc_lead16` = top 16 bits, `nc_cubic7` = cubic class of N mod 7 (= (1+3ω/(N))₃ by Theorem RA),
  `nc_ded5_bits`, `nc_ded5_sign` (6N·s(5,N): a fixed-argument Dedekind sum, a function of N mod 5 and the size of N),
  `nc_ms_a3` = 2[3/N]⁺ (fixed-argument modular symbol).
* **Menu (84) = 6 arguments × 14 statistics.** Arguments a(N): `sqrt` = ⌊√N⌋; `gold` = ⌊N(√5−1)/2⌋;
  `cbrt2` = ⌊N^{2/3}⌋; `ferm` = (⌊√N⌋+1)² − N; `euler` = 2^{(N−1)/2} mod N; `hash` = 1 + (SHA-256(str N) mod (N−1)).
  Per argument: `jac` = (a/N); Dedekind D = 6N·s(a,N): `dsign`, `dbits` = bitlen|D|, `dmod7`, `dmod5`;
  CF of a/N = [0; a₁..a_r]: `len` = r, `n_a1`, `maxbits`, `altsign`/`altbits` of Σ(−1)^{i+1}a_i (the Rademacher /
  Barkan–Hickerson–Knuth part of 12 s(a,N)); Euclid remainders r_i of (N, a): `rem1mod4`, `rem0mod3`;
  plus modular symbol of 11a, m = 2[a/N]⁺ = 2λ({∞, a/N}) ∈ Z (normalised [0]⁺ = 1): `ms_sign`, `ms_mod3`.
  If gcd(a, N) > 1 the argument's statistics are 0 and the nonunit is recorded (rule 3).
  Known in advance to be constant (degenerate, can never fire): `ferm_jac` (a ≡ x² mod N), `euler_jac`
  ((2/N)^{(N−1)/2} = 1), and the pure-size statistics of `sqrt` (first partial quotient ≈ √N) and `cbrt2` at some n.
* **PC (4), one-sided, reused from `src/eps/cf.py` (proved in `notes/eps_prereg.md` §3):** `pc2_T1024`, `pc3_T1024`
  (fire ⇒ ε = +1), `pcq2_T1024`, `pcq3_T1024` (fire ⇒ QR₂ = 1).

## 3. Predictions (registered)
* P1 (all menu statistics, both targets, every size): **null**, |z| < 6. Reason: (i) the class-field symbols at
  N-independent arguments are aggregates (Theorem RA); (ii) at the N-dependent arguments used here, the CF/Euclid data
  of a/N give congruences a·y_k ≡ ±r_k (mod N) with y_k·r_k ≍ N/a_{k+1}, never small, so the relation mechanism that
  produced the eps n = 32 CF hits (small Q_i of √N) is absent; (iii) for `hash`, the law of CF(a/N) for a uniform
  public a is asymptotically independent of the factorisation (Baladi–Vallée-type equidistribution; heuristic).
* P2: NC statistics null; `nc_N_mod8`, `nc_N_mod256` exactly degenerate.
* P3: PCs 0 violations at every size; PC2/PC3 fire with predicted z ≥ 8 at n ≤ 64 (eps measured rates); at n = 128
  they are expected never to fire.

## 4. Seeds and runs
`SEED_BASE` in `src/newfamily/harness.py`: main = 24.6·10¹² + n·10⁹, fresh = 24.7·10¹² + n·10⁹, pilot = 24.5·10¹²
(+ n·10⁹ for the pilot). Disjoint from all earlier project ranges (eps: 9.1–9.5·10¹²). Command:
`PYTHONPATH=. .venv/bin/python -m src.newfamily.harness --n 32 48 64 128 --range main --count 200000 --workers 6`.
Memory: ≤ 2 GB RSS; one heavy job at a time.

## 5. Test and decision rules
* Test: the eps stratified permutation z (`src/eps/analysis.py::strat_z`), strata N mod 2⁸. Hit iff |z| ≥ 6.0
  (S = 100 Bonferroni, as eps). Replication: a fresh-seed run of the hit size, same sign and |z| ≥ 4.
* Validity per size × target: V1 planted control (ρ = 0.05) within ±3 of its expected z; V2 all PCs 0 violations and,
  for n ≤ 64, at least one PC with z ≥ 8; V3 within-stratum shuffle mean z² ∈ [0.6, 1.5]. At n = 128 only V1 and V3
  apply (declared now).
* Conclusion per size × target: "refuted for this menu at |ρ| ≥ ρ_MDE" (80% power at |z| = 6) if valid and no
  confirmed hit; otherwise "claim not refuted" for each confirmed hit, followed by the eps audit (recompute from N in a
  fresh process; finer floors N mod 2¹² and (N mod 2⁸, N mod 105); at n = 32 the two-way cluster-robust variance
  over p and q, because n = 32 reuses primes, see `notes/eps_results.md` §5.3). The registered conclusion uses
  n ∈ {48, 64, 128}; n = 32 is reported but its z are anti-conservative.
