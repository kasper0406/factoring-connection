# Pre-registration: classE battery (agent `classE`, phase 25, 2026-10-23)

This note, `src/classE/{structure,fd,battery,harness,toy3,diag}.py` and `tests/test_classE.py` are committed
together **before** any instance of the `main` or `fresh` seed ranges, or any toy population, is labelled.

Runs before the freeze:
- unit tests (random small moduli; no target labels);
- exploratory identity checks on random moduli (scratch; no labels): the d = 3 congruences, Prop VIS, the
  √N decomposition;
- a **label-free** pilot on 3 000 instances of the `pilot` range at n = 48 and n = 128. It read only which
  statistics are constant or are functions of the floor strata, plus positive-control firing counts. One menu edit
  followed: `cf3_Q1mod4_odd` was constant and was replaced by `cf3_Qsq` and `cf3_Q1mod5`.

No target label was looked at.

## 1. Question, distribution, targets, floors
- **Instance law D_n.** As in eps and newfamily: N = pq from `balanced_semiprime(n, seed)`, where p < q are n/2-bit
  primes with the top two bits set, conditioned on N ≡ 1 (mod 4).
- **Targets.**
  - **ε** = (−1/p) on all of D_n.
  - **tb** = [8p ≥ 7q]. Since p/q ∈ (3/4, 1), this is the first non-constant bit of p/q.
- **Floors.**
  - ε: strata N mod 2⁸, as in eps and newfamily.
  - tb: strata = the top 12 bits of N. This is the registered primary floor.
    - The target is archimedean. The size of N predicts it with constant advantage (leakmap §0 item 5), so a floor
      of N mod 2⁸ alone would be met trivially by any size statistic.
    - tb against N mod 2⁸ strata is also computed. It is reported as a secondary, size-confounded view and is not
      used for any conclusion.
- **Question.** Do class-E functionals at structured arguments predict a target beyond its floor? The same question
  is asked of shifted d = 2 lattice-cone (Fourier–Dedekind) sums, reduced-form walk coefficients of disc −4N, and
  the infrastructure path of disc 12N.
- **Sizes.** n ∈ {48, 64, 128}, with 2·10⁵ instances per size.

## 2. Menu (130 statistics; order = `MENU` in `src/classE/battery.py`, frozen)
All are exact integers from `compute_stats(N)`, which is OPERATIONAL (it sees N only; AST isolation is tested).

- **NC (5), controls.** `nc_N_mod8`, `nc_N_mod256` (degenerate under the ε floor), `nc_bit8`, `nc_lead12`
  (degenerate under the tb floor), `nc_lead16`.
- **ARG (8 × 13 = 104), role `menu`.**
  - Arguments a(N):
    - `isqrt` = ⌊√N⌋;
    - `isqrt2` = ⌊√(2N)⌋;
    - `cbrt2` = ⌊N^{2/3}⌋;
    - `div3` = ⌊N/3⌋ − 1;
    - `fermres` = 2^{N−1} mod N;
    - `euler3` = 3^{(N−1)/2} mod N;
    - `invsqrt` = ⌊√N⌋⁻¹ mod N;
    - `hash2` = 1 + SHA-256("classE|N") mod (N − 1), a pseudo-random control.
  - Per argument, with a/N = [0; a₁..a_r]:
    - `len` = r, `n_a1`, `maxbits`;
    - `altsign` and `altbits` of Σ(−1)^{i+1}a_i;
    - `dmod7` = 6N·s(a, N) mod 7;
    - `phimod5` = the Rademacher Φ mod 5;
    - `ms_mod3`: 2[a/N]⁺ of 11a, mod 3;
    - `rem1mod4` = #{Euclid remainders ≡ 1 mod 4};
    - `midbits` = bitlength of the junction quotient a_{m+1}, where s_m ≤ √N < s_{m+1};
    - `prealt_sign` = the sign of the alternating sum over the first half;
    - `dr1_mod7` and `dr1_bits` of N²·s₁(a, 1; N), the shifted d = 2 Fourier–Dedekind sum (Dedekind–Rademacher
      type).
  - If gcd(a, N) > 1 the statistics are 0 and the nonunit is recorded (rule 3).
  - Every Euclid path is also scanned for a remainder sharing a factor with N (`path_hits`, rule 3; Prop VIS).
- **FD2 (6), role `menu`.** N²·s₁(a, b; N) for (a, b) ∈ {(isqrt, euler3), (cbrt2, fermres), (isqrt, invsqrt)}:
  the value mod 7 and its bitlength.
- **FW (7).** Disc −4N, prime form f of the least odd split prime ℓ, reduced forms f^j for j = 1..32.
  - Role `menu`: `fw_a1mod4`, `fw_a0mod3`, `fw_bitsum`, `fw_bneg`, `fw_cflen` (CF length of |b|/2a).
  - Role `menuG` (genus-expected): `fw_minbits`, `fw_sq` (a square leading coefficient).
  - Ambiguous forms on the walk are tested for a split (rule 3).
- **CF3 (5).** The CF of √(3N), 256 steps.
  - Role `menu`: `cf3_sumbits`, `cf3_Q1mod5`.
  - Role `menuG`: `cf3_n_a1`, `cf3_maxbits`, `cf3_Qsq`. These are small-Q statistics: by leakmap/eps §5 the small
    Q_i carry quadratic genus symbols.
- **PC (3).** Positive controls.
  - `pc2_T1024`, `pc3_T1024` (src/eps/cf.py; firing ⇒ ε = +1).
  - `pcA_fermat_T1024`: Fermat's method from ⌈√N⌉ succeeds within 1024 steps. Firing ⇒ q − p ≤ 2^{7.5}·N^{1/4}·O(1)
    ⇒ p/q > 7/8 at these sizes ⇒ tb = 1.

**Known degenerate (label-free pilot); these can never fire.**
- Constant at n = 48 and/or 128: `{isqrt, isqrt2, invsqrt}_{maxbits, altsign, altbits}`,
  `{isqrt2, cbrt2, div3}_prealt_sign`, `{cbrt2, div3}_altsign`, `cbrt2_maxbits`, `cbrt2_altbits` (n = 128),
  `invsqrt_midbits`, `cf3_Qsq` (n = 128).
- These are size facts: the huge first quotient ≈ √N, N^{1/3} or 3.
- Within-stratum constant under the tb floor at n = 128: `isqrt_dr1_bits`, `invsqrt_dr1_bits`.

## 3. Predictions (registered)
- **P1.** Every `menu` statistic, for both targets and every n ∈ {48, 64, 128}, is null: |z| < 6.
  - (i) Theory (notes/classE_results.md §1–2, Prop VIS): the Euclid path of a/N meets a multiple of p only when
    a ≡ ±r/s (mod q) with rs < q/p, which for balanced N means a ≡ ±1 (mod q). That is a group-order (M) event
    for `fermres`/`euler3` and size-impossible for the others.
  - (ii) The archimedean arguments reduce to CF prefixes of √(4N), √(8N) and N^{1/3}, plus a Euclid path of the
    pair (⌊√N⌋, N mod ⌊√N⌋) with N mod ⌊√N⌋ ≡ −uv. No mechanism links these to ε. For tb the size is absorbed by
    the floor.
  - (iii) The shifted d = 2 sums are class E (an aggregate congruence plus a CF-type integer).
  - (iv) The FW/CF3 `menu` statistics avoid smallness and squareness of the coefficients.
- **P2.** `menuG` statistics may hit ε at n = 48 only (family Q, via small/square coefficients). At n ≥ 64 they are
  null. A hit there is reported as Q, not as a new family.
- **P3.**
  - The PCs have 0 violations at every size.
  - At n = 48 and n = 64, pc2/pc3 fire with z ≥ 8 for ε and `pcA_fermat` fires with z ≥ 8 for tb.
  - At n = 128 none fire.
- **P4.**
  - `path_hits`: STRUCTURAL hits only on `fermres`/`euler3` at a Fermat/Euler-liar event at one prime; essentially
    0 expected at n ≥ 48.
  - ACCIDENTAL hits at the rate Σ_paths Σ_i (1/p + 1/q) ≈ 8·40·2/√N, i.e. ~0 at n ≥ 48.

## 4. Seeds and runs
- `SEED_BASE` in `src/classE/harness.py`: pilot = 25.51·10¹², main = 25.52·10¹², fresh = 25.53·10¹² (+ n·10⁹).
  These are disjoint from all earlier project ranges (eps 9.1–9.6·10¹², newfamily 24.5–24.7·10¹²).
- Command: `PYTHONPATH=. .venv/bin/python -m src.classE.harness --n 48 64 128 --range main --count 200000 --workers 4`.
- Memory ≤ 2 GB RSS (16 GB shared machine); one heavy job at a time; `memory_pressure` is checked first.

## 5. Test and decision rules (as newfamily)
- **Test.** The eps stratified permutation z (`src/eps/analysis.py::strat_z`), with the strata of §1. A hit is
  |z| ≥ 6.0.
- **Replication.** A fresh-seed run of the hit size, with the same sign and |z| ≥ 4.
- **Validity, per size × target.**
  - V1: the planted control (ρ = 0.05) lies within ±3 of its expected z.
  - V2: the PCs have 0 violations, and for n ≤ 64 at least one PC has z ≥ 8.
  - V3: the within-stratum shuffle mean z² over `menu` statistics lies in [0.6, 1.5].
  - At n = 128 only V1 and V3 apply.
- **Conclusion per size × target.**
  - "Refuted for this menu at |ρ| ≥ ρ_MDE" (80% power at |z| = 6), if the cell is valid and has no confirmed
    `menu` hit.
  - Otherwise, "claim not refuted" for each confirmed hit, followed by the eps audit:
    - recompute from N in a fresh process;
    - finer floors: N mod 2¹² for ε, top 16 bits for tb;
    - for tb, a hit must also survive the joint strata (top 8 bits, N mod 16).
- `menuG` hits are reported separately and classified as family Q (P2).

## 6. Toy test of the shifted d = 3 sums (exact O(N) evaluator; `src/classE/toy3.py`)
- **Populations.** All N = pq with p < q primes of n/2 bits with the top two bits set, N ≡ 1 (mod 4), for
  n ∈ {20, 22}. The populations are exhaustive.
- **Statistics (16).** S3 = N³·s₁(a, b, 1; N) and Z3 = ⌊N·s₁(a, b, 1; N)⌋, for
  (a, b) ∈ {(isqrt, cbrt2), (fermres, euler3), (isqrt, invsqrt), (hash2, isqrt)}. For each: S3 mod 7,
  bitlen|S3|, Z3 mod 7, sign Z3.
- **Targets.** ε and tb. Floors: N mod 2⁸ for ε; the top 6 bits for tb.
- **Prediction.** Null: |z| < 4.5 for all 64 tests.
- **Status: a sanity check only.**
  - There are only a few hundred to ≈ 1 100 instances, with heavy prime reuse, so the z are anti-conservative.
  - The MDE is |ρ| ≈ 0.2–0.35.
  - A |z| ≥ 4.5 would trigger a re-test at n = 24 (exhaustive), not a claim.
