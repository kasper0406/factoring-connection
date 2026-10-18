# Pre-registration: correlation battery for ε = (−1/p) and QR_2 on D_n (agent `epsbattery`, 2026-10-18)

NEXT_STEPS.md item 1. This note and the code it names (`src/eps/*`, `tests/test_eps.py`) are committed together, **before**
any battery instance (seed ranges `main`, `fresh`, `supp` below) has been generated. Before freezing, the runs were: the unit
tests (small primes and the `test` seed range); a cost benchmark on the separate `bench` range (statistics only; labels were
neither computed nor inspected); and one **pilot** of the full pipeline on 3 000 instances of `D_32` from the `test` range
(labels inspected; disclosed in §2a; it led to no menu change other than adding the §3 items PC4/PC5, which came from the
orchestrator). Any later change to `src/eps/{cf,algebra,forms,operational,analysis,harness}.py` is
listed as a deviation in `notes/eps_results.md`.

## 1. Question, distribution, targets

* `D_n`: `N = pq` from `balanced_semiprime(n, seed)` (`src/common/instances.py`: `p`, `q` independent primes of `n/2` bits
  with the top two bits set, `p ≠ q`), conditioned on `N ≡ 1 (mod 4)`. Then `p ≡ q (mod 4)`.
* Target 1: `ε = (−1/p) = (−1/q) ∈ {±1}`, on all of `D_n`; coded `y = [ε = +1]`.
* Target 2: `QR_2 = [2 is a square mod N] ∈ {0,1}`, on `D_n ∩ {N ≡ 1 (mod 8)}` (there `(2/p) = (2/q)`).
* Question: does a statistic of the menu (§2), computed from `N` alone, predict the target beyond the information floor
  (the best predictor from `N mod 2^k`)?
* Context: deciding ε for balanced semiprimes is Seurin's conjectured-hard "2-Φ/4-Hiding" problem (PKC 2014, eprint
  2013/256; `notes/codex_r6_eps_reduction.md`), whose prime windows and mixture weights differ from `D_n`. A positive
  result here would contradict that conjecture for a nearby distribution; a null is the expected outcome.
* Sizes `n ∈ {32, 48, 64}`; `10^5` instances of `D_n` per size; the QR_2 sample is the `N ≡ 1 (mod 8)` part of the same
  instances (expected ≈ `5·10^4`).

## 2. The statistic menu (94 statistics, frozen; order = `MENU` in `src/eps/operational.py`)

Conventions. All statistics are exact bounded integers computed by `compute_stats(N)` (OPERATIONAL: argument `N` only; it
never imports the instance generator or harness, never calls a factoring or primality routine; enforced by
`tests/test_eps.py::test_operational_modules_do_not_see_factors`). `jac(v)` = Jacobi symbol `(v mod N | N) ∈ {−1,0,1}`; before
every Jacobi symbol `gcd(v, N)` is computed and a value `> 1` is recorded as a nonunit (rule 3) and gives `0`.
`par(v)` = parity of the representative of `v` in `[0, N)`. `F_N = (1+X)^N − 1 − X^N`. `S_{r,M}(c) = Σ_j [X^{r+jM}]F_N · c^j`
= coefficient `r` of `F_N mod (N, X^M − c)` (Track A sections). A sketch `F_N mod (f, N)` has coefficients `a_0, a_1, …`
(Track B). Sums `jacsum`/`parsum` run over all `M` (or `deg f`) values. Roles: **NC** negative control (floor-type),
**menu**, **PC** positive control (§3).

**NC: floor-type negative controls (12).** `nc_N_mod8` = `N mod 8`; `nc_jac2_N` = `(2|N)`; `nc_N_mod256` = `N mod 256`;
`nc_bit8 … nc_bit11` = bits 8–11 of `N`; `nc_N_mod3_is1` = `[N ≡ 1 mod 3]`; `nc_jacN_5`, `nc_jacN_7` = `(N|5)`, `(N|7)`;
`nc_lead16` = top 16 bits of `N` (`N >> (n−16)`); `nc_chisum1000` = `Σ_{odd primes ℓ ≤ 1000} (−N|ℓ)` (a proxy for
`L(1, χ_{−4N})`). Functions of `N mod 2^8` (the first three) are exactly degenerate under the test of §5; bits 8–11 test for
residual 2-adic floor above `2^8`; the others are odd-modulus or size floor functions (0 bits about ε by
`notes/trackE_analysis.md` Lemma F1 plus Dirichlet equidistribution; **conjectured** null at these sizes).

**CF: continued-fraction prefix of √N (16 menu + 10 PC).** (Plus 2 MR/SS positive controls, §3.) `P_0 = 0, Q_0 = 1, a_0 = ⌊√N⌋`, `P_{i+1} = a_i Q_i − P_i`,
`Q_{i+1} = (N − P_{i+1}^2)/Q_i`, `a_{i+1} = ⌊(a_0 + P_{i+1})/Q_{i+1}⌋`. Identity (tested): `A_{i−1}^2 − N B_{i−1}^2 = (−1)^i Q_i`.
For `T ∈ {64, 1024}`, over `1 ≤ i ≤ T`:
`cf{T}_n_a1` = #{`a_i = 1`}; `cf{T}_n_a2` = #{`a_i = 2`}; `cf{T}_sumbits` = Σ bitlen(`a_i`); `cf{T}_maxbits` = max bitlen(`a_i`);
`cf{T}_n_Qeven` = #{`Q_i` even}; `cf{T}_n_sq_even` = #{even `i`: `Q_i` a square}; `cf{T}_n_Q1mod4_odd` = #{odd `i`: `Q_i ≡ 1 mod 4`};
`cf{T}_n_2sq_odd` = #{odd `i`: `Q_i = 2k^2`}.
Positive controls (indicators, §3): `pc1_T{T}`, `pc2_T{T}`, `pc3_T{T}` (for ε), `pcq2_T{T}`, `pcq3_T{T}` (for QR_2).
Factor base for PC3/PCQ3: `{−1} ∪ {primes ≤ 256}` (55 coordinates).

**A: Track A section values at fixed public (M, c) (26).**
`(1,2)`, `(1,3)`: `A1c{c}_jac = jac(F_N(c))`, `A1c{c}_par = par(F_N(c))` (4).
`(2,3)`: with `s_r = S_{r,2}(3)`: `A2c3_jac0 = jac(s_0)`, `A2c3_jac1 = jac(s_1)`, `A2c3_jacsum01 = jac(s_0+s_1)`,
`A2c3_jacdif01 = jac(s_0−s_1)`, `A2c3_par0`, `A2c3_par1` (6).
`(4,2)`: `A4c2_jac{r}`, `A4c2_par{r}`, `r = 0..3` (8). `(4,5)`: `A4c5_jacsum`, `A4c5_parsum` (2).
`(8,3)`: `A8c3_jacsum`, `A8c3_parsum`, `A8c3_jac0` (3). `(16,5)`: `A16c5_jacsum`, `A16c5_parsum`,
`A16c5_mod3is0 = [S_{0,16}(5) mod N ≡ 0 mod 3]` (3).

**B: Track B sketches, ranks, residues (18).** `f = X^3 − 2`, sketch `(a_0,a_1,a_2)`: `Bx3m2_jac{0,1,2}`, `Bx3m2_par0`,
`Bx3m2_mod3is0 = [a_0 ≡ 0 mod 3]`, `Bx3m2_mod5is0 = [a_0 ≡ 0 mod 5]` (residues of the representative), `Bx3m2_jacres` =
`jac(Res(f, S))` (determinant of multiplication by `S`, gcd-checked elimination; 0 on a nonunit pivot), `Bx3m2_rank` = rank over
`Z/N` of that 3×3 matrix (gcd-checked; −1 on a nonunit pivot). `f = Φ_5`: `Bphi5_jac0`, `Bphi5_jacsum`, `Bphi5_parsum`.
`f = Φ_12 = X^4 − X^2 + 1`: `Bphi12_jac0`, `Bphi12_jacsum`. `f = X^2 − X − 1`: `Bfib_jac0`, `Bfib_jac1`, `Bfib_par0`.
`G_1 = (1+X)^{N−1} − X^{N−1} mod (X^3 − 2, N)`: `BG1x3m2_jac0`, `BG1x3m2_jac1`.

**BM: Berlekamp–Massey lengths over Z/N, prefix of 64 terms (4).** gcd-checked BM (−1 on a nonunit discrepancy).
`bm_Fc`: sequence `F_N(c) mod N`, `c = 2..65`. `bm_coef64`: the 64 coefficients of `F_N mod (N, X^64 − 3)`.
`bm_euler`: `c^{(N−1)/2} mod N`, `c = 2..65`. `bm_coef64_jacC` = `jac` of the top coefficient `C_L` of the final connection
polynomial of `bm_coef64`.

**CL: reduced primitive forms of discriminant −4N with a ≤ 50 (6).** Forms `(a, b, c)`, `b^2 − 4ac = −4N`, `|b| ≤ a ≤ c`,
`b ≥ 0` if `|b| = a` or `a = c`, `gcd(a,b,c) = 1`. `cl_cnt` (all `a ≤ 50`), `cl_cnt_odd` (`a` odd), `cl_cnt_prime` (`a` prime),
`cl_cnt_le10`, `cl_cnt_le25`, `cl_min_a` = least `a ∈ [3, 50]` with a form (51 if none). (Checked against the analytic class
number for small `N`.) These are functions of `N mod lcm(4a)`, i.e. floor-type for odd moduli: **conjectured** null.

**Menu vetting before the freeze (bench range only).** The first draft had (M,c) = (2, −1), (4, −1) and `G_1 mod X^2+1`.
The benchmark showed `jac(s_0 − s_1)` and `jac(G_1[1])` constant. Reason (**proved**): in `Z/N[i]`, `(1+i)^2 = 2i`, so
`(1+i)^N = σ 2^{(N−1)/2}(1+i)` with `σ = (−1)^{(N−1)/4}`; hence `s_0 = s_1` and `G_1` is real. These items were replaced by
(2,3), (4,5) and `G_1 mod X^3−2` before the freeze. Still constant on the bench sample and **kept** (they are what the menu
asks for; they document that ranks and BM lengths are generic): `Bx3m2_rank` (= 3), `bm_Fc`, `bm_coef64`, `bm_euler` (= 32
at n = 48, 64). A rank or BM length below generic needs a coincidence mod p, at rate Θ(1/p) (Track B/D notes). Constant
statistics are reported as DEGENERATE and cannot fire. Also seen on the bench: one coefficient of the `Φ_5` sketch is
exactly `0 mod N`, at a position fixed by `N mod 5` (198/198 bench instances at n = 32; **experimentally supported**). This
makes `Bphi5_jac0` partly a floor function of `N mod 5`; it is kept unchanged.

## 2a. Pilot disclosure and a registered prediction for the CF group

Pilot (`test` range, n = 32, 3 000 instances, before the freeze): all PCs had 0 violations and matched `z_pred`; the planted
control matched; the within-stratum shuffle gave mean `z^2 = 0.90`. Four CF statistics at `T = 1024` had `z ≈ 3.4–5.3`
(`cf1024_n_Q1mod4_odd` 5.3, `cf1024_maxbits` 4.2, `cf1024_n_a1` 3.5, `cf1024_n_Qeven` 3.4); scaled to `10^5` instances these
would be hits at n = 32. **Registered prediction (conjectured mechanism):** at n = 32 the prefix `T = 1024` covers a
non-negligible part of the period of √N (the period closes within 1024 steps for ≈ 2.5% of bench instances), and
period-level structure carries ε (an odd period forces ε = +1; at the end of each period `a_ℓ = 2a_0`, `Q_ℓ = 1`, and with an
odd period the odd- and even-index `Q`'s are interchanged). So these hits are predicted to be the PC1 mechanism, and to
shrink with `n` (no hit at n = 64, where neither the period nor any PC is reachable). The diagnostic is registered in
§6 step 6: the z of a CF-group hit recomputed on the instances whose period does **not** close within 1024 steps
(`first_Q1 = 0`). A CF hit that survives there, or at n = 64, is **not** explained by this mechanism. Proper-factor gcd hits
in the pilot exceeded the null `Σ(1/p + 1/q)` (78 vs 27, mostly `A16c5` and the BM sequences); they are reported by source
with the null comparison (rule 3), classified ACCIDENTAL where they match the null and "excess, mechanism not analysed"
otherwise.

## 3. Positive controls (one-sided known signals; not leakage)

Let `c_p(x) ∈ F_2` be 1 iff `(x|p) = −1`. All proofs use `A_{i−1}^2 ≡ (−1)^i Q_i (mod N)` (identity of §2).

* **PC1 (frozen by NEXT_STEPS):** the period `ℓ` (least `i ≥ 1` with `Q_i = 1`) is odd and `ℓ ≤ T`. *Proof (proved):*
  `A_{ℓ−1}^2 − N B_{ℓ−1}^2 = −1`, so `−1` is a square mod `p` and mod `q`: `ε = +1`.
* **PC2:** some odd `i ≤ T` has `Q_i = k^2` with `gcd(k, N) = 1` (checked; automatic for balanced `N`, since
  `k ≤ (2√N)^{1/2} < p`). *Proof:* `(A_{i−1}k^{−1})^2 ≡ −1 (mod p)`, so `ε = +1`. (PC1 is the case `k = 1`.)
* **PC3:** `e_{−1}` lies in the F_2-span of the exponent-parity vectors of the factor-base-smooth `(−1)^i Q_i`, `i ≤ T`.
  *Proof:* every factor-base prime is coprime to `N` (checked: `gcd(N, ∏_{ℓ≤256} ℓ) = 1`, else PC3 is disabled), so
  `((−1)^i Q_i | p) = 1` gives `⟨v, c_p⟩ = 0` for each relation vector `v`; so `c_p ⊥ span`; `e_{−1} ∈ span` gives `c_p(−1) = 0`.
  This is CFRAC's linear algebra: where it fires routinely, the same relations factor `N` in `L_N[1/2]` time. It is a
  calibration of how much of `ε` the CF prefix exposes, not a candidate.
* **PCQ2 (for QR_2):** some even `i ≤ T` has `Q_i = 2k^2`, `gcd(k, N) = 1`. Then `2` is a square mod `p` and mod `q`.
* **PCQ3 (for QR_2):** `e_2 ∈` span. Then `c_p(2) = c_q(2) = 0`.
* **PC4 (MR; GPT-6.1-sol r6, `notes/codex_r6_primality_lit.md` §9, added before the freeze):** for some base
  `b ∈ [2, 65]`, with `N − 1 = 2^s d`, `d` odd, the Miller–Rabin chain `x_0 = b^d`, `x_i = x_{i−1}^2` has `x_i ≡ −1 (mod N)`
  for some `1 ≤ i ≤ s − 1`. Then `x_{i−1}^2 ≡ −1`, so `ε = +1` (**proved**). **PC5 (SS):** `b^{(N−1)/2} ≡ −1 (mod N)` for some
  `b ∈ [2, 65]`; this is `x_{s−1} = −1` (`s ≥ 2` on `D_n`), so PC5 ⊆ PC4. The r6 note gives the per-base rate as the
  Fermat-pass probability `F` times a constant in `[1/4, 1/3]` (PC4) or `1/4` (PC5, only if `v_2(p−1) ≠ v_2(q−1)`); `F` is
  exponentially small for random balanced `N`, so both are **predicted never to fire** at n ≥ 32 (0/575 on the bench). They
  are kept as zero-power controls: they test the never-wrong property only. Nontrivial square roots of 1 met in the chain are
  counted (`mr_split`) and never used.
* Known non-one-sided signal, kept as a menu statistic: `cf{T}_n_2sq_odd > 0` gives `(−2|p) = (−2|q) = 1`, which forces
  `N ≡ 1 (mod 8)` and `ε = (2|p)`, i.e. it exposes the **product** `ε·(2|p)`, not `ε` or `QR_2` alone. **Predicted:** null for
  each target separately (by the symmetry `ε ↔ −ε` of the remaining unknowns; **conjectured**).
* Bench firing rates (bench seeds, `≈ 190` instances per size, statistics only) at `T = 1024`: n = 32: PC1 1.5%, PC2 43%,
  PC3 45%, PCQ2 23%, PCQ3 24%; n = 48: PC1 0%, PC2 15%, PC3 44%; n = 64: PC1, PC2, PC3 0/199, PCQ2 = PCQ3 = 0.5%.
  Heuristic: PC2 fires at rate ≈ `#odd steps / N^{1/4}` (a random `Q ≈ √N` is a square w.p. `≈ N^{−1/4}`), PC1 much more rarely.

**Validity rule (per size and target).**
* V1: every positive control has **zero** violations (fires while the target is `ε = −1` / `QR_2 = 0`). One violation
  invalidates the run (a pipeline or theory bug).
* V2: for a one-sided indicator with firing frequency `f` and positive-class rate `π`, zero violations give a
  point-biserial correlation `ρ = √(f(1−π)/((1−f)π))`; the predicted z is `z_pred = ρ√(n−1)`. A PC is **powered** if
  `z_pred ≥ 8`. Every powered PC must have stratified `z_obs ≥ 6.0` (detected) and `|z_obs − z_pred| ≤ 3`.
* V3: the battery for a (size, target) is **validated** only if at least one of its own PCs is powered and passes V2.
  PC1 is reported separately, as NEXT_STEPS freezes it: if PC1 is unpowered at a size (expected at n = 48, 64), the note
  states that PC1 cannot validate that size and that validity there rests on PC2/PC3 (or PCQ2/PCQ3), as registered here.
  If no PC is powered, that (size, target) is reported **unvalidated** and its null results are not claimed.
* Planted pipeline control (analysis side; uses the label by construction): `x_pl = y` with probability 0.525 and `1−y`
  otherwise, from `random.Random(f"planted|{seed}")`. Expected `ρ = 0.05`; it must give `|z − 0.05·√neff| ≤ 3`. It checks the
  analysis code only.

## 4. Instances and seeds (disjoint from every earlier run)

Earlier runs used seeds `< 10^6` (configs in `experiments/configs/`). Ranges (`SEED_BASE` in `src/eps/harness.py`):
`bench = 9.1·10^12` (cost benchmark; never reused), `main = 9.2·10^12` (the battery), `fresh = 9.3·10^12` (re-runs of hits),
`supp = 9.4·10^12 + n·10^9` (supplementary curve), `test = 9.5·10^12` (unit tests). For size `n`, seeds `base + j`,
`j = 0, 1, …`; the sample is the first `10^5` seeds (in order) whose `N ≡ 1 (mod 4)`. Seeds are hashed with the family and
`n` (SHA-512 string seeding), so the same `j` gives unrelated instances at different `n`.

## 5. Floor and test statistic

**Floor baseline (reported, not used for the decision).** For `k = 2..12`, the majority-label-per-class predictor from
`N mod 2^k`, 2-fold cross-validated (even/odd positions), and a χ² test of the target against `N mod 2^8`. Expectation
(**conjectured** from Lemma F1 plus Dirichlet): accuracy `1/2 ± 0.002` for ε.

**Test statistic (decision).** Stratified permutation correlation (Mantel / CMH-correlation type). Strata `h = N mod 2^8`
(64 strata of ≈1560 on `D_n`; 32 strata of ≈1560 for QR_2). For statistic `x`, target `y`:
`U = Σ_h Σ_{i∈h} (x_i − x̄_h)(y_i − ȳ_h)`, `V = Σ_h S^h_{xx} S^h_{yy} / (n_h − 1)`, `z = U/√V`.
`V` is the exact variance of `U` when `y` is permuted within strata. Under H0 (`y ⟂ x | N mod 2^8`) `z → N(0,1)`.
**Why the floor cannot fool it:** any function of `N mod 2^k`, `k ≤ 8`, is constant inside every stratum, so it contributes
`U = V = 0` exactly (reported DEGENERATE, `z := 0`). For `x = g(N mod 2^8) + h(N)`, only the within-stratum variation of `h`
enters. A floor effect beyond `2^8` would show up in the negative controls `nc_bit8..11` (and is re-checked with `2^12`
strata for any hit). A linear residualisation on `N mod 2^k` would not have this exactness, which is why it was not chosen.
The unstratified point-biserial z is also reported, for comparison only. All bounded statistics; the Lindeberg ratio
`max (x_i − x̄_h)^2 / Σ (x_i − x̄_h)^2` is reported. Calibration: one within-stratum shuffle of `y` (numpy seed 0) must give
mean `z^2 ≈ 1` over the non-degenerate statistics.

**Threshold.** A **hit** is a menu or NC statistic with `|z| ≥ 6.0` at some size. Positive controls are not hits. With
`S = 100` nominal statistics per (size, target), the family-wise error is `≤ 100 · 1.97·10^{−9} ≈ 2·10^{−7}` (below the
two-sided 5σ level `5.7·10^{−7}`); over all `94 × 2 × 3 = 564` tests it is `≤ 1.1·10^{−6}`.

**Power analysis.** For a true within-stratum correlation `ρ`, `z ≈ ρ√neff`, `neff = Σ_h (n_h − 1)`. The minimum detectable
effect at 80% power is `ρ_MDE = (6.0 + 0.84)/√neff`: `≈ 0.0216` for ε (`neff ≈ 10^5`), `≈ 0.0306` for QR_2 (`neff ≈ 5·10^4`).
For a balanced binary statistic this is a prediction accuracy of `1/2 + ρ/2`: `51.1%` (ε), `51.5%` (QR_2). For a binary
statistic firing with frequency `f`, it is a shift `Δ = (ρ/2)√((1−f)/f)` of `P(y = 1 | x = 1)`. These are reported next
to every null. A degenerate statistic detects nothing.

## 6. Hits, leakage audit, replication, refutation criterion

For every hit:
1. **Static audit:** confirm the statistic's code path lies in the operational modules (AST test) and uses only `N`.
2. **Recompute from N alone:** in a fresh process importing only `src.eps.operational`, recompute the statistic for all
   instances from the stored `N`; values must match exactly.
3. **Floor sensitivity:** recompute z with strata `N mod 2^12`, and with strata `(N mod 2^8, N mod 105)`.
4. **Generator artefacts:** z of the statistic against `⌊16q/p⌋`, `nc_lead16`, seed parity.
5. **Fresh seeds:** re-run the full menu at that size on `10^5` instances from `fresh`. **Confirmed** iff same sign and
   `|z| ≥ 4.0` (one pre-specified test).
6. **Mechanism:** try to explain a confirmed hit by a one-sided identity of the §3 kind. For a CF-group hit, report its z
   on the instances whose period does not close within 1024 steps (§2a), and its z at the other sizes.
   A hit explained by a registered PC mechanism is reported as such and does not count against the refutation criterion
   below only if this diagnostic removes it (|z| < 4 on the non-closing instances); otherwise it counts as a menu hit.

**Refutation criterion (from NEXT_STEPS).** If no menu statistic other than the positive controls has a confirmed hit at
any of the three sizes, the claim "some statistic of this menu predicts ε (or QR_2) on D_n with `|ρ| ≥ ρ_MDE`" is refuted
for this menu and these sizes (at the sizes that are validated, V3). This excludes only effects of that size. It says
nothing about smaller correlations or about other statistics, and nothing about the reduction question (NEXT_STEPS 3).

## 7. Supplementary (descriptive): positive-control frequency versus n and T

`n ∈ {24, 32, 40, 48, 56, 64, 80, 96}`, `2·10^4` instances each from `supp`, CF only with `T_max = 16384`. Reported:
the firing frequency of PC1, PC2, PC3 (and PCQ2, PCQ3 on `N ≡ 1 mod 8`) at `T ∈ {64, 256, 1024, 4096, 16384}`, with
violation counts (must be 0). No hypothesis test. This quantifies how much of ε negative Pell and its prefix relatives
expose cheaply.

## 8. Cost (bench range, 16 concurrent workers) and budget

Final menu, mean per instance: 9.5 ms (n = 32), 10.1 ms (n = 48), 15.9 ms (n = 64); p95 ≤ 38 ms (an earlier bench of the
draft menu gave 16–27 ms under heavier contention). Main battery ≈ 3·10^5 × 12 ms ≈ 1 CPU-h; supplementary ≈ 1.3 CPU-h
(30 ms per instance at n = 96, T = 16384). No statistic was dropped for cost (largest group: BM, 6.5 ms at n = 64).

## 9. Commands (from the project root)

```
.venv/bin/python -m pytest -q tests/test_eps.py
.venv/bin/python -m src.eps.bench --seeds 400 --workers 16                    # done before the freeze
.venv/bin/python -m src.eps.harness --n 32 48 64 --range main --workers 30   # -> experiments/raw/eps_main_n{n}.jsonl.gz
.venv/bin/python -m src.eps.harness --n 24 32 40 48 56 64 80 96 --range supp --curve --count 20000 --workers 30
.venv/bin/python -m src.eps.analysis experiments/raw/eps_main_n32.jsonl.gz experiments/raw/eps_main_n48.jsonl.gz \
    experiments/raw/eps_main_n64.jsonl.gz --out experiments/analysis/eps_battery_summary.json
```
