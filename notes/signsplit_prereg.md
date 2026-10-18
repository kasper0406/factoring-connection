# signsplit: pre-registration (frozen before the experiments)

Agent `signsplit`, 2026-10-18. This file is committed before any run of
`src/signsplit/run.py main` on the instance set below.

**What ran before this file was written:**
- The exhaustive identity verification (`src/signsplit/verify_identity.py`, work-plan step 1, 0
  failures). This is not an experiment on the instance set.
- Smoke and debug runs at 24, 28 and 36 bits (seeds 0–2 and ≥ 100), outside the frozen set.
- Count-only model predictions on synthetic N = 2^n + 2^{n/2} + 1 (made coprime to 2310). These
  never touch a factor; they are the frozen predictions of §5.

## 1. Hypothesis under test (the user's claim, treated as unvalidated)

H: In Harvey–Hittmeir (HH 2022, Alg. 4.3) search with M = m_1⋯m_t·(2), one weight pair serves the
2^t residues with the same square class σ². The offsets satisfy τ(σ(ε)) ≡ Δ_0 + Σ ε_i Δ_i
(mod M²). A meet-in-the-middle over ε then tests the 2^t branches in about 2^{t/2} work, giving a
2^{t/2}-type saving over HH at matched success.

**Our prediction (from the accounting in notes/signsplit_results.md §3, made before the
experiment): H is refuted for this implementation.** No sign-split variant beats the HH baseline,
for these reasons:
- The literal per-group MITM (V2) destroys the shared baby list. It becomes per-target BSGS, with
  exponent 1/4.
- The only shared-list way to let one giant serve both signs is V3, the anti-lock. It costs lattice
  index m² per collapsed prime, and that loses to the HH lock by m^{2/5}/2^{3/5} per prime
  (1.024 at m = 3, 1.256 at m = 5).

## 2. Instances and sizes (frozen)

- `src.common.instances.balanced_semiprime(bits, seed)`, with bits ∈ {32, 40, 48, 56, 64} and
  seeds 0..19, for 100 instances in total.
- Public configuration: ratio bound q/p ≤ c = 2. Candidates for p range over [√(N/2), √N].
- Modulus menu: M ∈ {6, 30, 210, 2310}, i.e. t = 1, 2, 3, 4 odd primes, with 2 always included as
  in HH.

## 3. Variants (all operational: they see N, public parameters and a seed only)

| id | description | parameters |
|---|---|---|
| **V0** | HH Alg. 4.3 baseline. One Prop. 3.3 reduction and one giant per (σ, σ_0), with σ over all φ(M) units. Shared baby list β^{M² i}, 0 ≤ i ≤ 2λ, λ = ⌈4N^{1/2}I^{1/2}m_0^{−3/2}/M²⌉, I = M. Sort-and-match, then Alg. 4.1. Randomised base β (instead of HH's deterministic Lemma 2.5, as in cycloharvey). | (M, m_0) |
| **V1** | Shared-lattice HH: one reduction per (σ², σ_0). The 2^t giants are built by a Gray-code walk on the signed subset sum, with carries. Lists are identical to V0. | same as V0 |
| **V2** | Per-group sign MITM (literal claim). Left list: g·2^{t−t_R} elements per group. Right list: 2^{t_R}·⌈(2λ+1+c_hi−c_lo)/g⌉ elements per group, group-specific because Δ_R depends on (a,b). Product tree of the left lists, multipoint evaluation at the right lists, then gcds. | (M, m_0, t_R, g) |
| **V3** | Sign collapse (anti-lock) on a set A of odd primes: m \| a, m \| b, (a/m)N + (b/m)σ² ≡ 0 (mod m). One giant per (σ_L, ±σ_A), a shared baby list, and index I = M·∏_A m². | (M, m_0, A) |

**Parameter rules.** All are deterministic functions of N, in `src/signsplit/accounting.py`.
- V0's m_0 minimises its own count-only work (giant + baby). The grid is {1, 2, 3, 4} ∪
  {round(2^{k/8})}, followed by integer refinement ±3.
- **E1 (matched; "identical parameters except branch handling").** V1, V3(A={3}), V3(A=all odd)
  and V2 all run at V0's (M, m_0). For V2, t_R = max(1, ⌊t/2⌋) and g = 1.
- **E2 (own optimum).** V3 uses its own m_0 for each A. V2 minimises the upper-bound formula
  G·(2^{t−t_R}g + 2^{t_R}⌈(2λ+1+2t)/g⌉) over (m_0, t_R ∈ 1..t, g ∈ powers of 2).
- **E3 (best over menu).** For each instance, each variant takes its minimum work over every
  (M, A, parameter) record above.
- V2 runs operationally only when its list total is ≤ 300 000. Otherwise its exact list counts come
  from `v2_counts`, which performs the same lattice reductions.

## 4. Counters (frozen)

- **Primary:** `work` = total list elements (giants + babies; for V2, left + right). This is the
  quantity the 2^{t/2} claim is about, and HH's cost (4.3) is linear in it, up to polylog
  weights.
- **Secondary (raw):**
  - ring multiplications `rn_mul` (square-and-multiply counted per exponentiation) and
    `rn_pow`;
  - `gcd`, `inv`;
  - lattice reductions `lattice_red` and `gauss_iter`;
  - `poly_mul` and `poly_mul_dimsum` (product tree), `poly_rem`;
  - `exact_match`.
- **Success:** a proper factor is returned.
- **Classification (harness, using p):**
  - STRUCTURAL: the hit is the true branch (σ, σ_0) at the true baby index.
  - ACCIDENTAL: everything else, including `gcd_baby_order`, `gcd_base`, or a collision on another
    branch.
- **Coverage check (harness):** the true branch's u = aq + bp lies in its giant window.

## 5. Frozen predictions (count-only, synthetic N; best over the menu)

| n | best V0 | best V3 | best V2 (E2) | V3/V0 | V2/V0 |
|---|---|---|---|---|---|
| 32 | 97 | 97 | 272 | 1.00 | 2.8 |
| 40 | 267 | 269 | 1080 | 1.01 | 4.0 |
| 48 | 779 | 789 | 4308 | 1.01 | 5.5 |
| 56 | 2297 | 2349 | 17220 | 1.02 | 7.5 |
| 64 | 6797 | 6871 | 68868 | 1.01 | 10.1 |

- **P1.** V1 has exactly the same giant and baby counts as V0, and `lattice_red`(V1) =
  `lattice_red`(V0)/2^t.
- **P2.** At its own optimum (E2), V3(A = all odd) is worse than V0 at the same M by a factor that
  grows with M and n. At 64 bits and M = 2310 the predicted factor is about 120.
- **P3.**
  - V2's slope in n (log2 work against n) is about 1/4, against about 1/5 for V0.
  - V2/V0 ≥ 2 for every instance with n ≥ 40.
- **P4.**
  - Coverage holds in 100% of records.
  - Success is 100%, except for small-order bases (retry with up to 3 base seeds, reported).
  - At n ≥ 48 the STRUCTURAL share is ≥ 90%.

## 6. Decision rule for H (frozen)

Let R(n) be the median over instances of (best sign-split work, min over V2 and V3) / (best V0
work), in E3.
- **SUPPORTED** for this implementation if R(56) ≤ 0.5 and R(64) ≤ 0.5. This is a 2^{t/2}-type
  gain with slack ×2 for t ≥ 3.
- **REFUTED** for this implementation if R(56) ≥ 0.9 and R(64) ≥ 0.9, i.e. no ≥ 10% gain.
- **INCONCLUSIVE** otherwise.

V1's construction-cost savings (`rn_mul`, `lattice_red`) are reported separately. They are not
counted as support for H: they do not change any list size or the collision-detection cost.

A discrepancy between the counts in §5 and the instance medians of more than 10% is reported as a
discrepancy.

## 7. Commands

```
.venv/bin/python -m src.signsplit.verify_identity 2001
.venv/bin/python -m pytest -q tests/test_signsplit.py
.venv/bin/python -m src.signsplit.run main --procs 24
.venv/bin/python -m src.signsplit.run summarize
```
