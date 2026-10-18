# leanform: Lean 4 formalization of the Theorem D′ lemmas

Agent `leanform`, 2026-10-18. Code: `lean/` (Lake project; see `lean/README.md`). Source of
the statements: `notes/dcltheorem_results.md` §2–4, which remains authoritative.

## 0. Verdict (plain language)

**Lean has checked the algebraic and elementary number-theoretic core of Theorem D′.**
Lemmas 1, 2, 3, 4, 5, 6, 7, 8 and 10, plus the coverage criterion, are fully proved in Lean 4
with Mathlib. So is a capstone theorem that chains Lemmas 3, 4, 5 and 7.
- **No `sorry`. No `axiom` declaration.** `#print axioms` on all 31 main results shows only
  Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`).
- **The large sieve is not proved.** It enters as an explicit hypothesis `ArithmeticLargeSieve`
  (a `Prop` passed as an argument). Lemma 9 and the exceptional-set bound `E_b` are proved *from*
  it.
- **Not formalized:**
  - the analytic inputs (PNT in the progression 3 mod 4, `|S_n| ≍ P²/n²`);
  - the final asymptotic bookkeeping of §4(iii) and the cost analysis §4(ii);
  - the algorithm itself as a program. The capstone works with the giant's residues, not with
    executable code.

Status label for every "proved" row below: **proved (machine-checked)**. The large-sieve row is
**hypothesis (assumed, not proved)**.

## 1. Build

```sh
cd lean
lake exe cache get        # prebuilt Mathlib oleans (no Mathlib compilation)
lake build                # → "Build completed successfully"
lake env lean AxiomCheck.lean
```

- Toolchain `leanprover/lean4:v4.35.0`; Mathlib tag `v4.35.0` (commit `56a0058…`).
- Confirmed on 2026-10-18: `lake build` completes with no errors. The only warnings are Mathlib's
  generic "consider adding `module`" notices.
- `.lake/` (≈7 GB, mostly Mathlib oleans) is git-ignored by `lean/.gitignore`.

## 2. Lemma → Lean name → status

Notation as in the note: `u = x^M`, `γ = x^(M ℓ²)` with `M = L/ℓ²` and base `x` (= 2), and
`d = ord_q(x)`. All names are in namespace `FactoringLean`.

| note item | Lean name (file) | status | what exactly is proved |
|---|---|---|---|
| (a) exponent identity | `exponent_identity` (Exponents) | proved | `B + A(pq) − (Bp + Aq) = (p−1)(Aq − B)` over ℤ |
| (a) Fermat consequence | `pow_exponent_identity_zmod`, `pow_exponent_identity_modEq` (Exponents) | proved | `a^(B+AN) ≡ a^(Bp+Aq) (mod p)` for **every** `a` (`p ∤ a` not needed since `x^p = x` in `ZMod p`) |
| Fermat, integer exponents | `units_zpow_eq_of_sub_dvd` (Exponents) | proved | for a unit mod `p`: `(p−1) ∣ e₁ − e₂ ⇒ u^e₁ = u^e₂` |
| (b) coverage criterion | `coverage_criterion`, `coverage_criterion'`, `root_iff_of_modEq` (Coverage) | proved | prime `ℓ ∤ N`: `ℓ ∣ Bp² − AN ⇔ ℓ ∣ Bp − Aq`; root classes `p ≡ r` |
| (c) Lemma 1, class lock | `lemma1_class_lock` (Coverage) | proved | `ℓ` prime, `ℓ ∤ r`, `ℓ ∣ Br² − AN`, `ℓ² ∣ rs − N`, `ℓ² ∣ c − (Br + As)`, `ℓ ∣ p − r` ⇒ `ℓ² ∣ Bp + Aq − c` |
| Lemma 2, lower end (AM–GM) | `lemma2_lower` (SingleTarget) | proved | `⌊√(4ABN)⌋ ≤ Bp + Aq` |
| Lemma 2, upper end (convexity) | `lemma2_upper` (SingleTarget) | proved | `0 < a ≤ x ≤ b ⇒ Bx + AN/x ≤ max(F(a), F(b))` over ℝ |
| Lemma 2, indexing | `lemma2_index` (SingleTarget) | proved | `0 ≤ i* < m`, `g − (g_lo + j*m) = i*`, `0 ≤ j* ≤ (g_hi − g_lo)/m` |
| Lemma 2, window | `lemma2_window` (SingleTarget) | proved | `K = ℓ² > 0`, `Lo ≤ Kg + c ≤ Hi` ⇒ `⌈(Lo−c)/K⌉ ≤ g ≤ ⌊(Hi−c)/K⌋` |
| (d) Lemma 3, mod p | `lemma3_core`, `lemma3` (Exponents) | proved | `u^(B+AN−c) · γ^(−(g_lo+jm)) = γ^(g−g_lo−jm)` in `(ZMod p)ˣ`; no order hypothesis |
| (d) Lemma 4, mod q | `lemma4_core`, `lemma4`, `lemma4_exponent_difference` (Exponents) | proved | giant `= γ^(i*) · u^((A−B)(p−1))` in `(ZMod q)ˣ`; difference `≡ (A−B)(p−1) (mod q−1)` |
| Lemma 5, one-sided extraction | `lemma5_one_sided`, `not_full_iff` (Bands) | proved | `p ∣ G − γ^(i*)` and no baby collision mod `q` ⇒ `gcd(∏_{i<m}(G − γ^i), pq) = p` |
| (e) Lemma 6, four bands | `lemma6_clean_band`, `lemma6_clean_band_prime` (Bands) | proved | 4 pairwise disjoint prime sets, all `> Y`, and `0 < n < (Y+1)⁴` ⇒ some set has no prime dividing `n`; for prime `q ≤ Y⁴`, take `n = q − 1` |
| (f) Lemma 7, core | `coprime_dvd_mul_iff`, `coprime_dvd_mul_iff_int` (Bands) | proved | `gcd(e,d) = 1 ⇒ (d ∣ ex ⇔ d ∣ x)` |
| (f) Lemma 7, FULL criterion | `lemma7_core`, `lemma7_full_iff` (Bands) | proved | any group, `gcd(M, ord x) = 1`: giant FULL ⇔ `∃ δ ∈ [−i*, m−1−i*]`, `ord x ∣ D − ℓ²δ` |
| (f) cleanliness ⇒ coprimality | `clean_coprime`, `coprime_of_clean`, `lemma7_coprime_zmod` (Bands) | proved | no band prime divides `q−1` ⇒ `gcd(∏ℓ², q−1) = 1` ⇒ `gcd(M, ord_q x) = 1` |
| (g) Lemma 8, requested form | `lemma8` (SmallOrder) | proved | `#{odd primes q : ord_q(2) ≤ T} ≤ T²` |
| (g) Lemma 8, note's form | `lemma8_sharp` (SmallOrder) | proved (floor-log variant) | primes `> 2^b` (`b ≥ 1`) with `ord_q(2) ≤ T`: `2b·#s ≤ T(T+1)`. With `b = ⌊log₂ Q₀⌋` this is the note's `T(T+1)/(2 log₂ Q₀)` with `log₂` floored |
| (h) Lemma 10, class count | `card_le_of_pairwise_modEq` (CRTCount) | proved | a pairwise-congruent-mod-`M` subset of `[a, a+len)` has `≤ len/M + 1` elements |
| (h) Lemma 10, one triple (CRT) | `lemma10_pairwise_modEq`, `lemma10_single` (CRTCount) | proved | `ℓ` prime, `ℓ ∤ B`, `ℓ ∤ d`, `d > 0`: solutions of `ℓ ∣ Bp − Aq`, `d ∣ (A−B)(p−1) − ℓ²δ` in `[a, a+len)` number `≤ len/(ℓ·d/gcd(d, A−B)) + 1` |
| (h) Lemma 10, modulus size | `lemma10Modulus_ge` (CRTCount) | proved | `ℓd ≤ (ℓ·d/g′)·|A−B|` for `A ≠ B` (the note's `ℓd/g′ ≥ ℓd/H`) |
| (h) Lemma 10, union bound | `lemma10_union` (CRTCount) | proved | `|F| ≤ Σ_{(v,ℓ,δ)} (len/M_t + 1)` |
| Lemma 10, final real bound | — | not formalized | substituting `|V| ≤ H²`, `|𝒫_k| ≤ 2^kY`, `m̄_k`, `d > T` into the sum |
| Lemma 9, distinct classes | `same_class_dvd`, `weights_eq_of_same_class` (Coverage); `excludedResidue_spec`, `card_weights_le` (LargeSieve) | proved | primitive weights in `[1,H]²`, `H² < ℓ`, `ℓ ∤ q`: equal classes ⇒ equal weights; `|V| ≤ H²` |
| **Arithmetic large sieve** | `ArithmeticLargeSieve` (LargeSieve) | **HYPOTHESIS (assumed, not proved)** | prime-term form: `|U|·(1 + Σ_{ℓ∈𝒫} ω/(ℓ−ω)) ≤ N₀ + Q²` (Montgomery 1968; Montgomery–Vaughan 1973) |
| Lemma 9 | `lemma9_of_largeSieve` (LargeSieve) | proved **from the hypothesis** | `|U(q)| · |V| · Σ_{ℓ∈𝒫} 1/ℓ ≤ N₀ + Q²` |
| `E_b` exceptional set | `exceptional_Eb_of_largeSieve` (LargeSieve) | proved **from the hypothesis** | `(Σ_{q∈𝒬} |U(q)|) · J ≤ |𝒬|·(N₀ + Q²)` whenever `J ≤ |V| Σ_{ℓ∈𝒫(q)} 1/ℓ` |
| `Σ 1/ℓ ≥ c/n`, `|S_n| ≥ c₀P²/n²` | — | not formalized | needs PNT (in progressions); outside this scope |
| Capstone (Lemmas 3+4+7+5) | `single_target_extraction` (SingleTarget) | proved | `p ≠ q` primes, the giant's residues mod `p` and mod `q` as computed, `ℓ²g = Bp + Aq − c`, `g − g_lo = j*m + i*` with `0 ≤ i* < m`, clean band, and the FULL condition failing for all `δ ∈ [−i*, m−1−i*]` ⇒ `gcd(f(G), N) = p` |
| §4(iii) final assembly; §4(ii) cost | — | not formalized | asymptotic sums and the bit-complexity table |

**List of `sorry`s: none.** (`grep -rn sorry lean/FactoringLean` is empty, and `sorryAx` is
absent from every `#print axioms`.)
**List of axioms: none beyond Lean's standard three.** The large sieve is a hypothesis, not an
`axiom`.

## 3. Remarks on the formal statements

- **(a) is slightly stronger than the note.** In `ZMod p`, `x^p = x` holds for all `x`, so
  `a^(B+AN) ≡ a^(Bp+Aq) (mod p)` needs no `p ∤ a`. Lemmas 3 and 4 use integer (possibly
  negative) exponents, so they are stated in the unit group `(ZMod p)ˣ`. That matches the
  algorithm, which uses `γ^(−1)` mod `N`.
- **The large-sieve hypothesis is a true, weakened form of the theorem.**
  - Only the `d = 1` and prime `d` terms of `J` are kept. Since every term of `J` is ≥ 0, this
    weakens Montgomery's inequality.
  - `U` lies in `(M₀, M₀+N₀]`. Each `Ω ℓ` is a set of residues in `[0, ℓ)` with `|Ω ℓ| < ℓ`.
    `𝒫` is any finite set of primes `≤ Q`.
  - **Sanity check (experimentally supported):** I tested the exact formal statement by brute
    force on 20 000 random small instances (Q ≤ 12, N₀ ≤ 200, maximal U). There were no
    violations, and the largest LHS/RHS was 0.985. Script:
    `/tmp/claude-1000/-home-knielsen-code/4e032e48-3e06-5e52-8717-e3a7fde2f526/scratchpad/ls_check.py`
    (scratch; not committed).
- **Lemma 9 needs `H² < ℓ` for each band prime.** The note gets this from `ℓ > Y > H²`. The
  formal version also asks for `ℓ ∤ q`, which holds for `q > ℓ` prime.
- **The capstone takes the giant through its residues.** `G`, `γ` are integers whose images in
  `ZMod p`, `ZMod q` are the stated unit expressions; `p ≠ q`. This is exactly what the
  algorithm's computation mod `N` provides. The FULL failure is a hypothesis, because whether it
  holds is the exceptional-set question (`E_c`, Lemma 10).
- **Lemma 8's sharp form floors `log₂ Q₀`.** This loses at most a factor 2 when `Q₀ ≥ 2`.
- **No pre-registration.** This task runs no experiments. The only computation is the large-sieve
  sanity check above.

## 4. Review

An independent read-only Opus critic of the statements (faithfulness, vacuity, truth of the
large-sieve hypothesis) was launched. Its verdict had not arrived when this was committed:
**review pending**.
