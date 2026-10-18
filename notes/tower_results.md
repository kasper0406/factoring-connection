# tower128: what determines s = p + q mod 128 when N ≡ 1 (mod 8)

Agent `tower128`, 2026-10-18. NEXT_STEPS item 2. All code here is **diagnostic**: it knows `p` and `q` (rule 2).
Status labels follow AGENT_RULES rule 5. No independent critic pass was run on this note. The proofs are
cross-checked by exhaustive brute force (below), not by a second reviewer.

Notation: `p ≠ q` odd primes, `N = pq ≡ 1 (mod 8)` (equivalently `p ≡ q (mod 8)`), `ε = (−1/p) = (−1/q)`,
`p = ε + 4u`, `q = ε + 4v`, `s = p + q`. **QR₂** = "2 is a square mod N" = `[(2/p) = (2/q) = +1]`. In the JSONL files
`QR2 = 1` means QR₂ is true; REPORT.md's `b` is `1 − QR2`. Public key **K0 = (N mod 128, ε, QR₂)**. Hidden bit
**h = bit 2 of (uv mod 8)**. Since `uv mod 4` is a public function of K0 (REPORT §1.1), and `s = ε(N+1) − 16εuv`,
s mod 128 ⇔ h given K0.

## 0. Results at a glance

| # | Claim | Status | Evidence |
|---|---|---|---|
| T1 | `uv ≡ (5 − 2c)·a(c − a) (mod 8)`, with `p ≡ ε·5^a`, `N ≡ 5^c (mod 32)`. For `N ≡ 1 (mod 8)` this depends on `a` only mod 4. So **s mod 128 is a function of (N mod 128, ε, ψ(p))**, `ψ` the even quartic character mod 16 (`ψ = ω²`, `ω` the order-8 character mod 32). `ω(p)` itself is *not* needed. | proved | 0 failures, all 22 666 pairs `p<q<2000` with `p≡q (4)` |
| T1(iv) | Given K0, h is undetermined iff `ψ(p) = ψ(q)` (iff `p ≡ q (mod 16)`). Explicit formula for h in §1. | proved | 0 failures; 5 572 of 11 288 pairs (`N≡1 (8)`, `p<q<2000`) are in that case |
| T2 | For `N ≡ 1 (mod 8)` and `k ≥ 6`, `s mod 2^k` is a function of `(N mod 2^k, p mod 2^{k−3})` | proved | 0 ambiguous keys, `k = 6..11` |
| R1 | When QR₂ holds, `h = [2,2,p*] = [2,2,q*]` (Rédei symbol, Stevenhagen Def. 7.1) on the undetermined pairs | proved | follows from T1(iv) and Def. 7.1 |
| **T3** | `p ≡ q ≡ 7 (mod 8)` (ε = −1, QR₂): `r4 Cl(−4N) = r4 Cl(−8N) = 1` and **`r8 Cl(−4N) ⊕ r8 Cl(−8N) = [r ≡ 7 (mod 16)]`**, where `r ∈ {p,q}` is the prime with `(other/r) = +1`. Hence (C3) **`s ≡ −(N+1) + 64·[N ≡ 1 (16)]·(r8 Cl(−4N) ⊕ r8 Cl(−8N)) (mod 128)`**. | proved (proof uses Stevenhagen Thm 4.1, Def. 4.4, Def. 7.1, Thm 7.9, checked against the arXiv v3 text) | 0 failures on all 9 453 such pairs `11 ≤ p < q < 4000` |
| T4 | `p ≡ q ≡ 1 (mod 8)`, `(p/q) = 1`, `p ≡ q (16)`: with `z_k = [T_k ∈ Cl(−kN)^4]`, `(z4,z8) = (1,1) ⇒ h = 0`, `(1,0)` or `(0,1) ⇒ h = 1`, `(0,0)` undetermined | proved | 760 of 2 064 such pairs decided, 0 failures |
| T4′ | `p ≡ q ≡ 1 (mod 8)`, `(p/q) = −1`: `z4 ⊕ z8 = [N ≡ 9 (16)]` is public, so these depths carry nothing about h | proved | 4 217 pairs, 0 failures |
| E1 | No single tested discriminant (`−4kN`, `k ∈ {1,2,3,5,6,7}`; `−3N`, `−7N`) determines h. Its ranks, types and public-class depths carry ≤ 0.004 bits about h beyond K0, decaying with range. | experimentally supported | §2.3 |
| E2 | In the two QR₂-false cells (half of the undetermined pairs), every tested key, including all 8 discriminants' full 2-Sylow types and all public-class depths jointly, is at its cardinality-matched permutation null | experimentally supported | §2.2 |
| P5 | When `(2/p) = −1`, each Rédei symbol `[a,b,p*]` with `a, b` fixed is asymptotically independent of `ψ(p)` | proved (Chebotarev + a group-order argument) | §3.6 |
| C5 | No 2-power invariant (ranks, types, public-class depths) of these discriminants determines h when QR₂ is false | conjectured | E2, P5 |
| Q | `h` is unrelated to `(p/q)₄(q/p)₄`, `(2/p)₄(2/q)₄`, `(−4/p)₈(−4/q)₈`, `(p/q)₄` | proved for the structural part (h is a Dirichlet character value at p); experimentally supported (null) for the correlations | §3.8 |
| A | None of this is an algorithm | — | §4 |

Coverage at bound 4000 (37 030 pairs, `11 ≤ p < q < 4000`, `N ≡ 1 (mod 8)`). K0 leaves 18 390 pairs undetermined.
T3 resolves all 4 692 in the cell (ε = −1, QR₂). T4 resolves 760 of the 4 096 in (ε = +1, QR₂). That is 5 452 / 18 390
= **29.6%** of the undetermined pairs, all in the QR₂-true half. The QR₂-false half (9 602 pairs) is untouched.

**Correction to REPORT.md §1.1 / NEXT_STEPS item 2 / STATUS (for the orchestrator; not edited by me).** The text
"the missing quantity is the individual value at p of an order-8 Dirichlet character mod 32" is too strong. For
`N ≡ 1 (mod 8)`, only its square `ψ(p) = ω(p)²` (conductor 16, i.e. `p mod 16`) is needed (T1(iii)). The order-8 value
is needed only for `N ≡ 5 (mod 8)` (T1(v)). Also, "r4/r8 of Cl(−4N), Cl(−8N) do not continue the tower" is true
globally but false cell-wise: in the cell (ε = −1, QR₂) the XOR of the two 8-ranks determines s mod 128 exactly (T3).

## 1. Part 1: the exact hidden quantity (proved)

**2-adic characters.** `(Z/32)^* = {±1} × ⟨5⟩` with 5 of order 8. Write `p ≡ ε·5^a`, `q ≡ ε·5^b`, `N ≡ 5^c (mod 32)` with
`a, b, c ∈ Z/8` and `c = a + b`. The sign is ε because `(−1/5^a) = 1`. Let `ω` be the even order-8 character mod 32
with `ω(5) = ζ₈`, and `ψ = ω²` the even quartic character mod 16 (`ψ(5) = i`), so `ψ² = (2/·)`. Then `ω(p) = ζ₈^a`,
`ψ(p) = i^{a mod 4}`, `(2/p) = (−1)^a`. Only `ω(N), ψ(N), (2/N)` are public.

**Theorem T1.** Let `p ≡ q (mod 4)`.
1. `u ≡ ε·a(2a − 1) (mod 8)`.
2. `uv ≡ (5 − 2c)·a(c − a) (mod 8)`, hence `s ≡ ε(N+1) − 16ε(5 − 2c)a(c − a) (mod 128)`.
3. If `N ≡ 1 (mod 8)` (`c` even), `a(c − a) mod 8` depends on `a` only mod 4. So `s mod 128` is a function of
   `(N mod 128, ε, ψ(p))`, and `ω(p)` is not needed.
4. If `N ≡ 1 (mod 8)`, given K0 the bit `h` is undetermined iff `ψ(p) = ψ(q)`, iff `p ≡ q (mod 16)`. Then:
   * `N ≡ 1 (mod 16)` (QR₂ true): `h = [ψ(p) = −1] = [p ≡ ±7 (mod 16)]`;
   * `N ≡ 9 (mod 16)` (QR₂ false): `h = [ψ(p) = i] ⊕ [N ≡ 25 (mod 32)]`, where `[ψ(p) = i] = [p ≡ 5, 11 (mod 16)]`.

   In the other two cells h is public: `N ≡ 1 (16)` with QR₂ false gives `h = [N ≡ 17 (32)]`, and `N ≡ 9 (16)` with QR₂
   true gives `h = 0`.
5. If `N ≡ 5 (mod 8)` (`c` odd), `a(c−a)` changes by 4 when `a ↦ a + 4`, so there `ω(p)` (i.e. `p mod 32`) is needed.

*Proof.*
1. `5^a = (1+4)^a = 1 + 4a + 16·C(a,2) + 64(…)`, so `(5^a − 1)/4 ≡ a + 2a(a − 1) = a(2a − 1) (mod 8)`. This is
   well defined on `a mod 8`, since `(a+8)(2a+15) ≡ a(2a−1) (mod 8)`. Then `4u = p − ε ≡ ε(5^a − 1) (mod 32)`.
2. Put `x = a`, `y = b`. `uv ≡ xy(2x−1)(2y−1) = 4x²y² − 2xy(x+y) + xy` (ε² = 1). Here `4x²y² ≡ 4xy (mod 8)`, and
   `2xy(x+y) ≡ 2xyc (mod 16)` because `x + y ≡ c (mod 8)`. So `uv ≡ xy(5 − 2c)`, with `y ≡ c − a`.
3. For even `c`, `(a+4)(c−a−4) = a(c−a) + 4c − 8a − 16 ≡ a(c−a) (mod 8)`.
4. `a ↦ a+2` changes `a(c−a)` by `2c − 4a − 4 (mod 8)`. For odd `a` this is `2c`, nonzero iff `c ≡ 2 (4)`. For even `a`
   it is `2c − 4`, nonzero iff `c ≡ 0 (4)`. A nonzero change is exactly 4, and multiplying by the odd `5 − 2c` keeps
   it 4, so h flips. Given K0, `a mod 2` (QR₂) and `c` are known, so the hidden datum is `a mod 4`. h depends on it iff
   `2a ≡ c (mod 4)`, iff `a ≡ c − a (mod 4)`, iff `ψ(p) = ψ(q)`. Because ε and `ψ` determine `p mod 16`, this is
   `p ≡ q (mod 16)`. The explicit values follow by evaluating `(5−2c)a(c−a) mod 8` for `c ∈ {0,2,4,6}`.
5. As in 3, the change `4c` is ≡ 4 for odd `c`. ∎

**Lemma T2.** If `N ≡ 1 (mod 8)` and `k ≥ 6`, then `s mod 2^k = f_N(p)` with `f_N(x) = x + N x^{−1} mod 2^k`, and
`f_N(x') = f_N(x)` whenever `x' ≡ x (mod 2^{k−3})`. So `s mod 2^k` is a function of `(N mod 2^k, p mod 2^{k−3})`.
*Proof.* `f(x') − f(x) = (x' − x)(1 − N/(xx'))`. `x' ≡ x (mod 8)` gives `xx' ≡ x² ≡ 1 ≡ N (mod 8)`, so the second factor is
`≡ 0 (mod 8)` and the product is `≡ 0 (mod 2^{(k−3)+3})`. ∎ For `k = 7` this is `p mod 16 ⇔ (ε, ψ(p))`. In general the next
level needs the value at `p` of the order-`2^{k−5}` character mod `2^{k−3}`. This is the residue form of Prop. E2
(trackE).

**Verification** (`src/tower/formula.py`, `experiments/raw/tower128_formula.json`). All pairs `3 ≤ p < q < 2000`:

| Check | Pairs | Failures |
|---|---:|---:|
| index bookkeeping, T1(1), T1(2), `s mod 128` from T1(2) | 22 666 (`p ≡ q (4)`) | 0 |
| `(2/p) = (2/q)`; T1(4) rule for h; `ψ(p) = ψ(q) ⇔ p ≡ q (16)` | 11 288 (`N ≡ 1 (8)`) | 0 |
| T1(3): keys `(N mod 128, ε, ψ(p))` | 128 keys | 0 ambiguous |
| K0 alone | 64 keys | 32 ambiguous |
| T2, `k = 6, …, 11` | up to 17 547 keys | 0 ambiguous |
| T1(5): `N ≡ 5 (8)`, keys with `ψ(p)` vs with `ω(p)` | 128 / 256 keys | 128 / 0 ambiguous |

Undetermined-pair counts by cell (`N ≡ 1 (8)`, `p < q < 2000`): `N≡1 (16)` & QR₂: 2 605; `N≡9 (16)` & ¬QR₂: 2 967. The
other two cells (3 040 and 2 676 pairs) are determined.

## 2. Part 2: extended class-group keys (NEXT_STEPS item 2)

### 2.1 Setup

* **Pairs.** All primes `11 ≤ p < q < B` with `N ≡ 1 (mod 8)`: `B = 2000` gives 11 057 pairs; `B = 4000` gives 37 030.
  `p ≤ 7` is excluded so that `gcd(k, N) = 1` for every multiplier.
* **Discriminants.** `−4kN` for `k ∈ {1,2,3,5,6,7}`, i.e. `−4N, −8N, −12N, −20N, −24N, −28N`. Also `−3N` and `−7N`, the
  only `−kN` (`k ≤ 7`) that are ≡ 0, 1 (mod 4) when `N ≡ 1 (mod 8)`. Primitive forms. `−12N` and `−28N` are
  non-fundamental (conductor 2).
* **Invariants per D.** The exact 2-Sylow type of `Cl(D)` (all `r_{2^j}`). Also, for each public prime `ℓ | D/N` whose
  ramified form `(ℓ, b, c)` is primitive, its **2-depth** `max{j : [form] ∈ Cl(D)^{2^j}}`. Examples: `T₄ = (2,2,(N+1)/2)`
  for `−4N`, `T₈ = (2,0,N)` for `−8N`, `(3,0,N)` for `−12N`.
* **Method (`src/tower/classgroup.py`).**
  1. Enumerate all reduced primitive forms (numpy, exact int64), so `h = |Cl(D)|`.
  2. Write `h = 2^e m`. Adjoin `f^m` for uniformly random forms `f` (seeded) until the generated subgroup has order
     `2^e`. The stopping test is exact, so the 2-Sylow is exact whatever the draws.
  3. Count `|S[2^j]|`.

  Composition and reduction are imported from `src/repair/verify_r4.py` (not edited).
* **Validation.**
  * The form sets equal `verify_r4.reduced_forms` and the ranks `r2, r4, r8` equal `verify_r5.two_power_torsion` on 218
    discriminants (all `|D| < 3000` in a stride, plus `−164, −820`).
  * `Cl(−164) ≅ Z/8`.
  * 17 pytest tests (`tests/test_tower.py`).
  * Largest `h` met: 17 596 at bound 4000 (`|D| ≤ 28N < 4.4·10⁸`).
* **Cost.** On average 0.11 s (B = 2000) and 0.31 s (B = 4000) per pair for all 8 discriminants. Runs: 58 s and 565 s,
  20 workers.
  This is brute force, `Θ(|D|)` per discriminant: a diagnostic, not a method.
* **Proved redundancy.** `Cl(−12N)` and `Cl(−3N)`, and `Cl(−28N)` and `Cl(−7N)`, have isomorphic 2-Sylows. For an order
  of conductor 2, `h(O) = h(O_K)·(2 − (d_K/2))`, and `Cl(O) → Cl(O_K)` is onto with kernel of that order. Here
  `−3N ≡ 5 (8)` gives kernel order 3, and `−7N ≡ 1 (8)` gives order 1. Both are odd, and there are no extra units.
  Observed: 0 differences in 37 030 pairs. So only 6 of the 8 discriminants are distinct for this question.

### 2.2 Ambiguity tables (same convention as `verify_r5`)

Bound 2000 (11 057 pairs) and bound 4000 (37 030 pairs). An ambiguous key is one that occurs with two values of
`s mod 128`.

| Key | B = 2000: keys / ambiguous / pairs in them | B = 4000: pairs in ambiguous keys |
|---|---|---:|
| K0 | 64 / 32 / 5 459 | 18 390 |
| K0 + `(r4,r8,r16)` of `−4N` | 235 / 108 / 5 445 | — |
| K0 + `(r4,r8,r16)` of `−8N` | 240 / 111 / 5 439 | — |
| K0 + `(r4,r8,r16)` of `−20N` | 379 / 168 / 5 417 | — |
| K0 + `(r4,r8,r16)` of `−24N` | 375 / 174 / 5 431 | — |
| K0 + `(r4,r8,r16)` of `−3N` (≡ `−12N`) | 375 / 173 / 5 427 | — |
| K0 + `(r4,r8,r16)` of `−7N` (≡ `−28N`) | 387 / 169 / 5 430 | — |
| K0 + `(r4,r8)` of `−4N` and `−8N` | 279 / 85 / **3 972** | **13 557** |
| K0 + `(r4,r8,r16)` of all 8 | 8 323 / 322 / 1 095 | (overfit, see below) |
| K0 + full types + all public-class depths, all 8 | 9 623 / 174 / 648 | (overfit, see below) |

The last two keys take ~9 000 distinct values on 11 000 pairs, so most pairs are "determined" simply because their key
is unique. To separate information from cardinality, the extension part of each key was **permuted within K0 cells**.
This destroys any link to h but keeps the value multiset. The table compares the pairs left ambiguous with that null
(20 permutations; mean and minimum).

| Key (B = 4000), pairs left ambiguous / null mean | ε=−1, QR₂ (4 692) | ε=+1, QR₂ (4 096) | ε=−1, ¬QR₂ (4 839) | ε=+1, ¬QR₂ (4 763) |
|---|---|---|---|---|
| `(r4,r8)` of `−4N, −8N` | **0** / 4 692 | 3 955 / 4 085 | 4 839 / 4 839 | 4 763 / 4 763 |
| `(r4, depth T)` of `−4N, −8N` | **0** / 4 692 | **3 336** / 4 096 | 4 839 / 4 839 | 4 763 / 4 763 |
| full types + `T` depths, `−4N, −8N` | 0 / 4 527 | 3 066 / 3 460 | 4 813 / 4 817 | 4 755 / 4 751 |
| all 8 discriminants: `(r4,r8,r16)` | 0 / 1 792 | 1 723 / 1 754 | 2 844 / 2 879 | 2 180 / 2 263 |
| all 8: types + all depths | 0 / 476 | 1 372 / 1 390 | 836 / 870 | 1 453 / 1 491 |
| all except `−4N, −8N`: types + depths | 1 702 / 1 712 | 2 182 / 2 207 | 1 755 / 1 809 | 2 376 / 2 399 |

*Reading.*
* In (ε=−1, QR₂) the pair `−4N, −8N` removes all ambiguity, against a null that removes none (T3).
* In (ε=+1, QR₂) the gain over the null is about 400–500 pairs, concentrated in `(p/q) = 1` (T4). The last row shows
  that the other four discriminants add nothing beyond the null in any cell.
* In both ¬QR₂ cells every key sits at its null. The largest gap, 836 vs null min 812 (200 permutations), is not
  significant.

### 2.3 Conditional mutual information `I(h; X | K0)` per invariant

This is the plug-in estimate. The null is a within-cell permutation of `X` (200 permutations). "Excess" = observed − null
mean. `H(h | K0)` is 0.492 bits (B = 2000) and 0.496 bits (B = 4000). Full table:
`experiments/analysis/tower128_tables{,_b4000}.json` (key `cmi`, 63 non-constant invariants; 4 exceed their null p99 at
B = 2000 and 12 at B = 4000, all with excess ≤ 0.0035 bits). The largest excesses:

| Invariant | B = 2000 excess (obs / null p99) | B = 4000 excess (obs / null p99) |
|---|---|---|
| `−8N: (r4,r8,r16)` | 0.0035 (0.0097 / 0.0083) | 0.0019 (0.0038 / 0.0025) |
| `−8N: r8` | 0.0032 (0.0054 / 0.0032) | 0.0016 (0.0022 / 0.0010) |
| `−24N: r8` | 0.0018 (0.0051 / 0.0046) | 0.0012 (0.0023 / 0.0016) |
| `−20N: (r4,r8,r16)` | 0.0017 | 0.0012 |
| every other invariant | < 0.0025 | < 0.0012 |

Every single-discriminant invariant carries at most ~0.004 bits of the 0.49 bits, and the largest excess halves from
B = 2000 to B = 4000. Splitting by range inside the ¬QR₂ cell (ε=−1), the `−8N: r8` signal appears only at `q < 2000`
(0.0060 vs p99 0.0050). For `q ∈ [2000,3000)` it is 0.0011 (below the null mean 0.0021), and for `[3000,4000)` it is
0.0022 (p99 0.0035).

The QR₂-true part of the residual has a clean explanation from T3: `r8(−8N) = r8(−4N) ⊕ [r ≡ 7 (16)]`. So `r8(−8N)` is
informative about h alone only to the extent that `r8(−4N)` is biased. Observed: `P(r8(−4N) = 1) = 0.48–0.49` in that
cell. That is a finite-range bias of the 8-rank distribution, not a leak. Reading: no single-discriminant invariant
determines h, and none shows evidence of an asymptotic information leak (**experimentally supported**). The information
is in the *pair* `(−4N, −8N)`.

## 3. Part 3: theory (Rédei symbols, Stevenhagen arXiv:1806.06250)

### 3.1 Source statements used (checked against arXiv:1806.06250v3, 12 Apr 2021, Math. Proc. Camb. Phil. Soc. 2021)

All of these were read in the v3 PDF by this agent.
* **Thm 2.1:** `r2 = t − 1`, with `t` the number of prime discriminants in `D = ∏ p*_i`, `2* ∈ {−4, ±8}`.
* **(16), Thm 3.1:** Rédei matrix entries `ε_ij = (p*_i / p_j)` (Kronecker for `p_j = 2`), with diagonal = column sum.
  `r4 = r2 − rank R4`. The column of `p_j` is `γ(ϕ4([p_j]))`, so a zero column means `[p_j] ∈ 2C`.
* **Lemma 4.2:** a decomposition `D = d1 d2` is of the 2nd kind iff `(D/d_i / ℓ) = 1` for all primes `ℓ | d_i`.
* **Thm 4.1, (22), (23):** `r8 = r4 − rank R8`, with `R8` entries `ψ_i([m_j])`, `[m_j]` spanning `C[2] ∩ 2C`, and `ψ_i`
  quartic characters spanning `Ĉ[4]/Ĉ[2]` with `2ψ_i` the 2nd-kind genus characters. Also `[m] ∈ 4C ⇔ ϕ8([m]) = 0 ⇔`
  all `ψ_i([m]) = 0`.
* **Def 4.4:** `[d1,d2,m] = Art(m, H₂F/K) ∈ F₂`, for `F` the unramified cyclic quartic extension of `K = Q(√D)` through
  `Q(√d1,√d2)`.
* **Def 7.1, (35)–(38):** for `a` with `(a,−1)_ℓ = 1` everywhere, `[a,a,c] = ψ_a(c)`, with `ψ_2` "a character of order
  4 on (Z/16Z)^*", unique up to twisting by `χ_{−1}, χ_2`.
* **Def 7.6, Def 7.8, Lemma 7.7, Cor 8.2:** the general symbol is the Artin symbol in a minimally ramified `F_{a,b}`,
  independent of the choice. When one of `Δ(a), Δ(b)` is ≡ 1 (mod 8) and the other conditions hold, minimally
  ramified means unramified.
* **Thm 7.9 (Rédei reciprocity):** for `a, b, c` satisfying (1) `(a,b)_ℓ = (a,c)_ℓ = (b,c)_ℓ = 1` for all `ℓ ≤ ∞` and (2)
  `gcd(Δa, Δb, Δc) = 1`, the symbol is linear in each argument and invariant under all permutations.
* Example **(28)–(29):** for `p ≡ 1 (8)`, `r8 Cl(−4p) = 1` iff `p` splits completely in `Q(ζ₈, √(1+i))`, and
  `[−4,p,2] = [−1,2,p]`.

Linearity is used below only where all three symbols involved are verified to be defined.

### 3.2 R1: the hidden bit is a Rédei symbol (when QR₂ holds)

For `p ≡ ±1 (mod 8)`, `p* ≡ 1 (mod 8)`. The pair `(2, p*)` satisfies (37): `gcd(8, p) = 1`, `(2,2)_ℓ = (2,−1)_ℓ = 1`,
`(2,p*)_p = (2/p) = 1`, and `(2,p*)_2 = 1` since `p* ≡ 1 (8)`. So `[2,2,p*] = ψ₂(p*)`, which is well defined because
`χ_{−1}(p*) = χ_2(p*) = 1`. As `ψ₂(9) = ψ₂(5)² = −1`, `[2,2,p*] = [p* ≡ 9 (mod 16)] = [p ≡ ±7 (mod 16)] = [ψ(p) = −1]`.
With T1(4) this gives **`h = [2,2,p*] = [2,2,q*]`** on the undetermined QR₂-true pairs (**proved**). For a negative
argument `c = −p`, the even and odd minimal choices of `ψ₂` agree, because `c ≡ 1 (mod 4)`.

### 3.3 Theorem T3 (ε = −1, QR₂ true: p ≡ q ≡ 7 mod 8) — proved

Exactly one of `(q/p)`, `(p/q)` is `+1`, because `p ≡ q ≡ 3 (4)`. Call that prime `r` and the other `r'`, so
`(r'/r) = 1`.

*`D = −4N`.* Prime discriminants `−4, −r, −r'`, so `t = 3` and `r2 = 2`. Column of 2 in `R4`: `(−r/2) = (−r'/2) = +1`
because `−r ≡ 1 (8)`. It is zero, so `T₄ = [𝔭₂] ∈ 2C`. The columns of `r` and `r'` are both `(1, 0, 1)ᵀ`, from
rows `(−4, −r, −r')`. For column `r`: `(−4/r) = −1` and `(−r'/r) = (−1/r)(r'/r) = −1`, so the diagonal is 0. For column
`r'`: `(−4/r') = −1` and `(−r/r') = −(r/r') = (r'/r) = +1`, so the diagonal is 1. Hence `rank R4 = 1`, `r4 = 1`,
`ker R4 = ⟨e₂, e_r + e_{r'}⟩`, and `α(e_r + e_{r'}) = [(√−N)] = 1`. So `C[2] ∩ 2C = ⟨T₄⟩`, where `T₄ ≠ 1` because
`(2,2,(N+1)/2)` is reduced with `a = 2`.

The unique 2nd-kind decomposition is `(−r, 4r')`. Lemma 4.2: `(4r'/r) = 1`; `(−r/2) = 1`; `(−r/r') = (r'/r) = 1`. By
Thm 4.1, `r8 = 1 ⇔ T₄ ∈ 4C ⇔ A := [−r, 4r', 2] = 0`.

*`D = −8N`.* Prime discriminants `−8, −r, −r'`. The column of 2 is zero as before. `(−8/r) = (−2/r) = −1`, so
`r4 = 1`, and the same kernel computation gives `C[2] ∩ 2C = ⟨T₈⟩`: `α(e_r + e_{r'}) = [𝔭₂]` because
`(√−2N) = 𝔭₂𝔭_r𝔭_{r'}`. The decomposition `(−r, 8r')` is of the 2nd kind: `(8r'/r) = (2/r)(r'/r) = 1`, `(−r/2) = 1`,
`(−r/r') = 1`. So `r8 = 1 ⇔ T₈ ∈ 4C ⇔ A' := [−r, 8r', 2] = 0`.

*Linking the two.* In `Q*/Q*²`, `8r' = 2·r'·(square)` and `4r' = r'·(square)`. The symbols `[−r, r', 2]`, `[−r, 2, 2]` and
`[−r, 2r', 2]` are all defined (conditions (1), (2) of Thm 7.9):
* Hilbert symbols at `r`: `(−r, r')_r = (r'/r) = 1`, `(−r, 2)_r = (2/r) = 1`.
* At `r'`: `(−r, r')_{r'} = (−r/r') = 1`, `(r', 2)_{r'} = (2/r') = 1`.
* At 2: `(−r, r')_2 = 1` since `−r ≡ 1 (mod 4)`. Also `(2, x)_2 = (−1)^{(x²−1)/8} = 1` for `x ∈ {−r, r'}`, since
  `−r ≡ 1` and `r' ≡ 7 (mod 8)`, and `(2,2)_2 = 1`.
* At ∞: only `−r` is negative.
* Discriminants: `gcd(Δ(−r), Δ(r'), Δ(2)) = gcd(r, 4r', 8) = 1` and `gcd(r, 8, 8) = gcd(r, 8r', 8) = 1`.

Def 7.6: `gcd(Δ(−r), Δ(r'))` and `gcd(Δ(−r), Δ(2r'))` are 1 and `Δ(−r) = −r ≡ 1 (mod 8)`, so the minimally ramified
`F` is unramified. Therefore the Def 7.8 symbols with `c = 2` equal the Def 4.4 symbols `[d1, d2, 2]` at the ideal
`𝔭₂` of norm 2.

By Thm 7.9, `A' = [−r, r', 2] + [−r, 2, 2] = A + [2, 2, −r]`. By §3.2 applied to `c = −r`,
`[2, 2, −r] = [−r ≡ 9 (16)] = [r ≡ 7 (mod 16)]`. Hence
**`r8 Cl(−4N) ⊕ r8 Cl(−8N) = [T₄ ∈ C⁴] ⊕ [T₈ ∈ C⁴] = A ⊕ A' = [r ≡ 7 (mod 16)]`**. ∎

*Corollary C3.* On undetermined pairs `ψ(r) = ψ(r')`, so `h = [r ≡ 7 (16)]`. With `uv mod 4 = 0` for QR₂ true and
ε = −1, this gives `s ≡ −(N+1) + 64·[N ≡ 1 (16)]·(r8 Cl(−4N) ⊕ r8 Cl(−8N)) (mod 128)` for **all** `p ≡ q ≡ 7 (mod 8)`.
Verified: all 9 453 pairs `11 ≤ p < q < 4000`. The `r4 = 1` claim, `r8 = [T ∈ C⁴]`, the XOR rule and C3 each have
0 failures.

### 3.4 Theorems T4, T4′ (ε = +1, QR₂ true: p ≡ q ≡ 1 mod 8) — proved

*Case `λ = (p/q) = 1`.* For both `D = −4N` and `D = −8N`, every entry of `R4` is `+1`, so `R4 = 0`, `r4 = 2` and
`C[2] ∩ 2C = C[2] ∋ T`. The decompositions `(p, −4q)` and `(q, −4p)` (respectively `(p, −8q)`, `(q, −8p)`) are of the
2nd kind and give a basis of `Ĉ[4]/Ĉ[2]`. Thus `[T₄ ∈ 4C] ⇔ A = B = 0`, where `A = [p,−4q,2]` and `B = [q,−4p,2]`.
Likewise `[T₈ ∈ 4C] ⇔ A' = B' = 0`.

Linearity applies: `[p,−q,2]`, `[p,2,2]` and `[p,−2q,2]` are all defined, with `(−q,2)_2 = 1` because `−q ≡ 7 (8)`, and
`Δ(p) ≡ 1 (8)` gives unramified `F`. So `A' = A + [2,2,p]` and `B' = B + [2,2,q]`. On undetermined pairs
`[2,2,p] = [2,2,q] = h` (R1). Hence:
* `z4 = z8 = 1` ⇒ `A = B = 0` and `A' = B' = 0` ⇒ `h = 0`;
* exactly one of `z4, z8` is 1 ⇒ the other pair is `(h, h) ≠ (0,0)` ⇒ `h = 1`;
* `z4 = z8 = 0` is compatible with both values of h.

Under the heuristic that `(A, B, h)` are uniform and independent, the decided fraction is 3/8. Observed: 760 / 2 064 =
0.368, with 0 failures. The `r4 = 2` claim has 0 failures. *Corollary C4:* on decided pairs,
`s ≡ N + 1 + 64h (mod 128)` (0 failures).

*Case `λ = −1` (T4′).* Here `r4 = 1` for both discriminants, the only 2nd-kind decompositions are `(−4, pq)` and
`(−8, pq)`, and `C[2] ∩ 2C = ⟨T⟩`. Linearity gives `[−8, pq, 2] = [−1, pq, 2] + [2, 2, pq]`, so
`z4 ⊕ z8 = [2,2,N] = [N ≡ 9 (mod 16)]`, which is **public**. The ψ-dependence cancels in the product. 4 217 pairs,
0 failures.

### 3.5 Consistency checks against brute force (single primes, `src/tower/symbols.py`)

The same Rédei manipulations predict three single-prime statements:
* `p ≡ 1 (8)`: `r8 Cl(−4p) = [(−4/p)₈ = 1]` (Stevenhagen (28); Barrucand–Cohn);
* `p ≡ 1 (8)`: `r8 Cl(−8p) = [(2/p)₄ = 1]`, from `R8 = [p,−8,2] = [−1,2,p] + [2,2,p]`;
* `p ≡ 7 (8)`: `r8 Cl(−8p) = [p ≡ 15 (16)]`, from `R8 = [−p,2,2]`.

They also predict the prototype of T3/T4: `[T₄ ∈ C⁴] ⊕ [T₈ ∈ C⁴] = [p ≡ 9 (16)]` for `p ≡ 1 (8)`. The classical identity
**S1** `(−4/p)₈ = (2/p)₄·(−1)^{(p−1)/8}` is the Frobenius of `√(1+i) = 2^{1/4}ζ₁₆` over `Q(ζ₈)`. All were checked for
every prime `p < 20 000` in the relevant class (556 primes ≡ 1 (8); 565 primes ≡ 7 (8)), with 0 failures. This also
tests the dictionary "depth of `T` ↔ Rédei symbol at 2" against the brute-force class groups.

### 3.6 Why the QR₂-false half is not reached (P5 proved; C5 conjectured)

When QR₂ is false (`p ≡ q ≡ ±3 (mod 8)`), the undetermined pairs need `ψ(p) ∈ {i, −i}`: an element of **order 4** in
`Gal(Q(ζ₁₆)⁺/Q)`, with `p` inert in `Q(√2)`. The T3/T4 mechanism dies at three points:
1. `T ∉ 2C`, since the column of 2 in `R4` is nonzero (`(p*/2) = −1`), so `T` has depth 0;
2. `[2,2,p*]` is undefined, since `(2,p*)_p = (2/p) = −1`;
3. more generally, any Rédei symbol with `2` in one argument and `p` (to odd power) in another needs `(2·x/p) = 1`,
   which forces a compensating hidden factor `x` (e.g. `q*` with `(q*/p) = −1`).

**Proposition P5 (proved).** Let `a, b` be fixed (public), `(2/p) = −1`, and suppose `[a,b,p*]` is defined. Then, over
primes `p` in a fixed class mod 8, `[a,b,p*]` and `ψ(p)` are asymptotically independent.

*Proof.* A defined symbol needs `(a/p) = (b/p) = 1`, so `p` splits in `E = Q(√a,√b)`. Then `√2 ∉ E`, since `p` is inert
in `Q(√2)`. `[a,b,p*]` is the Frobenius of `p` in `F = F_{a,b}` (Def 7.8; it lies in `Gal(F/E)`, central of order 2).

`F ∩ Q(ζ₁₆)⁺` is `Q`, `Q(√2)` or `Q(ζ₁₆)⁺`. The last is impossible:
* if `a ≠ b`, `Gal(F/Q) ≅ D₄` has no cyclic quotient of order 4;
* if `a = b`, `F_{a,a} = Q(ζ₁₆)⁺` forces `a = 2`, which needs `(2/p) = 1`.

`√2 ∈ F` would make `F = E(√2)` abelian, contradicting the previous step, so the intersection is `Q`. By Chebotarev in
`F·Q(ζ₁₆)⁺`, the pair `(Frob_p|F, Frob_p|Q(ζ₁₆)⁺)` is asymptotically equidistributed on the product group. The same
holds for any finite family of such symbols: their Frobenius lies in the elementary abelian `Gal(L/E₁⋯E_k)`, whose image
in `Z/4` cannot contain `ψ(p)` of order 4, so again `L ∩ Q(ζ₁₆)⁺ ⊆ Q(√2)`. ∎

**Conjecture C5.** No 2-power invariant of the tested discriminants (ranks, types, public-class depths, or any
combination) determines h, or carries asymptotically positive information about it, in the QR₂-false cells.

P5 does not cover symbols in which `p` occurs in two arguments (quartic residue symbols `[p*,p*,c] = (c/p*)₄`), nor
16-rank data beyond Rédei symbols, nor the coupling through `q`. Evidence: E2 (§2.2) and the CMI table (§2.3).

**Why single discriminants fail when QR₂ holds** (proved form). In each discriminant, `ψ(p)` enters only through `A'`
(`= A + ψ`), or only through `A` and `B`. The other quantity, `A = [−r,4r',2] = [2,−r,r']` (the Frobenius of `r'` in
`F_{2,−r}`), is an independent hidden Rédei symbol. So one discriminant gives `ψ ⊕ (mask)`, and only the pair cancels
the mask. That the mask is asymptotically unbiased (so single-discriminant MI → 0) is the 8-rank equidistribution
expected from governing fields (Stevenhagen §9): **conjectured** here, consistent with §2.3.

**What `r4` can see (proved).** For every tested `D` (fundamental, or with the same 2-Sylow as a fundamental one, §2.1),
`R4` has entries `(ℓ*/p), (ℓ*/q)` for `ℓ* ∈ {−4, ±8, −3, 5, −7}`, entries `(p*/ℓ), (q*/ℓ)`, and `(p*/q), (q*/p)`. So
`r4` is a function of `(p mod 840, q mod 840, (p/q))`, which sees `p mod 16` only through `p mod 8`. By Dirichlet's
theorem (moduli `16·105·q` for fixed `q`), `p mod 16` is asymptotically independent of `(p mod 105, (p/q))` within a
class mod 8. So `r4` carries asymptotically no information about h given K0. Measured: `r4` excess ≤ 0.0006 bits.

### 3.7 Relation to quartic reciprocity and `(p/q)₄(q/p)₄`

* `h` is the value at `p` (equivalently at `q`, since `ψ(p) = ψ(q)` on undetermined pairs) of a **Dirichlet character
  of conductor 16**. It is not a mutual symbol between `p` and `q`.
* By S1, for `p ≡ 1 (8)`, `ψ(p) = (−4/p)₈·(2/p)₄`: a ratio of two residue symbols with fixed numerators at the hidden
  prime individually.
* Their products over the factors lose it. On undetermined pairs, `(−4/p)₈(−4/q)₈ = (2/p)₄(2/q)₄` exactly (proved:
  `ψ(p)ψ(q) = 1`; the two contingency tables are identical).
* Burde-type rational quartic reciprocity expresses `(p/q)₄(q/p)₄` through representations `p = a² + b²`,
  `q = c² + d²`. Its sign term `(−1)^{(q−1)/4}` sees only `q mod 8`. This statement is from memory: **UNVERIFIED**, no
  source read.
* Measured on the 570 pairs `p ≡ q ≡ 1 (8)`, `(p/q) = 1`, `p ≡ q (16)`, `q < 2000`. `I(h; ·|N mod 128)` for
  `(p/q)₄(q/p)₄`, `(2/p)₄(2/q)₄`, `(−4/p)₈(−4/q)₈` and `(p/q)₄` is 0.0014, 0.0116, 0.0116 and 0.0080 bits, against a
  permutation null mean of about 0.010. All are at the null, and h is not a function of any of them.

So s mod 128 is **not** governed by `(p/q)₄(q/p)₄`-type products. It is governed by 8-rank data of two discriminants
jointly: the Artin symbol of the prime above 2 in the 4-Hilbert class fields of `Q(√−N)` and `Q(√−2N)`.

## 4. Is any of this an algorithm? No.

* **T3 and T4 are structural identities, not procedures.** They say h equals `r8 Cl(−4N) ⊕ r8 Cl(−8N)` (or a function
  of `[T ∈ C⁴]`). Using them needs either the 8-ranks or a 4th-power membership test for the public classes `T₄, T₈`
  in groups of unknown order.
* **Via Rédei matrices.** Building `R4` and `R8` needs the prime discriminants `p*, q*`, i.e. the factors. The `R8`
  entries further need a solution of `x² − d1y² − d2z² = 0` with `d1, d2` built from the factors (Stevenhagen §5).
* **Via the class group.** No method is known to compute `r8` or "`T ∈ C⁴`" from an unfactored discriminant faster than
  computing the class group. Subexponential `L_{|D|}[1/2]` index-calculus methods (Hafner–McCurley; not re-read here)
  produce relations from which ambiguous forms, and hence the factors of `N`, follow anyway. The project's brute force
  is `Θ(|D|)` per discriminant.
* **The cell condition is itself hidden.** T3 applies when ε = −1 and QR₂ holds. Deciding ε is the SOS problem and
  deciding QR₂ is the 4-rank problem (REPORT §1.1).
* **No public `(N mod 2^k)` decision follows.** The information floor is unchanged.
* **The next level.** Exploratory, B = 2000, no null: within (ε=−1, QR₂), the full depths and types of `T₄, T₈` leave
  1 210 of 1 482 pairs ambiguous for `s mod 256`. So the depth of `T` does not by itself continue the tower to the next
  bit.

### Candidate card (rule 9)

```text
ID: TOWER-128-R8PAIR.  The bit of s mod 128 beyond (N mod 128, eps, QR_2) equals r8 Cl(-4N) XOR r8 Cl(-8N) when
    p = q = 7 (mod 8) (T3); partially via [T in C^4] when p = q = 1 (mod 8) (T4).
Target / input family: N = pq = 1 (mod 8); the eps = -1, QR_2 cell fully; eps = +1, QR_2, (p/q) = 1 for 3/8 of it.
Representation from N: the form class groups of -4N and -8N and the public classes T4 = (2,2,(N+1)/2), T8 = (2,0,N).
Extraction: 8-ranks, or 4th-power membership of T4, T8.
Why it exposes information: Redei linearity, [-r, 8r', 2] + [-r, 4r', 2] = [2, 2, -r] = psi_2(-r) (Stevenhagen Def 7.1, Thm 7.9).
Known lemmas: Stevenhagen Thm 3.1, Thm 4.1, Def 4.4, Def 7.1, Thm 7.9 (checked against v3).  New: T1-T4 (proved here).
Correctness: proved; exhaustive check 0 failures (9 453 + 2 064 + 4 217 pairs, q < 4000).
Cost in n: computing r8 / C^4-membership without the factors: no method below class-group computation (L[1/2]),
    which already factors N.  Brute force here: Theta(|D|).
Randomness: none needed for the identity.
Closest prior method: Redei-Reichardt 8-rank theory, Barrucand-Cohn / Kaplan single-prime 8-divisibility results.
    Difference: two multiplier discriminants combined so the hidden mask cancels.
Smallest informative experiment: done (sections 2-3).
Refutation: one pair p = q = 7 (mod 8) violating C3 (none found).  For usefulness: a poly-time r8 or C^4 test (none known).
Repair: none for an algorithm; the open question is the QR_2-false half (C5) and s mod 256.
```

## 5. Reproduction

All runs use `.venv/bin/python` from the project root. Config: `experiments/configs/tower128.json`.

```bash
.venv/bin/python -m src.tower.formula --bound 2000 --out experiments/raw/tower128_formula.json                 # Part 1, <1 s
.venv/bin/python -m src.tower.run_tower --bound 2000 --workers 20 --out experiments/raw/tower128_invariants.jsonl  # 58 s
.venv/bin/python -m src.tower.run_tower --bound 4000 --workers 20 --out experiments/raw/tower128_invariants_b4000.jsonl  # 565 s
gzip -9 experiments/raw/tower128_invariants.jsonl experiments/raw/tower128_invariants_b4000.jsonl     # as committed
.venv/bin/python -m src.tower.analyze --inp experiments/raw/tower128_invariants.jsonl.gz --out experiments/analysis/tower128_tables.json        # ~1 min
.venv/bin/python -m src.tower.analyze --inp experiments/raw/tower128_invariants_b4000.jsonl.gz --out experiments/analysis/tower128_tables_b4000.json  # ~4 min
.venv/bin/python -m src.tower.symbols --bound 2000 --single_bound 20000 --out experiments/raw/tower128_symbols.json  # S1-S3, 3 s
.venv/bin/python -m pytest -q tests/test_tower.py                                                              # 17 tests, ~25 s
```

| Kind | Files |
|---|---|
| Code | `src/tower/{formula,classgroup,run_tower,analyze,symbols}.py` |
| Raw | `experiments/raw/tower128_{formula,symbols}.json`, `experiments/raw/tower128_invariants{,_b4000}.jsonl.gz` (+ `.log`) |
| Tables | `experiments/analysis/tower128_tables{,_b4000}.json` (+ `.log`) |

Leakage: all of this is diagnostic code. Nothing inverts modulo `N`, so rule 3 does not arise.
