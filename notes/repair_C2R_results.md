# C-2R repair cycle: two-layer interval folding (results, agent `repair`, 2026-10-17)

Pre-registration: `notes/repair_C2R_prereg.md`, committed before any instance run (commit `a8920b7`).
Design source: `notes/codex_r3_critic.md` §4. Status labels follow AGENT_RULES §5. All arithmetic is exact
(integers mod p, numpy int64 with values < 2^31). Everything here is DIAGNOSTIC: p is known.

## 0. Verdict

| Question | Answer | Status |
|---|---|---|
| Does any non-character exponent reach `d ≤ 3L/8` on the 16 largest instances? | **No: 0/16** (0/32 overall). On all 16 the non-character arm gives no improvement at all over the one-layer `m+1`. | experimentally supported (exhaustive, exact) |
| Audit numbers (0/32, 14/32, 0/32) | **Reproduced exactly** by an independent implementation | experimentally supported |
| Pre-registered decision | **Reject** H-C2R, the order-independent constant-factor two-layer fold, for this model and range | refuted for this model |
| Frozen extension X1 (additive mix of two public branches) | **0/16** on the large instances (prediction 0/16). The 8 instances `p ≤ 409` pass degenerately, as pre-declared. | refuted for this model |
| Why `E ≡ (p+1)/2` compresses | First layer = signed identity `y ↦ χ(y)y`, image `A = H ∩ T` inside the interval; the second-layer reflection pairs squares | proved (mechanism); ratio `d/m → 3/4` experimentally supported |
| Non-character exponents have `o(L)` internal pairs? | **Proved** for the squaring exponents (`E ≡ 2`, and `E ≡ 2+(p−1)/2` for `p ≡ 1 mod 8`), for every `c`, exact or not. **False without exactness** for twisted identities (`(p−1) \| k(E−1)`, `k ≥ 3`): `Θ(L/k²)` pairs. **Open** for exact folds in general, but no exact non-character fold with ≥ 1 pair exists at any instance `p ≥ 4099`. | mixed; see §3 |

Strongest supported statement: *in the tested two-layer power-map model, every compression
beyond the one-layer degree that survives exactness comes from an exponent `E` with `(p−1) | 2(E−1)`,
i.e. an exponent that already encodes a multiple of `p−1` and factors `N` by itself.* This is proved for
`p ≡ 3 mod 4` as a mechanism and experimentally exhaustive on the 32 instances. As an asymptotic
classification it is conjectured (§3.5).

## 1. Implementation and validation

`src/repair/c2r_search.py` is my own implementation. It vectorizes over translations and checks exactness for every
pair-sum `c` of every gcd-2 exponent. Validation (`tests/test_repair_c2r.py`, plus a scratch run before the prereg):

- **Reflection criterion = brute force.** For every prime `11 ≤ p ≤ 83`, every gcd-2 `E` and every `v ∈ F_p`, exactness and `d`
  from the criterion equal direct evaluation of `Q(φ(x))` on all of `F_p`: 14 916 cases, 0 mismatches.
- **Pair sums cover every relevant translation.** For `p ≤ 83` the minimum exact `d` over *all* `(E, v)`, by brute force,
  equals the search result. No exponent with `gcd(E, p−1) ≠ 2` ever gives an exact fold with `d < m+1`, as the audit asserts:
  gcd 1 cannot beat `⌈L/2⌉`, and gcd ≥ 3 is never exact by O2.
- **Cross-check against the audit's code.** The audit's code is extracted from the note at test time, not copied. On all 143 primes
  101–997 and then on all 32 instances, the best `d` and the pass flags agree exactly.
- **Certificate** `p = 1031, L = 33, φ(x) = ((x+16)^{516} − 7)^2`: the zero set of `∏_{r∈R}(φ(x) − r)` over all 1031
  inputs is exactly `Z_33`. Also `φ(Z_33) = R = {1,4,9,16,25,36,49,64,81,196,441}` and `deg Q = 11`. Note `516 = (p+1)/2`.

## 2. Results (pre-registered runs)

Arms: 1 = public `E ≤ 64`, 2 = all `E`, 3 = all except `2E−2 ≡ 0 mod (p−1)`. `*` = passes `8d ≤ 3L`.
D1 = max internal pairs `r*(c)` over **all** `c` (non-character `E`), with the arg-max `E` and twist order `k`
(`k = 0`: no `k ≤ 12` with `(p−1) | k(E−1)`). D2 = the same over **exact** `c`. X1 = extension of §6 of the prereg
(min `d` over all `c`, exact or not, and whether an exact fold with `8d ≤ 3L` exists).

| p | p mod 4 | L | m+1 | 3L/8 | arm1 d | arm2 d (E) | arm3 d | D1 (E; k) | D2 | X1 min d | X1 exact pass |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 257 | 1 | 17 | 9 | 6.38 | 8 | 7 (70) | 7 | 2 (2; 0) | 2 | 3 | yes |
| 263 | 3 | 17 | 9 | 6.38 | 7 | 6* (132) | 7 | 3 (58; 0) | 2 | 4 | yes |
| 269 | 1 | 17 | 9 | 6.38 | 7 | 7 (22) | 7 | 3 (98; 0) | 2 | 4 | yes |
| 271 | 3 | 17 | 9 | 6.38 | 7 | 6* (136) | 7 | 3 (62; 0) | 2 | 4 | yes |
| 389 | 1 | 19 | 10 | 7.12 | 8 | 8 (22) | 8 | 2 (2; 0) | 2 | 5 | yes |
| 397 | 1 | 19 | 10 | 7.12 | 9 | 8 (74) | 8 | 2 (2; 0) | 2 | 4 | yes |
| 401 | 1 | 21 | 11 | 7.88 | 9 | 9 (54) | 9 | 3 (166; 0) | 2 | 4 | yes |
| 409 | 1 | 21 | 11 | 7.88 | 11 | 9 (370) | 9 | 3 (218; 0) | 2 | 5 | yes |
| 1031 | 3 | 33 | 17 | 12.38 | 17 | 11* (516) | 16 | 3 (34; 0) | 1 | 9 | no |
| 1033 | 1 | 33 | 17 | 12.38 | 17 | 16 (478) | 16 | 3 (26; 0) | 1 | 10 | no |
| 1039 | 3 | 33 | 17 | 12.38 | 17 | 12* (520) | 16 | 3 (16; 0) | 1 | 11 | no |
| 1049 | 1 | 33 | 17 | 12.38 | 17 | 16 (234) | 16 | 3 (30; 0) | 1 | 11 | no |
| 1543 | 3 | 39 | 20 | 14.62 | 20 | 14* (772) | 20 | 4 (1378; 0) | 0 | 14 | no |
| 1549 | 1 | 39 | 20 | 14.62 | 20 | 20 | 20 | 5 (1162; 4) | 0 | 14 | no |
| 1553 | 1 | 39 | 20 | 14.62 | 20 | 19 (1242) | 19 | 4 (790; 0) | 1 | 14 | no |
| 1559 | 3 | 39 | 20 | 14.62 | 20 | 18 (780) | 20 | 4 (178; 0) | 0 | 13 | no |
| 4099 | 3 | 65 | 33 | 24.38 | 33 | 23* (2050) | 33 | 4 (28; 0) | 0 | 25 | no |
| 4111 | 3 | 65 | 33 | 24.38 | 33 | 23* (2056) | 33 | 5 (422; 0) | 0 | 24 | no |
| 4127 | 3 | 65 | 33 | 24.38 | 33 | 23* (2064) | 33 | 4 (70; 0) | 0 | 26 | no |
| 4129 | 1 | 65 | 33 | 24.38 | 33 | 33 | 33 | 4 (98; 0) | 0 | 25 | no |
| 6151 | 3 | 79 | 40 | 29.62 | 40 | 28* (3076) | 40 | 4 (2; 0) | 0 | 28 | no |
| 6163 | 3 | 79 | 40 | 29.62 | 40 | 27* (3082) | 40 | 5 (2918; 0) | 0 | 32 | no |
| 6173 | 1 | 79 | 40 | 29.62 | 40 | 40 | 40 | 6 (4630; 4) | 0 | 33 | no |
| 6197 | 1 | 79 | 40 | 29.62 | 40 | 40 | 40 | 6 (1550; 4) | 0 | 32 | no |
| 16411 | 3 | 129 | 65 | 48.38 | 65 | 46* (8206) | 65 | 7 (13676; 6) | 0 | 57 | no |
| 16417 | 1 | 129 | 65 | 48.38 | 65 | 65 | 65 | 5 (7022; 0) | 0 | 56 | no |
| 16421 | 1 | 129 | 65 | 48.38 | 65 | 65 | 65 | 8 (4106; 4) | 0 | 55 | no |
| 16427 | 3 | 129 | 65 | 48.38 | 65 | 49 (8214) | 65 | 5 (356; 0) | 0 | 57 | no |
| 24593 | 1 | 157 | 79 | 58.88 | 79 | 79 | 79 | 6 (2; 0) | 0 | 70 | no |
| 24611 | 3 | 157 | 79 | 58.88 | 79 | 56* (12306) | 79 | 6 (2; 0) | 0 | 70 | no |
| 24623 | 3 | 157 | 79 | 58.88 | 79 | 54* (12312) | 79 | 6 (2; 0) | 0 | 70 | no |
| 24631 | 3 | 157 | 79 | 58.88 | 79 | 58* (12316) | 79 | 6 (2; 0) | 0 | 70 | no |

**Counts.** Arm 1: 0/32. Arm 2: 14/32. Arm 3: 0/32. These are identical to the audit. All 14 arm-2 passes use
`E = (p+1)/2` with `p ≡ 3 (mod 4)` (16 such instances; 1559 and 16427 miss by 3.4 and 0.6). On the 16 largest
instances arm 3 is exactly `m+1` everywhere. **D2 = 0 for every `p ≥ 4099`** (and for 1543, 1549, 1559; D2 = 1 at 1553): no exact non-character fold
has even one internal pair at the 16 largest instances. The non-exact maxima D1 are small (≤ 8) and come from `E = 2` (sums of two squares, e.g.
`c = 5525 = 5²·13·17`, 6 representations) or from twisted identities `k = 4, 6` (§3.4).

**X1.** 0/16 large instances have an exact fold with `8d ≤ 3L`. In two of them the non-exact minimum went below the threshold
(4111: 13 candidates; 6151: 477 candidates), and all failed exactness on `F_p`. At `p ≈ 24600` the best non-exact
`d = 70 > 58.9` comes from `(E1, E2) = (2, 4)`. That is a three-layer squaring fold whose pairs are sums of two squares (Thm C2R-3).
The 8 instances with `p ≤ 409` pass: degree-64 branches against `p ≈ 300` have fibres of size up to 64, which is the degeneracy the prereg excluded.

**Cost.** Main search 38 CPU-s; X1 4 899 CPU-s; theory and audit scripts < 300 CPU-s. Total ≈ 1.45 CPU-h of the 2 h cap.
No instance timed out.

## 3. Mathematics

Notation: `m = ⌊√p⌋//2`, `T = [−m, m]`, `H = {0} ∪ (F_p^*)^2`, `χ` the Legendre symbol, `A = A_E = {j^E : 0 ≤ j ≤ m}`
(`|A| = m+1` when `gcd(E, p−1) = 2`), `r*(c) = #{{a,b} ⊂ A : a ≠ b, a+b = c}`, `S_A(c) = Σ_{a∈A} χ(c−a)`.

### 3.1 Lemma C2R-1 (proved)

Let `gcd(E, p−1) = 2` and `c = −2v`.
(a) For every fold, exact or not, `d = |A| − r*(c)`. *Proof:* `(a+v)^2 = (b+v)^2 ⇔ a = b` or `a + b = c`, and given `c` each `a` has at most one partner. ∎
(b) If the fold is exact, then `4 r*(c) = |A| + S_A(c) + [c ∈ A] − 2[c/2 ∈ A]`.
*Proof:* exactness gives `{a ∈ A : c−a ∈ H} = {a ∈ A : c−a ∈ A}`. The right side has `2r*(c) + [c/2 ∈ A]` elements. The left side has
`Σ_{a∈A}(1+χ(c−a))/2 + [c ∈ A]/2` elements. ∎ (Checked on 20 exact random cases, 0 violations.)
**Consequences.** `d ≤ 3L/8 ⇔ r*(c) ≥ (2m+5)/8`. For exact folds this is equivalent to `S_A(c) ≥ 3/2 − [c∈A] + 2[c/2∈A]`.
So an exact fold beats the threshold iff the character sum over the reflected set `c − A` is *non-negative*. An exact fold with
`o(m)` pairs needs `S_A(c) = −(1 − o(1))|A|`, i.e. `c − A` lies almost entirely in the non-squares. Exactness, not compression, is the binding constraint.
For a "random-looking" `A` of size `≈ √p/2`, `|A ∩ (c − A)| = O(log p)`. Exactness then forces `c − a ∉ H` for essentially every `a`,
an event of probability `≈ 2^{−|A|}` per `c`.

### 3.2 Proposition C2R-2: why the character exponent compresses (proved)

Write `E = 2k`. Since `gcd(E, p−1) = 2`, `k` is coprime to `|H^*| = (p−1)/2`, so `y^E = π_k(y^2)` with `π_k : s ↦ s^k` a
**permutation of `H^*`**. Thus `A = π_k({j^2 : 0 ≤ j ≤ m})`: an unknown permutation of the set of small squares.
- `E ≡ (p+1)/2 (mod p−1)` has gcd 2 iff `p ≡ 3 (mod 4)`. Then `(p−1)/2` is odd, `k ≡ (p+1)/4` is the inverse of 2 mod `(p−1)/2`,
  and `π_k` is the **principal square root** on `H`. Equivalently `y^{(p+1)/2} = y^{(p−1)/2}·y = χ(y)·y`.
- Since `χ(−1) = −1`, exactly one of `±j` is a square for `j ≠ 0`, and `χ(j)j` is that one. Hence **`A = H ∩ T`**: the first layer folds
  `T` onto its own squares, *inside the interval* (217 primes `p ≡ 3 mod 4`, `p < 3000`: 0 violations).
- The second layer pairs `a` with `c − a`. For `1 ≤ c ≤ 2m`, `c − a ∈ [c−m, c+m]` does not wrap (`3m < p/2`), and `c − a ∈ T` iff `a ≥ c − m`.
  When `c − a ∈ T`, `c − a ∈ H ⇒ c − a ∈ A` automatically. So **exactness ⇔ `χ(c − a) = −1` for every square `a ∈ [−m, c−m−1]`**
  (about `c/2` boundary conditions; symmetric for `c < 0`). Also `r*(c) = #{a ∈ [c−m, m] : a < c/2, χ(a), χ(c−a) ∈ {0, 1}}`
  (217 primes, all `|c| ≤ 2m`: 0 violations; `c2r_theory.check_character_mechanism`).
- *Heuristic count.* If `χ` behaves like independent fair signs on intervals of length `≈ √p`, then `r*(c) ≈ (2m − c)/8`. So
  `d ≈ 3m/4 + 1 + c/8`, best at small `c`, where exactness costs only `O(c)` sign conditions. The threshold `3L/8 = 3m/4 + 3/8` sits at
  the mean, which is why the pass rate is about a coin flip at fixed size (14/16 here). Measured best exact `d/m` for `p ≡ 3 mod 4` from
  `2·10^4` to `3.6·10^6` (m = 70 … 942): 0.743, 0.787, 0.720, 0.731, 0.744, 0.745, 0.746, 0.728, 0.754, 0.742. That is consistent with `3/4`.
  `experimentally supported`. A proof would need `Σ_{a∈I} χ(a(c−a)) = o(|I|)` for intervals of length `≈ √p`. That is beyond Pólya–Vinogradov,
  and I have not checked whether a Burgess-type bound covers this quadratic argument.
- **The compression is the hidden order.** `(p−1) | 2E − 2`, so `gcd(a^{2E−2} − 1, N)` splits `N` (50/50 random semiprimes with both factors in `[10^5, 10^6]`;
  `c2r_theory.payoff_demo`). A constructor that supplies this `E` already factors.

### 3.3 Theorem C2R-3: squaring exponents give no constant-factor compression (proved)

(i) `E ≡ 2 (mod p−1)`, any odd prime `p`. For every `c`, exact or not, `r*(c) ≤ τ(c_0)`, where `c_0 ∈ [0, p)` represents `c`.
*Proof.* `A = {j^2}` with `j^2 ≤ m^2 ≤ p/4` as integers, so `a + b ∈ [0, p/2]` and `a + b ≡ c` forces the integer equation `j^2 + l^2 = c_0`.
The unordered solutions with `0 ≤ j < l` number at most `r_2(c_0)/4 ≤ τ(c_0)`, using `r_2(n) = 4(d_1(n) − d_3(n))`. ∎
(ii) `p ≡ 1 (mod 8)`, `E ≡ 2 + (p−1)/2`. Then `gcd(E, p−1) = 2` (with `p−1 = 8k`: `gcd(4k+2, 8k) = 2·gcd(2k+1, 4k) = 2`), and
`y^E = χ(y)y^2`. Here `A = {±j^2} ⊂ (−p/4, p/4)`, so sums lift to integers in `(−p/2, p/2)`. A pair is `j^2 + l^2 = ±c_0`
(`≤ τ`) or `j^2 − l^2 = c_0`, i.e. `(j−l)(j+l) = c_0`, which has `≤ τ(|c_0|)` solutions (`c_0 = 0` would force `j = l` with opposite characters, which is impossible).
Hence `r*(c) ≤ 2τ(|c_0|)`. ∎
**Corollary.** For these exponents every two-layer fold has `d ≥ m + 1 − 2·max_{n<p} τ(n) = L/2 − p^{O(1/log log p)}`, so there is no
constant-factor compression for large `p`, exact or not. The X1 branch `(2, 4)` is a three-layer squaring fold and falls under the same divisor-type
counting at its last layer. Check: all 779 primes `< 6000` (`E = 2`) and 186 with `p ≡ 1 mod 8` (case ii): 0 violations of the per-`c` bound, max `r*` = 4.

### 3.4 Twisted identities: a counterexample family to the unconditional `o(L)` statement (experimentally supported)

Let `(p−1) | k(E−1)` with `k ≥ 3` minimal and `gcd(E, p−1) = 2`. Then `y^E = θ(y)·y` with `θ = y^{E−1}` a character of order `k`, and
`θ(−1) = −1` (`E` even). For `ζ` in the image of `θ`, `ζ^{-1}(A ∩ ζT) ⊇ {u ∈ T : θ(u) = ζ}`. This set has density `≈ 1/k` in an interval, and
pairs inside it with `u + w = c'` give `r*(ζc') ≈ (2m − c')/(2k^2) ≈ m/k^2`. That is **linear in `L`**. Measured mean of `r*(ζc')` over all
`ζ` and `1 ≤ c' ≤ 10` (`c2r_theory.twisted_row`):

| k | m = 70 | 125 | 223 | 298 | 397 | 530 | 707 | 942 |
|---|---|---|---|---|---|---|---|---|
| 4 (mean / m/k²) | 4.88 / 4.38 | 7.97 / 7.81 | | 18.43 / 18.62 | | 33.62 / 33.12 | | 56.67 / 58.88 |
| 6 | 2.15 / 1.94 | 3.38 / 3.47 | | | | 15.08 / 14.72 | 18.63 / 19.64 | |
| 8 | 1.19 / 1.09 | | 3.12 / 3.48 | | 5.00 / 6.20 | | | 14.89 / 14.72 |

So the statement "non-character exponents have `o(L)` internal pairs" is **false without exactness**. The proof of linearity is open,
because Weil gives an error `O(√p log p)`, the same size as the main term at length `√p`. Two qualifications make this harmless for the
classification. First, **no twisted identity produced an exact fold with even one pair** (all `k ∈ {4, 6, 8}` rows, top 2 000 `c` each).
Exactness fails because `c − a` for `a` in the other rotated copies `ζ'T` lands in `H \ A` about half the time. Second, these exponents are
order-revealing as well: `(p−1) | k(E−1)`.

### 3.5 What remains open (precise)

**Conjecture C2R-X (classification, amended).** For every `ε > 0` there are `p_0, K` such that for `p > p_0`, an *exact* two-layer fold with both
gcds equal to 2 and at least `εL` internal pairs has `(p−1) | k(E−1)` for some `k ≤ K`. The audit's form, `k = 2` with `d ≤ 3L/8`, is the special case
that the data single out. By §3.4 one cannot drop exactness. By Lemma C2R-1(b) the conjecture is equivalent to a statement about translates
`c − A_E` whose character sum is non-negative while `c − A_E ∩ H ⊂ A_E`.
Status: `conjectured`. Evidence: D2 = 0 at every instance `p ≥ 4099` and in every twisted family tested. Proved cases: `E ≡ 2` and `E ≡ 2+(p−1)/2`
(Thm C2R-3, even without exactness). The character case is a proved mechanism (Prop. C2R-2). The missing tool is a bound on
`|A_E ∩ (c − A_E)|` or on `S_{A_E}(c)` for sets of size `√p`. Weil- or Burgess-type estimates do not reach that range for exponents of size `≈ p`.

## 4. Candidate card update (C-2R, after the repair cycle)

| Field | Update |
|---|---|
| **ID / mechanism** | C-2R: exact interval compression by two translated power maps (`gcd = 2, 2`) |
| **Status** | **Refuted for this model** (pre-registered rule met: 0/16 non-character passes; arm 3 = one-layer degree on all 16 large instances). The extension X1 (two-branch public mixing) is also refuted for this model (0/16). |
| **What compresses** | Only `E ≡ (p+1)/2`, `p ≡ 3 (mod 4)` (14/32 instances; `d/m ≈ 3/4` up to `m = 942`). It is the signed identity `y ↦ χ(y)y`, and it reveals `p−1` (proved). |
| **New lemma still required** | Conjecture C2R-X (§3.5): exactness + linear pairs ⇒ twisted identity. Proved for squaring exponents only. |
| **Cost in n** | Even the order-revealing fold gives only a constant factor (`3/4`) on a `√p`-degree object. No route to polynomial time. |
| **Closest prior** | Quadratic-character folding; Pollard `p−1` (the payoff identity `(p−1) \| 2E−2`). |
| **Refutation** | Done: pre-registered threshold not met by any order-independent exponent on 32/32 instances. |
| **Possible next repair** | None recommended in this family. Every route that compressed encoded a multiple of `p−1`. A different operation (not a power-map composition) would be needed. |

## 5. Reproduction

```bash
cd /home/knielsen/code/factoring-connection
.venv/bin/python -m pytest -q tests/test_repair_c2r.py                       # 6 tests (incl. audit-code cross-check, certificate)
.venv/bin/python -m src.repair.c2r_search --procs 32 --out experiments/raw/repair_C2R.jsonl        # ~10 s wall, 38 CPU-s
.venv/bin/python -m src.repair.c2r_mixing --procs 32 --out experiments/raw/repair_C2R_mixing.jsonl # ~13 min wall, 4 899 CPU-s
.venv/bin/python -m src.repair.c2r_theory --procs 16 --out experiments/raw/repair_C2R_theory.jsonl # ~5 s
.venv/bin/python -c "from src.repair.c2r_search import verify_certificate as v; print(v())"
```
Raw data: `experiments/raw/repair_C2R.jsonl`, `repair_C2R_mixing.jsonl` (+ `.log`), `repair_C2R_theory.jsonl`.
