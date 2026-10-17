# Track E: partial-information formats and lattice decoders. Analysis (agent trackE)

Date 2026-10-17. Status labels follow AGENT_RULES §5. Cards: `notes/trackE_candidates.md`.
Generated tables (all numbers below come from them): `experiments/analysis/trackE_tables.md`.
No novelty is claimed. The decoders are textbook Coppersmith / Howgrave-Graham. The reductions in §3 are
elementary, and they are written out here so that the audit can price other tracks' mechanisms.

**Primary source checked:** [May10] A. May, *Using LLL-Reduction for Solving RSA and Factorization
Problems: A Survey*, in *The LLL Algorithm* (2010), 315–348. I read the author's 34-page preprint
(archived copy of `crypto.rub.de/.../may/paper/lll.pdf`), §§2–5. The theorem numbers below are as
printed there:
* **Thm 1:** divisor `b ≥ N^β`, monic `f` of degree `δ`: all `x0` with `f(x0) ≡ 0 (mod b)` and `|x0| ≤ c N^{β²/δ}` can be found in time `O(c δ^5 log^9 N)`.
* **Thm 2 (Howgrave-Graham):** if `g(x0) ≡ 0 (mod b^m)`, `|x0| ≤ X` and `‖g(xX)‖ < b^m/√n`, then `g(x0) = 0` over Z.
* **Thm 3** (and its proof) gives the lattice `g_{i} = N^{m−i} f^i`, `h_j = x^j f^m`, `det = N^{δm(m+1)/2} X^{n(n−1)/2}`, `m = ⌈β²/(δε)⌉`, `t = ⌊δm(1/β − 1)⌋`, and the LLL bound `‖v‖ ≤ 2^{(n−1)/4} det^{1/n}`.
* **§5.1** (high bits known): `f(x) = p̃ + x`, `β = 1/2`, `δ = 1` gives `|p − p̃| ≤ N^{1/4}`.

Coppersmith 1997 Thms 4/5 (bivariate, high/low bits; `notes/prior_art_map.md` §7) give the same thresholds.

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| E1 | Affine-hint decoder `p ≡ r (mod M)`, `p ∈ [L,U]`. It succeeds whenever the exact inequality `2^{d(d−1)} N^{2m(m+1)} Y^{2d(d−1)} d^{2d} < L^{4md}` holds (Y = half-width of the unknown multiplier range). For balanced N with `m = (d−1)/2`, this is implied by `log2 Y ≤ ν(d−1)/(4d) − 3/2 − log2(d)/(d−1)`, where `ν = log2 N` (the real number, not the bit length). | **proved** ([May10] Thms 2–3) |
| E1-emp | Measured threshold, all formats, n = 64…512, d = 3…65: `k_min = n/4 − 2 + n/(4d) + δ`, with `δ = −0.05 ± 0.72` (low bits), `+0.48 ± 0.61` (random odd M) and `+0.05 ± 0.60` (high bits, using `n/4 − 1`). The rigorous guarantee is 1–2 bits above this, and the worst-case LLL factor costs ≤ 1 bit. | **experimentally supported** (570 (n,seed,fmt,d) cells) |
| E1-cost | Decode time grows steeply as the excess `c = k − n/4` shrinks, since `d ≈ n/(4(c+2))`. At n = 512 (LLL time at the empirical threshold): 0.003 s at c = 12 (d = 9), 0.13 s at c = 4 (d = 21), 2.6 s at c = 2 (d = 37), 12 s at c = 1 (d = 49), 52 s at c = 0 (d = 65). | **experimentally supported** |
| L | A list of T cosets costs `Σ decode(M_i)`. Wrong candidates cost the same as right ones and never yield a false factor. Cosets of a coarser modulus collapse. | **proved**, measured |
| N1 | Bit flips: the plain decoder fails if any flip lies in the used window (**proved**). Truncating to `2^{k'}` succeeds if the lowest flip is at a position `≥ k'` and the rigorous condition holds at `k'` (**proved**). The empirical form "iff lowest flip ≥ k' and k' ≥ k_min(d)" is **experimentally supported**. | experiment matches the prediction in every trial |
| N2 | Additive noise E on an interval hint costs `log2(1+2E/W)` bits of the `N^{1/4}` budget. Success iff tol ≥ noise and the widened width is below threshold. | **proved** (trivial); confirmed |
| N3 | Format (iii) with e wrong residues: a random-subset decoder succeeds per trial with probability `C(K−e,s)/C(K,s)`. That is exponential in `n/log B` at a constant error rate. | **proved** formula; measured |
| E2 | **φ(N) mod 2^k ⇔ s mod 2^k.** With `v = v_2(q−p)`, the roots of `x² − sx + N (mod 2^k)` are exactly the two cosets `p, q (mod 2^{k−v})` when `k > 2v`, and one coset mod `2^{⌈k/2⌉}` when `k ≤ 2v`. So `k ≥ k_min + v` bits and 2 decodes suffice. | **proved**; brute force (300 random pairs in a scratch check + 200 in the test suite), 240 decodes |
| E2' | For `N ≡ 1 (mod 4)`, the first bit of `s` not fixed by `N` (bit 2) is equivalent to "−1 is a square mod N", i.e. "N is a sum of two squares". | **proved** |
| E3 | **s mod M, M squarefree.** The root set is `∏{p,q mod ℓ_i}` and has `2^{ω'}` elements. The unordered pair `{p,q} mod M` removes the ambiguity: given s there are `2^{max(ω'−1,0)}` possible pairs, i.e. `max(ω'−1,0)` bits under the uniform-residue model. A Legendre label `(p/ℓ)` resolves ℓ iff `(N/ℓ) = −1`. No function of `(N mod M, s mod M)` resolves it. | **proved** (bit count: under the model); measured list `2^{ω'}`, decodes ≈ `2^{ω'}/2` |
| E4 | **p mod ℓ (consistent labels):** CRT, then E1. n = 512 needs the first ~26 odd primes. Legendre bits alone give a coset list of `∏(ℓ−1)/2` (2^105 at n = 512). | **proved**; measured |
| E5 | **N mod (p−1)** (Track B/D D3's linear complexity) is complete: the Pollard exponent `N − e`, with a 2-adic descent. | **proved**; 80/80 |
| E6 | Random known bits of p **and** q at density δ (low ≈ n/4 bits each), branch-and-prune, then E1. log2-width drifts by `1 − 2δ` per level, so `δ > 1/2` (≈ n/4 known bits in total) suffices typically. Heavy tails for δ < 2 − √2. | **experimentally supported** (n = 128, 256: δ ≥ 0.5 succeeded in 40/40) |
| F1 | **Information floor.** For `gcd(m,N) = 1`, a function of `(p mod m, q mod m)` is determined by `N mod m` iff it is constant on the fiber `{(u, c/u)}` (**proved**). Under the model where p mod m and q mod m are independent uniform units, its new information is `H(g(U, c/U))` (**proved under the model**; model checked by χ²). | exact census m ≤ 64 |
| F2 | For a prime ℓ: symmetric functions are exactly the functions of `s mod ℓ`, which has `(ℓ + (N/ℓ))/2` values. The only asymmetric content is a 1-bit label (**proved**). For squarefree m with ≥ 2 odd prime factors, `{p,q} mod m` is strictly finer than `s mod m` (**proved**, E3). The same holds for 2^k with k ≥ 3 (E2) and was observed for odd prime powers 9, 25, 27 (census only). It is false for m = 2ℓ, because (Z/2)^* is trivial. | see entries |
| X | Track mechanisms (§6): with their known algorithms and within the other tracks' restricted models, none yields a format bit beyond the floor at sub-N^{1/4} cost. Track B/D's D3 quantity is complete, but the known way to obtain it costs Θ(p). Codex E1/E2 masks carry zero information. | masks: **proved** (F1); cost statements: **conjectured** / model-conditional (they rest on known-algorithm costs and the other tracks' restricted obstructions, not on lower bounds) |

**Bottom line.** For the formats tested, the decoder side is measured precisely, and it is cheap. Bivariate/erasure decoders and CRT list decoding are untested. About `n/4 + n/(4d)`
bits of a *coset-shaped* hint (a residue of p, an interval, or φ(N) mod 2^k up to the 2-adic loss v)
suffice in dimension d. That costs milliseconds to tens of seconds at n = 512. Lists of polynomial size
are affordable, but errors inside the used window are not. Every quantity that is a function of `N mod m`
carries zero bits. The smallest unit of new information in the most natural format (φ(N) mod 2^k) is
already a classical open problem: E2' shows that its first bit is deciding whether N is a sum of two squares.

---

## 1. Decoders (a)/(b): theory, implementation, threshold

### 1.1 Construction (`src/trackE/coppersmith.py`, OPERATIONAL)

Hint `H = (M, r, L, U)`: "N has a divisor `b ≡ r (mod M)` with `L ≤ b ≤ U`". Write
`b = r + M(x_mid + y)`, `|y| ≤ Y`. First `gcd(M, N)` is checked, and a nonunit M is itself returned as a
factor (rule 3). Then `c = (r M^{-1} + x_mid) mod N` and `f(y) = y + c`, so `f(y0) ≡ M^{-1} b ≡ 0 (mod b)`.
The basis rows are the coefficient vectors of `N^{m−i} f(yY)^i` (i < m) and `(yY)^j f(yY)^m` (j < d − m).
We reduce with fplll (`fpylll` 0.6.4, (δ, η) = (0.99, 0.51); about 10× faster than flint's LLL at d = 37,
with identical output norms), take integer roots of the three shortest vectors with
`fmpz_poly.roots()`, and verify `1 < gcd(r + M(x_mid + y0), N) < N`.

Here `fmpz_poly.roots()` is factorization of an integer polynomial over Z, which is the standard final step
of [May10] §3, Step 4. It never factors N.

* (a) low bits / arbitrary modulus: `M = 2^k` or any M, with `[L,U] = [⌊√N⌋/2, ⌊√N⌋]`. This is the public
  balanced promise: the smaller factor of `balanced_semiprime` lies there.
* (b) interval / high bits: `M = 1`, `[L,U] = [A, A+W]`.

### 1.2 Rigorous parameter rule (proved)

[May10] Thm 2, plus the LLL bound (valid a fortiori for the (0.99, 0.51)-reduced basis, since
`(0.99 − 0.51²)^{−1/4} < 2^{1/4}`), plus `b ≥ L`, give a sufficient condition:
`2^{(d−1)/4} det^{1/d} < L^m/√d`. With `det = N^{m(m+1)/2} Y^{d(d−1)/2}`, raising to the power `4d`
gives the exact integer test `rigorous_ok`:

    2^{d(d−1)} · N^{2m(m+1)} · Y^{2d(d−1)} · d^{2d}  <  L^{4md}.

`choose_params` returns the smallest d (with near-optimal m) that satisfies it.

**Corollary (proved).** Take logs, use `L ≥ N^{1/2}/2` (up to the negligible `−1` from `isqrt`) and `m = (d−1)/2`, and write `ν = log2 N`:

    log2 Y_max ≥ ν(d−1)/(4d) − 3/2 − log2(d)/(d−1).

The critic's exact bisection check showed that this fails by up to (d−1)/(4d) bit if ν is replaced by the bit length n. It holds in all 45 tested cases with ν.

For low bits, `Y ≈ √N/2^{k+2}`, so `k ≥ n/4 + n/(4d) − 1/2 + log2(d)/(d−1)` suffices. The worst-case LLL
factor `2^{(d−1)/4}` contributes exactly 1/2 bit, as [May10] remarks after Thm 3. As `d → ∞` this tends
to `M ≥ N^{1/4}·O(1)`, which is [May10] Thm 1 with β = 1/2, δ = 1. Coming within c bits of the
asymptote needs `d ≈ n/(4c)`. Entries have at most `≈ m·n + t·log2 Y` bits, i.e. Θ(d·n) (measured column "entry bits").

### 1.3 Measured threshold (`exp_threshold`, 5 seeds × 3 formats × d-ladder, n ∈ {64,128,256,512})

Protocol. The oracle hands the decoder `p mod 2^k`, or `p mod M` (a public random odd k-bit M), or the top
k bits of p as an interval. For fixed d (and the heuristically optimal m), k is scanned down from a
sufficient value until two consecutive failures. `k_min` is the smallest k above which every attempted
k succeeded. The full table is in `trackE_tables.md` §Threshold. Excerpt (format "low", median over seeds):

| n | d=3 | d=5 | d=9 | d=15 | d=21 | d=37 | d=49 | d=65 | n/4−2 |
|---|---|---|---|---|---|---|---|---|---|
| 64 | 18 | 17 | 16 | 15 | 15 | – | – | – | 14 |
| 128 | 39 | 36 | 34 | 32 | 32 | 31 | – | – | 30 |
| 256 | 84 | 74 | 69 | 67 | 65 | 64 | 64 | – | 62 |
| 512 | 169 | 151 | 140 | 135 | 132 | 130 | 129 | 128 | 126 |
| 512, k rigorous | 171 | 154 | 142 | 137 | 134 | 131 | 131 | 130 | |
| 512, LLL s at k_min | 4e-5 | 3e-4 | 3e-3 | 0.027 | 0.13 | 2.6 | 12 | 52 | |

**Empirical law (experimentally supported).** Fit over all 570 cells with d ≥ 3:
`k_min − (n/4 − 2) − n/(4d)` has mean −0.05 (sd 0.72) for low bits, +0.48 (0.61) for odd M, and
+0.05 (0.60) for high bits (asymptote n/4 − 1 there). So in practice the threshold is
**`n/4 + n/(4d) − 2` known bits**.

The rigorous (guaranteed) k is 1–2 bits higher. Heuristic Hermite 1.02 and exact-SVP predictions coincide
with the rigorous one to within 1 bit, so LLL's worst-case factor is irrelevant here. Arbitrary odd moduli
behave like powers of two (+0.5 bit), and high bits behave like low bits. d = 2 (m = 1) needs ≈ n/2 − 3
bits, which is essentially trivial.

**Decode cost (price list; `exp_formats cost`, 3 seeds; full table in `trackE_tables.md`).** Low-bit
format at `k = n/4 + δ`, with the empirically minimal d from the ladder {2, 3, 4, 5, 7, 9, 11, 15, 21, 29,
37, 49} and the smallest *rigorous* d (`choose_params`). Median wall time per decode:

| n | δ = 2 | δ = 4 | δ = 8 | δ = 16 | δ = 32 |
|---|---|---|---|---|---|
| 128 | d=7, 0.5 ms (rig d=15, 7 ms) | d=5, 0.2 ms | d=3, 0.06 ms | 0.05 ms | 0.01 ms |
| 256 | d=21, 68 ms (rig d=27, 0.18 s) | d=11, 4 ms (rig 15) | d=7, 0.9 ms | 0.3 ms | 0.06 ms |
| 512 | d=37, 2.3 s (rig d=51, 17 s) | d=21, 0.13 s (rig 29, 0.69 s) | d=15, 30 ms | d=7, 1.3 ms | 0.6 ms |

Wrong candidates (random residue, same k and d) cost the same time (2.37 s against 2.32 s at n = 512,
δ = 2). In 60 runs a wrong candidate never produced a factor. The largest lattice entries are 16 000 bits
(n = 512, d = 51).

## 2. List tolerance (proved + measured)

* A format that yields a candidate set S for `p mod M` is priced as a coset list. Partition S into
  cosets of the coarsest modulus `M' | M` with `M' ≥ 2^{k_min(d)}` that is constant on each part. The cost
  is `#parts × t_dec(n, d)`. Examples: format (i) (2 cosets mod `2^{k−v}`, not `2^{v+1}` residues mod
  `2^k`), and erasures in high positions (truncate).
* The decoder cannot reject a wrong candidate early. LLL and root finding run to completion, and the
  measured wrong-candidate times equal the right-candidate times. A wrong candidate never yields a
  "factor", because every output is verified by divisibility. The decoder therefore cannot produce false
  positives, and the ACCIDENTAL-hit rate is 0. The only possible STRUCTURAL output is p itself, or q if the
  promise interval is widened to contain it.
* **Audit rule.** A mechanism producing T candidates for `p mod M` with `log2 M ≥ n/4 + c` factors in
  `T · t_dec(n, n/(4c))`. It is polynomial iff `T = poly(n)`. Equivalently, only `O(log n)` of the
  `log2 M` hint bits may be "missing".

## 3. Noise tolerance

### 3.1 Wrong low bits (proved; `exp_formats flips`)

If the hint residue `r' ≢ p (mod M)`, no `b ∈ [L,U]` with `b ≡ r' (mod M)` divides N, except possibly q.
q lies outside the balanced promise, so the decoder returns None. A flip at bit j leaves `r' ≡ p` correct
mod `2^j`. Hence the *truncation decoder* (use only the low `k'` bits) succeeds iff `lowest flip ≥ k'` and
`k' ≥ k_min(d)`.

Results (n = 128, d = 11 and n = 256, d = 15; 5 seeds × extra ∈ {0, 16, 64} × e ∈ {1,2,4} × 4 trials;
single-flip sweep at n = 256): the plain decoder succeeded 0 times whenever a flip was present. The
truncation decoder matched the prediction in every trial (`trackE_tables.md` §Bit-flip).

**Which errors are tolerated:**
* errors confined to positions ≥ `k_min(d)`: tolerated at zero extra cost (one decode at `2^{k_min}`);
* e unknown errors inside the needed window: enumerate `Σ_{i≤e} C(k', i)` patterns, which is
  polynomial only for e = O(1);
* known erasure positions: ×2 per erased bit. A contiguous erased block could use bivariate Herrmann–May
  "linear equations modulo divisors" (Asiacrypt 2008; their bound for two unknown blocks is
  `XY ≤ N^{0.207}`, against `N^{0.25}` for one block, which I derived from their general formula but did
  not re-verify in the paper): **not yet tested**;
* i.i.d. flips at rate ρ: `P(clean window) = (1−ρ)^{n/4}`, which is `1/poly` only for `ρ = O(log n / n)`.

A "smaller M" variant therefore helps exactly when the error support avoids a sub-modulus
`M' ≥ 2^{k_min}`. For residues mod many primes (format iii) this generalizes: any error pattern supported
outside some `M' | M` with `M' ≥ 2^{k_min}` is tolerated by decoding at M'. If the support is known this
is free; otherwise it is a search over sub-moduli (§3.3).

### 3.2 Additive noise on intervals (proved; `exp_formats highnoise`, n = 256, d = 15)

If the hint interval is off by E (unknown sign), widen by `tol ≥ E`. The decoding width becomes
`W + 2tol`, costing `log2(1 + 2tol/W)` bits. Observed: 100% success when tol ≥ |noise| and the width is
below the d = 15 threshold (log2 width ≤ 59–60 against n/4 − n/(4d) ≈ 59.7). 0% success whenever
tol < |noise|.

### 3.3 Wrong residues in format (iii) (`exp_formats crtnoise`, n = 128, d = 11)

The random-subset decoder succeeds per trial with probability `C(K−e, s)/C(K, s)`, where s is the number
of residues needed for `need` bits. Measured mean decodes (K = 40, e = 16): 35 (predicted ≈ 29). At a
constant error fraction and small primes this is exponential in `s ≈ n/(4 log2 B)`. Better decoders exist
in the CRT-code literature (Goldreich–Ron–Sudan 2000; Boneh, STOC 2000), but I did not test them.

## 4. Derived formats: exact reductions

### 4.1 Format (i): s = p + q mod 2^k, equivalently φ(N) = N + 1 − s mod 2^k

**Proposition E2 (proved).** Let `p ≠ q` be odd, `v = v_2(q − p) ≥ 1` and `s ≡ p + q (mod 2^k)`. Then
`x² − sx + N ≡ (x − p)(x − q) (mod 2^k)`, so x is a root iff `v_2(x−p) + v_2(x−q) ≥ k`. Write
`a = v_2(x − p)`, `b = v_2(x − q)`.
* If `a < v` then `b = a`. Hence `a ≥ ⌈k/2⌉`, which requires `v > ⌈k/2⌉`.
* If `a > v` then `b = v` and we need `a ≥ k − v`.
* If `a = v` then `b > v` (since `x − q = (x−p) + (p−q)` with both terms of valuation v), and we need
  `b ≥ k − v`.

Hence, for `k > 2v`, the root set is `{x ≡ p (mod 2^{k−v})} ∪ {x ≡ q (mod 2^{k−v})}`: two disjoint cosets,
`2^{v+1}` residues. For `k ≤ 2v` it is the single coset `{x ≡ p (mod 2^{⌈k/2⌉})}`, which contains q.
∎

Checked by brute force on 200 random odd pairs in the test suite (plus 300 in a preliminary scratch check) (`test_quad_root_cosets_bruteforce_and_v_formula`).

* **Decoder** (`formats.decode_sum_mod_2k`, operational). A 2-adic root tree for `x² − sx + N` computes
  the cosets from `(N, s, k)` alone. The tree is O(k) nodes per coset; it also reveals v. E1 is then run
  on each coset (≤ 2 decodes).
* **Bits needed:** `k ≥ k_min(d) + v`. Since `P(v = j) = 2^{−j}` for random odd p, q (heuristic,
  `E[v] = 2`), the cost over format (a) is two decodes and about 2 bits.
* **Experiment** (`sum2k`, n = 128/256, 20 seeds each): the coset modulus equals `2^{k−v}` in every case.
  Success occurs exactly when `k − v ≥ k_min`. At k = kmin_guess, the instances with v ≥ 3 failed and the
  instances with v ≤ 2 succeeded.
* **Prime powers ℓ^k, ℓ odd:** the same factorization `(x−p)(x−q)` gives exactly 2 roots when ℓ ∤ q − p.

**Proposition E2' (proved): the first bit of s beyond the floor.** Let `N ≡ 1 (mod 4)`. Then `p ≡ q (mod 4)`.
Write `p = ε + 4a`, `q = ε + 4b` with `ε ∈ {1, 3}`. Then `N ≡ 1 + 4(a + b) (mod 8)`, so `a + b (mod 2)` is
fixed by `N mod 8`, and `s = 2ε + 4(a+b)`. Therefore `s mod 8` is equivalent to `ε = p mod 4`, which is
equivalent to "−1 is a QR mod p and mod q", i.e. "−1 is a square mod N", i.e. "N is a sum of two squares".
(Brute-force check: given `N mod 2^K`, `s mod 4` takes 1 value and `s mod 8` takes 2 values.)

For `N ≡ 3 (mod 4)` the first free bit is bit 3 (also checked by brute force). So **φ(N) mod 8 for
N ≡ 1 (mod 4) is exactly the "sum of two squares" decision**, for which I know no polynomial-time method
that avoids factoring.

### 4.2 Format (ii): s mod M, M = ∏ ℓ_i squarefree, ℓ_i odd public primes not dividing N

**Proposition E3 (proved).** Modulo each ℓ, `x² − sx + N ≡ (x − p)(x − q)`, so the root set mod ℓ is
`{p mod ℓ, q mod ℓ}`. It has 2 elements iff `ℓ ∤ q − p`. By CRT the root set mod M has `2^{ω'}`
elements, `ω' = #{i : ℓ_i ∤ q − p}`. Exactly one element is `≡ p (mod M)` and one is `≡ q (mod M)`; the
rest are mixed. A mixed residue r has no divisor of N in `[L, U]` congruent to it (the only candidates are
p and q), so its decode returns None. Hence the list is `2^{ω'}` cosets with 1 success, and the expected
number of decodes in random order is `(2^{ω'}+1)/2`.

* **Exponential unless ω' = O(log n).** `M ≥ 2^{k_min}` from primes ≤ B forces `ω ≥ (n/4)/log2 B`. At
  n = 512, the first odd primes need ω ≈ 26, i.e. a 2^26-element list. One prime `ℓ ≥ 2^{k_min}`
  gives T = 2.
* **Which consistency information removes it (proved).** The `2^{ω'}` candidates correspond to the
  choices of root per ℓ. The list collapses to the 2 "pure" candidates iff the decoder knows, for every i,
  which root at ℓ_i belongs to the same prime as a fixed root at ℓ_1. That is exactly the unordered pair
  `{p mod M, q mod M}`, i.e. `ω' − 1` bits beyond `s mod M`. Equivalent forms:
  * `q − p mod M` up to one global sign (then `p ≡ (s ∓ t)/2`);
  * per-ℓ labels `g_ℓ(p mod ℓ)` that separate the two roots.
* **Lemma (proved).** If `uv ≡ c (mod ℓ)` and `u ≢ v`, then `(u/ℓ) ≠ (v/ℓ)` iff `(c/ℓ) = −1`, because
  `(u/ℓ)(v/ℓ) = (c/ℓ)`. So the Legendre bit `(p/ℓ)` is a perfect label exactly on the primes with
  `(N/ℓ) = −1`, about half of all primes. More generally, a character χ of order e | ℓ − 1 separates the
  roots iff `χ(u)² ≠ χ(c)`. Exhaustive test for ℓ ≤ 23 passes.
* **No label from (N mod M, s mod M) (proved).** For a prime ℓ_i | M with `ℓ_i ∤ q − p`, swap only the
  ℓ_i-component: replace `(p mod ℓ_i, q mod ℓ_i)` by `(q mod ℓ_i, p mod ℓ_i)` and keep the other
  components. This changes the labeled pair `(p mod M, q mod M)` but preserves both `N mod M` and
  `s mod M`. Every such per-prime swap is realised in the fiber, so no function of `(N mod M, s mod M)`
  can choose the root at ℓ_i.
* **Experiment** (`sumM`, n = 128, d = 7, 5 seeds per row). The list size was exactly `2^{ω'}` for
  ω' = 1…6, and the mean decodes to success were 1.6, 1.8, 4.4, 10.6/9.4, 14.6, 31.2/38.8, tracking
  `≈ 2^{ω'}/2`. With Legendre labels on primes with `(N/ℓ) = −1`, the list was 1 in all 39 runs.

### 4.3 Format (iii): p mod ℓ for many small primes, consistent labels

CRT gives `p mod ∏ℓ`, and then E1 decodes with T = 1. The smallest odd primes needed for
`n/4 + 4` bits are:

| n | 64 | 128 | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|---|---|
| #primes | 7 | 10 | 16 | 26 | 44 | 75 |
| largest prime | 19 | 31 | 59 | 103 | 197 | 383 |

(**proved** reduction; decodes verified in tests and in `crtnoise` with e = 0.)

### 4.4 Format (iv): one Legendre bit per prime

* **Information:** exactly 1 bit per prime (census: `H((p/ℓ) | N mod ℓ) = 1`). It is asymmetric.
* **Not coset-shaped (proved).** The set of residues mod `M = ∏ℓ` with prescribed symbols has
  `∏(ℓ−1)/2` elements. Feeding E1 needs `M ≥ 2^{k_min}`, so the list has `≈ log2 M − K` bits: 25, 52 and
  105 bits for n = 128, 256 and 512 (table in `trackE_tables.md`).
* **Information-theoretic decoding** needs ≈ `log2(U − L) ≈ n/2 − 1` bits. Toy n = 36: survivors in the
  promise interval roughly halve per prime; the true p always survives. Finding the survivors is the
  interval-CRT search problem of codex E1/E2 (Hittmeir's sieve), with no polynomial algorithm known.
  **refuted for this decoder model** as a stand-alone format.
* **Useful as labels** for (ii) on primes with `(N/ℓ) = −1` (Lemma above). Higher power-residue symbols
  of order e give `log2 e` bits and a list of size `(ℓ−1)/e` per prime.

### 4.5 Format (v): e = N mod (p − 1) = q mod (p − 1)

`E = N − e ≡ 0 (mod p − 1)`, so `a^E ≡ 1 (mod p)`.
* `E ≡ p − e (mod q − 1)` with `0 < p − e < q − 1`, so **`q − 1 ∤ E` always** (critic's observation).
  Hence `{a : a^E ≡ 1 mod q}` is a proper subgroup. A random base a gives `a^E ≢ 1 (mod q)` with
  probability ≥ 1/2, and then `gcd(a^E − 1, N) = p`.
* The code uses the fixed bases 2, 3, 5, 7, 11, 13 (plus the 2-part descent), so for the code the
  ≥ 1/2-per-base bound is heuristic. It never failed.

**proved** (standard); 80/80 instances at n = 64…512 (`pm1`). This is exactly the quantity that Track B/D's
D3 shows is the Berlekamp–Massey linear complexity mod p.

### 4.6 Format (vi): random known bits of p and q

`hs_branch` keeps all `(a, b) mod 2^i` with `ab ≡ N` that are consistent with the known bits. All branches
share the same known positions, so per level:
* both bits unknown: every branch doubles;
* one bit known: width unchanged;
* both bits known: each wrong branch survives with probability 1/2.

So log2-width is a walk with drift `(1−δ)² − δ² = 1 − 2δ`. The mean width grows like `(2 − 2δ + δ²/2)^i`,
which puts the mean threshold at `δ = 2 − √2`. Measured (5 seeds; `hsbits`; caps: width 2^20, 4096 decodes):
* n = 256, depth 69, δ ≥ 0.5: 20/20 succeeded, final width ≤ 16, ≤ 7 decodes, 73–105 known bits in
  total;
* δ = 0.45 and 0.4: 3/5 each; the other runs hit the decode or width caps (widths up to 1.1·10^5 and
  1.3·10^6);
* δ = 0.3: width > 10^6 in 5/5.

At n = 128 (depth 36) the transition is the same: δ ≥ 0.5 gives ≤ 21 decodes on average, δ = 0.3 gives
2/5.

So ≈ n/4 random bits spread over the low n/4 positions of p **and** q work about as well as n/4
contiguous low bits of p in the *typical* case. **experimentally supported.** For `1/2 < δ < 2 − √2`
the median width stays small but the *expected* width grows: heavy tails, rare blow-ups. Prior art: Heninger–Shacham, CRYPTO 2009
(branch-and-prune; their exact statement not re-checked). Codex E3 is the δ = 0 case (`2^{k−1}` states,
N^{1/4} at the needed depth).

## 5. The information floor

**Lemma F1 (proved).** Fix m with `gcd(m, N) = 1` and let `c = N mod m`. The pairs `(p mod m, q mod m)`
compatible with c form `F_c = {(u, c u^{−1}) : u ∈ (Z/m)^*}`, of size φ(m).
* A function g of `(p mod m, q mod m)` is a function of `N mod m` iff it is constant on every `F_c`. This
  is immediate: `N mod m = uv`.
* If `p mod m` and `q mod m` are independent and uniform on units (Dirichlet; empirically: χ²/dof ≤ 1.83
  over 34 cells from 4000 random 64-bit balanced semiprimes, mean 0.91), then conditional on c, u is uniform
  on `(Z/m)^*`. So the information g carries beyond N is exactly `H(g(U, c/U))`.

**Lemma F2 (proved).** Let `m = ℓ` be an odd prime. The unordered pair `{u, v}` is the root set of
`x² − sx + c`, so the swap-invariant functions are exactly the functions of `(c, s)`. The number of
values of s on `F_c` is `(ℓ + (c/ℓ))/2` (tested exhaustively for ℓ ≤ 17). So
`H(p mod ℓ) = log2(ℓ−1) = H(s) + H(label | s)` with `H(label | s) ≤ 1`.

For squarefree m with ≥ 2 odd prime factors (E3), and for 2^k with k ≥ 3 (E2), `{u, v} mod m` is strictly finer than `s mod m`. The census also shows this at m = 9, 12, 16, 25, 27; the odd prime powers are not proved. It is false for m = 2ℓ (census: "f(s)" holds at m = 6, 10, 14, 22, since (Z/2)^* is trivial). That gap is the label information of format (ii).

**Census** (`src/trackE/info_floor.py`, exact enumeration for m ∈ {3,…,47} ∪ {8,16,32,64,12,20,28,9,25,27};
`trackE_tables.md` §Information floor). Summary, prime ℓ:

| quantity (per ℓ) | new bits given N mod ℓ | swap-symmetric | carries new info? |
|---|---|---|---|
| N mod ℓ, (N/ℓ), any function of N mod ℓ | 0 | yes | **no** |
| Jacobi (a/N) | 0 (a function of N mod 4a; census m = 12, 20, 28) | yes | **no** |
| codex E1 mask `S_ℓ = {s : s² − 4N square}` | 0 (it is the set of values of s on F_c) | yes | **no** |
| s = p+q mod ℓ, φ(N) mod ℓ, {p,q} mod ℓ, p² + q² | ≈ log2 ℓ − 1 (2.67 of 3.58 at ℓ = 13) | yes | yes, needs labels across ℓ |
| (p − q)² mod ℓ | ≈ log2 ℓ − 2 | yes | yes |
| [ℓ \| φ(N)], [p ≡ q mod ℓ], lcm(ord p, ord q) | ≤ 1 (event of probability ≈ 2/ℓ, 1/ℓ) | yes | tiny |
| (p/ℓ) | 1 (a label iff (N/ℓ) = −1; symmetric iff (N/ℓ) = +1) | no | yes |
| [ℓ \| p − 1], ord_ℓ(p), e-th power class | `H(1/(ℓ−1))`, …, log2 e | no | yes |
| p mod ℓ | log2(ℓ − 1) | no | yes (full) |

Prime powers: at `m = 16` (N odd), `(p − q)²` and `p² + q²` are determined by `N mod 16`, s carries
exactly 1 bit, and the unordered pair carries 2–2.5 bits. This is consistent with E2/E2'.

**Reading the table.** A per-modulus mechanism only helps if its output depends on `(p mod m, q mod m)`
non-trivially on the fiber. Symmetric outputs (s, φ, orders of symmetric groups) need a consistent label
across moduli: `2^{ω'}`, or one large modulus. Asymmetric outputs need a coset shape (p mod ℓ, not a
character value). All the decoders together need ≈ `n/4` such bits.

## 6. Connection to the other tracks

For each mechanism: what partial information it produces, which format it could feed, and the shortfall.

* **Track A (binary sections).** A descent path reveals `t mod 2^k` for a support index
  `t ∈ S_q ∪ S_p = {qj} ∪ {pj}`. That feeds E1 only if `t = p` (or q), when it is format (a) with
  `k ≈ n/4 + n/(4d) − 2`: half the depth of a full descent.
  * *Shortfall 1 (C1):* at levels `k < log2(p−1)` both children are nonzero, so zero tests carry no
    information about which child contains p. "Contains index p" is itself the label we lack. For
    `t = qj` with unknown j, `t mod 2^k` says nothing about q mod 2^k (j ranges over all odd residues).
  * *Shortfall 2 (A1):* the level-k selection costs rank ≥ 2^k. At the depth E1 needs, that is
    `≥ 2^{n/4} = N^{1/4}`, no better than Pollard–Strassen/BGS.

  Pure classes (C3) give a factor directly.
* **Track C (prefix/interval unit predicate).** `U(K) = [K ≤ p]` by binary search gives format (b). E1
  needs width `≤ 2^{n/4 + n/(4d) − 1}` (smaller with d), so only ≈ n/4 adaptive calls are needed instead
  of ≈ n/2. Each call is still `Õ(√K) ≈ N^{1/4}` (BGS), so the exponent does not change.
  * Enumerating `N^{1/4}` windows of width `N^{1/4}` with one E1 decode each is the classical
    `N^{1/4+o(1)}` high-bits algorithm. Harvey's `N^{1/5}` already beats it.
  * Relaxation worth noting: E1 only needs a *coarse* predicate, correct whenever `|K − p| > N^{1/4−ε}`.
    Track C's obstructions (O1–O4) are stated for exact interval zero sets. Whether coarse predicates are
    cheaper in their models is **not yet tested**. I expect not: the first binary-search queries span
    `≈ √N`.
* **Track B/D (sketches, BM).**
  * B1/B4: one sketch's vanishing pattern depends on q only through `q mod L_α` with `L_α | p^d − 1`. The
    modulus depends on p, so this is not a residue of p modulo a public M and feeds none of (a)–(iv).
    A non-split observation carries `O(1/p)` bits (B6); a split is already a factor.
  * D3: the BM linear complexity mod p is `e_p = q mod (p−1) = N mod (p−1)`. That is **format (v), which
    is complete**: one exponentiation factors N (E5, 80/80). The shortfall is entirely in obtaining e_p
    (Θ(p) sequence terms). Note what this makes precise: Track B/D does not need the BM *run*, only its
    *output length*. Any `o(p)` method to learn the linear complexity of `c ↦ F_N(c) mod p`, or any
    multiple of `p−1`, suffices. This is the Pollard `p−1` interface again.
* **codex r1 E1 / r2 E2 (CRT discriminant masks).** `S_ℓ` is the image of `F_{N mod ℓ}` under `s`, a
  function of `N mod ℓ`, so by F1 it carries **zero** bits about the true s. Format (ii) needs the actual
  `s mod ℓ` (≈ log2 ℓ − 1 bits each, ≥ n/4 in total) *and* labels. Using masks alone requires isolating s
  in an interval of length `≈ √N`: ≈ n/2 masks plus the exponential interval-CRT search, which is exactly
  E1/E2's open obligation. E2's conversion `p'(s) = −p/(q−p)` additionally loses `log2(p/(q−p))` bits.
* **codex r2 E3 (low-bit carries).** This is format (vi) with δ = 0: `2^{k−1}` prefixes, and at
  `k ≈ n/4` that is `N^{1/4}` decodes. Any source of ≈ n/4 random bits of p and q (density > 1/2 in the
  low positions) would make it polynomial (§4.6).
* **codex r1 F1 / Williams / Pollard.** Order coincidences. A known multiple of `p ± 1` (or of
  `ord_p(a)`) is format-(v)-type and complete. F1's `d = gcd(p−ε_p, q−ε_q)` is a splitting probability,
  not a decodable hint.
* **codex r2 H2 (weighted sums).** The integer interval of `S = ap + bq` is format-(b)-like. Its masks
  are floor information (F1).

**Net.** The mask statement is **proved** (F1). The cost statements are **conjectured / model-conditional**: they use the costs of the known algorithms and the other tracks' restricted-model obstructions, not lower bounds. With that caveat, no mechanism in Tracks A–D or the codex rounds yields any bit of a
coset-shaped hint beyond the N-floor at cost below `N^{1/4}`. The two complete hints that appear (e_p,
pure classes) are already factors, at Θ(p)-type cost.

## 7. Most decisive open question

Every format that E1 can consume needs ≈ n/4 bits that are *not* functions of `N mod m` (F1). The
cheapest-looking one is φ(N) mod 2^k (E2: n/4 + v bits, two decodes). Its first free bit is already
"is N a sum of two squares?" (E2'). Concrete question:

> **Is there any polynomial-time R_N-algebra observable whose output distribution depends on
> `s mod 8 = (p + q) mod 8` beyond `N mod 8`, i.e. which decides whether −1 is a square mod N
> (N ≡ 1 mod 4) with advantage `1/poly(n)`?**

A positive answer is necessary for format (i) (bits beyond the floor are needed at every level). It is
also a clean single-bit target that the sketches of Tracks A/B could be tested against: does any sketch
statistic correlate with `[p ≡ 1 (mod 4)]` over instances with `N ≡ 1 (mod 4)`?

## 7a. Review

An independent read-only critic (Opus subagent, 2026-10-17) checked §§1–6 against the code with small
exact computations. Verdicts: the rigorous condition, the determinant, E2, E2', the Legendre lemma, the
counts in F2, the drift arithmetic and the leakage separation are CORRECT.

Fixes applied after the review:
* the corollary of §1.2 now uses ν = log2 N, not the bit length;
* the composite-m claim in F2 is restricted (it is false for m = 2ℓ);
* the per-prime swap proof replaces the global-swap argument in E3;
* the Pollard second case was vacuous and is removed;
* the heavy tails in §4.6 are noted;
* the model-dependence and model-conditionality labels in §0 X and §6 are added;
* the cost driver's rigorous mode now builds its hint through `oracles.py`.

The static leakage test is a substring audit, so it is weak. The operational modules were also checked by
reading them.

## 8. Reproduction

All commands are run from the project root, with `.venv` (python-flint 0.9.0; added `fpylll 0.6.4`,
`cysignals 1.13.1` for LLL, `pypdf` for reading [May10]).

```
.venv/bin/python -m pytest -q tests/test_trackE.py                            # 20 tests
.venv/bin/python -m src.trackE.exp_threshold experiments/configs/trackE_threshold.json
                                   # -> experiments/raw/trackE_threshold.jsonl (≈ 25 min, 16 workers)
TRACKE_WORKERS=8 .venv/bin/python -m src.trackE.exp_formats \
    pm1 sum2k legendre highnoise sumM crtnoise hsbits flips cost      # -> experiments/raw/trackE_<name>.jsonl
.venv/bin/python -m src.trackE.info_floor                             # -> experiments/raw/trackE_info_floor.jsonl
.venv/bin/python -m src.trackE.analyze                                # -> experiments/analysis/trackE_tables.md
```

Instances come from `balanced_semiprime(n, seed)`: seeds 0–4 for the threshold, 0–19 for sum2k/pm1, and
10000+ for the uniformity check. Auxiliary moduli are public: random odd M seeded by `(n, k)`, small
primes, and random primes seeded by `(n, seed, lbits)`.

Code map:

| path | role |
|---|---|
| `src/trackE/coppersmith.py` | OPERATIONAL decoder E1 |
| `src/trackE/formats.py` | OPERATIONAL reductions (i)–(vi) |
| `src/trackE/oracles.py` | HARNESS oracles |
| `src/trackE/info_floor.py` | census |
| `src/trackE/exp_*.py`, `src/trackE/analyze.py` | drivers |
| `tests/test_trackE.py` | tests, including a static leakage audit of the operational modules |
