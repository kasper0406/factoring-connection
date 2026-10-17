# Track E candidate cards (agent trackE, 2026-10-17)

Template: BRIEF §4. Status labels per AGENT_RULES §5. All formats are about the smaller factor
`p` of a balanced semiprime `N = pq` (`src/common/instances.balanced_semiprime`); `n = log2 N`.
"Decoder" = operational code in `src/trackE/{coppersmith,formats}.py` (sees `N` + the format object
only); "oracle" = harness code in `src/trackE/oracles.py` (knows `p, q`). Every card below is a
**diagnostic-oracle** result: it prices a format, it does NOT produce the format from `N`.
Analysis, proofs and tables: `notes/trackE_analysis.md`. No novelty is claimed for any decoder:
(a)/(b) are [May10] Theorems 1–3 (Coppersmith / Howgrave-Graham); the reductions are elementary.

---

### E-a — Affine-hint decoder (p ≡ r mod M, p ∈ [L, U])
- **Mechanism:** Howgrave-Graham lattice for `f(y) = y + (r M^{-1} + x_mid) mod N`, dimension `d = m + t`.
- **Target / family:** any `N` with a divisor `b ≥ L` in the hinted coset; balanced promise `L = ⌊√N⌋/2`.
- **Representation built from N:** the `d × d` basis of `N^{m-i} f^i`, `y^j f^m` (entries ≤ `m n + d log Y` bits).
- **Extraction:** LLL (fplll, (0.99,0.51)), integer roots of the short vector, `gcd(b, N)`.
- **Why it exposes a factor:** [May10] Thm 2 (Howgrave-Graham) + LLL bound; root ⇒ `b`.
- **Known lemmas:** [May10] Thms 1–3. **New lemma:** none.
- **Correctness:** **proved** whenever `2^{d(d−1)} N^{2m(m+1)} Y^{2d(d−1)} d^{2d} < L^{4md}` (exact test `rigorous_ok`); empirical thresholds 0–2 bits better than this guarantee (§1 tables).
- **Cost:** LLL in dimension `d`; to come within `c` bits of `N^{1/4}` needs `d ≈ n/(4c)`. Measured at n = 512 (LLL time at the empirical threshold): 0.003 s for k = n/4+12 (d=9), 0.13 s for n/4+4 (d=21), 2.6 s for n/4+2 (d=37), 12 s for n/4+1 (d=49), 52 s for k = n/4 (d=65).
- **Randomness:** deterministic.
- **Closest prior:** this IS the prior method ([May10] Thm 1; Coppersmith 1997 Thm 4/5 bivariate).
- **Experiment:** `exp_threshold` (n = 64…512, 3 hint shapes, 5 seeds, d ≤ 65). **Refuting behaviour:** none (it is a theorem); role here = calibrated price list.

### E-list — Candidate lists
- **Statement (proved):** a list of `T` cosets `(M_i, r_i)` costs `Σ_i decode(M_i)`; wrong candidates cost the same LLL (measured), never yield a false factor (divisibility check). A list that is a union of cosets of a coarser modulus `M' | M` with `M' ≥ N^{1/4+ε}` is decoded as one coset.
- **Audit price:** any mechanism outputting `T` candidates for `p mod M` (`M ≥ 2^{n/4+c}`) factors in `T · t_dec(n, c)`; polynomial iff `T = poly(n)`.

### E-noise — Errors in the hint
- **Bit flips in a low-bit hint:** plain decoder fails whenever any flip lies below `k` (proved: the coset is wrong; experiments 0 successes). Truncation to `M' = 2^{k'}` succeeds if no flip lies below `k'` and the rigorous condition holds at `k'` (proved). The empirical "iff, with k' ≥ k_min(d)" matched the prediction in all 360 trials (experimentally supported). Unknown flips inside the needed window: enumerate `Σ_{i≤e} C(k', i)` patterns (poly only for `e = O(1)`).
- **Additive noise on an interval hint:** absorbed by widening; costs `log2(1 + 2E/W)` bits of the `N^{1/4}` budget (experimentally exact).
- **Wrong residues in format (iii):** random-subset decoder; success per trial `C(K−e, s)/C(K, s)`; exponential at a constant error rate. CRT list-decoding (Goldreich–Ron–Sudan; Boneh STOC 2000) **not yet tested**.

### E-i — φ(N) mod 2^k  (≡ s = p+q mod 2^k)
- **Mechanism:** `x² − s x + N ≡ (x−p)(x−q) (mod 2^k)`; root set computed by a 2-adic root tree from `(N, s)`.
- **Exact ambiguity (proved, Prop. E2):** with `v = v_2(q−p) ≥ 1`: for `k > 2v` the roots are exactly the two cosets `p, q mod 2^{k−v}` (`2^{v+1}` residues); for `k ≤ 2v` one coset mod `2^{⌈k/2⌉}`.
- **Bits that suffice:** `k ≥ k_min(d) + v`; `E[v] = 2` (P(v = j) = 2^{−j} heuristically). List: 2 decodes.
- **Classification:** decodable in poly time (`n/4 + v + O(1)` bits). **Not known to be computable from N** (it is equivalent to `p + q mod 2^k`; no cheap source known).
- **Experiment:** `exp_formats sum2k`, 40 instances: success exactly when `k − v ≥ k_min`, coset modulus = `k − v` in all cases.

### E-ii — s mod M, M = ∏ ℓ_i squarefree (odd, public)
- **Exact ambiguity (proved, Prop. E3):** root set is `∏_i {p mod ℓ_i, q mod ℓ_i}` = `2^{ω'}` residues, `ω' = #{i : ℓ_i ∤ q−p}`; exactly one is `p mod M` (one is `q mod M`). Exponential when `M ≥ N^{1/4}` is built from small primes (`ω ≥ n/(4 log2 B)`).
- **Consistency information that removes it (proved):** the unordered pair `{p mod M, q mod M}` (equivalently `q − p mod M` up to one global sign) — exactly `ω' − 1` bits beyond `s mod M`. Per-prime labels suffice: e.g. the Legendre bit `(p/ℓ)` resolves `ℓ` iff `(N/ℓ) = −1` (proved).
- **Experiment:** `exp_formats sumM` (n=128): decodes to success ≈ `2^{ω'}/2` (1.6 … 38.8 for ω' = 1 … 6); with Legendre labels on `(N/ℓ)=−1` primes the list is 1 in all 39 runs.
- **Good special cases:** one prime `ℓ ≥ N^{1/4+ε}` (T = 2), or few large primes.

### E-iii — p mod ℓ for many small primes, consistent labels
- CRT ⇒ E-a with `M = ∏ ℓ`; needs `Σ log2 ℓ ≥ n/4 + c`; list size 1. E.g. `n = 512` needs the first ~30 odd primes (≤ ~130); see analysis §4.

### E-iv — one Legendre bit (p/ℓ) per prime
- **Information:** exactly 1 bit per prime (census; under the uniform-residue model). Not coset-shaped: coset list `∏(ℓ−1)/2`, i.e. `log2 M − K` bits of list for `M ≥ N^{1/4}` (n = 512: list ≈ 2^{100}). Isolating `p` information-theoretically needs ≈ `n/2` bits and an interval-CRT search (exhaustive in our toy; codex E1/E2 problem). **refuted for this decoder model** as a stand-alone format; **useful as labels** for E-ii.

### E-v — e = N mod (p−1)  (= q mod (p−1); Track B/D D3's linear complexity)
- `N − e` is a multiple of `p − 1` ⇒ Pollard exponent; `q−1 ∤ N−e` always (since `N−e ≡ p−e mod q−1` with `0 < p−e < q−1`), so a random base splits with probability ≥ 1/2. **proved**; the code's fixed bases succeeded on 80/80 instances (n = 64…512). Complete information in one number < p.

### E-vi — random known bits of p and q (density δ each, low t ≈ k_min bits) → branch-and-prune → E-a
- Width after level i: log2-width performs a walk with drift `1 − 2δ` (all branches share the known positions); typical success for `δ > 1/2` (mean-width threshold `2 − √2`). Measured n=256 (depth 69): for δ ≥ 0.5, 20/20 succeeded with final width ≤ 16 (73–105 known bits in total); δ = 0.4 and 0.45 succeeded 3/5 each; δ = 0.3 had width > 10^6. **experimentally supported**. Heavy tails for 1/2 < δ < 2 − √2. Prior: Heninger–Shacham (CRYPTO 2009) branch-and-prune; their exact theorem not re-checked by us.

### E-feed — what the other tracks' mechanisms could feed (analysis §6)
| Mechanism | Partial information it yields | Format fed | Shortfall |
|---|---|---|---|
| Track A descent along a support index `t` | `t mod 2^k`, `t ∈ {qj} ∪ {pj}` | E-a only if `t ∈ {p, q}` | path must follow `p` itself (the label problem: both children nonzero for k < log2 p, C1); level-k selection costs rank ≥ 2^k (A1) ⇒ ≥ N^{1/4} at the needed depth n/4 |
| Track C prefix/interval predicate | `p ∈ [A, A+W]` | E-b | halves the number of predicate calls (n/4 instead of n/2), each still Õ(N^{1/4}) |
| Track B/D sketch vanishing pattern (B1/B4) | event on `q mod L_α`, `L_α | p^d−1` | none (modulus depends on p) | O(1/p) bits per sketch |
| Track B/D D3 BM linear complexity | `e_p = N mod (p−1)` | E-v (complete) | obtaining it needs Θ(p) terms |
| codex E1/E2 CRT discriminant masks | `S_ℓ` = all values of `s` compatible with `N mod ℓ` | none | function of `N mod ℓ`: zero information (Lemma F1) |
| codex E3 low-bit carries | all `2^{k−1}` odd prefixes | E-a list | T = N^{1/4} at the needed depth |
| codex F1 / Pollard / Williams | multiple of `p ± 1` | E-v type (complete) | needs smoothness or luck |
