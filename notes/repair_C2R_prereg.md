# Pre-registration: C-2R two-layer interval folding (agent `repair`, 2026-10-17)

Written and committed BEFORE any of the runs below were made on the instance primes.
Source of the design: `notes/codex_r3_critic.md` §4 (GPT-6.1-sol r3 audit). Items 1–5 are
copied from the audit. Item 6 (extension X1) and item 7 (diagnostics) are mine.
They are frozen here.

**Honesty note.** This is a *replication* pre-registration, not a discovery one. The
critic already ran these searches and reported 0/32, 14/32 and 0/32 (arms 1, 2, 3). The
critic states that arm 3 was "an exploratory follow-up, not a pre-registered discovery". Before
writing this file, I ran my implementation only on non-instance primes: brute-force
validation for p ≤ 83, a cross-check against the critic's code for primes 101–997, and
timing at p = 25013 and p = 101.
The instance primes are listed in §4.

## 1. Hypothesis (audit §4, refutation threshold)

H-C2R: an exact two-layer fold whose exponents are constructible without hidden group-order
information can beat the optimal one-layer outer degree by a constant factor. Concretely,
some exponent **other than the character exponent** achieves `d ≤ 3L/8`.

Decision rule (copied): *"Reject the tested order-independent constant-factor mechanism if no
non-character exponent achieves d ≤ 3L/8 on the 16 largest specified instances."* A timeout
counts as inconclusive, never as a pass or a fail.

## 2. Model and search space (audit §4 "Exact experiment")

- `p` prime, `m = ⌊√p⌋ // 2`, `L = 2m + 1`, `Z_L = {0, −1, …, −(L−1)} mod p`, `T = Z_L + m = [−m, m]`.
- Fold: `φ(x) = ((x + m)^E + v)^2`, `Q(Y) = ∏_{w ∈ φ(Z_L)} (Y − w)`, outer degree `d = |φ(Z_L)|`.
- Exact: `{x ∈ F_p : Q(φ(x)) = 0} = Z_L`.
- Restriction (audit): first-layer gcd `gcd(E, p−1) = 2`, all even `E` in `[2, p−3]`, all
  `v ∈ F_p`. With `A = {j^E : 0 ≤ j ≤ m}`, `H = {0} ∪ (F_p^*)^2` and `c = −2v`:
  exactness ⇔ (`a ∈ A`, `c − a ∈ H` ⇒ `c − a ∈ A`), and `d = |A| − #{{a,b} ⊂ A : a ≠ b, a + b = c}`.
  Only `c` that are sums of two distinct elements of `A` can give `d < |A| = m + 1`,
  so enumerating pair sums covers every relevant translation.
- Baseline: the optimal one-layer degree `m + 1 = ⌈L/2⌉`.

## 3. Arms and threshold (copied)

Measure the minimum exact `d` per instance in three arms:
1. public exponents `2 ≤ E ≤ 64`;
2. all exponents;
3. all exponents excluding `2E − 2 ≡ 0 (mod p−1)` (the character exponent `E ≡ (p+1)/2`).

Per-instance success: `d ≤ 3L/8`, i.e. `8d ≤ 3L`. The critic reports 0/32, 14/32 and 0/32.

## 4. Instances (copied)

The first four primes at or above `2^j` and `3·2^{j−1}`, for `j = 8, 10, 12, 14`:

| j | ≥ 2^j | ≥ 3·2^{j−1} |
|---|---|---|
| 8 | 257, 263, 269, 271 | 389, 397, 401, 409 |
| 10 | 1031, 1033, 1039, 1049 | 1543, 1549, 1553, 1559 |
| 12 | 4099, 4111, 4127, 4129 | 6151, 6163, 6173, 6197 |
| 14 | 16411, 16417, 16421, 16427 | 24593, 24611, 24623, 24631 |

The "16 largest" are the `j = 12, 14` rows.

## 5. Budget and implementation

- Hard cap: 2 CPU-hours in total for everything in this file, on up to 32 cores.
- Per-instance CPU cap: 600 s for the main search and 1800 s for X1. Hitting a cap makes that instance inconclusive.
- Code: `src/repair/c2r_search.py` (my own implementation; the audit's code is used
  only for cross-checking, by extracting it from the audit note at test time),
  `src/repair/c2r_mixing.py`.
- Certificate check: `p = 1031`, `L = 33`, `φ(x) = ((x+16)^{516} − 7)^2`,
  `R = {1,4,9,16,25,36,49,64,81,196,441}`. Verify on all 1031 inputs.

## 6. Frozen extension X1: additive mixing of two public power branches

The audit's "possible repair" was: *"Permit additive mixing of branches, but state and freeze
that larger model before testing it."* My frozen interpretation follows.

- First layer: `ψ(y) = y^{E1} + β·y^{E2}`, where `E1 < E2` are both even with `2 ≤ E1 < E2 ≤ 64`
  (496 exponent pairs, all public), and `β ∈ F_p^*` (every value).
- Second layer: `φ(x) = (ψ(x+m) + v)^2`, `v ∈ F_p` (every value). `Q` and `d` are as in §2.
- Since `ψ` is even, `φ(Z_L) = {(w+v)^2 : w ∈ W}` with `W = ψ([0, m])`. With `c = −2v`,
  `d = |W| − #{{w,w'} ⊂ W : w ≠ w', w + w' = c}`. This is computed exactly for all
  `(E1, E2, β)` and the minimising `c`.
- Stage 2: for every `(E1, E2, β, c)` with `8d ≤ 3L`, test exactness by evaluating `φ` on all of `F_p`.
- Hypothesis H-X1: some public two-branch mix achieves an exact fold with `8d ≤ 3L` on at
  least one of the 16 largest instances. Decision: reject H-X1 for this model if 0/16.
  The 16 smaller instances (`p ≤ 1559`) are reported, but they do not enter the decision. There, `deg ψ ≤ 64`
  is comparable to `L` and to `p/4`, so many polynomial coincidences are expected. A pilot run at
  `p = 101` (not an instance) found exact folds with `d = 2`.
- My prediction, stated before running: 0/16 on the large instances. The `(E1,E2) = (2,4)`
  branch is a three-layer squaring fold whose pairs are sums of two squares; at p = 25013 its
  minimum `d` was 70, against a threshold of 59.6.

## 7. Pre-declared diagnostics (no pass/fail)

For every instance, over non-character exponents with gcd 2:
- D1: `max_c r*(c)` over **all** `c` (exactness not required), plus the arg-max `E` and its
  "twist order" `k = min{k ≤ 12 : (p−1) | k(E−1)}`, so that `y^E = χ(y)·y` with `χ` a
  character of order dividing `k`;
- D2: `max_c r*(c)` over **exact** `c`.
- For the character exponent: its best exact `d`, and the number of exact `c`.

Reproduction commands are in `notes/repair_C2R_results.md`.
