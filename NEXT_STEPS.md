# Next steps

Prioritized. Each item has a refutation criterion fixed in advance and an estimated cost. Context: `REPORT.md` §1, §5.
Freeze each item's parameters in a pre-registration note before running on its instances, as was done for C-2R.
Revised 2026-10-17 after the GPT-6.1-sol report review (`codex/r5_report_review.out`, items 1, 3, 7, 9, 11).

Throughout, `D_n` is the explicit challenge distribution: `N = pq` from `balanced_semiprime(n, seed)` (`src/common/instances.py`,
top-two-bit primes of `n/2` bits, `n` even, uniform seeds), conditioned on `N ≡ 1 (mod 4)`. On `D_n`, `ε = (−1/p) = (−1/q)`
is defined and both values occur. A predictor's advantage is measured against the best predictor that sees only
`N mod 2^k` (the "information floor" of `notes/trackE_analysis.md` §5, which is relative to `N mod m`, not to `N` itself).

## 1. Correlation battery for `ε` and `QR_2` with a positive control and a power analysis (cheap experiment)

**Question.** Does any statistic computable from `N` in polynomial time by the project's code predict `ε` on `D_n`, or
`QR_2 = [2 is a square mod N]` on `D_n` conditioned on `N ≡ 1 (mod 8)`, better than the floor?

**Menu (frozen before running).** Track A section values at fixed public `(M, c)`; Track B sketch coefficients, ranks and
residues mod small primes; Berlekamp–Massey lengths over a fixed prefix; integer-representative statistics: Jacobi symbols of
sketch values, continued-fraction **prefix** features of `√N` (first `T` partial quotients, `T ∈ {64, 1024}`), and counts of
reduced forms of `Cl(−4N)` with `a ≤ 50`. Complete continued-fraction periods are *not* polynomial-time in general and are
not in the menu as a candidate statistic.

**Positive control (frozen).** "The continued-fraction expansion of `√N` closes its period within `T` steps and the period
is odd." This event implies that `x² − N y² = −1` is solvable, hence that every odd prime factor is `≡ 1 (mod 4)`, hence
`ε = +1`. It is a known one-sided signal, not leakage. The battery must detect it, at the rate predicted by its measured
frequency, on the same instances. If it does not, the run is invalid.

**Power analysis (frozen).** With `10^5` instances per size and `S = 100` statistics, a two-sided `5σ` level
(`p ≈ 5.7·10^{−7}`) Bonferroni-corrected over `S` needs `|z| ≥ 6.0`. At 80% power the minimum detectable correlation is
`|ρ| ≈ (6.0 + 0.84)/√(10^5) ≈ 0.022`. Report this minimum detectable effect next to every null result.

**Run.** `10^5` instances of `D_n` at `n = 32, 48, 64`, seeds disjoint from all earlier runs. Leakage audit of every hit;
any hit re-run on fresh seeds.

**Refutation criterion.** If no menu statistic other than the positive control exceeds the corrected threshold at all three
sizes, the claim "some statistic of this menu predicts `ε` (or `QR_2`) with `|ρ| ≥` the detectable effect" is refuted for this
menu and these sizes. This excludes only effects of that size; it says nothing about smaller correlations or other statistics.

**Cost.** A few CPU-hours (CF prefixes and small-`a` form counts are `poly(n)` per instance); 1 agent-day including the
pre-registration and cost benchmarks.

## 2. Beyond `s mod 64`: what determines `s mod 128`? (cheap experiment)

**Status of the earlier version.** The planned mod-64 experiment was vacuous (r5 review item 3). For `N ≡ 1 (mod 8)`,
`uv mod 4 = b(1 + 2[N ≡ 1 mod 16])` with `b = [2 is a non-square mod N]`, so `(N mod 64, ε, QR_2)` already determines
`s mod 64`. Verified for all 3 418 pairs `p < q < 1000` with `N ≡ 1 (mod 8)` (`src/repair/verify_r5.py`, 0 failures).
`Cl(−N)` is dropped: `−N ≡ 3 (mod 4)` is not a discriminant.

**First result (2026-10-17, same script).** `s mod 128` (which needs `uv mod 8`) is **not** a function of
`(N mod 128, ε, QR_2)`: 32 of 64 keys ambiguous. It is also **not** a function after adding `r4` and `r8` of both `Cl(−4N)` and
`Cl(−8N)`: 73 of 248 keys ambiguous, covering 1 199 of 3 418 pairs. So 2-power ranks of these two class groups do not continue
the tower.

**Question.** Is `s mod 128` a function of `(N mod 128, ε, QR_2)` together with 2-power ranks of `Cl(−4kN)` for small public
`k` (`k ∈ {1, 2, 3, 5, 6, 7}`, all forms primitive)? Equivalently: does any public-discriminant class-group invariant
determine the individual value of an order-8 Dirichlet character mod 32 at `p` (only its product over `p, q` is public)?

**Attack.** Extend `src/repair/verify_r5.py` with the extra discriminants and `r16`; tabulate as before.

**Refutation criterion.** Two pairs with identical extended keys and different `s mod 128` refute determination by that key
set. If no ambiguity remains, state the functional form and attempt a proof via Rédei symbols (Stevenhagen arXiv:1806.06250).

**Cost.** About 5 CPU-minutes per added discriminant at bound 1000. Half an agent-day.

**Caveat.** Even a positive answer is not an algorithm: no way is known to compute `r4` or `r8` from an unfactored
discriminant faster than computing the class group.

## 3. Standard-model reduction question for `ε` (theory; long-open, low probability of resolution)

**Question (well-posed version).** Is there a randomized polynomial-time reduction that factors `N ← D_n` with probability
`1/poly(n)`, given oracle access to an algorithm `A` with `Pr_{N←D_n}[A(N) = ε(N)] ≥ 1/2 + δ`? The oracle may be called on
arbitrary integers, but its guarantee holds only on `D_n`. (The earlier "generic ring" phrasing was vacuous. With `N` public
and no hidden input, ring operations, equality and gcd tests can be simulated on canonical representatives, so a "generic"
decider is an ordinary one. GPT-6.1-sol r5 item 1.)

**Attack.** Literature first: Adleman–McCurley 1994 (O11), Jager–Schwenk (whose re-randomization over the input element is
exactly what is missing here), and any worst-case-to-average-case results for fixed-element residuosity. Then try the obvious
self-reductions (multipliers `kN`, related moduli) and record exactly where each fails.

**Refutation criteria (two distinct things).** (a) A specific extraction lemma ("from `A`, build a splitter with
probability `≥ δ/poly`") is refuted by an explicit oracle `A` on `D_n` with advantage `δ` for which the lemma's construction
provably fails. (b) The computational *separation* (`ε`-prediction easy, factoring hard) cannot be refuted or established
by such examples: it would need a hardness assumption or an unconditional result. An algorithm that predicts `ε` but never
outputs a factor does not refute equivalence, because a reduction may still use it.

**Cost.** 1–3 agent-days plus 2 GPT-6.1-sol rounds. Expected outcome: an annotated map of why known reductions fail.

## 4. Coarse prefix predicates in the Track C models (lemma)

**Question.** Coppersmith only needs a *coarse* predicate, correct whenever `|K − p| > N^{1/4−ε}`
(`notes/trackE_analysis.md` §6). Track C's obstructions O1–O2′ are stated for exact interval zero sets. Do they survive for
coarse predicates?

**Attack.** Restate O1 (fixed-dimension powers) and O2 (one power-map layer) with a "don't-care" window of width `W` around
`p`. Check exhaustively for small `p` (`p < 200`, all `W ≤ √p`) using the existing Berlekamp–Massey and fold-search code.

**Refutation criterion.** A construction in the O1 model with order `o(p)`, or in the O2 model with `deg Q < ⌈(L − W)/2⌉`,
refutes the "obstructions transfer" hypothesis. A proof with the window removed from the bound confirms it.

**Cost.** 1 agent-day of mathematics; exhaustive checks under 1 CPU-hour.

## 5. Conjecture C2R-X at larger p and the missing additive bound (lemma + experiment)

**Question.** Exact two-layer folds with `≥ εL` internal pairs come only from exponents with `(p−1) | k(E−1)`, `k ≤ K`
(`notes/repair_C2R_results.md` §3.5). Note that non-character exact compression does occur at small `p`: `p = 257, m = 8,
E = 70, v = 186` is exact with `d = 7 < 9` (`src/repair/verify_r5.py`), though it does not reach `d ≤ 3L/8`.

**Attack.** (a) Extend the exhaustive search to 16 primes near `2^16` and `2^17`, frozen in advance. The work per instance is
about `(p/2 exponents) × (m² pair sums)`, with early exit on the first exactness violation, so `p ≈ 2^20` is out of reach for
this code. (b) Seek a bound on `|A_E ∩ (c − A_E)|` for `A_E = {j^E : j ≤ √p/2}` via additive energy of short power images.

**Refutation criteria.** A finite example cannot refute C2R-X as quantified, since `p_0` and `K` are existential. What a
finite run can refute is the **stronger finite-range rule** "for `p ≥ 4099`, every exact fold with `≥ (2m+5)/8` internal pairs
has `(p−1) | 2(E−1)`": one exact non-character fold of that size refutes it. Refuting C2R-X itself needs, for some fixed
`ε > 0`, an unbounded family of exact folds with `≥ εL` internal pairs and `(p−1)/gcd(p−1, E−1) → ∞`. Run (a) reports the
observed maximum of `(p−1)/gcd(p−1, E−1)` among exact folds with linear pair counts, as evidence either way.

**Cost.** (a) Estimated 1–10 CPU-hours: the 32 current instances (`p ≤ 24 631`) took 38 CPU-s, and the cost grows roughly
like `p²` (about ×30 from `p ≈ 24 600` to `p ≈ 2^17`). Time one instance first. (b) Open-ended. Low priority: a proof would
only strengthen an obstruction.
