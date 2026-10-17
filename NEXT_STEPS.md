# Next steps

Prioritized. Each item has a refutation criterion fixed in advance and an estimated cost. Context: `REPORT.md` §1, §5.
Freeze each item's parameters in a pre-registration note before running on its instances, as was done for C-2R.

## 1. Generic-ring status of the first non-floor bit (decisive; theory)

**Question.** In the ring-plus-gcd model (opaque residues mod N; ring operations; equality and gcd tests), is deciding
`ε = (−1/p)` for semiprimes `N ≡ 1 (mod 4)` equivalent to factoring? Equivalently: bit 3 of `φ(N)`, or "N is a sum of two
squares".

**Attack.** Adapt Jager–Schwenk (eprint 2009/621 full version, Thms 2–3) from a random input element to the fixed element `−1`.
With a constant input the algorithm is a constant-only SLP with equality tests. Combine transcript freezing, the constant-only
SLP bound and the CRT-disagreement equivalence (`notes/codex_r2_generic_ring.md` §1, §3). Target lemma: a `T`-step generic
algorithm with advantage `δ` yields a factor with probability `≥ δ/poly(T)`. GPT-6.1-sol drafts blind; an independent Opus
critic audits.

**Refutation criterion.** The equivalence is refuted by a generic algorithm, on an explicit semiprime distribution, that
decides `ε` with advantage `≥ 1/poly(n)` and splits with probability `o(advantage)`. The proof attempt is abandoned if
transcript freezing provably loses more than `poly(T)` (exhibit an adaptive program on which every frozen transcript has
exponentially small mass).

**Cost.** 1–3 agent-days plus 2 GPT-6.1-sol rounds. No compute.

## 2. Does the bit tower continue as class-group invariants? (cheap experiment)

**Question.** GPT-6.1-sol r4 identifies bit 5 of `s = p + q` (for `N ≡ 1 mod 8`) with `r4 > 0` of `Cl(−8N)` or `Cl(−4N)`, and
calls any identification of later bits with 8-ranks UNVERIFIED. Is bit 6 of `s` (or `uv mod 4`) a function of
`(N mod 64, ε, QR_2, r4, r8)` of the two discriminants?

**Attack.** Extend `src/repair/verify_r4.py`: add `|G[8]|` counts (one more squaring per class) and tabulate `s mod 64` against
the tuple for all pairs `p < q < 1000` with `N ≡ 1 (mod 8)` (3 418 pairs). Repeat with `Cl(−4N)`, `Cl(−8N)`, `Cl(−N)` and the
public multipliers `−4kN` for small `k`.

**Refutation criterion.** Two pairs with identical tuples and different `s mod 64` refute "bit 6 is determined by these
invariants". Otherwise report the exact functional form and attempt a proof via Rédei symbols (Stevenhagen arXiv:1806.06250).

**Cost.** Under 1 CPU-hour (bound 1000 currently takes 45 s for `r4`). Half an agent-day.

**Caveat.** A positive answer is still not an algorithm: no way is known to compute `r4` or `r8` from an unfactored
discriminant faster than the class group.

## 3. Pre-registered correlation battery for `ε` and `QR_2` (cheap falsification)

**Question.** Does any statistic computable from `N` by the project's existing code correlate with `ε` (for `N ≡ 1 mod 4`) or
with `QR_2` (for `N ≡ 1 mod 8`)?

**Attack.** Freeze a menu before running: Track A section values at fixed public `(M, c)`; Track B sketch coefficients, ranks
and residues mod small primes; BM lengths over a fixed prefix; integer-representative statistics (continued-fraction period
features of `√N`; counts of reduced forms of `Cl(−4N)` with `a ≤ 50`). Run on `10^5` balanced semiprimes per size at 32, 48
and 64 bits, with seeds disjoint from all previous runs.

**Refutation criterion.** Fixed in advance: if no statistic reaches a Bonferroni-corrected `5σ` correlation at all three
sizes, the menu is refuted. Any hit is first checked for leakage, then re-run on fresh seeds.

**Cost.** A few CPU-hours; 1 agent-day including the pre-registration.

## 4. Coarse prefix predicates in the Track C models (lemma)

**Question.** Coppersmith only needs a *coarse* predicate, correct whenever `|K − p| > N^{1/4−ε}`
(`notes/trackE_analysis.md` §6). Track C's obstructions O1–O2′ are stated for exact interval zero sets. Do they survive for
coarse predicates?

**Attack.** Restate O1 (fixed-dimension powers) and O2 (one power-map layer) with a "don't-care" window of width `W` around
`p`. Check exhaustively for small `p` (`p < 200`, all `W ≤ √p`) using the existing Berlekamp–Massey and fold-search code.

**Refutation criterion.** A construction in the O1 model with order `o(p)` or in the O2 model with `deg Q < ⌈(L − W)/2⌉`
refutes the "obstructions transfer" hypothesis. A proof with the window removed from the bound confirms it.

**Cost.** 1 agent-day of mathematics; exhaustive checks under 1 CPU-hour.

## 5. Conjecture C2R-X at larger p and the missing additive bound (lemma + experiment)

**Question.** Exact two-layer folds with `≥ εL` internal pairs come only from exponents with `(p−1) | k(E−1)`, `k ≤ K`
(`notes/repair_C2R_results.md` §3.5).

**Attack.** (a) Extend the exhaustive search to 16 primes near `2^16` and `2^17`, frozen in advance. The work per instance is
about `(p/2 exponents) × (m² pair sums)`, with early exit on the first exactness violation, so `p ≈ 2^20` is out of reach for
this code. (b) Seek a bound on `|A_E ∩ (c − A_E)|` for `A_E = {j^E : j ≤ √p/2}` via additive energy of short power images.

**Refutation criterion.** One exact non-twisted fold with internal pairs `≥ (2m + 5)/8` at any `p ≥ 4099` refutes C2R-X in its
current form.

**Cost.** (a) Estimated 1–10 CPU-hours: the 32 current instances (`p ≤ 24 631`) took 38 CPU-s, and the cost grows roughly
like `p²` (about ×30 from `p ≈ 24 600` to `p ≈ 2^17`). Time one instance first. (b) Open-ended. Low priority: a proof would only strengthen an obstruction.
