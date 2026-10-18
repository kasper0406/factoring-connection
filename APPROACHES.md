# Approaches catalog

A catalog of everything tried in this project, from the Phase 1 tracks (2026-10-17) to Theorem D′
and the Phase 5 threads (2026-10-18). Each entry gives a one-line idea, a verdict with its
status label, the key numbers, and pointers to notes, code and commits.

**Conventions.**
- Status labels follow `notes/AGENT_RULES.md`: **proved**, **experimentally supported**,
  **conjectured**, **refuted for this model**, **open / not yet tested**, **UNVERIFIED**.
  "Refuted for this model" means refuted in the stated, restricted model, not in general.
- Commit SHAs are the *current* SHAs after the 2026-10-18 history rewrite. Notes often cite
  pre-registration commits by their *original* SHA; those appear as "orig `xxxxxxx`", and the
  full mapping is in [`COMMIT_SHA_MAP.md`](COMMIT_SHA_MAP.md).
- "codex rN" means a GPT-6.1-sol round run through the Codex CLI. Prompts and raw outputs are in
  `codex/`; curated copies are in `notes/codex_*.md`.
- Paths are relative to the repository root.

**Headline.** No polynomial-time or subexponential factorizer was found. The one positive
result of substance is **Theorem D′**: a deterministic algorithm that factors all but an
O(N^{−ε}) fraction of balanced semiprimes (p < q ≤ 2p) in N^{1/6+5ε/6}·polylog bit operations.
It is average-case only. A paper draft is in [`paper/`](paper/).

---

## Contents

- [0. Map by phase](#0-map-by-phase)
- [1. Walls and obstructions found](#1-walls-and-obstructions-found)
- [2. Refuted ideas](#2-refuted-ideas)
- [3. Positive results](#3-positive-results)
- [4. Open](#4-open)
- [5. Reviews and literature rounds (codex r1–r13)](#5-reviews-and-literature-rounds-codex-r1r13)
- [6. User ideas](#6-user-ideas)
- [7. Lessons](#7-lessons)

---

## 0. Map by phase

| Phase | Work item (agent) | Main entries below | Notes |
|---|---|---|---|
| 1 | Brief: is there a cheap operation on a compact representation of F_N = (1+X)^N − 1 − X^N that reveals a factor? Tracks A–F | — | `BRIEF.md` |
| 1 | Seed lemmas L1–L7 | P-1 | `notes/seed_lemmas.md` |
| 1 | Track A, binary coefficient sections (infraA) | W-1, W-2, R-1…R-5, P-2, O-1 | `notes/trackA_*.md` |
| 1 | Track B/D, quotient sketches, lifts, annihilators (sketchBD) | W-3, W-4, R-6…R-9, P-3 | `notes/trackBD_*.md` |
| 1 | Track C, factorial unit predicate (factC) | W-5…W-9, R-10…R-14, P-4 | `notes/trackC_*.md` |
| 1 | Track E, partial-information decoders (trackE) | W-10, P-5, P-6 | `notes/trackE_*.md` |
| 1 | Track F (no dedicated agent; codex cards) | R-15, O-2 | `notes/codex_r1_mechanisms.md`, `notes/codex_r4_sos_classgroup.md` |
| 1 | Repair of the r3 audit; C-2R cycle (repair) | W-11, R-16, O-3 | `notes/repair_*.md` |
| 1 | Final report | — | `REPORT.md`, `NEXT_STEPS.md` |
| 2 | ε/QR₂ correlation battery (epsbattery) | R-17, W-12 | `notes/eps_*.md` |
| 2 | s mod 128 tower (tower128) | P-7, O-4 | `notes/tower_results.md` |
| 2 | Primality-test leakage (primality) | W-13, W-14, R-18, P-8 | `notes/primality_*.md` |
| 2 | p−1 stage 1 + BSGS (pm1bsgs, user idea) | R-19 | `notes/pm1bsgs_analysis.md` |
| 2 | Cyclotomic Frobenius handles (cycloharvey) | W-15, R-20 | `notes/cycloharvey_*.md` |
| 3 | CRT sign-splitting (signsplit, user idea) | R-21, P-9 | `notes/signsplit_*.md` |
| 3 | Several moduli / multi-base (multimod) | W-16…W-18, R-22, P-10, O-5 | `notes/multimod_*.md` |
| 4 | Umans–Wang divisor covers (divcover) | W-19, P-11, O-6 | `notes/divcover_*.md` |
| 4 | Quasi-Riemann Hypothesis consequences (codex r11) | W-20 | `notes/codex_r11_qrh_factoring.md` |
| 4 | Theorem D′ (dcltheorem) | P-12 | `notes/dcltheorem_*.md` |
| 5 | Paper, sharper sieve, unbalanced ratios, worst case, below 1/6, Lean | O-7…O-11 | `STATUS.md` Phase 5 table |

---

## 1. Walls and obstructions found

Scoped statements of the form "in model M, X needs at least f". None is a factoring lower bound.

**W-1. Theorem A1, dimension conservation (Track A).**
- Idea: any rank-D algebra observable that isolates the residue class r mod M of the support
  of F_N must have large D.
- Verdict: **proved**, after correction. The original "D ≥ M" needs a nonzero index i ≥ M; in
  general D ≥ min(M, i+1). For odd balanced semiprimes, Cor. 3′ gives D ≥ min(M, p+1), flat in
  the level, so the per-level doubling does not survive. Counterexample to the original claim:
  N = 35, M = 16, r = 5.
- Numbers: exhaustive checks over F_2/F_3, 0 violations; the order-M recurrence A2 checked 45/45;
  the E1 mask corollary is tight in 38/38.
- Same linear-recurrence gap principle as Track C's O1, so one gap lemma blocks both the section
  and the prefix routes for finite-algebra evaluation.
- Pointers: `notes/trackA_analysis.md`, `notes/repair_adjudication.md`; `src/trackA/`,
  `src/circuit_sections.py`; `tests/test_trackA_sections.py`; commits `5b418c6`, `96cd690`,
  `13391f4`, STATUS correction `30fdbd8`.

**W-2. Section circuits grow exponentially (Track A, A-1/A-2).**
- Idea: descend to a support index by even/odd coefficient sections with zero testing.
- Verdict: **refuted for this model** (experimentally supported for gate-local sectioning; the
  exact selector is provably ×2 per level).
- Numbers: nodes ×3.6–3.7 per level (×3.0 with Karatsuba); 32-bit levels 0/3/6/9 have
  97 / 2 563 / 119 819 / 6 411 324 nodes. The exact quotient selector costs ×2.00 per level,
  Θ(N log N) total (4.7·10^7 products at 32-bit level 20). Zero testing is *not* the obstacle
  for odd balanced N (0 scalar false zeros in 44 110 sections; S_{0,4}(F_51) is a trap for
  other N).
- Pointers: `notes/trackA_analysis.md`, `notes/trackA_candidates.md`; `src/trackA/`; `5b418c6`.

**W-3. Single-sketch readouts succeed at Θ(1/p) (Track B).**
- Idea: reduce F_N modulo X^r − c (quotient sketches) and read a factor from content, resultant,
  rank or gcd.
- Verdict: **proved** exact criterion (B1: S ≡ F_m^{p^k} mod (p, f)), with per-sketch success
  = #(non-public F_p-roots of F_q)/p (B6, formula proved, Θ(1/p) experimentally supported).
  Correction (r3): "model M" itself has **no** proved lower bound; it contains Pollard p−1.
- Numbers: B1 0 violations for all composites < 5000; B3 column-gcd elimination 99 000 matrices,
  0 exceptions; rank splits fall ≈×4 per +4 bits, 0 at 40 bits.
- Pointers: `notes/trackBD_analysis.md`, `notes/trackBD_candidates.md`; `src/trackBD/`,
  `src/quotient_sketches.py`; `tests/test_trackBD.py`; `cc11920`, corrections `13391f4`.

**W-4. Lifted and divided derivatives carry nothing new (Track D, D1/D2).**
- Idea: lifts mod N^t, Hasse derivatives, F_N'/N, quotient folds.
- Verdict: **proved**: all are unit multiples of G_k = (1+X)^{N−k} − X^{N−k} or determined by
  the sketch; G_k(α) = 0 is a Fermat-type order condition.
- Pointers: `notes/trackBD_analysis.md`; `cc11920`.

**W-5. O1, C-finite / fixed-dimension recurrences (Track C).**
- Idea: evaluate the prefix/interval unit predicate "(K−1)! is a unit" by matrix powers.
- Verdict: **proved**: a prefix needs order R ≥ p+1, an interval R ≥ L+1 (tight). Correction:
  the polynomial-readout corollary bounds the lifted dimension C(R+D, D), not R.
- Pointers: `notes/trackC_analysis.md`; `src/trackC/`; `6c3a25f`, corrections `13391f4`.

**W-6. O2/O2′, power-map layers (Track C).**
- Verdict: **proved**: one power-map layer needs deg Q ≥ ⌈L/2⌉ (strengthened after r3);
  k layers need #roots(Q) ≥ L/∏gcd(e_i, p−1). O2′ does *not* exclude O(log L) layers with
  bounded gcds. "No prior formulation found by us" for O2.
- Numbers: 1 693 434 configurations, bound attained for all 986 (p, L).
- Pointers: `notes/trackC_analysis.md`.

**W-7. O3 degree (Track C).** **Proved**, corrected: forces degree ≥ L only, not ≥ p
(a Pochhammer of degree 10 at p = 101). `notes/trackC_analysis.md`, `notes/repair_adjudication.md`.

**W-8. O4 subgroup prefix (Track C).** **Proved** for x^K − 1 only: a long prefix forces
(p−1) | K (Pollard p−1). Weaker than Martin–Parvardi 2022 / Burthe, so not a new obstruction.
Max ratio 0.69 for p < 20 000. `notes/trackC_analysis.md`.

**W-9. Tal 2026 does not transfer (Track C audit).**
- Idea: use Tal's sub-square-root modular factorial (ECCC TR26-211) over Z/N.
- Verdict: **refuted for this model**: no step gives a sub-√ unit test with unknown factors;
  steps 1–14 need a known p, the BGS parts survive at √ cost. Checked against the PDF.
- Pointers: `notes/trackC_tal_audit.md`, `notes/codex_r1_tal_bgs.md`.

**W-10. Information floor for per-modulus data (Track E, F1/F2).**
- Idea: which functions of (p mod m, q mod m) carry information beyond N mod m?
- Verdict: **proved**: informative iff non-constant on the fibre {(u, c/u)}. N mod m, Jacobi
  symbols and the codex CRT discriminant masks carry 0 bits. The first informative bit of the
  cheapest decodable format (s mod 8) is "is N a sum of two squares" (ε), a known open problem.
- Pointers: `notes/trackE_analysis.md`; `src/trackE/info_floor.py`; `96cd690`.

**W-11. C-2R two-layer folding compresses only via the factor (repair cycle).**
- Idea: fold Z_L through φ = ((x+m)^E + v)² to beat the one-layer degree m+1 by a constant
  without hidden-order exponents.
- Verdict: **refuted for this model** (pre-registered). All 14 passes use E ≡ (p+1)/2, which
  encodes (p−1) | 2E−2, i.e. the factor (50/50 factor demo). Thm C2R-3 (**proved**): squaring
  exponents give no constant-factor compression.
- Numbers: arms 0/32, 14/32, 0/32; on the 16 largest instances (p ≥ 4099) no non-character
  exponent beats m+1; frozen extension X1 0/16. Later (r5): p = 257, E = 70 gives an exact fold
  d = 7 < 9, so "only the character exponent compresses" holds only at the 3L/8 threshold or
  on large instances.
- Pointers: `notes/repair_C2R_prereg.md` (prereg `795e681`, orig `a8920b7`),
  `notes/repair_C2R_results.md`; `src/repair/c2r_*.py`; `tests/test_repair_c2r.py`;
  `src/repair/verify_r5.py`.

**W-12. CF-prefix ε hits are the classical CFRAC relation family (epsbattery, post hoc).**
- The n = 32 hits are the A_{i−1}² ≡ (−1)^i Q_i relations on small Q_i (**conjectured**
  mechanism), decaying like a power of N and gone by n = 48. Positive control PC3 (CFRAC linear
  algebra) certifies 92% of ε = +1 at n = 64 with T = 16 384 and 0 at n = 96 (cost L[1/2]).
- Pointers: `notes/eps_results.md`.

**W-13. Primality tests leak only at rate ≈ 2g/√N (primality).**
- Idea: does a failed Fermat / MR / SS / Lucas / Frobenius / AKS test expose a factor?
- Verdict: **proved** CRT criterion (a factor appears iff the residue vanishes in exactly one
  component) and exact leak formulas; rate ≈ 2g/√N with g = gcd(p−1, q−1) (experimentally
  supported). AKS: Δ_N ≡ τ_p(Δ_q) mod p, with a coefficient leak confined to q − p ≲ r, the
  Fermat regime.
- Numbers: exhaustive over every base for all semiprimes < 20 000, 3-prime N < 20 000 and
  Carmichael ≤ 200 000, 0 failures; 32 bits 55 expected / 58 observed; 0 at 48–128 bits.
- Pointers: `notes/primality_analysis.md`, `notes/primality_explainer.md`,
  `experiments/analysis/primality_tables.md`; `src/primality/`; `tests/test_primality.py`;
  `caf656a` (orig `be7723d`), `69f5256`, `ce1e4c4`.

**W-14. N-derived exponents are q-driven (primality, PT-5).**
- **Proved**: gcd(N^k − 1, p^k − 1) = gcd(q^k − 1, p^k − 1); success ≈ p^{−deg Ψ}. Only smooth
  exponents (p±1, ECM), one-sided collisions (Harvey) or congruences of squares combine
  usefully, all prior art. Thm 10: witnesses are samplable, so a witness→factor map is itself a
  factoring algorithm. `notes/primality_analysis.md`.

**W-15. Theorem LB, the N^{1/5} floor with one modulus (cycloharvey).**
- Statement: any bounded-degree identity-model covering family, unguessed or with guesses
  modulo a single modulus, costs ≥ N^{1/5}·M^{−o(1)}.
- Verdict: **proved** in the model; one independent Opus review.
- Numbers: calibration gate exact LP 1/4, 1/3, 2/9, 1/5, 1/5 (BSGS, Lehman, Hittmeir, Harvey,
  HH); operational Harvey list = model count 112/112.
- Pointers: `notes/cycloharvey_results.md` §3.4; `src/cycloharvey/`; prereg `d6e5656`
  (orig `1d92a48`), results `c1b282a`.

**W-16. Theorem A, the N^{1/6} floor with any number of moduli (multimod).**
- Statement: in the enlarged one-base sumset model 𝔐⁺ (any moduli, any classes), every covering
  family costs ≥ c·N^{1/6}. Tight: DCL's targets meet it with equality.
- Verdict: **proved** in the model (Opus critic; r12: correct in the stipulated model).
- Pointers: `notes/multimod_results.md` §3.2; `paper/main.tex` §5; results `55fb688`.

**W-17. Theorems A′ and B: below 1/5 needs N^{Ω(1)} moduli (multimod).**
- A′: degree ≥ 2 targets never beat N^{1/5}. B: with J moduli, X ≥ c N^{1/5}(J N^{o(1)})^{−1/5}.
- Verdict: **proved** in the model (B's short counting proof supplied by r12).
- Pointers: `notes/multimod_results.md` §3.1–3.3.

**W-18. Theorem C: CRT is guess-neutral (multimod).**
- Class-lattice indices multiply over coprime moduli (**proved**; 20/20 HNF checks), so schemes
  covering *all* classes of a composite modulus stay neutral (heuristic consequence). Only
  *selective* covering escapes. This also explains why sign-splitting (R-21) and Gray-code
  bucketing cannot help. `notes/multimod_results.md` §3.4.

**W-19. He–Sahai at rank r, and a mass bound (divcover).**
- HS-r: v ≤ 1 + Δ(μΔ − 1) generalises HS Lemma 2.1 (**proved**); with a-priori Δ, μ bounds it
  gives nothing beyond the mass bound at rank ≥ 2. Lemma M (**proved**): beating N^{1/5} at
  α = β = 0.39 needs ≫ n^{0.61} GAP elements, each with a divisor ≥ exp(c n^{0.22}) from primes in
  (n^{0.78}, n].
- Pointers: `notes/divcover_results.md`; `src/divcover/`; `f48b36e` (orig `77d8b67`), `37d1929`.

**W-20. [QRH] gives no factoring exponent gain (codex r11).**
- Idea: what would a zero-free half-plane Re s > 11/12 (or 7/8) for all Dirichlet L-functions
  (the UNVERIFIED "Quasi-Riemann Hypothesis" preprint) imply?
- Verdict: strong auxiliary algorithms (Ankeny–Bach-type bounds, Miller φ/λ ⇔ factoring) but no
  improvement over the unconditional N^{1/5} and no deterministic subexponential method. Shanks
  under [QRH]: |D|^{6/25} or |D|^{4/17}, worse than HH.
- Pointers: `notes/codex_r11_qrh_factoring.md`; `69bb2d3`.

---

## 2. Refuted ideas

**R-1. A-3 fixed-budget selector (Track A).** Its "coin flip" was by construction. Downgraded
to **refuted for the implemented selector**. `notes/trackA_analysis.md`.

**R-2. A-4 parity descent via X^{2M} − s² (Track A).** Rank M + M, no saving. **Refuted for this
model.**

**R-3. A-5 pure-class gcd (Track A).** Θ((q−p) log N), first hit at 2M ∈ (q−p, 4(q−p)] in 25/25.
Dominated by Fermat. 338 structural vs 3 accidental hits (null ≈ 2.05).

**R-4. A-6 Hadamard masks; A-7 early stop + gcd scan (Track A).** **Refuted**; A-7 is ≈ √N.

**R-5. Track A C1–C4 as a factoring route.** The section-structure facts are proved (P-2) but
give only Fermat-like behaviour.

**R-6. B7 rank gaps over weights (Track B; codex r1 card B1).**
- **Refuted**: no uniform guarantee. P_gap = 0 for 65–86% of balanced instances; P_gap = 0 for
  101·113, 107·193, 1009·1481 at all r ≤ 16. The r1 quotient pilot found 123/10/2 factors at
  8/12/16 factor bits against a random-gcd null of 302.9/54.6/3.50. "r-independent" was not a
  theorem (N = 209).

**R-7. B8 Lucas-support gaps (Track B).** Occur when q − p ≲ r/2 (all 24 instances at ≥ 24 bits
have q − p ≤ 162). Dominated by Fermat.

**R-8. D3 Berlekamp–Massey over the weight c (Track D).** Linear complexity mod p = q mod (p−1);
**proved** bound: BM over Z/N hits a nonunit by index ≤ p + e_p − 1. Cost Θ(q−p): a weak
Fermat-gap method. The Δ^K zero predicate is not the prefix predicate (N = 391).

**R-9. M3/M3b, Q_N and Q_N − I (Track B/D).** 0 splits at ≥ 28 bits; the M3b necessity claim
was false (N = 247).

**R-10. M1 random-shift probes (Track C).** No information: hit rate 0.067 (true) vs 0.077
(false), null 0.062–0.077. **Refuted** for independent uniform shifts.

**R-11. M2/M2′ power sums (Track C).** The evaluator is circular (a unit test of (e+1)!).
**Refuted.**

**R-12. M3 truncated F_N as a factorial evaluator (Track C).** Equals U(K) (0 mismatches for
N < 600), cost linear in K, 11 000× BGS at 64 bits. **Refuted as a speedup.**

**R-13. M5 compressed discriminant (Track C; codex r1 card C1).** Õ(K) evaluator, worse than
BGS. **Refuted for this evaluator.**

**R-14. M6–M8 (Track C).** Fixed-dimension recurrences (by O1), power-map compositions (by
O2/O2′), group-law identities (only via (p−1) | K, i.e. Pollard p−1, by O4). **Refuted for these
models.**

**R-15. Track F torus exponent N − (D/N) (codex r1 card F1); CFRAC pilot (r2 card F2).**
Twist-only generalisation fails; the pilot gave 0/4 at 24/32 factor bits. Not independently
reproduced. `notes/codex_r1_mechanisms.md`, `notes/codex_r2_generic_ring.md`.

**R-16. C-2R as a constant-factor compressor.** See W-11. **Refuted for this model.**

**R-17. ε / QR₂ prediction from N alone (epsbattery).**
- Idea: does any of 94 statistics (Track A sections, Track B sketches, BM lengths, form counts,
  CF prefix of √N) predict ε or QR₂ beyond the N mod 2^8 floor?
- Verdict: **refuted for this menu** at n = 48 and 64 (|ρ| ≥ 0.022 for ε, 0.031 for QR₂). At
  n = 32, 9 ε and 5 QR₂ CF-prefix hits (|ρ| ≤ 0.094) are explained by W-12. The pre-registered
  period-closure explanation was itself refuted.
- Numbers: 10^5 instances per size; stratified permutation z, |z| ≥ 6; 7 positive controls with
  0 violations.
- Pointers: `notes/eps_prereg.md` (prereg `2b891b6`, orig `779a650`), `notes/eps_results.md`,
  `experiments/analysis/eps_battery_tables.md`; `src/eps/`; `tests/test_eps.py`; results `04d986c`.

**R-18. Iterated primality residues; ε from transcripts (primality, PT-2/PT-3).**
- Residues encode s or q − p as a discrete log (proved), so kangaroo/BSGS factor from them at
  slopes 0.24/0.24/0.23, the N^{1/4} Pollard 1978 baseline. Iterating (a^{f(s)}, deg f ≥ 2) turns
  into order finding (slope 0.52). MR/SS certify ε = +1 only by exhibiting √−1, at rate
  ≈ g²/(3N) (2^{−61} at 64 bits). Pre-registered 14-feature battery null at 32/48/64 bits (MDE
  0.029). **Refuted for this model.** `notes/primality_analysis.md`.

**R-19. p−1 stage 1 + BSGS on b^N ≡ b^x (pm1bsgs; user idea).**
- Verdict: the steps are correct, but b^N ≡ b^x ⇔ ord_p(b) | N − x, so this is order-finding
  BSGS with the window offset at N: the known p−1 FFT/BSGS continuation (Pollard 1974,
  Montgomery–Silverman 1990, Montgomery–Kruppa 2008). Prop D (**proved**): B1 < √X is dominated
  by B1 = √X, so composite-R coverage never helps at the optimum. **Refuted for this model** as a
  route below N^{1/4−o(1)}.
- Numbers: 73/73 budget-optimal choices have B1² ≥ X; 2 800 runs match the exact predictor; at
  80 bits, expected ops: hybrid 6.6e6, plain BSGS 3.1e8, rho 2.1e6.
- Pointers: `notes/pm1bsgs_analysis.md`; `src/pm1bsgs/`; `tests/test_pm1bsgs.py`; `7067c0e`.

**R-20. Cyclotomic Frobenius handles below N^{1/5} (cycloharvey).**
- Idea: in (Z/N)[X]/Φ_m, guess r ≡ p (mod m) and use σ_r as a Frobenius handle for
  higher-degree Taylor cancellation in Harvey-style search.
- Verdict: **refuted for this model**. D1 (**proved**): cyclotomic handles add no targets, since
  every F(p, q) already has the (Z/N)^* handle a^{F(1,N)}. With W-15, no bounded-degree family
  with one modulus beats N^{1/5}.
- Numbers: exact LLL search over all degree ≤ 3 exponents: slope 0.1998 vs Lehman 0.1997
  (40–80 bits); the cyclotomic version costs 10–26× the base ops.
- Pointers: `notes/cycloharvey_prereg.md`, `notes/cycloharvey_results.md`,
  `notes/codex_r7_cyclo_harvey.md`; `src/cycloharvey/`; `tests/test_cycloharvey.py`.

**R-21. CRT sign-splitting MITM for HH (signsplit; user idea).**
- Idea: the 2^t branches sharing a square class have offsets Δ0 + Σε_iΔ_i (mod M²); a
  meet-in-the-middle over signs should save 2^{t/2}.
- Verdict: the identity is **proved** (4.3·10^8 checks, 0 failures); the speedup is **refuted**
  by the pre-registered rule. Best sign-split / best HH list work = 1.021 (56 bits), 1.009
  (64 bits). A per-group MITM loses list sharing (slope 0.251 vs 0.196, 10× worse at 64 bits);
  shared-list sign collapse loses to the HH lock by m^{2/5}/2^{3/5} per prime; idealised cap
  2^{t/5}.
- Pointers: `notes/user_ideas/2026-10-18_crt_sign_splitting.md`, `notes/signsplit_prereg.md`
  (prereg `32772b1`, orig `a31a84f`), `notes/signsplit_results.md`; `src/signsplit/`;
  `tests/test_signsplit.py`; results `1468b67`, fix `9fb6efa`.

**R-22. Multi-base searches with independent or grouped bases (multimod, Thm E).** They reduce
coordinatewise to the 1-D bounds (**proved**). General known-relation lattices are open (O-5).
`notes/multimod_results.md` §3.7.

**R-23. Chirp / quadratic-in-k leaves (user idea; multimod §4).** Cheap generation saves only
Stage 2; the Stage-3 leaf count and the sumset bound are unchanged. **Refuted for this model**
as an exponent improvement. `notes/user_ideas/2026-10-18_chirp_quadratic_exponents.md`.

---

## 3. Positive results

**P-1. Seed lemmas L1–L7.** [X^p]F_N = N/p; k·a_k = 0; prefix ⇔ K ≤ p ⇔ (K−1)! a unit;
F_N' = 0; F_N ≡ 0 as a function iff N is prime or Carmichael; Lucas/Kummer support.
**Proved** (known facts; scaffolding). `notes/seed_lemmas.md`; `src/seed_checks.py`;
`tests/test_seed_checks.py`; `6c3a25f`.

**P-2. Track A section structure C1–C4.** All classes are hit for M ≤ p−1; at most 3 terms for
coprime M > p (217 772 cases); pure-class gcd criterion; no scalar false zeros for odd balanced
N. **Proved** (with r3 scope fixes). `notes/trackA_analysis.md`.

**P-3. Track B structure B1, B4, B5.** Exact sketch criterion; vanishing depends on q only via
q mod L_α; sporadic roots in anharmonic 6-orbits (classical; only the census is new). **Proved.**
`notes/trackBD_analysis.md`.

**P-4. Track C baselines.** Naive / Pollard–Strassen / BGS implemented and measured; slopes
(FFT model) PS 0.617, BGS 0.563 in K; end-to-end BGS factors 13/13, 64 bits in 29 s.
**Experimentally supported** (reference implementations). `notes/trackC_analysis.md`;
`src/trackC/`.

**P-5. Coppersmith / Howgrave-Graham decoder price list (Track E).**
- Rigorous integer test (**proved**, via May 2010 Thms 2–3); empirical threshold
  k_min ≈ n/4 − 2 + n/(4d) (**experimentally supported**, 570 cells, n = 64–512); decode 0.003 s
  (d = 9) to 52 s (d = 65) at n = 512.
- `notes/trackE_analysis.md`; `src/trackE/coppersmith.py`, `formats.py`; `tests/test_trackE.py`.

**P-6. Exact partial-information reductions (Track E).** φ(N) mod 2^k ⇔ s mod 2^k (two cosets,
k ≥ k_min + v₂(q−p), 2 decodes); s mod 8 ⇔ ε; s mod squarefree M has 2^ω ambiguity; Legendre
bits alone leave a 2^105 list at 512 bits; N mod (p−1) is complete. **Proved.**
`notes/trackE_analysis.md`.

**P-7. The s mod 128 tower (tower128).**
- T1 (**proved**): s mod 128 needs only ψ(p), the quartic character mod 16, not an order-8
  value (corrects REPORT §1.1). T3 (**proved**): for p ≡ q ≡ 7 (mod 8),
  s ≡ −(N+1) + 64·[N ≡ 1 (16)]·(r8 Cl(−4N) ⊕ r8 Cl(−8N)) (mod 128). T4 (**proved**) decides 3/8
  of the p ≡ q ≡ 1 (8) cell. Not an algorithm: r8 and fourth-power membership need the class
  group or the factors.
- Numbers: T1 0 failures on 22 666 pairs; T3 0 failures on 9 453 pairs (q < 4000); T3+T4
  resolve 29.6% of undetermined pairs, none in the QR₂-false half.
- Pointers: `notes/tower_results.md`; `src/tower/`; `tests/test_tower.py`; `99f56bc`, REPORT
  fold-in `7fbfd2e`, post-r9 `ca06f84`. Audited by codex r9 (P5's proof repaired).

**P-8. Exact leak formulas for primality tests (primality).** See W-13; the formulas themselves
are a positive, exhaustively verified result.

**P-9. Stage-2 savings for HH (signsplit; Gray-code user idea).** Shared lattice reduction plus
Gray code: 1.3–2.6× fewer ring multiplications than a fair fixed-base baseline (3–21× vs
naive). Stage-3 leaves are unchanged, so this is a constant. **Experimentally supported.**
`notes/signsplit_results.md` §6; `notes/user_ideas/2026-10-18_graycode_harvey_benchmark.md`.

**P-10. Theorem D: average-case N^{1/6+ε} coverage in the search model (multimod).**
- Idea (DCL, "divisor-cover Lehman"): Lehman/HH targets (Bp + Aq − c)/ℓ² with small fixed weights,
  over *all* primes ℓ near N^{1/6}, at the N-determined roots of Br² ≡ AN (mod ℓ); one shared baby
  list via γ = β^{∏ℓ²}.
- Verdict: **proved** in the search model: all but O(N^{−ε}) of balanced semiprimes are covered at
  cost N^{1/6+5ε/6+o(1)}. r12 found a repairable normalisation flaw and a substantive base-order
  gap, both closed in P-12.
- Numbers (pre-registered): model slope 0.171 vs HH 0.196 over 56–96 bits at ≈95% success
  (0.45× HH at 96 bits); full coverage slope 0.196; operational 128/128; verdict SUPPORTED only
  after a disclosed classifier fix.
- Pointers: `notes/multimod_prereg.md` (prereg `a1483b2`, orig `81dd28e`),
  `notes/multimod_results.md`; `src/multimod/`; `tests/test_multimod.py`; `55fb688`.

**P-11. Prime-only reformulation of Umans–Wang (divcover).**
- R1/R3 (**proved**): UW Thms 5.2/5.5 need only a *prime* cover (only primes in (n^{2δ}, n]
  after Pollard–Strassen), so He–Sahai's semiprime-based obstruction does not touch the
  integer-factoring use, even at rank 1. R2 (**proved** after repair): prefactoring can be
  replaced by a separation condition Z_p ≠ Z_q.
- Pointers: `notes/divcover_results.md`; reviews `codex/divcover_1.out`, `codex/divcover_2.out`.

**P-12. Theorem D′ and Corollary D″ (dcltheorem). The main result.**
- Statement (**proved**): for ε ∈ (0, 1/60] and n ≥ n₀(ε), the deterministic algorithm DCL-det
  sees only N, costs N^{1/6+5ε/6}·polylog bit operations on every input, outputs only a proper
  divisor or FAIL, and fails on O(N^{−ε}) of S_n = {pq : p < q ≤ 2p}, counted by pairs.
  Corollary D″: N^{1/6}·polylog with O((log N)^{−A}) exceptions.
- How: the collision mod p is exact (no order lemma); four disjoint bands guarantee one band
  coprime to q − 1, where a two-sided collision forces ord_q(2) | (A−B)(p−1) − ℓ²δ, rare by an
  elementary small-order count plus a CRT count; coverage by the large sieve.
- Reviews: GPT-6.1-sol dclt_1, an independent Opus critic, and the fresh r13 referee, all
  ACCEPT WITH FIXES (applied). Not yet checked by a human.
- Numbers (pre-registered): 0 wrong outputs and 0 lemma violations (240 lemma runs + 656
  operational runs; 3 002 clean covering targets); E slope 0.203 vs operational HH 0.226 over
  48–80 bits, **but HH is 2.3–4.9× cheaper at 40–80 bits**; full-coverage design slope 0.227 (no
  advantage); FULL ≤ 0.045% of covering targets.
- Scope: not worst case; q/p ≤ 2 only; prior art McKee–Pinch 1996 §5 is a direct precursor;
  "we found no prior formulation of this combination".
- Pointers: `notes/dcltheorem_prereg.md` (prereg `05d2bf5`, orig `3550e24`),
  `notes/dcltheorem_results.md`, `experiments/analysis/dcltheorem_summary.json`;
  `src/dcltheorem/`; `tests/test_dcltheorem.py`; theorem `bc1de32`, r13 fixes `539fdcf`;
  paper draft `paper/main.tex`.

---

## 4. Open

**O-1. E-1: binomial-specific SLPs for sections S_{r,2^k} (Track A).** No identity found (Lucas,
Chebyshev, roots-of-unity filter, holonomic, Frobenius). Became REPORT item 2 (a conditional
reduction with Agrawal–Biswas PIT). **Open.**

**O-2. Track F class-group cards (codex r4 F1–F3).** Square membership of a public class, 4-rank
positivity, ambiguous forms. All lack a polynomial extraction. The r4 claims were verified
exhaustively (6 901 pairs p < q < 1000, 0 failures; `src/repair/verify_r4.py`). **Open.**

**O-3. Conjecture C2R-X.** An exact fold with ≥ εL pairs forces (p−1) | k(E−1), k ≤ K.
**Conjectured.** NEXT_STEPS item 5 (test at p ≈ 2^16–2^17) is not done.

**O-4. Predicting ε on balanced N; the QR₂-false half of the tower.** Deciding ε is Seurin's
2-Φ/4-Hiding problem (conjectured hard); codex r6 found no reduction from factoring. Conjecture
C5: no 2-power class-group invariant decides the mod-128 bit when QR₂ is false. **Open.**

**O-5. Multi-base searches with non-split known relations (multimod Thm E(d)).** The 1-D
reduction is not proved there. **Open.**

**O-6. Low-rank GAP prime covers for Umans–Wang (divcover).** No construction found (best:
block-primorial subset sums, log H = Θ(n/log n)); exhaustive search over 198 cells (n ≤ 100,
ranks 1–3) follows the random model; next step: bound the target-prime radical of proper rank-2
GAPs. **Open.**

**O-7. Worst case for D′ (DCH).** Needs (a) the divisor-cover hypothesis DCH (every candidate p
has ℓ | Bp² − AN for small A, B and ℓ near N^{1/6}; **conjectured**, 0 uncovered in 25 exhaustive
instances at 32–48 bits), (b) large order mod the specific q on a clean band, (c) excluding
all-FULL pairs. Phase 5 agent `dchworst`.

**O-8. Sharper exceptional rate.** The large sieve gives n²/H²; the Poisson heuristic e^{−cH²/n}
would give N^{−ε} at polylog cost (**conjectured**). A two-variable (bilinear) sieve over q is
the candidate route. Phase 5 agent `sieve2d`.

**O-9. Unbalanced ratios q/p > 2.** DCL as specified costs worse than Harvey summed over all
ratios (heuristic). Phase 5 agent `unbalanced`.

**O-10. Below 1/6.** Blocked in the one-base sumset model by W-16; needs a non-sumset collision
mechanism and faster small-factor removal. Phase 5 agent `below16`.

**O-11. Lean 4 formalisation of D′'s key lemmas.** In progress (Phase 5 agent `leanform`,
`lean/`).

**O-12. Smaller items.** C-HHloc (HH restricted to a certified interval; untested); a one-sided
collision family inside a primality-style identity (untested); Elkies Z[p, q, t_p] targets
(cycloharvey §3.6, untested); NEXT_STEPS item 4, coarse prefix predicates in the Track C models
(untested).

---

## 5. Reviews and literature rounds (codex r1–r13)

| Round | Question | Main verdict | Pointers |
|---|---|---|---|
| r1 mechanisms | Six blind mechanism cards (A1 section PIT, B1 rank gaps, C1 discriminant, D1 lifted derivative, E1 CRT masks, F1 torus) | A1 ranked first (became E-1); others refuted or reduced (R-6, R-13, R-15; E1 masks later shown 0-bit by W-10) | `notes/codex_r1_mechanisms.md`, `codex/r1_mechanisms.out` |
| r1 prior art | Primary-source map | 39 entries; the reference base for later rounds | `notes/prior_art_map.md`, `notes/codex_r1_priorart.md` |
| r1 Tal/BGS | Does Tal 2026 transfer to unfactored N? | No (confirmed by factC, W-9) | `notes/codex_r1_tal_bgs.md` |
| r2 generic ring | Synthesis of generic-ring models | Useful invariant is CRT disagreement; preprocessing restriction essential; H2 (shared-list weighted collisions) ranked first, later realised as the cycloharvey/multimod line | `notes/codex_r2_generic_ring.md` |
| r3 critic | Adversarial audit of Tracks A, B/D, C | 36 WRONG/GAP items; the critic was right on 36/36 (27 changed, 9 scoped) | `notes/codex_r3_critic.md`, `notes/repair_adjudication.md` |
| r4 SOS / class groups | Track F: sum-of-two-squares bit and class groups | SOS decision ⇔ ε ≤ QR ≤ factoring; mod-32 correction = QR₂ ⇔ 4-rank; verified exhaustively | `notes/codex_r4_sos_classgroup.md` |
| r5 report review | Review of REPORT.md | 14 items, all accepted and fixed (e.g. mod-64 experiment vacuous; k_min relabelled) | `codex/r5_report_review.out` |
| r6 ε reduction | Reduction from factoring to an ε oracle? | None found; closest Seurin 2-Φ/4-Hiding | `notes/codex_r6_eps_reduction.md` |
| r6 primality lit | Why primality tests do not reveal factors | Exact CRT condition; all 13 testable claims verified by agent `primality` | `notes/codex_r6_primality_lit.md` |
| r7 cyclo Harvey | Cyclotomic handles for Harvey-style search | Handles valid; speedup UNVERIFIED; became cycloharvey (R-20) | `notes/codex_r7_cyclo_harvey.md` |
| r8 Coppersmith bridge | User's lattice ↔ exponent bridge | Claims 1–2 proved, 3 refuted; novelty refuted (Gao et al. 2512.19076); determinant bounds cannot certify < 1/5 | `notes/codex_r8_copp_bridge.md` |
| r9 tower audit | Audit of tower128 | T1–T4 correct; P5's Frobenius identification wrong, repaired | `notes/codex_r9_tower_audit.md` |
| r10 frontier | The N^{1/6} frontier and barriers | Harvey's 1/6 needs a full square-root speedup of Lehman; UW conditional; He–Sahai refutes the AP version | `notes/codex_r10_frontier.md` |
| r11 QRH | Consequences of [QRH] for factoring | No exponent gain (W-20) | `notes/codex_r11_qrh_factoring.md` |
| r12 DCL audit | Referee of multimod Thm D and Thms A–E | Coverage argument sound and repairable; base-order gap substantive; normalisation repaired; no prior almost-all result found | `notes/codex_r12_dcl_audit.md` |
| dclt_1 | Referee of D′ | Clean band closes the order gap; fixes (baby-length rule, T(T+1), cost terms) | `codex/dclt_1.out` |
| r13 final referee | Fresh referee of D′ and prior art | ACCEPT WITH FIXES; every proof step correct; finite-input output bug (returned N for prime N) fixed; McKee–Pinch §5 read; Carella 1308.2891 proof judged incomplete (referee's view) | `notes/codex_r13_final_referee.md` |
| divcover_1/2 | Blind attempt and review of prime-only covers | Prime-only reformulation proved; R2 repaired; no construction | `codex/divcover_1.out`, `codex/divcover_2.out` |

Independent read-only Opus critics also reviewed sketchBD (15 findings), trackE (7 fixes),
cycloharvey (LB scope fixes), multimod, divcover and dcltheorem. The review history up to Phase 1
is indexed in `notes/review_and_repairs.md`.

---

## 6. User ideas

| File | Idea | Verdict | Handled by |
|---|---|---|---|
| `notes/user_ideas/2026-10-18_mod_M2_lattice_T_MK.md` | Why the mod-M² lock (M) and Taylor centres (K) give T = MK = N^{1/5} | Matches HH Prop 3.3 bookkeeping; δ = 5/2 would give N^{1/7}, but δ_eff ≤ 3/2 → 1/5 | orchestrator; cycloharvey §1a |
| `…_schoof_semigroups_and_cross_difference.md` | Group-law vs coordinate-power semigroups; cross-difference law for separable BSGS | Separability is enforced by the model (proved); Elkies targets untested | cycloharvey §§3.5–3.6; multimod Thm E |
| `…_coppersmith_exponent_bridge.md` | Coppersmith lattice gives an exponent collision; m = 1 is Lehman | Claims 1–2 proved, 3 refuted; prior art Gao et al.; min over m ≥ 2 is exactly 1/5 | codex r8; cycloharvey §3.7 |
| `…_crt_sign_splitting.md` | MITM over CRT signs for a 2^{t/2} HH speedup | Identity proved; speedup refuted (R-21) | agent signsplit |
| `…_graycode_harvey_benchmark.md` | CRT fibre power tables + Gray code | Stage-2 constant only (P-9) | orchestrator; signsplit §6; multimod Thm C |
| `…_chirp_quadratic_exponents.md` | Quadratic-in-k giants at O(1) each | Stage-2 constant; refuted as exponent gain (R-23) | multimod §4 |
| `…_pm1_frobenius_bsgs.md` | p−1 stage 1 + BSGS covering composite R | Correct but known; refuted below N^{1/4} (R-19) | agent pm1bsgs |
| `…_paradigms_4_6.md` | a = b dichotomy, batch smoothness, mod-N² Fermat quotients, class-group BSGS, chirp-z Stage 3 | The mod-p² formula was wrong (corrected); each paradigm is prior art, at most a constant | orchestrator |

The cyclotomic-Frobenius idea entered through the codex r7 prompt rather than a user-ideas file
and became cycloharvey (R-20).

---

## 7. Lessons

1. **Most "compact representation" routes hit the same wall.** Finite-algebra evaluation of
   sections, prefixes and factorial predicates needs dimension or degree comparable to p (W-1,
   W-5–W-8). Every apparent signal reduced to a known method: Fermat, Pollard p−1, BSGS, or the
   CFRAC relations.
2. **The productive line was the deterministic collision-search literature.** Model lower bounds
   (W-15, W-16) said exactly where room was left (many moduli, selective covering), and the
   positive result (P-10, P-12) lives in that room. Mapping a model's floor before searching paid off.
3. **Independent adversarial review was decisive.** The r3 critic was right on 36/36 items; r12
   found the base-order gap that turned Theorem D into D′; r13 found a real output bug and the
   closest precursor (McKee–Pinch §5). Claims improved every time they were attacked.
4. **Pre-registration kept results honest.** Several pre-registered hypotheses failed (C-2R,
   sign-splitting, the ε explanation), and one verdict (multimod) passed only after a disclosed
   classifier fix. Writing the decision rule first made those calls unambiguous.
5. **Exponent wins are not practical wins.** D′ has the smaller measured exponent at 48–80 bits,
   and Harvey–Hittmeir is still 2.3–4.9× cheaper there.
6. **Prior art is closer than it looks.** Ingredients we first thought new were in HH 2022
   Prop. 3.3, McKee–Pinch 1996, Gao et al. 2025 and Hittmeir 2018. Claims are therefore limited
   to "no prior formulation found by us".
7. **Check the quantifiers.** "Almost all N", "the true factor of almost all pairs" and "every
   candidate for every N" are different statements, and confusing them was the most common
   review finding in the later phases.
