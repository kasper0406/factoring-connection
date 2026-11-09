# grdeg pre-registration (frozen before any experiment is run)

Agent `grdeg`, phase 5 restart (2026-11-09). Target: resbarrier open item #1, "degeneracy-driven G_r failure".
Code: `src/grdeg/`. Raw output: `experiments/raw/grdeg_E*.jsonl`. Exact integer / certified-ball arithmetic only;
floating point appears only in printed summaries. No operational code sees `p`, `q`, `phi`.

## 0. Background and the precise question (not under test; context)

resbarrier §2.6/§4.2 leaves open: families with rank `r ≥ 3` linked generators, **no rational side**, two nonlinear
sides. There the α ≥ 1/3 floor is closed only by the heuristic **G_r** (Schmidt fixed-degree point counts +
independence of the `r−1` linking conditions). If G_r fails, the aggregate ledger `log X ≈ x/d + (D/r)C` formally
permits L[1/4].

`rank3` (phase 17) and `rankdegen` (phase 18) already converted G_r into theorems for **general position**
(Δ_F ≠ 0: R1 forces `v_p(Δ_F) ≥ P(r,D)`, so `V_F ≤ C_r N^{−κ_r}`, κ_r = r(r−1)/(D(D−1)); R4 kills the §4.2
aggregate family for every N), and for **frame-coverable degenerate** types (frame–Hölder). What remains open:

- (a) **frame-less** degenerate modules (Δ_F = 0 with no `r+1` rows in general position): type **D** relative
  modules over a common subfield `k ≠ Q`, codex's finite-volume frame-less example, decoupled fibre products;
- (b) large boxes beyond side `B^{1/(r−1)}` (needs `LP_exc`);
- (d) the split side above threshold (needs `H_split`).

This agent targets **(a)**, the degeneracy direction the brief emphasises (CM/torsion, repeated factors,
Prym/Jacobian, subfield/tower). The question: can a frame-less type-D family actually be constructed and used to beat
L[1/3]?

## 1. Hypothesis to be developed (the obstruction), stated before coding

**Conjecture SL (subfield-linking obstruction).** A frame-less type-D linked rank-`r` module — whose candidate
module is, after Theorem-S reduction, defined over a common subfield `k` with `[k:Q] = e ≥ 2` (a Weil restriction of
a rank-`r/e` module over `k`) — can only be given concrete linking maps φ: O → Z/N, φ': O' → Z/N that agree on the
`k`-structure if one exhibits a ring homomorphism `O_k → Z/N`, i.e. a root mod N of the degree-`e` minimal
polynomial `μ_k` of a generator of `k`. For `N = pq`:

- (SL1) If `μ_k` has a root mod N that is **not** induced diagonally (same residue structure mod p and mod q), a
  gcd reveals a factor (ACCIDENTAL/STRUCTURAL split).
- (SL2) Producing any root of `μ_k` mod N is at least as hard as factoring in the worst case (`e = 2`: square roots
  mod N ⇔ factoring, Rabin 1979); for general `e` it is the polynomial-root-finding-mod-N problem, conjectured
  factoring-hard.

Hence the degeneracy that would correlate the `r−1` linking conditions (the only way, by R1, to force Δ_F = 0 and
escape the rank-r ledger) embeds a factoring-equivalent sub-problem in its own **construction**, so it is not a free
lunch. This is the rankdegen §5 conjecture made precise.

Status target: `proved` or `experimentally supported` for SL1; `assessment`/`conjectured` (with the Rabin reduction
cited) for SL2. Any theorem passes an independent Opus critic before being reported.

## 2. Experiments (frozen)

Instances are seeded (`seed = 20261109`). All norms are exact resultants `Res(f, A)`; Δ_F via certified `arb` balls
as in `src/rank3/forms.py` (reused by import).

### E1. Small-norm box excess for frame-less degenerate families with N^k | disc

**Goal.** Test whether a frame-less degenerate (type-D / subfield / CM) rank-`r` family produces box points with
**more / smaller** small-norm values than the rank-2 (resultant) prediction `C_2 N^{−2/(D(D−1))} X^{2/D}` would
allow — i.e. a measurable G_r violation at small sizes.

**Families built (exact, factorization of N known only to the harness that fixes the instance, never to counting
code):**
- (D1) Weil restriction over `k = Q(√D₀)`: `O = O' = Z[√D₀][θ]`, rank-`r` module `(1, θ, …)` with the same θ on both
  sides shifted by `N` on one side (the shifted-generator trick of rank3 §6 / rankdegen), chosen so Δ_F = 0 and no
  frame exists. `N` built from primes that **split** in `k` (so a root of `x²−D₀` mod N exists).
- (D2) CM / cyclotomic: `k = Q(ζ_m)`, module on `(1, ζ_m, ζ_m²)`; `N ≡ 1 (mod m)`.
- (D3) subfield tower: `K = k(α)` with `k` quadratic, module generators living in a rank-deficient `k`-subspace.
- Control: the matched **general-position** rank-`r` module (rank3 E1 style) and the rank-2 SNFS pair.

**Box.** side `E` chosen so the lattice error `rD(2E+1)^{r−1}` does **not** dwarf the box population; `N` of 30–60
bits; `X` over a dyadic grid.

**Predictions.**
- P1 (degeneracy really is frame-less): for D1–D3 as constructed, Δ_F = 0 and the number of frames (r+1 rows with
  every r independent) is 0. [kill if any D-family we build is secretly frame-coverable; then it is covered by
  rankdegen and reported as such.]
- P2 (no excess beyond the reduced module): the count of box points with `0 < |F| ≤ X` for a type-D family equals,
  up to the lattice error, the count for its Theorem-S **reduced** module of rank `r/e`, degree `D/e` — i.e. the
  degeneracy gives the *lower-rank* ledger, never a super-generic excess. Quantitatively, saturation
  `s = (count − error)_+ / main_reduced` stays `≤ 2` across the grid.
- P3 (the N-saving persists in the reduced module): the reduced module still satisfies `v_p(Δ_reduced) ≥ P(r/e, D/e)`
  at every `p | N`, so R1/R4 apply to it. [If P3 fails — a reduced module with Δ = 0 and genuine box excess beating
  rank 2 — that is a candidate G_r violation and is escalated.]

A P2 or P3 failure that survives scrutiny is a *refutation of SL* and a genuine sub-case to escalate.

### E2. Construction-cost / linking-map obstruction

**Goal.** Confirm SL1 operationally: that assembling the linking of a type-D family forces a root of `μ_k` mod N,
and that an *honestly constructed* one (without knowing p, q) either (i) cannot be produced below factoring, or (ii)
when produced, yields a split.

**Procedure (operational code sees only N and the public target k):**
- Attempt to build φ, φ' for D1–D3 using only `N` and `k` (no p, q): the step is "find a root of `μ_k` mod N".
- Instrument: does the only available construction call a square-root-mod-N / polynomial-root-mod-N routine? Record
  it as a charged sub-routine (NOTE per AGENT_RULES #2). Count ACCIDENTAL vs STRUCTURAL splits when a root is found
  by the harness (which *is* allowed p, q) and fed back — classify each `gcd(root-difference, N)` hit.

**Predictions.**
- P4: every D-family's linking construction requires a root of a degree-`≥2` polynomial mod N (no linear-only
  construction exists). [kill if a type-D linking is realisable with only linear mod-N arithmetic; that would be a
  real escape.]
- P5: for D1 with `N = pq` and a harness-supplied root `s` of `x²−D₀` mod N, the two CRT-inequivalent roots `±s`
  and the "mixed" roots `s'` give `gcd(s − s', N)` a proper factor in the mixed case and trivial otherwise — the
  expected Rabin behaviour, confirming SL1. Report STRUCTURAL-split rate vs the null random-gcd rate `≈ 1/p`.

## 3. Rules

- Pre-registration frozen; deviations disclosed in the results note.
- Reuse `src/rank3/forms.py` and `src/rankdegen/fh.py` by import; do not edit them.
- Operational/constructive code receives only `N`, `k`, public config, seed; never `p`, `q`, `phi`, and never calls a
  factorization routine. Any primality or root-mod-N call is named and charged.
- Status labels on every claim: proved / experimentally supported / conjectured / refuted for this model / UNVERIFIED.
- GPT-6.1-sol blind round `codex/grdeg_1`, review round `codex/grdeg_2`. Any new theorem passes a read-only Opus
  critic before being reported.
