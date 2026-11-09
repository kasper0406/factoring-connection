# grdeg: can the G_r (rank ≥ 3, no rational side) case beat L[1/3]? — degeneracy-driven G_r failure

Agent `grdeg`, phase 5 restart, 2026-11-09. Target: resbarrier open item #1. Revised after a GPT-6.1-sol blind round
(`codex/grdeg_1.out`) and a read-only Opus critic, which **converged** and together forced a major revision (§6, §7).

**Files.**
- Code: `src/grdeg/ledger.py` (exact exponent ledger, deliverable #1), `src/grdeg/construct.py` (frame-less type-D
  modules; reuses `src/rank3/forms.py` + `src/rankdegen/fh.py`), `src/grdeg/relres.py` (Proposition RR, the relative
  resultant), `src/grdeg/obstruction.py` (eigenvector lemma + Rabin remark), `src/grdeg/run_experiment.py`.
- Tests: `tests/test_grdeg.py`.
- Pre-registration: `notes/grdeg_prereg.md` (commit `b8c0209`, before any experiment).
- Raw data: `experiments/raw/grdeg_E{1,2}.jsonl`.
- Codex GPT-6.1-sol: blind `codex/grdeg_1.{md,out}`; review `codex/grdeg_2.{md,out}`.

**Labels.** proved / experimentally supported / conjectured / refuted for this model / assessment / KNOWN / UNVERIFIED.

---

## 0. Verdict

**Plain language. No. The rank-≥3 no-rational-side ("G_r") case does not beat L[1/3] in any mechanism we could
analyse, and — correcting the first draft — the reason the one open door (frame-less degenerate / subfield / Weil
restriction modules) stays shut is not that the linking is factoring-hard to construct (it often is NOT: both the
Opus critic and the codex blind round exhibit frame-less families with cheap, factoring-free linking maps). The real
obstruction is a HEIGHT penalty: a frame-less subfield degeneracy over a field `k` of degree `e` with a nonzero
relative resultant pays, through `N | Norm_{k/Q}(R)`, a generator-height floor `e·x/d`, which for the corresponding
rank-`2e` module is exactly the `~ r·x/D` penalty that restores α ≥ 1/3. The only escape from the penalty makes the
relative resultant vanish — collapsing both sides to the same form, which yields only trivial congruences of squares
(useless). This is Lemma RES (the ordinary resultant bound) relativized to `k`.**

Precisely:
1. **What G_r is and what must fail (deliverable #1, `ledger.py`) [proved, exponent arithmetic].** At the exponent
   level, write `d = x^δ`, `r = x^ρ`, `log B = x^α`; candidates `a_j ∈ {0,±1}` give `C = r log 3 ~ x^ρ`, and
   cost `= L[α]` forces `ρ = α`. The smoothness target is `log X = κ_gen·x + (D/r)C`. The achievable factoring
   exponent is `α = 1/(4−c)` where `κ_gen` has x-exponent `1 + c·ρ − δ`: **`c = 0`** (`κ_gen = 1/d`, independent of
   r) gives the aggregate **L[1/4]**; **`c = 1`** (`κ_gen = (r−1)/(D−1) ~ (r/2)/d`) gives **L[1/3]**
   (`aggregate_alpha`, `r4_alpha`, `crossover_floor`: 1/4, 4/15, 2/7, 4/13, 1/3). The claim that must fail is the
   *independence of the r−1 linking conditions*, failing toward **too many correlated families** so the N-saving
   exponent `c` drops below 1. The regime realising L[1/4]: `d ≍ x^{1/2}`, `r ≍ x^{1/4}`, `D = 2d`, box side 3.
   Codex reaches the identical regime and ledger independently.
2. **Which counting claim, sharpened (codex).** Two counting problems must be separated: **G-height** (counting linked
   modules of bounded height) and **G-count** (counting small-|F| candidates inside one module's box). G-height does
   *not* imply G-count — archimedean cancellation can give small values in a large-height module — so a genuine
   barrier needs the box-counting version. For *general position* this is rank3's volume bound; a G_r failure must
   evade it, which by rank3 R1 needs Δ_F = 0.
3. **Frame-less type-D box geometry (E1) [experimentally supported].** The Weil-restriction family D1
   (`K = K' = Q(√2,√3)`, rank-4 module over `k = Q(√2)`, one side N-shifted) is genuinely frame-less (Δ_F = 0, 0
   frames). Its small-norm box points lie only on the `a₃=a₄=0` sublattice where the two sides coincide (a rank-2
   general-position sub-family); genuinely higher-rank candidates pay `|F| ~ N^{Θ(1)}`. So the frame-less excess is a
   lower-rank-sublattice phenomenon (rankdegen LP_exc(ii)) — no new relation rank.
4. **The relative-resultant obstruction (Proposition RR, §3) [proved].** For a frame-less subfield degeneracy over
   `k` (degree `e`) with relative minimal polynomials `f, g` whose reductions share the linking root and relative
   resultant `R = Res_k(f, g) ≠ 0`: `N | Norm_{k/Q}(R)`, hence `h₁ + h₂ ≥ x/n − O(d) = e·x/d − O(d)`
   (`n = d/e`). For a rank-`r = 2e` module this restores `κ_gen ~ r·x/(2D)` and so α ≥ 1/3 in the box/height ledger.
   Verified on codex's examples (`relres.py`). If `R = 0` the two sides descend to the same form (repeated norms):
   products are squares already in `k`, giving trivial congruences of squares — useless (rankdegen Lemma CC;
   codex's `t³−t−1` example).
5. **True floor [assessment].** Within the box scope, the subfield frame-less case is re-closed at **1/3** by RR
   (height), *without* any factoring-hardness assumption. The genuinely-open sliver is unchanged and both reviewers
   agree on it: **non-subfield frame-less** configurations (rankdegen's finite-volume frame-less example; decoupled
   fibre products) and **large boxes beyond `B^{1/(r−1)}`** (LP_exc), where neither RR, nor frame–Hölder, nor the
   sublattice argument applies. These are the two doors still ajar.

**Probability assessment.** < 2% that the G_r case yields sub-L[1/3]. The single most promising remaining crack: a
*non-subfield* frame-less module (no `O_k`-action to invoke RR) whose box, larger than `B^{1/(r−1)}`, holds
small-|F| lattice points **in excess of volume** that are not confined to a lower-rank sublattice. This is rank3 §7(b)
/ rankdegen's LP_exc, restated; no construction is known (codex: "none of the constructions here meets that target").

---

## 1. What G_r is, and the exact ledger (deliverable #1)

**G_r (resbarrier §2.6).** (A) Schmidt's fixed-degree point counts applied at growing d, r to norm-form heights; (B)
the r−1 linking conditions (a common linear factor `ℓ ∈ P^{r−1}(Z/N)`) are equidistributed and **independent** mod N.

**Ledger (`ledger.py`) [proved arithmetic; endpoints heuristic].** With the parametrisation above,
`log X = κ_gen·x + (D/r)C`, entropy term `(D/r)C = 2d log 3 ~ x^δ`, and CEP smoothness `u log u ≤ C` (`u = log X/log B`)
give, after optimising δ, `α = 1/(4−c)` for `κ_gen`-exponent `1 + c·ρ − δ`. Endpoints `c=0 → 1/4`, `c=1 → 1/3`.

**Honest scope (Opus critic, folded in).**
- The *arithmetic* (`aggregate_alpha`, `r4_alpha`, `crossover_floor`, `ρ=α`, δ-optimisation) is **proved**.
- The *endpoints* α=1/4 and α=1/3 as factoring exponents are **heuristic**: the ledger form `log X = κ_gen·x + (D/r)C`
  is rank3's *volume* ledger, valid inside the box scope H5_r. The L[1/4] regime (`a_j ∈ {0,±1}`, side 3) has lattice
  error `rD·3^{r−1} ≫ B`, so H5_r **fails** and the volume ledger is not proved there; turning rank3 R4
  (which bounds generator heights) into box-norm sizes needs resbarrier's typical-size/balanced-profile heuristic.
  So both endpoints are conditional on that heuristic, as resbarrier §2.6 already flags.
- **"R1 pins c = 1, and c < 1 only if Δ_F = 0" is downgraded to an assessment.** It is the clean statement for
  general position, but (i) large-box lattice excess (LP_exc, rank3 §7(b)) is a *second* door to `c < 1` even with
  Δ_F ≠ 0, and (ii) frame-coverable Δ_F = 0 configurations exist and are also box-scoped. The precise relation is:
  rank3 R1/R4 and G_r **imply the same exponent** (both N-saving exponents grow like `r/d` up to a factor 2), not that
  one is literally "equivalent" to the other.

**G-height vs G-count (codex).** G_r's module-count (how many linked pairs of a given height) is weaker than what a
barrier needs, which is a box small-value count. The gap is archimedean cancellation. For general position rank3's
volume bound supplies the box count; a degenerate escape must defeat it.

## 2. Frame-less type-D box geometry (E1)

**Construction (`construct.py`).** D1 = Weil restriction `K = K' = Q(√2,√3)` (θ = √2+√3, min poly `x⁴−10x²+1`),
rank-4 `O_{Q(√2)}`-module `(2, 2√2, 2√3, 2√6)`, √3/√6 generators N-shifted on side 2 for genuine mod-N linking; plus
D2 (cyclotomic `Q(ζ₅)`) and a GP control.

**Findings (exact; `grdeg_E1.jsonl`).** D1 is genuinely **frame-less** (Δ_F = 0, 0 frames, 34/70 minors zero; P1
held). D2 and the GP control stay **framed** — a single-generator N-shift of a Vandermonde/CM module does not
degenerate it (so cyclotomic/CM structure per se did not create a frame-less degeneracy here; D1's degeneracy comes
from the two sides sharing their first two generators). **Small-norm confinement [P2 held]:** D1's small `|F|` occur
only on the `a₃=a₄=0` sublattice (sides coincide, W=W': a rank-2 GP sub-family, `test_D1_coincidence_sublattice_is_
rank2_gp`); `min|F|_genuine/N` = 6.3e22, 2.0e32, 5.3e41 at 23/34/44-bit N (grows with N).

**Disclosed deviation.** The pre-registration's P2/P3 described a Theorem-S reduction to a rank-`r/e` module with
perfect-square norms. That fits *within-side* subfield degeneracies; D1's degeneracy is **cross-side**, so its norms
are not perfect squares (the square test was wrong and removed). The operative, confirmed statement is the
sublattice-confinement one (rankdegen LP_exc(ii)) — a change of description, not of conclusion.

## 3. The relative-resultant obstruction (Proposition RR)

This is the corrected main result, replacing the first draft's factoring-hardness claim. It is codex's §2 argument,
made exact here and verified (`relres.py`); it is Lemma RES (snfsrig §1.1) relativized to `k`.

**Proposition RR [proved].** Let the frame-less degeneracy be a module over a common subfield `k = Q(w)`,
`[k:Q] = e ≥ 2`, realised on `K = k(s)`, `K' = k(s')` with relative minimal polynomials `f, g ∈ O_k[x]` of relative
degrees `n = n' = d/e`, and linking maps whose common reduction `ψ: O_k → Z/N` (surjective, `1 ↦ 1`) sends `f, g` to
polynomials over `Z/N` sharing the linking root `v`. Let `R = Res_k(f, g) ∈ O_k`.
- (RR1) `ψ(R) = Res_{Z/N}(ψf, ψg) = 0` (resultant commutes with the coefficient map; `ψf, ψg` share root `v`). So
  `R ∈ ker ψ`, an ideal of index `N`, whence **`N | Norm_{k/Q}(R)`**. [proved; `relres.psi_resultant_is_zero`,
  `relres.rr_quadratic`.]
- (RR2) By the resultant product formula / Hadamard, `|Norm_{k/Q}(R)| ≤ exp(n'h₁ + n h₂ + O(d))`, so
  `h₁ + h₂ ≥ x/n − O(d) = e·x/d − O(d)` (`h_i` = relative archimedean generator heights). [proved.]
- (RR3) A rank-`r = 2e` relative binary module therefore has N-saving exponent `≥ e·x/d = (r/2)·x/d`, i.e.
  `κ_gen`-exponent `c = 1` up to the factor-2 that separates SNFS from GNFS. In the box/height ledger (§1) this gives
  **α ≥ 1/3**, with the same typical-size caveat as §1. [proved (model).]
- (RR4) If `R = 0`: `f, g` share a factor over `k`, so the two sides descend to the same relative form; a product that
  is a square on both sides is already a square in `k`, and its two square roots have equal image up to sign —
  a trivial congruence of squares. Useless (rankdegen Lemma CC; codex's `k=Q(t), t³−t−1`, `F = G⁴` example). [proved.]

**Verification (`relres.py`).**
- codex's `N=65`, `k=Q(i)`, `f=x²−14`, `g=x²−79`: `R = 65²`, `Norm = 65²·274625`, `N | Norm` ✓; `ψ(R)=0 mod 65` ✓.
- codex's cheap-linking "escape" `k=Q(√−N)`, `f=x²−(1−N)`, `g=x²−(1−2N)`: `R = N²`, `Norm = N⁴`, `N | Norm` ✓. Here
  the linking maps `t↦0, s↦1, s'↦1` are written down in `poly(log N)` time with **no root extraction** — yet the
  generators `s, s' ~ √N` have height `x/2`, so the box norms are `~ N⁴` (confirmed scaling): the obstruction is
  **height, not construction time**, exactly as RR predicts and as the first draft got wrong.

**Relation to the (now-demoted) factoring-hardness remark.** The original SL claimed the linking is factoring-hard to
build. Both reviewers refuted this as a blanket claim. What survives, and is rigorous, is only:
- (SL0, non-circular form, Opus critic; `obstruction.eigenvector_lemma`) [proved]. If `φ` is a ring map and `θ_k`
  acts on the module by an integer matrix `C`, multiplicativity forces `C v = s v (mod N)` with `s = φ(θ_k)`, and
  `μ_k(C) = 0` gives `μ_k(s) ≡ 0 (mod N)`: a root of `μ_k` mod N is forced at a unit coordinate — *no assumed
  `O_k`-structure on `Z/N`*. (Verified: `Cv=sv`, `μ(s)=0`.)
- (SL1) [proved]. For `N=pq`, `Hom(O_k, Z/N) ↔ Hom(O_k,F_p)×Hom(O_k,F_q)` (CRT; roots of `μ_k` mod each prime, with
  the Dedekind caveat `p ∤ [O_k:Z[θ]]`); two maps differing at one prime give `gcd(s−s',N) ∈ {p,q}`
  (`rabin_split_harness`: STRUCTURAL by construction).
- (SL2) [conjectured, NOT proved]. Producing *one* root of a *fixed* `μ_k` mod N (e.g. a single `√2 mod N`) is
  **not** known to be factoring-equivalent — Rabin's reduction needs *random* squares, and a deterministic solver may
  always return the same CRT class. So the first draft's "for the quadratic subfield is factoring itself" is withdrawn;
  it is at most a conjectural hardness, and RR shows it is not needed anyway.

**Why RR subsumes SL.** If one tries to avoid RR's height penalty by taking small (`N^{o(1)}`) generators, then `R`
is small, so `N | Norm_{k/Q}(R)` forces `Norm(R) ≥ N` — impossible for small `R` unless `R = 0` (RR4, useless). So
small genuine linking is excluded by RR directly, with no appeal to the hardness of root extraction.

**Scope (honest; both reviewers).** RR covers degeneracies with a genuine common subfield `k` and nonzero relative
resultant. It does **not** cover: (a) non-subfield frame-less configurations (rankdegen's `(1,s,t,t²)`,`(1,s+23,t,t²)`
finite-volume example; decoupled fibre products), where there is no `O_k`-action and no relative resultant to charge;
(b) the **N-dependent subfield** / tower case, where `μ_k` is built by base-`m` expansion so that the root mod N is
free (height `N^{1/e}`): this is tower-NFS-like and, like SNFS, sits at L[1/3] by paying the subfield-generator height
`N^{1/e}` (consistent with RR's height reading), but it is a separate ledger not covered by RR's statement. We exhibit
no family beating L[1/3] in any of these.

## 4. Relation to prior notes

- rank3 R1/R4 (general position) + rankdegen frame–Hölder (frame-coverable Δ_F=0) + **grdeg Proposition RR
  (subfield frame-less Δ_F=0)** now cover all degenerate types *except* non-subfield frame-less and large boxes.
- RR is the subfield instance of snfsrig Lemma RES / resbarrier Lemma RB2, charged over `O_k`.
- Consistent with resbarrier §4.2's wish ("a structure such as CM/Galois symmetry making the r−1 conditions dependent
  without the geometric height growth of powers — we found none"): the subfield structure that creates the dependency
  is exactly what pays the relative-resultant height, by RR.

## 5. Reproduction

```bash
.venv/bin/python -m src.grdeg.ledger                        # deliverable #1 table
.venv/bin/python -m src.grdeg.construct                     # D1/D2/GP geometry
.venv/bin/python -m src.grdeg.relres                        # Proposition RR on codex's examples
.venv/bin/python -m src.grdeg.obstruction                   # SL0 eigenvector lemma + Rabin remark
.venv/bin/python -m src.grdeg.run_experiment all            # writes grdeg_E{1,2}.jsonl
.venv/bin/python -m pytest -q tests/test_grdeg.py
```

## 6. GPT-6.1-sol blind round (`codex/grdeg_1.out`, before/independent of our analysis)

Codex independently (a) derived the same L[1/4] ledger and regime (`d≍√(x/log x)`, `r≍x^{1/4}(log x)^{3/4}`,
`y≍√(x log x)`); (b) separated G-height from G-count; (c) built two explicit frame-less families with **cheap**
linking maps (refuting the first-draft factoring-hardness obstruction); (d) supplied the **relative-resultant**
obstruction now adopted as Proposition RR; (e) flagged that non-subfield frame-less families are not excluded by the
cited literature, and that Bean–Thunder/Thunder give volume/finite-type bounds whose passage to finite-box counts at
growing r, D needs extra input. All adopted.

## 7. Opus critic (read-only) and the revision it forced

The critic verified the ledger arithmetic and the Rabin/ledger numbers, and found (all applied):
- the ledger endpoints are **outside rank3's box scope** (H5_r fails for `a_j∈{0,±1}`), hence heuristic, not proved;
- "R1 pins c=1 / only Δ_F=0" overclaims — LP_exc large-box excess is a second door; "equivalent" → "implies the same
  exponent";
- **SL0 was circular** as phrased; replaced by the eigenvector lemma (non-circular, now `obstruction.eigenvector_lemma`);
- **SL2 (e=2) was a WRONG label**: one fixed `√D₀ mod N` is not known factoring-equivalent (Rabin needs random
  squares) — downgraded to conjectured and shown unnecessary once RR is in place;
- the §3 trichotomy was non-exhaustive — the **N-dependent subfield (tower)** case was missing (added to §3 scope);
- the "naive root search" E2 numbers are vacuous evidence (bound below √N; cost Θ(N) not Θ(√p)) — relabelled as an
  illustration only, not evidence.

Both reviewers independently reached the same correction: the obstruction is **height (relative resultant), not
construction hardness**, and the non-subfield frame-less / large-box case stays open.

## 8. Reviews pending
- GPT-6.1-sol review round on this revised note (`codex/grdeg_2.out`): §8 (to be folded in).
- Opus critic re-vet of Proposition RR: appended below once run.
