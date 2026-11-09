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

**Plain language. No sub-L[1/3] construction exists in any mechanism we could analyse — but this is "none found +
a partial rigorous barrier", not a theorem that the whole G_r case is closed. Correcting the first draft: the one open
door (frame-less degenerate / subfield / Weil-restriction modules) does not stay shut because the linking is
factoring-hard to build — it often is NOT (both reviewers exhibit frame-less families with cheap, factoring-free
linking). The real obstruction is a HEIGHT penalty: a frame-less subfield degeneracy over `k` of degree `e` with a
nonzero relative resultant pays, through `N | Norm_{k/Q}(R)` (Lemma RES relativized to `k`), a generator-height floor
`e·x/d = r·x/D` — which restores α ≥ 1/3 **only for relative-rank-2 (`r = 2e`) or relative-power-basis modules**. The
re-vet exposed the live gap: a type-D module over a *fixed small* subfield with relative rank `m ≥ 3` and non-power
generators gets only `κ_gen ~ (r/m)·x/d` from the relative resultant — the L[1/4] exponent — so it is genuinely open,
closable only by an unproved "relative R1" over `k` (blocked by the infinite volume of the acting units) or the
conjectural hardness of root-finding mod N. No construction for it is known; it is the sharpest remaining crack.**

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
   frames). Its small-norm box points lie only on the `a₃=a₄=0` sublattice where the two sides coincide — there
   `W = W' = 16(a₁²−2a₂²)²` (a *subfield restriction* to `Q(√2)` with repeated factor directions, **not** general
   position; codex), a lower-rank repeated-norm sub-family; genuinely higher-rank candidates pay `|F| ~ N^{Θ(1)}`
   (D1's true relative resultant is `Res(x²−3, (x−N)²−3) = N²(N²−12) ≠ 0`, so D1 is an RR height-penalty case). So the
   frame-less excess is a lower-rank-sublattice phenomenon (rankdegen LP_exc(ii)) — no new relation rank.
4. **The relative-resultant obstruction (Proposition RR, §3).** For a frame-less subfield degeneracy over `k`
   (degree `e`) with relative minimal polynomials `f, g` whose reductions share the linking root and relative
   resultant `R = Res_k(f, g) ≠ 0`: **`N | Norm_{k/Q}(R)`** [RR1, **proved**], hence `h₁ + h₂ ≥ x/n − d·log 2 = e·x/d
   − O(dn)` in the projective (Mahler) height `h_i` [RR2, **proved**]. If `R = 0` the two sides are `k`-isomorphic
   (type A), giving only trivial congruences of squares — useless [RR4, **proved**; rankdegen Lemma CC]. The step to a
   factoring exponent — *for a relative-binary (`r = 2e`) or relative-power-basis module* this height restores
   `κ_gen`-exponent `c = 1` and so α ≥ 1/3 [RR3, **heuristic/model**, same typical-size caveat as the ledger].
   Verified on codex's rational examples and a genuine `k = Q(i)` example where `O_k` does real work (`relres.py`).
5. **True floor [assessment], corrected after the Opus re-vet.** Within the box scope, the subfield frame-less case
   is re-closed at **1/3 only for relative rank 2 (or relative power bases)** by RR (height, no factoring-hardness).
   Three doors stay open, on which both reviewers agree:
   - **(i) Relative rank ≥ 3 over a fixed small subfield** (the new one the re-vet exposed): a type-D module over a
     fixed `k` (say `Q(√2)`, `e` constant) with relative rank `m ≥ 3` and *non-power* generators. There typical
     `log|F|` is the *max* of the generator heights, so RR alone gives only `κ_gen ~ e·x/d = (r/m)·x/d`, i.e. `c = 0`,
     the **L[1/4]** exponent. Closing it needs a *relative* R1 over `k` (an R1 analogue for the `m−1` relative linking
     conditions), which is unproved and obstructed by `V_F = ∞` (the units of `O_k` act), or the conjectural SL2.
   - **(ii) Non-subfield frame-less** configurations (rankdegen's finite-volume `(1,s,t,t²)`/`(1,s+23,t,t²)` example;
     decoupled fibre products) — no `O_k`-action, no relative resultant to charge.
   - **(iii) Large boxes beyond `B^{1/(r−1)}`** (LP_exc / rank3 §7(b)): small-|F| lattice points in excess of volume.

**Probability assessment.** < 2% that the G_r case yields sub-L[1/3]. The single most promising remaining crack is now
sharper — door (i): a **relative-rank-≥3 type-D family over a fixed small subfield with non-power relative
generators**, which is exactly a "relative R1" problem over `k`, blocked by the infinite volume (acting units). It is
the clean, concrete target: prove (or refute) that the `m−1` relative linking conditions force `v_𝔭(Δ^{rel}) ≥
P(m, D/e)` over `O_k`. No construction beating L[1/3] is known for any door (codex: "none of the constructions here
meets that target").

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

**Proposition RR.** Let the frame-less degeneracy be a module over a common subfield `k = Q(w)`, `[k:Q] = e ≥ 2`,
realised on `K = k(s)`, `K' = k(s')` with relative minimal polynomials `f, g ∈ O_k[x]`, monic (or with a unit linked
leading coordinate) of equal relative degree `n = n' = d/e`, and linking maps whose common reduction `ψ: O_k → Z/N`
(surjective, `1 ↦ 1`) sends `f, g` to polynomials over `Z/N` sharing the linking root `v`. Let `R = Res_k(f, g) ∈ O_k`,
`R ≠ 0`.
- (RR1) [**proved**] By the adjugate (Sylvester) identity `R = Af + Bg` in `O_k[x]`, applying `ψ` and evaluating at
  `v` gives `ψ(R) = 0` (so no leading-coefficient hypothesis is needed; `ψf, ψg` share root `v`). Thus
  `R ∈ ker ψ = RO_k ⊆ ker ψ`, an ideal of index `N` (any unital map onto `Z/N` is surjective), and
  `N = [O_k : ker ψ] | [O_k : RO_k] = |Norm_{k/Q}(R)|` (the index is `|det(×R)|`; needs only `R ≠ 0`). So
  **`N | Norm_{k/Q}(R)`**. [`relres.psi_resultant_is_zero`, `rr_quadratic`, `rr_over_Qi` (genuine `O_k`).]
- (RR2) [**proved**] By the resultant product formula, `|Norm_{k/Q}(R)| ≤ exp(n'h₁ + n h₂ + e·n·n'·log 2)` with
  `h_i` the **projective (Mahler) height** of the relative polynomial, summed over all `e` embeddings of `k` and
  including leading coefficients. Hence `h₁ + h₂ ≥ x/n − d·log 2 = e·x/d − O(dn)`.
- (RR3) [**heuristic / model; relative rank 2 (`r = 2e`) or relative power basis only**] For such a module the
  N-saving exponent is `≥ e·x/d = r·x/D` (exactly, since `r = 2e`, `D = 2d`), i.e. `κ_gen`-exponent `c = 1`, giving
  **α ≥ 1/3** in the box ledger (§1), with the same typical-size caveat. *This is a G-height statement; turning it
  into a box small-value count (G-count) is the model step, weakest here because type D is the `V_F = ∞` case.* For
  relative rank `m ≥ 3` with non-power generators, typical `log|F|` is the *max* (not sum) of generator heights, so RR
  gives only `(r/m)·x/d` (`c = 0`, L[1/4]); that sub-case is **open** (door (i), §0.5).
- (RR4) [**proved**] If `R = 0` with `f, g` irreducible over `k`, then `f = g` and there is a `k`-isomorphism
  `ι: K → K'` with `s ↦ s'`; the maps agree on `O_k` and on `s`, so `φ' ∘ ι = φ` — this is type A, and by rankdegen
  Lemma CC (indices coprime to `N`) a two-sided square gives `φ(γ) = ±φ'(ι γ)`, a trivial congruence. Useless.

**Verification (`relres.py`).**
- codex's `N=65`, `k=Q(i)`, `f=x²−14`, `g=x²−79`: `R = 65²`, `Norm = 65⁴ = 65·274625`, `N | Norm` ✓; `ψ(R)=0 mod 65` ✓.
  A genuine `k=Q(i)` case with non-rational `R = 436−152i`, `Norm = 213200`, `65 | Norm`, `ψ(R)=0` (`rr_over_Qi`,
  `test_rr_genuine_Qi_example_Ok_does_work`) confirms RR1 when `O_k` actually does work.
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

**Why RR subsumes SL (for relative rank 2 / power bases).** The clean statement is RR itself: total projective height
`h₁ + h₂ < x/n − d·log 2` forces `Norm_{k/Q}(R) < N`, so by RR1 `R = 0` (RR4, useless). (The naive "`N^{o(1)}`
generators ⇒ `Norm(R) < N`" is not quite right at growing degree — `Norm(R)` can reach `N` once the per-conjugate
height `ε ≥ 1/(2dn)` — but that `ε → 0` in the L[1/4] regime, and the correct threshold is the height one above.) So
genuine small linking is excluded by height, with no appeal to the hardness of root extraction. **This subsumption is
limited to relative rank 2 / power bases;** for relative rank `m ≥ 3` (door (i)) RR does not exclude it and SL2
(conjectural) would be needed.

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

**Opus re-vet of Proposition RR (second critic pass; all applied).** Verified RR1 (correct; needs only `R ≠ 0`; uses
the adjugate identity, no leading-coefficient hypothesis), RR2 (error term `O(dn)` not `O(d)`; `h_i` must be the
projective Mahler height over all embeddings incl. leading coefficients), RR4 (fixed: the collapse is via the
`k`-isomorphism `ι`/type A, not "square already in `k`"). **Forced downgrades:** (a) RR3 is a G-height result, not a
box count, so it is **heuristic/model**, not proved; (b) RR3 covers only **relative rank 2 / power bases** — a type-D
module over a fixed small `k` with relative rank `m ≥ 3` and non-power generators gets only `c = 0` (L[1/4]) from RR,
a **new third open door** (door (i)), closable only by an unproved *relative* R1 or the conjectural SL2. Also: SL0's
`eigenvector_lemma` is correct and rigorous (minimal-polynomial annihilation, non-circular) but elementary. Confirmed
the genuine `Q(i)` example (`O_k` does real work). All folded into §0.4–5 and §3.

## 8. GPT-6.1-sol review round (`codex/grdeg_2.out`) — converges with the Opus re-vet

Codex independently reached the same verdict: **RR1 is rigorous** (it gives the full index identity
`|Norm(R)| = [O_k:ker ψ][ker ψ:(R)] = N·[ker ψ:(R)]`, divisibility by all of `N` including prime powers, no
equidistribution assumption); the **categorical "No" is not a theorem**, "no construction found" is correct; and the
**main omission is growing relative rank** (door (i)). Additional fixes applied:
- RR2: the pre-division error is `O(enn')`, not `O(d)`; codex's witness `f=(T−1)ⁿ−2`, `g=(T+1)ⁿ−2` shows the
  resultant-norm vs height-sum gap is *quadratic* in `n`, so the Mahler-height definition (not coefficient 2-norm)
  must be used. The final bound `h₁+h₂ ≥ e·x/d − d·log 2` stands. Typo fixed: `Norm(65-example) = 65⁴`.
- RR3: the N-saving scale is `e/d = r/D` (the first draft's extra factor 2 removed).
- RR4: "a two-sided square is already a square in `k`" is **false** (codex: `14` is a square in `K=k(√14)` but not
  `k=Q(i)`); the correct argument is the `k`-isomorphism `ι` + agreeing reduction maps (one global sign), as now
  written. Codex also shows the order-index hypotheses of Lemma CC are indispensable (an `s=65√14`, split-map example
  yields a *precomputable* `gcd=13` — a free split, not an algorithm).
- **Concrete door-(i) witness (codex).** A non-binary subfield frame-less family: `N=65`, `k=Q(i)`, `K=k(s)` with
  `s³=14`, `K'=k(s')` with `s'³=79`, module `(1,i,s,is,s²,is²)` vs `(1,i,s',is',s'²,is'²)`, `r=6`, `D=12`, distinct
  norm forms `A³+14B³+196C³−42ABC` / `…79…`, genuinely frame-less (two complementary 3-blocks, no 7-row frame). RR1
  still holds (`R=−65³`, `Norm=65⁶`), but a two-generator resultant does **not** manufacture the missing factor `m=3`:
  RR charges `e`, not `r`. This is the concrete shape of the open crack.
- Scope: codex notes an N-dependent subfield is *not* intrinsically outside RR1 (RR1 does not need `k` fixed); what is
  extra is the base-field presentation/height ledger. Tower-NFS lower bounds are a finite-field-DLP literature and do
  not themselves bound integer factoring; door (i) is genuinely unresolved.

**Net effect of both reviews:** RR1 (divisibility) is a clean proved obstruction; RR2 (height) is proved with the
Mahler-height definition; RR3 (⇒ α≥1/3) is model/heuristic and only for relative rank 2 / power bases; RR4 is proved
with the `ι` argument and the order-index hypothesis. The subfield frame-less case is **not** fully closed — door (i)
(growing relative rank over a fixed small subfield) is the sharp open problem, a "relative R1 over `k`".

## 9. Status summary (labels after both reviews)

| claim | label |
|---|---|
| Ledger exponent arithmetic (α=1/(4−c), endpoints 1/4, 1/3) | proved (arithmetic); endpoints heuristic (outside H5_r) |
| RR1 (`N | Norm_{k/Q}(R)`) | **proved** |
| RR2 (height `h₁+h₂ ≥ e·x/d − O(dn)`) | **proved** (Mahler height) |
| RR3 (⇒ α ≥ 1/3) | heuristic/model, **relative rank 2 / power basis only** |
| RR4 (`R=0` ⇒ useless) | **proved** (ι + Lemma CC order-index hyp.) |
| SL0 eigenvector lemma | proved (elementary) |
| SL2 (one fixed √D₀ mod N factoring-hard) | conjectured, **not** proved, and not needed |
| E1 sublattice confinement | experimentally supported (sampled boxes, toy N) |
| door (i): relative rank ≥ 3 type-D over fixed small `k` | **OPEN** (would need a relative R1 or SL2) |
| doors (ii) non-subfield frame-less, (iii) large-box LP_exc | OPEN (unchanged from rankdegen) |
