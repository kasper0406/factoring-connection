# factoring-connection

AI-assisted research on the asymptotic complexity of integer factorization. Claude agents wrote the
code, proofs and notes; GPT-6.1-sol (Codex CLI) did literature sweeps and adversarial reviews.

> **Status:** the earlier main theorems are internal claims audited by AI referees, not human-refereed.
> The subexponential results share one chain of internal lemmas. The newer research-frontier notes
> explicitly distinguish reconstructed proofs, reported results, conditional reductions and open tasks.
> Treat these claims accordingly until an expert has checked them.

Notation: L_N[α, c] = exp((c + o(1))(log N)^α (log log N)^{1−α}).

## Context

- **Heuristic:** the number field sieve runs in L[1/3, 1.923] (MNFS 1.902). Its running time rests on
  unproven smoothness assumptions.
- **Rigorous record:** L[1/2, 1] (Lenstra–Pomerance 1992). Lee–Venkatesan (2018) rigorously find a
  congruence of squares in L[1/3], but cannot prove it is non-trivial ("fruitfulness"), so their result
  does not give rigorous factoring.
- **Deterministic record:** N^{1/5} (Harvey 2021; Harvey–Hittmeir 2022).

## Current research frontier (2026-11-26)

The latest non-sieve work connects multiplication-network counts, structured
square-weight collisions and quadratic-form minima. These are **partial results
and conditional reductions, not a new factoring algorithm**. In particular,
N^{1/8} is a conditional almost-all balanced-semiprime bound, not a replacement
for the existing N^{1/6} result or an improvement over heuristic L[1/3].

| Direction | Latest result / missing step | Details |
|---|---|---|
| **Canonical index (primary)** | One hidden quadratic-form minimum per prime/root branch; reported band coverage would yield N^{1/8}·polylog if that minimum were efficiently computable. Need independent coverage audit and native evaluator. | [Shared collision and canonical index](notes/shared_collision_nov26.md) |
| **Public ellipse count (primary)** | Modular lattice-point count mod 16 recovers the hidden squared ideal class at bounded cutoffs; exact principal contribution and class-square lifting required. No fast count known. | [Ellipse reduction and ambiguity](notes/shared_collision_nov26.md#4-public-ellipse-count-and-squared-class-decoding) |
| **Whole-disk parity** | Parity of the positive disk with up to two y-bit filters suffices to extract a factor. A Boolean decision diagram handles all norm levels simultaneously; intermediate carry diagrams still dominate. | [Parity reduction and benchmark](experiments/ellipse_disk_mod8/README.md), [which translation identities already expose a factor](experiments/ellipse_disk_mod8/translation.md) |
| **Centered form recurrence** | Ordinary signed counts mod 16 on explicitly changed forms suffice for factoring in O(log N) queries. The native traversal still enumerates primitive directions. A 35-bit control has one mixed hit among 8,192 inequivalent locally soluble forms. | [Recurrence, scope and checks](experiments/ellipse_disk_mod8/centered.md) |
| **Norm-N shell, all slopes** | A rational bisector and batched linear products distinguish the 8,192-slope fixtures. For Δ of size N^{1/2+o(1)}, B=1 costs N^{1/8+o(1)} including discriminant preparation. A counting bound rules out almost-all coverage by shallow ℓ=1 prefixes. | [Algorithm and coverage proof](experiments/ellipse_disk_mod8/bisector.md) |
| **Whole-prefix factor decision** | A public list of classes whose squares have a reduced minimum in the desired interval suffices to find a divisor or certify the interval empty. The current list builder scans the proposed minima. A batched verifier works for square and nonsquare minima. | [Decision, candidate budget and local obstruction](experiments/ellipse_disk_mod8/class_prefix.md) |
| **Selecting the prime class** | On the near-square ℓ=1 branch, any selector confined to reduced prime-class minima below N^{1/4−ε} has vanishing coverage. A rational-approximation shortlist found 14/96 fresh inputs. The shortlist misses typical-sized minima. | [Proof, blind benchmark and cutoff fixture](experiments/ellipse_disk_mod8/class_selection.md) |
| **Rectangles of norm levels and rows** | Discriminant projection and cancellation across square multiples reproduce the 47-bit residues 0 and 4 and recover a factor. The remaining parity work is concentrated on primitive nonsquare levels above B/4. The seed scan already costs at least B. | [Algorithm, cost and checked cutoffs](experiments/ellipse_disk_mod8/rectangle_rays.md) |
| **Primitive nonsquare band** | Below √Δ/2 each mixed lattice has at most one direction up to sign and scaling. Grouping possible levels in a smaller discriminant and verifying targeted roots reproduced the 581/582 transition. Two prime-discriminant controls isolate the remaining bulk-summation problem. | [Proof, comparisons and prime-discriminant fixtures](experiments/ellipse_disk_mod8/conductor_band.md) |
| **Weighted ellipse moment** | A branch with ℓ=1 is constructed from N alone. O(log N) weighted queries recover a factor. Two exact evaluators factor all 560 small balanced test pairs; both still take power-of-N work. | [End-to-end CLI and measured costs](experiments/ellipse_aggregate/README.md), [moment proof](experiments/ellipse_row_batches/README.md) |
| **Multiplication network (secondary)** | Normalized divisor-prefix count mod 3 suffices on semiprimes, but branch duplicates and intermediate state size must be handled. Cut-rank obstruction; masked-product local counts; reported G^{3/4+o(1)} short-list collision evaluator needs proof/cost audit. | [Network and grid](notes/shared_collision_nov26.md) |
| **QRP / hidden conductor (parked)** | Local-character oracle, conductor wheels, class-number and Jacobi-prefix identities; variable-endpoint parity factors mixed semiprimes. No efficient conductor decoder or aggregate known. | [QRP and interval targets](notes/qrp_conductor_nov26.md) |
| **Smooth decoding / binomial midpoint (paused)** | Tested jump and local-state mechanisms stalled; unrestricted adaptive decoding and a native single-index evaluator remain open. | [Priorities and negative results](notes/research_roadmap_nov26.md) |

Start here: [prioritized research roadmap](notes/research_roadmap_nov26.md) and
[experiment provenance and checks](experiments/nov26_consolidation/README.md).
The fresh exact checks include 12,697 reduced-form decodings, 700 grid-block
checks and 40 direct ellipse queries on the two ambiguity examples.
Original late-stage scratch files were not recovered; the new notes mark
reported asymptotic results separately from reconstructed proofs and checks.
The [global ellipse experiments](experiments/ellipse_global/README.md) give the finite
Fourier identity, the clipped-row reduction and exact class lifting with factored
discriminant. The [row experiments](experiments/ellipse_row_batches/README.md)
include 1,586 independently checked shells and 6,188 interval checks on the two
preserved collision inputs and a replacement cohort of 192 balanced semiprimes.
Neither these results nor the earlier model barriers settle general factoring.

## Top findings

**Rigorous randomized factoring** (Las Vegas, complete factorization of every N).

| Result | Bound | Audit | Where |
|---|---|---|---|
| **Theorem A** | L_N[1/2, c] for every c > 1/√2 (0.919 at degree 9) | 5 blind referees: accept with fixes | [`paper2/`](paper2/main.tex), [`thmA_proof`](notes/thmA_proof.md), [`thmA_audit`](notes/thmA_audit.md) |
| **Theorem RA + Cor. P** | L_N[1/2, o(1)] (recursive shape test) | 3 blind referees: accept with fixes | [`recurse`](notes/recurse_results.md), [`auditRA`](notes/auditRA_audit.md) |
| **Theorem U** | **L_N[2/5, c]** for a finite c | 4 blind referees: accept with fixes | [`lv5unif`](notes/lv5unif_results.md), [`auditU`](notes/auditU_audit.md) |
| **Theorem U-MF** (new) | **L_N[2/5, c]**, c ≈ 9.45 at C_γ = 0; self-contained, and **no unpublished input** (integer multiplier plus Design MF remove every BV-type estimate, Drappeau Lemme 5 and Harper's preprint) | GPT-6.1-sol adversarial: accept with minor fixes; Opus proof critic: accept with fixes (~80% each); literature critic applied (extends Lee–Venkatesan 2018 §6) | [`mfU_theorem`](notes/mfU_theorem.md), [`mfU`](notes/mfU_results.md) |
| **Theorem 2L-final** | **L_N[α] for every α ∈ (3/8, 2/5)**, and **L_N[3/8]** under a design constraint (two-level descent; Theorems AB/AB₂, RAT, U′, DC) | Blind audit [`audit2L`](notes/audit2L_audit.md): **accept with fixes** (window ~70%, 3/8 ~55%). The audit found and repaired a fatal-if-unfixed gap (RAT-2 parent filter, Lemma SPR-P′); one write-up gap is open (U′ under DW asserted). **Theorem 2L-MF** ([`mf2L`](notes/mf2L_results.md)) removes Drappeau Lemme 5 / Harper's preprint, using Harper Compositio 2016 Prop 2 and Drappeau CJM 2015 Props 2.9/2.13 (codex re-check 92%, relative to the internal chain) | [`thmAB`](notes/thmAB_results.md), [`c12rat`](notes/c12rat_results.md), [`c12alg`](notes/c12alg_results.md), [`pc3`](notes/pc3_results.md), [`uprime`](notes/uprime_results.md) |
| **Theorem DL** | discrete logs in F_p^*: L_p[1/2, c] for every c > 1/√2 (previous record L[1/2, √2], Pomerance 1987) | GPT-6.1-sol + 2 Opus critics | [`dlog`](notes/dlog_results.md) |
| **Theorem DL25** | discrete logs in F_p^*: **L_p[2/5, c]** | GPT-6.1-sol + 2 Opus critics | [`dlog25`](notes/dlog25_results.md) |

**Deterministic factoring** (exponential regime).

| Result | Bound | Where |
|---|---|---|
| **Theorem AV** | always factors, **average N^{1/6}·polylog** over balanced semiprimes, worst case N^{1/5} | [`avgtime`](notes/avgtime_results.md) |
| **Theorem D′** (+ 4th moment) | N^{1/6+5ε/12}, factors all but N^{−ε} of balanced semiprimes | [`dcltheorem`](notes/dcltheorem_results.md), [`sieve2d`](notes/sieve2d_results.md), [`paper/`](paper/main.tex) |
| **Theorem R** | the worst case reduces to a covering statement; coverage exponent c > 3/8 would give worst case below 1/5 | [`dchworst`](notes/dchworst_results.md) |

**Other results worth reusing:**
- #E(Z/N) for one fixed curve factors every good N (Theorem EXT, [`tracedata`](notes/tracedata_results.md)).
- Point counts of modular curves factor N: #X_0(N)(F_{ℓ²}) and #X_1(N)(F_ℓ) ([`geometry`](notes/geometry_results.md)).
- Dedekind sums and Salié sums mod N are factoring-equivalent ([`newfamily`](notes/newfamily_results.md)).
- Addition in the prime-exponent representation ⇔ factoring ([`factorbits`](notes/factorbits_results.md)).

## Key ideas (for anyone picking this up)

1. **Dixon descent removes fruitfulness.** Descending a *random square* rather than a fixed target makes
   the congruence non-trivial with probability ≥ 1/2. This is the step Lee–Venkatesan lacked; it is what
   costs the exponent (2/5 instead of 1/3).
2. **Relation supply from moments.** First and second moments come from smooth-modulus Bombieri–Vinogradov
   (Drappeau), divisor switching plus the large sieve, and an N·w multiplier that spreads forced
   residues. Batch smoothness replaces HECM, and recursion factors the split halves.
3. **Explicit-error, degree-uniform estimates** take one-level descent from 1/2 down to 2/5. The floor is
   the balance between split acceptance (log N)^{−2k} and norm exponent c/√k.
4. **Multi-level descent by union bound** (Prop UB, Prop DC). A child dies with small probability if
   pair-only conditional variances are small. This avoids joint control of one polynomial at ~d^{3/2}
   points, which earlier work (`nested`, `escapeR2`) took to be necessary. j levels give exponent
   (j+1)/(3j+2).
5. **Hybrid restricted BV** (Lemma 8.5H): B-smooth moduli, B₂-smooth values, normalised by the density of
   the actual modulus family.

Weakest link: the relation-supply chain (LD*, L2_fam, Lemma 8.5′). Lemma 8.5′ uses Drappeau 2015,
Lemme 5, whose proof cites Harper's unpublished preprint. **Theorem U-MF avoids all of this for L[2/5].** With an
integer multiplier, F(P₀) mod G(P₀) is uniform, so smoothness is an interval density (Lemma OP). Its external inputs
are all published, numbered theorems ([`mfU_theorem`](notes/mfU_theorem.md) is written for a human referee).
Theorem 2L-MF does the same for L[3/8], but through different published inputs. Only Theorems A/RA still use
Lemme 5.

## What blocks further progress (summary)

Full list, phase by phase: [`BARRIERS.md`](BARRIERS.md). These are model theorems, not lower bounds.

- **Rigorous below 2/5.** One-level designs stop at 2/5 (Prop R2). The two-level route reaches 3/8 (Theorem 2L-final, audited: accept with fixes); 4/11 is borderline.
  Depth ≥ 4 fails because the required precision falls below the natural covariance. Rigorous L[1/3]
  needs ERH plus a far stronger all-target supply (`nested` Theorem NR).
- **Heuristic below L[1/3] or MNFS 1.902.**
  - Resultant and area theorems: N² | disc.
  - Size–entropy law: every direct target is poly-time factoring in disguise ([`suffmap`](notes/suffmap_results.md)).
  - Constructing smooth numbers near kN costs L[2/3] per output for every constructor class tried
    (decoders, lattices, adaptive CRT, navigation).
- **Leaks without smoothness.** Every known way to make p and q behave differently is a group order
  (smoothness), an archimedean approximation (exponential) or a quadratic genus symbol.
  - Proved in formal models: pure algebra (Theorem PT), N-dependent exponents, AC⁰ on the bits of N,
    reciprocity invariants.
  - Searched without success: class groups, continued fractions, geometry, elliptic-curve traces.
  - Open only in narrow corners: deterministic depth-≥2 exponent gates; single-curve traces a_N.
- **Unexpected connections (phase 28).** A survey of about 300 connections across 40 fields (Langlands, K-theory,
  harmonic analysis, random matrices, knots, homotopy theory, tilings, operator algebras, ML, …) found none that gives
  faster factoring. They fail by seven recurring mechanisms ([`connections_results`](notes/connections_results.md)).
- **Factoring-equivalent invariants as new views (phase 29).** 55 invariants (φ, σ, class counts, K-groups, modular
  coefficients, class numbers, …): no native algorithm beats factoring. Native costs are N^{1/3}–N^{1/2} for divisor
  types and N^{1/5} for analytic ones; proved barriers hold only in restricted models (Prop TRI, LOC, Theorem LIN).
  The cleanest open view: one coefficient a_f(N) of a fixed non-CM cusp form ([`native_results`](notes/native_results.md)).
- **Deterministic.** Subexponential under ERH reduces to one constructive lemma, DSR, and every black-box
  derandomiser fails on it. The worst case below N^{1/5} needs coverage with λ < J^{1/4}.

## Open problems, ranked

1. **Close the audit's open item for 2L-final:** write out Hyp U′ under design DW; its leading-order checks pass
   ([`audit2L`](notes/audit2L_audit.md)).
2. **Merge the audit repair into one clean write-up:** 2L-final + audit2L's Lemma SPR-P′ + mf2L's published-input route,
   as a referee-ready note like [`mfU_theorem`](notes/mfU_theorem.md).
3. **Human refereeing of [`mfU_theorem`](notes/mfU_theorem.md).** It is the cleanest, self-contained candidate.
4. **NFS fruitfulness ⇒ rigorous L[1/3, 1.923].** Mixing descended rows into the matrix gives fruitfulness, but one such
   row costs as much as an L[1/3] descent ([`fruit`](notes/fruit_results.md), Theorem HYB). Sufficient: Conjecture F, or
   a rigorous L[1/3] single descent plus multi-row absorption.
5. **Three-level descent (4/11),** which is borderline by the same route ([`pc3`](notes/pc3_results.md) §7.5).
6. **Deterministic worst case below N^{1/5}.** Phase-discarding estimates provably miss by a saving of 7/25 vs 4/25
   ([`cov38b`](notes/cov38b_results.md) Theorem B). This needs phase-retaining multilinear estimates.

## Repository map

| Path | Contents |
|---|---|
| `paper2/` | Private draft preprint: Theorem A, RA/P, DL (54 pp, `tectonic paper2/main.tex`) |
| `paper/` | Draft paper: Theorem AV, D′, R, T1 (exponential regime) |
| `notes/` | One `*_results.md` per work item (verdict first), pre-registrations `*_prereg.md`, audits `*_audit.md` |
| `BARRIERS.md` | Detailed barrier map, phase by phase |
| `STATUS.md` | Chronological log of all phases and agent ownership |
| `APPROACHES.md` | Catalogue of everything tried, with verdicts |
| `src/`, `tests/` | Python code per work item; pytest suite |
| `codex/` | GPT-6.1-sol prompts and outputs; `codex/run.sh` wrapper |
| `lean/`, `leanA/` | Lean 4: D′ lemmas; conditional Theorems A/RA/P. Cost assembly and glue are proved; FM, LD* and 4 published theorems are hypotheses; checked by Comparator --paranoid and 6 kernels. See [`leanA_results`](notes/leanA_results.md) |
| `experiments/` | Configs, raw JSONL, summaries |

**Rules the agents followed:**
- exact arithmetic;
- operational code sees only N;
- a status label on every claim;
- pre-registration before experiments;
- every theorem reviewed by GPT-6.1-sol and an independent Opus critic.

## Running

```bash
python3.12 -m venv .venv && .venv/bin/pip install sympy gmpy2 python-flint numpy pytest
.venv/bin/python -m pytest -q          # full suite
```

Every results note ends with the commands that reproduce its experiments. The git history was
rewritten once to drop a 110 MB file; old pre-registration SHAs map via
[`COMMIT_SHA_MAP.md`](COMMIT_SHA_MAP.md).
