# factoring-connection

AI-assisted research on the asymptotic complexity of integer factorization. Claude agents wrote the
code, proofs and notes; GPT-6.1-sol (Codex CLI) did literature sweeps and adversarial reviews.

> **Status of every result below: proved internally and audited by blind AI referees, but NOT
> refereed by a human.** The subexponential results share one chain of internal lemmas. Treat them
> as claims, and the ideas as the main contribution, until an expert has checked that chain.

Notation: L_N[α, c] = exp((c + o(1))(log N)^α (log log N)^{1−α}).

## Context

- **Heuristic:** the number field sieve runs in L[1/3, 1.923] (MNFS 1.902). Its running time rests on
  unproven smoothness assumptions.
- **Rigorous record:** L[1/2, 1] (Lenstra–Pomerance 1992). Lee–Venkatesan (2018) rigorously find a
  congruence of squares in L[1/3], but cannot prove it is non-trivial ("fruitfulness"), so their result
  does not give rigorous factoring.
- **Deterministic record:** N^{1/5} (Harvey 2021; Harvey–Hittmeir 2022).

## Top findings

**Rigorous randomized factoring** (Las Vegas, complete factorization of every N).

| Result | Bound | Audit | Where |
|---|---|---|---|
| **Theorem A** | L_N[1/2, c] for every c > 1/√2 (0.919 at degree 9) | 5 blind referees: accept with fixes | [`paper2/`](paper2/main.tex), [`thmA_proof`](notes/thmA_proof.md), [`thmA_audit`](notes/thmA_audit.md) |
| **Theorem RA + Cor. P** | L_N[1/2, o(1)] (recursive shape test) | 3 blind referees: accept with fixes | [`recurse`](notes/recurse_results.md), [`auditRA`](notes/auditRA_audit.md) |
| **Theorem U** | **L_N[2/5, c]** for a finite c | 4 blind referees: accept with fixes | [`lv5unif`](notes/lv5unif_results.md), [`auditU`](notes/auditU_audit.md) |
| **Theorem 2L-DC** (conditional) | **L_N[α] for every α > 3/8**, assuming the algebraic block Theorem AB (written as a sketch); **L_N[3/8]** assuming AB₂. The rational blocks are proved (Theorem RAT) | GPT-6.1-sol + 2 Opus critics | [`multidescent`](notes/multidescent_results.md), [`uprime`](notes/uprime_results.md), [`pc3`](notes/pc3_results.md), [`c12rat`](notes/c12rat_results.md), [`c12alg`](notes/c12alg_results.md) |
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
Lemme 5, whose proof cites Harper's unpublished preprint.

## What blocks further progress (summary)

Full list, phase by phase: [`BARRIERS.md`](BARRIERS.md). These are model theorems, not lower bounds.

- **Rigorous below 2/5.** One-level designs stop at 2/5 (Prop R2). The two-level route now lacks only the algebraic block AB.
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

1. **Write out Theorem AB ⇒ rigorous L[α] for every α > 3/8.** Everything else in the two-level chain is proved
   internally: U′ ([`uprime`](notes/uprime_results.md)), decoupling ([`pc3`](notes/pc3_results.md)), the rational
   blocks RAT ([`c12rat`](notes/c12rat_results.md)), and the DW design, two-point smoothing and child bound
   ([`c12alg`](notes/c12alg_results.md)). AB₂, with its obligations O1–O3 and O5, would give L[3/8] itself.
2. **NFS fruitfulness ⇒ rigorous L[1/3, 1.923].** Prove Lemma M or Conjecture F ([`lemmaM`](notes/lemmaM_results.md),
   [`nfsfruit`](notes/nfsfruit_results.md)): Lee–Venkatesan's randomized NFS congruences are non-trivial with probability
   L[1/3, −o(1)]. Two phase-28 surveys independently ranked it first; feasibility with current tools is low.
3. **Three-level descent (4/11),** which is borderline by the same route ([`pc3`](notes/pc3_results.md), [`multidescent`](notes/multidescent_results.md) §3.5).
4. **Human refereeing** of the relation-supply chain: [`paper2/`](paper2/main.tex) §§5–6, Appendix A, and [`lv5unif`](notes/lv5unif_results.md).
5. **Deterministic worst case below N^{1/5}:** coverage with λ < J^{1/4} ([`coverage38`](notes/coverage38_results.md)).

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
