# factoring-connection

AI-assisted research on the asymptotic complexity of integer factorization. Claude agents wrote the
code, proofs and notes; GPT-6.1-sol (Codex CLI) did literature sweeps and adversarial reviews.

> **Status of every result below: proved internally and audited by blind AI referees, but NOT
> refereed by a human.** The subexponential results share one chain of internal lemmas. Treat them
> as claims until an expert has checked that chain.

## Top findings

Notation: L_N[α, c] = exp((c + o(1))(log N)^α (log log N)^{1−α}).

**Rigorous randomized factoring** (Las Vegas, complete factorization of every N). The previous rigorous
record is L[1/2, 1] (Lenstra–Pomerance 1992).

| Result | Bound | Audit | Where |
|---|---|---|---|
| **Theorem A** | L_N[1/2, c] for every c > 1/√2 (0.919 at degree 9) | 5 blind referees: accept with fixes | [`paper2/`](paper2/main.tex), [`thmA_proof`](notes/thmA_proof.md), [`thmA_audit`](notes/thmA_audit.md) |
| **Theorem RA + Cor. P** | L_N[1/2, o(1)] (recursive shape test) | 3 blind referees: accept with fixes | [`recurse`](notes/recurse_results.md), [`auditRA`](notes/auditRA_audit.md) |
| **Theorem U (lv5unif)** | **L_N[2/5, c]** for a finite c | 4 blind referees: accept with fixes | [`lv5unif`](notes/lv5unif_results.md), [`auditU`](notes/auditU_audit.md) |
| **Theorem DL** | discrete logs in F_p^*: L_p[1/2, c] for every c > 1/√2 (previous record L[1/2, √2], Pomerance 1987) | GPT-6.1-sol + 2 Opus critics | [`dlog`](notes/dlog_results.md) |
| **Theorem DL25** | discrete logs in F_p^*: **L_p[2/5, c]** (same constant as Theorem U) | GPT-6.1-sol + 2 Opus critics | [`dlog25`](notes/dlog25_results.md) |

The key ideas:
- **Dixon descent** of a random square removes NFS fruitfulness.
- **Relation supply** is proved by first and second moments: smooth-modulus Bombieri–Vinogradov, divisor switching plus the large sieve, and an N·w multiplier that spreads forced residues.
- **Batch smoothness** replaces HECM.
- **Recursion** factors the split halves using the algorithm itself.
- **Explicit-error, degree-uniform estimates** reach exponent 2/5.

Weakest link: the relation-supply chain, i.e. LD*, L2_fam and Lemma 8.5′. Lemma 8.5′ uses Drappeau 2015, Lemme 5, whose proof cites Harper's unpublished preprint.

**Deterministic factoring** (exponential regime; the worst-case record is N^{1/5}, Harvey–Hittmeir 2022).

| Result | Bound | Where |
|---|---|---|
| **Theorem AV** | always factors, **average N^{1/6}·polylog** over balanced semiprimes, worst case N^{1/5} | [`avgtime`](notes/avgtime_results.md) |
| **Theorem D′** (+ 4th moment) | N^{1/6+5ε/12}, factors all but N^{−ε} of balanced semiprimes | [`dcltheorem`](notes/dcltheorem_results.md), [`sieve2d`](notes/sieve2d_results.md), [`paper/`](paper/main.tex) |
| **Theorem R** | the worst case reduces to a covering statement; coverage exponent c > 3/8 would give worst case below 1/5 | [`dchworst`](notes/dchworst_results.md) |

**Barriers** (model theorems, not factoring lower bounds):
- T1/T2: local designs need N^{1/6}.
- Coverage cost ≥ N^{1/5}(λ⁴/J)^{1/5}.
- Ratio gap peaking at 1/5.
- Below 2/5 is blocked for every route tried:
  - nested descent (Prop R2);
  - exporting algebraic children ([`escapeR2`](notes/escapeR2_results.md) Cor EX);
  - averaging over N, which adds no parameters ([`almostall13`](notes/almostall13_results.md) Prop NB);
  - special forms ([`snfsrig`](notes/snfsrig_results.md)).
- Heuristic NFS: nothing beats MNFS 1.902 or exponent 1/3 ([`heurnfs`](notes/heurnfs_results.md)).
- Deterministic under GRH: fruitfulness is not the obstruction, deterministic relation supply is ([`detGRH`](notes/detGRH_results.md)).
- Inverse-sieve containers for the deterministic worst case fail: the information is spread too thin ([`inversesieve`](notes/inversesieve_results.md)).
- Below L[1/3]: no mechanism was found on any of the three routes tried:
  - integer "systematic equations" ([`systeq`](notes/systeq_results.md)); the exact PTE version is refuted by [`ptebarrier`](notes/ptebarrier_results.md) Lemma R / Theorem W, and the joint verdict is in [`ptesplit`](notes/ptesplit_results.md);
  - beating the resultant bound ([`resbarrier`](notes/resbarrier_results.md): N² | disc, area theorem);
  - non-smoothness invariants ([`nonsmooth`](notes/nonsmooth_results.md): FROB, OA). Open remainders, none with a candidate: non-exact compressions; pencils of several presentations; parameter corners; and resbarrier's rank ≥ 3 / no-rational-side case, closed only heuristically (floor L[1/4] even if that heuristic fails).
- Output-sensitive smooth-pair enumeration: no gain, MNFS stays at 1.902 ([`smoothenum`](notes/smoothenum_results.md)). The heuristic model was validated on real CADO-NFS runs ([`cadobench`](notes/cadobench_results.md)).
- Non-relational paradigms tried: analytic sums ([`analytic`](notes/analytic_results.md)), short circuits ([`slpcirc`](notes/slpcirc_results.md), Theorem SP), and the floor-sum ladder ([`floorsum`](notes/floorsum_results.md), [`floorhard`](notes/floorhard_results.md)). The ladder at degree j ≥ 3 would require factoring from fewer of p's top bits than Coppersmith needs, a known barrier. The general oracle is NP/#P-hard; top blocks stay open.
- Classical "Fourier side" ([`periodfind`](notes/periodfind_results.md)): for an observable g, g(a^x) has period-concentrated Fourier mass only if g predicts a large-order multiplicative character (Theorem T). Additive, ring-op and computable-character observables are provably flat. Smoothness observables are the only class with real peaks, which leads back to relation methods.
- Phase 18, where the limit bites for each open remainder:
  - algebraic preprocessing + readout: it would also solve prime-field discrete logs ([`prepread`](notes/prepread_results.md));
  - floor-sum ladder: exactly Coppersmith's known-bits barrier ([`topblock`](notes/topblock_results.md));
  - degenerate / large-box rank-r: area theorems plus one counting hypothesis ([`rankdegen`](notes/rankdegen_results.md));
  - product listing: Johnson-radius decoders gain only L^{o(1)} ([`pilist`](notes/pilist_results.md) Prop JR).
- Phase 19, the sufficient side: every target that would give α < 1/3 directly is as strong as polynomial-time factoring in disguise ([`suffmap`](notes/suffmap_results.md); size–entropy law, Theorem SS). The best realistic target improves only the constant: constructing L[1/3]-smooth numbers in short intervals near multiples of N would give L[1/3, 1.21] (Prop SI). S-units and abc add no barrier and no new target ([`abcsunit`](notes/abcsunit_results.md)). Lattice CVP helps only at approximation factor γ → 1, and Schnorr 2021 is quantitatively refuted ([`schnorrcvp`](notes/schnorrcvp_results.md)).
- Phase 20, Prop SI (constructing smooth numbers near kN). The 1.2114 implication is heuristically correct but needs explicit constructor requirements C1–C8 ([`sicheck`](notes/sicheck_results.md)). Every oblivious constructor costs L[2/3] per output: decoders gain nothing over the trivial divisor trick (Theorem J, [`sihard`](notes/sihard_results.md)), and meet-in-the-middle doesn't amortise. Prop SI is promise-CVP at factor 1+O(1/log N) in Micciancio's prime lattice ([`siconstruct`](notes/siconstruct_results.md) Prop CV).
- Representations ([`factorbits`](notes/factorbits_results.md)): addition in the prime-exponent representation is equivalent to factoring (a successor oracle factors N in popcount(N)−1 calls), mirroring Zech logarithms ⇔ discrete log. No natural single bit of a factor is known to be factoring-hard, because no computable map moves the secret factor. Theorem PT: methods that use N only as a modulus learn nothing about the factors except by splitting N.
- Phase 22:
  - [`zroute`](notes/zroute_results.md): deterministic subexponential factoring under ERH reduces to one constructive lemma, DSR (Theorem DZ). Everything else is proved: XK-Z, Prop PIN, rank forced deterministically.
  - [`selfred`](notes/selfred_results.md): no modulus-changing self-reduction exists in algebraic classes (SR1, SR2).
  - [`sinonobliv`](notes/sinonobliv_results.md): the Prop SI routes give nothing; one exception stays open, adaptive CRT decoding far below the Johnson radius.
- Phase 23, attacks on the three precise targets left open by phase 22. All three are negative:
  - [`dsr`](notes/dsr_results.md): DSR is not proved. Derandomising the JMV walk fails for every black-box tool (Theorem SB). Bach-type certificates fail because smooth sets have a Fourier-algebra norm of at least N^{1/2−o(1)} (Theorem FN). Short multiplier lists fail for infinitely many N (Prop BOX-∞).
  - [`ecount`](notes/ecount_results.md): #E(Z/N) for one fixed curve is NOT open. A single count factors N (Okamoto–Uchiyama 1998), and the agent made this rigorous outside a density-zero set. Dieulefait–Urroz needs either the factorisation of the count or counts on several curves.
  - [`crtdecode`](notes/crtdecode_results.md): no adaptive CRT decoder beats the divisor trick. Counting decoder calls cannot give a barrier, because the leaf computation is the sieve (Prop DT′). Barriers hold only for restricted classes.
- Phase 24 asked for a leak that is not splitting. Every known leak is one of three families: group order/smoothness, archimedean, or quadratic genus. All four agents found no fourth family and no cheaper smooth construction:
  - [`leakmap`](notes/leakmap_results.md): the three-family classification is proved for N-dependent exponent tests (NE1, matching Monier's liar rate), for AC⁰ on the bits of N (ε is unpredictable, unconditionally), and for pure trees plus Jacobi symbols. Integer representatives are too powerful a model for any such theorem. The strongest open model is "PRIM": retained, nested x^{E(N)} gates.
  - [`genusleak`](notes/genusleak_results.md): every route to an ambiguous form goes through a group order or a relation. Walk returns are smooth values (Lemma GIRTH), and a multiple of hR finds the form in poly time but costs |D|^{1/5}. Generic barrier: Theorem GB.
  - [`newfamily`](notes/newfamily_results.md): reciprocity invariants either depend only on N mod m (Theorem RA) or are factoring-equivalent. Dedekind sums (Prop DSF) and Salié sums (Prop SAL) are both factoring-equivalent. A pre-registered battery of 84 statistics found no signal.
  - [`smoothbuild`](notes/smoothbuild_results.md): constructing smooth numbers costs the same as sampling, index calculus or factoring. "Known factorisation" is no cheaper currency than smoothness (Lemmas KI/LP/AT). Open interface: hidden-root short-quotient navigation, which at L[1/3, 1.20] per call would give rigorous L[1/3, 1.80].
- Phase 25 worked on the open ends:
  - [`navigate`](notes/navigate_results.md): hidden-root short-quotient navigation is not a stepping stone. Without its height bound it is FACT ∧ DL (Prop H), and a solver for all N would give prime-field DL in L[1/3, 1.80] (Theorem DLB).
  - [`prim`](notes/prim_results.md): no leak in PRIM. The classification is proved for fresh-base chains of any depth, sparse gate sums (Theorem PS), two-layer gates off exceptional primes (Theorem D2) and bounded-height Jacobi tests (Theorem J1). It stays open at deterministic depth ≥ 2.
  - [`classE`](notes/classE_results.md): no leak. A CF convergent of a/N sees p only if a = uq ± 1 (Prop VIS), and for ⌊√N⌋/N the path is equidistributed up to denominators N^{7/12} (Theorem F-SQ).
  - [`tracedata`](notes/tracedata_results.md): #E(Z/N) for one fixed curve factors every good N (Theorem EXT, a degree-5 twist). One curve's traces (a_p, a_q) or a_N alone are still open; a_N for two CM curves factors N.
  - [`multidescent`](notes/multidescent_results.md): a two-level descent gives **L[3/8], conditional** on Hyp PC₃ and U′ (Theorem 2L). A union bound reduces it to three-point conditional statements (Prop UB), so a proof needs no control at d^{3/2} points. j levels would give (j+1)/(3j+2).
- Phase 26:
  - [`uprime`](notes/uprime_results.md): Hyp U′ is proved at 3/8 via a hybrid restricted BV (Lemma 8.5H; moduli B-smooth, values B₂-smooth), two-parameter level-1 moments and a mixed box scan.
  - [`pc3`](notes/pc3_results.md): PC₃ is stronger than needed. Decoupling (Prop DC) replaces it by two pair-only variance hypotheses C₁, C₂ (Theorem 2L-DC); both are still open.
  - [`geometry`](notes/geometry_results.md): no geometric shortcut. The area and genus of X_0(N) encode p+q, but every explicit model has size Θ(N); #X_0(N)(F_{ℓ²}) and #X_1(N)(F_ℓ) factor N, while #X_0(N)(F_ℓ) carries only QR bits.
- Fast factoring has no free consequences elsewhere in mathematics.

All of these are in [`APPROACHES.md`](APPROACHES.md).

## Open problems that would push further

- **Below L[2/5]:** prove the pair-only variance hypotheses C₁, C₂ ([`pc3`](notes/pc3_results.md) Theorem 2L-DC; U′ is done in [`uprime`](notes/uprime_results.md)) ⇒ rigorous L[3/8]. The main missing piece is a weighted BDH (Harper) + large sieve under the incidence measure.
- **Deterministic worst case below N^{1/5}:** coverage with λ < J^{1/4} ([`coverage38`](notes/coverage38_results.md)).
- **Human refereeing** of [`paper2/`](paper2/main.tex) §§5–6 and its Appendix A.

## Repository map

| Path | Contents |
|---|---|
| `paper2/` | Private draft preprint: Theorem A, RA/P, DL (54 pp, `tectonic paper2/main.tex`) |
| `paper/` | Draft paper: Theorem AV, D′, R, T1 (exponential regime) |
| `notes/` | One `*_results.md` per work item (verdict first), pre-registrations `*_prereg.md`, audits `*_audit.md` |
| `STATUS.md` | Chronological log of all phases and agent ownership |
| `APPROACHES.md` | Catalogue of everything tried, with verdicts |
| `src/`, `tests/` | Python code per work item; pytest suite |
| `codex/` | GPT-6.1-sol prompts and outputs; `codex/run.sh` wrapper |
| `lean/`, `leanA/` | Lean 4: D′ lemmas; conditional Theorems A/RA/P (cost assembly + glue proved; FM, LD* and 4 published theorems as hypotheses; checked by Comparator --paranoid, Nanoda, lean4lean) — see [`leanA_results`](notes/leanA_results.md) |
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
