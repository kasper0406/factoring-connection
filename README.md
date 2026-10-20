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
- Nested descent gives no route below 2/5 (Prop R2).
- Fast factoring has no free consequences elsewhere in mathematics.

All of these are in [`APPROACHES.md`](APPROACHES.md).

## Open problems that would push further

- **Below L[2/5]:** multi-level descent with "algebraic children" ([`nested`](notes/nested_results.md), [`recurse`](notes/recurse_results.md) §3).
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
| `lean/`, `leanA/` | Lean 4 formalizations: D′ lemmas; conditional Theorem A layer (in progress) |
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
