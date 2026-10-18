# factoring-connection

An exact computational research project on **classical deterministic integer factoring**. The
work was AI-assisted: Claude agents wrote the code, proofs and notes, and GPT-6.1-sol was used
for literature sweeps and adversarial review (see [How this was done](#how-this-was-done)).

The project started from one question: does a compact algebraic representation of
F_N(X) = (1+X)^N − 1 − X^N, or a related object, give a cheap operation that reveals a factor
of N? ([`BRIEF.md`](BRIEF.md)) That line closed negatively. The work then moved to the
deterministic collision-search methods of Lehman, Hittmeir, Harvey and Harvey–Hittmeir, and
produced the positive results below.

**No polynomial-time or subexponential factoring algorithm was found.**

## Headline results (balanced semiprimes; average case)

All results below concern S_n, the n-bit semiprimes N = pq with primes p < q ≤ 2p. Each is
**proved**, and each was reviewed by GPT-6.1-sol and an independent Claude Opus critic, with
"accept with fixes" and the fixes applied. **None has yet been checked by a human.**

| Result | Statement | Status |
|---|---|---|
| **Theorem AV** | A complete deterministic algorithm (it always outputs the full factorisation). Its worst case is Harvey–Hittmeir's N^{1/5}; its **average cost over uniform S_n is N^{1/6}·(log N)^{O(1)}** | proved ([`notes/avgtime_results.md`](notes/avgtime_results.md)) |
| **Theorem D′ + fourth moment** | One run of the search sees only N, is never wrong, costs N^{1/6+5ε/12}·polylog on every input, and fails on O(N^{−ε}) of S_n (ε ≤ 1/40). With polylog parameters: N^{1/6}·polylog with (log N)^{−A} exceptions | proved ([`notes/dcltheorem_results.md`](notes/dcltheorem_results.md), [`notes/sieve2d_results.md`](notes/sieve2d_results.md)) |
| **Theorem R** | The worst case reduces to coverage alone. A modified search factors every N ∈ S_n whose pair is caught by some band. So the pair-covering hypothesis DCH_pair ⇒ worst case N^{1/6+o(1)}, and any covering theorem with exponent c > 3/8 ⇒ worst case < 1/5 | proved; the hypotheses are **open** ([`notes/dchworst_results.md`](notes/dchworst_results.md)) |
| **T1/T2** | In a sumset collision-search cost model, residue-local degree-1 designs reach the true factor of a constant fraction of S_n only at cost ≥ N^{1/6}/polylog. The search attains this up to polylog | proved in the model ([`notes/below16_results.md`](notes/below16_results.md)) |

For comparison, the best proved deterministic exponent for *every* input is 1/5 (Harvey 2021;
Harvey–Hittmeir 2022). We are not aware of an earlier proved average-case or almost-all
deterministic exponent below 1/5 for balanced semiprimes. This is the outcome of a literature
search, not a proof of priority.

**Honest scope:**
- **Average case, not worst case.** The worst case is still N^{1/5}. Getting below it needs a
  covering theorem that is not proved.
- **Ratios.** The average-case theorem covers q/p ≤ 2. The almost-all statement extends to
  q/p ≤ N^{o(1)} (Theorem U, proved). For q/p = N^θ the best exponent we have is
  min((1+θ)/6, (1−θ)/4) (Theorem U_θ plus Pollard–Strassen). It peaks at 1/5 when q/p ≈ N^{1/5},
  and in a full-coverage model that peak is unavoidable (Proposition O′, a model theorem; see
  [`notes/unbalanced_results.md`](notes/unbalanced_results.md)).
- **Not faster in practice.** At 40–80 bits Harvey–Hittmeir is 2.3–4.9× cheaper in our cost
  proxy.
- **Built on prior art.** The key ingredients come from McKee–Pinch 1996 §5 (the divisibility
  event), Harvey–Hittmeir 2022 Prop. 3.3 (the lock), Harvey 2021 Alg. 4.1–4.2 (shared collision
  search, exact matches), Harvey–Hittmeir 2026 (order tools) and the large sieve.
- **Lower bounds are model theorems, not factoring lower bounds.** The full-candidate N^{1/6}
  covering floor (our Theorem A) is **pre-empted** by Ward–Maestro (iGent AI report,
  2 Sept 2026, Thm 3.6), whose model has no class locks. Pollard p − 1 already factors a positive
  proportion below N^{1/6} (Ward–Maestro Prop. 35.3); this lies outside our model.
- **Carella's claimed N^{1/6}** (arXiv:1308.2891) is not established: the §5 polynomial has
  content divisible by B, and after removing it the small-root condition fails.

Details:
- Paper draft: [`paper/main.tex`](paper/main.tex) (compiled with tectonic). It is a DRAFT, not
  for distribution. Review adjudication: [`notes/paper_review_r17.md`](notes/paper_review_r17.md).
- Also done: NFS fruitfulness, with no rigorous L[1/3] ([`notes/nfsfruit_results.md`](notes/nfsfruit_results.md)),
  and a Lean 4 formalisation of D′'s lemmas ([`notes/leanform_results.md`](notes/leanform_results.md)).

## Everything else that was tried

[`APPROACHES.md`](APPROACHES.md) catalogues every approach in the project. It covers:
- the Phase 1 tracks A–E;
- 13 GPT-6.1-sol rounds;
- the repair cycles;
- the ε battery, the s mod 128 tower and primality-test leakage;
- cyclotomic handles, sign-splitting and several moduli;
- divisor covers and Theorem D′;
- the user ideas.

Each entry has a verdict, its key numbers, and pointers to the notes, code and commits. The
entries are grouped into walls/obstructions, refuted ideas, positive results and open problems.

[`REPORT.md`](REPORT.md) is the Phase 1 final report, covering the F_N tracks, with
[`NEXT_STEPS.md`](NEXT_STEPS.md). It is kept as written, apart from corrections it records
itself. Later phases are logged in [`STATUS.md`](STATUS.md).

## Repository map

| Path | Contents |
|---|---|
| `paper/` | LaTeX draft: Theorem AV, D′ with the fourth-moment bound, Theorem R, T1 (`main.tex`, `refs.bib`) |
| `APPROACHES.md` | Catalogue of everything tried, with verdicts and pointers |
| `BRIEF.md` | The original research brief (mission, tracks A–F, audit rules) |
| `REPORT.md`, `NEXT_STEPS.md` | Phase 1 final report and its follow-up list |
| `STATUS.md` | Chronological log of all phases, agent ownership tables, and corrections |
| `references.md` | Primary sources consulted in Phase 1, with provenance codes |
| `notes/` | One results note per work item (`*_results.md`, `*_analysis.md`) and its pre-registration (`*_prereg.md`) |
| `notes/codex_*.md` | Curated GPT-6.1-sol rounds r1–r13 |
| `notes/prior_art_map.md` | Primary-source prior-art map |
| `notes/user_ideas/` | Ideas proposed by the human researcher, with assessments |
| `notes/AGENT_RULES.md` | Rules every agent followed |
| `src/` | Python code, one package per work item (`trackA`, `trackBD`, `trackC`, `trackE`, `repair`, `eps`, `tower`, `primality`, `pm1bsgs`, `cycloharvey`, `signsplit`, `multimod`, `divcover`, `dcltheorem`, …) plus shared infrastructure in `src/common/` |
| `src/dcltheorem/` | The Theorem D′ algorithm (`algorithm.py`), parameters, harness, Harvey–Hittmeir baseline, experiments |
| `tests/` | pytest suite |
| `codex/` | Raw GPT-6.1-sol prompts (`*.md`), outputs (`*.out`) and logs |
| `experiments/configs/`, `experiments/raw/`, `experiments/analysis/` | Configurations, raw JSONL results, aggregated tables and summaries |
| `references/` | Scans of McKee–Pinch 1996 (§5 is the direct precursor of Theorem D′) |
| `lean/` | Lean 4 formalisation of Theorem D′ lemmas (large sieve as an explicit hypothesis) |
| `COMMIT_SHA_MAP.md` | Old → new commit SHAs after the history rewrite (see below) |

**Leakage rule.** Operational code, meaning anything a candidate algorithm runs, receives only
N, public configuration and a seed. Code that knows the factors is named diagnostic, harness or
oracle, and operational code never imports it.

## Running the tests

The environment is a local virtualenv: Python 3.12 or later with sympy, gmpy2, python-flint,
numpy and pytest. fpylll is optional; the current machine runs without it, and the suite passes.

```bash
python3.12 -m venv .venv
.venv/bin/pip install sympy gmpy2 python-flint numpy pytest   # fpylll optional
.venv/bin/python -m pytest -q                         # full suite (about a minute)
.venv/bin/python -m pytest -q tests/test_dcltheorem.py
.venv/bin/python src/seed_checks.py                   # BRIEF §9 smoke test
```

Every results note ends with the exact commands that reproduce its experiments. For example,
Theorem D′'s experiments run as:

```bash
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments lemma   # also: config, hh, oper, exc, summarize
```

## History rewrite (commit SHAs)

On 2026-10-18 the git history was rewritten with `filter-branch` to remove one 110 MB raw file,
`experiments/raw/trackBD_tiny.jsonl`, which exceeds GitHub's file-size limit. Its gzipped copy
is kept.

Notes and `STATUS.md` cite pre-registration commits by their **original** SHAs. Look them up in
[`COMMIT_SHA_MAP.md`](COMMIT_SHA_MAP.md). For example, the Theorem D′ pre-registration was
`3550e24` and is now `05d2bf5`. `APPROACHES.md` uses current SHAs and gives the original SHA
next to each pre-registration commit.

## How this was done

- **The human researcher** set the brief and the research direction, proposed several ideas
  (`notes/user_ideas/`), and decided which threads to pursue.
- **An orchestrating Claude session** planned the work and reviewed the results.
- **Claude Opus 5.5 subagents** wrote all code, proofs, experiments and notes. They worked
  under written rules: exact arithmetic, no access to the factors in operational code, status
  labels on every claim, and pre-registration before experiments.
- **GPT-6.1-sol, via the Codex CLI,** did literature sweeps, blind attempts and adversarial
  referee reports.
- **Independent read-only Opus critics** reviewed each major result.

The review history matters for weighing the results:
- The r3 critic was right on all 36 of its findings.
- The r12 audit found the gap that turned Theorem D into Theorem D′.
- The r13 referee found a real output bug and the closest precursor.

Every claim carries a status label: *proved*, *experimentally supported*, *conjectured*,
*refuted for this model* or *UNVERIFIED*. No result here has yet been refereed by a human.
