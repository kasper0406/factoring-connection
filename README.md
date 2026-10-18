# factoring-connection

An exact computational research project on **classical deterministic integer factoring**. The
work was AI-assisted: Claude agents wrote the code, proofs and notes, and GPT-6.1-sol was used
for literature sweeps and adversarial review (see [How this was done](#how-this-was-done)).

The project started from one question: does a compact algebraic representation of
F_N(X) = (1+X)^N − 1 − X^N, or a related object, give a cheap operation that reveals a factor
of N? ([`BRIEF.md`](BRIEF.md)) That line closed negatively. The work then moved to the
deterministic collision-search methods of Lehman, Hittmeir, Harvey and Harvey–Hittmeir, and
produced one positive result.

**No polynomial-time or subexponential factoring algorithm was found.**

## Headline result (average case only)

**Theorem D′ (proved; reviewed by AI referees, not yet by a human).** Fix ε ∈ (0, 1/60]. For
all large n there is a deterministic algorithm with the following properties:
- it sees only an n-bit integer N;
- it uses N^{1/6+5ε/6}·(log N)^{O(1)} bit operations on every input;
- it outputs only a proper divisor of N or FAIL;
- it factors all but an O(N^{−ε}) fraction of the semiprimes N = pq with primes p < q ≤ 2p.

With polylogarithmic parameters the cost is N^{1/6}·polylog, with an O((log N)^{−A})
exceptional fraction (Corollary D″).

For comparison, the best proved deterministic exponent for *every* input is 1/5 (Harvey 2021;
Harvey–Hittmeir 2022).

**Honest scope:**
- **Average case only.** The worst case is open; it needs an unproved divisor-cover hypothesis
  (DCH) among other things.
- **Bounded ratio only.** Only q/p ≤ 2 is covered.
- **Not faster in practice.** At 40–80 bits the measured cost grows more slowly than that of
  a calibrated Harvey–Hittmeir implementation (slope 0.203 against 0.226), but
  Harvey–Hittmeir is still 2.3–4.9× cheaper at those sizes.
- **Built on prior art.** Every ingredient is known, and McKee–Pinch 1996 §5 is a direct
  precursor. We found no prior formulation of this combination.
- **A floor for full covering only.** In a one-base sumset collision-search model, any family
  that covers *every* candidate factor of a fixed N costs at least N^{1/6} (Theorem A, proved in
  the model). This is not a lower bound for factoring, and it does not apply to average-case
  algorithms such as Theorem D′. Whether one can go below 1/6 on almost all inputs is open.

Details:
- Paper draft: [`paper/main.tex`](paper/main.tex). It is a draft, not for distribution.
- Authoritative note: [`notes/dcltheorem_results.md`](notes/dcltheorem_results.md).
- Model lower bounds: [`notes/multimod_results.md`](notes/multimod_results.md).

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
| `paper/` | LaTeX draft of Theorem D′ and Corollary D″ (`main.tex`, `refs.bib`) |
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
| `lean/` | Lean 4 formalisation of some Theorem D′ lemmas (in progress) |
| `COMMIT_SHA_MAP.md` | Old → new commit SHAs after the history rewrite (see below) |

**Leakage rule.** Operational code, meaning anything a candidate algorithm runs, receives only
N, public configuration and a seed. Code that knows the factors is named diagnostic, harness or
oracle, and operational code never imports it.

## Running the tests

The environment is a local virtualenv: Python 3.13 with sympy, gmpy2, python-flint, numpy,
fpylll and pytest.

```bash
python3.13 -m venv .venv
.venv/bin/pip install sympy gmpy2 python-flint numpy fpylll pytest
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
