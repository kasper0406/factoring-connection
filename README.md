# factoring-connection

Exact computational research on whether a compact algebraic representation of `F_N = (1+X)^N − 1 − X^N` (or a related
object) yields a cheap operation that reveals a factor of `N`. Brief: `BRIEF.md`. Final findings: `REPORT.md`; next steps:
`NEXT_STEPS.md`; sources: `references.md`. No polynomial-time factorizer was found.

## Layout

| Path | Contents |
|---|---|
| `BRIEF.md` | The research brief (mission, tracks A–F, audit rules, deliverables) |
| `STATUS.md` | Chronological log, ownership table, and the corrections section after the r3 audit |
| `REPORT.md`, `NEXT_STEPS.md`, `references.md` | Final deliverables (BRIEF §7) |
| `notes/AGENT_RULES.md` | Rules every agent followed (exact arithmetic, no leakage, status labels) |
| `notes/seed_lemmas.md` | Proofs and checks of the brief's elementary lemmas (L1–L7) |
| `notes/trackA_*.md` | Binary coefficient sections: Theorem A1 (dimension conservation), C1–C4, selectors, masks |
| `notes/trackC_*.md` | Factorial unit predicate: BGS baselines, mechanisms M1–M8, obstructions O1–O4, Tal audit |
| `notes/trackBD_*.md` | Quotient sketches, lifts, annihilators: B1–B7, D1–D3 |
| `notes/trackE_*.md` | Coppersmith decoders, partial-information formats, information floor |
| `notes/repair_adjudication.md` | Item-by-item adjudication of the r3 critic audit (authoritative for corrected claims) |
| `notes/repair_C2R_prereg.md`, `notes/repair_C2R_results.md` | Pre-registered C-2R repair cycle and its results |
| `notes/codex_*.md`, `notes/prior_art_map.md` | Curated GPT-6.1-sol output (rounds 1–4, prior-art map); raw prompts and outputs in `codex/` |
| `notes/review_and_repairs.md` | Index of the review and repair history |
| `src/seed_checks.py` | BRIEF §9 smoke test plus extended checks |
| `src/common/` | Shared ring polynomials, seeded instance generator, metrics, baselines (`src/common/README.md`) |
| `src/circuit_sections.py`, `src/trackA/` | Track A code |
| `src/quotient_sketches.py`, `src/trackBD/` | Track B/D code |
| `src/trackC/` | Track C code |
| `src/trackE/` | Track E code (`coppersmith.py`, `formats.py` operational; `oracles.py` harness) |
| `src/repair/` | r3 audit reproductions, C-2R search, and `verify_r4.py` (exhaustive check of the r4 class-group claims) |
| `tests/` | pytest suite |
| `experiments/configs/`, `experiments/raw/`, `experiments/analysis/` | Configs, raw JSONL results, aggregated tables |

Operational code (anything a candidate algorithm runs) receives only `N`, public configuration and a seed. Code that knows
the factors is named diagnostic, harness or oracle, and is never imported by operational code.

## Running

```bash
cd /home/knielsen/code/factoring-connection
.venv/bin/python -m pytest -q                      # full test suite (~5 s)
.venv/bin/python src/seed_checks.py                 # BRIEF §9 smoke test
.venv/bin/python -m src.repair.verify_r4 --bound 300 --out experiments/raw/verify_r4.json
```

The environment is the project virtualenv `.venv` (Python 3.13; sympy, gmpy2, python-flint, numpy, fpylll, pytest).
Experiment commands for every track are collected in `REPORT.md` §3.3 and in the reproduction section of each track note.
The instance generator `balanced_semiprime` uses primes with the top two bits set, so `q/p < 4/3` for even bit lengths.
