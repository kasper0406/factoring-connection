# STATUS

**Project:** classical factoring connection (brief: `BRIEF.md`, received 2026-10-17)
**Orchestrator:** main Claude session. Coding is done by Opus 5.5 subagents. Critique, literature, and blind mechanism generation go to GPT-6.1-sol via `codex exec` (prompts and outputs live in `codex/`).

## Ownership (to avoid collisions)
| Area | Owner | Files |
|---|---|---|
| Shared infra: seed checks, ring polys, instance generator, metrics | agent `infraA` | `src/seed_checks.py`, `src/common/*`, `tests/test_seed*.py`, `tests/test_common*.py` |
| Track A: binary sections | agent `infraA` | `src/circuit_sections.py`, `src/trackA/*`, `notes/trackA_*.md`, `experiments/*/trackA*` |
| Track B + D: quotient sketches, lifts, pivots | agent `sketchBD` | `src/quotient_sketches.py`, `src/trackBD/*`, `notes/trackBD_*.md`, `experiments/*/trackBD*` |
| Track C: factorial unit predicate | agent `factC` | `src/trackC/*`, `notes/trackC_*.md`, `experiments/*/trackC*` |
| Prior art, mechanism generation, critique | codex GPT-6.1-sol | `codex/*` → curated into `notes/prior_art_map.md`, `notes/candidates_*.md` |

## Log
- 2026-10-17: workspace created; venv `.venv` (sympy 1.14, gmpy2 2.3.2, python-flint 0.9.0). Codex round 1 launched (mechanisms / prior art / Tal+BGS). Opus agents launched on Phase 1 and Tracks A, B/D, C.
- 2026-10-17: [infraA] shared infra READY: `src/common/{ring,instances,metrics,baselines}.py` + `src/common/README.md` (API, flint audit); `src/seed_checks.py` (BRIEF §9 + checks to N≤1000 incl. Lucas/Kummer support); 32 pytest tests pass (`.venv/bin/python -m pytest -q tests`).
- 2026-10-17: **factC done.** Built BGS / Pollard–Strassen baselines (slopes ~0.56 / 0.62 in K). The Tal audit was checked against the PDF: no sub-sqrt step survives over unfactored R_N. Proved restricted-model obstructions for the prefix/interval unit predicate: O1 (matrix powers need order ≥ L+1), O2/O2′ (power-map layers), O3 (degree), O4 (subgroup prefix ⇒ (p−1)|K). Orchestrator checked the O1/O2/O4 proofs and they hold. Escape requirements are in notes/trackC_analysis.md §3. Open question: are there polylog SLPs with interval zero sets mod p, uniform in p? (the mod-p form of Shub–Smale Thm 2). Test suite: 72 passed.
- 2026-10-17: **infraA done** (80 tests pass). Theorem A1 (dimension conservation; orchestrator checked the proof): any rank-D algebra observable that isolates residue class r mod M needs D ≥ M. This holds for nonlinear computation of f(x), for support-aware access, and for multilinear read-outs. It is the same linear-recurrence gap principle as Track C's O1, so **one gap lemma blocks both the section and prefix routes for all finite-algebra evaluation**. Experiments: section circuits grow ×3.7 per level; the best selector is exact but doubles its cost per level (exponential). There are zero scalar false zeros among 44k nonzero sections, so zero-testing is NOT the obstacle. Open escape E-1: binomial-specific SLPs for sections that do not go through algebra evaluation.
- 2026-10-17: **sketchBD done** (23 tests). B1 (proved): F_N ≡ F_m^{p^k} mod p; sketch roots depend on q only mod L_α = lcm(ord α, ord(1+α)); Frobenius coincidence is sufficient but not necessary; sporadic roots come in 6-orbits. Rank splits ⇔ pivot profiles differ. The Track D divided lifts are unit multiples of G_k = (1+X)^{N−k} − X^{N−k}, i.e. Fermat-type order tests. BM over c gives linear complexity q mod (p−1), a Fermat-gap method. Every single-sketch readout tested succeeds at Θ(1/p) per sketch, so model M is refuted. One independent Opus critic pass (15 findings, all fixed). NOTE: `balanced_semiprime` uses top-2-bit primes, so q/p < 4/3; this is an RSA-style restriction and is recorded as such.
- 2026-10-17: codex round 2 (generic-ring synthesis) curated → notes/codex_r2_generic_ring.md. trackE agent launched (Coppersmith decoders + info-floor table).
