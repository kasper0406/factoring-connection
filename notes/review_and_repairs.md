# Review and repair history (index)

Chronological index of every review, audit and repair in this project, with pointers. Created by agent `writer`,
2026-10-17. This file is append-only for later agents. The authoritative adjudication of the r3 audit is
`notes/repair_adjudication.md`.

| # | Date | Review / repair | Reviewer (independence) | Scope | Outcome | Where recorded |
|---|---|---|---|---|---|---|
| 1 | 2026-10-17 | GPT-6.1-sol round 1, internal critic and complexity-auditor passes | GPT-6.1-sol, same run (not independent) | Its own six mechanism cards (A1, B1, C1, D1, E1, F1) | B1 degree-2 rank guarantee refuted; C1 and D1 reduced to known predicates | `notes/codex_r1_mechanisms.md` §2; `codex/PHASE1_NOTES.md` |
| 2 | 2026-10-17 | Orchestrator proof check | main session | Theorem A1 (infraA); O1, O2, O4 (factC) | Proofs held as stated (A1's index hypothesis and O1's polynomial-readout corollary were later corrected, rows 6–7) | `STATUS.md` log |
| 3 | 2026-10-17 | factC check of GPT-6.1-sol r1 Tal/BGS claims against the PDFs | agent `factC` (independent of GPT-6.1-sol) | Tal ECCC TR26-211, BGS 2007 | All spot checks confirmed; HH adaptation (C-HHloc) not tested | `notes/trackC_tal_audit.md` §4 |
| 4 | 2026-10-17 | Independent read-only Opus critic on Track B/D | Opus subagent (independent) | `trackBD_analysis.md` §§2–5 and code | 15 findings, all incorporated (Euclid iff, squarefree hypothesis in B3, Model-M overclaim, non-equidistribution of `q mod (p−1)`, BM step, Stirling window, D1 hypotheses) | `notes/trackBD_analysis.md` header |
| 5 | 2026-10-17 | Independent read-only Opus critic on Track E | Opus subagent (independent) | `trackE_analysis.md` §§1–6 and code | Core results correct; 7 fixes (ν vs bit length, F2 restriction for `m = 2ℓ`, per-prime swap proof, Pollard case, heavy tails, model-conditional labels, oracle routing) | `notes/trackE_analysis.md` §7a |
| 6 | 2026-10-17 | GPT-6.1-sol round 3 adversarial critic audit | GPT-6.1-sol (independent of the track agents) | Tracks A, B/D, C proofs and claims; unified gap theorem; novelty calibration | 36 WRONG/GAP items raised | `notes/codex_r3_critic.md` |
| 7 | 2026-10-17 | Adjudication of the r3 audit, with exact reproductions | agent `repair` | All 36 items | **Critic right on 36/36** (27 claims changed, 9 scoped or clarified; wrong on 0). Notes corrected in place with dated lines; Corollary 3′ added; new D3 BM bound proved | `notes/repair_adjudication.md`; `src/repair/audit_track{A,C,BD}.py`; `tests/test_repair_audit.py` |
| 8 | 2026-10-17 | STATUS corrections | orchestrator | STATUS.md log lines repeating corrected claims | Corrections section added | `STATUS.md` (corrections list), commit `fcbcb6b` |
| 9 | 2026-10-17 | Pre-registered C-2R repair cycle | agent `repair` | Two-layer interval folding (design from r3 §4) | Pre-registered at `a8920b7`; audit numbers 0/32, 14/32, 0/32 reproduced; H-C2R rejected (0/16); frozen extension X1 rejected (0/16); Conjecture C2R-X stated | `notes/repair_C2R_prereg.md`, `notes/repair_C2R_results.md` |
| 10 | 2026-10-17 | Exhaustive verification of GPT-6.1-sol round-4 claims | agent `writer` (independent of GPT-6.1-sol) | `s = ε(N+1) − 16εuv`, mod-16/32/64 statements, 4-rank criterion, public class `T`, genus matrices, Milović Prop. 1 as quoted | All hold for every pair `p < q < 1000` in scope (6 901 pairs; 3 418 with `N ≡ 1 mod 8`); wrong-discriminant control agrees on only 46%. The genus-matrix proof itself was not re-derived | `src/repair/verify_r4.py`, `tests/test_repair_r4.py`, `experiments/raw/verify_r4.json`, `REPORT.md` §1.1 |

## Not independently reviewed (as of row 10)

* GPT-6.1-sol r2 theorems: exact exposure formulas (2)–(5), nonadaptive ring+gcd equivalence, constant-only SLP bound,
  restricted permanent barrier (`notes/codex_r2_generic_ring.md`).
* GPT-6.1-sol r1 F1 torus success formula and E1 mask density claims (`notes/codex_r1_mechanisms.md`).
* GPT-6.1-sol r4 genus-matrix proof of the 4-rank criterion, and its literature statements (Jager–Schwenk scope, Schoof,
  Rybalov, Milović) (`notes/codex_r4_sos_classgroup.md`).
* GPT-6.1-sol prior-art map entries not re-read by an agent (`references.md` marks provenance).
* C-HHloc (Harvey–Hittmeir restricted to a certified interval): not yet tested.

## Appended 2026-10-17 (agent `writer`): r5 report review

| # | Date | Review / repair | Reviewer (independence) | Scope | Outcome | Where recorded |
|---|---|---|---|---|---|---|
| 11 | 2026-10-17 | GPT-6.1-sol round 5 review of the first REPORT.md | GPT-6.1-sol (independent of `writer`); orchestrator agreed with all 14 items | REPORT.md, NEXT_STEPS.md | 14 items. Two computational items verified by `writer`: item 3 (`uv ≡ b(1+2[N≡1 mod 16]) mod 4` for `N ≡ 1 mod 8`, so `s mod 64` needs no class-group rank; 3 418 pairs `p<q<1000`, 0 failures) and item 10 (`p=257, m=8, E=70, v=186` is an exact non-character fold with `d = 7 < 9`; two implementations). New diagnostic: `s mod 128` is not determined by `r4`, `r8` of `Cl(−4N)`, `Cl(−8N)` (73/248 keys ambiguous). NEXT_STEPS rewritten for items 1, 3, 7, 9, 11; revised report text handed to the orchestrator | `codex/r5_report_review.out`; `src/repair/verify_r5.py`; `tests/test_repair_r5.py`; `experiments/raw/verify_r5.json`; `NEXT_STEPS.md` |
