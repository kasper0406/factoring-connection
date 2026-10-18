# Common instructions for Phase 5 agents

- Read notes/AGENT_RULES.md, STATUS.md (Phase 5 table; stay inside your files), and notes/dcltheorem_results.md (Theorem D′, which is authoritative). Also read notes/multimod_results.md, notes/codex_r12_dcl_audit.md and notes/codex_r13_final_referee.md. Skim notes/cycloharvey_results.md for the calibrated cost model.
- Reuse src/dcltheorem/, src/multimod/ and src/cycloharvey/ by import. Never edit other agents' files.
- Commit LOCALLY only, adding your own files by explicit path; never push. End every commit message with "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>". The only edit you may make to STATUS.md is to append one line, prefixed "- 2026-10-18: **<agent>**".
- **Pre-register before experiments:** notes/<agent>_prereg.md, committed first.
- **GPT-6.1-sol (Codex CLI)** is available for blind attempts, being stuck, or adversarial review. Write the prompt to codex/<agent>_<k>.md, then run it in the background (~10–40 min):
  `cd /home/knielsen/code/factoring-connection/codex && timeout 3000 codex exec -s read-only -m gpt-6.1-sol -c model_reasoning_effort="high" -c tools.web_search=true --skip-git-repo-check -o <agent>_<k>.out - < <agent>_<k>.md > <agent>_<k>.log 2>&1`
  Pass the prompt on stdin with `-` and do NOT add a second `< /dev/null`.
- **Claims:** label everything proved / experimentally supported / conjectured / refuted for this model / UNVERIFIED. Any claimed new theorem must pass an independent read-only Opus critic subagent before you report it.
- **Final message:** a summary of at most 350 words that opens with a plain-language verdict.

## Addendum (2026-10-18, phase 5 restart on a new machine)
- The project root is now **/Volumes/git/factoring-connection** (macOS). Read every `/home/knielsen/code/factoring-connection` path in the notes as this root. `.venv` was recreated: Python 3.12, sympy 1.14, gmpy2, python-flint 0.9.0, numpy, pytest. fpylll is not installed. The full suite passes (369 tests).
- macOS has no `timeout`. Run codex with `codex/run.sh <agent>_<k> [effort]`: the prompt goes in `codex/<agent>_<k>.md` and the outputs land in `codex/<agent>_<k>.out` and `.log`. The script uses gpt-6.1-sol with web search and a 80-minute alarm. The repo is read-only for codex, which writes scratch files to `codex/scratch/<name>/` (git-ignored). Launch it in the background, e.g. `(nohup codex/run.sh dchworst_2 xhigh >/dev/null 2>&1 &)`, then poll for the .out file.
- The first machine was lost mid-phase. The orchestrator already re-launched the `<agent>_1` blind prompts for sieve2d, unbalanced and dchworst; read `codex/<agent>_1.out` once it exists. Do not relaunch them.
- Pre-registrations committed before the loss (`notes/{sieve2d,unbalanced,dchworst}_prereg.md`) remain frozen and valid. Run exactly what they specify. If you must deviate, record it as a disclosed deviation.
