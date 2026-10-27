# Phase 28 survey brief (shared by conn_alg, conn_ana, conn_top, conn_cmp)

Goal: for each assigned field, look for **unexpected connections to integer factorization** that could lead to an
asymptotically faster algorithm (beating heuristic L[1/3, 1.902], or the rigorous L[2/5] / deterministic N^{1/5}
records). The user wants breadth and creativity, but verdicts must be honest.

Before starting, read: README.md, BARRIERS.md, notes/factorbits_results.md §3 (Theorem PT, fiber floor),
notes/leakmap_results.md §0 (the leak-family map: M group order/smoothness, A archimedean, Q quadratic genus),
notes/suffmap_results.md §0 (size–entropy law), notes/newfamily_results.md §0 and notes/geometry_results.md §0.
Use `grep -rl <keyword> notes/` to find prior coverage of each idea before analysing it.

For EACH field (aim for 3–6 concrete candidate connections per field):
1. State the connection concretely: what object is built from N alone, and what computable feature would carry
   factor information, or what algorithmic primitive would replace sieving or linear algebra.
2. Check prior coverage: in the repo (cite the note) and in the literature (WebSearch, or the web-fetch agent for
   key papers; cite them).
3. Verdict, one of:
   - RESTATEMENT (needs a factor to write down);
   - AGGREGATE (depends only on N mod m);
   - FACTORING-EQUIVALENT (or DL-equivalent);
   - REDUCES-TO-<M/A/Q or a named project barrier>;
   - KNOWN (already used in factoring, with reference);
   - OPEN-WORTH-A-LOOK.

   Give a 2–5 sentence justification. Be specific about WHERE it fails: the exact step that needs the factors,
   costs size N, or lands in a known family.
4. Mark anything surprising or new to the project.

Then rank your top 3 OPEN-WORTH-A-LOOK items, if any, with a precise first question or lemma for each. If none
survive, say so.

GPT-6.1-sol: one blind round EARLY asking for its own unexpected connections for your fields (codex/<agent>_1.md →
`(nohup codex/run.sh <agent>_1 xhigh >/dev/null 2>&1 &)`, then poll codex/<agent>_1.out). Merge its candidates into
your table, credited. One Opus read-only critic must review your OPEN-WORTH-A-LOOK items and any KNOWN attributions
before you finalise. Have it return its findings to you directly.

Rules:
- Survey work, so no pre-registered experiments are required. If you run any computation, keep it tiny (≤ 2 GB
  RSS, check `memory_pressure`) and leave no background processes.
- At most one sub-agent at a time. Use your own scratch directory, scratchpad/<agent>/.
- Deliverable: notes/<agent>.md, verdict first; then a table per field (connection | coverage | verdict |
  why/where it fails | refs); then the ranked shortlist.
- Append one STATUS.md line. Commit locally by explicit paths only; never push; never edit others' files. Commit
  messages end with "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>".
- Final message ≤ 300 words: verdict first, then the shortlist.
