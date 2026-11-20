(Read /Volumes/git/factoring-connection/codex/hatch_common.md and follow it.) HATCH 3, R / reuse.
Context: notes/lowerbound_results.md (R: extraction needs h distinct fresh successes in one block, log h ≥ (1−o(1))μ;
reuse budget; tree incidence T), notes/zroute_results.md (Prop PIN: only small-prime coordinates need pinning; XK-Z;
seedless PIN), notes/heurnfs_results.md (Assumption T, descent), notes/detGRH_results.md (Prop FL fruitfulness),
notes/fruit_results.md, notes/lowerbound_results.md batch discussion (dense-batch adversary A1), codex/lowerbound_2_*.
Questions:
 (a) Thin dependencies: is R a genuine linear-algebra necessity (rows > active columns) or can an algorithm extract a
     factor from far fewer than π(B) relations — e.g. relations concentrated on few columns, a Z-route PIN that needs
     only the ~2 log²N small primes pinned, large-prime cycles, quadratic-character/Schirokauer-map tricks, or a
     "planted" relation family whose dependency is known by construction? What is the true minimal row count, and
     what would a thin-dependency algorithm cost under H?
 (b) Shared descent DAGs: can a heavily shared DAG (one eliminator serving many parents) beat per-leaf pricing?
 (c) Reuse across related inputs: one witness serving many inputs (multipliers kN, N^j, factor-base sharing across
     k, many N in a structured family) — can per-target cost drop below L[1/3] for a SINGLE given N by manufacturing
     related auxiliary instances? Ledger each; give the best candidate, missing lemma and a toy test (python3).
Use only read-only shell commands plus python3 in your cwd. Final answer ≤ 2000 words, plain text.
