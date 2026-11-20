(Read /Volumes/git/factoring-connection/codex/equivtarget_common.md first and follow it.)

## Your focus (C): systematic harvest + re-scoring
1. Grep notes/ (e.g. rg -n -i "factoring-equivalent|FACT-equivalent|≡ FACT|FACT ∧ DL|DL-equivalent|in disguise|would
   (also )?solve|prime-field DL|F_p DL|ZPP\^FACT|Okamoto|class-group DL" notes) and build a table of every route that
   was dropped for being equivalent to FACT or DL (not the 55 phase-29 invariants). For each: target T, file:line,
   reduction direction and cost exponent, and T's native landscape (best known algorithm for T and for special
   families of T).
2. Re-score: which T have a known easy special regime or a transfer (see preamble list), and could the T-instances
   arising from N be steered there? Pay special attention to reductions that land in DL of a group OTHER than F_p^*
   (class groups Cl(−4N), (Z/N)^× subgroups, E(Z/N), Jacobians, F_{p^k}, unit groups, Zech/representation problems).
3. Work out the top 2–3 in detail with ledgers; run a toy test if one is cheap.
