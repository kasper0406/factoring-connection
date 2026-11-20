# hatch: attack the three escape hatches of the relation-model lower bound (common preamble)

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch only in your cwd).
User's only goal: factoring general N in L_N[α], α < 1/3 (heuristic fine; constants worthless).

Read first: notes/lowerbound_results.md (the model theorem: priced valuation-probe model; assumptions H conditional
Dickman, R quota, credit; profile λ(k) ≥ A√(x log k) ⇒ L[1/3,(4A/3)^{2/3}]; the three escape hatches), then the notes
named in your hatch below. Any sub-1/3 algorithm must break one of:
 1. H — a constructor/decoder finding smooth residuals faster than Dickman-rate testing (e.g. a beyond-Johnson CRT
    list decoder on N-public residue tables; Prop SI constructors).
 2. the profile — exceptionally small norms beyond the error budget at rank ≥ 3 (LP_exc^rel; type-D multi-cone unit
    cancellation; Theorem KH only bounds generator heights).
 3. R / reuse — thin dependencies (few pinned coordinates, PIN), heavily shared descent DAGs, one witness serving many
    related inputs.
Goal of this round: make a DENT — either a mechanism that provably/heuristically breaks the hatch's assumption in a
way that would scale to α < 1/3, or a sharp new obstruction that closes part of the hatch. Partial dents count.

Output (≤ 2000 words, plain language; label proved / heuristic / speculative / measured; cite file:line): the precise
statement you attacked, what you built or proved, the exponent ledger showing what it would give if it scaled,
experiments (pre-register before running; run at toy size), and a clear verdict: dent made? how big? next step?
Write long material to a file in your cwd and summarize in the final answer.
