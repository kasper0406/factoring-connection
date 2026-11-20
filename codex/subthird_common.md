# subthird: find a mechanism that beats L_N[1/3] for factoring (common preamble)

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit files there; write scratch work only in
your current working directory). The ONLY goal: an algorithm (heuristic is fine) for factoring general N in
L_N[α] with α < 1/3. Constant-factor improvements at α = 1/3 are worthless here. Analysis-only results count only if
they are a concrete step toward such an algorithm (or a clean kill of a candidate).

Read first: README.md, BARRIERS.md, notes/completeness_audit_results.md (§0, §5, §7), notes/suffmap_results.md
(Thm SS, the size–entropy law), notes/resbarrier_results.md (§2.6, §4.2, Lemma RES, Cor RB2, H5), notes/rank3_results.md,
notes/rankdegen_results.md, notes/midk_results.md §0, notes/grdeg_results.md (§0, Prop RR, the L[1/4] ledger),
notes/smoothbuild_results.md, notes/sihard_results.md (Thm J), notes/heurnfs_results.md. APPROACHES.md §0 indexes the rest.

Known facts you must not rediscover: every barrier in the repo is a MODEL theorem; the size–entropy law says α<1/3
needs candidates whose smoothness target is far below the resultant floor; rank-2 selection rules are closed by Cor RB2;
the aggregate rank-r≥3 linked-module family formally permits α = 1/4 if the r−1 linking conditions fail to be
independent "favourably" (G_r failure), and grdeg Prop RR closes relative-rank-2 / power-basis subfield cases; the
open crack is type-D modules over a fixed small subfield k with relative rank m ≥ 3 and non-power generators.
The repo has NO theorem bounding (i) adaptive constructors that build candidate and smoothness certificate jointly,
(ii) H5 sparse selection from boxes larger than B, (iii) compositional/hybrid methods across families.

Output discipline (≤ 2000 words, plain language):
- For each candidate mechanism: (1) one-paragraph description; (2) exponent ledger: candidate size log X, smoothness
  bound B = L[β], number of candidates, cost → α, with the x = log N bookkeeping explicit; (3) which repo barrier it
  must evade and the precise reason it is outside that barrier's model (file:line); (4) the single assumption it rests
  on, labelled proved / heuristic / speculative; (5) the cheapest falsification test (toy sizes) and its result if you ran it.
- Kill your own weak ideas explicitly ("killed by …"). Survivors only in the final ranking with P(α<1/3, heuristic).
- Pre-register any experiment (hypothesis + pass/fail criterion) in your scratch dir BEFORE running it.
