# smoothenum_2 — adversarial referee of notes/smoothenum_results.md

You are a hostile but fair referee (number theory / algorithms). Read-only repo at
/Volumes/git/factoring-connection; scratch in your cwd; python allowed; web search enabled.

Read notes/smoothenum_results.md in full, plus src/smoothenum/{target,barrier,mech,run_proto}.py and
notes/smoothenum_prereg.md. The note claims NO constant improvement, but it makes several barrier claims and
constant-table claims. Your job: break them.

For each item give verdict CORRECT / WRONG / GAP / OVERSTATED with a short justification and, if wrong, the fix.
1. Table 1 (T1–T8): recompute every constant independently (write your own minimiser). In particular T6
   (algebraic-listing infimum (4/9)^{1/3}) and T7/T8 (doubly listing → 0 / polynomial time): are the
   feasibility conditions right (δ < 6β_a; e → ∞ within the L[1/3] model; the degree-1 two-linear-forms
   argument with B = (log N)^κ and box N^K)? Is T5's σ-slope right?
2. Proposition LQ (§3.2): is the proof correct? Is (C u log u)^v the right bound for ρ(u−v)/ρ(u) uniformly
   for v ≤ v_box? Is assumption H_sel stated precisely enough, and is it the real crux? Construct, if you can,
   an algorithm inside the LQ access model (or a mild extension) that violates the bound — e.g. selecting
   nonempty lattices among M > area without paying one unit per M, or exploiting overlap of many lattices.
   Is the finite-size g* table (σ_eff, "76% of the rational penalty at 1024 bits") computed and caveated
   correctly (multiplicity, distinct yield)?
3. Prop GEN, Prop W (k-tree inequality (1+lg k)/k ≥ log(m/E)/log(Em)), Prop BC (Boneh decoding bound use,
   the S^{−c} hit rate with c ≥ 1/2 − o(1)): check each derivation. Is Boneh's bound quoted correctly
   (D, H, S, P)? Could a multivariate (two-variable a, b) Coppersmith/Boneh variant over the whole box do
   better than row-by-row?
4. Prop FT and §3.7 (adaptive special-q with bounded tails): correct? Any loophole (e.g. unbounded number of
   tail primes that are nonetheless cycle-eliminable at L^{o(1)} cost)?
5. Experiments (§4): are the registered predictions honestly evaluated against the prereg? Is the P2 failure
   correctly reported? Is the las comparison fair (regions differ; threads; cost units)? Any claim in the
   verdict not supported by the data?
6. Literature: is anything here (LQ-type barrier, factory per-N = rational-output program, two-sided special-q,
   Boneh CRT smooth-number search, Wagner) already published in a form the note should cite or credit?
7. Bottom line: list the fixes the note must make, ordered by severity.
Keep it under ~3000 words.
