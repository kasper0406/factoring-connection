You are an adversarial referee (computational and analytic number theory). Repository root: /Volumes/git/factoring-connection (read-only for you; scratch in your cwd or /tmp). Read notes/detGRH_results.md §§1–6 (the draft; §§7–10 are placeholders) and, for context, notes/thmA_proof.md, notes/nfsfruit_results.md §3, notes/lemmaM_results.md §1 and §5 (B1), notes/dlog_results.md §4 (Prop XK). Do NOT open codex/detGRH_1.out.

Your job is to find errors. For EACH numbered claim below give a verdict CORRECT / CORRECT WITH FIXES / GAP / WRONG, with the precise fix or a counterexample (check small cases by computation where possible; python with sympy is available in /Volumes/git/factoring-connection/.venv/bin/python).

1. Lemma MR (a) and (b), including uniqueness and the mod ±1 version; Corollary MR (the GRH step via Bach).
2. Proposition FL (i)–(iii): the isomorphism Phi_full, the criterion dim(V_c ∩ W') >= 2, the codimension bound c <= dim W' − 2.
3. Proposition CD (a), (b): certificate dependence; "y descends to a function of R iff G_B ∩ mu_2 ⊆ {±1}".
4. Theorem M: case (S) (every bullet, including the Miller-in-G_B claim and the B1 instance) and case (G) (honest certificates ⇒ y∘c not multiplicative mod ±1 ⇒ a fruitful cocycle exists). Is (G) correct as stated? Is it vacuous? Is Lemma SC right?
5. Proposition Z (Miller + Bach): the statement of Bach's bound used (every proper subgroup of (Z/m)^* misses a prime < 2 log^2 m under ERH — exact hypotheses/constants?), the proof that the subgroup H is proper and that every a not in H yields a factor, and the claim lambda(N) | [Gamma : Lambda_c].
6. Proposition EQ (equivalence of deterministic factoring with deterministic full-rank kernel-relation production, under ERH). Any hidden circularity or missing case (prime powers, small factors, N even)?
7. Consequences Z1–Z5, especially Z3 (no structural Q-rank obstruction) and Z4 (QS relations are not kernel relations; Z route does not apply to QS).
8. §5 table: is every randomness source of the Theorem A algorithm listed? Is each derandomization claim right (deterministic elimination cost and the slack arithmetic 3a*/c* = 0.879 in RA and 3a < (1+theta)/sqrt2 iff theta < 1/2; deterministic character primes under GRH; deterministic recursion; the Prop XK step that picks a random prime q)? Is "irreducible" justified for rows 1–4a?
9. Observation FS (family-size floor for the moment method) and Observation CS (the Hölder computation: error/main <= |E|^{-1/2} N^{1/2+o(1)} with the GRH sup bound, N^{3/4} without). Recompute CS. Are these fairly labelled?
10. The verdict: is anything claimed stronger than what is proved? Is Proposition EQ / the Z route known in the literature (give citations if so: Miller 1976, Bach 1984 "Discrete logarithms and factoring", McCurley, Pomerance, Shoup, Bach–Shallit book, Crandall–Pomerance)? Is a deterministic subexponential factoring algorithm under GRH known (for all or almost all N)? Use web search.

Write a markdown report with one section per claim and a final table of verdicts.
