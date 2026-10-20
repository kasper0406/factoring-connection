# lv5unif_2 — adversarial referee: Theorem U (rigorous L[2/5] factoring, conditional on the internal chain)

You are an adversarial referee (analytic number theory + algorithms). Work read-only in
/Volumes/git/factoring-connection; scratch files only in your working directory. Use web search to check
citations where useful.

Referee **notes/lv5unif_results.md** (the claim: Theorem U, a Las Vegas algorithm factoring every N in expected
time L_N[2/5, c], c finite, conditional on the project's internal fixed-degree chain). The base texts it modifies:
paper2/main.tex (the preprint of Theorem A; lemma labels are cited as in the .tex), notes/lem85_results.md
(Lemma 8.5′), notes/recurse_results.md (Lemma S″, Lemma F, Lemma OS, T_rec). Background: notes/thmA_proof.md,
notes/lv6_results.md. An earlier blind round (codex/lv5unif_1.out) reached the same verdict independently; do not
defer to it.

Your job is to find errors. In particular, attack:
1. **Uniformity in (d, k).** Is every constant in §§2–9 genuinely independent of d, k, N (polynomial factors in d
   allowed, since d ≤ log N)? Hunt for any remaining c^d, d!, C(2k,k), (log N)^{O(d)} or (log N)^{O(k)} factor that
   enters the *relative* moment discrepancy (E N)²/E[N(N−1)], rather than only the cost. Check §6 Step 7's power
   count and the claimed exact cancellation, with the common box radius V and the b-window ((1−1/d)V, V].
2. **Lemma Kern-U** (§3): is the kernel basis claim right, is the triangular count right, and is the hypothesis
   bb′ ≤ (2H+1)/d sufficient for the probability form?
3. **Lemma CS** (§2): the computable S via a Rankin value at an approximate saddle point, the use of
   Hildebrand–Tenenbaum's Theorem 1 (uniformity range!), the m bounds, and the circularity S ↔ X_F(S).
4. **Lemma 8.5′(R′)** (§5): does every step of notes/lem85_results.md really use only (R′1)–(R′4)? Is
   Drappeau's Lemme 5 applicable in this range (y = L_N[2/5], x = N^θ, θ ≍ D^{−1/5})?
5. **Lemma GM-U and FM-U** (§4): the layer argument with the new envelope X_F = (2d+1)MHV^d; the ϱ comparisons;
   the bad-mass step; #𝒫.
6. **LD*-U, L2-U, W-U, DS-U** (§6): especially the cushion change in R₀ and class (B)/(C) margins.
7. **Theorem K-U / Proposition select-U** (§8): Lemma S″'s hypotheses; the sharper (ii); is anything lost to 2^ω?
8. **LV6-U** (§9): the factorial-bit character primes; Stark's Lemma 8 constant; LO's hypotheses; the cost.
9. **Lemma F_α and the cost assembly** (§10): the recursion inequality (10.1), every term of P(x) and D(x),
   and whether the shape-test count and the split cost are correctly placed.
10. **Proposition MF** (§11): is the pointwise-plus-independence argument valid? Is the multi-field
    post-processing sketch sound (exceptional classes, (T))?
11. **The parameter choices** (§1, Lemma P): are all margins satisfiable simultaneously? Check
    src/lv5unif/params.py's margin list for missing comparisons.
12. **Logic of the overall claim**: does Theorem U follow from the listed ingredients? Is anything used at fixed d
    in the sources that is not re-checked here?

For each finding give: severity (FATAL / MAJOR / MINOR), location, the precise problem, and a fix if one exists.
End with a verdict on Theorem U (ACCEPT / ACCEPT WITH FIXES / REJECT) and the weakest point.
