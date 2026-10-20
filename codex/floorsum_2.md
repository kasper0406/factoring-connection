# floorsum_2: ADVERSARIAL review of notes/floorsum_results.md

You are a hostile but fair referee (computational number theory, lattice-point geometry, harmonic
analysis, lattice reduction) with web search. The repository at /Volumes/git/factoring-connection is
read-only; write scratch files only in your current directory. Read notes/floorsum_results.md in
full, notes/floorsum_prereg.md, src/floorsum/{farey,ladder_ops,lattice_block}.py, and
notes/analytic_results.md §3 (Theorem L_j). Your own blind report is codex/floorsum_1.out.

Attack every claim. For each numbered item give a verdict ACCEPT / ACCEPT WITH FIXES / REJECT, with
the precise defect and a fix. Check proofs line by line; check constants; try small counterexamples
by computation where possible (python3 with fractions is fine).

1. Lemma V (even j) and its use with Balog–Bárány: is S_b ⊆ G right, is the candidate exception
   right (including divisors), is the closed/open hull switch handled correctly?
2. Theorem T (convex-access model): is R_v = conv(L∖{v}) legitimate; does each LP query rule out at
   most one alternative; is the model statement precise; is the conclusion "Ω(N^{1/3} log N) for
   every even j" correctly scoped (the ladder needs F(P;0,H) per block — does it, or only the
   candidate-count difference)?
3. Prop PL: is "the lower hull of a linear staircase has O(log) vertices" correctly justified by
   Cook–Hartmann–Kannan–McDiarmid? Is the bound V_b/O(log) right?
4. Theorem O (Farey oracle): exactness (Lemma R, monotone runs, the split at the extremum) and the
   cost bound (window count, variation bound, |qs − p| < 1/Q). Look for inputs where the code is
   wrong or slow.
5. Prop F: Bernstein, the Remez inequality (Erdélyi 1992: exact statement and normalisation), the
   2K+2 interval count, the max-gap counting, the 0.099/√g − 3 optimisation, and the model's scope.
   Is the "≥ N^{1/4+1/(2(j+2))}" consequence correctly qualified?
6. Prop A (affine arc length invariance) and the claim that it "explains" why Euclid fails: is this a
   theorem or a heuristic? Fix the wording.
7. §6: the high-bits reformulation; the Coppersmith containment H ≤ N^{α²} for j = 2 (constants!);
   the j = 3 threshold α ≈ 0.276; the claim that beating (1/4)log₂N high bits is a standing open
   problem — is that accurate (cite precisely)? Is E5's lattice implementation correct (the exact
   sufficient inequality, root extraction)?
8. Prop D (D(N) via the ladder): correctness of the correction term for even and odd j; is the claim
   "below every known factorisation-free algorithm" accurate? Survey again: any exact D(N) or
   restricted divisor-sum algorithm below N^{1/3} without factoring?
9. The pre-registration handling: E2's failure, the post-hoc E2b/E4b, the E5 amendment with a
   disclosed peek. Is anything overclaimed?
10. MOST IMPORTANT: try hard to REFUTE the verdict. Is there an exact oracle for the ladder's actual
    queries (fixed N, shifts {0, −ρ−η}, dyadic subintervals) costing below H^{1/3} for j = 2, or below
    H^{1/10} for j = 4, by any means (arithmetic of N, lattice reduction, amortisation across blocks,
    precomputation shared by all blocks, a different decomposition)? Even a total cost below N^{1/3}
    for the whole ladder (amortised) would matter. Report the best idea you find with a cost analysis.

Deliverable: a referee report (≤ 3500 words) in your final message: overall verdict first, then the
numbered items.
