# rankdegen_1 — blind round: degenerate linked modules, large boxes, split side

You are a research mathematician (geometry of numbers, Diophantine approximation, algebraic number theory, NFS
complexity). Read-only access to /Volumes/git/factoring-connection; scratch in your cwd; python available
(/Volumes/git/factoring-connection/.venv/bin/python has sympy, python-flint, gmpy2, numpy). Web search enabled.
Label every claim PROVED / PROVED (model) / HEURISTIC / CONJECTURE / REFUTED. Rigour over speed; a precise negative
answer is valuable. Do not claim novelty.

## Background (read these)
- notes/rank3_results.md IN FULL (Theorem R1, Lemma R2, Corollary R3, Prop R4, R5, §7 open list).
- notes/resbarrier_results.md §1, §2.3, §2.6, §4.2.
- notes/ptesplit_results.md §2.4; notes/ptebarrier_results.md §6.3 (O3).
- src/rank3/forms.py (exact Δ_F etc.) if you want to compute.

Summary of the setting. Number fields K, K′ (degrees d, d′, D = d + d′), orders O, O′ with ring maps φ, φ′ to Z/N
(N squarefree). A linked rank-r module: r pairs (u_j, u′_j) ∈ O × O′ with φ(u_j) = φ′(u′_j). Candidates a ∈ Z^r in a
box |a_j| ≤ E; F(a) = N_K(Σa_j u_j)·N_K′(Σa_j u′_j) = ∏_{i≤D} L_i(a), decomposable of degree D in r variables.
Δ_F = ∏ over ordered r-tuples of distinct factors of det(L_{i1..ir}). Rank3 proved: if Δ_F ≠ 0 then N^{P(r,D)} | Δ_F
and so (Bean–Thunder) V_F := vol{|F| ≤ 1} ≤ C_r N^{−r(r−1)/(D(D−1))}; with Davenport's lemma (error rD(2E+1)^{r−1})
this gives, when the box error is ≤ B/4 (side ≲ B^{1/(r−1)}), a two-norm ledger
log X ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r), α ≥ 1/3, and constants increasing in r (rank 2 is best).

## Tasks
1. **Degenerate modules (Δ_F = 0).** Classify the linked modules in which some r of the D linear factors are
   dependent. At least: (i) modules lying in a proper subfield (all u_j, or all u′_j, in a subfield; F a power of a
   lower-degree norm form); (ii) repeated or conjugate rows, including linkings that hold in characteristic 0
   (cross rows proportional); (iii) modules whose linear factors span, on some r-subset, a space of dimension < r
   ("lower effective rank" on a sublattice; Schmidt's degenerate subspaces). For each type, prove a ledger (norm-size
   lower bound for many small-|F| candidates in a box) by reduction to lower rank or to a subfield, OR exhibit a
   degenerate family that beats the rank-2 ledger, and cost it in the L[1/3, c] size model of rank3 §4(b). Is there a
   general invariant (e.g. a nonzero semi-discriminant, the product of all NONZERO r×r minors, or minors of a
   maximal GP subconfiguration) that still carries the N-divisibility? Does "Δ_F = 0" with useful linking force the
   characteristic-0 coincidence that makes relations trivial (e.g. both sides' maps agreeing on a common subfield)?
2. **Large boxes, side > B^{1/(r−1)}.** Davenport's error (2E+1)^{r−1} then exceeds the factor base. Is there a
   lattice-point count for {a ∈ Z^r ∩ [−E,E]^r : |F(a)| ≤ X} (F decomposable, GP, N^{P} | Δ_F) of the form
   ≤ L^{o(1)}·V_F X^{r/D} + (lower-rank sublattice contributions), valid in the regime that matters (3-dim
   power-basis sieving at its size-model optimum e = 2β/3; tiny boxes a_j ∈ {0, ±1} with large r)? Tools:
   Davenport; Barroero–Widmer (lattice points in o-minimal definable families; error governed by volumes of
   projections); Thunder's counts for decomposable forms (Ann. Math. 2001; Compositio 2005, with m(F), a′(F));
   Schmidt subspace theorem for the degenerate directions. Prove the count, or show exactly where it fails with an
   explicit example (a GP linked module whose box has far more points with |F| ≤ X than volume predicts, in a way
   that would matter at the L[1/3] scale). Note: the excess must be compared to the factor base B, not to 1.
3. **Split side (ptesplit O3).** When the uncharged side's size log X₁ ≳ (r−1)x/(D−1), no bound on the charged norm
   alone follows from linking (rank3 §6 GP example: K = K′ = Q(ζ₅), generators (1, ζ₅+N, (ζ₅+N)²) vs (1, ζ₅, ζ₅²)).
   Decide whether this regime is genuinely exploitable by an algorithm — i.e. can the split side be split *by
   construction* at sizes far above the (r−1)x/(D−1) threshold while the charged side stays small — and cost it;
   or show it is excluded by another argument (e.g. counting of split/smooth values on the split side: a split
   element of an order of degree d with norm ≤ X₁ built from pieces of a factor base; Theorem E / Lemma ND
   constraints; the fact that a split-by-construction side is just a rational side in disguise).

Give the strongest PROVED statements you can, with complete proofs, the explicit examples, and a final verdict per
task. Under ~4000 words.
