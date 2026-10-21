# Adversarial round: break the negative analysis of "Prop SI constructors"

You are an adversarial referee (computational number theory, lattices, knapsack algorithms). The repository is
read-only for you; it is at /Volumes/git/factoring-connection. Read `notes/siconstruct_results.md` in full (and, as
needed, `notes/suffmap_results.md` §1.4 (Prop SI), `notes/schnorrcvp_results.md` §2 (Theorems P, G), and the code in
`src/siconstruct/`). You may write and run Python in your scratch directory. Label every finding FATAL / MAJOR /
MINOR, with a precise pointer (section, item) and a concrete fix.

The note claims that no constructor of y-smooth integers in I_k = [kN, kN + h] (y = L[1/3, κ], h = L[2/3, ν],
ν > 2/(3κ)) reaches L^{o(1)}, or even L[1/3, 0.757], per output, and supports this with:

1. **Prop CV** (§2.1): squarefree Prop SI outputs are CVP_{1+O(1/x)} answers in Adleman's ℓ2 prime lattice at
   scale C = kN/h, and conversely; λ₁² = (2+o(1))x so the target sits at ≈ λ₁/√2; the slack formula c(t); the
   translation γ − 1 = O(ln ln n/ln n) in dimension n.
2. **Lemma SR / Lemma RC** (§2.2): scale invariance, and "a constructed factor s plus a tested cofactor forces the
   cofactor ≥ N/(K(h+1))", hence L[2/3] per output.
3. **Lemma MA** (§2.3): N-oblivious two-list merging gives |A||B|(h+1)/N hits for any number K of multipliers served,
   so per output ≥ 2√(N/(R(h+1))); "amortisation over k is real but worthless".
4. §2.4 (k-lists, Wagner, representations), §3 table, §4 approach-by-approach analysis (including Boneh/Coppersmith:
   certified divisor d ≥ h, so Lemma RC applies), §5 total costs, and the experiments §6 (including the two disclosed
   FAILs P2c and P3a and their post-hoc corrections).

Tasks:
- (a) Check every proof in §2.1–§2.2 line by line. Are the inequalities right (signs, the role of v ≥ 1, the
  squarefree assumption, the 1/(4C²) term, the λ₁ lower bound proof, the O(√x) and O(1) factors in the converse)?
- (b) Attack the *scope* of Lemma RC and Lemma MA. Find a construction that escapes them: e.g. a constructed factor
  that is aimed *non-obliviously* (depending on N mod s) without solving the original problem; several constructed
  factors aimed jointly; exact or approximate counting of smooth numbers in intervals plus binary search; number-field
  norms; bivariate / whole-box Coppersmith; using k with special structure; recursion with unequal splits. For each,
  either show it reduces to an analysed case or give its cost.
- (c) Is the claim "one-sided is harder, not easier, than two-sided in oracle terms" correct and fairly stated?
- (d) Is the amortisation analysis right, including in the Prop SI regime Kh ≪ √N (merge cost independent of K)?
  Is there any way for one list (or a precomputed structure) to serve many k or many N so that the per-output cost
  drops below √(N/(Rh)) without exponential precomputation?
- (e) Check the cost table numbers (c(τ), the 0.757 threshold, 6^{1/3}, c(t) values, L[1/3, 2ν/(3κ)] gain, Boneh
  analysis), and whether the experiments' interpretation (toy h = N^{0.6+}) is honest.
- (f) Literature: anything that already states Prop CV, Lemma RC/MA, or a constructor for smooth numbers in
  [X, X + L_X[2/3]] at y = L_X[1/3]?

End with a list of required changes, most important first, and a one-paragraph verdict.
