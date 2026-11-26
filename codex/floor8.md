# floor8: does a 2D floor sum with a quadratic congruence N | Q_0(x,y) admit a fast mod-8 evaluation without splitting?

Precise target. Fix a nondegenerate binary quadratic form Q_0(x,y) = ax^2 + bxy + cy^2 (disc D = b^2 - 4ac,
D fixed, e.g. D = -23 as in the kern/coef work). For a box B = [0,X)^2 (or a hyperbola-type region), consider the
2D floor / lattice-point sum
    S(N) = #{ (x,y) in B : N | Q_0(x,y) }     (the count form),
and also a genuine floor variant S'(N) = sum_{(x,y) in B, N | Q_0} floor(L(x,y)) for a linear L.
Question: can S(N) mod 8 (resp. S'(N) mod 8) be computed in poly(log N) time WITHOUT "splitting the congruence" —
i.e. without factoring N and without decomposing {Q_0 ≡ 0 mod N} into its prime-ideal sublattices (which needs sqrt(D)
mod N / factoring)?

Context in the repo (/Volumes/git/factoring-connection; do NOT edit; scratch in your cwd): notes/kern_results.md and
notes/coef_results.md (u_2(N) ⟺ R(N) = #{s^2+st+6t^2=N} = 4 vs 8 is the factoring-hard bit; genus characters are the
computable abelian shadow; the finer 2-adic/mod-8 structure of representation numbers is the class-group 2-part),
notes/floorsum_results.md, notes/floorhard_results.md (general floor oracle NP/#P-hard; balanced top blocks open),
notes/analytic_results.md (floor-sum ladder), notes/topblock_results.md (Coppersmith known-bits), notes/classE_results.md
(Rademacher integer Phi(a/N), cotangent sums, Manin symbols mod small moduli), notes/leakmap_results.md (genus family Q).

Orchestrator's hypothesis to test (be adversarial, confirm or break it):
(H1) A nontrivial mod-8 SHADOW is fast: the quadratic-Gauss-sum root number / Weil index (an 8th root of unity,
     computable from Q_0 and N mod 8 by Landsberg–Schaar / Weil reciprocity) and the genus character (disc | N) are
     both poly(log N), no splitting. [likely PROVED/KNOWN]
(H2) The FULL count S(N) mod 8 is factoring-hard for generic N, because S(N) mod 8 contains the u_2-type / R(N)=4-vs-8
     bit = the class-group 2-part = QR-in-a-torus. [predicted]
(H3) The one real uncertainty is the genuine FLOOR variant S'(N) mod 8: does 2D Dedekind–Rademacher reciprocity (the
     floor, not the exponential sum) stay in the computable regime mod 8 (classE found Rademacher/cotangent invariants
     mod small moduli)? Decide this precisely.

Tasks (<= 2200 words; label proved / known / heuristic / measured / speculative; cite file:line / refs):
F1. Set up the orthogonality decomposition S(N) = (r/N) vol(B) + sum_{t != 0} (box factor) G_t(Q_0; N), state exactly
    which pieces are computable without splitting: the 2D quadratic Gauss sum G_t(Q_0; N) and its reciprocity (give the
    explicit evaluation / Weil index, and confirm it needs only N mod small powers of 2, (disc|N), gcd data — NOT a
    factorization). Where does r = #{Q_0 ≡ 0 mod N} enter, and exactly which part of r mod 8 is the genus (computable)
    shadow vs the class-group-2-part (hard)? Make H1 precise and H2 precise: reduce "S(N) mod 8 computable ⇒ ?" — does
    it compute R(N) mod 8 / u_2(N), and would that factor N or solve QR-in-the-torus (notes/kern_results.md)?
F2. The floor variant S'(N) and 2D Dedekind–Rademacher reciprocity: does the Euclidean/reciprocity reduction that makes
    1D linear floor sums O(log N) extend to the quadratic-congruence-constrained 2D sum mod 8? Give the reciprocity law
    if it exists (Dedekind–Rademacher, higher Dedekind sums à la Zagier, cotangent sums), the exact 2-adic denominators
    (why mod 8 specifically — the sawtooth values are in (1/2)Z, Dedekind sums have denominator dividing 12k, so mod-8
    control has specific obstructions), and whether evaluating it still needs the splitting. Confirm or break H3.
F3. Run a decisive pre-registered experiment (40–80 bit N = pq, D = -23 and one other D; exact ground truth from
    factors): compute S(N) mod 8 and S'(N) mod 8 directly, and compute every candidate fast formula (Gauss/Weil-index
    phase, genus character product, Dedekind–Rademacher reciprocity value, any closed form you derive). Measure exactly
    where fast formula = truth and where it diverges; in particular test whether (truth mod 8) determines / is
    determined by R(N) = 4 vs 8 (the hard bit), and whether any poly(log N) statistic predicts it beyond the genus
    shadow (held-out, planted-signal control, amplification not assumed).
F4. Verdict: precisely which mod-m (m | 8) evaluations are fast (proved), which contain the factoring-hard bit, the
    sharpest open sub-statement, and probability that a fast full mod-8 evaluation exists.
