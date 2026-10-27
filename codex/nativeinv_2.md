# Adversarial round: nativeinv (native cost of factoring-equivalent invariants)

You are an adversarial referee. The repository is at the parent of your working directory (read-only; your
scratch dir is writable). Read `../../../notes/nativeinv_results.md` (the note under review) and, as needed,
`../../../src/nativeinv/trichotomy.py`, `../../../notes/factorbits_results.md` §3.2 (Theorem PT),
`../../../notes/nonsmooth_results.md` §5.2–5.3, `../../../notes/geometry_results.md` §4. Your own blind report is
`../../nativeinv_1.out` (it was merged; check that it was represented fairly).

Attack, concretely:
1. **Prop TRI** (§3.1), **Prop LOC** (§3.2), **Prop AUT** (§3.3): find any error, missing hypothesis, or overclaim.
   Is Prop TRI already in the literature (e.g. Bach–Miller–Shallit 1986; Bach–Shallit; Kunihiro–Koyama;
   Dieulefait–Urroz) in this generality? Is Prop AUT known (Allouche–Shallit; Coons 2010; Bell–Coons–Hare)?
2. **The inventory** (§1): wrong reductions, wrong costs, wrong "best native" entries, wrong complexity of the
   generalisation, misassigned barriers, missing factoring-equivalent invariants (give at most 10 that are genuinely
   new to the table, each with a 2-line reduction), and duplicates.
3. **Deep dives** (§2): for each item, is the stated "forcing step" really forcing, or is there a compression we
   missed? In particular: (a) the SL₃(Z/N) polylog presentation and the π₂/identities-among-relations heuristic;
   (b) the gcd-sum/Schur-ring argument; (c) the lens-space oriented/unoriented homotopy statements;
   (d) whether any small-complex homology computation could output φ(N) (your blind report said nothing rules it
   out — try hard to construct one or to give an obstruction).
4. **The meta-section** (§3.4): is the "five compressions" taxonomy complete? Is the claim that an N^c lower bound
   in any model containing Theorem PT would give constant-free permanent lower bounds correctly stated and
   attributed (see `../../../notes/codex_r2_generic_ring.md` §3)?
5. **Shortlist** (§4): for S1 (exact summatory functions H₁(x) = Σ_{d≤x} d⌊x/d⌋, Φ(x), twisted F_χ(x)), S2
   (digit automata on semiprimes for fixed-modulus fiber functions) and S3 (square-invariant class-E functions =
   Euclid-computable lens-space homotopy invariants): is each genuinely open, non-NP-hard in generalisation, and not
   ruled out? Is there a known Õ(x^{1/3}) algorithm for Σ_{n≤x} σ(n) (give a reference) — if so, state the resulting
   native N^{1/3} bound. For S2, can you prove or refute: "no k-automatic g agrees with φ(N) mod 3 on all semiprimes
   N = pq, 3 < p < q"? For S3, can you prove the conjecture for Dedekind sums at least?
   Would you rank a different item higher?

Output: a markdown report with a verdict line (ACCEPT / MINOR REVISION / MAJOR REVISION), then numbered findings
labelled FATAL / MAJOR / MINOR, each with a concrete fix. Keep any computation tiny. Label claims proved / known
(with reference) / conjectured / heuristic.
