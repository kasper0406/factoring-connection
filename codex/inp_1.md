# inp round 1: SUPPOSE factoring is in P (classical, deterministic or randomized). Where must the algorithm live?

Treat as a working hypothesis (unverified rumours; do not assume any specific claimed method). Use the repo's barrier
map to locate the escape hatches a polynomial-time classical factoring algorithm MUST use.
Repo: /Volumes/git/factoring-connection (do NOT edit; scratch only in your cwd). Read: README.md, BARRIERS.md,
notes/lowerbound_results.md (α ≥ 1/3 model theorem; hatches H / LP_exc^rel / R), notes/suffmap_results.md (Thm SS
size–entropy law; "every sufficient target for α < 1/3 is poly-time factoring in disguise"), notes/factorbits_results.md
(Theorem PT: pure methods learn nothing except by splitting), notes/leakmap_results.md (three leak families M/A/Q; NE1;
AC0-ε; open: PRIM, unions), notes/prim_results.md (open at deterministic depth ≥ 2), notes/slpcirc_results.md,
notes/floorhard_results.md, notes/topblock_results.md (Coppersmith known-bits), notes/periodfind_results.md (Thm T),
notes/native_results.md (55 FACT invariants; open: single a_f(N)), notes/selfred_results.md (SR1/SR2),
notes/polyimpl_results.md (fast factoring has no external math consequences), notes/completeness_audit_results.md,
notes/equivtarget_results.md, notes/entropy_results.md (leaks decay as N^{−1/4}), notes/guide_results.md,
notes/cogentic_results.md, notes/multcopp_results.md.
Answer (≤ 1800 words, plain language):
P1. A poly-time algorithm is far below L[1/3], so it cannot be a relation method under assumption H with the model's
    credit rules. List every proved repo theorem / model it must violate or step outside of, and say for each whether
    escaping it is "cheap" (the theorem's model is narrow) or "deep" (would need something genuinely new).
P2. Rank the 3–6 most plausible LOCATIONS for such an algorithm (e.g. non-pure computation mixing modular arithmetic
    with integer representatives / floors (PRIM depth ≥ 2, the open part of leakmap); a fourth leak family beyond
    group-order / archimedean / genus; a native polynomial-time evaluation of a factoring-equivalent invariant such as a
    single non-CM coefficient a_f(N); an SLP / short-circuit for a factorial-type quantity (Pochhammer, τ-conjecture
    side); a constructive smooth-number / beyond-Johnson CRT decoder (violating H); lattice methods beyond
    Coppersmith's known-bits barrier; non-commutative or higher-dimensional group structures; something else). For
    each: what it would have to compute, which repo barrier is closest, why that barrier does not already exclude it,
    and a "signature" (an intermediate result that would be checkable and would indicate the route is real).
P3. Consistency checks: what would a poly-time algorithm imply that we could sanity-check (e.g. poly-time DL in prime
    fields? the repo's DL-equivalences; computing φ(N), σ(N), a_f(N); breaking specific conjectures like the
    τ-conjecture / permanent lower bounds via factorial circuits (Bürgisser), hardcore-bit results)? Which locations
    are ruled out because they would prove something believed false?
P4. Your ranking with rough subjective weights (conditional on factoring ∈ P), and one concrete probe per top
    location that the repo could run cheaply.
