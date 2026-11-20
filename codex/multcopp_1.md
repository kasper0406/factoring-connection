# multcopp: can the RH-remark's "multiplicity over Z" become a Coppersmith-style factoring algorithm below L[1/3]?

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch work only in your cwd).
The user's only goal is factoring general N in L_N[α] with α < 1/3 (heuristic is fine). Constant gains are worthless.

## The question
Coppersmith / Howgrave-Graham get their power from auxiliary polynomials vanishing to high order ("multiplicity")
at a hidden root modulo p (via powers N^j, shifts x^i f^j), then lattice reduction. That is a Stepanov-type
argument. Its barrier is the known-bits limit (≈ half the bits of p, i.e. n/4 bits of N, are needed). This repo's
floor-sum ladder (non-relational route) hits exactly that barrier: notes/topblock_results.md (Prop T),
notes/floorsum_results.md, notes/floorhard_results.md, notes/analytic_results.md; see also notes/codex_r8_copp_bridge.md,
notes/below16_results.md, notes/unbalanced_results.md, notes/completeness_audit_results.md §7 Q4 ("every leak
delivers O(1) bits; n/4 bits needed"), notes/suffmap_results.md (Thm SS).

A note (untrusted, unverified; text at /Volumes/git/factoring-connection/codex/scratch/remark.txt) claims a proof of
RH + Dirichlet GRH. Its key "new" idea: Stepanov's function-field method needs multiplicity, which over Z is "like
differentiating whole numbers"; they obtain multiplicity *natively automorphically* by chaining ~√X iterated theta
lifts U(X_0)→…→U(X_m) between totally definite unitary groups (π ↦ π ⊞ μ_j[t]), each lift "winning" one order of
multiplicity; a "Stepanov determinant"; multiplicity 2 on O(2,n) via a natively-automorphic determinant; divisibility
of iterated lifts proved via Chowla–Selberg, Weyl, Siegel mass; D_p ≡ 0 on a CM curve for p inert in K
(D the invariant derivation, D^a θ(0)Ω^a divisible by p^{a_p + a_p^3 + … − O(log a)}); "cost of each lift … one would
stop at 2/3 a priori".

The claimed Lean formalization is unpacked (sources only; DO NOT build or run anything from it) at
/private/tmp/claude-501/-Volumes-git/c1686473-142d-57cc-8202-223c6b26715f/scratchpad/leansrc/lean/src . Main.lean / Statements.lean give the top statements (RH and Dirichlet GRH for all moduli).
Large trees: Development/ (1 GB; e.g. Development/WeilBound, Development/AutomorphicArithmetic/ThetaCorrespondence),
PositivityTheory/, Formalism/, Supplements/Kernels/MainTheoremsAll.lean. Use rg/grep to locate the multiplicity
mechanism (e.g. rg -il 'ThetaLift|multiplicit|Stepanov|determinant|derivation' with -g '*.lean'); read the key
definitions/lemma statements, not whole files. It is huge: be targeted.

## Answer (≤ 1800 words, plain language; label proved / KNOWN / heuristic / speculative)
A1. What exactly is the multiplicity mechanism (from the remark and, where you can find it, the Lean statements: give
    file paths and lemma names)? What object vanishes to what order at what, and what is the analogue of the
    "auxiliary polynomial" and of the "point" being counted?
A2. Map it onto Coppersmith: in Coppersmith the multiplicity is p-adic (p^m | g(x0)) at a hidden root. Is the remark's
    multiplicity p-adic at a prime (cf. "D^a θ(0)Ω^a divisible by p^{…}" for p inert), archimedean, or spectral?
    Could it be made to vanish to high order at a HIDDEN p | N without knowing p (the essential Coppersmith feature)?
A3. Barrier check: does any version give more than Coppersmith's p^{1/2}-type lattice bound, i.e. recover p from
    fewer than ~n/4 known bits, or produce a polynomial-size auxiliary object with multiplicity m at p and small
    archimedean size beyond the Howgrave-Graham determinant tradeoff? Do the determinant/height ledger. Explain the
    "2/3 a priori" cost remark if you can.
A4. Cost: theta-lift chains of length ~√X — can any piece be computed in time L[α<1/3] for X related to N, or is it
    inherently a proof device? Is there a finite-dimensional shadow (e.g. theta series coefficients, CM values,
    Gauss-sum/Jacobi-sum identities, derivatives of theta at CM points) computable mod N that carries multiplicity at p?
A5. Verdict: concrete algorithm candidate(s) with ledger and the single missing lemma, or a clean kill (which repo
    barrier it lands in). P(α<1/3, heuristic). One cheap toy test (pre-register it before running, then run it if
    feasible at toy size).
