# buildbit: CONSTRUCTIVE attempt — assume a poly(log N) algorithm for the hidden-class bit exists, and build it

Your role is the **builder**, not the auditor. The user is ~90% sure a poly(log N) algorithm exists; the
orchestrator's prior is 10–15%. Your job is to produce the strongest constructive attempts. Kill a route only with a
concrete obstruction, and record what each attempt actually computes (its yield ledger). Repo is READ-ONLY for you;
work in your cwd. A sibling thread (codex/hidclass.md) does literature and hardness; don't duplicate it. Relevant
repo notes: notes/entropy_results.md (real CF(√7N) leak of a character of p, decaying ~N^{-1/4}), notes/kern_results.md
(cubic analogue ⇒ QR in a norm-one torus), notes/coef_results.md, notes/genusleak_results.md (ambiguous-form routes
reduce to group order or relations), notes/newfamily_results.md, notes/leakmap_results.md (leak families M/A/Q:
group order, archimedean, genus), BARRIERS.md.

## Targets, simplest first
- **T0 (Blum / sum of two squares).** N = pq ≡ 1 mod 4. Decide p ≡ q ≡ 1 mod 4 vs p ≡ q ≡ 3 mod 4. Equivalently
  r_2(N) = 4·Σ_{d|N} χ_{−4}(d) ∈ {16, 0}; equivalently χ_{−4}(p) + χ_{−4}(q) ∈ {2, −2}.
- **T1.** N = rq, r ≡ 1 mod 8, q ≡ 3 mod 8: decide r mod 16 ∈ {1, 9}. With A(z) = Σ_{n≥0} z^{n²},
  [r ≡ 9 mod 16] = [z^{N²}] A(z)(A(z^{128}) + A(z^{256})) mod 2 (verified).
- **T2.** N = pq ≡ 1 mod 7: both p, q split completely in the cubic field of conductor 7 (F = x³+x²−2x−1), or both inert.
An algorithm for any of these is a success. Note what generalises: T0 → T1 → full lift r mod 2^j (j ≈ ¼log₂N would
factor via Coppersmith).

## Routes to build (spend most effort on whichever looks most alive after a first pass)
R1. **Continued fractions / infrastructure of Q(√N), Q(√kN).** For T0: x² − Ny² = −1 can be solvable only if
    p ≡ q ≡ 1 mod 4; for p ≡ q ≡ 3 mod 4 one of px² − qy² = ±1, ±2 is solvable (the CF midpoint), and the period is even.
    The period is ~√N. Try to extract its parity or midpoint type without walking it: compact representations,
    distance arithmetic (Buchmann/Williams infrastructure), baby-step/giant-step on the cycle of reduced ideals,
    the ambiguous-ideal structure, the norm of the fundamental unit computed via something other than the regulator,
    or a short prefix or statistics of the CF of √(kN) for many chosen k aggregated across k (entropy found a per-k
    leak decaying like N^{-1/4}; can aggregating over poly many k, chosen adaptively, amplify it to constant advantage?
    Quantify the joint dependence, not just the marginal).
R2. **Class-group 2-structure of Q(√−N), Q(√−kN), Q(√N).** Genus theory gives the 2-rank; Rédei matrices give the
    4-rank from Legendre symbols (p/q), (−1/p), (2/p). Find a quantity that is computable in poly(log N) from N
    (random elements of the class group via forms of discriminant −4N, their orders' 2-parts? square roots of forms?
    composition with ambiguous forms? the Rédei symbol/spin?) and correlates with the hidden class. Ambiguous-form routes
    were killed in genusleak (reduce to group order or relations); state precisely how your attempt avoids that.
R3. **Algebraic tests in (Z/N)[i], (Z/N)[ζ_16], (Z/N)[x]/F.** Frobenius equals the Galois automorphism σ in the inert
    case and the identity in the split case, componentwise. Look for an identity in a, σ(a), N, and powers whose truth
    differs between the cases, without needing a group order (e.g. twisted Frobenius-pseudoprime tests, Lucas sequences
    with parameters tied to N, norms of N-th powers, Jacobi symbols of norms, the square class of a^{N±1}/σ(a)^…).
    kern reduced the cubic analogue to deciding squares in the norm-one torus: try to finish that step.
R4. **Your own route.** Analytic (twisted divisor sum Σ_{d|N} χ(d) via Voronoi / approximate functional equation /
    theta inversion of a single coefficient), CM curves (y² = x³ − x is supersingular at p ≡ 3 mod 4, so T0 is
    "both supersingular vs both ordinary"; is there an N-only supersingularity witness over Z/N?), quaternion orders,
    or anything else. One serious extra route is enough.

## Requirements
- For each attempt: an exact statement of what it computes in poly(log N), a proof or measurement of whether that
  quantity determines the target, and if not, the precise obstruction (which leak family, and whether it is a
  conservation law or a technique gap).
- **Experiments: pre-register first** (write prereg.md in your cwd with hypotheses, statistics, sample sizes and
  pass thresholds BEFORE running). 40–80-bit balanced semiprimes, ground truth from the factors, flag any gcd split.
  For a statistic with advantage, measure how it scales with bit size (constant vs decaying, and the exponent).
  A pass = advantage that is not decaying in N at 99% confidence across ≥3 bit sizes.
- If something passes, stop broadening: harden it (larger sizes, a proof of correctness, T0 → T1 → lift).
- Final: ≤ 2500 words. Per route: built / measured / killed (with obstruction), and the single most promising
  unfinished construction with a concrete next step. Label claims proved / known / heuristic / measured / speculative.
