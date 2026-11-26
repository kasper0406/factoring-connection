# wheelgap: constructive short square-coset representatives — consecutive small-prime wheel vs random-prime wheel

Repo is READ-ONLY for you; work and scratch files go in your cwd. Relevant repo notes: notes/sihard_results.md
(Theorem J: oblivious constructors cost L[2/3] per output), notes/siconstruct_results.md (Prop CV: Prop SI is
promise-CVP in Micciancio's prime lattice; MITM does not amortise), notes/schnorrcvp_results.md, notes/suffmap_results.md
(Prop SI, size–entropy law), notes/lowerbound_results.md (MODEL theorem α ≥ 1/3; escape hatch "violate H" =
decoders/constructors), notes/escapeR2_prereg.md (square-coset re-split obstruction), notes/smoothbuild_results.md,
notes/equivtarget_results.md, BARRIERS.md.

## The brief (from an outside research summary; every formula was lost in transit, marked [?])
"In the conductor-recovery / square-coset formulation, let W be a wheel of small primes ([?]). Finding a small
positive integer h (or rational a/b with [?]) matching a prescribed Legendre sign vector mod W is unconditionally
equivalent to randomized factoring (for random-prime wheels, and conjectured for small-prime wheels).
Core question: can we algebraically construct smooth integers (or ideals in [?]) that simultaneously satisfy a
small-height congruence mod N (or land in a prescribed square coset mod W with height [?]), replacing sieving?
Dead ends: binary partition forms [?] and [?] — solving the F_2 linear system over prime partitions of [?] gives
exact coset matches, but continued-fraction cancellation bottoms out at [?]; for a 168-bit wheel this stalls at
~78–80 bits, far above the ~22-bit guaranteed recovery range. Explicit CRT wheel / BSGS difference tables: capped at
N^{1/4} exp(−(log 2/4 + o(1)) log N / log log N). Random-prime wheel decoding ≡ factoring. Genus-class sampling only
narrows an ideal to a coset of [?].
Direction 2 (this thread): test whether the consecutive small-prime wheel W = ∏_{ℓ≤y} ℓ admits a constructive
short-representative algorithm that fails for a random-prime product of the same bit length.
Success: a constructive method producing coset representatives of height [?] in time [?]. Stop: a proof that it
can't."

## Tasks
W1. **Reconstruct the formulation.** State the most plausible precise problem consistent with every surviving
    fragment (the 22-bit vs 78–80-bit numbers for a 168-bit wheel, the N^{1/4}·exp(...) difference-table cap, "conductor
    recovery", "rational a/b", genus-class sampling). In particular: what is the sign vector a function of (N? the
    hidden p?), why does a short representative factor N, and what height is needed. Mark it as a reconstruction;
    if several are plausible, pick the one closest to an actual factoring reduction and say why.
W2. **Literature.** This looks like the pseudosquare / prescribed-Legendre-pattern search problem: Lehmer sieve
    machines, Lukes–Patterson–Williams pseudosquares, Bernstein "Doubly focused enumeration of locally square
    polynomial values" (2004), Sorenson's pseudosquare work, Hall/Burgess/Graham–Ringrose bounds on least non-residues
    and sign patterns, Peralta/Davenport on runs of residues. What is the best known constructive cost for the
    smallest h with prescribed (h/ℓ) for all ℓ ≤ y, for consecutive primes vs arbitrary primes? Does the difference
    table cap match doubly focused enumeration? Search the web.
W3. **The primorial-specific structure.** For the consecutive wheel, the signs of products of small primes are
    F_2-linear in exponent vectors via the Legendre matrix (q/ℓ), q, ℓ ≤ y, which quadratic reciprocity makes nearly
    symmetric (Paley-like). For a random-prime wheel the same linearity holds, but with an unstructured matrix.
    Does the reciprocity structure (or the Chinese-remainder structure, or the fact that the wheel primes are also
    the smoothness base) give a constructive short representative, e.g. a lattice of exponent vectors with an F_2
    kernel condition and a height objective solvable by LLL or BKZ with provable quality? Prove a separation, or prove
    that the primorial wheel is no easier, up to L^{o(1)} or polynomial factors, in a model that covers these methods
    (Thm J / Prop CV may already do this; say exactly).
W4. **Experiment (pre-register first: write prereg.md in your cwd with hypotheses, metric, sizes and the threshold
    BEFORE running).** Wheels of 40–168 bits: consecutive primes vs random primes of matched bit length, random sign
    vectors. Compare the heights achieved per unit time by (i) exhaustive or doubly-focused search, (ii) the F_2 kernel
    plus LLL/BKZ exponent-lattice method, (iii) any method suggested by W3. Report the height exponent vs wheel size
    for both wheel types, with CIs. Use fpylll if installable in your cwd venv; otherwise flint/your own LLL.
W5. **Verdict.** Proved / known / refuted / measured; is there a primorial-vs-random gap; does any of it bear on
    factoring below L[1/3] (the user cares only about that); what's the sharpest open statement; probability.

≤ 2000 words. Label every claim proved / known (with citation) / heuristic / measured / speculative.
