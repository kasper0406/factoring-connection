# inp round 2: orchestrator's pushback (Claude), then revise

Q1. Location 1 ("mixed modular/integer computation", 32%) is nearly vacuous as stated: reading integer
    representatives is computationally complete (leakmap says so), so EVERY algorithm is "mixed". What structural
    property would distinguish a plausible poly-time route inside it? Also, your probe (high bit of (2^N+3^N+r)^N mod N vs
    ⟨r, p⟩) is within what the entropy experiment already searched (depth-2 PRIM gates, Jacobi/bit readouts, GL parity
    targets, 20k tests, 600M semiprimes; only N^{−1/4}-decaying CF/genus leaks found — notes/entropy_results.md). Sharpen
    or downweight.
Q2. Location 2 (short matrix word that is the identity mod exactly one prime with inverse-poly probability, 18%):
    word-map theory (Larsen–Shalev, Larsen–Shalev–Tiep and successors on fibres of word maps in finite (quasi)simple
    groups) gives Pr_{g}[w(g) = 1] ≤ |G|^{−c} for nontrivial words in SL₂(F_p) / PSL₂, i.e. exponentially small, for
    random inputs. Does that kill the route, or does it survive with structured (non-random) inputs, fixed small
    integer matrices, or words of length poly(n) depending on N? Be precise.
Q3. Location 3 probe: why should [N+1−A](0,0) on 37a1 mod N with A = a_f(N) give a proper gcd? E(Z/N) ≅ E(F_p)×E(F_q)
    has order (p+1−a_p)(q+1−a_q), and a_N = a_p·a_q; N+1−a_N is not a group order of anything here. Give the mechanism
    or retract the probe (and say what recovery step from (N, a_f(N)) is actually plausible; cf. Bach–Charles using
    a_f(N) and a_f(N²)).
Q4. Missing locations? E.g. exact poly-time class numbers h(−4N) / regulators / L(1, χ_{−4N}) (genusleak: a multiple of
    hR finds an ambiguous form; nativeana precision wall); exact point counts of a fixed curve over Z/N (ecount: one
    count factors N); non-abelian/higher-genus Jacobians; "dequantized" period finding (Theorem T). Where do they rank?
Q5. Revised ranking with weights (conditional on factoring ∈ P), the single route you now think most likely, and the
    single best cheap probe the repo has NOT effectively run already. ≤ 1400 words.
