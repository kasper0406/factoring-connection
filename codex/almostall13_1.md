# almostall13_1: blind round. Below L[2/5] for ALMOST ALL N by averaging over the input?

You are a research referee and co-investigator. The repository /Volumes/git/factoring-connection is a research project on
rigorous subexponential integer factorization. It is read-only for you; write any scratch files only under
`codex/scratch/almostall13_1/`. Work blind: do NOT open any file named `notes/almostall13_*` or `codex/almostall13_*`
other than this prompt.

## Read first
- `notes/lv5unif_results.md` (Theorem U: claimed rigorous Las Vegas L_N[2/5, c] for EVERY N, conditional on the
  project's internal unrefereed chain; design (a‴): F_A = F̂_{g,w} + G·R_c with G = g₁x − g₀y, F_A(g₀, g₁) = N·w,
  c ∈ [−H, H]^d random, w a random prime) and `notes/auditU_audit.md`.
- `notes/nested_results.md` IN FULL (why multi-level descent below 2/5 is blocked worst-case): Prop SI (sparse image),
  Lemma FB, Lemma RW (algebraic components live in one field), Lemma OD (the per-field parameter budget d + 1),
  constrained resampling, Prop R / R2 (premises M1–M4), Prop BT, Lemma X_ERH and the L[1/3, Θ(δ)] character-prime cost,
  Theorem NR (ERH + Hypothesis AD + LV6′ ⇒ L[1/3]).
- `notes/recurse_results.md` §3 (size model, recursion Lemma F_α, Lemma S″).
- `notes/polyimpl_results.md` §1 and §7 D1 (transfer lemma: exceptional density o(n^{−2}) among n-bit integers ⇒ o(1)
  among balanced semiprimes).

## The idea to evaluate
Every attempt below 2/5 so far was worst-case over N. The blocker (nested) is that algebraic children have too few free
parameters within one field (Lemma OD, budget d + 1). Proposal: AVERAGE OVER N TOO. The base-m polynomial's coefficients
are essentially the base-m digits of N (or of N·w); for N uniform in [X, 2X] they are nearly uniform and independent. The
hope is that this adds ≈ d + 1 random parameters per field and couples the random field to the random input.

GOAL: a rigorous Las Vegas algorithm that factors all but an o(1) fraction of N in [X, 2X] (or of balanced semiprimes,
via the transfer lemma) in expected time L[1/3, c], or L[α] for some α < 2/5. Determine the best α achievable.

## Questions (answer each, with labels PROVED / KNOWN / MODEL / HEURISTIC / SPECULATIVE)
1. Re-examine nested's obstructions with N random. (a) Does Lemma OD's budget argument still bind when the field
   coefficients include N's digits? Count parameters exactly: how many genuinely new degrees of freedom does N add, in
   design (a‴) and in the plain base-m design without c? (b) Does Prop SI's sparsity disappear on average over N?
   (c) Which nested-descendability statement becomes an average over (N, field), and is it provable with Theorem U's
   τ = 0 machinery plus averaging over N (large sieve / Barban–Davenport–Halberstam-type averages with N as a
   modulus-like quantity)? Is there any place where N-averaging buys something genuine (e.g. Lemma W, rigour of the
   counting premise M1 of Prop R2)?
2. Multi-level descent for almost all N: size model (recurse §3.4), recursion (the halves u₀, v₀ and the values F(P),
   G(P) are factored recursively — what does an "almost all" guarantee for the algorithm imply for these structured
   sub-inputs, and what exceptional density is needed for the recursion to close?), transcript invariance (T), LV6
   at growing degree. Is ERH needed for the character-prime repair (nested: L[1/3, Θ(δ)] unconditionally), or does
   averaging over N help?
3. Quantifiers: state precisely what "almost all N" must mean (uniform fixed algorithm; the expected-time bound fails
   only on a small exceptional set; Las Vegas correctness always). State the transfer to balanced semiprimes with a
   quantified exceptional density.
4. Literature: are there almost-all-N (or average-over-input) rigorous subexponential factoring results below L[1/2]?
   Rigorous ECM-for-almost-all-primes statements? Analogous "average over input" rigorous NFS / discrete-log results
   (e.g. Lee–Venkatesan 2018, Kleinjung–Wesolowski, Pomerance 1987, anything else)? Cite precisely; say what you
   verified in a primary source.
5. Verdict: the best α you can support (rigorous, conditional with named hypotheses, or model only), with the exact
   missing lemma if it is not 2/5's machinery. If you find a genuine route below 2/5 for almost all N, write the proof
   in full.

Be adversarial with the idea itself: if N-averaging adds nothing, prove it as cleanly as possible; if it adds
something, quantify exactly what. Output a single report (markdown, ≤ 6000 words).
