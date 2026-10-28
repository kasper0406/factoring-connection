# fruit_1: blind round — fruitfulness of a rigorous randomized NFS (Lee–Venkatesan Conj 7.1)

You are a research mathematician (analytic/algorithmic number theory). Work blind: do not look for an answer in
this repo's notes beyond the files named below. The repo is read-only for you; write scratch files only in your
current directory (codex/scratch/fruit_1/). Rigour and honesty matter more than a positive result. Label every
claim: proved / heuristic / conjectured / refuted-in-model / UNVERIFIED (for literature you did not check).

## Background (read these files; paths relative to the repo root /Volumes/git/factoring-connection)
- notes/nfsfruit_results.md §§1, 3, 4 (Lee–Venkatesan reading; Lemma H: fruitfulness is a homomorphism
  Φ: K0 → μ₂(Z/N)/{±1} for every odd composite N; Prop G: Φ(S) = Σ_{i∈S} β_i for an explicit hidden per-pair
  function β; Prop R and the hierarchy F/E/M/I; the involution criterion).
- notes/lemmaM_results.md §§1, 4, 5 (Theorem DD: Dixon descent — a uniform square R = r² mod N, split R ≡ u/v by
  Gauss reduction, each prime of uv eliminated by a special-q relation; any dependency containing such a row is
  fruitful with probability exactly 1 − 2^{1−ω} given transcript invariance; averaging over N refuted as a route;
  unit twist; square-twist identity).
- paper2/main.tex §§1.4 ("The idea in one page"), 4 (Theorem DD), 2 (Algorithm SPLIT): the project's rigorous
  L[1/2, c<1] algorithm, all of whose rows are descended.
- notes/lv5unif_results.md Verdict (one-level descent at growing degree ⇒ L[2/5], conditional on an internal chain);
  notes/heurnfs_results.md §3.1–3.2 (descent cost models; HECM L[4/9] barrier; multi-level descent).

Lee–Venkatesan (JNT 187 (2018), arXiv 1805.08873) prove a randomized NFS finds x² ≡ y² (mod N) in rigorous
L[1/3, (64/9)^{1/3}]; non-triviality is proved only under their Conj 7.1 and for Blum N.

## Questions (answer each with a proof, a model barrier, or the exact failure point)
1. **Hybrid / partial descent.** Generate the bulk of relations by LV's L[1/3] sieve (no hidden randomness), and
   add h Dixon-descended rows (hidden roots). Over GF(2): exactly how many descended rows are needed so that the
   final dependency is fruitful with probability ≥ 1/2? Must the dependency use them, and can the algorithm force
   that without breaking transcript invariance? What role does the rank deficit of the LV matrix play? What does a
   single descended row cost, rigorously and heuristically, *into LV's own field and factor base* (degree ≍
   (log N/log log N)^{1/3}, factor base L[1/3])? Give the exponent of a k-level descent into such a field.
   Does the hybrid ever beat pure descent in exponent?
2. **Randomising the square-root step.** Can randomness in the algebraic square root (embedding choice, unit twist,
   Couveignes / Montgomery with random sign choices, choice of inert primes) make LV congruences non-trivial
   without descent? Prove or refute.
3. **Involution criterion.** Find an involution on the smooth-pair set X that preserves the square-class vector x
   and flips β on a set of measure L[1/3, −o(1)]: candidates N·w multipliers, Galois automorphisms, (a,b) ↦ (−a,b)
   type maps, integer multipliers of the resultant. Or prove a structural obstruction (what does x∘ι = x force?).
4. **Equidistribution.** Can the joint law of (x mod 2, β) on X be shown close to a product (in which sense?) at
   L[1/3] parameters by decoupling / Bombieri–Vinogradov-type tools? If not, what exactly fails?
5. **Barrier.** If nothing works, state the sharpest model statement explaining why, and the minimal hypothesis
   that would suffice for rigorous L[1/3] (with which constant).

Also: any literature 2018–2026 on rigorous fruitfulness of NFS/QS congruences, or rigorous individual-logarithm
descent at L[1/3] (cite precisely; mark unchecked items UNVERIFIED).

Output: a structured report, ≤ 4000 words, with numbered claims and labels. Small exact computations (Python in
your scratch dir) are welcome; report what you ran.
