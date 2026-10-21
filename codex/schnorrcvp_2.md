# Adversarial review: factoring via approximate CVP in prime-number lattices

You are a hostile but fair referee (computational number theory + lattice complexity). The repository is at
/Volumes/git/factoring-connection (read-only for you). Read `notes/schnorrcvp_results.md` in full, then the code it
cites (`src/schnorrcvp/*.py`) where needed, and the raw-data summaries in `experiments/raw/schnorrcvp_E*.jsonl` if
useful. Web search is available: check the primary literature where a claim depends on it.

Your job: find every error, overclaim, unjustified step or mislabelled status. For each issue give severity
(FATAL / MAJOR / MINOR), the exact location (section/claim), why it is wrong, and a concrete fix. Specifically attack:

1. **Theorem P (promise reduction)** and its rounding paragraph: is the proof correct (sinh inequalities, A ≤ Q,
   the choice C = e^{(L+A0)/2}/ω, the scan over C)? Is the "every γ-CVP answer" claim right? Is S ≥ L^2 enough?
2. **Lemma LA**: is c ↦ x(c) really a homomorphism on the GF(2) kernel; is the odd-index claim right; is Bach's GRH
   bound (primes < 2 ln^2 N generate (Z/N)^*) correctly stated for composite N?
3. **§2.3 distinctness via prime multipliers k > B** and the claimed **flaw in Adleman 1995 §2** (replacing V0 by
   V0+Vi+Vj does not change the lattice). Read Adleman's text (SchnorrGate repo, LATFAC.pdf) and decide whether the
   flaw is real or whether Adleman's argument intends something else.
4. **Theorem H / G / Prop A**: the claim "factoring ≤ CVP_{1+n^{-ε}} on prime lattices", the cost table as a function
   of γ (check the L[α, 2√(ατ)] derivation in §2.6), the claim that γ ≥ 2 / poly(n) / 2^{n^ε} gives nothing, and the
   adversarial-tightness Prop A.
5. **Prop NP**: is "factoring ≤_T GapCVP_γ for γ ≤ n^{c/loglog n}" correct as stated (deterministic? Karp vs Cook,
   promise problems)? Is "factoring ≤_R GapSVP_γ for every constant γ" correct (Khot 2005, Haviv–Regev)? Are the
   coAM / coNP placements (Goldreich–Goldwasser; Aharonov–Regev) stated with the right γ?
6. **§3 counting**: the exponent a(1−2/κ) + b(1−1/κ) − 1/κ; the thresholds κ = 3/2 and κ = 2; the claim that naive
   H_S over-estimates primitive relations (is that justified beyond the toy data?); the claim n_min grows like
   (ln N)^{3/2+o(1)}; whether E_tot ≈ 1e-52 / 1e-115 at Schnorr's parameters can be trusted (model, Ennola regime,
   bracket method in psi.py).
7. **§4 Schnorr 2021**: are the Lemma 3.1 and 0.917 criticisms correct as stated? Is "no algorithm can find
   relations that do not exist" a fair summary given it is heuristic? Is the BKZ-β paragraph honest?
8. **E0b cost model**: is it optimistic or pessimistic for the method; is comparing to crude L[1/2,1], L[1/3,1.923]
   misleading; is "worse than trial division at 800 bits" justified?
9. Anything else: novelty claims ("no prior statement found by us"), missing literature (e.g. analyses of Schnorr's
   method with explicit relation counts, Micciancio's work on the Schnorr–Adleman lattice, Bright 2014 abc), errors
   in the experiments' interpretation.

Write a structured report (≤ 3500 words). End with a table: issue | severity | fix.
