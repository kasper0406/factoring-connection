You are a research mathematician (analytic + computational number theory), working blind and independently. The
repository /Volumes/git/factoring-connection is read-only for you; write any scratch files under your current directory.

## Context

The project claims (internally, unrefereed) **Theorem U**: a Las Vegas algorithm that completely factors every integer N
in expected time L_N[2/5, c]. The best published rigorous bound for arbitrary N is Lenstra–Pomerance 1992, L[1/2, 1+o(1)].
Theorem U is in `notes/lv5unif_results.md` (§§1–10), resting on the fixed-degree preprint `paper2/main.tex` ([pp]),
`notes/lem85_results.md` (Lemma 8.5′, restricted BV for smooth numbers over moduli ℓ·(smooth)) and
`notes/recurse_results.md` (Lemma S″, Lemma F, Lemma OS, recursive shape test). A blind audit is in
`notes/auditU_audit.md`.

`notes/lv5unif_results.md` §11 proposes **Design MF**: one shared random rational side g = (g₀, g₁), and J_f independent
algebraic sides (w_j, c_j), j ≤ J_f, F_j := F̂_{g,w_j} + G·R_{c_j}, all with the same root m = g₀/g₁ mod N. A special-q
prime ℓ is MF-descendable if some candidate P in the ℓ-lattice has smooth rational cofactor G(P)/ℓ and F_j(P) is B-smooth
for SOME j. The claimed advantage: conditional on (g, ℓ), the fields are independent, so a pointwise one-candidate
first moment (Lemma GM-U) plus independence replaces the algebraic second moment (LD*, Kern, L2_fam, Lemma W). §11
states "Proposition MF [moment part proved; multi-field post-processing proved in outline]" with an extra height floor
log H ≥ 13(σ+1)dΛ log log N.

Blindness: do NOT read notes/mfU_*, codex/mfU_* (other than this prompt) or STATUS.md. Everything else in the repo and
the web is allowed. Treat every internal claim as unverified.

## Tasks

1. **Write the multi-field post-processing completely** (the step §11 only outlines). Concretely: given rows produced
   by Dixon descent where each special-q atom may come from a different field j, define the matrix (columns: shared
   rational primes/sign; per field j: parity, prime-ideal (r, s) columns, quadratic characters), the row target, and
   prove: (a) after character filtering, the remaining square obstructions define at most J_f F₂-linear functionals
   on the relation kernel (or say exactly what the correct count is and why); (b) with enough kernel dimension one can
   annihilate them and certify, with exact verification, that each field's product is a square; (c) the square roots
   in all J_f fields and the final congruence X² ≡ Y² (mod N), with one shared map α_j ↦ m; (d) transcript invariance
   for Dixon descent ([pp] Theorem `thm:DD`, Remark `rem:T`): Y must be a function of the squares R_i and of coins
   independent of the hidden roots; (e) explicit costs: rows, columns, characters (the character primes at growing
   degree have (2d)!·poly(log N) bits, lv5unif §9), square roots, linear algebra, as L[2/5]-exponents. Is
   J_f = (log N)^C/ρ_F fields affordable? Does the row target grow with J_f? Put the extra MF height floor and every
   per-field cost into an explicit cost ledger and give the resulting constant c as a formula.
2. **Hunt for failure modes of Design MF.** In particular: can rejected fields (irreducibility certificate, content,
   f′(m), anything else in SAMPLE_FIELD) absorb the pointwise smoothness mass ρ_F/(log N)^C? (Rejection probabilities
   like 2^{−d} or exp(−log H/(13 d log log H)) are far larger than ρ_F = exp(−Θ(D^{2/5} log log N)).) Is the
   "first good P₀" conditioning legitimate given that the success event involves fields' acceptance? Does simultaneous
   certification over J_f fields add an unaccounted cost? Is the rational supply γ_r ≥ (log N)^{−C} really proved
   from [pp] Proposition `prop:R` ("Theorem R") and the rational part of the first moment?
3. **Unrefereed inputs.** Drappeau, Compositio 151 (2015), Lemme 5 is stated without proof as a consequence of §3.3 of
   Harper's unpublished preprint arXiv:1208.5992. Where exactly does it enter Theorem U / Design MF (it enters via
   Lemma 8.5′, which makes most moduli t = |G(P)| "good")? Find, if possible, a route to a rigorous L[2/5] that avoids
   it: e.g. a different field law that makes the needed algebraic smoothness lower bound hold without any equidistribution
   of smooth numbers in progressions to a fixed large modulus; or a published BV-type theorem for smooth numbers with a
   relative saving in this range (y = L_x[1/2]) — Fouvry–Tenenbaum (1991/1996), Granville (1993), Harper (JNT 2012;
   Compositio 2016), de la Bretèche–Drappeau (JEMS 2020), others. For each candidate, state precisely what it gives and
   whether that suffices; if only a weaker input is available, which exponent α does the method then give?
4. **Assessment.** Is a self-contained, checkable proof of "rigorous L[2/5]" via Design MF realistic? What is the single
   weakest link?

## Output

Markdown. Start with a 5-line summary. Then full statements and proofs for task 1 (as a referee-ready lemma sequence),
then a table of failure modes (FATAL / MAJOR / MINOR, with fix), then task 3 with exact citations (author, venue,
year, theorem number) and whether you verified each, then task 4. Be concrete; do not pad.
