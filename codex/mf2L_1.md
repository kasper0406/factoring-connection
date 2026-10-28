You are a research mathematician (analytic + computational number theory), working blind and independently. The
repository /Volumes/git/factoring-connection is read-only for you; write any scratch files under your current directory.

## Context

The project has two internal, unrefereed rigorous-factoring results:

- **Theorem U-MF** (`notes/mfU_theorem.md`, self-contained): Las Vegas factoring in expected L_N[2/5]. Its key device
  (Lemma OP, §6.2) draws the multiplier w uniformly from all *integers* in (W, 2W], so the residue of the algebraic
  value F(P₀) modulo the rational value G(P₀) is uniform and the smoothness probability is an *interval density*.
  Design MF (one shared rational side g, J_f independent algebraic sides) removes the algebraic second moment; the
  rational second moment (Prop R2) needs only an UPPER bound with polylog loss, done with the elementary Lemma DS.
  Hence no Bombieri–Vinogradov (BV)-type input at all; in particular not Drappeau, Compositio 151 (2015), Lemme 5,
  whose proof is a pointer to §3.3 of Harper's unpublished arXiv:1208.5992.
- **Theorem 2L-final(α)** (`notes/thmAB_results.md` §0, §12; composition in `notes/c12alg_results.md` §5 and
  `notes/c12rat_results.md` §9): Las Vegas factoring in L_N[α] for α ∈ (3/8, 2/5), and L[3/8] under a constraint, via a
  TWO-level descent (level-1 parents with B₂-smooth algebraic values; each large prime q ∈ (B, B₂] of a parent's
  algebraic value is a "child" that must itself be descended in a box of its lattice). The proof needs every child to
  fail with probability ≲ 1/(number of children) (multidescent Prop UB, `notes/multidescent_results.md` §3.2), i.e.
  *second-order (covariance) precision* D^{−α} for pair correlations: rational blocks (`notes/c12rat_results.md`,
  Theorem RAT, through its Theorem WA) and algebraic blocks (`notes/thmAB_results.md`, Theorem AB₂, §9 Lemma EN-A/FAM-A
  + Theorem WA). Both use Drappeau's Lemme 5 ("Fact D") to bound character errors coming from the lattice condition
  x₁ ≡ ιx₂ (mod Δ′) on the pair values. Other inputs: `notes/pc3_results.md` (Prop DC), `notes/uprime_results.md`
  (Hyp U′, Lemma 8.5H hybrid BV), `notes/c12alg_results.md` (design DW: w uniform integer; Cor TP-q pair law).

Blindness: do NOT read notes/mf2L_*, notes/audit2L_*, codex/mf2L_* (other than this prompt), codex/audit2L_* or
STATUS.md. Everything else in the repo and the web is allowed. Treat every internal claim as unverified.

## Tasks

1. **Dependency map.** List every place in the two-level chain where Drappeau's Lemme 5, Theorem WA, Lemma 8.5H, or any
   BV / Barban–Davenport–Halberstam / large-sieve input is used (one-point child mean, rational pair covariances,
   algebraic pair covariances, parent counts, U′). For each: can an integer-multiplier / interval-density argument
   (Lemma OP-type) replace it, and how? Be precise about which statements are first-moment lower bounds (OP-type
   plausible), which are second-moment upper bounds with polylog loss (Lemma DS-type plausible), and which need
   second-order precision.
2. **Rational side.** The pair (G(P₁), G(P₂)) over g lies on {n₁ ≡ λn₂ (mod Δ)}, Δ = det(P₁, P₂). Can extra randomness
   (random shifts/multipliers on the rational side, randomising g, a second multiplier, averaging over pairs) make the
   rational residues uniform modulo the relevant moduli, so that pair counts become interval densities — without
   breaking F(g₀, g₁) ≡ 0 (mod N) / the common root m? Note that P₁ ≡ λP₂ (mod Δ) implies H(P₁) ≡ λ^e H(P₂) (mod Δ)
   for EVERY integral binary form H of degree e. What does averaging over pairs with the same Δ give (distribution of λ)?
3. **Algebraic side.** Under DW the pair law is uniform on L′ = {x₁ ≡ ιx₂ (mod Δ′)}, ι = λ^d, inside parallelograms.
   Can the covariance bound be proved from interval densities + Hildebrand–Tenenbaum short-interval asymptotics only,
   without characters? If not, what is the minimal residue-class statement needed (e.g. a variance over residue
   classes mod Δ′ for smooth numbers in short intervals, averaged over the family of Δ′), and at which moduli range
   (note children have Δ′ ≤ 2S₂² = y^{O(1)} with y = B, while rational child moduli contain the child prime q ≤ B₂)?
4. **Design changes.** Can an MF-type design (shared g, many independent algebraic sides; or other fresh randomness)
   remove the second-order precision requirement of the two-level chain, e.g. by converting per-child losses into
   search work? Note: algebraic children of field j must be descended in field j. Analyse the per-field joint event
   "parent accepted and all m ≍ D^α children good".
5. **Published inputs.** Give the minimal BV-type statement that the two-level chain still needs, and check whether it
   follows from a PUBLISHED result with full proof: Harper, Compositio 152 (2016); Harper, JNT 132 (2012);
   Soundararajan (2008); Granville, Acta Math. 170 (1993) and its Part II; Fouvry–Tenenbaum, Proc. LMS 63 (1991) and 72
   (1996); Drappeau's published theorems (Compositio 2015, other papers) with full proofs; de la Bretèche–Drappeau,
   JEMS 22 (2020); Hildebrand; others. For each: exact statement, range, error term, and whether it suffices at the
   needed precision (relative errors well below D^{−2α}, after the family losses).
6. **Assessment.** Is a "Theorem 2L-MF": rigorous L[α] for some α < 2/5 with only published inputs, realistic? If
   yes, give the argument; if no, isolate exactly where a BV-type input is unavoidable for this proof architecture.

## Output

Markdown. Start with a 5-line summary. Then tasks 1–6 with statements and proofs (or precise obstructions), exact
citations (author, venue, year, theorem number) and whether you verified each against the primary source. Be
concrete; do not pad.
