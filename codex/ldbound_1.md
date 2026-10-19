# ldbound_1: blind round on Hypothesis LD (lift dispersion), structured residues, and averaging

You are a research analytic number theorist acting as an independent, adversarial collaborator.
The repository /Volumes/git/factoring-connection is read-only for you (scratch: your cwd).
Rigour matters far more than speed. Label every claim PROVED / KNOWN (with exact citation and the
exact hypotheses/ranges of the cited theorem) / CONJECTURED / FALSE. Do not hide losses: L = L_N[1/2,1],
and "L^{o(1)}" means exp(o(sqrt(log N log log N))); N^{o(1)}, L^{epsilon} for fixed epsilon, and
fixed-epsilon losses are NOT L^{o(1)}.

## Read first
- /Volumes/git/factoring-connection/notes/pcest_results.md, in full. Especially §1 (design (a''_k),
  Lemma GB, Lemma 1.2: g1^d F_A(P) = N b^d mod G(P); Lemma 1.3), §2 (Theorem FM, Lemma 8.5'),
  §3 (Lemma Q, Lemma DS, Theorem R), §4 (Lemma Kern, Prop 4.2, Hypothesis LD verbatim, Theorem 4.3,
  §4.3 barriers B1', B2', B3', B5, Prop Z, the "plausible route" with pieces (i) and (ii), and the
  single-family form LV*).
- /Volumes/git/factoring-connection/notes/hypD_results.md §4 (RBDH', Prop G).
- /Volumes/git/factoring-connection/codex/pcest_1.out, pcest_2.out, pcest_3.out.
- Lee–Venkatesan arXiv 1805.08873 §5 and §8 (Lemma 8.5 and the majorant W(Q)); Harper arXiv 1208.5992
  (Theorems 1, 2, Props 1, 2); Harper JLMS 2025 (arXiv 2412.19644); de la Bretèche–Drappeau JEMS 2020;
  Drappeau; Pascadi 2025 — as needed.

## The problem
Theorem 4.3 of pcest shows: FM + Theorem R + Prop 4.2 + LD imply PC <= L^{o(1)}, hence D_k, hence
Theorem A (rigorous L[1/2, (1+3/(d+1))/sqrt 2 + o(1)] for every N, 0.919 at d = 9) modulo LV §6.
What Theorem 4.3 actually consumes is the bound

  (LD*)  E_ell sum_{P != P', ell | Delta} |G|^{-1} sum_{g in G cap L_ell(P), filtered} Y_P Y_P' (S^2/q) kappa_alg
            <= L^{o(1)} S^4 rho_r^2 ,

with kappa_alg as in pcest §4.2 (lift classes F^ + G j mod G|Delta|, |Delta| = ell q, q <= 2S^2).

## Questions (answer each; give full proofs for anything labelled PROVED)
Q1. Piece (ii): large values for the structured moving residues r0 = N (b/g1)^d mod G. For fixed modulus
    t = G, the fractions b/g1 come from a box of area ~ t; their representation counts are large (~K) only on
    a set of measure ~1/K^2 (Farey-type fractions b0/g0 with small b0, g0). Is there a proof that these
    structured residues carry no excess lift-L^2 mass at precision L^{o(1)}? If not, give a precise barrier
    (a construction or a model showing which tools cannot work).
Q2. Averaging. LD is averaged over ell, q, the field (g0, g1, c), and pairs. Does this averaging already
    spread the structured residues enough to reduce (ii) to an all-residue statement? Be careful about the
    coupling between the residue of P (which depends on g1 and n = G/ell) and the rational weight Y_P' of the
    partner point (smoothness of n' = G'(P')/ell, an affine function of g along the fibre G(P) = ell n).
    Exhibit exactly where any loss arises.
Q3. Piece (i): formulate and prove (or refute) an all-residue L^2 restricted Barban–Davenport–Halberstam
    upper bound, at precision L^{o(1)}, for the single family of moduli ell^2 n q (ell prime in J, n B-smooth
    in a dyadic range, q <= 2S^2), with values up to X_F and moduli up to ~N^{s+3/(4k)}. Track: imprimitive
    characters / non-reduced classes, small conductors (Siegel-type issues), the large-sieve diagonal, the
    sparsity of the family, and which averages (over ell? over n?) are essential.
Q4. Design changes. Consider changing the design (more randomness in f's coefficients, a different
    rational map, extra random parameters, ...) so that LD's structured residues disappear or become
    irrelevant. Any change must keep FM, Theorem R and Lemma Kern valid (re-check them) and keep the
    constant (1+theta)/sqrt 2 < 1; compute theta(d, k) and the infimum constant for any design you propose.
Q5. If you obtain a proof of LD* (for (a''_k) or for a modified design), write it out completely, with
    every cited range checked, and then attack it yourself as a referee would.

Finish with a short verdict: which of LD / LD* / (i) / (ii) you PROVED, for which design, at what cost in
the constant, and the smallest precise statement still missing (if any).
