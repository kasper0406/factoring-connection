# lowerbound: a unified relation-model theorem "α ≥ 1/3" (common preamble)

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch only in your cwd).

## Goal
No theorem exists saying "every relation-collection factoring method costs ≥ L_N[1/3]" in a BROAD model. The repo's
barrier map (BARRIERS.md, notes/completeness_audit_results.md §0/§5/§7) is a union of single-family model theorems in
disjoint models, with no compositional / amortized / adaptive theorem. Literature: Pomerance ICM 1994 and
Croot–Granville–Pemantle–Tetali (Annals 2012), Balister–Bollobás–Morris (Annals 2018) = random-squares thresholds
(L[1/2] level); Pomerance 2008 "size principle"; a 2025 Designs Codes Cryptogr. paper "On the complexity formulae of
the number field sieve and its variants" (doi 10.1007/s10623-025-01788-5) claims NFS and all variants are ≥ L(1/3)
(within NFS's own parametrization). An unconditional general lower bound is out of reach (it would separate P from
NP∩coNP), so the target is a MODEL theorem — as broad as possible, explicit about what it excludes.

Desiderata for the model M (state precisely):
- covers: QS/MPQS, NFS/GNFS, MNFS, SNFS/special forms, TNFS-like towers, lattice/line sieving, factorization factory /
  batch amortization over many N, large-prime variants, descent trees (heurnfs Assumption T), rank-r linked modules
  (rank3/rankdegen/grdeg; Theorem KH in notes/subthird_results.md), adaptive constructors that choose candidates
  using previous smoothness outcomes, compositions/hybrids of several families, Z-route (zroute) / square-root
  (F₂) dependencies;
- information primitive: what the algorithm learns about smoothness must come through (a) smoothness tests on
  candidates, priced, and (b) algebraic structure that is "N-public" — formalize this so that factoring-in-disguise
  primitives (oracles that already encode p) are excluded by a clean, checkable condition rather than ad hoc;
- arithmetic primitive: a size/height floor for candidates forced by N (resultant/discriminant/linking: resbarrier
  Lemma RES, Cor RB2, Theorem KH, suffmap Thm SS size–entropy law) — the theorem must say how the floor composes when
  candidates come from many families or adaptive choices;
- smoothness heuristic: the minimal, explicitly stated "randomness of smoothness" assumption (e.g. each tested
  candidate is B-smooth with probability ≤ ρ(u)^{1−o(1)} given everything the algorithm has computed so far, for
  u = log|candidate|/log B), and why it is plausible / where it is known to fail (special sets: exceptional small
  norms LP_exc^rel, beyond-Johnson decoders M2, smooth-number constructors Prop SI);
- conclusion: total cost ≥ L_N[1/3, c_M] for an explicit c_M > 0 (ideally matching MNFS's 1.902 or at least > 0).

Repo notes to use: notes/suffmap_results.md (Thm SS, size–entropy law, §1.2 model S and its lower bound), resbarrier,
rank3, rankdegen, grdeg, subthird (Thm KH, LP_exc^rel), heurnfs (Assumption T, Barriers), smoothenum (Prop LQ),
crtdecode, pilist (H_sel, Johnson), sihard (Thm J), smoothbuild, completeness_audit, equivtarget, multcopp.

Output (≤ 2200 words, plain language; label proved / heuristic / model-assumption): the model M (definitions),
the theorem statement, a proof (complete where you can; mark gaps precisely), the inclusion table (which known
algorithms/ideas are inside M — check NFS and MNFS reproduce α=1/3 as tight), the exclusion list (what M rules out by
assumption, and whether each exclusion is natural or a fudge), and the 3 weakest points of the argument.
