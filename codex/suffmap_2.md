# suffmap round 2 (ADVERSARIAL REFEREE): the map of sufficient conditions for faster factoring

You are an adversarial referee (computational number theory, smooth numbers, lattice algorithms, algebraic complexity).
Read the repository file `notes/suffmap_results.md` (agent suffmap) in full, and the code `src/suffmap/principle.py`,
`src/suffmap/listing.py`, `src/suffmap/test_suffmap.py` (you may run `.venv/bin/python -m pytest -q
src/suffmap/test_suffmap.py` from the repository root). You may consult the project notes it cites
(`notes/heurnfs_results.md`, `notes/smoothenum_results.md`, `notes/pilist_results.md`, `notes/resbarrier_results.md`,
`notes/slpcirc_results.md`, `notes/nonexact_results.md`) and the primary literature via web search. Work in the scratch
directory for computations. Your job is to BREAK the note: find wrong statements, mislabelled statuses
(PROVED / model-PROVED / KNOWN / HEURISTIC), quantifier errors, hidden costs, missing hypotheses, wrong citations, and any
sufficient condition that is not actually sufficient, or any classification (sufficient / necessary / unrelated) that
is wrong. Also say what is missing from the map (important sufficient conditions or high-profile problems not covered).

Focus in particular on:
1. Theorem SS (a)–(e): lower bound with assumption R; ω-independence claim; the power-law formula; the poly-time bound.
2. Theorem LS, Lemma SL, Prop PC (promise constructors: adversary bound 4H/h₀+1, tiling, and the "θ-equivalence"
   reading), Prop SI (single intervals near kN giving L[1/3, (16/9)^{1/3}]; the existence threshold; the claim that a
   promise constructor suffices given all-interval existence; part (iv) Bohr-set poly time; the heuristics on the
   residue and on oracle outputs). Is anything in Prop SI double-counted or too optimistic (e.g. are the L[1/3, κ]
   relations from distinct k independent; does "one call per interval" hide a cost)? Was the 1.2114 reduction published
   before (Boneh 2000, Schnorr, Pomerance, Lenstra, anyone)?
3. Prop DLT (modulus-generic listers ⇒ prime-field DL at the same exponent): correct? Does Prop SI then give an
   L_p[1/3, 1.2114] DL algorithm conditionally, and is the descent/individual-log step covered?
4. E4 (Pochhammer SLPs ⇒ P/poly) and the claim that the τ-conjecture contradicts the hypothesis; Bürgisser's exact
   theorem.
5. The tier table, the high-profile verdict table H1–H15, the "what is necessary" paragraph, and the ranking in §4.
   Do you agree that Prop SI is the most worth attacking and that the factorial/Pochhammer complexity is the unique
   "landmark either way"?

Output: numbered findings with severity FATAL / MAJOR / MINOR / NIT, exact location, the problem, a concrete fix.
Then an overall verdict (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION). ≤ 4000 words.
