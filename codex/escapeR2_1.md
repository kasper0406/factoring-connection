You are a research mathematician (analytic + algebraic number theory, rigorous index calculus). Work READ-ONLY on /Volumes/git/factoring-connection (your cwd is a writable scratch dir). Use web search freely. Rigour matters far more than speed. This is a BLIND round: produce your own analysis; do not look for or rely on any file named notes/escapeR2_*.

READ FIRST (in this order):
- README.md;
- notes/lv5unif_results.md (Theorem U: rigorous Las Vegas L[2/5] factoring via Dixon split + recursion + ONE-level special-q descent, no algebraic children; §1, §8, §10, §11 (shared-rational-side multi-field trick, Prop MF), §13 (why 2/5 is the floor of that design));
- notes/nested_results.md IN FULL (Prop SI: sparse image for a fixed field; Lemma FB; Lemma RW (row algebra); constrained resampling; Lemma OD (parameter budget d+1); Prop R (rational children exogenous via fresh fields); Prop R2 (designs without algebraic children have model floor 2/5); Prop BT; Theorem NR (ERH + Hyp AD ⇒ L[1/3]));
- notes/recurse_results.md §3 (size model; Lemma S″; multi-level descent reaches 1/3 in the model);
- skim notes/auditU_audit.md and APPROACHES.md (do not re-propose what is listed there without a new twist).

PROBLEM. Get a rigorous (unconditional, or under ERH as fallback) Las Vegas factoring algorithm for all N with expected time L_N[α, c] for some α < 2/5, ideally α = 1/3. The project's obstruction: below 2/5 the descent must be multi-level, and at level j ≥ 2 the prime ideals to be descended ("algebraic children") are factors of special-q algebraic values IN THE SAME FIELD; their descendability cannot be refreshed by averaging over that field (Prop SI, Lemma OD), while rational children can be exported to fresh fields (Prop R) but designs without algebraic children stall at 2/5 (Prop R2).

Analyse, with honest calculations, each of these escape families (and any better one you find):
1. EXPORT ALGEBRAIC CHILDREN to a fresh field:
   (a) replace a degree-one prime 𝔭 = (ℓ, r) of K = Q(α) by its norm ℓ viewed as a RATIONAL target, (ℓ) = ∏𝔭_i^{e_i}; what happens to the other primes above ℓ (Galois conjugates, Galois closure), and to the square-class / character bookkeeping (LV6-type post-processing: Lemma H, characters, Lemma SQ square roots) when the row algebra must cancel 𝔭 but only ℓ is available?
   (b) a second field K′ = Q(α′) sharing the root m mod N (MNFS / Coppersmith style: all fields share the rational side g = g₁x − g₀y) in which ℓ is re-descended; can the row algebra be rearranged so that each algebraic child of K is closed in K′ or via the rational side? (Note: an ideal 𝔭 of K is not an ideal of K′; ask what a "relation" that kills 𝔭 needs.)
   (c) the image under φ: Z[α] → Z/N of an S-unit or generator of 𝔭^h.
   For each: the exact square-class/relation bookkeeping, whether rigorous LV6-style post-processing survives with several fields, and whether the export is actually possible.
2. CLASS-GROUP DESCENT: express [𝔭] in Cl(K) via relations in the same field and descend only a principal part. Is rigorous class-group relation generation at L[1/3]-type sizes available (Biasse, Biasse–Fieker, de Boer–Pellet-Mary–Wesolowski arXiv 2512.01588 under ERH)? At which norm scale (dBPW sample at scale ≥ √|Δ_K|)? Does this reduce to the same Arakelov-scale L[1/2] barrier?
3. GALOIS / CHEBOTAREV: fields with large Galois group; conditioning on splitting types so that algebraic children of one field become "rational-like" across a family of conjugate fields.
4. CHANGE THE PROBLEM so no level-2 descent is needed:
   (a) several independent one-level descents whose products combine;
   (b) a different split law (Lemma S″ / Theorem S works for ANY target set) where the split halves are descended by recursion (factoring) rather than special-q;
   (c) factor-then-descend: recursion on intermediate targets of size N^{small}. Recursion makes factoring any public number cheap relative to T(N). Does "factor the special-q values recursively and RE-SPLIT them" (treat a smooth-but-large-prime special-q value as a new Dixon-type split target with fresh R's) bypass endogenous children? Note that a re-split needs the large prime to appear in a square relation mod N, not just to be factored.
5. For each approach: the exponent α it would give in the size model (recurse §3.3–3.4, nested Prop R2 premises M1–M4), and the EXACT missing rigorous lemma.

DELIVERABLE (Markdown): (i) a ranked table: approach | α in the model | missing lemma | known/open/false | your probability; (ii) for each approach a careful analysis with calculations; (iii) any approach that you believe gives a RIGOROUS α < 2/5 (even under ERH): state it as a theorem with a complete proof sketch, list every input with references, and then attack your own proof adversarially; (iv) any approach you can REFUTE (show it reduces to a known barrier, or is circular, or the model forbids it): give the refutation precisely. Label every claim PROVED / HEURISTIC / MODEL / KNOWN (ref) / SPECULATIVE. Do not inflate.
