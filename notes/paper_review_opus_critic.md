# Opus critic review of paper/main.tex (2026-10-18)

**Verdict: ACCEPT WITH FIXES.** No mathematical errors were found in the D′/D″ chain and nothing breaks compilation. The orchestrator transcribed this summary from the critic's report.

**Must fix**
1. **Quantifier overclaims** (r12 audit item 5). Theorem A bounds only full-candidate covering; D′ concerns the true factor, averaged over pairs.
   - (a) L798–799, the "Interpretation" remark: say "contain the true factor for almost every pair in S_n; whether they cover almost every candidate for fixed N is not proved; covering every candidate is Conj DCH".
   - (b) L182–183: say "floor for full-candidate covering in this model; it does not bound average-case algorithms".
   - (c) L972–973, "Below 1/6": add the same qualification.

**Should fix**
2. §4.1 model, L699: residue classes must be unit classes, r_τ ∈ (Z/M_τ)^×.
3. §4 per-target bounds:
   - (a) Add the e=1, K=1 case.
   - (b) Add window merging/rounding: either pay +C per original window with #windows ≤ ΣW ≤ X² = o(P), or add the cost hypothesis #windows ≤ X.
   - (c) L757: "If X = δP^{2/5}" → "If X ≤ δP^{2/5}".
4. Lemma lem:sieve, L536: define J_k := Σ_{ℓ∈P_k} |V|/(ℓ−|V|) ≥ c₁H²/n, and reconcile c′ with c₁.
5. Abstract L65–67: the slopes 0.203 vs 0.226 are over 48–80 bits; the 2.3–4.9× ratio is over 40–80 bits.
6. §6(a), L942: DCH's band and weights differ from DCL-det's. Write "for a variant of the algorithm with the band/weights of the conjecture, this would follow".
7. §5.1, L828–830: "This lies below the range of Cor polylog; at H ≍ √n the proved bound n²/H² exceeds 1 and is vacuous."
8. Carella paragraph (L301–302, L309–310): the descriptors "balanced N=pq", "arithmetic operations", "(lattice) mechanism", "worst-case" are not in our notes. Verify them against arXiv:1308.2891 or mark [sweep]/UNVERIFIED.
9. Acknowledgements L1002: write "shared base β^L (explicit example β=3, N=1333)", not 2^L.
10. Theorem B proof, L774–776: label the class-count bound a sketch (brute-force checked for composite M ≤ 105) and keep O(1)^{ω(M)} in the height count.

**Nits**
11. Wording:
   - L131: "Equivalently" → "In particular".
   - Abstract: write "primes ℓ ≡ 3 mod 4", and say the polylog exponent depends on A.
   - E_c, L618–620: generalize the exponent statement for ε ≤ 1/60.
   - Step 3, L365: notation.
   - R-d, L899: "within 0.005 of HH (pre-registered threshold)".
   - L1011: say "every proof step"; r13 marked the implementation WRONG and the §11 barrier a GAP.
   - L961: label "q/p ≤ N^{o(1)}" as heuristic.
12. Prior art:
   - L261–262: Hittmeir 2021 is a "time–space trade-off for Lehman's method"; mark [sweep] where unverified.
   - L219–220: either mark Lawrence UNVERIFIED in the text or drop the sentence promising markers.
   - L98: note the Carella claim.
13. Notation clashes:
   - P = √N (§4) vs P = 2^{n/2}.
   - d, k, H, e and c are each used for two things.
   - W_max is undefined when T_k = 0.
14. LaTeX: wrap the \ref in the titles at L410 and L673 in \texorpdfstring.
15. FYI, inconsistencies in the notes; the paper is correct in both cases:
   - dcltheorem: "896 runs" vs 656. The paper's 656 is right.
   - multimod: DCH candidate counts in §3.6 vs §6.2. The paper's §6.2 values are right.
