# polyimpl_2 — adversarial referee: consequence catalog, escape-signature table, new theorems

You are a hostile but fair referee (number theory + complexity). Repo at /Volumes/git/factoring-connection is
read-only (scratch: cwd). Read notes/polyimpl_results.md in full; it builds on your own blind round
codex/polyimpl_1.out and on notes/below16_results.md §2.1 (Theorem T1, the model and budget conventions),
notes/multimod_results.md (Thm A), notes/coverage38_results.md, notes/dchworst_results.md, and Ward–Maestro
(https://igent.ai/research/floors-not-bounds/factoring_report.pdf; Def 32.1, Conj 33.1, Prop 32.2, §35).

For EACH numbered item give a verdict CORRECT / CORRECT WITH FIXES / WRONG / GAP / OVERCLAIM, with the exact fix.

1. Theorem DR (§1): IW97 ⇒ deterministic L_N[1/2,O(1)] factoring for all N, by derandomizing Lenstra–Pomerance;
   and the contrapositive (a deterministic lower bound ⇒ E ⊆ i.o.-SIZE(2^{εm}) for all ε). Check the PRG
   parameters (seed O(log S), computable in poly(S), fooling size-S circuits), the role of Las Vegas truncation,
   uniformity, whether the hypothesis must hold for ALL large input lengths, and the quantifiers of the
   contrapositive. Is this stated anywhere in the literature (give a citation if so)?
2. Theorem CB (§4.3): T1-local collision search + (Z/N)^* B2-smooth-exponent annihilation (B2 ≤ (P/4)^{2/3},
   fixed base set A with Σ log|a| ≤ P^{1/3}/n^2, adaptive, no cost bound on Y1) fails on a positive proportion of
   S_n unless X ≥ cN^{1/6}/n^{4/3}; conditional on hypothesis F_iv (Fouvry-type positive proportion of primes in
   [aX,bX] with P+(p−1) > p^{2/3}). Attack: the box choice, the small-order removal count, ℓ_p ∥ p−1,
   ℓ_p | ord_p(a), the adaptivity argument in step 4 (is "every Y1 answer is 1 on G ⇒ runs as a fixed T1 design"
   valid given T1's locality/budget conventions?), whether a^E values used as T1 handles break T1's hypotheses,
   whether T1's bound applies to the sub-box G, and whether F_iv is actually known (Fouvry 1985, Baker–Harman,
   later work) in interval form — give the precise reference if so. Is the theorem new or trivially implied by
   WM §35 + T1?
3. Lemma KE and Theorem CS (your own §8, transcribed): recheck once more, especially Coll ≥ 1/R and the
   adaptivity of Model X after Model Y.
4. Lemma PF + Corollary PF (§3): the strip route e_strip = max(β,1−θ)/θ (θ > 1/2) and the claim that for all
   known short-interval theorems factoring is irrelevant; the region 2β < θ ≤ 2/3 where factoring could matter;
   the H4 condition β < 2/15 + δ. Any known theorem with θ > 2β that we missed (e.g. results for intervals of
   length x^β with β < 1/2, multiplicative/polynomial families)?
5. Prop LW: WM Conj 33.1 ⇒ least exposing Suyama index O(log p) for all but one prime per range; the comparison
   with least quadratic non-residues (Montgomery's GRH lower bound) — is the comparison fair?
6. The escape-signature table (§5) and the signature (S1)–(S4): any row with the wrong model, wrong
   "must break", or a wrong/missing witness? Is (S4) (super-constant depth; not TC^0 under Kalai) correctly
   stated? Any known algorithm that already has the full signature (classical, deterministic, rigorous,
   almost-all or worst case, below N^{1/6})?
7. Catalog §2: flag every wrong or overclaimed entry (H1–H5), incl. Kalai TC^0 implication, Sarnak quote use,
   Żrałek, Dieulefait–Urroz, Shallit–Shamir, Buchmann–Lenstra, class groups "not implied", Bach–Charles.
8. Experiment E1 (§8): is the generous exposure criterion P+(Φ_d(p)) ≤ B2 the right one for Theorem CB⁺?
   Is reporting "mixed" and the post-hoc reading honest? Is D6 "refuted for cyclotomic portfolios" justified?
9. Directions §7: for D1 (fixed-u deterministic QS almost-all, transfer from all integers to S_n), check the
   transfer lemma and the claimed lattice/AP reduction of the second moment; is the kill criterion right?
   Rank the directions yourself.
Output one markdown report, items 1–9, then a short list of required fixes.
