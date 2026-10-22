You are an adversarial referee (number theory / theoretical cryptography). Repo /Volumes/git/factoring-connection is
read-only for you; scratch goes in your current directory. Use web search to check literature. Be harsh and specific.

Referee notes/selfred_results.md (phase 22, agent selfred). Also read src/selfred/rigidity.py, src/selfred/reductions.py,
src/selfred/census.py, tests/test_selfred.py, notes/selfred_prereg.md and experiments/raw/selfred_e{1,2}.jsonl,
selfred_e1_posthoc.jsonl. Background: notes/factorbits_results.md §4.2 (Lemma PSR) and codex/selfred_1_report.md
(your own earlier blind report, from which Theorem H, Theorem L, the Bourgain corollary, SR2's rational exception and the
Fermat leak are taken).

For each item give: CORRECT / WRONG (with counterexample or the broken step) / NEEDS FIX (exact fix).
1. Theorem SR1(a): the Brun–Titchmarsh support bound, every case (ℓ ≥ y, √y ≤ ℓ < y, ℓ < √y), the hypothesis
   v_ℓ(G) > v_ℓ(c), the constant 4D/δ, and whether Montgomery–Vaughan's BT applies to a short interval (x0, x0+y]
   with a prime modulus ℓ < √y as used. Corollary SR1′'s averaging. Lemma DEG (finite differences). The "exact
   p-independence" refinement (is the Siegel–Walfisz step right? is the claim for ℓ | c with ℓ−1 > D right?).
2. Theorem SR2 and its exceptions; correctness and complexity of `integer_roots` / `sign_change_brackets` in
   rigidity.py (try to break it: repeated roots, roots at interval ends, huge coefficients, linear case); whether the
   operational code is leakage-free.
3. The RSA-vs-factoring explanation; SR3; the RSR literature paragraph (AFK §6 Question 1, FF Thm 3.1/Cor 3.3,
   Bogdanov–Trevisan, Beaver–Feigenbaum, Rabin per-modulus).
4. Status table §2.1: every row and reference, especially the claims about Shallit–Shamir, Woll, Landau,
   Buchmann–Lenstra Thms 1.3/6.12, Allender–Saks–Shparlinski Thm 1/Cor 2, Booker–Hiary–Keating, Adleman FOCS 1994
   C7–C9, the φ(N) mod ℓ (ℓ input) reduction (call count ≈ 0.35 n/ln n?), Lemma SIG, Lemma EUC.
5. Theorem L (GRH; m ≤ N^{1/2−η}): the Parseval/TV step, conditioning on N ≡ 1 mod 4, the smaller-factor ordering,
   and the "no secrecy threshold at N^{1/4}" claim. The Bourgain corollary (does Bourgain's theorem fix the top bit and
   allow r ≈ j prescribed low bits in [X, 2X] as used?). The KLM row; the Müllner row (is it honestly labelled?).
6. E1/E2: is the inequality A* ≤ 3TV₂/(1−2TV₂) correct as proved (derive it)? Are P1/P2 verdicts, the deviations
   D1–D4 and the post-hoc analysis reported honestly? Any bug in census.py (e.g. ordered vs unordered pairs, TV₂
   computation, U(h) computation)?
7. Overclaims anywhere in §0 (verdict) and §6. Anything missing that a referee would demand.
Deliver: a verdict (ACCEPT / MINOR / MAJOR REVISION) and a numbered list of required fixes.
