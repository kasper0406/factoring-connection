# detqs_2 — adversarial referee of notes/detqs_results.md

You are a hostile but fair referee (analytic number theory). Read ../../notes/detqs_results.md in full, plus
../../codex/detqs_1.out and ../../codex/scratch/detqs_1/fixed_u_moment_note.md (your own blind round), and
../../notes/detqs_prereg.md. Code is in ../../src/detqs/.

Your job: find every error, gap or overclaim. For each claimed result give a verdict
(CORRECT / CORRECT WITH FIXES / GAP / FALSE) with the exact fix.

Priorities:
1. **Proposition V / Theorem H1-S (§4).** This is the main claim: stratifying by N mod W (prime powers)
   and by archimedean slot, plus "averaging the main-term covariance over k", gives
   limsup (log B)² E(Z − M_σ)²/T² ≤ C r(W,J) with r → 0, hence H1 outside o(B/n²) integers.
   Attack especially: (a) step 4 — is Harper Theorem 1/2 applicable as used (moduli q = 2kW distinct,
   scales J/g, non-reduced classes, the cross terms with bounded multiplicity ν(W), the G-truncation)?
   (b) step 5 — do the marginals of the coupled residue law μ really equal the product-term marginals
   (including the per-t first moments, which are long-interval averages with t-dependent weights)?
   (c) step 6 — are (LD-a)–(LD-c) true at the stated precision from Ψ(x,y)=xρ(u)(1+O(1/log y)) plus
   inclusion–exclusion (note q′ can have many prime factors; Π(q′), σ(q′))? Is LD-c really "first order
   times σ(q′)", or does the O(1/log y) error of each Ψ break it?
   (d) step 7 — the claimed exact cancellation of the W-prime covariance after averaging k over residues
   mod W; the asymptotic independence of λ_k (depending on primes ∤ W of k) from k mod W; the
   truncation of D_W at valuation a; the archimedean slot giving J^{−2}.
   (e) is anything missing — e.g. the 2-adic/parity conditions w − v ≡ k² (mod 2k), the interaction
   between the W-part of k and the precision e_p = v_p(2kW), the dependence of boxes on k.
   If V is false as stated, give the correct statement; if it has a gap, say whether the gap is routine.
2. **Prop B** (L[1/2] barrier for residue-stratified Chebyshev): is the (log n)²/n² lower bound right?
3. **Lemma P, Prop S, Prop SZ, Theorem DQ** (cost N^{ωc+o(1)}, deterministic LSP elimination), the
   transfer, the "not average-time" remark.
4. **Experiments §8**: are the frozen verdicts and post-hoc labels honest? Any overclaim?
5. **Literature**: anything known that already proves H1 for almost all N for the QS polynomial, or a
   deterministic almost-all factoring result below N^{1/5}?

Output: numbered findings, each with verdict and fix; then a short overall verdict on what label each
result deserves. ≤ 2500 words.
