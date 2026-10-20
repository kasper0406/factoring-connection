# resbarrier_2 — adversarial referee of a generalized resultant barrier

You are a hostile but fair referee (number theory / computational number theory). Read-only repo at
/Volumes/git/factoring-connection; scratch in your cwd; web search enabled.

Read notes/resbarrier_results.md (the claims), src/resbarrier/{forms,ledger,e1}.py, tests/test_resbarrier.py,
notes/resbarrier_prereg.md, and for context codex/resbarrier_1.out (your own earlier blind round) and
notes/heurnfs_results.md §§1, 4.

Attack every claim labelled proved or proved (model). For each give: CORRECT / GAP (fixable; say how) / WRONG
(counterexample). In particular:
1. Lemma ND (§2.2): is "gcd([O_A:R], N) = 1 ⇒ ring-map extraction is useless" correct as stated? Are there
   extraction schemes in the stated model that escape it (e.g. maps from overorders that are not localisations,
   non-maximal overorders, elements whose square roots lie outside R[1/c])? Is the step "ψ1 = ψ2 on R[1/c] =
   O_A[1/c]" right?
2. Lemma C2 (§2.3): check the row/column interval argument, the use of Bean 1994 Thm 2 (A_F·|D_F|^{1/(n(n−1))} ≤
   3B(1/3,1/3), n ≥ 3), the disc identity disc(F1F2) = disc F1 disc F2 Res², the sublattice computation (c), the
   nonlinear-parametrisation claim (d) (N^e | Res for a shared degree-e factor mod p), and the derived ledger in (e)
   and Corollary RB2 constants (16k/9)^{1/3}. Is "SNFS is optimal among single-pair NFS for every N (H model)"
   justified? Is the k = 1 claim for ramified/class-group architectures justified (what exactly forces N | disc
   of the order the candidates live in)?
3. Theorem GAE (§2.4): check the proof and its restrictions; is it really unknown/folklore?
4. Lemma KL (§2.5): check, including the linking hypothesis and the claimed ledger consequence.
5. §2.6: the claim that rigorous constraints permit the aggregate ledger, the 1/4 floor, the G_r heuristic with
   Schmidt's counts (Schmidt 1993: X^{d(max(d,n)+1)} ≪ #{deg-d points of P^n, H ≤ X} ≪ X^{d(d+n)}), the
   "r − 1 independent linking conditions", and the expected-count exponent N^{(d+1)/d − (r−1)}.
6. §4.1 Kirchner reading and the "nonlinear congruence disc(f) ≡ 0 mod N" construction obstacle; §4.2.
7. The escape table §3, especially (d) non-commutative (strong approximation / Eichler class number), (h) Schnorr.
8. The MNFS RES-only floor 1.7114 (heurnfs mnfs(kr=1, ka=0.5)): is the modelling (kr → 1, ka ≥ 1/2) right?
Finally: is any sub-L[1/3] mechanism hiding in a gap you found? List required fixes in priority order. ≤ 3000 words.
