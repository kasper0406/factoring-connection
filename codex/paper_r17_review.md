You are an adversarial journal referee (target: Mathematics of Computation / ANTS) for the draft
paper /Volumes/git/factoring-connection/paper/main.tex (bibliography paper/refs.bib; compiled PDF
paper/main.pdf). The repository is read-only for you; write scratch files only under your
current working directory (codex/scratch/paper_r17_review/).

The paper claims, for balanced semiprimes S_n = {pq : p<q≤2p primes, n bits}:
- Theorem 1.1 (AV): a complete deterministic algorithm LAD, worst case = Harvey–Hittmeir's bound,
  average cost over uniform S_n of N^{1/6}(log N)^{c+1} log log N.
- Theorem 1.2: DCL_{n,H} (never wrong) with cost N^{1/6+5ε/12+o(1)} and O(N^{−ε}) exceptions for
  ε ≤ 1/40; and N^{1/6}polylog with (log N)^{−max(4a,1+2a)} exceptions.
- Theorem 1.3/7.1 (R): a modified algorithm DCL-wc factors every covered pair; Corollary 1.4:
  DCH_pair ⇒ worst case N^{1/6+o(1)}; coverage with c > 3/8 ⇒ worst case < 1/5.
- Theorem 1.5/8.2 (T1): a true-factor floor C n² X^{3/2} N^{−1/4} for residue-local degree-1
  sumset designs; Proposition 8.3: DCL lies in the model.
- Section 9 is an explicitly marked PLACEHOLDER (unbalanced ratios); do not referee its content
  beyond checking it is labelled as a cost analysis, not a theorem.

Authoritative project notes (read them to check faithfulness, but referee the PAPER):
notes/avgtime_results.md, notes/dcltheorem_results.md, notes/sieve2d_results.md,
notes/dchworst_results.md, notes/below16_results.md, notes/multimod_results.md.
Prior art read by the authors: Ward–Maestro, "Floors for deterministic integer factoring…"
(2 Sept 2026, https://igent.ai/research/floors-not-bounds/factoring_report.pdf); Carella
arXiv:1308.2891v3 §5; McKee–Pinch 1996 §5 (scans in references/mckee_pinch_1996/).

Please produce a full adversarial referee report:
1. For every theorem, lemma and proposition: CORRECT / GAP / WRONG, with the exact line or
   equation, and a concrete fix. Check every exponent computation (e.g. E_c exponents in the
   proof of Thm 1.2(a) and in §6 (d); the ladder reach-set inclusions; the H-range (4.1) and
   scale separations (3.2); the fourth-moment error terms and their normalisation in Thm 5.3;
   Lemma 4.8's composite-moduli large sieve; Lemma 6.1; Theorem 7.1 (iii) cost).
2. Quantifiers and populations: is every "almost all"/"average"/"worst case"/"every candidate"
   statement exact? Is anything overclaimed (abstract, intro, novelty list)? Is the
   distinction between the true-factor task and full-candidate covering kept everywhere?
3. Bit-complexity accounting: anything missing or wrong (detector, exponent trees, AKS,
   remainder trees over Z/N, exploration in DCL-wc)?
4. Prior art and attribution: verify (with web search) the claims about McKee–Pinch §5, Harvey
   2021 Alg. 4.1/4.2, Harvey–Hittmeir 2022 Prop. 3.3, Harvey–Hittmeir 2026 (arXiv 2601.11131)
   Lemmas 2.1/2.3/2.4/§3.1, Ward–Maestro (Def. 2.1, Lemma 3.2, Thms 3.4, 3.6, 9.4, Props 14.2,
   35.3), Friedlander 1989, Carella §5's normalisation defect as stated in the paper,
   Umans–Wang, He–Sahai, Gao et al. 2512.19076, Peled–Yadin–Yehudayoff 2011, AKS. Flag any
   misattribution or bibliographic error in refs.bib. Search for any further prior work that
   pre-empts Theorem 1.1, 1.2, 1.3 or 1.5 (average/almost-all deterministic factoring below
   1/5 for balanced semiprimes; true-factor lower bounds in collision models).
5. Presentation for a Math. Comp. submission: structure, what to cut or move to appendices,
   missing definitions, notation clashes, LaTeX issues.
6. Final verdict: ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT, and the 5 most
   important fixes in priority order.
Be precise and terse; number your findings.
