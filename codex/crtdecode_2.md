# crtdecode round 2 (ADVERSARIAL): referee notes/crtdecode_results.md (v1, commit 5cb24e0)

You are a hostile but fair referee (coding theory, lattice algorithms, analytic number theory, complexity). The repo
/Volumes/git/factoring-connection is read-only for you; write scratch files in your current directory. Python:
/Volumes/git/factoring-connection/.venv/bin/python (sympy, gmpy2, python-flint, numpy). Web search is on.

Read, in this order:
- notes/crtdecode_results.md (the note under review), notes/crtdecode_prereg.md (pre-registration, commit 7d07797);
- src/crtdecode/*.py, tests/test_crtdecode.py, experiments/raw/crtdecode_E{1,2}.jsonl;
- context: notes/sinonobliv_results.md §2 (Lemma IC, Props JC/AT, Prop RD), notes/sihard_results.md §5.2 (Theorem J),
  your own blind round codex/crtdecode_1.out.

Tasks.
1. **Proofs.** Check every PROVED item: Prop F1, F2 (a)(b)(c) — especially the claim that the code-specific capacity
   threshold is ν = 2/(3κ) and that the uniform q-ary formula gives ν = 1/κ; F3 (t/k² → κ/ν², ln n/(k ln k) → 3κ²/ν,
   = 1/2 and 1 at κ*, ν*); F4 (minimum-weight reading); Prop DT (correctness of the tree, the exact leaf rule, the
   depth bounds, worst case and expectation); Lemma NU; Lemma GC; Prop PT (the sup over d ≥ h/2, the use of Prop JC's
   case D > V, the scope remarks); Prop SG (Λ_J, the injection, Ψ(X,(ln X)^A) = X^{1−1/A+o(1)}, A → 4/3, uniformity;
   is the weighted case covered? is "full agreement" the right reading of "pure decoding"?); Prop PNP; Lemma SC.
   Give counterexamples where something fails.
2. **Overclaims.** Is "Prop SI sits exactly at the code-specific capacity" justified, and is "capacity" the right word
   for an average-list threshold of a fixed code? Is "the natural adaptive model cannot carry a barrier" stated with
   the right scope? Are the labels (PROVED / MODEL / assessment) right? Is anything presented as new that is standard?
3. **Literature.** Check the citations and claims in §1.3–§1.5 and §3.2: Cheng–Wan (what exactly do they prove, which
   regime, which venue/year), Kiayias–Yung, Naor–Pinkas, Bleichenbacher–Nguyen (noisy CRT), Gandikota–Ghazi–Grigorescu,
   Guruswami–Vardy, BKR, CHHS, Erickson 1999 (r-linear decision trees for r-SUM; does the transfer statement make
   sense?), Kane–Lovett–Moran, Dyer's knapsack FPRAS, de Bruijn's Ψ(x, log^A x), BGJT 2014. Re-check BCPZZ Thm 4.1 and
   Jeronimo against the arXiv texts.
4. **Experiments.** Reproduce E1 and E2 (`python -m src.crtdecode.run_exp E1`, `E2`, `python -m src.crtdecode.analyze`;
   ~1 min total, run them into your scratch dir only if you can redirect output; otherwise re-score from the committed
   JSONL), re-score P1–P5 against the pre-registration, check the decoders for leakage or accounting bugs (cost
   counting, budget checks, the PH threshold, the information decoder using X only through the oracle), and check the
   E1/E2 tables in the note.
5. **Missed constructions.** Is there an adaptive decoder, in any model you consider natural, that the note's ladder
   (§2.2) misses and that plausibly beats the divisor trick by more than L^{o(1)}? Be concrete.

Output: a referee report with a verdict (ACCEPT / MINOR / MAJOR REVISION), then a numbered table of findings
(severity FATAL/MAJOR/MINOR/NIT, location, finding, suggested fix). ≤ 3500 words.
