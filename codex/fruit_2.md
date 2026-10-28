# fruit_2: adversarial referee of notes/fruit_results.md

You are an adversarial referee (analytic/algorithmic number theory). The repo /Volumes/git/factoring-connection is
read-only for you; scratch files go in your cwd. Read notes/fruit_results.md in full, plus what it cites:
paper2/main.tex §4 (Theorem DD) and §2 (SPLIT), notes/nfsfruit_results.md §§3–4, notes/lemmaM_results.md §§1–5,
src/fruit/theory.py, src/fruit/linalg.py, src/fruit/nfs_hybrid.py.

For EACH claim labelled proved (Lemma HY, Lemma DD′, Theorem HYB and "Choosing S", Lemma ABS-dom, Prop T, Prop SR (a)–(d),
Prop IV, Prop IV-LV and its consequences, Cor IV-TV) give a verdict CORRECT / CORRECT WITH FIXES / WRONG, with the exact
flaw and a fix. Try hard to break:
1. Theorem HYB: is everything except the hidden roots really F-measurable in an LV-type algorithm (characters, the
   exceptional-set pairing, rejection of T ∩ D = ∅, LV's non-monic f_d factors, the f'(α) factor)? Is the claim
   "for every N with ω ≥ 2" right (prime powers? small factors? non-unit R?)
2. Prop IV-LV: check the counting (squarefree part divides Δ; one m per (d,u); the box; signs; LV's exact m-interval
   and box; does LV's x actually contain rational parities and sign?). Is the conclusion "x injective on X" right, and
   does it kill criterion (I) and fibre mixing η* as claimed? Any loophole (e.g. involutions defined on batches, on
   (field, pair) jointly, or allowing x(ιz) − x(z) in a fixed subspace)?
3. The descent ladders (k+2)/(3k+3) (hybrid into an LV field) and (k+1)/(3k+2) (pure): check the model derivation,
   and the claim that the descended row must live in LV's field/factor base (is there a two-field or MNFS escape?).
4. Prop HX(i): does SD₁ really give rigorous L[1/3] by pure descent without any LV fruitfulness? What else is needed?
5. Prop T: is the formal statement correct and non-vacuous? Is the "distributional vs hidden-root" dichotomy honest?
6. Prop SR(d): Rabin-type argument. Any subtlety?
7. Are any of these already in the literature (Dixon, Pomerance, Buhler–Lenstra–Pomerance, Couveignes, Montgomery,
   Lee–Venkatesan, Seysen, NFS-DL descent: Gordon, Joux–Lercier, Commeine–Semaev, Barbulescu, Guillevic)? Cite precisely
   or mark UNVERIFIED.
8. Overclaims or underclaims in the Verdict.

Output ≤ 3000 words: numbered findings with severity (FATAL / MAJOR / MINOR) and concrete fixes.
