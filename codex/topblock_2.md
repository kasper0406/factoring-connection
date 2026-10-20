# topblock_2: ADVERSARIAL REFEREE ROUND on notes/topblock_results.md

You are an adversarial referee (lattice reduction / Coppersmith methods, computational number
theory) with web search. The repository at /Volumes/git/factoring-connection is read-only context.
Write scratch files only under your current directory. Your blind answer is in codex/topblock_1.out.

Read notes/topblock_results.md IN FULL, notes/topblock_prereg.md, and src/topblock/{core,lattices}.py.
Raw data are in experiments/raw/topblock_*.jsonl.

Your job is to BREAK the claims. For each item, give a verdict: CORRECT / GAP / WRONG / OVERCLAIMED.
Back each with a proof, a counterexample, or an exact reference. You may run small exact computations
(.venv/bin/python has sympy, flint, fpylll).

1. **Prop T (transfer theorem), §3.1.**
   - Is (x+h)σ(w_b) ∈ D_b·(N, x+h) correct for the coefficientwise-reduced g_b and shifted centres
     x_b = x − Δ_b?
   - Does every spanning polynomial land in ∏D_b^{m2_b}·(N, x+h)^{m1+T}? Check the claimed degree
     bound d and the norm transfer constant (3 + 3√2)^d·(1 + X/x)^T.
   - Is the claim that it covers "every shift strategy" (polynomials of the form monomial × N^α(x+h)^s
     × ∏w^t D^{m2−t} with α + s ≥ m1) correct? What natural lattices are NOT covered: shifted windows
     for near-divisors, unreduced/other representatives, rational centres, polynomials vanishing only
     mod D (not mod p), sublattice tricks, non-Howgrave-Graham certificates?
   - Is the exponent-level consequence ("o(1) loss") correctly stated when d grows with log N?
2. **The F-artifact, §3.2.** Is the HG-shortness boundary of F·N^{m1−1}D^{m2−1} right? Is the
   "negative per-row contribution = F" explanation right, or only suggestive?
3. **Theorem EQ and Props KB / LOC / IV, §1.**
   - Is the corrected Prop LOC bound proved (including ν₀ ≥ 2 blocks and odd j)?
   - Is "the window oracle is never easier than promised search" justified, given the candidate
     factor?
   - Is "the ladder's target is reachable through some block primitive iff promised search costs
     N^{o(1)}" correctly stated?
   - Is Lemma TS correct and correctly interpreted?
4. **Lemma MI and §2.3.** Is the basic lattice equal to the full truncated ideal lattice? Are the AASW,
   CHHS, Cohn–Heninger and Herrmann–May statements accurate? In particular: does CHHS (ASIACRYPT 2016)
   also treat the unknown-divisor case b ≥ N^β? What exactly are AASW's lattice L_d and Assumption 1?
   Does our observation λ₁ ≪ GH (beyond the Howgrave-Graham threshold, for the univariate ideal
   lattice) contradict their Assumption 1, or is it outside its regime?
5. **§3.3 HM split.** Is the "full truncated ideal in total degree" claim right? Is the adjacency
   argument right?
6. **Experiments.** Do the data (experiments/raw/topblock_x3–x6.jsonl) support each pre-registered
   prediction as reported? Look for design flaws, e.g. comparators that make P3.1 tautological, or
   scan-protocol artefacts.
7. **Verdict and barrier statement, §0.** Is anything overclaimed? Give replacement text for any
   sentence you would change. Is there ANY mechanism by which the ladder's window structure could
   still beat partial key exposure that the note fails to exclude? If so, state it precisely.

Return an itemised report with verdicts, then a list of required fixes. Under ~3500 words.
