You are an independent expert (computational + analytic number theory) making a BLIND attempt at a research question.
The repository at /Volumes/git/factoring-connection is read-only for you; write any scratch scripts under your current
directory. Rigour matters far more than speed.

## Background (internal, unrefereed results of this project)

- `notes/dlog_results.md` and `paper2/main.tex` §16 (`\section{Application: discrete logarithms ...}`, label sec:dlog):
  **Theorem DL**, a Las Vegas L_p[1/2, c] discrete-log algorithm in F_p^* for every c > 1/√2, obtained from the relation
  supply of the project's factoring Theorem A with N := p, plus **Proposition XK** (exact integer dependencies: exact
  cancellation of rational and prime-ideal exponents, units removed by archimedean logs + LLL + the Blanksby–Montgomery
  gap, hidden exponents T_i = g^{a_i}h^{b_i}; no factorisation of p − 1, no Schirokauer maps). Fixed degree d.
- `notes/lv5unif_results.md` (audited in `notes/auditU_audit.md`): **Theorem U**, Las Vegas factoring of every N in
  L_N[2/5, c], obtained by re-running the relation-supply chain at growing degree d ≍ D^{1/5}, k ≍ D^{2/5}
  (D = log N/log log N), B = L[2/5, a], special-q primes L[3/5], with all precision losses polylogarithmic ("τ = 0"),
  and a recursive shape test (factor the Gauss-split halves |u₀|, |v₀| ≤ √(2N) by the algorithm itself; Lemma F_α).
- `notes/recurse_results.md` (Lemma S″, Lemma F, §3.3 size model), `notes/thmA_proof.md` (roadmap of Theorem A).

You may treat Theorem U's internal chain (FM-U, LD*-U, Kern-U, Lemma 8.5′(R′), CS, D-U, K-U, select-U) as given
**for composite or prime N**, but check that every use of N's arithmetic survives N := p.

**Blindness rules.** Do NOT read STATUS.md, any file named `notes/dlog25_*`, `src/dlog25/*`, `tests/test_dlog25*`, or any
`codex/dlog25_*` file other than this prompt.

## Question

Is there a rigorous (conditional only on the same internal chain as Theorem U) Las Vegas algorithm computing discrete
logarithms in F_p^* for every prime p in expected time L_p[2/5, c] with c finite ("Theorem DL25")? If not, what is the best
rigorous α you can justify?

## What to do

1. **Transfer of Theorem U's design to DL.** With N := p, check each ingredient at growing degree: FM-U, LD*-U, Kern-U,
   Lemma 8.5′(R′), Lemma CS, D-U, K-U (split counting for T = g^a h^b uniform on F_p^*: what replaces Lemma S″, which is
   stated for R uniform on the squares mod N?), select-U (rounds now succeed with probability φ(n)/n, not 1/2), field
   rejections. Where does N's arithmetic enter? Any factor exits?
2. **The shape test and recursion.** In DL the halves u₀, v₀ must be factored, not logged. Can the shape test call
   Theorem U's factoring algorithm as a subroutine, and what does that do to the cost recursion (is there still a
   fixpoint, or is the outer level non-recursive)? Is there any meaningful recursion on logarithms?
3. **Proposition XK at growing degree d ≍ D^{1/5}.** Archimedean unit logs (conjugate lower bounds, precision), the
   Blanksby–Montgomery gap (≈ 1/(d² log d)), the 2^{O(d²)} or d! factors in the relation-size bound, LLL in dimension ≤ 2d,
   root isolation of f, r-maximal orders at bad primes r | f_d·disc f, torsion μ_K (d odd), exact rational solving (Dixon
   lifting, Hadamard bounds on entry size), hidden-exponent uniformity. Is the total cost L^{3a + o(1)} with o(1) in the
   L_p[2/5] sense (i.e. exp(o((log p)^{2/5}(log log p)^{3/5})))? Does the cubic linear algebra term bind at the
   optimum of the L[2/5] cost formula (the fixed-d note says "the 3a term takes over" in some L[1/2] regimes)?
4. **Individual-log descent.** In this design every relation is a descended randomized target. Check the target
   randomization and the hidden-label independence when the shape test is itself a randomized recursive algorithm.
   Is a precomputation / per-log split possible at L[2/5] (and at what cost), or does every log cost the full run?
5. **Assemble** the theorem with exact parameters and give the cost constant as a formula in (a, κ, d₀, C_γ) and in the
   factoring constant used for the shape test; compare with Theorem U's c.
6. **Literature.** Best rigorous prime-field DL bounds known (Pomerance 1987 L[1/2, √2]; anything better, for all p, almost
   all p, or under GRH?). Extension fields F_{p^n} in medium characteristic: what is known rigorously (Bender–Pomerance,
   Kleinjung–Wesolowski, others), and would this transfer plausibly give anything there? Use web search; cite precisely
   and say what you actually read.

## Output

Markdown. Start with a one-line verdict (YES L[2/5] / YES WITH CAVEATS / ONLY α = … / NO). Then: the algorithm and the
theorem as you would state it; a table of ingredients (transfers verbatim / modified / replaced / fails) with reasons; the
cost formula; the weakest link; answers to items 1–6; any counterexample or gap you find. Be concrete; do not pad.
