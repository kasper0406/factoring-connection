You are a blind, adversarial referee (analytic number theory: smooth numbers, Bombieri–Vinogradov/BDH-type estimates,
lattices; and computational number theory). The repository at /Volumes/git/factoring-connection is read-only for you;
write any scratch scripts under your current directory.

## The claim to referee

`notes/audit2L_blind_claims.md` (an excerpt of a chain of internal notes, with review sections, verdicts,
probabilities and reviewer attributions removed) claims **Theorem 2L-final**: a rigorous Las Vegas algorithm that
completely factors every integer N in expected time L_N[α, c(α)] for every fixed α ∈ (3/8, 2/5), and also at α = 3/8
under a design constraint a > 2c₂. This would be a major result: the best published rigorous bound for arbitrary N is
Lenstra–Pomerance 1992, L[1/2, 1+o(1)]. The chain is conditional on an internal, unrefereed one-level chain giving
rigorous L[2/5] ("Theorem U": `paper2/main.tex` and `notes/lv5unif_results.md`). **Treat Theorem U's chain as given**
(it has been audited separately); your job is the two-level layer on top of it: Hyp U′ in design DW, Prop UB, Prop
DC/DC_e, Lemma η, Prop PB, Lemma VD, Theorem RAT (rational two-point blocks), Theorem AB/AB₂ (algebraic two-point
blocks), Prop CLB, and the assembly/cost (excerpt A §13).

**Blindness rules.** Do NOT read STATUS.md, README.md, REPORT.md, BARRIERS.md, NEXT_STEPS.md, any notes/*audit*.md
(except notes/audit2L_blind_claims.md), notes/mf2L_*, or any codex/* file other than this prompt. You MAY read the
source notes the excerpt names (notes/{thmAB,c12rat,c12alg,pc3,uprime,multidescent,lv5unif,nested,escapeR2,lem85}_results.md,
paper2/main.tex) and src/{thmAB,c12rat,c12alg,uprime}/*.py, **but not their review sections ("Reviews"), their
"Verdict"/"0. Verdict" sections, or any probability estimate**, before you have written your verdict. The web is
allowed (Drappeau Compositio 2015 = arXiv 1307.7554; Harper arXiv 1208.5992; Hildebrand JNT 1986;
Hildebrand–Tenenbaum 1986 and the 1993 survey; de la Bretèche–Tenenbaum 2005; anything else).

## What to do

Try hard to break the claim. A single polylogarithmic loss in the wrong place is fatal: the two-level precision budget
needs the incidence-weighted relative variances A (child counts) and B (parent counts) to be ≤ D^{−α}λ^C and
≤ D^{7α−3}λ^C in the open window (D = log N/log log N, λ = log log N), and ≪ D^{−3/8} at α = 3/8, while the
"first-order" error of every smooth-number asymptotic used is ≍ D^{−α}. So the proofs claim **second-order** bounds
D^{−2α}λ^C via a Cauchy–Schwarz cell-covariance device. Check in particular:

1. **The cell-covariance device** (RAT Prop PD, AB Prop PD-A): exact cell factorisation (Lemmas TYPE, CP, CL, CF, FIB),
   the covariance identity with the atom ∗, and whether the marginal masses really equal cell densities on interior
   cells (fibre-length weighting!). Is every first-order quantity really entering only squared? Is any polylog factor
   multiplying a merely D-polynomially small quantity? Check Lemma OPV / OPV-A (one-point second moments ≤ D^{−2α}λ^C,
   including the local factor, prime powers, β-drift, the archimedean factor, near-zero and large-gcd exceptional
   sets) and the lower bound p₁ ≥ c₁/2 (OPV-A step 6).
2. **Character errors** (Theorem WA, Lemmas EN/EN-A, FAM/FAM-A, DOM, SPR, MULT): is Drappeau's Lemme 5 used within
   its hypotheses (conductor ≤ x^{η_D}, (log x)^{c_D} ≤ y, its exact normalisation and the 1/φ(k)-weighting), is the
   domination quantity 𝔇 really e^{O(D^θλ²)} with θ = 2α − 1/2 < α (children) and θ_p = max(1−2α, 2−5α) < 3−7α
   (parents), and is the comparison with the label weights (max_q vs Σ_q over a dyadic window; modulus k containing q
   or not) legitimate? Are the pointwise weight bounds (DOM(i)) and spread bounds (SPR) true?
3. **Thin pairs and skewed labels** (Lemmas TF, LS, PUB, PU-A; Prop SK; AB §10): fractions, losses, and whether the
   per-point bound for skewed labels really avoids any cross-box comparison.
4. **Parents** (AB §§9, 11; RAT RAT-2): separate thresholds, Lemma W-P (P1)–(P3), the hard filter and π_P-centring,
   the root-conditioned acceptability subtraction, and the parent "weight bridge" to RAT (AB §13 step 1).
5. **The assembly**: Prop DC/DC_e (incl. environment cells), Lemma VD, Lemma η (W/μ₁ ≤ 3u₂ with no polylog loss),
   Prop PB, Prop UB, U′ under DW (c12alg §2.7, including the "routine range restatement" and, at 3/8, the
   leading-order re-check at η₀ = D^{−3/8}), Prop CLB, field rejections, split acceptance, rows/LV6/fruitfulness,
   recursion (Lemma F_α), and the cost formula c(α). Is the window really 3/8 < α < 2/5, is a > 2c₂ the right
   endpoint constraint, and are the exponent gaps (½ − α for children, 2 − 5α for parents) correct?
6. **Sanity.** If this works, why is it not known? Is there any known obstruction to two-level descent at L[3/8] with
   rigorous analysis?

## Output

Write your report in markdown (it is captured automatically). Start with a one-line verdict for each of
(i) α ∈ (3/8, 2/5) and (ii) α = 3/8: ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT. Then a table of issues
(FATAL / MAJOR / MINOR, location, why, whether it changes the exponent or only the constant, fix if any), then the
single **weakest link**, then answers to items 1–6. Be concrete: quote the line, give the counterexample or the exact
estimate that fails. Do not pad; do not accept "verbatim as …" without checking the cited proof yourself.
