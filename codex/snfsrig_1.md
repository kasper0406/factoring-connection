# snfsrig_1 — blind round: rigorous sub-L[2/5] factoring for special-form N?

You are a research mathematician (analytic number theory and the number field sieve). Work independently. The repository at
/Volumes/git/factoring-connection is read-only for you; write any scratch files only into your current directory.
Rigour matters more than speed. Label every claim: proved / heuristic (size model) / conjectured / open / literature (with
exact citation).

## Background you must read first (repository paths)

- notes/lv5unif_results.md — Theorem U: a claimed rigorous Las Vegas L_N[2/5, c] factoring algorithm for every N
  (internal, unrefereed). Read §1 (design: random field A = (g₀, g₁, w, c); F_A = F̂_{g,w} + G·R_c with F̂(g₀,g₁) = N·w;
  special-q primes ℓ ≈ N^{1/(2k)}; box V ≈ √ℓ·S), §3 (Kern-U), §4 (FM-U, GM-U), §6 (LD*-U, Lemma W), §13 (why not
  below 2/5).
- notes/recurse_results.md §3.3–3.4 — the size model: (R1) split acceptance forces k = O(D^α), D = log N/log log N;
  (R2) smoothness; (R3) design floors θ = 2μ + h + (d−1)ν with h ≥ max(2ν, μ − (d−2)ν, 1 − dμ + ν), giving
  θ ≥ 3/(d+1) (window+level) and θ ≥ (d+1)/(4k) (Kern); (R4) box. (R1)+(R3) give α ≥ 2/5. §3.4: multi-level descent
  in one field reaches 1/3 in the size model.
- notes/nested_results.md §§3–5 — Prop SI (sparse image: fixed field ⇒ no density-only lemma), Lemma FB, Prop R,
  **Prop R2** (model: designs without algebraic children have α ≥ 2/5; premise (M1) is the counting heuristic
  (d+1)log‖f‖ + 2 log‖g‖ ≳ log N), Lemma OD.
- notes/ldbound_results.md "Verdict" and §3 — why the N·w multiplier is needed (the forced residue
  F(P) ≡ N(b/g₁)^d mod G(P) is spread by w via Brun–Titchmarsh, Lemma W).
- notes/thmA_proof.md §1 — the field law 𝒜‴ and the role of each random parameter.

## The question

Special NFS: for N = r^e − s with small r, s (Cunningham-type numbers), there is a FIXED polynomial pair
f(x) = r^{j}x^d − s (or similar), g(x) = x − m, m = r^{(e+j)/d}, with ‖f‖ = N^{o(1)} and ‖g‖ ≈ N^{1/d}. Heuristic SNFS is
L[1/3, (32/9)^{1/3}]. Smaller norms relax the design floors (R3); but the polynomial is fixed, so the field averaging
that all the project's relation-supply proofs (FM, LD*, Kern, Lemma W, Theorem R for the rational side) rely on is lost.

Goal: a rigorous Las Vegas algorithm with exponent α < 2/5 (ideally 1/3) for a natural infinite special family of N;
failing that, the best α and the precise obstruction.

Please address, with proofs where you claim them:

1. **Size model.** Redo recurse §3.3–3.4 for small-height f (‖f‖ = N^{o(1)}, degree d growing). Does the floor in (R3)
   move? Which α does one-level descent reach? Which does multi-level descent reach? **Recompute Prop R2** (designs
   without algebraic children) when (M1) fails because ‖f‖ = N^{o(1)}, allowing each descent node to choose its own
   degree d among the special polynomials available for the same N (for N = r^e − s, every d ≤ log N works).
2. **Recovering randomness.** For each device decide which of FM / LD* / Kern / Lemma W / Theorem R survives, and why:
   (a) a randomized rational side g = g₁x − g₀ sharing the root mod N (Buhler–Lenstra–Pomerance base-(m₁, m₂));
       twisting the root by a small multiplier;
   (b) multipliers kN or N·w with w random — can r^e − s times w keep a small-coefficient polynomial in a NEW field?
   (c) twists f(x+t), u^d f(x/u), Möbius transforms;
   (d) families of SNFS polynomials for the same N (several (d, j), several representations);
   (e) random special-q choices (averaging over ℓ only);
   (f) small perturbations f + c(x)·g(x) with c of height H.
   Is there any proof that the polynomials sharing a known root with tiny height are essentially unique (so that any
   field randomisation costs height ≥ N^{1/d})? What is the cleanest rigorous form?
3. **Special families with natural rigour.** N = r^e ± 1 with cyclotomic structure (algebraic factors Φ_d(r),
   Aurifeuillian factors). What does this buy rigorously (exponent, or only constant)? Any other family where a
   rigorous L[α < 2/5] is natural (e.g. N with known small-height polynomial representation; N whose prime factors
   have special structure)?
4. **Literature** (search the web): rigorous or conditional results for SNFS, special-form factoring, rigorous index
   calculus in F_p for special p (Semaev 2002 "special primes"? Gordon 1993; Joux–Lercier 2003; SNFS-DL), and any
   rigorous L[1/3] (or sub-1/2) result for special forms, in factoring or prime-field DL. Give exact citations and say
   what is actually proved.
5. **Verdict.** Best rigorous α you can justify for a natural infinite special family, and the precise obstruction
   (state it as a mathematical statement that would suffice, e.g. a pointwise smooth-values estimate for a fixed binary
   form at L-scale, and compare with the best known results on friable values of binary forms / polynomials).

Write your answer as a structured report (≤ 4000 words), verdict first.
