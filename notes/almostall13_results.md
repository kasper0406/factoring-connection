# almostall13: averaging over N to supply parameters for algebraic children — results

Agent `almostall13`, phase 14, 2026-10-20.

Files:
- This note; pre-registration `notes/almostall13_prereg.md` (commit `246ac18`); code frozen at `a94e61c`.
- Code: `src/almostall13/rank.py` (E1, exact), `harness_e2.py`, `analyze_e2.py` (E2, harness/diagnostic; imports
  `src/nested/toy_e1.py` unchanged), `bounds.py` (exact helpers for Lemmas M1-AA and MP and for the NR-AA constant).
  Tests: `tests/test_almostall13.py`.
- Data: `experiments/raw/almostall13_e1.json`, `almostall13_e2.jsonl`, `almostall13_e2_summary.json`; config
  `experiments/configs/almostall13_e2.json`.
- Codex: `codex/almostall13_1.{md,out}` (blind round), `codex/almostall13_2.*` (adversarial referee, §9).

Base texts: nested (`notes/nested_results.md`, cited as [nested]), lv5unif/auditU (Theorem U, [U]), recurse §3 ([rec]),
polyimpl §1 and §7 (transfer lemma). Everything about the project's own algorithms rests on the internal, unrefereed
chain. Labels: **proved**, **proved conditional on …**, **model computation**, **heuristic**, **experimentally
supported**, **open**.

## Verdict

**Plain language.** We found no way for averaging over the input N to supply the missing parameters, and this note
proves no improvement from it. Assuming the internal chain, **the best bound established for almost all N remains
L[2/5]**, from Theorem U, which already covers every N. We refute the specific idea that N's digits add about d + 1
free parameters per field. This is a refutation of that parameter-counting argument, not an impossibility theorem: it
does not exclude every benefit a different field law might bring. The reason is short, for the fixed-polynomial
special-q descent of [nested] §8:
- Once the field (F, G, w) is fixed, N is fixed too (N = F(g₀, g₁)/w).
- The descent inside a field never looks at N, except through gcd checks (Lemma DF).
- So averaging over N is averaging over fields under a different law on the same space of polynomial pairs.
- N's digits move only one coordinate of F, its value T = F(g₀, g₁) at the common root, and the multiplier w already
  randomises that coordinate. In the plain base-m design N's digits are a *substitute* for the random coefficients c,
  not an addition to them.
- Nested's per-field budget of d + 1 parameters is a property of a degree-d form, whatever its law.

Designs that use arithmetic modulo N inside the descent are outside this scope: escapeR2's re-split (Prop RS there),
for example, Gauss-reduces w·s² mod N. They have their own obstructions (notes/escapeR2_results.md).

Precisely:
1. **Prop NB (no budget gain) [proved].** For any fixed-polynomial special-q design with N | F(g₀, g₁), in the
   contracting regime where every scanned value has 0 < |value| < N, the following hold.
   - Every descent statistic is a function of the field alone, up to factor exits (Lemma DF).
   - An average over (N, field) equals a field average under the mixture law μ̄ (pushforward).
   - At s fixed points, the values have rank ≤ min(s, d + 1) under every joint law, with law-independent integer
     relations when s > d + 1.
   - In design (a‴), under explicit nondegeneracy assumptions, varying N adds 0 affine dimensions beyond w, and 1
     beyond fixed (N, w).
   - Lemma OD, Prop SI, the information version and the constrained-resampling budget hold verbatim (Cors OD-AA,
     SI-AA). E1 confirms the dimension counts exactly.
2. **What averaging over N does buy [proved]. None of it moves α.**
   - (a) Lemma W becomes unnecessary: N itself spreads the forced residues, Pr ≤ 1/t + 1/X (Remark W-AA).
   - (b) **Lemma M1-AA.** For all N outside a set of density X^{−ε}, every polynomial pair with a common root mod N
     obeys the counting premise (M1) of [nested] Prop R2. The constant is ½ for all N and 1 − o(1) for B-rough N.
     Here d ≤ (log X)^{0.9}. So R2's counting premise is *proved* for almost all N. R2 itself stays a model
     computation, because its premises M2–M4 are models. SNFS-type pairs at degree d exist only on density
     ≤ X^{−1/2+2/d+o(1)}. If anything, the almost-all setting is the less favourable one for the size model.
   - (c) An average-case AD implies AD for most N, by Markov (Lemma AV).
3. **Recursion does not tolerate an "o(1)" exceptional set (Prop HR) [proved].** The recursive sub-inputs live in thin
   sets of density e^{−Θ(E)}: the accepted halves u₀, v₀, which are products of k primes in a window, and the
   algebraic values F(P). An exceptional *measure* bound of density o(1), or even e^{−cE} for small c, is consistent
   with all of them being exceptional. A recursive almost-all design therefore needs a *hereditary* hypothesis. With
   one (closure under divisors as well as under descent values and split halves), **Theorem NR-AA** holds [proved
   conditional on Hypothesis HAD + LV6′; no ERH]: a Las Vegas algorithm, alternating with FACTOR_U, completely
   factors every N.
   - Its expected time is L[2/5] for all N.
   - For N in the tame set 𝒦 it is L_N[1/3, c] for every c > c_NR-AA.
   - The transfer to semiprimes needs ε_𝒦 = o(n^{−2}).
   - HAD follows from nested's worst-case AD, so NR-AA also removes ERH from nested's Theorem NR.
4. **ERH and characters [proved, routine].** ERH is not needed. Unconditional factorial-bit character primes ([U] §9)
   cost L[1/3, C_χδ] at d = δD′, which only raises the constant. Averaging over N gives no growing-degree family
   zero-free region (none is known).
5. **Lemma MP and the union-bound route [proved; method barrier, not impossibility].**
   - Moment information accurate to relative error η cannot certify Pr[dead] < η.
   - So a uniform per-target union bound over e^{b₁Q} targets, run from approximate moments alone, needs
     smooth-correlation asymptotics with relative error ≈ e^{−b₁Q} = L[2/3]^{−Θ(1)}. Methods that use extra
     structure (support constraints, exact identities) are not covered.
   - Known L-scale errors are of the shape exp(−cu/log²u) + y^{−c} on average over moduli, with u ≍ D′. That is far
     too weak. None of this depends on N.
6. **Literature (§7).**
   - No rigorous almost-all-N factoring below L[1/2] was found. The explicit rigorous almost-all-N complete
     factorisation result is Hafner–McCurley 1989 (pp. 544–545): exp((log N)^θ) for θ > 5/6, far above L[1/2].
   - The almost-all-prime ECM bounds (Pomerance 1987 Thms 2.1/B′ with exponent 6/7, re-used in LP92 §7 Thms 7.3–7.4;
     HM Thms 3.1–3.2) are much weaker than L_p[1/2].
   - The worst-case rigorous elliptic-type bound is Lenstra–Pila–Pomerance's L_p[2/3] (1993/2002), for every p.
   - Lee–Venkatesan average over the algorithm's coins, not over N.
   - Kleinjung–Wesolowski is a worst-case result for every small-characteristic field and uses Frobenius structure
     that integers lack.
7. **Experiments.**
   - E1: all four predictions hold exactly (affine dimension d with (N, w) fixed; d + 1 otherwise, so N adds nothing
     beyond w; 0 kernel violations).
   - E2 (toy, d = 3, 44-bit N, 72 + 72 fields): all four predictions hold.
     - Q1: random-N versus fixed-N mean algebraic descendability 0.285 vs 0.270, CI [−0.028, 0.057].
     - Q2: between-field SD ratio 1.14. N-averaging does *not* reduce field heterogeneity.
     - Q3: no statistically detected N-group effect within FIX in this toy (permutation p = 0.28).
     - Q4: nested's Jensen gap persists, J_pool = 1.23 / 1.30, and per-field rates explain it, J_field = 0.90 / 0.91.
8. **Reviews (§9).**
   - GPT-6.1-sol's blind round (`codex/almostall13_1`) reached the same verdict independently: "the supported exponent
     remains 2/5 … no extra coefficient rank beyond w, c".
   - The adversarial round (`codex/almostall13_2`, MAJOR REVISION; one FATAL issue, the heredity gap in NR-AA) and
     two independent Opus critics (proof: no FATAL, three MAJOR in NR-AA; literature: attribution corrections) were
     all applied.

---

## 0. Setting, notation, quantifiers

t := log N, E := t^{1/3}(log t)^{2/3}, Q := t^{2/3}(log t)^{1/3}, D′ := (t/log t)^{1/3}; L_N[1/3, c] = e^{(c+o(1))E}.
B = e^{bE} throughout §§2–3 (at α = 1/3).

A **fixed-polynomial special-q design** (the family of [nested] §8 and [U]) samples, on input N, a field datum
A = (F, G, w). Here:
- F ∈ Z[x, y]_d is a binary form;
- G = g₁x − g₀y with gcd(g₀, g₁) = 1;
- F(g₀, g₁) = N·w. This is design (a‴); w = 1 gives design (a″_k).

The design then descends targets by scanning special-q boxes of that fixed pair. More generally only N | F(g₀, g₁) ≠ 0
is used.
- Public parameters (d, B, box radii, M, H, …) are taken to depend only on the bit length of N. This is a design
  convention, needed to integrate over N at fixed g; [U] defines M, H, W, X_F from N itself.
- Descent (DESC), targets 𝒯(A), boxes and "good points" are as in [nested] §8.
- Designs that do arithmetic modulo N inside the descent (escapeR2's re-split, Prop RS there) are outside this family.

**Definition AA (almost-all Las Vegas bound).** A single uniform algorithm 𝒜, fixed once and for all, *factors almost
all N in expected time L[α, c]* if the following hold.
- (i) On every input it halts with probability 1 and every output is a certified complete factorisation.
- (ii) There are η(X), ε(X) → 0 such that for every large X,
  #{N ∈ [X, 2X) : E_ω T_𝒜(N) > exp((c + η(X))Λ_α(X))} ≤ ε(X)·X.

Remarks:
- The exceptional set E_X need not be computable, and 𝒜 does not know whether N ∈ E_X.
- Exceptions concern expected time only, never correctness.
- Running 𝒜 in alternation with a worst-case algorithm (here FACTOR_U of [U]) keeps (ii) and makes the expected time
  finite for every N.
- A high-probability statement jointly over (N, ω) does not by itself imply (ii) (codex).

**Lemma TR (transfer) [proved; polyimpl §7 D1, detqs §1].** Let S_n be the n-bit balanced semiprimes, |S_n| ≍ 2^n/n².
- If the bad set in [2^{n−1}, 2^n) has size ≤ ε_n2^{n−1}, then the bad fraction of S_n is O(ε_n n²). So
  ε_n = o(n^{−2}) gives "all but o(1) of balanced semiprimes". Density merely o(1) can contain every balanced
  semiprime.
- An almost-all bound is not a mean bound. With an L[2/5, c_U] fallback, the worst-case envelope bounds the mean over
  [X, 2X) if ε(X) ≤ L_X[α, c]/L_X[2/5, c_U]. This is a sufficient condition.

Theorem U satisfies Definition AA with E_X = ∅ at α = 2/5. That is the best almost-all bound established here.

## 1. Averaging over N is averaging over fields (Prop NB)

**Lemma DF (descent depends only on the field) [proved, in the contracting regime].** Assume that every value scanned
by DESC satisfies 0 < |G(P)|, |F(P)| < N. This holds for large N, since log|G(P)| = O(Q) and log|F(P)| ≤ τ_F t + O(Q)
with τ_F = δb₁ < 1 ([nested] §8). Let DESC°(A, τ) be DESC with D_A° := g₁·f_d·disc(f) instead of
D_A = N·g₁·f_d·disc(f), and without gcd checks against N. Then for every N with N | F(g₀, g₁) and every coin sequence of
the recursive calls, DESC(N; A, τ) either makes a factor exit or returns exactly DESC°(A, τ): the same tree, the same
success indicator. The target sets satisfy 𝒯°(A) ∖ 𝒯(A) ⊆ {primes dividing N}.

*Proof.* Inspection of [nested] §8.
- The box comes from the Gauss-reduced basis of the target lattice, which is fixed by (q, G) or by (q, r) with
  F(r, 1) ≡ 0 mod q.
- The values are G(P) and F(P), and their canonical factorisations are functions of the values (and of the coins of
  the recursive calls).
- The good test compares prime factors with B, with e^{λz} and with D_A.

N enters only in two places:
- through gcd checks, where by the size assumption a nontrivial gcd is a proper divisor of N, hence a factor exit
  (gcd = N cannot occur; codex's example with F(17, 1) = 15 = N is excluded by |F(P)| < N);
- through the primes p | N in D_A. If such a p divides a value, the gcd check has already exited.

Removing N from D_A can only add targets, namely primes dividing N. Rational children sent to fresh fields
([nested] Prop R) are covered by applying this to each field. ∎

**Proposition NB (no budget gain) [proved].** Let ν be any law on N (e.g. uniform on [X, 2X)), μ_N the design's field
law, and μ̄ := Σ_N ν(N)μ_N. In general μ̄ is not the uniform law on pairs of given heights.
- (i) *Pushforward.* For every event ℰ that is a function of DESC° runs in one field,
  Pr_{N∼ν, A∼μ_N}[ℰ(A)] = Pr_{A∼μ̄}[ℰ(A)], up to factor exits.
- (ii) *N is a function of the field.* Given (F, G, w), N = F(g₀, g₁)/w. If only N | F(g₀, g₁) is used, N is one of at
  most τ(F(g₀, g₁)) divisors. DESC° does not use N at all.
- (iii) *The budget is intrinsic.* For points P₁, …, P_s and any joint law of (N, A), the vector (F(P_i))_i equals
  E_P·(f₀, …, f_d)^T, with E_P[i][j] = a_i^j b_i^{d−j}. So it lies a.s. in a Z-module of rank ≤ min(s, d + 1) that
  does not depend on the law. For s > d + 1, the left kernel of E_P (integer relations λ depending only on the points)
  annihilates it a.s. Jointly with (G(P_i))_i the rank is ≤ min(s, d + 1) + min(s, 2). These are upper bounds; a law
  may have smaller support.
- (iv) *What N adds.* Fix g. Then {F : F(g₀, g₁) = T} = F̂_{g,T} + G·Z[x, y]_{d−1}. In design (a‴),
  F = F̂_{g,Nw} + G·R_c, so (N, w, c) ↦ F factors through (T, c) with T = Nw, and N changes only the coordinate T.

  Assume the points are pairwise non-proportional, G(P_i) ≠ 0, the c-box spans Z^d, and T takes at least two values.
  Then the affine dimension of the support of (F(P_i))_{i≤s} is:
  - min(s, d) with (N, g, w) fixed;
  - min(s, d + 1) as soon as w or N varies.

  So N adds 0 dimensions beyond w, and 1 beyond fixed (N, w). In the plain base-m design (no c, no w), the d + 1 digits
  of N *are* the coefficients: averaging over N supplies at most the d + 1 directions that (c, w) together supply in
  (a‴); c alone supplies d.

*Proof.*
- (i): Lemma DF and the definition of a mixture.
- (ii): trivial.
- (iii): F(P_i) = Σ_j f_j a_i^j b_i^{d−j}, so the rank is ≤ rank E_P; the G-part has 2 coefficients.
- (iv): if F(g₀, g₁) = F̂(g₀, g₁) = T, then F − F̂ vanishes at the primitive point (g₀ : g₁), so G divides it in
  Q[x, y] and, by Gauss's lemma, in Z[x, y]. In `make_f`, F̂ = gbe(Nw, …) depends on (N, w) only through Nw.
  - Under the stated assumptions the c-part G(P_i)R_c(P_i) has rank min(s, d): binary forms of degree d − 1 evaluated
    at non-proportional points have Vandermonde-type rank, and G(P_i) ≠ 0.
  - The T-direction is transverse to it whenever s ≥ d + 1: evaluation of degree-d forms at d + 1 distinct points is
    injective, and F̂_{g,T} − F̂_{g,T′} ∉ G·Z[x, y]_{d−1} for T ≠ T′. ∎

**Corollary OD-AA [proved; the children count is nested's typical-factorisation model].** [nested] Lemma OD holds
verbatim for the joint (N, field) law: with s > d + 1 fixed points the values satisfy fixed linear relations.
- Constrained resampling F + G·∏_{i≤k}L_{P_i}·H, with deg H = d − 1 − k, has d − k free directions during a run, because
  N is fixed by the field. Releasing the evaluation constraint across runs gives at most one more, d + 1 − k.
- In the typical-factorisation model ([nested] §5) a node has ≈ d/(2λ) large algebraic children, and that count is
  unchanged. So the model's depth-1 overdetermination for λ < 1/2 persists.

**Corollary SI-AA [proved].** [nested] Prop SI(i) is a fixed-field statement and holds pointwise in (N, A): the family
𝒮_{N,A} := 𝓘 ∖ Φ_{N,A}(Ω) has density 1 − o(1) and is hit with probability 0 for every (N, A), so also on average.
SI(ii), full support for *field-independent* families, is the only averaging route, and averaging over N is a special
case of it (NB(i)). SI(iii) is unchanged: the good family of an algebraic child is "descendable in A", a function of A.

**Corollary Info-AA [proved].** Given g and w, after d + 1 revealed values at points with independent evaluation rows,
F is determined, and with it N = F(g₀, g₁)/w. So beyond depth 1 every needed statement is pointwise in the field,
exactly as in [nested] §5's information version.

*Answer to task 1(a).* The idea "N's digits add ≈ d + 1 parameters" is false in design (a‴): c already randomises the d
directions of the fibre {F(g₀, g₁) = T}, w randomises T, and N moves only T. It is true in the plain base-m design, but
there the digits are the whole budget. Across r fields sharing g, N contributes one *shared* direction: with fixed
multipliers w_j the joint coefficient space has dimension ≤ rd + 1. With independent w_j, as in (a‴), it is r(d + 1),
but the extra directions come from the w_j, not from N (codex referee: rank 8 for two cubic fields).

## 2. What averaging over N does buy

**Remark W-AA (Lemma W is replaced by N) [proved for the congruence; proved with the stated changes for its use].**
Fix g, P = (a, b) and t := |G(P)| with gcd(b·g₁·w, t) = 1 ([pp] Lemma coprime). Since a ≡ g₀bg₁^{−1} (mod t) and F̂ is
homogeneous of degree d,
  F̂(P) ≡ (bg₁^{−1})^d·F̂(g₀, g₁) = (bg₁^{−1})^d·Nw (mod t),
and F_A(P) ≡ F̂(P) (mod t) as well. The map N ↦ (bg₁^{−1})^d w·N is a bijection mod t, so for N uniform on X
consecutive integers,
  Pr_N[F̂(P) ≡ r (mod t)] ≤ 1/t + 1/X,  hence  E_N 𝔔₀ ≤ (1/t + 1/X)·𝔔(t|Δ|) ≤ ((1 + t/X)/φ(t))·𝔔(t|Δ|)
in LD*-U Step 4. Step 2's AM–GM needs only marginal residues, so this substitution is legitimate.

This is Lemma W-U's conclusion, with C_W = 1 + o(1) instead of 24/η₀, and with no Brun–Titchmarsh. Conditions:
- the public parameters must be functions of X (§0), and N-dependent rejection indicators must be dropped before taking
  the upper bound;
- the algorithm splits the B-rough part N′, not N; for N′ uniform on the rough integers of an interval a trivial count
  gives C_W = O(log B).

To conclude "most N", apply Markov to the per-N second moment, then Paley–Zygmund with FM-U's per-N first-moment lower
bound. The exceptional density is 1/K with γ ≥ (log N)^{−C}/K, i.e. a (log N)^{−C}-small fraction. Consequences:
- This is a simplification of the proof of *top-level* relation supply, not a new exponent.
- It is a weaker theorem than Theorem U, which covers every N.
- It cannot be used in the recursive shape-test calls, whose inputs are not uniform (Prop HR).

**Lemma M1-AA (the counting premise of R2 holds for almost all N) [proved].** Call (f, g) a pair if:
- f ∈ Z[x, y] is a form of degree ≤ d;
- g = g₁x − g₀y with gcd(g₀, g₁) = 1;
- R := f(g₀, g₁) ≠ 0.

N *admits* the pair if N | R (a common root modulo N). Let ε > 0 and X be large.
- (a) The number of N ∈ [X, 2X) that admit some pair with d ≤ (log X)^{0.9} and
  (d + 1)log(2‖f‖ + 1) + 2 log(2‖g‖ + 1) < (½ − ε)log X is ≤ X^{1−ε+o(1)}.
- (b) For fixed β ∈ [0, 1) and d ≤ (log X)^β, the constant ½ improves to 1 − σ for any σ with (1 + β)(1 − σ) < 1, i.e.
  σ > β/(1 + β). This gives ¾ − ε at β = 1/3.
- (c) (codex) If N is B-rough (no prime factor ≤ B = e^{bE}) and d = o(log B), e.g. d ≍ D′, then the constant is
  1 − ε for every fixed ε > 0, hence 1 − o(1) by a diagonal choice of ε.

*Proof.*
- (a) Rankin with exponent ½: #{n ∈ [X, 2X) : n | R} ≤ Σ_{n|R}(2X/n)^{1/2} ≤ (2X)^{1/2}∏_{p|R}(1 − p^{−1/2})^{−1}.
  - The product over ω(R) primes is largest for the first ω(R) primes. Their primorial is ≤ |R|, so they are all
    ≤ y = O(L), where L := log|R|; and Σ_{p≤y}p^{−1/2} = O(y^{1/2}/log y).
  - For a violating pair, ‖f‖, ‖g‖ < X^{1/2} and |R| ≤ (d + 1)‖f‖‖g‖^d. So L ≤ (d + 2)log X ≤ (log X)^{1.9+o(1)}, and
    the window count is ≤ X^{1/2+o(1)}.
  - Bin the pairs by d and by dyadic heights H_f, H_g ≤ X^{1/2}: O((log X)^{2.9}) bins. Each bin has
    ≤ (2H_f + 1)^{d+1}(2H_g + 1)² ≤ X^{1/2−ε+o(1)} pairs; the dyadic rounding costs (d + 3)log 2 = o(log X).
  - Summing gives X^{1−ε+o(1)}.
- (b) The same with exponent σ: one needs y^{1−σ} = o(log X·log y) with y ≍ d log X.
- (c) N | R and N B-rough give N | R_{>B}, the B-rough part.
  - The number of its divisors is ≤ 2^{Ω(R_{>B})} ≤ 2^{L/log B}.
  - With L ≤ (d + 2)log X this is exp(O(d log X/log B)) = X^{o(1)} when d = o(log B).
  - At d = δD′ and log B = bE: d/log B ≍ D′/E = 1/log t → 0.
  - Repeat (a) with window count X^{o(1)}. ∎

*Consequences.*
- [nested] Prop R2's premise (M1), that usable pairs satisfy (d + 1)log H_f + 2 log H_g ≳ t, was a heuristic for
  worst-case N. Prop R2 noted that the rigorous constraint alone would let special N escape, e.g. f = x³ − 2. For
  almost all N, (M1) is a theorem: with constant ½ for every N (density X^{−ε+o(1)}), and 1 − o(1) for B-rough N.
- For the B-rough cofactor N′ = N/s of an arbitrary N (s the B-smooth part), apply (c) at scale X/s and truncate at
  s ≤ X^{1/2}.
  - Inputs with s > X^{1/2} have density ≤ exp(−(1 + o(1))u log u) with u = t/(2bE) ≍ D′²: super-polylogarithmically
    small, though not X^{−ε}.
  - The total exceptional density for preprocessed cofactors is X^{−ε+o(1)} + e^{−Θ(D′² log t)}, which is o(n^{−2}).
- R2's window argument uses (M1) only up to constants, so it goes through unchanged for almost all N. R2 stays a
  **model computation** (premises M2–M4). In that model, designs without algebraic children are still stuck at 2/5.
- SNFS-type pairs at degree d (‖f‖ = O(1), ‖g‖ ≤ CX^{1/d}) have |R| = O(d·C^d·X). So each is admitted by O(d·C^d) values
  of N, and such N have density ≤ X^{−1+2/d+o(1)} for d = o(log X) (proof critic). Averaging over N removes the
  structured inputs where polynomial selection could be lucky: it is the less favourable setting for sizes.
- `bounds.py` checks the Rankin window bound on toy inputs (tests).

**Lemma AV (average ⇒ most N) [proved, Markov].** If Pr_{N∼ν, A∼μ_N}[bad] ≤ δ, then ν{N : Pr_{μ_N}[bad] > ½} ≤ 2δ.
So "AD-avg(δ)", the average-case form of [nested] Hypothesis AD, implies AD for all N outside a set of density 2δ.
A lower bound on E_N p_N alone allows a positive fraction of N with p_N = 0 (codex).

## 3. Multi-level descent for almost all N

### 3.1 Size model and the needed supply statement

The size model of [rec] §3.4 and [nested] §8 uses typical sizes (log‖f‖ ≍ Q, values N^{τ_F}, τ_F = δb₁). Averaging over
N changes none of them, and by Lemma M1-AA no luckier pairs exist for almost all N. So the model still reaches L[1/3, c]
with O(log t) levels and e^{O(log²t)} nodes [model computation, as in nested].

The needed supply statement becomes AD-avg. By Prop NB this is AD° for the N-free mixture law μ̄: "a random pair (F, G)
drawn from μ̄ has, with probability ≥ 1 − δ, a good point in every special-q box". It is cosmetically nicer than AD,
because no conditioning on N | F(g₀, g₁) is needed. It meets the same method barriers:
- **component-level independence moment methods** meet Lemma OD (Cor OD-AA);
- **a uniform per-target union bound from approximate moments** meets Lemma MP.

**Lemma MP (moment precision) [proved, trivial].** Let X ≥ 0 be a random variable and η ∈ (0, 1). The mixture X′ := 0
with probability η, else X, using a coin independent of X, satisfies Pr[X′ = 0] ≥ η and E[X′^k] = (1 − η)E[X^k] for all
k ≥ 1. The same holds for conditional moments given any σ-algebra 𝒢 for which the coin is independent of σ(X, 𝒢).
Hence an argument whose only input is moment estimates accurate to relative error η (conditional on such 𝒢) cannot
certify Pr[X = 0] < η. (`bounds.mp_mixture`.)

*Consequence for AD [method barrier for one proof strategy, not impossibility].*
- A uniform per-target union bound over the e^{(b₁+o(1))Q} targets of a field needs Pr_A[τ dead] ≤ δe^{−b₁Q} per
  target, even for AD-avg(δ).
- If it is run from approximate moments of the good-point count alone, Lemma MP says those moments must be known to
  relative error ≈ e^{−b₁Q} = L[2/3]^{−Θ(1)}.
- Arguments that use additional structure (support constraints, exact identities) are not covered.
- *Heuristic order count.* With mean count e^{(2s−φ/3)E} ([nested] §8), a 2j-th central moment with j ≍ b₁D′/(2s − φ/3)
  is needed. That involves 2j-point configurations, against the budget d + 1 = δD′ + 1.
  - At nested's example parameters, 2s − φ/3 = 5/48, so j ≈ 9.6D′, far above 0.5D′.
  - Much larger boxes (s ≳ 1.7; codex gives s > 40/23 in one convention) bring it within budget at a higher L[1/3]
    constant.
  - So OD is not the only obstacle to this route; precision is.
- The best known relative errors for smooth numbers in progressions at L-scale have the shape
  exp(−cu/log²(u+1)) + y^{−c}, on average over moduli (Harper, arXiv 1208.5992, Thms 1–2; Fouvry–Tenenbaum; ε(N) in
  [U] §4).
  - For y = e^{λz} children (u ≍ D′) this is e^{−o(E)}.
  - For B-smoothness (u ≍ D′²) it is e^{−ω(E)}.
  - Either way it is far above e^{−b₁Q}, and the y^{−c} term alone is e^{−O(E)}.
- Conditioning on local factors does not change this (Lemma MP is conditional).
- None of this depends on N. It refines [nested] §6's union-bound paragraph.

### 3.2 Recursion and quantifiers

**Proposition HR (almost-all bounds do not compose; hereditary sufficiency).**
- (a) [proved, trivial; measure form] Let 𝒮 ⊂ [Y, 2Y) be any set of size ≤ εY. The statement "the time bound fails only
  on a set of density ≤ ε" is consistent with failure exactly on 𝒮.
  - As an algorithm: for a base algorithm 𝒜₀ meeting the bound everywhere and any 𝒮 decidable within the bound (e.g.
    from the computed factorisation), the variant that pads to FACTOR_U's time on 𝒮 meets Definition AA with density ε.
  - So an outer analysis that needs the recursive sub-inputs to lie in a set of density σ ≤ ε at scale Y, and is fed
    only "𝒜 has exceptional density ε", cannot conclude anything: the needed sub-inputs may be exactly the exceptions.
- (b) [proved] In NEST ([nested] §8) the needed recursive sub-inputs include two kinds, each of density e^{−Θ(E)} in
  its dyadic window:
  - accepted halves: |u₀|, |v₀| ∈ 𝒰_A, the products of k = ⌈t/(2b₁Q)⌉ + 1 ≍ D′ distinct primes of (Y, 2Y], of density
    ≤ (C/log N)^{k}·e^{O(k)} = e^{−Θ(E)};
  - algebraic values |F(P)| at good points: max(B, q^λ)-smooth numbers with u = O(D′) at every level, so of density
    e^{−Θ(E)}.

  So "o(1)", and even "e^{−cE}" for c below the relevant constants, is vacuous for the recursion.
- (c) [proved; codex's dovetail form; a sufficient condition] Suppose each recursive query at scale y is dominated by
  K(y) times the uniform measure, conditionally on history.
  - *Dovetailing* the uncapped fast algorithm with a universal fallback of expected cost U(y) gives expected cost
    O(T_f(y)) + O(K(y)ε(y)U(y)).
  - A single clocked attempt followed by the fallback would add Θ(θU) for the clock-failure probability θ (Markov),
    giving O(rT_f) + O((Kε + θ)U).
  - With U = L_y[2/5, c_U] at y = N^β, closing the recursion at L_N[1/3] this way needs
    log(1/ε(y)) ≥ c_UΛ_{2/5}(y) + log K + log M − O(E): a super-polynomially small exceptional density *and* a
    domination hypothesis on structured queries.
  - Neither is supplied by averaging over N.

**Hypothesis HAD(𝒦) (hereditary AD).** Fix the parameters (δ, λ, b, b₁, s) and field law 𝔉 of [nested] §8 and a set
𝒦 ⊆ Z_{≥2} ("tame inputs"). Put n_U(m) := exp((log m)^{4/5}). A point P of a target of a field for N is *𝒦-good* if:
- it is good in the sense of [nested] §8; and
- each of |G(P)| and |F(P)| is ≤ n_U(N) or lies in 𝒦.

AD_𝒦(A) means that every τ ∈ 𝒯(A) has a 𝒦-good point. HAD(𝒦) asserts the following.
- (H0) *Divisor closure.* For every m ∈ 𝒦, every divisor m′ | m with m′ > n_U(m) lies in 𝒦. This covers the B-rough
  part after preprocessing, perfect-power bases and the factor-tree nodes after a split or factor exit.
- For all large N ∈ 𝒦:
  - (H1) Pr_{A∼𝔉(N)}[factor exit, or A admissible and AD_𝒦(A)] ≥ ½;
  - (H2) for every admissible A, |𝒰_A ∩ 𝒦| ≥ |𝒰_A|/2, where 𝒰_A is the split set of [nested] §8's proof (products of
    k distinct primes of (Y, 2Y] not dividing D_A).

Its exceptional density is ε_𝒦(X) := |[X, 2X) ∖ 𝒦|/X. Relations to other hypotheses:
- Nested's worst-case AD implies HAD(Z_{≥2}), and (H0) and (H2) are then empty.
- For a proper 𝒦, HAD(𝒦) and AD are *incomparable*: 𝒦-good is pointwise stronger than good, while (H1) is required only
  on 𝒦.
- (H0)–(H2) concern thin sets, which by Prop HR(a) no density bound controls. Lemma AV gives only the 𝔉-average part of
  (H1), and only for all N outside density 2δ, from AD-avg(δ).

**Theorem NR-AA [proved conditional on HAD(𝒦), LV6′ ([nested] §8, with the Hensel prime and characters of item 3) and
Theorem U (internal chain); no ERH].** Run the algorithm NEST-AA[c′, K] in alternation with FACTOR_U, where K ≥ 1 is a
hard-coded constant. Its existence is shown in the proof; it is not effective. NEST-AA[c′, K] is NEST ([nested] §8),
inside the recurse wrapper (preprocessing, perfect powers, factor tree), with these changes:
1. Every recursive factorisation of an integer m > n_U(N) is a call of NEST-AA[c′, K] itself, alternated with trial
   division, with clock Kβ(m), β(m) := exp(c′Λ_{1/3}(m)).
   - It is repeated up to r := ⌈E²⌉ times with fresh coins. If all r calls time out, the point is treated as not
     good, or the row as rejected.
   - Integers m ≤ n_U(N) go to FACTOR_U, at cost L_m[2/5]^{O(1)} = exp(t^{8/25+o(1)}) = e^{o(E)}. This needs
     Theorem U; LP92 at this size would cost exp(t^{2/5+o(1)}).
2. A DESC call that returns FAIL in a field discards that field, as in NEST.
3. The character columns use the unconditional factorial-bit character primes of [U] §9 instead of Lemma X_ERH:
   s_χ = ⌈(2d)!(log N + log Δ̄ + d)⁴⌉ bits, Lemma E's ≤ 1 exceptional class, three square tests. The Hensel prime is
   one of the character primes, as in [nested].

Then the following hold.
- The algorithm is Las Vegas and completely factors every N in expected time ≤ 2·E T_U(N) = L_N[2/5, c_U + o(1)].
- For every c′ > c_NR-AA there is a K such that every N ∈ 𝒦 is factored in expected time ≤ Kβ(N) = L_N[1/3, c′],
  where
  c_NR-AA := max(2b, b + C_χδ, (b + 1/b₁)/(1 − 2^{−1/3}), (b + 2s)/(1 − τ_F^{1/3})), τ_F = δb₁ < 1.
  Here C_χ is an absolute constant covering the polynomial cost of s_χ-bit arithmetic: character sampling and
  evaluation, and Lemma SQ's Hensel lifting modulo an s_χ-bit prime.
- If ε_𝒦(X) = o((log X)^{−2}), then by Lemma TR the L[1/3] bound holds for all but o(1) of balanced semiprimes.

*Proof (changes to [nested]'s proof of Theorem NR).*
- **Correctness and termination.** Outputs are certified (AKS, gcds), and timeouts only reject candidates. The alternation gives the L[2/5] bound for every N.
- **Base case and the constant K.** Because every call alternates with trial division, E T(m) ≤ 2√m + O(1) for all m.
  - Let n₂(c′) be the threshold beyond which the recursion inequality below closes. Its o(1)s do not depend on K,
    because K enters linearly.
  - Choose K ≥ 8·max_{m<n₂}(2√m + O(1))/β(m). Then E T(m) ≤ Kβ(m)/4 for every m < n₂.
- **Strong induction over 𝒦.** Claim: E T(m) ≤ Kβ(m)/4 for all m ∈ 𝒦. Assume it below N ∈ 𝒦, N ≥ n₂.
  - By Markov a clocked call on m ∈ 𝒦 succeeds with probability ≥ ¾, so all r calls fail with probability ≤ 4^{−E²}.
  - A call on m ∉ 𝒦 costs at most rKβ(m) whatever happens.
  - *Factor tree (recurse §1 accounting).* Complete factorisation of N runs SPLIT on ≤ 2log₂N factor-tree nodes. By
    (H0) every node above n_U(N) lies in 𝒦, so HAD applies to each node's SPLIT. Nodes ≤ n_U(N) go to FACTOR_U at
    e^{o(E)}. Hence E T(N) ≤ 2log₂N·max_node E T_SPLIT(node) + e^{o(E)}.
  - The induction hypothesis is used only for SPLIT's own sub-calls, which are on integers ≤ √(2N) (halves) and
    ≤ N^{τ_F+o(1)} (values). It is never used for the factor-tree nodes, whose sizes need not shrink by a constant
    factor in the exponent.
- **Fields, per LV6 round.**
  - In a field with AD_𝒦(A), each target visited in one round has a first 𝒦-good point, whose two clocked calls both
    succeed except with probability ≤ 2·4^{−E²}.
  - One round visits ≤ e^{(b+o(1))E} targets (M_rows rows × e^{O(log²t)}-node trees, Lemma TS). So a round fails
    through timeouts with probability ≤ e^{(b+o(1))E}·4^{−E²} = o(1).
  - Over the expected ≤ 2 + o(1) rounds (Lemma T_multi) and, with (H1), the expected ≤ 2 + o(1) fields, the costs add
    up as in NR.
  - Work per row attempt is bounded by the clocks in good and bad fields alike. The factor r = E² = e^{o(E)} is
    absorbed.
- **Split.** Lemma S″ ([rec] §3.2) applied to S′ := 𝒰_A ∩ 𝒦 gives Pr ≥ |S′|²/(2φ(N)).
  - S′ has |S′| ≥ |𝒰_A|/2 ≥ C(|𝒫|, k)/2 ≥ C(|𝒫|, k)/(2k), so S″'s size hypothesis holds directly.
  - The condition |𝒫| ≥ 4k³2^ω holds as in NR, with |𝒫| = e^{Θ(Q)} and 2^ω ≤ e^{O(t/E)}.
  - The acceptance exponent is unchanged, e^{−(1/b₁+o(1))E}.
- **(T).** Timeouts and repetitions are functions of sub-inputs and fresh coins, which are independent of the
  hidden square roots. So Lemma T_multi holds and each LV6 round is fruitful with probability ≥ ½.
- **Characters.** [U] §9's argument is degree-uniform.
  - LO's unconditional effective Chebotarev needs log x ≥ 10n_L(log d_L)² = O(d³t²), and s_χ ≫ that.
  - Lemma E holds for odd d.
  - s_χ ≫ b₁Q, so the character primes exceed every atom prime.
  - At d = δD′, log s_χ = 2d log d + O(d + log t) = (2δ/3 + o(1))E. The convergence is slow: 0.17 / 0.20 / 0.24 of
    the limit 1/3 at δ = ½ for t = 10⁴ / 10⁶ / 10⁹ (proof critic).
  - So sampling, evaluation and Hensel lifting modulo s_χ-bit primes cost L[1/3, C_χδ]. This matches codex, [nested]
    §7's Θ(δ) correction, and LV Remark 6.7's 4δ/3 (whose δ normalisation we did not match to ours).
  - SQ's sign enumeration 2^{d−1} = e^{o(E)} is part of LV6′.
- **Recursion.** As in NR, with the clocked sub-calls, for each node m of N's factor tree (m ≤ N):
  T_SPLIT(m) ≤ K·[e^{(max(2b, b+C_χδ)+o(1))E} + e^{(b+1/b₁+o(1))E}·rβ(√(2N)) + e^{(b+2s+o(1))E}·rβ(N^{τ_F+o(1)})].
  - Since β(N^τ) = e^{c′τ^{1/3}E(1+o(1))} and c′ > c_NR-AA, the exponent max(2b, b + C_χδ, b + 1/b₁ + 2^{−1/3}c′,
    b + 2s + τ_F^{1/3}c′) is strictly below c′.
  - With the factor 2log₂N = e^{o(E)} from the factor tree, E T(N) ≤ K·e^{(c′−η+o(1))E} ≤ Kβ(N)/4 for N ≥ n₂(c′), for
    some η = η(c′) > 0. The o(1)s do not depend on K or N. ∎

*Remarks.*
- With 𝒦 = Z_{≥2}, NR-AA is nested's Theorem NR without ERH, at the cost of the additive C_χδ (and with Theorem U used
  for small sub-inputs).
- For a proper 𝒦, HAD is not "AD with a weaker quantifier": it adds the hereditary conditions (H0), (H2) and
  𝒦-goodness. Its substance, the same-field supply at L-scale, is untouched by averaging over N (Prop NB, Lemma MP,
  Cor OD-AA). So NR-AA is a correct bookkeeping of what an almost-all L[1/3] would need, not progress toward it.

### 3.3 LV6 at growing degree and ERH (task 2, last question)

- ERH is not needed for a finite constant at α = 1/3 (Theorem NR-AA, item 3). It would only remove the additive C_χδ
  and allow polynomial-bit characters.
- Averaging over N would need a zero-free region or an exceptional-zero density bound for the family {K_A(√h)} of
  fields of degree 2d, with d → ∞, on average over the field law.
  - The known family results (Pierce–Turnage-Butterbaugh–Wood, Invent. Math. 2020; Lemke Oliver–Thorner–Zaman, Invent. Math. 235 (2024))
    are for fixed degree and fixed Galois group.
  - We found nothing at growing degree.
  - In any case, Prop NB makes this a field average with an N-independent law.

## 4. Answers to the task, item by item

1.
   - (a) OD binds. N adds 0 directions beyond w (design a‴), or is the budget itself (plain base-m) (NB(iv), E1).
   - (b) SI's sparsity persists pointwise and on average (Cor SI-AA).
   - (c) AD becomes AD-avg = AD° under the N-free mixture law μ̄ (NB(i)). Theorem U's τ = 0 machinery gives
     single-box first and second moments only. An average over N acts as BDH-type spreading of N mod t, which replaces
     Lemma W (Remark W-AA) and nothing more. The union bound needs precision e^{−b₁Q} (Lemma MP), and component
     methods meet OD.
   - Genuine gains: Remark W-AA and Lemma M1-AA (rigour of R2's M1).
2.
   - The size model is unchanged.
   - Recursion needs hereditary control (Prop HR) and is closed under HAD (Theorem NR-AA).
   - (T) holds (budgets keep it).
   - LV6 at growing degree is unconditional with L[1/3, C_χδ] characters. ERH is not needed; averaging does not help.
3. Definition AA, Lemma TR (ε_n = o(n^{−2})), and the mean-versus-almost-all remark (§0).
4. Literature: §7.
5. **Best α established: 2/5** (Theorem U, every N, internal chain).
   - Conditional: 1/3 on 𝒦 under HAD + LV6′ + Theorem U.
   - Model: 1/3.
   - We know no tool by which averaging over N reaches an exponent in [1/3, 2/5). This is not an impossibility claim.

## 5. Experiment E1 (pre-registered `246ac18`; exact)

d ∈ {3, 5}; s ∈ {d − 1, d, d + 1, d + 2, 2d + 2}; K = 3(d + 3) + 10 samples per regime; 40-bit N.

| (d, s) | a_F(A) | a_F(B) | a_F(C) | a_F(D) | a_FG(D) | kernel dim | violations |
|---|---|---|---|---|---|---|---|
| (3, 2) | 2 | 2 | 2 | 2 | 4 | 0 | 0 |
| (3, 3) | 3 | 3 | 3 | 3 | 5 | 0 | 0 |
| (3, 4) | 3 | 4 | 4 | 4 | 6 | 0 | 0 |
| (3, 5) | 3 | 4 | 4 | 4 | 6 | 1 | 0 |
| (3, 8) | 3 | 4 | 4 | 4 | 6 | 4 | 0 |
| (5, 4) | 4 | 4 | 4 | 4 | 6 | 0 | 0 |
| (5, 5) | 5 | 5 | 5 | 5 | 7 | 0 | 0 |
| (5, 6) | 5 | 6 | 6 | 6 | 8 | 0 | 0 |
| (5, 7) | 5 | 6 | 6 | 6 | 8 | 1 | 0 |
| (5, 12) | 5 | 6 | 6 | 6 | 8 | 6 | 0 |

Regimes: (A) N, g, w fixed; (B) N, g fixed; (C) g fixed, N varies; (D) all vary.

R1–R4 all hold: a_F = min(s, d) in (A) and min(s, d + 1) in (B), (C) and (D); a_F(C) = a_F(B); a_FG = min(s, d + 1) +
min(s, 2); 0 kernel violations. This is an exact check of Prop NB(iii)–(iv) (code correctness, not new evidence beyond
the algebra).

## 6. Experiment E2 (pre-registered `246ac18`; toy)

Setup:
- the toy of [nested] E1: d = 3, 44-bit balanced semiprimes, B = 2^{10}, Y = 2^{14}, S = 32;
- FIX: 6 fixed N × 12 fields; RND: 72 N × 1 field;
- per field: 150 uniform algebraic and 150 uniform rational targets, and ≤ 80 parent nodes with matched baselines;
- 609 s total.

| prediction | result | verdict |
|---|---|---|
| Q1 \|mean_RND − mean_FIX\| ≤ 0.05, CI ⊂ [−0.08, 0.08] | 0.2852 vs 0.2705; diff +0.015, CI [−0.028, +0.057] | holds |
| Q2 SD_btw ratio ∈ [0.67, 1.5] | 0.1057 / 0.0930 = 1.14 | holds |
| Q3 ANOVA over the 6 fixed N, p > 0.01 | F = 1.30, permutation p = 0.276 | holds |
| Q4 J_pool > 1.2 and J_field ∈ [0.75, 1.33], both regimes | FIX 1.229 / 0.901; RND 1.301 / 0.912 (5509 / 5459 nodes) | holds |

Descriptive:
- Per-field algebraic rates range 0.09–0.55 (FIX) and 0.08–0.61 (RND).
- Rational rates 0.108 vs 0.122; SD ratio 1.35; ANOVA p = 0.049.
- The correlation of the algebraic rate with the root count of f modulo primes ≤ 50 is 0.37 (FIX) and 0.46 (RND): local
  factors carry much of the heterogeneity.

Interpretation (fixed in advance): consistent with Prop NB. Averaging over N acts only through the field law and leaves
the field-level heterogeneity intact (nested §6's Θ(1) relative variance). Q4 is J_field slightly below 1 in both
regimes, inside its band. Toy sizes support nothing rigorous.

## 7. Literature

Sources: codex almostall13_1 §4 (primary sources read where stated) and the Opus literature critic (§9), who read LP92,
Hafner–McCurley, Pomerance 1987, Lenstra 1987, LV, KW, Barbulescu–Jouve, PTBW, LPP I/II and Harper as primary texts.
- **No rigorous almost-all-N (or balanced-semiprime) complete factorisation below L[1/2] was found**, nor L[1/2, c < 1]
  on average over N. This is a search result, not a nonexistence theorem.
  - The worst-case record is Lenstra–Pomerance, JAMS 5 (1992): the unnumbered Theorem of §1, a consequence of
    Thm 10.5, L[1/2, 1 + o(1)].
  - LP92's remark after Thm 10.3 notes that Algorithm 10.1's time is *equal* to L[1/2, 1 + o(1)], and that changing y
    does not improve it.
- **The explicit rigorous almost-all-N result.** Hafner–McCurley, J. Algorithms 10 (1989), pp. 544–545: almost all
  n ≤ x are completely factored in expected time exp((log x)^θ) for θ > 5/6. They note that Pomerance's L[1/2, √2] is
  better.
- **ECM for almost all primes.**
  - Lenstra, Ann. Math. 126 (1987): Cor (2.8) is rigorous given the Hasse-interval smoothness probability; the
    L_p[1/2, √2] running time, Conj (2.10), rests on an unproved short-interval smoothness conjecture.
  - Pomerance 1987 (Discrete Algorithms and Complexity, 119–143), Thms 2.1 and B′ (via Friedlander–Lagarias): every
    prime ≤ v outside an exceptional set of size O(x·exp(−½(log x)^{1/6})) is found in time
    O((log N)⁴·exp(2(log v)^{6/7})). LP92 §7 (Thms 7.3–7.4) re-uses this as a smoothness *recogniser*.
  - Hafner–McCurley Thms 3.1–3.2 extend this to every θ > 5/6.
  - Barbulescu–Jouve (arXiv 2212.11724) give ECM smoothness asymptotics for a fixed curve over primes under
    Elliott–Halberstam-type hypotheses (Conj 1.1/4.1, Hypothesis 1).
  - **Worst case:** Lenstra–Pila–Pomerance, "A hyperelliptic smoothness test I", Phil. Trans. R. Soc. A 345 (1993),
    Thm 1.1/Cor 1.2, and II, Proc. LMS (3) 84 (2002): rigorous L_p[2/3, c₀]·(log n)² for every p.
  - All of these give at best L_N[2/3]-type or worse time for complete factorisation of typical N, because the
    second-largest prime factor is typically N^{Θ(1)}.
- **Averaging in rigorous NFS/DL.**
  - Lee–Venkatesan, JNT 187 (2018) (arXiv 1805.08873):
    - Thm 2.1: a possibly trivial congruence in L[1/3, (64/9)^{1/3}], for any n;
    - Thm 2.3: semiprimes with p ≡ q ≡ 3 mod 4, under Conj 7.1;
    - Remark 6.7: characters cost L(1/3, 4δ/3) unconditionally, polynomial under GRH. Their δ normalisation was not
      matched to our d = δD′ (the constant C_χ is uncomputed anyway).
    - The expectation is over the algorithm's own randomness, "for any n". They do not average over n.
  - Pomerance 1987: worst case L[1/2, √2]; its ECM step averages over primes via Thm B′.
  - Kleinjung–Wesolowski, JAMS 35 (2022), Thm 1.1: worst case over every F_{p^n}, (pn)^{2log₂n+O(1)}, from Frobenius
    structure absent for integers.
  - Neither we nor the critic found an input-averaged rigorous NFS or index-calculus theorem (e.g. for almost all p)
    for factoring or for prime-field DL.
- **Smooth numbers in progressions.** The exp(−cu/log²u) error shape is Harper (arXiv 1208.5992, Thms 1–2). It is a
  Bombieri–Vinogradov-type average over moduli, with error ≪ Ψ(x,y)²(e^{−cu/log²(u+1)} + y^{−c}) + Ψ(x,y)Q; the same
  shape appears in Fouvry–Tenenbaum. Either way it is e^{−o(E)} at our scale, far above e^{−b₁Q}.
- **Smooth values of random polynomials and forms** (Barbulescu–Lachand, Math. Comp. 2017; BBDT 2012; Martin; Dartyge):
  fixed degree, or y a power of x; no L-scale joint-correlation results ([nested] §10).
- **Family zero-free regions.** All known results are for fixed degree:
  - Pierce–Turnage-Butterbaugh–Wood, Invent. Math. 2020: the parameters depend on n;
  - Lemke Oliver–Thorner–Zaman, Invent. Math. 235 (2024) 893–971;
  - Chen An, arXiv 2004.14410.

  None is uniform as the degree → ∞. The Lagarias–Odlyzko label: [U] and lv6 cite "Theorem 1.3", PTBW "Corollary 1.3",
  and it is unverified in the printed volume. Only the content is used.

## 8. Candidate card

- **Mechanism.** Average over the input N, in addition to the field, to supply the free parameters that algebraic
  components lack ([nested] Lemma OD), for an almost-all L[1/3] (or α < 2/5).
- **Model.** Theorem NR / Theorem U architecture, design (a‴) and plain base-m, d ≍ D′.
- **Result.**
  - **Refuted for this model:** the claim that N adds ≈ d + 1 parameters; it adds none beyond w (Prop NB; E1, E2).
  - Proved: Prop NB, Cors OD/SI/Info-AA, Remark W-AA, Lemma M1-AA, Lemma AV, Lemma MP, Prop HR.
  - Proved conditional on HAD + LV6′: Theorem NR-AA (no ERH).
- **Failure scope.** The parameter-supply mechanism, for fixed-polynomial special-q designs (the [nested] §8 family:
  multi-field, constrained resampling, plain base-m, design a‴) in which N | F(g₀, g₁) and the descent inside a field
  depends only on the field (Lemma DF, contracting regime). Designs doing arithmetic mod N inside the descent (escapeR2
  re-split, Prop RS) are outside this scope and have their own obstructions.
- **Not excluded.**
  - Benefits of a different field law that are not parameter counts (this note proves no improvement, but no
    impossibility either).
  - A genuine L-scale joint-smoothness theorem for random polynomial pairs. That would be new analysis, independent of
    N.
- **Next step.** The bottleneck is a precision problem (Lemma MP) or a same-field correlation problem (OD) for a random
  pair (F, G), and N-free. Averaging over N is not the missing ingredient.

## 9. Reviews

All valid findings were applied in the text above.

**GPT-6.1-sol blind round `codex/almostall13_1`** (xhigh, before my draft).
- Same verdict: 2/5 is the supported exponent, and averaging over N adds no coefficient rank beyond (w, c).
- Contributions: the rank table; the Lemma W substitute; the counting form of M1; the B-rough divisor bound (used as
  Lemma M1-AA(c)); the dovetail fallback analysis (Prop HR(c)); unconditional characters at L[1/3, O(δ)]; the
  quantified Definition AA; Hafner–McCurley.

**GPT-6.1-sol adversarial referee `codex/almostall13_2`** (xhigh). Recommendation: MAJOR REVISION. Fixes applied:
1. *FATAL.* NR-AA lacked heredity for the factor-tree inputs (preprocessing cofactor, perfect-power bases, divisors after
   a split). Added (H0), divisor closure.
2. *MAJOR.* Lemma DF must exclude gcd = N. It is now restricted to the contracting regime 0 < |value| < N. The
   target-set inclusion was reversed.
3. *MAJOR.* The scope claimed "all designs". It is now restricted to fixed-polynomial special-q designs, with escapeR2's
   re-split acknowledged; and μ̄ is not uniform.
4. *MAJOR.* NB(iv) needs nondegeneracy assumptions, now stated. The "rd + 1" multi-field bound is false with
   independent w_j; it is now restricted to fixed w_j.
5. OD-AA's off-by-one: d − k during a run, d + 1 − k across runs. Info-AA needs "given g and w".
6. W-AA: the displayed inequality, X-dependent parameters, Markov on the per-N second moment.
7. M1-AA: the β range, the d-dependent SNFS density, and B-rough inputs versus cofactors.
8. *MAJOR.* MP's scope narrowed: the coin must be independent of 𝒢; the barrier is for uniform per-target union bounds
   from approximate moments. The moment order j ≈ 9.6D′ at nested's parameters exceeds d + 1.
9. *MAJOR.* HR(c): a clocked fallback adds a θU term; it is now dovetailed.
10. HR(a) put in measure form; HR(b)'s blanket G-value claim dropped.
11. Bound per LV6 round.
12. *MAJOR.* The headline is narrowed to "the best bound established here; no improvement proved; not an impossibility".
13. Literature: LP92 7.2 is an algorithm; Barbulescu–Jouve's published version is 2024 (venue not checked by us).
14. Q3 wording: "no statistically detected effect".
- It reproduced E1 and all E2 summary values exactly from the raw files.

**Opus proof critic (read-only, independent).** No FATAL issue; three MAJOR issues, all in NR-AA. Fixes applied:
1. The repetition count r = ⌈E⌉ was too small for b > ln 4. It is now r = ⌈E²⌉ (I had changed this before the report
   arrived).
2. The base case for sub-inputs below the ineffective threshold: added the hard-coded constant K and alternation with
   trial division.
3. 𝒦 was not closed under the integers NEST actually runs on: (H0).
4. The dependency on Theorem U is now in NR-AA's label.
5. C_χ now includes SQ lifting modulo the s_χ-bit Hensel prime.
6. W-AA: rough N′ gives C_W = O(log B); the exceptional fraction is only (log N)^{−C}; top-level relation supply only.
7. M1-AA consequences: truncation for preprocessed cofactors, and unified B.
8. A better SNFS density, X^{−1+2/d+o(1)}.
9. DF: "for every coin sequence".
10. Labels: OD-AA's children count is a model; HAD and AD are incomparable; the Verdict's overclaim; TR "only if" →
    "if"; HR(c) is sufficient, not necessary.
11. HR(a): stated relative to a base algorithm.
12. MP: 2j-point configurations; u ≍ D′² for B-smoothness.
13. The d − k inconsistency, and "rd + 1".

I also found myself that the factor-tree bookkeeping in NR-AA's induction must go through SPLIT costs, as in recurse
§1, not through E T of the divisors. This is now fixed.

**Opus literature critic (read-only, independent; primary sources).** Fixes applied:
- The 6/7 bound and the exceptional set are re-attributed to Pomerance 1987 (Thms 2.1, B′), re-used in LP92
  Thms 7.3–7.4.
- LP92's main theorem is the §1 Theorem via Thm 10.5, with the remark after 10.3.
- Hafner–McCurley pp. 544–545 is cited as the explicit almost-all-N complete-factorisation result (θ > 5/6).
- Added LPP I/II (rigorous L_p[2/3] for every p).
- KW is worst-case, not "average over input".
- Barbulescu–Jouve moved to ECM, with its hypotheses stated.
- Venue added for LOTZ, plus Chen An's fixed-degree log-free zero-density result.
- Harper's BV-type error shape includes the y^{−c} term.
- The Lagarias–Odlyzko theorem label is unverified (content only).

## 10. Reproduction

```
.venv/bin/python -m src.almostall13.rank experiments/raw/almostall13_e1.json
.venv/bin/python -m src.almostall13.harness_e2 experiments/configs/almostall13_e2.json
.venv/bin/python -m src.almostall13.analyze_e2 experiments/raw/almostall13_e2.jsonl experiments/raw/almostall13_e2_summary.json
.venv/bin/python -m pytest -q tests/test_almostall13.py
```
