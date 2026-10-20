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

**Plain language.** Averaging over the input N does not supply the missing parameters, so it does not get below
L[2/5]. The best exponent we can support for almost all N is **2/5**, from Theorem U, which already covers every N
(conditional on the internal chain). Nothing below 2/5 comes from this route, rigorously or under ERH. The reason is
short:
- Once the field is fixed, N is fixed too (N = F(g₀, g₁)/w).
- The descent inside a field never looks at N.
- So averaging over N is just averaging over fields under a different law on the same (d + 3)-dimensional space of
  polynomial pairs.
- N's digits move only one coordinate of F, its value T = F(g₀, g₁) at the common root, and the multiplier w already
  randomises that coordinate. In the plain base-m design N's digits are a *substitute* for the random coefficients c,
  not an addition to them.
- Nested's per-field budget of d + 1 parameters is a property of a degree-d form, whatever its law.

Precisely:
1. **Prop NB (no budget gain) [proved].** For any design with N | F(g₀, g₁), the following hold.
   - Every descent statistic is a function of the field alone, up to factor exits (Lemma DF).
   - An average over (N, field) equals a field average under the mixture law μ̄ (pushforward).
   - At s fixed points, the values have rank ≤ min(s, d + 1) under every joint law, with law-independent integer
     relations when s > d + 1.
   - In design (a‴), varying N adds exactly 0 affine dimensions beyond w, and 1 beyond fixed (N, w).
   - Lemma OD, Prop SI, the information version and the constrained-resampling budget hold verbatim (Cors OD-AA,
     SI-AA). E1 confirms the dimension counts exactly.
2. **What averaging over N does buy [proved]. None of it moves α.**
   - (a) Lemma W becomes unnecessary: N itself spreads the forced residues, Pr ≤ 1/t + 1/X (Remark W-AA).
   - (b) **Lemma M1-AA.** For all N outside a set of density X^{−ε}, every polynomial pair with a common root mod N
     obeys the counting premise (M1) of [nested] Prop R2. The constant is ½ for all N and 1 − o(1) for B-rough N.
     Here d ≤ (log X)^{0.9}. So R2's model floor of 2/5 for designs without algebraic children applies to almost all
     N with M1 *proved*. The SNFS-type escapes exist only on density X^{−1/2+o(1)}. If anything, the almost-all
     setting is the less favourable one for the size model.
   - (c) An average-case AD implies AD for most N, by Markov (Lemma AV).
3. **Recursion does not tolerate an "o(1)" exceptional set (Prop HR) [proved].** The recursive sub-inputs live in thin
   sets of density L^{−Θ(1)}: the accepted halves u₀, v₀, which are products of k primes in a window, and the good
   smooth values F(P), G(P). An exceptional set of density o(1), or even of density L^{−C} for small C, can contain
   all of them. A recursive almost-all design therefore needs a *hereditary* hypothesis. With one, **Theorem NR-AA**
   holds [proved conditional on Hypothesis HAD + LV6′; no ERH]: a Las Vegas algorithm, alternating with FACTOR_U,
   completely factors every N.
   - Its expected time is L[2/5] for all N.
   - For N in the tame set 𝒦 it is L_N[1/3, c] for every c > c_NR-AA.
   - The transfer to semiprimes needs ε_𝒦 = o(n^{−2}).
   - HAD follows from nested's worst-case AD, so NR-AA also removes ERH from nested's Theorem NR.
4. **ERH and characters [proved, routine].** ERH is not needed. Unconditional factorial-bit character primes ([U] §9)
   cost L[1/3, C_χδ] at d = δD′, which only raises the constant. Averaging over N gives no growing-degree family
   zero-free region (none is known).
5. **Lemma MP and the union-bound route [proved; method barrier, not impossibility].**
   - Moment information accurate to relative error η cannot certify Pr[dead] < η.
   - So proving AD, or average AD, by a union bound over e^{b₁Q} targets needs smooth-correlation asymptotics with
     relative error e^{−b₁Q} = L[2/3]^{−Θ(1)}.
   - Known L-scale errors are exp(−cu/log²u) with u ≍ D′, which is far too weak. This does not depend on N.
6. **Literature (§7).**
   - No rigorous almost-all-N factoring below L[1/2] was found. The rigorous almost-all-prime ECM results are much
     weaker than L_p[1/2]: Hafner–McCurley 1989 give exp((log v)^θ) for θ > 5/6, and LP92 §7 gives
     exp((log p)^{6/7}).
   - Lee–Venkatesan average over the algorithm's coins, not over N.
   - The "average over input" rigorous DL results (Kleinjung–Wesolowski) cover every small-characteristic field and use
     Frobenius structure that integers lack.
7. **Experiments.**
   - E1: all four predictions hold exactly (affine dimension d with (N, w) fixed; d + 1 otherwise, so N adds nothing
     beyond w; 0 kernel violations).
   - E2 (toy, d = 3, 44-bit N, 72 + 72 fields): all four predictions hold.
     - Q1: random-N versus fixed-N mean algebraic descendability 0.285 vs 0.270, CI [−0.028, 0.057].
     - Q2: between-field SD ratio 1.14. N-averaging does *not* reduce field heterogeneity.
     - Q3: no N effect within FIX (permutation p = 0.28).
     - Q4: nested's Jensen gap persists, J_pool = 1.23 / 1.30, and per-field rates explain it, J_field = 0.90 / 0.91.
8. GPT-6.1-sol's blind round (`codex/almostall13_1`) reached the same verdict independently: "the supported exponent
   remains 2/5 … input averaging supplies residue averaging and the counting bound … no extra coefficient rank beyond
   w, c". It contributed the B-rough form of M1-AA and the dovetail form of Prop HR.

---

## 0. Setting, notation, quantifiers

t := log N, E := t^{1/3}(log t)^{2/3}, Q := t^{2/3}(log t)^{1/3}, D′ := (t/log t)^{1/3}; L_N[1/3, c] = e^{(c+o(1))E}.

An **NFS-type design** samples, on input N, a field datum A = (F, G, w). Here:
- F ∈ Z[x, y]_d is a binary form;
- G = g₁x − g₀y with gcd(g₀, g₁) = 1;
- F(g₀, g₁) = N·w. This is design (a‴) of [U] and [nested] §8; w = 1 gives design (a″_k).

More generally only N | F(g₀, g₁) ≠ 0 is used. Public parameters (d, B, box radii, …) depend only on the bit length of
N. Descent (DESC), targets 𝒯(A), boxes and "good points" are as in [nested] §8.

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

**Lemma TR (transfer) [proved; polyimpl §7 D1, detqs §1].** Let S_n be the n-bit balanced semiprimes, |S_n| ≍ 2^n/n².
If the bad set in [2^{n−1}, 2^n) has size ≤ ε_n2^{n−1}, then the bad fraction of S_n is O(ε_n n²). So ε_n = o(n^{−2})
gives "all but o(1) of balanced semiprimes"; density merely o(1) can contain every balanced semiprime. An almost-all
bound is not a mean bound: with an L[2/5, c_U] fallback, the mean over [X, 2X) is controlled only if
ε(X) ≤ L_X[α, c]/L_X[2/5, c_U].

Theorem U satisfies Definition AA with E_X = ∅ at α = 2/5. That is our best almost-all exponent.

## 1. Averaging over N is averaging over fields (Prop NB)

**Lemma DF (descent depends only on the field) [proved].** Let DESC°(A, τ) be DESC with D_A° := g₁·f_d·disc(f)
instead of D_A = N·g₁·f_d·disc(f), and without gcd checks against N. For every N with N | F(g₀, g₁), DESC(N; A, τ)
either makes a factor exit or returns exactly DESC°(A, τ): the same tree, the same success indicator. Likewise
𝒯(A) ∖ 𝒯°(A) consists only of primes dividing N.
*Proof.* Inspection of [nested] §8. The box comes from the Gauss-reduced basis of the target lattice, which is fixed by
(q, G) or by (q, r) with F(r, 1) ≡ 0 mod q. The values are G(P) and F(P). The canonical factorisations are functions of
the values. The good test compares prime factors with B and e^{λz} and with D_A. N enters only through gcd checks,
where any nontrivial gcd is a factor exit, and through the primes p | N in D_A. If such a p divides a value, the gcd
check has already exited. Rational children sent to fresh fields ([nested] Prop R) are covered by applying this to each
field. ∎

**Proposition NB (no budget gain) [proved].** Let ν be any law on N (e.g. uniform on [X, 2X)), μ_N the design's field
law, and μ̄ := Σ_N ν(N)μ_N.
- (i) *Pushforward.* For every event ℰ that is a function of DESC° runs in one field,
  Pr_{N∼ν, A∼μ_N}[ℰ(A)] = Pr_{A∼μ̄}[ℰ(A)], up to factor exits.
- (ii) *N is a function of the field.* Given A, N = F(g₀, g₁)/w. If only N | F(g₀, g₁) is used, N is one of at most
  τ(F(g₀, g₁)) divisors. DESC° does not use N at all.
- (iii) *The budget is intrinsic.* For points P₁, …, P_s and any joint law of (N, A), the vector (F(P_i))_i equals
  E_P·(f₀, …, f_d)^T, with E_P[i][j] = a_i^j b_i^{d−j}. So it lies a.s. in a Z-module of rank ≤ min(s, d + 1) that
  does not depend on the law. For s > d + 1, the left kernel of E_P (integer relations λ depending only on the points)
  annihilates it a.s. Jointly with (G(P_i))_i the rank is ≤ min(s, d + 1) + min(s, 2).
- (iv) *Exactly what N adds.* Fix g. Then {F : F(g₀, g₁) = T} = F̂_{g,T} + G·Z[x, y]_{d−1}. In design (a‴),
  F = F̂_{g,Nw} + G·R_c, so (N, w, c) ↦ F factors through (T, c) with T = Nw, and N changes only the coordinate T. The
  affine dimension of the support of (F(P_i))_{i≤s} is:
  - min(s, d) with (N, g, w) fixed;
  - min(s, d + 1) as soon as w or N varies.

  So N adds 0 dimensions beyond w, and 1 beyond fixed (N, w). In the plain base-m design (no c, no w), the d + 1 digits
  of N *are* the coefficients: averaging over N supplies at most the d + 1 directions that c supplies in (a‴).

*Proof.*
- (i): Lemma DF and the definition of a mixture.
- (ii): trivial.
- (iii): F(P_i) = Σ_j f_j a_i^j b_i^{d−j}, so the rank is ≤ rank E_P; the G-part has 2 coefficients.
- (iv): if F(g₀, g₁) = F̂(g₀, g₁) = T, then F − F̂ vanishes at the primitive point (g₀ : g₁), so G divides it in
  Q[x, y] and, by Gauss's lemma, in Z[x, y]. In `make_f`, F̂ = gbe(Nw, …) depends on (N, w) only through Nw. The affine
  dimensions follow from rank min(s, d) of the c-part (R_c has d coefficients, multiplied by the nonzero G(P_i)) plus
  one transverse direction T. ∎

**Corollary OD-AA [proved].** [nested] Lemma OD holds verbatim for the joint (N, field) law. With s > d + 1 fixed points
the values satisfy fixed linear relations. The per-component budget of [nested] §4 (constrained resampling
F + G·∏L_{P_i}·H keeps d + 1 − k directions) is also unchanged during a run, because N is fixed by the field. Releasing
the evaluation constraint across runs gives at most one more direction (codex: d + 1 − k instead of d − k). The
typical-factorisation count ≈ d/(2λ) of children is unchanged, so depth-1 overdetermination for λ < 1/2 persists.

**Corollary SI-AA [proved].** [nested] Prop SI(i) is a fixed-field statement and holds pointwise in (N, A): the family
𝒮_{N,A} := 𝓘 ∖ Φ_{N,A}(Ω) has density 1 − o(1) and is hit with probability 0 for every (N, A), so also on average.
SI(ii), full support for *field-independent* families, is the only averaging route, and averaging over N is a special
case of it (NB(i)). SI(iii) is unchanged: the good family of an algebraic child is "descendable in A", a function of A.

**Corollary Info-AA [proved].** Given g, after d + 1 revealed values at points with independent evaluation rows, F is
determined, and with it N = F(g₀, g₁)/w. So beyond depth 1 every needed statement is pointwise in the field, exactly as
in [nested] §5's information version.

*Answer to task 1(a).* The idea "N's digits add ≈ d + 1 parameters" is false in design (a‴): c already randomises the d
directions of the fibre {F(g₀, g₁) = T}, w randomises T, and N moves only T. It is true in the plain base-m design, but
there the digits are the whole budget. The joint coefficient space of r fields that share g and N has dimension
≤ rd + 1, not r(d + 1) (codex).

## 2. What averaging over N does buy

**Remark W-AA (Lemma W is replaced by N) [proved; the substitution in LD*-U Step 4 is line by line].** Fix g, P = (a, b)
and t := |G(P)| with gcd(b·g₁·w, t) = 1. Since a ≡ g₀bg₁^{−1} (mod t) and F̂ is homogeneous of degree d,
F̂(P) ≡ (bg₁^{−1})^d·F̂(g₀, g₁) = (bg₁^{−1})^d·Nw (mod t). For N uniform on X consecutive integers,
Pr_N[F̂(P) ≡ r (mod t)] ≤ 1/t + 1/X.

This is Lemma W-U's conclusion, with C_W ≈ 1 instead of 24/η₀, and with no Brun–Titchmarsh and no N^{η₀}-margin for w.
So in design (a″_k) (w = 1) the second moment LD*-U holds on average over N, and by Markov (Lemma AV) Theorem U's
relation supply holds for all but a small fraction of N without the w-device.
- This is a simplification of the proof, not a new exponent.
- It is a *weaker* theorem than Theorem U, which covers every N.
- Caveat (codex): field acceptance and earlier relations can bias N mod t, so conditional uses need justification. In
  LD*-U the average is unconditional, so that is fine there.

**Lemma M1-AA (the counting premise of R2 holds for almost all N) [proved].** Call (f, g) a pair if:
- f ∈ Z[x, y] is a form of degree ≤ d;
- g = g₁x − g₀y with gcd(g₀, g₁) = 1;
- R := f(g₀, g₁) ≠ 0.

N *admits* the pair if N | R (a common root modulo N). Let ε > 0 and X be large.
- (a) The number of N ∈ [X, 2X) that admit some pair with d ≤ (log X)^{0.9} and
  (d + 1)log(2‖f‖ + 1) + 2 log(2‖g‖ + 1) < (½ − ε)log X is ≤ X^{1−ε+o(1)}.
- (b) For d ≤ (log X)^β the constant ½ improves to 1 − σ for any σ > β/(1 + β) (¾ − ε at β = 1/3, 1 − ε at
  β = o(1)).
- (c) (codex) If N is B-rough (no prime factor ≤ B = e^{bΛ}, as after preprocessing) and d = o(log B) (in particular d ≍ D′ with
  B = e^{bE}), then the constant is 1 − o(1).

*Proof.*
- (a) Rankin with exponent ½: #{n ∈ [X, 2X) : n | R} ≤ Σ_{n|R}(2X/n)^{1/2} ≤ (2X)^{1/2}∏_{p|R}(1 − p^{−1/2})^{−1}.
  - The product over ω(R) primes is largest for the first ω(R) primes. They are all ≤ y ≍ L := log|R|, and
    Σ_{p≤y}p^{−1/2} = O(y^{1/2}/log y).
  - For a violating pair, ‖f‖, ‖g‖ < X^{1/2} and |R| ≤ (d + 1)‖f‖‖g‖^d, so L ≤ (d + 2)log X ≤ (log X)^{1.9+o(1)}, and
    the window count is ≤ X^{1/2+o(1)}.
  - Bin the pairs by d and by dyadic heights H_f, H_g ≤ X^{1/2}: O((log X)^{2.9}) bins. Each bin has
    ≤ (2H_f + 1)^{d+1}(2H_g + 1)² ≤ X^{1/2−ε+o(1)} pairs (the dyadic rounding costs (d + 3)log 2 = o(log X)).
  - Summing gives X^{1−ε+o(1)}.
- (b) The same with exponent σ: one needs y^{1−σ} = o(log X·log y) with y ≍ d log X.
- (c) N | R and N B-rough give N | R_{>B}, the B-rough part. The number of its divisors is
  ≤ 2^{Ω(R_{>B})} ≤ 2^{L/log B}. With L ≤ (d + 2)log X this is exp(O(d log X/log B)) = X^{o(1)} when d = o(log B). At
  d = δD′ and B = e^{bE}: d/log B ≍ D′/E = 1/log t → 0. Repeat (a) with window count X^{o(1)}. ∎

*Consequences.*
- [nested] Prop R2's premise (M1), that usable pairs satisfy (d + 1)log H_f + 2 log H_g ≳ t, was a heuristic for
  worst-case N. Prop R2 noted that the rigorous constraint alone would let special N escape, e.g. f = x³ − 2. For
  almost all N, (M1) is a theorem: with constant 1 − o(1) after preprocessing, and ½ unconditionally.
- R2's window argument uses (M1) only up to constants, so it goes through unchanged for almost all N. In the size model
  (premises M2–M4 unchanged), designs without algebraic children are still stuck at 2/5.
- SNFS-type pairs (‖f‖ = O(1), ‖g‖ ≈ N^{1/d}) exist only on a set of density ≤ X^{−1/2+2/d+o(1)}. Averaging over N
  removes the structured inputs where polynomial selection could be lucky: it is the less favourable setting for sizes.
- `bounds.py` checks the window bound exactly on toy inputs (tests).

**Lemma AV (average ⇒ most N) [proved, Markov].** If Pr_{N∼ν, A∼μ_N}[bad] ≤ δ, then ν{N : Pr_{μ_N}[bad] > ½} ≤ 2δ.
So "AD-avg(δ)", the average-case form of [nested] Hypothesis AD, implies AD for all N outside a set of density 2δ.
A lower bound on E_N p_N alone allows a positive fraction of N with p_N = 0 (codex).

## 3. Multi-level descent for almost all N

### 3.1 Size model and the needed supply statement

The size model of [rec] §3.4 and [nested] §8 uses typical sizes (log‖f‖ ≍ Q, values N^{τ_F}, τ_F = δb₁). Averaging over
N changes none of them, and by Lemma M1-AA no luckier pairs exist for almost all N. So the model still reaches L[1/3, c]
with O(log t) levels and e^{O(log²t)} nodes [model computation, as in nested]. The needed supply statement becomes
AD-avg. By Prop NB this is AD° for the N-free mixture law μ̄: "a random polynomial pair (F, G) of the given heights has,
with probability ≥ 1 − δ, a good point in every special-q box". It is cosmetically nicer than AD, because no
conditioning on N | F(g₀, g₁) is needed. It has exactly the same method barriers:
- **component-level moment methods** meet Lemma OD (Cor OD-AA);
- **the union bound over single boxes** meets Lemma MP.

**Lemma MP (moment precision) [proved, trivial].** Let X ≥ 0 be a random variable and η ∈ (0, 1). The mixture X′ := 0
with probability η, else X (independently), satisfies Pr[X′ = 0] ≥ η and E[X′^k] = (1 − η)E[X^k] for all k ≥ 1. The
same holds for conditional moments given any σ-algebra (local factors, history). Hence no argument that uses only moment
estimates accurate to relative error η can certify Pr[X = 0] < η. (`bounds.mp_mixture`.)

*Consequence for AD [method barrier].*
- A union bound over the e^{(b₁+o(1))Q} targets of a field needs Pr_A[τ dead] ≤ δe^{−b₁Q} per target, even for
  AD-avg(δ). By Lemma MP, the moments of the good-point count must then be known to relative error e^{−b₁Q} =
  L[2/3]^{−Θ(1)}.
- The moment order needed is j ≍ b₁D′/(2s − φ/3), using the heuristic mean e^{(2s−φ/3)E} of [nested] §8. This is
  within the budget d + 1 = δD′ + 1 if the boxes are large, so OD does *not* block this route.
- Precision does. The best relative errors known for smooth numbers in progressions at this scale are of the form
  exp(−cu/log²u) with u ≍ D′ (e.g. ε(N) in [U] §4, from Drappeau/Harper-type estimates). That is e^{−o(E)}, far above
  e^{−b₁Q}.
- Conditioning on local factors does not help: Lemma MP is conditional.
- None of this depends on N. This refines [nested] §6's union-bound paragraph. It is a barrier for moment methods, not
  an impossibility result.

### 3.2 Recursion and quantifiers

**Proposition HR (almost-all bounds do not compose; hereditary sufficiency).**
- (a) [proved, trivial] Let 𝒮 ⊂ [Y, 2Y) be any set of size ≤ εY. There is an algorithm 𝒜_𝒮 (fast off 𝒮, running
  FACTOR_U on 𝒮) that meets Definition AA with exceptional density ε and is slow on all of 𝒮. Suppose the outer
  analysis needs the recursive sub-inputs to lie in a set of density σ ≤ ε at scale Y, and is fed only "𝒜 has
  exceptional density ε". Then it cannot conclude anything: 𝒮 may be exactly the set of needed sub-inputs.
- (b) [proved] In NEST ([nested] §8) the needed sub-inputs are of two kinds, both of density L^{−Θ(1)}:
  - accepted halves: |u₀|, |v₀| ∈ 𝒰_A, the products of k ≍ D′ distinct primes of (Y, 2Y], of density
    ≈ (2/log N)^k = e^{−Θ(E)} in their window;
  - good values |F(P)|, |G(P)|: smooth numbers of density ρ = e^{−Θ(E)}.

  So "o(1)", and even "L^{−C}" for C below the relevant densities, is vacuous for the recursion.
- (c) [proved, codex's dovetail form] Suppose instead that each recursive query at scale y is dominated by K(y) times
  the uniform measure, conditionally on history. Then a clocked fast attempt plus a universal fallback of expected cost
  U(y) costs O(T_f(y)) + O(K(y)ε(y)U(y)). With U = L_y[2/5, c_U] at y = N^β, closing the recursion at L_N[1/3] needs
  log(1/ε(y)) ≥ c_UΛ_{2/5}(y) + log K + log M − O(E): a super-polynomially small exceptional density *and* a domination
  hypothesis on structured queries. Neither is supplied by averaging over N.

**Hypothesis HAD(𝒦) (hereditary AD).** Fix the parameters (δ, λ, b, b₁, s) and field law 𝔉 of [nested] §8 and a set
𝒦 ⊆ Z_{≥2} ("tame inputs"). Put n_U(N) := exp(t^{4/5}). A point P of a target is *𝒦-good* if:
- it is good in the sense of [nested] §8; and
- each of |G(P)| and |F(P)| is ≤ n_U(N) or lies in 𝒦.

AD_𝒦(A) means that every τ ∈ 𝒯(A) has a 𝒦-good point. HAD(𝒦) asserts, for all large N ∈ 𝒦:
- (H1) Pr_{A∼𝔉(N)}[factor exit, or A admissible and AD_𝒦(A)] ≥ ½;
- (H2) for every admissible A, |𝒰_A ∩ 𝒦| ≥ |𝒰_A|/2, where 𝒰_A is the split set of [nested] §8's proof (products of
  k distinct primes of (Y, 2Y] not dividing D_A).

Its exceptional density is ε_𝒦(X) := |[X, 2X) ∖ 𝒦|/X. Nested's worst-case AD implies HAD(Z_{≥2}), and Lemma AV gives
(H1) for all N outside density 2δ from AD-avg(δ). But (H1)'s "𝒦-good" and (H2) concern thin sets, which Prop HR(a)
says no density bound controls.

**Theorem NR-AA [proved conditional on HAD(𝒦) and LV6′ ([nested] §8); no ERH].** Run the algorithm NEST-AA[c′] in
alternation with FACTOR_U. NEST-AA[c′] is NEST ([nested] §8) with these changes:
1. Every recursive factorisation of an integer m > n_U(N) is a call of NEST-AA[c′] itself with clock
   β(m) := exp(c′Λ_{1/3}(m)), repeated up to r := ⌈E⌉ times with fresh coins. If all r calls time out, the point is
   treated as not good, or the row as rejected. Integers m ≤ n_U(N) go to FACTOR_U, at cost
   L_m[2/5]^{O(1)} = exp(t^{8/25+o(1)}) = e^{o(E)}.
2. A DESC call that returns FAIL in a field discards that field, as in NEST.
3. The character columns use the unconditional factorial-bit character primes of [U] §9 (s_χ = ⌈(2d)!(log N + log Δ̄ +
   d)⁴⌉ bits, Lemma E's ≤ 1 exceptional class, three square tests) instead of Lemma X_ERH. The Hensel prime is one of
   the character primes, as in [nested].

Then the following hold.
- The algorithm is Las Vegas and completely factors every N in expected time ≤ 2·E T_U(N) = L_N[2/5, c_U + o(1)].
- For every c′ > c_NR-AA, every N ∈ 𝒦 is factored in expected time L_N[1/3, c′], where
  c_NR-AA := max(2b, b + C_χδ, (b + 1/b₁)/(1 − 2^{−1/3}), (b + 2s)/(1 − τ_F^{1/3})), τ_F = δb₁ < 1,
  and C_χ is an absolute constant from the polynomial cost of s_χ-bit arithmetic.
- If ε_𝒦(X) = o((log X)^{−2}), then by Lemma TR the L[1/3] bound holds for all but o(1) of balanced semiprimes.

*Proof (changes to [nested]'s proof of Theorem NR).*
- **Correctness and termination.** Outputs are certified (AKS, gcds), and timeouts only reject candidates. The
  alternation gives the L[2/5] bound for every N.
- **Strong induction over 𝒦.** Claim: E T(m) ≤ β(m)/4 for all m ∈ 𝒦 with m ≥ n₀(c′). Assume it below N ∈ 𝒦.
  - By Markov a clocked call on m ∈ 𝒦 succeeds with probability ≥ ¾, so all r calls fail with probability ≤ 4^{−E}.
  - A call on m ∉ 𝒦 costs at most rβ(m) whatever happens.
  - Sub-inputs below n₀ are handled by alternation with trial division.
- **Fields.** In a field with AD_𝒦(A), each visited target's first 𝒦-good point survives its two clocked calls except
  with probability ≤ 2·4^{−E}. A run visits ≤ e^{(b+o(1))E} targets (rows × the e^{O(log²t)}-node trees of Lemma TS),
  so the field fails with probability e^{−Ω(E)}. With (H1), the expected number of fields is ≤ 2 + o(1). Work per row
  attempt is bounded by budgets in good and bad fields alike, as in NR.
- **Split.** Lemma S″ ([rec] §3.2) applied to S′ := 𝒰_A ∩ 𝒦 with |S′| ≥ |𝒰_A|/2 (H2) gives Pr ≥ |S′|²/(2φ(N)).
  - The shared-prime step uses k·C(|G| − 1, k − 1) ≤ (2k³/|G|)|𝒰_A| ≤ (4k³/|G|)|S′|, so the hypothesis becomes
    |G| ≥ 8k³2^ω.
  - This holds as in NR: |G| = e^{Θ(Q)} ≥ e^{O(t/E)}k³.
  - The acceptance exponent is unchanged, e^{−(1/b₁+o(1))E}.
- **(T).** Timeouts and repetitions are functions of sub-inputs and fresh coins, which are independent of the hidden
  square roots. So Lemma T_multi holds and each LV6 round is fruitful with probability ≥ ½.
- **Characters.** [U] §9's argument is degree-uniform.
  - LO Theorem 1.3 needs log x ≥ 10n_L(log d_L)² = O(d³t²), and s_χ ≫ that.
  - Lemma E holds for odd d.
  - s_χ ≫ b₁Q, so the character primes exceed every atom prime.
  - At d = δD′, log s_χ = 2d log d + O(d + log t) = (2δ/3 + o(1))E, so sampling and evaluation cost L[1/3, C_χδ]
    (codex; [nested] §7's Θ(δ) correction; LV Remark 6.7's 4δ/3).
  - SQ's sign enumeration 2^{d−1} = e^{o(E)} is part of LV6′.
- **Recursion.** As in NR, with the clocked sub-calls:
  T(N) ≤ e^{(max(2b, b+C_χδ)+o(1))E} + e^{(b+1/b₁+o(1))E}·rβ(√(2N)) + e^{(b+2s+o(1))E}·rβ(N^{τ_F+o(1)}).
  Since β(N^τ) = e^{c′τ^{1/3}E(1+o(1))} and c′ > c_NR-AA, the exponent max(2b, b + C_χδ, b + 1/b₁ + 2^{−1/3}c′,
  b + 2s + τ_F^{1/3}c′) is strictly below c′. So T(N) ≤ β(N)/4 for N ≥ n₀(c′), with uniform o(1)s. ∎

*Remarks.*
- With 𝒦 = Z_{≥2}, NR-AA is nested's Theorem NR without ERH, at the cost of the additive C_χδ.
- HAD is weaker than AD only in its quantifier. Its substance, the same-field supply at L-scale, is untouched by
  averaging over N (Prop NB, Lemma MP, Cor OD-AA). So NR-AA is a correct bookkeeping of what an almost-all L[1/3]
  would need, not progress toward it.

### 3.3 LV6 at growing degree and ERH (task 2, last question)

- ERH is not needed for a finite constant at α = 1/3 (Theorem NR-AA, item 3). It would only remove the additive C_χδ
  and allow polynomial-bit characters.
- Averaging over N would need a zero-free region or an exceptional-zero density bound for the family {K_A(√h)} of
  fields of degree 2d, with d → ∞, on average over the field law.
  - The known family results (Pierce–Turnage-Butterbaugh–Wood, Invent. Math. 2020; Lemke Oliver–Thorner–Zaman 2024)
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
5. **Best α: 2/5** (Theorem U, every N, internal chain). Conditional: 1/3 on 𝒦 under HAD + LV6′. Model: 1/3. No
   exponent in (1/3, 2/5) is reachable by averaging over N with any tool we know.

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

Sources: codex almostall13_1 §4 (primary sources read where stated), a literature-search subagent, and the Opus
literature critic (§9).
- **No rigorous almost-all-N (or balanced-semiprime) complete factorisation below L[1/2] was found**, nor
  L[1/2, c < 1] on average over N. The worst-case record is Lenstra–Pomerance, JAMS 5 (1992), L[1/2, 1 + o(1)]. This is
  a search result, not a nonexistence theorem.
- **ECM for almost all primes.**
  - Lenstra, Ann. Math. 126 (1987), Cor 2.8/Conj 2.10: the L_p[1/2, √2] bound needs a conjectural Hasse-interval
    smoothness estimate.
  - Hafner–McCurley, J. Algorithms 10 (1989), Thms 3.1–3.2: rigorous for all but x·exp(−½(log x)^{1/6}) primes, in
    time exp(O((log v)^θ)) for θ > 5/6.
  - Lenstra–Pomerance 1992 §7 (Thms 7.1–7.4): an almost-all-p ECM bound of the form exp((log p)^{6/7+o(1)}) (per the
    literature search; see the critic's verification in §9).
  - Both are far weaker than L_p[1/2]. Complete factorisation of almost all N must find the second-largest prime factor,
    which is typically N^{Θ(1)}, so neither gives anything below L_N[1/2] for complete factorisation.
- **Averaging in rigorous NFS/DL.**
  - Lee–Venkatesan, JNT 187 (2018): Thm 2.1, a possibly trivial congruence in L[1/3, (64/9)^{1/3}] for every N;
    Thm 2.3, for p ≡ q ≡ 3 mod 4 under their Conjecture 7.1. They average over the algorithm's random polynomial
    choices, not over N.
  - Pomerance 1987: worst case, L[1/2, √2]; its averaging is over auxiliary ECM inputs.
  - Kleinjung–Wesolowski, JAMS 35 (2022): every F_{p^n} with small p, from Frobenius structure absent for integers.
  - Barbulescu–Jouve (arXiv 2212.11724): conditional on Elliott–Halberstam-type hypotheses.
- **Smooth values of random polynomials and forms** (Barbulescu–Lachand, Math. Comp. 2017; BBDT 2012; Martin; Dartyge):
  fixed degree, or y a power of x; no L-scale joint-correlation results ([nested] §10).
- **Family zero-free regions:** fixed degree only (§3.3).

## 8. Candidate card

- **Mechanism.** Average over the input N, in addition to the field, to supply the free parameters that algebraic
  components lack ([nested] Lemma OD), for an almost-all L[1/3] (or α < 2/5).
- **Model.** Theorem NR / Theorem U architecture, design (a‴) and plain base-m, d ≍ D′.
- **Result.**
  - **Refuted for this model:** N adds no parameter beyond w (Prop NB; E1, E2).
  - Proved: Prop NB, Cors OD/SI/Info-AA, Remark W-AA, Lemma M1-AA, Lemma AV, Lemma MP, Prop HR.
  - Proved conditional on HAD + LV6′: Theorem NR-AA (no ERH).
- **Failure scope.** Every design in which N | F(g₀, g₁) and the descent inside a field depends only on the field.
  That is all NFS-type designs we know, including multi-field, constrained resampling and plain base-m.
- **Not excluded.**
  - Designs where N enters the descent other than through the field. None is known.
  - A genuine L-scale joint-smoothness theorem for random polynomial pairs. That would be new analysis, independent of
    N.
- **Next step.** The bottleneck is a precision problem (Lemma MP) or a same-field correlation problem (OD) for a random
  pair (F, G), and N-free. Averaging over N is not the missing ingredient.

## 9. Reviews

(Filled in after the referee and critic rounds; see below.)

## 10. Reproduction

```
.venv/bin/python -m src.almostall13.rank experiments/raw/almostall13_e1.json
.venv/bin/python -m src.almostall13.harness_e2 experiments/configs/almostall13_e2.json
.venv/bin/python -m src.almostall13.analyze_e2 experiments/raw/almostall13_e2.jsonl experiments/raw/almostall13_e2_summary.json
.venv/bin/python -m pytest -q tests/test_almostall13.py
```
