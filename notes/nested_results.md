# nested: ERH and nested descendability ("endogenous children") — results

Agent `nested`, phase 12, 2026-10-20.

Files:
- This note; pre-registration `notes/nested_prereg.md` (commit `7c17ffb`, Addendum A `e432e59`).
- Code: `src/nested/toy_e1.py`, `harness_e1.py`, `analyze_e1.py` (E1, diagnostic), `src/nested/barrier.py` (exact
  helper computations for §§2, 5, 6, 8). Tests: `tests/test_nested.py`.
- Data: `experiments/raw/nested_e1.jsonl`, `nested_e1_summary.json`, `nested_e1_calib.jsonl`; config
  `experiments/configs/nested_e1.json`.
- Codex: `codex/nested_1.{md,out,log}` (blind round; full report copied to `notes/nested_codex1_report.md`),
  `codex/nested_2.*` (adversarial referee, §12).

Everything about the project's own algorithm rests on Theorem A (`notes/thmA_proof.md`) and recurse
(`notes/recurse_results.md`), which are internal and unrefereed. External results are used as stated in their
sources; dBPW (arXiv 2512.01588v2) was read in the original (Thms 1.1, 6.4, 7.1 with its proof, 9.5, 13.2, 17.7,
§14, Defs 2.10–2.11). Labels: **proved**, **proved conditional on …**, **model computation**, **heuristic**,
**experimentally supported**, **open**.

## Verdict

**Plain language.** We did not find an ERH-conditional proof that the multi-level descent works, and so no rigorous
or ERH-conditional L[1/3], nor any fixed exponent below 1/2. No such result is in the literature either. What we
can say precisely:
- The de Boer–Pellet-Mary–Wesolowski sampler does not fit this job. Up to a factor e^{±1/4} its success
  probability is the density of the target family among ideals of norm up to about √|Δ_K| (a rereading of their own
  identities). For the fields we construct that norm is typically a power of N; we do not prove a lower bound for
  the field discriminant, and special fields such as x^d − 2 show none holds in general. Under the usual
  smooth-ideal heuristic, this sampler therefore costs L[1/2] per descent step and gives no contraction.
- The obstruction is not ERH-shaped. In the standard size model, a descent below L[2/5] has algebraic children:
  prime ideals of the field it works in. The row algebra closes such a child only by an element of the same field, so
  a whole "algebraic component" of the tree shares one field. Every argument we know that averages over the field
  then needs joint control of one polynomial's values at more points than the field has free coefficients. No such
  tool exists, with or without ERH. This is a barrier for the known methods, not an impossibility theorem.
- Under ERH plus two explicit hypotheses the L[1/3] algorithm is rigorous (Theorem NR, conditional):
  - AD, a strong all-target relation-supply statement in reduced-basis boxes;
  - LV6′, lv6's remaining lemmas at growing degree with stated costs, plus an admissible-field sampler.
  The tree, costs, recursion, transcript invariance and row algebra are proved here. ERH is used only for the
  character columns.
- **Best α under ERH from this work: α = 1/2 (no improvement).**

Precisely:
1. **dBPW extracted (§1).**
   - Any degree n ≥ 2; ERH for the Hecke L-functions of K; an LLL basis of O_K, or of an order R (§9.5) at
     discriminant [O_K:R]²Δ_K.
   - One-sided: β ∈ 𝔞 with (β)𝔞^{−1} ∈ 𝒮·𝒮_B, with probability ≥ δ_𝒮[r^n]/3 − ε for an *arbitrary* family 𝒮.
     δ_𝒮 is the local density among *all integral ideals* of norm ≈ r^n, and
     r^n ≥ 48^n b^{2n²/b} n^{7n/2} |Δ_K|^{3/2} (Thm 9.5's normalisation).
   - One sampling step; Thm 17.7 repeats it until the cofactor is L_{|Δ|}(1/2)-smooth.
   - Worst case over the input ideal; full-dimensional; no rational side.
2. **Theorem AF [proved, given an Arakelov-uniform input; a rereading of dBPW (34)–(48)].**
   - For Arakelov-uniform sampling in a box of radius ρ, the probability that the cofactor lies in 𝒮 is, up to
     e^{±1/4}, E_v[δ̃_𝒮(ρ^n e^{−v})], where v ~ Gamma(r₁ + r₂ − 1, 1) and δ̃ is the cumulative density among ideals of
     norm ≤ ρ^n e^{−v}. Here ρ^n ≥ (8n²)^n|Δ_K|^{1/2}.
   - This is a comparison, not a law of the cofactor norm. In the uniform ideal-count model the log-norm deficit is
     Gamma(r₁ + r₂, 1) (codex nested_2).
   - Corollary AF [heuristic step]: in fields with log|Δ_K| ≍ log N, Arakelov-type steps cost L[1/2] at best under
     the smooth-ideal model. This is dBPW's own L_{|Δ|}(1/2).
3. **Prop SI, the sparse-image dilemma [proved].**
   - For a fixed field the candidate values of a node form a sparse set. An adversarial good family of density
     1 − o(1) misses all of them, so no lemma that uses density alone can hold.
   - Averaging over the field restores full support, but only for families that do not depend on the field.
   - The existing arbitrary-family tools (Theorem K/S″, LP92 Thm 6.1 as used in LP §10, dBPW) are all counting plus
     a sampler that is uniform on a group. The size of that group sets the value size: N for Dixon, √|Δ_K| for
     Arakelov.
   - This does not show that the actual, recursively defined good sets are adversarial.
4. **Fresh fields (§4) [proved / structural / model].**
   - Rational children are integers and can be sent to an independent field (Prop R). A joint supply statement for
     the targets sharing that field is still needed.
   - Algebraic children are ideals of their field. The row algebra closes them only with elements of the same field
     (Lemma RW). Lemma FB excludes the tower transfer.
   - Constrained polynomial resampling F + G·∏L_{P_i}·H (codex nested_2) can refresh a field while keeping up to
     d − 1 retained points. That is the same parameter budget, now applied per algebraic component.
   - Prop R2 [model, explicit premises]: in the balanced-height, typical-size, random-integer model, designs without
     algebraic children do not beat the one-level 2/5.
5. **Lemma OD and the parameter budget [proved linear algebra; children counts are a typical-factorisation model].**
   - A degree-d field has d + 1 free coefficients.
   - In the typical-factorisation model a node with children ≤ q^λ has ≈ d/(2λ) large algebraic children. A thin
     tree whose values are mostly B-smooth is not excluded, but it is costlier.
   - Then depth 1 with λ < 1/2, depth 2 at any λ < 1, and any pair of depth-1 configurations (as a second moment
     needs) exceed the budget.
   - Independence-over-the-field moment methods (FM, Kern, LD*) then do not apply directly. One would need L-scale
     smooth-correlation estimates for overdetermined affine systems, which do not exist.
6. **Precision (§6) [model computation].**
   - In the naive backtracking model (an upper bound for that strategy only), per-level good-set densities
     ≥ (log N)^{−C} over O(log log N) levels give exp(O(E log log N)) = L[1/3 + o(1)].
   - A bounded geometric-mean density gives L[1/3, O(1)].
7. **What ERH buys (§7) [proved from GRH-Chebotarev].**
   - Lemma X_ERH: LV6's character columns at growing degree d, with poly(t)-bit character primes and no exceptional
     class. This closes the Lemma-X part of the LV6 growing-degree gap in recurse §3.5; the rest is Hypothesis LV6′.
   - Correction to recurse §3.5: unconditionally the Stark-threshold repair costs L[1/3, Θ(δ)] at α = 1/3, not
     L^{o(1)}. LV Remark 6.7 already asserts both points, without proof.
   - The adaptive one-sided lemma (codex nested_1) handles endogenous *one-sided* children under ERH, at the
     Theorem-AF norm scale.
   - ERH does not supply the needed two-sided supply theorem.
8. **Theorem NR [proved conditional on ERH, Hypothesis AD and Hypothesis LV6′].**
   - Hypothesis AD: with probability ≥ 1/2 over an admissible field law, every rational or degree-one target of
     norm q ∈ (B, e^{b₁Q}] has a point in its reduced-basis box (≤ (2e^{sE} + 1)² lattice points) whose two
     cofactors are max(B, q^λ)-smooth, with the large primes avoiding D_A.
   - Conclusion: a Las Vegas algorithm factors every N in expected L_N[1/3, c_NR + o(1)], with an explicit finite
     c_NR.
   - FM, LD*, Theorem K and the 2^ω issue all drop out.
   - AD is a *sufficient strengthening* of the visited-target supply that is really needed, not an equivalent of it.
     It is the NFS-DL descent heuristic in reduced-basis form; with sup-norm boxes it would fail on irregular
     lattices (both reviewers).
   - It is a pointwise statement about one field, for which no technique is known.
9. **E1 [experimentally supported, toy, pre-registered; tests sup-norm boxes, not the reduced-basis AD].**
   - Setup: 12 toy fields (d = 3, 44-bit N), 4800 parent nodes, 14333 algebraic children.
   - P1: endogenous children are one-step descendable at the same rate as matched uniform targets: 0.2957 vs
     0.2944, CI [−0.008, 0.012]. P2 and P4 hold.
   - P3 **failed**: its ratio was 1.71 against the pooled rate, because descendability varies strongly between fields
     (matched-uniform algebraic rates per field 0.110–0.587). Against each field's own rates the aggregate is 1.05
     (exploratory), consistent with a product model within a field.
10. **Literature (§10).**
    - No rigorous or ERH-conditional L[α < 1/2] factoring is known; LP92 L[1/2, 1] stands.
    - Lee–Venkatesan 2018 prove L[1/3] *congruences*; nontriviality only for p ≡ q ≡ 3 mod 4, under their
      Conjecture 7.1.
    - Rigorous L(1/3)-or-better discrete logs exist only in small or fixed characteristic (Kleinjung–Wesolowski
      2022, quasi-polynomial).
    - Class groups: dBPW rigorous (ERH) L(1/2). Jacobians: all L(1/3) results are heuristic.

---

## 0. Notation and setting

t := log N; E := t^{1/3}(log t)^{2/3}; Q := t^{2/3}(log t)^{1/3}; D′ := Q/E = (t/log t)^{1/3}. So E = D′ log t,
Q = D′² log t, D′Q = t, and L_N[1/3, c] = e^{(c+o(1))E}, L_N[2/3, c] = e^{(c+o(1))Q}.

The architecture is that of Theorem A (`thmA_proof.md` §3) with recursion (recurse §1, §3.1): rows r² ≡ R, a Gauss
split R ≡ u₀/v₀, complete factorisation of |u₀|, |v₀| by a recursive call, descent of every large prime factor by
degree-one special-q relations in one field A = (g₀, g₁, w, c) with rational side G = g₁x − g₀y and algebraic side
F_A of degree d, LV6 linear algebra, fruitfulness by Theorem DD. The size model of recurse §3.4 (d = δD′, factor
base B = e^{bE}, box inflation e^{sE}, children ≤ e^{λz} for a target of size e^z, O(log t) levels) is the target.

A *target* is a rational prime q (lattice {P : q | G(P)}) or a degree-one prime ideal 𝔮 = (q, r) of K = Q(α)
(lattice {(a, b) : a ≡ rb mod q}). The *children* of a box point P for a target are the rational primes > B of the
rational cofactor and the degree-one prime ideals > B of the algebraic cofactor. Rational children are integers;
algebraic children are ideals of K.

## 1. What dBPW prove (read in the original, arXiv 2512.01588v2, 19 Feb 2026)

**Setting and hypotheses.**
- K any number field, degree n ≥ 2 (no upper bound), discriminant Δ_K, residue ρ_K of ζ_K at 1.
- Input: an LLL-reduced basis of O_K (Thms 1.1, 9.5, 13.2). §9.5 (Lemma 9.14, Remark 9.15) gives a variant for an
  order R ⊂ O_K with the conductor in the modulus and radius |Δ(R)|^{3/(2n)}, |Δ(R)| = [O_K : R]²|Δ_K| (found by the
  codex blind round; I checked that §9.5 exists but did not audit it). So computing O_K is not an unavoidable
  obstruction; the discriminant scale only gets worse.
- ERH = the Riemann hypothesis for Hecke L-functions of K (§2.2). It enters through the random-walk theorem
  (Thm 6.4 / Cor 6.5: walk on Pic⁰_{K𝔪} with N ≈ 7n + 2log(1/ε) + log|Pic⁰| steps by uniform prime ideals of norm
  ≤ B = Õ(n²(n² (log log 1/ε)² + (log |Δ_K|N𝔪)²))), through prime-ideal sampling (Lemma 5.4), and through
  smooth-ideal density estimates in Part 2. Theorem 7.1 (the averaging statement) is unconditional.
- **Density notion (Def 2.11).** δ_𝒮[x] = min_{t ∈ [x/e^n, x]} |{𝔟 ∈ 𝒮 : N𝔟 ≤ t}|/(ρ_K t): the proportion of 𝒮
  among *all integral ideals* of norm ≤ t, for t in a window of ratio e^n below x. It is not a relative density
  among primes, nor among values of any polynomial.

**Theorem 7.1 (unconditional).** If the Arakelov ray class of the input is uniform on a coset of G ≤ Pic⁰_{K𝔪} and
ρ ≥ 8n²Γ_K|Δ_K|^{1/(2n)}N(𝔪)^{1/n} (Γ_K ≤ |Δ_K|^{1/n}), then for α uniform in Exp(𝔞) ∩ ρB_∞,
E Pr[(α)Exp(−𝔞) ∈ 𝒮 | ray condition] ≥ (N𝔪₀/φ(𝔪₀))·[Pic⁰:G]·δ_𝒮[ρ^n]/3.

**Theorem 9.5 (ERH).** For *every* input ideal 𝔟 (worst case) and any family 𝒮 (with the class condition
[d⁰(𝒮)] ⊆ G − [d⁰(𝔟)]), Algorithm 2 outputs β ∈ 𝔟 with (β)𝔟^{−1} ∈ 𝒮·𝒮_B^G, β ∈ τK_{𝔪,1},
|N(β)| ≤ N(𝔟)·B^{N}·r^n, with probability ≥ (N𝔪₀/φ(𝔪₀))[Pic⁰:G]δ_𝒮[r^n]/3 − ε, where
r = 48ω·b^{2n/b}·n^{7/2}·|Δ_K|^{3/(2n)}·N(𝔪₀)^{1/n}; time poly(b^b, log|Δ_K|, size(𝔟), log 1/ε, log N𝔪, index,
…). The success event is (β)𝔟^{−1} = 𝔠·∏_{i ≤ N}𝔭_i with known walk primes 𝔭_i and N𝔠 ≤ r^n.
Thm 1.1 is the specialisation G = Pic⁰, 𝔪 = O_K, except that the printed Thm 1.1 puts an extra N(𝔞)^{1/n} into r
(and Thm 9.5 instead puts N(𝔟) into the norm bound); the two printings are not literally consistent in v2, and we
use Thm 9.5's normalisation. Thm 13.2 is the version with success ≥ (N𝔪₀/φ(𝔪₀))δ_𝒮[r^n]/6. (For G ≠ Pic⁰ the
walk bound B carries a factor [Pic⁰:G]² and the walk length a term −log[Pic⁰:G] + log N𝔪.)

**§14 (Lemma 14.1, ERH: the prime-ideal count uses Lagarias–Odlyzko; needs B ≥ C(n + log|Δ|)^{2+ε},
A ≤ B/(4 log B), x ≥ Be^n).** A lower bound, combinatorial apart from the prime count, δ_{𝒮_{A,B}}[x] ≥ (4 log B)^{1−u}u^{−u}/(ρ_K B) for ideals composed of
primes with norm in (A, B] (counting products, as in LP92 Thm 6.1). Codex nested_1 gives the analogous bound for the
monoid of an *arbitrary* set of m ≥ δy/log y prime ideals of norm ≤ y: δ ≥ (δ/(u log y))^u/(ρ_K y).

**Theorem 17.7 (ERH).** If S contains all prime ideals of norm ≤ B_max = poly(L_{|Δ|}(1/2), L_{n^n}(2/3),
ρ_K^{cut}), Algorithm 3 writes any 𝔞 coprime to 𝔪₀ as αO_K·∏_{𝔭∈S}𝔭^{−v_𝔭} in expected time
poly(L_{|Δ|}(1/2), L_{n^n}(2/3), size, ρ_K^{cut}, …). One sampling step, repeated until the cofactor is smooth.

**Answers to the task's questions.**
- *Degree range:* unrestricted. *Density:* among all ideals of norm ≈ r^n (local density, Def 2.11).
- *One-sided vs two-sided:* strictly one-sided. The output β ∈ 𝔞 is a general element of a rank-n lattice. The NFS
  set {a − bα} is a rank-2 subgroup of K, not an ideal (an ideal contains cO_K for some c ≠ 0, so has rank n), so it
  cannot be an input; conditioning dBPW's output on lying in it destroys the bound (codex: O(H²) of Θ(H^n) points).
- *Rational side:* none. ψ(β) = Σβ_i m^i mod N has no small rational counterpart for a general β.
- *Level structure:* one step. The cofactor bound N𝔠 ≤ r^n does not shrink with N𝔞 (the N(𝔞) cancels), so iterating
  gives no contraction.
- *Cofactor scale for NFS fields.* For a degree-n polynomial with log‖f‖ = O(Q) (the size-model fields),
  log|disc f| = O(n log‖f‖ + n log n) = O(t) (Hadamard on the Sylvester matrix), so log|Δ_K| = O(t). The typical
  value is log|Δ_K| ≍ t (the polynomial discriminant is ≍ ‖f‖^{2n−2}); we do not prove a lower bound for the field
  discriminant of the random fields. dBPW's certified cofactor scale is then N^{3ζ/2+o(1)} when log|Δ_K| = ζt
  (codex (14), with b = 2 the extra terms are o(t) at n ≍ D′).

## 2. Theorem AF: the cofactor law of Arakelov-uniform sampling

**Theorem AF [proved; it is dBPW's (46)–(48), restated as a law of the cofactor norm].** Let K have
r* := r₁ + r₂ − 1 ≥ 1 (any field but Q and imaginary quadratic ones), trivial modulus. Let the Arakelov class of 𝔞 ∈ Div⁰_K be uniform on Pic⁰_K, let
ρ ≥ 8n²Γ_K|Δ_K|^{1/(2n)}, and let α be uniform in Exp(𝔞) ∩ ρB_∞. For any set 𝒮 of integral ideals put
δ̃_𝒮(x) := |𝒮(x)|/(ρ_K x), X := ρ^n. Then

  e^{−1/4}·𝔼_v[δ̃_𝒮(Xe^{−v})·1_{v ≤ log X}] ≤ Pr[(α)·Exp(−𝔞) ∈ 𝒮] ≤ e^{1/4}·𝔼_v[δ̃_𝒮(Xe^{−v})·1_{v ≤ log X}],

with v ~ Gamma(r*, 1) (density v^{r*−1}e^{−v}/(r*−1)!).

*Proof.* For a single ideal 𝔠, p_{𝔞,𝔠} := Pr_α[(α)Exp(−𝔞) = 𝔠] equals |Exp(𝔞 + d(𝔠))^× ∩ ρB_∞| divided by
|Exp(𝔞) ∩ ρB_∞| (dBPW Lemma 7.4 with trivial modulus; Exp(·)^× is the set of generators). The denominator lies in
[e^{−1/4}, e^{1/4}]·2^{r₁}(2π)^{r₂}ρ^n/√|Δ_K| for every 𝔞 (dBPW (34), from Lemma 2.8, using the hypothesis on ρ).
The numerator depends only on the Arakelov class, and its expectation over a uniform class is exactly
w·Vol(S_{n log ρ − log N𝔠})/(h_K·Vol(T)) = w(n log ρ − log N𝔠)_+^{r*}/(r*!·h_K R_K) (dBPW (38)–(42), Lemma 2.4:
the expected number of points of a uniformly shifted log-unit lattice in a simplex is volume over covolume; these are
equalities). With the class number formula ρ_K = 2^{r₁}(2π)^{r₂}h_K R_K/(w√|Δ_K|) this gives
E_𝔞 p_{𝔞,𝔠} ∈ [e^{−1/4}, e^{1/4}]·C(N𝔠)/(ρ_K X), C(N) := (log(X/N))_+^{r*}/r*! (dBPW (44), (46)). Summing over
𝔠 ∈ 𝒮 and applying Abel summation as in dBPW (47):
Σ_{𝔠∈𝒮}C(N𝔠) = ∫_1^X |𝒮(t)|(log(X/t))^{r*−1}/((r*−1)!·t) dt; the substitution t = Xe^{−v} turns
Σ_{𝔠∈𝒮} C(N𝔠)/(ρ_K X) into 𝔼_v[δ̃_𝒮(Xe^{−v})1_{v≤log X}]. ∎

**Remarks.**
- dBPW state (46) two-sided but use only the lower half (with δ_𝒮[X] ≤ δ̃_𝒮(t) for t ∈ [X/e^n, X], and
  Pr[v ≤ n] ≥ 1/2). The upper half shows their bound is tight up to e^{±1/4} and the Gamma weighting: *the method
  succeeds as often as the cumulative density of 𝒮 among ideals of norm ≤ Xe^{−v}, averaged over v.* No cleverness in choosing 𝒮 or ρ changes the norm scale X ≥ (8n²)^n√|Δ_K|.
- For 𝒮 = {N𝔠 ≤ T}: Pr[N(cofactor) ≤ T] ≤ e^{1/4}(n log ρ)^{r*}·|{𝔠 : N𝔠 ≤ T}|/(r*!·ρ_K·ρ^n) (from
  C ≤ (n log ρ)^{r*}/r*!).
- **AF is a comparison, not a law of the cofactor norm (codex nested_2, issue 1).** The actual law weights each ideal
  𝔠 by C(N𝔠)/(ρ_K X). In the uniform ideal-count model |𝒮(u)| = ρ_K u, the log-norm deficit log(X/N𝔠) has the
  Gamma(r* + 1, 1) density, so Pr[N𝔠 ≤ Xe^{−V}] ≈ e^{−V}Σ_{j≤r*}V^j/j! (`barrier.af_weight`). The auxiliary
  variable v ~ Gamma(r*, 1) in AF is not that deficit; the extra exponential comes from the cumulative cutoff. This
  changes nothing at L-scale.
- The theorem is about the averaged (Haar) step; dBPW's ERH random walk makes the actual input ε-close to it, so the
  same identity holds for Algorithm 2 up to ±ε.

**Corollary AF [heuristic in its quantitative step; three explicit premises].** Premises: (i) log|Δ_K| ≥ ζt, an
*assumption* (the field discriminant can be far below the polynomial discriminant because of the index; x^d − 2 with
N | m^d − 2 has |Δ_K| ≤ d^d2^{d−1}, so no lower bound holds for all usable fields); (ii) the smooth-ideal density
*upper bound* below, which dBPW do not prove (they prove lower bounds); (iii) the sampler is Arakelov-uniform *and*
samples uniformly in a box satisfying the radius hypothesis. Under these, every such step succeeds for a family 𝒮 with probability ≈ the local density of 𝒮 at
norm Xe^{−v}, X ≥ |Δ_K|^{1/2}. Moving the effective norm down by e^{−V} costs Gamma weight ≈ e^{−V}V^{r*}/r*!, while
for 𝒮 = 𝒮_y (y-smooth ideals) the density gains only a factor ≈ e^{V log u/log y} with log u/log y = o(1); so the
weighted average is governed by norms X·e^{−O(n log n)} = N^{ζ/2+o(1)} (n ≤ t^{1/3+o(1)}). If smooth ideals of norm X have density
≤ u^{−(1+o(1))u} (u = log X/log y), which is the standard heuristic and the value dBPW's own analysis uses, then
making such a cofactor e^{bΛ_α}-smooth costs exp((ζ(1−α)/(2b) + o(1))t^{1−α}(log t)^α) = L[1 − α]. Balancing
against a factor base L[α] forces α ≥ 1/2, and iterating the step gives no contraction (the cofactor scale does not
depend on the input norm). This is exactly dBPW's L_{|Δ|}(1/2) (codex nested_1 §3(c) independently). It is a statement
about this sampler, not about all possible samplers or designs.

**Remark (one-sided rows; codex nested_1 (21)–(22)).** With the evaluation map ψ: R → Z/N of kernel 𝔫 (N𝔫 = N), the
ray condition β ≡ R = r² (mod 𝔫) gives valid one-sided rows r² = ψ(β), with no rational side needed. But the modulus
puts N into the radius: r^n ⊇ |Δ|^{3/2}·N (N√|Δ| for the Haar statement). This direct implementation is "Dixon's
method in K", at sampling scale ≥ N√|Δ|, hence L[1/2] at best under the smoothness model; these are scales and
envelopes, not pointwise lower bounds on every output norm.

## 3. Prop SI: the sparse-image dilemma

**Proposition SI [proved].** Fix a node: a finite candidate set Ω (the box), a field A, and the value map
Φ_A: Ω → 𝓘_X (ideals, or integers, of norm ≤ X; |Ω| ≤ W). Let 𝒮(A) ⊆ 𝓘_X be any "good family" allowed to depend
on A.
- (i) *Fixed field.* If 𝒮(A) := 𝓘_X ∖ Φ_A(Ω), then Pr_{P∈Ω}[Φ_A(P) ∈ 𝒮(A)] = 0, while the relative density of
  𝒮(A) is ≥ 1 − W/|𝓘_X|, and its local density satisfies δ_{𝒮(A)}[x] ≥ δ_{all}[x] − e^nW/(ρ_K x) (codex (19)).
  Hence no lemma of the form "for every family of relative density ≥ η the node hits it with probability ≥ f(η) > 0"
  holds for a fixed field when W = o(|𝓘_X|).
- (ii) *Averaging over the field.* If 𝒮 does not depend on A and, for each P, the law of Φ_A(P) over A has density
  ≥ c with respect to the uniform law on some set 𝓘′ (a full-support sampler on 𝓘′), then
  E_A Pr_P[Φ_A(P) ∈ 𝒮] ≥ c·|𝒮 ∩ 𝓘′|/|𝓘′|.
- (iii) In nested descent the good family of a child is "descendable in the field the child lives in". For an
  algebraic child that field is A itself, so (ii) does not apply and (i) shows that density alone cannot suffice.

*Proof.* (i): Φ_A(Ω) has at most W elements; delete them from every cumulative count. (ii): linearity. (iii): an
algebraic child is an ideal of K_A (Lemma FB below). ∎

**The principle behind every existing arbitrary-family tool.** Each is counting plus a sampler that is uniform on a
group; the group's size fixes the value size.
- Theorem K / Lemma S″ (recurse §3.2; codex (26)): R uniform on (Z/N)^{×2}; every coprime same-class pair (u, v) is hit
  with probability 2^ω/φ(N); value size |u₀v₀| ≈ N.
- LP92 Thm 6.1 (p. 499, read in the original): a pure count of integers ≤ x composed of primes from an arbitrary set
  𝒫, under y ≤ exp((log x)^{1/2}(log log x)^η) and two mass conditions on 𝒫. LP use it with forms uniform over a class
  group (§§5, 10), at value size ≈ √|Δ|. Its y-range already excludes x = e^{Θ(t)}, y = e^{Θ(Q)} (codex).
- dBPW: Arakelov-uniform input, value size ≥ √|Δ_K| (Theorem AF).
In an NFS node the candidate set has W = e^{O(E)} points while the values have size X = e^{Θ(z d)} with z ≥ bE, so
W/X → 0 and only route (ii) is available.

**Class obstruction (codex (18)).** If the good primes all lie in a proper subgroup H of Cl(K) and [𝔮] ∉ H, no relation
(β) = 𝔮∏𝔭^{e} exists. Relative prime density cannot see this; including all small primes removes this particular
obstruction but proves nothing two-sided.

## 4. Fresh fields: what can and cannot be made exogenous

**Lemma FB (algebraic children are field-bound) [proved].** Let K₁ ⊂ K₂ with [K₂ : K₁] = e ≥ 2, K₂ = K₁(γ), and let
𝔮 ⊂ O_{K₁} be a prime of residue field F_q such that γ − c ∉ 𝔮O_{K₂} for every c ∈ Z. If a, b ∈ Z and
a − bγ ∈ 𝔮O_{K₂}, then q | a and q | b. *Proof.* O_{K₂}/𝔮O_{K₂} is an F_q-algebra in which 1 and γ̄ are F_q-linearly
independent (otherwise γ ≡ c mod 𝔮O_{K₂} for some c ∈ Z), so ā − b̄γ̄ = 0 forces ā = b̄ = 0. ∎ So degree-one
elements of an extension cannot carry the whole extended ideal 𝔮O_{K₂}; they can contain one prime 𝔔 | 𝔮, which
leaves the other primes above 𝔮 uncancelled (codex (23): passing to the norm q loses the conjugates). In unrelated
fields the ideal 𝔮 of K₁ has no meaning at all. So algebraic children have to be descended in their own field as far
as we can see; this last sentence is structural (no transfer construction is known to us), not a theorem.
(The hypothesis on γ fails only for 𝔮 dividing a fixed index-type quantity; it is not needed for the conclusion that
the ideal is field-bound, only for the tower variant.)

**Prop R (rational children can be made exogenous) [proved, structural].** A rational child is an integer. A relation
in field A_j (rational side G_j(P)/g₁^{(j)} ≡ ψ_j(a − bα_j) mod N) and one in an unrelated field A_{j+1} can be
multiplied, and the final linear algebra has one block of algebraic columns per field and shared rational columns
(Coppersmith-type multiple-field bookkeeping; each field needs its own LV6 square root). So a design in which level
j + 1 uses a field independent of levels ≤ j sees its rational targets as exogenous: conditional on the history, the
target q is fixed and the field A_{j+1} is fresh, so route (ii) of Prop SI is *available* for q; the required
value-distribution minorisation and a simultaneous supply statement for all targets sharing A_{j+1} remain
hypotheses, and the number of field blocks must be charged. Two caveats (codex
nested_1 §3(d)): all targets at one level share A_{j+1}, so a *joint* statement over the m targets of a node is still
needed; and the bound must hold conditionally on the history that selected q.

**Constrained polynomial resampling (codex nested_2, issue 3) [proved algebra].** Keep a parent point P₀ = (a₀, b₀)
and put L₀(x, y) := b₀x − a₀y. For any binary form H of degree d − 2, F_H := F + G·L₀·H satisfies F_H(P₀) = F(P₀) and
F_H(g₀, g₁) = F(g₀, g₁) ≡ 0 (mod N). So the parent's relation, its residue identity and its children's labels
(p, a₀/b₀ mod p) survive in the new field K_H, whose remaining coefficients are fresh randomness. With k retained points
one uses ∏L_{P_i} and H of degree d − 1 − k. This does make algebraic labels exogenous relative to H. But Lemma RW closes
an algebraic target only by an element *of the field in which it was opened*. So all nodes linked by algebraic edges
(an *algebraic component*: an algebraic target, the element that opened it, the elements that close it, and so on,
down to rational targets) must live in one field, and every point of the component must be retained. The budget
d + 1 − k therefore applies to a whole algebraic component, which is the parameter budget of Lemma OD again. Rational
targets cut components (Prop R). This route is not excluded by our arguments for components with fewer than d points.

**Prop R2 (no gain without algebraic children) [model computation with explicit premises].** Premises:
- (M1) the *counting heuristic*: usable pairs (f, g) with a common root mod N satisfy (d + 1)log H_f + 2 log H_g ≳ t
  (H_f = ‖f‖, H_g = ‖g‖). The only rigorous constraint is H_f H_g^d ≥ N/(d + 1), and with it alone the top level
  could avoid algebraic children (d = 2, H_f = O(1), H_g ≈ √N). The first draft's "‖f‖ ≥ N^{1/(d+1)}" was wrong; f = x³ − 2
  has a root mod suitable N with height 2 (both reviewers).
- (M2) typical sizes on regular lattices: the rational cofactor has log-size ≈ log H_g + sE − z/2, the algebraic value
  ≈ log H_f + (d/2)z.
- (M3) random-integer smoothness: log-size X is e^y-smooth with probability exp(−(1 + o(1))(X/y)log(X/y)).
- (M4) the split forces level-1 targets z₁ ≍ t^{1−α}(log t)^α.
Each node may use its own field (any d, H_f, H_g), but no node may have algebraic children (algebraic values B-smooth).
Write Λ := log B and U for the largest u with u log u ≤ Λ (so U ≍ Λ/log t); a node costs L[α] only if each of its
smoothness events has u ≤ U.
- *Window.* Put A := Λ²U²/t and W := t/(8U). If log H_g ≤ t/4, (M1) gives log H_f ≥ t/(2(d + 1)), so
  log H_f + (d/2)z ≳ √(tz), and B-smoothness within budget forces z ≲ A. If log H_g > t/4 and z ≤ t/4, the rational
  cofactor has log-size ≥ t/8, so its child bound y satisfies y ≥ t/(8U) = W. So targets with z ∈ (A, W) are
  infeasible.
- *No jump over the window.* For a feasible node above W, requiring *all* children below e^A means that a log-size
  ≥ t/8 cofactor is e^A-smooth, i.e. u ≥ t/(8A) = t²/(8Λ²U²). That exceeds U exactly when the window is non-empty.
  So, within budget, every node above the window has a child of log-size ≥ A, which is either in the window
  (infeasible) or above it. Leaves are ≤ Λ < A (as α > 1/3). So a finite tree of cost L[α] cannot exist.
- The window (A, W) is non-empty iff Λ²U³ ≪ t², i.e. 5α < 2.
So in this model designs without algebraic children have α ≥ 2/5, the one-level floor. Not covered: structured or rare
polynomial selection outside (M1), atypically small values, designs that keep large columns for recombination,
absorb bad primes in H, or change the search dimension. None of these is known to give a faster algorithm.

## 5. Prop OD: overdetermination of same-field configurations

**Lemma OD [proved].** Fix (g₀, g₁, w), so F_A(P) = F̂(P) + G(P)·⟨c, v(P)⟩ with v(P) = (a^{d−1}, a^{d−2}b, …, b^{d−1})
and c ∈ Z^d. For s > d points P₁, …, P_s with G(P_i) ≠ 0 there are integers λ_i, not all 0, such that
Σ_i λ_i F_A(P_i)/G(P_i) = Σ_i λ_i F̂(P_i)/G(P_i) for every c. Over the full law (g, w, c) the vector
(F_A(P_i))_{i≤s} ranges over the integer points of an image of a (d + 3)-dimensional family, so for s > d + 3 its
law is supported on a variety of positive codimension. *Proof.* s vectors in Q^d are linearly dependent; the second
statement is a dimension count. ∎ (Cleaner, codex nested_2: a degree-d binary form has d + 1 coefficients, so at
s > d + 1 fixed points the evaluation vector satisfies a fixed linear relation, whatever the parametrisation.)

**Children count [typical-factorisation model, not a consequence of the descent definitions].** The cap p ≤ e^{λz}
only gives Ω(C_large) ≥ log(C_large)/(λz) for the part C_large of the algebraic cofactor above B, counted with
multiplicity. There is no lower bound for log C_large: factor-base primes may carry the mass, repeated primes give one
target, and memoisation shares targets (codex nested_2, issue 5). The first draft's "m ≥ d/(2λ)" was wrong as a
deduction. In the typical-factorisation model (log|F(P)| ≈ (d/2)z + log‖f‖, and a random e^{λz}-smooth value has
most of its mass in ≈ u distinct primes near the top), the node has m ≈ d/(2λ) large algebraic children
(`barrier.min_children`). Thin trees, whose values are mostly B-smooth, are not excluded; they cost more.
Terminology: depth 1 = a node with its children; depth 2 = with grandchildren. Hence, in this model:
- λ < 1/2 (in particular any design with exponent contraction, as below 2/5): m > d, so already depth 1 has more
  values than the field has parameters;
- λ ∈ [1/2, 1) (constant contraction at α = 1/3): depth 1 (≈ d/(2λ) + 1 values) fits, but a pair of depth-1
  configurations (needed for a second moment) does not, and depth 2 does not (`barrier.overdetermined`).

**Consequence (method barrier).** FM, Kern and LD* work by fixing the points and using that the values are
independent-ish arithmetic progressions in c. For depth 1 with λ < 1/2, for depth 2 at any λ, or for the second
moment of a depth-1 configuration, the values are affinely dependent in c; a proof would need an *L-scale smooth
correlation estimate for an overdetermined system of affine forms with huge coefficients*. The literature has
nothing of the kind: Harper (Compositio, arXiv 1408.1662) handles a + b = c for y ≥ (log x)^K; Matthiesen–Wang (ANT
2025) need Ψ(N, y) > N/log⁸N; Lachand needs fixed u; BBDT's binary-form results need y ≥ x^{d−2+ε} (x^{1/√e+ε}
for cubics). This is a barrier to the method, not an impossibility theorem: one random object can support many
correlated searches (codex).

**Information version [proved, trivial].** After d successful nodes with independent v(P_i) and their values revealed,
c is determined. Any later node's success is a deterministic function of (A, target). So beyond depth 1 every needed
statement is pointwise in the field, i.e. about smooth values of one fixed binary form in special-q lattices at
L-scale; even primes in such values are known only for special forms (Heath-Brown x³ + 2y³; Maynard's incomplete
norm forms need n ≥ 4k, i.e. at least 3n/4 of the n variables; n ≥ 22k/7 for Q(θ^{1/n})).
Barbulescu–Lachand (Math. Comp. 86, 2017) study friable values of NFS binary forms, again at power-type y.

## 6. Precision: backtracking and union bounds

**Prop BT [model computation].** In a descent whose good sets have density p_j at level j, a node with m children
must retry until all children succeed; with independent children the expected work multiplies by
m·p_{j+1}^{−m} per level. Over O(log t) levels with m ≍ D′ (α = 1/3):
log(overhead) = Σ_j (log m + m log(1/p_j)).
- p_j ≥ (log N)^{−C}: log(overhead) = O(D′·(log t)²) = O(E·log t) (attained when p_j ≍ t^{−C} throughout), so the
  model bound is exp(O(E log log N)) =
  L[1/3]^{O(log log N)}, i.e. exp((log N)^{1/3}(log log N)^{5/3+o(1)}), not L[1/3, c] (`barrier.backtrack_log_cost`;
  codex (b): a t^{−c} loss per level becomes e^{−O(log² t)}, and raised to 2k it costs e^{O(E log t)}).
- bounded average of log(1/p_j) (geometric-mean density bounded below): log(overhead) = O(D′ log t) = O(E), so
  L[1/3, O(1)]. Fully correlated children would cost p_j^{−1} instead of p_j^{−m}; the model assumes independence.
So with the nested lemma requested in the task ("relative density ≥ (log N)^{−C}") the naive backtracking strategy
gives at most L[1/3 + o(1)]; to reach L[1/3, c] by it one needs constant densities, or AD (p_j = 1). A smarter
global selection argument is not excluded by this model computation.

**Union-bound route to AD.** To get "every target is good" from field moments by a union bound one needs
Pr_A[τ dead] ≤ e^{−b₁Q}/(d + 1) per target. The particular *unconditioned* Chebyshev/Paley–Zygmund estimate gives
only Pr[dead] ≤ Var/E², which is Θ(1) here because the count of good points has Θ(1) relative variance across
fields (local factors such as root counts of f mod small primes; E1's per-field rates 0.11–0.59 are a toy
illustration, not an asymptotic theorem). Large relative variance does not itself imply a large dead probability
(codex: Z = Mc_A with c_A ∈ {1, 2}), and conditioning on local factors, or higher moments, are not excluded by this
remark. We know no way to carry them out at the required precision; joint 2k-point smoothness asymptotics are not
available at L-scale.

## 7. What ERH does buy

**Lemma X_ERH [proved from GRH-Chebotarev].** Assume GRH for ζ_K and for the Dedekind zeta functions of the fields
L_h = K(√h), h ∈ H ∖ {1} (lv6 §3.2–3.3 notation; [L_h : Q] = 2d, |Δ_{L_h}| ≤ Δ̄ with log Δ̄ = O(d·t)). Sample
first-degree primes 𝔭 = (p, s) of K with N𝔭 ∈ (x/2, x] and p ∤ D as in lv6 (rejecting p | D, i.e. at most d·ω(D) prime ideals),
with x the least integer such that √x ≥ C₀(d + log Δ̄)²(log Δ̄)² and log x ≥ b₁Q + 1. Then for every h ≠ 1, Pr_𝔭[χ_𝔭(h) = −1] ≥ 1/(5d) for large N, and with
t_χ := ⌈5d(ln|H| + 2 ln t)⌉ characters the joint kernel meets H only in {1} with probability ≥ 1 − t^{−2}. No
exceptional class (lv6 Lemma E) occurs.

*Proof.* Lagarias–Odlyzko 1977, Thm 1.1 (GRH version): for L/K Galois with group G and class C,
|π_C(x) − (|C|/|G|)Li(x)| ≤ c₁((|C|/|G|)x^{1/2}log(d_L x^{n_L}) + log d_L), c₁ absolute. Apply to L_h/K, G = C₂,
C = {σ}: the number of primes of K with Frobenius σ and norm in (x/2, x] is
≥ ½(Li(x) − Li(x/2)) − 2c₁(√x(log Δ̄ + 2d log x) + log Δ̄) ≥ (x/(4 log x))(1 − o(1)) for the stated x (C₀ large).
Subtract the ≤ d√x primes of degree ≥ 2 and the ≤ log₂Δ̄ ramified ones, both o(x/log x). Upper bound: first-degree
primes with norm in (x/2, x] number ≤ d(π(x) − π(x/2)) ≤ 1.1dx/log x (prime number theorem). The ratio is
≥ (1 − o(1))/(4.4d) ≥ 1/(5d).
Independence and a union bound over |H| classes give the kernel statement; under GRH no h has an exceptional real
zero, so H₀ = {1}. The condition log x ≥ b₁Q + 1 makes every character prime larger than every prime > B that divides an
atom of a multi-level row (such primes are targets, of norm ≤ e^{b₁Q}), and the rejection of p | D keeps character
primes away from D; primes ≤ B are below x. Sampling and evaluation cost poly(log x) = poly(t). ∎

The character primes have poly(t) bits (x ≥ e^{b₁Q} because of the atom condition; without it √x polylogarithmic in N
would do), so LV6's character step works at any growing degree under ERH.
**Correction to recurse §3.5:** without ERH, the Stark threshold ν₁ ≍ 1/((2d)!log Δ̄) forces log x ≳ (2d)!·log Δ̄,
i.e. characters of (2d)!-bit size; at α = 1/3 with d = δD′ this is e^{Θ(δ)E}·poly, an L[1/3, Θ(δ)] factor (still
finite), not L^{o(1)} as recurse §3.5 says (that statement is right for α > 1/3). LV's own §6 (unconditional,
Stark–Brauer–Siegel with an exceptional field set) is the unconditional alternative; we did not transfer it.
**Prior statement (literature critic):** LV Remark 6.7 (p. 31) already says that their characters cost
L(1/3, 4δ/3 + o(1)) unconditionally and that "conditional on GRH these L_n(1/3) bounds become polynomial in log n".
Lemma X_ERH is the proof of that remark in lv6's setting, and our Θ(δ) correction matches their 4δ/3.

**Can ERH be weakened or averaged?** For Lemma X_ERH: yes in principle (LV §6, unconditional). Family-averaged
zero-free regions (Pierce–Turnage-Butterbaugh–Wood, Invent. Math. 2020; Lemke Oliver–Thorner–Zaman 2024) are for
families of *fixed* degree and fixed Galois group; nothing is known for growing-degree families such as ours. For
nested descent ERH is not the missing input at all (§§3–5), so weakening it changes nothing there. For dBPW the
Hecke ERH drives the random walk; averaging over our random fields would need the spectral gap for growing-degree
Arakelov groups on average, which is open, and in any case Theorem AF caps the outcome at L[1/2].

**The adaptive one-sided lemma (codex nested_1 §2) [proved under ERH, from Thm 9.5].** If the target ideal and the
family 𝒮 are arbitrary functions of the history and δ_𝒮[r^n] ≥ η for every history, then with fresh coins and
ε ≤ η/6 each attempt succeeds with probability ≥ η/6 conditionally on the history. This genuinely handles endogenous
*one-sided* children; its price is the Theorem-AF norm scale.

## 8. Theorem NR: ERH + AD ⇒ L[1/3]

Parameters: constants δ, λ ∈ (0, 1), b, b₁, s > 0 with τ_F := δb₁ < 1. d := the least odd integer ≥ δD′;
B := ⌊e^{bE}⌋. Field law 𝔉: one uniform sampler for all N, running in expected time e^{o(E)} and outputting either
a proper factor of N (a *factor exit*) or an *admissible* field: gcd(g₁·f_d·disc(f), N) = 1, ψ(g′(β)) a unit, and A = (f, g₀, g₁) with f irreducible of degree d
(with a certificate), F(g₀, g₁) ≡ 0 (mod N), gcd(g₀, g₁) = 1, and log‖f‖_∞, log max(g₀, g₁) ≤ κQ for a constant
κ (e.g. design (a‴) with M = ⌊N^{1/(d+1)}⌋, w ≈ M, H ≈ M). D_A := N·g₁·f_d·disc(f).

*Targets* 𝒯(A): rational primes q and degree-one prime ideals (q, r) (r a simple root of f mod q) with
B < q ≤ e^{b₁Q} and q ∤ D_A. *Box* of a target of norm q = e^z: let (u₁, u₂) be the Gauss-reduced basis of its
lattice (deterministic); the box is {P = iu₁ + ju₂ : |i|, |j| ≤ ⌊e^{sE}⌋, P = (a, b) with gcd(a, b) = 1,
G(P)F(P) ≠ 0}, normalised to b ≥ 0. It has ≤ (2e^{sE} + 1)² points, and since ‖u₁‖‖u₂‖ ≤ (2/√3)q,
‖P‖_∞ ≤ e^{sE}(‖u₁‖ + ‖u₂‖) ≤ 3e^{sE}q (≈ e^{sE}√q for regular lattices). *Why not a sup-norm box (proof critic
issue 1):* for an irregular lattice (a primitive u₁ with ‖u₁‖ < √q·e^{−sE}/√2, e.g. u₁ = (0, 1) for q | g₀) the box
|a|, b ≤ √q·e^{sE} meets only the line Zu₁, so it has one primitive point, which is good with probability bounded
away from 1; there are e^{Θ(E)} such targets per field, so AD would fail under its own heuristic. *Good point:* every prime p dividing the
target's two cofactors (G(P)/q or G(P); F(P)/q or F(P)) satisfies p ≤ B, or p ≤ e^{λz} and p ∤ D_A.

**Hypothesis AD(δ, λ, b, b₁, s; 𝔉).** For all large N, Pr_𝔉[factor exit, or A admissible and every τ ∈ 𝒯(A) has a
good point] ≥ 1/2. (Admissibility is part of the event: AD-good raw fields that fail NEST's gcd checks do not count;
codex nested_2, issue 8.) AD is a *sufficient strengthening* of the visited-target supply that the descent needs.
Cofactors are taken in absolute value and "good" means max(B, q^λ)-smooth with the large primes avoiding D_A.

**Hypothesis LV6′.** lv6's Lemmas C, H, SQ and IR (good primes, residual square-class group, square root in Z[β],
irreducibility certificate) hold at degree d ≍ D′ with cost e^{o(E)} (polynomial bounds in d and t, except SQ's sign
enumeration 2^{d−1} = e^{Θ(D′)} = e^{o(E)}), with SQ's Hensel prime taken as below. (UNVERIFIED bookkeeping; the
only growing-degree gap identified in recurse §3.5 was Lemma X, which Lemma X_ERH replaces.)

**Theorem NR [proved conditional on ERH, AD and LV6′].** There is a Las Vegas algorithm that outputs the complete
factorisation of every N in expected time L_N[1/3, c + o(1)] for every
c > c_NR := max(2b, (b + 1/b₁)/(1 − 2^{−1/3}), (b + 2s)/(1 − τ_F^{1/3})).

*Algorithm NEST* (run inside the recurse wrapper: preprocessing, trial division to B, perfect powers, AKS, and
alternation with trial division instead of an N₀ branch, recurse §1.1 item 3 / auditRA):
1. Sample A ∼ 𝔉 (certificates; gcd checks of g₁, f_d, disc f, content with N; a nontrivial gcd is a factor exit).
2. Rows. Repeat: r uniform in [1, N−1] (gcd check), R := r², Gauss-reduce to (u₀, v₀); completely factor |u₀|, |v₀|
   by recursive calls with fresh coins; accept iff every prime factor p satisfies p ≤ B or (p ≤ e^{b₁Q}, p ∤ D_A).
   For each accepted row, run DESC on every prime factor > B; if any DESC returns FAIL, discard A and all rows and go
   to step 1.
   DESC(τ): scan the box of τ in a fixed order; for each P, gcd-check G(P), F(P) with N, completely factor both by
   recursive calls; at the first good point, run DESC on every child (prime > B of the cofactors; algebraic children
   are (p, a·b^{−1} mod p)) and return the tree; if the box has no good point return FAIL.
3. When M_rows := #columns + 12 rows are collected, run LV6 with Lemma X_ERH (columns: sign, v_p mod 2 for p ≤ B,
   parity column, (p, s) for p ≤ B, p ∤ D, t_χ characters), the rows being given by Lemma RW below; compute
   X = ∏_{i∈S} r_i, Y, and gcd(X − Y, N). **Hensel prime:** one of the sampled character primes p (> e^{b₁Q}, p ∤ D, so Γ is a unit at an unramified prime;
finding "the least" such prime would have no proved cost bound, codex nested_2 issue 14). (lv6 takes
   p ∈ (B, 2B], and its Lemma SQ needs Γ to be a unit at p; multi-level atoms contain target primes up to e^{b₁Q},
   squared in Γ = ∏Z_i⁺Z_i⁻, so the lv6 choice can fail; proof critic issue 2). On a trivial gcd discard the rows and collect fresh ones (same field).

**Lemma RW (row identity for multi-level trees) [proved].** Let T := f_d (written t in thmA §3), β := Tα,
z_P := Ta − bβ ∈ Z[β], so ψ(z_P) ≡ T·G(P)/g₁ (mod N). In this lemma read every "t" below as T. Give each tree node v an exponent ε_v ∈ Z: a root prime ℓ gets +v_ℓ(u₀) or
−v_ℓ(v₀); a child of v gets −ε_v times its valuation in v's *rational value* (n_v := |G(P_v)|/q for a rational
target, |G(P_v)| for an algebraic one) or in v's *element* z_{P_v} (algebraic children). Then
  r² ≡ q_row·ψ(Z_row) (mod N),  Z_row := ∏_v z_{P_v}^{ε_v},
with q_row = ±(g₁/t)^{E_row}·(a B-smooth rational number), E_row := Σ_v ε_v over the expanded nodes (its parity is
thmA's single parity column), and the ideal of Z_row supported on first-degree primes of norm ≤ B and on primes
dividing D (the latter handled by lv6's group H). Every algebraic child is a *simple* degree-one prime: p ∤ D gives
p ∤ [O_K : Z[α]], and with gcd(a, b) = 1, p ∤ b the only prime above p containing a − bα is (p, a·b^{−1}). *Proof.* Induction on the number of expanded nodes, with the invariant
"r² ≡ (product of open targets and closed factors) mod N". Expanding a rational target q with exponent ε uses
q = σn_v^{−1}G(P) ≡ σg₁t^{−1}ψ(z_P)/n_v, which removes q^ε, adds z_P^ε to Z and n_v^{−ε} to the rational side; the
algebraic primes of z_P^ε become open with exponents ε·v_𝔭(z_P), and the rational primes of n_v^{−ε} become open.
Expanding an algebraic target 𝔮 with exponent ε_w multiplies Z by z_{P′}^{ε_w} and the rational side by
(tG(P′)/g₁)^{−ε_w}; since P′ is good, v_𝔮(z_{P′}) = 1 (otherwise q would divide the cofactor, and q > e^{λz}), so
the ε_w = −ε_parent·v_𝔮(z_P) exactly cancels 𝔮 in the ideal of Z. A good point's children do not divide D, so all
non-factor-base primes are opened and later closed; the tree is finite (below), so at the end only factor-base and
D-primes remain. ∎ The exponents satisfy |ε_v| ≤ (C log N)^{depth+1} = e^{O((log log N)²)} (every valuation is
O(log N)). Negative exponents: write Z_i = A_i/B_i with A_i, B_i ∈ Z[β]; LV6 is run on A_iB_i = Z_iB_i², and on a
dependency S the accumulated ∏_{i∈S}B_i is used as an explicit square correction when forming Y (replacing Z_i by
A_iB_i *without* that correction would change the residue row). Characters are unaffected because character primes
divide no atom (Lemma X_ERH).

**Lemma TS (tree size) [proved].** In any field, a DESC tree has depth ≤ ⌈log(b₁Q/(bE))/log(1/λ)⌉ + 1 = O(log t)
(children of a target e^z are ≤ e^{λz}, and targets are > B) and each node has ≤ 2·(log|F(P)| + log|G(P)|)/log B
= O(t/E) children, so ≤ e^{O((log t)²)} = e^{o(E)} nodes. The box has ≤ (2e^{sE} + 1)² points. (codex (30).)

**Lemma T_multi (transcript invariance for adaptive multi-level trees) [proved].** Every exposed quantity of NEST
(field, R_i, (u₀, v₀), canonical factorisations returned by recursive calls, all box scans, all failed and
successful DESC calls and their order, field restarts, columns, characters, kernel vectors, Hensel prime, Y) is a
function of (R_i) and coins independent of the hidden roots: DESC is a deterministic function of (A, target) given
the canonical factorisations, and recursive calls use fresh coins (recurse Lemmas OS, T_rec). Hence, as in Theorem
DD, X/Y is uniform on μ₂(Z/N) given the transcript, and Pr[1 < gcd(X − Y, N) < N | 𝒯] = 1 − 2^{1−ω(N)} ≥ 1/2.
(codex (34)–(35) gives the same argument; the dependency must be chosen before the hidden roots are used, and after a
trivial gcd fresh rows are drawn, as in lv6 §3.7.) ∎

*Proof of Theorem NR.*
- **Correctness.** Every output is certified (AKS, gcds), so the algorithm is Las Vegas; correctness never uses AD.
- **Field restarts.** By AD a field is good with probability ≥ 1/2; in a good field DESC never returns FAIL (every
  target has a good point, and its children are again targets). A bad field either fails (restart) or produces valid
  rows. By Lemma TS the work of every row attempt is bounded by the same quantity in good and bad fields, so the
  expected number of fields is ≤ 2 and the expected cost is at most twice the good-field bound.
- **Split acceptance.** Let k := ⌈t/(2b₁Q)⌉ + 1, X := ⌊√(N/2)⌋/2, Y := ⌊X^{1/k}/2⌋ (so B < Y < 2Y ≤ e^{b₁Q}), and
  𝒫 := primes in (Y, 2Y] not dividing D_A (all but O(t/Q) of them). After trial division to B, 2^{ω(N)} ≤ e^{O(t/E)}
  while |𝒫| ≥ Y/(3 log Y) = e^{(1+o(1))Q·(b₁ − o(1))} ≥ 4k³2^ω. Lemma S″ (recurse §3.2; codex (26)) with 𝒰 := all
  products of k distinct primes of 𝒫 (≤ (2Y)^k ≤ X, 2X² < N) gives Pr[accept] ≥ |𝒰|²/(2φ(N)) ≥ (|𝒫|/k)^{2k}/(2N)
  ≥ exp(−2k log(6k log Y) − O(k)) = e^{−(1/b₁ + o(1))E}, using 2k log(6k log Y) = (D′/b₁)(log t)(1 + o(1)).
  (Lemma S″ is used without the lower end of its window, a harmless generalisation; codex nested_2 issue 11.)
  Each attempt costs two recursive calls on integers ≤ √(2N).
- **Descent.** Per accepted row: e^{o(E)} nodes (Lemma TS), each scanning ≤ (2e^{sE} + 1)² points, each point
  costing two recursive factorisations of integers ≤ max(|F(P)|, |G(P)|) ≤ N^{τ_F + o(1)}, since with
  ‖P‖_∞ ≤ 3e^{sE}q: log|F(P)| ≤ log(d + 1) + κQ + d(z + sE + log 3) ≤ δD′(b₁Q + sE) + O(Q) = τ_F·t + O(Q) and
  log|G(P)| ≤ κQ + z + sE + O(1) = O(Q).
- **Linear algebra.** #columns ≤ (d + 2)π(B) + O(d²t) = e^{(b+o(1))E}; row weight e^{o(E)}; Wiedemann/KS
  e^{(2b+o(1))E}; Lemma SQ's Hensel lifting to e^{(b+o(1))E} digits (exponents ≤ e^{O(log² t)}): e^{(b+o(1))E}·poly.
  Expected number of LV6 rounds ≤ 2 (Lemma T_multi).
- **Recursion.** With F(n) the worst expected cost of complete factorisation of integers ≤ n, the above gives
  F(N) ≤ e^{(2b+o(1))E} + e^{(b+1/b₁+o(1))E}F(√(2N)) + e^{(b+2s+o(1))E}F(N^{τ_F+o(1)}).
  Since Λ_{1/3}(n^τ) = (τ^{1/3} + o(1))Λ_{1/3}(n), strong induction (recurse §3.1; codex (32)) gives
  F(N) ≤ L_N[1/3, c] for every c > c_NR. ∎

**Remarks on Theorem NR.**
- **What AD replaces.** FM, LD*, Lemma 8.5′, Kern, Prop PZ and Theorem K (every target is good, so no density γ is
  needed and the 2^ω loss is absent by Lemma S″). The precision problem of recurse §3.3 disappears. ERH is used only in
  Lemma X_ERH (and could be replaced by LV §6).
- **AD is the standard heuristic (reduced-basis form).** Under the usual "values behave like random integers of
  their size" model, using the worst-case sizes log|F(P)| ≤ κ_FQ + d(z + sE) and log|G(P)/q|, log|G(P)| ≤ m₀Q + sE + O(z),
  every target (regular or not) has ≈ e^{(2s − φ/3 + o(1))E} good points with
  φ = (κ_F + m₀)/(λb) + δ/λ + δs/(λb) (the δ/λ term is δ/(2λ) for regular lattices, recurse §3.4 / codex (29)),
  so the failure probability exp(−e^{Ω(E)}) per target beats the union bound over e^{O(Q)} targets whenever
  2s > φ/3 (`barrier.ad_heuristic_ok`). Example: δ = 1/2, b₁ = 1 (τ_F = 1/2), λ = 1/2, b = 4, κ_F = 4, m₀ = 2 (design
  (a‴) with M ≈ H ≈ N^{1/(d+1)}, so log M ≈ Q/δ), s = 3/4 satisfies the heuristic and gives c_NR ≈ 26.7
  (unoptimised; `barrier.nr_constant`).
- **AD is exactly what §§3–5 say we cannot prove.** It is a pointwise-in-the-field statement about smooth values of
  one binary form of degree d ≍ D′ on special-q lattices at L-scale. Averaging over 𝔉 cannot certify it by moments
  (§6: Θ(1) relative variance across fields; overdetermination for joint statements).
- AD for *every* target is stronger than needed (only targets that occur matter), but any weaker endogenous version
  runs into Prop SI(i).

## 9. Experiment E1 (pre-registered `7c17ffb`, Addendum A `e432e59` before the main run)

Toy endogeneity test: d = 3, 44-bit balanced semiprimes (seeds 1–12), M = 2^{20}, B = 2^{10}, Y = 2^{14}, S = 32,
400 parent nodes per instance (rational ℓ ≈ 2^{22}, first box point with Y-smooth cofactors). D1 = one-step
descendability (box with both cofactors B-smooth). Matched uniform baseline per child (same type, norm within a
factor 1.25); independent field A′ for rational children. 145 s total.

| prediction | result | verdict |
|---|---|---|
| P1 algebraic children, same field: CI(endo − unif) ⊂ [−0.08, 0.08] | 14333 children; 0.2957 vs 0.2944; diff +0.0013, CI [−0.0082, +0.0118] | holds |
| P2 rational children, same field | 1383; 0.1837 vs 0.1540; +0.030, CI [+0.003, +0.057] | holds |
| P4 rational children, independent field A′ | 1383; 0.2198 vs 0.1757; +0.044, CI [+0.016, +0.073] | holds |
| P3 all-children-good / ∏ pooled per-type rate ∈ [0.75, 1.33] | 4571 nodes; 0.0462 vs 0.0270; ratio 1.71 | **fails** |

- **P3 failure, diagnosed (exploratory, not pre-registered).** The per-field matched-uniform algebraic rates range
  from 0.110 to 0.587 across the 12 fields (the endogenous rates from 0.108 to 0.612), so the pooled rate understates
  E[p^m] (Jensen). With each node's expectation computed from its own field's rates, the aggregate ratio is 1.05:
  consistent, in aggregate, with a product model within a field (one aggregate statistic does not establish
  independence). Between fields descendability varies by a factor ≈ 5. The parent bootstrap gives no uncertainty
  over the field law (12 field clusters).
- **P2/P4 positive bias.** Endogenous rational children are slightly *more* descendable than matched uniform ones, in
  the same field and in an independent field alike, so it is not endogeneity: children are skewed to the small end
  of the ±25% window (mean log p 8.115 vs 8.133 for their matched draws).
- **Descriptive.** Children per node: 1: 229, 2: 903, 3: 1671, 4: 1393, 5: 529, 6: 74, 7: 1 (d = 3 has 6 field
  parameters, so the toy is too small to show overdetermination). Dead fractions fall with S (seed 1, algebraic:
  0.76 / 0.605 / 0.45 at S = 16/32/64). 18 degenerate children.
- **Interpretation (fixed in advance).** No endogeneity penalty at toy size, as the heuristic predicts; this supports
  nothing rigorous. The large between-field variance is the empirical face of §6.

## 10. Literature (literature subagent and codex nested_1, primary sources where marked)

- **Factoring.** Rigorous record L[1/2, 1 + o(1)]: Lenstra–Pomerance, JAMS 5 (1992) 483–516 (read; Thm 6.1 p. 499
  read in the original scan). Earlier: Pomerance 1987 L[1/2, √2]; Vallée L[1/2, √(4/3)]. Nothing below 1/2,
  rigorous or GRH-conditional, 1990–2026 (searches). Lee–Venkatesan, JNT 187 (2018), arXiv 1805.08873 (read): Thm 2.1
  unconditional L[1/3, (64/9)^{1/3}] for a congruence x² ≡ y² that "may be trivial" (Remark 2.2); Thm 2.3
  nontriviality only for N = pq, p ≡ q ≡ 3 (mod 4), under their Conjecture 7.1. Their relation collection uses HECM
  (lv6 §6.5).
- **Discrete logs.** Rigorous L[1/2, √2] for F_p and F_{2^n} (Pomerance 1987). Kleinjung–Wesolowski, JAMS 35 (2022),
  arXiv 1906.10668: rigorous (pn)^{2log₂n+O(1)} in F_{p^n}, hence rigorous L(α) for p = L_{p^n}(α), α < 1/2 (their
  corollary); also Granger–Kleinjung–Zumbrägel (Trans. AMS 2018), Lido (arXiv 2206.10327). The mechanism (Frobenius
  model; each elimination step is a point count on an explicit curve, Weil bounds; smoothness of polynomials is exact
  counting) has no analogue for integers. No rigorous L(α < 1/2) in medium/large characteristic.
- **Class groups.** Hafner–McCurley 1989 (ERH, imaginary quadratic, L(1/2, √2)); dBPW 2025/26 (ERH, all fields,
  L_{|Δ|}(1/2)); Biasse (Math. Comp. 2014, arXiv 0912.1927; arXiv 1204.1292) L(1/3) for special families, heuristic;
  Biasse–Fieker (LMS JCM 2014) GRH + heuristics.
- **Jacobians.** Enge–Gaudry–Thomé (arXiv 0905.2177): L(1/3), heuristic (Heuristics 5–6 cover relation collection
  and the special-Q descent). Rigorous results are L(1/2)-type (Enge 2002, Enge–Gaudry 2002, Couveignes 2001, Hess
  2004) or fixed-genus powers (Diem, Math. Comp. 2011).
- **Smooth values / linear systems.** BBDT (Comment. Math. Helv. 2012): y ≥ x^{g−2+ε} for binary forms of degree g;
  Tenenbaum–Wu (2003): L-type smoothness only for binary *quadratic* forms; one-variable polynomials of degree ≥ 2:
  no L-type results (de la Bretèche–Drappeau). Harper (Compositio) a + b = c for y ≥ (log x)^K; Matthiesen–Wang
  (ANT 2025) Ψ > N/log⁸N; Lachand (Acta Arith. 2017) fixed u.
- **GRH and smooth numbers in progressions.** Pointwise asymptotics are known for q ≤ y^{4√e−δ} (Soundararajan 2008;
  Harper JNT 2012), on average up to x^{1/2+} (Harper; Drappeau; Pascadi). A pointwise GRH statement for q ≫ y was
  not found (the literature agent suggests one may follow from Harper's Prop. 1 + GRH; unchecked, not used).
- **Family-averaged GRH.** Pierce–Turnage-Butterbaugh–Wood (Invent. Math. 2020), Lemke Oliver–Thorner–Zaman (2024):
  fixed degree, fixed Galois group only.
- **Theorem AF is not new as an identity:** dBPW (46) already states the two-sided bracket [e^{−1/4}, e^{1/4}], and
  (48) writes the weight with the upper incomplete Gamma function (their Remark 7.2 notes the 1/3 can be pushed toward
  1). What we add is reading it as an expectation over v ~ Gamma(r*, 1), i.e. as a law of the cofactor norm, and the
  consequence for NFS fields. Earlier one-sided relatives: de Boer–Ducas–Pellet-Mary–Wesolowski (CRYPTO 2020), de
  Boer's thesis (2022).
- Prop SI: no prior statement found (it is elementary and may be folklore); Prop SI(i) is also in codex nested_1
  ((19)), found independently.
- Theorem NR: no prior statement with a descent hypothesis found. Closest analogues: LV Thm 2.3 (a hypothesis ⇒
  L[1/3] factoring, for p ≡ q ≡ 3 mod 4 only), EGT's main theorem (Heuristics 5–6, including descent, ⇒ L(1/3) in
  Jacobians), Buhler–Lenstra–Pomerance 1993 (L[1/3] under explicit heuristics). An ETNFS record paper (arXiv
  2212.04999) calls termination of the descent an open problem.
- GKZ (Trans. AMS 2018) is rigorous only for fields admitting a suitable model; the KW L(α) consequence is stated in
  their introduction, not as a numbered corollary. de la Bretèche–Drappeau: JEMS 2020. Tenenbaum–Wu treat X² + Y²;
  the extension to all binary quadratic forms is BBDT's remark.

## 11. Answers to the task, item by item

1. dBPW: §1. Two-sided values with a shared rational side are outside its scope; its level structure is a single
   sampling step; its cofactor scale is ≥ √|Δ_K| (Theorem AF), typically N^{Θ(1)} for our fields (log|Δ_K| = O(t) is proved; a matching lower bound is not).
2. Designs: (a) ideal-class descent on the algebraic side with a randomised rational side: the rational side is not
   independent of the field once ψ is fixed (g₀ ≡ mg₁ mod N; codex (20)), and the algebraic step is dBPW-type
   (Theorem AF); (b) inductive good sets: well defined (codex (17)) but their density is exactly the open statement,
   and (log N)^{−C} per level gives only L[1/3 + o(1)] in the naive backtracking model (Prop BT); (c) dBPW at each
   level: no contraction (Cor AF, heuristic step); (d) fresh fields: rational children yes (Prop R), algebraic
   children only within their algebraic component (Lemma RW; constrained resampling, Lemma FB), and designs without
   algebraic children stay at 2/5 in the model (Prop R2).
3. Costs: tree e^{O(log² t)} nodes, recursion closes (Theorem NR's proof); the size model is right, the supply is
   not proved.
4. (T): Lemma T_multi, proved for adaptive multi-level trees with restarts.
5. ERH: needed for LV6 characters at growing degree (Lemma X_ERH), avoidable via LV §6; not the missing input for
   nested descent; no family-averaged replacement exists for growing degree.
6. Literature: §10. L[1/3] under ERH for factoring is not known.

**Best α under ERH: 1/2.** (Unconditional L[1/2, 1] by LP92; internally L[1/2, o(1)] by RA + P, unrefereed.)
**Sharpest barrier:** a quantitative supply of short two-sided special-q relations whose children lie in the
recursively defined good sets *of the same field* (the same algebraic component). Hypothesis AD is a sufficient
strengthening of it (Theorem NR). Density-only lemmas cannot give it for a fixed field (Prop SI); independence-over-
the-field moment methods meet the parameter budget (Lemma OD, typical-factorisation model); Arakelov-uniform sampling
works at the wrong norm scale (Theorem AF, heuristic step); and ERH, as far as we can see, does not supply it. None of
this is an impossibility theorem.

## 12. Reviews

All valid findings were applied in the text above.

**GPT-6.1-sol blind round `codex/nested_1`** (xhigh; report copied to `notes/nested_codex1_report.md`).
- No ERH-only α < 1/2. Its sharpest barrier is the same as ours: the supply of short two-sided relations whose
  children lie in completed-certificate good sets of the same field.
- It proved, independently of us:
  - the adaptive one-sided dBPW lemma (§7);
  - the arbitrary-prime-ideal monoid density bound;
  - the sparse-family counterexample (our Prop SI(i));
  - uniform split counting (= Lemma S″);
  - the completed-tree size bound;
  - the recursion lemma and (T) for adaptive trees.
- It found the dBPW order variant (§9.5), so computing O_K is not an obstruction (our first idea was that it was).
- It found the one-sided ray-condition rows, and the class obstruction.

**Opus proof critic (read-only, independent).**
- Verdicts:
  - Theorem AF: CORRECT.
  - SI, FB, OD (linear algebra), BT arithmetic, the recurse §3.5 correction: CORRECT.
  - X_ERH: CORRECT WITH FIXES.
  - Theorem NR: implication CORRECT WITH FIXES, but **AD as first stated was defective**.
- Fixes applied:
  1. *MAJOR.* Sup-norm boxes contain a single primitive point for irregular lattices, so AD would fail under its own
     heuristic. Fixed by switching to reduced-basis boxes. This gives τ_F = δb₁ instead of δb₁/2, a new φ, and a new
     example with c_NR ≈ 26.7.
  2. *MAJOR.* The Hensel prime in (B, 2B] can divide Γ. Fixed by taking it among the character primes.
  3. *MAJOR.* The Verdict omitted LV6′.
  4. *MAJOR.* R2's premise ‖f‖ ≥ N^{1/(d+1)} is false. Rewritten with the two-constraint window argument.
  5. *MINOR.* Box count, the Chebotarev ratio 1/(4.4d), poly(t)-bit character primes, LV6′ cost e^{o(E)},
     children-count wording, FB scope, BT phrased as a naive-strategy bound, 𝔉 sampling cost, AF novelty, and the
     cumulative-density wording.

**Opus literature critic (read-only, independent).** Most attributions CONFIRMED. Fixes applied:
- AF is dBPW (46)+(48) reread, not a new identity.
- dBPW Lemma 14.1 is ERH (it uses Lagarias–Odlyzko).
- Maynard's condition is n ≥ 4k.
- LV Remark 6.7 already asserts Lemma X_ERH's conclusion and the Θ(δ) correction.
- The printed radius of dBPW Thm 1.1 is inconsistent with Thm 9.5's.
- Closest analogues of Theorem NR: LV Thm 2.3, EGT, BLP 1993.
- Added references: ETNFS (arXiv 2212.04999), Barbulescu–Lachand, GKZ's model restriction.

**GPT-6.1-sol adversarial round `codex/nested_2`** (xhigh; report copied to `notes/nested_codex2_report.md`).
- Recommendation: MAJOR REVISION; no FATAL issue. It reviewed both the first draft and the revision made after the
  critics.
- It confirms: AF's comparison; SI; FB's local lemma; the BT arithmetic; X_ERH after fixes; RW, TS and T_multi; the
  split, the sizes, τ_F = δb₁ and c_NR (it reproduces 26.66); and the E1 numbers exactly, including the honestly
  reported P3 failure.
- Fixes applied:
  1. AF is a comparison, not a cofactor law. In the uniform model the deficit is Gamma(r* + 1); `barrier.af_weight`
     was corrected accordingly.
  2. Corollary AF now has three explicit premises, including the x^d − 2 counterexample to any universal
     discriminant lower bound.
  3. **Constrained polynomial resampling F + G·L₀·H** (a new design), discussed in §4. It is limited by the
     per-component parameter budget, because Lemma RW closes algebraic targets only within the field that opened
     them.
  4. R2 is restated as a model computation with premises (M1)–(M4). The "at most U pieces" step is replaced by a
     cost argument for jumping the window.
  5. The children-count lower bound was a wrong deduction. It is now a typical-factorisation model with the thin-tree
     caveat.
  6. The BT necessity language is removed, and the union-bound paragraph is narrowed (relative variance ≠ dead
     probability).
  7. X_ERH: GRH for ζ_K stated; d·ω(D) rejected ideals.
  8. NR: an admissible-field sampler, with admissibility part of the AD event; a Hensel prime reused from the
     character primes; T := f_d in RW; the exponent bound (C log N)^{depth+1}; the explicit denominator-square
     correction; the split constant 6.
  9. E1: rate labels (matched-uniform 0.110–0.587), the independence wording, and a field-cluster caveat.
  10. Verdict and Candidate card rescoped. AD is a sufficient strengthening, not an equivalent. "ERH does nothing"
      became "ERH does not supply". There are no impossibility claims.
- Rejected: none.

## 13. Candidate card

- **Mechanism.** ERH-conditional nested descent for the recursive Dixon-descent factoring algorithm, via dBPW sampling,
  fresh fields, or inductive good sets.
- **Model.** Theorem A / recurse architecture, one field of degree δD′, O(log log N) levels.
- **Result.** No rigorous α < 1/2 under ERH. Proved: Theorem AF (dBPW's success probability ≈ density among ideals of
  norm ≲ √|Δ_K|, a rereading of their (46)–(48)), Prop SI (sparse-image dilemma), Lemma FB, Lemma OD, Lemma X_ERH,
  Lemmas RW/TS/T_multi, and Theorem NR (ERH + AD + LV6′ + admissible sampler ⇒ L[1/3, c_NR]). Model computations:
  Prop R2, children counts, Prop BT.
- **Failure scope.** The known proof methods for designs whose algebraic components exceed the parameter budget
  (in the standard size model, every design below 2/5); the dBPW sampler under the smooth-ideal heuristic;
  independence-over-the-field moment methods beyond depth 1. Not excluded: conditional or local-factor moments,
  constrained resampling for small components, global recombination, higher-dimensional searches.
- **Next step.** Either a pointwise (fixed-field) L-scale smoothness theorem for binary forms on special-q lattices,
  or a non-descent source of hidden randomness at L[1/3]. Neither is in sight.
