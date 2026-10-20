# escapeR2: escaping the algebraic-children barrier below L[2/5] — results

Agent `escapeR2`, phase 13, 2026-10-20.

Files:
- This note. Pre-registration `notes/escapeR2_prereg.md` (commit `4624c3b`; code frozen at `42312c2`).
- Code: `src/escapeR2/transfer.py` (E1, parity-transfer spaces), `sizemodel.py` (E2, exact-rational LP for the size
  model), `resplit.py` (E3, harness), `run_experiments.py`. Tests: `tests/test_escapeR2.py`.
- Data: `experiments/raw/escapeR2_e{1,2,3}.json`.
- Codex: `codex/escapeR2_1.{md,out}` (blind round on the whole problem, before my draft), `codex/escapeR2_2.*`
  (adversarial referee, §12). The orchestrator's wild round `codex/r20_wild.out` items #1, #4, #7–#9 are evaluated in
  §§5, 7.

Labels: **proved** (a complete argument here, from published theorems or the project's internal results as cited),
**proved conditional on …**, **model** (exponent arithmetic under the stated premises, as in recurse §3.3 and nested
Prop R2), **heuristic**, **experimentally checked** (exact computation, toy), **open**. Everything about the project's
own algorithm is internal and unrefereed.

---

## Verdict

**Plain language.** None of the proposed escapes gives a rigorous factoring algorithm below L[2/5], and neither does
ERH. We tried norm export, fresh shared-root fields, class-group principalisation, Galois and Kummer transport,
re-splitting and factor-then-descend. In every one we inspected, a specific step fails, and the failure is now
written down precisely. Three kinds of finding come out of it:
- **Proved obstructions.** In the standard certificate form (per-field squares), an algebraic child can only be
  removed by elements of its own field, unless the fields are tied together by a common subfield structure (Lemma TI).
  The Kummer/odd-isogeny transport of r20 #4 lands back in the original field for the generic Galois group (Prop KT).
  The counting lemma behind the Dixon split does not extend to arbitrary square cosets (Prop RS).
- **A sharp model picture.** Inside one field, descent below 2/5 needs joint control of more values of one polynomial
  than it has coefficients. Components of O(d) points give nothing (Prop R3). In the chain sub-model, d^{C} points with
  any C > 1 give an exponent below 2/5: 3/8 at C = 3/2, tending to 1/3 as C → ∞.
- **A weaker sufficient hypothesis.** Pooling fields (Theorem NR-q) reduces "rigorous L[α]" to a per-target,
  field-averaged success probability of only L^{−O(1)}. It must hold on a polylog-dense set of split primes.

This is a catalogue of failed implementations plus one precise sufficient target. It is **not** a proof that every
route reduces to that target: conditional randomness in tied fields, directed recombination with non-diffuse labels,
and structured image certificates remain separate open alternatives (§9; codex escapeR2_2 §9).

**Precisely.**
1. **Lemma TI (cross-field parity transfer) [proved].** Let K_1, …, K_J lie in a Galois field L, and let ∏Γ_j be a square
   in L with Γ_j ∈ K_j. At every prime ℓ unramified in L, the parity vector of Γ_1 at the primes of K_1 above ℓ is
   pulled back from k_1 := K_1 ∩ K_2⋯K_J. **Corollary EX:** if k_1 = Q (e.g. K_1 primitive and not contained in the
   compositum of the others), then other fields, rational integers and norms cannot remove an algebraic child of a
   degree-one element of K_1 at such an ℓ. The certificate must contain further elements of K_1 in primes above ℓ. This
   blocks families 1(a), 1(b) and the conjugate/splitting-type versions of 3 in their stated forms (§§2–3, 5). Children
   ℓ dividing another field's discriminant are excluded from the statement; they are O(J·d·t/Λ) primes (§2.1). E1
   checks TI exactly on 72 824 group configurations: 0 violations.
2. **Tied fields do transfer parity [proved example].** A prime of Q(√ab) is ≡ mod 2 a sum of primes of Q(√a) and Q(√b)
   (E1: τ = 1). Whether such transfer can be made *useful* is open, i.e. realised by short elements with a compatible
   evaluation map and with fresh randomness (b can be fresh given a). We do not exclude it (§2.3).
3. **Kummer/odd-isogeny transport (r20 #4) [proved, for f of odd degree d ≥ 7 with Gal(f) = S_d and deg f′ ∈ {d, d+1}].**
   If J′[2] ≅ J[2] for J = Jac(y² = f), J′ = Jac(y² = f′), then the étale algebra of f′ is K or K × Q (Prop KT; Fact (F)
   proved from Bochert's theorem). So 2-torsion-preserving transport through hyperelliptic models returns to the same
   field. KT does not rule out other cohomological representations, other Galois groups, or even d (§5.2).
4. **Principalisation and class groups (families 1(c), 2) [proved reduction + model].** Removing 𝔭 modulo squares
   requires θ with (θ) = 𝔭·(factor base)·(square) (codex escapeR2_1 (5.1)). Under ERH for the Hecke L-functions of K, θ exists once
   the factor base contains all prime ideals of norm < 12 log²|Δ_K| (Bach). The only general rigorous way to *find* θ is at the
   Arakelov scale ≥ √|Δ_K| (dBPW; nested Thm AF). ψ(θ) is then a generic residue, whose only treatment we know is a
   re-split modulo N. A re-split puts the targets back at the top scale (§4, §6). For fields with explicit ψ, |Δ_K| = N^{Θ(1)}
   heuristically, so this route is L[1/2]. A structured way of producing θ with small evaluated image is a separate
   open alternative.
5. **Re-splitting (family 4(b,c)) [proved obstruction + proved bound for uniform masks].**
   - The square-corrected re-split w ≡ u/(v s²) is valid.
   - Lemma S″ does not extend from the square class to a coset wQ_N: pigeonholing the window primes by square class
     gives target sets with zero pairs (Prop RS(ii); E3).
   - For uniform masks, Pr[max(|u|, |v|) ≤ H] ≤ (2H + 1)²/|Q_N| (codex escapeR2_2), so masked re-splits do not contract.
   - Recursion supplies factorisations and tests. It is not itself a modular bridge (§6).
6. **Budget-respecting in-field descent stays at 2/5 [model] (Prop R3).** Under nested's premises M1–M4 and the
   typical-factorisation model, the floor is 2/5 when every algebraic component has O(d) points. O(d) is what the known
   moment methods control. This extends Prop R2 from "no algebraic children" to "budget-respecting algebraic children"
   (§7.1).
7. **Superlinear budgets in the chain sub-model [model; E2 exact].** Allow components of up to d^{C} points. The
   sub-model floor is α*(C) = min_{r ≤ C} max((r+2)/(3r+5), [C < r+1]·(r+2−C)/(2r+5−2C)):
   - (3 − C)/(7 − 2C) on [1, 3/2], so it falls below 2/5 for every C > 1;
   - 3/8 on [3/2, 5/2], 7/19 at C = 13/5, 4/11 on [8/3, 11/3], 5/14 at C = 4, 11/32 at C = 10;
   - plateaus (r+2)/(3r+5) on [r+1 − 1/(r+1), r+2 − 1/(r+1)], so it is above 1/3 for every finite C.

   This is an *upper bound* on the general model floor: one field, B-smooth rational sides, uniform child sizes.
   Prop R3 is the only general lower bound (C = 1 + o(1)) (§7.2). Codex escapeR2_1 found the C = ∞ level thresholds
   independently.
8. **Theorem NR-q [proved conditional on Hyp T_q, Hyp LV6′ at degree D^δ, character columns via ERH or the
   factorial-bit repair, and lv5unif's multi-field post-processing (outline)].** Suppose that for a set of split primes q
   of density ≥ (log N)^{−C₀}, Pr_A[the in-field search for q completes] ≥ L[α, −c]. Then factoring runs in L[α, c′].
   - The split uses a padded modulus N_* = 9^{k+1}N; codex escapeR2_2 found and repaired a window error in the draft.
   - This is weaker than nested's AD: averaged over fields, pointwise in the target, with only L^{−O(1)} success.
   - Consequence for r20 #1: its two-level moment ratio is needed only up to L^{O(1)} **fieldwise, for each q in a
     polylog-dense set**, rather than up to (log N)^{O(1)}. Jointly averaged moments do not suffice (codex
     escapeR2_2 §8 counterexample). The target density itself is still needed at polylog precision (§§7.3, 8).
9. **ERH.** ERH supplies none of the missing inputs. It removes character costs (Lemma X_ERH) and gives Bach's bound,
   but supplies neither in-field joint smoothness nor small evaluated images.
10. **Literature (§11).** We found no rigorous or ERH-conditional L[α < 1/2] factoring algorithm, no rigorous L(1/3)
    class-group or principal-ideal algorithm, and no smooth-number result for overdetermined linear systems at L-scale.
    Lemma TI is elementary and close to the norm-relation and Brauer-relation literature. We make no novelty claim.

### Per-approach exponent table

α is the size-model exponent the approach would give if its missing lemma were proved. "Blocked" means blocked in the
certificate class of §1 (per-field square certificates, or squares in a common Galois field with a supplied
compatible evaluation), not for every conceivable algorithm.

| # | approach | α in the model | status | missing lemma |
|---|---|---|---|---|
| 1a | replace 𝔭 by its norm ℓ (Galois conjugates / closure) | none: complement ≥ as hard | **blocked** (Cor EX, codex §3) | — |
| 1b | second field K′ sharing m; re-descend ℓ there | 2/5 (it is Prop MF) | cross-field closing **blocked** when K ∩ K′ = Q (TI) | — |
| 1b′ | tied fields in one Galois L (norm relations) | open | transfer exists (example TB); usefulness unproved | short transfer elements with compatible ψ and fresh randomness |
| 1c | ψ-image of a generator / S-unit | ≥ 2/5 if treated by re-split | model | contracting image certificate (codex (5.2)) |
| 2 | class-group descent (dBPW, Biasse–Fieker) | 1/2 rigorous (ERH); 1/3 only under descent heuristics | KNOWN + model | rigorous L_{|Δ|}(1/3) PIP with controlled ψ-images, or small-|Δ| fields with explicit ψ |
| 3 | Galois/Chebotarev; Kummer transport (r20 #4); Gassmann (r20 #7) | no gain in the inspected forms | blocked for 2-torsion-preserving transport, S_d, odd d ≥ 7 (KT); r20 #7 kill | a transport not induced by 2-torsion, or a non-S_d family with useful transfer |
| 4a | independent one-level descents; diffuse recombination (r20 #9) | 2/5 (R2, R3) | model; r20 #9 killed under diffuse labels | directed recombination with seeded rank surplus (open) |
| 4b | split law with halves re-split | no contraction for uniform masks | S″ fails on cosets (RS(ii)); uniform-mask bound (RS(iii)) | coset equidistribution + a contracting bridge |
| 4c | factor-then-descend | needs a modular bridge | structural | a costed modular bridge |
| 5a | in-field multi-level, components O(d) points | 2/5 | model (Prop R3) | — |
| 5b | in-field multi-level, components d^{C} points | chain sub-model: 7/18 at C = 5/4, 3/8 at C = 3/2, → 1/3 | model (E2) | Hyp LS_C (targetwise moments + search) |
| 5c | r20 #1: moments of complete two-level certificates | 3/8 | open; ratio needed only to L^{O(1)} fieldwise, targetwise (§7.3) | Hyp T_q at 3/8 |
| 5d | full multi-level descent | 1/3 | conditional (NR: AD; NR-q: T_q) | Hyp T_q at 1/3 (with a search analysis for growing depth) |

---

## 1. Setting and the certificate class

Notation as in nested §0 and lv5unif §1: t = log N, D = t/log t, Λ_α = D^α log t, L[α, c] = e^{(c+o(1))Λ_α}. A *field*
is A = (f, g) with F(g₀, g₁) ≡ 0 (mod N), K = Q(α), ψ: Z[α] → Z/N, α ↦ m = g₀/g₁. A *row* is
r² ≡ a·∏_j ψ_j(Z_j) (mod N) with a ∈ Q^× and Z_j ∈ K_j^× (codex escapeR2_1 (2.1)). A dependency makes a = a₀² and
Z_j ∈ K_j^{×2} for every j. The congruence is then completed by codex's square-root identity (2.5) with
denominator-square corrections. That identity is a direct check and is correct; it is lv5unif §11's multi-field
post-processing written as a formula.

**A more general multi-field square certificate** asks instead that ∏_j Γ_j be a square in a Galois field L
containing all K_j. To use it one also needs a ring homomorphism ψ_L on a suitable order of L, restricting to every
ψ_j, defined at the square root.
- *Such a ψ_L need not exist* (both reviewers). The map x_j ↦ m on Z[x_1,…,x_J]/(f_j(x_j)) descends to the order
  generated by the α_j only if the extra relations among the α_j vanish at (m, …, m).
- *Counterexample* (codex escapeR2_2): N = 119, α_1 = √2, α_2 = −√2, m = 11 (11² ≡ 2). Each field separately allows
  α_j ↦ 11, but α_1 + α_2 = 0 while 11 + 11 ≢ 0.
- In the V₄ triple Q(√ab), Q(√a), Q(√b), forcing all three generators to m forces m ≡ ±1 (proof critic).
- *When it does exist:* the compatibility is automatic, away from finitely many primes, when K_1 ⊗ … ⊗ K_J is a field
  (linearly disjoint fields).
- *How we use the class.* We treat it as a *superset* of the usable certificates: a ψ_L must be supplied, and the
  per-field certificate (Γ_j ∈ K_j^{×2}) is the special case actually used by NR-q and lv5unif. Lemma TI is a
  statement about valuations only, so it applies to the whole superset. Its obstructions therefore hold a fortiori.

An **algebraic child** of a degree-one element z_P = a − bα (gcd(a, b) = 1) is a prime 𝔭 = (ℓ, a/b mod ℓ) with
ℓ > B, ℓ ∤ b·f_d and ℓ ∤ [O_{K,(ℓ)} : Z_{(ℓ)}[α]] (localised index, so non-monic f is fine). To *remove* it means
to reach a certificate in which v_𝔭 is even.

---

## 2. Lemma TI: where cancellation can happen

### 2.1 Statement and proof

**Lemma TI [proved].** Let L/Q be Galois with group G, K_j = L^{H_j} (1 ≤ j ≤ J), and Γ_j ∈ K_j^× with ∏_j Γ_j ∈ L^{×2}.
Let ℓ be a rational prime unramified in L. For each j let x_j be the function 𝔭 ↦ v_𝔭(Γ_j) mod 2 on the primes of K_j
above ℓ. Put H′ := ∩_{j≥2} H_j (fixed field K_2⋯K_J), H₁⁺ := ⟨H_1, H′⟩ (fixed field k_1 := K_1 ∩ K_2⋯K_J). Then x_1
is constant on the fibres of the restriction map {primes of K_1 above ℓ} → {primes of k_1 above ℓ}.

*Proof.* Fix a prime 𝔓₀ of L above ℓ with decomposition group D. The primes of L above ℓ are σ𝔓₀ (σ ∈ G), and
σ𝔓₀ = τ𝔓₀ iff σD = τD. The prime of K_j below σ𝔓₀ is σ𝔓₀ ∩ K_j. Moreover σ𝔓₀ and τ𝔓₀ lie over the same prime of K_j
iff τ ∈ H_jσD, because H_j acts transitively on the primes of L above a given prime of K_j. Since ℓ is unramified in
L, v_{σ𝔓₀}(Γ_j) = v_{σ𝔓₀∩K_j}(Γ_j) for Γ_j ∈ K_j. Squareness in L gives, for every σ ∈ G,
  Σ_j x_j(σ𝔓₀ ∩ K_j) ≡ 0 (mod 2).
Write F_j(σ) := x_j(σ𝔓₀ ∩ K_j). Then F_j(hσ) = F_j(σ) for h ∈ H_j, since h fixes K_j. So F_1 = Σ_{j≥2}F_j is
invariant under left multiplication by H_1 and by H′, hence by H₁⁺. Two primes σ𝔓₀ ∩ K_1 and τ𝔓₀ ∩ K_1 lie over the
same prime of k_1 iff τ ∈ H₁⁺σD. So x_1 takes equal values on them. ∎

*Remarks.*
- *Excluded primes (codex escapeR2_2 #2, escapeR2_3 (e)).* TI needs ℓ unramified in L. Take L to be the compositum of
  the Galois closures of the K_j; that suffices for TI. Then the excluded ℓ divide some disc(f_j), j ≥ 1, and not only
  the source field's own exceptional primes. For a larger L, its extra ramification must be charged. There are ≤ Σ_j log|disc f_j|/log B =
  O(J·d·t/Λ) such ℓ > B. They are fixed once the fields are fixed. A descent can declare them non-targets, as nested
  excludes q | D_A, at the price of rejecting points whose values they divide. Parity transfer at ramified ℓ is not
  analysed here.
- The two-field case is sharp: Im_1 ∩ Im_2 consists of the functions left-invariant under ⟨H_1, H_2⟩, i.e. exactly the
  image of the ideals of K_1 ∩ K_2. This is a one-line consequence of the proof; E1 P3 confirms it, with τ₂ = 0 in all
  6 883 two-field cases.
- For three or more fields, the sum of the images can contain more than the separate pairwise intersections; that is
  why k_1 is defined with the compositum (E1(b)). The bound in TI is an inclusion, not an equality.

### 2.2 Corollary EX: export of algebraic children is impossible between unrelated fields

**Corollary EX [proved].** Assume k_1 = K_1 ∩ K_2⋯K_J = Q. Let Γ_1 = ∏ z_i^{e_i} be a product of degree-one elements
z_i = a_i − b_iα_1 (gcd(a_i, b_i) = 1). Let 𝔭 = (ℓ, r) be a prime of degree one with ℓ unramified in L and ℓ ∤ f_d. Since
ℓ ∤ f_d, α_1 is ℓ-integral. Assume ℓ ∤ [O_{K_1,(ℓ)} : Z_{(ℓ)}[α_1]] (localised index; this handles non-monic f), and
v_𝔭(Γ_1) odd. Then:
- (i) every prime of K_1 above ℓ has odd valuation in Γ_1;
- (ii) all primes of K_1 above ℓ have degree one, and each of them contains some z_i with e_i odd;
- (iii) in particular Γ_1 contains, besides the element that opened 𝔭, at least one further z_i lying in a prime above
  ℓ, unless K_1 = Q.

*Proof.* (i) is TI with k_1 = Q: there is one prime of Q above ℓ, so x_1 is constant, and it equals 1 at 𝔭.
(ii) If ℓ ∤ b_i, then z_i lies in at most one prime above ℓ, namely (ℓ, a_i/b_i), which has degree one. If ℓ | b_i,
then ℓ ∤ a_i and z_i is an ℓ-unit (its norm F(a_i, b_i)/f_d is ≡ a_i^d ≢ 0 mod ℓ). A prime 𝔓 of degree ≥ 2 above ℓ therefore has
v_𝔓(Γ_1) = 0, contradicting (i). So every prime above ℓ has degree one, and for each such prime the sum Σ e_i v(z_i)
is odd, so some z_i with e_i odd lies in it. (iii) If [K_1 : Q] = d > 1, then (ii) gives d ≥ 2 distinct primes
above ℓ, each containing an element of Γ_1. ∎

So norm export (replacing 𝔭 by ℓ and its conjugate complement), re-descending ℓ in a fresh field, and any Galois or
Chebotarev conditioning that keeps K_1 ∩ (other fields) = Q all leave the same obligation: elements of K_1 in primes
above ℓ. That obligation is an in-field descent of ℓ's primes. Its lattices {P : a ≡ r′b (mod ℓ)} are fixed by K_1, and
its values are values of the same form F_1. The norm route is in fact worse: it needs d such elements, one in each
prime above ℓ, instead of one in 𝔭 (codex escapeR2_1 §3 proves the d > 1 case of (i)–(ii) directly).

**When does k_1 = Q hold?** If f_1 has Galois group S_d with d ≥ 3, then K_1 is primitive (the point stabiliser is
maximal), so k_1 ∈ {Q, K_1}. Thus k_1 = Q unless K_1 ⊂ K_2⋯K_J. For the independent random fields of lv5unif's Design MF
(shared g, independent (w_j, c_j)), inclusion is a closed condition we do not expect. We do not prove a probability
bound. The inclusion K_1 ⊂ K_2⋯K_J is decidable, but we do not cost the test. Sharing the rational side does not, by
itself, create a common subfield.

### 2.3 Tied fields: parity transfer exists; usefulness is open

**Example TB [proved; E1(b)].** In L = Q(√a, √b) with G = V₄, let ℓ split completely, and label the four primes of L
by sign pairs (ε_a, ε_b) of √a, √b mod ℓ. A prime of Q(√a) is a set {(ε, ±)}, a prime of Q(√b) is a set {(±, ε′)},
and a prime of Q(√ab) is a set {(ε, ε′), (−ε, −ε′)}. Then {(ε, ε′), (ε, −ε′)} + {(ε, ε′), (−ε, ε′)} = {(ε, −ε′), (−ε, ε′)} mod 2.
So a prime of Q(√ab) is ≡ a sum of a prime of Q(√a) and a prime of Q(√b). Pairwise there is no transfer (τ₂ = 0);
jointly there is (τ = 1). This is the parity shadow of the norm relation (with denominator 2, in the sense of Biasse–Fieker–Hofmann–Page) 2 = (1 + σ) + (1 + τ) − σ(1 + στ)
in Q[V₄] (§11).

E1's census shows the same thing is common in small groups. Among 65 940 subgroup triples with D ranging over all
cyclic subgroups, τ > 0 in 40 149. In every one of them k_1 ≠ Q: among the 7 948 triples with k_1 = Q, τ > 0 in 0.

*What this does and does not show.*
- *First obstruction: an evaluation map.* Transfer needs a common Galois field L in which K_1 meets the compositum of
  the others non-trivially. Parity transfer by itself is not a certificate: it also needs a compatible ψ_L. In the
  V₄ triple, all generators mapping to the same m forces m ≡ ±1 (§1). With ψ(√a) = s_a and ψ(√b) = s_b chosen
  separately, each subfield needs its own small rational representation of s_a, s_b and s_a s_b for degree-one sieving.
- *Our first draft overstated the second obstruction.* It said such fields bring no independent randomness. That is
  false: fixing a and choosing a fresh b gives Q(√a), Q(√b), Q(√ab) with exactly the V₄ structure, while b remains free
  given a (codex escapeR2_2 §9).
- **[open]** Can parity transfer between tied fields be realised by short elements, with a compatible evaluation and
  enough conditional randomness for a moment argument? We do not exclude it, and we know no construction. The
  multiquadratic scale problem of §4 is one obstacle: discriminants ≥ √N per quadratic factor.

---

## 3. Family 1: exporting an algebraic child

**1(a) Norm export [refuted].** In parity coordinates at ℓ, multiplying by the rational integer ℓ adds the all-ones
vector (1, …, 1) over the primes of K above ℓ, while removing 𝔭 needs (1, 0, …, 0). For d > 1 these differ, and the
complement must be closed in K (Cor EX). Codex escapeR2_1 §3 adds three points:
- for odd d, the norm parity of the complement is d − 1 ≡ 0, so a rational norm check can falsely report success;
- N(z) ≡ (a − bm)·C(a, b, m) (mod N), so descending N(z) does not descend ψ(z);
- square norm does not imply square element (2 + √3).
With ramification ((ℓ) = 𝔭^d, d odd) the substitution would be valid. But ramified ℓ divide the discriminant, so
there are finitely many per field and they cannot carry endogenous children.

**1(b) Re-descend ℓ in a field K′ sharing the root m [refuted as an export; = Prop MF as a design].** A relation
ℓ ≡ b·ψ′(W) in K′ is valid and can be multiplied in. But TI (K ∩ K′ = Q) shows it does not change the parity of 𝔭
in K. The legitimate cross-field identity ψ(a − bα) = ψ′(a − bα′) moves a whole atom to another field, with its
entire factorisation (codex §4). Used at the parent's point, this is lv5unif's Prop MF: it replaces the parent's
algebraic value by another field's value at the same point. It gives the one-level 2/5 (lv5unif §11), not a closing.

**1(c) ψ-image of a generator or S-unit.** See §4: it reduces to class-group principalisation plus an image problem.

---

## 4. Family 2: class-group descent and principalisation

**Lemma PR (principalisation modulo squares) [proved; codex escapeR2_1 (5.1)].** An element θ with
(θ) = 𝔭·𝔞·𝔟² (𝔞 supported on the factor base 𝒮) exists iff [𝔭] ∈ 2Cl(K) + ⟨[𝔰] : 𝔰 ∈ 𝒮⟩. Replacing a row's Z by Zθ^{−1}
and a by a·ψ(θ) preserves the row and removes 𝔭 mod 2 (codex (2.4)).

*Under ERH (Hecke L-functions), existence is free.* Bach (Math. Comp. 55, 1990) shows that Cl(O_K) is generated by the prime
ideals of norm < 12 log²|Δ_K|. These are prime ideals of every residue degree, including those above index primes,
and the hypothesis is ERH for the Hecke L-functions of the class-group characters (i.e. RH for the Hilbert class
field's zeta function), not RH for ζ_K alone. Grenié–Molteni (Math. Comp. 87, 2018) give the degree-one, (4 + o(1))
log² version. The NFS factor base of nested Lemma RW has only degree-one primes of Z[α], so 𝒮 must be
enlarged to all prime ideals of O_K of norm ≤ 12 log²|Δ_K|. That is cheap, since |Δ_K| ≤ |disc f| and log|disc f| = O(d t).
Then every 𝔭 satisfies the condition of Lemma PR, and the question is never existence. It is (i) finding θ, and (ii) handling ψ(θ).

*(i) Finding θ.* The only rigorous general method is dBPW's sampler (ERH for Hecke L-functions). Their class-group
theorem costs poly·L_{|Δ_K|}(1/2)·L_{n^n}(2/3)·min(ρ_K, L_{|Δ_K|}(2/3 + o(1))); the L_{n^n}(2/3) factor is negligible
against L_N[1/2] for d = D^δ with δ < 3/4 (literature critic). The cofactor norm is ≥ √|Δ_K| (nested Theorem AF). For fields carrying an explicit ψ, the counting premise M1 gives a polynomial
discriminant ≈ H_f^{2d−2} ≈ N^{2(d−1)/(d+1)}. Whether the field discriminant matches is [heuristic], and nested's
Cor AF premise (i) records that no universal lower bound holds (x^d − 2). Under the smooth-ideal model this costs
L_N[1/2] (codex (6.1)). The L(1/3) class-group algorithms all rest on smoothness and relation-generation heuristics,
which are the descent heuristics again: Biasse (Math. Comp. 2014; arXiv 1204.1292) and Gélin–Joux (2016), for special
families. Biasse–Fieker (GRH + heuristics) gives L(2/3 + ε) in general, or L(1/2) for n ≤ (log|Δ|)^{3/4−ε}; the draft
wrongly grouped it with L(1/3) (§11).

*(ii) ψ(θ).* θ is a general element of the ideal lattice, not a degree-one element. ψ(θ) = Σ θ_i m^i mod N has no
small rational representation. It is a generic residue. The only known way to write it over smooth integers modulo N
is a re-split (§6). That returns targets of size ≥ L[1 − α] in an arbitrary square class (Prop RS), so this route
collapses to the top of the descent [model]. Small degree-one θ would be an in-field closing, i.e. the original
problem.

*Small-discriminant fields (codex row 4).* If fields with log|Δ_K| = O(Λ_{2/3}) and explicit ψ existed for every N,
then L_{|Δ|}(1/2) would be L_N[1/3]. But an explicit ψ on a field of small discriminant means an explicit root modulo N
of a polynomial of small height. For a fixed small field that is a root-finding problem modulo a composite. For
x² − c with fixed small c, finding a root modulo N is not known to be easier than factoring (no method below the NFS is
known). Rabin's equivalence with factoring is for roots of random residues. For fixed c, factoring follows only from two
roots that do not differ by sign. The NFS avoids this only by building f from N, which forces large height.
Example: quadratic fields Q(√(a² − kN)) have explicit ψ(√D) = a, but |D| ≳ √N unless one is lucky (probability
≈ |D|/√N per trial), so the Arakelov scale is ≥ N^{1/4} and the method is L[1/2] (it is MPQS-like). [heuristic.] Codex
found the same and leaves it as a 5% speculative row; we agree it is not excluded, but we know no construction.

**Dimension cost of putting large primes into LV6's residual group [proved; codex (8.1)].** dim_{F₂} H_S ≥ |S| + r₁ + r₂.
So absorbing private large primes into the exceptional support costs one kernel dimension each. There is no free
compression.

---

## 5. Family 3: Galois, Chebotarev, Kummer and Gassmann transport

### 5.1 Galois structure and splitting types

- A typical S_d field has Aut(K) = 1. Conjugate fields Q(α_i) are isomorphic *pointed* fields only if ψ is transported
  with them. Their evaluation maps send α_i to other roots of f modulo N, which have no small rational representation.
  With representation, they are not fresh samples (codex §7).
- If K is Galois and ℓ splits completely, (ℓ) = ∏_σ σ𝔭. Conjugation permutes coordinates and does not identify them
  (codex §7; TI with H_1 = 1).
- Conditioning on splitting types changes which primes appear, not where they must be closed. Cor EX(ii) shows that
  once a degree-one 𝔭 has odd exponent, ℓ must split into degree-one primes and every one of them needs an element of
  K_1. Chebotarev controls how often that happens; it supplies no element.
- **Odd-degree extensions are no square-class sink [proved; codex §7].** K^×/K^{×2} → L^×/L^{×2} is injective for odd
  [L : K]. Even-degree extensions need √ψ(z) modulo N, which is exactly the hidden information.

### 5.2 Kummer transport through odd isogenies (r20 #4)

r20 #4 proposes to transport the norm-kernel part of an algebraic square class through Kummer/2-descent
correspondences between Jac(y² = f) and Jac(y² = f′), using odd-degree isogenies to keep 2-torsion. Its own
adversarial verdict already noted that rational odd isogenies are rare and that isogenous varieties need not be
Jacobians. We add a structural kill for the generic case.

**Fact (F) [proved, from Bochert's theorem].** For d ≥ 7, every subgroup H of S_d of index ≤ d + 1 is S_d, A_d or a
point stabiliser. (The reviewers disagreed on whether Dixon–Mortimer Thm 5.2B covers d = 7, 8: codex read it as
d ≥ 7, while the literature critic quotes "n > 8" via Babai, arXiv 1512.03547. So we give a direct argument.)
- (i) *H intransitive.* Take an orbit of size k ≤ d/2. Then H ≤ S_k × S_{d−k} has index ≥ C(d, k), and C(d, k) > d + 1
  for 2 ≤ k ≤ d/2 (d ≥ 4). So k = 1, H ≤ S_{d−1}, and the index bound forces H = S_{d−1}.
- (ii) *H transitive and imprimitive,* with blocks of size b, 2 ≤ b ≤ d/2. The index is ≥ the number of such block
  systems, d!/(b!^{d/b}(d/b)!). This is ≥ 35 for d = 8 and > d + 1 for every composite d ≥ 8; d = 7 has no such H.
- (iii) *H primitive, not containing A_d.* Bochert's theorem gives index ≥ ⌊(d + 1)/2⌋! ≥ 24 > d + 1.
- (iv) Otherwise H ⊇ A_d.

The bound d ≥ 7 is sharp: PGL₂(F₅) has index 6 in S_5, and S_6 has a transitive S_5. For d ≥ 9 this is also
Dixon–Mortimer Thm 5.2B with r = 2.

**Proposition KT [proved; uses Fact (F), i.e. Bochert's theorem].** Let f ∈ Z[x] be squarefree of odd degree d ≥ 7 with Galois group S_d,
K = Q[x]/(f), J = Jac(y² = f). Let f′ be squarefree of degree d or d + 1, J′ = Jac(y² = f′) of the same genus, and
suppose J′[2] ≅ J[2] as Gal(Q̄/Q)-modules. This holds whenever J′ and J are related by an isogeny of odd degree
defined over Q. Then Q[x]/(f′) ≅ K if deg f′ = d, and ≅ K × Q if deg f′ = d + 1.

*Proof.* Let X, X′ be the root sets. For odd |X|, J[2] ≅ F₂[X]⁰ (sum-zero vectors), and for even |X′|,
J′[2] ≅ F₂[X′]⁰/⟨𝟙⟩ (Schaefer 1995 §2; Stoll 2001). An odd-degree isogeny has odd-order kernel, so it is injective on
2-torsion, hence an isomorphism J[2] → J′[2] of Galois modules. The field Q(J[2]) is the splitting field L of f,
because S_X acts faithfully on F₂[X]⁰ for |X| ≥ 3. Likewise Q(J′[2]) is the splitting field of f′, because S_{X′}
acts faithfully on F₂[X′]⁰/⟨𝟙⟩ for |X′| ≥ 5 (repaired after both reviewers):
- take π ≠ 1, x with π(x) = y ≠ x, and z ∉ {x, y};
- then π(e_x + e_z) = e_y + e_{π(z)} is a weight-2 vector different from e_x + e_z, since {x, z} = {y, π(z)} would
  force z = y;
- it is also different from e_x + e_z + 𝟙, which has weight |X′| − 2 ≠ 2;
- so π acts non-trivially. The same argument without the 𝟙 case gives the odd case.

The module descriptions are standard: Schaefer (JNT 51, 1995) for odd degree; Stoll (Acta Arith. 98, 2001) and Poonen–Schaefer (Crelle 488, 1997) for even degree; see also Morgan (Proc. LMS 131, 2025, e70066). So both splitting fields equal L, and G = S_d acts faithfully on X′ with |X′| ≤ d + 1. By
(F), each orbit has size 1, 2 or d, since 2d > d + 1. If no orbit has size d, then A_d fixes X′ pointwise,
contradicting faithfulness. So there is exactly one orbit of size d, with point stabiliser S_{d−1}, and it is
isomorphic to X as a G-set (d ≠ 6, all index-d subgroups are conjugate). The remainder has size 0 or 1. The étale
algebra of f′ is therefore K, or K × Q. ∎

So any transport that preserves 2-torsion through hyperelliptic models lands in a conjugate of K itself.
- *A new representation of the same field* is a different generator α′ with ψ(α′) ≡ g′₀/g′₁.
- *Withdrawn.* Our draft claimed, by a heuristic count, that every such α′ has Arakelov-scale height. That is false
  (codex escapeR2_2 #6). α′ = α + 1 has H_{f(x−1)} ≤ (d+1)2^d H_f and ψ(α′) = (g₀ + g₁)/g₁.
- *Möbius changes.* For α′ = (pα + q)/(rα + s) with small p, q, r, s, a − bα′ = ((as − bq) − (bp − ar)α)/(rα + s). So
  degree-one elements in α′ are degree-one elements in α, divided by one fixed element, under a unimodular-up-to-scalar
  change of (a, b). They give no new values.
- *Other generators.* Generators that are not Möbius images of α give elements of degree ≥ 2 in α. We prove nothing about
  their heights.

*Scope (narrowed after codex escapeR2_2 #7).*
- KT identifies the étale algebra of any hyperelliptic model whose 2-torsion is isomorphic to J[2], for Gal(f) = S_d,
  d odd, d ≥ 7, and deg f′ ∈ {d, d + 1}. That kills r20 #4's odd-isogeny formulation for generic NFS fields.
- It does not classify all useful representations or transports. Excluded from the statement:
  - other Galois groups, d < 7, even d;
  - cohomological transports not induced by a 2-torsion isomorphism;
  - non-hyperelliptic carriers.
- For exotic groups (Gassmann-type), r20 #7's kill covers the transports that send embedding vectors to embedding
  vectors. Mixing transports there land in tied fields (§2.3), whose usefulness is open.

### 5.3 Gassmann / torus transport (r20 #7)

We checked r20's kill and agree. Over a common Galois closure the character lattices are permutation lattices, and an
odd-degree isogeny is invertible mod 2. A basis vector mapping to a basis vector propagates along the orbit to a
G-set bijection, hence isomorphic fields. This is a special case of the phenomenon TI captures (transfer only through
tied structure).

---

## 6. Family 4: re-split, factor-then-descend, independent descents

**Prop RS (re-splitting modulo N).**
- (i) [proved; codex (10.1)] For a unit w and fresh s, R := ws² and a Gauss-reduced u ≡ Rv give w ≡ u/(v s²). The
  original target is preserved with an explicit square correction. Non-units give gcd exits.
- (ii) [proved; construction after codex escapeR2_2 #8] Lemma S″'s bound Pr ≥ |S|²/(2φ(N)) holds for R uniform on the
  squares Q_N. For R uniform on a coset wQ_N the relevant count is Σ_c n_c n_{wc} (codex (10.2)), and it can vanish.
  - *Construction:* partition the primes of any window (coprime to N) into their h = 2^{ω(N)} square classes, and let G
    be a largest class. Then |G| ≥ #window/h, every k-product from G lies in one class, so every ratio of two lies in
    Q_N.
  - For w ∉ Q_N ∪ −Q_N (such w exist for ω ≥ 2), no signed pair exists. For large windows, |G| ≥ 4k³h and S″'s product
    window are available too, so this is a full-hypothesis counterexample. No prime-distribution theorem is needed.
  - *E3* (5 seeds, 48-bit N; an illustration with smaller adversarial G, 59–70 primes < 128) found 0 positive pairs for
    non-square w, against 263 040–640 312 for w = 1 on the same S.
  - *Disclosure (codex escapeR2_2 #9).* E3 drew w only with (w/p₁) = −1. In seeds 2, 3 and 5 the class of w is (−1, 1).
    For seed 3, whose primes are 3 and 1 mod 4, that is the class of −1. There signed pairs exist and the
    absolute-value hit probability is positive, although 0 hits occurred in 2 000 draws. The positive-pair count of R1
    is unaffected; the intended obstruction needs w ∉ ±Q_N.
  - The control G = all window primes gives 4.02·10⁷ pairs ≥ |S|²/(2h) = 2.06·10⁷.
- (iii) [proved for uniform masks; codex escapeR2_2 §4] Let R be uniform on wQ_N. The pairs (u, v) with |u|, |v| ≤ H
  account for at most (2H + 1)² residue ratios, so Pr[max(|u|, |v|) ≤ H] ≤ (2H + 1)²/|Q_N| ≤ (2H + 1)²·4^{ω}/N, uniformly
  in w. For H = N^{1/2−ε} this is N^{−2ε+o(1)} after preprocessing (ω = o(log N)). So L[α, O(1)] masked trials do not
  contract the value size. *[model]* Translating "values ≈ √N" into "targets ≥ L[1 − α]" uses the smoothness model.
  Deliberately constructed identities (e.g. s = 1 and small w) are not covered by this sampling-law statement.
- *Missing lemma:* class equidistribution of k-products of the descendable primes in arbitrary cosets of Q_N, plus a
  contracting bridge not of uniform-mask type. (i) and (ii) were found independently by codex.

**What recursion does and does not supply [structural; reworded after both reviewers].**
- The recursive call returns the complete factorisation of a public integer M < N.
- A congruence modulo N needs a multiplicative relation. Factorisations enter such relations only through a separate
  modular bridge, i.e. an identity x ≡ a·M, x ≡ M/s², x ≡ u/v, … (mod N) that is constructed and costed independently.
- The bridges used in this project are special-q lattice values (M = G(P) ≡ g₁ψ(z_P)) and Gauss splits / re-splits
  (values ≈ √N; RS(iii)). Recursion adds no new bridge; it makes testing cheap.
- We do not claim these are the only possible bridges. Codex §10 also notes that a square congruence modulo M says
  nothing modulo N.

**Independent one-level descents [model].** If every descent closes fully at one level, combining them changes
constants only (R2, R3). Global recombination of partial rows needs repeated large labels. Under diffuse labels the
expected number of repeats is ≪ M²/V = o(1) (codex §8, r20 #9). We agree with r20's kill under diffuse statistics.
*Directed* recombination is not covered: rows deliberately constrained to share large ideals, with a seeded rank-surplus
theorem for the resulting graph. It is a separate open alternative, and cyclic components escape R3's unrolling
(codex escapeR2_2 §9).

---

## 7. In-field descent and the component budget

### 7.1 Prop R3: budget-respecting components stay at 2/5

**Premises.** nested Prop R2's (M1)–(M4), together with the typical-factorisation model of nested §5: the algebraic value
of a node, of log-size X, is the sum of the log-sizes of its children plus a B-smooth part of log-size ≤ ΛU. Here
Λ = log B, U is the largest u with u log u ≤ Λ, and a node costs L[α] only if each smoothness event has u ≤ U. Nodes may
use their own fields (any d, H_f, H_g); algebraic children are closed in the parent's field. An **algebraic
component** is the set of nodes linked by algebraic edges, all in one field. *Budget:* every component has at most
κ(d − 2) points, for a constant κ ≥ 2. κ = 2 covers the d + 1 points that the known moment methods control (nested
Lemma OD), for d ≥ 5. Components are trees: one closing point per child, as in RW. A point that closes several children
may also be shared (an acyclic component, i.e. a DAG). Then the path-weighted unrolling below double-counts shared
descendants. Instead, order the n ≤ κ(d − 2) distinct nodes topologically and apply the same edge recurrence node by
node: X_v ≤ (d/(d − 2))(1 + 2/(d − 2))^{n−1}ΛU ≤ 2e^{2κ}ΛU for d ≥ 4 (codex escapeR2_3 (b)). This is again O_κ(ΛU), and
the conclusion below is unchanged. Cyclic components are not covered.

**Prop R3 [model computation].** Under these premises, α ≥ 2/5.

*Proof.* Claim: every node v satisfies X_v := log H_f + (d/2)z_v ≤ 2(1 + 2κ)ΛU for d ≥ 4, where z_v is the target's
log-size. The value of v consists of: the target itself (only if v is an algebraic node, log-size z_v ≤ 2X_v/d), the
children (log-sizes y_w), and a B-smooth part (≤ ΛU in the model). So X_v ≤ z_v + ΛU + Σ_w y_w, hence
X_v ≤ c(ΛU + Σ_w y_w) with c := d/(d − 2). (This corrects the draft, which dropped z_v. Codex escapeR2_2 #11 showed the
draft's constant 1 + 2κ fails: a d = 5 chain has full norm logarithms 5/3, 25/9, 95/27 > 3. The corrected constant
2(1 + 2κ) = 6 covers it.) The closing node w of a child of log-size y_w has X_w ≥ (d/2)y_w, so
y_w ≤ (2/d)X_w. Hence X_v ≤ cΛU + (2/(d − 2))Σ_{w child of v} X_w. Unrolling over the component below v gives
X_v ≤ cΛU·Σ_u (2/(d − 2))^{depth(u)}. Every node u below v contributes ≤ 2/(d − 2), and there are ≤ κ(d − 2) of
them, so X_v ≤ cΛU(1 + 2κ) ≤ 2(1 + 2κ)ΛU. This is R2's "algebraic feasibility gives log H_f + (d/2)z = O(ΛU)", with
the constant 2(1 + 2κ) in place of 1. For d = 3 the weights are 2, and the bounded node count κ gives
X_v ≤ 3·κ2^{κ}ΛU = O_κ(ΛU) instead. R2's window argument then applies verbatim with A replaced by O_κ(A):
- if log H_g ≤ t/4, then log H_f ≥ t/(2(d + 1)), so by AM–GM X_v ≥ √(tz/2) and z ≤ 8(1 + 2κ)²Λ²U²/t ≍ A;
- if z ≫ A, the rational cofactor has log-size ≥ t/8, so a rational child has log-size ≥ W = t/(8U);
- children in the window (A, W) are infeasible, and jumping below e^A costs u ≥ t/(8A) ≫ U for α < 2/5;
- the top targets z₁ ≍ D^{1−α} log t lie above A ≍ D^{4α−1} log t exactly when α < 2/5.
So a finite tree of cost L[α] cannot exist for α < 2/5. ∎

So nested's Lemma OD and Prop R2 combine into one statement: **with any jointly controllable component size O(d),
algebraic children buy nothing.**

### 7.2 Superlinear budgets: the chain sub-model (E2)

Chain sub-model (`sizemodel.py`). One field, d = D^δ, base-m heights log H_f = log H_g = t/(d + 1). Root target
ζ₀ = 1 − α; algebraic levels ζ_1 ≥ … ≥ ζ_r; value v_j = max(1 − δ, δ + ζ_j) with v_j ≤ α + ζ_{j+1}; terminal
v_r ≤ 2α; B-smooth rational sides (1 − δ ≤ 2α, ζ_j ≤ 2α); component size Σ_j max(0, v_j − ζ_{j+1}) ≤ Cδ (exponents of D).

*Closed form.* When v_j = δ + ζ_j, the component exponent telescopes to rδ + ζ₀ − ζ_r. It is minimised at
ζ_r = 2α − δ. The chain needs δ ≤ ((r + 3)α − 1)/(r + 1), and the rational sides need δ ≥ 1 − 2α. This gives
α ≥ (r + 2)/(3r + 5) without a budget (codex escapeR2_1 §1: (r_c + 1)/(3r_c + 2) with r_c = r + 1). With a budget,
(r + 1 − C)δ ≤ 3α − 1, so α ≥ (r + 2 − C)/(2r + 5 − 2C) when r + 1 > C. Moreover fixed depth r ≥ 1 needs C ≥ r:
monotone targets give v_j − ζ_{j+1} ≥ δ, so the component exponent is ≥ rδ (codex escapeR2_2 #13). The LP enforces
this implicitly; the draft's helper `closed_budget` did not (it returned 3/7 at r = 2, C = 1, where the LP is
infeasible). It was fixed after review; the global minimum is unaffected. Hence
  α_r(C) = max((r+2)/(3r+5), (r+2−C)/(2r+5−2C)) for r ≤ C < r + 1;  (r+2)/(3r+5) for C ≥ r + 1;  infeasible for C < r,
  α*(C) = min_r α_r(C).
The plateau at value (r+2)/(3r+5) is exactly [r+1 − 1/(r+1), r+2 − 1/(r+1)] (both reviewers); level r is also optimal on the descending segment just before it. Between plateaus α* interpolates; for
example α*(13/5) = 7/19 ≈ 0.368 (proof critic).

**E2 (exact-rational simplex; all pre-registered predictions hold).**

| C (component ≤ d^C points) | α*(C) | levels r | closed form |
|---|---|---|---|
| 1 | **2/5** | 0 | 2/5 (= R3) |
| 5/4 | 7/18 ≈ 0.389 | 1 | (3 − C)/(7 − 2C) |
| 3/2, 2, 5/2 | **3/8** | 1 | 3/8 |
| 13/5 (post-review check) | 7/19 ≈ 0.368 | 2 | (r+2−C)/(2r+5−2C) |
| 8/3, 3 | 4/11 ≈ 0.364 | 2 | 4/11 |
| 4 | 5/14 ≈ 0.357 | 3 | 5/14 |
| 6 | 7/20 = 0.35 | 5 | 7/20 |
| 10 | 11/32 ≈ 0.344 | 9 | 11/32 |
| ∞ | → 1/3 | → ∞ | (r + 2)/(3r + 5) |

The LP agrees with the closed form to 10⁻⁹ at every grid point. Each closed form is feasible exactly and infeasible
10⁻⁶ below. At α = 3/8 (δ = 1/4, ζ_1 = 1/2) the first-level value has exponent v₀ = 7/8. Its children number
D^{3/8} = d^{3/2}: the "3/8 vs 1/4" count of codex §11 and r20 #1.

*Reading.*
- These are floors of the **chain sub-model**: one field, B-smooth rational sides at the algebraic levels, base-m
  heights, uniform child sizes, no sharing. Each is achievable in the size model, so it is an *upper bound* on the
  general model floor for budget d^C.
- The only general *lower bound* we have is Prop R3, for O(d)-point components.
- Within the sub-model, any fixed C > 1 gives an exponent below 2/5, and no finite C reaches 1/3.
- The *component size* is the D-exponent of the typical point count (distinct, typically sized children). It is a
  model count, not one forced by smoothness alone.

The in-field input this needs is:

**Hypothesis LS_C (joint smooth supply; stated as moments plus search) [open].** For a degree-d field law, and for each
target q in a set of split primes of density ≥ (log N)^{−C₀}, let Z_{A,q} count the complete certificates for q (of
≤ d^C points, in the prescribed boxes). Two conditions are needed:
- (a) *targetwise moments:* E_A[Z_{A,q}²] ≤ L^{c}·(E_A Z_{A,q})² with E_A Z_{A,q} > 0;
- (b) *search:* a deterministic search, of cost L[α, O(1)] in every field, finds a certificate whenever Z_{A,q} ≥ 1.
Codex escapeR2_2 #14: a lower bound on an averaged count alone controls neither overlaps nor search. (a) and (b)
together imply Hyp T_q (§8). (b) holds at bounded depth (§7.3).

For s > d + 1 the forms are linearly dependent (nested Lemma OD). This is a "linear equations in smooth numbers"
problem of growing complexity (s ≈ d^C forms in d variables), with coefficients of size up to exp(Λ_{7/8}), at
smoothness u ≍ D^α. The nearest known results are Harper's a + b = c (one equation, three variables, coefficients 1;
y ≥ (log x)^K; Compositio 2016) and Matthiesen–Wang (finite-complexity systems, near-full y range per their arXiv
abstract; u bounded or growing extremely slowly). Neither reaches it (§11). (The draft's "Ψ(x, y) ≥ x/(log x)^8", inherited
from nested §10, does not match their abstract and is dropped.)

### 7.3 r20 #1, evaluated adversarially

r20 #1 proposes first and second moments of Z_{A,q} = #complete two-level certificates for q, at α = 3/8, d ≍ D^{1/4},
with saddle-point-tilted weights and sparse transference, and requires polylogarithmic precision.

- **The exponent is right [model; E2].** Two levels at 3/8 need C = 3/2. Smaller budgets already give exponents
  below 2/5 in the chain sub-model (7/18 at C = 5/4), so 3/8 is not the "first" step. It is the first plateau, the
  least C at which a single algebraic level is unconstrained by the budget.
- **The moment ratio can be relaxed to L^{O(1)}, fieldwise and targetwise [proved conditional, §8].** With a pool of
  J_f = L^{O(1)} fields (lv5unif Prop MF's mechanism), it suffices that, for every q in a set of density
  ≥ (log N)^{−C₀}, E_A Z² ≤ L^{c}(E_A Z)². By Paley–Zygmund, Pr_A[Z ≥ 1] ≥ (E_A Z)²/E_A Z² ≥ L^{−c}. This needs no
  independence between certificates. Then:
  - The pool covers half of these q with probability ≥ 1/2.
  - The split sees a target set of polylog density, and the γ^{2k} amplification hits only that density, not the moment
    ratio.
  - *The quantifier matters (codex escapeR2_2 §8, proof critic #3).* Moments averaged jointly over (A, q) do not
    suffice. Take Z_{A,q} = 1 on a target subset of density ε = e^{−cΛ} and 0 elsewhere, for every A. Then
    E Z = ε and E Z²/(E Z)² = 1/ε: both relaxed joint inequalities hold. But the split pays ε^{2k} = exp(−Θ(kΛ)), outside
    L[α]. r20 #1 stated its moments with q uniform in the shell; the relaxation applies only to the targetwise
    version, and the target *density* is still needed at polylog precision.
- **Search at bounded depth [proved; codex escapeR2_2 §8].** A recursive feasibility search finds a certificate whenever
  one exists: for each candidate point, check each child independently. Its cost satisfies T_h ≤ W·(work per
  candidate + m·T_{h−1}), with W = L[α, O(1)] candidates per box and m = e^{o(Λ)} children per candidate. That is
  L[α, O(h)] at fixed depth h. (nested's Lemma TS bounds successful trees only; failed searches need this recurrence.)
- **Consequence for r20's kill test.** r20 would kill the transference route if its majorant lost exp(Ω(d)) in the
  relative moment ratio. After pooling, a targetwise relative loss of exp(O(s)) with s = d^{3/2} = D^{3/8} is
  exp(o(Λ_{3/8})), and even exp(O(s log log N)) is L[3/8, O(1)]. Neither is raised to the power 2k. So local-factor
  losses of e^{O(1)} per form, which a transference majorant typically incurs, no longer kill the route on that count.
- **What remains [open].**
  - A lower bound, within L^{O(1)} of the main term, for the number of c at which s ≈ d^{3/2} dependent affine forms are
    simultaneously smooth with prescribed divisibility.
  - A matching targetwise second moment.
  - The main term is tiny: ∏_i ρ_i ≈ exp(−s·u log u) = L[3/4]^{−Θ(1)} per configuration. A relative loss exp(O(s)) is
    affordable, but additive errors must be below that scale, which no existing circle-method or transference estimate
    gives at L-scale.
  - An averaged first moment alone also does not control rare fields carrying very many certificates; that is what the
    second moment is for.
- **Verdict.** It survives as a precise sufficient analytic target. We agree with r20's 3% and do not raise it: the
  relaxation removes one precision objection, not the core estimate.

---

## 8. Theorem NR-q: the weakest hypothesis we know that suffices

**Setting.** Fix α ∈ [1/3, 1/2), parameters as in nested §8 (α = 1/3) or the chain sub-model of §7.2, and a field law
𝔉, sampled in time e^{o(Λ_α)}. Its admissible fields have degree d ≍ D^δ, with δ < α, or α = 1/3 and nested's
d = δ′D^{1/3}. In the boundary case the character and LV6′ costs are L[1/3, Θ(δ′)] rather than e^{o(Λ)}, which only
changes c′. Fix also a deterministic in-field
descent procedure 𝔇, for example nested's DESC (first good point, recursing) or the bounded-depth feasibility search of
§7.3. Requirements on 𝔇:
- on every admissible field and every target, its total work per call, including all backtracking and failed
  branches, is ≤ L[α, s′] for a constant s′, excluding oracle calls;
- it calls the recursive factoring oracle only on integers ≤ N^{τ_F}, τ_F < 1;
- it outputs FAIL, or a factor exit, or an RW tree (nested Lemma RW) whose leaves lie in the factor base.

(The draft required "e^{o(Λ)} tree nodes", which excludes the bounded-depth search; proof critic #6.) For a rational
prime q put p_q := Pr_{A∼𝔉}[a factor exit occurs, or A is admissible and 𝔇_A(q) returns a tree]. Factor exits count as
success, as in lv5unif Prop MF.

**Hypothesis T_q(α; 𝔉, 𝔇, c, C₀).** Let k := ⌈κD^α⌉, Y := ⌊N^{1/(2k)}⌋, and J := primes in (Y, 2Y] not dividing N.
Then #{q ∈ J : p_q ≥ e^{−cΛ_α}} ≥ (log N)^{−C₀}·#J.

**Theorem NR-q [proved conditional on T_q; on nested's Hypothesis LV6′ restated at degree d ≍ D^δ with overheads
e^{o(Λ_α)}; on character columns via Lemma X_ERH (ERH) or lv5unif's factorial-bit repair (cost e^{O(d log d)},
which is e^{o(Λ_α)} for δ < α, and L[1/3, Θ(δ)] at α = 1/3 with d = δD^{1/3}); and on the multi-field post-processing of
lv5unif §11 (proved there in outline; codex escapeR2_1 (2.5) gives the square-root identity)].** There is a Las Vegas
algorithm that completely factors every N in expected time L_N[α, c′] for a finite c′ = c′(α, c, C₀, κ, s′, …).

*Proof.*
1. **Pool.** Sample J_f := ⌈4e^{cΛ_α}⌉ fields A_1, …, A_{J_f} ∼ 𝔉 independently (factor exits end the run). Let
   G_pool := {q ∈ J : ∃j, 𝔇_{A_j}(q) returns a tree or a factor exit} (the same disjunction as in p_q; codex
   escapeR2_3 (a)). A factor exit met while descending a split target ends the run successfully. For a good q (p_q ≥ e^{−cΛ_α}), Pr[q ∉ G_pool] ≤ (1 − p_q)^{J_f} ≤ e^{−4}.
   By linearity the expected number of missed good q is ≤ e^{−4}·#good, with no independence between different q
   needed. By Markov, with probability ≥ 1 − 2e^{−4} at least half the good q are in G_pool, so |G_pool| ≥ γ#J with
   γ := (log N)^{−C₀}/2. Call such a pool *good*. G_pool is never enumerated; membership of a given q is tested by running
   𝔇 in the pool order.
2. **Split, on a padded modulus (repair by codex escapeR2_2 #15).** With Y = ⌊N^{1/(2k)}⌋, every k-product of primes in
   (Y, 2Y] exceeds about √N. It therefore cannot be the Gauss-shortest vector modulo N, so Lemma S″ cannot be applied
   modulo N as the draft did. Instead:
   - After preprocessing (3 ∤ N), put N_* := 9^{k+1}N, draw r uniform among the units modulo N_* (non-units give gcd
     exits or rejections), and Gauss-reduce R = r² modulo N_*.
   - All k-products of G_pool lie in (Y^k, 2^kY^k]. One of the k dyadic windows inside contains ≥ C(|G_pool|, k)/k of
     them, and its upper end X satisfies 2X² ≤ 2·4^kN < N_*.
   - Lemma S″ applies modulo N_*: its square-class index is 2^{ω(N)+1} ≤ |G_pool|/(4k³), by lv5unif (P6) at general α.
     Acceptance is ≥ (γ/(4 log N))^{2k}·e^{−O(k)}/poly(k) = L[α, −2κ(1 + C₀) − o(1)], since k = o(Λ_α).
   - The row is read modulo N. Halves are ≤ √(2N_*) = N^{1/2+o(1)}, which preserves the recursive contraction.
   - By CRT, given R mod N_*, the reduction r mod N is uniform among the square roots of R mod N. The hidden-root
     argument survives.
3. **Descent.** For each of the 2k primes of an accepted split, run 𝔇_{A_1}, 𝔇_{A_2}, … in order until success. That is
   at most J_f·L[α, s′] work, plus oracle calls on integers ≤ N^{τ_F}. By nested Lemma RW each tree gives
   ℓ ≡ (rational)·ψ_j(Z) with Z supported on factor-base primes of A_j.
4. **Rows and columns.** Columns: the shared rational ones, plus for each pool field its parity, (p, r) and character
   columns with t_χ^{(j)} = ⌈5d(ln|H_j| + 2 ln log N + ln 4J_f)⌉. That is ≤ J_f·L[α, b + o(1)] columns in total. Collect
   M_rows := #columns + J_f + 12 rows.
5. **Post-processing.** lv5unif §11, written for this pool: per-field Lemma X (ERH, or the factorial-bit repair at
   cost e^{O(d log d)}); a J_f × (J_f + 2) bit matrix from Algorithm SQ to kill the residual characters; exact
   verification; codex (2.5) for Y. The cost is polynomial in J_f, M_rows and the expanded row size, i.e. L[α, O(1)].
6. **Fruitfulness.** 𝔇 is deterministic and the pool is chosen before the Dixon draws. Every exposed quantity is a
   function of the R_i and coins independent of the hidden roots (nested Lemma T_multi, with the pool as part of the
   coins). So each LV6 round succeeds with probability ≥ 1/2.
7. **Restarts (doubling schedule).** Let T* be the good-pool expected time *until a verified proper factor*, including
   the inner LV6 loop: character resampling, and fresh rows after a trivial gcd (each round succeeds with probability
   ≥ 1/2 by step 6). T* is not computable in advance (c, C₀ are not explicit). Run rounds i = 1, 2, …: in round i, sample a fresh pool and run steps 2–6 with a work budget 2^i,
   aborting on overrun.
   - Once 2^i ≥ 4T*, a round succeeds with probability ≥ p₀ := (1 − 2e^{−4})·3/4 > 0.7: the pool is good, and by Markov
     the cost is within the budget.
   - The total expected work is ≤ Σ_{i ≤ i₀} 2^i + Σ_{j ≥ 0} 2^{i₀+j}(1 − p₀)^j = O(T*), because 2(1 − p₀) < 1.
   - (The draft's "probability ≥ 1/4" would not make this sum converge.)
8. **Recursion.** F(N) ≤ L[α, p₀] + L[α, δ₀]F(N^{1/2+o(1)}) + L[α, δ₁]F(N^{τ_F + o(1)}) with finite exponents. Strong
   induction with Λ_α(N^τ) = (τ^α + o(1))Λ_α(N) (nested's proof of NR; lv5unif Lemma F_α for the first two terms)
   gives L[α, c′]. ∎

*Relations to known hypotheses.*
- AD (nested §8) ⇒ T_q with p_q ≥ 1/2 − o(1) for every q ∈ J, using DESC (factor exits are in both events; the q | D_A
  are o(#J)). So NR-q contains NR, at α = 1/3 with nested's parameters.
- T_q is field-averaged and target-pointwise. That is the form Prop SI(ii) allows: the root target q is exogenous to
  the pool. The endogeneity now lives entirely inside one field's tree, i.e. in Hypothesis LS_C.
- A sufficient condition for T_q is r20 #1's moment pair at L^{O(1)} precision (§7.3), when 𝔇 is exhaustive at
  bounded depth.
- *What T_q is not:* it is not proved for any α < 2/5, unconditionally or under ERH. At α = 2/5 with one level it
  follows from lv5unif's Theorem D-U, in the stronger, density form.

---

## 9. The missing lemmas, collected (a catalogue, not a reduction)

These are separate sufficient inputs. We do not show that the others are stronger than item 1. As both reviewers
stress, the failed implementations in §§3–6 do not prove that every route reduces to item 1.

1. **Hyp LS_C → T_q** (in-field joint smooth supply: targetwise moments plus search; §§7.2, 8). This is the most
   precise target. At C = 3/2 it gives 3/8 in the chain sub-model, and as C → ∞ it gives → 1/3.
2. **Contracting image certificate** (codex (5.2)): for θ from principalisation, ψ(θ) ≡ w²·a·∏q_ν (mod N) with
   q_ν ≤ q^λ, λ < 1. We know no way to produce this other than by a re-split (§4).
3. **Coset re-split** (Prop RS): class equidistribution of products of the actual descendable primes in every coset,
   plus a contracting bridge not of uniform-mask type. False for adversarial sets; open for the actual ones.
4. **Useful transfer between tied fields** (§2.3, §5): fields with parity transfer (k_1 ≠ Q), a compatible evaluation
   map, short transfer elements and enough conditional randomness (e.g. Q(√a), Q(√b), Q(√ab) with b fresh). No
   construction is known. KT excludes the hyperelliptic 2-torsion-preserving version for S_d-fields of odd degree ≥ 7.
5. **Directed recombination** (§6): rows constrained to share large ideals, with a seeded rank-surplus theorem.
6. **Small-discriminant fields with explicit ψ** (§4): these would make L_{|Δ|}(1/2) = L_N[1/3]. For a small fixed field
   this is root-finding modulo composites, which is not known to be easier than factoring.

---

## 10. Experiments (pre-registered, `notes/escapeR2_prereg.md`, commit `4624c3b`; code frozen at `42312c2`)

Command: `.venv/bin/python -m src.escapeR2.run_experiments experiments/raw` (44 s).

| | prediction | result |
|---|---|---|
| E1 P1 | dim(Im_1 ∩ ΣIm_j) ≤ dim F₂^{H₁⁺\G/D} always; τ = 0 for disjoint fields | **0 violations** in 51 product-group cases (S_3², S_3×S_4, S_3³, all D) and 65 940 census triples; 0 of 7 948 triples with k_1 = Q have τ > 0 |
| E1 P2 | V₄: τ = 1 | **τ = 1** (pairwise τ₂ = 0, 0, 0) |
| E1 P3 | two fields: τ₂ = 0 | **0 violations** in 6 832 census pair cases + 51 product-group pairs |
| E1 P4 | descriptive | τ > 0 in 40 149 / 65 940 triples (S_4: 32 080/49 140; D_4 456/840; A_4 330/672; C₂³ 4 074/8 736; S_3×C₂ 3 209/6 552), all with k_1 ≠ Q |
| E2 Q1 | (r + 2)/(3r + 5) | **holds** for r = 0…8: feasible at the value, infeasible 10⁻⁶ below |
| E2 Q2 | α*(1) = 2/5 | **holds** (infeasible below 2/5 for all r ≤ 8) |
| E2 Q3 | α*(3/2) = 3/8; (r + 2)/(3r + 5) for C ≥ r + 1; > 1/3 | **holds** at all 10 grid C; LP = closed form to 10⁻⁹ |
| E2 Q4 | children exponent 3/8 at α = 3/8 | **3/8** (C needed: 3/2) |
| E3 R1 | 0 pairs, 0 hits for non-square w | **0 and 0** on all 5 seeds |
| E3 R2 | x = 1 pairs ≥ |S|²/(2h) | **4.02·10⁷ ≥ 2.06·10⁷** (hypothesis |G| ≥ 4k³h holds: 255 ≥ 128) |

*Disclosed deviations.*
- E1(c)'s pair-case count is 6 832 including both orders; the prereg did not fix a number.
- E3 used 48-bit moduli for all seeds (within the registered 40–48 range). The product window is [Y′/2, Y′] with
  Y′ = 2^{⌊(n−2)/2⌋} (so 2Y′² < N, as Lemma S″ needs); the prereg said "the top window" without fixing Y′.
- E3 (codex escapeR2_2 #8–9): the adversarial G (59–70 primes) are below S″'s size threshold 4k³h = 128, so E3 is an
  illustration. The full-hypothesis counterexample is the pigeonhole construction of Prop RS(ii). w was drawn with
  (w/p₁) = −1 only; in seed 3 it lies in −Q_N, so the zero-hit record there is luck, not the obstruction.
- Post-review code change: `sizemodel.closed_budget` now returns "infeasible" for r > C (codex #13), and
  `tests/test_escapeR2.py` gained `test_closed_budget_matches_lp_off_grid` (α*(13/5) = 7/19). E2's registered outputs
  are unchanged: the LP never used the helper, and the global minima agree.

E1–E3 check algebra and arithmetic only. They say nothing about relation supply.

---

## 11. Literature

From the Opus literature critic (re-spawned instance). Sources read first-hand where stated in §12; the rest is
flagged.
- **Norm relations.**
  - Lemma TI's double-coset parity form was not found in the literature. It is elementary and close to standard
    material; **no novelty claim.**
  - Biasse, Fieker, Hofmann, Page, "Norm relations and computational problems in number fields", J. London Math. Soc.
    (2) 105 (2022) 2373–2414 (arXiv 2002.12332).
  - Kani, Rosen, Math. Ann. 284 (1989) 307–327, and J. Number Theory 46 (1994) 230–254.
  - Kuroda's class number formula is the classical V₄ instance (Lemmermeyer, Acta Arith. 66 (1994); from memory). The
    corresponding Brauer relation is [G/1] − Σ[G/C_i] + 2[G/G] = 0.
  - Multiquadratics: Bauch, Bernstein, de Valence, Lange, van Vredendaal, EUROCRYPT 2017, LNCS 10210, 27–59; Biasse,
    van Vredendaal, ANTS XIII, Open Book Ser. 2 (2019) 103–118.
- **Class groups.**
  - dBPW, arXiv 2512.01588: ERH for Hecke L-functions; poly·L_{|Δ|}(1/2)·L_{n^n}(2/3)·min(ρ_K, L_{|Δ|}(2/3+o(1)));
    rigorous.
  - Hafner–McCurley, J. AMS 2 (1989): imaginary quadratic, ERH.
  - Biasse, Math. Comp. 83 (2014) 2005–2031 (arXiv 0912.1927), and arXiv 1204.1292 (DLP / principality): heuristic
    L(1/3), special families.
  - Biasse–Fieker, LMS JCM 17A (2014) 385–403: GRH plus heuristics; L(2/3+ε), or L(1/2) for n ≤ (log|Δ|)^{3/4−ε}.
  - Gélin–Joux, LMS JCM 19A (2016) 315–331: heuristic L(1/3).
  - Bach, Math. Comp. 55 (1990) 355–380: under ERH for Hecke L-functions of class-group characters, Cl(O_K) is
    generated by primes of norm < 12 log²|Δ_K|. Grenié–Molteni, Math. Comp. 87 (2018) 2483–2511: degree-one primes,
    (4 + o(1)) log².
  - **No rigorous L(1/3) class-group or PIP algorithm was found, 2023–2026 included.**
- **Smooth numbers in linear systems.**
  - Harper, Compositio 152 (2016) 1121–1158: a + b = c, y ≥ log^C x.
  - Matthiesen–Wang, Algebra & Number Theory 19 (2025) 1881–1946: finite-complexity linear forms.
  - Lachand, Acta Arith. 181 (2017) 97–109: fixed u.
  - Balog–Blomer–Dartyge–Tenenbaum, Comment. Math. Helv. 87 (2012) 639–667.
  - **None has a growing number of forms, or u at L-scale.**
- **2-descent.** Schaefer, J. Number Theory 51 (1995) 219–232 (odd degree); Stoll, Acta Arith. 98 (2001) 245–277; Poonen–Schaefer, Crelle 488 (1997) 141–188 (from memory); see also Morgan, Proc. LMS 131 (2025) e70066.
- **Permutation groups.** Bochert's theorem (index ≥ ⌊(n+1)/2⌋! for primitive groups not containing A_n), as cited in
  Dixon–Mortimer and by Babai (arXiv 1512.03547). Dixon–Mortimer Thm 5.2B. Fact (F) is proved directly in §5.2.
- **Rabin**, MIT/LCS/TR-212 (1979): random-residue square roots are equivalent to factoring. For fixed c, nothing is
  known.
- **Factoring.**
  - No rigorous or ERH-conditional algorithm below L[1/2] for arbitrary N was found (2024–2026 searched).
  - Lenstra–Pomerance, J. AMS 5 (1992), L[1/2, 1 + o(1)], remains the record found.
  - Umans–Wang (arXiv 2511.10851) is conditional and deterministic-exponential; it is irrelevant here.
  - No prior formulation of an overdetermination / parameter-budget obstruction for rigorous NFS descent was found. The
    search was a single targeted one, so this is low confidence.

---

## 12. Reviews

All valid findings were applied in the text above. Nothing was rejected outright; two items were applied in narrowed
form (noted below).

**GPT-6.1-sol blind round `codex/escapeR2_1` (xhigh; the whole problem, before my draft).**
- *Verdict:* no rigorous α < 2/5, unconditionally or under ERH; the same per-family failures.
- *Contributions adopted (credited in place):*
  - the parity invariant (2.2) and the multi-field square root (2.5);
  - the norm-substitution failures and odd-degree norm parity (§3);
  - mod-square principalisation (5.1) and the contracting-image lemma (5.2);
  - the S-unit dimension bound (8.1);
  - the square-corrected re-split and the coset pair count (10.1)–(10.2);
  - the level thresholds (r+1)/(3r+2);
  - the 3/8 two-level table with its D^{3/8}-vs-D^{1/4} overdetermination count.
- Its ranking put "fixed-field completed descent at 3/8" first, which is what §7 sharpens.

**Orchestrator's wild round `codex/r20_wild`.**
- #1 (complete-certificate moments) is evaluated in §7.3. Its precision requirement is relaxed to L^{O(1)}, fieldwise
  and targetwise; the core estimate stays open.
- #4 (Kummer transport) is evaluated in §5.2, with a structural kill of its 2-torsion-preserving hyperelliptic form for
  S_d-fields.
- #7–#9 kills checked and agreed (§§5.3, 6).

**GPT-6.1-sol adversarial referee `codex/escapeR2_2` (xhigh).** TI CORRECT; EX, KT, RS, R3, chain model, NR-q and
precision CORRECT WITH FIXES; exhaustive Verdict WRONG. 19 issues:

| # | sev. | issue | adjudication |
|---|---|---|---|
| 1 | MAJOR | §1: a common evaluation on the compositum need not exist (√2, −√2 counterexample) | accepted: §1 rewritten; the general class is a superset, ψ_L must be supplied |
| 2 | MINOR | TI excludes ℓ ramified via other fields | accepted: excluded set stated, O(J d t/Λ) primes |
| 3 | MINOR | non-monic index / norm normalisation in EX | accepted: localised index, F(a,b)/f_d |
| 4 | MAJOR | tied fields may keep conditional randomness | accepted: §2.3 rewritten as an open question |
| 5 | MINOR | KT faithfulness via one transposition | accepted: general argument |
| 6 | MAJOR | generator-height count false (α + 1) | accepted: claim withdrawn; Möbius remark instead |
| 7 | MAJOR | KT does not kill all transports | accepted: scope narrowed in Verdict, table, §5.2 |
| 8 | MINOR | Chebotarev density in a growing window; E3 sets below S″ threshold | accepted: pigeonhole construction; E3 relabelled illustration |
| 9 | MINOR | E3 seed 3 has w ∈ −Q_N | accepted: disclosed |
| 10 | MAJOR | "every re-split" / bridge classification unsupported | accepted: RS(iii) restricted to uniform masks (codex's proved bound); recursion paragraph restated |
| 11 | MAJOR | R3 recurrence drops the incoming target; constant fails (d = 5 chain) | accepted: recurrence corrected (my own pre-review fix coincided), constant 2(1+2κ) |
| 12 | MINOR | R3 weights at d = 3 | accepted: O_κ bound for d = 3 |
| 13 | MINOR | C ≥ r omitted in closed form / `closed_budget` | accepted: formula and code fixed, test added |
| 14 | MAJOR | LS_C as an averaged count controls neither overlaps nor search | accepted: LS_C restated as targetwise moments + search |
| 15 | MAJOR | NR-q split window invalid (k-products > √N) | accepted: codex's padded modulus N_* = 9^{k+1}N adopted |
| 16 | MINOR | uniform costs and implementable runtime cap | accepted: doubling schedule, with corrected convergence arithmetic; stopping time = verified factor incl. LV6 retries (escapeR2_3) |
| 17 | MAJOR | jointly averaged moments can concentrate | accepted: targetwise quantifier made explicit, with codex's counterexample |
| 18 | MINOR | TS bounds successful trees only | accepted: search recurrence T_h ≤ W(work + mT_{h−1}) |
| 19 | MAJOR | failed implementations promoted to an exhaustive reduction | accepted: Verdict and §9 rewritten as a catalogue |

**GPT-6.1-sol confirmation round `codex/escapeR2_3` (high; on commit e6a1adb).**
- Status of the 19 issues: 15 FIXED; #2, #3, #4 and #16 PARTLY FIXED. The residuals are now applied: the ramification set
  restricted to the compositum of Galois closures; the localised index in §1; softened §2.2 wording; the stopping time.
- NR-q: the padded split passes. Remaining bookkeeping, now applied:
  - the factor-exit disjunction in G_pool;
  - sampling units modulo N_*;
  - T* counted until a verified factor;
  - admission of the α = 1/3 boundary case.
- R3: the tree proof passes. The DAG remark was replaced by the topological-order bound 2e^{2κ}ΛU.
- Chain envelope and 115 exact checks pass. A plateau wording slip was fixed.
- The Verdict no longer claims an exhaustive reduction.
- Bach: its "RH for ζ_K alone" finding referred to e6a1adb; it was already fixed to ERH for Hecke L-functions in 71701e4.
- *Overall:* NR-q CORRECT WITH FIXES (now applied). No FATAL defect. No unconditional improvement below L[2/5].

**Opus proof critic (read-only, independent).** TI, EX, TB CORRECT; KT, PR/Bach, R3, NR-q CORRECT WITH FIXES; the
recursion claim WRONG as worded; §7.2's general α*(C) equality WRONG (grid values right). No FATAL.

| # | sev. | issue | adjudication |
|---|---|---|---|
| 1 | MAJOR | §1 extension of ψ backwards (x² − 2 and V₄ counterexamples) | accepted (= codex #1) |
| 2 | MAJOR | "recursion cannot replace descent [proved]" contradicts RS(i) | accepted: relabelled structural, restated |
| 3 | MAJOR | precision relaxation only targetwise; the density is still polylog | accepted (= codex #17); Verdict wording fixed |
| 4 | MINOR | α*(C) plateau formula false off the grid: α*(13/5) = 7/19 | accepted: exact lower-envelope formula, plateaus [r+1−1/(r+1), r+2−1/(r+1)], new table row and test |
| 5 | MINOR | R3: z_v omitted; κ ≥ 1 misses d + 1; tree assumption | accepted: κ ≥ 2, DAG remark |
| 6 | MINOR | NR-q's node bound excludes exhaustive search | accepted: total-work condition |
| 7 | MINOR | α = 1/3 excluded; factor exits in p_q; AD ⇒ p_q ≥ 1/2 − o(1) | accepted |
| 8 | MINOR | Y window, ambiguous reference, LV6′ at D^δ, budget, ERH in Verdict | accepted (window via codex's padded modulus) |
| 9 | MINOR | chain-model floors are upper bounds on the model floor | accepted: labelled throughout |
| 10 | MINOR | PZ gives L^{−c}, not L^{−2c}; "cheapest step" inaccurate | accepted |
| 11 | MINOR | KT even faithfulness; (F) at d = 7, 8; Verdict hypotheses; height wording | accepted: (F) now proved directly from Bochert's theorem for d ≥ 7 (codex read 5.2B as d ≥ 7, the literature critic as n > 8); height claim withdrawn (codex #6) |
| 12 | MINOR | Bach needs all prime ideals and GRH(ζ_K) | accepted |
| 13 | MINOR | index undefined for non-monic f | accepted (= codex #3) |
| 14 | MINOR | Rabin attribution | accepted: "not known to be easier than factoring" |
| 15 | MINOR | table rows 3 and 4b overstated | accepted |

**Opus literature critic (read-only, independent; re-spawned after the first instance was lost to an API limit).**

*Sources.* Read first-hand: the arXiv abstracts of 0912.1927, 1204.1292, 1408.1662, 2002.12332, 2211.16892 and
1609.08872; part of the dBPW v2 HTML; Crossref metadata; Wesolowski 1807.01561 on Bach; Umans–Wang. Not read:
Dixon–Mortimer, Bach, Schaefer, Stoll, Morgan and Biasse–Fieker in full.

*Corrections applied:*
- Fact (F)'s citation: 5.2B may assume n > 8. We now prove (F) directly from Bochert's theorem for all d ≥ 7. The
  critic and codex disagreed on 5.2B's range, and the direct proof makes that disagreement moot.
- An even permutation is needed for faithfulness; the general argument of §5.2 covers it.
- Schaefer covers the odd case only; Stoll and Poonen–Schaefer cover the even case; Morgan is a "see also".
- Bach's hypothesis is ERH for Hecke L-functions; Grenié–Molteni added.
- dBPW's L_{n^n}(2/3) factor.
- Biasse–Fieker is not L(1/3); Biasse's arXiv numbers correctly assigned; Gélin–Joux added.
- The Matthiesen–Wang range: the nested-inherited "Ψ ≥ x/log⁸x" dropped.
- Rabin is for random residues; fixed c is open.
- The V₄ identity is a norm relation (denominator 2), not literally a Brauer relation; Kuroda added.

*Confirmed:*
- no rigorous L(1/3) class-group or PIP algorithm;
- no smooth-number result for growing systems at L-scale;
- no rigorous or ERH factoring below L[1/2];
- no prior statement of TI's parity form (no novelty claimed);
- no prior budget-obstruction formulation found (low confidence).

---

## 13. Candidate card

- **Mechanism.** Escape nested Prop R2 in one of three ways:
  - export algebraic children (norms, fresh shared-root fields, generators, Galois/Kummer transport);
  - class-group principalisation;
  - change the problem (re-splits, recursion, independent descents).
- **Model.** The Theorem U / nested architecture; the certificate class of §1; the size model of recurse §3 and nested
  R2.
- **Result.** No rigorous or ERH-conditional α < 2/5.
  - Proved:
    - Lemma TI and Cor EX: no export across fields with k_1 = Q, at primes unramified in the common field;
    - Example TB: tied fields transfer parity; usefulness open;
    - Prop KT: 2-torsion-preserving hyperelliptic transport returns to K, for S_d with odd d ≥ 7;
    - Lemma PR;
    - Prop RS(i)–(iii): RS(iii) for uniform masks.
  - Model:
    - Prop R3 (O(d)-point components ⇒ 2/5);
    - chain-sub-model floors α*(C) (E2): 7/18 at C = 5/4, 3/8 at C = 3/2, → 1/3.
  - Conditional: Theorem NR-q (T_q ⇒ L[α]), which relaxes r20 #1's moment ratio to L^{O(1)} fieldwise and
    targetwise.
- **Failure scope.** Square certificates built from field elements; designs under (M1)–(M4); the inspected
  implementations. Not covered:
  - useful transfer between tied fields;
  - directed recombination;
  - structured image certificates;
  - fields with explicit ψ and small discriminant;
  - non-square certificates;
  - structured polynomial selection.
- **Next step.** The smallest case of Hyp LS_{3/2}: two-level certificates at α = 3/8. It needs a targetwise
  L^{O(1)}-precision second moment, i.e. any lower bound for smooth points of an overdetermined affine system of d^{3/2}
  forms in d variables at L-scale.

## Reproduction

- `.venv/bin/python -m src.escapeR2.run_experiments experiments/raw` (E1–E3, 44 s).
- `.venv/bin/python -m pytest -q tests/test_escapeR2.py`.
