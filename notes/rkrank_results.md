# rkrank: Hypothesis RK (GF(2) rank of NFS relations) and rigorous small-square samplers

Agent `rkrank`, phase 8, 2026-10-19.
- Pre-registration: `notes/rkrank_prereg.md`, committed in `7a2d0fe` before the E1 run.
- Code: `src/rkrank/`. Tests: `tests/test_rkrank.py` (6 tests).
- Raw data: `experiments/raw/rkrank_e1.jsonl`, plus the supplementary run
  `experiments/raw/rkrank_e1_p5supp.jsonl` (disclosed deviation, §5).
- GPT-6.1-sol rounds: `codex/rkrank_1.{md,out}` (blind, lever A), `codex/rkrank_2.{md,out}` (blind,
  lever B) and `codex/rkrank_3.{md,out}` (adversarial review of the first draft).
- Reviews: two independent read-only Opus critics, one on proofs and one on literature and novelty. All
  three reviews have been adjudicated and applied (§7).
- Notation: **K_V** is the sampler's quasi-uniformity ratio (hypD's κ). **κ** is reserved for the
  codimension exponent in Theorem V_κ.

## Verdict

**Neither lever succeeded. We did not prove Hypothesis RK. We found no rigorous small-square sampler
with θ < 2/3, either in the literature or by our own attempts. Lenstra–Pomerance's L[1/2, 1+o(1)]
(1992) stands, and none of the four searches (two GPT-6.1-sol, two Opus) found a published rigorous
c < 1.** Design V remains conditional at c = √(2/3).

What remains is a sharper map of the problem: what RK reduces to, which weaker form already suffices,
which obstructions are forced, and where the natural proof templates fail. Pre-registered toy
experiments support the sparse form of RK ("W = Λ/2 on the well-covered primes"). They say nothing about
the asymptotic regime.

**Lever A (RK).** Labels in brackets.
- **Theorem V_κ [proved implication; routine cost balance].** RK may be weakened to codimension D^κ, where
  D is the number of factor-base columns. With sampler exponent θ and LV size exponent η = 3/(d+1), the
  cost constant is c = √(2·max{η, θ/(2−κ)}).
  - For Vallée (θ = 2/3) and d ≥ 9, c < 1 iff **κ < 2/3**.
  - So deficiency D^{2/3−ε} would already beat Lenstra–Pomerance, conditionally on Design V's other
    inputs. RK as originally stated is κ = 0. We found no prior formulation of this weakening.
- **Theorem U [proved; routine coupon-collector argument].** Fix a forced space G′ given in advance, e.g.
  the algebraic identities of Prop C. Suppose every dual space V with V ∩ G′ = 0 and dim V > r is escaped
  by a random relation with probability ≥ dim V/T. Then T(1 + ln D) relations reach deficiency
  ≤ r + dim G′, up to a constant-probability Markov loss.
  - The hypothesis ("expansion of the smooth set X") is the real content. It is open for NFS.
  - It is the k-graded analogue of de Boer–Ducas–Fehr's evenly-distributed generation lemma.
- **Prop G′ [proved, existential; classical, Buhler–Lenstra–Pomerance 1993 §8].** Raw matrix membership
  and genuine certificates differ by at most dim H_K ≤ r₁ + r₂ + log₂ h_K = O_d(log N) dimensions, plus
  bookkeeping that is asserted but not written out. So r + O(log N) + 1 sampler rows force a genuine
  dependency with nonzero sampler projection. Efficient extraction still needs LV §6.
- **Prop C [proved; elementary plus Weil].** This classifies the quadratic-character identities
  χ(a − bm)·ψ(a − bα) ≡ 1 that hold on **all** coprime pairs.
  - Away from primes ℓ ≤ 4(d+1)² and ℓ | 2f_d·disc f, they are exactly the ω characters at p | N.
  - The excluded primes contribute O_d(log N) dimensions.
  - There are no archimedean identities on all of Z².
  - **Box-level obstructions are separate.** In LV's non-monic setting, e_sign + e_{|S|} ∈ R^⊥ always,
    because a − bm < 0 on the box. A real place whose sign is constant 1 on the box adds at most r₁ more.
  - Prop C gives **no quantitative RK bound**. An obstruction not explained by an algebraic identity
    might still hold on the whole box, or on X.
- **Lemma S [proved; switching argument from codex rkrank_1].** For uniform B-smooth integers, n ↦ p·n
  switching gives exactly the escape rate Theorem U needs (the rational-side model). No such move is
  known for the joint kept law of one NFS field.
- **Barriers [proved, for the stated templates].**
  - (B-F) A per-functional union bound averaged over ≤ N^{O(1)} fields cannot reach the 2^{−D} scale. A
    spread construction (codex) has every fixed functional fail in ≤ 1 field while every field has
    deficiency D/2.
  - (B-M) Certifying span through collision probabilities of k-fold sums needs k ≥ D/L^{o(1)}.
  - (B-C) Coverage does not imply span.
  - These rule out templates only, not proofs of RK in general.
- **Experiment E1 [experimentally supported; pre-registered; toy monic cubic NFS at 64/80 bits].**
  - W ⊆ Λ/2 in 72/72 records.
  - W ⊇ Λ/2 on primes occurring ≥ 3 times (c₃ = 0) in all 32 records at ≥ 1.5·cols. It fails in 5
    records at 1.0–1.25·cols.
  - 13 fixed Liouville-type and random dense functionals are balanced: δ ∈ [0.45, 0.55] in 20/20
    instances, using the supplementary run (§5).
  - Two sharper predictions failed narrowly. P1 said peeling explains every rare-prime obstruction; P2b
    said c₁ = 0 at ≥ 2·cols. The cause was ≤ 5 rare-prime obstructions, never one on well-covered primes.

**Lever B (samplers).**
- **Prop Q [proved implications; classical at s = 0].** Take a fibre-symmetric sampler of exponent θ and
  per-call cost L^{s+o(1)}.
  - **Dixon-style: c = s + √(2θ)**, with no rank hypothesis. It is < 1 iff s + √(2θ) < 1, so θ < 1/2 alone
    would beat LP.
  - **Design V** (under RK and its inputs): c = max{√(2η), (s + √(s² + 4θ))/2}. At s = 0 this is √θ, as hypD
    claimed, given d + 1 ≥ 6/θ.
- **Theorem L (linearisation barrier) [proved, restricted model].** Take a sampler whose random cell has
  uniformly small integer quadratic part and covers each point with probability ≤ λ·|cell|/N.
  - (i) For **rank-one cells** x₁ + vτ, v² ≡ σ, the expected yield is ≤ λ2^ω(2h+1)(2√h+1)/N, for every
    θ < 2/3. So N^{1−3θ/2−o(1)} cells are needed per output.
  - (ii) For **cells of any rank**, the same holds for θ < 1/2, by a rank-collapse lemma.
  - For 1/2 ≤ θ < 2/3, cells of rank ≥ 2 that do not collapse exist (e.g. w₂ ≈ √N). The barrier is open
    there.
  - The standard Coppersmith/Howgrave-Graham determinant guarantee for (x₁ + t)² − z mod N requires
    T²h < N, which is incompatible with Th ≥ N^{1−o(1)}. That guarantee is unavailable, not provably
    useless.
- **Literature.** Post-1991 progress on Vallée's sampler is in uniformity, not exponent.
  - Coron (CRYPTO 2002): wider intervals, statistically close to uniform.
  - Gentry (CRYPTO 2004 §4): exactly uniform at width ≈ 8N^{2/3}, per codex rkrank_2/3. The literature
    critic could confirm only quasi-uniformity from secondary sources, so this is unverified by us.
  - No θ < 2/3 was found.

---

## 1. Setting and notation

These follow hypD §5.
- Field A = (m, f), with f = f̂_{N,m} + (x − m)R, f(m) = N, fixed odd d and B = L[1/2, a].
- A relation is a coprime (a, b) in a box of size L^{O(1)} with a − bm and F(a, b) both B-smooth.
- Its row consists of:
  - rational parities q(a, b) ∈ F₂^{D_rat}: the sign and v_p(a − bm) for p ≤ B;
  - a high block z(a, b): degree-one ideal parities, real-embedding signs, LV characters, and the
    |S|-parity column when f is non-monic.
- R is the row space and W = {u : (u, 0) ∈ R}.
- Λ/2 = {u : ∏p^{u_p} ∈ (Z/N)^{×2}}.
- Gen = (Λ/2)^⊥, spanned by the ω functionals w_i(p) = [(p/p_i) = −1].
- U is the target space: Λ/2 restricted to a set of well-covered primes, or the span of the sampler
  supports.
- Deficiency means codim_U(W ∩ U).
- X = X_{m,f} is the set of smooth pairs in one field's box.

## 2. Lever A: Hypothesis RK

### 2.1 Duality [proved; elementary]
- W^⊥ = proj_rat(R^⊥) and codim_U(W ∩ U) = dim res_U(W^⊥).
- *Proof.* A rational functional w annihilates W iff some (w, y) annihilates R. ∎
- Also dim W = rank R − rank(z-block). The code asserts this (`harness.level_record`).
- So a **dual obstruction** is a pair (w, y) with w(q(a, b)) = y(z(a, b)) for every relation.

### 2.2 Prop G′: raw versus genuine certificates [proved, existential; classical]
1. Let E₀ be the combinations S whose ideal parities and chosen character columns vanish.
2. Then ∏_S(a − bα) has a square ideal, so it defines a class
   ρ(S) ∈ H_K := {x ∈ K^× : (x) a square of an ideal}/K^{×2}.
3. H_K fits into 0 → O_K^×/O_K^{×2} → H_K → Cl(K)[2] → 0, so dim H_K ≤ r₁ + r₂ + log₂ h_K = O_d(log N).
   This is Buhler–Lenstra–Pomerance 1993 §8; van Leeuwen's thesis restates it. Codex rkrank_1 rederived
   it.
4. Order, index, denominator and non-monic bookkeeping are asserted to add O(log N). This is not written
   out.
5. Genuine certificates form ker ρ, so dim W_raw − dim W_true ≤ dim H_K + O(log N).

**Consequence.** If the raw deficiency on U is ≤ r, then t > r + dim H_K + O(log N) sampler rows supported
in U give a genuine dependency with nonzero sampler projection.
- It is fruitful with probability ≥ 1/(2K_V), provided dependency selection does not inspect the sampler
  roots.
- Efficient extraction needs LV §6.

### 2.3 Theorem U: uniform rank growth [proved; routine]
**Setup.** Z₁, Z₂, … are i.i.d. with law μ on F₂^s. G′ ⊆ Ann(μ) is a subspace fixed in advance.
Δ_M = span(Z₁..Z_M)^⊥ ⊇ Ann(μ) ⊇ G′.

**Hypothesis (Exp_T(G′)).** For every subspace V with V ∩ G′ = 0 and dim V = k > r,
Pr_μ[Z ∉ V^⊥] ≥ k/T.

**Theorem U.** Under Exp_T(G′):
- E min{M : dim Δ_M − dim G′ ≤ r} ≤ T(1 + ln s);
- in particular dim Ann(μ) ≤ dim G′ + r;
- codim_U(W ∩ U) ≤ dim res_U(proj_rat G′) + r after M relations, except with probability
  ≤ T(1 + ln s)/M.

*Proof.*
1. Write Δ_M = G′ ⊕ V_M. Escaping V_M lowers dim Δ by one, with conditional probability ≥ dim V_M/T.
2. The hypothesis covers every eligible V, so adaptivity is harmless.
3. Sum the geometric holding times and apply Markov's inequality. ∎

*Remarks.*
- With G′ = Ann(μ), the theorem is only about finite-sample convergence to span(X) and says nothing
  about RK. With G′ the algebraic space of Prop C, Exp_T itself excludes "accidental" obstructions of X:
  their escape probability is 0. So Exp_T(G′_alg) is the full content of RK. It is a deterministic
  expansion property of the finite set X.
- The natural scale is T = D·L^{o(1)}: rare primes force a k-dependent rate. The frequency ≈ u/B of a
  prime p ≈ B in a relation is heuristic for the joint law.
- This is the k-graded analogue of de Boer–Ducas–Fehr (Eurocrypt 2020, Def. 4, Lemma 5: p-evenly
  distributed samples generate a lattice). de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588, §18; ERH)
  use the same kind of condition. That template's constant escape probability cannot hold for NFS rows.

### 2.4 Fourier identities and the collision barrier [proved; codex rkrank_1, checked]
- E 2^{s − rank R_M} = Σ_w ((1 + μ̂(w))/2)^M and Pr[S_k = S′_k] = 2^{−s}Σ_w μ̂(w)^{2k}.
- One row has support ≤ J = L^{O(1)}, so Pr[S_k = S′_k] ≥ J^{−k}.
- **(B-M):** certifying span by showing that k-fold sums are nearly uniform (collision ≈ 2^{−s+r}) needs
  k ≥ (s − r)log 2/log J = D/L^{o(1)}. Pair-correlation-level information cannot do it *through this
  route*.
- This is not a general lower bound: a law uniform on a basis spans in s steps with poor mixing.

### 2.5 Prop C: algebraic identities and box obstructions [proved]
**Definition.** An **algebraic identity** is a pair (ψ₁, ψ₂) such that:
- ψ₁ is a quadratic Dirichlet character, of any conductor;
- ψ₂ is a quadratic character of (O_K/𝔮)^×;
- ψ₁(a − bm)ψ₂(a − bα) = 1 for all coprime (a, b) with units as arguments.

Every dual functional on the smooth groups is the restriction of such a pair (Kummer plus Chebotarev),
generally of enormous conductor. Prop C says which pairs hold identically.

**Prop C.** Let ℓ be a prime with ℓ > 4(d+1)² and ℓ ∤ 2f_d·disc(f). Then the ℓ-component of an
algebraic identity is trivial unless ℓ | N. At ℓ | N it is a power of (χ_ℓ, χ_{𝔓_ℓ}),
where 𝔓_ℓ = (ℓ, α − m).
- The identities therefore span ≤ ω + O_d(#{ℓ ≤ 4(d+1)²} + ω(f_d disc f)) = ω + O_d(log N) dimensions.
- There are none at the archimedean places on all of Z².

*Proof (fixes from the reviews folded in).*
1. **CRT and lifting.** The identity is a function of (a, b) mod ∏ℓ^{e_ℓ}. Every primitive local tuple
   lifts to a globally coprime integer pair: choose a ≡ a_ℓ, b ≡ b_ℓ by CRT, then adjust b by a multiple
   of ∏ℓ^{e_ℓ} to make gcd(a, b) = 1 (Dirichlet). Setting the other components to (1, 0), where every
   factor is 1, shows each local component is ≡ 1.
2. **Local form.** For odd ℓ ∤ f_d disc f, the ℓ-component is χ_ℓ(G(a, b)), where
   G = (a − bm)^{ε₀}∏_{𝔩|ℓ}N_𝔩(a − bα)^{ε_𝔩} is a binary form over F_ℓ, using (x/𝔩) = χ_ℓ(N(x)). The
   factors N_𝔩(a − bα) are the distinct irreducible factors of the homogenised f mod ℓ. They are coprime
   to a − bm unless a degree-one 𝔩 has root m, i.e. ℓ | f(m) = N.
3. **Odd degree.** Replace (a, b) by (λa, λb) with λ a non-residue. This is again an admissible local
   pair, and χ_ℓ(G) changes sign at every unit point, contradicting the identity. No Weil bound is needed.
4. **Even degree, ℓ ∤ N, G ≠ const.** G is squarefree and non-square. Up to d + 1 points of P¹(F_ℓ) are
   excluded, namely the zeros of every argument required to be a unit. If the identity held,
   Σ_{P¹} χ_ℓ(G) ≥ ℓ − 2d − 1. Weil gives |Σ| ≤ (deg G − 1)√ℓ + 1 ≤ d√ℓ + 1. These are incompatible for
   ℓ > 4(d+1)², since then √ℓ(√ℓ − d) > 2(d+1)(d+2) > 2d + 2.
5. **ℓ | N.** Apply step 4 to G/(a − bm)². Only ε₀ = ε_{𝔓_ℓ} survives. (If ℓ | N and ℓ | disc f, then
   gcd(N, disc f) factors N.)
6. Each excluded prime contributes O_d(1) dimensions: ≤ d + 2 for odd ℓ, ≤ 2d + 2 at ℓ = 2 because of
   principal units. ∎

**Box-level obstructions [proved].** These are identities on the box only.
- **Non-monic f (LV).** The |S|-parity column is constant 1 and a − bm < 0 on the box, so
  e_sign + e_{|S|} ∈ R^⊥ and W ⊆ {u_sign = 0} always.
- **Real places.** A real-sign column that is constant **1** on the rows gives e_sign + e_{rs_j} ∈ R^⊥.
  A column constant 0 gives nothing for the rational side.
- Total: ≤ r₁ + 1 dimensions. Design V uses positive targets, or pairs of negative ones.
- **Correction to hypD §5.** hypD attributed a sign obstruction to a "real root near m". In LV the sign
  obstruction comes from the |S| column regardless.
  - A real root within m/4 of m requires H ≥ c_d·m. Diagnostic only, not a proof: with degree-9 LV-style
    fields and R-height m^η, it occurred in 0/20 trials for η ≤ 1 and 20/20 for η ≥ 1.5. This is
    consistent with a displacement N/f′(m) ≈ m/(1 + H/m).
- **Monic toy (E1).** There is no |S| column. P4 found no sign obstruction in any record, and every
  constant real column was constant 0.

**Scope.** Prop C classifies identities on all of Z². It gives no bound on obstructions that hold on the
finite box or on X, and that is exactly what RK is about.

### 2.6 Lemma S: the rational model satisfies Exp_T [proved; switching from codex rkrank_1]
**Lemma S.** Let n be uniform among B-smooth integers in [1, x], with x = N^{c+o(1)} for fixed c > 0 and
B = L[1/2, a] for fixed a > 0. Let V be a k-dimensional space of prime-parity functionals. Then
  Pr[n ∉ V^⊥] ≥ kΨ(x/B, B)/((1 + log₂x)Ψ(x, B)) ≥ k/(D·L^{o(1)}),
uniformly in V and k.

*Proof.*
1. Pick primes p₁..p_k ≤ B whose evaluation vectors form a basis of V*.
2. For smooth n ≤ x/B, the integers n, p₁n, …, p_kn have distinct evaluation vectors, so at most one of
   them is annihilated.
3. Each y ≤ x has ≤ 1 + ω(y) such preimages.
4. Canfield–Erdős–Pomerance gives Ψ(x/B, B)/Ψ(x, B) = B^{−1}L^{−o(1)}. ∎

This is in the spirit of the standard smooth-number switching identities (cf. Croot–Granville–Pemantle–
Tetali 2012). We found no statement of this F₂ escape bound.

**The joint law.**
- Nothing guarantees a bounded-congestion move (a, b) ↦ (a′, b′) with a′ − b′m = p(a − bm) that stays in
  the box and keeps F smooth. Such pairs may exist by accident.
- F(a, b) ≡ Nb^d (mod a − bm) couples the two sides.
- **A bounded-congestion switching for the joint kept law would prove Exp_T, and hence RK.** By Theorem
  V_κ, it would suffice to have it for dual spaces of dimension > D^{2/3−ε}.

### 2.7 Barriers [proved for the stated templates]
- **(B-F) Naive field averaging.**
  - With ≤ N^{O(1)} fields drawn uniformly, Pr_A[w ∈ W_A^⊥] is 0 or ≥ N^{−O(1)} for each rational
    functional w. A per-functional union bound over 2^{D_rat} functionals therefore needs all but
    2^r·N^{O(1)} of them to be obstructions in no field.
  - Spread example (codex rkrank_1): rows of field A are uniform on E_A^⊥ for a spread {E_A} of
    D/2-dimensional subspaces. Every fixed nonzero functional is an obstruction in ≤ 1 field, yet every
    field has deficiency D/2.
  - This kills the template. It does **not** show that averaging cannot prove RK. Aggregate counts,
    classification, or extra sampling randomness can evade it.
- **Many fields.** Congruences from different fields multiply mod N; only the algebraic square root must
  stay within one field. So W = Σ_i W_i and W^⊥ = ∩W_i^⊥.
  - L^{o(1)} fields help under a uniform per-field contraction of the surviving dual.
  - Fixed-functional survival bounds alone could need many more fields. That is not shown to be
    necessary.
- **(B-M)** is in §2.4.
- **(B-C) Coverage does not imply span.** Rows (e_i, e_i) cover every coordinate arbitrarily often yet give
  W = 0 (codex).
- **(B-L) Random exponents [assessment].**
  - Hafner–McCurley-type generation needs relation vectors that are quasi-uniform on a coset given the
    observable. One NFS row has support L^{O(1)}, which is much smaller than |Λ/2|. Compact random-exponent
    constructions are not excluded.
  - A twisted Vallée sampler (c_e x² mod N small, c_e a product of k factor-base primes) plausibly
    satisfies Exp_T [sketch, unproved]. But each of its rows costs a sampler smoothness test, so spanning
    with it is Vallée's own c = √(2θ).

### 2.8 Theorem V_κ: weaker RK suffices [proved implication]
**Assumptions.**
- Design V's inputs: LV supply at fixed odd d, LV §6, and a sampler of exponent θ with L^{o(1)} cost per
  call and pointwise ratio K_V.
- RK weakened to: raw deficiency on U ≤ D^κ = L^{κa+o(1)}.
- Supply, rank and postprocessing hold **jointly** on a set of fields of measure L^{−o(1)} under the
  law actually reached.
- **Sampler rows must be supported in U.** If U excludes poorly covered primes, the sampler is
  conditioned on avoiding them. This costs L^{−o(1)} provided the excluded primes have bounded
  reciprocal-sum weight; that is an assumption.

**Algorithm.** Collect D^κ + O(log N) + 1 sampler rows (§2.2). The costs, with η = 3/(d+1), are:
- linear algebra: L^{2a};
- LV rows: L^{a+η/(2a)};
- sampler rows: L^{κa+θ/(2a)}.

**Optimum.**
- c(κ) = √(2 max{η, θ/(2−κ)}), at a = √(θ/(2(2−κ))). For κ ≤ 1 this crossing is the minimiser, since the
  sampler term decreases in a.
- κ = 0 recovers Theorem V's √θ. κ = 1 is Vallée/Dixon √(2θ): D + 1 rows are always dependent.
- **c < 1 iff κ < 2 − 2θ** (given η ≤ θ/(2−κ)). For θ = 2/3 this means κ < 2/3, with d ≥ 9 for every
  κ ∈ [0, 1].
- The κ = 0 shape is the classical index-calculus split between precomputation and per-target work
  (Pomerance 1987; Enge–Gaudry 2002). We found no prior use of "codimension D^κ" as a weakened span
  hypothesis.

## 3. Lever B: small-square samplers

### 3.1 Prop Q: consequences of a better sampler [proved implications]
**Definition.** A sampler of exponent θ:
- outputs x with |x² mod N| ≤ N^{θ+o(1)};
- satisfies Pr[εx] ≥ Pr[x]/K_V for all ε ∈ μ₂;
- outputs z that is B-smooth with probability ≥ ρ(log h/log B)L^{−o(1)}.

**Inputs.**
- Fruitfulness ≥ 1/(2K_V) for root-blind dependency selection.
- Smooth-square counting via products of two same-class smooth numbers ≤ √h (Dixon 1981; Pomerance
  1982). Restricting to Ω ≤ K = o(log L) keeps the multiplicities 2^{2K} = L^{o(1)}.
- Trial division to exp(√log N), which makes 2^ω = L^{o(1)}.
- Bernstein's batch smooth-part algorithm (2004 draft) for L^{o(1)} amortised smoothness tests.

**Costs.**
- **Dixon-style:** min_a max{2a, a + s + θ/(2a)} = **s + √(2θ)**, at a = √(θ/2). It is < 1 iff
  s + √(2θ) < 1.
  - This is codex rkrank_2's formula. Our first draft's "correction" s + √(s² + 2θ) was wrong: it balanced
    the two terms past the minimum of the collection term. Proof critic and codex rkrank_3 both caught it.
- **Design V:** min_a max{2a, a + η/(2a), s + θ/(2a)} = max{√(2η), (s + √(s² + 4θ))/2}. It is < 1 iff
  θ + s < 1 when η does not bind.
- At s = 0 these give √(2θ) and √θ. hypD's "θ_V < 2/3 ⇒ √θ_V" is correct, given RK, supply and
  d + 1 ≥ 6/θ.
- That θ < 1/2 alone beats c = 1 is immediate and folklore-level: class-group methods are √|Δ|-size
  objects.

### 3.2 Why 2/3, and Theorem L [proved for a restricted model]
**Mechanism (Vallée; reconstructed by us and by codex rkrank_2).**
- (x₁ + t)² ≡ x₁² + 2x₁t + t².
- Absorbing t² ≤ h leaves an affine lattice of determinant N in a T × h box with T ≤ √h.
- Productive cells need Th ≳ N, i.e. h ≳ N^{2/3}.
- Vallée's Farey dissection controls skewed cells exactly down to this threshold (Vallée §5.1, per
  codex).

**Rank-collapse lemma [proved; folklore-level linear algebra].** Setup:
- w ∈ (Z/N)^r and the symmetric box 𝔅 = ∏[−T_i, T_i] with T_i ≥ 1;
- s_ij = w_iw_j mod N, centred;
- q(t) = Σs_ij t_it_j;
- assume |q| ≤ h on 𝔅 and (5/4)h² < N, which forces **θ < 1/2**.

Then, unless a gcd with N factors N, w ≡ v·u with u ∈ Z^r primitive, v² ≡ σ ∈ Z \ {0}, and
q(t) = σ(u·t)². So the cell x₁ + w·𝔅 has ≤ 2√h + 1 distinct points.

*Proof.*
1. 4s_ij = q(e_i + e_j) − q(e_i − e_j), so |s_ij| ≤ h/2, and |s_ii| ≤ h.
2. Every 2×2 minor is ≡ 0 (mod N), because S ≡ ww^T, and has absolute value ≤ (5/4)h² < N. So every
   minor vanishes and S = σuu^T.
3. Choose i with u_i, w_i units, and put v = w_i/u_i. Then w_j ≡ vu_j.
4. τ_max = Σ|u_i|T_i is attained, so |σ|τ_max² ≤ h. ∎

The hypothesis is needed: for √N < h, e.g. w = (1, ⌈√N⌉), integral forms with nonzero minors kN exist.
This particular rank-two cell has fewer points than a rank-one cell, but no general bound is proved for
1/2 ≤ θ < 2/3 (codex rkrank_3 and the proof critic).

**Theorem L [proved].** Setup:
- the sampler picks a random cell 𝒞;
- (i) either 𝒞 is rank-one, x₁ + vτ with v² ≡ σ and |σ|τ² ≤ h (so |𝒞| ≤ 2√h + 1), for any θ < 2/3;
- (ii) or 𝒞 has arbitrary rank under the lemma's hypotheses, for θ < 1/2;
- obliviousness: every y ∈ Z/N has Pr[y ∈ 𝒞] ≤ λ(2√h+1)/N.

Examples of obliviousness: a centre with density ≤ λ/N and directions independent of the centre, or a
partition into cells selected by a point of density ≤ λ/N. Vallée selects her Farey cells this way; her
cells themselves are curved, not linearisable.

Then E#(𝒞 ∩ S_h) ≤ λ·2^ω(2h+1)(2√h+1)/N, where S_h is the set of units x with |x² mod N| ≤ h. Hence
N^{1−3θ/2−o(1)} cells are needed per output when λ, 2^ω = L^{o(1)}.

*Proof.* E# = Σ_{y∈S_h}Pr[y ∈ 𝒞] and |S_h| ≤ 2^ω(2h + 1) for squarefree N. ∎

*Scope.*
- Without obliviousness, cells can be aimed at known small squares. Codex rkrank_3 gives an explicit
  counterexample to the first draft's centre-only density hypothesis.
- The theorem does not cover arithmetic choice of productive cells, exact curved regions, or general
  lattices.
- No unconditional bound θ > 0 is possible without a factoring lower bound: polynomial-time factoring
  gives θ = 0.
- Codex rkrank_2 proves the interval-cell special case.

**Coppersmith [proved for the construction].**
- Take the shift lattice N^{m−j}t^iF^j, 0 ≤ i ≤ 2(m − j), for F = (x₁ + t)² − z.
- Its dimension is (m+1)² and det = (NT)^A h^C, with A = m(m+1)(4m+5)/6 and C = m(m+1)(2m+1)/6.
- The determinant (Howgrave-Graham/LLL) sufficient condition is T^{(4m+5)/(2m+1)}h < N, which tends to
  T²h < N.
- This is incompatible with Th ≥ N^{1−o(1)}. So the standard guarantee is unavailable. That does not
  prove the lattice useless.
- This is a routine Jochemsz–May count, derived independently by us and by codex rkrank_2.

### 3.3 Other routes [assessment unless marked]
- **Ratio forms** x² ≡ u/v are equivalent to (vx)² ≡ uv. [proved]
- **Class groups (LP, Seysen, Hafner–McCurley).** These are θ = 1/2 objects: reduced forms have
  A ≤ √(|Δ|/3).
  - [proved, codex rkrank_2] Ideals of norm ≤ |Δ|^θ with θ < 1/2 reach only O(A log A) of
    |Δ|^{1/2+o(1)} classes, so uniform class sampling cannot go below θ = 1/2.
- **Elkies (ANTS 2000):** tangent subdivision reproduces the curvature threshold.
- **Schnorr (ePrint 2021/933), Ajtai–Kumar–Sivakumar, Boneh's CRT decoding:** none supplies rigorous
  smoothness, rank and fibre symmetry together.
- **Hittmeir, Harvey, Harvey–Hittmeir:** rigorous but exponential.
- **Square identities.** If P(t)² ≡ Q(t)² identically, a fibre-symmetric family factors N directly.

### 3.4 Literature (literature subagent, codex rkrank_1–3, literature critic)
- **Vallée.**
  - Math. Comp. 56 (1991), no. 194, 823–849; STOC 1989.
  - Target |x² mod N| ≤ 4N^{2/3}, Farey order ≈ N^{1/3}/4.
  - Secondary source (Gentry, US patent 7,957,525): quasi-uniformity constants l₁ = 1/5, l₂ = 4.
  - We could not open the primary PDF. Codex read it via ResearchGate (rkrank_2, hypD_2).
  - LP 1992 cites the resulting L[1/2, √(4/3)].
- **Coron**, CRYPTO 2002, LNCS 2442, 613–626: width O(N^{2/3+ε}), statistical distance ≈ N^{−3ε/13}
  (secondary source).
- **Gentry**, CRYPTO 2004, LNCS 3152, 179–200 (§4): exact uniformity at ≈ 8N^{2/3}, per codex. Not
  confirmed by the literature critic; treat as unverified.
- **Stange**, arXiv 2211.06821v2.
  - Hyp. 3.1 (heuristic): saturation of a one-dimensional restriction of an integer relation lattice.
  - Thm 3.2 is relative to a relation oracle; Thm 2.1 uses ERH.
  - NFS appears in the abstract, the introduction and remark (9). There is no GF(2) statement.
  - "1 − 1/ζ(c+1)" is consistently misstated (six occurrences) for 1/ζ(c+1) (codex; literature critic).
- **Lattice generation.**
  - Hafner–McCurley, JAMS 2 (1989) 837–850 (ERH; √2).
  - Biasse–Erukulangara 2021 (3/√8, GRH).
  - de Boer–Ducas–Fehr, Eurocrypt 2020, Lemma 5.
  - de Boer–Pellet-Mary–Wesolowski, arXiv 2512.01588 (ERH).
  - Enge–Gaudry, Acta Arith. 2002.
- **Squares and rank.**
  - Croot–Granville–Pemantle–Tetali, Annals 175 (2012).
  - Balister–Bollobás–Morris, Annals 188 (2018): first-dependency threshold, with the 2-core, for uniform
    integers. E1's 2-core peeling follows their language.
  - Coja-Oghlan et al. and Cooper–Frieze–Pegden: random sparse matrix models.
  - **No GF(2) span or saturation theorem for sieve relations was found.**
- **Rigorous c < 1:** none. LP 1992 remains the benchmark (Harvey 2020 calls it the best known rigorous
  probabilistic bound).
  - Lee–Venkatesan (J. Number Theory 187 (2018)) prove only congruences of squares unconditionally.
  - Citers checked: 115 of LP 1992 and the citers of LV (literature subagent; counts not re-verified).

## 4. What would close the gap
1. **A bounded-congestion switching for the joint kept law of one LV field.** That is Exp_T(G′_alg) with
   T = D·L^{o(1)}. By V_κ it is needed only for dual spaces of dimension > D^{2/3−ε}.
   - A possible lever, not worked out: L^{o(1)} fields sharing m (common rational parts) with independent
     R. Congruences from different fields can be multiplied (§2.7), so switching a pair *between*
     fields is not ruled out by the bookkeeping.
2. **A non-oblivious, fibre-symmetric sampler below 2/3**, or Theorem L for rank ≥ 2 cells at
   1/2 ≤ θ < 2/3. That case is genuinely open, not just unproved.

## 5. Experiment E1 (pre-registered; exact; toy)
- Commands are in `notes/rkrank_prereg.md`. Analysis:
  `.venv/bin/python -m src.rkrank.analyze experiments/raw/rkrank_e1.jsonl`.
- 20 instances (blum/tri × 64/80 bits × seeds 0–4), 72 level-records.
- Runtime: 15–40 s at 64 bits, 5–13 min at 80 bits. 0 accidental gcd events.

**Disclosures.**
1. *Limitation (as registered).* Boxes ran out before 3.0·cols in 16 of 20 instances, and unreached
   levels are skipped.
   - Levels reached: 3.0 in 4 instances, 2.0 in 8, 1.5 in 4, 1.25 in 4.
   - So the "levels ≥ 2.0" predictions rest on 16 records.
2. *Deviation (found by codex rkrank_3 and the proof critic).*
   - The frozen harness computed the P5 functionals only at level 3.0. The registration says "top level
     reached", so the frozen run tested P5 in only 4/20 instances.
   - The supplementary script `src/rkrank/p5supp.py` regenerates the identical relation sets and
     evaluates the same 13 functionals at the top level reached:
     `.venv/bin/python -m src.rkrank.p5supp experiments/configs/rkrank_e1.json experiments/raw/rkrank_e1_p5supp.jsonl 10`.
   - It draws its random dense functionals from the same RNG stream as the frozen harness, so in
     the instances that reached 3.0 they coincide (checked: 4/4 identical).
   - It also checks the independence assumption behind excess_peel.

| prediction | registered | result |
|---|---|---|
| P0 W ⊆ Λ/2 | all records | 72/72 — **holds** |
| P1a excess_peel ≤ 4 (levels ≥ 1.5) | all records | 30/32 — **fails** (5 in tri 64 s1 and tri 80 s1, both at 1.5) |
| P1b excess_peel = 0 (levels ≥ 2.0) | ≥ 80% | 12/16 = 75% — **fails narrowly** (others: 1, 1, 2 at 2.0; 1 at 3.0) |
| P2a c₃ = 0 (levels ≥ 1.5); c₃ ≤ 2 | ≥ 90%; all | 32/32 — **holds** |
| P2b c₁ = 0 (levels ≥ 2.0) | ≥ 80% | 12/16 = 75% — **fails narrowly** |
| P3 small squares on P₃ in W (levels ≥ 1.5) | ≥ 0.95 | 1141/1141 observations, 514/514 with z < 0 — **holds** (see caveat) |
| P4 sign obstruction ⟺ constant-1 real column | all records | 72/72 — **holds** (vacuously: no sign obstruction and no constant-1 column occurred) |
| P5 δ ∈ [0.45, 0.55], 13 functionals, top level reached | all instances | frozen run: 4/4 tested; supplementary: 20/20 at the top level reached (all 13 values in [0.468, 0.525]) — **holds** |

**Caveats on reading the table.**
- c₃ ≠ 0 in five records at lower levels: four at 1.0 (blum 64 s1, blum 80 s3, tri 80 s1, tri 80 s3)
  and one at 1.25 (tri 80 s3).
- For z > 0, P3 is implied by c₃ = 0. For z < 0 it needs one W-vector with the sign bit set and support
  in P₃. The 1141 observations reuse 1010 sampled vectors across nested prefixes, so they are not
  independent trials.
- P5 tests 13 fixed functionals. Random dense functionals are balanced almost regardless of structure,
  so P5 shows only that *these* are balanced. Exact dense obstructions would already appear in dim W^⊥;
  adaptive near-obstructions are not tested.
- The independence check (supplementary): rank(Gen ∪ {e_p : p peeled}) = ω + n_peel in 72/72 level-records.
- c₁ = dim W^⊥ − ω − #{occ = 0} in every record. This is close to an identity, not a finding.
- Descriptive: c₂ = 0 in 12/16 records at 1.5, 11/12 at 2.0 and 4/4 at 3.0. The excess over the peeling
  predictor falls with M: 4–23 at level 1.0, 0–5 at 1.5, 0–2 at 2.0, 0–1 at 3.0.

**Reading [experimentally supported; toy only].**
- In these monic cubic toys, at ≥ 1.5·cols, the observed dual obstructions are the ω genuine characters
  plus functionals involving primes that occur at most twice. Mostly these are primes outside the
  algebraic 2-core.
- No obstruction involving only well-covered primes was seen (P2a), and no sign obstruction arose (P4).
- The failed sharp predictions (P1, P2b) failed through ≤ 5 rare-prime obstructions that escape the
  simple peeling predictor.
- This is the sparse form of RK at toy scale: raw matrix, cubic, small B. It says nothing about d ≥ 9,
  L[1/2]-size factor bases, the LV coefficient law, or genuine certificates.

## 6. Candidate cards
```text
ID: RK-U  uniform rank growth (Theorem U) for one LV field's smooth set.
Proved: Theorem U (Exp_T(G') => deficiency <= r + dim G' after T(1+ln s) relations; routine);
  Lemma S (rational model satisfies Exp_T with T = D L^{o(1)}); Prop G' (raw vs genuine <= O(log N);
  classical, BLP §8); Prop C (algebraic identities on all of Z^2 = Gen + O_d(log N); box: <= r1+1,
  automatic sign obstruction for non-monic f); Theorem V_kappa (deficiency D^kappa, kappa < 2-2theta, suffices).
New lemma required: Exp_T(G'_alg) for the joint kept law (bounded-congestion joint switching), needed only
  for dual spaces of dimension > D^{2/3-eps}.
Barriers (templates only): naive field averaging (B-F), collision certificates (B-M), coverage != span (B-C).
Status: OPEN. E1 toy: obstructions = Gen + rare-prime functionals; none on well-covered primes at >= 1.5 cols.

ID: VS-L  small-square samplers below 2/3.
Proved: Prop Q (Dixon s+sqrt(2 theta); Design V max{sqrt(2 eta),(s+sqrt(s^2+4theta))/2});
  rank-collapse lemma (theta<1/2); Theorem L (oblivious rank-one cells: theta >= 2/3 - o(1);
  any rank: theta >= 1/2 - o(1)); Coppersmith determinant guarantee needs T^2 h < N.
Status: refuted for this model (oblivious rank-one linearisable cells; oblivious any-rank cells below 1/2;
  natural Coppersmith guarantee); open for rank >= 2 at 1/2 <= theta < 2/3 and for arithmetic cell choice.
```

## 7. Reviews and adjudication
**GPT-6.1-sol blind rounds.**
- rkrank_1 (A) contributed: the duality, Prop G′ (rederived), Lemma S, the Fourier identities and the
  collision barrier, the spread example, the (e_i, e_i) example, the check of the V_κ formula, the Stange
  ζ issue, and the "bounded-congestion switching" route.
- rkrank_2 (B) contributed: Prop Q (with the correct Dixon formula), the Vallée mechanism, the
  Coppersmith exact determinant, the interval-cell barrier, class-group coverage, and the Coron/Gentry
  literature.

**Adversarial round codex rkrank_3 and two Opus critics (proof; literature/novelty).** All findings were
accepted.
- **WRONG, fixed.**
  - The Dixon "correction" (codex and proof critic): restored to s + √(2θ).
  - Theorem L with a centre-only density hypothesis: codex counterexample with N = 1009·1013, h = 100.
    Restated with a per-point coverage hypothesis and split into (i)/(ii).
  - The rank-collapse lemma: it checked only principal minors, used the wrong bound, and is valid only
    for θ < 1/2. Repaired via q(e_i ± e_j), with its scope stated.
  - E1 reporting: P5 4/4 rather than 20/20 (disclosed deviation plus supplementary run); the c₃
    statements; P3 dependence; c₁ near-identity.
- **CORRECT WITH FIXES, applied.**
  - Theorem U now has a G′ fixed in advance. With G′ = Ann(μ) it says nothing about RK.
  - Prop C: CRT lifting, odd degree by non-residue scaling, excluded points in the Weil count, ℓ = 2 dims,
    and the "accident of X" wording withdrawn. Its scope is all of Z².
  - Box obstructions: constant-1 vs constant-0; the automatic |S| sign obstruction in LV; the "root near
    m" claim relabelled as a diagnostic.
  - Lemma S: parameter range, and "no guaranteed move" wording.
  - B-F/B-M/B-L narrowed to templates. B-L's twisted-sampler claim downgraded to a sketch.
  - "Dependencies must share a field" was wrong: congruences multiply.
  - V_κ: sampler support in U and joint field conditions.
  - Coppersmith: "guarantee unavailable", not "fails".
  - Prop G′ bookkeeping marked as asserted.
  - κ notation split (K_V).
- **Re-crediting (literature critic).**
  - Prop G′ to Buhler–Lenstra–Pomerance §8.
  - The square-class product trick to Dixon/Pomerance.
  - Theorem U and the rank-collapse lemma: routine.
  - V_κ's κ = 0 shape: the index-calculus split.
  - Bernstein cited as a draft.
  - Gentry "exactly uniform" marked unverified.
  - Stange: NFS in the abstract; ζ misstatement systematic.
  - CGPT/BBM credited for the 2-core.
- Novelty statements are limited to "no prior formulation found by us": V_κ's weakening, Prop C's
  packaging, Lemma S's F₂ bound, Theorem L's model.
