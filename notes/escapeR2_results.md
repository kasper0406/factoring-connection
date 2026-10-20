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
re-splitting and factor-then-descend. Each one fails at an identifiable step: it either cannot remove an algebraic
child at all, or it hands the problem back in an equivalent form. What we can now say is narrower than before:
- the only way to remove an algebraic child is inside its own field;
- inside a field, the whole difficulty is one analytic statement: a lower bound for the number of field coefficient
  vectors at which **more values of one polynomial than it has coefficients** are simultaneously smooth;
- in the size model, every polynomial excess helps. If such joint control held for d^{C} points with any fixed
  C > 1, the exponent would drop below 2/5. C = 3/2 gives 3/8, and C → ∞ approaches 1/3.
- after a pooling trick, this statement is needed only up to a factor L^{O(1)}, not to polylogarithmic precision.

So the problem is not in the algebra, the precision bookkeeping or ERH. It is one open estimate about smooth
values of overdetermined linear systems.

**Precisely.**
1. **Lemma TI (cross-field parity transfer) [proved].** Let K_1, …, K_J lie in a Galois field L, and let ∏Γ_j be a square
   in L with Γ_j ∈ K_j. At every prime ℓ unramified in L, the parity vector of Γ_1 at the primes of K_1 above ℓ is
   pulled back from the subfield k_1 := K_1 ∩ K_2⋯K_J. **Corollary EX:** if k_1 = Q (for example, K_1 primitive and not
   contained in the compositum of the other fields), no combination of other fields, rational integers or norms can
   remove an algebraic child of a degree-one element of K_1. The certificate must contain other elements of K_1 lying
   in primes above the same ℓ. This kills families 1(a), 1(b) and 3 in their stated forms (§§2–3, 5). E1 checks it
   exactly on 72 824 group configurations (0 violations).
2. **Tied fields do transfer, but are not fresh [proved example + structural].** If k_1 ≠ Q, transfer happens. For
   example, a prime of Q(√ab) is congruent mod 2 to a sum of primes of Q(√a) and Q(√b) (E1: τ = 1). Such fields are
   subfields of one Galois field L, so they bring no independent randomness (§2.3).
3. **Kummer / odd-isogeny transport (r20 #4) [proved for S_d-fields, d ≥ 7, given fact (F)].** If J′[2] ≅ J[2] (in particular if J′
   is odd-degree isogenous to J = Jac(y² = f)), and Gal(f) = S_d, then the étale algebra of f′ is K or K × Q. The
   transported class lands back in K (Prop KT, §5.2). A different generator of the same field needs a min-polynomial
   of height ≈ N (heuristic count, §5.2), so the transport gives no fresh representation.
4. **Principalisation and class groups (families 1(c), 2) [proved reduction + model].** Removing 𝔭 modulo squares
   requires θ with (θ) = 𝔭·(factor-base)·(square) (codex escapeR2_1 (5.1)). Under ERH such θ always exist (Bach), but
   only at the Arakelov scale ≥ √|Δ_K| (nested Thm AF). There ψ(θ) is a generic residue. Its only known treatment is
   a re-split modulo N, which resets targets to the top scale L[1 − α] (§4, §6). Rigorous class-group algorithms are
   L_{|Δ|}(1/2) (dBPW, ERH). For fields with explicit ψ, |Δ_K| = N^{Θ(1)} heuristically, so this route gives L[1/2].
5. **Re-splitting (family 4(b,c)) [proved obstruction + model].** The square-corrected re-split w ≡ u/(v s²) is valid,
   but Lemma S″ does not extend from the square class to a coset wQ. An adversarial target set gives zero pairs
   (Prop RS; E3: 0 pairs on 5 toy moduli). Any re-split also returns targets at scale ≥ L[1 − α]. Recursion supplies
   factorisations, never relations modulo N (§6).
6. **Budget-respecting in-field descent stays at 2/5 [model] (Prop R3).** Under nested's premises M1–M4 and the
   typical-factorisation model, suppose every algebraic component has O(d) points, which is what the known moment
   method controls. Then the floor is 2/5. This extends Prop R2 from "no algebraic children" to "budget-respecting
   algebraic children" (§7.1).
7. **The budget exponent is the whole story [model; E2 exact].** Allow components of up to d^{C} points (chain
   sub-model). The floor is α*(C) = (3 − C)/(7 − 2C) for 1 ≤ C ≤ 3/2, and α*(C) = (r + 2)/(3r + 5) once C ≥ r + 1 (7/18
   at C = 5/4; 3/8 at C = 3/2; 4/11 at C = 8/3; 5/14 at C = 4; 11/32 at C = 10). It is below 2/5 for every C > 1 and
   above 1/3 for every finite C (§7.2). Codex escapeR2_1 found the C = ∞ level thresholds independently.
8. **Theorem NR-q [proved conditional on Hyp T_q, LV6′ and the multi-field post-processing of lv5unif Prop MF
   (outline)].** Suppose a positive-density set of split primes q (density ≥ (log N)^{−C₀}) each satisfies:
   Pr_A[the in-field descent of q completes] ≥ L[α, −c]. Then factoring runs in L[α, c′]. This is weaker than nested's
   AD in three ways: it is averaged over the field, pointwise in the target, and needs only L^{−O(1)} success (§8). Consequence for r20 #1:
   its two-level moment lemma is needed only up to L^{O(1)}, not (log N)^{O(1)}. Losses of exp(O(#forms)) are
   affordable, so the "kill test" in r20 #1 does not kill it (§7.3).
9. **ERH.** ERH is not the missing input anywhere. It removes character costs (Lemma X_ERH) and makes Bach's generation
   bound available, but supplies neither in-field joint smoothness nor small evaluated images.
10. **Literature (§11).** No rigorous or ERH-conditional L[α < 1/2] factoring was found. No rigorous L(1/3) class-group
    or principal-ideal algorithm was found. No smooth-number result for overdetermined linear systems at L-scale was
    found. Lemma TI is elementary and close to the norm-relation and Brauer-relation literature. We make no novelty claim.

### Per-approach exponent table

α is the size-model exponent the approach would give if its missing lemma were proved. "Refuted" means refuted in
the certificate class of §1 (products of field elements that are squares in a common Galois field), not for every
conceivable algorithm.

| # | approach | α in the model | status | exact missing lemma |
|---|---|---|---|---|
| 1a | replace 𝔭 by its norm ℓ (with Galois conjugates / closure) | none: the complement is ≥ as hard | **refuted** (Cor EX, codex §3) | — |
| 1b | second field K′ sharing m; re-descend ℓ there | 2/5 (it is Prop MF) | cross-field closing **refuted** when K ∩ K′ = Q (TI) | — |
| 1b′ | tied fields inside one Galois L (norm relations) | no gain: the component lives in L, same budget | proved example; structural | budget lemma for L |
| 1c | ψ-image of a generator / S-unit | ≥ 2/5 (collapses to a re-split) | model | contracting image certificate (codex (5.2)) |
| 2 | class-group descent (dBPW, Biasse–Fieker) | 1/2 rigorous (ERH); 1/3 only under the same descent heuristics | KNOWN + model | rigorous L_{|Δ|}(1/3) PIP with controlled ψ-images, or small-|Δ| fields with explicit ψ |
| 3 | Galois/Chebotarev, conjugate families, Kummer transport (r20 #4), Gassmann (r20 #7) | no gain | **refuted** for S_d-fields (TI, KT; r20 #7) | a non-S_d family with transfer *and* freshness (none known) |
| 4a | independent one-level descents; global recombination (r20 #9) | 2/5 (R2, R3) | model; r20 #9 killed under diffuse labels | concentrated-label reservoirs |
| 4b | split law with halves re-split recursively | no contraction (scale reset) | **refuted as stated** (RS) | coset equidistribution + a contracting re-split |
| 4c | factor-then-descend | = special-q descent | structural (recursion gives no relations mod N) | — |
| 5a | in-field multi-level, components O(d) points | 2/5 | model (Prop R3) | — |
| 5b | in-field multi-level, components d^{C} points | α*(C): 3/8 at C = 3/2 → 1/3 | model (E2) | **Hyp LS_C** (joint smooth supply for d^{C}-point configurations) |
| 5c | r20 #1: moments of complete two-level certificates | 3/8 (needs C = 3/2) | open; precision relaxed to L^{O(1)} (§7.3) | Hyp T_q at 3/8 |
| 5d | full multi-level descent | 1/3 | conditional (Thm NR: AD; Thm NR-q: T_q) | Hyp T_q at 1/3 |

---

## 1. Setting and the certificate class

Notation as in nested §0 and lv5unif §1: t = log N, D = t/log t, Λ_α = D^α log t, L[α, c] = e^{(c+o(1))Λ_α}. A *field*
is A = (f, g) with F(g₀, g₁) ≡ 0 (mod N), K = Q(α), ψ: Z[α] → Z/N, α ↦ m = g₀/g₁. A *row* is
r² ≡ a·∏_j ψ_j(Z_j) (mod N) with a ∈ Q^× and Z_j ∈ K_j^× (codex escapeR2_1 (2.1)). A dependency makes a = a₀² and
Z_j ∈ K_j^{×2} for every j. The congruence is then completed by codex's square-root identity (2.5) with
denominator-square corrections. That identity is a direct check and is correct; it is lv5unif §11's multi-field
post-processing written as a formula.

**The most general multi-field square certificate** asks instead that ∏_j Γ_j be a square in a Galois field L
containing all K_j, with ψ extended to O_L. Such an extension exists on Z[α_1, …, α_J] (all α_j ↦ m). It is a ring
homomorphism because Z[α_1,…,α_J] is a quotient of Z[x_1,…,x_J]/(f_1(x_1),…,f_J(x_J)) and every f_j(m) ≡ 0. The
per-field certificate is the special case Γ_j ∈ K_j^{×2}. Lemma TI applies to the general version.

An **algebraic child** of a degree-one element z_P = a − bα (gcd(a, b) = 1) is a prime 𝔭 = (ℓ, a/b mod ℓ) with
ℓ > B, ℓ ∤ b·f_d·[O_K : Z[α]]. To *remove* it means to reach a certificate in which v_𝔭 is even.

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
- Ramified primes and primes dividing indices are finitely many per field. They are exactly the primes that LV6
  places in its residual group H, handled by characters. TI is about large unramified primes.
- The two-field case is sharp: Im_1 ∩ Im_2 is exactly the image of the ideals of K_1 ∩ K_2 (E1 P3: τ₂ = 0 in all
  6 883 two-field cases).
- For three or more fields, the sum of the images can contain more than the separate pairwise intersections; that is
  why k_1 is defined with the compositum (E1(b)). The bound in TI is an inclusion, not an equality.

### 2.2 Corollary EX: export of algebraic children is impossible between unrelated fields

**Corollary EX [proved].** Assume k_1 = K_1 ∩ K_2⋯K_J = Q. Let Γ_1 = ∏ z_i^{e_i} be a product of degree-one elements
z_i = a_i − b_iα_1 (gcd(a_i, b_i) = 1). Let 𝔭 = (ℓ, r) be a prime of degree one with ℓ unramified in L, ℓ ∤ f_d and
ℓ ∤ [O_{K_1}:Z[α_1]], and with v_𝔭(Γ_1) odd. Then:
- (i) every prime of K_1 above ℓ has odd valuation in Γ_1;
- (ii) all primes of K_1 above ℓ have degree one, and each of them contains some z_i with e_i odd;
- (iii) in particular Γ_1 contains, besides the element that opened 𝔭, at least one further z_i lying in a prime above
  ℓ, unless K_1 = Q.

*Proof.* (i) is TI with k_1 = Q: there is one prime of Q above ℓ, so x_1 is constant, and it equals 1 at 𝔭.
(ii) If ℓ ∤ b_i, then z_i lies in at most one prime above ℓ, namely (ℓ, a_i/b_i), which has degree one. If ℓ | b_i,
then ℓ ∤ a_i and z_i is an ℓ-unit (its norm F(a_i, b_i) is ≡ f_d a_i^d ≢ 0 mod ℓ). A prime 𝔓 of degree ≥ 2 above ℓ therefore has
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
bound; the algorithm can test it exactly, since K_1 ⊂ K_2⋯K_J is decidable. Sharing the rational side ties the fields
only through Q.

### 2.3 Tied fields: transfer exists, freshness does not

**Example TB [proved; E1(b)].** In L = Q(√a, √b) with G = V₄, let ℓ split completely, and label the four primes of L
by sign pairs (ε_a, ε_b) of √a, √b mod ℓ. A prime of Q(√a) is a set {(ε, ±)}, a prime of Q(√b) is a set {(±, ε′)},
and a prime of Q(√ab) is a set {(ε, ε′), (−ε, −ε′)}. Then {(ε, ε′), (ε, −ε′)} + {(ε, ε′), (−ε, ε′)} = {(ε, −ε′), (−ε, ε′)} mod 2.
So a prime of Q(√ab) is ≡ a sum of a prime of Q(√a) and a prime of Q(√b). Pairwise there is no transfer (τ₂ = 0);
jointly there is (τ = 1). This is the parity shadow of the Brauer/norm relation 2 = (1 + σ) + (1 + τ) − σ(1 + στ)
in Q[V₄] (§11).

E1's census shows the same thing is common in small groups. Among 65 940 subgroup triples with D ranging over all
cyclic subgroups, τ > 0 in 40 149. In every one of them k_1 ≠ Q: among the 7 948 triples with k_1 = Q, τ > 0 in 0.

*Why this is no escape [structural].* Transfer needs a common Galois field L in which K_1 is not disjoint from the
others. All such fields are subfields of one object L, and none of them is chosen
independently of K_1. The descent problem moves into L. Its children are ideals of L, closed by elements of L, whose
parameter budget is that of L's defining data. This is the same situation as one field. We do not claim more; a
construction with many subfields *and* independent randomness per subfield would be new. The multiquadratic fields
of §4 show the scale problem such constructions meet.

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

*Under ERH existence is free.* Bach's bound: Cl(K) is generated by primes of norm ≤ 12 log²|Δ_K|. Once the factor base
contains those primes (they are below B), every 𝔭 satisfies the condition of Lemma PR. So the question is never
existence. It is (i) finding θ, and (ii) handling ψ(θ).

*(i) Finding θ.* The only rigorous general method is dBPW's sampler (ERH), at cost L_{|Δ_K|}(1/2), with cofactor norm
≥ √|Δ_K| (nested Theorem AF). For fields carrying an explicit ψ, the counting premise M1 gives a polynomial
discriminant ≈ H_f^{2d−2} ≈ N^{2(d−1)/(d+1)}. Whether the field discriminant matches is [heuristic], and nested's
Cor AF premise (i) records that no universal lower bound holds (x^d − 2). Under the smooth-ideal model this costs
L_N[1/2] (codex (6.1)). Biasse's and Biasse–Fieker's L(1/3) class-group algorithms rest on smoothness and
relation-generation heuristics, which are the descent heuristics again (§11).

*(ii) ψ(θ).* θ is a general element of the ideal lattice, not a degree-one element. ψ(θ) = Σ θ_i m^i mod N has no
small rational representation. It is a generic residue. The only known way to write it over smooth integers modulo N
is a re-split (§6). That returns targets of size ≥ L[1 − α] in an arbitrary square class (Prop RS), so this route
collapses to the top of the descent [model]. Small degree-one θ would be an in-field closing, i.e. the original
problem.

*Small-discriminant fields (codex row 4).* If fields with log|Δ_K| = O(Λ_{2/3}) and explicit ψ existed for every N,
then L_{|Δ|}(1/2) would be L_N[1/3]. But an explicit ψ on a field of small discriminant means an explicit root modulo N
of a polynomial of small height. For a fixed small field that is a root-finding problem modulo a composite: for
x² − c it is equivalent to factoring (Rabin). The NFS avoids this only by building f from N, which forces large height.
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

**Fact (F) [KNOWN; Dixon–Mortimer, *Permutation Groups*, §5.2 (Thm 5.2B); flagged for the literature critic].** For
d ≥ 7, every subgroup of S_d of index ≤ d + 1 is S_d, A_d or a point stabiliser.

**Proposition KT [proved, given (F)].** Let f ∈ Z[x] be squarefree of odd degree d ≥ 7 with Galois group S_d,
K = Q[x]/(f), J = Jac(y² = f). Let f′ be squarefree of degree d or d + 1, J′ = Jac(y² = f′) of the same genus, and
suppose J′[2] ≅ J[2] as Gal(Q̄/Q)-modules. This holds whenever J′ and J are related by an isogeny of odd degree
defined over Q. Then Q[x]/(f′) ≅ K if deg f′ = d, and ≅ K × Q if deg f′ = d + 1.

*Proof.* Let X, X′ be the root sets. For odd |X|, J[2] ≅ F₂[X]⁰ (sum-zero vectors), and for even |X′|,
J′[2] ≅ F₂[X′]⁰/⟨𝟙⟩ (Schaefer 1995 §2; Stoll 2001). An odd-degree isogeny has odd-order kernel, so it is injective on
2-torsion, hence an isomorphism J[2] → J′[2] of Galois modules. The field Q(J[2]) is the splitting field L of f,
because S_X acts faithfully on F₂[X]⁰ for |X| ≥ 3. Likewise Q(J′[2]) is the splitting field of f′, because S_{X′}
acts faithfully on F₂[X′]⁰/⟨𝟙⟩ for |X′| ≥ 5: a transposition (12) sends e_1 + e_3 to e_2 + e_3, and the difference
e_1 + e_2 is not in ⟨𝟙⟩. So both splitting fields equal L, and G = S_d acts faithfully on X′ with |X′| ≤ d + 1. By
(F), each orbit has size 1, 2 or d, since 2d > d + 1. If no orbit has size d, then A_d fixes X′ pointwise,
contradicting faithfulness. So there is exactly one orbit of size d, with point stabiliser S_{d−1}, and it is
isomorphic to X as a G-set (d ≠ 6, all index-d subgroups are conjugate). The remainder has size 0 or 1. The étale
algebra of f′ is therefore K, or K × Q. ∎

So any transport that preserves 2-torsion lands in a conjugate of K itself. Within K, a "new representation" means
a different generator α′ with ψ(α′) ≡ g′₀/g′₁. Rigorously |F′(g′₀, g′₁)| ≥ N, so H_{f′}·H_{g′}^d ≥ N/(d + 1), the
same constraint as for α. Heuristically, generators of O_K with min-polynomial height ≤ H number about H/√|Δ_K|, so
the next one after α has height ≈ √|Δ_K| ≈ N^{(d−1)/(d+1)}. Each also needs ψ(α′) to be a small fraction modulo N,
with probability ≈ H_{g′}²/N. [heuristic count.] So the transported representation has Arakelov-scale heights.

*Scope.* KT needs Gal(f) = S_d and d ≥ 7. For exotic Galois groups with non-conjugate subfields carrying isomorphic
2-torsion modules (Gassmann-type), r20 #7 shows that transport mapping embedding vectors to embedding vectors forces
isomorphic fields. Mixing transports remain conceivable, but such fields are tied (§2.3) and hence not fresh. We do
not rule out a cohomological transport that is not induced by a 2-torsion isomorphism. r20 #4 must then name its
mechanism.

### 5.3 Gassmann / torus transport (r20 #7)

We checked r20's kill and agree. Over a common Galois closure the character lattices are permutation lattices, and an
odd-degree isogeny is invertible mod 2. A basis vector mapping to a basis vector propagates along the orbit to a
G-set bijection, hence isomorphic fields. This is a special case of the phenomenon TI captures (transfer only through
tied structure).

---

## 6. Family 4: re-split, factor-then-descend, independent descents

**Prop RS (re-splitting modulo N).**
- (i) [proved; codex (10.1)] For a unit w and fresh s, R := ws² and a Gauss-reduced u ≡ Rv give w ≡ u/(v s²). The
  original target is preserved with an explicit square correction.
- (ii) [proved; E3] Lemma S″'s lower bound Pr ≥ |S|²/(2φ(N)) holds for R uniform on the squares Q_N. For R uniform on
  a coset wQ_N, the relevant pair count is Σ_c n_c n_{wc} (codex (10.2)), and it can be 0. Take the target set
  G ⊂ {primes that are squares mod every p | N}; then every product of G-primes is in Q_N, and for w ∉ Q_N no pair
  exists. Such G have relative density 2^{−ω(N)} among primes (Chebotarev), which can exceed (log N)^{−C} when ω is
  small. E3 (5 seeds, 48-bit N): 0 pairs and 0 hits in 2 000 draws for non-square w, against 263 040–640 312 pairs for
  w = 1 on the same S. The control G = all window primes gives 4.02·10⁷ pairs ≥ |S|²/(2h) = 2.06·10⁷.
- (iii) [model] Every re-split is a 2-dimensional Gauss split modulo N, so its values are ≈ √N whatever the incoming
  target size. Targets return to scale Y with u-cost log(1/Pr) ≈ (t/(2 log Y))·log(t/log Y), so L[α] cost forces
  Y ≥ L[1 − α]. A re-split never contracts.
- *Missing lemma:* class equidistribution of k-products of the descendable primes in arbitrary cosets of Q_N, and a
  contracting replacement for (iii). Neither is in sight. (i) and (ii) were found independently by codex.

**Recursion cannot replace descent [proved, structural].** The recursive call returns the complete factorisation of a
public integer M < N. A congruence modulo N needs a multiplicative relation x ≡ ∏p_i^{e_i} (mod N). A factorisation of
M gives one only if M ≡ x (mod N) by construction, and that is exactly a special-q relation (M = G(P), x = g₁ψ(z_P)).
Factoring M = ℓ + jN directly would need M ≥ N, which recursion does not cover. Codex §10 also observes that a square
congruence modulo M says nothing modulo N. So 4(c) is the special-q design, with factoring used only for testing,
and 4(b)'s "halves descended by recursion" means re-splits, i.e. (iii).

**Independent one-level descents [model].** If every descent closes fully at one level, combining them changes
constants only (R2, R3). Global recombination of partial rows needs repeated large labels. Under diffuse labels the
expected number of repeats is ≪ M²/V = o(1) (codex §8, r20 #9). We agree with r20's kill under diffuse statistics.

---

## 7. In-field descent: the budget is the whole problem

### 7.1 Prop R3: budget-respecting components stay at 2/5

**Premises.** nested Prop R2's (M1)–(M4), together with the typical-factorisation model of nested §5: the algebraic value
of a node, of log-size X, is the sum of the log-sizes of its children plus a B-smooth part of log-size ≤ ΛU. Here
Λ = log B, U is the largest u with u log u ≤ Λ, and a node costs L[α] only if each smoothness event has u ≤ U. Nodes may
use their own fields (any d, H_f, H_g); algebraic children are closed in the parent's field. An **algebraic
component** is the set of nodes linked by algebraic edges, all in one field. *Budget:* every component has at most
κ(d − 2) points, for a constant κ ≥ 1. The known moment methods need ≤ d + 1 (nested Lemma OD).

**Prop R3 [model computation].** Under these premises, α ≥ 2/5.

*Proof.* Claim: every node v satisfies X_v := log H_f + (d/2)z_v ≤ (1 + 2κ)ΛU, where z_v is the target's log-size
(for an algebraic target, X_v − z_v is the cofactor). The closing node w of a child of log-size y_w has value
≥ log H_f + (d/2 − 1)y_w. So y_w ≤ 2X_w/(d − 2). Write X_v ≤ ΛU + Σ_{w child of v} y_w ≤ ΛU + (2/(d − 2))Σ_w X_w.
Unrolling over the component, rooted at v, gives X_v ≤ ΛU·Σ_u (2/(d − 2))^{depth(u)} ≤ ΛU(1 + 2κ(d − 2)/(d − 2)) =
(1 + 2κ)ΛU, since every node below v has weight ≤ 2/(d − 2) and there are ≤ κ(d − 2) of them. This is R2's
"algebraic feasibility gives log H_f + (d/2)z = O(ΛU)", with the constant 1 + 2κ in place of 1, and R2's window
argument applies verbatim with A replaced by (1 + 2κ)²A:
- if log H_g ≤ t/4, then log H_f ≥ t/(2(d + 1)), so √(tz)/2 ≤ (1 + 2κ)ΛU and z ≲ A;
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
(r + 1 − C)δ ≤ 3α − 1, so α ≥ (r + 2 − C)/(2r + 5 − 2C) when r + 1 > C.

**E2 (exact-rational simplex; all pre-registered predictions hold).**

| C (component ≤ d^C points) | α*(C) | levels r | closed form |
|---|---|---|---|
| 1 | **2/5** | 0 | 2/5 (= R3) |
| 5/4 | 7/18 ≈ 0.389 | 1 | (3 − C)/(7 − 2C) |
| 3/2, 2, 5/2 | **3/8** | 1 | 3/8 |
| 8/3, 3 | 4/11 ≈ 0.364 | 2 | 4/11 |
| 4 | 5/14 ≈ 0.357 | 3 | 5/14 |
| 6 | 7/20 = 0.35 | 5 | 7/20 |
| 10 | 11/32 ≈ 0.344 | 9 | 11/32 |
| ∞ | → 1/3 | → ∞ | (r + 2)/(3r + 5) |

The LP agrees with the closed form to 10⁻⁹ at every grid point. Each closed form is feasible exactly and infeasible
10⁻⁶ below. At α = 3/8 (δ = 1/4, ζ_1 = 1/2) the first-level value has exponent v₀ = 7/8. Its children number
D^{3/8} = d^{3/2}: the "3/8 vs 1/4" count of codex §11 and r20 #1.

*Reading.* Any fixed C > 1 gives an exponent below 2/5 in the model, and no finite C reaches 1/3. The single missing
input is:

**Hypothesis LS_C (joint smooth supply, informal) [open].** For a degree-d field law and configurations of s ≤ d^C
points arising as complete descent certificates, the number of field coefficient vectors c ∈ [−H, H]^d at which all s
affine forms F_c(P_i) have the prescribed large divisors and B-smooth cofactors is ≥ L^{−O(1)} times its heuristic
value, on average over configurations.

For s > d + 1 the forms are linearly dependent (nested Lemma OD). This is a "linear equations in smooth numbers"
problem of growing complexity (s ≈ d^C forms in d variables), with coefficients of size up to exp(Λ_{7/8}), at
smoothness u ≍ D^α. The nearest known results are Harper's a + b = c (one equation, three variables, coefficients 1;
y ≥ (log x)^K) and Matthiesen–Wang (finite complexity, Ψ(x, y) ≥ x/(log x)^8). Neither reaches it (§11).

### 7.3 r20 #1, evaluated adversarially

r20 #1 proposes first and second moments of Z_{A,q} = #complete two-level certificates for q, at α = 3/8, d ≍ D^{1/4},
with saddle-point-tilted weights and sparse transference, and requires polylogarithmic precision.

- **The exponent is right [model; E2].** Two levels at 3/8 need C = 3/2. They are the cheapest step below 2/5.
- **The precision requirement can be relaxed to L^{O(1)} [proved conditional, §8].** With a pool of J_f = L^{O(1)}
  fields (lv5unif Prop MF's mechanism), one needs, for a (log N)^{−C₀} fraction of q,
  E_A Z ≥ L^{−c} and E_A Z² ≤ L^{c}(E_A Z)². By Paley–Zygmund, Pr_A[Z ≥ 1] ≥ L^{−2c}. At bounded depth the
  exhaustive search finds a certificate whenever Z ≥ 1, at cost L^{O(1)}. The split then sees a target set of
  polylog density, and the γ^{2k} amplification hits only the pool's coverage, not the moment ratio.
- **Consequence for r20's kill test.** r20 would kill the transference route if its majorant lost exp(Ω(d)) in the
  relative moment ratio. After the pool trick, a loss of exp(O(s)) with s = d^{3/2} = D^{3/8} is exp(o(Λ_{3/8})), and
  even exp(O(s log log N)) is L[3/8, O(1)]. Both are affordable. So local-factor losses of e^{O(1)} per form, which a
  transference majorant typically incurs, no longer kill the route.
- **What remains [open].** A lower bound of the right order L^{−O(1)} times the main term for the number of c at which
  s ≈ d^{3/2} dependent affine forms are simultaneously smooth with prescribed divisibility. The main term itself is
  tiny: ∏_i ρ_i ≈ exp(−s·u log u) = L[3/4]^{−Θ(1)} per configuration. So any additive error must be smaller than that,
  which no existing circle-method or transference estimate gives at L-scale.
- **Verdict.** It survives as the right analytic target and is the cheapest one. We agree with r20's 3% and do not
  raise it: the relaxation removes the precision objection, not the core estimate.

---

## 8. Theorem NR-q: the weakest hypothesis we know that suffices

**Setting.** Fix α ∈ (1/3, 1/2), parameters as in nested §8 or §7.2, a field law 𝔉 (admissible fields of degree
d ≍ D^δ, sampled in time e^{o(Λ_α)}), and a deterministic in-field descent procedure 𝔇. Examples: nested's DESC (first
good point, recursing), or exhaustive search at bounded depth. On every admissible field and every target, 𝔇 must make
at most e^{o(Λ_α)} tree nodes and scan boxes of ≤ L[α, 2s] points (nested Lemma TS). It calls the recursive factoring
oracle on integers ≤ N^{τ_F}, τ_F < 1, and outputs either FAIL or an RW tree (nested Lemma RW) whose leaves lie in the
factor base. For a rational prime q put p_q := Pr_{A∼𝔉}[A admissible and 𝔇_A(q) ≠ FAIL].

**Hypothesis T_q(α; 𝔉, 𝔇, c, C₀).** With J = primes in (Y, 2Y], Y = N^{1/(2k)}, k = ⌈κD^α⌉:
#{q ∈ J : p_q ≥ e^{−cΛ_α}} ≥ (log N)^{−C₀}·#J.

**Theorem NR-q [proved conditional on T_q, nested's Hypothesis LV6′ for each field, character columns via Lemma
X_ERH (ERH) or lv5unif's factorial-bit repair, and the multi-field post-processing of lv5unif §11 (proved there in
outline; codex escapeR2_1 (2.5) gives the square-root identity)].** There is a Las Vegas algorithm that completely
factors every N in expected time L_N[α, c′] for a finite c′ = c′(α, c, C₀, κ, …).

*Proof.*
1. **Pool.** Sample J_f := ⌈4e^{cΛ_α}⌉ fields A_1, …, A_{J_f} ∼ 𝔉 independently (factor exits end the run). Let
   G_pool := {q ∈ J : ∃j, 𝔇_{A_j}(q) ≠ FAIL}. For a good q (p_q ≥ e^{−cΛ_α}), Pr[q ∉ G_pool] ≤ (1 − p_q)^{J_f} ≤ e^{−4}. By
   Markov, with probability ≥ 1 − 2e^{−4} ≥ 1/2 at least half the good q are in G_pool, so |G_pool| ≥ γ#J with
   γ := (log N)^{−C₀}/2. Call such a pool *good*.
2. **Split.** G_pool is a fixed set before any Dixon draw, so Lemma S″ (recurse §3.2) applies to it with the window
   bookkeeping of lv5unif Theorem K-U (|G_pool| ≥ 4k³2^{ω} by lv5unif (P6)). An attempt is accepted with probability
   ≥ (γ/(4 log N))^{2k}/(32k²) = L[α, −2κ(1 + C₀) − o(1)].
3. **Descent.** For each of the 2k primes of an accepted split, run 𝔇_{A_j} for j = 1, 2, … until success: at most
   J_f·e^{o(Λ)}·L[α, 2s] work, plus oracle calls on integers ≤ N^{τ_F}. By nested Lemma RW each tree gives
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
7. **Restarts.** A bad pool is detected only by running time. Run each pool with a budget of 4× the good-pool expected
   cost, and resample the pool on overrun. The expected number of pools is ≤ 2.
8. **Recursion.** F(N) ≤ L[α, p₀] + L[α, δ₀]F(√(2N)) + L[α, δ₁]F(N^{τ_F + o(1)}) with finite exponents, so Lemma F_α
   (lv5unif §10; for the third term, nested's strong induction with Λ_α(N^τ) = τ^α Λ_α(N)) gives L[α, c′]. ∎

*Relations to known hypotheses.*
- AD (nested §8) ⇒ T_q with p_q ≥ 1/2 for every q, using DESC. So NR-q contains NR.
- T_q is field-averaged and target-pointwise. That is the form Prop SI(ii) allows: the root target q is exogenous to
  the pool. The endogeneity now lives entirely inside one field's tree, i.e. in Hypothesis LS_C.
- A sufficient condition for T_q is r20 #1's moment pair at L^{O(1)} precision (§7.3), when 𝔇 is exhaustive at
  bounded depth.
- *What T_q is not:* it is not proved for any α < 2/5, unconditionally or under ERH. At α = 2/5 with one level it
  follows from lv5unif's Theorem D-U, in the stronger, density form.

---

## 9. The missing lemmas, collected

1. **Hyp LS_C / T_q** (in-field joint smooth supply; §§7.2, 8). This is *the* missing lemma: every other route
   reduces to it or to something stronger. At C = 3/2 it gives 3/8; as C → ∞ it gives → 1/3.
2. **Contracting image certificate** (codex (5.2)): for θ from principalisation, ψ(θ) ≡ w²·a·∏q_ν (mod N) with
   q_ν ≤ q^λ, λ < 1. It is not known how to produce this without a re-split (§4).
3. **Coset re-split** (Prop RS): class equidistribution of products of descendable primes in every coset, plus
   contraction. False for adversarial sets; open for the actual ones; no contraction in any case.
4. **A fresh-and-tied construction** (§2.3, §5): fields that share a Galois structure rich enough for parity transfer
   (TI with k_1 ≠ Q) while carrying independent randomness per subfield. No candidate is known. KT excludes the
   Kummer/odd-isogeny version for S_d-fields.
5. **Small-discriminant fields with explicit ψ** (§4): would make L_{|Δ|}(1/2) = L_N[1/3]; equivalent, for small fixed
   fields, to root-finding modulo composites.

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

E1–E3 check algebra and arithmetic only. They say nothing about relation supply.

---

## 11. Literature

(To be completed by the Opus literature critic; preliminary.)
- **Norm relations / Brauer relations.** Cross-subfield transfer (Example TB) is the parity shadow of norm relations in
  Q[G]. Those are used to compute units, class groups and S-units of a Galois field from its subfields: Biasse, Fieker,
  Hofmann and Page, "Norm relations and computational problems in number fields", J. LMS (2022). For multiquadratic
  fields: Bauch, Bernstein, de Valence, Lange and van Vredendaal, Eurocrypt 2017; Biasse and van Vredendaal, ANTS 2018.
  Kani and Rosen relate idempotent relations to class groups and Jacobians. Lemma TI is elementary; we do not claim it
  is new.
- **Class groups.** dBPW (arXiv 2512.01588, ERH, L_{|Δ|}(1/2), rigorous). Biasse (arXiv 0912.1927, 1204.1292) and
  Biasse–Fieker (LMS JCM 2014): L(1/3) for families, under smoothness and relation-generation heuristics. Bach (Math.
  Comp. 1990): ERH generators of norm ≤ 12 log²|Δ|.
- **Smooth numbers in linear systems.** Harper (Compositio; arXiv 1408.1662), a + b = c. Matthiesen–Wang (arXiv
  2211.16892). Lachand. None covers growing complexity at L-scale (nested §5).
- **Hyperelliptic 2-descent.** Schaefer (JNT 1995), Stoll (Acta Arith. 2001).
- **Permutation groups.** Dixon–Mortimer, *Permutation Groups*, GTM 163, §5.2.

---

## 12. Reviews

(Pending: codex escapeR2_2 adversarial referee; Opus proof critic; Opus literature critic.)

---

## 13. Candidate card

- **Mechanism.** Escape nested Prop R2 by exporting algebraic children (norms, fresh shared-root fields, generators,
  Galois/Kummer transport), by class-group principalisation, or by changing the problem (re-splits, recursion,
  independent descents).
- **Model.** Theorem U / nested architecture; certificate class of §1; size model of recurse §3 and nested R2.
- **Result.** No rigorous or ERH-conditional α < 2/5.
  - Proved: Lemma TI and Cor EX (export impossible across fields with k_1 = Q); Example TB (tied fields transfer);
    Prop KT (Kummer/odd-isogeny transport returns to K for S_d, d ≥ 7, given (F)); Lemma PR and Prop RS(i, ii).
  - Model: Prop R3 (O(d)-point components ⇒ 2/5); α*(C) (E2).
  - Conditional: Theorem NR-q (T_q ⇒ L[α]).
- **Failure scope.** Certificates built from field elements squared in a common Galois field; designs under (M1)–(M4).
  Outside: non-square certificates, structured polynomial selection, fields with explicit ψ and small discriminant.
- **Next step.** The smallest case of Hyp LS_{3/2}: two-level certificates at α = 3/8 with L^{O(1)} precision. Find
  any lower bound for smooth points of an overdetermined affine system of d^{3/2} forms in d variables at L-scale.
