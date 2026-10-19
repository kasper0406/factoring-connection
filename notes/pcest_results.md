# pcest: the second moment PC for bounded-k Dixon descent

Agent `pcest`, phase 8, 2026-10-19.
- Pre-registration: `notes/pcest_prereg.md` (commit `8aa8898`, before E1 ran).
- Code: `src/pcest/` (`geometry.py`, `params.py`, `toy.py`, `e1.py`). Tests: `tests/test_pcest.py`.
- Raw data: `experiments/raw/pcest_e1.jsonl`, summary `experiments/raw/pcest_e1_summary.json`.
- GPT-6.1-sol: blind round `codex/pcest_1.{md,out}`; adversarial round `codex/pcest_2.{md,out}` (§7).
- Opus critics: one on proofs, one on literature and novelty (§7).

## Verdict

**We did not obtain a complete theorem. Lenstra–Pomerance's rigorous L[1/2, 1+o(1)] still stands.**
Of the three missing inputs that hypD left for the fixed-k design, two are now proved, after a change
of design. One precisely stated estimate remains open.

In plain terms:
- The *rational* half of the pair correlation (PC_rat) is the problem hypD considered hardest. It is a
  binary correlation of smooth values along a line, of "smooth twins n, n + h" type.
- It becomes easy once the rational polynomial is randomized in two dimensions, g = g₁x − g₀ with
  (g₀, g₁) random. Pairs of rational values then fill a two-dimensional lattice region instead of a
  line.
- After that change, an elementary divisor-switching argument bounds the pair count within a
  polylogarithmic factor. No sieve or equidistribution theorem is needed.
- GPT-6.1-sol found the same route independently in its blind round.
- The *algebraic* half (PC_alg) is still open. It needs smooth numbers in arithmetic progressions to
  very large moduli to be equidistributed among the lifts of a moving class, in an L² sense, under
  rational-smoothness weights. We found no published theorem that gives this, and we explain below why the known
  tools fall short.

Labelled results:
1. **Design (a″_k) and Lemma GB [proved; elementary].**
   - Choose (g₀, g₁) uniform among coprime pairs in [M/2, M]², with M = ⌊N^{1/(d+1)}⌋.
   - The generalized base-(g₀/g₁) expansion gives f̂ with F̂(g₀, g₁) = N and ‖f̂‖ = O_d(M).
   - Then f = f̂ + (g₁x − g₀)R.
   - There is no Kleinjung residue constraint on g₀ mod g₁: the leading coefficient is allowed to
     vary instead.
   - All norm sizes, and hence θ and c, are unchanged.
2. **Theorem FM [proved, regime (B), design (a″_k)].** This is the first moment, with hypD §4.1′
   items 1–5 written out line by line (§2).
   - Lemma 8.5′ is re-proved with every cited range checked.
   - The rational side now needs no Bombieri–Vinogradov input at all. This also removes hypD's
     constraint k ≥ (3/4)(d+1).
3. **Theorem R (PC_rat) [proved, design (a″_k)].** E[N_r(N_r − 1)] ≤ (log N)^{O(1)}·(E N_r)² for the
   rational-only count (§3). It rests on two elementary lemmas:
   - Lemma Q: the pair image lattice has index q = |Δ|/ℓ ≤ 2S² and graph form.
   - Lemma DS: an upper bound for smooth numbers in a short interval and a progression to a small
     modulus, uniform in the residue.
4. **Lemma Kern [proved].** If the R-height satisfies H ≥ V² (V = √ℓ·S, the box size), then the
   coefficient law of the pair (F, F′) is dominated pointwise by C_d × (uniform measure on
   Γ ∩ Φ(C_H)). Against the bounding box [−T, T]² the factor becomes C_d·S²/q (Prop 4.2) (§4.1).
   - This is the "2D kernel lemma" with codex's H ≫ V² condition.
   - Its cost: the h-floor rises from 1/(4k) to 1/(2k). That binds only for k < d(d+1)/4, so **the
     infimum constant is unchanged**.
   - Example: d = 9, k = 20 gives θ = 0.325 (c = 0.937) instead of 0.3125 (c = 0.928). From k = 23
     onward, θ is the same as before.
5. **What remains: Hypothesis LD [open] (§4.2).** This is a restricted, rational-smoothness-weighted
   Barban–Davenport–Halberstam-type upper bound, at precision L^{o(1)}, for smooth numbers in the
   lift classes F̂ + Gj mod G|Δ|.
   - It is equivalent in substance to hypD's RBDH′ and to codex's statement (6), but stated as
     exactly what the proof consumes.
   - The *conclusions* of the max-residue (L¹) Bombieri–Vinogradov theorems cannot imply it (Prop B5).
     Naive uses of Harper's all-residue BDH and of the large sieve lose polynomial factors (§4.3).
   - These are barriers to naive applications, not to the method. §4.3 gives a plausible route: an L²
     version of LV's character-majorant argument, plus a large-values bound for structured residues.
6. **Theorem A [conditional].** Assume Hypothesis LD and the LV §6 adaptation (square roots and
   characters at our heights, inherited from lemmaM Corollary DD). Then every N is factored in
   expected time L_N[1/2, (1 + θ(d, k))/√2 + o(1)].
   - Infimum over k: c(d) = (1 + 3/(d+1))/√2.
   - d = 7: 0.9723. d = 9: 0.9192. d = 11: 0.8839. d = 19: 0.8132.
   - c < 1 holds exactly for odd d ≥ 7 (d = 5 gives 1.0607). At d = 7 it needs k ≥ 13; at d = 9,
     k ≥ 12 (with the old 1/(4k) floor: 11 and 11).
   - HECM's rigour (LPP Part III) is the same caveat as in LV and hypD (§5).

Compared with hypD, the conditions for (a′_k) were FM (outline) + PC_rat (open) + PC_alg (open) +
LV §6. **They are now LD + LV §6.**

**Experiments (pre-registered E1, toy, exact).** Setup: 8 instances at 56–64 bits, d = 3, and 9 600
(field, ℓ) cells per design.
- PC_rat(2D) = 1.007–1.013 (P1 holds).
- PC(2D) = 1.06–1.28 (P2 holds).
- |PC_rat(2D) − PC_rat(1D)| ≤ 0.003 (P3 holds).
- E N(2D)/E N(1D) = 0.77–0.94 (P4 holds). The 2D rational values are slightly larger, a constant
  factor.
- P5: 0 identity violations and 160/160 Lemma Q checks. The harness also flagged 3 degenerate points
  with G = 0 in irregular 1D lattices (m ≡ 0 mod ℓ). These are not identity failures (disclosed, §6).

The toy cannot test any asymptotic claim; it checks the code and shows that the 2D design behaves like
the 1D one in practice.

### Per-design table (L[1/2, c]; "inputs" = what a proof still needs)

| design | c | status | missing inputs |
|---|---|---|---|
| Lenstra–Pomerance 1992 | 1 | proved (prior art) | — |
| Vallée 1991 (θ = 2/3) | √(4/3) ≈ 1.155 | proved (prior art) | — |
| hypD (a′_k): rational side x − m, LP shape test | (1+θ)/√2, θ → 3/(d+1) (d = 9: 0.920) | conditional | PC_rat (binary, open), PC_alg (open), LV §6; FM by the §2 method with Harper BV on the rational side (k ≥ 3(d+1)/4): outline only, not written out |
| **(a″_k): 2D-randomized rational side, LP shape test** | **(1+θ)/√2**, θ(d,k) = 2s + (d−1)/(4k) + max(1/(2k), s − (d−2)/(4k)); infimum (1+3/(d+1))/√2: d = 7: 0.972; d = 9: 0.919; d = 11: 0.884; d = 19: 0.813 | **conditional on LD + LV §6** | **Hypothesis LD** (restricted lift-dispersion upper bound, §4.2); LV §6 adaptation; FM, PC_rat and the kernel transfer are proved |
| (a″_k), recursive shape test | √((1+√2)θ) (d = 9: 0.851) | conditional on LD at every scale | as above, at all scales |
| (a′_k)/(a″_k) at Horn 1 (k ≍ (log N)^{1/4}, HECM shape test) | √(6/(d+1)) infimum (d = 9: 0.775) | conditional | LD at precision exp(−o(Λ/k)) (B3); explicit-error LV §5. Theorem R's polylog loss would suffice there, but FM and LD would need explicit errors |
| two randomized algebraic sides (codex) | (1 + 4/(d+1))/√2 (limit) | conditional | analogue of LD with special-root conditioning; completely-split target law |
| Design V (hypD §5) | √(2/3) ≈ 0.816 | conditional on RK | unchanged by this note |

---

## 1. Design (a″_k): a two-dimensional random rational side

Everything not mentioned here is as in hypD §1–§3: Dixon descent, Theorem DD, Theorems S/S′ and K, the
field-selection rule and the LP shape test. Those results treat the field A abstractly and act only on
R and A, so they apply verbatim to the field law below.

**Notation.**
- L = L_N[1/2, 1] and Λ = (log N log log N)^{1/2}.
- s = 1/(d+1), M = ⌊N^s⌋, B = L[1/2, a], ϱ(x, y) = Ψ(x, y)/x, and 𝒮_B = the B-smooth integers (sign
  ignored).
- ℓ is uniform among the primes of J = (2^{j−1}, 2^j], with 2^{kj} ≍ √N, so ℓ = N^{1/(2k)+o(1)}.
- S = L^σ with σ > 0 fixed (S² will be T_SQ up to L^{o(1)}), and V = ⌊√ℓ·S⌋.
- H = N^h.

**Preprocessing [standard].** Trial division up to B, at cost B·(log N)^{O(1)} = L^{a+o(1)}. This is
below the final cost L^{c}, since c ≥ 2a. Afterwards N has no prime factor ≤ B, and is not a prime
power. This replaces Lemma Ω for everything below; 2^ω = L^{o(1)} still holds.

**Field law 𝒜.**
- (g₀, g₁) is uniform on 𝒢 := {(g₀, g₁) ∈ Z² ∩ [M/2, M]² : gcd(g₀, g₁) = 1}, and G(x, y) := g₁x − g₀y.
- c is uniform on C_H := [−H, H]^d ∩ Z^d, and R_c(x, y) := Σ_{i<d} c_i x^{d−1−i}y^i.
- F_A := F̂_g + G·R_c, with F̂_g from Lemma GB. For P = (a, b):
    F_A(P) = F̂_g(P) + G(P)·⟨c, v(P)⟩,  v(P) := (a^{d−1−i}b^i)_{0≤i<d}.
- If gcd(g₁, N) > 1 or gcd(G(P), N) > 1, a factor of N has been found. Otherwise put
  ψ: Z[x] → Z/N, x ↦ g₀g₁^{−1}. It is well defined on Z[x]/(f): g₁^d f(g₀/g₁) = F(g₀, g₁) = N ≡ 0
  (mod N).
- F_A(P) = 0 is excluded from the counts.

**Lemma GB (generalized base expansion) [proved].** Let g₁ ≥ 2 and gcd(g₀, g₁) = 1. Define
T_d := N, and for i = d, d−1, …, 1:
  a_i ≡ T_i g₀^{−i} (mod g₁), |a_i| ≤ g₁/2,  T_{i−1} := (T_i − a_i g₀^i)/g₁,
and finally a_0 := T_0. Then F̂(g₀, g₁) := Σ_i a_i g₀^i g₁^{d−i} = N. If g₀, g₁ ∈ [M/2, M], then
|a_0| ≤ 2^{d+1}M + d2^{d−1}M, so ‖f̂‖_∞ = O_d(M).
*Proof.*
1. T_i ≡ a_i g₀^i (mod g₁) by the choice of a_i, so T_{i−1} ∈ Z.
2. By induction, T_i = Σ_{j≤i} a_j g₀^j g₁^{i−j}. At i = d this is the claim.
3. For the bound, a_0 = (N − Σ_{i≥1} a_i g₀^i g₁^{d−i})/g₁^d. Here N/g₁^d ≤ 2^d N/M^d ≤ 2^{d+1}M, and
   each |a_i|(g₀/g₁)^i ≤ (g₁/2)·2^i ≤ 2^{i−1}M. ∎
(`src/pcest/geometry.py: gbe`; tested.)

*Remarks.*
- **Prior art.** The non-monic linear side g₁x − g₀ with a base-(g₀, g₁) expansion is Buhler–Lenstra–
  Pomerance's base-(m₁, m₂) variant (LNM 1554 (1993) 50–94; Murphy, thesis, ANU 1999, §2.3.2). Lemma GB
  is that construction with explicit bounds. The design is not new. What we did not find elsewhere is
  its use to randomize the rational side in a rigorous pair-correlation argument (Theorem R).
- **Kleinjung constraint.** Kleinjung (Math. Comp. 75, 2006) prescribes a small leading coefficient
  a_d, and therefore needs a_d g₀^d ≡ N (mod g₁). Here a_d is *determined* by (g₀, g₁), so (g₀, g₁)
  can be sampled freely subject only to coprimality. Codex makes the same point.
- **Sizes.** The price is a_d ≈ g₁ rather than a small leading coefficient. That is irrelevant at the
  exponent level, because LV's f̂ already has coefficients ≈ M.

**Lemma 1.2 [proved].** For every P = (a, b), g₁^d F_A(P) ≡ N b^d (mod G(P)).
*Proof.* g₁·P − b·(g₀, g₁) = (G(P), 0). Since F_A is an integral form, g₁^d F_A(P) = F_A(g₁P) ≡
F_A(b·(g₀, g₁)) = b^d F_A(g₀, g₁) = b^d N, because G(g₀, g₁) = 0. ∎

**Relations.** If G(P) = ±ℓn with n ∈ 𝒮_B and F_A(P) ∈ 𝒮_B, then
  ψ(a − bα) = G(P)/g₁ = ±ℓn/g₁ (mod N).
This is hypD's descent relation, multiplied by the per-field unit g₁^{−1}.
- One extra column, the parity of the number of relations, makes g₁ a square in every dependency.
  This is the same device LV and Buhler–Lenstra–Pomerance use for f_d.
- **Irreducibility.** LV Lemma 5.1 transfers. Its Turán-sieve count uses only the affine structure
  f ∈ f̂ + (x − m)·Z[x]_{d−1}, and the condition "g(m) ≡ n (mod r)" for the reductions g mod r. With
  g₁x − g₀ in place of x − m, the condition becomes ĝ(g₀, g₁) ≡ N (mod r) for r ∤ g₁N. Claim 5.2
  counts monic irreducibles with a prescribed value at an affine point m, here m = g₀g₁^{−1} mod r
  (r ∤ g₁). It is unchanged.
  [transfer of a published proof, checked at the level of its statements]
  - a_d ≠ 0 in Lemma GB, because a_d ≡ N g₀^{−d} (mod g₁) and gcd(g₁, N) = 1.
- A reducible f either gives F₁(g₀, g₁)F₂(g₀, g₁) = N with both factors ≠ ±1, which is a factorization,
  or is rare.
- **Quantified (Opus proof critic).** LV's bound transplants as Pr(f reducible) ≤ H^{−1/3+o(1)}. Field
  selection keeps a field with law ∝ q_A, which amplifies rare events by at most 1/μ. So
  Pr(kept field reducible) ≤ H^{−1/3+o(1)}/μ = N^{−h/3+o(1)}. The degree-drop event
  a_d + g₁c₀ = 0 has probability ≤ 1/(2H + 1). Reducible fields are rejected.
  - Status: transfer of LV Lemma 5.1, checked at the level of statements. It is listed among the
    assumptions of Theorem A.

**Counted candidates and the capped search.**
- Put 𝒫_ℓ := {(a, b) : b prime, V/2 < b ≤ V, |a| ≤ V, gcd(a, b) = 1}, and
  Λ_ℓ(g) := {P : ℓ | G(P)}, which has index ℓ when ℓ ∤ g₁.
- The search enumerates Λ_ℓ(g) ∩ [−V, V]² from a Gauss-reduced basis, tests membership in 𝒫_ℓ, and
  tests G(P)/ℓ and F_A(P) for B-smoothness.
  - Primality of b uses a deterministic polynomial-time test (AKS).
  - Smoothness is tested by HECM, as in LV Cor. 3.19.
- N_ℓ(A) := #{P ∈ 𝒫_ℓ ∩ Λ_ℓ(g) : |G(P)| ≥ MV/32, G(P)/ℓ ∈ 𝒮_B, F_A(P) ∈ 𝒮_B}. Each such P is a
  certified descent relation. ℓ² ∤ G(P) holds automatically, because ℓ > B.
  - The **filter |G(P)| ≥ MV/32** was added after codex pcest_2 (§7.2). Without it, candidates with
    |G| as small as ℓ have algebraic values smaller by a polynomial factor. Their smoothness density
    can then exceed ρ_F by L^{Θ(1)}, which breaks Theorem 4.3's comparison. FM only ever uses
    |G| ∈ [MV/32, 2MV], so the filter costs nothing there.
  - A capped randomized smoothness test is repeated O(log N) times, so that a smooth value is
    recognized with probability ≥ 1 − N^{−2}. This costs a polylogarithmic factor.
- **Lemma 1.3 (regularity) [proved].** Suppose λ₁(Λ_ℓ(g)) < √ℓ/(2S), attained at w. Then
  Λ_ℓ(g) ∩ 𝒫_ℓ = ∅. Otherwise #(Λ_ℓ(g) ∩ [−V, V]²) ≤ C_abs·S² for an absolute constant C_abs (e.g. 100).
  *Proof.*
  - w is primitive in Z²: its content e satisfies e ≤ |w| < ℓ, so gcd(e, ℓ) = 1 and w/e ∈ Λ_ℓ(g); by
    minimality e = 1.
  - So w is primitive in Λ_ℓ(g) (index ℓ), and the lattice is the union of the lines jw′ + Zw, j ∈ Z,
    which are parallel to w at mutual distance ℓ/|w| > 2V.
  - Every point of [−V, V]² has norm ≤ √2V < 2V, so it lies on the central line Zw. The only primitive
    points there are ±w, and |w| < V/2 < b for every point of 𝒫_ℓ.
  - In the regular case take a reduced basis (λ₁ ≤ λ₂, λ₁λ₂ ≤ (2/√3)ℓ). The translates of the
    fundamental cell F at the lattice points of K = [−V, V]² are disjoint and lie in K + F. So the count
    is ≤ (area(K) + perimeter(K)·diam F + π diam(F)²)/ℓ ≤ 4S² + O(V/λ₁) + O(λ₂/λ₁) ≤ C_abs S², using
    λ₁ ≥ √ℓ/(2S), V/λ₁ ≤ 2S² and λ₂/λ₁ ≤ (2/√3)ℓ/λ₁² ≤ 8S²/√3. ∎
- So the cap T_SQ := C_abs S²·(log N)^{O(1)} is never binding on a lattice that contains a candidate.
  Irregular lattices need no separate exclusion: they contribute nothing.

**Sizes [proved arithmetic].** For P ∈ 𝒫_ℓ and g ∈ 𝒢:
- |G(P)| ≤ 2MV, so |G/ℓ| ≤ X_r := 2MV/ℓ = N^{s−1/(4k)}L^{O(1)};
- |F̂_g(P)| = O_d(MV^d);
- |G(P)⟨c, v(P)⟩| ≤ 2dMHV^d, so |F_A(P)| ≤ X_F := C_d MHV^d = N^{s+h+d/(4k)}L^{O(1)}.

With ρ_r := ϱ(X_r, B) and ρ_F := ϱ(X_F, B), these are exactly hypD's exponents, so
θ = 2s + h + (d−1)/(4k) as in `src/pcest/params.py`.

## 2. FM in full (regime (B), design (a″_k))

### 2.1 Cited results, with their exact ranges

Throughout, x = N^{ξ+o(1)} with ξ = s + h + d/(4k) fixed, y = B = L[1/2, a], and u = log x/log y ≍
(log N/log log N)^{1/2}. So (log x)^K < y < x^{1/log log x} for any fixed K and large N.

- **(H) Harper, arXiv 1208.5992; LV Facts 8.8–8.12.** The constants are absolute and effective.
  - Fact 8.11 (Harper Prop. 2): for 0 < η < 1/80, y ≤ x^{9/10}, x^η ≤ Q ≤ √x and M₀ = x^η,
      Σ_{M₀≤r≤Q} Σ_{χ*_r} Σ_{t∈[M₀,Q]} φ(t)^{−1} Σ_{χ_t induced by χ*_r} |Ψ(x, y; χ_t)|
        ≪ log^{7/2}x · √Ψ(x, y)·(Q + x^{1/2−η} log²x).
- **(LV8.5) LV Lemma 8.5 and its proof.** Fix ε > 0. For (log x)^K < y < x^{1/log log x},
  x^ε ≤ Q ≤ √Ψ(x, y), ω₀ → ∞ and ω₀ = y^{O(1)}, the majorant
    W(Q) := Σ_{t∈[Q/ω₀, Q], t y-smooth} φ(t)^{−1} Σ_{χ_t ≠ χ₀} |Ψ(x, y; χ_t)|
  (LV eq. (8.2)) satisfies W(Q) ≪ Ψ(x, y)ϱ(Q, y)(H(u)^{−c} + y^{−c}) + Q√Ψ(x, y) log^{7/2}x, where
  H(u) = e^{u/log²(u+1)}.
  - LV's proof bounds W itself, not only the max-residue sum, which W majorizes.
  - The constant is effective; Siegel and Page characters are treated in LV's G₁ part.
- **(T) Tenenbaum (LV Fact 3.16, Cor. 3.17).** If ω(t) ≤ y^c/log(1 + u), then
    Ψ_t(x, y) = (φ(t)/t)Ψ(x, y)(1 + O(log(1+u) log(1+ω(t))/log y))(1 + O(ω(t)/y)).
  In our range the error is o(1) uniformly for t ≤ x, since ω(t) ≤ log x.
- **(Hi) Hildebrand.**
  - LV Fact 3.13: for log y ≥ (log log x)^{5/3+ε} and 1 ≤ z ≤ y^{5/12},
      Ψ(x(1 + 1/z), y) − Ψ(x, y) = (Ψ(x, y)/z)(1 + O(log(u+1)/log y)).
  - Ψ(2x, y) ≤ 2Ψ(x, y) (LV ref. [19, Thm 4]). Hence ϱ(z₂) ≤ 2ϱ(z₁) whenever z₂ ≥ z₁.
- **(C8.16) LV Claim 8.16.** If u → ∞ and c ≥ 0 is fixed, then ϱ(x, y) = ϱ(xy^c, y)(log x)^{O(⌈c⌉)}.
  So ϱ changes by (log N)^{O(1)} = L^{o(1)} when its argument moves by a factor y^{O(1)} = L^{O(1)}.
- **(LV5) LV Lemmas 5.13, 5.18, 5.20.** These concern only the linear map φ_P(c) = ⟨c, v(P)⟩. They do
  not involve m, so they transplant verbatim (§2.3, step 2).

### 2.2 Lemma 8.5′ (restricted BV for moduli ℓ·(smooth)) [proved; regime (B)]

**Statement.**
- Let 0 < q < ξ/2 with q > 1/(2k) + η₁, for some fixed η₁ > 0, and Q = N^{q+o(1)}.
- Let Y = 2^{j−1} and 𝒯_Q := {t = ℓs ∈ (Q/2, Q] : ℓ ∈ (Y, 2Y] prime, P⁺(s) ≤ B}.
- Let δ_𝒯 := ϱ(Q/Y, B)/log Y. Then #𝒯_Q ≍ Qδ_𝒯.
- Conclusion: Σ_{t∈𝒯_Q} max_{(r,t)=1} |Ψ(x, B; t, r) − Ψ_t(x, B)/φ(t)| ≪ Ψ(x, B)δ_𝒯·ε_N, where
  ε_N := H(u)^{−c} + B^{−c} + x^{−η₂} for some fixed η₂ > 0.

*Proof.*
- **Setup.**
  - #𝒯_Q = Σ_ℓ (Ψ(Q/ℓ) − Ψ(Q/2ℓ)).
  - Fact 3.15 gives Ψ(z/2) = 2^{−α}Ψ(z)(1 + o(1)) with α ≥ 1/2, so each difference is ≍ Ψ(Q/ℓ).
  - Summing gives #𝒯_Q ≍ Qδ_𝒯.
  - Since ℓ > B ≥ P⁺(s), ℓ ∤ s and φ(t) = (ℓ − 1)φ(s). Every character mod t is χ_ℓχ_s, and t
    determines (ℓ, s).
  - The left side is at most W′ := Σ_{t∈𝒯_Q} φ(t)^{−1} Σ_{χ≠χ₀ mod t} |Ψ(x, B; χ)|.
- **Part 1: χ_ℓ principal, χ_s ≠ χ₀.**
  - Every B-smooth n is coprime to ℓ, so Ψ(x, B; χ₀^{(ℓ)}χ_s) = Ψ(x, B; χ_s).
  - So this part equals Σ_ℓ (ℓ−1)^{−1} Σ_{s∈(Q/2ℓ, Q/ℓ], smooth} φ(s)^{−1} Σ_{χ_s≠χ₀} |Ψ(x, B; χ_s)|.
  - By positivity this is ≤ Σ_ℓ (ℓ−1)^{−1} W(Q/ℓ), with ω₀ = log x.
  - (LV8.5) applies with Q′ = Q/ℓ = N^{q−1/(2k)+o(1)} ≥ x^ε, where ε := η₁/(2ξ), and Q′ ≤ √Ψ.
  - Hence Part 1 ≪ Ψ(H(u)^{−c} + B^{−c})·Σ_ℓ ϱ(Q/ℓ)/ℓ + Q√Ψ log^{7/2}x·Σ_ℓ ℓ^{−2}.
  - The first sum is ≤ 2ϱ(Q/2Y)Σ_ℓ ℓ^{−1} ≪ δ_𝒯, using (Hi) and Mertens on (Y, 2Y].
  - The second term is ≪ (Q/Y)√Ψ log^{7/2}x ≤ Ψδ_𝒯x^{−η₂}. This holds because Q ≤ x^{1/2−(ξ/2−q)/2ξ}
    and Ψ, δ_𝒯 ≥ x^{−o(1)}.
- **Part 2: χ_ℓ non-principal.**
  - The conductor is ℓ·r′ ≥ Y = N^{1/(2k)+o(1)}.
  - Choose η := min(1/100, 1/(4kξ)). Then Y ≥ x^η = M₀, so all these characters have conductor in
    [M₀, Q], with moduli t ∈ (Q/2, Q] ⊂ [M₀, Q].
  - By positivity, Part 2 is at most the left side of Fact 8.11, which is
      ≪ log^{7/2}x √Ψ (Q + x^{1/2−η} log²x) ≤ Ψδ_𝒯 x^{−η₂},
    using q < ξ/2 and Ψ, δ_𝒯 ≥ x^{−o(1)}.
- The principal character of t contributes exactly Ψ_t(x)/φ(t), which is the subtracted main term. ∎

*Remarks.*
- Unlike hypD §4.1′, no conductor range [Y, M₀] arises in regime (B), so neither Fact 8.9 nor Fact 8.10
  is needed for χ_ℓ ≠ χ₀.
- Inherited from LV's proof: the Page-character step absorbs (log x)^{O(1)} into exp(−cu/log²u), which
  needs u/log²u ≫ log log x. That holds in regime (B), where u ≍ (log N/log log N)^{1/2}, but not at
  the edge of LV's stated range. We use Lemma 8.5′ only in regime (B).
- In the setup, Fact 3.15 is used with α → 1. Ψ(2x, y) ≤ 2Ψ(x, y) can equally be taken from Fact 3.15,
  up to a factor 1 + o(1).
- ε_N → 0 uniformly in N, and ε_N·(log L)^{O(1)} → 0, since H(u)^{−c} = exp(−c(log N)^{1/2}/
  (log log N)^{5/2+o(1)}).

### 2.3 Theorem FM [proved]

**Statement.**
- Fix odd d ≥ 3, k > (d+1)/4, a > 0, σ > 0 and η₀ > 0, and let h ≥ max(1/(4k), s − (d−2)/(4k)) + η₀.
- Then for all large N (preprocessed), with an implied constant depending only on (d, k, a, σ, h, η₀):
    E_{ℓ∈J, A∼𝒜} N_ℓ(A) ≥ S²ρ_rρ_F/(log N)^{O(1)}.

*Proof.* Write Pr_A for g uniform on 𝒢 and c uniform on C_H, independently. Then
  E N = π_J^{−1} Σ_{ℓ∈J} Σ_{P∈𝒫_ℓ} Σ_{n∈𝒮_B} (|𝒢_P(n)|/|𝒢|)·E_{g∈𝒢_P(n)} Pr_c[F_A(P) ∈ 𝒮_B],
where 𝒢_P(n) := {g ∈ 𝒢 : G_g(P) = ℓn}. Throughout, (log)^{O(1)} denotes (log N)^{O(1)}.

**Step 1 (fibres and the induced law; hypD item 2).** Fix ℓ and P = (a, b) ∈ 𝒫_ℓ, and put V′ :=
max(|a|, b) ∈ (V/2, V].
- (i) The g with G_g(P) = ℓn lie on a line with direction (a, b) and point spacing |P|₂ ≥ V/2. The
  square has chords of length ≤ M/√2. So |𝒢_P(n)| ≤ √2M/V + 1 for every n.
- (ii) Consider the continuous model (g₀, g₁) uniform on [M/2, M]².
  - G = g₁a − g₀b is a sum of independent uniforms on intervals of lengths |a|M/2 and bM/2.
  - Its density is ≥ 1/(MV′) on a central interval of length ≥ MV′/2 ≥ MV/4.
  - Removing |G| < MV/32 leaves an interval I*_P of G-values, of length ≥ MV/16, on which
    |G| ∈ [MV/32, 2MV] and G has constant sign.
  - For t ∈ I*_P the chord {G = t} has length ≥ (M²/4)|P|₂/(MV′) ≥ M/4, so it carries ≥ M/(8V)
    lattice points for large N.
  - Put I_P := I*_P/ℓ, a set of n-values of length ≥ MV/(16ℓ).
- (iii) **Coprimality.** The integer points of the line g = g* + j(a, b) on which e | gcd(g₀, g₁) form
  one residue class of j mod e if e | ℓn, and none otherwise, because gcd(a, b) = 1. By Möbius
  inversion, for n ∈ I_P,
    |𝒢_P(n)| ≥ (φ(ℓn)/ℓn)·M/(8V) − 2^{ω(ℓn)} ≥ M/(C V log log N),
  since 2^{ω(ℓn)} ≤ N^{o(1)} ≪ M/V.

**Step 2 (the algebraic side; hypD items 3, 4).** Fix ℓ, P, n ∈ I_P ∩ 𝒮_B and g ∈ 𝒢_P(n). Put
t := |G(P)| = ℓ|n| and r := F̂_g(P) mod t.
- (a) **(r, t) = 1, so the moving residue is reduced (item 3).**
  - gcd(G, b) = gcd(g₁, b), and g₁a ≡ ℓn (mod b). Since b is a prime > B ≥ P⁺(n) and b ≠ ℓ (b ≤ V < ℓ),
    gcd(g₁, b) = 1.
  - gcd(G, g₁) = gcd(g₀b, g₁) = 1.
  - Lemma 1.2 then gives gcd(r, t) = gcd(N, t) = 1, by preprocessing and ℓ ∤ N.
  - The residue r moves with g and b. It enters only through max_r below.
- (b) **Kernel and ϑ (LV 5.13, 5.18; item 4).**
  - LV's induction in Lemma 5.13 uses only gcd(a, b) = 1. Exchanging the roles of a and b when
    |a| > b, and absorbing signs into the x_i, it gives S_P ⊂ [−2V′, 2V′]^d on which φ_P is a
    bijection onto an interval of length V′^{d−1}. It also gives the sub-kernels ν_l with
    φ_Pν_l = Unif(I(l)).
  - As in Lemma 5.18, with μ = Unif(C_H) and ν = μ ⋆ (⋆_i ν_{l_i}), l_i = |a|^{d−1−i}b^i:
    ‖μ − ν‖₁ ≤ 2·Pr_μ[dist(c, ∂C_H) ≤ 4dV′] ≤ 8d²V/H.
  - Also φ_Pν = ϑ, a convolution of d uniforms on intervals of lengths 2Hl_i, shifted by
    |M₀| ≤ dV^{d−1}.
  - Its density is ≤ (2H max_i l_i)^{−1} = (2HV′^{d−1})^{−1}, and its support lies in
    [M₀ − T′, M₀ + T′] with T′ ≤ dHV^{d−1}.
- (c) **Truncation (LV 5.20; item 4).** Fix ω := L^{C₁}, where C₁ exceeds the exponent of 1/ρ_F.
  - Put Y_ω := (2H + 1)V′^{d−1}/ω. Let ϑ′ be LV's **flattening** of ϑ on the central part:
    ϑ′(z) = ϑ(max(z, Y_ω)) for z ≥ 0 and ϑ(min(z, −Y_ω)) for z < 0. Nothing is deleted. Discrete
    interval lengths are 2H + 1 rather than 2H; this changes only constants.
  - Then ‖ϑ − ϑ′‖₁ ≤ 2/ω.
  - ϑ′ is monotone away from 0, since Y_ω > |M₀| when H > dωV. So ϑ′ = Σ_{y≥Y_ω} (W_y Unif[0, y] +
    W_{−y} Unif[−y, 0)), with W ≥ 0 and |1 − ΣW| ≤ 2/ω.
- (d) **One uniform piece.**
  - Let y ∈ [Y_ω, 2T′]. The values F̂ + G·m, m ∈ [0, y], form a progression ≡ r (mod t).
  - **Sign correction (codex pcest_2).** If G and F̂ have opposite signs, the progression can cross 0;
    its far endpoint then has absolute value ty − |F̂|, not ty. In every case, all terms except at most
    2(|F̂|/t + 1) of them have one sign and absolute value in (|F̂|, ty − |F̂|]. This interval contains
    (ty/4, ty/2], because |F̂| ≤ ty/4.
  - Let X_F^max := C_d MH V_max^d with V_max := ⌊√(2Y)·S⌋, so it is common to all ℓ ∈ J. Let
    𝒳 := {X_F^max ω^{−2}(1 + 1/2)^i} ∩ [X_F^max ω^{−2}, 16dX_F^max].
    - 𝒳 has O(log ω) = O(Λ) elements.
    - It is independent of t, P, ℓ and the residue.
    - It contains a point x₁ ∈ (ty/4, ty/2]. Consecutive points differ by the factor 3/2 < 2, and
      ty/4 ≥ tY_ω/4 ≥ X_F^max ω^{−2} for large N, since tY_ω ≥ X_F^max/(C_d 2^{2d}ω).
  - Call t **good** if, for all x ∈ 𝒳 and all reduced residues r′,
    |Ψ(x, B; t, r′) − Ψ_t(x, B)/φ(t)| ≤ Ψ_t(x, B)/(2φ(t)).
  - For good t, (T) and (Hi) give
      #{smooth terms} ≥ Ψ(x₁; t, ±r) − Ψ(|F̂|; t, ±r) − 2(|F̂|/t + 1)
                      ≥ Ψ(x₁)/(4t) − 3(|F̂|/t + 1) ≥ Ψ(ty)/(16t) − 3(|F̂|/t + 1).
      (The progression contains every integer ≡ ±r (mod t) with absolute value in (|F̂|, x₁].)
  - Dividing by y + 1 and using ϱ(ty) ≥ ϱ(16dX_F^max)/2 ≥ ρ_F/(log)^{O(1)} (by (Hi) and (C8.16)) gives
    Pr ≥ ρ_F/(log)^{O(1)} − O(ωV/H). The error term is polynomially small, since |F̂|/(tY_ω) ≪ ωV/H.
- Combining (b)–(d): **if t is good, then Pr_c[F_A(P) ∈ 𝒮_B] ≥ ρ_F/(log)^{O(1)} − 8d²V/H − 4/ω ≥
  ρ_F/(log)^{O(1)}.**
  - Here V/H ≤ N^{−η₀}L^{O(1)} by the h-floor 1/(4k) + η₀, and 1/ω = L^{−C₁}.
  - This is hypD item 4: additive errors ≪ ρ_F rather than o(1), with truncation ω ≥ L^C.

**Step 3 (bad moduli are rare; hypD item 1).**
- If t is not good, then for some x ∈ 𝒳 and some r′, |E(x; t, r′)| > Ψ_t(x)/(2φ(t)) ≥ Ψ(x)/(4t), by (T).
- So, by Lemma 8.5′ applied at each x ∈ 𝒳 and each of the O(1) dyadic Q covering
  [MV/32, 2MV]·L^{o(1)}:
    #{bad t ∈ 𝒯_Q} ≤ Σ_{x∈𝒳} (4Q/Ψ(x))·CΨ(x)δ_𝒯ε_N ≪ (log ω)·ε_N·#𝒯_Q.
- The hypotheses of Lemma 8.5′ hold:
  - q = s + 1/(4k) + o(1) > 1/(2k) + η₁, because k > (d+1)/4;
  - q < ξ/2 ⟺ h > s − (d−2)/(4k), the level-½ condition, with margin η₀.
- The goodness exponent is a constant 1/2. LV's Definition 5.12 exponent 1 + o(1) is not used, so the
  N^{o(1)} loss that the blind round flagged does not arise.

**Step 4 (assembly; hypD items 2 and 5).**
- *Good mass.* By Step 1(ii)–(iii) and (Hi) with z ≤ 64:
    Σ_{n∈I_P∩𝒮_B} |𝒢_P(n)| ≥ (M/(CV log log N))·Ψ(I_P) ≥ M²ρ_r/(Cℓ log log N).
- *Bad mass.* By Step 1(i), the n with ℓn bad carry total weight ≤ (2M/V)·#{bad n ≤ 2MV/ℓ}. Averaged
  over ℓ, Step 3 gives
    π_J^{−1} Σ_ℓ #{bad n for ℓ} ≤ (log ω)ε_N·π_J^{−1} Σ_ℓ #{n ≤ 2MV/ℓ : n ∈ 𝒮_B} ≪ (log ω)ε_N MVρ_r/ℓ.
  This is o(1/log log N) times the good mass, uniformly in N.
- *Total.* Summing over P ∈ 𝒫_ℓ, with #𝒫_ℓ ≍ V²/log V and |𝒢| ≍ M², and inserting Step 2:
    E N ≥ (log)^{−O(1)}·(V²/ℓ)·ρ_rρ_F = S²ρ_rρ_F/(log N)^{O(1)}.
- *Item 5.* Small b cannot occur, since b > V/2. Irregular lattices contain no candidates (Lemma 1.3).
  The finite-interval error of sampling g in a square rather than mod ℓ is inside Step 1, because the
  fibres are counted exactly. ∎

*Remarks.*
- **What was needed beyond hypD §4.1′.**
  - Restricting b to primes in (V/2, V]. This makes every moving residue reduced and removes small b.
  - The 2D field. This makes the rational side a plain interval count, (Hi), with no BV.
  - A residue-free grid 𝒳 for goodness.
- **The 1D design.** The same proof works for (a′_k), with Step 1 replaced by Harper's Theorem 1 for
  the dense moduli family {b}. That needs √ℓS ≤ (bM/ℓ)^{1/2−η}, i.e. k ≥ (3/4)(d+1)(1 + O(η)). We do
  not write it out, since (a″_k) is the design used below.
- **Precision.** The loss is (log N)^{O(1)}, not merely L^{o(1)}. The losses are 1/log V from prime b,
  log log factors from φ(ℓn)/ℓn, and the (log)^{O(1)} of (C8.16).
  - At Horn 1 (k ≍ (log N)^{1/4}), the number of grid points, ε_N and Lemma 8.5′ change. That regime
    still needs the explicit-error rewrite hypD describes.

## 3. PC_rat: Theorem R

Throughout this section Y_P := 1[G(P)/ℓ ∈ 𝒮_B] is the rational event for a candidate P. The rational
counts are
  R1 := E_{ℓ,g} Σ_{P∈𝒫_ℓ} 1[P ∈ Λ_ℓ(g)]Y_P,
  R2 := E_{ℓ,g} Σ_{P≠P′∈𝒫_ℓ} 1[P, P′ ∈ Λ_ℓ(g)]Y_PY_{P′}.

**Lemma Q (pair image lattice) [proved; tested].**
- Let P = (a, b) ≠ P′ = (a′, b′) be in 𝒫_ℓ with Δ := ab′ − a′b ≡ 0 (mod ℓ). Δ ≠ 0, since b, b′ > 0 and
  both points are primitive.
- Then L_ℓ(P) := {g ∈ Z² : ℓ | G_g(P)} = L_ℓ(P′).
- The map g ↦ (n, n′) := (G_g(P)/ℓ, G_g(P′)/ℓ) is injective on L_ℓ(P), and its image is
    Γ_q = {(n, n′) ∈ Z² : n′ ≡ λn (mod q)},  q := |Δ|/ℓ ≤ 2S²,  λ :≡ b′b^{−1} (mod q).
- Since b is a prime > V/2 > 2S² ≥ q, b is invertible mod q, and λ is a unit.

*Proof.*
1. g ↦ (G_g(P), G_g(P′)) is the integer matrix [[a, −b], [a′, −b′]] acting on (g₁, g₀). Its determinant
   is −Δ ≠ 0, so the map is injective with image of index |Δ|.
2. L_ℓ(P) has index ℓ in Z². ℓ | Δ forces ℓ | G_g(P′) whenever ℓ | G_g(P), so L_ℓ(P) = L_ℓ(P′), and its
   image lies in ℓZ² with index ℓ|Δ| in Z². Dividing by ℓ gives a lattice of index ℓ|Δ|/ℓ² = q.
3. b·G_g(P′) − b′·G_g(P) = g₁(ba′ − b′a) = −g₁Δ, so bn′ − b′n = ∓g₁q ≡ 0 (mod q). The image therefore
   lies in Γ_q, which also has index q, so the two are equal. ∎
- This matches codex pcest_1 (3). The index is q, not ℓq: the special-q condition is absorbed by the
  division by ℓ.

**Lemma DS (divisor switching: smooth numbers in a short interval and a progression) [proved;
elementary].**
- Setup: let C₀ ≥ 1 be fixed, y → ∞ with y ≥ (log X)^K, u = log X/log y → ∞, and 1 ≤ q ≤ y^{C₀}. Let
  I ⊂ [X y^{−C₀}, X] be an interval of length A ≥ X y^{−C₀}.
- Conclusion: for every r ∈ Z,
    #{n ∈ I : n ≡ r (mod q), P⁺(n) ≤ y} ≤ (A/q)·ϱ(X, y)·(log X)^{O_{C₀}(1)}.

*Proof.*
1. Let e := gcd(r, q). If e has a prime factor > y, the count is 0. Otherwise write n = en₁. Then
   n₁ ∈ I/e has length A/e and n₁ ≡ r/e (mod q′) with q′ = q/e and gcd(r/e, q′) = 1, so n₁ is coprime
   to q′.
2. Put D := min(A/q, X/(e y^{C₀})). Then D ≤ min(I/e).
3. Every y-smooth n₁ ≥ D has a divisor d with D/y < d ≤ D: take the largest product of its smallest
   prime factors that is ≤ D.
4. Write n₁ = dm. Then gcd(d, q′) = 1, m ≡ (r/e)d^{−1} (mod q′), and m lies in an interval of length
   A/(ed). So there are ≤ A/(dq) + 1 ≤ 2A/(dq) choices of m, since A/(Dq) ≥ 1.
5. Hence the count is ≤ (2A/q) Σ_{D/y<d≤D, P⁺(d)≤y} 1/d.
6. By dyadic blocks and (Hi), the sum is ≤ 8 log₂(2y)·ϱ(D/y, y).
7. Finally X/(D/y) ≤ y^{2C₀+1}, so (C8.16) gives ϱ(D/y, y) ≤ ϱ(X, y)(log X)^{O(C₀)}. ∎
*Attribution.*
- The case q = 1 follows from Hildebrand's subadditivity Ψ(x + z, y) ≤ Ψ(x, y) + Ψ(z, y) (Quart. J.
  Math. 36 (1985) 57–69) together with (C8.16).
- Choosing a divisor in (D/y, D] is the standard flexible factorization of smooth numbers (Friedlander,
  Acta Arith. 39 (1981); Vaughan's iterative method).
- Codex pcest_2 points to Harper, Compositio 152 (2016), §2.1, "Smooth Numbers Result 3", as a stronger
  published progression version, including non-coprime residues. We have not checked its exact ranges
  against ours, so we use our own proof above.
- Lemma DS is therefore a progression extension of known short-interval bounds, not a new result.
- Codex's blind round states the same bound (its (1)).

**Theorem R (rational pair correlation) [proved].** In design (a″_k), with B = L[1/2, a] and S = L^σ
satisfying S² ≥ ρ_r^{−2}:
  R2 ≤ (log N)^{O(1)}·R1²  and  R1 ≥ S²ρ_r/(log N)^{O(1)}.

*Proof.*
- **The lower bound on R1** is Steps 1 and 4 of §2.3 without the algebraic factor.
- **The pair count.** Fix ℓ and a pair (P, P′) as in Lemma Q. Then
    #{g ∈ 𝒢 ∩ L_ℓ(P) : Y_PY_{P′} = 1} ≤ #{(n, n′) ∈ Γ_q ∩ Π : n, n′ ∈ 𝒮_B}.
  Here Π := ℓ^{−1}[[a, −b], [a′, −b′]]·[M/2, M]² is a parallelogram of area M²q/(4ℓ), contained in
  [−X_r, X_r]².
- **Shape of Π.** One side is (M/2ℓ)(b, b′), of length ≍ MV/ℓ and with direction in
  [atan(1/2), atan 2], because b, b′ ∈ (V/2, V]. So Π has width ≥ area/(2·long side) ≥ Mq/(8√2V) ≥
  w₀/12, where w₀ := M/V.
- **Covering.** Cover Π by the grid squares of side w₀ that meet it. There are ≤ C·area(Π)/w₀² = C·qS²
  of them.
- **A square away from the axes.** Let 𝒬 = I₁ × I₂ have |n|, |n′| ≥ X_r S^{−2} throughout.
  - Lemma DS applies to each coordinate (sign symmetric; I_i ⊂ [X_rS^{−2}, X_r];
    A = w₀ = X_r/(2S²) ≥ X_r y^{−C₀}; q ≤ 2S² ≤ y^{C₀}).
  - Hence #(Γ_q ∩ 𝒬 ∩ 𝒮_B²) = Σ_{r mod q} #{n ∈ I₁ : n ≡ r}·#{n′ ∈ I₂ : n′ ≡ λr}
    ≤ q·((w₀/q)ρ_r(log)^{O(1)})² = (w₀²/q)ρ_r²(log)^{O(1)}.
  - Summing over the ≤ CqS² squares gives ≤ CS²w₀²ρ_r²(log)^{O(1)} = C(M²/ℓ)ρ_r²(log)^{O(1)}.
- **Near the axes.** For a square that meets the band |n| < X_rS^{−2}, count only its points that lie
  in Π (the original points). They have |n| < 2X_rS^{−2}, since w₀ = X_r/(2S²). Each n-value has
  ≤ √2M/V + 1 preimages g (Step 1(i)). So there are at most 8MX_r/(VS²) = 16M²/(ℓS²) such points; the
  same holds for n′.
  - For the filtered count of §1 (|G| ≥ MV/32, i.e. |n| ≥ X_r/64), this band is empty.
  - Squares crossing the boundary of [−X_r, X_r]² are handled by applying Lemma DS with X = 2X_r.
  - With the choice S² = L^ε/(ρ_rρ_F) of §5, and since ρ_F ≤ ρ_r, we get S^{−2} ≤ ρ_r², so this is
    absorbed. In general Theorem R needs S² ≥ ρ_r^{−2}.
- **Per pair**, the count is ≤ (M²/ℓ)ρ_r²(log N)^{O(1)}. This bound is independent of q.
- **Counting pairs.** For fixed P, the P′ with ℓ | Δ are the points of Λ_ℓ(g) ∩ 𝒫_ℓ for g = (a, b),
  i.e. g₀/g₁ ≡ a/b. This lattice contains P ∈ 𝒫_ℓ, so it is regular. There are
  ≤ 22S² + 1 of them (Lemma 1.3). Hence
    R2(ℓ) ≤ #𝒫_ℓ·23S²·|𝒢|^{−1}(M²/ℓ)ρ_r²(log)^{O(1)} ≍ (V²/log V)(S²/ℓ)ρ_r²(log)^{O(1)}
          = S⁴ρ_r²(log)^{O(1)}.
- **Conclusion.** By the lower bound, R1(ℓ) ≥ S²ρ_r/(log)^{O(1)} for every ℓ ∈ J. Averaging over ℓ
  proves the theorem. ∎

*Remarks.*
1. **Why it works in 2D and not in 1D.** In (a′_k), the rational pair (s, s′) of a pair moves on the
   *line* b′s − bs′ = Δ/ℓ as m varies (hypD §4.2).
   - Divisor switching then sends the problem to (e, e′) solutions of a single linear equation in a
     box. That count is a Diophantine 0/1 quantity with no area main term, and the method gives
     nothing.
   - Averaging over (b, b′, h) does not repair this: the averaged weight is concentrated on pairs with
     good rational approximations b′/b to s′/s.
   - With (g₀, g₁) random, the pair fills a 2D lattice region of area ≫ perimeter·(lattice scale),
     which is all the argument needs. **This is the "key simplification" the prompt asked about. It is
     not a sieve bound: sieve and Rankin bounds fail at density L^{−c}, see below. It is divisor
     switching, made possible by the 2D randomization.**
2. **Sieve and Rankin alternatives (prompt item 3, checked).**
   - Shiu and Nair–Tenenbaum upper bounds for Σ F(n)F(n + h) give only ∏_{B<p≤X}(1 − 2/p) ≍
     (log B/log X)², a logarithmic saving against the needed ρ².
   - A Rankin tilt n^{α−1} that recovers the Dickman exponent cannot be localized to an interval or a
     line: extended to all n it is dominated by small n, and a class mod q receives mass ≥ 1 > L(χ₀)/q
     once q > e^{Cu}.
   - Matomäki–Radziwiłł–Tao-type averaged correlations have o(1) or log-power errors, useless at
     density L^{−c}.
   - Pascadi's Cor. 7.2 gives exponent 13/8 − ε < 2. Codex notes that its coefficient condition can
     hold for large fixed k, but the exponent still falls short.
   - So, as far as we can tell, no sieve-type upper bound with loss L^{o(1)} is available in 1D; codex
     reached the same conclusion.
3. **Precision.** The loss is polylogarithmic, so Theorem R also meets Horn 1's precision requirement
   exp(−o(Λ/k)) for k ≍ (log N)^{1/4}. The remaining Horn-1 obstacles are FM's explicit errors and LD.

## 4. PC_alg: the kernel transfer, the remaining hypothesis, and why known tools fall short

### 4.1 Lemma Kern (the 2D kernel lemma) [proved; tested]

**Setup.** Let d ≥ 3, and let P ≠ P′ be primitive with |P|, |P′| ∈ [V/2, √2V] and Δ = det(P, P′) ≠ 0.
Put v = v(P), v′ = v(P′), and Φ: Z^d → Z², c ↦ (⟨c, v⟩, ⟨c, v′⟩). Write P′ ≡ κP (mod Δ); κ exists
because P is primitive.

**Statement.**
- (i) Γ := Φ(Z^d) has index |Δ| and equals {(X, X′) : X′ ≡ κ^{d−1}X (mod Δ)}.
- (ii) ker Φ = Q_{PP′}·Z[x, y]_{d−3}, where Q_{PP′} = (bx − ay)(b′x − a′y). This kernel is saturated, and
  its monomial basis has ‖·‖_∞ ≤ 2V² on 𝒫_ℓ, where |a|, |b| ≤ V. In general the bound is
  ≤ 4V² for |P| ≤ √2V, and then H ≥ 2V² is required; only constants change.
- (iii) If H ≥ V², then for every z ∈ Z²,
    Pr_{c∈C_H}[Φ(c) = z] ≤ C_d·|Δ|/(H²‖v ∧ v′‖) ≤ C′_d/(H²V^{2d−4}).

*Proof.*
- **(ii).** ⟨c, v(P)⟩ = C(P) for the form C = Σ c_i x^{d−1−i}y^i. C(P) = 0 with P primitive gives
  (bx − ay) | C over Z, by Gauss's lemma. Then C(P′) = 0 and det(P, P′) ≠ 0 give (b′x − a′y) | C/(bx − ay).
  So ker Φ = Q·Z[x, y]_{d−3}. A multiple of Q that is integral is an integral multiple, again by Gauss's
  lemma, so the kernel is saturated.
- **(i).** The index of Φ(Z^d) is the gcd of the 2×2 minors of [v; v′], which is |Δ| (hypD Prop G(i),
  (I6); also codex pcest_1, via a unimodular change sending P to (1, 0)). Since P′ ≡ κP (mod Δ),
  v′ ≡ κ^{d−1}v (mod Δ), so Γ lies in the displayed lattice, which also has index |Δ|. Hence they are
  equal.
- **(iii).**
  - A fibre is (c₀ + K) ∩ C_H, with K = ker Φ of rank d − 2 and det K = ‖v ∧ v′‖/|Δ|. Here we use that
    a saturated lattice and its orthogonal complement in Z^d have the same determinant, and that the row
    lattice of [v; v′] has index |Δ| in its saturation.
  - The translates c + F, where F is the fundamental cell of the basis in (ii), are disjoint. They lie
    in the section of the cube [−(2d−3)H, (2d−3)H]^d by a (d−2)-flat, whose volume is ≤ C_dH^{d−2}.
  - So the fibre has ≤ C_dH^{d−2}|Δ|/‖v ∧ v′‖ points. Divide by (2H + 1)^d.
  - For the second inequality: every 2×2 minor of [v; v′] equals Δ·W_{ij}(P, P′), with W bihomogeneous
    of bidegree (d−2, d−2). For example W_{01} = (aa′)^{d−2}.
  - The vector W(P, P′) of all divided minors never vanishes on (S¹)² (repaired after codex pcest_2).
    - Off the diagonal and antidiagonal, this is injectivity of the projective Veronese map: v(P) ∦ v(P′)
      when P ∦ P′, and Δ ≠ 0.
    - On the diagonal, W(P, P) is the wedge of v(P) with its tangent ∂_u v(P) divided by det(P, u). It
      is nonzero because the rational normal curve is nondegenerate. For example, at P = (1, 0) the
      adjacent divided minor W_{01} equals 1. The antidiagonal is the same up to sign.
    - W is continuous, so by compactness ‖v ∧ v′‖ ≥ c_d|Δ|(|P||P′|)^{d−2}. ∎
    - **Effective version for candidates (Opus proof critic).** The minor (d−2, d−1) equals
      (bb′)^{d−2}Δ. On 𝒫_ℓ we have b, b′ > V/2, so ‖v ∧ v′‖ ≥ |Δ|(V/2)^{2d−4} directly, without any
      compactness argument. This is all that Prop 4.2 uses.
    - The exact identity det K·|Δ| = ‖v ∧ v′‖ and the bound with c_d > 0.05 (d ≤ 7) are checked exactly
      in `tests/test_pcest.py`.

*Cost [proved arithmetic].*
- H ≥ V² = ℓS² means h > 1/(2k), which raises the h-floor from LV's 1/(4k) (FM) to 1/(2k).
- Then θ = 2s + (d−1)/(4k) + max(1/(2k), s − (d−2)/(4k)). The floor binds only for k < d(d+1)/4,
  against (d² − 1)/4 before.
- **The infimum over k is unchanged.**
- Finite-k examples: d = 9, k = 20 gives θ = 0.325 instead of 0.3125, so c = 0.937 instead of 0.928. At
  d = 9 the minimal k with c < 1 moves from 11 to 12 (`src/pcest/params.py`). Codex's (5) agrees.
- If one also wanted the moduli G|Δ to stay below |F|^{1/2} (not needed below), the condition would
  become h > s + (6 − d)/(4k) (codex).

### 4.2 Reduction of PC to one arithmetic hypothesis

**Notation.** For a pair (P, P′) ∈ 𝒫_ℓ² with ℓ | Δ, and g ∈ 𝒢 ∩ L_ℓ(P):
- G = G_g(P), G′ = G_g(P′), and T := dHV^{d−1};
- the full windows I := F̂_g(P) + G·[−T, T] and I′ likewise;
- the lift classes r_j := F̂_g(P) + Gj (mod G|Δ|) and r′_j := F̂_g(P′) + G′κ^{d−1}j (mod G′|Δ|), for
  j mod |Δ|.

**The normalized lift correlation.**
  κ_alg(ℓ, P, P′, g) := |Δ|·Σ_{j mod |Δ|} Ψ_I(G|Δ|, r_j)·Ψ_{I′}(G′|Δ|, r′_j) / (4T²·ϱ(T|G|, B)·ϱ(T|G′|, B)),
where Ψ_I(m, r) := #{F ∈ I ∩ 𝒮_B : F ≡ r (mod m)}. If the smooth numbers of each window were spread
evenly over the |Δ| lift classes, κ_alg would be ≍ 1.

**Prop 4.2 [proved].** Assume h ≥ 1/(2k) + η₀. For each such (ℓ, P, P′, g),
  Pr_c[F_A(P), F_A(P′) ∈ 𝒮_B] ≤ C_d·(V²/|Δ|)·ϱ(T|G|)ϱ(T|G′|)·κ_alg.
*Proof.*
- By Lemma Kern(iii), Pr_c ≤ C′_d(H²V^{2d−4})^{−1}·#{z ∈ Γ ∩ [−T, T]² : F̂ + Gz₁ and F̂′ + G′z₂ smooth}.
- By Lemma Kern(i), that count is Σ_j Ψ_I(G|Δ|, r_j)Ψ_{I′}(G′|Δ|, r′_j).
- Insert the definition of κ_alg and T² = d²H²V^{2d−2}. ∎

Here V²/|Δ| = S²/q is comparable to (not equal to) the bounding-box factor 1/sin∠(P, P′).
- If κ_alg were ≤ L^{o(1)}, the factor would cost only O(log S) on average. For fixed P and q, the
  equation det(P, P′) = ±ℓq has O(1) solutions P′ in the box (they lie on two lines parallel to P,
  spaced by |P| ≥ V/2).
- So Σ over pairs of (S²/q)·(Theorem R's uniform per-pair bound) loses only a harmonic sum
  Σ_{q≤2S²} 1/q.
- LD keeps the product (S²/q)·κ_alg; the two factors are not separated.

**Hypothesis LD (lift dispersion under rational-smoothness weights) [open].** Uniformly in N:
  E_ℓ Σ_{P≠P′, ℓ|Δ} |𝒢|^{−1} Σ_{g∈𝒢∩L_ℓ(P)} Y_PY_{P′}·(S²/q)·κ_alg ≤ L^{o(1)}·E_ℓ Σ_{P≠P′, ℓ|Δ} |𝒢|^{−1} Σ_g Y_PY_{P′}.
- It is an upper bound, at precision L^{o(1)}, for a restricted second moment of smooth numbers in the
  lift classes of moving residues to moduli G|Δ| = ℓ²nq.
- The weights are the actual rational-smooth multiplicities. Theorem R shows that the right side is
  ≤ (log N)^{O(1)}S⁴ρ_r².
- LD is implied by a suitably weighted form of hypD's RBDH′ (an o(1) relative variance), but it is
  *weaker*: it is an upper bound within L^{o(1)}. It is codex pcest_1's (6) with the kernel transfer
  carried out. We do not claim an equivalence.

**Theorem 4.3 [proved implication; repaired after codex pcest_2].**
- **Hypotheses.** Those of Theorem FM, with h ≥ 1/(2k) + η₀ and S² := (log N)^{C₃}/(ρ_rρ_F), where C₃ is
  larger than the polylog exponent lost in Theorem FM. Let N_ℓ(A) be the **filtered** count of §1
  (|G(P)| ≥ MV/32).
- **Conclusions.**
  - (i) FM + Theorem R + Prop 4.2 + LD ⟹ PC := E[N(N−1)]/(E N)² ≤ L^{o(1)}.
  - (ii) Hence, by Prop PZ (E N ≥ 1), D_k(γ) holds with γ ≥ L^{−o(1)}.
*Proof.*
- E[N(N−1)] ≤ E_ℓ Σ_{P≠P′, ℓ|Δ} |𝒢|^{−1} Σ_{g: |G|,|G′| ≥ MV/32} Y_PY_{P′} Pr_c[·].
- By the filter, T|G| ∈ [TMV/32, 2TMV], so ϱ(T|G|) ≤ 64ϱ(TMV) ≤ ρ_F(log N)^{O(1)} by (Hi) and (C8.16).
  - The unfiltered version of this step was wrong: for |G| ≈ ℓ, the ratio X_F/(T|G|) is a positive
    power of N.
- By positivity, LD (stated over all g) bounds the filtered sum.
- So, by Prop 4.2, LD and Theorem R, E[N(N−1)] ≤ C_dρ_F²·L^{o(1)}·S⁴ρ_r²(log N)^{O(1)}.
- By Theorem FM, which only counts |G| ∈ [MV/32, 2MV] and so applies to the filtered count,
  (E N)² ≥ S⁴ρ_r²ρ_F²/(log)^{O(1)} and E N ≥ 1.
- Theorem R's condition S² ≥ ρ_r^{−2} holds for large N, since ρ_F/ρ_r = L^{−c+o(1)} for some c > 0. ∎
- (The first version took S² = L^ε/(ρ_rρ_F), giving L^{O(ε)} losses that are not L^{o(1)}. The
  polylogarithmic choice of S² fixes this.)

### 4.3 Why LD is not available, and what each known tool gives

All items below are statements about what the cited tools give by themselves. They are not lower bounds
on LD, and they do not exclude proofs that combine these tools with further arithmetic information.

- **(B1′) Naive use of Harper's BDH (all residues, arbitrary residue) loses the factor G.** For one
  modulus m = G|Δ|:
  - LD needs the |Δ| lift classes, with main-term mass |Δ|(Ψ/m)² = Ψ²/(G·m).
  - Harper's Theorem 2 (Q ≤ Ψ) controls Σ_{m≤Q} Σ_{all reduced r}|E(m, r)|², a global sum whose
    per-modulus scale is Ψ²/m.
  - Harper's asymptotic BDH for general sequences (J. London Math. Soc. 112 (2025) e70293,
    arXiv 2412.19644, Cor. 1) confirms that the all-residue variance is dominated by its diagonal QΨ.
    It has no restricted-moduli or restricted-residue form either.
  - Positivity cannot recover the factor G = N^{s+1/(4k)+o(1)}, and e^{−cu/log²u} cannot compensate a
    polynomial factor (codex pcest_1, Q2). The modulus family is also sparse, with density
    ≈ ρ_r/(ℓ log Y), so restricted versions are needed anyway, as in barrier B1.
- **(B2′) Divisor switching needs small moduli.** Lemma DS works because q ≤ y^{O(1)}.
  - Applied to F in a progression mod G ≈ N^{s+1/(4k)}, the complementary-factor count is only
    nontrivial when the divisor cutoff is D ≲ X_F/G.
  - The divisor sum is then evaluated at density ϱ(X_F/(Gy)) instead of ϱ(X_F), a loss of L^{Θ(1)}.
  - The q-part of the coupling (modulus q ≤ 2S²) is small. But it sits inside an AP to the large
    modulus G, so the same obstruction applies.
- **(B3′) The arithmetic large sieve loses 1/(ρ_Fδ) [assessment].** Use the large sieve in its standard
  form, over primitive characters, or over primes p ≤ Q with Σ_p p Σ_w |Z(p, w) − Z/p|² ≤ (x + Q²)Z for
  general moduli via induced-character multiplicities. It is insensitive to the sparseness of smooth
  numbers.
  - Its diagonal term x·Ψ exceeds the needed scale by the factor 1/(ρ_Fδ) = L^{Θ(1)}, where δ is the
    density of the modulus family.
  - The same phenomenon is why Fouvry–Tenenbaum's large-sieve BV works only for y ≥ exp(log^{2/3+ε}x)
    (Harper's introduction).
- **Prop B5 (an L¹ saving alone does not give L² control) [proved; elementary; scope narrowed after
  codex pcest_2].**
  - The max-residue estimates we know are L¹ bounds over moduli of β − 1, with saving
    ε_N = H(u)^{−c} + …. These are Harper Thm 1, LV Lemma 8.5 and our Lemma 8.5′, and, with fixed
    residues, de la Bretèche–Drappeau, Drappeau and Pascadi.
  - The proof of LD from β-information needs an L² bound: κ_alg ≤ β(G|Δ)·β(G′|Δ)-type products, and
    neither factor is controlled pointwise.
  - Example: let β = T with probability ε_N/T, and β = 1 otherwise, with T ≤ 1/ρ_F (the trivial range).
    The L¹ excess is ε_N, but the L² excess is ε_NT, which reaches ε_N/ρ_F = L^{Θ(1)}.
  - So no inequality of the form "L¹ saving ⟹ L² control" holds at these densities in general.
  - A proof of LD that goes through max-class ratios needs, in addition, a *large-values* estimate: the
    moduli with an overloaded lift class must have weighted measure ≪ L^{o(1)}T^{−2}, for T up to
    1/ρ_F.
  - This does not show that every LD proof must take this form.
- **Prop Z (the ℓ-part is one shared random variable) [proved; observation].**
  - For P ∈ Λ_ℓ(g), P ≡ (b/g₁)(g₀, g₁) (mod ℓ). Hence ⟨c, v(P)⟩ ≡ (b/g₁)^{d−1}R_c(g₀, g₁) (mod ℓ).
  - So the ℓ-adic lift of every candidate F_A(P) mod ℓ|G| is determined by one field constant
    ζ := R_c(g₀, g₁)g₁^{1−d} (mod ℓ), through ζb^{d−1}.
  - Conditioning PZ on ζ removes the ℓ-part of the coupling. But it then needs anti-concentration of
    E[N | ζ] in ζ. That is again a lift-equidistribution statement of LD's type [assessment; not a
    proved equivalence].
  - It does localize the difficulty: equidistribution of smooth numbers in an AP mod ℓn over the ℓ
    lifts mod ℓ²n, jointly with a mod-q coupling.
- **A plausible route to LD (Opus proof critic; not carried out) [assessment].**
  - B1′ and B5 treat the residue as arbitrary. But by Lemma 1.2 the moving residue is explicit:
    r₀ ≡ N b^d g₁^{−d} (mod G).
  - At a fixed modulus t, there are about t/log V pairs (P, g), and the fractions b/g₁ mod t come from a
    box of area about t. So the coarse classes are spread over a coset of the d-th powers, of index
    ≤ d^{ω(t)} = L^{o(1)} on average.
  - LD would then follow from two pieces:
    - (i) an all-residue L² restricted BDH for the moduli ℓ²nq, with LV-type relative saving ϱ(Q). This
      is the L² analogue of LV's majorant argument, Σ φ^{−1} Σ |Ψ(χ)|², which LV say their method
      "naturally provides", and Lemma 8.5′'s character factorization should carry over;
    - (ii) a large-values bound for the structured residues r ≡ N(b₀/g₀)^d with small b₀, g₀, where
      the representation count is ≈ K on a set of measure ≈ 1/K. This is where a real obstruction
      remains.
  - Heuristically LD is plausible. The obvious correlation, shared divisibility by p | q, costs
    ∏_{p|q}(1 + p^{1−2α}) = L^{o(1)}. E1's PC_alg of 1.05–1.27 is consistent with it.
- **Two algebraic sides (codex).** Independent fields do not decorrelate two evaluations in one field.
  The formal constant worsens to (1 + 4/(d+1))/√2, and an analogue of LD is still needed.

**Sufficient single-family form [outline, not a proved reduction].**
- By Cauchy–Schwarz over j, by bounding the moving class by its maximum over classes mod G (which
  depends only on n), and by Lemma DS in the n′ variable (Π's long side makes an angle in
  [atan ½, atan 2] with the axes), LD would follow from LV*:
    Σ_{n∈(X′, 2X′], n∈𝒮_B} max_{(r,ℓn)=1} |Δ|·Σ_{j mod |Δ|} Ψ_I(ℓn|Δ|, r + ℓnj)² / (2Tϱ_F)² ≤ L^{o(1)}·Ψ(2X′, B)
  on average over ℓ ∈ J and |Δ| = ℓq with weight 1/q, for X′ ∈ [X_rS^{−2}, X_r].
- LV* is a "max over coarse classes, L² over the fine lifts" statement for the single family {ℓ²nq}.
  We have not written this reduction out in full.

## 5. Assembly: the sharpest conditional theorem

**Assumption (LV6) [open; inherited from lemmaM Corollary DD and hypD Corollary K].** LV §6
(characters, square roots, the residual square-class kernel) has to be adapted to:
- pure-descent rows from pairs of height ≤ V = N^{1/(4k)}L^{O(1)};
- the non-monic pair (f, g₁x − g₀);
- the one extra g₁-parity column.
The adaptation must keep the cost (rows)^{2+o(1)} and the transcript invariance (T) of Theorem DD.
LV's Theorem 2.6 is written for heights L(1/3). The adaptation has not been written out, by lemmaM,
hypD or us.

**Theorem A [conditional on LD and LV6].**
- A possibly cheaper alternative to A′, not checked: Bernstein's deterministic batch smoothness
  algorithm, run over batches of ≥ B/polylog candidates gathered across R-trials and fields, might
  remove HECM at no cost to the exponent (Opus proof critic).
- **Statement.**
  - Fix odd d ≥ 3, k > (d+1)/4, η₀ > 0, and let θ := θ(d, k) = 2s + (d−1)/(4k) + max(1/(2k),
    s − (d−2)/(4k)) with θ < 1.
  - Let h be the corresponding height plus η₀, let a := θ/√2, and B = L[1/2, a].
  - Assume LD for design (a″_k) at these parameters, uniformly in N, and LV6.
  - Then every N is factored in expected time L_N[1/2, (1+θ)/√2 + O(η₀) + o(1)].
- **Consequences.** Infimum over k: (1 + 3/(d+1))/√2, i.e. 0.9723 (d = 7), 0.9192 (d = 9), 0.8839
  (d = 11) and 0.8132 (d = 19). The constant is < 1 iff d ≥ 7, with k ≥ 13 (d = 7) or k ≥ 12 (d = 9).
- **Further explicit assumptions (codex pcest_2).**
  - (HECM) LPP's running-time theorem for the hyperelliptic curve method, as used by LV Cor. 3.19. The
    published proof trail is incomplete (LPP II defers details to an unpublished Part III); Theorem A′
    removes this assumption.
  - Reducible f are rejected. The irreducibility transfer omits the O(log N) primes dividing g₁N; the
    exceptional probability is polynomially small.
  - The g₁ column records the parity of the number of atomic special-q relations in an assembled row,
    counted with multiplicity.
- **Ingredients and their status.**
  - Proved: Lemma GB, Theorem FM, Theorem R, Lemma Kern, Prop 4.2, Theorem 4.3.
  - From hypD and lemmaM: Theorems DD, S′, K and Corollary K's algorithm.
  - Prior art: Lenstra–Pomerance as the shape test, and HECM for the B-smoothness tests (caveat below).

*Proof.*
1. Preprocess by trial division to B.
2. Theorem 4.3 gives D_k(γ) with γ ≥ L^{−o(1)}.
3. hypD Theorem K(i)–(iii) and its field-selection rule give pure-descent rows at expected cost
   (T_shape + T_SQ)·L^{o(1)} each, with T_shape = L_N[1/2, 1/√2] (LP on √N-size integers) and
   T_SQ = S²L^{o(1)} = L^{θ/(2a)+o(1)}.
4. B·L^{o(1)} rows, LV6 and Theorem DD give a nontrivial factor with probability ≥ 1/2 per dependency.
5. Cost: L[1/2, a + max(1/√2, θ/(2a))] + L[1/2, 2a]. At a = θ/√2 this is L[1/2, (1+θ)/√2], since
   θ ≤ 1. ∎

**Theorem A′ (no HECM) [conditional on LD and LV6].**
- Replace HECM by the deterministic Pollard–Strassen method, which finds all prime factors ≤ B of m in
  B^{1/2+o(1)}(log m)^{O(1)} time (Pollard, Proc. Camb. Phil. Soc. 76 (1974); Strassen, Jber. DMV 78
  (1976); Costa–Harvey, Math. Comp. 83 (2014)).
- Each smoothness test then costs L^{a/2}, and
    c_PS = min_a max(2a, a + max(1/√2, θ/(2a) + a/2)).
- For θ ≤ 1/2 the optimum is at a² − √2a + θ = 0, giving c_PS = 1/√2 + (√2 − √(2 − 4θ))/2.
- θ → 3/(d+1) gives d = 7: 1.0607; d = 9: 0.9670 (a = 0.2599); d = 11: 0.9142; d = 19: 0.8226. These
  were corrected after codex pcest_2; the first version printed 0.9671 and 0.9144 from a grid search.
- So **the HECM caveat (LPP Part III unpublished) can be removed at a cost, and c < 1 survives for
  d ≥ 9**.
- A possibly cheaper route, not checked: Pomerance's (1987) rigorous ECM detection of *most* smooth
  numbers. Whether a "most smooth values" guarantee suffices for FM and PZ (literature critic) is open. LV6 might itself use HECM; if so, the same replacement applies there at L^{o(1)}-many tests
  per row. We did not check LV §6 for this.

**What exactly is missing for a theorem beating Lenstra–Pomerance [the sharpest form we can state].**
1. **LD** (§4.2): a restricted, rational-smoothness-weighted BDH upper bound at precision L^{o(1)} for
   smooth numbers in the |Δ| lift classes of a reduced residue to moduli ℓ²nq. Here n is B-smooth,
   q ≤ 2S², the moduli are ≍ N^{s+3/(4k)}L^{O(1)}, and the values are ≍ N^{s+h+d/(4k)}.
2. **LV6** (bookkeeping of a known method at new heights; not expected to need a new idea, but not
   written).

Within this project's proved results, nothing else is missing. FM, the binary rational correlation,
the kernel transfer, the split counting (Theorem S′), the complete-target-list reduction (Theorem K)
and fruitfulness (DD) are proved here or in hypD/lemmaM.

Caveats:
- these are internal, unrefereed results;
- the HECM running-time theorem (Theorem A only; A′ avoids it);
- the Lemma 5.1 irreducibility transfer, checked at the level of statements.

**Has anyone beaten c = 1 rigorously (1992–2026)?** To our knowledge, no. Codex pcest_1 and pcest_2 and
the Opus literature critic (§7) all searched. Neither found a rigorous algorithm factoring every N in L[1/2, c] with c < 1. The
rigorous record remains Lenstra–Pomerance 1992 (JAMS 5, 483–516).
- Related results are either heuristic (NFS; Stange 2022) or concern other problems:
  - Kleinjung–Wesolowski: discrete logarithms in finite fields of fixed characteristic (JAMS 35, 2022).
  - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588): ERH, class and unit groups.
- Under GRH the best rigorous constant is also c = 1: Seysen and A. K. Lenstra's class-group method,
  whose GRH dependence Lenstra–Pomerance removed. We know of no c < 1 under a standard hypothesis.
- Lee–Venkatesan (J. Number Theory 187 (2018) 92–159) produce congruences in L[1/3], but their
  factoring conclusion is
  conditional on their Conjecture 7.1.
- This is a search result, not a certification.

## 6. Experiment E1 (pre-registered; exact; toy)

Command (as registered):
`.venv/bin/python -m src.pcest.e1 experiments/configs/pcest_e1.json experiments/raw/pcest_e1.jsonl 8`,
then `--analyze` (summary in `experiments/raw/pcest_e1_summary.json`).
- Setup: 8 instances (semi/tri × 56/64 bits × seeds 0, 1), 40 fields × 30 primes ℓ ∈ (2¹¹, 2¹²] per
  design, d = 3, S = 16, B_rat = 256, B_alg = 16384, H = 1024.

| prediction | registered | result |
|---|---|---|
| P1 PC_rat(2D) | [0.8, 1.5] in ≥ 7/8 | 1.007–1.013 in 8/8 — holds |
| P2 PC(2D) | [0.8, 2.0] in ≥ 7/8 | 1.059–1.282 in 8/8 — holds |
| P3 \|PC_rat(2D) − PC_rat(1D)\| | ≤ 0.25 in ≥ 7/8 | ≤ 0.003 in 8/8 — holds |
| P4 E N(2D)/E N(1D) | [0.7, 1.4] in ≥ 7/8 | 0.766–0.938 in 8/8 — holds |
| P5 exact identities | 0 violations | 0 identity violations; 160/160 Lemma Q checks pass (see disclosure) |
| P6 (reported) | — | γ(2D) 0.68–0.89 vs γ(1D) 0.75–0.92; PC_alg(2D) 1.05–1.27 vs PC_alg(1D) 1.06–1.14 |

**Disclosure (P5).** The harness counts "G = 0 or ℓ ∤ G" as a violation, and flagged 1 event in each of
3 one-dimensional records.
- In each case m ≡ 0 (mod ℓ). Then the lattice {a ≡ 0 (mod ℓ)} is irregular (λ₁ = 1), and its box
  contains the root point (m, 1), where G = a − m b = 0.
- These are degenerate points of irregular lattices, which the theory excludes (Lemma 1.3), and the
  code skips them. They are not failures of ℓ | G or of Lemma 1.2.
- Factor discoveries: none were possible, by construction. Every prime factor of the instances
  (≥ 2¹⁸) exceeds ℓ and both smoothness bounds, so gcd(G, N) = gcd(F, N) = 1 for every tested smooth
  value. STRUCTURAL and ACCIDENTAL rates are both 0, with null expectation 0.

**Reading.**
- The 2D design yields about 14% fewer pairs. Its rational values g₁a − g₀b are larger by a constant
  factor than a − bm.
- Its correlations are indistinguishable from the 1D design's at toy size.
- PC_rat ≈ 1.01 reflects large E N_rat (130–190) at these bounds. The toy is not informative about the
  asymptotic difficulty of PC_rat in 1D, which is a question of provability, not of size.
- PC_alg ≈ 1.1–1.3 is consistent with hypD E2 (1.08–1.21). This is weak support that LD is true; it is
  not evidence for any asymptotic statement.

## 7. Reviews

**Status.**
- GPT-6.1-sol adversarial round (`codex/pcest_2.out`): done, all findings applied (§7.2).
- Opus literature critic: done, all fixes applied (§7.3).
- Opus proof critic: done, all fixes applied (§7.4). It reviewed the version before the pcest_2 repairs.
  Its one major proof finding, the small-|G| gap in Theorem 4.3, coincides with codex's and was already
  repaired.
- GPT-6.1-sol pcest_3: adversarial pass on the *repaired* Theorem 4.3 and Theorem R (§7.5).

### 7.2 GPT-6.1-sol adversarial round (pcest_2): adjudication (all accepted)
- **Verdicts.**
  - CORRECT: Lemma GB, Lemma 1.2, Lemma 8.5′ (the majorant W use is confirmed), Lemma Q, Lemma DS,
    Prop 4.2, Lemma Kern (i)–(ii), the det K identity, and the minimal-k numbers.
  - CORRECT WITH FIXES: Lemma 1.3, Theorem FM, Theorem R, Lemma Kern (iii), Prop B5, Theorems A and A′.
- **Theorem 4.3: WRONG as first stated.** For |G| ≈ ℓ the bound ϱ(T|G|) ≤ ρ_F·polylog fails by a factor
  L^{Θ(1)}. Repaired by the filter |G(P)| ≥ MV/32 on counted candidates; FM is unaffected. The L^{O(ε)}
  loss is repaired by S² = (log N)^{C₃}/(ρ_rρ_F).
- **Other fixes applied.**
  - Primitivity of w (Lemma 1.3).
  - LV's truncation is a flattening.
  - The sign case in FM Step 2(d).
  - An ℓ-independent grid with O(Λ) points.
  - The compactness argument via the full divided-minor vector.
  - The q-density argument (O(1) solutions of det(P, P′) = ±ℓq).
  - LD is weaker than RBDH′, not equivalent to it.
  - B2′, B3′, B5 and Z wording.
  - A′ constants: 0.9670 and 0.9142, with a closed form.
  - The explicit HECM, reducibility, g₁-column and repeated-test assumptions.
- **Literature (codex).** Lemma DS is prior art in substance (Harper, Compositio 2016, Smooth Numbers
  Result 3; Friedlander 1981). No earlier rigorous use of the 2D rational randomization was found.

### 7.4 Opus proof critic: adjudication (all accepted)
- **Verdicts.**
  - CORRECT: Lemma GB, Lemmas 1.2 and 1.3, Lemma 8.5′ (both parts; "LV's proof bounds W piece by
    piece"), Lemma Q, Lemma DS, Lemma Kern (i)–(iii) (300 random pairs checked independently), Prop 4.2,
    Prop Z, and the Theorem A arithmetic.
  - CORRECT WITH FIXES: Theorem FM, Theorem R, Theorem A′ rounding, and the parameter table.
  - Theorem 4.3: a major but repairable gap at small |G|, which was already repaired.
- **"Does LD + proved pieces give PC and Theorem A?"** Yes, after the repair. Two inputs were used
  silently: irreducibility in the kept field (now quantified in §1) and HECM (listed; A′ removes it).
  The critic found no L[1/2]-size loss.
- **Fixes applied.**
  - Irreducibility quantified, including the kept-field amplification by 1/μ and the degree drop.
  - B1′ and B5 reworded as barriers to naive applications, and the critic's route to LD added (§4.3).
  - Kern: the summary wording, the 4V² constant, and an effective wedge bound via the minor (bb′)^{d−2}Δ.
  - FM Step 2(d): 2(|F̂|/t + 1), and F = 0 excluded.
  - Theorem R: X = 2X_r, and the Λ_ℓ notation.
  - ψ well defined, and a_d ≠ 0.
  - The regime-(B) caveat for LV's Page step.
  - The table relabel for (a′_k) FM (outline only).
  - The Bernstein batch-smoothness remark.

### 7.3 Opus literature critic: adjudication (all accepted)
- Unconditional and GRH-conditional c < 1: not found. Wording now reads "to our knowledge".
- BLP 1993 / Murphy §2.3.2 credited for the g₁x − g₀ design.
- Hildebrand 1985 and divisor switching credited for Lemma DS.
- Harper's JLMS 2025 BDH asymptotics added to B1′.
- Claim 5.2 wording (affine point).
- "Nothing else is missing" qualified.
- B5's "every published" narrowed.
- LV journal reference, Kleinjung–Wesolowski description and Pollard–Strassen references added.
- The ECM "most smooth numbers" route noted as unchecked.

### 7.1 GPT-6.1-sol blind round (pcest_1): adjudication
- **Independent agreement.** Codex independently proposed the balanced rational randomization, Lemma GB
  (the generalized base expansion with A_d free), the elementary progression bound (its (1) = our Lemma
  DS), the index-q pair lattice (its (3) = our Lemma Q), and the proof of PC_rat (its (4) = our
  Theorem R). It also proposed the kernel transfer with cost η > max(2ν, s − (d−2)ν) (its (5) = our
  Lemma Kern and cost). It left the algebraic weighted bound (6) open. **Accepted; all consistent with
  §1–§4.**
- **Corrections taken from codex.**
  - "Never applicable for fixed k" for Pascadi Cor. 7.2 is false. The coefficient condition can hold
    for large fixed k, but the exponent 13/8 still falls short (§3, Remark 2).
  - The e^{O(u₂)} factor in de la Bretèche–Drappeau is harmless in upper bounds and fatal only against
    a relative saving.
  - Restricting the enlarged moduli to below |F|^{1/2} would need h > s + (6−d)/(4k). We do not need
    this.
- **Codex's alternatives.**
  - Two algebraic sides: (1 + 4/(d+1))/√2, worse (table).
  - Thinning: preserves the ratio, so it does not help.
  - Design V: unchanged.

## 8. Candidate cards

```text
ID: KD-2 bounded-k Dixon descent, design (a''_k): 2D-randomized rational side g = g1 x - g0.
Target: every N (after trial division to B).
Operation: as KD-1 (hypD), with (g0, g1) uniform coprime in [M/2, M]^2, f = gbe_g + g R, candidates
  P = (a, b), b prime in (V/2, V], |a| <= V, in the special-q lattice of l.
Proved: Lemma GB; Theorem FM (first moment, full proof); Theorem R (PC_rat, polylog loss, elementary
  divisor switching); Lemma Kern (H >= V^2 transfer); Prop 4.2; Theorem 4.3 (LD => D_k).
New lemma required: Hypothesis LD (restricted lift-dispersion BDH upper bound, precision L^{o(1)});
  LV6 (LV section 6 at height N^{1/(4k)}).
Cost: L[1/2, (1+theta)/sqrt 2] conditional; inf over k (1 + 3/(d+1))/sqrt 2 (d = 9: 0.919); with
  Pollard-Strassen instead of HECM: d = 9: 0.967.
Status: conditional on LD + LV6. E1 toy: PC_rat 1.01, PC 1.06-1.28, 2D ~ 1D.

ID: R-1 divisor switching for smooth pairs in 2D lattice regions (Lemma DS + Theorem R).
Status: proved. Scope: moduli <= y^{O(1)}, regions of width >= X y^{-O(1)}; fails in 1D (no area main term)
  and for large moduli (B2').

ID: B5 L^1 -> L^2 barrier for lift dispersion.
Status: proved (elementary). Any LD proof needs a large-values estimate for overloaded lift classes.
```
