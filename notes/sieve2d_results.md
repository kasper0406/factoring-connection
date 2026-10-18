# sieve2d: a sharper exceptional-set bound for Theorem D′ (two-variable fourth moment)

Agent `sieve2d`, 2026-10-18 (restarted run; the first run's machine was lost after the
pre-registration commit `0e47a2d`). Pre-registration: `notes/sieve2d_prereg.md` (frozen,
unchanged). Code: `src/sieve2d/{coverage,sampling,experiments,analyze,bound,explore}.py`.
Tests: `tests/test_sieve2d.py`, `tests/test_sieve2d_bound.py`. Raw data:
`experiments/raw/sieve2d_{exh,samp,explore}.jsonl`. Summary:
`experiments/analysis/sieve2d_summary.json`. Blind GPT-6.1-sol attempt: `codex/sieve2d_1.out`
(reached the same theorem independently, §7). Notation is that of
`notes/dcltheorem_results.md` (Theorem D′, authoritative).

## 0. Verdict (plain language)

- **Goal (a), an e^{−cλ} (Poisson-quality) bound, is NOT achieved.** I also show it cannot
  come from any argument that only uses the joint behaviour of a bounded number of band
  primes (§4). Any such proof needs genuinely new arithmetic input.
- **Goal (b): proved, a λ^{−2} bound.** The two-variable fourth moment, averaged over both
  p and q (Theorem 1), gives
  - the uncovered integer fraction of a box is ≤ 3/λ² + 1/λ³ + (explicit errors), and
    the errors are o(λ^{−2}) in the D′ regime as n → ∞;
  - consequently |E_b|/|S_n| ≤ C·n⁴/H⁴ (Theorem 2), instead of D′'s C·n²/H².
- **Effect on the factoring theorem (Theorem 3, proved).** The algorithm DCL-det_{n,ε} of D′,
  completely unchanged, fails on only O(N^{−2ε}) of S_n, not O(N^{−ε}), for ε ≤ 1/80.
  - Equivalently: N^{−ε} exceptions at cost N^{1/6+5ε/12}·polylog, instead of N^{1/6+5ε/6}.
  - Equivalently: cost N^{1/6+δ} buys exceptions N^{−12δ/5}, instead of N^{−6δ/5}.
  - In the polylog regime (H = n^{1+a}) exceptions drop from (log N)^{−2a} to
    (log N)^{−max(4a, 1+2a)}. This uses an analysis-only T = ⌊2^{n/5}⌋ (Corollary 3′).
    The Opus critic caught the original T.
- **Also proved, by a standard argument (Proposition 4).** Sieving p by the class 0 modulo small
  primes, in the one-variable large sieve, gives the baseline O(1/λ) = O(n/H²). This was
  found independently by GPT-6.1-sol (sieve2d_1, r16).
- **Barrier (goal (c)), in two parts.**
  - *Proved, citing Peled–Yadin–Yehudayoff 2011:* exact knowledge of all joint statistics of ≤ J
    band primes is consistent with an uncovered mass ≥ c_J·λ^{−⌊J/2⌋} (Proposition 5).
  - *Precise computation, not a theorem about all methods:* this error analysis of averaging
    over (p,q) controls those statistics for J ≤ 4 (proved), and J = 5 is sketched. At J = 6
    the per-lattice O(1) error *bounds* add up to ≈ X²/H², far above the needed X²λ^{−3}.
  - Exploratory data (§5.4) show that the true 6-prime discrepancies cancel. So this is a
    barrier for the proof technique, not for the truth.
  - So λ^{−2} is where *this* level-type analysis stops. It is **not** proved that the box
    itself determines nothing beyond ≤5-prime statistics.
  - The fundamental lemma needs a level exponent s ≳ κ ≈ |V|/2. The box offers s ≤ 6.
- **Experiments (pre-registered, run exactly as frozen).**
  - P1, P2, P3, P4, P6 **hold**. Uncovered counts match the independent-residue model to
    pooled 1.00001 over 202 cells; the decay is exponential, with fitted slope −0.9995 in λ_eff.
  - **P5 FAILED**, in 34 integer cells (0 prime cells). The integer-pair 4th moment is far
    off the model (ratios 0.007–1.27), while every prime cell lies in [0.988, 1.015].
  - The diagnosis (post hoc, §5.3) is a heavy-tailed rare atom in the integer model, not a
    failure of the main term.
  - Exploratory sixth moments match the model (§5.4). So the j = 6 obstruction is a proof
    barrier, not a visible effect in the data.

## 1. Setting

- X ≥ 1 is an integer and the box is 𝔅₀ = (0, X]² ∩ Z².
- The band 𝒫 is a nonempty set of primes in (Y, 2Y].
- V is a nonempty set of coprime (A,B) ∈ [1,H]². D′ uses V(H), which also requires A ≠ B;
  that is not needed here.
- s = Σ_{ℓ∈𝒫} 1/ℓ and λ = |V|·s. |𝒫| is the number of band primes.
- For x = (a,b) and v = (A,B), write F_v(x) = Ba − Ab = det(x, v).
  - z_ℓ(x) = #{v ∈ V : ℓ | F_v(x)}.
  - W(x) = #{ℓ ∈ 𝒫 : z_ℓ(x) > 0}.
  - x is *uncovered* iff W(x) = 0. For a prime pair this is exactly D′'s condition
    "no (v,ℓ) with ℓ | Bp − Aq".
- ℒ = ∪_{v∈V} R·v is the union of the lines through the weights, and 𝔅 = 𝔅₀ ∖ ℒ.
  - Points of ℒ are covered by every band prime.
  - No prime pair p ≠ q, with p, q > H, lies on ℒ: Bp = Aq forces p | A ≤ H.
- A *configuration* of order j is c = ((ℓ₁,v₁), …, (ℓ_j,v_j)) with **distinct** ℓ_i ∈ 𝒫
  and v_i ∈ V. Its lattice is Λ_c = {x ∈ Z² : ℓ_i | F_{v_i}(x) ∀i}, and d_c = ∏ℓ_i.
- The model: Z = Σ_ℓ Z_ℓ with independent Z_ℓ ~ Bernoulli(|V|/ℓ). Then
  E(Z)_j = Σ_c 1/d_c, where (w)_j = w(w−1)⋯(w−j+1).

## 2. Theorem 1 (two-variable fourth moment) — **proved**

**Theorem 1.** Assume H ≥ 2, Y ≥ 16H², HX < Y⁴ and λ ≥ 1. Then

  #{x ∈ 𝔅 : W(x) = 0} ≤ X²(3/λ² + 1/λ³) + 2HX + λ^{−4}·Σ_{j=1}^{4} κ_j·Err_j,

with κ₁ = 15λ³, κ₂ = 25λ², κ₃ = 10λ, κ₄ = 1. Here Err_j = A_j + B_j + C_j + D_j + E_j + F_j,
where (⌈·⌉, ⌊·⌋ applied to j/2):

| term | value |
|---|---|
| A_j | (\|𝒫\|·\|V\|)^j |
| B_j | 28.3·X·(2Y)^{j/2}·(2\|V\|)^j |
| C_j | 2HX·Σ_{r=⌈j/2⌉}^{j} C(j,r) λ^r \|𝒫\|^{j−r} |
| D_j | X²\|V\|s⁴ if j = 4, else 0 |
| E_j | 6HX·Σ_{r∈I_j} C(j,r)λ^r\|𝒫\|^{j−r} + 2HX\|V\|·Σ_{r∈I_j} C(j,r) s^{j−r}(\|𝒫\|\|V\|)^r, with I_j = [max(0,j−3), ⌊j/2⌋] |
| F_j | (3\|V\|)^j·X²·s/Y |

(`src/sieve2d/bound.py::err_terms`; `tests/test_sieve2d_bound.py` checks every displayed
inequality exhaustively on small boxes, and every per-configuration claim of §2.2–2.3 on
all configurations of order ≤ 4 for small bands.)

**Corollary 1 (asymptotic form; proved).** Under the same hypotheses, with absolute O(·),

  #{x ∈ 𝔅 : W = 0}/X² ≤ 3/λ² + 1/λ³ + O( |V|^{−3} + HY³/(|V|³X) + Y²/(Xs⁴) + 1/(Ys³)
                                          + Y⁴/X² + (H/X)(1 + Y²/|V|³) ).

The terms come from rows D, E (first part), B, F, A and E (second part, j = 1), using
|𝒫| ≤ 2Ys and λ ≤ |𝒫|. The latter holds because λ = |V|s ≤ (Y/16)s < |𝒫|/16.
- Rows C and E (second part, j ≥ 2) give terms like H·Y²/(X|V|²) and (H/X)(2Y)^r|V|^{1+r−j},
  for example 4HY²/(X|V|) at j = 4, r = 2, and 2HY/X at j = 2, r = 1. These are absorbed by
  max(HY³/(|V|³X), Y²/(Xs⁴)):
  - if |V|² ≤ Y, the first term covers them;
  - otherwise |V| > 4H and s ≤ 1, and the second does.
- (Absorption line added after the Opus critic.)

### 2.1 Elementary facts

- **(F1)** If ℓ ∈ 𝒫 does not divide both a and b, then z_ℓ(x) ≤ 1.
  - Suppose ℓ | F_v(x) and ℓ | F_{v′}(x) with v ≠ v′. Then ℓ | (AB′−A′B)·a and ℓ | (AB′−A′B)·b.
  - Here 0 < |AB′ − A′B| < H² < ℓ: distinct primitive positive vectors are not parallel.
  - So ℓ | a and ℓ | b.
- **(F2)** If x ∈ 𝔅, every F_v(x) is a nonzero integer with |F_v(x)| < HX < Y⁴.
  - Indeed −AX < Ba − Ab < BX.
  - So F_v(x) has at most 3 prime factors in (Y,2Y], and Σ_ℓ z_ℓ(x) ≤ 3|V|.
- **(F3)** Λ_c has index d_c.
  - Since gcd(A_i,B_i) = 1, x ↦ F_{v_i}(x) mod ℓ_i is onto Z/ℓ_i, with kernel of index ℓ_i.
  - The ℓ_i are distinct, so the CRT applies.
- **(F4)** For v ∈ V, put r_v(c) = #{i : v_i ≠ v} and t_v(c) = ∏_{i: v_i≠v} ℓ_i. Then
  Λ_c ∩ Rv = Z·t_v v.
  - σv ∈ Λ_c (σ ∈ Z) iff ℓ_i | σ(B_iA − A_iB) for all i.
  - That holds automatically when v_i = v. Otherwise it holds iff ℓ_i | σ, since
    0 < |B_iA − A_iB| < ℓ_i.
- **(F5)** The shortest vector u_c of Λ_c satisfies |u_c| ≤ (2/√3)^{1/2}·√d_c ≤ 1.0746√d_c
  (Hermite's constant in dimension 2).
  - u_c is primitive in Λ_c, so Λ_c has a basis (u_c, w) with |det(u_c, w)| = d_c.

### 2.2 Lemma L (lattice points in a box) — proved

For a lattice Λ ⊂ R² of covolume d with shortest vector length ρ, and any convex K of
diameter D,

  |#(Λ ∩ K) − area(K)/d| ≤ 1 + (4/√3)·D/ρ.

For K = [0,X]² (or (0,X]²), D = √2X and the bound is 1 + 3.266X/ρ.

*Proof.*
- Take a Lagrange–Gauss reduced basis b₁, b₂, with |b₁| = ρ ≤ |b₂| and angle in [60°, 120°].
  Then d ≥ (√3/2)|b₁||b₂|.
- In coordinates x = σb₁ + τb₂, Λ becomes Z². K becomes a convex set K′ of area area(K)/d.
- σ(x) = det(x,b₂)/det(b₁,b₂) varies over K by at most D|b₂|/d ≤ (2/√3)D/ρ. Likewise τ
  varies by at most D|b₁|/d ≤ (2/√3)D/ρ.
- Count Z² ∩ K′ column by column.
  - Each column σ ∈ Z meets K′ in an interval of length h(σ), containing between h(σ) − 1
    and h(σ) + 1 integers.
  - There are at most w_σ + 1 nonempty columns, where w_σ is the width of K′ in σ.
  - h is concave on its support, hence unimodal. So |Σ_σ h(σ) − area(K′)| ≤ max h ≤ w_τ.
- Hence |# − area(K′)| ≤ w_σ + w_τ + 1. ∎

(`test_lattice_lemma_random`: 300 random lattices.)

### 2.3 Lemma C (one configuration) — proved

Let c have order j ≤ 4, d = d_c and u = u_c. Exactly one of the following holds.

**Case (i): u ∉ ℒ.** Then

  |#(Λ_c ∩ 𝔅) − X²/d| ≤ 1 + 3.266X/|u| + Σ_{v′: r_{v′}(c) ≥ ⌈j/2⌉} X/(t_{v′}(c)·max(A′,B′)).

**Case (ii): u = ±t_v(c)·v for one v ∈ V.** Write r = r_v(c) and g = d/t_v, the product of
the j − r primes labelled v. Then r ≤ ⌊j/2⌋, and:

- if j − r ≥ 4, then Λ_c ∩ 𝔅 = ∅;
- if j − r ≤ 3, then

  |#(Λ_c ∩ 𝔅) − X²/d| ≤ 2√2X/(t_v|v|) + 2HX/g + Σ_{v′≠v} X/(t_{v′}·max(A′,B′)).

In both cases every v′ in the last sum (v′ not parallel to u) has r_{v′}(c) ≥ ⌈j/2⌉.

*Proof.*

- **(a) Directions not parallel to u.** Let v′ ∈ V with t_{v′}v′ ∦ u.
  - t_{v′}v′ ∈ Λ_c is independent of u. So |t_{v′}v′|·|u| ≥ d, i.e.
    |t_{v′}v′| ≥ √d/1.0746 ≥ 0.93·Y^{j/2}.
  - But |t_{v′}v′| ≤ (2Y)^{r_{v′}}·√2H.
  - If r_{v′} ≤ ⌈j/2⌉ − 1, then 2^{r}·1.53H ≥ 0.93·Y^{j/2−r} with j/2 − r ≥ 1/2. Checking
    j = 1, …, 4 against Y ≥ 16H² rules this out: the tightest case is j = 3, r = 1, which
    needs √Y ≤ 3.3H, while √Y ≥ 4H.
  - Points of Λ_c ∩ 𝔅₀ on Rv′ are σt_{v′}v′ with 1 ≤ σ ≤ X/(t_{v′}max(A′,B′)). This gives
    the last sums. (In case (i) no v′ is parallel to u.)
- **(b) Case (i).** Use Lemma L for #(Λ_c ∩ 𝔅₀), then subtract the points on ℒ, which (a)
  bounds.
- **(c) Case (ii), r ≤ ⌊j/2⌋.**
  - |u| = t_v|v| ≥ Y^r and |u| ≤ 1.0746(2Y)^{j/2}.
  - So r > j/2 would force Y^{1/2} ≤ 1.0746·2^{j/2} ≤ 4.3, contradicting Y ≥ 64.
- **(d) Case (ii), line structure.**
  - With a basis (u, w), det(u, x) ∈ dZ for x ∈ Λ_c. So det(v, x) = det(u,x)/(±t_v) ∈ gZ.
  - Each line {det(v,·) = mg} meets Λ_c in a coset of Zu, which is nonempty since x = mw lies
    on it.
  - On 𝔅₀, |det(v,x)| = |F_v(x)| < HX.
  - If j − r ≥ 4, then g > Y⁴ > HX, so only m = 0 occurs. That is the line Rv ⊂ ℒ, removed.
  - If j − r ≤ 3, let h(σ) be the length of the chord {y ∈ [0,X]² : det(v,y) = σ}.
    - h is concave on its support, max h ≤ √2X, and ∫h dσ = |v|X² (coarea; |∇ det(v,·)| = |v|).
    - Unimodality gives |g·Σ_{m∈Z} h(mg) − ∫h| ≤ g·max h.
    - Each nonempty line carries h(mg)/|u| ± 1 lattice points, and there are at most
      2HX/g values m ≠ 0 with |mg| < HX.
    - Dropping m = 0 (the removed line Rv) costs ≤ h(0)/|u|.
    - Since g|u| = d|v|, we get |#(Λ_c ∩ 𝔅₀ ∖ Rv) − X²/d| ≤ (max h + h(0))/|u| + 2HX/g.
  - Subtracting the other ℒ-lines as in (a) completes the proof. ∎

### 2.4 Summation over configurations — proved

Let Σ_c run over the configurations of order j. Then

  Σ_c |#(Λ_c ∩ 𝔅) − X²/d_c| ≤ A_j + B_j + C_j + D_j + E_j.

- **A_j:** the "+1" of case (i), summed over at most (|𝒫||V|)^j configurations.
- **B_j:** in case (i), c is counted at its own shortest vector u = u_c ∉ ℒ, with
  |u| ≤ R_j := 1.0746(2Y)^{j/2}.
  - For such u, each F_v(u) is nonzero with |F_v(u)| ≤ √2H·R_j ≤ 6.1HY² < Y³ (Y > 6.1H).
    So F_v(u) has ≤ 2 band prime factors.
  - Hence u satisfies ≤ 2|V| events (ℓ, v), and lies in Λ_c for ≤ (2|V|)^j configurations.
  - With Σ_{0<|u|≤R} 1/|u| ≤ 8R (there are 8m vectors with |u|_∞ = m, each with |u| ≥ m):
    Σ_{case (i)} 3.266X/|u_c| ≤ 3.266X·8R_j(2|V|)^j ≤ B_j.
- **C_j:** swap the sum over c with the sum over v′.
  - For fixed v′ and r, Σ_{c: r_{v′}(c)=r} 1/t_{v′}(c) ≤ C(j,r)·λ^r·|𝒫|^{j−r}. Choose the r
    positions; each of them contributes Σ_{ℓ,v} 1/ℓ = λ, and each other position contributes
    a factor ≤ |𝒫|. Dropping distinctness only enlarges the sum.
  - Σ_{v′} 1/max(A′,B′) ≤ Σ_{m≤H}(2m−1)/m < 2H.
- **D_j:** case (ii) with j − r ≥ 4 occurs only for j = 4 and r = 0 (all four labels equal).
  The error is X²/d, and Σ_v Σ X²/(ℓ₁⋯ℓ₄) ≤ |V|X²s⁴.
- **E_j:** case (ii) with j − r ≤ 3, so r ∈ I_j.
  - Σ_{c: r_v=r} 1/t_v ≤ C(j,r)λ^r|𝒫|^{j−r} and Σ_{c: r_v=r} 1/g ≤ C(j,r)s^{j−r}(|𝒫||V|)^r.
  - Σ_v 2√2X/|v| ≤ 2√2X·2H ≤ 6HX.

### 2.5 Factorial moments and assembly — proved

- Put F_j(x) := Σ_c 1[x ∈ Λ_c] = Σ_{ℓ₁,…,ℓ_j distinct} ∏ z_{ℓ_i}(x).
  - By (F1), F_j(x) = (W(x))_j unless x ∈ 𝒟 := {x : ∃ℓ ∈ 𝒫, ℓ | a, ℓ | b}.
  - On 𝔅 both lie in [0, (Σ_ℓ z_ℓ)^j] ⊂ [0, (3|V|)^j] by (F2).
  - #(𝒟 ∩ 𝔅₀) ≤ X²Σ_ℓ ℓ^{−2} ≤ X²s/Y.
  - Hence |Σ_{𝔅}(W)_j − Σ_c #(Λ_c∩𝔅)| ≤ F_j (the table entry), and
    |Σ_𝔅 (W)_j − X²E(Z)_j| ≤ Err_j.
- For j = 0: |#𝔅 − X²| = #(ℒ ∩ 𝔅₀) ≤ Σ_v X/max(A,B) < 2HX.
- Expand in falling factorials, using Stirling numbers:
  (w−λ)⁴ = (w)₄ + (6−4λ)(w)₃ + (7−12λ+6λ²)(w)₂ + (1−4λ+6λ²−4λ³)(w)₁ + λ⁴.
  - For λ ≥ 1 the coefficients are bounded by κ_j.
  - The same identity holds for Z, and E(Z−λ)⁴ = Σ_ℓ μ₄(Z_ℓ) + 3Σ_{ℓ≠ℓ′}σ_ℓ²σ_ℓ′² ≤ λ + 3λ².
    A Bernoulli(π) variable has μ₄ = π(1−π)(1−3π+3π²) ≤ π and σ² ≤ π.
- Every uncovered x ∈ 𝔅 contributes (0−λ)⁴ = λ⁴. So

  λ⁴·#{W=0} ≤ Σ_𝔅 (W−λ)⁴ ≤ X²(3λ²+λ) + λ⁴·2HX + Σ_{j=1}^{4} κ_j Err_j. ∎

**Remarks.**
- (1) The repeated-direction (same-slope) degeneracy is essential. Exact equidistribution is
  false for 4 equal labels: the count is 0, not X²/d (case (ii), j − r ≥ 4). Removing ℒ, as
  the orchestrator/r16 warned, is what makes it harmless (row D).
- (2) The weight-free statistic Q = mean of (1 − M/λ)⁴ of the pre-registration (P5) is this
  fourth moment.

## 3. Consequences for Theorem D′ — proved

**Theorem 2 (coverage exceptional set).**
- Setting: the D′ design with any H satisfying 2 ≤ H ≤ N^{1/15}/C₀; n ≥ n₀; Y = ⌈2^{n/6}H^{−1/3}⌉;
  four blocks 𝒫_k ⊂ (2^kY, 2^{k+1}Y] of primes ≡ 3 (mod 4).
- E_b is the set of pairs of S_n uncovered by their clean block k(q).
- Then

  |E_b|/|S_n| ≤ C·n²·(λ₀^{−2} + H^{−6}) + C·n⁶H^{1/3}N^{−1/6} ≤ C′·n⁴/H⁴ + C′·n⁶H^{1/3}N^{−1/6}.

*Proof.*
- E_b ⊂ ∪_{k=0}^{3} {(p,q) ∈ S_n : W_k(p,q) = 0}.
- Prime pairs of S_n lie in 𝔅 with X = ⌈√2·2^{n/2}⌉. They avoid ℒ (§1).
- Apply Theorem 1 to each block, with Y_k = 2^kY ≥ Y.
  - Hypotheses: Y ≥ 16H² ⟸ 2^{n/6} ≥ 16H^{7/3}, and HX < Y⁴ ⟸ H^{7/3} < 2^{n/6}/√2. Both
    hold for H ≤ N^{1/15}/C₀.
  - λ_k ≥ c₁H²/n by D′ Lemma 9. Take H ≥ C√n so that λ_k ≥ 1; otherwise the claim is trivial.
- Use Corollary 1 at this scale:
  - Y_k³ ≤ 512Y³ and Y³H ≤ 2P, so HY_k³/(|V|³X) = O(|V|^{−3}) = O(H^{−6});
  - Y_k²/(Xs⁴) = O(n⁴P^{−1/3}H^{−2/3});
  - 1/(Y_k s³) = O(n³P^{−1/3}H^{1/3});
  - the remaining terms are smaller.
- So #{uncovered in block k} ≤ X²[3λ_k^{−2} + λ_k^{−3} + O(H^{−6} + n⁴P^{−1/3}H^{1/3})].
- Divide by |S_n| ≥ c₀P²/n². ∎

**The explicit constants are poor.** `bound.scale_terms(n, H, k)` evaluates every
normalised term at the D′ scale, and the error terms beat the main term only for large n.
The dominant pieces are rows B₄ (shortest-vector sum) and F₄ (crude bound on 𝒟).

| n | H | block | λ | 3/λ² + 1/λ³ | Σ errors (normalised) | D′ bound n²/H² (no constant) |
|---|---|---|---|---|---|---|
| 96 | 32 | 0 | 21.7 | 6.5e−3 | 1.4e3 | 9.0 |
| 200 | 2^{10} | 0 | 1.1e4 | 2.7e−8 | 3.2e−2 | 3.8e−2 |
| 300 | 2^{15} | 0 | 7.3e6 | 5.7e−14 | 1.8e−6 | 8.4e−5 |
| 400 | 2^{20} | 0 | 5.6e9 | 9.7e−20 | 1.3e−10 | 1.5e−7 |

So Theorem 2 is an asymptotic improvement. At cryptographic-scale n the second term
(n⁶H^{1/3}N^{−1/6}-type) dominates the proved bound. Even then the proved bound already beats
D′'s n²/H² at n ≥ 300 in this table. Sharpening rows B and F, for example with smooth
weights and an average over the shortest vectors, is routine but not done.

**Proposition 4 (one-variable large sieve with a prime presieve; proved, standard).**
|E_b|/|S_n| ≤ C/(1 + λ_min) = O(n/H²).

*Proof.*
- Fix a prime q and a block k. Let U(q) be the set of *prime* p ∈ (P/2, P] uncovered by block k.
- Put Q = ⌊√P⌋ and Z = ⌊Q/2^{k+1}Y⌋. Then 1 ≤ Z < Y and ln Z ≥ (1/6 − o(1))·ln P.
- U(q) avoids the |V| classes AqB^{−1} mod each ℓ ∈ 𝒫_k, and the class 0 mod each prime r ≤ Z,
  since p > P/2 > Z.
- The arithmetic large sieve (Montgomery 1968; Montgomery–Vaughan 1973) gives
  |U(q)| ≤ (P/2 + 1 + Q²)/J.
- In J, keep d = m and d = mℓ, with m ≤ Z squarefree and ℓ ∈ 𝒫_k.
  - These d are ≤ Q and pairwise distinct, because the primes of m are < Y < ℓ.
  - So J ≥ (Σ_{m≤Z} μ²(m)/φ(m))·(1 + Σ_ℓ |V|/(ℓ−|V|)) ≥ (1+λ_k)·ln Z.
- Summing over the ≤ 3P/n primes q gives |E_b| ≤ C P²/(n²(1+λ)). ∎

This is the r16 / sieve2d_1 refinement. I verified it line by line. It removes D′'s n-loss but
not the power of λ.

**Theorem 3 (cost/exception trade-off; proved).** Let ε′ ∈ (0, 1/80] and n ≥ n₀(ε′). The
algorithm DCL-det_{n,ε′} of D′ §2, unchanged, has:
- H = ⌈n·2^{ε′n/2}⌉;
- cost N^{1/6+5ε′/6}(log N)^{O(1)};
- at most C·2^{−2ε′n}·|S_n| FAIL pairs.

Restated with ε = 2ε′: **N^{−ε} exceptions at cost N^{1/6+5ε/12}(log N)^{O(1)}, for
ε ≤ 1/40.** Equivalently, cost N^{1/6+δ+o(1)} with an exceptional fraction O(N^{−12δ/5}),
for δ ≤ 1/96.

*Proof.* Follow D′ §4(iii), with T = ⌊2^{(1/4−ε′)n}/n⌋. T enters only the analysis.
- **E_a.** By D′ Lemma 8, |E_a| ≤ (T(T+1)/(n−1))·(P/2+1). Dividing by |S_n| ≥ c₀P²/n² gives
  O(n·T²/P) = O(2^{−2ε′n}/n).
- **E_b.** By Theorem 2, n⁴/H⁴ ≤ 2^{−2ε′n}, and n⁶H^{1/3}N^{−1/6} is far smaller.
- **E_c.** D′ bounds it by C·n·(H^{14/3}P^{1/3}/T + H^{10/3}/P^{1/3}).
  - H^{14/3} ≤ (2n)^{14/3}2^{7ε′n/3} and P^{1/3}/T ≤ 2n·2^{(ε′−1/12)n}. So
    H^{14/3}P^{1/3}/T ≤ (2n)^{6}·2^{(10ε′/3 − 1/12)n}.
  - We need 10ε′/3 − 1/12 < −2ε′, i.e. 16ε′/3 < 1/12, i.e. ε′ < 1/64. For ε′ ≤ 1/80 the
    margin is 1/12 − 16/240 = 1/60. So this term is O(2^{−2ε′n − n/60}·n^{O(1)}).
  - The second term is O(2^{(5ε′/3−1/6)n}).
- **Scale separations.** D′ needs Y > H² and Y⁴ > √2P. Theorem 1 needs Y ≥ 16H² and
  HX < Y⁴. All hold for ε′ ≤ 1/80.
- Cost and correctness are D′(i)–(ii) verbatim. ∎

**Corollary 3′ (polylog regime).** Take H = ⌈n^{1+a}⌉ and the analysis-only T = ⌊2^{n/5}⌋.
- Cost: N^{1/6}(log N)^{O(1+a)}.
- E_a = O(nT²/P) = O(n·2^{−n/10}), and E_c = n^{O(1)}·2^{−n/30}.
- By Proposition 4 and Theorem 2, E_b = O(min(n^{−1−2a}, n^{−4a})).
- So the exceptional fraction is O(n^{−max(1+2a, 4a)}). D″ had O(n^{−2a}).
- *Correction (Opus critic):* the first draft kept D″'s T = 2^{n/4}n^{−1−a}. Its E_a term
  O(n^{−1−2a}) is not o(n^{−4a}) for a ≥ 1/2.

## 4. Barrier analysis

### 4.1 An information barrier (Proposition 5; proved, using PYY 2011)

**Proposition 5.**
- Let 𝒫 be any band and π_ℓ ∈ (0, 1/2] with max π_ℓ ≤ 2 min π_ℓ and Σπ_ℓ = λ.
- For every fixed J ≥ 2 there are c_J > 0 and λ_J such that for λ ≥ λ_J there exist
  J-wise independent Bernoulli variables (X_ℓ)_{ℓ∈𝒫} with P(X_ℓ = 1) = π_ℓ and

  P(X_ℓ = 0 for all ℓ) ≥ c_J·λ^{−⌊J/2⌋}.

*Proof.*
- Let π* = max π_ℓ, n′ = |𝒫| and μ = n′π* ≤ 2λ.
- Peled–Yadin–Yehudayoff, "The maximal probability that k-wise independent bits are all 1",
  Random Struct. Alg. 38 (2011) 502–525 (arXiv 0801.0059), handle even k:
  - Let M(n′,k,p) be the maximum of P(all ones) over k-wise independent bits with marginal p.
  - Corollary 1.2(1): for fixed even k and N := n′p(1−p) − 1 > N₀, M ≥ c^k·M̃ with
    M̃ = p^{n′}/P(Bin(n′,1−p) ≤ k/2). The PYY upper bound is M ≤ M̃.
  - Corollary 1.3(1) gives directly M ≥ (c/√k)·(pk/(2e(1−p)n′))^{k/2}.
  - *Correction (Opus critic):* the first draft misquoted Theorem 1.1 as (c/k)·M̃. Its
    actual factor exp(−c₂k·V(N/k)) is not constant, but for fixed k Corollary 1.2(1) gives
    the constant c^k.
  - Their (1.1), from Benjamini–Gurel-Gurevich–Peled, handles odd k: M(n′,k+1,p) = p·M(n′−1,k,p).
- Take p = 1 − π*, so N ≈ μ → ∞. Then P(Bin ≤ k/2) ≤ p^{n′}Σ_{i≤k/2}(2μ)^i/i!, so
  M ≥ c_k μ^{−k/2}.
- Complementing the bits gives J-wise independent Y_ℓ with P(Y_ℓ=1) = π* and
  P(all Y = 0) ≥ c_J μ^{−⌊J/2⌋}.
- Thin independently: X_ℓ = Y_ℓ·η_ℓ with η_ℓ ~ Bernoulli(π_ℓ/π*). Any ≤ J of the X_ℓ are
  independent with the right marginals, and P(all X = 0) ≥ P(all Y = 0). ∎

For prime pairs, attach to each X_ℓ = 1 an independent uniform label v ∈ V. All labelled
(ℓ,v)-event statistics of order ≤ J then match the independent-residue model exactly, while
the uncovered mass is ≥ c_J λ^{−⌊J/2⌋}.

So **no argument whose only input is the joint distribution of ≤ J band-prime events (moments,
Selberg/Bonferroni majorants, any polynomial of degree ≤ J in the event indicators) can prove
o(λ^{−⌊J/2⌋})**, even with that input known exactly.

The matching upper side is Selberg's. The degree-2m square majorant gives
[Σ_{i≤m} e_i(π_ℓ/(1−π_ℓ))]^{−1}. In the Poisson limit max π_ℓ → 0 this is
(Σ_{i≤m} λ^i/i!)^{−1} ~ m!/λ^m, the Christoffel number of the Poisson law, i.e. PYY's M̃.

*Heuristic:* e^{−cλ} by this route needs J ≍ λ, i.e. lattices of index Y^{cλ}, astronomically
beyond the box volume X² ≈ Y⁶H². Proposition 5 is stated for fixed J; the growing-J version
would need PYY's growing-degree bounds (their Corollary 1.2(2)–(3)) and is not checked here.

### 4.2 Which J the box gives (precise computation; not a theorem about all methods)

- Theorem 1's machinery works for any j. The normalised error of order j, divided by X²λ^j, is
  dominated by:
  - B_j: ≈ (2Y)^{j/2}2^j/(Xs^j);
  - A_j: (|𝒫|/s)^j/X² ≈ (2Y)^j/X².
- **j = 4:** proved (Theorem 1). The normalised error is O(N^{−c} + |V|^{−3}).
- **j = 5:** *sketch, not written out.* Rows A and B are N^{−c}. Configurations with four
  repeated labels add a further normalised O(|V|^{−3}) (pointed out by GPT-6.1-sol
  sieve2d_2). So the accuracy is not uniformly N^{−c} in the polylog regime.
  - By Proposition 5 (J = 5 gives ⌊5/2⌋ = 2), the fifth moment cannot improve λ^{−2}.
  - The same holds for odd-degree majorants: tiny masses far out absorb odd moments.
- **j = 6:** Y⁶ ≈ X²/(2H²), so a typical 6-prime lattice has ≈ 2H² points. The two error
  sources are then:
  - A₆/X² ≈ 64Y⁶/X² ≈ 32/H². Each lattice count is an integer, so an O(1) error per lattice
    cannot be removed lattice by lattice. There are ≈ (λY)⁶·(const) lattices.
  - B₆/X² ≈ 10⁴/(H s⁶) ≈ n⁶/H.
  - Both must be ≪ λ^{−3} ≈ n³/H⁶ for a λ^{−3} bound. They are larger by factors ≥ H⁴/n³.
- So the displayed absolute-discrepancy estimates do not establish a sixth-moment
  improvement. Large upper bounds on errors do not prove that cancellation is impossible.
- A sixth-moment proof would need the **signed** discrepancies of ≈ (λY)⁶ lattices, each of
  size Θ(1), to cancel to o(n³/H⁴) per lattice on average.
- **Where those discrepancies live (heuristic, not proved).**
  - With smooth weights, a generic 6-prime lattice has negligible error. Its dual vectors
    have length ≈ H/P.
  - All the error sits on *elongated* lattices, those with a vector u of length ≤ d/X ≈ Y³/H.
    These are a ≈ 1/H² fraction, and each has count 0 or ≈ 2H².
  - u ∈ Λ_c says exactly that the events of c hold at the point u. So counting elongated
    lattices is the coverage problem again, at the smaller scale X/H², with a worse
    level-to-volume ratio.
  - The signs depend on where one lattice line falls relative to the box. Their cancellation
    is an equidistribution statement of the same type at the next scale.
  - This is the regress that level methods cannot break.
- GPT-6.1-sol (sieve2d_1, its (24)) gets the same j = 6 failure with the same factor.

**Fake distribution realised.** Proposition 5 gives a distribution with exact ≤5-prime
statistics and uncovered mass ≍ λ^{−2}.
- So λ^{−2} is sharp for every method whose input is ≤5-prime statistics. That covers
  Theorem 1's method, which obtains them to accuracy N^{−c}.
- It is **not** proved that the box determines nothing beyond ≤5-prime statistics. The
  failure at j = 6 is a failure of this error analysis, and §5.4 suggests the truth is much
  better.

### 4.3 Standard sieve tools, quantified

- **Fundamental lemma / β-sieve** (**heuristic** quantification: no uniform growing-dimension
  fundamental-lemma theorem is invoked).
  - The sieve dimension is κ ≈ |V|/2 ≈ 0.3H²: ω(ℓ) = |V| on a band of log-length log 2, and
    half of all primes are ≡ 3 (mod 4).
  - An error factor e^{−s} needs s = log D/log z ≳ κ, which is ≳ λ·log Y.
  - The box offers s ≤ log(X²)/log(2Y) ≈ 6 (1D: ≈ 3).
  - The upper-bound sieve below the sieving limit gives at best Selberg's 1/G(√D). With the
    *effective* level, D such that (|𝒫||V|)^{ω(d)} × O(1) ≤ X²λ^{−ω(d)/2}, that is ω(d) ≤ 5,
    this heuristically gives λ^{−2} again (**heuristic**, not proved).
- **Janson / Suen.** They need a product space or a dependency graph. Here all events depend on
  one point (p,q), and Proposition 5 shows that low-order correlations alone cannot imply the
  Janson conclusion.
- **BV / BDH.** Theorem 1 never needs primes in progressions; integers suffice. The obstacle at
  j = 6 is the lattice discrepancy itself, not the prime restriction. Moduli d ≈ Y⁶ ≈ P²/H² are
  in any case far beyond BV's √P range in one variable.

### 4.4 Primes versus integers

- **Experimentally they agree** (P4 ratio 1.00004).
- **Theorem 1 loses n²** by bounding prime pairs by all integer points. A Selberg presieve
  inside the fourth moment should give C/λ² for S_n (not done; **conjectured**). It changes
  only polylog factors and none of the exponents above.
- **Any e^{−cλ} proof must use a restriction that excludes the following dilates.**
  Primality, or even coprimality of a and b, would do. (Wording narrowed after GPT-6.1-sol
  sieve2d_2.)
  - The dilates x = t(H, H+1) with gcd(t, ∏𝒫) = 1 are uncovered integer points. Indeed
    F_v(x) = t(BH − A(H+1)) with 0 < |BH − A(H+1)| ≤ H² < Y.
  - There are ≫ P/H of them in the balanced region.
  - So an integer analogue of e^{−cλ} is false once λ ≫ log P (noted by sieve2d_1; checked).
  - Prime pairs are never such dilates.

## 5. Experiments (pre-registered; notes/sieve2d_prereg.md, commit 0e47a2d)

- All runs were exactly as frozen: same configs and seeds. `config` re-wrote an identical
  `experiments/configs/sieve2d.json` (git diff empty).
- E1: 620 cells, n = 20–30, prime and integer. E2/E3: 228 cells, 10⁶ pairs each.
- The runs took 9 min on 10 cores.

### 5.1 E2, prime pairs, block 0 (obs = uncovered of 10⁶; E = model expectation)

| n | H | Y | λ | λ_eff | obs | E | var/model | m4/model |
|---|---|---|---|---|---|---|---|---|
| 48 | 5 | 256 | 1.04 | 1.07 | 342283 | 342497 | 0.998 | 0.995 |
| 48 | 8 | 128 | 2.99 | 3.45 | 32033 | 31831 | 1.001 | 1.004 |
| 48 | 10 | 128 | 4.42 | 5.51 | 4032 | 4040 | 1.003 | 1.007 |
| 64 | 10 | 813 | 3.20 | 3.30 | 36886 | 37010 | 0.998 | 0.997 |
| 64 | 14 | 813 | 6.51 | 6.91 | 1063 | 999.6 | 1.000 | 0.999 |
| 64 | 16 | 813 | 8.16 | 8.80 | 166 | 150.4 | 1.003 | 1.004 |
| 64 | 20 | 813 | 13.1 | 14.9 | 0 | 0.34 | 1.001 | 1.002 |
| 80 | 12 | 5161 | 3.50 | 3.52 | 29941 | 29657 | 1.000 | 1.001 |
| 80 | 16 | 5161 | 6.14 | 6.21 | 2027 | 2017 | 0.999 | 0.998 |
| 80 | 20 | 5161 | 9.87 | 10.05 | 47 | 43.4 | 1.001 | 0.998 |
| 96 | 12 | 32768 | 2.91 | 2.91 | 54348 | 54369 | 1.004 | 1.015 |
| 96 | 16 | 32768 | 5.11 | 5.12 | 5938 | 6000 | 1.000 | 1.006 |
| 96 | 20 | 32768 | 8.21 | 8.23 | 236 | 265.7 | 0.999 | 0.998 |
| 96 | 24 | 32768 | 11.57 | 11.62 | 5 | 9.0 | 1.000 | 0.998 |

The full tables (all 228 E2 cells, all 620 E1 cells) are `table_E1`/`table_E2` in the
summary JSON.

### 5.2 Frozen predictions

| pred. | rule (frozen) | observed | verdict |
|---|---|---|---|
| P1 | E2 cells with E ≥ 30: obs/E ∈ [0.85, 1.15]; pooled ∈ [0.95, 1.05] | 202 cells, range [0.877, 1.104], pooled 1.00001 (prime 0.99999, int 1.00003) | **holds** |
| P1 (E1, descriptive) | obs/E ∈ [0.67, 1.5] | 6/558 cells outside (0.44–0.665, all n ≤ 30 prime thm cells, E = 32–297); pooled 0.996 | fails (descriptive only) |
| P2 | block 0, E < 5 pooled: Σobs ≤ ΣE + 3√ΣE + 3; every λ ≥ 6 cell ≤ 0.1·3/λ² | 1 ≤ 9.50 (ΣE = 2.13, 14 cells); 42/42 cells below 0.1·3/λ² | **holds** |
| P3 | weighted LSQ slope of log(obs/pairs) vs λ_eff ∈ [−1.15, −0.85]; vs λ ∈ [−1.3, −0.85] | −0.9995 (λ_eff), −1.023 (λ); 43 cells | **holds** |
| P4 | pooled int/prime ∈ [0.9, 1.1] | 1.00004 (101 cell pairs) | **holds** |
| P5 | every λ ≥ 1 cell: var ratio ∈ [0.9, 1.1], m4 ratio ∈ [0.8, 1.2]; int Q ratio ∈ [0.8, 1.2] | var 0.972–1.034 (all pass); m4: **34 integer cells fail** (ratios 0.007–1.27), all 68 prime cells pass (0.988–1.015) | **FAILS** |
| P6 | n = 30, block 0, f ∈ [0.01, 0.5]: dispersion D ∈ [0.5, 2] | D = 0.56, 0.61, 0.68, 0.79 (H = 4–7, 1932 q each) | **holds** |

**Decisions (frozen rules).**
- P1 ∧ P2 ∧ P3: "the uncovered fraction ≈ the independent-residue model ≈ e^{−λ}" is
  **experimentally supported at 48–96 bits**.
- P4: "restricting to primes vs integers changes the rate by < 10%" holds.
- P5 fails, so the frozen rule does **not** let me claim "the 4th-moment main term is
  numerically accurate". For prime pairs it is accurate to ±1.5% in all 68 cells, but that is
  a post-hoc reading.
- P6 holds, with under-dispersion: per-q counts vary less than binomially.

### 5.3 Why P5 failed (post hoc, exploratory X1; not pre-registered)

- In the integer model, x ≡ (0,0) (mod ℓ) has probability 1/ℓ² and gives M_ℓ = |V|.
- With 10⁶ samples the number of such points is O(1)–O(10) per cell. Yet the model's fourth
  moment contains |V|⁴Σ_ℓ ℓ^{−2}, which is dominated by them. So the sample m4 is a
  heavy-tailed Poisson draw: mostly far below the model, occasionally above.
- This is the set 𝒟 of §2.5. The proof also treats it separately and bounds it crudely (row F).
- X1 recomputes, on the identical samples, m4 with 𝒟 removed against the model conditioned on
  x ≢ (0,0) mod every ℓ.

Result, over all 68 integer cells with λ ≥ 1 (`experiments/raw/sieve2d_explore.jsonl`):
- **off 𝒟**, m4/model ∈ [0.986, 1.007] and m6/model ∈ [0.960, 1.021]. That includes all 34
  cells that failed P5, whose off-𝒟 ratios lie in [0.986, 1.005].
- **#𝒟 per cell** is 0–439, and 0–25 in the failing cells.
- **Share of the sample's 4th-moment sum carried by 𝒟** ranges up to 0.999.
- **Reproduction:** on the frozen population the diagnostic reproduces the frozen m4 exactly
  (relative difference 0).

So the main-term identity is accurate to ~1% for integer pairs too, once the atom is
separated (**experimentally supported**, post hoc). The frozen P5 verdict stays **FAILED**.

*Disclosed deviation inside the exploratory script.* The first X1 run forgot to remove the
lines ℒ from integer samples, which the frozen E2 had removed. Its ratios reached 34× (m4)
and 4·10⁴× (m6). It was re-run with ℒ removed. Only the re-run is in the raw file.
- The pre-registration did not foresee the atom. The prediction was ill-posed for integer
  cells. It is still recorded as a **failure**.

### 5.4 Sixth moments (exploratory X2; not pre-registered)

Question: is the j = 6 obstruction (§4.2) visible, i.e. does the true 6th central moment of M
exceed the model?

Result (`sieve2d_explore.jsonl`, prime cells with λ ≥ 1, block 0, n = 48–96, 10⁶ pairs):
- **m6/model ∈ [0.961, 1.038]** over 33 cells, and m4/model ∈ [0.988, 1.015].
- **In the 12 cells with λ³/H² ≥ 1** (up to 18.6, at n = 80, H = 32) m6/model ∈ [0.990, 1.010].
- **Comparison with the proof's bound.** The j = 6 per-lattice bound summed without
  cancellation would allow a deviation of ≈ 32λ⁶/H², i.e. ≈ 2λ³/H² ≈ 37 times the model
  sixth moment (≈ 15λ³) in the largest cell. The data show ≤ 1%.
- So the signed 6-prime lattice discrepancies **do** cancel in reality, by at least 3 orders of
  magnitude in these cells (**experimentally supported**, exploratory).
- The j = 6 obstruction of §4.2 is therefore a barrier for this proof technique, not a
  feature of the true distribution. This is consistent with P1–P3 (e^{−λ} decay).

## 6. Reproduction

```
.venv/bin/python -m pytest -q tests/test_sieve2d.py tests/test_sieve2d_bound.py
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments config     # identical config
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments exh        # E1 (~2.5 min)
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments samp       # E2/E3 (~6 min)
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments summarize  # P1–P6
PYTHONPATH=. .venv/bin/python -m src.sieve2d.explore                # X1, X2 (exploratory)
```

## 7. Review and provenance

- **Blind GPT-6.1-sol attempt** (`codex/sieve2d_1.out`; it did not see these notes).
  - It independently proved the same results: the n⁴/H⁴ fourth-moment bound with explicit
    errors, the cost N^{1/6+5ε/12}, the O(1/λ) prime presieve, the PYY fake distribution, and
    the j = 6 failure.
  - It also supplied the integer-dilate example (§4.4) and the literature: Daniel 1999 (Crelle
    507) Lemma 3.3 and Shute 2025 (JIMJ) §3 for averaged planar lattice counting over shortest
    vectors; Brady, "Poisson imitators".
  - My proof (§2) was written in parallel and differs in its bookkeeping: an explicit case split
    on the shortest vector, explicit constants, and per-configuration lemmas tested by code.
    It borrows its shortest-vector averaging (row B) from that answer.
- **Orchestrator / r16 input** (prime presieve, the same-slope degeneracy): verified and
  incorporated (Proposition 4; Remark 1 in §2.5).
- **Independent Opus critic** (read-only; own scripts: 2,750 random configurations including
  case (ii) with HX near Y⁴, 2·10⁵ lattices for Lemma L, symbolic falling factorials,
  exhaustive boxes).
  - Verdict: **ACCEPT WITH FIXES.** Theorems 1–3 and Propositions 4–5 stand as proved.
  - Fixes, all applied:
    - (1) Corollary 3′ needed a different analysis-only T. D″'s T gives E_a = O(n^{−1−2a}),
      which is not o(n^{−4a}) for a ≥ 1/2.
    - (2) The PYY citation was corrected to Corollary 1.2(1)/1.3(1); Theorem 1.1 has a
      nonconstant factor.
    - (3) An absorption line was added to Corollary 1.
    - (4) §0 and §4.2–4.3 wording was softened (the j = 6 failure concerns this analysis, not
      the box; the effective-level Selberg claim is heuristic).
    - (5) Placeholders filled.
- **Adversarial GPT-6.1-sol referee** (`codex/sieve2d_2.out`, xhigh; 20,868 configurations and
  72 boxes, no violations).
  - Verdict: **ACCEPT WITH FIXES.** "The fourth-moment theorem and exponential factoring
    improvement survive."
  - Fixes, all applied:
    - the same Corollary 3′ gap (it proposed T = 2^{n/4}n^{−1−2a}; either choice works);
    - narrower barrier scope ("the displayed absolute-discrepancy estimates do not establish a
      sixth-moment improvement");
    - j = 5 downgraded to a sketch, with the extra O(|V|^{−3}) noted;
    - the fundamental-lemma quantification marked heuristic;
    - the Christoffel formula qualified by max π_ℓ → 0;
    - "primality must matter" narrowed to "a restriction excluding the dilates is necessary".
- **Labels.**
  - Proved: Theorems 1–3, Corollaries 1, 3′, Propositions 4, 5.
  - Experimentally supported: e^{−λ} model agreement (P1–P3), primes ≈ integers (P4), and
    6th-moment cancellation (X2, exploratory).
  - Failed as pre-registered: P5.
  - Conjectured: an e^{−cλ} bound for S_n; C/λ² for S_n via a presieved fourth moment.
  - Heuristic: §4.2's self-similarity description and §4.3's quantifications.
- **Novelty.** None claimed. The method is standard (lattice counting averaged over shortest
  vectors; Daniel, Shute). No prior formulation of this growing-H Farey-band bound was found
  by us or by GPT-6.1-sol.
