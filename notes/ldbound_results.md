# ldbound: Hypothesis LD, via a w-randomized design

Agent `ldbound`, phase 9, 2026-10-19.
- Code: `src/ldbound/` (`params.py`: exponents and constants; `checks.py`: exact identity checks).
  Tests: `tests/test_ldbound.py`.
- GPT-6.1-sol: blind round `codex/ldbound_1.{md,out}`; adversarial rounds `codex/ldbound_2.*` (§8).
- Opus critics: one on proofs, one on literature (§8).
- No prediction-testing experiment was run, so there is no pre-registration.

## Verdict

**DRAFT — reviews pending (§8). Labels below are the author's until the reviews are applied.**

**For a slightly modified design, the LD input of Theorem A is proved.** The modification adds one
random parameter to the field, and its infimum constant is the same as before.

In plain terms:
- **Why LD was hard.** In design (a″_k), the residue class of an algebraic value F_A(P) modulo the
  large rational value G(P) is *forced*: r₀ ≡ N(b/g₁)^d (mod G) (pcest Lemma 1.2).
  - Given the rational data, this residue is a point mass.
  - It is also coupled to the partner point's rational smoothness through g₁.
  - So the proof needed a large-values estimate for smooth numbers in a deterministic, structured
    family of residues. Nobody has such an estimate.
- **The fix.** Choose a random prime w ≈ G·N^{η₀}, independent of everything else, and build f with
  F(g₀, g₁) = N·w instead of N.
  - This is still an NFS pair: g₁^d f(g₀/g₁) ≡ 0 (mod N).
  - The forced residue becomes r₀ ≡ N·w·(b/g₁)^d. Its law over w is spread, by Brun–Titchmarsh, within
    a constant factor of uniform over the reduced classes mod G. This holds whatever the rational data
    are.
  - The structured-residue problem (critic piece (ii)) therefore disappears. Only an *all-residue*
    statement (piece (i)) remains.
- **Piece (i) is elementary at precision L^{o(1)}.** Because only an upper bound within L^{o(1)} is
  needed, the all-residue L² Barban–Davenport–Halberstam-type bound for the family ℓ²nq follows from
  two tools:
  - divisor switching (pcest Lemma DS) for conductors up to B^{C₀};
  - the additive large sieve for larger conductors. For conductors divisible by ℓ this uses the
    average over the special-q primes ℓ.
- **The cost.** w sits in the constant coefficient a₀ ≈ M·w. It is invisible in the algebraic norm as
  long as w ≤ H·N^{−η₀/2}, which needs h ≥ s + 1/(4k) + η₀.
  - This raises θ from 3s + 1/(4k) to 3s + d/(4k) at finite k.
  - The infimum over k is unchanged: c = (1 + 3/(d+1))/√2, which is 0.9192 at d = 9.
  - c < 1 at d = 9 now needs k ≥ 20 (before: 12). At d = 11 it needs k ≥ 17.

Labelled results:
1. **Design (a‴_k) and Lemma GB_w [proved; elementary]** (§1). FM, Theorem R, Lemma Kern, Lemma 1.3
   and Prop 4.2 are re-checked for (a‴_k) (§2).
2. **Lemma W (residue spreading) [proved; Brun–Titchmarsh]** (§3).
3. **Theorem BDH_fam (all-residue L² restricted BDH upper bound for moduli ℓ²nq) [proved; DS +
   large sieve]** (§4). This is critic piece (i), in the form the proof consumes.
4. **Theorem LD* (the input of Theorem 4.3) for (a‴_k) [proved]** (§5): E[N(N−1)] ≤ L^{o(1)}·S⁴ρ_r²ρ_F².
   With FM: PC ≤ L^{o(1)}, so D_k(γ) holds with γ ≥ L^{−o(1)}.
5. **Theorem A‴ [conditional on LV6 only, plus pcest's explicit caveats (HECM, or A′ without it; the
   Lemma 5.1 irreducibility transfer)]** (§6). Every N is factored in expected time
   L_N[1/2, (1 + θ‴(d,k))/√2 + O(η₀) + o(1)], with θ‴ = 3/(d+1) + d/(4k).
6. **Design (a″_k) itself: LD remains open** (§7). We give:
   - the exact missing correlation statement (Prop 7.1);
   - a proof that the averaging available in (a″_k) cannot reduce (ii) to (i) without a correlation
     input (Prop 7.2);
   - a barrier for the multiplier method: any residue randomizer needs support ≥ t^{1−o(1)} (Prop 7.3).
     So h ≥ s + 1/(4k) − o(1) is forced within this method.

**Is Theorem A's LD input now proved?**
- For design (a‴_k): **yes**, in the form Theorem 4.3 consumes (LD*), pending the reviews in §8.
- For design (a″_k): no, but (a″_k) is no longer needed. The constant's infimum is the same.

**What remains for Theorem A‴:** LV6 (agent `lv6`), now for Res(f, g₁x − g₀) = ±N·w rather than ±N.
The square-root step uses only f(g₀/g₁) ≡ 0 (mod N) and the size of ‖f‖ ≍ MH, which is unchanged.

### Constants (L[1/2, c], Lenstra–Pomerance shape test; `src/ldbound/params.py`)

| d | (a″_k): min k for c < 1 | (a‴_k): min k for c < 1 | (a‴_k), k = 30 | k = 50 | k = 100 | k = 1000 | inf_k |
|---|---|---|---|---|---|---|---|
| 7 | 13 | 45 | 1.0135 | 0.9970 | 0.9846 | 0.9735 | 0.9723 |
| 9 | 12 | 20 | 0.9723 | 0.9511 | 0.9351 | 0.9208 | 0.9192 |
| 11 | 13 | 17 | 0.9487 | 0.9228 | 0.9033 | 0.8858 | 0.8839 |
| 13 | 13 | 17 | 0.9352 | 0.9046 | 0.8816 | 0.8609 | 0.8586 |
| 19 | 16 | 18 | 0.9251 | 0.8803 | 0.8468 | 0.8165 | 0.8132 |

With Pollard–Strassen in place of HECM (Theorem A′):
- d = 9: 0.9929 at k = 100, 0.9695 at k = 1000, limit 0.9670;
- d = 11: 0.9033 → 0.9142 (limit; k = 100 gives 0.9425);
- d = 19: limit 0.8226.

---

## 1. Design (a‴_k)

Notation is pcest's (§1 there):
- L = L_N[1/2, 1], Λ = (log N log log N)^{1/2};
- s = 1/(d+1), M = ⌊N^s⌋, B = L[1/2, a], ϱ(x) = Ψ(x, B)/x;
- ℓ uniform among the primes of J = (Y, 2Y], 2Y = 2^j, with 2^{kj} ≍ √N; π_J := #(primes in J);
- S = L^{σ+o(1)}, V = ⌊√ℓ·S⌋, V_max := ⌊√(2Y)·S⌋;
- 𝒢 = coprime pairs in [M/2, M]², G_g(x, y) = g₁x − g₀y;
- 𝒫_ℓ = {(a, b): b prime in (V/2, V], |a| ≤ V, gcd(a, b) = 1}, Λ_ℓ(g) = {P : ℓ | G_g(P)};
- trial division to B as preprocessing.

Fix d odd ≥ 3, k > (d+1)/4, a > 0, a compact σ-range, η₀ ∈ (0, 1/(8k)).

**New parameters.**
- W := ⌈N^{s + 1/(4k) + η₀/2}⌉.
- **w** is uniform among the primes in (W, 2W], independent of (ℓ, g, c).
- h ≥ s + 1/(4k) + η₀ and H = ⌊N^h⌋. This implies pcest's floors 1/(2k) + η₀ (Kern), 1/(4k) + η₀ and
  s − (d−2)/(4k) + η₀ (FM). The first implication uses s > 1/(4k), i.e. k > (d+1)/4.

**Field law 𝒜‴.** (g, w, c) are independent: g uniform on 𝒢, w as above, c uniform on C_H = [−H, H]^d ∩ Z^d.
- F_A := F̂_{g,w} + G_g·R_c, where F̂_{g,w} is Lemma GB applied to the integer N·w.
- So F_A(g₀, g₁) = N·w and f(g₀/g₁) ≡ 0 (mod N) after clearing g₁^d.
- The map ψ: x ↦ g₀g₁^{−1} (mod N) is well defined on Z[x]/(f), because g₁^d f(g₀/g₁) = Nw ≡ 0.
- **Gcd checks.** If gcd(w, N) > 1, gcd(g₁, N) > 1 or gcd(G(P), N) > 1, a factor of N has been found.
  Since w is prime, the first case means w | N.
- Reducible f are rejected, as in pcest §1.
- Everything else is as in (a″_k): candidates, special-q lattices, the filter |G(P)| ≥ MV/32, the capped
  search, the counted quantity N_ℓ(A), the LP shape test and Theorem K's field selection.
- A uniform random prime in (W, 2W] is sampled by drawing uniform integers and testing primality with
  AKS. The expected number of trials is O(log W).

**Lemma GB_w [proved].** Apply pcest Lemma GB to T_d := Nw. Then:
- |a_i| ≤ g₁/2 for 1 ≤ i ≤ d;
- |a₀| ≤ 2^{d+2}Mw + d2^{d−1}M for large N;
- a_d ≢ 0 (mod g₁).
*Proof.*
- The recursion and identity are GB's, verbatim.
- For the bound, a₀ = (Nw − Σ_{i≥1} a_i g₀^i g₁^{d−i})/g₁^d. Here Nw/g₁^d ≤ 2^d Nw/M^d ≤ 2^{d+1}Mw·(N/M^{d+1})
  ≤ 2^{d+2}Mw, and each |a_i|(g₀/g₁)^i ≤ 2^{i−1}M.
- a_d ≡ Nw g₀^{−d} (mod g₁), and gcd(g₁, Nw) = 1 because w > M ≥ g₁ is prime and gcd(g₁, N) = 1. ∎
(Tested: `field_w`, `test_field_w_resultant_and_lemma12`.)

**Sizes [proved arithmetic].** For P ∈ 𝒫_ℓ:
- |F̂_{g,w}(P)| ≤ (dM/2 + |a₀|)V^d ≤ C_d·MwV^d ≤ C_d·MHV^d·N^{−η₀/2}, using 2W ≤ 3N^{s+1/(4k)+η₀/2} and
  H ≥ N^{s+1/(4k)+η₀}/2;
- |G(P)⟨c, v(P)⟩| ≤ 2dMHV^d.
- So |F_A(P)| ≤ X_F := C′_d MHV^d, the same order as in (a″_k).
- The rational side is untouched: |G/ℓ| ≤ X_r := 2MV/ℓ.
- Hence ρ_r, ρ_F and θ = 2s + h + (d−1)/(4k) are as before, with the new h.
- At h = s + 1/(4k) + η₀ this gives **θ‴ = 3s + d/(4k) + η₀**.

**Lemma 1.2_w [proved].** g₁^d F_A(P) ≡ Nw·b^d (mod G(P)). Same proof as pcest Lemma 1.2. (Tested.)

**Coprimality [proved].** Let P ∈ 𝒫_ℓ ∩ Λ_ℓ(g) with gcd(G(P), N) = 1, and put t := |G(P)| ≤ 2MV.
- gcd(b, t) = gcd(g₁, t) = 1, by pcest FM Step 2(a).
- gcd(w, t) = 1, since w > W > 2MV ≥ t and w is prime.
- Hence the moving residue r₀ := F̂_{g,w}(P) mod t satisfies
    r₀ ≡ c_P·w (mod t),  where c_P := N b^d g₁^{−d} ∈ (Z/t)^×.   (1.1)
- **c_P depends only on (ℓ, P, g), not on w.**

## 2. Re-check of the proved pcest results for (a‴_k)

- **Lemma 1.3, Lemma Q, Lemma DS, Theorem R, Lemma Kern, Prop 4.2.** These involve only the rational
  form G_g, the candidate set and the linear map c ↦ (⟨c, v⟩, ⟨c, v′⟩). None of them involves N, w or
  F̂. Lemma Kern needs H ≥ V², which holds. **Unchanged** [proved there].
- **Theorem FM for (a‴_k) [proved].** pcest's proof applies for each fixed admissible w, with constants
  uniform in w; averaging over w then gives FM for 𝒜‴. Line by line:
  - *Step 1 (fibres).* Rational only. Unchanged.
  - *Step 2(a).* Use (1.1) in place of gcd(r, t) = gcd(N, t) = 1.
  - *Steps 2(b), 2(c).* These concern c only. Unchanged; they need h > 1/(4k) + η₀ and H > dωV, both of
    which hold.
  - *Step 2(d).*
    - The sign argument needs |F̂| ≤ ty/4 for y ≥ Y_ω = (2H+1)V′^{d−1}/ω. Now |F̂|/t ≤ 32C_d wV^{d−1}
      (t ≥ MV/32), and wV^{d−1}/Y_ω ≤ 4ωw/H ≤ L^{O(1)}N^{−η₀/2} → 0.
    - The additive error 3(|F̂|/t + 1)/(y + 1) is ≤ C·ωw/H = N^{−η₀/2}L^{O(1)}. This is polynomially
      small against ρ_F, as before (pcest had ωV/H there).
    - The grid 𝒳 uses X_F^max = C′_d MH V_max^d. It is still independent of t, P, ℓ, w and the residue.
  - *Step 3 (Lemma 8.5′).* The moduli are t ∈ [MV/32, 2MV], so q = s + 1/(4k) + o(1), and
    x = X_F = N^{ξ+o(1)} with ξ = s + h + d/(4k).
    - q > 1/(2k) + η₁ holds, since k > (d+1)/4.
    - q < ξ/2 ⟺ h > s − (d−2)/(4k). This holds with margin (d−1)/(4k) + η₀.
    - Goodness is a statement about all reduced residues, so it does not see w.
  - *Step 4 (assembly).* Unchanged.
  - Conclusion: E_{ℓ,A∼𝒜‴} N_ℓ(A) ≥ S²ρ_rρ_F/(log N)^{O(1)}.
- **Irreducibility.** The LV Lemma 5.1 transfer reads exactly as in pcest §1, with the affine family
  f̂_{g,w} + (g₁x − g₀)Z[x]_{d−1} and the value condition ĝ(g₀, g₁) ≡ Nw (mod r) for r ∤ g₁Nw.
  - Status unchanged: checked at the level of statements, and listed as an assumption.
  - A reducible f may now split as F₁(g₀, g₁)F₂(g₀, g₁) = Nw without factoring N; such f are rejected
    anyway, and their probability is ≤ H^{−1/3+o(1)}.
- **Relations.** For P = (a, b) with G(P) = ±ℓn and both values smooth, ψ(a − bα) = G(P)/g₁ (mod N),
  exactly as before.

## 3. Lemma W (the forced residue becomes spread)

**Cited.**
- **(BT)** Montgomery–Vaughan (Mathematika 20 (1973) 119–134, Thm 2): for 1 ≤ q < x and (a, q) = 1,
  π(x; q, a) ≤ 2x/(φ(q) log(x/q)).
- **(PNT)** π(2W) − π(W) ≥ W/(2 log W) for W ≥ W₀.

**Lemma W [proved].** For every t ∈ [1, 2MV_max] and every a ∈ Z,
  Pr_w[w ≡ a (mod t)] ≤ C_W/φ(t),  with C_W := 24/η₀,
for all N ≥ N₀(d, k, σ, η₀).

*Proof.*
- If (a, t) > 1 the probability is 0, since w > t is prime.
- Otherwise #{w ∈ (W, 2W] prime : w ≡ a (mod t)} ≤ π(2W; t, a) ≤ 4W/(φ(t) log(2W/t)).
- Divide by W/(2 log W). This gives ≤ 8 log W/(φ(t) log(2W/t)).
- Now log W ≤ (s + 1/(4k) + η₀) log N ≤ log N. Since 2MV_max ≤ N^{s+1/(4k)}L^{O(1)},
  log(2W/t) ≥ (η₀/2) log N − O(Λ) ≥ (η₀/3) log N. ∎
(Exact toy check: `bt_ratio`, `test_brun_titchmarsh_ratio`.)

**Corollary W.** Let Φ: Z/t → [0, ∞) be any function. Then, by (1.1),
  E_w Φ(r₀) ≤ (C_W/φ(t)) Σ_{c ∈ (Z/t)^×} Φ(c).
This holds for any (ℓ, P, g) with gcd(G(P), N) = 1. ∎

This replaces the large-values problem (critic piece (ii)). The bound is uniform in the rational data,
so the rational weights can be summed separately (§5, Step 5).

## 4. Theorem BDH_fam (critic piece (i), all-residue, at precision L^{o(1)})

**Notation.**
- X* := X_F^max = C′_d MH V_max^d. It is common to all ℓ ∈ J and dominates every window below.
- 𝒮* := 𝒮_B ∩ [−X*, X*] ∖ {0}, with #𝒮* ≤ 3Ψ(X*, B).
- S(α) := Σ_{F∈𝒮*} e(αF).
- For m ≥ 1: 𝔔(m) := #{(F, F′) ∈ 𝒮*² : F ≡ F′ (mod m)}.
- X′₊ := 2MV_max/Y. For ℓ ∈ J: 𝒩_ℓ := {n ∈ Z : P⁺(|n|) ≤ B, ℓ|n| ∈ [MV/32, 2MV]}. Then |n| ≤ X′₊ on 𝒩_ℓ.

**Additive Parseval [proved; exact].**
  𝔔(m) = (1/m) Σ_{h mod m} |S(h/m)|² = (1/m) Σ_{m′|m} f(m′),  f(m′) := Σ_{a mod m′, (a,m′)=1} |S(a/m′)|²,
since each h/m is a unique reduced fraction a/m′ with m′ | m. The term m′ = 1 is |S(0)|²/m = (#𝒮*)²/m.
(Integer form via Ramanujan sums: `pair_count_ramanujan`, tested.)

**Cited.**
- **(LS)** The additive large sieve (Montgomery–Vaughan, J. London Math. Soc. 8 (1974); Selberg): for
  any a_F supported on an interval of K integers and Q ≥ 1,
    Σ_{m′≤Q} Σ_{(a,m′)=1} |Σ_F a_F e(aF/m′)|² ≤ (K + Q²) Σ|a_F|².
- **(DS)** pcest Lemma DS, for every residue (non-reduced included), modulus ≤ B^{C₀}, and intervals of
  length ≥ XB^{−C₀} inside [XB^{−C₀}, X].
- **(Hi), (C8.16)** as in pcest §2.1:
  - ϱ(z₂) ≤ 2ϱ(z₁) for z₂ ≥ z₁;
  - ϱ(x) = ϱ(xB^c)(log x)^{O(⌈c⌉)} for fixed c ≥ 0.

Fix C₀ with B^{C₀} ≥ 64S²·L/(ϱ(X′₊)ϱ(X*)) for large N. This is possible because S², 1/ϱ(X′₊) and
1/ϱ(X*) are L^{O(1)}. Put R₀ := B^{C₀}.

**Theorem BDH_fam [proved].** Uniformly for 1 ≤ q ≤ 2S²,
  (***)  E_ℓ Σ_{n∈𝒩_ℓ} Σ_{m′ | ℓ²|n|q, m′>1} f(m′) ≤ τ(q)(log N)^{O(1)}·Ψ(X′₊, B)·Ψ(X*, B)².
Consequently, since #𝒩_ℓ ≤ 3Ψ(X′₊),
  E_ℓ Σ_{n∈𝒩_ℓ} 𝔔(ℓ²|n|q)/φ(ℓ|n|) ≤ τ(q)(log N)^{O(1)}·Ψ(X′₊)Ψ(X*)²/(Y q (MV)²).   (4.1)

*Proof.*
- **Rearrangement.** Write the left side of (***) as Σ_{m′>1} f(m′)·ν(m′), where
  ν(m′) := E_ℓ #{n ∈ 𝒩_ℓ : m′ | ℓ²|n|q}. Recall ℓ > max(B, 2S²) ≥ max(P⁺(n), q), so ℓ ∤ |n|q.
  Split the m′ into three classes.
- **(A) ℓ ∤ m′ and 1 < m′ ≤ R₀.**
  - *Bound on f.* f(m′) ≤ Σ_{a mod m′} |S(a/m′)|² = m′ Σ_{c mod m′} N_c², with N_c := #{F ∈ 𝒮* : F ≡ c (mod m′)}.
  - For each c, DS (for the positive and the reflected negative parts of [X*R₀^{−1}, X*]), together
    with the trivial count ≤ X*R₀^{−1}/m′ + 1 on |F| < X*R₀^{−1}, gives N_c ≤ (X*/m′)ϱ(X*)(log N)^{O(1)}.
    This uses R₀^{−1} ≤ ϱ(X*) and m′ ≤ X*ϱ(X*).
  - Hence f(m′) ≤ m′·max_c N_c·#𝒮* ≤ Ψ(X*)²(log N)^{O(1)}.
  - *Bound on ν.* Here m′ | |n|q. Put m₂ := gcd(m′, q) and m₁ := m′/m₂, so m₁ | n. Then
    ν(m′) ≤ 2Ψ(X′₊/m₁) ≤ 2(X′₊/m₁)·2ϱ(X′₊R₀^{−1}) ≤ (Ψ(X′₊)/m₁)(log N)^{O(C₀)}, by (Hi) and (C8.16).
  - The map m′ ↦ (m₁, m₂) is injective, with m₂ | q. Hence Σ_{m′≤R₀} 1/m₁ ≤ τ(q)(1 + log R₀).
  - So (A) ≤ τ(q)(log N)^{O(1)}Ψ(X′₊)Ψ(X*)².
- **(B) ℓ ∤ m′ and m′ > R₀.**
  - Necessarily m′ ≤ |n|q ≤ X′₊q.
  - Take dyadic blocks m′ ∈ (R, 2R], R = R₀2^i. On a block, ν(m′) ≤ #{1 ≤ |n| ≤ X′₊ : m₁ | n} ≤ 2X′₊/m₁ ≤
    2X′₊q/R, since m₁ ≥ m′/q.
  - By (LS) with K = 2X* + 1 and Q = 2R: Σ_{m′∈(R,2R]} f(m′) ≤ (3X* + 4R²)·3Ψ(X*).
  - Summing over blocks:
      (B) ≤ Σ_R (2X′₊q/R)(3X* + 4R²)·3Ψ(X*) ≤ 36X′₊qX*Ψ(X*)/R₀ + 48X′₊qΨ(X*)·(4X′₊q).
  - First term: 36X′₊qX*/R₀ ≤ Ψ(X′₊)Ψ(X*), because R₀ ≥ 36q/(ϱ(X′₊)ϱ(X*)).
  - Second term: 192X′₊q² ≤ ϱ(X′₊)ϱ(X*)X* holds for large N. Indeed X*/X′₊ ≥ c_d HV^{d−1}Y ≥ N^{h}, and
    q²/(ϱϱ) = L^{O(1)}.
  - So (B) ≤ 2Ψ(X′₊)Ψ(X*)².
- **(C) ℓ | m′.**
  - Write m′ = ℓ^i m″, with i ∈ {1, 2} and m″ | |n|q. The prime ℓ is determined by m′: it is m′'s unique
    prime factor > max(B, 2S²).
  - So only one ℓ contributes to ν(m′), and ν(m′) = π_J^{−1}#{n ∈ 𝒩_ℓ : m″ | |n|q}.
  - Two bounds hold: ν(m′) ≤ π_J^{−1}·3Ψ(X′₊), and, with m″ ≥ m′/ℓ² > R/(4Y²) on a block (R, 2R],
    ν(m′) ≤ π_J^{−1}·2X′₊q/m″ ≤ π_J^{−1}·8X′₊qY²/R.
  - The blocks run over Y/2 ≤ R ≤ R_max := 4Y²X′₊q, which is O(log N) blocks.
  - By (LS):
      (C) ≤ π_J^{−1} Σ_R [3Ψ(X′₊)·3X*·3Ψ(X*) + (8X′₊qY²/R)·4R²·3Ψ(X*)]
          ≤ π_J^{−1}[27 log N·Ψ(X′₊)X*Ψ(X*) + 192X′₊qY²R_maxΨ(X*)].
  - First term: X*/π_J ≤ Ψ(X*)/log² N, because π_J ≥ Y/(2 log Y) = N^{1/(2k)+o(1)} while 1/ϱ(X*) = L^{O(1)}.
  - Second term: π_J^{−1}·768X′₊²q²Y⁴ ≤ Ψ(X′₊)Ψ(X*), provided X′₊q²Y³ log Y ≤ c·ϱ(X′₊)ϱ(X*)X*. Since
    X*/X′₊ ≥ c_d HV^{d−1}Y, this reduces to Y²q² log Y ≤ c·ϱϱ·HV^{d−1}.
    - In exponents: 1/k < h + (d−1)/(4k), i.e. h > (5−d)/(4k).
    - This holds with polynomial margin, since h > s + 1/(4k) > 1/(2k) ≥ (5−d)/(4k) for d ≥ 3.
  - So (C) ≤ 2Ψ(X′₊)Ψ(X*)².
- **(4.1) from (***).**
  - For n ∈ 𝒩_ℓ, m := ℓ²|n|q ≥ ℓqMV/32 and φ(ℓ|n|) ≥ c ℓ|n|/log log N ≥ c′MV/log log N.
  - By the Parseval identity, 𝔔(m)/φ(ℓ|n|) ≤ (C log log N/(Yq(MV)²))·[(#𝒮*)² + Σ_{m′|m, m′>1} f(m′)].
  - Sum over n, average over ℓ, and use #𝒩_ℓ ≤ 3Ψ(X′₊) and (***). ∎

*Remarks.*
- **What is new, and what is not.** No equidistribution (no Siegel–Walfisz input) is needed. LD only
  asks for an upper bound within L^{o(1)}, so each range can be handled crudely:
  - conductors ≤ B^{C₀} = L^{O(1)} with a uniform pointwise bound (DS), for any fixed C₀;
  - all larger conductors with the large sieve, whose diagonal X*·Ψ is beaten either by the sparsity
    ν(m′) ≲ X′₊q/m′ of the moduli divisible by m′ (class B), or by the average over the N^{1/(2k)} primes
    ℓ (class C).
- **The ℓ-average is essential in (C).** For one fixed ℓ, the characters to the modulus ℓ² alone give
  Σ_{a mod ℓ²} |S(a/ℓ²)|² = ℓ²·Σ_c N_c², and the trivial bound for it is X*Ψ(X*) ≫ Ψ(X*)².
  - Proving the single-modulus bound Σ_c N_c² ≲ L^{o(1)}Ψ²/ℓ² for one modulus ℓ² = N^{1/k} would need
    smooth numbers to be L²-equidistributed modulo one modulus of size B^{Θ(u)}.
  - Unconditionally this is out of reach. For example, the primes ≤ B could all lie in a subgroup of
    (Z/ℓ²)^× of index up to 1/ϱ(ℓ², B) = L^{Θ(1)}: Burgess-type bounds on the least non-residue are far
    above B = ℓ^{o(1)}.
  - The average over ℓ ∈ J removes this.
- **The family's sparsity is used only through ν(m′) ≤ 2X′₊q/m₁.** This is a divisibility count. No
  smoothness of the moduli is needed beyond (A).
- **Prop B5 is not contradicted.** BDH_fam is an L² statement: all residue classes at once, weighted by
  their own mass. It never passes through max-class ratios.

## 5. Theorem LD* for (a‴_k)

**Theorem LD* [proved].**
- **Hypotheses.** Those of Theorem FM for (a‴_k) (§2), and S² = (log N)^{C₃}/(ρ_rρ_F) as in pcest
  Theorem 4.3. Let N_ℓ(A) be the filtered count.
- **Conclusion.**
    E[N(N−1)] ≤ τ*(S)(log N)^{O(1)}·S⁴ρ_r²ρ_F²,  τ*(S) := max_{q≤2S²} τ(q) = L^{o(1)}.
- Equivalently, in pcest's normalization (Prop 4.2 and ϱ(T|G|) ≥ ρ_F/(log N)^{O(1)} under the filter):
    E_{ℓ,w} Σ_{P≠P′, ℓ|Δ} |𝒢|^{−1} Σ_{g, filtered} Y_PY_{P′}(S²/q)κ_alg ≤ L^{o(1)}·S⁴ρ_r².
- With Theorem FM: PC ≤ L^{o(1)}, and D_k(γ) holds with γ ≥ L^{−o(1)} (Prop PZ).

*Proof.* Write T := dHV^{d−1}, t := |G_g(P)|, t′ := |G_g(P′)| and q := |Δ|/ℓ ∈ [1, 2S²]. Sums over pairs run
over ordered pairs P ≠ P′ in 𝒫_ℓ ∩ Λ_ℓ(g) with ℓ | Δ and both |G|, |G′| ∈ [MV/32, 2MV].

- **Step 0.** E[N(N−1)] = E_{ℓ,g,w} Σ_{P≠P′} Y_PY_{P′}·Pr_c[F_A(P), F_A(P′) ∈ 𝒮_B].
  - Pairs with ℓ ∤ Δ contribute nothing (pcest Lemma Q, step 0).
  - Every quantity except F̂_{g,w} is independent of w.
- **Step 1 (kernel).**
  - By Lemma Kern(iii) with the effective wedge bound ‖v ∧ v′‖ ≥ |Δ|(V/2)^{2d−4}, Pr_c[Φ(c) = z] ≤
    C_d/(H²V^{2d−4}) for every z.
  - Φ(C_H) ⊂ [−T, T]², and Γ = Φ(Z^d) = {(z₁, z₂) : z₂ ≡ κ^{d−1}z₁ (mod Δ)} (Kern(i)). Hence
      Pr_c[both smooth] ≤ (C_d/(H²V^{2d−4}))·Z,  Z := Σ_{j mod |Δ|} A_j(P)A_{κ^{d−1}j}(P′),
    where A_j(P) := #{z ∈ [−T, T] : z ≡ j (mod |Δ|), F̂(P) + G(P)z ∈ 𝒮_B}.
- **Step 2 (AM–GM).**
  - κ is a unit mod Δ (P′ is primitive), so j ↦ κ^{d−1}j is a bijection of Z/|Δ|.
  - Hence Z ≤ ½(Q(P) + Q(P′)), where Q(P) := Σ_j A_j(P)². Q(P) depends on (P, g, w, |Δ|) only.
- **Step 3 (enlargement to a residue-only quantity).**
  - z ↦ F = F̂(P) + G(P)z is injective, F ≡ r₀ (mod t), and z ≡ z′ (mod |Δ|) ⟺ F ≡ F′ (mod t|Δ|).
  - Every such F lies in [−X*, X*], by §1 sizes. Hence
      Q(P) ≤ 𝔔₀(t, |Δ|; r₀) := #{(F, F′) ∈ 𝒮*² : F ≡ F′ ≡ r₀ (mod t), F ≡ F′ (mod t|Δ|)}.
  - 𝔔₀ depends on w only through r₀ = F̂_{g,w}(P) mod t.
- **Step 4 (average over w).**
  - By Corollary W, E_w 𝔔₀(t, |Δ|; r₀) ≤ (C_W/φ(t)) Σ_{c∈(Z/t)^×} 𝔔₀(t, |Δ|; c).
  - The classes c partition the pairs counted, so this is ≤ (C_W/φ(t))·𝔔(t|Δ|).
  - Likewise for P′, with t′.
  - The pair sum is symmetric under (P, P′) ↦ (P′, P): Δ ↦ −Δ, and the filter and the Y's are symmetric.
  - Therefore
      E[N(N−1)] ≤ (C_dC_W/(H²V^{2d−4}))·E_{ℓ,g} Σ_{P≠P′} Y_PY_{P′}·𝔔(t ℓq)/φ(t).   (5.1)
- **Step 5 (the rational weights; fibres).**
  - Fix ℓ, P = (a, b) and P′ with |Δ| = ℓq. The summand of (5.1) depends on g only through n := G_g(P)/ℓ,
    and Y_P = 1 ⟺ n ∈ 𝒩_ℓ.
  - For fixed n, the g ∈ 𝒢 with G_g(P) = ℓn lie on the line g = g* + j·(a, b) in (g₀, g₁)-coordinates.
    The direction is primitive and |(a, b)| ≥ b > V/2, so j ranges over ≤ √2M/V + 1 integers.
  - Along it, G_g(P′) changes by ba′ − ab′ = −Δ per step. So n′ := G_g(P′)/ℓ runs through an arithmetic
    progression of difference q, inside an interval of length A := q(√2M/V + 1).
  - The filter puts |n′| ∈ [X_r/64, X_r], with X_r := 2MV/ℓ.
  - So for each sign the n′-values lie in an interval I_± ⊂ [X_r/64, X_r] with |I_±| ≤ A, in a fixed
    class mod q.
  - Apply DS with X = X_r and modulus q ≤ 2S² ≤ B^{C₀}. Enlarge I_± to an interval of length
    max(|I_±|, X_rB^{−C₀}) inside [X_rB^{−C₀}, X_r]; positivity allows this. Since X_rB^{−C₀}/q ≤ M/V,
    this gives
      #{g ∈ 𝒢 : G_g(P) = ℓn, Y_{P′} = 1, P′ filtered} ≤ U := (3M/V)ρ_r(log N)^{O(1)},
    uniformly in ℓ, P, P′ and n.
  - Hence Σ_g Y_PY_{P′}𝔔(tℓq)/φ(t) ≤ U·Σ_{n∈𝒩_ℓ} 𝔔(ℓ²|n|q)/φ(ℓ|n|).
- **Step 6 (partners).**
  - For fixed P and q, at most 12 points P′ ∈ 𝒫_ℓ have det(P, P′) = ±ℓq.
  - They lie on two lines parallel to P, on which consecutive integer points differ by P, with |P| ≥ V/2.
    Each line meets [−V, V]² in a segment of length ≤ 2√2V.
- **Step 7 (assembly).**
  - By (5.1), Steps 5–6 and (4.1), using |𝒢| ≥ cM², #𝒫_ℓ ≤ CV_max²/log V and V_max/V ≤ 2:
      E[N(N−1)] ≤ C·(H²V^{2d−4})^{−1}·M^{−2}·V²·(M/V)ρ_r·Σ_{q≤2S²} τ(q)(log N)^{O(1)}Ψ(X′₊)Ψ(X*)²/(Yq(MV)²).
  - Here Ψ(X′₊) ≤ (2MV/Y)ρ_r(log N)^{O(1)} and Ψ(X*) ≤ C MHV^d ρ_F(log N)^{O(1)}, by (Hi) and (C8.16), and
    Σ_q τ(q)/q ≤ (log S)².
  - Collecting powers:
      V²·(M/V)·(MV/Y)·(MHV^d)²/(H²V^{2d−4}·M²·Y·(MV)²) = V⁴/Y² ≍ S⁴.
  - Hence E[N(N−1)] ≤ τ*(S)(log N)^{O(1)}S⁴ρ_r²ρ_F².
  - **The normalized form.** (S²/q)κ_alg = Z/(4d²H²V^{2d−4}ϱ(T|G|)ϱ(T|G′|)), and ϱ(T|G|) ≥ ρ_F/(log N)^{O(1)}
    under the filter. Steps 2–7, applied to Z, give the displayed bound.
  - **PC.** By Theorem FM, (E N)² ≥ S⁴ρ_r²ρ_F²/(log N)^{O(1)} and E N ≥ 1. Hence PC ≤ τ*(S)(log N)^{O(1)} =
    L^{o(1)}, since τ(q) ≤ exp(O(log q/log log q)) and log q = O(Λ). ∎

*Remarks.*
- **Theorem R is not used.** The rational pair correlation is re-derived inside Step 5, by the same
  divisor-switching argument, because the algebraic factor no longer depends on the fibre variable.
- **No large-values input appears.** The only place where the residue enters is Step 4, and there it is
  integrated out over w before any arithmetic estimate is applied.
- **Uniformity.** All constants depend only on (d, k, a, the σ-range, η₀). Every "for large N" is a
  polynomial margin (`src/ldbound/params.py: margins`) or an L^{O(1)} comparison.
- **LD versus LD\*.** pcest's literal LD has right side L^{o(1)} × (rational pair count). We prove the
  bound Theorem 4.3 consumes, L^{o(1)}S⁴ρ_r². The two agree if the rational pair count is
  ≥ S⁴ρ_r²L^{−o(1)}; this lower bound is not needed and we do not prove it.

## 6. Theorem A‴

**Theorem A‴ [conditional on LV6; plus pcest's explicit caveats].**
- **Hypotheses.** Fix odd d ≥ 3, k > (d+1)/4 and η₀ > 0, and put θ‴ := 3/(d+1) + d/(4k) + η₀ ≤ 1.
  Let h := 1/(d+1) + 1/(4k) + η₀, a := θ‴/√2 and B = L[1/2, a].
- Assume:
  - LV6 for design (a‴_k): LV §6 adapted to pure-descent rows at height V, the non-monic pair
    (f, g₁x − g₀) with Res = ±Nw, and the g₁-parity column;
  - HECM's running-time theorem, as in pcest (Theorem A′ removes this assumption, at the constants of
    the table);
  - the LV Lemma 5.1 irreducibility transfer, at the level of statements.
- **Conclusion.** Every N is factored in expected time L_N[1/2, (1+θ‴)/√2 + o(1)].
- Infimum over k: (1 + 3/(d+1))/√2, i.e. d = 9: 0.9192, d = 11: 0.8839, d = 19: 0.8132.
- c < 1 holds for d = 9 from k = 20 on, for d = 11 and 13 from k = 17, and for d = 7 from k = 45.

*Proof.* pcest Theorem A's proof (steps 1–5) applies verbatim, with design (a‴_k):
- Theorem 4.3(ii) is replaced by Theorem LD* + FM + Prop PZ;
- T_SQ = S²L^{o(1)} = L^{θ‴/(2a)+o(1)}.
- The new field parameter w costs O(log W) primality tests per field, i.e. (log N)^{O(1)}. ∎

## 7. Design (a″_k): what remains open, and why the averaging does not suffice

**Prop 7.1 (the exact missing statement for (a″_k)) [proved reduction; the statement itself open].**
Steps 0–3 and 5–7 of §5 hold for (a″_k) verbatim; take w = 1, and X* of the old size. Step 4 is the only
step that uses w. Hence LD* for (a″_k) follows from the following estimate, with c_P = N b^d g₁^{−d} and
g₁ = g₁(j) along the fibre:

    E_ℓ Σ_{P, q} Σ_{n∈𝒩_ℓ} Σ_{j : Y_{P′}(j)=1} 𝔔₀(ℓ|n|, ℓq; c_{P,j})
        ≤ L^{o(1)}·E_ℓ Σ_{P,q} Σ_{n∈𝒩_ℓ} U·𝔔(ℓ²|n|q)/φ(ℓ|n|).

This is a *correlation* between two objects along the fibre:
- the smoothness of the affine sequence n′(j) = n′₀ ∓ qj;
- the lift-L² mass 𝔔₀ at the residues c_P(j) = N b^d (g₁* + jb)^{−d} (mod ℓn), a "Möbius progression"
  of d-th powers.

**Prop 7.2 (what the averaging in (a″_k) gives) [proved].**
- (a) For fixed rational data (ℓ, P, P′, g), the residue r₀ is deterministic.
  - So in (a″_k), no average that keeps the rational weights Y_PY_{P′} inside can replace Step 4 by an
    average over residues. The only randomness left beyond the rational data is c, and c does not move
    r₀.
- (b) Dropping the partner weight (Y_{P′} ≤ 1) decouples the residue from the partner's smoothness, but
  costs exactly the factor 1/ρ_r.
  - In Step 5 the fibre count becomes ≤ √2M/V + 1 instead of U = (M/V)ρ_r(log N)^{O(1)}.
  - Every other step is unchanged at best, so the bound obtained exceeds S⁴ρ_r²ρ_F² by 1/ρ_r = L^{Θ(1)}.
    No estimate on the residue side can recover this loss.
- (c) Keeping Y_{P′} requires the correlation of Prop 7.1. Two tools fail:
  - Lemma DS cannot see residues mod t.
  - The large sieve controls Σ_{all classes}, not a ρ_r-sparse, j-dependent subfamily of classes.

  We know no tool, published or here, that gives Prop 7.1 [assessment].

**Prop 7.3 (barrier for the multiplier method: the randomizer must be as large as the modulus)
[proved; elementary].**
- **Hypotheses.** Let a design randomize the forced residue by F(g₀, g₁) = N·w, with w drawn from a law
  μ that is independent of the rational data. Suppose Step 4 is to hold with loss λ for *every* nonnegative
  Φ on Z/t, i.e. E_w Φ(c_P w) ≤ (λ/φ(t))Σ_c Φ(c).
- **Conclusions.**
  - μ(w ≡ a (mod t)) ≤ λ/φ(t) for all a.
  - Hence #supp(μ) ≥ φ(t)/λ, so max|w| ≥ φ(t)/(2λ) = t^{1−o(1)} when λ = L^{o(1)}.
  - F(g₀, g₁) = Nw with g₀, g₁ ≤ M forces max_i |a_i| ≥ Nw/((d+1)M^d) ≥ Mw/(d+1).
  - By equivalence of norms on binary forms of degree d, restricted to the region |a| ≤ V, V/2 < b ≤ V,
    this gives sup_P |F̂(P)| ≥ c_d MwV^d.
  - FM's window argument (Step 2(d)) needs |F̂(P)| ≪ t·HV^{d−1}. The algebraic value bound X_F used by
    FM and LD* is ≥ sup_P |F̂(P)|. So keeping X_F = O(MHV^d) requires H ≥ t^{1−o(1)}, i.e.
    h ≥ s + 1/(4k) − o(1).
- *Proof.* Take Φ = indicator of the class c_P a. The bound on #supp follows. The size statement follows
  from Σ a_i g₀^i g₁^{d−i} = Nw and g₀, g₁ ≤ M. ∎
- **Consequence [assessment: a barrier to this method, not to LD].** Within "independent residue
  randomization + all-residue BDH", the θ cost
  (d−1)/(4k) of (a‴_k) over (a″_k) is necessary. It vanishes as k → ∞.

**Status of the critic's route.**
- (i): proved, in the all-residue form consumed (Theorem BDH_fam).
- (ii): not proved for (a″_k). It is made unnecessary by the design change (Lemma W).
- LV*: not needed.

## 8. Reviews

*(pending)*

## 9. Candidate cards

```text
ID: KD-3 bounded-k Dixon descent, design (a'''_k): (a''_k) plus a random prime multiplier w ~ N^{s+1/(4k)+eta}
  of the resultant, F(g0, g1) = N w.
Target: every N (after trial division to B).
Proved: Lemma GB_w; FM, Theorem R, Lemma Kern re-checked; Lemma W (Brun-Titchmarsh spreading of the
  forced residue N w (b/g1)^d mod G); Theorem BDH_fam (all-residue L^2 restricted BDH upper bound for
  moduli l^2 n q, precision L^{o(1)}, via Lemma DS + additive large sieve + average over l); Theorem LD*.
Conditional on: LV6 (LV section 6 at height N^{1/(4k)}, Res = N w), HECM (or Pollard-Strassen, A'),
  irreducibility transfer.
Cost: L[1/2, (1 + 3/(d+1) + d/(4k))/sqrt 2]; inf over k (1 + 3/(d+1))/sqrt 2 (d = 9: 0.919);
  c < 1 for d = 9 from k = 20.
Status: DRAFT pending reviews.

ID: B-mult (barrier) multiplier randomization needs support >= t^{1-o(1)}, hence h >= s + 1/(4k).
Status: proved (elementary). Scope: designs that randomize the forced residue through Res(f, g) = N w.
```
