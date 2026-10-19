# Theorem A in design (a‴_k; μ): statement, algorithm and proof roadmap

Agent `thmA`, phase 10, 2026-10-19. This document assembles a claimed theorem from five internal,
unrefereed project notes. It is a roadmap that cites them; it does not re-prove every lemma. Each step below
names the exact result it uses and the place where that result is proved, and it states the conditions under
which that result is applied in the final design. Assessments of correctness belong in `notes/thmA_audit.md`,
not here.

Sources (all in `notes/`):
- **[DD]** `lemmaM_results.md` §1: Theorem DD and remarks R1–R4.
- **[hypD]** `hypD_results.md`: §2 Lemma Ω, Theorems S, S′; §3 Theorem K(i)–(iii), Corollary K; §4 Prop PZ.
- **[pcest]** `pcest_results.md`: §1 design, Lemmas GB, 1.2, 1.3, sizes; §2.1 cited facts; §2.2 Lemma 8.5′;
  §2.3 Theorem FM; §3 Lemma DS; §4.1 Lemma Kern.
- **[lv6]** `lv6_results.md`: §1 rows and (1.1); §2 Prop LV6 and Algorithm LV6; §3.1–3.7 Lemmas C, H, E, X,
  linear algebra, SQ, congruence, (T); §4 field events and Lemma IR; §6.2 Lemma BS, Prop BS; §7 cost.
- **[ld]** `ldbound_results.md`: §1 design (a‴_k; μ), Lemma GB_w, sizes, Lemma 1.2_w, coprimality (1.1) and
  order of operations; §2 re-checks; §3 Lemma W; §4 Theorem L2_fam; §5 Theorem LD*; §6 Theorem A‴; §6.1
  transfer of lv6.

Notation: L = L_N[1/2, 1] = exp((log N log log N)^{1/2}); L[1/2, c] = L^c. "L^{o(1)}" means exp(o((log N log log
N)^{1/2})) uniformly in N; fixed powers L^ε and N^{o(1)} are not L^{o(1)}. ϱ(x) = ϱ(x, B) = Ψ(x, B)/x.

---

## 1. Parameters

Fixed (independent of N):
- an odd degree **d ≥ 3**;
- an integer **k > (d + 1)/4**, and ν := 1/(4k);
- a rational base exponent **μ** with ν + η₀ ≤ μ < 1;
- a slack **η₀ > 0**.

Derived exponents:
- **h-floor**: h_*(μ) := max(2ν [Kern], μ − (d − 2)ν [FM level ½], 1 − dμ + ν [window]); h := h_*(μ) + η₀.
- **θ** := (μ − ν) + (μ + h + dν) = 2μ + h + (d − 1)ν. The two summands are the exponents of the rational
  cofactor G(P)/ℓ and of the algebraic value F_A(P).
- **a** := θ/√2. Require **θ < 1** (θ ≤ 1 suffices for every size side-condition; θ < 1 for the final
  inequality 2a < (1 + θ)/√2).

N-dependent quantities (N is the integer after preprocessing, §3):
- **B** := ⌊L_N[1/2, a]⌋, the smoothness bound and trial-division bound.
- **X** := ⌊√(N/2)⌋. Let j be maximal with 2^{kj} ≤ X, Y := 2^{j−1}, and **J** := (Y, 2Y] ∩ primes, so
  ℓ ∈ J means ℓ = N^{2ν+o(1)}, and every product of k primes of J is ≤ X.
- **S**: a computable function of N with S = L^{θ/(4a)+o(1)}; for example S² := ⌈L^{θ/(2a)}⌉·⌈log N⌉^{C₃}
  (§6, item U5). The theory needs only S = L^{σ+o(1)} with σ in a fixed compact set, S² ≥ ρ_r^{−2}L^{c}
  for some fixed c > 0, and S² = L^{o(1)}/(ρ_rρ_F); here ρ_r := ϱ(X_r), ρ_F := ϱ(X_F) with X_r, X_F below.
  The roadmap's quoted form S² = (log N)^{C₃}/(ρ_rρ_F) ([pcest] Thm 4.3, [ld] §5) is one such choice.
- **V_ℓ** := ⌊√ℓ·S⌋ for ℓ ∈ J; envelopes V₋ := ⌊√Y·S⌋ ≤ V_ℓ ≤ V_max := ⌊√(2Y)·S⌋.
- **M** := ⌊N^μ⌋; **H** := ⌊N^h⌋; **W** := ⌈N^{μ+ν+η₀/2}⌉.
- **Sizes** ([ld] §1 "Sizes"): X_r := 2MV/ℓ = N^{μ−ν}L^{O(1)}; X_F ≍ X* := C′_d·M·H·V_max^d =
  N^{μ+h+dν}L^{O(1)}.
- **Row target** M_rows := C_col + 12, C_col := 2 + π(B) + #{(r, s): r ≤ B prime, r ∤ D, f(s) ≡ 0 mod r} + t_χ
  ≤ (d + 1)π(B) + 2 + t_χ, with t_χ := ⌈5d(ln|H|_bound + 2 ln log N)⌉ = O_d(log N) ([lv6] §2, §3.3).

**w distribution.** w is uniform among the primes in (W, 2W]. Since W > 2M·V_max for large N, w exceeds
every rational modulus t = |G(P)| ≤ 2MV, every ℓ ∈ J, B and M ≥ g₁. Since μ + ν + η₀/2 < 1 (θ ≤ 1 gives
μ + ν + η₀/2 ≤ ½), w < N.

**Field law 𝒜‴** ([ld] §1). Independent: (g₀, g₁) uniform on 𝒢 := {(g₀, g₁) ∈ [M/2, M]² ∩ Z² : gcd = 1};
w as above; c uniform on C_H := [−H, H]^d ∩ Z^d. Then
- F̂_{g,w} := Lemma GB applied to the integer N·w, so F̂(g₀, g₁) = Nw ([ld] Lemma GB_w);
- G(x, y) := g₁x − g₀y; R_c(x, y) := Σ_{i<d} c_i x^{d−1−i}y^i; F_A := F̂_{g,w} + G·R_c; f(x) := F_A(x, 1).
- So F_A(g₀, g₁) = Nw, f_d = a_d + g₁c₀, and ψ: x ↦ g₀g₁^{−1} (mod N) kills f.

**Candidates** ([pcest] §1, [ld] §1). For ℓ ∈ J:
- 𝒫_ℓ := {P = (a, b) : b prime, V_ℓ/2 < b ≤ V_ℓ, |a| ≤ V_ℓ, gcd(a, b) = 1};
- Λ_ℓ(g) := {P : ℓ | G(P)};
- the **filtered count** N_ℓ(A) := #{P ∈ 𝒫_ℓ ∩ Λ_ℓ(g) : MV_ℓ/32 ≤ |G(P)| ≤ 2MV_ℓ, G(P)/ℓ B-smooth,
  F_A(P) ≠ 0 B-smooth}, and N_ℓ(A) := 0 if ℓ | N;
- **G_A** := {ℓ ∈ J : N_ℓ(A) ≥ 1}, g_A := |G_A|/π_J, γ := E_{𝒜‴} g_A.

**The regime in which c < 1.** c := (1 + θ)/√2 < 1 ⟺ (1 + θ)² < 2 ⟺ θ < √2 − 1 = 0.41421….
Minimising θ over μ (exactly, `src/thmA/params.py`, independent of `src/ldbound/params.py` and agreeing with
it on every tested (d, k)):
- for k < (d² + 1)/4: μ_opt = (1 − ν)/d, θ = 2/d + ((d² + d − 2)/d)ν (Kern and window floors bind);
- for k ≥ (d² + 1)/4: μ_opt = (1 + (d − 1)ν)/(d + 1), θ = 3/(d + 1) + ((4d − 2)/(d + 1))ν (level and window
  floors bind).

| d | min k with c < 1 (η₀ ↓ 0) | θ there | c there | η₀ must be < | c at k = 50 | k = 1000 | inf over k, η₀ |
|---|---|---|---|---|---|---|---|
| 7 | 21 | 139/336 | 0.99963 | 0.00052 | 0.98376 | 0.97285 | 0.97227 |
| 9 | **13** | **16/39** (μ = 17/156, h = 1/26 + η₀) | **0.99720** | **0.00396** | 0.93126 | 0.91984 | **0.91924** |
| 11 | 13 | 9/22 | 0.99638 | 0.00512 | 0.89626 | 0.88450 | 0.88388 |
| 19 | 17 | 257/646 | 0.98842 | 0.01638 | 0.85188 | 0.81383 | 0.81317 |

At d = 9, k = 12 the minimum is θ = 23/54 = 0.4259 > √2 − 1, so k ≥ 13 is necessary and sufficient
((1 + 16/39)² = 3025/1521 < 2 = 3042/1521). The infimum (1 + 3/(d + 1))/√2 over (k, η₀) is not attained.

## 2. Statement

**Theorem A (design (a‴_k; μ)) [claimed; assembled from the sources].** Fix (d, k, μ, η₀) as in §1 with θ < 1.
There is a randomized Las Vegas algorithm that, on input any integer N ≥ 2, outputs the complete prime
factorisation of N, and whose expected running time is at most
  L_N[1/2, (1 + θ)/√2 + o(1)],
where the o(1) → 0 as N → ∞ at a rate depending only on (d, k, μ, η₀).

Remarks on the statement.
- **(S1) Infimum, not o(1).** For each ε > 0 there are fixed parameters with (1 + θ)/√2 < 0.91924 + ε at d = 9
  (and < 1/√2 + ε as d grows). The form "L[1/2, (1 + 3/(d + 1))/√2 + o(1)] at fixed d" would require letting
  k = k(N) → ∞, and the sources do not carry out that diagonalisation: their N₀ and polylog exponents depend on
  k, and k enters as (log N)^{O(k)} and C(2k, k). At d = 9 the smallest admissible fixed k gives c ≈ 0.9972.
- **(S2)** Everything is unconditional apart from the cited published theorems (§7). The internal results are
  unrefereed.

## 3. The algorithm

All primality tests are AKS (deterministic polynomial time). "gcd check x" means: compute g = gcd(x, N); if
1 < g < N output g (a proper factor; recursion continues on g and N/g); if g = N reject the current object.

```
FACTOR(N):
  if N = 1: return []
  if N is prime (AKS): return [N]
  if N = m^e for some e ≥ 2 (integer roots, e ≤ log2 N): return FACTOR(m) repeated e times
  trial-divide N by all primes ≤ B(N); recurse on each prime found; let N' be the cofactor
  if N' ∈ {1} or N' prime or N' a perfect power: handle as above
  if N' < N₀(d,k,μ,η₀): factor N' by trial division (O(1) time)
  g := SPLIT(N')           # 1 < g < N', g | N'
  return FACTOR(g) ++ FACTOR(N'/g)

SPLIT(N):    # N > N₀, odd, every prime factor > B, not a prime power (so ω(N) ≥ 2)
  # Phase F: field selection ([hypD] Corollary K, "retain the first field whose attempt succeeds")
  loop:
    A := SAMPLE_FIELD(N)                       # may exit with a proper factor
    if A = REJECT: continue
    row := DESCEND_ONCE(A, N)                  # may exit with a proper factor
    if row ≠ FAIL: keep A; rows := [row]; break
  # Phase R + L
  loop:
    while |rows| < M_rows: row := DESCEND_ONCE(A, N); if row ≠ FAIL: append row
    (S, Y) := LV6(A, rows)                     # [lv6] Algorithm LV6 steps 1–6, Las Vegas
    X := ∏_{i ∈ S} r_i mod N
    g := gcd(X − Y, N); if 1 < g < N: return g
    rows := []                                 # trivial gcd: fresh hidden roots, same field ([lv6] §3.7)

SAMPLE_FIELD(N):
  draw (g0, g1) uniform in [M/2, M]^2 until gcd(g0, g1) = 1
  draw w uniform among primes in (W, 2W] (uniform integers + AKS)
  draw c uniform in [−H, H]^d
  gcd checks: g1, w, then f_d = a_d + g1 c0 (0 < |f_d| < N)
  f := F̂_{g,w}(x, 1) + (g1 x − g0) R_c(x, 1)        # Lemma GB on N·w
  e := content(f);  e | N w;  if 1 < gcd(e, N) < N: output it;  if e = w: return REJECT   # e ∈ {1, w} otherwise
  if no prime r ∈ (r0(d), ½ log H] gives an irreducible degree-d reduction f mod r: return REJECT   # Lemma IR
  m := g0 g1^{-1} mod N;  gcd check f'(m) (with the factor t^{d−2}); if N | f'(m): return REJECT
  return A = (g0, g1, w, c, f)

DESCEND_ONCE(A, N):
  draw r uniform in [1, N − 1]; gcd check r            # a non-unit r is a factor (gcd = N impossible)
  R := r^2 mod N
  (u0, v0) := shortest nonzero vector of {(u, v) ∈ Z^2 : u ≡ R v (mod N)} (Gauss reduction)
  factor |u0| and |v0| completely with Lenstra–Pomerance (1992)          # the shape test, T_shape
  certify every output prime by AKS                                       # Las Vegas (complexity critic C1)
  gcd checks u0, v0
  require: |u0|, |v0| each a product of at most k DISTINCT primes, all in J;  else return FAIL
  for each prime ℓ | u0 v0:  P_ℓ := SQ_SEARCH(A, ℓ)  (cached per (A, ℓ));  if FAIL: return FAIL
  return row (r, R, sign(u0 v0), {ℓ: ±1}, {P_ℓ with factorisations of G(P_ℓ)/ℓ and F_A(P_ℓ)})

SQ_SEARCH(A, ℓ):           # deterministic function of (A, ℓ); [pcest] §1, Lemma 1.3; [lv6] Prop BS
  gcd check ℓ
  Gauss-reduce a basis of Λ_ℓ(g); if λ1 < √ℓ/(2S): return FAIL      # irregular: no candidates (Lemma 1.3)
  enumerate Λ_ℓ(g) ∩ [−V, V]^2 (≤ C_abs S^2 points)
  keep P ∈ 𝒫_ℓ (b prime by AKS), G(P) ≠ 0, F_A(P) ≠ 0, MV/32 ≤ |G(P)| ≤ 2MV; gcd checks G(P), F_A(P)
  batch-compute the B-smooth parts of all |G(P)/ℓ| and |F_A(P)| (Bernstein product/remainder tree)
  return the first P in a fixed order with both values B-smooth (trial-divide its two values), else FAIL

LV6(A, rows):              # [lv6] §2 Algorithm LV6, transferred by [ld] §6.1
  columns: sign; v_p mod 2 for p ≤ B; ONE parity column E_i = Σ_P e_iP mod 2;
           (r, s) for r ≤ B, r ∤ D := f_d·disc f, s a root of f mod r; t_χ quadratic characters at
           uniform first-degree primes of norm in (x/2, x], x := the least integer with log(x/2) ≥ max(20d (log Δ̄)^2, (log N)^2), so log x = O_d((log N)^2)
  two uniform left-kernel vectors s1, s2 (Wiedemann / Kaltofen–Saunders, verified; retry until independent)
  for S in (s1, s2, s1+s2): Γ_S := g'(β)^2 ∏_P z_P^{|e_P(S)|};  η := SQ(Γ_S) (Hensel at the least odd prime p ∈ (B, 2B] with p ∤ D);
                            if η found: break      (if none: resample characters and kernel vectors)
  Y := √q_S · t^{−E(S)/2} · ψ(η) · (ψ(γ⁻_S) ψ(g'(β)))^{−1} mod N
  return (S, Y)
```

Row algebra ([lv6] §1, with the sign made explicit). With t := f_d, β := tα, z_P := ta_P − b_Pβ ∈ Z[β], and
G(P_ℓ) = σ_ℓ ℓ n_ℓ (σ_ℓ = ±1, n_ℓ > 0 B-smooth), ψ(z_P) = tG(P)/g₁, hence ℓ ≡ σ_ℓ g₁ t^{−1}ψ(z_P)/n_ℓ (mod N).
For a row with R ≡ u₀/v₀ and e_ℓ = +1 for ℓ | u₀, −1 for ℓ | v₀:
  r² ≡ q · t^{−E} · ∏_ℓ ψ(z_{P_ℓ})^{e_ℓ} (mod N),  q := sign(u₀v₀)·∏_ℓ (σ_ℓ g₁/n_ℓ)^{e_ℓ} ∈ Q^×,  E := Σ e_ℓ.
(The factor sign(u₀v₀) is implicit in [lv6] (1.1), which writes u, v as positive products; see §6, U6.)

## 4. Proof roadmap

Each step lists: what is used, where it is proved, the hypotheses as they must hold in (a‴_k; μ), and what the
next step consumes. "Proved in X" means proved in that internal note.

**Step 0 — preprocessing and the size of ω.**
- FACTOR removes prime factors ≤ B, perfect powers and primes. The surviving N has every prime factor > B
  and ω(N) ≥ 2. Then ω(N) ≤ log N/log B, so 2^{ω} ≤ exp(O((log N/log log N)^{1/2})) = L^{o(1)}
  (cf. [hypD] Lemma Ω; [pcest] §1 "Preprocessing").
- Cost L^{a+o(1)}. Consumed by Theorem S′ (the factor 2^{−ω}), by Lemma 8.5′/FM ((N, t) = 1), by
  coprimality (1.1), and by DD (μ₂(Z/N) ≅ {±1}^ω for odd N, repeated prime factors allowed).

**Step 1 — field law and algebraic identities.** [ld] §1: Lemma GB_w (F̂(g₀, g₁) = Nw, |a_i| ≤ g₁/2 for
i ≥ 1, |a₀| ≤ 2^d Nw/M^d + d2^{d−1}M, a_d ≢ 0 mod g₁), the sizes paragraph (|F̂(P)| ≤ C_d MHV^d N^{−η₀/2} by
the window floor; all values in [−X*, X*]), Lemma 1.2_w (g₁^d F_A(P) ≡ Nw b^d mod G(P)), and the
coprimality statement (1.1): **for ℓ ∤ N and P ∈ 𝒫_ℓ ∩ Λ_ℓ(g) with G(P)/ℓ B-smooth, t := |G(P)| satisfies
(t, bg₁Nw) = 1 and F̂_{g,w}(P) ≡ c_P·w (mod t) with c_P := Nb^dg₁^{−d} ∈ (Z/t)^× independent of w.**
Consumed by FM Step 2(a) (reduced moving residue) and by LD* Step 4 (Lemma W).

**Step 2 — the special-q search recognises exactly N_ℓ(A).** [pcest] Lemma 1.3 (an irregular lattice,
λ₁ < √ℓ/(2S), has no candidate; a regular one has ≤ C_abs S² box points), [lv6] Lemma BS / Prop BS (batch
smooth parts per complete box, cost Õ(B + S² log N), deterministic in (A, ℓ)), [ld] §6.1 (row "Prop BS":
unchanged). Hence SQ_SEARCH succeeds iff N_ℓ(A) ≥ 1, so the set G_A defined in §1 is the set of descendable
primes, as Theorem K requires. Cost T_SQ := Õ(B + S² log N) = L^{θ/(2a)+o(1)} (since B ≤ S² when θ ≤ 1).

**Step 3 — first moment FM.** [pcest] §2.3 Theorem FM, with [pcest] §2.2 Lemma 8.5′, re-checked for
(a‴_k; μ) in [ld] §2 (proof applied for each fixed w, uniformly in w, then averaged).
- Conclusion consumed: **E_{ℓ∈J, A∼𝒜‴} N_ℓ(A) ≥ S²ρ_rρ_F/(log N)^{O(1)}**, averaged over ℓ as well as A. The
  bad-modulus step does not give a lower bound for each individual ℓ, and none is needed (codex thmA_1 §5,
  thmA_2 §8).
- Tenenbaum's hypothesis is ω(t) ≤ y^{c/log(1+u)}, not y^c/log(1+u) as quoted in [pcest] §2.1. It is applied to
  the B-smooth part |n| of t = ℓ|n|, with Ψ_t(x, B) = Ψ_{|n|}(x, B), since ℓ > B. This changes φ(t)/t by a
  factor (1 − 1/ℓ). The range is ample (codex thmA_1 §2; Opus analytic critic F3).
- Hypotheses in this design: x = X_F with exponent ξ = μ + h + dν; moduli t ∈ [MV/32, 2MV] with exponent
  q₀ = μ + ν; (i) q₀ > 2ν + η₀/2 ⟸ μ ≥ ν + η₀; (ii) q₀ < ξ/2 ⟺ h > μ − (d − 2)ν (level floor);
  (iii) Kern/LV-kernel: h > ν and H > dωV with ω = L^{C₁} (Kern floor); (iv) |F̂|/(tY_ω) ≪ ωN^{−η₀/2}
  (window floor); (v) the moving residue is reduced (Step 1); (vi) M/V = N^{μ−ν+o(1)} dominates the Möbius
  error 2^{ω(ℓn)} in Step 1(iii) of FM.
- Structure (re-derived): Σ_{P∈𝒫_ℓ} ≈ V²/log V candidates; for each, the rational fibre count gives weight
  ≈ ρ_r/ℓ (Step 1 of FM: |𝒢_P(n)| ≍ M/V on an n-interval of length ≍ MV/ℓ); for a good modulus t the
  algebraic progression F̂ + G·⟨c, v(P)⟩ has smooth density ≥ ρ_F/(log N)^{O(1)} (Step 2, via LV Lemmas 5.13,
  5.18, 5.20); bad moduli carry an o(1/log log N) share by Lemma 8.5′ (Step 3). Product:
  (V²/ℓ)ρ_rρ_F/(log N)^{O(1)} = S²ρ_rρ_F/(log N)^{O(1)}.

**Step 4 — second moment LD*.** [ld] §5 Theorem LD*, using [pcest] §4.1 Lemma Kern (i), (iii) (with the
effective wedge bound ‖v ∧ v′‖ ≥ |Δ|(bb′)^{d−2}), [ld] §3 Lemma W and Corollary W (Brun–Titchmarsh:
Pr_w[w ≡ a mod t] ≤ (24/η₀)/φ(t) for t ≤ 2MV_max), [ld] §4 Theorem L2_fam, and [pcest] §3 Lemma DS.
- Conclusion consumed: **E_{ℓ,𝒜‴}[N(N − 1)] ≤ (log N)^{O(1)}S⁴ρ_r²ρ_F²**, uniformly for every S = L^{σ+o(1)}
  with σ in a fixed compact set.
  - [ld] §5 states its hypothesis with S² = (log N)^{C₃}/(ρ_rρ_F), but Steps 0–7 use S only through
    B^{C₀} ≥ 64S²L/(ϱϱ), 2S² < Y and S < √ℓ.
  - [ld] Step 0 should read "≤", because of the filters.
- Chain of LD* (re-derived scaling): Pr_c[both F_A(P), F_A(P′) smooth] ≤ C_d(H²V^{2d−4})^{−1}·Z (Kern) →
  Z ≤ ½(Q(P) + Q(P′)) (AM–GM; κ a unit mod Δ) → Q(P) ≤ 𝔔₀(t, |Δ|; r₀), a function of r₀ = c_Pw mod t only
  (window enlargement to [−X*, X*], bounded scale since tT ≍ X*) → E_w 𝔔₀ ≤ (C_W/φ(t))𝔔(t|Δ|) (Corollary W) →
  the algebraic factor depends on g only through n = G(P)/ℓ, so the partner's rational indicator is summed
  along the fibre by Lemma DS: ≤ (M/V)ρ_r(log N)^{O(1)} → ≤ 12 partners per (P, q) → Theorem L2_fam (4.1):
  E_ℓΣ_n 𝔔(ℓ²|n|q)/φ(ℓ|n|) ≤ τ(q)(log N)^{O(1)}Ψ(X′₊)Ψ(X*)²/(Yq(MV)²) → Σ_q τ(q)/q ≪ log² S. Powers:
  (H²V^{2d−4})^{−1}·M^{−2}·V²·(M/V)ρ_r·(MV/Y)ρ_r·(MHV^d)²ρ_F²/(Y(MV)²) = (V⁴/Y²)ρ_r²ρ_F² ≍ S⁴ρ_r²ρ_F².
- Hypotheses in this design: Kern needs H ≥ V² (Kern floor); L2_fam needs R₀ := B^{C₀} ≥ 64S²L/(ϱ(X′₊)ϱ(X*))
  and h > (5 − d)ν (implied by the Kern floor); Lemma W needs log(W/t) ≥ (η₀/3)log N (W ≥ N^{η₀/2}·2MV_max)
  and log W ≤ 2 log N.

**Step 5 — Paley–Zygmund ⇒ D_k.** [hypD] §4 Prop PZ, in the filtered form of [pcest] Thm 4.3 (with LD* in
place of LD + Theorem R + Prop 4.2, as [ld] §5 and §6 state):
  γ = Pr_{A,ℓ}[N_ℓ(A) ≥ 1] ≥ (E N)²/(E[N(N−1)] + E N).
Put m := S²ρ_rρ_F. With Steps 3–4, E N ≥ m/P and E[N(N−1)] ≤ Pm² for a fixed polylogarithm P, so
γ ≥ 1/(P³ + P/m). For the implicit choice m = (log N)^{C₃} this is ≥ 1/(log N)^{O(1)}. For the computable choice
of §1, m = L^{o(1)}, and m is in fact smaller than every inverse power of log N because of de Bruijn's
second-order term. So only **γ ≥ L^{−o(1)}** holds there, and it suffices (codex thmA_1 §§14–15). Consumed: **D_k(γ): γ ≥ L^{−o(1)} uniformly in N.** Theorem R ([pcest] §3) is not on
the critical path.

**Step 6 — rejections and factor exits.** [lv6] §4 (Lemma IR, N | f′(m)), [ld] §1 "Order of operations" and
§6.1 (content e = w; gcd(f_d, N) is a factor exit since 0 < |f_d| < N; reducible f rejected via the Lemma IR
certificate).
- The moments are proved under the raw product law 𝒜‴. Rejection event 𝓡 = {no Lemma-IR certificate} ∪
  {N | f′(m)} ∪ {e = w} has Pr(𝓡) ≤ exp(−log H/(13d log log H)) + O(1/H + 1/W) = L^{−ω(1)} = o(γ).
- Hence E[g_A·1_{𝓡ᶜ}] ≥ γ − Pr(𝓡) ≥ γ/2 under the raw law.
- **Factor exits (repair, codex thmA_2 §5).** This inequality does not by itself survive conditioning on fields
  that neither reject nor produce a factor, and a special-q search may exit with a factor on an earlier
  candidate. Absorbing repair: for one raw field trial, let δ := Pr[factor exit] and b := E q_A := Pr[row
  delivered], and p := δ + b.
  - The counting of Step 7, applied to "row or factor exit", gives p ≥ L^{−o(1)}.
  - A retained row-field is reached with probability b/p, and its law is ∝ q_A. So its contribution to the
    expected cost of C rows is ≤ (b/p)·C·C_trial/b = C·C_trial/p.
  - Selection costs ≤ C_trial/p. Further factor exits only shorten the run. The exponent is unchanged.

**Step 7 — bounded-k splits: from D_k to rows.** [hypD] §2 Theorems S, S′ and §3 Theorem K(i)–(iii), field
selection of Corollary K.
- K(i): Pr_A[g_A ≥ γ/2] ≥ γ/2.
- K(ii) via S′ (needs |G_A| ≥ 4k², G_A ⊂ J coprime to N, R uniform on Q_N, products ≤ X with 2X² ≤ N):
  Pr_R[R ≡ u/v, u, v ∈ S_A] ≥ c_k(g_A/(4 log N))^{2k}2^{−ω}, c_k > 0.
- K(iii): accept iff the primitive vector (u₀, v₀) has both sides products of ≤ k distinct primes of G_A
  (divisor-closed). Pr[right shape]/Pr[accept] ≤ (log N)^{O(k)}g_A^{−2k}2^{2ω} = L^{o(1)}.
- Field selection: per-attempt success q_A; μ_sel := E q_A ≥ (γ/4)·c_k(γ/(16 log N))^{2k}2^{−ω} = L^{−o(1)}; the
  kept field has law ∝ q_A, so E[1/q_A | kept] = Pr(q_A > 0)/μ_sel ≤ 1/μ_sel ([hypD] Cor. K, item 1).
- Consumed: in the kept field, the expected cost of one row is
  (1/q_A)(T_shape + Pr[shape]·2k·T_SQ) with E[1/q_A] ≤ L^{o(1)}, i.e. (T_shape + T_SQ)·L^{o(1)},
  where T_shape = L_{√N}[1/2, 1 + o(1)] = L_N[1/2, 1/√2 + o(1)] (Lenstra–Pomerance 1992, complete factorisation
  of |u₀|, |v₀| ≤ √(2N)).
- The rows are i.i.d. given the field. Accepted roots are not uniform on (Z/N)^× (acceptance filters R), but
  conditional on all certificates each r_i is uniform on its fibre r_i·μ₂, independently. That is exactly what
  DD uses ([lv6] §2's "uniform on (Z/N)^×" hypothesis should read this way; codex thmA_2 §13).
- [hypD] K(ii) as displayed uses the fourth-power form of Theorem S; the 2k exponent above is the Theorem S′
  form. In K(iii)'s shape-count, replace S_J by T_J := ∪_{j≤k}{products of j distinct J-primes}, with
  |T_J| ≤ C_kπ_J^k (codex thmA_2 §4).

**Step 8 — LV6: a congruence of squares from M_rows rows.** [lv6] §2 Prop LV6 with §3.1–3.6 (Lemma C: good
primes; Lemma H: |H| ≤ 2^{dω(D)+d}h_K, log|H| = O_d(log N); Lemma E: ≤ 1 exceptional class per field (Stark
Lemma 8 + Landau); Lemma X: characters separate H ∖ (E ∪ {0}) w.p. ≥ 1 − (log N)^{−2} (Lagarias–Odlyzko);
§3.4: three tests s₁, s₂, s₁ + s₂ suffice; Lemma SQ: exact square root in Z[β] via g′(β)O_K ⊂ Z[β]; §3.6: the
congruence), transferred to (a‴_k; μ) in [ld] §6.1 (ψ well defined since N | Nw; degree drop impossible;
content e ∈ {1, w} after factor exits; discriminant, index and class-number bounds N^{O(1)} because
‖f‖ ≪ MH < N; nothing requires (w, f_d disc f [O_K:Z[β]]) = 1).
- Consumed: a Las Vegas procedure, expected time (M_rows)²(log N)^{O_d(1)} + L^{a+o(1)}, returning S ≠ ∅ and Y
  with X_S² ≡ Y² (mod N), where S and Y are functions of the field, the certificates and fresh coins only.

**Step 9 — fruitfulness: Theorem DD and (T).** [DD] Theorem DD; [lv6] §3.7.
- (T): every exposed quantity (field (g₀, g₁, w, c), the R_i, splits, shape-test outputs, special-q atoms,
  columns, characters, kernel vectors, Hensel prime, sign branch, η, the choice among s₁, s₂, s₁ + s₂, Y) is a
  function of R_i = r_i² and coins independent of the r_i. The hidden roots enter only X_S.
- Then X_S/Y = ε_{i₀}·Z with ε_{i₀} uniform on μ₂ ≅ {±1}^ω independent of everything else, so
  Pr[1 < gcd(X_S − Y, N) < N | 𝒯] = 1 − 2^{1−ω} ≥ 1/2.
- After a trivial gcd, fresh rows are drawn (no adaptive reuse of the same rows), so the number of LV6 rounds is
  geometric with mean ≤ 2, independently of the field.

**Step 10 — cost.** [lv6] §7 "Cost assembly", [ld] §6.
- Rows: M_rows = (d + 1)π(B) + O(log N) = B·L^{o(1)}.
- Per row: (T_shape + T_SQ)·L^{o(1)} with T_shape = L^{1/√2+o(1)} and T_SQ = Õ(B + S²) = L^{θ/(2a)+o(1)}.
- Field selection L^{o(1)} attempts, each T_shape + 2k·T_SQ + polylog (corrected after codex thmA_2 §17); LV6 L^{2a+o(1)}; preprocessing L^{a+o(1)}.
- Total: L^{a + max(1/√2, θ/(2a)) + o(1)} + L^{2a + o(1)}. Minimised at a = θ/√2: both relation terms are
  L^{(1+θ)/√2}, and 2a = √2θ < (1 + θ)/√2 iff θ < 1. For a < θ/√2 the term a + θ/(2a) is larger (it decreases
  on (0, √(θ/2)] and θ/√2 < √(θ/2) for θ < 1); for a > θ/√2 the term a + 1/√2 is larger.
- Expected time: Wald's identity over the geometric number of field attempts, of rows per field
  (E[1/q_A | kept] ≤ 1/μ_sel), of LV6 rounds, and the ≤ log₂N recursive SPLIT calls on factors of N (each
  ≤ L_N[1/2, c + o(1)]).

## 5. Exponent ledger (re-derived for this document)

| quantity | exponent of N (zero slack) | source |
|---|---|---|
| special-q prime ℓ | 2ν | §1, J |
| box radius V | ν (+ L^{O(1)}) | V = √ℓ S |
| rational modulus t = \|G(P)\| | μ + ν | \|G\| ∈ [MV/32, 2MV] |
| rational cofactor G/ℓ | μ − ν | X_r |
| constant coefficient a₀ | max(μ, 1 + μ + ν − dμ) | Lemma GB_w with w ≈ N^{μ+ν} |
| algebraic value F_A(P) | μ + h + dν | X* ≍ MHV^d, needs a₀V^d ≤ MHV^d·N^{−η₀/2} ⟺ window floor |
| pair modulus t\|Δ\| | ≤ μ + 3ν | \|Δ\| = ℓq, q ≤ 2S² |
| θ | 2μ + h + (d − 1)ν | sum of the two smooth sizes |
| log(1/ρ_r) + log(1/ρ_F) | (θ/(2a))·(log N log log N)^{1/2}(1 + o(1)) | u log u with u = ξ log N/log B |

## 6. Uniformity and bookkeeping (precisions supplied by this document)

- **U1 (small factors).** Trial division to B(N) is part of FACTOR. Every later estimate is applied to the
  cofactor N′ with its own B(N′) ≤ B(N); N′ has no prime factor ≤ B(N) ≥ B(N′). If N′ < N₀ it is factored in
  O(1) time. L_{N′} ≤ L_N.
- **U2 (prime powers, perfect powers).** Removed before SPLIT; N′ = p^e is detected by integer roots plus AKS.
  Non-squarefree N′ with ω ≥ 2 (e.g. p²q) needs no special treatment: DD uses |μ₂| = 2^ω, and FM/LD*/LV6
  use only (N, t) = 1 and N | F(g₀, g₁).
- **U3 (N₀).** Every "for N ≥ N₀" in the sources depends only on (d, k, μ, η₀, the σ-range). The effective
  constants of Stark, Landau and Lagarias–Odlyzko make N₀ effective but huge ([lv6] §7).
- **U4 (o(1) terms).** The L^{o(1)} factors are (log N)^{O(k)}, 2^{O(ω)}, γ^{−O(k)}, C(2k, k), and the Wiedemann
  and Hensel polylogs. All are L^{o(1)} uniformly in N, because k is fixed and γ ≥ L^{−o(1)} (Step 5).
  - With the computable S, γ is not ≥ 1/(log N)^{O(1)} (corrected after codex thmA_1 §15).
  - Pr(𝓡) = exp(−Ω(log N/log log N)) is still o(γ).
- **U5 (S is computable).** Ψ is not computed by the algorithm. Any computable S with S = L^{θ/(4a)+o(1)} works:
  FM gives E N = L^{o(1)} (possibly < 1, by an L^{o(1)} factor), LD*'s bound is uniform over S = L^{σ+o(1)},
  and Prop PZ then gives γ ≥ L^{−o(1)} without E N ≥ 1. The cost T_SQ = S²L^{o(1)} is unchanged in exponent.
- **U6 (sign of the split).** The Gauss-reduced (u₀, v₀) can have u₀v₀ < 0; sign(u₀v₀) belongs in q_i and is
  killed by the sign column.
- **U7 (expected time versus probability).** All loops are Las Vegas with verified outputs (kernel vectors,
  η² = Γ, gcd), so only expected costs need bounding. LP's algorithm is used through its expected-time bound;
  no cap is needed.

## 7. External inputs

Published: Lenstra–Pomerance 1992 (JAMS 5; complete factorisation in expected L_n[1/2, 1 + o(1)]); Lee–Venkatesan
(JNT 187, 2018; arXiv 1805.08873) Lemma 8.5 and Lemmas 5.13, 5.18, 5.20, with Harper's Facts 8.8–8.12
(arXiv 1208.5992); Tenenbaum and Hildebrand on Ψ (LV Facts 3.13, 3.15–3.17, Claim 8.16); Montgomery–Vaughan
1973 (large sieve; Brun–Titchmarsh); Lagarias–Odlyzko 1977 (effective Chebotarev); Stark 1974 Lemma 8; Landau
(via Davenport Ch. 14); Lenstra Bull. AMS 1992 Thm 6.5 (or the elementary fallback of [lv6] §3.2);
Wiedemann 1986 / Kaltofen–Saunders 1991; AKS; Bernstein's smooth parts (elementary version proved in [lv6]
§6.2); standard finite-field factoring and root finding; Gauss's lemma; Buhler–Lenstra–Pomerance 1993 for the
NFS post-processing architecture.

Internal (unrefereed): [DD], [hypD] §§2–4, [pcest] §§1–4.1, [lv6] §§1–4, 6.2, 7, [ld] §§1–6.1.
