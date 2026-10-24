# uprime: Hyp U′ (mixed-smoothness re-proof of Theorem U's moments at the 3/8 parameters) — results

Agent `uprime`, phase 26, 2026-10-24. Draft v1 (before the adversarial round and the Opus critics).

Files:
- This note; pre-registration `notes/uprime_prereg.md` (commit `f5a9bd7`, code frozen there).
- Code: `src/uprime/params.py` (exact leading-order bookkeeping of every margin used here, 27 checks; σ₁, σ₂, the
  constraint (SQ), cost brackets), `hybrid_bv.py` (E1), `toy.py` + `experiments.py` (E2), `analyze.py`,
  `posthoc.py` (E2 post hoc, not pre-registered). Tests: `tests/test_uprime.py`.
- Data: `experiments/raw/uprime_e1.jsonl`, `uprime_e2.jsonl`, summaries `uprime_e{1,2}_summary.json`,
  `uprime_e2_posthoc.json`.
- Codex: `codex/uprime_1.{md,out}` (blind round, launched before any analysis was written down).

Labels: **proved** (from the project's internal, unrefereed chain plus the published results cited — the same status
as Theorem U), **proof plan** (every step named, conditions checked, not written out), **open**, **heuristic**,
**experimentally supported** (toy), **UNVERIFIED**.

Notation as lv5unif §1 and multidescent §3: λ := log log N, D := log N/λ, Λ := D^{3/8}λ, L[3/8, c] = e^{(c+o(1))Λ},
y₁ := B = e^{aΛ}, y₂ := B₂ = e^{b₂D^{1/2}λ}, ρ_i(z) := Ψ(z, y_i)/z.

---

## 0. Verdict

**Plain language.** Multidescent's two-level L[3/8] theorem needed two unproved inputs, PC₃ and U′. U′ asked us to
re-prove Theorem U's moment estimates when the level-1 algebraic value is allowed to be smooth with respect to a
much larger bound B₂ ≈ L[1/2] (its large prime factors are the children) while the rational side stays B-smooth.
**The four re-proof items listed for U′ all go through at α = 3/8**, conditional on the project's internal chain
exactly as Theorem U is. The core item, a "hybrid" restricted Bombieri–Vinogradov estimate with B-smooth modulus
cofactors and B₂-smooth values, follows by re-running the project's proof of Lemma 8.5′ with two smoothness
parameters. Only one step changes, and it produces one new condition, which holds with a margin of D^{1/4}. Two small
design changes are needed: a larger polynomial height (found by the codex blind round) and one inequality between the
constants, a·b₂ > 3d₀/(32κ). **One bullet of U′ is not a re-proof of anything in Theorem U and remains open:** the
statement that a child's box, *conditionally* on its parent having been accepted, still contains about the expected
number of good points. That is a new two-point lower bound (parent point and child point). We state it precisely
(Hyp CM), give a proof plan whose analytic inputs are of the same kind as Lemma 8.5′, and check its level conditions.
It is not proved. Codex, working blind, reached the same verdict, including the same constraint. So Theorem 2L now
rests on PC₃ and CM (which is the first-order part of PC₃), plus the internal chain; U′'s level-1 part is no longer a
hypothesis. No exponent above 3/8 is needed for U′: 3/8 is delivered.

Precisely:
1. **Lemma 8.5H (hybrid restricted BV) [proved].** For moduli t = ℓs with ℓ ∈ (Y, 2Y] prime, ℓ > y₂, and s y₁-smooth,
   Σ_t max_{(r,t)=1}|Ψ(x, y₂; t, r) − Ψ_t(x, y₂)/φ(t)| ≤ C·Ψ(x, y₂)·δ_𝒯·ε₂, δ_𝒯 = ρ₁(Q/Y)/log Y (the density of the
   *actual* modulus family), ε₂ = exp(−cu₂/log²(u₂ + 1)) + y₂^{−1/2}, u₂ = log x/log y₂, under (H1)–(H5) of §2. The one
   new condition is (H4): the small-conductor transfer factor R₁^{1−α(Q′, y₁)}, which was polylogarithmic in
   Lemma 8.5′, is now exp(O((log y₂/log y₁)log log x)); Drappeau's saving exp(−cu₂/log²u₂) must beat it. At 3/8 the
   loss is exp(O(D^{1/8}λ)) and the saving exp(−cD^{3/8}/λ²).
2. **Level 1 with two smoothness parameters [proved].** GM-U′, FM-U′, DS-U′, L2-U′, LD*-U′, D-U′ (§3):
   μ₁ := E N₁ ≥ S₁²ρ₁(X_r)ρ₂(X_F)/(d(log N)^{C}) ≥ (log N)^{−C}, E[N₁²] ≤ (log N)^C(μ₁² + μ₁), γ₁ ≥ (log N)^{−C}.
3. **Mixed scan (Lemma CS′) and child radii [proved].** The first S₁ = 2^j with S²R_B(X_r(S))R_{B₂}(X_F(S)) ≥ 1 gives
   m₁ = S₁²ρ₁ρ₂ ∈ [(log N)^{−2C_R}/4, 64] and S₁² = L[3/8, σ₁], **σ₁ = 3/(8ad₀) + 3d₀/(32κb₂)**. (Lemma CS at a single
   y = B would give an L[1/2]-scale box.) R_{B₂} is computed without enumerating primes up to B₂ (cost B^{1+o(1)}).
   Child boxes: a per-window scan gives S₂(Q₂) with nominal main term m₂ ≥ exp((log log N)²) and
   W₂ = L[3/8, σ₂], σ₂ = (3/(8a))(3/d₀ + d₀b₂) (conservative geometry).
4. **Squarefree large part and discriminant [proved].** Union bound: Pr_c[F(P) B₂-smooth but not acceptable] ≤
   (log N)²/B; subtracting it from GM-U′ needs **(SQ) a·b₂ > 3d₀/(32κ)** (strict). The AP-mod-tp² alternative is not
   needed.
5. **Recursion, search, Lemma IR [proved].** Value factorisations by FACTOR_U (Theorem U) on values ≤ X_F cost
   L_N[7/20]^{O(1)} = L[3/8]^{o(1)} each; Lemmas OS/W apply; the bounded-depth search costs W₁(1 + m̄W₂) =
   L[3/8, σ₁ + σ₂]; Lemma IR's rejection is exp(−Θ(D^{1/2})) = o(γ₂). Cost: c₂L = max(2a, a + s + max(a, σ₁ + σ₂),
   (a + s)/(1 − 2^{−3/8})), s = 2κ(1 + C_{γ₂}).
6. **Hyp CM (conditional child mean) [open; proof plan, §7].** Outside an incidence set of weight ≤ μ₁/log N,
   E[N_𝔮 | E] ≥ m₂(q)/(log N)^{C}. It is a two-point lower bound: a rational pair (G(P)/ℓ, G(P′)) in an AP of modulus
   Δ = det(P, P′) along the parent's fibre, and an algebraic pair whose child value runs over an AP of modulus
   |G(P′)|·|Δ|/q. The plan uses (i) an all-moduli smooth BV (Drappeau Lemme 5 + the large sieve) averaged over the
   child window, (ii) a Lemma-SM-type BV for the product family {s·k : s B-smooth}, (iii) a two-point version of
   one-point smoothing, (iv) multiplicity bounds that may lose exp(O(D^{1/4}λ)), which ε₂-type savings absorb, (v) a
   Markov upgrade. The height floor (H) puts both pairs below level 1/2. Codex's warning that an unconditional error
   ε divided by Pr(E) is exponentially large is answered in §7.3: the plan never divides by Pr(E).
7. **Achievable α.** U′ (items 1–5) holds at 3/8; nothing between 3/8 and 2/5 is needed. CM is equally open on the
   whole interval (codex agrees).
8. **Experiments (pre-registered; toy; 9/9 predictions held).** E1: y₁-smooth-cofactor moduli are not worse than all
   moduli (hybrid sum ratio 0.976–1.0; variance ratio 1.05–1.08; no bad modulus). E2: rational/algebraic smoothness
   multiply (1.041 [1.027, 1.055]); the squarefree/discriminant removal is 0.39%; the conditional child mean exceeds
   the fresh-field null by 1.21 [1.17, 1.25], which post hoc is entirely field quality: against a same-field null the
   ratio is 1.005 [0.975, 1.035].

---

## 1. Setting

**Design 2L′.** Theorem U's design (lv5unif §1) at d := least odd ≥ max(3, d₀D^{1/4}), k := max(⌈κD^{3/8}⌉, d + 1),
ν = 1/(4k), μ := μ_opt = (1 − ν)/d (since k < (d² + 1)/4), η₀ := 1/⌈D^{1/2}⌉, B := ⌊e^{aΛ}⌋, B₂ := ⌊e^{b₂D^{1/2}λ}⌋,
Y = N^{1/(2k)+o(1)} (special-q primes ℓ = L[5/8]), with two changes:
- **(H) height floor**: h := μ + 4ν + η₀ (codex uprime_1: needed so that every child value lies above level 1/2 of
  its modulus; §4.2, §7). Since μ ≍ D^{−1/4} and dν ≍ D^{−1/8}, θ_F := μ + h + dν = (d₀/(4κ))D^{−1/8}(1 + o(1)) is
  unchanged at leading order, and h ≥ h_*(μ) + η₀ (all three of lv5unif's floors hold: 2ν, μ − (d − 2)ν < 0,
  1 − dμ + ν = 2ν).
- **Acc_ℓ(P)** (multidescent §3.2) with one addition: P ∈ 𝒫 ∩ Λ_ℓ(g), |G(P)| ≥ MV₁/32, G(P)/ℓ ∈ 𝒮_B, and
  F(P) = s·q₁⋯q_m with s B-smooth and q_i distinct primes in (B, B₂] dividing none of b·f_d·disc(f)·**N**. (A child
  q | N is a factor exit through gcd(F(P), N); excluding it from Acc only simplifies the moments.)

Sizes (leading order; `params.sizes`): log Y ∼ D^{5/8}λ/(2κ), log X_r ∼ D^{3/4}λ/d₀ (rational cofactor),
log X_F ∼ (d₀/(4κ))D^{7/8}λ (algebraic value), log H ∼ log M ∼ D^{3/4}λ/d₀, u_r := log X_r/log B ∼ D^{3/8}/(ad₀),
u_{F2} := log X_F/log B₂ ∼ (d₀/(4κb₂))D^{3/8}, log u ∼ (3/8)λ. Child cap m̄ := log X_F/log B ∼ (d₀/(4κa))D^{1/2}.

**What fails if one smoothness parameter is used** (multidescent §3.3, codex multidescent_2 #1), confirmed:
- y = B₂ in Lemma 8.5′ controls the family of moduli ℓs with s B₂-smooth, of density ρ₂(Q/Y)/log Y. Restricted to the
  actual family (s B-smooth) the relative error becomes ε·ρ₂(X_r)/ρ₁(X_r) = e^{−cD^{3/8}/λ²}·e^{Θ(D^{3/8}λ)} → ∞.
- y = B for the values: u = log X_F/log B ≍ D^{1/2}, and (R′3) u log(4u log y) ≤ K₀ log y fails (ratio ≍ D^{1/8}).
  Independently, B-smooth values at size X_F cost L[1/2] (that is why the children exist).

## 2. Lemma 8.5H: the hybrid restricted BV [proved]

Notation of paper2 Appendix A: E(z, y; q, a) := Ψ(z, y; q, a) − Ψ_q(z, y)/φ(q); Δ_{z,y}(q) := φ(q)^{−1}Σ_{χ≠χ₀}|Ψ(z, y; χ)|;
H(u) := exp(u/log²(u + 1)); α(z, y) the Hildebrand–Tenenbaum saddle point; δ_D, η_D, c_D, C_D the constants of
Drappeau's Lemme 5 (Fact D), C_dBT that of Fact dBT, C_α that of (α2). E1–E8, Lemma LS and Corollary LS (= lv5unif's
LS′) are the regime-free lemmas of paper2 Appendix A, used as proved there.

**Statement.** Let 2 ≤ y₁ ≤ y₂ ≤ x, and Y, Q > 0. Put u₂ := log x/log y₂ and, for z ≤ x, u₁(z) := log z/log y₁. Let
J_Y := {ℓ prime : Y < ℓ ≤ 2Y}, 𝒯_Q := {t = ℓs ∈ (Q/2, Q] : ℓ ∈ J_Y, P⁺(s) ≤ y₁}, δ_𝒯 := ρ₁(Q/Y)/log Y. Fix K₀ ≥ 1,
K := 2K₀ + 6, and U₀ ≥ 8 such that Hildebrand–Tenenbaum's Theorem 3 error O(1/u + log y/y) is ≤ 1/100 for
u ≥ U₀ − 1. Assume, for x ≥ x₀(K₀):
- (H1) log y₁ ≥ max{4(log log x)², (8K(1 + C_α)/δ_D)·log log x}, and (log x)^{c_D} ≤ y₂;
- (H2) u₂ ≥ max{(log log x)⁴, 2K/η_D};
- (H3) u₂ log(4u₂ log y₂) ≤ K₀ log y₂ and u₁(Q) log(4u₁(Q) log y₁) ≤ K₀ log y₂;
- (H4) K(log y₂/log y₁)(log log x + C_α log log y₁) + log(60C_dBT²C_D log log(x + 16)) ≤ δ_D u₂/(8 log²(u₂ + 1));
- (H5) Y > y₂, Y ≥ y₂^K, and Γ := min{log(Q/Y), log(√x/Q)} ≥ (2K₀ + U₀)log y₂ + 4 log log x.

Then, with C, c > 0 depending only on K₀ and the absolute constants above:
- (i) Qδ_𝒯/10 ≤ #𝒯_Q ≤ 3Qδ_𝒯;
- (ii) Σ_{t∈𝒯_Q} max_{(r,t)=1}|E(x, y₂; t, r)| ≤ C·Ψ(x, y₂)·δ_𝒯·ε₂, ε₂ := exp(−cu₂/log²(u₂ + 1)) + y₂^{−1/2}.

At y₁ = y₂ this is Lemma 8.5′(R′) (with (H4) implied by (R′1)–(R′3)).

**Proof.** We follow lv5unif §5 (Steps A–E), marking every change.

*Step A (saddle points), now at both parameters.* For 2y_i ≤ z ≤ x, (α2) gives 1 − α(z, y_i) ≤
(log(u log y_i) + C_α log log y_i)/log y_i with u log y_i = log z ≤ log x, so by (H1) 1 − α(z, y_i) ≤ 1/(4 log log x)
≤ 1/100 for i = 1, 2 (y₂ ≥ y₁). In particular α(z, y_i) ≥ 3/4. Moreover, by the second part of (H1),
  K(1 − α(z, y₁)) ≤ K(1 + C_α)log log x/log y₁ ≤ δ_D/8.   (A′)

*Step B (density floors), now separately.* E7 at y₂ and z = x: ρ₂(x) ≥ y₂^{−1}exp(−u₂ log(4u₂ log y₂)) ≥ y₂^{−1−K₀}
by (H3). E7 at y₁ and 2y₁ ≤ z ≤ Q: ρ₁(z) ≥ y₁^{−1}exp(−u₁(z)log(4u₁(z)log y₁)) ≥ y₂^{−1−K₀}, because
u ↦ u log(4u log y₁) is increasing, u₁(z) ≤ u₁(Q), and y₁ ≤ y₂. Hence for 2y₁ ≤ Q′ ≤ Q:
  Ψ(x, y₂)ρ₁(Q′) ≥ x·y₂^{−2−2K₀}.   (B′)

*Step C (Lemma SM-H: y₁-smooth moduli, y₂-smooth values, relative saving).* Let 2y₁ ≤ Q′ ≤ √x·y₂^{−(2K₀+U₀)}(log x)^{−4}
and 𝒮₁ := {s ∈ (Q′/2, Q′] : P⁺(s) ≤ y₁}. Claim: Σ_{s∈𝒮₁}Δ_{x,y₂}(s) ≤ C·Ψ(x, y₂)ρ₁(Q′)·ε₂.
By E1 this bounds Σ_s max_a|E(x, y₂; s, a)|. Split the characters by conductor r: (a) 1 < r ≤ R₁ := y₂^K; (b) r > R₁.
- (a) **[the only changed step].** E3 is applied with the value parameter y₂: the Möbius variable d runs over divisors
  of s_{y₂} = ∏_{p|s, p≤y₂}p = rad(s), because s is y₁-smooth and y₁ ≤ y₂; so d is squarefree, y₁-smooth, d ≤ Q′,
  (d, r) = 1. Exchanging sums, the modulus weights are Σ_{s∈𝒮₁, rd|s}1/φ(s), and E8 **with the modulus parameter y₁**
  gives ≤ 10C_dBT log log(Q′ + 16)·ρ₁(Q′)·(rd)^{1−β₁}/(φ(r)φ(d)), β₁ := α(Q′, y₁). (If rd is not y₁-smooth the
  sum is empty and the bound is trivial.) Bound r^{1−β₁} ≤ R₁^{1−β₁}. The d-sum is
  Σ_d d^{1−β₁}φ(d)^{−1}|Ψ(x/d, y₂; χ*)| summed over primitive χ*; by Fact D at (x/d, y₂, Q := R₁, A = 1/2) — its
  hypotheses hold since x/d ≥ √x, (log x)^{c_D} ≤ y₂ by (H1), R₁ = y₂^K ≤ x^{η_D/2} ≤ (x/d)^{η_D} by (H2) — and Fact dBT
  at y₂ (Ψ(x/d, y₂) ≤ C_dBTΨ(x, y₂)d^{−β₂}, β₂ := α(x, y₂)), and E6 (Σ_d d^{1−β₁−β₂}/φ(d) < 6, as β₁, β₂ ≥ 3/4):
    (a) ≤ 60C_dBT²C_D log log(x + 16)·R₁^{1−β₁}·{H(u₂/2)^{−δ_D} + y₂^{−δ_D}}·Ψ(x, y₂)ρ₁(Q′).
  *The transfer factor.* By (α2) at (Q′, y₁), log R₁^{1−β₁} = K log y₂·(1 − β₁) ≤ K(log y₂/log y₁)(log log x +
  C_α log log y₁). In Lemma 8.5′ (y₁ = y₂) this was (log x)^{O(1)}; here it is exp(O((log y₂/log y₁)log log x)).
  By (H4) and δ_D(u₂/2)/log²(u₂/2 + 1) ≥ δ_Du₂/(2 log²(u₂ + 1)), the H-term contributes ≤ exp(−3δ_Du₂/(8 log²(u₂+1))).
  By (A′), R₁^{1−β₁} ≤ y₂^{δ_D/8}, so the y₂-term contributes ≤ y₂^{−3δ_D/4} (the log log prefactor is ≤ y₂^{δ_D/8} for
  x ≥ x₀ by (H1)), and by (H3) (u₂ log u₂ ≤ K₀ log y₂) y₂^{−3δ_D/4} ≤ exp(−3δ_Du₂/(4K₀ log²(u₂ + 1))). So
    (a) ≤ C·Ψ(x, y₂)ρ₁(Q′)·exp(−cu₂/log²(u₂ + 1)), c := δ_D·min(3/8, 3/(4K₀)).
- (b) Unchanged in form: Corollary LS with value parameter y₂, z₀ = x, Q₀ = Q′, R₀ = R₁, 𝒳(s) = {χ mod s : cond χ > R₁}
  for s ∈ 𝒮₁ (any subset of moduli is allowed). Its hypotheses hold since x/Q′ ≥ √x ≥ max(4, y₂^{3/2}). It gives
  (b) ≤ 1.5·10³π(y₂)(log x)²(y₂^{3/2}x/R₁ + y₂^{1/2}x^{1/2}Q′). Relative to (B′) this is
  ≤ 1.5·10³(log x)²(y₂^{9/2+2K₀−K} + y₂^{7/2+2K₀}Q′/√x) ≤ 1.5·10³(log x)²(y₂^{−3/2} + y₂^{7/2−U₀}(log x)^{−4}) ≤ y₂^{−1/2}.

*Step D (count of 𝒯_Q), at y₁.* As lv5unif Step D with y = y₁: ℓ > y₂ ≥ y₁ ≥ P⁺(s), so t determines (ℓ, s);
z := Q/(2ℓ) ≥ e^Γ/4 ≥ y₂^{U₀−1} ≥ y₁^{U₀−1}, α(z, y₁) ≥ 0.99 (Step A), and HT Theorem 3 at y₁ gives (i).

*Step E (the error sum).* E2 (lift averaging) needs ℓ > y₂ (every y₂-smooth n is coprime to ℓ) — this is (H5).
- Part 1 (χ_ℓ principal): Σ_ℓ(ℓ − 1)^{−1}Σ_{s∈(Q/2ℓ,Q/ℓ], P⁺(s)≤y₁}max_a|E(x, y₂; s, a)|. Step C applies with Q′ = Q/ℓ
  (2y₁ ≤ Q/(2Y) by Γ; Q/ℓ ≤ √x e^{−Γ}). Since Ψ(·, y₁) is nondecreasing, ρ₁(Q/ℓ) ≤ 2ρ₁(Q/Y); Σ_{ℓ∈J_Y}1/(ℓ − 1) ≤ 2.6/log Y.
  So Part 1 ≤ 5.2C·Ψ(x, y₂)δ_𝒯ε₂.
- Part 2 (ℓ | cond χ, so cond χ > Y): Corollary LS with z₀ = x, Q₀ = Q, R₀ = Y gives
  ≤ 1.5·10³π(y₂)(log x)²(y₂^{3/2}x/Y + y₂^{1/2}x^{1/2}Q). By (B′) with Q′ = Q/Y, Ψ(x, y₂)δ_𝒯 ≥ x y₂^{−2−2K₀}/log Y, and
  by (H5) Part 2 ≤ Ψ(x, y₂)δ_𝒯·y₂^{−1/2}.
- The principal character of t contributes the main term Ψ_t/φ(t) exactly. ∎

**Corollary bad₂.** Call t *good₂* if for every x in the grid 𝒳 (§3) and every (r, t) = 1,
|Ψ(x, B₂; t, r) − Ψ_t(x, B₂)/φ(t)| ≤ Ψ_t(x, B₂)/(2φ(t)). In the setting of Lemma 8.5H with x ranging over 𝒳 and
y₁ = B, y₂ = B₂, #{t ∈ 𝒯_Q bad₂} ≤ 30C|𝒳|ε₂·#𝒯_Q.
*Proof.* As paper2 Corollary `cor:bad`: Tenenbaum's Ψ_q theorem (Fact Ten) at y = B₂ — its range ω(t) ≤ log N ≤
B₂^{c/log(1+u₂)} holds and its error O(log(1 + u₂)log(1 + ω)/log B₂) + O(ω/B₂) → 0 — gives
Ψ_t(x, B₂)/(2φ(t)) ≥ Ψ(x, B₂)/(3t) ≥ Ψ(x, B₂)/(3Q); Markov and (i). ∎

**Where the design sits** (`params.margins`, all hold for every (a, d₀, κ, b₂)). x ∈ 𝒳 has log x = θ_F log N(1 + o(1)),
u₂ ∼ (d₀/(4κb₂))D^{3/8}; Q ∈ [MV₁/32, 2MV₁], log(Q/Y) ∼ log X_r ∼ D^{3/4}λ/d₀. (H1): log B = aD^{3/8}λ ≫ λ². (H2):
D^{3/8} ≫ λ⁴. (H3): u₂ log u₂ ≍ D^{3/8}λ ≪ log B₂ ≍ D^{1/2}λ, and u_r log u_r ≍ D^{3/8}λ ≪ D^{1/2}λ. **(H4):** LHS
≍ (b₂/a)D^{1/8}λ, RHS ≍ D^{3/8}/λ² — margin D^{1/4}/λ³. (H5): log Y ≍ D^{5/8}λ ≫ log B₂; Γ ≥ min(log X_r,
(1/2)θ_F log N − log Q) ≍ D^{3/4}λ ≫ log B₂ (the level gap is automatic at 3/8: it does not use η₀). The saving
ε₂ = exp(−Θ(D^{3/8}/λ²)) beats every polylogarithm, as FM needs. *General form:* with y₁ = L_x[β₁], y₂ = L_x[β₂]
relative to the value size x, (H4) requires β₂ − β₁ < 1 − β₂; at 3/8, β₁ = 3/7, β₂ = 4/7.

## 3. Level 1 with two smoothness parameters [proved]

The grid is lv5unif's: ϖ := B²S₁², 𝒳 := {X_Fϖ^{−2}(3/2)^i} ∩ [X_Fϖ^{−2}, 16dX_F], X_F := (2d + 1)MHV₁^d, V₁ = ⌊√(2Y)S₁⌋,
X_r := 2MV₁/Y. ρ_r := ρ₁(X_r), ρ_{F2} := ρ₂(X_F).

**Lemma U1₂ (local behaviour at y₂) [proved].** For z₁ with log z₁/log B₂ ≥ D^{1/4} and z₁ ≤ z₂ ≤ z₁B₂^C:
ρ₂(z₂) ≤ 2ρ₂(z₁) and ρ₂(z₁) ≤ (log N)^{C_U(C+1)}ρ₂(z₂). *Proof.* HT Theorem 3 in ≤ C + 1 steps with c ≤ y = B₂;
(1 − α(z, B₂))log B₂ ≤ 2 log log N. ∎ Every shift used below is by a factor ≤ ϖ² = B⁴S₁⁴ ≤ B₂ or ≤ B₂² — one or two
steps — and in fact (1 − α₂)log ϖ² = O(D^{−1/8}λ) → 0, so these shifts cost 1 + o(1) at y₂.

**Lemma GM-U′ [proved].** Hypotheses of lv5unif GM-U, with t good₂. Then Pr_c[F_A(P) ∈ 𝒮_{B₂}] ≥ ρ_{F2}/(log N)^{C_g}.
*Proof.* paper2 `lem:goodmod` / lv5unif GM-U line by line with y = B₂ for the values: the small-layer bound, the
envelope (|F̂(P)| ≤ MHV^dN^{−η₀/4}; Lemma `lem:valuesizes` needs h ≥ 1 − dμ + ν + η₀, true under (H)), the choice of
x₁ ∈ 𝒳 ∩ (ty′/4, ty′/2], goodness₂ and Fact Ten at y = B₂ give Ψ(x₁; t, r′) ≥ Ψ(x₁, B₂)/(3t), and U1₂ gives
ρ₂(x₁) ≥ ρ_{F2}/(log N)^{2C_U}. The two error terms 33eϖN^{−η₀/4} and 3d²V/(2H + 1) are ≤ ρ_{F2}B^{−1}, since
ρ_{F2} = L[3/8]^{−O(1)} while N^{−η₀/4} = L[1/2]^{−Θ(1)} and V/H ≤ N^{−μ}. ∎

**Lemma SQF (local exclusions) [proved].** Let (ℓ, P, g, w) satisfy the hypotheses of GM-U′ (so t = |G(P)| = ℓn,
n B-smooth, gcd(t, bg₁Nw) = 1; b prime > (1 − 1/d)V₁ > B₂), and let c be uniform on C_H.
- (a) For a prime p ∈ (B, B₂]: Pr_c[p² | F_A(P)] ≤ 1/p² + 1/(2H + 1).
- (b) For a prime q ∈ (B, B₂], q ∤ N: Pr_c[q | F_A(P) and q | f_d·disc(f)] ≤ 2e(d − 1)/q².
- (c) For a prime q ∈ (B, B₂], q | N: Pr_c[q | F_A(P)] ≤ 1/q + 1/(2H + 1).
Hence Pr_c[F_A(P) ∈ 𝒮_{B₂} ∖ 𝒜cc] ≤ (log N)²/B for N ≥ N₀, where 𝒜cc is the algebraic part of Acc.

*Proof.* (a) p ∤ G(P) (its prime factors are ≤ B or equal to ℓ > B₂). F_A(P) = F̂(P) + G(P)⟨c, v(P)⟩, so p² | F_A(P)
fixes ⟨c, v(P)⟩ modulo p². The coefficient of c_{d−1} in ⟨c, v(P)⟩ is b^{d−1}, a unit mod p. For fixed c₀, …, c_{d−2},
c_{d−1} ∈ [−H, H] lies in one class mod p²: ≤ (2H + 1)/p² + 1 values. (c) The same with q in place of p².
(b) Since q ≤ B₂ ≤ (2H + 1)/d, every class of C_H mod q has between ⌊(2H+1)/q⌋^d and ⌈(2H+1)/q⌉^d elements, so the law
of c mod q is ≤ (1 + q/(2H + 1))^d ≤ e times uniform. Under the uniform law, F_A mod q is uniform on
𝒜_q := {F ∈ F_q[x, y]_d : F(g₀, g₁) = Nw}: forms vanishing at the point (g₀ : g₁) are exactly G·F_q[x, y]_{d−1}
(G ≢ 0 since gcd(g₀, g₁) = 1), and Nw ≢ 0 (q ∤ N, q ≠ w). P ≢ λ(g₀, g₁) mod q since q ∤ G(P); so F ↦ F(P) is a
nonconstant affine function on 𝒜_q and Pr[F(P) = 0] = 1/q. Given F(P) = 0, F = L_P·H with L_P := bx − ay and H uniform
on the affine hyperplane ℋ := {H ∈ F_q[x, y]_{d−1} : H(g₀, g₁) = κ₁}, κ₁ := Nw/L_P(g₀, g₁) ≠ 0 (L_P(g₀, g₁) = −G(P)).
Now f_d = F(1, 0) = b·H(1, 0), and if f_d ≢ 0 then q | disc(f) iff the binary form F has a repeated linear factor over
F̄_q, i.e. iff L_P | H or H has a repeated factor. So the event is contained in {H(1, 0) = 0} ∪ {H(P) = 0} ∪
{Disc(H) = 0}. The first two are nonconstant affine conditions on ℋ (or empty, when (1 : 0) or P equals (g₀ : g₁)
projectively), probability ≤ 1/q each. Disc(H) is a polynomial of degree 2(d − 2) in the coefficients, not identically
zero on ℋ when q > d (take H = κ∏_{i<d}(x − r_iy) with distinct r_i ∈ F_q ∖ {g₀/g₁} and κ fixing H(g₀, g₁) = κ₁), so
Schwartz–Zippel gives ≤ 2(d − 2)/q. Total ≤ e·(1/q)·(2(d − 1)/q). ∎ (Sum: Σ_{p>B}(1 + 2e(d−1))/p² + (log N/log B)(1/B) +
2π(B₂)/(2H + 1) ≤ (log N)²/B.)

**Theorem FM-U′ [proved].** μ₁ = E_{ℓ,A}N₁(ℓ, A) ≥ S₁²ρ_rρ_{F2}/(d(log N)^{C_FM}), provided (SQ) holds.
*Proof.* lv5unif FM-U / paper2 `thm:FM` with: (1) for each counted (ℓ, P, n, g) with t good₂ and every w,
Pr_c[F_A(P) ∈ 𝒜cc] ≥ ρ_{F2}/(log N)^{C_g} − (log N)²/B ≥ ρ_{F2}/(2(log N)^{C_g}) (GM-U′, Lemma SQF), where the last step
is **(SQ)**: log(1/ρ_{F2}) = u_{F2}log u_{F2}(1 + o(1)) = (3d₀/(32κb₂))Λ(1 + o(1)) < aΛ = log B. (2) Good rational mass
unchanged (B-smooth n, Lemma FIB, Lemma U1 at B). (3) Bad mass: the moduli t = ℓ|n| lie in 𝒯_Q for the 7 dyadic Q
covering [MV₁/32, 2MV₁], with y₁ = B; Corollary bad₂ gives π_J^{−1}Σ_ℓ𝔟_ℓ ≤ C|𝒳|ε₂(MV₁/Y)ρ_r(log N)^C, which is
(log N)^{−C′} times the good mass for every C′, since |𝒳| = O(Λ) and ε₂ = e^{−Θ(D^{3/8}/λ²)}. ∎

**If (SQ) fails** (a·b₂ < 3d₀/(32κ) strictly), the subtraction (log N)²/B exceeds ρ_{F2} by e^{Θ(Λ)} and the union
bound gives nothing; that is a failure of the elementary estimate, not evidence that acceptable values are scarce
(codex). The alternative, B₂-smooth counts in APs mod t·p² for p ∈ (B, B₂], would need moduli t·p² — still in a
hybrid family (p² is y₂-smooth but not y₁-smooth), so it would need a three-parameter version; not pursued, as (SQ)
costs only a constraint on constants.

**Lemma DS-U′ [proved].** For fixed C₀, X with log X/log B₂ ≥ D^{1/4}, 1 ≤ q ≤ B^{C₀}, an interval I ⊂ [XB^{−C₀}, X] of
length A ≥ XB^{−C₀}: #{n ∈ I : n ≡ r (q), P⁺(n) ≤ B₂} ≤ (A/q)ρ₂(X)(log N)^{C}.
*Proof.* paper2 `lem:DS` with divisors d₁ ∈ (D′/B₂, D′] (y = B₂); the dyadic sum has ≤ log₂(2B₂) ≤ log N blocks;
X/(D′/B₂) ≤ B₂B^{2C₀} ≤ B₂², so U1₂ with C = 2 applies. ∎

**Theorem L2-U′ [proved].** paper2 `thm:L2` / lv5unif L2-U with 𝒮* := 𝒮_{B₂} ∩ [−X_F, X_F] ∖ {0} and the rational
family 𝒩_ℓ (B-smooth n) unchanged:
  E_ℓΣ_{n∈𝒩_ℓ}𝔔(ℓ²|n|q)/φ(ℓ|n|) ≤ τ(q)(log N)^{C_L}Ψ(X_r, B)Ψ(X_F, B₂)²/(Yq(MV₁)²), uniformly in 1 ≤ q ≤ 4S₁².
*Proof.* R₀ := B^{C₀} ≥ 64S₁²B/(ρ_rρ_{F2}) with C₀ ≤ 2σ₁/a + 3 fixed. Class (A): N_c ≤ (X_F/m′)ρ_{F2}(log N)^C by DS-U′
(and R₀^{−1} ≤ ρ_{F2}); multiplicities mult(m′) ≤ (Ψ(X_r, B)/m₁)(log N)^{O(C₀)} by Lemma U1 at B (unchanged).
Class (B): large sieve with Σ|a_F|² = #𝒮* ≤ 2Ψ(X_F, B₂); ratios q/(R₀ρ_rρ_{F2}) ≤ 1 and X_rq²/(X_Fρ_rρ_{F2}) =
N^{−h−(d+1)ν}B^{O(1)} = o(1). Class (C): ratios log N/(π_Jρ_{F2}) and X_rq²Y⁴/(π_JX_Fρ_rρ_{F2}) = N^{(5−d)ν−h}B^{O(1)},
both o(1). ∎

**Theorem LD*-U′ [proved].** E_{ℓ,A}[N₁(N₁ − 1)] ≤ d·η₀^{−1}(log N)^{C_LD}S₁⁴ρ_r²ρ_{F2}².
*Proof.* An upper bound may drop the squarefree, discriminant and q ∤ N conditions: Acc ⊂ {rational filter, F(P) ∈ 𝒮_{B₂}}.
Then lv5unif LD*-U Steps 0–7 verbatim (Kern-U with bb′ ≤ V₁² ≤ (2H + 1)/d; W-U; Step 5 with DS-U at B; the pair bound
q ≤ 4S₁²), with L2-U′ in Step 7. The power count V⁴/Y² ≤ 4S₁⁴ and the constant (2d + 1)²/d ≤ 9d are unchanged. ∎

**Theorem D-U′ [proved].** γ₁ := Pr[N₁ ≥ 1] ≥ (log N)^{−C_γ}, and E[N₁²] ≤ (log N)^C(μ₁² + μ₁).
*Proof.* Paley–Zygmund with Lemma CS′ (§4.1): m₁ ≥ (log N)^{−2C_R}/4 and m₁ ≤ μ₁d(log N)^{C_FM}. ∎

These are exactly the level-1 moments Prop UB consumes (μ₁ and E[N₁²]; E[N₂²] ≤ E[N₁²]).

## 4. Box sizes

### 4.1 Lemma CS′ (mixed scan) [proved]

For y ∈ {B, B₂} and z with log z/log y ≥ D^{1/4}, let α̂_y(z) solve Σ_{p≤y}log p/(p^α − 1) = log z and
R_y(z) := z^{α̂−1}∏_{p≤y}(1 − p^{−α̂})^{−1}; R*_y(z) the exact value. As in Lemma CS, ρ_y(z) ≤ R*_y(z) ≤ (log N)^{C_R}ρ_y(z)
(Rankin above; Hildebrand–Tenenbaum Theorem 1 at the saddle point below, φ₂ ≍ log z log y).

*Computation without primes up to B₂.* Enumerating p ≤ B₂ costs L[1/2]. Instead: sum exactly over p ≤ B, and replace
the sums over p ∈ (B, B₂] of log p/(p^α − 1) and −log(1 − p^{−α}) by the integrals ∫_B^{B₂}dt/(t^α − 1) and
∫_B^{B₂}−log(1 − t^{−α})dt/log t. With a de la Vallée Poussin-type explicit error |θ(t) − t| ≤ C₁t exp(−c₂√log t)
(explicit constants exist; our citation of a specific explicit form is UNVERIFIED, see §10), integration by parts bounds
both errors by B₂^{1−α}(log B₂)²exp(−c₂√log B) = exp(O(λ) − c₂D^{3/16}√(aλ)) → 0 uniformly near the saddle point (here
B₂^{1−α} ≤ (log N)^{2}). Hence α̂ is found to accuracy (log N)^{−3} and R_{B₂} to relative accuracy 2 at cost
B^{1+o(1)} (codex uprime_1 suggests the cut at exp((log log N)⁴) instead of B; either works).

**Lemma CS′.** The scan S = 2^j, j = 0, 1, …, stopping at the first S with S²R_B(X_r(S))R_{B₂}(X_F(S)) ≥ 1, stops at some
S₁ with S₁² ≤ e^{(σ₁+o(1))Λ}, σ₁ := 3/(8ad₀) + 3d₀/(32κb₂), and m₁ := S₁²ρ_rρ_{F2} ∈ [(log N)^{−2C_R}/4, 64].
*Proof.* Lemma CS verbatim with the two Rankin factors at different y. Existence: for S² ≤ e^{2σ₁Λ},
d log S = O(D^{5/8}λ) = o(log X_F), so log X_r and log X_F keep their leading terms, and CEP (valid: u ≤ log z/(2 log log z))
gives log(R_BR_{B₂}) = −(u_r log u_r + u_{F2}log u_{F2})(1 + o(1)) = −σ₁Λ(1 + o(1)). Overshoot: X(S/2) ≤ X(S) ≤ 2^{d+1}X(S/2)
≤ BX(S/2), and Lemma U1(i) at B and U1₂(i) at B₂ give m₁ ≤ 64. ∎

*Single-y comparison.* Lemma CS with both factors at B gives log(1/(ρ_B(X_r)ρ_B(X_F))) ∼ (d₀/(8κa))D^{1/2}λ (u ∼ D^{1/2},
log u ∼ λ/2): an L[1/2]-scale box (`params.sigma1_single_y`), as multidescent said.

### 4.2 Child boxes and S₂ [proved]

For a child 𝔮 = (q, r) of P (r = a/b mod q), Λ_𝔮 := {P′ : a′ ≡ rb′ (mod q)}. Two geometries:
- **Conservative (default; codex uprime_1).** Box(𝔮) := primitive points, mod ±1, of {ie₁ + je₂ : |i|, |j| ≤ S₂} for a
  Gauss-reduced basis e₁, e₂ of Λ_𝔮. Since |e₂| ≤ q, every box point has |P′| ≤ V₂ := 2Q₂S₂ for q ∈ (Q₂, 2Q₂], and there
  are ≥ c S₂² of them. Envelopes: X_{r2} := 2MV₂, X_{F2} := (2d + 1)MHV₂^d/Q₂ ≥ |F(P′)|/q.
- *Square (sharper; needs CM's label-law bound).* Λ_𝔮 ∩ [−V₂, V₂]² with V₂ := ⌊√(2Q₂)S₂⌋. Regular lattices
  (λ₁ ≥ √q/(2S₂)) have ≍ S₂² points; irregular ones have ≤ 2 primitive points (Lemma `lem:reg`), so such children are
  dead. The number of residues r mod q with an irregular lattice is ≤ #{primitive v : |v|² < q/(4S₂²)} ≤ q/S₂², i.e. a
  fraction ≤ 1/S₂² of labels; that this is also a negligible fraction *of incidences* needs an upper bound on the
  children's label law (Lemma CL type), which belongs to CM (§7).

**Scan S₂(Q₂)** (conservative geometry): the first 2^j with S²R_B(X_{r2}(S))R_B(X_{F2}(S)) ≥ 4(log N)^{2C_R}·
exp((log log N)²). Then m₂(Q₂) := S₂²ρ₁(X_{r2})ρ₁(X_{F2}) ≥ exp((log log N)²) and S₂² = e^{(σ₂(Q₂)+o(1))Λ} with
  σ₂(Q₂) = (3/(8a))(3/d₀ + d₀b(Q₂)) ≤ σ₂ := (3/(8a))(3/d₀ + d₀b₂),  b(Q₂) := log(2Q₂)/(D^{1/2}λ),
using log X_{r2} ∼ D^{3/4}λ/d₀ and, under (H), log X_{F2} ∼ (2/d₀ + d₀b(Q₂))D^{3/4}λ (square geometry: d₀b/2). One scan
per dyadic window that actually occurs (≤ log₂B₂ of them), each B^{1+o(1)}. The threshold exp((log log N)²) is
super-polylogarithmic but L^{o(1)}: it is chosen so that m₂ beats the unknown polylogarithmic loss of CM whatever its
exponent (the algorithm cannot know C_CM).

**Search sizes.** W₁ := #(Λ_ℓ(g) ∩ [−V₁, V₁]²) ≤ 400S₁² = L[3/8, σ₁]; W₂ := max (2S₂ + 1)² = L[3/8, σ₂].

**What this does and does not give.** m₂ is the *nominal* (product) main term of a child box. That the conditional
mean E[N_𝔮 | E] is ≥ m₂/(log N)^C is Hyp CM (§7) — not a consequence of anything in this section (codex agrees).

## 5. Squarefree large part: constants

Item 3 of multidescent's U′ asked for the constants of the union-bound option. They are in Lemma SQF and FM-U′:
- per prime p ∈ (B, B₂]: Pr_c[p² | F(P)] ≤ p^{−2} + (2H + 1)^{−1} (needs B₂² ≤ (2H + 1)/d; true, log H ≍ D^{3/4}λ);
- discriminant/leading coefficient: ≤ 2e(d − 1)q^{−2} (codex's Schwartz–Zippel on D_A of degree 2d − 1 gives
  (2d − 1)q^{−2}·e, the same up to a constant);
- children dividing N: ≤ (log N/log B)(1/B + 1/(2H + 1));
- total ≤ (log N)²/B; usable iff **a > 3d₀/(32κb₂)**, i.e. log B beats log(1/ρ_{F2}) at leading order. Since σ₁ contains
  3d₀/(32κb₂) as well, (SQ) says: the factor-base exponent a must exceed the algebraic half of σ₁.
- (SQ) is satisfiable with room (e.g. a = 3/5, d₀ = 1, κ = 1/8, b₂ = 2: 3d₀/(32κb₂) = 3/8 < 3/5). Increasing b₂ relaxes
  (SQ) and lowers σ₁, but raises σ₂ (larger children); the cost optimum (`params.optimise_2L`) balances these.

## 6. Recursion, search and Lemma IR [proved]

- **Value factorisations.** Batch smoothness (Bernstein) to B decides the rational side and the B-part of F(P); the
  cofactor c₀ := |F(P)|/(B-smooth part), 1 ≤ c₀ ≤ X_F, is factored by FACTOR_U (Theorem U) with fresh coins, and the
  canonical factorisation is checked against Acc. No sieving or trial division to B₂. Expected cost per call
  ≤ L_{X_F}[2/5, c_U + o(1)] = exp(O((θ_F log N)^{2/5}λ^{3/5})) = exp(O(D^{7/20}λ)) = L[3/8]^{o(1)}
  (7/20 < 3/8). (Codex's alternative — self-recursion with a three-term recurrence F(x) ≤ L^p + L^{δ₀}F(x^{1/2+o(1)}) +
  L^{δ₁}F(x^{τ_x}), τ_x → 0 — also closes; using Theorem U as a black box is simpler.)
- **Lemma OS** (canonical outputs ⇒ the outer transcript's law is independent of the oracle) and **Lemma W** (Wald for
  adaptive fresh-coin sub-calls) apply verbatim; hence (T) and Theorem DD's fruitfulness (multidescent Theorem 2L (v)).
- **Search (escapeR2 §7.3 recurrence at depth 2).** For each of the ≤ W₁ level-1 candidates: batch work plus ≤ 1
  FACTOR_U call; for each accepted point and each of its ≤ m̄ children: a child scan (B^{1+o(1)}, cached per window)
  and ≤ W₂ candidates with batch work. Deterministic count ≤ W₁(1 + m̄W₂) candidate examinations; expected cost
  T_SQ ≤ W₁(1 + m̄W₂)·L^{o(1)} = L[3/8, σ₁ + σ₂ + o(1)], uniformly in the field and the target, including failed searches.
- **Lemma IR at γ₂.** With (H), log H/(13d log log H) ∼ D^{1/2}/(13d₀²·(3/4)) → the rejection probability
  exp(−Θ(D^{1/2})) + 3/(2H + 1) + 1/W is o(γ₂) for γ₂ ≥ (log N)^{−C}. As codex notes, this is used *globally* (as in
  select-U: the set 𝒢* = {A ∉ 𝓡 : δ_A ≥ γ₂/4} has measure ≥ γ₂/4), not conditionally on a rare incidence E. Lemma IR's
  primes r ≤ (1/2)log H are < B, and Q_IR = H^{1/2+o(1)} ≤ (2H + 1)/d.
- **Cost.** Lemma F_α at α = 3/8 with (10.1) of lv5unif, T_SQ as above: c₂L = max(2a, a + s + max(a, σ₁ + σ₂),
  (a + s)/(1 − 2^{−3/8})), s = 2κ(1 + C_{γ₂}); 1/(1 − 2^{−3/8}) ∈ [4.3900, 4.3901] (exact bracket via x⁸ vs 8).
  Illustration (rational grid, (SQ) imposed; C_{γ₂} not computable): c ≤ 4.23, 4.61, 5.35 for C_{γ₂} = 0, 1, 5.

## 7. The conditional child mean: Hyp CM [open; proof plan]

### 7.1 Statement and why it is not a level-1 statement

For an incidence (ℓ, P, q) let E = E_{ℓ,P,q} (multidescent §3.2), ν(ℓ, P, q) := Pr[ℓ, E] (total mass W ≤ m̄μ₁), and
  E[N_𝔮 | E] = Σ_{P′∈Box(𝔮)}Pr[E ∧ Good(P′)]/Pr[E].   (CM)

**Hyp CM.** There is C_CM such that Σ_{(ℓ,P,q)}ν(ℓ, P, q)·1[E[N_𝔮 | E] < m₂(q)(log N)^{−C_CM}] ≤ μ₁/log N.

With m₂ ≥ exp((log log N)²) this gives Σν·min(1, 1/E[N_𝔮|E]) ≤ 2μ₁/log N, i.e. the "1/E" part of PC₃'s ϕ is
negligible; the rest of PC₃ is the covariance (precision) part. **So CM is the first-order part of PC₃** and the only
part of multidescent's U′ not proved here. The numerator is a probability of a *two-point* event (P and P′), which
none of FM/LD*/GM bound from below: FM is one-point, LD* is a two-point *upper* bound.

### 7.2 Proof plan

(CM1) *Denominator, pointwise upper bound [outline; standard].* Pr[E] ≤ (log N)^C·Pr_g[rational filter at P]·q^{−β₂}ρ_{F2},
β₂ := α(X_F, B₂), pointwise in (ℓ, P, q): one-point smoothing for c, Lemma W-U for w (sum over all reduced residues mod t),
and #{F ∈ [−X_F, X_F] : q | F, F ∈ 𝒮_{B₂}} = 2Ψ(X_F/q, B₂) ≤ Cq^{−β₂}Ψ(X_F, B₂) (HT Theorem 3; q ≤ B₂). This is the
algebraic analogue of multidescent's Lemma CL.

(CM2) *Rational pair.* Fix (ℓ, P, n) with n = G(P)/ℓ in I_P ∩ 𝒮_B; on the fibre g = g* + jP (≥ M/(CV₁ log log N) values
of j, Lemma FIB), G(P′) = G_{g*}(P′) − jΔ, Δ := det(P, P′) = qk. Restrict Box(𝔮) to transverse points
|Δ| ≥ |P||P′|/4 (a positive proportion). Then G(P′) runs over an AP of difference |Δ| whose range has length ≍ MV₂, the
value scale: so "Δ good" (max over residues, B-smooth values at scale X_{r2}, grid as in def:good) gives the count of
smooth G(P′) along every fibre, uniformly in n and the residue, with polylog loss.
  *Δ-goodness.* Moduli qk, k ≤ K := 2V₁V₂/q. Since q > B, E2 (lift averaging at the prime q) splits the error into
  E(x, B; k, ·)/(q − 1) — controlled by an all-moduli smooth BV over k ≤ K (Drappeau Lemme 5 is itself a sum over all
  moduli ≤ Q of all characters of conductor ≤ x^{η}; large conductors by Corollary LS with R₁ = L[3/8, C] ≤ x^{η};
  level K ≤ x^{o(1)}) — plus characters with q | cond. For q ≥ L[3/8, C] the latter are handled by the large sieve with
  R₀ = q. For q ∈ (B, L[3/8, C]) (small children) the family {qk} is sparse (density 1/q) and its q-divisible small
  conductors are not controlled for a single q; **averaged over q in a dyadic window** the family {qk : q ∈ (Q₂, 2Q₂],
  k ≤ K} has density ≍ 1/log Q₂ and Drappeau's Lemme 5 (all moduli) applies directly. So Δ-goodness holds for all but an
  ε-fraction of (q, k) per window, ε = e^{−Θ(D^{3/8}/λ²)}.
  *Multiplicity.* Each (q, k) is hit by ≤ C′V₁V₂ pairs (P, P′) with |P′| ≤ V₂ (for each P′, the P with det(P, P′) = qk lie
  on a line with ≤ 2√2V₁/|P′| + 1 points in the box), against an average ≍ V₁V₂/(d log V₁) (b prime): ratio ≤ (log N)².
  With the pointwise bound (CM1), the ν-weight of incidences meeting a bad Δ is ≤ (log N)^C ε·(total) — the bad set is a
  set of *moduli*, deterministic given (P, P′), so ε is never divided by Pr[E] (§7.3).

(CM3) *Algebraic pair.* Given (g, w), c ↦ Φ(c) = (R_c(P), R_c(P′)) ∈ Γ_Φ = {z₂ ≡ κ^{d−1}z₁ (mod Δ)} (paper2 `lem:kern`(i)).
  *Two-point smoothing [to be proved].* Law(Φ(c)) is within TV O(d²V₂²/H) of a mixture of uniform measures on Γ_Φ ∩ K_m,
  K_m convex and containing translates of boxes [−y₁′, y₁′] × [−y₂′, y₂′], y_i′ ≍ HV_i^{d−1}/ϖ. (Sketch: Σ c_iw_i with
  w_i = (v_i(P), v_i(P′)); two coordinates with large entries give a parallelogram; Kern-U's triangular kernel basis,
  entries ≤ bb′ ≤ V₁V₂, fills the cosets of their sublattice in Γ_Φ with shifts of size ≤ dV₁V₂ ≪ H.)
  In a large layer, for every z₁ with F(P) = F̂(P) + G(P)z₁ acceptable and q | F(P), the admissible z₂ form
  {κ^{d−1}z₁ + jΔ} in an interval of length ≍ y₂′ ≫ |Δ|; then F(P′)/q runs over an AP of difference
  m₁₂ := |G(P′)|·|Δ|/q with ≫ 1 terms over the value range ≍ X_{F2}. So "m₁₂ good" (B-smooth values at scale X_{F2}, max
  over residues) gives the count for every z₁. *Level:* m₁₂ ≤ 4MV₁V₂²/q and X_{F2} ≥ c MHV₂^d/q, so m₁₂²·N^{η₀} ≤ X_{F2}
  ⟸ H ≥ 16MV₁²N^{η₀} — exactly what (H) supplies (h − μ = 4ν + η₀ > 2ν + 2·log S₁/log N + η₀).
  *Residues.* Mod |G(P′)| the class is reduced (Lemma `lem:coprime` for P′, restricted to gcd(b′, g₁G(P′)) = 1, a
  constant-factor loss); mod Δ/q it is κ^d·F(P)/q, reduced unless a prime of Δ/q divides F(P); a prime ≤ B is a local
  factor; a prime in (B, B₂] is another child q₂ of P with q₂ | det(P, P′), in which case P′ is dead (q₂ | F(P′)/q) — a
  fraction ≤ m̄/B of the box.
  *m₁₂-goodness.* Moduli m₁₂ = s·k with s = |G(P′)| B-smooth and k = |Δ|/q arbitrary: a Lemma-SM-type BV for the
  product family {s·k : P⁺(s) ≤ B, k ≤ K}, with E3/E8 modified so that the Möbius variable and the conductor may divide
  k (weights Σ_{s,k : rd | sk}1/φ(sk) ≲ ρ₁(Q_s)·Σ_{r=r₁r₂}r₁^{1−β}τ(r)/(φ(r₁)φ(r₂))); the transfer factor is again
  R₁^{1−β₁} (as in (H4)); large conductors by Corollary LS. [Not written; the E8 variant is the one genuinely new
  computation.]
  *Multiplicity, the costly step.* The conditional multiplicity of a modulus value s = |G(P′)| is
  #{g : G(P′) = s, G(P)/ℓ ∈ 𝒮_B}: smooth numbers in an AP of difference |Δ| ≈ L[5/8]·L[1/2] — far beyond DS-U's
  q ≤ B^{C₀}. Divisor switching still gives the main term times Δ^{1−α(X_r, B)} = exp(O((log Δ/log B)λ)) = exp(O(D^{1/4}λ)).
  This loss is super-polylogarithmic, but it multiplies only the bad set, of relative size ε = e^{−Θ(D^{3/8}/λ²)}:
  D^{3/8}/λ² ≫ D^{1/4}λ.

(CM4) *Assembly and Markov.* Good mass ≥ (main)/(log N)^C, bad mass ≤ (main)·ε·e^{O(D^{1/4}λ)}; so the ν-weighted fraction
of incidences whose box has bad fraction > ε^{1/2} is ≤ ε^{1/2}e^{O(D^{1/4}λ)} ≪ 1/(m̄ log N); for the others, (CM2) and
(CM3) at good moduli give Σ_{P′}Pr[E ∧ Good(P′)] ≥ Pr_g[rational]·q^{−β₂}ρ_{F2}·m₂/(log N)^C, and (CM1) gives Pr[E] ≤
the same main term times (log N)^C. That is CM. Small children are handled through window averages, then a Markov step
over q within the window using the upper bound of the label law (Lemma CL type).

### 7.3 Codex's objection, and where the plan could fail

Codex (uprime_1 §5): "discarding bad moduli before conditioning is insufficient: an unconditional error ε can become
ε/Pr(E), exponentially large." This is right for the naive argument (bound the unconditional probability of a bad event,
then condition on E, of probability e^{−Θ(Λ)}). The plan avoids it: badness is a property of *moduli* (Δ = det(P, P′),
deterministic given the points; m₁₂ = |G(P′)|·k), the bad set is counted in the incidence measure itself with
multiplicities bounded pointwise (CM1 and the e^{O(D^{1/4}λ)} divisor-switching bound), and only ratios of masses in that
same measure are formed. What can fail:
- the two-point smoothing lemma (layers must be large in *both* directions, y₂′ ≫ |Δ|; true by (H), but the 2D
  log-concave/coset-filling argument is not written);
- the product-family E8 variant (the conductor and Möbius variables may divide the unrestricted k; the sum over
  r₁ | k of φ(r₁)-weighted primitive characters must still carry 1/φ(r₁) — it does when k is averaged, not for fixed k);
- the uniformity of Drappeau's Lemme 5 for the all-moduli sum at level K = x^{o(1)} (fine as stated: Q ≤ x);
- bookkeeping of the children's label law at small q (needs Lemma CL's algebraic analogue, outlined in multidescent §2.2,
  UNVERIFIED there).
None of these is an identified obstruction. Our estimate that CM is provable along these lines by a focused effort:
50–70%. Codex (blind): unproved, "not a demonstrated counterexample"; it did not assess the plan.

**Coordination with pc3.** CM is the 2-point (first-order) part of PC₃'s weighted form. If pc3 proves PC₃ via
"3-form count / (1- and 2-form counts) = 1 + O(η)", the 2-form counts it normalises by are CM's numerators, and their
lower bound is CM. No pc3 notes existed at the time of writing; we did not edit their files.

## 8. Achievable α

U′ items 1–5 hold at α = 3/8 (§§2–6), so no retreat into (3/8, 2/5) is needed for them. Codex's general-α check (one
level of algebraic children, d ≍ D^{1−2α}, log B₂ = b₂D^{2−4α}λ) shows the hybrid transfer loss is
O(D^{2−5α}λ) = o(D^α/log²D) for every α > 1/3, σ₁(α) = α/(ad₀) + αd₀/(4κb₂), and (SQ) becomes ab₂ > αd₀/(4κ). CM is open
on the whole interval, so moving α up does not help with what remains. (Deeper descents j ≥ 3 would need hybrid BVs
with three or more smoothness parameters at the deeper levels; not examined.)

## 9. Experiments (pre-registered, `notes/uprime_prereg.md`, commit `f5a9bd7`)

Toy checks only.

### 9.1 E1: hybrid BV at toy scale (x = 2^28, y₂ = 200, Ψ(x, y₂) = 4 885 929; 20 s, 0.25 GB)

| Q | y₁ | #F₁ | #F_all | bad F₁ | H(F₁)/H(F_all) | med V ratio | med z(F₁) |
|---|---|---|---|---|---|---|---|
| 2^12 | 7 | 209 | 239 | 0 | 0.998 | 1.062 | 0.31 |
| 2^12 | 13 | 239 | 239 | 0 | 1.000 | 1.000 | 0.31 |
| 2^13 | 7 | 301 | 479 | 0 | 0.976 | 1.078 | 0.40 |
| 2^13 | 13 | 387 | 479 | 0 | 0.990 | 1.060 | 0.39 |
| 2^14 | 7 | 456 | 958 | 0 | 0.976 | 1.076 | 0.46 |
| 2^14 | 13 | 617 | 958 | 0 | 0.978 | 1.054 | 0.46 |

**P1a–P1d all held.** Disclosure: in the cell (2^12, 13) every cofactor s ≤ 16 is 13-smooth, so F₁ = F_all and the
cell is trivially 1 (not noticed at registration). The variance ratio V ≈ 0.12–0.25 is strongly sub-Poisson for all
families (smooth numbers are more evenly spread over reduced classes than random sets); y₁-smooth cofactors are
5–8% higher than all moduli, consistently, but inside the band. Reading: smooth modulus cofactors are not special at
this scale, as the hybrid lemma requires asymptotically.

### 9.2 E2: mixed acceptance and the conditional child mean (3 600 boxes, 1.15·10⁶ points, 7 493 children; 63 s)

| statistic | value | band | verdict |
|---|---|---|---|
| P2a product rule: #both/Σ(#rat·#alg2/#pts) | 1.041 [1.027, 1.055] (field clusters) | [0.8, 1.3] | held |
| P2b removed (square + disc)/#both | 0.0039 (19 + 24 of 10 933; 0 for q | N) | ≤ 0.02 | held |
| P2c R_CM conditional/fresh-field null | 1.208 [1.167, 1.253] (parents); [0.97, 1.46] (fields) | ∈ [0.8, 1.5], CI ⊂ [0.6, 2.0] | held |
| P2d irregular child lattices | 0 of 7 493 | ≤ 0.0052 | held |
| P2e Pr[N_𝔮 = 0]/e^{−mean} | 1.056 (0.394 vs 0.374) | [0.8, 1.3] | held |

*Post hoc (not pre-registered; `src/uprime/posthoc.py`).* The registered null drew labels from fresh fields. Fields
contribute children in proportion to how many level-1 points they accept, so better fields are over-represented
(size bias; cf. multidescent E3's height/Jensen effect). Against a **same-field** null (a random matched label of the
parent's own field), R_CM_within = 1.005 [0.975, 1.035] (parents), [0.973, 1.044] (fields), per-field 0.90–1.20. So, at
toy scale, conditioning on the parent's acceptance neither depletes nor enriches its children's boxes beyond field
quality — the behaviour CM asserts, with polylog room to spare.

## 10. Literature

Inputs used (all as in Theorem U unless marked): Drappeau, Compos. Math. 151 (2015) 828–862, Lemme 5 (small conductors;
proof refers to Harper arXiv 1208.5992 §3.3) and Lemme 1(i) (from de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005),
Thm 2.4(i)); Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Thms 1, 2(i), 3; Tenenbaum's Ψ_q theorem; Iwaniec–Kowalski
Thm 7.13; Canfield–Erdős–Pomerance 1983; Schwartz–Zippel (new here, elementary); an explicit prime number theorem
(new here, for Lemma CS′'s computation; e.g. Rosser–Schoenfeld 1962/1975 or Trudgian, Ramanujan J. 2016 — exact form
UNVERIFIED, only the shape C₁t·exp(−c₂√log t) with explicit constants is used). Hybrid statements with different
smoothness for moduli and values: Granville (1993) and Fouvry–Tenenbaum (1996) prove BV for smooth numbers over all
moduli; Harper (1208.5992) and Drappeau 2015 are single-parameter; we found no prior "smooth moduli / smoother values"
restricted BV, but this is a search result to be checked by the literature critic. Lenstra–Pomerance 1992 remains the
best rigorous bound for arbitrary N we know (L[1/2, 1 + o(1)]).

## 11. Reviews

- **GPT-6.1-sol blind round `codex/uprime_1`** (xhigh, launched before this note): verdict "PROVED EXCEPT the
  conditional child first moment". Independently: the hybrid BV with the same new transfer condition (its (T) = our
  (H4)); the same σ₁; Rankin at B₂ without enumerating primes to B₂; the constraint ab₂ > 3d₀/(32κ); the conservative
  child geometry (reduced-basis boxes can have physical radius ≍ q) and the **height floor (H)** for small children's
  level-1/2 conditions (adopted here); IR used globally; the two exact obstacles of CM (AP of modulus Δ with unrestricted
  Δ/q; child AP of modulus |G(P′)||Δ|/q) and the ε/Pr(E) warning (answered in §7.3); general-α extension.
- Codex adversarial round, Opus proof critic, Opus literature critic: pending (this is v1).

## 12. Candidate card

- **Mechanism.** Theorem U's relation-supply proof with two smoothness bounds: B for moduli/rational sides, B₂ = L[1/2]
  for the level-1 algebraic value whose large primes are children.
- **Result.** Lemma 8.5H (hybrid restricted BV, proved), FM-U′/LD*-U′/D-U′, Lemma CS′, Lemma SQF, child scans,
  recursion and IR at 3/8 (proved, internal chain); constraint (SQ) and height floor (H). Hyp CM open (proof plan).
- **Failure scope.** None found for U′'s level-1 part. CM is the remaining first-order input; PC₃'s covariance part is
  pc3's.
- **Next step.** Write CM's two-point smoothing lemma and the product-family E8 variant; then CM, then PC₃'s covariance.

## 13. Reproduction

```
.venv/bin/python -m pytest -q tests/test_uprime.py
.venv/bin/python -m src.uprime.hybrid_bv experiments/raw/uprime_e1.jsonl
.venv/bin/python -m src.uprime.analyze E1 experiments/raw/uprime_e1.jsonl
.venv/bin/python -m src.uprime.experiments experiments/raw/uprime_e2.jsonl 6 1
.venv/bin/python -m src.uprime.analyze E2 experiments/raw/uprime_e2.jsonl
.venv/bin/python -m src.uprime.posthoc experiments/raw/uprime_e2.jsonl experiments/raw/uprime_e2_posthoc.json
.venv/bin/python -c "from src.uprime import params as P; from fractions import Fraction as F; print(P.margins(F(3,5),1,F(1,8),2))"
```
