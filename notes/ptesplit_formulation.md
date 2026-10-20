# ptesplit: exact formulation of route 1's open escape (systeq §6)

Agent `ptesplit`, phase 16, 2026-10-20. This is the shared formulation that `ptebarrier` and `ptecons` read.
Labels: PROVED (proof here), HEURISTIC (assumption named), KNOWN (cited).

Notation: x = log N, b = log B, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}), B = L[α].

## F1. Objects

**Relation.** A pair (u, v) of nonzero integers with u ≡ v (mod N) and u ≠ v. Write u − v = kN, k ≠ 0.
- Codex aggregate bound [PROVED, systeq §4]: max(|u|, |v|) ≥ N/2.
- u is *split by construction*: u = ∏_{j=1}^M π_j(m) for pieces π_j from a factor base 𝒫 of size |𝒫| ≤ L[α]
  (small integers, or values b m − a of linear forms, or elements of an order mapping to Z/N).
- v must be B-smooth. Smoothness of a random integer of log-size ℓ costs ρ(ℓ/b) = exp(−(ℓ/b) log(ℓ/b)(1+o(1))).

**Presentation.** P = (O, φ): an order O = Z[X]/(F̃) (F̃ ∈ Z[X] monic, degree d, height H; or any order of a number
field, or Z itself when d = 1) and a ring map φ : O → A to a Z/N-algebra A (usually A = Z/N, φ(X) = m with N | F̃(m)).
The *free* presentation is A = (Z/N)[X]/(F) with F̃ a monic lift (systeq Lemma EX).

**Split member.** P(X) = ∏_{j=1}^M (b_j X − a_j), with (a_j, b_j) ∈ [−B, B]² (the monic case b_j = 1 is systeq's).
Its *reduction* is r_P = P mod F̃ ∈ Z[X], deg r_P < d, height h. Its relation is

  ∏_j (b_j m − a_j) ≡ r_P(m) (mod N),  v = r_P(m), |v| ≲ d·h·m^{d−1}.

**PTE fibre.** For K = M − k: 𝓕_K(h) = {multisets a ∈ [−B, B]^M : e_i(a) = h_i, i ≤ K}
(equivalently p_i(a) = h_i, i ≤ K, by Newton in characteristic 0). Members of one fibre differ by
P − P′ of degree < k. Ideal PTE of size M is the case K = M − 1 (k = 1).

**Family.** 𝓕 is a set of split members, all computed from N and public data in time T.

## F2. Two distinct mechanisms by which N enters (these must not be conflated)

- **(A) Exact presentation (one congruence compresses the family).** Each member satisfies P ≡ r_P mod F̃ *exactly in
  Z[X]*; N enters only through φ (N | F̃(m)). This is the literal BGJT shape: X^q ≡ h₀/h₁ in the presentation.
- **(B) Seeded additive compression (the PTE shape).** One seed relation P₁(m) ≡ v₁ (mod N), found by any means; every
  P in the fibre of P₁ gives P(m) ≡ v₁ + (P − P₁)(m), an exact integer identity on top of the seed. All members share
  the same k = (u − v)/N. Relations per seed ≤ |fibre|.
- **(C) Non-exact free presentation.** P − r_P ∈ (N, F̃) with (P − r_P) mod F̃ ≠ 0, so ‖·‖ ≥ N (EX). This is
  systeq's G_aux case: d conditions mod N.

## F3. Targets (cost model of smoothenum/heurnfs: LA quadratic in |𝒫|, relations needed ≈ |𝒫|)

Algorithm cost L[α'] needs at least |𝒫| ≈ L[α] useful relations at total cost L[α'], α' = α.

**(T-v) Size of the small side.** v must be L[α]-smooth with probability L[α]^{−O(1)}:
  log|v| ≤ c·x^{2α}(log x)^{1−2α}.
For α = 1/3, 1/4 and 1/10 this means log|v| ≲ x^{2/3}, x^{1/2} and x^{1/5}. "Independent of x" is stronger than needed.

**(T-u) Size of the split side.** |u| ≥ N/2, so Σ_j log|π_j(m)| ≥ x − 1. With pieces of log-size λ: M ≥ x/λ.

**(T-N) Number.** Σ over the members of the per-member success probability (v smooth) must reach |𝒫| ≈ L[α].

**Window lemma (HEURISTIC sizes; PROVED inequality).** In mechanism (B) with pieces m − a_j (λ = log m > b):
(P − P₁)(m) = Σ_{l<k} Δe_{M−l} m^l with |Δe_{M−l}| ≤ 2·C(M, l)·B^{M−l}. For *generic* Δe, (T-v) and (T-u) force
  (k−1)λ ≤ x^{2α} and (M−k+1)b ≤ x^{2α}, hence λ ≥ x^{1−α} and k − 1 ≤ x^{3α−1}.
So for α < 1/3 the fibre must be an **ideal-PTE fibre (k = 1)**: M ≈ x^α, B = L[α], pieces of log-size x^{1−α}.
Escaping k = 1 needs Δe much smaller than generic (structured families). At α = 1/3, k = O(1).

**What would beat L[1/3].** The cheapest sufficient package, in mechanism (B), is all of:
- (i) ideal-PTE fibres with M ≈ x^α → ∞, entries ≤ B = L[α], and |fibre| ≥ L[α, c] for some c > 0 (the trivial bound is
  2B + 1);
- (ii) seeds: about |𝒫|/|fibre| seed relations of the special shape P₁(m) ≡ small (mod N), found at total cost L[α];
- (iii) the v-sides v₁ + c are L[α]-smooth at rate L[α]^{−O(1)}.

In mechanism (A), what is needed is many small r with r(θ) split into linear elements of the box, and r(m) smooth
(see F4).

For L[1/4] and L[1/10], put α = 1/4 or 1/10 in the above: M ≈ x^{1/4} or x^{1/10}, k = 1.

## F4. Two facts that the other agents can use

**Lemma NC (norm compression) [PROVED].** Let F̃ ∈ Z[X] be monic and G a monic irreducible factor of degree δ with
N | G(m). If no such G exists, then some gcd(N, G_i(m)) is a proper factor of N. Let θ be a root of G and
P = ∏_j(b_jX − a_j) ≡ r (mod F̃) in Z[X]. Then P(θ) = r(θ) in Z[θ], so

  ∏_j |Ĝ(a_j, b_j)| = |N_{Q(θ)/Q}(r(θ))| ≤ ‖r‖₁^δ · ‖F̃‖₂^{d−1},

where Ĝ(a, b) = b^δ G(a/b) ∈ Z. For δ ≥ 2, Ĝ(a, b) ≠ 0.
- Proof of the bound: |r(θ_i)| ≤ ‖r‖₁ max(1, |θ_i|)^{deg r}, and ∏_i max(1, |θ_i|) = M(G) ≤ M(F̃) ≤ ‖F̃‖₂ (Landau).
- Hence at most log₂(‖r‖₁^δ‖F̃‖₂^{d−1}) pieces are non-units. Every other piece is a unit bθ − a of Z[θ]:
  - at most 2δ monic ones;
  - O(δ) primitive ones for δ ≥ 3 (Thue; Bombieri–Schmidt 1987 count);
  - O(log B) in the box for δ = 2.
- Members with the same non-unit part differ by exact unit relations, which are uninformative (EX).
- So every mechanism-(A) member is an NFS relation for the pair (G, X − m), sieving element r(θ) of degree ≤ d − 1,
  with an element (not ideal) factorisation. It is subject to heurnfs B1/B2 (α ≥ 1/3 under G).
- Function-field contrast: in BGJT the identity holds in the finite ring F_{q²}[X]/(I), where every piece X − a is a
  unit, so the norm carries no size.

**Prop NF (N-free families are useless for almost all N) [PROVED].** Let S be any set of integers of bit-size ≤ T,
fixed before N. The number of balanced semiprimes N = pq ∈ [X, 2X] (p, q > √X/2) with some s ∈ S, |v| ≤ V, v ≠ s,
N | s − v is at most |S|(2V+1)(2T/log₂(√X/2))². That is a fraction e^{−x+o(x)} when |S|, V, T = e^{o(x)}.
- Proof: s − v ≠ 0 has at most T/log₂(√X/2) prime factors above √X/2, so at most that number squared such divisors pq.
- Consequence: N must enter the construction of the family (seed, or presentation).

## F5. Heuristic count (generic fibres)

Multisets: C(2B+M, M) ≈ (2B)^M/M!. Prescribing e₁..e_K costs ∏_{i≤K} 2C(M,i)B^i ≈ B^{K(K+1)/2}.
- The expected size of a generic fibre is ≈ (2B)^M / (M!·∏_{i≤K} 2C(M,i)B^i). It is < 1 as soon as K(K+1)/2 > M + o(M).
- For k = 1 that means M ≥ 3 (up to constants).
- Pair count: the expected number of nondiagonal pairs is ≈ (2B)^{2M}/(M!²·B^{K(K+1)/2}), which is → 0 for K ≥ 4 when k = 1.
- **Dense fibres are therefore non-generic for every M ≥ 5 with k ≤ M/2**: any family must come from structure.
- Trivial L^∞ bound [PROVED]: |𝓕_K(h)| ≤ C(2B + k, k). Choose any k-submultiset; the complement is determined by
  Newton's identities.
