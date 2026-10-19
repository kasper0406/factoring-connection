# lv6: the LV §6 adaptation (LV6) for design (a″_k), and the HECM caveat

Agent `lv6`, phase 9, 2026-10-19.
- Pre-registration: `notes/lv6_prereg.md` (commit `1297e8e`, before E1 ran).
- Code: `src/lv6/` (`numfield.py`, `sqrt.py`, `chars.py`, `pipeline.py`, `e1.py`). Tests: `tests/test_lv6.py`.
- Raw data: `experiments/raw/lv6_e1.jsonl`, summary `experiments/raw/lv6_e1_summary.json`.
- GPT-6.1-sol: blind round `codex/lv6_1.{md,out}`; adversarial round `codex/lv6_2.{md,out}` (§9).
- Opus critics: one on proofs, one on literature (§9).

## Verdict

**LV6 is proved (Proposition LV6, §2–§3). The HECM caveat is removed (§6). If Hypothesis LD were proved, Theorem A
would hold with c = (1+θ)/√2 (infimum (1 + 3/(d+1))/√2, i.e. 0.9192 at d = 9), and every other ingredient would be
a proved result or a published theorem (§7).** LD is still open (agent `ldbound`). Nothing here is refereed.

In plain terms:
- LV's §6 was written for degree d ≍ (log N)^{1/3} and heights L(1/3). Our degree is **fixed** (d = 9), and that
  makes everything easier. Every field-dependent quantity (discriminants, the index, the class number, the number
  of bad primes) is N^{O(1)}. The obstruction group H therefore has dimension O(log N), and O(log N) character
  columns suffice.
- The height V = N^{1/(4k)} enters in one place only: the bit size of the element whose square root is taken,
  O(rows·log N) bits instead of LV's O(rows·(log N)^{1/3}). Square roots cost quasi-linear time in that size, so
  this is harmless.
- The square root needs no inert prime and no lattice reduction. We Hensel-lift at **one** unramified prime
  p > B and try the 2^{t−1} ≤ 2^{d−1} sign patterns. Each candidate is checked exactly (η² = Γ). Fixed d makes
  this rigorous and cheap.
- Siegel zeros are handled as in LV (Stark plus Landau), but at fixed d they cost even less: **at most one**
  nonzero square class can escape the characters, for every field. Testing the three dependencies s₁, s₂,
  s₁ + s₂ then always finds a square. (Codex found this sharpening; LV had O(log log N) exceptional classes and an
  exceptional set of fields.)
- The non-monic pair (f, g₁x − g₀) needs **one** parity column. It makes both g₁^E and f_d^E squares. Bad primes
  (dividing f_d·disc f) get **no** columns; they are absorbed into H.
- The whole post-processing is a function of the certificates (R_i, splits, special-q relations, field) and fresh
  coins. It never sees a hidden root r_i, so Theorem DD's transcript invariance (T) holds.
- Cost: linear algebra L^{2a+o(1)} = (rows)^{2+o(1)}. Characters, square root and the rational square root cost
  L^{a+o(1)}.

**HECM.** Theorem A does not need it:
- (i) **Batch smoothness [proved here; elementary].** Every special-q search enumerates a complete box of ≤ C_abs·S²
  candidates (pcest Lemma 1.3). Testing all of them at once with Bernstein's smooth-part algorithm costs
  Õ(B + S²·log N) per search, deterministically, and recognizes *every* smooth value. Rows × B is B² = the
  linear-algebra cost, so the exponent never changes. At Theorem A's optimum B ≤ S² anyway.
- (ii) **Unconditional HECM [proved modulo genus2rh's internal Props A/B].** genus2rh's separation argument plus
  LPP I Thm 1.3 (refereed, 1993) and LPP II Thm 1.1 (refereed, 2002) gives a smoothness test for all B-smooth
  inputs at cost L_B[2/3, O(1)]·(log N)^{O(1)} = L^{o(1)}. We write out the bookkeeping that genus2rh §6.3 left as
  a sketch (§6.3). Route (i) does not depend on it.
- (iii) **Lenstra–Pomerance 1992 does not use HECM.** It uses ECM with Pomerance's 1987 "recognizable smooth
  numbers" (LP §7, read from the paper). So the shape test was always HECM-free.
- **Headline:** Theorem A (0.919 at d = 9) is the honest conditional headline. Theorem A′ (Pollard–Strassen, 0.967)
  is superseded.

**LV6 needs no smoothness tests at all.** The rows arrive factored. The only new integer tests are primality tests
of the character primes and of the Hensel prime.

**Changes to pcest's theorems: none to FM, R, Kern or 4.3.** Theorem A's *assembly* gains four field-rejection
events (§4). Their total probability is ≤ 1/(2H+1) + exp(−Ω(log H/log log H)), so γ drops by at most a factor 2.
Lemma IR (§4) is an elementary proof of the irreducibility bound we need. It replaces the "LV Lemma 5.1 transfer,
checked at the level of statements" caveat.

**Experiment E1 (pre-registered, exact, toy; §8).** 16 instances, d ∈ {3, 5}, 36–48-bit N, up to 6000 relations.
- 116 certified squares: 0 congruence violations (P1).
- 92 non-squares, all certified (P2).
- With 40 characters every tested full-kernel vector was a square in 16/16 instances (P3), and the combination
  search found a square in 16/16 (P4).
- Without characters, the parity-only kernel is mostly non-square: character rank 1–16, square fraction ≤ 1/2 in
  16/16 (P5). So characters are needed, and the bad-prime classes are real.
- The toy checks correctness only. It says nothing about asymptotics.

Labels used below: **proved** (proof in this note), **cited** (exact reference, hypotheses checked at our
parameters), **proved (internal)** (proved in another project note, unrefereed), **open**.

---

## 1. Setting: certificates of a pure-descent row

Notation follows pcest §1 and hypD §1–§3.
- d ≥ 3 is odd and **fixed**; k is fixed; s = 1/(d+1); M = ⌊N^s⌋; H = N^h; V = ⌊√ℓ·S⌋ = N^{1/(4k)}L^{O(1)}.
- B = L[1/2, a]. After preprocessing (trial division to B, perfect-power test), N has no prime factor ≤ B.
- The field is A = (g₀, g₁, c), f = F̂_g + (g₁x − g₀)R_c, with f = Σ f_i x^i, F(a, b) = Σ f_i a^i b^{d−i}, and α a
  root of f.
- Write t := f_d and β := tα. The monic polynomial g(y) := t^{d−1}f(y/t) = Σ f_i t^{d−1−i} y^i has root β.
  Let Z[β] = Z[y]/(g) ⊂ O_K and D := t·disc(f).
  - disc(g) = t^{(d−1)(d−2)} disc(f), so the primes dividing disc(g) are exactly those dividing D (d ≥ 3).
- ψ: Z[β] → Z/N, β ↦ t·g₀g₁^{−1}. It is a ring homomorphism, because g(t g₀/g₁) = t^{d−1} g₁^{−d} F(g₀, g₁)
  = t^{d−1} g₁^{−d} N ≡ 0.
- For an atom P = (a, b) (a special-q relation), put z_P := t a − bβ ∈ Z[β]. Then
  Norm_{K/Q}(z_P) = t^{d−1}F(a, b) and ψ(z_P) = t·G(P)/g₁, where G(P) = g₁a − g₀b = σ_P ℓ n_P, with σ_P = ±1 and
  n_P > 0 B-smooth.

**Rows.** Row i comes from a uniform r_i, R_i = r_i², and an accepted split R_i ≡ u/v with u, v products of
distinct primes ℓ. Each ℓ has one certified atom P_ℓ, and e_{iP} = +1 if ℓ | u, −1 if ℓ | v. Since
ℓ ≡ σ_P g₁ ψ(a − bα)/n_P (mod N),

  r_i² ≡ q_i · ∏_P ψ(a − bα)^{e_{iP}} = q_i · t^{−E_i} ∏_P ψ(z_P)^{e_{iP}}   (mod N),
  q_i := ∏_P (σ_P g₁/n_P)^{e_{iP}} ∈ Q^×,  E_i := Σ_P e_{iP}.   (1.1)

A row is a DD descent relation with X_i = r_i (W_i = 1). Its certificate is c_i = (q_i, (z_P, e_{iP})_P). The same
atom may occur in several rows (cached special-q relations). For a set S of rows, write e_P(S) := Σ_{i∈S} e_{iP}
∈ Z and E(S) := Σ_P e_P(S).

## 2. Proposition LV6

**Proposition LV6 [proved, given the cited theorems listed in §3].**
- **Hypotheses.**
  - Fix odd d ≥ 3 and k.
  - N > N₀(d, k) has been preprocessed.
  - The field passes the checks of §4: gcd(g₁ t, N) = 1, f irreducible and primitive, and ψ(g′(β)) is a unit
    mod N; otherwise a factor of N was output or the field was rejected.
  - We are given M_rows pure-descent rows as in §1, with M_rows ≥ C_col + 12. Here
    C_col := π(B) + #{(r, s)} + 2 + t_χ ≤ (d+1)π(B) + 2 + t_χ, and t_χ = O_d(log N) is the character count of §3.3.
- **Algorithm.** There is a Las Vegas algorithm. In expected time (M_rows)^{2}·(log N)^{O_d(1)}, it outputs a
  nonempty set S of rows and Y ∈ (Z/N)^× with X_S² ≡ Y² (mod N), where X_S = ∏_{i∈S} r_i.
- **Invariance.** S and Y are functions of the certificates, the field and the algorithm's own coins only. So the
  exposed transcript satisfies Theorem DD's hypothesis (T).
- **Consequence.** Pr[1 < gcd(X_S − Y, N) < N] = 1 − 2^{1−ω(N)} ≥ 1/2 (lemmaM Theorem DD).

**Algorithm LV6.**
1. *Columns* (§3.1). The columns are: the sign of q_i; for each prime p ≤ B, the parity of v_p(q_i); one parity
   column E_i mod 2; for each prime r ≤ B with r ∤ D and each root s of f mod r, the parity of
   Σ_P e_{iP}·[a_P ≡ b_P s mod r]·v_r(F(a_P, b_P)); and t_χ character columns (§3.3).
2. *Characters* (§3.3). Compute x from the explicit discriminant bound. Sample t_χ first-degree primes (r_j, s_j)
   of K independently and uniformly among those with norm in (x/2, x]. Column j holds
   Σ_P e_{iP}·χ_j(a_P − b_P α) mod 2, where χ_j(a − bα) = ((a − b s_j)/r_j).
3. *Linear algebra* (§3.4). Sample two left-kernel vectors s₁, s₂ of the M_rows × C_col matrix over F₂, verify
   them, and retry until they are linearly independent.
4. *Square test* (§3.5). For S ∈ (s₁, s₂, s₁ + s₂), in this order, run Algorithm SQ on
   Γ_S := g′(β)²·∏_P z_P^{|e_P(S)|}. Take the first S for which it returns η ∈ Z[β] with η² = Γ_S.
5. *Congruence* (§3.6). Set Y := √q_S · t^{−E(S)/2} · ψ(η) · (ψ(γ⁻_S)·ψ(g′(β)))^{−1} (mod N), where
   q_S = ∏_{i∈S} q_i, √q_S is computed by halving exponents, and γ⁻_S := ∏_{e_P(S)<0} z_P^{|e_P(S)|}.
6. If step 3 or step 4 fails (probability ≤ 2/(log N)² + 2^{−9} per attempt; §3.3, §3.4), resample the characters
   and kernel vectors.
7. **After a trivial gcd, restart with fresh hidden roots** (new rows). Do not try another dependency of the same
   rows that was chosen in view of the gcd outcome (§3.7).

## 3. Proof

### 3.1 Columns: Lemma C (good primes)

**Lemma C [proved].** Let r be a prime with r ∤ D, and let (a, b) be coprime with F(a, b) ≠ 0. Then:
- (i) if r | b, then r ∤ F(a, b), and a − bα is a unit at every prime of O_K above r;
- (ii) if r ∤ b, the only prime 𝔭 of O_K above r with ord_𝔭(a − bα) ≠ 0 is the first-degree prime
  𝔭 = (r, α − ab^{−1}). It occurs only if r | F(a, b), and then ord_𝔭(a − bα) = v_r(F(a, b)).

*Proof.*
- Since r ∤ t, α is r-integral. Z_(r)[α] has discriminant disc(f)/t^{2d−2} up to an r-unit, which is prime to r.
  So Z_(r)[α] is the maximal order at r, and by Dedekind–Kummer the primes above r are (r, φ_j(α)), one for each
  irreducible factor φ_j of f mod r, all unramified.
- (i) If r | b, then r ∤ a, and F(a, b) ≡ t a^d ≢ 0 (mod r). Also a − bα ≡ a (mod 𝔭) for every 𝔭 | r.
- (ii) a − bα ∈ (r, φ_j(α)) iff φ_j(ab^{−1}) ≡ 0, i.e. iff φ_j = x − ab^{−1}, which is a first-degree prime 𝔭.
  - Norm_{K/Q}(a − bα) = F(a, b)/t, and 𝔭 is the only prime above r dividing a − bα.
  - Since f(𝔭/r) = 1, ord_𝔭(a − bα) = v_r(F(a, b)/t) = v_r(F(a, b)). ∎
- Primes r > B with r ∤ D contribute nothing: r ∤ F(a, b) because F(a, b) is B-smooth.
- Hence, for any S whose algebraic (r, s)-columns sum to zero mod 2:
  ord_𝔓(γ_S) ≡ 0 (mod 2) for every prime 𝔓 of O_K not above a prime dividing D, where γ_S := ∏_P (a_P − b_Pα)^{e_P(S)}.
- Bad primes need no columns. Neither the index [O_K : Z[β]] nor the factorization of D is ever computed.
- (Tested: `test_good_prime_columns_lemma_C`.)

**Rational side.** If the sign, the p-parity (p ≤ B) and the E-parity columns of S sum to zero, then
q_S = g₁^{E(S)}·∏_P n_P^{−e_P(S)} (sign +1) is a square in Q^×, with explicit square root
√q_S = g₁^{E/2} ∏_{p≤B} p^{−(Σ_P e_P v_p(n_P))/2}.
- The **one** parity column also gives E(S) even, which makes t^{−E(S)} a square.
- So γ_S = t^{−E(S)}·∏_P z_P^{e_P(S)} is a square in K iff γ′_S := ∏_P z_P^{e_P(S)} is.
- This is the "g₁ column" of pcest §5. It counts atomic relations with multiplicity, and it handles f_d as well.

### 3.2 The residual square-class group: Lemma H

Let S_D be the set of primes of O_K above rational primes dividing D. Define
  V_D := {γ ∈ K^× : ord_𝔓(γ) even for all 𝔓 ∉ S_D},   H := V_D/K^{×2}.
By §3.1, every S in the kernel of the non-character columns has γ_S ∈ V_D. Write φ(S) := [γ_S] ∈ H. Then φ is
F₂-linear on the kernel, because γ_{S+S′} = γ_Sγ_{S′}/γ_{S∩S′}².

**Lemma H [proved; standard].** |H| ≤ 2^{|S_D|}·|Cl(K)[2]|·2^{r₁+r₂} ≤ 2^{d·ω(D) + d}·h_K. With the bound of
Lenstra (Bull. AMS 26 (1992), Thm 6.5) applied to M_K ≤ √|Δ_K| ≤ √|disc g|, this gives
  log₂|H| ≤ d·log₂|D| + d + ½log₂|disc g| + (d−1)·log₂(d + log₂|disc g|) = O_d(log N).
*Proof.*
- Put W := {γ : γO_K = 𝔞²} and Y := O_K^×K^{×2}. Then V_D ⊃ W ⊃ Y ⊃ K^{×2}.
- The map V_D → ⊕_{𝔓∈S_D} Z/2, γ ↦ (ord_𝔓 γ mod 2), has kernel W.
- The map W → Cl(K)[2], γ ↦ [𝔞], is well defined (the ideal group is torsion-free) and has kernel Y.
- Y/K^{×2} ≅ O_K^×/O_K^{×2} has order 2^{r₁+r₂} ≤ 2^d (Dirichlet; roots of unity ±1 for odd d).
- |S_D| ≤ d·ω(D) ≤ d·log₂|D|.
- Sizes. ‖f‖_∞ ≤ C_d·g₁H = N^{O(1)}, so |disc f| ≤ d^{2d}‖f‖^{2d−2} and log|D|, log|disc g| = O_d(log N). ∎
- **Compared with LV (Claim 6.5).**
  - LV use BLP's [V : W] ≤ [O_K : Z[f_dα]]. We use the cruder 2^{|S_D|}, which needs no l-functions on the
    order, no projective-root bookkeeping, and no index computation.
  - Units, class-group 2-torsion, the index and the non-monic leading coefficient are all inside H. None of them
    is assumed trivial.
  - The cost is a constant factor in the dimension: O(log N) either way.
- (`numfield.class_group_log2_bound` computes this bound; E1 reports 271–604 bits for the toy fields.)

### 3.3 Character columns: Lemma X

For 0 ≠ h ∈ H, let L_h := K(√h), so [L_h : Q] = 2d.

**Discriminant bound [proved].** |Δ_{L_h}| ≤ Δ̄ := |disc g|²·8^d·|D|^d.
- Δ_{L_h} = Δ_K²·N_{K/Q}(𝔡_{L_h/K}).
- Take a representative γ ∈ V_D. At 𝔓 ∤ 2D, ord_𝔓γ is even, so after removing a local square L_h/K is unramified
  at 𝔓.
- At 𝔓 | D, 𝔓 ∤ 2: ord_𝔓 𝔡 ≤ 1. At 𝔓 | 2: ord_𝔓 𝔡 ≤ 2·ord_𝔓(2) + 1.
- So N(𝔡) ≤ 8^d·∏_{𝔓|D} N𝔓 ≤ 8^d |D|^d, and |Δ_K| ≤ |disc g|. Then log Δ̄ = O_d(log N).
- No small representative of h is needed. LV's minimal-norm representative argument is not used.

**Exceptional classes.** Let c_L > 0 be the absolute constant in Landau's theorem below. Put
  ν₁ := min( 1/(4·(2d)!·log Δ̄), c_L/(2 log Δ̄) ),
and let E := {0 ≠ h ∈ H : ζ_{L_h} has a real zero β with 1 − ν₁ ≤ β < 1}.

**Lemma E [proved from cited results]. |E| ≤ 1, for every field.**
- **Cited inputs.**
  - Stark (Invent. Math. 23 (1974), Lemma 8): if a number field L has a real zero β of ζ_L with
    1 − (c(L) log|Δ_L|)^{−1} ≤ β < 1, where c(L) = 4[L:Q]! for non-normal L, then ζ_F(β) = 0 for some quadratic
    F ⊂ L.
  - Landau (1918; Davenport, *Multiplicative Number Theory*, Ch. 14; LV Fact 6.24): there is an effective absolute
    c_L > 0 such that, for distinct real primitive characters χ₁, χ₂ of conductors q₁, q₂, L(s, χ₁)L(s, χ₂) has at
    most one real zero in (1 − c_L/log(q₁q₂), 1).
- *Proof.*
  - Let h ∈ E with zero β. Since ν₁ ≤ 1/(4(2d)! log Δ_{L_h}), Stark gives a quadratic F_h ⊂ L_h with
    ζ_{F_h}(β) = ζ(β)L(β, χ_{F_h}) = 0. Since ζ(β) ≠ 0 on (0, 1), L(β, χ_{F_h}) = 0.
  - [L_h : Q] = 2d with d odd, so L_h has at most one quadratic subfield: two would generate a biquadratic subfield
    of degree 4 ∤ 2d. K has none, so L_h = K·F_h.
  - h ↦ L_h is injective (Kummer theory), hence so is h ↦ F_h.
  - |Δ_{F_h}| ≤ |Δ_{L_h}|^{1/d} ≤ Δ̄. Two distinct h, h′ ∈ E would give distinct real primitive characters with
    zeros ≥ 1 − ν₁ ≥ 1 − c_L/log(|Δ_{F_h}||Δ_{F_{h′}}|). Landau forbids this. ∎
  - (Codex lv6_1 found the bound |E| ≤ 1. LV's dyadic-block argument gives O(log log N), plus an exceptional set
    of fields. Neither is needed at fixed d.)

**Sampler.**
- Let log x := ⌈max(20d·(log Δ̄)², (log N)²)⌉.
- Draw r uniform in (x/2, x]. If r is prime (AKS), compute the set Z_r of simple roots of f mod r (Las Vegas root
  finding). With probability |Z_r|/d, output (r, s) with s uniform in Z_r; otherwise redraw.
- The output is uniform over first-degree primes 𝔭 = (r, α − s) of K with N𝔭 ∈ (x/2, x]. All of these have
  r > Δ̄ ≥ |D| and r > B.
- The expected number of draws is O(d log x). Each draw costs (log x)^{O(1)} = (log N)^{O(1)}.

**Lemma X [proved from cited results].** For every field and every h ∈ H ∖ (E ∪ {0}),
  Pr_𝔭[χ_𝔭(h) = −1] ≥ 1/(5d)  (N ≥ N₀(d)).
Consequently, with t_χ := ⌈5d·(ln |H|_bound + 2 ln log N)⌉ independent characters,
  Pr[H ∩ ⋂_j ker χ_j ⊄ E ∪ {0}] ≤ (log N)^{−2},
and on the complement H₀ := H ∩ ⋂_j ker χ_j is a subgroup of order ≤ 2.
- **Cited input.** Lagarias–Odlyzko (in *Algebraic Number Fields*, Durham 1975, Academic Press 1977, Thm 1.3).
  - Setting: L/K Galois with group G, C a conjugacy class, and x ≥ exp(10 n_L (log d_L)²).
  - Statement: |π_C(x) − (|C|/|G|)Li(x)| ≤ (|C|/|G|)Li(x^{β₀}) + c₁x·exp(−c₂ n_L^{−1/2}(log x)^{1/2}), with c₁, c₂
    absolute and effective.
  - The β₀ term is present only if ζ_L has a real zero β₀ ∈ [1 − (4 log d_L)^{−1}, 1).
- *Proof.*
  - **Setup.**
    - Apply LO to L_h/K, G = C₂, C = {σ}. Here n_L = 2d, d_L ≤ Δ̄, and log x ≥ 20d(log Δ̄)².
    - For 𝔭 ∤ 2D unramified, χ_𝔭(h) = −1 iff Frob_𝔭 = σ. The value is the Legendre symbol of a 𝔭-unit
      representative; for atomic products it is ∏_P ((a_P − b_P s)/r)^{e_P}.
  - **Lower bound on inert first-degree primes with norm in (x/2, x].**
    - Main term: ½(Li(x) − Li(x/2)) ≥ x/(4 log x).
    - The β₀ term, if present, has β₀ < 1 − ν₁ because h ∉ E. It is at most 2x^{1−ν₁} ≤ 2x·exp(−ν₁(log N)²), which
      is o(x/log x) since ν₁ ≥ c_d/log N.
    - The error term is 2c₁x·exp(−c₂ log N/√(2d)) = o(x/log x).
    - Primes of degree ≥ 2 number ≤ d√x, and ramified primes ≤ log₂ Δ̄. Both are negligible.
    - Total: ≥ x/(4 log x)·(1 − o(1)).
  - **Upper bound on all first-degree primes with norm in (x/2, x]:** ≤ d(π(x) − π(x/2)) ≤ 1.1·d·x/log x.
  - The ratio is ≥ 1/(5d).
  - **The kernel.** Independence gives Pr[h ∈ ⋂ ker χ_j] ≤ (1 − 1/(5d))^{t_χ} for each nonexceptional h ≠ 0. A
    union bound over |H| classes gives the stated probability.
  - **The group.** H₀ is a subgroup contained in E ∪ {0}, so |H₀| ≤ 2. ∎

**Cost.** t_χ = O_d(log N) characters. Each is sampled in (log N)^{O(1)} expected time and evaluated on all atoms in
M_rows·k·(log N)^{O(1)}.
- Under GRH, characters at primes of size (log N)^{O(1)} would do (cf. LV Claim 6.21). We do not use GRH.

### 3.4 Linear algebra

- The matrix has M_rows rows and C_col ≤ M_rows − 12 columns. Its density is w ≤ M_rows·(2k·log₂N + t_χ + 2)
  nonzero entries.
- A uniformly random left-kernel vector is obtained by Wiedemann's method:
  - It runs over F_{2^e} with 2^e ≥ M_rows⁴ and is applied to the transpose.
  - A singular-system solver with preconditioning is analysed rigorously in Wiedemann (IEEE Trans. IT 32 (1986)
    54–62) and Kaltofen–Saunders (AAECC-9, LNCS 539 (1991), Thms 2–4). LP 1992 §10 uses the same method.
  - "x − A⁻(Ax)" gives a uniform element of the kernel, given the solver's coins. Coordinate expansion over an F₂
    basis of F_{2^e} returns an F₂ kernel vector.
- Each kernel vector is verified, so the method is Las Vegas.
- Expected cost O(M_rows·w + M_rows²)·(log N)^{O(1)} = (M_rows)^{2}(log N)^{O(1)}.
- dim ker ≥ 12. So two uniform samples are independent with probability ≥ 1 − 2^{−11}, and the third test vector
  s₁ + s₂ costs nothing.
- **Why three tests suffice [proved].** φ(s₁), φ(s₂) ∈ H₀, and |H₀| ≤ 2. Either one of them is 0, or both equal the
  unique nonzero element, and then φ(s₁ + s₂) = 0. In every case some S ∈ {s₁, s₂, s₁ + s₂} has γ′_S ∈ K^{×2},
  and S ≠ ∅.

### 3.5 Square root in K: Algorithm SQ and Lemma SQ

**Algorithm SQ.**
- **Input.**
  - Γ ∈ Z[β], given as g′(β)² times a product of atoms z_P (with multiplicities).
  - An integer A with this property: if Γ = η² for some η ∈ Z[β], then ‖η‖_∞ ≤ A.
- **Steps.**
  1. Let p be the least prime > B with p ∤ D. It is found among the O(log N) primes after B, since ω(D) = O(log N).
     Then g mod p is squarefree, and Γ is a unit at p, because Norm(Γ) = ±disc(g)²·∏(t^{d−1}F(a_P, b_P))^{|e_P|}
     has only prime factors ≤ B or dividing D.
  2. Factor ḡ = φ₁⋯φ_τ over F_p (τ ≤ d), Las Vegas. In each F_p[y]/φ_j, test whether Γ is a square (Euler) and take
     a square root ρ_j (Tonelli–Shanks in F_{p^{deg φ_j}}, Las Vegas). If some component is a non-square, output
     NONSQUARE.
  3. For each of the 2^{τ−1} sign vectors ε with ε₁ = +1:
     - CRT the values ε_jρ_j to ρ_ε ∈ Z[β]/p.
     - Newton-lift the inverse square root, z ↦ z + z(1 − Γz²)/2 mod p^{2^i}, to the first precision p^e > 2A.
     - Put η_ε := Γz mod p^e with centred coefficients. If η_ε² = Γ exactly, output η_ε.
  4. Output NONSQUARE.

**Lemma SQ [proved].**
- (i) **Correctness.**
  - Any output η satisfies η² = Γ exactly.
  - If Γ ∈ K^{×2}, then Γ = η² with η ∈ Z[β], and SQ outputs ±η.
- (ii) **The bound A.** Let Mh := ⌊‖g‖₂⌋ + 1 (≥ Mahler(g) ≥ max_k max(1, |β_k|)) and g_d := d‖g‖₁Mh^{d−1}. Then
  A := d·⌈U^{1/2}⌉·2^{d−1}Mh·g_d^{d−1}, U := g_d²·∏_P (‖z_P‖₁Mh^{d−1})^{|e_P|},
  is valid. log₂A ≤ ½·Σ_P|e_P|·log₂(‖z_P‖₁Mh^{d−1}) + O(d² log(d‖g‖₁Mh)) = O_d(Σ|e_P|·log N).
- (iii) **Cost.** The cost is O_d(𝖬(log A)·log log A) + (log N)^{O(1)} bit operations, where 𝖬 is integer
  multiplication time. This is quasi-linear in the size of Γ.

*Proof.*
- **(i).**
  - Write Γ = g′(β)²·γ with γ := ∏ z_P^{|e_P|} ∈ Z[β].
  - If Γ = δ², δ ∈ K, then (δ/g′(β))² = γ ∈ O_K, so δ/g′(β) ∈ O_K.
  - g′(β)O_K ⊂ Z[β]: indeed O_K ⊂ O_K^∨ ⊂ Z[β]^∨ = g′(β)^{−1}Z[β] (Euler's lemma; Z[β] ⊂ O_K). Hence
    η := δ ∈ Z[β].
  - p ∤ disc(g), so Z_p[β] ≅ ∏_j Z_p[y]/(G_j) is a product of unramified local rings, with G_j the Hensel lifts of
    the φ_j.
  - η mod p equals ±ρ_ε for exactly one class ε with ε₁ = +1.
  - For p odd and Γ a unit, a unit square root with prescribed residue lifts uniquely (Hensel). So
    η_ε ≡ ±η (mod p^e). Since p^e > 2A ≥ 2‖η‖_∞, centring recovers ±η exactly.
  - Conversely, the check η_ε² = Γ is exact. NONSQUARE at step 2 is correct because η ∈ O_K would reduce to a
    square in O_K/𝔓 = F_p[y]/φ_j (p is prime to the index).
- **(ii).**
  - |σ_k(z)| ≤ ‖z‖₁Mh^{d−1} and |g′(β_k)| ≤ g_d, so |σ_k(η)|² ≤ U.
  - Coefficients come from the inverse Vandermonde: (V^{−1})_{jk} = [y^j](g(y)/(y − β_k))/g′(β_k).
    - The numerator coefficients are ≤ 2^{d−1}Mh.
    - ∏_k |g′(β_k)| = |disc g| ≥ 1 gives |g′(β_k)| ≥ g_d^{−(d−1)}.
- **(iii).**
  - Γ is computed by a product tree.
  - Newton doubles the precision. Each step is O(d²) multiplications of ≤ log₂(p^e) = O(log A) bits.
  - The sign enumeration multiplies by 2^{d−1} = O(1). ∎
- (E1: all 92 NONSQUARE outputs were independently certified, and all 116 squares were verified. The ratio
  log₂‖η‖/log₂A was 0.08–0.22, so A is conservative but polynomially tight.)
- **Prior methods.**
  - Couveignes (LNM 1554, 1993) and Thomé (2012) need a prime inert in K, or a sign fixed by the norm for odd d.
  - Montgomery's and Nguyen's methods are lattice-based and heuristic.
  - The 2^{τ−1}-branch enumeration avoids both issues, and is affordable because d is fixed.
  - LV cite Montgomery "in time O(M²)" without a rigorous analysis. Our Lemma SQ replaces that step.

### 3.6 The congruence and ψ

- Let S be the chosen dependency, e_P := e_P(S), E := E(S) (even), and γ⁻ := ∏_{e_P<0} z_P^{|e_P|}.
- ψ is a ring homomorphism on Z[β], so ψ(η)² = ψ(Γ_S) = ψ(g′(β))²·∏_P ψ(z_P)^{|e_P|}.
- Hence, with Y as in step 5,
  Y² = q_S·t^{−E}·∏_P ψ(z_P)^{|e_P|}/ψ(γ⁻)² = q_S·t^{−E}·∏_P ψ(z_P)^{e_P} ≡ ∏_{i∈S} r_i² (mod N) by (1.1).
- **Units.**
  - ψ(z_P) = tG(P)/g₁ is a unit: G(P) = σℓn, with n B-smooth, and ℓ ∤ N, otherwise a factor was found.
  - ψ(g′(β)) = t^{d−2}f′(g₀/g₁) is a unit by the field check (§4).
  - √q_S involves only g₁ and primes ≤ B.
- No map on K^× is needed. Every identity is a ring identity in Z[β] or in Q, followed by reduction mod N.
- E1 checks this exactly: 0 violations in 116 cases.
- **Final rational square root.** √q_S is computed from the recorded exponent vectors (halving), with
  M_rows·k·(log N)^{O(1)} operations. No square root mod N is ever extracted.

### 3.7 Transcript invariance (T)

Theorem DD needs the exposed transcript 𝒯 to be invariant under r_i ↦ ε_ir_i. Every quantity computed by LV6 is
a function of (field, R_i, splits, atoms, certificates) and fresh coins:
- the columns, the characters (r_j, s_j), the kernel vectors, the Hensel prime p, the branch ε, η, the choice
  among s₁, s₂, s₁ + s₂, and Y.

Field rejections (§4) depend on the field only. None of these reads r_i, a Jacobi symbol of r_i, or any gcd
outcome. So (T) holds, and DD gives Pr[nontrivial | 𝒯] = 1 − 2^{1−ω}.
- **Retry policy.**
  - After a trivial gcd we restart with fresh rows, as in lemmaM's remark R1.
  - (Optional refinement.) Reusing the rows with a new dependency S′ is also exact, provided S′ ∉ {0, S}:
    conditioned on the outcome for S, ∏_{S′}ε_i is still uniform on μ₂. But a policy that adapts the *choice* of
    S′ must keep S′ outside the span of all previously tested S. We do not need this.

### 3.8 Cost and the five requirements (a)–(e)

| item | where | cost / status |
|---|---|---|
| (a) heights V = N^{1/(4k)}L^{O(1)} | only in log A = O(k·M_rows·log N) | Lemma SQ quasi-linear; no box enumeration in post-processing |
| (b) non-monic (f, g₁x − g₀) | β = f_dα, monic g, ψ(β) = f_d g₀/g₁; bad primes in H | Lemmas C, H, SQ |
| (c) the g₁/f_d parity column | one column E_i mod 2 | §3.1 |
| (d) (rows)^{2+o(1)} | Wiedemann L^{2a+o(1)}; characters, SQ (3 calls), √q_S at L^{a+o(1)} | §3.3–§3.6 |
| (e) transcript invariance (T) | all steps are functions of certificates and coins | §3.7 |

## 4. Field-level events, Lemma IR, and what changes in pcest

The kept field of hypD Corollary K has law ∝ q_A, which amplifies the probability of a field event by at most
1/μ, with μ ≥ (γ/2)^{2k+1}L^{−o(1)} = L^{−o(1)}. We instead **reject** bad fields before any descent attempt. The
field law then becomes 𝒜 conditioned on the good event 𝒢. Since g_A ≤ 1, E[g_A | 𝒢] ≥ γ − Pr[𝒢ᶜ] ≥ γ/2. So D_k
survives with γ/2, and FM, Theorem R, Lemma Kern, Prop 4.2 and Theorem 4.3, all stated for 𝒜, need **no change**.

| event | bound | treatment |
|---|---|---|
| gcd(g₁, N) ≠ 1 | — | proper factor (0 < g₁ < N) |
| degree drop f_d = 0 | **impossible**: f_d ≡ a_d ≡ N g₀^{−d} ≢ 0 (mod g₁) (codex) | — |
| gcd(f_d, N) ≠ 1 | — | proper factor (0 < \|f_d\| < N) |
| content(f) > 1 | content divides F(g₀, g₁) = N and is < N | proper factor |
| f reducible | ≤ exp(−log H/(13d log log H)) (Lemma IR) | reject; a factorization F₁F₂ = N with both ≠ ±1 is output |
| N \| f′(g₀/g₁) | ≤ 1/(2H+1) | reject (a proper gcd is a factor) |
| bad primes, index, units, Cl[2] | no rarity needed | inside H (Lemma H) |
| exceptional zero | ≤ 1 class, every field | Lemma E, three tests |
| characters miss | ≤ (log N)^{−2} per attempt | resample |
| Wiedemann / finite-field roots | Las Vegas | retry |

**N | f′(m) [proved].** Put m = g₀/g₁ mod N, so g₁m ≡ g₀.
- f′ = F̂′ + g₁R + (g₁x − g₀)R′, so f′(m) ≡ f̂′(m) + g₁R_c(m) (mod N).
- R_c(m) = Σ_i c_i m^{d−1−i} has coefficient 1 at c_{d−1}, and g₁ is a unit.
- So for fixed c₀..c_{d−2}, at most one value of c_{d−1} ∈ [−H, H] (since 2H + 1 < N) gives N | f′(m). ∎

**Lemma IR (irreducibility) [proved; elementary].** Fix g with gcd(g₁, N) = 1 and N free of primes ≤ B. For c
uniform in C_H,
  Pr[f reducible over Q] ≤ e·(1 − 1/(4d))^{|𝒫|},
where 𝒫 is the set of primes in (r₀(d), ½ log H]. Hence the bound is ≤ exp(−log H/(13 d log log H)) for large N.
*Proof.*
- Let Q := ∏_{r∈𝒫} r ≤ H^{1/2+o(1)} ≤ 2H/d.
- **Near-uniformity.** Each residue class of c mod Q has between ⌊(2H+1)/Q⌋^d and ⌈(2H+1)/Q⌉^d elements. So the law
  of c mod Q is ≤ (1 + Q/(2H+1))^d ≤ e times uniform.
- **The map mod r.** For r ∈ 𝒫 (r ∤ N since r < B), c ↦ f mod r is an affine bijection from F_r^d onto:
  - if r ∤ g₁: all polynomials of degree ≤ d with f(m_r) = v_r := N g₁^{−d} ≠ 0, where m_r = g₀/g₁;
  - if r | g₁: all polynomials with leading coefficient a_d ≢ 0 and arbitrary lower coefficients. Here
    (g₁x − g₀) ≡ −g₀ is a unit.
- **Counting irreducibles.**
  - Monic irreducible Q of degree d with Q(0) = w ≠ 0 correspond to elements of F_{r^d} of degree d and norm
    (−1)^d w. There are ≥ ((r^d − 1)/(r − 1) − 2r^{d/2})/d ≥ r^{d−1}/(2d) of them for r ≥ 16.
  - Translating x ↦ x + m_r and summing over the r − 1 leading coefficients, the irreducible fraction of the affine
    space is ≥ 1/(4d). In the case r | g₁ it is ≥ 1/(2d).
- **Conclusion.** An irreducible reduction of degree d, with r ∤ f_d, forces f irreducible over Q. By CRT the
  reductions are independent under the uniform law. Finally |𝒫| ≥ log H/(3 log log H). ∎
- **Comparison.**
  - pcest quoted LV's H^{−1/3+o(1)}, via a transfer of LV Lemma 5.1, checked at the level of statements. Codex
    sketches the same order, O(H^{−1/3} log H), by a variance argument.
  - Lemma IR is weaker but fully proved, and it is enough: it is L^{−ω(1)}, while D_k only needs Pr[𝒢ᶜ] ≤ γ/2.

**What LV6 changes in pcest's theorems.**
- **Nothing in FM, Theorem R, Lemma Kern, Prop 4.2 or Theorem 4.3.**
- In Theorem A's assembly:
  - (1) reject bad fields (table above); this halves γ at worst;
  - (2) Lemma IR replaces the LV 5.1 transfer caveat;
  - (3) the row count is M_rows = C_col + 12 with C_col ≤ (d+1)π(B) + O(log N);
  - (4) FM's search uses batch smoothness (§6.2), and FM's count N_ℓ(A) is then recognized exactly, not with
    probability 1 − N^{−2}.

## 5. Where LV §6 is used and where it is not

| LV §6 item | LV (d ≍ (log N)^{1/3}, heights L(1/3)) | here (fixed d, height N^{1/(4k)}) |
|---|---|---|
| |H| bound | Claim 6.5 via BLP Prop 7.4, √\|Δ_g\|(log n)^{O(d)} | Lemma H, 2^{dω(D)+d}h_K = N^{O(1)} |
| characters | Lemma 6.6: L(1/3)-size primes, sampler Υ, probability ≥ 1/(2d) | Lemma X: exp(O((log N)²))-size first-degree primes, uniform sampler, ≥ 1/(5d) |
| Siegel zeros | Lemma 6.10: ≤ (4/3) log log n classes, except an L(2/3)^{−…} set of f | Lemma E: ≤ 1 class, every f |
| pigeonhole | ℓ = (4/3) log log n repetitions, (ℓ choose 2) tests | 3 tests |
| square root | Montgomery "O(M²)", not analysed | Lemma SQ, quasi-linear, rigorous |
| index / non-monic | f′(α)², sketched; "v ∈ Z[α]" for non-monic α | Z[β], β = f_dα, g′(β)O_K ⊂ Z[β]; bad primes in H |
| transcript invariance | not considered (LV need Conj. 7.1) | §3.7 |

Our proof uses LV §6's architecture: random large characters, Stark plus Landau, a pigeonhole on the residual
kernel. It is otherwise self-contained.

## 6. The HECM caveat

### 6.1 What Theorem A needs
- **Smoothness tests occur only inside special-q searches.** These test G(P)/ℓ and F_A(P) for the ≤ C_abs·S²
  points of one lattice box. They also occur in the preprocessing, which is trial division.
- **The shape test** factors u₀, v₀ ≈ √N completely, using LP 1992. LV6 needs no smoothness tests.
- FM counts *every* doubly smooth candidate. A test that recognizes only "most" smooth numbers would require
  redoing FM (and LD) for a restricted set.

### 6.2 Route 1: batch smoothness (Prop BS) [proved]

**Lemma BS (Bernstein's smooth parts) [proved; elementary].** Let 𝒫 = ∏_{p≤B} p, let x ≥ 2, and let 2^e ≥ log₂x.
- The B-smooth part of x is gcd(x, (𝒫 mod x)^{2^e} mod x).
- Given x₁, …, x_n, all smooth parts can be computed in O(𝖬(b) log b) bit operations, b = log₂𝒫 + Σ log₂x_i.
- *Proof.*
  - Write x = s·u, with s B-smooth and gcd(u, 𝒫) = 1. Every p^j ∥ s has j ≤ log₂x ≤ 2^e, so s | 𝒫^{2^e}, while
    gcd(u, 𝒫^{2^e}) = 1.
  - Cost: a product tree of the x_i, one reduction of 𝒫 modulo their product, and a remainder tree. Each level
    costs O(𝖬(b)) (standard simultaneous reduction: von zur Gathen–Gerhard, *Modern Computer Algebra*, §10.1–10.3).
    Then e = O(log log x) squarings mod x_i, and a gcd.
  - This is Bernstein, "How to find smooth parts of integers" (2004), Algorithm 2.1 / Thms 2.2–2.3. We use only
    the elementary version proved here. ∎

**Prop BS [proved].** Replace the per-candidate HECM test in the special-q search (pcest §1) as follows:
1. enumerate the whole box (≤ C_abs·S² points, or stop at once if the lattice is irregular, Lemma 1.3);
2. batch-test all values G(P)/ℓ and F_A(P);
3. output the first candidate, in a fixed order, that passes the filters and is doubly smooth;
4. trial-divide that candidate's two values by the primes ≤ B.

Properties:
- The search is deterministic given (A, ℓ). It recognizes every counted candidate, so N_ℓ(A) is exactly pcest's
  count, and FM, R, LD and 4.3 apply unchanged.
- Its cost is Õ(B + S²·log N).
- It is a function of (A, ℓ), so (T) is unaffected.

**Cost in Theorem A.** Per row, the cost changes from (T_shape + T_SQ) to (T_shape + T_SQ + B)·L^{o(1)}. Since
#rows·B = B^{2+o(1)} is the linear-algebra cost already present, **the exponent never changes**. At Theorem A's
optimum, a = θ/√2 < 1/√2 and T_SQ = L^{1/√2+o(1)} ≥ B in any case.

### 6.3 Route 2: unconditional HECM (Prop HS) [proved modulo genus2rh Props A/B (internal)]

**Ingredients.**
- **LPP I, Thm 1.3** (Phil. Trans. R. Soc. A 345 (1993) 397–408; statement checked against the page image in
  `codex/scratch/genus2rh_1/lpp1_theorem13.png`): with c₁ = 1980000^{1/3}, x ≥ x₃ (effective),
  z = L_x[2/3, c₁] and x^{3/5} ≤ y ≤ x, the number of z-smooth integers in [x, x + y] is
  ≥ y·exp(−(log x)^{1/3}(log log x)^{8/3}).
- **LPP II, Thm 1.1** (Proc. LMS 84 (2002)), as quoted in genus2rh §2.
- **genus2rh**: Lemma 1, Lemmas 2–3, Prop A (separation at a second prime q), Lemma 4 and Prop B. These were
  reviewed by codex and an Opus critic, with no fatal error, and are unrefereed. They hold for any modulus n with
  two distinct prime factors p ≤ v and q > 8100; repeated factors are allowed.

**Prop HS.** There is a randomized algorithm SMOOTH(n, B) with these properties:
- It never declares a non-B-smooth n smooth: it outputs a factorization certified by AKS.
- It factors every B-smooth n ≤ N^{O(1)} with probability ≥ 1/2.
- Its expected time is ≤ L_B[2/3, 2^{2/3}c₁ + o(1)]·(log N)^{O(1)} = L^{o(1)}.

*Proof (bookkeeping).*
1. **Small primes.** Trial division by primes ≤ 8100 (and ≤ x₃^{1/2}).
2. **Recursive splitting.** Maintain a list of unfactored parts. Remove perfect powers by integer roots and prime
   parts by AKS. On a composite part n′ that is not a prime power, run stages j = 14, …, ⌈log₂B⌉, with
   v_j = 2^j, w_j = ⌈L_{v_j²}[2/3, c₁]⌉ and
   h_j = ⌈κ·exp((log v_j²)^{1/3}(log log v_j²)^{8/3})·λ(v_j)·λ(N^{O(1)})·log(8 log₂N)⌉
   curves of Algorithm C with parameters (v_j, w_j).
3. **One adequate stage succeeds.** If n′ is B-smooth, its least prime p ≤ B. Let j* be the first stage with
   v_{j*} ≥ p, and let q be another prime factor of n′.
   - Take x = p² − ½p^{3/2} and y = p^{3/2} ≥ x^{3/5}. Since w_{j*} ≥ L_x[2/3, c₁], LPP I Thm 1.3 gives
     S_w(p) ≥ p^{3/2}exp(−(log p²)^{1/3}(log log p²)^{8/3}) ≥ 2p.
   - Then Prop B gives a per-curve success probability π* ≥ exp(−(log p²)^{1/3}(log log p²)^{8/3})/(C·λ(p)λ(q)).
   - With h_{j*}π* ≥ log(8 log₂N), stage j* fails with probability ≤ 1/(8 log₂N).
   - Primes p ≤ max(8100, x₃^{1/2}) were removed in step 1.
4. **All splits.** There are ≤ log₂n ≤ O(log N) splits. A union bound gives overall failure ≤ 1/8 + (Las Vegas
   overhead).
5. **Cost.** Stage j costs h_j·w_j·O(log v_j)(log N)² (genus2rh §3.1). The largest stage dominates, giving
   L_{B²}[2/3, c₁]·exp(O((log B)^{1/3+o(1)}))·(log N)^{O(1)} = L_B[2/3, 2^{2/3}c₁ + o(1)](log N)^{O(1)}.
6. **Truncation.** Truncating at twice this expectation (Markov) yields the stated ≥ 1/2. ∎

**Status.** This writes out genus2rh §6.3 ("V4, sketch") for the smoothness-test use. It closes the gap that
LPP's unpublished Part III left, for Theorem A's purposes. Its weak point is that genus2rh's Props A/B are internal
and unrefereed. Route 1 does not depend on them, so Theorem A uses Route 1.

### 6.4 Other routes, checked
- **Pomerance 1987 / Lenstra–Pomerance §7 (ECM on "recognizable" smooth numbers).**
  - LP Theorem 7.1 (from Lenstra 1987 Cor. 2.8) handles numbers with all prime factors in a set 𝒮 of primes with
    many smooth numbers in (p − √p, p + √p).
  - LP §6 shows that recognizable smooth numbers are a fair fraction of the smooth numbers needed *by LP's
    sampler*.
  - For Theorem A this would require FM (a lower bound in moving progressions) for 𝒮-restricted smooth numbers.
    Theorem R and LD restrict by positivity, so they are unaffected. Not needed.
- **Pollard–Strassen / Costa–Harvey per-candidate tests (A′):** B^{1/2} per test. Superseded by Prop BS.
- **Lenstra–Pomerance 1992 itself (the shape test).** It uses ECM with recognizable smooth numbers (LP §7,
  Theorem 7.1, Algorithm 7.2; checked from the paper) and Wiedemann's linear algebra. **It does not use HECM.** So
  Theorem A′'s "removal of HECM" was correct, but it was not needed.

### 6.5 Verdict
- **Theorem A, at c = (1+θ)/√2, is the honest conditional headline, with no HECM assumption.** Its infimum is
  0.9192 at d = 9, approached as k → ∞ and not attained at any finite k.
- Theorem A′'s 0.967 is superseded.
- LV6 needs no smoothness test.

## 7. Theorem A, assembled (as it would stand if LD were proved)

**Theorem A [conditional on Hypothesis LD only].**
- **Parameters.**
  - Fix odd d ≥ 7, an integer k > (d+1)/4 and η₀ > 0.
  - Put s = 1/(d+1), h = max(1/(2k), s − (d−2)/(4k)) + η₀ and θ = 2s + (d−1)/(4k) + h.
  - Assume θ < 1. Let a = θ/√2, B = L_N[1/2, a], and S² = (log N)^{C₃}/(ρ_rρ_F).
- **Hypothesis.** Assume pcest's Hypothesis LD for design (a″_k) at these parameters, uniformly in N.
- **Conclusion.**
  - There is a Las Vegas algorithm that splits every composite N that is not a prime power. Its expected time is
    L_N[1/2, (1+θ)/√2 + o(1)].
  - Recursion (with AKS) factors N completely within the same bound.
  - The infimum over k of the constant is (1 + 3/(d+1))/√2: d = 7: 0.9723; d = 9: 0.9192; d = 11: 0.8839;
    d = 19: 0.8132.
  - c < 1 needs d ≥ 7. It needs k ≥ 13 at d = 7 and k ≥ 12 at d = 9.

**Checklist (every ingredient, with its status).**

| # | ingredient | status | source |
|---|---|---|---|
| 1 | preprocessing: trial division to B, perfect powers, AKS | standard | — |
| 2 | field law (a″_k), Lemma GB, Lemma 1.2, Lemma 1.3 | proved (internal), both critics | pcest §1 |
| 3 | field rejections: content, N \| f′(m), reducibility (Lemma IR); γ → γ/2 | proved here | §4 |
| 4 | Theorem FM (uses Harper, LV §8, Tenenbaum, Hildebrand) | proved (internal) + cited | pcest §2 |
| 5 | Theorem R (Lemmas Q, DS) | proved (internal) | pcest §3 |
| 6 | Lemma Kern, Prop 4.2 | proved (internal) | pcest §4.1–4.2 |
| 7 | **Hypothesis LD** | **OPEN** | pcest §4.2; agent `ldbound` |
| 8 | Theorem 4.3: FM + R + 4.2 + LD ⟹ PC ≤ L^{o(1)} ⟹ D_k(γ) | proved implication (internal) | pcest §4.2 |
| 9 | Theorems S, S′, K(i)–(iii), Corollary K field selection | proved (internal) | hypD §2–§3 |
| 10 | shape test: Lenstra–Pomerance 1992 (ECM + Pomerance 1987; no HECM) | published | LP 1992 Thm, §7, §10 |
| 11 | special-q smoothness: Prop BS (batch; Lemma BS) | proved here | §6.2 |
| 12 | Theorem DD with transcript invariance (T) | proved (internal) | lemmaM §1 |
| 13 | **LV6**: Prop LV6 (Lemmas C, H, E, X, SQ; Wiedemann) | proved here, from cited LO 1977 Thm 1.3, Stark 1974 Lemma 8, Landau, Lenstra 1992 Thm 6.5, Wiedemann 1986 / Kaltofen–Saunders 1991 | §2–§3 |
| 14 | cost assembly | proved arithmetic | below |

**Cost assembly.**
- Rows: M_rows = (d+1)π(B) + O(log N) = B·L^{o(1)}.
- Per row (Theorem K(iii), Corollary K): (T_shape + T_SQ + B)·L^{o(1)}, with T_shape = L[1/2, 1/√2] and
  T_SQ = L[1/2, θ/(2a)]. The field-selection overhead is L^{o(1)}, since μ ≥ L^{−o(1)}.
- LV6: L[1/2, 2a + o(1)].
- Total: L[1/2, a + max(1/√2, θ/(2a))] + L[1/2, 2a]. At a = θ/√2 this is L[1/2, (1+θ)/√2], because θ < 1.
- Splitting succeeds with probability ≥ 1/2 per run, so the expected number of runs is ≤ 2.

**Caveats.**
- All internal results are unrefereed.
- LD is open.
- The constant is an infimum over fixed k.
- Effective constants (Stark, Landau, LO c₁, c₂, LPP x₃) make "N ≥ N₀" effective but enormous.

## 8. Experiment E1 (pre-registered; exact; toy)

**Commands** (as registered):
```
.venv/bin/python -m src.lv6.e1 experiments/configs/lv6_e1.json experiments/raw/lv6_e1.jsonl
.venv/bin/python -m src.lv6.e1 --analyze experiments/raw/lv6_e1.jsonl > experiments/raw/lv6_e1_summary.json
```
The run took 216 s in total. Calibration used seeds 100 and 101 only (disclosed in the prereg).

| prediction | registered | result |
|---|---|---|
| P1 congruence | 0 violations, no ψ(g′) non-unit | 0 / 116 certified squares; 0 events — **holds** |
| P2 non-squares certified | 0 UNCERTIFIED | 0 / 92 (46 residue non-square at p, 46 character witness) — **holds** |
| P3 full kernel squares | ≥ 15/16 eligible; all 6 squares in all but ≤ 1 | 16/16 eligible; 16/16 all square — **holds** |
| P4 combination search | finds a square in every eligible instance | 16/16 (always the first vector) — **holds** |
| P5 characters needed | rank ≥ 1 in ≥ 15/16; parity-only square fraction ≤ ½ in ≥ 14/16 | rank 1–16 (16/16); fraction 0 in 14, ½ and ⅙ in 2 — **holds** |
| P6 (reported) | — | gcd nontrivial 57/116 (plain NFS congruences; no DD claim); log₂‖η‖/log₂A = 0.08–0.22; τ_p ∈ {1, 2, 3}, ≤ 4 sign tries; bad-prime incidences 0–10 158 per instance |

**Summary-file key names.** The summary keys differ from the prereg labels:
- `P3_all_full_square` is P3.
- `P4_combo_found` is P4.
- `P6_par_frac_le_half` is P5's second half.
- `P5_char_rank_le_bound` is a trivial bound check. The registered P5 rank condition is read from the per-instance
  `char_rank`: min = 1.

**Reading.**
- Bad primes (dividing f_d·disc f) are frequent in the toy fields: up to 10⁴ incidences.
- Giving them no columns did not hurt: H absorbs them, and the characters detected classes of rank up to 16.
- E1 cannot test Lemma E (no Siegel zeros are expected at any size) or any cost claim.
- **Factor discoveries.** There were 0 ACCIDENTAL gcd hits in relation collection; the null rate is tiny, since all
  factors are ≥ 2^{17} > B. The 57 nontrivial final gcds are what the plain NFS congruence gives, at ≈ ½ for
  semiprimes (heuristic). They are STRUCTURAL in the classical sense, not a DD statement.

## 9. Reviews

### 9.1 GPT-6.1-sol blind round (lv6_1): adjudication
- **Independent agreement.** It reached the same architecture without seeing our draft:
  - the bad-prime set S with H_S ⊂ 2^{|S|+d}h_K;
  - LO + Stark + Landau with explicit discriminant bounds;
  - Hensel lifting at a squarefree prime with ≤ 2^d branches and exact verification;
  - g′(β)O_K ⊂ Z[β] via the trace dual;
  - one parity column for both g₁ and f_d;
  - batch smoothness removing HECM at no cost to the exponent;
  - the genus2rh + LPP I Thm 1.3 completion of HECM;
  - LP 1992 being HECM-free;
  - headline A.
- **Taken from codex.**
  - |E| ≤ 1 (Landau, with all conductors ≤ N^{O(1)}), and the three-test pigeonhole.
  - Degree drop is impossible.
  - Content of f divides N.
  - The affine argument for N | f′(m).
  - The Kaltofen–Saunders citation.
- **Differences.**
  - Codex clears denominators with a norm multiplier D₀. We use Γ = g′(β)²γ⁺γ⁻. These are equivalent.
  - Codex's irreducibility bound is O(H^{−1/3} log H), sketched. We prove the weaker Lemma IR in full, which
    suffices.
  - Codex takes x = 2^{(log N)⁴}. We take log x = max(20d(log Δ̄)², (log N)²).
- Pending: adversarial round lv6_2 and the two Opus critics (§9.2–§9.4).

## 10. Candidate card

```text
ID: LV6-1 fixed-degree NFS post-processing for Dixon-descent rows (design a''_k).
Operation: columns = sign, rational primes <= B, ONE parity column (g1 and f_d), first-degree (r,s) with
  r !| f_d disc f; O(log N) quadratic characters at exp(O((log N)^2))-size first-degree primes; Wiedemann;
  three square tests by Hensel lifting at one unramified prime with 2^(tau-1) branches; Y via psi on Z[beta].
Proved: Lemmas C, H, E (<= 1 exceptional class), X, SQ; Prop LV6 (cost rows^2 polylog; transcript invariance).
Status: proved (internal, unrefereed); E1 toy: 0 violations, 0 uncertified.

ID: BS-1 batch smoothness inside complete special-q boxes (removes HECM from Theorem A).
Status: proved (elementary; Bernstein 2004). Cost O~(B + S^2 log N) per search; exponent unchanged.

ID: HS-1 unconditional HECM smoothness test (genus2rh Props A/B + LPP I Thm 1.3 + LPP II Thm 1.1).
Status: proved modulo genus2rh internal results; cost L_B[2/3, 2^(2/3) c1 + o(1)] (log N)^O(1).
```
