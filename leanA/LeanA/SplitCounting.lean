import LeanA.Defs

/-!
# (c) Split counting: [pp] Theorems thm:S, thm:Sprime and the combinatorial core of thm:K

`R` is uniform on the unit squares `Q_N` (`unitSquares N`). `SplitEvent N S R` is the event
`u ≡ R v (mod N)` for some `u, v ∈ S`.

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset Real

/-- MAIN. `[(ℤ/N)ˣ : Q_N] = |μ₂(ℤ/N)|` (squaring is a homomorphism with kernel `μ₂`), so for odd `N`
the index is `2^{ω(N)}`. -/
theorem card_unitSquares_mul (N : ℕ) [NeZero N] :
    (unitSquares N).card * (muTwo N).card = Fintype.card (ZMod N)ˣ := by
  set f : (ZMod N)ˣ →* (ZMod N)ˣ := powMonoidHom 2
  have h1 : (f.ker : Set (ZMod N)ˣ) = ↑(muTwo N) := by
    ext u; simp [muTwo, f]
  have h2 : (f.range : Set (ZMod N)ˣ) = ↑(unitSquares N) := by
    ext u; simp [unitSquares, f]
  have := Subgroup.card_ker_mul_card_range f
  rw [← Nat.card_eq_fintype_card, ← this, mul_comm]
  congr 1
  · rw [← SetLike.coe_sort_coe, h1, Nat.card_coe_set_eq, Set.ncard_coe_finset]
  · rw [← SetLike.coe_sort_coe, h2, Nat.card_coe_set_eq, Set.ncard_coe_finset]

/-! ## Square roots of unity: `|μ₂| ≤ 2^ω(N)` for odd `N` -/

lemma ordProj_dvd_sq_sub (N : ℕ) (hodd : Odd N) (p : ℕ) (hp : p ∈ N.primeFactors) (a : ℤ)
    (ha : (N : ℤ) ∣ a ^ 2 - 1) :
    ((p ^ N.factorization p : ℕ) : ℤ) ∣ a - 1 ∨ ((p ^ N.factorization p : ℕ) : ℤ) ∣ a + 1 := by
  have hpp : p.Prime := Nat.prime_of_mem_primeFactors hp
  have hpi : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hpp
  have hq : ((p ^ N.factorization p : ℕ) : ℤ) ∣ (a - 1) * (a + 1) := by
    have : ((p ^ N.factorization p : ℕ) : ℤ) ∣ (N : ℤ) :=
      Int.natCast_dvd_natCast.2 (Nat.ordProj_dvd N p)
    exact this.trans (by convert ha using 1; ring)
  have hp2 : p ≠ 2 := by
    rintro rfl
    have : 2 ∣ N := Nat.dvd_of_mem_primeFactors hp
    exact (Nat.not_even_iff_odd.2 hodd) (even_iff_two_dvd.2 this)
  push_cast at hq ⊢
  by_cases h1 : (p : ℤ) ∣ a - 1
  · left
    have h2 : ¬ (p : ℤ) ∣ a + 1 := by
      intro h2
      have : (p : ℤ) ∣ 2 := by
        have := dvd_sub h2 h1; convert this using 1; ring
      have : p ∣ 2 := by exact_mod_cast this
      exact hp2 ((Nat.prime_dvd_prime_iff_eq hpp Nat.prime_two).1 this)
    have hc : IsCoprime ((p : ℤ) ^ N.factorization p) (a + 1) :=
      ((Prime.coprime_iff_not_dvd hpi).2 h2).pow_left
    exact hc.dvd_of_dvd_mul_right hq
  · right
    have hc : IsCoprime ((p : ℤ) ^ N.factorization p) (a - 1) :=
      ((Prime.coprime_iff_not_dvd hpi).2 h1).pow_left
    exact hc.dvd_of_dvd_mul_left hq

lemma not_both (N : ℕ) (hodd : Odd N) (p : ℕ) (hp : p ∈ N.primeFactors) (a : ℤ) :
    ¬ (((p ^ N.factorization p : ℕ) : ℤ) ∣ a - 1 ∧ ((p ^ N.factorization p : ℕ) : ℤ) ∣ a + 1) := by
  rintro ⟨h1, h2⟩
  have hpp : p.Prime := Nat.prime_of_mem_primeFactors hp
  have hpq : (p : ℤ) ∣ ((p ^ N.factorization p : ℕ) : ℤ) := by
    have : 1 ≤ N.factorization p :=
      (Nat.Prime.dvd_iff_one_le_factorization hpp (Nat.mem_primeFactors.1 hp).2.2).1
        (Nat.dvd_of_mem_primeFactors hp) |> fun h => h
    exact_mod_cast dvd_pow_self p (by omega)
  have : (p : ℤ) ∣ 2 := by
    have := dvd_sub (hpq.trans h2) (hpq.trans h1); convert this using 1; ring
  have : p ∣ 2 := by exact_mod_cast this
  have hp2 : p = 2 := (Nat.prime_dvd_prime_iff_eq hpp Nat.prime_two).1 this
  subst hp2
  exact (Nat.not_even_iff_odd.2 hodd) (even_iff_two_dvd.2 (Nat.dvd_of_mem_primeFactors hp))

lemma card_muTwo_le (N : ℕ) [NeZero N] (hodd : Odd N) : (muTwo N).card ≤ 2 ^ omega N := by
  classical
  let A : (ZMod N)ˣ → ℤ := fun u => ((u : ZMod N).val : ℤ)
  have hA : ∀ u, ((A u : ℤ) : ZMod N) = (u : ZMod N) := by intro u; simp [A]
  let F : (ZMod N)ˣ → Finset ℕ := fun u =>
    N.primeFactors.filter (fun p => ((p ^ N.factorization p : ℕ) : ℤ) ∣ A u - 1)
  have hsq : ∀ u ∈ muTwo N, (N : ℤ) ∣ A u ^ 2 - 1 := by
    intro u hu
    simp only [muTwo, mem_filter, mem_univ, true_and] at hu
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    push_cast
    rw [hA, sub_eq_zero]
    have := congrArg (fun w : (ZMod N)ˣ => (w : ZMod N)) hu
    simpa using this
  have hinj : Set.InjOn F (muTwo N : Set (ZMod N)ˣ) := by
    intro u hu w hw huw
    have key : ∀ p ∈ N.primeFactors, ((p ^ N.factorization p : ℕ) : ℤ) ∣ A u - A w := by
      intro p hp
      have hmem := congrArg (fun T => p ∈ T) huw
      simp only [F, mem_filter, hp, true_and, eq_iff_iff] at hmem
      rcases ordProj_dvd_sq_sub N hodd p hp (A u) (hsq u hu) with h | h
      · have := dvd_sub h (hmem.1 h); convert this using 1; ring
      · rcases ordProj_dvd_sq_sub N hodd p hp (A w) (hsq w hw) with h' | h'
        · exact absurd ⟨hmem.2 h', h⟩ (not_both N hodd p hp (A u))
        · have := dvd_sub h h'; convert this using 1; ring
    have hN : (N : ℤ) ∣ A u - A w := by
      have hprod : ∏ p ∈ N.primeFactors, ((p ^ N.factorization p : ℕ) : ℤ) = (N : ℤ) := by
        rw [← Nat.cast_prod]
        congr 1
        exact (Nat.prod_primeFactors_pow_factorization (NeZero.ne N)).symm
      rw [← hprod]
      refine Finset.prod_dvd_of_coprime ?_ key
      intro p hp q hq hpq
      simp only [Function.onFun]
      rw [Nat.isCoprime_iff_coprime]
      exact Nat.Coprime.pow _ _ ((Nat.coprime_primes (Nat.prime_of_mem_primeFactors hp)
        (Nat.prime_of_mem_primeFactors hq)).2 hpq)
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd] at hN
    push_cast at hN
    rw [hA, hA, sub_eq_zero] at hN
    exact Units.ext hN
  calc (muTwo N).card ≤ N.primeFactors.powerset.card := by
        refine Finset.card_le_card_of_injOn F ?_ hinj
        intro u _
        simp only [F, coe_powerset, Set.mem_preimage, Set.mem_powerset_iff, coe_filter]
        intro x hx; exact hx.1
    _ = 2 ^ omega N := by rw [card_powerset]; rfl

/-! ## The moment method shared by Theorems S and S′ -/

/-- `n(R) = #{(u, v) ∈ S² : u ≡ R v}`. -/
def nR (N : ℕ) [NeZero N] (S : Finset ℕ) (R : (ZMod N)ˣ) : ℕ :=
  ((S ×ˢ S).filter (fun p => (p.1 : ZMod N) = (R : ZMod N) * (p.2 : ZMod N))).card

/-- The number of proportional pairs `((u, v), (u', v')) ∈ (S²)²`, `u v' = u' v`. -/
def quadCount (S : Finset ℕ) : ℕ :=
  (((S ×ˢ S) ×ˢ (S ×ˢ S)).filter (fun q => q.1.1 * q.2.2 = q.2.1 * q.1.2)).card

/-- The subgroup of unit squares. -/
def sqSub (N : ℕ) [NeZero N] : Subgroup (ZMod N)ˣ := (powMonoidHom 2 : (ZMod N)ˣ →* (ZMod N)ˣ).range

lemma mem_unitSquares_iff (N : ℕ) [NeZero N] (R : (ZMod N)ˣ) :
    R ∈ unitSquares N ↔ R ∈ sqSub N := by
  simp [unitSquares, sqSub]

lemma card_sqSub (N : ℕ) [NeZero N] : Nat.card (sqSub N) = (unitSquares N).card := by
  have h2 : ((sqSub N) : Set (ZMod N)ˣ) = ↑(unitSquares N) := by
    ext u; simp [unitSquares, sqSub]
  rw [← SetLike.coe_sort_coe, h2, Nat.card_coe_set_eq, Set.ncard_coe_finset]

/-- The class of `u` in `(ℤ/N)ˣ / Q_N` (junk value `1` if `u` is not coprime to `N`). -/
def sqClass (N : ℕ) [NeZero N] (u : ℕ) : (ZMod N)ˣ ⧸ sqSub N :=
  if h : u.Coprime N then QuotientGroup.mk (ZMod.unitOfCoprime u h) else 1

lemma sqClass_of_coprime (N : ℕ) [NeZero N] (u : ℕ) (h : u.Coprime N) :
    sqClass N u = QuotientGroup.mk (ZMod.unitOfCoprime u h) := by
  unfold sqClass
  split_ifs
  · rfl

lemma second_moment (N : ℕ) [NeZero N] (S : Finset ℕ) (hS : ∀ s ∈ S, s.Coprime N)
    (hprod : ∀ a ∈ S, ∀ b ∈ S, a * b < N) :
    ∑ R ∈ unitSquares N, (nR N S R) ^ 2 ≤ quadCount S := by
  classical
  set T := (S ×ˢ S) ×ˢ (S ×ˢ S)
  let c : (ZMod N)ˣ → ℕ × ℕ → Prop := fun R p => (p.1 : ZMod N) = (R : ZMod N) * (p.2 : ZMod N)
  have hsq : ∀ R, (nR N S R) ^ 2 = (T.filter (fun q => c R q.1 ∧ c R q.2)).card := by
    intro R
    rw [sq, nR, ← card_product, filter_product]
  simp_rw [hsq, card_filter]
  rw [Finset.sum_comm, quadCount, card_filter]
  refine Finset.sum_le_sum (fun q hq => ?_)
  rw [← card_filter]
  simp only [T, mem_product] at hq
  obtain ⟨⟨hu, hv⟩, hu', hv'⟩ := hq
  split_ifs with heq
  · refine Finset.card_le_one.2 (fun R hR R' hR' => ?_)
    simp only [mem_filter, c] at hR hR'
    have hvu : IsUnit (q.1.2 : ZMod N) := (ZMod.unitOfCoprime _ (hS _ hv)).isUnit
    apply Units.ext
    have := hR.2.1.symm.trans hR'.2.1
    exact (hvu.mul_left_inj).1 this
  · rw [Nat.le_zero, card_eq_zero, filter_eq_empty_iff]
    rintro R - ⟨h1, h2⟩
    apply heq
    have hz : ((q.1.1 * q.2.2 : ℕ) : ZMod N) = ((q.2.1 * q.1.2 : ℕ) : ZMod N) := by
      push_cast
      simp only [c] at h1 h2
      rw [h1, h2]; ring
    rw [ZMod.natCast_eq_natCast_iff', Nat.mod_eq_of_lt (hprod _ hu _ hv'),
      Nat.mod_eq_of_lt (hprod _ hu' _ hv)] at hz
    exact hz

lemma first_moment (N : ℕ) [NeZero N] (S : Finset ℕ) (hS : ∀ s ∈ S, s.Coprime N) :
    S.card ^ 2 ≤ Fintype.card ((ZMod N)ˣ ⧸ sqSub N) * ∑ R ∈ unitSquares N, nR N S R := by
  classical
  set G := (ZMod N)ˣ ⧸ sqSub N
  let fib : G → Finset ℕ := fun x => S.filter (fun u => sqClass N u = x)
  have h1 : S.card = ∑ x : G, (fib x).card :=
    card_eq_sum_card_fiberwise (fun _ _ => mem_univ _)
  have h2 : ((S ×ˢ S).filter (fun p => sqClass N p.1 = sqClass N p.2)).card =
      ∑ x : G, (fib x).card ^ 2 := by
    rw [card_eq_sum_card_fiberwise (f := fun p => sqClass N p.1) (t := univ)
      (fun _ _ => mem_univ _)]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    rw [sq, ← card_product]
    congr 1
    ext ⟨u, v⟩
    simp only [fib, mem_filter, mem_product]
    constructor
    · rintro ⟨⟨⟨hu, hv⟩, he⟩, hx⟩; exact ⟨⟨hu, hx⟩, hv, he ▸ hx⟩
    · rintro ⟨⟨hu, hx⟩, hv, hx'⟩; exact ⟨⟨⟨hu, hv⟩, hx.trans hx'.symm⟩, hx⟩
  have h3 : ((S ×ˢ S).filter (fun p => sqClass N p.1 = sqClass N p.2)).card ≤
      ∑ R ∈ unitSquares N, nR N S R := by
    simp_rw [nR, card_filter]
    rw [Finset.sum_comm]
    refine Finset.sum_le_sum (fun p hp => ?_)
    rw [← card_filter]
    split_ifs with he
    · rw [Nat.one_le_iff_ne_zero, ← pos_iff_ne_zero, card_pos]
      simp only [mem_product] at hp
      obtain ⟨hu, hv⟩ := hp
      set U := ZMod.unitOfCoprime p.1 (hS _ hu)
      set V := ZMod.unitOfCoprime p.2 (hS _ hv)
      have he' : (QuotientGroup.mk U : G) = QuotientGroup.mk V := by
        rw [sqClass_of_coprime N _ (hS _ hu), sqClass_of_coprime N _ (hS _ hv)] at he
        exact he
      rw [QuotientGroup.eq] at he'
      refine ⟨U * V⁻¹, ?_⟩
      simp only [mem_filter]
      constructor
      · rw [mem_unitSquares_iff]
        have := (sqSub N).inv_mem he'
        convert this using 1
        rw [mul_inv_rev, inv_inv, mul_comm]
      · have : (U : ZMod N) = ((U * V⁻¹ : (ZMod N)ˣ) : ZMod N) * (V : ZMod N) := by
          rw [← Units.val_mul, inv_mul_cancel_right]
        simpa [U, V] using this
    · exact Nat.zero_le _
  calc S.card ^ 2 = (∑ x : G, (fib x).card) ^ 2 := by rw [h1]
    _ ≤ (univ : Finset G).card * ∑ x : G, (fib x).card ^ 2 := sq_sum_le_card_mul_sum_sq
    _ = Fintype.card G * ((S ×ˢ S).filter (fun p => sqClass N p.1 = sqClass N p.2)).card := by
        rw [h2, card_univ]
    _ ≤ Fintype.card G * ∑ R ∈ unitSquares N, nR N S R := Nat.mul_le_mul_left _ h3

lemma card_quot_mul (N : ℕ) [NeZero N] :
    Fintype.card ((ZMod N)ˣ ⧸ sqSub N) * (unitSquares N).card = Nat.totient N := by
  rw [← card_sqSub, ← ZMod.card_units_eq_totient, ← Nat.card_eq_fintype_card,
    ← Nat.card_eq_fintype_card, ← Subgroup.index_eq_card, Subgroup.index_mul_card]

lemma unitSquares_card_pos (N : ℕ) [NeZero N] : 0 < (unitSquares N).card :=
  card_pos.2 ⟨1, mem_image.2 ⟨1, mem_univ _, one_pow _⟩⟩

lemma card_quot_le (N : ℕ) [NeZero N] (hodd : Odd N) :
    Fintype.card ((ZMod N)ˣ ⧸ sqSub N) ≤ 2 ^ omega N := by
  have h1 := card_quot_mul N
  have h2 := card_unitSquares_mul N
  rw [ZMod.card_units_eq_totient, ← h1, mul_comm] at h2
  have := Nat.eq_of_mul_eq_mul_right (unitSquares_card_pos N) h2
  rw [← this]
  exact card_muTwo_le N hodd

/-- The moment bound: `|S|⁴ ≤ 2^ω φ(N) · #{proportional pairs} · Pr[E_S]`. -/
theorem moment_bound (N : ℕ) [NeZero N] (hodd : Odd N) (S : Finset ℕ)
    (hS : ∀ s ∈ S, s.Coprime N) (hprod : ∀ a ∈ S, ∀ b ∈ S, a * b < N) :
    (S.card : ℝ) ^ 4 ≤ 2 ^ omega N * (Nat.totient N : ℝ) * (quadCount S : ℝ) *
      prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by
  classical
  set Q := unitSquares N
  set F := Q.filter (fun R : (ZMod N)ˣ => SplitEvent N S (R : ZMod N))
  set idx := Fintype.card ((ZMod N)ˣ ⧸ sqSub N)
  have hidx := card_quot_mul N
  have hidxle := card_quot_le N hodd
  have hQ := unitSquares_card_pos N
  have h1 := first_moment N S hS
  have h2 := second_moment N S hS hprod
  have hzero : ∀ R : (ZMod N)ˣ, R ∈ Q → R ∉ F → nR N S R = 0 := by
    intro R hR hRF
    rw [nR, card_eq_zero, filter_eq_empty_iff]
    rintro ⟨u, v⟩ huv h
    apply hRF
    simp only [mem_product] at huv
    exact mem_filter.2 ⟨hR, u, huv.1, v, huv.2, h⟩
  have hsumF : ∑ R ∈ Q, nR N S R = ∑ R ∈ F, nR N S R :=
    (Finset.sum_subset (filter_subset _ _) hzero).symm
  have hCS : (∑ R ∈ F, nR N S R) ^ 2 ≤ F.card * quadCount S := by
    calc (∑ R ∈ F, nR N S R) ^ 2 ≤ F.card * ∑ R ∈ F, nR N S R ^ 2 := sq_sum_le_card_mul_sum_sq
      _ ≤ F.card * ∑ R ∈ Q, nR N S R ^ 2 :=
          Nat.mul_le_mul_left _ (Finset.sum_le_sum_of_subset (filter_subset _ _))
      _ ≤ F.card * quadCount S := Nat.mul_le_mul_left _ h2
  have key : S.card ^ 4 ≤ idx ^ 2 * (F.card * quadCount S) := by
    calc S.card ^ 4 = (S.card ^ 2) ^ 2 := by ring
      _ ≤ (idx * ∑ R ∈ Q, nR N S R) ^ 2 := Nat.pow_le_pow_left h1 2
      _ = idx ^ 2 * (∑ R ∈ F, nR N S R) ^ 2 := by rw [hsumF]; ring
      _ ≤ idx ^ 2 * (F.card * quadCount S) := Nat.mul_le_mul_left _ hCS
  have key2 : S.card ^ 4 ≤ 2 ^ omega N * idx * (F.card * quadCount S) := by
    calc S.card ^ 4 ≤ idx ^ 2 * (F.card * quadCount S) := key
      _ = idx * idx * (F.card * quadCount S) := by ring
      _ ≤ 2 ^ omega N * idx * (F.card * quadCount S) :=
          Nat.mul_le_mul_right _ (Nat.mul_le_mul_right _ hidxle)
  have hprob : prob Q (fun R => SplitEvent N S (R : ZMod N)) = (F.card : ℝ) / Q.card := rfl
  rw [hprob, ← hidx]
  have hQr : (0 : ℝ) < Q.card := by exact_mod_cast hQ
  have : (2 : ℝ) ^ omega N * ((idx * Q.card : ℕ) : ℝ) * (quadCount S : ℝ) * ((F.card : ℝ) / Q.card)
      = ((2 ^ omega N * idx * (F.card * quadCount S) : ℕ) : ℝ) := by
    push_cast
    field_simp
  rw [this]
  exact_mod_cast key2

/-! ## Theorem S: the gcd bound on proportional pairs -/

lemma quadCount_eq (S : Finset ℕ) :
    quadCount S = ∑ p ∈ S ×ˢ S, ((S ×ˢ S).filter (fun p' => p.1 * p'.2 = p'.1 * p.2)).card := by
  rw [quadCount, card_filter, sum_product]
  refine sum_congr rfl (fun p _ => ?_)
  rw [card_filter]

/-- The partners of `(u, v)` in a window `[Y'/2, Y']`: at most `2 gcd(u, v)`. -/
lemma partners_le_two_gcd (S : Finset ℕ) (Y' : ℝ) (hS : ∀ s ∈ S, Y' / 2 ≤ (s : ℝ) ∧ (s : ℝ) ≤ Y')
    (hpos : ∀ s ∈ S, 0 < s) (u v : ℕ) (hu : u ∈ S) (hu0 : 0 < u) :
    ((S ×ˢ S).filter (fun p' => u * p'.2 = p'.1 * v)).card ≤ 2 * Nat.gcd u v := by
  classical
  set g := Nat.gcd u v
  have hg : 0 < g := Nat.gcd_pos_of_pos_left _ hu0
  set u₁ := u / g
  set v₁ := v / g
  have hu₁ : u₁ * g = u := Nat.div_mul_cancel (Nat.gcd_dvd_left u v)
  have hv₁ : v₁ * g = v := Nat.div_mul_cancel (Nat.gcd_dvd_right u v)
  have hcop : Nat.Coprime u₁ v₁ := Nat.coprime_div_gcd_div_gcd hg
  set M := ⌊Y'⌋₊
  have himg : ((S ×ˢ S).filter (fun p' => u * p'.2 = p'.1 * v)).image Prod.fst ⊆
      (Ioc 0 M).filter (fun x => u₁ ∣ x) := by
    intro x hx
    simp only [mem_image, mem_filter, mem_product] at hx
    obtain ⟨⟨a, b⟩, ⟨⟨ha, hb⟩, hab⟩, rfl⟩ := hx
    simp only at hab ⊢
    simp only [mem_filter, mem_Ioc]
    refine ⟨⟨hpos a ha, Nat.le_floor (hS a ha).2⟩, ?_⟩
    · have : u₁ * b = a * v₁ := by
        have h2 : (u₁ * b) * g = (a * v₁) * g := by
          calc (u₁ * b) * g = (u₁ * g) * b := by ring
            _ = a * (v₁ * g) := by rw [hu₁, hv₁, hab]
            _ = (a * v₁) * g := by ring
        exact Nat.eq_of_mul_eq_mul_right hg h2
      have : u₁ ∣ a * v₁ := ⟨b, this.symm⟩
      exact hcop.dvd_of_dvd_mul_right this
  have hinj : Set.InjOn Prod.fst
      (((S ×ˢ S).filter (fun p' => u * p'.2 = p'.1 * v)) : Set (ℕ × ℕ)) := by
    rintro ⟨a, b⟩ hab ⟨a', b'⟩ hab' (h : a = a')
    simp only [coe_filter, mem_product] at hab hab'
    subst h
    have : u * b = u * b' := by rw [hab.2, hab'.2]
    rw [Nat.eq_of_mul_eq_mul_left hu0 this]
  rw [← card_image_of_injOn hinj]
  refine (card_le_card himg).trans ?_
  rw [Nat.Ioc_filter_dvd_card_eq_div]
  apply Nat.div_le_of_le_mul
  have hu1pos : 0 < u₁ := Nat.div_pos (Nat.gcd_le_left v hu0) hg
  have : M ≤ 2 * u := by
    have h1 : (M : ℝ) ≤ Y' := Nat.floor_le (by linarith [(hS u hu).1, (hS u hu).2])
    have h2 : Y' ≤ 2 * (u : ℝ) := by linarith [(hS u hu).1]
    have : (M : ℝ) ≤ ((2 * u : ℕ) : ℝ) := by push_cast; linarith
    exact_mod_cast this
  calc M ≤ 2 * u := this
    _ = u₁ * (2 * g) := by rw [← hu₁]; ring

/-- `∑_{u, v ≤ M} gcd(u, v) ≤ M² (1 + log M)`. -/
lemma sum_gcd_le (M : ℕ) :
    ((∑ u ∈ Icc 1 M, ∑ v ∈ Icc 1 M, Nat.gcd u v : ℕ) : ℝ) ≤ (M : ℝ) ^ 2 * (1 + Real.log M) := by
  classical
  set I := Icc 1 M
  have h1 : ∀ u ∈ I, ∀ v ∈ I, Nat.gcd u v ≤ ∑ d ∈ I, if d ∣ u ∧ d ∣ v then d else 0 := by
    intro u hu v hv
    simp only [I, mem_Icc] at hu hv
    rw [← sum_filter]
    calc Nat.gcd u v = ∑ d ∈ (Nat.gcd u v).divisors, Nat.totient d := (Nat.sum_totient _).symm
      _ ≤ ∑ d ∈ (Nat.gcd u v).divisors, d := sum_le_sum (fun d _ => Nat.totient_le d)
      _ ≤ ∑ d ∈ I.filter (fun d => d ∣ u ∧ d ∣ v), d := by
          refine sum_le_sum_of_subset (fun d hd => ?_)
          have hg0 : Nat.gcd u v ≠ 0 := (Nat.gcd_pos_of_pos_left _ (by omega)).ne'
          rw [Nat.mem_divisors] at hd
          have hdu : d ∣ u := hd.1.trans (Nat.gcd_dvd_left u v)
          have hdv : d ∣ v := hd.1.trans (Nat.gcd_dvd_right u v)
          simp only [I, mem_filter, mem_Icc]
          refine ⟨⟨Nat.pos_of_dvd_of_pos hdu (by omega), (Nat.le_of_dvd (by omega) hdu).trans hu.2⟩,
            hdu, hdv⟩
  have h2 : ∑ u ∈ I, ∑ v ∈ I, Nat.gcd u v ≤ ∑ d ∈ I, d * (M / d) ^ 2 := by
    calc ∑ u ∈ I, ∑ v ∈ I, Nat.gcd u v
        ≤ ∑ u ∈ I, ∑ v ∈ I, ∑ d ∈ I, if d ∣ u ∧ d ∣ v then d else 0 :=
          sum_le_sum (fun u hu => sum_le_sum (fun v hv => h1 u hu v hv))
      _ = ∑ u ∈ I, ∑ d ∈ I, ∑ v ∈ I, if d ∣ u ∧ d ∣ v then d else 0 :=
          sum_congr rfl (fun u _ => sum_comm)
      _ = ∑ d ∈ I, ∑ u ∈ I, ∑ v ∈ I, if d ∣ u ∧ d ∣ v then d else 0 := sum_comm
      _ = ∑ d ∈ I, d * (M / d) ^ 2 := by
          refine sum_congr rfl (fun d _ => ?_)
          rw [← sum_product' (f := fun u v => if d ∣ u ∧ d ∣ v then d else 0), ← sum_filter,
            sum_const, smul_eq_mul]
          rw [show (I ×ˢ I).filter (fun x : ℕ × ℕ => d ∣ x.1 ∧ d ∣ x.2) =
              (I.filter (fun x => d ∣ x)) ×ˢ (I.filter (fun x => d ∣ x)) from by
                ext ⟨a, b⟩; simp only [mem_filter, mem_product]; tauto, card_product]
          have : I.filter (fun x => d ∣ x) = (Ioc 0 M).filter (fun x => d ∣ x) := by
            congr 1
          rw [this, Nat.Ioc_filter_dvd_card_eq_div]; ring
  have h3 : ((∑ d ∈ I, d * (M / d) ^ 2 : ℕ) : ℝ) ≤ (M : ℝ) ^ 2 * (harmonic M : ℝ) := by
    rw [harmonic_eq_sum_Icc]
    push_cast
    rw [mul_sum]
    refine sum_le_sum (fun d hd => ?_)
    have hd1 : (1 : ℝ) ≤ d := by
      simp only [I, mem_Icc] at hd; exact_mod_cast hd.1
    have hdiv : ((M / d : ℕ) : ℝ) ≤ (M : ℝ) / d := Nat.cast_div_le
    have h0 : (0 : ℝ) ≤ ((M / d : ℕ) : ℝ) := Nat.cast_nonneg _
    calc (d : ℝ) * ((M / d : ℕ) : ℝ) ^ 2 ≤ (d : ℝ) * ((M : ℝ) / d) ^ 2 := by gcongr
      _ = (M : ℝ) ^ 2 * (d : ℝ)⁻¹ := by field_simp
  have h4 := harmonic_le_one_add_log M
  have hM0 : (0 : ℝ) ≤ (M : ℝ) ^ 2 := by positivity
  calc ((∑ u ∈ I, ∑ v ∈ I, Nat.gcd u v : ℕ) : ℝ) ≤ ((∑ d ∈ I, d * (M / d) ^ 2 : ℕ) : ℝ) := by
        exact_mod_cast h2
    _ ≤ (M : ℝ) ^ 2 * (harmonic M : ℝ) := h3
    _ ≤ (M : ℝ) ^ 2 * (1 + Real.log M) := by gcongr

/-- MAIN. [pp] Theorem thm:S (split counting; ratio form of the Dixon–Pomerance lemma), probability
bound. -/
theorem split_counting (N : ℕ) [NeZero N] (hodd : Odd N) (Y' : ℝ) (hY : 2 ≤ Y')
    (hYN : 2 * Y' ^ 2 < N) (S : Finset ℕ)
    (hS : ∀ s ∈ S, Y' / 2 ≤ (s : ℝ) ∧ (s : ℝ) ≤ Y' ∧ Nat.Coprime s N) :
    (S.card : ℝ) ^ 4 / (2 * 2 ^ omega N * (Nat.totient N : ℝ) * Y' ^ 2 * (1 + Real.log Y')) ≤
      prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by
  have hpos : ∀ s ∈ S, 0 < s := fun s hs => by
    have := (hS s hs).1
    exact_mod_cast (by linarith : (0 : ℝ) < s)
  have hcop : ∀ s ∈ S, s.Coprime N := fun s hs => (hS s hs).2.2
  have hYN' : Y' ^ 2 < N := by nlinarith [sq_nonneg Y']
  have hprod : ∀ a ∈ S, ∀ b ∈ S, a * b < N := by
    intro a ha b hb
    have h1 := (hS a ha).2.1
    have h2 := (hS b hb).2.1
    have h0a : (0 : ℝ) ≤ a := Nat.cast_nonneg _
    have h0b : (0 : ℝ) ≤ b := Nat.cast_nonneg _
    have : ((a * b : ℕ) : ℝ) < N := by
      push_cast
      calc (a : ℝ) * b ≤ Y' * Y' := mul_le_mul h1 h2 h0b (by linarith)
        _ = Y' ^ 2 := by ring
        _ < N := hYN'
    exact_mod_cast this
  have hmom := moment_bound N hodd S hcop hprod
  set M := ⌊Y'⌋₊
  have hM : (M : ℝ) ≤ Y' := Nat.floor_le (by linarith)
  have hM2 : 2 ≤ M := Nat.le_floor (by exact_mod_cast hY)
  have hSI : S ⊆ Icc 1 M := fun s hs => mem_Icc.2 ⟨hpos s hs, Nat.le_floor (hS s hs).2.1⟩
  have hquad : (quadCount S : ℝ) ≤ 2 * Y' ^ 2 * (1 + Real.log Y') := by
    have h1 : quadCount S ≤ ∑ p ∈ S ×ˢ S, 2 * Nat.gcd p.1 p.2 := by
      rw [quadCount_eq]
      exact sum_le_sum (fun p hp => partners_le_two_gcd S Y'
        (fun s hs => ⟨(hS s hs).1, (hS s hs).2.1⟩) hpos p.1 p.2 (mem_product.1 hp).1
        (hpos _ (mem_product.1 hp).1))
    have h2 : ∑ p ∈ S ×ˢ S, 2 * Nat.gcd p.1 p.2 ≤
        2 * ∑ u ∈ Icc 1 M, ∑ v ∈ Icc 1 M, Nat.gcd u v := by
      rw [← mul_sum, sum_product]
      refine Nat.mul_le_mul_left _ ?_
      calc ∑ u ∈ S, ∑ v ∈ S, Nat.gcd u v ≤ ∑ u ∈ S, ∑ v ∈ Icc 1 M, Nat.gcd u v :=
            sum_le_sum (fun u _ => sum_le_sum_of_subset hSI)
        _ ≤ ∑ u ∈ Icc 1 M, ∑ v ∈ Icc 1 M, Nat.gcd u v := sum_le_sum_of_subset hSI
    have h3 := sum_gcd_le M
    have hM1 : (1 : ℝ) ≤ M := by exact_mod_cast (by omega : 1 ≤ M)
    have hlog : Real.log M ≤ Real.log Y' := Real.log_le_log (by linarith) hM
    have hlog0 : 0 ≤ Real.log M := Real.log_nonneg hM1
    have h4 : (M : ℝ) ^ 2 * (1 + Real.log M) ≤ Y' ^ 2 * (1 + Real.log Y') := by
      apply mul_le_mul _ (by linarith) (by linarith) (by positivity)
      exact pow_le_pow_left₀ (by linarith) hM 2
    have h12 : (quadCount S : ℝ) ≤ 2 * ((∑ u ∈ Icc 1 M, ∑ v ∈ Icc 1 M, Nat.gcd u v : ℕ) : ℝ) := by
      exact_mod_cast h1.trans h2
    linarith
  have hprob0 : 0 ≤ prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) :=
    div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _)
  have hphi : (0 : ℝ) < Nat.totient N := by
    exact_mod_cast Nat.totient_pos.2 (NeZero.pos N)
  have hlogY : 0 < 1 + Real.log Y' := by
    have := Real.log_nonneg (by linarith : (1 : ℝ) ≤ Y'); linarith
  rw [div_le_iff₀ (by positivity)]
  calc (S.card : ℝ) ^ 4 ≤ 2 ^ omega N * (Nat.totient N : ℝ) * (quadCount S : ℝ) *
        prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := hmom
    _ ≤ 2 ^ omega N * (Nat.totient N : ℝ) * (2 * Y' ^ 2 * (1 + Real.log Y')) *
        prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by gcongr
    _ = _ := by ring

/-- MAIN. [pp] Theorem thm:S, second part: if `2Y'² < N`, `u, v ∈ [1, Y']` and `u ≡ R v`, then
`(u, v)` is an integer multiple of every shortest nonzero vector of `Λ_R`. -/
theorem split_multiple_of_shortest (N : ℕ) (Y' : ℝ) (hYN : 2 * Y' ^ 2 < N) (R : ZMod N)
    (u v : ℕ) (hu : 1 ≤ u) (hv : 1 ≤ v) (huY : (u : ℝ) ≤ Y') (hvY : (v : ℝ) ≤ Y')
    (huv : (u : ZMod N) = R * (v : ZMod N)) (x y : ℤ) (hxy : IsShortestVec N R x y) :
    ∃ j : ℤ, (u : ℤ) = j * x ∧ (v : ℤ) = j * y := by
  obtain ⟨hlat, hne, hmin⟩ := hxy
  have hNpos : 0 < N := by
    have : (0 : ℝ) < N := by nlinarith [sq_nonneg Y']
    exact_mod_cast this
  have : NeZero N := ⟨hNpos.ne'⟩
  have huvlat : InSplitLattice N R (u : ℤ) (v : ℤ) := by
    simp only [InSplitLattice, Int.cast_natCast]; exact huv
  have hne' : ((u : ℤ), (v : ℤ)) ≠ (0, 0) := by
    intro h; simp only [Prod.mk.injEq] at h; omega
  have hmin1 := hmin u v huvlat hne'
  -- the determinant `u y - x v` vanishes
  have hdet : (u : ℤ) * y - x * v = 0 := by
    set D := (u : ℤ) * y - x * v
    have hD : (D : ZMod N) = 0 := by
      simp only [InSplitLattice] at hlat huvlat
      simp only [D]; push_cast
      rw [hlat, huv]; ring
    rw [ZMod.intCast_zmod_eq_zero_iff_dvd] at hD
    refine Int.eq_zero_of_abs_lt_dvd hD ?_
    refine abs_lt_of_sq_lt_sq ?_ (by positivity)
    have hCS : D ^ 2 ≤ ((u : ℤ) ^ 2 + (v : ℤ) ^ 2) * (x ^ 2 + y ^ 2) := by
      have : ((u : ℤ) ^ 2 + (v : ℤ) ^ 2) * (x ^ 2 + y ^ 2) - D ^ 2 = ((u : ℤ) * x + v * y) ^ 2 := by
        simp only [D]; ring
      nlinarith [sq_nonneg ((u : ℤ) * x + v * y)]
    have h2 : D ^ 2 ≤ ((u : ℤ) ^ 2 + (v : ℤ) ^ 2) ^ 2 := by
      refine hCS.trans ?_
      rw [sq ((u : ℤ) ^ 2 + (v : ℤ) ^ 2)]
      exact mul_le_mul_of_nonneg_left hmin1 (by positivity)
    have h3 : (((u : ℤ) ^ 2 + (v : ℤ) ^ 2 : ℤ) : ℝ) ≤ 2 * Y' ^ 2 := by
      push_cast
      have hu0 : (0 : ℝ) ≤ u := Nat.cast_nonneg _
      have hv0 : (0 : ℝ) ≤ v := Nat.cast_nonneg _
      nlinarith
    have h4 : ((D ^ 2 : ℤ) : ℝ) < ((N : ℤ) ^ 2 : ℤ) := by
      have h5 : ((D ^ 2 : ℤ) : ℝ) ≤ ((((u : ℤ) ^ 2 + (v : ℤ) ^ 2) ^ 2 : ℤ) : ℝ) := by
        exact_mod_cast h2
      have h6 : (0 : ℝ) ≤ (((u : ℤ) ^ 2 + (v : ℤ) ^ 2 : ℤ) : ℝ) := by positivity
      push_cast at h5 h6 h3 ⊢
      have : ((u : ℝ) ^ 2 + (v : ℝ) ^ 2) ^ 2 ≤ (2 * Y' ^ 2) ^ 2 := pow_le_pow_left₀ h6 h3 2
      have : (2 * Y' ^ 2) ^ 2 < (N : ℝ) ^ 2 := pow_lt_pow_left₀ hYN (by positivity) (by norm_num)
      linarith
    exact_mod_cast h4
  -- primitive part of `(x, y)`
  set g : ℤ := (Int.gcd x y : ℤ)
  have hg : 0 < g := by
    have : 0 < Int.gcd x y := Int.gcd_pos_iff.2 (by
      by_contra h; push Not at h; exact hne (by rw [h.1, h.2]))
    show (0 : ℤ) < (Int.gcd x y : ℤ)
    exact_mod_cast this
  set x₁ := x / g
  set y₁ := y / g
  have hx₁ : x₁ * g = x := Int.ediv_mul_cancel (Int.gcd_dvd_left x y)
  have hy₁ : y₁ * g = y := Int.ediv_mul_cancel (Int.gcd_dvd_right x y)
  have hcop : IsCoprime x₁ y₁ := by
    rw [Int.isCoprime_iff_gcd_eq_one]
    exact Int.gcd_div_gcd_div_gcd (Int.gcd_pos_iff.2 (by
      by_contra h; push Not at h; exact hne (by rw [h.1, h.2])))
  have h1 : (u : ℤ) * y₁ = x₁ * v := by
    have : ((u : ℤ) * y₁) * g = (x₁ * v) * g := by
      calc ((u : ℤ) * y₁) * g = (u : ℤ) * (y₁ * g) := by ring
        _ = x * v := by rw [hy₁]; linarith
        _ = (x₁ * v) * g := by rw [← hx₁]; ring
    exact mul_right_cancel₀ hg.ne' this
  obtain ⟨a, b, hab⟩ := hcop
  set m : ℤ := a * u + b * v
  have hmx : m * x₁ = u := by simp only [m]; linear_combination (u : ℤ) * hab - b * h1
  have hmy : m * y₁ = v := by simp only [m]; linear_combination (v : ℤ) * hab + a * h1
  -- `g ∣ m` by minimality
  set r := m % g
  set q := m / g
  have hmr : r + g * q = m := Int.emod_add_mul_ediv m g
  have hr0 : 0 ≤ r := Int.emod_nonneg m hg.ne'
  have hrg : r < g := Int.emod_lt_of_pos m hg
  have hr : r = 0 := by
    by_contra hr
    have hrpos : 0 < r := lt_of_le_of_ne hr0 (Ne.symm hr)
    have hlat' : InSplitLattice N R (r * x₁) (r * y₁) := by
      have ex : r * x₁ = (u : ℤ) - q * x := by rw [← hmx, ← hx₁, ← hmr]; ring
      have ey : r * y₁ = (v : ℤ) - q * y := by rw [← hmy, ← hy₁, ← hmr]; ring
      simp only [InSplitLattice] at hlat huvlat ⊢
      rw [ex, ey]; push_cast
      simp only [Int.cast_natCast] at huvlat
      rw [huvlat, hlat]; ring
    have hx1y1 : x₁ ^ 2 + y₁ ^ 2 > 0 := by
      by_contra h
      have hx0 : x₁ = 0 := by nlinarith [sq_nonneg x₁, sq_nonneg y₁]
      have hy0 : y₁ = 0 := by nlinarith [sq_nonneg x₁, sq_nonneg y₁]
      apply hne; rw [← hx₁, ← hy₁, hx0, hy0]; simp
    have hne'' : (r * x₁, r * y₁) ≠ (0, 0) := by
      intro h; simp only [Prod.mk.injEq, mul_eq_zero] at h
      rcases h with ⟨h1 | h1, h2 | h2⟩ <;> first | omega | (rw [h1, h2] at hx1y1; simp at hx1y1)
    have := hmin _ _ hlat' hne''
    rw [← hx₁, ← hy₁] at this
    have h7 : g ^ 2 * (x₁ ^ 2 + y₁ ^ 2) ≤ r ^ 2 * (x₁ ^ 2 + y₁ ^ 2) := by nlinarith
    have h8 : g ^ 2 ≤ r ^ 2 := le_of_mul_le_mul_right h7 hx1y1
    nlinarith
  refine ⟨q, ?_, ?_⟩
  · rw [← hmx, ← hx₁, ← hmr, hr]; ring
  · rw [← hmy, ← hy₁, ← hmr, hr]; ring

/-! ## Theorem S′: binomial estimates -/

lemma choose_ratio_aux (n k : ℕ) (hkn : k ≤ n) :
    ∀ i, i ≤ k → ((n.choose i * n.choose (k - i) : ℕ) : ℝ) * ((n - k : ℕ) : ℝ) ^ i ≤
      ((k.choose i * n.choose k : ℕ) : ℝ) * (n : ℝ) ^ i := by
  intro i
  induction i with
  | zero => intro _; simp
  | succ i ih =>
    intro hi
    have ih := ih (by omega)
    have e1 := Nat.choose_succ_right_eq n i
    have e2 := Nat.choose_succ_right_eq n (k - i - 1)
    have e3 := Nat.choose_succ_right_eq k i
    rw [show k - i - 1 + 1 = k - i by omega] at e2
    rw [show k - (i + 1) = k - i - 1 by omega]
    have e1' : ((n.choose (i + 1) : ℕ) : ℝ) * ((i : ℝ) + 1) = (n.choose i : ℝ) * ((n : ℝ) - i) := by
      have := congrArg (fun t : ℕ => (t : ℝ)) e1
      push_cast [Nat.cast_sub (show i ≤ n by omega)] at this
      exact this
    have e2' : ((n.choose (k - i) : ℕ) : ℝ) * ((k : ℝ) - i) =
        (n.choose (k - i - 1) : ℝ) * ((n : ℝ) - k + i + 1) := by
      have := congrArg (fun t : ℕ => (t : ℝ)) e2
      push_cast [Nat.cast_sub (show i ≤ k by omega), Nat.cast_sub (show k - i - 1 ≤ n by omega),
        Nat.cast_sub (show 1 ≤ k - i by omega)] at this
      rw [this]; ring
    have e3' : ((k.choose (i + 1) : ℕ) : ℝ) * ((i : ℝ) + 1) = (k.choose i : ℝ) * ((k : ℝ) - i) := by
      have := congrArg (fun t : ℕ => (t : ℝ)) e3
      push_cast [Nat.cast_sub (show i ≤ k by omega)] at this
      exact this
    push_cast [Nat.cast_sub hkn] at ih ⊢
    have hik : (i : ℝ) + 1 ≤ k := by exact_mod_cast hi
    have hkn' : (k : ℝ) ≤ n := by exact_mod_cast hkn
    have hpos : (0 : ℝ) < ((i : ℝ) + 1) * ((n : ℝ) - k + i + 1) := by positivity
    rw [← mul_le_mul_iff_of_pos_right hpos]
    have hnk : (0 : ℝ) ≤ (n : ℝ) - k := by linarith
    have hi0 : (0 : ℝ) ≤ i := Nat.cast_nonneg _
    -- rewrite both sides using the identities
    have L : (n.choose (i + 1) : ℝ) * (n.choose (k - i - 1) : ℝ) * ((n : ℝ) - k) ^ (i + 1) *
        (((i : ℝ) + 1) * ((n : ℝ) - k + i + 1)) =
        ((n.choose i : ℝ) * (n.choose (k - i) : ℝ) * ((n : ℝ) - k) ^ i) *
          (((n : ℝ) - i) * ((k : ℝ) - i) * ((n : ℝ) - k)) := by
      have : (n.choose (i + 1) : ℝ) * (n.choose (k - i - 1) : ℝ) * ((n : ℝ) - k) ^ (i + 1) *
          (((i : ℝ) + 1) * ((n : ℝ) - k + i + 1)) =
          ((n.choose (i + 1) : ℝ) * ((i : ℝ) + 1)) *
            ((n.choose (k - i - 1) : ℝ) * ((n : ℝ) - k + i + 1)) * ((n : ℝ) - k) ^ (i + 1) := by ring
      rw [this, e1', ← e2']; ring
    have R : (k.choose (i + 1) : ℝ) * (n.choose k : ℝ) * (n : ℝ) ^ (i + 1) *
        (((i : ℝ) + 1) * ((n : ℝ) - k + i + 1)) =
        ((k.choose i : ℝ) * (n.choose k : ℝ) * (n : ℝ) ^ i) *
          (((k : ℝ) - i) * (n : ℝ) * ((n : ℝ) - k + i + 1)) := by
      have : (k.choose (i + 1) : ℝ) * (n.choose k : ℝ) * (n : ℝ) ^ (i + 1) *
          (((i : ℝ) + 1) * ((n : ℝ) - k + i + 1)) =
          ((k.choose (i + 1) : ℝ) * ((i : ℝ) + 1)) * (n.choose k : ℝ) * (n : ℝ) ^ (i + 1) *
            ((n : ℝ) - k + i + 1) := by ring
      rw [this, e3']; ring
    rw [L, R]
    have hfac : ((n : ℝ) - i) * ((k : ℝ) - i) * ((n : ℝ) - k) ≤
        ((k : ℝ) - i) * (n : ℝ) * ((n : ℝ) - k + i + 1) := by
      have hki : (0 : ℝ) ≤ (k : ℝ) - i := by linarith
      have : ((n : ℝ) - i) * ((n : ℝ) - k) ≤ (n : ℝ) * ((n : ℝ) - k + i + 1) := by nlinarith
      nlinarith
    have hA0 : (0 : ℝ) ≤ (n.choose i : ℝ) * (n.choose (k - i) : ℝ) * ((n : ℝ) - k) ^ i := by
      positivity
    have hf0 : (0 : ℝ) ≤ ((n : ℝ) - i) * ((k : ℝ) - i) * ((n : ℝ) - k) := by
      have : (0 : ℝ) ≤ (n : ℝ) - i := by linarith
      have : (0 : ℝ) ≤ (k : ℝ) - i := by linarith
      positivity
    calc ((n.choose i : ℝ) * (n.choose (k - i) : ℝ) * ((n : ℝ) - k) ^ i) *
          (((n : ℝ) - i) * ((k : ℝ) - i) * ((n : ℝ) - k))
        ≤ ((k.choose i : ℝ) * (n.choose k : ℝ) * (n : ℝ) ^ i) *
          (((n : ℝ) - i) * ((k : ℝ) - i) * ((n : ℝ) - k)) := mul_le_mul_of_nonneg_right ih hf0
      _ ≤ ((k.choose i : ℝ) * (n.choose k : ℝ) * (n : ℝ) ^ i) *
          (((k : ℝ) - i) * (n : ℝ) * ((n : ℝ) - k + i + 1)) :=
          mul_le_mul_of_nonneg_left hfac (by positivity)

/-- `C(n, i) C(n, k − i) ≤ 2 C(k, i) C(n, k)` for `n ≥ 4k²`. -/
lemma choose_mul_choose_le (n k i : ℕ) (hk : 1 ≤ k) (hn : 4 * k ^ 2 ≤ n) (hi : i ≤ k) :
    n.choose i * n.choose (k - i) ≤ 2 * (k.choose i * n.choose k) := by
  have hkn : k ≤ n := by nlinarith
  have h1 := choose_ratio_aux n k hkn i hi
  have hn' : (4 : ℝ) * (k : ℝ) ^ 2 ≤ n := by exact_mod_cast hn
  have hk' : (1 : ℝ) ≤ k := by exact_mod_cast hk
  have hnpos : (0 : ℝ) < n := by nlinarith
  have hnk : ((n - k : ℕ) : ℝ) = (n : ℝ) - k := Nat.cast_sub hkn
  have hnkpos : (0 : ℝ) < (n : ℝ) - k := by nlinarith
  -- Bernoulli: `(1 - k/n)^i ≥ 1 - i k / n ≥ 3/4`
  have hb : (3 / 4 : ℝ) * (n : ℝ) ^ i ≤ ((n : ℝ) - k) ^ i := by
    have hB := one_add_mul_le_pow (show (-2 : ℝ) ≤ -(k : ℝ) / n by
      rw [le_div_iff₀ hnpos]; nlinarith) i
    have hik : (i : ℝ) ≤ k := by exact_mod_cast hi
    have h34 : (3 / 4 : ℝ) ≤ 1 + (i : ℝ) * (-(k : ℝ) / n) := by
      have e : (i : ℝ) * (-(k : ℝ) / n) = -((i : ℝ) * k / n) := by ring
      rw [e]
      have : (i : ℝ) * k / n ≤ 1 / 4 := by
        rw [div_le_iff₀ hnpos]
        have hk0 : (0 : ℝ) ≤ k := by linarith
        nlinarith
      linarith
    have : (1 + -(k : ℝ) / n) = ((n : ℝ) - k) / n := by field_simp; ring
    rw [this, div_pow] at hB
    have hnpow : (0 : ℝ) < (n : ℝ) ^ i := by positivity
    have := h34.trans hB
    rw [le_div_iff₀ hnpow] at this
    exact this
  rw [hnk] at h1
  have hpow : (0 : ℝ) < ((n : ℝ) - k) ^ i := by positivity
  have : ((n.choose i * n.choose (k - i) : ℕ) : ℝ) ≤ ((2 * (k.choose i * n.choose k) : ℕ) : ℝ) := by
    have h0 : (0 : ℝ) ≤ ((k.choose i * n.choose k : ℕ) : ℝ) := Nat.cast_nonneg _
    have h2 : ((n.choose i * n.choose (k - i) : ℕ) : ℝ) * ((n : ℝ) - k) ^ i ≤
        ((2 * (k.choose i * n.choose k) : ℕ) : ℝ) * ((n : ℝ) - k) ^ i := by
      push_cast at h0 h1 ⊢
      nlinarith
    exact le_of_mul_le_mul_right h2 hpow
  exact_mod_cast this


/-! ## Theorem S′: counting proportional pairs of `k`-products -/

lemma kProd_spec (D : Finset ℕ) (k : ℕ) (hD : ∀ p ∈ D, p.Prime) {s : ℕ} (hs : s ∈ kProducts D k) :
    s.primeFactors ∈ powersetCard k D ∧ ∏ p ∈ s.primeFactors, p = s ∧ 0 < s := by
  obtain ⟨T, hT, rfl⟩ := mem_image.1 hs
  have hTp : ∀ p ∈ T, p.Prime := fun p hp => hD p ((mem_powersetCard.1 hT).1 hp)
  rw [Nat.primeFactors_prod hTp]
  exact ⟨hT, rfl, prod_pos (fun p hp => (hTp p hp).pos)⟩

lemma kProd_partners (D : Finset ℕ) (k : ℕ) (hD : ∀ p ∈ D, p.Prime) (S : Finset ℕ)
    (hS : S ⊆ kProducts D k) (u v : ℕ) (hu : u ∈ S) (hv : v ∈ S) :
    ((S ×ˢ S).filter (fun p' => u * p'.2 = p'.1 * v)).card ≤
      D.card.choose (u.primeFactors ∩ v.primeFactors).card := by
  classical
  obtain ⟨huP, huprod, hu0⟩ := kProd_spec D k hD (hS hu)
  obtain ⟨hvP, hvprod, hv0⟩ := kProd_spec D k hD (hS hv)
  set C := u.primeFactors \ v.primeFactors
  set i := (u.primeFactors ∩ v.primeFactors).card
  have hCcard : C.card + i = k := by
    rw [card_sdiff_add_card_inter]; exact (mem_powersetCard.1 huP).2
  have hsub : ∀ p' ∈ (S ×ˢ S).filter (fun p' => u * p'.2 = p'.1 * v), C ⊆ p'.1.primeFactors := by
    intro p' hp' q hq
    simp only [mem_filter, mem_product] at hp'
    obtain ⟨⟨ha, hb⟩, hab⟩ := hp'
    have ha0 := (kProd_spec D k hD (hS ha)).2.2
    simp only [C, mem_sdiff] at hq
    have hqp : q.Prime := Nat.prime_of_mem_primeFactors hq.1
    have hqu : q ∣ u := Nat.dvd_of_mem_primeFactors hq.1
    have : q ∣ p'.1 * v := hab ▸ dvd_mul_of_dvd_left hqu _
    rcases (Nat.Prime.dvd_mul hqp).1 this with h | h
    · exact Nat.mem_primeFactors.2 ⟨hqp, h, ha0.ne'⟩
    · exact absurd (Nat.mem_primeFactors.2 ⟨hqp, h, hv0.ne'⟩) hq.2
  calc ((S ×ˢ S).filter (fun p' => u * p'.2 = p'.1 * v)).card ≤ (powersetCard i D).card := by
        refine card_le_card_of_injOn (fun p' => p'.1.primeFactors \ C) ?_ ?_
        · intro p' hp'
          rw [mem_coe] at hp' ⊢
          have hp'' := hp'
          simp only [mem_filter, mem_product] at hp''
          obtain ⟨hP, -, -⟩ := kProd_spec D k hD (hS hp''.1.1)
          rw [mem_powersetCard] at hP ⊢
          refine ⟨sdiff_subset.trans hP.1, ?_⟩
          rw [card_sdiff_of_subset (hsub p' hp'), hP.2]
          omega
        · intro p₁ hp₁ p₂ hp₂ heq
          rw [mem_coe] at hp₁ hp₂
          have e : p₁.1.primeFactors = p₂.1.primeFactors := by
            rw [← sdiff_union_of_subset (hsub p₁ hp₁), ← sdiff_union_of_subset (hsub p₂ hp₂)]
            simp only at heq
            rw [heq]
          have hp₁' := hp₁
          have hp₂' := hp₂
          simp only [mem_filter, mem_product] at hp₁' hp₂'
          have e1 : p₁.1 = p₂.1 := by
            rw [← (kProd_spec D k hD (hS hp₁'.1.1)).2.1, ← (kProd_spec D k hD (hS hp₂'.1.1)).2.1, e]
          have e2 : p₁.2 = p₂.2 := by
            have : u * p₁.2 = u * p₂.2 := by rw [hp₁'.2, hp₂'.2, e1]
            exact Nat.eq_of_mul_eq_mul_left hu0 this
          exact Prod.ext e1 e2
    _ = D.card.choose i := card_powersetCard _ _


lemma card_inter_fiber_le (D : Finset ℕ) (k : ℕ) (A : Finset ℕ) (hA : A ∈ powersetCard k D) (i : ℕ) :
    ((powersetCard k D).filter (fun B => (A ∩ B).card = i)).card ≤
      k.choose i * (D.card - k).choose (k - i) := by
  classical
  rw [mem_powersetCard] at hA
  have hDA : (D \ A).card = D.card - k := by rw [card_sdiff_of_subset hA.1, hA.2]
  calc ((powersetCard k D).filter (fun B => (A ∩ B).card = i)).card
      ≤ (powersetCard i A ×ˢ powersetCard (k - i) (D \ A)).card := by
        refine card_le_card_of_injOn (fun B => (A ∩ B, B \ A)) ?_ ?_
        · intro B hB
          rw [mem_coe, mem_filter, mem_powersetCard] at hB
          rw [mem_coe, mem_product, mem_powersetCard, mem_powersetCard]
          refine ⟨⟨inter_subset_left, hB.2⟩, ?_, ?_⟩
          · exact sdiff_subset_sdiff hB.1.1 (le_refl _)
          · show (B \ A).card = k - i
            have := card_sdiff_add_card_inter B A
            rw [inter_comm, hB.2, hB.1.2] at this
            omega
        · intro B₁ _ B₂ _ h
          simp only [Prod.mk.injEq] at h
          ext x
          constructor
          · intro hx
            by_cases hxA : x ∈ A
            · have : x ∈ A ∩ B₁ := mem_inter.2 ⟨hxA, hx⟩
              rw [h.1] at this; exact (mem_inter.1 this).2
            · have : x ∈ B₁ \ A := mem_sdiff.2 ⟨hx, hxA⟩
              rw [h.2] at this; exact (mem_sdiff.1 this).1
          · intro hx
            by_cases hxA : x ∈ A
            · have : x ∈ A ∩ B₂ := mem_inter.2 ⟨hxA, hx⟩
              rw [← h.1] at this; exact (mem_inter.1 this).2
            · have : x ∈ B₂ \ A := mem_sdiff.2 ⟨hx, hxA⟩
              rw [← h.2] at this; exact (mem_sdiff.1 this).1
    _ = k.choose i * (D.card - k).choose (k - i) := by
        rw [card_product, card_powersetCard, card_powersetCard, hA.2, hDA]

lemma sum_inter_choose_le (D : Finset ℕ) (k : ℕ) (hk : 1 ≤ k) (hn : 4 * k ^ 2 ≤ D.card)
    (A : Finset ℕ) (hA : A ∈ powersetCard k D) :
    ∑ B ∈ powersetCard k D, D.card.choose (A ∩ B).card ≤
      2 * (2 * k).choose k * D.card.choose k := by
  classical
  set n := D.card
  have hAk : A.card = k := (mem_powersetCard.1 hA).2
  have hmaps : ∀ B ∈ powersetCard k D, (A ∩ B).card ∈ range (k + 1) := by
    intro B _
    rw [mem_range]
    have := card_le_card (inter_subset_left : A ∩ B ⊆ A)
    omega
  rw [← sum_fiberwise_of_maps_to hmaps]
  calc ∑ i ∈ range (k + 1), ∑ B ∈ (powersetCard k D).filter (fun B => (A ∩ B).card = i),
        n.choose (A ∩ B).card
      = ∑ i ∈ range (k + 1), ((powersetCard k D).filter (fun B => (A ∩ B).card = i)).card *
          n.choose i := by
        refine sum_congr rfl (fun i _ => ?_)
        rw [sum_congr rfl (g := fun _ => n.choose i), sum_const, smul_eq_mul]
        intro B hB
        rw [(mem_filter.1 hB).2]
    _ ≤ ∑ i ∈ range (k + 1), k.choose i * (2 * (k.choose i * n.choose k)) := by
        refine sum_le_sum (fun i hi => ?_)
        have hik : i ≤ k := Nat.lt_succ_iff.1 (mem_range.1 hi)
        calc ((powersetCard k D).filter (fun B => (A ∩ B).card = i)).card * n.choose i
            ≤ (k.choose i * (n - k).choose (k - i)) * n.choose i :=
              Nat.mul_le_mul_right _ (card_inter_fiber_le D k A hA i)
          _ ≤ (k.choose i * n.choose (k - i)) * n.choose i :=
              Nat.mul_le_mul_right _ (Nat.mul_le_mul_left _ (Nat.choose_le_choose _ (by omega)))
          _ = k.choose i * (n.choose i * n.choose (k - i)) := by ring
          _ ≤ k.choose i * (2 * (k.choose i * n.choose k)) :=
              Nat.mul_le_mul_left _ (choose_mul_choose_le n k i hk hn hik)
    _ = 2 * n.choose k * ∑ i ∈ range (k + 1), k.choose i ^ 2 := by
        rw [mul_sum]; refine sum_congr rfl (fun i _ => ?_); ring
    _ = 2 * (2 * k).choose k * n.choose k := by rw [Nat.sum_range_choose_sq]; ring

lemma quadCount_kprod_le (D : Finset ℕ) (k : ℕ) (hk : 1 ≤ k) (hD : ∀ p ∈ D, p.Prime)
    (hn : 4 * k ^ 2 ≤ D.card) (S : Finset ℕ) (hS : S ⊆ kProducts D k) :
    quadCount S ≤ 2 * (2 * k).choose k * D.card.choose k ^ 2 := by
  classical
  set n := D.card
  set P := powersetCard k D
  set f : Finset ℕ → Finset ℕ → ℕ := fun A B => n.choose (A ∩ B).card
  have hinj : Set.InjOn Nat.primeFactors (S : Set ℕ) := by
    intro a ha b hb h
    rw [mem_coe] at ha hb
    rw [← (kProd_spec D k hD (hS ha)).2.1, ← (kProd_spec D k hD (hS hb)).2.1, h]
  have himg : S.image Nat.primeFactors ⊆ P := by
    intro A hA
    obtain ⟨s, hs, rfl⟩ := mem_image.1 hA
    exact (kProd_spec D k hD (hS hs)).1
  calc quadCount S
      = ∑ p ∈ S ×ˢ S, ((S ×ˢ S).filter (fun p' => p.1 * p'.2 = p'.1 * p.2)).card := quadCount_eq S
    _ ≤ ∑ p ∈ S ×ˢ S, f p.1.primeFactors p.2.primeFactors := by
        refine sum_le_sum (fun p hp => ?_)
        rw [mem_product] at hp
        exact kProd_partners D k hD S hS p.1 p.2 hp.1 hp.2
    _ = ∑ A ∈ S.image Nat.primeFactors, ∑ B ∈ S.image Nat.primeFactors, f A B := by
        rw [sum_product, sum_image hinj]
        refine sum_congr rfl (fun u _ => ?_)
        rw [sum_image hinj]
    _ ≤ ∑ A ∈ P, ∑ B ∈ P, f A B := by
        refine (sum_le_sum_of_subset himg).trans' ?_
        exact sum_le_sum (fun A _ => sum_le_sum_of_subset himg)
    _ ≤ ∑ A ∈ P, 2 * (2 * k).choose k * n.choose k :=
        sum_le_sum (fun A hA => sum_inter_choose_le D k hk hn A hA)
    _ = 2 * (2 * k).choose k * n.choose k ^ 2 := by
        rw [sum_const, card_powersetCard, smul_eq_mul]; ring


/-- MAIN. [pp] Theorem thm:Sprime (k-products). `D` is a set of primes in a dyadic interval
`(Z, 2Z]`, all coprime to `N`, with `|D| ≥ 4k²`; `S` is a set of products of exactly `k` distinct
primes of `D` in the window `[Y'/2, Y']`, `2Y'² < N`, and `|S| ≥ C(|D|, k)/(2k)`. -/
theorem split_counting_kprod (N : ℕ) [NeZero N] (hodd : Odd N) (k : ℕ) (hk : 1 ≤ k) (Z : ℕ)
    (D : Finset ℕ) (hD : ∀ p ∈ D, p.Prime ∧ Z < p ∧ p ≤ 2 * Z ∧ Nat.Coprime p N)
    (hDcard : 4 * k ^ 2 ≤ D.card) (Y' : ℝ) (hYN : 2 * Y' ^ 2 < N) (S : Finset ℕ)
    (hSsub : S ⊆ kProducts D k) (hSwin : ∀ s ∈ S, Y' / 2 ≤ (s : ℝ) ∧ (s : ℝ) ≤ Y')
    (hScard : (Nat.choose D.card k : ℝ) / (2 * k) ≤ S.card) :
    (S.card : ℝ) ^ 2 / (cK k * 2 ^ omega N * (Nat.totient N : ℝ)) ≤
      prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by
  have hDp : ∀ p ∈ D, p.Prime := fun p hp => (hD p hp).1
  have hcop : ∀ s ∈ S, s.Coprime N := by
    intro s hs
    obtain ⟨T, hT, rfl⟩ := mem_image.1 (hSsub hs)
    exact Nat.Coprime.prod_left (fun p hp => (hD p ((mem_powersetCard.1 hT).1 hp)).2.2.2)
  have hYN' : Y' ^ 2 < N := by nlinarith [sq_nonneg Y']
  have hprod : ∀ a ∈ S, ∀ b ∈ S, a * b < N := by
    intro a ha b hb
    have h1 := (hSwin a ha).2
    have h2 := (hSwin b hb).2
    have h0a : (0 : ℝ) ≤ a := Nat.cast_nonneg _
    have h0b : (0 : ℝ) ≤ b := Nat.cast_nonneg _
    have : ((a * b : ℕ) : ℝ) < N := by
      push_cast
      calc (a : ℝ) * b ≤ Y' * Y' := mul_le_mul h1 h2 h0b (h0a.trans h1)
        _ = Y' ^ 2 := by ring
        _ < N := hYN'
    exact_mod_cast this
  have hmom := moment_bound N hodd S hcop hprod
  have hquad : (quadCount S : ℝ) ≤ 2 * ((2 * k).choose k : ℝ) * (D.card.choose k : ℝ) ^ 2 := by
    exact_mod_cast quadCount_kprod_le D k hk hDp hDcard S hSsub
  have hk' : (0 : ℝ) < k := by exact_mod_cast hk
  have hC0 : (0 : ℝ) ≤ (D.card.choose k : ℝ) / (2 * k) := by positivity
  have hck : 2 * ((2 * k).choose k : ℝ) * (D.card.choose k : ℝ) ^ 2 ≤ cK k * (S.card : ℝ) ^ 2 := by
    have h1 : ((D.card.choose k : ℝ) / (2 * k)) ^ 2 ≤ (S.card : ℝ) ^ 2 :=
      pow_le_pow_left₀ hC0 hScard 2
    have h2 : 2 * ((2 * k).choose k : ℝ) * (D.card.choose k : ℝ) ^ 2 =
        8 * (k : ℝ) ^ 2 * ((2 * k).choose k : ℝ) * ((D.card.choose k : ℝ) / (2 * k)) ^ 2 := by
      field_simp; ring
    rw [h2, cK]
    have : (0 : ℝ) ≤ 8 * (k : ℝ) ^ 2 * ((2 * k).choose k : ℝ) := by positivity
    exact mul_le_mul_of_nonneg_left h1 this
  have hprob0 : 0 ≤ prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) :=
    div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _)
  have hphi : (0 : ℝ) < Nat.totient N := by
    exact_mod_cast Nat.totient_pos.2 (NeZero.pos N)
  have hcK : 0 < cK k := by
    rw [cK]
    have : (0 : ℝ) < ((2 * k).choose k : ℝ) := by exact_mod_cast Nat.choose_pos (by omega)
    positivity
  rcases Nat.eq_zero_or_pos S.card with h0 | hpos
  · rw [h0]; simp only [Nat.cast_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true,
      zero_pow, zero_div]; exact hprob0
  have hS2 : (0 : ℝ) < (S.card : ℝ) ^ 2 := by positivity
  have key : (S.card : ℝ) ^ 2 * (S.card : ℝ) ^ 2 ≤
      (2 ^ omega N * (Nat.totient N : ℝ) * cK k *
        prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N))) * (S.card : ℝ) ^ 2 := by
    calc (S.card : ℝ) ^ 2 * (S.card : ℝ) ^ 2 = (S.card : ℝ) ^ 4 := by ring
      _ ≤ 2 ^ omega N * (Nat.totient N : ℝ) * (quadCount S : ℝ) *
          prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := hmom
      _ ≤ 2 ^ omega N * (Nat.totient N : ℝ) * (cK k * (S.card : ℝ) ^ 2) *
          prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by
          gcongr; exact hquad.trans hck
      _ = _ := by ring
  have := le_of_mul_le_mul_right key hS2
  rw [div_le_iff₀ (by positivity)]
  linarith

/-! ## Theorem K: helpers -/

lemma pow_div_le_choose (n k : ℕ) (hkn : k ≤ n) : ((n : ℝ) / k) ^ k ≤ n.choose k := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp
  have hterm : ∀ i ∈ range k, n * (k - i) ≤ k * (n - i) := by
    intro i hi
    have hik : i < k := mem_range.1 hi
    obtain ⟨a, rfl⟩ := Nat.exists_eq_add_of_le hik.le
    obtain ⟨b, rfl⟩ := Nat.exists_eq_add_of_le hkn
    rw [show i + a + b - i = a + b by omega, Nat.add_sub_cancel_left]
    nlinarith
  have h1 := Finset.prod_le_prod hterm
  rw [prod_mul_distrib, prod_mul_distrib, prod_const, prod_const, card_range,
    ← Nat.descFactorial_eq_prod_range, ← Nat.descFactorial_eq_prod_range, Nat.descFactorial_self,
    Nat.descFactorial_eq_factorial_mul_choose] at h1
  have h2 : n ^ k ≤ k ^ k * n.choose k := by
    have : n ^ k * k.factorial ≤ (k ^ k * n.choose k) * k.factorial := by
      calc n ^ k * k.factorial ≤ k ^ k * (k.factorial * n.choose k) := h1
        _ = (k ^ k * n.choose k) * k.factorial := by ring
    exact Nat.le_of_mul_le_mul_right this (Nat.factorial_pos k)
  have hkpos : (0 : ℝ) < (k : ℝ) ^ k := by positivity
  rw [div_pow, div_le_iff₀ hkpos]
  have : ((n ^ k : ℕ) : ℝ) ≤ ((k ^ k * n.choose k : ℕ) : ℝ) := by exact_mod_cast h2
  push_cast at this
  linarith

lemma card_kProducts (D : Finset ℕ) (k : ℕ) (hD : ∀ p ∈ D, p.Prime) :
    (kProducts D k).card = D.card.choose k := by
  rw [kProducts, card_image_of_injOn, card_powersetCard]
  intro A hA B hB h
  rw [mem_coe, mem_powersetCard] at hA hB
  have hA' : ∀ p ∈ A, p.Prime := fun p hp => hD p (hA.1 hp)
  have hB' : ∀ p ∈ B, p.Prime := fun p hp => hD p (hB.1 hp)
  rw [← Nat.primeFactors_prod hA', ← Nat.primeFactors_prod hB']
  simp only at h
  rw [h]

lemma kProd_bounds (D : Finset ℕ) (k Y : ℕ) (hD : ∀ p ∈ D, Y < p ∧ p ≤ 2 * Y) {s : ℕ}
    (hs : s ∈ kProducts D k) : (Y + 1) ^ k ≤ s ∧ s ≤ (2 * Y) ^ k := by
  obtain ⟨T, hT, rfl⟩ := mem_image.1 hs
  rw [mem_powersetCard] at hT
  constructor
  · rw [← hT.2]; exact Finset.pow_card_le_prod _ _ _ (fun p hp => (hD p (hT.1 hp)).1)
  · rw [← hT.2]; exact Finset.prod_le_pow_card _ _ _ (fun p hp => (hD p (hT.1 hp)).2)

lemma isProdAtMostK_of_dvd (D : Finset ℕ) (k : ℕ) (hD : ∀ p ∈ D, p.Prime) {a d : ℕ}
    (ha : a ∈ kProducts D k) (hd : d ∣ a) : IsProdAtMostK D k d := by
  obtain ⟨hP, hprod, hpos⟩ := kProd_spec D k hD ha
  rw [mem_powersetCard] at hP
  have hsq : Squarefree a := by
    rw [← hprod]
    refine Finset.squarefree_prod_of_pairwise_isCoprime (fun p hp q hq hpq => ?_)
      (fun p hp => (Nat.prime_of_mem_primeFactors hp).squarefree)
    simp only [Function.onFun, ← Nat.coprime_iff_isRelPrime]
    exact (Nat.coprime_primes (Nat.prime_of_mem_primeFactors hp)
      (Nat.prime_of_mem_primeFactors hq)).2 hpq
  have hsub : d.primeFactors ⊆ a.primeFactors := Nat.primeFactors_mono hd hpos.ne'
  refine ⟨d.primeFactors, hsub.trans hP.1, ?_, ?_⟩
  · rw [← hP.2]; exact card_le_card hsub
  · exact (Nat.prod_primeFactors_of_squarefree (hsq.squarefree_of_dvd hd)).symm

lemma card_sq_fiber (N : ℕ) [NeZero N] (F : Finset (ZMod N)ˣ) (hF : F ⊆ unitSquares N) :
    (muTwo N).card * F.card ≤ (univ.filter (fun u : (ZMod N)ˣ => u ^ 2 ∈ F)).card := by
  classical
  have hfib := card_eq_sum_card_fiberwise (s := univ.filter (fun u : (ZMod N)ˣ => u ^ 2 ∈ F))
    (f := fun u => u ^ 2) (t := F) (fun u hu => (mem_filter.1 hu).2)
  rw [hfib, mul_comm, ← smul_eq_mul, ← sum_const]
  refine sum_le_sum (fun R hR => ?_)
  obtain ⟨w, -, hw⟩ := mem_image.1 (hF hR)
  refine card_le_card_of_injOn (fun ε => ε * w) ?_ ?_
  · intro ε hε
    rw [mem_coe, muTwo, mem_filter] at hε
    rw [mem_coe, mem_filter, mem_filter]
    have : (ε * w) ^ 2 = R := by rw [mul_pow, hε.2, one_mul, hw]
    exact ⟨⟨mem_univ _, this ▸ hR⟩, this⟩
  · intro a _ b _ h
    exact mul_right_cancel h

lemma window_exists (s X k : ℕ) (hk : 1 ≤ k) (hs1 : 1 ≤ s) (hsX : s ≤ X) (hX : X ≤ s * 4 ^ k) :
    ∃ i < 2 * k, s * 2 ^ i ≤ X ∧ X ≤ s * 2 ^ (i + 1) := by
  classical
  have hex : ∃ i, X ≤ s * 2 ^ (i + 1) := ⟨X, by
    have := Nat.lt_two_pow_self (n := X + 1)
    nlinarith⟩
  refine ⟨Nat.find hex, ?_, ?_, Nat.find_spec hex⟩
  · have : Nat.find hex ≤ 2 * k - 1 := Nat.find_min' hex (by
      rw [show 2 * k - 1 + 1 = 2 * k by omega, pow_mul]; norm_num; exact hX)
    omega
  · rcases Nat.eq_zero_or_pos (Nat.find hex) with h0 | hpos
    · rw [h0]; simpa using hsX
    · have := Nat.find_min hex (show Nat.find hex - 1 < Nat.find hex by omega)
      rw [show Nat.find hex - 1 + 1 = Nat.find hex by omega] at this
      omega


/-- MAIN. [pp] Theorem thm:K, combinatorial core. Let `N` be odd, `X` with `2X² < N`, and `Y` with
`(2Y)^k ≤ X < (4Y)^k` (as for `J = (Y, 2Y]` in [pp] §sec:params). Let `D ⊆ J` be primes coprime to
`N` with `|D| ≥ 4k²`. For `r` uniform on `[1, N−1]`, with probability at least
`(|D|/k)^{2k} / (4k² c_k 2^{ω(N)} N)` the root `r` is a unit and every shortest nonzero vector
`(u₀, v₀)` of `Λ_{r²}` has `u₀v₀ ≠ 0` and both `|u₀|`, `|v₀|` products of at most `k` distinct
primes of `D` (so the shape test of DESCEND_ONCE passes with targets in `D`). -/
theorem complete_descent_prob (N k X Y : ℕ) [NeZero N] (hodd : Odd N) (hk : 1 ≤ k)
    (hXN : 2 * X ^ 2 < N) (hYX : (2 * Y) ^ k ≤ X) (hXY : X < (4 * Y) ^ k)
    (D : Finset ℕ) (hD : ∀ p ∈ D, p.Prime ∧ Y < p ∧ p ≤ 2 * Y ∧ Nat.Coprime p N)
    (hDcard : 4 * k ^ 2 ≤ D.card) :
    ((D.card : ℝ) / k) ^ (2 * k) / (4 * (k : ℝ) ^ 2 * cK k * 2 ^ omega N * N) ≤
      prob (Icc 1 (N - 1)) (fun r : ℕ => IsUnit (r : ZMod N) ∧
        ∀ x y : ℤ, IsShortestVec N ((r : ZMod N) ^ 2) x y →
          x * y ≠ 0 ∧ IsProdAtMostK D k x.natAbs ∧ IsProdAtMostK D k y.natAbs) := by
  classical
  have hDp : ∀ p ∈ D, p.Prime := fun p hp => (hD p hp).1
  set n := D.card
  set K := kProducts D k
  have hkn : k ≤ n := by nlinarith
  have hY1 : 1 ≤ Y := by
    obtain ⟨p, hp⟩ : D.Nonempty := card_pos.1 (by nlinarith)
    have := hD p hp
    have := this.1.two_le
    omega
  -- every `k`-product lies in one of the `2k` dyadic windows below `X`
  have hwin : ∀ s ∈ K, ∃ i < 2 * k, s * 2 ^ i ≤ X ∧ X ≤ s * 2 ^ (i + 1) := by
    intro s hs
    obtain ⟨h1, h2⟩ := kProd_bounds D k Y (fun p hp => ⟨(hD p hp).2.1, (hD p hp).2.2.1⟩) hs
    apply window_exists s X k hk
    · exact le_trans (Nat.one_le_pow _ _ (by omega)) h1
    · exact h2.trans hYX
    · have : Y ^ k ≤ s := le_trans (Nat.pow_le_pow_left (by omega) k) h1
      calc X ≤ (4 * Y) ^ k := hXY.le
        _ = 4 ^ k * Y ^ k := mul_pow _ _ _
        _ ≤ 4 ^ k * s := Nat.mul_le_mul_left _ this
        _ = s * 4 ^ k := mul_comm _ _
  set Sf : ℕ → Finset ℕ := fun i => K.filter (fun s => s * 2 ^ i ≤ X ∧ X ≤ s * 2 ^ (i + 1))
  have hcover : K.card ≤ ∑ i ∈ range (2 * k), (Sf i).card := by
    refine le_trans (card_le_card ?_) card_biUnion_le
    intro s hs
    obtain ⟨i, hi, hw⟩ := hwin s hs
    exact mem_biUnion.2 ⟨i, mem_range.2 hi, mem_filter.2 ⟨hs, hw⟩⟩
  have hKcard : K.card = n.choose k := card_kProducts D k hDp
  have hk' : (0 : ℝ) < k := by exact_mod_cast hk
  obtain ⟨i, -, hSi⟩ : ∃ i ∈ range (2 * k), (n.choose k : ℝ) / (2 * k) ≤ (Sf i).card := by
    refine exists_le_of_sum_le (by rw [nonempty_range_iff]; omega) ?_
    rw [sum_const, card_range, nsmul_eq_mul]
    calc ((2 * k : ℕ) : ℝ) * ((n.choose k : ℝ) / (2 * k)) = n.choose k := by
          push_cast; field_simp
      _ = (K.card : ℝ) := by rw [hKcard]
      _ ≤ ∑ i ∈ range (2 * k), ((Sf i).card : ℝ) := by exact_mod_cast hcover
  set S := Sf i
  set Y' : ℝ := (X : ℝ) / 2 ^ i
  have h2i : (1 : ℝ) ≤ 2 ^ i := one_le_pow₀ (by norm_num)
  have h2ipos : (0 : ℝ) < 2 ^ i := by positivity
  have hX0 : (0 : ℝ) ≤ X := Nat.cast_nonneg _
  have hY'X : Y' ≤ X := div_le_self hX0 h2i
  have hY'0 : 0 ≤ Y' := by positivity
  have hYN : 2 * Y' ^ 2 < N := by
    have : Y' ^ 2 ≤ (X : ℝ) ^ 2 := pow_le_pow_left₀ hY'0 hY'X 2
    have : (2 * X ^ 2 : ℝ) < N := by exact_mod_cast hXN
    linarith
  have hSK : S ⊆ kProducts D k := filter_subset _ _
  have hSwin : ∀ s ∈ S, Y' / 2 ≤ (s : ℝ) ∧ (s : ℝ) ≤ Y' := by
    intro s hs
    obtain ⟨-, h1, h2⟩ := mem_filter.1 hs
    have h1' : (s : ℝ) * 2 ^ i ≤ X := by exact_mod_cast h1
    have h2' : (X : ℝ) ≤ s * 2 ^ (i + 1) := by exact_mod_cast h2
    constructor
    · rw [div_div, ← pow_succ, div_le_iff₀ (by positivity)]; exact h2'
    · rw [le_div_iff₀ h2ipos]; exact h1'
  have hmain := split_counting_kprod N hodd k hk Y D hD hDcard Y' hYN S hSK hSwin hSi
  -- from `R` uniform on `Q_N` to `r` uniform on `[1, N-1]`
  have hX2 : 2 ≤ X := by
    have : 2 ≤ (2 * Y) ^ k := le_trans (by omega) (Nat.le_self_pow (by omega) _)
    omega
  have hN1 : 1 < N := by nlinarith
  have : Fact (1 < N) := ⟨hN1⟩
  set Q := unitSquares N
  set F := Q.filter (fun R : (ZMod N)ˣ => SplitEvent N S (R : ZMod N))
  set U := univ.filter (fun u : (ZMod N)ˣ => u ^ 2 ∈ F)
  have hU : (muTwo N).card * F.card ≤ U.card := card_sq_fiber N F (filter_subset _ _)
  set E := fun r : ℕ => IsUnit (r : ZMod N) ∧ ∀ x y : ℤ, IsShortestVec N ((r : ZMod N) ^ 2) x y →
          x * y ≠ 0 ∧ IsProdAtMostK D k x.natAbs ∧ IsProdAtMostK D k y.natAbs
  have hUE : U.card ≤ ((Icc 1 (N - 1)).filter E).card := by
    refine card_le_card_of_injOn (fun u => (u : ZMod N).val) ?_ ?_
    · intro u hu
      rw [mem_coe] at hu ⊢
      have hu2 : u ^ 2 ∈ F := (mem_filter.1 hu).2
      have hval : (((u : ZMod N).val : ℕ) : ZMod N) = (u : ZMod N) := ZMod.natCast_zmod_val _
      rw [mem_filter, mem_Icc]
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · rw [Nat.one_le_iff_ne_zero]
        intro h
        exact u.ne_zero ((ZMod.val_eq_zero _).1 h)
      · show (u : ZMod N).val ≤ N - 1
        have := ZMod.val_lt (u : ZMod N); omega
      · show E (u : ZMod N).val
        simp only [E]
        rw [hval]
        refine ⟨u.isUnit, ?_⟩
        intro x y hxy
        obtain ⟨-, a, ha, b, hb, hab⟩ := mem_filter.1 hu2
        have hR : (u : ZMod N) ^ 2 = ((u ^ 2 : (ZMod N)ˣ) : ZMod N) := (Units.val_pow_eq_pow_val u 2).symm
        rw [hR] at hxy
        have haK : a ∈ kProducts D k := hSK ha
        have hbK : b ∈ kProducts D k := hSK hb
        have ha1 : 1 ≤ a := (kProd_spec D k hDp haK).2.2
        have hb1 : 1 ≤ b := (kProd_spec D k hDp hbK).2.2
        obtain ⟨j, hj1, hj2⟩ := split_multiple_of_shortest N Y' hYN _ a b ha1 hb1
          (hSwin a ha).2 (hSwin b hb).2 hab x y hxy
        have hxa : a = j.natAbs * x.natAbs := by
          have := congrArg Int.natAbs hj1
          rw [Int.natAbs_mul] at this
          simpa using this
        have hyb : b = j.natAbs * y.natAbs := by
          have := congrArg Int.natAbs hj2
          rw [Int.natAbs_mul] at this
          simpa using this
        refine ⟨?_, isProdAtMostK_of_dvd D k hDp haK (Dvd.intro_left _ hxa.symm),
          isProdAtMostK_of_dvd D k hDp hbK (Dvd.intro_left _ hyb.symm)⟩
        intro h
        rcases mul_eq_zero.1 h with h | h
        · rw [h] at hxa; simp at hxa; omega
        · rw [h] at hyb; simp at hyb; omega
    · intro u _ w _ h
      exact Units.ext (ZMod.val_injective N h)
  -- the final chain of inequalities
  have hphi : (0 : ℝ) < Nat.totient N := by exact_mod_cast Nat.totient_pos.2 (NeZero.pos N)
  have hQpos : (0 : ℝ) < Q.card := by exact_mod_cast unitSquares_card_pos N
  have hmuQ : (Q.card : ℝ) * (muTwo N).card = Nat.totient N := by
    have := card_unitSquares_mul N
    rw [ZMod.card_units_eq_totient] at this
    exact_mod_cast this
  have hcK : 0 < cK k := by
    rw [cK]
    have : (0 : ℝ) < ((2 * k).choose k : ℝ) := by exact_mod_cast Nat.choose_pos (by omega)
    positivity
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (NeZero.pos N)
  have hN1r : (0 : ℝ) < ((N - 1 : ℕ) : ℝ) := by exact_mod_cast (by omega : 0 < N - 1)
  have hN1le : ((N - 1 : ℕ) : ℝ) ≤ N := by exact_mod_cast Nat.sub_le N 1
  have hprob : prob (Icc 1 (N - 1)) E = (((Icc 1 (N - 1)).filter E).card : ℝ) / ((N - 1 : ℕ) : ℝ) := by
    rw [prob, Nat.card_Icc, show N - 1 + 1 - 1 = N - 1 by omega]
    convert rfl
  have hprobQ : prob Q (fun R => SplitEvent N S (R : ZMod N)) = (F.card : ℝ) / Q.card := rfl
  rw [hprobQ] at hmain
  have hS0 : ((n : ℝ) / k) ^ k / (2 * k) ≤ S.card := by
    refine le_trans ?_ hSi
    exact div_le_div_of_nonneg_right (pow_div_le_choose n k hkn) (by positivity)
  have hD' : (0 : ℝ) < cK k * 2 ^ omega N * N := by positivity
  calc ((n : ℝ) / k) ^ (2 * k) / (4 * (k : ℝ) ^ 2 * cK k * 2 ^ omega N * N)
      = (((n : ℝ) / k) ^ k / (2 * k)) ^ 2 / (cK k * 2 ^ omega N * N) := by
        rw [pow_mul']; field_simp; ring
    _ ≤ (S.card : ℝ) ^ 2 / (cK k * 2 ^ omega N * N) := by
        gcongr
    _ = ((S.card : ℝ) ^ 2 / (cK k * 2 ^ omega N * (Nat.totient N : ℝ))) * Nat.totient N / N := by
        field_simp
    _ ≤ ((F.card : ℝ) / Q.card) * Nat.totient N / N := by gcongr
    _ = ((muTwo N).card * F.card : ℝ) / N := by
        rw [← hmuQ]; field_simp
    _ ≤ (((Icc 1 (N - 1)).filter E).card : ℝ) / N := by
        gcongr
        exact_mod_cast hU.trans hUE
    _ ≤ (((Icc 1 (N - 1)).filter E).card : ℝ) / ((N - 1 : ℕ) : ℝ) :=
        div_le_div_of_nonneg_left (Nat.cast_nonneg _) hN1r hN1le
    _ = prob (Icc 1 (N - 1)) E := hprob.symm

end FactoringA
