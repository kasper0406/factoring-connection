# snfsrig: special forms (SNFS) and rigorous sub-L[2/5] factoring — results

Agent `snfsrig`, phase 14, 2026-10-20.

**Files.**
- This note.
  - No prediction-testing experiment was run, so there is no pre-registration (`notes/snfsrig_prereg.md` does not
    exist).
  - The code only checks proved lemmas and does exact exponent arithmetic.
- Code (tests in `tests/test_snfsrig.py`):
  - `src/snfsrig/sizemodel.py`: exact size model (§1).
  - `src/snfsrig/rigidity.py`: checkers for Lemmas RES/HR/RS/PH, Prop DI and the value-set count of Prop RB (§2).
    It is diagnostic.
  - `src/snfsrig/cyclo.py`: cyclotomic and Aurifeuillian pieces of r^e ± 1 (§3). It is operational-safe and never
    factors N.
- Codex:
  - `codex/snfsrig_1.{md,out}`: blind round, launched before my derivations were written.
  - `codex/snfsrig_2.{md,out}`: adversarial referee of the first draft (§6).

**Base texts.**
- "[U]": Theorem U and its design (`notes/lv5unif_results.md`). Everything the project calls "Theorem U" is internal
  and unrefereed.
- "[rec]": the size model (`notes/recurse_results.md` §3.1–3.4).
- "[nest]": Props SI/R/R2 and Theorem NR (`notes/nested_results.md`).
- "[ld]": the w multiplier (`notes/ldbound_results.md`).
- The field law: `notes/thmA_proof.md` §1.
- The HECM barrier: `notes/lemmaM_results.md` (c).
- "[g2rh]": RH-conditional L_p[1/2, 2] (`notes/genus2rh_results.md`).

**Labels.**
- **proved**: a mathematical proof is given here, from published results.
- **proved (model arithmetic)**: an exact consequence of the stated size-model premises, which are themselves
  heuristic.
- **heuristic**.
- **open**.
- **assessment**: a judgement, not a theorem.
- **literature**, with the access level.

**Scope.** "Special family" means **N = r^e − s with r ≥ 2 and s ≠ 0 fixed, e → ∞**, without any promise about the
prime factors of N. Families defined by their factors (close primes, smooth p − 1, N = p^r q) are rigorously easy and
are a different question (§3.2).

---

## 0. Verdict

**Plain language.**
- Special forms do not give a rigorous algorithm below L[2/5] with known tools. We can say precisely where it breaks.
- The tiny SNFS polynomial does help the *size model*. It removes the 2/5 barrier for designs without algebraic
  children, so in the model rational-only multi-level descent reaches L[1/3].
- Rigorously it does not help, for two reasons.
  - **(O1)** Every relation-supply proof in the project averages over random fields.
    - Re-randomising the special field either reuses the values of the fixed pair (scaled, on a transformed region)
      or a small fixed family of pairs, which we prove for every device we examined, or costs polynomial height about
      N^{1/d}, which turns SNFS back into GNFS.
    - The polynomials sharing the special root obey this dichotomy as a theorem (Lemma HR).
    - What remains is a smoothness estimate for the values of one fixed binary form at L-scale. Nothing of the kind is
      known.
    - We also show that proofs using only *linear* statistics of smooth numbers (the class 𝒦 of Prop RB) cannot supply
      it for any small family of fixed fields.
  - **(O2)** The Dixon-split architecture must factor the split halves u₀, v₀ completely.
    - These are arbitrary integers of size √N, not special ones, so the recursion that removes recognition costs in
      general cannot stay inside the family.
    - Unconditionally, none of the rigorous recognisers available (Theorem U recursion, HECM, the Pomerance-type
      averaged ECM) has a proved bound below what α = 2/5 already allows.
    - Under RH, the project's genus-2 method [g2rh] removes (O2). Then (O1) alone remains.
- **Best rigorous α for the family: 2/5** [assessment], from Theorem U, which applies to every N and is conditional on
  the internal chain.
- The special structure improves only the constant:
  - **L[2/5, c_U·ε_e^{2/5} + o(1)]** for r^e ± 1 via cyclotomic pieces, with ε_e = φ(e)/e (resp. φ(2e)/e);
  - **L[2/5, o(1)] along N = r^e − 1 with e running over the primorials**;
  - from published results only (Lenstra–Pomerance 1992 on the pieces), **L[1/2, o(1)]** along the same sequence.

**Results.**

1. **Size model (§1).**
   - **Lemma RES [proved].** If coprime f, g ∈ Z[x] have a common root mod N, then
     N ≤ |Res(f, g)| ≤ ‖f‖₂^{deg g}‖g‖₂^{deg f}. So every pair with a linear side has coefficient exponent
     θ_coeff ≥ 1/d, and SNFS attains this.
   - **Prop UF [proved (model arithmetic), within the Dixon-split one-level architecture of [rec] (R1)–(R3)].**
     - If θ_coeff ≥ c·d^{−p}, the one-level floor is (1 + p)/(3p + 2).
     - RES lets us apply this with p = 1, so the floor is 2/5 for every pair with a linear side, special or not.
     - The (R3) floor changes only in its constant: θ ≈ 1/√k for SNFS, against √2/√k for heuristic GNFS and for
       Theorem U's design at 2/5.
   - Hypothetical constants. Suppose SNFS had Theorem U's supply. The split halves would still be factored by
     Theorem U (Prop REC), giving c ≈ 3.34 / 3.85 / 4.79 / 6.19 against Theorem U's 3.51 / 4.05 / 5.04 / 6.50
     (C_γ = 0 / 1 / 5 / 20, grid upper bounds). That is about 5% lower, and **not a theorem**.
   - Multi-level descent in one field: 1/3 in the model, unchanged from GNFS.
   - **Prop R2-S [proved (model arithmetic)].** Nested's Prop R2 recomputed with its premise (M1) replaced by the SNFS
     supply of special polynomials.
     - Choosing a special degree d(z) ≍ Ξ/z at each node, rational-only descent contracts targets by t^{1−3α}.
     - So the no-algebraic-children floor drops from 2/5 to **1/3**, attained with O(log log N) levels.
     - This is exactly where small norms shift R2. Codex's blind round found the same.
2. **Randomness (§2).**
   - Proved:
     - Lemma HR (height rigidity at the special root);
     - Lemma PH (a perturbation f₀ + (x − m)c costs height ≥ (m − 1)‖c‖ − ‖f₀‖);
     - Lemma RS (a rational side sharing the root with height < N/(m + 1) is a multiple of x − m);
     - Prop DI (w·N, kN, root twists or scalings, shifts, Möbius maps, and the other j for the same d act by
       (F, G)(P) = (c₁F₀(TP), c₂G₀(TP)), T ∈ M₂(Z));
     - Prop RB (a relativisation barrier for any family of N^{o(ν)} fixed pairs, for a listed class 𝒦 of *linear*
       statistical inputs: deleting the at most N^{(2+o(1))ν} queried rational values leaves those inputs valid with
       negligible extra error, and makes supply 0).
   - Quadratic inputs such as L2_fam and BDH are not covered by RB.
   - **Device table (§2.6).** No device we examined that keeps the SNFS norms keeps any of FM, LD*, Kern, Lemma W or
     Theorem R. The devices that restore them are exactly Theorem U's design (heights ≥ N^{1/d}).
   - Not proved: that *every* small-height pair with a known root is special or a twist of a special one. That is the
     (M1) counting heuristic. HR proves it only at the special root.
3. **Recognition (Prop REC, §1.5) [proved accounting, unconditional].**
   - Each accepted row needs the complete factorisation of u₀, v₀ ≈ √N, whose prime factors are of size L[1 − α].
   - A recogniser with a proved L_p[β] bound fits the budget iff α ≥ β/(1 + β).
   - The bounds we have:
     - HECM: β = 2/3, giving 2/5. Its algorithmic part III is unpublished.
     - The explicit Pomerance-type ECM (Lenstra–Pomerance 1992 §7): exp(O((log y)^{6/7})), giving 6/13.
     - Theorem U recursion: L_N[2/5].
   - Under RH, [g2rh] gives β = 1/2, so α ≥ 1/3.
4. **Special families (§3).**
   - **Prop CY [proved, conditional on Theorem U].** L_N[2/5, c_Uε_e^{2/5} + o(1)].
   - **Corollary CY-P.** L[2/5, o(1)] along primorial e. The o(1) is ≍ (log log log N)^{−2/5}, which decays extremely
     slowly.
   - **Corollary CY-LP [proved, published inputs only].** L[1/2, ε_e^{1/2} + o(1)].
   - Aurifeuillian splitting gives a further factor 2^{−2/5} (resp. 2^{−1/2}).
5. **What the model identifies as missing (§4).**
   - (O1) Hypothesis RT: pointwise rational-target supply for the explicit Kummer forms x^d − s·r^j at L-scale, with
     the large cofactor primes in a recognisable set.
   - (O2) A recogniser for that set.
   - Characters: ERH for Lemma X_ERH, or LV §6 unconditionally at cost L[1/3, Θ(δ)].
   - Multi-field LV6 bookkeeping.
   - This is a list of missing inputs, not a proved implication.
   - The known friable-value results (BBDT, Lachand, Barbulescu–Lachand, Mounier) and prescribed-digit results
     (Bourgain, Swaenepoel) do not establish RT, even for its rational half.
6. **Literature (§5).**
   - No rigorous or GRH-conditional sub-1/2 result was found for factoring special-form N or for discrete logs in prime
     fields with special p.
   - SNFS (LLMP 1990: "this paper does not contain a rigorous mathematical result"), Gordon, Semaev, Joux–Lercier and
     Joux–Pierrot are heuristic. Lee–Venkatesan is for random polynomials.
   - The only rigorous L[α < 1/2] for a natural family that we found is Kleinjung–Wesolowski's, for F_{p^n} with
     p = L(α). It concerns extension fields, not integers or prime fields.

**Weakest points.**
- The barriers are method barriers, not impossibility theorems:
  - RB covers proofs from a listed class of statistical inputs;
  - REC covers the Dixon-split architecture with the proved recognisers we know;
  - DI covers the listed devices.
- Everything positive at α = 2/5 rests on Theorem U's internal chain.

---

## 1. Size model

**Notation.**
- t = log N, D = t/log t, Λ_α = t^α(log t)^{1−α}, B = e^{aΛ_α}.
- U ≍ D^α is the largest u with u log u ≤ Λ_α.
- Ξ_α := Λ_αU ≍ t^{2α}(log t)^{1−2α}.
- In the random-integer model, a value of log-size X is B-smooth at cost L[α]^{O(1)} iff X = O(Ξ_α) (CEP).

**The special polynomials.** N = r^e − s with r, s fixed and gcd(r, N) = 1. For each d ≥ 2 put q_d := ⌈e/d⌉,
m_d := r^{q_d} and j_d := dq_d − e ∈ [0, d). Let

  f_d(x) := x^d − s·r^{j_d}.

Then
- f_d(m_d) = r^{j_d}N ≡ 0 (mod N);
- log‖f_d‖ ≤ d log r + log|s|;
- log m_d ≤ t/d + log r + O(1/N).

Codex snfsrig_1 corrected the sign of j in my prompt.

**Eisenstein family.** 𝔈 := {r^e − s : e ≥ 1} for fixed r ≥ 2 and a fixed prime s ∤ r.
- Every f_d is Eisenstein at s, hence irreducible (codex's suggestion).
- For Cunningham numbers r^e ± 1 one uses x^d ∓ r^{j} with gcd(j, d) = 1 and r not a perfect power. These are
  irreducible by Capelli's theorem (x^d − a is reducible only if a is a p-th power for some p | d, or 4 | d and
  a ∈ −4Q⁴; neither can happen for a = ±r^j with gcd(j, d) = 1).

### 1.1 Lemma RES and the coefficient floor [proved]

**Lemma RES.** Let f, g ∈ Z[x] with deg f = d ≥ 1 and deg g = δ ≥ 1. Suppose they are coprime in Q[x] and
f(m) ≡ g(m) ≡ 0 (mod N) for some integer m. Then N | Res(f, g), Res(f, g) ≠ 0, and
N ≤ |Res(f, g)| ≤ ‖f‖₂^δ‖g‖₂^d.

*Proof.*
- The adjugate of the Sylvester matrix gives A, B ∈ Z[x] with Af + Bg = Res(f, g) (Lang, *Algebra*, IV §8).
  Evaluating at m gives Res ≡ 0 (mod N).
- Res ≠ 0 because f and g have no common root in Q̄.
- The Sylvester matrix consists of δ shifted copies of f's coefficient vector and d copies of g's. Hadamard's
  inequality on its rows gives |Res| ≤ ‖f‖₂^δ‖g‖₂^d. ∎

**Corollary (linear side).** Let g = g₁x − g₀ and let f be irreducible of degree d ≥ 2. Then ‖f‖₂‖g‖₂^d ≥ N. Write
a := log‖f‖₂ ≥ 0 and b := log‖g‖₂, so a + db ≥ t. Then a + b ≥ t/d + a(1 − 1/d) ≥ t/d. So
**θ_coeff := (log‖f‖₂ + log‖g‖₂)/t ≥ 1/d** for every NFS pair with a linear side.
- SNFS attains equality up to O(d log r/t).
- RES also covers non-monic f, projective roots (via homogeneous resultants) and several rational sides pairwise
  (codex snfsrig_2).
- Skewed boxes change the weighted heights but not the resultant obstruction.

`rigidity.check_res_floor` checks the divisibility and the Hadamard bound exactly on SNFS pairs.

### 1.2 One level [proved (model arithmetic)]

**Architecture.** The project's Dixon split + one-level special-q descent ([rec] §3.3). The premises are:
- (R1) Split. The halves u₀, v₀ are products of k primes of size N^{1/(2k)}, with k = O(D^α).
- (R2) Smoothness. θ = O(D^{2α−1}), in the random-integer model.
- (R3) Norm product. With box radius V ≈ √ℓ·S, the rational value is ≈ ‖g‖V/ℓ and the algebraic value is ≈ ‖f‖V^d.
  So θ = θ_coeff + (d − 1)ν with ν = 1/(4k), up to O(d log S/t).

For SNFS, θ_SNFS = 1/d + (d − 1)/(4k). Theorem U's randomisation floors (h ≥ 2ν, level, window) are absent because
nothing is randomised.

**Prop UF.** Suppose θ_coeff ≥ c·d^{−p} with c > 0 and p ≥ 0. Then α ≥ (1 + p)/(3p + 2), and α ≥ 1/2 when p = 0.
By the Corollary of RES the hypothesis holds with p = 1 for every pair with a linear side. So **α ≥ 2/5 for every
one-level design of this architecture, special or not.**

*Proof.*
- Put d ≍ D^x and k ≍ D^y.
- Applying (R2) to the two terms of θ gives px ≥ 1 − 2α and y − x ≥ 1 − 2α. (R1) gives y ≤ α.
- Eliminating x and y: (1 − 2α)(1 + p)/p ≤ α, i.e. α ≥ (1 + p)/(3p + 2), which is decreasing in p.
- See `sizemodel.one_level_floor`; the tests include a brute-force rational LP. ∎

The scope is the stated architecture. "One level" means Dixon-split primes descended in one special-q step.
Ordinary SNFS relation collection, with no split, is a different architecture (codex snfsrig_2).

**Constants at α = 2/5.** Put d = d₀D^{1/5} and k = κD^{2/5}, and θ₀ := lim θD^{1/5} (`sizemodel.theta0`).

| design | θ₀ | min over d₀ |
|---|---|---|
| SNFS | 1/d₀ + d₀/(4κ) | 1/√κ |
| heuristic GNFS (base m, both heights N^{1/(d+1)}) | 2/d₀ + d₀/(4κ) | √(2/κ) |
| Theorem U | 3/d₀ if κ ≥ d₀²/4, else 2/d₀ + d₀/(4κ) | √(2/κ), attained in the Kern regime d₀ = √(8κ) |

**Hypothetical constants.** These assume that SNFS had Theorem U's supply γ ≥ (log N)^{−C_γ}, which is exactly what
§2 shows we cannot prove.
- By Prop REC the split halves are factored by Theorem U, so the recursion term becomes a + s + 2^{−2/5}c_U.
- The cost formula is c = max(2a, a + s + max(a, σ), a + s + 2^{−2/5}c_U), with s = 2κ(1 + C_γ) and σ = 2θ₀/(5a)
  (`sizemodel.optimise_nonrec`, exact rational grid).

| C_γ | Theorem U (reproduced) | SNFS with U's supply, U-factored halves | (draft: SNFS self-recursive — inconsistent with REC) |
|---|---|---|---|
| 0 | 3.510 | 3.335 | 3.056 |
| 1 | 4.047 | 3.847 | 3.510 |
| 5 | 5.038 | 4.788 | 4.378 |
| 20 | 6.498 | 6.189 | 5.658 |

The last column is my first draft. The Opus proof critic pointed out that it is inconsistent with Prop REC.

### 1.3 Multi-level descent in one field [model arithmetic, as [rec] §3.4]

Take d = D^δ. A level-j prime of log-size t^{β_j} has values of log-sizes t^{1−δ} (rational, ≈ log m) and t^{δ+β_j}
(algebraic). These exponents are the same for GNFS and SNFS heights; only the constants differ.
- So the level count ⌈(1 − 2α)/(3α − 1)⌉ (`sizemodel.multilevel_levels`) is unchanged.
- So is the floor 1/3, attained with O(log log N) levels at constant contraction.
- The obstruction there is nested's: algebraic children are field-bound ([nest] Lemmas FB and OD).

### 1.4 Prop R2-S: designs without algebraic children [proved (model arithmetic)]

Nested's Prop R2 premises:
- (M1) (d + 1)log H_f + 2 log H_g ≳ t;
- (M2) typical sizes;
- (M3) random-integer smoothness;
- (M4) split targets z₁ ≍ t^{1−α}(log t)^α.

Under them, designs whose algebraic values are all B-smooth have α ≥ 2/5. Nested itself notes that (M1) is only a
counting heuristic, and that small-height f break it.

**Premise (M1-S).** For every 2 ≤ d ≤ D there is an irreducible pair with log H_f ≤ d log r + log|s| = O(d) and
log H_g ≤ t/d + O(1). This holds for the family f_d, m_d with r and s fixed; for 𝔈 every d works. (M2)–(M4) are
kept, and each node may choose its own d.

**Prop R2-S.** In this model:
- (i) for α > 1/3 a rational-only descent tree exists with ⌈(1 − 2α)/(3α − 1)⌉ levels;
- (ii) at α = 1/3 such a tree exists with O(log log N) levels at constant contraction;
- (iii) for α < 1/3, the child-smoothness threshold needed to make a *typical* node affordable exceeds the target
  (y/z → ∞), so this model gives no finite tree.

So the floor is **1/3** instead of 2/5.

*Derivation (model arithmetic).*
- **Node set-up.** A node has a rational target q = e^z. It searches a reduced special-q box of radius e^{σ_bΛ}√q in
  the field of degree d. This holds for regular lattices; the worst case is treated below. (σ_b is the box exponent,
  renamed from s, which clashed with s in r^e − s.)
- **Algebraic side.** The algebraic value has log-size d(z/2 + σ_bΛ) + O(d). B-smoothness within budget needs this to
  be ≤ φ_FΞ, so take d(z) ≈ φ_FΞ/(z/2 + σ_bΛ), rounded to an integer.
- **Rational side.** The rational cofactor has log-size t/d − z/2 + σ_bΛ ≈ t(z/2 + σ_bΛ)/(φ_FΞ). Its children are
  bounded by y ≈ (cofactor)/(φ_RU), so y/z ≈ t(1/2 + σ_bΛ/z)/(φ_Fφ_RΞU).
- **Contraction.** ΞU ≍ t^{3α}(log t)^{1−3α}, so the ratio is ≍ t^{1−3α}(log t)^{3α−1}
  (`sizemodel.r2s_contraction_exponent`).
- **(iii).** d ≤ φ_FΞ/(z/2) forces a rational cofactor ≥ tz/(2φ_FΞ), hence y/z → ∞. This is the same argument as codex
  snfsrig_1. It is a statement about typical-size nodes, not about every candidate.
- **(ii), α = 1/3.** Put E = Λ_{1/3}, Q = t^{2/3}(log t)^{1/3}, D′ = Q/E, so t = D′Q.
  - Take B = e^{bE}, u_F = φ_FD′, u_R = φ_RD′ and top targets z₁ = b₁Q.
  - Then d(z) ≈ φ_FbQ/(z/2 + σ_bE) and λ = (1/2 + σ_b/b)/(φ_Fφ_Rb).
  - Example: b = 4, σ_b = 3/4, φ_F = φ_R = 1, b₁ = 1. Then λ = 11/64. The unrounded top degree tends to 8 from below,
    so the actual degree is 7 (codex snfsrig_2), and rounding costs a factor ≤ 8/7 in λ.
  - The good-point margin is 2σ_b − (φ_F + φ_R)/3 = 5/6 > 0.
- **Worst-case reduced-basis boxes.** Here ‖P‖ ≤ 3e^{σ_bE}q, as in [nest] §8. The algebraic term is d(z + σ_bE),
  which gives the *continuous* bound λ = (1 + σ_b/b)/(φ_Fφ_Rb) = 19/64 and an unrounded top degree tending to 4 from
  below (`r2s_alpha_third(regular=False)`).
  - Rounding down to d = 3 multiplies the rational cofactor, and hence λ, by at most d*/⌊d*⌋ ≤ 4/3. So λ ≤ 19/48 + ε < 1
    (codex snfsrig_3).
  - In general, rounding costs at most a factor d*/(d* − 1) at the top, where d* ≥ 2 is needed, and 1 + o(1) further
    down.
  - The good-point margin 2σ_b − (φ_F + φ_R)/3 = 5/6 is computed in this note's accounting for both lattice settings.
    The Opus critic's estimate of ≈ 1/2 under nested's AD formula was not derived here; we claim only that the margin
    is positive in this accounting.
- **Tree size.** Depth O(log D′) = O(log log N) and fan-out O(t/E) give e^{O((log log N)²)} = L^{o(1)} nodes
  ([nest] Lemma TS). The degree rises from O(1) at the top to ≤ 2φ_FD′ at the leaves, and the coefficient term
  d log r = O(D′) = o(E) is negligible. ∎

**Relation to GNFS (corrected after both reviewers).** Prop R2 excludes trees that satisfy (M1) at every node. HR
(§2.1) shows that polynomials *sharing the special root* either are multiples of f₀ or have height ≳ N^{1/d}. So in the
model, a sub-2/5 rational-only tree must contain nodes whose pairs violate (M1).
- For N ∈ 𝔈 the only such pairs we know are the special ones (and the generator changes of §2.4).
- That *every* (M1)-violating pair with a known root is of this kind is **not proved**. It is the (M1) counting
  heuristic. My draft's exhaustive classification is withdrawn.
- Hybrid trees are not excluded.
- Codex's blind round reached the same 1/3 ("the former forbidden window disappears").

### 1.5 Recognition: Prop REC [proved accounting; unconditional]

Each accepted row of the Dixon-split architecture requires the complete factorisation (the shape certificate) of
the Gauss-reduced halves u₀, v₀ of a random square.
- These are arbitrary integers ≤ √(2N) whose prime factors lie in the target set ⊂ (Y, 2Y], Y = N^{1/(2k)} =
  L_N[1 − α].
- Every intermediate node of a multi-level tree likewise has to factor rational cofactors with prime factors up to
  L[1 − α].
- These integers are not special, so the special-family algorithm cannot be called on them. That would be circular
  (codex snfsrig_1 §1, "algorithmic qualification").

**Prop REC.** A recogniser with a proved expected cost L_p[β, O(1)] per prime p found, for all p in a target set of
density (log N)^{−O(1)}, meets the budget L_N[α]^{O(1)} iff β(1 − α) ≤ α, i.e. **α ≥ β/(1 + β)**. The bounds
available:

| recogniser | proved bound | α threshold |
|---|---|---|
| Theorem U recursion on u₀, v₀ (internal) | L_{√(2N)}[2/5, c_U] = L_N[2/5, 2^{−2/5}c_U] | 2/5 |
| HECM (Lenstra–Pila–Pomerance) | L_p[2/3, c₀] (part III unpublished, §5) | 2/5 (lemmaM (c)) |
| averaged ECM (Pomerance 1987; explicit form Lenstra–Pomerance 1992 §7) | exp(O((log y)^{6/7})), recognisable primes only | 6/13 |
| trial division, Pollard–Strassen, batch smoothness | Y^{Ω(1)}, or total cost ≥ π(Y) | 1 − α ≤ α only for α ≥ 1/2 |
| **under RH:** genus-2 HECM [g2rh] (RH + Younis, unrefereed) | L_p[1/2, 2 + o(1)] for every p | **1/3** |

So, **unconditionally**, none of these guarantees establishes recognition within L[α] for α < 2/5. These are the
proved upper bounds of the known methods; they are not lower bounds for recognition (codex snfsrig_2).

**Target restriction does not help cheaply.** Lemma S″ accepts any target set G ⊂ J. But with a set of relative
density δ, the counting guarantee drops by ≈ δ^{2k} with k ≍ D^α, so the guarantee needs δ ≥ (log N)^{−O(1)}. For
example, the primes in J with p − 1 being L[α]-smooth number at most Ψ(2Y, B). That is a relative density exp(−(1 +
o(1))u log u) with u ≍ D^{1−2α}, so δ^{2k} = exp(−Θ(D^{1−α}log D)) is not L[α]^{O(1)}. This is an upper bound for the
guarantee, not a proof that other target sets fail.

*Proof of the table.* The first row uses L_{√x}[α, c] = L_x[α, 2^{−α}c + o(1)] ([rec] §3.1). The other rows insert
y = 2Y = L_N[1 − α] into the cited bounds. ∎

---

## 2. Recovering randomness: what survives

The relation-supply proofs use field randomness in four places ([U] §§4–6, [ld], thmA §1) [assessment from reading
those proofs]:
- **(R) Theorem R / FIB.** For fixed P, G_g(P) = g₁a − g₀b is near-uniform on an interval as g ranges over the 2D box
  𝒢 = [M/2, M]². Each value n has mass |𝒢_P(n)|/|𝒢| ≤ 16/(MV). This gives rational smoothness.
- **(GM) GM-U, one-point smoothing.** For fixed P, F_A(P) = F̂(P) + G(P)⟨c, v(P)⟩ ranges over an AP mod G(P) as c ranges
  over [−H, H]^d. This gives algebraic smoothness via Lemma 8.5′.
- **(K) Kern-U.** For P ≠ P′, the pair (F_A(P), F_A(P′)) is near-uniform over c (a dimension-free fibre count). This is
  LD*'s Step 1.
- **(W) Lemma W.** F̂(g₀, g₁) = Nw with w a random prime spreads the forced residue F_A(P) ≡ wN(b/g₁)^d mod G(P). This
  is LD*'s Step 4.

### 2.1 Lemma HR (height rigidity at the special root) [proved]

Let f₀ ∈ Z[x] be irreducible of degree d with f₀(m) ≡ 0 (mod N). Let f ∈ Z[x] have degree δ ≥ 1 with f(m) ≡ 0
(mod N), and suppose f₀ ∤ f in Q[x]. Then ‖f‖₂ ≥ (N/‖f₀‖₂^δ)^{1/d}.

*Proof.* f₀ is irreducible and does not divide f, so the two are coprime. Apply Lemma RES. ∎

**Consequences for N ∈ 𝔈 with f₀ = f_d and δ ≤ d.** We get ‖f‖₂ ≥ N^{1/d}/(√2·r^d|s|). This is N^{(1−o(1))/d} when
d² log r + d log|s| = o(t). So any polynomial with the root m_d that is not a multiple of f_d, including any other
defining polynomial of K_d, has height ≥ N^{(1−o(1))/d}. Its algebraic values then have *typical* size ≥ N^{1/d}V^d.

**Scope (codex).** This is uniqueness at a fixed known root and bounded degree. The other special degrees use other
roots m_{d′}; they are available, and they are among the polynomially many pairs covered by Prop RB.

**Exhaustive check.** `rigidity.exhaustive_small_height` enumerates all f of degree ≤ 3 and height ≤ 3 with f(m) ≡ 0.
- For N = 2^{30} − 3 it finds only ±f₀.
- For N = 509, where HR guarantees nothing, non-multiples appear, all consistent with the bound.

### 2.2 Lemma PH (perturbation height) [proved; codex snfsrig_1 §2]

For m > 1 and 0 ≠ c ∈ Z[x], ‖(x − m)c‖_∞ ≥ (m − 1)‖c‖_∞. Hence ‖f₀ + (x − m)c‖_∞ ≥ (m − 1)‖c‖_∞ − ‖f₀‖_∞.

*Proof.* Pick i with |c_i| = ‖c‖_∞. The x^i-coefficient of (x − m)c is c_{i−1} − mc_i (with c_{−1} := 0), so its
absolute value is ≥ m|c_i| − |c_{i−1}| ≥ (m − 1)‖c‖_∞. ∎

So GM/Kern-type perturbations of the special field cost height ≈ H·N^{1/d} at once. The resulting pair satisfies (M1):
(d + 1)t/d + 2t/d ≥ t.

### 2.3 Lemma RS (rational-side rigidity) [proved]

Let 0 ≤ m < N and M ≥ 1 with (m + 1)M < N. Let (g₀, g₁) ∈ Z² satisfy g₀ ≡ mg₁ (mod N) and |g₀|, |g₁| ≤ M. Then
g₀ = mg₁. So the only primitive pairs are ±(m, 1), and G_g(P) = ±(a − bm).

*Proof.* N divides mg₁ − g₀, and |mg₁ − g₀| ≤ (m + 1)M < N. ∎

**Consequences.**
- The 2D-random rational side (R) with the special root would need M ≥ N/(m + 1) ≈ N^{1−1/d}. That gives rational
  values ≥ N^{1−1/d}V/ℓ and θ ≥ 1 − 1/d.
- A random root m′ (BLP base-(m₁, m₂)) restores (R) but discards f₀; that is GNFS.
- (My draft's remark on Cunningham roots mζ^i was wrong and is deleted. Those roots are powers of r and give the other
  special fields, which RB covers; both reviewers.)

### 2.4 Prop DI (device invariance for the listed devices) [proved]

Let F₀(a, b) := b^d f₀(a/b) and G₀(a, b) := a − mb. For each device below, the pair (F, G) satisfies
(F, G)(P) = (c₁F₀(TP), c₂G₀(TP)) for all P ∈ Z², with T ∈ M₂(Z), det T ≠ 0, and c_i ∈ Q^×. The field is Q(α₀)
via the corresponding Möbius image of α₀. Compatibility with ψ needs the determinants and denominators to be units
mod N; a gcd check handles them.

| device | T | (c₁, c₂) |
|---|---|---|
| multiplier N·w, kN (f = wf₀, same root) | id | (w, 1) |
| root twist / scaling by u (root um, f = u^d f₀(x/u)) | diag(1, u) | (1, 1) |
| shift f₀(x + t) (root m − t) | [[1, t], [0, 1]] | (1, 1) |
| Möbius (γx + δ)^d f₀((αx + β)/(γx + δ)), G := G₀∘T | [[α, β], [γ, δ]] | (1, 1) |
| same d, exponent j + dh: r^{j+dh}N = (r^h m)^d − s r^{j+dh} | diag(1, r^h) | (1, 1) |

*Proof.* Each identity is the substitution x = TP in the homogenised forms. The table was corrected after codex
snfsrig_2 and the Opus critic: the last row has T = diag(1, r^h), c₁ = 1, and only j′ ≡ j (mod d) with j′ ≥ j. Other
residues of j may give different fields (they need not: for r = 4, s = 3, d = 2, x² − 3 and x² − 12 define
the same field; codex snfsrig_3). `rigidity.device_identities` checks the first four rows symbolically. ∎

**Consequences.**
- On a box Ω, a device's values are those of the fixed pair (F₀, G₀) on TΩ (scaled). A non-unimodular T changes the
  sampled sublattice, and a shift changes the region.
- For a prime w > B every nonzero value of wF₀ is non-smooth. After removing the content nothing changes (codex).
  Hence N·w gives no new randomness, and Lemma W has nothing to act on.
- **Not exhaustive.** A generator change such as x⁵ − s ↦ x⁵ − s² (root m², α ↦ α²) keeps the field but is *not*
  Möbius. Its degree-one values are degree-two elements in α₀ (codex snfsrig_2). Such changes, and all other
  polynomially many special pairs, are covered by Prop RB, not by DI.

### 2.5 Prop RB (relativisation barrier for small fixed families) [proved, for a listed class of inputs]

**Setting.**
- A one-level design 𝔇 is given by Φ fixed pairs (f_i, g_i), with devices allowed, and special-q targets ℓ ∈ J.
  Each target is searched in a box of ≤ W lattice points; [U] has W ≤ 400S².
- For 𝒮 ⊂ Z, let γ_𝔇(𝒮) be the descendability density of [U] §1 computed with 𝒮 in place of the B-smooth integers
  𝒮_B. The filter keeps the rational cofactor in [X_r/C, X_r].
- Let 𝒱 be the set of all rational cofactors G_i(P)/ℓ over all (i, ℓ, P). Then |𝒱| ≤ Φπ_JW ≤ N^{(2+o(1))ν} whenever
  ΦW = N^{o(ν)}.

**Input class 𝒦.** The class consists of the following statements, applied to the integer set 𝒮 at sizes x ≥ X_r/C.
- (ii) linear upper bounds #(𝒮 ∩ I) ≤ U(I) (DS-U-type counts in progressions; Rankin).
- (iii) linear two-sided estimates Σ_{I∈ℐ} w_I|#(𝒮 ∩ I) − main(I)| ≤ err, with weights |w_I| ≤ 1 and a *fixed* main(I).
  Examples: HT local estimates, Lemma U1, Ψ in a single AP, and BV-type sums with Σ over moduli t ≤ Q₀ of max over
  residues, as in Lemma 8.5′(R′). Here |ℐ| counts one slot per modulus, so |ℐ| ≤ Q₀.
- **Not in 𝒦** (codex snfsrig_3):
  - quadratic or pair statements such as L2_fam and Barban–Davenport–Halberstam sums over all residues, which have
    ≍ Q₀² slots, and Fourier moments, which are not monotone under deletion;
  - multiplicative closure of 𝒮, which deletion destroys.
- So RB does **not** say that LD*'s full input set fails. It says that the linear inputs in 𝒦 cannot carry a supply
  proof on their own.

**Prop RB.**
- (i) γ_𝔇(𝒮_B ∖ 𝒱) = 0.
- (ii) Every type-(ii) input holds for 𝒮′ := 𝒮_B ∖ 𝒱.
- (iii) Every type-(iii) input holds for 𝒮′ with err replaced by err + |𝒱|·|ℐ|.

Hence a proof of γ_𝔇 ≥ γ₀ > 0 from inputs in 𝒦, valid with some slack (each type-(iii) input holding with error
err − |𝒱|·|ℐ|), would also prove γ_𝔇(𝒮′) ≥ γ₀. That contradicts (i). Such a proof therefore cannot exist.

*Proof.* (i): every success requires some G_i(P)/ℓ ∈ 𝒮, and 𝒱 ∩ 𝒮′ = ∅. (ii): monotonicity. (iii): for each I,
|#(𝒮_B ∩ I) − #(𝒮′ ∩ I)| ≤ |𝒱|. ∎

**SNFS numerology [proved arithmetic].**
- *Sizes.* Φ = (log N)^{O(1)} covers the K_d, the generator changes and devices with entries (log N)^{O(1)}. With
  π_J ≤ 2Y = O(N^{2ν}) and W = L^{O(1)}, we get |𝒱| ≤ N^{(2+o(1))ν}. The values have size x ≥ X_r/C = N^{1/d−ν−o(1)}.
- *Single sets.* Single intervals or APs of expected size ≥ xρ/B^{O(1)} absorb |𝒱| when 3ν < 1/d.
- *Moduli sums.* BV-type sums over Q₀ ≤ x^{1/2} moduli (one max-over-residues slot each) have budget ≍ εΨ(x) with
  ε = L^{−O(1)}. The deletion error is ≤ |𝒱|x^{1/2}. With x ≥ N^{1/d−ν−o(1)} the margin is
  N^{(1/d − 5ν)/2 − o(1)}, which beats L^{O(1)} when **5ν < 1/d** with room, e.g. ν ≤ c/d with c < 1/5 fixed. At the 2/5
  design νd → 0, so the margin is N^{(1−o(1))/(2d)}.
- *Regime.* This holds at the 2/5 design (d ≍ D^{1/5}, ν ≍ D^{−2/5}) and at every fixed d once k > 5d/4.
- *Not checked.* We did not check every individual input of FM/LD* in its exact indexed form (codex snfsrig_2,
  snfsrig_3). The claim is for the class 𝒦 as defined.

**Contrast with Theorem U [consistency check, corrected].** In [U] the rational value G_g(P) of a fixed P has point
mass ≤ 16/(MV) over g (FIB), and the algebraic value has point mass ≤ 1/((2H + 1)V′^{d−1}) over c (SM1). So deleting
any K values changes FM's main term by at most K·#𝒫·16/(MV) on the rational side. That is negligible for
K ≤ N^{μ+ν}ρ_r·L^{−O(1)}, and similarly for K ≤ N^{h+(d−1)ν−o(1)} on the algebraic side. The same statistical inputs
therefore do suffice there. What matters is quantitative spreading of the queried values, not support alone. This is
[nest] Prop SI(i)/(ii), applied to the plain smooth set.

**What a proof would need instead.** It would need form-specific, non-relativising information about Norm(a − bα₀) in
a Kummer field, or about a − b·r^q in base r. Examples are type I/II sums for incomplete norm forms (Heath-Brown,
Maynard), BBDT's methods, and prescribed-digit methods (Bourgain). §4 and §5 compare the ranges.

### 2.6 Device table (what survives)

"✗": the ingredient's hypotheses fail. "✓": it holds.

| device | (R) | (GM) | (K) | (W) | FM | LD* | keeps SNFS norms? |
|---|---|---|---|---|---|---|---|
| (a) random rational side, special root | ✗ (RS) | — | — | — | ✗ | ✗ | only with g = ±(x − m) |
| (a′) random root m′ (BLP) | ✓ | needs new f | — | — | GNFS | GNFS | **no** |
| (a″) root twist / scaling | ✗ (DI) | ✗ | ✗ | ✗ | ✗ | ✗ | yes |
| (b) kN, N·w | ✗ | ✗ | ✗ | ✗ (nothing to spread; w > B kills smoothness) | ✗ | ✗ | yes |
| (c) f(x + t), Möbius | ✗ (DI) | ✗ | ✗ | ✗ | ✗ | ✗ | heights grow like t^d |
| (d) the polynomially many special pairs (all d, generator changes) | ✗ (RB) | ✗ | ✗ | ✗ | ✗ (RB) | ✗ | yes |
| (e) random special-q only | ✗ (RB: union over ℓ) | ✗ | ✗ | ✗ | ✗ (RB) | ✗ | yes |
| (f) f₀ + (x − m)c, fixed g | ✗ (RS) | ✓ | ✓ | ✗ (Res fixed = r^jN) | ✗ | ✗ | **no** (PH) |
| (f′) F̂_{g,w} + G·R_c with random g and w | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **no**: this is Theorem U's design |

Summary:
- [proved, for the listed devices] Every device in this table that keeps ‖f‖ = N^{o(1/d)} fails every ingredient: by
  DI its values are scaled values of (F₀, G₀) on a transformed region, or its family is a small fixed family to which
  RB applies (for the linear inputs 𝒦). Every device that restores the ingredients pays height
  ≥ N^{(1−o(1))/d} (HR, PH) or a GNFS root (RS).
- [heuristic, (M1)] Unlisted small-height pairs with known roots do not exist beyond polynomially many.

---

## 3. Special families with natural rigour

### 3.1 Prop CY (cyclotomic pieces) [proved, conditional on Theorem U]

- r^e − 1 = ∏_{a|e}Φ_a(r) and r^e + 1 = ∏_{a|2e, a∤e}Φ_a(r).
- By Möbius inversion, log Φ_a(r) = φ(a)log r + Σ_{b|a}μ(a/b)log(1 − r^{−b}). The sum is bounded in absolute value by
  Σ_{b≥1} −log(1 − 2^{−b}) < 1.25. So **log Φ_a(r) ≤ φ(a)log r + 1.25** (critic's fix; my draft used
  (r + 1)^{φ(a)}, whose ratio log(r + 1)/log r could push ε above 1).
- For a | e, φ(a) ≤ φ(e). So every piece has log-size ≤ ε_e·t + 1.25 with **ε_e := φ(e)/e** (resp. φ(2e)/e for
  r^e + 1). Since t = e log r − O(1), we have ε_e ≤ 1 + O(1/t).

**Prop CY.** Compute the pieces with `cyclo.algebraic_pieces` (polynomially many integer operations; N is never
factored). Then factor each completely with Theorem U's FACTOR, which includes perfect-power detection. This factors
N = r^e ± 1 in expected time L_N[2/5, c_U·ε_e^{2/5} + o(1)].

*Proof.*
- There are τ(e) ≤ e = O(t) pieces.
- Theorem U gives, for every η > 0, a constant C_η with T_U(n) ≤ C_η + L_n[2/5, c_U + η].
- L_n[2/5, c] is increasing in n. With log n ≤ ε_et + 1.25 ≤ (1 + o(1))t:
  L_n[2/5, c] ≤ exp(c(ε_et + 1.25)^{2/5}(log 2t)^{3/5}) = L_N[2/5, cε_e^{2/5} + o(1)].
- Pieces sharing a prime (only primes dividing e can) are merged at the end. Correctness follows from Theorem U's AKS
  certification. ∎

**Corollary CY-P.** Along e = p_k# (primorials), Mertens gives ε_e = ∏_{p≤p_k}(1 − 1/p) ≍ 1/log p_k ≍ 1/log log e.
So N_e = r^e − 1 factors in L_{N_e}[2/5, o(1)], with o(1) ≍ (log log log N)^{−2/5}, which decays extremely slowly. This is an
infinite natural family with a strictly better rigorous *constant* and the same α.
`cyclo.largest_piece_log_ratio` gives 0.31, 0.234 and 0.208 for e = 30, 210 and 2310 (r = 2).

**Corollary CY-LP [proved; published inputs only].** The same reduction with Lenstra–Pomerance 1992 (L_n[1/2, 1 + o(1)])
on the pieces gives L_N[1/2, ε_e^{1/2} + o(1)]. Hence L[1/2, o(1)] along primorial e. This does not depend on the
project's internal chain (suggested by the Opus critic).

**Aurifeuillian factors.** Examples: 2^{4k+2} + 1 = (2^{2k+1} − 2^{k+1} + 1)(2^{2k+1} + 2^{k+1} + 1); in general
Brent 1993 and Stevenhagen 1987. Where they apply, Φ_{2e}(r) splits into two halves of equal size, which multiplies the
constant by 2^{−2/5} (2^{−1/2} for LP92). `cyclo.special_pieces` implements r = 2.

**No exponent from this route.**
- φ(e) ≥ e/(e^γ log log e + 3/log log e) for e ≥ 3 (Rosser–Schoenfeld 1962, Thm 15, with 2.50637 rounded up), so the primitive piece always has log-size ≥ t/O(log log t).
- So applying a general algorithm to the cyclotomic and bounded-degree Aurifeuillian pieces cannot lower the fixed
  exponent through this size calculation.
- This is not a factoring lower bound. Φ_e(r) is itself a special number; it is just not one our tools exploit.
- For p | Φ_e(r) with p ∤ e one knows e | p − 1, which is only log e bits of information.

### 3.2 Families defined by their factors (codex snfsrig_1 §4; literature subagent) [proved or literature; for completeness]

These families are rigorously sub-2/5. That is why the verdict is restricted to special *forms* without factor
promises.
- Semiprimes with q − p ≤ (log p)². One Fermat step solves them in polynomial time, and there are infinitely many
  (by the PNT).
- Squarefree N with every p − 1 being L[β]-smooth. Pollard p − 1 with random bases gives expected L[β].
- N = p^r q with r ≥ (log p)^{1/3+δ}. Boneh–Durfee–Howgrave-Graham (Coppersmith/LLL) is deterministic and rigorous and
  gives L[α < 1/2] (§5). The time bound exp(O(log p/r)) is recalled from the abstract and is not verified by us.

---

## 4. What the model identifies as missing

Take the Eisenstein family 𝔈 (r, s fixed). §1.4 and §1.5 identify the following inputs as missing for a rigorous
L[1/3] in the Dixon-split + rational-only-descent architecture. **This is a list of missing inputs, not a proved
implication** (codex snfsrig_2 #7: the draft's "would suffice" had a compatibility gap, fixed below).

**(O1) Hypothesis RT (rational-target supply in the Kummer fields) [open].** Fix a set 𝒢 of primes (the recognisable
targets).
- *Constants.* There are constants (b, b₁, σ_b, λ, φ_F).
- *Targets.* The statement is for all large N ∈ 𝔈 and every prime q ∈ 𝒢 ∩ (e^{bE}, e^{b₁Q}] not dividing rsN.
- *Box.* With d = d(log q) as in §1.4, take the reduced-basis box of radius e^{σ_bE} of {(a, b) : a ≡ bm_d (mod q)}.
- *Good point.* The box contains a primitive (a, b) such that:
  - F_d(a, b) = a^d − sr^{j_d}b^d is e^{bE}-smooth;
  - every prime factor of (a − bm_d)/q is ≤ e^{bE}, or lies in 𝒢 ∩ (e^{bE}, q^λ].

The second condition is the **compatibility condition**: RT's large children must again be recognisable (codex
snfsrig_2 #7, FATAL to the draft's sufficiency claim). It is vacuous when 𝒢 = all primes.

RT is the rational-target-only, special-field case of [nest] Hypothesis AD. In the random-integer model it holds with
room when 𝒢 has density 1 − o(1): the good-point margin is 5/6 for regular lattices and ≈ 1/2 in the worst case, and
the union bound runs over e^{O(Q)} targets.

**(O2) Recognition of 𝒢.** A rigorous procedure that finds every prime factor p ∈ 𝒢 of its input in expected time
L_p[1/2, O(1)], for 𝒢 of density ≥ (log x)^{−O(1)} in every dyadic range up to L_N[2/3].
- Unconditionally this is open (§5).
- **Under RH, [g2rh] gives it with 𝒢 = all primes.** Then the compatibility condition is vacuous, and only (O1)
  remains.

**Also needed.**
- Character columns at growing degree: ERH via [nest] Lemma X_ERH, or LV §6 unconditionally at cost L[1/3, Θ(δ)].
- LV6 with one algebraic block per K_d ([nest] Prop R).
- Multi-field square roots and the row algebra (Lemma RW, written in [nest] for a single field).

These are UNVERIFIED bookkeeping. Note that [nest] Theorem NR cannot simply be specialised to a point-mass field law,
because its recursive calls need the algorithm for all integers, which is (O2) again.

**How far (O1) is from known results.**
- **Algebraic half.**
  - RT needs F_d(a, b) = a^d − sr^jb^d, of degree d growing up to O(D′), to be y-smooth with
    log y ≍ (log|F|)^{1/2+o(1)}, so that u → ∞. And it needs this *pointwise* in special-q boxes.
  - BBDT 2012 prove Ψ_F(x, y) ≍ x² for fixed degree-g forms only for y ≥ x^{g−2+ε} (x^{1/√e+ε} for cubics), with x the
    variable size: power-size y and bounded u.
  - Barbulescu–Lachand 2017, Thm 4.2, go down to y ≥ exp((log log x)^{5/3+ε}), but only for fixed primitive imaginary
    *quadratic* forms of fundamental discriminant, in form-adapted regions, and without joint rational smoothness.
  - None of these results establishes the required uniform statement.
- **Rational half.**
  - The values a − b·m_d with |a|, |b| ≤ V lie in the integers below m_dV whose residue mod m_d = r^{q_d} lies in
    [0, V] ∪ [r^{q_d} − V, r^{q_d}). Up to endpoint corrections, a block of base-r digits in positions [log_r V, q_d)
    is constant (all 0 or all r − 1).
  - The prescribed proportion is (log m − log V)/(log m + log V) → 1. At the one-level design this is (1 − dν)/(1 + dν);
    at α = 1/3 it is ≈ (D′ − φ_Fb)/(D′ + φ_Fb).
  - RT tests (a − bm)/q, and division by q destroys the block. The relevant problem is therefore *multiples of q* with a
    prescribed block whose cofactor is smooth.
  - Bourgain (2013, 2015) and Swaenepoel (2020) prescribe only a small positive proportion of the digits of primes.
  - None of these establishes even the rational half. Prop RB shows that statistics of smooth numbers alone cannot
    replace it.

## 5. Literature

Sources:
- the literature subagent (web; access level marked: [read] primary, [abs] abstract or publisher page only,
  [sec] secondary);
- codex snfsrig_1 §5.

The Opus literature critic's check of these citations is still outstanding; see §6.1.

- **SNFS is heuristic.**
  - Lenstra–Lenstra–Manasse–Pollard, "The number field sieve", STOC 1990, 564–572 [read, preprint]. They give
    L[1/3, (32/9)^{1/3}] for r^e ± s under a random-integer heuristic, and say verbatim: "Unfortunately we are unable to
    give a rigorous proof … this paper does not contain a rigorous mathematical result."
  - LNM 1554 (1993): the expanded version, likewise heuristic [sec].
  - No later rigorous SNFS analysis was found.
- **Lee–Venkatesan**, J. Number Theory 187 (2018) 92–159 [read, arXiv 1805.08873].
  - Thm 2.1: congruences of squares (possibly trivial) in L[1/3, (64/9)^{1/3}].
  - Thm 2.3: non-triviality only for p ≡ q ≡ 3 (mod 4), under their Conjecture 7.1.
  - The polynomials are random with height L[2/3, κ]. This is exactly the randomisation SNFS gives up, so special forms
    are not covered.
- **Special-prime discrete logs: all heuristic.**
  - Gordon, SIAM J. Discrete Math. 6 (1993) 124–138 [abs]: heuristic L_p[1/3, 3^{2/3}].
  - Semaev, Math. Comp. 71 (2002) 363–377 [abs, + codex read]: primes ab^t + c, heuristic L_p[1/3, (32/9)^{1/3}].
  - Joux–Lercier, Math. Comp. 72 (2003) 953–967 [sec, + codex read].
  - Joux–Pierrot, Pairing 2013 (LNCS 8365) [abs].
  - Fried–Gaudry–Heninger–Thomé, Eurocrypt 2017 [abs]: experimental.
  - Schirokauer, Math. Comp. 79 (2010) 583–602, low-weight integers [abs]: conjectural.
  - Schirokauer, MSRI Publ. 44 (2008) survey [read]: these times are "conjectured".
- **Rigorous index calculus.**
  - Prime fields: Pomerance 1987 gives L_p[1/2, √2] [sec].
  - F_{p^n}: Lovorn Bender–Pomerance 1998 give L[1/2, √2] [sec].
  - Kleinjung–Wesolowski, JAMS 35 (2022) 581–624 [read]: rigorous quasi-polynomial time in fixed characteristic, and
    L_{p^n}(α + o(1)) when p = L_{p^n}(α), α ∈ (0, 1/2). This is the only rigorous L[α < 1/2] for a natural infinite
    family that we found. It concerns extension fields of small characteristic, not prime fields or integers.
  - No rigorous L[1/3] for prime fields, special or not, was found.
- **Rigorous factoring of special N.**
  - Families with structured *factors* are easy: Boneh–Durfee–Howgrave-Graham, CRYPTO '99 [abs], factor N = p^r q in
    polynomial time when r ≈ log p.
  - Families with structured *form* get only algebraic splittings: Brent, Math. Comp. 61 (1993) 131–149, and
    Stevenhagen 1987 [sec], i.e. cyclotomic and Aurifeuillian pieces.
  - For Cunningham numbers nothing rigorous beats the general L[1/2, 1 + o(1)] (Lenstra–Pomerance 1992), apart from
    this project's internal Theorem U.
- **Friable values of fixed forms** (the (O1) gap).
  - Balog–Blomer–Dartyge–Tenenbaum, Comment. Math. Helv. 87 (2012) 639–667 [read]:
    - Thm 1: Ψ_F(x, y) ≍ x² for y ≥ x^{g−2+ε} (degree g);
    - Thm 2: y ≥ x^{1/√e+ε} for cubics.
  - Dartyge–Martin–Tenenbaum 2001 [sec]: power-size thresholds.
  - Lachand (thesis; Ann. Inst. Fourier 68 (2018)) [sec/abs]: degree 2 down to exp((log x)^{2/5+ε}); degree 3 with an
    error that is non-trivial only near power-size y.
  - Mounier, Acta Arith. 220 (2025) [read]: effective versions, still power-size.
  - Barbulescu–Lachand, Math. Comp. 86 (2017) 397–418, Thm 4.2 [read; confirmed]: for primitive irreducible imaginary
    quadratic forms of fundamental discriminant, smooth-value asymptotics uniformly for
    exp((log log x)^{5/3+ε}) ≤ B ≤ x(log x)^{−κ}. Their Thm 4.1 counts smooth ideals in any number field. Fixed degree 2
    only.
  - **For degree ≥ 3 nothing is proved at y = exp((log x)^c) with c < 1.** RT needs growing degree, special-q lattices,
    joint smoothness with the rational form, and pointwise boxes.
- **Prescribed digits.** Bourgain, Israel J. Math. 194 (2013) and 206 (2015) [abs]: primes with up to cn prescribed
  binary digits, for a small unspecified c > 0. Swaenepoel, Proc. LMS (2020) [sec]: any base, an explicit proportion
  c(g) < 1. The SNFS rational side needs a proportion → 1, for smooth numbers rather than primes.
- **Recognition ((O2), Prop REC).**
  - Lenstra, Ann. Math. 126 (1987) 649–673 [read]:
    - unconditional: Cor. 2.8, success probability ≥ 1 − e^{−h f(w)/(c log v)}, where f(w) is the density of
      w-smooth integers in the Hasse interval of p;
    - conditional: the L_p[1/2, √2] running time (Conjecture 2.10), which needs smooth numbers in √p-intervals.
  - Pomerance 1987, Thm 2.1 [sec; text not obtained]: a rigorous ECM smoothness procedure that recognises all but an
    exceptional set (Granville, MSRI 44 (2008) §2.6, [read]). We could not verify the exceptional-set bound, nor whether
    it covers primes of size L[2/3] with w = L_p[1/2]. It is the natural candidate for (O2) with a target set
    G = J ∖ (bad primes). Codex snfsrig_2 reports
    that the explicit version in Lenstra–Pomerance 1992 §7 runs in (log n)^{O(1)}exp(O((log y)^{6/7})) on
    recognisable y-smooth integers, with a large recognisable prime set. That is β = 6/7, not 1/2.
    **Unconditionally (O2) remains open.**
  - **Under RH:** the project's [g2rh] Theorem R (RH + Younis 2024, unrefereed) gives L_p[1/2, 2 + o(1)]·(log N)^{O(1)}
    for the least prime factor p of any N. That is (O2) with 𝒢 = all primes, internal.
  - Lenstra–Pila–Pomerance, Phil. Trans. R. Soc. A 345 (1993) 397–408 and Proc. LMS 84 (2002) 105–146 [read]: HECM
    finds p in L_p[2/3, c₀](log n)². The final assembly is deferred to a part III that Pomerance's slides call
    "to be written". **The rigorous L_y[2/3] HECM bound used in Prop REC(b) and lemmaM (c) is therefore announced but
    not completely in print.** That makes REC's conclusion stronger: even fewer proved recognisers exist.
  - Short intervals: Matomäki–Radziwiłł 2016 (x^ε-smooth in [x, x + c√x]) [abs]; Croot 2007 [abs]; Hildebrand–Tenenbaum
    almost-all-intervals (exp((log x)^{5/6+ε})) [sec]. None gives L[1/2]-smooth numbers in every Hasse interval.
- **Explicit statements that rigorous NFS is open.**
  - LLMP 1990 (above).
  - Lee–Venkatesan 2018: "Even assuming standard conjectures (e.g. GRH), there is no analysis that any substantial part
    of the NFS will halt."
  - Pomerance and Schirokauer, MSRI 44 (2008).
- **Search conclusion (not a novelty claim).** We found no prior statement of Lemma HR/RS/PH as obstructions to
  randomising SNFS, of Prop RB, or of Prop R2-S. Lemma RES is classical (resultant/Hadamard). The cyclotomic reduction
  is classical; its L[2/5, o(1)] consequence depends on Theorem U.

---

---

## 6. Reviews

Every valid fix is applied above. Records: codex `codex/snfsrig_2.out`; the Opus critics' reports are summarised here.

**GPT-6.1-sol blind round (`codex/snfsrig_1`, before my draft).** It reached the same verdict independently:
- no rigorous α < 2/5; Theorem U is the fallback;
- one level 2/5; multi-level 1/3;
- "Prop R2's 2/5 barrier disappears" (rational-only per-node degree gives 1/3).

It contributed:
- the j-sign correction, the Eisenstein family, Lemma PH, the RES-based uniqueness at a fixed root (= HR) and the RS
  threshold;
- the **recognition qualification** (u₀, v₀ are not special), which became Prop REC;
- the factor-defined fast families (§3.2);
- the Barbulescu–Lachand exception and Mounier 2025.

**GPT-6.1-sol adversarial referee (`codex/snfsrig_2`).** Verdict: "GAP — major revision". The model calculations
survive; the universal rigidity, the necessity of recognition and the relativisation conclusions exceeded their proofs,
and the RT + (O2) sufficiency lacked a compatibility hypothesis. All points are accepted and applied:
1. DI's last row was wrong; it is now diag(1, r^h) with c₁ = 1. "Same values" became "scaled values on the transformed
   region". The Minkowski remark is deleted. The unit/gcd condition for ψ is added.
2. UF: "RES forces p ≤ 1" became "RES lets us apply UF with p = 1". The scope is restricted to the Dixon-split
   architecture. Theorem U's θ formula in §1.2 is fixed (the draft's 3/(d+1) + (d+1)/(4k) would give √3, not √2).
3. R2-S: the top degree is 7, not 8, with slack after rounding. Below 1/3 the claim is now about typical nodes. The
   dilemma paragraph's exhaustive classification is withdrawn.
4. REC: "guarantees do not establish", not lower bounds. The general threshold α ≥ β/(1 + β) is added. Pomerance's
   averaged ECM (explicit exp(O((log y)^{6/7})), LP92 §7) is included. The p − 1 density is now via the upper bound
   Ψ(2Y, B).
5. RB: it is now stated as a meta-theorem over a listed input class 𝒦, with slack. The bound is
   |𝒱| ≤ N^{(2+o(1))ν} under ΦW = N^{o(ν)}. The rational-side spreading in U is explained (the contrast had mismatched
   marginals). The check is not claimed for every indexed input.
6. CY: ε_e = φ(e)/e via the Möbius bound. "No exponent gain possible" is restricted to this route.
7. §4 (FATAL to the draft's sufficiency claim): RT now requires large children in 𝒢. The section is relabelled "missing
   inputs, not an implication". The reduced-basis envelope is treated, s is renamed σ_b, and the digit description is
   made an approximation with the division-by-q caveat.
8. The family is now r, s fixed and e → ∞. DI is not exhaustive (the generator change x⁵ − s → x⁵ − s² goes to RB). The
   verdict is restricted to forms without factor promises.

**GPT-6.1-sol re-check (`codex/snfsrig_3`, of the revised note).** Verdict: "GAP — further revision, principally for
RB's application"; no new FATAL; the draft's FATAL compatibility gap is FIXED.
- Items FIXED: 1, 2, 4, 6 and 7.
- Items PARTLY fixed: 3, 5 and 8. All are now applied:
  - worst-case rounding: d = 3, so λ ≤ 19/48 + ε;
  - the unsupported "≈ 1/2" margin removed;
  - 𝒦 restricted to *linear* inputs. L2_fam and BDH are quadratic and are excluded; BV slots are counted per modulus;
    the margin is made explicit;
  - the false "different fields" remark corrected (x² − 3 and x² − 12);
  - the summaries in §§0, 2.6 and 7 narrowed to match the propositions;
  - the column reference fixed.
- It reproduced 3.335 / 3.847 / 4.788 / 6.189. It confirmed the RH route via [g2rh] Theorem R and Corollary R2, with
  Younis's interval length x^{3/4}.

**Opus proof critic (independent, read-only).** Verdict: ACCEPT WITH FIXES, no FATAL. It re-derived the arithmetic
(Hadamard bound, (1+p)/(3p+2), λ = 11/64, margin 5/6, 5ν < 1/d, Φ_a(r) bounds) and ran the tests. Fixes applied:
- MAJOR 1: rigidity over-claimed beyond the known root. This coincides with codex #3/#8.
- MAJOR 2: **the RH route via [g2rh]** gives (O2) with 𝒢 = all primes, so under RH only (O1) remains. It also asked for
  REC to be labelled unconditional.
- MAJOR 3: the "would suffice" overstatement, the compatibility of 𝒢, and the missing character input.
- MINOR 1–13:
  - Theorem U's θ formula;
  - the **hypothetical constants must use U-factored halves**: 3.34/3.85/4.79/6.19, not 3.06/…;
  - the ε_e definition (it could exceed 1);
  - the false Cunningham-root remark;
  - the DI duplicate row and last row;
  - HR wording, plus the d log|s| condition;
  - RB: size filter, rational-only count, K range, "six results";
  - REC: batch smoothness is a total cost; Theorem U is internal; LPP III is unwritten;
  - the worst-case lattice margin ≈ 1/2;
  - fixed r, s;
  - wording ("nobody can prove", "universal", assessments);
  - BBDT's cubic exponent, which the literature subagent and codex both read as 1/√e (kept);
  - **Corollary CY-LP** (its suggestion).

**Opus literature critic (independent, read-only).** See §6.1.

### 6.1 Opus literature critic

Launched on the revised note but not yet reported when this note was committed. Its fixes, if any, are outstanding.
Until it reports, the literature items marked [abs]/[sec] in §5 are unverified by a second reader.

---

## 7. Candidate card

- **Mechanism.** Fixed SNFS polynomials f_d = x^d − sr^j, with the degree chosen per node, inside the project's
  Dixon-split + descent architecture. The randomness sources tried were w, kN, twists, root twists, several (d, j),
  generator changes and special-q.
- **Model.** The [rec]/[nest] size model and Theorem U's supply framework.
- **Result.**
  - Size model: one level 2/5 (Prop UF), multi-level 1/3, rational-only 1/3 (Prop R2-S; GNFS has 2/5).
  - Rigorous: best α = 2/5 (Theorem U, internal).
  - Constant: c_Uε_e^{2/5} for r^e ± 1 (Prop CY), hence L[2/5, o(1)] along primorial e.
  - Published inputs only: L[1/2, o(1)] along primorial e (CY-LP).
- **Failure scope.**
  - (O1) Fixed fields. HR, RS and PH: re-randomising at the special root keeps the field or costs height N^{1/d}. DI:
    the listed devices give scaled values of the fixed pair on transformed regions. RB: proofs from linear statistical
    inputs (class 𝒦) fail for any family of N^{o(ν)} fixed pairs.
  - (O2) Recognition (Prop REC). Unconditionally, the proved recognisers give only 2/5 (U, HECM) or 6/13 (averaged
    ECM). Under RH, [g2rh] gives 1/3.
- **Most promising open question.** A pointwise friable-values estimate for the Kummer forms a^d − cb^d on special-q
  lattices at L-scale (Hypothesis RT). Even the rational half alone would be new: multiples of q with a long prescribed
  base-r digit block and a smooth cofactor. Under RH this is the *only* missing analytic input for L[1/3] on 𝔈, apart
  from the multi-field LV6 bookkeeping.
