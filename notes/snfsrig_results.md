# snfsrig: special forms (SNFS) and rigorous sub-L[2/5] factoring — results

Agent `snfsrig`, phase 14, 2026-10-20.

Files:
- This note. No prediction-testing experiment was run, so there is no pre-registration (`notes/snfsrig_prereg.md`
  does not exist). The code only checks proved lemmas and does exact exponent arithmetic.
- Code: `src/snfsrig/sizemodel.py` (exact size model, §1), `src/snfsrig/rigidity.py` (checkers for Lemmas RES/HR/RS/PH,
  Prop DI and the value-set count of Prop RB, §2; diagnostic), `src/snfsrig/cyclo.py` (cyclotomic and Aurifeuillian
  pieces of r^e ± 1, §3; operational-safe, never factors N). Tests: `tests/test_snfsrig.py`.
- Codex: `codex/snfsrig_1.{md,out}` (blind round, launched before my derivations were written), `codex/snfsrig_2.*`
  (adversarial referee of this note, §6).

Base texts: Theorem U and its design (`notes/lv5unif_results.md`, "[U]"), the size model (`notes/recurse_results.md`
§3.3–3.4, "[rec]"), Props SI/R/R2 and Theorem NR (`notes/nested_results.md`, "[nest]"), the w multiplier
(`notes/ldbound_results.md`, "[ld]"), the field law (`notes/thmA_proof.md` §1), HECM barrier (`notes/lemmaM_results.md`
(c)). Everything the project calls "Theorem U" is internal and unrefereed.

Labels: **proved** (a mathematical proof is given here, from published results); **proved (model arithmetic)**
(an exact consequence of the stated size-model premises, which are heuristic); **heuristic**; **open**;
**literature** (with access level).

---

## 0. Verdict

**Plain language.** No. Special forms do not give a rigorous algorithm below L[2/5] with the tools we have, and we can
say exactly why. The tiny SNFS polynomial helps the *size model* a lot: it removes the 2/5 barrier for designs without
algebraic children, so in the model rational-only multi-level descent reaches L[1/3]. Rigorously it does not help,
for two independent reasons. (O1) Every relation-supply proof in the project averages over random fields. For a
special N, every way of re-randomising the field either keeps the same field and the same set of values, or costs
polynomial height ≥ N^{1/d} and so turns SNFS back into GNFS. What is left would need smoothness estimates for the
values of one fixed binary form at L-scale, which nobody can prove. (O2) The Dixon-split architecture must also fully
factor the split halves u₀, v₀. These are arbitrary integers of size √N, not special ones, so the recursion that removes
the recognition cost in general cannot stay inside the special family. With the rigorous recognition tools we have
(the general algorithm, Theorem U, or the hyperelliptic method HECM), this cost alone already forces α ≥ 2/5.
The best rigorous result for a natural special family is therefore still **α = 2/5** (Theorem U applies to every N,
conditional on the internal chain). The special structure improves only the constant: L[2/5, c_U·ε^{2/5}] for
r^e ± 1 via cyclotomic pieces, where ε = φ(e)/e, and hence **L[2/5, o(1)] along N = r^e − 1 with e running over
primorials**.

Precisely:

1. **Size model (§1).**
   - **Lemma RES [proved].** Two coprime f, g ∈ Z[x] with a common root mod N satisfy N ≤ |Res(f, g)| ≤ ‖f‖₂^{deg g}‖g‖₂^{deg f}.
     So every pair with a linear side has coefficient exponent θ_coeff ≥ 1/d. SNFS attains this minimum.
   - **Prop UF [proved (model arithmetic)].** One-level descent has α ≥ 2/5 for every polynomial pair with a linear side,
     special or not. If θ_coeff ≍ d^{−p}, the floor is (1 + p)/(3p + 2). RES forces p ≤ 1, so the floor is 2/5. The
     (R3) floor moves only in its constant, θ ≥ 1/d + (d−1)/(4k) instead of 3/(d+1) + …, so the minimum over d is
     θ ≈ 1/√k for SNFS against √2/√k for both heuristic GNFS and Theorem U's design at 2/5.
   - Hypothetical constants. Suppose SNFS had Theorem U's supply. Then Theorem U's cost formula would give c ≈ 3.06 / 3.51 /
     4.38 / 5.66 against Theorem U's 3.51 / 4.05 / 5.04 / 6.50, for C_γ = 0 / 1 / 5 / 20. This is about 13% lower and is
     *not* a theorem.
   - Multi-level descent in one field: 1/3 in the model, the same as GNFS ([rec] §3.4).
   - **Prop R2-S [proved (model arithmetic)].** Recomputing nested's Prop R2 without its premise (M1): for SNFS N, choose
     at each node the degree d(z) ≍ Ξ/z of a special polynomial. Then designs *without algebraic children* contract
     targets by t^{1−3α}. The floor drops from 2/5 to **1/3**, attained with O(log log N) levels at constant
     contraction. This is exactly where small norms shift R2. GPT-6.1-sol's blind round found the same.
2. **Randomness (§2).** Six results, all **proved**:
   - Lemma HR: a polynomial with the special root mod N is either a multiple of f₀ or has height ≥ (N/‖f₀‖^{δ})^{1/d}.
   - Lemma PH: a perturbation f₀ + (x − m)c costs height ≥ (m − 1)‖c‖ − ‖f₀‖.
   - Lemma RS: a rational side sharing the root with height < N/(m + 1) is a multiple of x − m.
   - Prop DI: w·N, kN, root twists, f(x + t), u^d f(x/u) and Möbius maps give (F, G)(P) = (c₁F₀(TP), c₂G₀(TP)) with
     T ∈ M₂(Z). These are the same field and the same values, so a prime w > B even makes every value non-smooth.
   - Prop RB, a relativisation barrier for fixed field families. A fixed family of Φ fields queries at most
     2Φπ_J(2S+1)²·O(1) = N^{2ν+o(1)} values. Deleting them from the smooth set changes every count used by FM/LD*/8.5′/L2_fam
     by a negligible amount, and makes supply exactly 0. So no proof that uses only statistics of smooth numbers can
     give supply for SNFS. Theorem U escapes because its field law has full support.

   Device table (§2.6): **none** of FM, LD*, Kern, Lemma W or Theorem R survives a device that keeps the SNFS norms.
   The ones that survive (f + cg with fresh g and w) are exactly Theorem U's design, whose norms do not see the special form.
3. **Recognition (§1.5, Prop REC) [proved accounting, given the tools].**
   - Each accepted row needs the complete factorisation of u₀, v₀ ≈ √N. These have prime factors of size L[1 − α] and are
     not special.
   - Rigorous HECM costs L_N[2(1−α)/3], which is ≤ L[α] iff α ≥ 2/5 (lemmaM (c)).
   - Recursion through Theorem U costs L_N[2/5, 2^{−2/5}c_U].
   - So sub-2/5 for a special family also needs a rigorous ECM-speed recogniser for √N-size integers with L[1 − α]
     factors, at least on a target set of polylogarithmic density; this would suffice down to α = 1/3. §5 says what the
     literature gives.
4. **Special families (§3).**
   - **Prop CY [proved, conditional on Theorem U for the pieces].** N = r^e ± 1 factors in L_N[2/5, c_U·ε_e^{2/5} + o(1)],
     where ε_e = φ(e)/e (respectively φ(2e)/e).
   - **Corollary CY-P.** Along e = primorial this is L[2/5, o(1)].
   - Aurifeuillian splitting gives a further factor 2^{−2/5} where it applies.
   - Families defined by their *factors* are trivially easier: close primes are polynomial time, and p − 1 smooth gives
     L[β] (codex). They are not special *forms* of N.
5. **Precise obstruction (§4).** Two statements would suffice for a rigorous L[1/3] on the Eisenstein family
   N = r^e − s (s prime, s ∤ r):
   - (O1) Hypothesis RT, a pointwise rational-target supply for the explicit Kummer forms x^d − sr^j at L-scale;
   - (O2) ECM-speed rigorous recognition.

   On (O1), the best friable-value results for fixed binary forms of degree ≥ 3 need power-size y = x^{c}. Barbulescu–Lachand reach
   y = exp((log log X)^{5/3+ε}), but only for fixed *quadratic* (imaginary) forms (§5). The SNFS rational side alone,
   a − b·r^j, consists of integers with a (1 − o(1))-proportion of prescribed base-r digits. Bourgain-type results
   prescribe only a small positive proportion, so even the rational half of (O1) is out of reach.
6. **Literature (§5).** No rigorous or GRH-conditional sub-1/2 result for factoring special-form N or for prime-field DL
   with special p was found. SNFS (LLMP 1990/93), Gordon 1993, Semaev 2002, Joux–Lercier 2003 and Joux–Pierrot 2013 are
   heuristic. Lee–Venkatesan 2018 is for random polynomials, not fixed special pairs.

**Weakest points.**
- The barrier results are method barriers, not impossibility theorems.
  - Prop RB excludes only proofs that use statistics of the smooth set.
  - Prop REC concerns the Dixon-split architecture with the recognition tools we have.
- Everything positive about α = 2/5 rests on Theorem U's internal chain.

---

## 1. Size model

Notation: t = log N, D = t/log t, Λ_α = t^α(log t)^{1−α}, B = e^{aΛ_α}, U ≍ D^α the largest u with u log u ≤ Λ_α,
Ξ_α := Λ_αU ≍ t^{2α}(log t)^{1−2α}. A value of log-size X is B-smooth at cost L[α]^{O(1)} in the random-integer model
iff X = O(Ξ_α) (CEP). The special family: N = r^e − s, r ≥ 2, s ≠ 0 fixed (or |s| = N^{o(1)}), gcd(r, N) = 1. For
each d ≥ 2 put q_d := ⌈e/d⌉, m_d := r^{q_d}, j_d := dq_d − e ∈ [0, d), and

  f_d(x) := x^d − s·r^{j_d},  so f_d(m_d) = r^{j_d}N ≡ 0 (mod N),  log‖f_d‖ ≤ d log r + log|s|,  log m_d ≤ t/d + log r.

(codex snfsrig_1 corrected the sign of j in my prompt.) **Eisenstein family.** If s is a prime with s ∤ r, then f_d
is Eisenstein at s for every d, hence irreducible. We call 𝔈 := {r^e − s : s prime, s ∤ r} the Eisenstein family;
it is the natural infinite family below (codex's suggestion). For Cunningham numbers r^e ± 1 one uses x^d ∓ r^{j}
with gcd(j, d) = 1 and r not a perfect power, which is irreducible by Capelli's theorem when 4 ∤ d.

### 1.1 Lemma RES and the coefficient floor [proved]

**Lemma RES.** Let f, g ∈ Z[x], deg f = d ≥ 1, deg g = δ ≥ 1, coprime in Q[x], and f(m) ≡ g(m) ≡ 0 (mod N) for some
integer m. Then N | Res(f, g), Res(f, g) ≠ 0, and N ≤ |Res(f, g)| ≤ ‖f‖₂^δ‖g‖₂^d.

*Proof.* The adjugate of the Sylvester matrix gives A, B ∈ Z[x] with A f + B g = Res(f, g) (e.g. Lang, *Algebra*,
IV §8). Evaluating at m gives Res ≡ 0 (mod N). Res ≠ 0 because f and g have no common root in Q̄. Hadamard's
inequality applied to the rows of the Sylvester matrix (δ shifted copies of f's coefficient vector and d of g's) gives
|Res| ≤ ‖f‖₂^δ‖g‖₂^d. ∎

**Corollary (linear side).** If g = g₁x − g₀ and f is irreducible of degree d ≥ 2, then ‖f‖₂‖g‖₂^d ≥ N. With
a := log‖f‖₂ and b := log‖g‖₂ we get a + db ≥ t, hence a + b ≥ t/d + a(1 − 1/d) ≥ t/d. So **θ_coeff :=
(log‖f‖ + log‖g‖)/t ≥ 1/d − O(log d/t)** for every NFS pair with a linear side. SNFS attains equality up to
O(d log r/t).

`rigidity.check_res_floor` checks the divisibility and the Hadamard bound exactly on SNFS pairs (tests).

### 1.2 One level [proved (model arithmetic)]

The premises are [rec] (R1)–(R4), the random-integer model.
- **(R1) Split.** The halves u₀ and v₀ are products of k primes of size N^{1/(2k)}, with k = O(D^α).
- **(R2) Smoothness.** θ = O(D^{2α−1}).
- **(R3) Norm-product exponent.** θ = θ_coeff + (d − 1)ν with ν = 1/(4k). This holds for the box V ≈ √ℓ·S, whose
  rational value is ‖g‖V/ℓ and whose algebraic value is ‖f‖V^d.

For SNFS, θ_SNFS = 1/d + (d − 1)/(4k) + O(d log S/t): the randomisation floors of Theorem U's design (h ≥ 2ν, level,
window) are absent because nothing is randomised.

**Prop UF (universal one-level floor).** Suppose θ_coeff ≥ c·d^{−p} with c > 0 and p ≥ 0. Then α ≥ (1 + p)/(3p + 2)
(and 1/2 for p = 0). By the Corollary of RES, every pair with a linear side has p ≤ 1, so **α ≥ 2/5 for every
one-level design, special or not.**

*Proof.* Write d ≍ D^x and k ≍ D^y. Then (R2) applied to the two terms gives px ≥ 1 − 2α and y − x ≥ 1 − 2α, and (R1) gives y ≤ α.
Eliminating x and y gives (1 − 2α)(1 + p)/p ≤ α, i.e. α ≥ (1 + p)/(3p + 2). At p = 1 this is 2/5. `sizemodel.one_level_floor`, and a brute-force rational LP in
the tests. ∎

So the (R3) floor moves only in its constant. At fixed k the minimum over d is:
- SNFS: θ_min = 1/√k + O(1/k), attained at d ≈ 2√k;
- heuristic GNFS (base m, both heights N^{1/(d+1)}): √2/√k, at d ≈ √(8k);
- Theorem U's design: 3/(d+1) + (d+1)/(4k), whose minimum over d is again √2/√k − in its Kern regime it coincides with
  heuristic GNFS at 2/5 (`sizemodel.theta0`).

Hence θ₀ := lim θD^{1/5} = 1/d₀ + d₀/(4κ) for SNFS, against 2/d₀ + d₀/(4κ) for both others.

**Hypothetical constant.** Insert SNFS's θ₀ into Theorem U's cost formula (c = max(2a, a + s + max(a, σ),
(a + s)/(1 − 2^{−2/5})), s = 2κ(1 + C_γ), σ = 2θ₀/(5a); exact rational grid, `sizemodel.optimise`):

| C_γ | Theorem U (reproduced) | SNFS, *if* it had U's supply |
|---|---|---|
| 0 | 3.510 | 3.056 |
| 1 | 4.047 | 3.510 |
| 5 | 5.038 | 4.378 |
| 20 | 6.498 | 5.658 |

This is not a theorem: the SNFS supply is exactly what §2 shows we cannot prove.

### 1.3 Multi-level descent in one field [model arithmetic, as [rec] §3.4]

With d = D^δ, a level-j prime of log-size t^{β_j} has values of log-sizes t^{1−δ} (rational, ≈ log m) and t^{δ+β_j}
(algebraic). These exponents are the same for GNFS and SNFS heights, because only the constants differ. So the level
count ⌈(1 − 2α)/(3α − 1)⌉ (`sizemodel.multilevel_levels`) and the floor 1/3 (attained with O(log log N) levels at
constant contraction) are unchanged. The obstruction there is nested's: algebraic children are field-bound ([nest]
Lemma FB, Lemma OD).

### 1.4 Prop R2-S: designs without algebraic children [proved (model arithmetic)]

Nested's Prop R2 premises are (M1) (d + 1)log H_f + 2 log H_g ≳ t, (M2) typical sizes, (M3) random-integer
smoothness, and (M4) split targets z₁ ≍ t^{1−α}(log t)^α. Under them, designs whose algebraic values are all B-smooth
have α ≥ 2/5. Nested itself notes that (M1) is only a counting heuristic, and that small-height f breaks it.

**Premise (M1-S).** For every degree 2 ≤ d ≤ D there is an irreducible pair with log H_f ≤ d log r + log|s| =
O(d) and log H_g ≤ t/d + O(1) (the family f_d, m_d above; for 𝔈 every d works). (M2)–(M4) are kept, and each node
may choose its own d.

**Prop R2-S.**
(i) For α > 1/3, a rational-only descent tree exists with a bounded number of levels, ⌈(1 − 2α)/(3α − 1)⌉.
(ii) At α = 1/3 it exists with O(log log N) levels, whose targets contract by a constant factor.
(iii) For α < 1/3 every node expands its targets, so no finite tree exists.
So the floor is **1/3** instead of 2/5.

*Derivation (model arithmetic).*
- **Algebraic side.** A node with rational target q = e^z uses a reduced special-q box of radius e^{sΛ}√q in the field
  of degree d. Its algebraic value has log-size d(z/2 + sΛ) + O(d). B-smoothness within budget needs this to be
  ≤ φ_FΞ, so choose d(z) := ⌊φ_FΞ/(z/2 + sΛ)⌋.
- **Rational cofactor.** Its log-size is t/d − z/2 + sΛ ≈ t(z/2 + sΛ)/(φ_FΞ). For it to be e^y-smooth within budget,
  y ≥ (log-size)/(φ_RU), so y/z ≈ t(1/2 + sΛ/z)/(φ_Fφ_RΞU).
- **Contraction.** Since ΞU ≍ t^{3α}(log t)^{1−3α}, the ratio is ≍ t^{1−3α}(log t)^{3α−1}. That is a polynomial
  contraction for α > 1/3, a constant for α = 1/3 and an expansion for α < 1/3 (`sizemodel.r2s_contraction_exponent`).
  For (iii), the reverse inequality d ≤ φ_FΞ/(z/2) forces a rational size ≥ tz/(2φ_FΞ), hence y/z → ∞ (codex
  snfsrig_1 makes the same argument).
- **At α = 1/3.** Write E = Λ_{1/3}, Q = t^{2/3}(log t)^{1/3}, D′ = Q/E, so t = D′Q.
  - Parameters: B = e^{bE}, box radius e^{sE}√q, u_F = φ_FD′, u_R = φ_RD′, top targets z₁ = b₁Q.
  - Each node then has d(z) = φ_FbQ/(z/2 + sE) ∈ [2φ_Fb/b₁, 2φ_FD′] and children y ≤ λz, with
    λ = (1/2 + s/b)/(φ_Fφ_Rb).
  - Example: b = 4, s = 3/4, φ_F = φ_R = 1, b₁ = 1. Then λ = 11/64, d₁ = 8, and the good-point margin
    2s − (φ_F + φ_R)/3 = 5/6 > 0 (`sizemodel.r2s_alpha_third`).
  - Tree: depth O(log D′) = O(log log N), fan-out O(t/E), so e^{O((log log N)²)} = L^{o(1)} nodes ([nest] Lemma TS).
- The top degree is O(1) and the leaf degree is O(D′): the degree rises down the tree. ∎

**Why GNFS cannot do this.** By RES, a field sharing a root with a small rational side is either the special field
or has log H_f ≥ (t − δ log‖f₀‖)/d (Lemma HR, §2). Fields with log H_f ≳ t/d obey (M1) and fall under Prop R2. So
in the model **every node of a sub-2/5 rational-only tree whose target is large must use a field of height
N^{o(1/d)}**, which for N ∈ 𝔈 means a special field, i.e. one of O(D²) fixed fields. This is the dilemma of §2.
Codex's blind round found the same 1/3 ("the former forbidden window disappears").

### 1.5 Recognition: Prop REC [proved accounting, given the tools]

In the Dixon-split architecture, each accepted row requires the complete factorisation (equivalently, the shape
certificate) of the Gauss-reduced halves u₀, v₀ of a random square. These are arbitrary integers ≤ √(2N), with all
prime factors in the target set ⊂ (Y, 2Y], Y = N^{1/(2k)} = L_N[1 − α]. Every intermediate node of a multi-level tree
similarly has to factor rational cofactors with prime factors up to L[1 − α].

**Prop REC.** Let α < 2/5, and suppose the algorithm uses only these rigorous recognisers: (a) a recursive call of
a general factoring algorithm, of which the best rigorous one available is Theorem U; (b) LPP's hyperelliptic method
(HECM); (c) trial division, Pollard–Strassen and batch smoothness. Then none of them has a *proved* cost bound of
L_N[α]^{O(1)} for certifying the split. (These are upper bounds for the known methods, not lower bounds for
recognition.)
- (a) costs L_{√(2N)}[2/5, c_U] = L_N[2/5, 2^{−2/5}c_U + o(1)] ≫ L_N[α]^{O(1)}.
- (b) costs L_Y[2/3]^{O(1)} = L_N[2(1 − α)/3]^{O(1)}, and 2(1 − α)/3 > α ⟺ α < 2/5 (lemmaM (c)).
- (c) costs Y^{Ω(1)} or π(2Y)^{Ω(1)} = L_N[1 − α]^{Ω(1)}.
- Restricting the target set does not help (Lemma S″). Any target set G ⊂ J of relative density δ multiplies the
  acceptance by about δ^{2k}, and k ≍ D^α. So G must have density (log N)^{−O(1)}. The primes p ∈ J with p − 1
  L[α]-smooth (Pollard p − 1) have density exp(−(1 + o(1))u log u) with u ≍ D^{1−2α}. Then δ^{2k} = exp(−Θ(D^{1−α} log D)),
  which is not L[α]^{O(1)}.

*Proof.* (a) is Theorem U at √(2N) with L_{√x}[α, c] = L_x[α, 2^{−α}c + o(1)] ([rec] §3.1). (b) is LPP's bound with
y = 2Y. (c) is direct. ∎

So in the general case recursion is what removes the recognition cost ([rec] §3.1). It cannot be used inside a
special family, because u₀ and v₀ are not special (codex snfsrig_1 §1, "algorithmic qualification"). A rigorous
recogniser at ECM speed, L_p[1/2, O(1)], for all p in a target set of polylogarithmic density would cost L_N[(1 − α)/2]
and suffice down to α = 1/3. What is known is in §5.

---

## 2. Recovering randomness: what survives

The relation-supply proofs use field randomness in exactly four places ([U] §§4–6, [ld], thmA §1):
- (R) **Theorem R / FIB**: for fixed P, G_g(P) = g₁a − g₀b is near-uniform on an interval as g ranges over the 2D box
  𝒢 = [M/2, M]². This gives rational smoothness.
- (GM) **GM-U, one-point smoothing**: for fixed P, F_A(P) = F̂(P) + G(P)⟨c, v(P)⟩ ranges over an AP mod G(P) as c
  ranges over [−H, H]^d. This gives algebraic smoothness via Lemma 8.5′.
- (K) **Kern-U**: for P ≠ P′, (F_A(P), F_A(P′)) is near-uniform over c (dimension-free fibre count). This is LD*'s
  Step 1.
- (W) **Lemma W**: F̂(g₀, g₁) = Nw with w a random prime spreads the forced residue F_A(P) ≡ wN(b/g₁)^d mod G(P).
  This is LD*'s Step 4.

### 2.1 Lemma HR (height rigidity) [proved]

Let f₀ ∈ Z[x] be irreducible of degree d with f₀(m) ≡ 0 (mod N). Let f ∈ Z[x] have degree δ ≥ 1 with f(m) ≡ 0
(mod N), and suppose f₀ ∤ f in Q[x]. Then ‖f‖₂ ≥ (N/‖f₀‖₂^δ)^{1/d}.

*Proof.* f₀ is irreducible and f₀ ∤ f, so f₀ and f are coprime. Apply Lemma RES. ∎

For N ∈ 𝔈 and f₀ = f_d, with δ ≤ d: ‖f‖₂ ≥ N^{1/d}/(√2·r^{d}|s|). So **any field with the special root other than
K_d itself has height ≥ N^{(1−o(1))/d}**, as long as d² log r = o(t). The algebraic values of such a field have size
≥ N^{1/d}V^d. Codex adds the scope: this is uniqueness at a fixed root and bounded degree. Other special degrees use
other roots m_{d′}, and they *are* available. They are the O(D²) fields of the special family, and Prop RB covers them.
`rigidity.exhaustive_small_height` enumerates all f of degree ≤ 3 and height ≤ 3 with f(m) ≡ 0 for N = 2^{30} − 3 and
finds only ±f₀. For N = 509, where HR guarantees nothing, non-multiples appear, all consistent with the bound.

### 2.2 Lemma PH (perturbation height) [proved; codex snfsrig_1 §2]

For m > 1 and 0 ≠ c ∈ Z[x], ‖(x − m)c‖_∞ ≥ (m − 1)‖c‖_∞. Hence ‖f₀ + (x − m)c‖_∞ ≥ (m − 1)‖c‖_∞ − ‖f₀‖_∞.

*Proof.* Pick i with |c_i| = ‖c‖_∞. The x^i-coefficient of (x − m)c is c_{i−1} − mc_i (with c_{−1} := 0), and
|c_{i−1} − mc_i| ≥ m|c_i| − |c_{i−1}| ≥ (m − 1)‖c‖_∞. ∎

So GM/Kern-type perturbations of the special field cost height ≈ H·N^{1/d} at once. At H = 1 the algebraic norm is
already that of a GNFS pair with rational side N^{1/d}, and (M1) holds: (d + 1)t/d + 2t/d ≥ t.

### 2.3 Lemma RS (rational-side rigidity) [proved]

Let 0 ≤ m < N, M ≥ 1 with (m + 1)M < N, and let (g₀, g₁) ∈ Z² satisfy g₀ ≡ mg₁ (mod N) and |g₀|, |g₁| ≤ M.
Then g₀ = mg₁. So the only primitive pairs are ±(m, 1), and G_g(P) = ±(a − bm).

*Proof.* |mg₁ − g₀| ≤ (m + 1)M < N, and N divides it. ∎

So the 2D-random rational side (R) with the special root needs M ≥ N/(m + 1) ≈ N^{1−1/d}, which gives rational values
≥ N^{1−1/d}V/ℓ and θ ≥ 1 − 1/d. A random root m′ (BLP base-(m₁, m₂)) restores (R) but discards f₀: that is GNFS. For
Cunningham N | r^e − 1 with d | e the other roots mζ^i (ζ = r^{e/d}) are known, but they are generic residues
mod N. Their root lattices have shortest vectors ≍ √N by Minkowski, so this does not help either (remark, not
needed).

### 2.4 Prop DI (device invariance) [proved]

Let F₀(a, b) := b^d f₀(a/b) and G₀(a, b) := a − mb. For each device below there are T ∈ M₂(Z) with det T ≠ 0 and
c₁, c₂ ∈ Q^× such that the device's pair (F, G) satisfies (F, G)(P) = (c₁F₀(TP), c₂G₀(TP)) for all P ∈ Z², and the
device's field is isomorphic to Q(α₀) compatibly with ψ:

| device | T | c₁, c₂ |
|---|---|---|
| multiplier N·w, kN (f = w f₀, the same root) | id | w, 1 |
| root twist λ (root λm, f = λ^d f₀(x/λ)) | diag(1, λ) | 1, 1 |
| shift f₀(x + t) (root m − t) | [[1, t], [0, 1]] | 1, 1 |
| scaling u^d f₀(x/u) (root um) | diag(1, u) | 1, 1 |
| Möbius (γx + δ)^d f₀((αx + β)/(γx + δ)) | [[α, β], [γ, δ]] | 1, (by G₀∘T) |
| N r^i = m′^d − s r^{i′} (other j, same d) | diag(1, r^{i}) | r-power, 1 |

*Proof.* Each identity is the substitution x = TP in the homogenised forms (exact symbolic checks:
`rigidity.device_identities`). The fields are Q(α₀) with the corresponding Möbius image of α₀. ∎

So the value set of a device on a box Ω is the value set of the fixed pair (F₀, G₀) on TΩ ⊂ Z². For a prime w > B,
every nonzero value of wF₀ is non-smooth; after removing content nothing changes (codex). Hence N·w gives no new
randomness, and Lemma W has nothing to act on. A small-coefficient polynomial for r^e − s times w in a *new* field
would need a new known root; by HR it cannot share m.

### 2.5 Prop RB (relativisation barrier for fixed field families) [proved]

**Setting.** A one-level design 𝔇 is given by Φ fixed pairs (f_i, g_i) (devices allowed) and the special-q targets
ℓ ∈ J. Each target is searched in a box of ≤ W points of its lattice; [U] has W ≤ 400S². For a set 𝒮 ⊂ Z let
γ_𝔇(𝒮) be the descendability density of [U] §1 computed with 𝒮 in place of the B-smooth integers 𝒮_B. Let 𝒱 be
the set of all values G_i(P)/ℓ over all (i, ℓ, P). Then |𝒱| ≤ Φπ_JW.

**Prop RB.** (i) γ_𝔇(𝒮_B ∖ 𝒱) = 0.
(ii) Every input of the form "#(𝒮 ∩ I) is bounded above by …" (DS-U, L2_fam, Rankin, Brun–Titchmarsh-type) holds for
𝒮′ := 𝒮_B ∖ 𝒱 because 𝒮′ ⊂ 𝒮_B.
(iii) Every input of the form "Σ_{I∈ℐ}|#(𝒮 ∩ I) − main(I)| ≤ err" (HT local estimates, Lemma U1, Lemma 8.5′(R′),
Tenenbaum's count in APs) holds for 𝒮′ with err + |𝒱|·|ℐ|.

Hence a proof of γ_𝔇 > 0 that uses only inputs of type (ii), and inputs of type (iii) whose error budgets exceed
|𝒱|·|ℐ|, would also prove γ_𝔇(𝒮′) > 0, contradicting (i). Such a proof therefore cannot exist. In the SNFS regime the
budgets do exceed |𝒱|·|ℐ| (next paragraph).

*Proof.* (i): every success requires some G_i(P)/ℓ ∈ 𝒮, and 𝒱 ∩ 𝒮′ = ∅. (ii): monotonicity. (iii): for each I,
|#(𝒮_B ∩ I) − #(𝒮′ ∩ I)| ≤ |𝒱|. ∎

**The SNFS numbers [proved arithmetic].**
- *Size of 𝒱.* Φ = (log N)^{O(1)} (the fields K_d, d ≤ D, together with devices whose T have entries
  (log N)^{O(1)}). Also π_J = N^{2ν+o(1)} and W = L^{O(1)}, so |𝒱| = N^{2ν+o(1)}.
- *Size of the counted values.* The rational values have size X_r = N^{1/d−ν+o(1)}.
- *Inputs on intervals or single APs.* These count sets I of integers of size x ≥ X_r with expected size
  ≥ xρ/B^{O(1)} and budget ε·(main), ε = L^{−O(1)}. Their deletion error |𝒱| is negligible because 1/d − ν > 2ν.
- *Barban–Davenport–Halberstam / Bombieri–Vinogradov-type sums* (Lemma 8.5′). These range over |ℐ| ≤ x^{1/2} classes,
  with budget ≍ ε·Ψ(x). The deletion error is ≤ |𝒱|·x^{1/2}, which is ≤ ε·xρ iff |𝒱| ≤ ερx^{1/2}. At x ≥ X_r this holds
  iff 2ν < (1/d − ν)/2, i.e. **5ν < 1/d (d < 4k/5)**.
- *Regime.* This holds throughout the regime that matters: at the 2/5 design d ≍ D^{1/5} and ν ≍ D^{−2/5}, and for every
  fixed d once k > 5d/4. Algebraic values are larger still.

*Contrast with Theorem U [consistency check].* There the union over the field law of the values F_A(P) is an AP of
length ≍ HV^{d−1}. The point masses are ≤ 1/((2H + 1)V′^{d−1}) ([U] Lemma SM1). Deleting K values changes FM by
≤ K·#𝒫·max mass, which is negligible for any K = N^{O(ν)}. So the same statistical inputs do suffice. The barrier is
the small support of a fixed family, exactly Prop SI(i)/(ii) of [nest] but for the plain smooth set.

**What a proof would need instead.** It would need non-relativising, form-specific information: arithmetic of the
specific values (Norm(a − bα₀) in a Kummer field, a − b·r^j in base r), such as type I/II sums for incomplete norm forms
(Heath-Brown, Maynard), BBDT's methods for binary forms, or prescribed-digit methods (Bourgain). §4–§5 compare the
required ranges with what these give.

### 2.6 Device table (what survives)

"✗" = the ingredient's hypotheses fail. "✓" = it holds, but only by giving up the SNFS norms.

| device | (R) Theorem R | (GM) | (K) Kern | (W) | FM | LD* | keeps SNFS norms? |
|---|---|---|---|---|---|---|---|
| (a) random rational side, special root | ✗ (RS: only ±(x − m) below N^{1−1/d}) | — | — | — | ✗ | ✗ | yes, but then fixed g |
| (a′) random root m′ (BLP) | ✓ | needs a new f | — | — | GNFS | GNFS | **no** (base-m′ heights) |
| (a″) root twist λm | ✗ (DI) | ✗ | ✗ | ✗ | ✗ | ✗ | yes |
| (b) kN, N·w | ✗ | ✗ | ✗ | ✗ (nothing to spread; w > B kills smoothness) | ✗ | ✗ | yes |
| (c) f(x + t), u^d f(x/u), Möbius | ✗ (DI) | ✗ | ✗ | ✗ | ✗ | ✗ | yes (heights grow like t^d, u^d) |
| (d) the O(D²) special fields | ✗ (RB) | ✗ | ✗ | ✗ | ✗ | ✗ | yes |
| (e) random special-q only | ✗ (RB: union over ℓ ⊂ values on a disc) | ✗ | ✗ | ✗ | ✗ | ✗ | yes |
| (f) f₀ + (x − m)c, fixed g | ✗ | ✓ | ✓ | ✗ (Res fixed = r^jN) | ✗ | ✗ | **no** (PH: height ≥ (m−1)‖c‖) |
| (f′) f₀ + G·c with random g, w | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **no**: this is Theorem U's design; f̂_{g,w} no longer sees f₀ |

So, **proved**: every device that keeps ‖f‖ = N^{o(1/d)} keeps the field and the value set (HR, RS, DI). Its union of
values is N^{2ν+o(1)}-sparse, and relativising supply proofs fail (RB). Every device that restores a full-support field
law pays height ≥ N^{1/d} (HR, PH) and lands in Theorem U's design, which is GNFS-like.

---

## 3. Special families with natural rigour

### 3.1 Prop CY (cyclotomic pieces) [proved, conditional on Theorem U]

- Write r^e − 1 = ∏_{a|e}Φ_a(r) and r^e + 1 = ∏_{a|2e, a∤e}Φ_a(r).
- Each piece satisfies Φ_a(r) ≤ (r + 1)^{φ(a)}, since |r − ζ| ≤ r + 1, and φ(a) ≤ φ(e) for a | e.
- So every piece has log-size ≤ ε_e·log N·(1 + O(1/log N)), with ε_e := φ(e)log(r + 1)/(e log r). For r^e + 1 use
  φ(2e) instead of φ(e).

**Prop CY.** Completely factor each piece with Theorem U's algorithm, after computing the pieces (polynomially many
integer operations, `cyclo.algebraic_pieces`, which never factors N). This factors N = r^e ± 1 in expected time
L_N[2/5, c_U·ε_e^{2/5} + o(1)].

*Proof.*
- There are τ(e) ≤ e = O(log N) pieces.
- Theorem U gives, for every η > 0, a constant C_η with T_U(n) ≤ C_η + L_n[2/5, c_U + η].
- L_n[2/5, c] is increasing in n, and log Φ_a(r) ≤ ε_e t. Hence L_{Φ_a(r)}[2/5, c] ≤ exp(c(ε_et)^{2/5}(log t)^{3/5}) =
  L_N[2/5, cε_e^{2/5}].
- Pieces that share primes (those dividing e) only need a final merge.
- Correctness follows from Theorem U and AKS certification. ∎

**Corollary CY-P.** Along e = p_k# (primorials), ε_e ≤ (1 + o(1))log(r + 1)/(log r)·∏_{p≤p_k}(1 − 1/p) → 0 (Mertens),
so N_e = r^e − 1 factors in L_{N_e}[2/5, o(1)]. This is an infinite natural family with a strictly better rigorous
*constant*, but the same α. `cyclo.largest_piece_log_ratio` gives 0.31, 0.234 and 0.208 for e = 30, 210 and 2310
(r = 2).

**Aurifeuillian factors.** Where they exist (e.g. 2^{4k+2} + 1 = (2^{2k+1} − 2^{k+1} + 1)(2^{2k+1} + 2^{k+1} + 1);
in general Brent 1993 and Stevenhagen 1987), they split Φ_{2e}(r) into two halves of equal size. That multiplies the
constant by 2^{−2/5} for those (r, e) (`cyclo.special_pieces`, r = 2). Since φ(e) ≥ e/(e^γ log log e + 3/log log e),
the primitive piece always has log-size ≥ t/O(log log t). So **no exponent gain is possible from algebraic factors**:
the primitive part is a general-looking integer of size N^{Ω(1/log log log N)}.

### 3.2 Families defined by their factors (codex snfsrig_1 §4) [proved; recorded for completeness]

These are trivially fast, but they are not special *forms*.
- Semiprimes with q − p ≤ (log p)²: one Fermat step, polynomial time. There are infinitely many by the prime number
  theorem.
- Squarefree N with every p − 1 being L[β]-smooth: Pollard p − 1 with random bases gives expected L[β]. A Cunningham
  representation implies nothing of the kind. (For p | Φ_e(r) one only knows e | p − 1, which is log e bits.)

---

## 4. The precise obstruction (what would suffice)

Take the Eisenstein family 𝔈. The size model of §1.4 identifies two sufficient inputs for a rigorous L[1/3] in the
Dixon-split + rational-only-descent architecture.

**(O1) Hypothesis RT (rational-target supply in the Kummer fields) [open].** There are constants (b, b₁, s, λ, φ_F)
such that, for all large N = r^e − s ∈ 𝔈 and every prime q ∈ (e^{bE}, e^{b₁Q}] not dividing rsN, the following holds.
With d = d(log q) as in §1.4, the reduced-basis box of radius e^{sE} of {(a, b) : a ≡ bm_d (mod q)} contains a primitive
(a, b) such that
- F_d(a, b) = a^d − sr^{j_d}b^d is e^{bE}-smooth, and
- (a − bm_d)/q is q^λ-smooth (large primes coprime to rsN).

This is the rational-target-only, special-field case of [nest] Hypothesis AD. AD itself is not needed here (no
algebraic targets). In the random-integer model RT holds with room (§1.4: good-point margin 2s − (φ_F + φ_R)/3 > 0,
union bound over e^{O(Q)} targets).

**(O2) ECM-speed recognition [status: see §5].** A rigorous procedure that, for all p in some set G of primes of
density ≥ (log x)^{−C} in every dyadic range (x, 2x] with x ≤ L_N[2/3], finds every factor p ∈ G of its input in
expected time L_p[1/2, O(1)].

**Statement [model, not proved].** RT + (O2) + LV6 at degree O(D′) in each Kummer field, with multi-field row algebra
(one algebraic block per K_d, [nest] Prop R), give a Las Vegas L[1/3, c] algorithm for 𝔈. I have not written this
implication out. It is a sufficiency claim at the level of [nest] Theorem NR's proof. Theorem NR's single-field
argument and transcript invariance carry over; the multi-field LV6 square roots and character columns are UNVERIFIED
bookkeeping. Note that [nest] Theorem NR cannot simply be specialised to a point-mass field law: its recursive calls
need the algorithm for *all* integers. That is (O2) again.

**How far (O1) is from known results.**
- *Algebraic half.* RT needs F_d(a, b) = a^d − sr^jb^d, of degree d growing up to D′, to be y-smooth with
  log y ≍ (log|F|)^{1/2+o(1)}: u → ∞ like (log X)^{1/2}. And it needs this *pointwise* in each special-q box.
  BBDT 2012 give positive proportions for fixed degree-d binary forms only at y = x^{d−2+ε} (x^{1/√e+ε} for cubics),
  where x is the variable size: power-size y, bounded u. Barbulescu–Lachand 2017, Thm 4.2, reach
  y ≥ exp((log log X)^{5/3+ε}), but for fixed imaginary *quadratic* forms of fundamental discriminant, unrestricted
  boxes, and without joint smoothness with a rational form (codex). Nothing is known at growing degree or in special-q
  lattices.
- *Rational half.* The values a − b·r^{q_d} with |a|, |b| ≤ V are exactly the integers ≤ r^{q_d}V whose base-r
  digits in positions [log_r V, q_d) are all 0 or all r − 1 (those of a mod r^{q_d}). That is a proportion
  (1 − dν)/(1 + dν) → 1 of their digits. Prescribed-digit theorems for primes (Bourgain 2013/2015, Swaenepoel 2020)
  allow only a small positive proportion of prescribed digits. Even the rational half of RT, for smooth numbers instead
  of primes, is far beyond them. Also, by Prop RB, statistics of smooth numbers cannot replace it.

---

## 5. Literature

(To be completed from the literature subagent; see below.)

---

## 6. Reviews

(To be completed.)

---

## 7. Candidate card

- **Mechanism.** SNFS fixed polynomials f_d = x^d − sr^j (degree chosen per node) inside the project's
  Dixon-split + descent architecture. Randomness from w, kN, twists, root twists, several (d, j), special-q.
- **Model.** [rec]/[nest] size model; Theorem U's supply framework.
- **Result.**
  - Size model: one-level 2/5 (universal, Prop UF), multi-level 1/3, rational-only 1/3 (Prop R2-S; GNFS: 2/5).
  - Rigorous: the best α is 2/5 (Theorem U); constant ε_e^{2/5} for r^e ± 1 (Prop CY), o(1) along primorial e.
- **Failure scope.**
  - (O1) Fixed fields: HR/RS/PH/DI show that randomisation keeps the value set or costs height N^{1/d}.
    RB: relativising supply proofs fail for any fixed family of N^{o(ν)} fields.
  - (O2) Recognition (Prop REC): the split halves are not special. HECM or Theorem U recursion gives 2/5.
- **Most promising open question.** A pointwise friable-values estimate for the Kummer forms a^d − cb^d on special-q
  lattices at L-scale (Hypothesis RT), even only for the rational half a − b·r^j (prescribed-digit smooth numbers).
