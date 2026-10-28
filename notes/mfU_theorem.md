# Theorem U-MF: a Las Vegas factoring algorithm in expected time L_N[2/5]

Agent `mfU`, phase 30, 2026-10-28. This document is self-contained. It is meant to be read by an analytic number theorist
without the project's other notes. The review log, the audit of the earlier route and the history are in
`notes/mfU_results.md`.

**Status.** Every internal statement below has a complete proof here. External inputs are published, numbered theorems,
listed in §1.2 with exact citations. **Nothing has been refereed by a human.** All checks so far are AI reviews
(GPT-6.1-sol and Opus critics; `notes/mfU_results.md` §6). Read this as a careful preprint draft, not as an
established result.

---

## 0. Statement and overview

Notation: λ := log log N, D := log N/λ, Λ := D^{2/5}λ = (log N)^{2/5}(log log N)^{3/5}, and L_N[2/5, c] := e^{cΛ}.

> **Theorem U-MF.** Fix rationals a, κ, d₀ > 0. The algorithm FACTOR_MF = FACTOR_MF[a, κ, d₀] of §3 is Las Vegas.
> On every input N ≥ 2 it returns the complete prime factorisation of N (every prime is either a sieved prime ≤ B or
> AKS-certified). Its expected running time is **at most** L_N[2/5, c_MF + o(1)], where
>
>   c_MF = max( 2(a + σ_F),  a + σ_F + s + max(a, σ_r, σ_F),  (a + σ_F + s)/(1 − 2^{−2/5}) ) = (a + σ_F + s)/(1 − 2^{−2/5}),
>
>   σ_r = 2/(5a d₀),  σ_F = (2/(5a))(1/d₀ + d₀/(4κ)),  s = 2κ(1 + C_γ),
>
> and C_γ is a constant depending only on (a, κ, d₀) and on the implied constants of the cited theorems (§10.3). The
> o(1) is uniform in N.

The best rigorous bound for factoring arbitrary N that we found in the literature is Lenstra–Pomerance, J. Amer. Math.
Soc. 5 (1992), 483–516: L_N[1/2, 1 + o(1)]. This is a search result, not a novelty claim.

**Closest prior work: Lee–Venkatesan**, *Rigorous analysis of a randomised number field sieve*, J. Number Theory 187 (2018),
92–159 (arXiv:1805.08873). Their Thm 2.1 is an unconditional rigorous L[1/3, (64/9)^{1/3}] bound for finding a congruence of
squares, which may be trivial. Their Thm 2.3 gives fruitfulness only for semiprimes with both primes ≡ 3 (mod 4), and it is
conditional on their Conjecture 7.1 (a character-independence hypothesis). Theorem U-MF obtains fruitfulness for every N
without any such hypothesis (Dixon descent, Theorem DD), at the cost of exponent 2/5. We use two ingredients from LV:
1. their randomised polynomial f = f̂ + (a large random multiple of x − my), with uniform coefficients ((3.1), §4), which
   is the c-randomisation of our Lemmas SM1/OP;
2. the Stark–Landau analysis of exceptional square classes (§6), which our Lemma 8.2 extends.
As far as our literature critic checked, neither the integer multiplier w (resultant N·w) used to make the forced residue
uniform, nor norm-parity columns detecting the rational exceptional class, appears in prior work. Buhler–Lenstra–Pomerance
§5 calls a square norm "necessary … but very far from sufficient".

**Ideas.**
1. *Dixon descent* (§7). Draw a random square R = r² mod N. Reduce the lattice {(u, v) : u ≡ Rv (mod N)} to a shortest
   vector (u₀, v₀), and ask that |u₀| and |v₀| be products of at most k distinct primes from J = (Y, 2Y], with
   Y = N^{1/(2k)+o(1)}. Each such prime ℓ ("special-q") is then written as a product of B-smooth objects by one number
   field sieve relation. The hidden root r never enters any computation except the final gcd. So the final congruence of
   squares is nontrivial with probability ≥ 1/2 (Theorem DD), which gives *fruitfulness* without any character
   independence hypothesis.
2. *Relation supply* (§§5–6). A special-q prime ℓ must be "descendable" with probability γ ≥ (log N)^{−C}. The power
   matters: a split needs 2k ≍ D^{2/5} descendable primes simultaneously, so γ enters as γ^{2k}, and any loss of the form
   2^d or (log N)^d in γ would be fatal.
3. *Design MF* (one shared rational side, J_f independent algebraic sides). The rational side g = (g₀, g₁) is random.
   For each ℓ it picks the first candidate P₀ in the ℓ-lattice whose rational value is smooth. There are J_f ≈ 1/ρ_F
   independent algebraic polynomials F_j, all with the root m = g₀/g₁ mod N. Conditionally on the rational data the
   F_j(P₀) are independent, and each is smooth with probability ≥ ρ_F/61 (Lemma OP). So some j succeeds with probability
   ≥ 1/2. **No algebraic second moment is needed.** The rational supply needs only a first and a second moment over
   (g, ℓ) (Theorem R), which is elementary given Hildebrand–Tenenbaum.
4. *Integer multipliers* (Lemma OP). The resultant of F_j with the rational form is N·w_j, with w_j a uniform *integer*
   in (W, 2W] and W ≥ |G(P₀)|·N^{η₀/3}. Then the residue class of F_j(P₀) modulo G(P₀) is uniformly distributed, so the
   smoothness probability is an *interval* density. **No equidistribution of smooth numbers in progressions is needed.**
   In particular no Bombieri–Vinogradov-type theorem is used, and so no input whose proof rests on Harper's unpublished
   preprint arXiv:1208.5992. This is the main difference from the project's earlier Theorem U.
5. *Simultaneous certification* (§8). After quadratic characters, each field's residual square obstruction is either
   trivial or the class of one rational number. It is detected by the norm. So π(B) + 1 shared "norm-parity" columns
   annihilate all J_f residual obstructions at once, and one kernel vector suffices (Proposition PP).
6. *Recursion* (§9). The shape test factors |u₀|, |v₀| ≤ √(2N) by calling FACTOR_MF itself. Lemma F_α solves the
   recursion.

**Why 2/5 and not lower.** The descent is one level: norms of size N^{θ} with θ ≍ D^{−1/5}, k ≍ D^{2/5} special-q primes.
Balancing the split acceptance (log N)^{−2k} against the smoothness of N^{θ}-size values at B = L[α] gives α ≥ 2/5 for
this geometry. Remark 11.4 has the arithmetic. Nothing here goes below 2/5.

---

## 1. Notation and published inputs

### 1.1 Notation

- P⁺(n) is the largest prime factor (P⁺(1) = 1). 𝒮_B := {n ∈ Z ∖ {0} : P⁺(|n|) ≤ B}. "Smooth" means B-smooth.
- Ψ(x) := #{1 ≤ n ≤ x : P⁺(n) ≤ B}, ϱ(x) := Ψ(x)/x; Ψ(x, y) for other bounds y.
- ω(n): number of distinct prime factors. φ: Euler's function.
- μ_N := {ε ∈ (Z/N)^× : ε² = 1}; Q_N := ((Z/N)^×)². For odd N, μ_N ≅ {±1}^{ω(N)} and [(Z/N)^× : Q_N] = 2^{ω(N)}.
- A *gcd check* of an integer x computes gcd(x, N). If it is a proper divisor, the current procedure stops and returns it (a
  *factor exit*). If it is N, the current object is rejected.
- All parameters with real exponents (B, H, η₀, J_f, …) are computed by a fixed procedure that returns a value within a
  factor 2 of the stated real number, or an integer within ±1 of it (e.g. via rational approximations of log). Every
  statement below tolerates this slack; we suppress it.
- C, c denote positive constants depending only on (a, κ, d₀) and on the absolute constants of §1.2; they may change from
  line to line. They **never** depend on d, k or N. "For N ≥ N₀" means: for N larger than such a constant.
- z is **in range** if B^{D^{1/3}} ≤ z ≤ N². Then u_z := log z/log B satisfies D^{1/3} ≤ u_z ≤ D^{3/5}/a, and
  u_z ≤ log z/(2 log log z) for N ≥ N₀.

### 1.2 Published inputs (used only through these statements)

| Label | Statement used | Source (verified as stated in §6) |
|---|---|---|
| **HT1** | Ψ(x, y) = x^{α}ζ(α, y)/(α√(2πφ₂(α, y)))·(1 + O(1/u + log y/y)) uniformly for x ≥ y ≥ 2, at the saddle point α = α(x, y) (root of Σ_{p≤y} log p/(p^α − 1) = log x), with φ₂(α, y) = Σ_{p≤y}p^α(log p)²/(p^α − 1)². | A. Hildebrand, G. Tenenbaum, *On integers free of large prime factors*, Trans. Amer. Math. Soc. 296 (1986), 265–290, Theorem 1 |
| **HT2** | α(x, y) = log(1 + y/log x)/log y·(1 + O(log log(1 + y)/log y)) uniformly for x ≥ y ≥ 2. | ibid., Theorem 2 (part "(i)" UNVERIFIED against print) |
| **HT3** | Ψ(cx, y) = c^{α(x,y)}Ψ(x, y)(1 + O(1/u + log y/y)) uniformly for x ≥ y ≥ 2, 1 ≤ c ≤ y. | ibid., Theorem 3 |
| **CEP** | For fixed ε > 0, Ψ(x, x^{1/u}) = x·u^{−u(1+o(1))} uniformly for u → ∞, u ≤ (1 − ε)log x/log log x. | E. R. Canfield, P. Erdős, C. Pomerance, *On a problem of Oppenheim concerning "factorisatio numerorum"*, J. Number Theory 17 (1983), 1–28, the unnumbered Corollary after Theorem 3.1 (§3): for ε > 0, 3 < u < (1 − ε)log x/log₂x, Ψ(x, x^{1/u}) = x·exp{−u(log u + log₂u − 1 + (log₂u − 1)/log u + E)}, |E| ≤ c_ε(log₂u/log u)² |
| **PNT** | π(x) > x/log x (x ≥ 17); π(x) < 1.25506x/log x (x > 1); π(x) − π(x/2) > 3x/(10 log x) (x ≥ 41); and π(x) − π(x/2) ≤ 0.6x/log x for x ≥ x₀. | J. B. Rosser, L. Schoenfeld, *Approximate formulas for some functions of prime numbers*, Illinois J. Math. 6 (1962), 64–94, Theorem 1 and Corollaries 1–3 |
| **Stark** | If a number field L of degree n has ζ_L(β₀) = 0 with β₀ real, 1 − (c(L) log|Δ_L|)^{−1} ≤ β₀ < 1, where c(L) = 4·n! in general, then ζ_F(β₀) = 0 for some quadratic subfield F ⊂ L. | H. M. Stark, *Some effective cases of the Brauer–Siegel theorem*, Invent. Math. 23 (1974), 135–152, Lemma 8 |
| **Landau** | There is an effective absolute c_L > 0 such that for distinct real primitive characters χ₁, χ₂ of conductors q₁, q₂, L(s, χ₁)L(s, χ₂) has at most one real zero, counted with multiplicity, in (1 − c_L/log(q₁q₂), 1) (so χ₁, χ₂ cannot both vanish at one β there). | H. Davenport, *Multiplicative Number Theory*, 3rd ed., GTM 74, Ch. 14; stated with multiplicity in H. L. Montgomery, R. C. Vaughan, *Multiplicative Number Theory I*, CUP 2007, Ch. 11 |
| **LO** | For L/K Galois with group G, conjugacy class 𝒞, n_L = [L : Q], d_L = |Δ_L| and x ≥ exp(10n_L(log d_L)²): |π_𝒞(x) − (|𝒞|/|G|)Li(x)| ≤ (|𝒞|/|G|)Li(x^{β₀}) + c₁x·exp(−c₂n_L^{−1/2}(log x)^{1/2}), c₁, c₂ absolute and effective; the β₀-term only if ζ_L has a real zero β₀ ∈ [1 − (4 log d_L)^{−1}, 1). | J. C. Lagarias, A. M. Odlyzko, *Effective versions of the Chebotarev density theorem*, in: Algebraic Number Fields (Durham 1975), Academic Press 1977, 409–464, Theorem 1.3 (number UNVERIFIED against print; Pierce–Turnage-Butterbaugh–Wood call it Corollary 1.3; content confirmed) |
| **AKS** | Primality is decidable deterministically in time (log n)^{O(1)}. | M. Agrawal, N. Kayal, N. Saxena, Ann. of Math. 160 (2004), 781–793 |
| **CZ** | Factoring a degree-d polynomial over F_r: Las Vegas, expected (d log r)^{O(1)}. | D. G. Cantor, H. Zassenhaus, Math. Comp. 36 (1981), 587–592 |
| **KS** | A uniformly random element of the kernel of a sparse M × M matrix over F_q (q ≥ M⁴) is found Las Vegas in expected Õ(M·nnz + M²) field operations. | E. Kaltofen, B. D. Saunders, *On Wiedemann's method of solving sparse linear systems*, AAECC-9, LNCS 539 (1991), 29–38, Theorems 2–4. Their rank computation (Thm 3) is Monte Carlo; Las Vegas behaviour here comes from the exact verification in §8.3, and uniformity on the F₂-kernel from our own coset/projection argument (§8.3 step 2), not from KS Thm 4 |
| **Gauss** | Gauss's lemma (products of primitive integral polynomials are primitive); Dedekind–Kummer; Kummer theory of quadratic extensions; Minkowski's bound; Gauss reduction of rank-2 lattices. | standard |

Standard algebraic and elementary inputs (textbook; used in Lemmas FIB, 8.1, 8.4, 8.5): Dirichlet's unit theorem
(|𝒪_K^×/𝒪_K^{×2}| = 2^{r₁+r₂}), the trace-dual identity Z[β]^∨ = g′(β)^{−1}Z[β] for a monogenic order (Euler), φ(t)/t ≫ 1/log log t,
ω(t) ≪ log t/log log t, and the Mahler-measure and discriminant bounds |disc f| ≤ d^{2d−1}‖f‖_∞^{2d−2}. Lower-order inputs, all standard: Stanley's log-concavity of convolutions (R. P. Stanley, Ann. NY Acad. Sci. 576
(1989), 500–535, Prop. 2 and §2; proposition number UNVERIFIED), the remainder-tree and product-tree algorithms (J. von zur Gathen, J. Gerhard,
*Modern Computer Algebra*, §§10.1–10.3), and Tonelli–Shanks.

**Not used** (in contrast with the project's Theorem U): Drappeau's Lemme 5 and therefore Harper's arXiv:1208.5992; any
Bombieri–Vinogradov or large-sieve estimate; Brun–Titchmarsh; Tenenbaum's Ψ_q theorem; Lenstra–Pomerance (the shape test
is recursive); the hyperelliptic curve method.

---

## 2. Parameters

The input of SPLIT is an odd N with no prime factor ≤ B that is not a prime power. Fixed rationals a, κ, d₀ > 0.

- d := the least odd integer ≥ max(3, d₀D^{1/5}); k := max(⌈κD^{2/5}⌉, d + 1); ν := 1/(4k); μ := 1/d.
- η₀ := 1/⌈D^{1/2}⌉.
- B := ⌊e^{aΛ}⌋.
- X := ⌊√(N/2)⌋; j₀ maximal with 2^{kj₀} ≤ X; Y := 2^{j₀−1}; J := {ℓ prime : Y < ℓ ≤ 2Y}; π_J := #J.
- σ̄ := (2/(5a))·(2/d₀ + d₀/(4κ)) (a rational constant; σ̄ = σ_r + σ_F).
- M := ⌊N^{μ}⌋; W := ⌈N^{μ+ν+η₀/2}⌉; **H := ⌊N^{ν+η₀}⌋·⌈exp(30(σ̄ + 1)·d·Λ·λ)⌉** (the second factor is the *MF height floor*).
- S: the first S = 2^i, i = 1, 2, …, with S²·R(X_r(S)) ≥ 1, where V(S) := ⌊√(2Y)·S⌋, X_r(S) := 2M·V(S)/Y and R is
  the computable Rankin value of Lemma CS. Then **V := V(S)**, X_r := X_r(S), **X_F := (2d + 1)MHV^d**,
  ρ_r := ϱ(X_r), ρ_F := ϱ(X_F).
- **J_f := ⌈(log N)^{λ₃}/R(X_F)⌉** with λ₃ := max(1, log log log N), the number of algebraic sides.
- Candidates: 𝒫 := {P = (a, b) ∈ Z² : b prime, V/2 < b ≤ V, |a| ≤ V, gcd(a, b) = 1}, in a fixed order.

Sizes (leading order; D^{1/5}θ → θ₀ below): with θ := 2μ + h + (d − 1)ν, h := log H/log N,
- log X_r = (μ − ν)log N + O(Λ) = (1/d₀)D^{4/5}λ(1 + o(1));
- log X_F = (μ + h + dν)log N + d log(√2S) + O(log d) = (1/d₀ + d₀/(4κ))D^{4/5}λ(1 + o(1));
- the MF height floor contributes 30(σ̄ + 1)dΛλ ≍ D^{3/5}λ² to log H, and d log S ≍ D^{3/5}λ to log X_F; both are o(D^{4/5}λ).

**Lemma P (margins) [proved; leading-order comparisons, encoded in `src/mfU/params.py::margins`].** For N ≥ N₀:
- (P1) B^{C} ≤ Y, 32S² < Y (so b > V/2 > 4S²), B < V/2, V < ℓ for ℓ ∈ J; and V² ≥ 7Y.
- (P2) N^{η₀/4} ≥ max(B^{C}, S^{C}, 2^{Cd}); explicitly used: 2^{d+12}dB ≤ B², 2^{d+9}B ≤ N^{η₀/4}, 2^{d+5} ≤ N^{η₀/4}, N^{η₀/2} ≥ 10⁵d²/ρ_F.
- (P3) 2MV·N^{η₀/3} ≤ W and 2W < N; M/V ≥ N^{η₀/2}; H/V ≥ N^{η₀/2}; M(H + 1) < N; 2H + 1 < N; 2X² < N.
- (P4) d log d ≤ Λ/C; d log S + log J_f = O(Λ·D^{1/5}) = o(D^{4/5}λ) (lower order than log X_F); (log N)^{λ₃} = e^{o(Λ)}.
- (P5) log Y ≥ ω(N) log 2 + C log log N (using ω(N) ≤ log N/log B).
- (P6) log H/(24 d log log H) ≥ 1.25(σ̄ + 1)Λ, and 2^{(log H)/2} ≤ W·N^{−1/(2d)}.
- (P7) X_r^{1/2} and X_F are in range, and u_{X_F} ≤ log X_F/(2 log log X_F).
- (P8) log(1/ρ_r) = σ_rΛ(1 + o(1)), log(1/ρ_F) = σ_FΛ(1 + o(1)), with σ_r, σ_F as in Theorem U-MF.

*Proof.* Each is a comparison of the leading terms of two expressions of the form c(log N)^p(log log N)^q. (P1): log Y ≍
D^{3/5}λ against log B, log S ≍ Λ = D^{2/5}λ; V² ≥ (2√(2Y) − 1)² ≥ 7Y as S ≥ 2. (P2): η₀ log N ≍ D^{1/2}λ ≫ Λ ≫ d. (P3):
log W − log(2MV) = (η₀/2)log N − log S − O(1); h − ν = η₀ + (MF floor) ≥ η₀; μ + h < 1 and ν + η₀ < 1/2. (P4):
d log d ≍ D^{1/5}λ; d log S ≍ D^{3/5}λ. (P5): log Y ≍ D^{3/5}λ/(2κ) against log N/log B ≍ D^{3/5}/a, a factor λ apart.
(P6): log H ≥ 30(σ̄ + 1)dΛλ and log log H ≤ λ; and (log H)/2 · log 2 ≍ D^{3/5}λ² ≪ log N/(2d) ≍ D^{4/5}λ. (P7):
u ≍ D^{2/5}. (P8): u_{X_r} = (1/(ad₀))D^{2/5}(1 + o(1)), log u = (2/5)λ(1 + o(1)), and CEP; similarly for X_F. ∎

*Where the asymptotics start.* The N^{η₀}-margins beat L[2/5]-scale quantities only by a factor D^{1/10} in the exponent,
so N₀ is astronomically large (for N^{η₀/4} ≥ B^{10}, D ≳ 10^{14}). **The algorithm never uses N₀**: FACTOR_MF runs SPLIT
and trial division in alternation (§3), so the theorem is a statement about the asymptotic expected time only.

---

## 3. The algorithm

**FACTOR_MF(N)** (returns the sorted list of primes with multiplicity).
1. If N = 1 return []. If N is prime (AKS), return [N]. If N = m^e with e ≥ 2 (integer roots, e ≤ log₂N), return
   FACTOR_MF(m) with multiplicities times e.
2. Trial-divide by all primes ≤ B(N). Let N′ be the cofactor. If N′ = 1 return; if N′ is prime or a perfect power,
   treat it as in step 1.
3. Run SPLIT(N′) and trial division of N′ **in alternation, one bit operation each** (a bit operation of SPLIT's whole nested computation, recursive shape tests included), until one returns a proper factor g.
4. Return the sorted merge of FACTOR_MF(g), FACTOR_MF(N′/g) and the primes of step 2.

**SPLIT(N)** (N odd, all prime factors > B, not a prime power).
0. Compute the parameters of §2 (including the scans for S and J_f). If a structural prerequisite fails (J = ∅,
   M < 4, H < 1, 2X² ≥ N, or S² ≥ Y), return FAIL forever (this happens only for N < N₀; trial division then finishes).
1. *Field selection.* Repeat: A ← SAMPLE_FIELD_MF; row ← DESCEND_ONCE(A); if a row is returned, keep A and go to 2.
   (Any factor exit returns the factor.)
2. *Rows.* Compute M_rows(A) (§8.2). Call DESCEND_ONCE(A) until M_rows rows are held.
3. *Post-processing.* (𝒮, 𝖸) ← POSTPROCESS_MF(A, rows). Compute 𝖷 := ∏_{i∈𝒮} r_i mod N and g := gcd(𝖷 − 𝖸, N).
   If 1 < g < N return g.
4. Otherwise discard all rows and go to step 2 (same A, fresh hidden roots).

**SAMPLE_FIELD_MF.** Returns A = (g; (w_j, c_j, f̃_j, acc_j)_{j ≤ J_f}).
1. Draw (g₀, g₁) uniform in [M/2, M]² ∩ Z² until gcd(g₀, g₁) = 1. Gcd-check g₁. Put G(x, y) := g₁x − g₀y and
   m := g₀g₁^{−1} mod N.
2. For j = 1, …, J_f independently:
   a. Draw w_j **uniform among the integers** in (W, 2W] and c_j uniform in C_H := [−H, H]^d ∩ Z^d. Gcd-check w_j.
   b. F_j := F̂_{g,w_j} + G·R_{c_j} (Definition 4.1), f_j(x) := F_j(x, 1). If the leading coefficient f_{j,d} = 0, set
      acc_j := false and continue. Gcd-check f_{j,d}.
   c. Let e_j be the content of f_j; gcd-check e_j; put **f̃_j := f_j/e_j** and F̃_j := F_j/e_j.
   d. *Irreducibility certificate:* if no prime r ∈ (16, ½ log H] has f̃_j mod r of degree d and irreducible over F_r,
      set acc_j := false and continue.
   e. Gcd-check f̃_{j,d}^{d−2}f̃_j′(m); if N divides it, acc_j := false and continue.
   f. acc_j := true.
3. Return A.

**DESCEND_ONCE(A).**
1. Draw r uniform in [1, N − 1]; gcd-check r. Put R := r² mod N.
2. Let (u₀, v₀) be a shortest nonzero vector of Λ_R := {(u, v) ∈ Z² : u ≡ Rv (mod N)} (Gauss reduction, fixed
   tie rule). Gcd-check u₀ and v₀.
3. *Shape test:* run FACTOR_MF on |u₀| and on |v₀| with fresh coins. Require that each is a product of at most k
   **distinct** primes, all in J. Otherwise return FAIL.
4. For each prime ℓ | u₀v₀ (in increasing order): gcd-check ℓ and call SQ_SEARCH_MF(A, ℓ) (cached per (A, ℓ)). If one
   fails, return FAIL.
5. Return the row (r, R, sgn(u₀v₀), (e_ℓ)_ℓ, (P_ℓ, j_ℓ)_ℓ), with e_ℓ = +1 if ℓ | u₀ and −1 if ℓ | v₀, and the
   factorisations of n_{P_ℓ} := |G(P_ℓ)|/ℓ and F̃_{j_ℓ}(P_ℓ).

**SQ_SEARCH_MF(A, ℓ)** (a deterministic function of (A, ℓ)).
1. Gauss-reduce Λ_ℓ(g) := {P ∈ Z² : ℓ | G(P)}. If its first minimum λ₁ < ℓ/(2V), return FAIL.
2. Enumerate Λ_ℓ(g) ∩ [−V, V]² (≤ 400S² points, Lemma 4.4). Keep the P ∈ 𝒫 (b prime by AKS) with |G(P)| ≥ MV/32.
   Compute the B-smooth parts of all |G(P)|/ℓ in one batch (Lemma 9.1). Let P₀ be the **first** kept P with
   G(P)/ℓ ∈ 𝒮_B. If there is none, return FAIL. Gcd-check G(P₀).
3. For every j with acc_j, compute F̃_j(P₀) and gcd-check it; compute the smooth parts of all of them in one batch.
   Return (P₀, j) for the least j with F̃_j(P₀) ∈ 𝒮_B, with the factorisations of |G(P₀)|/ℓ and F̃_j(P₀) (trial
   division by primes ≤ B). If there is none, return FAIL.

**POSTPROCESS_MF(A, rows).** §8.3.

*Coins and transcripts.* Every computation except the final 𝖷 in SPLIT step 3 is a function of the squares R_i, the
field A, the shape-test outputs (canonical factorisations) and coins independent of the roots r_i. This is checked in
Lemma 7.3.

---

## 4. Elementary lemmas: field law, sizes, regularity

**Definition 4.1 (base expansion; cf. Buhler–Lenstra–Pomerance 1993 §3; the random multiple G·R_c with uniform c is Lee–Venkatesan's randomised polynomial, JNT 187 (2018), (3.1) and §4).** Let g₀, g₁ ∈ [M/2, M] with gcd(g₀, g₁) = 1, and
let w ≥ 1 be an integer. Put T_d := Nw and, for i = d, d − 1, …, 1: a_i is the integer in (−g₁/2, g₁/2] with
a_i ≡ T_ig₀^{−i} (mod g₁), and T_{i−1} := (T_i − a_ig₀^i)/g₁. Put a₀ := T₀ and F̂_{g,w}(x, y) := Σ_{i=0}^{d} a_ix^iy^{d−i}.
For c ∈ Z^d put R_c(x, y) := Σ_{i<d}c_ix^{d−1−i}y^i and F := F̂_{g,w} + G·R_c.

**Lemma 4.2 (base expansion) [proved].** F̂_{g,w}(g₀, g₁) = Nw, hence F(g₀, g₁) = Nw; |a_i| ≤ g₁/2 for i ≥ 1 and
|a₀| ≤ 2^dNw/M^d + d2^{d−1}M. If gcd(g₁, N) = 1, then f(m) ≡ 0 (mod N) for m = g₀g₁^{−1} and f(x) := F(x, 1).
*Proof.* T_i ≡ a_ig₀^i (mod g₁), so T_{i−1} ∈ Z and inductively T_i = Σ_{l≤i}a_lg₀^lg₁^{i−l}; at i = d this is F̂(g₀, g₁).
F − F̂ = G·R_c vanishes at (g₀, g₁). For a₀: g₁^da₀ = Nw − Σ_{i≥1}a_ig₀^ig₁^{d−i} and |a_i|g₀^ig₁^{d−i} ≤ (M/2)·M^{d−1}·2^{i}…
precisely |a₀| ≤ Nw/g₁^d + Σ_{i≥1}(g₁/2)(g₀/g₁)^i ≤ 2^dNw/M^d + d2^{d−1}M, since g₀/g₁ ≤ 2. Finally
f(m) ≡ g₁^{−d}F(g₀, g₁) = g₁^{−d}Nw ≡ 0. ∎

**Lemma 4.3 (sizes) [proved].** For P = (a, b) with |a|, |b| ≤ V and every w ∈ (W, 2W]:
|G(P)| ≤ 2MV, |F̂_{g,w}(P)| ≤ MHV^dN^{−η₀/4}, and |F(P)| ≤ X_F for every c ∈ C_H. Also ‖f‖_∞ ≤ 3M(H + 1) < N.
*Proof.* |G(P)| ≤ |g₁a| + |g₀b| ≤ 2MV. Since M^d ≥ N(1 − dN^{−1/d}) ≥ N/2 (as M = ⌊N^{1/d}⌋), |a₀| ≤ 2^{d+1}·2W + d2^{d−1}M
≤ 2^{d+3}N^{μ+ν+η₀/2}. Then |F̂(P)| ≤ (dM/2 + |a₀|)V^d ≤ 2^{d+4}N^{μ+ν+η₀/2}V^d ≤ MHV^dN^{−η₀/4}, because
MH ≥ N^{μ+ν+η₀}/2 and 2^{d+5} ≤ N^{η₀/4} (P2). Next |G(P)R_c(P)| ≤ 2MV·dHV^{d−1}, so |F(P)| ≤ (2d + 1)MHV^d = X_F.
The coefficients of G·R_c are ≤ 2MH, and |a₀| ≤ 2^{d+3}N^{μ+ν+η₀/2} ≤ MH. ∎

**Lemma 4.4 (regularity, with one box for all ℓ) [proved; cf. Pollard 1993].** Let ℓ ∈ J and gcd(g₀, g₁) = 1. If
λ₁(Λ_ℓ(g)) < ℓ/(2V), then Λ_ℓ(g) ∩ 𝒫 = ∅. Otherwise #(Λ_ℓ(g) ∩ [−V, V]²) ≤ 400S².
*Proof.* Λ_ℓ(g) has index ℓ (g ≢ 0 mod ℓ, since gcd(g₀, g₁) = 1). A shortest vector w₀ is primitive in the lattice (its content is
< ℓ, hence a unit mod ℓ, and w₀/content would be shorter). The lattice is the union of lines parallel to w₀ at spacing
ℓ/|w₀| > 2V; every point of [−V, V]² has norm ≤ √2V < 2V, so the box meets only the line Zw₀, whose primitive points are ±w₀,
and |w₀| < ℓ/(2V) ≤ Y/V ≤ V/7 < V/2 < b for (a, b) ∈ 𝒫 (P1). In the regular case take a reduced basis λ₁ ≤ λ₂,
λ₁λ₂ ≤ (2/√3)ℓ; fundamental parallelograms at the box points are disjoint and lie in a square of side 2V + 4λ₂, so the count
is ≤ (8V² + 32λ₂²)/ℓ with λ₂² ≤ (4/3)ℓ²/λ₁² ≤ (16/3)V²; this is ≤ 179V²/ℓ ≤ 358S² since V² ≤ 2YS² and ℓ > Y. ∎

**Lemma 4.5 (homogeneity and the forced residue) [proved].** For every P, g₁^dF(P) ≡ Nw·b^d (mod G(P)). If moreover
P = (a, b) ∈ 𝒫 ∩ Λ_ℓ(g), ℓ ∤ N, and n := G(P)/ℓ ∈ 𝒮_B (N having no prime factor ≤ B), then t := |G(P)| satisfies
gcd(t, Ng₁b) = 1 and

  F̂_{g,w}(P) ≡ c_P·w (mod t),  c_P := Nb^dg₁^{−d} ∈ (Z/t)^×, for every integer w,

where c_P depends on (ℓ, P, g) but not on w.
*Proof.* g₁P − b(g₀, g₁) = (G(P), 0) and F has integer coefficients, so g₁^dF(P) = F(g₁P) ≡ F(bg₀, bg₁) = b^dNw. The same
holds for F̂ (take c = 0). b is a prime with B < V/2 < b ≤ V < ℓ (P1), so b ∤ ℓ|n| = t. A prime p | gcd(g₁, t) divides
g₀b, hence b (gcd(g₀, g₁) = 1), impossible. gcd(t, N) = 1 because n is B-smooth and ℓ ∤ N. ∎

**Lemma 4.6 (counting) [proved].** #𝒫 ≥ V²/(4 log V) and #𝒫 ≤ 3V². (PNT: > 3V/(10 log V) primes b ∈ (V/2, V]; for each,
≥ 2V − 3 values of a ∈ [−V, V] are coprime to b.) Also π_J ≥ Y/(2 log Y), and at most log N/log Y ≤ 2k + 1 primes of J divide N. ∎

---

## 5. Smooth numbers

**Lemma U1 (local behaviour of ϱ) [proved].** Fix C ≥ 1. For N ≥ N₀ and z₁ ≤ z₂ ≤ z₁B^C with z₁ in range:
(i) ϱ(z₂) ≤ 2ϱ(z₁); (ii) ϱ(z₁) ≤ (log N)^{C_U(C+1)}ϱ(z₂), C_U absolute; (iii) for fixed Z ≥ 1,
Ψ(z₁(1 + 1/Z)) − Ψ(z₁) ≥ Ψ(z₁)/(3Z); (iv) log(1/ϱ(z)) = u_z log u_z(1 + o(1)).
*Proof.* By HT2, (1 − α(z, B))log B = log log z − log(1 + log z/B) + O(log log B) ≤ λ + C log log B, so
B^{1−α} ≤ (log N)^{C_U}; and α < 1 (the defining sum at α = 1 is O(log B) < log z). For 1 ≤ c ≤ B, HT3 gives
ϱ(cz)/ϱ(z) = c^{α−1}(1 + O(1/u_z)), which lies in [B^{α−1}(1 − O(1/u)), 1 + O(1/u)]. Write z₂/z₁ as a product of at most
C + 1 factors in [1, B] (and apply HT3 at each intermediate point, which is in range); this gives (i) and (ii) because
1/u ≤ D^{−1/3} → 0 uniformly. (iii): Ψ(z₁(1 + 1/Z)) − Ψ(z₁) = Ψ(z₁)((1 + 1/Z)^α(1 + O(1/u)) − 1) ≥ Ψ(z₁)(α/(Z + 1) −
O(1/u)) ≥ Ψ(z₁)/(3Z) as α → 1. (iv): CEP, applicable by (P7). ∎

**Lemma CS (computable box) [proved].** For z in range let α̂(z) be the root of Σ_{p≤B}log p/(p^α − 1) = log z, computed by
bisection to absolute accuracy (log N)^{−3}, and R(z) := z^{α̂−1}∏_{p≤B}(1 − p^{−α̂})^{−1}, computed to relative accuracy 2. Then
R*(z)/2 ≤ R(z) ≤ 2R*(z), where R* is the exact value at α̂, and ϱ(z) ≤ R*(z) ≤ (log N)^{C_R}ϱ(z) (C_R absolute).
The scan of §2 stops at some S with S² ≤ e^{(σ_r+o(1))Λ}, and at that S

  m_r := S²ρ_r ∈ [(log N)^{−C_R}/4, 16].

Moreover R(X_F)^{−1} ∈ [ρ_F^{−1}(log N)^{−C_R}/2, 2ρ_F^{−1}], so J_f = ρ_F^{−1}·e^{o(Λ)} = L^{σ_F+o(1)} and
**J_f ≥ 244/ρ_F** for N ≥ N₀. Computing S and J_f costs B^{1+o(1)}.
*Proof.* *Rankin:* for α ∈ (0, 1], Ψ(z, B) ≤ Σ_{P⁺(n)≤B}(z/n)^α = z^α∏_{p≤B}(1 − p^{−α})^{−1}. *Lower bound:* by HT1, at the saddle
point α the Rankin value exceeds Ψ(z, B) by the factor α√(2πφ₂)(1 + o(1)) ≤ log N (φ₂ ≍ log z log B). The function
g(α′) := log(z^{α′}ζ(α′, B)) is convex with g′(α) = 0 and g″ = φ₂; for |α̂ − α| ≤ (log N)^{−3} every term of φ₂ changes by a
factor ≤ e³, so g(α̂) − g(α) ≤ e³φ₂(α)(log N)^{−6}/2 ≤ 1. Bisection is valid because the defining sum decreases in α.
*Existence:* for 4 ≤ S² ≤ e^{2σ_rΛ}, log X_r(S) = (μ − ν)log N + O(Λ), so by U1(iv) and the R*/ϱ bounds,
log R(X_r(S)) = −σ_rΛ(1 + o(1)) uniformly; at S² = e^{(σ_r+ε)Λ} the product S²R ≥ e^{(ε+o(1))Λ} ≥ 1. *m_r:* at the
chosen S, 1 ≤ S²R(X_r) ≤ 2S²R*(X_r) ≤ 2(log N)^{C_R}S²ρ_r, so m_r ≥ (log N)^{−C_R}/2. If S ≥ 4, the scan did not stop at
S/2, so S²/4 < 1/R(X_r(S/2)) ≤ 2/R*(X_r(S/2)) ≤ 2/ϱ(X_r(S/2)) ≤ 4/ρ_r by U1(i) (X_r(S)/3 ≤ X_r(S/2) ≤ X_r(S)); so
m_r < 16. If S = 2, m_r = 4ρ_r < 16. *J_f:* R(X_F) ≤ 2R*(X_F) ≤ 2(log N)^{C_R}ρ_F gives J_f ≥ (log N)^{λ₃−C_R}/(2ρ_F) ≥ 244/ρ_F
once (log N)^{λ₃−C_R} ≥ 488; and R(X_F) ≥ R*/2 ≥ ρ_F/2 gives J_f ≤ 2(log N)^{λ₃}/ρ_F + 1. *Cost:* each R costs
π(B)(log N)^{O(1)}; the scan makes O(Λ) evaluations. ∎

**Lemma DS (smooth numbers in short progressions) [proved].** Fix C₀ ≥ 1. Let X ≥ B^{D^{1/3}+2C₀+1} be in range (true for X = X_r, X_F with a large margin), 1 ≤ q ≤ B^{C₀}, and
I ⊂ [XB^{−C₀}, X] an interval of length A ≥ XB^{−C₀} (so A/q ≥ XB^{−2C₀} ≥ 1, X being in range). Then for every r,
#{n ∈ I : n ≡ r (mod q), P⁺(n) ≤ B} ≤ (A/q)ϱ(X)(log N)^{C_U(2C₀+2)+2}.
*Proof.* Let e := gcd(r, q). If e has a prime factor > B the count is 0. Otherwise n = en₁ with n₁ ∈ I/e,
n₁ ≡ r/e (mod q′), q′ := q/e, gcd(n₁, q′) = 1. Put Δ₀ := min(A/q, X/(eB^{C₀})) ≤ min(I/e). Every smooth n₁ ≥ Δ₀ has a divisor
d₁ ∈ (Δ₀/B, Δ₀] (the largest product of its smallest prime factors not exceeding Δ₀). Writing n₁ = d₁m: (d₁, q′) = 1,
m ≡ (r/e)d₁^{−1} (mod q′) and m lies in an interval of length A/(ed₁), so there are ≤ A/(d₁q) + 1 ≤ 2A/(d₁q) choices
(d₁ ≤ A/q). Hence the count is ≤ (2A/q)Σ_{Δ₀/B<d₁≤Δ₀, P⁺(d₁)≤B}1/d₁ ≤ (2A/q)·8 log₂(2B)·ϱ(Δ₀/B) (dyadic blocks and U1(i)).
Since Δ₀ ≥ XB^{−2C₀}, X/(Δ₀/B) ≤ B^{2C₀+1}, and U1(ii) gives ϱ(Δ₀/B) ≤ (log N)^{C_U(2C₀+2)}ϱ(X). ∎

---

## 6. Relation supply

### 6.1 The rational side

For ℓ ∈ J and g, put Y_P(g) := 1 if P ∈ 𝒫 ∩ Λ_ℓ(g), |G_g(P)| ≥ MV/32 and G_g(P)/ℓ ∈ 𝒮_B; else 0. Put
N_r(g, ℓ) := Σ_P Y_P(g) if ℓ ∤ N, and N_r := 0 if ℓ | N. Let g be uniform on 𝒢 := {(g₀, g₁) ∈ [M/2, M]² ∩ Z² : gcd = 1}
(|𝒢| ≍ M²) and ℓ uniform on J, independent.

**Lemma FIB (rational fibres) [proved].** Fix ℓ ∈ J and P = (a, b) ∈ 𝒫, and put 𝒢_P(n) := {g ∈ 𝒢 : G_g(P) = ℓn}.
(i) |𝒢_P(n)| ≤ 2M/V for every n. (ii) There is a set I_P of consecutive integers of length ≥ MV/(16ℓ) − 1 on which ℓn has
constant sign, MV/32 ≤ ℓ|n| ≤ 2MV and |𝒢_P(n)| ≥ M/(C_FV log log N).
*Proof.* (i) G_g(P) = g₁a − g₀b is constant along the primitive direction (a, b) of length ≥ b > V/2, and the square has
diameter M/√2. (ii) For (g₀, g₁) uniform on the real square [M/2, M]², G = ag₁ − bg₀ is a sum of independent uniforms on
intervals of lengths α′ ≤ β′ = V′M/2 (V′ := max(|a|, b) ∈ (V/2, V]); its density is ≥ 1/(2β′) on an interval of length
β′, and the density at s equals (chord length of {G = s})/(|(a, b)|M²/4), so there chords have length ≥ M/4. Removing
(−MV/32, MV/32) leaves a subinterval of length ≥ MV/16 of constant sign; I_P := the integers n with ℓn in it. For such n the
fibre line meets the closed square in J₀ ≥ ⌊(M/4)/|(a, b)|⌋ ≥ M/(8V) consecutive lattice points g* + i(a, b).
*Coprimality:* choose x, y with xa + yb = 1; the unimodular map g ↦ (xg₀ + yg₁, G) shows that e | gcd(g₀, g₁) iff e | ℓn and
i ≡ i₀(e) (mod e). By Möbius over e | ℓ|n|, at least J₀φ(ℓ|n|)/(ℓ|n|) − 2^{ω(ℓ|n|)} of the J₀ points are primitive. Finally
φ(t)/t ≥ c/log log t, and log 2^{ω(t)} ≤ C log t/log log t = O(μ log N/λ) = o(log(M/V)) since log(M/V) ≥ (μ/2)log N. ∎

**Proposition R1 (rational first moment) [proved].** E_{g,ℓ}N_r ≥ c·m_r/(log V·log log N).
*Proof.* For ℓ ∤ N, E_gN_r ≥ |𝒢|^{−1}Σ_{P∈𝒫}Σ_{n∈I_P∩𝒮_B}|𝒢_P(n)| (the filters hold on I_P). After reflection
I_P = (z₁, z₂] with z₂ − z₁ =: λ_P ≥ MV/(16ℓ) − 1 and z₂ ≤ 2MV/ℓ ≤ 33λ_P. If z₁ ≤ z₂/2, U1(iii) (Z = 1) gives ≥ Ψ(z₂/2)/3
smooth integers in I_P; otherwise z₂ ≥ z₁(1 + 1/33) and U1(iii) (Z = 33) gives ≥ Ψ(z₁)/99. In both cases, as all these
scales are within a factor 128 of X_r and below it, U1(i) gives ≥ cλ_Pρ_r smooth integers. With FIB(ii), Lemma 4.6,
|𝒢| ≤ M² and V² ≥ ℓS²(1 − o(1)): E_gN_r ≥ c(V²/log V)(MV/ℓ)ρ_r(M/(Vλ))/M² ≥ cS²ρ_r/(log V·λ). Averaging over ℓ loses
the factor 1 − (2k + 1)/π_J. ∎

**Proposition R2 (rational second moment; "Theorem R") [proved].** E_{g,ℓ}[N_r(N_r − 1)] ≤ (log N)^{C_2}·m_r².
*Proof.* Fix ℓ ∈ J (ℓ ∤ N) and P ≠ P′ in 𝒫 with ℓ | Δ := ab′ − a′b (a necessary condition for both to lie in Λ_ℓ(g): the
two points are annihilated mod ℓ by the nonzero form g₁x − g₀y). Then q := |Δ|/ℓ ≤ 2V²/ℓ ≤ 4S², and b ≠ b′ (if b = b′,
ℓ | (a − a′) with |a − a′| ≤ 2V < ℓ forces P = P′).
- *Injectivity and image.* On {g : ℓ | G_g(P)} the map g ↦ (n, n′) := (G_g(P)/ℓ, G_g(P′)/ℓ) is injective (its matrix
  ((−b, a), (−b′, a′)) has determinant Δ ≠ 0), with image of index q in Z². From bG(P′) − b′G(P) = −g₁Δ we get
  bn′ − b′n ≡ 0 (mod q). Since b is a prime > V/2 > 4S² ≥ q (P1), b is a unit mod q, so the image is exactly
  Γ_q := {(n, n′) : n′ ≡ λn (mod q)}, λ := b′b^{−1}, which also has index q.
- *Region.* The image of the real square [M/2, M]² is a parallelogram Π of area q·M²/(4ℓ) with sides of length ≤ √2VM/(2ℓ).
  Admissible points have |n|, |n′| ∈ [X_r/128, X_r] (filter: |G| ≥ MV/32 ⇒ |n| ≥ MV/(64Y) = X_r/128).
- *Covering.* Cover Π by axis-parallel squares of side w := M/(2V). Their number is ≤ area/w² + 4·perimeter/w + 4 ≤
  qV²/ℓ + 16√2V²/ℓ + 4 ≤ CqS². For each square I₁ × I₂ meeting the admissible set, replace I_i by an interval I_i′ of
  length w with I_i ∩ ±[X_r/128, X_r] ⊂ I_i′ ⊂ ±[X_r/256, X_r] (possible as w = X_rY/(4V²) ≤ X_r(1 + o(1))/(8S²) ≤ X_r/256, since S² ≥ 33 for N ≥ N₀ by Lemma CS).
  Apply Lemma DS with X = X_r, modulus q ≤ 4S² ≤ B^{C₀} and length w ≥ X_rB^{−C₀}, where C₀ := ⌈σ̄/a⌉ + 2
  (X_r/w ≤ 8S² ≤ B^{C₀}); for negative intervals apply it to −n. The number of points of Γ_q ∩ (I₁′ × I₂′) with both
  coordinates smooth is ≤ Σ_{ρ mod q}#{n ∈ I₁′ smooth, n ≡ ρ}·#{n′ ∈ I₂′ smooth, n′ ≡ λρ} ≤ q·((w/q)ρ_r(log N)^{C})² =
  (w²/q)ρ_r²(log N)^{2C}.
- So the number of g ∈ 𝒢 with Y_P(g) = Y_{P′}(g) = 1 is ≤ CqS²·(w²/q)ρ_r²(log N)^{2C} = CS²M²ρ_r²(log N)^{2C}/(4V²) ≤
  C(M²/ℓ)ρ_r²(log N)^{2C}, as S²/V² ≤ (1 + o(1))/(2Y) ≤ 1/ℓ.
- *Partners.* For fixed P, the P′ with ℓ | det(P, P′) lie in the lattice Λ_ℓ((a, b)) (the form bx′ − ay′), which contains the
  candidate P; by Lemma 4.4 (applied with g = (a, b); b < ℓ) it is regular and has ≤ 400S² points in [−V, V]².
- *Sum.* E_g[N_r(N_r − 1)] ≤ |𝒢|^{−1}·#𝒫·400S²·C(M²/ℓ)ρ_r²(log N)^{2C} ≤ C·3V²S²ρ_r²(log N)^{2C}/ℓ ≤ (log N)^{C₂}S⁴ρ_r². ∎

**Theorem R0 (rational supply) [proved].** γ_r := Pr_{g,ℓ}[ℓ ∤ N and N_r(g, ℓ) ≥ 1] ≥ (log N)^{−C_r}.
*Proof.* Paley–Zygmund for the nonnegative integer variable N_r under the product law:
γ_r ≥ (EN_r)²/(E[N_r(N_r − 1)] + EN_r). The right side is increasing in EN_r; with EN_r ≥ m_r/P₁, P₁ := C log V log log N
≤ (log N)², and E[N_r(N_r − 1)] ≤ P₂m_r², P₂ := (log N)^{C₂}: γ_r ≥ 1/(P₁²P₂ + P₁/m_r) ≥ (log N)^{−C_r} by Lemma CS. ∎

### 6.2 One algebraic side at one point

**Lemma SM1 (one-point smoothing) [proved].** Let P = (a, b) with gcd(a, b) = 1, ab ≠ 0, V′ := max(|a|, |b|), l_i :=
|a|^{d−1−i}|b|^i, and H ≥ 3dV′. For c uniform on C_H there is a probability measure ϑ = ϑ_P on Z with
(a) ‖Law(⟨c, v(P)⟩) − ϑ‖_TV ≤ 3d²V′/(2H + 1), v(P) := (a^{d−1−i}b^i)_{i<d};
(b) ϑ is the law of Σ_iY_i, with Y_i independent and uniform on integer intervals of lengths (2H + 1)l_i;
(c) ϑ(z) is non-increasing in |z − M₀| for some M₀ ∈ ½Z with |M₀| ≤ dV′^{d−1};
(d) max_z ϑ(z) ≤ 1/((2H + 1)V′^{d−1}).
*Proof.* The vectors u_i := be_i − ae_{i+1} (0 ≤ i ≤ d − 2) lie in K := ker⟨·, v⟩ ∩ Z^d; the maximal minor omitting column i
is ±a^{d−1−i}b^i, whose gcd over i is 1, so the u_i are a basis of K, and a fundamental cell has ℓ^∞-radius ≤ |a| + |b| ≤ 2V′.
For |n| ≤ ‖v‖₂ the real point nv/‖v‖₂² (ℓ^∞-norm ≤ 1) lies on {⟨x, v⟩ = n}; rounding it through the cell gives an integer
x with ⟨x, v⟩ = n and ‖x‖_∞ ≤ 1 + 2V′. As l_i ≤ V′^{d−1} ≤ ‖v‖₂, fix for each i a set S_i of l_i vectors with
‖x‖_∞ ≤ 3V′ on which ⟨·, v⟩ is a bijection onto {0, …, l_i − 1}. With x^{(i)} uniform on S_i, independent, put
c′ := c + Σ_ix^{(i)}; then ⟨c′, v⟩ = Σ_iY_i with Y_i := c_iv_i + ⟨x^{(i)}, v⟩ uniform on an interval of length (2H + 1)l_i
(since ±c_i is uniform on [−H, H]); the Y_i are independent. Each shift moves the uniform law on C_H by ≤ 3V′ in each of d
coordinates, so ‖Law c − Law c′‖_TV ≤ d·3dV′/(2H + 1); apply ⟨·, v⟩. (c): uniform laws on intervals are log-concave and
symmetric; convolutions of log-concave sequences without internal zeros are log-concave (Stanley 1989), hence unimodal;
symmetry gives monotonicity in |z − M₀| with M₀ = Σ(l_i − 1)/2. (d): a convolution's maximum is at most that of any
factor; take l_i = V′^{d−1}. ∎

By (c), ϑ = Σ_mW_m·Unif([M₀ − y_m, M₀ + y_m] ∩ Z) with W_m ≥ 0, ΣW_m = 1 (layer-cake decomposition).

**Lemma OP (one point, integer multiplier) [proved].** Let ℓ ∈ J, ℓ ∤ N, g ∈ 𝒢 and P ∈ 𝒫 with Y_P(g) = 1. Let w be
uniform on the integers of (W, 2W] and c uniform on C_H, independent. Then

  p(g, ℓ, P) := Pr_{w,c}[F_{g,w,c}(P) ∈ 𝒮_B] ≥ ρ_F/61.

*Proof.* Put t := |G(P)| ∈ [MV/32, 2MV], σ := sgn G(P), V′ := max(|a|, b) ∈ (V/2, V], ϖ := B and
Y_ϖ := (2H + 1)V′^{d−1}/ϖ. For fixed w write F̂_w := F̂_{g,w}(P). Since F(P) = F̂_w + G(P)⟨c, v(P)⟩ and H ≥ 3dV, SM1 gives
Pr_c[F(P) ∈ 𝒮_B] ≥ Pr_{Z∼ϑ}[F̂_w + G(P)Z ∈ 𝒮_B] − 3d²V/(2H + 1).
- *Small layers* (y_m < Y_ϖ) are supported in |z − M₀| < Y_ϖ, of ϑ-mass ≤ (2Y_ϖ + 1)/((2H + 1)V′^{d−1}) ≤ 3/ϖ by SM1(d).
- *A large layer* (y := y_m ≥ Y_ϖ, y′ := ⌊y/2⌋). As |M₀| ≤ dV′^{d−1} ≤ Y_ϖ/2 (P3: H ≥ dϖ), the layer's interval contains
  [0, y′] and [−y′, 0]. Put x₁ := ty′/4. By Lemma 4.3, |F̂_w| ≤ MHV^dN^{−η₀/4}, while
  ty′ ≥ (MV/32)(Y_ϖ/2 − 1) ≥ MHV^d/(2^{d+7}ϖ); hence |F̂_w| ≤ x₁ **for every w** by (P2). For z ∈ [0, y′] the values
  F̂_w + G(P)z run through all integers ≡ F̂_w (mod t) between F̂_w and F̂_w + σty′, so they include every n with
  σn ∈ (x₁, 2x₁] and n ≡ F̂_w (mod t). By Lemma 4.5, F̂_w ≡ c_Pw (mod t) with c_P a unit, and since W is an integer,
  Pr_w[c_Pw ≡ n (mod t)] = Pr_w[w ≡ c_P^{−1}n (mod t)] ≥ ⌊W/t⌋/W ≥ 1/t − 1/W for every n. Therefore
  E_w#{z ∈ [0, y′] : F̂_w + G(P)z ∈ 𝒮_B} ≥ (1/t − 1/W)(Ψ(2x₁) − Ψ(x₁)) ≥ (1 − t/W)x₁ϱ(x₁)/(3t),
  by U1(iii) with Z = 1. The same holds for z ∈ [−y′, 0] (values −σn). Now x₁ ≤ t·d(H + 1)V^{d−1}/4 ≤ X_F and
  X_F/x₁ ≤ 2^{d+12}dϖ ≤ B² (P2), so U1(i) gives ϱ(x₁) ≥ ρ_F/2. Hence the expected number of smooth values in the layer is
  ≥ 2(1 − t/W)(y′/4)(ρ_F/2)/3 = (1 − t/W)y′ρ_F/12, and the expected proportion is ≥ (1 − t/W)ρ_F·y′/(12(2y + 1)) ≥
  (1 − N^{−η₀/3})ρ_F/60, because y′/(2y + 1) ≥ (y − 1)/(2(2y + 1)) ≥ 1/5 for y ≥ Y_ϖ ≥ 9.
- *Total.* Summing over the large layers (total weight ≥ 1 − 3/ϖ),
  p ≥ (1 − 3/ϖ)(1 − N^{−η₀/3})ρ_F/60 − 3d²V/(2H + 1) ≥ ρ_F/61, because ϖ = B, and H/V ≥ N^{η₀/2} ≥ 10⁵d²/ρ_F (P2, P8).
  (The average over w commutes with the layer decomposition, because ϑ depends only on P and H.) ∎

*Remark.* Lemma OP is the only algebraic smoothness estimate in the proof. It uses Lemma 4.5 (a unit c_P) and the
uniformity of w modulo t, which holds because w is uniform on an interval of integers of length W ≥ tN^{η₀/3}. With w
uniform on *primes* (the design of the earlier Theorem U) only Brun–Titchmarsh's upper bound is available, and a lower bound
then requires that t be a modulus at which smooth numbers are equidistributed ("good"); that is where the restricted
Bombieri–Vinogradov Lemma 8.5′, Drappeau's Lemme 5 and Harper's preprint entered.

### 6.3 Rejections of one algebraic side

**Lemma REJ [proved].** Fix g ∈ 𝒢 with gcd(g₁, N) = 1. For (w, c) uniform as above, the probability that SAMPLE_FIELD_MF
step 2 sets acc_j = false (with no factor exit) is

  ≤ 2/(2H + 1) + e·exp(−log H/(24d log log H)) + 2^{−|𝒫_IR|} + 2^{|𝒫_IR|}/W ≤ exp(−(σ̄ + 1)Λ) ≤ ρ_F/244

for N ≥ N₀, where 𝒫_IR := primes in (16, ½ log H].
*Proof.* *Leading coefficient (b).* f_d = a_d + g₁c₀, where a_d does not depend on c; so f_d = 0 for at most one c₀:
probability ≤ 1/(2H + 1). (If f_d ≠ 0 then 0 < |f_d| < N, so its gcd check exits or passes.)
*Content (c)* is never a rejection: e | F(g₀, g₁) = Nw and e ≤ |f_d| < N, so after its gcd check gcd(e, N) = 1 and f̃ = f/e is a
primitive integral polynomial with f̃(m) ≡ 0 (mod N) and F̃(g₀, g₁) = Nw/e.
*Derivative (e).* Mod N, f′(m) ≡ f̂′(m) + g₁R_c(m, 1) since g₁m − g₀ ≡ 0, and R_c(m, 1) has coefficient 1 at c_{d−1}; g₁ is a
unit, so for fixed c₀, …, c_{d−2} at most one c_{d−1} ∈ [−H, H] gives N | f′(m) (2H + 1 < N). As f̃′(m) = f′(m)/e with e a unit
and f̃_d^{d−2} a unit after (b), the probability is ≤ 1/(2H + 1).
*Irreducibility (d).* Let 𝒫_w := {r ∈ 𝒫_IR : r | w} and 𝒫′ := 𝒫_IR ∖ 𝒫_w. Fix w. Q′ := ∏_{r∈𝒫′}r ≤ e^{(1/2)log H(1+o(1))} ≤
(2H + 1)/d, so every residue class mod Q′ of C_H has between ⌊(2H + 1)/Q′⌋^d and ⌈(2H + 1)/Q′⌉^d elements, and the law
of c mod Q′ is at most (1 + 1/d)^d ≤ e times uniform; under the uniform law the reductions mod distinct r are independent
(CRT). Fix r ∈ 𝒫′: r < B, so r ∤ N, and r ∤ w, so r ∤ Nw. If r ∤ g₁, c ↦ f mod r is an affine bijection from F_r^d onto the
polynomials of degree ≤ d with value Nwg₁^{−d} ≠ 0 at m_r := g₀g₁^{−1} (multiplication by g₁x − g₀ maps degree ≤ d − 1
bijectively onto the polynomials vanishing at m_r). If r | g₁, then g₁x − g₀ ≡ −g₀ is a unit and f mod r ranges over all
polynomials with leading coefficient a_d ≢ 0 (a_dg₀^d ≡ Nw). Monic irreducible polynomials of degree d over F_r with
prescribed nonzero value at a point number ≥ ((r^d − 1)/(r − 1) − 2r^{d/2})/d ≥ r^{d−1}/(2d) (elements of degree d of F_{r^d}
with prescribed norm, translated); summing over the r − 1 leading coefficients, the irreducible fraction is ≥ (r − 1)/(2dr) ≥
1/(4d). Hence Pr_c[no certificate | w] ≤ e(1 − 1/(4d))^{|𝒫′|}. A certificate for f is one for f̃ (same reduction up to the unit
e mod r, as r ∤ e: r | e would give r | w). If |𝒫′| ≥ |𝒫_IR|/2 this is ≤ e·exp(−|𝒫_IR|/(8d)) ≤ e·exp(−log H/(24d log log H))
(PNT: |𝒫_IR| ≥ log H/(3 log log H)). Otherwise ≥ |𝒫_IR|/2 primes of 𝒫_IR divide w; the probability of that over w is
≤ Σ_{T⊂𝒫_IR, |T|=⌈|𝒫_IR|/2⌉}(1/∏_Tr + 1/W) ≤ 2^{|𝒫_IR|}(16^{−|𝒫_IR|/2} + 1/W). *Certificate ⇒ irreducible:* a factorisation of
f̃ over Z reduces mod r to one with positive degrees, as r ∤ f̃_d.
*Size.* By (P6), the second term is ≤ e·exp(−1.25(σ̄ + 1)Λ) and the last is ≤ N^{−1/(2d)}; 2^{−|𝒫_IR|} and 2/(2H + 1) are
exp(−Ω(log H/log log H)) and N^{−h}; all are ≤ exp(−(σ̄ + 1)Λ)/4 for N ≥ N₀. Finally ρ_F ≥ exp(−(σ_F + o(1))Λ) with
σ_F < σ̄. ∎

### 6.4 The multi-field supply

For a field A = (g; (w_j, c_j)_j) and ℓ ∈ J with ℓ ∤ N, let P₀(g, ℓ) be the first P ∈ 𝒫 in the fixed order with Y_P(g) = 1
(undefined if N_r(g, ℓ) = 0). Call ℓ **MF-descendable in A** if P₀ exists and some j with acc_j has F̃_j(P₀) ∈ 𝒮_B. Let
𝒟_A be the set of such ℓ, δ_A := |𝒟_A|/π_J, and let 𝔛_A be the event that SAMPLE_FIELD_MF(A) makes a factor exit.

**Proposition MF [proved].** Under the raw law of (A, ℓ) (A from SAMPLE_FIELD_MF's coins, ℓ uniform on J, independent):

  γ_MF := Pr[𝔛_A or ℓ ∈ 𝒟_A] ≥ γ_r/2.

*Proof.* Condition on (g, ℓ). If gcd(g₁, N) > 1, 𝔛_A holds. Otherwise suppose N_r(g, ℓ) ≥ 1 (and ℓ ∤ N); P₀ is then a
function of (g, ℓ). For each j let E_j := {F_{g,w_j,c_j}(P₀) ∈ 𝒮_B and acc_j} ∪ {the steps for j make a factor exit}. The
pairs (w_j, c_j) are independent of each other and of (g, ℓ), so the E_j are conditionally independent. By Lemmas OP and
REJ, Pr[E_j | g, ℓ] ≥ Pr[F_j(P₀) ∈ 𝒮_B] − Pr[side j rejected without an exit] ≥ ρ_F/61 − ρ_F/244 ≥ ρ_F/82. So
Pr[no E_j | g, ℓ] ≤ (1 − ρ_F/82)^{J_f} ≤ exp(−J_fρ_F/82) ≤ exp(−244/82) < 1/2 by Lemma CS. If some E_j holds, either 𝔛_A holds or F_j(P₀) is smooth with acc_j, and then F̃_j(P₀) = F_j(P₀)/e_j is smooth, so
ℓ ∈ 𝒟_A. Averaging over (g, ℓ) with N_r ≥ 1, ℓ ∤ N gives ≥ γ_r/2. ∎

*Why no algebraic second moment is needed.* The bound concerns one candidate P₀ fixed by the rational data; independence
across the J_f algebraic sides turns the one-point probability into a probability of success. The counterexample to naive
"multi-field amplification" — all fields succeeding on the same target set — is excluded because the lower bound p(g, ℓ, P₀)
holds pointwise for every (g, ℓ) with N_r ≥ 1.

*Why rejected sides cannot absorb the mass.* The pointwise mass is only ≍ ρ_F = e^{−Θ(D^{2/5}λ)}. Without the MF height
floor, the irreducibility rejection e^{−log H/(13d log log H)} would be e^{−Θ(D^{2/5})} ≫ ρ_F (one factor λ short), and all of
the smooth mass could lie in rejected sides. The floor makes the rejection e^{−1.25(σ̄+1)Λ} ≪ ρ_F. Content is never a
rejection (primitive part), because its probability can be ≍ 2^{−d} ≫ ρ_F when w is even.

---

## 7. Splits, fruitfulness and field selection

**Lemma S″ (split counting without 2^{ω}) [proved].** Let N be odd, h₂ := 2^{ω(N)} = [(Z/N)^× : Q_N] and k ≥ 1. Let 𝒟 be a set
of primes coprime to N with |𝒟| ≥ 4k³h₂, and 𝒮 a set of products of exactly k distinct primes of 𝒟 lying in a window
[Y′/2, Y′] with 2Y′² < N and |𝒮| ≥ C(|𝒟|, k)/(2k). For R uniform on Q_N let E be the event that the shortest nonzero vectors
of Λ_R are ±(u, v) with u, v ∈ 𝒮 coprime. Then Pr[E] ≥ |𝒮|²/(2φ(N)).
*Proof.* (1) For a class c of (Z/N)^×/Q_N put 𝒮_c := 𝒮 ∩ c; uv^{−1} ∈ Q_N iff u, v lie in the same class, and Cauchy–Schwarz
over ≤ h₂ classes gives Σ_c|𝒮_c|² ≥ |𝒮|²/h₂. (2) For u ∈ 𝒮, at most k·C(|𝒟| − 1, k − 1) = (k²/|𝒟|)C(|𝒟|, k) ≤ (2k³/|𝒟|)|𝒮|
elements of 𝒮 share a prime with u. (3) So the number of coprime same-class pairs is ≥ |𝒮|²(1/h₂ − 2k³/|𝒟|) ≥ |𝒮|²/(2h₂).
(4) Map such a pair to R := uv^{−1} ∈ Q_N. If (u, v), (u′, v′) give the same R then uv′ ≡ u′v and |uv′ − u′v| < Y′² < N, so
uv′ = u′v and, both being coprime and positive, (u, v) = (u′, v′). (5) For such R and any shortest s ∈ Λ_R,
|det(s, (u, v))| ≤ ‖s‖‖(u, v)‖ ≤ 2Y′² < N; determinants of lattice vectors are multiples of N, so (u, v) = ms, and m = ±1
since (u, v) is primitive. (6) Pr[E] ≥ (|𝒮|²/(2h₂))/|Q_N| = |𝒮|²/(2φ(N)). ∎

**Theorem K (complete descents in one field) [proved].** Let A be a field returned by SAMPLE_FIELD_MF (no factor exit) with
δ_A ≥ γ_r/4. Then one call of DESCEND_ONCE(A) returns a row or makes a factor exit with probability
≥ (δ_A/(4 log N))^{2k}/(32k²).
*Proof.* |𝒟_A| ≥ δ_Aπ_J ≥ 4k³2^{ω(N)} by (P5) and Theorem R0. Products of k distinct primes of 𝒟_A lie in
(Y^k, (2Y)^k], covered by the k windows [Y′/2, Y′], Y′ = 2^{i+1}Y^k (0 ≤ i < k), all with Y′ ≤ (2Y)^k ≤ X and so 2Y′² < N;
one holds a set 𝒮 of ≥ C(|𝒟_A|, k)/k of them. Given that r is a unit, R is uniform on Q_N, so by Lemma S″
Pr_r[E] ≥ (φ(N)/(N − 1))|𝒮|²/(2φ(N)) ≥ |𝒮|²/(2N). With |𝒮| ≥ (|𝒟_A|/k)^k/k, |𝒟_A| ≥ δ_AY/(2 log Y), Y^{2k} ≥ X²/16^k ≥ N/(4·16^k)
and 2k log Y ≤ log N, this is ≥ (δ_A/(4 log N))^{2k}/(32k²). On E, unless a gcd check exits, |u₀| and |v₀| are products of k
distinct primes of 𝒟_A ⊂ J, so the shape test passes (FACTOR_MF returns the exact factorisation), and for each ℓ, SQ_SEARCH_MF
returns an atom (Proposition 9.2), because ℓ ∈ 𝒟_A. ∎

**Theorem DD (Dixon descent; Rabin/Dixon fibre symmetry) [proved].** Let N be odd with ω(N) ≥ 2, r₁, r₂, … independent uniform
on [1, N − 1], R_i := r_i² mod N, and ξ randomness independent of (r_i). Let ℱ := σ((R_i)_i, ξ). On an event ℰ ∈ ℱ let an
ℱ-measurable nonempty finite set 𝒮 of indices with every R_i (i ∈ 𝒮) a unit, and an ℱ-measurable 𝖸 ∈ (Z/N)^×, be given with
𝖷² ≡ 𝖸² (mod N), 𝖷 := ∏_{i∈𝒮}r_i. Then on ℰ, conditionally on ℱ, 𝖷/𝖸 is uniform on μ_N, and
Pr[1 < gcd(𝖷 − 𝖸, N) < N | ℱ] = 1 − 2^{1−ω(N)} ≥ 1/2.
*Proof.* Fix a root ρ(R) of each unit square R. Given R_i (a unit), r_i = ε_iρ(R_i) with ε_i uniform on μ_N (each unit square
has |μ_N| roots and r_i is uniform), independently of each other and of ξ. On ℰ, 𝖷/𝖸 = ∏_{i∈𝒮}ε_i·(∏ρ(R_i)/𝖸), the second
factor ℱ-measurable and in μ_N; a product of independent uniform elements of a finite group over a nonempty index set is
uniform. For ζ ∈ μ_N, gcd(ζ − 1, N) = ∏_{ζ≡1 (p_j^{e_j})}p_j^{e_j} (ζ − 1 is a unit mod p^e where ζ ≡ −1, N odd), proper unless
ζ = ±1. ∎

**Lemma 7.3 (transcript invariance) [proved].** In every round of SPLIT, the hypothesis of Theorem DD holds with 𝒮, 𝖸 the
output of POSTPROCESS_MF.
*Proof.* Everything before 𝖷 is a function of the R_i and of coins independent of the r_i: the field A and its checks; Gauss
reduction of Λ_{R_i}; the gcd checks (r_i is a unit iff R_i is); the shape test, whose input |u₀|, |v₀| is a function of R_i, whose
coins are fresh, and whose output is the canonical factorisation (so its whole internal transcript, including its own hidden
roots, is a function of (|u₀|, |v₀|, its coins) — conditional invariance under r_i ↦ εr_i is what is needed, not
independence); the special-q searches and their cache, which depend on (A, ℓ) only; and in POSTPROCESS_MF the row target, the
columns, the character primes, the kernel vector, δ_j(S), the Hensel primes and sign vectors of SQ, and 𝖸. Previous rounds'
roots go into ξ; they are independent of the new ones. Field retention biases the law of A and of the accepted R_i, not
of the fibre signs ε_i. The roots must be uniform on all of [1, N − 1] (drawing from [2, N − 2] breaks it, e.g. N = 15,
R = 1). ∎

**Lemma 7.4 (reverse Markov).** If Z ∈ [0, 1] and EZ ≥ β then Pr[Z ≥ β/2] ≥ β/2. ∎

**Proposition select (field selection) [proved].** Let a *trial* be SAMPLE_FIELD_MF followed (if no exit) by DESCEND_ONCE.
Let p be the probability that a trial ends in a factor exit or a row, C_trial a uniform bound on the expected cost of a trial and
of one DESCEND_ONCE call conditional on the past, and T_post a uniform bound on the expected cost of one POSTPROCESS_MF call.
Then (i) p ≥ (γ_r/4)(γ_r/(16 log N))^{2k}/(32k²), and (ii) E[cost of SPLIT] ≤ (2M̄_rows + 1)C_trial/p + 2T_post, where
M̄_rows is a deterministic upper bound for M_rows(A).
*Proof.* (i) Put Z_A := 1 on 𝔛_A and Z_A := δ_A otherwise. By Proposition MF, EZ_A ≥ γ_r/2, so by Lemma 7.4 the set of A with
Z_A ≥ γ_r/4 has probability ≥ γ_r/4. For such A either SAMPLE_FIELD_MF exits, or δ_A ≥ γ_r/4 and Theorem K applies.
(ii) Trials are i.i.d. until the first factor exit or row: their number is geometric with mean 1/p, and by Wald selection costs
≤ C_trial/p. Let b ≤ p be the probability of a row and ς_A the per-call row probability of DESCEND_ONCE(A). The kept field A*
has the size-biased law ς_A Pr(dA)/b, so E[1/ς_{A*}] = Pr[ς_A > 0]/b ≤ 1/b. Given A* = A, a round collects ≤ M̄_rows rows in
≤ M̄_rows/ς_A expected calls and then runs POSTPROCESS_MF; it succeeds with conditional probability ≥ 1/2 (Theorem DD and
Lemma 7.3), whatever happened before; so the expected cost given A* = A is ≤ 2(M̄_rowsC_trial/ς_A + T_post). Hence
E[cost] ≤ C_trial/p + (b/p)·2(M̄_rowsC_trial/b + T_post) ≤ (2M̄_rows + 1)C_trial/p + 2T_post. ∎

---

## 8. Simultaneous square certification for J_f fields

Fix a field A kept by SPLIT. Let 𝒥 := {j : acc_j}. For j ∈ 𝒥: f̃_j is primitive and irreducible, f̃_{j,d} and
f̃_{j,d}^{d−2}f̃_j′(m) are units mod N, and f̃_j(m) ≡ 0 (mod N). Let α_j be a root of f̃_j, K_j := Q(α_j) (degree d, odd),
β_j := f̃_{j,d}α_j, g_j(y) := f̃_{j,d}^{d−1}f̃_j(y/f̃_{j,d}) the monic minimal polynomial of β_j, D_j := f̃_{j,d}disc(f̃_j), and

  ψ_j : Z[β_j] → Z/N,  β_j ↦ f̃_{j,d}m,

a ring homomorphism (g_j(f̃_{j,d}m) = f̃_{j,d}^{d−1}f̃_j(m) ≡ 0). For an atom (P, ℓ, j), P = (a, b), put
z_P := f̃_{j,d}a − bβ_j = f̃_{j,d}(a − bα_j). Then Norm_{K_j/Q}(a − bα_j) = F̃_j(P)/f̃_{j,d} and
ψ_j(z_P) = f̃_{j,d}(a − bm) = f̃_{j,d}G(P)/g₁ (mod N).

### 8.1 The row identity

Write G(P) = σ_Pℓn_P with σ_P = ±1, n_P > 0 smooth. Then ℓ ≡ σ_Pg₁f̃_{j,d}^{−1}ψ_j(z_P)/n_P (mod N). Row i comes from
R_i = r_i² ≡ u₀/v₀ = sgn(u₀v₀)∏_ℓℓ^{e_{iℓ}}, so

  r_i² ≡ 𝔮_i·∏_ℓ f̃_{j_ℓ,d}^{−e_{iℓ}}ψ_{j_ℓ}(z_{P_ℓ})^{e_{iℓ}} (mod N),  𝔮_i := sgn(u₀v₀)∏_ℓ(σ_{P_ℓ}g₁/n_{P_ℓ})^{e_{iℓ}} ∈ Q^×.  (8.1)

For a set 𝒮 of rows, e_P(𝒮) := Σ_{i∈𝒮}e_{iP} ∈ Z (the same atom may occur in several rows), E_j(𝒮) := Σ_{P in field j}e_P(𝒮),
E(𝒮) := Σ_jE_j(𝒮).

### 8.2 The matrix

One row per row of SPLIT; columns over F₂:
- (a) the sign of 𝔮_i; (b) for each prime p ≤ B: Σ_ℓe_{iℓ}v_p(n_{P_ℓ});
- (a′) **norm sign:** Σ_ℓe_{iℓ}[F̃_{j_ℓ}(P_ℓ) < 0]; (b′) **norm parity:** for each prime p ≤ B: Σ_ℓe_{iℓ}v_p(F̃_{j_ℓ}(P_ℓ));
- for each j ∈ 𝒥: (c_j) the parity Σ_{ℓ: j_ℓ=j}e_{iℓ}; (d_j) for each prime r ≤ B with r ∤ D_j and each root s of f̃_j mod r:
  Σ_{ℓ: j_ℓ=j}e_{iℓ}[a_{P_ℓ} ≡ b_{P_ℓ}s (r)]v_r(F̃_j(P_ℓ)); (e_j) t_χ^{(j)} character columns (§8.4).

C_col := 2 + 2π(B) + Σ_{j∈𝒥}(1 + #{(r, s)}_j + t_χ^{(j)}), and **M_rows := C_col + 12**. Since #{(r, s)}_j ≤ dπ(B) and
t_χ^{(j)} = O(d² log N) (§8.4), M_rows ≤ M̄_rows := 2 + 2π(B) + J_f(1 + dπ(B) + Cd² log N) + 12 = L^{a+σ_F+o(1)}.
The kernel 𝒦 of all columns has dimension ≥ 12.

### 8.3 POSTPROCESS_MF(A, rows)

1. For each j ∈ 𝒥 touched by some row, sample character primes 𝔭^{(j)}_1, …, 𝔭^{(j)}_{t_χ^{(j)}} (§8.4); build the matrix.
2. Sample a uniform element S of 𝒦 (KS over F_{2^e}, 2^e ≥ M_rows⁴: draw x uniform in F₂^{M_rows}, give the solver only the
   right-hand side A^⊤x and independent coins, obtain y ∈ F_{2^e}^{M_rows} with A^⊤y = A^⊤x, and return S := x − π(y), where
   π: F_{2^e} → F₂ is a fixed F₂-linear map with π(1) = 1, applied coordinatewise; since A^⊤ has entries in F₂,
   A^⊤π(y) = π(A^⊤y) = A^⊤x, and given A^⊤x, x is uniform on a coset of 𝒦 independent of π(y), so S is uniform on 𝒦;
   verify). If S = 0, repeat step 2.
3. For each field j touched by S: compute δ_j(S) := the squarefree kernel, with sign, of ∏_{P in field j}F̃_j(P)^{e_P(S)}
   (from the stored factorisations), and run SQ_j on Γ′_{j,S} := δ_j(S)·g_j′(β_j)²∏_{P in field j}z_P^{|e_P(S)|} (§8.5). If
   some SQ_j returns NONSQUARE, go back to step 1 with fresh character primes.
4. Return 𝒮 := S and

   𝖸 := √(𝔮_S/∏_jδ_j(S)) · ∏_j f̃_{j,d}^{−E_j(S)/2}·ψ_j(η_j)·(ψ_j(γ⁻_{j,S})ψ_j(g_j′(β_j)))^{−1} mod N,

   with 𝔮_S := ∏_{i∈S}𝔮_i, γ⁻_{j,S} := ∏_{P in field j, e_P(S)<0}z_P^{|e_P(S)|}, η_j the output of SQ_j, and the products over
   touched j.

### 8.4 Characters at growing degree

For j ∈ 𝒥 (drop the index): S_D := primes of 𝒪_K above primes dividing D; V_D := {γ ∈ K^× : ord_𝔓γ even for all
𝔓 ∉ S_D}; ℋ := V_D/K^{×2}. For 0 ≠ 𝔮 ∈ ℋ, L_𝔮 := K(√𝔮), of degree 2d. Δ̄ := |disc g|²8^d|D|^d.

**Lemma 8.1 (sizes) [proved].** |Δ_{L_𝔮}| ≤ Δ̄; a quadratic subfield F of L_𝔮 has |Δ_F| ≤ Δ̄^{1/d}. |ℋ| ≤ 2^{dω(D)+d}h_K ≤
|ℋ|_bd := 2^d|D|^dX₀(1 + log X₀)^{d−1}, X₀ := ⌈√|disc g|⌉, and log|ℋ|_bd = O(d log N), log Δ̄ = O(d log N).
*Proof.* Δ_L = Δ_K²Norm(𝔡_{L/K}); L/K is unramified outside 2D; ord_𝔓𝔡 ≤ 1 at odd 𝔓 | D and ≤ 2ord_𝔓2 + 1 at 𝔓 | 2, so
Norm 𝔡 ≤ 8^d|D|^d; |Δ_K| ≤ |disc g|. For F ⊂ L, |Δ_F|^{[L:F]} | Δ_L. For |ℋ|: V_D ⊃ W₀ := {γ : γ𝒪_K = 𝔞²} ⊃ Y₀ := 𝒪_K^×K^{×2};
V_D/W₀ ↪ (Z/2)^{S_D}, W₀/Y₀ ↪ Cl(K)[2], |Y₀/K^{×2}| = 2^{r₁+r₂} ≤ 2^d (only ±1 are roots of unity, d odd); each class contains an
integral ideal of norm ≤ √|Δ_K| (Minkowski, crudely), and there are ≤ X(1 + log X)^{d−1} integral ideals of norm ≤ X. Sizes:
‖f̃‖_∞ ≤ 3M(H + 1) (Lemma 4.3), so log|disc f̃| ≤ (2d − 2)log(d‖f̃‖) + O(d log d) = O(d(μ + h)log N + d log d) = O(log N), as
d(μ + h) = O(1); log|disc g| = (d − 1)(d − 2)log|f̃_d| + log|disc f̃| = O(d log N). ∎

Put ν₁ := min(1/(4(2d)! log Δ̄), c_L/(2 log Δ̄)) and 𝖤 := {0 ≠ 𝔮 ∈ ℋ : ζ_{L_𝔮} has a real zero in [1 − ν₁, 1)}.

**Lemma 8.2 (the exceptional classes are rational) [proved; an extension of Lee–Venkatesan 2018 §6].** The argument (Stark's Lemma 8; d odd, so K has no quadratic subfield; F the unique quadratic subfield of L_𝔮; L_𝔮 = KF; injectivity; Landau) is that of Lee–Venkatesan, JNT 187 (2018), §6 (Lemma 6.10, Claim 6.25, Facts 6.22, 6.24, Remark 6.27), who bound the exceptional classes by (4/3) log log n and remove them by brute force. New here: |𝖤| ≤ 1 at the factorial threshold, and the observation that each exceptional class is *rational* (Kummer: 𝔮 ≡ δ_F), which lets the norm-parity columns clear it. |𝖤| ≤ 1, and every 𝔮 ∈ 𝖤 is the class of a rational number:
𝔮 = [δ] with δ ∈ Q^×.
*Proof.* Let 𝔮 ∈ 𝖤 with zero β₀. Since 1 − ν₁ ≥ 1 − 1/(4(2d)! log|Δ_{L_𝔮}|), Stark (c(L) = 4·(2d)!) gives a quadratic F ⊂ L_𝔮
with ζ_F(β₀) = 0, so L(β₀, χ_F) = 0. As [L_𝔮 : Q] = 2d with d odd, L_𝔮 has at most one quadratic subfield (two would generate
a biquadratic subfield, 4 ∤ 2d) and K has none; so L_𝔮 = KF, i.e. K(√𝔮) = K(√δ_F). By Kummer theory 𝔮 ≡ δ_F (mod K^{×2}).
If 𝔮 ≠ 𝔮′ were both in 𝖤, then F_𝔮 ≠ F_𝔮′ (Kummer: 𝔮 ↦ L_𝔮 is injective), and χ_{F_𝔮}, χ_{F_𝔮′} would be distinct real
primitive characters of conductors ≤ Δ̄^{1/d} with real zeros ≥ 1 − ν₁ > 1 − c_L/log(Δ̄^{2/d}), contradicting Landau. ∎

**Character sampler.** Let s := ⌈(2d)!·(log N + log Δ̄ + d)⁴⌉ and x := e^{s}. Draw r uniform in (x/2, x]; if r is prime (AKS),
find the roots Z_r of f̃ mod r (CZ); with probability |Z_r|/d output (r, s′) with s′ uniform in Z_r, else redraw. The output
is uniform over the degree-one primes 𝔭 = (r, α − s′) of norm in (x/2, x]; r > Δ̄ > |D|, r > B, and r exceeds every
|F̃(P)|, so no atom lies in 𝔭. The expected number of draws is O(ds): by LO for K/K (trivial group, n_L = d) there are
≥ x/(3 log x) prime ideals of norm in (x/2, x] — a real zero of ζ_K in [1 − 1/(4d! log|Δ_K|), 1) would, by Stark, give a
quadratic subfield of K (none), and a zero below contributes ≤ x^{1−1/(4d! log Δ̄)} = o(x/log x) since s ≥ (2d)!(log Δ̄)⁴;
those of degree ≥ 2 number ≤ d√x.
For such 𝔭, χ_𝔭(𝔮) := the Legendre symbol mod 𝔭 of γπ^{−ord_𝔭γ} (γ ∈ V_D a representative, π a uniformiser at 𝔭; ord_𝔭γ is even)
is a well-defined quadratic character of ℋ, and χ_𝔭(a − bα) = ((a − bs′)/r) on atoms. The character columns are
Σ_ℓe_{iℓ}(1 − χ_𝔭(a_{P_ℓ} − b_{P_ℓ}α))/2 mod 2 over the atoms of field j.

**Lemma 8.3 (characters separate) [proved].** For every 𝔮 ∈ ℋ ∖ (𝖤 ∪ {0}), Pr_𝔭[χ_𝔭(𝔮) = −1] ≥ 1/(5d) for N ≥ N₀. With
t_χ^{(j)} := ⌈5d(ln|ℋ_j|_bd + ln(4J_f))⌉ = O(d² log N) independent samples, the probability that
ℋ_{j,0} := ℋ_j ∩ ⋂_t ker χ_{𝔭_t^{(j)}} is not contained in 𝖤_j ∪ {0} is ≤ 1/(4J_f); otherwise |ℋ_{j,0}| ≤ 2.
*Proof.* Apply LO to L_𝔮/K, G = C₂, 𝒞 = {σ}, n_L = 2d, d_L ≤ Δ̄. LO's range condition: s ≥ (2d)!(log Δ̄)⁴ ≥ 20d(log Δ̄)². For 𝔭 ∤ 2D,
χ_𝔭(𝔮) = −1 iff Frob_𝔭 = σ. The main term gives ½(Li(x) − Li(x/2)) ≥ (x/(4 log x))(1 − o(1)) such prime ideals with norm in
(x/2, x]. A real zero, if any, is < 1 − ν₁ (𝔮 ∉ 𝖤) and contributes ≤ 2x^{1−ν₁} = 2x·e^{−ν₁s} = o(x/log x), since
ν₁s ≥ (log N + log Δ̄)⁴/(4 log Δ̄)·min(1, 2(2d)!c_L) ≫ log s = O(d log d + log log N). The error term is
2c₁x·exp(−c₂√(s/(2d))) = o(x/log x). Prime ideals of degree ≥ 2 (≤ d√x) and ramified ones are negligible. There are ≤ d(π(x) −
π(x/2)) ≤ 0.6dx/log x degree-one primes in all, so the ratio is ≥ 1/(5d). Independence: Pr[𝔮 ∈ ℋ_{j,0}] ≤ (1 − 1/(5d))^{t_χ}; a
union bound over ℋ_j gives ≤ |ℋ_j|_bd·e^{−t_χ/(5d)} ≤ 1/(4J_f). Finally ℋ_{j,0} is a subgroup contained in 𝖤_j ∪ {0}, which
has ≤ 2 elements. ∎

*Cost.* Each character prime has s = exp(O(d log d))(log N)^{O(1)} bits; sampling and each Legendre symbol cost s^{O(1)}.
Field j needs t_χ^{(j)} primes (O(ds) draws each) and one Legendre symbol per (character, atom of field j). In total
≤ (J_f + 2kM_rows)·max_jt_χ^{(j)}·s^{O(1)} = M̄_rows·e^{O(d log d)} = L^{a+σ_F+o(1)}, by (P4).

### 8.5 Square roots

**Algorithm SQ_j(Γ)** (after Buhler–Lenstra–Pomerance 1993 §9). Input Γ = δ·g′(β)²∏_Pz_P^{|e_P|} ∈ Z[β], δ ∈ Z squarefree and
B-smooth. (1) Let p be the least prime in (B, 2B] with p ∤ D (the primes of (B, 2B] are sieved once per SPLIT; at most
log|D|/log B of them divide D). Then g mod p is squarefree and Γ is a unit at p (its norm has prime factors ≤ B or dividing
D). (2) Factor g mod p = φ₁⋯φ_τ (CZ); in each F_p[y]/(φ_i) test whether Γ is a square and compute a root ρ_i (Tonelli–Shanks); if
some component is a non-square output NONSQUARE. (3) For each of the 2^{τ−1} sign vectors ε with ε₁ = +1: CRT-combine ε_iρ_i to
ρ_ε ∈ Z[β]/p, Newton-lift the inverse square root z ↦ z + z(1 − Γz²)/2 from ρ_ε^{−1} to precision p^e > 2A_sq, put
η_ε := Γz mod p^e (centred); if η_ε² = Γ exactly, output η_ε. (4) Output NONSQUARE.

**Lemma 8.4 (square roots) [proved].** Let Mh := ⌊‖g‖₂⌋ + 1, g_d := d‖g‖₁Mh^{d−1}, U := |δ|g_d²∏_P(‖z_P‖₁Mh^{d−1})^{|e_P|},
A_sq := d⌈U^{1/2}⌉2^{d−1}Mh·g_d^{d−1}. (i) Any output η satisfies η² = Γ; if Γ ∈ K^{×2}, then Γ = η² with η ∈ Z[β],
‖η‖_∞ ≤ A_sq, and SQ outputs ±η. So SQ decides exactly whether Γ is a square. (ii) log A_sq = O((1 + Σ_P|e_P|)d log N + d² log N),
and SQ runs in expected 2^{d}·Õ(log A_sq)·(log N)^{O(1)} bit operations.
*Proof.* (i) If Γ = η², then (η/g′(β))² = δ∏z_P^{|e_P|} ∈ 𝒪_K, so η/g′(β) ∈ 𝒪_K, and g′(β)𝒪_K ⊂ Z[β] (𝒪_K ⊂ 𝒪_K^∨ ⊂ Z[β]^∨ =
g′(β)^{−1}Z[β], Euler). As p ∤ disc g, Z_p[β] ≅ ∏_iZ_p[y]/(G_i) with unramified factors lifting the φ_i; the reduction of η is ±ρ_ε for
exactly one ε with ε₁ = +1, a unit square root with prescribed residue lifts uniquely (p odd), so η_ε ≡ ±η (mod p^e), and centring
recovers ±η since p^e > 2A_sq. NONSQUARE in step 2 is correct since η would reduce to a square in each residue field. Bound:
|σ_k(z)| ≤ ‖z‖₁Mh^{d−1}, |g′(β_k)| ≤ g_d, so |σ_k(η)| ≤ U^{1/2}; recover coefficients with the inverse Vandermonde matrix,
whose numerators are coefficients of ∏_{i≠k}(y − β_i), bounded by 2^{d−1}Mahler(g) ≤ 2^{d−1}Mh, and |g′(β_k)| ≥ g_d^{−(d−1)} since
∏_k|g′(β_k)| = |disc g| ≥ 1. (ii) |δ| ≤ ∏|F̃(P)|^{|e_P|}; product tree, Newton iteration with O(d²) multiplications per doubling,
2^{τ−1} ≤ 2^{d−1} sign vectors. ∎

### 8.6 The certification theorem

**Lemma 8.5 (residual obstructions are rational) [proved].** Let S ∈ 𝒦 and j ∈ 𝒥. Then E_j(S) is even, γ_{j,S} :=
∏_{P in field j}(a_P − b_Pα_j)^{e_P(S)} lies in V_{D_j}, and φ_j(S) := [γ_{j,S}] ∈ ℋ_{j,0}. If ℋ_{j,0} ⊂ 𝖤_j ∪ {0}, then
δ_j(S)·γ_{j,S} ∈ K_j^{×2}.
*Proof.* Columns (c_j) give E_j(S) even. For r ∤ D_j, r ≤ B: r ∤ f̃_{j,d}, so α_j is r-integral and Z_{(r)}[α_j] has discriminant
disc(f̃_j)f̃_{j,d}^{−(2d−2)}, an r-unit; it is the maximal order at r and (Dedekind–Kummer) the primes above r correspond to the
irreducible factors of f̃_j mod r. For coprime (a, b): if r | b, then F̃_j(a, b) ≡ f̃_{j,d}a^d ≢ 0 and a − bα_j is a unit above r;
if r ∤ b, the only prime above r with ord(a − bα_j) ≠ 0 is the degree-one prime (r, α_j − ab^{−1}), with
ord = v_r(F̃_j(a, b)). So columns (d_j) give even valuations at all primes not above D_j and ≤ B; primes > B with r ∤ D_j do
not occur (F̃_j(P) is smooth). Hence γ_{j,S} ∈ V_{D_j}; φ_j is F₂-linear (γ_{S+S′} = γ_Sγ_{S′}/γ_{S∩S′}²), and the character
columns give φ_j(S) ∈ ℋ_{j,0}. If ℋ_{j,0} ⊂ 𝖤_j ∪ {0}, Lemma 8.2 gives γ_{j,S} = δξ² with δ ∈ Q^×, ξ ∈ K_j^× (δ = 1 if
φ_j(S) = 0). Then Norm(γ_{j,S}) = δ^dNorm(ξ)² ≡ δ (mod Q^{×2}) as d is odd. On the other hand
Norm(γ_{j,S}) = ∏_P(F̃_j(P)/f̃_{j,d})^{e_P(S)} ≡ ∏_PF̃_j(P)^{e_P(S)} (E_j(S) even), whose squarefree kernel is δ_j(S). So
δ_j(S) ≡ δ (mod Q^{×2}) and δ_j(S)γ_{j,S} = δδ_j(S)ξ² ∈ K_j^{×2}. ∎

*Remark (the J_f-functional form).* On the event ℋ_{j,0} ⊂ 𝖤_j ∪ {0} for all j, b_j(S) := [φ_j(S) ≠ 0] is an F₂-linear functional
on 𝒦, so the residual obstructions define at most |𝒥| ≤ J_f linear functionals, and a common zero S ≠ 0 exists once
dim 𝒦 ≥ J_f + 1. This is the route of the earlier outline (lv5unif §11): it needs J_f + 2 kernel vectors and a J_f × (J_f + 2)
bit matrix, i.e. J_f extra linear-algebra solves. Lemma 8.5 shows the obstructions are moreover *rational* classes detected by the
norm, so the π(B) + 1 shared columns (a′), (b′) annihilate all of them simultaneously, and one kernel vector suffices.

**Proposition PP (post-processing) [proved].** Let A be a kept field and M_rows rows be given. POSTPROCESS_MF is Las Vegas,
returns a nonempty S and 𝖸 ∈ (Z/N)^× with ∏_{i∈S}r_i² ≡ 𝖸² (mod N), and every quantity it computes is a function of A, of
the rows' certificates (R_i, u₀, v₀, atoms, factorisations) and of its own coins. Its expected cost is
T_post ≤ M̄_rows²·L^{o(1)} + B^{1+o(1)} = L^{2(a+σ_F)+o(1)}.
*Proof.* *Kernel.* Each row has ≤ 2k atoms, and each atom contributes to ≤ 4 log₂N + 4 + t_χ^{(j)} columns, so the matrix
has nnz ≤ M_rows·2k(4 log₂N + 4 + max_jt_χ^{(j)}) = M_rows·L^{o(1)} nonzero entries; KS returns a uniform kernel element (given A^⊤x, x is uniform on a coset of the
kernel and the solver sees only A^⊤x), verified, in Õ(M_rows²·k log N + M_rows·Σt_χ) operations; dim 𝒦 ≥ 12 makes S ≠ 0 with
probability ≥ 1 − 2^{−12}.
*Squares.* With probability ≥ 3/4 (Lemma 8.3, union over j), ℋ_{j,0} ⊂ 𝖤_j ∪ {0} for every j; then by Lemma 8.5,
δ_j(S)γ_{j,S} ∈ K_j^{×2}. Since γ′_{j,S} := ∏_Pz_P^{e_P(S)} = f̃_{j,d}^{E_j(S)}γ_{j,S} with E_j(S) even,
Γ′_{j,S} = δ_j(S)g_j′(β_j)²γ′_{j,S}(γ⁻_{j,S})² is a square, and SQ_j returns η_j (Lemma 8.4). A failed round has probability
≤ 1/4 + 2^{−12}; expected O(1) rounds.
*Rational parts.* Columns (a), (b) and E(S) = Σ_jE_j(S) even make 𝔮_S = g₁^{E(S)}∏_Pn_P^{−e_P(S)}·(sign) a positive rational square,
with explicit root g₁^{E(S)/2}∏_{p≤B}p^{−½Σ_Pe_P(S)v_p(n_P)}. Columns (a′), (b′) make ∏_PF̃_{j_P}(P)^{e_P(S)} a positive square;
since δ_j(S) ≡ ∏_{P in j}F̃_j(P)^{e_P(S)} (mod Q^{×2}) and the δ_j are squarefree, ∏_jδ_j(S) is a positive square. So
𝔮_S/∏_jδ_j(S) is a positive rational square whose root is read off the exponent vectors.
*Congruence.* ψ_j is a ring homomorphism, so ψ_j(η_j)² = δ_j(S)ψ_j(g_j′(β_j))²∏_{P in j}ψ_j(z_P)^{|e_P(S)|}. Hence
𝖸² = 𝔮_S∏_jf̃_{j,d}^{−E_j(S)}∏_{P in j}ψ_j(z_P)^{|e_P(S)|}/ψ_j(γ⁻_{j,S})² = 𝔮_S∏_jf̃_{j,d}^{−E_j(S)}∏_{P in j}ψ_j(z_P)^{e_P(S)}
≡ ∏_{i∈S}r_i² by (8.1). Everything inverted is a unit mod N: ψ_j(z_P) = f̃_{j,d}G(P)/g₁ with G(P) = ±ℓn_P coprime to N;
ψ_j(g_j′(β_j)) = f̃_{j,d}^{d−2}f̃_j′(m) (checked); δ_j(S) has only prime factors ≤ B; g₁, f̃_{j,d} are checked.
*Transcript.* No step reads a root r_i.
*Cost.* Kernel: M̄_rows²L^{o(1)}. Characters: (2kM_rows + J_f)·max_jt_χ^{(j)}·e^{O(d log d)} = M̄_rows·L^{o(1)} (§8.4). (r, s)-columns: J_f·π(B)·poly(d, log B) (root
finding mod r ≤ B, done once per kept field). SQ: Σ_j over touched fields of 2^{d}Õ(log A_sq^{(j)})(log N)^{O(1)} with
Σ_jΣ_P|e_P| ≤ 2kM_rows, i.e. ≤ M̄_rows·L^{o(1)}; the Hensel-prime sieve B^{1+o(1)}. ∎

---

## 9. Batch smoothness, the special-q search, recursion

**Lemma 9.1 (smooth parts; Bernstein) [proved].** Let 𝒫_B := ∏_{p≤B}p, x ≥ 2, 2^e ≥ log₂x. The B-smooth part of x is
gcd(x, (𝒫_B mod x)^{2^e} mod x). For x₁, …, x_n all smooth parts cost O(𝖬(b)log b) bit operations, b := log₂𝒫_B + Σlog₂x_i.
*Proof.* x = s·u with s smooth, gcd(u, 𝒫_B) = 1; every p^i ∥ s has i ≤ log₂x ≤ 2^e, so s | 𝒫_B^{2^e} and gcd(u, 𝒫_B^{2^e}) = 1.
Product tree of the x_i, reduce 𝒫_B modulo the product, remainder tree (von zur Gathen–Gerhard §§10.1–10.3); then e
squarings and a gcd per x_i. ∎

**Proposition 9.2 (the special-q search) [proved].** SQ_SEARCH_MF(A, ℓ) is a deterministic function of (A, ℓ) and costs
Õ(B + S²log N + J_f·d log N). If ℓ ∤ N, it either makes a factor exit or returns an atom (P₀, j) iff ℓ ∈ 𝒟_A, and every
returned atom has P₀ ∈ 𝒫 ∩ Λ_ℓ(g), |G(P₀)| ≥ MV/32, G(P₀)/ℓ ∈ 𝒮_B, acc_j and F̃_j(P₀) ∈ 𝒮_B.
*Proof.* If some P ∈ 𝒫 lies in Λ_ℓ(g), step 1 does not fail (Lemma 4.4); step 2 enumerates ≤ 400S² points from the reduced basis
and finds exactly the first P with Y_P = 1, i.e. P₀(g, ℓ). Step 3 tests exactly the accepted j. Smooth parts are exact
(Lemma 9.1, all values O(log N) bits by Lemma 4.3). ∎

**Lemma OS (oracle substitution) [proved].** Let SPLIT^O denote SPLIT with shape test O, where O on input n returns the canonical
(sorted) factorisation of n, terminates a.s., and uses fresh coins. The joint law of the outer logical transcript (outer coins
and every value computed outside O) is the same for every such O.
*Proof.* The outer computation is a deterministic function of outer coins and O's return values, and every return value is the
factorisation of its input, a function of the input alone. Induction over the steps. ∎

**Lemma W (Wald for adaptive sub-calls) [proved].** If a procedure makes sub-calls i = 1, 2, …, the decision to make call i
and its input are measurable w.r.t. the history before it, call i uses fresh coins, and E[C_i | history, input] ≤ F̄, then
E[Σ_{i≤τ}C_i] ≤ F̄·E[τ]. *Proof.* EΣC_i1{τ ≥ i} = ΣE[1{τ ≥ i}E[C_i | ·]] ≤ F̄ΣPr[τ ≥ i]. ∎

**Lemma F_α (solving the recursion) [proved].** Let 0 < α < 1, Λ_α(x) := (log x)^α(max(1, log log x))^{1−α}, L_x := e^{Λ_α(x)}.
If F: [2, ∞) → [0, ∞) is nondecreasing and finite, and F(x) ≤ L_x^{p+ε(x)} + L_x^{δ+ε(x)}F(√(2x)) for x ≥ x₀ with ε(x) → 0, then
F(x) ≤ L_x^{c+o(1)}, c := max(p, δ/(1 − 2^{−α})).
*Proof.* For y = √(2x): (log y)^α ≤ 2^{−α}(log x)^α + O((log x)^{α−1}) and log log y ≤ log log x, so Λ_α(y) ≤ 2^{−α}Λ_α(x) + o(1)
and L_y^b ≤ e·L_x^{2^{−α}b} for x ≥ x₁(b). Fix η > 0, η′ := η(1 − 2^{−α})/4, and X ≥ max(x₀, x₁(c + η), 4) with ε ≤ η′,
L_x^{η−η′} ≥ 2 and L_x^{3η′} ≥ 2e on [X, ∞); K := max(1, F(X)). The bound F ≤ KL_x^{c+η} holds on [2, X]. If it holds on [2, Y′],
Y′ ≥ X, then for x ∈ (Y′, Y′²/2]: F(x) ≤ L_x^{p+η′} + Ke·L_x^{δ+η′+2^{−α}(c+η)}, and δ ≤ c(1 − 2^{−α}) makes the second exponent
≤ c + η − 3η′, so both terms are ≤ (K/2)L_x^{c+η}. ∎

---

## 10. Cost ledger and proof of Theorem U-MF

### 10.1 Ledger (expected costs, exponents of L_N[2/5])

| Item | Bound | Exponent |
|---|---|---|
| preprocessing: trial division to B, S-scan, J_f (Lemma CS) | B^{1+o(1)} | a |
| SAMPLE_FIELD_MF: J_f sides, each (d log N)^{O(1)} (IR via CZ on r ≤ ½ log H) | J_f·L^{o(1)} | σ_F |
| SQ_SEARCH_MF (Prop 9.2) | Õ(B + S² log N + J_f d log N) | t_q := max(a, σ_r, σ_F) |
| one DESCEND_ONCE, excluding shape tests: 2k searches + Gauss + gcds | L^{t_q+o(1)} | t_q |
| shape tests per DESCEND_ONCE | 2F(√(2N)) | (recursive) |
| C_trial | 2F(√(2N)) + L^{t_q+o(1)} | |
| 1/p (Prop select (i), γ_r ≥ (log N)^{−C_r}) | exp((2k + 1)((1 + C_r)λ + O(1)) + O(log k)) | s = 2κ(1 + C_r) |
| M̄_rows (§8.2) | J_f(1 + dπ(B) + O(d² log N)) + 2π(B) + 14 | a + σ_F |
| T_post (Prop PP): one kernel vector, characters, SQ in the touched fields | M̄_rows²L^{o(1)} + B^{1+o(1)} | 2(a + σ_F) |
| **the MF height floor** | multiplies H, hence X_F, by e^{30(σ̄+1)dΛλ}; changes log(1/ρ_F) by ≤ (30(σ̄+1)dΛλ/log B)·log u_F(1 + o(1)) = O(D^{1/5}λ²) = o(Λ) | +0 |
| character primes (factorial bits) | e^{O(d log d)} per operation | +0 (d log d = o(Λ)) |
| SQ sign vectors | 2^{d−1} per field | +0 |
| J_f's slack factor (log N)^{λ₃} | e^{λλ₃} | +0 |

So C_γ := C_r in the statement (we write γ := γ_r/2 for the per-prime success probability).

### 10.2 The recursion

Let F(x) := sup_{2≤n≤x}E[cost of FACTOR_MF(n)]. *Finiteness:* because of the alternation with trial division, every node of
the factor tree costs ≤ 2(√n + O(1)) + (cost of preprocessing), so F(x) < ∞. *Accounting:* FACTOR_MF(N) builds a factor tree
with ≤ 2 log₂N nodes; each node n | N costs L_n^{a+o(1)} for AKS, perfect powers and trial division, plus (at most twice, by the
alternation) one SPLIT on its cofactor. Inside SPLIT, by Proposition select (ii) with the ledger, Lemma W for the shape tests (their
inputs are ≤ √(2N), their coins fresh) and Lemma OS (the law of SPLIT's outer transcript does not depend on which exact
factoring oracle is used, so the bounds proved for an ideal oracle hold for O = FACTOR_MF), for x ≥ x₀:

  F(x) ≤ 2 log₂x·[P(x) + D(x)F(√(2x))],  P(x) = L_x^{max(2(a+σ_F), a+σ_F+s+t_q)+o(1)},  D(x) = 8(2M̄_rows + 1)/p = L_x^{a+σ_F+s+o(1)}.  (10.1)

(The factor 8 collects the alternation (2), the two shape tests per trial (2) and Proposition select (ii) (2).)

(Suprema over n ≤ x of bounds L_n^{q+ε(n)} with ε → 0 are again L_x^{q+o(1)}.)

**Proof of Theorem U-MF.** Lemma F_α with α = 2/5, p = max(2(a + σ_F), a + σ_F + s + t_q) and δ = a + σ_F + s, applied to (10.1)
(the 2 log₂x is L^{o(1)}), gives F(x) ≤ L_x[2/5, c_MF + o(1)]. *Correctness:* every output prime is AKS-certified, and termination
is almost sure because the expectation is finite. *Las Vegas:* no output is ever wrong; only the running time is random. ∎

### 10.3 The constant

The third term always dominates (codex/mfU_2): max(a, σ_r, σ_F) ≤ a + σ_F, so the first two terms are ≤ 2(a + σ_F) + s
< (a + σ_F + s)/(1 − 2^{−2/5}), as 1/(1 − 2^{−2/5}) > 4.1. Hence c_MF = (a + σ_F + s)/(1 − 2^{−2/5}).


C_γ = C_r is assembled in Theorem R0 from C_R (Lemma CS, from HT1), the constant C_U of Lemma U1 (from HT2, HT3), C₀ = ⌈σ̄/a⌉ + 2
in Lemma DS, and the constant of Lemma FIB; it depends only on (a, κ, d₀) and absolute constants. These come from published
theorems with effective but unspecified implied constants, so c_MF is finite but not numerically known. With 2^{−2/5} bracketed by
exact rationals (`src/mfU/params.py`), the illustrative optima of c_MF over a rational grid are recorded in `notes/mfU_results.md`
§4.

---

## 11. Dependency graph and inputs

### 11.1 Dependency graph (→ means "is used by")

```
HT1 ─────────────→ Lemma CS ──────────────┐
HT2, HT3, CEP ───→ Lemma U1 ─→ Lemma DS ─→ Prop R2 ─┐
                       │                             ├→ Theorem R0 (γ_r) ─┐
PNT, Lemma 4.6 ──→ Prop R1 ←── Lemma FIB ────────────┘                    │
Lemma 4.4 (reg) ─→ Prop R2, Prop 9.2                                      │
Lemmas 4.2, 4.3, 4.5 ─→ Lemma OP ←── Lemma SM1 (Stanley), Lemma U1        │
Lemma REJ (PNT, CRT) ──────────────┐                                      │
Lemma OP ──────────────────────────┴→ Proposition MF ←────────────────────┘
Proposition MF + Lemma S″ ─→ Theorem K ─→ Proposition select
Stark, Landau ─→ Lemma 8.2 ─┐
LO, Lemma 8.1 ─→ Lemma 8.3 ─┼→ Lemma 8.5 ─→ Proposition PP ←── Lemma 8.4 (SQ), KS
Dedekind–Kummer, Gauss ─────┘
Theorem DD + Lemma 7.3 (+ Lemma OS) ─→ Proposition select (ii)
Lemma 9.1 ─→ Proposition 9.2 ─→ Theorem K, ledger
Proposition select + Prop PP + Lemmas OS, W, F_α + AKS ─→ Theorem U-MF
```

### 11.2 Where unrefereed or unpublished inputs enter

1. **Every internal lemma** (§§2–10) is unrefereed: it has been checked only by the author agent and by AI reviewers. The ones a
   referee should read first are Lemma OP (the only algebraic smoothness estimate), Lemma REJ (the height floor, integer
   multipliers in the irreducibility certificate), Proposition R2 (the rational second moment), Lemma 8.5 with Proposition PP
   (simultaneous certification), and Lemma 7.3 (transcript invariance with the recursive shape test).
2. **No external unpublished input.** In particular Drappeau's Lemme 5 (Compositio 151 (2015)), which is stated there without
   proof as a consequence of §3.3 of Harper's unpublished arXiv:1208.5992, is **not used**. Neither is any other
   Bombieri–Vinogradov-type estimate.
3. **External published inputs whose labels we could not check against the printed source:** Lagarias–Odlyzko's theorem
   number (content checked against two restatements; §6 of the results note). CEP's label (unnumbered Corollary after Thm 3.1), HT Thm 2 "(i)", Stanley's Prop. 2 are likewise UNVERIFIED as labels. The rest were checked as stated (HT Thms 1, 3 via the HT 1993 survey; Rosser–Schoenfeld; Stark Lemma 8 verbatim including F ⊂ L in the non-normal case; AKS; CZ).
4. **Effectivity.** All cited constants are effective (HT, CEP, Stark, Landau, LO); no Siegel-type ineffective bound is used. The
   algorithm uses no analytic constant at all (alternation with trial division; J_f uses the slack factor (log N)^{λ₃}).

### 11.3 Variants

- *Prime multipliers and good moduli.* If w is uniform on primes (as in the project's Theorem U), Lemma OP must be replaced by a
  lower bound at "good" moduli t (smooth numbers equidistributed mod t), plus a restricted Bombieri–Vinogradov estimate showing
  that most t are good; that brings back Drappeau's Lemme 5 and Harper's preprint. The integer multiplier makes this unnecessary.
- *Single field.* The project's Theorem U uses one algebraic side and an algebraic second moment (LD*, Kern-U, L2_fam, Lemma W).
  With integer multipliers its first moment would also follow from Lemma OP and avoid Lemma 8.5′, but its second moment is not
  re-proved here.

### 11.4 Why the exponent is 2/5 [heuristic arithmetic for this design; not a lower bound for factoring]

The split needs (log N)^{−2k}-scale acceptance, so cost L[α] forces k = O(D^{α}); the total norm exponent (rational plus
algebraic, θ = 2/d + dν + o(·); the algebraic part alone is 1/d + dν) is ≥ c/√k at the optimum d ≍ √k and smoothness of N^{θ}-size values at B = L[α] costs exp(u log u) with u = θ log N/log B, which is L[α]^{O(1)}
only if θ = O(D^{2α−1}). Together k ≳ D^{2(1−2α)} and k = O(D^{α}), so α ≥ 2/5. Theorem U-MF attains it.
