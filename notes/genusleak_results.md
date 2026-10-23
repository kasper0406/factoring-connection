# genusleak: family (Q) without smoothness (phase 24)

Agent `genusleak`, 2026-10-23. Pre-registration: `notes/genusleak_prereg.md` (commit fdef700). Code:
`src/genusleak/` (operational: `forms.py`, `cf.py`, `cayley.py`, `sqrtred.py`, `genus.py`; harness: `harness.py`).
Tests: `tests/test_genusleak.py`. Raw data: `experiments/raw/genusleak_E{1..5}.jsonl`.

Throughout, N = pq with p < q odd primes, k ≥ 1 squarefree with gcd(k, N) = 1, D = −4kN (imaginary case) or
D = 4kN (real case), and h = h(D). "Hidden 2-torsion" means a class of Cl(D)[2] outside the subgroup A_pub generated
by the classes that k alone determines (Lemma 0).

## 0. Verdict (plain language)

**Negative. We found no subexponential route to a hidden ambiguous class that avoids a group order, a relation, or
a smoothness event. Every route examined either leads back to one of those or is exponential.**

- The genus characters are not free. Deciding whether a class is a square (whether it lies in the principal genus)
  is *equivalent* to quadratic residuosity mod N (Prop GEN). That problem is not known to be easy, nor known to be as
  hard as factoring.
- **Square roots in Cl(D) are factoring-equivalent** (Prop SQ, known).
- **Random-walk statistics on the Cayley graph of small prime forms carry information only through relations.**
  - Every closed walk is a smooth value x² + kNy² ≥ kN (Lemma GIRTH). Walks shorter than log(kN)/log ℓ_max carry no
    information at all.
  - The hidden genus term in the return probability needs ≥ δ²|H|/4 walks to detect (Prop WALK). Its eigenvalue is
    the mean of (ℓ/p) over the generators, a QR-type quantity.
- **Infrastructure.**
  - Any exact multiple of the regulator (hR is enough; h need not be separated from R) gives a factor in polynomial
    time, by a 2-adic sweep (Prop HR).
  - An analytic approximation of hR to absolute precision Δ buys exactly a √Δ search. The Euler-product route
    therefore stops at |D|^{1/5} under GRH, as Lenstra and Schoof found.
  - SQUFOF's square forms are the 2-torsion classes' shadows in the principal cycle. They have density ≍ D^{−1/4}, and
    giant steps cannot raise that density (E2).
- **Isogeny volcanoes** fail at construction. A curve with Frobenius discriminant −4kN·f² needs either a
  degree-h(O) ≥ N^{1/2+o(1)} class polynomial or about √(kN) random trials. Even with the curve, ambiguous classes act
  by isogenies of the unknown degree p.
- **Ring class groups and multiplier pencils** only add more groups whose orders are hidden. The kernels have known
  order, but that order is useless against the hidden odd part.
- **Barrier (Theorem GB).** In a generic model of a hidden-order group, a q-query algorithm outputs a hidden 2-torsion
  element with probability O(q³/X) when the odd order has a prime factor P that is uniform-like in [X, 2X].
  - Any successful generic algorithm hands you a relation (2e ∈ Λ, e ∉ Λ). In rank 1 that is a multiple of an
    element's odd order: the route *is* a group-order problem.
  - The non-generic features that escape the model are listed in §6. Every known escape is one of three kinds:
    - a multiplicative coincidence of reduced-form norms: smooth (relations, i.e. (M)), square (SQUFOF, N^{1/4}), or a
      collision;
    - an analytic order approximation (N^{1/5}, GRH);
    - a direct divisor search ((A)).
  - The open class is precisely: non-generic uses of the integers of reduced forms that are none of these.
- **Structural remark.** h(−4kN) is an Euler product of *public* local factors (D/ℓ). The hidden characters (ℓ/p)
  never enter h individually; they enter only the group *structure*. It is the public order that unlocks the hidden
  2-torsion. So, as far as we can see, family (Q) reaches the factors through (M)'s door.

Nothing here is a claimed novelty. Labels: PROVED (proof here), KNOWN (reference), OBSERVATION, CONJECTURED.

## 1. Background facts used (KNOWN)

- **(F1) Ambiguous reduced forms** (elementary; Buell, *Binary Quadratic Forms* (1989) Ch. 3; Cox, *Primes of the
  form x² + ny²*, §3).
  - A reduced form f = (a, b, c) has [f]² = 1 iff f ~ (a, −b, c). By uniqueness of reduced forms this happens iff
    b = 0, a = b or a = c.
  - So the reduced form of each class in Cl(D)[2] is ambiguous; there is exactly one per class.
  - The reduced ambiguous forms are enumerable from the divisors of D (tested: `genus.ambiguous_reduced_forms` against
    PARI's 2-rank).
- **(F2) Gauss genus theory** (Cox Thm 3.15, Thm 6.1; Gauss DA §§229–287).
  - Cl/Cl² ≅ (Z/2)^{μ−1}, cut out by the assigned characters.
  - For D = −4kN these are (·/p), (·/q), (·/ℓ) for odd ℓ | k, and the 2-adic characters among χ₋₄, χ₈, χ₋₈ (which
    ones depends on kN mod 8).
  - The principal genus equals Cl² (Gauss's duplication theorem).
- **(F3) The class number formula.**
  - h(D) = (w/2π)√|D|·L(1, χ_D) for D < 0. This holds for non-fundamental D too, with χ_D the Kronecker symbol.
  - For D > 0: hR = ½√D·L(1, χ_D).
- **(F4) Generation** (Bach, Math. Comp. 55 (1990) 355–380; GRH). The classes of primes ≤ 6 log²|D| generate Cl(D)
  for fundamental D.

## 2. Route 1: 2-Sylow access

### 2.1 Lemma 0: what an ambiguous form gives (PROVED; elementary; unit-tested)
Let N > 12k² and q > 4k. Every reduced ambiguous form of D = −4kN is one of two kinds:
- it has first coefficient a | 4k (a *public* form; its class lies in A_pub, the subgroup generated by the classes of
  the prime ideals above 2k); or
- one of a, b, 2a ± b has gcd with N equal to p or q.

The class [𝔭] of the prime above p is not in A_pub. So A_pub ≠ Cl[2] and at least half of Cl[2] is useful.

*Proof.*
1. **Forms (a, 0, c).** Here ac = kN and a ≤ c, so a ≤ √(kN) < N. Hence gcd(a, N) ∈ {1, p, q}, and gcd = 1 means
   a | k.
2. **Forms (a, a, c).** Here a(4c − a) = 4kN and 4c − a ≥ 3a, so a ≤ √(4kN/3) < N. Again gcd(a, N) ∈ {1, p, q}, and
   gcd = 1 means a | 4k.
3. **Forms (a, b, a), 0 < b < a.** Put d₁ = 2a − b and d₂ = 2a + b. Then d₁d₂ = 4kN and d₁ ≤ d₂ ≤ 3d₁.
   - If gcd(d₁, N) = 1, then d₁ | 4k, so d₂ ≥ N > 12k ≥ 3d₁, which is impossible.
   - If N | d₁, then d₁d₂ ≥ N² > 4kN, also impossible.
   - So gcd(d₁, N) ∈ {p, q}.
4. **[𝔭] ∉ A_pub** (stated for the maximal order; for the non-maximal case read "ideal" as "form"). A relation
   [𝔭] = ∏[𝔩_i] over distinct ramified 𝔩_i above primes dividing 2k (these classes have order ≤ 2) would make the
   primitive ideal 𝔭∏𝔩_i principal, with norm p·d, d | 2k. Its generator x + y√−kN satisfies x² + kNy² = pd.
   - If y ≠ 0 this forces kN ≤ 2kp, i.e. q ≤ 2.
   - If y = 0 it forces x² = pd with p ∤ d.

   Both are impossible.
5. **Classes outside A_pub give p.** By F1 the reduced form of a class outside A_pub is ambiguous. If it were public
   (a | 4k), the class would be that of an ideal of norm a, a product of ramified primes above 2k, hence in A_pub.

∎ (Checked for k ∈ {1, 2, 3, 5, 6, 7, 10, 11} on six semiprimes; the count of reduced ambiguous forms matches
2^{2-rank} from PARI.)

### 2.2 Prop SQ: square roots in Cl(D) ⇔ factoring (KNOWN; precise version)
- **FACT ≤ SQRT** (randomized; GRH used only for generation).
  - Take D = −4kN fundamental: k = 1 if N ≡ 1 (mod 4), k = 2 otherwise. Let S be the split primes ≤ 6 log²|D|.
  - Sample x = ∏ f_ℓ^{e_ℓ} with e_ℓ uniform in [0, N²). Then x is within statistical distance |S|·h/N² of uniform on
    ⟨S⟩ = Cl (F4).
  - Give y = x² to an arbitrary root oracle, which returns some z with z² = y. Then t = z·x⁻¹ ∈ Cl[2].
  - Given y, x is uniform on a coset of Cl[2], and z depends on y only. So t is uniform on Cl[2], and
    Pr[t ∉ A_pub] ≥ ½ − |S|h/N².
  - Lemma 0 then gives p.
  - Experiment E3 used a deterministic oracle, blind to x: the lexicographically least root. Success rates were
    0.485 and 0.4875 per query (n = 24, 32; predicted ½), with 389/389 correct factors.
- **SQRT ≤ FACT.** With p and q known, membership in Cl² is decided by genus characters. A square root is then
  computed in polynomial time by Gauss's duplication/ternary-form method. References: Lagarias, "Worst-case complexity
  bounds for algorithms in the theory of integral quadratic forms", J. Algorithms 1 (1980); Simon, Math. Comp. 74
  (2005) — the exact statements are **UNVERIFIED** (literature critic).
- The equivalence is folklore in imaginary-quadratic cryptography (Buchmann–Hamdy survey, 2001; Biehl–Buchmann–
  Hamdy–Meyer, DCC 25 (2002)) — UNVERIFIED wording. **So "take a square root to land in Cl[2]" is factoring, not a
  route.**

### 2.3 Prop GEN: genus characters ≡ quadratic residuosity mod N (PROVED; randomized reductions; one direction uses
GRH for primes in progressions)
Let GEN-SQ be: given a form f of discriminant D = −4kN, decide whether [f] ∈ Cl² (the principal genus, F2). Let
QR_N be: given a with (a/N) = 1, decide whether a is a square mod N.

- **GEN-SQ ≤ QR_N.**
  - Evaluate f at random small (x, y) to get a represented m coprime to 2kN (a nontrivial gcd gives a factor).
  - The public assigned characters of m (2-adic characters and (m/ℓ) for ℓ | k) and (m/N) are computable.
  - [f] ∈ Cl² iff all of these are 1 and QR_N(m) holds, because (m/p) = (m/q) = 1 ⇔ (m/N) = 1 ∧ (m/p) = 1.
- **QR_N ≤ GEN-SQ** (`genus.qr_via_principal_genus`; unit-tested against a brute-force "is a square in Cl" oracle).
  1. Given a, use CRT to pick a prime m ≡ a (mod N) with m ≡ 1 (mod 8k), m < 8kN³.
  2. Then (D/m) = (−1/m)(k/m)(N/m) = 1·1·(m/N) = (a/N) = 1, by quadratic and Jacobi reciprocity using m ≡ 1 (mod 8)
     and m ≡ 1 (mod ℓ) for ℓ | k. So the prime form f_m exists and is computable (Tonelli–Shanks mod the prime m).
  3. All public characters of m equal their value at 1, namely 1. The hidden ones are (m/p) = (a/p) and
     (m/q) = (a/q) = (a/p).
  4. Hence [f_m] ∈ Cl² ⇔ (a/p) = 1 ⇔ a is a square mod N.
  5. Under GRH a random m in that progression below (8kN)³ is prime with probability ≍ 1/log N.

∎

- **Consequence.** "Genus characters are hidden" is exactly the QR_N assumption. It is not circular in the strict
  sense: QR_N is not known to be factoring-hard (Adleman–McCurley 1994, open problem O11; see
  `notes/codex_r4_sos_classgroup.md` §1).
- But no subexponential non-relational algorithm for QR_N is known either.
- Even a GEN-SQ oracle is not known to yield Cl[2] elements. Testing squareness does not produce a root.

### 2.4 Lemma GIRTH and Prop WALK: Cayley-graph statistics are relation statistics (PROVED)
Let S be a set of split primes ℓ ∤ D with prime forms f_ℓ (ideals 𝔩, conjugates 𝔩̄ = f_ℓ⁻¹). Let H = ⟨S⟩ and
Λ_S = {e ∈ Z^S : ∏ f_ℓ^{e_ℓ} = 1}.

**Lemma GIRTH.**
- (i) For e ≠ 0 in Λ_S, ∏ℓ^{|e_ℓ|} = x² + kNy² with gcd(x, y) = 1 and y ≠ 0. So ∏ℓ^{|e_ℓ|} ≥ kN and
  ‖e‖₁ ≥ g₀ := ⌈log(kN)/log ℓ_max⌉.
- (ii) Conversely, every primitive S-smooth value of x² + kNy² with y ≠ 0 gives an element of Λ_S.

*Proof.*
- 𝔞_e = ∏_{e>0}𝔩^{e}∏_{e<0}𝔩̄^{|e|} is primitive (for each ℓ only one of 𝔩, 𝔩̄ occurs) and has norm ∏ℓ^{|e|}.
- [𝔞_e] = 1 iff 𝔞_e = (α) with α = x + y√−kN ∈ Z[√−kN].
- If y = 0, then (x) is primitive only for |x| = 1, i.e. e = 0.
- For (ii): (x + y√−kN) is primitive of S-smooth norm, so it factors into 𝔩's and 𝔩̄'s.

∎

So the closed-walk count W_t = Σ_{λ∈Λ_S} #{walks in Z^S of length t ending at λ} is a function of Λ_S ∩ {‖·‖₁ ≤ t}.
In particular W_t = W_t^free for t < g₀. **Every return or collision observed by a walk statistic is a smooth value
of the principal form. That is the Hafner–McCurley/Seysen relation source, i.e. family (M).** E1 confirms
W_t = W_t^free below g₀ in 60/60 instances. The first excess appears at g₀ + {0..4}.

**Prop WALK (return-frequency statistics; PROVED).**
- The spectral expansion is P_t = W_t/(2|S|)^t = (1/|H|)Σ_{ψ∈Ĥ} λ_ψ^t, with λ_ψ = |S|⁻¹Σ_ℓ Re ψ(f_ℓ).
- The hidden genus characters restrict to χ_p(f_ℓ) = (ℓ/p), so λ_{χ_p} = |S|⁻¹Σ_{ℓ∈S}(ℓ/p). This is a QR-type
  aggregate.
- Let G_t be the contribution of the hidden characters that are distinct on H from the public ones (deviation D1,
  §7). For even t, every term is ≥ 0, so b := P_t − G_t ≥ 1/|H| and G_t ≤ 2/|H|.
- Bernoulli(P_t) and Bernoulli(b) have squared Hellinger distance ≤ (P_t − b)²/(P_t + b) ≤ 2/|H|. M independent
  walks therefore distinguish "with" from "without" the hidden term with advantage ≤ √(4M/|H|).
- Advantage δ needs M ≥ δ²|H|/4 ≈ δ²·N^{1/2}.
- For collision counts (pairs of walks), the signal-to-noise heuristic gives M ≳ √(|H|/2). This is the birthday
  bound, consistent with Theorem GB.
- E1, for the variance criterion M = P_t/G_t²: the minimum over t of M/|H| was ≥ 1.0 (as-registered) and ≥ 2.0
  (D1) in every instance. Medians were 1.4·10³, 4.1·10³, 1.3·10⁴ at n = 20, 24, 28. The predicted upper bound of 100
  (P2) was **not met**: in the median instance, detection is far *harder* than the proved floor.

### 2.5 4-rank via Rédei
- Rédei matrices are indexed by the prime-discriminant factors (·/p) and (·/q). Building them needs the factors.
- Without the factors, the 4-rank still decides a fixed-element QR bit (codex_r4 F2). For N ≡ 1 (mod 8):
  [2 is a square mod N] ⇔ r₄(Cl(−8N)) > 0 if p ≡ q ≡ 1 (mod 4), and ⇔ r₄(Cl(−4N)) > 0 if p ≡ q ≡ 3 (mod 4).
  - We re-verified this independently: 92/92 random pairs p, q < 400, using PARI `quadclassunit`
    (`scratchpad/genusleak/r4check.py`).
- So computing r₄ is at least as hard as that QR instance. It is also not a route to Cl[2]: r₄ describes Cl[4]/Cl[2],
  and Cl[2]'s *size* is already public.

### 2.6 Route-1 verdict
Every way into Cl[2] that we found is one of three things:
- a root (factoring-equivalent);
- a group-order computation (an odd multiple kills the odd part, then square down: Shanks; Schnorr–Lenstra);
- a relation. Walks and spectra are relations by GIRTH.

Genus characters are QR_N-hard and do not by themselves produce Cl[2] elements.

## 3. Route 2: infrastructure (D = 4kN)

### 3.1 Midpoint (KNOWN; E4)
- When N(ε) = +1, the principal cycle is symmetric under conjugation. The reduced principal ideal at distance R/2 is
  ambiguous, with a | 4kN, and a ∤ 4k gives p.
- In the continued-fraction picture this is the symmetry point of the even period of √(kN) (Legendre; Lehmer–Powers
  1931).
- E4: at k = 1 the midpoint splits N in 72% (n = 28) and 64% (n = 32) of instances; some k ≤ 11 splits 97% / 98%.

### 3.2 Prop HR: any multiple of R suffices (PROVED, given standard infrastructure arithmetic)
- **Statement.** Suppose X = M·R is known to absolute precision 0.1, for an unknown integer M ≥ 1. Then the
  ambiguous ideal at R/2 is reached in poly(log X) operations.
- **Proof.**
  - For j = 1, 2, …, ⌊log₂ X⌋ + 1, compute the reduced principal ideal nearest distance X/2^j. Shanks–Lenstra
    giant steps track distances to O(log X) bits: Lenstra 1982; Buchmann–Williams; Jacobson–Williams, *Solving the
    Pell Equation* (2009), Ch. 10–11.
  - Inspect its O(1) neighbours for an ambiguous form.
  - At j = v₂(M) + 1 the distance X/2^j = (odd)·R/2 ≡ R/2 (mod R).

  ∎
- **Consequence.** Separating h from R is unnecessary: hR itself (F3) is such an X. Exact hR ⇒ factor with constant
  probability per k (E4).
- This is the classical "regulator ⇒ factoring" reduction (Lenstra 1982; Schoof 1982; used by Hallgren, JACM 54
  (2007)). Our 2-adic sweep is only a convenience; no novelty claimed.

### 3.3 What a precision-Δ approximation of hR buys (PROVED as an upper bound; KNOWN in substance)
- Given |A − hR| ≤ Δ, find a distance in [A − Δ, A + Δ] at which the principal ideal recurs. Any such distance is a
  multiple of R, which is enough by Prop HR.
- BSGS (baby steps of total distance w around 0, giant steps of length w) finds one in O(w + min(Δ, R)/w) steps,
  i.e. O(√min(Δ, R)).
- The imaginary case is the same with R := 1: BSGS for an element's order in the window, then an odd multiple, then
  square down.
- **Euler products.**
  - Under GRH the truncated product to x has |log L(1,χ) − log L_x| ≪ x^{−1/2}log²(x|D|) (Bach, Math. Comp. 1995,
    "Improved approximations for Euler products" — exact constant UNVERIFIED). So Δ(x) ≈ √|D|·x^{−1/2}.
  - The total cost x + |D|^{1/4}x^{−1/4} is minimised at x = |D|^{1/5}. That is the Lenstra/Schoof exponent (also
    `notes/analytic_results.md`, `notes/below16_results.md` §5).
  - E5 measures the fitted slope d log Δ/d log x ∈ {−0.52, −0.53, −0.56} and log C_min/log|D| ∈ {0.20, 0.18, 0.18}
    at n = 32, 40, 48.
- **Exponent 1/5 is the floor of this architecture.**
  - A subexponential factoring algorithm along this route needs Δ ≤ L[½] in subexponential time, i.e. relative
    precision |D|^{−1/2+o(1)} for L(1, χ).
  - Euler products need x ≈ |D|^{1−o(1)} for that: tail variance ≍ 1/(x log x) (below16 §5.3; heuristic).
  - We know no other non-relational subexponential approximation of L(1, χ_D) to that precision.

### 3.4 Square forms and giant steps (OBSERVATION + E2)
- A proper square Q_i = r² in the principal cycle is 𝔯² with N𝔯 = r ≲ D^{1/4} and [𝔯] ∈ Cl⁺[2] ∖ {1}.
- So SQUFOF finds 2-torsion through a *multiplicative self-coincidence* of a reduced norm (a = r·r). The density of
  such points per cycle position is ≍ |Cl⁺[2]|·D^{−1/4} (heuristic; Gower–Wagstaff 2008).
- Giant steps reach chosen distances. The distance of 𝔯² is 2·dist(𝔯), where dist(𝔯) is measured in the cycle of a
  hidden class, so giant steps cannot aim at squares.
- In the generic infrastructure model below, jumping to random distances costs 1/ρ forms, exactly like walking.
- **E2.**
  - Median i₁·ρ (first proper square × density): 0.45, 0.54, 1.00 at n = 28, 32, 36. P4 (∈ [0.5, 2]) is **met at
    32 and 36, not met at 28**: sequential walking is, if anything, slightly better than random jumps.
  - Fitted slope of log₂ median i₁: 0.38 (P5, ∈ [0.18, 0.32], **not met**).
  - Post hoc: the slope of log₂ median(1/ρ) is 0.24, consistent with N^{1/4}. The excess in i₁ at small n comes from
    i₁·ρ drifting upward with n.
- CFRAC is the relational completion of SQUFOF: products of Q_i that form squares.

### 3.5 Are Buchmann/Vollmer subexponential regulator algorithms relational? (KNOWN; references UNVERIFIED in detail)
Yes. They all collect relations among prime ideals of smooth norm, then do lattice reduction / HNF (Hafner–McCurley
style):
- Buchmann, "A subexponential algorithm for the determination of class groups and regulators of algebraic number
  fields", Sém. Théorie des Nombres Paris 1988–89 (1990);
- Vollmer, ANTS-IV (LNCS 1838, 2000) and ANTS-V (LNCS 2369, 2002);
- Biasse–Jacobson, ANTS-IX (2010).

Smooth norms put them in (M).

### 3.6 Prop GI: generic infrastructure barrier (PROOF SKETCH; same argument as Theorem GB)
- **Model.** The algorithm holds handles of reduced principal ideals with known approximate distances. It may
  compose them (distances add, error ≤ c) and compare handles. The regulator R is hidden, with density ≤ 1/X on
  [X, 2X], or ≤ 1/Δ on an analytic window of width Δ.
- **Collisions.** A collision means some known difference d = δ₁ − δ₂ (|d| ≤ 2^L) is within c of mR. The set of R
  for which a fixed d does this has measure ≤ Σ_{m ≤ |d|/X} 2c/m ≤ 2c(1 + L ln 2).
- **Bound.** The rest is as in Theorem GB: Pr[success] ≤ (q + 2)²·2c(1 + L)/X with L ≤ q + O(1), so q = Ω(X^{1/3}).
  BSGS achieves √X.
- **Reading.** Without relations and without square coincidences, the infrastructure is a hidden-period problem. R
  is the order of R/RZ.

## 4. Route 3: isogeny volcanoes / CM (failure points)

### 4.1 Construction
- An ordinary E/F_ℓ with Z[π] of discriminant t² − 4ℓ = −4kN f² needs ℓ ≥ kNf² (4ℓ = t² + 4kNf²).
- The curves with trace t number H(4ℓ − t²) ≈ f√(kN)·polylog (Deuring; Kronecker–Hurwitz class number), out of
  ≈ 2ℓ curves in total. So random curves plus SEA need ≈ 2ℓ/H ≥ √(kN)/polylog trials, which is worse than trial
  division.
- Fixing ℓ does not enlarge the target. At most 6 traces t satisfy (t/2)² + kNf² = ℓ.
- The CM method needs the class polynomial H_O mod ℓ.
  - deg H_O = h(O) ≥ h(O_K)/3 = |d_K|^{1/2+o(1)} by Siegel (ineffective), with K = Q(√−kN).
  - Computing it costs Õ(|D|): Enge, Math. Comp. 78 (2009); Sutherland, Math. Comp. 80 (2011).
- We know of no method that produces any curve whose endomorphism algebra is K in time |d_K|^{o(1)}. Producing one
  is the well-known "CM bottleneck"; see e.g. Bröker–Stevenhagen 2007 on prime-order curves, where |D| must be small.
- **Failure point:** construction is exponential.

### 4.2 Even given the curve
- **Crater walks.** Small split ℓ′-isogenies realise the class-group *action*: a hidden-order group action, weaker
  than a black-box group. The crater length is ord([𝔩′]), an order computation.
- **The ambiguous classes.** [𝔭] acts by an isogeny of degree p, with kernel E[p] ∩ ker(π − t/2), since
  (π − t/2)² = −kNf². Computing it needs p.
- **Choosing t ≡ 2 (mod 2N).** This makes ker(√−kN) ∩ E[N] F_ℓ-rational: E(F_ℓ) contains Z/N. But a point of order N
  does not separate p from q.
- **Group structure.** E(F_ℓ) ≅ O/(π − 1) as O-modules (H. W. Lenstra Jr., J. Number Theory 56 (1996)). Its order
  #E = (t/2 − 1)² + kN f² is public.
- **No leak found.**

## 5. Route 4: non-maximal orders and multiplier pencils

### 5.1 Ring class groups
- 1 → K_m → Cl(Dm²) → Cl(D) → 1, with |K_m| = m∏_{ℓ|m}(1 − (D/ℓ)/ℓ)/[O*: O_m*] (Cox Thm 7.24). The order of K_m is
  public.
- A known multiple of |K_m| kills K_m, not the hidden odd part of Cl(D).
- In Theorem GB the Z/P-component survives unless P | m∏(…), i.e. unless the algorithm "guessed" a multiple of P.
  That is an L-bit integer divisible by P, which has probability ≤ Lρ/log X.
- No gain.

### 5.2 Multiplier pencil
- The hidden data common to all Cl(−4kN) is the character (·/p), not a group element. The orders h(−4kN) are
  *public* Euler products in k (F3), and they vary independently enough to be separate lottery tickets: Schnorr–
  Lenstra's multiplier trick, i.e. (M).
- **The compositum L = Q(√−kN, √−k′N) ⊃ Q(√kk′).**
  - Since p ∤ kk′, p is unramified in Q(√kk′). So 𝔭_k O_L = 𝔭_{k′}O_L (both are the product of the primes of L above
    p), and [𝔭_k], [𝔭_{k′}] extend to the same class.
  - Finding such a pair is a cross-group collision in Cl(L). There, equality testing is a principal-ideal problem in a
    quartic field (subexponential only relationally), over a group of order ≈ h_k h_{k′}h(kk′) (Kuroda's class
    number formula; constant UNVERIFIED).
- Many discriminants give more (M)-type chances, not a non-smooth genus character.

## 6. Route 5: the barrier, and the non-generic escapes

### 6.1 Theorem GB (generic hidden-order barrier; PROVED in the model)
**Model.**
- G = G′ ⊕ Z/P is finite abelian. P is an odd prime drawn from a law π with min supp π ≥ X and max point mass ρ.
  P ∤ |G′|.
- A ≤ G′[2] is public, given with a basis.
- Generators g_i = (g_i′, u_i), i ≤ s. The u_i are uniform and independent in Z/P, and g_i′ is arbitrary — the
  algorithm may even know G′ completely.
- Random injective encoding σ: G → Σ. The oracle returns σ(x − y). The inputs are σ(0), σ(g_i), and σ(a_j) for
  j ≤ r.
- Success: output σ(t) with t ∈ G[2] ∖ A.

**Theorem GB.** Any algorithm making q oracle calls succeeds with probability

  ≤ ((q + s + r + 1)² + 1)·(1/X + ρ(q + 1)/log₂X) + (q + 1)|G|/|Σ|.

For P uniform on the primes of [X, 2X] this is O(q³/X). So q = Ω(X^{1/3}).

*Proof (Shoup-style simulation; the order is unknown, so divisibility by P is bounded by counting).*
1. **Formal vectors.** Track, for every handle, a formal vector (e, c) ∈ Z^s × F₂^r. Inputs are unit vectors and
   oracle answers are differences, so |e| ≤ 2^{q+1}.
2. **The simulation.** Simulate with fresh random strings, identifying two handles only when their formal vectors
   are equal. Since A has a known basis, e = e′ and c ≠ c′ really do give different elements.
3. **When simulation and reality differ.** They can diverge only if some pair with d = e − e′ ≠ 0 has d·u ≡ 0
   (mod P).
   - If d ≢ 0 (mod P), this has probability 1/P ≤ 1/X over u.
   - If d ≡ 0 (mod P), then P divides a nonzero integer of absolute value ≤ 2^{q+2}. That integer has
     ≤ (q + 2)/log₂X prime factors ≥ X, so the probability over P is ≤ ρ(q + 2)/log₂X.
   - In the simulation the pairs are independent of (u, P), so a union bound over pairs applies.
4. **Success in the simulation.** The output is either a fresh string (probability ≤ |G|/|Σ|) or a handle (e, c).
   For t ∈ G[2] ∖ A we need e ≠ 0 (otherwise t = c·a ∈ A) and e·u ≡ 0 (mod P), because P is odd and t has order
   ≤ 2. That has the same probability bound.

∎

**Corollary GB-X (extraction; any G, no distribution).**
- Suppose a generic algorithm outputs a known combination t = e·g + c·a ∈ G[2] ∖ A. Then e·π_odd(g) = 0 and
  e·g ≠ c·a, so 2e is a nonzero relation of G.
- For s = 1, e is a nonzero multiple of the odd order of g.
- Conversely, an odd multiple M′ of ord_odd(g) gives g^{M′} ∈ G₂; squaring down reaches G[2].

So in rank 1, "hidden 2-torsion in ⟨g⟩" and "a multiple of ord_odd(g)" are interchangeable. In general rank, success
hands over a relation. **In this model, family (Q) is a group-order / relation problem.** Several multipliers k give
several groups; the union bound adds their success probabilities. A small P⁺(h_odd) for some k is exactly the
Schnorr–Lenstra smoothness event, i.e. (M).

**Scope.**
- The bound O(q³/X) is weaker than the √X of BSGS. The q³ comes from exponents divisible by many primes in [X, 2X].
  We conjecture q = Θ̃(√X).
- Generic lower bounds for root extraction in groups of hidden order: Damgård–Koprowski, EUROCRYPT 2002. Related:
  Rotem–Segev–Shahaf, EUROCRYPT 2020; the low-order assumption, Boneh–Bünz–Fisch, CRYPTO 2019; Belabas–Kleinjung–
  Sanso–Wesolowski, ePrint 2020/1310. Exact statements UNVERIFIED (literature critic).
- Real class groups are not drawn from π. Fitting the model to them is a Cohen–Lenstra-type heuristic, and the model
  is broken by the non-generic escapes below.

### 6.2 Which non-generic features escape, and to where

| Feature of class groups | In the generic model? | Known escape | Lands in |
|---|---|---|---|
| Canonical reduced forms (equality) | yes (encoding) | — | generic |
| Integer norm a of a reduced form: smoothness test | no | relations (Seysen 1987; Hafner–McCurley 1989; Buchmann; Vollmer; Biasse); L[½] under GRH | (M) |
| Norm a: squareness test | no | SQUFOF (2-torsion shadows, density D^{−1/4}) | (Q), N^{1/4} |
| Norm a: equality between forms / factorisation over earlier norms | no | collisions / large-prime variation | generic / (M) |
| gcd(a, N) on random forms | no | hits with probability ≍ 1/p | accidental |
| Explicit forms of prescribed prime norm | no | sampling ideals; QR_N ≤ GEN-SQ (Prop GEN) | QR_N |
| Genus characters | no | ≡ QR_N | QR_N (open) |
| Explicit list of ambiguous forms (= divisors of D) | no | divisor search | (A) / trial division |
| Distances in the infrastructure | partly (GI) | R = period; any multiple of R ⇒ factor (HR) | group order of R/RZ |
| Analytic class number formula | no | hR to precision Δ ⇒ √Δ search; |D|^{1/5} (GRH) | order approximation |
| CM curves / isogenies | as a group action | construction exponential | — |
| Conductors / multipliers / composita | yes (more groups) | more lottery tickets | (M) |

**Open class.** Algorithms that use the integer coefficients of reduced forms (or CF partial quotients, or
distances) in ways that are not multiplicative coincidences (smooth, square, equal), not divisor searches, and not
order approximations.
- A pure-integer statistic of the (a_i) along the principal cycle is an example. eps_results §5 found such statistics
  correlating with ε, but only through ((−1)^iQ_i | p) admissibility, i.e. the QR_N-type genus data of Prop GEN, not
  through Cl[2] elements.
- We have no barrier for this class. A barrier there would be a factoring lower bound.

## 7. Experiments (pre-registered; outcomes)

| P | Statement | Outcome |
|---|---|---|
| P1 | W_t = W_t^free for t < g₀ | **met** (60/60) |
| P2 | min M_ret/\|H\| ≥ ¼ (proved) and median ≤ 100 | lower bound **met** (min 1.0); upper **not met** (medians 1.4e3, 4.1e3, 1.3e4) |
| P3 | min M_coll/√\|H\| ≥ 1/√2 (proved) and median ≤ 50 | lower bound **met** (min 1.41); upper **not met** (53, 90, 160) |
| P4 | median i₁ρ ∈ [0.5, 2] | met at n = 32 (0.54) and 36 (1.00); **not met** at 28 (0.45) |
| P5 | slope of log₂ median i₁ ∈ [0.18, 0.32] | **not met** (0.38); post hoc, the slope of 1/ρ is 0.24 |
| P6 | sqrt-oracle success 0.5 ± 0.08 | **met** (0.485, 0.4875; H = Cl in 20/20) |
| P7 | claimed factors correct | **met** (389/389) |
| P8 | midpoint splits at k = 1 in [0.2, 0.8] | **met** (0.72, 0.64) |
| P9 | some k ≤ 11 splits ≥ 0.9 | **met** (0.97, 0.98) |
| P6′ | Euler-product error slope ∈ [−0.75, −0.25] | **met** (−0.52, −0.53, −0.56) |
| P7′ | log C_min / log\|D\| ∈ [0.15, 0.30] | **met** (0.20, 0.18, 0.18) |

**Deviation D1 (disclosed).** The registered G_t double-counts when a hidden character coincides on S with a public
one: χ_p|_S trivial (the H = h/2 cases) or χ_p|_S = χ₋₄|_S. In those cases the "hidden" term is the public one. We
report both versions. The proved lower bounds hold for both. No instance changed status. Instances with no distinct
hidden character on H: 1, 1, 2 at n = 20, 24, 28.

**Reproduction.** `.venv/bin/python -m src.genusleak.harness E1` (likewise E2–E5); `.venv/bin/python -m pytest -q
tests/test_genusleak.py`. Wall times: E1 6 s, E2 18 s, E3 9 s, E4 2 s, E5 10 s.

## 8. Candidate cards (rule 9)
- **GL-1 Cayley spectral detection.**
  - Target: λ_{χ_p} = mean (ℓ/p).
  - Missing lemma: detecting it below |H| samples.
  - Failure: Lemma GIRTH + Prop WALK give ≥ δ²|H|/4 walks for returns, and returns are relations.
  - Scope: refuted for this model (return/collision statistics).
- **GL-2 Infrastructure jump to R/2.**
  - Works given any multiple of R (Prop HR).
  - Failure: obtaining that multiple is the hidden-period problem (Prop GI), or relations, or an analytic window
    (|D|^{1/5}).
- **GL-3 Square roots / genus characters.** Factoring-equivalent (SQ) / QR_N-equivalent (GEN). No route.
- **GL-4 CM volcano.**
  - Failure: curve construction (√(kN) trials, or a degree N^{1/2+o(1)} class polynomial).
  - Also: the ambiguous isogeny has the unknown degree p.
- **GL-5 Multiplier pencil / compositum.** Only more hidden-order groups (GB union bound). Cross-group collisions
  live in Cl(L), which is larger.

## 9. Most decisive open question
Is there a polynomial-time computable statistic of the reduced forms in the principal cycle of √(kN) (or of the
reduced forms reached by composition), *not* of the form "multiplicative coincidence among norms", whose correlation
with a hidden 2-torsion class (rather than with the QR_N-type characters (ℓ/p)) is non-negligible?

By Prop GEN, correlation with (ℓ/p) alone buys at most QR_N information. By Theorem GB, any generic construction of
the class itself is a relation.

## 10. Reviews
(See below; filled in after the codex rounds and the Opus critics.)
