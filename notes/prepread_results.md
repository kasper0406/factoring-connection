# prepread: high-degree algebraic preprocessing, characters, and dense short SLPs (phase 18)

Agent `prepread`, 2026-10-20.
- Pre-registration: `notes/prepread_prereg.md` (commit `fdb1328`). Amendment A1 (E4b cost criterion) was
  committed before any experiment ran.
- Code: `src/prepread/{slp,spectra,dlred,compsearch,run_experiments,analyze}.py` and `src/prepread/compsearch.c`.
  Tests: `tests/test_prepread.py`.
- Raw data: `experiments/raw/prepread_e{1,1b,2,3,4,5,5b}.jsonl`. E1b and E5b are post-hoc.
- GPT-6.1-sol:
  - blind round `codex/prepread_1.{md,out}` (full report in `codex/scratch/prepread_1/report.md`,
    git-ignored);
  - adversarial round `codex/prepread_2.{md,out}`, verdict MAJOR REVISION. All fixes are applied (§9).
- Two read-only Opus critics, one on proofs and one on literature. Adjudication in §9.
- Labels:
  - **KNOWN** (with reference);
  - **PROVED** (proof here; "no prior formulation found by us" is the strongest novelty claim);
  - **CONDITIONAL**;
  - **HEURISTIC**;
  - **CONJECTURED**;
  - **experimentally supported**;
  - **REFUTED for this model**;
  - **UNVERIFIED**.

**Notation.**
- N = pq is balanced, with p, q ∈ (X, 2X]. e(x) = exp(2πix). π*(X) = π(2X) − π(X).
- For Φ ∈ Z[y] (or an SLP), a prime p, a multiplier c ∈ F_p and a multiplicative character χ of F_p^*:

  A_p(χ, c) = (1/(p−1)) Σ_{y∈F_p^*} χ(y)·e(cΦ(y)/p).

- For the composite observable Ψ(y) = e(hΦ(y)/N), CRT gives e(hΦ(y)/N) = e(hq̄Φ(y)/p)·e(hp̄Φ(y)/q), with
  q̄ = q^{−1} mod p and p̄ = p^{−1} mod q. Hence Ψ̂(χ_pχ_q) = A_p(χ̄_p, hq̄)·A_q(χ̄_q, hp̄).

---

## 0. Verdict (plain language)

**No mechanism. The two remainders behave very differently, and the difference tracks the order of the character
that the readout needs.** This is a framing, not a proved complexity dichotomy.

**O1 (high-degree preprocessing, additive readout).** A working O1 observable would also be a prime-field
discrete-log algorithm.
- Theorem O1-DL (§2.4, PROVED; routine Fourier/list-decoding methodology): an efficiently computable phase that
  correlates with a character of order d at a known prime p yields discrete logarithms modulo d.
- Corollary O1-N: an efficiently constructible O1 observable with the heavy large-order peaks that periodfind's
  Theorem T requires would give polynomial-time discrete logarithms in F_p^* for at least a μ/2 fraction of the
  primes in the range, on top of factoring.
- Prime-field DL is not known to follow from factoring. So an O1 factoring method would have to solve a second
  problem that is believed hard (L_p[1/3] by the best known methods). This is a CONDITIONAL barrier. It is not a
  separation.
- Unconditionally, Weil bounds the correlation only for cyclic exponent width below √p (§2.1). That needs SLP size
  only ½·log₂ p, so the unconditional bound is vacuous for short SLPs.
- A barrier from "small SLP plus large character order" alone is false. Pohlig–Hellman is a short SLP (§2.3). An
  N-independent two-input program of size O(log m + e²) has correlation ≥ 0.95 with a character of order 2^e at
  every prime with p − 1 = 2^e·u, u | m odd, u ≥ 3, for a p-dependent multiplier c. That is the p − 1
  (order-annihilation) regime.
- O1-DL needs the order to contain a prime ℓ whose DL is hard. With near-full order the relevant prime is
  P⁺(p − 1).

**The dense-torus question has no per-prime obstruction.** For unrestricted prime-specific polynomials, root
density 1/2 at a known prime is trivial (x^{(p−1)/2} − 1). The difficulty is uniformity in p.
- **Galois theory alone cannot give a barrier** (Prop GAL-0, PROVED).
  - On every interval there are polynomials with trivial, cyclic and symmetric Galois groups whose root density
    lies in (1/6, 1/3] at every prime of the interval.
  - The fixed-point proportion of the group is the wrong statistic (§3.1).
  - Under GRH, effective Chebotarev does exclude dense roots, but only for polynomials of degree ≲ X^{3/2} with few
    irreducible factors, small discriminants and few bad primes in the interval (Prop CS, §3.2, corrected after
    review). Short SLPs exceed all of these at size ≈ 1.5·log₂ X.
  - Abelian versus non-abelian is the wrong dividing line. ECM's division-polynomial fields have GL₂ images
    (§3.4).
- **Galois theory plus the arithmetic of the roots gives a restricted inverse theorem.** Theorem MIN (§4.4,
  PROVED; no prior formulation found by us; checked by codex and the proof critic) says:
  - take compositions of monic integer polynomials from a fixed finite set (e.g. x² + c with |c| ≤ 2), of degree at
    most 2^Δ·X;
  - then the proportion of primes p ∈ (X, 2X] at which root density ≥ μ tends to 0, uniformly over such
    compositions.
  - The proof combines Chebotarev on the finitely many outer prefixes, Burnside, Markov and Northcott. Roots of
    such compositions have bounded house.
  - So a dense composition needs excess depth Δ → ∞.
  - Polylogarithmic depth stays open (Conjecture QC, §4.6). Codex's fixed-prefix and outer-constant obstructions
    (§4.2–4.3) cover further subclasses.
- **Experiments agree.**
  - Exhaustive search over all 5^k quadratic compositions (E2): no word reaches a minimum fibre above k across 8
    primes.
  - Iterates of x² + 1 follow the Haar wreath-product model even far beyond the Chebotarev range (E3).
  - No high-degree random SLP correlates with large-prime-order characters beyond 2× the random baseline, except
    a monomial at a prime with smooth p − 1, which is order annihilation (E1, E1b).

**Where the limit bites** (§5):
- For O1, at DL in the part of F_p^* detected by the character, i.e. at P⁺(p − 1).
- For the dense-torus question, nowhere per prime; only through uniformity in p. There slpcirc's Prop BB makes a
  general barrier permanent-hard, and GAL-0 rules out a group-theoretic shortcut.

**Subjective assessment.** Probability that either remainder yields sub-L[1/3] factoring: < 1%.

---

## 1. The common frame: one two-parameter family of sums

**Prop UF (the two faces). PROVED, trivial.** For Φ ∈ Z[y] and a prime p,

  #{y ∈ F_p^* : Φ(y) ≡ 0}/(p − 1) = E_{c∈F_p} A_p(χ₀, c),   χ₀ = trivial character.

*Proof.* E_{c∈F_p} e(cv/p) = [v ≡ 0]. ∎

So both remainders ask for large values of the same family A_p(χ, c).
- **slpcirc (model U, W_U).** The trivial character, averaged over c, must lie in [μ, 1 − μ] at almost every
  prime. The readout is a gcd.
- **periodfind O1.** A single (χ, c) needs |A_p(χ, c)| ≥ τ, with χ of order ≥ p/poly, at both p and q
  (Corollary O1-N). The readout is a Fourier peak, so only large-order characters help (Theorem T).
- **Small orders.** They reveal only [d | r]-type residuosity bits (periodfind §1.2, §2).

The character order distinguishes the observables formally. The computational difficulty depends on the order's
factorisation, on access to the multiplier and frequency, and on uniformity (codex).

---

## 2. O1: high-degree preprocessing with an additive readout

### 2.1 The Weil threshold. KNOWN (Weil; Castro–Moreno; Fu–Wan) + PROVED (width reduction)

*Bound.* Let χ be nontrivial and let f ∈ F_p[y, y^{−1}] be a nonconstant Laurent polynomial with exponents in
[−a, b], a, b ≥ 0, a, b < p. Then

  |Σ_{y∈F_p^*} χ(y)·e(f(y)/p)| ≤ (a + b)√p.

- This is Weil's bound for L_χ ⊗ L_ψ(f) on G_m. The Swan conductors are at most a at 0 and b at ∞, and H²_c = 0.
- References: Castro–Moreno, Proc. AMS 128 (2000) 2529–2537. The conductor form for rational f is Fu–Wan,
  Quart. J. Math. 65 (2014) 1195–1211, Theorem 4.6 of arXiv:1303.3650v3 (checked by codex against v3).
- Brute force (proof critic): p ≤ 127, a, b ≤ 6, all nontrivial χ. The maximum ratio to the bound is 1.0000,
  attained by Gauss sums.

*Cyclic width (PROVED).* Only Φ mod (y^{p−1} − 1) matters. Lift its nonzero exponents to integers in a window
[−a, b] that minimises a + b. Call the minimum the cyclic width w_p(Φ). It is the length of the shortest arc of
Z/(p − 1) containing the support and 0 (`slp.weil_width`). Then

  |A_p(χ, c)| ≤ w_p(Φ)·√p/(p − 1)  for every nontrivial χ and every c ≠ 0.

**Consequences.**
- **Threshold.** The bound is trivial exactly when w_p ≥ (p − 1)/√p, i.e. **w ≈ √p**. A τ-correlation needs
  w_p ≥ τ(p − 1)/√p.
- **Composite.** For χ_p, χ_q nontrivial, |Ψ̂(χ)| ≤ w_p w_q √N/((p − 1)(q − 1)) ≈ w_p w_q/√N.
- **Short SLPs.** w_p ≤ deg Φ ≤ 2^s, so the bound binds only when s < ½·log₂ p − log₂(1/τ). It gives no typical
  value of w_p. In E1b, 131 of the 153 member-prime pairs whose width was computed had w_p ≥ √p (§6).
- **Conductor.** Conductor, not degree, is the right complexity (codex). y^{p−2} = 1/y has width 1. Finding a
  low-conductor representative need not preserve SLP size.
- **Compositions (HEURISTIC remark).** For even χ (χ(−1) = 1; for odd χ the sum vanishes by x ↔ −x),
  Σ_x χ(x)e(G(x² + c)/p) = Σ_u [λ(u − c) + λη(u − c)]e(G(u)/p) with λ² = χ. A k-layer composition is therefore
  a sum against a trace function of rank 2^k. Its conductor grows like the degree, so a sheaf bound is again
  vacuous above √p.

### 2.2 Monomials. PROVED from a KNOWN theorem (codex)

Let d = gcd(m, p − 1), H = (F_p^*)^m and L = |H|.
- If χ is nontrivial on ker(y ↦ y^m), then A_p = 0.
- Otherwise |A_p(χ, c)| = (d/(p − 1))·|Σ_{z∈H} θ(z)e(cz/p)|, with ord χ ≤ L. Expanding the square and applying
  Bourgain–Glibichuk–Konyagin (J. London Math. Soc. 73 (2006) 380–398, Thm 6) gives

    |A_p(χ, c)| ≪ (L^{−1} + p^{−η(ε)})^{1/2}  for ord χ ≥ p^ε.

  A checked statement of the input with an explicit implied constant is Kowalski's notes on sums over small
  subgroups, Thm 1.1 (codex).

So monomial preprocessing of any degree has power cancellation against characters of polynomially large order.
Concentration at small L is order annihilation (periodfind §1.3f). E1's largest non-control ratio, 6.6, is
exactly this case: x^{lcm(1..30)} at p = 13841, where p − 1 = 2⁴·5·173 and the image has 173 elements.

### 2.3 Prop PH (Pohlig–Hellman is a short SLP). PROVED

**The program.**
- Fix e ≥ 1 and an odd integer m. The program has inputs y and w (w is an auxiliary random input) and uses ring
  operations and halving. Halving is division by the unit 2, gcd-checked as in slpcirc §1.4.
- Set z = y^m and ζ = w^m.
- For i = 0..e − 1: s_i = z_i^{2^{e−1−i}} ∈ {±1}, b_i = (1 − s_i)/2, and
  z_{i+1} = z_i·((1 + s_i) + (1 − s_i)ζ_i^{−1})/2, with ζ_0^{−1} = ζ^{2^e−1} and ζ_{i+1}^{−1} = (ζ_i^{−1})².
- Output L = Σ b_i 2^i. Size O(log m + e²). Every exponent and constant is a fixed integer.

**Claim.**
- Let p − 1 = 2^e·u with u odd and u | m. Note that u | m is stronger than u being B-smooth when m is the odd part
  of lcm(1..B), because prime powers matter.
- Let w be a quadratic non-residue (probability 1/2). Then ζ generates μ_{2^e}, and L(y) = log_ζ(y^m) ∈ [0, 2^e).
- Let χ(y) = e(L(y)/2^e), a character of order 2^e. Take c = round(p/2^e), which depends on p. Then
  |A_p(χ̄, c)| = |sin(π2^eδ)/(2^e sin πδ)| with |δ| = |c/p − 2^{−e}| ≤ 1/(2p). This is ≥ sinc(π/(2u)) ≥ 0.95 for
  u ≥ 3.

*Proof.*
- y ↦ y^m maps F_p^* onto μ_{2^e} because u | m.
- ζ generates because w^{m2^{e−1}} = (w/p)^{m/u} = −1.
- Peeling works because s_i = (−1)^{b_i} when z_i ∈ ⟨ζ^{2^i}⟩.
- L is uniform on [0, 2^e), so A_p(χ̄, c) is a geometric sum. ∎

Checks:
- unit test at p = 3·2^{12} + 1: correlation 0.99 with 274 operations;
- proof critic: 37 primes p = u·2^e + 1 with e = 4..13, minimum correlation 0.994;
- sign convention: codex's counterexample p = 13 shows that it is χ̄, not χ, that correlates.

**Reading.**
- Codex's prime-specific version uses y^{(p−1)/R}, a field constant ζ and a Boolean rounding circuit. It gives
  correlation ≥ cos(π/p) with order R = 2^b, and R ≥ p^{0.19} infinitely often (Xylouris).
- So an unconditional barrier "SLP size s and order ≥ R ⇒ small correlation" is false, even for programs whose
  arithmetic is N-independent.
- PH supplies examples. It does not classify all exceptions.
- The next theorem says what any example must be: a discrete-log algorithm for the part of the group that the
  character detects.

### 2.4 Theorem O1-DL (correlation ⇒ discrete logarithm). PROVED (routine)

This is established Fourier/list-decoding methodology applied to our observable (codex), the classical analogue of
the Fourier step in Shor's DL algorithm.

**Statement.**
- *Setting.* Let p be prime, n = p − 1, and f: F_p^* → {|z| ≤ 1} computable in time T. Suppose |⟨f, χ⟩| ≥ τ for
  some character χ.
- *Output.* Given p, a generator g of F_p^* and a target y, the algorithm outputs a list of size O(τ^{−2}) of pairs
  (j′, x′ mod n/gcd(j′, n)). The list contains the pair for the frequency j of χ (χ(g) = e(j/n)) with
  probability ≥ 1 − δ.
- *Cost.* The SFT search costs poly(log p, 1/τ, log(1/δ))·T. The phase measurements cost
  O(τ^{−2}·log p·log(log p/δ)) evaluations of f per list entry.

**Proof.**
1. *Frequency list.* F(t) = f(g^t) on Z/n has F̂(j) = ⟨f, χ⟩. Since g is a generator, restriction classes are
   singletons (Theorem T with i = 1). The sparse Fourier transform over Z_n with known n returns a list of O(τ^{−2})
   candidates containing every τ-heavy frequency (Akavia–Goldwasser–Safra, FOCS 2003, Thm 5; corrected citation,
   codex). Mansour's algorithm assumes a power-of-two modulus, and Akavia 2009 Thm 4 depends on the Fourier ℓ¹
   norm.
2. *One phase.* For any z, E_s f(zg^s)e(−js/n) = F̂(j)·e(j·log_g z/n). Estimate this and F̂(j) to additive accuracy
   cτ. Their quotient determines θ(z) = j·log_g z/n mod 1 to within 1/16 turn, using O(τ^{−2} log(1/δ′)) samples.
3. *Doubling.* Measure θ_k ≈ frac(2^k θ) at z = y^{2^k}, for k = 0..K with K = ⌈log₂ n⌉ + 3. Go backwards: of the
   two half-angle lifts of θ_{k+1}, one has half its error and the other is half a turn away. So 1/16-accurate
   measurements select the right one (it suffices that 2ε + e < 1/2; proof critic). After K steps, round to get
   j·log_g y mod n. With a = gcd(j, n) this gives log_g y ≡ (jx/a)(j/a)^{−1} mod n/a.
4. *Confidence.* A union bound over the SFT and the K + 2 estimates, with δ′ = δ/(K + 3). ∎

**Notes.**
- *Full log.* Cofactor completion solves yg^{−x₀} = (g^{n/a})^t, 0 ≤ t < a, by BSGS in O(√a) operations.
  Self-testing on y = g^r with known r identifies the list entry that works.
- *No generator given.* Retry random bases, complete full logs and verify. A partial log alone does not certify
  success. For a prescribed base b, recover verified logs of b and y to a common generator and divide.
- *Quotient versus subgroup (codex).* The algorithm gives DL in the quotient F_p^*/ker χ. That is DL in an order-ℓ
  subgroup only when ℓ ∥ ord χ. If ℓ² | p − 1, an order-ℓ character can be trivial on the order-ℓ subgroup. The
  near-full-order corollary below avoids this by bounded cofactor completion.
- Code: `src/prepread/dlred.py` (the operational part, given the frequency). Tests:
  `test_dlog_reduction_planted`, `test_dlog_reduction_noncoprime_frequency`.
- Codex proved the same per-prime statement in its blind round and labelled the resulting barrier CONJECTURED.
  Ours adds the N-level corollary.

**Corollary O1-N (an O1 factoring observable is a prime-field DL algorithm). PROVED.**

*Hypotheses* (made explicit after the adversarial round):
- (i) The programs Φ_X of size s are N-independent and constructible in poly(log X) time, or given as advice.
  With advice the conclusion is nonuniform.
- (ii) The multiplier N ↦ h(N; ρ) is computable in polynomial time with coins ρ.
- (iii) Uniform bounds i ≤ I and κ ≤ K hold.
- (iv) With probability ≥ μ over a uniform balanced pair (p, q) ∈ (X, 2X]² and the coins ρ, the observable
  Ψ(y) = e(hΦ_X(y)/N) has, for some a ∈ (Z/N)^* of index ≤ I, a τ-heavy peak of order ≥ r/K.

*Conclusion.* For at least a μ/2 fraction of primes p ∈ (X, 2X], full discrete logarithms in F_p^* are computable
in expected time poly(s, log X, 1/τ, I, K, 1/μ) plus the cost of h.

*Proof.*
1. *Correlation.* Theorem T (T2) gives χ = χ_pχ_q with |Ψ̂(χ)| ≥ √τ/I and ord χ ≥ r/K. By the CRT factorisation,
   |A_p(χ̄_p, c_p)| ≥ √τ/I with c_p = hq̄ mod p, since the q-factor is ≤ 1.
2. *Order.* (p − 1)(q − 1)/(KI) ≤ r/K ≤ ord χ ≤ ord χ_p·(q − 1), so ord χ_p ≥ (p − 1)/(KI).
3. *Averaging.* At least μ/2 of the primes p have at least a μ/2 fraction of good (partner, coins) choices.
   Reverse averaging even gives μ/(2 − μ) (codex).
4. *Simulation at p.* Sample q′ ≠ p uniformly among the primes of (X, 2X] by rejection (this uses a primality
   test). Fix coins ρ and set c = h(pq′; ρ)·q′^{−1} mod p. This needs neither the factor nor the witness a. Run
   O1-DL on f = e(cΦ_X/p).
5. *Completion.* Keep the list entries with gcd(j′, p − 1) ≤ KI, complete the cofactor by BSGS, and verify
   g^x = y. Each trial succeeds with probability ≥ μ/2. ∎

### 2.5 The barrier, stated three ways

1. **CONDITIONAL (from O1-DL).**
   - Suppose a family of size-s programs has τ-correlations with near-full-order characters (order ≥ (p − 1)/K)
     at a positive fraction of primes.
   - Then it yields prime-field DL there in time poly(s, log p, 1/τ, K).
   - So any hardness of prime-field DL beyond poly(s, 1/τ) excludes it. The best known algorithms take
     L_p[1/3, 1.92].
   - The first draft displayed an explicit square-root trade-off τ ≤ Õ((s log² p/T_DL)^{1/2}). It ignored the
     poly(1/τ) cost of the SFT search and was deleted (codex; proof critic). It holds only when the frequency is
     supplied as advice.
2. **Not an order-only statement** (Prop PH). Smooth or annihilated orders must be excluded. The hardness enters
   through the large prime factors of the order.
3. **Unconditional** (Weil, §2.1; BGK for monomials, §2.2).
   - τ ≤ w_p√p/(p − 1). This is vacuous for s ≥ ½·log₂ p.
   - An unconditional bound τ ≤ poly(s)·p^{−c} for all size-s SLPs and hard-order characters would be a
     superpolylogarithmic algebraic-circuit lower bound for weak DL predictors.
   - Coppersmith–Shparlinski (J. Cryptology 13 (2000) 339–360) leave the SLP complexity of DL open. The known
     lower bounds concern degree, sparsity or linear complexity: Coppersmith–Shparlinski Thms 1 and 3;
     Lange–Winterhof, Acta Arith. 101 (2002); Meidl–Winterhof, IEEE-IT 47 (2001).

**What remains open in O1.** O1-DL shows that a short SLP with large hard-order correlation is a polylog
algebraic-circuit weak DL predictor for F_p^*. We know of no reason to expect one, and no proof that none exists.

---

## 3. The Galois framing of the dense-torus question

### 3.1 Dictionary, and the right statistic. KNOWN / PROVED (routine)

- **Fixed points.** For F ∈ Z[x] and p ∤ lc(F)·disc(rad F),
  Z_p := #{x ∈ F_p : F(x) = 0} = fix(Frob_p) on the distinct roots. This holds also when deg F ≥ p.
- **Burnside and Markov.** For G acting with t orbits, E_G fix = t, and Pr_G[fix ≥ M] ≤ t/M.
- **What model U asks.** A specific element, Frob_p, must fix ≥ μp roots for every p ∈ (X, 2X], with
  deg F ≥ μX. That is a statement about the Frobenius sequence at small primes, not about the group.
- **FPP is the wrong statistic.** The fixed-point proportion FPP(G) = Pr[fix > 0] only says whether a root exists.
  - For transitive actions on D points the mean fixed-point count is 1, and FPP can be anywhere in [1/D, 1 − 1/D].
    AGL₁(F_D) attains 1 − 1/D (proof critic and codex; the draft's "1 − e^{−1}" was only the S_D value).
  - So "FPP ∈ [1/polylog, 1 − 1/polylog]" is plentiful and says nothing about density.
  - The relevant set {σ : fix(σ) ≥ μp} has Haar measure ≤ t/(μp).

### 3.2 Prop CS (what effective Chebotarev can and cannot do). CONDITIONAL on GRH. Corrected after review

The draft claimed that Chebotarev is "silent by a factor √X" whenever deg F ≥ μX. Codex showed this is false
twice over:
- bad primes were not accounted for (F = x(x − M) with M = ∏_{X<p≤2X} p);
- the target was wrong. Dense roots need Σ_p Z_p ≳ μX·π*(X) ≈ μX²/log X, not the main term t·π*(X).

**Corrected statement.**
- Let F have irreducible factors F_i with fields K_i of degree n_i, Σ n_i = deg rad F =: D.
- Let B_X be the set of primes in (X, 2X] dividing lc(F)·disc(rad F).
- Under GRH (Lagarias–Odlyzko 1977, Thm 1.1, the prime ideal theorem form):

    Σ_{p∈(X,2X]∖B_X} Z_p(F) ≤ t·π*(X) + O( √X·Σ_i (log d_{K_i} + n_i log X) ).

- **Consequence.** If Z_p ≥ μp at all but one prime and |B_X| ≤ π*(X)/2, then one of the following holds:
  - t ≳ μX (many irreducible factors: the Pochhammer route);
  - D ≳ μX^{3/2}/log² X;
  - Σ log d_{K_i} ≳ μX^{3/2}/log X.
- **Escapes.** W_U is possible only with many factors, degree ≳ X^{3/2}, huge discriminants, or bad primes at
  half the interval.
- **Reading.** Short SLPs reach degree X^{3/2} at size ≈ 1.5·log₂ X. Compositions of depth ≈ log₂ X have
  discriminants with log ≈ D log D and possibly many bad primes in the interval. So CS is a genuine but weak
  conditional filter. MIN (§4.4) is the unconditional statement for the composition class.
- The unconditional range in Lagarias–Odlyzko Thm 1.3 (x ≥ exp(10 n_K (log d_K)²)) is a sufficient range for that
  theorem, not a limit on all methods (codex).

### 3.3 Prop GAL-0 (no Galois-only barrier). PROVED

Codex found the S_D case independently. Fix X and D = ⌊X/3⌋. Then there are monic F ∈ Z[x] of degree D with
Z_p(F) = D, i.e. δ_p = D/p ∈ (1/6, 1/3], at every prime p ∈ (X, 2X]. Their Galois groups can be:

- **(a) Trivial.** F = ∏_{j≤D}(x − j). The roots stay distinct mod p because D < p.
- **(b) S_D.** F = ∏_{j≤D}(x − j) + M·h(x), with M = ∏_{X<p≤2X} p and deg h < D.
  - Choose h by CRT with three auxiliary primes ℓ₁, ℓ₂, ℓ₃ ∉ (X, 2X], ℓ_i ∤ M, and ℓ₂, ℓ₃ > D:
    - F is Eisenstein at ℓ₁ (congruences mod ℓ₁²);
    - F is squarefree of type (D − 1, 1) mod ℓ₂;
    - F is of type (2, 1, …, 1) with distinct factors mod ℓ₃.
  - Eisenstein gives transitivity, the (D − 1)-cycle gives double transitivity, and conjugating the transposition
    gives all transpositions. So the group is S_D.
- **(c) Cyclic of order D.**
  - Take a prime ℓ > 2X that splits completely in Q(ζ_D, p^{1/D} : p ∈ (X, 2X]). Such ℓ exist by Chebotarev.
  - Then ℓ ≡ 1 (mod D) and every p ∈ (X, 2X] is a D-th power mod ℓ. So every such p splits completely in the
    degree-D subfield K ⊂ Q(ζ_ℓ).
  - Choose θ ∈ O_K by CRT with D distinct residues at the D primes above each p (possible since D < p). Then θ is
    primitive, and its minimal polynomial has D distinct roots mod each p.
  - Codex and the proof critic both checked this case.

Codex's variant prescribes ⌊p/2⌋ linear factors at each p, giving δ_p ≈ 1/2 with group S_D.

**Consequence.**
- A barrier must use the SLP structure.
- Any barrier covering the trivial-group case must cover the unit multiples 2^{e(n)}∏_{a≤n}(x − a). By Bürgisser,
  Comput. Complexity 18 (2009), preprint Thm 4.1(2), that would make τ(Per) superpolynomial (slpcirc Prop BB(b)).
- So a barrier "short SLP + Galois group in 𝒢 ⇒ no W_U" is permanent-hard whenever 𝒢 contains the trivial
  group. Hardness of the unscaled monic Pochhammer program alone does not follow from the quoted construction
  (codex).

### 3.4 Abelian, non-abelian, and "pinned Frobenius"

**Prop AB (abelian ⇒ cyclotomic containment). PROVED (Kronecker–Weber).**
- If F has abelian splitting field of conductor m, its roots lie in Q(ζ_m) = Q(G_m[m]).
- For p ∤ m·disc(rad F), Z_p(F) = Σ_χ n_χ χ(p). Here χ runs over Dirichlet characters mod m, the n_χ ≥ 0 are
  integers, Σ n_χ = deg rad F, and n_{χ₀} = t.
- Each irreducible factor has a regular quotient action, so it contributes its full degree or nothing.
- p − 1, Williams p + 1 and Bach–Shallit Φ_k (Gauss periods) are of this form.
- Prop AB establishes cyclotomic containment and dependence on p mod m. It does not turn every such polynomial into
  an efficient torus algorithm, and it gives no size barrier: GAL-0(c) is abelian.

**Non-abelian does not mean "not group order".**
- Division polynomials of a non-CM curve have Galois image ⊆ GL₂(Z/M)/±1.
- For odd M, away from bad primes, Z_p(ψ_M) = (|E(F_p)[M]| + |E^{tw}(F_p)[M]| − 2)/2 (codex; corrected).
- The fixed-space sizes depend on the full Frobenius action, not only on its characteristic polynomial: identity
  and a nontrivial unipotent have the same characteristic polynomial.
- This is ECM (nonsmooth Prop OA). So abelian versus non-abelian is the wrong dividing line.

**CONJECTURED (Frobenius pinning).** Every polylog N-independent SLP with W_U falls into one of three classes:
- Frob_p acting on torsion of a commutative group scheme (G_m-tori, abelian varieties). Multiples [M] have SLP
  size O(log M) via the group law;
- dynamics with a pinned critical structure;
- the explicit-root route (trivial group).

Kronecker–Weber and Ritt–Eremenko motivate this. They do not prove it. Prop BB shows that the third class cannot
be excluded without a permanent lower bound.

---

## 4. Quadratic compositions: toward the inverse theorem

**Setting.**
- g_i(x) = x² + c_i with c_i in a fixed finite set C ⊂ Z (C = {−2..2} in the experiments), and
  F = g_k ∘ ⋯ ∘ g_1.
- For odd p, set T_k = {0} and T_{j−1} = g_j^{−1}(T_j). Then Z_p(F) = |T_0|, and for nonempty T_j

    |T_{j−1}| = Σ_{u∈T_j} (1 + η(u − c_j)) = |T_j|(1 + β_j),

  with η the Legendre symbol (η(0) = 0).
- In Galois terms, β_j is the "non-swap" frequency of Frob_p at the fixed vertices of level j of the preimage
  tree.

### 4.1 Lemma DEF (deficit). PROVED, elementary

**Statement.** If Z_p(F) ≥ μp > 0 with k = log₂(μp) + Δ, then Σ_j (1 − β_j) ≤ 2 ln 2·Δ ≈ 1.39Δ.

*Proof.* Σ_j (1 − log₂(1 + β_j)) = k − log₂ Z_p ≤ Δ. The function −log₂(1 + β) is convex, and its tangent line
at β = 1 gives 1 − log₂(1 + β) ≥ (1 − β)/(2 ln 2). ∎

**Reading.** Near-minimal-depth dense compositions need almost every level to be almost entirely
"square-shifted". At a single prime, nothing can close this: the T_j are roots of univariate polynomials, and
character sums over such 0-dimensional sets have no square-root-cancellation theory.

### 4.2 Prop OUT (outer constant). PROVED (codex), trivial

A root of F mod p gives a root of g_k. If −c_k is not a square in Z, there are no roots at the primes with
(−c_k/p) = −1. That is half the primes asymptotically, uniformly over the finite alphabet. So W_U forces the
outermost map to be x² − s². OUT excludes about half the primes; it does not prove QC's o(1) conclusion.

### 4.3 Prop PFX (fixed-prefix towers). PROVED (codex)

**Statement.**
- Let P_1, P_2, … be a nested tower: P_{m+1} = P_m ∘ g for maps g from the alphabet.
- Let F_X be compositions whose depth-m outer prefix equals P_m for every m ≤ k(X), with k(X) → ∞.
- If FPP(Gal(P_m)) → 0, then the proportion of primes of (X, 2X] where F_X has any root tends to 0.

Iterates of x² + c with full wreath image qualify: Stoll 1992, as stated in Jones 2013 (Publ. Math. Besançon),
Thm 2.3, covers x² + 1.

*Proof.* A root of F_X maps to a root of P_m. Fixed-depth Chebotarev bounds the proportion by FPP(G_m) + o(1).
Let m → ∞. ∎

For the full binary wreath product, q_{n+1} = q_n − q_n²/2 and q_n ~ 2/n (Odoni 1985). For bounded degree the FPP
is O(1/n) (codex). E3 tests this Haar model in the silent regime.

### 4.4 Theorem MIN (no dense compositions of near-minimal degree). PROVED; no prior formulation found by us

**Statement.**
- Let 𝒢 be a finite set of monic polynomials in Z[x] of degree ≥ 2, and fix μ ∈ (0, 1] and Δ ≥ 0.
- Then

    sup_{F} #{p ∈ (X, 2X] : Z_p(F) ≥ μp} / π*(X)  →  0  as X → ∞,

  where F ranges over all compositions g_k ∘ ⋯ ∘ g_1 with g_i ∈ 𝒢 and deg F ≤ 2^Δ·X.
- In particular W_U(μ) fails for all such F once X is large. Since W_U needs deg F ≥ μX, **a dense composition
  needs degree ≥ X·2^{Δ(X)} with Δ(X) → ∞.**

The first draft concluded only "fails at a positive proportion κ ≈ 1/|G_P|". Both the proof critic and codex
pointed out that Markov in step 5 gives the "almost all primes" form.

**Proof.**
1. *Invariant disc.* Let degrees be ≤ D₀ and lower coefficients be ≤ H in absolute value, and set R = D₀H + 2
   (codex's simplification). For |z| = r > R, |g(z)| ≥ r^{d−1}(r − dH) > r. So every preimage of the closed
   R-disc lies inside it, and by induction from 0 every complex root of every composition has all conjugates in
   |z| ≤ R. For x² + c with |c| ≤ H one may take R = (1 + √(1 + 4H))/2.
2. *Northcott count.* Let N(d) be the number of monic irreducible integer polynomials of degree ≤ d with all roots
   in |z| ≤ R. N(d) is finite, since the coefficients are bounded by binom(d, i)R^i. For R = 2 and d = 1..4,
   N = 5, 31, 263 and ≈ 6113 (proof critic).
3. *Prefix transfer.* Write F = P ∘ Q with P the outer m layers. Q is monic of degree deg F/deg P, so each root
   of P mod p has ≤ deg Q preimages. Hence Z_p(F) ≥ μp > μX implies
   #{y ∈ F_p : P(y) = 0} ≥ μX·deg P/deg F ≥ a·deg P, with a := μ2^{−Δ}.
4. *Few factors.* Fix ε > 0. Pick d with d + 1 > 2/(aε), then m with 2^m > 2N(d)/(aε). Every m-layer prefix P has
   deg P ≥ 2^m, and its number t_P of distinct irreducible factors satisfies
   t_P ≤ N(d) + deg(rad P)/(d + 1) < εa·deg P. Multiplicities only enlarge deg P (codex).
5. *Markov + Chebotarev on finitely many prefixes.* For each of the finitely many m-layer words P over 𝒢, G_P acts
   on the distinct roots of rad P.
   - By Burnside and Markov, the set S_P = {σ : fix(σ) ≥ a·deg P} has Haar measure ≤ t_P/(a deg P) < ε.
   - By Chebotarev for the fixed field (dyadic form, from π_C(x) ~ (|C|/|G|) li(x)), the primes of (X, 2X] with
     Frob_p ∈ S_P have proportion < ε + o(1).
   - At every other prime with p ∤ lc·disc(rad P), P has < a·deg P roots, so Z_p(F) < μp by step 3.
   - The bad primes are finitely many.
6. *Uniformity.* There are finitely many prefixes, and compositions with fewer than m layers have
   Z_p ≤ deg F ≤ D₀^{m−1} < μX for large X. So for X ≥ X₀(ε) every admissible F has density ≥ μ at < 2ε of the
   primes. ∎

**Remarks.**
- *Where it stops being group theory.* The non-group-theoretic step is Northcott: roots of compositions over 𝒢
  have bounded house, so a prefix has only boundedly many low-degree irreducible factors. This is exactly what
  GAL-0 shows Galois theory alone cannot supply. The Pochhammer product, whose roots 1..D have unbounded height, is
  the excluded case.
- *Counterexample attempts* (codex; proof critic). Powers of x², fully split prefixes (finitely many integer roots
  in [−R, R]), x² − 1, PCF behaviour, repeated roots and Chebyshev/power families all satisfy the low-degree-factor
  bound.
  - The proof critic enumerated all 5^m prefixes for m ≤ 7. The maximum t_P is m + 1, attained only by
    (x² − 1)∘(x² − 2)^{∘(m−1)}.
  - Its empirical root frequencies at m = 5, 6 never exceeded the Markov bound beyond noise.
- *Scope.* Fixed μ and Δ are essential.
  - A diagonal argument allows some arbitrarily slowly growing Δ(X).
  - Under GRH, log N(d) = O(d²) suggests Δ ≈ c·log log log X. That is a sketch: the splitting-field estimates are
    not supplied (UNVERIFIED).
  - The method does not reach polylogarithmic depth, or μ = (log X)^{−B}.
- *General SLPs.* MIN is for compositions. A general SLP has no outer prefix of bounded degree.
- *Literature* (codex). Related results have different hypotheses:
  - Doyle–Healey–Hindes–Jones (arXiv 2108.11233), Thm 1.5: assumes large arboreal image;
  - Ferraguti (Proc. AMS 2018), Thm 2.3: stable primes;
  - Benedetto–Juul and Looper: surjective arboreal representations.

  None states MIN.

### 4.5 Random coefficients. PROVED, trivial

- If c_1..c_k are independent and uniform mod p, then E_c Z_p(F) = 1 exactly, because E_c η(u − c) = 0 at each
  layer. So Pr_c[Z_p ≥ μp] ≤ 1/(μp).
- For integer c_i uniform in [0, H) this becomes an upper bound E Z_p ≤ (1 + p/H)^k.
- So generic large coefficients are excluded at every single prime. The interesting case is small fixed constants.

### 4.6 What remains: Conjecture QC (codex's restricted target, endorsed)

**Conjecture QC.** For fixed finite C ⊂ Z and A, B > 0, uniformly over words of length ≤ (log X)^A,

  #{p ∈ (X, 2X] : δ_p(F) ≥ (log X)^{−B}} = o(π*(X)).

**What is proved**, for fixed μ, which does not settle QC's shrinking density:
- OUT: the outer constant is −s²;
- PFX: fixed prefix towers with vanishing FPP fail;
- MIN: the excess depth must → ∞;
- DEF: near-minimal depth needs almost-total squareness at each prime.

**Not proved.** That the remaining words are only "intersective" ones: words whose prefixes have roots at almost
every prime, e.g. those carrying an integer backbone 0 ← y_1 ← y_2 ← … with y_j² + c = y_{j−1} and |y_j| ≤ 2.
That is a description of what E2's maximisers look like, not a theorem.

---

## 5. Where the limit bites (summary)

| Question | Per-prime status | Uniform status | What is proved here |
|---|---|---|---|
| O1, χ of near-full order (hard prime factor ℓ) | correlation τ ⇒ DL mod ord χ in poly(T, log p, 1/τ) (O1-DL) | efficient O1 factoring ⇒ full prime-field DL at ≥ μ/2 of primes (O1-N, hypotheses (i)–(iv)) | CONDITIONAL barrier on prime-field DL; Weil only for width < √p |
| O1, χ of smooth/annihilated order | easy (PH) | arithmetic N-independent when p − 1 = 2^e·u with u \| m; multiplier p-dependent | examples only; this is the OA regime |
| O1, monomial Φ | power saving for ord χ ≥ p^ε (BGK) | — | §2.2 |
| dense-torus, unrestricted | trivial for prime-specific polynomials | needs uniformity | no Galois-only barrier (GAL-0); a general barrier is permanent-hard (BB) |
| dense-torus, GRH filter | — | W_U ⇒ many factors, deg ≳ X^{3/2}, huge discriminants, or bad primes at half the interval | CS (CONDITIONAL) |
| compositions, outer −c non-square | — | no roots at ≈ ½ of primes | OUT |
| compositions, nested prefix tower with FPP → 0 | — | roots at o(1) of primes | PFX |
| compositions, degree ≤ 2^Δ X | — | density ≥ μ at o(1) of primes, uniformly | **MIN** |
| compositions, polylog depth | — | open | Conjecture QC; E2 and E3 consistent with it |

---

## 6. Experiments (pre-registered in `notes/prepread_prereg.md`; post-hoc runs labelled)

Floating point appears only in E1/E1b spectra (float64 FFT; relative error ~1e-12, orders of magnitude below every
comparison). All root counts are exact. Peak RSS ≤ 1.2 GB (E1, at p ≈ 2^{15.8}). One heavy job ran at a time.

### 6.1 E1: large-prime-order correlations (6 primes: 4 in (2^13, 2^14) with all c; 51157 and 56779 with 256 c)

- **Baseline.** B_p = 0.040, 0.039, 0.038, 0.036 (small primes) and 0.0179, 0.0173 (large primes). These are
  medians of the maximum over c and j ∈ J_p for random maps.
- **H1** (nondegenerate F1–F3 members, image ≥ p/4): 300 member-prime pairs, **0 with ratio > 2** (maximum
  1.43). **Holds.**
- **H2** (Weil sanity, width < √p): 270 checks, **0 violations**. The maximum M/bound is exactly 1.0, attained by
  linear Φ (Gauss sums). **Holds.**
- **H3** (planted DL control): ratio 25–58 at every prime (M = 1.00). **Holds** (≥ 20).
- **Kill criterion** (non-control ratio ≥ 10 at ≥ half the primes): **not triggered.** The largest non-control
  ratio is 6.6, for x^{lcm(1..30)} at p = 13841 (image 173, p − 1 = 2⁴·5·173). That is order annihilation,
  explained in §2.2. The member is "degenerate" by the pre-registered image criterion. Its ratio is 1.76 at
  11657 and below 1.5 at the other primes.
- **Degenerate members.** All F2 compositions have small images (image of f^{∘k} ≈ 2p/k), so they were outside
  H1 as pre-registered. Descriptively, they have ratios 1.3–1.57 at every prime, all with width ≥ √p.
- **Disclosed weakness of the pre-registered F1 generator.** Uniformly random operands produced
  low-degree programs: all 360 F1 member-prime pairs have width < √p, and 90 are constant. So F1 tests the Weil
  regime, not O1. The post-hoc E1b fixes this.

### 6.2 E1b (POST-HOC): 30 high-degree chain SLPs (sizes 24, 32, 40; degrees 2^9–2^24), same primes

This run was added after the F1 generator weakness was seen. Its generator, `slp.chain_slp`, always operates on the
previous node and squares often. The statistic and baseline are those of E1.

| | value |
|---|---|
| member-prime pairs | 180 (12 constant at some prime) |
| width computed | 153 (formal degree < p − 1, or small prime) |
| of which w_p ≥ √p | 131 |
| median ratio M_p/B_p, w_p ≥ √p | 1.01 |
| max ratio, w_p ≥ √p | 1.44 |
| max ratio, all pairs | 1.97 (E1b_s24_2, w_p = 40 at every prime: low effective degree, Weil regime) |
| ratio ≥ 10 (kill criterion analogue) | 0 |
| Weil checks (w_p < √p) | 14, 0 violations |

**Reading.** Dense SLPs of degree up to 2^24, far beyond p, behave like random maps against large-prime-order
characters. Every observed ratio lies within 2× of the random baseline. This is experimental support for the
O1-DL picture at toy sizes. It is not evidence about cleverly designed SLPs.

### 6.3 E2: exhaustive search over all 5^k quadratic compositions (C = {−2..2}; C helper cross-checked against brute force in tests)

**Per prime** (k = b + 2): M* = max over all 5^k words of the fibre.

| b | primes | M* | k |
|---|---|---|---|
| 7 | 139, 179, 211, 233 | 43, 34, 36, 54 | 9 |
| 8 | 337, 353, 383, 479 | 97, 102, 164, 86 | 10 |
| 9 | 787, 911, 971, 991 | 70, 92, 68, 108 | 11 |
| 10 | 1031, 1093, 1511, 1697 | 102, 92, 82, 168 | 12 |

- **H4a** (M* ≤ max(2k², 4·2^{max v₂(p ± 1)})): **holds at all 16 primes.** The largest ratio is 164/512 at
  p = 383, where p + 1 = 2^7·3.
- **H4b** (median at b = 10 is ≤ 3× the median at b = 7): 97/39.5 = 2.46. **Holds.** The medians are not
  monotone (99.5 at b = 8, 81 at b = 9). A density-bounded mechanism would give ≈ 8×. Densities M*/p fall from
  0.17–0.31 at b = 7 to 0.05–0.10 at b = 10.

**Uniform** (8 primes jointly): U* = max over words of the minimum fibre.

| b | k | U* | nodes |
|---|---|---|---|
| 9 | 11 | 8 | 1.2·10^7 |
| 9 | 13 | 10 | 3.1·10^8 |
| 10 | 12 | 12 | 6.1·10^7 |
| 10 | 14 | 14 | 1.5·10^9 |

- **H5** (U* ≤ k): **holds**, with equality at b = 10. The maximum minimum density is ≈ 14/1500 ≈ 0.009.
- **All argmax words have outermost constant −1** (c_k = −1), as OUT requires. They contain long runs of 0 and −2
  layers, i.e. power/Chebyshev pieces.
- *Disclosed:* one timing run before the pre-registered runs (k = 11, p = 1031, M* = 88) was a development check
  outside the pre-registered grid.

### 6.4 E3: x² + 1 iterates at all 4642 primes in (2^14, 2^16), k = 1..24

| k | 1 | 2 | 4 | 8 | 12 | 16 | 20 | 24 |
|---|---|---|---|---|---|---|---|---|
| fraction with a root | 0.499 | 0.375 | 0.255 | 0.163 | 0.127 | 0.104 | 0.085 | 0.076 |
| Haar wreath q_k | 0.500 | 0.375 | 0.258 | 0.164 | 0.121 | 0.096 | 0.080 | 0.069 |
| mean #roots | 1.00 | 1.00 | 1.00 | 1.05 | 1.09 | 1.04 | 1.05 | 1.05 |

- **H6 holds for all 24 depths** (tolerance 3σ + 0.01). The maximum fibre over all primes and k ≤ 24 is 84.
- *Disclosed nuance.* For k ≥ 11 the observed fraction exceeds q_k by 0.005–0.007, about 1.5–2σ. This is
  systematic, and the pre-registered 0.01 slack absorbs it. A plausible source is finite-p effects: near short
  cycles the preimage tree is not a free branching process. We did not investigate further.
- **Reading.** For k ≥ 15 the degree 2^k exceeds p, the regime where Chebotarev is silent. The Frobenius at these
  small primes still looks Haar-distributed on the tree: critical branching, mean 1, survival ~2/k. This is
  evidence for QC's heuristic, not a proof.

### 6.5 E4: the DL reduction in operation

The reduction code sees only p, g, a target y and oracle access to f.

| prime | (a) planted: recovered | queries / recovery | (b) best E1 member (ratio) | (b) recovered | queries / recovery |
|---|---|---|---|---|---|
| 11657 | 20/20 | 3 618 | Dickson 2^16 (1.38) | 20/20 | 1.17·10^6 |
| 12161 | 20/20 | 3 618 | Dickson 2^10 (1.43) | 20/20 | 1.17·10^6 |
| 13249 | 20/20 | 3 618 | Dickson 2^16 (1.42) | 20/20 | 1.25·10^6 |
| 13841 | 20/20 | 3 618 | Dickson 2^16 (1.42) | 20/20 | 1.36·10^6 |
| 51157 | 20/20 | 4 020 | Dickson 2^10 (1.40) | 20/20 | 6.42·10^6 |
| 56779 | 20/20 | 4 020 | Dickson 2^16 (1.43) | 20/20 | 6.50·10^6 |

- **H7 (amended A1) holds.**
  - (a) 120/120 recoveries with ≤ 4 020 ≤ 10^4·⌈log₂ p⌉ queries.
  - (b) 120/120 recoveries, but with ≈ 100·p queries per recovery, worse than exhaustive search.
- The original H7(b) prediction (0/20) was mis-specified and was amended before running. The estimator is
  unbiased for any f with F̂(j) ≠ 0. A flat observable makes the reduction cost Θ(1/|F̂|²) ≈ Θ(p), exactly as
  O1-DL predicts.
- In (b) the (c, j) pair came from E1's exhaustive search over c and the full FFT (HARNESS; cost p²). A real
  attack would pay that too.

### 6.6 E5: random-SLP root-density census; E5b (POST-HOC) with chain SLPs

- **E5** (2000 pre-registered random SLPs of size 20, 12 primes in (2^13, 2^14)).
  - 171 programs have min density ≥ 1/64. **All 171 are the zero polynomial** (identically 0 at all 12 primes),
    which H8 excludes. **H8 holds.**
  - Among the 1829 nonzero programs the largest minimum density is 3.3·10^{−4}, i.e. ≈ 3 roots per prime.
  - *Disclosed.* The same generator weakness as F1 applies: no nonzero program has degree ≥ 2^13.
- **E5b (post-hoc)** (2000 chain SLPs of size 40, 1949 of them of degree ≥ 2^13, median degree 2^{19.1}).
  - 33 programs are identically zero.
  - Among the high-degree nonzero programs, the largest minimum density over the 12 primes is 5.9·10^{−4} (≈ 5–8
    roots per prime). The median density is 2.1·10^{−4} (≈ 2 roots).
  - This is the Chebotarev-like "about one root per irreducible factor" picture of slpcirc E2 (H2), now at degrees
    2^13–2^24 far beyond p.

---

## 7. Candidate cards (BRIEF §4 template, condensed)

- **C1. O1 with dense high-degree SLPs** (compositions, chain SLPs, sparse binomials).
  - Mechanism: algebraic preprocessing creates character-like structure.
  - Results: E1/E1b ratios ≤ 2 except OA. O1-DL: any success is a DL algorithm.
  - **REFUTED for this model experimentally** at the tested sizes; **CONDITIONALLY barred** in general.
- **C2. O1 with PH-type programs.**
  - Correlation ≈ 1 at primes with p − 1 = 2^e·u, u | m.
  - **Reduced to OA (p − 1 method).**
- **C3. O1 with monomials.** BGK power saving. **REFUTED (proved)** for ord χ ≥ p^ε.
- **C4. Dense-torus via quadratic compositions.**
  - OUT, PFX, MIN (proved).
  - E2: the uniform maximum fibre is ≤ k. E3: Haar-like Frobenius.
  - **REFUTED for near-minimal degree (proved).** Polylog depth: **open**, no candidate.
- **C5. Dense-torus via Galois-designed polynomials** (CRT; trivial, cyclic or S_D group).
  - These exist (GAL-0), but no short SLP is known.
  - **Reduced** to the Pochhammer/τ question (Prop BB).
- **C6. Gleason-type parameter polynomials H_{j+1} = xH_j² + 1** (codex blind round).
  - About one root per prime in exact counts up to X = 2048.
  - **Killed (heuristic + data).**

---

## 8. Reproduction

```
.venv/bin/python -m pytest -q tests/test_prepread.py
.venv/bin/python -m src.prepread.run_experiments e1    # ~28 min, ≤ 1.2 GB
.venv/bin/python -m src.prepread.run_experiments e1b   # post-hoc
.venv/bin/python -m src.prepread.run_experiments e2    # ~2 min (C helper compiled on the fly with cc -O2)
.venv/bin/python -m src.prepread.run_experiments e3    # < 1 s
.venv/bin/python -m src.prepread.run_experiments e4    # needs e1 output
.venv/bin/python -m src.prepread.run_experiments e5
.venv/bin/python -m src.prepread.run_experiments e5b   # post-hoc
.venv/bin/python -m src.prepread.analyze e1 e1:e1b
```

## 9. Review record

REVIEW_PLACEHOLDER
