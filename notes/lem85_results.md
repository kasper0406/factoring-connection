# lem85: a self-contained proof of Lemma 8.5′ (restricted Bombieri–Vinogradov for moduli ℓ·(smooth))

Agent `lem85`, phase 11, 2026-10-20.
- Code (exact checks of the elementary identities only): `src/lem85/identities.py`; tests `tests/test_lem85.py`.
- GPT-6.1-sol adversarial referee: `codex/lem85_1.{md,out}`. Opus critics (proof; citations): §10.
- No experiments, so no pre-registration.

## Verdict

**Lemma 8.5′ is now proved (status: proved), in the generality Theorem A needs. Part 1, the
referees' weakest link, is closed in the sense the brief asked for.** No step uses the inside of
another paper's proof. Harper's preprint arXiv 1208.5992 is no longer used as a direct input. It appears only as the
provenance of one published lemma, below. Lee–Venkatesan's
Lemma 8.5 is no longer used, neither its statement nor its proof.

In plain terms:
- **Part 1 needed no character majorant.** The old proof imported LV's majorant W from inside LV's
  proof. That import was unnecessary. For t = ℓs with ℓ a prime > y, the "χ_ℓ principal" part of the
  error mod t is *exactly* the error mod s divided by ℓ − 1 (lift averaging, §3 E2). What is then
  needed is a restricted Bombieri–Vinogradov theorem over smooth moduli s, with a *relative* saving
  (Lemma SM, §5). We prove Lemma SM ourselves.
- **Large conductors** (conductor ≥ y^K, K fixed) are handled by our Lemma LS (§4). It uses only the
  classical multiplicative large sieve (Iwaniec–Kowalski Thm 7.13), a truncated Perron kernel
  bound (proved in §2), and the smooth-number factorisation n = mk. At y = L[1/2] every smooth-number density is
  y^{−O(1)} = x^{−o(1)}, so crude losses of y^{O(1)} are affordable. This replaces Harper's
  Proposition 2 (LV Fact 8.11).
- **Small conductors** (2 ≤ r ≤ y^K) use one published, numbered statement: **Drappeau, Compositio
  Math. 151 (2015), Lemme 5**, restated as Pascadi, Compositio Math. 161 (2025), Lemma 3.2. It is an
  averaged bound for smooth-number character sums of conductor ≤ x^η. We apply it only to
  *primitive* characters, which makes it independent of the moduli. Möbius inversion plus de la
  Bretèche–Tenenbaum's local bound Ψ(x/d, y) ≪ Ψ(x, y)d^{−α} transfer it to smooth moduli with
  a relative saving.
  - The transfer costs only R₁^{1−α} = (log x)^{O(1)}. That is because log R₁ ≍ log y and
    (1 − α) log y ≍ log(u log y). The saving H(u)^{−δ} absorbs this loss because u/log²u ≫ log log x.
  - Both facts are specific to y = L[1/2, a]. This regime is exactly where the argument works.
- **Exceptional (Siegel/Page) zeros.** Our own steps use no information about zeros. Drappeau's
  Lemme 5 is used with A = 1/2, and its implied constant is stated to be effective for A < 1. Its
  (log x)^{−A} factor is discarded, so Siegel's theorem is never needed.
  - **Every bound holds uniformly for N ≥ N₀.** N₀ is effective provided Lemme 5's constants
    η, c, δ are. Drappeau does not state this explicitly; we read Harper's source as giving it.
  - The algorithm uses no analytic constant except through thmA's "N′ < N₀" branch. That branch can
    be removed (§7).

**What remains external (stated precisely, §2).**
- The published input is Drappeau's Lemme 5. In both published sources its proof is a reference to
  the computations of §3.3 of Harper's *unpublished* preprint 1208.5992. Drappeau writes "il découle
  des calculs de la section 3.3 de [Har12b]"; Pascadi writes "follows from the work of Harper in
  [21, Section 3.3]".
- So the chain now rests on a refereed, published, numbered statement, as the brief allows. Our
  argument nowhere enters the preprint, but the preprint is still behind that statement. We could not
  remove it without re-proving Harper's Theorem 3 and Proposition 1. That is a zero-density and
  zero-free-region analysis of Ψ(x, y; χ), with Page's theorem, and we did not attempt it.
- **Harper 1208.5992 was never published.** zbMATH (queried 2026-10-20) lists it as "Preprint,
  arXiv:1208.5992 (2012)", and it is cited as a preprint in 2019–2025 papers.

**Reviews (§10).**
- GPT-6.1-sol adversarial referee: ACCEPT WITH FIXES (4 MINOR).
- Opus proof critic: CORRECT WITH MINOR FIXES.
- Opus citation critic: every citation confirmed, or confirmed after a correction.
- None found anything FATAL or MAJOR, and all fixes are applied. All three name (D) as the weakest
  remaining link.

Everything else cited is a classical published statement, checked against the source where marked:
Hildebrand–Tenenbaum 1986, Drappeau 2015 Lemme 1(i) (from de la Bretèche–Tenenbaum 2005), Iwaniec–Kowalski
Thm 7.13 and Rosser–Schoenfeld 1962.

---

## 1. Statement

Notation.
- P⁺(n) and P⁻(n) are the largest and smallest prime factors of n, with P⁺(1) = 1 and P⁻(1) = +∞.
- S(z, y) is the set of y-smooth n ≤ z, Ψ(z, y) = #S(z, y), and ϱ(z, y) := Ψ(z, y)/z.
- For a Dirichlet character χ of any modulus, Ψ(z, y; χ) := Σ_{n∈S(z,y)} χ(n). Here χ(n) = 0 when
  (n, mod χ) > 1.
- For (a, q) = 1: Ψ(z, y; q, a) := #{n ∈ S(z, y) : n ≡ a (mod q)}, Ψ_q(z, y) := #{n ∈ S(z, y) :
  (n, q) = 1}, and E(z, y; q, a) := Ψ(z, y; q, a) − Ψ_q(z, y)/φ(q).
- u := log x/log y, and H(u) := exp(u/log²(u + 1)) (Drappeau's normalisation).
- α(z, y) is the Hildebrand–Tenenbaum saddle point: the solution of Σ_{p≤y} log p/(p^α − 1) = log z.
- L = L_N[1/2, 1] = exp((log N log log N)^{1/2}).

**Regime (R).** Fix a > 0, a compact Ξ ⊂ (0, ∞) with ξ_max := max Ξ, and λ ∈ (0, 1/4).
- N → ∞, y := ⌊L_N[1/2, a]⌋ (this is B), and x is real with log x/log N ∈ Ξ.
- K := 6 + ⌈2ξ_max/a²⌉ and R₁ := y^K.
- ε_N := exp(−c₀ (log N)^{1/2}/(log log N)^{5/2}). Here c₀ := δ_D·min Ξ/(2a) > 0, and δ_D is the
  constant of Drappeau's Lemme 5 (§2, D).

**Lemma 8.5′.** In regime (R), let Y, Q be reals with y^K ≤ Y and x^λ·Y ≤ Q ≤ x^{1/2−λ}. Put
- J := {ℓ prime : Y < ℓ ≤ 2Y};
- 𝒯_Q := {t = ℓs ∈ (Q/2, Q] : ℓ ∈ J, P⁺(s) ≤ y};
- δ_𝒯 := ϱ(Q/Y, y)/log Y.

Then there are N₀, C > 0, depending only on (a, Ξ, λ) and on the constants of the results in §2, all
effective, such that for N ≥ N₀:
- (i) Qδ_𝒯/10 ≤ #𝒯_Q ≤ 3Qδ_𝒯;
- (ii) Σ_{t∈𝒯_Q} max_{(r,t)=1} |Ψ(x, y; t, r) − Ψ_t(x, y)/φ(t)| ≤ C·Ψ(x, y)·δ_𝒯·ε_N.

*Remarks on the statement.*
- This contains pcest §2.2's Lemma 8.5′ (regime (B), which is what Theorem A uses), with ε_N in
  place of H(u)^{−c} + B^{−c} + x^{−η₂}.
  - In regime (R), H(u)^{−c} = exp(−c′(log N)^{1/2}/(log log N)^{5/2}(1 + o(1))), so the two forms of
    ε_N are of the same type.
  - FM Step 3 (pcest §2.3) needs only (log ω)·ε_N·log log N → 0 with log ω = O((log N log log N)^{1/2}).
    This holds.
- The condition y^K ≤ Y also covers hypD's regime (A), Y = exp((log N)^{3/4}). Regime (A) is not used
  by Theorem A.
- Interface with the final design (thmA_proof §1):
  - Y = 2^{j−1} = N^{2ν+o(1)}, Q = N^{μ+ν+o(1)}, Q/Y = N^{μ−ν+o(1)} with μ − ν ≥ η₀, and
    x = N^{ξ+o(1)} with ξ = μ + h + dν. Also ξ/2 − (μ + ν) ≥ η₀/2 (the FM level-½ floor).
  - So the hypotheses hold with λ := η₀/(4ξ_max), Ξ a small closed neighbourhood of ξ, and a = θ/√2.
    They hold uniformly over FM's grid 𝒳 and its O(1) dyadic Q, because ω = L^{C₁} changes log x
    by o(log N).

---

## 2. Cited results (exact statements; only these are used)

All are published. "Checked" means we read the statement in the source, or in an author copy whose
numbering matches the published version.

- **(D) Drappeau, *Théorèmes de type Fouvry–Iwaniec pour les entiers friables*, Compos. Math. 151
  (2015), no. 5, 828–862, Lemme 5.**
  - Statement: there exist constants η_D, c_D, δ_D > 0 such that whenever (log x)^{c_D} ≤ y ≤ x and
    Q ≤ x, for every fixed A > 0,
      Σ_{q≤Q} φ(q)^{−1} Σ_{χ (mod q), 1<cond(χ)≤x^{η_D}} |Σ_{n∈S(x,y)} χ(n)|
        ≪_A Ψ(x, y){H(u)^{−δ_D}(log x)^{−A} + y^{−δ_D}}.
  - "La constante implicite est effective si A < 1."
  - Checked in arXiv 1307.7554v3 and in the author's copy (bffi-friable.pdf, numbering identical).
  - The same statement appears as **Pascadi, *Smooth numbers in arithmetic progressions to large
    moduli*, Compos. Math. 161 (2025), no. 8, 1923–1974, Lemma 3.2**. Pascadi gives the range
    (log x)^C ≤ y ≤ x, Q ≤ x, conductor ≤ x^ε, and "effective if A < 1". Checked in arXiv
    2304.11696v3. Pascadi cites it as "[7, Lemme 5]", where [7] is the Compositio version.
  - **Provenance.** Drappeau writes only "il découle des calculs de la section 3.3 de [Har12b]"; Pascadi
    writes "follows from the work of Harper in [21, Section 3.3]".
    - *Our reading*, not either author's statement: the relevant computation is Harper's eq. (3.1),
      the effective bound Σ_{1<r≤x^η} Σ*_{χ*} Σ_{q≤Q} φ(q)^{−1} Σ_{χ induced by χ*} |Ψ(x, y; χ)| ≪
      Ψ(x, y)(e^{−cu/log²(u+1)} + y^{−c}), with an ineffective (log x)^{−A} strengthening asserted after it.
    - Harper's literal context is that of his Theorems 1–2 (Q ≤ √Ψ(x, y), log^K x ≤ y). Drappeau's
      "Q ≤ x" is a legitimate reading of the computation, since it uses only Σ_{s≤Q/r} 1/φ(s) ≪ log x.
      The "effective if A < 1" form is Drappeau's addition.
    - **Our use lies inside even Harper's literal context.** We take Q = R₁ = y^K ≤ √Ψ(x/d, y), and we
      discard the (log x)^{−A} factor entirely, so Harper's effective form (3.1) is all that is needed.
      We take **A = 1/2** only to quote Drappeau's effective statement.
- **(dBT) Drappeau 2015 (as in (D)), Lemme 1(i).** (α is the saddle point. The citation critic
  confirmed this from the proof of dBT's (2·19), via their Lemme 5.1.)
  - Statement: "Lorsque 2 ⩽ y ⩽ x et 1 ⩽ d ⩽ x, on a Ψ(x/d, y) ≪ Ψ(x, y)/d^α", where α = α(x, y) is
    the saddle point of Drappeau's (1.6) (= HT86 (3.1)). Checked in arXiv v3 and the author's copy.
  - Drappeau derives it from de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), 139–202,
    Théorème 2.4(i) with m = 1. For m = 1, ω(m) ≪ 1 holds, which allows the full range 1 ≤ d ≤ x.
    - The dBT author copy we read says it includes corrections relative to the published version
      ("Nous incluons ici certaines corrections par rapport à la version publiée"). So we anchor the
      citation to Drappeau's published, numbered Lemme 1(i), and give dBT only as its source.
  - We write C_dBT for the implied constant (absolute).
- **(HT) Hildebrand–Tenenbaum, *On integers free of large prime factors*, Trans. AMS 296 (1986),
  265–290.** Checked (scan).
  - Theorem 2(i), eq. (2.4): α(x, y) = (log(1 + y/log x)/log y)·(1 + O(log log(1 + y)/log y)),
    uniformly for x ≥ y ≥ 2.
  - Theorem 3: Ψ(cx, y) = Ψ(x, y)c^{α(x,y)}(1 + O(1/u + log y/y)), uniformly for x ≥ y ≥ 2 and
    1 ≤ c ≤ y.
- **(LS) Iwaniec–Kowalski, *Analytic Number Theory* (AMS Colloq. Publ. 53, 2004), Theorem 7.13.**
  For complex a_n with M < n ≤ M + N:
    Σ_{q≤Q} (q/φ(q)) Σ*_{χ (mod q)} |Σ_n a_n χ(n)|² ≤ (Q² + N − 1) Σ|a_n|²,
  where Σ* is over primitive characters. Quoted verbatim as Drappeau's Lemme 6, and as Pascadi's
  Lemma 3.3 in arXiv v3, which is **Lemma 3.4 in the published version** (codex lem85_1, finding 12).
  We apply it with Q = ⌊2R⌋, an integer.
- **(P) Truncated Perron kernel — proved here, not cited.** For c > 0, T > 0 and w > 0 with w ≠ 1:
    |(2πi)^{−1}∫_{c−iT}^{c+iT} w^s ds/s − 1_{w>1}| ≤ w^c/(πT|log w|).
  *Proof.*
  - The integrand w^s/s is meromorphic, with a single simple pole at s = 0 of residue 1.
  - If w > 1, close the segment [c − iT, c + iT] by the rectangle with vertical side Re s = −σ₀,
    σ₀ → ∞. The left side contributes ≤ w^{−σ₀}·2T/σ₀ → 0. Each horizontal side contributes
    ≤ (2π)^{−1}∫_{−∞}^{c} w^σ dσ/T = w^c/(2πT log w). The residue theorem gives the claim.
  - If w < 1, close to the right instead. There is no pole, and the same bound holds with |log w|. ∎
  - (This is the standard lemma of Davenport, *Multiplicative Number Theory*, ch. 17, and of
    Montgomery–Vaughan, *Multiplicative Number Theory I*, Cor. 5.3. We prove it so that nothing
    unnumbered is cited.)
- **(RS) Rosser–Schoenfeld, Illinois J. Math. 6 (1962), 64–94.** Numbering as commonly cited, not
  re-checked; any Chebyshev-strength substitute changes only constants.
  - Corollary 1: π(z) < 1.25506 z/log z for z > 1.
  - Corollary 3: π(2z) − π(z) > 3z/(5 log z) for z ≥ 20.5.
  - Theorem 15: n/φ(n) < e^γ log log n + 2.50637/log log n for n ≥ 3, except n = 223092870, which
    needs 2.51. We use the 2.51 form. These statements were confirmed from secondary sources; the
    "Corollary 3" label was not.

Consequences used, valid in regime (R) for N ≥ N₀. Write β(z, y) := α(z, y), the saddle point.
- **(α1)** 3/4 ≤ β(z, y) ≤ 1 for 2y ≤ z ≤ x.
- **(α2)** 1 − β(z, y) ≤ (log(u log y) + C_α log log y)/log y for 2y ≤ z ≤ x, with C_α absolute.
- Proof of (α2): by (HT) Theorem 2(i), α ≥ ᾱ(1 − C loglog y/log y) with ᾱ := log(1 + y/log z)/log y.
  Also log(1 + y/log z) ≥ log y − log log z, and log z ≤ log x = u log y.
- (α1) follows from (α2) in regime (R). The upper bound α ≤ 1 holds because y > log z, so that
  Σ_{p≤y} log p/(p − 1) > log z (Mertens) and the defining equation forces α < 1. Only α ≥ 3/4 is
  actually used.

---

## 3. Elementary lemmas (all proved here)

**E1 (orthogonality majorant).** For (a, q) = 1, E(z, y; q, a) = φ(q)^{−1} Σ_{χ≠χ₀ (mod q)} χ̄(a)Ψ(z, y; χ).
Hence max_{(a,q)=1} |E(z, y; q, a)| ≤ Δ(q) := φ(q)^{−1} Σ_{χ≠χ₀} |Ψ(z, y; χ)|.
*Proof.* Orthogonality, together with the fact that Ψ(z, y; χ₀) = Ψ_q(z, y). ∎

**E2 (lift averaging).** Let ℓ be a prime with ℓ > y, let s ≥ 1 with ℓ ∤ s, and put t = ℓs. For every
r with (r, t) = 1:
  E(z, y; t, r) = E(z, y; s, r)/(ℓ − 1) + Rem(t, r),
  where |Rem(t, r)| ≤ Δ_ℓ(t) := φ(t)^{−1} Σ_{χ (mod t), ℓ | cond χ} |Ψ(z, y; χ)|.

*Proof.* (ℤ/t)^* = (ℤ/ℓ)^* × (ℤ/s)^*, and every χ (mod t) factors as χ_ℓχ_s. Also ℓ | cond χ if and
only if χ_ℓ ≠ χ₀^{(ℓ)}.
- Every y-smooth n is coprime to ℓ. So for χ_ℓ = χ₀^{(ℓ)} we get Ψ(z, y; χ₀^{(ℓ)}χ_s) = Ψ(z, y; χ_s).
  Also Ψ_t(z, y) = Ψ_s(z, y), and φ(t) = (ℓ − 1)φ(s).
- In E1 for modulus t, the terms with χ_ℓ = χ₀^{(ℓ)} and χ_s ≠ χ₀^{(s)} sum to
  ((ℓ−1)φ(s))^{−1} Σ_{χ_s≠χ₀} χ̄_s(r)Ψ(z, y; χ_s) = E(z, y; s, r)/(ℓ − 1).
- The remaining terms have ℓ | cond χ.

Equivalently, without characters: averaging E(z, y; t, ·) over the ℓ − 1 lifts of r (mod s) to
(ℤ/t)^* gives exactly E(z, y; s, r)/(ℓ − 1). This is checked exactly in `tests/test_lem85.py`. ∎

**E3 (coprimality by Möbius).** Let χ (mod q) be induced by a primitive χ* (mod r), and put
q_y := ∏_{p|q, p≤y} p. Then
  Ψ(z, y; χ) = Σ_{d | q_y, (d,r)=1} μ(d)χ*(d)Ψ(z/d, y; χ*).

*Proof.* For y-smooth n, (n, q) = 1 if and only if (n, q_y) = 1. Expand 1_{(n,q_y)=1} = Σ_{d|(n,q_y)} μ(d).
Then n = dn′ with n′ ∈ S(z/d, y), because d is y-smooth. Terms with (d, r) > 1 vanish. Checked
exactly for real characters in the tests. ∎

**E4 (factorisation).** Let M₁ ≥ 1. Every n ∈ S(z, y) with n > M₁ has exactly one factorisation n = mk
with:
- M₁ < m ≤ P⁺(m)M₁ (that is, m/P⁺(m) ≤ M₁ < m);
- P⁺(m) ≤ P⁻(k);
- k ∈ S(z, y).

Moreover k ≤ z/M₁.

*Proof.* Write n = p₁p₂⋯p_J with p₁ ≤ ⋯ ≤ p_J, and let m = p₁⋯p_j with j minimal such that the product
exceeds M₁. Conversely, if (m, k) satisfies the three conditions, the primes of m come first in the
sorted list of n's primes (ties at P⁺(m) = P⁻(k) are harmless). Then m/P⁺(m) ≤ M₁ < m forces j to be
the number of prime factors of m. Checked exhaustively on toy ranges in the tests. ∎

**E5.** Σ_{j≤X} 1/φ(j) ≤ 2(1 + log X) for X ≥ 1.
*Proof.* j/φ(j) = Σ_{e|j} μ²(e)/φ(e). So the sum is at most (1 + log X)·Σ_e μ²(e)/(eφ(e)), and
Σ_e μ²(e)/(eφ(e)) = ζ(2)ζ(3)/ζ(6) < 1.95. ∎

**E6 (constants).** Over squarefree d:
- Σ_d 1/(dφ(d)) < 1.95;
- Σ_d 1/(φ(d)d^{1/2}) ≤ ∏_p (1 + 2p^{−3/2}) ≤ exp(2Σ_p p^{−3/2}) < e^{1.7} < 6, since Σ_p p^{−3/2} < 0.85;
- if β, β′ ≥ 3/4, then Σ_d d^{1−β−β′}/φ(d) ≤ ∏_p (1 + p^{−1/2}/(p − 1)) < 6.

Also n/φ(n) ≤ 5 log log(n + 16) for all n ≥ 1, by (RS) Theorem 15 for n ≥ 16 and by inspection for
n < 16.

**E7 (density lower bound).** For z ≥ 2y with u_z := log z/log y ≤ u and N ≥ N₀:
  ϱ(z, y) ≥ y^{−1}exp(−u log(4u log y)).
In regime (R), u log(4u log y) = (ξ/a²)(1 + o(1)) log y. Hence ϱ(z, y) ≥ y^{−1−ξ_max/a²−o(1)}.

*Proof.* Let k := ⌊u_z⌋ and 𝒫 := primes in (y/2, y]. By (RS) Cor. 3, #𝒫 ≥ 3y/(10 log y) ≥ k. The
products of k distinct primes of 𝒫 are distinct, y-smooth, and ≤ y^k ≤ z. So
  Ψ(z, y) ≥ C(#𝒫, k) ≥ (#𝒫/k)^k ≥ (z/y)(3/(10k log y))^k.
Then use k ≤ u. In regime (R), log y = a(log N log log N)^{1/2}(1 + o(1)) and
u = (ξ/a)(log N/log log N)^{1/2}(1 + o(1)). ∎

**E8 (weights over smooth moduli).** Let 2y ≤ Q′ ≤ x, r, d ≥ 1 with (r, d) = 1, and β′ := β(Q′, y).
Then
  W(r, d) := Σ_{Q′/2<s≤Q′, P⁺(s)≤y, rd|s} 1/φ(s) ≤ 10 C_dBT log log(Q′ + 16)·ϱ(Q′, y)·(rd)^{1−β′}/(φ(r)φ(d)).

*Proof.*
- Write s = rdj, so that 1/φ(s) ≤ 1/(φ(rd)φ(j)) and φ(rd) = φ(r)φ(d). Here j is y-smooth and
  j ∈ (J/2, J] with J := Q′/(rd); the sum is empty if J < 1.
- By E6, Σ_j 1/φ(j) ≤ 5 log log(Q′ + 16)·(2/J)Ψ(J, y) = 10 log log(Q′ + 16)·ϱ(J, y).
- By (dBT) at (Q′, y) with d := rd: ϱ(J, y) = Ψ(Q′/(rd), y)·rd/Q′ ≤ C_dBT ϱ(Q′, y)(rd)^{1−β′}. ∎

---

## 4. Lemma LS (large conductors, from the classical large sieve)

**Lemma LS.** Let y ≥ 2, z ≥ max(4, y^{3/2}) and R ≥ 1. Then
  Σ_{R<r≤2R} (r/φ(r)) Σ*_{χ (mod r)} |Ψ(z, y; χ)| ≤ 30π(y) log(2zR)·(y^{3/2}z + y^{1/2}z^{1/2}R²).

*Proof.* Write w_χ := r/φ(r) for χ primitive of conductor r ∈ (R, 2R]. Then Σ_χ w_χ ≤ Σ_{r≤2R} r ≤ 4R².
Put M₁ := z^{1/2}. By E4,
  Ψ(z, y; χ) = Σ_{n≤M₁, n∈S(y)} χ(n) + Σ_{p≤y} Σ_{m,k: mk≤z} a^{(p)}_m b^{(p)}_k χ(m)χ(k),
where:
- a^{(p)}_m := 1 if P⁺(m) = p and M₁ < m ≤ pM₁, and 0 otherwise;
- b^{(p)}_k := 1 if P⁻(k) ≥ p, P⁺(k) ≤ y and k ≤ z/M₁, and 0 otherwise.

Neither sequence depends on χ, and #supp a^{(p)} ≤ yM₁, #supp b^{(p)} ≤ z/M₁.

*Short part.* By Cauchy–Schwarz and (LS) with Q = 2R:
  Σ_χ w_χ |Σ_{n≤M₁}| ≤ (4R²)^{1/2}((4R² + M₁)M₁)^{1/2} ≤ 2Rz^{1/2} + 4R²z^{1/4}.

*Bilinear part.*
- Let z′ := ⌊z⌋ + 1/2, c := 1/log(2z) and T := (2zR)^4. For integers m, k ≥ 1, 1_{mk≤z} = 1_{z′/(mk)>1}.
- On the support, mk ≤ (yM₁)(z/M₁) = yz and |z′ − mk| ≥ 1/2. Hence |log(z′/mk)| ≥ |z′−mk|/max(z′, mk)
  ≥ 1/(4yz), and (z′/mk)^c ≤ (2z)^c = e.
- By (P), each pair (m, k) contributes an error ≤ 4eyz/(πT) < 4eyz/T. There are ≤ yz pairs per p, at most y
  primes p, and total character weight ≤ 4R². So the total error over all χ and p is
  < 16e y³z²R²/T ≤ e·y³/(z²R²) ≤ 3.
- The main term is (2πi)^{−1}∫_{c−iT}^{c+iT} A_p(χ, s)B_p(χ, s) z′^s ds/s, where
  A_p(χ, s) = Σ_m a^{(p)}_m χ(m)m^{−s} and B_p is defined similarly.
- For each fixed t, Cauchy–Schwarz and (LS) with coefficients a_m m^{−c−it} (of modulus ≤ 1) give
    Σ_χ w_χ |A_p(χ, c+it)B_p(χ, c+it)| ≤ (yM₁ + 4R²)^{1/2}(z/M₁ + 4R²)^{1/2}(yM₁)^{1/2}(z/M₁)^{1/2}
      ≤ (yz^{1/2} + 4R²)(yz)^{1/2} = y^{3/2}z + 4y^{1/2}z^{1/2}R².
- Also ∫_{−T}^{T} dt/|c + it| ≤ 2/c + 2 log T ≤ 10 log(2zR) and z′^c ≤ e. So the bilinear part is
  ≤ 5π(y) log(2zR)(y^{3/2}z + 4y^{1/2}z^{1/2}R²) + 3.

Adding the two parts, and using 2Rz^{1/2} + 4R²z^{1/4} + 3 ≤ 9y^{1/2}z^{1/2}R² + y^{3/2}z, gives the
bound. ∎

**Corollary LS′ (moduli with smooth cofactor, any sub-family).** Let z₀ ≥ max(4, y^{3/2}), and let
1 ≤ R₀ ≤ Q₀ with z₀/Q₀ ≥ max(4, y^{3/2}). For each modulus q ≤ Q₀ let 𝒳(q) be any set of characters
mod q whose conductors lie in (R₀, Q₀]. Then
  Σ_{q≤Q₀} φ(q)^{−1} Σ_{χ∈𝒳(q)} |Ψ(z₀, y; χ)| ≤ 1.5·10³ π(y)(log z₀)²·(y^{3/2}z₀/R₀ + y^{1/2}z₀^{1/2}Q₀).

*Proof.*
- Group (q, χ) by the inducing primitive χ* (mod r) and apply E3: |Ψ(z₀, y; χ)| ≤ Σ_{d|q_y, (d,r)=1}
  |Ψ(z₀/d, y; χ*)|, with d squarefree and d ≤ Q₀.
- Swap the sums. By E5 and φ(rdj) ≥ φ(r)φ(d)φ(j), Σ_{q≤Q₀, rd|q} 1/φ(q) ≤ 2(1 + log Q₀)/(φ(r)φ(d)).
- On each dyadic block r ∈ (R, 2R] with R₀ ≤ R ≤ Q₀, use 1/φ(r) ≤ (r/φ(r))/R and Lemma LS at
  z = z₀/d ≥ z₀/Q₀. This gives 30π(y)·2 log z₀·(y^{3/2}z/R + y^{1/2}z^{1/2}R).
  (Here log(2zR) ≤ log(2z₀Q₀) ≤ 2 log z₀, because Q₀ ≤ z₀/4.)
- Sum over the blocks R = 2^iR₀ (Σ1/R ≤ 2/R₀, ΣR ≤ 2Q₀), then over d with weights 1/φ(d), using
  Σ_d 1/(dφ(d)) < 2 and Σ_d 1/(φ(d)d^{1/2}) < 6 (E6). This gives
    ≤ 2(1 + log Q₀)·60π(y) log z₀·(4y^{3/2}z₀/R₀ + 12y^{1/2}z₀^{1/2}Q₀),
  and 1 + log Q₀ ≤ log z₀ because eQ₀ ≤ z₀. ∎

---

## 5. Lemma SM (smooth moduli, relative saving)

**Lemma SM.** In regime (R), let 2y ≤ Q′ ≤ x^{1/2−λ}, and let 𝒮 := {s ∈ (Q′/2, Q′] : P⁺(s) ≤ y}. Then
for N ≥ N₀:
  Σ_{s∈𝒮} max_{(a,s)=1} |E(x, y; s, a)| ≤ Σ_{s∈𝒮} Δ(s) ≤ C_SM·Ψ(x, y)·ϱ(Q′, y)·ε_N.

*Proof.* The first inequality is E1. Split each Δ(s) according to the conductor r of χ: case (a) is
1 < r ≤ R₁ = y^K, and case (b) is r > R₁.

**(a) Small conductors.**
- Apply E3 (s is y-smooth, so s_y = rad s) and swap the sums. W(r, d) = 0 unless rd ≤ Q′, so only
  d ≤ Q′ occurs. With W(r, d) from E8:
    (a) ≤ Σ_{1<r≤R₁} Σ*_{χ* (mod r)} Σ_{d≤Q′ sqfree, P⁺(d)≤y, (d,r)=1} |Ψ(x/d, y; χ*)|·W(r, d)
        ≤ 10 C_dBT log log(Q′+16)·ϱ(Q′, y)·R₁^{1−β′}·Σ_{d≤Q′} (d^{1−β′}/φ(d))·D(x/d),
    where D(z) := Σ_{1<r≤R₁} φ(r)^{−1} Σ*_{χ (mod r)} |Ψ(z, y; χ)|.
  The cutoff d ≤ Q′ is kept through the next two steps, (D) and (dBT). Only the final Euler sum is
  extended to all squarefree d, by positivity. (Fix from codex lem85_1 #2 and the Opus proof critic #5.)
- *Apply (D) to the primitive characters only.* Take q = r ≤ R₁ and χ = χ* primitive. These pairs
  are among the terms of Lemme 5's left side with Q := R₁. Take A = 1/2 and the argument z = x/d.
  - The hypotheses hold: x/d ≥ x/Q′ ≥ x^{1/2}, so (log z)^{c_D} ≤ y ≤ z and R₁ = y^K =
    N^{o(1)} ≤ z^{η_D} for N ≥ N₀.
  - Hence D(x/d) ≤ C_D Ψ(x/d, y)(H(u_d)^{−δ_D} + y^{−δ_D}), where u_d := log(x/d)/log y ≥ u/2
    (because d ≤ Q′ ≤ x^{1/2}).
  - H is increasing on [4, ∞), so H(u_d) ≥ H(u/2).
- By (dBT) (d ≤ Q′ ≤ x), Ψ(x/d, y) ≤ C_dBT Ψ(x, y)d^{−β}, with β := β(x, y). By (α1) and E6,
  extending to all squarefree d, Σ_{d≤Q′} d^{1−β′−β}/φ(d) < 6.
- By (α2), R₁^{1−β′} = exp(K(1 − β′) log y) ≤ (u log y)^K (log y)^{KC_α} ≤ (log x)^{K(1+C_α)}.
- Therefore
    (a) ≤ 60 C_dBT² C_D·log log x·(log x)^{K(1+C_α)}·(e^{−δ_D (u/2)/log²(u/2+1)} + y^{−δ_D})·Ψ(x, y)ϱ(Q′, y).
- In regime (R), u = (ξ/a)(log N/log log N)^{1/2}(1 + o(1)) and log²(u + 1) = (1/4)(log log N)²(1 + o(1)).
  So δ_D(u/2)/log²(u/2+1) = (2δ_Dξ/a)(log N)^{1/2}/(log log N)^{5/2}·(1 + o(1)).
  - This beats the (log x)^{O(1)} = exp(O(log log N)) factor and the threshold 2c₀ = δ_D min Ξ/a.
  - y^{−δ_D} = exp(−δ_D a(log N log log N)^{1/2}) is smaller still.
- Hence (a) ≤ Ψ(x, y)ϱ(Q′, y)ε_N for N ≥ N₀.

**(b) Large conductors.**
- If Q′ ≤ R₁, case (b) is empty: every conductor divides s ≤ Q′. (Fix from codex #3 and Opus #6.)
  Otherwise apply Corollary LS′ with z₀ = x, Q₀ = Q′, R₀ = R₁ and 𝒳(s) = {χ (mod s) : cond χ > R₁}, keeping
  only the moduli s ∈ 𝒮 (positivity). This gives
    (b) ≤ 1.5·10³π(y)(log x)²(y^{3/2}x/y^K + y^{1/2}x^{1/2}Q′).
- By E7, Ψ(x, y)ϱ(Q′, y) ≥ x·y^{−2−2ξ_max/a²−o(1)}.
  - The first term is ≤ Ψϱ(Q′)·y^{5/2−K+2+2ξ_max/a²+o(1)} ≤ Ψϱ(Q′)y^{−1}, because
    K = 6 + ⌈2ξ_max/a²⌉.
  - The second is ≤ Ψϱ(Q′)·x^{−λ}y^{O(1)}, because Q′ ≤ x^{1/2−λ}.
  - Both are ≤ Ψϱ(Q′)ε_N for N ≥ N₀. ∎

*Why this works exactly at y = L[1/2].*
- The large sieve needs R₁ ≥ 1/(ϱ(x)ϱ(Q′)) up to y^{O(1)}, that is, log R₁ ≳ u log u.
- The small-conductor transfer costs R₁^{1−β} = exp(log R₁·log(u log y)/log y).
- Both are compatible with log R₁ ≍ log y precisely when u log u ≍ log y, i.e. y = L^{Θ(1)}.
- For much smaller y one would need the medium-conductor zero-density analysis of Harper/LV; we
  avoid it.

---

## 6. Proof of Lemma 8.5′

*(i) Counting 𝒯_Q.*
- Since ℓ > y ≥ P⁺(s), each t ∈ 𝒯_Q determines (ℓ, s), and ℓ ∤ s. So
  #𝒯_Q = Σ_{ℓ∈J} (Ψ(Q/ℓ, y) − Ψ(Q/2ℓ, y)).
- Upper bound: #𝒯_Q ≤ #J·Ψ(Q/Y, y) ≤ (1.256·2Y/log 2Y)(Q/Y)ϱ(Q/Y, y) ≤ 3Qδ_𝒯, by (RS) Cor. 1.
- Lower bound:
  - Put z := Q/(2ℓ) ≥ x^λ/4. Then u_z → ∞ and α(z, y) → 1, by (HT) Theorem 2(i).
  - (HT) Theorem 3 with c = 2 gives Ψ(2z, y) − Ψ(z, y) ≥ 0.8Ψ(z, y). With c = 4 it gives
    Ψ(Q/Y, y) ≤ 4.1Ψ(Q/4Y, y) ≤ 4.1Ψ(z, y).
  - So each ℓ contributes ≥ 0.19(Q/Y)ϱ(Q/Y, y).
  - #J ≥ 3Y/(5 log Y), by (RS) Cor. 3. Hence #𝒯_Q ≥ 0.11Qδ_𝒯.

*(ii) The error sum.* By E2, for t = ℓs ∈ 𝒯_Q,
  max_{(r,t)=1} |E(x, y; t, r)| ≤ max_{(a,s)=1} |E(x, y; s, a)|/(ℓ − 1) + Δ_ℓ(t).

**Part 1 (χ_ℓ principal).**
- Σ_t max_a |E(x, y; s, a)|/(ℓ − 1) = Σ_{ℓ∈J} (ℓ − 1)^{−1} Σ_{s∈(Q/2ℓ, Q/ℓ], P⁺(s)≤y} max_a |E(x, y; s, a)|.
- Apply Lemma SM with Q′ = Q/ℓ. Its hypotheses hold: Q/ℓ ≥ Q/(2Y) ≥ x^λ/2 ≥ 2y and
  Q/ℓ < Q ≤ x^{1/2−λ}.
- This gives ≤ C_SM Ψ(x, y)ε_N Σ_ℓ ϱ(Q/ℓ, y)/(ℓ − 1).
- Ψ(·, y) is nondecreasing and ℓ ∈ (Y, 2Y], so ϱ(Q/ℓ, y) = Ψ(Q/ℓ, y)ℓ/Q ≤ Ψ(Q/Y, y)·2Y/Q =
  2ϱ(Q/Y, y). Also Σ_{ℓ∈J} 1/(ℓ − 1) ≤ #J/Y ≤ 2.6/log Y, by (RS) Cor. 1.
- So **Part 1 ≤ 5.2 C_SM·Ψ(x, y)δ_𝒯ε_N.**

**Part 2 (ℓ | conductor).**
- Every χ counted in Δ_ℓ(t) has conductor r with ℓ | r, so r > Y.
- Apply Corollary LS′ with z₀ = x, Q₀ = Q, R₀ = Y and 𝒳(t) = {χ : ℓ_t | cond χ} for t ∈ 𝒯_Q (empty
  for other q). This gives Part 2 ≤ 1.5·10³π(y)(log x)²(y^{3/2}x/Y + y^{1/2}x^{1/2}Q).
- By E7, Ψ(x, y)δ_𝒯 ≥ x·y^{−2−2ξ_max/a²−o(1)}/log Y. With Y ≥ y^K and Q ≤ x^{1/2−λ} this gives
  **Part 2 ≤ Ψ(x, y)δ_𝒯·y^{−1/2} ≤ Ψ(x, y)δ_𝒯ε_N** for N ≥ N₀.

Adding Parts 1 and 2 proves (ii) with C := 5.2 C_SM + 1. ∎

---

## 7. Exceptional zeros, effectivity, uniformity in N

- **Where a Siegel or Page zero could enter.** Only in case (a) of Lemma SM, the primitive characters
  of conductor ≤ y^K. Lemma LS, Corollary LS′ and E1–E8 use no information about L-function zeros.
- **How it is handled.** Through (D) with A = 1/2. Drappeau states that its implied constant is
  effective for A < 1.
  - We describe the source computation (Harper §3.3) here; this is not a step of our proof.
    - By Page's theorem, at most one exceptional real primitive character χ_bad has small conductor.
    - For it, Harper bounds Ψ(x, y; χ_bad) by an Euler-product (saddle-point) argument. For
      |t| < 1/(2 log y), that argument uses the favourable sign of the exceptional-zero term in
      Σ Λ(n)χ(n). It uses **no lower bound on 1 − β**, which is why it is effective. It yields the
      e^{−cu/log²(u+1)} saving.
    - In regime (R), u ≫ (log log x)², so only this Euler-product branch is relevant. Harper's other
      branch, y > x^{1/(log log x)²}, rests on Fouvry–Tenenbaum and an unwritten check, and is not
      needed.
  - Siegel's ineffective theorem enters only the (log x)^{−A} strengthening, which we discard.
- **Uniformity.** Every bound above holds for all N ≥ N₀(a, Ξ, λ), uniformly in x, Y, Q in the stated
  ranges. There is no exceptional set of N and no exceptional set of moduli: the conductor of a
  potential exceptional character is absorbed in (a) like any other. The dependence on the
  constants of (D), (dBT) and (HT) is explicit.
- **Effectivity.**
  - Drappeau states effectivity only for the *implied constant* of (D) when A < 1. He does not state
    it for η_D, c_D, δ_D.
  - In Harper's source, c and K are "fixed and effective constants" (his Theorem 1), and η, δ are
    absolute constants fixed explicitly in §3.3 from Theorem 3, Proposition 1 and the log-free
    zero-density estimate.
  - So all our constants are effective, provided (D)'s η_D, c_D, δ_D are. We read Harper as giving
    this, but neither published source states it. (Correction from the Opus critics.)
  - The constants of (dBT), (HT), (LS), (P) and (RS) are effective. Note that δ_D enters c₀ and hence
    ε_N.
- **Does the algorithm need them?** Only through one threshold, which can be removed.
  - Lemma 8.5′ enters Theorem A only through FM, which is a lower bound on an expectation in the
    running-time analysis.
  - The algorithm's parameters (B, S, M, H, W, J, the grid) are explicit functions of N and of the
    fixed (d, k, μ, η₀). C and c₀ are not used by the algorithm.
  - **However**, thmA_proof §3 (FACTOR) branches on "N′ < N₀", and U3 collects the sources' N₀ into
    that threshold. So the explicit algorithm does use N₀ (Opus proof critic #10).
  - If every N₀ is effective, as we believe (see above), nothing changes. If one were ineffective,
    Theorem A as an existence statement survives, but the algorithm could not be written down
    explicitly.
  - *Robust fix, for the thmA/preprint owners:* run trial division (or any unconditional method) and
    SPLIT in alternation, step by step, with no N₀ branch. This costs at most a factor 2 in expected
    time for large N, and the algorithm then uses no analytic constant at all.

---

## 8. Literature status of the inputs (checked 2026-10-20)

- **Harper, arXiv 1208.5992: unpublished.** zbMATH's author profile (API query) lists it as "Preprint,
  arXiv:1208.5992 [math.NT] (2012)". The following still cite it as a preprint:
  - Drappeau 2015 [Har12b];
  - Pascadi 2025 [21];
  - later papers found by web search (not itemised).

  Harper's published related work is a different paper: J. Number Theory 132 (2012), 182–199, which
  treats individual moduli q ≤ y^{4√e−δ} and does not contain the needed averaged statements. Harper's
  JLMS 2025 BDH paper is about general sequences. **No longer used by us.**
- **Lee–Venkatesan, *Rigorous analysis of a randomised number field sieve*, J. Number Theory 187
  (2018), 92–159** (zbMATH). We did not compare its Lemma 8.5 numbering with the arXiv version.
  **No longer used by us.**
- **What replaced what.**

  | old input (pcest §2.1–2.2) | status | replaced by |
  |---|---|---|
  | majorant W from *inside* LV Lemma 8.5's proof (Part 1) | proof-internal | E2 (exact identity) + Lemma SM |
  | LV Lemma 8.5 statement | published, but proof gap at the top of its y-range (Page step) | Lemma SM (our proof, our regime) |
  | Harper Prop. 2 = LV Fact 8.11 (Part 2, large conductors) | preprint | Lemma LS / LS′ (IK 7.13 + Perron) |
  | Harper Thm 3, Prop. 1, p. 15–16 facts = LV Facts 8.8–8.10, 8.12; Page's theorem | preprint / proof-internal | (D) Drappeau Lemme 5 = Pascadi Lemma 3.2 (published, numbered) |
  | LV Claims 8.14, 8.16 | proof-internal | E6–E8 with (dBT), (HT) |

- **Ruled out as substitutes for (D).**
  - de la Bretèche–Drappeau JEMS 2020 Thm 2.1 carries e^{O(u₂)} (hypD §8).
  - Fouvry–Tenenbaum 1996 needs y ≥ exp((log x)^{2/3+ε}).
  - Harper JNT 2012 is an individual-modulus result with q ≤ y^{4√e−δ}. Granville's Acta Math. 1993
    paper has individual-modulus results for q < √y (saving 1/log y), and an averaged BV result with
    saving (log x)^{−A}, but only for Q ≤ min{√x/log^B x, y^{C(A) log log y/log log log y}}. That Q range is
    far below our Q′, and the saving is not relative to smooth moduli. (Corrected after the citation
    critic.)
  - A key point: *no unrestricted Bombieri–Vinogradov statement over all moduli can give the relative
    saving over smooth moduli*. Its total modulus weight is ≍ log x, while smooth moduli carry weight
    ϱ(Q′). That is why (D) is applied to primitive characters and transferred by E3 and E8, not
    applied to the moduli s directly.

## 9. What a referee should check (our own list)
1. E2: the identity, and the condition ℓ > y with ℓ ∤ s.
2. E3/E8 bookkeeping: (d, r) = 1 makes rd | s with φ(rd) = φ(r)φ(d); the sum over d is over all
   squarefree y-smooth d.
3. (a): applying (D) to the sub-family {q = r, χ primitive}, with Q := R₁ and argument x/d; the range
   conditions; A = 1/2.
4. (α2) and R₁^{1−β′} ≤ (log x)^{O(1)}. This is the step where y = L[1/2] is used.
5. Lemma LS: the Perron constants, and that (LS) is applied with coefficient sequences independent of χ.
6. E7 and the choice of K.
7. Parameter interface with thmA_proof §1 (§1 Remarks).

## 10. Reviews (all three returned "accept with fixes"; every fix applied)

**GPT-6.1-sol adversarial referee (`codex/lem85_1.out`, xhigh): ACCEPT WITH FIXES.**
- No FATAL or MAJOR findings.
- It independently confirmed:
  - the quotation and application of Drappeau's Lemme 5;
  - the transfer weights (no hidden 2^{ω(s)} loss);
  - the absorption of the transfer cost: log H(u/2) has coefficient ≥ 4c₀, against an O(log log N) loss;
  - E4, the Perron and bilinear steps (constants 30 and 1440 ≤ 1500);
  - E7 and the choice of K, with exponent ≤ −3/2;
  - E2, the counting in (i), the FM interface, and effectivity.
- Its four MINOR findings, all applied:

| # | finding | fix |
|---|---|---|
| 2 | §5(a): the displayed d-sum dropped the cutoff d ≤ Q′, which (D) and u_d ≥ u/2 need | cutoff kept through (D) and (dBT); extended only in the final Euler sum |
| 3 | §5(b): LS′ needs R₀ ≤ Q₀, which fails if Q′ < y^K | case (b) is empty then |
| 12 | Pascadi's large-sieve lemma is Lemma 3.4 in print (3.3 in arXiv v3) | corrected |
| 13 | "only published numbered statements" overstated: Davenport's Perron lemma is unnumbered and was unchecked; the dBT author copy contains corrections to the printed version | Perron kernel bound now proved in §2 (P); (dBT) anchored to Drappeau's published Lemme 1(i), with β = α the saddle point |

**Opus proof critic (read-only): CORRECT WITH (MINOR) FIXES.**
- It recomputed every constant numerically, with primes up to 10⁷:
  - ∏(1 + 1/(p(p−1))) = 1.9436;
  - Σp^{−3/2} = 0.84953;
  - ∏(1 + 2p^{−3/2}) = 4.24;
  - ∏(1 + p^{−1/2}/(p−1)) = 3.07;
  - H is increasing from u ≈ 3.92.
- It rechecked the Perron total e·y³/(z²R²), the E7 asymptotics, the margin in (a) (a factor 4 in
  the exponent over ε_N), and the FM interface (#bad ≤ 40C|𝒳|ε_N·#𝒯_Q).
- Findings and fixes:
  - #5 and #6 are the same as codex #2 and #3; applied.
  - #10: thmA's algorithm branches on N′ < N₀. We corrected §7 and proposed a fix for the thmA and
    preprint owners: alternate trial division with SPLIT, with no N₀ branch.
  - #11: Drappeau states effectivity only for the implied constant, not for η, c, δ. §7 and the Verdict
    are corrected. The §7 paraphrase of Harper's exceptional-character step was also corrected.
  - #12: the Verdict said Harper's preprint was "no longer cited at all"; this is reworded. pcest §2
    and thmA_proof §7 still list LV 8.5 and Harper as inputs. Those files belong to other owners and
    are not edited here; they should now cite this note.
- Weakest link named: (D). Next: the (log x)^{O(K)} transfer cost, which is fine but specific to regime (R).

**Opus citation critic (read-only): every citation CONFIRMED, or confirmed after a correction.**
- Confirmed verbatim or exactly: Drappeau Lemme 5 (statement and numbering, also via Pascadi's
  "[7, Lemme 5]"), Lemme 1(i), Lemme 6 and H(u); Pascadi Lemma 3.2 (arXiv v3); dBT Théorème 2.4(i)
  (α is the saddle point, from its proof); HT86 Theorem 2(i) and Theorem 3; the zbMATH entries for LV,
  Drappeau, Pascadi and Harper. Harper's preprint is unpublished: confirmed.
- Corrections, all applied:
  - the "(3.1) and the remark after it" pinpoint is our reading, not the authors';
  - effectivity of η, c, δ is not stated by Drappeau;
  - §7: "independent of the location of the zero" becomes "uses no lower bound on 1 − β (Page plus a
    sign argument)";
  - unidentified "2019/2023 papers" removed;
  - Granville 1993 also has an averaged BV, but only for tiny Q;
  - RS Theorem 15's exact form, and that the "Corollary 3" label is unverified.
- It also noted, and §2 now records, that our use of (D) lies inside even Harper's *literal* context:
  Q = y^K ≤ √Ψ, the (log x)^{−A} factor discarded, and only the Euler-product branch, since
  u ≫ (log log x)².
- Not verified by anyone (books or paywalls): IK Thm 7.13 from the book itself (confirmed through two
  published verbatim quotations); RS numbering; de la Bretèche–Drappeau JEMS 2020 Thm 2.1 (only used
  to rule it out).

## 11. Checks
`.venv/bin/python -m pytest -q tests/test_lem85.py`: exact checks of E2 (lift averaging, 4 parameter
sets, all residues), E3 (Möbius identity, real characters, 4 sets), E4 (uniqueness of n = mk on 3
ranges) and E5. These are correctness checks of elementary identities only.
