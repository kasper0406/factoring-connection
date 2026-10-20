# dlog: rigorous discrete logarithms in F_p^* from the Theorem-A machinery

Agent `dlog`, phase 11, 2026-10-20.
- Files: `src/dlog/{params,pipeline,harness_e1,analyze}.py`, `tests/test_dlog.py`.
- Pre-registration: `notes/dlog_prereg.md` (commit `dfc7f26`).
- Data: `experiments/raw/dlog_e1{.jsonl,_summary.json}`, `experiments/raw/dlog_e2s.jsonl` (supplementary).
- Codex: `codex/dlog_1.{md,out}` (blind), `codex/dlog_2.{md,out}` (adversarial referee).
- Related: `codex/r19_transfer.out` §2, which reached the same route independently.

Labels: **proved** (proof here) · **proved from A** (proved here, conditional on the internal unrefereed relation-supply
chain of Theorem A, Steps 1–7 of `notes/thmA_proof.md`) · **known** (published, with reference) · **conditional** ·
**conjectured** · **experimentally supported**.

## Verdict

In plain terms: yes, the transfer works. Theorem A's relation supply, with the prime p in the role of N, gives an
index-calculus algorithm for discrete logarithms in F_p^*. Its constant is below Pomerance's √2 for every
optimised tuple (θ < 8/9), and below 1 at degree d = 9. As with Theorem A, this rests on the internal, unrefereed supply chain. The DL-specific part is new and
proved here, and it needs no unproved hypothesis.

1. **Theorem DL** (proved from A). Fix any admissible Theorem-A parameter tuple (d, k, μ, η₀) with θ < 1.
   - There is a Las Vegas algorithm that, given a prime p, a primitive root g and h ∈ F_p^*, computes log_g h mod
     (p − 1) in expected time L_p[1/2, c_DL(θ) + o(1)], uniformly in p.
   - The constant is c_DL(θ) = min_{a>0} max(a + 1/√2, a + θ/(2a), 3a). This equals Theorem A's (1 + θ)/√2 whenever
     θ ≤ 1/2, which covers every tuple with c < 1.
   - Exact values are zero-slack infima (η₀ ↓ 0). With fixed slack η₀ > 0 add η₀/√2 in the (1 + θ)/√2 branch.
     - d = 9, k = 13: 55√2/78 + η₀/√2 = 0.99720 + η₀/√2, below 1 for η₀ < 0.0039571.
     - Infimum at d = 9: 13√2/20 = 0.91924.
     - d = 5: 3√2/4 = 1.0607.
     - d = 3: 3√3/4 = 1.2990, where the linear algebra binds.
     - Infimum over all fixed tuples: 1/√2, not attained.
   - The same bound holds for an arbitrary base g with the promise h ∈ ⟨g⟩, via auxiliary random bases and restarts
     (Corollary DL′, following codex dlog_1 §4).
   - No GRH, no factorisation of p − 1, no class-group or unit-group computation, no Schirokauer maps and no
     power-residue characters are used.
2. **The DL-specific obstruction is removed by exact dependencies (Proposition XK, proved).**
   - The linear algebra is done over Z, not modulo the prime divisors ℓ of p − 1. Rational and prime-ideal exponents
     must cancel exactly, and units are removed by archimedean logarithms, LLL and the Blanksby–Montgomery gap.
   - Then ∏T_i^{x_i} is ±1 exactly. One primitive dependency gives A + X·B ≡ λ (mod p − 1), with B uniform mod p − 1
     because of the hidden exponents (Lemma HX).
   - All prime powers of p − 1 are handled at once, and the cost is L^{3a+o(1)}.
   - This quantified cost (post-processing exponent C_LA = 3) is what keeps the d = 9 constants. Codex dlog_1 shows that a post-processing exponent
     C_LA preserves them iff C_LA ≤ 1 + 1/θ (55/16 at k = 13, 13/3 for the infimum). The black-box route through
     Ge / Lenstra–Silverberg Thm 1.11 has an unquantified exponent.
3. **Literature (known).**
   - The best rigorous prime-field DL bound we could find remains Pomerance 1987: L_p[1/2, √2] precomputation plus
     L_p[1/2, 1/√2] per logarithm.
   - We found nothing below √2 for all primes or for almost all primes, unconditionally or under GRH.
   - Theorem DL would therefore be the first, but only as a reduction to the internal chain. No novelty claim beyond
     "no prior formulation found by us".
4. **What does not transfer.**
   - lv6's square-class machinery does not extend uniformly to ℓ-th powers for large ℓ | p − 1 (§6). The Kummer degree
     grows with ℓ, character values need an order-ℓ DL to be linearised, and Schirokauer maps need a Leopoldt/Wieferich
     condition that is not known to hold for a prescribed ℓ.
   - Within the cost model analysed (§7), the DL setting does not lower the shape-test term. The constant is the
     factoring constant, not better. This is not a lower bound on algorithms.
5. **Independent convergence.** Two blind GPT-6.1-sol rounds (r19 §2, dlog_1 §4) reached the exact-dependency route
   via Ge's algorithm, independently of my derivation. Neither quantified the post-processing exponent; dlog_1 named that
   as the missing "fast exact-dependence lemma". Proposition XK supplies it with C_LA = 3.
6. **Recursion.**
   - The shape test is integer factorisation of |u₀|, |v₀| ≤ √(2p). Any factoring constant τ (in L_N units) that agent
     `recurse` proves for that subroutine therefore transfers verbatim: c_DL = min_a max(a + τ/√2, a + θ/(2a), 3a).
   - In DL the "halves" cannot be handled by recursive *logarithms*. Their logs live in the same group F_p^*, so there
     is no size reduction.
   - Not developed here (owned by `recurse`).
7. **Weakest links.**
   - (i) The internal supply chain, especially Lemma 8.5′ Part 1 (agent `lem85`).
   - (ii) The transfer of Steps 1–7 to prime N = p (§3). It is mechanical, and both codex rounds checked it.
   - (iii) The precision bookkeeping of the unit-elimination step (Lemma XK4).
8. **Toy runs (experimentally supported, correctness only).**
   - 36/36 pre-registered instances (30–50-bit p) were solved correctly, with 0 row, kernel, torsion or answer
     violations.
   - 3465 certified-torsion dependencies all had W = ±1 and passed 3 auxiliary-prime checks. 702 of them had a
     nontrivial algebraic part.
   - 0 of 707 non-torsion kernel vectors with a nontrivial algebraic part passed. This shows the unit step is necessary
     and effective.
   - The hidden-exponent law held for ℓ = 2, 3, 5, 7 (all p ≥ 0.29).
   - One prediction failed. P6b ("full unit rank seen in ≥ 9/12 E2 runs") came out 5/12, because the pair box often
     produced no nontrivial γ. This is a design shortfall, not a violation (§9).

---

## 1. Literature

| problem | best rigorous bound found | hypotheses | reference |
|---|---|---|---|
| DL in F_p^* (p prime), all p | L_p[1/2, √2] precomputation, L_p[1/2, 1/√2] per logarithm (expected) | none (ECM-based smoothness test, Wiedemann) | Pomerance, *Fast, rigorous factorization and discrete logarithm algorithms*, Discrete Algorithms and Complexity (1987) 119–143, §1 (read from the scan: "L(q)^{√2+o(1)} for a preprocessing stage and … L(q)^{√1/2+o(1)} for the actual discrete logarithm", q prime or 2^n) |
| same, earlier | L[1/2, 2√2] (Adleman's sketch), √(5/2) with Pomerance's tools | — | Pomerance 1987 §1 |
| same, heuristic | L_p[1/3, (64/9)^{1/3}] (NFS; Gordon 1993: 3^{2/3}) | heuristic, not GRH-conditional | Gordon, SIAM J. Discrete Math. 6 (1993); Schirokauer 2008 survey §2.8 |
| F_{p^n}, p ≤ n^{o(n)} | L_q[1/2, √2] | none | Lovorn Bender–Pomerance, AMS/IP Stud. Adv. Math. 7 (1998) 221–232, §1. Also L(q)^{√(8/3)} for p < n^{4n/3}, and p^{2+o(1)} above that |
| F_{p^n}, n = 2 | rigorous subexponential (Lovorn's thesis; constant 3/2 commonly quoted, not verified) | none | Lovorn 1992 thesis, cited in Bender–Pomerance §1 |
| F_{p^n}, any p, n | (pn)^{2 log₂ n + O(1)}; quasi-polynomial for fixed p; L_Q[α + o(1)] for p = L_Q[α] (a fixed constant needs log p ≍ t^α(log t)^{−α}, t = log Q; codex dlog_1 §1(b)) | none | Kleinjung–Wesolowski, JAMS 35 (2022) 581–624, Thm 1.1 (arXiv 1906.10668) |
| imag. quadratic class group structure | L[1/2, √2] | ERH | Hafner–McCurley, JAMS 2 (1989) |
| imag. quadratic DL / class group | L[1/2, 3√2/4] | GRH | Vollmer, ANTS IV (2000), LNCS 1838; Biasse–Erukulangara, AMC 17 (2023) for structure |
| imag. quadratic DL, cyclic class group of known order | L[1/2, 1] | ERH | Enge–Gaudry, Acta Arith. 102 (2002), Thm 1 (per codex dlog_1 §1(c)) |
| arbitrary number fields: class/unit group | L_{|Δ|}[1/2, O(1)] (fixed degree) | ERH | de Boer–Pellet-Mary–Wesolowski, arXiv 2512.01588 (2025), Thm 1.3; finite-field DL explicitly not investigated |
| hyperelliptic Jacobians, g ≥ ϑ log q | L_{q^g}[1/2, (5/√6)(√(1+3/(2ϑ)) + √(3/(2ϑ)))] | none | Enge, Math. Comp. 71 (2002), Thm 1.1 |
| same, cyclic of known order, g/log q → ∞ | L_{q^g}[1/2, √2] | none | Enge–Gaudry, Acta Arith. 102 (2002), Thm 1 |
| fixed-genus curves | Õ(q^{2−2/g}) (exponential in log q^g) | none | Diem, Math. Comp. 80 (2011), Thm 1 |
| exact multiplicative relations in a Q-algebra | deterministic polynomial time (exponent not stated) | none | Ge 1993 thesis (fields); Lenstra–Silverberg, arXiv 1509.08843, Thm 1.11(i) |

- **known (search result).** No rigorous prime-field DL with constant below √2, for all p or for almost all p, under
  any standard hypothesis. Searched by me (WebSearch, the Galbraith book ch. 15, the KW introduction, the
  Bender–Pomerance introduction) and by codex dlog_1 §1 and r19 §2.
- This is a search result, not a theorem. Factoring has LP 1992's L[1/2, 1], but its class-group method has no
  known DL analogue.

## 2. Setting and the algorithm

- L = L_p[1/2, 1].
- Parameters (d, k, μ, η₀), θ, B = L^a, S, J, V_ℓ, M, H, W are exactly those of `notes/thmA_proof.md` §1, with N := p.
  The only change is that a is a free fixed constant in (0, √(θ/2)], not necessarily θ/√2. Steps 1–7 are proved for
  any such fixed a (proof critic, issue 2). Three facts are used:
  - T_SQ = L^{max(a, θ/(2a))+o(1)}, with S = L^{θ/(4a)+o(1)}. The resulting extra row term 2a is dominated by 3a;
  - L2_fam's free constant C₀ can be chosen with aC₀ > θ/a + 1;
  - the Ψ-estimates are uniform for a in a compact set.
  The choice of a matters only for θ > 1/2 (the d = 3 row), since a = θ/√2 otherwise.
- n := p − 1.
- The input is p prime ≥ p₀(d, k, μ, η₀), a primitive root g and h ∈ F_p^*. Write X := log_g h.

```
DLOG(p, g, h):
  A := field selection exactly as thmA SPLIT Phase F, with the target law below  # SAMPLE_FIELD incl. rejections
  loop (rounds):
    rows := []
    while |rows| < C̄(A) + d:          # computable bound C̄ := (d + 1)π(B) + 1 ≥ C_col
      draw a, b uniform in Z/n;  T := g^a h^b mod p        # T uniform on F_p^*
      row := DESCEND_ONCE(A, p) with R := T                 # Gauss split T = u0/v0, LP shape test, special-q
      if row ≠ FAIL: store (row, a, b)
    (x*, λ) := EXACT_DEPENDENCE(A, rows)                    # Proposition XK; x* primitive in Z^M
    A_ := Σ x*_i a_i,  B_ := Σ x*_i b_i  (mod n)
    if gcd(B_, n) = 1: X := (λ − A_)·B_^{−1} mod n; if g^X = h return X
```

Row algebra (thmA §3 with N := p). Notation for one row:
- T ≡ u₀/v₀ (mod p), with u₀v₀ ≠ 0 and |u₀|, |v₀| ≤ √(2p).
- Every large prime ℓ | u₀v₀ has a special-q atom P_ℓ = (a_ℓ, b_ℓ) with G(P_ℓ) = σ_ℓ ℓ n_ℓ, n_ℓ B-smooth.
- F(P_ℓ) ≠ 0 is B-smooth, and gcd(a_ℓ, b_ℓ) = 1.
- Let m := g₀g₁^{−1} mod p and ψ(a − bα) := a − bm. Then G(P) = g₁(a − bm) − b(g₀ − g₁m) ≡ g₁ψ(a − bα) (mod p), since g₀ ≡ g₁m.

Hence
  T ≡ ε · ∏_{q≤B} q^{e_q} · g₁^{E} · ∏_ℓ ψ(a_ℓ − b_ℓα)^{e_ℓ}  (mod p),        (2.1)
where:
- ε ∈ {±1};
- e_ℓ = +1 for ℓ | u₀ and −1 for ℓ | v₀;
- E := Σ_ℓ e_ℓ;
- e_q = v_q(u₀) − v_q(v₀) − Σ_ℓ e_ℓ v_q(n_ℓ).

With the relaxed toy shape, the B-smooth parts of u₀, v₀ also enter e_q. The theorem's shape has none, which changes
nothing.

## 3. Transfer of Theorem A, Steps 0–7, to N = p (proved from A)

Each step of `notes/thmA_proof.md` §4 is listed with how it transfers. N's arithmetic enters the sources only through:
(N, t) = 1 for t ≤ 2MV; N | F(g₀, g₁); trial division / 2^ω; and the law of R (checked by grep of
pcest/ldbound/lv6/hypD and by reading the statements).

| step | content | N = p |
|---|---|---|
| 0 | preprocessing, 2^ω = L^{o(1)} | Not needed. ω = 1, p > B. |
| 1 | Lemma GB_w, sizes, Lemma 1.2_w, coprimality (1.1) | Verbatim. GB applies to any integer Nw. (t, p) = 1 since 0 < t ≤ 2MV < p. |
| 2 | SQ_SEARCH recognises N_ℓ(A) | Verbatim. The gcd checks never exit, since all checked quantities are nonzero and < p in absolute value: \|G\| ≤ 2MV, \|F_A\| ≤ X* = p^{μ+h+dν+o(1)} < p, and 0 < \|f_d\| < p. |
| 3 | FM (Lemma 8.5′) | Verbatim. Uses only (N, t) = 1 and the field law. |
| 4 | LD* (Kern, Lemma W, L2_fam, DS) | Verbatim. Same reason. |
| 5 | PZ ⇒ D_k(γ), γ ≥ L^{−o(1)} | Verbatim. |
| 6 | rejections, factor exits | Rejections are verbatim, keeping IR, e = w, and p \| f′(m). There are no factor exits, so the absorbing bookkeeping is trivial (δ = 0). |
| 7 | Theorems S′, K, field selection | R uniform on Q_N is replaced by T uniform on F_p^*: Lemma S_p below. The 2^{−ω} factors become 1. |
| 8 | LV6 (characters, Hensel root) | **Replaced** by Proposition XK (§4). |
| 9 | DD + (T) | **Replaced** by Lemma HX (§4.6). |
| 10 | cost | §5. The LV6 term L^{2a} becomes L^{3a}. |

**Lemma S_p (proved; Theorem S′ with h = 1).** Let G be a set of primes in a dyadic interval with |G| ≥ 4k². Let S
be the products of exactly k distinct primes of G in a window [Y/2, Y] with 2Y² < p and |S| ≥ C(|G|, k)/(2k). If T is
uniform on F_p^*, then
  Pr[∃ u, v ∈ S : u ≡ Tv (mod p)] ≥ |S|²/(c_k(p − 1)),   c_k = 8k²C(2k, k).

*Proof.* Put n(T) := #{(u, v) ∈ S² : u ≡ Tv}.
- Σ_T n(T) = |S|² exactly. Each pair determines T = uv^{−1}, and no square-class condition is imposed.
- The second-moment bound Σ n(T)² ≤ c_k|S|² is hypD Theorem S′'s count, steps 1–4, verbatim. It uses only
  |uv′ − u′v| < Y² < p ⇒ uv′ = u′v and the shared-prime multiplicity count.
- Cauchy–Schwarz gives #{T : n(T) > 0} ≥ |S|⁴/(c_k|S|²). Divide by p − 1. ∎

Consequently Theorem K(ii) holds with the 2^{−ω} removed. K(iii) (the Gauss-reduced primitive vector is the split,
since 2Y² < p) and the divisor-closed acceptance rule hold verbatim. The kept field has E[1/q_A | kept] ≤ 1/μ_sel with
μ_sel ≥ L^{−o(1)} (thmA Step 7).

**Lemma T_p (target law; proved).**
- Since g generates F_p^*, (a, b) ↦ (log_g T, b) = (a + Xb, b) is a bijection of (Z/n)². So T is uniform, and given
  T the label b is uniform.
- Every accepted row, every field decision and every later computation except the final A_, B_ is a function of the
  T's and of coins independent of the labels (a_i, b_i). Within a round j, the rows' labels b_i are therefore i.i.d.
  uniform on Z/n conditional on 𝒯_j: all T's, all coins, and all labels and extraction outcomes of earlier rounds. 𝒯_j excludes round j's labels and extraction outcome.
  - The decision to start round j + 1 depends on round j's labels, through gcd(B_, n). So the labels are not uniform
    given the transcript of all rounds, and HX is stated per round (proof critic, issue 7).
  - Each round draws fresh labels, and b is independent of T, which is stronger than conditional uniformity. Hence
    the rounds succeed independently, each with probability ≥ φ(n)/n.
- This is thmA's (T), with r ↦ εr replaced by b ↦ b + c, a ↦ a − Xc. ∎

## 4. Proposition XK: exact dependencies (proved)

Fix a kept field A (f irreducible, content 1, 0 < |f_d| < p, m a simple root of f mod p) and rows i = 1..M with
labels (a_i, b_i). Notation:
- K = Q(α), f(α) = 0, t := f_d, β := tα.
- 𝔞 := O_K + αO_K, a fractional ideal.
- For an atom P = (a, b) with gcd(a, b) = 1: 𝔍_P := (a − bα)𝔞^{−1}.

**Lemma XK0 (ideal bookkeeping; standard).**
- (i) 𝔍_P is integral, and N(𝔍_P) = |F(a, b)|.
- (ii) If r ∤ t·disc(f), the prime ideals of O_K dividing 𝔍_P above r are first-degree, (r, s) with s ≡ ab^{−1}
  (mod r), and v_{(r,s)}(𝔍_P) = v_r(F(a, b)).
- (iii) For each of the O(d log p) primes r ≤ B with r | t·disc(f), an r-maximal order, the primes 𝔯 | r and the
  valuations v_𝔯(𝔍_P) are computable in time poly(d, log r, log p).

*Proof.*
- (i) Integrality: a − bα ∈ 𝔞. For the norm: f is primitive, so Gauss's content lemma for Dedekind domains (over the
  Galois closure, content(f) = (t)·∏_σ(1, σα)) gives N(𝔞) = |t|^{−1}. Then N(𝔍_P) = |N_{K/Q}(a − bα)|·|t| =
  |F(a, b)|.
- (ii) For r ∤ t, α is r-integral and 𝔞 is trivial at r. For r ∤ t·disc f, Z_(r)[α] is the r-maximal order and Dedekind–Kummer applies. If r | F(a, b)
  then r ∤ b (otherwise r | ta^d, so r | t), so a − bα lies only in (r, α − ab^{−1}), with valuation
  v_r(N(a − bα)) = v_r(F(a, b)).
- (iii) Buchmann–Lenstra (1994) / Cohen GTM 138 Alg. 6.1.8 (Round 2 at a known prime), 6.2.9 (decomposition),
  4.8.17 (valuations). Then v_𝔯(𝔍_P) = v_𝔯(z_P) − min(v_𝔯(t), v_𝔯(β)), with z_P = ta − bβ ∈ Z[β] and
  t𝔞 = (t, β). Finding the bad r: compute t·disc(f) mod r for r ≤ B, in Õ(B) time. ∎

**Columns.** One column per rational prime q ≤ B (entry e_q), one column "g₁" (entry E_i), and one column per prime
ideal 𝔯 above r ≤ B that occurs (entry Σ_ℓ e_ℓ v_𝔯(𝔍_{P_ℓ})). Hence C_col ≤ (d + 1)π(B) + 1. Signs are not columns.
Let E ∈ Z^{M×C_col} be the exact exponent matrix and
  Λ₁ := {x ∈ Z^M : xᵀE = 0},  γ_x := ∏_P (a_P − b_Pα)^{y_P(x)},  y_P(x) := Σ_i x_i Σ_{ℓ: P_ℓ = P} e_{iℓ}.
(Atoms are cached per (A, ℓ), so P is a function of ℓ.)

**Lemma XK1 (exactness ⇒ units).** For x ∈ Λ₁, γ_x ∈ O_K^×.
- *Proof.* (γ_x) = ∏_P 𝔍_P^{y_P}·𝔞^{Σ_P y_P}.
  - Σ_P y_P = Σ_i x_i E_i = 0 (the g₁ column).
  - For every 𝔯 above r ≤ B, v_𝔯 = Σ_i x_i E_{i,𝔯} = 0.
  - Every prime dividing some 𝔍_P lies above some r ≤ B, because N(𝔍_P) = |F(P)| is B-smooth, and every such 𝔯 is a
    column.
  - So (γ_x) = O_K. ∎

**Lemma XK2 (the reduction map).** Let 𝔭 be the prime of O_K corresponding to the root m: the unique prime above p
with α ≡ m (mod 𝔭). Other simple roots of f mod p give other degree-1 primes, which are not used. It exists, has
residue field F_p and is unique for this m, because p ∤ t and m is a simple root mod p (the rejection p | f′(m) of
thmA SAMPLE_FIELD). Reduction on the local unit group, red_𝔭 : O_{K,𝔭}^× → F_p^×, is a homomorphism with red_𝔭(a − bα) = a − bm = ψ(a − bα) ≠ 0. For x ∈ Λ₁:
  W(x) := ∏_i T_i^{x_i} ≡ ±red_𝔭(γ_x)  (mod p).                                          (4.1)
- *Proof.*
  - Each a − bα is a 𝔭-unit. 𝔍_P has norm |F(P)| < p, and v_𝔭(𝔞) = 0 since p ∤ t.
  - Take (2.1) to the power x_i and multiply. The q-, g₁- and sign-factors have exponents xᵀE = 0 in the rational
    columns (the signs give ±1).
  - The ψ-factors give ∏_P ψ(a_P − b_Pα)^{y_P} = red_𝔭(γ_x). ∎

**Lemma XK3 (lattices).**
- Λ₁ is saturated (it is the kernel of an integer matrix), and rank Λ₁ ≥ M − C_col.
- The map u: Λ₁ → O_K^×/μ_K, x ↦ [γ_x], is a homomorphism to a free group of rank r_U := r₁ + r₂ − 1 ≤ d − 1.
- Λ₂ := ker u = {x ∈ Λ₁ : γ_x ∈ μ_K} is saturated in Z^M (if jx ∈ Λ₂ then γ_x^j ∈ μ_K, so γ_x ∈ μ_K), and
  rank Λ₂ ≥ M − C_col − r_U ≥ 1 when M ≥ C_col + d.
- Every primitive vector of Λ₂ is primitive in Z^M. ∎

**Lemma XK4 (unit elimination; proved).** Given s := r_U + 1 ≤ d independent x₁, …, x_s ∈ Λ₁ with entries of at
most β_x bits, one finds a nonzero n ∈ Z^s with Σ_j n_j x_j ∈ Λ₂ in time Õ_d(M(β_x + log p) + (β_x + log p)²).

*Proof.*
- **Log vectors.** Let Log: K^× → R^{r₁+r₂}, Log(z) = (log|σ_i z|)_{real} ∪ (2 log|σ_j z|)_{complex pairs}, and
  v_j := Log(γ_{x_j}) = Σ_P y_P(x_j)·Log(a_P − b_Pα). These are logs of units, so they lie in the unit lattice
  Λ_U = Log(O_K^×), of rank r_U.
- **Gap (Blanksby–Montgomery 1971).** A non-torsion algebraic integer of degree d′ ≤ d has Mahler measure
  > 1 + 1/(52d′ log 6d′). For a unit z of K, log M(z) ≤ d·‖Log z‖_∞. So every nonzero w ∈ Λ_U has
  ‖w‖_∞ ≥ β_d := log(1 + 1/(52d log 6d))/d > 0. By Kronecker, Log z = 0 iff z ∈ μ_K.
- **Relations exist and are small.**
  - Rel := {n ∈ Z^s : Σ n_j v_j = 0} has rank ≥ s − r_U ≥ 1.
  - Fix an LLL-reduced basis u_1..u_{r_U} of Λ_U (only for the proof). Its Gram–Schmidt norms are ≥ 2^{−r_U/2}β_d.
  - So the integer coordinates c_{ji} of v_j satisfy |c_{ji}| ≤ 2^{O(d²)}‖v_j‖/β_d.
  - Cramer on an r_U × (r_U + 1) integer matrix gives a nonzero n₀ ∈ Rel with
    ‖n₀‖ ≤ R := max(1, d!·(2^{O(d²)} max_j‖v_j‖/β_d)^{d}).
    If the coordinate matrix has rank ρ_c < r_U, apply the minor argument at rank ρ_c. If all v_j = 0, any e_j is a
    relation and R = 1 bounds it. That is the E1 case (codex dlog_2 #3).
  - ‖v_j‖ ≤ 2^{β_x}·M·2k·max_P‖Log(a_P − b_Pα)‖, and the last factor is O_d(log p) from both sides:
    - Upper: |σ_i(a − bα)| ≤ U := |a| + |b|·max|σα| = p^{O(1)}.
    - Lower: |σ_i(a − bα)| = |N(a − bα)|/∏_{j≠i}|σ_j(a − bα)| ≥ (|F(a,b)|/|t|)/U^{d−1} ≥ p^{−O(d)}.
    So log R = O_d(β_x + log p). The same lower bound shows that absolute precision κ + β_x + O(d log p) bits on the
    roots of f gives the relative precision needed for the logs.
- **LLL.** Let b_j := (e_j, ⌊2^κ v_j⌉) ∈ Z^{s+r₁+r₂}, where ⌊2^κ v_j⌉ is any integer vector within 1 of 2^κ v_j in
  each coordinate (computed from approximations).
  - If n ∈ Rel, then ‖Σ n_j b_j‖ ≤ (1 + s√d)‖n‖.
  - If n ∉ Rel, then ‖Σ n_j b_j‖ ≥ 2^κβ_d − s√d‖n‖.
  - LLL's first vector b has ‖b‖ ≤ 2^{(s−1)/2}(1 + s√d)R =: R′, and its n-part has ‖n‖ ≤ R′.
  - Choose 2^κβ_d > (1 + s√d)R′. Then b's n-part lies in Rel, since otherwise ‖b‖ > R′. And n ≠ 0, since n = 0 would
    make b = 0.
  - So κ = O_d(β_x + log p).
- **Precision.**
  - Roots of f to κ + O(β_x + log p) bits: Õ(d²κ) (Schönhage / Pan).
  - Each log|a − bσα| to that precision: Õ(κ).
  - The sums Σ_P y_P·Log(·) over ≤ 2kM atoms: Õ(M(κ + β_x)).
  - LLL in fixed dimension s + r₁ + r₂ ≤ 2d with O(κ)-bit entries: Õ(κ²), using a quadratic-in-bit-length LLL such as
    Nguyen–Stehlé L² (SIAM J. Comput. 39, 2009).
  - Precision ledger (codex dlog_2 #3): absolute precision P ≥ κ + β_x + log₂M + O_{d,k}(1) on each atom logarithm
    suffices for the rounding, the accumulation and the certification.
- **Certification.** In the algorithm, accept n only if the ball |Σ n_j v_j|_∞ < β_d/2 is certified. A true relation
  passes, by the precision. A non-relation cannot pass, by the gap. ∎

**Lemma XK5 (torsion ⇒ known logarithm).** For x ∈ Λ₂, W(x) ∈ μ_{2w}(F_p), where w := |μ_K|. For d odd, μ_K = {±1},
so W(x) = ±1. Then λ := log_g W(x) ∈ {0, n/2} is read off. In general λ = j·n/e with e = gcd(2w, n), and j < e is
found by trying.
- *Proof.* By (4.1), W(x) = ±red_𝔭(ζ) with ζ ∈ μ_K, and red_𝔭 maps μ_K into μ_w(F_p). ∎

**Lemma HX (hidden exponents; proved).** Let x* ∈ Z^M be primitive, i.e. gcd of its entries is 1, and a function of
the round's transcript 𝒯_j (Lemma T_p). Then B_ := Σ x*_i b_i is uniform on Z/n given 𝒯_j, so Pr[gcd(B_, n) = 1 | 𝒯_j] = φ(n)/n ≥
1/(1 + log₂ n).
- If gcd(B_, n) = 1 and A_ := Σ x*_i a_i, then A_ + X·B_ ≡ λ (mod n), so X = (λ − A_)B_^{−1} mod n. This is exact
  for every prime power of n at once, with no factorisation of n.
- *Proof.*
  - By Lemma T_p the b_i are i.i.d. uniform given 𝒯.
  - b ↦ b·x* mod n is a surjective homomorphism, since gcd(x*) = 1 and so a Bezout vector exists. Hence the image of
    the uniform law is uniform.
  - log_g W(x*) = Σ x*_i(a_i + Xb_i) = A_ + X·B_ (mod n), by T_i = g^{a_i}h^{b_i}. ∎
  - The bound φ(n)/n ≥ 1/(1 + log₂ n) is codex dlog_1's elementary telescoping. Mertens gives ≫ 1/log log n.

**Algorithm EXACT_DEPENDENCE and Proposition XK (proved).** On M ≥ C̄ + d ≥ C_col + d rows:
1. Compute E exactly (Lemma XK0).
2. Pick a random prime q ∈ [Q, 2Q] with Q := M·(log p)³. Compute by Gaussian elimination mod q a rank profile
   (I, J), |I| = |J| = ρ.
3. For d rows i ∉ I, solve yᵀE_{I,J} = −E_{i,J} over Q by Dixon's p-adic lifting (Dixon 1982), and clear
   denominators: x_i := D·(y, e_i).
4. Verify x_iᵀE = 0 exactly on all columns. If it fails, then rank_Q E > ρ: restart step 2.
5. Apply Lemma XK4 to s = r_U + 1 of the x_i, giving n ≠ 0.
6. Put x* := Σ n_j x_j / gcd. Then x* ∈ Λ₂ is primitive (Lemma XK3).
7. Compute W(x*), λ (Lemma XK5) and A_, B_ (Lemma HX).

**Proposition XK.** Steps 1–7 are Las Vegas. They output a primitive x* ∈ Λ₂ and λ with log_g W(x*) = λ in expected
time L^{3a+o(1)} when M = BL^{o(1)} and B = L^a.

*Proof.*
- **Correctness.** Lemmas XK1–XK5.
- **Step 2.**
  - The ρ_Q × ρ_Q minor of a fixed basis has |det| ≤ ∏‖rows‖ ≤ (O(k log p))^{1.5ρ}. Each row has ≤ 2k(1 + log₂p)
    nonzeros, each O(k log p).
  - So at most O(ρ log log p) primes divide it. [Q, 2Q] contains ≫ M(log p)²/log M ≫ ρ(log p) primes, so
    Pr[rank mod q < rank_Q] ≤ 1/2.
  - A nonsingular minor mod q is nonsingular over Q, so ρ ≤ rank_Q.
  - Step 4 does *not* certify ρ = rank_Q (codex dlog_2 #4). The chosen d rows may already lie in the rational span of
    the rows in I.
  - That does no harm. If verification passes, the d vectors are valid kernel vectors, independent through their free
    coordinates e_i, which is all XK4 needs. If it fails, resample q.
  - A q that preserves a fixed maximal minor (Pr ≥ 1/2) guarantees that verification passes, so the expected number of
    restarts is ≤ 2.
  - Cost O(M·C_col·ρ) operations mod q = L^{3a+o(1)}.
- **Step 3.**
  - Dixon: one Q-inverse mod a prime (ρ³), then O(ρ log(ρ‖E‖)) lifting steps of O(ρ²) each, for each of d right-hand
    sides. Total Õ(ρ³) = L^{3a+o(1)}.
  - The outputs have β_x = O(ρ log log p) = L^{a+o(1)} bits (Cramer/Hadamard).
- **Steps 4–7.**
  - Step 4: nnz(E)·β_x = L^{2a+o(1)}.
  - Step 5: Lemma XK4 with β_x = L^{a+o(1)} costs L^{2a+o(1)}.
  - Step 7: M modular exponentiations with L^{a+o(1)}-bit exponents cost L^{2a+o(1)}.
  - Lemma XK0(iii): Õ(B) + poly(log p).
- ∎

Remarks.
- (i) **C_LA = 3** comes from dense elimination. A sparse rank profile plus sparse rational solving would give smaller
  C_LA, but rigorous sparse integer solvers (EGGSV 2006) rely on an unproved block-projection conjecture, so we do not
  use them. Only C_LA ≤ 1 + 1/θ matters for the constants.
- (ii) **Black-box alternative.** Ge / Lenstra–Silverberg Thm 1.11(i) computes generators of the multiplicative
  relation lattice of the explicit row elements in polynomial time with unspecified exponent. This is the route of r19
  §2 and dlog_1 §4. It gives Theorem DL with 3a replaced by Ca, and C unquantified.
- (iii) **Factoring.** The same device gives a character-free replacement for LV6 in Theorem A (also noted in r19 §2).
  - Take Λ₂′ := {x ∈ Λ₂ : γ_x ∈ μ_K², sign product +1}, a sublattice of index ≤ 4.
  - Its reduction mod 2 has dimension ≥ rank Λ₂ − 2, so with M ≥ C_col + d + 2 there is x with some odd x_i.
  - Then (∏r_i^{x_i})² ≡ ψ(ζ′)² with ζ′ ∈ μ_K computable, and DD applies.
  - For composite N, ψ must be well defined on the group generated by the a − bα. It is, because Z[α] ≅ Z[x]/(f) for
    primitive f (Gauss's lemma), and one localises at the elements whose ψ-images are units, all of which are
    gcd-checked.
  - Cost L^{3a} instead of L^{2a}. That is harmless for θ ≤ 1/2, i.e. whenever c < 1.
  - It would remove Lagarias–Odlyzko, Stark/Landau and the Hensel root from Theorem A.
  - Not worked out further here (lv6's area).
  - **Algorithmic sketch: existence proved, cost not established.** XK outputs one primitive torsion relation. The
    factoring variant needs enough independent torsion relations, the saturation of their span, and the two sign
    constraints imposed within budget (codex dlog_2 #12).

## 5. Theorem DL, cost and constants (proved from A)

**Theorem DL.** Fix admissible (d, k, μ, η₀) with θ < 1. DLOG is Las Vegas and has expected time
  L_p[1/2, c_DL(θ) + o(1)],  c_DL(θ) := min_{a>0} max(a + 1/√2, a + θ/(2a), 3a),
uniformly for primes p ≥ p₀(d, k, μ, η₀). For θ ≤ 1/2, c_DL = (1 + θ)/√2, attained at a = θ/√2.

*Proof.*
- **Rows.** With B = L^a and the same per-row cost as Theorem A's Step 7 (via §3 and Lemma S_p), C_col + d =
  BL^{o(1)} rows cost L^{a + max(1/√2, θ/(2a)) + o(1)}. Field selection costs L^{max(1/√2, θ/(2a)) + o(1)}.
- **Post-processing** costs L^{3a+o(1)} (Proposition XK).
- **Rounds** are geometric with success ≥ φ(n)/n (Lemma HX), so the expected number is ≤ 1 + log₂ n = L^{o(1)}.
  Fresh rows are drawn each round, so the rounds are independent. The field is kept.
- **Total.** Wald's identity gives L^{max(a + 1/√2, a + θ/(2a), 3a) + o(1)}.
- **Minimisation.** Let f₁ = a + 1/√2, f₂ = a + θ/(2a), f₃ = 3a.
  - At a = θ/√2: f₁ = f₂ = (1 + θ)/√2, and f₃ = 3θ/√2 ≤ (1 + θ)/√2 iff θ ≤ 1/2.
  - min_a max(f₁, f₂) = (1 + θ)/√2 for θ < 1 (thmA Step 10), so no a does better.
  - For θ > 1/2, the minimum is where f₂ = f₃, at a = √θ/2, giving c = 3√θ/2. Check: f₁ = √θ/2 + 1/√2 ≤ 3√θ/2 iff
    θ ≥ 1/2.
  - `src/dlog/params.py` computes all of this exactly with sympy. ∎

So c_DL(θ) = (1 + θ)/√2 for θ ≤ 1/2, and 3√θ/2 for 1/2 < θ < 1. It is < √2 iff θ < 8/9, since 3√θ/2 < √2
iff θ < 8/9. This includes every optimised tuple: the minimal θ at d = 3 is 3/4 + O(1/k). Admissible but unoptimised
tuples with 8/9 ≤ θ < 1 do not beat √2.

| d | k | θ | c_DL (exact) | ≈ | Theorem A | LA binds? |
|---|---|---|---|---|---|---|
| 3 | ∞ (inf) | 3/4 | 3√3/4 | 1.29904 | 1.23744 | yes |
| 5 | ∞ (inf) | 1/2 | 3√2/4 | 1.06066 | 1.06066 | boundary |
| 7 | 21 (min for c < 1) | 139/336 | 475√2/672 | 0.99963 | same | no |
| 9 | **13** (min) | **16/39** | **55√2/78** | **0.99720** | same | no |
| 9 | 50 | 317/1000 | 1317√2/2000 | 0.93126 | same | no |
| 9 | ∞ (inf) | 3/10 | 13√2/20 | 0.91924 | same | no |
| 11 | ∞ (inf) | 1/4 | 5√2/8 | 0.88388 | same | no |
| 19 | ∞ (inf) | 3/20 | 23√2/40 | 0.81317 | same | no |
| d → ∞ | | → 0 | → 1/√2 | 0.70711 | same | no |

- All rows are zero-slack values (η₀ ↓ 0); with slack add η₀/√2 where (1 + θ)/√2 applies. Rows "∞" are also infima over k, not attained. At d = 9, k = 13 one needs η₀ < √2 − 1 − 16/39 ≈ 0.003957
  for c < 1.
- Every tuple with Theorem-A constant < 1 has θ < √2 − 1 < 1/2, so the cubic linear algebra never binds there.
  This is tested in `tests/test_dlog.py`.

**Corollary DL′ (arbitrary base; proved from A, following codex dlog_1 §4).** For g, h ∈ F_p^* with h ∈ ⟨g⟩, log_g h
mod ord(g) can be computed in the same expected time.
- *Proof.* Draw a random auxiliary base t.
  - A primitive root occurs with probability φ(n)/n ≥ 1/(1 + log₂ n).
  - Run DLOG(p, t, g) and DLOG(p, t, h) under a time budget that is doubled each stage, with c·(1 + log₂ n) fresh t's
    per stage. The fixed c is chosen so that a stage whose budget exceeds 4× the expected primitive-root time fails
    with probability ≤ 1/4.
  - For primitive t and budget ≥ 4·E[time], Markov's inequality gives each of the two runs failure probability
    ≤ 1/4, so both succeed with probability ≥ 1/2.
  - Since stage failure is ≤ 1/4 < 1/2 once the budget is large enough, the doubling series converges.
  - Every output is verified (t^u = g, t^v = h). Then solve u·x ≡ v (mod n).
    - If t is primitive, this is consistent because h ∈ ⟨g⟩.
    - If t is non-primitive and both runs still return verified exponents, it may be inconsistent. If gcd(u, n) ∤ v,
      discard and continue.
    - Any solution x is checked against g^x = h before output.
  - For a non-primitive t, DLOG may fail to terminate within budget, but it never returns an unverified answer.
- Expected total: L^{c_DL + o(1)}. No factorisation of p − 1 is needed.
- Alternatively, factor p − 1 with Theorem A at L^{(1+θ)/√2} ≤ L^{c_DL}. That route uses Theorem A's Steps 8–10
  (LV6, DD) as well, so it depends on more of the internal chain than Theorem DL does. ∎

**Uniformity.**
- p₀ and the o(1) depend only on (d, k, μ, η₀), as in thmA §6 U3–U4.
- Primes p < p₀ take O(1) time.
- Nothing depends on the factorisation of p − 1.

## 6. The mod-ℓ route: what does and does not transfer

This answers the task's questions on Z/ℓ linear algebra, Schirokauer maps and ℓ-th-power classes. Let ℓ | p − 1 and
H_ℓ := {γ ∈ K^× : (γ) = 𝔟^ℓ (away from S)}/K^{×ℓ}. There is an exact sequence
0 → O_S^×/ℓ → H_ℓ → Cl(O_S)[ℓ] → 0 (codex dlog_1 §3 (7)). Its dimension is
≤ |S| + r_U + 1 + log h_K/log ℓ.

- **(a) Characters (lv6 transfer). The lv6 method's cost bounds degrade for large ℓ (a cost obstruction for this method, not an impossibility result); partial transfer for small ℓ.**
  - The ℓ-th power residue characters at degree-1 primes q ≡ 1 (mod ℓ) need effective Chebotarev in K(μ_ℓ, γ^{1/ℓ}).
    That extension has degree ≤ dℓ(ℓ − 1), so lv6 Lemma X's range log x ≫ n_L(log d_L)² becomes ≫ ℓ⁶ log²p.
  - Linearising μ_ℓ-valued symbols into Z/ℓ needs an order-ℓ DL in F_q^* for every character value. BSGS gives an Õ(√ℓ) upper bound, but no uniformly cheap linearisation for large ℓ is supplied: each of ≈ BL^{o(1)} values needs a second DL problem in a group of order ≥ ℓ.
    - Hence the lv6-style transfer is affordable for ℓ ≤ L^{o(1)}, but such ℓ are already handled by
    Pohlig–Hellman–BSGS in Õ(e√ℓ).
  - Exceptional zeros need care. Odd-order Kummer characters are non-real, but real quadratic characters from the
    cyclotomic base K(μ_ℓ) still enter the prime sampling. That would need a separate argument (codex dlog_2 #10), and
    it is not supplied: **unverified**.
  - For ℓ = 2, lv6 transfers verbatim and gives X mod 2. For 2^e and large ℓ it does not.
- **(b) Schirokauer maps. Conditional.**
  - They are linear with values in (Z/ℓ)^d, so they cannot be injective on H_ℓ when dim H_ℓ > d.
  - The S-unit valuations must therefore be handled by the ideal columns, and the maps used only on ordinary units of
    the good-prime domain. Example: K = Q, S = {2, 3}, ℓ = 5. The Fermat-quotient map kills 18 = 2·3², which is not a
    fifth power (codex dlog_2 #9).
  - On that domain, sufficient conditions are ℓ ∤ h_K (no class-group obstruction) and full rank of the Schirokauer
    map on O_K^×/O_K^{×ℓ} (a Leopoldt/Wieferich-type condition), as in Schirokauer's own exposition (2008, §2.2).
  - Neither can be proved for a prescribed ℓ and a random field in the family. Even "infinitely many non-Wieferich
    primes" is open for 2 ∈ Q. So this route gives only a conditional theorem.
- **(c) Prime powers ℓ^e.**
  - The mod-ℓ route gives only F_p^×/F_p^{×ℓ} information. Pohlig–Hellman lifting needs ℓ^e-th power classes or
    roots.
  - Exact dependencies (Prop. XK) give all of Z/(p − 1) at once.
- **(d) The answer to "how does lv6's square-class kernel analysis transfer".**
  - Its structure transfers: the dimension bound for H_ℓ is the analogue of lv6 Lemma H.
  - Its algorithms do not, for large ℓ.
  - The exact route makes H_ℓ irrelevant. It never forms a mod-ℓ kernel, so the class group cannot enter: only
    (γ) = (1) is used.
  - The toy run shows the integral-versus-mod-2 rank discrepancy that a GF(2) kernel would have. It reflects Smith-invariant divisibility, not specifically class-group torsion. dim ker(E mod 2) − dim ker_Z(E)
    was 1 in 6/36 E1/E2 runs (P7).

## 7. DL-specific questions

- **Is the split counting easier? Yes, but only by L^{o(1)}.** Lemma S_p removes 2^{−ω} and the square-class
  characters: Σ n(T) = |S|² exactly.
- **Can choosing the target's randomisation lower the shape cost? Not within the cost model below.** This optimises an assigned cost expression for the described implementations. It is not a lower bound on algorithms (codex dlog_2 #11).
  - Hybrid shapes: u₀, v₀ = (B-smooth part of size p^{1/2−ς}) × (cofactor ≤ p^ς made of J-primes). These need
    rigorous smoothness extraction (Pomerance's ECM theorem), and then LP on the cofactor at L^{√ς}.
  - The per-row exponent is f(ς) = (1 − 2ς)/(2a) + √ς. f is concave on [0, 1/2], so its minimum is at an endpoint:
    f(1/2) = 1/√2 versus f(0) = 1/(2a).
  - f(1/2) < f(0) iff a < 1/√2, which holds at a = θ/√2.
  - So, within this family and cost model, the pure k-product shape of Theorem A is the best choice.
  - Targets T = g^a h^b·(chosen v) do not help: v·T is uniform and ≈ p, not ≈ √p.
- **Recursion.** The DL "halves" are integers whose logarithms are as hard as any log in F_p^*, so there is no
  recursion on logs. The only recursion is the factoring subroutine.
  - With a factoring constant τ for integers ≤ √(2p), c_DL(τ) = min_a max(a + τ/√2, a + θ/(2a), 3a). With quadratic
    post-processing it would be min_a max(a + τ/√2, a + θ/(2a), 2a).
  - Exactly: with s := τ/√2,
    - c_DL(s) = s + θ/(2s) (at a = θ/(2s)) if s ≥ √θ;
    - c_DL(s) = 3√θ/2 (at a = √θ/2, where the cubic term binds) if s < √θ.
    So if `recurse` proves τ arbitrarily small, c_DL → 3√θ/2, which → 0 as θ → 0 with d → ∞. This is r19 §2's
    "every positive constant", and it needs no quadratic post-processing.
  - Agent `recurse` owns the fixpoint and its validity. Nothing here depends on it.
- **Two randomized number fields (codex dlog_1 §2 (6)).** Not needed. The design keeps Theorem A's rational side.
  dlog_1 shows that two independent fields would also work: E|G_{A₊} ∩ G_{A₋}|/|J| ≥ γ².
- **Cheaper rows. Conditional, not pursued.**
  - Mixing cheap NFS pair relations (no hidden labels) with few descended targets would need a spanning lemma: the
    target rows' exponent vectors must lie in the saturated span of the pair rows. This is rkrank's RK-type
    obstruction.
  - At fixed d this could lower the row term. The present analysis still retains an L^{1/√2} shape-test term for the
    single target descent without recursion.

## 8. Which Theorem-A lemmas transfer verbatim

- **Verbatim:**
  - Lemma GB_w, sizes, Lemma 1.2_w, coprimality (1.1);
  - Lemma 1.3 / Prop BS (SQ_SEARCH);
  - Lemma 8.5′, Theorem FM;
  - Lemma Kern, Lemma W / Cor. W, Theorem L2_fam, Lemma DS, Theorem LD*;
  - Prop PZ / Thm 4.3;
  - Lemma IR, the content and f′(m) rejections;
  - Theorem K(i), (iii), and field selection (Cor. K).
- **Modified:**
  - Theorem S′ → Lemma S_p (h = 1);
  - K(ii) without 2^{−ω};
  - (T) → Lemma T_p.
- **Replaced:**
  - LV6 (Lemmas C, H, E, X, SQ, congruence) → Proposition XK;
  - Theorem DD → Lemma HX;
  - cost 2a → 3a.
- **Dropped:** preprocessing and Lemma Ω; factor-exit bookkeeping.

## 9. Toy experiments (pre-registered, `notes/dlog_prereg.md`)

Pipeline: d = 3, μ = 1/3, thmA toy settings. Deviations are as pre-registered: toy shape, bad-prime atoms discarded,
Miller–Rabin, flint nullspace.

| prediction | result |
|---|---|
| P1 correctness | E1 **24/24**, E2 **12/12** (and X equal to the secret, 36/36) |
| P2 row identity | **0** violations (target and pair rows) |
| P3 torsion | **0** violations in 3465 certified-torsion dependencies (W = ±1 and 3 auxiliary primes). 702 had γ ≠ 1 |
| P4 exactness / units | **0** kernel violations, **0** non-unit kernel vectors. unit_rank_implied ≤ r₁ + r₂ − 1 in 36/36 |
| P5 hidden-exponent law (E1) | ℓ = 2: 210/398 = 0.528 (p = 0.29); ℓ = 3: 84/269 = 0.312 (p = 0.48); ℓ = 5: 12/68 (p = 0.76); ℓ = 7: 13/120 (p = 0.36). **Holds** |
| P6a non-torsion γ pass auxiliary ±1 | **0/707** (≤ 5% predicted). **Holds** |
| P6b full unit rank in ≥ 9/12 E2 runs | **5/12. FAILS.** The 7 failing runs had no kernel vector with γ ≠ 1, since too few pair rows came from the 300-box. In all 5 runs with γ ≠ 1, unit rank = r₁ + r₂ − 1 = 1. A design shortfall, not a violation |
| P7 integral vs mod-2 rank discrepancy, dim ker(E mod 2) − dim ker_Z(E) | 0 in 30 runs, 1 in 6 runs (descriptive) |

- In E1 (no pair rows), no kernel vector has γ ≠ 1. At toy sizes each atom has a private prime ideal, so exact
  cancellation forces y_P = 0. E1 therefore tests the hidden-exponent extraction and exactness. E2 tests unit
  elimination.
- **Supplementary E2s** (post-prereg, disclosed: pair box 600, 18 primes, `experiments/configs/dlog_e2s.json`):
  - 18/18 correct, 0 violations;
  - 0/2611 non-torsion γ passed;
  - full unit rank in all 9 runs with γ ≠ 1, including one totally real cubic field with unit rank 2.
- The toy runs LLL on the whole returned kernel, not on the fixed-size subset of XK4, so it does not validate XK's
  asymptotic cost.
- Commands are in the prereg. Analysis: `.venv/bin/python -m src.dlog.analyze experiments/raw/dlog_e1.jsonl`.

## 10. Reviews

- **GPT-6.1-sol dlog_1 (blind, xhigh).** It independently found exact dependencies via Ge's algorithm, Lemma S_p, the
  hidden-label uniformity, the generator-free corollary, and the threshold C_LA ≤ 1 + 1/θ. Its literature table was
  merged into §1. Its one gap ("fast exact-dependence lemma" unquantified) is Proposition XK here.
- **r19 §2** reached the same route through Lenstra–Silverberg Thm 1.11, with "every positive constant" conditional
  on recursion.
- Adversarial referee (codex dlog_2), Opus proof critic and Opus literature critic: see §11, added after the reviews.

## 11. Review log

(Filled in after the reviews.)
