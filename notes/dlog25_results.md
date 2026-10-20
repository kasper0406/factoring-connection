# dlog25: rigorous L_p[2/5] discrete logarithms in F_p^* (Theorem DL × Theorem U)

Agent `dlog25`, phase 13, 2026-10-20.
- Files: `src/dlog25/{params,toy_sp,toy_xk}.py`, `tests/test_dlog25.py` (8 tests).
- Pre-registration: `notes/dlog25_prereg.md` (commit `70317f3`).
- Data: `experiments/raw/dlog25_e1.jsonl`, `experiments/raw/dlog25_e2.jsonl`, `experiments/raw/dlog25_params.txt`.
- Codex: `codex/dlog25_1.{md,out}` (blind), `codex/dlog25_2.{md,out}` (adversarial referee).

Base texts: `notes/dlog_results.md` ("[DL]", Theorem DL, Proposition XK at fixed d), `notes/lv5unif_results.md` ("[U]",
Theorem U), `notes/auditU_audit.md`, `notes/recurse_results.md` ("[rec]"), `paper2/main.tex` ("[pp]", §16 = [DL]).

Labels: **proved** (proof here, from published results) · **proved from U** (proved here, conditional on exactly the
internal, unrefereed chain on which Theorem U rests, [U] §10 "Dependencies", and on Theorem U itself) · **known**
(published, with reference) · **conjectured** · **experimentally supported**.

## Verdict

**In plain terms: yes, conditionally, in the same sense as Theorem U.** Theorem U's growing-degree relation supply
with N := p, plus Theorem U itself as the shape-test subroutine, plus Proposition XK re-checked at growing degree, gives
a Las Vegas algorithm for discrete logarithms in F_p^* in expected time L_p[2/5, c] for every prime p, with a finite c.
It is exactly as strong as Theorem U, which no outside referee has read. The DL-specific additions are proved here
and need nothing new.

1. **Theorem DL25 (proved from U).** Fix rationals a, κ, d₀ > 0 (outer tuple τ) and a′, κ′, d₀′ > 0 (shape-test
   tuple τ′). Algorithm DLOG25[τ; τ′] (§2) is Las Vegas. Given a prime p, a primitive root g and h ∈ F_p^*, it returns
   log_g h mod (p − 1) in expected time L_p[2/5, c_DL + o(1)], where
     c_DL(τ; τ′) = max( 3a,  a + s + max(a, σ),  a + s + 2^{−2/5}·c_U(τ′) ),
     s = 2κ(1 + C_γ(τ)),  σ = 2θ₀(d₀, κ)/(5a),
   and c_U(τ′) = max(2a′, a′ + s′ + max(a′, σ′), (a′ + s′)/(1 − 2^{−2/5})) is Theorem U's constant ([U] (10.3)).
   - No factorisation of p − 1, no Schirokauer maps, no class or unit group, no GRH.
   - Arbitrary base g with h ∈ ⟨g⟩: same bound (Corollary DL25′, as [DL] Cor. DL′).
   - Hence also L_p[α, o(1)] discrete logarithms for every α > 2/5.
2. **The DL constant equals the factoring constant (Lemma C25, proved; arithmetic checked exactly in
   `src/dlog25/params.py`).**
   - For every tuple τ, using the same tuple at both levels: **c_DL(τ; τ) = c_U(τ) exactly.** On this diagonal the
     cubic linear algebra never binds, for any τ:
     c_U(τ) ≥ (a + s)/(1 − 2^{−2/5}) > 4.129a > 3a (observation of codex dlog25_1).
   - For all τ, τ′: c_DL(τ; τ′) ≥ min(c_U(τ), c_U(τ′)). So DL is never cheaper than Theorem U in this cost model, and
     inf c_DL = inf c_U.
   - Grid values (one grid for both U and DL, `dlog25.params.default_grids`, not [U]'s default): for C_γ = 0, 1, 5, 20 the upper brackets are 3.5103, 4.0469, 5.0384, 6.4975,
     identical to Theorem U's; 3a = 1.65, 1.95, 2.40, 3.45 there.
   - Like Theorem U's, these numbers are illustrations: C_γ is a finite but uncomputed constant.
3. **Every ingredient transfers (§3, §4).**
   - The supply chain (FM-U, Lemma 8.5′(R′), LD*-U, Kern-U, W-U, L2-U, DS-U, CS, D-U, Lemma IR and rejections) holds
     verbatim for prime N = p. N's arithmetic enters only through (N, t) = 1 for t ≤ 2MV, N | F(g₀, g₁), ω(N) and the
     target law; there are no factor exits.
   - Split counting: Lemma S_p″ (elementary; a counting form of the classical uniqueness of rational number reconstruction, cf. Wang 1981 and von zur Gathen–Gerhard, *Modern Computer Algebra*, §5.10) replaces Lemma S″. For T uniform on F_p^*, Pr ≥ |S|²/(2(p − 1))
     once |G| ≥ 4k³. No square classes, no 2^ω.
   - Proposition XK at growing degree d ≍ D^{1/5} (Proposition XK-U): every d-dependence is polynomial. The
     Blanksby–Montgomery gap is ≍ 1/(d² log d); conjugates lie in [p^{−d}, p]; the relation bound is
     log₂R = O(dβ_x + d² + d log(kMd log p)) bits; LLL runs in dimension ≤ 2d; ≤ O(log p) bad primes need local
     maximal orders; d is odd, so μ_K = {±1}. Total L^{3a+o(1)} in the L_p[2/5] sense.
4. **Recursion.** There is no recursion on logarithms: the halves u₀, v₀ live in the same group F_p^*. The halves are
   *factored*, by FACTOR_U (Theorem U) as a black-box subroutine. So the outer level is not self-recursive, and
   Lemma F_α is used only inside Theorem U. The constant is nonetheless the same (item 2).
5. **Individual logarithms.** Every row is a descended randomised target T_i = g^{a_i}h^{b_i}, using the one-level
   descent geometry of Theorem U. So each logarithm costs a full run, L_p[2/5, c_DL]. A precomputation / cheap per-log
   split is **not established**: it would need the projection index of the torsion lattice Λ₂ onto the new rows to be
   coprime to p − 1, a unit/class-index condition we cannot control (§7).
6. **Literature (search result, §1).** The best rigorous prime-field DL bound we found is still Pomerance 1987,
   L_p[1/2, √2] (precomputation) and L_p[1/2, 1/√2] per logarithm; nothing below α = 1/2 for prime fields, for all p
   or almost all p, unconditional or under GRH. Splitting randomised targets T ≡ u/v with |u|, |v| ≈ √p and
   descending the halves is the classical *heuristic* individual-logarithm step (Blake–Fuji-Hara–Mullin–Vanstone 1984,
   cf. Galbraith, *Mathematics of Public Key Cryptography*, Ex. 15.5.10; Coppersmith–Odlyzko–Schroeppel, Algorithmica 1
   (1986) 1–15; Joux–Lercier, Math. Comp. 72 (2003) 953–967; Commeine–Semaev, PKC 2006, LNCS 3958, 174–190), and
   Pomerance 1987 already tests smoothness with a factoring subroutine (ECM). We found no rigorous analysis of the
   split-and-descend step and no prior formulation of an L_p[2/5] DL algorithm; this is a search result, not a priority
   claim. For F_{p^n}, the rigorous results below L_Q[1/2] that we found are in
   small and lower-medium characteristic, p = L_Q[β] with β < 1/2 (Kleinjung–Wesolowski: quasi-polynomial for fixed p,
   L_Q[β + o(1)] for p = L_Q[β]; Lido 2026 for p < n). For p ≥ L_Q[1/2] we found nothing rigorous even at L_Q[1/2]. Here "medium/large characteristic" means p ≥ L_Q[1/2]. **Extension fields with
   p ≥ L_Q[1/2]: the transfer does not give anything as it stands** (§8). Even the injectivity step of the split counting
   transfers only up to units (heuristically, for fixed n); the rational side becomes a number field, and Lemma 8.5′ would be
   needed for smooth ideals in residue classes modulo ideals, which is not available.
7. **Toy checks (experimentally supported; correctness only; §9).** E1: the split identity
   #{T : E(T)} = #{coprime pairs} held exactly in 8/8 exhaustive cases, and the lower bound held in every case where
   its hypothesis held. E2: Lemma XK4-U with the *explicit* κ and R certified a nonzero torsion relation in 15/15 runs at
   d = 3, 5, with 0 torsion violations. At d = 7 the toy pair box produced too few kernel vectors (0/4 runs, descriptive).

**Weakest links.** (i) Everything in Theorem U's dependency list ([U] §10; weakest: γ ≥ (log N)^{−C_γ}, i.e. LD* with
Lemma 8.5′(R′) and Drappeau's Lemme 5; auditU). (ii) Theorem U itself as a subroutine. (iii) Lemma XK4-U's bookkeeping
at growing degree (§4), the only DL-specific estimate not inherited from [DL] (its method, Log + LLL with a
Mahler-measure gap, is classical: [DL] §9b).

---

## 1. Literature

| problem | best rigorous bound found | hypotheses | reference |
|---|---|---|---|
| DL in F_p^*, all primes p | L_p[1/2, √2] precomputation, L_p[1/2, 1/√2] per logarithm | none | Pomerance, *Fast, rigorous factorization and discrete logarithm algorithms*, in Discrete Algorithms and Complexity (Academic Press 1987) 119–143 ([DL] §1 read the scan); same constant for one logarithm including precomputation (group order p − 1 known): Enge–Gaudry, Acta Arith. 102 (2002) 83–103, §5 Example 1 (read by the literature critic) |
| DL in F_p^*, split-and-descend of randomised targets | heuristic individual-log step only; no rigorous analysis found | — | Blake–Fuji-Hara–Mullin–Vanstone 1984 (via Galbraith Ex. 15.5.10, "does not lower the complexity" at L[1/2]); Coppersmith–Odlyzko–Schroeppel 1986; Joux–Lercier 2003; Commeine–Semaev 2006 |
| same | nothing with α < 1/2 or with constant < √2 found (all p, almost all p, GRH) | — | search: this note, [DL] §1, codex dlog_1, dlog25_1 |
| DL in F_p^*, heuristic | L_p[1/3, (64/9)^{1/3}] (Schirokauer 1993; Gordon 1993: 3^{2/3}) | heuristic | Gordon, SIAM J. Discrete Math. 6 (1993) 124–138; Schirokauer, Phil. Trans. R. Soc. A 345 (1993) 409–423 |
| F_{p^n}, p ≤ n^{o(n)} | L_Q[1/2, √2] | none | Lovorn Bender–Pomerance, AMS/IP Stud. Adv. Math. 7 (1998) 221–232 |
| F_{p^n}, any p, n | (pn)^{2 log₂ n + O(1)}; quasi-polynomial for fixed p; L_Q[α + o(1)] when p = L_Q[α] (stated by KW for α ∈ (0, 1/2), §1 after Thm 1.1) | none | Kleinjung–Wesolowski, J. Amer. Math. Soc. 35 (2022) 581–624, Thm 1.1 (arXiv 1906.10668) |
| F_{p^n}, p < n | (log Q)^{O(log log Q)} expected | none | Lido, Finite Fields Appl. 111 (2026) 102753 (arXiv 2206.10327), Thm 1.1; also Granger–Kleinjung–Zumbrägel, Trans. AMS 370 (2018) 3129–3145, for fields with a suitable model |
| F_{p^n}, medium/large characteristic (here: p ≥ L_Q[1/2]) | nothing at or below L_Q[1/2] found: Bender–Pomerance needs p ≤ n^{o(n)}; KW gives L_Q[α + o(1)] for p = L_Q[α], 1/2 ≤ α < 1; nothing subexponential found for fixed n ≥ 3, p → ∞ (n = 2: Lovorn 1992 thesis, unverified) | — | search (this note and the Opus literature critic) |
| factoring, all N | L_N[1/2, 1 + o(1)]; a possibly trivial congruence of squares in L_N[1/3, (64/9)^{1/3} + o(1)] | none (non-triviality only for N = pq, p ≡ q ≡ 3 mod 4, under their Conj. 7.1) | Lenstra–Pomerance, JAMS 5 (1992) 483–516; Lee–Venkatesan, J. Number Theory 187 (2018) 92–159, Thms 2.1, 2.3 (they state in §1 that they do not treat DL) |

Rows on class groups, hyperelliptic Jacobians and Ge / Lenstra–Silverberg exact relations are as in [DL] §1. The
"none found" rows are search results, not theorems.

## 2. Parameters and the algorithm DLOG25

Notation as in [U] §1 with N := p: λ = log log p, D = log p/λ, Λ = D^{2/5}λ, L^c := L_p[2/5, c] = e^{cΛ}; "L^{o(1)}"
means exp(o(Λ)) uniformly in p. n := p − 1.

**Parameters.** For the outer tuple τ = (a, κ, d₀): d, k, ν, η₀, μ, h, θ, B, X, j, Y, J, M, H, W, S, V, 𝒫 and the filter
are exactly as in [U] §1 with N := p (d odd; common V; b-window ((1 − 1/d)V, V]; X_F = (2d + 1)MHV^d). The shape-test
subroutine is FACTOR_U[τ′] for a second tuple τ′ = (a′, κ′, d₀′) (τ′ = τ is allowed).

**Row target.** M_DL := C̄ + d, C̄ := (d + 1)π(B) + 1 (a computable bound on the number of columns, §4).

```
DLOG25[τ; τ'](p, g, h):                      # run in alternation with exhaustive search x = 0, 1, ... (g^x = h)
  A := SELECT: repeat A := SAMPLE_FIELD (as [U]/[pp], N := p); if rejected continue;
               T := g^a h^b with (a, b) uniform in (Z/n)^2;  if DESCEND_ONCE(A, T) returns a row: keep A, stop
  loop (rounds, same field A):
    rows := []
    while |rows| < M_DL:
      draw (a_i, b_i) uniform in (Z/n)^2;  T_i := g^{a_i} h^{b_i} mod p
      row := DESCEND_ONCE(A, T_i)           # Gauss split T_i = u0/v0 (|u0|,|v0| <= sqrt(2p));
                                            # shape test: FACTOR_U[tau'](|u0|), FACTOR_U[tau'](|v0|) with fresh coins,
                                            #   accept iff each is a product of <= k distinct primes, all in J;
                                            # SQ_SEARCH(A, l) for each l | u0 v0 (cached per (A, l))
      if row != FAIL: store (row, a_i, b_i)
    (x*, lambda) := EXACT_DEPENDENCE_U(A, rows)          # Proposition XK-U (§4)
    A* := sum x*_i a_i,  B* := sum x*_i b_i  (mod n)
    if gcd(B*, n) = 1: X := (lambda - A*) B*^{-1} mod n;  if g^X = h: return X
```

Row identity ([DL] (2.1), with [U]'s shape, i.e. no B-smooth part in u₀, v₀): with G(P_ℓ) = σ_ℓℓn_ℓ and
G(P) ≡ g₁ψ(a − bα) (mod p),
  T ≡ ε ∏_{q≤B} q^{e_q} · g₁^{E} · ∏_ℓ ψ(a_ℓ − b_ℓα)^{e_ℓ} (mod p),  e_q = −Σ_ℓ e_ℓ v_q(n_ℓ),  E = Σ_ℓ e_ℓ.   (2.1)

## 3. Transfer of Theorem U's supply chain to N := p (proved from U)

**Where N's arithmetic enters.** Checked against [U] §§1–8 and the [pp] lemmas they cite. The prime N enters only
through: (N, t) = 1 for 0 < t ≤ 2MV; N | F(g₀, g₁) (= Nw); ω(N) (Lemma S″'s 2^ω, (P6), FM-U's primes of J dividing N);
the gcd checks (factor exits); the law of the target R; and three size or trivial conditions, all satisfied by a large
prime p: N has no prime factor ≤ B ([pp] Lemma coprime), 2X² < N (N odd), and the f′(m) rejection modulo N. This is a
reading claim over a long chain (the weakest DL-specific step, per the Opus proof critic); we found no other use.

| [U] ingredient | N := p |
|---|---|
| Lemma P (P1)–(P10), U1, CS | verbatim: they are statements about log N, d, k, B, Y, S only; (P6) is trivial (ω = 1) |
| FM-U, GM-U, FIB, SM1 | verbatim. No prime of J divides p (J ⊂ (Y, 2Y], 2Y < p). (t, p) = 1 since t ≤ 2MV < p; the residue of GM-U is a unit for the same reason |
| Lemma 8.5′(R′) | verbatim: it is a statement about Ψ(x, y; t, r) with no N |
| W-U, DS-U, L2-U, LD*-U, Kern-U | verbatim: W < p, so w ≠ p; everything else is N-free |
| D-U (γ ≥ (log p)^{−C_γ}) | verbatim, with the same C_γ = C(a, κ, d₀) (the theorem holds for every N ≥ N₀, in particular for primes) |
| Lemma IR, rejections | verbatim; the derivative rejection p ∣ f′(m) remains a rejection |
| gcd checks / factor exits | never exit: g₁, w, f_d, the content, G(P), F_A(P), u₀, v₀ and ℓ are nonzero and of absolute value < p (|F_A(P)| ≤ X_F < p because log X_F = (μ + h + dν) log p·(1 + o(1)) and μ + h + dν ≤ θ → 0) |
| Lemma S″ (R uniform on squares) | **replaced** by Lemma S_p″ below (T uniform on F_p^*) |
| K-U, select-U | **modified**: K-U with Lemma S_p″; select-U with rounds of success ≥ φ(n)/n (Lemma Sel-DL) |
| LV6-U, Theorem DD, (T), T_rec | **replaced** by Proposition XK-U (§4) and Lemmas T_p-U, HX |
| recursion (Lemma F_α) | **not used at the outer level**; FACTOR_U is called as a subroutine (§5) |

**Lemma S_p″ (proved; elementary).** Let p be prime, G a set of primes < p with |G| ≥ 4k³, and S the set of products of
exactly k distinct primes of G lying in [Y/2, Y], with 2Y² < p and |S| ≥ C(|G|, k)/(2k). Let T be uniform on F_p^*, and
E(T) the event that the Gauss-reduced primitive vector of Λ_T := {(u, v) ∈ Z² : u ≡ Tv (mod p)} is ±(u, v) with
u, v ∈ S and gcd(u, v) = 1. Then
  #{T : E(T)} = P′ := #{(u, v) ∈ S² : gcd(u, v) = 1} ≥ |S|²/2,  so  Pr[E] ≥ |S|²/(2(p − 1)).
*Proof.*
- **Shared primes** ([rec] Lemma S″ step 2, verbatim): for u ∈ S, the number of v ∈ S sharing a prime with u is
  ≤ k·C(|G| − 1, k − 1) = (k²/|G|)C(|G|, k) ≤ (2k³/|G|)|S|. So P′ ≥ |S|²(1 − 2k³/|G|) ≥ |S|²/2.
- **Injectivity.** Map a coprime pair to T := uv^{−1} mod p (v < p is a unit). If (u, v), (u′, v′) give the same T, then
  uv′ ≡ u′v and |uv′ − u′v| < Y² < p, so uv′ = u′v; coprimality and positivity give (u, v) = (u′, v′).
- **Primitive vector** ([rec] step 5): for such T, the shortest vector s of Λ_T has |det(s, (u, v))| ≤ ‖s‖‖(u, v)‖ ≤ 2Y²
  < p, and every determinant of two lattice vectors is ≡ 0 (mod p), so (u, v) = m·s with m = ±1 by primitivity. So E(T)
  holds for the T of every coprime pair. Conversely E(T) produces a coprime pair with T = uv^{−1}. Hence #{T : E} = P′.
- Divide by |F_p^*| = p − 1. ∎
*Remark.* Compared with Lemma S″ there are no square classes (T is uniform on the whole group), so h = 1. E1 checks
the identity exhaustively (§9).

**Theorem K-U for DL (proved from U).** Let A be a kept field with δ_A ≥ γ/4. One call of DESCEND_ONCE(A, T), T uniform,
returns a row with probability ≥ (δ_A/(4 log p))^{2k}/(32k²).
*Proof.* [U] Theorem K-U verbatim with Lemma S_p″ in place of Lemma S″: |𝒟_A| ≥ γπ_J/4 ≥ 4k³ by (P6); one of the
k windows holds |𝒮| ≥ C(|𝒟_A|, k)/k ≥ C(|𝒟_A|, k)/(2k) products, with 2Y′² < p; Pr ≥ |𝒮|²/(2(p − 1)) ≥ |𝒮|²/(2p); the
rest of the computation is [U]'s. On E the shape test passes (FACTOR_U returns the true factorisations, which consist
of exactly k distinct primes of 𝒟_A ⊂ J on each side), and every special-q search returns an atom ([pp] Prop BS). ∎

**Lemma T_p-U (hidden labels; proved).** (The transcript below is that of the *untruncated core* algorithm after
target construction: drawing (a_i, b_i) and computing T_i are not functions of 𝒯_j, and the exhaustive-search wrapper is
analysed by domination, T_wrapped ≤ 2T_core + O(1), never by conditioning on the wrapper still running; codex dlog25_2.) (a, b) ↦ (log_g T, b) = (a + Xb, b) is a bijection of (Z/n)², so T is uniform
on F_p^* and, given T, b is uniform on Z/n and independent of T. Let 𝒯_j be: the field, all targets of round j and their
acceptance outcomes, every coin of every FACTOR_U call, SQ_SEARCH and EXACT_DEPENDENCE_U in round j, and everything from
earlier rounds; it excludes the label coins of round j and step "A*, B*". Every computation of round j before that step
is a function of 𝒯_j: FACTOR_U is Las Vegas and its *output* (the factorisation of |u₀|, |v₀|) is a function of T alone,
its running time a function of T and its own coins; acceptance is a function of the factorisations, the cached atoms and
the SQ coins. Hence, conditional on 𝒯_j, the labels b_i of the accepted rows of round j are i.i.d. uniform on Z/n. ∎
(This is [DL] Lemma T_p with the shape-test coins added to 𝒯_j; the recursive structure inside FACTOR_U is irrelevant,
because only its output enters the row and that output is deterministic.)

**Lemma Sel-DL (cost of selection and rounds; proved).** Let p_sel be the probability that one trial (SAMPLE_FIELD +
DESCEND_ONCE) returns a row, C_trial a bound on the conditional expected cost of one trial or one DESCEND_ONCE given the
past, and T_post a bound on the expected cost of EXACT_DEPENDENCE_U. Then
  E[cost of DLOG25] ≤ (1 + log₂ n)·((M_DL + 1)·C_trial/p_sel + T_post) + T_S,
where T_S = B^{1+o(1)} is the cost of S (Lemma CS), and p_sel ≥ (γ/4)(γ/(16 log p))^{2k}/(32k²).
*Proof.* p_sel: [U] select-U(i) with K-U for DL; there are no factor exits. Selection: Wald, ≤ C_trial/p_sel. Kept field:
size-biased law, E[1/ς_{A*}] ≤ 1/p_sel ([U] select-U(ii), verbatim). A round in field A costs ≤ M_DL·C_trial/ς_A + T_post;
C_trial includes the O(log³ p) cost of drawing (a_i, b_i) and computing T_i, paid by every draw. By Lemmas T_p-U and HX, each round succeeds with conditional probability
≥ φ(n)/n ≥ 1/(1 + log₂ n) given all earlier rounds, so the number of rounds is dominated by a geometric variable with
mean ≤ 1 + log₂ n; Wald's identity (each round's conditional expected cost is bounded as above given A) gives the
bound. ∎

## 4. Proposition XK at growing degree (Proposition XK-U; proved)

Fix a kept field A: f irreducible of odd degree d with content 1, 0 < |f_d| < p, m a simple root of f mod p, and
‖f‖_∞ ≤ 3MH ([pp] Lemma GB_w with (P2)). Rows i = 1..M, M = M_DL, with labels. Notation of [DL] §4: K = Q(α), t = f_d,
𝔞 = O_K + αO_K, 𝔍_P = (a − bα)𝔞^{−1}, Λ₁, γ_x, y_P(x), Λ₂.

**What is degree-free.** Lemmas XK0(i)–(ii), XK1, XK2, XK3, XK5 and HX of [DL] contain no constant depending on d:
their proofs are pure algebra (Gauss's content lemma, Dedekind–Kummer, the reduction map at the degree-one prime
𝔭 = (p, α − m), saturation of kernels, μ_K = {±1} for d odd). Their numerical hypotheses hold in the U design:
- |F_A(P)| ≤ X_F < p and p ∤ t, so every a − bα is a 𝔭-unit (N𝔍_P < p);
- |G(P)| ≤ 2MV < p, g₁ < p, and all q ≤ B are < p, so (2.1) holds in F_p^*;
- d is odd ([U] §1: d is the least odd integer ≥ max(3, d₀D^{1/5})), so W(x) = ±1 on Λ₂.
Columns: one per q ≤ B, one for g₁, and one per prime ideal above some r ≤ B that occurs; at most d per r, so
C_col ≤ (d + 1)π(B) + 1 = C̄ = L^{a+o(1)} (d ≤ log p). Rank Λ₂ ≥ M − C_col − r_U ≥ 1 for M = C̄ + d.

**Lemma XK0-U (bad primes; proved).** The primes r ≤ B with r | t·disc(f) number ≤ log₂|t·disc f| = O(log p), because
log|disc f| ≤ d log d + (2d − 2) log(√(d + 1)·3MH) = O(d(μ + h) log p + d log d) = O(log p) (d(μ + h) = O(dθ) = O(1)).
For each, an r-maximal order, the primes above r and the valuations v_𝔯(𝔍_P) are computed in expected time
poly(d, log p) (an overorder maximal at r in deterministic polynomial time: Buchmann–Lenstra, J. Théor. Nombres
Bordeaux 6 (1994) 221–260, Thm 1.2 / Alg. 6.1 and Prop. 6.2, a Round-2-type iteration, cf. Cohen GTM 138 Alg. 6.1.8;
since f need not be monic, the input order is Z[tα] (or the standard order of f), as in [DL]'s Z_(r)[α] remark;
decomposition and valuations: Cohen §§6.2, 4.8.3 (algorithm numbers 6.2.9, 4.8.17 not verified against the book); the only randomness is Cantor–Zassenhaus factoring over F_r, and r ≤ B; the r-adic
precision is bounded by v_r of the discriminant of the actual starting order, e.g. disc Z[tα] = t^{(d−1)(d−2)}disc f, whose
logarithm is O(d log p); codex dlog25_2). Finding them: t·disc f mod r for all r ≤ B, Õ(dB). Total, including the
valuations of the ≤ 2kM atom occurrences at the bad primes: ≤ 2kM·O(log p)·poly(d, log p) + Õ(dB) = L^{a+o(1)}. ∎

**Lemma XK4-U (unit elimination at growing degree; proved).** Let s := r₁ + r₂ ≤ d, r_U := s − 1, and x₁, …, x_s ∈ Λ₁
independent with entries of at most β_x bits. One finds a nonzero n ∈ Z^s with Σ n_jx_j ∈ Λ₂ in time
poly(d)·Õ(M·(P + β_x) + κ²) with κ, P = O(dβ_x + d² + d·log(kMd log p)).
*Proof (every d-dependence explicit).*
- **Conjugate bounds.** For an atom P = (a, b) (|a|, |b| ≤ V, F(a, b) ≠ 0) and every embedding σ:
  |σ(a − bα)| ≤ U := V(2 + 3MH) < p (Cauchy: |σα| ≤ 1 + ‖f‖_∞/|t| ≤ 1 + 3MH; and 8MHV ≤ X_F < p), and
  |σ(a − bα)| = |N(a − bα)|/∏_{σ′≠σ}|σ′(a − bα)| ≥ (|F(a, b)|/|t|)·U^{−(d−1)} ≥ p^{−d}. So every coordinate of
  Log(a − bα) has absolute value ≤ 2d log p, and ‖Log(a − bα)‖₂ ≤ 2d^{3/2} log p.
- **Gap.** Blanksby–Montgomery (Acta Arith. 18 (1971) 355–369): a nonzero algebraic integer of degree d′ that is not a
  root of unity has Mahler measure > 1 + 1/(52d′ log(6d′)). For a unit z ∈ K of degree d′ | d,
  log M(z) = (d′/d)Σ_σ log⁺|σz| ≤ d′·‖Log z‖_∞ ≤ d‖Log z‖_∞ (complex coordinates carry the factor 2). So every
  nonzero w ∈ Λ_U := Log(O_K^×) has ‖w‖₂ ≥ ‖w‖_∞ ≥ β_d := log(1 + 1/(52d log 6d))/d ≥ 1/(105d² log 6d), and Log z = 0
  iff z ∈ μ_K (Kronecker).
- **Small relations exist.** v_j := Log γ_{x_j} = Σ_P y_P(x_j) Log(a_P − b_Pα) ∈ Λ_U, with
  ‖v_j‖₂ ≤ 2^{β_x}·2kM·2d^{3/2} log p (≤ 2kM atom occurrences, each |coefficient| ≤ 2^{β_x}). Fix (for the proof only) an
  LLL-reduced basis u₁, …, u_{r_U} of Λ_U: ‖u_i*‖ ≥ 2^{−(i−1)/2}λ₁(Λ_U) ≥ 2^{−r_U/2}β_d. Back-substitution with
  |μ_{ij}| ≤ 1/2 gives integer coordinates |c_{ji}| ≤ (3/2)^{r_U}·2^{r_U/2}‖v_j‖₂/β_d =: C₁ (the fixed-d text [DL] wrote
  2^{O(d²)}; 2^{O(d)} suffices, and either is polynomial in bit length). The r_U × s integer matrix (c_{ji}) has a
  nonzero kernel vector n₀ with ‖n₀‖_∞ ≤ R := max(1, r_U!·C₁^{r_U}) (Cramer on a maximal nonsingular minor; if all v_j = 0,
  n₀ = e₁). Thus
    log₂R ≤ r_U log₂ r_U + r_U(r_U log₂(3/2) + r_U/2 + β_x + log₂(4kMd^{3/2} log p) + log₂(105d² log 6d))
          = O(dβ_x + d² + d log(kMd log p)).
- **LLL with the explicit precision.** b_j := (e_j, ⌊2^κ v_j⌉) ∈ Z^{s + s}, coordinates within 1 of 2^κv_j. If n ∈ Rel,
  ‖Σn_jb_j‖ ≤ (1 + s√s)‖n‖; if n ∉ Rel, ‖Σn_jb_j‖ ≥ 2^κβ_d − s√s‖n‖. Since R bounds ‖n₀‖_∞, ‖n₀‖₂ ≤ √s·R, so the
  lattice contains a vector of norm ≤ √s(1 + s√s)R and LLL's first vector has norm ≤ R′ := 2^{(s−1)/2}√s(1 + s√s)R; its
  coefficient vector has norm ≤ R′. With 2^κβ_d > (1 + s√s)R′ it lies in Rel, and it is nonzero. So
  κ := ⌈log₂((1 + s√s)R′/β_d)⌉ + 1 = log₂R + O(d + log d) + log₂(1/β_d) = O(dβ_x + d² + d log(kMd log p)).
  (Corrected after the Opus proof critic, #2: the first draft used R as a Euclidean bound. The toy's
  `toy_xk.explicit_kappa` uses the draft's R′ without √s but adds 2 to the ceiling, so it gives
  2^κβ_d ≥ 4(1 + s√s)²2^{(s−1)/2}R ≥ (1 + s√s)R′ whenever √s ≤ 4; every E2 run has s ≤ 4, so E2 used a valid κ.)
- **Precision ledger.** Each coordinate of each Log(a_P − b_Pα) to absolute precision 2^{−P}, P := κ + β_x +
  log₂(2kM) + O(log d), suffices for rounding, accumulation and certification ([DL] Lemma XK4 ledger, with the sum over
  ≤ 2kM occurrences and s ≤ d coordinates). By the conjugate lower bound, absolute precision P + log₂V + d log₂p + O(1)
  on the roots of f gives this (relative error of σ(a − bα) ≤ V·ε·p^d). Root approximation of a degree-d integer
  polynomial with O(log p)-bit coefficients to b bits costs Õ(d²(b + d log p + d)), near-linear in b (Pan 2002,
  Thm 2.1.1 and Cor. 2.1.2, after scaling the roots into the unit disk; codex dlog25_2); with b = P + log₂V + d log₂p + O(1)
  this is L^{a+o(1)} (A. Schönhage, *The fundamental theorem of algebra in terms of computational
  complexity*, preliminary report, Math. Inst. Univ. Tübingen, 1982; Pan, J. Symbolic Comput. 33 (2002) 701–733). The roots must be *certified* (isolated), to separate real from complex roots and to pair the
  complex embeddings consistently across atoms: f is squarefree (irreducible), its root separation is
  ≥ 2^{−O(d log d + d log‖f‖)} = 2^{−O(d log p)} (Mahler), and P exceeds this, so isolating enclosures (e.g. by
  Newton–Kantorovich or Gershgorin tests on the approximations) cost poly(d, log p, P) (Opus proof critic, #3).
- **Cost.** Roots: Õ(d²(P + d log p)). Logs of ≤ 2kM atoms in s coordinates: d·2kM·Õ(P). The sums: Õ(dM(P + β_x)).
  LLL (Nguyen–Stehlé L², SIAM J. Comput. 39 (2009) 874–903) in dimension ≤ 2d with O(κ)-bit entries: poly(d)·Õ(κ²).
- **Certification.** Accept n only if ‖Σn_jv_j‖_∞ < β_d/2 is certified by ball arithmetic: a true relation passes by the
  precision, a non-relation cannot pass by the gap. ∎

**Proposition XK-U (proved).** EXACT_DEPENDENCE_U = [DL] EXACT_DEPENDENCE with Lemmas XK0-U, XK4-U and M = M_DL is Las
Vegas, outputs a primitive x* ∈ Λ₂ and λ = log_g W(x*) ∈ {0, n/2}, and has expected cost L^{3a+o(1)} (L = L_p[2/5]).
*Proof.* Correctness: [DL] Lemmas XK1–XK5 (degree-free) with XK0-U, XK4-U. Cost, step by step:
- E exactly: each row has ≤ 2k(1 + log₂(2MV)) + 1 + 2k·log₂X_F ≤ 5k log₂ p nonzero entries, each of absolute value
  ≤ 2k log₂ p (k ≤ log p); the valuations are read off the factorisations of G(P) and F(P) (already computed by the
  batch smoothness of SQ_SEARCH) plus XK0-U. Cost M·poly(k, d, log p) + L^{a+o(1)}.
- Rank profile modulo a random prime q ∈ [Q, 2Q], Q = M(log p)³: a fixed maximal nonsingular minor has
  log|det| ≤ ρ·log(√(5k log p)·2k log p) = O(ρ log log p) (Hadamard; k, d ≤ log p), so ≤ O(ρ log log p/log Q) primes in
  [Q, 2Q] divide it, against ≫ Q/log Q primes there: Pr[the minor survives] ≥ 1/2. Dense elimination: O(M·C̄·ρ) = L^{3a+o(1)}.
- Dixon p-adic lifting for d right-hand sides: Õ(dρ³) = L^{3a+o(1)}; outputs β_x = O(ρ log log p) = L^{a+o(1)} bits.
- Verification x_iᵀE = 0: nnz(E)·d·β_x = L^{2a+o(1)}; a failure resamples q (expected ≤ 2 restarts, [DL]).
- Lemma XK4-U with this β_x: κ, P = O(dβ_x + d² + d log(kMd log p)) = L^{a+o(1)}, cost L^{2a+o(1)}.
- W(x*), A*, B*: M exponentiations with L^{a+o(1)}-bit exponents, L^{2a+o(1)}.
Every factor poly(d, k, log p) is (log p)^{O(1)} = L^{o(1)} (in fact exp(O(log log p)) = exp(o(Λ))). ∎

So the answer to "does the cost stay L^{o(1)} relative to L[2/5]": yes. The growing degree costs only poly(d) factors
and O(d)-fold larger bit lengths (κ ≍ dβ_x instead of β_x), never an exponent. The cubic linear algebra L^{3a} is the
only post-processing term that is not L^{2a}; its effect is computed in §6.

## 5. Shape test, recursion, and why DL has no log recursion

- **The shape test is factoring.** DESCEND_ONCE needs the complete factorisations of |u₀|, |v₀| ≤ √(2p). DLOG25 calls
  FACTOR_U[τ′] on each, with fresh coins. Theorem U's proof bounds the sup form F_U(x) := sup_{2≤m≤x} E[cost of
  FACTOR_U[τ′](m)] ([U] §10: Lemma F_α is applied to this F), so for every integer m ≤ √(2p), E[cost of FACTOR_U[τ′](m)] ≤
  F_U(√(2p)) ≤ L_{√(2p)}[2/5, c_U(τ′) + o(1)] ≤ e·L_p[2/5, 2^{−2/5}c_U(τ′) + o(1)] (the scaling step of [U] Lemma F_α:
  Λ(√(2p)) ≤ 2^{−2/5}Λ(p) + o(1)). This bound is uniform in m, and conditional on the past (fresh coins), so it is a valid
  C_trial term in Lemma Sel-DL ([rec] Lemma W covers adaptively chosen sub-calls).
- **No recursion on logarithms.** The logarithms of the halves u₀, v₀ are logarithms in the same group F_p^*, so a
  recursive DL call would be on the same p: no size reduction. The only recursion is inside FACTOR_U (Theorem U's own,
  via Lemma F_α). DLOG25 itself is one level.
- **Why the constant is still the factoring constant.** See Lemma C25 (§6): a one-level outer loop around a factoring
  subroutine with constant c* costs ≥ min(c_U(τ), c*) for the same reason that the fixed point of Lemma F_α is c_U.

## 6. Theorem DL25 and the cost (proved from U)

**Theorem DL25.** Fix rationals a, κ, d₀, a′, κ′, d₀′ > 0. DLOG25[τ; τ′] is Las Vegas, and for every prime p and
primitive root g it returns log_g h mod (p − 1) in expected time L_p[2/5, c_DL(τ; τ′) + o(1)], with
  c_DL(τ; τ′) = max( 3a,  a + s + max(a, σ),  a + s + 2^{−2/5}c_U(τ′) ),  s = 2κ(1 + C_γ(τ)),  σ = 2θ₀(d₀, κ)/(5a). (6.1)
*Proof.*
- **Correctness.** Every output is verified (g^X = h); by Lemmas XK1–XK5, HX and T_p-U, A* + XB* ≡ λ (mod n) for the
  true X, so a round with gcd(B*, n) = 1 returns X. The alternation with exhaustive search makes the algorithm
  terminate for small p (where structural guards may fail) and costs a factor ≤ 2.
- **Cost.** Lemma Sel-DL with:
  - 1/p_sel ≤ 4γ^{−1}(16 log p/γ)^{2k}·32k² = L^{s+o(1)} (γ ≥ (log p)^{−C_γ}, 2k log log p = 2κΛ(1 + o(1)); [U] (10.1));
  - M_DL = (d + 1)π(B) + 1 + d = L^{a+o(1)};
  - C_trial ≤ 2F_U(√(2p)) + 2k·T_SQ + poly(d, log p), T_SQ = L^{max(a,σ)+o(1)} ([U] §10), so
    C_trial ≤ L^{max(2^{−2/5}c_U(τ′), max(a, σ))+o(1)};
  - T_post = L^{3a+o(1)} (Proposition XK-U); T_S = L^{a+o(1)}; 1 + log₂ n = L^{o(1)}.
  Total: L^{max(a + s + 2^{−2/5}c_U(τ′), a + s + max(a, σ), 3a) + o(1)}. ∎

**Corollary DL25′ (arbitrary base; proved from U).** For g, h ∈ F_p^* with h ∈ ⟨g⟩, log_g h mod ord(g) in the same
expected time: [DL] Corollary DL′ verbatim (random auxiliary bases t, primitive with probability ≥ 1/(1 + log₂ n),
doubling budgets with Markov, every output verified, inconsistent congruences discarded).

**Corollary.** For every α > 2/5, DL in F_p^* in expected L_p[α, o(1)]: L_p[2/5, c] = L_p[α, o(1)].

**Lemma C25 (the DL constant is the factoring constant; proved).** Write T₁(τ) := a + s + max(a, σ) and
c_U(τ) = max(2a, T₁(τ), (a + s)/(1 − 2^{−2/5})).
- (i) c_DL(τ; τ) = c_U(τ) for every τ; in particular the term 3a never binds when the same tuple is used at both
  levels. (Off the diagonal it can bind: for C_γ = 20, τ = (3, 1/100, 3/10), τ′ = (23/20, 1/100, 3/10) the terms are
  9, 6.42, 8.344, so 3a binds; codex dlog25_2. This does not affect (iii).)
- (ii) For all τ, τ′: c_DL(τ; τ′) ≥ min(c_U(τ), c_U(τ′)).
- (iii) Hence inf_{τ,τ′} c_DL = inf_τ c_U.
*Proof.* (i) ≤: (a + s) ≤ (1 − 2^{−2/5})c_U(τ), so a + s + 2^{−2/5}c_U(τ) ≤ c_U(τ); T₁ ≤ c_U(τ); and
3a < a/(1 − 2^{−2/5}) ≤ (a + s)/(1 − 2^{−2/5}) ≤ c_U(τ), since 1/(1 − 2^{−2/5}) > 4.129 > 3 and s ≥ 0 (this observation is
from codex dlog25_1). ≥: by (ii) with τ′ = τ. (ii): let x := c_DL(τ; τ′) and suppose x < c_U(τ′). Then
a + s ≤ x − 2^{−2/5}c_U(τ′) < (1 − 2^{−2/5})x, so (a + s)/(1 − 2^{−2/5}) < x; also T₁ ≤ x and 2a < 3a ≤ x; so
c_U(τ) ≤ x. (iii): (i) and (ii). ∎

**Numbers** (`.venv/bin/python -m src.dlog25.params`, exact rationals, 2^{±2/5} bracketed; the same grid for U and DL,
`dlog25.params.default_grids`, which differs from [U]'s default grid: κ ∈ {i/200 ≤ 1/5} ∪ {2^{−j}}, a ≤ 2 — this is why
C_γ = 0 gives 3.51035 here against [U]'s 3.517; the
subroutine uses the grid-optimal τ′ for the same C_γ; the reported value is the upper bracket):

| C_γ | inf c_U (grid, upper) | c_DL (grid, upper) | outer τ = (a, κ, d₀) | 3a | binding term |
|---|---|---|---|---|---|
| 0 | 3.51035 | 3.51035 | (11/20, 3/20, 23/20) | 1.65 | shape test (special-q within 10^{−5}) |
| 1 | 4.04692 | 4.04692 | (13/20, 2/25, 4/5) | 1.95 | special-q |
| 5 | 5.03838 | 5.03838 | (4/5, 7/200, 1/2) | 2.40 | shape test |
| 20 | 6.49754 | 6.49754 | (23/20, 1/100, 3/10) | 3.45 | special-q |

So the linear-algebra term ("the 3a term takes over" at L[1/2] for θ > 1/2, [DL] §5) never binds at L[2/5] on the
diagonal τ′ = τ, which already attains inf c_DL (Lemma C25): the recursion already makes c_U ≥ 4.129a. At L[1/2] there is no recursion factor 1/(1 − 2^{−1/2})
in the fixed-shape-test constant, which is why 3a could bind there. As in [U], the table entries are grid values and C_γ
is not computed.

## 7. Individual logarithms and the target randomisation

- **Target randomisation.** T_i = g^{a_i}h^{b_i} with (a_i, b_i) uniform is uniform on F_p^* and carries an independent
  uniform hidden label b_i (Lemma T_p-U). Theorem K-U's split law needs only "T uniform on F_p^*" (Lemma S_p″). The
  descent geometry is Theorem U's: one Gauss split into halves ≤ √(2p), k ≍ D^{2/5} primes of size L[3/5] per half,
  each descended by one special-q search into B = L[2/5]-smooth rational and algebraic values.
- **Every logarithm costs a full run.** In this design there is no separate "individual-log phase": the target h enters
  every row. A precomputation/descent split (rows with b = 0 precomputed once; one new row per logarithm) would need,
  for the new coordinate(s), the projection of the saturated torsion lattice Λ₂ to have index coprime to n; that index is
  governed by unit and class indices of K and cannot be controlled (this is exactly what Schirokauer maps or virtual logs
  address, at the cost of factoring n and the conditions of [DL] §6). With ≥ C̄ + d fresh h-rows the precomputation is
  useless. **Not established.** At L[2/5] this matters less than at L[1/2]: both phases would be L[2/5, Θ(1)].
- What can be shared between logarithms: the field. SELECT may use targets g^a only (no h), so the kept field, the
  factor-base tables and the special-q cache are target-independent and reusable (codex dlog25_1); the per-log bound
  stays L_p[2/5, c_DL].

## 8. Extension fields F_{p^n} in medium characteristic

- **Known (rigorous).** Bender–Pomerance L_Q[1/2, √2] for p ≤ n^{o(n)}; Kleinjung–Wesolowski (pn)^{2 log₂n + O(1)}. For
  Q = p^n and p = L_Q[β, c] with fixed 0 < β < 1 this is exp((2c(1 − β)/log 2 + o(1))(log Q)^β(log log Q)^{2−β})
  (codex dlog25_1), i.e. L_Q[β + o(1)] but not L_Q[β, C] for a fixed C; it beats any L_Q[2/5, C] only for β < 2/5, and
  is below L_Q[1/2] only for β < 1/2. For p ≥ L_Q[1/2] (medium and large characteristic) we found nothing rigorous
  at or below L_Q[1/2]: Bender–Pomerance needs p ≤ n^{o(n)}, and for fixed n ≥ 3, p → ∞ nothing subexponential
  (Enge–Gaudry p. 97 remark; n = 2: Lovorn's thesis, unverified) (search result).
- **Does the transfer give anything? Not as it stands.**
  - What plausibly transfers (heuristic; fixed n): one step of the Dixon split. Represent F_{p^n} = Z[x]/(p, φ), φ
    monic of degree n irreducible mod p, and let Λ_T := {(u, v) ∈ Z[x]_{<n}² : u ≡ Tv mod (p, φ)}, of rank 2n and
    covolume p^n. The analogue of Lemma S_p″'s *injectivity* holds when the coefficients of u, v are ≤ Y with
    N_{Q(φ)}(uv′ − u′v) < p^n, since uv′ − u′v ∈ (p, φ(θ)) forces the norm to be ≥ p^n unless it is 0. The
    "injectivity" only identifies the *ratio* u/v ∈ Q(φ)^×: unit multiples (εu, εv) give the same ratio (e.g. in
    Q(√2), (3 + √2, 5 + 2√2) and ε = 3 + 2√2; codex dlog25_2), so unit normalisation and principality issues remain. The
    *primitive-vector* step does not transfer as is: Λ_T contains θ^i·(u, v), so there is no unique shortest vector
    (module-lattice reduction would be needed), and coprimality needs ideal-theoretic care since Z[θ] need not be a UFD
    (Opus proof critic, #7).
  - What fails: the halves are elements of a degree-n number field Q(φ), not integers. The rational side of the relation
    supply becomes an algebraic side, and FM/LD* need Lemma 8.5′ (restricted BV for smooth numbers, Drappeau Lemme 5) for
    **smooth ideals in residue classes modulo ideals of Q(φ)**, which is not available. (The shape test itself is not
    an obstacle: an explicit ideal of norm ≈ Q^{1/2} can be factored by factoring its integer norm with Theorem U and
    decomposing locally; codex dlog25_2.) The outstanding obstacles are the smooth-ideal distribution and the field
    construction. The pair (Q(φ), K_A) would also need a common-root construction F(g₀, g₁) =
    p·w over Q(φ) (a "conjugation"-type pair of fields), whose field law is not covered by Lemma IR / FM-U.
  - Two simpler objections (codex dlog25_1): putting N := p^n in Theorem U's design gives arithmetic in Z/p^n, not in
    F_{p^n}; and integer ratios u₀/v₀ reduce into the prime field, so they cannot represent uniform targets of F_{p^n}.
    A number-field reduction onto F_{p^n} needs a prime of residue degree n, i.e. field degree ≥ n, which for n ≫ d
    exceeds this design's degree budget.
  - So in medium characteristic this is a research programme, not a corollary. Not pursued.

## 9. Toy experiments (pre-registered, `notes/dlog25_prereg.md`, commit 70317f3)

Commands as in the prereg; E1 ≈ 23 s, E2 ≈ 9 s.

| prediction | result |
|---|---|
| P1 E1 identity #{T : E(T)} = P′ | **8/8** exhaustive cases (p up to 4 000 037, k ∈ {1, 2, 3}) |
| P2 E1 bound where \|G\| ≥ 4k³ | holds in both such cases (k = 1): Pr[E]/(|S|²/(2(p − 1))) = 1.91, 1.95. In the 6 cases without the hypothesis the identity still holds; the ratio ranged 0–1.63 (two cases had |S| ≤ 2 and P′ = 0) |
| P3 E2 certification at explicit κ | **15/15** ok runs (d = 3: 8, d = 5: 7) certified a nonzero relation with the LLL first vector |
| P4 E2 torsion | **0** violations: every certified x* gave ψ_q(γ) = ±1 with a consistent sign at 5 auxiliary 40-bit primes |
| P5 E2 units/rank | 0 non-unit kernel vectors; unit_rank_implied = r_U in 15/15 (the s = r_U + 1 vectors always spanned the full unit rank image); 0 non-relation LLL rows passed the auxiliary test |
| P6 coverage | d = 3: 8/8, d = 5: 7/8 (one field had 165 relations for 181 columns). d = 7 (descriptive): 0/4, kernels of dimension 1–3 < s |
| P7 descriptive | κ = 31–55 at d = 3, 55–95 at d = 5; log₂R = 13–36 (d = 3) and 34–73 (d = 5), while log₂‖n‖ ≤ 8.9: the explicit bound is very loose but valid |

- E1 is a check of an exact identity; its hypothesis |G| ≥ 4k³ is only reachable at k = 1 at exhaustive sizes.
- E2 uses NFS-type pair elements of small fields, not descended rows (disclosed), and toy d ≤ 5 only. It checks that
  the *explicit* κ of Lemma XK4-U is sufficient; it does not test the asymptotic cost. Seeds 1 and 2 at d = 3 drew the
  same field (identical records; disclosed, not a failure).

## 10. Reviews

**GPT-6.1-sol blind round `codex/dlog25_1` (xhigh; before it saw any draft).** Verdict: "YES (L[2/5]), conditional on
the stated internal chain; the growing-degree XK step introduces no larger exponent." It independently derived the
same algorithm, the same coprime-pair replacement for Lemma S″ (bound (3) = Lemma S_p″), the same cost formula
c_DL = max(3a, a + s + max(a, σ, 2^{−2/5}c_F)), c_DL = c_U at the same tuple, XK at growing degree with
log R = O(dβ_x + d log M + d log(kd log p) + d³) (it kept the 2^{O(d²)} coordinate bound of [DL]; we use 2^{O(d)}), the
per-round hidden-label conditioning, and the same negative answers on per-log precomputation and medium characteristic.
Adopted from it: (a) c_U ≥ (a + s)/(1 − 2^{−2/5}) > 4.129a > 3a, which upgrades Lemma C25 from "inf c_DL ∈ [c*, 3c*/2],
= c* on the grid" to c_DL(τ; τ) = c_U(τ) for every τ; (b) the precise KW exponent for p = L_Q[β, c]; (c) the two
elementary extension-field objections in §8; (d) the remark that p ∤ disc f is not needed globally (only a simple root
m), which is already how SAMPLE_FIELD is written; (e) the reusable target-independent field (selected with powers of g
only, §7). It read Buchmann–Lenstra 1994 Thm 1.2, Alg. 6.1/Prop. 6.2 (local maximal orders in polynomial time) and
Dixon 1982.

**Opus proof critic (read-only, independent; reviewed the pre-merge draft): ACCEPT WITH FIXES, no FATAL/MAJOR.** It
re-checked §3 (where N enters; Lemma S_p″; K-U; T_p-U; Sel-DL), §4 line by line (Cauchy and conjugate bounds, BM gap
with log M(z) ≤ d′‖Log z‖_∞, LLL coordinate bound, Cramer, κ, ledger, Hadamard/random prime, Dixon, μ_K), the scaling
of the subroutine cost, (6.1) and Lemma C25; reproduced the §6 table; found c_U/(3a) ≥ 1.382 at all 211 680 grid points;
and re-ran the S_p″ identity in 5 cases. Fixes applied: (1) 3a never binds for any tuple (same as codex's observation,
already merged); (2) R is an ℓ^∞ bound, so R′ gains √s (fixed in XK4-U; the toy κ remains valid for s ≤ 16); (3) certified
root isolation via the root-separation bound; (4) per-atom bad-prime valuations in XK0-U's total; (5) the grid differs
from [U]'s default (explains 3.51035 vs 3.517); (6), (7) scope wording (lower-medium characteristic; only the injectivity
step transfers to F_{p^n}, heuristically); (8) three size conditions on N listed in §3; (9) exponentiation cost in
C_trial, sup form of F_U; (10) test count. Weakest DL-specific step in its view: the §3 reading claim, then XK4-U.

**Opus literature critic (read-only, independent): OK WITH FIXES.** It read Enge–Gaudry 2002 (§1, §5 Ex. 1, p. 97),
Galbraith ch. 15, Lee–Venkatesan pp. 1–3, Buchmann–Lenstra 1994 (Thms 1.1–1.4, Alg. 6.1, Prop. 6.2), Kleinjung–Wesolowski
Thm 1.1, Lido Thm 1.1 and the dBPW abstract; Pomerance 1987 and Bender–Pomerance were image-only scans. It found nothing
beating Pomerance 1987 for prime fields. Fixes applied: three novelty overclaims (Lemma S_p″ is a counting form of
rational-reconstruction uniqueness; XK4-U's method is classical; "new here" → "proved here, no prior formulation found
by us"); the split-and-descend step is the classical heuristic individual-log step (BFMV 1984, COS 1986, Joux–Lercier
2003, Commeine–Semaev 2006), now cited, with no rigorous analysis found; Lee–Venkatesan row; Enge–Gaudry confirmation of
√2; KW's α-range; for p ≥ L_Q[1/2] nothing rigorous even at L_Q[1/2]; Lido and GKZ rows; Gordon's constant 3^{2/3};
Buchmann–Lenstra pinpoint and the non-monic input order; Cohen algorithm numbers 6.2.9/4.8.17 flagged unverified;
Schönhage 1982 full reference.

**GPT-6.1-sol adversarial referee `codex/dlog25_2` (xhigh): ACCEPT WITH FIXES, "no surviving FATAL or MAJOR gap".** It
confirmed the transfer to N = p, the size-biased selection and Wald accounting, Lemma S_p″'s exact identity and K-U's
windows, the hidden labels with the randomized subroutine (round success exactly φ(n)/n given the core transcript), XK-U
at growing degree (no hidden exp(d)-type cost; R is never enumerated), (6.1), Lemma C25 and the table, and the scope of
§§7–8; it re-ran the 8 tests and its own scratch checks. Fixes applied (all MINOR): √s in R′ (as the proof critic);
near-linear root-approximation cost via Pan Thm 2.1.1/Cor. 2.1.2; r-adic precision from the starting order's
discriminant (O(d log p)); the hidden-label transcript is that of the untruncated core after target construction, with
the exhaustive-search wrapper handled by domination; "3a never binds" only on the diagonal τ′ = τ (off-diagonal
counterexample recorded); in §8, unit multiples defeat pair uniqueness and the shape test is not an obstacle; the
characteristic ranges stated quantitatively. Weakest new link in its view: XK4-U's quantified implementation.

| reviewer | scope | verdict | main fixes (all applied) |
|---|---|---|---|
| GPT-6.1-sol dlog25_1 (blind, xhigh) | whole question | YES L[2/5], conditional on U's chain | c_U > 4.129a > 3a (diagonal identity), KW exponent, extension-field objections |
| Opus proof critic | §§3–8 line by line | ACCEPT WITH FIXES, no FATAL/MAJOR | √s in R′, certified roots, bad-prime valuations, grid wording, scope wording |
| Opus literature critic | §1, citations, novelty | OK WITH FIXES | 3 novelty overclaims removed; heuristic split-and-descend prior art; Lee–Venkatesan, Lido, GKZ, Enge–Gaudry rows |
| GPT-6.1-sol dlog25_2 (adversarial, xhigh) | whole note | ACCEPT WITH FIXES, no FATAL/MAJOR | Pan's near-linear root cost, starting-order discriminant, core transcript, diagonal qualifier |

None of the reviewers certifies Theorem U or its chain; Theorem DL25 stays labelled **proved from U**.

## 11. Dependencies of Theorem DL25

- Internal (unrefereed): Theorem U and its whole dependency list ([U] §10), used twice: the supply chain with N := p,
  and FACTOR_U as the shape-test subroutine; [DL] Lemmas XK0–XK5, HX, T_p, Corollary DL′ (reviewed at fixed d, all
  degree-free except XK0(iii) and XK4, re-proved here); [rec] Lemma W (adaptive sub-calls).
- Proved here (new to this project; no prior formulation of these statements found by us; methods partly classical,
  see §1 and [DL] §9b): Lemma S_p″, K-U for DL, Lemma T_p-U, Lemma Sel-DL, Lemma XK0-U, Lemma XK4-U, Proposition XK-U, Theorem DL25,
  Corollary DL25′, Lemma C25.
- Published: Blanksby–Montgomery 1971; Nguyen–Stehlé 2009 (L²); Schönhage 1982 / Pan 2002 (roots); Dixon 1982
  (p-adic lifting); Buchmann–Lenstra 1994 and Cohen GTM 138 (local maximal orders, decomposition, valuations);
  Kronecker's theorem; AKS; everything Theorem U cites.
