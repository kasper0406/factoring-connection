# Track B + D candidate cards (agent sketchBD, 2026-10-17)

Template: BRIEF §4. Proofs and data: `notes/trackBD_analysis.md` (section numbers below refer to it).
No novelty is claimed for anything here.

---

## BD-1 — Single quotient sketch, coefficient / content gcd
- **ID and one-sentence mechanism:** compute S = F_N mod (f, N) for small public f and output gcd(N, coefficient) / gcd(N, content).
- **Target / family:** any composite N; analysed for balanced semiprimes, prime powers, multi-prime, Carmichael.
- **Representation:** r coefficients in R_N, O(n) multiplications in R_N[X]/(f) by repeated squaring.
- **Extraction:** coefficient gcds; content gcd.
- **Why a factor appears:** content gcd ∋ p ⇔ every root α of f mod p satisfies (1+α)^q = 1+α^q (Thm B1, Cor B2.1); other single-coefficient zeros are accidental (rate = uniform-residue null) or Lucas-support gaps (BD-6).
- **Known lemmas:** Frobenius in characteristic p; CRT.
- **New lemma required:** a family of f making "all roots of f mod p are roots of F_q" likely. None is possible for fixed f: each root is an independent Θ(1/p^d) event (B4/B6).
- **Correctness:** proved criterion (B1). Prime powers: always succeed (F_1 = 0 ⇒ p | every coefficient).
- **Cost:** O(n·M(r)) per sketch; success Θ(r/p) per sketch at best (null level). Exponential.
- **Randomness:** weights c random; per-instance rate is exactly computable (census).
- **Closest prior method:** AKS identity test (detects compositeness, no factor); Fermat-test gcds.
- **Smallest experiment:** done (§6.1, 6.2, 6.7).
- **Refutation:** observed coefficient-hit rates equal the null r(1/p+1/q) at 24–40 bits; structural hits only via public roots (W), rare FROB, sporadic orbits.
- **Status:** refuted for model M (§5).
- **Repair:** use the resultant / Euclid readout (BD-2), which sees single components.

## BD-2 — Resultant / Euclid / inversion failure of the sketch
- **Mechanism:** run Euclid on (f, S) over Z/N with gcd-checked leading coefficients (incl. the final constant); equivalently Res(f, S) and the inverse of S in R_N[X]/(f).
- **Why:** a factor is returned iff the remainder-degree sequences differ; forced when deg gcd(f, F_q^p) mod p ≠ deg gcd(f, F_p^q) mod q, otherwise accidental (B3). Res(f,S) ≡ Res(f,F_q) (mod p).
- **New lemma required:** a public f whose roots mod p hit roots of F_q with probability ≫ deg(f)/p.
- **Correctness:** proved; implementation bug (missing constant-remainder gcd) found and fixed, regression-tested.
- **Cost:** O(r) Euclid steps; success per sketch = (# non-public F_p-roots of F_q among roots of f)/p ≈ 1/p·E|A_p|.
- **Closest prior:** AKS + gcd; Pollard p−1 / Bach–Shallit (group-order coincidence) for the FROB part.
- **Refutation:** §6.7: SPLIT counts match the A_p-based prediction at 16–32 bits; per-sketch rate ~1e-5 at 32 bits.
- **Status:** refuted for model M.
- **Repair:** none within sketches of F_N; the exponent would have to change (Pollard-type), which is outside the track.

## BD-3 — Weighted sketch families X^r − c and rank-split linear algebra (M1, M2)
- **Mechanism:** multiplication-by-sketch matrix / Hankel matrix of projected coefficients; column-gcd elimination finds a factor iff pivot profiles differ (Prop B3).
- **Representation:** r×r matrices over R_N.
- **Why:** rank_p = r − deg gcd(f, F_q^p) (Cor B2.3).
- **New lemma required:** a distribution over public (r, c) with P(rank gap) ≥ 1/poly(n).
- **Correctness / data:** exact rank distributions over all weights (§6.5): P_gap = XOR of (#weights x^r, x ∈ A_l)/l, **independent of r ∈ [2,16]**, zero for 65–86 % of balanced instances, ≈ 2.3–2.6/p on average. Codex's r = 2 zeros extend to all r ≤ 16 for the same N.
- **Cost:** O(r³) per probe; expected Θ(p) probes.
- **Closest prior:** codex card B1 (degree-two rank gaps); AKS.
- **Status:** refuted for model M; "vary the degree" repair refuted for the binomial family (B6).

## BD-4 — Berlekamp-type matrices Q_N (g ↦ g(X^N) mod f) and Q_N − I
- **Mechanism:** rank / kernel of the composition operator by X^N.
- **Why:** rank Q_N mod p = #{distinct α^q} over roots α; split needs ord(α_i/α_j) = q (q | p^D − 1); Q_N − I kernel jumps need α^{pq−1} = 1 (Fermat) or point coincidences.
- **Correctness:** rank formula proved (Vandermonde at the points α^N); public-order roots (Φ_m, X^r − c with c of small order) give the same integer matrix on both sides — never split.
- **Data:** 0 splits at ≥ 28 bits in 19 800 matrices each (§6.4).
- **Status:** refuted for model M.

## BD-5 — Krylov / BM over the weight c (sequence c ↦ F_N(c)) and the finite-difference window
- **Mechanism:** BM over Z/N on s_c = F_N(c), c = 0..T−1 with gcd on discrepancies; equivalently Δ^K F_N(0) = (K+1)!·S(N,K+1).
- **Why:** mod p, s_c is a polynomial of degree e_p − 1, e_p = q mod (p−1); mod q of degree p − 1 (Prop D3). LC_p = e_p (for e_p ≥ 3); for the balanced generator e_p = q − p + 1, so this is a weak Fermat-gap method. Window: (K+1)!·S(N,K+1) (K ≥ 1) ≡ 0 mod p for all K ≥ e_p; proper nonunit for e_p ≤ K ≤ p−1 up to accidental q-side zeros.
- **New lemma required:** evaluate Δ^K F_N(0) = Σ_{i≤K}(−1)^{K−i}C(K,i)[(1+i)^N − i^N] mod N for one K ≈ √N in o(K) (and better than Õ(√K) to beat Pollard–Strassen).
- **Correctness:** LC and window proved; first structural nonunit observed at step 2e_p (58) or 2e_p − 1 (1); 31/90 runs hit an accidental nonunit first (§6.6).
- **Cost:** Θ(e_p) = Θ(q − p) terms (Fermat's method needs only ≈ (q−p)²/(8√N) steps). Exponential for typical balanced N.
- **Closest prior:** brief §2.4 prefix predicate / factorial unit test (Track C), BGS.
- **Status:** mechanism proved, cost exponential, missing primitive; low priority (sibling of the factorial window).
- **Other sequences:** c ↦ [X^1](F_N mod X²−c) has LC = p on both sides (no structure; 28/30 hits accidental at ~p/2 steps); m ↦ [X^0](F_N mod X^m − 1) has maximal LC (all hits accidental).

## BD-6 — Lucas-support gaps (SUPPORT class)
- **Mechanism:** a coefficient of F_N mod (X^r − c) vanishes mod p because no exponent in the Lucas support of F_N mod p = F_q(X^p) folds onto it.
- **Why:** for q = p + δ < 2p the support is p·([1,δ] ∪ [p, p+δ−1]); it misses residues mod r when 2δ ≲ r.
- **Data:** every SUPPORT instance at ≥ 24 bits has q − p ≤ 162 (24 instances, §6.7).
- **Closest prior:** Fermat's difference-of-squares method (finds p in O(1) steps when q − p ≤ N^{1/4}): strictly stronger.
- **Status:** proved mechanism; dominated by Fermat; refuted as a general method.

## BD-7 — Lifted derivative and divided Hasse derivatives (Track D)
- **Mechanism:** compute (1+X)^N mod (N², f), form F_N'/N, D^[k]F_N/N, quotient folds; take gcds / ranks.
- **Why:** these quantities (k below the least prime; r, c units) are unit multiples of G_k = (1+X)^{N−k} − X^{N−k} or are determined by S (Prop D1). G_k(α) = 0 ⇔ (1+1/α)^{N−k} = 1: a Fermat-type order condition.
- **Data:** scalar G_1 split probability (g−1)(p+q−2(g−1))/(pq) confirmed (codex claim (c)); non-public root count g − 2, mean τ(p−1) − 2; nothing non-public when g = 2.
- **Cost:** O(n) per evaluation; success Θ(τ(p−1)/p). Exponential.
- **Closest prior:** Fermat-witness gcd(a^{N−1} − 1, N); Miller–Rabin factoring of Carmichael numbers; Grantham's Frobenius test (extension-ring Fermat).
- **Status:** refuted for model M; the lift carries no information beyond G_k and S (proved).

## BD-8 — Annihilators of projected coefficient sequences in a fixed quotient (Wiedemann/M2)
- **Mechanism:** BM on u_k = ℓ(X^k S mod f): annihilator = f / gcd(f, S) generically.
- **Why / status:** same vanishing-pattern criterion as BD-2/BD-3, plus accidental ℓ-degeneracy ≈ r/p. Refuted for model M.

## Most promising (but negative-leaning) lead from this track
None of BD-1…BD-8 survives. The precise obstruction: in model M, success ⇔ a coincidence in the power map
x ↦ x^{q mod (p^d−1)} at public algebraic points (proved); empirically such coincidences occur with
probability ≈ 1/p^d per condition (conjectured), except Frobenius coincidences, which are Pollard-type order conditions. A revision must change the
operation (the exponent or the group), not the readout.
