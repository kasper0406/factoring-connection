# tracedata pre-registration (phase 25)

Agent `tracedata`, 2026-10-23. Committed before any experiment below is run. Exact integer arithmetic throughout;
floats only for likelihoods/display. Code: `src/tracedata/`; raw output `experiments/raw/tracedata_e{1..5}.jsonl`.
Local data (a_p, #E(F_p), group structure) come from PARI (cypari2) in the HARNESS only. OPERATIONAL code
(`src/tracedata/reductions.py`, `src/tracedata/decode.py`) receives N and the stated data only, never p, q, and never
calls factorisation; an AST test enforces this.

Notation: N = pq, p < q; E fixed; a_r = r + 1 − #E(F_r); a_N = a_p a_q; M = #E(Z/N) = (p+1−a_p)(q+1−a_q);
s = p + q; S = q a_p² + p a_q².

## E1. Exact reductions (checks of statements to be proved in notes/tracedata_results.md)
Instances: `balanced_semiprime(n, seed)`, seeds 0–29, n ∈ {32, 64, 128}, skipping N with gcd(N, 6Δ_E) > 1.
Curves: C1 = 11a1 (E(Q)_tors = Z/5), C2 = 37a1, C3: y² = x³ + 2x + 3 (Z/2), CM1: y² = x³ − x (32a2, CM by Z[i]),
T7 = 26b1 (Z/7).
- (a) BC2: input (N, a_N, a_{N²}); roots of Z² − SZ + N a_N² with S = a_N² + N − a_{N²}; gcd with N.
  **P1a:** success on every instance with (a_p, a_q) ≠ (0, 0); failure exactly when a_p = a_q = 0.
- (b) BC3: input (N, a_N, a_{N³}). **P1b:** success iff a_N ≠ 0.
- (c) MIX: input (N, {a_p, a_q} unordered, M). **P1c:** success 100% (deterministic quadratic root-finding).
- (d) CM2: CM1 only, N with p ≡ q ≡ 1 (mod 4) (seeds drawn until 30 such N per n): input
  (N, a_N(f), a_N(g)) with f = 32a2's newform and g = η(4z)⁶ (weight 3, a_p(g) = a_p(f)² − 2p). **P1d:** 100%.
- (e) TOR: C1 (m = 5), T7 (m = 7), C3 (m = 2): s mod m recovered from (N, a_N) via a_N ≡ (p+1)(q+1) (mod m).
  **P1e:** 100% agreement with the true s mod m.
- (f) CM-QR: CM1, N ≡ 1 (mod 4): [a_N ≠ 0] = [(−1/p) = 1]. **P1f:** 100% agreement.

## E2. Haar-model constants (exact enumeration over GL₂(F_ℓ), ℓ prime 3–47; no primes involved)
Model: g_p, g_q independent Haar on GL₂(F_ℓ); p mod ℓ = det g_p; t = tr g.
- I_ord(ℓ) := I(det g_p ; (t_p, t_q) | det g_p det g_q) in bits.
  **P2a:** ℓ²·I_ord(ℓ) ∈ [1.0, 1.9] for all ℓ ≥ 11, and |ℓ²·I_ord(ℓ) − 1/ln 2| ≤ 0.15 for ℓ ≥ 23.
  **P2b:** Σ_{3 ≤ ℓ ≤ 47} I_ord(ℓ) ≤ 1 bit.
- I_prod(ℓ) := I({det g_p, det g_q} ; t_p t_q | det g_p det g_q) (information about s mod ℓ).
  **P2c:** ℓ³·I_prod(ℓ) ∈ [0.1, 20] for all ℓ ≥ 11 (scaling check; exploratory constant).
- PT constants: for F ∈ {M mod ℓ (ℓ = 3, 5, 7), a_N mod ℓ (ℓ = 3, 5, 7), r_N(ψ₃) (model PGL₂(F₃) ≅ S₄ on lines),
  r_N(f) (S₃ cubic)}: compute 1 − max_b Pr[F = b] (unconditioned model). **P2d:** each ≥ 0.2.

## E3. Multi-curve Galois decoding (HEURISTIC claim: trace numbers of many fixed curves determine p mod ℓ)
Curve list (fixed, N-independent): y² = x³ + Ax + B, 1 ≤ A, B ≤ 200, ordered by (max(A, B), A, B); index i = 0, 1, ….
All primes used are ≥ 2^26 > |4A³ + 27B²|, so every curve has good reduction (no accidental gcds possible).
Decoder (OPERATIONAL): maximum likelihood over d ∈ F_ℓ^×, with Haar likelihood P(t | d) ∝ ℓ + χ(t² − 4d), χ the
Legendre symbol (χ(0) = 0); for products, P(u | d) = Σ_{t t' ≡ u} P(t | d) P(t' | N/d).
Haar simulation: same K, ℓ, 2000 trials per cell, traces drawn from the model.
- **E3a (ordered).** 200 random 30-bit primes (seed 1). K ∈ {250, 1000, 4000, 10000}; ℓ ∈ {3, 5, 7, 11, 13, 17,
  19, 23, 29, 31}. Input: (a_p(E_i) mod ℓ)_{i < K}. Target p mod ℓ.
  **P3a:** |real accuracy − Haar-sim accuracy| ≤ 0.10 in ≥ 38 of the 40 cells.
  **P3b (null control):** traces of a different random prime p′: accuracy within 3 SE of 1/(ℓ−1) in ≥ 38/40 cells.
- **E3b (products).** 100 N = pq, p, q random 30-bit (seed 2). K ∈ {1000, 10000, 40000}; ℓ ∈ {3, 5, 7, 11, 13}.
  Input (N, a_N(E_i) mod ℓ). Target s mod ℓ. **P3c:** |real − Haar-sim| ≤ 0.12 in ≥ 13 of 15 cells.
  **P3d (null):** a_N from an unrelated N′: accuracy within 3 SE of the model's chance level in ≥ 13/15 cells.
- **E3c (end-to-end, ordered).** 30 N = pq, p, q random 28-bit (seed 3), K = 40000, ℓ ∈ {3, …, 29} (odd primes),
  list-decoding: top-2 residues per ℓ, all 2⁹ CRT combinations with p odd tested by divisibility of N.
  Input: N and (a_p(E_i))_{i<K} (one prime's traces). **P3e:** success ≥ 24/30.
- **E3d (end-to-end, products only).** 30 N as in E3c (seed 4), K = 40000, input N and (a_N(E_i))_{i<K};
  top-3 candidates for s mod ℓ per ℓ ∈ {3, …, 29}, s ≡ 0 (mod 2), all CRT combinations tested by solving
  x² − sx + N = 0 over Z. **P3f:** success rate within ±0.25 of the Haar-simulated end-to-end rate (200 simulated N).
These are toy sizes; any success is a check of the model, not evidence of an efficient factoring algorithm (the
traces are given by the harness).

## E4. Fixed-ℓ counts on a prime census
Primes: all primes in [2^15, 2^16) with the top two bits set (as ecount E4); ordered pairs p ≠ q. Curves: C2 = 37a1
and C5: y² = x³ − x + 1 (conductor 92). If the harness finds either curve non-surjective mod 3 (a Frobenius-class
census missing a GL₂(F₃) class), that curve is reported but excluded from P4a.
For ℓ = 3: Bayes advantage for predicting (p/3) on (N/3) = 1 from (a) (N mod 3, a_N mod 3), (b) r_N(ψ₃).
**P4a:** census advantage within ±0.03 of the Haar/S₄ model value, for both items and both curves.
**P4b:** r_N(ψ₃) ∈ {1, 16} ⇒ p ≡ 1 (mod 3), and r_p(ψ₃) = 2 ⇒ p ≡ 2 (mod 3): 0 exceptions.

## E5. X_E descent
Instances: the ecount census X_E member (11a1, p = 54443, q = 62819) and cx1's example (y² = x³ + 27331x + 2,
N = 27331). Algorithm: TW* plus descent over the primes ℓ ≤ 1000 dividing M (scalars M/ℓ^j), random x₀, ≤ 20 trials.
**P5:** both split.

## Labels
Results will be labelled PROVED / HEURISTIC / experimentally supported / REFUTED (prediction) / OPEN. Deviations from
this file will be disclosed.
