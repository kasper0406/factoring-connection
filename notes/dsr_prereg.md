# dsr pre-registration (phase 23, 2026-10-22)

Frozen before any coverage or Fourier statistic was computed. The only computation done before this file was a
density calibration (`delta_bal` at a ∈ {0.5, …, 0.8}, 20–24 bits, seed 0) to choose B so that δ is between 0.03 and
0.15; it computed no statistic used below.

**Purpose.** Toy-scale calibration only. Nothing here can establish or refute Hypothesis DSR, which is a statement
about all N at L-scale. The experiments measure (E1) whether smooth sets have the small Fourier-algebra norm that
Bach-type certificates need (zroute Prop FC(c)), and (E2) how explicit deterministic multiplier lists, including a
deterministic "walk", compare in the WORST CASE over targets with random lists and with the independence prediction.

**Disclosed deviation from AGENT_RULES rule 1 (exact arithmetic).** E1's quantity Σ|Ŝ(χ)| is irrational; it is
computed with float64 FFTs, with Parseval's identity checked to relative error ≤ 1e-9 in every run (a run failing the
check is reported as failed). E2 is exact integer arithmetic throughout.

## Common definitions

- L_M[1/2, a] = exp(a·sqrt(ln M · ln ln M)); B = ⌊L_M[1/2, a]⌋; P = primes ≤ B; G = (Z/M)^×, φ = |G|.
- **S_can(a)** = {x ∈ G : the least positive residue of x is B-smooth}. Primary set, a = 0.8.
- **S_bal(a)** = {x ∈ G : x ≡ u·v^{-1} with u, v B-smooth, u ≠ 0, v > 0, u² + v² < M}. Secondary set, a = 0.5. (For
  such x the vector (u, v) is the unique shortest vector of {(u, v): u ≡ xv}, up to sign, so this is the "Gauss split
  smooth and short" set.)
- δ = |S|/φ (measured, not predicted).
- Instances: `src/common/instances.py` balanced_semiprime(bits, seed), prime_instance(bits, seed).

## E1: Fourier-algebra norm of smooth sets (approach 2; Prop FC(c))

- Moduli: primes at 16, 18, 20, 22 bits (seeds 0, 1, 2); balanced semiprimes at 20, 22 bits (seeds 0, 1, 2). 18 moduli.
- Sets: S_can(0.8), S_bal(0.5). Controls: 3 uniformly random subsets of G of the same size (seed = 1000·bits + seed + k).
- Statistics: A′(S) = φ^{-1} Σ_{χ≠χ₀} |Ŝ(χ)|; β(S) = max_{χ≠χ₀} |Ŝ(χ)|; r_A = A′(S)/mean A′(random);
  r_∞ = β(S)/mean β(random). Harness uses the factorisation (discrete-log tables); labelled diagnostic.
- **P1a.** r_A ∈ [0.7, 1.3] for every modulus and both sets (smooth sets are Fourier-generic in L¹; no small-norm,
  coset-like structure).
- **P1b.** Over the four prime sizes, the least-squares slope of ln A′(S_can) against ln |S_can| is in [0.4, 0.6]
  (A′ grows like |S|^{1/2}, as for random sets).
- **P1c (directional, low confidence).** r_∞ > 1 for S_can in ≥ 75% of moduli (small-order characters bias smooth
  sets more than random sets).

## E2: worst-case coverage of all targets by explicit multiplier lists (approach 1)

- Moduli: balanced semiprimes at 20, 22, 24 bits, seeds 0–9 (30 instances).
- A multiplier list W = (w₁, w₂, …) covers target a ∈ G at step s if a·w_s ∈ S. m*(N, W) := least m such that every
  a ∈ G is covered by some w_s, s ≤ m (equivalently S·{w₁..w_m}^{-1} = G). Cap M_cap = ⌈50·ln φ/δ⌉; censored values
  are recorded as such.
- Lists (w must be B-smooth; t := ⌈ln N/ln B⌉ + 1):
  - **R_G** (control, not a legal DSR output): i.i.d. uniform elements of G.
  - **R_walk** (control; the JMV randomized algorithm): product of t i.i.d. uniform primes from P.
  - **D_box** (deterministic): the B-smooth integers 1, 2, 3, … in increasing order.
  - **D_slide** (deterministic "walk"): labels ℓ_i = ⌊iφ_g⌋ mod |P| with φ_g = (1+√5)/2 computed exactly as
    (i + isqrt(5i²)) // 2; w_s = ∏_{i=s}^{s+t−1} p_{ℓ_i}, s = 1, 2, … (a sliding window along one explicit walk).
- Set S_can(0.8) for all four lists; S_bal(0.5) for R_G and D_box only.
- Secondary statistic (the actual DSR targets): for each q ∈ P, the first s with q^H·w_s ∈ S (H = 8⌈log₂N⌉), and
  T* = max_q. Also **D_orbit**: first j ≥ 0 with q^{H+j} ∈ S (legal per the orbit reformulation in dsr_results).
- Prediction for i.i.d. covering: Pr[a uncovered after m] = (1−δ)^m exactly for R_G, so
  E m* ≈ (ln φ + 0.577)/(−ln(1−δ)). Define the normalised value ν := m*·(−ln(1−δ))/ln φ.
- **P2a (calibration).** ν(R_G) ∈ [0.8, 1.5] in ≥ 90% of (instance, set) pairs.
- **P2b.** Median over instances of m*(R_walk)/m*(R_G) ∈ [0.7, 1.5].
- **P2c.** D_box and D_slide reach full coverage within M_cap in every instance.
- **P2d.** Median of m*(D_box)/m*(R_G) ≥ 1.5 for S_can (structured multipliers are worse in the worst case; zroute E4).
- **P2e.** Median of m*(D_slide)/m*(R_G) ∈ [0.7, 1.5] (null: an explicit pseudo-random walk behaves like a random one).
- **P2f (scaling).** The median ν(D_box) for S_can at 24 bits is ≤ 2 × its median at 20 bits.
- Descriptive only: T* per list; D_orbit's max_q first hit; uncovered counts at m = ⌈ln φ/δ⌉.

## Limits stated in advance

- Sizes are 2^16–2^24; L-scale behaviour cannot be inferred. A pass of P2c/P2f is consistent with DSR's heuristic,
  not evidence for a proof; a failure would show a toy-scale worst-case defect of that list only.
- Memory: ≤ 2 GB RSS; one job at a time; `memory_pressure` checked before each run.
