# nativemod: pre-registration (phase 29)

Agent `nativemod`, 2026-10-27. Written and committed before any experiment below is run. Scratch prototypes run so far
(not experiments, not recorded): PARI `mfheckemat` on M_2(Γ₀(11)) at n = 35 (charpoly (x − 48)(x + 2)); a search for
the genus partner of the level-11 quaternary form (x² + xy + 3y²) ⊕ (z² + zw + 3w²), which found
L2 = [[4,−1,−1,2],[−1,4,2,−1],[−1,2,4,0],[2,−1,0,4]] with θ_{L1} − θ_{L2} = 4·f_{11a1} through q³⁰.
The codex blind round `codex/nativemod_1.out` has been read.

Notation: N = pq, p < q odd primes; instances from `src/common/instances.balanced_semiprime(n, seed)`.
All experiments are correctness checks or a null leak search; none is a factoring algorithm.

## E1 — exact identities (harness; PARI via cypari2; all predictions: 0 mismatches)

- **E1a (Eichler–Selberg Hurwitz moments).** C_k(N) = Σ_{t² < 4N} P_k(t, N) H(4N − t²),
  P_k = [z^{k−2}](1 − tz + Nz²)^{−1}. Population: every odd semiprime N = pq ≤ 4000 with p ≥ 3 (both coprime to 6
  not required). Predictions:
  - k ∈ {4, 6, 8, 10, 14}: C_k(N) = −2(1 + p^{k−1}); k = 2: C_2(N) = 2(N + q).
  - k ∈ {12, 16, 18, 20, 22, 26}: C_k(N) = −2 Tr(T_N | S_k(SL₂Z)) − 2(1 + p^{k−1}), Tr from PARI `mftraceform`
    / `mfheckemat`. (Restricted to N ≤ 1500 for k ≥ 12 if PARI is slow; disclosed if so.)
- **E1b (fixed-level Hecke matrices).** M ∈ {11, 37, 43, 53}; 25 semiprimes N ≤ 10⁴ coprime to M per level
  (seeded). On the full space M_2(Γ₀(M)) (PARI `mfinit([M,2],4)`):
  (i) T_N = T_p T_q; (ii) σ(N) is an eigenvalue of T_N; (iii) T_N² − T_{N²} − N·I = p·T_{q²} + q·T_{p²}.
- **E1c (theta series at level 11 and class-number-one quaternary lattices).** For all n ≤ 3000:
  (i) a_{11a1}(n) = (r_{L1}(n) − r_{L2}(n))/4;
  (ii) for 11 ∤ n: (r_{L1}(n)/|O(L1)| + r_{L2}(n)/|O(L2)|) / (1/|O(L1)| + 1/|O(L2)|) = (12/5)·σ(n);
  (iii) odd n: r_{Z⁴}(n) = 8σ(n); r_{D4}(n) = 24σ(n); for 3 ∤ n: r_{A2⊕A2}(n) = 12σ(n).
  If (ii)'s constant or (iii)'s A2⊕A2/D4 constants differ from the stated ones, that is a prediction failure of
  the constant, recorded as such; the linear relation itself is what is tested.
- **E1d (Merel's Heilbronn set, descriptive only, no prediction).** #X_n for X_n = {(a b; c d) : ad − bc = n,
  a > b ≥ 0, d > c ≥ 0}, n ≤ 600, against (12 log 2/π²)·σ(n)·log n.

## E2 — operational reductions (toy correctness; prediction: 100% success, 0 accidental)

Operational code sees only N and the stated datum; the datum is produced by the harness from (p, q).
n ∈ {64, 128, 256} (E2b: n ∈ {64, 128}), 40 instances per cell, seeds 0..39.
- **E2a (Prop EIS-k).** Datum: a(N) for E_k(χ, χ′), a(n) = Σ_{d|n} χ(n/d)χ′(d)d^{k−1}, for
  (k, χ, χ′) ∈ {(2,1,1), (4,1,1), (12,1,1), (3,1,χ₋₄), (3,χ₋₄,1), (2,χ₋₄,χ₋₄)}. Prediction: 240/240 per family
  factored deterministically by the quadratic-in-p^{k−1} method over the character patterns consistent with N.
- **E2b (CRT trace).** Datum: t_N = CRT(a_p mod p, a_q mod q) for E ∈ {37a1, 11a1, y² = x³ − x}. Decoder:
  univariate Howgrave-Graham/Coppersmith for the small root x = a_q of x − t ≡ 0 (mod q), q ≥ √N. Prediction:
  80/80 per curve. Null control: t′ uniform in Z/N, 40 per cell: 0 factors.
- **E2c (Frobenius class pair).** Datum: for primes ℓ ≤ L (∏ℓ > 2√N + 2N), the unordered pair {p mod ℓ, q mod ℓ}
  (what det of the two Frobenius classes of ρ̄_{f,ℓ} gives in weight 2). Prediction: 240/240 (unit test only).

## E3 — null leak search in the open PRIM fragment (extension components + Jacobi of deterministic gate values)

Algebra A_ℓ = (Z/N)[x]/(ψ_ℓ/ℓ) (monic division polynomial, Schoof's algebra; its components at p are F_{p^f}),
E = 37a1 (y² + y = x³ − x, full mod-ℓ image), ℓ ∈ {3, 5}. Population: balanced semiprimes, n ∈ {32, 48},
20 000 per (n, ℓ), seeds 0..19999, excluding gcd(N, 6·37·ℓ·disc ψ_ℓ) > 1 (trial division, counted).
Operational statistics (N only): y = x^N mod (ψ_ℓ, N);
- S1 = (Tr_{A/(Z/N)}(y) | N), S2 = (Tr(x·y) | N), S3 = (Tr(y²) | N) (Jacobi symbols; 0 counts as a split check),
- S4 = [ψ_ℓ(y) = 0 in A] (element test),
- every computed Z/N value is gcd-checked; splits classified STRUCTURAL (explained: value ≡ 0 at exactly one prime
  through a stated mechanism) or ACCIDENTAL.
Targets (harness): T1 = a_N mod ℓ; T2 = [a_N ≡ 0 mod ℓ]; T3 = ε = (−1/p) on N ≡ 1 (mod 4); T4 = (p/ℓ) on
(N/ℓ) = 1.
Test: for each (statistic ∈ {S1, S2, S3}, target), a stratified association z-score (Cochran–Mantel–Haenszel over
strata N mod 4ℓ, binary statistic vs binary target; T1 as the indicator [a_N ≡ c] for each c, max |z| reported with
the multiplicity count). Positive control: J_CRT = (c | N) with c built by CRT so that (c|p) = target (harness), for
T2 at n = 32.
Predictions:
- P3a: every |z| < 4 for S1–S3 (null).
- P3b: S4 = 1 in ≤ 0.1% of instances.
- P3c: ACCIDENTAL splits ≤ 3 per cell; STRUCTURAL splits: none expected.
- P3d: positive control |z| > 10.
A violation of P3a would be reported as a candidate leak and re-tested on fresh seeds 20000..39999 before any claim.

## Resources
Each job alone; peak RSS ≤ 2 GB (checked with /usr/bin/time -l); no background processes left.
