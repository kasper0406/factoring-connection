# cuspcoef: pre-registration (phase 30)

Agent `cuspcoef`, 2026-10-28. Written and committed BEFORE any experiment below is run. Deviations will be disclosed in
`notes/cuspcoef_results.md`.

Question. One coefficient a_f(N) of a fixed non-CM cusp form at N = pq: does it (or a few such coefficients) factor N
(R), can it be computed natively below L[1/3] (A), or is there a sharp model barrier (B)?

Working hypotheses going in (to be proved, refuted or reworded in the results):
- H-B (barrier): in the Lang–Trotter independent-trace model, any datum that is a fixed function of N and of the
  trace pairs (a_p(f_i), a_q(f_i))_{i <= K} ("determinant-free data") changes the success probability of any
  factoring algorithm by at most a factor C^K, with C = O(B^2) and B the spread constant of the trace law
  (B = O(log p loglog p) for the Gekeler law). So one coefficient is worth O(log log N) bits in the model.
- H-R (many curves): the products a_N(E_i) of K = O(n^2 log n) fixed non-CM curves, TOGETHER WITH the factorisations
  of these products (integers of about n/2 bits), factor N (heuristic: Chebotarev independence across a family).
  Without the factorisations, the congruence content is quadratic-residue bits only (tracedata Theorem SC).

Code: OPERATIONAL `src/cuspcoef/decode.py` (inputs: N, the public curve list, the list of a_N values; no PARI, no
factorisation of N, no p, q). HARNESS `src/cuspcoef/harness.py`, `src/cuspcoef/run_experiments.py` (PARI via cypari2,
nbthreads = 1). Raw data in `experiments/raw/cuspcoef_<exp>.jsonl`. Peak RSS target ≤ 2 GB; jobs run one at a time.

Instances: p, q independent uniform random primes in [X, 2X), X = 2^{n/2 - 1}, p ≠ q, seeded (`random.Random(seed)`),
seeds 1000 + i.

Curve family F (fixed, public): y^2 = x^3 + Ax + B with A ∈ [−160, 160] \ {0}, B ∈ [1, 160], 4A^3 + 27B^2 ≠ 0,
ordered by (max(|A|, B), A, B); keep the first curve for each j-invariant (removes quadratic twists, whose a_N differ
only by (d/N)); drop the 13 CM j-invariants. F_K = first K curves.

## E1. Products of K curves → factoring (exact toy end-to-end)
- n ∈ {32, 40, 48}; 20 instances per n; K ∈ {2000, 8000, 32000}.
- Harness supplies A_i = a_p(E_i)·a_q(E_i) (PARI `ellap`). Operational decoder gets (N, F_K, (A_i)).
- Decoder:
  1. Drop curve i if gcd(Δ_i, N) ≠ 1 (a proper gcd is recorded as an ACCIDENTAL split, reported separately; the
     structural decoder continues without it). Drop A_i = 0.
  2. "Half-size factoring oracle": factor |A_i| completely (toy: smallest-prime-factor sieve; |A_i| ≤ 2^{n/2+2}).
     Candidate pairs C_i = {(x, y) ∈ Z^2 : xy = A_i, |x|, |y| ≤ R}, R = ⌊2√(2X)⌋.
  3. For each odd prime ℓ ≤ L_n (L_32 = 17, L_40 = 19, L_48 = 23; the least L with ∏_{3≤ℓ≤L} ℓ > 2X), log-likelihood
     of d = p mod ℓ ∈ F_ℓ^*: Σ_i log Σ_{(x,y)∈C_i} n_ℓ(x, d)·n_ℓ(y, N/d), with n_ℓ(t, d) = ℓ(ℓ + χ(t^2 − 4d)) the
     GL2(F_ℓ) class count (Haar/Chebotarev model). Convert to a posterior on s mod ℓ (s = d + N/d), uniform prior on d.
  4. Best-first list decoding over the per-ℓ rankings (≤ 3000 combinations), CRT with s ≡ 0 (mod 2), all lifts
     s ∈ [2√N, 4X], accept if s^2 − 4N is a square and the root divides N.
- Controls (K = 32000 only, same instances):
  - C1 (null lift = the model simulator): A'_i = x'·y' with x', y' independent uniform integers in [−R, R]; same decoder.
  - C2 (congruence-only): decoder uses only A_i mod ℓ (likelihood Σ_{tt' ≡ A_i} n_ℓ(t, d) n_ℓ(t', N/d)); no factorisation.
- Predictions:
  - P1a: K = 32000 succeeds on ≥ 17/20 instances at n = 32 and n = 40, and ≥ 14/20 at n = 48.
  - P1b: success is non-decreasing in K at each n (tolerance 1 instance).
  - P1c: C1 succeeds on ≤ 1/20 at each n.
  - P1d: C2 succeeds on ≤ 2/20 at each n.
  - P1e: at K = 32000, per-ℓ top-1 accuracy of s mod ℓ ≥ 0.8 for every ℓ ≤ 13, pooled over n.
- Reported: success, per-ℓ top-1 accuracy and rank of the truth, number of combinations tried, mean |C_i|, accidental
  gcd splits (expected O(K/p) per instance at n = 32), wall time.

## E2. Can the half-size factorisation be avoided by selecting easy products? (scaling)
- n ∈ {64, 128, 192, 256}; 600 (instance, curve) samples per n, curve drawn uniformly from F_32000, instance as above.
- Acceptance event Acc(A): A ≠ 0 and A/(its B-smooth part) ∈ {1} ∪ {primes}, B = n^2 (then the complete
  factorisation, hence all candidate pairs, is known in poly(n) time). Real A = a_p a_q (PARI) vs simulator
  A' = x'y' (x', y' uniform in the Hasse ranges [−2√p, 2√p], [−2√q, 2√q]; 10^5 samples per n).
- Predictions:
  - P2a: real acceptance rate within a factor 2 of the simulator rate at each n where the simulator rate ≥ 0.02
    (binomial noise allowed: |difference| ≤ 3 SE also counts as agreement).
  - P2b: simulator rate is strictly decreasing in n and ≤ 0.02 at n = 256 (super-polynomial decay, since one of the two
    factors must be B-smooth; Dickman ρ(u) with u ≈ (n/4)·ln 2/ln B).

## E3. Spread constant of the Gekeler law (model constant for H-B)
- For 10 random primes p at each of 16, 20, 24 bits: B_emp(p) = max_{|t| < 2√p} √p·H(4p − t^2)/(2p)
  (H = Hurwitz class number, PARI `qfbhclassno`; Σ_t H(4p − t^2) = 2p is checked for every p).
- Predictions: P3a: the Kronecker–Hurwitz sum check holds for all 30 primes. P3b: max B_emp ≤ 4 over all 30 primes,
  and below the explicit bound (1/π)(2 + ln 4p)(e^γ ln ln 4p + 1) for every p.

## E4. 37a1 image checks (input to the scalar hypothesis of Theorem SC)
- For every prime ℓ ≤ 47: the set of (a_p mod ℓ, p mod ℓ) over good primes p < 2·10^5, p ≠ ℓ.
- Predictions: P4a: all ℓ(ℓ − 1) pairs (t, d), d ≠ 0, occur for every 5 ≤ ℓ ≤ 47 (with Dickson this proves surjectivity
  mod ℓ, given det surjective). P4b: for ℓ = 3 all 6 pairs occur. P4c: a_p odd ⇒ (37/p) = +1 holds for every p
  (0 violations; the Serre-curve entanglement Q(√37) ⊂ Q(E[2])), and the converse fails for some p.

Order: E4, E3, E1, E2. Code changes after this commit that affect outcomes will be disclosed.
