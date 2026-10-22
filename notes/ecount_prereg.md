# ecount pre-registration (phase 23)

Agent `ecount`, 2026-10-22. This file is committed before any experiment below is run. All experiments use exact
integer arithmetic; floats appear only for display. Code lives in `src/ecount/`, raw output in
`experiments/raw/ecount_e{1,2,3,4}.jsonl`, and the harness computes #E(F_p) and group structures with PARI (cypari2
2.2.4, added to `.venv`). The PARI values are cross-checked against naive Legendre-sum counts in the tests.

Notation:
- N = pq with p < q primes.
- E: y² = f(x) = x³ + ax + b is a FIXED curve, and E′ is its quadratic twist.
- A_r = #E(F_r) and A′_r = #E′(F_r) = 2r + 2 − A_r.
- M = #E(Z/N) := A_p·A_q.
- r_r is the number of roots of f mod r.
- H′_r = E′(F_r)[M] (the M-torsion of the twist).

## Claim under test (to be proved in notes/ecount_results.md): Theorem TW, one-curve count ⇒ factoring

**Algorithm TW(N, a, b, M)** (OPERATIONAL). Its only inputs are N, the curve and M.
- Pick x₀ uniform in Z/N and set δ = f(x₀) mod N. If 1 < gcd(δ, N) < N, output the factor.
- If the Jacobi symbol (δ/N) ≠ −1, this trial fails.
- Otherwise P = (δx₀, δ²) is a point of E^(δ): Y² = X³ + aδ²X + bδ³ over Z/N. Compute [M]P with the affine group law
  over Z/N, splitting cases by gcds; every non-invertible denominator with 1 < gcd < N is output as a factor.

**Claimed bound.** The per-trial split probability is ≥ ρ(p, q), where

  ρ(p, q) = [(A_p − 1 − r_p)(A′_q − |H′_q|) + (A′_p − |H′_p|)(A_q − 1 − r_q)] / (4pq).

So TW fails on every trial iff exp(E′(F_p)) | M and exp(E′(F_q)) | M. Call this the exceptional set Exc(E). Outside
Exc(E), ρ ≥ 1/8 − O(p^{−1/2}).

## E1. TW on random balanced semiprimes

**Instances.**
- Generator: `src/common/instances.balanced_semiprime` with seeds 0–39.
- Sizes: n ∈ {32, 64, 96, 128}.
- Curves:
  - C1 = 11a1 in short form, y² = x³ − 13392x − 1080432;
  - C2 = 37a1, y² = x³ − 1296x + 11664;
  - C3: y² = x³ + 2x + 3;
  - CM1: y² = x³ − x;
  - CM2: y² = x³ + 1.
- Instances are resampled (with seed suffix) until gcd(N, 6Δ·cond) = 1.
- Up to 64 TW trials per instance. The trial RNG is `random.Random(f"ecount-e1|{curve}|{n}|{seed}")`.

**Controls.**
- (i) Wrong M: run TW with M + 2 on the same instances. This is the null; any split there is classified ACCIDENTAL.
- (ii) Designed negative: CM1 on p ≡ q ≡ 3 (mod 4). Here a_p = a_q = 0, so E ≅ E′ in count. We use 20 instances per n
  from a seeded sampler that resamples until p ≡ q ≡ 3 (mod 4). TW must fail, and SR2 (M = (p+1)(q+1)) must split.

**Predictions.**
- P1: on non-CM curves, every instance outside Exc(E) splits within 64 trials. We expect Exc to contain none of the
  480 instances.
- P2: the pooled empirical per-trial split rate is ≥ the pooled exact ρ (harness, PARI ellgroup), and within
  ±3 binomial standard errors of ρ + bonus, where bonus = Pr[f(x₀) ≡ 0 at exactly one prime]. The algorithm may split
  slightly more often than ρ through intermediate steps. That excess is counted, not predicted.
- P3: the wrong-M control has 0 splits, except for gcd(δ, N) hits whose rate matches the bonus term.
- P4: CM1 with p ≡ q ≡ 3 (mod 4) gives 0 TW splits, and SR2 splits 100%.

## E2. Exact census of the exceptional set Exc(E) and of ρ

**Population.**
- S_b = all primes in [3·2^{b−2}, 2^b) with b ∈ {12, 14, 16, 18}.
- All unordered pairs p < q in S_b with good reduction, for curves C1, C2, C3, CM1 and CM2.

**Per-pair data.** The membership of (p, q) in Exc(E), exactly (via E′(F_r) ≅ Z/m₁ × Z/m₂ from PARI), and ρ(p, q).

**Predictions.**
- P5: for the non-CM curves, |Exc ∩ S_b²|/|S_b²| ≤ 10^{−3} for every b, and it is non-increasing in b.
- P6: the fraction of pairs with ρ < 1/8 − 4/√p is ≤ |Exc|/pairs + 10^{−3}. In other words, outside Exc the lower
  bound 1/8 essentially holds, and median ρ ≈ 1/2.
- P7: for CM1, Exc contains every pair with p ≡ q ≡ 3 (mod 4) (they are supersingular at both primes). For CM2, Exc
  contains every pair with p ≡ q ≡ 2 (mod 3). For both CM curves, the remaining pairs obey P5.
- Descriptive: the composition of Exc by type (both-supersingular / one-supersingular / ordinary-ordinary).

## E3. Dieulefait–Urroz (arXiv:1911.11004 v1) checked on toy examples

All D–U algorithms are OPERATIONAL. Their inputs are N, the curve(s), the count(s) and, for (c), the factorisation of
M from the harness (an oracle, as in their Game 1).

- (a) **Affine twist reduction (their Thm 4, affine case).** Instances: n ∈ {32, 64, 128}, 20 seeds, curves C1–C3, and
  d chosen as the least d > 1 with (d/N) = −1 (case 2), plus d with (d/p) = (d/q) = −1 (case 1; d chosen by the
  harness). Prediction P8: 100% split in both cases.
- (b) **Doubly supersingular counterexample to their Thm 1 / Thm 5 (affine).** CM1 with p ≡ q ≡ 3 (mod 4): for every
  d, E*(Z/N) = E_d*(Z/N) = N. Prediction P9: the affine counts equal N in 20/20 instances. No algorithm can factor N
  from them, because they are functions of N alone.
- (c) **Malleability theorem (their Thm 6).** The input is N together with the factorisation of M (harness oracle).
  The algorithm enumerates the divisors u of M in [√N/4, 4√N] and runs Coppersmith (Howgrave-Graham lattice,
  python-flint LLL) for a factor within 2√u + 2 of u, on n ∈ {48, 64, 96}, 20 seeds, curves C1–C3.
  - Prediction P10: 100% split.
  - Descriptive: the distribution of d(M) and of the number of divisors tried.
- (d) **Projective twist reduction.** We implement D–U's Lemma 3 / Thm 4 as written, and our repaired version.
  - The repair: in case 2 ((d/N) = −1), the reduced fraction (E + E_d)/(E − E_d) yields P/k with k | a_p, then
    Coppersmith on p ≡ −1 (mod P/k). In case 1 we use TW.
  - Prediction P11: the as-written projective step cannot be executed without (p+1)(q+1), so it is not run. The
    repaired version splits 100% (C1–C3, n ∈ {32, 64, 96}, 20 seeds).

## E4. Partial information: what M mod ℓ and root counts would reveal (exact census, HARNESS)

**Population.** All ordered pairs (p, q), p ≠ q, from S_16 with good reduction, for curves C1, C2 and C3.

**Quantities.**
- (i) The Bayes advantage for predicting the fixed-base QR bit β_ℓ = (ℓ*/p), on pairs with (N/ℓ) = 1, from
  (N mod ℓ, M mod ℓ), for ℓ ∈ {3, 5, 7}.
- (ii) The advantage on β_Δ = (Δ_f/p), on pairs with (Δ_f/N) = 1, from M mod 2. Equivalently, from the root count
  r_p·r_q of f mod N, and from #E(Z/N)[2].
- (iii) The same advantages predicted by Chebotarev for independent Frobenius classes in the mod-ℓ image (GL₂(F_ℓ),
  or the actual image for C1 at ℓ = 5, which has rational 5-torsion) and in S₃ for f.

**Predictions.**
- P12: the empirical advantages match the Chebotarev predictions within 0.01 (finite-population effects at b = 16 are
  O(10^{−2})).
- P13: the advantages are nonzero for every ℓ with surjective image. They are 0 for C1 at ℓ = 5, where M ≡ 0 always.
- P14: M mod 2 has positive advantage on β_Δ, and the advantage equals Pr[M odd | (Δ/N) = 1]/2, exactly as derived
  in the notes.

## Not pre-registered
- Unit tests of the group law, TW, D–U and Coppersmith.
- Any post-hoc analysis, which will be labelled as such.
