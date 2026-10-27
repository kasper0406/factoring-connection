# nativeinv pre-registration (phase 29)

Agent `nativeinv`, 2026-10-27. Frozen before any code in `src/nativeinv/` runs.

The phase-29 task is mostly an inventory and an analysis. The only computations are small correctness checks of
identities and of one operational decoder. None of them is evidence about asymptotic cost (AGENT_RULES rule 8).

## E1. Exact identities by brute force (HARNESS; N = pq, small odd primes)
Each check uses N ∈ {15, 21, 33, 35} where feasible (group orders ≤ ~10⁵ elements; ≤ 2 GB RSS).
Prediction: every identity holds exactly for every N tested.
1. k(SL₂(Z/N)) (conjugacy classes, by orbit enumeration) = (p + 4)(q + 4) = Σ_{t ∈ Z/N} 3^{ω(gcd(t² − 4, N))}.
2. k(GL₂(Z/N)) = (p² − 1)(q² − 1) = Σ_{t ∈ Z/N, d ∈ (Z/N)^*} 2^{ω(gcd(t² − 4d, N))}  (N ∈ {15, 21} only).
3. k(Aff(Z/N)) = N, and Σ of squared irreducible degrees over the predicted degree multiset
   {1^{φ(N)}, (p−1)^{q−1}, (q−1)^{p−1}, φ(N)^1} equals |Aff(Z/N)| = Nφ(N); the number of predicted irreducibles is N.
4. k(Heis(Z/N)) (upper unitriangular 3×3 over Z/N) = (p² + p − 1)(q² + q − 1).
5. Pillai P(N) = Σ_g gcd(g, N) = #{(x, y) ∈ (Z/N)²: xy ≡ 0} = (2p − 1)(2q − 1).
6. Number of subgroups of (Z/N)² = (p + 3)(q + 3).

## E2. Trichotomy decoder (OPERATIONAL decoder; harness supplies f(N))
Quasi-polynomial local data F(r) = Σ_i c_i(r mod m) r^i, f(pq) = F(p)F(q). The operational decoder sees only N, the
value f(N) and the public description of F. It enumerates the residue pairs (u, v) with uv ≡ N (mod m), forms the
integer polynomial in p, finds its integer roots exactly, and tests divisibility.
- FACT-tier local data (two distinct degrees at the true residues): r − 1, r + 1, 2r − 1, r + 4, r² − 1, r² + r − 1,
  r + 3, r − χ₄(r), r³ − 1, (r − 1)² (φ, ψ/σ, Pillai, k(SL₂), k(GL₂)=J₂, k(Heis), subgroups of Z_r², unit circle,
  #K₅, a non-squarefree-type example).
- FIBER-tier local data (monomial of equal degree at both residues): 1 + χ₄(r) (r₂-type), 1 + (−3/r) (ν₃-type),
  r·χ₄(r) + r (= r(1 + χ₄(r))).
- AGG-tier: χ₄(r)·r² (computable as χ₄(N)N²).
Instances: `src/common/instances.balanced_semiprime(n, seed)`, n ∈ {32, 64, 128}, seeds 0..39 each (120 instances).
Predictions:
- P1: FACT tier: the decoder returns {p, q} on 120/120 instances for every FACT-tier F.
- P2: FIBER tier: whenever the true residue pair is monomial-equal, the decoder returns no split; any split is classified
  ACCIDENTAL. Predicted: 0 splits. Where a FIBER-tier F happens to have a mixed-degree pair (e.g. r(1 + χ₄(r)) vanishes
  identically on r ≡ 3 mod 4, so the pair is "zero"), the decoder reports the tier it found; no split predicted.
- P3: AGG tier: f(N) equals χ₄(N)N² on 120/120 (no factor needed).

Raw output: `experiments/raw/nativeinv_e2.jsonl`. Command: `.venv/bin/python -m src.nativeinv.run_experiments`.
Any deviation will be recorded in the results note as a disclosed deviation.
