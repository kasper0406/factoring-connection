# sicheck round 2 (ADVERSARIAL): attack the validation verdict on "Prop SI"

You are an adversarial referee (computational number theory, smooth numbers, index calculus). Your job is to find
errors, overclaims, missing hypotheses and weak experimental inferences. Be concrete: for each finding give severity
(FATAL / MAJOR / MINOR / NIT), the exact location, why it is wrong, and the fix. Label claims PROVED / KNOWN (cite) /
HEURISTIC / OPEN / UNVERIFIED. You may run small exact Python checks (≤ 10 min, ≤ 1 GB) in the current directory.

Read these files (read-only; the repository root is /Volumes/git/factoring-connection):
- notes/sicheck_results.md (the document under review; §4 end-to-end toy may still be a placeholder)
- notes/sicheck_prereg.md (registered predictions and Addendum A1)
- src/sicheck/core.py, model.py, posthoc.py, score.py, run_experiment.py, diag.py; tests/test_sicheck.py
- experiments/raw/sicheck_E1_score.txt, sicheck_E1_posthoc.txt, sicheck_E3_score.txt, sicheck_E3b_score.txt
- For context only: notes/suffmap_results.md §1.4 (Prop SI as stated), notes/nfsfruit_results.md §3–4,
  notes/abcsunit_results.md §3.2.

Attack in particular:
1. **Lemma I, Lemma CP, Lemma D** (§2.1). Is Lemma CP's proof right (saddle point, the bound y^δ ≈ u log u, the
   split at p ≤ x², the count of large prime factors)? Is the local-factor model stated honestly, and is there a
   correct citation for Ψ_q(x, y) ≈ ∏_{p|q}(1 − p^{−α})Ψ(x, y) in a range covering this use? Is Lemma D's
   P(ℓ | u) ≈ ℓ^{−α} for the *chosen* u justified for each constructor? Does the "fresh rate ≥ 0.6·ρ" survive?
2. **The E1 inference.** The pre-registered P1/P4 failed; the note explains the failures post hoc by inheritance and
   coprimality, and calls the refined H_S "experimentally supported". Is the coprime null p_cop(v; u) a fair,
   non-circular null (it conditions on u, which the constructor chose)? Are the trend claims (deficit shrinking with
   size) confounded by y, K and h varying across cells? Is "10–69% of the naive H_S rate" computed correctly (PH3
   divides fresh relations by Σ p_null over *all* outputs)?
3. **Prop B and Cor. B** (§3.2–3.3). Check the proof, including ℓ = 2, the sign of v, the case gcd(V, N) > 1, and
   whether β is well defined independently of how dependencies are lifted. Is Cor. B(a)'s "QS factor base" remark
   right (Jacobi vs Legendre, N ≡ 1 mod 4 for Blum N, reciprocity for odd ℓ)? Is beta0 a legitimate "constructor
   seeing only kN and N" (it tests smoothness of candidate v's)?
4. **§3.4.** Is the claim "no hidden-root guarantee is available within the I_k interface" proved, or only argued? Is
   the cost estimate L[2/3, 2/(3κ)] for one hidden-root relation right? Is Prop R′ correctly transcribed from
   nfsfruit Prop R (sampling with/without replacement, the 2^{π(y)} union bound, what x must contain)?
5. **The verdict and the requirement list C1–C7** (§5). Is anything missing or wrong (e.g. the oracle domain, the
   promise vs mean density, uniformity in k, smoothness tests, large primes, what "τ < 0.757" means)? Is "correct
   heuristic implication" an honest summary? Is anything overclaimed in §0?
6. **Code.** Any bug in the sieve's no-miss threshold, the constructors, primitive-pair dedup, Z = W/V, the F₂
   elimination, or the scoring (e.g. P2's parity model, the KS against a truncated exponential, R5's prefix statistic)
   that would change a reported number?

End with a table "finding | severity | fix" and a one-paragraph overall assessment.
