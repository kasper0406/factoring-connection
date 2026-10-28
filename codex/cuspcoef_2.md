You are an adversarial referee (GPT-6.1-sol) for agent `cuspcoef` (phase 30) of a research project on classical
integer factorization. Repo: /Volumes/git/factoring-connection (read-only for you); cwd is a writable scratch dir;
python with sympy/gmpy2/python-flint/cypari2 is at /Volumes/git/factoring-connection/.venv/bin/python. Keep memory
< 2 GB, no background processes. Use web search to verify literature claims and give exact citations.

READ: notes/cuspcoef_results.md (the claims under review), notes/cuspcoef_prereg.md, src/cuspcoef/*.py,
tests/test_cuspcoef.py, experiments/raw/cuspcoef_*.jsonl. Background: notes/tracedata_results.md (§2.1 DEC, §2.6 SC,
§2.8 MC), notes/nativemod_results.md (§1.3, §3.3, §4, §5). Your own blind round is codex/cuspcoef_1.out (it is merged
into §8 of the note); an earlier, killed attempt of that round left scratch checks in codex/scratch/cuspcoef_1/
(trace_check.py, elliptic_trace_check.py, audit.log: PARI's qfbhclassno factors its argument internally).

STATE OF THE NOTE. It is v1 plus the blind-round merge (commit ad25398: the trace-formula route is now described as a
native N^{1/2+o(1)} FACTORING method, since its class-number sums output s = p + q). E2 is now complete (4 lines in cuspcoef_e2.jsonl; real/sim
acceptance 0.943/0.913, 0.400/0.335, 0.045/0.053, 0.005/0.0059 at n = 64, 128, 192, 256; 600 real, 1e5 simulated
samples each). Issues ALREADY KNOWN to the authors (do not just repeat them; say whether the proposed fix is right and
push further):
 (k1) DFB's right-hand side is "an algorithm given N alone" only if D~ = Phi(N, sigma, omega) is efficiently samplable;
      for the datum "factorisation of a_p a_q" the authors propose Kalai/Bach generation of random integers with known
      factorisation. They also propose a time corollary: since the bound holds for every N, expected time drops by at
      most a factor C^K.
 (k2) The joint spread of K forms is NOT B^K for isogenous curves or twists. The public family F (one curve per j)
      contains 34 isogenous pairs among its first 32000 curves, e.g. (A, B) = (1, 2) and (-19, 30).
 (k3) The Kronecker–Hurwitz relation sum_t H(4N - t^2) = 2 sigma(N) - sum_{d|N} min(d, N/d) = 2(q + N) for N = pq, so
      the class-number sums of the "native N^{1/2+o(1)}" Eichler–Selberg route reveal q at the same cost; and N^{1/2} is
      above deterministic factoring (Pollard–Strassen N^{1/4}, Harvey N^{1/5}), after which a_f(N) costs poly(n).
 (k4) cx1's MODEL bound (< 1.187K bits of congruence information, all odd l) is to be added as a complement to SC/SCX.

TASK: try hard to break every claim. For each item give a verdict (CORRECT / CORRECT WITH FIXES / WRONG / OVERCLAIM)
with a concrete reason, counterexample or computation.
1. Theorem DFB (§2.1): statement, proof, constant C = 2^{k+3}B^2, the corollary on K forms. Is the model LT(B)
   coherent? Is "determinant-free" defined precisely enough (e.g. does Phi get N but not p, q; can Phi use the
   factorisation of a_p a_q; is the ordered pair really allowed)? Is there any known reduction that DOES fit the
   hypotheses (that would refute the claim "every known reduction violates them")? Check the table in §2.3 row by row.
2. §2.2: the bound on H(D) and B_Gek (class number formula for non-maximal orders; L(1,chi) <= 3 + ln D; Robin); the
   claim B_E <= [GL2(Zhat):G] B_Gek in the Lang–Trotter law; the Serre-curve claim for 37a1; weight-12 remarks.
3. §2.4 DFB(m): statement and proof; the remarks on group-order coupling (is there in fact a reduction from a_N plus
   computations in E(Z/N) or E((Z/N)[t]/(g)) that factors N? try to find one).
4. Prop SCX (§3): openness of Lambda; the bound log2[Zhat^x : Lambda]; (iii) scalars in a Serre curve's image;
   joint scalars for K curves.
5. Prop RM (§4) and E1/E2: decoder correctness, leakage (operational code must not use p, q), the heuristic, the
   K = O(n^2 log n) claim, the cost claim L_N[1/3, (64/9)^{1/3} 2^{-1/3}], the E2 argument that easy-product selection
   is super-polynomially rare, and the registered-prediction bookkeeping (§6).
6. §5 native table and §0 item 3: the Eichler–Selberg cancellation (attribution, cost, and (k3): is "native" even a
   meaningful notion above the cost of factoring N outright?); your own 1.187K-bit bound (k4): check its proof once
   more, in particular the step from per-trace chi^2 to mutual information, independence across l for a FIXED curve
   (finite entanglement, e.g. Q(sqrt 37) in Q(E[2]) for 37a1), and whether it is stated conditionally on N mod l;
   Remark BC-native; Remark JS (verify the formula
   sum_{x mod N^2} (f(x)/N^2) = N (p - r_p)(q - r_q) and the factoring step).
7. Reconciliation: recompute every number quoted in §0, §2.2, §3, §4, §6 from experiments/raw/cuspcoef_*.jsonl and
   flag any mismatch (a previous draft of a report contained fabricated numbers; be strict).
8. Anything missing: a better reduction (R) or native algorithm (A) than what is stated; a sharper barrier.
Output: a numbered list of findings, most severe first, then a one-paragraph overall verdict
(ACCEPT / MINOR REVISION / MAJOR REVISION). <= 2000 words.
