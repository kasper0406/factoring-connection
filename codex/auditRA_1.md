You are an independent, sceptical referee. The repository root is /Volumes/git/factoring-connection (read-only for you; scratch files go in your current directory or /tmp). Your job is to find out whether a claimed result is TRUE. A negative or "major gap" finding is just as valuable as a positive one; do not be generous.

A research project claims (internal, UNREFEREED notes):

THEOREM RA. Fix odd d >= 3, k, mu, eta0 admissible for the project's "Theorem A" (theta = theta(d,k,mu,eta0) < 1). Replace Theorem A's Lenstra–Pomerance shape test (complete factorisation of the two Gauss-reduced halves |u0|, |v0| <= sqrt(2N) of R = r^2 mod N) by a recursive call of the algorithm itself on |u0| and |v0|, with fresh coins, and take the smoothness bound B = L_N[1/2, a] with a > 0 a FREE fixed parameter (Theorem A used a = theta/sqrt2). Then for every fixed a > 0 the algorithm is Las Vegas, completely factors every N, and runs in expected time L_N[1/2, c(a) + o(1)], c(a) = max(a + theta/(2a), (2+sqrt2) a). Optimum c* = sqrt((1+sqrt2) theta): 0.851 at d = 9 (as an infimum over k, eta0), tending to 0 as d -> infinity.

COROLLARY P. A single Las Vegas algorithm (a restart schedule over a computable sequence of parameter sets with c*_j -> 0) factors every N in expected time exp(o(sqrt(log N log log N))) = L_N[1/2, o(1)], at an ineffective rate.

These would beat Lenstra–Pomerance's rigorous L[1/2, 1] (JAMS 1992) dramatically. Extraordinary claims need extraordinary scrutiny.

READ:
  - notes/auditRA_blind_claims.md: the claims, algorithm, proofs of RA (Lemmas OS, W, F, T_rec) and Corollary P. This is a verbatim excerpt.
  - notes/thmA_proof.md: Theorem A's statement, algorithm and proof roadmap (the base on which RA is built).
  - The source notes as needed: notes/pcest_results.md (FM, Lemma DS, Theorem R, Kern), notes/ldbound_results.md (Lemma W, L2_fam, LD*), notes/hypD_results.md (Lemma Omega, Theorems S, S', K, Prop PZ), notes/lv6_results.md (LV6, Prop BS), notes/lemmaM_results.md (Theorem DD, section 1), notes/lem85_results.md (Lemma 8.5').
  - src/thmA/params.py and src/recurse/constants.py (exact constants), if useful.

BLINDNESS RULES (important): earlier reviewers have commented on these notes; form your own judgement.
  - Do NOT open notes/recurse_results.md, notes/recurse_codex1_report.md, notes/thmA_audit.md, notes/auditRA_audit.md, STATUS.md, REPORT.md, README.md, or anything in codex/ other than this prompt.
  - In the source notes IGNORE every "Verdict" section, every "Reviews" section, every "[proved]" / "CORRECT" / "passed critics" label and every candidate card. Treat every claim as unverified.
  - You may take Theorem A's proof (at a = theta/sqrt2, with the Lenstra–Pomerance oracle) as GIVEN, but you must check whether each source estimate it uses really holds for EVERY fixed a > 0 and whether RA's use of it is legitimate.
  - Web search is available: use it for the literature questions.

CHECK IN DEPTH:
 1. a-freedom. Does every estimate used (FM, Theorem R, Kern, Lemma W, L2_fam, LD*, Theorem 4.3 / Prop PZ, S/S', K, Lemma Omega, Lemma 8.5', LV6, Prop BS) hold for every fixed a > 0, not only a = theta/sqrt2? Look at the parameter constraints linking a, B, S^2, V, H, the h-floor, constants such as C0 (B^{C0} >= 64 S^2 L / (rho rho)), C1, K = 6 + ceil(2 xi/a^2), and any "for N >= N0" or level-of-distribution range (e.g. y >= (log x)^K, u/log^2 u >> log log x, conductor <= y^K <= x^eta). Do small a (B tiny, S huge) or large a (B huge, S tiny) break anything? Re-derive the cost ledger for general a.
 2. The recursion. The expected-time recurrence F(x) <= 2 log2 x [P(x) + 2 D(x) F(sqrt(2x))], Lemma W (Wald), Lemma F (solving it; check the induction and constants carefully), the strong-induction finiteness, inputs <= sqrt(2N), the factor-tree accounting, the number of shape tests per relation, whether the halves can be 1, primes, prime powers or tiny, and whether failures/restarts of the inner calls are accounted for. Is it legitimate that Theorem A's bound on the number of DESCEND_ONCE calls is independent of the oracle?
 3. Transcript invariance (T) of Theorem DD under recursion (Lemma T_rec), and Lemma OS (oracle substitution): is "swapping in any exact factoring oracle leaves the law of the outer transcript unchanged" valid when the oracle's running time is random? Are there caps, timeouts, races or N0 branches (e.g. "if N' < N0: trial division", with N0 perhaps not computable) that could break it? Could the recursion's own randomness or running time leak information correlated with the hidden roots r?
 4. Corollary P. Is a restart schedule over a growing sequence of parameter sets a legitimate "single algorithm with expected time L[1/2, o(1)]"? Check the Markov/phase computation, whether the parameter sequence is admissible and computable for all j (does eta0_j -> 0 need to be imposed? does theta_j -> 0?), whether the fixed-parameter thresholds N0(d_j, ...) break the argument, what exactly is ineffective, and whether finite E[T_j] suffices (the "E[T log T]" issue for live dovetailing).
 5. Sanity and literature. If rigorous L[1/2, o(1)] factoring followed this easily from a rigorous relation supply, why has it not been noticed? Is there a known (published or folklore) obstruction saying rigorous methods relying on smoothness of sqrt(N)-size auxiliary numbers cannot beat L[1/2, 1] (Pomerance's analyses of random-squares-type methods; Croot–Granville–Pemantle–Tetali; Lenstra–Pomerance p. 509)? Does RA evade it legitimately? Compare with the heuristic fixed-degree NFS constant (what L[1/2, c(d)] does heuristic NFS with fixed degree d give?), to check that RA does not claim more than the heuristic. Search for prior art on recursive use of a factoring algorithm on auxiliary numbers.
 6. Re-derive c(a), a*, c*, the fixed-depth iteration c_{j+1} = (theta/c_j + c_j)/sqrt2, the threshold c* < 1 iff theta < sqrt2 - 1, and c* at d = 7, 9, 11, 13, 19 (infimum theta = 3/(d+1)) and at d = 9, k = 13, with exact arithmetic where possible.

Write a plain markdown report. For every item give CORRECT / CORRECT WITH FIXES / GAP / WRONG with the exact location and a precise reason. Distinguish (i) a real mathematical gap, (ii) an imprecision with an obvious fix, (iii) an inherited Theorem A issue.

END your report with exactly these lines:
  VERDICT RA: one of ACCEPT / ACCEPT WITH FIXES / MAJOR GAP / FALSE (conditional on Theorem A at a = theta/sqrt2 being correct)
  VERDICT P: one of ACCEPT / ACCEPT WITH FIXES / MAJOR GAP / FALSE (same condition)
  WEAKEST LINK: the single step most likely to be wrong, and why.
