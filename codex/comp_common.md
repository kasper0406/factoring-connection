# comp: the compensation lemma (common preamble)

Context: notes/dd_results.md (just written), notes/chase_results.md (Lemma 2), experiments/chase/ch_L2c/proof.txt
(threshold inequality (1)–(2); "necessary combination" ~lines 42–51), experiments/dd/codex/proof.txt (cancellation
identity G_a = V_a + √2 W_a; the bound u(λ + log Cm) ≤ H + D log(Cm); the half-embedding construction with all capped
rows scalar), experiments/dd/gemini/proofs.txt (the (1 − (m−2)/n)-fraction construction; the identity
w_{2+j}(α_i) = ∏_{l≠i}(σ(r_j) − α_l) and the ρ-inflation Q ≥ Ω((1−m/n) D m² log m); the low-H family C1b with Q < 0
but log|F| ≥ 1.24 log N), notes/subthird_results.md (Theorem KH ideal-relative; LP_exc^rel).
Notation: Q = Σ_i log t_i − log(N𝔠N𝔠′); H = Σ_i log M_i − log J (ideal-relative envelope); T = 2⌈2n/s⌉+1 (KH);
balance: fixed k, n ~ √(x/log x), m ~ x^{1/4}(log x)^{3/4}, ternary box, cap log Y ~ √(x log x).
Target (COMPENSATION LEMMA): for reduced linked O_k-modules with H = O(x/T), show that if log(M_i/t_i) ≥ λ at u ≍ D
embeddings with λ ≍ m log m, then EITHER Q ≥ c_k x/T (so every capped candidate is confined to ≤ D specified proper
subspaces — then charge their parity rank via chase Lemma 2's fresh-rank bound) OR every non-scalar capped candidate
has log|F| ≥ (1 − o(1))·(KH floor) (so N_ind is tiny). I.e. bound the contribution of the compensating conjugates to
Σ log t_i. Natural tools: the per-row resultant integrality (product over 2n conjugates ≥ 1) combined with the
per-embedding thresholds across a full ternary basis (the basis attaining t_i), Hadamard/Cramer at the compensating
places, KH on V^h restricted to the loss-aligned sub-configuration, and Gemini's ρ-inflation identity generalized
beyond linear poles.
You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch only in your cwd).
Output ≤ 1800 words, plain language, labels proved / heuristic / measured / speculative, file:line cites; long proofs
to a file in your cwd; pre-register experiments. Verdict: proved / refuted (with the counterexample family and its
yield ledger) / partial (what remains) + next step.
