# inversesieve_2 — adversarial referee (GPT-6.1-sol, read-only)

Referee `notes/inversesieve_results.md` (repo root = cwd's parent's parent; read-only). Also read
`notes/inversesieve_prereg.md`, `src/inversesieve/*.py`, `tests/test_inversesieve.py`,
`experiments/analysis/inversesieve_summary.json`, and your own blind round `codex/inversesieve_1.out`.

Be adversarial. For each labelled claim (Theorem S (a)-(f), Prop C, Prop E, Prop M, Theorem L and
Corollaries L1-L3, Prop D, Cor. D2, Prop Q, Prop I table, every experimental verdict in §7):
1. Is the proof correct and complete? Give counterexamples or missing hypotheses (e.g. Theorem L: cell
   alignment, boundary cells, gcd(x, M_i), sub-rule count, the Δ_i ≥ 2M_i reduction, Theorem-G charging,
   Theorem N's hypotheses in L1; Cor. L3's arithmetic; Prop D's inequality and the J^{4/5} conversion;
   Cor. D2's cost claim).
2. Are the exponents right (6/31, 7/124, 69/124, 55/62, 4/31, k ≥ 7/2, 3/124)?
3. Is any claim overstated relative to its scope (e.g. "kill", "equivalent", "no better than brute-force
   Coppersmith")? Is the verdict justified?
4. Experiments: do the frozen analyses match the pre-registration; are the disclosed deviations and the
   exploratory follow-ups (matched controls, fresh-seed replication, held-out deficit) honestly labelled;
   is the "U(N) = {p} in all 88 large-J worst pairs, selection bias" reading right?
5. Anything you can verify independently with a short exact script (write in your cwd), do so.
Output: a table claim | verdict (CORRECT / FIX / WRONG) | required fix, then an overall verdict
(ACCEPT / ACCEPT WITH FIXES / REJECT). ≤ 2500 words.
