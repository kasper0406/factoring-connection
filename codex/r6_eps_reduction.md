Research task (use web search, primary sources; label claims proved / conjectured / UNVERIFIED). Project context: /home/knielsen/code/factoring-connection (read REPORT.md §1.1 and §5, and NEXT_STEPS.md item 3; read-only).

Setting: D_n = balanced semiprimes N = pq with N ≡ 1 (mod 4), and ε = (−1/p) = (−1/q). Deciding ε is the sum-of-two-squares (SOS) decision. Question: is there a randomized polynomial-time reduction that factors N ← D_n with probability 1/poly(n), given an oracle A with Pr_{N←D_n}[A(N) = ε(N)] ≥ 1/2 + δ? The oracle may be queried on arbitrary integers, but it is guaranteed only on D_n.

Tasks:
1. Literature map with exact statements: Adleman–McCurley 1994 (open problem O11 and its neighbours); any result on the complexity of deciding SOS or of fixed-element quadratic residuosity (−1, 2, −2) for composite N; worst-case/average-case results; Rabin 1979 (square roots ⇔ factoring); results that use oracles on related moduli (kN, N + small, N^2, N·r for a known prime r); "Ireland–Rosen"-style genus theory; and the relation between deciding ε and computing a square root of −1 mod N. Note that SOS search ⇔ a sqrt of −1, but two independent sqrts factor N.
2. The reduction map. For each natural attempt, state exactly where it fails, with a derivation or explicit counterexample. Attempts:
   (a) query A on kN for small primes k, or on k·N with k ≡ 1 (mod 4) of known factorization;
   (b) query on N·r^2;
   (c) use A to decide (−1/·) for related products and combine via multiplicativity;
   (d) Jager–Schwenk-style re-randomization (why it needs a hidden input element);
   (e) a "search-to-decision" reduction that obtains a square root of −1 using the decision oracle on varied moduli;
   (f) quadratic forms: using A to decide ε for discriminants related to N, then the class-group / ambiguous-form route.
   Note: the oracle only promises correctness on D_n, i.e. on balanced semiprimes ≡ 1 mod 4. Queries outside D_n get arbitrary answers. This severely limits (a)–(c); quantify it.
3. Is there a relativized or black-box separation (an oracle world where ε is easy and factoring hard), or a known barrier? If none is known, say so.
4. Give your best assessment of what is plausibly provable in 1–3 agent-days. Propose ONE precise lemma that a research agent should attempt, plus the experiment that would test it.
Output Markdown, concise.
