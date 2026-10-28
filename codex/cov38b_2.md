You are an adversarial referee (analytic number theory, Fourier analysis on Z/q, geometry of numbers). Your job is to BREAK the claims, not to agree. You have read-only access to /Volumes/git/factoring-connection; write scratch scripts in your current directory (python3 / numpy; /Volumes/git/factoring-connection/.venv/bin/python has numpy, sympy, gmpy2). Keep runs under ~10 minutes and ~2 GB.

Read notes/cov38b_results.md (the claims under review), with background notes/coverage38_results.md (Theorem U, Prop L, section 3.1) and notes/dchworst_results.md (Theorems R, W). Your own blind round on the same questions is codex/cov38b_1.out; several of its results were integrated (Lemma F2 = your eq. (16)-(18); the moment correction; the Kerr-Shparlinski-Wu-Xi row). Check that they were integrated CORRECTLY and with correct attribution.

Referee in particular:
1. Prop T (target, reciprocal form, uniform form, constants) and Prop 0 (no intermediate exponent; averaged forms; the 26th-moment and 50th-moment numbers; the claim that even the 4th centred moment over a is a Kloosterman-type problem at this length).
2. Lemma F2 as written (Moebius decomposition of 1_S, the per-m bound min((h/m)A, K(1+log h)), the divisor count, the assembly to (K + h log q) q^{o(1)}). Is the q^{o(1)} honest, and are the hypotheses (2hK < q) enough?
3. Theorem B: (i) sufficiency via Fejer sums (supp G inside the coverage-implying set? W(G) <= H/4?), (ii) necessity for EVERY real g <= 0 off S_out (the definition of a valid S, the identity used, the t = 0 term), (iii) the exponent bookkeeping delta* = (1+2c)/(7-2c) vs the Parseval cap 1/(7-2c), with log factors. Is "for every c > 0" justified? Is the claimed exactness ("iff") justified up to q^{o(1)}?
4. Lemma F (energy route), Lemma F1, Remark B2, Prop B' (exact energy 2J^2-J when 16Y^3<q; Paley-Zygmund), Lemma M, Lemma S.
5. Theorem TP (two-pronged): union design, Lehman fraction (X/P^{2/5})^{5/4} J_L^{1/4}, linear families, the (x, floor(alpha x)) survivor, and the wording of what is assessment vs theorem.
6. The tables in sections 2.2 and 2.4: any wrong reference, range, or overclaim ("every tool", "no tool gives positivity"). Is there ANY known positivity (every-residue) theorem at total product size q^{1+o(1)} with one factor a prime in a short interval that the note missed?
7. The experiments section 5 against experiments/raw/cov38b_e*.jsonl and experiments/analysis/cov38b_summary.json (recompute a few numbers), and the pre-registration notes/cov38b_prereg.md (was anything run that deviates from it?).

Output Markdown: verdict first (ACCEPT / ACCEPT WITH FIXES / REJECT), then a numbered list of issues labelled MAJOR / MINOR with exact locations and concrete fixes, then anything new you found (with PROVED / UNVERIFIED labels).
