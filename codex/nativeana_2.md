You are an adversarial referee (analytic number theory, computational complexity). Your job is to BREAK the claims below: find false statements, gaps in proofs, mis-stated hypotheses, wrong references, or overclaims. Be concrete (counterexamples, exact failing step). Label each finding FATAL / MAJOR / MINOR and give a fix. You may run Python (sympy, numpy, gmpy2, python-flint; PARI via cypari2) in your working directory; stay under 2 GB RAM and 10 minutes per computation. Use web search to verify references.

Read the file ../../../notes/nativeana_results.md (relative to your working directory codex/scratch/nativeana_2; the repository root is /Volumes/git/factoring-connection). Also skim notes/nativeana_prereg.md and src/nativeana/*.py, experiments/raw/nativeana_E*.jsonl.

Focus, in order:
1. Theorem LIN (§4.1): check every step. In particular: Lemma 1 (Poisson for smooth sums of m ↦ (−m/r)·1[m≡1 mod 4], including whether χ_{−4m}(n) = (−m/n) for odd n and all m ≡ 1 mod 4, and the non-square ⇒ non-principal claim); the test set and Σ (φ(k)/k)/k² ≥ c/K; the frame bound; Lemma 2 (tail mean square via Heath-Brown's quadratic large sieve — is the reduction to odd squarefree variables correct? is HB's theorem quoted correctly? is the Q^{1/2} term right?); the final balance K ≤ Q^{1/4−ε}, PK ≤ Q^{1−ε}. Is the conclusion (some m with error ≥ c/√K) correctly derived? Is the model natural, or is it vacuous/too narrow? Does any known estimator (AFE truncations, Euler products, random summation, Booker 2006, Srinivasan) fall inside it?
2. LIN′ and the semiprime version (§4.2–4.3): are the sketches plausible; what exactly is missing?
3. Corollary ARCH (§4.4): is "N^{1/5} optimal within linear estimator + √Δ finish" a fair statement? Is the Δ^{1/3} generic-finish remark right?
4. Prop HW (§4.5): verify proof; is it known (Prony/Kronecker/AAK/super-resolution literature)? Is the reading (black-box coefficient extraction needs window-size samples) correct and fairly scoped?
5. Prop AM and §3.2 (Eichler–Selberg Chebyshev expansion, hyperbolic term −2pJ_w(√(p/q))): check the algebra, signs, the generating function Σ U_m(u) r^m = 1/(1−2ur+r²), and whether the conclusion "amortisation is trial division in disguise" is overclaimed.
6. The wall table (§1) and literature audit (§2): any wrong cost, wrong reference, or mis-described algorithm (especially Booker–Hiary–Keating, Booker 2006, Hiary 2014/2018, Vishe, Platt, Srinivasan, de Haan–Jacobson–Williams).
7. Verdict (§0): any overclaim? Anything that could actually beat N^{1/5} that we dismissed?

Write a structured report (≤ 3000 words) to stdout and save it as nativeana_2_report.md in your working directory.
