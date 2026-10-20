You are a senior research mathematician in computational and analytic number theory (index calculus, NFS, algebraic number theory, rigorous algorithm analysis). Work read-only on /Volumes/git/factoring-connection; your cwd is a writable scratch dir. Use web search freely and cite exact statements (paper, theorem number, year).

CONTEXT. Read notes/thmA_proof.md first (Theorem A: a claimed rigorous Las Vegas factoring algorithm in L_N[1/2, (1+θ)/√2 + o(1)], with infimum 0.919 at degree 9 and → 1/√2 as the degree grows). Its sources are notes/{lemmaM,hypD,pcest,lv6,ldbound}_results.md, and notes/thmA_audit.md is an internal audit. Treat all of these as internal and unrefereed.

Theorem A has two parts that look reusable for discrete logarithms:
- relation supply for special-q descent of RANDOM targets (Theorems S/S′, K, FM, Lemma Kern, Lemma W, Theorem L2_fam, LD*, Theorem 4.3/PZ): a uniform target R mod N is split R ≡ u/v with u, v ≈ √N, the shape test factors u, v with Lenstra–Pomerance, and every prime factor ℓ ≈ N^{1/(2k)} is descended by one special-q step in a randomized degree-d field with F(g₀,g₁) = N·w;
- NFS post-processing (lv6: GF(2) kernel with quadratic characters, Lagarias–Odlyzko, Stark/Landau exceptional class, square root in Z[β]).
The factoring-specific part (fruitfulness via hidden square roots, Theorem DD) is not needed for discrete logarithms.

TASK (blind round; nobody has worked this out for you).

1. LITERATURE. What are the best RIGOROUS (unconditional, and separately GRH-conditional) bounds for:
   (a) discrete logarithms in F_p^* (p prime) — Pomerance 1987 gives L_p[1/2, √2]; has anything rigorous beaten √2 for prime fields, for all p or for almost all p?
   (b) F_{p^n}, including small/fixed characteristic (Kleinjung–Wesolowski quasi-polynomial, 2019/2021/JAMS 2022), Lovorn Bender–Pomerance, and any rigorous L[1/2] or L[1/3] results in medium/large characteristic;
   (c) class groups of imaginary quadratic fields (Hafner–McCurley 1989 under GRH, Vollmer, others);
   (d) Jacobians of hyperelliptic curves (Enge, Enge–Gaudry, Heß, Diem), with exact constants and hypotheses.
   Give exact statements and constants. Say explicitly whether any rigorous prime-field DL beats L_p[1/2, √2].

2. DESIGN. Design a rigorous NFS-type index-calculus algorithm for DL in F_p^* that reuses Theorem A's relation supply with p in the role of N (two randomized sides, F(g₀,g₁) = p·w, special-q descent of uniform targets g^a h^b). Address in detail:
   - Which steps of thmA_proof Steps 1–7 transfer verbatim when N = p is prime, and which need changes (e.g. Theorem S for R uniform on F_p^* instead of on squares)?
   - The linear algebra: modulo the prime factors ℓ of p − 1, including small ℓ and prime powers ℓ^e, or otherwise. Must p − 1 be factored?
   - The number-field obstruction: units and class group modulo ℓ-th powers, Schirokauer maps, ℓ-th power residue characters, virtual logarithms. What can be made RIGOROUS at fixed degree d, for ℓ small and for ℓ large (ℓ up to p/2)? Does lv6's square-class analysis transfer to ℓ-th power classes? Is there a way to avoid the obstruction entirely?
   - How the answer log_g h is extracted (fruitfulness analogue), and whether it needs a generator g.
   - The total cost and the constant c, as a function of θ, with exact arithmetic for d = 9, k = 13 and the infimum; any new cost term (e.g. from linear algebra) that can bind.
   - Uniformity in p.
3. VERDICT. Does this give a rigorous L_p[1/2, c] with c < √2 (or c < 1)? If not, what exactly is missing (state the missing lemma precisely)? If yes, state the theorem precisely, with every hypothesis, and identify the weakest step.
4. RECURSION. The shape test factors u, v ≈ √p. Comment briefly on how a recursive shape test (factor the halves with the algorithm itself) would affect the DL constant; do not develop it in depth.

Be adversarial with your own proposals. Label every claim KNOWN (with reference) / PROVED (give the proof) / HEURISTIC / SPECULATIVE. Output Markdown, verdict first.
