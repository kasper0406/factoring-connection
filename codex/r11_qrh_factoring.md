Research question (use web search; primary sources; label proved / conditional / conjectured / UNVERIFIED; be rigorous and adversarial).

Hypothesis [QRH]: OpenAI preprint "The Quasi-Riemann Hypothesis" (Oct 5 2026, github.com/openai/math, preprints/The-Quasi-Riemann-Hypothesis-October-5-2026). It claims that every Dirichlet L-function, including ζ, and every finite-order Hecke L-function over Q(√−3), has no zeros for Re s > 11/12. Hence there are no Landau–Siegel zeros. A companion paper claims Re s > 7/8. The paper lists consequences: least quadratic nonresidue ≪ (log p)^{O(1)}; deterministic Tonelli–Shanks; a deterministic Miller test; the effective bound h(D) ≫ √|D|/log log|D|. Treat [QRH] as an unverified hypothesis. Do NOT assume full GRH.

Question: what does [QRH] imply for integer factorization and the deterministic algorithms around it? Be quantitative: a zero-free half-plane Re s > 1−δ, with δ = 1/12 (or 1/8), gives polylog bounds with exponents depending on δ, not GRH's (log)^2. Work out each item, with explicit exponents where possible.

1. List the GRH/ERH-conditional results in computational number theory that bear on factoring. For each, decide whether [QRH] suffices, and with what loss. Cover at least:
   - Miller 1976 (computing φ or λ ⇔ factoring, deterministically, under ERH);
   - Bach's bounds (least nonresidue, generators of (Z/n)*, prime ideals generating class groups);
   - Shanks' SPAR / class-group factoring at N^{1/5} under GRH;
   - Schoof's class-group computations;
   - Hafner–McCurley's subexponential class group computation (GRH);
   - Buchmann's algorithms;
   - Lenstra–Pomerance (unconditional, randomized; does [QRH] help derandomize?);
   - Seysen;
   - Vallée;
   - deterministic factoring of special forms (e.g. Adleman–Huang / Lenstra results conditional on GRH);
   - Harvey/HH deterministic factoring (which step uses large-order elements; does [QRH] simplify it or change the exponent?);
   - Oznovich–Volk / HH 2026 order-finding;
   - deterministic primality certificates.
2. Key question: does [QRH] give any NEW unconditional deterministic factoring bound, either an exponent below 1/5 or deterministic subexponential (e.g. derandomizing Hafner–McCurley or Lenstra–Pomerance class-group relation collection)? Identify precisely what randomness those algorithms use: random group elements, random forms, smoothness of random elements. Then decide whether polylog generator bounds plus a zero-free half-plane suffice to replace it by a deterministic enumeration with provable smoothness yields. Smooth-number counts in short intervals or for norms from quadratic fields need what L-function input? Note that [QRH] covers Dirichlet L-functions and Hecke over Q(√−3), not Hecke L-functions over Q(√−N).
3. The deterministic quadratic-sieve missing lemma: B-smooth values among (⌈√N⌉+k)² − N for k < H, with B, H = L_N[1/2]. Does any zero-free-region input (Dirichlet, or Hecke over Q(√N)) imply the needed smooth-value count plus a useful dependency? What else is missing?
4. Output a ranked list of concrete theorems one could try to PROVE conditional on [QRH], each improving a known unconditional factoring-related bound. For each, give the precise statement, the proof strategy, and where [QRH] enters. Flag anything that would need Hecke L-functions beyond Q(√−3).
Output concise Markdown.
