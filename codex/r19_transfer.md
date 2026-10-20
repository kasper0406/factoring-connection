You are a senior research mathematician in analytic and computational number theory. Work read-only on /Volumes/git/factoring-connection; your cwd is a writable scratch dir. Use web search freely.

READ FIRST:
- notes/thmA_proof.md: Theorem A. As far as an internal 5-referee audit can tell, it gives a rigorous Las Vegas factoring algorithm for every N in L_N[1/2, c], for every c > 1/√2. That would beat Lenstra–Pomerance's L[1/2, 1].
- notes/thmA_audit.md.
- The source notes: notes/{lemmaM,hypD,pcest,lv6,ldbound}_results.md.
- APPROACHES.md.

THE NEW DEVICES developed there:
1. Dixon descent: a random square r² is split as u/v and descended through special-q. It removes the NFS fruitfulness problem, because the hidden root stays hidden by transcript invariance.
2. The N·w multiplier: build f with F(g₀,g₁) = N·w, w a random prime, so that forced residues mod moving moduli become equidistributed (Brun–Titchmarsh).
3. A 2D-randomized rational side, g₁x − g₀ with (g₀,g₁) random (Buhler–Lenstra–Pomerance base-(m₁,m₂)). It turns binary smooth correlations into lattice problems.
4. Divisor switching plus the large sieve, giving all-residue L² BDH-type upper bounds at precision L^{o(1)}; Paley–Zygmund then gives existence.
5. Bernstein batch smoothness, which replaces HECM rigorously at no exponent cost.
6. Split counting by character positivity (Theorem S), and bounded-k splits (Theorem K).
7. (Proposed, not yet checked) recursive self-application. The shape test factors the halves u, v ≈ √N with the algorithm itself rather than with Lenstra–Pomerance, which would remove the 1/√2 floor.

TASK. Find where else these devices, alone or combined, give PROVABLE speedups. For each candidate, give:
- the target, with the best known rigorous bound and exact reference;
- the mechanism;
- the claimed new bound, with a feasibility calculation;
- the exact missing lemma;
- the first step;
- your probability of success and the significance.

Candidates to consider, not exhaustive:
- (a) Recursion (device 7): check the fixed-d fixpoint constant and the growing-degree limit. Could the HECM-based L[2/5] barrier of notes/lemmaM_results.md §2(c) and the one-level-descent limits now give a rigorous L[α] with α < 1/2? What exactly is the descent bottleneck?
- (b) Rigorous discrete logarithms in F_p^* (Pomerance 1987, L[1/2, √2]), in F_{p^n}, and in class groups and Jacobians. Has anyone beaten √2 rigorously for prime fields?
- (c) Rigorous class group / unit group computation for imaginary quadratic fields WITHOUT GRH (Hafner–McCurley is under GRH).
- (d) Rigorous L[1/2] computation of the multiplicative order, or of φ(N).
- (e) Rigorous "smooth relation" problems in other algorithms (ECM-type almost-all results).
- (f) The project's deterministic threads. Does the N·w multiplier or 2D randomization help derandomize, e.g. the deterministic QS H2 obstacle in notes/detqs_results.md, or worst-case coverage in notes/coverage38_results.md?
- (g) Anything else.

Be adversarial with your own proposals. Kill those that fail, and rank the survivors.

Output Markdown with a ranked table first. Label claims KNOWN (ref) / PROVED / HEURISTIC / SPECULATIVE.
