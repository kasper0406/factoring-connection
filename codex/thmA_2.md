You are an independent, sceptical referee. The repository root is /Volumes/git/factoring-connection (read-only for you; scratch files go in your current directory or /tmp). A research project claims:

THEOREM A. A rigorous randomized algorithm factors EVERY integer N in expected time L_N[1/2, (1+theta)/sqrt2 + o(1)], with theta = theta(d,k,mu,eta0) < 1 for fixed odd d, fixed k, base exponent mu and slack eta0; at d = 9 the infimum over (k, eta0) is 0.9192 (k >= 13 needed for c < 1). This would be the first rigorous improvement on Lenstra–Pomerance's L[1/2, 1 + o(1)] (JAMS 1992).

The claim is assembled from internal, UNREFEREED notes. Read FIRST notes/thmA_proof.md (the statement, algorithm and a roadmap citing each lemma), then the underlying notes it cites:
  notes/lemmaM_results.md (§1 Theorem DD), notes/hypD_results.md (§§2-4), notes/pcest_results.md (§§1-4.1),
  notes/lv6_results.md (§§1-4, 6.2, 7), notes/ldbound_results.md (§§1-6.1).

BLINDNESS RULES (important): earlier reviewers have already commented on these notes. You must form your own judgement. Therefore:
  - IGNORE the "Verdict" section at the top of each note, every "Reviews" section (lemmaM §9, hypD §9, pcest §7, lv6 §9, ldbound §8), every "[proved]" / "passed both critics" / "CORRECT" label, and every candidate card. Treat every claim as unverified.
  - Do NOT open anything in codex/ other than this prompt, do NOT open notes/*review*, notes/thmA_audit.md (if it exists), STATUS.md, REPORT.md or README.md.
  - You may read source code under src/ and the primary literature (web search is available; use it for the statements of Lee–Venkatesan arXiv 1805.08873, Harper arXiv 1208.5992, Lenstra–Pomerance 1992, Montgomery–Vaughan 1973, Lagarias–Odlyzko, Stark 1974, etc.).

Write your report as plain text/markdown. For every result you check give a verdict CORRECT / CORRECT WITH FIXES / GAP / WRONG with the exact location and a precise reason (a counterexample, a parameter computation, or a quoted hypothesis that fails). Distinguish carefully between (i) a real mathematical gap, (ii) an imprecision with an obvious fix, (iii) a cited external theorem whose hypotheses you could not confirm. Re-derive the key estimates yourself rather than trusting the notes' own checks. Look especially at INTERFACES: does each lemma's hypothesis actually hold in the FINAL design (a'''_k; mu) (F(g0,g1) = N*w, w a random prime, M = N^mu free), and is its conclusion exactly what the next step consumes? Look for any hidden L[1/2]-size loss (a factor L^{epsilon} or N^{o(1)} that is not L^{o(1)}), and any non-uniformity in N.

END your report with:
  WHOLE-THEOREM VERDICT: one of ACCEPT / ACCEPT WITH FIXES / MAJOR GAP / FALSE (for Theorem A as stated in notes/thmA_proof.md §2).
  WEAKEST LINK: the single step most likely to be wrong, and why.

YOUR HALF: the ALGEBRAIC and ALGORITHMIC half. Check in depth:
  1. Theorem DD (lemmaM §1) and the transcript-invariance hypothesis (T) for the WHOLE final algorithm in notes/thmA_proof.md §3: field selection, w, g0, g1, c, the Gauss split, the Lenstra–Pomerance shape test, cached special-q atoms, batch smoothness, LV6 columns/characters/kernel sampling/square root, the choice among s1, s2, s1+s2, restarts after a trivial gcd. Does ANY exposed quantity or decision depend on a hidden root r rather than on R = r^2? Is the probability exactly 1 - 2^{1-omega}?
  2. Theorems S and S' (hypD §2) and Theorem K(i)-(iii) with the field-selection rule of Corollary K (hypD §3), which were written for earlier designs (1-D field x - m, then (a''_k)). Are they valid verbatim with F(g0,g1) = N*w, free M = N^mu, the filtered candidate set with b prime, and G_A defined via the capped deterministic search? Check the expected-cost bookkeeping (E[1/q_A | kept] <= 1/mu), the 'at most k distinct primes' acceptance rule, Pr[shape]/Pr[success] = L^{o(1)}, and the claim that the shape test via Lenstra–Pomerance costs L_N[1/2, 1/sqrt2 + o(1)] per trial (check what LP 1992 actually proves: complete factorization? expected time for every n?).
  3. Prop LV6 (lv6 §§2-3: Lemmas C, H, E, X, linear algebra, Lemma SQ, the congruence) and its transfer to (a'''_k; mu) (ldbound §6.1): the row identity (lv6 (1.1)) including signs, the parity column, content e | N w, w | f_d or w | disc f, the Stark/Landau |E| <= 1 argument, the Lagarias–Odlyzko character sampler, the Hensel square root.
  4. Field-level events and Lemma IR (lv6 §4; ldbound §1 order of operations and §6.1): rejection vs factor exits, and the claim E[g_A 1_{not rejected}] >= gamma/2.
  5. Prop BS (lv6 §6.2): does the batch test recognise exactly the counted candidates N_l(A)? Cost.
  6. Cost assembly (thmA_proof §4 Step 10, lv6 §7, ldbound §6): re-derive min_a max(2a, a + 1/sqrt2, a + theta/(2a)) = (1+theta)/sqrt2; look for any hidden L[1/2] loss (rows, field selection, LV6, characters, Hensel prime, product trees, recursion).
  7. Uniformity in N: small prime factors, prime powers, perfect powers, non-squarefree N, omega(N) large, the cofactor after trial division, N0, effectivity, and whether "o(1)" vs "infimum over k" is stated correctly (thmA_proof §2, remark S1).
