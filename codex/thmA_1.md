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

YOUR HALF: the ANALYTIC half. Check in depth:
  1. Theorem FM (pcest §2.3, with Lemma 8.5' in pcest §2.2 and the cited LV/Harper/Tenenbaum/Hildebrand facts in pcest §2.1) as re-checked for (a'''_k; mu) in ldbound §2. Check that LV Lemma 8.5 and Harper's Fact 8.11 (LV numbering) really say what is quoted, with the ranges used (x = N^{mu+h+d nu}, moduli t = l*s in [MV/32, 2MV], Q/l = N^{mu-nu}, B = L[1/2,a]). Check the goodness/grid argument, the induced law on t, the fibre/Möbius count, the truncation and sign arguments, and the bad-moduli bookkeeping.
  2. Lemma W and Corollary W (ldbound §3): the Brun–Titchmarsh form, the range t <= 2 M V_max < W, the claim that c_P does not depend on w.
  3. Theorem L2_fam (ldbound §4): the additive Parseval identity, the partition (A)/(B)/(C), each bound, the use of the large sieve and of the average over l, every margin.
  4. Theorem LD* (ldbound §5) Steps 0-7, including Lemma Kern (pcest §4.1) and Lemma DS (pcest §3). Re-derive the final power count V^4/Y^2 = S^4.
  5. Prop PZ (hypD §4) in the filtered form, and whether FM + LD* give D_k: gamma >= L^{-o(1)} uniformly in N (thmA_proof §4 Step 5, and §6 item U5 on the computability of S).
  6. Theorem R (pcest §3) is claimed NOT to be needed on the critical path; confirm or refute.
Also re-derive the exponent theta(d,k,mu) and the claim that k = 13 is minimal at d = 9 (you may run python; src/thmA/params.py has an exact implementation you should check independently).
