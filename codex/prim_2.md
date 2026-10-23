You are an adversarial referee (research mathematician / theoretical computer scientist). Repo:
/Volumes/git/factoring-connection (readable; your cwd is a scratch dir you may write to; python:
/Volumes/git/factoring-connection/.venv/bin/python with sympy, gmpy2, numpy). Be harsh and precise. Use web search to
check literature claims. Label each finding OK / WRONG / GAP / OVERCLAIM / MISCITED / MINOR, with a concrete fix.

TARGET: notes/prim_results.md (agent `prim`, v1), which builds on notes/leakmap_results.md (PT-X, NE1, PT-J),
notes/factorbits_results.md §3 (PT), notes/slpcirc_results.md §4 (Theorem SP, Kelley), and codex/prim_1.out (your own
earlier blind round; results credited there). Code: src/prim/*.py, tests/test_prim.py. Pre-registrations:
notes/prim_prereg.md and notes/prim_posthoc_prereg.md. Data summaries: experiments/analysis/prim_e1_*.json.

CHECK IN PARTICULAR (try to break each; give counterexamples where possible, numerically if useful):
1. Theorem U (§2): the Delta-form, floor gates, the "<= 2 primes per residue" counting, the domination inequality and
   the constant 8/pi(I) after PT-X. Is the reference kappa(p) legitimately a function of p alone in adaptive trees?
2. Lemma P (N-gates are permutations) and its use.
3. Theorem R (§3, NEW): modal decoupling Pr[correct or split] <= B + 4s + 2 Pr[D_mod and m_p = m_q] for programs with
   fresh coins from Z/N. Check every step: conditional independence of tau_p, tau_q given (p, q) when coins are
   shared across nodes and nodes are adaptive; the inequality s(p,q) >= 1/2 when modal transcripts differ; the
   PT-X step with distinct pairs; whether "transcript" must include the coin-dependent path; p = q.
   Then Corollary R-chain: the reduction from modal deviation to a single node with local zero-probability >= 1/(2h),
   the inverse-fiber bound, and the final estimate O(h k alpha tau (2hD)^{1/k} X^{-1/k} ln X).
4. Theorem PS (§4) built on your Lemma KS (12): Brun-Titchmarsh use, the dens(V_L) term, whether "removed cosets are
   cosets of a single subgroup" is right, and the claimed range of t.
5. Theorem D2 (§5, NEW): |Z(p)| <= gcd(alpha,d) tau(d) (n/d + d (n/o)(sqrt p + 2)) + n/o for nodes
   A (y^{ar+b} + c)^{alpha r + beta} - B with A, B monomials in gate outputs. Scrutinise the character expansion, the
   reduction to Sum_{h in <y^a>} psi(h) chi(y^b h + c), the Weil application (non-degeneracy, number of roots, the
   w = 0 term), the count |S_d|, and the corollaries (sieve for divisors of p-1 in [p^delta, p^{1/4}], GRH order
   input, "o(1) not negligible"). Is the extension to inner u = M1 (M2 + c) correct? Is there a known result that
   already gives D2 (e.g. via Bourgain-type sum-product bounds for Zech logarithms, or a better bound for the Zech
   hyperbola ind(u) ind(u+c) = t mod p-1)?
6. Lemmas FZ (frozen exponents) and CL (closure), and the claim that both reduce to (M)/(A) events.
7. Theorem J1 (§6): the Fourier bound with (2^{k-1}+1)s, degenerate/non-degenerate split, parity nu public, and the
   min(1, .) error terms.
8. §7 depth hierarchy: is "first break for deterministic bases at shift depth 2" justified, or does D2's method (or
   another) handle ((g^r + c1)^r + c2)^r? Conversely, is depth 2 really "mostly classified" given the exceptional set?
9. §9 experiments: are the pass/fail statements faithful to the registrations? Is the second-stage held-out test
   (b = 16) legitimate given that the post-hoc classes were developed on b <= 14? Any statistical problems (e.g.
   clumping of zeros in cosets inflating variance; the q-world model including PIN residues; dyadic variance
   negative)? Spot-check numbers against experiments/analysis/*.json.
10. Verdict (§0): any overclaim? Is "within PRIM every non-split leak we can analyse is (M), (A) or (Q)" accurate given
    what is proved vs experimental?
11. Try once more, briefly, to construct a PRIM leak beyond the fiber floor without splitting in any of the OPEN
    fragments listed in §0 (deterministic shift-depth-2 towers, exceptional primes, exponential-degree post-processing,
    deterministic gate Jacobi, compressed Jacobi). Report any candidate with a numerical check.

OUTPUT: a numbered list of findings with labels and fixes, then an overall recommendation (ACCEPT / ACCEPT WITH FIXES /
MAJOR REVISION) and the 3 most important changes.
