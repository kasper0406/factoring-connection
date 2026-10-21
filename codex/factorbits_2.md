You are an adversarial referee (theoretical cryptography / computational number theory). Read
/Volumes/git/factoring-connection/notes/factorbits_results.md (read-only repo; scratch in your cwd). Code is in
src/factorbits/ (equiv.py, pure.py, gl.py), tests in tests/test_factorbits.py, raw data
experiments/raw/factorbits_e1.jsonl. Your own blind report is codex/factorbits_1.out.

Your job: try hard to BREAK the claims. For each item below say CORRECT / CORRECT WITH FIXES (give the fix) /
WRONG (give a counterexample), and check literature where relevant (web search allowed; cite precisely).

1. Theorem PT and its Partition lemma (§3.2). Is the model ("pure" trees: integer SLPs with constants independent of N,
   zero tests mod N, gcd instrumentation, outputs function of the leaf, any size, adaptive, randomized with N-independent
   coins) stated precisely enough that the lemma is true? Edge cases: division/inversion nodes, equality tests, square N
   (X = Y), primes dividing constants, tests whose value is 0 over Z, randomized trees, side information N mod m (PT-m),
   and trees that reuse ring elements across nodes (DAG SLPs). Is the inequality
   Pr[correct or split] ≤ max_b Pr[F=b] + P_split + η² right, and is "P_split ≥ η" right? Is the "advantage δ ⇒
   P_split ≥ δ/2" corollary right? Check Corollary PT2 (ε cross bound) line by line, including the identity
   ‖α−β‖² = 2P_split(B) − 2P_split(A) and the use of the same tree T₁ on population B. Is the generalisation to two fiber
   points stated correctly? Is the list of included algorithms right (Pollard rho? It uses gcd(x_i − x_j, N) with x_i a
   polynomial iterate: pure? Pollard–Strassen? ECM with random curves?). Is anything like PT already in the literature
   (generic ring model: Aggarwal–Maurer, Jager–Schwenk, Maurer–Wolf, Boneh–Lipton, Damgård–Koprowski; Lipton 1994
   straight-line complexity and factoring; Shoup's generic lower bounds)? Give the closest prior statement.

2. Lemma PSR (§4.2) — proof and hypothesis; is "randomized polynomial maps are covered coefficient-wise" right?

3. Prop GLM (§2.2): f(M) = parity(P⁺(M) AND M); M = N·k, k = r·N^{-1} mod 2^L. Check the claim P⁺(M) = q except for
   O(1/L) of r, the advantage accounting, and that the GL list decoding applies (r ranges over {0,1}^L uniformly; the
   predictor is evaluated on M = N·k(r), and Rackoff's pairwise-independent version queries r_J ⊕ e_i). Any prior
   formulation of embedding GL randomness into a cofactor?

4. §1.2: the call-count table. In particular (a) the randomized O(log log N) short-interval recursion (Huxley θ > 7/12:
   is the density statement for primes in [X − X^θ, X] correct and enough for an expected-polynomial Las Vegas search?);
   (b) the GRH claim for the one-call successor reduction (random t ≤ N·polylog N, tN − 1 prime with probability ≈ 1/ln):
   what exactly does GRH give for π(x; q, a) with q = N and x = N^{1+o(1)}? Is x ≥ q^{2+ε} needed? Fix the statement.

5. §1.5 partial forms: the quasi-polynomial bound for fixed-m exponents (recursion T(b) ≤ b + b·T(b/m)) — correct?
   Can you make it polynomial? The P⁺-only row. The approximate-addition row ("no harder than factoring one
   (n − log₂E)-bit number").

6. §2.3 status table: anything mislabelled (KNOWN vs OPEN), any missing known equivalence (e.g. results showing
   squarefreeness or ω is as hard as factoring under some reduction; Landau/Woll content; Shallit–Shamir), anything
   about p mod 4 / ε known?

7. §3.1/3.3: the trichotomy summary. Is "every algorithm faster than N^{1/4} is non-pure" right? (Pollard–Strassen is
   pure at N^{1/4}; is there any pure method below N^{1/4}?)

Deliverable: an itemised referee report (≤ 3000 words) with an overall verdict (ACCEPT / MINOR / MAJOR REVISION) and a
numbered list of required fixes.
