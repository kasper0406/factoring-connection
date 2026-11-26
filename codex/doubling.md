# doubling: try to evaluate the central-binomial zero test mod N fast (builder role)

You are the **builder**. The user wants this attempted seriously; they are ~90% sure the hidden-class bit has a
poly(log N) algorithm. Repo READ-ONLY; work in your cwd. Read notes/buildbit_results.md and
experiments/buildbit/ledger.md first (your previous session, 01a21773-2898-7433-9eac-dd535aa3988d, built the carry test).

## Object
B(N) = C((N−1)/2, (N−1)/4) mod N for N = pq ≡ 1 mod 4 (= C(2m, m), m = ⌊N/4⌋).
T0: B ≡ 0 mod N ⟺ p ≡ q ≡ 3 mod 4. Your proposed next step was a dyadic doubling rule for the zero image of
U_n = C(2n,n)/4^n, 2nU_n = (2n−1)U_{n−1}.

## Calibration from the orchestrator (check it; don't just accept it)
- **General index ⇒ factoring.** For n < p, p | C(2n, n) ⟺ p ≤ 2n. So a zero/gcd test of C(2n, n) mod N at
  ARBITRARY n gives p by binary search. Any doubling rule that works at all n is therefore poly-time factoring, and
  it would also evaluate factorial-type products mod N (best known ≈ n^{1/2}: Bostan–Gaudry–Schost, Strassen).
- **So the only non-factoring-equivalent hope is the special index (N−1)/4**, where the index is tied to N.
  Lucas in the 1,1 case gives B ≡ C_q · C((p−1)/2, (p−1)/4) ≡ C_q · 2a_p mod p. Here p = a_p² + b_p² (Gauss), and
  C_q = C((q−1)/2, (q−1)/4) is a cross term: an integer reduced mod p. Symmetrically mod q.

## Attempts (spend effort where it looks alive; record a yield ledger per attempt)
D1. **Special-index identities.** Can B(N) mod N be rewritten through something N-native? Candidates: Hasse
    coefficient of y² = x³ − x at "N"; the CM value E_{N−1}(i) (or the E_4^{(N−1)/4} coefficient of E_{N−1}); the
    Gross–Koblitz / Jacobi-sum form of Gauss's congruence; Dwork-type and supercongruences for C(2n, n) along n ≡ −1/4;
    a hypergeometric (₂F₁(1/4, 3/4; 1; ·)-type) truncation at N. Does any of these have an evaluation cheaper than the
    length-N recurrence, and does it keep the zero/nonzero distinction mod N?
D2. **Doubling for the zero image only.** Over F_s the zero set of C(2n, n) is "some base-s digit of n ≥ s/2".
    Look for a mod-N quantity Z(n) with a composition law Z(2n), Z(2n+1) ← Z(n) (or any poly(log) addition chain)
    whose vanishing mod p tracks that set at n = (N−1)/4. Examples: polynomial-matrix products mod N with a twist by
    N, Frobenius-free q-analogues (Gaussian binomials at roots of unity mod N), Lucas-sequence analogues, Carlitz-style
    digit sums. Prove each candidate either works or collapses to factoring or to an M-family order.
D3. **Amortisation.** Harvey's average-polynomial-time remainder trees amortise over many primes. Can a single N be
    embedded in a family (kN, N + k·M, ...) where the batch cost per member is poly(log N) and one member's answer
    gives T0 for N?
D4. **Experiment (pre-register first in prereg.md, with hypotheses, sizes and pass threshold, BEFORE running).**
    For 32–64-bit N, compute B(N) mod N exactly (BGS or plain), and test every cheap public function of it, or of your
    D1/D2 candidates, for agreement with T0. Also test whether (B/N) or B mod small numbers correlates with the 1,1
    case's hidden data (a_p, a_q). Ground truth from the factors; flag gcd splits.

## Output
≤ 2000 words. Per attempt: built / measured / killed, with the precise obstruction. If you find a working
poly(log N) evaluation at the special index, stop broadening and harden it (proof, 80+-bit tests, the lift to T1).
Label every claim proved / known / heuristic / measured / speculative.
