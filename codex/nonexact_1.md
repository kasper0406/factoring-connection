# nonexact round 1 (BLIND): non-exact / evaluated compressions of split polynomials modulo N — barrier or construction

You are a research mathematician (analytic number theory, geometry of numbers, lattice algorithms) collaborating on a
speculative project about asymptotically faster integer factorization. Rigour matters more than speed and honesty more
than a positive result. Label every claim PROVED (give the proof), KNOWN (cite precisely), HEURISTIC (name the
assumption), or OPEN. Do not read other files of this repository; work from this prompt, your knowledge and web
search. Scratch computation (Python; exact integer arithmetic) goes in the current directory.

## Notation and background

x = log N, L[β] = exp(Θ(x^β (log x)^{1−β})). N = pq balanced semiprime. A multiset A of M integers in [−B, B] has
the split polynomial P_A = ∏_{a∈A}(X − a). B = L[β], β < 1/3.

**Target (the "systematic-equation escape").** A family 𝓕 of pairs (u, v), u ≡ v (mod N), u ≠ v, with
(a) u = ∏ of ≳ x/log B pieces, each ≤ B^{O(1)}, by an identity ("smooth by construction");
(b) v small: log|v| = O(x^{2β}(log x)^{1−2β}) (then v is L[β]-smooth with probability L[β]^{−O(1)});
(c) |𝓕| ≥ L[β] with β < 1/3, computable from N alone at amortised cost L^{o(1)} per member.
In small-characteristic finite fields this exists (BGJT/Joux/GKZ/Kleinjung–Wesolowski): one presentation
X^q ≡ h₀/h₁ makes ≈ q³ split polynomials reduce to degree ≤ 3.

**Presentations.** F ∈ Z[X] monic of degree d, with φ: Z[X]/(F) → Z/N.
- View B: φ(X) = μ ∈ Z/N with F(μ) ≡ 0. Factor base {μ − a : |a| ≤ B} plus small primes.
- View A: φ(X) = m ∈ Z with N | F(m). Pieces are the integers m − a (or, with m = 0, small integers/primes).

**A member** is (A, h, Q, r): h, Q ∈ Z[X] multipliers (deg ≤ j₁, j₂, height ≤ H_Q), r ≠ 0 the small side
(deg < k, height ≤ H). Its relation is φ(h)φ(P_A) = φ(r).
- *exact*: hP_A = FQ + r in Z[X];
- *non-exact*: hP_A ≡ r only mod the ideal (N, F) of Z[X];
- *evaluated* (View A): only h(m)P_A(m) ≡ r(m) (mod N), r small in VALUE at m, not necessarily in height.

**Proved results you may use (exact case is closed).**
- Lemma R (rigidity). P, P' monic degree M, split over Z, roots in [−B, B]; Q, Q' ≠ 0 of degree ≤ j;
  G = QP − Q'P' of degree < K, height ≤ h. If M ≥ K + 2⌊log₂(2^{K−1}√K h B^{K−1})⌋ + 4j + 3 then G = 0.
  Proof: remove the gcd; Mahler-measure bound |G̃(a)| ≤ Δ at integer points; at a root a of U, |U'(a)| = ∏|a − a'|
  ≤ Δ so at most log₂Δ roots of U' are at distance ≥ 2 (discreteness of Z in R); vanishing moments of the
  difference measure force support ≥ M − K + 2 (Vandermonde); count supports.
- Theorem W. For exact members on any F (any height): in that regime hP_A/r is the same rational function for all
  members — one presentation gives one relation.
- Lemma EX. If u ≡ v mod (N, F) (F monic) and R₀ = (u − v) mod F ≠ 0 then N | R₀ coefficientwise, so ‖R₀‖ ≥ N.
- Prop NX. In a PTE fibre (members share e₁..e_{M−k}), d ≥ k, non-exact members have pairwise distinct reduced
  defects (in Lemma R's regime). BUT injective defects do not make the congruence conditions independent: e.g.
  P = X²(X² − 49)², P' = (X² − 64)(X² − 25)(X² − 9), P − P' = 14400, N = 25: both members compatible because
  25 | 14400. For k = 1 such an inherited compatibility is itself a relation (14400 = ∏(b − a) is 2B-smooth and
  ≡ small mod N).
- Kleiman: (M − 1)! | P − P' for an ideal PTE pair. Fibres with one free coefficient have ≤ c B^{1/2}(log B)^κ
  members (Binyamini–Cluckers–Novikov).
- A cost-model fact: an output-sensitive enumerator of {split P : P(m) ≡ small mod N} with output ≥ L[β] at cost
  |output|·L^{o(1)} would suffice for sub-L[1/3] (and a "doubly-smooth" lister would give heuristic polynomial
  time). Schnorr/Adleman prime-number-lattice factoring is a known failed instance of this idea.

## The open remainder O1 (your task)

Compressions that hold only modulo N (or modulo an N-dependent ideal), or remainders small only in VALUE at m.
Lemma R and Theorem W refute EXACT compressions in their regime; O1 is the rest.

1. **Formalise O1** precisely. Distinguish (i) "holds mod (N, F)" from (ii) "small in value at m". For each, show
   exactly which step of Lemma R / Theorem W breaks and why (quantitatively, in the L-range M ≍ x/log B).
2. **Barrier attempt.**
   - Lattice view: for fixed split P, {(Q, r) : QP ≡ r mod (N, F)} is a lattice; compute its determinant,
     Minkowski / Gaussian-heuristic short-vector sizes, and a rigorous Davenport / successive-minima count.
   - Charge the split structure honestly: P = ∏(X − a_j) with small roots is a NONLINEAR (thin) condition on the
     coefficients. Can Bombieri–Pila / the determinant method / dimension growth count split polynomials in a lattice
     coset? If not, say exactly why not.
   - Is ANY counting barrier possible for (ii)? (Is the number of members abundant, ≈ S·V/N? Can one prove that?)
     If counting cannot work, what is the right model for a barrier (a computational/query model), and what can be
     proved in it?
   - Inheritance: codex-type examples show members can share compatibility through one identity P − P' ≡ small.
     For k ≥ 2, can you show that an inherited compatibility is itself (or yields) relations of the target type?
3. **Construction attempt.** Use lattice reduction to FIND non-exact compressions: given N and F (or μ, m), search for
   split P with small roots (and multipliers Q) with QP ≡ small mod (N, F). Analyse quantitatively what residual
   size is achievable vs the needed one, and how many independent relations one lattice yields. In particular
   analyse "linearise the last s roots through e₁..e_s and solve a CVP" and meet-in-the-middle variants; give the
   cost per relation as a function of (N, B, V, s). Run small exact experiments if useful (report them as
   UNVERIFIED scratch results).
4. **Verdict**: barrier, construction, or a precise residual open statement. Is O1 a genuinely separate escape, or
   is it equivalent to (contained in / containing) the general problem of output-sensitive listing of
   split-by-construction integers with small residue mod N? Give your probability that O1 yields sub-L[1/3].

Be concrete; prefer short exact statements with proofs over surveys. Cite literature precisely (authors, venue,
year, theorem number) and flag anything you could not verify.
