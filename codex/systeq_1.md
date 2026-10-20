# systeq round 1 (BLIND): an integer analogue of the systematic equation x^q − x = ∏(x − a)?

You are a research referee and collaborator on a speculative project about asymptotically faster
integer factorization. Be bold in exploring and ruthless in killing ideas. Honesty matters more than
a positive result. Do not read other files in this repository; work from this prompt, your own
knowledge and web search. Put any scratch computation in the current directory.

## Background (the wall)

Notation: x = log N, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}).

Smoothness-based factoring with candidates of size X and smoothness bound L[α] needs
log X ≲ x^{2α} (up to log factors) for the smoothness probability to be L[α]^{-O(1)}. NFS ties
candidates to N via a common root mod N of two integer polynomials f, g, so N | Res(f, g), and
Hadamard gives N ≤ ‖f‖^{deg g}‖g‖^{deg f}. With a search box of entropy C this forces
log(norm product) ≥ 2√(κ x C) and hence α ≥ 1/3. Every polynomial-selection variant we know (MNFS,
tower/TNFS, higher sieving dimension, multiple rational sides) obeys this ledger for general N.

In small-characteristic finite fields the wall was broken: Joux 2013 (L[1/4]),
Barbulescu–Gaudry–Joux–Thomé 2014 (quasi-polynomial), Granger–Kleinjung–Zumbrägel 2014/2018
(on the powers of 2 / rigorous for a family), Kleinjung–Wesolowski 2019 (arXiv 1906.10668,
rigorous quasi-polynomial for fixed characteristic). The engine is the systematic equation
x^q − x = ∏_{a∈F_q}(x − a), plus Frobenius and the subfield F_q ⊂ F_{q^k}: the representation
F_{q^{k}} = F_q[x]/(I) with h₁(x^q) ... (x^q ≡ h₀/h₁), the PGL₂(F_q) homographic action on
x^q − x giving q³ relations among linear factors at once, and degree-2 elimination (or KW's
"traps"/ Frobenius-twisted descent on an elliptic-curve-type representation).

## Tasks

1. **Axiomatise the descent engine.** State the minimal abstract properties that the small-char
   algorithms use (for BGJT/GKZ and separately for KW): e.g. (E1) a ring R with a cheap
   endomorphism φ whose fixed ring is a "small" subring k; (E2) an identity in R[X] that splits
   completely into linear factors over k, of degree comparable to |k|; (E3) a group action
   (PGL₂(k)) producing |k|³ (or more) relations from one identity; (E4) a representation where
   φ(X) equals a low-degree rational function of X, so the identity becomes a low-degree relation
   in the target field; (E5) a degree-reduction (elimination) step, and the size measure that it
   decreases. Be precise about which are essential for quasi-polynomial vs L[1/4].

2. **Integer analogues.** For each candidate below, decide whether it provides an analogue of each
   axiom, compute (or estimate) the number and size of smooth-by-construction numbers it produces
   at the scale needed, and — crucially — how those numbers are *related to N* (i.e. what
   congruence mod N they carry). Then test against the wall: does it yield relations among numbers
   of log-size x^{2α} with α < 1/3 carrying non-trivial information mod N, at density enough to
   produce ≥ (factor-base size) relations in L[α] time?
   - cyclotomic identities X^n − 1 = ∏_{d|n} Φ_d(X) evaluated at integers or at X ≡ ... mod N;
     Aurifeuillian factorizations;
   - products of consecutive integers / binomial coefficients (smooth by Kummer/Legendre);
     n! mod N (Pollard–Strassen; why are factorial methods exponential?);
   - Chebyshev/Dickson polynomials, D_m ∘ D_n = D_{mn} (a "Frobenius-iterate" analogue);
   - elliptic divisibility sequences / division polynomials; Lucas sequences;
   - norm-form identities in cyclotomic fields (Gauss/Jacobi sums, Stickelberger);
   - p-adic / Witt-vector lifts of X^q − X (Teichmüller: X^{p^k}… ), and why a ring map
     Z/p^k → Z/N is impossible;
   - group-ring identities; anything else you think of (e.g. the "field with one element" heuristic
     F_1-adic: X^n − 1 as the F_{1^n} systematic equation; Λ-rings / Adams operations ψ^n as integer
     Frobenius lifts; Borger's Λ-ring F_1).
   For the Λ-ring / Adams-operation idea specifically: Z[X] with ψ^p(X) = X^p is a Frobenius lift
   for every p simultaneously. Does the identity X^p ≡ X mod p (Fermat) give anything mod N?

3. **Most promising mechanism.** Develop the best candidate into an algorithm sketch with a
   heuristic cost exponent, or else state why each fails as a precise barrier of the form "any
   integer identity of type X produces ≤ ... relations mod N among numbers of size ≤ ...". A good
   generic barrier would be: a family of numbers smooth by construction is a family of values of
   products of small-height linear forms; a congruence mod N among them is a common-root condition
   and is subject to Res ≥ N. Is that right? Can the "Frobenius" (an endomorphism of Z/N that is
   not the identity and is cheaply computable) exist at all for Z/N? (Note: the only ring
   endomorphism of Z/N is the identity; but of (Z/N)[X]/(F) there are many — e.g. X ↦ X^N is not
   a ring endomorphism unless N is prime. What about the ring automorphisms of (Z/N)[ζ_m] given by
   ζ ↦ ζ^a? Those exist for every N. Do they give anything beyond the cyclotomic/Gauss-sum
   identities?)

4. **Kill test.** Propose one cheap, exact, toy experiment (N of 40–80 bits) that would measure
   whether the best candidate's relation yield beats the plain size-bound prediction (i.e. random
   integers of the same size). Give a precise predicted number.

Label every claim KNOWN (with a citation), PROVED (give the proof), HEURISTIC (name the assumption),
or SPECULATIVE. Give honest probabilities. Keep the answer under ~4000 words.
