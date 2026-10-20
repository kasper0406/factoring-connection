# ptebarrier round 1 (BLIND): prove a counting barrier for dense split families (multigrade Prouhet–Tarry–Escott with N-dependence)

You are a research mathematician (analytic number theory / Diophantine geometry) acting as a collaborator on a
speculative project about asymptotically faster integer factorization. Your job here is the BARRIER side:
prove upper bounds, rigorously, and state exactly where a proof stops. Honesty matters more than a positive
result. Do not read other files in this repository; work from this prompt, your knowledge and web search. Put
scratch computation (Python is fine) in the current directory.

## Background

Notation: x = log N, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}); L[α] means L[α, c] for some c > 0.

Smoothness-based factoring needs relations u ≡ v (mod N), u ≠ v, both sides B-smooth, B = L[α]. NFS ties
candidates to N via a common root mod N of two integer polynomials, which forces α ≥ 1/3 (product of the two
norms ≥ N, and the search box must be ≥ L[1/3] to get enough relations).

In small-characteristic finite fields the L[1/3] wall was broken (Joux 2013; Barbulescu–Gaudry–Joux–Thomé 2014;
Granger–Kleinjung–Zumbrägel; Kleinjung–Wesolowski 2019) by the systematic equation
X^q − X = ∏_{a∈F_q}(X − a) together with a presentation X^q ≡ h₀(X)/h₁(X) (deg h_i ≤ 2): the
PGL₂(F_q)-orbit of X^q − X is ≈ q³ fully split polynomials of degree q+1 which all reduce, modulo the
presentation, to polynomials of degree ≤ 3. One fixed congruence compresses a dense family of
smooth-by-construction objects to small ones.

A previous agent proved the following integer obstructions (take them as given):
- additive Frobenius-type maps on Z/N-algebras force λ(N) | s − 1 (so exhibiting one factors N), and a
  restricted version (X+b)^s ∈ span(1, X^s) for s+1 constants is impossible in rank ≥ 4;
- a product ∏_{c∈S}(z − c) vanishing on all scalars needs |S| ≥ p;
- (Descartes) span(X^{n+1}, X^n, X, 1) has no nondegenerate member split over R for n ≥ 4;
- (Lemma EX) in a presentation (Z/N)[X]/(F) with monic lift F̃ ∈ Z[X], a congruence u ≡ v mod (N, F̃) whose
  difference reduced mod F̃ has all coefficients of absolute value < N is an EXACT identity in Z[X]/(F̃); exact
  identities in a characteristic-0 ring carry no information about N (the relation lattice they generate is
  that of the ring itself, independent of the prime factors of N);
- if integers U ≡ V (mod N), U ≠ V, then max(|U|, |V|) ≥ N/2.

**The open escape (stated by that agent).** A family 𝓕 of pairs (u, v), u ≡ v (mod N), u ≠ v, with
(a) u a product of ≳ log N / log B pieces, each ≤ B^{O(1)}, by an identity (smooth by construction);
(b) v small: log-size o(x^{2/3}), or independent of x;
(c) |𝓕| ≥ L[α], α < 1/3, computable without knowing the factors of N.
Concretely: many split polynomials P(X) = ∏_{j=1}^M (X − a_j), a_j small integers, that share a common
reduction to a small polynomial modulo some presentation (N, F). If all P in the family agree in e₁..e_{M−k}
(elementary symmetric functions of the roots), then P − P′ has degree < k. Over Z this is a multigrade PTE-type
system: many integer multisets of size M, entries in [−B, B], with prescribed e₁..e_{s}, s = M − k
(equivalently, by Newton, prescribed power sums p₁..p_s), with N-dependence on top. Example of one member:
(X² − 4)(X² − 1)² = X⁶ − 6X⁴ + 9X² − 4 ≡ X + 3 modulo X⁶ − 6X⁴ + 9X² − X − 7 (exact over Z).

## Tasks

1. **Model.** Write down the most general precise model you can for this escape in which a theorem can be
   proved: presentation data (N; an integer or residue m, or a monic F of degree d and height H_F with a root
   mod N; or a number-field order), the family (size M, entry bound B, number k of free low coefficients, height
   bound H_v of the small side), and what counts as a relation and as *informative* (not an exact identity,
   cf. EX). Then state which parameter ranges would give factoring below L[1/3] in the usual L-cost model
   (smoothness probability u^{−u}, LA quadratic in the factor base). Is it enough to bound the number of
   multisets with prescribed e₁..e_s, or must the small side's height also be prescribed (note: the low
   coefficients e_{M−k+1}, …, e_M are the LARGE ones, up to ~(MB)^M)?

2. **Heuristic count.** Expected number of multisets of size M with entries in [−B, B] and prescribed
   e₁..e_s, as a function of M, s, B (and the height window for the remaining coefficients). Find the
   threshold where it drops below 1, and compare to what (c) needs.

3. **Vinogradov / PTE.** Let R(c) = #{a ∈ [1,B]^M : Σ a_i^j = c_j, 1 ≤ j ≤ s}. Then Σ_c R(c)² = J_{M,s}(B),
   the Vinogradov mean value. Work out precisely what Bourgain–Demeter–Guth (and Wooley) give for max_c R(c)
   and for the number of values c with large R(c) (Markov / Cauchy–Schwarz / higher moments), and — crucially —
   whether the implied constants are uniform when M and s grow with N (here M ≈ x / log B → ∞). If they are not,
   find what IS available uniformly (elementary bounds, Newton-identity injection, explicit Vinogradov
   constants, large-sieve / polynomial-method bounds). Is there an elementary bound of the form
   max_c #fiber ≤ (C B)^{M−s} · (something decaying in s) that is uniform?

4. **Determinant method.** Count integer points of height ≤ B on V_c = {e₁ = c₁, …, e_s = c_s} ⊂ A^M
   (dimension M − s, degree s!) via Bombieri–Pila, Heath-Brown, Salberger (dimension growth; Castryck–Cluckers–
   Dittmann–Nguyen uniform versions with polynomial dependence on the degree). Which is strongest in this range,
   and are any of them uniform in M, s?

5. **N-dependence.** The reduction is modulo a presentation depending on N with N ≫ B. Does the congruence
   version make solutions rarer, or can it be satisfied for free by small exact identities? Combine with EX:
   show (or refute) that any family whose members' congruences all follow from one fixed small-height
   presentation is exact (uninformative), and quantify how many independent N-dependent conditions each
   informative member must pay, as a function of the heights.

6. **Theorem.** State and PROVE the strongest clean theorem you can of the form: "in model X, the number of
   informative relations produced by any dense split family with parameters in the range that would beat
   L[1/3] is at most ..., hence route-1 families of type X cannot give α < 1/3". If a full proof is impossible,
   prove the strongest partial result and state precisely the missing step. Be explicit about which
   statements are counting (existence) bounds and which are cost (finding) bounds: a counting bound cannot
   kill an escape whose objects exist in abundance but are hard to find.

7. **Literature.** Is any such bound already in the PTE / Vinogradov / Tarry–Escott literature (Hua, Wright,
   Borwein–Ingalls, Alpers–Tijdeman, Caley, Wooley "Vinogradov systems", Prendiville, Brandes–Hughes, counting
   PTE solutions of size M via determinant method)? Give citations you are confident of.

Label every claim KNOWN (citation), PROVED (give the proof), HEURISTIC (name the assumption) or SPECULATIVE.
Give an honest probability that this escape yields sub-L[1/3]. Keep it under ~5000 words.
