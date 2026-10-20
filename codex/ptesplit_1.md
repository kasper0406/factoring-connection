# ptesplit round 1 (BLIND): dense N-dependent split families (multigrade Prouhet–Tarry–Escott) — construction or barrier?

You are a research referee and collaborator on a speculative project about asymptotically faster
integer factorization. Be bold in exploring and ruthless in killing ideas. Honesty matters more than
a positive result. Do not read other files in this repository; work from this prompt, your own
knowledge and web search. Put any scratch computation (Python is fine) in the current directory.

## Background

Notation: x = log N, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}).

Smoothness-based factoring needs relations u ≡ v (mod N), u ≠ v, both sides B-smooth with B = L[α].
Random candidates of log-size X are B-smooth with probability u^{−u}, u = X/log B. NFS ties
candidates to N through a common root mod N of two integer polynomials, so N | Res(f, g); this gives
log(norm product) ≥ 2√(2 x log E) for search box E and hence α ≥ 1/3.

In small-characteristic finite fields the wall was broken (Joux 2013, BGJT 2014, GKZ, Kleinjung–
Wesolowski 2019): the systematic equation X^q − X = ∏_{a∈F_q}(X − a) with a presentation where
X^q ≡ h₀(X)/h₁(X) (deg h_i ≤ 2) makes the PGL₂(F_q)-orbit of X^q − X (≈ q³ fully split
polynomials) all reduce, modulo the presentation, to polynomials of degree ≤ 3. Thus one fixed
congruence compresses a dense family of split (= smooth by construction) objects to small ones.

A previous agent ("systeq") proved these integer obstructions:
- (F′) restricted Frobenius compatibility (X+b)^s ∈ span(1, X^s) in a rank ≥ 4 Z/N-algebra for s+1
  constants with unit differences is impossible when all primes of N exceed s;
- (CF) a systematic product ∏_{c∈S}(z − c) vanishing on all scalars needs |S| ≥ p;
- (DES) by Descartes' rule, V_n = span(X^{n+1}, X^n, X, 1) has no nondegenerate member split over R
  for n ≥ 4;
- (EX) in a free presentation A = (Z/N)[X]/(F), any congruence u ≡ v mod (N, F̃) whose difference
  reduced mod F̃ has coefficients < N is an exact identity over Z[X]/(F̃) and carries no
  information about N; imposing a non-exact compression costs d conditions mod N (heuristically
  factor-base norms ≥ N^{d/(d+1)}).
- Barrier for the trivial version: products of small integers that are below N and congruent mod N
  are equal; if U ≡ V (mod N), U ≠ V, then max(|U|, |V|) ≥ N/2.

**The open target.** A family 𝓕 of pairs (u, v) with u ≡ v (mod N), u ≠ v, such that
- (a) u is a product of ≳ log N / log B pieces, each ≤ B^{O(1)}, by an identity (smooth by construction);
- (b) v is small: log-size o(x^{2/3}), or even independent of x;
- (c) |𝓕| ≥ L[α] with α < 1/3, and the family is computable without knowing the factors.

Concretely: many split polynomials P(X) = ∏_{j=1}^M (X − a_j), a_j small integers, that share a common
reduction to a small (low-degree, small-height) polynomial modulo some presentation (N, F) — e.g.
P(m) ≡ r(m) (mod N) for an integer m with F(m) ≡ 0 (mod N). If F has small degree d, and all P in the
family agree in their top M − k elementary symmetric functions e₁..e_{M−k}, then the differences
P − P′ have degree < k. Over Z this is a multigrade Prouhet–Tarry–Escott (PTE) problem. Example:
the Dickson-type identity (X² − 4)(X² − 1)² = X⁶ − 6X⁴ + 9X² − 4, which is ≡ X + 3 modulo
X⁶ − 6X⁴ + 9X² − X − 7 — a split polynomial congruent to a small one by an identity.

## Tasks

1. **Formalise exactly.** Give a precise formulation allowing general presentations (N, F) (F of any
   degree d, possibly several variables / a number field), the number k of free low coefficients, the
   multiset size M, entry bound B, and the residual side size. Using the standard L-cost model
   (smoothness probability u^{−u}, linear algebra quadratic in factor-base size), state what family
   size and parameters would beat L[1/3]; L[1/4]; L[1/10]. Is condition (a)–(c) correctly stated, or
   does the residual side v need to be smooth as well (a "doubly smooth" family), and how does the
   value P(m) at a root m mod N enter?

2. **Counting / barrier.** Multisets of size M with entries ≤ B number ~ B^M / M!. Prescribing
   e₁..e_{M−k} (with e_i up to size ~ (MB)^i) costs heuristically ∏ (MB)^i. Find when dense families
   are generically impossible. Then try to PROVE an upper bound on the number of multisets
   (a₁..a_M) ∈ [−B, B]^M sharing prescribed e₁..e_{M−k} (equivalently power sums p₁..p_{M−k}).
   Candidate tools: Vinogradov's mean value theorem (Bourgain–Demeter–Guth 2016; Wooley's efficient
   congruencing, cubic case 2016), Bombieri–Pila / Heath-Brown determinant method / Salberger
   dimension growth. Determine precisely whether the Vinogradov main conjecture (now a theorem)
   bounds the number of PTE-type solutions in the needed range, and whether this yields a clean
   barrier: e.g. "for every presentation of degree d and height H, the number of B-smooth-by-
   construction split polynomials with a common reduction of degree < k is ≤ ...". Include the
   N-dependence: does requiring the common reduction modulo an N-dependent presentation cut the count
   further, or does it reduce to an exact identity over Z (as in EX)? Distinguish a family where all
   members are congruent to the SAME small v from one where each member is congruent to its own small
   v_P (which is what a relation needs).

3. **Constructions.** Look for structured parametric families beating the generic count: Frolov /
   Tarry–Escott symmetric constructions, Prouhet–Thue–Morse, Chebyshev/Dickson compositions
   D_m ∘ D_n (e.g. D_n(X) − c split for special c), norm forms, cyclotomic and Galois-orbit families
   (the closest analogue of PGL₂ orbits), symmetric-group actions, Borwein's results. For each,
   compute the density, the residual size and check (a)–(c). In particular: the family
   {D_n(X) − D_n(t)} = ∏ (X − (ζ^j τ + ζ^{−j} τ^{−1})) splits over a cyclotomic/real field, not Z; and
   the polynomials whose values at m mod N are ≡ small. Does any integer family reach the
   F_q density q³ at "degree" q?

4. **Kill test.** Propose a small exact experiment (M ≤ 10, entries ≤ 30–60) that counts multisets
   sharing e₁..e_{M−k} (or sharing a reduction modulo a random presentation of N) and compares to the
   heuristic; give predicted numbers.

5. **Literature.** PTE (Borwein, Ingalls, Alpers–Tijdeman, Caley, Wright's conjecture, Hua), PTE and
   Vinogradov systems, Bombieri–Pila counts on such varieties, any use of PTE in factoring or
   cryptanalysis.

Label every claim KNOWN (citation), PROVED (give the proof), HEURISTIC (name the assumption) or
SPECULATIVE. Give honest probabilities that this route yields sub-L[1/3]. Keep it under ~4500 words.
