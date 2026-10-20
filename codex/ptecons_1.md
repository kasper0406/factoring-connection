# ptecons round 1 (BLIND): CONSTRUCT dense N-dependent split families, or show why each structured attempt fails

You are a research collaborator on a speculative project about asymptotically faster integer
factorization. Be inventive first, then ruthless with every candidate. Honesty matters more than a
positive result. Do not read other files in this repository; work from this prompt, your own
knowledge and web search. Put any scratch computation (Python is fine; sympy is available) in the
current directory. Keep each computation small (< 2 GB RAM, < 10 min).

## Background

Notation: x = log N, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}).

Smoothness-based factoring needs relations u ≡ v (mod N), u ≠ v, both sides B-smooth. NFS ties
candidates to N via a common root of two integer polynomials, so N | Res(f, g); this forces α ≥ 1/3.

In small-characteristic finite fields the wall was broken (Joux 2013; Barbulescu–Gaudry–Joux–Thomé
2014; Granger–Kleinjung–Zumbrägel; Kleinjung–Wesolowski 2019). Engine: the systematic equation
X^q − X = ∏_{a∈F_q}(X − a), plus a presentation F_{q^{2k}} = F_{q²}[X]/(I) with X^q ≡ h₀/h₁
(deg h_i ≤ 2). Because Frobenius is additive, every member of the PGL₂(F_{q²})-orbit of X^q − X —
(aX+b)^q(cX+d) − (aX+b)(cX+d)^q, about q³ + q distinct squarefree fully split polynomials — reduces
modulo I to a polynomial of degree ≤ 3. One fixed congruence compresses a dense family of split
(= smooth by construction) objects to small objects.

A previous agent proved these integer obstructions (take them as given):
- (F′) restricted Frobenius compatibility, (X + b)^s ∈ span(1, X^s) for s + 1 constants b with unit
  differences in a rank-≥4 Z/N-algebra, is impossible when all primes of N exceed s;
- (CF) a systematic product ∏_{c∈S}(z − c) vanishing on all scalars of a Z/N-algebra needs |S| ≥ p_min;
- (DES) by Descartes' rule, V_n = span(X^{n+1}, X^n, X, 1) has no nondegenerate member split over R
  for n ≥ 4;
- (EX) in a free presentation A = (Z/N)[X]/(F), a congruence u ≡ v mod (N, F̃) whose difference
  reduced mod F̃ has coefficients < N is an exact identity over Z[X]/(F̃) and carries no information
  about N; imposing a non-exact compression costs d conditions mod N (heuristically factor-base norms
  ≥ N^{d/(d+1)}, worse than ordinary polynomial selection).
- If U ≡ V (mod N), U ≠ V, then max(|U|, |V|) ≥ N/2; products of small integers that are below N and
  congruent mod N are equal.
- The F_1 analogue (μ_n ∪ {0}) has |PGL₂| = 2n: orbits of size O(1), not q³.

**The open target.** A family 𝓕 of pairs (u, v), u ≡ v (mod N), u ≠ v, with
- (a) u a product of ≳ log N / log B pieces, each ≤ B^{O(1)} (or of norm ≤ B^{O(1)} in a fixed small
  number field), by an identity;
- (b) v small: log-size o(x^{2/3}), ideally independent of x;
- (c) |𝓕| ≥ L[α] with α < 1/3 (or at least ≥ the factor-base size at a B giving sub-L[1/3]), and 𝓕
  computable without knowing the factors, at cost ≈ |𝓕|·L[o(1)].

Concretely: dense families of split polynomials ∏_j (X − a_j) with small a_j (in Z, or in a fixed
small ring of integers) sharing a common reduction to small polynomials modulo an N-dependent
presentation (N, F). Example (exact over Z, not N-dependent): the Dickson identity
D₆(X) − 2 = (X² − 4)(X² − 1)² = X⁶ − 6X⁴ + 9X² − 4 ≡ X + 3 (mod X⁶ − 6X⁴ + 9X² − X − 7).

## Tasks — CONSTRUCTION side

For each candidate below: write the family explicitly; compute (i) its density (number of members
with entries ≤ A, as a function of A, degree M, field degree), (ii) the residual size |v|, (iii) how
N enters (which N, which presentation, at what cost), and (iv) check (a)–(c) and each barrier F′, CF,
DES, EX. Be explicit about whether the many relations of a family are *genuinely distinct* (distinct
exponent vectors carrying new information mod N) or collapse to one relation (e.g. Galois/automorphism
images of one relation, or differences that are exact integer identities).

1. **Dickson/Chebyshev compositions.** D_n(x + 1/x) = x^n + x^{−n}, D_m∘D_n = D_{mn}. D_n(X) − D_n(t)
   splits over Q(ζ_n + ζ_n^{−1}, τ) with τ + 1/τ = t. When do many shifts D_n(X) − c split over a
   small number field with small roots, and can they share a common reduction modulo some F depending
   on N? Generalise the D₆ example; classify integer t for which D_n(X) − D_n(t) splits over Z.
2. **Prouhet–Tarry–Escott.** Frolov/Tarry-type parametric PTE families, ideal solutions of size ≤ 12,
   ideal chains, Prouhet–Thue–Morse. A PTE class {multisets with equal e₁..e_{M−t}} gives split
   polynomials differing by degree < t. Can such classes be made to reduce compatibly modulo an
   N-dependent F, or aligned so that their common top part T satisfies T(m) ≡ small (mod N) for a
   root m of F? What is the cost of that alignment? How long can a class (chain) be, as a function of
   M, t and the entry bound A, for structured constructions vs. generic counting?
3. **Cyclotomic and Galois orbits.** ∏(X − ζ^a) over Z[ζ_n]; orbits under (Z/n)^*, Galois groups,
   or the normaliser of μ_n acting on roots — the characteristic-0 analogue of PGL₂(F_q)-orbits.
   Use N with a convenient reduction (N | Φ_n(m) Cunningham-type; N with a known small-norm element;
   N = norm of something). Gauss/Jacobi-sum identities (Stickelberger). Do orbits give distinct
   relations, or Galois conjugates of one relation?
4. **Elliptic analogues.** Division polynomials, multiplication-by-m, Lattès maps; preimage families
   split over torsion fields. Does anything replace additivity of Frobenius?
5. **Group-ring / norm constructions.** Norms from Z[G], symmetric-group actions on root multisets,
   resultants of split polynomials.
6. **Lattice/LLL families.** Families where N enters through a lattice: e.g. polynomials whose
   coefficient vectors lie in a lattice built from N (polynomials vanishing at m mod N), intersected
   with the variety of split polynomials; Schnorr/Adleman-type relation lattices. What is the expected
   number of short lattice points that are split polynomials, versus random?
7. **Your own ideas.** Anything that supplies "one N-dependent alignment, many relations".

Then:
- Pick the best one or two candidates. Give an algorithm sketch and a heuristic L-cost analysis
  (random-integer smoothness ρ(u) ≈ u^{−u}, linear algebra quadratic in factor-base size). State
  precisely which parameter would have to be achievable for sub-L[1/3], and whether it is.
- Name the single strongest obstruction you see to ANY construction of this type, and say whether it
  is provable (cite tools: Vinogradov MVT / Bourgain–Demeter–Guth, determinant method, Bilu–Tichy /
  Siegel for polynomials with many integer-root shifts, Newton identities).
- Suggest one small exact experiment (≤ 10 min Python) that would measure the yield of the best
  candidate against a random control. If cheap, run it and report numbers.

Label every claim PROVED / KNOWN (with citation) / HEURISTIC / CONJECTURED / UNVERIFIED. End with a
verdict paragraph and a probability estimate that some construction of this type beats L[1/3].
