# pilist round 1 (BLIND): the product-listing primitive behind several "below L[1/3]" open remainders

You are a research mathematician (analytic number theory, geometry of numbers, coding theory / list decoding,
lattice algorithms) collaborating on a speculative project about asymptotically faster classical integer
factorization. Rigour matters more than speed and honesty more than a positive result. Label every claim PROVED
(give the proof), KNOWN (cite precisely: authors, venue, year, theorem number), HEURISTIC (name the assumption),
model-PROVED (exact consequence of a named model), or OPEN. Do not read other files of this repository; work from
this prompt, your knowledge and web search. Scratch computation (Python, exact integer arithmetic) goes in the
current directory. Final answer: a structured report (≤ 6000 words) with numbered sections matching the tasks.

## Notation and background

x = log N, L[α, c] = exp((c+o(1)) x^α (log x)^{1−α}); L[c] means L[1/3, c]. N = pq a balanced semiprime.

**Heuristic H (NFS cost model).** Values of size exp(s x^{2/3}(log x)^{1/3}) are L[β]-smooth with probability
L[−s/(3β)]; joint smoothness events multiply; a smoothness test costs L^{o(1)}; sparse linear algebra on C columns
costs C^{2+o(1)}; dependencies are fruitful. Base-m NFS of degree d = δ(x/log x)^{1/3}, box |a|,|b| ≤ E = L[e],
rational factor base L[β], algebraic L[β_a], V = L[v] fields (MNFS). Rational penalty r = 1/(3βδ), algebraic
a = (δe + 1/δ)/(3β_a), columns s = max(β, β_a + v); relations need 2e + v − r − a ≥ s; LA 2s; box sieve 2e.
Known constants: GNFS 1.9230, MNFS (Coppersmith) 1.9019, factory per-N 1.6386.
"Listing" targets computed in that model (treat as given): output-sensitive listing of the rational-smooth pairs
gives 1.6386; partial listing at cost E²·L[−θr] gives c(θ) with initial slope ≈ −0.21 per unit θ (every fixed θ>0
helps; θ = 0, i.e. savings L^{o(1)}, gives nothing); output-sensitive algebraic-side listing gives infimum
(4/9)^{1/3} = 0.7631; output-sensitive doubly-smooth listing gives infimum 0 (and, with early stopping and
polylog factor bases, heuristic polynomial time — a quantified Schnorr/Adleman observation).

**The primitive (one name for several open problems).** "Π": cheaply LIST products of factor-base pieces that are
congruent to small values (or to other such products) modulo N, or split polynomials in an N-lattice, at cost below
one unit per output. Instances:
- (Π-nonexact) Factor base 𝔛 of ≤ L[β] units mod N computable from N (small primes; m − a for a computed m with
  |a| ≤ L[β]). Output L[β]^{1+o(1)} distinct multisets A with ∏_{a∈A} x_a ≡ v (mod N), where either
  0 < |v| ≤ e^{O(ℓ)} with ℓ = x^{2β}(log x)^{1−2β} (one-sided), or v itself split over 𝔛 (doubly split); u ≠ v as
  integers; F₂-rank of the outputs over the primes ≥ (1−o(1))·#columns; total time L[β]^{1+o(1)} including any
  presentation. Natural catalogues are NOT scarce: S distinct split products give ≥ S²/N − S congruent pairs, and
  under GRH the one-sided count is (1+o(1))S(2V+1)/N. So any obstruction must be about the cost of FINDING members.
- (Π-split) Output-sensitive enumeration of {P = ∏_{l≤M}(X − a_l), |a_l| ≤ B} ∩ {P(m) ≡ small (mod N)} at cost
  |output|·L^{o(1)}, output ≥ L[α], α < 1/3, with log|v| ≤ x^γ, γ < 2/3, or v smooth by identity.
- (Π-NFS) Output-sensitive or partial listing of rational-smooth / algebraic-smooth / doubly-smooth pairs (a, b)
  in the NFS box (targets above).

**The unproved crux.** A "lattice-query model": the algorithm may set up Λ_{M,R} = {(a, b) : a ≡ Rb mod M} for M
with known factorisation into factor-base primes (cost ≥ 1 even if empty), enumerate its box points (cost 1 each),
test points (L^{o(1)}). Under conditional cofactor randomness and the assumption H_sel ("for M > U = box area,
E[#new box points of Λ_M | history] ≤ C·U/M", i.e. nonempty large-modulus lattices cannot be identified for less
than one unit each), the gain over the box sieve is ≤ (C u log u)^{log U/log B} = L^{o(1)}. H_sel is ASSUMED, not
proved. Boneh's CRT list decoding (STOC 2000; JCSS 64 (2002), Thm 3.1/4.1: in poly time find every |t| < H with
gcd(f(t), P) ≥ S provided log H ≲ (log S)²/(D log P)) is the one known batch identifier; applied row by row it gains
only L^{o(1)}. Not excluded: whole-box multivariate decoders; adaptive special-q elimination where each relation has
Θ(u) "tail" primes in (L[η], L[β]] (u = log(value)/log B) — with o(u) tail primes it is provably no better than MNFS
with factor base L[η].

**PTE / pencils background.** P_A = ∏_{a∈A}(X − a), A a multiset of M integers in [−B, B]. Lemma R (proved):
if P, P′ are such, Q, Q′ ≠ 0 have degree ≤ j, G = QP − Q′P′ has degree < K and height ≤ h, and
M ≥ K + 2⌊log₂(2^{K−1}√K·h·B^{K−1})⌋ + 4j + 3, then G = 0. (Proof: remove the gcd; a Mahler-measure bound gives
|G̃(a)| ≤ Δ at integer points; at a root a of U, |U′(a)| = ∏|a − a′| ≤ Δ, so at most log₂Δ roots of U′ lie at distance
≥ 2 from a; vanishing moments of the difference measure force support ≥ M − K + 2; count supports.) Consequence:
one presentation F (any height) with exact members hP_A = FQ + r (small h, Q, r) yields one relation. For k = 1
fibres (members Q₀ − t, all split, roots in [−B, B], M ≥ 4) there are ≤ 2c·B^{1/2}(log B)^κ + 4 members (root-pair
curve (Q₀(x) − Q₀(y))/(x − y) = 0 plus the Binyamini–Cluckers–Novikov uniform bound c·d²·B^{1/d}(log B)^κ for
integral points of height ≤ B on an integral affine plane curve of degree d, arXiv 2302.04209 Thm 2).
Open item O2: linear systems of ≥ 2 unrelated presentations, members P_A = F₀ + ν_A F₁ + r_A (φ(F₀) = φ(F₁) = 0,
r_A small, deg r_A < k). Each ν gives ≤ 1 member by Lemma R. "Missing lemma: rigidity for small integer three-term
combinations c₀P₀ + c₁P₁ + c₂P₂ = G of split polynomials with G small and of low degree. Only the obvious pencils,
such as (X + ν)G, are seen to degenerate."

## Tasks

1. **Canonical statement.** State ONE canonical product-listing problem Π that subsumes Π-nonexact, Π-split and the
   Π-NFS listing primitives (say exactly how each is an instance, and which direction each containment goes). Give
   the exact per-output cost (as an L-exponent, with the window size and factor-base size as parameters) needed for:
   (a) beating MNFS 1.9019 at exponent 1/3; (b) any exponent α < 1/3; (c) exponent 1/4. Give formulas, not just
   words; check them numerically.

2. **Barrier attempts.** Try to PROVE the crux (H_sel, or a lower bound for listing smooth congruent products) in a
   natural model:
   - (a) a lattice-oracle / generic model; be explicit about why the model is or is not natural, given that all
     the rational-side roots are reductions of one integer m (R_q = m mod q);
   - (b) a "Coppersmith/LLL with bounded dimension" model, or more generally any decoder that certifies outputs as
     roots of small polynomials in the ideal (f, P)^k (Howgrave-Graham), including soft-decision / weighted CRT list
     decoding (Guruswami–Sahai–Sudan) and prime powers;
   - (c) a generic model of CRT/lattice manipulations.
   Compare with known lower bounds for lattice enumeration (e.g. pruned enumeration lower bounds), with the
   Coppersmith-bound tightness examples, and with hardness results for decoding Reed–Solomon codes beyond the
   Johnson radius (e.g. Cheng–Wan; Guruswami–Vardy; Gandikota–Ghazi–Grigorescu). Decide: is Boneh's L^{o(1)} a
   barrier (for which class) or an artifact (of what)? Quantify the expected number of useful outputs per decoder call.

3. **Construction attempts.** For each, compute the heuristic cost in the model above:
   - whole-box multivariate decoders (multivariate Coppersmith, Herrmann–May linear bounds, Jochemsz–May);
   - Guruswami–Sudan / Parvaresh–Vardy / folded or multiplicity-style list decoding of CRT (number-field) codes
     (Goldreich–Ron–Sudan; Boneh; Guruswami–Sahai–Sudan; any beyond-Johnson CRT decoder you know of) applied to
     finding many smooth values at once; what property would a beyond-Johnson decoder need, and would it suffice?
   - adaptive special-q elimination with Θ(u) tail primes: can adaptivity make the tail hypergraph close up with
     fewer relations than the number of tail primes, or make its elimination cheaper than sparse LA on L[β] columns?
     Try to prove a lower bound on the number of relations needed for a single dependency in a model where each
     relation has a bounded number of adaptively chosen ("forced") tail primes and Θ(u) random ones.

4. **PTE O2 pencils.** Attempt the three-term rigidity lemma generalising Lemma R, or find a counterexample family.
   In particular: is the naive statement ("c₀P₀ + c₁P₁ + c₂P₂ = G with small c_i, small low-degree G forces a
   degenerate pencil") true? If not, what is the right replacement (e.g. a count of split members of a pencil
   F₀ + tF₁, or of a linear system F₀ + νF₁ + (deg < k))? Try to extend the k = 1 fibre argument to general pencils.

5. **Experiments.** Propose 2–3 small exact experiments (with falsifiable predictions) that would test your claims
   at toy scale on a 16 GB machine (≤ 2 GB RSS). You may run small scratch checks yourself.

6. **Verdict.** Where exactly does the L[1/3] limit bite for Π? Probability (subjective) that Π_β is solvable for some
   β < 1/3.
