# prepread_1: blind round. High-degree algebraic preprocessing, characters, and dense short SLPs

You are a research mathematician (analytic number theory, arithmetic dynamics, Galois theory, algebraic
complexity) with web search. The repository at /Volumes/git/factoring-connection is read-only context (an
AI-assisted research project on CLASSICAL integer factoring). Background, if wanted:
notes/periodfind_results.md (§1.1 Theorem T, §1.3 Prop W and GR-F, §6 open item O1),
notes/slpcirc_results.md (§1.3 model U and W_U, §4 Theorem SP, §6 open question),
notes/nonsmooth_results.md §5.1 (Prop OA), §5.2 (Prop FROB). Cite them; do not re-derive them. Be bold
generating, ruthless killing; honesty beats a positive result. Label every claim PROVED / KNOWN (exact
reference: author, year, venue, theorem number) / HEURISTIC / CONJECTURED / UNVERIFIED.

## Setting

N = pq balanced, G = (Z/N)^*, a ∈ G random, r = ord(a). Theorem T (periodfind): the Fourier
coefficients of F(t) = g(a^t) over Z/r are sums of the multiplicative-character coefficients ĝ(χ) over the
χ restricting to a given character of ⟨a⟩; a 1/poly-heavy peak of order ≥ r/poly needs g to correlate
(≥ 1/poly) with a multiplicative character of large order. Additive observables e(h y/N) are flat
(Korobov / Gauss sums), polynomial phases e(f(y)/N) of degree D are flat with |F̂(j)| ≤ D²√N/r (Weil).

Model U (slpcirc): an N-independent constant-free integer SLP F_X of size s, evaluated mod N at a random
input, then gcd with N. It factors all balanced N with prime factors in (X, 2X] in poly time if s =
polylog X and the root density δ_p(F_X) ∈ [1/polylog, 1 − 1/polylog] for (almost) every p ∈ (X, 2X].
Theorem SP: if F is t-sparse in a torus coordinate, density ≥ λ forces every large prime factor of the
group order to divide the exponent-difference product (group-order method in disguise). Open: F dense in
every torus coordinate. Proving that no such F exists implies constant-free permanent lower bounds
(Bürgisser 2009), so only restricted barriers are possible.

## Tasks

1. **O1: preprocessing + additive readout.** Ψ(y) = e(h·Φ(y)/N) with Φ ∈ Z[y] of degree D computed by an
   SLP of size s (so D ≤ 2^s, possibly D ≫ p). By CRT and Theorem T the question per prime is: how large
   can the mixed sum S(χ, Φ) = Σ_{y∈F_p^*} χ(y) e(Φ(y)/p) be, relative to p, for χ of large order?
   (a) Weil/Deligne/Castro–Moreno: |S| ≤ (D + something)√p for Φ not of the form c + (Artin–Schreier
   degenerate) — state the exact bound and hypotheses (also for rational Φ). At which D does it become
   trivial (D ≈ √p)? Note Φ mod p only matters as a function on F_p, so reduce Φ modulo y^p − y: what is
   the right "effective degree"?
   (b) Can a polynomial with a short SLP (size polylog p) correlate ≥ 1/polylog with a multiplicative
   character of order ≥ p^{Ω(1)}, i.e. |S(χ, Φ)| ≥ p/polylog? Constructions or obstructions. Example to
   analyse: Φ = y^m (Gauss periods / Korobov), Φ with e(Φ(y)/p) ≈ χ(y) by design (how short can an SLP be
   whose additive phase imitates a multiplicative character on F_p^*? A discrete-log interpolation
   polynomial has degree p − 2 and is believed to need large circuits — what is known: Shparlinski,
   Coppersmith–Shparlinski "on polynomial approximation of the discrete logarithm", Lange–Winterhof,
   Meidl–Winterhof lower bounds on degree/sparsity/circuit size of DL-interpolating polynomials?).
   Note the uniformity issue: Φ must work for all p in the range, with a single integer program. Is there
   any barrier of the shape "SLP size s ⇒ |S(χ, Φ)| ≤ f(s)·p^{1−c} for ord χ ≥ R", perhaps conditional, or
   is that as hard as circuit lower bounds for DL? Is there a reduction showing a large correlation would
   itself give a DL/order algorithm (so O1 is equivalent to a known hard problem)?
2. **Dense short SLPs and Galois theory.** For F ∈ Z[x] squarefree and p ∤ disc F, #roots mod p =
   #fixed points of Frob_p on the roots, so Chebotarev gives density E[fix] = #orbits-of-Gal = number of
   irreducible factors — but only "on average over p" and only once p ≫ (disc, degree), whereas model U
   needs deg F ≳ δp (p < deg F, Chebotarev-silent regime). (a) Make precise what Galois theory / effective
   Chebotarev can and cannot say when deg F ≥ p. (b) Iterates of x² + c: arboreal Galois representations
   (Odoni 1985, Stoll 1992, Jones 2008/2013/2014 survey, Juul–Kurlberg–Madhu–Tucker 2016 on proportion of
   periodic points / Galois images; Boston–Jones; Ferraguti) give the fraction of p for which f^{∘n} has a
   root mod p → 0 (fixed-point proportion of iterated wreath products → 0). State exact theorems. Do these
   give, or nearly give, an inverse theorem for quadratic compositions (x² + c_k)∘⋯∘(x² + c₁) with fixed
   c_i (root density mod p at polylog depth k ≈ log p)? What happens in the Chebotarev-silent regime
   k ≥ log₂ p, where Pollard-rho/functional-graph heuristics (Flynn–Garton, Bach, Burnette–Schmutz, random
   mapping models) apply instead? (c) Which finite permutation groups (degree D, transitive or with few
   orbits) have fixed-point proportion (fraction of elements with a fixed point / mean number of fixed
   points) bounded in [1/polylog, 1 − 1/polylog] AND arise as Galois groups of polynomials with polylog
   SLPs of exponential degree? Abelian (cyclotomic/CM: x^M − 1, division polynomials) = group-order
   methods. Are there non-abelian candidates (e.g. GL₂(Z/m) on torsion of a non-CM curve is "group order"
   too; wreath products from dynamics have proportion → 0; Fried's exceptional/Chebyshev; monodromy of
   compositions; derangement theory: Fulman–Guralnick on derangements in classical groups; Łuczak–Pyber)?
   Is "a Frobenius element with many fixed points uniformly in p" only possible via a group-order
   mechanism? Formulate the cleanest conjecture/barrier you can, and prove whatever restricted version you
   can (e.g. for compositions of bounded-degree maps; for polynomials whose Galois group is an iterated
   wreath product).
3. **Kill or keep.** Give the most promising concrete N-independent polylog-SLP family that is dense in
   every torus coordinate and not a group-order method, with its predicted root density mod p, and the
   strongest argument against it. Suggest 1–2 small exact experiments that would discriminate.

Deliver a structured report (≤ 4000 words) with the labelled claims, exact references, and proofs for
anything PROVED. Write scratch files only under your working directory.
