# Adversarial referee round: classical Fourier-side period finding for a^x mod N

You are a hostile but fair referee (number theory + TCS). Below are claims from a research note. For EACH claim:
- verdict: CORRECT / CORRECT WITH FIXES / FALSE (counterexample) / OVERCLAIMED;
- a check of the proof, plus small-case computations where useful (you may run Python in your scratch dir);
- literature: is it KNOWN? Cite precisely, especially Akavia (2008 thesis; CRYPTO 2009; COLT 2010), Kushilevitz–Mansour 1993,
  Goldreich–Levin 1989, Gilbert–Muthukrishnan–Strauss, Aggarwal–Maurer 2009, Jager–Schwenk 2009, Cleve 2004,
  Van den Nest 2011, Schwarz–Van den Nest 2013, Aaronson "BQP and the PH" (Fourier Fishing), Donoho–Stark 1989,
  Tao 2005 uncertainty principle for Z/p, Granville–Soundararajan (pretentious characters), Bourgain–Glibichuk–Konyagin.
Be concrete. Report in ≤ 3500 words. Then add a section "strongest objection to the verdict".

Notation: N = pq odd balanced semiprime, G = (Z/N)^*, a ∈ G, H = ⟨a⟩, r = |H|, r_p = ord_p(a), r_q = ord_q(a),
g = gcd(r_p, r_q), i = [G:H] = φ(N)/r. e(x) = exp(2πix). For g: G → C, ĝ(χ) = E_{y∈G} g(y) χ̄(y). F_g(t) = g(a^t),
F̂(j) = E_{t∈Z/r} F_g(t) e(−jt/r). η_j is the character of H with η_j(a) = e(j/r).

**Verdict of the note (plain language).** There is no classical Fourier-side mechanism. Sparse Fourier learning on
g(a^x) works only if g correlates with a multiplicative character of large order. Every computable class we could
analyse is either provably flat, CRT-local (and then it factors N directly at least as fast), of bounded character
order, or (smoothness observables) limited by its density. Linear "interference" filters need ≥ r/2 evaluations to
isolate one peak.

**Claim T (transfer identity; PROVED).** F̂(j) = Σ_{χ ∈ Ĝ : χ|_H = η_j} ĝ(χ). Hence |F̂(j)|² ≤ i·Σ_{χ|_H=η_j}|ĝ(χ)|².
The mass on peaks of order ≥ d is ≤ i·(mass of ĝ on characters of order ≥ d). A τ-heavy peak (|F̂(j)|² ≥ τ, |g| ≤ 1)
of order d forces some χ of order ≥ d with |ĝ(χ)| ≥ √τ/i.
("Only characters concentrate": an exact identity, not a model statement.)

**Claim W (additive/Gauss class; PROVED).** If g(y) = b_0 + Σ_{h∈G} b_h e(hy/N), then for j ≠ 0,
|F̂(j)| ≤ (Σ_{h≠0}|b_h|)·√N/r. Same bound for g(y) = ((y+c)/N) (c ∈ G; Jacobi sums), with no exceptional frequency.
For the interval/MSB bit [y < N/2], |F̂(j)| = O(√N log N / r) for balanced p, q (the nonunit-h terms included).
Numerically the bound √N/r is attained to 3 digits as the max over j (r ~ 10^6–10^7).

**Claim HNP (why Akavia's HNP success does not transfer; explanation).** In the HNP the hidden map x ↦ αx is an
additive automorphism of Z/p, and MSB has O(log p) additive ℓ¹ mass, so MSB(αx) has heavy coefficients at αk, k small.
Here the hidden map t ↦ a^t lands multiplicatively, so SFT on t sees the multiplicative spectrum of g (Claim T). By
Pólya–Vinogradov/Gauss, MSB is multiplicatively flat. The observable must be concentrated in the dual of the group
where the hidden homomorphism lands.

**Claim GR-F (generic-ring Fourier barrier; PROVED, packaging no prior formulation found by us).** Let Ψ be computed
from y by any algebraic decision tree over Z/N: ring operations, constants from Z/N, branching on "v = 0 mod N",
any size, output a function of the leaf with |Ψ| ≤ 1. Extension rings (Z/N)[x]/(f), tori, elliptic-curve formulas are
included. Let ε = Pr_t[some test on the run on a^t is 0 mod exactly one of p, q], and say such a test "splits"; then
gcd(v, N) is a factor. Then Σ_{j : (r/g) ∤ j} |F̂(j)|² ≤ 8ε. The frequencies in (r/g)Z reveal only divisors of g | gcd(p−1, q−1).
Proof: off the split set S, the N-run equals the p-run and the q-run, so F = F_p + e_p = F_q + e_q with F_p r_p-periodic,
F_q r_q-periodic and ‖e‖² ≤ 4ε. F̂_p lives on (r/r_p)Z, F̂_q on (r/r_q)Z, and their intersection is (r/g)Z.
Corollary: sparse Fourier on a generic observable is never better than evaluating it at ~1/ε random orbit points and
taking gcds.

**Claim DB (density bound; PROVED, trivial).** If g = 1_E (centred), every nonzero relative peak satisfies
|F̂(j)|²/‖F − mean‖² ≤ ρ/(1−ρ), ρ = orbit density of E. So rare events (smoothness at small B) cannot give
1/poly-heavy peaks.

**Claim SM (smoothness observable; HEURISTIC + data).** For g = [y is B-smooth] (integer rep in [1, N)), a
pretentious character χ with χ(ℓ) ≈ 1 for all ℓ ≤ B exists (pigeonhole) with max_ℓ ‖arg χ(ℓ)/2π‖ ≤ φ^{−1/π(B)}.
It makes peaks heavy, up to the density bound ρ, once φ^{−1/π(B)}·Ω(y) ≪ 1. Data (B = 50): the max peak is 40–340×
the random baseline (ln r)/r at n = 22–26 bits, ≈ 0.3–0.5·ρ. For SFT one needs ρ ≥ 1/poly, so u ≲ log n/log log n,
so π(B) is huge, so no coherent χ exists (heuristic, equidistribution of (χ(ℓ))_ℓ). Hence "dual index calculus" is
dominated by ordinary (primal) relation collection. Framing: classical relation methods sample the hidden lattice L
(primal); Shor samples L* (dual).

**Claim UP (linear interference; PROVED).** Let w be supported on m points of Z/r. Suppose that for every F with
prescribed flat magnitudes |F̂(j)| = c (j ≠ j0, worst-case phases) the filter output Σ_x w(x)F(x) is dominated by the
j0 term, i.e. |ŵ(j0)| > Σ_{j≠j0}|ŵ(j)| (unnormalised ŵ). Then m > r/2.
Proof: |ŵ(j0)| ≤ ‖w‖₁ ≤ √m‖w‖₂ and ‖ŵ‖₁ ≥ ‖ŵ‖₂²/‖ŵ‖_∞ ≥ r‖w‖₂/√m. So 2√m‖w‖₂ ≥ r‖w‖₂/√m.
Interpretation: building interference by evaluating at strategic x is a linear filter; isolating one period peak costs
Θ(r), versus √r for birthday (a quadratic statistic). Data: estimating the strongest peak of e(ha^t/N) from m random
samples needs m ≈ r/ln r.

**Claim BASIS (compressed representation; PROVED).** For any orthonormal basis {β_k} of the second register, Shor's
state is Σ_k √w_k |φ_k⟩|β_k⟩ with φ_k ∝ F_{β̄_k}. The classical strategy "sample k ∝ w_k, then SFT φ_k" succeeds iff
the branches are Fourier-concentrated, i.e. (Claim T) iff the β_k correlate with characters. The Schmidt rank across
x|y is r in every basis. In the character basis the state is Σ_j |tone_j⟩|η_j⟩, a superposition of r x-product states
(maximally compressed), but the basis change has entries χ(y), i.e. discrete logs. So "compression in a
multiplicatively adapted basis" is equivalent to the concentrating-g question.

**Claim Q1 (partial information).**
(a) For Blum N, deciding v₂(r) ⟺ quadratic residuosity (poly-time equivalent).
(b) An oracle for [3 | ord_N(a)] on random a decides, for N ≡ 1 mod 3, whether p ≡ q ≡ 1 mod 3, i.e. whether −3 is
a square mod N (the fixed-base QR bit; no reduction from factoring known).
(c) a^{N+1} = a^{p+q}, giving an interval DL of width O(√N), i.e. N^{1/4}; iterating gives a^{N^k} = a^{(s−1)^k}.
Cheon-type speedups need a known divisor of the exponent-group order (Z/r)^*, which is hidden.
(d) The only reciprocity-computable characters of G are generated by the Jacobi symbol. Power-residue symbols of
rational integers at (N) are trivial, and order-ℓ symbols need an element π ∈ Z[ζ_ℓ] of norm N. Every computable
character has poly(n)-smooth order, so it reveals divisibility of r by small primes only.

**Claim Q3 (special a).** Special elements (small a, Lucas, tori, cyclotomic rings, a^{N^e} cycling) fall under
order annihilation: smooth order or a {p,q} coincidence (e.g. r_p | N^e − 1 ⟺ q^e ≡ 1 mod r_p). The formal
1 + √N in (Z/N)[t]/(t²) has known order N.

**Claim PF (barrier statement).** Combining T, GR-F, W, DB and UP: any classical algorithm that learns r by
sparse-Fourier/heavy-hitter queries to g(a^x), for g CRT-local (any size), additive/Weil-class, Jacobi-generated,
or a rare event of density ρ, finds no nontrivial τ-heavy peak unless τ ≤ max(8ε_split, B²N/r², ρ/(1−ρ)). In the
CRT-local case the algorithm factors N by gcd in O(1/ε_split) evaluations anyway. Every heavy-hitter algorithm with
an output-list guarantee costs ≥ 1/τ in these regimes.

Attack hardest: GR-F (is the "(r/g)Z reveals only divisors of g" reading right? adaptivity? randomised trees?
outputs that depend on computed values rather than just the leaf?), UP (is the worst-case-phase hypothesis
meaningful?), SM, and Q1(d).
