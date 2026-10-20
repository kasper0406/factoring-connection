# analytic_2: ADVERSARIAL referee round. Break these claims

You are a hostile referee: an expert in analytic number theory, geometry of numbers and complexity theory, with
web search. The repository /Volumes/git/factoring-connection is read-only context. Relevant files:
notes/analytic_prereg.md, src/analytic/{hiary,hull,cubic,ladder}.py, your own earlier blind report
codex/scratch/analytic_1/analytic_1_report.md (you may disagree with it), notes/nonsmooth_results.md.
Only classical algorithms matter.

For EACH claim: (a) is it true as stated? Give a counterexample or the exact gap if not. (b) Is it already
in the literature? Give exact references. (c) Is the label (PROVED / HEURISTIC / model statement) right?
(d) Give the corrected statement. Be concrete, check the algebra yourself, and run small exact computations
in your scratch dir where helpful.

## C1. Why Hiary's exponent is 1/3: the integer-hull / curvature law

Hiary's Algorithm 1 (arXiv 1408.2608) uses blocks with H = ⌊(17N)^{-1/3} x⌋, i.e. ≈ 0.891·N^{1/3} blocks per
dyadic range [X, 2X], independent of X, for 0.214·N^{1/3} ln N blocks in total. The integer hull
conv{(x,y) ∈ Z²: xy ≥ N} has Θ(N^{1/3} log N) vertices (Balog–Bárány arXiv 2602.06897;
Alcántara–Blanco–Criado–Santos arXiv 2501.19193). Every divisor (d, N/d) is a vertex. Measured (exact, brute
force, n = 24..40 bits): 0.48–0.49·N^{1/3} hull vertices in every full dyadic window in [N^{1/3}, √N], and a
max/min window ratio ≤ 1.03. Claimed interpretation: 1/3 is the Jarník/Andrews curvature scale. A convex arc of
length L and radius of curvature ρ carries ≍ L ρ^{-1/3} hull vertices, and for xy = N near √N, L ≍ ρ ≍ √N.
Likewise for the circle x² + y² = N, where lattice points on it factor N when p ≡ q ≡ 1 mod 4: ≍ N^{1/3} vertices
(Balog–Bárány 1991 for discs). Claim: "hull-walk" algorithms are Ω(N^{1/3} log N) (PROVED, trivially, as output
size), and the 2-D conic geometry explains why Lehman, Hiary, McKee and Rubinstein all land at 1/3 (HEURISTIC /
taxonomy). Attack the taxonomy part especially: does it really explain Lehman and Rubinstein, or is it
numerology?

## C2. Reduction L_j: the floor-sum ladder gives N^{1/(j+2)}, not N^{1/(j+1)} or N^{1/(2j+2)}

Model: an oracle returns exact Σ_{h1≤h<h2} ⌊P(h)⌋ for rational degree-j P. Algorithm (src/analytic/ladder.py):
B = N^{1/(j+2)}, blocks [x, x+H) with H = ⌊c x/B⌋, x from B/c to √N (trial division below). P = degree-j Taylor
polynomial of f(h) = N/(x+h). The remainder R(h) = N(−h)^{j+1}/(x^{j+1}(x+h)) has constant sign and
|R| ≤ ρ = N H^{j+1}/x^{j+2} ≤ c^{j+1} B/x < 1. If x+h | N then {P(h)} ∈ [1−ρ, 1) (j odd, or h = 0), or
{P(h)} ∈ [0, ρ] (j even). The candidate set is counted by a difference of two floor sums, located by bisection,
and each candidate is tested by one division. Candidates d satisfy |N mod± d| ≤ 2ρd ≈ 2c^{j+1}B =: Y, so the
total over all blocks is ≤ Σ_{|k|≤Y} τ(N−k) ≪ Y log N by Shiu's Brun–Titchmarsh theorem for τ in short
intervals (Y ≥ N^ε). Total: O(N^{1/(j+2)} log N) blocks and O(N^{1/(j+2)} log² N) oracle calls. PROVED. At j = 1
(Euclid floor sums) this is an N^{1/3} algorithm (a cousin of Hiary). j = 2 → N^{1/4}, j = 3 → N^{1/5}, j = 4 →
N^{1/6}. Your blind report's construction demanded precision 1/(8M) and got N^{1/(j+1)}; nonsmooth's table says
N^{1/(2j+2)}. Which is right? Check Shiu's theorem's exact hypotheses (Shiu, J. reine angew. Math. 313 (1980)),
the sign/endpoint conventions, and whether the candidate bound needs any equidistribution. Also judge the
claimed model lower bound: "any method that sees the block only through a degree-j polynomial approximation,
pays 1 per block and 1 per ambiguous h, costs ≥ c_j N^{1/(j+2)} per dyadic range near √N", via Chebyshev's
bound E_j(f) ≥ (b−a)^{j+1} min|f^{(j+1)}| / (2^{2j+1}(j+1)!) balanced against blocks. HEURISTIC or provable?

## C3. Band-limited divisor detectors need bandwidth ≍ x

Claim (PROVED): let F(u) = Σ_{|k|≤K} c_k e(ku), F(0) = 1, with |F(r/d)| ≤ ε for all d ∈ I ⊂ [x, 2x] and
r ≢ 0 mod d, where ε < 1/(2W) and W = |I|, so that Σ_{d∈I} F(N/d) rounds to the divisor count. Then for any
d ∈ I with K < d/2: sampling F at r/d has no aliasing, and the DFT of the samples must vanish at the d − 2K − 1
frequencies with K < |k| < d − K. Parseval gives (d−2K−1)/d² ≤ (d−1)ε²/d, so 2K+1 ≥ d − d(d−1)ε². With
W ≥ √(2x) this forces K ≥ 3x/8 − 1. Hence evaluating Σ_k c_k Σ_{d∈I} e(kN/d) frequency by frequency costs ≳ x
sums, more than trial division over I. Check the algebra and the hypotheses.

## C4. Truncated Ramanujan / singular-series expansions carry nothing below the least prime factor

If f(N) = Σ_c a_c(N) c_c(N) with a_c depending on N only through |N| (archimedean size), then every term with
c < P⁻(N) has c_c(N) = μ(c), so the truncation at c < P⁻(N) is the same function of |N| as for a prime. For r_k
(k ≡ 0 mod 4), the odd-q singular-series terms are ε_q^k q^{−k/2} c_q(N), so this applies. For odd k, they are
Jacobi-symbol weighted, so information lives in L(·, χ_N) (class b) or N mod q (floor). Correct?

## C5. Complete cubic Weyl sums mod N factor N (and the degree-2 / degree-3 dividing line)

S(a,N) = Σ_{x mod N} e(a x³/N). For N = pq: S(a,N) = S(a q̄, p) S(a p̄, q). S(b,p) = 0 if p ≡ 2 mod 3, and
S(b,p) = 2√p cos(θ_p − 2πk(b)/3) if p ≡ 1 mod 3, where G(χ) = √p e^{iθ_p} and G(χ)³ = p J(χ,χ), N(J) = p. From
S(a^i b^j, N), i, j ∈ {0,1,2}, a rank-one 3×3 array 4√N u_α v_β up to an unknown F_3-linear relabeling
(48 options): z_u = Σ u_k ω^k ∝ e^{iθ_p}, z_u³ ∝ π_p (up to conjugation/units). Then γ = z_u³ z_v³ and
γ' = z_u³ conj(z_v³), rescaled to norm N and rounded in Z[ω], have gcd of norm p. Implemented in
src/analytic/cubic.py; it factored 100/100 N (14–22 bits, p ≡ q ≡ 1 mod 3; 95 structural, 5 accidental gcd).
With uniform noise δ on each S value it still succeeds ≥ 50% for δ ≤ 3 (16–22 bits), i.e. it needs absolute
precision O(1) on a quantity of size ≈ 4√N. Consequences claimed:
(i) Zero/non-zero of S(1,N) for N ≡ 1 mod 3 decides whether p ≡ q ≡ 1 or p ≡ q ≡ 2 mod 3, i.e. whether −3
is a QR mod N (the fixed-base QR bit), PROVED.
(ii) A polylog algorithm approximating complete cubic sums mod N to O(1) absolute error would factor every
N = pq with p ≡ q ≡ 1 mod 3 (a positive proportion), PROVED. So Hiary's cubic-sum theorem cannot extend from
|c| ≤ K^{-2} (μ ≤ 1) to the complete-sum range c = a/N, K = N (μ = 2) with polylog cost unless such N are
easy to factor.
(iii) Degree 2 complete sums mod N are Jacobi symbols (class c), so the reciprocity / Weil-representation wall sits
exactly between degree 2 and 3. Same for degree k with gcd(k, p−1), gcd(k, q−1) > 1 (Jacobi sums in Z[ζ_k]).
Is this known? (van Dam–Seroussi quant-ph/0207131 relate Gauss-sum estimation to discrete log; Kummer's
problem; anything on composite-modulus cubic sums and factoring?) Is (ii) correctly scoped? Could the precision
requirement be cut to N^{c} (c > 0) by a lattice decoder, e.g. Coppersmith on the norm-N conic?

## C6. Dequantisation verdict

Shor's amplitude A(k) = (1/M) Σ_{x<M, a^x ≡ y} e(kx/M) = (1/(MN)) Σ_h Σ_x e(kx/M + h(a^x − y)/N). Claims:
(i) the phase class {αx + h a^x/N} is closed under Weyl differencing (Δ_u: h → h(a^u − 1)), so the A-process never
lowers "degree", and there is no stationary-phase B-process (KNOWN-ish / HEURISTIC);
(ii) for the bit-split x = 2^s x_hi + x_lo, the matrix [e(h a^{2^s x_hi + x_lo}/N)] has rank
min(2^s, 2^{m−s}, r)-ish generically, so any MPS / tensor contraction along bit order has bond dimension
≳ min(2^{m/2}, r) (HEURISTIC; prove or refute the generic rank claim via characters of the cyclic group ⟨a⟩);
(iii) the QFT itself is classically simulable (Griffiths–Niu semiclassical QFT; Yoran–Short 2007; Browne 2007;
Aharonov–Landau–Makowsky 2006), so all hardness is in modular exponentiation;
(iv) black-box query bounds (Cleve 2004; Chakraborty–García-Soriano–Matsliah 2011) do not apply to explicit
modular multiplication.
Find errors, missing references, or an actual classical extraction I have missed.

## C7. Subexponential feasibility

The only analytic route to L[1/2] found is the uniform-degree ladder: N^{1/(j+2)} · T_j, with
T_j = exp(O(j log j)) polylog and j ≍ √(log N / log log N). The minimal OPEN primitive is the exact quadratic
floor sum Σ_{h<H} ⌊αh + βh² + γ⌋ in the "full-curvature short" regime βH² ≍ 1 (rational coefficients with
denominators ≈ x³). There the lattice points under the arc have ≍ H^{1/3} hull vertices, so hull walking costs
H^{1/3} per block and recovers exactly N^{1/3}. Polylog would give N^{1/4} (ties Pollard–Strassen and
Coppersmith); beating N^{1/5} needs j ≥ 4. Complete quadratic floor sums Σ_{x<p} ⌊a x²/p⌋ give h(−p) via
Dirichlet's formula Σ_{r<p} r (r/p) = −p h(−p) for p ≡ 3 mod 4, p > 3 (class b), but that is a different
regime. Is the short regime known to be hard or easy? Is there any reciprocity / Ostrowski-type recursion for
quadratic floor sums (cf. Hardy–Littlewood lattice points under a parabola, Bombieri–Iwaniec, "quadratic
Dedekind sums")? Also: Darren Li, arXiv 2609.29949, "Computing the divisor summatory function in N^{9/28+o(1)}".
Does it exist, what does it actually do, and could a restricted-cutoff variant locate divisors below N^{1/3}
without the Lenstra–Pomerance step?

Deliver a Markdown report: one section per claim with verdict ACCEPT / ACCEPT WITH FIXES / REJECT, the fix,
and references. Then list anything important that is missing.
