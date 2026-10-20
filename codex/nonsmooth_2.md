# nonsmooth_2: adversarial referee of route-3 barrier claims

You are an adversarial referee (number theory + complexity), with web search. Repository
/Volumes/git/factoring-connection is read-only context. Your earlier blind report for this route is
codex/nonsmooth_1.out; read it. The claims below are agent `nonsmooth`'s draft. Your job: BREAK them.
For each claim: is it true as stated? Is the proof complete? Is the scope stated correctly (what it
does NOT cover)? Is it already in the literature (give references; it is fine and expected if these
are folklore, we only want correct attribution)? Then: is there any non-smooth paradigm, or any
variant of the killed ones, that these claims wrongly kill or simply miss?

Notation: N = pq, p < q ≤ 2p unless stated. "Operational" = sees N only.

## Claim OA (order-annihilation dichotomy)
Let G be a smooth commutative group scheme of finite type over Z[1/M0] with connected fibres. An
"order-annihilation" method computes [m]x for x ∈ G(Z/N), m an integer it can write down, and
detects x ↦ 0 modulo a prime p | N by a gcd. Write #G(F_p) = p^e · m_p with p ∤ m_p.
(a) By Chevalley's theorem and Lang's theorem (H^1(F_p, L) = 0 for connected L), with L = T × U the
linear part (commutative, perfect field) and A the abelian quotient,
#G(F_p) = #T(F_p) · p^{dim U} · #A(F_p), and #T(F_p) = |det(p − F | X*(T))| is prime to p and
≥ (p−1)^{dim T}.
(b) If ℓ is a prime dividing m_p with ℓ ∤ m, then Pr_x[[m]x ≡ 0 mod p] ≤ 1/ℓ.
(c) Hence with m = N^k·S, S composed of primes ≤ B, success with probability > 1/B needs m_p to be
B-smooth (Pollard p−1, Williams p+1, Bach–Shallit cyclotomic, ECM, ...). The only "smoothness-free"
case m_p = 1 forces dim T = 0 (for p > 2) and A(F_p) a p-group ("anomalous"); and splitting needs
the same not to hold at q.
(d) Coverage heuristics: a fixed non-CM curve over Q has a_p = 1 for ≍ √X/log X primes p ≤ X
(Lang–Trotter, conjectural), so "anomalous ECM" with K curves covers ≈ K/√p of primes: cost ≍ √p
= N^{1/4} for constant coverage (no better than rho). CM curves with D ≤ Y (Cheng 4p−1) cover
≍ √(Y/p) of primes at cost ≥ Σ_{D≤Y} h(−D) ≍ Y^{3/2}. For dim A = g ≥ 2, P_A(1) = p^j is rarer
(heuristically ≍ p^{−(?)} per variety; please supply a correct heuristic exponent or say UNKNOWN).

## Claim FROB (Frobenius operators over Z/N are factoring-hard)
Let B be an étale (Z/N)-algebra given by structure constants, and Φ_B its componentwise Frobenius
(Φ_B ⊗ F_p = Frob_p for each p | N), a (Z/N)-algebra automorphism.
(a) For B = (Z/N)[y]/(y^2 − d) with (d/N) = −1: Φ(y) = εy, ε = CRT((d/p), (d/q)), and
gcd(ε − 1, N) is a proper factor. (Deterministic: least such d.)
(b) For B = (Z/N)[x]/(f), f monic squarefree mod N with (disc f / N) = −1: over F_p, dim of the
fixed algebra ker(Φ − 1) equals the number r_p of irreducible factors of f mod p, and
(−1)^{deg f − r_p} = (disc f / p). So r_p ≢ r_q (mod 2), the (Z/N)-module ker(Φ_B − 1) is not free,
and Gaussian elimination / Howell form on Φ_B − 1 over Z/N meets a nonunit pivot: a factor.
(c) For B ⊇ (Z/N)[ζ_ℓ] with p ≢ q (mod ℓ): Φ(ζ) − ζ^{p mod ℓ} vanishes in one CRT component only;
trying the ℓ − 1 exponents factors N.
(d) Consequence (our reading): Schoof/SEA, Dokchitser resolvents, the Frobenius-identification
step of Edixhoven–Couveignes/Bruin, and Kedlaya/Lauder/Satoh-type Frobenius lifts all consume Φ on
such algebras at prime characteristic, where it is x ↦ x^p; an algorithm producing the same
operator at composite N factors N. The substitute x ↦ x^N is not a ring endomorphism.
(Is (d) correctly scoped? E.g. Schoof uses (Z/N)[x]/ψ_ℓ; does its Frobenius determine (det = p mod
ℓ) and is (b) or (c) applicable there? Lenstra 1985 "Galois theory and primality testing" and
Grantham's Frobenius pseudoprimes are the nearest literature we know.)

## Claim EIS (the Hurwitz aggregate is σ(N) in disguise)
C(N) := Σ_{t^2<4N} H(4N − t^2) = Σ_{ab=N} max(a, b) = 2N + 2q (Kronecker–Hurwitz; Gross–Zagier
(4.4)). Your nonsmooth_1 survivor: approximate C(N) within N^{1/4−ε} in L[1/4], then Coppersmith.
(i) Writing H(|Δ|) = (√|Δ|/π) L(1, Δ) and expanding L(1, Δ_t), Δ_t = t^2 − 4N, as Euler products: for
ℓ ∤ 4N the average over t mod ℓ of the local factor is an N-independent constant; for ℓ = p,
Δ_t ≡ t^2 (mod p) is a nonzero square for t ≢ 0, so the local factor is (1 − 1/p)^{-1} for a 1 − 1/p
fraction of t. This bias, times the main term ≈ 2N, is the 2q. Hence any approximation that uses
Euler factors / local densities only at primes ℓ ≤ X < p ("singular series") has error ≍ 2q ≍ √N,
and capturing it needs the local factor at the unknown p (or all primes up to p).
(ii) In weight 2, level 1 there are no cusp forms, so the Eichler–Selberg/Kronecker–Hurwitz identity
has no spectral term: C(N) has no representation other than (class numbers) = (divisor sum).
(iii) Hence the survivor is "approximate σ(N) to within N^{1/4}", i.e. the Coppersmith hint
itself; killed as a mechanism (not as a computational problem). Is this kill valid? Is there a
representation of C(N) (theta lifts, Zagier's weight-3/2 mock modular Eisenstein series,
Shintani zeta functions, Cohen–Zagier, Rademacher-type exact formulas for coefficients of mock
modular forms) that evades (i)? NB: Rademacher-type exact series for Fourier coefficients of
(mock) modular forms involve Kloosterman sums mod c for c up to ... — would truncating such a
series give C(N) to N^{1/4} cheaply, and if not, exactly why not?

## Claim RANK (minimal primitive)
For every invariant in the equivalence map whose prime values are easy (φ, σ, ψ, r_4, #E, #T(Z/N),
τ, a_f, dim S_k(Γ0(N))), its value mod a single small prime ℓ already determines the unordered pair
{p, q} mod ℓ in the reducible (Eisenstein) cases, hence the fixed-base quadratic-residuosity bit
QR_N(ℓ*) when (N/ℓ) = 1; and the computable-by-reciprocity part is only χ(N) for Dirichlet
characters (function of N mod conductor). Is "fixed-base QR mod N" (is −1, −3, ℓ* a square mod N)
known to be easier than, equivalent to, or independent of factoring? References?

## Also
- Our verdict: no smoothness-free general-purpose route below L[1/3] survives; the non-smooth
  methods that exist (Fermat/Lehman, Pollard–Strassen/BGS, Hiary N^{1/3}, Coppersmith N^{1/4},
  Harvey–Hittmeir N^{1/5}, SQUFOF N^{1/4}, Cheng/anomalous CM, Shor) are exponential, special-
  purpose, or quantum. Is there any classical non-smooth method with a proved or heuristic
  subexponential bound for general N that we are missing (e.g. from the generic-ring, real
  quadratic infrastructure, Witt-vector, p-adic, or lattice literature)?
- Label every point TRUE / FALSE / NEEDS-FIX / UNVERIFIED, with a one-line reason, then give the
  fixes. Be concise.
