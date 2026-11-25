# coef: compute u_ℓ(N) = N-th coefficient of ζ(s)ζ(s−1)/L(E,s) mod ℓ without factoring N

Your own target from dig A (this session). The user believes there are good reasons (rumours) that classical
poly-time factoring exists and that we are missing something obvious. Work hard on THIS target; be creative, but be
precise. Repo (do NOT edit; scratch in your cwd).

Orchestrator's framing to test (Claude): Jacobi symbols are computable from N because Dirichlet characters are
completely multiplicative AND given by congruence conditions (reciprocity), so χ(N) = χ(p)χ(q) is read off N.
u_ℓ(N) = det(1−g_p)det(1−g_q) mod ℓ depends on Frobenius CLASSES of p, q in Gal(Q(E[ℓ])/Q) ⊆ GL₂(F_ℓ) — non-abelian;
from N alone one gets only abelian shadows (det g_p · det g_q = N mod ℓ). Expansion:
u_ℓ(N) = (1+p)(1+q) − a_p(1+q) − a_q(1+p) + a_p a_q mod ℓ, i.e. it needs mixed terms beyond a_N = a_p a_q.
So a native evaluator would be a "computable non-abelian reciprocity at composite index".
Tasks (≤ 2200 words; label proved / known / heuristic / speculative):
C1. Check and sharpen the framing. Which quantities attached to (g_p, g_q) ARE computable from N (abelian characters of
    Gal; anything else? e.g. via the Weil pairing, via traces in auxiliary representations Sym^k, ⊗, via Rankin–
    Selberg L(f×g) coefficients at N, via base change to quadratic fields where some classes become abelian, via
    CM where Gal is (virtually) abelian so a_p is a Hecke character value — for CM curves is u_ℓ(N) computable from N
    by reciprocity in the CM field? why not? (needs the local norm-p element, i.e. representing p by a form)).
    Is there a group-theoretic proof that no function of the abelian shadows determines u_ℓ(N) mod ℓ?
C2. Small cases exhaustively: ℓ = 2 (M mod 2 = cubic root-count condition at p and q; Gal ⊆ S₃): is M mod 2 computable
    from N for a fixed S₃ cubic? It is a function of the splitting types (Frobenius classes in S₃) of p and q; the
    abelian quotient S₃ → C₂ gives (disc/N). Decide whether the pair of classes, or just the needed event
    "both classes are 3-cycles", is determined by computable data, and whether computing it would factor N
    (a decision version: "is M odd?" — by PT/symmetry arguments, if it were computable, could you split N?). Then ℓ = 3.
C3. Mine the recent results for any lever: (i) Cogentic deterministic L[3/4] paper (codex/scratch/cogentic/
    deterministic_factoring_l34.txt; biquadratic cocycle positivity, Walsh–Hadamard nonnegativity forcing QNRs,
    algSplit from a multiple of p−1, Gauss cycles / short principal-ideal generators); (ii) Cogentic L[1/3]
    (codex/scratch/cogentic/integer_factorization.txt); (iii) the alpo.ge Lean mechanism (notes/multcopp_results.md:
    CM derivatives D^a θ(0) mod N vanish mod p for p inert, i.e. blindly detect supersingular reduction — a Frobenius
    property! — but only at order a ≥ p); (iv) OpenAI FFT paper (codex/scratch/openai_fft/paper.txt); (v) Masullo GRH.
    For each: does it give ANY handle on non-abelian Frobenius data at N, or on #E(Z/N), or on a_N?
C4. Attempt constructions: anything that evaluates a non-abelian class function at the pair (g_p, g_q) from N — e.g.
    via a composite-level modular form whose N-th coefficient is computable (Eisenstein-type, theta series of
    quadratic forms of discriminant related to E?), via the curve over the ring Z/N with an auxiliary algebra
    (EXT-style) where some operation's success probability depends on u_ℓ(N), via sampling: is u_ℓ(N) mod ℓ
    estimable from statistics of random-point operations mod N (e.g. frequency of [m]P hitting a denominator gcd),
    via ℓ-adic / Tate-module lifting mod N^k, via isogeny walks mod N. Run cheap pre-registered experiments where
    sensible (40–80 bits, exact ground truth from the factors), including a sampling-based estimator.
C5. Verdict: proved obstruction(s), the sharpest remaining open statement, and an honest probability.
