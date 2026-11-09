# Adversarial review (midk_2)

You are an adversarial referee. Attack the following claims about the "intermediate-k" factoring
barrier (cost L_N[α] = exp((c+o(1))(log N)^α(log log N)^{1−α}); heuristic record α = 1/3). Be
concrete; find errors, missing cases, or overclaims. Use exact arithmetic.

**Claim MAP.** In the NFS/index-calculus relation model, exactly three families carry a discrete
parameter k whose intermediate values are blocked only by a counting/independence hypothesis (not a
proved area/discriminant/Coppersmith wall), and **only one can beat L[1/3]**:
- (A) rank-r linked modules (two degree-d number fields, r ≪ d linked generators), blocked at
  intermediate rank 3 ≤ r ≤ d by a counting+independence hypothesis G_r (Schmidt point counts +
  independence of the r−1 linking congruences mod N). Favourable failure → L[1/4].
- (B) the multi-level descent ladder, cost α_j = (j+1)/(3j+2), blocked by pair-correlation C₁/C₂;
  claim: α_j > 1/3 for every finite j, inf 1/3 only as j→∞, so B cannot beat or reach 1/3.
- (C) the floor-sum ladder, deterministic metric N^{1/(j+2)} ≫ L[1/3], not a route to L[1/3].

**Claim LEDGER.** For (A): per-candidate log-norm log X = κ·x/d + d·c (x=log N, c=log(E√r) ≥ 1,
κ = height cost of linking). The rigorous 2-plane resultant gives κ ≥ 1, so AM–GM gives
log X ≥ 2√(κ c x). Hence (i) fallback κ=1 ⇒ α ≥ 1/4 rigorously, attained at d≈√(x/log x),
r≈x^{1/4}(log x)^{3/4}, and **no (d,r) gives α < 1/4**; (ii) under G_r, κ = r−1, ⇒ α = 1/3 with r=2
optimal; the G_r/KL height floors are non-decreasing in r.

**Claim SPL-INT.** A favourable G_r failure needs the r−1 linking conditions dependent at height
N^{1/d}. The only sources are (I) a global char-0 relation (= ring power/subfield, pays geometric
height) or (II) a dependence realised only mod N via a non-scalar Z/N intertwiner T. **Lemma SPL-INT:**
fraction-free Gaussian elimination of (T − λI) over Z/N either returns a global Z/N eigenvector or
hits a nonzero non-unit pivot whose gcd with N splits N. Hence an exploitable favourable failure is
factoring-equivalent to exhibit. (Witness: the cyclotomic automorphism-graph family links all r
generators at height 1 but the two sides are Galois-conjugate, so norms and modular square roots
coincide — no gcd.)

**Claim 2ND-ORDER.** The minor-product invariant Δ = ∏_{|I|=r} M_I² ∈ Z∖{0} with N^{2·C(D−2,r−2)} | Δ
(unramified linking) + Hadamard gives a *deterministic* height penalty h₂ ≥ (r−1)x/(D−1) whenever all
maximal embedding minors are nonzero. This re-imposes α ≥ 1/3 across the whole NONSINGULAR intermediate
range (no independence heuristic). The hole is exactly the SINGULAR case (some M_I = 0, i.e.
frame-less / rank-deficient), floored at L[1/4], where no construction is known.

**Questions.** (1) Is Claim MAP's trichotomy complete — is there a fourth intermediate-k family whose
favourable failure beats 1/3? (2) Is the 1/4 floor in Claim LEDGER truly a floor, or can a cleverer
ledger (unequal degrees, non-balanced profiles, p-adic/special-q sublattices, several fields, box side
> B) push below 1/4? (3) Is Lemma SPL-INT correct and does it actually cover every exploitable
favourable failure, or is there a dependence source it misses (e.g. a global non-power relation, or a
nonlinear/non-matrix intertwiner)? (4) Does Claim 2ND-ORDER's minor-product penalty really hold at
growing (d, r), and does the singular hole genuinely survive, or does a union bound over deleted
minors / lower sub-ranks close it? Cite literature (Schmidt 1993, Bean–Thunder 1996, Evertse–Győry,
Davenport, Thunder 2001, Coppersmith) where decisive.
