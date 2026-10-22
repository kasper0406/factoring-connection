# sinonobliv round 1 (blind): non-oblivious constructors for smooth integers in short intervals near kN

You are a research mathematician (analytic number theory, lattice algorithms, factoring). Work rigorously, label every
claim PROVED / HEURISTIC / MODEL / OPEN / UNVERIFIED, and do not inflate. The repo at
/Volumes/git/factoring-connection is read-only for you; scratch files go in your current directory. You may run
python (`/Volumes/git/factoring-connection/.venv/bin/python` has sympy, gmpy2, python-flint, numpy) for exact checks.
Use web search for literature (Friedlander–Lagarias 1987, Cohn–Heninger ANTS X 2013 "Approximate common divisors via
lattices", Herrmann–May, Boneh 2000/2002 "Finding smooth integers in short intervals using CRT decoding",
Sarkar–Maitra ePrint 2010/146, Balog / Harman / Granville survey 2008 MSRI 44 §4–5, Coppersmith's "Modifications to the
NFS", SNFS/MNFS constants).

## The problem (Prop SI)

Notation: x = ln N (or ln X for the size X ≈ kN), F_α = x^α (ln x)^{1−α}, L[α, c] = exp((c + o(1))F_α).
Factor base y = L[1/3, κ]; window h = L[2/3, ν] with ν > 2/(3κ); intervals I_k = [kN, kN + h], k ≤ L[1/3, C].

**Prop SI** (model-proved under heuristics H_S on the residue and R): if one can construct a y-smooth u ∈ I_k (with
its factorisation) at cost L^{o(1)} per output, then u ≡ v = u − kN (mod N) with v ≤ h tested for smoothness gives
relations, and N factors in L[1/3, (16/9)^{1/3}] = L[1/3, 1.2114] (κ* = (2/9)^{1/3}, ν* = 2/(3κ*) = 1.1006).
With per-output cost L[1/3, τ] the constant is min_κ max(κ + τ + 2/(9κ²), 2κ); it beats MNFS (1.902) iff τ < 0.757.
The baseline (interval sieve) costs 1/ρ(u₀) = L[2/3, 2/(3κ)] per output, u₀ = x/ln y.

Known results from earlier work (take as given; they are about OBLIVIOUS methods):
- **Theorem J** (proved, averaged over the window position): any *single-integer* Johnson-radius certificate decoder
  (Coppersmith / Howgrave-Graham / Boneh / GSS weighted CRT decoding; certified set {z ∈ R : Σ_q min(v_q(z), e_q) w_q
  log q ≥ √(Q log|R|)}, Q = Σ e_q w_q² log q), with modulus P and weights fixed in advance (or from a repertoire of
  size L[2/3]^{o(1)}), certifies on average ≤ L[2/3, −2/(3κ) + o(1)] smooth integers per call. The divisor trick
  z = d·⌈x/d⌉ (d ≈ h a power of 2, cofactor tested) attains this. Herrmann–May batching over k ≤ h is covered.
- **Scale invariance** (Lemma SR): if u = d·j ∈ [X, X + h] then d lies in a window of the same relative width h/X.
  An independently placed factor s ≲ Kh buys only the smoothness of a cofactor of size ≥ N/(Kh) (Lemma RC).
- **Fixed two-list product tables** (MITM, all k at once): ≥ √(X/3)/(h + 1) cost per output (exponential).
- **Prop CV**: squarefree Prop SI ⇔ promise CVP_{1+O(1/x)} in Adleman/Schnorr's prime lattice of dimension π(y).
- **Consequences**: a modulus-generic constructor (works for every modulus) would also give F_p discrete logs in
  L[1/3, ≤ 1.232]; all-interval *existence* at this scale implies Vinogradov's least-QNR conjecture for p ≡ 3 (mod 4).
  So only promise constructors are realistic.
- Exact *split/compression* families (polynomial identities producing products of linear factors near a target) are
  limited to ≤ 1 relation per presentation by rigidity results (Lemma R / Theorem W: exact identities h_A P_A = F Q_A +
  r_A on one presentation F give one rational function). Any gain must come from *smoothness of the factors*, not
  from exact splitting.

## The open routes. For EACH: construct something and compute its per-output cost τ (as L[1/3, τ], or say it is
## L[2/3, ·] / exponential), OR prove a barrier in a clearly stated model.

1. **Adaptive choice of moduli/lattices.** Theorem J assumes parameters fixed in advance. Can an adaptive sequence of
   Coppersmith/Boneh calls beat the divisor trick, where each call's modulus depends on previous outputs (e.g. building
   a smooth number digit by digit or prime by prime with lattice steps: choose d₁ | z, then condition on z ≡ 0 mod d₁
   and look for d₂ in the induced progression, etc.)? Formalise a natural adaptive model (decision tree over
   divisibility / gcd queries, or a sequence of certificate calls whose parameters depend on earlier answers) and give
   its lower bound, or find a gain. Think about whether the information-theoretic count (adaptive tree with L^{o(1)}
   leaves) is the right measure, and whether averaging over x still gives a bound for adaptive strategies.
2. **Friedlander–Lagarias-type correlated identities.** Polynomial identities that put many integers near a target
   with structured, correlated factorisations: products of consecutive integers t(t+1)…(t+r−1); F(t) with forced
   factorisation like t² − a² = (t − a)(t + a), t(t+1); Aurifeuillian factorisations; (m − a)(m + a) = m² − a².
   Can a family place many y-smooth integers in [kN, kN + h]? Compute the density gain over random (smoothness
   probability of the structured value vs a random integer of the same size) AND the placement cost (how many members
   fall in the window). Note: the product of r factors of size X^{1/r} each is smooth with probability ρ(u/r)^r ≈
   ρ(u)·r^{u}, but how many members land in a window of relative width h/X? Is there any identity family where both the
   gain and the hit count work at L-scale? What is the actual Friedlander–Lagarias 1987 construction (x^{1/2}-smooth in
   intervals of length x^{1/4+ε}), and how far does that method go toward y = L[1/3], h = L[2/3]?
3. **Multi-sample Cohn–Heninger decoders.** Approximate common divisors: m samples a_i = q_i·d + r_i with |r_i| < R,
   d ≥ X^β unknown, recoverable for R up to X^{β^{(m+1)/m}}. Here the windows I_{k_1}, …, I_{k_m} (or shifted windows)
   are decoded jointly for a common smooth divisor. Is there a gain beyond Theorem J's single-integer bound, in the
   averaged count of certified smooth outputs? Note that a common divisor d of integers in m windows is a much rarer
   event; compute the expected number of certified smooth outputs per call averaged over N.
4. **Smooth norm or polynomial values near (kN)^{1/d}.** Find u = F(t) or u = Norm(α) ≈ kN with F of degree d
   (t ≈ (kN)^{1/d}) or α in a number field of degree d, whose value is a product of small norms / small factors, so
   smooth with higher probability. E.g. u = t^d − a with a small? u = Norm(a − bθ)? u = ∏(t − r_i)? The window forces
   F(t) ∈ [kN, kN + h]: how many (t, parameters) hit, and how smooth is the value? Cost this in the NFS-style sampling
   model (size–entropy law: total log-size λ of tested integers vs log of the number of candidates C; cost
   e^{C*}, C*² = 2λ log λ), and check whether it reproduces SNFS-like constants (1.526) or MNFS / GNFS ones, or does
   better than 1.2114-route needs.
5. **The softer target.** Report the per-output τ for every route; τ < 0.757 already beats MNFS.

Also: "uncertified decoder outputs" — are there decoders whose uncertified roots are biased toward smooth numbers?

## Output format
- Verdict first (plain language).
- Per-route table: route, best construction found, per-output cost (τ or exponent), status label, the model in which
  any barrier holds.
- Then derivations, with every inequality written out. Exact constants where possible.
- Literature: what Friedlander–Lagarias, Cohn–Heninger, Herrmann–May, Boneh, Sarkar–Maitra actually prove, with the
  exact statements you checked (say UNVERIFIED if not checked).
- Plausibility (%) that some route yields τ < 0.757.
