# crtdecode round 1 (blind): adaptive CRT decoding far below the Johnson radius, for smooth integers in short intervals

You are a research mathematician (coding theory, lattice algorithms, analytic number theory, query/decision-tree
complexity). Work rigorously, label every claim PROVED / HEURISTIC / MODEL / OPEN / UNVERIFIED, and do not inflate.
The repo at /Volumes/git/factoring-connection is read-only for you; scratch files go in your current directory. You may
run python (`/Volumes/git/factoring-connection/.venv/bin/python` has sympy, gmpy2, python-flint, numpy) for exact
checks. Use web search for literature: Goldreich–Ron–Sudan (IEEE Trans. IT 46, 2000) "Chinese remaindering with
errors"; Boneh (STOC 2000; JCSS 64 (2002)) "Finding smooth integers in short intervals using CRT decoding";
Guruswami–Sahai–Sudan (FOCS 2000) soft-decision CRT decoding; list-decoding capacity for CRT codes (Guruswami's
thesis; Guruswami 2003 number-field codes; Cohn–Heninger ANTS X 2013 "Approximate common divisors via lattices";
Chinburg–Hemenway–Heninger–Scherr ASIACRYPT 2016); Brakensiek–Chen–Putterman–Zhang–Zheng, arXiv 2609.08005 (2026),
"Algorithmic list decoding of Reed–Solomon codes up to capacity"; Jeronimo arXiv 2609.05870.

## Notation and target

x = ln X with X ≈ kN; F_α = x^α (ln x)^{1−α}; L[α, c] = exp((c + o(1))F_α). Factor base y = L[1/3, κ]; window
h = L[2/3, ν] with ν slightly above 2/(3κ); u₀ = x/ln y, so the density of y-smooth integers near X is
ρ(u₀) = L[2/3, −2/(3κ)]. Optimum κ* = (2/9)^{1/3} ≈ 0.606.

**Prop SI.** Given X = kN (and h, y), output a y-smooth u ∈ [X, X + h] with its factorisation. If this costs
L^{o(1)} per output, N factors in L[1/3, 1.2114] (heuristic). With cost L[1/3, τ] per output the factoring constant is
min_κ max(κ + τ + 2/(9κ²), 2κ); it beats MNFS (1.9019) iff τ < 0.757. Equivalently, a decoder costing π(y)^a per window
(L[1/3, aκ]) beats MNFS iff a < 1.1416 (after re-optimising κ; a = 1 gives 6^{1/3} = 1.817).
The baseline (interval sieve, or the "divisor trick" z = d⌈X/d⌉ with d ≈ h a power of 2 and the cofactor tested)
costs 1/ρ(u₀) = L[2/3, 2/(3κ)] per output. Promise version: the window is assumed to contain a smooth integer.

**CRT view.** Let P = ∏_{q ≤ y} q^{e_q} with q^{e_q} ≤ 3X < q^{e_q+1}. The candidate offsets j ∈ [0, h] form a
CRT code c(j) = (j mod q^e)_{q,e}; the received word is (−X mod q^e). The weighted agreement (weights ln q) of
c(j) with the received word is A(j) = ln(y-smooth part of X + j) = ln gcd(X + j, P). Prop SI asks for j with
A(j) ≥ ln(X + j) ≈ x. Total weight ln P ≍ π(y)·x = L[1/3, κ]; message length ln(h+1) = νF_{2/3}.
Known decoders (Coppersmith/Howgrave-Graham/Boneh/GSS) certify j iff A(j) ≥ √(ln(h+1)·ln P) (the Johnson bound),
which here is ≫ x by an L[1/3] factor.

Known results you may take as given:
- **Theorem J** (proved, averaged over X): every single-integer Johnson-class certificate decoder (any y-smooth
  modulus, any weights, any lattice dimension) whose parameters are fixed in advance certifies on average
  ≤ ρ(u₀)^{1−o(1)} smooth integers per call. The divisor trick attains this. Fixed adaptive TREES of such calls with M
  distinct nodes give ≤ M·ρ(u₀)^{1−o(1)}.
- **Information is cheap (Lemma IC):** the residues of X modulo all q^{e_q}, q ≤ y, cost L[1/3, κ] bit operations and
  determine exactly which positions in the window are smooth. So no query/information barrier can exist; any barrier
  must be computational.
- **Prop CV:** squarefree Prop SI ⇔ promise CVP_{1+O(1/x)} in Adleman/Schnorr/Micciancio's prime-number lattice of
  dimension π(y), on a centred target.
- Existence for EVERY window at this scale would imply Vinogradov's least-QNR conjecture (so only promise or average
  versions are realistic).
- The explicit-prime-list version (smooth over a given set S in an interval) is NP-hard under randomized reductions
  (DIVISOR / subset product). With S = all primes ≤ y nothing is known.

## Questions (answer each; be quantitative)

1. **Formalise** the exact CRT list-decoding problem Prop SI needs (code, coordinates, weights, rate, agreement
   fraction, list size — worst case vs average over X — and the interval constraint). Place it on the map:
   agreement fraction vs the Johnson radius and vs list-decoding capacity for CRT codes. Compare precisely with GRS,
   Boneh, GSS, the capacity results for CRT/number-field codes, and the 2026 plain-RS capacity decoder (what parameter
   of that decoder blows up here, and why it would not transfer to CRT codes even if it did not blow up).

2. **Adaptive decoders.** Smooth numbers have multiplicative structure: u ≡ 0 mod many q means u is a multiple of
   their product M, so the decoder can reduce to "find a multiple of a large smooth M in the window whose cofactor is
   smooth" (the divisor trick, nested). Define a natural adaptive model: a decision tree whose internal nodes are
   queries (e.g. "give me the multiples of M in [X, X+h]" or "give me ⌈X/M⌉ mod q", or one Coppersmith/Boneh call
   with modulus M), chosen as a function of all previous answers, each costing 1 (or poly(log X)). Can ADAPTIVELY
   chosen sequences of moduli beat the divisor trick's ρ(u₀)^{−1} per output? Analyse carefully; note that X itself is
   known to the algorithm, so "adaptivity" must be measured in a model that charges for computation, not for
   information.

3. **Either a lower bound in a natural adaptive model, or a construction with a < 1.1416.** Candidate frameworks:
   query complexity with a divisibility-test oracle (but the residues are free — make the model non-trivial);
   lattice-based adaptive decoders (sequences of Coppersmith calls whose moduli depend on earlier lattice outputs);
   reductions to/from subset-product / subset-sum / CVP; generic-group or algebraic-decision-tree models; any known
   lower bound for CRT/RS list decoding beyond Johnson (e.g. Ben-Sasson–Kopparty–Radhakrishnan list-size bounds,
   Guruswami–Rudra limits for RS, NP-hardness of maximum-likelihood RS decoding by Guruswami–Vardy, bounded-distance
   decoding hardness, Gandikota–Ghazi–Grigorescu). Which of these transfer, and what exactly do they say here?
   Be explicit about what is PROVED, in which model, and what remains OPEN. If you find a construction, compute its
   exponent a (or τ) honestly, including all L[2/3] factors.

4. **Toy experiment design.** Propose 2–3 exact toy experiments (40–64-bit X, y ≈ 2^8–2^12, h ≈ 2^10–2^16) comparing
   adaptive decoders against the divisor trick on "smooth outputs per unit cost", with a cost measure that cannot be
   gamed by toy-size artefacts (e.g. h ≥ √X at toy size). Say what result would surprise you.

5. Give a subjective probability that a decoder with a < 1.1416 exists, and the single most promising mechanism.

Write your answer as a self-contained report (≤ 4000 words), sections 1–5.
