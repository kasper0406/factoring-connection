# Blind round: factoring via approximate CVP/SVP in "prime number lattices"

You are a research mathematician (computational number theory + lattice algorithms). Work independently and rigorously.
Label every claim PROVED (with proof), KNOWN (with a primary citation you have checked; give section/theorem numbers),
HEURISTIC (name the assumption), or UNVERIFIED. Do not trust secondary summaries; this area contains a famous erroneous
claim (Schnorr 2021, ePrint 2021/933 and later versions). Web search is enabled: please read the primary sources
(Schnorr 1991 "Factoring integers and computing discrete logarithms via diophantine approximation", EUROCRYPT'91 /
Adv. Comput. Complexity Theory 1993; Adleman 1995 "Factoring and lattice reduction" (manuscript); Schnorr 2021;
Ducas 2021 SchnorrGate repository and write-up; Ritter–Rössner 1997; Aboumrad–Widdows–Kaushik arXiv 2308.07804;
Yan et al. arXiv 2212.12372 and Khattar–Yosri arXiv 2307.09651; any 2022–2026 analyses, e.g. "lattice-based factoring
revisited" style papers, Aono–Nguyen-type estimates). Quote exact parameters when you rely on them.

## Setting

N odd composite, not a prime power, n-th prime p_n, factor base P = {p_1..p_n} (possibly p_1 = −1). Schnorr's prime
number lattice: basis rows b_i = (√f(i)·e_i, ⌈C·ln p_i⌋) for i ≤ n (f(i) a weight, e.g. f(i)=1 or a random
permutation of 1..n, or ln p_i), target t = (0,…,0, ⌈C·ln N⌋), scale C = N^c. A lattice vector b = Σ e_i b_i with
e ∈ Z^n gives u = ∏_{e_i>0} p_i^{e_i}, v = ∏_{e_i<0} p_i^{−e_i}, and one hopes |u − vN| is p_n-smooth, which gives a
"fac-relation" u ≡ (u−vN) mod N with both sides p_n-smooth. Need ≈ n+1 (or π(p_n)+small) relations, then linear
algebra mod 2 and a gcd.

## Questions (answer each, numbered)

1. **Exact reduction.** State as a precise theorem (or explain why none exists): for which approximation factor γ for
   CVP (or which SVP / BDD / list-CVP quality) in which explicit lattice (dimension n, weights, scale C) a single oracle
   call yields one relation, with |u − vN| ≤ ? Derive the exact identity linking |u − vN| to the distance
   ‖b − t‖ (the last coordinate error and the weighted ℓ2 or ℓ1 norm of e). How small must |u − vN| be to be
   p_n-smooth with useful probability (use Dickman ρ / Canfield–Erdős–Pomerance)? Then: how many oracle calls/relations
   are needed, how do you get *distinct, independent* relations from a deterministic oracle (randomised weights?), and
   what is the total factoring cost as a function of the oracle factor γ(n)? Find the threshold γ(n) (as a function
   of n and N) at which the method reaches (a) polynomial time, (b) L_N[α] for α < 1/3, (c) L[1/2]. Be explicit about
   the counting: how many lattice vectors exist in the relevant ball (Gaussian heuristic / exact count of p_n-smooth
   u,v with |u−vN| small), and is a *close* vector likely to be a *smooth-residual* vector at all?

2. **Rigour.** Which steps are proved, which rest on Gaussian heuristic, smoothness-of-residual heuristics,
   independence heuristics? State the cleanest provable statement you can (e.g. "an exact CVP oracle for lattices of
   dimension n=(log N)^c gives factoring in time …"), including whether an *exact* CVP oracle in Schnorr's lattice
   provably returns relations at all.

3. **Schnorr 2021.** Where exactly does the claimed polynomial-time (or "destroys RSA") argument fail? Compare the
   number of close vectors that actually give smooth relations with the claim. What does the approach really cost
   with BKZ-β (as a function of β, n, log N), including the root-Hermite factor δ_β ≈ (β/(2πe))^{1/(2(β−1))} and BKZ
   cost 2^{0.292β} (sieving) or β^{β/8} (enumeration)? Is it ever better than L[1/2]? Than the quadratic sieve? Note
   an elementary fact: for nonzero integers u, v, w with u − vN = w and gcd(uv, N) = 1, what lower bound on |u·w|
   holds, and what does that imply for the smoothness probability of a relation (both u and w must be smooth)?

4. **Consequences.** If approx-CVP with factor γ(n) = poly(n), or 2^{n^ε}, or 2^{n/log n}, were solvable in poly(n)
   time on these specific lattices, what would factoring cost (heuristically, and provably)? Conversely, does faster
   factoring imply anything for lattice problems (SVP/CVP/LWE) in these or general lattices? Be careful about
   worst-case vs. these structured instances. Give a clean statement "factoring ≤ γ-CVP on prime lattices with
   γ = …" or explain precisely why no clean reduction exists. Is the reduction direction ever known the other way
   (e.g. lattice problems on prime lattices being easy / hard relative to factoring)?

5. **Variants.** Is there a variant that is plausibly better: NFS algebraic side (CVP in a lattice encoding
   ln|norms| of ideals), two lattices, using divisor classes mod N (u ≡ vN + w with w restricted), Coppersmith-style
   small root formulations, list decoding, Micciancio–Voulgaris / sieving on these specific lattices, or exploiting
   that the prime number lattice is "almost" a scaled Z^n? Cost each one honestly in L-notation.

Write a self-contained report (≤ 4000 words), with formulas. End with a table: claim | status | source.
