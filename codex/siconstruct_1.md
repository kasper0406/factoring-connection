# Blind round: constructing smooth integers in short intervals near multiples of N

You are a research mathematician (computational number theory, lattice algorithms, knapsack / generalized-birthday
algorithms). Work independently and rigorously. Label every claim PROVED (with proof), KNOWN (with a primary
citation you have checked, with theorem/section numbers), HEURISTIC (name the assumption), or UNVERIFIED. Web search
is enabled; use primary sources. Your working directory is a scratch directory; you may write and run small Python
scripts there (exact integer arithmetic) to test claims at toy size.

## Notation

x = ln N, F_α = x^α (ln x)^{1−α}, L[α, c] = exp((c + o(1)) F_α). Dickman: a uniformly random integer near X is
y-smooth with probability ρ(u)^{1+o(1)}, u = ln X/ln y, −ln ρ(u) = (1+o(1)) u ln u.

## The target (call it Prop SI)

Factor base y = L[1/3, κ] with κ near (2/9)^{1/3} ≈ 0.606. Window h = L[2/3, ν] with ν slightly above 2/(3κ) ≈ 1.10,
so that each interval I_k = [kN, kN + h] contains on average L[2/3, ν − 2/(3κ)] → ∞ y-smooth integers. Multipliers
k ≤ L[1/3, C] may be CHOSEN freely. Wanted: an algorithm that outputs a y-smooth u_k ∈ I_k (with its factorisation)
in time L^{o(1)} = L[1/3, o(1)] per output, for at least a 1/L^{o(1)} fraction of the k it tries. ("Promise
constructor": it may fail on intervals with no smooth element.) About L[1/3, κ + ν/(3κ)] outputs are needed in
total, because each relation u ≡ v := u − kN (mod N) is useful only when v ∈ [0, h] is also y-smooth (probability
L[1/3, −ν/(3κ)], assuming v behaves randomly). With linear algebra L[1/3, 2κ], such a constructor would factor N in
L[1/3, (16/9)^{1/3} + ε] ≈ L[1/3, 1.2114], far below the number field sieve. Even a per-output cost L[1/3, τ] with
τ < ~0.75 would beat MNFS (1.9019).

The baseline: sieving the interval costs h per interval and yields L[2/3, ν − 2/(3κ)] smooth numbers, i.e.
L[2/3, 2/(3κ)] per output. That is far too slow.

## Questions (answer each, numbered)

1. **Log-knapsack / lattice.** u = ∏ p_i^{e_i} ∈ I_k ⟺ Σ e_i ln p_i ∈ [ln kN, ln(kN + h)], a knapsack with relative
   precision h/(kN). Relate this precisely to CVP in a prime-number lattice (Schnorr 1991 / Adleman 1995 style,
   dimension π(y) = L[1/3]). What approximation factor γ is needed, and is the one-sided problem (v = 1, k free)
   easier or harder than Schnorr's two-sided problem (u, v both smooth, |u − vN| small)? Be quantitative. What would
   BKZ-β or any known CVP algorithm cost here?
2. **Meet-in-the-middle / multiplicative birthday.** u = a·b with a, b y-smooth from sorted lists A, B (sorted by
   ln). Compute the number of hits as a function of |A|, |B|, h, N, K (the number of multipliers k served at once).
   Crucially: does one sorted list amortise over many k (or many N-related targets) so that the cost PER OUTPUT
   becomes L^{o(1)}? Give the exact amortised per-output cost. Then consider k-fold products (Wagner's k-tree in log
   space, Schroeppel–Shamir, Howgrave-Graham–Joux representations, BCJ). Do the many factorisations of a smooth number
   (representations) help?
3. **Joint choice of (k, u).** We only need some y-smooth u with u mod N ∈ [0, h]. Is there a way to exploit the free
   k (or free multipliers / cofactors) that beats single-interval search? Be careful: an algorithm that listed y-smooth
   u with u mod N ≤ L[2/3] at L[α] per output for α < 1/3, or with u mod N ≤ y, would give a sub-L[1/3] or polynomial
   factoring algorithm; if your method seems to do that, look for the bug.
4. **CRT / Coppersmith / Boneh smooth-in-interval decoding** (Boneh STOC 2000 "Finding smooth integers in short
   intervals using CRT decoding"; Howgrave-Graham; Guruswami–Sahai–Sudan). What is the best per-output cost in this
   regime? Note the window h = L[2/3] is large relative to y.
5. **Anything else**: continued fractions, Farey/Stern–Brocot, products of consecutive integers, binomial
   coefficients, ratios of smooth numbers, special-form k, the base-m expansion of N, recursion / self-reduction,
   choosing k so that kN has special structure. Try hard to find a genuinely new construction. If you believe no
   L^{o(1)} (or L[1/3, < 0.75]) constructor exists within some broad class of methods, state the class and prove (or
   give the heuristic argument for) the barrier, e.g. a lower bound on the size of the part of u that must be tested
   for smoothness rather than constructed.
6. **Literature.** Is "find a y-smooth integer in [X, X + h]" with y = L_X[1/3], h = L_X[2/3] studied anywhere
   (algorithms or hardness)? Is the reduction "such a constructor ⇒ factoring in L[1/3, 1.2114]" in the literature?
   Is this just MPQS / lattice sieving / a known method in disguise?

Give a final table: approach | per-output cost (heuristic, L-notation) | meets L^{o(1)}? | status label. Then a
one-paragraph verdict.
