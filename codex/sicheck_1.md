# sicheck round 1 (BLIND): validate a conditional L[1/3, 1.2114] factoring reduction ("Prop SI")

You are a research mathematician (computational number theory, smooth numbers, index calculus, sieve algorithms)
reviewing a speculative claim in a project on asymptotically faster classical integer factorization. Scepticism and
precision matter more than a positive result. Label every claim PROVED (give the proof), KNOWN (cite precisely; write
UNVERIFIED if you could not check), HEURISTIC (name the assumption), model-PROVED (exact consequence of a named model),
or OPEN. Do not read other files of this repository; work from this prompt, your knowledge and web search. Scratch
computation (Python, exact integer arithmetic; you may run small experiments up to a few minutes and < 1 GB RAM) goes in
the current directory. Final answer: a structured report (≤ 6000 words) with numbered sections matching the tasks, ending
with a list "requirements the constructor must satisfy for the claim to hold".

## Notation

x = log N, F_α = x^α (log x)^{1−α}, L[α, c] = exp((c + o(1)) F_α). Single N (an RSA-type semiprime p·q is the main
case), RAM model. Records: GNFS 1.923, MNFS 1.9019, SNFS (special N) 1.526, factorization factory 1.6386 per N after
precomputation.

## The claim (Prop SI)

- Factor base: primes ≤ y, y = L[1/3, κ]. Intervals I_k = [kN, kN + h], h = L[2/3, ν], k = 1, 2, ….
- Hypothesis: an oracle ("constructor"), given kN (and h, y), returns in L^{o(1)} time a y-smooth u_k ∈ I_k, with its
  factorization. (Only the I_k are queried; the mean number of y-smooth integers in I_k is L[2/3, ν − 2/(3κ)], so
  ν > 2/(3κ) is needed for existence.)
- Each output gives the congruence u_k ≡ v_k (mod N) with v_k = u_k − kN ∈ [0, h]. When v_k is also y-smooth this is a
  relation; F₂-dependencies among relations give X² ≡ Y² (mod N).
- Claimed cost: L[1/3, max(κ + ν/(3κ), 2κ)], infimum (16/9)^{1/3} ≈ 1.2114 at κ³ = 2/9, ν ↓ 2/(3κ). Claimed also: with
  "single intervals" of this kind, no choice of smoothness level L[α] gives an exponent below 1/3.
- Stated assumptions: H_S (v_k is "random", so y-smooth with probability ρ(log h / log y)) and R (enough relations
  give fruitful, full-rank dependencies).

## Tasks

1. **Re-derive independently.** Why must v_k be smooth? Oracle calls per relation, number of relations, linear algebra;
   the optimisation (exact arithmetic); the claim that single intervals cannot beat exponent 1/3 (and in exactly which
   model). Is the relation u − v ≡ 0 with u, v both smooth the same kind of relation as Dixon/QS (where one side is a
   square x² with x a free random element)? If not, what kind is it, and what does that change? Is 1.2114 consistent
   with the general "size principle" (cost governed by the log-size of the parts that must be smooth by chance)? Why does
   the claim beat SNFS (1.526) and the factory (1.6386) without special form, and exactly what is the oracle giving away
   for free (what would it cost to implement it by sieving)?
2. **H_S for deterministic constructors.** If u_k is a deterministic function of kN (examples: the smallest y-smooth
   number ≥ kN; the largest ≤ kN + h; a uniformly random one in I_k; a "cofactor scan" that tries smooth u₀ ≤ h in
   increasing order and returns u₀·⌈kN/u₀⌉ when that cofactor is smooth), is v_k still "random-like"? Identify every
   correlation between v_k and u_k or k you can (local/congruence structure, size distribution, gaps between smooth
   numbers), say whether each matters at the L-scale or only at toy scale, and what restriction on the set of k used is
   needed (e.g. what goes wrong if k is restricted to odd k, to primes, or if k may be as large as N).
3. **R: fruitfulness and rank for deterministic relation sources.** Can all dependencies be trivial (X ≡ ±Y)? Give the
   precise algebraic criterion for a set of such relations to yield a nontrivial square root of 1 (in terms of the
   relation lattice Λ_N = {e ∈ Z^{π(y)} : ∏ℓ^{e_ℓ} ≡ 1 mod N} and/or a "fruitfulness functional" on the F₂-dependency
   space). Is that functional computable from N alone for some N (e.g. Blum N), and what does that imply (both
   degeneracy examples and possible fixes)? Construct, if possible, an explicit constructor that sees only kN and N yet
   makes every dependency trivial. Which randomisations (random k, random choice among smooth numbers in I_k, random
   multipliers, Dixon-style hidden square roots, …) restore R, what exactly do they reduce R to, and what do they cost?
   Is a Dixon/Rabin-style hidden-root argument available within this oracle interface at all?
4. **Duplicates and the domain of k.** Can different k give the same relation (after removing common factors)? How
   should k be chosen?
5. **Verdict.** Is 1.2114 a correct heuristic implication? List exactly what the constructor must satisfy.
