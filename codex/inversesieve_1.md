# inversesieve_1 — blind round (GPT-6.1-sol, read-only on the repo)

You are a research referee and collaborator on a project on deterministic integer factorization
(repository: the current repo root, read-only; scratch writes go to your cwd). Read first:
`codex/r20_wild.out` §2 (the idea under test), `notes/dchworst_results.md` §§1–5 (Theorems W, R,
barriers B1–B4), `notes/coverage38_results.md` §§0–5 (Theorem U, Prop L, Lehman dichotomy, Prop P),
`notes/below16_results.md` §0 (T1/T2 local-selector floors) and `notes/specialfam_results.md` §0
(Theorem G, Coppersmith containers).

## Setting

N = pq, p < q ≤ 2p primes, P = √N. A DCL design is a band 𝒫 of J primes ℓ ≡ 3 (mod 4) in (Y, 2Y]
and a weight set V of primitive (A, B), height ≤ H (H² < Y/2 so classes are distinct). The pair is
**covered** iff ℓ | Bp − Aq for some (A,B) ∈ V, ℓ ∈ 𝒫. Theorem R extracts p from any covered pair in
the list cost. The design cost is X ≈ HP/Y² + JH²; with λ := |V|·Σ 1/ℓ ≍ JH²/Y = J^ϑ the balanced
cost is X = N^{(1+ϑ)/(6+ϑ)}; ϑ = 1/5 gives N^{6/31}. Theorem U says no design in the family has
worst-case coverage at cost < N^{1/5} unless λ < J^{1/4}; no one can prove coverage at λ = J^{1/5}.

## The idea to evaluate (r20 #2)

Run the design at ϑ = 1/5 and, instead of proving every pair covered, prove an **inverse / structure
theorem**: if the true pair is uncovered by all band primes, the forbidden-residue pattern forces
structure that can be turned, from N alone in N^{6/31+o(1)} work, into ≤ N^{6/31+o(1)} explicit
Coppersmith containers (intervals of length N^{1/4−ε}, or residue classes p ≡ r mod M with
M ≥ N^{1/4+ε}, or Theorem-G polynomial boxes of volume ≤ N^{12/31}) one of which contains p. Kill
tests: an argument that merely rediscovers pq = N is useless; a container needing N^{1/5} work to
enumerate is useless.

Per prime: uncovered ⇔ t_ℓ = p q^{-1} mod ℓ ∉ F_H(ℓ) := {A B^{-1} mod ℓ : (A,B) ∈ V}, equivalently
p mod ℓ ∉ R_ℓ(N) := {x : B x² ≡ A N (mod ℓ) for some (A,B) ∈ V}. U(N) := {x ∈ [P/2, P] : x mod ℓ ∉
R_ℓ(N) for all ℓ ∈ 𝒫} is the set of all integer candidates left uncovered.

## What I want from you (work blind: do not search the repo for an `inversesieve` note — there is none yet)

1. Is there ANY inverse-sieve theorem in the literature (Green–Harper arXiv 1311.6176; Helfgott–
   Venkatesh; Walsh's inverse sieve in high dimensions; Croot–Elsholtz; Gallagher's larger sieve;
   "large sieve sharp only for quadratic sets") whose hypotheses are met by U(N) or by the single
   point p in the regime ρ_ℓ := |R_ℓ(N)|/ℓ ≍ λ/J → 0 (few classes removed per prime), |U(N)| tiny?
   State hypotheses precisely and say yes/no for each, with reasons.
2. Information/locality: the uncovered event has probability ≈ e^{−λ} = e^{−J^{1/5}}, but each prime
   contributes only ρ_ℓ ≈ λ/J. Coppersmith can absorb only a modulus M ≤ P (≤ log P / log Y band
   primes). Prove or refute: any container family that depends on p only through p mod M with M a
   product of band primes (CRT boxes) plus a size window needs ≥ N^{1/4}(1−o(1)) containers.
   Extend, if you can, to Theorem-G polynomial boxes or to containers built from any bounded number
   of band primes. Is there a cleaner, stronger "diffuseness" barrier?
3. Under the Poisson/random model at ϑ = 1/5, what is E Σ_N |U(N)| and E #{uncovered pairs}? What
   does that mean for the logical content of the inverse theorem (is it equivalent to a coverage
   theorem for pairs outside the containers)? Compare quantitatively with Theorem D′ / large-sieve
   bounds: |U(N)| ≤ (P + 4Y²)/λ. What exceptional-set bound would an inverse theorem imply, and is
   it beyond the large sieve?
4. Try hard to find a NON-trivial positive structure statement: something provable about uncovered
   pairs that is stronger than pq = N and does not need N^{1/5} to enumerate (e.g. Diophantine
   structure of q/p, continued fractions, Lehman-easiness, a polynomial-value container, CRT
   forcing families). Give your best candidate and its kill-test status. If you think none exists,
   give the sharpest precise reason (e.g. a pseudorandomness statement with respect to every
   container family enumerable in N^{1/5}, in a clearly defined model), and say what you can prove.
5. Suggest the 3–5 most informative small-n experiments (n = 24–40 bits) on exact U(N), with
   concrete statistics and null models, that would distinguish "structured" from "pseudorandom".

Label every claim PROVED (with proof), KNOWN (with reference), HEURISTIC, or UNVERIFIED. Be
adversarial and quantitative. Plain text / markdown output, ≤ 4000 words.
