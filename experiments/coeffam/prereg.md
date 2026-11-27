# coeffam pre-registration (2026-11-27, before any coeffam experiment)

Goal: a poly(log N) evaluator for one coefficient at N = pq from a uniform FAMILY that yields > O(1) bits
about {p, q}. Semiprimes: balanced, p < q odd primes, sampled with the factorisation in the HARNESS only.
Evaluators under test take N (and the family parameter) only.

## E1. Level-1 mod-2 family {Delta^k mod 2 : k odd}
Delta = sum_{m odd} q^{m^2} mod 2, so a(N, Delta^k) = parity of #{N = sum_{i in bits(k)} 2^i x_i^2, x_i odd > 0}.
Nicolas-Serre height h(k) = n3(k) + n5(k) (bit 2j-1 -> n3 += 2^{j-1}, bit 2j -> n5 += 2^{j-1}).
- P1a: h(k) <= 1  =>  a(pq, Delta^k) = 0 for all good pq.
- P1b: h(k) = 2 (k in {7, 9, 17}) => a(pq, Delta^k) is a function of (p mod 8, q mod 8); on the public
  stratum N = 1 mod 8, Delta^9 decides whether p = 3 mod 4 (a Blum-integer-recognition bit).
- P1c: some k < 2^12 has a(pq, Delta^k) NOT a function of (p mod 2^6, q mod 2^6) (non-abelian 2-adic data,
  e.g. quartic character of 2 / x^2+32y^2).
- Ledger: conditional entropy of the vector (a(pq, Delta^k))_{k<K} given public data (N mod 2^J), as a
  function of K; and the subset of that entropy explained by (p mod 2^J, q mod 2^J).
- Evaluator cost: best factor-blind evaluator implemented (enumeration over the last variable), time vs log N.

## E2. Binary theta family x^2 + D y^2 (D varying, incl. D = 2^m), cost/information trade-off
Enumeration evaluator costs ~ sqrt(N/D). Prediction P2: the fraction of N on which the coefficient differs
from its public-data baseline is <= C sqrt(N/D)/sqrt(N) * polylog, i.e. the cheap members (D > N^{1-eps})
carry ~N^{-(1-eps)/2} bits; a polylog member gives o(1) bits.

## E3. Exponent-N Hasse-Witt / Dwork truncation (direction 2)
HW_N(f) := [x^{N-1}] f(x)^{(N-1)/2} mod N for f = x^3 + a x + b.  P3: HW_N = H_p * C_q (mod p), with
H_p the Hasse invariant (= a_p mod p) and C_q := [x^{q-1}] f^{(q-1)/2} in Z (a non-Frobenius integer).
Same for the Legendre-family truncation F_{N-1}(lambda) = sum_{k<N} binom(2k,k)^2 (lambda/16)^k.
Native cost: O(N) direct, O(N^{1/2}) baby-step/giant-step (holonomic). Kill if no sub-power method.

## E4. Own idea: sums-of-three-squares / class-number family via orbit (Smith-type) recursion
r_3(N) = 12 H(4N) (N = 1,2 mod 4).  P4a: r_3(N) mod 48 is determined by r_2(N) and r_{x^2+2y^2}(N)
(points fixed by subgroups of the order-48 hyperoctahedral group), so its hidden bits are (-1/p), (-2/p).
P4b: h(-4N) mod 8 (Redei) is a function of (-1/p), (2/p), (p/q); family {h(-4Nm) mod 8 : m <= M} carries
~ pi(M) quadratic-character bits of p. Native evaluator: none polylog expected (class group L[1/2]).

## Kill criterion (all)
Every evaluator found costs a power of N (or L[1/2]) or needs a factorisation-dependent local datum
(Frobenius / Legendre symbol at p).  Positive claim requires a factor-blind implementation given only N,
checked on >= 35-bit inputs.
