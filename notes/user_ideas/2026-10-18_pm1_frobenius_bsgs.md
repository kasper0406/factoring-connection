# User-supplied idea (2026-10-18): "Pollard p−1 + Frobenius-BSGS hybrid"

Pasted by the user, verbatim in substance:

- Run Pollard p−1 stage 1: b = a^E mod N, E = lcm(1..B). If p−1 = S·R with S | E (S = B-smooth part), then ord_p(b) | R = (p−1)/S.
- Since R | p−1, N = pq ≡ q (mod R), so b^N ≡ b^{N mod R} ≡ b^x (mod p) for some x ∈ [0, R).
- BSGS (or Pollard–Strassen multipoint evaluation) over x ∈ [0, √N/Σ], looking for gcd(b^N − b^x, N) = p. Cost O(N^{1/4}/√Σ) when S ≥ Σ, regardless of whether R is prime or composite; standard stage 2 only handles R prime ≤ B2.
- (Earlier part of the note: p+q ≡ N+1 (mod R) is trivial; the Lucas torus analogue gives nothing new by itself.)
