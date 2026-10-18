# User-supplied idea (2026-10-18, unvalidated): CRT sign-splitting of residue branches in Harvey-style search

- Take M = m_1⋯m_t, a product of small distinct primes. A weight pair (a,b) with Nba^{-1} a QR mod every m_i serves 2^t residues r(ε) = CRT(ε_i r_i).
- The mod-M² offset V(r) = a r + b s + λ(N − r s), with s ≡ N r^{-1} and λ ≡ b r^{-1}, satisfies V(r(ε)) ≡ Σ ε_i Δ_i (mod M²), a signed subset sum. So meet-in-the-middle can test the 2^t branches in about 2^{t/2} steps.
- Orchestrator check (by hand): V is representative-independent given a r² ≡ bN (mod M); V mod m_i² depends only on r mod m_i; and V_i(−r_i) ≡ −V_i(r_i). So the separability identity holds.
- Orchestrator assessment: the gain is capped by 2^{ω(M)/2} = exp(O(log M / log log M)) = N^{o(1)}. It attacks only the per-branch "+T" setup term, so it cannot change the 1/5 exponent. It could still give N^{1/5}·exp(−Ω(log N / log log N)), which would be better than Harvey–Hittmeir 2022's (log log N)^{3/5} saving. Unverified; possible prior art is Hittmeir 2018-style CRT sieving.
