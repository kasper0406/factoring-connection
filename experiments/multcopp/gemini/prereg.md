# Pre-registration: `multcopp` Toy Test (CM Invariant-Derivation Multiplicity vs. Coppersmith Barrier)

Frozen before writing or running `toy_test.py`.

## Setup
Consider the CM elliptic curve $E: w^2 = 1 - s^4$ over $\mathbb{Z}[1/2]$ (CM by $\mathbb{Z}[i]$, $K = \mathbb{Q}(i)$) with invariant derivation $D = w \partial_s - 2s^3 \partial_w$ (equivalently $s'' = -2s^3$, $s(0)=0, s'(0)=1$) and exact integer jet values $E_a = D^a s(0) \in \mathbb{Z}$ ($a = 4k+1$).

## Pre-registered Predictions & Pass/Fail Criteria

1. **P1 (Remark item 2 — $p$-adic multiplicity at inert primes vs. $a < p$ threshold):**
   For all inert primes $p \in \{7, 11, 19, 23, 31\}$ ($p \equiv 3 \pmod 4$) and all $a = 4k+1 \le 401$:
   - **P1a:** $p \mid E_a$ holds for **100%** of $a = 4k+1 \ge p$.
   - **P1b:** $v_p(E_a) = v_p(a!) - \lfloor (a-1)/(p^2-1) \rfloor$? More conservatively, test the two-sided asymptotic slope:
     for every tested $p$, $v_p(E_a) \ge \sum_{j \ge 0} \lfloor a / p^{2j+1} \rfloor - \lfloor \log_p a \rfloor - 1$ holds for **100%** of $a = 4k+1 \le 401$, while at $a = 401$ for $p=7$, $v_7(E_{401})$ is strictly less than $v_7(401!)$ (confirming loss of the $\sum_{j \ge 1} a/p^{2j}$ formal-log denominators and retention of $\frac{p}{p^2-1} a$).

2. **P2 (Split-prime non-vanishing via $D^q \equiv A_q D \pmod q$):**
   For all split primes $q \in \{5, 13, 17, 29, 37\}$ ($q \equiv 1 \pmod 4$):
   - Let $A_q = E_q \bmod q$. Then $A_q \equiv (-1)^{(q-1)/4} \binom{(q-1)/2}{(q-1)/4} \not\equiv 0 \pmod q$, and $E_{r + k(q-1)} \equiv A_q^k E_r \pmod q$ holds with **0 violations** across all $r \in \{1, 5, \dots, q-4\}$ and $k \ge 0$ with $r + k(q-1) \le 401$. In particular, $q \nmid E_{1 + k(q-1)}$ for **100%** of $k$.

3. **P3 (Archimedean Bost-slope ledger vs. single-prime Coppersmith requirement):**
   - **P3a:** $\frac{1}{a}\bigl(\log |E_a| - (a \log a - a)\bigr) \to -\log \varpi \approx -0.270812$ (within $0.02$ at $a = 401$), matching Lean's `phF(a) - a * logOm`.
   - **P3b:** For every single inert prime $p \in \{7, 11, 19, 23, 31\}$ and $1 < a = 4k+1 \le 401$, the Coppersmith single-prime ratio $\frac{v_p(E_a)\log p}{\log |E_a|} < 0.15 < 1$ (cannot force vanishing at a single hidden prime).

4. **P4 (Mod-$N$ hidden-factor test & the $a \ge p \asymp \sqrt{N}$ wall):**
   Generate 30 deterministic balanced semiprimes $N = pq$ (14–20 bits, seed 42) with $\left(\frac{-1}{N}\right) = -1$ ($p \equiv 3 \pmod 4$ inert, $q \equiv 1 \pmod 4$ split):
   - **P4a (Below-$p$ failure at Coppersmith scale $N^{1/4}$):** At $a_{\text{low}} = 4\lfloor N^{1/4}/4 \rfloor + 1 < p$, $\gcd(E_{a_{\text{low}}} \bmod N, N) = p$ succeeds on $\le 2/30$ instances (accidental $O(1/p)$ hits only; no systematic multiplicity).
   - **P4b (Above-$p$ split/inert separation at $\sqrt{N}$ scale):** At $a_{\text{high}} = 4\lceil \sqrt{2N}/4 \rceil + 1 > p$, $p \mid E_{a_{\text{high}}}$ holds on **30/30** instances, and $\gcd(E_{a_{\text{high}}} \bmod N, N) = p$ recovers $p$ on $\ge 25/30$ instances — but requires jet order $a_{\text{high}} > p \asymp \sqrt{N}$ (exponential cost $\ge \tilde{\Omega}(N^{1/2})$).
