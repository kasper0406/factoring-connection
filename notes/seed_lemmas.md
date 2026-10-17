# Seed lemmas: independent proofs and checks (agent infraA)

Notation: `N > 1`, `R_N = Z/NZ`, `F_N(X) = (1+X)^N - 1 - X^N`, `a_k = C(N,k) mod N`.
Since the endpoint terms cancel exactly, `F_N = sum_{0<k<N} a_k X^k` in `R_N[X]`.
`p` = least prime factor of `N`. For a prime `l`, `j <=_l m` means every base-`l`
digit of `j` is at most the corresponding digit of `m` ("digit-dominated").

Verification command (all checks below, all composites `N <= 1000`, 2.3 s):

```
.venv/bin/python src/seed_checks.py            # BRIEF §9 smoke test + check_extended(1000)
.venv/bin/python -m pytest -q tests/test_seed_checks.py
```

Output: `Seed checks passed for 153 small composites.` and
`{'composites': 831, 'semiprimes_pq': 288, 'prime_powers': 25, 'three_factor': 135, 'scalar_zero_composites': [561]}`.

---

## L1 (§2.1) First surviving coefficient is `N/p` at index `p` — **proved**, checked N ≤ 1000

*Proof.* For `1 <= k < p`, every integer in `[1,k]` has all prime factors `< p`,
so `k!` is a unit mod `N`. The integer identity `k! C(N,k) = N(N-1)...(N-k+1)`
has right side `≡ 0 (mod N)`; multiplying by `(k!)^{-1}` gives `a_k = 0`.
For `k = p`: `C(N,p) = (N/p) C(N-1,p-1)` (integer identity). Also
`(p-1)! C(N-1,p-1) = prod_{j=1}^{p-1} (N-j) ≡ (-1)^{p-1} (p-1)! (mod N)`, and
`(p-1)!` is a unit, so `C(N-1,p-1) = (-1)^{p-1} + N t` for an integer `t`.
Hence `C(N,p) = (-1)^{p-1} N/p + (N/p) N t ≡ (-1)^{p-1} N/p (mod N)`. For odd `p`
this is `N/p`; for `p = 2`, `-N/2 ≡ N/2`. As `0 < N/p < N`, `a_p = N/p ≠ 0`. ∎
(Valid for all composite `N`, including prime powers; `p < N` so the removed
`X^N` is irrelevant.)

## L2 (§2.2) `k a_k = 0`; any support index gives a factor — **proved**, checked

*Proof.* `k C(N,k) = N C(N-1,k-1)`, so `k a_k ≡ 0 (mod N)`. If `gcd(k,N) = 1`,
`k` is a unit and `a_k = 0`. For `0 < k < N`, `gcd(k,N) ≠ N`. So `a_k ≠ 0`
implies `1 < gcd(k,N) < N`. ∎

## L3 (§2.4) Prefix test ⇔ `K <= p` ⇔ `(K-1)!` is a unit — **proved**, checked for all `1 <= K <= N`, N ≤ 1000

*Proof.* `F_N ≡ 0 (mod N, X^K)` iff `a_k = 0` for all `0 < k < K` (the constant
coefficient is 0). For composite `N`, by L1 the least support index is `p`, so
this holds iff `K <= p`. And `gcd(N,(K-1)!) = 1` iff no prime `<= K-1` divides
`N` iff `p > K-1` iff `K <= p`. ∎ (For prime `N`, `F_N = 0` and the factorial
side holds iff `K <= N`, so the equivalence holds for `K <= N`.)

## L4 (§2.5) `F_N' = 0` in `R_N[X]`; Hasse derivatives keep coefficients — **proved**, checked

*Proof.* Formal differentiation is coefficientwise: the coefficient of `X^{k-1}`
in `F_N'` is `k a_k`, zero by L2; the removed endpoints contribute `N X^{N-1}`
and `-N X^{N-1}`, both 0 mod `N`. So every ordinary derivative of order
`>= 1` vanishes. The `k`-th Hasse derivative at 0 is
`D^{[k]}F_N(0) = [X^k]F_N = a_k` (it is defined by the integer binomials
`C(i,k)` without division), so it retains the coefficients, but querying one
order `k` at a time is subject to §2.3 (a gcd with `k` answers it). ∎

## L5 (§2.6) The 561 scalar-zero trap, and its exact scope — **proved**, checked

*Claim.* `F_N(a) = 0` in `R_N` for every `a ∈ R_N` iff `a^N = a` for all `a`
iff `N` is prime or a Carmichael number (squarefree, `l-1 | N-1` for all `l | N`).

*Proof.* Let `g(a) = a^N - a` (a function on `R_N`). Then
`F_N(a) = (1+a)^N - 1 - a^N = g(a+1) - g(a)`. So `F_N ≡ 0` as a function iff `g`
is invariant under `a ↦ a+1`, i.e. constant, i.e. `g ≡ g(0) = 0`. The condition
`a^N ≡ a` for all `a` is Korselt's criterion (squarefree, and `l - 1 | N - 1`
for each prime `l | N`), whose composite solutions are exactly the Carmichael
numbers. For `561 = 3·11·17`: `2, 10, 16 | 560`. Meanwhile `F_561 ≠ 0` formally
(L1: `a_3 = 187`). ∎
Check: the only composite `N <= 1000` with `F_N` vanishing on all of `R_N` is
561 (`scalar_zero_composites = [561]`), consistent with the claim.

*Consequence for zero testing.* Scalar evaluation in `R_N` is not a sound
identity test for any Carmichael `N`. More generally, evaluating at elements of
an `R_N`-algebra `A` sees `F_N` only through the ideal of polynomials vanishing on `A`;
L5 shows this ideal contains `F_N` for `A = R_N`, `N` Carmichael.

## L6 (§2.6) Support lower bound: `|supp F_N| >= p-1` for `N = pq`, `p < q` — **proved**, checked

*Proof.* Mod `q`, Frobenius gives `(1+X)^{pq} = ((1+X)^q)^p ≡ (1+X^q)^p`, so
`F_N ≡ sum_{j=1}^{p-1} C(p,j) X^{qj} (mod q)`. Each `C(p,j)` (`0<j<p`) is a
positive integer whose prime factors are `<= p < q`, hence nonzero mod `q`, so
`a_{qj} ≠ 0` for the `p-1` indices `qj`. ∎ (Exact count in L7.)

## L7 Lucas description of the support — **proved**, checked for all N ≤ 1000

**(a) Squarefree prime.** Let `l` be a prime with `l || N`, `N = l m`, `l ∤ m`.
Then for `0 < k < N`:
`C(N,k) ≢ 0 (mod l)` iff `k = l j` with `j <=_l m`, and then `C(N, l j) ≡ C(m, j) (mod l)`.

*Proof.* The base-`l` digits of `N` are `(0, m_0, m_1, ...)`. Lucas:
`C(N,k) ≡ prod_i C(n_i, k_i) (mod l)` with `C(n,k) = 0` for `k > n`, and each
`C(n_i,k_i)` with `k_i <= n_i < l` is nonzero mod `l`. So the product is nonzero
iff `k_0 <= 0` and `k_{i+1} <= m_i` for all `i`, i.e. `k = l j`, `j <=_l m`;
dropping the `C(0,0)=1` factor gives `C(N,lj) ≡ C(m,j)`. ∎

**(b) `N = pq`, `p < q` primes.** By CRT, `a_k ≠ 0` iff `C(N,k) ≢ 0` mod `p` or mod `q`.
* mod-`p` support: `S_p = {p j : 0 < j < q, j <=_p q}`;
* mod-`q` support: `q` dominates exactly `0..p` in base `q` (single digit `p`), so `S_q = {q j : 1 <= j <= p-1}`;
* `S_p ∩ S_q = ∅` (an index in both is divisible by `pq = N`).

Hence `supp F_N = S_p ⊔ S_q` and **`|supp F_N| = (p-1) + prod_i (q_i + 1) - 2`** where
`q_i` are the base-`p` digits of `q`. (Checked equal for all 288 such `N <= 1000`.)
For `p < q < p^2`: `q = q_1 p + q_0` with `1 <= q_1 <= p-1`, `1 <= q_0 <= p-1`
(`q_0 ≠ 0` since `p ∤ q`), so `|S_p| = (q_1+1)(q_0+1) - 2`. For **same-bit-length**
primes (`p < q < 2p`, the balanced family) `q_1 = 1`, `q_0 = q - p`, hence
`S_p = {p j : 1 <= j <= q-p} ∪ {p j : p <= j <= q-1}` (the two "blocks"),
`|S_p| = 2(q-p)` and **`|supp F_N| = 2q - p - 1`** (e.g. `N = 143`: 14 indices).
So the support has size `Θ(p)`, density `≈ 2/p` in `[1, N)`: exponential size, exponentially small density.

**(c) Prime powers `N = l^e`.** For `0 < k < l^e` with `t = v_l(k)`:
`v_l(C(l^e,k)) = e - t`. *Proof:* `C(l^e,k) = (l^e/k) C(l^e - 1, k - 1)` and
`l^e - 1` has all base-`l` digits `l-1`, so by Lucas `C(l^e-1,k-1) ≢ 0 (mod l)`;
the factor `l^e/k` has valuation `e - t`. ∎ Hence `a_k ≠ 0` iff `t >= 1`:
**`supp F_{l^e} = {l, 2l, ..., l^e - l}`** (`l^{e-1} - 1` terms; checked for all 25
prime powers `<= 1000`).

**(d) General `N`, Kummer form.** `v_l(C(N,k))` = number of carries when adding
`k` and `N-k` in base `l` (Kummer). With `N = prod l^{e_l}`,
**`a_k ≠ 0` iff there is a prime `l` with `carries_l(k, N-k) < e_l`.**
(Kummer + CRT; checked as `support == kummer_support(N)` for all 831 composites.)
For squarefree `N` this specialises to `supp F_N = ⋃_{l | N} {l j : j <=_l N/l}`
(interior indices), which for three distinct primes `p<q<r` is the union of
three Lucas sets (not disjoint in general: `k` divisible by `pq` can lie in two
of them). Checked for all 135 three-prime `N <= 1000`.

---
Code: `src/seed_checks.py` (`check_extended`, `kummer_support`,
`lucas_support_squarefree_part`, `dominated`, `carries`). Reference code only
(trial division, explicit binomials).
