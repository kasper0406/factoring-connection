# Why a primality test does not tell you the factors (explainer)

*Agent `primality`, 2026-10-18. The proofs, exact formulas, data and reproduction commands are in
`notes/primality_analysis.md`.*

## 1. What a primality test actually checks

Take Fermat's test. Pick a base `a` and compute `y = a^(N−1) mod N`. If `N` is prime, `y = 1`. If `y ≠ 1`, then `N` is
certainly composite.

Now suppose `N = p·q`. By the Chinese remainder theorem, arithmetic mod `N` is arithmetic mod `p` and mod `q` side by
side. The number `y` is really the pair `(y mod p, y mod q)`. "`y ≠ 1`" means **at least one** of the two components is
not 1.

A factor would come from `gcd(y − 1, N)`. That gcd equals `p` exactly when `y ≡ 1 (mod p)` but `y ≢ 1 (mod q)`, i.e. when
the test fails in **one** component and not the other.

The test needs only "fails somewhere" (an OR). A factor needs "fails here but not there" (an XOR). For a random base the
failure is almost always in both components at once. Then the pair `(≠1, ≠1)` carries no visible label saying which
prime is which.

When does a component give 1? `a^(N−1) ≡ 1 (mod p)` exactly when the multiplicative order of `a` mod `p` divides `N − 1`.
That order divides `p − 1`, and `N − 1 ≡ q − 1 (mod p − 1)`. So it must divide `g = gcd(p−1, q−1)`, which is usually a
small number. A random base has this property with probability `g/(p−1)`. So the factor-revealing event has probability
exactly

`P = g/(p−1) + g/(q−1) − 2g²/((p−1)(q−1)) ≈ 2g/√N`.

* For a 32-bit `N` this is about `1/8 000`.
* For a 2048-bit RSA modulus it is about `2^(−1020)`.

The same structure holds for Miller–Rabin, Solovay–Strassen, Lucas and Frobenius tests. Each has its own group in place of
"numbers mod p": `p − 1` becomes `p + 1` for Lucas, and so on. Every one leaks only through such a **group-order
coincidence**. We derived and checked the exact probabilities for each test against every base of every semiprime below
20 000, and they match.

## 2. Worked example: N = 10403 = 101 · 103

Here `g = gcd(100, 102) = 2`.
* **Base `a = 2`.** `2^10402 mod 10403 = 9296`, so `N` is composite. But `gcd(9295, 10403) = 1`: no factor. The order of 2
  is 100 mod 101 and 51 mod 103; neither divides 10402. Both components fail.
* **Base `a = 100`.** Here `100 ≡ −1 (mod 101)`, which has order 2, and 2 divides 10402. Mod 103 the order is 17, which does
  not divide 10402. So `100^10402 ≡ 1 (mod 101)` only, and `gcd(100^10402 − 1, 10403) = 101`. A factor!

The exact probability of a lucky base here is `33/850 ≈ 3.9%`. It is large only because `N` is tiny.

## 3. Worked example: N = 561 = 3 · 11 · 17 (a Carmichael number)

* **Fermat.** For 561, `a^560 ≡ 1` for **every** base coprime to 561. Fermat's test is fooled completely: there is no
  failure and no factor.
* **Miller–Rabin.** It writes `560 = 2^4 · 35` and walks up the ladder `a^35, a^70, a^140, a^280, a^560`. For `a = 2` this
  is `263, 166, 67, 1, 1`. The value 67 squares to 1 but is not ±1, so the components reached 1 at *different* steps.
  Then `gcd(67 − 1, 561) = 33 = 3 · 11`. Here the coincidence is guaranteed (every component order divides 560). The only
  question is whether the components hit 1 at the same step, and they do not for 31 of every 32 bases.

So MR famously factors Carmichael numbers; on ordinary semiprimes the same mechanism fires with probability
`≈ (2/3)·g²/N`.

## 4. The factors *are* in the output, as a discrete logarithm

Since `φ(N) = N + 1 − (p + q)` and `a^φ(N) = 1`:

`a^(N−1) = a^(p+q−2) (mod N)`.

The failed test's output is `a` raised to the **unknown** exponent `s − 2`, where `s = p + q` lies between `2√N` and
`3√N`. Finding the exponent is a discrete logarithm in an interval of width about `√N`. Pollard's kangaroo (or
baby-step/giant-step) does that in about `√(√N) = N^(1/4)` steps. Then `p` and `q` are the roots of `X² − sX + N`.

For `N = 10403`, `a = 2`: from `y = 9296` alone the kangaroo found the exponent 202 in 65 multiplications. So
`s = 204` and `204² − 4·10403 = 4`, giving `101 · 103`.

We ran this on random semiprimes from 16 to 68 bits. Every run succeeded, and the cost grows like `N^0.24`. Lucas tests
are similar: their output is `z^(±(q−p)+2)`, which hides `q − p` instead, also at cost `N^(1/4)`.

So the honest answer to "why no insight?" is: **the insight is there, but encrypted as a discrete logarithm**. The cheapest
known generic way to decrypt it costs `N^(1/4)`: exponential in the number of digits. This is Pollard's 1978 idea and the
starting point of the best deterministic methods (Hittmeir; Harvey's `N^(1/5)`). It is a baseline, not a breakthrough.

## 5. Can we iterate on the residue, or combine checks?

* **Iterating.** Raising `y` to `N − 1` again gives `a^((s−2)²)`, and other combinations give `a^f(s)` for a known
  polynomial `f`. The unknown is still the single number `s`, but now in a much wider exponent range.
  * The kangaroo trick needs `f` linear, so this is strictly worse.
  * Measured: cost grows like `N^0.52`, and the squared form never returned `s` itself. It only produced
    order-finding collisions, a different and costlier route.
  * Combining the Fermat and Lucas outputs gives `s` and `q − p`, but these determine each other
    (`s² − (q−p)² = 4N`). We found no way to make the two searches help each other.
* **Better exponents or groups than `a^(N−1)`.** Any exponent built from `N` (such as `N+1`, `N²−1`, `N⁶−1`), in any group
  (quadratic tori, extension fields, elliptic curves), reduces modulo the hidden group order to a polynomial in `q`.
  * A leak then needs a coincidence between two independent primes, probability about `p^(−1)` per unit of group size.
    Bigger groups make it worse: about `2^(−186)` at 64 bits in `F_(p^6)`.
  * The exponents that *do* work are the highly composite ones, `lcm(1..B)`. They succeed when the hidden group order is
    smooth: Pollard `p−1`, Williams `p+1` and, re-rolling the group, Lenstra's elliptic-curve method.
  * Mixing in `N − 1` adds nothing.
  * The two ways known to genuinely "combine" are prior art:
    * Harvey's weighted collisions `a^(Aq+Bp) ≡ a^(AN+B) (mod p)`. They hold mod `p` but generally not mod `q`, so they
      create exactly the one-sided coincidence a gcd needs. Cost `N^(1/5)`.
    * Congruence-of-squares sieves (QS, NFS). They combine thousands of partial facts into one asymmetric event and are
      the only subexponential methods.

## 6. Partial information, and the "sum of two squares" bit

The project's open question is whether one can cheaply tell if `N` is a sum of two squares (`ε = +1`, both primes
`≡ 1 mod 4`).
* Miller–Rabin occasionally sees `a^(m·2^i) ≡ −1` with `i ≥ 1`. Then the previous value is a square root of −1 mod `N`.
  That *proves* `ε = +1`, and it even gives `N = u² + v²`.
  * But this happens with probability about `g²/(3N)` per base: `6·10⁻⁷` at 32 bits and `6·10⁻¹⁵` at 64 bits.
  * A continued-fraction search for `√N` certifies the same bit far more often for the same effort, and that search
    costs about `N^(1/4)` too.
* A pre-registered test of 14 transcript features (32–64 bits, about 20 000 numbers per size) found **no** correlation with
  `ε` down to `|ρ| ≈ 0.03`. Its two controls (one natural, one planted) were detected strongly.
* There is a structural reason. Reading 2-adic digits of `s` out of `a^(s−2)` needs elements of known 2-power order. Mod `N`
  the only such elements are `±1`, square roots of −1 (whose existence *is* the answer `ε = +1`), or elements that already
  factor `N`.

## 7. Why it has to be this way

A Fermat or Miller–Rabin witness can be produced from `N` alone: pick a random `a`, at least 3/4 of them work. So if any
fast procedure turned witnesses into factors, picking a random `a` and running it would itself be a fast factoring
algorithm. Recursion (split, then test and split the pieces) would then factor every number completely.

Compositeness proofs that are easy to *find* therefore cannot be easy to *decode* into factors, unless factoring is easy.
The factor `p` is also a short proof of compositeness, but nobody knows how to find it quickly. A Fermat witness proves "some
prime component misbehaves" without saying which, and the CRT structure is exactly what lets it do that.
