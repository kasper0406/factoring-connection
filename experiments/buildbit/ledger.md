# Hidden-class bit yield ledger

No public endpoint passed the preregistered non-decay test. The constructive results are a one-sided T0 algorithm with verifiable cross-multiplier certificates and a uniform reduction from full T1 factor recovery to one Boolean central-binomial test. The latter works with a factor-backed oracle. Its factor-free zero test remains open.

## R1. Signed relations across continued fractions

**Built and proved.** For convergent numerator and denominator H and K just before complete denominator Q at step j of √(kN), H²−kNK²=(−1)^j Q. For each N the algorithm collects 16 denominators at each of six fixed and six adaptively selected multipliers. Every chosen Q has a modular square root of its signed value. Split these positive integers into pairwise coprime atoms using gcds. For an atom that is a nonsquare, use its exponent parity as a column over F₂. Ignore square atoms. Coprimality makes this an exact rational square-class matrix even when the atoms are composite. Its nullspace contains every subset whose product is a square.

For a null vector with product v² and odd sign, the product of the corresponding H divided by v modulo N squares to −1. Gcds treat noninvertible denominators and incompatible modular square roots as explicit splits. Arbitrarily many polynomially bounded steps and multipliers preserve polynomial bit complexity with arbitrary-precision arithmetic. The native word-size helper is sized for the registered measurements.

**Measured.** Each size has 2,048 products per class, with both factors of half the stated bit length and the product of exactly the stated bit length. In the table, all certificate columns count split-clean inputs with both primes 1 modulo 4. The denominator column is the number of such inputs. The final column includes both classes. “Sum” counts the total individual-multiplier hits, “union” counts the inputs those hits covered.

| Bits | Denominator | First six pooled | Individual sum | Individual union | Multiple individual hits | Full pool | Proper gcd splits |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 1700 | 101 | 177 | 174 | 3 | 318 | 635 |
| 52 | 2012 | 12 | 25 | 25 | 0 | 26 | 57 |
| 64 | 2045 | 2 | 3 | 3 | 0 | 3 | 7 |
| 80 | 2048 | 1 | 1 | 1 | 0 | 1 | 0 |

One joint certificate for N=647442828917 uses Q=517873 at (k,j)=(1,1) and Q=4660857=9·517873 at (47,16). Its verified root is 467763698354. At 80 bits the one certificate used a single Q at (2,1). A fair fallback gives balanced advantage equal to one quarter of clean sensitivity. Simultaneous 99% intervals, in percentage points, are 4.676 [3.886,5.545], 0.3231 [0.1475,0.6037], 0.03667 [0.00144,0.1803], and 0.01221 [0.0000031,0.1310]. The effective endpoint exponent for advantage proportional to N^(−β) is 0.215 [0.122,0.520]. Only one 80-bit hit contributes to that estimate. All certificates were consistent with T0.

**Obstruction.** The extra joint yields are 144, 1, 0 and 0. Thus independent repetitions would misdescribe the small-size gain. At larger sizes this bank lacks square relations of odd sign, an M/Q relation-supply gap. Prescribing Q₁=Q₊u² after seeing a positive signed Q₊=H² mod N is already equivalent to constructing a root of −1 when Hu is invertible. Indeed a₀²=−Q₊u² mod N. Conversely choose a representative a₀ of a root of that congruence sufficiently large, and k=(a₀²+Q₊u²)/N realizes the prescribed first denominator. This addresses deliberate matching without ruling out other adaptive statistics. Even access to the parity of the entire period only gives a one-sided T0 test. For example 221=13·17 has period six.

## R2. A public ambiguous class with only the target genus characters

**Built and proved using Gauss genus theory.** Set d=4N−1, D=−Nd and f=(N,N,(5N−1)/4). This primitive form is ambiguous and represents m=f(1,2)=8N−1, coprime to D. Since D is odd its assigned characters are the Legendre symbols at the odd primes dividing D, evaluated on m. Every prime dividing d sees m=1. At p and q the value is (−1/p) and (−1/q). The principal genus is the square subgroup, so [f] is a square exactly in the 1,1 case. This includes nonfundamental D. Construction uses only integer arithmetic on N. Enumeration of the square subgroup verified all 125 small pairs specified in the follow-up, including 26 nonfundamental discriminants. [Genus theorem](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/jlms.12921)

**Yield and obstruction.** The form is already in hand. Its order divides two in both classes. The product of its two remaining characters is always +1, an exact Q-family cancellation. For this pencil the question whether m is a square modulo |D| reduces by CRT to whether −1 is a square modulo N, since m=1 modulo d. Computing the single class's square membership by other means remains a technique gap. No hidden involution or class number is needed to construct the query.

## R3. A scalar quartic for the cubic norm-one square test

**Built and proved.** In the conductor-seven cubic algebra take b=a³/Norm(a) for a uniform unit a. Form T=Tr(b), S=Tr(b⁻¹), and H(X)=X⁴−2TX²−8X+T²−4S. Over an algebraic closure there are four norm-one choices of square roots of the three conjugates of b. Their traces are the roots of H. For a choice with trace u and inverse trace v, u²−2v=T and v²−2u=S. When H is separable its trace distinguishes all four choices, so a local scalar root is equivalent to a local norm-one square root. This gives an explicit univariate test without an algebraic group order.

At an inert prime the torus has odd order s²+s+1 and there is exactly one root. At a split prime the three square-class bits of a are independent. The map to b projects them onto the two-dimensional even-parity subspace, uniformly. All four choices are rational if b is square, otherwise none are. The probability that b is square at both split primes is 1/16. The four-root discriminant is always a square since the sign changes form a Klein four group and conjugates are cyclically permuted.

**Measured.** At 40, 60, 80 bits the factor-backed square predicate succeeded on 512/512 inert inputs each time and on 39, 29, 32 of 512 split inputs. Publicly computed Res(H,X^N−X) was a unit on all 3,072 inputs. There were no degenerate quartics or proper gcds. Its Jacobi sign gave balanced advantages −0.20, +0.39, +2.15 percentage points. All simultaneous 99% intervals contain zero, so a decay exponent is unidentifiable. For a nonzero root the replacement tests ord(α) dividing pq−1, an M-family condition. The discriminant loses the distinction identically. The gap for the full torus method is recognizing a scalar root modulo N.

## R4. One central-binomial zero test through the entire lift

**Proved.** Put m=⌊hN/d⌋ and B=binom(2m,m). Writing N=st, v=ht mod d and e=hN mod d gives

    m = s·⌊ht/d⌋ + (sv−e)/d.

For d=2^j≥4, h odd, and s>d, the low base-s digit causes a carry in m+m exactly when v>d/2. When 1≤h≤d/4 and t<2s the upper digit alone cannot cause a carry. Kummer's theorem identifies divisibility of B by s with at least one carry. For unbalanced inputs a lower-half pair still prevents divisibility by the larger prime, and an upper-half pair forces both primes to divide. Small primes up to the fixed denominator can be tried directly. [Binomial carries](https://www.cecm.sfu.ca/organics/papers/granville/)

- **T0.** Take (h,d)=(1,4). Then N divides B exactly in the 3,3 case. Up to sign B is the coefficient of x^(N−1) in (x³−x)^((N−1)/2).
- **T1.** Use d=16 with h=1 if N=3 mod 16 and h=3 if N=11 mod 16. Multiplication by h puts both factor residues in the lower half for r=1 mod 16 and both in the upper half for r=9 mod 16. A single zero answer decides T1. Balanced inputs produce only zero or a unit.
- **T2.** Take (h,d)=(1,7) on balanced inputs. The same digit computation uses e=1. The split residue pairs (1,1),(6,6) give respectively a unit and zero. Both inert pairs (2,4),(3,5) give a proper gcd. All 3,072 carry-oracle checks agreed, with 1,536 oracle splits. Access to the coefficient or its gcd suffices here.

**Uniform T1 lift.** Suppose a=r mod d/2 is known for d=2^j, j≥4. Set c=N/a² mod d, so c=3 mod 8. The lower candidate a pairs with b=N/a mod d. Changing lifts adds d/2 to both candidates. Find an odd h<d/4 for which ha,hb are in the same half. To construct it, let L=d/8. By successive squaring, c^L=1+d/2 mod d. The upper-half bit of c^i XOR the parity of i differs at i=0 and i=L. Bisect this disagreement to adjacent indices i,i+1. The upper-half bits there agree. Start with h=c^i/a modulo d and fold by sign and by h ↦ d/2−h into (0,d/4). These transformations preserve agreement for two odd inputs. The direction at a is public. One Boolean divisibility query determines the actual lift.

The balanced implementation recovered r on 768/768 fresh 40-, 60- and 80-bit inputs with 16, 26 and 36 factor-backed queries respectively, followed by enumeration of a constant-size interval of candidates. No query produced a proper binomial gcd. In the general version, try the two possible lifts of each factor as exact divisors before every query. A miss ensures both primes exceed d. This proves full recovery for arbitrary T1 factor ratios in O(log N) Boolean queries with polynomial additional work. All 300 small pairs passed, including 168 unbalanced pairs. The general oracle discarded the proper-gcd information on 52 queries. This is a uniform coefficient reduction and gives no reduction from full lifting to a fixed T0 or fixed T1 oracle.

**Evaluation gap, measured.** The code computes h and m in polynomial bit time. The large-input Kummer oracle uses the factors. Direct small tests used exact integer binomials. The local valuation reads the floor of a scaled hidden cofactor, an A-type quantity. The public coefficient fold mod X^L−1 aliases the desired term with its neighbors. All 2,304 folds at L=8,16,32 returned units, including all 1,152 required-zero observations. One counterexample for every L is N=753071038451, r=945961, q=796091. The fold is refuted as an exact replacement. General native coefficient evaluation remains a technique gap.

**Single next construction, speculative.** For Uₙ=binom(2n,n)/4ⁿ retain the endpoint relation 2nUₙ=(2n−1)Uₙ₋₁ over Z/N, allowing arbitrary starting value. Over F_s the possible terminal values at m are forced to zero exactly in the upper half of the low base-s digit. This follows from the general formal solution (1−4z)^((s−1)/2)V(z^s) of (1−4z)U′=2U. On the aligned inputs even the support of this one coordinate decides the bit. Compose the two-endpoint linear relations on dyadic intervals, retaining the midpoint coordinate. At T0 the reflection n ↦ (N+1)/2−n interchanges the two recurrence coefficients up to sign. The outstanding step is a shift/doubling formula for the midpoint image ideal that avoids visiting every singular index. Ordinary interval-by-interval elimination takes exponentially many leaves.

## Reproduction

From this scratch directory, compile basis.cpp as a shared C++17 library named basis.so. With the repository virtualenv Python use samples.py, cf.py separately at 40,52,64,80, cubic_carry.py, summarize.py, lift.py, genus.py, general_lift.py and check_all.py in that order. The original registration is prereg.md, SHA-256 8bf2d71e3cf75f89125c175c2ea9fefff860661f91b15ec01701f77147fac429. The three subsequent proof-validation registrations are separate files. summary.json contains raw denominators and simultaneous intervals. The other counts files record oracle and validation outcomes. Each experimental row, selected multiplier list and extracted first root is retained. All 41,581 distinct sampled factors passed PARI's primality proof.
