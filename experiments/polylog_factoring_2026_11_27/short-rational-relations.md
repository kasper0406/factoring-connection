# What short rational relations can evaluate

For a unit target t modulo N, consider **all** fractions u/v with u≡tv mod N, gcd(v,N)=1 and 0<|uv|≤N. With a prime bound L polynomial in the bit length of N, we can decide in polynomial bit time whether an odd number of these fractions multiply to a signed L-smooth integer times a rational square. When they do, the evaluator returns a unit w and a signed L-smooth integer s with t w²≡s mod N. This is a public-input algorithm. The probability that it returns a relation on a random square target is the outstanding condition for using it in the existing all-input Dixon reduction.

## Exact construction

Reduce each fraction first. Its common divisor is a unit modulo N, and its removal changes the product uv by a square. Write a=(tv−u)/N. Then a/v approximates t/N within 1/v². The continued-fraction determinant bound puts (a,v) in a pair of adjacent convergent vectors with both integer coefficients of absolute value at most three. Enumerating those combinations, filtering by height and checking unit denominators gives the complete list in O(log N) candidates up to an absolute constant. The code includes the determinant argument and a brute-force comparison of complete lists. [Enumeration](check_rational_height_shell.py)

Strip the base primes from each |uv| and keep the remaining positive integers. Split overlapping cofactors by gcd until the surviving nonsquare atoms are pairwise coprime. Each split by g>1 reduces the log of the product of the work items by at least log g. Keeping binary row-incidence vectors through the splits gives a polynomial-size parity matrix. An odd subset has square cofactor product precisely when its row sum vanishes. Eliminate over the field with two elements and check the parity of the support of each kernel basis vector. If any kernel vector has odd support, at least one basis vector does.

For a selected subset of k rows write ∏uv=s z², stripping the square cofactor outside the base, with s signed and smooth. Since t^k(∏v)²≡s z², the known unit w=t^((k−1)/2)∏v/z mod N gives the promised relation. The public entry point is `evaluate_target(n, target, base_limit)`. It returns the row count, w and s, or None if this entire family has no witness. With O(log N) rows the integer products have polynomial bit length. Repetitions of one fraction add only squares. [Constructor and exact subset check](check_same_target_square_relations.py)

For a positive example, take N=247, t=23 and L=3. The fractions −17/10, −5/32 and −1/204 have rough parts 85, 5 and 17. The product of these three rough parts is 85². No single eligible fraction for this target has a square rough part. After absorbing the base square as well, the three rows certify 23·116²≡−3 mod 247.

## Matched random-square targets

The samples, seed strings and base cutoff L=b² are the same as in the [constructor baseline](https://github.com/kasper0406/factoring-connection/commit/ffa890d). Each cell counts targets out of 4,096, with 1,024 per modulus and four balanced semiprimes per bit size. The final two columns allow square cofactors.

| Bits b | Ordinary convergents, both sides smooth | All eligible fractions, both sides smooth | Any one eligible fraction | Any odd subset |
|---:|---:|---:|---:|---:|
| 32 | 3,892 | 4,022 | 4,023 | 4,023 |
| 48 | 2,056 | 2,479 | 2,480 | 2,480 |
| 64 | 488 | 632 | 632 | 632 |
| 80 | 73 | 93 | 93 | 93 |
| 96 | 12 | 17 | 17 | 17 |

The final column is an exact per-cohort ceiling for the specified one-target construction. It is a finite measurement. The earlier single-fraction counting bound does not apply directly to a product of several rows, whose height can exceed N. Enumeration matched direct search at all 6,856 unit targets for 3≤N≤150. The cofactor solver matched exhaustive subset search on 10,006 small lists. The public entry point also matched 1,166 end-to-end brute-force decisions over all unit targets for 2≤N≤43 at cutoffs 1 and 3, including verification of the returned multipliers.

## Which cancellations can move between independent targets?

Let N=pq for distinct odd primes. Draw m independent uniform unit squares modulo N. Let C=16N(1+ln N)/φ(N). For integer L≥1, the probability that a prime greater than L appears in any short fraction for each of two different target draws is at most

    binom(m,2) C²/L.

Indeed, for a fixed prime ℓ there are at most 4(N/ℓ)(1+ln N) signed fractions within the height bound with ℓ dividing either coordinate. For example, with ℓ dividing u, count all |u|=ℓa and v>0 with av≤N/ℓ, then allow both signs and both coordinates. A unit fraction fixes one target. There are φ(N)/4 eligible square targets. This gives C/ℓ for the occurrence probability at one draw, even when the fraction is selected after examining the target. Independence permits squaring the probability for two distinct draws. Sum over pairs and use Σ_{j>L}1/j²≤1/L. No independence between fractions of one target is assumed.

Consequently, suppose a batch using only these fractions produces, with probability δ, a dependency with a nonzero vector of target-use parities, cancelling all non-base prime parities. Choose L above the original base bound and at least 2 binom(m,2) C²/δ. Except on an event of probability at most δ/2, every prime above L in a dependency must cancel *within* its target. Some target occurs oddly and must pass the exact one-target test at cutoff L. A uniformly chosen index passes with probability at least δ/(2m). For m polynomial in log N and inverse-polynomial δ, the enlarged cutoff is still polynomial in log N (one may use N/φ(N)≤2 for distinct odd primes). This connects batch success to a measurable one-target condition at a **different** base cutoff. The bound allows full adaptive choice among each target's short fractions. It assumes independent target draws and it leaves open dependencies using each target an even number of times.

## Dependencies using public targets only

Given row atoms, record a parity vector for the target identities as well as one for the atoms. Write A for the map from row selections to atom parity, and B for the map to target parity. The dimension of B(ker A) is rank([A;B])−rank(A), where brackets denote vertical stacking. Each independent target vector retains fresh private square-root signs in the Dixon reduction. A vector in ker A with Bc=0 instead yields a test from the public targets alone. Because each t occurs an even number of times and ∏uv is a positive square, compare its integer square root to (∏v)∏t^(number of uses/2) modulo N and take a gcd.

An experiment used the ordinary Euclidean fractions for 1,024 targets per modulus. It accepted square rough parts directly and let identical nonsquare rough parts cancel in pairs, treating each whole cofactor as one atom. It performed all elimination separately per modulus. On the four 48-bit moduli there were 901 independent zero-target-vector dependencies for the square targets and 925 for the uniform-unit control. Of the first 40 such basis vectors per modulus, 75/160 and 78/160 respectively gave proper factors, splitting all four moduli in each cohort. No saved target roots entered those gcds. At 64, 80 and 96 bits the same budget produced zero dependencies, even allowing those cofactor matches. This measures one batch policy and budget; it does not establish a lower bound on public-target constructions. [Reproduction and gcd verifier](check_cross_target_dependencies.py)
