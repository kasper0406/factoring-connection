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

## Does a known exponent make the evaluator cheaper?

A known smooth representation can make a collected row usable in a gcd. The [earlier degree-one construction](../../notes/zroute_results.md#32-seedless-and-field-free-designs) already establishes this. Here we measure the effect on the *same targets and fractions*. Take t=2^e mod N with e a uniformly sampled (b+64)-bit string, where b is the bit length of N. A Euclidean row u≡tv gives

    (v · 2^ceil(e/2))² ≡ uv · 2^(e mod 2) mod N.

The stored exponent uses the existing base column for prime 2 and a modular exponentiation to compute the witness. In the control, withhold e from elimination and extraction. The collector still sees precisely the same t, cofactor matches and fractions. A row counts only after cancelling its opaque cofactor or recognizing the remaining cofactor as a square. For the random-square comparison use t=x² and save the sampled x as the witness. These are finite-size benchmarks of relation supply, with seeded pseudorandom draws. The table counts binary dependencies, then verified proper gcds among the first forty basis vectors per modulus. All linear algebra is separate per modulus before summing.

| Input bits | Targets per modulus | Target source | Rows after cofactor matching (four moduli) | Dependencies without witness | Dependencies with witness | Proper gcds / checks with witness |
|---:|---:|---|---:|---:|---:|---:|
| 64 | 1,024 | Powers of two | 667 | 0 | 0 | 0 / 0 |
| 64 | 2,048 | Powers of two | 1,618 | 0 | 0 | 0 / 0 |
| 64 | 4,096 | Powers of two | 4,152 | 0 | 1,920 | 77 / 160 |
| 64 | 4,096 | Uniform unit squares | 4,210 | 0 | 1,978 | 85 / 160 |
| 64 | 8,192 | Powers of two | 11,283 | 0 | 9,031 | 77 / 160 |
| 80 | 8,192 | Powers of two | 813 | 0 | 0 | 0 / 0 |
| 96 | 8,192 | Powers of two | 76 | 0 | 0 | 0 / 0 |

At 64 bits the 4,096-target power and square samples both factored all four moduli with their witnesses. For arbitrary public unit targets at 8,192 per modulus we checked 1,233,037 convergents and obtained 11,233 rows after cofactor matching. All 8,973 dependency dimensions used nonzero target parities. No public gcd candidate survived. These larger pools show exactly where the implementation benefited from knowing the target construction. They provide no improved asymptotic factoring bound. [Power and square comparison](check_known_power_targets.py) · [Public-unit control](check_public_pool_scaling.py)

### Generating near-uniform targets together with labels

There is a public polynomial-time construction for every N≥2 when the labels may use arbitrary unit generators. Put G=(Z/NZ)* and b=bitlength(N). Draw k=b+2λ independent uniform units g_i. For each target draw a fresh uniform binary word a and compute t=∏g_i^a_i modulo N. Store a. This supplies an exact representation over a fixed dictionary; the generators themselves can be large. Rejection sampling the units uses no factorization and takes polynomial expected time since N/φ(N)≤b+1.

Over the random dictionary, the expected total variation distance of the conditional target law from uniform on G is at most 2^(−λ−1). For any two different words, choose a position on which they differ. Conditioning on all the other generators leaves the ratio of their products uniform on G. Thus the expected collision probability for two random words is 2^(−k)+(1−2^(−k))/|G|. Cauchy–Schwarz on the squared deviations of the conditional probabilities from 1/|G| gives the distance bound. By coupling, m fresh targets from a shared dictionary, together with the dictionary, are within m·2^(−λ−1) of that dictionary with m independent uniform units. This statement compares target values; revealing a label alongside its target explicitly certifies their relation.

In the finite control take k=b+64 (λ=32), and compare the same targets, rows and cofactor matches with the binary words retained or withheld. Given t=∏g_i^a_i and u≡tv, a selection with even total word in every generator column and square positive ∏uv yields the public congruence between √(∏uv) and (∏v)∏g_i^(total exponent_i/2). The benchmark checks that congruence before taking each gcd.

| Bits | Targets per modulus | Projected rows, four moduli | Dependencies with words withheld | Dependencies with words retained | Verified proper gcds |
|---:|---:|---:|---:|---:|---:|
| 64 | 4,096 | 4,177 | 0 | 1,428 | 83 of 160, splitting all four moduli |
| 80 | 8,192 | 854 | 0 | 0 | No candidate |

These labels use arbitrary units sampled once per modulus. All the 64-bit splits here were extracted from dependencies using the public words; no factor was encountered in generating either dictionary cohort. For 64 bits the joint distance bound for the 4,096 target *values* is at most 2^(−21). The protocol constructs the target before testing it. It supplies no representation of an externally prescribed target. [Binary-dictionary experiment](check_labelled_dictionary.py)

### What large-cofactor matches can contribute on independent targets

For any fixed polynomial budget m≤b^a, consider *all* fractions 0<|uv|≤N with unit denominator for each target. On independent uniform unit targets, the following two events together have probability smaller than any fixed inverse polynomial in b, asymptotically:

1. Some target has a single fraction whose |uv| has squarefree kernel supported on primes ≤b².
2. Two different draws have fractions with the same nonsquare cofactor after stripping all primes ≤b² from |uv|. Repetitions on two independently sampled but coincident residues are included.

To see this, fix a desired exponent c>0. For a uniform unit target, the prime-occurrence bound in the earlier collision proof applies with C=4N(1+ln N)/φ(N)=O(b²). Choose L=b^K with fixed K>2a+4+c and K>2. The probability of sharing a prime >L between any two draws is O(m²b⁴/L). An identical cofactor having any prime >L is covered by this bound. Otherwise that cofactor is L-smooth. In that case, or for an already-square rough part in the first event, at least one short fraction has L-smooth squarefree kernel. By the [short rational counting bound](https://github.com/kasper0406/factoring-connection/commit/df10c6e) there are at most N^(1−min(1/2,1/K)+o(1)) such unit residues; the integer count there is independent of N's factorization. Divide by φ(N)≥N/(b+1) and take a union bound over m draws. The result is o(b^(−c)). For example, to see the totient bound, the i-th distinct prime dividing N is at least i+1, so ∏p/(p−1)≤∏(i+1)/i≤b+1.

For dictionary targets the same conclusion holds with λ=b by the joint distance bound. Hence even polynomial-time generation of well-distributed labelled targets leaves both direct hits and chance equal-cofactor matches across independent targets too rare for this policy. Equal cofactors produced by different fractions of the **same** target remain possible. Deliberately correlating distinct targets or using a different representation family also falls outside the bound. The proof does not assert uniformity for powers of the fixed base 2.

### A single random generator

There is also a single-row coverage bound for a uniformly random unit base g and an exponent sampled uniformly from a range of M≥N consecutive integers. For a fixed semiprime N=pq write A=p−1 and B=q−1. Relative to local primitive roots, let u,v be the independent uniform discrete logarithms of g (used only in the proof). Set d=gcd(A,u) and f=gcd(B,v). The index of the subgroup generated by g in the unit group is

    d f gcd(A/d,B/f) ≤ d f gcd(A,B).

The expected index at this modulus is at most τ(A)τ(B)gcd(A,B), since averaging gcd(A,u) costs at most τ(A), the divisor count. Over distinct primes p,q∈[X,2X] the expected index is O(log⁸ X). Indeed, sum the displayed upper bound with the divisor weights over all positive A,B≤2X. Expanding gcd(A,B)=Σ_(d|A,B) φ(d) and using τ(dr)≤τ(d)τ(r) bounds this sum by O(X² log² X)Σ_(d≤2X)τ(d)²/d. The last sum is O(log⁴ X) by τ(d)²≤τ₄(d), where τ₄ counts ordered quadruple factorizations. Divide by the order of X²/log² X distinct prime pairs from the prime number theorem.

Each possible value of g^e has probability at most 2/ord(g). Thus for a base of at most (log N)^A₀ primes with fixed A₀>1, a direct short rational representation smooth up to a square has probability at most N^(−min(1/2,1/A₀)+o(1)) on this balanced-prime, random-generator ensemble (N is of order X²). This bounds direct representations only. It allows the factor base to depend on N and g when its size bound is uniform; it must be fixed before e is drawn.
