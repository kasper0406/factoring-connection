# Complete quadratic weighting of the fixed phases

Write $M=N^2$, $e_m(x)=\exp(2\pi i x/m)$, and $C_P(t)=N^{-1}\sum_{x\bmod M}e_M(tP(x))$ for the literal cubic or quartic. At unit scales the CRT stationary formula gives an extension $W(t)=\sum_f w_f e_M(ft)$ defined on all scales. Its values at nonunits can differ from $C_P(t)$. The formulas below specify the domain and normalization of each average.

## A native evaluator for projected literal sums

Take the outer quadratic coefficient $Na$ with $a$ a unit modulo $N$. For the literal raw sum define
$$J=\frac1{N\sqrt N}\sum_{t\bmod M} C_P(t)e_M(Nat^2+Bt).$$
Writing $t=u+Nv$ projects onto $N\mid P(x)+B$. The surviving kernel is the complete quadratic Gauss sum modulo $N$ with linear coefficient $(P(x)+B)/N$. Put $d=(4a)^{-1}\bmod N$ and $\eta=\epsilon_N(a\mid N)$, where $\epsilon_N$ is one for $N\equiv1\pmod4$ and $i$ otherwise. On writing $x=x_0+Nw$ the exact identity is
$$J=\eta\sum_{P(x_0)+B\equiv0\ (N)}\frac1N\sum_{w\bmod N}e_N\big(-d(\beta_0+wP'(x_0))^2\big),\qquad \beta_0=\frac{P(x_0)+B}{N}.$$
This expression is used for the proof and the factor-supplied control.

For squarefree $N=pq$, a term with derivative nonzero at either prime has magnitude at most $1/\sqrt{\min(p,q)}$. There are at most $\deg(P)^2$ terms. Terms with zero derivative at both primes are exactly the stationary extension under the same quadratic weight. Test the gcds of $B+f$ with $N$ for each public local label $f$. A proper gcd yields a divisor. Otherwise the projection retains a sole public diagonal when $B\equiv-f\pmod N$, or gives zero. For the diagonal compute its weight and the one standard Gauss phase via the Jacobi symbol. Trial division through $(2\deg(P)^2/\varepsilon)^2$ either gives a factor or bounds the discarded terms by $\varepsilon/2$. Including the numerical phase, this gives polynomial bit complexity in $\log N$ and $1/\varepsilon$ with no supplied factors. [Classical Gauss evaluations](https://mathweb.ucsd.edu/~revans/contents.html) · [Stationary phase](https://doi.org/10.4064/aa-80-1-1-48)

If all $B+f$ are units, each quotient root is simple at every component. Then the inner change of variable is bijective and the two Gauss factors cancel their phases. The exact value is $J=N^{-1/2}\#\{x_0:P(x_0)+B=0\pmod N\}$, also for repeated prime powers. This identifies the additional root-count information at precision $O(N^{-1/2})$. For instance, if the offset is constructed with a public factorization of the shifted polynomial into distinct linear factors modulo $N$, and all differences are units, the count is exactly $\deg(P)^2$ on two-component inputs. A general critical offset with repeated factors needs a separate calculation, since its discarded terms need not be small.

## What changes for unit outer coefficients

For $A$ a unit and $D=(4A)^{-1}\bmod M$, the stationary extension obeys
$$\frac1N\sum_{t\bmod M} W(t)e_M(At^2+Bt)=e_M(-DB^2)\sum_f w_f e_M(-Df^2-2DBf).$$
If the sum uses only $t$ coprime to $N$, insert the exact mask $\gcd(B+f,N)=1$. The block $t=u+Nv$ survives at the unique stationary $u\equiv-(B+f)/(2A)\pmod N$, which proves the mask. Test the gcds for each local label. A proper gcd splits the input, and otherwise a global collision deletes that label at every component. In particular, the unit-scale raw cubic at $B\equiv\pm2$ and raw quartic at $B\equiv1$ leave a single public line.

For a spectrum satisfying $f^2+u f+v=0\pmod M$, the *all-scale extension* transform is a point evaluation at $(u-2B)D$ times the public phase $e_M((v-B^2)D)$. Raw cubic and quartic have $(u,v)=(0,-4)$ and $(1,0)$ respectively. A generic unit-coefficient weighted evaluator of the extension therefore retains the original magnitude problem.

For the unit-scale squared magnitude, local labels are $0,\pm h$, with $h=4$ for cubic and $h=1$ for quartic. Set $B=-t_0/(2D)$ for a desired unit scale $t_0$. When the mask retains all labels, the corrected transform as a function of $D$ is a sequence of at most $R=2^k$ unit-circle modes for $k$ prime-power components. Indeed $f^2$ has labels $0$ or $h^2$ componentwise. At $D=0$ the sequence gives the squared magnitude. For fixed $k$, choose $2R$ samples $D=jD_0$ at one of $D_0=1,2R+1,2R+2$. The three sets are disjoint. For large $N$, at most two values of $D$ have a global label collision, a proper gcd yields a split, and small $N=O(R^2)$ can be trial-divided in time polynomial in $R$.

Extrapolation is stable without separating the modes. A unit-mode sequence $y_j$ of rank at most $R$ has a backward characteristic recurrence $y_j+\sum_{i=1}^R b_i y_{j+i}=0$ with coefficient sum at most $2^R$, padding the roots if needed. With samples $\tilde y_1,\ldots,\tilde y_{2R}$ at absolute error $\delta$, fit $b$ on the $R$ observed rows using least squares with penalty $\delta^2\|b\|_2^2$. The true recurrence gives an objective at most $(R+1)4^R\delta^2$. For the fitted $b$, the true residual sequence is itself of rank at most $R$. Its first $R$ residuals have size $O(R2^R\delta)$, so its unobserved residual at zero has size $O(R4^R\delta)$ by the true recurrence. The estimate $-\sum_i b_i\tilde y_i$ has the same error bound. Exact Gaussian-rational linear algebra charges finite precision and avoids unstable root fitting. Eight samples suffice for two components. A fixed inverse-polynomial point accuracy requires inverse-polynomial chirp accuracy in that case.

## Literal all-scale quartic square

The actual $|C_P(t)|^2$ includes nonunit contributions. For the quartic and a unit $A$, at offset $B=0$ its exact transform is
$$T(A)=\frac1N\sum_{t\bmod M}|C_P(t)|^2 e_M(At^2)=\prod_{q=p^e\mid\mid N}\big(C_q+4e_M(-D e_q)\big),\qquad C_{p^e}=p^{e-1}\big((5e-2)p-5e-\chi_p(-1)\big),$$
where $e_q$ is the CRT idempotent for the $q^2$ component. Here $C_q$ counts pairs with equal quartic value modulo $q$.

Locally at $p^{2e}$, the sector $p^e\mid t$ gives $C_q$ by summing over the lower-modulus scales and using additive orthogonality. At valuation $h<e$ the complete quartic reduces to modulus $p^{2e-h}$. Stationary phase gives squared magnitude $p^h$ times a constant five plus two opposite harmonics with weights $2$ for even $h$ and $2\chi_p(-2)$ for odd $h$. The constants average to zero with the restricted unit quadratic weight. For $h>0$ each harmonic has unit linear coefficient and a quadratic coefficient divisible by $p$, so cancels on the upper digit. At $h=0$ the two harmonics give the equal Gauss value $4e_M(-D e_q)$. This proves the transform over repeated factors as well.

For the count, the equality factors as $(x-y)(x+y)(x^2+y^2-2)=0$. Modulo $p$ the two lines and conic give $3p-5-\chi_p(-1)$ points, with five transverse crossings and no triple crossing. Each smooth point lifts to $p^{e-1}$ points modulo $p^e$. At a crossing the two functions are local coordinates. Count $p^e\mid uv$ with $p\mid u,v$. The zero $u$ contributes $p^{e-1}$ choices of $v$. Each valuation $1\le v_p(u)<e$ contributes $(p-1)p^{e-1}$ pairs. Combining the five crossing counts with the smooth points yields the displayed $C_q$.

With two components $P,Q$, five values at $D=jt_0$, $1\le j\le5$, recover the original quartic point magnitude. Put $x=e_M(-t_0 e_P)$, $y=e_M(-t_0 e_Q)$, and $z=xy=e_M(-t_0)$, which is public. Subtract $16z^j$ and take consecutive differences, scaled by $\sqrt N$. The four differences $h_j$ satisfy $h_{j+2}-(x+y)h_{j+1}+zh_j=0$. From samples at absolute error $\delta\sqrt N$, estimate $x+y$ by the ratio using the larger measured denominator $h_2,h_3$ when it exceeds $\delta^{2/3}$. The slope error is $O(\delta^{1/3})$. If both are below threshold, the recurrence bounds $|h_0|$ by $5|h_2|+2|h_3|$. Since $-\Re h_0=4[C_Q(1-\Re x)+C_P(1-\Re y)]/\sqrt N$ and $\max(C_P,C_Q)\ge\sqrt N$, one of $x,y$ is within $O(\delta^{1/3})$ of one. The other is then close to $z$. In both branches the magnitude estimate has error $O(\delta^{1/3})$, including imbalanced and repeated components. The large unknown constant $C_P C_Q$ cancels in the differences.

At the same zero offset the *unit-only* squared quartic transform is public given the component count: the mask retains only labels $\pm1$ at every component, so the value is $4^k e_M(-D)$. These observations leave the native unit-coefficient literal evaluator and the original point magnitude open. They separate projections where the useful phase has been removed from weighted averages that retain it.

The diagnostic matched 24 literal projected chirps, nine literal squared chirps, 12 unfactored coarse projections, six noisy eight-sample extrapolations, and nine noisy five-chirp point evaluations. [Evaluator and checks](quadratic_chirp.py)

```sh
python experiments/polylog_factoring_2026_11_27/quadratic_chirp.py
```

## Regular quartic offsets on all odd moduli

At gcd(N,B(B−1))=1 all roots of x⁴−2x²+B are simple modulo each prime
component. Thus the exact normalized projection is the root count divided
by √N, as above. This yields a uniform native additive evaluator even
when neither the number of components nor their exponents is bounded.
Strip primes through 64; count roots modulo each stripped prime, using
unique simple lifts for every exponent. Write c for the product of these
counts, m for the remaining cofactor and S=N/m. The number of roots on
m is at most 4^ω(m)≤m^(1/3), so |J|≤c/(√S m^(1/6)).
If this is at most 1/(2R), return zero. Otherwise m<(2Rc)^6:
enumerate the roots on m and return the exact count divided by √N at
sufficient numerical precision. Since c≤4^π(64), this costs polynomial
bit time in log N and R, with an additional polynomial input-reading
cost for B. The elementary enumeration has a large worst-case constant.

## Exactly critical quartic offset B=1

Write d=(4a)^(-1) mod N and η=ε_N(a|N). In the projected identity
sum over all x mod N² and put z=x²−1. The projection condition is
N|z² and the phase is e_{N³}(−d z⁴), with normalizing factor η/N.
Locally on p^e exactly dividing N the condition forces p|z. Each such
z mod p^(2e) has two preimages, by simple lifting of x=±1. The local
character has unit coefficient c_p=−d((N/p^e)^3)^(-1) mod p^(3e).
Thus the local factor, with the multiplicity two separated, is
p^(−e) times the sum over z mod p^(2e) of valuation ≥ceil(e/2).

For v_p(z)=h with 4h<3e, the sum reduces to a character on the units
of modulus p^(3e−4h). There are enough digits since h≥ceil(e/2).
For modulus exponent at least two the last-digit sum vanishes, as the
derivative 4c_p u³ is a unit. The region 4h≥3e contributes p^floor(e/4).
The only remaining layer satisfies 4h=3e−1, possible precisely at
e=3 mod 4. Including the zero-residue term from the deep region, the
local factor then becomes p^floor(e/4) ∑_{u mod p} e_p(c_p u⁴).

In particular if no e is 3 mod 4 the exact projection is
J=η 2^ω(N) ∏p^floor(e/4). On a cubefree input of known component count
this critical offset is therefore easy. It also explains why using the
squarefree-semiprime diagonal prediction unchanged on arbitrary powers
would be incorrect.

For odd n=∏p^e take a=B=1 on moduli n² and n⁴. Both Gauss prefactors
are one and neither query has an exceptional local exponent. The values
are integers: J(n⁴)=2^ω(n)n and J(n²)=2^ω(n)∏p^floor(e/2).
Hence two absolute-error-<1/2 evaluations determine ω(n) and the
largest s such that s² divides n exactly. This can split mixed-multiplicity
inputs; on squarefree semiprimes it gives ω=2 and s=1. The pair of
queries is not a complete factoring reduction.
[Regular evaluator and critical power decoder](quartic_projected_powers.py)

## A fixed projected cubic value on fourth powers is factoring-equivalent

Use the branch's literal raw cubic P(x)=3x−x³ with a=1 and B=2. Write
J(M) for its normalized projection at odd modulus M. For a base v
coprime to 6 and prime-power components q=p^e exactly dividing v,
$$
             J(v^4)=\prod_{q\parallel v\ (q=p^e)}(q+q^{-2}).
$$
The product is over maximal prime powers only. In particular this is real.
The evaluation promise below is a numerical estimate on every such input;
an optional return of an arbitrary proper split of the *padded* modulus
is not the same oracle.

**Local calculation.** The shifted polynomial is Q(x)=(2−x)(x+1)^2.
Use the projected sum η/M times the sum over x mod M² with M|Q(x)
of e_{M³}(−dQ(x)²), where d=4^(-1) mod M. Here η=1 since M is a fourth
power. At p^h exactly dividing M, h=4e. CRT gives a unit local character
coefficient, regardless of the other prime powers. The simple root class
x=2 mod p maps by Q bijectively onto multiples of p modulo p^(2h),
by simple lifting. Restriction to p^h|Q leaves a complete quadratic
Gauss sum divided by p^h, equal to p^(−h/2) since h is even.
In the other class put y=x+1 and z=y sqrt(1−y/3), taking the square
root equal to one mod p. This is bijective on multiples of p: the square
root lifts uniquely, and the map lifts one digit at a time with derivative
one mod p. Here Q=3z². The quartic-character layers with 4v_p(z)<3h
cancel as in the quartic calculation above, since 3h is divisible by four
and the remaining character modulus cannot be p¹. The deep region gives
p^(h/4) after normalization. Adding the two contributions proves the formula.

**Extracting exact arithmetic from unit additive error.** Remove factors
2,3,5 from the target n. Choose u a power of 5 exceeding 4n², and
query the single public base nu at moduli (nu)^(4t). Multiplying the
answer by n^(2t)/(u^t+u^(−2t)) and rounding recovers exactly
H_t=∏(1+q^(3t)) for the unknown components q of n. The scaled error
is less than 1/4 even when the original error is one. For n=pq a distinct
semiprime, **one query suffices**: subtract n³+1 from H₁ to get p³+q³.
Solve s³−3ns=p³+q³ by integer bisection on [ceil(2√n),n+1], where
the left side is strictly increasing, then recover p,q from their sum.

**Arbitrary component count.** Put b=bitlength(n) and
L_t=log H_t−3t log n. If P_j=∑q^(−3j), then
L_t=∑_{h≥1} (−1)^(h+1)P_{th}/h by the log series. Its Dirichlet inverse
has coefficients c_m=(1/m)∑_{2^v|m}2^v μ(m/2^v). Indeed the numerator
sequence (−1)^(m+1) is the convolution of the constant-one sequence
with the sequence supported on 1,2 having values 1,−2. Inverting those
two sequences gives the formula. Thus P_j=∑_{m≥1} c_m L_{jm}.
Absolute convergence follows from q≥7 and |c_m|<2.

Compute j=1,...,b, truncating m at b+5. The omitted absolute error
is bounded by 2b ∑_{m>b+5} 2^(−3jm), smaller than 1/(8n^(3j)).
Rounding n^(3j)P_j recovers the exact integer power sum of the unknown
integers (n/q)^3. The implementation computes the two integer logarithms
with correctly rounded Decimal ln, then uses rational arithmetic for their
difference and weighted sum. The precision of 3bj+3b+30 decimal digits
leaves more than the required margin: the integer logs are O(bj(b+5))
in magnitude, the subtraction multiplies the log n error by at most
3j(b+5), and at most b+5 terms of coefficient magnitude <2 are added.
Newton identities recover the monic polynomial with these roots padded
by zeros to degree b. Deterministic rational polynomial factorization
and integer root extraction recover the prime powers and their prime bases.
All queried exponents t≤b(b+5), the moduli have O(b³) bits, and all
arithmetic and oracle counts are polynomial in b. Conversely the factors
of the queried base give the rational projection by the local formula.

The numerical decoder was checked with signed unit errors and factor-generated
controls including squarefree, mixed-multiplicity and prime-power inputs.
The local identity was separately checked by literal Fourier evaluation
at 5⁴ and by direct projected lifting sums at 5⁴, 7⁴ and 11⁴. The remaining
algorithmic target is an unfactored unit-additive evaluator at this one
critical offset on fourth powers.
[Decoder](cubic_critical_factoring.py)

## Unpadded projected cubic on squares: a local character mask

Let n be coprime to 6, M=n², B=2+Mv with v a unit mod n. Write t=−dv² mod M
with d=(4a)^(-1) mod M, and let E_q be the CRT idempotent equal to one
mod q² and zero at all other components, where q=p^e exactly divides n.
The exact literal cubic projection is
$$
 J(M,a,2+Mv)=\prod_q\left[e_M(t E_q)+1+
                    \left(\frac{-3v}{p}\right)+\frac1q\right].
$$
This uses the projected raw cubic, without an outer magnitude or extension
at nonunit scales. The sign mask and unit scale can be set independently:
for desired unit u take v=−3u and a=−v²/(4t) mod M. Then (−3v|p)=(u|p).

For the calculation let m=n/q. The CRT coefficient at p in the projected
phase is −d m^(−6) mod q⁶. On the simple class x=2 mod p the shifted
polynomial still maps bijectively onto multiples of p, contributing 1/q
by the quadratic Gauss sum at the even exponent 2e. At x=−1 mod p use
the same local coordinate z as above, so Q=3z². The constraint q²|Q+Mv
requires q|z. Put z=q m w locally. Reducing to the true period leaves
q^(−1) times the complete sum modulo q² of the phase with polynomial
−d m^(−2)(3w²+v)². Summing the top e digits leaves evaluations at
stationary residues modulo q. They are w=0 and the roots of 3w²=−v;
these classes are disjoint and simple as roots of the derivative since
p≥5 and v is a unit. The first phase is e_M(t E_q) and the others are
one. The outer Gauss prefactor is one on M=n². This proves the formula.

For exactly two distinct prime-power components, choose a random unit u.
With probability 1/2 its two Legendre symbols differ. Subtract the masked
projection from the all-positive projection at the same t. On a mixed mask
this difference is exactly 2(x+2+1/q) with x=e_M(t E_q) at the surviving
positive component. Thus (difference)/2−2 estimates x within 1/q+1/32
when each complex answer has absolute error ≤1/32. Since q≥5 this is
below sin(π/12). The closest of the four cardinal directions gives the
true phase in turns within 1/8+1/24=1/6. Query t=2^j mod M through
j=bitlength(M). From the last quarter-turn estimate work backward,
choosing the half-angle closest to the preceding observation. The correct
branch is within 1/6+1/12=1/4 turn of it, strictly, so the other half-angle
loses. The final error is <1/(6·2^j); multiply by M and round to recover
E_q. Its gcd with n is a verified proper divisor. Nonunits among the
sampled masks split immediately. Failed masks are harmless because the
decoder verifies every output. Constant expected trials and O(log n)
calls per trial suffice; all calls in a trial have the same unpadded M.
On a squarefree semiprime one may preselect Jacobi-negative masks.

For quadratic residuosity on a distinct semiprime, a Jacobi-negative unit
is already a nonsquare. For a Jacobi-positive unit the two local signs
agree. In the positive case the masked and unmasked values coincide.
In the negative case their difference divided by two is
2+1/p+1/q+x+y at scale one, with x,y replaced by x²,y² at scale two.
Since cos θ+cos 2θ≥−9/8, at least one real difference is ≥7/4. Four
constant-accuracy queries distinguish the cases. [Decoders](cubic_shifted_projection.py)

Direct quotient-lift sums checked the identity on prime, composite and
repeated bases. Noisy factor-supplied controls separately exercised
phase rounding and character masking. The native complex evaluator on
these shifted critical inputs remains open.

## Unpadded splitting with arbitrarily many components

The same projected cubic on M=n², B=2+Mv now suffices for randomized
expected polynomial-time factoring on arbitrary inputs, with absolute
complex error at most 1/256. Per trial all oracle calls use that one
unpadded modulus. Strip twos, test primality and handle perfect powers
with ordinary integer roots. On each remaining composite first trial-divide
through 2b for b=bitlength(n). A nonunit encountered during random sampling
immediately gives a proper split. Below, condition on uniform units.

Let k be the component count, bounded by the public K=floor(b/floor(log₂(2b))).
At scales t=1,2 write x_i=e_M(tE_i), c_i=1/q_i. For a mask u the local
factor is x_i+c_i+1+χ_i(u), where χ_i is the prime's Legendre symbol.
A positive sign has magnitude ≥1+c_i, a negative sign ≥1−c_i. Thus the
whole magnitude exceeds 1−∑c_i≥1/2 after trial division. Its log magnitude
is exactly a constant plus ∑β_i(t)χ_i(u), with
β_i(t)=(1/2)log(|x_i+2+c_i|/|x_i+c_i|). All these coefficients are
nonnegative. For each i at least one of the two scales has Re x_i≥−1/2:
if cos θ<−1/2 then cos 2θ>−1/2. There the squared magnitudes differ by
4(1+c_i+Re x_i)≥2 and |x_i+c_i|≤1.2, so β_i>0.21.
The specified absolute error perturbs log magnitudes by less than 0.01.

Choose s independent random unit generators with 8(K+1)²≤2^s<16(K+1)².
By CRT their Legendre signs give each component an independent uniform
s-bit signature. With probability above 7/8 all signatures are nonzero and
distinct, by a union bound. Evaluate both scales on the 2^s subset products
of the generators. In the normalized Walsh transform, the coefficient at
each nonzero signature is exactly its β_i, and unoccupied coefficients
are zero. Each coefficient's error is at most the maximum input error.
Select signatures with coefficient >1/8 at either scale and keep the
stronger scale per signature. After multiplying all subset products by a
fresh unit w, its selected coefficients become β_i χ_i(w), so their signs
decode the full vector of local characters. Neither prime values nor
discrete logarithms are used. With k+4 fresh vectors the chance of not
spanning the binary k-space is <1/16 by the union bound over nonzero
annihilating linear forms. Binary row elimination while multiplying the
associated units constructs a single mask negative at just one component.

Finally let A be the all-positive complex value and C the value with this
single sign flipped. At any unit t, A=H(x_i+c_i+2) and C=H(x_i+c_i),
where |H|≥1. Then 2C/(A−C)=x_i+c_i. Using perturbed values changes this
ratio by at most (2+2·1.2)/(256(2−2/256)); including c_i≤1/5 keeps the
complex phase error below sin(π/12). The existing quarter-turn decoder at
dyadic scales recovers E_i, and its gcd with n is verified. A failure of
signature separation or span costs at most one bounded trial; an accidental
verified split is equally usable. Recursing on proper divisors terminates
in O(b) nodes. At most O(K³+b) evaluations are needed per trial, with
polynomial bit work including the transforms. Input log magnitudes are
O(b), and the strict numerical margins permit polynomial precision.

No magnitude-only oracle is assumed for the final extraction. The numerical
checks used independent noisy factor-supplied values for three to five
components, including repeated exponents, and sampled only units to exercise
the mask decoder rather than incidental gcds.
[Multi-component decoder](cubic_mask_spectrum.py)

### Relative magnitudes suffice, even with unknown gain at each scale

For the unpadded shifted cubic, the oracle may return only positive
magnitudes, within relative error 1/4096. It may multiply all answers at
a fixed scale t by any common positive gain, provided its outputs have
polynomial descriptions. This weaker oracle still suffices for randomized
expected polynomial-time factoring. Trial-divide through max(4096,2b)
first. Each log-magnitude error is <0.001. A scale-dependent additive
constant in the logs changes only the zero Walsh coefficient, which the
mask decoder discards. Thus it recovers the same singleton flip with the
earlier constant success probability.

At a chosen component write c=1/q<1/4096, and x=e_M(tE_q) with real
part z. For the all-positive and flipped magnitudes the squared ratio is
R²=[5+4z+2c(2+z)+c²]/[1+2cz+c²]. The common gain cancels. R² differs
from 5+4z by at most (24c+10c²)/(1−c)², and the measured squared ratio
has relative error below 0.001 since each magnitude has relative error
≤1/4096. Consequently (measured R²−5)/4 estimates z within 0.01.
Clipping to [−1,1] preserves the error bound.

To resolve signs, first collect cosines at all dyadic scales used in the
quarter-turn decoder. One dyadic reference r has measured absolute cosine
≤0.9. Indeed if the true sine stays below 1/2 in absolute value, the
distance of each phase in turns to the nearest half-integer is <1/12.
It doubles on the next step as long as this continues. Its initial
nonzero distance is at least 1/(2M), so this cannot persist through
bitlength(M) doublings. At a true absolute sine ≥1/2 the observed absolute
cosine is ≤√3/2+0.01<0.9.

Choose the orientation in which the reference sine is positive. For each
dyadic t use an additional cosine at t+r, or t−r if t+r vanishes mod n.
A proper gcd of a proposed scale with n already splits; both choices
cannot vanish mod n. With sign s=±1, derive the oriented sine as
s(cos(t)cos(r)−cos(t+sr))/sin(r), writing angles in scale notation.
The estimated denominator is at least √0.19 and the true one exceeds
0.41; their difference is at most 2.4·0.01 by rationalizing the square
roots. The numerator error is at most 3·0.01, so the sine error is below
0.13. Together with the cosine this gives complex error below sin(π/12).
Quarter-turn rounding recovers E_q or −E_q, with identical gcd. The
extra real queries cost O(b) on top of the mask learning. Factoring
conversely evaluates the unscaled magnitude by the local product formula.

Checks used only real outputs on three to four large components, with
signed endpoint relative errors and varying gains. A separate noisy
cosine check covered small moduli and repeated factors. These checks use
factor-supplied oracle controls; they do not implement a native estimator.
[Magnitude mode](cubic_mask_spectrum.py) · [Cosine decoder](cubic_shifted_projection.py)

### Recovering normalized magnitudes from a comparison bit

Let b be the bit length of n, with gcd(n,6)=1. At unit scale t and
unit mask u, write R=|J_n(t,1)|/|J_n(t,u)|. Each flipped prime power
q contributes |x+2+1/q|/|x+1/q| for a unit complex number x. This
lies between 1 and 4 by squaring for the lower bound and the triangle
inequality for the upper bound (q≥5). Thus log₂ R lies in [0,2b].

Append L=2b+14 known prime powers q coprime to each other and n, each
at least 10⁸ L. With s their product, use the padded base ns and scale
s²t + Σ k_q(ns/q)². In the local product this preserves every original
phase and sets each new phase to exp(2πi k_q/q²). For a desired gain g
choose a phase with cosine (g²−5)/4, so that |x+2|=g. Round its frequency
to a unit modulo q² by moving at most 1.5 lattice steps. Adding 1/q to
the numerator and denominator perturbs their magnitudes by at most 1/q
each, and rounding changes x by at most 3π/q². Since both ideal
magnitudes are at least one, the log₂ gain error is below 12/q.

Use 2b gains of 2 and one gain of 2^(2^(-i)) for each i from 1 to 14.
CRT masks activate any selected gains against an all-negative padding
baseline. Comparing the original mask u with baseline mask 1 on these
two padded inputs tests whether the selected log₂ gain exceeds log₂ R.
Binary search on a grid of width 2^(-14) takes O(log(b+1)+14) queries.
Even if replies are arbitrary within log₂ distance 2^(-16) of a tie,
each retained interval expanded at both ends by 2^(-16)+12·10^(-8)
contains the true value. Its final midpoint therefore gives 1/R within
relative error 1/4096. This normalizes all masks by the same scale gain,
as required by the magnitude reduction above.

The padding has polynomial bit length. Only O(b) small primes can be
excluded by n, and the prime powers can be built by repeated multiplication.
The comparison oracle must work on the padded inputs. Constructing such
a comparator from the unfactored input remains open.
[Comparison calibration implementation](cubic_magnitude_comparison.py)

### Two fixed scales decide quadratic residuosity

For a unit u modulo n coprime to 30, put R(t)=|J_n(t,1)|/|J_n(t,u)|.
For squares, R(1)=R(2)=1. Otherwise pick any prime where the local
Legendre sign is negative. Write a=1/q for its reciprocal prime power.
Its squared gain is 1+4(1+a+c)/(1+a²+2ac), with c the cosine of its
phase. This increases with c since 1−2a−a²>0, and is at least 3 at
c=−1/2. At least one of c and 2c²−1 is ≥−1/2, so scale 1 or 2
gives R(t)≥√3. All other local gains are at least one.

A single extra factor 5 sets the comparison threshold. Give it local
phase exp(2πi·12/25), negative baseline mask 2 and positive candidate
mask 1. Its gain G is between 1.50 and 1.53 by the displayed formula.
The candidate/baseline magnitude ratio is G/R(t). Both comparisons
are positive for squares, and at least one is negative for nonsquares,
even with arbitrary replies for absolute log₂ ratios below 0.17
(since log₂(√3/1.53)>0.17). This directly embeds quadratic residuosity
on Jacobi-positive semiprime inputs into two constant-gap comparisons.
It supplies an explicit test of any proposed native comparison algorithm.

### What the two-scale test counts

For a promised product of two distinct odd primes coprime to 30, the
Jacobi symbol and the two-scale square test together return the exact
number of negative local signs. Jacobi −1 means one. With Jacobi +1
the square test distinguishes zero from two. Repeated exponents change
the Jacobi parity rule, so this shortcut requires the stated promise.

For unrestricted n there is a stronger pairwise query. Let A(u) be the
set of prime divisors where u is a nonsquare, and form
D(t)=R_u(t)R_v(t)/R_{uv}(t). Each prime in A(u)∩A(v) contributes its
squared local gain, and all other primes cancel. Therefore disjoint sets
give D(1)=D(2)=1, whereas an intersection gives max(D(1),D(2))≥3.
The calibrated magnitude oracle above supplies all three normalized
ratios at each scale. Relative errors of 1/4096 leave ample room for
the threshold 2, including when the component sets are large.

This gives a randomized count for arbitrary unit u without learning
individual signs. For uniform independent unit v, disjointness occurs
with probability 2^(−|A(u)|) by CRT. Direct frequency estimation thus
costs exponentially many samples in the support size. For supports
of size O(log log n) it costs polynomially many in log n. The two
fixed scales themselves do not yield a general deterministic count by
this sampling argument.

### Counting dense supports by averaging over unit scales

For odd n of bit length b, trial divide through 16b and count negative
Legendre signs on the removed primes directly. On the cofactor sample
t uniformly from units modulo its square. CRT makes the phases for distinct
prime powers independent. At q=p^e the phase is uniform on primitive
q²-th roots. Put a=1/q, A=2+a, d=q² and h=d/p. The mean natural log
gain for a flipped sign is exactly

    log A + [log(1+A^(-d)) - log(1+A^(-h))
             - log(1+a^d) + log(1+a^h)] / (d-h).

Indeed the product over primitive roots is the product over all d-th
roots divided by the product over h-th roots, both orders odd. Thus the
absolute correction to log A is at most 4·2^(-h)/(d-h). Since
log A−log 2 ≤ 1/(2q), after conversion to base two and summing over
at most b surviving components, the mean log₂ ratio differs from the
unknown support size by less than 1/16. Each local log₂ gain lies in
[0,2]. For m independent scales and k>0 negative components, Hoeffding
bounds the probability that the sample average differs from its mean by
1/4 or more by 2 exp(−m/(32k)). For k=0 the true ratio is always one.
A relative magnitude error of 1/4096 adds less than 0.001 to the
log₂ average even if all errors have the same sign. Rounding recovers
the integer count. Taking m=32b(s+1) gives failure below 2^(-s).
The acceptance probability for sampling units is at least 1−b/(16b)
by a union bound after trial division. The method still needs the
normalized magnitude oracle, or its calibrated comparison simulation.

### A factor-free product for the semiprime count

Suppose N=pq for distinct odd primes p≤q and u is a unit. Compute its
Jacobi symbol by reciprocity. The negative-symbol branch has count one.
On the positive branch put L=floor(sqrt(N)/2) and form

    D = gcd(N, product_{j=0}^L (j²-u)).

The count is two if D=1 and zero otherwise. In the square case a pair
of opposite roots modulo p has a representative between 0 and (p−1)/2.
In the nonsquare case neither prime can divide any term. D=N is a
valid positive witness even though it does not split the semiprime.

A fast product tree gives a deterministic baseline
(more precisely, a fourth-root bound in N). Choose m=ceil(sqrt(L+1))
and build the monic polynomial F(X)=product_{j<m} ((X+j)²-u)
over Z/NZ. Evaluate it at 0,m,2m,... for the complete blocks by a
monic multipoint remainder tree. Multiply the values and finish the
last partial block directly. Monic division works over this composite
ring. Fast polynomial arithmetic gives soft-O(N^(1/4) poly(log N))
bit operations. This is exponential in bit length and gives no
asymptotic improvement over classical deterministic factoring bounds.

Merely shortening the interval to poly(log N) loses the guarantee.
For example take u=4^(-1) mod N. For each prime r dividing N its
least nonnegative root is (r−1)/2. Even on uniform random squares,
any fixed K candidate residues hit a root modulo either prime with
probability at most 2K/(p−1)+2K/(q−1), by a union bound. This bound
concerns candidate lists chosen independently of u. The example does
not rule out adaptive algebraic use of u or another product identity.

### Rational witnesses and the fourth-root scale

For a semiprime N=pq, p<q, search for a congruence a²=u b² modulo
either factor with 1≤a,b≤H. A nonunit among 1,...,H gives a factor
directly and must be handled before interpreting a product collision.
Otherwise any collision certifies a local square. On Jacobi-positive
inputs this certifies a square at both primes.

Every square modulo p has such a witness for H=floor(sqrt(p)).
For a root r, the H+1 residues 0,r,...,Hr are distinct modulo p.
Some adjacent cyclic gap a is at most floor(p/(H+1))≤H. The
corresponding indices differ by a nonzero signed integer of magnitude
b≤H, giving a=±br modulo p. Therefore H=floor(N^(1/4)) suffices
without knowing p. With F(X)=product_{a=1}^H(X−a²), evaluate F
at u,4u,...,H²u modulo N and take a gcd of the product with N.
Product and monic remainder trees use soft-O(H poly(log N)) bits of
work. The earlier difficult example u=1/4 now has witness (1,2).

This change reduces witness height but leaves the fourth-root work
bound. Indeed for a uniform square unit u and H<p the chance of
any collision is at most 2H²/(p−1)+2H²/(q−1). Each pair determines
one square value locally; the bound follows by a union over pairs and
primes. For balanced semiprimes and H=o(N^(1/4)) this probability
vanishes. This limits exhaustive rational witness searches by height;
it leaves open algorithms exploiting the structured product without
locating a small witness.

The local height guarantee is exactly sharp for every odd prime. Put
h=floor(sqrt(p)) and take u=h². For 1≤a,b<h, both hb−a and hb+a are
strictly between 0 and h²≤p. Neither p nor q can divide either factor
of h²b²−a², so the same example also gives a semiprime lower bound.
The witness (h,1) attains height h.

### Exact algebraic degree of the squared magnitude

For a unit mask and unit scale on an odd base n coprime to 6, let
Y=|J_n(t,u)|². Its degree over the rationals is exactly

    D = ∏_{q=p^e || n} (φ(q²)/2) = n φ(n) / 2^ω(n).

Indeed the local factor is f_q=1+c²+c(x+x⁻¹), with rational
c∈{1/q,2+1/q} and x a primitive q²-th root. Its stabilizer in the
local cyclotomic Galois group is precisely {1,−1}. To check that the
product introduces no extra stabilizer, suppose σ fixes Y. Then
σ(f_q)/f_q lies both in the q² cyclotomic field and in the compositum
of the remaining cyclotomic fields. Their conductors are coprime, so
their intersection is rational (equivalently use the CRT product of
Galois groups). This rational ratio has finite multiplicative order
by iteration of σ and is positive since every conjugate is positive.
It is therefore one. Orbit-stabilizer gives the formula.

This also applies to the squared ratio for two masks, with the product
restricted to primes where their signs differ. Each nonconstant local
factor is an invertible fractional linear function of x+x⁻¹, or its
reciprocal: the two constants c and c+2 satisfy c(c+2)≠1 for c=1/q,
so the numerator and denominator are not proportional. Positivity and
the field-intersection argument carry over.

Even the degree alone would suffice for factoring. For a semiprime
with distinct factors, φ(n)=4D/n and their sum is n+1−φ(n).
More generally if b is the bit length, 2^b D is a known multiple of
φ(n) with O(b) bits. The standard randomized reduction from a multiple
of the group exponent to factoring applies (strip powers of two,
try random bases and take gcds at nontrivial square roots of one;
detect prime powers separately). This use of the degree does not
require printing the polynomial. Printing its dense coefficient list
already costs exponentially many bits on prime inputs. High degree
alone gives no comparable bound for arithmetic circuits or approximation.

### One normalized trace suffices on distinct-prime semiprimes

Let T(n) be the mean of |J_n(t,1)|² over unit scales t modulo n².
This is the normalized rational Galois trace. CRT makes the primitive
local phases independent. Each phase has mean zero because its order
is a nonsquarefree prime power (subtract the sum over nonunits from
the full geometric sum). Hence

    T(n) = ∏_{q=p^e || n} [1+(2+1/q)²].

This number has a denominator dividing n² and an O(log n)-bit
numerator, since every local factor is less than 6 for q≥5.
For n=pr with distinct primes and s=p+r, expansion gives

    n² T(n) = 25n² + 20ns + 5s² + 6n + 4s + 1.

Consecutive candidate integer sums give trace values separated by
(20n+10s+9)/n² > 20/n. Thus a rational estimate with absolute
error ≤10/n suffices. Binary search the increasing quadratic on
[0,n+1], select the nearest value, and verify the roots of
X²−sX+n by an integer square root. This needs only O(log n)
precision bits and polynomial decoding work.

A direct average of independent unit-scale samples has a different
cost. For each local factor Y=1+c²+c(x+x⁻¹), Var(Y)=2c² since
the means of x and x² both vanish. Independence gives
Var(Y₁Y₂) ≥ Var(Y₁) E[Y₂]² ≥ 200 for the all-positive semiprime
case. This variance is also O(1). Thus achieving root-mean-square
error 10/n with the plain sample mean takes Θ(n²) samples. This
calculation only constrains that estimator. A symbolic uniform-scale
summation that computes this trace to the stated error from n alone
would factor distinct-prime semiprimes.

### Masked normalized traces give an all-input splitting reduction

Let T(n,u) be the same unit-scale mean square with a unit mask u.
As above n²T(n,u) is integral. Its local numerator at q=p^e is
A=5q²+4q+1 for a positive sign and B=q²+1 for a negative sign.
The gain g(q)=A/B=5+4(q−1)/(q²+1) is strictly decreasing on
[5,∞) and lies strictly between 5 and 6.

Choose k independent uniform units with H=2^k the first power of two
at least 2 bitlength(n). Give each p its k-bit signature of negative
signs. CRT makes these signatures uniform and independent. For each
subset x query the numerator at the corresponding product of masks.
Take its multiplicative Walsh transform, replacing each additive
butterfly by (a·b,a/b). At each nonzero index h the result equals

    ( ∏_{q with signature h} g(q) )^(H/2).

Indeed the logarithm of each local numerator is a constant plus
one half of log(g(q)) times its signature character. A singleton
bucket yields g(q)^(H/2), which can be inverted by exact rational
binary search. Empty and multiple buckets fall outside (5^(H/2),
6^(H/2)). For a fixed component, the probability its signature is
zero or shared is at most ω(n)/H≤1/2. A successful trial therefore
splits n with probability at least one half when it has two or more
distinct prime divisors. A nonunit sampled while choosing masks
already supplies a divisor. Recurse on verified splits, stripping
2 and 3 and handling primes and perfect powers by standard tests.

There are O(log n) trace queries per trial. Transformed rational
numerators and denominators have O((log n)²) bits. Absolute trace
error strictly less than 1/(2n²) suffices by rounding to the known
denominator. This is a randomized polynomial-time equivalence between
factoring and computing a rational scale average for arbitrary masks.
[Exact trace splitter](cubic_masked_trace.py)

### A fixed congruence system whose count is the masked trace numerator

Over Z/nZ take variables t, y₁, x₁, y₂, x₂. Impose t²=t and, for
i=1,2, impose the three equations

    yᵢ(yᵢ²−u)=0,      (1−t)yᵢ=0,      t(yᵢ²−u)xᵢ=0.

For a unit mask u and a prime power q=p^e exactly dividing n, t=0
forces both y's to vanish and gives q² solutions. At t=1, the
factors y and y²−u are coprime at p. Hence each y either vanishes
(with x=0) or satisfies y²=u (with q free choices of x).
Simple-root lifting gives r=1+(u|p) square roots. The two independent
pairs give (1+rq)² solutions on this branch. Thus the local total is
q²+(1+rq)², exactly the local trace numerator for either sign.
CRT proves the global count equals n²T(n,u), including repeated powers.

This reformulation uses five variables, seven equations and total degree
at most four, independently of the bit length. One can eliminate both
x variables by gcds: sum over idempotents t the square of the sum of
gcd(t(y²−u),n) over admissible y. The supplied reference implementation
scans residues, so it does not have a polynomial bit-time bound.
The all-input trace reduction above makes exact counting for this one
system randomized-polynomial-time equivalent to factoring. Showing a
generic fixed-dimension lattice-point algorithm applies would require
handling these nonlinear modular constraints and the binary modulus.

### A fixed cubic with no mask suffices for factoring

A simpler count bypasses the masked trace decoder altogether. Let C(m)
count the roots of x(x−1)²=0 modulo m². At each q=p^e exactly dividing
m the root classes are disjoint: x=0 modulo q², or x=1 modulo q.
These contribute 1+q. Hence C(m)=∏(1+q), the unitary divisor sum.
For K a power of two define

    P = ∏_{j=1,2,4,...,K/2} C(n^j)/n^j
      = ∏_{q || n, q prime power} (1−q^(−K))/(1−q^(−1)).

Thus n/P = A / ∏(1−q^(−K)), where A=∏(q−1) is an integer.
With b=bitlength(n) and K the first power of two ≥3b+2, the denominator
is at least 1−b·2^(−K). Since A≤n<2^b, n/P differs from A by less
than 1/2. Round exactly using rational arithmetic. Each p−1 divides
its own p^e−1, so nA is a multiple of φ(n). The usual randomized
factoring reduction from an exponent multiple applies, including to
non-squarefree inputs. On squarefree n, A itself equals φ(n).

This uses O(log b) exact counts of a single fixed monic cubic, on
moduli of O(b²) bits, with no mask and no random choices before the
standard exponent-based splitter. Conversely the factorization gives
C immediately. The bottleneck is counting solutions with different
local choices without knowing the prime powers. This result gives
a simpler equivalent counting target; it does not compute the count
in polynomial bit time. [Decoder](cubic_masked_trace.py)

### The mask-one trace alone suffices on arbitrary inputs

Write b=bitlength(n), K=16(b+1), F(x)=1+4x/5+x²/5. For n coprime
to 6 the fixed-mask trace T(n^j,1) equals B times the product of
F(q^(−j)) over maximal prime powers q, with B=5^ω(n). Rounding
T(n^K,1) recovers B: the excess is at most
5^b (exp(b·5^(−K))−1). Now choose c_j so that through degree K

    Σ c_j log F(x^j) = log(1−x).

If f_i is the i-th log F coefficient, f_1=4/5, and
c_j=(−1/j−Σ_{d|j,d<j} c_d f_{j/d})/f_1. The two reciprocal
roots of F have modulus 1/√5, giving |f_i|≤2·5^(−i/2)/i.
Induction gives |c_j|≤2j²: the divisor sum is bounded by
2j² Σ_{i≥2}|f_i|/i² < 2j²/4, and c_1=−5/4.
For degree m>K the absolute coefficient error is at most 3m³ by
summing over divisors. At x≤1/5 the resulting total log error over
all components is at most b Σ_{m>K}3m³5^(−m)<1/(16n)
(the successive terms have ratio below 1/2 for m≥K).
Therefore n exp(Σ c_j log(T(n^j,1)/B)) rounds to ∏(q−1).

Absolute trace errors ≤1/(128 n K³) still suffice: the baseline rounds
identically, all traces exceed one, and the total propagated log error
is at most 4K³/(128nK³). Elementary functions can be evaluated with
polynomial precision or interval bounds. The queried bases have O(b²)
bits; only O(b) absolute precision bits are required. Multiplication
of the recovered integer by n supplies a multiple of φ(n) as above.
The numerical implementation uses guard bits. No mask selection,
phase data, or factorization is an input to the decoder.

### Exact dense-interval baseline

For H consecutive bases ending at X, a segmented factor sieve gives
each mask-one trace exactly. Generate primes through floor(sqrt(X))
and remove the maximal power q of each visiting prime from each eligible
base, accumulating 5q²+4q+1. Each residual above one is prime since two
unvisited primes would have product exceeding X. Divide the accumulated
numerator by the square of the base after processing the residual.
The empty product at n=1 is one.

The prime sieve costs O(sqrt(X) log X) elementary steps as a loose bound.
Across primes, visits to the interval are at most H Σ_{p≤√X} 1/p+√X;
all successful divisions together cost O(H log X) steps. Thus total bit
work is (H+√X) poly(log X), with polynomial bit cost per output when
H≥√X and O(H log X+√X) bits of working memory. For a single requested
base, this implementation still generates the sieve through √X. The
power bases used by the trace decoder are sparse, so the dense-interval
bound does not give that decoder polynomial total work.

### Precision threshold for the fixed-mask average

At relative error 1/4 the mask-one trace oracle is deterministic
polynomial-time equivalent to computing the distinct prime count ω(n).
For the forward direction a single query at n² suffices. With
F(x)=1+(4x+x²)/5 and 0≤x≤1, log F(x)≤x. The distinct maximal
prime powers q≥5 give Σ q^(−2) ≤ Σ_{i≥5} i^(−2) ≤ 1/4 by
integration. Thus T(n²,1)/5^ω lies in [1,exp(1/4)]. Even after
25% multiplicative error the value lies strictly between
2·5^(ω−1) and 2·5^ω. Rational threshold comparisons recover ω.

Conversely start with 5^ω and multiply in the exact F(1/q) factors
for all primes ≤bD found by trial division, where b=bitlength(n)
and D is a positive integer. The missing log correction is at most
b/(bD)=1/D. The resulting underestimate has relative error
≤1−exp(−1/D)≤1/D. The work is polynomial in b and D. In
particular even inverse-polynomial relative error is available on a
promised distinct-prime semiprime (ω=2) without an oracle. This
does not supply the finer O(1/n) accuracy of the semiprime decoder.
For general integers, 2 and 3 can be stripped before querying.

### Mixed masks avoid the prime-gap loss in approximate factor recovery

For n=pr with distinct primes p<r, draw a random unit mask with Jacobi
symbol −1. Its positive local sign belongs to p with probability 1/2.
In that case the mean squared magnitude equals H(r), where

    H(z)=(5+4z/n+z²/n²)(1+1/z²).

Throughout [√n,n], H'(z)=4/n+2z/n²−10/z³−4/(nz²)≥3/n for n≥144.
Consequently an absolute trace estimate within δ locates the larger
prime to distance ≤1+nδ/3 by monotone integer binary search, including
when the estimate lies outside the range (clip to the endpoints).
The bound stays uniform even for almost equal primes. In comparison,
inverting their sum first introduces a derivative proportional to the
inverse prime gap.

For fixed 0<ε<1/4, δ≤n^(−3/4−ε) suffices for polynomial-time
factoring of every distinct-prime semiprime, given this masked average.
Use the standard univariate Coppersmith small-root theorem for a hidden
divisor ≥n^(1/2), applied to the monic linear polynomial hint+X.
The unknown root r−hint has magnitude ≤1+n^(1/4−ε)/3, within the
n^(1/4−ε/2) bound for sufficiently large n. Test candidate roots by
gcd and retry with independent masks. Sampling also exposes a factor
immediately if the mask is a nonunit. This uses an asymmetric mask,
so it does not give the same accuracy threshold for the mask-one trace.
The theorem supplies the lattice step; the implementation here only
computes the hint.

The mixed-trace wrapper now passes its interval to the repository's
trackE Coppersmith decoder with the rigorous parameter check enabled.
Tests perturbed factor-supplied traces in both directions at 82, 130
and 202 bits. Each run returned a proper divisor with one FLINT LLL
call and one integer-polynomial root call. The 202-bit hint errors were
−2147483648 and +2147483649. These runs test decoding only.

### Nine consecutive point values suffice for a semiprime trace

Strip primes ≤9 first. For distinct-prime semiprimes, expand the local
squared-magnitude product at t=j: the sequence is a sum of nine
exponentials zᵃ (at exponent a=j), each with strictly positive weight.
The nodes are n²-th roots, distinct by CRT since each component
frequency index is −1,0,1. The weight at z=1 is the desired trace.
The other nodes have average zero over unit scales because at least one
primitive local phase remains. All scales j=1,...,9 are units.

Form the ten-by-ten symmetric Toeplitz matrix with these nine samples
as off-diagonals and zero diagonal. The true moment matrix is PSD of
rank nine, so its missing common diagonal is minus the least eigenvalue
of the zero-diagonal matrix. With pointwise absolute error η the latter
matrix changes by at most 9η in operator norm. Restore the estimated
diagonal, retain the leading nine-by-nine block, and add 20ηI. The
result B exceeds the true block in PSD order and differs by at most 38η.
Then 1/(1ᵀB⁻¹1) is the minimum of vᵀBv subject to Σvᵢ=1. The true
objective equals the weighted sum of |P(z)|² for the coefficient
vector v. The polynomial ∏_{z≠1}(X−z)/(1−z) attains the exact trace.
It has real coefficients by conjugate pairing, and coefficient 1-norm
at most (n²)^8 since |1−z|≥4/n². Hence the regularized minimum lies
between the trace and trace+38η n^32 in exact arithmetic. This deliberately
loose bound still costs only O(log n) extra precision bits. In fixed
dimension inverse conditioning is polynomial in 1/η after regularization;
working precision can be increased for numerical roundoff. The supplied
implementation is numerical, not an interval certificate.

In particular the previous semiprime mask-one decoder needs only nine
absolutely calibrated squared point values of the same mask, each with
O(log n) bits of absolute precision. Computing those values natively
remains the required evaluator step.

### Five complex samples split any two-component base directly

Let n have two distinct prime-power components, both coprime to 30,
and put M=n². Their unit-scale local phases x,y are M-th roots
with xy=z=exp(2πi/M) by the CRT idempotent sum. For mask one the
sample J(j)=(x^j+a)(y^j+b) has a,b>2. Subtract z^j and difference
adjacent residuals to obtain D_j=b(x−1)x^j+a(y−1)y^j. It follows that

    D_(j+1)+z D_(j−1)=(x+y)D_j.

Five consecutive samples give both candidate divisors D₂ and D₃.
Since D₃−yD₂=b(x−1)x²(x−y), the larger absolute value is at least
16/M²: distinct M-th roots have chord distance at least 4/M.
With absolute complex sample errors η≤1/(100M⁴), each D error is
at most 2η. Using the larger observed divisor gives root-sum error
at most 8η/(16/M²−2η) < η M². Each root w of the resulting quadratic
has |w|<4 and satisfies |(w−x)(w−y)|≤4ηM². It lies within 1/(5M)
of a true local phase, enough to round its argument to the correct
M-th-root index. A gcd of this index with n gives a proper prime-power
factor. This proof allows repeated exponents on both components.
Small primes dividing 30 can be removed before use. The implementation
uses guarded numerical arithmetic and verifies the divisor by gcd.

### A single calibrated complex sample determines a two-component split

For n=qr with q,r its two maximal prime powers, put a=2+1/q,
b=2+1/r, s=a+b=4+(q+r)/n and P=ab=2s−4+1/n. Let
h=exp(−πi/n²)=c+id, z=exp(2πi/n²), and U+iV=h(J_n(1,1)−z).
The two rotated local phases hx and hy are conjugates because xy=z.
Thus U−Pc=s cos θ and V−Pd=±(a−b) sin θ. With
D=(s−4)²−4/n=(a−b)², the following quartic vanishes at the true s:

    Q(s)=D(U−Pc)²+s²(V−Pd)²−s²D.

The leading coefficient is 4(c²+d²)−1=3. This matters even when
the roots cluster. Use rational approximations to U,V,c,d accurate
to 2/(10¹² n⁴) per coordinate (sample error ≤1/(10¹² n⁴)
plus guard precision suffices). At the true s, s<5, P<5, D<1,
and the two centered coordinates are bounded by 5 and 1 respectively.
Perturbing the displayed squares gives |Q_approx(s)|<τ=10^(−6)n^(−4).
Its leading coefficient remains above 2. Wherever |Q_approx|≤τ,
the distance to some complex root is at most (τ/2)^(1/4) by
factorization. Hence each connected real component of this sublevel
set has length at most 8(τ/2)^(1/4)<1/n. It is bounded and its
ends solve Q_approx=±τ. After substituting k=n(s−4), isolate these
real boundary roots with rational intervals of width 1/8. Testing
nearby integers finds k=q+r; an integer discriminant and division
verify the split. Fixed-degree rational root isolation costs polynomial
bit time. This proof requires the absolute complex normalization of J.

### A direct root sum for the calibrated sample

Fix scale and mask one. Scan w modulo n subject to
w(w−1)(w+1)(w−2)=0 mod n. At primes ≥5 the four labels are
pairwise distinct, hence each maximal prime power has exactly four
local choices. Write g=gcd(w,n), h=gcd(w−2,n) and let E be the
CRT idempotent one modulo g² and zero modulo (n/g)², taking E=0
for g=1. The sample equals the sum of exp(2πiE/n²)/h over these w.
Expanding the established stationary-phase product proves the identity:
label zero selects the local phase, ±1 contribute unit weight and
label two contributes 1/q. All quantities in each summand are obtained
by integer arithmetic on n,w. The scan costs n poly(log n) at O(log n)
precision. Every root outside the four public residues already reveals
a proper divisor via one of the four gcds. Explicitly enumerating the
mixed roots therefore entails factor recovery.

For n=5,7,25,35,55 the compressed sum agreed within 10^(−55) with
a separate evaluation scanning x modulo n⁴ in the literal cubic
projection. At 35 and 55 those literal values also passed the one-sample
decoder, yielding 5 in both cases. This checks the complex normalization
against the original sum without using factor-supplied point values.

### Fixed-mask magnitude at constant relative precision on all bases

Suppose a deterministic oracle returns |J_n(t,1)| within 1% relative
error on unit scales. Trial divide small inputs (n≤4096) and strip 2,3.
For n with at least two maximal prime powers, put M=n² and b=bitlength(n).
The formal log magnitude f(t) on the whole cyclic group Z/M is the sum
of log|c_q+exp(2πitE_q/M)|, c_q=2+1/q. It lies between 0 and 2b.
At frequency E_q its normalized Fourier coefficient has magnitude >0.22:
the leading coefficient is 1/(2c_q)≥1/4.4, and the power series for
log(1+z/c_q) has coefficients bounded by 2^(−j)/j. Aliases from
other positive or negative harmonics begin at j≥q²−1≥24. Other local
components cannot contribute there, by the disjoint CRT coordinates.
The gcd of E_q with n equals n/q.

Use log of the oracle output as the observed function. If gcd(t,n) is
proper, return that factor immediately. If n divides t, use zero.
For analysis extend the observed function on proper-gcd inputs by the
true f. Any algorithmic query there already succeeds. The uniform
Fourier coefficient perturbation is bounded by −log(0.99)+2b/n<0.02
for n>4096. Thus each E_q has observed coefficient magnitude >0.20.
This argument assumes absolute calibration across scales.

Here is an explicit polynomial interval-energy search that does not
scan M coefficients. For a consecutive frequency interval I write
K(y)=Σ_{k∈I}exp(2πiky/M), computable by a geometric sum. For uniform
independent t,y, the expectation of g(t)g(t+y)Re K(y) is the Fourier
energy of I (normalized transform, real observed function g). Importance
sample y with mass w(y)/W, where w(0)=M and for 2^j≤min(y,M−y)<2^(j+1)
set w(y)=M/2^j. W is the sum of w over all residues. There are O(log M)
shells, uniformly samplable within each shell. The geometric-sum bound
gives |K(y)|≤w(y), and W/M=O(log M). Multiply each sampled product
by W/(M w(y)). Its magnitude is O(b³), so Hoeffding gives additive
0.005 accuracy in polynomially many samples, including the logarithmic
confidence factor. Trigonometric values need only polynomial working
precision, since geometric-sum denominators are Ω(1/M) away from y=0.

Bisect [0,M−1], retaining children with estimated energy ≥0.025.
Every heavy fundamental survives since its energy exceeds 0.04.
Each retained interval has true energy ≥0.020 on the successful event;
Parseval bounds the count at each depth by (2b+1)²/0.020. Cap the
search at this bound per level, restarting if exceeded; a union bound
over the resulting polynomial query budget ensures success with
probability at least 3/4 per run. At singletons test gcds against n.
Primes and perfect powers can be handled separately, and recursion
gives randomized expected polynomial factoring. This is a reduction
using only mask one, variable unit scales, and coarse magnitudes.
It concerns pointwise values; the fixed-mask scale average loses
the location of the frequency. No native pointwise estimator is supplied.

A vectorized finite-size check of interval-energy retrieval now samples
only unit scales from the magnitude callback. To exercise the frequency
search it ignores all sampled gcd factors, fills missing log values with
the pilot mean, and centers the samples by that mean. This different
missing-data convention is an empirical test, not the proof above.
With 32,768 correlated sample pairs per depth, fixed seed 2026, and
alternating ±1% relative magnitude perturbations, the returned frequency
and gcd were (265226,101) at n=10403, (33175177,391) at n=7429, and
(15856203,121) at n=5929. The oracle controls used the local factors;
the search used no full Fourier transform.
[Finite-size retrieval implementation](cubic_scale_spectrum.py)

### Comparing two scales at mask one suffices

Let a comparator on unit scales give the sign of f(t)−f(s), with
arbitrary fixed answers when the absolute natural-log difference is ≤0.01.
Strip primes ≤16b first (b the original bit length). Return any proper
gcd encountered in sampling. For analysis only, extend comparisons on
proper-gcd inputs by the correct sign; pairs involving a multiple of n
can receive a zero score. Use an independent fair argument swap and
negate the reversed result. The resulting mean pair score A(t,s) is
antisymmetric and bounded by one. Define g(t)=E_s A(t,s) on uniform
scales in Z/n². Then E g=0 and

    E[f(t)g(t)] = (1/2) E[(f(t)−f(s)) A(t,s)].

For one local log factor the range is ≤log 3, while the variance is
≥2(0.22)² by its two fundamental coefficients. Its independent-copy
absolute difference therefore has mean ≥4(0.22)²/log 3>0.17.
Summing independent CRT components cannot decrease this mean by
conditional Jensen. Charging 0.02 for arbitrary near-tie signs and
at most 4b/n for pairs involving multiples of n leaves covariance
>0.06 on remaining composite bases (or trial divide n≤4096).

The sum of absolute nonconstant Fourier coefficients of f is ≤b log 2,
by its log power series. A local harmonic jE_q always yields a proper
frequency gcd unless q divides j. The total nonconstant Fourier mass
from such bad harmonics is ≤b·2^(1−min q)<0.001 after trial division.
Parseval and |ĝ|≤1 therefore imply a frequency with proper gcd and
|ĝ|≥0.05/b. The interval-energy search above applies with energy
tolerance a fixed small multiple of 1/b² and polynomial sampling cost.
Each occurrence of g in a sample product uses a fresh comparator pair
with independent reference and swap, preserving its expected energy
even at zero shift. Verify every frequency gcd and restart on failure.
This proves a randomized polynomial reduction from factoring to a
one-bit comparison across scales of a fixed mask with no padding.

The finite-size unit-only comparison experiment used one randomized
comparison per score, with false returned throughout the 0.01 dead
band. At the previous fixed budget and threshold it recovered 101
from 10403 and 121 from 5929; the run on 7429 returned no divisor.
That empirical cutoff is larger than the worst-case theorem permits.

### The second scale can always be twice the first

Put h(t)=f(t)−f(2t) on Z/M. Doubling permutes this odd cyclic
group, so h has zero mean. At E_q, the contribution from f(2t)
requires a local harmonic 2j≡±1 mod q². All such harmonics have
|j|≥(q²−1)/2≥12, with combined absolute coefficient bounded by
Σ_{j≥12}2^(−j)/j<0.001. Other CRT components cannot contribute.
Thus |ĥ(E_q)|>0.219 and E|h|≥|ĥ(E_q)|. The sign comparison g(t)
between t and 2t obeys E[hg]≥E|h|−0.02 apart from unavailable
scales. As before, return a factor on a proper gcd, extend those
inputs correctly for analysis, and set g=0 on multiples of n.
The last change costs at most 2b/n. For n>4096 the correlation
is therefore greater than 0.19.

The nonconstant Fourier coefficient 1-norm of h is ≤2b log 2.
After stripping primes ≤16b, the coefficient mass at nonzero
frequencies failing the gcd test is at most 2b·2^(1−min q);
doubling preserves the gcd. By Parseval, at least one proper-gcd
frequency of g has magnitude >0.13/b. The same polynomial
interval-energy search now consumes direct bits and requires no
reference-scale sampling. All oracle calls use units and mask one.
This reduces arbitrary factoring to deciding the direction of change
under doubling, allowing arbitrary fixed replies in the 0.01 natural-log
dead band.

With the same finite-size budget and seed as the random-reference
test, the doubling-bit experiment returned (frequency,gcd) pairs
(265226,101), (23787012,323), and (3440633,49) at 10403, 7429,
and 5929, respectively. All near-ties returned false. These tests
used factor-generated bits and disabled early factor recovery on
nonunit scales.

### A threshold list in the public scale needs many intervals

Take distinct primes 101≤p<q, n=pq, M=n². Extend any correct doubling
bit to signs g on the whole scale circle, with arbitrary choices on
nonunits and inside the 0.01 natural-log dead band. As above E|h|>0.219.
Each local difference is bounded by log 3, so |h|≤2 log 3.
Replacing the exact sign on the exceptional inputs gives

    E[hg] ≥ 0.219 − 0.02 − 4 log(3)(1/p+1/q) > 0.11.

Every nonconstant frequency of h is a multiple of p² or q². Its
distance from zero modulo M is therefore ≥p². This remains true
under **any invertible affine reordering** t=au+v modulo M; such
reordering also preserves the Fourier coefficient 1-norm, ≤4 log 2
for h. If g has K cyclic sign changes in that ordering, its total
variation is 2K. Taking a first difference in the normalized Fourier
transform bounds |ĝ(k)| by 2K/(M |1−exp(2πik/M)|), at most
K/(2p²) on the support of h. Parseval now gives K>0.079p².
For balanced primes with q≤2p this implies more than 0.039n changes.
Thus even an input-dependent invertible stride cannot yield a
poly(log n)-length interval list for the promised comparisons. This
is a representation bound; an arithmetic program may oscillate
frequently at low cost.

Direct small-input checks counted mandatory cyclic transitions after
discarding all nonunits and dead-band entries, giving the minimum over
all assignments to those entries. In natural order the counts were
7926, 3938, and 18402 for n=143,323,899. At stride n+1 they were
1214, 43926, and 366058. These are illustrative cases below the
prime cutoff in the uniform proof.

### Squaring inside the subgroup of scales congruent to one

Take t=1+nk modulo n², k modulo n. Then t²=1+n(2k) modulo n²,
so a single prescribed comparison between t and its square stays in
the same subgroup. No unit tests or unavailable scales occur. Strip
primes below 17 and handle single prime powers separately. On Z/n
the log magnitude F(k) has each local term
log|c_q+α_q exp(2πi k E_q/n)|, with public c_q given q and
an unknown unit phase α_q. The local character has exact order q.
The shift α_q preserves the magnitudes of the individual power-series
coefficients before aliasing.

For H(k)=F(k)−F(2k), the Fourier coefficient at E_q mod n has
magnitude >0.22: its leading term has modulus ≥1/4.4, the other
harmonics at that index in F start at q−1, and the harmonics in
F(2k) start at (q−1)/2≥8. Their combined tail is at most
Σ_{j≥16}2^(−j)/j + Σ_{j≥8}2^(−j)/j <0.002. Distinct local
components have disjoint nonzero supports. In this shorter group
every nonzero local harmonic has proper gcd with n, because each
is a multiple of n/q and the only multiple having gcd n is zero.
H has zero mean, so there is no bad-frequency subtraction.

The comparison bit G (encoded ±1), even with arbitrary choices when
|H|≤0.01, satisfies E[HG]≥E|H|−0.02>0.20. Since the nonconstant
Fourier 1-norm of H is at most 2b log 2, Parseval gives a proper-gcd
frequency of G of magnitude >0.14/b. The earlier cyclic interval
search applies on n points. This proves equivalence to expected
polynomial factoring for the restricted squaring comparison oracle.
No claim of a native implementation of that oracle is made.

The finite-size experiment, with factor-generated comparisons and the
same sample budget, recovered frequency/divisor pairs (103,103),
(646,323), and (2058,49) at n=10403,7429,5929 respectively.

### Fixed shifted Jacobi probes for the lift comparison

For distinct primes 17≤p<q and n=pq consider a predictor whose entire
input about k is m ternary symbols (k+a_i | n), with fixed offsets.
Allow arbitrary postprocessing depending on n. Remove duplicate offsets
mod n; if a pair collides at just one prime, their difference already
splits n by gcd. Otherwise the offsets are distinct at both primes.

Write the comparison log difference H(k)=h_p(k mod p)+h_q(k mod q).
Both local differences have mean zero and absolute value ≤log 3.
Restrict the postprocessing to the Boolean cube {−1,+1}^m and use
its multilinear extension on the full cube [−1,1]^m. This extension
is bounded by one by interpolation. Its Walsh coefficient 1-norm is
at most 2^(m/2) by Parseval and Cauchy–Schwarz. A nonempty monomial
on a subset S is the product of two independent prime-field characters
of the polynomial ∏_{i∈S}(k+a_i). The standard Weil bound for a
quadratic character of a squarefree polynomial bounds its normalized
mean at prime r by (|S|−1)/√r. In the correlation with h_p the q
coordinate is averaged, and conversely for h_q. The empty monomial
has zero correlation. On a zero-symbol input the true postprocessing
may differ from its multilinear extension by two. A union bound then
gives the following bound on the absolute correlation with H:

    log(3) 2^(m/2) (m−1)(1/√p+1/√q) + 4 log(3) m(1/p+1/q).

The promised comparison requires correlation at least E|H|−0.02>0.20.
In particular any such realization for balanced large primes requires
m ≥ log₂ p − O(log log p). This only constrains a fixed list of
linear Jacobi probes, and a logarithmic count would still fit the
desired runtime.

An exhaustive check grouped k by its full ternary symbol vector at
offsets 0 through m−1. Summing the absolute within-group sums of H
and dividing by n gives the best possible correlation over all tables.
For n=10403 the optima for m=1,3,6 were <0.000001, 0.044044, 0.139384
against a required 0.644228. For n=256027 they were <0.000001,
0.005579, 0.012504 against 0.549090. The local factors were used to
score predictors, and no training or random sampling was involved.

### Adaptive Jacobi trees and the cost of decoding their answers

Consider an explicit ternary decision tree with S internal nodes and
depth d. Each node selects a fixed shifted Jacobi probe, and its edges
encode −1,0,+1. Offsets can depend on the path. Duplicate probes on
a path can be removed, and two offsets on a path congruent at exactly
one prime already split n by gcd of their difference.

Ignore zero edges temporarily and use the binary subtree. Its signed
output has a multilinear extension given by a sum over at most S+1
leaves of signed path indicators ∏(1±X_i)/2. Each path expansion has
coefficient 1-norm at most one and degree at most d. At ternary
inputs the sum is still bounded by one: interpret a zero as a fair
binary branch. If any of the S probes is zero, the original tree may
differ by two; the union of these events has probability at most
S(1/p+1/q). Applying the preceding CRT/Weil estimate term by term
bounds the absolute correlation by

    log(3)(S+1)d(1/√p+1/√q) + 4 log(3) S(1/p+1/q).

Consequently for balanced large primes (S+1)d=Ω(√p) is necessary
for a valid comparator in this explicit-tree model. This gives a
size barrier even though only one root-to-leaf path runs per input.
The argument does not cover trees generated on demand by succinct
arithmetic rules or circuits with shared computation.

There is also an upper bound for distinguishing inputs alone. For two
different scales k,l, pick a prime r where they differ. For a uniform
shift a the two nonzero local Legendre signs have product +1 on
(r−3)/2 choices and −1 on (r−1)/2; this follows from the exact
quadratic character sum −1. At the other prime s both symbols are
nonzero with probability at least 1−2/s. Therefore the global ternary
symbols differ with probability at least (1−2/s)(r−3)/(2r)>1/3 for
p,q≥17. For m independent shifts the chance of any indistinguishable
pair is at most n(n−1)(2/3)^m/2. Thus O(log n) probes can encode
every k injectively. This alone supplies no efficient rule for the
desired labels. Exhaustive checks also found injectivity for the first
23 consecutive offsets at n=10403 and the first 31 at n=256027;
these counts refer to the full ternary vectors, without target labels.

### One lifted magnitude against a public power of two

The comparison need not involve another scale. On a cofactor n with
at least two prime-power components, take b=bitlength(n) and strip
primes ≤32b (using the pre-stripping bit length is also sufficient).
Use F(k)=log|J_n(1+nk,1)| as above. For r=ω(n), the constant
coefficient of F−r log 2 has absolute value <0.02. Indeed each local
mean differs from log(2+1/q) by at most Σ_{j≥q}2^(−j)/j (only
multiples of q survive), while log(1+1/(2q))≤1/(2q).
The sum is bounded by b/(2 min q)+b·2^(1−min q)<0.02.

For the correct threshold exponent r the fundamental coefficient
of F−r log 2 at E_q mod n has magnitude >0.22 by the earlier alias
bound. Its mean absolute value is therefore >0.22. For any threshold
bit G∈{−1,1} correct outside a 0.01 natural-log dead band, it follows
that E[(F−r log 2)G]>0.20. The constant mode contributes at most
0.02. All remaining frequencies of F have proper gcd with n, with
total coefficient 1-norm at most b log 2. Hence some proper-gcd
frequency of G has magnitude >0.18/(b log 2). The cyclic interval
search applies. Try exponents 2,...,b with capped polynomial searches
and verify returned divisors; repeat a sweep if necessary. One trial
at the correct exponent succeeds with constant probability, irrespective
of the outcomes at other exponents. For a promised two-component base
only the public threshold 4 is required after the small-prime step.
This reduction presupposes absolute magnitude calibration.

The experimental sweep started at exponent two without receiving the
component count. With seed 2026 and false answers throughout the dead
band it returned (exponent,frequency,divisor) = (2,5150,103),
(3,323,323), (2,2058,49) for the same three inputs. The factor-generated
oracle controls used the original unstripped inputs in these examples.

### A public-angle predictor for the semiprime threshold

For n=pq and a lift t=1+nk the two local phases x,y have known
product exp(2πit/n²). Consider the predictor that replaces 2+1/p and
2+1/q by 2, then models one phase as uniform on the circle conditional
on this product. Put C=cos(πt/n²) and D=Re(x exp(−πit/n²)). In this
model the squared magnitude minus 16 is exactly
16D²+40CD+16C²−7, with D distributed as the cosine of a uniform angle.
Thus the model probability of exceeding four is one minus the
arcsine-law mass between (−5C±sqrt(9C²+7))/4, clipped to [−1,1].
This is a factor-free heuristic whose conditional distribution is a
model assumption.

The majority prediction is positive for C²<1/34 or C²>7/16. To see
the inner transition, when both roots are interior and straddle zero,
the intervening angular width equals π/2 when their squared sum equals
one. This gives 17C²/4+7/8=1. After the outer root clips at an endpoint,
the other crosses zero at C²=7/16. Symmetry handles negative C.
Therefore the bit changes sign at most four times per cycle. By the
first-difference bound on Z/n its normalized Fourier amplitude at any
proper-gcd index is ≤2/min(p,q). This decays exponentially in bit
length for balanced factors, too small for the polynomial-budget
frequency search.

Enumerating all lifts and scoring only outside the 0.01 natural-log
dead band gave 54.976% accuracy at n=10403 and 55.097% at n=256027,
versus 51.931% and 51.493% for always predicting positive. The
largest proper-gcd Fourier amplitudes were 0.007977 and 0.001290.
Thus improvement over the majority baseline in these examples did
not reproduce the factor-bearing signal required by the reduction.

### Nonlinear polynomial preprocessing of a public-angle threshold

Let n=pq be a distinct-prime semiprime. Let B be a sign sequence on
Z/n with at most K cyclic changes, and define g(k)=B(P(k) mod n)
for an integral polynomial P of degree 2≤d<min(p,q) with unit leading
coefficient. A nonunit leading coefficient can be checked by gcd, or
discarded when it vanishes mod n. For every nonzero target frequency j,

    |ĝ(j)| ≤ K(1+log n)[(d−1)²/√n + 1/p + 1/q].

Indeed expand B in normalized additive Fourier modes a. Mode zero has
zero contribution at j≠0. For unit a the normalized complete sum with
phase aP(k)−jk factors by CRT into two prime-field additive character
sums, each bounded by (d−1)/√r by Weil: the linear term cannot cancel
the leading term. Bound each nonunit mode trivially by one. The
first-difference bound |Bhat(a)|≤K/(2 dist(a,nZ)) gives total nonzero
mass at most K(1+log n); summing over nonzero multiples of p or q
separately gives at most K(1+log n)(1/p+1/q). This proves the claim.

Thus on balanced semiprimes even the largest individual coefficient
is exponentially small in bit length for polylogarithmic K and d.
This excludes this class of modular preprocessing from meeting the
inverse-polynomial spectral requirement. The argument counts polynomial
degree, not circuit size; iterated squaring can reach much higher
degree cheaply.

For the four-transition public-angle predictor, exhaustive FFTs gave
the following largest proper-gcd amplitudes:

| n | identity | k² | k³+3k+7 |
|---:|---:|---:|---:|
| 10403 | 0.0079771 | 0.0289617 | 0.0205733 |
| 256027 | 0.0012903 | 0.0049959 | 0.0069722 |

The changes on these two inputs illustrate why a modest finite-size
gain alone is insufficient evidence for a polynomial-time reduction.

### Repeated squaring: an invariant independent of formal degree

Let n=pq be squarefree and define g_s(k)=B(k^(2^s) mod n), for any
fixed real function B. Under uniform CRT inputs X,Y, the total squared
Fourier amplitude at proper-gcd frequencies is exactly

    Var(E[g_s | X]) + Var(E[g_s | Y]).

Indeed the nonconstant characters on each single coordinate form an
orthonormal basis for its centered conditional expectation. In particular
this sum is unchanged by any separate bijections on X and Y.

At prime r the power map of exponent e has one zero preimage of zero
and has gcd(e,r−1) preimages of each point in the nonzero image, the
unique subgroup of that index. Two exponents having the same gcd with
r−1 therefore give the same map up to a permutation of the input
(matching fibers over each output). Consequently the proper-gcd energy
is exactly the same for all s≥max(v₂(p−1),v₂(q−1)). This applies to
any B, with no requirement that the individual Fourier coefficients
remain unchanged.

For p,q≡3 mod 4 this saturation starts at the first squaring. For a
K-transition predictor the degree-two bound gives each nonzero
coefficient after that first step amplitude at most
A=K(1+log n)(1/√n+1/p+1/q). There are p+q−2 proper-gcd frequencies.
At every subsequent depth each such coefficient thus has amplitude
at most √(p+q−2) A by the invariant. On balanced Blum semiprimes this
is O(K log(n)/n^(1/4)), exponentially small in bit length for
polylogarithmic K, uniformly over all squaring depths.

Using the public-angle threshold as B, the computed proper-gcd energy
at depths 1 through 9 stayed at 0.0173125452 for 103·107 and
0.0035651609 for 503·523. For 97·113 the energy rose from
0.0169629122 to 0.3753766464 at depth five and stayed there through
depth nine, as predicted by v₂(96)=5 and v₂(112)=4. This result is
specific to a single modular power followed by the predictor.

### Arbitrary modular powers on safe-prime semiprimes

Let n=pq with p=2r+1 and q=2s+1 distinct safe primes, r,s odd.
For any supplied positive exponent e, write u for its odd part.
If r or s divides e, a uniform unit a gives a proper divisor via
gcd(a^u−1,n) with probability **exactly one half**. At a covered
prime the equation a^u=1 selects exactly half the units, since
gcd(u,2r)=r (or gcd(u,2s)=s). The outcome at the other prime is
independent by CRT. Exactly one success thus has probability one half
regardless of that other prime's selection probability. This test uses
only n and e, and has polynomial bit cost when e has polynomial bit
length. It requires no knowledge of which case applies.

If neither odd half-order divides e, the e-th power map has the
same local fiber profiles as identity for odd e and squaring for
even e. The preceding conditional-variance invariant therefore gives
exactly the identity or first-square total proper-gcd energy for any
fixed output function B. For a K-transition threshold list and balanced
safe primes, the amplitude of every factor-bearing frequency is bounded
by O(K log(n)/n^(1/4)) uniformly in e, as above (for identity use
the first-difference bound and sum over the p+q−2 frequencies).
Thus any exponent achieving inverse-polynomial factor-bearing amplitude
beyond this bound already gives a separate factoring shortcut. This
is a per-instance statement and assumes no infinitude of safe primes.

On n=47·59 all 2668 units were checked. Exponents 1 and 5 had identical
proper-gcd energy 0.00104350 for the public-angle predictor. Exponents
2 and 10 both gave 0.03365720. For 23, 58 and 10672 the direct gcd
split exactly 1334 units in each case, even though the predictor for
10672 had zero proper-gcd energy. The calculation supplies no method
to construct an exponent containing a hidden half-order.

### Random affine output thresholds versus direct circuit collisions

For a distinct-prime semiprime n=pq let A(k) be any computable map
acting separately on CRT coordinates (in particular any fixed modular
arithmetic circuit). Let B be a fixed sign sequence. Choose an
independent uniform unit a and uniform offset v and read B(aA(k)+v).
Write W_* for the squared normalized Fourier mass of B on unit
frequencies, W_p for frequencies of order p, and W_q for those of
order q. Let c_r be the probability that two independent uniform
inputs to the local map A_r have the same output, and put
β_r=(r c_r−1)/(r−1). Then the expected total energy at frequencies
with proper gcd is exactly

    W_* (β_p+β_q−2β_pβ_q) + W_p (1−β_p) + W_q (1−β_q).

To see this, expand B in additive characters. Averaging over v
removes cross terms. For one character depending on both output
coordinates, let z_r be the mean of its local phase. The two
conditional variances sum to
(1−|z_p|²)|z_q|²+(1−|z_q|²)|z_p|². The two local multipliers from
a are independent. Orthogonality gives mean |z_r|² over nonzero
local multipliers equal to β_r. For a single-coordinate character,
the other local phase is constant. This yields the displayed formula.

Directly sampling independent x,y and taking gcd(A(x)−A(y),n)
succeeds with probability D=c_p+c_q−2c_p c_q. Since
|β_r−c_r|≤1/r and the function x+y−2xy is 1-Lipschitz in
each coordinate on the unit square, the expected energy is at most
D+1/p+1/q+W_p+W_q. For a list with K cyclic sign changes the
first-difference bound gives W_p+W_q≤(π²/12)K²(1/p²+1/q²).
Consequently on balanced semiprimes, an inverse-polynomial probability
of an inverse-polynomial Fourier amplitude under random affine output
thresholding implies inverse-polynomial direct collision success for
polylogarithmic K. No degree or iteration-depth bound on A is needed.
This averaged statement leaves deliberately chosen exceptional output
transforms open.

At n=77, enumerating all 4620 affine output transformations for a
one-interval indicator encoded ±1 reproduced the formula to 10^(-12)
for A=k, k², and (k²+1)²+1. Their mean energies were 0.02428740,
0.19885529 and 0.35945353; direct split probabilities were 0.20779221,
0.34677011 and 0.46753247.

### Squarefree products and the repeated-prime distinction

For squarefree n with any number of prime factors the random-affine
energy identity has a short product form. For a source frequency a,
put u_r=1 when r divides a and u_r=β_r otherwise, with β_r defined
above for the local circuit at each r|n. Its mean fraction of energy
at proper-gcd frequencies is

    1 − ∏ u_r − ∏(1−u_r).

After translation averaging, sum this against |Bhat(a)|². To prove
it, take the tensor product of the local Fourier expansions. In each
coordinate the zero-frequency fraction is the squared phase mean,
and its complement is the total nonzero fraction. Subtract the
all-zero and all-nonzero cases and average independent unit multipliers.
For uniform input pairs the direct splitting probability is the same
expression evaluated at the local collision probabilities c_r.
For unit a the two expressions differ by at most Σ1/r: the product
expression is the probability of a mixed outcome for independent
Bernoulli variables, so couple each coordinate at cost |β_r−c_r|.
For nonzero nonunit source frequencies bound the mass of a K-transition
threshold by (π²K²/12)Σ1/r² via a union bound. Thus the earlier
collision diagnostic extends to squarefree n whenever these explicit
correction terms are small enough. At n=105 exhaustive affine averages
for the same three circuits were 0.05297052, 0.49573221 and 0.81734240,
matching the product identity to within 10^(−12).

A different phenomenon occurs on repeated prime powers. The readily
computable entire function A_n(k)=k^n mod n has **least additive period
rad(n)**. For q=p^e||n, if k≡l mod p and both are units then the
ratio lies in the kernel of reduction modulo p, a group of order
p^(e−1) dividing n. Hence k^n=l^n mod q. For nonunits both powers
vanish mod q. This proves rad(n) is a period. Conversely any period d
satisfies d^n=A_n(0)=0 mod n, forcing every prime divisor of n to
divide d. In particular every scalar readout of A_n has Fourier support
on multiples of n/rad(n). If n is nonsquarefree every nonzero supported
frequency has proper gcd. The mass can still be spread over rad(n)−1
frequencies or the readout can be constant. No polynomial classical
algorithm for detecting the hidden period follows from support alone.
Enumeration confirmed the exact periods 35,77,35,385 at n=175,539,
245,385. Squarefree inputs give only the trivial full period.

### Exact gcd law for two power-map outputs

For A(x)=x^n mod n, q=p^e||n and d_p=gcd(n,p−1), the condition
p|A(x)−A(y) already implies q|A(x)−A(y). Indeed nonunits at p
map to zero modulo q. Unit outputs have order dividing p−1 since
their (p−1)-th powers have exponent n(p−1), divisible by φ(q).
Reduction mod p is injective on this prime-to-p torsion. Thus every
output gcd is a unitary divisor, selecting whole prime powers.

For independent uniform x,y the probability of selecting q is
[1+(p−1)d_p]/p² by counting zero pairs and equal nonzero p-field
power images. For uniform x and y=x+1 the probability is (d_p−1)/p:
away from 0 and −1 the ratio (x+1)/x runs over the nonzero residues
except 1. For either experiment these events are independent over
primes. With their probabilities denoted by c_p, a proper split has
probability exactly 1−∏c_p−∏(1−c_p).

If all d_p=1 the gcd is just the product of all p^e with p dividing
x−y. One can obtain the very same answer from gcd(n,gcd(x−y,n)^b),
b=bitlength(n), without evaluating A. This includes n=p^e q^f for
distinct odd primes p<q≤2p: neither divides the other's predecessor
(q−1 cannot equal p), so both residue power maps are permutations.
Adjacent inputs never split these bases. Repeated powers can therefore
give a genuine shorter period without improving this collision test.

Exhaustive checks on n=35,55,175,245,539,775 covered all input pairs.
At 175 and 245 independent pairs split with probability 2/7 and
adjacent pairs never split. At 55, where 5 divides 11−1, adjacent
pairs split with probability 4/11. These rates are for the indicated
experiments and do not bound other uses of the map.

### Factorization-free reference table on all lifted scales

Each term (E,h) of the four-label root scan contributes
exp(2πiE/n²)/h times exp(2πikE/n) at scale 1+nk. Accumulate the
first multiplier in bin E mod n and apply a positive-sign unnormalized
DFT (n times the numerical inverse FFT). One scan and one length-n
transform produce all n complex samples in O(n poly(log n)) work at
O(log n) bits of target precision with appropriate guard bits. The
experimental table uses ordinary double precision. It is a batch
reference; the linear-size setup is exponential in input bit length.

Feeding only comparisons from this table into the lifted-squaring
decoder gave (frequency,divisor)=(103,103) at n=10403, (646,323)
at n=7429, and (2058,49) at n=5929. No supplied factorization or
root-scan gcd was passed to the decoder. A separate comparison with
the local factor formula over every lifted scale had maximum absolute
discrepancy below 1.5·10^(-14) on each input. Constructing the table
does encounter factor-revealing roots, so this cannot establish an
advantage over directly taking their gcds.

### Short random-stride windows of calibrated complex lift values

For a promised n=pq with distinct primes and mask one, consider the
mean of J_n(1+n(k+jh),1) over 0≤j<L. Choose the stride h uniformly
among units. This is a smoothed target: compute 4 plus the mean of
the public phases exp(2πi(1+n(k+jh))/n²) by a geometric sum. This
estimate takes polynomial bit work in log n, log L and requested
precision and does not receive the factors.

For 0<ε,δ≤1, strip primes ≤16/ε first; finding one already splits
the promised semiprime. Otherwise the exact constant term
(2+1/p)(2+1/q) differs from 4 by <ε/2. The two remaining mixed
characters on the lift index have orders p and q, and their weights
sum to <5. A random unit stride makes each local increment uniform
on the nonzero residues of that order. Except with probability ≤δ,
both increments are at distance at least δ/6 from an integral turn:
the probability for one prime r is at most 2(δ/6)r/(r−1)≤δ/2.
For such a stride each mixed window mean has modulus at most
1/(2L·δ/6) by the geometric-sum bound. Thus L≥30/(εδ) gives
absolute error ≤ε simultaneously at every start k. For fixed ε,δ,
the sufficient window length is independent of n. This supplies no
pointwise estimator or guarantee for an adversarial stride. In fact
the returned formula has only the public constant and public frequency
one as functions of the start, so the earlier spectral decoder gets
no proper-gcd frequency from it directly.

A numerical check at n=5003·7001 used L=4096, ε=δ=0.1 and 1000
random unit strides with random starts (seed 2). One error exceeded
0.1; the median was about 0.001095 and the maximum 1.04945 rounded
up. The comparison used the factored exact geometric-sum formula.

### Stable point recovery from windows all of the same length

For any complex sequence f on Z/n let A_s(k) average the L values
starting at k with step s. For every stride h there is an exact identity

    f(k) = A_(Lh)(k) + Σ_{j=1}^{L−1} [A_(jh)(k) − A_(jh)(k+jh)].

The bracket telescopes to [f(k)−f(k+Ljh)]/L. Summing the brackets
cancels all but the desired point from the first mean. It costs
2L−1 queries, and absolute window errors η amplify by at most
(2L−1)η. This order in L is necessary for arbitrary sequences:
if all queried progressions have L distinct points, an impulse of
height Lη at k changes every mean by at most η. A zero answer is
then compatible with both that impulse and the zero sequence.
No analogous necessity claim is made for the restricted cubic family.

For unit-stride-only evaluators, gcd(j,n) for 1≤j≤L either supplies
a split or certifies the multipliers are units when L<n. If all
but a δ fraction of unit strides have the error guarantee at every
start, a uniform h makes all L strides good with probability at least
1−Lδ. Hence inverse-polynomial accuracy in polylogarithmic windows
can suffice for constant pointwise accuracy at an appropriately
small stride failure rate. The preceding public geometric estimator
cannot be substituted at this precision: applying the identity to
its outputs gives exactly 4+exp(2πi(1+nk)/n²), still omitting both
hidden components.

The formula passed 2540 exact rational tests, including wraparound
and degenerate strides. Additional cubic checks reconstructed every
lifted sample at n=35,175,385 from length-four windows of the
factorization-free reference table.

### Four window queries per point-difference trial

For real f on Z/n, take gcd(L,n)=1. For a uniform anchor w and each
endpoint u choose s=(w−u)/L mod n. Then

    L[A_s(u)−A_s(u+s)] = f(u)−f(w).

Use the same w at both endpoints and subtract. Four window calls
determine f(u)−f(v), with error ≤4Lη if each window has error ≤η.
Each stride is marginally uniform, although the strides are dependent.
If at most a δ fraction of strides are unreliable and every other
stride works at all starts, the trial succeeds with probability at
least 1−2δ. For δ<1/4, independent anchors and a median reduce the
failure probability to exp(−2m(1/2−2δ)²) with m trials. The number
of windows per trial is independent of L. This proof concerns means
of the target sequence itself. All lifts 1+nk are units regardless
of the stride between their indices. If the window oracle restricts
strides to units, count nonunits as additional bad strides and return
a factor on any proper stride gcd, or count the trial as failed.

Apply this to f(k)=log|J_n(1+nk,1)|, u=k, v=2k. At η≤0.002/L
the difference error is ≤0.008, sufficient outside the earlier
0.01 comparison dead band. Choose the repetition count for a union
bound over the polynomial Fourier search. Thus the all-input reduction
can tolerate a fixed bad-stride fraction bounded away from 1/4.
The earlier native complex-window estimate is for a different observable
and supplies neither these log means nor the required accuracy.
A test used the factorization-free n=385 table with L=8, replaced
38 entire strides by adversarial signed values of magnitude 1000,
and perturbed every good window by ±0.002/L. With 15 anchor trials
per pair all 385 differences met the 0.008 error bound (seed 73).

### How far a constant log-window baseline gets

Let n=pq for distinct primes both greater than max(17,L); smaller
primes can be trial-divided when L is polynomial in the bit length.
For a uniform unit stride, the public estimate log 4 is within
(24/δ+2)/L of every length-L log-magnitude mean, except for a
fraction δ of strides. For each local order r and positive harmonic j
not divisible by r, its increment in turns is uniform over the nonzero
r-grid. The probability it lies within δ/(12j²) of an integer is at
most 3δ/(12j²). Union over both primes and all j costs less than δ.
The geometric mean multiplier for the remaining strides has magnitude
at most 6j²/(Lδ). Multiplying by local harmonic weights ≤2^(−j)/j
and summing gives 24/(Lδ) for the two components. Harmonics divisible
by r contribute at most 2^(1−r) per component; together with the
constant-term correction Σ log(1+1/(2r)) this costs at most 2/L.

The order 1/L conceals a decisive constant. At a fundamental of order
p the pointwise coefficient exceeds 0.22 in magnitude. A length-L
mean multiplies it by |sin(πLa/p)|/(L|sin(πa/p)|), where a is a
uniform nonzero residue when the stride is uniform over units. Since
L is invertible mod p, at least half the strides put La mod p in
[p/6,5p/6], giving magnitude >0.11/L. Subtracting any estimate
constant in the start, even one chosen separately for every stride,
does not change that coefficient. Its magnitude lower-bounds the
worst-start error. Thus no such predictor can meet the 0.002/L
window target on the requisite majority of strides.

For n=101·103, 200 sampled unit strides (seed 2026) and all starts,
the median worst-start errors of log 4, multiplied by L, were
1.5156, 1.7101, 2.0754 at L=8,32,64. Even the optimal constant for
each sampled stride (midrange of the actual windows) exceeded
0.11/L on every sampled stride at all three lengths. Values came
from the factorization-free lifted reference table.

### Average-case log windows suffice, even with no wholly reliable stride

For each remaining n suppose the length-L log mean is within 0.002/L
on at least 99.9% of uniformly random (start,stride) pairs in (Z/n)²,
with gcd(L,n)=1. The bad set can be fixed adversarially. In the
four-window construction use uniform k and independent uniform anchor
w, with endpoints k and 2k. Each of the four queried pairs is uniform:
doubling is bijective, (u,w)→(u,(w−u)/L) is bijective, and shifting
the start by the stride is a shear. Hence the probability of any bad
window is at most 0.004. Only one anchor is needed per comparison.

Write H=F(k)−F(2k) as before, after removing primes below 17, and
σ²=E H². The local differences are independent, centered, bounded
by log 3, and each has variance ≥2·0.21² from its conjugate fundamental
coefficients. Therefore E H⁴ ≤3σ⁴+(log 3)²σ² <17σ⁴. Interpolating
the first, second and fourth moments gives E|H|≥σ/√17. On good
window queries the estimated difference is within 0.008; its sign
loses at most 0.016 in the correlation with H. On the remaining
event Cauchy–Schwarz bounds the loss by 2σ√0.004. Thus for the
random-anchor sign G,

    E[HG] ≥ σ/√17 − 0.016 − 2σ√0.004 > 0.06σ,

using σ≥√2·0.21. The earlier Fourier argument applies to the
conditional mean g(k)=E[G|k]. Its interval-energy estimator queries
two independent anchors for the product, including at equal points.
This average-case criterion is per modulus; it provides no allowance
for discarding hard moduli.

With factorization-free n-point reference tables, length eight and
seed 2026, the comparator was tested with error ±0.002/L except on
(k−s mod n)<floor(n/1000), where the oracle returned 1000. Every
stride had bad starts. The empirical interval search returned factors
103,323,49 at n=10403,7429,5929 respectively.

### Using average-case windows at an unfactored coprime multiple

Let M=nm with m coprime to n, M coprime to 6 and L invertible
modulo M. For spectral input k mod n take a uniform lift u=k+nz
mod M. Choose v by the CRT conditions v=2k mod n and v=u mod m.
Use the shared uniform anchor modulo M to estimate F_M(u)−F_M(v).
The logarithmic contributions of *all* prime powers of m cancel,
since each depends on the lifted index only through its local residue.
No prime factors of m need be supplied. At each q||n (a maximal
prime power) the two remaining arguments depend on k and 2k through
the same local unit rotation and primitive q-th-root character: the
global lift is 1+Mu and the local projector is one modulo q².
Changes to the coefficient caused by m are unit multipliers, so the
fundamental coefficient, moment bounds, and single-coordinate Fourier
support arguments above continue to apply on Z/n.

For uniform k, u and v are individually uniform mod M (CRT and
invertibility of 2 mod n). The preceding four joint-uniform window
queries and 0.004 union bound therefore still apply, now at M. The
interval search runs on n, so the verified divisor splits the original
target. If an efficiently sampled coprime padding of polynomial bit
length gives a good window oracle on M with inverse-polynomial
probability for each composite target, cap every attempt, verify
divisors and retry fresh paddings. Bad moduli may return arbitrary
bounded-time answers. This is a conditional weakening of the per-base
requirement; high success density among unrestricted integers does not
imply high density among multiples of each prescribed n.

Using the same corrupted-window experiment with factorization-free
reference tables, (target,padding,verified divisor) was (10403,35,101),
(7429,25,323), and (5929,65,121). The unpadded regression at 10403
returned 103. Only the integer padding was passed to the search.

### Window recovery on an arbitrary multiple, including powers

For any M=dn coprime to 6 choose one uniform z mod d at the start of
the search. Compare u=z+dk and v=z+d(2k mod n), using the same
four-window anchor construction modulo M. If p^a||n and p^e||M,
the local lift phase before restriction is a unit twist times a
primitive character of order p^e in the lift index. Since v_p(d)=e−a,
its restriction in k has exact order p^a and is supported on that
one CRT coordinate of n. For p dividing M but not n the restriction
is constant and cancels. The local additive offset in the cubic
product is now 2+1/p^e; it obeys the earlier bounds on the fundamental,
aliases and local range. The independent centered local differences and
proper-gcd spectral support therefore apply to every fixed z after the
same small-prime preprocessing of the target. No decomposition of d
or splitting via gcd(d,n) is needed.

To transfer the 99.9%-good window promise, average first over z and
uniform k. Both u and v are then individually uniform modulo M.
With the independent uniform anchor each requested window pair is
uniform as before, assuming gcd(L,M)=1. The mean conditional
probability of any bad call is ≤0.004. Markov gives probability at
least 1/3 of an offset with conditional failure probability ≤0.006.
On each such offset, uniformly over k and the fresh anchor, the
correlation lower bound is σ/√17−0.016−2σ√0.006 > 0.03σ,
using σ≥√2·0.21. Thus a bounded search followed by verification
and fresh-offset retries gives the same polynomial reduction. The
offset must remain fixed within each frequency search.

With linear-size reference tables and the earlier corrupted length-eight
windows, targets and bases (323,323²), (7429,17·7429) and (323,85·323)
yielded verified divisors 17,323,17 respectively. The first test is
the important case where gcd(d,n)=n conveys no split. The unpadded
10403 regression returned 103. The oracle is still required to be
usable on the selected multiple, whose bit length must be polynomial
in the target length.

### Average-case obstruction when one start coordinate is missing

Let q≥17 be a maximal prime-power component of n and gcd(L,n)=1.
Even a predictor given the full stride and every start coordinate
except k mod q must fail the 0.002/L tolerance on more than 39% of
uniform (start,stride) pairs. It may use arbitrary computation on
the available data. This concerns predictors whose output, conditional
on those data, does not depend on the missing coordinate.

Grant the exact local window contributions of every other component
and subtract them. For fixed remaining inputs the estimate for the
q-window is now constant in its starting index i. Write X(i) for
the one-component log magnitude and s for the stride reduced mod q.
If the predictions at i and i+s both meet tolerance, telescoping
implies |X(i)−X(i+Ls)|≤0.004. Here the auxiliary start coordinates
are kept fixed in both tests; only the omitted coordinate is shifted.
The two test inputs have the same marginal uniform distribution.
Moreover (i,i+Ls) is uniform on independent residues mod q since
L is invertible. Conditioning on all auxiliary start coordinates is
harmless; in the final probability bound s is averaged, not fixed.

For X=log|c+exp(iθ)| with c=2+1/q and θ on a rotated q-grid, any
interval of log width 0.008 imposes a cosine interval of width at most
(c+1)²(1−exp(−0.016))/(2c)<0.04. The largest circular probability
of any cosine interval of this width is arccos(0.96)/π (attained
at an endpoint of the cosine range). Its preimage is at most two arcs,
so the grid probability is at most arccos(0.96)/π+2/q<0.21.
Conditioning on the independent first value gives probability <0.21
that the two log values differ by at most 0.004. By a union bound
the single-window failure probability exceeds (1−0.21)/2. This
also rules out any stride-only baseline under the joint average-case
criterion; the proof does not assume that a stride works at all starts.

Exhaustive finite checks used the direct lift table, L=8 and the
optimal predictor for each (stride, visible-start-coordinates) class:
sort the true values and take the most populated interval of width
0.004/L. Its success rates were 9025/104329 for n=323 with q=17
hidden, 10693/152881 for n=391 with q=23 hidden, and 148225/1713481
for n=1309 with q=17 hidden.

### Losing one digit of the start at a repeated prime

Suppose q=p^e exactly divides n, with p≥17 and e≥2, and gcd(L,n)=1.
Give the predictor the full stride s and the start k modulo n/p. Under
the uniform distribution on (k,s), its failure probability δ at tolerance
0.002/L satisfies

    δ ≥ (2/3 − 1/q)(1/3 − 1/p)/2 > 0.09.

Indeed, subtract the window at k+s from that at k and multiply by L.
The exact difference telescopes to F(k)−F(k+Ls). Both queries have
uniform marginals, so the predicted difference has error at most 0.004
except with probability at most 2δ. Fix s and k modulo n/p. The
prediction and every other prime-power contribution are now fixed.
For the remaining term write X(θ)=log|c+exp(iθ)|, c=2+1/q. Across
the p possible starts, θ runs over a rotated p-grid. The angle increment
β for Ls is fixed in this class and uniform on the q-grid over s, since
both L and the local frequency are units.

At least 2/3−1/q of strides have |sin(β/2)|≥1/2. For any such stride,
each of sin(θ+β/2)≥1/2 and sin(θ+β/2)≤−1/2 holds on at least 1/3−1/p
of the starts in the class, by counting grid points in a single arc.
The difference of squared moduli is 4c sin(θ+β/2) sin(β/2), so it
has opposite signs on the two sets and magnitude at least c. Taking
half the logarithm preserves the sign and gives magnitude at least
c/[2(c+1)^2]>0.10 by the mean value theorem. Any fixed prediction
therefore misses one set by more than 0.004. Averaging proves the bound.
Randomized predictors obey the same bound by fixing their random seed.

In particular rad(n) divides n/p here. The period of k^n modulo n
therefore prevents a predictor using only that start statistic from
meeting the 99.9% target, even with unrestricted stride information.
For L=8, exhaustive direct-table checks maximizing the number of hits
in each visible class gave 7736/83521 successes at n=17², 11158/130321
at n=19², and 193900/2088025 at n=5·17² (hiding the top digit at 17
in the last case).

### A single Jacobi correlation separates the repeated-prime case

For odd n>1 let C(n) be the unnormalized sum over k modulo n of the
Jacobi symbol (k(k+1) | n). It obeys the exact identity

    C(n) = ∏_{p^e || n} p^(e−1) · b(p,e),
    b(p,e) = −1 for odd e, and p−2 for even e.

For even e the local symbol is one except at the two roots modulo p.
For odd e the Legendre sum is −1: the number of pairs (k,y) with
y²=k(k+1) is p−1, since (2k+1−2y)(2k+1+2y)=1 parametrizes them by
one nonzero field element. Each k modulo p lifts p^(e−1) times, and
CRT multiplies the local sums. Thus 0<|C(n)|<n. On squarefree n the
value is ±1. Otherwise n/rad(n)>1 divides C(n), so gcd(n,C(n))
is a proper factor even if cross-prime cancellation enlarges the gcd.

An exact polynomial-time correlation evaluator together with a
promise-squarefree factoring algorithm would suffice on all integers.
Remove powers of two, apply the gcd, recurse on a split, and invoke
the promise algorithm only when the gcd is one. The recursion tree
has O(log n) nodes, counting prime factors with multiplicity.
This is a conditional reduction: direct summation still uses n terms.
For n=p²q with distinct odd primes, C(n)=−p(p−2), hence
1+sqrt(1−C(n))=p. Even an unnormalized additive error strictly less
than p−5/4 suffices to round the resulting estimate of p correctly
by the spacing of squares around (p−1)². For balanced p and q this
allows normalized error of order n^(−2/3), still exponentially small
in input bit length. Direct sums verified both the product identity
and the gcd criterion on all 251 odd integers from 3 through 503.

### Squaring the input and varying the quadratic

Even the fixed C has information on squarefree semiprimes when queried
at a square. If n=pq for distinct odd primes, C(n²)/n=(p−2)(q−2).
Thus p+q=(n+4−C(n²)/n)/2, which recovers both factors from one call.
This is an exact-query reduction, with a modulus of twice the bit length.

For arbitrary odd n one can instead vary the discriminant. Define
A(n,d) as the count of x modulo n with gcd(x²−d,n)=1. Equivalently it
is the sum of the Jacobi symbol (x²−d | n²) over x modulo n² divided
by n. For unit d, let S be the primes dividing n at which d is a
quadratic residue. Counting the two roots at each such prime gives

    A(n,d)/n = ∏_{p in S} (p−2)/p.

Its reduced denominator h=n/gcd(n,A) divides n and contains no prime
outside S. If the largest prime P dividing n belongs to S, its
denominator factor cannot cancel, since every numerator factor is
positive and less than P. If another prime dividing n lies outside S,
h is a proper divisor. For uniform unit d the Legendre signs are
independent by CRT. Thus for n with at least two distinct prime
divisors, imposing these two signs at the largest and smallest primes
gives a split with probability at least 1/4, regardless of exponents.
Sample d uniformly from 1 through n−1 and take any nontrivial gcd
immediately. This preserves the unconditional lower bound. Standard
primality and perfect-power tests handle prime powers, and recursion
gives expected polynomial time given a polynomial exact count oracle.
Conversely the factors give the count in polynomial time for unit d.
For nonunit d each prime dividing d contributes the density (p−1)/p
in place of the two-root or no-root density.

Direct counts over every unit d gave proper reduced denominators on
6/8, 18/24, 27/36, 30/40, 42/48, 90/120, and 210/240 cases at
n=15,45,63,75,105,225,385. The count is still the unknown operation:
testing a given x with a gcd takes polynomial time, while enumerating
all x takes n tests. No evaluation algorithm for the original window
or for this exact count follows from the oracle equivalence.

### Approximate quadratic counts on two distinct odd primes

Let n=pq with p<q. Choose a uniform d of Jacobi symbol −1 by rejection
sampling (a nontrivial input gcd already splits n). Exactly one prime
admits roots. In the orientation where d is a square modulo p, there
are 2q excluded residues modulo n. Hence A=n−2q. The two orientations
are equally likely by CRT, and neither the selection nor the decoder
needs the factorization. Given an approximation B, set H to the nearest
integer to (n−B)/2. The monic linear polynomial x+H has a root q−H
modulo the unknown divisor q≥sqrt(n).

The univariate Coppersmith theorem for unknown divisors of size at
least n^β recovers roots up to n^(β²−ε) for a monic linear polynomial
and fixed ε>0. Taking β=1/2, any additive count error ≤n^(1/4−ε)
suffices for all sufficiently large n (0<ε<1/4), including rounding.
Check gcd(n,H+x) for each returned integer root and repeat on another
d if necessary. An always-accurate oracle needs at most two selected
queries in expectation. This applies also to unbalanced semiprimes.
In normalized density the sufficient additive error is n^(−3/4−ε).
This bound concerns the two-prime promise; the exact denominator
argument for general n does not give the same noise tolerance.

Direct example: n=1009·2017=2035153 and d=5 give A=2031119 by counting
all residues. Perturbing A upwards by 6 gives H=2014. For X=5, LLL
on the coefficient rows of n, H+x, x(H+x), with column i scaled by
X^i starting at i=0, returned x²−7x+12 after unscaling the first row.
Its roots 3 and 4 gave gcds 2017 and 1009. This small example tests
the decoding step; the asymptotic exponent uses the general theorem.

### Approximate counts suffice for arbitrary multiplicities

The n^(1/4−ε) absolute count tolerance extends to all odd composites
with at least two distinct prime divisors. Put R=rad(n), let S be
the primes at which a random unit d is a square, and put

    b(S) = (n/R) ∏_{p outside S} p,
    A(n,d) = b(S) ∏_{p in S}(p−2).

Then b(S) divides both n and A. Complementary sign patterns satisfy
b(S)b(Sᶜ)=n²/R≥n. Consequently at least one member of each pair has
b≥sqrt(n). To ensure the gcd is proper, restrict to pairs whose two
largest primes P>Q have opposite signs. In either member let r be the
largest prime in S. It equals P or Q. Its exponent in b is one short
of its exponent in n, and none of the positive factors p−2 with p in
S reaches r. Thus gcd(n,A)<n in both members. In fact this properness
argument applies to every nonempty S; it only uses the largest prime
*in S*. Opposite signs simply guarantee both members are nonempty.

More sharply, discard only the complementary pair consisting of the
empty and full sets. Each remaining pair contributes at least one
good pattern. With t≥2 distinct primes, a uniform unit d therefore
gives sqrt(n)≤gcd(n,A)<n with probability at least 1/2−2^(−t)≥1/4.
Independence and uniformity of signs follow from CRT. This is an
existence and probability argument; the decoder never forms S or R.

Round an approximate count to H. Apply the unknown-divisor univariate
small-root algorithm to H+x modulo n with β=1/2. On a good pattern,
x=A−H is a root modulo gcd(n,A). For fixed 0<ε<1/4 the count error
n^(1/4−ε) plus rounding fits n^(1/4−ε/2) for sufficiently large n.
Test gcd(n,H+x) for every returned root, accepting only proper factors.
Sample d uniformly in 1,...,n−1; nonunits split immediately. An
always-accurate polynomial-time evaluator gives at most four trials
in expectation per split. An evaluator exceeding tolerance on a
fraction δ<1/4 of uniform unit queries still gives unit success at
least 1/4−δ; a constant gap suffices. Remove powers of two, test
primality and perfect powers, and recurse to obtain complete factors.

Direct enumeration of all x and all unit d gave large proper gcds
on 36/48, 108/144, 120/240, 360/720, 270/480 queries respectively at
105, 315, 385, 1001, 1155. The theorem relaxes the accuracy demanded
of the evaluator; it does not supply the count.

### Precision under powering the input modulus

Write D(n,d)=A(n,d)/n. For M=n^k and d'≡d mod n, D(M,d')=D(n,d)
exactly: the coprimality test is identical and every residue has
M/n lifts. This gives a precision reduction without any lattice step.
Suppose an evaluator runs in polynomial time in the bit length L of
M and achieves normalized additive error 2^(−L^α) for fixed α>0.
For b=bitlength(n) choose k=(b+2)^ceil(1/α), so L≥k and L^α≥b+2
for n≥2. Constructing M takes polynomial time in b, and rounding n
times the resulting density recovers A(n,d) exactly. A uniform unit
query modulo M projects to a uniform unit modulo n, so the earlier
random-query split also tolerates a constant failure fraction below
its good-pattern probability. This reduction requires the precision
guarantee on prime-power lifts of the input, even when n is squarefree.

Conversely any fixed inverse polynomial normalized accuracy in b
already has a deterministic polynomial-time evaluator. For target
b^(−c), trial-divide through T=ceil(2 b^(c+1)) and multiply the
exact local densities for the small prime divisors. For odd p the
root count modulo p is 1, 2 or 0 according as d vanishes, is a
nonzero square, or is a nonsquare; at 2 there is one root.
Omitting larger primes overestimates D by at most the sum of their
root densities, bounded by 2b/T since n has at most b prime divisors.
The product is rational of polynomial bit length. Thus this coarse
evaluator works for every d. Polynomial input powering cannot turn
its guaranteed inverse-polynomial error into the <1/(2n) rounding
precision. Intermediate rates, for example 2^(−(log b)^2), are left
open by these two arguments. Exhaustive enumeration checked the
powering identity at k=2 for every residue of 15,21,25,35,45,77.

### Resolving individual primes from intermediate precision

Let n be odd, b its bit length, and assume a pointwise additive density
error η≤1/[16(2b+1)(B+1)²], for an integer cutoff B≥3. A randomized
polynomial-size batch recovers every prime factor ≤B with probability
at least 3/4. Choose m with 4(b+1)²≤2^m<8(b+1)², sample m independent
uniform units g_j, and query at all subset products indexed by binary
vectors x. For each p|n let v_p record the nonresidue bits of the g_j
modulo p. By CRT these vectors are independent and uniform. A union
bound over pairs and zero vectors gives probability <1/4 of any
collision or zero. Units can be sampled by rejection with expected
polynomial cost (the unit density is at least the density bound below).

For a unit d, each local density is either 1 or 1−2/p. Thus the
log density on this binary cube is a constant plus a sum of Walsh
characters (−1)^(v_p·x), each with coefficient w(p)=log(1−2/p)/2.
On the no-collision event the normalized Walsh transform isolates
each coefficient. Sorting the distinct odd primes gives p_i≥2i+1,
so D≥∏_{i=1}^t (2i−1)/(2i+1)=1/(2t+1)≥1/(2b+1).
Taking logs therefore incurs error at most 2(2b+1)η, and averaging
with signs does not increase that bound. Budget another equal amount
for numerical precision. Since w'(z)=1/[z(z−2)], candidate integer
values from 3 through B have spacing greater than 1/(B+1)².
Find the nearest candidate by binary search; for a true p≤B it
is uniquely correct. Check primality and divisibility before keeping
any candidate and divide repeatedly for multiplicity. These checks
exclude spurious answers; repeating reduces the miss probability.
The arithmetic costs are polynomial in b and log(1/η) for B≤n.

For the concrete intermediate rate η=2^(−(log₂ b)²), polynomial
powering of the modulus boosts the constant in the exponent arbitrarily:
the bit length of n^k is Θ(kb). Consequently a polynomial-time
pointwise evaluator at this rate would extract every prime up to
2^(C(log₂ b)²) for any fixed C with high probability. Conversely,
given all primes through 2b/η, multiplying their exact local factors
gives absolute error ≤η by the earlier omitted-prime bound. This
locates the intermediate target at small-factor extraction with a
quasipolynomial cutoff, and leaves balanced large-prime inputs open.

A direct enumeration at n=38675=5²·7·13·17 used 128 subset products,
independent uniform additive perturbations of size ≤10^(−6), and
recovered 5,7,13,17 from the transformed log densities (seed 41).

### Random-input errors in the small-factor decoder

Let t be the number of distinct odd prime divisors and let T≥t be a
known upper bound (the bit length b works). Take a power of two K
with 4T≤K<8T. Use log₂ K random unit generators as above, but start
the subset products at an additional independent uniform unit u.
Every queried argument is now individually uniform, including the
empty-subset query. Translation multiplies the coefficient of each
prime by its Legendre sign at u, preserving its absolute value.

For any fixed target prime the probability its binary vector is zero
or coincides with another prime's vector is at most t/K≤1/4. Other
primes may collide freely. If the evaluator violates the earlier
accuracy tolerance on fraction δ≤1/(32T) of units, the probability
of any bad value in the batch is at most Kδ≤1/4, even for correlated
query locations. Hence with probability at least 1/2 the absolute
transformed coefficient identifies this target by the same spacing
bound. Verify primality and divisibility for every decoded candidate.
With independent batches a fixed target is missed with probability
at most 2^(−r) after r batches; r=ceil(log₂(T/ρ)) suffices for all
targets simultaneously with probability ≥1−ρ. This argument permits
arbitrary errors on a fixed bad set and requires no knowledge of its
locations. It does not give constant error tolerance for unbounded t.

For an arithmetic check on n=38675, mark every 256th element of the
ordered unit list as bad (90/23040 units) and replace its log density
by −10. Using the local product formula for other entries, six batches
of 32 translated products recovered all four primes in 100/100 seeded
trials (seed 61). This checks the decoder, not an implementation of
the density evaluator.

### Euler witness split in the certified density interval

If trial division leaves an odd cofactor R with at most two prime
factors counted with multiplicity, a negative Jacobi query excludes
the square case. Test gcd(R,u±1) for u=d^((R−1)/2) modulo R. A proper
gcd gives both primes, so the local density is exact. At R=1009·2017
and d=5 the gcd with u−1 is 1009 and the density is 1007/1009.

For R=pq and a uniform unit query, put h=gcd(p−1,q−1). A split by
this test requires u=±1 at some prime, equivalently d^(R−1)=1 there.
In the cyclic unit group modulo p this has probability h/(p−1),
and analogously at q. Hence the unconditional split probability is
at most h/(p−1)+h/(q−1); conditioning on negative Jacobi costs at most
a factor of two. Small h on balanced inputs limits this particular
witness test even with repeated random queries.

### Deterministic queries using Schemmel totients

Let S_k(n) count starts of k consecutive units modulo n, with S_0=n.
These are the Schemmel totient functions. At a prime p such a run
excludes min(k,p) starting residues, hence its local count at p^e
is p^(e−1) max(p−k,0). Let b be the original input bit length.
Remove all prime factors ≤b by trial division. If a cofactor m>1
remains, then throughout k=0,...,b one observes the single polynomial

    S_k(m) = (m/rad(m)) ∏_{p|m} (p−k).

Its degree is at most b. Exact interpolation and deterministic rational
polynomial factorization recover all its roots; division of m recovers
multiplicities. Coefficient bit lengths are polynomial (the sum of
their absolute values is at most m·2^b). This deterministic reduction
uses at most b nontrivial count queries of polynomial run length.
Conversely the prime factors give every query by the local formula.
It is a different evaluator family from the quadratic density: the
tested product has k linear factors. Direct enumeration followed by
interpolation for 7429 and 5491 returned respectively
−(k−17)(k−19)(k−23) and 17(k−17)(k−19).

There is a simpler distinction for binary-encoded long runs: S_k(n)>0
exactly when k is smaller than the least prime factor of n. Thus a
zero test with cost polynomial in log n and log k factors deterministically
by binary search and recursion. Equivalently that test decides whether
gcd(k!,n)=1. This identifies the cost hidden by a succinct specification
of the run; iterating over k residues would not meet that bound.

### Rounding noisy short-run counts before interpolation

Fix any polynomially bounded integer H≥1 as a function of the original
bit length b. Trial-divide through H+b and let m be the cofactor.
For H≤k≤H+b each surviving prime exceeds k. The count polynomial
P(k) has integer coefficients and P(0)=m, hence S_k(m)≡m mod k.
Each estimate Y with |Y−S_k(m)|<k/2 gives the exact count by

    m + k · round((Y−m)/k).

The b+1 counts therefore recover the factor polynomial even with
adversarial raw errors strictly below H/2 throughout. This short-run
reduction still needs normalized accuracy of order H/m: choosing a
polynomial H only saves O(log b) precision bits for a b-bit cofactor.

Corrupted locations can also be tolerated. For s arbitrary bad
queries, use b+2s+1 consecutive k starting at H, increasing the trial
cutoff to the largest k. After rounding, solve over the rationals
for a monic E of degree s and a polynomial Q of degree at most b+s
such that Q(k)=E(k) times the rounded value at every query. Such a
solution exists by taking E zero on the errors (and adding factors
to reach degree s). Every solution has Q=EP since Q−EP has at least
b+s+1 good zeros. Exact division recovers P. This is the standard
Berlekamp–Welch interpolation argument with integer count rounding.
Runtime is polynomial for polynomial s,H and polynomial-bit oracle
answers. Enumeration at 101·103, k=20,...,34, with alternating
additive errors ±9, recovered (k−101)(k−103).

### A divisor from one leading coefficient of a quadratic run count

For unit d, count starts x modulo m for which every x²−d j² is a
unit for 1≤j≤k. Call this Q(k). If the least prime of m exceeds 2k,
the excluded residues at a split prime p are the 2k distinct elements
±j√d; at an inert prime there are none. Write S for the primes where
d is a square and t=|S|. CRT gives

    Q(k) = c ∏_{p in S}(p−2k),     c=m/∏_{p in S}p.

For any positive integer H with 2tH below the least prime, the t-th
forward difference on 0,H,...,tH equals (−2H)^t t! c. If each raw
count has absolute error less than H^t t!/2, normalization and
nearest-integer rounding recover c: the absolute binomial weights
sum to 2^t. The degree need not be supplied. Try orders 1 through
bitlength(m) and accept only candidates that divide m properly.
Only the true order needs to obey the precision and prime-size
conditions. With polynomial H one may first strip primes through
2H times the input bit length. A mixed sign pattern guarantees
1<c<m, and uniform unit d gives such a pattern with probability
1−2^(1−ω(m)) for at least two distinct primes. The stated tolerance
depends on the realized t; for mixed signs on a semiprime t is one.

Direct enumeration at m=7429=17·19·23, d=30 and H=3 produced
7429,3289,805. The second difference 1656 divided by 72 gives 23.
Even with raw errors (0,−8,+8), this rounds to 23. This reduction
uses a family of joint counts. Its k=1 member is the earlier
quadratic count; separate marginal evaluations do not give these
joint counts.

### Binary exponents as a proposed replacement for long root lists

Consider a balanced safe-prime semiprime n=pq, p=2r+1, q=2s+1,
with r,s distinct odd primes. For unit d the binomial x^E−d has
either zero or gcd(E,p−1) roots mod p, by cyclicity, and likewise
at q. Unless r or s divides E, the total excluded fraction on
uniform x modulo n is at most 2/p+2/q. The same bound sums over any
fixed list of such binomials, without requiring independent root sets.
Thus polynomially many tests without a resonant exponent still have
exponentially small excluded density as a function of bit length.

Conversely if either r or s divides E, write M for its odd part.
For a uniform unit a the chance of a^M=1 mod p is gcd(M,2r)/(2r).
At a prime where its half-order divides E this equals 1/2. The two
prime events are independent by CRT, so exactly one holds with
probability precisely 1/2, even if both half-orders divide E.
Consequently gcd(a^M−1,n) splits directly, without a count oracle.
No factorization is needed to strip the twos or verify a split.
At p=47,q=59, enumerating all 2668 units gave 1334 splits each
for E=46,58,1334, and 102 for E=16. This is a restriction on
amplification by binomials and uniform sampling, not on arbitrary
polynomial circuits or exact evaluation of tiny density differences.

### Dickson trace tests and an exponent-halving split

The Dickson polynomial V_E with parameter one satisfies
V_E(z+z^{-1})=z^E+z^{-E}. For odd p the fraction of x modulo p
satisfying V_E(x)=2 is exactly

    a_p(E) = [gcd(E,p−1)+gcd(E,p+1)]/(2p).

Indeed z lies in the split multiplicative group or the norm-one
subgroup over the quadratic extension. The trace equals 2 after
powering exactly when z^E=1. In each cyclic group the trace quotient
pairs inverse elements, except the fixed points (1 and, for even E,
−1). The two groups share exactly those fixed kernel points, so
their exceptional terms cancel in the combined count.

Suppose n=pq, with distinct odd primes. If max(a_p(E),a_q(E))≥δ
for 0<δ≤1/4, test all indices from the odd part of E up to E by
successive doubling. At the odd initial index both fractions are
at most 1/2, since each gcd is at most half the even group order.
Each doubling increases either fraction by at most a factor of two.
At the first index with maximum ≥δ, that maximum A is thus at most
1/2 (including if the first index is the odd part). Write B for the
other fraction. A uniform x modulo n yields a proper gcd of n with
V_index(x)−2 with probability A+B−2AB≥A≥δ by CRT. One may use
the same random x throughout a pass: the successful index is fixed
by n,E,δ. Compute the initial trace by binary Lucas powering, then
iterate V_{2j}=V_j²−2. Each pass uses O(log E) modular operations.

This allows p+1 smoothness as in Williams' p+1 method. More generally,
if the full exponent excludes fraction ε of residues (root at either
prime), then the maximum local fraction is ≥ε/2. A detectable
inverse-polynomial exclusion rate therefore suffices to split directly
on these semiprimes. Constructing suitable polynomial-bit exponents
on arbitrary inputs remains open here. Enumeration checked the local
formula for every odd prime below 101 and every E from 1 to 128.

### Bit-length-only exponent lists miss almost all balanced prime pairs

Let L=log X and fix A,B,D>0. Consider any list of positive integer
exponents of total binary length ≤L^A, chosen independently of the
primes in [X,2X]. Except for an O(L^(−D)) fraction of primes p in
this interval, **every** listed exponent E has both gcd(E,p−1)/(p−1)
and gcd(E,p+1)/(p+1) smaller than L^(−B). Constants are uniform over
the list. No assumption on safe primes or shifted-prime smoothness
is needed.

To prove it, take C>max(2,B,A+D+2) and y=L^C. The exponent product
has at most L^A distinct prime divisors. For each such divisor ℓ>y,
at most O(X/y+1) integers p in the interval have ℓ dividing p−1 or
p+1, by counting progressions. Remove their union. Next remove p
for which either neighbor is y-smooth. For σ=1−1/(2C), Rankin's
bound for the number of y-smooth positive integers ≤3X is

    (3X)^σ ∏_{ℓ≤y prime}(1−ℓ^(−σ))^(−1)
      ≤ (3X)^σ exp(O_C(y^(1−σ))) = X^(σ+o(1)).

Indeed the Euler product bounds the sum of m^(−σ) over smooth m;
its logarithm is bounded by a constant times the sum of j^(−σ)
over integers 2≤j≤y. Dividing the two exceptional counts by the
number of primes in the interval, asymptotic to X/log X, gives
O(L^(−D)). For each survivor, each neighbor has a prime divisor
above y absent from every E. Its gcd ratio is therefore at most
1/y < L^(−B), as claimed.

The exact trace root formula is a weighted average of the two gcd
ratios, so it obeys the same bound. For two random distinct primes
in this interval the exception probability grows by at most a factor
of two. Their product has only constantly many possible bit lengths;
taking the union of the corresponding precomputed lists changes only
the length constant. Any fixed polynomial number of sampling tests
can be absorbed by increasing B. Random lists independent of the
input primes obey the same distributional result by conditioning.
This rules out an input-size-only library as a universal source of
detectable trace roots. It does not address exponents computed from
the actual integer n, nonuniform choices of x, or other evaluators.

### Exponents given by polynomial functions of the public modulus

The preceding distributional obstruction extends to E=|F(pq)|>0.
Fix a list of nonzero integer polynomials chosen per bit length, with
sum of degrees, number of entries and total bit length of their
contents all ≤L^A, L=log X. The content is the positive gcd of the
coefficients. For random distinct primes p,q in [X,2X], any fixed
B,D>0, the probability that any listed E has a Dickson trace root
fraction ≥L^(−B) at either prime is O(L^(−D)). Polynomials with
zero value on an input supply no positive exponent there.

Indeed apply the preceding smooth-number/progression exclusion to
the product of the contents. Outside the exceptional p set, choose
one prime ℓ>y=L^C from each of p−1,p+1, neither dividing any
content. Conditional on p, each F(pq) modulo such ℓ is F(q) or
F(−q) respectively. These are nonzero polynomials over the field,
so at most the sum of degrees of the entire list of residue classes
of q can make any value zero. Each class has O(X/ℓ+1) integer
representatives in [X,2X]. Dividing by the number of eligible primes
q gives conditional exception probability O(L^(A+1)(1/y+1/X))
per witness. Choose C>max(2,B,A+D+2). Outside both exceptions, each
neighbor of p has a missing divisor >y for every exponent. Both gcd
ratios are then ≤1/y. Swap p and q and union bound for the claim.
The proof permits selecting a sublist using n, and taking any positive
divisor of a nonzero listed value cannot increase a gcd. Hence it
includes (n±1)/2, products of polynomially many short shifts, and
halving schedules for these exponents. It does not cover general
polynomial-time functions of n, such as input-dependent modular
exponentiation producing an exponent outside the list.

### Exponents centered on the Fermat estimate of the totient

Take E_c=|n+1−2 ceil(sqrt(n))+c|>0 for offsets in any preselected
list of at most L^A integers per bit length, L=log X. The preceding
negligible-density conclusion still holds on random distinct primes
p,q in [X,2X], even though rounding the square root puts E_c outside
the polynomial family. This includes testing all polynomially many
successive Fermat estimates.

Here are the extra counting details. Discard |p−q|<X/L^K. Counting
integers q for fixed p and dividing by the number of interval primes
costs O(L^(1−K)). As above, outside a negligible exceptional set of
p, choose a prime ℓ>L^C dividing each neighbor, separately. For
ℓ|p−1 use the plus sign and for ℓ|p+1 the minus sign in

    f(q) = ±q + 1 − 2 sqrt(pq) + c.

If ℓ divides E_c then f(q) is within 2 of a multiple of ℓ. Its
range has length O(X), independent of the offset size. The minus
branch has derivative of magnitude at least one. On each remaining
side of p the plus branch is monotone with derivative magnitude
|1−sqrt(p/q)|=|q−p|/[sqrt(q)(sqrt(q)+sqrt(p))] ≥1/(4L^K).
Thus each multiple admits O(L^K) integer preimages within the
width-four window, and only O(X/ℓ+1) multiples need be considered.
After summing over offsets and normalizing by prime count the bound
is O(L^(A+K+1)(L^(−C)+1/X)). For any prescribed B,D>0 choose
K>D+2, then C>A+K+B+D+3. Using the smooth-number estimate above
and swapping p,q gives exception probability O(L^(−D)) for any
local trace hit rate ≥L^(−B). This permits arbitrary selection among
listed offsets but does not cover an unrestricted offset computed
from n. When ceil(sqrt(n))=(p+q)/2 the zero-offset guess is exactly
φ(n); in that case the first Fermat square test already recovers the
factors via ((p+q)/2)²−n=((q−p)/2)².

### Random elliptic curves with a polynomial-bit exponent budget

In the uniform-curve, uniform-point model, even a positive exponent
list depending arbitrarily on n=pq has negligible probability of a
local annihilation on random balanced primes, if picked before the
curve and of total bit length ≤L^A, L=log X. Here p,q are distinct
uniform primes in [X,2X]. Locally draw a uniform nonsingular short
Weierstrass model over F_p and then a uniform point on it. For every
D>0 the probability that any listed scalar kills the point is
O(L^(−D)). Polynomially many trials follow by union bound. This
statement targets the specified scalar test and sampling law.

Write M for the curve order. Hasse places it within 2√p of p+1.
For ordinary traces the Deuring class-number formula bounds the
probability of any one order at fixed p by O(L²/√X). Indeed each
isomorphism class of short models is weighted proportionally to
1/#Aut, and the resulting mass for a trace is bounded by a constant
times H(4p−t²)/p. The quadratic order class-number bound
h(Δ)=O(√|Δ| log |Δ|) and summation over square divisors give
H(d)=O(√d log²(2d)). For p>3 the only possible supersingular trace
is zero; handle its order p+1 separately below.

Take y=L^C. On averaging over p, the probability that an ordinary
order is y-smooth is at most O(L³ Ψ(3X,y)/X): each candidate order
occurs in Hasse intervals of at most O(√X) integer centers, there
are asymptotically X/L eligible primes, and the individual-order
probabilities obey the bound above. The preceding Rankin estimate
makes this negligible for fixed C>2. The supersingular case needs
care because the exponent list can depend on p through n: use the
standard trace-zero part of Deuring's formula, whose mass is bounded
by O((h(−4p)+h(−p))/p)=O(L/√X), omitting invalid discriminants.
Thus all trace-zero draws may simply be charged as exceptional.

Condition now on the primes and their preselected list. Its product
has at most L^A distinct prime divisors above y. For each such ℓ,
the Hasse interval contains O(√X/y+1) multiples. The ordinary
single-order mass bound makes the probability of meeting any of
these divisors O(L^(A+2)(1/y+1/√X)). Otherwise a nonsmooth M has
a prime divisor ℓ>y absent from every exponent. Multiplication by
each scalar is invertible on the ℓ-primary subgroup, so annihilation
requires the point to project to zero there, of probability ≤1/ℓ.
The same necessary event applies simultaneously to every listed
scalar. Choose C>A+D+3, and sum the bounds for the two hidden primes.
The argument leaves curve-dependent exponent choices, specially
constructed curves or points, and other factoring signals open.

### Point-first sampling and intermediate affine denominators

The preceding negligible-probability estimate transfers to sampling
uniform a,x,y modulo n and setting b=y²−x³−ax. At a prime p>3 each
short model with its marked affine point has exactly one generating
triple. The p singular coefficient pairs have at most 2p affine
points each, so singular reduction costs at most 2/p. Conditional
on nonsingularity, the sampler is uniform over model/affine-point
pairs. Compare to the preceding model-first law (including infinity):
the likelihood ratio for an affine pair on a model of order M is
M divided by the mean affine size over all nonsingular models.
Hasse bounds this by (p+1+2√p)/(p−2√p), uniformly O(1).
Thus local rare-event bounds transfer without assuming that the
two samplers are identical.

A chain starting from one marked point P has formal integer labels
i for [i]P. At an ordinary affine addition of labels i,j, equality
of x coordinates means [i]P=±[j]P, hence [i−j]P=O or [i+j]P=O.
At a doubling, a zero denominator means [2i]P=O. Infinity at an
intermediate nonzero label is itself an annihilation. Exclude
identities of formal integers: adding P to itself uses doubling,
and adding opposite labels gives infinity identically. A proper
denominator gcd must therefore be preceded or witnessed by a
nonzero local annihilator from these labels and their signed sums.
For s additions/doublings starting with label 1, magnitudes are at
most 2^s; the O(s) tested integers have total O(s²+s) bits.

The same first-event proof permits branching only on these equality,
infinity and denominator tests. Fix n and independent scheduling
coins, symbolically follow the path on a point of infinite order,
and form its label list before the curve is drawn. Up to the first
unexpected local relation the actual path agrees. The previous
bound applies directly to that list; there is no union over the
whole branch tree. This covers polynomial-length single-point
affine addition chains under the stated sampler on random balanced
prime pairs. Choices using the numerical coefficients or coordinates,
other curve families, and general factorization algorithms are outside
this argument.

### Coordinate-dependent selection among many scalar chains

Let L=log X, W=sqrt(L log L). Under the previous sampling assumptions,
allow O(W) bits of decisions derived arbitrarily from numerical
coordinates or model coefficients, even adaptively, solely to select
the scalar addition schedule. Formal equality/infinity branches need
no extra bits before the first unexpected local relation. Conditional
on n and independent coins, enumerate all coordinate decision strings
and collect the candidate labels on their formal paths. For chains
of polynomial length their union has total binary length ≤exp(KW)
for a fixed K. The union need only exist in the proof. Conditioning
on the chosen string itself would incorrectly bias the curve law.

For this enlarged menu the probability of any local annihilation
is exp(−Ω_K(W)). Set y=exp(aW) with a>K+1 and take
σ=1−hW/L, h=1/(4a). Rankin's Euler product now has logarithm at most
O(y^(1−σ)·sum_{ℓ≤y prime}1/ℓ)=O(L^(1/4) log log y)=o(W)
by the prime reciprocal sum bound and σ tending to one. Consequently
Ψ(3X,y)/X ≤ exp(−hW+o(W)). Averaging over Hasse intervals and using
the single-order mass bound as above costs only polynomial factors
in L. For a prime divisor ℓ>y of any listed scalar, the probability
of an ordinary order divisible by ℓ is O(L²(1/y+1/√X)).
Multiplying by exp(KW) still gives exp(−Ω_K(W)). Outside these two
exceptions and trace zero, the order has a prime divisor >y missing
from the whole menu. One point can be killed by a listed scalar only
if its corresponding primary component is zero, of probability ≤1/y
simultaneously for all paths. The sampler comparison and union over
the two primes preserve the bound. For b-bit balanced inputs W is
Θ(sqrt(b log b)). This counts steering bits only; gcds of arbitrary
coordinate expressions, a coordinate-selected initial curve/point,
and unrestricted numerical algorithms are separate models.

### Gcd tests on polynomials in the initially sampled coordinates

Fix any n=pq for distinct odd primes. The point-first sampler starts
with three independent uniform residues (a,x,y); substitute
b=y²−x³−ax before measuring the total degree d of a proposed test F.
Coefficients may depend on n. Suppose its reductions are both
nonzero as formal polynomials. Schwartz–Zippel on the three free
coordinates bounds the chance of a proper gcd by d/p+d/q. A menu
with sum of degrees T has bound T(1/p+1/q) even if coordinates
select the test; union bound before conditioning on the selection.
Conditioning on nonsingularity at both primes changes the bound by
at most 1/[(1−2/p)(1−2/q)] for p,q>3. Rational tests with invertible
denominators use their numerators; denominator failures are counted
as separate tests. For balanced primes and total menu degree
exp(o(log p)) this contribution decays as p^(−1+o(1)).
Unlike the earlier scalar-order estimate this holds for each pair.

There is an important exceptional case. If F vanishes formally at
both primes the gcd is always n. If it vanishes formally at exactly
one, this *is* a factoring construction: the gcd of n and all its
expanded coefficients is precisely that prime. One need not expand
a polynomial-size circuit when d is polynomial in log n. Assuming
d<min(p,q), evaluate on the grid {0,...,d}³, taking the gcd of n
and all grid outputs. A polynomial of individual degrees ≤d over
either field vanishes on that grid iff it is the zero polynomial
(by univariate interpolation in each variable). The grid gcd
therefore gives the exceptional prime deterministically with at
most (d+1)³ circuit evaluations. Small primes ≤d can first be
trial-divided. This is a decoder for a supplied circuit; it does
not construct one with differing polynomial identities.

For ordinary low-degree tests this bound can be added to the earlier
first-event chain estimate along the formal paths. Any test that is
identically zero on both sides must take its deterministic zero
branch there. Compact circuits with very high formal degree, or a
menu too large for the degree sum bound, remain outside this result.

### A randomized extension-ring decoder for a compact identity

Let a division-free integer arithmetic circuit compute F of total
degree at most D. Suppose that at distinct prime divisors p,q of n
it is the zero formal polynomial modulo p and nonzero modulo q.
This promise suffices to split n in randomized polynomial time in
the circuit size, log n and log(D+2), with no factor or field
construction supplied. In particular exponential degree is allowed:
a size-s binary arithmetic circuit on variables and constants has
degree at most 2^s.

Take even k≥16 with 2^(k/2)≥2 max(1,D). Sample a uniform monic
degree-k polynomial H over Z/n and evaluate the circuit on independent
uniform elements of (Z/n)[T]/H. Take the gcd of n and all k output
coefficients in the monomial basis. It contains p. To see why it
often omits q, reduce H modulo q. For any k/2<j≤k, each irreducible
monic polynomial of degree j divides a uniform H with probability
q^(−j). All these events over j and irreducibles are disjoint. There
are at least q^j/(2j) such irreducibles for each j: at least half the
elements of F_(q^j) have full degree, by bounding the union of proper
subfields by the sum of q^i up to floor(j/2), and each orbit has j
elements. Thus the probability of some such factor is at least
Σ 1/(2j) ≥1/4. Project onto the corresponding field for the proof
only. The circuit inputs project to independent uniform field values,
and the chance of a spurious zero is ≤D/q^j≤1/2 by Schwartz–Zippel.
Consequently a proper gcd occurs with probability at least 1/8 per
trial. Monic reduction uses no ring inversions. Repeating until a
verified split has constant expected trial count under the promise.
The proof permits arbitrary multiplicities in n.

A fixture with n=35 and F=((X+1)^5−X^5−1)(X^7−X) illustrates the
distinction from a function identity: every scalar evaluation modulo
35 vanishes. But evaluating at T with T²+1=0 gives 10T−10, whose
coefficient gcd with 35 is 5. The fixture intentionally uses known
prime exponents. The construction of a separating circuit from n
alone remains the central unsolved step.

### Constructing a separator by short finite differences of powers

For positive M use F_r(X)=Σ_{i=0}^r (−1)^(r−i) binom(r,i)(X+i)^M.
For any prime p>r≥1 it is formally zero mod p exactly when the
base-p digit sum s_p(M) is less than r. To prove this, the coefficient
of X^(M−j) is binom(M,j) times the r-th difference of i^j at zero.
By Lucas, the first factor can be nonzero only when every digit of j
is at most the corresponding digit of M. At such j the values i^j
on F_p equal i^s_p(j) (also at zero, using the constant monomial for
j=0). Thus if s_p(M)<r every coefficient vanishes. Otherwise choose
a digit subvector of M of sum exactly r. Its j gives a nonzero Lucas
coefficient and a moment r! nonzero mod p. This proves nonvanishing
as a polynomial, without requiring a nonzero base-field evaluation.

The circuit has O(r log(M+1)) arithmetic gates via repeated squaring.
So a polynomial-bit M computable from n whose digit sums at two
hidden primes straddle a small r yields a constructible separator.
One can scan polynomially many r with bounded decoder trials, accepting
only verified factors, rather than waiting indefinitely on a bad r.

For n=pq with p<q<2p and M=n the sums are 1+q−p at p and p at q.
Hence r=q−p+2 works when r<p. At n=323, r=4, the values F_r(0),
F_r(1), F_r(2) modulo n are 85,289,102; each gcd is 17. For this
specific M a polynomial-size explicit difference on well-separated
balanced primes cannot produce a one-prime identity: r below both
digit sums gives nonzero formal polynomials on both sides. When the
gap is polynomial in log n, Fermat factoring is already efficient.
A different exponent construction or more compact difference circuit
would be needed to turn this sufficient criterion into an all-n method.

### Small digit sums already yield a direct multiplicative split

For M≥r≥2 form E=∏_{s=1}^{r−1}(M−s). If any odd p|n has s_p(M)<r,
then p−1 divides E since M≡s_p(M) modulo p−1. This only requires
one small digit sum, irrespective of all other digit sums. The bound
M≥r loses no separating finite difference: below that degree the
r-th difference vanishes over the integers. E has O(r log(M+1)) bits.

Write E=2^v u with u odd. For a uniform unit a test gcd(n,a^(2^j u)−1)
for 0≤j≤v via one initial modular power and repeated squaring.
At j=v₂(p−1)−1 the exponent contains the full odd part of p−1
and exactly half its two-part. The event of getting 1 mod p has
probability exactly 1/2 by cyclicity. For any other prime q|n, its
corresponding event is independent, so with probability exactly 1/2
the two events differ. Then the gcd is proper regardless of prime
multiplicities. No knowledge of p or the successful index is needed.
At n=M=323, r=4 this gives E=33075840; exactly 144 of 288 units split
at exponent 2067240. Thus a successful short digit-sum exponent search
would already solve the splitting problem without polynomial identities.

For balanced semiprimes, if M=F(n) comes from the previously bounded
polynomial lists and r is polynomial in bit length, apply the
polynomial-exponent obstruction to the expanded menu F(n)−s. Provided
M≥r, none of the relevant evaluations is zero; discard identically
zero candidates. For the expanded theorem its coefficient-content
budget must still be polynomial, as it is for polynomial-size dense
integer coefficient encodings. A small digit sum would give gcd ratio
one at p−1, so these lists do not give the promised small digit sum
on a nonnegligible fraction of random balanced pairs. This restricts
the choice of M, not arbitrary high-order difference circuits.

### Binomial prefixes expose the least prime, with a truncation cost

Let G_k=gcd(n,binom(n,1),...,binom(n,k)) for 1≤k≤n. Exactly

    G_k = n / gcd(n,lcm(1,...,k)).

For p^e exactly dividing n, j binom(n,j)=n binom(n−1,j−1) gives
v_p(binom(n,j))≥e−v_p(j). Set a=min(e,floor(log_p k)).
At j=p^a the second binomial is a unit mod p by Lucas: the bottom
a digits of n−1 are all p−1. Thus the prefix minimum valuation,
including n itself, is e−a. This proves the identity. Over F_p the
order of vanishing at zero of (1+X)^n−1 is instead p^e by Frobenius.
These two thresholds differ on repeated prime powers: G_k first
loses a factor p at k=p and loses the last one at k=p^e.

In particular G_k=n is exactly the predicate that the least prime
of n exceeds k. A binary-index prefix oracle giving only this bit
would permit deterministic factoring by binary search and recursion.
Ordinary powering modulo (n,X^(k+1)) needs a degree-k representation.
Moreover a nilpotent element z of a finite-dimensional algebra over
a field has z^h=0 for h at most the dimension (use its minimal
polynomial). Hence on squarefree n=pq with p<q, any substitution
into algebras of dimension ≤p at both primes using nilpotent inputs
makes (1+z)^n−1 vanish at both. Random field extensions have no
nonzero nilpotents. This rules out direct low-dimensional nilpotent
substitution as a shortcut to the truncation; it is not a general
lower bound on computed prefix predicates. Explicit coefficients
verified the gcd identity for all 45,149 pairs 2≤n≤300, 1≤k≤n.

### A single-coefficient decision oracle replaces a prefix oracle

For composite n let p be its least prime. For 0≤k<2p with k<n set
B_k=(−1)^k binom(n−1,k). The product formula B_k=∏_{j=1}^k(1−n/j)
can be reduced modulo n: denominators coprime to n are invertible,
and every other j in this range is itself a prime divisor s of n.
For j=s the factor is an integer. Thus

    B_k ≡ 1 − Σ_{s|n prime, s≤k} n/s   (mod n).

All cross terms vanish because distinct selected primes have product
dividing n. For any selected s^e exactly dividing n, the term n/s
has valuation e−1 at s, whereas every other summand is divisible by
s^e. Hence B_k=1 modulo n precisely when k<p throughout this
interval, including when n has repeated prime factors.

Starting at 1, double k until the equality test fails. The first
doubling reaching p lies below 2p and below n, so it cannot miss.
Binary search within the final bracket recovers p with O(log p)
single-coefficient bit queries. Primality testing and recursion give
a deterministic reduction for general n. Nothing about monotonicity
at k≥2p is required. At n=323 and k=17, the signed coefficient
equals 1−19 mod 323; the single deviation has gcd 19 with n.
Using direct binomial coefficients to answer queries, exhaustive
decoding returned the least prime for all 4330 composites ≤5000.

By comparison a row-n coefficient with gcd(k,n)=1 always vanishes
mod n, by k binom(n,k)=n binom(n−1,k−1). The shift to row n−1 puts
the signal at indices that need not themselves share a divisor with
n. Computing the decision bit in poly(log n) at binary indices
remains open here; forming the ordinary coefficient is too large
in balanced cases.

### A fixed central coefficient splits every balanced odd semiprime of ratio below two

For n=pq, distinct odd primes p<q, row n−1 in base q has digits
(p−1,q−1), high first. Lucas's theorem makes every row coefficient a
unit modulo q. In base p the low digit of n−1 is p−1 and the upper
digits represent q−1. For uniform k<n the quotient floor(k/p) is
uniform below q, so the chance gcd(n,binom(n−1,k))=p is precisely
1−∏(a_i+1)/q for base-p digits a_i of q−1. All other gcds equal 1.
When q<2p this reduces to (2p−q)/q; arbitrary errors on a δ fraction
of uniformly queried residues lower this success guarantee by at most δ.

More strongly, for q<2p use just k=(n−1)/2. Its quotient on division
by p is (q−1)/2, strictly below p and strictly above q−p−1. But q−1
has low base-p digit q−p−1. Lucas therefore forces a zero modulo p
at that central index. The evaluator need only supply the residue of
one central coefficient, and the gcd deterministically returns p.
Equivalently the p-adic carry at the central index occurs exactly
when some base-p digit of q−1 is odd: absence of carries means
q−1 is twice (q−1)/2 digitwise. The appended low digit p−1 is even.

Direct checks covered all prime pairs 3≤p<100, p<q<2p. For n=323
and 527 the random-index split counts were respectively 255/323 and
51/527, but the central index split both. For n=799=17·47 the central
index missed, in agreement with q−1=(2,12) in base 17. Computing the
central residue in poly(log n) is the unresolved evaluation task;
the one-index conclusion alone gives no fast way to do so.

### Variable-index central coefficient unit test is factoring-equivalent

Allow the binary index k to vary independently of odd composite n.
Only the bit U(n,k)=[gcd(n,binom(2k,k))=1] is needed. For p the least
prime divisor and 0≤k<p, every denominator factor in (2k)!/(k!)² is
a unit, so U is true exactly for 2k<p. Double k starting from 1;
the first power of two at least (p+1)/2 is still ≤p−1 (since the
preceding power is ≤(p−1)/2). Hence the first false query brackets
the change point, without relying on monotonicity after k reaches p.
Bisect that bracket to recover its first false index h=(p+1)/2
and output 2h−1. This uses O(log p) decisions and k<p≤sqrt(n).
Strip powers of two, test primality and recurse for general inputs.

Conversely the prime factors suffice to compute U for *any* binary k
by checking for a carry in k+k in each prime base (Kummer); this is
polynomial in the two input lengths, including when primes repeat.
Thus the decision problem itself is deterministically factoring-
equivalent. This converse does not assert that full coefficient
residues are easy from a factorization. Direct integer coefficients
verified the decoder on all 7738 odd composites through 20000.
For the odd composite 799, where the fixed central query failed,
the first false variable-index query is k=9, recovering 17.

### What the RH ramp would have to evaluate for factoring

Section I.1 of the Alpöge–Claude zeta exposition uses the ramp-weighted
von Mangoldt sum S(x)=Σ Λ(m)(x/m−1)_+. Insert the principal Dirichlet
character modulo N in this sum to get S₀(x;N). Subtraction isolates
an exact factor signal:

    S(x)−S₀(x;N) = Σ_{p|N} Σ_{a≥1} log(p) (x/p^a−1)_+.

If P is the least prime of N, the difference vanishes through x=P.
At x=P+1 it equals log(P)/P, since every other prime power with
positive weight would have to be at most P. Positivity for every
larger x gives a monotone decision bit and binary search for P.
This identity is unconditional. The exposition's bound on the global
remainder on long intervals supplies neither the restricted sum nor
a certified sign for this subtraction. Even two independent two-sided
absolute error bounds would need combined error below half the first
signal for simple threshold decoding there. This is a concrete
precision target for that subtraction, not a lower bound on other
algorithms or on estimates whose errors cancel.
Source: https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf

### A single macroscopic ramp query on balanced semiprimes

For n=pq, 3≤p<q<2p, choose x=n/2. Now q<x<p², so precisely the
first powers of p and q contribute. The ramp difference is

    R=(q log p + p log q)/2 − log n.

This determines s=p+q stably even near p=q. Put z=√n, L=log z,
p=z exp(−u), q=z exp(u). Across candidate sums one has
s=2z cosh u and 0≤u≤(log 2)/2. The smooth interpolation of R is
z(L cosh u−u sinh u)−log n. Its derivative with respect to s is
(L−1−u coth u)/2, including the continuous endpoint limit.
Since u coth u ≤ cosh u <2 on this interval, consecutive integer
sums give values separated by more than (L−3)/2. In particular
for L≥8 a single approximation of absolute error ≤1 recovers the
exact integer sum by monotone search and nearest-value decoding,
leaving constant slack for numerical error. Search integer sums
between ceil(2√n) and floor(3√(n/2)), computing the two quadratic
roots and logs at polynomial bit precision. Check the discriminant
and multiply the factors. Smaller n can be covered by finite trial
division. At n=1009·2017 both injected errors −1 and +1 decoded
the sum 3026 numerically (the sufficient size cutoff is conservative).
The ramp input was generated from known factors for this decoder test.

This changes the required precision compared with detecting onset.
It requires unit absolute accuracy on the *difference* at a specified
macroscopic scale; the global analytic remainder estimate does not
give that accuracy. The clean two-term expression assumes q<2p.

### A unit-error ramp query suffices on every distinct-prime semiprime

Choose x=n^4. For n=pq define A=log(p)/(p−1)+log(q)/(q−1).
For each prime l, summing the geometric series past its last power
l^h≤x shows that xA−R(x) has local contribution
h log l + x log(l)/(l^h(l−1)). This is positive and at most
log x+2 log l. Thus one raw answer with error ≤1 estimates A with
error ≤(2 log x+2 log n+1)/x. No assumptions about q/p are needed.

For candidate integer k in [2,√n] compute the same slope using real
arguments k,n/k. It decreases strictly: set g(t)=t/(exp(t)−1),
L=(log n)/2 and u=L−log k. The slope is g(L−u)+g(L+u).
On log 2≤t≤log n, g''(t)≥c/n for an absolute c>0. Indeed
g''=e^t[(t−2)e^t+t+2]/(e^t−1)^3 is positive since
t coth(t/2)>2. For t≥3 it is at least exp(−t), and the remaining
compact interval has a positive minimum. The symmetric sum has
derivative zero at u=0 and second derivative ≥2c/n. Adjacent k's
are spaced by at least 1/√n in u, so adjacent candidate slopes
differ by at least c/n² even at the most crowded end. This dominates
O(log n/n^4) for sufficiently large n. Binary search for the nearest
slope recovers p; rounding errors can be bounded with O(log n)
bits and remaining candidates checked by division. Small inputs are
finite. This is still a conditional decoder for the raw ramp difference,
not a way to compute that difference. Tests using 90 decimal digits
and raw errors ±1 recovered p on 3·1000003, 101·1000003, 1009·1013.

### Twisted ramps separate arbitrary prime supports at a public scale

Remove the factor 2. For odd n and a unit d, multiply each summand
of the principal-character ramp difference by the Jacobi symbol (d/m)
at its odd summation index m. Equivalently its value at x is

    R_d(x) = Σ_{p|n} log(p) Σ_{a≥1} χ_p(d)^a (x/p^a−1)_+.

This is an extra evaluator parameter; d=1 is the untwisted case.
At every x≥n the slope R_d(x)/x differs from its infinite-series
limit by at most [ω(n)log x+2 log n]/x, by absolute values in the
geometric tail and the finite constant terms. The limit is

    Σ_{p|n} log(p) [1+pχ_p(d)]/(p²−1).

Thus a prime's coefficient is w(p)=p log(p)/(p²−1). For real v≥3
the negative derivative of w is
[(v²+1)log v−(v²−1)]/(v²−1)² ≥ (log 3−1)/v².
Consequently the candidate integer weights from 3 to n have gaps
at least (log 3−1)/n². Unit absolute error in R_d(n^4) is ample,
including the finite-tail error, for sufficiently large n.

As in the binary-cube density decoder, choose independent uniform
units g_j and evaluate d at their subset products. With b=bitlength(n)
use a power-of-two table of size between 4(b+1)² and 8(b+1)².
Legendre sign vectors for different primes are independent uniform
by CRT. With probability at least 3/4 all are nonzero and distinct.
The normalized Walsh transform then yields one w(p) per occupied
nonzero frequency, to the same uniform absolute error as the input
slopes. Nearest-weight search, primality and divisibility checks,
and repeated division recover the complete factorization on that
event. Repeat if the verified prime powers leave a cofactor. The
logarithm precision required is polynomial in b; bounded small n
can be handled directly.

More generally any computable estimates with raw error O(x^(1−δ))
uniformly over these twists for a fixed δ>0 suffice at x=n^K with
fixed K>2/δ and K>2. This is because the error after division by x
is O(n^(−Kδ)). A theorem bounding a remainder around a main term
that itself contains the unknown slope does not supply such estimates.
The evaluator must actually return the twisted difference or a
computable approximation to it. Given the factors the direct finite
prime-power sums already do so in polynomial time at these scales.

### One Gaussian ring class number supplies all splitting exponents

Specialize the finite Euler product to the fixed character χ=χ_{−4}.
For odd n>1 the standard ring class-number formula in Q(i) gives

    E=2 h(−4n²)=∏_{p^e || n} p^(e−1)(p−χ(p)).

Here h counts primitive proper form classes. The maximal order has
class number one, and its unit index over the conductor-n order is
two. The formula makes E a simultaneous multiple of each p−χ(p).
Only one class-number call is required for the entire factorization.

For a current divisor m of n, sample a uniform z=a+bi modulo m.
If gcd(a²+b²,m) is proper, use it. If it equals m, retry. For
invertible norm form w=z/conjugate(z). At every p|m the local w
is uniform in a cyclic group of order p−χ(p). In the split algebra
the map on unit pairs is (s,t)↦(s/t,t/s); in the inert field it
is z↦z^(1−p), with image of size p+1. CRT gives independence.
Write E=2^v u with u odd and test gcd(m,real(w^j−1),imag(w^j−1))
at j=u,2u,...,E. For a chosen p the index j=u·2^(v₂(p−χ(p))−1)
makes the identity event have probability 1/2. At any distinct q
the event is independent, so exactly one occurs with probability
1/2. This gives a verified proper gcd regardless of multiplicities.
Conditioning on unit norm yields uniform unit reductions; at any odd
prime the ill-norm probability is at most 2/p, so even unconditional
sampling has a constant probability of avoiding the retry case
(gcd=m) and either splitting immediately or entering this test.
The same E works down the recursion tree. Prime powers and powers
of two are resolved separately. E has polynomial bit length since
each local factor is at most 2p^e.

Independent enumeration of primitive reduced positive forms at all
100 odd conductors 3 through 201 agreed with the product. This is
a consequence of the classical formula, not a fast algorithm for
evaluating the class number. The finite Euler product at s=1 is a
different observable from its logarithmic derivative in the ramp.

### Exact Gauss reduction for the norm-one count

For every odd n>1, 2h(−4n²) equals the number C(n) of ordered
residue pairs with x²+y²≡1 mod n. There is a direct character-sum
derivation. Write e_n(z)=exp(2πiz/n), G(t,n)=Σ_x e_n(tx²). Then

    C(n) = (1/n) Σ_t e_n(−t) G(t,n)²
         = Σ_{m|n} (n/m) χ_{−4}(m) μ(m).

Group t by its additive-character order m: t=(n/m)u for a unit u
mod m. The Gauss sum scales by n/m and its primitive square is
χ_{−4}(m)m independently of u, by the odd-modulus quadratic Gauss
evaluation. The remaining sum of e_m(−u) is μ(m), including m=1
with contribution n. Multiplication of local factors recovers
the ring class-number formula. Replacing 1 by any unit right-hand
side gives the identical result since multiplication permutes units
in each Ramanujan sum.

Thus all phase dependence within each order has an explicit sum.
For a semiprime n=pq the four contributions are n, −χ(p)q, −χ(q)p,
and χ(n). In particular ignoring the two proper-order contributions
gives an exactly quantified error; they are the only hidden terms
once squarefreeness and two primes are promised. Grouping by order
still requires locating divisors or evaluating the aggregate by some
other means. The equivalent sum Σ_{m=1}^n [m|n](n/m)χ(m)μ(m) is not
a polylog evaluation procedure. Direct square-residue histograms
verified the identity for all 150 odd moduli from 3 through 301.

### What one threshold bit of a twisted conic count reveals

For odd semiprime n=pq with p<q and unit d, let C_d count solutions
of x²−d y²=1 modulo n. The split and nonsplit norm fibers give
C_d=(p−a)(q−b), where a=(d|p), b=(d|q). Hence C_d−n=−aq−bp+ab.
Since q≥p+2, the last two terms have magnitude at most p+1<q:
C_d<n holds exactly when a=1. On the Jacobi-positive domain,
this single bit is precisely the quadratic residuosity decision
modulo n. This requires varying d; a fixed Gaussian conic supplies
only the character value at −1.

Conversely fix a Jacobi-negative reference u, found by sampling units.
On Jacobi-positive inputs the residuosity oracle gives a directly.
On a Jacobi-negative d, du is Jacobi-positive and a(d)=a(du)a(u).
Thus one orientation bit a(u) plus the ordinary residuosity oracle
simulates the full threshold oracle; the opposite guess implements
the character at q instead. For a proposed factoring reduction with
verifiable outputs, run both simulations and verify any returned
factor, eliminating the advice. Consequently a reduction from
factoring these semiprimes to this threshold oracle would also be a
randomized reduction to quadratic residuosity. No such reduction is
proved here. Direct square histograms checked the threshold on every
unit for 15,35,77,143,323.

### Deterministic polynomial interpolation from a fixed split conic

Let P(M)=#{(x,y) mod M: xy=0 mod M}=Σ_x gcd(x,M), the classical
Pillai gcd-sum function. The count is multiplicative by CRT. At a
prime power p^a the mean P(p^a)/p^a equals 1+a(1−1/p): the zero
residue contributes 1 and each exact valuation 0 through a−1
contributes 1−1/p. Consequently at M=n^k the mean is the polynomial

    F(k)=∏_{p^e||n} (1+k e(1−1/p)).

For b=bitlength(n), remove primes ≤b by trial division and apply
the formula to the remaining cofactor m. With k=0,...,b, interpolate
m F(k), an integer-coefficient polynomial since rad(m) divides m.
Errors strictly below 1/(2m) in each mean can be removed by rounding
m times the estimate. Each root is −p/[e(p−1)], already reduced
because p>b≥e; hence its numerator in absolute value is the hidden
prime. Standard deterministic rational polynomial factorization and
repeated division recover the complete decomposition. The coefficients
have polynomial bit length: their absolute sum is bounded by
m(1+b)^b. Queries have O(b²)-bit moduli. Factoring conversely gives
P exactly by the product formula. This is a deterministic reduction
using a single equation and an ordinary modulus as input, with no
choice of twist or unknown field.

Direct summation and interpolation gave (4k+5)(4k+3)/15 at n=45
and (6k+7)(4k+5)(2k+3)/105 at n=105. On squarefree semiprimes
already P(n)=(2p−1)(2q−1)=4n−2(p+q)+1, so a single count suffices.
Direct evaluation still takes a sum over M residues; the easy local
gcd term does not compute the sum in polynomial bit complexity.

### Sampling the split conic and a slow conditional-update chain

Write P(n) for the zero-product count. Under the uniform solution
law, testing gcds of both coordinates splits with exact probability

    1 − (2φ(n)+1)/P(n).

The only failures are (0,0), (unit,0) and (0,unit). For odd composites
the success probability is at least 1/4: if n has ≥2 distinct primes,
P(n)/n≥(1+2/3)² and 2φ(n)+1≤2n−1, giving success >1/4.
For n=p^e, e≥2, the x-coordinate alone has an intermediate valuation
with probability (e−1)(1−1/p)/(1+e(1−1/p))≥2/7. Hence any polynomial
sampler within total variation ≤1/8 gives a constant-probability
verified split; independent retries and recursion suffice. Conversely
known factors permit exact uniform sampling locally: choose the
x-valuation using weights 1−1/p for 0,...,e−1 and 1 for e, choose
uniform x of that valuation, then uniform y in its annihilator, and
combine by CRT. Direct enumeration checked the failure counts for
9,15,35,77,323.

The simple Gibbs chain, which refreshes one coordinate uniformly given
the other, has the uniform solution law as a stationary distribution.
Nevertheless, on n=pq its escape probability per update from the
failure set is at most (p+q−2)/n: fixing a unit forces the other
coordinate to zero; fixing zero gives a uniform residue, only p+q−2
of which expose a proper gcd. Thus after t updates from (1,0), the
chain's split probability is ≤t(p+q−2)/n, whereas at stationarity
it tends to 1/2 on growing balanced semiprimes. The same escape
bound applies from stationarity conditioned on the failure set, a
warm start of asymptotic density ratio two. Reaching TV ≤1/8 on
balanced semiprimes requires Ω(√n) such updates. This gives a
concrete failure mode for an otherwise cheap exact conditional sampler.

### Unconditioned multiplicative convolution at a binary-encoded time

Let Y_k be the product modulo n of k independent uniform residues.
For any prime p|n the event p∤Y_k has probability s_p=(1−1/p)^k;
these events are independent across distinct primes by CRT. For odd
n with at least two distinct primes, take p its least prime and k
the first power of two with s_p≤1/2. It exists by k<2p and the
preceding survival exceeds 1/2, hence 1/4<s_p≤1/2. For any other
prime q, exactly one of p,q divides the product with probability
s_p+(1−2s_p)s_q≥s_p>1/4. This guarantees a proper gcd irrespective
of the other primes or their multiplicities.

A sampler taking n and binary k and returning a residue within 1/8
total variation of Y_k therefore suffices: try dyadic lengths up
through the first power of two ≥n, verify outputs by gcd, and repeat
with fresh randomness. There are O(log n) calls per pass and no need
to locate the crossing in advance. Handle prime powers and powers
of two separately. Exact small-modulus convolution confirmed the
split bound at n=15,35,45,105,323.

The target is unconditional product sampling, but the binary length
of k matters. The explicit product takes k draws; repeated squaring
of one sampled value never introduces a new prime divisor of n.
More generally products and positive powers of t sampled residues
can expose a proper gcd on balanced pq only if an input was divisible
by p or q, an event of probability ≤t(1/p+1/q) for uniform draws.
This explains why a multiplication-only reuse of polynomially many
samples cannot implement the required distribution at the useful
length. It does not rule out additive operations or other samplers.

### What the additive Fourier trace of a long product preserves

Write Y_j for the product of j independent uniform residues mod n,
Y_0=1. Conditional averaging over the last residue gives exactly

    E exp(2πitY_k/n) = Pr[n/gcd(t,n) divides Y_(k−1)].

All unit frequencies therefore measure just z=Pr[Y_(k−1)=0].
At k=3 this scalar is the earlier Pillai count divided by n².
In particular the character trace at fixed t=1,k=3 would already
suffice for deterministic factoring via the previous power-modulus
interpolation, given enough precision to round the required counts.

There is a distinct distribution with exactly the same unit-frequency
coefficients: return zero with probability z and an independent uniform
residue otherwise. Its proper-gcd probability for n=pq is only
(1−z)(p+q−2)/n. At the useful dyadic k the actual Y_k has split
probability >1/4, so on growing balanced semiprimes the two laws
differ in total variation by at least 1/4−o(1). For n=323,k=16,
z=(1−(16/17)^15)(1−(18/19)^15); actual and surrogate split rates
are respectively about 0.481 and 0.070. A nonzero nonunit frequency
can itself be checked for a proper gcd before evaluating the trace.

This concerns what matching the coefficients guarantees about a
sampler. It does not rule out decoding the number z using its known
arithmetic dependence on n, as the k=3 identity illustrates.

### Ten-percent relative cubic estimates recover ω on arbitrary inputs

Let F(M)=P(M)/M, equivalently M times the normalized cubic Fourier
average. For b=bitlength(n) and k=16b the ratio F(n^(2k))/F(n^k)
is within a known factor of 2^ω(n). Indeed for each p^e||n set
a=e(1−1/p)≥1/2. Its contribution divided by 2 equals
1−1/[2(1+ka)]≥1−1/k. Since ω(n)≤b, the full ratio lies between
(15/16)2^ω and 2^ω. If both positive values have relative error
at most ten percent, their observed ratio divided by 2^ω is between
(15/16)(9/11) and 11/9, strictly inside (1/√2,√2). Round the base-2
logarithm to obtain exactly the number of distinct primes, without
a squarefree promise. The powered inputs still have O(b²) bits.
Finite checks using local formulas verified both extreme perturbations
for every integer from 2 through 10000. This recovers a count, not
the factors.

On promised squarefree inputs the reverse reduction holds at the
original modulus. Strip primes through 16b, multiplying their exact
local factors. If the remaining cofactor has t primes, replace its
F by 2^t: its true value divided by 2^t is ∏(1−1/(2p))≥1−1/32,
so the relative overestimate is at most 1/31. To decode t from a
ten-percent oracle using squarefree inputs only, query it on the
stripped cofactor directly and round log₂ of its answer; the known
factors are counted separately. Thus relative approximation on
squarefree moduli alone and computing ω on that promise are
deterministically interreducible.

Precision convention matters: P(M)/M²≤τ(M)/M≤2/√M by divisor
pairing. Returning zero for the ordinary cubic average therefore
already meets any fixed inverse-polynomial *additive* accuracy in
bitlength(M) for all sufficiently large M. Such accuracy gives no
relative guarantee for the two powered queries.

### Ten-percent relative twisted cubic estimates count split primes

For odd n and gcd(d,n)=1 consider the normalized cubic trace of
z(x²−dy²) on three uniform residues mod M=n^k. Averaging over z
makes it the probability of x²=dy² mod M. Multiply this trace by M
to obtain the solution count divided by M, denoted F_d(k). At an
odd prime p^e||n where d is a square, a square root lifts to p^(ke);
the invertible change of variables (x−√d y,x+√d y) gives the same
local factor 1+ke(1−1/p) as xy=0. At a nonsquare prime the valuation
of x²−dy² is twice min(v_p(x),v_p(y)) until it reaches ke: after
removing the common power of p the residue cannot vanish. For even
k this forces both coordinates to be multiples of p^(ke/2), giving
exactly p^(ke) solutions and normalized local factor 1.

Thus at even powers only the square primes occur in the product.
Apply the previous ratio argument with k=16 bitlength(n): ten-percent
relative estimates of F_d(2k) and F_d(k) recover the exact number
t_d of distinct p|n where d is a square. For d=1 this gives ω(n);
equality t_d=ω(n) decides quadratic residuosity modulo n for unit d
by lifting and CRT, including repeated factors. The d=1 denominator
and numerator can be reused across all d. These are single cubic
phase evaluations with an explicit parameter, a stronger evaluator
than access just to the fixed xyz phase. Interpolating the full values
would reveal the primes, but the ratio test alone only delivers their
count in each sign class. Direct square histograms checked the
local-product prediction at M=n² for all units 1≤d<n and odd n≤35.

### Converse for coarse twisted estimates on semiprimes

Promise n=pq for distinct odd primes, and use even k. QR decision
on Jacobi-positive units suffices to approximate F_d(k) uniformly
within 10 percent: a negative Jacobi symbol means exactly one split
prime; a positive symbol means two or zero according to QR. If t
is this count, return (1+k)^t after first trying trial divisors
through 23 (if a divisor is found compute the local formula exactly).
For every remaining prime the ratio of the true split factor to 1+k
is 1−k/[(1+k)p]≥1−1/p. Thus the relative overestimate is at most
(1−1/23)^(−2)−1<0.1 independently of k. For inverse-polynomial
relative tolerance ε≤1, trial division through ceil(8/ε) gives the
same simulation within ε. Evaluator input here is the known base n
and even exponent k; the simulated queries are on powers of n.

Conversely a single k=2 query decides QR on a Jacobi-positive unit.
Two nonsquare signs give F_d(2)=1; two square signs give a value
at least (3−2/3)^2. Threshold 2 distinguishes even with ten-percent
relative error. Thus on this promise, producing estimates to this
guarantee and deciding QR are deterministically interreducible.
This statement concerns algorithms that must work with any estimate
within tolerance, not hidden information in a particular implementation's
extra digits. It does not apply to odd-power relative estimates.
The local formulas at n=29·31 passed the simulation error check
for all units and k=2,4,16,128, including both endpoints of the
allowed estimation error for the single-query decision.

### Odd-power relative estimates and guessing a factor magnitude

For distinct odd primes n=pq and an odd power k, the nonsquare local
factor is 1/p (both coordinates must be divisible by p^((k+1)/2)).
In the mixed-sign case let s denote the split prime. The entire
normalized answer has the particularly simple form

    F_d(k) = ((1+k)s−k)/n.

It measures s to roughly the same relative accuracy as the oracle,
independently of the size of k. For two nonsquares F=1/n; for two
squares use the earlier product. In particular increasing k does
not amplify the relative precision of a mixed-sign query.

For a robust polynomial-time *factoring decoder* restricted to
queries on powers of its given semiprime, odd-power queries at
relative error ε do not evade the reduction to a QR decision oracle
when 1/ε is polynomial in the bit length. To see this, trial-divide
through ceil(8/ε) as before and take 0<ε≤1/2. Choose a random unit
u with negative Jacobi symbol and try both signs for its character
at the smaller prime. QR on Jacobi-positive units then simulates
both labeled local characters for all d: for negative Jacobi d use
du and the guessed sign. Nonunit d can be gcd-checked.

Enumerate the grid 2(1+ε/8)^j up to √n, assigning the smaller prime
the grid estimate h and the larger n/h. One grid point is within
a multiplicative factor 1+ε/8 of the smaller prime, and both assigned
sizes then satisfy that same bound. There are O(log n/ε) guesses.
For mixed odd queries plug the appropriate guessed size into the
linear formula. Its relative error is at most (ε/8)s/(s−1)≤ε/4.
For even queries or two positive signs use the previous simulator
returning (1+k)^t. For two negative signs and odd k return 1/n.
The correct grid and orientation jointly give valid answers on the
whole adaptive path. Run the time-bounded decoder for each guess,
retaining only verified factors. Independent randomness gives the
promised success probability on the correct path. This is a
search reduction with verification; it does not yield a standalone
odd-power estimator from QR with a certified choice of guess.

As a precision benchmark on balanced p<q<2p, ε≤n^(−1/4−δ) for
fixed 0<δ<1/4 instead gives a factor via a Jacobi-negative query:
recover an approximation to s from the linear formula, with error
O(n^(1/4−δ)); with probability 1/2 the split prime is the larger.
The usual unknown-divisor small-root bound then applies. The regime
between inverse-polynomial relative precision in log n and this
sufficient precision remains open here.

### A tunable reference phase converts the gcd sum into a sign query

Let M>1, A=P(M)/M² be the xyz trace and Q=16M². With a fair Boolean
selector, use the xyz phase exp(2πi xyz/M) on one branch and a
constant phase exp(±2πit/Q) on the other, with independent fair sign.
The resulting real average is (A+cos(2πt/Q))/2. All phase angles
can use the single known denominator Q since M divides Q. This is
a gated phase family with two Boolean inputs, not the previous
unmodified cubic evaluator.

At t=Q/4 this average is positive; at Q/2 it is negative because
A<1. Bisect on the sign to adjacent indices enclosing its zero,
accepting equality directly. Adjacent cosine values differ by at
most 2π/Q, hence either endpoint yields P(M) by rounding minus its
cosine times M²: the error is at most π/8<1/2. There is constant
rounding slack for certified polynomial-bit cosine arithmetic.
O(log M) sign queries suffice. Any signed relative-error guarantee
strictly below one on the combined value gives this sign. The
earlier deterministic interpolation then factors arbitrary N using
M of polynomial bit length. This reduction needs no relative oracle
for odd-power twists and uses no Coppersmith step.

The precision requirement has moved into the sign query: near
cancellation, separately estimated branches with ordinary relative
errors do not control the sign of their sum. Testing the bisection
with exact gcd sums and high-precision reference phases recovered
P(M) at M=9,35,45,105,323,1022117; this exercises the decoding
and supplies no implementation of the combined sign oracle.

### Sign margins and the uniform-reference baseline

Write P=P(M), A=P/M² and use the combined gate of the preceding
section with Q=64M². To test P>h for an integer 0≤h<M², target
cos(2πt/Q)=−(h+1/2)/M² on the decreasing quarter circle. Choose an
integer t within one grid step of the real index (approximate the
inverse cosine index to within 1/4 and round). The cosine discrepancy
is ≤2π/Q. Thus the gate sign gives the exact comparison and its
absolute mean is at least (1/2−π/32)/(2M²)>1/(5M²). Ordinary
integer binary search recovers P with O(log M) promised-gap queries.
The inverse cosine and grid computation use polynomially many bits,
including at the endpoint approached to within 1/(2M²). Direct
tests checked every search path for all 499 moduli 2 through 500.

For a uniform t over the full grid, predicting the sign of the
reference cosine alone errs only on the arcs −A≤cos≤0, with tie
conventions affecting O(1/Q). Their total fraction is arcsin(A)/π
up to O(1/Q) grid discrepancy. Since A≤τ(M)/M≤2/√M, the failure
rate is O(M^(−1/2)), exponentially small in bit length. This
baseline uses no cubic evaluation and its sign is determined by
integer quadrant comparisons. Thus very high uniform-reference
accuracy by itself does not certify the adaptive comparisons above.
This concerns the gated reformulation and does not obstruct different
query distributions or the original average-case window reduction.

### Cheap power traces versus the shear's invariant-space dimension

The permutation U(a,x)=(a,x+a) on pairs modulo n has exactly P(n)
cycles, because translation by a has gcd(a,n) orbits. Equivalently
P(n) is the multiplicity of eigenvalue one of the permutation matrix.
For any binary k its power maps (a,x) to (a,x+ka); its trace equals
n gcd(k,n). Hence every power trace is exactly evaluable in polynomial
bit complexity. The projector onto the invariant space is the mean
of the first n powers, whose trace recovers P(n) but requires the
same sum as before. Explicit permutation orbits verified the count
at every modulus 2 through 64.

Restricted filter obstruction: the block a=1 already contains every
n-th root of unity. A trigonometric polynomial T of degree d whose
values on this grid approximate the indicator of 1 uniformly within
1/4 must have d=Ω(n). Indeed Bernstein gives sup|T'|≤d sup|T|
in the angular variable. If d<n/(2π), the nearest grid point
bounds sup|T| by 5/4+(πd/n)sup|T|, hence by 5/2. Values at angles
0 and 2π/n differ by at least 1/2, implying d≥n/(10π). Larger d
already satisfies the claim. This excludes direct low-degree uniform
projector approximation, including for a constant operator-norm
tolerance. Degree alone does not bound circuit complexity: modular
powers at binary indices are cheap here. Nor does this address
filters targeting only the trace through cancellation. Those
distinctions matter for any attempted fast summation.

### Exact short projector circuit and a determinant with small residue outputs

For S_j=I+U+...+U^(j−1), form (U^j,S_j) by binary recursion:
S_(2j)=(I+U^j)S_j and S_(2j+1)=S_(2j)+U^(2j). Thus S_n/n is
already an exact projector described by O(log n) matrix operations.
The large degree lower bound does not explain the remaining cost:
taking the trace of this succinct matrix expression is the hard step.

For odd n every cycle length d divides n and is odd. On a d-cycle
det(I+U)=1−(−1)^d=2, by its characteristic polynomial X^d−1.
Consequently det(I+U)=2^P(n). One need not ask a determinant oracle
for an exponentially long output. Given a prime r≥2, reduce modulo
2^r−1; the canonical residue is exactly the single bit 2^(P(n) mod r).
Its bit index gives P(n) modulo r, including r=2. The first 2b
primes for b=bitlength(n) have product exceeding n² and largest
member O(b log(2b)); ordinary CRT then recovers P(n) exactly.
Each determinant residue uses polynomially many bits. Together with
the earlier interpolation at powers of the odd cofactor, this gives
a deterministic reduction for arbitrary input after stripping twos.
Conversely known factors give P directly, hence all these residues.
Direct integer determinants of the n individual translation blocks
and CRT decoding were checked for all nine odd n from 3 through 19.

The number of rows is n² and even enumerating the n blocks costs
more than polynomial in b. This isolates a structured modular
determinant target; no generic determinant bound on an explicit
matrix supplies the needed fast evaluation.

### One 3-adic determinant residue per gcd sum

For odd n>1 choose the smallest t with 2·3^(t−1)>n². A single
residue D=det(I+U) mod 3^t recovers all of P(n). The order of 2
modulo 3^j is 2·3^(j−1): odd powers cannot equal 1 mod 3, and
v₃(4^m−1)=1+v₃(m) by the lifting-the-exponent identity. Start the
exponent at 0 or 1 according to D mod 3. Given the exponent e
modulo h=2·3^(j−1), exactly one of e,e+h,e+2h matches D modulo
3^(j+1) after raising 2 to that exponent. Lift through j=t−1.
Since P(n)≤n² the resulting exponent is exact. The modulus has
O(log n) bits and the decoder uses O(log n) small modular powers.
This is an elementary smooth-order discrete logarithm; there is
no generic hard discrete log in the decoding step.

At n=pq for two distinct odd primes use P(n)=4n−2(p+q)+1.
One may even select 2·3^(t−1)>n+1: the exponent residue gives
P(n) modulo 2·3^(t−1), hence p+q modulo 3^(t−1). As p+q
≤3+n/3 < (n+1)/2 for n>15 (with n=15 handled directly), this
already fixes the sum. The general-input interpolation reduction
still queries multiple powered moduli; only the recovery of each
individual gcd sum has become a single call. All 499 odd n from
3 through 999 passed the lifting test on residues generated from
direct gcd sums. No fast determinant computation is assumed in
these tests.

### Partial determinant digits suffice on balanced semiprimes

Fix 0<ε<1/8 and n=pq with distinct odd p<q<2p. A constant list
(depending on ε) of determinant residues, each of bit length at most
(1/4+ε)log₂ n+O_ε(1), suffices for deterministic polynomial-time
recovery. This reduces per-query precision, not necessarily total
output length. Pick J>4/ε fixed distinct odd primes ℓ. Trial gcd
them with n. Set T=ceil((1/4+ε)log_ℓ n), and obtain the determinant
mod ℓ^(T+c), with c=v_ℓ(2^h−1), h=ord_ℓ(2), fixed constants.
Raising the residue to h gives (2^h)^P. Its order modulo ℓ^(c+j)
is ℓ^j by LTE, so ℓ-way digit lifting recovers P modulo ℓ^T.
Thus s=p+q=(4n+1−P)/2 is known modulo ℓ^T.

Write a=v_ℓ(q−p). For at least one listed prime ℓ^a≤n^(ε/4)
since the product of these prime powers divides q−p<n. For this
prime 2a<T. The discriminant s²−4n modulo ℓ^T has valuation 2a.
Divide it by ℓ^(2a) and Hensel-lift its two unit square roots modulo
ℓ^(T−2a); multiply back by ℓ^a. This recovers ±(q−p) modulo
B=ℓ^(T−a) without enumerating the repeated roots. Combining with
s and dividing by 2 yields both factor residues mod B. Skip primes
with zero discriminant residue. All operations use known data.

Here B≥n^(1/4+3ε/4), q<√(2n), so q=r+Bx with |x|≤√(2n)/B
for its candidate residue r. The monic linear polynomial
X+r B^(−1) mod n has this root modulo the unknown divisor q≥√n.
For sufficiently large n the bound is below n^(1/4−ε/2), within
the standard unknown-divisor Coppersmith theorem. Try both roots
at every usable prime and verify by gcd. Bounded n can be checked
directly. This argument assumes access to the specified determinant
digits; their efficient computation from the shear is still open.

### A quadratic shear isolates fixed points in small 3-adic precision

For odd semiprime n=pq, p<q, select a unit d of Jacobi symbol −1.
Let V_d fix x,y and translate z by x²−d y², all modulo n. It has
n(2s−1) fixed points, where s is whichever prime sees d as a square,
by the local split/inert counts. A cycle of length j in det(I+3V_d)
contributes 1−(−3)^j. All cycle lengths divide n. If the least prime
factor is ≥T+1 then every nontrivial cycle contributes 1 modulo
3^(T+1); otherwise trial division through T already splits n.
Handle 3 separately. Consequently the determinant at that modulus
is 4^[n(2s−1)]. Ternary exponent lifting gives n(2s−1) mod B=3^T;
2n is invertible there. The recovered s residue is a unit, so
n/s mod B gives the other factor residue as well, without knowing
the orientation. There is no quadratic discriminant to lift.

For any fixed 0<ε<1/8 take T=ceil((1/4+ε)log_3 n). One residue
now suffices even for unbalanced distinct primes. To spell out the
unknown-divisor bound, let β=log_n q∈[1/2,1]. From its residue r
the root x=(q−r)/B has size ≤n^(β−1/4−ε). Try a fixed grid of
lower bounds β₀≤β with β−β₀≤ε/8, starting at 1/2. Since
β−β²≤1/4 and β²−β₀²≤ε/4, this x is below n^(β₀²−ε/2)
for large n. Apply the monic linear small-root theorem modulo the
unknown divisor q≥n^β₀ to X+r/B (inverse of B taken mod n).
Try both residues and verify factors. Sampling a Jacobi-negative
unit costs constant expected trials or immediately gives a gcd.
The matrix query has changed from the two-coordinate linear shear
to this three-coordinate quadratic shear. No polylog determinant
algorithm follows. At n=35 and modulus 81 direct multiplication
of translation-cycle determinants agreed with the fixed-point
expression for d=2,3.

### Recovering all prime powers with polynomially many determinant queries

Let L be the bit length of N. Strip primes up to T=2L+2 by trial
division, and call the remainder n. For any unit d, let C(d) count
pairs with x²=d y² modulo n. The cycle argument above returns
n C(d) modulo 3^T from one determinant modulo 3^(T+1). Divide by
the unit n there. Since C(d)≤n²<3^T this gives the exact count.

Choose m independent uniform units u_i with H=2^m the smallest
power of two above 4(L+1)^2. Query C at every subset product u^t
for t in the binary m-cube. For each prime power p^a exactly
dividing n let v_p record which u_i are nonsquares modulo p.
These vectors are independent uniform m-bit strings by CRT.
With probability greater than 1/2 they are distinct and nonzero,
by a union bound over collisions and zero vectors.

For nonzero v form the rational product of C(u^t) raised to
the signs (-1)^(v dot t). At p^a the two local counts are
S=p^(a-1)((a+1)p-a) for squares and I=p^(2 floor(a/2)) for
nonsquares. For S, change coordinates to the two linear factors
and count pairs whose product vanishes by valuation. For I, both
coordinates must be divisible by p^ceil(a/2). Character
orthogonality cancels every local factor except those with v_p=v,
each of which contributes (S/I)^(H/2). Thus on the distinct,
nonzero event each occupied slot gives the exact rational R=S/I
by integer roots of numerator and denominator. Empty slots give 1.

Try each a from 1 through L. For odd a the candidate prime is
(R+a)/(a+1), and for even a it is a/(a+1-R). Skip undefined
values and accept only integral primes dividing n. Extract their
full multiplicities by division. Repeat with fresh units until
the verified factors exhaust n. On each round the success event
has constant probability. All intermediate integers have O(HL)
bits, and even the naive product loops take polynomial time in L.
Uniform units can be obtained by rejection after the trial division
(the chance of a nonunit is at most L/T). This reduces general
factoring to exact evaluation of the same succinct determinant on
arbitrary unit parameters with O(L)-bit precision. It leaves the
runtime of that evaluation open.

### Reverse reduction at arbitrary output precision

Given the factorization of odd n and a unit d, the determinant modulo
3^K is deterministic polynomial time in log n and K, including when
small primes divide n. For 1≤j<K put g=gcd(n,j). The jth iterate fixes
exactly F_j=n g² C_{n/g}(d) triples, since its translation vanishes
precisely when x²=d y² modulo n/g. Each pair there lifts g² ways and
z is free. Compute C_{n/g} by the preceding local counts (C_1=1).
In increasing order obtain the number h_j of cycles of length j by
subtracting sum(i h_i for proper divisors i of j) from F_j and dividing
by j. Multiply (1−(−3)^j)^h_j modulo 3^K over j<K. All longer cycles
contribute 1. The counts are at most n³; modular powering, small
divisor enumeration, and the local Legendre tests have polynomial cost.

Together the two reductions give randomized polynomial-time Turing
equivalence between factoring and evaluating this determinant family
for odd n, unit d, and polynomial output precision. For a factor-free
approach it already suffices to compute the single pair count C_n(d)
on the rough odd moduli of the forward reduction; higher iterates
carry no additional information at the precision used there.

### Fixed d=1 suffices deterministically by querying powers of the input

For a semiprime n=pq with distinct odd primes, C_n(1)=4n−2(p+q)+1,
so one exact count gives p+q and the integer quadratic recovers p,q.
More generally, let L be the bit length of the original input. Strip
primes ≤4(L+1)^2 by trial division. For the remaining odd n>1 query
the fixed shear d=1 at moduli n^k for k=1,...,L. All prime divisors
are above the precision cutoff T_k=2 bitlength(n^k)+2 so the earlier
cycle isolation and ternary lifting recover C_{n^k}(1) exactly.

Normalize by n^k. The local split count gives the rational polynomial
P(k)=C_{n^k}(1)/n^k=∏_{p^a || n}(1+k a(p−1)/p).
Use P(0)=1 and rational interpolation; the degree is at most L.
Factor P over the rationals by deterministic polynomial-time polynomial
factorization. For every root r form −1/r=a(p−1)/p. Its reduced
denominator is p because p>L≥a. Verify by division and remove the
whole prime power. This avoids searching for Jacobi-negative units
and works for arbitrary multiplicities. The query moduli have at
most L² bits. Values and interpolated coefficients have polynomial
bit length (the expanded product has common denominator dividing n
and numerator coefficients bounded by n(1+L)^L).
The local formula also identifies the evaluator target with the
classical gcd sum: (u,v)=(x−y,x+y) is invertible at odd moduli,
and uv=0 has exactly gcd(u,n) choices of v for each u.

### Linear output precision for all fixed-shear power queries

Set T to the least integer with 3^T>N² for the original input N;
strip primes ≤max(T,L+1,3), where L is the original bit length.
For n>1 remaining and 1≤k≤L put m=n^k and Q(k)=n P(k).
This is integral because the denominator of P divides rad(n).
It is bounded by n²: locally 1+k a(p−1)/p ≤1+k a ≤(1+k)^a
≤p^a. Thus the least nonnegative residue of Q(k) modulo 3^T
determines it exactly, independently of k.

Query det(I+3V_1) at m modulo 3^(T+1). Ternary lifting gives
m G(m) modulo 3^T, because every nontrivial cycle has length >T.
Multiplying this residue by n/m² modulo 3^T gives Q(k) directly;
m is a unit. Interpolate over the rationals on the recovered
integers, with Q(0)=n. This avoids any precision loss from factorials
in modular interpolation. Each answer takes at most 2 log₂ N+O(1)
bits. Even the number of positive evaluation points can be reduced:
with B=max(T,L+1,3), take D=floor(log_B n) by integer comparisons.
The degree ω(n) is at most D since every remaining prime exceeds B.
Thus D queries suffice and the largest query modulus has O(LD) bits.

### Small-prime traces suffice if the determinant coefficient can vary

For L=bitlength(N), strip primes through L+1 as before, obtaining n.
For each m=n^k used in interpolation, ask for det(I+ℓV_1) mod ℓ² at
successive primes ℓ not dividing n. The linear coefficient of the
determinant is the trace, so the returned residue is
1+ℓ m G(m) mod ℓ². This identity holds regardless of cycle lengths.
Subtract one, divide by ℓ and multiply by n/m² in the field F_ℓ.
The resulting residue is Q(k) mod ℓ. Stop when the product of
the distinct primes exceeds n². Ordinary CRT then gives Q(k)
exactly by the preceding bound. Use the same prime list for all k.

The skipped primes have product at most n. By the Chebyshev bound
sum_{p≤x} log p = Θ(x), primes up to O(L) suffice, so each answer
uses O(log L) bits, with O(L/log L) primes for each k. The oracle
here permits varying the coefficient ℓ; the preceding fixed-3
result does not require that flexibility. Equivalently a uniform
algorithm for G(m) mod ℓ in time polynomial in log m and ℓ over
these inputs suffices for deterministic polynomial-time factoring.
This uses the first-order determinant identity only; an efficient
method for the modular trace itself is still needed.

### A divisibility decision oracle suffices

After the same trial division, write the slopes of P as A_i/p_i
where A_i=a_i(p_i−1). They are nonzero and distinct over Q since
their reduced denominators are distinct. For a prime ℓ coprime to n,
ask just whether ℓ divides G(n^k), for k=1,...,ℓ−1. The positive
answers are precisely the nonzero roots of P modulo ℓ. No root is
zero since P(0)=1. If all slopes stay distinct and nonzero, the
root set R determines the full polynomial modulo ℓ as
∏_{r in R}(1−X/r). Otherwise it has fewer distinct roots.

One can choose a sufficient search range without knowing the
factorization. All exceptional primes, including primes dividing n,
divide the nonzero integer
E=n ∏ A_i ∏_{i<j}|A_i p_j−A_j p_i|.
With L=bitlength(N), at most L slopes, A_i≤L N and each difference
at most 2 L N², this integer has O(L³) bits. Enumerate consecutive
primes until their product exceeds
n (L N)^L (2 L N²)^(L(L−1)/2) n².
The product of the nonexceptional primes then exceeds n². The
Chebyshev bound gives maximum prime O(L³), hence polynomial total
query size and time outside the oracle. Discard primes dividing n.
Keep exactly the primes attaining the maximum observed root count:
this equals the true degree since at least one nonexceptional prime
was tested. For those primes reconstruct P modulo ℓ from R, then
CRT-reconstruct nP(j) for j=0,...,deg P as above. Interpolation
and rational polynomial factorization finish the deterministic reduction.

This changes the query model: each answer is a single bit, while
the number and size of queries can exceed the modular-value version.
For squarefree semiprimes, the positive answer has a direct
interpretation for k≠−1 modulo ℓ: one of the prime factors lies in
the class k/(k+1) modulo ℓ. For k=−1 there is no root in that case.

### Direct cost of a semiprime residue-class decision by progression products

Given n=pq with distinct primes, ℓ coprime to n, and a unit class r,
the smaller prime must be in {r,n/r mod ℓ} whenever either prime
is in r. In both classes take positive representatives ≤floor(sqrt n).
A product over all these integers with gcd 1 against n certifies no.
A nonunit product gives a factor: if a block gcd is n, inspect
the entries of that block to find a nonunit. Each entry is <n.
Then check both recovered primes modulo ℓ to decide the actual bit.
This last check matters: composite multiples also trigger the gcd.

For M representatives per progression, take b=ceil(sqrt M), skipping
empty progressions. If the first term is h, form the polynomial
F(X)=∏_{i=0}^{b−1}(X+iℓ) over Z/nZ. Evaluate at h+jbℓ for all
full blocks by a monic subproduct/remainder tree. A final partial
block costs O(b) multiplications. Using fast polynomial arithmetic
and testing the gcd of each block value costs softly O(b+1) modular
operations (up to polylogarithmic bit costs). This is the usual
baby-step/giant-step factorial technique. M≤1+sqrt(n)/ℓ. Hence
for ℓ polynomial in log n the resulting worst-case bound still
contains n^(1/4)/sqrt(ℓ). This calculation is an upper bound for
one approach and gives no general lower bound for the predicate.

### The query prime 3 and the fixed discriminant −3

For n=pq with distinct primes >3, the oracle bit at ℓ=3,k=1 is
[3 divides G(n)]. It is true exactly when p or q is 2 modulo 3,
since G(n)=(2p−1)(2q−1). Equivalently it is true exactly when
−3 has no square root modulo n: modulo a prime >3, a square root
of −3 is equivalent by the quadratic formula to a root of X²+X+1,
hence to an element of order 3. Apply CRT. If n=2 modulo 3 the
bit is always true. For n=1 modulo 3 the factors have equal
residue classes; both possibilities give Jacobi symbol (−3/n)=+1.

For example 55 and 91 are both 1 modulo 3. The bit is true for 55
(factors 5,11) and false for 91 (19²=−3 modulo 91). All other k
modulo 3 give the constant false answer for these semiprimes:
the local factor (k+1)p−k is nonzero for k=0,−1 modulo 3.
Thus querying more powers at this same prime adds no information.
A test based solely on reducing X^n modulo X²+X+1 similarly
depends only on n mod 3 by X³=1. This identifies a specific
quadratic-residuosity decision subproblem of the uniform evaluator;
the equivalence here does not establish a hardness result for the
fixed discriminant alone.

### Testing a single residue class by lifting to a larger modulus

For odd n=pq, p<q, let 1<ℓ≤n be coprime to n and r a unit mod ℓ.
On a positive instance q is in one of the two known classes r,
n/r mod ℓ. Fix a small rational ε>0 and take M=ℓ·2^t with t the
least nonnegative integer for which M≥n^(1/4+ε). Enumerate the
2^t lifts h in [0,M) of each class. The genuine lift h=q mod M
has x=(q−h)/M≤n^(β−1/4−ε), β=log_n q≥1/2. Modulo q it solves
the monic linear polynomial X+h/M with coefficients taken mod n.
As in the earlier precision argument, a constant grid of lower
bounds β₀ on β and β−β²≤1/4 put this root within the
unknown-divisor Coppersmith bound n^(β₀²−ε/2). Verify candidate
divisors and check the original bit; exhaustion certifies a negative
answer. This handles unbalanced semiprimes as well.

For fixed ε, this costs O(1+n^(1/4+ε)/ℓ) calls of polynomial bit
complexity. In particular the decision costs polynomial time for
ℓ≥n^(1/4+ε). For ℓ=n^α with fixed 0<α<1/4, choosing ε<α/2
improves the exponent of the preceding progression-product bound.
The small primes in the general reduction are only poly(log n);
this lift enumeration still has exponentially many candidates there.

### Removing the epsilon in the known-residue bound

A direct LLL estimate suffices for M≥16 n^(1/4), with gcd(M,n)=1.
Let q be the larger factor, h its least residue, β=log_n q, and
x=(q−h)/M. Check x=0 directly. Otherwise try dyadic bounds X up to
2n; one satisfies x≤X<2x≤n^(β²)/8 by β−β²≤1/4. Put
d=4 bitlength(n), and try integers a=1,...,d−1. One equals floor(dβ).
Use the coefficient rows of n^max(a−i,0) f(Z)^i for 0≤i<d with
column j scaled by X^j, where f(Z)=Z+(h/M mod n). At x every row
polynomial vanishes mod q^a. The triangular determinant is
n^(a(a+1)/2) X^(d(d−1)/2). LLL gives a nonzero polynomial H with
|H(x)| ≤ sqrt(d) 2^((d−1)/4) det^(1/d) by Cauchy–Schwarz.
This is strictly below q^a whenever

log₂ X < [(2aβ−a(a+1)/d)/(d−1)] log₂ n − 1/2 − log₂(d)/(d−1).

Indeed the bracket is at least β²−2/(d−1), by completing the
square and |a−dβ|<1. With this d the total loss on the right is
less than 3 bits, as required. Hence H(x)=0 in the integers.
Recover its integer roots, check h+Mx by division, and try both
possible residue classes. Dimension, entry bit lengths, grid size
and root extraction cost are polynomial in log n. This is a use of
the usual univariate small-root lattice, not a new lattice algorithm.

Doubling the query modulus until M≥16 n^(1/4) therefore costs
O(1+n^(1/4)/ℓ) lifts, each of polynomial bit complexity. It also
sharpens the earlier Jacobi-negative shear precision: 3^T≥16 n^(1/4)
suffices, so a single determinant answer needs (log₂ n)/4+O(1)
bits for arbitrary distinct odd semiprimes after the small-prime check.

### Keeping d fixed at quarter-bit precision by choosing the base

For n=pq with distinct odd primes let S=p+q. The fixed shear d=1
has trace F=n(4n+1−2S). For successive odd primes ℓ test gcd(ℓ,n)
and probe det(I+ℓV_1) mod ℓ². The returned value is 1+ℓF modulo
ℓ². Thus we can compute S mod ℓ, then test S²−4n mod ℓ. Stop at
the first nonzero discriminant. Every rejected prime divides p−q,
so the primorial bound guarantees success by ℓ=O(log n) unless a
factor was already found. No factoring oracle is needed to select
the base given the low-precision determinant probes.

Take the least T with M=ℓ^T≥16 n^(1/4); trial-divide through T.
If no factor is found, all nontrivial cycles vanish from the
determinant modulo ℓ^(T+1), leaving (1+ℓ)^F. For odd ℓ its
exponent is uniquely recoverable mod M by digit lifting, since
v_ℓ((1+ℓ)^(ℓ^j)−1)=j+1. This gives S mod M. Enumerate the two
simple roots of X²−SX+n mod ℓ and lift each uniquely to M. These
are the two factor residues, to which the preceding lattice argument
applies. The large query costs at most (log₂ n)/4+2 log₂ ℓ+O(1)
answer bits, in addition to the probes. The combined probe output is
O(log n) bits by the primorial bound. The coefficient of the matrix
varies with ℓ; the shear itself is fixed throughout.

### Reusing primes where the two factors coincide

Keep the product B of the odd probe primes where the discriminant
vanishes, starting at 1. Each such probe gives both factor residues:
p=q=S/2 mod ℓ. CRT gives a single common residue mod B, with no
choice of orientations. Stop immediately if B≥H=16 n^(1/4).
Otherwise, at the first separating prime ℓ use the least T≥1 with
B ℓ^T≥H. Recover the two simple roots mod ℓ^T as above (or use
the probe itself for T=1), then combine each by CRT with the same
residue mod B. The lattice argument applies to these two classes.
As before, trial division through T precedes the high precision query.

At termination the final combined modulus M lies in [H,H ℓ_last).
For an archived transcript, encode the index of the last prime and
whether it separated, then pack the local S residues in mixed radix.
Use one digit per collision prime and one residue mod ℓ^T at the
separating prime: its probe digit is already included. These residues
can be decoded back to the full determinant answers by the trace
formulas and exponentiation. Since ℓ_last=O(log n), the packed
transcript has at most (log₂ n)/4+O(log log n) bits, including
metadata, on paths reaching lattice recovery. Gcd or trial-division
exits already produce a factor. This bound concerns packed replies;
the actual adaptive evaluator still receives the low precision probes
and may receive a second request at the same prime.

### A single deterministic base-3 query when the discriminant is singular

Keep both d=1 and the determinant coefficient 3 fixed. For n=pq
with distinct odd primes other than 3, let h be least with
3^h≥16 n^(1/4), and v=v_3(p−q). First try the at most two square
roots of n modulo 3^h in the known-residue recovery algorithm.
These roots are computable by simple Hensel lifting (n is a unit).
If v≥h this already recovers a factor because q²=n mod 3^h.

After that attempt fails we know v<h. Put K=2h−1, trial-divide
through K, and query det(I+3V_1) modulo 3^(K+1). As above recover
S=p+q modulo 3^K from the fixed-point exponent. The discriminant
Δ=S²−4n has valuation exactly 2v<K, so the residue reveals v.
Divide Δ by 3^(2v) and take its two unit square roots modulo
3^(K−2v). Multiplication by 3^v recovers ±(q−p) modulo 3^(K−v).
Hence (S±(q−p))/2 gives the two factor classes modulo at least 3^h.
This avoids enumerating the exponentially many apparent quadratic
roots at the original precision. The single answer uses at most
(log₂ n)/2+O(1) bits. If a bound v≤w<h is known, K=h+w suffices
by exactly the same argument. In particular the extra precision over
the quarter-bit leading term is at most w ternary digits. This is a
sufficient precision bound, not a lower bound.

### What a Dirichlet GRH bound buys for base selection

For n=pq with distinct odd primes, the Jacobi character (·/n)
has a proper kernel on the units. Under GRH for Dirichlet L-functions,
the Ankeny–Bach least-nonmember bound yields a unit d with
(d/n)=−1 by deterministic scanning in O((log n)^2) candidates.
Gcd and Jacobi computations use only n. The earlier split/inert
shear then takes one base-3 determinant query of quarter-bit
leading precision. GRH enters only in bounding this search.

There is also a fixed-shear version. Use the quadratic character
u ↦ (n/u) on odd units modulo 4n (well-defined by quadratic
reciprocity). It is nontrivial: choose u=1 mod 4 with (u/n)=−1
by CRT. A small nonmember from the same bound has an odd prime
divisor ℓ with (n/ℓ)=−1 by multiplicativity. Equivalently scan odd
primes and test Euler's criterion, taking any gcd factor directly.
Then p and q cannot coincide modulo ℓ, since their product would
be a square. This chooses a separating base without determinant
probes. With d=1 a single answer of (log₂ n)/4+O(log log n) bits
suffices under the same hypothesis. Neither application bounds
the time for computing the determinant. The relevant analytic
assumption concerns Dirichlet L-functions. See Eric Bach, Explicit
bounds for primality testing and related problems, Math. Comp. 55
(1990), 355–380, https://doi.org/10.1090/S0025-5718-1990-1023756-8 .

### Deterministic separation on the original modulus under Dirichlet GRH

In the arbitrary-input Walsh reduction, Dirichlet GRH allows every
determinant query to use the same remaining modulus n and removes the
random choice of units. Strip small primes as in the exact-count
reduction. List the units u_j up to the uniform Ankeny–Bach bound
O((log n)^2), testing gcds and splitting recursively on a proper gcd.
One may alternatively strip all primes up to this bound first, using
the original input bit length for a uniform cutoff. For each p|n let
v_p be its vector of nonsquare bits on this list. Every vector is
nonzero and every pair distinct: the Legendre character at p and
the product of the characters at any distinct p,q are nontrivial
on the units modulo n by CRT. Apply the least-nonmember bound to
each kernel. Multiplicity does not affect this argument.

To avoid exponentially many subset products on this long list,
write V_p(X)=sum_j v_p[j] X^j over F_2 with indices starting at 0.
For m list entries and L input bits, choose a binary field of
size H>m(L+1)^2, rounding up to a power of two. A defining
irreducible polynomial can be found by exhaustive trial in time
polynomial in H. Enumerate evaluation points z in this field.
Some z makes all V_p(z) nonzero and distinct: the nonzero polynomials
V_p and V_p−V_q together exclude fewer than H roots. The vectors
themselves need never be known. In a fixed binary basis form one
unit per coordinate b by multiplying exactly those u_j for which
z^j has bit b set. These log₂ H units have prime signatures V_p(z).
For each z run the earlier exact-count Walsh extraction on its H
subset products. Verify candidate factors by primality and division,
and stop on a complete factorization. The good z guarantees success
within H iterations. The number of queries and all ordinary bit
costs are polynomial in L. Compared with deterministic interpolation
at n^k, this result keeps the input modulus fixed at the cost of a
GRH assumption and more queries. The determinant runtime is open.

### Falling factorials realize the interpolation values at fixed modulus

Let A_k(x)=x(x−1)...(x−k+1), A_0=1. Use the shear sending
z to z+y A_k(x), fixing x,y modulo n. Its pair count is exactly
Q(k)=sum_{x mod n} gcd(A_k(x),n). If p^a exactly divides n and
k<p, the k linear factors have disjoint zero classes modulo p.
Away from all those classes the gcd is 1. Within the class of i,
all factors except x−i are units, so the excess over 1 summed
on this class is G(p^a)−p^a. Consequently the local count is
p^a+k a(p−1)p^(a−1). CRT gives precisely the integer nP(k) used
in the power-query interpolation above, now on the original ring.

Concretely strip primes through max(L+1,T,3) with 3^T>N² and
T minimal, L=bitlength(N). For each k up to the earlier degree
bound D, det(I+3V_k) mod 3^(T+1) yields the exponent nQ(k)
modulo 3^T. Dividing by the unit n gives Q(k) exactly since
Q(k)≤n². Interpolate and extract factors as above. This is
deterministic and unconditional as an oracle reduction, with
fixed modulus and fixed determinant coefficient. The shear has
degree k+1 and an O(k)-multiplication description. It changes
the family beyond the earlier quadratic shears. Since k! is
coprime to n, gcd(A_k(x),n)=gcd(binomial(x,k),n) for canonical
representatives x≥0 (binomial(x,k)=0 for x<k). This gives a
concrete binomial sum whose uniform efficient evaluation would
suffice. No bound for computing that sum in poly(log n) follows.

### One count using idempotents and coefficient packing

Let L=bitlength(N), s=2L+1, K=2^s. For any remaining modulus n|N
use variables e,x,b_1,...,b_s and the quadratic equations
e(e−1)=0, (1−e)x=0, b_i(b_i−1)=0, (1−e)b_i=0 for every i.
At each prime power u=p^a exactly dividing n, the only idempotents
are 0 and 1. If e=0 all other variables are forced to zero. If e=1
then x is free and the s idempotents are independent. The local
count is 1+K u. Thus the global count C is F(K), where
F(X)=∏_{p^a || n}(1+p^a X). Every coefficient is less than K
since their sum is ∏(1+p^a)≤2^L n<K. Base-K digits give the exact
coefficients, and the negative reciprocal roots are precisely the
maximal prime powers. Integer root extraction over exponents 1,...,L
and primality testing recover their primes. Digit extraction and
rational polynomial factoring have polynomial bit cost.

For a single determinant query choose T=s(L+1)+1 and strip primes
through T before constructing the system. Introduce one translation
coordinate per equation (r=2+2s), shifted by its left-hand side.
The permutation has n^r C fixed points; every other cycle has length
dividing n and greater than T. From det(I+3V) modulo 3^(T+1), lift
the exponent of 4 and divide by the unit n^r modulo 3^T. The
result is exactly C because C<K^(L+1)<3^T. Thus this enlarged
family of vector translations needs only one O(L²)-bit answer for
general inputs. The price is O(L) coordinates and equations;
this is not a single query to the earlier three-coordinate shear.

### Packing the binomial counts into a first-order determinant

After stripping primes ≤L+1 (L the original bit length), take D any
known upper bound on the number of remaining distinct prime factors
with D≤L. Skip the case n=1. On D labeled copies of (Z/nZ)^3 use
the falling-factorial shears V_k for 1≤k≤D. Choose B a power of
two strictly above n³, and ask for the determinant of the block
diagonal matrix with blocks I+B^(D+k−1)V_k modulo B^(2D).
All products of two perturbation entries vanish modulo B^(2D).
The full determinant therefore equals
1+B^D sum_{k=1}^D B^(k−1) nQ(k) modulo B^(2D).
Each trace is at most n³<B, so the digits are exact. Subtract one,
divide by B^D, read the D digits and divide each by n. This is
ordinary integer division and needs no inverse mod B. Interpolate
with Q(0)=n. No trial division for cycle isolation or exponent
lifting is needed, only the small-prime check for interpolation.
The answer takes O(D log n) bits. This explicitly requests a
weighted block diagonal matrix; an evaluator restricted to one
unweighted shear cannot be assumed to accept that input. Indeed
the encoded answer is computable in polynomial time from the D
separate counts, so packing does not itself accelerate evaluation.

### Exact hit rate for direct sampling of the binomial sum

For n=pq with distinct primes and 1≤k<min(p,q), sample x uniformly
mod n. The conditions that x mod p and x mod q fall among 0,...,k−1
are independent. Thus gcd(binomial(x,k),n) takes values 1,p,q,n
on respectively (p−k)(q−k), k(q−k), k(p−k), k² inputs.
Among the k² inputs with gcd n, k have the same shift at both primes
(x=0,...,k−1). The other k(k−1) have distinct shifts and yield a
proper divisor by gcd(x−i,n) for an appropriate i. Altogether the
chance that any of the k individual shifts gives a proper gcd is
ρ=k(p+q−k−1)/n. The exact success probability after t independent
samples is 1−(1−ρ)^t. Querying all lower orders at the same x does
not enlarge the set of shifts. For balanced semiprimes p,q=Θ(√n)
and k polynomial in log n, constant success by this procedure takes
Θ(√n/k) samples. Also, with probability ((p−k)(q−k)/n)^t the
ordinary sample mean of the gcd values equals 1, while its target
mean is (1+k(1−1/p))(1+k(1−1/q)). Consequently polynomially many
uniform samples miss even a factor-two relative approximation with
probability tending to one for balanced growing semiprimes. This
conclusion concerns the raw sample mean and direct shifted-gcd tests.

### Divisor grouping and what the extreme coefficients reveal

For the rough modulus n and 1≤k<min_{p|n}p, expand the gcd by
gcd(t,n)=sum_{d|n, d|t} φ(d). The falling factorial has exactly
k^ω(d) roots modulo d: at each prime power dividing d exactly one
linear factor must vanish to the full prime-power precision.
Consequently Q(k)=sum_{d|n}(n/d)φ(d) k^ω(d), including d=1.
This groups the polynomial coefficients by number of distinct
primes in a divisor. Enumerating those divisors already requires
access to the factor structure.

Two consequences give simpler targets than all coefficients. Write
r=ω(n) and a_p=v_p(n). The leading coefficient is φ(n)∏ a_p,
an integer multiple of the group exponent λ(n) with polynomial bit
length. Given this coefficient alone, standard randomized order
splitting factors n in expected polynomial time: try uniform bases,
check gcd first, exponentiate by the odd part of the multiple and
square successively, testing gcds with n when encountering a
nontrivial square root of 1. Recurse on splits and handle prime
powers by integer roots and primality testing. For two distinct odd
prime divisors, the 2-parts of the orders of independent uniform
local units have different sizes with probability at least 1/2,
which suffices for splitting; the same multiple works recursively.

The linear coefficient divided by n is sum_p a_p(1−1/p).
Its reduced denominator is exactly rad(n) once p>a_p for all p:
after clearing the squarefree common denominator, exactly the term
−a_p rad(n)/p survives modulo p. This identifies the squarefree
kernel even with repeated primes, but for squarefree n this
denominator merely returns n. Neither coefficient formula supplies
an efficient evaluation procedure on unfactored inputs.

### A central-binomial decision predicate on binary-size inputs

Consider the single bit B_n(j)=[gcd(binomial(2j,j),n)>1]. Define
H_n(m)=[gcd(m!,n)>1]. For m=2j+b with b in {0,1}, the factorial
identity (2j)!=binomial(2j,j)(j!)² gives H_n(m) as the OR of
H_n(j), B_n(j), and (when b=1) [gcd(m,n)>1]. Start at H_n(0)=false.
This computes H_n(m) by scanning the binary prefixes in O(log(m+1))
oracle queries and ordinary gcds. For n composite, H_n is a monotone
threshold whose first true index is the least prime factor. Binary
search up to floor(sqrt n) finds it with O(log² n) central-binomial
bit queries. Remove its multiplicity, test primality and recurse.
All queried indices satisfy 2j≤sqrt n at that splitting step.

Conversely, given the factors, B_n(j) is tested by checking whether
j+j has a carry in base p for any p dividing n (Kummer's theorem,
also immediate by subtracting factorial valuations). Thus this
particular predicate, on binary inputs, is deterministically
polynomial-time Turing equivalent to factoring. Its output is just
a coprimality bit, including when the gcd would be all of n.
The recurrence avoids needing the numerical residue of the binomial.
Computing its full integer expansion can take Θ(j) bits and is not
a polylog implementation of the predicate.

### Logarithmically many promised binomial queries for a full factorization

Let p be the least prime divisor of the current composite n. At
j<p, the denominator (j!)² is invertible modulo n. Thus B_n(j)
from the preceding section is exactly [p≤2j]. More generally, for
m=2j+1 test B_n(j) or gcd(m,n)>1 to decide [p≤m] under the same
promise. This covers odd cutoffs without exceeding floor(sqrt n).

Start with a=1<p. While 2a<floor(sqrt n), probe at j=a.
A negative answer doubles the certified lower bound a; a positive
answer brackets p in (a,2a]. If no positive answer was needed,
the final bracket is (a,floor(sqrt n)] with upper endpoint ≤2a.
Binary-search this bracket using the cutoff test for m at its
midpoint. Throughout, floor(m/2)≤a for the initial bracket lower
endpoint a, so every query meets j<p. The doubling stage takes
O(log p) calls and the bracket width is at most p, giving O(log p)
calls for this split in total. Remove the full multiplicity of p
and repeat on the remainder, stopping on a prime or 1. Summing
log p over the distinct extracted primes gives O(log N) one-bit
queries for the entire factorization. All 2j≤floor(sqrt n).
The unresolved case is the promised bit itself: it asks whether a
prime divisor first appears by 2j given that none appears by j.

### Baseline cost and a limitation of symbolic interval products

Under the promise p_min(n)>j the query is exactly the nonunit test
for ∏_{h=1}^j(j+h) modulo n. A block polynomial with b=ceil(sqrt j)
consecutive factors evaluated at O(b) equally spaced starts, plus a
final partial block, computes this product in softly O(sqrt j)
ring operations by fast monic multipoint evaluation. Alternatively
factor n by Harvey's deterministic n^(1/5+o(1)) bit algorithm and
compare its least prime with 2j. The two upper bounds cross at
j=n^(2/5) at the exponent level. For j=n^α, 0<α≤1/2 fixed, the
combined bit bound is n^(min(α/2,1/5)+o(1)). These are benchmarks,
not lower bounds. See https://arxiv.org/abs/2010.05450 and the
refinement https://arxiv.org/abs/2105.11105 .

A possible shortcut would compute P_j(X)=∏_{h=1}^j(X+h) itself by
a short constant-free division-free straight-line program over Z.
It has j distinct integer roots. The Shub–Smale τ conjecture bounds
the number of distinct integer roots of a nonzero integer polynomial
by a fixed polynomial in its straight-line complexity. Under that
conjecture the proposed symbolic program needs j^Ω(1) operations.
The argument does not apply to the single specialized input X=j,
modulus-specific algorithms, or algorithms returning only the gcd bit.

### A polylogarithmic decision window below the square root

For odd composite n in the promised query, put t=2j≤√n and assume
t³≥n. A negative answer requires every prime factor to exceed t;
then n must have precisely two prime factors counted with multiplicity.
If n=pq with t<p≤q, the integer a=(p+q)/2 lies between ceil(√n)
and floor((n+t²)/(2t)), since x+n/x decreases up to √n.
Enumerate just these a and test whether a²−n is a square b².
Return negative if a−b and a+b are both primes greater than t;
return positive if the range is exhausted. Deterministic primality
tests make this valid even when the input is not known to be
semiprime. Includes p=q. Even inputs and primes are handled directly.

The interval has at most 1+(√n−t)²/(2t) candidates. Thus for each
fixed c≥0, if 0≤√n−2j≤n^(1/4)(log n)^c the promised decision is
deterministic polynomial time. For sufficiently large n this window
satisfies t³≥n. Exhaustion can certify a positive bit without
identifying the factor. The window does not cover general thresholds
that occur while bracketing the least prime divisor.

### Batching the Fermat decision window with exact randomized answers

Keep t³≥n and the previous interval [A,B] for a, of size W.
An empty interval certifies positive. Handle prime powers directly by
integer roots and primality tests. Otherwise n has at least two
distinct prime divisors. Choose g uniformly in 1,...,n−1, using a
proper gcd immediately. For a unit g, seek any a in [A,B] satisfying
(g²)^a=g^(n+1) mod n by ordinary interval baby-step/giant-step:
store the first ceil(sqrt W) powers and shift the target by inverse
giant steps, offsetting by A. Keep the smallest baby index for
duplicate values and check the last block's endpoint. If no match
exists the answer is positive: on every negative instance the true
a=(p+q)/2 gives exponent n+1−2a=φ(n). If there is a match, its
positive exponent E=n+1−2a kills g. Strip powers of 2 from E and
use repeated squaring with gcds against n to extract a proper split.
If extraction fails, retry with fresh g. Any proper split suffices
to classify the bit: it is negative iff both parts are primes >t.

For any two distinct odd primes dividing n, the local orders of a
uniform unit have independent 2-adic valuations. In a cyclic group
of even order no single such valuation has probability above 1/2.
Hence with probability at least 1/2 the two valuations differ.
On this event *every* killing exponent of g yields a split, since
at some squaring step one component is 1 and the other is not
(take gcd at each step). Thus each trial gives an exact decision
with probability at least 1/2 on either kind of input, and the
expected number of trials is at most two. Expected time is
O((1+sqrt W) poly(log n)), space O((1+sqrt W) poly(log n)) bits.
This replaces W individual square tests by a square-root search;
the polylogarithmic window still requires W=poly(log n).

### Unique matches and deterministic bases for the batched test

In the interval search keep the two smallest baby indices for each
group value. For each giant block collect up to two legal matches
within its endpoint; stop on a second distinct match. This
distinguishes zero, one and multiple matches in the same square-root
bound. A unique match is decisive without order splitting: test
a²−n for a square and verify that both resulting factors are primes
greater than t. Every negative instance supplies its true average,
so failure of that check certifies positive. In particular a unit g
with ord_n(g²)≥W always settles the bit in one search.

For a deterministic guarantee under Dirichlet GRH, handle prime
powers first, then scan the small units from the Ankeny–Bach bound.
To justify the scan pick any two distinct odd prime divisors p,q
in the proof. Write h_p=v₂(p−1), h_q=v₂(q−1). If h_p>h_q,
any quadratic nonresidue at p has local order of 2-adic valuation
h_p, greater than any local valuation at q. Interchange p,q for
the reverse inequality. For equal h, any unit with opposite Legendre
symbols at p and q gives unequal valuations. In each case the good
units include the complement of the kernel of one nontrivial
quadratic character modulo n. Under Dirichlet GRH a witness is
therefore reached in O(log² n) candidates. The algorithm need not
know that character: use gcd, interval search, and, for multiple
matches, the earlier order-splitting test on any one match.
A witness either finds no match or settles the bit, even in the
multiple-match case. This gives deterministic
O((1+sqrt W) poly(log n)) time conditionally. The unique-match
certificate itself needs no analytic assumption.

### Unconditional deterministic batching by simultaneous bases

For W≤n^(2/5) take all prime bases ≤y=ceil((log₂ n)^8),
handling proper gcds immediately and stopping the list before n.
Perform one interval baby-step/giant-step search on tuples of powers
for these bases, using sorted tuple keys for deterministic lookup.
On a negative instance n=pq (handle squares beforehand), the true
average matches simultaneously. Two matches at distance d<W imply
g^(2d)=1 modulo n for every base and hence every product of them.

There is an explicit certificate that this cannot occur. Take
m=floor(log_y(n−1)), and s the number of prime bases. The binomial
coefficient binom(s,m) counts distinct products of m distinct bases,
all ≤n−1 and therefore distinct modulo n. On a squarefree semiprime
at most (2d)^2 residues solve X^(2d)=1 by the root bound separately
in each prime field. Consequently if binom(s,m)>4W² then *multiple*
simultaneous matches certify positive. Zero and unique matches use
the previous exact decisions. Only the cardinality is computed;
the products are never enumerated.

For sufficiently large n the certificate holds throughout this range:
s=Ω(y/log y), m=(1+o(1))log n/log y and binom(s,m)≥(s/m)^m
give n^(7/8−o(1)) products, whereas 4W²≤4n^(4/5).
Use integer comparisons (for example y=bitlength(n)^8); when the
certificate fails fall back to the exact Fermat scan. The exceptions
in this range occur only for bounded n, without needing a numerical
cutoff for the proof. Time in this range is O((1+sqrt W) poly(log n)).
For larger W simply factor n by the Harvey–Hittmeir deterministic
O(n^(1/5) poly(log n)) algorithm and classify directly; its bound
is within O(sqrt W poly(log n)) there. Thus the deterministic
square-root-in-W upper bound holds unconditionally on the previous
t³≥n domain. The many-base proof can have substantial overhead.
It does not yield polylog time for general W.

### Sieving the baby indices by the Fermat discriminant

A true average a satisfies a²−n=b² over the integers. For an odd
prime ℓ not dividing n the permitted a mod ℓ form a set of size
(ℓ+χ_ℓ(n))/2, with zero allowed as a square. Indeed
(a−b)(a+b)=n gives ℓ−1 ordered solutions; the 1+χ_ℓ(n)
solutions with b=0 are counted only once on projecting to a.
Use consecutive odd primes, taking any proper gcd directly. For the
product M of the selected primes the CRT gives R admissible classes,
the product of their local sizes. Generate those classes directly
by CRT. For a window [A,A+W), turn them into offsets s=(a−A) mod M.
Store the tuple powers at the R offsets, and giant-step by M over
ceil(W/M) blocks, rejecting offsets outside the final endpoint.
The previous simultaneous-base certificates apply to this subset:
every true negative still supplies its true average. Keep the two
least offsets at each duplicate baby key as before.

Stop at the first M with MR≥W. Then W/M≤R, and MR≤W times the
square of the last prime. If k primes are used, R/M≤(2/3)^k.
By the primorial estimates the last prime is O(log W) and
k=Ω(log W/log log W) as W grows: M≥sqrt W, and the previous
primorial was below W. Thus R+ceil(W/M) is at most
sqrt(W) exp(−Ω(log W/log log W)) up to polynomial logarithmic
factors. CRT generation, powering and sorted-key lookup add only
poly(log n) overhead per entry. This bound holds deterministically
in the simultaneous-base range W≤n^(2/5), and in expected time
for the single random-base method throughout its original domain.
It saves a subpower factor while leaving exponent 1/2 in W.

### Limits of the odd-modulus discriminant sieve

Let R(M) count the classes a mod M for which a²−n is a square,
including zero, with M odd and coprime to n. For M=p^e let
z=1+(n|p). The exact count is
((p−1−z)/2)p^(e−1) + z [1 + sum_{1≤j,2j<e} φ(p^(e−2j))/2].
The first term lifts the classes with nonzero square discriminant.
On each of the z singular classes the map a ↦ a²−n is bijective
onto multiples of p modulo p^e by simple-root lifting in a.
The bracket counts square multiples of p by their even valuation.

For p≥5 the first term alone gives density at least 1/5. For p=3
and z=0 it gives 1/3. For p=3 and z=2, the two singular classes
give density at least 2/27: at e≥3 require a²−n=9 mod 27 so the
remaining unit after division by 9 is a square modulo all higher
powers of 3. At e≤2 the zero discriminant already suffices.
CRT now gives R(M)/M ≥ (2/27)^ω(M). Since M exceeds or equals
the product of the first ω(M) primes, ω(M)=O(log M/log log M).
Consequently R(M)≥M exp(−O(log M/log log M)).

For the specific search that enumerates every admissible baby class
and every giant block, the entry budget R(M)+ceil(W/M) is therefore
at least sqrt(W) exp(−O(log W/log log W)) for 1≤M≤W. This follows
by taking the geometric mean of the two terms. For M>W, building
the entire baby table already costs W^(1−o(1)). Together with the
preceding construction this pins down the scale of savings from
choosing odd wheel moduli, including prime powers. The bound concerns
full residue tables; implicit representations and interval-specific
pruning require a separate analysis.

### A parity floor-sum target for the promised interval

Let F(m,j,ℓ,r) be the parity of the sum of floor(m/x) over
j<x≤2j with x=r mod ℓ. XOR the answers at m=n and m=n−1.
Their difference counts, modulo two, the divisors of n in that
slice: the two floors differ exactly when x divides n. Under
p_min(n)>j, every divisor in the interval is prime, since a product
of two primes greater than j exceeds 2j for j≥1. For semiprimes
with 2j≤√n there is at most one such divisor (also for squares).
Hence in this case the unfiltered parity, ℓ=1, suffices.

For general inputs let L=bitlength(n) and S be the set of prime
divisors in the interval. If S is nonempty, fix one p in S in the
proof. A prime ℓ failing to isolate p must divide the product of
|p−q| for q in S other than p. That product is less than or equal
to n^L, including the empty product. Enumerate consecutive primes
until their product exceeds n^L; some ℓ isolates p and its residue
class returns odd parity. Try all residues at all these primes.
An empty S always returns zero. The largest prime tested is O(L²)
by the primorial bound, so the number and sizes of calls are
polynomial in L. For example n=5·7·11, j=4 has two hits at 5 and 7:
unfiltered parity vanishes and the two classes modulo 3 separate them.

Combined with the promised binomial search this is a deterministic
reduction of general factoring to parity of truncated reciprocal
floor sums on small-modulus progressions. It requires no factorization
for the reduction steps. A polylog evaluator for these parity sums
is still missing; enumerating the terms would cost up to Θ(j).

### Unfiltered interval parity factors every squarefree input

Let H(n,t) denote the parity of the number of divisors in [1,t].
To obtain it from unfiltered dyadic floor-sum parities, put j=floor(t/2).
XOR H(n,j) with the dyadic floor-sum difference at n and n−1 and,
when t is odd, the directly tested bit [t divides n]. Use H(n,0)=0
and handle t=1 directly. Each prefix therefore needs O(log(t+1))
unfiltered queries, even when j exceeds the least prime factor.

For nonsquare n, divisor pairing shows that [2,floor(sqrt n)]
contains precisely τ(n)/2−1 divisors. If 4 divides τ(n), this is
odd. Its parity is directly observable, so the algorithm needs no
prior knowledge of τ(n). Bisect this interval, querying child parity
from prefix differences and keeping an odd child. The final singleton
is a verified proper divisor. For squarefree composite n, τ(n)=2^r
with r≥2, and both sides of every split are squarefree. Recursive
splitting with a primality test therefore factors any squarefree n
using polynomially many unfiltered calls. Per split the bound is
O(log² n) oracle calls. There is no progression restriction on x.
The oracle must allow dyadic intervals up to sqrt n, including ones
containing composite divisors.

For general n, exact square roots handle squares. The same odd-child
argument succeeds on nonsquares with 4 dividing τ(n). If it fails,
n has precisely one odd prime exponent, congruent to 1 mod 4, all
other exponents even. This follows from v₂(τ(n))=1. For example
45 has the two divisors 3 and 5 in the initial interval, so its
parity is even. The previous progression-filtered reduction still
covers these cases.

### A few low bits remove cancellation on all promised inputs

Use the least-factor search's original promise j<p_min(n). Every
divisor in (j,2j] is prime. If there are h of them, their product
divides n and is at least (j+1)^h. Thus h≤H=floor(log_{j+1} n),
computable by integer comparisons. Put b=ceil(log₂(H+1)).
Two unfiltered sums of floor(m/x) on (j,2j], reduced mod 2^b
at m=n and n−1, recover h *exactly* by modular subtraction.
No modulus filters, power inputs, or squarefree assumption are
needed. The j=0 case in an odd-cutoff test is trivial. The inputs
to the evaluator must include m=n−1 even if that integer itself
has a small prime factor; the promise is on n.

For L=bitlength(N), each answer uses at most ceil(log₂(L+1))
bits, and the preceding doubling-and-bisection reduction needs
O(L) pairs of answers for a full factorization. For j≥n^α with
fixed α>0 the required precision is constant. An arbitrary fixed
precision cannot recover h by this direct subtraction on all scales:
choose 2^b distinct primes in (j,2j] for sufficiently large j
and take their product n (b≥2). This satisfies 2j≤sqrt n and
has exactly 2^b hits, giving zero modulo 2^b. The prime number
theorem guarantees enough primes for each fixed b. This only
limits the single difference as a detection rule.

### A restricted exact floor-sum baseline via rational slopes

For m≥1 and 1≤j≤sqrt(m), the dyadic sum admits an
O(min(j,m^(1/3)) poly(log m)) bit-operation baseline. Below twice
the ceiling cube root use direct division. Otherwise take blocks
[c,c+w) inside (j,2j] of length w≤h=floor(j/ceil(m^(1/3))).
Continued fractions provide integers a,q with 1≤q≤h and
|a−qm/c²|≤1/h. Put f(x)=ax+qm/x. For x=c+d in the block,
f(x)−f(c)=(a−qm/c²)d+qm d²/(c²(c+d)).
Since c≥j, d<h, q≤h and mh³≤j³, this lies between −1 and 2.

Set C=floor(f(c))−1. The baseline sum of floor((C−ax)/q) is
an affine floor sum computable by the Euclidean algorithm. To
correct it consider s=C+1,...,C+3. Add the count of block integers
satisfying ax=s modulo q and ax²−sx+qm≥0. Each such x adds the
ordinate y=(s−ax)/q between the baseline and m/x. The quadratic
inequality cuts the block into at most two intervals, whose integer
endpoints and congruence counts take poly(log m) bit operations,
also for linear or constant cases. There are O(m^(1/3)) blocks.
Taking the sum modulo a small power of two changes the accumulator
size but leaves this block count.

Richard Sladkey gives an O(m^(1/3))-time algorithm for the full
divisor summatory function by successive approximation:
https://arxiv.org/abs/1206.3369 . The full sum modulo two equals
floor(sqrt m) modulo two by pairing (x,y) with (y,x); truncation
to j<x≤2j prevents that cancellation in general.

### Adaptive straight-line blocks still need cube-root many pieces

Take j=floor(sqrt(m)/4) and consider covers of (j,2j] by consecutive
integer blocks. On each block allow an arbitrary integer a and positive
integer q, requiring the oscillation of f(x)=ax+qm/x at its integer
points to be at most a fixed constant C. This models the preceding
bounded correction-level method. Even adaptive covers require Ω_C(m^(1/3))
blocks for sufficiently large m.

Write H for the endpoint distance of a long block. Here f''(x) is
bounded above and below by positive absolute constants times q/sqrt(m).
The centered second difference at equally spaced integer points,
with spacing floor(H/2), is at most 2C and at least a constant
times qH²/sqrt(m). Hence q=O_C(sqrt(m)/H²). The endpoint secant
slope has magnitude at most C/H. By the mean value theorem and the
upper curvature bound, every x in the block therefore satisfies
|m/x²−a/q|=O_C(1/(qH)).

Set H₀=K m^(1/6), with fixed K sufficiently large in terms of C.
Every point of every block of length at least H₀ (endpoint distance)
has slope m/x² near one of the fractions just described with
q≤Q=O_C(sqrt(m)/H₀²). In the bounded slope range there are O(q)
possible numerators per q. The total length of the resulting slope
intervals is O_C(Q/H₀)=O_C(K^(−3)). The inverse slope map has
derivative O(sqrt(m)), so these intervals cover at most
O_C(sqrt(m)/K³)+O(Q²) integer x, allowing one rounding error per
interval. Here Q²=O_C(m^(1/3)), negligible compared with j.
Choose K so that long blocks cover at most half the integers.
Covering the rest in blocks of at most H₀+1 points takes
Ω_C(sqrt(m)/H₀)=Ω_C(m^(1/3)) pieces, matching the upper bound.

This bound concerns enumeration of constant-oscillation affine blocks.
It does not cover algorithms combining corrections across blocks,
allowing many correction levels with implicit evaluation, or computing
only the modular difference of the two sums.

### Sparse additive Fourier approximation of the floor-sum difference

Write A_j(m)=sum_{j<x≤2j} [x divides m]. Expand each divisibility
indicator as the average of exp(2πiam/x) over a modulo x. For every
reduced frequency with denominator d in (j,2j], the coefficient is
exactly 1/d: no second multiple of d occurs in the band. Their
total squared coefficient mass is
sum_{j<d≤2j} φ(d)/d² = (6/π²)log 2 + o(1), about 0.4214,
by partial summation of sum_{d≤u} φ(d)=3u²/π²+O(u log u).

This gives a finite-input obstruction. Use the dictionary of all
reduced rational frequencies modulo 1 of denominator at most 2j.
On any consecutive T=j³ integers its normalized Gram matrix has
diagonal one and off-diagonal row sums in absolute value
O(j² log j/T). Indeed the circular separation is at least 1/(4j²),
and the exponential-sum bound at distance δ is at most 1/(2Tδ).
Sum over neighboring frequencies in both directions. Hence mean
squared error of any coefficient approximation is at least
(1−o(1)) times the squared coefficient error. Retaining k frequencies
with arbitrary weights can remove at most k/j² of the indicated
mass from omitted coordinates. Uniform absolute error <1/2 thus
requires k≥((6/π²)log 2−1/4−o(1))j². One can take m starting
at 4j², so all these input integers have O(log j) bits.

This bound concerns a single linear Fourier approximation, with
frequencies and coefficients fixed across the full interval of m.
It leaves open grouped or nonlinear evaluations, frequencies outside
the dictionary, and approximations on the promised rough inputs alone.

### Fourier size on the rough-input CRT family

For fixed j consider one linear Fourier sum, frequencies of reduced
denominator at most 2j, required to approximate A_j(m) within 1/2
(strictly) on all composite m with p_min(m)>j and 4j²≤m. Its minimum
number of frequencies is Θ(j²/log j). The upper bound expands only
the indicators for primes j<p≤2j, since no other divisor can occur.

For the lower bound, in lcm(1,...,2j) fix the residue 1 at every
prime-power component belonging to a prime ≤j. Leave each band
prime p free. CRT realizes the full product of these coordinates.
Impose divisibility by one further prime in (2j,4j) and take a
representative exceeding both that prime and 4j². This gives
promised composite inputs with O(j) bits by the usual Chebyshev
bound on the logarithm of the lcm.

On this family every dictionary frequency is constant or depends
on exactly one free coordinate (a denominator containing p must
equal p). Take the real part, at most doubling the support. Write
the error as a constant plus local real errors e_p, each relative
to the point indicator of residue zero. Uniform error <1/2 implies
sum_p osc(e_p)<1: independent coordinates realize all local maxima
simultaneously and likewise all minima. If k_p nonconstant modes
are present at p, Parseval gives variance at least
(p−1−k_p)/p², since every missing nonconstant coefficient of the
indicator is 1/p. The oscillation is at least the standard deviation.
Thus a prime with fewer than (p−1)/2 modes consumes Ω(1/sqrt(j))
of the total oscillation. Only O(sqrt(j)) primes can do so.
The prime number theorem gives Θ(j/log j) band primes, and the
others each need Ω(j) terms. This proves the bound.

The CRT witnesses can have Θ(j) bits. The proof leaves open rough
inputs of O(log j) bits, grouped evaluations and decision-only tests.

### A single Fourier threshold for the yes/no query

On the preceding rough CRT family, the same Θ(j²/log j) bound
holds for a real linear Fourier score followed by one threshold
deciding whether *any* band prime divides m. Allow arbitrary real
weights and arbitrarily small separation margins. On the finite
product domain a correct threshold strictly separates the maximum
no-hit score from the minimum hit score, adjusting the threshold
between them if necessary.

Write the score as a constant plus sum_p f_p(m mod p). Subtract
from each f_p its minimum on the nonzero residues, absorbing the
offset in the constant. Let b_p be its resulting maximum on nonzero
residues, D_p=f_p(0), and B=sum_p b_p. The largest all-miss score
uses the maxima simultaneously. A single hit at p with all other
coordinates minimizing has score D_p plus the same constant.
Thus D_p>B for every p. Replace each local score by f_p/D_p.
It now gives exactly 1 on a hit and a value in [0,b_p/D_p] on
a miss, using the same nonconstant frequencies. The total excess
over the true count is in [0,E] where E=sum_p b_p/D_p<1 (also
when B=0). Subtracting E/2 gives uniform error <1/2. The preceding
Parseval bound on this CRT family therefore applies. Conversely,
thresholding the exact prime-indicator sum gives the upper bound.

The rescaling proves existence for the counting lower-bound argument;
it need not be computed by the evaluator. The restriction is one
fixed additive score and one threshold. Nonlinear combinations
of separate prime tests, input-dependent scores and the polynomial-
in-j input range are separate cases.

### Nonlinear residue queries and a single gcd witness

On the same rough CRT family, allowing arbitrary nonlinear postprocessing
of residue queries with modulus ≤2j drops the deterministic worst-case
query count to exactly π(2j)−π(j), when the input is accessible only
through these queries. Test each band prime for the upper bound.
On an all-miss path, any unqueried band prime can have its coordinate
changed to zero without affecting the transcript: among moduli ≤2j
only that prime itself can see its coordinate. This works for adaptive
queries too. This query model excludes direct arithmetic on the input.

For a different model let a fixed positive integer U(j) serve as a
gcd witness for every promised composite input m: gcd(m,U)>1 iff
some band prime divides m. Every band prime p must divide U. Indeed
choose one prime q in (2j,4j) with pq≥4j², or simply choose q in
(4j,8j). The negative input q² forces q not to divide U; the positive
input pq forces p to divide U. Both inputs have O(log j) bits.
Similarly every prime q>2j is forbidden in U by q². These conditions
are sufficient, since factors ≤j cannot divide promised inputs.
The smallest positive witness is consequently the band primorial.
Its binary length is Θ(j) by the prime number theorem. Moreover its
gcd with m gives exactly the squarefree product of the band hits.

The central binomial coefficient is another witness: each band prime
has valuation one and larger primes have valuation zero. The explicit
length lower bound leaves open a uniform short division-free arithmetic
circuit for any admissible witness. Such a circuit could be evaluated
gate by gate modulo m and would suffice without materializing U.

### Limits of products of short-base power differences as gcd witnesses

Consider a positive witness of the previous type that is a product of
terms a^k−1 or a^k+1, with a≥2, k≥1 and log a≤(log j)^C for fixed C.
Each term must have all its prime divisors ≤2j. Set K=k for the
minus sign and K=2k for the plus sign. For K>6, Zsigmondy's primitive
divisor theorem for a^K−1 supplies a prime of multiplicative order
K for a. For K=2k that prime divides a^k+1. Its size is at least
K+1, hence K<2j. Small exceptional exponents cause no asymptotic
change.

Group the band primes p dividing a term by d=ord_p(a), a divisor of
K. Each group has at most d log(a)/log(j) members by taking their
product inside a^d−1. It also has at most 1+j/d members since all
are 1 mod d. Splitting d at sqrt(j), each group thus has size
j^(1/2+o(1)). Uniformly for K≤2j the divisor bound τ(K)=j^o(1)
as an upper bound gives at most j^(1/2+o(1)) hits per term.
The full band contains Θ(j/log j) primes. Hence every witness of
this form needs at least j^(1/2−o(1)) distinct terms. Raising terms
to powers or multiplying by j-smooth padding cannot add band hits.
The argument applies to universal, input-independent positive gcd
witnesses and bases of the stated bit length. General circuits can
use additions between products or succinct bases with much larger
integer bit length.

### One-sided prime coverage when extra divisors are allowed

Suppose a product of terms a^k±1 (a≥2, k≥1) only has to be divisible
by every prime in (j,2j]. Permit arbitrarily large extra prime
divisors. Measure input size S by the sum of bit lengths of the
listed bases and exponents. On infinitely many dyadic bands, coverage
alone forces S=Ω(j^(1/3)).

Fouvry's shifted-prime result gives a positive proportion of primes p
with a prime r|p−1 exceeding p^(2/3). Consequently infinitely many
dyadic bands contain Ω(j/log j) such primes: otherwise summing over
bands contradicts the positive proportion. Fix one such band. A hit
by one term has d=ord_p(a) dividing K=2k. If r divides K, there
are O(log(2K)/log j) possibilities for r>j^(2/3) in that term.
Each supplies at most 1+j/r possible p in the band. If r does not
divide K, then d divides (p−1)/r and hence d≤(2j)^(1/3).
All hits in this latter case divide the product of a^i−1 for
1≤i≤floor((2j)^(1/3)). Its logarithm is O(j^(2/3) log a),
so it accounts for O(j^(2/3) log a/log j) distinct band primes.
Summing both bounds gives O(j^(2/3) S/log j) special hits; comparison
with the required number proves the claim. A term divisible by p
cannot have p dividing its base, so all orders used here exist.

This also applies if an integer witness divides such a product:
exact divisions can remove primes but cannot cover a missing prime.
The lower bound charges explicit integer descriptions; arbitrary
straight-line programs can reuse gates or succinctly generate large
bases. Source for the shifted-prime bound: Fouvry, Théorème de
Brun-Titchmarsh; application au théorème de Fermat (1985); see also
https://mathoverflow.net/questions/458013/ .

### Input-dependent order exponent for semiprimes with shared group factors

For n=pq with distinct odd primes set M=lcm(1,...,B) and E=(n−1)M.
Write h=gcd(p−1,q−1), u=(p−1)/h and v=(q−1)/h. Since
gcd(n−1,p−1)=h=gcd(n−1,q−1), p−1 divides E exactly when u
divides M, and likewise for v and q. The algorithm uses only n and B.
If either condition holds then one random trial splits n with probability
at least 1/2. Take a uniform base in [1,n−1], checking its gcd, then
write E=2^s e with e odd. Test gcd(g^(2^i e)−1,n) for i=0,...,s
by modular powering followed by squaring.

Indeed, if only one local group is annihilated by E, the fraction
of units in the other killed by E is gcd(E,q−1)/(q−1)≤1/2
(interchange primes if needed). The final gcd succeeds otherwise.
If both groups are annihilated, the independent local order
2-valuations differ with probability at least 1/2 and some stage
separates them. Nonunits sampled here give a proper factor directly.
For B polynomial in log n, the exponent and all arithmetic have
polynomial bit cost. The guarantee requires a residual cofactor
dividing M; it does not require all prime factors of p−1 or q−1
to be small.

For n=6533=47·139 and B=3, both local group orders contain 23,
and their residual cofactors after removing h=46 are 1 and 3.
Base 2 gives gcd(2^9798−1,n)=n but already at the odd exponent
4899 the gcd is 47. The ≥1/2 bound was also checked by enumerating
every base on all 56 eligible pairs of distinct odd primes below 48
with B=3. This is a special case of order-based factoring; no
condition covering arbitrary semiprimes is established.

### Signed residual cofactors with a Jacobi-selected exponent

For n=pq as above choose a uniform unit D, put J=(D/n) and use
E=(n−J) lcm(1,...,B). Work in (Z/nZ)[w]/(w²−D).
At p let σ=(D/p), and at q let τ=(D/q), so J=στ. For a uniform
unit z the quotient g=z/conj(z) is uniform in the norm-one group:
locally the split algebra gives (u/v,v/u) of order group size p−1,
while in the field case conjugation is Frobenius and the power map
z↦z^(1−p) has image of size p+1. Both local groups are cyclic and
the two samples are independent by CRT. Sample z by drawing its two
coefficients and testing its norm for invertibility. Return any proper
gcd. For fixed D the acceptance probability is bounded below by an
absolute constant even in the split algebras. Resample z at fixed D
on rejection so that conditioning does not bias the sign probabilities.

Set h=gcd(p−σ,q−τ). Since n−στ=σ(q−τ) mod (p−σ),
gcd(n−J,p−σ)=h and likewise at q. Therefore E kills the entire
local group at p iff (p−σ)/h divides lcm(1,...,B). If either
residual cofactor divides this lcm, the earlier success proof applies:
split at the final exponent if just one local element is killed,
or along the repeated-squaring chain when order valuations differ.
For an element a+bw the identity test returns gcd(n,a−1,b).
The probability of splitting at the chosen signs is at least 1/2.

If *any* of the four sign pairs meets the criterion, a uniform D
selects that pair with probability 1/4, giving at most eight outer
trials in expectation (or earlier success from a nonunit). This
extends the special-input guarantee while keeping the exponent at
O(log n+B) bits. The arithmetic is the norm-one group used in
p+1 factoring. For n=6439=47·137 and B=3 the unsigned residuals
23 and 68 both fail. For D=3 the local signs are + and −, the
orders are 46 and 138, and h=46. With z=1+w the identity gcd at
the odd exponent 2415 gives 137 by direct quadratic-ring powering.

### A batched second stage for a surviving signed order

Keep the signed exponent E above, with M=lcm(1,...,B). On the
chosen sign pair the image group at p after powering by E has
order R_p=u/gcd(u,M), where u=(p−σ)/h. This follows by dividing
p−σ and n−J by their exact gcd h. Define R_q similarly. If
min(R_p,R_q)≤S, there is a square-root-cost search at bound S.

Write v=g^E, b=ceil(sqrt S), and T_i=v^i+v^(−i), a scalar since
conjugation inverts v. Build the monic polynomial
P(X)=∏_{0≤i<b}(X−T_i) mod n. Evaluate at T_(jb), 1≤j≤b, by
monic multipoint evaluation. Generate both trace lists by the
second-order trace recurrence. If a local order d of v is ≤S,
write d=jb−i in these ranges. The corresponding evaluation has
nonunit gcd with n. An individual nonunit difference of traces can
be found by scanning the b babies at one nonunit evaluation.

Take a proper gcd directly. If the individual difference vanishes
at both primes, try the two positive multipliers jb−i and jb+i.
Locally equality of traces implies equality or inversion of the
corresponding powers: in the field model they are roots of the
same quadratic, and in the split model use the first coordinate.
Thus at least one multiplier times E kills the local g at one
prime. If it kills only one, the identity gcd splits n; if it
kills both, the 2-adic squaring chain succeeds whenever the two
original local order valuations differ. This latter event alone
has probability at least 1/2 for uniform g and suffices for any
chosen nonunit trace difference. Including random D, any qualifying
sign pair gives at most eight expected outer trials. With B polylog n
each trial costs (1+sqrt S) poly(log n,log(S+1)) bit operations
by fast polynomial arithmetic. A polylog S still gives expected
polylog factoring for this restricted family. For such S one could
also enlarge stage one's lcm bound to BS: the entire residual u
then divides lcm(1,...,BS) since each prime power in u is ≤BS.
The second stage reduces the dependence on the search bound.

### A simultaneous obstruction to all four signed orders

Choose four distinct primes r₁,...,r₄ in [R,2R] for large R.
On these four moduli prescribe residues (1,−1,2,2) for p and
(2,2,1,−1) for q. CRT gives two distinct reduced classes modulo
their product. Linnik's least-prime theorem gives primes in both
classes of size R^O(1); hence for n=pq we have R≥n^c for an
absolute c>0 for sufficiently large R. Both primes tend to infinity.

For every σ,τ in {1,−1}, one of r₁,r₂ divides p−σ, and q is 2
modulo that prime. Thus n−στ is σ(2−τ) there and is nonzero.
Symmetrically one of r₃,r₄ divides q−τ and does not divide n−στ.
For any B<R both private primes survive E=(n−στ)lcm(1,...,B).
Both residual group sizes consequently exceed or equal R for all
four sign pairs at once. For every fixed polylogarithmic B,S these
bounds beat S on infinitely many n.

There is also a probability bound for the specific second-stage
test. After sampling a uniform norm-one element and powering by E,
each local image element is uniform in its cyclic image group.
The chance its order omits the designated prime r is at most 1/r.
With probability at least 1−2/R both local orders have their
designated prime divisor. For b=ceil(sqrt S), if b²+b<R,
no positive multiplier jb±i tested in either stage is divisible by
either designated prime. Thus no trace hit or identity gcd occurs
on this event. This bound conditions on sampling units, as in the
signed method. Incidental nonunits during sampling cost O(1/p+1/q)
in probability per draw. The claim applies to these exponents and
this bounded collision search, not to other factoring algorithms.

At (r₁,r₂,r₃,r₄)=(11,13,17,19) the least positive prime
representatives are p=1064933 and q=393109. For B=10 the four
residual pairs, in (++,+−,−+,−−) order, are (266233,32759),
(266233,39311), (481,799), and (59163,39311).

### Simultaneous avoidance of length-indexed polynomial exponents

For each bit length L allow a predetermined list of nonzero integer
polynomials P(X), with list size, degrees and coefficient bit lengths
bounded by fixed polynomials in L. Even trying all |P(n)| as
exponents, times an lcm cutoff B=poly(log n), leaves both signed
local groups with polynomially large residuals on infinitely many
semiprimes. This includes lists containing products or powers
provided the resulting degree and coefficient bounds hold.

Indeed Linnik's theorem supplies an absolute constant C such that
the construction using four primes in [R,2R] always has output
bit length ≤C log R for large R. Pool the lists for every integer
length up to this upper bound before choosing the classes. For each
polynomial choose a nonzero coefficient. Their product has only
poly(log R) bits in absolute value, so four distinct primes in the
band can be chosen avoiding it. Over each selected field the total
number of roots of all polynomials in the pool is bounded by their
total degree, poly(log R). Choose x_i distinct from 0 and 1 outside
the root union in the i-th field. Prescribe p residues
(1,−1,x₃,−x₄) and q residues (x₁,−x₂,1,−1). The resulting primes
are distinct, and n has residue x_i in each field. No listed
polynomial value vanishes there. Each signed local order therefore
retains its designated prime ≥R≥n^c for an absolute c>0. Products
of tested exponents preserve this avoidance. With uniform group
samples the preceding small-order collision probability bound still
applies; polynomially many trials cannot restore constant success
through that collision mechanism.

For the twelve polynomials X^k±1, 1≤k≤6, the four primes
31,37,41,43 and residues x_i=3,2,2,2 yield p=3192877,
q=2790871. Direct modular evaluation confirmed all 48 nonvanishings.
The argument allows adaptive selection *within* the length-indexed
list. It does not cover polynomial coefficients computed from the
value of n, other local groups, or arbitrary modular queries.

### A cubic local group outside the signed first-stage criterion

For n=pq with distinct odd primes sample a uniform monic cubic f
mod n and a uniform element z of A=(Z/nZ)[X]/f. Compute N(z)
as the determinant of multiplication by z. A proper gcd of N(z)
with n gives a split; retry if n divides N(z). Otherwise form
g=z³/N(z) and test the gcd with n of all coefficients of g^E−1,
where E=(n³−1)lcm(1,...,B).

At a prime p where f is irreducible the unit group has order
(p−1)(p²+p+1), and the field norm raises z to p²+p+1. The map
g=z^(-(p−1)(p+2)) thus has cyclic image of order
H(p)=(p²+p+1)/gcd(3,p+2). A uniform nonzero z gives a uniform
image. The probability of irreducibility is (1−p^(-2))/3 by
counting degree-three elements up to Frobenius conjugacy. The
local cubics and elements at p and q are independent by CRT.
Hence both cubics are irreducible and both z components nonzero
with probability bounded below by an absolute constant. The algorithm
does not need to recognize irreducibility at either hidden prime.

If exactly one of H(p),H(q) divides E, then conditional on that
event of irreducibility and nonzero sampling, one component is
always killed and the other survives with probability at least
1/2. Thus fresh trials yield expected polylogarithmic factoring
for this class when B is polylogarithmic. This changes the local
groups used in the preceding polynomial-exponent obstruction.

At n=61597=31·1987 and B=3 none of the four signed first-stage
pairs has a fully annihilated local group. The cubic image orders
are 331 and 1316719, and only 331 divides E. The cubic
X³+X+12 is irreducible at both primes (no local roots). For z=X
the quotient g is 1+X/12. Direct powering modulo n and f returned
g^E=7875+51243X+54994X²; the coefficient gcd of g^E−1 with n
is 31. This example compares the first-stage guarantees only.

### Simultaneous obstruction for the degree-d norm-quotient extension

Consider precisely the extension using g=z^d/Norm(z), a uniform monic
degree-d algebra, and E=(n^d−1)lcm(1,...,B). There are infinitely
many odd semiprimes on which every degree 2≤d≤c sqrt(log n) has a
large surviving local prime at both factors for every squarefree
splitting pattern. Here c>0 is an absolute constant and B may be
any fixed polynomial in log n.

Construction: take D=floor(log R) and T=lcm(1,...,D)=exp(O(D)).
Choose a prime t between 2R and 4R. Use 2D distinct residues between
ceil(t/2) and t−1 modulo t, also requiring residue 1 modulo T.
Linnik's theorem gives 2D distinct auxiliary primes r>R of size
R^O(1), each 1 modulo every d≤D. Designate one r for each local
side and each index 1,...,D. At the auxiliary prime for side p,
index k, prescribe p of multiplicative order k and prescribe pq
to be a primitive root modulo r (this fixes q). Reverse roles for
side q. CRT and Linnik again give p,q with log(pq)=O(D log R).
They are distinct: at the index-one modulus the designated prime
is 1 and the other is a primitive root. Both exceed R by their
own index-one conditions. Thus D≥c sqrt(log n), R grows faster
than any fixed polynomial in log n, and no auxiliary r divides
n^d−1 for d≤D.

For an irreducible degree-d algebra at p, the index-d auxiliary r
divides p^d−1 but not p−1. Take an element z of order r. Its norm
lies in the scalar group and hence equals 1; z^d retains order r.
For a squarefree reducible algebra use the index-one r dividing
p−1. In one field component take z to be a scalar of order r and
in every other component set z=1. In the quotient z^d/Norm(z)
the ratio between these scalar coordinates still has order r.
Thus the image order is divisible by r in either case, independently
at both primes. A uniform unit maps uniformly to the image. Since
r does not divide E, the probability of the identity at either
local side is at most 2/R by the union bound. Uniform monic f
has a repeated factor with probability 1/p locally for d≥2;
sampling a nonunit z costs at most d/p locally for squarefree f.
These exceptional events also have vanishing probability over
poly(log n) fresh trials. This limits the terminal coefficient-gcd
test of g^E−1 for the specified exponent. Other ways to exploit
the algebra, and larger extension degrees, are outside the proof.

### Partial component hits and eigenvalue collisions

The scalar gcd of Norm(g^E−1) with n can succeed when the gcd of
all coefficients is 1. For n=143, f=X(X−1)(X−2), choose z with
values (2,3,5) at the three roots (interpolation is valid since 2
is a unit). At E=n³−1 the components of g^E−1 are (22,41,48)
modulo n. Their common gcd with n is 1 and their product gives 11.

On the preceding CRT family, however, even this partial-hit test
and the discriminant of the characteristic polynomial of multiplication
by g^E have combined probability O(d²/R) of vanishing locally,
conditioned on squarefree f and a uniform unit. To see this at p,
split the local algebra into fields, whose embeddings label the d
eigenvalues. The units are uniform on the product. For one field
component only, the preceding index-d argument bounds the chance
of identity by 1/R. For multiple components, test an eigenvalue
in one component while multiplying z in another component of degree
s by a scalar of order r>R dividing p−1 (the index-one auxiliary).
The eigenvalue of g^E changes by the scalar to power −sE, coprime
to r. Uniformity on each such orbit bounds the hit probability by
1/r. The same argument works for a prescribed target scalar.

For equality of two eigenvalues the common norm denominator cancels.
Between distinct field components, scale just one by the order-r
scalar; the ratio changes by its dE-th power. Within one field
of degree k, the two embeddings differ by a Frobenius power i
with 1≤i<k. Equality requires z_component^(dE(p^i−1))=1, after
applying an inverse Frobenius if necessary. The index-k auxiliary
prime divides p^k−1, divides none of d,E,p^i−1, and therefore gives
probability ≤1/R. Sum over d eigenvalues and their unordered pairs,
then repeat at q. Nonunits and nonsquarefree sampling have the same
exception bounds as above. This controls the norm and characteristic-
polynomial discriminant tests only; other coefficient identities
can behave differently.

### Scalar rescaling bounds a Frobenius-style norm test

Take distinct odd primes p<q and n=pq. In any fixed commutative
unital algebra free of rank d over Z/nZ, choose a uniform unit z
and a fixed base-linear algebra automorphism σ. Extract a candidate
factor by gcd(n,Norm(z^n−σ(z))). This includes a cyclotomic
automorphism chosen using n when it is defined. Put h=gcd(p−1,q−1),
x=h/(p−1), y=h/(q−1). The probability of a proper factor is at
most d(x+y), even for nonreduced algebras.

Indeed the norm is a multiplication determinant. Modulo p it
vanishes precisely if the difference belongs to one of the at most
d maximal ideals of the local algebra. Fix such an ideal and a unit
z. On the scalar orbit cz for c in F_p^*, the equality in its
residue field requires c^(n−1)=σ(z)/z^n in that field. Both sides
used in the ratio are units. There are at most gcd(n−1,p−1)=h
solutions among p−1 possible c. Uniform units are invariant under
this scaling. Union bounds over the ideals and the two primes prove
the estimate without requiring the automorphism to preserve each
ideal.

For comparison, a uniform scalar unit a extracts a proper factor
by gcd(a^(n−1)−1,n) with probability x(1−y)+y(1−x)≥x/2,
since y≤1/2 and y≤x. Thus the algebra norm test's per-trial
probability is at most 4d times the scalar test's. With d polynomial
in log n, inverse-polynomial success of the former implies the same
for the latter on the same inputs. The bound assumes σ and the algebra
are fixed before sampling z; norm-normalized or deliberately
structured samples need a different argument.

### An infinite separation from uniform sampling in the Gaussian algebra

Take any prime p=3 mod 4. Require a prime q to satisfy q=1 mod
4(p+1) and q=−1 mod (p−1). The moduli have gcd 2, so CRT gives
a reduced class; Linnik gives such a q of size p^O(1). Also
q≥4(p+1)+1. Thus n=pq gives an infinite family with p≥n^c
for some absolute c>0, and gcd(p−1,q−1)=2.

In the same rank-two algebra (Z/nZ)[i], i²=−1, use the same
automorphism (conjugation) and the same scalar output
gcd(n,Norm(z^n−conj(z))). With z a uniform unit, the preceding
scalar-rescaling bound gives success at most
2(2/(p−1)+2/(q−1))=O(1/p). Instead take a uniform unit w and
feed z=w/conj(w) into that test. Its local components are independent
uniform norm-one elements, of cyclic orders p+1 at the inert prime
and q−1 at the split prime. Now the test matches locally exactly
when z^(n+1)=1. At the split prime a hit in one component is a hit
in both, since the two coordinates are mutual inverses. Because
q−1 is a multiple of p+1, gcd(n+1,p+1)=p+1 and
gcd(n+1,q−1)=p+1. Thus conditional on a sampled unit, the exact
probability that the output is p is 1−(p+1)/(q−1)≥3/4.
Nonunit sampling is handled by the norm gcd as above. All steps
use the fixed Gaussian algebra and n alone, with no smoothness
bound or sign search. This is the mixed-sign mechanism specialized
to give a matched comparison between two input distributions.

### Structured samplers with nonzero scalar weight

Let n=pq with distinct odd primes. In the preceding automorphism
test suppose the input unit w has a distribution invariant under
scalar units (uniform suffices) and the map producing z satisfies
F(cw)=c^t F(w) for scalar units c and a known nonzero integer t.
Assume F outputs units and the algebra, F and automorphism are fixed
before w is sampled. Put E=|t|(n−1). For r=p,q set
x_r=gcd(E,r−1)/(r−1). A local component hit after rescaling requires
c^(t(n−1)) to attain one fixed nonzero residue-field value. The
same orbit argument bounds the probability of an algebra norm split
by d(x_p+x_q) for rank d.

A scalar trial with exponent E and its repeated-squaring chain
splits with probability at least max(x_p,x_q)/2. Indeed if either
fraction is below one then the smaller is at most 1/2, so the final
gcd alone works with at least the stated probability. If both equal
one, the unequal local 2-order valuations suffice as before. Thus
the rank-d test gains at most a factor 4d over this scalar trial,
even if |t| itself is large. Efficient comparison requires E to
have polynomial bit length and be available to the algorithm.
This applies in particular to w^a conj(w)^b of nonzero weight a+b.
At weight zero that monomial is a power of w/conj(w); scalar
rescaling fixes it, so the hypothesis yielding a useful exponent
disappears. The infinite Gaussian separation above is precisely
of this latter form.

### Bounded-degree reparameterizations after norm normalization

Fix the quadratic algebra, a unit discriminant at both odd primes,
and the tested positive exponent E. A uniform norm-one input at
r=p,q belongs to a cyclic group of size m_r=r−(D|r). Its local
identity probability after powering is x_r=gcd(E,m_r)/m_r.
Now pass the input through a rational self-map of the norm-one
torus, defined on its sampled points, whose local reduction at each
prime is nonconstant and has degree at most δ on the projective
completion. Over an algebraic closure this completion is a projective
line. Each target point has at most δ distinct preimages, by the
degree bound on the numerator of the fiber equation (or by using
homogeneous coordinates at infinity). Thus the new local probability
is at most δ x_r. This counts the final identity test per original
uniform draw; exceptional denominator gcds are separate tests.

For the terminal proper-factor event, the unmodified probability is
P=x_p(1−x_q)+x_q(1−x_p) by CRT. If both x_r=1, neither method
can split at the terminal test. Otherwise min(x_p,x_q)≤1/2 and
P≥max(x_p,x_q)/2. The reparameterized success probability is
therefore at most 4δP by the union bound. In particular a map of
polylogarithmic degree cannot promote negligible terminal success
to inverse-polynomial success. This is a comparison after passing
to norm-one inputs; the initial normalization from all algebra units
is not a self-map of this one-dimensional group.

The local degree qualifier matters for succinct constructions. The
map u↦u^k gives exactly gcd(kE,m_r)/m_r, so this particular
high-degree map is already captured by substituting kE for E.
Compact high-degree nonmonomial maps and maps degenerating to a
constant at just one prime are outside the stated comparison.

### Extracting a formal one-prime degeneration from a succinct circuit

Let n=pq with distinct primes, and suppose division-free arithmetic
circuits over Z/nZ describe P,Q of degree at most Δ. For a rational
map constant modulo p and nonconstant modulo q (with both local
denominators nonzero polynomials), the two-variable circuit
H(X,Y)=P(X)Q(Y)−P(Y)Q(X) is identically zero modulo p and nonzero
modulo q. This is a formal distinction: agreement on every base-field
point does not by itself imply a polynomial identity.

Choose k=max(2,ceil(log₂(8 max(1,Δ)))). Sample a uniform monic f of
degree k over Z/nZ and evaluate H at two independent uniform elements
of (Z/nZ)[T]/f. Take the gcd of n with all k output coefficients.
They are always zero modulo p. With probability at least 1/(2k)
the reduction of f modulo q is irreducible: the elements of degree
less than k in F_(q^k) are in proper subfields, with union bounded
by sum_{i=1}^{floor(k/2)} q^i≤q^k/2 for k≥2. Each degree-k irreducible
accounts for k elements. Conditional on irreducibility, Schwartz–Zippel
bounds the zero probability by 2Δ/q^k≤1/4. Thus a single trial
gives p with probability at least 3/(8k). No irreducibility test
or factorization is needed. With s arithmetic gates on degree-one
inputs one may use Δ≤2^s, so expected bit complexity is polynomial
in s and log n even if the expanded coefficient list is exponential.

As a diagnostic, at n=35 the polynomial 5(X^7−X) vanishes at every
scalar modulo n. In the quotient by T²+1 its value at T is −10T,
whose coefficient gcd is 5. The known factor in this toy polynomial
makes it a check on extraction only. A uniform construction of
suitable circuits from unfactored n remains open.
