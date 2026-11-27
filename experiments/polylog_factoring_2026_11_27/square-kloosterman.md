# Factoring from square modulus Kloosterman approximation

For an odd input with exactly two distinct prime divisors, constant-precision approximation of an ordinary Kloosterman sum modulo $N^2$ suffices for randomized polynomial-time factoring. Only its upper base-$N$ digit varies. This gives an alternative to the Salié target with no Jacobi twist. [Salié approximation](salie-quarter.md)

Let
$$K_M(1,b)=\sum_{x\in(\mathbb Z/M)^\times} e_M(x+b/x),\qquad e_M(t)=e^{2\pi it/M}.$$
For $M=N^2$ and unit $c$, write $x=x_0+Ny$. Modulo $N^2$,
$$x+c^2/x=x_0+c^2/x_0+Ny(1-c^2/x_0^2).$$
Summing over $y\bmod N$ gives zero unless $x_0^2=c^2\pmod N$, in which case it gives $N$ times the first phase. Each stationary $x_0$ has a unique lift $r$ with $r^2=c^2\pmod{N^2}$, by the unit derivative $2x_0$. At the lifted value the phase is $2r$. Consequently
$$\frac{K_{N^2}(1,c^2)}N=\sum_{r^2=c^2\bmod N^2}e_{N^2}(2r).$$

## Vary just one digit

Let $u$ be a mixed root modulo $N$ and $U$ its lift modulo $N^2$. Choose $c=1+Nj$ with $j\bmod N$. The public root pair $\pm c$ contributes $2\cos(4\pi c/N^2)$. The remaining pair contributes
$$2\cos(4\pi U/N^2+4\pi uj/N).$$
Thus the unknown component is a sinusoid in $j$ with frequency $2u\bmod N$, even though all the scales have the same lower digit.

Suppose a public algorithm approximates $\operatorname{Re}K_{N^2}(1,(1+Nj)^2)/N$ to RMS error at most $1/2$ over uniform $j$. Subtract the public term, divide by two and clip to $[-1,1]$. The resulting signal has RMS error at most $1/4$ from a unit-amplitude cosine. Its correlation with that cosine is at least $1/4$ by Cauchy–Schwarz and the exact mean square $1/2$ of the cosine. This forces a Fourier coefficient at $2u$ of magnitude at least $1/4$ despite the unknown phase. The public sampled Fourier search finds the frequency. Testing its difference from 2 by gcd gives a proper divisor. The sum identity and lifting also apply when $N$ has repeated factors.

## What the block cancellation saves

The lift sum reduces $N^2$ terms to a scan over at most $N$ residue blocks, with a polynomial-bit-time stationarity test per block. Scanning only for the unknown stationary residues is exactly a search for the mixed roots modulo $N$. Once the two public stationary residues have been removed, a uniform block hits either unknown block with probability $2/N$. Independent uniform block sampling therefore does not provide the needed precision efficiently. An evaluator that aggregates or locates stationary blocks without this scan remains open.

## Checks

Direct Kloosterman enumeration verified the identity and the two upper-digit frequencies on all 179 upper digits across five small moduli, including two repeated-factor inputs. The diagnostic also verified a proper gcd from both recovered frequencies. Run `python check_square_kloosterman.py`.

## One fixed sum at higher precision

For odd n with exactly two maximal prime powers, no parameter variation
is needed if K_(n²)(1,1)/n is available within absolute error n^(−4).
Put M=n². The four square roots of one mod M are ±1 and ±U for one
mixed root U. The stationary-phase identity above gives

    K_M(1,1)/n = 2 cos(4π/M) + 2 cos(4πU/M).

Subtract the known cosine pair and divide by two. The result estimates
cos(2πj/M), where j is the folded residue of 2U, within 1/(2M²).
Clip to [−1,1], take arccos, and round in units of 2π/M. Moving
π/M from a true angle within [0,π] changes its cosine by at least
1−cos(π/M)≥2/M²; clipping can only improve the error. Thus the
rounded index is exactly ±2U mod M. Halving modulo M and taking the
gcd of the resulting root minus one with n yields a proper divisor.
This works with repeated exponents and costs polynomially many bit
operations given the sample. For the unnormalized sum the sufficient
absolute error is n^(−3).

Direct summation over units, without supplied factors, verified the
decoder on 15,21,35,45,63,77,99,143 with both endpoint perturbations.
The hard step is evaluation at binary modulus in polynomial bit time;
the single-sum decoder does not accelerate that step.
[Decoder](cubic_masked_trace.py)

## Folded-angle square: inverse-square precision suffices

Keep M=n² and let j be the folded index ±2U from the preceding
section, 0≤j≤M/2. Put w=min(2j,M−2j). Then w≡±4U mod M,
0≤w≤M/2, and |cos(2πj/M)|=cos(πw/M). In particular
L=(w²−16)/M is an integer. A mixed root cannot give w≤4 since
M>16 and the only possible square there congruent to 16 is 16,
which corresponds to a public root.

For a normalized-sum estimate within B/M, subtract and halve as above,
take absolute value, and clip to [0,1], obtaining c within B/(2M)
of the true absolute cosine. Approximate L by
M(arccos(c)/π)²−16/M. On [0,1] the function (arccos(c)/π)² is
1/π-Lipschitz: its derivative magnitude is 2θ/(π² sin θ)≤1/π
for θ in [0,π/2], with a finite endpoint limit. Thus the estimate
is within B/(2π) of the true integer, even at either cosine extremum.
For B=1 nearest-integer rounding suffices. More generally enumerate
O(1+B) nearby integers, testing whether 16+ML is an integer square.
At the correct L, gcd(w−4,n) is a proper divisor since n is odd
and w represents a mixed choice of the signs of four. All candidates
are checked by gcd. Polynomial B in the bit length preserves polynomial
decoding time and costs O(log n) bits of working precision.

This uses only one K_(n²)(1,1) value, now with sufficient absolute
error B/n in the unnormalized sum. The implementation passed 32
endpoint checks using literal sums on the eight preceding inputs with
B=1 and 32. It supplies no faster way of computing the sum.

## One random upper digit and constant absolute error in the raw sum

Choose c=1+nk with k uniform modulo n. If a raw K_(n²)(1,c²)
estimate has absolute error ≤1/16, subtract the public cosine as
above to estimate cos(2πj/M) within 1/(32n), with folded j=±2cU.
Keep estimates with absolute value ≤3/4. On acceptance the entire
arccos error interval lies in [−0.76,0.76] for n≥17, so its derivative
has magnitude below 2. The rounded angle index J then satisfies
|j−J|≤ceil(n/64)+1 = X. The true angle traverses an equally spaced
n-point grid, since 2U is a unit. At least n/3−2 grid points have
absolute cosine ≤1/2, ensuring constant acceptance probability.

The small correction x=j−J solves the monic quadratic
f(x)=(J+x)²−4c² mod M. Use the lattice of coefficients after
substitution x→Xx in x^i f(x)^l M^(m−l), for i=0,1 and 0≤l≤m.
Its dimension D=2m+2 and determinant Δ=M^(m(m+1))X^(D(D−1)/2)
follow from triangular leading terms. Each polynomial vanishes modulo
M^m at the correction. Choose m until
2^(D(D−1)) Δ^4 D^(2D)<M^(4mD). The LLL length bound times √D
then makes the first reduced polynomial vanish over the integers.
This loop needs O(log M) iterations: after taking the 4D-th root
and dividing by M^m, the bound is
√D M^(1/4) (√2 X/n)^(m+1/2), and X/n≤1/8 for n≥17.
Extract integer roots, check the congruence and test gcd(J+x−2c,n).
Polynomial reductions of coefficients of the initial monic f preserve
all these conditions. Small bases can be handled by direct division.

This achieves expected O(1) oracle queries at constant absolute error
in the unnormalized sum on two-component inputs. The parameter is
randomized and the lattice dimension can grow with input length.
A controlled test with endpoint errors ±1/16 on eight random lifts
per input recovered splits on 7,5,5,4 trials at n=35,99,10403,
137446031423, respectively. Other trials were rejected by the angle
check. These tests used the stationary formula with supplied components.

## Two scales tolerate square-root absolute error in the raw sums

Let H be the largest power of two ≤√n. Sample c=1+nk as before
and query the raw sums at square parameters c² and (Hc)² mod M.
Each raw error may now be H/64. Subtracting the public term leaves
cosine error ≤H/(128n). Accept when both observed cosine magnitudes
are ≤0.95. Then each arccos divided by 2π estimates its folded angle
in turns to within H/(128n): the arccos derivative is below 4 on
the intervening segment. Both true angles individually traverse uniform
n-point grids (H is a unit); the union bound on the arcs where absolute
cosine exceeds 0.94 gives constant joint acceptance for large n,
without assuming independence. Small bases can be trial divided.

Let A,B be the two estimates and x the first true folded angle.
For one sign σ, Hx is an integer plus the true signed fine angle.
Rounding HA−σB identifies this integer because
(H+1)H/(128n)<1/2. For that sign, (round(HA−σB)+σB)/H differs
from x by at most 1/(128n). Multiply by M and round: one of the
two centers is within the preceding certified quadratic lattice radius
of the index. Test both signs and verify any divisor. Each trial takes
two queries, with O(1) expected trials. The two square parameters need
not have the same lower base-n digit.

Endpoint-perturbed factor-supplied checks over twelve trials per input
gave 9,8,6,7 verified splits for n=35,99,10403,137446031423
respectively. At the largest input each raw error magnitude was 4096.

## Adjustable radix with a shared reference

For a power-of-two radix R≥2, a single trial needs at most
1+2 ceil(log_R(n²)) raw estimates, each within n/(128R).
Choose c=1+nk uniformly as above. A normalized, public-term-subtracted
sample at multiplier m estimates cos(mθ) to δ=1/(256R), where
θ=4πcU/n². Accept the reference if its estimated absolute cosine is
≤3/4. Its sine magnitude is then bounded away from zero. Orient θ so
the reference sine is positive (this chooses the global sign of U).

At m=R^i also sample m+1; if that multiplier produces a proper gcd
with n return the split, and if it vanishes modulo n use m−1, a unit.
With d=±1 the signed sine estimate is
d(C_m C_1−C_(m+d))/sqrt(1−C_1²). Here C denotes measured cosines.
The numerator error is ≤4δ, denominator error ≤2δ, and denominator
≥√7/4. Thus the reconstructed complex point differs from exp(imθ)
by <11δ. Its argument error in turns is <4δ (the point lies inside
a disk of radius 11δ about the true point). The initial folded reference
has angle error <4δ as well. If the previous angle estimate has error
≤4δ/(m/R), multiplying by m incurs error ≤4Rδ. The new angle's
integer winding is therefore unambiguous since 4(R+1)δ<1/2. Dividing
by m gives error ≤4δ/m, with no branching. At m≥n², round the index
and verify the divisor via gcd. The reference succeeds with constant
probability as in the single-scale trial.

Choosing R comparable to log n gives O(log n/log log n) queries with
raw error Θ(n/log n); fixed R gives O(log n) queries and Θ(n) raw
error. The implementation uses no LLL. Direct finite Kloosterman sums
at bases 35,99,143 were tested with endpoint perturbations and the
permitted early-gcd shortcut. A factor-supplied synthetic control at
47·59 also tested the full radix-8 winding reconstruction.
