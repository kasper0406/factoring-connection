# Factoring from N alone and evaluating the weighted disk in one pass

The entry point takes a balanced product of two distinct odd primes. It constructs one public principal form, evaluates mixed-orbit moments, and returns two verified factors. Both available evaluators take time polynomial in the *value* of N in the worst-case bounds. A polylogarithmic evaluation of the same moment on the public cutoffs would instead give polylogarithmic factoring for this input family. The existence of the reduction supplies no algorithm with that cost.

## Public branch construction

Let N=pq, p<q≤2p, and set ℓ=1, c=floor(2√N), Δ=4N−c². When Δ≤4 decrease c by 1 and recompute Δ. If gcd(c,N)>1 it is already a proper divisor since 0<c<N. Otherwise Δ is coprime to N. Its negative is an order discriminant and the vector α=((c−b)/2,1), b=Δ mod 2, is primitive of norm N in the principal form Q. The two public roots are ±c mod N.

The weighted-moment proof only needs this primitive norm-N vector, splitting at p and q, units ±1, and its first search cap A=floor(√(Δ/3)) below p. It does not require ℓ prime or coverage of a grid of weights. Before adjustment Δ<4√N. If adjusted, Δ≤4+2c_old−1<4√N+4. Since √N≤√2 p, 4√N+4<3p² for every odd p≥3. Thus A<p, including the small cases.

The earlier prime-ℓ assumption still belongs to the separate canonical-grid reduction. Here ℓ=1 is justified directly by the moment proof. No auxiliary prime search, modular square root, class lift, supplied norm index, or input factor is needed.

## Why a moment gives a divisor

Write w=2x+by and L=(cy−ℓw)/2. The moment W(B) sums L² over points with 0<Q≤NB and L>0, modulo N. Subtract the known principal contribution and multiply by c⁻² mod N. The resulting S(B) sums one y² per mixed conjugation/sign orbit. At the first mixed shell a, S(a) is the ordinary integer y² when the shortest direction is unique. A reduced-basis identity recovers the two y² values in the tied case. The positivity and non-wrapping argument locates this shell in O(log N) radial queries. The details, including ties and endpoint behavior, are in [the moment reduction](../ellipse_row_batches/README.md#one-global-second-moment-suffices-to-factor).

For a true candidate, w²=4Na−Δy² yields a mixed point and gcd((cy−ℓw)/2,N) returns a proper factor. All divisions and square tests are verified. Once a mixed point or its y coordinate is known, extraction has polynomial bit complexity.

For example take N=11009. The public branch has c=209, Δ=355 and b=1. At B=1, the raw moment is 9952 mod N. Subtract 209² and divide by 209² modulo N to obtain S=64. Then y=8, w²=44036−355·64=21316=146². The row is L=(209·8−146)/2=763 and gcd(763,11009)=109. Division gives the other factor 101.

## Aggregating all norm indices

A disk through B is equivalently the positive-quadrant matches

    w² ≡ −Δy² (mod 4N),       w²+Δy² ≤ 4NB.

The congruence forces the required parity w≡by (mod 2), and j=(w²+Δy²)/(4N) is an integer. Both coordinate axes are empty away from the origin in this range: a zero coordinate would force N to divide the other, contradicting X,Y<N for B≤A<p and q≥5. Each positive pair accounts for four signed vectors and one conjugation/sign orbit. The two signed directions with positive-or-zero representative L contribute, together,

    ((cy+ℓw)² + (cy−ℓw)²)/4 ≡ c² y² (mod N).

Accordingly a **scalar count and a sum of y²** suffice. They also give the unweighted signed count modulo 16. In the decoder range the principal orbit on shell k² has y=kℓ. Subtract ℓ²·sum(k², 1≤k≤√B) to obtain S(B) directly.

The join indexes y≤Y=floor(2√(NB/Δ)) by the residue −Δy² mod 4N. It probes with w² for w≤X=floor(2√(NB)). Each bucket is sorted by Δy². As w increases, a suffix beyond the ellipse boundary can be removed permanently. Accepted matches are grouped by j, producing a sparse list of prefix moments for every cutoff through B in one pass. The implementation builds geometrically increasing capacities on demand and reuses completed prefixes during binary search. Its arithmetic never factors N or Δ.

With ordinary hash-table operations the cost of a build is Õ(X+Y+h), where h counts accepted matches. For each y the quadratic equation in w has at most four roots modulo squarefree N, including nonunits. Since X<2N in this regime, there are at most eight candidate w per y, and h=O(Y). Memory is O(Y) stored integers up to logarithmic factors. There is a deterministic alternative to the hash-table complexity assumption: sort both residue lists and merge equal keys, at the cost of storing X+Y entries. Importantly, the implemented join always performs **X probes**, even for an empty mixed disk. At B=1, X=floor(2√N). At the public cap B=O(N^(1/4)), the join upper bound is Õ(N^(5/8)). Removing the shell loop has moved the cost to the long coordinate list. The logarithmic number of decoder calls does not bound the work inside a call.

The shell implementation uses the earlier geometric bound. Its conservative bound summed to B is Õ(B + min(√(N/Δ) B^(3/2), (N/√Δ)^(1/3) B^(4/3))). Neither implementation achieves the sought speedup.

## Why a polynomial product loses the signal

One tempting replacement for the long-side scan is P_X(t)=∏_{w=1}^X(t−w²), evaluated at t=−Δy² mod 4N. For this full rectangle, P_X(t)=0 mod 4N **for every y**, even when t is absent from the w list. Indeed the public root c/ℓ supplies a square at each prime p,q. A residue root up to sign at either prime has a representative in 1..q (allow zero via the prime itself), and X≥floor(2√N)≥q for balanced inputs. Thus each prime divides the product, possibly at different w. Modulo 4 the parity can be satisfied by some w as well. Local vanishing loses the requirement of a single common w in the interval.

At B=1 we computed the actual products for 4,071 targets across the 194 public branches. All were zero modulo 4N, including all 3,772 targets without any common w in the *rectangle*, before ellipse clipping. This rules out that product-vanishing test. It does not rule out aligned, blocked or differently weighted aggregate tests.

The harmonic coefficient formulation in the [row note](../ellipse_row_batches/README.md) remains exact. Bruin's fast Hecke-operator construction for arbitrary composite indices uses their prime factorizations, and the level dependence also matters for varying Δ. Substituting that routine for an unfactored index Nj would assume the missing information. [Bruin, §5.3–5.4](https://pub.math.leidenuniv.nl/~bruinpj/coefficients.pdf)

## Reproduction and measured costs

From the repository root run:

    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_aggregate/factor.py 11009
    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_aggregate/factor.py 11009 --method join
    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_aggregate/benchmark.py --output results-rerun.json

The input promise matches the mathematical reduction. The code does not certify the promise with a factoring or primality oracle. The CLI checks every returned factor by exact division.

* **Complete small run.** Both evaluators recovered the correct factors on every one of the 560 balanced odd-prime pairs with q≤241. Nine required the discriminant adjustment. None used the gcd shortcut. On the 194 public branches from the prior experiment, the join matched independent shell-based counts and moments at all 1,586 disk cutoffs. Those 194 comprise the two preserved collision cases and a replacement sample of 192.
* **End-to-end timing.** There are three seeded inputs in each scaling tier, chosen without filtering for easy shell minima. Only N reaches each timed factor function. The table includes branch construction, all queries and verification. Timings are single runs per ellipse engine, and medians of five runs for simple odd trial division per input. The table reports the median of the three inputs in each tier. The prime sieve used to construct fixtures is outside all timed algorithms.

| Actual bit lengths | Geometric shells | All-shell join | Odd trial division |
|---|---:|---:|---:|
| 19–20 | 0.00046 s | 0.00114 s | 0.0000075 s |
| 23–24 | 0.00435 s | 0.0116 s | 0.000039 s |
| 27–28 | 0.0170 s | 0.0551 s | 0.000181 s |
| 31–32 | 0.105 s | 0.304 s | 0.000644 s |
| 35–36 | 0.504 s | 2.430 s | 0.002714 s |

On the 36-bit input N=38,942,974,393, both engines found 153,247 × 254,119. There were 19 scalar moment queries. Geometric shells used 1,745,681 square tests, and the join used 20,607,475 long-list probes including geometric rebuilds. Previous 49-bit shell timings used a supplied shell and cannot serve as timings for the complete factorer.

**Open question.** Compute these prefix moments while preserving which two local roots belong to the same short integer representative, without a polynomial-size coordinate scan. The decoder, branch generation and factor verification already run in polynomial bit time apart from their oracle calls. No polylogarithmic evaluator, complexity lower bound for arbitrary evaluators, or new competitive factoring exponent follows from these experiments.
