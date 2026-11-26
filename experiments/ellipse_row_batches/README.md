# Batching the rows of the modular ellipse

26 November 2026. This experiment continues the ellipse work on [factoring-connection, base 9e70736](https://github.com/kasper0406/factoring-connection/blob/9e70736b4abaacff226c93e4fdddaf5e6fcad5f5/notes/shared_collision_nov26.md). It includes an exact public evaluator for both whole disks and clipped intervals, a bound for its running time, and a reduction from factoring to one weighted radial observable. The polylogarithmic evaluator and the conditional almost-all exponent remain open.

## Inputs and notation

Use the principal form of discriminant −Δ, with Δ>4:

    b = Δ mod 2,     Q(x,y) = x²+bxy+(Δ+b²)y²/4.
    4N = c²+Δℓ²,    gcd(N,ℓΔ)=1,    N=pq,    p<q≤2p.
    α = ((c−bℓ)/2, ℓ),    L(x,y) = det(α,(x,y)).

The eligible branches of the repository use positive c and ℓ. A shell has Q=Nj. The disk is the union of shells 1≤j≤B. The raw radial result counts both signs modulo 16. The clipped result counts points with L in an inclusive interval of positive integers modulo 8. Setting R to the raw count and A to the positive half count gives R=2A+2⌊√B⌋. This last term consists of the multiples of α. Assume B<min(p,q) when using the four-line decomposition. For the factoring reduction it suffices that floor(√(Δ/3))<p.

## An exact native batch over a row interval

Set w=2x+by. Each shell consists of the integer solutions of

    w²+Δy² = H = 4Nj,    w ≡ by (mod 2),    2L=cy−ℓw.

No modular square roots, divisors of N, or factorization of Δ enter either of the following evaluators.

**Project the interval.** From the two equations,

    (2Ny−cL)² = ℓ²(4N²j−ΔL²).

For l≤L≤h, intersect the interval with [1,floor(2N√(j/Δ))]. The functions cL and √(4N²j−ΔL²) have opposite monotonicities there. Their endpoint extrema give two enclosing integer y intervals, one for each sign. Merge them, test H−Δy² for squareness, then check the original L inequality exactly. Outward-rounded integer square roots preserve endpoint ties. This batches an arbitrarily positioned interval and needs at most O(1+√(Nj/Δ)) square tests per shell. For narrow intervals the projected range can be much smaller.

**Cover curved arcs with integer lines.** In the first quadrant write either coordinate as z=f(x)=√((H−Ax²)/C), with (A,C) equal to (Δ,1) or (1,Δ). Use each graph only where Ax²≤H/2. The two sectors cover the quadrant, including their common boundary. On this sector,

    |f''(x)| ≤ 3A/√(CH).

Take K=max(1,floor((√(CH)/(4A))^(1/3))). For K=1 directly check the coordinate rows in that sector. Otherwise start a horizontal block of K integer x values at x₀. A continued-fraction convergent u/v to the signed derivative has 1≤v≤K and |v f'(x₀)−u|≤9/(8K). The code computes a rational approximation to the derivative with error ≤1/(8K²) and takes its last convergent of denominator at most K. All operations are on integers. Taylor's formula gives, throughout this block,

    |v f(x)−u x − (v f(x₀)−u x₀)| ≤ 9/8 + 3/8.

Consequently every integer point belongs to one of five integer lines vz−ux=k centered on floor(vf(x₀)−ux₀). For each line let d=Av²+Cu². A solution requires H d−ACk² to be a nonnegative square, and then x=(−Cuk ± v√(Hd−ACk²))/d. Verify divisibility, block membership, and z≥0. Restore signs, parity and the requested L interval. Sector-boundary points are counted once.

The number of blocks per sector is O(1+(H/√(AC))^(1/3)). This follows by dividing its horizontal width O(√(H/A)) by K; if the unrounded K is smaller than 2 the direct scan satisfies the same bound. Each block uses a constant number of exact quadratic intersections and polylogarithmic bit work. With the public choice of the cheaper estimated method, a full shell or clipped shell costs

    Õ(1 + min{ √(Nj/Δ), (Nj/√Δ)^(1/3) }).

Here Õ suppresses factors polynomial in the bit length. A disk costs the sum over j≤B, at most Õ(B + min{√(N/Δ) B^(3/2), (N/√Δ)^(1/3) B^(4/3)}). The implementation compares the projected row count with the number of direct rows and line tests in the arc cover. This proxy is a scheduling choice, not a prediction of elapsed time. The arc branch enumerates all points on a shell before applying a clip. Related geometric background is [Cilleruelo and Córdoba, lattice points on short elliptic arcs](https://matematicas.uam.es/~antonio.cordoba/documentos/articulos/Latticepointsonellipses.pdf).

## What the product screen buys

For j≤floor(√(Δ/3)), the public positive L row on the exact shell exists only if j=k², and is k cℓ. Every mixed row shares exactly one prime with N. For a block of d rows starting after u, compute F_d(u)=∏_{i=1}^d(u+i) mod N. If this is a unit, the entire block has zero mixed contribution. All full blocks use a common product polynomial and multipoint evaluation modulo the composite N. The final partial block is multiplied directly. A surviving block falls back to the geometric interval projection; no divisor obtained from a gcd is returned or used to choose a mixed root. Polynomial block products are a classical construction in deterministic factoring. [Strassen background](https://doi.org/10.1090/mcom/3313)

For a prefix 1..T with T<N, exactly floor(T/p)+floor(T/q) indices are nonunits. A block of length at most p contains at most two of them. Thus an interval-by-interval unit screen must handle at least half that many surviving blocks. At full shell height T≈2N√(j/Δ), this is already of order √(Nj/Δ) on balanced inputs. If the block length exceeds p, every full block contains a nonunit. These statements concern this particular rejection certificate.

| First mixed shells, 194 branches summed | Blocks | Surviving blocks | Surviving blocks with no mixed point | Projected square tests |
|---|---:|---:|---:|---:|
| Block width floor(N^(1/4)) | 76,651 | 11,607 | 11,215 | 2,279 |
| Block width floor(√N) | 6,127 | 6,076 | 5,689 | 14,384 |
| Whole-interval projection, no screen | 194 intervals | — | — | 8,847 |

The first row also requires 76,478 polynomial values. A reduction in square tests alone overstates its benefit.

## One global second moment suffices to factor

Define a full radial observable using the signed vectors:

    W(B) = (1/2) Σ L(x,y)²  mod N,   0<Q(x,y)≤NB, N|Q(x,y).

Equivalently sum over L>0. The factor 2 is invertible mod N and rows with L=0 contribute zero. W uses a larger output modulus and a different weight from the original unweighted query.

Let a denote the first mixed norm after division by N. By reduction a≤A=floor(√(Δ/3)). Since A<Δ/4, at any B≤A the public positive rows are precisely k cℓ for 1≤k≤floor(√B). Subtract their moments and normalize:

    S(B) = c^(−2) [ W(B) − (cℓ)² Σ_{k≤√B} k² ]  mod N.

Consider the orbit consisting of a mixed vector, its conjugate, and both negatives on a shell j. It has four members: a mixed vector cannot equal its conjugate up to sign, since conjugation exchanges the two distinct mixed lines and they are disjoint at these cutoffs. Its two positive squared L values add to

    2Njℓ² + (c²−Δℓ²)y²/2  ≡ c² y²  (mod N).

Thus S(B) is the sum of y², once per mixed conjugate/sign orbit. It is zero before a. At B=a its least nonnegative residue is the actual positive integer sum, as follows.

* **One shortest pair in one mixed lattice.** There is one conjugate/sign orbit across the two mixed lattices. The moment is y², with 0<y²≤4Na/Δ<N.
* **Two shortest pairs.** Their vectors form a reduced basis of one mixed lattice, whose scaled form is [a,b₁,a]. Here 0<|b₁|<a; a tied square or hexagonal primitive form would have Δ=4 or 3. Three shortest pairs are possible only in the hexagonal case. These assertions also follow directly from |b₁|≤a≤d in a reduced form by testing vectors of norm a. Let y₁,y₂ be the physical y coordinates of the two basis vectors. If T has these vectors as columns and G₀ is the Gram matrix of Q, then TᵀG₀T=N G for the Gram matrix G of [a,b₁,a]. Taking the lower-right entry of TG⁻¹Tᵀ=N G₀⁻¹ gives

      a(y₁²+y₂²)−b₁y₁y₂=N.

  Consequently S=y₁²+y₂²≤N/(a−|b₁|/2)<N. Since b₁²=4a²−Δ is known, (y₁y₂)²=(aS−N)²/(4a²−Δ). The two squared coordinates are the integer roots of the quadratic with sum S and this product. All candidates are checked exactly.

To find a with W alone, scan B up to min(7,A). If a larger search is necessary, double from 8 with endpoints capped at A, then binary-search the first nonzero S. For a≥8 and a≤B<2a there are at most three orbits: distinct antipodal pairs with norm <2a occupy distinct nonzero parity classes in one mixed lattice. Each contributes at most 4NB/Δ, and their total is strictly less than 24Na/Δ≤8N/a≤N. It cannot vanish modulo N. Thus the search uses O(log N) weighted radial queries, including exceptional ties, without access to a radial count oracle.

Finally try the one-orbit square test and, when admissible, the tied quadratic. For a candidate positive y, compute w=√(4Na−Δy²) exactly, check the parity, and test gcd((cy−ℓw)/2,N). For a true mixed vector this is a proper divisor. This proves the reduction; the measured native evaluators still enumerate shells. No class lifts are needed for this observable.

The usual theta series provides another statement of the target. On Q=Nj the quadratic harmonic weight (2x+by)²−Δy² is congruent to −2Δy² modulo N. Its coefficient at Nj therefore carries the same physical-coordinate second moment after accounting for multiplicities and the known principal terms. Efficient composite-index coefficient extraction remains unresolved. [Weighted theta background](https://people.math.harvard.edu/~elkies/M272.19/oct16.pdf)

## Periodic row batches with output modulo 16

One can retain the small output modulus by changing the query to a residue-filtered count. Use the ordinary radial queries to get a. On the exact shell Q=Na, ask for the signed count with L≡r mod s, for all residues of small primes s. A shell response is also obtainable by differencing two disk responses. Subtract the public points in each bucket. Since there are at most eight signed mixed points, each bucket count modulo 16 is now exact. Its r²-weighted sum gives the signed mixed second moment modulo s.

Chinese-remainder the moments until the product of the auxiliary primes exceeds 8T², where T=floor(2N√(a/Δ)). The signed moment is then recovered as an ordinary integer. Divide by 2 and normalize by c² modulo N to recover S. Using all primes up to O(log N) suffices by the growth of the primorial, giving O(log² N) scalar bucket queries. The largest auxiliary prime in the 194-branch run was 29. These queries have an additional arithmetic-progression filter; no implementation of them from only the ordinary radial count is implied.

## Benchmarks and reproduction

The original late-stage 192-input sample remains unavailable. For continuity, public_cases.json contains exactly the two preserved collision cases and the separately generated 192 cases from the preceding ellipse experiment. It contains only N,ℓ,c,Δ. All 1,586 shells up to the public search caps matched an independent x-coordinate enumeration. On the first shells, 6,188 clipped projection checks included endpoints immediately before, at, and after occupied rows. Both full-disk and clipped-disk residue and moment queries were checked on every branch.

All 194 branches yielded a proper divisor via the global moment with at most 11 radial queries. The 192 single-orbit and 2 tied branches all passed. Ten additional deliberately selected tie controls also passed. The collision examples give S=1849=43² for N=11009 and S=1444=38² for N=28829. At N=8927, a=9 and S=1037 is not a square. The tied equation gives {y₁²,y₂²}={196,841} and the extraction returns 113.

Two series of four larger eligible inputs compare per-shell costs. One keeps Δ=20 and j=1. The other uses growing Δ and a supplied first mixed index selected with known fixture factors. These are shell timings, excluding the search for that index. On the 49-bit fixed-Δ input, projection used 7,798,922 square tests in about 2.20 s, while arcs used 731,400 line tests in about 0.52 s. On the 49-bit growing-Δ input at j=3,173, projection took about 0.10 s and arcs about 0.68 s. The default evaluator selects projection on the latter case. Timings are single-run comparisons; the checked operation counts are in scaling.json.

From this directory:

    uv run --no-project --with python-flint==0.9.0 python row_batches.py --branch 11009 11 204 20 --cutoff 5
    uv run --no-project --with python-flint==0.9.0 python row_batches.py --branch 11009 11 204 20 --cutoff 1 --rows 3924 3924
    uv run --no-project --with python-flint==0.9.0 python row_batches.py --branch 8927 7 141 323
    uv run --no-project --with python-flint==0.9.0 python row_batches.py --output results-rerun.json
    uv run --no-project --with python-flint==0.9.0 python scaling.py --output scaling-rerun.json

The first command gives the raw signed residue and positive second moment. The second clips to one determinant row. Omitting the cutoff runs the public moment-based factor extraction on a promised eligible branch. The technical question left by the experiments is how to aggregate over many norm indices or extract the harmonic coefficients without locating individual integer points.
