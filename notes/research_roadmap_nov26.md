# Research roadmap — 2026-11-26

Goal: classical factoring of general semiprimes, ideally polynomial in log N
or below L_N[1/3]. Do not conflate an improvement in N^c deterministic
complexity, an improved NFS constant, a new reduction, and a new general
factoring algorithm. The internally audited prior results and their
hypotheses remain in [STATUS](../STATUS.md).

## 1. Canonical index and public modular ellipse count

Read [shared_collision_nov26](shared_collision_nov26.md) first.
For a public eligible branch (N,ell,c), Delta=(4N-c^2)/ell^2 and
Q0 the principal form of discriminant -Delta, evaluate
    #{(x,y): 0<Q0(x,y)<=NB, Q0(x,y)=0 mod N} mod 16
without the factorization or mixed CRT root. Charge failed branches and all
preprocessing. The integer-count reduction recovers the squared hidden
class, up to inversion, by O(log N) queries at bounded cutoffs; class-square
lifting remains separate. The principal term must be counted exactly at
endpoints, not always replaced by floor(sqrt(B)). Fresh small validation is
in experiments/nov26_consolidation/check_identities.py.

**Priority for polylogarithmic factoring.** Park norm-level scans and class
catalogues. The primitive band now has an exact formula that sums across its
levels, but producing its inputs still costs a power of $N$ on growing
discriminants. Reopen this particular approach on evidence for an evaluator
that obtains the hidden character information directly from the public input.
A [one-bit test at fixed discriminant](../experiments/ellipse_disk_mod8/cubic59.md)
now has a quartic formulation and explicit controls. The canonical-index
coverage audit remains relevant to the separate conditional $N^{1/8}$ goal.
[Decision, counterexample and reopening test](../experiments/ellipse_disk_mod8/class_trace.md)

**Weighted route.** With ℓ=1 and c chosen next to 2√N, the squared-coordinate
moment gives a direct factor from N alone. Both the geometric evaluator and
a residue-list join run through verification. Their costs remain powers of N.
The immediate target is computing a sharp harmonic prefix at multiples of N
without shell enumeration or a length-√(NB) scan. A full-range product test
loses local root alignment. See [implemented factorer, measurements and obstruction](../experiments/ellipse_aggregate/README.md).
This route bypasses class lifting and the canonical-grid coverage hypothesis.

**Whole-disk parity route.** The first shell can be located and decoded using
only parity of positive-quadrant disk counts with at most two y-bit filters.
These are stronger queries than the original unmasked count. The native
Boolean decision diagram shares all j≤B in one carry relation and can eliminate
coordinates over F₂. Its construction is expensive even when the final result
is tiny. Pursue a different intermediate representation that cancels histories
before a common w² fixes them; measure peak state and query construction.
[Proof, implementation and controls](../experiments/ellipse_disk_mod8/README.md).
The [translation and shear audit](../experiments/ellipse_disk_mod8/translation.md)
shows that constructing a short exact linear cancellation of translates of the completed
mixed congruence already gives a factor. Ordinary unit shears lose every primitive
mixed point on the original cone. Splitting both coordinates by low bits shrinks the disk but creates Ω(√N) nonempty children with distinct translated kernels already at B=1 in the literal recursion. Deprioritize fixed translation identities and
decoder tuning. First test an exact recurrence at the required disk center, allowing
a change of form. Charge the boundary and all child calls. Stress-test it immediately
before and at the first mixed shell, then check that it supports the intended decoder.
The [centered projective recurrence](../experiments/ellipse_disk_mod8/centered.md)
now achieves that exact transformation. Its unmasked mod-16 decoder uses O(log N)
queries on changed forms, including nonprincipal classes. The native traversal
takes Θ(NB/√Δ) nodes in the growing principal-branch regime. A stronger next
test is to evaluate a whole group of sibling slopes without visiting them
individually. At B=1 the 35-bit and 36-bit controls each give 8,192 inequivalent
forms representing N locally everywhere. Only the former group contains a
mixed global representation. Account for building the group and testing its
residue on both fixtures.
The [bisector evaluator and coverage bound](../experiments/ellipse_disk_mod8/bisector.md)
now solve the opposite-outcome B=1 groups. Exact reflection turns a mixed
norm-N point into a short linear nonunit search, shared across all slopes.
Discriminant factorization and all product-tree work are charged. For the
ℓ=1 near-square family, the fraction with any mixed point through B is at most
P^{-1/2+o(1)} B^{3/2} for primes of size P. Prioritize a global aggregate
at growing nonsquare cutoffs; testing more nearby public points only at B=1
cannot give almost-all coverage by mixed witnesses.
The [class-prefix decision and measured list budget](../experiments/ellipse_disk_mod8/class_prefix.md)
provide a second way to handle nonsquare minima. Lift all possible squared classes
in a proposed interval, then batch linear products to decide whether a listed
form represents a prime divisor of N. The full-cap version factors without a
supplied index, but generating its current list costs at least one visit per
possible minimum. The immediate research question is whether public c can select
or test a block of classes cheaply. Measure list construction, hit rate, and
the sum of short coordinate ranges on held-out inputs. For the ℓ=1 integer
basis formulation, both unrestricted local solubility and continuous balanced
feasibility hold for every primitive candidate; a selector must combine integer
constraints with the archimedean bounds. A cheaper selector would address the
power-exponent track. The supplied-class verifier itself still has a power-of-N
bound, so a polylogarithmic result would require a further idea.
The [public middle-coefficient shortlist test](../experiments/ellipse_disk_mod8/class_selection.md)
found 14 of 96 fresh inputs with small rational weights, versus 13 for equally many
low indices and 96 for the full class list. For primes in one dyadic interval,
selecting only prime-class minima below $N^{1/4-\epsilon}$ on this branch
has vanishing coverage. Move the next evaluator test to a growing nonsquare
cutoff. At the 47-bit fixture the signed punctured-disk residues at 581 and
582 are 0 and 4. Test bulk aggregation with all preprocessing and boundary work
charged. A selector that reaches typical-sized minima remains a separate possibility.
The [rectangle evaluator](../experiments/ellipse_disk_mod8/rectangle_rays.md)
now reproduces that cutoff pair and factors the 47-bit input from clipped queries.
Scaling a primitive point couples multiple norm levels; its coefficient can vanish
before the final square test. At 582 the filtered parity evaluator still runs
16,546 of its 18,042 square tests on primitive candidates in nonsquare levels
above B/4, where the multiplier is one. Prioritize a joint test over these
progressions that can also avoid listing at least B seed residues.
The [upper-band experiment](../experiments/ellipse_disk_mod8/conductor_band.md)
shows why the 581/582 transition is hard to resolve by cancellation. Below
√Δ/2 each mixed lattice has at most one direction. Public conductor descent
produced 24 possible norm levels at 582; verifying them by targeted modular
roots was cheaper than selecting the prime class. A counting bound rules out a
polynomial conductor saving on almost all near-square ℓ=1 inputs.
The [class-character experiment](../experiments/ellipse_disk_mod8/class_trace.md)
then reproduced 578/579 with one trace of degree 155 and 280/281 with five
traces of total degree 59. Both computations still required the full class
catalogue and a norm-$N$ representation scan. All lower-order traces together
miss the relative sign on two actual inputs at the same discriminant.
For a new public readout, test the order-three bit at $\Delta=59$ before
building any growing class fields. The substitution of $N$ for the prime
Frobenius exponent gave the same gcd degree on 128 cases with both bit values.
The [fixed-discriminant follow-up](../experiments/ellipse_disk_mod8/cubic59.md)
connects the cubic residue symbol and the earlier norm-one group tests.
For $N=t^2+3t+17$ the quartic $X^4-2tX^2-8X+t-5$ has one root at bit 1
and zero or sixteen at bit 0. All 320 controls agree. The public cubic symbol
cancels, the tested pair of shifted norm-Jacobi symbols has vanishing advantage,
and substituting $N$ in the prime halving exponent missed all 221 unique roots.
The [level-59 modular-symbol evaluator](../experiments/ellipse_disk_mod8/manin59.md)
computes this parity from $N$ alone using a sum over every $b/N$. It agreed on
zero, one, and sixteen roots. The faster selection of seven endpoints lost the
bit on an explicit matched triple. The next evaluator test is a joint recurrence
for the parity of continued-fraction edge visits over all $b$ at the branch
denominator, with all intermediate states charged.
The [growing-discriminant test](../experiments/ellipse_disk_mod8/growing59.md)
shows exactly which part of a class all unramified cubic bits determine. At
$-4027$ two or three factor-backed oracle answers decoded 288 of 288 hidden
classes up to inversion. At $-1579$ all such answers collide at two different
minima. The direct nonnegative-coefficient near-square catalogue also misses
almost all balanced inputs asymptotically. A general reduction needs both
branch coverage and a way to distinguish classes lost modulo cubes.
For the independent variable-square reflection reduction, see [floor8](floor8_results.md).

Next tasks for the original unweighted route:
- Reconstruct and independently audit the canonical band-coverage proof
  (large sieve + Gauss-sum estimate for small Delta), including small-order
  exceptions, parity, nonfundamental discriminants and square-weight
  distinctness. Verify the conditional N^{1/8} exponent.
- Try a reciprocity/floor-sum or mod-16 lattice count with a provably shrinking
  argument; use the known two unmixed norm-N roots but do not supply a mixed
  root or secretly factor N. Validate against exact counts at multiple B.
- Determine whether class lifts can be enumerated in the needed time given
  public Delta, including the 2-primary ambiguity; verify each minimum
  by the square-discriminant test. Do not assume C has odd order.
- Compare runtime of the complete construction with D'' (deterministic,
  almost all balanced pairs), and report the failure fraction separately.

A fast ellipse counter alone is progress, but does not establish the
canonical-minimum evaluator or the conditional factoring exponent.

## 2. Secondary: multiplication network / structured collision search

Inspect the normalized mod-3 contraction of Rudich–Rousseau and the exact
branch de-duplication. Check cost over two successive compositions, not
only the last scalar. Low/high cut rank is exponential in the half-width;
symbolic affine-line representations escape that cut bound but require
a real composition algorithm. Reconstruct the G^{3/4+o(1)} short-target
multipoint bound and benchmark the actual evaluation including setup.
A shortlist must retain genuine one-sided collisions at non-negligible
frequency, and have length <G^{1/2-epsilon}. Inert-prime filtering alone
has not been shown to meet that threshold. A complete-period modular level
count is not the desired cutoff oracle.

## 3. Paused: smooth-relation construction

For constant improvements at L[1/3], keep the beyond-Johnson CRT decoder
target (cost exponent a<1.1416 in the repo model, with all the availability,
fresh-yield and fruitfulness conditions). The new two-donor/two-generation
jump bound reportedly allows N-dependent selection and free donor divisors
but still has vanishing probability on the specified average distribution;
the original proof and exact range of window sizes are unavailable.
Three-plus generations and conditioning on nonempty windows were not
closed. Do not treat that result as an arbitrary adaptive lower bound.
Only reopen for a concrete mechanism reusing failed sub-threshold matches
or locating collisions cheaper than generating them. Compare work per
fresh, fruitful relation with a matched sieve, not raw smooth counts.
See [crtdecode](crtdecode_results.md), [sicheck](sicheck_results.md),
[lowerbound](lowerbound_results.md), [comp](comp_results.md).

## 4. Parked: binomial midpoint, hidden characters and conductor

The new midpoint audit reportedly finds the first noninvertible coefficient
at distance ceil(p/4) from m=(N-1)/4 for N=pq=1 mod 4, p<q; shorter local
neighborhoods allow arbitrary midpoint value, and a gcd at the first
singular step factors. The reflection trace on the full solution space
contains the bit but the alternating kernel/cokernel trace cancels it.
Original late scratch proofs missing; keep scope restricted and rederive
before citation. Existing [doubling](doubling_results.md) covers whole-section
obstructions and the fixed CM evaluator. No cheap uniform variable-index
test was found.

The [QRP/conductor note](qrp_conductor_nov26.md) records the complete local
character, prefix-parity factoring reduction, padded class-number tests
and fixed-window arithmetic. Reopen only for a shrinking-argument native
evaluator or an explicit conductor decoder exploiting r|N. Merely rewriting
one bit via another representation is insufficient.

## Research hygiene

- Label elementary proofs, externally conditional reductions, measured results,
  conjectures, and conversation-reconstructed reports separately.
- Preserve original preregistrations and data; distinguish fresh verification
  from reproducing an earlier run. Do not publish inferred raw data.
- Cost every operation, including construction, empty calls, auxiliary prime
  selection, class-group lifting, false positives, and verification.
- Keep the global goal below L[1/3] distinct from improvements in exponential
  deterministic complexity. No new sub-L[1/3] algorithm is established here.
- All late-stage asymptotic claims are internally reported, not independently
  human-refereed; original scratch material is missing for some. See
  [provenance](../experiments/nov26_consolidation/README.md).
