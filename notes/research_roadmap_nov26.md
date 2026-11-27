# Research roadmap — 2026-11-26

Goal: classical factoring of general semiprimes, ideally polynomial in log N
or below L_N[1/3]. Do not conflate an improvement in N^c deterministic
complexity, an improved NFS constant, a new reduction, and a new general
factoring algorithm. The internally audited prior results and their
hypotheses remain in [STATUS](../STATUS.md).

The [November 27 checkpoint](../experiments/polylog_factoring_2026_11_27/README.md) records follow-up reductions, reference evaluators, and predictor tests. The open objective is a native polynomial-bit-time evaluator.

**November 27 triage.** An independent audit verified the uniform reduction to one-sided gap comparisons for a clipped modular hyperbola, including repeated factors, recursion, bit widths and a conservative $O(n^4)$ comparison budget per constant-success round. Exact arbitrary affine half-window counts are #P-hard, with $\oplus$P-hard parity, even when the modulus factorization is supplied. Future work should therefore target the special unscaled or doubling-only windows, or the relaxed one-sided gap. Generic exact affine isolation is no longer a priority. [Gap audit](../experiments/polylog_factoring_2026_11_27/single-quadratic-audit.md) · [Affine hardness](../experiments/polylog_factoring_2026_11_27/affine-root-hardness.md)

On the balanced fixed-cap branch, literal floor-sum reciprocity has $\Theta(N^{1/2}\Delta^{1/4})$ terms and the standard theta, Hecke and finite Poisson transforms retain either the high coefficient or the hidden slopes. Genus characters are constant on the queried odd powers. These results narrow the open mechanism to a transform that aggregates the sharp boundary and absolute coordinate jointly, or a special recurrence for the fixed-exponent query family. [Floor audit](../experiments/polylog_factoring_2026_11_27/fixed-cap-carries.md) · [Transform audit](../experiments/polylog_factoring_2026_11_27/theta-hecke-transform.md) · [Cyclic-class audit](../experiments/polylog_factoring_2026_11_27/cyclic-class-evaluator-audit.md)

## Active question for classical polynomial time

Let $n=\lceil\log_2(N+1)\rceil$. We have polynomial-time reductions to explicit sums and verified divisors on finite inputs. Every public evaluator tested so far still takes more than polynomial work in $n$ in its stated worst-case bound.

| Priority | What the reduction already gives | What a new algorithm must compute |
|---|---|---|
| Balanced semiprimes | On one public near-square branch, the first absolute-coordinate moment modulo $N$ at a single cutoff gives a verified divisor. Five even moments at the same cutoff also recover the minimum. | Compute $\sum w\bmod N$ over positive matches in the whole disk. The square-residue join scans a power of $N$. Check a formula at the prescribed cutoff itself, where a common recurrence for all prefixes is unnecessary. [Single-moment proof](../experiments/ellipse_aggregate/weighted_recurrence_single.md) · [Five even moments](../experiments/ellipse_aggregate/weighted_recurrence.md) |
| General integers | A randomized reduction to a one-sided clipped modular-hyperbola gap comparison covers repeated factors and arbitrary odd composites. A deterministic route uses $O(n)$ signed counts mod 4 per split. Even factors and prime-power leaves are explicit. | Evaluate the gap in polynomial bit time, or evaluate the special unscaled root window. Exact arbitrary affine root-window counts are #P-hard even with supplied factors. Tangent blocks plus exact Euclidean floor sums give a public root-prefix control with $\Theta(N^{2/3})$ blocks at central cutoffs. [Gap audit](../experiments/polylog_factoring_2026_11_27/single-quadratic-audit.md) · [Affine hardness](../experiments/polylog_factoring_2026_11_27/affine-root-hardness.md) · [Prefix control](../experiments/ellipse_disk_mod8/native_dyadic.md) |
| Independent arithmetic control | For $p<q<2p$, $B=\lfloor\sqrt{pq}\rfloor$ gives $\gcd(\binom{2B}{B},pq)=q$. | Extract that middle coefficient modulo $N$ at the one prescribed index. All 128 tested unweighted cyclic folds miss the divisor on two fixtures. [Valuations, rank bounds and fixtures](direct_split_poly.md) |

First invest in computing a residue on unfactored input. Charge every arithmetic operation needed to build the representation and to answer the query, including empty blocks. A smaller geometric subproblem makes progress only if the number and size of all child calls also give a polynomial bit bound. On the balanced branch, check the two collision inputs and a growing cap. On general inputs, require the same evaluator to handle repeated factors, where points on the horizontal axis can change the parity. Keep factor-supplied calculations as independent controls.

**Next experiment.** Work directly on the clipped modular-hyperbola contrast after the change of variables $xv=c_s^2+r+N^2q$. Seek a transform that aggregates many consecutive quotient values $q$ while preserving the signed complement cancellation to additive error below $N/2$. Test the unscaled half-window first as the smallest exact falsification target. Reject a construction if it enumerates $q$, divisors, lower blocks, upper digits or a power-sized Fourier spectrum. [Reduction and required resolution](../experiments/polylog_factoring_2026_11_27/single-quadratic-audit.md)

The balanced fixed-cap experiment remains a secondary control. Integer carries and the ellipse boundary can couple the two prime congruences, but literal floor-square reciprocity retains $\Theta(N^{1/2}\Delta^{1/4})$ terms. A future formula must contract exact square boundaries across quotient cells rather than evaluate adjacent cumulative ellipse counts separately. [Carry audit](../experiments/polylog_factoring_2026_11_27/fixed-cap-carries.md)

A candidate normalized residue $z$ can be checked immediately by taking $t=\lfloor\sqrt{z^2\bmod N}\rfloor$, checking the square and verifying a divisor by gcd. This also makes weaker randomized mechanisms testable. For a uniform unit $z$ on a balanced $N=pq$, exactly $2(\lfloor\sqrt N\rfloor-1)$ of the $\varphi(N)$ units give a split by this test. Indeed, each positive $t\le\lfloor\sqrt N\rfloor$ except $p$ is a unit, its ordinary square is distinct and below $N$, and two of the four CRT signs give a split. Uniform guessing needs order $\sqrt N$ trials on average. A randomized proposal targeting polynomial expected time needs an inverse-polynomial success guarantee in the **bit length**, together with a polynomial cost per trial. The fixed-cap residue supplies a deterministic successful $z$ or an immediate gcd.

**Why the fixed-discriminant bit moved down.** Suppose a deterministic reduction uses at most $t$ yes/no oracle answers and polynomial additional work to return a verifiable divisor. One can run every answer transcript, verify each claimed divisor, and retain a correct split in $2^t\operatorname{poly}(n)$ time. If only the true path has a proved polynomial runtime, apply that bound to each simulated path. A single bit, or even $O(\log n)$ adaptive bits total, can therefore be eliminated at polynomial cost. The current $-59$ bit also has a sparse input promise. This observation concerns verified search with that total answer budget. A polynomial number of bits, or a residue modulo $N$ carrying $O(n)$ bits, has a different budget. [Endpoint experiment](../experiments/ellipse_disk_mod8/manin59.md) · [Growing-discriminant collisions and coverage](../experiments/ellipse_disk_mod8/growing59.md)

The older canonical-grid target has a separate payoff. Its conditional almost-all $N^{1/8}\operatorname{polylog}N$ bound still has an exponential factor in $n$. Keep the coverage and class-lifting work under that objective.

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

**Earlier weighted prefix route.** With ℓ=1 and c chosen next to 2√N, the squared-coordinate
moment gives a direct factor from N alone. Both the geometric evaluator and
a residue-list join run through verification. Their costs remain powers of N.
The fixed-cap reduction above removes the need for adaptive prefixes. A full-range product test
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
bit on an explicit matched triple. A joint recurrence for continued-fraction edge
visits would improve this benchmark. The terminal-edge multiset has an elementary
formula, and keeping its reduced-denominator weight recovers the known gcd sum.
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
