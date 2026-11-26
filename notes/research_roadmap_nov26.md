# Research roadmap — 2026-11-26

Goal: classical factoring of general semiprimes, ideally polynomial in log N
or below L_N[1/3]. Do not conflate an improvement in N^c deterministic
complexity, an improved NFS constant, a new reduction, and a new general
factoring algorithm. The internally audited prior results and their
hypotheses remain in [STATUS](../STATUS.md).

## 1. Primary: canonical index / public modular ellipse count

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

Next tasks (in order):
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
