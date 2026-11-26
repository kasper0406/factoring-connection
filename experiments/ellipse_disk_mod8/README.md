# Whole-disk factoring with masked parity

**Result.** For a balanced semiprime, O(log² N) one-bit answers suffice to recover a verified factor. Each answer is the parity of a positive-quadrant disk count with zero, one or two specified bits of y required to be 1. This strengthens the earlier [periodic buckets](../ellipse_row_batches/README.md) and removes the modulo-N moment from the interface. An unmasked signed-vector count by itself does not supply the masked queries.

**Implementation result.** A carry-based binary decision diagram computes these queries without scanning w, y, or the successive norm levels. It works on the test cases. The measured intermediate diagrams dominate the final diagrams, and no improved running-time exponent was obtained.

## Guidance from existing work

* The concurrent [floor8 result](../../notes/floor8_results.md) uses reflection on variable-size squares to read the parity of multiples of the hidden factors, and verifies extraction by bisection. It provides a second reason why an unsplit small-modulus count can contain factor information. Its complete-period simplifications do not evaluate a clipped ellipse. On our split/split branch the full-period congruence count modulo 8 is already public: (2p−1)(2q−1) is 1 if N≡1 mod 4 and 5 if N≡3 mod 4. The clipped disk is essential here.
* The [multiplication-network note](../../notes/shared_collision_nov26.md#1-multiplication-counting-modulo-3) requires a unique assignment per representation. Our auxiliary j is uniquely determined by w,y; each carry is deterministic before partial summation. We count positive-quadrant orbits because taking the parity of a full signed count would immediately lose information.
* The [all-shell join](../ellipse_aggregate/README.md) already removes the per-j loop. It uses a long coordinate list. Its full-interval product test loses the requirement that both local congruences hold at the same integer coordinate. Here the integer relation and its disk bound are enforced jointly.
* The [coefficient note](../../notes/coef_results.md) treats a parity identity and a fixed-form automaton question. The [modular-form audit](../../notes/nativemod_results.md) discusses the missing local prime data at composite indices. Neither supplies a native oracle for this variable-discriminant, bit-filtered disk. The [analytic note](../../notes/nativeana_results.md) treats restricted linear extraction models. Its barriers do not apply to every nonlinear disk algorithm.
* The [clipped-row and weighted reductions](../ellipse_row_batches/README.md) already isolate the first mixed shell. The change here is which information one can extract using tiny counts.

## Disk and parity reduction

Use the public branch ℓ=1 constructed [from N alone](../ellipse_aggregate/README.md). The proof also covers the earlier eligible branches in the tests. Set b=Δ mod 2, w=|2x+by| and take y positive. For B≤A=floor(√(Δ/3))<min(p,q), define a query to count

    w>0, y>0, w²+Δy² ≡ 0 mod 4N, w²+Δy² ≤ 4NB.

Optionally require one or two specified binary digits of y to be 1. Return only the parity. The displayed congruence enforces w≡by mod 2. The coordinate axes contain no counted positive norm in this range. One positive pair thus stands for one sign/conjugation orbit. For the unfiltered disk, multiplying the positive count modulo 4 by four gives the original signed count modulo 16.

**Public subtraction.** In the decoder range B<Δ/4, the positive principal pairs are (kc,kℓ) for k=1,…,floor(√B). Subtract their count modulo 2. For testing bit i on a disk prefix, their parity is

    Σ_{k=1}^{floor(√B)} floor(kℓ / 2^i) mod 2.

The usual Euclidean floor-sum algorithm evaluates this in polynomial bit time. At a single shell there is at most one public pair to subtract.

**Locate the first mixed shell using only parity.** Write a for the mixed reduced minimum. For B<2a, the number of mixed positive-quadrant pairs through B is at most three. In either mixed lattice two non-antipodal vectors below 2a cannot occupy the same nonzero parity class, by the half-sum and half-difference argument. Conjugation exchanges the two mixed lattices, so passage to positive quadrant retains one orbit per antipodal pair in one of them.

If there are two mixed orbits, their positive y coordinates are distinct, even when the two norms differ. To see this, project each orbit into the same one of the two mixed CRT lines. For purposes of this proof only, that line specifies a single slope w/y modulo N. For a fixed |y| the possible positive w values therefore have residues ±t y, for that slope t. Yet for B≤A we have w≤2√(NB)<N/2 for p≥17, since B<p and q>16. At most one positive w in this range has either residue. No slope is computed. If 16A≥N, the CLI handles the input by bounded trial division: the balance and cap inequalities then force p<17 and N<512. Equality of y would now force the same positive pair.

The total parity detects 1 or 3 orbits. If it is even, query successive single-bit masks. For two distinct y's some bit has odd count, while no mask has odd count for zero orbits. These tests therefore detect whether anything is present over the whole interval B<2a. Doubling, capped at A, and binary search locate a in O(log N) cutoff evaluations and O(log² N) parity queries.

**Read out a factor.** The first shell has one or two positive mixed pairs. Three shortest directions are possible only for a hexagonal form and would give Δ=3 for a primitive form. Difference each masked disk answer at a and a−1 and subtract the principal shell if present. With one mixed pair the surviving bits give y. With two, they give the bitwise XOR of the two y's. Pick a set bit i of this XOR; exactly one of the two y's has that bit set. Asking for bits i AND k for every k recovers all bits of that one coordinate. Only O(log N) further queries are needed.

In either case take the integer square root w of 4Na−Δy². The determinant (cy−ℓw)/2 has a proper gcd with N because this mixed line agrees with the public line modulo exactly one prime. Check the proposed divisor by exact division. There is no implicit access to a mixed root, weighted moment, point coordinates or factor inside the query algorithm.

For example the recorded branch (N,ℓ,c,Δ)=(8927,7,141,323) has a tie at a=9. The mixed positive y coordinates are 14 and 29. The single-bit answers give XOR 19. Its low bit selects the odd coordinate, so the two-bit answers reconstruct y=29. Then w=223 and gcd((141·29−7·223)/2,8927)=79. The other factor is 113.

**Alternative mod-8 decoder.** With a mod-8 positive-quadrant oracle, the unmasked disk suffices for locating a and detecting whether it is tied. At the shell a, single-bit mask counts are exactly 0, 1 or 2. Summing these with bit weights gives s=y₁+y₂ in the tied case. In a mixed reduced basis, let d=√(4a²−Δ) be the absolute middle coefficient. The earlier Gram identity gives y₁y₂=(a s²−N)/(2a±d). Try both signs, solve the quadratic, and verify. These modulo-eight disk queries require O(log N) calls. Their single-bit restriction still adds information beyond an unmasked query.

## Global evaluator and its measured bottleneck

The relation w²+Δy²=4Nj is compiled in increasing order of binary columns, with j represented by O(log B) Boolean variables and constrained through a comparator. We never instantiate a shell per value of j. Conditions on every reachable signed carry are reduced ordered binary decision diagrams. Each column adds its square monomials and linear terms and checks parity before dividing the carry by two. The final carry must vanish. A single diagram represents a range of cutoffs; memory and time for its construction are charged.

The comparison modes are:

1. Retain the Boolean relation, then sum its assignments modulo 8. Three skipped free digits annihilate a branch in this ring.
2. After the last column using a particular w bit, add its low and high cofactors at fixed carry in Z/8Z, using three binary planes and discarding the fourth carry.
3. Use F₂ throughout, XORing the two cofactors. The stronger decoder above consumes these one-bit answers directly.
4. Fix each requested cutoff and also eliminate bits of j as soon as their last arithmetic use is past. This gives up reuse across unlike cutoffs. The experiment charges every rebuilt diagram and reports the sum of cache lookups.

The low/high addition happens during construction in modes 2–4. In the first three modes the remaining coordinates y,j determine at most one positive w. In fact for the implemented column order, if any two partially eliminated w assignments survived on the very same remaining assignment and carry, all future arithmetic would be identical. The contribution left by the eliminated variables is solely the difference of the two positive squares, which cannot vanish. This explains why merging final outputs is cheap and why pure w elimination provides little early arithmetic cancellation. The recorded overlap counter checks this explicitly. For mode 4 the cutoff is fixed; a small final relation still does not bound its construction cost.

**Separate rank calculation.** In the complete-period matrix with one axis w mod N, the other y mod N and entry 1 if w²=−Δy² mod N, the rank over F₂ is exactly (p+1)(q+1)/4. The public unit square root of −Δ makes the two sets of square residues equal. Grouping equal residues gives one disjoint all-ones block per square. Conversely, selecting one representative per residue exhibits an identity minor. This only rules out a low-rank decomposition across the full w|y coordinate cut. It gives no lower bound for the *clipped scalar*, mixed variable cuts, or other arithmetic representations.

## Verification and use

The public sample contains the two surviving historical collision inputs and the replacement 192-semiprime sample. On all 194 branches and ten extra, deliberately selected tied-minimum branches, both decoders returned verified factors. For the 194 branches alone, 12,754 distinct cutoff/mask combinations, including unmasked cases, zero cutoffs and every individual bit, were compared with direct independent area counts in both rings. The ten extra branches added 853 combinations. The parity decoder also exercised two-bit masks on the tied cases. At N=11009 on the historical Δ=20 branch, larger cutoffs 16, 32 and 64 yielded positive counts 22, 46 and 90, returning 6, 6 and 2 modulo 8.

On the 23-bit N=6,494,011 the smallest mixed norm was 28. Early elimination of w gave an 80-node final representation at the largest built disk, while construction reached about 81,000 live nodes. Retaining the full relation gave 150 final nodes and a similar peak. All three shared-cutoff variants performed about ten million BDD cache lookups. Fixing each cutoff and projecting j as well used ten builds and 31.3 million cache lookups, with the same peak and no overlapping cofactors on this input. The per-input data also include a three-order comparison. These are diagnostic resource counts, not a fitted exponent. The existing geometric factorer is already faster on this input.

From the repository root:

    uv run --no-project --with dd==0.6.0 --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/decode.py 6494011 --parity-masks
    uv run --no-project --with dd==0.6.0 --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/decode.py 438751 --parity-masks --project-quotient
    uv run --no-project --with dd==0.6.0 --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/benchmark.py --output results-rerun.json

The CUDD backend is provided by the dd wheel. Elapsed time and cache counters can vary between runs. The CLI constructs a branch from N alone; tests of supplied branches maintain continuity with previous work.

The remaining target is computing a *single masked parity of a clipped congruence* with small intermediate state. A useful next test is to sum pairs of coordinate bits before low-column equality fixes their complete carries and retain only a signed or parity-valued boundary condition. A candidate must prove that configurations it merges have the same effect on all future high-bit disk bounds. Reducing the final Boolean relation alone does not meet that test.

The [translation and reflection test](translation.md) gives an exact cancellation law for batches of shifted copies of the mixed weight. It narrows the next recurrence search to the clipped scalar or another representation.

The [centered form recurrence](centered.md) gives a second interface: unmasked signed mod-16 queries on changed forms. All 204 factoring controls passed. Its native geometric evaluation still branches over individual primitive directions.
