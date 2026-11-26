# Counting the primitive nonsquare band across directions

**Decision.** Keep the modular-root lookup for rough discriminants. Park the search for pairs of actual mixed directions to cancel below $\sqrt{\Delta}/2$. A determinant bound shows that a second direction cannot occur there. The class-group calculation gives a cheap way to propose norm levels on the $53^2$ fixture. Verifying those levels directly was much cheaper than deciding membership in a list of prime classes. Both methods still have power-of-$N$ work.

We evaluated the [upper band left by the rectangle experiment](rectangle_rays.md#next-test). Its output is the parity of positive base points at nonsquare levels $B/4<t\le B$. Base points use the previous experiment's content convention. This bit is one contribution to the full disk count. The direct verifier also produces the exact band count. All setup uses public branch data. Modular square roots in the direct verifier are taken modulo prime powers in the public discriminant.

## How much cross-direction cancellation is possible?

Let $F$ be a mixed reduced form of discriminant $-\Delta$. The repository identifies its nonzero vectors up to negation with the positive mixed orbits, pairing the two inverse mixed classes by conjugation. For integer vectors $u,v$, let $b_F(u,v)=F(u+v)-F(u)-F(v)$. The two-by-two Gram determinant gives
$$\Delta\det(u,v)^2=4F(u)F(v)-b_F(u,v)^2.$$
If both norms are at most $B$ and $4B^2<\Delta$, every such pair is collinear. At the entire decoder cap $B\le\sqrt{\Delta/3}$ every independent pair has determinant $\pm1$. There are then at most three directions up to sign. Once two directions form a basis, a third must have both coefficients $\pm1$. The two possible additional directions have determinant 2 with each other. In particular the 581/582 benchmark asks whether a single possible mixed direction is occupied. Pairing distinct true directions inside its mixed lattice cannot prune the empty candidates.

## Count in a smaller discriminant

Write $C$ for the hidden prime class, so the mixed form belongs to $C^2$ up to inversion. The principal contribution through the public cap occurs at square levels and drops out of this query. For even discriminant use $(L,D)=(N,\Delta/4)$ and for odd discriminant use $(L,D)=(4N,\Delta)$. A physical point obeys $x^2+D y^2=Lt$. Put $E=\prod_{p^2\mid D}p$, with one factor per prime. A base vector in a mixed form has no common coordinate divisor prime to $E$ greater than 1.

To see why the content test transfers, the mixed lattice has index $N$. A prime shared by both physical coordinates cannot divide $N$ in this range, since it would force its square to divide $Nt$ in the original integral form. The lattice change therefore preserves divisibility by primes relevant to a removable multiplier. On the odd branch, after dividing $w,y$ by a removable $h$, the equation $w_0^2+\Delta y_0^2=4Nt_0$ restores their equal parity. Thus the integral coordinate $(w-y)/2$ also divides by $h$. This gives the same base convention in physical, mixed and reduced coordinates.

Take the largest odd prime $f$ such that $f^2\mid\Delta$, or $f=1$ if none exists. Write $\delta=\Delta/f^2$. We use one prime descent so its incident sublattices are ordinary lines in $\mathbb F_f^2$. Enlarge the quadratic order from discriminant $-\Delta$ to $-\delta$. Extension is a homomorphism on ideal classes, so the parent of a square class is again a square. Enumerate primitive reduced classes at $-\delta$ up to inversion, square them by Gauss composition, and retain distinct parent squares. This costs one leading-coefficient probe through $\lfloor\sqrt{\delta/3}\rfloor$, including modular roots and all compositions. The orders of discriminant $-3$ and $-4$ each have one parent class. [Class extension for a prime-index overorder, slide 17](https://math.mit.edu/~drew/PANTSXX.pdf)

A primitive child form is a restriction of its parent to an index-$f$ sublattice. In elementary terms, a primitive child modulo odd $f$ has rank one. Put its radical in the first coordinate. Then $f$ divides the first two coefficients, the third is a unit, and the discriminant identity forces $f^2$ to divide the first coefficient. Dividing that basis vector by $f$ gives an integral primitive parent. Conversely an index-$f$ sublattice is the inverse image of a projective line. Keep its restriction if its coefficients are primitive. Canonical reduction up to inversion handles equivalent restrictions; their counts agree.

For a parent $[a,b,d]$, enumerate vectors up to sign by the exact row bound $(2au+bv)^2\le4aB-\delta v^2$. Keep nonsquare $t=F(u,v)$ in the band whose coordinate content is supported on $E$. Passing to a child of index $f$ preserves this last test at every prime outside $E$, since $f\mid E$. A nonzero vector modulo $f$ belongs to exactly one child line. A zero vector belongs to every line. These statements group vectors of different levels and directions before generating their physical progressions.

We compared two uses of this construction.

* **Select levels, then verify points.** From the parent vectors collect a set of possible levels. If the vector is nonzero modulo $f$ and $f\mid t$, its line is isotropic and its restriction has all coefficients divisible by $f$, so discard it. Retain a zero vector's level conservatively. Every physical hit lies on the resulting list. For each distinct proposed target compute the roots of $x^2=Lt$ modulo the prime powers in $D$. The existing rectangle routine combines roots by CRT modulo $D$, checks content and tests $(Lt-x^2)/D$ for an integer square along each progression. No factor of $N$ is used or needed. This method removes the old seed scan of width at least $B$. Its parent catalogue has a separate setup cost.
* **Select classes, then verify membership.** Within each parent toggle the parity for each vector's incident line; a vector on every line toggles a common bit. Reduce primitive restrictions of the odd lines. The existing conic and class-duplication routines lift the selected child squares, and the existing batched verifier tests whether a listed class represents a proper divisor of $N$. Exactly the hidden prime class or its inverse can represent such a divisor. An odd answer yields a factor witness as a byproduct. With $f=1$ the catalogue already supplies the square preimages. This route incorporates cancellations within classes before norm verification. [Lift and membership proofs](class_prefix.md)

## Results and controls

| Public $N$ | $\Delta$ | $B$ | Band counts at $B-1,B$ | Proposed levels at $B$ | Physical square tests, old / proposed |
|---:|---:|---:|---:|---:|---:|
| 140,000,556,000,551 | 5,820,248 | 582 | 0, 1 | 24 | 16,546 / 15,210 |
| 138,468,065,327 | 1,410,683 | 579 | 0, 1 | 77 | 1,694 / 1,463 |
| 39,866,375,420,117 | 1,639,243 | 281 | 0, 1 | 19 | 4,299 / 3,919 |

The last two inputs are prime-discriminant stress controls selected from the existing held-out cohort for a nonsquare, coprime, unique-direction first hit. One first norm is the prime 281, the other is $579=3\cdot193$. On these inputs the old rectangle projection scanned the whole prime modulus to establish its residue table. A simpler control scans the norm band and asks for roots only at the proposed residues. A second control also rejects levels having an inert prime divisor before asking for roots.

| First norm | Original rectangle | All levels, targeted roots | Inert sieve first | Parent levels, targeted roots | Class membership |
|---:|---:|---:|---:|---:|---:|
| 582 | 8.6 ms | 11.1 ms | 13.9 ms | 6.8 ms | 524 ms |
| 579 | 151 ms | 19.0 ms | 11.5 ms | 9.1 ms | 63 ms |
| 281 | 178 ms | 10.8 ms | 5.3 ms | 6.0 ms | 50 ms |

Timings are medians of five same-host runs with fresh evaluator objects and rotated order. They include setup, public factorizations, root lookups and verification. The two targeted-root controls and the parent-level method use the same final verifier. Most of the improvement on prime discriminants came from removing the complete residue-table scan. The simpler inert-first control was faster at 281.

At 582 the parent discriminant is $\delta=2072$. The class route probed 26 parent leading coefficients, found three parent square classes, inspected 121 nonzero parent vectors, and found 51 eligible vectors. It selected 35 candidate child squares before the square-class tests. The final 160 prime classes required 149,938 polynomial evaluations for membership. These evaluations have a different cost from integer square tests. On the level route, the 121 parent vectors yielded just 24 distinct proposed levels. The raw record accounts for the 46 targeted modular-root calls and all remaining physical positions.

The checks include 605 band queries on 204 repository branches, including the two historical collision inputs and the tie fixtures. Eight queries had positive even counts; a parity-zero result can therefore coexist with physical hits. Independent signed-form enumeration checked 5,586 class labels in 996 lists through discriminant 1,000, of which 180 lists used descent. Ninety-two labels had more than one independent direction. The timed cases were checked against independent reduced-form counts using their fixture factors only on the control side. [Implementation](conductor_band.py) · [Independent checks](conductor_band_check.py) · [Raw input-level results and timing medians](conductor_band_results.json)

## Why this does not remove a power of N

For distinct primes in $[P,2P]$ and the near-square $\ell=1$ branch, $c=O(P)$ and $0<\Delta=O(P)$. If an integer $f\ge F\ge2$ has $f^2\mid\Delta$, write $\Delta=f^2d$. For each $c,f$ there are $O(P/f^2)$ choices of $d$, and each triple fixes $N=(c^2+f^2d)/4$. Summing over $c$ and $f$ gives $O(P^2/F)$ distinct semiprimes at most. There are $\asymp P^2/(\log P)^2$ prime pairs, so the affected fraction is $O((\log P)^2/F)$. Likewise, $\Delta\le P^{1-\epsilon}$ accounts for a fraction $O(P^{-\epsilon}(\log P)^2)$. For every fixed exponent saving, removing even the full square conductor therefore leaves this implementation's leading-coefficient scan at $N^{1/4-o(1)}$ on almost all such pairs, in the sense of these fixed-$\epsilon$ bounds. In the 96 held-out inputs 79 admitted no odd prime descent at all. The bound concerns the near-square branch and this catalogue construction. [Prime number theorem](https://mathworld.wolfram.com/PrimeNumberTheorem.html)

## Next blocker

Use the prime-discriminant controls at 578/579 and 280/281 to test a bulk calculation on the surviving levels. No conductor shrinkage is available. For an odd prime discriminant magnitude $\Delta>3$ the class group has odd order: the only ambiguous reduced class is principal. Indeed $b=0$ is impossible, $|b|=a$ forces $a\mid\Delta$, and $a=d$ would factor $\Delta=(2a-b)(2a+b)$ with the reduced bounds. Hence every class is a square. Every odd split prime $r$ in the band yields a primitive represented norm via $[r,b,(b^2+\Delta)/(4r)]$, so a catalogue of all geometrically possible levels includes it. At 579 this accounts for 43 of 77 proposed levels. At 281 all 19 proposals are split primes. This is an instance count, with no uniform prime-distribution claim in a growing discriminant.

The test for the next idea is to sum over those possible levels without individually constructing their progressions or scanning the full parent class catalogue. Use both controls so that a mechanism for prime target levels alone cannot pass the composite first hit. Charge preprocessing and compare against targeted roots with the inert sieve. A route confined to cancellations within one occupied mixed lattice has no independent pair to exploit at any of the three benchmark transitions.

From the repository root with python-flint, NumPy and SymPy installed, run

    python experiments/ellipse_disk_mod8/conductor_band.py 140000556000551 582
    python experiments/ellipse_disk_mod8/conductor_band.py 138468065327 579

The class-membership variant uses PARI/GP. The full check also imports the earlier dd-based test module.

    python experiments/ellipse_disk_mod8/conductor_band.py 140000556000551 582 --class-membership --gp /path/to/gp
    python experiments/ellipse_disk_mod8/conductor_band_check.py --gp /path/to/gp --output rerun.json
