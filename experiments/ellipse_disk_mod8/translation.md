# Can reflection or translation cancel the hard part of the ellipse?

**Decision.** Keep the disk as a specification for a factoring oracle. Deprioritize searches for short, fixed translation identities for its completed mixed congruence. Constructing any such exact cancellation already gives the factors. A reciprocity at the particular cutoff could still succeed. Its first milestone is a public transformation with a cheaper boundary computation and a bound on the total work of all children.

The [previous experiment](README.md) provides a decoder and a whole-disk implementation. The [square-reflection experiment](../../notes/floor8_results.md) gives a reason to test pairing. Here we determine exactly what a family of translated or centrally reflected copies can cancel. The [earlier Fourier identity](../ellipse_global/README.md#fourier-transform-and-the-2-adic-test) supplies the completed mixed weight used below.

## The congruence and the observable

For an eligible branch write $N=pq$ with distinct odd primes $p<q\le 2p$. Define public split coordinates modulo $N$ by

    u = ℓw − cy,       v = ℓw + cy.

The coordinate change is invertible modulo $N$. Since $c^2+\Delta \ell ^2=4N$, the congruence in odd characteristic is $uv=0$. The physical candidate set still enforces $w\equiv \Delta y\pmod 2$, positivity, the disk inequality, and any requested y-bit mask.

Define two subgroups of $(\mathbb Z/N)^2$, using the unknown factors **only in the proof and validation**,

    H = {(u,v) : p divides u, q divides v},
    J = {(u,v) : q divides u, p divides v}.

Both have size $N$, their intersection is the origin, and any coset of H meets any coset of J in exactly one point. Write $M=1_H+1_J$ as an integer periodic function. On the disk through $B<p$, it counts exactly the mixed contribution. Intersections of public and mixed lines require a prime s dividing both physical coordinates, forcing $Ns$ to divide the form value. The origin is excluded. Outside the disk this completion assigns multiplicity two at the origin and can include points lying on public lines. It matters which global function is being paired.

## What a batch of translates cancels

Let $S$ be distinct public offsets modulo $N$. Add their translates of $M$ in characteristic two. Let $a,b$ be the numbers of H-cosets and J-cosets respectively that meet S an odd number of times. Then the number of anchors with odd output is exactly

    a(N−b) + b(N−a).

The selected H-cosets and J-cosets are two families of parallel lines. There are $Na+Nb$ incidences and $ab$ crossings counted twice. This counts the entire period for any S, including masked and irregular sets. Its proof uses no assumption of randomness.

* For $k=|S|$ offsets with distinct images in both quotients, the support is $2k(N-k)$. The only cancellations are the isolated crossings. An H- or J-coset collision between two different offsets gives a difference in H or J. At least one of its coordinates then has a proper gcd with N, since the difference is nonzero modulo N. For a public progression of length $k\le N$ whose step has two unit coordinates, both quotient maps are injective throughout.
* Every nonempty, identically vanishing batch with $k<N$ factors N. The formula forces $a=b=0$, since $a=b=N$ needs at least N offsets. Fix any offset. Its H-coset is even, so some other offset shares that coset. Testing both coordinate differences against N for every other offset takes at most $2(k-1)$ gcds.
* The same implication holds for a linear identity with coefficients modulo 8 or 16. First combine equal offsets and discard zero coefficients. Divide integer representatives of all coefficients by their minimum common power of two, then reduce the identity modulo two. At least one odd coefficient survives and the preceding extraction applies to those offsets.
* Conversely, given p and q choose nonzero $h=(p,q)\in H$, $j=(q,p)\in J$. The four-term difference $(1-T_h)(1-T_j)M$ vanishes over the integers. Four offsets suffice and fewer cannot. For a null parity identity below N, the occupied cosets in each family must all have even occupancy. With two distinct offsets they cannot collide in both families, and with one or three at least one coset is odd. Thus a polylogarithmic-size explicit translation certificate is equivalent to factoring in polynomial bit time. The strict $k<N$ matters. Summing N successive unit-step translates vanishes in parity without given factors, and expanding that identity costs N terms.
* Any nonzero parity residual has support at least $2(N-1)$ on the full period. For $1\le a,b<N$, the bilinear formula takes its minimum on the four corners of this rectangle. For $a$ or $b$ equal to 0 or N, the parity relation $a\equiv b\equiv k\pmod 2$ gives the same bound. A single translate attains it. A large support can still have a compact symbolic description.

A central reflection at offset $h$ replaces $M(z)$ by $M(h-z)=M(z-h)$, since M is even. The only center yielding equality for every anchor is $h=0\pmod N$. At that center the positive-quadrant geometry has already identified the signed counterparts. More generally, a real affine isometry of the whole centered ellipse that bijects its physical integer lattice fixes the center. The principal form's unique shortest pair is horizontal because a nonzero y has norm at least $\Delta /4>1$. The isometry must preserve the horizontal axis and its orthogonal complement. Each positive orbit is therefore fixed after folding signs. This statement concerns a single disk and its lattice; the square comparison uses different domains.

**Connection to the sharp cutoff.** Let m be the first mixed norm divided by N. In any positive disk through $2B\le m$, the difference d of two distinct candidate offsets has $Q(d)<2NB\le Nm$, interpreting d in physical form coordinates. The strict inequality follows from the positive inner product in $w,y$. A difference in H or J would have norm a positive multiple of N. Below $Nm$ it cannot be primitive mixed by minimality, nor singular since $m<p$. Hence the quotient projections of the entire disk, or any mask of it, are injective up to this cutoff. Even a disk whose mixed count at the required anchor is zero can have a large nonzero response over translated anchors. This distinction is precisely where an anchor-specific recurrence could escape the result.

## A centered recurrence that really shrinks the disk

Consider the elementary recurrence at the requested center. Split each physical integer coordinate x,y by its low bit, change variables in the disk, and repeat. Each child has half the radial scale. After t steps, let d be the stride, with d = 2^t. A child corresponds to a residue pair r modulo d. Because d is invertible modulo N and both mixed subgroups are closed under scaling by d, its remaining periodic weight is a translate of the original by the public offset r/d in split coordinates. Different residue pairs give different translated kernels for d < N, by invertibility of the change of coordinates and the trivial translation stabilizer proved above.

At B = 1 on the public ℓ = 1 branch, let L be the integer part of √(4N/(Δ+9)) and take the largest power of two d ≤ L. Every child at this depth is geometrically nonempty. The square 1 ≤ x,y ≤ d supplies representatives of all d² residue pairs, with w = 2x + (Δ mod 2)y ≤ 3d and w² + Δy² ≤ 4N. The branch satisfies Δ < 4√N + 4, so the literal recurrence keeping these calls separate creates Ω(√N) nonempty children even at the first cutoff. The count in any given child may of course be zero. The test checks each representative and each translated-kernel offset using only N.

For N = 38,942,974,393 the stride was 512, giving 262,144 distinct nonempty children at B = 1. This example shows why shrinking the geometric radius alone is insufficient. It leaves open symbolic evaluation of many children together, early certificates for a clipped child, and other choices of recurrence.

## A direct test of a Euclidean shear

In split coordinates the unit-determinant shear $(u,v)\mapsto(u,v-tu)$ keeps the primitive part of H on the cone only when $q\mid t$. Indeed $v=0\pmod q$ and u is a unit there. The primitive part of J survives only when $p\mid t$. A unit coefficient therefore retains no primitive mixed point on the same cone. Keeping only one family already gives a proper gcd from t.

More generally, an invertible matrix preserving the full split cone must permute the two coordinate axes locally at each prime. Its columns are the images of the axes. A linear map that chooses different permutations at the two primes has a matrix entry divisible by just one prime. If both choices agree, it is diagonal or anti-diagonal modulo N. A Euclidean swap which **mutates the form** escapes the fixed-cone argument. It still needs a descent measure and a way to pay for the transformed cutoff.

A useful control for any claimed evaluator covering arbitrary polygons in split coordinates is the triangle $1\le u<v\le T$. For $\lceil\sqrt N\rceil\le T\le\lfloor\sqrt{2N}\rfloor<2p$, the literal congruence count is exactly $[q\le T]$. Within that triangle the only possible solution is $(p,q)$. Its parity locates q by bisection. This control needs no mask, curved boundary, or preliminary class computation. An ellipse-specific algorithm need not handle the triangle.

## Checks and research priorities

The [standard-library script](translation_check.py) and [results](translation_results.json) reproduce the following.

| Check | Result |
|---|---|
| Full periods for N = 15, 21, 35 | Tested every nonzero reflection offset, every shear parameter, unit-step progressions through length N, and a four-term factored identity at every anchor. |
| 194 public branches plus 10 deliberate ties | All 1,473 sharp-cutoff and mask queries agreed between direct mixed counts at the origin and factor-backed quotient parities over 5,304,161 candidate coordinates in total. |
| Early prefix $2B\le m$ | All 297 nonempty tested sets of offsets had injective projections into both quotients, as predicted. For the full set of tested cutoffs, 822 of 1,158 nonempty filters had injective projections. |
| Unmasked first mixed shell | 192 of the 194 public branches had one mixed orbit, the other two had two. The ten deliberate ties all had two. This sample is the repo's two preserved collision cases and 192 replacements. |
| Two-coordinate dyadic descent | Explicit physical representatives and distinct kernel offsets checked for every child at the selected depth on five N-only branches, from 64 to 262,144 children. |
| Triangle control | Direct enumeration agreed with the threshold at all 12,043 tested cutoffs across the 204 branches. |

Run from the repository root:

    python experiments/ellipse_disk_mod8/translation_check.py

The fixture factors and first norms enter only the controls, including the choice of cutoffs. These are algebraic checks, with no runtime fit. The translation theorem assumes an identity of the completed mixed weight at **every** anchor. A sum of implicitly many terms, a nonlinear identity, an alternative completion, or a formula specialized to the single clipped scalar can evade it. The theorem gives no general complexity lower bound.

The next experiment worth funding is a recurrence for the parity at the actual anchor, initially with one sharp cutoff and no y-mask. Require an exact transformation, a shrinking parameter, and a bound on the sum of work across recursive calls, including construction of any exceptional boundary. Check it at the first nonzero mixed shell and its predecessor, where losing one point flips the answer in the singleton cases. Then test one- and two-bit masks for the direct factor decoder, or show compatibility with the original unmasked class route. Further decoder or class-lift benchmarks do not address this uncertainty. A claim of improving on the public deterministic exponent needs a complete factoring cost below the [Harvey–Hittmeir bound](https://doi.org/10.1090/mcom/3708) for the stated input family, with any restrictions or failure sets explicit.
