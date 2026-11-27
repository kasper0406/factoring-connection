# Classical factoring research checkpoint

The objective is a classical factoring algorithm with polynomial bit complexity, or progress toward a general algorithm below the number field sieve scale. The work here establishes stronger oracle reductions, a faster evaluator for one fixed scalar, and barriers for specific representations. The native polynomial-bit-time evaluation step remains open. [Original November 26 roadmap](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/notes/research_roadmap_nov26.md)

## Working ledger for the quarter-parity evaluator

**Current target.** Given an unfactored odd modulus N and an aligned binary interval I in its lower half, compute the parity of the solutions to x² ≡ 1 (mod N) in I in polynomial bit time. A general evaluator remains open.

**Where the search stalled.** The product, character, correlation, and run-count reductions all leave a global quantity to evaluate. The next benchmark is a public computation on a wide interval where the current evaluator returns unresolved. Report the interval width, bit-operation cost, and assumptions about N. A proved improvement on a specified family is worth recording even before it covers all N. For the CM proposal the corresponding benchmark is an explicit auxiliary value with certified divisibility and size, together with its public construction cost.

| Approach | What we know | What would advance it |
|---|---|---|
| Wide intervals with a small-prime power dividing N | For roots of one, extract a prime power s from N=sm by trial division through the bit length. Two progressions reduce total exceptional length up to s⌊√m⌋ to short quadratic lattices modulo the remaining cofactor. [Implementation](check_lower_half_prefix.py) · [Width proof and checks](single-quadratic-count.md#wide-intervals-via-an-extracted-prime-power) | The width gain is √s. Balanced semiprimes with both primes above the trial bound receive no gain. |
| Balanced-semiprime control for wide queries | The [public small-multiplier baseline](check_lehman_window_baseline.py) tests 256 integer square roots after computing the 32 starting points, then counts the four roots if it finds a split. For k=ij the matching center is ip+jq; its excess over 2√(kN) is (ip−jq)²/(ip+jq+2√(kN)). At p sizes 16, 32, 48 bits, all 100 draws per size with q=nextprime(⌊3p/2⌋) resolved the second-quarter query. Controls drawing q near a random point between 1.35p and 1.65p resolved 11, 0, 0 of 100 respectively. All answers agreed with CRT. | This is a known [Fermat/Lehman control](../entropy/round1/gates.cpp), with multiplier one included here. Benchmark future interval methods on the misses to measure gain beyond finding an existing square witness. |
| Second-quarter mask | A [dyadic rational sign filter](single-quadratic-count.md#a-compact-rational-filter-for-the-second-quarter) of degree O(log² N) gives total error at most 1/4 over the roots, with no guard band. | Evaluate the rationally weighted root sum from N without enumerating roots. |
| Property sieve for balanced factor pairs | The existing [hyperbolic sieve analysis](../../notes/signsplit_results.md) tests a candidate sum S by requiring S²−4N to be a square modulo each auxiliary prime. At N=2,989,601,243 under the promise p≤q≤2p, parity and the interval from 2√N through 3√(N/2) give 3,316 candidates. Sieving by the odd primes through 31 leaves three, including the true sum 112,716. A [clipped CRT wheel](check_factor_sum_wheel.py) visits 86 individual sums after building residue classes and certifies the lower half empty with 43 probes. At N=118464083331556819902753234409, the same fixed tests leave at least 46,495,088,640 sums in the balanced range: each complete period of 200,560,490,130 integers contains exactly 223,534,080 allowed even sums by CRT. | A [common randomized Jacobi sieve](single-quadratic-count.md#a-common-jacobi-sieve-for-all-candidate-factor-sums) separates every nonsquare discriminant up to N with high probability using O(log² N) public moduli. A [single-modulus block certificate](single-quadratic-count.md#contiguous-blocks-rejected-by-a-single-small-modulus) needs at least W/(L−1) contiguous blocks when the moduli are at most L. On the empty benchmark half, [eight rejection moments plus the interval size](single-quadratic-count.md#certifying-an-empty-interval-from-rejection-moments) certify zero survivors. The large-modulus moments were enumerated. [Small-prime CRT moments](single-quadratic-count.md#computing-fixed-order-moments-from-small-prime-filters) are polynomial-time at fixed order and certify this half using the ten odd primes through 31. [Order nine is necessary and sufficient](single-quadratic-count.md#exact-moment-threshold-on-the-ten-prime-benchmark) for this benchmark from aggregate moments alone. A [degree lower bound](single-quadratic-count.md#a-degree-lower-bound-for-small-prime-interval-certificates) rules out fixed-order nonnegative polynomial certificates for power-width intervals with polylogarithmic primes. Uniform small-prime separation on balanced odd composites implies a [least quadratic nonresidue bound](single-quadratic-count.md#a-least-nonresidue-problem-inside-uniform-small-prime-separation). An [output-size bound](single-quadratic-count.md#the-size-of-an-explicit-survivor-list-on-general-inputs) rules out polynomial explicit lists of all genuine balanced sums on general squarefree inputs. A [two-input balancing reduction](single-quadratic-count.md#two-odd-inputs-suffice-to-balance-an-arbitrary-split) extends exact balanced interval occupancy on odd composites to general factoring. A [single odd squarefree multiplier](single-quadratic-count.md#balancing-with-one-squarefree-odd-multiplier) also preserves squarefree inputs. [Ordinary Fermat search](single-quadratic-count.md#fermat-search-after-primorial-balancing) after choosing among polynomially many such multipliers misses almost all balanced semiprimes in B^(1−ε) midpoint tests for factors of scale B. A [short midpoint interval](single-quadratic-count.md#certifying-a-short-interval-of-fermat-midpoints) can be finished by an unknown-divisor lattice using its exact transformed width. Efficiently locate one surviving sum in a specified interval as the initial range grows; full-period CRT densities alone do not give that algorithm. |
| Factor-sum properties plus lattice recovery | A proposed sum residue modulo prime ℓ gives at most two residues for the larger factor q via z²−sz+N. In each progression write q=c+ℓy, center the interval, then use the existing [certified linear divisor lattice](../../src/nativemod/reductions.py) on y+c/ℓ modulo the unknown q. The [hybrid check](check_sum_residue_lattice.py) returned a factor on all 100 larger-control inputs given the planted residue modulo 17,179,869,263, and rejected a locally admissible false residue on every input. The full balanced range was used. | Produce a short list containing the actual sum residue using only N. Across a full prime modulus the square-discriminant test alone retains (ℓ+(N/ℓ))/2 residues by pairing z with N/z. The supplied residue was side information. The [low-bit decoder](check_sum_low_bits.py) also handles s mod 2^k. For an even representative s=2h, solve (q−h)²≡h²−N mod 2^k. A zero right side gives one progression of stride 2^ceil(k/2). For nonzero valuation 2v, divide out 2^(2v) and lift at most four odd square roots, yielding progressions of stride 2^(k−v); odd valuation is impossible. On the 100 unplanted 48-bit-prime controls, the first successful precision for the fixed rank-9 decoder ranged from 25 to 34 supplied low bits, median 26. This measures decoding given true bits. [Padding by 15 when needed](single-quadratic-count.md#avoiding-precision-loss-when-decoding-low-sum-bits) preserves balance and ensures four full-stride classes for every genuine sum, eliminating valuation loss. On semiprimes 3 mod 4, the [first ambiguous sum bit](single-quadratic-count.md#the-first-ambiguous-sum-bit-as-quadratic-residuosity) decides whether a specified one of ±2 is a square modulo N. Equivalently [4 is a fourth power](single-quadratic-count.md#a-single-quartic-test-for-the-first-ambiguous-bit) in the N+1 sum branch, connecting the bit to the existing quartic audit. The [norm conic](single-quadratic-count.md#where-the-norm-conic-first-detects-the-missing-bit) passes every squarefree modulus on either branch; primitive solubility modulo N² distinguishes the branches. With opposite local signs the [affine conic count](single-quadratic-count.md#counting-the-affine-conic-with-opposite-local-signs) equals twice one hidden prime minus one. [Every nonzero Fourier coefficient](single-quadratic-count.md#the-conic-fourier-transform-away-from-its-origin) is efficiently evaluable; finding one with exceptional value already yields a factor. [Fractional frequencies on the coordinate axis](single-quadratic-count.md#exact-fractional-frequencies-on-a-coordinate-axis) reduce exactly to a geometric sum over multiples of a hidden prime. The [norm-one variant](single-quadratic-count.md#uniform-sampling-on-the-norm-one-conic-and-the-accuracy-needed) permits exact uniform sampling; its cardinality minus (N−1) gives the signed prime difference. [Two samples encode that difference](single-quadratic-count.md#an-interval-logarithm-for-the-signed-prime-difference) as a unique discrete logarithm in its signed half-window with probability at least 6/π². A [known conic point order](single-quadratic-count.md#a-known-order-on-both-components-shortens-the-interval-logarithm) either splits the input by a local order mismatch or reduces the candidate count by the square of that order, up to rounding. [Primitive complex roots of unity](single-quadratic-count.md#reducing-a-characteristic-zero-torsion-point) retain their order at every prime coprime to it, ruling out order mismatch by direct reduction. [Interval bins of any primitive linear projection](single-quadratic-count.md#interval-bias-of-sampled-norm-one-coordinates) on fresh samples have negligible distinguishing advantage with polylogarithmically many bins and draws. The [same-point Jacobi label](single-quadratic-count.md#the-square-class-bit-of-a-sampled-norm-one-point) gives the product of the two local square-class signs; halving a point is equivalent to taking a modular square root of its sampler norm. [Varying the conic](single-quadratic-count.md#prescribing-a-norm-by-choosing-the-conic) embeds any unit norm in at most fifteen expected trials. [Halving uniform squares on the single circle](single-quadratic-count.md#halving-uniform-squares-on-the-fixed-circle) already suffices for general factoring, including repeated prime factors. On semiprimes, [unique odd-order halving](single-quadratic-count.md#forcing-a-split-from-a-unique-odd-order-half) suffices via randomized conic selection and negative-label planting. A [coverage bound for affine exponents](single-quadratic-count.md#coverage-of-affine-exponents-for-unique-halving) rules out polynomial-size shared menus of such powering formulas on typical balanced semiprimes. A [joint real-bin and low-bit bound](single-quadratic-count.md#combining-an-archimedean-bin-with-the-low-product-bits) limits short-list prediction even when the bin partition depends on the observed product residue. For the 97-bit example, an [exact projection certificate](check_sum_bit_projection.py) shows all 1,398,104 admissible 26-bit endings still have a lift passing the odd-prime sieve through 31 within the search interval. [Counting proof](single-quadratic-count.md#projecting-the-odd-prime-sum-sieve-onto-low-bits) |
| Shared low bits tested from the product alone | Under p<q≤2p, the [shared-suffix test](check_shared_suffix_split.py) chooses the first power 2^k whose progression width meets the rank-9 lattice certificate. If p≡q mod 2^k then q²≡N mod 2^k. Lift its at most four odd square roots from N and search those progressions for q. A missed split excludes the hypothesized suffix under the balance promise. In 100 planted 48-bit-prime cases the primes shared 32 bits while retaining separated sizes. The test used 25–26 bits computed from N and split all 100, with zero overlap with the fixed Lehman control. It split zero of the 100 original unplanted cases. | This standard modular-square-root and partial-factor combination handles a specific relation between p and q. No general shared-suffix property has been established. [Sparse signed-gap audit](check_sparse_suffix.py): choose e=1 for N≡1 mod 4 and e=−1 otherwise, so δ=q−ep is divisible by four. For a proposed δ modulo 2^k, solve (q−δ/2)²≡eN+(δ/2)²; the right side is odd, giving at most four progressions. Enumerating at most two signed powers costs O(k²) gap proposals. Factor-assisted scoring of the 100 unplanted controls at each size found one match at 32-bit p and zero at 48-bit p. For the matching modulus an N-only replay found a proper factor after 330 distinct proposals. At 48 bits the minimum observed weight was five. Treat this as a special-case baseline. In the [local null model](check_sparse_suffix_null.py), x is uniform odd modulo 2^k and the partner is N/x. The gap x−eN/x is exactly uniform on δ≡1−eN mod 8 for k≥3: each such δ gives four odd roots of (x−δ/2)²=eN+(δ/2)² modulo 2^k, since its right side is 1 mod 8. The expected two-term matches on the benchmark moduli were 1.28 and 0.0063, versus 1 and 0 observed. This comparison gives no evidence of enrichment over local guessing. |
| Curve coefficients and finite-state paths | [Higher-genus coefficient matrices](cartier-paths.md) reverse the two hidden-prime factors in the two reductions. These paths read hidden prime labels. They give no digit automaton on N. Truncating a permutation to half its states can make one order survive and the other vanish. Three small curve degrees split 73/110 toy balanced inputs with exact coefficients. | A public cubic powering control obeyed the order identity but split 0/200 larger controls via entries or determinant. For full-matrix rank the [minimum local collision degree](cartier-paths.md#minimum-degree-for-a-local-rank-drop-in-power-substitution) is max(2,ord_q(p)) at p, excluding all small squarefree quotients on those controls. Find a fast recurrence also producing a one-sided zero. In this example the weights cost large-index binomial evaluations. A [modified adjacent-binomial family](cartier-paths.md#guaranteed-witness-menu-from-adjacent-binomial-carries) guarantees a splitting witness or elementary preprocessing exit for every distinct balanced semiprime in O(log² N) public tests; evaluating the binomials remains open at polynomial bit cost. For mixed residues modulo four, the [neighbor of the Gaussian CM middle coefficient](cartier-paths.md#the-coefficient-beside-the-gaussian-cm-middle-term) is the d=4 specialization and singles out the prime congruent to one. [Formal multiplication](cartier-paths.md#formal-multiplication-on-the-gaussian-cm-curve) first becomes nonzero at that prime index; detecting it by explicit truncation requires power-size precision. [Scalar evaluations at multiples of N](cartier-paths.md#scalar-evaluations-at-multiples-of-the-modulus) admit a public logarithm and suppress the degree-p term through precision N^p. [Multiplication by N on actual points](cartier-paths.md#multiplication-on-actual-points-in-degrees-one-and-two) is bijective through extension degree two on balanced inputs; degree three can have a one-sided kernel. [Small-extension coverage](cartier-paths.md#how-often-small-extensions-can-contain-n-torsion) is negligible on typical balanced semiprimes even allowing adaptive curve choice with rational two-torsion. The [trace-one exception](cartier-paths.md#the-trace-one-exception-already-has-an-operational-cm-splitter) was already implemented by the Cheng track; even H widely spaced discriminants cover only O(log B · √(H/B)) of prime pairs at scale B. Given the [actual fundamental discriminant](cartier-paths.md#what-an-input-specific-trace-one-cm-field-would-reveal), short square-factor enumeration already recovers a typical prime without constructing a class polynomial. A [joint square-class hint](cartier-paths.md#joint-square-class-gives-a-uniform-balanced-decoder) for both factors gives a two-trial decoder for every balanced pair. With variable shifts, [binary local signs with up to one eighth corrupted](cartier-paths.md#a-noisy-sign-decoder-for-a-family-of-shifted-products) suffice for CRT recovery without balance. A [held-out small-prime valuation check](cartier-paths.md#do-cheap-coefficient-valuations-predict-the-path-type) scored 50.17% on the balanced degree-five path label. |
| Split versus inert quadratic algebra | Jacobi −1 forces opposite splitting. [Rank tests](single-quadratic-count.md#which-ranks-can-see-a-split-quadratic-algebra) invariant under algebra isomorphism and scalar extension have the same local ranks as a split control. For marked polynomial matrices of entry degree d and equal generic rank r, unequal local ranks occur with probability at most rd/(p−1)+rd/(q−1) over a Jacobi-negative parameter. | [Scalar-square decision](single-quadratic-count.md#scalar-squares-in-the-quadratic-algebra-and-the-hidden-character) does distinguish the split-prime character and is decision-equivalent to ordinary quadratic residuosity on Jacobi-positive inputs. Find a ground-field-sensitive operation of polynomial bit cost; public compact high-degree functions remain outside this bound. |
| Primality-test residues | The [October analysis](../../notes/primality_analysis.md) already compares test rejection with asymmetric coefficient divisibility. On the 200 larger missed balanced controls, [small AKS-style checks](check_aks_window_controls.py) rejected every input. Neither individual coefficient gcds nor scalar Fermat gcds split any, for shifts 1,2,3 or the floored second-quarter endpoints and midpoint proxy ⌊3N/8⌋. Cyclic degrees were 3,5,7,11; this was a diagnostic, not a full AKS run. | Use a necessary property of factor pairs to eliminate candidates even when no gcd splits immediately. |
| Complex reciprocal phase | A [normalized Poisson kernel](analytic-divisor-window.md) equals one on divisors and at most 1/5 elsewhere. Raising it to O(log W) gives an interval sum accurate to 1/4 for W integer positions. | A direct contour shift crosses Θ(√N) poles in a balanced window. Grouping a related kernel by product error bounds leakage below 1/4 but leaves the target count as its central term. A fast collective evaluation remains open. |
| Shared linear basis for the divisor kernel | [Rank bound](analytic-divisor-window.md#rank-required-by-a-shared-linear-expansion) of at least B on (B,2B] under total absolute pointwise error below one, even with arbitrary input coefficients. | Applies to pointwise compression. For one hidden divisor, all interval totals admit a [polynomial-size shared basis](analytic-divisor-window.md#a-shared-basis-for-single-divisor-interval-totals). Computing its row from N remains open. |
| Prefix descent | At most half the bit length in parity queries, then lattice recovery. Verified on 676 odd moduli using exhaustive queries. [Lattice backend check](bounded-quadratic-roots.md#lattice-backend-check) [Decoder](check_lower_half_prefix.py) | Evaluate a wide query from N and I. |
| Inverse pairs | Swapping the two coordinates gives exactly the same parity. Separate prime-factor counts lose their overlap. [Identity](single-quadratic-count.md#swap-symmetry-gives-an-exact-interval-parity-target) · [Counterexample](quadratic-chirp.md#separate-local-sign-counts-lose-prefix-parity) | [Endpoint rectangles](single-quadratic-count.md#empty-corners-in-the-inverse-pair-box) can be pruned exactly, but a linear number of interior rows survives. [Explicit rational-center exclusions](single-quadratic-count.md#coverage-of-rational-center-root-exclusions) require Ω(N) distinct centers to cover a fixed fraction. Compute the joint parity without scanning the interval. |
| Lean and CM | The Gaussian-binomial identity sums over subspaces. A one-sided CM divisibility relation with enough lifts would yield a root or factor. [Lean comparison](single-quadratic-count.md#comparison-with-the-finite-subspace-identity-in-grh-lean) · [Extractor](check_auxiliary_relation.py) | Construct the relation from public inputs and prove its divisibility. |
| Tensor contraction and masks | The root predicate has a small circuit. Dense elimination is expensive. A weak mask predictor already suffices on distinct odd semiprimes. [Circuit](quadratic-chirp.md#a-public-tensor-network-for-the-prefix-query) · [Prediction criterion](quadratic-chirp.md#a-zero-test-using-a-quadratic-boolean-phase) | Give a compact contraction or a predictor with a proved advantage at each modulus. |
| Exponents and square lifts | Several public exponent formulas have negligible mean coverage. A sufficiently accurate quotient interval enables lattice recovery. [Exponent bounds](bounded-quadratic-roots.md#fixed-polynomial-exponents-on-separated-balanced-primes) · [Quotient accuracy](bounded-quadratic-roots.md#square-root-cutoff-for-approximate-quotients) | Prove coverage for a public construction outside the tested families. |

<details>
<summary>Full attempt ledger, including individual tests and sources</summary>

| Approach tried | Established result | What would justify reopening it |
|---|---|---|
| ePrint 2026/3020 | The authors report conditional deterministic subexponential discrete logarithms on known prime-field inputs. [Scope and contraction comparison](quadratic-chirp.md#comparison-with-eprint-20263020) | A public O(b²)-gate network counts prefix roots with unique witnesses. [Construction](quadratic-chirp.md#a-public-tensor-network-for-the-prefix-query) Its schoolbook square has a [clique minor on the free input bits](quadratic-chirp.md#width-of-dense-elimination-for-the-schoolbook-square), forcing exponential dense messages on wide queries. On distinct odd semiprimes the final Y=1 indicator has rank ≤4 at every cut. A [public exact tensor-train constructor](quadratic-chirp.md#constructing-the-compressed-root-indicator-on-semiprimes) is equivalent to deterministic splitting under the stated bit bounds. [Approximate trains](quadratic-chirp.md#approximation-in-squared-error-suffices-with-a-margin) also suffice if total squared error beats the two-public-root baseline by an inverse-polynomial margin δ. [Subtracting the public roots](quadratic-chirp.md#prefix-search-after-subtracting-the-public-roots) reduces the prefix budget to O(b/δ). With K explicit Walsh characters the [optimal squared-error improvement](quadratic-chirp.md#exact-walsh-approximation-error-for-the-mixed-root-pair) is min(2,4K/2^b). However [2(b−1) targeted signed sums](quadratic-chirp.md#reading-a-mixed-root-with-weight-two-walsh-queries) of weights one and two suffice to read the odd mixed root. A [weak callable predictor on uniform masks](quadratic-chirp.md#a-weak-random-mask-predictor-suffices) for the parity-selected moment also suffices by Goldreich–Levin decoding. For [single moments](quadratic-chirp.md#predicting-a-single-signed-moment), inverse-polynomial signed correlation on nonzero masks suffices without identifying those masks. With a [low-bit-gated quadratic phase](quadratic-chirp.md#a-zero-test-using-a-quadratic-boolean-phase), predicting zero for the full four-root sum with any inverse-polynomial advantage suffices. |
| Odd-root prefix parity | For any odd N with at least two distinct primes, the number of odd roots of 1 after removing 1 is odd. [One parity query per digit](quadratic-chirp.md#odd-root-prefix-parity-on-general-odd-moduli) deterministically gives a proper factor, including for repeated prime factors. | The [earlier lower-half argument](../ellipse_disk_mod8/native_dyadic.md) removes the odd-root filter. The [aligned-prefix decoder](check_lower_half_prefix.py) maintains odd corrected parity in its current block clipped to [0,(N+1)/2). If the midpoint reaches the cutoff, only the left child can contribute. Otherwise the left child is a complete dyadic block, so a single ordinary prefix query selects the next child after correcting for root 1. The [bounded quadratic solver](bounded-quadratic-roots.md#square-root-cutoff-for-approximate-quotients) finishes once the clipped interval contains at most floor(√N) integers, in at most eight blocks of span at most max(1,floor(√N/8)). Halving from size 2^(b−1) gives at most floor(b/2) prefix queries. On the 676 qualifying inputs below 2000, total oracle calls fell from 6055 to 3073. At N=133, [two aligned prefixes](quadratic-chirp.md#separate-local-sign-counts-lose-prefix-parity) have identical per-prime ±1 counts and opposite joint parity. The [existing swap identity](single-quadratic-count.md#swap-symmetry-gives-an-exact-interval-parity-target) gives the same oracle by counting units in the block whose inverses stay there. End-to-end enumeration also split all 676 test inputs in 3073 calls using this formulation. Computing a wide-block parity without scanning remains open. |
| Quadratic floor sums | The [earlier floor-sum audit](../../notes/floorsum_results.md) already gives the adjacent-offset identity. On each aligned block sum floor((x²−1)/N)−floor((x²−2)/N) modulo two. Direct evaluation reproduced all 676 splits and 3073 queries, including blocks containing zero. | A method that aggregates this fixed-target difference without scanning the block. Interchanging x and the quotient k enumerates the same witnesses x²=1+kN. Hardness results for arbitrary targets do not settle the fixed-target parity problem. |
| Removing a known square scale | Multiplication pulls the target interval back to exactly min(h,N−h,w,N−w) cyclic runs. [Proof and small-scale reduction](single-quadratic-count.md#interval-cost-of-removing-a-known-unit-scale) | For wide intervals and large signed scales, a shortcut must preserve the root count without expanding the entire selector into intervals. |
| Reflection and short exceptional sets | Compare membership of x and N−x on the lower half. If either the agreement or disagreement set has size at most √N poly(log N), the lattice solver evaluates the wide interval exactly. [Criterion and checks](single-quadratic-count.md#wide-intervals-with-a-short-exceptional-set-after-reflection) | The [implemented evaluator](check_lower_half_prefix.py) uses the √N cutoff. It resolved 74,702 of 94,910 small cyclic queries exactly; the other 20,208 returned unresolved. In the 676-input split check, remaining oracle calls fell from 3073 to 2332. At quarter width its [exact coverage over uniform cyclic starts](single-quadratic-count.md#frequency-of-short-reflection-exceptions) is (2 floor(√N)+1)/N for large odd N. A [pointwise combination](single-quadratic-count.md#combining-reflection-evaluable-selectors-at-several-scales) still needs order √N windows for a balanced prefix, even with multiple unit scales. A [polynomial map preserving all roots](single-quadratic-count.md#two-evaluations-classify-polynomial-maps-on-the-roots-of-one) either acts as a constant ±1, identity or negation, or yields a factor at ±1. The [quadratic prefix product](single-quadratic-count.md#quadratic-products-lose-prefix-parity-at-bounded-precision) loses parity even with polynomially many low coefficients modulo N^k at polynomial k. For the full marked product, the first surviving degree modulo a distinct odd semiprime is twice its smallest prime factor. [One full-product evaluation](single-quadratic-count.md#a-single-value-of-the-full-product-splits-a-semiprime) at a Jacobi-selected point returns a factor by gcd. [Prime-power lifting](single-quadratic-count.md#full-product-evaluation-given-an-arbitrary-factorization) gives the reverse reduction on all moduli. [Public normalization](single-quadratic-count.md#normalizing-the-product-gives-the-local-quadratic-characters) leaves the CRT vector of local Legendre signs. One exponent implements it on all units exactly when the prime divisors have matching v₂(p−1). For a sparse polynomial or rational sign formula, [the exponents alone suffice to factor even with error below the reciprocal of the total term count](single-quadratic-count.md#the-exponents-of-a-sparse-sign-polynomial-suffice-to-factor). On distinct odd semiprimes, [weak random-mask prediction](single-quadratic-count.md#weak-predictions-on-jacobi-negative-inputs) of the even mixed root also suffices by Goldreich–Levin decoding. The [ordinary low bit](single-quadratic-count.md#the-low-bit-selects-one-hidden-prime-character) gives one hidden prime’s Legendre character. Tracking which character is selected at successive powers of two supplies the earlier binary-prefix decoder. Its [exact unit-shift correlation](single-quadratic-count.md#exact-shifted-low-bit-correlation-returns-the-complementary-factor) gives the complementary prime; independent sign products have bias of order 1/p. The [unsigned support count](single-quadratic-count.md#counting-consecutive-units-already-contains-the-factors) is affine in the totient for semiprimes; logarithmic-length Schemmel totient queries recover arbitrary factorizations by interpolation. At even lengths, parity alone locates the least prime divisor by binary search; the fixed-point test reduces to factorial coprimality. Efficient public evaluation and an aggregate parity identity remain open. |
| Lean subspace identity | The Gaussian-binomial sum collapses to a kernel dimension. Our quadratic fiber fails the affine midpoint test; a binary prefix can contain exactly three roots. [Comparison](single-quadratic-count.md#comparison-with-the-finite-subspace-identity-in-grh-lean) | A public small kernel testing whether any root has prefix 0000 suffices on distinct odd semiprimes. [Existence reduction](single-quadratic-count.md#existence-alone-on-prefix-0000-suffices-for-a-semiprime-split) · [Same-modulus cost](single-quadratic-count.md#choosing-the-zero-prefix-length-with-more-prime-factors). With auxiliary moduli and variable prefix length, O(log N) decisions recover a root for arbitrary odd N coprime to 3. [Digit reduction](single-quadratic-count.md#recovering-a-root-by-changing-the-query-modulus) |
| Divisibility versus size; proposed CM valuation | Chebyshev's binomial argument bounds a product of primes by a binomial coefficient. The Howgrave–Graham step forces an integer auxiliary value to vanish when its absolute value is smaller than its known divisor. The Lean kernel identity above is an exact subspace sum. | CM proposal remains unspecified. To test it, give a public auxiliary family, a proved divisor at the unknown root, and a coefficient-size bound on the proposed interval. Compare the certified interval at polynomial construction cost. In [the one-curve exponent](https://github.com/lsernau/grh-lean/blob/72fead102f663d5c7a407964a13dc0cf9d977207/grh/src/New/Derivatives/FormalGroup/Exponent.lean), leastFlatDegree is p+1 for inert p and 2 for the ramified p=d. The inert coefficient at p+1 is a unit in the normalized correction series. This is a support threshold for that component. Choosing d to hit the ramified case makes gcd(d,N) available directly. The [derivative divisibility bound](https://github.com/lsernau/grh-lean/blob/72fead102f663d5c7a407964a13dc0cf9d977207/grh/src/New/Derivatives/Vacuum/Divisibility.lean) concerns scaled jets of CM standard data; transferring it to integer auxiliary values at modular roots remains open. For f(T)=(A+T)²−y and known N, first check whether the normalized candidate belongs to ℤ[T,f(T)/N]. If so, [CHHS, Theorem 2](https://arxiv.org/pdf/1605.08065) already precludes a nonzero candidate uniformly smaller than one on a real interval centered at zero of radius greater than 2√N. For odd N, unit squares, and min(m,1+deg H) at most the smallest prime factor, an integer polynomial divisible by N^m on every integer lift necessarily belongs to (N,f)^m. Even one opposite sign pair suffices given s=min(m,1+deg H) consecutive lifts in each class. At any multiplicity, an explicit H outside the standard ideal with sufficient lifts of just one sign yields a root or factor by [coefficient extraction](bounded-quadratic-roots.md#extracting-a-one-sided-relation-at-arbitrary-multiplicity). [Circuit reduction](bounded-quadratic-roots.md#evaluating-a-compressed-auxiliary-polynomial) covers compact division-free H at polynomial multiplicity, even with exponential degree. [Division gates](bounded-quadratic-roots.md#division-gates-and-the-quadratic-norm) admit unit inversion; a nonunit with nonzero base residue gives a root or factor. [Extractor](check_auxiliary_relation.py) checked on 120 planted root relations, 120 planted factor relations and 48 standard-ideal controls. [Proof](bounded-quadratic-roots.md#two-opposite-root-classes-already-suffice) For standard polynomial-degree integer shifts even an integer-grid bound has a [finite-difference obstruction](../ellipse_aggregate/native_cap.md#small-root-certificates-on-the-whole-disk). |
| Public power denominator | For D=(X+c)^(N−1)−1, X²=y and uniform c, the norm hits either prime of N=pq with probability at most 2 gcd(p−1,q−1)(1/p+1/q). [Analysis and checks](bounded-quadratic-roots.md#a-public-power-denominator) | For separated balanced prime ranges the [mean hit rate](bounded-quadratic-roots.md#average-coverage-for-separated-balanced-primes) is O(log³ N/√N). On uniform unit squares a polynomial-size list of predetermined shifts also has negligible mean coverage. Even [target-dependent rational shifts](bounded-quadratic-roots.md#target-dependent-rational-shifts) of polylogarithmic degree have negligible average coverage under the stated denominator hypothesis. [Monomial shifts](bounded-quadratic-roots.md#monomial-shifts-of-arbitrary-degree) have negligible mean coverage at arbitrary degree by counting distinct solution ratios. For an [affine power shift](bounded-quadratic-roots.md#the-exceptional-value-for-an-affine-power-shift), gcd((a y^k)²−y,N) captures all exceptional hits; a=1 maximizes its root-or-factor success on uniform unit squares. A [bounded missing cofactor](bounded-quadratic-roots.md#completing-an-exponent-that-covers-one-hidden-prime) at one prime gives a direct split with probability at least 1/2 after exponent completion. A [fixed polynomial exponent](bounded-quadratic-roots.md#fixed-polynomial-exponents-on-separated-balanced-primes) has mean scalar hit rate at most N^(−1/2+o(1)); the linked earlier worst-case result covers larger length-indexed menus. The [rounded square-root estimate for φ(N)](bounded-quadratic-roots.md#estimating-the-exponent-from-the-square-root-of-n) also has mean coverage at most N^(−1/2+o(1)) on separated balanced primes. For primes drawn from the [same interval](bounded-quadratic-roots.md#balanced-primes-in-the-same-interval), the best mean exponent among these offsets is N^(−1/4+o(1)), attained by offset −1 for group-order coverage. For that offset, [removing exact totient hits](bounded-quadratic-roots.md#removing-the-exact-totient-hits) leaves mean coverage N^(−1/2+o(1)). Constructing a general exponent remains open. |
| Ternary decoder with the existing short-root solver | For a surviving class up to sign modulo D, write x=r+Dt for each sign representative r. If floor((B-1)/D) is within the [lattice solver's bound](bounded-quadratic-roots.md), its monic quadratic routine recovers the root directly for any unit square y by dividing (r+Dt)²−y by D² modulo N. Here B is the decoder's power-of-two cutoff and D is coprime to N. | The [decoder](check_prefix_root_recovery.py) now uses this cutoff. On all 181 unit squares across its six test moduli, total prefix queries fell from 1760 to 1309 with the lattice cutoff, then to 1152 by inferring the third child after two negative answers. The three signed children partition the occupied parent, so each refinement needs at most two queries. The first wide prefix query still needs an evaluator. Revisit if a method covers those wide intervals without a factorization. |
| Binary root recovery with even query moduli | For any odd N>1 and any square y, a root below B exists in ±1 mod 4: of h and N−h one is odd (take 0≤h<N). For an occupied signed class ±a mod D, D a power of two at least 4, its two children are ±a and ±(a+D) mod 2D. Test the first by imposing x²=a² mod 4D alongside x²=y mod N using CRT and the same cutoff B. For odd x,a one of x−a,x+a is exactly 2 mod 4, so their product vanishes mod 4D exactly when the other is divisible by 2D. A negative answer selects the second child. | One query per refinement, followed by the existing short-root solver. Requires an evaluator on even moduli 4DN. [Binary decoder](check_prefix_root_recovery.py) verified on all 6176 square residues for odd N from 3 through 255, including nonunits. The wide-query evaluator remains open. |
| Probabilistic binary branch choice | Inverse-polynomial total probability on occupied leaves suffices for verified randomized recovery. A single occupied path of inverse-polynomial probability also permits deterministic search by pruning low-probability prefixes. [Bounds](single-quadratic-count.md#a-probability-bound-sufficient-for-binary-root-recovery) | Square blinding reduces the unit case to mean success at each fixed N. A K-valued target summary gives mean success at most O(K/N^(3/5)) on distinct odd semiprimes with the current cutoff. [Counting argument](single-quadratic-count.md#how-much-target-information-must-reach-the-leaf-selector) A sufficiently accurate public predictor remains open. |
| Skipping initial ternary digits using extra representatives | At N=77 and cutoff B=128, exhaustive enumeration of 0≤x<128 shows that the classes ±1, ±2, ±4 modulo 9 respectively miss the unit squares 4, 16, 9. Thus no first child chosen from N alone always works. | Enlarging the cutoff to DN makes any prescribed class modulo D=3ʲ feasible by CRT, but leaves exactly N candidates in that progression. This shortcut alone gives no interval shrinkage. Even (N, y mod 9) is insufficient. For N=205, B=256, targets 4 and 139 agree modulo 9. Their complete root lists below B are {2,43,162,203,207,248} and {37,78,127,168,242}, giving disjoint admissible branches {±2,±4} and {±1}. A choice using the full y remains open. |
| Clipped inverse involution | Exact binary nullity formula. Signed edges give one kernel basis vector per root over F₃, still one column per candidate. [Signed construction](single-quadratic-count.md#removing-two-cycle-kernel-vectors-for-the-existence-test) The unsigned binary version has one column per interval element; squaring recovers the prior inverse-pair count. [Derivation](single-quadratic-count.md#kernel-dimension-for-a-clipped-inverse-involution) | A query-dependent kernel with efficiently computable entries. A two-dimensional cosine determinant gives the fixed-cutoff sign, but its determinant has rational degree φ(N)/4. [Semiprime determinant](single-quadratic-count.md#a-two-dimensional-quarter-determinant-over-a-large-number-field) · [General odd case](single-quadratic-count.md#extending-the-quarter-determinant-to-any-odd-modulus) · [Recovering the pair count from traces](single-quadratic-count.md#recovering-the-missing-zeroth-trace-from-a-pair-count-bound). Shared binary dot products across all unit squares and half-window translations require full column rank. [Rank proof](single-quadratic-count.md#rank-of-the-translated-window-answers-across-square-inputs). Per-query random error permits a short vector; computing its bucket parities is open. [Randomized summary](single-quadratic-count.md#randomized-linear-summaries-for-translated-windows) |
| Window and modulus symmetries | Binary prefixes handle the cutoff. Translation solves unit-scale quarter parity when four divides the modulus; doubling an odd modulus preserves its hard bit. [Results](single-quadratic-count.md#binary-prefixes-isolate-the-interval-condition) | For distinct odd semiprimes it suffices to evaluate root parity on the single prefix 000. [Three-bit reduction](single-quadratic-count.md#a-single-public-binary-prefix-suffices-on-distinct-odd-semiprimes) A balanced split of the free bits has full rank for the pointwise phase; a single requested coefficient has rank at most two. [Split proof](single-quadratic-count.md#splitting-the-free-bits-of-the-000-prefix) Computing its parity remains open. |
| Salié sums and prediction | Exact and noisy decoders exist. A square-only sign predictor with inverse-polynomial advantage suffices on distinct odd semiprimes. Public direct summation still scans the units. [Sign reduction](salie-quarter.md#a-square-only-sign-predictor-needs-any-nonnegligible-advantage) | An estimator with a runtime or advantage guarantee on unfactored inputs. |
| Public features | Interval tables, inverse bins, sparse phase fits and fixed shifted Jacobi probes have model-specific bounds. Held-out square-only scoring needs no factors. An inverse-polynomial advantage on the sampled-root label suffices for a split at any odd modulus with at least two distinct prime factors. [Reduction](salie-quarter.md#distinct-prime-dependence-of-two-sign-reductions). The advantage itself is impossible on families with ω(N) growing faster than log log N. [Information limit](salie-quarter.md#information-limit-for-the-sampled-root-label-at-many-prime-factors) [Scoring method](salie-quarter.md#public-labels-for-square-only-learning-and-held-out-scoring) | A feature outside the recorded bound with reproducible held-out advantage. |
| Individual square witnesses | The best fixed L-level list for distinct odd semiprimes has hit rate Θ(√(L/N)) for 1≤L≤N. Choosing among K such lists gives at most O(√(KL/N)) for KL≤N. Affine selectors and bounded-degree polynomial lift menus have negligible coverage under the stated size bounds. [Coverage proofs](salie-quarter.md#coverage-of-searching-integer-square-lifts) | A selector outside those families with a coverage argument. |
| Approximate square-lift quotient | A nonzero root has a representative in the upper half by changing sign if needed. Given a predicted interval [a,b] containing its quotient k=(x²−y)/N, restrict x to the upper half and to y+aN≤x²≤y+bN. The integer endpoints differ by at most b−a since their sum is at least N. The [shifted quadratic solver](bounded_quadratic_roots.py) therefore suffices for b−a≤X with the current X of order N^(2/5). Checked containment on all 5858 upper-half representatives for N=35,77,205,1001,10403 using intervals supplied from the true quotient. | The quotient routine now splits wider intervals into square-root-scale blocks and adjusts lattice rank per block. The [variable-rank bound](bounded-quadratic-roots.md#square-root-cutoff-for-approximate-quotients) permits width √N times a fixed polylogarithmic factor in polynomial bit time. Open: predict such an interval from (N,y) with inverse-polylogarithmic coverage on uniform unit squares at each modulus. |

</details>

For each proposed experiment, check the closest row and its linked proof first. Record the changed assumption, the result, and what remains open in that row. Add a row when the mathematical mechanism changes. Detailed proofs and runs belong in the linked notes.

## Priorities

1. **General input: evaluate one clipped root count.** Use the parity of mixed roots of $u^2=1\pmod N$ in a translated window of length $(N-1)/2$. It suffices for deterministic bisection on any odd modulus with at least two distinct prime divisors. The same task is realized by signed disk counts modulo 4 at one cutoff and one discriminant. The forms are powers of a public class of known order, and query exponents take polynomial bit time to construct. The simple one-row formulation is the smaller test of a new evaluator. The class-group formulation supplies an alternative if the evaluator can exploit the order structure. [Fixed window](fixed-window-split.md) · [Common discriminant](fixed-discriminant-split.md) · [Cyclic class encoding](cyclic-class-queries.md)
2. **Balanced input: retain the fixed absolute moment as a control.** The implemented hybrid uses roots modulo the public discriminant and tangent batching. It computes the exact scalar with expected $N^{2/5+o(1)}$ bit cost. The 36-bit control has trace 69,529,082 and divisor 254,119. The bound concerns the scalar and remains larger than the complexity of established general factoring algorithms. Further geometric tuning is a reference-evaluator track. [Progression evaluator](discriminant-progressions.md) · [Tangent hybrid](tangent-progressions.md)
3. **Test partial information before requiring exact output.** The quarter-window reduction splits inputs with exactly two distinct prime divisors without adaptive queries. A lattice finishes from a bounded interval or progression, so a fixed-precision evaluator can target just enough bits to narrow a root to $N^{1/2-\varepsilon}$ candidates. Uniformly chosen windows of this size hit mixed roots with exponentially small probability. A candidate must explain how it obtains the partial coordinate from the public input. [Quarter queries](nonadaptive-quarter.md) · [Small-root decoder](bounded-quadratic-roots.md)

A direct lattice-block evaluator now computes central shifted-window parity in $\sqrt N\operatorname{poly}(\log N)$ bit operations. It improves the direct reference exponent for this scalar and still covers a power-of-$N$ number of blocks. [Lattice block evaluator](lattice-window.md)

For any number of distinct odd prime divisors, nonadaptive affine isolation takes $O(\log N)$ weighted queries, or $O((\log N)^3)$ unweighted modulo-4 queries, with constant success probability. These are alternatives when a candidate evaluator can process a batch more cheaply than adaptive queries. [Affine isolation](affine-pair-bound.md) · [Unweighted reconstruction](affine-parity.md)

On inputs with two distinct prime divisors, an inverse-polynomial prediction advantage on the scaled quarter-window bit also suffices. Fourier localization replaces exact bit reconstruction in that setting. [Weak prediction reduction](quarter-fourier.md)

## Approximate counts and analytic targets

The general odd-input reduction also works with a one-bit gap test on ordinary quadratic counts in a fixed-dimensional domain. The inverse phase modulo $N^2$ gives exact periods for all nonstationary blocks. A lower-input interval isolates a mixed root, and a second affine interval keeps the nonstationary variation small. On a good trial, the gap bit has a constant-size Fourier coefficient at the hidden frequency. Repeated prime factors are included. [Input isolation](square-input-isolation.md) · [Affine interval filter](affine-interval-variance.md) · [One-bit comparison](count-comparison.md) · [Single quadratic count](single-quadratic-count.md)

For inputs with exactly two prime divisors, a constant-precision normalized Salié or square-modulus Kloosterman evaluation would also suffice. One signed phase-interval sum per scale preserves the hidden frequency. Single-term sampling has high variance both before and after grouping by phase, so this route needs aggregation across distinct phase values. [Salié sum](salie-quarter.md) · [Square-modulus sum](square-kloosterman.md) · [Signed interval](salie-phase-bins.md) · [Grouped variance](salie-grouped-variance.md)

The Fourier decoder reuses one sample pool across its interval search. It recovered planted noisy signals on inputs up to 36 bits. The public predictor tests found no clear signal for the specified adaptive exponent rules or for the carry and Euclidean feature models on held-out moduli. [Decoder and calibration](quarter-fourier.md) · [Exponent screen](adaptive-exponent-screen.md) · [Held-out feature screen](carry-feature-screen.md)

## Native evaluator audit for the complete sums

Distinguish the summation modulus $M$ from the integer $N$ to factor. The older repo studies pure cubic Weyl sums $\sum_xe_N(ax^3)$ at ordinary modulus $N$. Its nine-value reduction for split semiprimes assumes access to approximate sum values. The small-scale target here uses the fixed phases $P_3(x)=3x-x^3$ and $P_4(x)=x^4-2x^2$ at $M=N^2$. Classical stationary formulas give direct local evaluations for known prime powers. Hiary's fast cubic method for shorter, small-cubic-coefficient sums does not cover the coefficients at the full period length. [Older cubic reduction](../../notes/analytic_results.md#4-reduction-c3-complete-cubic-weyl-sums-mod-n) · [Prime-power stationary phase](https://doi.org/10.4064/aa-80-1-1-48) · [Hiary's cubic method](https://annals.math.princeton.edu/wp-content/uploads/annals-v174-n2-p04-p.pdf)

For an arbitrary odd square $M=B^2$ and unit scale $t$, writing $x=a+By$ gives a direct evaluation from $M$ alone. Sum complete upper-digit blocks first. Exactly those $a\bmod B$ with $P'(a)=0\pmod B$ contribute, each with normalized phase $e_M(tP(a))$. Scan these derivative residues and cache a histogram of exact phase exponents. No factorization or lifting is needed. For the cubic simplified derivative test remove the public prime 3. At ordinary nonsquare composite moduli the implementation uses literal symmetry-paired summation. [Native phase evaluator](native_phase_evaluator.py)

If $M=p^k$ with $k\ge2$, exact perfect-power recognition and primality give the prime publicly. For even $k$ the local critical points $\pm1$ for cubic and $0,\pm1$ for quartic directly give the normalized terms. For odd $k=2r+1$ and $r\ge1$, Taylor expansion in a stationary class has only the quadratic term after the constant, with coefficient $p^{2r}tP''(a)/2$. The remaining local Gauss factor after normalization is $\epsilon_p(tP''(a)/2\mid p)$. The critical roots and their values remain public. This is a polynomial-bit evaluator on these nondegenerate perfect-power branches. The general odd square case still uses $B$ derivative checks; the lattice cache above retains a better worst-case exponent.

The evaluator matched 60 literal sums across ordinary moduli, square moduli, and prime powers. It used two cached phases for both polynomials on a public 90-bit prime power. On $N=11009$, $M=N^2$, the derivative scan used 11,009 exact checks and produced four weighted phases in approximately 0.001–0.002 seconds of setup. The earlier Python lattice cache on the same control took about 1.21 seconds on this run. The scan also prepared four phases at $N=471311$ in about 0.068 seconds. These are measured small-input comparisons. [Evaluator and tests](native_phase_evaluator.py)

On inputs with two prime-power components, averages of the magnitude and squared magnitude over long scale windows now have native polynomial-bit evaluators. The windows can lie entirely within the small-scale range when they have length $N/\operatorname{poly}(\log N)$ or more. An unusually low mixed line is found by an ordinary quadratic equation and gives a factor. Otherwise a public diagonal sum is accurate to inverse-polynomial additive error on those windows. [Window evaluator and complexity](coarse-phase-window.md)

Multiplying the same window by an exponential at any supplied frequency gives a native test of a translated band. Certified quadratic lattice blocks recover a factor if that band contains a mixed line. On balanced primes, even all rational centers of polynomially bounded denominator cover only an exponentially small fraction of inputs. Selecting a productive center from the arithmetic of $N$ remains open. [Modulated window and coverage bound](modulated-phase-window.md)

Scales $t_0+Nj$ fold those bands to modulus $N$. A native window now needs only $\sqrt N/\operatorname{poly}(\log N)$ progression indices for polynomial-bit evaluation, with scales outside the earlier $2N$ range. Complete folded coverage still costs $\Theta(\sqrt N)$ lattice blocks. The center audit records the older geometric shortlist, predictor screens, and scope of arithmetic centers from near-rational factor ratios. [Bandwidth and center audit](bandwidth-center-audit.md)

Complete quadratic weighting gives a separate native family when the outer quadratic coefficient is a multiple of $N$. At inverse-polynomial accuracy it evaluates the actual all-scale raw cubic or quartic on squarefree semiprimes, or finds a factor. Unit-coefficient weighted sums can preserve the point magnitude. Five centered literal squared-quartic chirps recover that magnitude even with repeated and imbalanced components. [Quadratic weighting](quadratic-chirp.md)

### Variance after subtracting the public quartic terms

For distinct odd primes $N=pq$, subtract $1+4e_M(-t)$ from the normalized quartic using the Gaussian completion below. The remaining contribution $R(t)$ is the average of $Z_s=e_M(-s^2/(4t)-s-t)(g_M(s)-N)$ over uniform $s\bmod M$, with $Z_s=0$ outside $1<\gcd(s,N)<N$. Its magnitude is at most four by the CRT product. At gcd $p$, the quadratic Gauss factor is $N\sqrt p$ times a signed fourth root of unity, with signs balanced. At gcd $p^2$ it is $Np$. The corresponding probabilities are $(p-1)(1-1/q)/p^2$ and $(1-1/q)/p^2$. The symmetric $q$ cases give the exact second moment
$$\mathbb E|Z_s|^2=4N^2(1-1/p)(1-1/q).$$
Thus independent uniform summand sampling needs order $N^2$ samples for constant squared error on balanced semiprimes. Optimal one-term importance sampling still has second moment at least $(\mathbb E|Z_s|)^2=\Omega(N^{3/2})$ there. This bound follows from Cauchy–Schwarz and the gcd-$p$ terms alone. Both bounds apply at small unit scales because the scale changes only the phase of each summand. They concern independent sampling after subtracting the public terms. Analytic aggregation across the remaining indices is open. The formulas matched five direct controls. [Sampling diagnostic](quartic_sampling_check.py)

## Degree-eight benchmark and one reciprocal quotient range

The smooth upper-digit moment uses $P_8(t)=B_8(t)/B_8(0)=1-20h^2-40h^3-30h^4$, where $h=t(1-t)$. Average $P_8(z/N)$ over the lifts of each unit lower digit $a$ in the isolation window and subtract the two public stationary blocks. For $d=\gcd(a^2-c_0^2,N)$ the upper digits run through $s+d\mathbb Z$ with multiplicity $d$. Bernoulli multiplication gives the exact block average $(d/N)^8P_8(s/d)$. The remaining stationary block on a singleton trial contributes $P_8(((z_0+2uj)\bmod N)/N)$. [Benchmark](bernoulli_moment.py) · [Isolation](affine-pair-bound.md)

The local sum of eighth-power weights at $p^e$ is at most $2/(1-p^{-7})$ by the valuation counts in the input-isolation proof. Write $Q=2^{\omega(N)}$. Subtract the $Q$ stationary weights of one and average over the translated window. The expected supremum of the remaining signal is at most $(wQ/N)/1457$, using $\sum_{p\ge3}p^{-7}\le1/1458$. At the isolation width $wQ/N\le 1/8192$, Markov's inequality leaves a constant probability of a singleton with remainder at most $1/64$. A singleton in a width-one window has zero remainder. The Fourier series of $P_8$ gives $\|P_8-\cos(2\pi\cdot)\|_\infty<1/96$. With a $1/64$ comparison gap about zero, the bit can disagree with the cosine sign only on $|\cos|<1/24$. This union of two arcs occupies at most $N/24+2$ digits. The clean sign has coefficient $2\cos(\pi/(2N))/(N\sin(\pi/N))$ at frequency $2u$. The corrupted coefficient remains above $1/4$ for $N\ge15$. The query uses a fixed number of quadratic counts because a $z^k$ weight can be counted by adding $k$ variables in $[0,z)$. [Local valuation counts](square-input-isolation.md) · [Fourier recovery](quarter-fourier.md)

For the reciprocal benchmark let $M=N^2$, $A=(c_0+Nj)^2\bmod M$, and $I=[L,H)$ be one interval of centered denominators with the same nearest quotient $q$. Inverse reciprocity rewrites
$$e_M(A/x)=e_{|x|}(-\operatorname{sgn}(x)A\overline M)\,e^{2\pi i A/(Mx)}.$$
Omitting the smooth correction over the centered range costs $O(\log N/N)$ in the normalized sum. Within the bin, $r=M-q|x|$ runs in a progression and a second inverse reciprocity step replaces the reciprocal modulus by $|r|$ when $q$ is a unit. Recursing on each remaining modulus individually retains the cost of visiting its denominator. [Reciprocity checks and evaluator](reciprocal_range.py)

The prototype evaluates the whole bin at once in the original coordinate. For each unit lower digit $a\in W$, put $b(a)=1-c_0^2/a^2\bmod N$ and $Y_a=\{y:L\le a+Ny<H\}$. The normalized $j$-signal has coefficient
$$ e_M(a+c_0^2/a)\,\frac1N\sum_{y\in Y_a} e_N(b(a)y) $$
at the distinct frequency $2c_0/a\bmod N$. This gives the exact scan of $w$ lower digits and removes the upper-lift enumeration. For $b\ne0$, the finite geometric sum divided by $N$ has magnitude at most $1/(2|b|_N)$. Each derivative value has at most $Q$ preimages since it specifies one unit square. By orthogonality in $j$, truncating to $|b|_N\le K$ thus has mean squared error at most $wQ/(2NK)$ over uniform translations of $W$ and uniform $j$, at the actual bin endpoints. The stored spectrum serves every query digit. A fixed macroscopic interval retains the stationary coefficient with magnitude $|Y_a|/N=|I|/M+O(1/N)$.

The truncated evaluator generates the retained digits by solving the $2K+1$ public congruences $(1-b)a^2=c_0^2\pmod N$, skipping nonunit coefficients. The same certified windows of size $\Theta(\sqrt N)$ used for mixed-root parity solve these general monic quadratics. Its preprocessing costs $O((K+1)(1+w/\sqrt N)\operatorname{poly}(\log N))$ bit operations, including empty blocks. Each subsequent digit uses only the retained coefficients. Fixed accuracy permits a fixed $K$ at the isolation width. Central windows still require exponentially many blocks in the bit length. [Certified block solver](lattice-window.md)

The degree-eight control recovered a verified divisor on six synthetic singleton windows, including repeated factors. On the $N=1575$ fixture it used 22 progression weights in place of 34,650 direct lift weights per query. The bin diagnostic matched 90 literal sums against both reciprocity steps, checked the spectrum by Parseval, and compared the truncated solver with a full lower-digit scan. At $N=471311$ and width 3682, 688 lattice blocks retained one coefficient, with RMS truncation error $0.004342$. These figures measure work units and approximation error on controls. [Moment benchmark](bernoulli_moment.py) · [Bin benchmark](reciprocal_range.py)

On the 36-bit control $N=38942974393$, the positive $q=2$ bin at start 12391793379 and width 131072 used 96 lattice blocks and retained one coefficient at $K=8$. The lower-digit scan took about 1.09 seconds and the truncated evaluation about 0.10 seconds in a run with FLINT. The normalized RMS truncation error was $5.61\times 10^{-6}$. The fixture uses known factors to place a singleton window. The evaluators take only the public window parameters. SymPy's Python LLL backend failed an internal assertion on one basis in this run. To use the FLINT backend, run the following from this directory.

```sh
uv run --no-project --with sympy --with python-flint==0.9.0 python reciprocal_range.py 38942974393 --start 12391793379 --width 131072
```

## Recovering a factor from the full stationary magnitude

Let $M=N^2=\prod_i q_i$, where $q_i$ are the complete prime-power components, and let $e_i$ be the CRT idempotents modulo $M$. For unit $t$, stationary phase gives the full normalized Kloosterman sum
$$F(t)=K_M(1,t^2)/N=\sum_{u^2=1\bmod M}e_M(2ut)=\prod_i2\cos(4\pi e_i t/M).$$
The product extends the signal to all $t$. A proper $\gcd(t,N)$ already splits the input. At the $1/N$ fraction where $t\equiv0\pmod N$, using zero in the transformed signal changes each Fourier coefficient by at most $O(\log N/N)$. [Stationary-phase identity](square-kloosterman.md)

Use the smooth transform $H(t)=\log(1+|F(t)|)$. On uniform $t\bmod M$, put $C_i=\cos(8\pi e_it/M)$ and $X_i=\frac12\log(2+2C_i)$. The $X_i$ are independent, and $\mathbb E X_i=\log 2/q_i$ by taking absolute values in the product of $1+e_{q_i}(r)$ over $r\bmod q_i$. Their variances are bounded below and their centered fourth moments bounded above by absolute constants. In fact each of $C_i\ge1/2$ and $C_i\le-1/2$ contains a constant fraction of the odd grid. Also $|2\cos(2\pi r/q_i)|\ge2/q_i$, and counting grid points near the two cosine zeros gives $\Pr(X_i<-s)\le(3/2)e^{-s}$ for $s\ge0$. These estimates are uniform in the prime powers.

Fix $i$ and write $S=\sum_{k\ne i}X_k$. If $Y=S-\mathbb ES$ and $\sigma^2=\mathbb EY^2$, independence and the fourth-moment bounds give $\mathbb EY^4\le C\sigma^4$ for an absolute $C$. Hölder gives $\mathbb E|Y|\ge\sigma/\sqrt C$. Since $\mathbb E Y=0$, Cauchy–Schwarz then gives $\Pr(Y>0)\ge1/(4C)$. Hence $\Pr(S\ge0)$ is bounded below. For fixed $S$, $H=\log(1+e^{S+X_i})$ is increasing in $C_i$. When $S\ge0$, its values on the two bands differ by at least $\log((1+\sqrt3)/2)$. The independent-copy formula for covariance now gives $\mathbb E[H C_i]\ge c$ for an absolute $c>0$. This is a real Fourier coefficient of $H$ at $4e_i$, and $\gcd(4e_i,N)=N/\sqrt{q_i}$ is proper whenever the input has at least two distinct prime divisors.

The transform is 1-Lipschitz in the trace magnitude and is $O(\log N)$ in size. A sufficiently small constant RMS error in the normalized magnitude therefore leaves an inverse-polynomial heavy coefficient after rescaling and clipping. The sampled Fourier decoder works directly on uniform $t$ and can verify factors by gcd. Even one-bit magnitude comparisons suffice. Choose a public uniform threshold on a constant-mesh grid for $H$ and test the corresponding threshold $e^h-1$ on $|F|$. The mean answer approximates rescaled $H$, and a constant raw additive gap perturbs it only over an interval of thresholds of constant length. No lower-digit isolation window is needed. With factors supplied the complete trace is a product of $\omega(N)$ cosines and takes polynomial bit time at the required precision. Computing its magnitude without factors remains open. [Fourier recovery](quarter-fourier.md)

The diagnostic verified the product identity on seven synthetic grids up to 50,625 points and checked 18 literal square-modulus Kloosterman sums. Quantizing the raw trace to a grid of step $1/64$ preserved the predicted nonunit coefficients. Dense Fourier recovery using the quantized values found a proper divisor on all six grids where it was run. Known prime powers generated the controls. [Diagnostic](stationary_log.py)

### A fixed comparison and an order-only comparison

Two weaker one-bit promises suffice. At uniform unit $t$ compare $|F(t)|$ to $1$, with either answer allowed when $|\log|F(t)||\le\varepsilon$ for a sufficiently small absolute $\varepsilon>0$. Alternatively compare $|F(2t)|$ to $|F(t)|$, allowing either answer when the ratio lies in $[e^{-\varepsilon},e^{\varepsilon}]$. Away from the respective gap the answer is $1$ on the larger side and $0$ on the smaller side. Here $F=K_{N^2}(1,t^2)/N$ on units as above. Neither comparator needs varying thresholds. [Comparison diagnostic](stationary_threshold.py)

Let $k=\omega(N)\ge2$. For the fixed comparison set $S=\sum_i X_i=\log|F|$ and $\mu=\mathbb E S$. Since $q_i\ge p_i^2$, $0\le\mu=\log2\sum_i1/q_i=O(1)$. Write $Y=S-\mu$. The moment bounds above give $\operatorname{Var}Y\asymp k$ and $\mathbb E Y^4=O(k^2)$, so $\mathbb E Y_-=\mathbb E|Y|/2\gg\sqrt k$. For any admissible bit $B$, the gap promise gives
$$\mathbb E[YB]=\mathbb E[(\mu-S)(1-B)]\ge \mathbb E Y_- -\mu-2\varepsilon.$$
This is $\Omega(\sqrt k)$ for large $k$. For bounded $k$, on the event that every local $|2\cos|\le1/2$, the bit is zero and the contribution to the right-hand expectation is uniformly positive. Each such local band has probability bounded below because the odd grid has size at least $9$. The only possible negative contribution is at most $\varepsilon$ from the gap. Choosing $\varepsilon$ small enough gives the same bound for bounded $k$.

The fixed threshold also allows the explicit gap $\varepsilon=1/4$ after trial division through 31. Then $q_i\ge37^2$ and $\mu\le(\log2)(\pi^2/8-1)<1/6$. Improve the lower bound on the covariance to $\mathbb E[(-\varepsilon-S)_+]-\varepsilon\Pr(\mu\le S\le\varepsilon)$. By Jensen the first term is at least $\mathbb E[(-a-X_i)_+]$ for any $i$, with $a=\varepsilon+\sum_{j\ne i}\mathbb E X_j\le5/12$. The continuous log-cosine lower tail at $u\ge0$ is $(2/\pi)\arcsin(e^{-u}/2)\ge e^{-u}/\pi$. The grid changes that probability by at most $2/q_i$. Integrate over $u\in[a,a+3]$ to get a lower bound $e^{-a}(1-e^{-3})/\pi-6/q_i>0.19$. The continuous log-cosine density is increasing up to $\log2$, so any interval of width $1/4$ has mass at most $(2/\pi)\arccos(e^{-1/4})<1/2$. With at most $4/q_i$ grid discrepancy, conditioning on the other components bounds the covariance loss by less than $0.127$. This uniform positive bound and the large-$k$ bound together give $\Omega(\sqrt k)$. In raw normalized units the smaller symmetric gap $[0.8,1.2]$ is contained in the proved band. For the original unnormalized sum a comparison with a gap of $\pm N/5$ around $N$ therefore suffices.

For the order comparison take $S=\sum_i(X_i(2t)-X_i(t))$ and $\mu=0$. Doubling permutes each odd grid, so the local differences are centered. Their fourth moments are uniformly bounded and their variances uniformly positive. For large $q_i$, two fixed angle sectors on which the differences are separated give the variance bound. On the remaining finite grids the differences are nonconstant, since the local cosine magnitudes at grid indices $1$ and $2$ differ for odd $q_i\ge9$. In this case $\mathbb E[SB]\ge\mathbb E S_+-\varepsilon=\Omega(\sqrt k)$ for sufficiently small absolute $\varepsilon$. These bounds allow adversarial choices at every ambiguous input.

For either comparison, changing the bit at an additional fraction $\delta$ of uniform inputs changes the total covariance by at most $\sqrt{\delta\operatorname{Var}S}=O(\sqrt{\delta k})$. Thus a sufficiently small absolute error rate away from the gap is also allowed.

In either case some centered local variable has covariance $\Omega(1/\sqrt k)$ with the bit. Its harmonics can be truncated uniformly on the odd grids. For $X_i$ use the first $D$ terms of $\Re\sum_{h\ge1}(-1)^{h-1}e_M(4h e_i t)/h$. With $Z=|2\cos(4\pi e_i t/M)|$, the Abel tail has magnitude at most $2/((D+1)Z)$. On $Z\ge1/D$, the bound $\Pr(Z\le z)\le(3/2)z$ gives $\mathbb E[Z^{-2}\mathbf 1_{Z\ge1/D}]=O(D)$. On the remaining event the fourth-moment tail estimate for $\log Z$, or direct integration of its exponential tail, gives squared error $O(\log^2(D+2)/D)$. The same error order holds for the difference under doubling, since doubling permutes the grid. Thus a polynomial choice of $D$ in $k$ leaves covariance $\Omega(1/\sqrt k)$ with the truncation. Its Fourier coefficient absolute sum is $O(\log(D+2))$. At least one nonzero multiple of $4e_i$ therefore has bit coefficient $\Omega(1/(\sqrt k\log(k+2)))$.

Apply the sampled heavy-frequency decoder to the bounded bit. At any properly nonunit query $t$ return the gcd immediately. Substituting zero at the multiples of $N$ can alter a coefficient by at most $1/N$. For each recovered $f$, test $\gcd(N,M/\gcd(M,f))$. A nonzero local frequency has conductor a power of one prime, so this test gives a proper split even when $f$ itself is a multiple of $N$. The same argument applies to repeated factors. Computing either promised comparison from public input remains open. [Sampled Fourier decoder](quarter-fourier.md)

With logarithmic gap $1/64$, an assignment minimizing the total covariance inside the gap still yielded a split on all ten dense comparison runs. The strongest factor-revealing frequencies had magnitudes between $0.1227$ and $0.1913$. Independent CRT samples also checked total and local covariances on five- and twelve-prime controls. The fixtures synthesize comparison bits from known factors. [Reproduction](stationary_threshold.py)

### Evaluating a supplied frequency and retaining the upper digit

For the unmodified $L(t)=\log|F(t)|$, any supplied nonzero Fourier frequency $f$ admits an evaluation without a supplied factorization. Compute its conductor $d=M/\gcd(M,f)$. If $d$ is a prime power, exact power and primality tests give its base $p$ and successive divisions give the full component $q=p^{v_p(M)}$. Otherwise the coefficient is zero, since each local logarithm has frequencies in a single prime-power component.

For $q$ in the support, write $r\equiv f\cdot4^{-1}\pmod q$ and $s=\min(r,q-r)$ with canonical $1\le r<q$. Aliasing the log-cosine harmonic series gives
$$\widehat L(f)=\frac{(-1)^{s-1}}2\int_0^1\frac{x^{s-1}-x^{q-s-1}}{1+x^q}\,dx.$$
This is nonzero with the displayed sign because $q$ is odd and $s<q-s$. Substitute $u=x^q$ and expand $(1+u)^{-1}=\sum_{j\ge0}2^{-j-1}(1-u)^j$. The integrals of each power are rational beta integrals. Truncating after $J$ terms gives a one-sided absolute error bounded by $2^{-J-1}(1/s-1/(q-s))$. The bit cost is polynomial in $\log M$ and the requested precision using polynomial-time primality. The executable uses exact rational arithmetic and documents its primality backends. It checked 43,984 frequencies against full-grid transforms. [Dual evaluator](stationary_dual.py)

The nonzero support contains exactly $\sum_i(q_i-1)$ frequencies across the components. For a balanced semiprime and $M=N^2$, uniform probes require order $N$ samples to reach this support. Each hit itself gives a proper gcd via its conductor. The pointwise result leaves the bulk search over frequencies open.

There is a direct loss from taking $t$ only modulo $N$. If $P_i=\sqrt{q_i}=p_i^{v_{p_i}(N)}$, averaging $X_i(t)=\log|1+e_{q_i}(4a_it)|$ over the lifts of one residue modulo $N$ gives $P_i^{-1}\log|1+e_{P_i}(4a_it)|$, where $a_i=(M/q_i)^{-1}\pmod{q_i}$. This follows by multiplying the $P_i$ factors along that local lift orbit. The conditional centered log has $L^2$ norm $O(\sum_i1/P_i)$, and the same bound holds for its difference under doubling. Thus any bounded predictor depending only on $t\bmod N$ has covariance at most this size with the centered sum used in the comparison proof. For balanced semiprimes the bound goes to zero, while an admissible bit has covariance bounded below. The lift digit must be retained by a predictor for this target.

### Complete cubic and quartic sums

There is an exact polynomial-phase alternative to the inverse sum. If $\gcd(3t,N)=1$ and $M=N^2$, then
$$K_M(1,t^2)=\sum_{y\bmod M}e_M(t(3y-y^3)).$$
Writing $y=a+Nz$ cancels all blocks except $a^2=1\pmod N$. Each stationary $a$ has a unique lift $u^2=1\pmod M$ and polynomial phase $2tu$. This is the same stationary trace. The fixed comparison at normalized magnitude $1$ therefore also applies to this complete cubic sum. Small factors of $3$ can be detected before calling the evaluator. [Polynomial-phase diagnostic](stationary_cubic.py)

The quartic $y^4-2y^2$ gives a separate target for every odd $N$. For unit $t$,
$$G(t)=\frac1N\sum_{y\bmod M}e_M(t(y^4-2y^2))=\prod_i(1+2e_M(-e_i t)).$$
Its derivative is $4y(y-1)(y+1)$. At each prime power the three stationary lifts $0,1,-1$ give the phase values $0,-1,-1$. In particular $|G(t)|\ge1$, and
$$\log|G(t)|=k\log2+\Re\sum_i\sum_{h\ge1}\frac{(-1)^{h-1}}{h\,2^h}e_M(h e_i t).$$
The coefficient at $e_i$ has magnitude at least $0.24$. The leading term is $1/4$, and all aliases are at index at least $q_i-1\ge8$. After clamping a coarse magnitude estimate below at one, its logarithm is 1-Lipschitz. The same bounded Fourier decoder and conductor gcd give a split from sufficiently small constant additive error, with no clipping window or logarithmic singularity. The exact polynomial diagnostics checked 30 cubic and 45 quartic sums on odd fixtures, including repeated factors. Literal complete sums still scale with the modulus. [Cubic and quartic checks](stationary_cubic.py)

Quadratic completion identifies where the quartic loses cancellation. Write $g_M(s)=\sum_y e_M(sy^2)$. Completing the square in the Fourier transform of $e_M(tz^2)$ gives
$$G(t)=1+\frac{e_M(-t)}M\sum_{s\bmod M}e_M(-s^2/(4t)-s)\big(g_M(s)-N\big).$$
Here the denominator $4t$ is inverted modulo $M$. Unit $s$ give $g_M(s)=N$ because $M$ is an odd square. The terms $s=Nj$ contribute exactly $2^{\omega(N)}e_M(-t)$ by orthogonality on $y^2=1\pmod N$. Every remaining nonzero term has a proper $\gcd(s,N)$. On balanced two-prime inputs these indices have density $O(N^{-1/2})$. This rules out a direct uniform sampler over the remaining terms as a polynomial-bit way to access them. An aggregate that evaluates their contribution implicitly remains open. The diagnostic verified 30 decompositions against the literal quartic sums. [Gaussian completion checks](stationary_cubic.py)

### Precompute the stationary phases from one square-congruence scan

For the quartic, a critical residue $y$ has $y^2=s$, where $s$ is an idempotent. Over a local component, $s=0$ has one critical root and $s=1$ has two. The phase value is $-s$. Enumerate roots of $u^2=1\pmod N$ by the certified square-root-width lattice blocks, lift them to $N^2$, and set $s=(1+u)/2$. The multiplicity of $-s$ is the number of lower idempotents $a$ with $sa=a\pmod N$. This equals the number of subsets of the local components on which $s=1$. The cubic phases are the same lifted roots multiplied by two.

This produces a reference from $N$ with no supplied prime factors. Its preprocessing uses $O(\sqrt N\operatorname{poly}(\log N))$ bit work for the lattice blocks and $O(4^{\omega(N)}\operatorname{poly}(\log N))$ for incidence tests. The second term is $N^{o(1)}$ since the components contain distinct primes. Each query then uses $2^{\omega(N)}$ cached weighted phases. The tests checked roots and multiplicities on eight inputs, along with 24 quartic and 15 cubic sums by literal evaluation. At $N=11009$, the reference used 787 lattice blocks and 16 incidence tests in setup, followed by four weighted phases per query. The literal symmetric calculation uses 60,599,041 phases per query. The cache reuses the same work across scales, and the preprocessing still has exponential cost in $\log N$. [Cached reference](stationary_window_reference.py)

### Limits of lower-degree pieces and linear updates

Approximating the literal phase by a lower-degree polynomial on a long unit-step progression has a direct obstruction. More generally, on a progression of unit stride modulo $M$, the degree-$d$ leading coefficient is $b/M$ with $b$ a unit. Suppose an arbitrary real polynomial of degree below $d$ matches the phase modulo one to error $\varepsilon<2^{-d-1}$ on a consecutive index interval of length $L$. Write the phase difference as $k_x+r_x$ with $k_x\in\mathbb Z$ and $|r_x|\le\varepsilon$. Each $\Delta^d k_x$ must be the unique integer within $2^d\varepsilon$ of $d!b/M$. Thus $\Delta^d r_x$ is constant and has magnitude at least $1/M$ for $M>d!$. For a stride $h$ with $dh<L$, $\Delta_h^d r_x=h^d\Delta^d r_x$, yet $|\Delta_h^d r_x|\le2^d\varepsilon$. Consequently $L=O_d(1+(\varepsilon M)^{1/d})$. A polynomial number of such pieces cannot cover the full period. The argument concerns pointwise approximation modulo one.

Even the short nonlinear forward-difference update gives a long exact linear recurrence for the literal term stream. For squarefree $N=pq$ with $p,q\nmid d$, write $M=N^2$ and $f(y)=e_M(tP(y))$ for either degree-$d$ phase and unit scale $t$. In its length-$M$ transform, the coefficient at $k_0+Nj$ is
$$\widehat f(k_0+Nj)=\frac1N\sum_{a\bmod N:\,tP'(a)=k_0\bmod N}e_M(tP(a)-k_0a)e_N(-ja).$$
Each derivative fiber contains at most $R=(d-1)^2$ points by CRT. Hence every coefficient has magnitude at most $R/N$. By Parseval, at least $M/R^2$ coefficients are nonzero. Any approximation supported on $r$ frequencies of the same period has squared RMS error at least $1-rR^2/M$. An exact constant-coefficient recurrence for the period also needs order at least $M/R^2$, since its characteristic polynomial must vanish at every supported root of unity. These bounds apply even to the small-scale queries. The diagnostic checked 40 full-grid transforms against the fiber formula, with maximum discrepancy $1.25\times10^{-16}$. The smallest observed support fraction was $0.2746$. [Frequency diagnostic](polynomial_stream_limits.py)

### Two prime components need only scales up to $2N$

Assume $N$ has exactly two distinct odd prime divisors, allowing repeated powers, and write $M=N^2$. Small $N$ and the factor 3 can be removed by trial division. For the quartic, let $e$ be one primitive CRT idempotent modulo $M$. Reflect the angle of $v=2e-1\pmod{2M}$ to choose $0<v<M$. Then $v^2=1\pmod M$ and $v$ is mixed. For $c=\cos(\pi t/M)$ and $w=\cos(\pi vt/M)$, the magnitude of the normalized quartic is
$$|G(t)|^2=(5c+4w)^2+9\sin^2(\pi t/M).$$
Indeed, the two local phases multiply to $e_M(-t)$. Rotating the expanded product by $e^{\pi i t/M}$ gives real part $5c+4w$ and imaginary part $-3\sin(\pi t/M)$. At $t\le2N$ this real part is bounded away from zero for $N\ge 32$. A magnitude oracle with absolute error at most $10^{-7}$ therefore gives $w$ to error below $10^{-5}$ using the positive square root.

The cubic magnitude has the same readout with $v=4(2e-1)\pmod{2M}$ reflected to $0<v<M$. Now $v^2=16\pmod M$ and the signed normalized trace is $2\cos(4\pi t/M)+2\cos(\pi vt/M)$. The second cosine is approximated by half the magnitude minus $\cos(4\pi t/M)$. The possible negative sign costs at most $2(1-\cos(4\pi t/M))\le 64\pi^2/N^2$ on the same scales. Thus the same accuracy bound holds for $N\ge 2^{14}$ after a constant-size trial-division precheck.

For either readout, start with $t=1$. If the first cosine is near $-1$, reflect $v$ to $u=M-v$ and multiply each readout by $(-1)^t$. Otherwise take $u=v$. If the resulting first cosine is above $0.9$, double the scale until the estimated cosine crosses $0.9$. This produces an anchor $s$ with $\theta=\pi su/M$ having sine bounded below by $0.4$. The first crossing keeps the positive angle below $1$ in the doubling case. Any mixed $u$ satisfies $\min(u,M-u)\ge N$: for a smaller representative its square minus $1$ or $16$ could vanish modulo $M$ only by integer equality with the public root. In particular $s\le N/2$ after the precheck.

At scale $s2^j$ query also $s(2^j+1)$ or $s(2^j-1)$. A proper gcd of any query scale with $N$ already suffices. The two shifted scales cannot both be multiples of odd $N$. The cosine addition law gives the sine of $2^j\theta$ to constant error using the anchor, so the circle angle is determined to error below $\rho=10^{-3}$. Dividing by $2^j$ and choosing the lift closest to the preceding estimate is unambiguous since $3\rho<\pi$. The error in $\theta$ after this step is at most $\rho/2^j$. Stop at the first $s2^j\ge N/2$. The largest queried scale is below $2N$, and the inferred $u$ has error at most $2\rho N/\pi<N/64$. The certified quadratic lattice for modulus $M$ recovers a root from this short interval, and the gcd of $u-1$ or $u-4$ with $N$ verifies the split. All public calculations and the $O(\log N)$ oracle calls are deterministic. [Quadratic lattice](lattice-window.md)

The diagnostic checked the phase intervals on 54 synthetic controls with bounded noise, including endpoint stress phases and repeated prime powers. The largest queried scale was $1.179N$, and the largest coordinate error was $8.18\times10^{-8}N$. It exercises the numerical reconstruction separately from the lattice finish. Evaluating the polynomial sums at these scales in polynomial bit time remains open. [Phase interval checks](small_scale_phase_check.py)

## One quadratic box count with no clipping

There is also a complete count target using idempotents, with no arcs or windows. Count $0\le s,y<N$ satisfying $s(s-1)=sy=0\pmod N$. In a component $q_i=p_i^{e_i}$, the selector $s=0$ allows $q_i$ choices for $y$, while $s=1$ forces $y=0$. The count is $\prod_i(q_i+1)$. For two components it yields $q_1+q_2$ by subtracting $N+1$, and a quadratic equation recovers both. Two bounded quotient coordinates let one pack the congruences into one equation with quadratic part of rank two. Additive error below $1/2$ is sufficient by rounding.

The same construction gives one count for general $N$ when the number of coordinates may grow. Add $r$ pairs $(z_j,w_j)$ with $w_j=s z_j$ and $z_j(w_j-1)=0\pmod N$. On $s=0$ the pairs vanish; on $s=1$ they give $A=2^r$ choices in place of one. The count becomes $P(A)=\prod_i(A+q_i)$. Choose $r=2\operatorname{bit\_length}(N)$, so $A>2^{\omega(N)}N$ exceeds every coefficient of $P$. Its base-$A$ digits give the exact polynomial, and integer polynomial roots give the components. For packing, bind each polynomial residue to a unique integer quotient. A base greater than twice the largest possible residual makes a single weighted sum of the residuals vanish only when they all vanish. The box, coefficients and count have polynomial bit length. This gives a deterministic one-query reduction to an ordinary count with one quadratic equation. A polynomial-bit evaluator for the count is open. [Count construction and checks](idempotent_count.py)

### Absolute error, relative error, and small residues

For two components, a power query preserves the four-coordinate count and permits much larger additive error. Write $N=ab$, $S=a+b$. At modulus $L=N^m$, subtract $L+1$ from the count to get $D_m(S)=a^m+b^m$. This has the public recurrence $D_0=2$, $D_1=S$, $D_j=SD_{j-1}-ND_{j-2}$. On real $S\ge2\sqrt N$ its derivative is
$$D_m'(S)=m\sum_{i=0}^{m-1}a^i b^{m-1-i}\ge m^2 N^{(m-1)/2}.$$
Pair terms in reverse order and apply AM–GM. The endpoint follows by continuity. Thus absolute error less than $(m^2/2)N^{(m-1)/2}$ allows monotone binary search for the unique nearest integer $S$ and exact quadratic reconstruction. Any fixed promise $L^\alpha$ with $\alpha<1/2$ is enough by taking a fixed $m>1/(1-2\alpha)$. This applies to equal-size, imbalanced, and repeated-prime components. The same excess power sum can be counted on modulus $N$ by excluding public selectors $s=0,1$ and using $m$ independent coordinates $y_i$ with $sy_i=0\pmod N$. The packed quadratic part has rank two because all quadratic monomials have a factor $s$.

When the count evaluator allows additional components, known padding permits any fixed sublinear absolute error $L^\alpha$, $\alpha<1$. For odd $N$ use $L=NP$ with $P=2^b$ public. For any previously chosen $A$, the new count is $(P+A)P_N(A)$, where $P_N(A)=\prod_i(q_i+A)$. Choose $b$ so $\alpha\operatorname{bit\_length}(N)-(1-\alpha)b<-1$. Dividing the noisy result by $P+A$ gives error below $1/2$, so rounding recovers the exact count. The padding keeps the query bit length polynomial and allows the existing deterministic radix decoder for general $N$.

Relative error behaves differently for both transforms. Padding preserves relative error. For the pair power count, write $S=2\sqrt N\cosh\gamma$ and $D_m=2N^{m/2}\cosh(m\gamma)$. Then $D_m'/D_m\le m^2/S$ by concavity of $\tanh$. Adjacent candidate sums thus have relative spacing at most $O(m^2/\sqrt N)$ for polynomial $m$. In fact ordinary relative approximation to $P_N(1)$ is easy. Trial-divide and strip every component with underlying prime at most $B=\lceil2\operatorname{bit\_length}(N)/\varepsilon\rceil$. For unresolved cofactor $R$, use $R$ in place of $P_R(1)$. The missing factor is $\prod_{q\mid\mid R}(1+1/q)\le\exp(\operatorname{bit\_length}(N)/(B+1))\le1+\varepsilon$ for $0<\varepsilon\le1$. This gives a certified interval in polynomial time in $\log N$ and $1/\varepsilon$. The diagnostic passed 42 power-query checks, 36 padded queries, and 27 relative intervals. [Precision checks](count_precision.py)

At the square-root absolute-error boundary for power queries, a fixed component ratio $b/a\ge1+\delta$ permits $m=O_\delta(\log N)$. The derivative is at least $mb^{m-1}$, so the spacing dominates $N^{m/2}$ for a sufficiently large constant in $m$. The same lower bound applies to the lower adjacent step with a smaller constant imbalance for large $N$. In the near-balanced regime the public surrogate $L+1+2\sqrt L$ can instead meet the entire error budget. If $|\gamma|\le c_0/m_{\max}$, its error is at most $2\sqrt L(\cosh c_0-1)$ at every exponent up to $m_{\max}$. Thus such outcomes can be reproduced using only $N$ in that regime. The numerical checks include two imbalanced decodes and five near-balanced surrogate controls. [Boundary checks](count_precision.py)

Small modular answers also suffice for deterministic decoding. With $K=\lfloor\log_2N\rfloor$, choose auxiliary primes $\ell$ where $2^0,\ldots,2^{K-1}$ are distinct. Excluding prime divisors of $\prod_{j=1}^{K-1}(2^j-1)$ discards at most $O(K^2)$ primes by bit length, so polynomially many small candidates give enough retained primes. At the first retained prime, interpolate $(P_N(X)-N)/X$ from $K$ count residues at $X=2^r$ and learn the degree $k$ by monicity. At each later prime, use $k-1$ points for $(P_N(X)-N-X^k)/X$. CRT with a product exceeding $2^K N$ identifies every coefficient exactly. The eight-component fixture used 84 oracle calls, with largest auxiliary prime 67. The diagnostic uses synthetic counts and independent selector scans as controls. [Residue decoder](idempotent_residues.py)

### One-bit answers and smaller boxes

Arbitrary positive $A$ is available in the same selector gadget. Use $r=\lceil\log_2 A\rceil$ binary coordinates on $s=1$. If $A=2^r$, all strings qualify. Otherwise append Boolean states for the predicate that the processed lower bits encode an integer below the corresponding prefix of $A$. Starting at zero, one step is $b'=(1-z)b$ at a zero bit of $A$ and $b'=1-z+zb$ at a one bit. These are quadratic equations. Fix the final state to one. At $s=0$ all bits are zero and automatically meet every positive threshold. At $s=1$ exactly $A$ strings qualify. Each state is unique, and bounded quotients pack the system into one quadratic equality. This gives $P_N(A)=\prod_i(q_i+A)$ with at most $6r+4$ coordinates. [Bounded selector construction](idempotent_zero_bits.py)

For a small prime $\ell$, ask only whether $P_N(A)$ is zero modulo $\ell$, scanning $A=1,\ldots,\ell$. This gives the set of distinct component residues. If $\ell$ does not divide $\prod_{i<j}(q_i-q_j)$, the set has $k=\omega(N)$ elements and determines the monic polynomial modulo $\ell$. For $n=\operatorname{bit\_length}(N)$ the difference product excludes fewer than $n^3$ primes. Scan the first $n^3+2n+1$ primes, retain those with the largest observed set size, and combine their coefficients by CRT. Their product exceeds $2^nN$, a bound on each coefficient. The maximum set size gives the degree without a supplied factor count. Every query uses one bit and $O(\log n)$ coordinates because the auxiliary primes are polynomial in $n$. The reference decoder can stop earlier on a verified prime-power decomposition. Six synthetic decodes and 30 explicit selector counts passed. The eight-component control used 568 bits with largest auxiliary prime 67. [Zero-bit diagnostic](idempotent_zero_bits.py)

For two odd components $N=PQ$, one fixed half-selector interval yields a deterministic parity-only decoder with four coordinates. The two mixed selectors sum to $N+1$, and exactly one lies in $[2,(N-1)/2]$. Its count with $0\le y\le B$ is $1+\lfloor B/d\rfloor$ for one component $d$. Query at $B=2^j$ and subtract the public one modulo two. Each parity is the next binary digit of $1/d$. At the largest $K=2^j\le N-1$, recover $h=\lfloor K/d\rfloor$. The other component $N/d$ belongs to $[Nh/K,N(h+1)/K)$, of width at most two, so direct divisor checks finish. In the full interval $[2,N-1]$ the count is instead $2+h_j$ with $h_j=\lfloor2^j/P\rfloor+\lfloor2^j/Q\rfloor$. A residue modulo three gives the next carry-free digit sum $h_j-2h_{j-1}\in\{0,1,2\}$. This determines $P+Q$ to an interval of width at most four at the same maximum bound. Both deterministic decoders use $O(\log N)$ queries. [Selector digit checks](selector_digits.py)

The full interval also supports parity alone with a sampled Fourier decoder. Its sign as a function of uniform $B\bmod N$ is $(-1)^{\lfloor B/P\rfloor+\lfloor B/Q\rfloor}=(-1)^{(B\bmod P)+(B\bmod Q)}$. Over odd modulus $d$, the local parity wave has transform magnitude $1/(d|\cos(\pi h/d)|)$, at least $2/\pi$ at $h=(d\pm1)/2$. CRT gives heavy coefficients of size at least $(2/\pi)^2$ at frequencies $(\pm P\pm Q)/2$. Recover candidates with the sampled decoder and test their doubled centered magnitudes as either a sum or a difference. The parity at the public scales $B=\lfloor N/2^j\rfloor$ also gives the exact bits of $P\mathbin{\rm XOR}Q$ by integer division. The diagnostic checked ten bounded-digit cases, including repeated prime powers, and five unrestricted parity spectra. [Spectral checks](selector_digits.py) · [Sampled decoder](quarter-fourier.md)

The complete count has much coarser native residues at fixed bases. Modulo four, strip the known power of two; the odd part contributes zero if it has at least two distinct primes and otherwise its value follows by prime-power testing. For odd $3\nmid N$, the complete count is nonzero modulo three exactly when all components satisfy $q_i\equiv1\pmod3$. Equivalently, every prime $p\equiv2\pmod3$ has even exponent, the norm criterion for $x^2-xy+y^2$. On a squarefree input the zero bit is the fixed-discriminant residuosity obstruction for $-3$. These statements concern the unshifted complete count. The varying-bound parity construction still needs its own evaluator.

### Native zero-bit evaluation at large auxiliary primes

For one supplied $A$ and auxiliary prime $\ell$, the bit $P_N(A)\equiv0\pmod\ell$ asks whether a complete prime-power component belongs to $r=-A\bmod\ell$. Strip any public $\ell$-power first. If that component does not already make the bit zero, continue on the coprime residual. There $r=0$ gives no hit. A component at most $\sqrt N$ can be found in the progression $d\equiv r\pmod\ell$. A larger component has a unitary complement at most $\sqrt N$ in $d\equiv N/r\pmod\ell$. Test $d=1$ directly. The verified candidate on the appropriate side must itself be a prime power, and $\gcd(d,N/d)=1$ ensures a complete component. [Native large-prime bit](large_mod_residue_evaluator.py)

For a progression in the divisor band $[2^j,2^{j+1})$, use $n=\operatorname{bit\_length}(N)$ and block coordinate $0\le x\le X=\lfloor 2^{j^2/n}/16\rfloor$. A divisor in this progression has the form $a+\ell x$, hence $f(x)=x+a\ell^{-1}\pmod N$ vanishes modulo that divisor. Use the $n$ rows $N^{j-i}f^i$, $0\le i<j$, and $x^i f^j$, $0\le i<n-j$, after substitution $x\mapsto Xx$. Each row is divisible by $d^j$ at a true root. The determinant is $N^{j(j+1)/2}X^{n(n-1)/2}$. Dividing the LLL evaluation bound by $(2^j)^j$ gives at most
$$\sqrt n\,2^{(n-1)/4+j/2-j^2/(2n)}16^{-(n-1)/2}<1.$$
Thus all true roots are ordinary integer roots of the reduced polynomial. Shorter blocks are checked directly. The integer exponent floor used by the executable is still within a constant factor of the stated width. Checking all roots, including the primality checks on the appropriate side, costs polynomial bit time per block.

There are $O(1+2^{j-j^2/n}/\ell)$ blocks in a band. Since $j-j^2/n\le n/4$ up to the square-root cutoff, the full cost per bit is $\widetilde O(1+N^{1/4}/\ell)$, with polynomial bit overhead and no supplied factors. This is polynomial bit time for a *specified* bit at $\ell\ge N^{1/4}$. It applies to repeated prime powers and to arbitrary numbers of components. The diagnostic verified 451 synthetic bits, including true large components found through their small complements and false composite-complement hits. [Lattice certificate and controls](large_mod_residue_evaluator.py)

On squarefree semiprimes a separate Pollard–Strassen progression evaluator costs $\widetilde O(1+N^{1/4}/\sqrt\ell)$. Scaling a progression by $\ell^{-1}$ modulo $N$ makes it contiguous for the existing product and remainder tree. Nonunit hits yield the smaller prime, whose residue verifies the answer. A progression can also be certified empty of nonunits without extracting a component. The code checked 399 bits, including 239 negative decisions without factor recovery. [Progression reference](sparse_residue_evaluator.py)

The asymptotic tradeoff is in the total query count. The earlier zero-bit decoder scans all $\ell$ possible values of $A$ at each retained auxiliary prime. With the lattice evaluator this costs $\widetilde O(\ell+N^{1/4})$ per prime, even when one bit is cheap. On two-component inputs and $\ell\nmid N$, the unstructured equality-query model has about $\ell/2$ possible unordered residue pairs $\{u,N/u\}$. Each point query tests just one pair. A random search for any nonzero-support position takes order $\ell$ queries in that model. A factorization algorithm must exploit more arithmetic than this scan, or evaluate at small primes without the $N^{1/4}/\ell$ block cost.

For two components, a dense bit at $\ell\ge N^{1/4}$ gives a separate route to bulk extraction. An inverse-polynomial advantage in predicting the centered quarter sign of $P_N(A)\bmod\ell$ suffices for sampled Fourier recovery of the component sum modulo $\ell$, followed by the native lattice. On a balanced semiprime distribution, predictors based only on $N\bmod\ell$ have exponentially small average advantage. Even an arbitrary short summary of the full input before seeing $A$ needs about $\frac14\log_2N$ bits to support inverse-polynomial advantage. [Prediction reduction and congruence bound](dense-count-prediction.md)

Adding independently sampled prime moduli gives negligible prediction advantage even when their combined CRT modulus with $\ell$ is $N^{1-\eta}$ for fixed $\eta>0$ and the target shift is fixed independently of the extra moduli. The same comparison covers a fixed vector of shifts and successive independent residues. [Several public congruences](congruence-summary.md)

### Shifted power-of-two bits and a binary search for a component

Keep the odd input $N$ fixed and let a query ask whether $2^b$ divides $P_N(A)$ at a positive shift $A$. For two complete components $N=PQ$ and an odd $A$, write $u=v_2(P+A)$, $v=v_2(Q+A)$, and $c=v_2(N-A^2)$, with $c=\infty$ when $N=A^2$. The identity
$$N-A^2=(P+A)(Q+A)-A((P+A)+(Q+A))$$
gives $c=\min(u,v)$ when $u\ne v$ and $c\ge u$ when $u=v$. For $m=\max(u,v)$, the count valuation is therefore $m+\min(m,c)$. The query is equivalent to a component in $-A\bmod 2^h$ at the public depth
$$h=\max(\lceil b/2\rceil,b-c).$$
For $c=\infty$, use $h=\lceil b/2\rceil$. An even $A$ gives an odd count. The progression lattice above uses only that the step is coprime to $N$, so it evaluates this bit from $N,A,b$ in $(1+N^{1/4}/2^h)\operatorname{poly}(\log N+\log(A+1)+b)$ bit operations. In particular, $b\ge 2\lceil n/4\rceil$ always suffices for polynomial bit time, where $n=\operatorname{bit\_length}(N)$. The test also applies at smaller $b$ when $b-c\ge\lceil n/4\rceil$.

In the other direction, to test a class $r\bmod 2^h$, take $A=-r\bmod 2^h$ and query at $b=h+\min(h,c)$. Starting from the occupied odd class, one query per level chooses an occupied child. After $\lceil n/4\rceil-1$ queries, the residue has enough bits for the native lattice to return a verified component. The box uses $O(n)$ coordinates via the arbitrary-shift construction above. The decoder never visits all classes at its terminal depth. Native execution still pays a power of $N$ on the shallow levels. For example, on squarefree two-component inputs with $N\equiv1\pmod4$, the level-two class $3\bmod4$ is occupied exactly when $-1$ fails to be a square modulo $N$.

For any number of odd components, use an upper bound $K$ on their number, such as $n$. Test the class $r\bmod 2^h$ by taking one shift from each lift of $-r$ modulo $2^{h+s}$, where $2^s>K$ and $s$ is minimal. If the class is empty all shifts have the same count valuation. If it is occupied, there are both occupied and empty leaves among the lifts. Climb from an empty leaf to the first occupied ancestor. Its two children give a pair of shifts whose valuations differ. All components outside the ancestor give equal contributions, and all components inside it contribute an extra bit of valuation only on the occupied side. The lower valuation is below $K(h+s)$, so capping at that value still detects variation. Binary search for one capped valuation and compare each remaining lift at that threshold and the next one. A class test uses $O(K+\log(K(h+s)))$ zero bits. The binary search uses $O(n^2)$ bits per split with $K=n$, then uses the same native lattice. Recurse on verified divisors for complete factorization. Every shift and answer modulus has polynomial bit length. This gives a single fixed answer prime for the oracle reduction. Native evaluation of the shallow queries for general inputs remains open. [Power-of-two evaluator and checks](prime_power_count_bits.py)

### The four-coordinate count with general odd inputs

The unamplified count $E=\prod_i(p_i^{e_i}+1)$ itself suffices for a randomized split on any odd nonsquare input with at least two distinct primes. Choose a component $p^e$ with odd $e$, so $p+1\mid E$. In the quadratic ring modulo $N$, sample public $D,a,b$ and form $x=(a+b\sqrt D)/(a-b\sqrt D)$ when its denominator has invertible norm. A proper gcd of the norm already splits $N$. On the event that $D$ is a nonsquare modulo $p$ and $(a,b)$ is nonzero there, the reduction of $x$ is uniform in the cyclic group of order $p+1$. CRT makes this component independent of any other prime component. Write $E=2^j E_0$ with $E_0$ odd. In the chain $x^{E_0},x^{2E_0},\ldots,x^E$, the first identity time modulo $p$ has each atom of probability at most $1/2$. It therefore differs from the identity time modulo a second prime with probability at least $1/2$, including the case where the second component never reaches the identity. A gcd of the two coefficients of a chain value minus one with $N$ splits at that point. Any nonunit norm at other primes still gives an early gcd when the norm is a unit at $p$. The event at $p$ has absolute constant probability. Exact square-root preprocessing reduces the all-even exponent case to this one, and primality and perfect-power tests handle leaves. This uses the same rank-two count and a randomized group calculation, with one count per split.

## Native unweighted counts and joint phase moments

The number of stationary residues can be evaluated from $N$ without listing the roots. First remove any maximal perfect power, which preserves the distinct prime factors. For a composite base $B$, strip all prime powers with prime at most $\lfloor B^{1/3}\rfloor$ using the repository's bounded Pollard–Strassen product and remainder trees. The remaining cofactor has at most two prime factors counted with multiplicity. Primality and an exact square test distinguish an empty cofactor, a prime, a prime square, and a product of two distinct primes. In the last case the cofactor can stay intact. This gives $2^{\omega(N)}$ idempotents, $2^{\omega(N)}$ cubic critical roots when $\gcd(N,6)=1$, and $3^{\omega(N)}$ quartic critical roots when $N$ is odd. [Product and remainder trees](https://doi.org/10.1137/S0097539704443793)

The bounded search costs $\widetilde O(B^{1/6})$ bit operations, including location of the small prime factors by gcd. The reference checked direct enumeration for 512 moduli and six larger controls covering prime powers and unsplit cofactors. At the 70-bit input $N=1050809297549059047257$ it left the input intact and certified four idempotents and nine quartic critical roots with 3,189 block gcds. The code uses proved primality in its classification. It computes counts of stationary points, so it leaves their separate phase weights unresolved. [Unweighted evaluator](stationary_count.py)

These counts suffice for a broad class of exact joint moments. For odd $N$, and with $3\nmid N$ when using the cubic, take uniform unit $t\bmod N^2$. Form a joint moment of normalized cubic values $F(m_jt)$ and squared quartic magnitudes $|G(m_jt)|^2$, with each multiplier a unit and fixed nonnegative integral exponents $r_j$. Each prime component $q_i=p_i^{e_i}$ contributes the same Laurent polynomial $H(z)$ evaluated at a local unit phase with denominator $q_i^2$. A cubic factor is $z^{m_j}+z^{-m_j}$, and a squared quartic factor is $5+2(z^{m_j}+z^{-m_j})$. The absolute degree of the joint polynomial is at most $D=\sum_jm_jr_j$. If $D<\min_i q_i^2/p_i$, the only retained local term is the constant coefficient. The average over units of $e_{q_i^2}(h\,t)$ is zero at nonzero $|h|<q_i^2/p_i$. CRT independence then gives the exact joint moment $([z^0]H)^{\omega(N)}$. In particular these moments carry the unweighted component count and no phase location data. The diagnostic matched 25 joint moments and two cases at the first alias boundary. [Moment controls](phase_moment_control.py)

Run the native count with the dependencies used by the existing polynomial product code.

```sh
uv run --no-project --with sympy --with python-flint==0.9.0 --with gmpy2 python experiments/polylog_factoring_2026_11_27/stationary_count.py
```

## Factor-supplied range control and hardness of arbitrary windows

The retained reciprocal spectrum can be evaluated with $O((K+1)2^{\lceil\omega(N)/2\rceil}\operatorname{poly}(\log N))$ work per digit at fixed accuracy and a fixed number of denominator intervals, given the prime powers of $N$. For a retained slope $b$, square roots of $(1-b)^{-1}$ modulo the prime-power components of $N^2$ give CRT sign choices. Lift a candidate lower digit to $r=a+Nh$ with $(1-b)r^2=c^2\pmod{N^2}$. Then the phase at $a$ is $(2-2b)r+ba\pmod{N^2}$ by subtracting the lift change $Nbh$ from the phase at $r$. Its exponential factors across two half lists, with a single wrap correction $e_N(-b)$ when the lower digits sum to at least $N$. Sort the right half list by lower digit and take prefix sums of phase weights. Each left choice queries shifted right intervals where the geometric lift weight is constant. The ordering is reused between digits. Square roots and lifting are included in the setup. The input factorization is supplied. [CRT evaluator](crt_range.py)

The implementation passed 352 local-root checks and 420 comparisons against the retained lower-digit scan, including repeated factors. On $N=35336848261$ with eight distinct primes, a ten-million-digit window and $K=32$ gave five retained slope classes. The reference used 81 phase evaluations and 160 range searches per tested digit after setup. The fixture window was placed around a known mixed root. Run `python crt_range.py 11 13 17 19 23 29 31 37 --start 9632322254 --width 10000000 --cutoff 32 --digit 0 --digit 17` to reproduce it. [Reproduction and counters](crt_range.py)

Exact clipped root counting is #P-hard on arbitrary windows, even given the factorization and a unit whose square is the target. For positive subset weights $w_i$ and target $T_0$, write $T=2T_0-\sum_i w_i$ and choose $L=16(1+\sum_i w_i+|T|)$. Choose powers $q_i$ of distinct small odd primes, all larger than $2^{m+8}mL$ for $m$ weights, and take $N=\prod_i q_i$. Choose a unit $k_i\bmod q_i$ with $|k_i/q_i-w_i/L|\le2/q_i$ by rounding and adjusting one unit. The CRT sign choices for the known square with components $(N/q_i)k_i$ all lie at circle positions $\sum_i\sigma_i k_i/q_i\bmod1$. A cyclic window of radius $1/(4L)$ centered at $T/L$ selects exactly the target subset sums. Every parameter has polynomial bit length.

The same construction gives hardness of constant-additive-error evaluation of the real clipped stationary sum. Choose $k'_i\in[1,q_i)$ with $(N/q_i)k'_i\equiv k_i\pmod{q_i}$ and set $V_i=(N/q_i)^2k'_i$. The roots of $(\sum_i V_i)^2$ modulo $N^2$ are the signed sums of $V_i$. Their lower digits are unchanged and each phase $e_{N^2}(2r)$ is within $4\pi\sum_i1/q_i$ of one. The sum of all phase errors is below $1/8$. Approximating the clipped real sum to $1/4$ and rounding recovers the exact subset count. Two shifted half-width queries with error at most $1/8$ also suffice. Their difference removes a copy of the target window shifted by $(N-1)/2$, where there are no roots. The reduction uses adversarial windows with an unbounded number of distinct primes. It leaves the random isolation promise and the unmasked magnitude target open. The diagnostic checked 49 encodings, with moduli up to 375 bits. [Encoding diagnostic](stationary_hardness.py)

## Acceptance test for an evaluator

Use the unfactored input and public window parameters. Compare the exact residue with a direct row scan on both sides of a transition. For the first check, use the two semiprime controls 11,009 and 38,942,974,393, then add repeated-factor inputs such as 315 and 1,575. Include queries with no mixed roots, queries containing multiple mixed roots, and shifts for which the public-root or axis correction changes the raw count. Charge preprocessing, empty blocks, endpoint corrections and any auxiliary modulus construction. After the finite check, give a bound in the bit length for the entire split, including every query the decoder needs.

Full-period totals, local CRT counts taken independently, and representations supplied after finding the hidden root are correctness controls. A successful evaluator must account for which two local residues occur in the same clipped integer interval. On the class formulation, require counts only at the powers actually queried by the reduction. A full-orbit evaluator solves a stronger task.

## Restricted approaches already checked

| Approach | Established restriction |
|---|---|
| Explicit separation of the cutoff into functions of the two prime residues | The unfiltered row can have rank $\Omega(\sqrt N)$. The mixed filtered table can still have rank $\Theta(N^{1/8})$, even when one direction suffices for a closed-form sum. [Full boundary](fixed-cutoff-crt-rank.md) · [Filtered table](filtered-crt-rank.md) |
| Stationary linear state updates | The shift-indexed sequence requires $\Omega(\sqrt N)$ dimension on balanced semiprimes. The full odd-class orbit needs order greater than $L/4$. Read-once binary matrix products can also have exponential width on general inputs. [Shift](window-threshold.md) · [Class power](class-exponent-period.md) · [Binary width](binary-prefix-rank.md) |
| Uniform-input modular samplers | Low-degree rational maps have small success under the square decoder. On independent CRT outputs a pairwise gcd captures rare inverse-polynomial success. Equality branches can be instrumented with gcd probes. These arguments leave integer-representative operations open. [Rational maps](uniform-rational-sampling.md) · [Power maps](power-sampling.md) · [Independent outputs](independent-sampling.md) · [Equality branches](arithmetic-branching.md) |
| Universal shift certificates for mixed roots | Removing the public roots leaves the low-degree ideal of vanishing certificates unchanged. The argument is for every integer lift and leaves interval-specific certificates open. [Proof](mixed-root-certificates.md) |
| Fixed finite families of unweighted cyclic binomial folds | Infinitely many balanced semiprimes make all individual folds units although the exact middle coefficient splits. The proof uses bounded gaps in a fixed progression and does not address growing width families or combined outputs. [Proof](fixed-cyclic-folds.md) |

## Complex phase of a cyclotomic polynomial

For distinct primes $p,q$ and $N=pq$, take $z=\exp(i\pi/N)$. The identity
$\Phi_N(z)=(1-z^N)(1-z)/((1-z^p)(1-z^q))$ gives principal argument
$\pi(N+1-p-q)/(2N)$: each factor $1-e^{it}$ has argument $t/2-\pi/2$ for $0<t\le\pi$.
Thus rounding $2N\arg(\Phi_N(z))/\pi$ recovers the totient, then the factor sum.
The modulus is at least $4/\pi^2$ by $2x/\pi\le\sin x\le x$ on $[0,\pi/2]$.
Absolute error $1/(10N)$ suffices, so precision costs only $O(\log N)$ bits.
This is an evaluator equivalence for semiprimes. With factors the four powers are cheap.
From $N$ alone, the formula for $z^N-1$ includes roots of every order dividing
$N$; extracting the primitive-root contribution in polynomial bit time is the open step here.

### Constant absolute precision via conjugate samples

For each auxiliary prime $r\ge17$ coprime to $N=pq$ (distinct primes), put
$z_k=\exp(2\pi i k/r)$. Use indices for which both $k$ and $Nk$ modulo $r$
lie in $[r/8,7r/8]$. The four-chord formula bounds
$|\Phi_N(z_k)|\ge\sin^2(\pi/8)>0.14$ there. Reciprocity gives
$\Phi_N(z_k)/\overline{\Phi_N(z_k)}=z_k^{\varphi(N)}$.
Given complex values within 0.01 in absolute distance, their ratios to their
conjugates err by at most $0.02/(0.14-0.01)<0.16$.
Enumerate the residues $a$ modulo $r$ and retain those whose predicted ratios
$z_k^a$ lie within 0.2 of every observed ratio. For any wrong $a$, multiplication
by $a-\varphi(N)$ permutes residues. There is a $k$ for which all three of
$k,Nk,(a-\varphi(N))k$ are in the indicated interval modulo $r$: each
condition excludes at most $r/4+1$ indices and $3(r/4+1)<r$.
At this index the two predicted ratios differ by at least $2\sin(\pi/8)>0.76$.
Hence the surviving residue is unique. CRT over primes with product exceeding
$N$ recovers the totient. The prime number theorem allows $r=O(\log N)$
asymptotically, and an auxiliary prime dividing $N$ already splits it.
This reduces the numerical specification of a sufficient evaluator to constant
absolute accuracy on small-order roots. It still requires a way to evaluate the
primitive contribution from $N$ alone. In particular $N$ modulo one probe
order cannot generally give the local answer: $33$ and $305$ are both $16$
modulo $17$, but their totients are $20$ and $240$, distinct modulo $17$.

### Moving the probe inside the disk

Write $N=pq$ with odd primes $p<q$, and put $a=R^p$ for $0<R<1$.
Throughout $|z|\le R$ the four-chord identity gives
$|\Phi_N(z)-(1-z)|\le 2(2a+2a^2)/(1-a)^2$ by expanding the numerator
of the difference. At $R=\exp(-7/p)$ this is below 0.004. Thus the same
linear polynomial already meets the 0.01 tolerance throughout that disk.
At $R=\exp(-1/p)$, however, $z=-R$ gives
$(1+R)-\Phi_N(-R)=(1+R)(a+R^q+aR^q-R^N)/((1+a)(1+R^q))$,
at least $(e^{-1}-e^{-q})/4>0.08$.
The transition for a fixed absolute signal therefore occurs at radial distance
of order $1/p$ from the unit circle. At the origin the first coefficient where
$\Phi_N$ differs from $1-z$ is the coefficient of $z^p$, equal to 1.
This excludes extracting the first deviation by reading a polynomial-in-bit-length
Taylor prefix for balanced semiprimes. The argument gives no lower bound on
evaluation by a different representation.

### Magnitude-only samples with 1% relative error

For $N=pq$ with distinct primes, let $r>5$ be a prime coprime to $N$.
Use every nontrivial $r$th root $z$ and a 1%-relative estimate of
$|\Phi_N(z)|$. Divide $|(1-z^N)(1-z)|^2$ by the squared estimate.
The exact quotient is $D(z)$ for the Laurent polynomial
$D(X)=(2-X^p-X^{-p})(2-X^q-X^{-q})$. Since $0\le D(z)\le16$,
the absolute error is at most $16(0.99^{-2}-1)<0.33$. Supply $D(1)=0$.
The normalized inverse discrete Fourier transform now gives every integer
coefficient $c_j$ of $D$ modulo $X^r-1$ within 0.33, so rounding is exact.

The sixth moment obeys
$\sum_{j=0}^{r-1} j^6 c_j=60N^2(p^2+q^2)\pmod r$.
Indeed exponents may be reduced modulo $r$ before taking their sixth powers.
Before reduction, the moment is the sixth derivative of $D(e^t)$ at zero:
each factor begins $-p^2t^2-p^4t^4/12$ or its counterpart for $q$.
Divide modulo $r$ by $60N^2$. CRT with product exceeding $N^2$ gives
$p^2+q^2$ as an integer and adding $2N$ gives $(p+q)^2$.
The prime sizes and transform lengths are polynomial in the input bit length.
At these probes $|\Phi_N(z)|\ge4/r^2$ by the chord bound, so ordinary
polynomial-bit arithmetic suffices for the decoder. The outstanding evaluation
problem is to obtain those real magnitudes from $N$ alone.

### Uniform sampling of the nonprimitive-root correction

The product over all nonidentity $N$th roots has public magnitude
$B=|(z^N-1)/(z-1)|$. To obtain $\log|\Phi_N(z)|$, subtract from $\log B$
the sum of $\log|z-e^{2\pi i j/N}|$ over $1\le j<N$ with $\gcd(j,N)>1$.
For distinct primes $N=pq$, a uniform index hits this correction with exact
probability $h=(p+q-2)/(N-1)$. Independent uniform sampling of the weighted
indicator gives zero correction on the event of no hit, of probability $(1-h)^m$
after $m$ draws. The first hit itself factors $N$.

This causes actual magnitude error on an infinite family. At $z=e^{2\pi i/5}$
take distinct primes $p,q\equiv2\pmod5$. Then $B=1$, whereas the
four-chord formula gives $|\Phi_N(z)|=(\sin(\pi/5)/\sin(2\pi/5))^2
=(3-\sqrt5)/2$. For comparable such primes (available in dyadic intervals
by the prime number theorem in arithmetic progressions), achieving 1% relative
accuracy with probability at least $2/3$ requires at least
$\log(3)/[-\log(1-h)]=\Theta(\sqrt N)$ draws for this estimator.
This only addresses the estimator that samples the correction term by term.

### When is a prime-order magnitude exactly one?

For distinct primes $p,q$, let $N=pq$ and let $r\ge7$ be prime to $N$ and prime.
At any one primitive $r$th root $z$ the following are equivalent:
$|\Phi_N(z)|=1$, and at least one of $p,q$ is $\pm1\bmod r$.
Sufficiency follows by cancellation of chord lengths. For necessity, the equality
$\Phi_N(z)\Phi_N(z^{-1})=1$ propagates to all conjugate roots. The denominator
spectrum $D$ from the sixth-moment calculation above then equals the numerator
spectrum obtained with exponents $1,N$. Their moments give
$p^2+q^2\equiv1+N^2\pmod r$, hence $(p^2-1)(q^2-1)\equiv0\pmod r$.

This is an exact local factor-property test if the evaluator is available.
Even 1% relative error at every conjugate suffices to resolve it by the
rounded Fourier decoder. For each fixed public residue $n\ne0\bmod r$,
both outcomes occur among semiprimes: choose factor residues $1,n$ for equality,
or choose $a$ outside $\{1,-1,n,-n\}$ and use $a,n/a$ for inequality.
Dirichlet's theorem supplies distinct primes in either pair of classes.
At the smaller probe $r=5$, $N\equiv\pm2$ forces one factor to be $\pm1$,
so all four magnitudes are then exactly one from public data.

### How much do all the unit-magnitude flags reveal?

Choose two distinct primes uniformly without replacement from $[X,2X]$.
Let $B$ record the preceding yes/no test for every prime $7\le r\le y$,
where $y=(\log X)^A$ for fixed $A>0$. Then its Shannon entropy in bits is
$H(B)=O(\log y)=O_A(\log\log X)$. Indeed Brun–Titchmarsh bounds the
number of primes in either class $\pm1\bmod r$ in the interval by
$O(X/(r\log X))$ uniformly in this range. The prime number theorem and a
union bound give $\Pr(B_r=1)=O(1/r)$. Binary entropy is therefore
$O((\log r)/r)$, including small primes by enlarging the constant.
Entropy subadditivity and $\sum_{r\le y,\ r\text{ prime}}(\log r)/r
=O(\log y)$ finish the proof without any independence assumption.

Algorithmic consequence on this distribution: suppose an algorithm factors
each input in polynomial time given its true flag vector. Sample another
independent pair of primes in the same interval, compute its flags, run the
algorithm with these guessed flags using the polynomial time cutoff, and verify
by division. Uniform prime samples take expected polynomial time by rejection
sampling with a primality test. For at least a $1-\epsilon$ fraction of inputs,
the probability mass of their true vector is at least $2^{-H(B)/\epsilon}$,
by Markov applied to the surprisal. For fixed $\epsilon>0$, polynomially
many independent guesses therefore suffice for constant success on each
of those inputs. This removes the flags at polynomial overhead for a typical
balanced input in this model. It gives no worst-case simulation, and it does
not cover the full magnitudes or predicates with variable target residues.

### Choosing a short list of probe points from public data

Fix a prime $r\ge7$ coprime to $N$. For each nonzero residue $a$ set
$b=N/a\bmod r$. The candidate denominator spectrum at $z_k=e^{2\pi i k/r}$
is $D_a(k)=|1-z_k^a|^2|1-z_k^b|^2$. Group candidates by the unordered
pair $\{a^2,b^2\}$ (with multiplicity); members have identical spectra.
Distinct groups have distinct sums of squares since their products are fixed.
The sixth-moment identity above thus makes their integer coefficient vectors
different modulo $X^r-1$. Parseval gives
$r^{-1}\sum_k|D_a(k)-D_c(k)|^2\ge1$ for distinct groups.
Since each spectrum lies in $[0,16]$, at least $(1-0.9^2)/256$ of all
indices have difference exceeding 0.9.

A deterministic greedy choice now separates every pair of groups using
$O(\log r)$ points. From the remaining pairs pick an index distinguishing
the most pairs, and remove those it distinguishes. Compute the public spectra
within 0.01 and count a difference above 0.85 as distinguishing.
Every difference above 0.9 is counted, and every counted difference exceeds
0.83 in truth. Averaging removes a fixed fraction of the outstanding pairs
each time. Index zero never distinguishes a pair. The candidate table and
the greedy search have polynomial cost in $r$. One-percent relative magnitude
estimates at just the selected points give denominator estimates within 0.33,
so matching the approximate table to tolerance 0.35 recovers the unique group
(allow another 0.005 for decoder arithmetic). This determines whether *either*
factor is $\pm a\bmod r$ for every target $a$ simultaneously. The public
selection uses only $N\bmod r$; evaluation at the selected points is still open.

### One log-magnitude component decides restricted quadratic residuosity

Take odd distinct primes $p,q$, $N=pq$, and a probe prime $r\equiv1\pmod4$
coprime to $N$. Put $\chi(k)=(k/r)$, $z_k=e^{2\pi i k/r}$, and
$C_r=\sum_{k=1}^{r-1}\chi(k)\log|1-z_k|$. The weighted log-magnitude sum is
$$\sum_{k=1}^{r-1}\chi(k)\log|\Phi_N(z_k)|
=C_r(1-\chi(p))(1-\chi(q)).$$
This follows from the four chords by permuting $k$ by each exponent.
The public constant is nonzero: the Abel-limit Fourier series for
$\log|1-z_k|$, followed by the real quadratic Gauss sum, gives
$C_r=-\sqrt r\,L(1,\chi)$, and Dirichlet nonvanishing applies.
When $(r/N)=+1$, reciprocity says the two hidden signs agree.
The sum is zero precisely when $r$ is a square modulo $N$, and $4C_r$
otherwise. For Jacobi symbol $-1$ the sum is publicly zero.

At $r=5$ a single magnitude already suffices on the Jacobi-positive promise.
The two positive chord lengths are $u=2\sin(\pi/5)$ and $v=2\sin(2\pi/5)$.
Both local symbols positive gives magnitude one; both negative gives
$u^2/v^2=(3-\sqrt5)/2$, since $N\equiv\pm1\pmod5$. Thus even a fixed probe
with 0.01 absolute error decides whether 5 is a square modulo a Jacobi-positive
odd semiprime. This concerns specified small prime radicands. It supplies
neither a reduction from arbitrary radicands nor a factoring algorithm.

### A bounded real norm equation for the fifth-root reading

Under the preceding Jacobi-positive promise with $r=5$, the magnitude is one
exactly when $4N+5y^2=t^2$ for some integer $1\le y\le\sqrt{N/5}$.
In that case there are exactly two such positive $y$, and finding both splits $N$.
Indeed the ring of integers of $\mathbb Q(\sqrt5)$ has class number one
(Minkowski bound $\sqrt5/2<2$). A split prime gives two prime ideals of
norm that prime; an inert prime cannot occur to odd exponent in an element
norm. In the split case the four ideals of norm $N$ each have a totally positive
generator $\alpha$ of norm $N$, adjusting signs with $-1$ and the norm-minus-one
unit $\phi=(1+\sqrt5)/2$. Scaling by powers of the fundamental totally positive
unit $\phi^2$ puts $\alpha/\sqrt N$ between $1/\phi$ and $\phi$ uniquely.
Equality at the endpoints would make $\sqrt N$ an element of the field, impossible
here. Write $\alpha=(t+y\sqrt5)/2$. Then
$|y|=|\alpha-N/\alpha|/\sqrt5<\sqrt{N/5}$ and conjugate ideals pair
opposite signs. Conversely a positive square test gives a generator since
$t,y$ have the same parity; the bound places it in this same interval.

Each witness gives a root $t/y$ of 5 modulo $N$: a prime dividing both $y$
and $N$ would also divide $t$, contradicting squarefreeness. Choosing the
conjugate ideal above a prime flips that local root's sign. Thus the two
positive witnesses yield roots differing in sign at exactly one prime, and
the gcd of their difference with $N$ splits it. For $N=341$ the roots are
37 and 180, with difference gcd 11. Scanning the reduced interval uses
$O(\sqrt N)$ square tests. Unit reduction bounds the search but supplies no
polynomial-bit-time method for finding its two witnesses.

### Cost of a full-period modular sieve for the norm witnesses

Put $Y=\lfloor\sqrt{N/5}\rfloor$. For any prime $r\nmid10N$, let $A_r$
count the residues $y$ for which $4N+5y^2$ is a square modulo $r$, including
zero. With $\chi$ the Legendre symbol modulo $r$,
$$2A_r=r-\chi(5)+1+\chi(-20N).$$
Indeed $\sum_y\chi(4N+5y^2)=-\chi(5)$: for nonzero $c$, counting
$(u-y)(u+y)=c$ gives $\sum_y\chi(y^2+c)=-1$. The polynomial has
$1+\chi(-20N)$ zeros, which require the extra unit weight in $2A_r$.
In particular $A_r/r\ge1/3$.

For any squarefree $M\le Y$ coprime to $10N$, the interval contains
$\lfloor Y/M\rfloor$ full periods. CRT gives at least
$(Y/2)3^{-\omega(M)}$ survivors. Since the product of $k$ distinct primes
is at least $k!$, uniformly $\omega(M)=O(\log Y/\log\log Y)$.
Thus even choosing the primes using $N$ leaves $Y^{1-o(1)}$ candidates for
explicit enumeration as long as the combined period fits the search interval.
For larger combined periods the full-period argument gives no bound.
Checked the exact local count for every nonzero input residue at every eligible
odd prime below 100.

### Simultaneous isolation of the norm witnesses

Assuming GRH for Dirichlet L-functions, all primes $r\le C(\log N)^2$
coprime to $10N$ together remove every false candidate in the reduced interval,
for an absolute constant $C$. For $D=4N+5y^2\le5N$ nonsquare, the
quadratic character of $\mathbb Q(\sqrt D)$ induces a nonprincipal
character on the units modulo $m=40ND$ (its conductor divides $4D$).
Bach's GRH bound for generators of the unit group gives a prime
$r\nmid m$ of size $O(\log^2 m)$ outside its kernel. Thus $(D/r)=-1$.
The same cutoff works simultaneously for all $y$, and a square always passes.
The allowed residue tables can be built from $N$ in polynomial bit time.

There is also an unconditional randomized isolation statement with larger
moduli. Draw independent odd positive $u\le100N$ uniformly and reject a candidate
when the Jacobi symbol $(D/u)=-1$. For each nonsquare $D$, in every full
period modulo $4D$ exactly half the units give $-1$ by quadratic reciprocity
and nonprincipality. Full periods cover at least half the sampling range,
so a draw rejects this $D$ with probability $\Omega(1/\log\log N)$
by the standard lower bound on $\varphi(4D)/(4D)$. Squares never give $-1$.
A union bound over at most $\sqrt N$ candidates shows that
$O(\log N\log\log N)$ draws with a sufficiently large constant isolate
all witnesses with probability at least $0.99$. These Jacobi predicates
are compact to evaluate on a given candidate, but their full residue tables
may be large. Neither isolation statement gives an algorithm to list the
surviving integers in polynomial bit time.

### A fixed menu of real norm equations covers every distinct-prime semiprime

Strip primes at most 5. For $N=pq$ the Legendre symbols at $p,q$ give two
linear functionals on the three-dimensional squareclass space generated by
$2,3,5$ over $\mathbb F_2$ (encode $+1$ by zero). Their joint kernel has
dimension at least one. Thus some $d\in\{2,3,5,6,10,15,30\}$ splits at both
primes, meaning $d$ is a nonzero square modulo both. One may try the entire
fixed menu; knowing the factors is unnecessary to construct it.

For that $d$, in the ring of integers of $\mathbb Q(\sqrt d)$ choose one
prime ideal over each of $p,q$. Their product $I$ has norm $N$.
Minkowski's bound gives an integral ideal $J$ in its inverse class of norm
$c\le\sqrt{\Delta_d}/2\le\sqrt{120}/2<6$. Thus $IJ=(\alpha)$,
with $|\operatorname{Norm}(\alpha)|=cN$ and $1\le c\le5$. Both primes
are coprime to $c$. Multiplication by a fixed field unit balances the two
real embeddings to size $O_d(\sqrt N)$: multiply by a suitable power of a
unit whose first embedding exceeds one until that embedding has magnitude
between $\sqrt{cN}$ and a constant times $\sqrt{cN}$.
Writing $2\alpha=t+y\sqrt d$ therefore gives
$t^2-dy^2=\pm4cN$ with $|t|,|y|=O(\sqrt N)$ uniformly in the fixed menu.
Since $p,q$ divide the right side only once, $y$ is invertible modulo $N$.
At the chosen prime ideals $-t/y$ is the reduction of $\sqrt d$.
Changing the ideal above just $p$ flips the root sign there and preserves it
at $q$, even if the multiplier $c$ changes. A gcd of the difference of the
two recovered roots splits $N$. Therefore enumerating bounded witnesses for
these finitely many conics suffices for every such semiprime. The coverage and
the small multipliers are unconditional; finding the witnesses efficiently
is still open.

### What two-dimensional lattice reduction needs

For $N=pq$ with distinct primes greater than 5, fix a menu entry $d$ and
a root $u^2=d\bmod N$. The lattice $L_u=\langle(N,0),(u,1)\rangle$
has index $N$. Gauss reduction in the positive metric $t^2+dy^2$ returns
a shortest nonzero vector in polynomial bit time using integer scalar products.
The scaled covolume is $\sqrt d\,N$, so the planar Hermite bound gives
$0<|t^2-dy^2|\le t^2+dy^2\le(2/\sqrt3)\sqrt d\,N<7N$.
The strict lower bound uses that $d$ is nonsquare. Hence the norm is $kN$
with $1\le|k|\le6$. A prime dividing $y$ and $N$ would also divide $t$
and then divide $k$, impossible here. The short vector therefore returns
exactly the original root $t/y=u\bmod N$.

There is a converse obstruction for any proposed additive lattice on which
the norm vanishes identically modulo $N$. If it contains $(s,v)$ with
$v$ invertible modulo $N$, polarization against any $(t,y)$ gives
$st=dvy\bmod N$. Since $s/v$ is a unit square root of $d$, this forces
$t=(s/v)y\bmod N$ throughout the lattice. In particular one such lattice
cannot combine two incompatible root choices. If a proposed basis has
second coordinates $y_1,y_2$, compute $g=\gcd(y_1,y_2,N)$. A proper
nonunit gcd factors $N$; $g=1$ gives a vector with second coordinate 1
modulo $N$ by Bezout; $g=N$ forces all coordinates to vanish modulo $N$
by squarefreeness. Thus constructing a nontrivial norm-divisible lattice
already supplies a root or factor before reduction. This applies only to
additive lattices with this exact divisibility property.

### Rational secants from one norm witness

Start at an integral point $(t,y)$ of $t^2-dy^2=kN$ with $\gcd(2dky,N)=1$,
and put $u=t/y\bmod N$. Intersect the conic again with the line through the
point of integer slope $m$ in the $(y,t)$ plane. Write the second point as
$(T/Z,Y/Z)$ with $Z=m^2-d$, $Y=y(m^2+d)-2mt$, $T=mY+(t-my)Z$.
Reduce the homogeneous triple by its common gcd. For each $s=p,q$, the new
ratio $T/Y$ after reduction is $+u\bmod s$ if $m=u\bmod s$ and $-u$
otherwise. Indeed $yY=(my-t)^2-kN$. Off the exceptional residue,
$Y=y(m-u)^2\ne0\bmod s$ and the formulas give $T=-uY\bmod s$.
On the exceptional residue $Y$ has valuation exactly one, $s\mid Z$,
and $s^2\mid T-mY$. Cancellation therefore gives $T/Y=m\bmod s$.

Consequently a secant gives an incompatible root exactly when
$1<\gcd(m-u,N)<N$. For uniform integer slopes modulo $N$ this has probability
$(p+q-2)/N$. The gcd already splits $N$ before the intersection.
Rational slopes $a/b$ with $b$ invertible modulo $N$ have the same criterion
after replacing $m$ by $a/b$ locally (clear denominators by the unit $b^2$).
A reduced denominator sharing just one prime already splits $N$; if divisible
by $N$, the slope is vertical at both primes and gives the opposite signs
at both by the same homogeneous formulas. This accounts for one-step rational
parametrization from a known integral witness.

### Three input-selected fields suffice; seven fixed fields are sharp

Compute the Jacobi signs at $N$ on the basis $2,3,5$. Their kernel in the
three-dimensional squareclass space has dimension at least two. Choose a
two-dimensional plane in the kernel. Its three nonzero classes are the whole
input-specific menu: both hidden primes have identical Legendre signs on
this plane, and a linear functional on a plane has a nonzero kernel vector.
Explicitly, keep all three Jacobi-positive menu entries when there are three;
when all seven are positive, use $2,3,6$. This preserves the previous
small norm multipliers. On the all-positive Jacobi pattern even two entries
can both have negative local signs, so Jacobi data on this basis alone
cannot guarantee a shorter list of its entries.

For comparison, no fixed list of six nonsquare rational squareclasses covers
all semiprimes. In their binary span a uniformly chosen nonzero linear
functional vanishes at each listed class with probability strictly below
one half. Some such functional thus vanishes on at most two listed classes.
A second functional can take value one on both remaining classes, since
two distinct nonzero binary vectors are independent. Their common kernel
misses the list. Any prescribed sign functional occurs at infinitely many
unramified primes by Chebotarev in the multiquadratic extension attached to
a basis of the span. Choose two distinct primes realizing the two functionals.
The same realization justifies the sharpness assertion on the all-positive
Jacobi pattern by taking the two local functionals equal. These bounds
concern field selection, not the cost of producing norm witnesses.

### Extending field selection beyond two prime factors

Let $k\ge2$ bound the number of distinct primes in an odd nonsquare $N$.
Use the first $k+1$ primes as independent squareclasses, first taking gcds
with $N$. In the public Jacobi kernel choose any $k$-dimensional subspace.
On it the hidden Legendre functionals have rank at most $k-1$: the sum of
the rows for odd prime exponents vanishes there and the sum involves at
least one row. Hence some nonzero class among these $2^k-1$ entries splits
at every prime. Its square roots lift to all prime powers by Hensel's lemma.
This covers repeated factors as well. Perfect powers can be detected separately.

Write $P$ for the product of the first $k+1$ primes. Every radicand $d$ is
at most $P$. For each CRT root $u$ the preceding planar lattice gives a
witness $t^2-dy^2=jN$ with $0<|j|\le2\sqrt{P/3}$ and
$t^2+dy^2\le2\sqrt{P/3}\,N$. If $\gcd(y,N)>1$ it already splits $N$,
unless $N$ divides both coordinates; that is impossible for $N>2\sqrt{P/3}$.
Otherwise $t/y$ returns $u$. Two CRT roots differing at just one prime
give a proper factor when $N$ has at least two distinct primes.
For $k\log(k+2)=O(\log\log N)$ the primorial bound gives
$P=(\log N)^{O(1)}$. The complete list of equations is then polynomial
in the input bit length and may be generated given the bound $k$, without
knowing the individual primes or exponents. Finding witnesses is still open.

The original seven fields can all fail on a product of three primes:
for $N=7\cdot11\cdot13$, none of their radicands is a square at all three.
This is a coverage limitation distinct from the cost of solving a selected conic.

### Sharp ambiguity even with the full residue modulo the basis primes

For the first $k+1$ prime generators put $M=8\prod_{i=2}^{k+1}p_i$ and
$V=\mathbb F_2^{k+1}$. On units modulo $M$ the vector of Legendre signs on
these generators defines a surjective homomorphism to $V^*$: choose the
class modulo 8 equal to 1 or 5 and prescribe the quadratic character at each
odd generator by CRT and reciprocity. Fix one unit class $C$ giving a nonzero
sign functional $J$. Every nonzero $v\in\ker J$ occurs as the **unique**
common split class for some squarefree $N\equiv C\pmod M$ with $k$ factors.

Indeed $v^\perp$ is a $k$-dimensional space containing $J$. It has a basis
whose sum is $J$, by applying an invertible linear map to any basis and
its nonzero sum. Choose unit residue classes realizing the first $k-1$
basis rows and take the final class to make their product $C$. Its signs
then give the final row. Linnik's least-prime bound supplies a prime of
size $M^{O(1)}$ in each class. Their independent rows make the primes
distinct. Their product has precisely the claimed one-dimensional common
kernel and has bit length $O(k^2\log(k+2))$.

Hence a field selector restricted to $k$ and this full public residue must
retain all $2^k-1$ possibilities in the chosen squareclass space. Even if
the selector also knows the bit length, some length needs at least
$(2^k-1)/O(k^2\log(k+2))$ entries by pigeonholing these inputs.
This lower bound does not constrain use of other information in $N$ or
radicands outside the specified span.

### One input-dependent conic for every odd non-prime-power input

After square and gcd checks, set $a=\lfloor\sqrt N\rfloor$,
$D=N-a^2>0$, with $\gcd(a,N)=1$. The public conic
$T^2+DY^2=NZ^2$ has the point $(a,1,1)$. Here $-D$ is a known unit
square modulo $N$ regardless of the number of distinct prime divisors.
For any full prime power $h=p^e\parallel N$, $h<N$, there is a primitive
point of polynomial bit length with $\gcd(T-aY,N)=h$.
For existence take $m=a+h$, $Z=m^2+D$, $Y=h^2-N$, $T=mY-hZ$ and
divide the triple by its common gcd. Direct substitution gives the equation.
Before division $v_p(Y)=e$, $h\mid Z$, and $h^2\mid T-mY$, so
after division $T=aY\bmod h$ and $Y$ is a unit there.
At every other prime power dividing $N$, $Y=h^2\bmod q^f$ is a unit and
$T=-aY\bmod q^f$ by substitution. These congruences give the stated gcd.
The raw coordinates have magnitude $O(N^3)$. Conversely any primitive point
whose difference $T-aY$ has a proper gcd with $N$ gives a factor immediately.
Thus a single positive definite binary form with a variable square on the
right removes the field-selection cost. Producing a point with the required
gcd remains the factoring search itself. Example: $N=77$, $D=13$,
$(T,Y,Z)=(-149,-2,17)$ gives gcd 7.

### A smaller denominator on the single public conic

Write $a=\lfloor\sqrt N\rfloor$, $D=N-a^2>0$, $\gcd(a,N)=1$ as above.
For every odd proper prime divisor $p>2\sqrt{D/3}$ there is a point with
$1\le Z\le2\sqrt{D/3}$ and $\gcd(T-aY,N)=N/p$, even if $p^2\mid N$.
For existence take a shortest nonzero vector $(t,y)$ in
$t=-ay\bmod p$, with squared length $t^2+Dy^2$. The planar Hermite bound
gives $t^2+Dy^2=cp$ with $1\le c\le2\sqrt{D/3}<p$.
Thus $p$ divides neither $t$ nor $y$. Put $s=\sqrt{-D}$ and set
$T+Ys=(a+s)(t+ys)^2/p$. Both coefficients are integral because
$(a+s)(t+ys)\equiv y(a+s)(s-a)=-Ny\equiv0\bmod p$.
Taking norms gives $T^2+DY^2=Nc^2$, and direct subtraction gives
$T-aY=-2(N/p)ty$. Hence $Z=c$ has the claimed properties.
If some common coordinate factor also divides $N$, its gcd with $N$ already
splits the input (it is at most $c<p$); otherwise primitive reduction preserves
the extracted divisor. Since $D\le2a$, the denominator bound is $O(N^{1/4})$.
The short-vector lattice in the proof uses $p$. Bounding $Z$ alone gives
no procedure to find $T,Y$ from $N$ within polynomial bit time.

### Harmonic theta coefficient for the denominator-one circle

For the preceding $N=77$ example an integral point is $(-5,-2,1)$:
$25+13\cdot4=77$ and $\gcd(-5+16,77)=11$.

For $N=a^2+1=pq$ odd with distinct primes, define the fourth harmonic
moment $H_N=\frac18\sum_{x^2+y^2=N}\operatorname{Re}((x+iy)^4)$, with the
real part taken after the fourth power. Unique factorization in the Gaussian
integers gives exactly two representations up to sign changes and interchange:
each prime splits, and conjugating just one prime factor gives the other orbit.
There are eight points in each orbit. One representative is $(a,1)$; call the
other $(t,y)$ with $t>y>0$. Since the summand equals $N^2-8x^2y^2$,
$$ (ty)^2=(2N^2-8a^2-H_N)/8. $$
Square roots then give $t+y$ and $t-y$ from $N\pm2ty$, hence both coordinates.
The ratios $t/y$ and $a$ modulo $N$ agree at precisely one prime up to a
global sign, by the two Gaussian-factor choices; $\gcd(t-ay,N)$ splits $N$.
Thus an exact single coefficient of the harmonic theta series suffices here.
The ordinary circle count is always 16 on this family. Extracting the harmonic
coefficient at a composite index in polynomial bit time is the missing step;
the direct sum enumerates the circle points. This observation is restricted
to the stated input family.

### Coarse harmonic readings at several degrees

On the same two-representation family use the normalized sum
$\frac18\sum_{x^2+y^2=N}\operatorname{Re}(((x+iy)/\sqrt N)^{4m})$.
After subtracting the publicly known orbit it gives $c_m=\cos(m\theta)$,
where $\theta=4\arctan(y/t)\in(0,\pi)$. Allow total absolute error 0.001
per residual reading. Only $O(\log N)$ readings with $m=O(\sqrt N)$ suffice.
The known orbit at binary-encoded $m$ is computable by repeated squaring
with $O(\log N)$ bits of working precision.

To see this, find a power of two $h$ with observed $|c_h|\le0.8$ by doubling.
The distance of $\theta$ to $\{0,\pi\}$ is $\Omega(1/\sqrt N)$ because
$y\ge1$ and $t-y\ge1$ (the latter gives
$\tan((\pi-\theta)/4)=(t-y)/(t+y)$). Doubling that endpoint distance
until it first reaches $\pi/4$ gives an angle with absolute cosine at most
$1/\sqrt2$, including the case already between $\pi/4$ and $\pi/2$.
Thus $h=O(\sqrt N)$. Choose the global sign of the angle so $\sin(h\theta)>0$.
For $k=1,2,4,\dots$ estimate $\sin(k\theta)$ by
$(c_{|k-h|}-c_{k+h})/(2\sqrt{1-c_h^2})$ using the noisy readings.
The true and estimated denominating square roots exceed 0.59 and differ by
less than 0.002, so the sine error is below 0.006. Together with $c_k$ this
determines $k\theta$ modulo $2\pi$ within 0.02 radians. At each doubling
lift the measured angle to the multiple-of-$2\pi$ branch nearest twice the
previous lifted angle; the accumulated discrepancy before choosing the branch
is at most 0.06 radians. Dividing the final lift by $k\ge\sqrt N$ gives
angle error at most $0.02/\sqrt N$. Fold to $[0,\pi]$, divide by four,
and round the radius-$\sqrt N$ Cartesian coordinates to recover $t,y$.
This trades coefficient precision for growing harmonic degree. The unknown
operation remains evaluating the normalized sums from $N$ in polynomial bit time.

### Precision barrier for real heat samples, even on circle radii

Consider a black-box radial sum $F(t)=\sum_n b_n e^{-tn}$ for real $t\ge0$.
Suppose only that the shell weights are bounded and supported on sums of two
squares. At $N=a^2+1$ there are two such nonnegative spectra, with target
coefficients differing by 1, whose heat readings differ by at most
$4^{-2m}$ uniformly in $t$, where $m=\lfloor a/32\rfloor\ge1$.
To construct their difference put alternating binomial weights
$(-1)^{j+m}\binom{2m}{j+m}/\binom{2m}{m}$ at the distinct radii
$(a+j)^2+1$, $-m\le j\le m$, and split into positive and negative parts.

For $f(x)=\exp(-t((a+x)^2+1))$ the resulting heat difference is the
normalized $2m$th forward difference at $-m$. On the complex disk of
radius $a/4$ about any real $x\in[-m,m]$ one has
$\operatorname{Re}((a+z)^2+1)>0$, hence $|f(z)|\le1$. Cauchy's bound
and the repeated unit-interval integral formula for finite differences give
an absolute difference at most $(2m)!(4/a)^{2m}/\binom{2m}{m}
\le(8m/a)^{2m}\le4^{-2m}$.
Thus even adaptive real heat queries with independent worst-case absolute
errors of half this size can have identical transcripts on these two spectra.
Recovering the target uniformly to error below $1/2$ needs
$\Omega(\sqrt N)$ bits of absolute precision. This is a black-box bound
with the actual sum-of-two-squares support restriction. The constructed weights
are not asserted to be the harmonic theta coefficients; arithmetic relations
between those coefficients and complex-time queries fall outside the bound.

### Variance of uniform complex-phase extraction on the true circle weights

Let $b_n$ be the normalized degree-$4m$ angular reading on shell $n$,
with axis points and diagonal points removed. Those points do not contribute
on the target family of odd nonsquare semiprimes. Set $F(z)=\sum_{n\ge1}b_nz^n$.
For uniform $U$ on the unit circle and any $0<\rho<1$, the complex estimator
$F(\rho U)U^{-N}/\rho^N$ has mean $b_N$. Orthogonality gives complex
mean-square error $\sum_{n\ne N}|b_n|^2\rho^{2(n-N)}$ for one sample.

Set $\theta=4\arctan(1/2)$. For each prime $l\equiv3\pmod4$ every
representation of $5l^2$ or $25l^2$ scales from 5 or 25, since $l$ divides
both coordinates. The sole off-axis orbits at 5 and 25 are generated by
$2+i$ and $(2+i)^2=3+4i$ respectively. Consequently the two coefficients
are $\cos(m\theta)$ and $\cos(2m\theta)$, at least one of absolute
value at least $1/2$ by the double-angle identity. For $l<\sqrt N/5$
both indices are below $N$, making their radial weights at least one.
All these indices are distinct. The prime number theorem in progressions
gives a complex mean-square error of $\Omega(\sqrt N/\log N)$, uniformly
in the harmonic degree and radius. Independent uniform phases averaged
arithmetically divide this error by the sample count. Thus even with exact
evaluations and axes and diagonals removed, this particular estimator needs
$\Omega(\sqrt N/\log N)$ samples for fixed complex RMS accuracy.
This variance bound does not establish a tail bound or cover structured
quadrature, other estimators, or additional arithmetic subtractions.

### Randomly rotated Fourier grids after removing imprimitive points

On the odd squarefree target family $N=a^2+1$, retaining only points
with coprime coordinates preserves all target points. Let $b_j(n)$ denote
the resulting normalized angular coefficient at degree $4j$.
Average the complex Fourier estimator over $Q$ equally spaced phases on
a common radius with a uniform random initial rotation. The geometric sum
eliminates all indices except $n\equiv N\pmod Q$. Orthogonality in the
rotation leaves variance $\sum_{n\ne N,\ n\equiv N\,(Q)}
|b_j(n)|^2\rho^{2(n-N)}$.

For the two degrees indexed by $m,2m$, take primes $p<N$ congruent to $N$
modulo $L=\mathrm{lcm}(4,Q)$, with $\gcd(Q,N)=1$. Each is $1\bmod4$,
has one primitive orbit, and contributes $\cos^2 u+\cos^2(2u)\ge7/16$
to the sum of variances, since $z+(2z-1)^2\ge7/16$ for $z=\cos^2 u$.
Both radial weights are at least one, even with different radii for the
two degrees. For $Q\le(\log N)^A$ with fixed $A$, Siegel-Walfisz gives
$\gg_A N/(\varphi(L)\log N)$ such primes uniformly for large $N$.
Averaging $B$ independent rotations of this grid for each degree divides
the sum of errors by $B$. Thus fixed joint complex RMS accuracy requires
$BQ=\Omega_A(N/\log N)$ phase evaluations up to a constant factor, since
$\varphi(L)\le4Q$. This bound concerns arithmetic means on a common grid
size at the two degrees and complex mean-square error. It neither assumes
noise in a theta evaluation nor covers arbitrary deterministic quadrature.

### Direct circle evaluator by rational tangents

[circle_tangent_arcs.py](circle_tangent_arcs.py) lists all representations
in $0\le x\le y$ from $N$ alone in $O(N^{1/3}\operatorname{polylog}N)$ bit
operations. Let $R=\sqrt N$, take block width $H\asymp R^{1/3}$ with
$H\le R^{1/3}$, and write $f(x)=\sqrt{N-x^2}$. At each block start $b$
the last rational convergent $p/q$ to $b/\lfloor f(b)\rfloor$ with $q\le H$
has $|q b/\lfloor f(b)\rfloor-p|\le1/H$. Replacing the floor by the true
value changes the slope by at most $2/R$ since $b\le f(b)$ and
$\lfloor f(b)\rfloor\ge1$ (equivalently use $b\le\lfloor f(b)\rfloor$).
On the sector $|f''|\le2\sqrt2/R$. Taylor's theorem bounds the variation
of $px+qf(x)$ from its starting value by $1+2H^2/R+\sqrt2 H^3/R<5$.
Only constantly many integer levels need testing. On a level $px+qy=L$,
$(qx-py)^2=(p^2+q^2)N-L^2$ gives both possible points by an exact
square test and two divisions. There are $O(R/H)$ blocks.
The list suffices for angular evaluation, or directly for the gcd split on
the preceding two-representation family. This improves the coordinate-scan
bound for this evaluator. It remains exponential in input bit length.

### Complexity baseline for circle enumeration and angular evaluation

For $N=a^2+1=pq$ with distinct odd primes, finding the other orbit and
factoring are deterministically equivalent up to polynomial-in-bit-length
postprocessing. Given $p,q$, use CRT for $u=a\bmod p$, $u=-a\bmod q$.
Gauss-reduce the basis $(N,0),(u,1)$. Its shortest nonzero vector $(t,y)$
has squared norm at most $2N/\sqrt3<2N$, and the norm is divisible by $N$,
so it equals $N$. Squarefreeness gives $\gcd(y,N)=1$. Hence the original
gcd decoder recovers $p$. The other direction already follows from any
representative of the second orbit. With both orbits known the normalized
binary-degree angular sums to fixed precision take polynomial bit time by
repeated squaring. In particular the deterministic exponent-one-fifth
factoring algorithm gives an $N^{1/5+o(1)}$ upper bound for this evaluator.
The direct arc algorithm improves on coordinate scanning but is asymptotically
slower than this existing factoring-first route.
[Harvey, Theorem 1.1](https://arxiv.org/abs/2010.05450)

On unrestricted $N$, listing all circle points cannot take polynomial time
in the bit length of $N$ because of output size. For the product of the
first $k$ primes equal to 1 modulo 4 there are $2^{k-1}$ symmetry orbits
by Gaussian factorization, while $\log N=O(k\log k)$ by the prime number
theorem in progressions. This obstruction concerns enumeration only.
For example, on this squarefree family, given a Gaussian prime
$\pi_j$ above each rational factor $p_j$, the signed and ordered angular
sum divided by four equals $\prod_j 2\operatorname{Re}
((\pi_j/\sqrt{p_j})^{4m})$ by independent conjugate choices. It can be
evaluated to fixed absolute precision with polynomially many bit operations
using guard bits for the at-most-$2^k$ product amplification, without listing
the exponentially many orbits.

### Fixed-degree, one-query amplification on unrestricted indices

Define $C(n)=\frac14\sum_{x^2+y^2=n}\operatorname{Re}((x+iy)^4)$
and $A(n)=C(n)/n^2$. The fourfold unit symmetry makes $C(n)$ integral.
Gaussian factorization gives $A(uv)=A(u)A(v)$ for coprime indices,
including vanishing coefficients. For a split prime $r=x^2+y^2$,
$A(r)=2(1-8x^2y^2/r^2)$. A fixed positive proportion of primes
have $A(r)\ge3/2$, by Hecke's prime number theorem for Gaussian primes
in fixed angular sectors (take a narrow sector next to the real axis).

Given $n$, take distinct such primes coprime to $n$ until the exact rational
gain $G=\prod A(r)>n^2$, and let $M$ be their product. This needs
$O(\log(2n))$ primes. At most $\log_2 n$ primes are skipped for dividing
$n$; the fixed-sector prime count ensures scanning to $O((\log(2n))^2)$
asymptotically suffices. Primality, circle enumeration for those small primes,
and all rational comparisons therefore cost polynomial bit time.
The query $nM$ has $O(\log(2n)\log\log(3n))$ bits. One reading of $A(nM)$
with absolute error at most $1/4$ recovers $C(n)$ by multiplying by $n^2/G$
and rounding. This works also at zero or negative coefficients. On the
two-representation family the earlier exact angular moment is $H_N=C(N)/2$.

The query uses the sum divided by four, with no division by the number of
representations. The growing magnitude on the padded input is essential;
a fixed relative error or a constant-error orbit average does not obtain
this rounding guarantee. The evaluator must accept the padded indices,
which generally leave the original semiprime family.

### Two fixed-degree readings separate the hidden local denominators

**Prior result.** The two-index elimination and denominator gcd below
instantiate the Bach–Charles proof strategy (Section 2, equations 1–3).
The repo previously recorded the non-CM density-one theorem and implemented
the weight-two reduction. The Gaussian calculation here checks the odd-weight
CM instance explicitly, including its exact inert-prime exception. The earlier
padding lemma addresses the absolute precision of an unrestricted evaluator.
No composite-index coefficient algorithm follows from either reduction.
Reopen this route on evidence for an evaluator, rather than another change
of eigenform yielding the same elimination.
[Paper](https://arxiv.org/abs/0708.1192) ·
[Prior survey](../../notes/polyimpl_results.md) ·
[Composite-index evaluator audit](../../notes/nativemod_results.md) ·
[Weight-two implementation](../../src/tracedata/reductions.py)

Let $N=pq$ for distinct odd primes, with no supplied circle point. Write
$v=A(N)$ and $w=A(N^2)$ using the preceding normalization.
For a split prime $p=x^2+y^2$, $a_p=A(p)=2-16x^2y^2/p^2$ has reduced
denominator exactly $p^2$. The three Gaussian-factor choices at $p^2$
give $A(p^2)=a_p^2-1$. For an inert prime the corresponding readings
are $A(p)=0$, $A(p^2)=1$.

If both primes split, put $X=a_p^2$, $Y=a_q^2$. Multiplicativity yields
$XY=v^2$ and $X+Y=v^2+1-w$. Solve this rational quadratic using an exact
square root of its discriminant. The reduced denominators of its roots are
$p^4$ and $q^4$; a gcd with $N$ splits it. In the mixed case $v=0$ and
$w=a_s^2-1$ for the split prime $s$, so $\gcd(N,\operatorname{den}(w))=s$.
If both primes are inert then $(v,w)=(0,1)$, and more generally
$A(N^j)$ is zero for odd $j$ and one for even $j$. Consequently this particular
power sequence carries no factor information in that case.

The amplification lemma recovers each of the two exact readings with one
fixed-absolute-error query at a padded input. Coverage now includes every
odd distinct-prime semiprime with at least one prime equal to 1 modulo 4.
The missing algorithm is for the angular evaluator itself.

### What a low-bit coefficient evaluator would actually compute

For odd $n$ put $R=r_2(n)/4$, $s=1$ for a square and zero otherwise,
and $h=(R-s)/2$. There are $h$ off-axis orbits of size eight. Their coordinate
parities are opposite. The integer fourth harmonic coefficient therefore satisfies
$$C(n)=n^2R-16\sum_{\rm off\text{-}axis\ orbits} x^2y^2
\equiv n^2 R-64\,1_{n\equiv5\bmod8}\,h\pmod{256}.$$
Indeed $x^2y^2$ is 4 modulo 16 if the even coordinate is 2 modulo 4
and zero otherwise. The first case occurs exactly for $n=5\bmod8$.
In particular $C(n)=n^2 R\bmod64$; the bound is sharp already at $n=5$.
Conversely $C(n)=(n^2-32I)R+32Is\bmod256$ for that indicator $I$.
The coefficient of $R$ is odd and invertible. Thus these eight bits are
efficiently interconvertible with the ordinary representation count modulo
256 from $n$ alone; they contain no additional angular information.

Even a modulo-8 evaluator has a specific arithmetic burden. For $n=pq=1\bmod4$
with distinct odd primes, $C(n)\bmod8$ is 4 if both primes split and zero
if both are inert, by the two-squares theorem. This decides whether $-1$
is a square modulo $n$ on that promise. Verified the full modulo-256 identity
on all 5,000 positive odd inputs below 10,000 by direct enumeration.
This checks a low-precision modular starting point; it supplies no general
polynomial-time composite-index evaluator.

### The ninth bit detects simultaneous quadratic residuosity of 2

Suppose $N=pq\equiv1\pmod8$ with distinct $p,q\equiv1\pmod4$. Then
$$C(N)\equiv4N^2+256\,1_{p\equiv q\equiv5\bmod8}\pmod{512}.$$
For $p=x^2+y^2$ with even $y$, write $C(p)=2p^2-64u_p$ where
$u_p=(xy/2)^2$. This square is 0 modulo 4 when $p=1\bmod8$ and 1 when
$p=5\bmod8$. Multiply the local formulas and use $p^2=q^2=1\bmod4$
and the equality of the two primes' residues modulo 8. By the supplement
to quadratic reciprocity, the zero correction is precisely the case where
2 has a square root modulo $N$. On this family $r_2(N)/4=4$ identically.

For every fixed $b\ge9$, both outcomes occur with $N=1\bmod 2^b$.
Choose the two prime residues modulo $2^b$ to be $1,1$ or $5,5^{-1}$;
Dirichlet's theorem supplies infinitely many distinct choices.
Thus the eight-bit count equivalence is sharp for any rule using only
the count and a fixed-length binary suffix of $N$. This gives a test of
a proposed evaluator, not an efficient construction of that bit.

### Direct evaluation on a finite complex grid

[circle_fourier_grid.py](circle_fourier_grid.py) takes an unfactored $N$
and computes the fourth angular coefficient. Truncate each one-dimensional
theta polynomial $T_j(z)=\sum_{|x|\le\sqrt N}x^j z^{x^2}$ for
$j=0,2,4$, with constant term 1 for $T_0$. Expanding the fourth power gives
$C(N)=[z^N](T_0 T_4-3T_2^2)/2$. This polynomial has degree at most $2N$.
On the $(N+1)$st roots of unity the Fourier bin of index $N$ has no other
representative in $[0,2N]$. Averaging its demodulated samples therefore
gives the desired integer exactly in exact arithmetic.

There is no exponential precision cost on this grid. For example, evaluate
each unit phase (including the demodulator) to absolute error
$\epsilon=(10000(N+1)^4)^{-1}$ and round subsequent operations to a fixed
absolute grid of that spacing. The theta sums have magnitude at most
$3(N+1)^{(j+1)/2}$, and termwise phase errors and additive roundoffs
sum to at most $10\epsilon(N+1)^{(j+1)/2}$ for each. The product,
demodulation and averaging errors are consequently below $1/4$. Both
integer and fractional bit requirements are $O(\log N)$. The reference
uses mpmath guard bits and conjugate symmetry; its rounding was compared
with exact coordinate sums on 69 inputs. It uses $\Theta(N)$ grid points
and $O(N^{3/2})$ phase terms; FFT convolution can compute the same bin with
$O(N\log N)$ arithmetic operations. This is a direct evaluation baseline,
slower than the tangent-arc method. The grid count is a cost of this
construction, not a general lower bound for complex evaluation.

## Reproduce the principal results

Run these commands in this directory. The lattice and tangent programs use SymPy. The direct disk enumerators are controls with power-of-$N$ running time.

```sh
python fixed_discriminant_split.py
python cyclic_class_queries.py
python affine_parity.py
python bounded_quadratic_roots.py
python bernoulli_moment.py
python reciprocal_range.py
python stationary_log.py
python stationary_threshold.py
python stationary_dual.py
python stationary_cubic.py
python native_phase_evaluator.py
python quartic_sampling_check.py
python stationary_window_reference.py
python polynomial_stream_limits.py
python small_scale_phase_check.py
python idempotent_count.py
python count_precision.py
python idempotent_residues.py
python idempotent_zero_bits.py
python selector_digits.py
python phase_moment_control.py
python crt_range.py
python stationary_hardness.py
python tangent_progressions.py 38942974393
```
