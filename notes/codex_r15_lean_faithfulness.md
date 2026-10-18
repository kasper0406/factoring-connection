> Source: GPT-6.1-sol r15 Lean faithfulness review, 2026-10-18 (run against the pre-fix commit f77e05d).

The algebraic core is faithful. I found **no substantive cast/subtraction error, globally vacuous main theorem, or circular extraction proof**. The qualifications are quantitative: Lemma 8 loses the exact logarithm, Lemmas 9–10 formalize only parts of their paper statements, and `ArithmeticLargeSieve` is a valid prime-term corollary of the classical theorem.

Below, **STRONGER** means a valid generalization. **WEAKER** identifies a lost paper-level conclusion. Related helper names share rows; every theorem in the seven files is covered.

| Lemma | Lean name | Verdict | Explanation |
|---|---|---|---|
| Exponent identity | `exponent_identity` | FAITHFUL | Exact identity over ℤ; subtraction preserves signs. |
| Fermat consequence | `pow_exponent_identity_zmod`, `pow_exponent_identity_modEq` | STRONGER | Correctly permits every base, including zero. Natural exponents and `a^p = a` justify this extension. |
| Fermat for signed exponents | `units_zpow_eq_of_sub_dvd` | FAITHFUL | Uses units, integer exponents, and divisibility by `(p : ℤ) - 1`. Negative exponents have their intended meaning. |
| Coverage criterion | `coverage_criterion` | STRONGER | Only `ℓ ∤ p` is needed; the note assumes `ℓ ∤ N`. Cancellation is legitimate. |
| Coverage with note’s hypothesis | `coverage_criterion'` | FAITHFUL | Exactly supplies cancellation from `ℓ ∤ pq`. |
| Root representative invariance | `root_iff_of_modEq` | FAITHFUL | Correctly transfers the quadratic congruence between representatives. It does not implement root enumeration. |
| Lemma 1 | `lemma1_class_lock` | FAITHFUL | Root, inverse-lift, representative, and coverage congruences are explicit. Conclusion is integer divisibility by `ℓ²`, avoiding any quotient-rounding ambiguity. |
| Distinct classes: determinant | `same_class_dvd` | FAITHFUL | Correct determinant congruence; the necessary cancellation hypothesis `ℓ ∤ q` is present. |
| Distinct primitive weights | `weights_eq_of_same_class` | FAITHFUL | Positivity, height bounds, both coprimality conditions, and `H² < ℓ` are present. These exclude proportional but unequal weights. |
| Lemma 2: lower bound | `lemma2_lower` | FAITHFUL | `Nat.sqrt (4ABN)` is precisely the required floored square root. No truncated subtraction occurs. |
| Lemma 2: upper bound | `lemma2_upper` | STRONGER | Valid even for negative `B`: the linear term does not affect convexity. Positivity of denominators follows from `0 < a ≤ x ≤ b`. |
| Lemma 2: window conversion | `lemma2_window` | FAITHFUL | For positive `K`, Lean’s integer division is floor division. `-((c-Lo)/K)` correctly expresses the ceiling, including negative numerators. |
| Lemma 2: indexing | `lemma2_index` | FAITHFUL | Euclidean remainder gives `0 ≤ i* < m`. Nonnegativity and the upper bound for `j*` are correctly conditional on window membership. |
| Lemma 3 | `lemma3_core`, `lemma3` | FAITHFUL | Derives the collision from the class-lock equation and Fermat. No order assumption is hidden. Arbitrary units validly generalize base 2. |
| Lemma 4 | `lemma4_core`, `lemma4` | FAITHFUL | Correct discrepancy `(A-B)(p-1)`, with the correct sign and modulus. |
| Lemma 4: exponent difference | `lemma4_exponent_difference` | FAITHFUL | Exact integer divisibility by `q-1`; no natural subtraction. |
| Lemma 5 | `lemma5_one_sided` | FAITHFUL | Distinct prime factors, one collision modulo `p`, and absence of every baby collision modulo `q` imply exactly gcd `p`. Integer-product evaluation preserves the modular meaning. |
| FULL/divisibility translation | `not_full_iff` | FAITHFUL | Exact translation of modular equality into divisibility of an integer difference. |
| Lemma 6: four bands | `lemma6_clean_band` | STRONGER | The larger admissible bound `n < (Y+1)^4` is valid because natural band primes exceeding `Y` are at least `Y+1`. Positivity of `n` is present. |
| Lemma 6: prime specialization | `lemma6_clean_band_prime` | FAITHFUL | The paper’s hypotheses imply `q ≤ Y⁴`. Primality ensures `q-1 > 0`, so natural subtraction is harmless. Empty bands being clean is legitimate. |
| Coprime cancellation | `coprime_dvd_mul_iff`, `coprime_dvd_mul_iff_int` | FAITHFUL | Correct cancellation in ℕ and ℤ. |
| Cleanliness to coprimality | `clean_coprime`, `coprime_of_clean` | FAITHFUL | Correctly transfers cleanliness through the product, its divisors, and divisors of `q-1`. |
| Lemma 7: group core | `lemma7_core` | FAITHFUL | Correct order-divisibility equivalence after cancelling `M`. Valid for arbitrary groups; even `orderOf = 0` has consistent integer-divisibility semantics. |
| Lemma 7: FULL equivalence | `lemma7_full_iff` | FAITHFUL | Uses the **exact** shifted interval `[-i*, m-1-i*]`. Omitting `0 ≤ i* < m` is harmless for this translation equivalence; those hypotheses are needed for containment in `(-m,m)`. |
| Lemma 7: modular specialization | `lemma7_coprime_zmod` | FAITHFUL | Correctly uses unit order dividing `q-1`. The natural subtraction is safe because `q` is prime. |
| Prime-product counting helper | `pow_card_le_of_primes_dvd` | FAITHFUL | Distinct primes come from a `Finset`; nonzero `n` is explicitly required. |
| Small order implies Mersenne divisor | `dvd_two_pow_sub_one_of_orderOf_le` | FAITHFUL | Excludes `q=2` and establishes positive order before choosing `k`. Subtraction in `2^k-1` is safe. |
| Mersenne prime-divisor count | `card_mul_lt_of_dvd_mersenne` | FAITHFUL | Correct integer bound `b·card < k`. Its weak threshold `2^b ≤ q` already suffices. |
| Mersenne union bound | `card_le_sum_mersenne` | FAITHFUL | Correctly handles primes appearing for multiple exponents through an upper bound. |
| Gauss sum | `gauss_sum_Icc` | FAITHFUL | Exact natural-number identity, including `T=0`. |
| Lemma 8: “sharp” form | `lemma8_sharp` | WEAKER | Proves `2b·card ≤ T(T+1)` under `q > 2^b`. Taking `b = floor(log₂ Q₀)` loses the paper’s exact real-log denominator. |
| Lemma 8: coarse form | `lemma8` | WEAKER | Correct `card ≤ T²` bound for odd primes, but lacks the paper’s logarithmic saving from the lower threshold. |
| Residue-class counting | `card_le_of_pairwise_modEq` | FAITHFUL | Natural subtraction is protected by interval membership. `M=0` is harmless: congruence modulo zero means equality, hence at most one element. |
| Lemma 10: CRT and single triple | `lemma10_pairwise_modEq`, `lemma10_single` | FAITHFUL | Correct modulus `ℓ·(d/gcd(d,A-B))`; all required cancellation and coprimality hypotheses are present. |
| Lemma 10: modulus lower bound | `lemma10Modulus_ge` | FAITHFUL | Uses `A ≠ B` and integer absolute value. Since the gcd divides `d`, natural division is exact. |
| Lemma 10: union | `lemma10_union` | WEAKER | Faithful union-bound core, but does not instantiate the weight/band/shift family or derive the paper’s final bound using `H,Y,T,m̄`. |
| Excluded residue specification | `excludedResidue_spec` | FAITHFUL | `ℓ ∤ B` makes the modular inverse legitimate; the representative is canonical and satisfies the intended congruence. |
| Weight cardinality | `card_weights_le` | FAITHFUL | Correct upper bound for positive weights in `[1,H]²`. |
| Classical large-sieve input | `ArithmeticLargeSieve` | WEAKER | Correct prime-term inequality; omits composite squarefree terms of the full classical denominator. Detailed assessment below. |
| Lemma 9 | `lemma9_of_largeSieve` | WEAKER | Faithfully proves the conditional product bound with `|V|Σ1/ℓ`. It does not formalize the full `J`, the reciprocal-sum estimate, or the lower bound on `|V|`. |
| Exceptional-set summation | `exceptional_Eb_of_largeSieve` | WEAKER | Correct conditional summation. It leaves positivity and the useful lower bound for `J`, prime counts, and population normalization external. |
| Single-target extraction | `single_target_extraction` | FAITHFUL | Meaningfully chains Lemmas 3, 4, 7, and 5 under residue-formula and arithmetic hypotheses. It is not an executable-algorithm correctness theorem. |

The auxiliary definitions `excludedResidue`, `uncovered`, `lemma10Set`, and `lemma10Modulus` also have the intended meanings. In particular, inspection of the elaborated types confirms that `lemma10Set` uses `(p : ℤ) - 1`, and the capstone uses `(m : ℤ) - 1 - istar`.

**The large-sieve hypothesis is mathematically sound, but is not the full classical statement.** At [LargeSieve.lean:23](/home/knielsen/code/factoring-connection/lean/FactoringLean/LargeSieve.lean:23), it asserts
\[
|U|\left(1+\sum_{\ell\in\mathrm{Ps}}
 \frac{\omega(\ell)}{\ell-\omega(\ell)}\right)\le N_0+Q^2.
\]
The canonical representatives ensure that cardinalities count distinct residue classes; `ω(ℓ)<ℓ` ensures strictly positive real denominators.

For \(Q\ge1\), the classical arithmetic form instead has the squarefree-modulus sum
\[
J(Q)=
\sum_{\substack{1\le d\le Q\\d\ \mathrm{squarefree}\\
                 \text{prime factors of }d\in\mathrm{Ps}}}
\prod_{\ell\mid d}\frac{\omega(\ell)}{\ell-\omega(\ell)}.
\]
Keeping just \(d=1\) and prime \(d\) gives the formal hypothesis. Thus it is a valid weaker bound, not an accidentally stronger assumption. See Montgomery–Vaughan, [Theorem 19.13, printed p. 165](https://personal.science.psu.edu/rcv4/571s25/montgomery-vaughanII.pdf#page=177).

One can also justify the exact formal inequality directly: apply the additive large sieve to frequency zero and all nonzero fractions \(a/\ell\). They are \(Q^{-2}\)-separated; zero contributes \(|U|^2\), and each prime contributes at least \(\omega/(\ell-\omega)|U|^2\). The upper bound is \((N_0+Q^2)|U|\). This uses the [Montgomery–Vaughan inequality, Theorem 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf). The cases `N₀=0` or `Q=0` are harmless elementary degeneracies.

For the intended sufficiently large **single dyadic band**, products of two band primes already exceed \(Q\). Consequently, no supported composite terms occur: the truncation loses nothing there.

**The capstone is substantive and satisfiable.** At [SingleTarget.lean:91](/home/knielsen/code/factoring-connection/lean/FactoringLean/SingleTarget.lean:91), its residue hypotheses specify how to compute the giant; they do not assume the resulting baby collision. The proof then:

1. Derives the collision modulo `p` using Lemma 3.
2. Derives the discrepancy modulo `q` using Lemma 4.
3. Uses Lemma 7 and arithmetic nondivisibility to exclude every baby collision modulo `q`.
4. Applies Lemma 5.

For a concrete nonvacuous instance, take base 2, `p=11`, `q=13`, `M=1`, `ℓ=5`, `A=2`, `B=1`, `c=12`, `g=glo=1`, `jstar=istar=0`, and `m=1`. Then `γ=54` and `G=23` modulo 143. The order modulo 13 is 12, the discrepancy is 10, and
\[
\gcd(G-\gamma^0,143)=\gcd(22,143)=11.
\]
Also, although `A ≠ B` is not explicit, the capstone’s non-FULL hypothesis implies it: `δ=0` lies in the shift interval and would contradict that hypothesis if `A=B`.

Its boundary is clear: the class-lock equation, coprimality, and residue formulas are supplied as hypotheses. It does not establish that the implementation produced these values or visited this giant. The file header’s claim to chain Lemma 1 is therefore inaccurate.

Concrete fixes:

1. **Clarify the large-sieve name and scope.** Rename it `ArithmeticLargeSievePrimeTerms`, or define the full squarefree form and derive the current corollary. Keep its explicitly conditional status.
2. **Add the exact real-log Lemma 8.** Retain `lemma8_sharp` as an integer-threshold corollary; its current “as in the note” description overstates exact agreement.
3. **Assemble the note-level window theorem.** Connect coverage, computed `Lo/Hi`, Lemma 1, and the indexing results. The current Lemma 2 components are sound but separate.
4. **Finish the quantitative Lemma 10 corollary.** Instantiate weights, band primes, and `|δ|<m̄`; apply the modulus bound and cardinality estimates. Handle the paper’s closed interval explicitly—the current counting theorem uses a half-open interval.
5. **Add positive-denominator coverage corollaries.** Require `0<J` when turning the product inequalities into cardinality bounds. Keep the analytic estimates as clearly identified external inputs.
6. **Add a note-level extraction wrapper.** Instantiate base 2, `L=∏ℓ²`, `M=L/ℓ²`, computed residues, and the selected indices; derive `hg` and `hclean`. Correct the header meanwhile.
7. **Narrow the README/results claim.** Say “algebraic core, window components, floor-log small-order bound, and CRT/union counting core”; reserve “fully formalized Lemmas 8–10” for their complete paper-level conclusions.