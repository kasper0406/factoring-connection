# Can short matrix products expose the first nonunit?

**Candidate and verdict.** For distinct primes \(p<q\), let \(N=pq\), \(B=\lfloor\sqrt N\rfloor\), and \(n=N.\mathrm{bit\_length}()\). The interval-product reduction gives \(\gcd(B!,N)=p\), since \(p\le B<q\). A weaker output \(u(B)\) with \(u(B)=0\pmod p\) and \(u(B)\ne0\pmod q\) also suffices, with one gcd and one exact division. We test two linear representations for computing all cutoffs up to \(B\) from public \(N\). The first fails even for a zero test, and the second fails for exact factorial residues. A uniform constructor for the weaker binary output remains open. The underlying gcd reduction and the question of short factorial straight-line programs are established prior work. [Track C](trackC_analysis.md) · [Straight-line programs, §§1.5 and 4.6](slpcirc_results.md) · [Sufficiency map, E4](suffmap_results.md) · [Lipton, *Factoring and Factorials*](https://rjlipton.com/2009/02/23/factoring-and-factorials/)

## Fixed linear recurrences for the zero test

**Proved.** Let \(B\ge p\). Suppose \(u(k)=L A^k R\) over \(\mathbb F_p\), for a square \(d\)-dimensional matrix \(A\), row \(L\) and column \(R\). Require \(u(p-1)\ne0\) and \(u(k)=0\) for \(p\le k\le B\). Then
\[
d\ \ge\ s=\min(p,B-p+2).
\]
Indeed \(H_{ij}=u(p-s+i+j)\) for \(0\le i,j<s\) lies in the specified range \(0,\ldots,B\). Its entries on the antidiagonal equal \(u(p-1)\), and every entry beyond that antidiagonal vanishes. Reversing the columns gives a triangular matrix with nonzero diagonal. Hence \(\operatorname{rank} H=s\). But \(H_{ij}=(LA^{p-s+i})(A^jR)\) has rank at most \(d\). If the transitions and readouts have period \(h\), a cyclic clock of \(h\) blocks lifts them to one constant matrix of dimension \(hd\), so \(hd\ge s\). The conclusion applies to matrices over \(\mathbb Z/N\mathbb Z\) by reduction modulo \(p\), even if the chosen matrices depend on \(N\) and \(B\).

**Refuted for this model.** Powering an explicit, polylogarithmic-dimensional matrix would use polynomial bit time, since it takes \(O(\log B)\) matrix multiplications. For any fixed \(0<\eta<1\), on inputs with \((1+\eta)p\le q<2p\), the bound is \(\Omega_\eta(p)\). For twin primes \(B=p\) and the bound is only two. The claim concerns representations of *all* the queried sequence values by a common recurrence. It gives no lower bound for an algorithm required to return just the single output at \(B\), or for the cost of an implicit high-dimensional operation.

## Binary matrix products for exact factorial residues

**Proved.** Read the padded bits of \(k\) from high to low, or from low to high. Multiply one matrix for each bit, allowing its entries and dimensions to depend on the layer, its digit, \(N\), and \(B\). Attach a fixed input and output vector. If the scalar output is \(k!\pmod p\) for \(0\le k\le B\), the width at a cut between consecutive bits is at least
\[
m=2^{\lfloor\log_2\sqrt{\min(p,B+1)}\rfloor}.
\]
Use the cut just above the low \(\log_2 m\) bits. Fix the unused high bits at zero and set \(k=mj+i\) for \(0\le i,j<m\). Splitting the product at the cut factors the \(m\)-by-\(m\) table through its width. The true table \(K_{ij}=(mj+i)!\) has determinant
\[
\det K =
\left(\prod_{j=0}^{m-1}(mj)!\right)
\prod_{0\le a<b<m}m(b-a)\ne0\pmod p .
\]
Indeed every column factor \((mj)!\) is nonzero. What remains in row \(i\) evaluates the monic degree-\(i\) polynomial \(\prod_{t=1}^i(X+t)\) at the distinct points \(X=mj\). The unit triangular change to the monomial basis preserves the Vandermonde determinant. Thus the width is at least \(m\), which exceeds \(\sqrt p/2\) when \(B\ge p\).

**Refuted for this model.** Explicit binary linear representations of exact factorials require exponential width in \(n\) on balanced semiprimes. This leaves variable reordering, repeated passes, nonlinear arithmetic, a single cutoff, and zero-only outputs open. The stated rank arguments do not lower-bound general arithmetic circuits for factorials.

## A two-state comparator and its construction cost

**Proved.** For a known prime threshold \(r\), the binary function \([k<r]\) has width two. Start with \((l,e)=(0,1)\) indicating that the prefixes so far are less or equal. Reading input bit \(x\) and border bit \(b\) from high to low gives
\[
l'=l+[x<b]e,\qquad e'=[x=b]e .
\]
The output is \(l\). Applying this separately at \(p,q\) and combining the border bits by CRT would supply one gcd-revealing output at \(B\). The missing public input is the sequence of combined border bits. For the input digit one the layer matrix is \(\operatorname{diag}(1,b)\). At any position where the two prime bits differ, the determinant's gcd with \(N\) is a proper factor. Invertible changes of basis at either end of a layer multiply the determinant by units and preserve this gcd.

**Open.** A factor-blind polynomial-time compiler of a zero-only digit program would suffice even if its output below each local prime were an arbitrary nonzero residue. For the single-query reduction it would suffice to separate the two primes just at \(B\). Width two with factor-supplied bits establishes representability only. A different construction needs public transitions and a proof of distinct local vanishing at the actual cutoff. A proposal supplying CRT-packed threshold bits or conjugates of the explicit comparator must already produce a proper divisor at a layer. Neither lower bound above obstructs a different zero-only digit program.

## A single central coefficient that factors every balanced input

**Proved.** At the one public cutoff \(B=\lfloor\sqrt N\rfloor\), on odd primes \(p<q<2p\),
\[
\gcd\!\left(\binom{2B}{B},N\right)=q .
\]
Since \(B<q\le 2B<2q\), the numerator \((2B)!\) contains exactly one multiple of \(q\), and neither denominator factorial contains one. For \(p\), the inequalities \(p\le B<\sqrt{2}\,p<3p/2\) place \(p,2p\) in the numerator and exactly one \(p\) in each denominator. All these multiples are below \(p^2\) because \(p\ge3\). Thus the central coefficient has \(q\)-valuation one and \(p\)-valuation zero. This changes the repository's central-binomial observable from \(\binom{2N}{N}\) to \(\binom{2B}{B}\), reduces the polynomial degree from \(2N\) to \(2B\), and extends its documented guarantee from \(q<3p/2\) to \(q<2p\). [Native invariant V4](nativeinv_results.md) · [Family 54](newfamily_results.md)

**Refuted for the tested readout.** Compute \((1+x)^{2B}\bmod(N,x^d-1)\) by binary powering, then take the coefficient at \(B\bmod d\). For \(d=1,2\) this returns \(2^{2B}\) and \(2^{2B-1}\) respectively, always units for odd \(N\). More strongly, on \(N=503\cdot937=471311\) *every* width \(1\le d\le128\) gives a unit at that position. The exact central coefficient has gcd 937. A second control, \(1009\cdot1877\), also gives a unit at every one of these widths. The harness verifies the quotient coefficient exactly by summing the original integer binomial coefficients whose exponents agree modulo \(d\). On other fixtures individual widths occasionally reveal a proper divisor. These two counterexamples reject a guarantee from the tested bounded family of unweighted middle-position folds. They leave larger, weighted or adaptively combined projections open.

**Priority.** Repeated squaring computes the circuit for \((1+x)^{2B}\) using \(O(\log N)\) gates. Computing the one middle coefficient modulo \(N\), or any residue with the same local zero pattern, is the unresolved step. The integer coefficient already has \(\Theta(B)\) bits. A proposal should specify how it obtains this local zero pattern without materializing the coefficient, traversing \(B\) entries, or relying on known local characteristic. The degree reduction alone does not justify prioritizing the central coefficient over a mechanism with an independently efficient public extraction operation.

## Exact checks

**Experimentally supported.** The factor-aware harness uses five fixed prime pairs, including a twin-prime control. It computes the public reference product by \(B\) modular multiplications and one gcd. It independently takes ranks over \(\mathbb F_p\) for factorial values and for a different nonzero labelling before the same threshold. It verifies each binary determinant directly and exhausts the packed comparator at every input of the chosen bit length. All five final gcds and all five comparator layer divisors were structural. These are deterministic identities, with no random-gcd sampling.

| \((p,q)\) | \(B\) | Unary rank, both labellings | Binary exact rank |
|---|---:|---:|---:|
| (101, 103) | 101 | 2 | 8 |
| (101, 179) | 134 | 35 | 8 |
| (251, 467) | 342 | 93 | 8 |
| (503, 937) | 686 | 185 | 16 |
| (1009, 1877) | 1376 | 369 | 16 |

With python-flint 0.9.0 installed, reproduce from the repository root using the command below. The matrices are dense verification objects with respectively \(s^2\) and \(m^2\) field entries. The factor-aware comparator is a diagnostic control. [Diagnostic code](../experiments/direct_split_poly/check.py) · [Raw records](../experiments/direct_split_poly/results.jsonl)

```sh
python experiments/direct_split_poly/check.py | diff experiments/direct_split_poly/results.jsonl -
```
