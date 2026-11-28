# Uniform integer lattices of prescribed determinant

Given positive integers \(D,d\), sample a full-rank sublattice of \(\mathbb Z^d\) uniformly among all those of index **exactly \(D\)**. The output is its upper-triangular column Hermite normal form. The algorithm uses uniform integers and gcds, takes expected polynomial bit time in \(d,\log(2D)\), and averages fewer than 2.295 diagonal proposals. It requires no prime searches. At \(D=91,d=3\), seed zero returns the rows \((7,2,4),(0,1,0),(0,0,13)\). [Sampler](determinant_lattices.py)

The fixed-index ensemble is classical. Zou gives its exact size, and Goldstein and Mayer study its equidistribution. Hu et al. use a determinant **upper bound** in their sampling algorithm and use asymptotic diagonal distributions. Here the input fixes a possibly composite determinant, and the output law is exact. The repository's earlier intermediate-lattice sampler allows all indices between two fixed endpoint lattices. Its parity dichotomy concerns a different law. [Zou](https://arxiv.org/abs/math/0610684) · [Goldstein–Mayer](https://doi.org/10.1515/form.2003.009) · [Hu et al., pp. 2–4](https://eprint.iacr.org/2020/208.pdf) · [Intermediate lattices](sublattice-sampling.md)

## Sampling the diagonal

Write \(a_1,\ldots,a_d\) for positive diagonal entries of product \(D\). Each entry above \(a_i\)'s row diagonal lies in \([0,a_i)\), so this diagonal occurs in \(\prod_i a_i^{d-i}\) Hermite forms. Eliminating \(a_1\) gives multiplicity \(D^{d-1}\prod_{s=1}^{d-1}a_{s+1}^{-s}\).

For each \(s=1,\ldots,d-1\), draw \(s\) independent uniform residues modulo \(D^2\) and set \(a_{s+1}\) to their joint gcd with \(D^2\). Accept if \(\prod_{i=2}^d a_i\mid D\), then set \(a_1=D/\prod_{i=2}^d a_i\). Restart on rejection. Draw each remaining matrix entry independently and uniformly in its specified range.

For any accepted candidate, the valuation of \(a_{s+1}\) at \(p^e\parallel D\) is at most \(e\), strictly below the gcd's saturation value \(2e\). Thus the probability of proposing a particular admissible \(a_{s+1}=a\) is \(a^{-s}\prod_{p\mid D}(1-p^{-s})\). Independence across rows gives exactly the diagonal weights above. Dimension one and determinant one give the unique legal matrix in each case.

Let \(A_d(D)\) count the index-\(D\) lattices. It is multiplicative, and

\[
A_d(p^m)={m+d-1\brack d-1}_p
=p^{m(d-1)}\prod_{s=1}^{d-1}\frac{1-p^{-m-s}}{1-p^{-s}}.
\]

Summing the diagonal proposal probabilities and substituting this identity gives exact acceptance

\[
\prod_{p^e\parallel D}\prod_{s=1}^{d-1}(1-p^{-e-s})
\;\ge\;\prod_{j=2}^\infty\zeta(j)^{-1}.
\]

Hence the expected number of proposals is at most \(\prod_{j=2}^\infty\zeta(j)\simeq2.29486\), uniformly even when the dimension grows. A trial uses at most \(d(d-1)/2\) random residues of \(O(\log(2D))\) bits and as many gcds. The sampler also checks divisibility of the remaining determinant before each assignment. Constructing the final matrix takes \(O(d^2)\) entries. [Fixed-index counting formula](https://arxiv.org/abs/math/0610684)

## What happens at determinant \(N^k\)

Let \(N=pq\) be a product of distinct primes and let \(k\ge1\). CRT identifies global index-\(N^k\) lattices with independent choices at \(p,q\), each of local index \(p^k,q^k\) respectively. This follows, for example, by viewing the quotient as a quotient of \((\mathbb Z/N^k)^d\).

**The HNF diagonal.** For \(k=1\), define \(h_d(x)=\sum_{j=0}^{d-1}x^{-j}\). At each prime \(p\), its diagonal position is \(j+1\) with probability \(p^{-j}/h_d(p)\). The diagonal yields a factor exactly when the two prime positions differ, with probability \(1-h_d(pq)/(h_d(p)h_d(q))\le1/p+1/q\).

For every \(k\), the local probability that the full prime power \(p^k\) lies in the first diagonal entry is \(p^{k(d-1)}/A_d(p^k)\ge\prod_{i=1}^\infty(1-p^{-i})\). The probability it lands anywhere else in whole or in part is at most \(\sum_{i\ge1}p^{-i}=1/(p-1)\). Consequently the entire HNF is within \(1/(p-1)+1/(q-1)\) in total variation of a public draw whose diagonal is \((N^k,1,\ldots,1)\) and whose remaining first-row entries are independent uniform residues modulo \(N^k\). This is a comparison of the two sampling laws. Every subsequent randomized operation obeys the same total-variation bound.

**Reduction and saturation.** Among index-\(p^m\) lattices the number with cyclic quotient is the number of primitive row vectors modulo \(p^m\) divided by \(\varphi(p^m)\), since unit multiples give the same kernel. This gives the exact cyclic probability

\[
C_d(p,m)=
\frac{\prod_{i=2}^{d}(1-p^{-i})}
     {\prod_{i=m+1}^{m+d-1}(1-p^{-i})}
\;\ge\;1-\frac1{p(p-1)} .
\]

If the two local quotients at \(p,q\) are cyclic, both local Smith exponent lists are \((0,\ldots,0,k)\). At every depth \(j\ge0\) the reduction of \(N^{-j}(L\cap N^j\mathbb Z^d)\) therefore has dimension \(d-1\) at both primes for \(j<k\), and \(d\) at both for \(j\ge k\). The chance of any discrepancy in this entire rank filtration is at most \(1/[p(p-1)]+1/[q(q-1)]\). In particular, dimension three typically gives the same plane rank at both primes at every depth below \(k\). This bound is uniform in \(d,k\). Writing the determinant explicitly costs \(O(k\log N)\) bits.

More generally, reduction of a uniform index-\(p^m\) lattice has codimension at least \(r\) only if it lies in some intermediate lattice with elementary quotient of order \(p^r\). Counting incidences gives, for \(1\le r\le\min(m,d)\),

\[
\Pr(\operatorname{codim}\ge r)
\le {d\brack r}_p\,\frac{A_d(p^{m-r})}{A_d(p^m)}
\le\left[\prod_{i\ge1}(1-2^{-i})\right]^{-1}p^{-r(r-1)} .
\]

Indeed, the first factor counts possible intermediate lattices and the second counts the relative number of sublattices in each. The Gaussian coefficient is at most \(p^{r(d-r)}/\prod_{i\ge1}(1-2^{-i})\), and the ratio of the two index counts is at most \(p^{-r(d-1)}\). Growing the determinant power while seeking the two middle ranks in higher odd dimension thus leaves the smaller ranks very rare.

## Checks

The diagnostic replayed every one of 1,888 Hermite forms over 36 determinant/dimension pairs through the public sampler. Direct gcd distributions verified equal per-matrix masses and the acceptance formula. In dimension three, all 4,550 lattices of index \(36=6^2\) were checked, of which 1,183 have different ranks modulo 2 and 3. For index \(5^2\), 31 of 806 lattices gave rank one modulo 5 and 775 gave rank two. Direct congruence enumeration checked every saturation depth through \(k=2\) at both primes of a cyclic sample. Five further public draws reached dimension 40 and 521-bit determinant. [Runnable checks](check_determinant_lattices.py)
