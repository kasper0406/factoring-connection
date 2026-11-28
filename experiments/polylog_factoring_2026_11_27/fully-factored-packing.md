# A prime conductor and a factored congruence modulus

For every fixed mixed input \(N,x\), one can construct infinitely many packed sphere instances for which **both the squared radius and the congruence modulus come with their complete prime factorizations**. The conductor in the squared radius is a single prime. An unconditional theorem proves termination for each fixed input. Its constants depend on \(N\), so a polynomial bound for this simultaneous-prime search remains open. Allowing composite coordinates with known prime factors gives a different construction with a uniform estimate. [Factored coordinate sampler](factored-coordinate-packing.md) The preceding two-prime construction supplies a factored squared radius in expected polynomial bit time for short \(x\). [Two-prime construction](factored-gaussian-packing.md) · [Prime-conductor generator](fully_factored_packing.py)

## The construction works on either side of the diagonal

As before, \(N\) is a product of two distinct odd primes, \(A=x^2+1\), and \((\frac{-A}{N})=-1\). Choose distinct odd primes \(B,q\) coprime to \(2AN\) with \(B>2N\), and test that \(T=B^2+4N^2q^2\) is prime and coprime to \(A\). Take \(t=2q,\ C=Nt\). Then
\[
 Q=2ABt=4ABq,\qquad P=B^2-C^2-2xBC,\qquad M=A T^2.
\]
The prime factorizations of \(Q\) and \(M\) are now explicit after trial division of the short \(A\). The queried modulus is \(Q\ell\), where the test prime \(\ell\) is public and prime to \(Q\). The numbers \(P,Q\) are coprime. At primes dividing \(2Bt\) this follows from \(\gcd(B,C)=1\) and the opposite parity of \(B,C\). At a prime dividing \(A\), the congruence \(P\equiv(B-xC)^2\) shows that a common factor would divide \(T\).

A negative \(P\) causes no difficulty. For a potential \(W\equiv Pu\bmod Q\) of squared norm \(AT^2\), set \(v=(Pu-W)/Q\) and \(D=u\cdot v+xN\). Comparing norms gives \(ABt\mid D\), and the identity \(Q D=A(B^2-C^2)-u\cdot W\) gives
\[
 -\frac{N^2t}{B}\le D\le\frac{B}{t}.
\]
The bounds lie strictly between \(-ABt\) and \(ABt\) whenever \(AB^2>N^2\) and \(At^2>1\), as here. Hence \(D=0\) and \(v\cdot v=N^2\). The same 72 marked frames and factor-revealing residue predicates apply. Giving the complete prime factorization of the radius and modulus takes no hidden inputs once \(B,q,T\) are certified. [Bijection, primitivity, and CRT readout](fixed-field-packing.md)

## Fixed-input existence and the size question

Take \(H=1,2,4,\ldots\). The implementation draws \(B,q\) uniformly among the primes in \([2NH,4NH]\) and \([16H,32H]\), respectively. A successful first box has \(q\le31\) and \(T<(16+4\cdot31^2)N^2\), an absolute constant times \(N^2\). If it fails, the search widens both intervals. In the first box, at least one of \(19,23,31\) is coprime to \(N\), and all three are coprime to \(A=x^2+1\).

Green and Sawhney's Theorem 1.1 counts primes \(B,q,B^2+nq^2\) for each fixed \(n\equiv0,4\bmod6\), with any smooth compactly supported weight on the two input coordinates. Here \(n=4N^2\) satisfies the congruence. With \(N\) fixed, take a nonnegative smooth weight supported inside \((2N,4N)\times(16,32)\) with positive integral. At scale \(H\) the asymptotic gives \(\gg_N H^2/(\log H)^3\) actual prime pairs in the box for all sufficiently large \(H\). Removing prime-power terms is included in the paper's derivation of the box count. For large \(H\), the excluded prime factors of \(2AN\) cannot occur. The prime number theorem in each marginal interval gives a conditional success probability \(\gg_N 1/\log H\) per uniformly drawn prime pair. With \(8\,\operatorname{bitlength}(2NH)^2\) independent trials at each height, the randomized search therefore terminates almost surely. Deterministic exhaustive search over successive boxes also terminates. [Green–Sawhney, Theorem 1.1 and its remarks](https://arxiv.org/pdf/2410.04189#page=1)

The Euler product in their equation (1.1) gives
\[
 \frac{\kappa_{4N^2}}{\kappa_4}
 =\prod_{\substack{r\mid N\\r\equiv1\bmod4}}\frac{r-1}{r-3}\ \ge1 .
\]
The obstacle to a uniform running-time estimate is their error \(O_{W,n}(H^2(\log\log H)^2/\log H)\) and the height from which the main term dominates. The associated number field stays \(\mathbb Q(i)\). The ray modulus in their §8.4 still contains \(2n\). A later generalization also permits arbitrary dependence on the fixed linear forms and the choice of field generator. Neither stated estimate gives a bit-length bound for varying \(N\). [Green–Sawhney, §8.4](https://arxiv.org/pdf/2410.04189#page=42) · [Achenjang–Woo, Theorem A and §2](https://arxiv.org/html/2611.16062v1)

One testable sufficient hypothesis is an inverse-polynomial lower bound in \(\log N\) for the fraction of successful prime pairs in the **first** box for each sufficiently large balanced mixed input. That would make the complete-factorization construction at size \(T=O(N^2)\) run in expected polynomial bit time with an exact polynomial-time primality test. The implementation uses deterministic Miller–Rabin in its proved range and Lucas certificates after factoring the candidate minus one above that range. This certificate generation itself has no polynomial worst-case guarantee.

## Finite checks

All 32 seeded balanced inputs at 16, 24, 32 and 40 bits certified a prime conductor in the first box. Their median was 10 sampled prime pairs and their maximum was 62. A separate complete enumeration at \(N=251\cdot257=64{,}507,\ x=0\) found 4,157 prime conductors among all 52,995 first-box prime pairs. The 64- and 96-bit fixtures also passed in the first box. The latter produced a certified 203-bit conductor and a 406-bit squared radius.

| \(N\) | \(x\) | \(B\) | \(q\) | Prime \(T\) | \(Q\) | Scalar queries | Recovered factor |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 15 | 0 | 37 | 23 | 477,469 | 3,404 | 120 | 3 |
| 65 | 4 | 211 | 19 | 6,145,421 | 272,612 | 528 | 5 |
| 77 | 1 | 163 | 17 | 6,880,493 | 22,168 | 238 | 7 |

Each row used direct integer-square tests on the requested scalar sphere residues and checked both supplied factorizations by multiplication and primality certification. These results concern the cost of constructing the instances. The cost of answering a specified residue query on growing inputs remains open. [Reproduce all checks](check_fully_factored_packing.py)
