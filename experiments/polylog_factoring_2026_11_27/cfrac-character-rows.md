# What a continued-fraction relation reveals before a square dependency

On the 64-bit semiprime $N=13069950186793286509$, the continued fraction of $\sqrt N$ found the first two independent relations after 2,465 steps. They certify that **3 and $-5$ are squares at both hidden primes**. There was no matrix dependency at that point. The certificates reduce the possible factor residues modulo 840 from 192 to 48. The ensuing deterministic Harvey–Hittmeir search used 1,189 baby steps and 34 giants, compared with 2,559 and 21 for all 192 unit classes at the same modulus. Both searches recovered 3,545,884,609. The 2,465 collection steps are additional work and have a different cost from the search steps. [Extractor](cfrac_character_rows.py) · [Completion algorithm](hh_residue_lists.py) · [Exact checks](check_cfrac_character_rows.py)

## Individual labels from public congruences

Fix distinct odd base primes $\ell$ and let $m=8\prod\ell$. Associate one binary sign coordinate to each of $-1,2$ and $\ell^*=(-1)^{(\ell-1)/2}\ell$. For a prime $p\nmid m$, the coordinate is one when the corresponding Legendre symbol at $p$ is $-1$. Quadratic reciprocity gives $(\ell^*/p)=(p/\ell)$. Thus the sign word depends only on $p\bmod m$. The character map from units modulo $m$ onto sign words is surjective: modulo 8 the four units give the four combinations for $-1,2$, and each odd prime contributes equally many squares and nonsquares.

For a parity row $w$, let $D(w)$ be the signed product of the corresponding base elements. A public unit $x$ with $x^2\equiv D(w)\pmod N$ imposes $w\cdot\chi(p)=w\cdot\chi(q)=0$ over $\mathbb F_2$. At rank $h$ the admitted unit classes consequently number exactly $S=\varphi(m)/2^h$. They can be enumerated with work polynomial in the input bit length per class, plus $O(\sum\ell)$ local residues and their arithmetic. Solve the triangular row equations for the sign words, then combine local square and nonsquare classes by CRT. These classes contain both prime residues. For the example the two roots modulo $N$ are 11642448869718210396 for 3 and 6446772307084315212 for $-5$.

The Jacobi symbols of $N$ alone give the word $\delta=\chi(p)+\chi(q)$. Every proposed word $v$ for the smaller prime is compatible with the partner $v+\delta$. Even if either prime will do, a nonzero $\delta$ merely partitions all sign words into pairs: the smallest cover contains one word per pair. When $\delta=0$, every word is still possible. Square witnesses impose equations on each prime individually. [Prior study of hidden genus characters](../../notes/genusleak_results.md)

## Relation collection and its measured cost

For the $i$th convergent $a_i/b_i$ of $\sqrt N$, starting at $i=0$, the continued-fraction recurrence gives $a_i^2-Nb_i^2=(-1)^{i+1}Q_{i+1}$ with positive complete denominator $Q_{i+1}$. Keep $a_i$ modulo $N$, strip the factor-base primes and absorb any remaining integer square. Two incomplete rows with the same positive remaining cofactor also yield a complete row by multiplication. In either case the collector retains a verifying modular square root. Gaussian elimination multiplies the roots and removes the square from overlapping base coordinates. A dependency gives a square root of one and is tested immediately for a proper gcd. The traditional CFRAC algorithm uses these dependencies to factor. Here independent rows are usable *before* that stage. [Morrison and Brillhart's CFRAC method](https://doi.org/10.1090/S0025-5718-1975-0371800-5) · [Earlier repository reduction using completed square dependencies](../../notes/smoothbuild_results.md) · [Earlier smooth square relations for a specified target](short-rational-relations.md)

Each row below uses 64 seeded balanced semiprimes of exactly $b$ bits, with the same factor base $3,5,7$. The extractor has only $N$ and a budget of $b^2$ convergents. It stops at two independent rows, a direct proper gcd or exhaustion. Known factors generate the test inputs and independently verify the returned rows.

| $b$ | Ended with rank 0 | Rank 1 | Rank 2 | Direct factor | Total convergents tested |
|---:|---:|---:|---:|---:|---:|
| 32 | 0 | 0 | 26 | 38 | 6,886 |
| 48 | 12 | 9 | 19 | 24 | 85,879 |
| 64 | 48 | 12 | 2 | 2 | 252,308 |
| 80 | 64 | 0 | 0 | 0 | 409,600 |

## Conductor limit on binary characters

For *any* modulus $m$, the group of quadratic Dirichlet characters has binary dimension at most $2+\omega_{\rm odd}(m)$. Each odd prime power contributes one coordinate, and the entire power of two contributes at most two. Thus even perfect independent labels from arbitrary sources can thin the units by at most $4\,2^{\omega_{\rm odd}(m)}$. If $t$ distinct odd primes divide $m$, their product is at least $t!$, giving $t=O(\log m/\log\log m)$. For $m=N^{O(1)}$ the maximum saving in class count is $N^{o(1)}$. Also $m/\varphi(m)\le2^{\omega(m)}=N^{o(1)}$. Therefore the term $N^{1/5}(S/m)^{3/5}$ in this list algorithm remains $N^{1/5-o(1)}$ for the full set of character-compatible classes. Producing or reading a very large explicit list can cost more. These restrictions still yield exact partial information on individual inputs. The 80-bit cohort illustrates the separate difficulty of finding even a few witnesses at fixed conductor and budget. [List completion bound](hh-residue-lists.md)
