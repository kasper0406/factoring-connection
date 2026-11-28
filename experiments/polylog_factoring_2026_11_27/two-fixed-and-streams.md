# Two fixed multipliers give polynomial-size lists on typical inputs

Let $P,Q$ be independent uniform odd residues modulo $2^k$, and fix $h=3$ or $h=5$. Observe the low product $n=PQ\bmod 2^k$ and the low $k$ bits of the two words

$$
 P\mathbin{\&}Q,\qquad (hP)\mathbin{\&}(hQ).
$$

Let $F_k$ count all compatible ordered pairs of odd residues. A single word $P\mathbin{\&}Q$ leaves $2^{(0.18872\ldots+o(1))k}$ pairs on typical inputs. With the second word, **the distributions of $F_k$ are tight**. For any prescribed error probability there is a finite list bound independent of $k$. Moreover, the ordinary exact bit-by-bit enumerator visits only polynomially many nodes with probability tending to one. Both claims still hold if any fixed number $r$ of bottom bit positions is omitted from either or both words. The probability here averages over the two inputs, including their product. [Single-stream result](shared-and-fiber-entropy.md)

More precisely, for $r\le t\le k/2$ and $t\ge1$ there is an absolute constant $K_h$ with

$$
 \Pr(F_k>2\cdot4^t)\le K_h(9/10)^t. \tag{1}
$$

This bound concerns exact lists for individual observations, with the probability of each observation proportional to its number of preimages. It gives a polynomial runtime on a density-one set. The worst-case list and the mean runtime have separate quantifiers.

## The finite carry calculation

In binary multiplication by $h$, a bit $a$ and an incoming carry $c\in\{0,\ldots,h-1\}$ give output $(c+ha)\bmod2$ and next carry $\lfloor(c+ha)/2\rfloor$. To compare the two observed words for two ordered input pairs, use one carry for each of their four coordinates. Form an $h^4$ by $h^4$ nonnegative integer matrix $T$. From carry state $(c_1,c_2,c_3,c_4)$ try each $(a_1,a_2,a_3,a_4)\in\{0,1\}^4$. Put $u_i=c_i+ha_i$. Increment the entry for the four next carries when

$$
 a_1a_2=a_3a_4,\qquad
 (u_1\bmod2)(u_2\bmod2)=(u_3\bmod2)(u_4\bmod2).
$$

Consequently $(T^b\mathbf1)_s$ counts all compatible quadruples of length-$b$ words starting at state $s$. Integer arithmetic gives the following certificates. In the last three columns $v=T^7\mathbf1$.

| $h$ | Short row-sum certificate | $\min v$ | $\max v$ | $\max_s(5(Tv)_s-36v_s)$ |
|---:|---|---:|---:|---:|
| 3 | $\max T^3\mathbf1=508<8^3$ | 584,368 | 1,319,652 | $-16,424$ |
| 5 | $\max T^5\mathbf1=28{,}744<8^5$ | 410,192 | 1,390,424 | $-724,264$ |

The short certificate already gives growth strictly below $8^b$ by grouping consecutive blocks. The vector certificate gives a common explicit bound. Positivity of $T$ as an order-preserving map and $Tv\le(36/5)v$ imply, for every $b\ge0$,

$$
 \max_s(T^b\mathbf1)_s\le
 \frac{\max v}{\min v}(36/5)^b. \tag{2}
$$

The checker constructs every transition from the two equalities above. [Exact matrix and direct finite checks](check_two_fixed_and_streams.py)

## Where an alternative can first differ

Fix $1\le d<k$ above any omitted positions and set $b=\min(d,k-d)$. Consider a potential alternative whose first coordinate agrees with $P$ in its low $d$ bits and flips bit $d$. The equality of the odd products forces its second coordinate to agree with $Q$ through the same $d$ bits and flip the same bit. Indeed multiplication by an odd unit preserves the 2-adic valuation of the difference in the inverse-partner map.

Call the common low representatives $x=P\bmod2^d$, $y=Q\bmod2^d$. Write the following $b$ bits of the actual pair as $u,v$ and of the alternative as $u',v'$. Then the necessary product condition modulo $2^{d+b}$ becomes

$$
 v'-v\equiv c(u'-u)\pmod{2^b},
 \qquad c=-y/x\bmod2^b .
$$

The term of order $2^{2d}$ disappears since $b\le d$. As $x,y$ range over independent uniform odd inputs, $c$ is uniform over the $2^{b-1}$ odd residues, independently of the actual high words $u,v$. The odd differences $u'-u,v'-v$ prescribe exactly one $c$. The incoming multiplier carries depend on $x,y$ and can correlate with $c$. For an upper bound, allow all $h^2$ common incoming carry pairs for the two paths. For each pair of carries, (2) bounds how many high quadruples can meet both observed streams. Dropping restrictions on that count only enlarges it.

Each specified $u,v$ has probability $4^{-b}$, and the required slope has probability $2^{-(b-1)}$. Summing these probabilities bounds the expected number of matching alternative **blocks** with the prescribed first difference by

$$
 h^2\max_s(T^b\mathbf1)_s\,4^{-b}\,2^{-(b-1)}
 \le 2h^2\frac{\max v}{\min v}(9/10)^b. \tag{3}
$$

Every full alternative with first difference $d$ requires one of these blocks. By (3), the probability of any alternative whose first difference lies in $[t,k-t]$ is at most $C_h(9/10)^t$. Indeed each value of $\min(d,k-d)$ occurs at most twice and the resulting geometric series converges.

For a fixed observation, call a leaf good if no other leaf first differs from it in $[t,k-t]$. All good leaves sharing a low-$t$ first coordinate must agree through bit $k-t$ inclusive. Their remaining top bits give at most $2^t$ choices, and there are at most $2^t$ low words. Thus at most $4^t$ leaves are good. Conditional on the observation, the actual ordered pair is uniform on its fiber. Whenever $F_k>2\cdot4^t$, the actual leaf has conditional probability greater than $1/2$ of having an alternative in the middle. This proves (1). The argument works with the bottom $r$ positions omitted because (3) used only the fully available bits at and above $d\ge r$.

## Exponential worst-case fibers for every fixed odd menu

For any finite collection of fixed positive odd multipliers, choose $k=2s$ with all multipliers less than $2^s$. Take $x=1+2^s a$ and $y=1-2^s a$ modulo $2^{2s}$ for every odd $a\bmod2^s$. Their product is always $1\bmod2^{2s}$. For any multiplier $h$ in the collection, the low $s$ bits of $hx,hy$ both equal $h$. Their upper words are $ha,-ha\bmod2^s$. For odd $z$, the words $z$ and $-z\bmod2^s$ have bitwise AND equal to one, since all their bits above the common least set bit are complementary. Consequently every observed stream equals $h+2^s$ on the entire family. The same low product and menu values therefore admit at least $2^{s-1}$ ordered pairs. This is a worst-case statement about modular prefixes. The density-one bound and the planted example concern different input distributions.

## Exact search and the binomial observations

At depth $j$, try both lifts of the smaller-factor candidate residue, compute the odd inverse partner using $n\bmod 2^{j+1}$, and compare the next available bits of both ANDs. This visits exactly the compatible prefixes. Equation (1) gives a terminal-list tail $O_h(S^{-\eta})$ for list cutoff $S$, with $\eta=\log(10/9)/\log4$. The estimate extends to cutoffs beyond $2^j$ by triviality and to the remaining constant range of cutoffs by increasing the implicit constant for fixed $r$. Applying it separately to the at most $k$ depth marginals shows, for example, that every list has at most $k^{14}$ elements with probability tending to one, because $14\eta>1$. This deliberately loose exponent follows directly from the integer certificate. It does not estimate typical practical list sizes.

For two distinct primes in $[R,2R]$ and $2^k\le R^{1-\delta}$, with fixed $\delta>0$, Brun–Titchmarsh and the prime number theorem bound the probability of each odd residue pair by a constant times its uniform-model probability. Thus the two list estimates transfer to this prime sample. Ordering the primes preserves the observations and list size. One can take $k=\lfloor(\log_2N)/4\rfloor+O(1)$, since only a bounded number of such depths occur in the dyadic sample. At that depth each residue takes polynomial time to complete by the balanced affine divisor algorithm. [Prime-counting transfer](shared-and-fiber-entropy.md#enumeration-and-the-balanced-prime-model) · [Affine recovery](coupled-rank-windows.md)

For $N=pq$ with distinct odd primes $p<q<2p$, the required bits have explicit central-binomial indices. To read position $j$, put $D=2^{j+1}$ and $t=\lfloor hN/D\rfloor$. If $D\ge4h$ and $D/2<p$, the answer to $N\mid\binom{2t}{t}$ is exactly bit $j$ of $(hp)\mathbin{\&}(hq)$. In fact $2t<p^2$ since $q<2p$. For a prime $s\mid N$ put $v=N/s$, and take the odd representative $a=hv\bmod D$. The units digit of $t$ in base $s$ is $\lfloor as/D\rfloor$. Doubling this digit carries exactly when $a>D/2$, since its distance from $s/2$ before flooring is at least $s/D>1/2$. The inequality $2t<p^2$ eliminates higher base-$s$ carries. This is exactly the valuation criterion from the earlier adaptive binomial reduction. At quarter depth the upper size condition holds for sufficiently large balanced inputs. For $h=3$ begin at $j=3$, for $h=5$ at $j=4$, and for $h=1$ at $j=1$. Omitting those bottom bits is covered by (1). This gives two fixed, nonadaptive streams with polynomial postprocessing on a density-one set of inputs. Computing their Boolean values from $N$ remains open. [Earlier varying-multiplier construction](adaptive-binomial-lift.md)

As a finite check, independently grouping complete residue pairs through nine bits verified 212,080 fiber/scale cases for $h=3$ and 214,178 for $h=5$. At nine bits each multiplier had a **maximum** ordered list of 32 with both full words, and 40 when the bottom two positions of both words were omitted. A separate low-slope/high-word calculation counted exactly the same 5,466 and 5,458 first-difference pairs as direct grouping on the tested shorter blocks. These finite counts and the asymptotic tail concern different scales. [Reproduce all checks](check_two_fixed_and_streams.py)
