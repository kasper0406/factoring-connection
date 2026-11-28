# Binary multiplication carry parity and public factor tests

Let $N=pq$ for odd primes $p<q\le 2p$, let $b$ be its bit length, and write
$U_a(T)$ for the binary digit polynomial of $a$. A polynomial-size menu of
proposed carry parity words gives a deterministic factor test. A fixed
number of bit errors per word can be corrected by search. A proposed word
can even split $N$ when it differs substantially from the true word.

## Which information suffices

Write $s_j=[T^j]U_pU_q$ over the integers, $n_j=[T^j]U_N$, and $c_j$ for
the carry entering column $j$. Schoolbook multiplication gives
$s_j+c_j=n_j+2c_{j+1}$ and $c_0=c_b=0$. Modulo two this becomes

$$
U_p U_q=U_N+H,\qquad H(T)=\sum_{j=1}^{b-1}(c_j\bmod 2)T^j
\quad\text{in }\mathbf F_2[T].
$$

Enumerate all words $E$ within Hamming distance $k$ of a proposed word,
retaining the possible constant bit $E_0=0$. Factor $F=U_N+E$ over
$\mathbf F_2$ by squarefree decomposition and Berlekamp's algorithm.
The lower prime lies in $[\lceil\sqrt{N/2}\rceil,\lfloor\sqrt N\rfloor]$.
For each monic divisor of $F$ of an eligible degree, read its coefficients
as an integer and test that integer against the interval and $N$ by division.
This verifies every reported split. [Berlekamp, 1967](https://doi.org/10.1002/j.1538-7305.1967.tb03174.x)

Single irreducibles can miss the answer. For example,
$1+T^2=(1+T)^2$ in $\mathbf F_2[T]$. The coefficient strings encode the
prime 5 and the integer 3 respectively, whereas $3^2=9$.
Binary evaluation therefore requires recombining irreducibles in the
polynomial ring. A Hensel lift to higher powers of two would require a
target polynomial specified to that precision.

If $F=\prod_i f_i^{e_i}$ and $d_i=\deg f_i$, the number of eligible divisors
is the sum over eligible degrees of the coefficients of
$\prod_i(1+Y^{d_i}+\cdots+Y^{e_i d_i})$. A dynamic program computes the
counts for each suffix of the factor list. Skip $F$ when its eligible count
exceeds a public budget $B$. Otherwise, enumerate only exponent choices
with a feasible suffix. For $L$ templates there are at most
$L\sum_{t=0}^k\binom bt$ polynomials, each requiring polynomial time to count
and at most $B\,\mathrm{poly}(b)$ time to enumerate. With fixed $k$,
polynomial $B,L$ and a polynomial-time list generator, the complete
procedure runs in polynomial time on every input.

The precise success criterion on a non-skipped word is $U_p\mid H+E$, by
subtracting $U_pU_q=U_N+H$. Matching $H$ itself suffices. A different word
can have the same residue modulo $U_p$, in which case the complementary
polynomial need not be $U_q$.

## Cost of the cap on random prime pairs

For a uniformly random monic degree $d\ge1$ binary polynomial with
constant term one, the expected number of monic divisors is exactly
$(d+3)/2$. The constant and whole-polynomial divisors contribute one each.
For $1\le j<d$ there are $2^{j-1}$ possible degree-$j$ divisors with
constant term one and $2^{d-j-1}$ possible complements. Dividing by the
$2^{d-1}$ possible inputs gives $1/2$ per intermediate degree.

For uniform distinct primes in $[R,2R]$, put $b\asymp\log R$. One prime
digit polynomial has expected divisor count $O(b^2)$. Indeed, the interval
meets at most two whole binary lengths. The preceding identity bounds the
sum of counts over its odd integers by $O(Rb)$, and the prime number
theorem supplies $\gg R/b$ primes. For two independent primes the product
of their counts consequently has mean $O(b^4)$. Conditioning on
distinctness costs a bounded factor. Ordering the pair preserves this
symmetric count. Divisor counts satisfy
$\tau(U_pU_q)\le\tau(U_p)\tau(U_q)$, including shared factors.

Take $B=b^{4+A}$ for any fixed $A>0$. By Markov's inequality the
probability that the *true* corrected polynomial exceeds the budget is
$O(b^{-A})$. If an $N$-dependent menu covers the true word within radius
$k$ with probability $\delta$, the public algorithm splits with probability
at least $\delta-O(b^{-A})$. The subtraction holds for randomized menus
and for adversarial dependence between coverage and divisor counts. It
bounds an unconditional loss over prime pairs. Incorrect templates are
subject to the same public cap.

## Prevalence of a fixed parity word

On two independent uniformly random odd words, reveal successive pairs of
factor bits. At a free position $j\ge1$, all terms except the two endpoint
bits in column $j$ are already known. For a determined integer $a$,

$$
c_{j+1}=\left\lfloor(a+p_j+q_j)/2\right\rfloor.
$$

The sum of the fresh bits is 0, 1 or 2 with probabilities $1/4,1/2,1/4$.
For each $a\bmod4$, either specified parity for $c_{j+1}$ is missed with
conditional probability at least $1/4$. By successive coupling with fresh
uniform random variables, the number of misses in $m$ free columns
stochastically dominates a binomial variable with parameters $m,1/4$.
Taking all odd words below the least power of two above $2R$ gives
$m=\log_2 R+O(1)$. Conditioning both samples to be distinct primes in
$[R,2R]$ costs $O((\log R)^2)$. For fixed $k$, the chance of a true parity
word within $k$ errors of a fixed template is at most

$$
O((\log R)^2)\sum_{t=0}^k\binom mt4^{-t}(3/4)^{m-t}
=R^{-\log_2(4/3)+o(1)} .
$$

A union bound preserves the exponent for polynomially many fixed templates.
Templates computed from the complete $N$ can correlate with unrevealed
factor bits. The estimate also leaves open hits on a different word
having the same residue modulo $U_p$.

A template may even read the whole low half of the product. Take independent
uniform odd $p,q<2^{m+1}$, and reveal $u=pq\bmod 2^{m+1}$. Conditional on
$u$, the odd word $p$ is uniform and $q=u/p\bmod 2^{m+1}$. At step $j$,
$1\le j\le m$, the low $j$ factor bits and $u$ determine $p_j\oplus q_j$
by the column equation. When this XOR is zero, the fresh fair bit $p_j$
chooses between endpoint sums 0 and 2. The next carry parity is then fair
even given all previous factor bits and $u$. Unconditionally the number of
these ambiguous positions is $\mathrm{Bin}(m,1/2)$ because the original
factor bits were independent. For any cutoff $K$ and any word chosen from
$u$, the chance of at most $k$ misses is therefore at most

$$
\Pr\{\mathrm{Bin}(m,1/2)<K\}+\Pr\{\mathrm{Bin}(K,1/2)\le k\}.
$$

For the second term, take the first $K$ ambiguous positions. Each is
identified before its fresh fair bit is revealed. If fewer than $K$ occur,
pad the sequence of errors with independent fair coins. The resulting $K$
errors are independent fair bits, without conditioning on the number of
ambiguous positions. Taking $K=\lfloor m/4\rfloor$ and fixed $k$ gives
$2^{-(1-H_2(1/4))m+o(m)}$, where $H_2$ is binary entropy. A menu of
polynomial size and the restriction to balanced prime pairs preserve the
exponent. A further $s$ bits of arbitrary information multiply this bound
by at most $2^s$ through enumeration of the possible messages. Again, this
bounds agreement with the true word.

## Reproducible checks and prior work

Run `python check_carry_parity.py`. Independent long division checks all
255 odd monic polynomials of degrees 1 through 8, including repeated
factors. Across 1,024 pairs of odd words with five free bits, the exact
histogram of odd carries in the five tested columns is
$147,211,254,222,127,63$. All 78 unordered prime pairs in $[64,128)$
split when given their measured words. The public experiment always starts
from the zero template and uses only $N$ during decoding.

| Error radius | True words covered out of 78 | Verified public splits out of 78 |
|---:|---:|---:|
| 0 | 0 | 1 |
| 1 | 1 | 16 |
| 2 | 9 | 60 |

For $851927=881\cdot967$, only incoming carries at positions 6 and 19 have
odd value. Twelve other nonzero carries have value two. The public
two-error run tests 113 masks and 62 eligible divisors before returning
881. This instance was selected by its scored carry pattern.
[Decoder](carry_parity.py) · [Exact checks](check_carry_parity.py)

A separate seeded experiment drew 1,024 unplanted products at each prime
bit length. In order of radii 0, 1 and 2, the uncapped split counts were
$(2,17,114)$ at 12 bits, $(0,2,15)$ at 16 bits, $(0,0,2)$ at 20 bits,
and $(0,0,0)$ at each of 24, 32, 48 and 64 bits. True-word coverage was
$(0,2,4)$ at 12 bits and zero at all larger sizes. Using just $N$, the
capped decoder reproduced all 131 inputs that the offline residue test
identified as positive at any radius. A separate unfiltered replay on 19
inputs gave no splits. For instance, $1870182893=35449\cdot52757$ splits at
radius one with true carry weight 17. A separate bound on these actual
zero-template splits includes the extra residue matches.
[Full proof and cohort table](carry-syndrome-bounds.md) ·
[Seeded check](check_carry_syndromes.py)

Multiplication carry polynomials occur explicitly in Devevey, Sakzad,
Stehlé and Steinfeld's work on integer polynomial learning with errors.
The polynomial reformulation of US 12,034,849 also studies carry corrections
in integer factorization, via transformations of quadratic equations.
Binary convolution height studies the overlap counts $s_j$. The earlier
repository experiment used modular-reduction carries as public features
for a quarter-interval sign. Here the unknowns are the parities of the
hidden multiplication carries, and the public test enumerates eligible
binary-polynomial divisors. [Devevey et al., Lemma 2](https://eprint.iacr.org/2021/277) ·
[Patent, polynomial reformulation](https://patents.justia.com/patent/12034849) ·
[Binary convolution height](https://doi.org/10.1080/27684830.2026.2727395) ·
[Earlier carry-feature screen](carry-feature-screen.md)
