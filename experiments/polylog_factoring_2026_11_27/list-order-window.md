# Factoring from unlabelled local residue lists

Suppose one divisor $r$ of $N$ lies in $[\sqrt{N/2},\sqrt{2N}]$. For each of several pairwise coprime moduli $M_i$ coprime to $N$, a supplied list of $k_i$ residues contains $r\bmod M_i$. One lattice can recover $r$ without identifying the true entry in any list. For fixed $\epsilon>0$ and sufficiently large $N$, a sufficient condition is
$$
\prod_i M_i^{1/k_i}>N^{1/4+\epsilon}.
$$
Every list must cover the same divisor. More generally, a supplied monic polynomial $F_i$ of degree $k_i$ can specify the constraint $F_i(r)\equiv0\pmod{M_i}$. For a list, use $F_i(T)=\prod_{a\text{ in list }i}(T-a)$. The final check tests literal membership as well, since a composite modulus can have additional polynomial roots. [Decoder](list_order_window.py)

## One coefficient basis for all constraints

Write $r=c+x$ with $|x|\le X$ and a public lower bound $L\le r$. In degrees below $D$, intersect the ideal powers $(N,T+c)^m$ and $(M_i,F_i(T+c))^{t_i}$ in $\mathbb Z[T]$. At $x=r-c$ every resulting polynomial is divisible in value by $r^m\prod_i M_i^{t_i}$. For one coordinate $(M,F,k,t)$ and degree $j$, set $u=\min(t,\lfloor j/k\rfloor)$. A row of degree $j$ is given locally by
$$
M^{t-u} T^{j-ku} F(T+c)^u.
$$
Its leading coefficient $M^{t-u}$ is minimal. Indeed, expansion in powers of the monic $F(T+c)$ has remainders of degree below $k$, and the coefficient of power $v<t$ must be divisible by $M^{t-v}$. The modulus-$N$ coordinate has $F(T)=T$ and $k=1$.

To build a global row, multiply the local leading coefficients. Assemble the monic parts coefficientwise by CRT modulo the local $M^t$, and multiply by that common leading coefficient. The resulting row attains the required local divisibilities. Triangular elimination in descending degree shows that the rows span the intersection. For $m<D$ and $k_it_i\le D$, scaling column $j$ by $X^j$ gives determinant
$$
\Delta=X^{D(D-1)/2}N^{m(m+1)/2}\prod_i M_i^{k_i t_i(t_i+1)/2}.
$$
The exact sufficient test is $\sqrt D\,2^{(D-1)/4}\Delta^{1/D}<L^m\prod_i M_i^{t_i}$. LLL and Cauchy–Schwarz then make every qualifying offset an integer root of a single nonzero polynomial. The decoder verifies the row norm, isolates its integer roots and tests the resulting divisors.

At each $D$, take $m=\lfloor D\log_N L\rfloor$ clamped to $[0,D-1]$ and $t_i=\lfloor(D-1)/k_i\rfloor$. The log-radius threshold approaches $(\log L)^2/\log N+\sum_i(\log M_i)/k_i$. The rounding and LLL losses are $O(1+(\log N+\sum_i\log M_i)/D)$ uniformly in the degrees: for a full local modulus the exact exponent is $[2Dt_i-k_it_i(t_i+1)]/[D(D-1)]$, which differs from $1/k_i$ by $O(1/D)$, including $k_i\ge D$. Thus a dimension polynomial in the input size suffices for the stated fixed-power margin. For balanced divisors $L\asymp\sqrt N$ and the entire box has $X=O(\sqrt N)$.

The degree-one case is the earlier weighted order-hint calculation with full local agreement. Classical CRT soft decoding gives weighted equality tests at supplied coordinates. Coppersmith's modular-root theorem applied after CRT gives the lists-only radius $(\prod_i M_i)^{1/k}$ at common degree $k$. Here the constraint $r\mid N$ contributes a further $N^{1/4}$ in the determinant bound. This comparison assumes the lists or monic equations have been supplied. An interval written by its endpoints still costs its cardinality when expanded as a list. [Earlier weighted factor hints](factoring-with-an-order-hint.md) · [Repository CRT-decoding audit](../../notes/crtdecode_results.md) · [Guruswami–Sahai–Sudan, Theorem 4](https://people.csail.mit.edu/madhu/papers/2000/gss-conf.pdf) · [Boneh on Coppersmith and simultaneous polynomial equations](https://www.ams.org/notices/199902/boneh.pdf)

## When some lists are entirely wrong

A single public lower bound $H$ on the product of the matching moduli suffices. Their locations can be unknown. Set $M=\prod_i M_i$ and pass `min_agreement=H`. The polynomial interface recovers all divisors $r$ in the box for which $G(r)=\prod_i\gcd(M_i,F_i(r))\ge H$. The residue-list interface checks the stricter score $\prod_{i:\ r\bmod M_i\ \text{is listed}}M_i\ge H$ after extracting the roots. For composite $M_i$, a partially shared factor or a polynomial root assembled from different listed residues can increase $G(r)$.

Assign every list the same multiplicity $t$. Each term of $(M_i,F_i(T+c))^t$ at $r-c$ is divisible by $\gcd(M_i,F_i(r))^t$. Consequently the same lattice and determinant certify all qualifying offsets as soon as
$$
\sqrt D\,2^{(D-1)/4}\Delta^{1/D}<L^m H^t .
$$
When this finite inequality holds, the returned tuple is the entire set of balanced divisors in the public box that meet the chosen score. An empty tuple therefore certifies that no such divisor exists. If the search reaches the dimension cap without a certificate, the implementation raises an error and the predicate remains undecided.

For lists of common degree $k$, let $t/D$ approach $\log H/(k\log M)$. Completing the square in $t/D$ gives the sufficient range, with a fixed positive margin in units of $\log N$,
$$
\log X<\frac{(\log L)^2}{\log N}+\frac{(\log H)^2}{k\log M}.
$$
For unequal degrees put $B=\prod_i M_i^{k_i}$ and $k_{\max}=\max_i k_i$. The same calculation gives the second term $\max_{0\le u\le1/k_{\max}}(2u\log H-u^2\log B)$, which is at least $(\log H)^2/(k_{\max}\log M)$. The implementation chooses $t=\min(\lfloor D\log_B H\rfloor,\lfloor D/k_{\max}\rfloor)$ by integer comparisons and checks the finite determinant inequality. When full agreement is requested, it retains the separate degree-specific multiplicities.

The precedent is the ideal interpolation of Guruswami, Sahai and Sudan for CRT decoding with errors. Their Theorem 3 uses an ideal for each proposed residue and charges its multiplicity in the determinant. For a list of distinct residues modulo a prime, the corresponding ideals are pairwise comaximal and their intersection equals $(M_i,F_i)$. This gives the local list-recovery construction directly. The hidden divisor coordinate contributes $(\log L)^2/\log N$ to the range here. [Guruswami–Sahai–Sudan, §3.3](https://people.csail.mit.edu/madhu/papers/2000/gss-conf.pdf#page=5)

For example, take binary lists, a full balanced search interval and $\log M=(1+o(1))\log N$. Agreement on more than $1/\sqrt2$ of the log-modulus weight by a fixed margin suffices. This is a concrete test for a proposed source of hints. Once the determinant test passes, the lists and $N$ alone suffice to verify whether a hidden factor achieves the requested agreement.

Two planted trials deliberately omitted the true residue from 5 of 16 and 8 of 24 binary lists. Each modulus exceeded 256, so the public bounds $H=2^{88}$ and $H=2^{128}$ followed just from the promised number of matches. These bounds are powers of two and the matching-modulus products are odd. Only a numerical lower bound is supplied. The decoder used dimensions 7 and 9. In the second trial it recovered $1162342947127129$ from $N=1763389057895639888195926065613$. Requiring every list to match returned no divisor on either trial. An exhaustive test of all 256 four-list patterns on $29\cdot41$ compared every returned set with trial division. It certified 25 empty sets, 110 singletons and 121 sets containing both factors. [Reproducible checks](check_list_order_window.py)

## Finite check

Four planted balanced semiprimes used 4, 8, 16 and 32 binary lists. Each list held the actual residue and an independently drawn alternative, shuffled. The corresponding lattice dimensions were 7, 5, 5 and 5. On the 325-bit input the 32 lists describe 4,294,967,296 CRT combinations; the five-dimensional calculation returned the factor $5421621943611810540134878558792359926713223424479$. In every fixture, even the centered balanced box had radius greater than the lists-only $\sqrt{\prod_i M_i}$ threshold. Separate controls checked the ideal powers, exact determinants, mixed degrees and the composite-modulus membership check. [Reproducible checks](check_list_order_window.py)

## Can a local rule improve its reach by predicting a list?

Let $p,q$ be distinct uniform primes in $[R,2R)$, with $X=R^2$. For each prime $\ell\le(\log R)^A$, $A>1$ fixed, let a preselected rule return any list from $N\bmod\ell$. Give the list weight $\log\ell/k$ when it contains $p\bmod\ell$, where $k$ is its number of distinct entries. Give empty or incorrect lists zero weight. Even counting all correct coordinates with hindsight, the probability of total weight at least $t\log N$ is at most $X^{-t/A+o(1)}$ for fixed $0<t<1/2$. Single-entry tables already attain this exponent for the score. The bound also applies if either factor can be the target.

To see this, replace each nonempty list independently by one uniform member, sampling once for each possible local input. Empty lists can be filled arbitrarily. On a fixed pair the expected single-entry agreement weight is at least the original list score. Its maximum $M=\sum_{\ell\le(\log R)^A}\log\ell$ is polynomial in $\log R$. If the list score reaches $t\log N$, the probability that the rounded table reaches $(t-\eta)\log N$ is at least $\eta\log N/M$ for $0<\eta<t$, by splitting the expectation at that threshold. Apply the previously proved upper bound to each realized single-entry table, average, and let $\eta$ approach zero. The loss from rounding is polynomial in the bit length and does not affect the exponent. This argument allows lists to vary in size at every coordinate. It pinpoints the missing source for the successful planted decoder: lists using additional information about the full input. [Single-entry theorem and attaining construction](local-order-predictions.md)
