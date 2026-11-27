# Modulated window averages and the cost of locating a band

Take the same two-component input and cubic or quartic phase as in the [long-window evaluator](coarse-phase-window.md). For a supplied center $c\bmod M$, $M=N^2$, the evaluator now supports the mean of $e_M(-ct)$ times the normalized magnitude or squared magnitude over the unit scales in a window of length $L$. The center can be anywhere on the frequency circle. At $L\ge N/\operatorname{poly}(\log N)$ and inverse-polynomial additive error the total cost is polynomial in the input bit length.

## Public evaluator at any supplied center

Keep the same finite local Laurent list and mass bound as before. Each mixed line has different local labels $a,b$ and satisfies $(f-a)(f-b)=0\pmod M$. The geometric kernel now depends on the centered distance $|f-c|_M$. Use the old omission bound with this distance, and search the cyclic band of radius $B=\min(\lfloor M/2\rfloor,\lceil WrM/L\rceil)$ where $r$ is the inverse error parameter. The unit-density correction is unchanged because the multiplier has modulus one. The short-window branch still searches the whole circle.

For each pair of labels, put $h=(a+b)/2\pmod M$. The target is $(f-h)^2=h^2-ab\pmod M$. The interval $[c-B-h,c+B-h]$ can be searched as an ordinary integer interval, using periodicity to cover a cyclic wrap. The existing monic quadratic lattice uses disjoint blocks of at most $1+\lfloor N/8\rfloor$ integers at modulus $M$. It certifies all modular roots independently of the starting point. Each verified root with $1<\gcd(f-a,N)<N$ returns a component after the trial-division cutoff. If all bands are empty of mixed lines, return the public diagonal sums with the modulation included. [Certified square-root blocks](lattice-window.md)

There are $O(1+B/N)$ blocks per pair of labels. For bounded precision costs this removes the quadratic dependence on $N/L$ of the quotient search. More precisely, squared magnitudes take $O(1+rN/L)$ block work and unsquared magnitudes take $O(r^2(1+rN/L))$. Both have a polynomial factor in $\log N$ and the bit lengths of the window and precision parameters. For zero-centered windows, the code can still use quotient enumeration when its cost is within a polynomial factor of block search. The trigonometric kernels need the complex geometric sum after modulation. The checks cover a negative center, wraparound, two planted mixed bands, and an empty band.

## Coverage of simple centers

A uniform integer center hits a mixed line with probability at most $(2B+1)/M$. At inverse-polynomial error and $L\ge N/\operatorname{poly}(\log N)$, even a polynomial number of uniform centers and local labels gives total hit probability at most $\operatorname{poly}(\log N)/N$ for any fixed two-component input.

There is a stronger statement for sweeping small rational centers, even if their numerators and denominators are chosen from the full input. Let $p,q$ be distinct primes in $[X,5X/4]$, and fix unequal labels $a,b$ with magnitude polynomial in $\log X$. Suppose $f\equiv a\pmod{p^2}$ and $f\equiv b\pmod{q^2}$ lies within $NT$ of some $jM/k$, where $T\ge1$ and $1\le k\le K$. Choose the lift of $f$ closest to this center and set
$$U=\frac{kf-jM-ka}{p^2},\qquad V=\frac{kf-jM-kb}{q^2}.$$
Both are integers, $|U|,|V|=O(kT)$, and $Up^2-Vq^2=C=k(b-a)$. For $K$ polynomial in $\log X$ and large $X$, $0<|C|<X^2$. This forces $U,V$ to be nonzero with the same sign, even when $T$ grows faster than a polynomial in the bit length.

Fix $U,V$. The ratio $p/q$ differs from $\sqrt{V/U}$ by $O(|C|/(|U|X^2))$. Distinct prime pairs give distinct reduced fractions spaced by at least $1/(5X/4)^2$. There are therefore $O(1+|C|/|U|)$ such ordered pairs. For each $|U|=A$, the ratio bound allows only $O(A)$ choices of $|V|$. Summing gives $O(k^2T(T+|a-b|))$ pairs for fixed $k$, uniformly over all numerators $j$. Summing $k\le K$ gives
$$ O\big(K^3T(T+|a-b|)\big) $$
exceptions among $\Theta(X^2/\log^2 X)$ ordered balanced prime pairs. A polynomial-size local label set and polynomial tolerances retain exponentially small coverage in the bit length. Rounding the center or displacing it by $O(N\operatorname{poly}(\log N))$ only enlarges $T$ by a polynomial amount. With polynomially many labels and denominators, inverse-polynomial coverage needs $T$ of order at least $X/\operatorname{poly}(\log X)$. At that radius the block search visits a power of $N$ blocks when certifying an empty band.

The exceptional inputs have a separate factor search. Enumerate the small coefficients $U,V$ and test the discriminant of $Uz^2-Cz-VN^2=0$ at $z=p^2$. Integer square and divisor checks verify the resulting factors. The recovered components also let us check the band condition directly. This costs polynomial bit time for the specified ranges and needs no supplied numerator.

For a fixed real scaled center $\gamma$ there is also a weaker bound without a denominator restriction. Fix $q$ and write $e=p^2(p^{-2}\bmod q^2)$. Then the mixed residue satisfies $(f-a)/M=(b-a)(p^{-2}\bmod q^2)/q^2\pmod1$. The residues $p^{-2}\bmod q^2$ are distinct as $p$ ranges over the prime interval, and multiplication by the nonzero $b-a$ is invertible. An arc of width $O(T/N)$ contains $O(T+1)$ of these grid positions, including the $O(|a|/M)$ correction. A polynomial list of such scaled centers independent of $p$ thus has coverage at most $\operatorname{poly}(\log X)/X$ on this family.

The remaining problem is a center-selection rule that uses the arithmetic of $N$ beyond these families and hits useful frequencies with inverse-polynomial probability. If no band contains a mixed line, the explicit public diagonal is already a legal answer at the requested accuracy.

[Modulated evaluator and tests](coarse_phase_window.py)
