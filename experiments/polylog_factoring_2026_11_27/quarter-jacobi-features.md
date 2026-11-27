# Polynomial Jacobi features for the scaled quarter bit

A predictor built from a small number of Jacobi symbols of low-degree polynomials still has exponentially small advantage on balanced semiprimes. This extends the single Jacobi feature test. The feature polynomials and the classifier may even depend on the factors. The result concerns the hidden-pair bit and fixed feature lists. [Earlier baseline](quarter-periodic-features.md) · [Weak prediction target](quarter-fourier.md)

Let $N=pq$ with distinct odd primes $p<q\le2p$. On each of at most $S$ integer intervals, let the predictor select its expected sign from the scale's residue modulo a period $d\le B<p$ and the $s$ features $(P_i(c)/N)$. The feature polynomials are fixed for the modulus and have degree at most $D$. If a polynomial reduces identically to zero at either prime, its feature is constant zero and can be incorporated in the labels. For the remaining polynomials, the probability that any feature is zero on a uniform unit is $O(sD/p)$ by the finite-field root bound.

On each arithmetic progression cell, expand the label on signs in the Boolean Fourier basis. The sum of absolute values of the coefficients is at most $2^{s/2}$, since their squared sum is at most one. Extend the expansion at zero by averaging the label over the missing signs. This changes the predictor only on the exceptional scales already counted. Each term in the expansion is the Jacobi symbol of the product $P$ of a subset of the feature polynomials, of degree at most $sD$.

There are two cases for such a product.

* If it is a constant times a square over both prime fields, the Jacobi feature is constant outside its zero set. Replace it there by the constant and charge the same root bound for the change. The remaining correlation with the quarter bit is bounded by the periodic-feature discrepancy.
* Otherwise at least one prime-field reduction has a factor of odd multiplicity. For that prime $r$, the Weil bound for a quadratic multiplicative character with an additive twist gives
  $$\left|\sum_{c\bmod r}(P(c)/r)e^{2\pi i kc/r}\right|=O((1+sD)\sqrt r)$$
  uniformly over $k$. One can remove square factors first and account for their zeros. At the other prime use the trivial bound. By CRT, every normalized additive Fourier coefficient of the Jacobi feature is therefore $O((1+sD)/\sqrt p)$. The indicator of the progression cell and the hidden quarter bit each have Fourier absolute sum $O(\log N)$ by the geometric-sum bound. Fourier convolution gives correlation $O((1+sD)\log^2 N/\sqrt p)$ for the normalized sum on the cell.

Summing the feature expansions over at most $SB$ cells and converting from uniform residues to uniform units yields the advantage bound
$$O\left(\frac{SB\,2^{s/2}(1+sD)\log^2 N}{N^{1/4}}+\frac{sD}{\sqrt N}\right).$$
In particular polynomial $S,B,D$ in the bit length and $s=o(\log N)$ give exponentially small advantage. The bound does not cover arbitrary high-degree maps, many adaptive tests, or predictors using other global information.

The diagnostic computes the optimal in-sample label for the interval quartile, the scale modulo 3, and four fixed polynomial Jacobi symbols. On moduli 1,147, 11,009 and 50,429 its exact advantages were approximately 0.231, 0.064 and 0.032. These are factor-supplied, finite controls rather than evidence for the asymptotic claim. Run `python check_quarter_jacobi_features.py`.
