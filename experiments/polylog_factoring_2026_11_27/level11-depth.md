# What a full quaternion frame reveals about one prime component

For a public projective line modulo a balanced semiprime, the level-11 ideal class was already a computable predictor candidate. Its two possible values can have exactly the same frequencies on every row and column. We now keep the change-of-basis matrix modulo a selectable power of two. There is a transition at modulus of order $N^{1/6}$ in how much fixing a prime component changes the distribution of this public output. Exact equality of two outputs can also be analyzed across the entire range of moduli. [Earlier class-bit experiment](level11-sampling.md) · [Frame construction and level-44 coefficient](level11-marking.md) · [Equal-frame gcds](level11-collisions.md)

## A public experiment with a hidden conditioning event

Take $N=pq$ for distinct odd primes other than 11, and fix the public splitting of the quaternion order modulo $N$. A uniform line $L$ modulo $N$ consists of independent uniform lines $L_p,L_q$ by the Chinese remainder theorem. There are $p+1$ and $q+1$ possibilities, respectively. Fix $M=2^k$ with $k\ge2$. The public function of $L$ is its ideal class together with its transport matrix over $\mathbb Z/M$, up to the units of the reference ideal. It takes expected polynomial time in $\log N+k$ to set up and evaluate. All probability statements below hold for **each** fixed splitting. [Implementation](level11_marking.py)

Write $\nu$ for the law of the output, $\nu_i$ for its conditional law when $L_p=i$, and $u$ for the uniform law on its determinant sector. In this sector the determinant on the principal class is $N$ modulo $M$, and on the second class it is $N/3$ modulo $M$. The sector contains

$$S=|\mathrm{SL}_2(\mathbb Z/M)|(1/4+1/6)=5M^3/16$$

states. Indeed the determinant fibers have $3M^3/4$ matrices each, by reduction modulo 2 and the three-dimensional lift at each step. The two reference unit groups have orders four and six and act freely. Their reductions are faithful already modulo 4. For example, a nonidentity finite-order matrix congruent to the identity modulo 4 would have a nonidentity power of prime order in the same subgroup. If $v_2(A-I)=h\ge2$, then squaring increases this valuation by exactly one, and taking an odd power retains it. Thus all marked states carry equal mass in every sector.

## Bounds above and below the transition

Put $\rho_r=2\sqrt r/(r+1)$. The normalized Hecke transition for a prime $r$ sends sector uniformity to sector uniformity and contracts centered $L^2$ norms by at most $\rho_r$. For this statement use the supersingular isogeny graph in characteristic 11 with full level $M$. Its determinant sectors are the Weil-pairing fibers. The adjacency operator is normal, and its restriction to functions whose sum is zero on each fiber has spectral radius at most $2\sqrt r$. Extending a centered density by zero from one fiber and restricting to the destination gives the same bound for a single step, also for the adjoint. The common factor converting between the full and single-fiber $L^2$ norms cancels. This uses the known spectral theorem. [Codogni and Lido, §§2.2–2.3 and Theorem 2.3.8](https://doi.org/10.1016/j.jnt.2026.02.006)

For fixed $i$, the $q+1$ choices in the row run through the $q$-neighbors of the intermediate marked ideal of norm $p$. The full law starts at one vertex and takes a uniform $p$-step followed by a uniform $q$-step. These assertions follow by localizing the line ideal at $p$ and $q$. The two prime choices are independent, and changing the local basis permutes the neighbors. A point mass has centered density norm $\sqrt{S-1}$ in its sector. Therefore

$$\lVert\nu_i/u-1\rVert_{2,u}\le\sqrt{S-1}\rho_q,\qquad
\lVert\nu/u-1\rVert_{2,u}\le\sqrt{S-1}\rho_p\rho_q.$$

In particular the total variation distance between **every** row and the full law is at most $\tfrac12\sqrt{S-1}\rho_q(1+\rho_p)$. Swapping $p,q$ gives the column bound. Thus when $M^3=o(q)$ and $q$ grows, every row becomes indistinguishable from the full law in total variation.

In the other direction let $A_i$ be the support of row $i$, and set $a_i=|A_i|/S$. A row has at most $q+1$ outputs, so $a_i\le(q+1)/S$. Applying Cauchy–Schwarz to the centered indicator of $A_i$ gives

$$\nu(A_i)\le a_i+\rho_p\rho_q\sqrt{(S-1)a_i(1-a_i)}
\le \frac{q+1}{S}+\frac4{\sqrt p}.$$

Since $\nu_i(A_i)=1$, the distance from **every** row to the full law tends to one when $M^3/q\to\infty$ and $p\to\infty$. This conclusion still holds when the full law itself occupies a small fraction of the sector. For balanced inputs $p\le q\le2p$ the transition scale in both directions is $M\asymp N^{1/6}$. The bounds specify the two sides of that scale and leave the constant-scale profile open.

For example, set $M$ to the first power of two at or above $N^{1/4}$. Each row and column then has distance at least $1-O(N^{-1/4})$ from the full law, whereas the full law is within $O(N^{-1/8})$ of sector uniformity. There is also an information bound. Writing $X$ for the public output, its Shannon entropy in bits is at least its collision entropy, so

$$I(X;L_p)\ge \log_2\frac{S}{q+1}-\log_2\big(1+(S-1)\rho_p^2\rho_q^2\big)
\ge \tfrac14\log_2 N-O(1).$$

Here $H(X\mid L_p)\le\log_2(q+1)$. The same statement holds for $L_q$. For a **fixed** $N$ these bits concern a random line reduced at a fixed prime. Estimating a hidden numerical prime from the observations is a different computational task.

## Which information would actually give a split?

The support test attaining the lower bound uses a row indexed by the hidden prime. A supplied intermediate mark of norm $p$ at $M\asymp N^{1/4}$ would already recover a balanced factor. For the principal starting ideal and reference norm $d\in\{1,3\}$, undo the known matrix splitting and multiply by $d$. The resulting four residues are the coordinates of an integral quaternion $\beta$ with norm $dp$. Positive definiteness bounds its coordinates by $O(N^{1/4})$. Enumerate their constant number of lifts, divide each norm by $d$ when integral, and test for a proper divisor of $N$. This applies to either prime. Given the factors one can also compute the appropriate intermediate mark directly, using the reduction of the same splitting. The evaluator actually implemented on $N$ produces the final mark of norm $N$. [Lift control](level11_frame_lift.py) · [Executable checks](check_level11_frame_lift.py)

Another public test draws two independent lines, retains distinct lines with equal marks and takes the gcd of their coordinate determinant with $N$. This gcd is a proper divisor precisely when the lines agree in exactly one prime component. The collision bound covers this selection rule for every $M=2^k$, $k\ge2$. In particular, equality of marks for distinct lines **precludes** a shared component as soon as $M>2\sqrt{\max(p,q)}$. The large statistical distance above therefore supplies no equality-based shortcut at that precision. [Proof and bounds at lower precision](level11-collisions.md)

## All 600 lines at $N=19\cdot29$

The third column measures how uniform the output is over the entire sector. The fourth column is the *least* full-law mass excluded by any one fixed-19-line row. The last two columns count ordered pairs of distinct lines.

| $M$ | States observed / sector size | Full TV to uniform | Minimum mass excluded by a row | Equal-mark pairs | Pairs giving a factor |
|---:|---:|---:|---:|---:|---:|
| 4 | 20 / 20 | 1/20 | 16/75 | 17,640 | 1,200 |
| 8 | 156 / 160 | 13/75 | 11/15 | 2,072 | 96 |
| 16 | 494 / 1,280 | 393/640 | 553/600 | 228 | 0 |
| 32 | 578 / 10,240 | 4831/5120 | 113/120 | 44 | 0 |

The diagnostic uses the known factors solely to construct the grid, group the public labels and verify the results. It checks the one-step and two-step collision bounds, the support inequality and the injectivity criterion using exact arithmetic. Run `python experiments/polylog_factoring_2026_11_27/check_level11_depth.py`. [Diagnostic](check_level11_depth.py)
