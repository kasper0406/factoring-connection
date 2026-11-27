# A square root block evaluator for mixed root parity

The shifted mixed-root parity admits a direct public evaluator using $O((1+w/\sqrt N)\operatorname{poly}(\log N))$ bit operations for an interval of width $w$. This improves the exponent of the roadmap's cubic tangent blocks for central intervals. It still takes exponential time in the bit length on those intervals, and the bound does not improve on general-purpose factoring algorithms. [Earlier tangent evaluator](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_disk_mod8/native_dyadic.md)

## Cover by certified short blocks

For a block starting at $s$, put $f(x)=(s+x)^2-1\bmod N$, written as a monic integer quadratic. Take $X=\lfloor\sqrt N/8\rfloor$ and cover the interval by disjoint blocks with at most $X+1$ integers. Handle $X\le1$ by individual tests. Use the monic shift lattice with multiplicity $m$ on any larger block. Its LLL and evaluation bound relative to $N^m$ is
$$\sqrt{2m+1}\,2^{m/2}X^m N^{-m^2/(2m+1)}.$$
Taking $m=n=N.\mathrm{bit\_length}()$ and $X\le\sqrt N/8$ makes the ratio at most $\sqrt{2m+1}\,2^{-9m/4}<1$, using $N^{m/(4m+2)}\le2^{m/4}$. Thus every modular root in the block is an integer root of a nonzero polynomial of degree at most $2m$. Integer polynomial factorization or exact real-root isolation recovers the integer candidates, and direct evaluation discards spurious roots and the known roots $\pm1$. [Lattice construction](bounded-quadratic-roots.md)

The rank is $O(n)$ and the entries have $O(n^2)$ bits, so reduction and root extraction take polynomial bit time per block with polynomial-time exact algorithms. Empty blocks are charged. The code chooses the least multiplicity passing the exact integer version of the bound. It validates the resulting auxiliary polynomial by its scaled coefficient sum. No factors of $N$ are supplied to the evaluator. Using it in the fixed-window bisection gives a complete split within the same $\sqrt N\operatorname{poly}(\log N)$ scale.

## Checks

The implementation matched direct parity on 198 seeded intervals and verified splits on three repeated-factor inputs. Four planted blocks checked the polynomial stage, including one at $N=38{,}942{,}974{,}393$ with width 24,668 and multiplicity four. Those fixture factors placed the interval and were absent from the evaluation. An end-to-end run at $N=1{,}147$ used 1,035 blocks and returned 31. These are correctness checks. The larger fixture exercises one block, and no large central-window speedup was measured.

Run `python lattice_window.py 11009 --start 100` or `python lattice_window.py 1147 --split` with SymPy available.
