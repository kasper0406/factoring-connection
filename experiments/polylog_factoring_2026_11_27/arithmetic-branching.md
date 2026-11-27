# Equality branches in arithmetic samplers

Adaptive modular arithmetic with equality tests can correlate the two hidden CRT residues. For this model, that correlation has a quantitative cost. Insert a gcd probe at each equality test. If the probes seldom split the input, the sampler follows one transcript with high probability, and its output conditioned on that transcript has independent local residues.

## Model

Fix $N=pq$ with distinct odd primes. Inputs are independent uniform residues, or any inputs whose complete local random tapes are independent at $p$ and $q$. On a fixed transcript, all operations act separately on the two CRT components. Examples include modular arithmetic and public exponentiation. A branch tests whether a residue is zero modulo $N$. Equality of two residues uses their difference. A division can use the same probe for its denominator and continue with inversion only on units. Take a finite execution tree, and charge a gcd per probe.

Instrument each test of $w$ by computing $\gcd(w,N)$. A proper gcd returns a divisor. The continuing branches correspond to both local values being zero or both being nonzero. At leaves, probe the output in the same way, treating a doubly zero output as a failed decode with dummy output 1. The remaining outputs are units.

## Concentration of the transcript

Condition on reaching a node without a split. Every preceding continuing branch constrained the two local tapes by a product event. Consequently their conditional distributions remain independent. If the local probabilities of zero at this test are $a,b$, the two continuing branches have probabilities $ab$ and $(1-a)(1-b)$. The split probability is $a(1-b)+(1-a)b$. The arithmetic-geometric mean inequality gives
$$a(1-b)+(1-a)b\ge 2\sqrt{ab(1-a)(1-b)}\ge2\min\{ab,(1-a)(1-b)\}.$$

For the proof, follow the larger continuing branch at each node. If $r$ is the probability that any probe returns a proper divisor, the total probability of departing this path into a smaller continuing branch is at most $r/2$. Indeed, sum the node inequalities weighted by the probability of reaching each node on the chosen path. The events that a probe splits at those nodes are disjoint. The resulting terminal transcript has probability $m\ge1-3r/2$, and its unit output has a product distribution. The path is used only in the proof. No computation of branch probabilities is required. When the lower bound on $m$ is positive, a continuing path always exists through the chosen nodes.

## Implication for square decoding

Let $s$ be the square-decoder success probability for the original sampler. If $r\ge s/3$, the probes already split with probability at least $s/3$. Otherwise the original successful mass on the chosen continuing path is at least $s-3r/2>s/2$, so the conditional decoder success probability there is at least $s/(2m)$. The [independent sampling bound](independent-sampling.md) supplies two cases.

* If the pairwise gcd bound holds on that path, two instrumented runs both reach its leaf with probability $m^2$. Taking the gcd of the difference of their outputs thus splits with unconditional probability at least $m^2(s/(2m))^2/4=s^2/16$.
* If the mode case holds, a unit output with a successful square decode occurs with unconditional probability above $m/4\ge1/8$.

The same public routine works across the cases. Probe the equality tests in two runs, decode each available output, and try the gcd of their difference. This gives a test for proposed ways of correlating the local residues. Representative-dependent branching such as floors, order comparisons, bit extraction and sampling from a short integer interval lies outside the model. The reduction does not bound the success probability of such proposals.
