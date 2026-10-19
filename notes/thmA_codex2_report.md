Independent referee report: algebraic and algorithmic half of Theorem A
19 October 2026

I find no fatal algebraic or algorithmic obstruction to the final design. Several statements do not apply literally as written, but I give elementary repairs below that preserve the claimed exponent. In particular, the final hidden-root argument works, the randomized multiplier does not break the number-field congruence, and the square-root computation has acceptable bit complexity.

Scope: I read thmA_proof.md first, followed by the specified mathematical portions of lemmaM_results.md, hypD_results.md, pcest_results.md, lv6_results.md, and ldbound_results.md. I did not consult the excluded review sections, candidate cards, audit files, or earlier referee reports. My verdict concerns the requested algebraic/algorithmic audit and its interfaces. The first- and second-moment arguments are examined here for the inputs and outputs consumed by the algorithm, rather than presented as a separate, exhaustive analytic referee report. Acceptance below should be read with that division of responsibility.

Locations below refer to files in /Volumes/git/factoring-connection/notes/. Line numbers identify the version read in this audit.

1. Theorem DD and the exact fruitfulness probability — CORRECT.

Location: lemmaM_results.md §1, lines 82–111; thmA_proof.md §3 and §4 Step 9.

For odd N = product p_j^{e_j}, each unit square has exactly 2^omega(N) roots. Conditional on R_i, write r_i = r_i^0 epsilon_i, with independent uniform epsilon_i in mu_2(Z/N). If the exposed transcript and the nonempty dependency S are invariant under changing the epsilon_i, fixing every epsilon except one indexed by S leaves X_S/Y uniform in mu_2. Exactly two of these elements give trivial gcds: the globally constant signs +1 and -1. For every other sign pattern,

  gcd(X_S - Y, N) = product of p_j^{e_j} over the +1 components.

Thus the conditional probability at a reached final gcd is exactly 1 - 2^{1-omega(N)}, including for nonsquarefree N. This is a statement about that final gcd; the complete algorithm also has earlier factor exits, so its overall success probability need not equal this number.

The final sampling range [1,N-1] is essential to exactness. The earlier [2,N-2] range mentioned at lemmaM line 93 does not have the claimed exact conditional probability. For example, at N=15 and R=1 it retains roots 4 and 11 and excludes 1 and 14; both retained roots yield proper factors against Y=1. The final pseudocode has repaired this issue.

2. Transcript invariance for the whole final algorithm — CORRECT.

Location: thmA_proof.md §3, lines 121–173, and §4 Step 9; lv6_results.md §3.7.

I checked the following dependencies explicitly:

* Field construction, including g0, g1, w and c, uses independent field coins. Keeping a field after a successful descent biases the field and R distribution, but does not distinguish roots in a fixed square fibre.
* Gauss reduction, ties resolved independently of r, depends on R. Complete factorization of the resulting coordinates, the distinct-prime shape decision, and LP's random coins are invariant on that fibre.
* SQ_SEARCH, its candidate ordering, batch smoothness, and cached atoms depend on (A,ell). A cache shared between rows does not expose their hidden signs.
* Rational and ideal columns, characters, kernel vectors, the Hensel prime, local square-root branches, exact square tests, the choice among s1,s2,s1+s2, and Y depend only on certificates and fresh coins.
* Gcd checks on field and certificate data are invariant. The check on r itself has invariant unit/nonunit status; a nonunit draw terminates with a factor.
* The roots are first used non-invariantly when X_S is formed for the final gcd. After a trivial gcd the algorithm discards all rows. Old gcd outcomes can bias the retained field, but do not bias new independent signs inside new square fibres.

Consequently there is no hidden-root leak in the specified algorithm. The retained field need not be independent of the first accepted R; fibre invariance, which is what DD actually uses, survives that conditioning.

3. Theorems S and S' — CORRECT.

Location: hypD_results.md §2, lines 161–213; application at thmA_proof.md §4 Step 7.

Let h2 = 2^omega(N). Character orthogonality gives at least |S|^2/h2 ordered pairs with square ratio. For pairs in [Y/2,Y]^2, equality of ratios modulo N implies equality over Q because |uv'-u'v|<N. The total ratio multiplicity is bounded by

  2 sum_{g<=Y} g floor(Y/g)^2 <= 2Y^2(1+log Y).

Cauchy–Schwarz and |Q_N|=phi(N)/h2 give precisely the lower bound in S.

For k-products, classify a pair by the number i of common primes. Its primitive ratio can have at most binomial(n,i) partners, where n=|G|. Summing the resulting bounds gives O_k(|S|^2), since sum_i binomial(k,i)^2 = binomial(2k,k). For n>=4k^2 the binomial comparisons used in S' hold uniformly. Hence

  Pr(square ratio represented) >= |S|^2/(C_k h2 phi(N)).

The stronger condition 2Y^2<N is needed when turning existence into the Euclidean shortest-vector test. If a target vector has coordinates bounded by Y, its norm is at most sqrt(2)Y; the shortest vector has no larger norm. A nonzero determinant between them would be a multiple of N but have absolute value at most 2Y^2. Thus they are proportional. In the final design 2 floor(sqrt(N/2))^2<N because N is odd, so the boundary causes no problem.

These results do not depend on the field model, F(g0,g1), or the candidate's prime-b restriction.

4. Theorem K and the final descendable-prime set — CORRECT WITH FIXES.

Location: hypD_results.md §3, lines 221–268; pcest_results.md §1, lines 186–215; thmA_proof.md §4 Steps 2 and 7.

K treats G_A abstractly. Its counting argument transfers to free M=N^mu and F(g0,g1)=Nw once G_A is the set recognized by the actual search. The prime-b and |G| filters have to be included in that definition; they cannot be added after proving a moment bound for a different set. The final definitions do include them.

At most k distinct primes means squarefree products, including the empty product 1. This rule is divisor-closed. When the primitive part of a represented k-product ratio is taken, common primes disappear and at most k primes remain on each side. Requiring exactly k would invalidate this implication. Coordinate signs must be handled separately from membership of the absolute values.

If g_A>=gamma/2 and gamma=L^{-o(1)}, then |G_A|>=4k^2 for all sufficiently large N: pi_J=N^{2nu+o(1)} dominates the subexponential loss. S' then gives

  q_A >= C_k g_A^{2k} (log N)^{-O(k)} 2^{-omega(N)}

for the mathematical row event, with early factor exits counted as successful termination where necessary.

The old displayed K(ii) and its proof still use the weaker fourth-power version of S. The final roadmap uses the valid S' strengthening to the second power. The standalone theorem should be updated instead of being cited as if its display already said this.

There is also a small error in the old shape-count justification at hypD line 267: all accepted absolute coordinates belong to the UNION of products of j distinct primes, 0<=j<=k, rather than necessarily the exactly-k set S_J. Use

  T_J = union_{j=0}^k {products of j distinct primes of J}.

Then |T_J|<=C_k pi_J^k. At most 2|T_J|^2 signed ratios can have the right shape, and division by |Q_N| supplies the needed shape upper bound. Combining it with the preceding success lower bound proves the asserted polylogarithmic times g_A^{-2k} times 2^{2omega} ratio on good fields. This ratio is L^{o(1)} on good fields; it need not be so on every field that can be retained. The field-selection averaging, not a pointwise good-field assertion for every retained field, supplies the global expected bound.

5. Retaining the first successful field: harmonic mean and factor exits — CORRECT WITH FIXES.

Location: hypD_results.md Corollary K, lines 277–283; lv6_results.md §4, lines 440–452; ldbound_results.md §1, lines 152–167; thmA_proof.md §4 Steps 6–7.

Without competing factor exits, if q_A is the probability that one attempt returns a row and m_row=E q_A, the retained-field distribution is q_A dP(A)/m_row. Therefore

  E[1/q_A | retained] = Pr(q_A>0)/m_row <= 1/m_row.

This identity is correct. No lower bound on every retained field's q_A is required. For C further rows, a uniform bound C_trial on the expected cost of an attempt gives conditional cost at most C C_trial/q_A, even when the cost of an attempt correlates with its outcome.

The treatment of factor exits in the assembled proof needs more care. From E[g_A 1_{not rejected}]>=gamma/2 one cannot, without another argument, infer the same bound after conditioning on fields that neither reject nor factor N. Factor-producing fields can carry some of that mass. Likewise, a special-q search may find a factor on an earlier nonsmooth candidate instead of returning its later counted atom. Thus the lower bound for row probability alone does not follow verbatim from the displayed gamma/2 inequality.

Here is a complete repair that requires no new distribution theorem. For one raw field trial, let delta be the probability of returning a factor, b=E q_A the probability of delivering a row, and p=delta+b the total probability of useful termination. A mathematical successful descent either supplies a row or is interrupted by a successful factor exit. Applying the good-field counting argument with these two outcomes combined gives p>=L^{-o(1)}. Conditional on reaching a retained row-field, its law is still proportional to q_A. It is reached with probability b/p. Consequently its contribution to the expected cost of C subsequent rows is at most

  (b/p) (C C_trial/b) = C C_trial/p.

The initial selection costs at most C_trial/p, including its field-sampling cost. Fresh-row rounds have final-gcd success probability at least 1/2, so they add at most a factor 2. Additional factor exits during collection only shorten the run. This proves the required unconditional expected bound without claiming that the nonterminating-field law preserves gamma.

This is a local probabilistic bookkeeping defect with an explicit elementary repair, not an unresolved relation-supply gap.

6. Lenstra–Pomerance shape test — CORRECT.

Location: thmA_proof.md §3 line 151 and §4 Step 7, lines 262–263.

The external result really is an expected-time algorithm for COMPLETE prime factorization of every input integer. It is not an average-over-inputs theorem, a smoothness recognizer, or a theorem conditional on the factor sizes. The main theorem on the first page states this explicitly. See Lenstra–Pomerance, JAMS 5 (1992), “A rigorous time bound for factoring integers”:
https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf

A shortest vector in the determinant-N Gauss lattice has coordinate sizes O(sqrt N). Uniformly for inputs n<=sqrt(2N), LP's bound becomes

  exp((1+o(1)) sqrt(log n log log n))
    <= L_N[1/2,1/sqrt2+o(1)].

Complete factorization permits exact checking of the prime interval and squarefreeness. A cap on LP is unnecessary for an expected-time analysis. Its random coins do not leak r because its integer inputs depend on R alone.

7. Base expansion, final sizes, and coprimality — CORRECT.

Location: pcest_results.md §1 Lemmas GB and 1.2; ldbound_results.md §1, lines 171–201; thmA_proof.md §4 Step 1.

Applying the base expansion to Nw gives the exact identity Fhat(g0,g1)=Nw. Its coefficient estimates give

  |Fhat(P)| <<_d (M + Nw/M^d) V^d.

Since w=N^{mu+nu+eta0/2+o(1)} and h>=1-d mu+nu+eta0,

  Nw/M^{d+1} <= N^{h-eta0/2+o(1)}.

Thus the offset is smaller than the random-value scale MHV^d by a fixed polynomial margin. This is the margin actually needed in the first-moment sign/window argument; equality F(g0,g1)=N is unnecessary.

For a rational-smooth candidate, G(P)=+-ell n with n B-smooth. The prime b satisfies B<b<ell, so b divides neither ell nor n. From primitivity of g it follows that (g1,G(P))=1. Preprocessing gives (N,G(P))=1 when ell does not divide N. Also w>2MVmax>=|G(P)|, so (w,G(P))=1. Hence

  Fhat(P) == N b^d g1^{-d} w (mod |G(P)|)

is indeed a reduced residue with coefficient independent of w. The rational-smoothness hypothesis is indispensable here; the statement does not hold for arbitrary candidates before that filter.

The sizes consumed downstream are X_r=N^{mu-nu}L^{O(1)} and X_F=N^{mu+h+d nu}L^{O(1)}. Their exponent sum is the stated theta.

8. Moment-to-algorithm interface — CORRECT WITH FIXES.

Location: pcest_results.md §§2–4.1; ldbound_results.md §§2–5; thmA_proof.md §4 Steps 3–5 and §6 U5.

The final FM hypotheses have fixed margins:

  q0-2nu = mu-nu >= eta0;
  xi/2-q0 = [h-mu+(d-2)nu]/2 >= eta0/2;
  H/V^2 = N^{h-2nu}L^{-O(1)} >= N^{eta0}L^{-O(1)}.

The affine offset estimate in item 7 is uniform in w. Goodness is quantified uniformly over reduced residues, so fixing w and then averaging it is legitimate. The published LV Lemma 8.5 has Q<=sqrt(Psi(x,y)), effective constants, and a restricted smooth-modulus error with a density factor. It is not ordinary BV over all moduli. Those distinctions are respected in pcest's decomposition. Harper's large-conductor term includes both Q and x^{1/2-eta} log^2 x; the fixed level margin absorbs both. Sources checked: LV Lemma 8.5 and Fact 8.11, and Harper Proposition 2:
https://arxiv.org/pdf/1805.08873
https://arxiv.org/pdf/1208.5992

The kernel transfer also has the required final margin. Its exact index is |det(P,P')|, and on the prime-b candidate set the last adjacent minor is det(P,P') (bb')^{d-2}. Consequently its fibre bound is O_d(1/(H^2 V^{2d-4})), without an unquantified compactness constant.

In LD*, averaging over w is done after replacing the translated windows by a common window. The resulting pair-count function then depends on g through n=G(P)/ell. This order is essential to applying divisor switching to the partner along its g-fibre. The sizes give W/t=N^{eta0/2}L^{-O(1)}, so Brun–Titchmarsh costs only a fixed constant depending on eta0. The large-sieve family argument's most restrictive polynomial comparison is

  N^{(5-d)nu-h} L^{O(1)} = o(1),

which follows from h>=2nu+eta0 for every d>=3. Its kernel and fibre factors multiply to V^4/Y^2 times rho_r^2 rho_F^2, namely S^4 rho_r^2 rho_F^2. I found no extra power of L in this interface.

Two precisions are needed:

* FM's proved conclusion is averaged over ell as well as A. thmA Step 3 writes E_A N_ell as though the lower bound held for each ell. The bad-modulus argument proves only the averaged assertion. Use E_{A,ell}; this is exactly what PZ and gamma consume.
* A computable S with the stated leading L-exponent need not make E N>=1 or gamma a reciprocal polylogarithm. Put T=S^2 rho_r rho_F=L^{o(1)}. FM and LD* give gamma>=T/[(log N)^C(1+T)]>=L^{-o(1)}, which suffices for fixed k. U5 uses this correct logic, but U4's stronger blanket gamma>=1/(log N)^{O(1)} must be weakened for that choice of S.

I have not independently reconstructed all of LV's small-conductor analysis or Harper's proof. These are identified external inputs, not newly established results of this report. I found no failed range hypothesis at this interface.

9. Signed row identity and the parity column — CORRECT in the final roadmap; CORRECT WITH FIXES for standalone LV6.

Location: lv6_results.md (1.1), lines 100–106, and §3.1; thmA_proof.md §3, lines 176–180 and §6 U6.

With t=f_d, beta=t alpha and z_P=t a-b beta,

  psi(z_P)=t G(P)/g1,
  ell == sigma_P g1 t^{-1} psi(z_P)/n_P (mod N).

Thus a signed primitive split gives

  r_i^2 == q_i t^{-E_i} product_P psi(z_P)^{e_iP},
  q_i=sign(u_i v_i) product_P (sigma_P g1/n_P)^{e_iP}.

The split sign cannot be omitted. For example, when R=-1 is a unit square, the shortest vector can be (1,-1); its empty atomic product has q=-1, not q=1. The final roadmap includes this sign. The old LV6 formula assumes positive u,v and requires the same explicit correction.

The E-parity column makes E even and handles both the un-factored g1 contribution to q and t^{-E}. Small-prime parity columns then make q a positive rational square. There is no need to factor the large-prime part of g1, nor to add a w-parity column.

10. Lemmas C and H, and the Nw transfer — CORRECT.

Location: lv6_results.md §§3.1–3.2; ldbound_results.md §6.1, lines 521–532.

Away from t disc(f), the local order is maximal and unramified. A primitive linear atom a-b alpha can be divisible only at the degree-one prime corresponding to a/b, with valuation v_r(F(a,b)); if r divides b, it is a unit. Since F(a,b) is B-smooth, good primes larger than B contribute no valuations. Thus the non-character kernel maps into the group of square classes with even valuations outside primes over D=t disc(f).

The exact-sequence argument gives

  |H| <= 2^{d omega(D)} |Cl(K)[2]| 2^{r1+r2}.

The elementary ideal-count fallback bounds h_K by sqrt(|Delta_K|) times a fixed power of log|Delta_K|, up to a d-dependent constant. Since ||f||=N^{O(1)} and d is fixed, log|H|=O_d(log N). Computing the maximal order, class group, or factorization of D is unnecessary.

The map psi on Z[beta] remains valid because

  g(t g0/g1)=t^{d-1}g1^{-d}Nw == 0 (mod N).

If w divides t, disc(f), or the order index, it is simply among the bad rational primes supported in H. After gcd(w,N)=1, this does not obstruct reduction modulo N. No unjustified coprimality with w is being used.

11. Lemma E: at most one exceptional class — CORRECT, with the external-input distinction below.

Location: lv6_results.md §3.3, lines 215–248.

The relative quadratic discriminant has support only over 2D. At odd primes its exponent is at most 1; at primes over 2 it is at most 2 ord_P(2)+1. Therefore

  |Delta_{K(sqrt h)}| <= |disc g|^2 8^d |D|^d = Delta_bar.

Stark's required zero interval for degree 2d is satisfied by nu1<=1/[4(2d)! log Delta_bar]. Since d is odd, K has no quadratic subfield, and a degree-2d extension has at most one quadratic subfield: two would generate a degree-four subfield, while 4 does not divide 2d. If its quadratic subfield is F, the extension must be KF. Distinct square classes therefore give distinct quadratic F. The discriminant tower gives |Delta_F|^d dividing Delta_{K(sqrt h)}. The chosen Landau interval consequently forbids two exceptional classes, including a common zero counted in two distinct L-functions.

I verified the exact Stark and Landau formulations in LV Facts 6.22 and 6.24, pp. 36 of the arXiv pagination. The original Stark scan was not retrievable through the available access, and I did not independently inspect Davenport's original presentation. This is a limitation of original-source verification, not a failure of a hypothesis: the explicit versions reproduced in the primary LV research paper suffice for the application above.
https://arxiv.org/pdf/1805.08873

12. Lemma X and the character sampler — CORRECT WITH FIXES.

Location: lv6_results.md §3.3, lines 256–299; thmA_proof.md §3 lines 167–168.

The sampler's root-count acceptance probability |Z_r|/d gives every pair (r,s) the same probability before conditioning on output. It therefore samples uniformly over degree-one prime ideals in the interval. Here r is much larger than every atomic norm, so the Legendre symbols evaluated on atoms are nonzero.

For a nonexceptional quadratic extension, Lagarias–Odlyzko applies with n_L=2d and log(x/2)>=20d(log Delta_bar)^2. Its ordinary error is o(x/log x). Any real-zero term has 1-beta>=nu1>>_d 1/log N, and log x>> (log N)^2 makes that term o(x/log x) as well. Removing prime ideals of degree at least two costs O_d(sqrt x). Dividing the surviving inert degree-one count by the upper bound d pi(x) yields the conservative separation probability 1/(5d). A union bound over |H| with the prescribed t_chi leaves a subgroup of order at most two.

The primary LO Theorem 1.3 really has the threshold exp(10 n_L (log Delta_L)^2) and the stated exceptional-zero term:
https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf

Two presentation fixes: specify, for example, the least integer x satisfying the displayed lower bound. An arbitrary enormous x also satisfies it and need not have polylogarithmic bit length. With that choice, log x=O_d((log N)^2), and sampling/root-finding is polynomial in log N. Also, absence of a quadratic subfield excludes a zero in Stark's narrower interval, not necessarily all of LO's wider interval. A possible zero outside the narrower interval is bounded by exactly the decay calculation just given. The sampler's claimed cost still follows.

13. Kernel sampling, three dependencies, and the accepted-root hypothesis — CORRECT WITH FIXES.

Location: lv6_results.md §2 lines 118–121 and §§3.2–3.4, especially lines 321–327; thmA_proof.md §4 Steps 7–8.

The matrix has at least twelve dimensions of left kernel. For two independent uniform samples in a binary kernel of dimension r, the probability of linear independence is (1-2^{-r})(1-2^{1-r}). If the surviving square-class group has order at most two, linearity of the square-class map implies that one of s1,s2,s1+s2 maps to zero; all three vectors are nonzero after the independence check. This does not assume that the rows themselves are independent as vectors.

Kaltofen–Saunders Theorems 2–4 do provide the required black-box singular-system solver and uniform solution-manifold construction with essentially quadratic field-operation cost. Padding the transpose and using fast Toeplitz products preserves that bound:
https://kaltofen.math.ncsu.edu/bibliography/91/KaSa91.pdf

For the binary sampling construction, explicitly take the coefficient of the basis element 1 when projecting x-y from F_{2^e}. Conditional on A^T x, binary x is uniform on an affine kernel coset, while the solver sees only that right-hand side and its own coins. Subtracting the fixed 1-coordinate of y therefore gives a uniform binary kernel vector. The phrase “one fixed coordinate” should identify this coordinate; an arbitrary other coordinate has no x contribution.

A separate literal hypothesis mismatch occurs at lv6 §2 line 120. Accepted roots are not uniform over the entire unit group: accepting a row filters their squares. What is true, and what the proof needs, is independence and uniformity inside the square fibres conditional on the certificates. With the retained field fixed, the first retained row and subsequent successful rows have the same accepted-row law. Replace the overly strong unconditional-uniformity hypothesis by this fibre statement. DD then applies exactly as in items 1–2. The algebraic part of LV6 does not need an unconditioned root distribution at all.

14. Lemma SQ, the Hensel prime, and the final congruence — CORRECT.

Location: lv6_results.md §§3.5–3.6; thmA_proof.md §3 lines 170–172.

For Gamma=g'(beta)^2 product z_P^{|e_P|}, a square root delta in K satisfies delta/g'(beta) integral because its square is integral. Trace duality gives

  O_K subset O_K^vee subset Z[beta]^vee = g'(beta)^{-1} Z[beta],

so delta lies in Z[beta]. This justifies exact coefficient recovery in the polynomial order; it does not assume that order maximal.

The norm of Gamma has only primes <=B or dividing D. The prime p in (B,2B] avoiding D is therefore an odd unramified unit prime for Gamma. There are far more primes in that interval than the O(log N) possible prime divisors of D, and a sieve finds p in L^{a+o(1)} time. Factoring g mod p and enumerating the 2^{tau-1} sign choices is harmless because tau<=d is fixed. A true square root must reduce to one of those choices up to a global sign and lifts uniquely; exact verification rejects false positives.

The inverse-Vandermonde bound is valid: each derivative at a conjugate is bounded above by g_d, their product has absolute value |disc g|>=1, and therefore each derivative is bounded below by g_d^{-(d-1)}. This proves the stated coefficient bound A and

  log A = O_d((1+sum |e_P|) log N) = O_{d,k}(M_rows log N).

Product-tree construction and Newton lifting thus have L^{a+o(1)} bit cost, up to fixed powers of log N, rather than an unaccounted N^{o(1)} cost.

Finally, for gamma_minus=product_{e_P<0} z_P^{|e_P|}, the proposed Y satisfies

  Y^2 = q_S t^{-E} product psi(z_P)^{|e_P|}/psi(gamma_minus)^2
      = q_S t^{-E} product psi(z_P)^{e_P}
      = product_{i in S} r_i^2 (mod N).

The denominators are units: g1 and t were checked; G(P)=+-ell n is coprime to N; and psi(g'(beta))=t^{d-2}f'(g0/g1) was checked. This reasoning is unchanged for Nw and for repeated prime factors of N.

15. Field rejection and Lemma IR — CORRECT WITH FIXES to the conditioning inference only.

Location: lv6_results.md §4; ldbound_results.md §1 “Order of operations” and §6.1; thmA_proof.md §4 Step 6.

The irreducibility certificate argument is uniform in the final g and w. For every small sieve prime r, r does not divide Nw. If r does not divide g1, the affine family is exactly the degree-at-most-d polynomials with a prescribed nonzero value at g0/g1. If r divides g1, its leading coefficient is fixed nonzero and all lower coefficients vary freely. Counting full-degree elements in each nonzero norm fibre of F_{r^d} gives an irreducible fraction at least 1/(4d), uniformly in that prescribed value.

For Q the product of the sieve primes up to (log H)/2, Q=H^{1/2+o(1)}. The coefficient-cube distribution modulo Q is at most e times uniform. CRT and the irreducible fractions bound the probability of NO certificate, not just reducibility, by

  exp(-log H/(13d log log H)).

Thus an irreducible polynomial lacking the certificate is correctly included in the rare rejection event. The test costs a polynomial in log N.

After gcd(g1,N)=1, the leading coefficient is nonzero modulo g1 because Nw is a unit there. Also |f_d|<<MH<N when theta<1. Consequently a nonunit f_d yields a proper factor, not a degree-drop rejection. Content e divides Nw and satisfies e<=|f_d|<N. If it does not give a factor of N, e is 1 or w. In the latter case the field has zero counted algebraic-smooth values because w>B. The elementary upper bound Pr(e=w)<=1/w+1/(2H+1) is sufficient regardless.

Fixing all coefficients except c_{d-1}, the derivative condition N|f'(m) fixes at most one value in [-H,H], since g1 is a unit and 2H+1<N. All true rejection probabilities are therefore L^{-omega(1)}. Hence, under the RAW product law,

  E[g_A 1_{not rejected}] >= gamma-Pr(rejected) >= gamma/2.

This inequality is valid. Its subsequent use after deleting factor exits requires the repair in item 5. Rejecting before applying the moment estimates to the raw law would indeed be invalid; the final order of operations avoids that mistake.

16. Prop BS and exact recognition of the counted candidates — CORRECT WITH FIXES to the exit wording.

Location: pcest_results.md Lemma 1.3; lv6_results.md §6.2; thmA_proof.md SQ_SEARCH and §4 Step 2.

An irregular lattice has no prime-b primitive candidate: noncentral lines have distance greater than 2V from the origin, and the central line contains only two primitive points, both too short to have b>V/2. In a regular lattice, reduced-basis counting bounds the whole box by O(S^2). This argument is independent of M, w and the field degree.

Let P_B be the product of primes <=B. For a nonzero positive input x and 2^e>=log_2 x,

  gcd(x, P_B^{2^e})

is its full B-smooth part. Repeated squaring after reducing P_B modulo x gives this gcd without forming the enormous power. Product and remainder trees share the initial reductions among all O(S^2) values. With each value having O(log N) bits, the total cost is soft-O(B+S^2 log N). Trial division of the chosen two smooth values costs soft-O(B log N), already included at this scale. Zero values must be excluded first, and 1 handled as smooth; the final search does so.

Subject to no earlier proper-factor exit, the deterministic search succeeds exactly when N_ell(A)>=1. The strict phrase “succeeds iff” must include the successful-termination alternative when an earlier nonsmooth candidate yields a factor. A counted candidate itself has both gcds equal to 1, by rational smoothness, preprocessing, and ell not dividing N. There is no loss of a counted row through an ordinary rejection. Factor exits are handled as in item 5.

17. Cost assembly and optimization — CORRECT WITH FIXES to the field-selection line.

Location: thmA_proof.md §4 Step 10; lv6_results.md §7; ldbound_results.md §6.

The row target is B L^{o(1)}, with O(log N) character columns. Matrix-vector products have B (log N)^{O(1)} cost. Kernel sampling therefore costs L^{2a+o(1)}; character evaluation and each of the constantly many square-root tests cost L^{a+o(1)}. A relation's worst attempt cost is

  L^{max(1/sqrt2, a, theta/(2a))+o(1)}.

After the field-selection calculation, the total bound is

  L^{max(2a, a+1/sqrt2, a+theta/(2a))+o(1)}.

Let a*=theta/sqrt2. At a*, the last two terms are both (1+theta)/sqrt2 and 2a*=sqrt2 theta is smaller when theta<1. For a>a*, a+1/sqrt2 is larger. For 0<a<a*, the function a+theta/(2a) is strictly decreasing up to sqrt(theta/2), and a*<sqrt(theta/2) for theta<1. It is therefore larger than its value at a*. This independently proves the asserted global minimum.

At this optimum S^2=L^{1/sqrt2+o(1)} and B=L^{theta/sqrt2+o(1)}, so batch preprocessing fits and the row cost has exponent 1/sqrt2. No hidden fixed L-power was found in the algebraic postprocessing, Hensel precision, prime sampling, caching, or recursion.

thmA line 289 says each field-selection attempt costs T_shape+polylog. Such an attempt can invoke up to 2k special-q searches. Replace this by T_shape+2k T_SQ+polylog. Its exponent is 1/sqrt2 at the optimum and is below the relation-collection exponent, so the final bound is unchanged. The same uniform bound is what item 5 uses.

18. Uniformity in N and numerical constants — CORRECT WITH FIXES already identified for S and root sampling.

Location: thmA_proof.md §§1–2 and §6; lv6_results.md §7 effectivity.

Trial division to B leaves

  omega(N) <= log N/log B,
  log(2^omega) = O(sqrt(log N/log log N)) = o(sqrt(log N log log N)).

Thus all 2^{O(omega)} losses are uniformly L^{o(1)}. Testing primality and perfect powers removes prime powers before SPLIT; a composite not a perfect power has at least two distinct prime divisors. General nonsquarefree N needs no further assumption, as checked in DD and the ring congruence.

After trial division the cofactor uses its own parameters. Its prime factors exceed the previous B and hence its own, smaller B. There are at most O(log N) nontrivial recursive calls. Bounds uniform in the input integer, including LP's bound, can be applied to all smaller inputs using a monotone asymptotic envelope; this costs only another polylogarithm.

The cutoffs and constants may depend on the fixed parameters. The original-source effectivity confirmed for LO, and the effective LV restricted-modulus statement, do not introduce a family of exceptional N. The character construction uses an explicit discriminant upper bound rather than an unknown exceptional character. Fixed small inputs may be handled by an encoded finite cutoff and trial division. No quantitative claim here is uniform as k grows.

I independently minimized the maximum of the three affine h-floors, using their pairwise intersections, rather than importing the parameter scripts. The optimum is

  mu=(1-nu)/d, theta=2/d+(d^2+d-2)nu/d

below k=(d^2+1)/4, and

  mu=(1+(d-1)nu)/(d+1), theta=3/(d+1)+(4d-2)nu/(d+1)

above that transition, before adding eta0 to theta. At d=9,k=13 this gives mu=17/156, h=1/26+eta0, theta=16/39+eta0, and

  c=0.9972018709041053 + eta0/sqrt2.

At k=12 the minimum theta is 23/54, so it cannot beat 1. The slack allowance at k=13 is eta0<sqrt2-1-16/39=0.0039571521166849. The d=9 infimum is

  (1+3/10)/sqrt2 = 0.9192388155425117.

It is not achieved at any finite fixed k and positive eta0. Remark S1 correctly distinguishes this infimum from a single fixed-d algorithm with that exact leading constant plus o(1). No growing-k uniformity has been proved in the cited notes. The analogous fixed-degree numerical entries checked for d=7,11,19 also agree.

The corrections in this report are local: use the absorbing factor-exit calculation, formulate accepted roots by their fibres, reconcile the standalone K and LV6 statements with the final design, and correct the indicated averaged, cost, and computability statements. None requires a new analytic estimate or changes the displayed leading constant. My audit does not supply a separate exhaustive certification of every analytic argument behind FM and LD*.

WHOLE-THEOREM VERDICT: ACCEPT WITH FIXES
WEAKEST LINK: thmA_proof.md §4 Steps 6–7, the passage from raw-law relation supply to field selection with factor exits. The displayed gamma/2 inequality alone does not justify the claimed row-probability bound after removing successful factor-producing events. The absorbing-outcome calculation in item 5 repairs this interface without changing the exponent.
