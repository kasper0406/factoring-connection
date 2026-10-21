# sihard round 1 (BLIND): hardness of constructing smooth numbers in short intervals near kN

You are a research mathematician (computational number theory, smooth numbers, lattice / Coppersmith methods, CRT list
decoding, complexity theory: NP-hardness, TFNP, fine-grained) collaborating on a speculative project about
asymptotically faster classical integer factorization. Rigour matters more than speed, and honesty more than a positive
result. Label every claim PROVED (give the proof), KNOWN (cite precisely: authors, venue, year, theorem/section number —
check with web search; write UNVERIFIED if you could not check), HEURISTIC (name the assumption), model-PROVED (exact
consequence of a named model), CONJECTURED or OPEN. Do not read other files of this repository; work from this prompt,
your knowledge and web search. Scratch computation (Python, exact integer/rational arithmetic) goes in the current
directory. Final answer: a structured report (≤ 7000 words) with numbered sections matching the tasks, ending with a
one-paragraph "sharpest barrier statement".

## Notation and the target

n = log N, F_α = n^α (log n)^{1−α}, L[α, c] = exp((c + o(1)) F_α). Records: GNFS L[1/3, 1.923], MNFS L[1/3, 1.9019],
SNFS (special N) 1.526. ρ = Dickman function.

**Prop SI (the target, a sufficient condition found earlier in this project).** Factor base y = L[1/3, κ]; intervals
I_k = [kN, kN + h], h = L[2/3, ν], k = 1, 2, …, L[1/3, O(1)]. Suppose an algorithm ("constructor") outputs, in
L^{o(1)} time per call, a y-smooth u_k ∈ I_k (it may be a promise constructor, valid only when I_k contains a y-smooth
integer). The mean number of y-smooth integers in I_k is L[2/3, ν − 2/(3κ)], so ν > 2/(3κ) is needed. Each u_k gives a
relation u_k ≡ v_k (mod N) with v_k = u_k − kN ∈ [0, h]; v_k is y-smooth with probability L[1/3, −ν/(3κ)] (heuristic:
the constructor does not bias v_k). Collect L[1/3, κ] relations, sparse linear algebra L[1/3, 2κ]. Cost
L[1/3, max(κ + ν/(3κ), 2κ)], optimised at κ³ = 2/9, ν ↓ 2/(3κ): **L[1/3, (16/9)^{1/3} = 1.2114]**. If each call
costs L[1/3, τ], the program beats MNFS iff τ < 0.757 (at κ = (4/9)^{1/3}). Equivalently: the constructor solves a
one-sided Schnorr-type inhomogeneous approximation problem Σ e_p log p ∈ [log kN, log(kN + h)], e ≥ 0, over the
L[1/3]-dimensional prime lattice.

Your job is the HARDNESS side: is there any reason to believe such a constructor cannot exist, at least in natural
models, and what would it imply?

## Tasks

1. **Reductions and consequences.**
   (a) Note that a constructor that works for every input integer N is, with k = 1, a constructor for y-smooth
   integers in an arbitrary interval [x, x + L_x[2/3, ν]] at y = L_x[1/3, κ]. What does such a modulus-generic
   constructor imply beyond factoring? Work out carefully: prime-field discrete logarithms (relation collection with
   N → p; linear algebra mod p − 1; **individual logarithms**: initial splitting ("booting") and descent of medium
   primes — can the constructor be used in the descent, e.g. by querying it with a different modulus M ≈ jp/q? What
   is the resulting L[1/3] constant, and does the booting step (ECM-based, standard constant 3^{1/3}?) dominate?);
   class groups of imaginary quadratic fields; discrete logs in other groups; deterministic construction of smooth
   numbers in short intervals at other scales (which scales follow trivially, which do not).
   (b) Conversely: is Prop SI implied by fast factoring (a factoring oracle), or by anything known? Does a factoring
   oracle help construct smooth numbers in short intervals? Formalise a query model in which it provably does not.
2. **Literature status of "construct a y-smooth integer in a given interval [x, x + h]".** Survey: trivial sampling /
   sieving costs; Boneh, "Finding smooth integers in short intervals using CRT decoding" (STOC 2000; JCSS 64 (2002)) —
   exact theorem, parameter range, and what it gives at the Prop SI scale (x ≈ N·L[1/3], h = L[2/3], y = L[1/3]);
   Coppersmith / Howgrave-Graham "approximate integer common divisors" and Guruswami–Sahai–Sudan weighted CRT decoding;
   Granville's survey "Smooth numbers: computational number theory and beyond" (MSRI 44, 2008), especially any section
   on finding smooth numbers in short intervals; existence results (Hildebrand, Friedlander–Granville, Croot,
   Matomäki–Radziwiłł, and any 2020s results); deterministic construction results; any explicit statement that this is
   an open problem. Any hardness results (NP-hardness, relation to knapsack / subset sum / subset product) for finding
   smooth numbers in intervals, or for related problems (divisors in intervals, Manders–Adleman, subset product,
   Garey–Johnson)?
3. **Model lower bounds.** In natural algorithm models, prove that the constructor needs ≥ L[1/3, c] per output for
   some c > 0 (or more, e.g. L[2/3, ·]) at the Prop SI scale, and say whether the 1.2114 target then collapses to
   ≥ 1.9019. Candidate models:
   (a) sampling / sieving (the algorithm examines integers and tests them);
   (b) lattice / Coppersmith–Howgrave-Graham–Boneh certificate decoders: a call on a region [x, x + V] with modulus
   P (a product of factor-base prime powers) outputs all z with gcd(z, P) ≥ S, guaranteed when
   log V ≤ (log S)²/log P (unweighted Johnson-type bound), or the GSS weighted analogue. Bound the expected number
   of *y-smooth* certified outputs per call, at V ≤ h = L[2/3, ν], for every choice of P, in a random-divisibility or
   averaged-over-x model. Does any P beat the trivial "z = d·⌈x/d⌉ with d smooth" trick by more than L[2/3]^{o(1)}?
   Batching across many k (bivariate regions)?
   (c) product / meet-in-the-middle / generalized-birthday (Wagner k-tree) over factor-base products, including
   amortisation over many k (the target set is then the Bohr set {u : u mod N ∈ [0, h]}).
   (d) the actual Schnorr log-lattice CVP route (approximation factor and dimension needed).
4. **NP-hardness / TFNP.** Is "given (x, h, y), find a y-smooth integer in [x, x + h]" (or the decision version) NP-hard
   for some ranges of (y, h)? Consider: very short h (poly(log x)) — NP ∩ coNP?; ranges where every interval provably
   contains a smooth number (TFNP, Megiddo–Papadimitriou); the Prop SI scale (subexponential sieve time; existence only
   conjectured); a variant where the allowed primes form an explicit given set S (subset product / exact cover); bounded
   exponent variants. State exactly which ranges are provably NP-hard, which are provably not NP-hard unless
   NP = coNP (or ETH fails), and which are open, and where Prop SI's scale falls.
5. **Barrier statement.** The sharpest honest statement: either "the constructor, in models X, Y, Z, costs ≥ …", or
   "no barrier is known, and here is why it is plausible or implausible". Give your probability that an L[1/3, τ]
   constructor with τ < 0.757 exists.
