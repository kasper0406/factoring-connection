# suffmap round 1 (BLIND): the map of SUFFICIENT conditions for asymptotically faster classical factoring

You are a research mathematician (computational number theory, smooth numbers, algebraic number theory, lattice
algorithms, coding theory, algebraic complexity) collaborating on a speculative project about asymptotically faster
classical integer factorization. Rigour matters more than speed, and honesty more than a positive result. Label every
claim PROVED (give the proof), KNOWN (cite precisely: authors, venue, year, theorem/section number — check with web
search; write UNVERIFIED if you could not check), HEURISTIC (name the assumption), model-PROVED (exact consequence of a
named model), or OPEN. Do not read other files of this repository; work from this prompt, your knowledge and web
search. Scratch computation (Python, exact integer/rational arithmetic) goes in the current directory. Final answer: a
structured report (≤ 7000 words) with numbered sections matching the tasks, ending with a table
"condition | gives | relation to known problem | believed | status of the implication".

## Notation and background

x = log N, F_α = x^α (log x)^{1−α}, L[α, c] = exp((c + o(1)) F_α). Single N, RAM model, no precomputation unless
stated. Classical records: GNFS L[1/3, (64/9)^{1/3} = 1.923], Coppersmith MNFS L[1/3, 1.9019], SNFS (special N)
L[1/3, (32/9)^{1/3} = 1.526]; Coppersmith's factorization factory 1.6386 per N after an L[1/3, 2.0069]
precomputation. Rigorous record L[1/2, 1] (Lenstra–Pomerance 1992).

**Standard heuristic H.** An examined candidate whose smoothness integers have total log-size λ is B-smooth with
probability ρ(λ/log B)·L^{o(1)} (Canfield–Erdős–Pomerance); joint events multiply; one smoothness test costs L^{o(1)};
sparse linear algebra on K columns costs K^{2+o(1)}; dependencies are fruitful (split N with probability ≥ 1/2).

Literature is full of *barriers* (necessary conditions: "any faster method of type X would need Y"). We want the
**sufficient** side: concrete mathematical statements whose truth (plus at most the standard heuristic H) would
**directly** give faster factoring, ideally tied to recognised open problems — and, conversely, a precise
classification of famous open problems that are often mentioned near factoring but are actually unrelated.

## Tasks

1. **Smoothness-size principle as a heuristic theorem.** State and prove (within H) a precise two-sided statement:
   (a) sufficiency: a candidate generator supplying enough candidates of total log-size λ (each in L^{o(1)} time, with a
   multiplicative relation mod N) gives factoring in exp((c+o(1))·√(λ log λ)) with an explicit c; in L-notation, with
   λ = s·F_{2α}, give the exact constant as a function of (α, s), and the minimum number of candidates required;
   (b) necessity within the sampling model; (c) the role of the *size–entropy law* λ(C) (how the candidate size grows
   with the log of the number of candidates C), with a formula that reads off the exponent α from λ(C) (e.g. if
   λ(C) ≍ x^a C^b) and the constant (check that it reproduces GNFS 1.923 and SNFS 1.526 from their norm ledgers);
   (d) the exact σ needed (candidates of log-size ≍ x^σ) for α = 1/4, 1/5, 1/10, quasi-polynomial and polynomial
   time; (e) what changes for *listing* (output-sensitive) rather than *sampling* methods, i.e. when the cost is the
   number of outputs, not the number of examined candidates: one-sided (one side smooth by construction, the other
   small and tested) versus doubly-smooth listing. Who first stated each piece (Pomerance's L(X)^{√2} analysis, etc.)?

2. **Sufficient structural properties.** List every structural property you know or can derive that would yield
   candidates of the needed size together with a congruence to N, and for each give: the exact statement; what α or
   L[1/3] constant it gives (with derivation); whether it is implied by / equivalent to / implies a recognised problem
   or conjecture; whether that problem is believed true or false, and why. At least:
   - polynomial pairs (f, g) with a common root mod N and coefficient heights below the generic (Gaussian-heuristic)
     bound; what constant each height profile gives; can heights alone ever change the exponent 1/3? Is the generic
     bound a theorem for almost all N?
   - Kirchner 2018 ("NFS factorization: new hopes", Prikl. Diskr. Mat. Suppl. 11): number fields of special shape
     (ramified at N, small height, Galois subfield); what exactly is conditionally claimed;
   - number fields / orders of small discriminant with a ring map to Z/N;
   - explicit smooth numbers near multiples of N (u ≡ v mod N with u smooth by construction and |v| small);
   - doubly-smooth listers (u ≡ v mod N, both y-smooth), and their NFS analogue (pairs (a, b) with both norms smooth);
   - output-sensitive enumeration of smooth values of a − bm (rational side) or of F(a, b) (algebraic side);
   - higher-rank / multivariate candidate families, "systematic equations" (integer analogue of x^q − x = ∏(x − a)
     from quasi-polynomial finite-field DL), dense Prouhet–Tarry–Escott-type families;
   - anything from algebraic complexity (short straight-line programs for n!, Strassen 1976, Lipton 1994, Shub–Smale
     τ-conjecture, Bürgisser 2009 on the permanent and τ(n!)) and from complexity classes (NP ∩ coNP, TFNP classes).
   - group-order methods (ECM-type): what property of a family of algebraic groups over Z/N would give L_p[α], α < 1/2?

3. **High-profile candidates: decide precisely sufficient / necessary / unrelated** (and if "sufficient", for what: an
   exponent below 1/3, polynomial time, or only an L[1/3] constant):
   - smooth numbers in short intervals or sparse sets (Matomäki–Radziwiłł 2016, Harper, "every interval [x, x + x^ε]
     contains an x^ε-smooth number"), including the open problem of *deterministically constructing* a y-smooth
     number in a given interval [X, X + h] — at what (X, h, y) would a polynomial-time construction matter, and what
     would it give? (Hint to examine, not to assume: can a search oracle be turned into an output-sensitive lister?)
   - polynomial-time deterministic construction of primes (Polymath4, Tao–Croot–Helfgott 2012) or of smooth numbers;
   - sum-product phenomena and expansion (Cayley graphs of (Z/N)^* generated by small primes; navigation/path finding);
   - existence of number fields with small discriminant and a prescribed root mod N;
   - explicit S-unit or abc statements (Evertse-type counts; effective abc); algorithmic S-unit equation solving;
   - approximate CVP in Schnorr–Adleman prime-number lattices (only classify; another agent works this out in detail);
   - Coppersmith / partial-key-exposure bounds (factoring from a fraction of the bits of p);
   - list decoding of Chinese-remainder codes beyond the Johnson radius (Goldreich–Ron–Sudan, Boneh 2000,
     Guruswami–Sahai–Sudan; the 2026 plain-RS capacity result of Brakensiek–Chen–Putterman–Zhang–Zheng if you can find
     it); sub-quadratic sparse linear algebra over F₂ (vs Wiedemann/Lanczos, Peng–Vempala 2021).

4. **Rank** the sufficient targets by directness × plausibility × profile. Which is most worth attacking? Which would
   be a landmark *either way* (a proof or a disproof would be important)?

5. Give 3–6 small exact numeric sanity checks of your key formulas (state inputs and outputs, e.g. the constants that
   fall out of the size–entropy formula for Dixon, QS, GNFS, SNFS, and a table of σ against α).
