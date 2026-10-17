# Track C — audit of Tal (ECCC TR26-211) for an unknown-factor modulus

Owner: factC. Date 2026-10-17.

**Sources consulted.** I read these myself:
- Y. Tal, *Computing Modular Factorials Below the Square-Root Barrier*, ECCC TR26-211. This is the 69-page PDF from https://eccc.weizmann.ac.il/report/2026/211/download, dated September 2026 and posted 25–26 Sep 2026. I read §§1–5, §6.1–6.2, §8 and App. H. I read the other appendices only where cited below.
- Bostan–Gaudry–Schost, SIAM J. Comput. 36(6) 2007, in the author-hosted PDF. I read Lemma 3, Lemma 4, Thm 5, Thm 8, Lemma 12, Lemma 13 and Thm 11.
- `codex/r1_tal_bgs.out` (GPT-6.1-sol). I checked its key claims against the PDFs (see §4).

**Notation.** Below, `p` is Tal's *known* prime and `q` is Tal's auxiliary conductor. For our problem `N = P·Q`, the factors are unknown.

## 1. Tal's results in one paragraph

Thm 1.1 (and Thm 1.2 mod `p²`) assume `q | p−1` and compute `n! mod p` in `Õ(q^c + sqrt(p/q))`. Thm 1.3 assumes `q | Φ_r(p)` for a fixed odd prime power `r`, and works mod `p^k` in `Õ_r(k^{O(1)}[q^{c_r} + (Φ_r(p)/q)^{1/(2φ(r))}])`. The exponents `c` and `c_r` are not made explicit (§8). The saving comes from recovering factorial **ratios at special arguments** from Jacobi/Gauss sums. Those sums are reconstructed from their Stickelberger ideal plus a polarization (Lenstra–Silverberg, Thm 2.7), and the root of unity is then fixed (van Wamelen, Lemma 2.8). What remains is short interval products done with BGS (Lemma 2.6). Prop. 4.3 / App. H is a separate result: a uniform `sqrt(log p/log log p)` saving with no divisor assumption, via modular halving plus BGS Thm 15.

## 2. Step-by-step table

Columns: what each step needs, and whether it survives over `Z/NZ` with unknown factors, treating failed inversions as outputs.

| # | Step (Tal ref.) | Needs | Survives over Z/NZ? |
|---|---|---|---|
| 1 | Choose `q` with `q \| p−1` (or `q \| Φ_r(p)`, `q \| p+1`); Prop. 1.6 picks it from the factorisation of `p−1` (resp. `Φ_r(p)`) | **known p**, factorisation of `p−1` / `Φ_r(p)` | **No.** `q \| N−1` or `q \| Φ_r(N)` says nothing about `P−1` or `Q−1`. Knowing a `q \| P−1` is Pollard-`p−1` information. |
| 2 | Representation of `F_p` or `F_{p^r}` and an element `η` of exact order `q`. Deterministic setup costs `p^{1/4+o(1)}` (power-nonresidue search, l. 668/1684); Las Vegas uses `x^{(p^r−1)/q}` | **known p, group order** `p^r−1` | **No.** `x^{(N−1)/q}` mod N does not have order `q` in either CRT component unless `q \| P−1`. You can test componentwise order with `gcd(η^{q/s}−1, N)`, but a pass does not create the field. |
| 3 | Character `χ_m` defined by `ρ(χ_m(x)) = x^{−(p^f−1)/m}` (§2.2) | **character of known order**, exponent `(p^f−1)/m` | **No.** Over Z/NZ the only characters computable from N alone are powers of the Jacobi symbol (order 2). |
| 4 | Gauss/Jacobi identities (Prop. 2.1, eqs. (9), (10), (11)); norm `J·J̄ = p`; App. A uses `1_{y=0} = 1 − y^{p−1}` (l. 2659) and `Σ_x x^d` over `F_p` | **field `F_p`, Fermat, character orthogonality** | **No.** These are false identities mod N, not merely missing inverses. Example: `1 − 3^{14} ≡ 7 (mod 15)` although `3 ≠ 0`. |
| 5 | Stickelberger ideal `∏ p_{m,t}^{c_t}` (§2.4, eq. (15)/(16), Prop. 2.4) | **known p**: the primes `p_{m,t} = (p, ζ_m − η^t)` | **No.** It needs the primes of `Q(ζ_m)` above `p`. Primes above `N` would be the ideals `(N, …)`, which are not prime, and factoring them is factoring N. |
| 6 | Polarized generator recovery (Lenstra–Silverberg, Thm 2.7) | CM-order, ideal `I`, `w = p^f` | The algorithm itself is generic (poly time in its input). It has no meaningful input here, because step 5 does not exist. |
| 7 | Root-of-unity orientation (Lemma 2.8, van Wamelen Thm 9) | Sums at proper divisors of the conductor, the congruence `Σ j a_j ≡ … (mod m)` | Generic integer arithmetic. Nothing to orient without step 6. |
| 8 | Local reduction `ρ_{p^k}` of a product representation (Lemma 2.5): clear denominators, **remove the known `p`-adic valuation**, Hensel-lift | **known uniformizer `p`, local prime, valuations** | **No.** Removing a valuation is not an inversion that could fail informatively. |
| 9 | Decoder (2): `(aK)! ≡ (−1)^{a−1} ρ(J_a)(K!)^a (mod p)` with `K = (p−1)/q` | **known p** (both the block size `K` and the identity) | **No.** No analogous congruence mod N; `(p−1)/q` is unknown. |
| 10 | Multinomial decoder (Prop. 5.1), via first-order Gross–Koblitz (13)/(14) and base-`p` digits | **base-p digits, Wilson, valuation `d`** | **No.** Base-N digits do not encode the hidden-prime digits. |
| 11 | Cyclotomic lattice rounding (Lemma 5.2: LLL in dim `φ(r)`, rotation by `ζ_r`) | integers `p`, `D \| Φ_r(p)` | Survives as an **integer** algorithm with `N` in place of `p`. But the lattice point has no decoder (step 10 fails), so it is useless. |
| 12 | Wilson complementation `n! ≡ (−1)^{m+1}(m!)^{−1}` (§4.2) | `(p−1)! ≡ −1` | **No.** `(N−1)! ≡ 0 (mod N)` for composite N > 4. Only the trivial reflection `j ↦ N−j` survives, and it preserves unit status. |
| 13 | Mod-`p²` harmonic correction (Prop. 6.1, Lemma 6.2: `H_{aK} − H_K = Σ (η^{at} − η^t) Q_t`, `Q_t = ((1−η̃^t)^{p−1} − 1)/p`) | **known p**, Fermat quotient, lifted roots | **No.** |
| 14 | Shifted factorial polynomial `P_M(Z) = ∏(j + pZ) mod p^k`, degree `< k` (§2.7, §7) | **known p** (the coefficient of `Z^i` is divisible by `p^i`) | **No** as stated. The analog `∏(j + N·Z) mod N^k` exists, but it is a statement about `N`-adic lifts of an interval product, not a way to shorten it. |
| 15 | BGS interval product (Lemma 2.6) and batched products with shifts (Lemma H.2) | ring with monic division, unit small integers | **Yes.** These are exactly our baselines (b) and (c). Cost is `Õ(sqrt(total length))`; see trackC_analysis §1. |
| 16 | Modular halving (eq. (20), App. H): `∏_{a<j≤b} j = 2^{b−a} ∏(evens/2) ∏((odds+p)/2)` | `2` a unit | **Survives algebraically for odd N, but gives nothing** (proof below). |
| 17 | BGS Thm 15 (selected indices of a 2×2 polynomial-matrix recurrence) | units among `O(sqrt H)` small integers | **Yes.** It is generic, with gcd-checked inverses and sqrt cost. |

**Why step 16 gives nothing over Z/NZ.** `proved` (elementary). After `t` halving rounds, an interval `[A, A+L)` with `L ≪ N` maps to `2^t` arithmetic sub-progressions. They are rescaled to sit near `(A + kN)/2^t` for `k = 0..2^t−1`, with pairwise separation about `N/2^t ≫ L`. The map is a bijection on residues and multiplies by a unit, so unit status is preserved and the total length stays exactly `L`.

Tal's saving comes from routing branches into one *shared* prefix `[1, H]` with `H ≈ p/log p`, whose values are batch-evaluated once (Lemma H.1). That requires the product length to be comparable to the modulus, so that branches wrap around into the prefix. For our intervals (`L ≤ sqrt N ≪ N`), no two branches overlap and none falls into a common short prefix except the `k = 0` branch (the multiples of `2^t`). The work is therefore unchanged.

**Sub-step question.** Does any step compute a product of linear factors `∏_{j∈S}(x − j)` mod N, or a factorial-type quantity, below sqrt cost without knowing the prime?

**No** (`refuted for this model`). Every sub-square-root step is one of:
- reconstructing a prime-specific factorial ratio (steps 2–10 and 13), or
- a log-factor rearrangement whose gain needs product length ≍ modulus (step 16).

The only generic product routines (steps 15 and 17) are square-root in the length. This agrees with `codex/r1_tal_bgs.out` §3.

## 3. What would be needed to transfer anything

1. **A source of factorial ratios mod N.** The special values `(aK)!`, `C(A−1, B−1)` and so on come from evaluating Gauss sums at a prime above `p`. A CRT-compatible analog would need a ring hom `Z[ζ_m] → Z/NZ` whose kernel lies above *both* `P` and `Q` with compatible characters. That needs `m | P−1` and `m | Q−1` simultaneously, plus knowledge of those characters, which is order information about `(Z/NZ)^*`. So the transfer is blocked at step 2 already. It is the classical "known group order" obstruction of `p−1`-type methods, in a new guise.
2. **A divisor `q` of `p−1` of size `p^{δ}`.** For random RSA primes this exists only with smoothness-type probability, and it is unknown. Tal's §1.3 density analysis is over primes `p` that are *given*.
3. Tal's own remark §1.4: the CRT direction only helps when the factors are known. No converse (fast factoring ⇒ fast `n! mod p`) is known.

## 4. Check of `codex/r1_tal_bgs.out` (= curated `notes/codex_r1_tal_bgs.md`) against the PDFs

Orchestrator-requested spot checks were done directly in the Tal PDF text: Prop. 6.1/Lemma 6.2, App. A indicator identity (l. 2659), the `p^{1/4+o(1)}` setup (l. 668/1684), Lemma H.1/H.2, and §1.4 (Lipton/Shub–Smale/Cheng remarks). All matched.

| codex claim | Check | Status |
|---|---|---|
| Thm 1.1–1.3 bounds, `c` and `c_r` not explicit | §1.1, §8 | confirmed |
| Deterministic setup `p^{1/4+o(1)}` via power nonresidues | l. 668, proof of Thm 1.1 (l. 1684) | confirmed |
| App. A uses `1_{y=0} = 1 − y^{p−1}` and `Σ x^d` | l. 2659 | confirmed |
| Prop. 6.1 (mod `p²` with `1 + (ap/q)(H_{aK} − H_K)`) and Lemma 6.2 (roots-of-unity filter, `q−1` exponentiations) | §6.1–6.2 | confirmed |
| Lemma H.1 (≤ `H^{3/8}` endpoints, BGS Thm 15 with ε = 1/10), Lemma H.2 (`O((s+1+sqrt U)L^6)`) | App. H | confirmed |
| BGS Thm 5 condition: `β, 2..d, α−dβ..α+dβ` units; Thm 8 / Lemma 13 / Lemma 12 as described | BGS pp. 1784–1791 | confirmed |
| Modular halving "routing survives for odd N, but at scale N" | my analysis above | confirmed, sharpened: zero saving for `L ≪ N` |
| Codex Idea 1: HH restricted to a certified interval near `sqrt N`, cost `Õ(W^{3/5} N^{−1/10})` | HH not read by me | **not yet tested.** It tests location of the *prime* `P`, not arbitrary multiples. It is not a unit-only mechanism. Recorded as card C-HHloc. |
| Harvey/HH details (§5 of codex) | not read by me | UNVERIFIED by factC |

## 5. Bottom line

`proved`/`refuted for this model`: no step of Tal's construction yields a sub-square-root unit test over `Z/NZ` with unknown factors.

The necessary ingredients are field identities in a known characteristic and characters of known order dividing `p−1` (or `Φ_r(p)`). The failure is not a missing inversion that could leak a gcd. The identities themselves are false mod N.

The generic parts that survive (BGS shifting, Thm 15, lattice rounding, halving) are either square-root cost or give no saving for `L ≪ N`.
