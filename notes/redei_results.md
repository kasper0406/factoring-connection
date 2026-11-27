# redei: 2-power class-group bits of discriminants attached to N (2026-11-27)

Agent `redei`. Prereg `experiments/redei/prereg.md` (eaf9bfa). Code/data `experiments/redei/`. Codex: `codex/redei_1.md`
-> `codex/redei_1_codex.out` (gpt-6.2-astra, session 01a21e3a-9a6b-74f3-879f-876dbc9758db).
All ledger code is diagnostic (knows p, q). No factor-blind positive result exists, so there is no 35-bit+ positive test.

**Verdict.** We found no native poly(log N) 2-power class-group bit beyond genus. The 4-rank family is a QR_N oracle,
and the 8-rank family is a hidden-component QR oracle. Neither has a known factoring decoder. The one family
that provably gives a factor is Kronecker–Hurwitz, and it has 4√N members, each with a non-native class number. Nothing
here goes below L[1/3].

## 1. Bit ledger (which hidden symbols each invariant carries)

| Invariant, family | Hidden content | Bits | Status |
|---|---|---|---|
| r4 of a single Cl(−4N) or Cl(−8N), or h mod 8 | f(ε, (2/p), (p/q)), e.g. [2 is a square mod N] (genusleak §2.5); for N≡1 (8), [r4 Cl(−4N)=0] = a_N(Δ⁹ mod 2) = [p≡q≡3 (8)] (codex Q1/Q4) | ≤ 2 | proved; O(1) bits, so brute-forceable |
| r4 of Q(√−kN), k ∈ {1, r, rs}·{1, 2}, primes r ≤ B | σ = ((−1/p), (2/p), (q/p), (r/p) for r ≤ B), determined up to the p↔q swap | π(B)+2 | Rédei formula proved; E1 gives 600/600 agreement with PARI. With pairs rs, E2b gives exactly H−1 bits in 60/60 cases. Singletons alone give 4.1/6, 8.5/12 and 18.3/27 bits at B = 10, 30, 100 (E2) |
| the same family with a prime multiplier a ≈ 2⁴⁰ | (a/p) for any a. **The 4-rank family is a QR_N oracle**: ≤ follows from Rédei, ≥ from E2c (39/40) | full | proved (≤); experimentally supported (≥) |
| does the fingerprint (r/p), r ≤ B, pin down p? | Information-theoretically yes: B ≤ (log N)² under GRH (codex Q2). E6: it is injective on all primes in [2^b, 2^{b+1}] after 20, 22, 26, 33 characters at b = 14, 16, 18, 20, which is about 2 log₂ #primes | ~log N bits | no decoder is known. QR ⇒ factoring is the open Gap–Rabin problem (codex Q2) |
| r8, Rédei symbols [a,b,p] | Artin bit of p in a D4/C4 field. Sketch: this equals χ_𝔭(x+y√a) on elements whose norm is a square mod p, i.e. a restricted Legendre-PRF oracle with hidden key s = √a mod p. Recovering s factors N via gcd(s²−a, N). Also p mod 16 in some cells (tower T3) | O(1) per (a,b) | reduction: conjectured (sketch). Published Legendre-PRF key recovery costs p^{Ω(1)} (from memory) |
| r8, r16 and full 2-Sylow types of 8 discriminants −kN → p mod 32 | none detected. Held-out accuracy on s=p+q bit 7: r4+r8 0.739 vs null 0.718; r16 0.743 vs 0.733; full 0.775 vs 0.770. The control bit 6 (p mod 16) shows signal, 0.831 vs 0.744 (E3) | 0 | experimentally supported. Codex gives a certified pair (367·751 and 383·991, both 161 mod 256) with identical 2-parts of Cl(−4N) and Cl(−8N) that differ mod 32. It also proves that finitely many fixed-slot Rédei symbols have governing exponent ≤ 4, so they cannot reach Q(ζ₃₂) |
| Coppersmith-usable residues | Needs p mod M with M > N^{1/4}. 2-power data gives only the 2-Sylow projection of p in (Z/M)*, about 1 bit per r instead of log₂ r, and p mod 2^j only up to j = 4 | O(1) | proved for fixed-slot symbols; growing twists are open |
| **Kronecker–Hurwitz**: Σ_{t²<4N} H(4N−t²) = 2(N+q) | H mod 2^k over the family gives q mod 2^{k−1}, which factors N at k ≈ log₂N/4 (Coppersmith). Mod 8 it already gives the Blum bit (−1/q) | full | proved (classical); 40/40 exact checks (E4). The family has 4√N members and about 15% are nonzero mod 8. Each term needs the class group of 4N−t² |

Note on E3: the first run sampled N with replacement and showed spurious signal at p mod 32. The cause was train/test
duplicates being enriched among key hits. The rerun with distinct N is the one reported above.

## 2. Native routes: where each one needs the factors

| Route | Failure point | Native cost |
|---|---|---|
| Composition on public forms (ℓ, b, c) of small norm | They generate Cl under GRH, but [(ℓ,b,c)] ∈ Cl² ⇔ (ℓ/p)=(ℓ/q)=1, which is QR. The 2-Sylow structure needs relations or h | L[1/2] (Hafner–McCurley), or N^{1/4} |
| Square root / 2-descent of the public T=(2,2,(N+1)/2) | **Proved:** if B² = T with B=(a,·,·), gcd(a,N)=1, then T(x,y)=a² gives (2x+y)²+Ny² = 2a², so (2x+y)/a is √2 mod N. Two roots differing by a hidden ambiguous class factor N. The 8-rank also needs χ_p(B)=(a/p), which is QR again | needs √ mod N |
| Rédei reciprocity / governing fields | [a,b,N] = [a,b,p][a,b,q] is public, the Jacobi analogue. [a,b,p] needs the prime above p in Q(√a), i.e. √a mod p. [a,b,p] = [a,p,b] needs the field F_{a,p}, i.e. the conic x²−ay² = pz² | needs √a mod p |
| h mod 2^k via Σ_{a<\|D\|/2} χ(a) | Gauss's lemma linearises (a/N) only mod 2, so it reaches #{χ=−1} mod 2, i.e. h mod 4 (genus). Bit 3 and above are nonlinear in the Gauss counts | Θ(N); exact h costs N^{1/4} by BSGS |
| Hurwitz / r₃ | r₃(N) = 12h(−4N), so h mod 2^k ⇔ r₃ mod 2^{k+2}, which is a lattice-point count mod 2^k (the roadmap's open counting core) | power of N |
| Nicolas–Serre mod-2 Hecke | a_N(Δ^e) = parity of positive odd solutions of N = Σ_{i∈e} 2^i m_i². Δ⁹ gives the hidden 4-rank bit. No fixed e is both native and non-public (codex Q1) | binary-form count mod 2, which needs Cornacchia (√ mod N) |
| Kronecker–Hurwitz, Eichler–Selberg mod 2^k | The trace formula sums over all t² < 4N. Restricting t to a progression mod M moves the work into a trace on level-M² forms | ≥ N^{1/2}/M terms plus a trace of cost ≥ M |

## 3. Own ideas (both killed)
- **Short Kronecker–Hurwitz window** (E5/E5b). The centre window |t| ≤ 40 is null (perm p=0.21). The edge window t ≈ 2√N
  does leak q mod 4 (perm p<0.005), but the whole leak comes from square discriminants 4N−t² = (2s)². Those are just
  sum-of-two-squares representations N = (t/2)²+s² found by search: 6H ≡ 2 (4) occurred 85 times with q≡1 and 0 times
  with q≡3 (saved run, `e5b_edge.out.json`; an earlier run gave 87/0). Held-out accuracy overall was 0.512 vs 0.467 base
  (0.516 vs 0.511 in the earlier run). This is Fermat-type search, costing √N in general.
- **Bounded norm search for r4>0 on N≡1 (8)** (codex Q4). Here [r4>0] ⇔ N = x²+y² or N = x²−2y². A cutoff y ≤ log₂N
  gives 87 TP / 0 FP / 342 FN on all pairs 101 ≤ p < q ≤ 499, and it can recognise only O(√X log X) inputs up to X.
  Reproduced via `experiments/redei/codex1_check.gp`.

## 4. Sharpest open step
To decode the hidden key s = √a mod p from the Rédei-symbol oracle χ_𝔭(x+ys), with N's factors and p unknown, in
poly(log N) time. Even with a native 8-rank evaluator, the decoding problem is Legendre-PRF key recovery, and the
known attacks cost p^{Ω(1)}. For 4-ranks the open step is Gap–Rabin, i.e. QR ⇒ factoring. In both cases the native
evaluator itself is also missing.
