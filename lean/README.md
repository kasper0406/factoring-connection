# FactoringLean: Lean 4 formalization of the Theorem D′ lemmas

This Lake project machine-checks parts of Theorem D′: the algebraic core (Lemmas 1, 3–7 and
the coverage criterion), the components of the window lemma (Lemma 2), a small-order bound with a
floored log (Lemma 8, partial), and the core counting step of the CRT/union argument (Lemma 10,
partial). Lemma 9 is proved only from a large-sieve hypothesis. See Theorem D′ (`../notes/dcltheorem_results.md`, §3–4). The status table (lemma, Lean name,
proved / hypothesis / not formalized) is in `../notes/leanform_results.md`.

- Toolchain: `leanprover/lean4:v4.35.0` (`lean-toolchain`).
- Dependency: Mathlib, tag `v4.35.0` (`lakefile.toml`; pinned commit in `lake-manifest.json`).
- No `sorry` and no `axiom` declarations. The arithmetic large sieve is **not** proved here; its prime-term weakening
  (Montgomery–Vaughan, MNT I Thm 19.13, keeping only `d = 1` and prime `d`) is assumed.
  It is the `Prop` `FactoringLean.ArithmeticLargeSievePrimeTerms`, and the two theorems that use it take
  it as an explicit hypothesis `(hLS : ArithmeticLargeSievePrimeTerms)`.

## Build

```sh
cd lean
lake exe cache get      # download prebuilt Mathlib .olean files (do not compile Mathlib)
lake build              # builds FactoringLean (expected: "Build completed successfully")
lake env lean AxiomCheck.lean   # prints the axioms of the 34 main results
```

`elan` lives in `~/.elan/bin`; put it on `PATH` first. `AxiomCheck.lean` should list only
`propext`, `Classical.choice` and `Quot.sound` for every result. `sorryAx` must never appear.
The `.lake/` build directory is git-ignored (`lean/.gitignore`).

## Files

| file | contents |
|---|---|
| `FactoringLean/Exponents.lean` | (a) exponent identity and Fermat consequence; Lemma 3 (mod p); Lemma 4 (mod q) |
| `FactoringLean/Coverage.lean` | (b) coverage criterion; (c) Lemma 1 class lock; distinct excluded classes (for Lemma 9) |
| `FactoringLean/Bands.lean` | Lemma 5 one-sided extraction; (e) Lemma 6 four-band pigeonhole; (f) Lemma 7 FULL criterion |
| `FactoringLean/SmallOrder.lean` | (g) Lemma 8: few primes with small `ord_q(2)` |
| `FactoringLean/CRTCount.lean` | (h) Lemma 10: residue-class count, CRT, union bound |
| `FactoringLean/LargeSieve.lean` | large-sieve **hypothesis**; Lemma 9 and the `E_b` bound deduced from it |
| `FactoringLean/SingleTarget.lean` | Lemma 2 (window); single-target extraction chaining Lemmas 3, 4, 5, 7; wrapper from Lemmas 1, 2 |
| `AxiomCheck.lean` | `#print axioms` for the main results |
