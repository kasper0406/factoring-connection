# FactoringLean: Lean 4 formalization of the Theorem D′ lemmas

This Lake project machine-checks the algebraic and elementary number-theoretic lemmas of
Theorem D′ (`../notes/dcltheorem_results.md`, §3–4). The status table (lemma, Lean name,
proved / hypothesis / not formalized) is in `../notes/leanform_results.md`.

- Toolchain: `leanprover/lean4:v4.35.0` (`lean-toolchain`).
- Dependency: Mathlib, tag `v4.35.0` (`lakefile.toml`; pinned commit in `lake-manifest.json`).
- No `sorry` and no `axiom` declarations. The arithmetic large sieve is **not** proved here.
  It is the `Prop` `FactoringLean.ArithmeticLargeSieve`, and the two theorems that use it take
  it as an explicit hypothesis `(hLS : ArithmeticLargeSieve)`.

## Build

```sh
cd lean
lake exe cache get      # download prebuilt Mathlib .olean files (do not compile Mathlib)
lake build              # builds FactoringLean (expected: "Build completed successfully")
lake env lean AxiomCheck.lean   # prints the axioms of the 31 main results
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
| `FactoringLean/SingleTarget.lean` | Lemma 2 (window); single-target extraction chaining Lemmas 3, 4, 5, 7 |
| `AxiomCheck.lean` | `#print axioms` for the main results |
