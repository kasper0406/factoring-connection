#!/bin/bash
# Sandboxed paranoid check of leanA inside a Linux container (Docker Desktop on macOS).
#
#   docker build -t leana-checker leanA/docker      # once
#   leanA/docker/check.sh                           # from the repository root
#
# The container gets a fresh copy of the sources (no .lake from the host). Mathlib's .olean files come
# from `lake exe cache get` (the Mathlib cache, trusted as in the Comparator README); `LeanA` and
# `Challenge` are built from source. Then:
#   1. `lake env lean AxiomCheck.lean`        (#print axioms for all MAIN theorems)
#   2. `lake comparator --paranoid`           (Comparator, built into Lean 4.35: statement matching
#                                              against Challenge.lean, axiom whitelist, and replay of the
#                                              solution through Lean's kernel, leanchecker-paranoid,
#                                              lean4lean, nanoda, con-leche and con-ron, inside a
#                                              bubblewrap sandbox with an empty network namespace)
# Docker's own seccomp/AppArmor/masked-/proc confinement is relaxed only so that bubblewrap can create
# its namespaces; the judged code runs inside bubblewrap.
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
out="${1:-$here/docker/last_check.log}"
docker volume create leana-mathlib-cache >/dev/null
docker volume create leana-lake >/dev/null
docker run --rm \
  --security-opt seccomp=unconfined --security-opt apparmor=unconfined \
  --security-opt systempaths=unconfined \
  -v "$here":/src:ro -v leana-mathlib-cache:/home/checker/.cache -v leana-lake:/home/checker/leanA/.lake \
  leana-checker bash -c '
    set -euo pipefail
    mkdir -p ~/leanA && cd /src && tar --exclude=./.lake --exclude=./docker/last_check.log -cf - . | (cd ~/leanA && tar xf -)
    cd ~/leanA
    rm -rf .lake/build   # packages (Mathlib checkout) persist in a volume; our own build is always fresh
    export PATH="$(dirname "$(elan which lean)")":$PATH   # bundled nanoda_bin, lean4lean, ... on PATH
    echo "== toolchain"; lean --version; lake --version; which nanoda_bin
    echo "== mathlib cache"; lake exe cache get 2>&1 | tail -3
    echo "== build"; lake build 2>&1 | grep -v "^✔\|^⚠\|imports .Mathlib. from package\|declaration uses .sorry." | tail -20 || true
    echo "== sorry scan (proof layer)"; grep -rn "sorry" LeanA/ LeanA.lean || echo "no sorry in LeanA/"
    echo "== axiom scan (proof layer)"; grep -rnE "^\s*(axiom|opaque|unsafe|@\[implemented_by|@\[extern)" LeanA/ || echo "no axiom/opaque/unsafe/implemented_by/extern in LeanA/"
    echo "== #print axioms"; (lake env lean AxiomCheck.lean 2>&1 || true) | sort | uniq -c
    echo "== lake comparator --paranoid"; if lake comparator --paranoid; then echo "COMPARATOR_EXIT=0"; else echo "COMPARATOR_EXIT=$?"; fi
  ' 2>&1 | tee "$out"
