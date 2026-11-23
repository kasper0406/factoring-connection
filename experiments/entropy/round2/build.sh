#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
cxx="${CXX:-clang++}"
gmp="${GMP_PREFIX:-/opt/homebrew/opt/gmp}"
ssl="${OPENSSL_PREFIX:-/opt/homebrew/opt/openssl}"
args=(-O3 -std=c++17 -Wno-deprecated-literal-operator -Wno-deprecated-declarations -I"$gmp/include" -I"$ssl/include")
for name in budget nodes identity engine; do
    "$cxx" "${args[@]}" "$name.cpp" -L"$gmp/lib" -lgmpxx -lgmp -o "$name"
done
"$cxx" "${args[@]}" sampler.cpp -L"$ssl/lib" -lcrypto -o sampler
"$cxx" "${args[@]}" locked.cpp -L"$gmp/lib" -L"$ssl/lib" -lgmpxx -lgmp -lcrypto -o locked
for name in matched dyad; do "$cxx" "${args[@]}" -shared -fPIC "$name.cpp" -o "$name.dylib"; done
./engine names > variants.txt
