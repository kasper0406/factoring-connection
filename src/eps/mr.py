"""Miller-Rabin / Solovay-Strassen certificates of eps = +1.  OPERATIONAL: input N only.

Source: GPT-6.1-sol r6 (notes/codex_r6_primality_lit.md section 9), relayed by the orchestrator before the freeze.
Write N - 1 = 2^s d, d odd (s >= 2 on D_n).  For base b: x_0 = b^d, x_i = x_{i-1}^2 mod N.
  PC4 (MR):  x_i = -1 (mod N) for some 1 <= i <= s-1  =>  x_{i-1} is a square root of -1 mod N  =>  eps = +1.
  PC5 (SS):  b^{(N-1)/2} = -1 (mod N)  (this is x_{s-1} = -1, so PC5 implies PC4 for the same base).
A nontrivial square root of 1 met in the chain (x_i = 1, x_{i-1} != +-1) gives a factor gcd(x_{i-1} - 1, N);
it is recorded (rule 3), never used.
"""
from __future__ import annotations

from math import gcd
from typing import Dict

BASES = range(2, 66)


def mr_features(N: int) -> Dict[str, int]:
    s, d = 0, N - 1
    while d % 2 == 0:
        d //= 2
        s += 1
    pc4 = pc5 = 0
    splits = 0
    for b in BASES:
        x = pow(b, d, N)
        for i in range(1, s):
            y = x * x % N
            if y == N - 1:
                pc4 = 1
                if i == s - 1:
                    pc5 = 1
            if y == 1 and x not in (1, N - 1):
                if 1 < gcd(x - 1, N) < N:
                    splits += 1
            x = y
    return {"pc4_mr": pc4, "pc5_ss": pc5, "mr_split": splits}
