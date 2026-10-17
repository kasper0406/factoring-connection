"""Track A (e): masked observables of F_N via quotient algebras.  OPERATIONAL.

Each function returns the masked information about F_N that the cheapest known
algebra in the isolation model provides, with its rank, so costs can be
compared with the lower bounds of Theorem A1 / Corollary E1:

  prefix(N, K)           coefficients a_0..a_{K-1}     algebra R_N[X]/(X^K)        rank K
  interval(N, a, b)      coefficients a_a..a_{b-1}     algebra R_N[X]/(X^b)        rank b
  residue_class(N,r,M,c) S_{r,M}(c)                    algebra R_N[X]/(X^M - c)    rank M
  bit_pattern(N,j,c)     sum over classes r with bit j of r set, of S_{r,2^{j+1}}(c)
                         = masked sketch (P (.) m)(c^{1/2^{j+1}}-type)   rank 2^{j+1}
The mask for bit j has zero runs of length 2^j, so by the gap lemma any
isolating observable needs rank >= 2^j + 1; prefix/interval need rank >= the
last nontrivial index + 1 over a field (Fitting form, squarefree N).
"""
from __future__ import annotations

from typing import List, Optional

from ..common.metrics import Metrics
from ..common.ring import FlintQuotient, power_of_linear
from .sections import section_values


def prefix(N: int, K: int, m: Optional[Metrics] = None) -> List[int]:
    Q = FlintQuotient(N, [0] * K + [1], m)
    v = Q.coeffs(power_of_linear(Q, N))
    v[0] = (v[0] - 1) % N
    return v


def interval(N: int, a: int, b: int, m: Optional[Metrics] = None) -> List[int]:
    return prefix(N, b, m)[a:b]


def residue_class(N: int, r: int, M: int, c: int, m: Optional[Metrics] = None) -> int:
    return section_values(N, M, c, m)[r]


def bit_pattern(N: int, j: int, c: int, m: Optional[Metrics] = None) -> int:
    M = 1 << (j + 1)
    vals = section_values(N, M, c, m)
    return sum(vals[r] for r in range(M) if (r >> j) & 1) % N
