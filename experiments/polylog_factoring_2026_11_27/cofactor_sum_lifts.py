"""Exact factor-sum precision and a signed Fermat search on dyadic progressions."""

from math import gcd, isqrt


def state_sum_class(n, endpoint, bits):
    """Maximal power-of-two congruence shared by x+n/x for odd x in the state."""
    if n < 1 or n % 2 == 0 or endpoint % 2 == 0 or bits < 1:
        raise ValueError("expected an odd positive product, an odd residue and positive precision")
    modulus = 1 << bits
    endpoint %= modulus
    common = gcd(endpoint*endpoint-n, 2*modulus)
    if common < modulus:
        stride = modulus*common
    elif common == modulus:
        stride = 2*modulus**2
    else:
        stride = modulus**2
    residue = (endpoint+n*pow(endpoint, -1, stride)) % stride
    return residue, stride


def _square_residues(value, levels):
    """Successively lift square roots modulo 1, 4, 16, ... without restarting."""
    roots, modulus = [0], 1
    for _ in range(levels):
        yield modulus, roots
        lifted_modulus = 4*modulus
        roots = [r for root in roots for r in range(root, lifted_modulus, modulus)
                 if (r*r-value) % lifted_modulus == 0]
        if not roots:
            return
        modulus = lifted_modulus


def signed_fermat_factor(n, budget=1):
    """Try budget centers in each signed square-root progression using n alone."""
    if n < 3 or n % 2 == 0 or budget < 1:
        raise ValueError("expected an odd input >=3 and a positive budget")
    for sign in (1, -1):
        lower = isqrt(n-1)+1 if sign == 1 else 0
        for modulus, roots in _square_residues(sign*n, n.bit_length()):
            for root in roots:
                first = lower+(root-lower) % modulus
                for i in range(budget):
                    center = first+i*modulus
                    if center > (n+1)//2:
                        break
                    delta = center*center-sign*n
                    other = isqrt(delta)
                    if other*other == delta:
                        divisor = abs(center-other)
                        if 1 < divisor < n:
                            return divisor
    return None
