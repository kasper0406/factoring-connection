"""Public arithmetic probes on the branch 4N = c² + 59."""

import json
from argparse import ArgumentParser
from math import gcd, isqrt
from random import Random

from flint import fmpz, fmpz_mod_poly_ctx


def public_parameters(n):
    c = isqrt(4 * n)
    if c * c + 59 != 4 * n or gcd(n, 6 * 59) != 1 or c % 2 != 1:
        raise ValueError(
            "expected an odd canonical branch of discriminant -59 coprime to 6 and 59"
        )
    # beta = u + v*omega, with omega² + omega + 1 = 0.
    u, v = (c + 3) // 2, 3
    beta_norm = u * u - u * v + v * v
    t = (c - 3) // 2
    assert beta_norm == t * t + 3 * t + 9 == n - 8 and n % 9 == 8
    return {
        "c": c,
        "beta": (u, v),
        "beta_norm": beta_norm,
        "shanks_parameter": t,
        # Coefficients in ascending order of the trace equation for a norm-one
        # square root of the Shanks generator. Entries are simplified modulo N.
        "halving_quartic": [t - 5, -8, -2 * t, 0, 1],
        # By cubic reciprocity the unit and ramified supplements are 1 at
        # N = -1 mod 9, and N reduces to the cube 8 modulo the prime-to-3 part of beta.
        "cubic_jacobi_exponent": 0,
    }


class CyclicCubic:
    """Shanks algebra with its order-three automorphism over Z/nZ."""

    def __init__(self, n, t):
        self.n = n
        self.ring = fmpz_mod_poly_ctx(n)
        z = self.ring.gen()
        self.f = z**3 - t * z**2 - (t + 3) * z - 1
        self.sigma_z = -(z**2) + t * z + t + 2  # -(z+1)/z, since z is a unit

    def conjugates(self, a):
        b = a.compose_mod(self.sigma_z, self.f)
        return b, b.compose_mod(self.sigma_z, self.f)

    def norm(self, a):
        b, d = self.conjugates(a)
        return self.norm_from_conjugates(a, b, d)

    def norm_from_conjugates(self, a, b, d):
        norm = ((a * b % self.f) * d) % self.f
        if norm.degree() > 0:
            raise ArithmeticError("norm failed to descend to the base ring")
        return int(norm[0])

    def sample_from(self, a):
        """Return sigma(a)/a for a unit; otherwise return gcd(N, Norm(a))."""
        b, d = self.conjugates(a)
        norm = self.norm_from_conjugates(a, b, d)
        divisor = gcd(norm, self.n)
        if divisor != 1:
            return None, divisor
        return (((b * b % self.f) * d % self.f) * pow(norm, -1, self.n)), None

    def shifted_symbols(self, u):
        norms = [self.norm(1 + u), self.norm(1 + (u * u % self.f))]
        signs = [fmpz(v).jacobi(self.n) for v in norms]
        divisors = {gcd(v, self.n) for v in norms}
        return tuple(map(int, signs)), sorted(d for d in divisors if 1 < d < self.n)

    def public_half_probe(self, quartic):
        """Try the inert-prime halving exponent with the public modulus."""
        n = self.n
        z = self.ring.gen()
        v = z.pow_mod((n * n + n + 2) // 2, self.f)
        v1, v2 = self.conjugates(v)
        trace = int((v + v1 + v2)[0])
        error = (v * v % self.f) - z
        # A prime dividing every coefficient is a component where v is a half.
        half_gcd = gcd(n, *(int(error[i]) for i in range(3)))
        root_gcd = gcd(n, int(self.ring(quartic)(trace)))
        return {"trace": trace, "half_gcd": half_gcd, "quartic_gcd": root_gcd}


def public_probe(n, draws, seed):
    if draws < 0:
        raise ValueError("draw count must be nonnegative")
    params = public_parameters(n)
    algebra = CyclicCubic(n, params["shanks_parameter"])
    half = algebra.public_half_probe(params["halving_quartic"])
    rng = Random(seed)
    pairs = []
    factors = {d for d in (half["half_gcd"], half["quartic_gcd"]) if 1 < d < n}
    for _ in range(draws):
        a = algebra.ring([rng.randrange(n) for _ in range(3)])
        u, divisor = algebra.sample_from(a)
        if divisor is not None:
            if 1 < divisor < n:
                factors.add(divisor)
            continue
        signs, divisors = algebra.shifted_symbols(u)
        pairs.append(signs)
        factors.update(divisors)
    return {
        "N": n,
        **params,
        "public_half": half,
        "attempts": draws,
        "shifted_symbols": pairs,
        "gcd_factors": sorted(factors),
    }


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("N", type=int)
    p.add_argument("--draws", type=int, default=8)
    p.add_argument("--seed", type=int, default=59)
    args = p.parse_args()
    print(json.dumps(public_probe(args.N, args.draws, args.seed), indent=2))
