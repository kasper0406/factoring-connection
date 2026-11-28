"""Find balanced divisors with unlabelled polynomial or residue-list constraints."""

from math import gcd, prod

from sympy import Poly, Symbol

if __package__:
    from .affine_prime_window import _integer_roots, lll_reduce
    from .weighted_order_hint import _box, _multiplicity
else:
    from affine_prime_window import _integer_roots, lll_reduce
    from weighted_order_hint import _box, _multiplicity


T = Symbol('T')


def _parameters(n, lower, hints, radius, max_dimension, min_agreement=None):
    if min_agreement is not None:
        degree_modulus = prod(mod**f.degree() for mod, f in hints)
        largest_degree = max(f.degree() for _, f in hints)
    for d in range(2, max_dimension+1):
        m = _multiplicity(n, lower, d)
        if min_agreement is None:
            weights = [(d-1)//f.degree() for _, f in hints]
            covered = prod(mod**t for (mod, _), t in zip(hints, weights))
        else:
            t = min(_multiplicity(degree_modulus, min_agreement, d), d//largest_degree)
            weights = [t]*len(hints)
            covered = min_agreement**t
        bound = lower**m * covered
        det4 = radius**(2*d*(d-1)) * n**(2*m*(m+1))
        det4 *= prod(mod**(2*f.degree()*t*(t+1)) for (mod,f),t in zip(hints,weights))
        if (1 << (d*(d-1))) * d**(2*d) * det4 < bound**(4*d):
            return d, (m, *weights), bound
    raise ValueError('no list-polynomial lattice certificate through the dimension cap')


def _basis(equations, weights, center, radius, dimension):
    active = [(mod,f,t) for (mod,f),t in zip(equations,weights) if t]
    modulus = prod(mod**t for mod,_,t in active)
    coordinates = []
    for mod,f,t in active:
        local = mod**t
        cofactor = modulus//local
        projector = cofactor*pow(cofactor, -1, local)
        shifted = f.shift(center).trunc(local)
        powers = [[1]]
        power = Poly(1,T,domain='ZZ')
        for _ in range(t):
            power = (power*shifted).trunc(local)
            powers.append(list(map(int,reversed(power.all_coeffs()))))
        coordinates.append((mod,f.degree(),t,projector,powers))

    rows = []
    for j in range(dimension):
        leading = 1
        residues = [0]*j
        for mod,k,t,projector,powers in coordinates:
            u = min(t, j//k)
            leading *= mod**(t-u)
            shift = j-k*u
            for i, coeff in enumerate(powers[u][:-1]):
                residues[shift+i] += projector*coeff
        row = [leading*(v % modulus)*radius**i for i,v in enumerate(residues)]
        row.extend([leading*radius**j] + [0]*(dimension-j-1))
        for i in reversed(range(j)):
            pivot = rows[i][i]
            quotient = (row[i]+pivot//2)//pivot
            if quotient:
                for h in range(i+1):
                    row[h] -= quotient*rows[i][h]
        rows.append(row)
    return rows


def polynomial_list_divisors(n, hints, *, center=0, radius=None, max_dimension=None,
                            min_agreement=None, stats=None):
    """All balanced divisors r in the box satisfying F_i(r)=0 mod M_i.

    Each hint is (M_i, coefficients in increasing degree). Polynomials must be
    monic modulo M_i, of positive degree. All M_i and n must be pairwise coprime.
    With min_agreement=H, returns all qualifying divisors with
    product_i gcd(M_i,F_i(r)) >= H. No matching coordinates are provided.
    Raises ValueError when no determinant certificate is available at the cap.
    """
    product = 1
    equations = []
    for mod, coeffs in hints:
        coeffs = tuple(coeffs)
        if (not isinstance(mod,int) or mod < 2 or gcd(product*n,mod) != 1
                or len(coeffs) < 2 or any(not isinstance(c,int) for c in coeffs)
                or coeffs[-1] % mod != 1):
            raise ValueError('expected pairwise-coprime moduli and positive-degree monic polynomials')
        product *= mod
        equations.append((mod, Poly.from_list([v % mod for v in reversed(coeffs)],T,domain='ZZ')))
    lower, upper, radius = _box(n, product, center, radius)
    agreement = product if min_agreement is None else min_agreement
    if not isinstance(agreement, int) or not 1 <= agreement <= product:
        raise ValueError('min_agreement must be an integer between one and the modulus product')
    if max_dimension is None:
        max_dimension = max(8, n.bit_length()+product.bit_length()+2)
    if not isinstance(max_dimension,int) or max_dimension < 2:
        raise ValueError('dimension cap must be an integer of at least two')
    if stats is not None:
        stats.update(dimension=0, multiplicities=(0,)*(1+len(equations)), evaluation_bound_squared=0)
    if lower > upper:
        return ()
    if lower == upper:
        candidates = [lower]
    else:
        d, weights, bound = _parameters(n, lower, equations, radius, max_dimension,
                                        None if agreement == product else agreement)
        basis = _basis([(n,Poly(T,T)),*equations],weights,center,radius,d)
        row = list(map(int,lll_reduce(basis)[0]))
        norm_bound = d*sum(v*v for v in row)
        if not 0 < norm_bound < bound**2 or any(v % radius**i for i,v in enumerate(row)):
            raise ArithmeticError('LLL row failed its integer certificate')
        f = Poly.from_list([v//radius**i for i,v in enumerate(row)][::-1],T,domain='ZZ')
        candidates = [center+t for t in _integer_roots(f,radius)]
        if stats is not None:
            stats.update(dimension=d, multiplicities=weights, evaluation_bound_squared=norm_bound)
    return tuple(sorted(r for r in candidates if lower <= r <= upper and n % r == 0
                        and prod(gcd(mod, int(f.eval(r))) for mod, f in equations) >= agreement))


def residue_list_divisors(n, hints, *, min_agreement=None, **kwargs):
    """Require the product of moduli with literal list membership to reach H.

    Omitting min_agreement requires membership in every supplied nonempty list.
    """
    lists, equations = [], []
    for mod, values in hints:
        if not isinstance(mod,int) or mod < 2:
            raise ValueError('expected integer moduli of at least two')
        values = tuple(values)
        if not values or any(not isinstance(v,int) for v in values):
            raise ValueError('each list must contain one or more integer residues')
        residues = frozenset(v % mod for v in values)
        f = Poly(1,T,domain='ZZ')
        for a in sorted(residues):
            f = (f*Poly(T-a,T)).trunc(mod)
        equations.append((mod, [int(c) % mod for c in reversed(f.all_coeffs())]))
        lists.append((mod,residues))
    roots = polynomial_list_divisors(n, equations, min_agreement=min_agreement, **kwargs)
    agreement = prod(mod for mod, _ in lists) if min_agreement is None else min_agreement
    return tuple(r for r in roots
                 if prod(mod for mod, residues in lists if r % mod in residues) >= agreement)
