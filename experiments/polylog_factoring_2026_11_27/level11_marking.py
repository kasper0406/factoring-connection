"""Classify a line ideal with a projective mod-2 or an integral mod-2^k frame."""

from math import gcd
from random import SystemRandom

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

if __package__:
    from .level11_automorphism import _product, _conjugate
    from .level11_sampler import _inverse_splitting, _line_basis, _shell_points
else:
    from level11_automorphism import _product, _conjugate
    from level11_sampler import _inverse_splitting, _line_basis, _shell_points


# At 3 use omega=diag(1,0), j=[[0,1],[-1,0]]. The line (1,1) imposes
# x+y+z = z+w-x = 0 mod 3. The norm-3 points in O have coordinate image
# one of the two axes, so this right ideal J represents the six-unit class.
_J = ((1, -1, 0, 1), (0, -1, 1, -1), (0, 3, 0, 0), (0, 0, 0, 3))
_BAR_J = tuple(_conjugate(v) for v in _J)
_LINES = ((1, 0), (0, 1), (1, 1))


def _witness_mark(point, principal):
    # In O/2 the matrices for omega and j are [[0,1],[1,1]] and
    # [[0,1],[1,0]]. The three nonzero vectors are the projective lines.
    x, y, z, w = point
    a, b, c, d = x+w, y+z, y+z+w, x+y+w
    images = [_LINES.index(((a*u+b*v) % 2, (c*u+d*v) % 2)) for u, v in _LINES]
    if principal:
        # The projective units of O fix (1,1), the unique principal 2-neighbor.
        return images[2]
    # The projective units of the left order of J form the even permutations.
    inversions = sum(images[i] > images[j] for i in range(3) for j in range(i+1, 3))
    return 3 + inversions % 2


def _class_witnesses(n, inverse, line, base):
    """Return the class and all integral numerators of its reference transports."""
    basis = _line_basis(n, inverse, line)
    if base == "nonprincipal":
        # gamma is in J and has norm 6. J/gamma O has size 4, so gamma
        # generates J at every odd prime. K = gamma I + n J agrees with
        # gamma I at r|n and J elsewhere, and retains the fixed 2-frame.
        gamma = _J[0]
        products = [_product(gamma, tuple(map(int, basis.col(k)))) for k in range(4)]
        basis = hermite_normal_form(Matrix(products).T.row_join(n*Matrix(_J).T))
        n *= 3  # absolute norm of the relative norm-n subideal of J
    points = _shell_points(basis, n)
    if points:
        if len(points) != 4:
            raise ArithmeticError("expected four generators of a principal ideal")
        return True, points

    # I*bar(J) is exactly {beta : beta J ⊂ 3I}. At norm 3n the
    # inclusion has equal indices, so (beta/3) J = I. This holds at
    # primes dividing both n and 3 as well. The factor 3 is a unit mod 2.
    columns = [tuple(map(int, basis.col(k))) for k in range(4)]
    transport = hermite_normal_form(Matrix([_product(a, b) for a in columns for b in _BAR_J]).T)
    witnesses = _shell_points(transport, 3*n)
    if len(witnesses) != 6:
        raise ArithmeticError("expected the six class equivalences of the nonprincipal ideal")
    return False, witnesses


def _marked_line(n, inverse, line, base="principal"):
    principal, points = _class_witnesses(n, inverse, line, base)
    return _witness_mark(points[0], principal)


def _setup(n, rng, base):
    if n < 1 or gcd(n, 22) != 1:
        raise ValueError("expected a positive odd norm coprime to 11")
    if base not in ("principal", "nonprincipal"):
        raise ValueError("expected a principal or nonprincipal base")
    return _inverse_splitting(n, SystemRandom() if rng is None else rng)


def _line(n, u, v):
    if gcd(u, v, n) != 1:
        raise ValueError("the coordinates must generate the unit ideal modulo n")
    return u % n, v % n


def marked_line_classifier(n, rng=None, *, base="principal"):
    """Return a function (u,v) -> one of five marked classes for odd n coprime to 11.

    Values 0,1,2 give the distinguished line of the principal class in
    the order (1,0), (0,1), (1,1). In the six-unit class, 3 denotes even
    and 4 odd transport from the fixed norm-three representative.
    The principal base starts at state 2, the nonprincipal base at state 3.
    """
    inverse = _setup(n, rng, base)

    def classify(u, v):
        return _marked_line(n, inverse, _line(n, u, v), base)

    return classify


def framed_line_classifier(n, k, rng=None, *, base="principal"):
    """Retain the integral 2^k-frame of a relative norm-n line ideal, k >= 2.

    The result is (target class, row-major matrix) where the matrix is the
    lexicographically least transport in its right coset of reference units.
    For absolute ideal norm n_abs and target representative norm d in {1,3},
    the determinant is n_abs/d mod 2^k. Here n_abs equals n from the
    principal base and 3n from the nonprincipal base.
    """
    if not isinstance(k, int) or k < 2:
        raise ValueError("expected a marking exponent of at least two")
    inverse = _setup(n, rng, base)
    modulus = 1 << k
    # The 2-adic splitting lifts the same root at every level. Its matrix
    # here takes four order coordinates to four row-major matrix entries.
    split_inverse = _inverse_splitting(modulus, rng)
    forward = split_inverse.adjugate()*pow(int(split_inverse.det()), -1, modulus)
    rows = [tuple(int(a) % modulus for a in row) for row in forward.tolist()]
    third = pow(3, -1, modulus)

    def classify(u, v):
        principal, witnesses = _class_witnesses(n, inverse, _line(n, u, v), base)
        # Each transport is beta/d. Enumerating beta lists exactly its
        # right translates by units of the reference left order, four or
        # six elements irrespective of k. No matrix group is enumerated.
        scale = 1 if principal else third
        frame = min(tuple(scale*sum(a*b for a,b in zip(row, beta)) % modulus
                          for row in rows) for beta in witnesses)
        return ("principal" if principal else "nonprincipal"), frame

    return classify
