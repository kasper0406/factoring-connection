"""A Metropolis walk uniform at stationarity on all subgroups of (Z/n)^d."""

from random import SystemRandom

from sympy import Matrix, zeros
from sympy.matrices.normalforms import hermite_normal_form

if __package__:
    from .idempotent_sampler import bernoulli
else:
    from idempotent_sampler import bernoulli


def proposal(n, lattice, vector, *, intersect):
    """Return the HNF after adding a vector or restricting by its dot product.

    Input columns span the inverse image in Z^d of the current subgroup,
    so this full-rank lattice must contain n Z^d.
    """
    d = lattice.rows
    if n < 1 or d <= 0 or lattice.cols != d or len(vector) != d:
        raise ValueError("expected a positive modulus and matching lattice and vector")
    index = abs(int(lattice.det()))
    v = Matrix(vector)
    if not intersect:
        return hermite_normal_form(lattice.row_join(v), D=index)
    # The augmented lattice is the graph of the character modulo n. In
    # column HNF exactly the first d columns span its section at height 0.
    graph = lattice.row_join(zeros(d, 1)).col_join((v.T*lattice).row_join(Matrix([[n]])))
    return hermite_normal_form(graph, D=n*index)[:d,:d]


def acceptance(n, dimension, old_index, new_index, *, intersect):
    """Integer numerator and denominator of the Metropolis acceptance probability."""
    ambient = n**dimension
    indices = old_index * new_index
    numerator, denominator = (ambient, indices) if intersect else (indices, ambient)
    return min(numerator, denominator), denominator


def subgroup_step(n, lattice, rng=None):
    """One factor-free update of the inverse-image lattice in column HNF."""
    rng = SystemRandom() if rng is None else rng
    intersect = bool(rng.randrange(2))
    vector = [rng.randrange(n) for _ in range(lattice.rows)]
    candidate = proposal(n, lattice, vector, intersect=intersect)
    if candidate == lattice:
        return lattice
    a, b = acceptance(n, lattice.rows, abs(int(lattice.det())), abs(int(candidate.det())),
                      intersect=intersect)
    return candidate if bernoulli(a, b, rng) else lattice
