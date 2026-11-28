"""Use factorizations of smaller signed low words to produce candidate factor residues."""

from math import isqrt

from sympy import divisors

if __package__:
    from .hh_low_residue import hh_low_factor
    from .hh_residue_lists import factor_from_residue_list
else:
    from hh_low_residue import hh_low_factor
    from hh_residue_lists import factor_from_residue_list


def centered(value, modulus):
    return (value+modulus//2) % modulus-modulus//2


def word_modulus(n):
    return 1 << max(2, 5*n.bit_length()//24)


def short_product_classes(n, modulus, budget=0):
    """All odd residues whose centered inverse partner has product carry in [-B,B]."""
    if n < 1 or n % 2 == 0 or modulus < 4 or modulus & (modulus-1) or budget < 0:
        raise ValueError("expected an odd product, a dyadic word >=4, and a nonnegative budget")
    head = centered(n, modulus)
    classes = set()
    for carry in range(-budget, budget+1):
        word = head+carry*modulus
        for divisor in divisors(abs(word), generator=True):
            x = int(divisor)
            if 2*x < modulus and 2*abs(word//x) < modulus:
                classes.add(x % modulus)
                classes.add(-x % modulus)
    return sorted(classes)


def factor_from_short_products(n, budget=0, *, stats=None):
    """Try the low-product promise from n alone on a balanced odd semiprime."""
    modulus = word_modulus(n)
    classes = short_product_classes(n, modulus, budget)
    work = {} if stats is None else stats
    work.update(modulus=modulus, smaller_integers=2*budget+1)
    return factor_from_residue_list(n, modulus, classes, stats=work)


def near_product_classes(n, modulus, carry_budget, center_budget):
    """Find suffixes whose unsigned product has a square in the specified Fermat range."""
    if (n < 1 or n % 2 == 0 or modulus < 4 or modulus & (modulus-1)
            or carry_budget < 0 or center_budget < 1):
        raise ValueError("expected an odd input, a dyadic word and valid search budgets")
    head = centered(n, modulus)
    classes = set()
    for carry in range(-carry_budget, carry_budget+1):
        word = abs(head+carry*modulus)
        lower = isqrt(word-1)+1
        for midpoint in range(lower, lower+center_budget):
            excess = midpoint*midpoint-word
            root = isqrt(excess)
            if root*root != excess:
                continue
            small, large = midpoint-root, midpoint+root
            if 2*large < modulus:
                classes.update(d % modulus for d in (small, -small, large, -large))
    return sorted(classes)


def near_word_modulus(n):
    # A fixed margin keeps the affine lattices small on executable controls.
    return 1 << (12+(n.bit_length()+3)//4)


def factor_from_near_products(n, carry_budget=1, center_budget=1, *, stats=None):
    """Split a balanced semiprime whose short suffix product has close absolute factors."""
    m = near_word_modulus(n)
    classes = near_product_classes(n, m, carry_budget, center_budget)
    work = {} if stats is None else stats
    work.update(modulus=m, classes=len(classes), affine_windows=0)
    for residue in classes:
        local = {}
        divisor = hh_low_factor(n, m, residue, stats=local)
        work["affine_windows"] += local["intervals"]
        if divisor is not None:
            return divisor
    return None
