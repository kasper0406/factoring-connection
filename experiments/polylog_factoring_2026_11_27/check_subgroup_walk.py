"""Exact transition, reversibility, and first-hit checks for the subgroup walk."""

from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random

from sympy import Matrix, eye, ones, totient

if __package__:
    from .check_sublattice_sampler import hnfs, profile_of
    from .subgroup_walk import acceptance, proposal, subgroup_step
else:
    from check_sublattice_sampler import hnfs, profile_of
    from subgroup_walk import acceptance, proposal, subgroup_step


def rank_row(p, q, d, r, s, *, raw=False):
    """Rank-pair transition law on F_p^d times F_q^d."""
    row = Counter()
    old = p**(d-r)*q**(d-s)
    for intersect in (False, True):
        stay_p = Fraction(1, p**(r if intersect else d-r))
        stay_q = Fraction(1, q**(s if intersect else d-s))
        direction = -1 if intersect else 1
        for i,j in product((0,1), repeat=2):
            mass = ((1-stay_p) if i else stay_p) * ((1-stay_q) if j else stay_q) / 2
            if not mass:
                continue
            u, v = r+direction*i, s+direction*j
            new = p**(d-u)*q**(d-v)
            prob = Fraction(1) if raw else Fraction(*acceptance(p*q,d,old,new,intersect=intersect))
            row[u,v] += mass*prob
            row[r,s] += mass*(1-prob)
    assert sum(row.values()) == 1
    return +row


def finite_walk(n, d, pq=None):
    matrices = {tuple(h): h for h in hnfs(n,d)}
    indices = {key: abs(int(h.det())) for key,h in matrices.items()}
    ambient = n**d
    if pq:
        labels = {key: tuple(profile_of(h,p,1)[0] for p in pq) for key,h in matrices.items()}
    transition = {}
    expanded = 0
    for key,h in matrices.items():
        row = Counter()
        for intersect in (False, True):
            counts = Counter(tuple(proposal(n,h,v,intersect=intersect))
                             for v in product(range(n), repeat=d))
            for other,count in counts.items():
                old, new = indices[key], indices[other]
                if other != key:
                    quotient = max(old,new)//min(old,new)
                    expected = old*int(totient(quotient)) if intersect else ambient//old*int(totient(quotient))
                    assert count == expected
                proposal_mass = Fraction(count, 2*ambient)
                accept = Fraction(*acceptance(n,d,old,new,intersect=intersect))
                row[other] += proposal_mass*accept
                row[key] += proposal_mass*(1-accept)
            expanded += ambient
        assert sum(row.values()) == 1
        transition[key] = row
        if pq:
            ranked = Counter()
            for other,prob in row.items():
                ranked[labels[other]] += prob
            assert +ranked == rank_row(*pq,d,*labels[key])
    for key,row in transition.items():
        for other,prob in row.items():
            assert transition[other][key] == prob
    return len(matrices), expanded


def first_hit(p, q, d):
    epsilon = Fraction(1,p**((d+1)//2)) + Fraction(1,q**((d+1)//2))
    rows = [rank_row(p,q,d,r,r) for r in range(d+1)]
    for row in rows:
        escape = sum(prob for (u,v),prob in row.items() if u != v)
        assert escape <= epsilon
    aligned = Matrix(d+1, d+1, lambda r,s: rows[r][s,s])
    expectations = (eye(d+1)-aligned).inv() * ones(d+1,1)
    assert min(expectations) >= 1/epsilon
    assert max(expectations) <= 9*d + 40/epsilon
    return expectations[d//2]


def run():
    groups = proposals = 0
    for n,d,pq in ((4,1,None), (12,1,None), (4,2,None), (9,2,None),
                   (6,1,(2,3)), (6,2,(2,3))):
        count, moves = finite_walk(n,d,pq)
        groups += count
        proposals += moves
    laws = 0
    for p,q in ((2,3), (3,5), (17,19), (101,103)):
        for d in (1,3,5,7):
            expected = first_hit(p,q,d)
            laws += 1
            if (p,q,d) == (101,103,3):
                print(f"central initial rank at N={p*q}, d={d}: E[T]={float(expected):.3f}")
    raw = rank_row(3,5,3,1,1,raw=True)[1,0]
    accepted = rank_row(3,5,3,1,1)[1,0]
    assert raw == Fraction(2,15) and accepted == raw/75

    n, h, rng = 35, 35*eye(3), Random(106)
    for _ in range(32):
        h = subgroup_step(n,h,rng)
        assert all(x.is_Integer for x in n*h.inv())
    print(f"{groups} subgroup states, {proposals} raw proposals, exact detailed balance")
    print(f"{laws} exact odd-rank first-hit systems and 32 public HNF steps")


if __name__ == "__main__":
    run()
