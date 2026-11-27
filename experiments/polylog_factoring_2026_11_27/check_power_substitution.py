"""Check ordered Frobenius products for the cubic T^3-T-1."""

from math import gcd
from sympy import Matrix
from check_lehman_window_baseline import cohorts


def power_matrix(index, n):
    def multiply(a, b):
        c = [0]*5
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                c[i+j] += x*y
        for k in (4,3):
            c[k-3] += c[k]
            c[k-2] += c[k]
        return [v%n for v in c[:3]]

    t = [0,1,0]
    v = [1,0,0]
    for bit in bin(index)[2:]:
        v = multiply(v,v)
        if bit == '1':
            v = multiply(v,t)
    return Matrix([[1,0,0],v,multiply(v,v)]).T


if __name__ == '__main__':
    for size,planted,cases in cohorts():
        if planted or size == 16:
            continue
        hits = 0
        for p,q in cases:
            n = p*q
            c = power_matrix(n,n)
            hits += any(1 < gcd(int(v),n) < n for v in [*c,c.det()])
            a,b = power_matrix(p,n),power_matrix(q,n)
            assert all(v%p==0 for v in c-a*b)
            assert all(v%q==0 for v in c-b*a)
        print(size,len(cases),hits)
