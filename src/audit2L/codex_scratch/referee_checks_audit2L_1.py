"""Independent finite and exponent checks for the blind two-level referee report."""
from fractions import Fraction as F
from math import gcd


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def small_part(n, bound):
    n, out = abs(n), 1
    for p in range(2, bound + 1):
        if all(p % d for d in range(2, int(p**0.5) + 1)):
            while n % p == 0:
                out *= p
                n //= p
    return out


def smooth(n, bound):
    return small_part(n, bound) == abs(n)


print('RAT EXC(c): a mid-prime contribution to s_i is missed by a small-part-only argument')
p1, p2, g = (1, 1), (1, 210), (1, 23)
n1, n2 = det(p1, g), det(p2, g)
K, P0, S0, y, q = abs(det(p1, p2)), 5, 7, 17, 19
delta_mid = 11
assert all(gcd(*v) == 1 for v in [p1, p2, g])
assert K == q * delta_mid and q > y
assert smooth(n1, y) and smooth(n2, y)
assert small_part(n1, P0) <= S0 and small_part(n2, P0) <= S0
assert small_part(n1, P0) * delta_mid > S0
assert small_part(n2, P0) * delta_mid > S0
print(dict(P1=p1, P2=p2, g=g, values=(n1, n2), K=K,
           small_parts=(small_part(n1, P0), small_part(n2, P0)),
           s=(small_part(n1, P0)*delta_mid, small_part(n2, P0)*delta_mid)))

print('Parent PRIM: p > y dividing g can give two nonzero parent indicators even when p does not divide K')
p1, p2, g = (4, 5), (3, 7), (23, 23)
K, y, q, ell, Y, V, M = abs(det(p1, p2)), 7, 13, 23, 20, 7, 23
vals = (det(p1, g), det(p2, g))
assert K == q and K % ell != 0 and Y < ell <= 2*Y
assert all(smooth(n // ell, y) for n in vals)
assert all(F(M*V, 32) <= abs(n) <= 2*M*V for n in vals)
assert gcd(*g) == ell > y
assert p1[0] * pow(p1[1], -1, q) % q == p2[0] * pow(p2[1], -1, q) % q
print(dict(P1=p1, P2=p2, g=g, values=vals, K=K, content=gcd(*g)))

print('Exponents: delta, theta, theta_p, correct budget margins, and effect of ONE extra log N in B')
for al in [F(3, 8), F(7, 18), F(39, 100), F(399, 1000)]:
    delta = 3 - 7*al
    th = max(3*al-1, 1-2*al, 2*al-F(1,2))
    thp = max(1-2*al, 2-5*al)
    assert th < al and thp < delta
    assert 3*al-F(1,2) > 4*al-1
    assert 2-5*al > 0
    assert delta-2*al < 0
    assert -al < 0
    print(dict(alpha=str(al), delta=str(delta), theta=str(th), theta_p=str(thp),
               A_over_eta_power=str(delta-2*al),
               sqrt_AB_over_eta_power=str(-al),
               sqrt_AB_over_eta_with_extra_logN_power=str(F(1,2)-al)))
