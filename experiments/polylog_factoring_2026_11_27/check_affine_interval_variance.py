"""Check the affine sign covariance and the prime-power variance factor."""

from fractions import Fraction
from math import gcd


def sign(n, value):
    return 1 if value % n <= (n-1)//2 else -1


def check_covariance(n):
    divisors = [d for d in range(1,n+1) if n%d==0]
    for d in divisors:
        total = sum(sign(n,b)*sign(n,a*d+b) for a in range(n) for b in range(n))
        assert total == d*d
    return len(divisors)


def check_local(p,e):
    modulus = p**e
    units = [x for x in range(1,modulus) if gcd(x,modulus)==1]
    weights = {x:gcd(x*x-1,modulus) for x in units}
    direct = sum(weights[x]*weights[y]*gcd(x-y,modulus)**2 for x in units for y in units)
    by_classes = 0
    for k in range(e+1):
        divisor = p**k
        masses = [0]*divisor
        for x,value in weights.items():
            masses[x%divisor] += value
        jordan = 1 if k==0 else divisor**2-(divisor//p)**2
        by_classes += jordan*sum(value**2 for value in masses)
    assert direct == by_classes
    exact = Fraction(direct,modulus**4)
    assert exact <= 2+Fraction(10,p)
    print(p,e,f"{float(exact):.8f}")


if __name__ == "__main__":
    checked = sum(check_covariance(n) for n in (9,15,45,105,225))
    print(f"{checked} covariance checks")
    for p,e in [(3,1),(3,4),(3,6),(5,1),(5,4),(7,3),(101,1)]:
        check_local(p,e)
