"""Factor-blind exact region counts under a positive definite binary form, by integer-hull
(Stern-Brocot / continued-fraction) descent with Pick-type edge bookkeeping.

count_box(b, c, T, X) = #{(x,y): 0<=x,y<=X, x^2+bxy+cy^2 <= T},  b,c >= 0.
The accumulator can be reduced mod 2^k (mod_bits); geometry must stay exact.
Each hull vertex costs O(log T) isqrt-based line tests. Stats count vertices and line tests.
"""
from math import isqrt


class Region:
    """K' = {(x,y): L<=x<=R, y<0 or x^2+bxy+cy^2<=T}: convex and down-closed for x,y>=0."""

    def __init__(self, b, c, T, L, R):
        self.b, self.c, self.T, self.L, self.R = b, c, T, L, R
        self.tests = 0

    def F(self, x, y):
        return x*x + self.b*x*y + self.c*y*y

    def inside(self, x, y):
        return self.L <= x <= self.R and (y < 0 or self.F(x, y) <= self.T)

    def top(self, x):
        """max y>=0 with F(x,y)<=T (assumes F(x,0)<=T)."""
        b, c, T = self.b, self.c, self.T
        disc = b*b*x*x - 4*c*(x*x - T)
        y = (isqrt(disc) - b*x)//(2*c)
        while self.F(x, y+1) <= T:
            y += 1
        while y > 0 and self.F(x, y) > T:
            y -= 1
        return y

    def first_on_line(self, px, py, dx, dy, tmin):
        """smallest integer t>=tmin with (px+t*dx, py-t*dy) in K', dx>=1; None if none.
        K' convex => the feasible t form an interval; a quadratic gives it exactly."""
        self.tests += 1
        # x-range
        tmax = (self.R - px)//dx
        tlo = max(tmin, -((px - self.L)//dx))  # ceil((L-px)/dx)
        if tlo > tmax:
            return None
        # quadratic q(t) = F(px+t dx, py - t dy) - T = A t^2 + B t + C
        b, c = self.b, self.c
        ex, ey = dx, -dy
        A = ex*ex + b*ex*ey + c*ey*ey
        B = 2*px*ex + b*(px*ey + py*ex) + 2*c*py*ey
        C = self.F(px, py) - self.T

        def ok(t):
            return self.inside(px + t*dx, py - t*dy)

        if ok(tlo):
            return tlo
        # Feasible set is an interval; tlo infeasible -> it lies to the right.
        # Candidates: the left root of q, or the first t with y<0.
        cands = []
        D = B*B - 4*A*C
        if D >= 0:
            s = isqrt(D)
            t = (-B - s)//(2*A)  # near left root
            for u in (t - 1, t, t + 1, t + 2):
                if u >= tlo:
                    cands.append(u)
        if dy > 0:
            t = py//dy + 1  # first t with py - t*dy < 0
            cands.append(max(t, tlo))
        best = None
        for u in sorted(set(cands)):
            if u <= tmax and ok(u):
                best = u
                break
        return best

    def last_on_line(self, px, py, dx, dy):
        """largest t>=0 with point in K' given t=0 is in K'."""
        self.tests += 1
        tmax = (self.R - px)//dx
        b, c = self.b, self.c
        ex, ey = dx, -dy
        A = ex*ex + b*ex*ey + c*ey*ey
        B = 2*px*ex + b*(px*ey + py*ex) + 2*c*py*ey
        C = self.F(px, py) - self.T
        D = B*B - 4*A*C
        s = isqrt(max(D, 0))
        t = (-B + s)//(2*A) + 2
        t = min(t, tmax)
        # beyond y<0 everything (within x range) is inside
        if dy > 0 and py - tmax*dy < 0:
            return tmax
        while t > 0 and not self.inside(px + t*dx, py - t*dy):
            t -= 1
        return t


def column_sum(b, c, T, L, R, mod_bits=None, stats=None):
    """sum_{x=L}^{R} Y(x), Y(x)=max{y>=0: F(x,y)<=T}; requires F(R,0)<=T and L<=R."""
    K = Region(b, c, T, L, R)
    M = (1 << mod_bits) if mod_bits else None
    x, y = L, K.top(L)
    total = 0
    verts = 0
    while x < R:
        verts += 1
        b0 = y - K.top(x + 1)
        lo, hi = (1, b0), (1, b0 - 1)        # (dx, down)
        while True:
            j = K.first_on_line(x + hi[0], y - hi[1], lo[0], lo[1], 1)
            if j is None:
                break
            lo, hi = (j*lo[0] + hi[0], j*lo[1] + hi[1]), ((j-1)*lo[0] + hi[0], (j-1)*lo[1] + hi[1])
        a, d = lo
        k = K.last_on_line(x, y, a, d)
        assert k >= 1
        total += k*a*y - (a*d*k*(k-1)//2 + k*(a-1)*(d+1)//2)
        if M:
            total %= M
        x, y = x + k*a, y - k*d
    total += y
    if stats is not None:
        stats['vertices'] = stats.get('vertices', 0) + verts
        stats['tests'] = stats.get('tests', 0) + K.tests
    return total % M if M else total


def count_box(b, c, T, X, mod_bits=None, stats=None):
    """#{0<=x,y<=X : x^2+bxy+cy^2<=T}."""
    if T < 0:
        return 0
    xmax = isqrt(T)               # F(x,0)=x^2<=T
    R = min(X, xmax)
    # x1 = first x with top(x) <= X, i.e. F(x, X+1) > T
    # F(x,X+1) = x^2 + b(X+1)x + c(X+1)^2 increasing in x>=0
    Y1 = X + 1
    C0 = c*Y1*Y1 - T
    if C0 > 0:
        x1 = 0
    else:
        disc = b*b*Y1*Y1 - 4*C0
        x1 = (isqrt(disc) - b*Y1)//2
        while x1*x1 + b*Y1*x1 + C0 <= 0:
            x1 += 1
        while x1 > 0 and (x1-1)**2 + b*Y1*(x1-1) + C0 > 0:
            x1 -= 1
    x1 = min(x1, R + 1)
    total = x1*(X + 1)
    if x1 <= R:
        total += column_sum(b, c, T, x1, R, mod_bits, stats) + (R - x1 + 1)
    if mod_bits:
        total %= 1 << mod_bits
    return total


def brute_box(b, c, T, X):
    n = 0
    for x in range(X + 1):
        for y in range(X + 1):
            if x*x + b*x*y + c*y*y <= T:
                n += 1
    return n


def short_square(N, b, c, mod_bits=3, stats=None):
    """floor8 query S(N) = #{0<=x,y<X : N | x^2+bxy+cy^2}, X = isqrt(N)+1, mod 2^mod_bits.
    Box is [0, X-1]^2; values F<=maxF; S = 1 + sum_j (C(jN) - C(jN-1))."""
    Xb = isqrt(N)                     # inclusive bound X-1
    maxF = (1 + b + c)*Xb*Xb
    J = maxF//N
    M = 1 << mod_bits
    s = 1
    for j in range(1, J + 1):
        s += count_box(b, c, j*N, Xb, mod_bits + 1, stats) - count_box(b, c, j*N - 1, Xb, mod_bits + 1, stats)
    return s % M
