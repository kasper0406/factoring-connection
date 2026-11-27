"""E5 (own idea, HARNESS): the whole mod-2 level-1 family at N is one polynomial P_N(X) = sum_k a(N,Delta^k) X^k
= J_N(X) mod 2 (Faber / Hecke correspondence T_N on the j-line).  Evaluate at the ordinary point j=1 (E: y^2+xy=x^3+1
over F_2, CM by O_K, K=Q(sqrt-7)) and test the predicted CM-trace identity
   P_N(1) = Tr H_{-7N^2} + (1+chi(p)) Tr H_{-7q^2} + (1+chi(q)) Tr H_{-7p^2} + (1+chi(p))(1+chi(q))  (mod 2),
which for p,q unramified in K collapses to P_N(1) = Tr H_{-7N^2} (mod 2)."""
import json
import cypari2
from sympy import primerange, factorint
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
M = 301; MASK = (1 << (M + 1)) - 1
Dl = 0
r, m = 0, 1
while m * m <= M: r ^= 1 << (m * m); m += 2
Delta = r
def mul(A, B):
    out = 0
    while B:
        lb = B & -B; out ^= A << (lb.bit_length() - 1); B ^= lb
    return out & MASK
# P_n(1) = [q^n] sum_{k>=1} Delta^k
S, Pw = 0, 1
for k in range(1, M + 1):
    Pw = mul(Pw, Delta); S ^= Pw
    if Pw == 0: break
def tr(D):
    H = pari.polclass(D); deg = int(pari.poldegree(H)); return -int(pari.polcoef(H, deg - 1))
chi = lambda l: int(pari.kronecker(-7, l))
res = []
for n in range(3, M + 1, 2):
    f = factorint(n)
    if any(l == 7 for l in f) or any(e > 1 for e in f.values()) or len(f) > 2: continue
    P1 = (S >> n) & 1
    if len(f) == 1:
        l = n; pred = (tr(-7 * l * l) + (1 + chi(l))) % 2
    else:
        p, q = sorted(f)
        pred = (tr(-7 * n * n) + (1 + chi(p)) * tr(-7 * q * q) + (1 + chi(q)) * tr(-7 * p * p)
                + (1 + chi(p)) * (1 + chi(q))) % 2
    res.append(dict(n=n, P_n_at_1=P1, predicted=pred, collapsed=tr(-7 * n * n) % 2))
ok = sum(r["P_n_at_1"] == r["predicted"] for r in res)
ok2 = sum(r["P_n_at_1"] == r["collapsed"] for r in res if len(factorint(r["n"])) == 2)
print("tested", len(res), "identity holds", ok, "| semiprimes: collapsed form holds", ok2, "of",
      sum(len(factorint(r['n'])) == 2 for r in res), "| P=1 count", sum(r["P_n_at_1"] for r in res))
json.dump(dict(tested=len(res), holds=ok, rows=res), open("/Volumes/git/fc-coeffam/experiments/coeffam/e5_faber_cm.json", "w"))
