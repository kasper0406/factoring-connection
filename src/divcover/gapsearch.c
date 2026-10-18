/*
 * divcover exact GAP-cover search (diagnostic/combinatorial; no factoring of N).
 *
 * Finds the minimal height H*(n, shape) of a positive rank-r generalized
 * arithmetic progression
 *     Q = { b + sum_k x_k c_k : 0 <= x_k < L_k },   b >= 1, c_k >= 1,
 * (height = max element = b + sum_k (L_k - 1) c_k) such that every target d
 * divides some element of Q.  Targets: mode 0 = all primes p <= n,
 * mode 1 = all d in [2, n] (full n-divisor property).
 * Every GAP has a representation with all c_k >= 0 (flip x_k -> L_k-1-x_k);
 * c_k = 0 lowers the rank, so c_k >= 1 loses nothing for "rank <= r".
 *
 * Usage: gapsearch n mode Hmax L1 [L2 [L3]]
 * Output (one line, key=value): best height (0 if none <= Hmax), b, c..., and
 * the number of parameter vectors examined.
 *
 * Exhaustive: every (b, c) with height < best is either examined or excluded
 * by an exact residue test, so the reported best is the exact minimum
 * (or "none" if no GAP of height <= Hmax exists).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXT 4096
static int n, mode, r;
static long Hmax;
static int L[3];
static int nt;            /* number of targets */
static int tgt[MAXT];     /* targets, most restrictive first (descending) */
static unsigned char *mask[MAXT]; /* mask[t][res] = 1 iff res = -<x,c> mod d for some x */
static int maskok[MAXT];
static long cc[3];
static long nparams = 0;

static int is_prime(int x) {
    if (x < 2) return 0;
    for (int q = 2; q * q <= x; q++) if (x % q == 0) return 0;
    return 1;
}

static void build_mask(int t) {
    int d = tgt[t];
    unsigned char *m = mask[t];
    memset(m, 0, d);
    long c0 = cc[0] % d, c1 = r > 1 ? cc[1] % d : 0, c2 = r > 2 ? cc[2] % d : 0;
    int l1 = r > 1 ? L[1] : 1, l2 = r > 2 ? L[2] : 1;
    for (int z = 0; z < l2; z++) {
        long s2 = (z * c2) % d;
        for (int y = 0; y < l1; y++) {
            long s1 = (s2 + y * c1) % d;
            long s = s1;
            for (int x = 0; x < L[0]; x++) {
                m[(d - s) % d] = 1;
                s += c0; if (s >= d) s -= d;
            }
        }
    }
    maskok[t] = 1;
}

static inline int good(int t, long b) {
    if (!maskok[t]) build_mask(t);
    return mask[t][b % tgt[t]];
}

/* For current c: smallest b in [1, bmax] covering all targets, or 0. */
static long *reslist = NULL;
static long minimal_b(long bmax) {
    if (bmax < 1) return 0;
    for (int t = 0; t < nt; t++) maskok[t] = 0;
    /* CRT pre-filter on the two most restrictive coprime targets */
    int t0 = 0, t1 = -1;
    for (int t = 1; t < nt; t++) {
        int a = tgt[t0], bb = tgt[t], g = a, h = bb;
        while (h) { int tmp = g % h; g = h; h = tmp; }
        if (g == 1) { t1 = t; break; }
    }
    build_mask(t0);
    long M = tgt[t0];
    long nres = 0;
    if (t1 >= 0) {
        build_mask(t1);
        M = (long)tgt[t0] * tgt[t1];
        for (long x = 0; x < M; x++)
            if (mask[t0][x % tgt[t0]] && mask[t1][x % tgt[t1]]) reslist[nres++] = x;
    } else {
        for (long x = 0; x < M; x++) if (mask[t0][x]) reslist[nres++] = x;
    }
    if (nres == 0) return 0;
    for (long base = 0; base <= bmax; base += M) {
        for (long i = 0; i < nres; i++) {
            long b = base + reslist[i];
            if (b < 1) continue;
            if (b > bmax) return 0;
            int ok = 1;
            for (int t = 0; t < nt; t++) {
                if (t == t0 || t == t1) continue;
                if (!good(t, b)) { ok = 0; break; }
            }
            if (ok) return b;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: gapsearch n mode Hmax L1 [L2 [L3]]\n"); return 2; }
    n = atoi(argv[1]); mode = atoi(argv[2]); Hmax = atol(argv[3]);
    r = argc - 4; if (r > 3) r = 3;
    for (int k = 0; k < r; k++) L[k] = atoi(argv[4 + k]);
    for (int k = 0; k < r; k++) if (L[k] < 2) { fprintf(stderr, "L_k >= 2 required\n"); return 2; }
    /* targets descending; primes first in prime mode, all d in full mode */
    nt = 0;
    for (int d = n; d >= 2; d--) {
        if (mode == 0 && !is_prime(d)) continue;
        if (nt >= MAXT) { fprintf(stderr, "too many targets\n"); return 2; }
        tgt[nt++] = d;
    }
    /* full mode: put primes before composites within descending order? keep
       descending but move large primes first: sort key = (is_prime desc, d desc) */
    if (mode == 1) {
        int tmp[MAXT], k = 0;
        for (int t = 0; t < nt; t++) if (is_prime(tgt[t])) tmp[k++] = tgt[t];
        for (int t = 0; t < nt; t++) if (!is_prime(tgt[t])) tmp[k++] = tgt[t];
        memcpy(tgt, tmp, sizeof(int) * nt);
    }
    for (int t = 0; t < nt; t++) mask[t] = malloc(tgt[t] + 1);
    reslist = malloc(sizeof(long) * (long)(n + 1) * (n + 1));
    long best = Hmax + 1, bb = 0, bc[3] = {0, 0, 0};
    long w0 = L[0] - 1, w1 = r > 1 ? L[1] - 1 : 0, w2 = r > 2 ? L[2] - 1 : 0;
    /* symmetry: if two sides are equal we may order their generators */
    if (r == 1) {
        for (cc[0] = 1; w0 * cc[0] + 1 < best; cc[0]++) {
            nparams++;
            long b = minimal_b(best - 1 - w0 * cc[0]);
            if (b) { best = b + w0 * cc[0]; bb = b; bc[0] = cc[0]; }
        }
    } else if (r == 2) {
        for (cc[0] = 1; w0 * cc[0] + w1 + 1 < best; cc[0]++) {
            for (cc[1] = 1; w0 * cc[0] + w1 * cc[1] + 1 < best; cc[1]++) {
                if (L[0] == L[1] && cc[1] < cc[0]) continue;
                nparams++;
                long b = minimal_b(best - 1 - w0 * cc[0] - w1 * cc[1]);
                if (b) { best = b + w0 * cc[0] + w1 * cc[1]; bb = b; bc[0] = cc[0]; bc[1] = cc[1]; }
            }
        }
    } else {
        for (cc[0] = 1; w0 * cc[0] + w1 + w2 + 1 < best; cc[0]++)
          for (cc[1] = 1; w0 * cc[0] + w1 * cc[1] + w2 + 1 < best; cc[1]++) {
            if (L[0] == L[1] && cc[1] < cc[0]) continue;
            for (cc[2] = 1; w0 * cc[0] + w1 * cc[1] + w2 * cc[2] + 1 < best; cc[2]++) {
                if (L[1] == L[2] && cc[2] < cc[1]) continue;
                nparams++;
                long b = minimal_b(best - 1 - w0 * cc[0] - w1 * cc[1] - w2 * cc[2]);
                if (b) { best = b + w0 * cc[0] + w1 * cc[1] + w2 * cc[2]; bb = b;
                         bc[0] = cc[0]; bc[1] = cc[1]; bc[2] = cc[2]; }
            }
          }
    }
    if (best > Hmax) printf("best=0 b=0 c=0,0,0 nparams=%ld\n", nparams);
    else printf("best=%ld b=%ld c=%ld,%ld,%ld nparams=%ld\n", best, bb, bc[0], bc[1], bc[2], nparams);
    return 0;
}
