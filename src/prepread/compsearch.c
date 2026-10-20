/* Exhaustive branch-and-bound over quadratic compositions f = (x^2+c_k) o ... o (x^2+c_1), c_i in CS,
 * maximising min over the given primes of #{x in F_p : f(x) = 0}. Exact integer arithmetic.
 * Usage: compsearch k p1 [p2 ...]      Output: one JSON object on stdout.
 * The fibre is propagated backwards from {0}: T_{i-1} = {x : x^2 + c_i in T_i}; preimages of
 * distinct points are disjoint, so no deduplication is needed. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXP 8
#define MAXK 20
#define MAXARG 20
static const int CS[5] = {-2, -1, 0, 1, 2};
static int np_, K, P[MAXP];
static int *SQ[MAXP];               /* SQ[i][w] = a square root of w mod P[i], or -1 */
static int *BUF[MAXK + 1][MAXP];
static int LEN[MAXK + 1][MAXP];
static int seq[MAXK];
static long best = -1, nties = 0, nodes = 0;
static int argmax[MAXARG][MAXK], nargmax = 0;

static void dfs(int d) {
    nodes++;
    if (d == K) {
        long v = LEN[d][0];
        for (int i = 1; i < np_; i++) if (LEN[d][i] < v) v = LEN[d][i];
        if (v > best) { best = v; nties = 0; nargmax = 0; }
        if (v == best) {
            nties++;
            if (nargmax < MAXARG) { memcpy(argmax[nargmax], seq, sizeof(int) * K); nargmax++; }
        }
        return;
    }
    for (int ci = 0; ci < 5; ci++) {
        int c = CS[ci];
        long bound = -1;
        int dead = 0;
        for (int i = 0; i < np_; i++) {
            int p = P[i], n = 0;
            int *src = BUF[d][i], *dst = BUF[d + 1][i];
            for (int a = 0; a < LEN[d][i]; a++) {
                int w = ((src[a] - c) % p + p) % p;
                if (w == 0) dst[n++] = 0;
                else { int r = SQ[i][w]; if (r >= 0) { dst[n++] = r; dst[n++] = p - r; } }
            }
            LEN[d + 1][i] = n;
            if (n == 0) { dead = 1; }
            long b = (long)n << (K - d - 1);
            if (b > p) b = p;
            if (bound < 0 || b < bound) bound = b;
        }
        seq[d] = c;   /* seq[d] is the constant of layer K - d (applied d-th from the top) */
        if (dead && best >= 0 && 0 < best) continue;
        if (bound < best) continue;
        dfs(d + 1);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
    K = atoi(argv[1]); np_ = argc - 2;
    if (K > MAXK || np_ > MAXP) return 2;
    for (int i = 0; i < np_; i++) {
        int p = P[i] = atoi(argv[2 + i]);
        SQ[i] = malloc(sizeof(int) * p);
        for (int w = 0; w < p; w++) SQ[i][w] = -1;
        for (long x = 0; x <= (p - 1) / 2; x++) SQ[i][(x * x) % p] = (int)x;
        for (int d = 0; d <= K; d++) BUF[d][i] = malloc(sizeof(int) * (p + 2));
        BUF[0][i][0] = 0; LEN[0][i] = 1;
    }
    dfs(0);
    printf("{\"k\": %d, \"primes\": [", K);
    for (int i = 0; i < np_; i++) printf("%s%d", i ? ", " : "", P[i]);
    printf("], \"best\": %ld, \"nties\": %ld, \"nodes\": %ld, \"argmax\": [", best, nties, nodes);
    for (int a = 0; a < nargmax; a++) {
        printf("%s[", a ? ", " : "");
        /* print as (c_1, ..., c_k): c_1 is the innermost layer = seq[K-1] */
        for (int t = 0; t < K; t++) printf("%s%d", t ? ", " : "", argmax[a][K - 1 - t]);
        printf("]");
    }
    printf("]}\n");
    return 0;
}
