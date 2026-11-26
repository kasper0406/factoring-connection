// Factor-free operations on complete denominators of the registered 80-bit bank.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>
using Word = uint64_t;

extern "C" void select_multipliers(const Word *den, int *order) {
    Word score[18]{};
    bool used[18]{};
    for (int round = 0; round < 12; ++round) {
        int best = round;
        if (round >= 6) {
            best = -1;
            for (int k = 6; k < 18; ++k)
                if (!used[k] && (best < 0 || score[k] > score[best])) best = k;
        }
        order[round] = best;
        used[best] = true;
        for (int k = 6; k < 18; ++k) {
            if (used[k]) continue;
            for (int j = 0; j < 16; ++j)
                score[k] = std::max(score[k], std::gcd(den[16*k],den[16*best+j]));
        }
    }
}

// Replace intersecting factors a,b with gcd(a,b), a/g, b/g.
// Deduplication leaves the multiplicative span intact. No primality assumption.
static void insert_atom(std::vector<Word>& atoms, Word value) {
    std::vector<Word> pending{value};
    while (!pending.empty()) {
        Word v = pending.back();
        pending.pop_back();
        if (v == 1) continue;
        bool handled = false;
        for (size_t i = 0; i < atoms.size(); ++i) {
            Word b = atoms[i], g = std::gcd(v,b);
            if (g == 1) continue;
            handled = true;
            if (g == b) {
                pending.push_back(v/g);
            } else {
                atoms[i] = atoms.back();
                atoms.pop_back();
                pending.push_back(b/g);
                pending.push_back(v/g);
                pending.push_back(g);
            }
            break;
        }
        if (!handled) atoms.push_back(v);
    }
}

extern "C" int square_columns(const Word *den, int n, Word *out, int words) {
    // A square coprime atom contributes no square class. Other coprime
    // atoms are independent, even when composite.
    const Word small[]{2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    std::vector<Word> rest(n), atoms;
    for (int i=0; i<n; ++i) {
        Word v = den[i];
        for (int j=0;j<25;++j) {
            while (v % small[j] == 0) {
                out[i*words] ^= Word(1) << j;
                v /= small[j];
            }
        }
        rest[i] = v;
        insert_atom(atoms,v);
    }
    int col = 25;
    for (Word b: atoms) {
        Word s = Word(std::sqrt((long double)b));
        while (s+1 <= b/(s+1)) ++s;
        while (s > b/s) --s;
        if (s*s == b) continue;
        if (col >= words*64) return -1;
        for (int i=0;i<n;++i) {
            while (rest[i] % b == 0) {
                out[i*words+col/64] ^= Word(1) << (col%64);
                rest[i] /= b;
            }
        }
        ++col;
    }
    return col;
}
