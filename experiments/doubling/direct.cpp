#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

struct Sample {
    uint64_t modulus, index, residue;
    void include(uint64_t prime) {
        if (residue == 0 || prime > 2 * index) return;
        for (uint64_t power = prime; power <= 2 * index; power *= prime) {
            uint64_t low = index % power;
            if (2 * low >= power) residue = residue * prime % modulus;
            if (power > 2 * index / prime) break;
        }
    }
};

int main(int argc, char **argv) {
    std::vector<Sample> samples;
    uint64_t last = 0;
    for (int i = 1; i < argc; ++i) {
        uint64_t n = std::strtoull(argv[i], nullptr, 10);
        if (n < 5 || n > UINT32_MAX || n % 4 != 1) return 1;
        samples.push_back({n, n / 4, 1});
        last = std::max(last, 2 * (n / 4));
    }
    std::vector<uint64_t> base;
    for (uint64_t j = 3; j * j <= last; j += 2) {
        bool prime = true;
        for (uint64_t p : base) {
            if (p * p > j) break;
            if (j % p == 0) { prime = false; break; }
        }
        if (prime) base.push_back(j);
    }
    for (auto &s : samples) s.include(2);
    constexpr uint64_t size = 1 << 18;
    for (uint64_t start = 3; start <= last; start += 2 * size) {
        uint64_t end = std::min(last + 1, start + 2 * size);
        std::vector<uint8_t> marked((end - start + 1) / 2);
        for (uint64_t p : base) {
            uint64_t k = std::max(p * p, ((start + p - 1) / p) * p);
            if (k % 2 == 0) k += p;
            for (; k < end; k += 2 * p) marked[(k - start) / 2] = 1;
        }
        for (uint64_t i = 0; i < marked.size(); ++i) {
            if (!marked[i]) for (auto &s : samples) s.include(start + 2 * i);
        }
    }
    for (const auto &s : samples) std::cout << s.modulus << ' ' << s.residue << '\n';
}
