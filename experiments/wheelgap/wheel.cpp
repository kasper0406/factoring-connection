#include <algorithm>
#include <cstdint>
#include <ctime>
#include <unordered_map>
#include <vector>

struct Result { int a, b, prefix; };
using Table = std::vector<std::vector<signed char>>;
static Table characters(const int *p, int k) {
    Table table(k);
    for (int i = 0; i < k; ++i) {
        table[i].assign(p[i], 1);
        table[i][0] = -1;
        for (int x = 1; x <= p[i]/2; ++x)
            table[i][int64_t(x)*x % p[i]] = 0;
    }
    return table;
}

extern "C" void exact(const int *p, int k, uint32_t target, int rational, Result *out) {
    const auto start = std::clock();
    const auto expired = [&]() {return double(std::clock()-start)/CLOCKS_PER_SEC >= .25;};
    const auto ch = characters(p,k);
    *out = {0,0,0};
    if (rational) {
        std::unordered_map<uint32_t,int> first;
        for (int n=1; n <= (1<<19); ++n) {
            if (!(n % 1024) && expired()) break;
            uint32_t sig = 0;
            bool unit = true;
            for (int i=0; i<k; ++i) {
                const int bit=ch[i][n % p[i]];
                if (bit < 0) {unit=false; break;}
                sig |= uint32_t(bit)<<i;
            }
            out->prefix=n;
            if (!unit) continue;
            if (!target) { *out={1,1,n}; return; }
            auto it=first.find(sig^target);
            if (it!=first.end()) { *out={n,it->second,n}; return; }
            first.emplace(sig,n);
        }
    } else {
        int step=1, split=0;
        while (split<k && step <= 30000/p[split]) step*=p[split++];
        std::vector<int> residues;
        const auto matches = [&](int n,int lo,int hi) {
            for (int i=lo;i<hi;++i)
                if (ch[i][n%p[i]] != int((target>>i)&1)) return false;
            return true;
        };
        for (int r=1;r<=step;++r) if (matches(r,0,split)) residues.push_back(r);
        for (int base=0;base<(1<<22);base+=step) {
            if (expired()) break;
            for (const int r : residues) {
                int n=base+r;
                if (n>(1<<22)) break;
                if (matches(n,split,k)) { *out={n,1,n}; return; }
            }
            out->prefix=std::min(base+step,1<<22);
        }
    }
}
