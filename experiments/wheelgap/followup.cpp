#include "wheel.cpp"
#include <iostream>

extern "C" uint64_t first_integer(const int *p,int k,uint32_t target,uint64_t bound) {
    const auto ch=characters(p,k);
    std::vector<int> residues{0};
    int step=1,i=0;
    for (;i<k && step<=100000000/p[i] && uint64_t(step)*p[i]<=bound; ++i) {
        std::vector<int> crt;
        for (int block=0;block<p[i];++block) {
            for (const int r:residues) {
                const int x=block*step+r;
                if (ch[i][x%p[i]]==int((target>>i)&1)) crt.push_back(x);
            }
        }
        step*=p[i];
        residues.swap(crt);
    }
    for (uint64_t base=0;base<=bound;base+=step) {
        // Starting strictly above zero works also when no CRT components fit.
        for (auto r:residues) {
            const uint64_t x=base+r;
            if (x>bound) return 0;
            if (!x) continue;
            int j=i;
            for (;j<k;++j)
                if (ch[j][x%p[j]]!=int((target>>j)&1)) break;
            if (j==k) return x;
        }
    }
    return 0;
}
extern "C" void first_rational(const int *p,int k,uint32_t target,int *a,int *b) {
    if (!target) { *a=1; *b=1; return; }
    const auto ch=characters(p,k);
    std::unordered_map<uint32_t,int> first;
    for (int n=1;n<(1<<24);++n) {
        uint32_t sig=0;
        int i=0;
        for (;i<k;++i) {
            const int s=ch[i][n%p[i]];
            if (s<0) break;
            sig |= uint32_t(s)<<i;
        }
        if (i!=k) continue;
        auto it=first.find(sig^target);
        if (it!=first.end()) { *a=n; *b=it->second; return; }
        first.emplace(sig,n);
    }
    *a=*b=0;
}
