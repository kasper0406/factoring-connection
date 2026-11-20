// Enumerate the simultaneous inner and outer boxes. A half-sum strip at the
// tightest place is a necessary condition for the outer box.
#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <vector>
#include <cassert>

using I = int64_t;
using U = uint64_t;
using Wide = __int128_t;
constexpr U prime = (U(1)<<61)-1, root2=U(1)<<31;
using Point = std::array<I,2>;
using Digits = std::array<int,16>;
using Row = std::array<U,8>;

static U multiply(U a,U b) {
    __uint128_t v=__uint128_t(a)*b;
    U x=U(v>>61)+(U(v)&prime);
    return x>=prime?x-prime:x;
}
static U inverse(U a) {
    U v=1;
    for(U e=prime-2;e;e>>=1,a=multiply(a,a))
        if(e&1) v=multiply(v,a);
    return v;
}
struct Basis {
    int m, rank=0;
    Row pivots[8]{};
    int* output;
    void add(const Digits& digits) {
        if(rank==m) return;
        Row row{};
        for(int j=0;j<m;++j) {
            I t=digits[2*j]+I(root2)*digits[2*j+1];
            row[j]=t<0?I(prime)+t:t;
        }
        for(int j=0;j<m;++j) {
            U a=row[j];
            if(!a) continue;
            if(!pivots[j][j]) {
                U scale=inverse(a);
                for(int k=j;k<m;++k) pivots[j][k]=multiply(scale,row[k]);
                std::copy_n(digits.begin(),2*m,output+rank*2*m);
                ++rank;
                return;
            }
            for(int k=j+1;k<m;++k) {
                U v=multiply(a,pivots[j][k]);
                row[k]=(row[k]>=v?row[k]:row[k]+prime)-v;
            }
        }
    }
};
static std::vector<Digits> make_digits(int m) {
    int count=1;
    for(int j=0;j<m;++j) count*=3;
    std::vector<Digits> out(count);
    for(int i=0;i<count;++i) {
        int t=i;
        for(int j=m-1;j>=0;--j) {out[i][j]=t%3-1;t/=3;}
    }
    return out;
}
struct Place {
    I inner,outer;
    Wide inner2,outer2;
    std::vector<Point> left,right;
};
static std::vector<Point> sums(int m,const I* weights,const std::vector<Digits>& ds) {
    std::vector<Point> result(ds.size());
    for(size_t i=0;i<ds.size();++i) for(int j=0;j<m;++j)
        for(int k=0;k<2;++k) result[i][k]+=ds[i][j]*weights[2*j+k];
    return result;
}
extern "C" void ranks(int m,int count,const I* data,const I* bounds,
                       int* lower,int* upper,int* result) {
    assert(m>=3 && m<=8 && count>0);
    const auto ds=make_digits(m);
    Basis inside{m,0,{},lower}, outside{m,0,{},upper};
    std::vector<Place> places;
    for(int i=0;i<count;++i) {
        const I* row=data+i*4*m;
        I inner=bounds[2*i],outer=bounds[2*i+1];
        I cap=0;
        for(int j=0;j<4*m;++j) cap+=std::abs(row[j]);
        if(inner>=cap) continue;
        places.push_back({inner,outer,Wide(inner)*inner,Wide(outer)*outer,
                          sums(m,row,ds),sums(m,row+2*m,ds)});
    }
    auto accept=[&](int a,int b) {
        bool certain=true;
        for(const auto& p:places) {
            I x=p.left[a][0]+p.right[b][0],y=p.left[a][1]+p.right[b][1];
            Wide dist=Wide(x)*x+Wide(y)*y;
            if(dist>p.outer2) return;
            if(dist>p.inner2) certain=false;
        }
        Digits digits=ds[a];
        std::copy_n(ds[b].begin(),m,digits.begin()+m);
        outside.add(digits);
        if(certain) inside.add(digits);
    };
    int zero=(ds.size()-1)/2;
    // Try the coordinate lines first, including small scalar letters.
    int powers[16];powers[2*m-1]=1;
    for(int j=2*m-1;j>0;--j) powers[j-1]=powers[j]*3;
    for(int j=0;j<m;++j) for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b) {
        int code=zero*int(ds.size())+zero+a*powers[2*j]+b*powers[2*j+1];
        accept(code/int(ds.size()),code%int(ds.size()));
    }
    if(inside.rank<m) {
        assert(!places.empty());
        std::sort(places.begin(),places.end(),[](const Place& a,const Place& b){return a.outer<b.outer;});
        const auto& p=places[0];
        std::vector<int> order(ds.size());
        std::iota(order.begin(),order.end(),0);
        I ranges[2]={};
        for(auto z:p.right) for(int k=0;k<2;++k) ranges[k]=std::max(ranges[k],std::abs(z[k]));
        int axis=ranges[1]>ranges[0];
        std::sort(order.begin(),order.end(),[&](int a,int b){return p.right[a][axis]<p.right[b][axis];});
        for(int a=0;a<=zero && inside.rank<m;++a) {
            I target=-p.left[a][axis];
            auto first=std::lower_bound(order.begin(),order.end(),target-p.outer,
                [&](int b,I value){return p.right[b][axis]<value;});
            for(auto it=first;it!=order.end() && p.right[*it][axis]<=target+p.outer && inside.rank<m;++it)
                accept(a,*it);
        }
    }
    result[0]=inside.rank;result[1]=outside.rank;
}
