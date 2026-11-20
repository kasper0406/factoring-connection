// Exact closest pair from two half-boxes, subject to a nonzero linear functional.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
using I = int64_t;
using W = __int128_t;
using V = std::array<I,2>;
struct Point { V z; int code; };
struct Node { V lo, hi, tag; int begin, end, left=-1, right=-1; bool uniform; };
static W sq(I x) { return W(x)*x; }
static std::vector<Point> half(int digits, const I* data) {
    std::vector<Point> pts{{{0,0},0}};
    for (int i=0; i<digits; ++i) {
        std::vector<Point> out;
        out.reserve(3*pts.size());
        for (const auto &p:pts) for (int v=-1;v<=1;++v)
            out.push_back({{p.z[0]+v*data[2*i],p.z[1]+v*data[2*i+1]},3*p.code+v+1});
        pts.swap(out);
    }
    return pts;
}
struct Search {
    int l,r;
    std::vector<Point> a,b;
    std::vector<Node> tree;
    std::vector<Point> ta,tb;
    W best;
    int ia,ib;
    Search(int m,const I* data): l(2*(m/2)),r(2*m-l),a(half(l,data)),b(half(r,data+2*l)) { build(0,b.size()); }
    int build(int s,int e) {
        Node n{};
        n.begin=s; n.end=e; n.lo=n.hi=b[s].z;
        for(int j=s+1;j<e;++j) for(int k=0;k<2;++k) {
            n.lo[k]=std::min(n.lo[k],b[j].z[k]);
            n.hi[k]=std::max(n.hi[k],b[j].z[k]);
        }
        int id=tree.size();
        tree.push_back(n);
        if(e-s>8) {
            int k=(n.hi[1]-n.lo[1]>n.hi[0]-n.lo[0]);
            int mid=(s+e)/2;
            std::nth_element(b.begin()+s,b.begin()+mid,b.begin()+e,
                [k](const Point& u,const Point& v){return u.z[k]<v.z[k];});
            int left=build(s,mid),right=build(mid,e);
            tree[id].left=left;tree[id].right=right;
        }
        return id;
    }
    void label(const I* data) {
        ta=half(l,data);tb=half(r,data+2*l);
        for(auto it=tree.rbegin();it!=tree.rend();++it) {
            Node& n=*it;
            if(n.left>=0) {
                auto &x=tree[n.left],&y=tree[n.right];
                n.tag=x.tag;
                n.uniform=x.uniform && y.uniform && x.tag==y.tag;
            } else {
                n.tag=tb[b[n.begin].code].z;n.uniform=true;
                for(int j=n.begin+1;j<n.end;++j) if(tb[b[j].code].z!=n.tag){n.uniform=false;break;}
            }
        }
    }
    static W bound(const Node& n,V z) {
        W d=0;
        for(int k=0;k<2;++k) {
            if(z[k]<n.lo[k]) d+=sq(n.lo[k]-z[k]);
            if(z[k]>n.hi[k]) d+=sq(z[k]-n.hi[k]);
        }
        return d;
    }
    void visit(int i,const Point& p,V z,V tag) {
        const auto &n=tree[i];
        if(n.uniform && n.tag==tag) return;
        if(n.left>=0) {
            int s=n.left,t=n.right;
            W ds=bound(tree[s],z),dt=bound(tree[t],z);
            if(dt<ds){std::swap(ds,dt);std::swap(s,t);}
            if(ds<best) visit(s,p,z,tag);
            if(dt<best) visit(t,p,z,tag);
        } else for(int j=n.begin;j<n.end;++j) {
            const auto &q=b[j];
            W d=sq(q.z[0]-z[0])+sq(q.z[1]-z[1]);
            if(d<best && tb[q.code].z!=tag) {best=d;ia=p.code;ib=q.code;}
        }
    }
    void solve(const I* weights,int* digits) {
        label(weights);
        best=W(1)<<126;
        for(const auto &p:a) {
            V z={-p.z[0],-p.z[1]},tag={-ta[p.code].z[0],-ta[p.code].z[1]};
            if(bound(tree[0],z)<best) visit(0,p,z,tag);
        }
        assert(best<(W(1)<<126));
        for(int j=l-1;j>=0;--j){digits[j]=ia%3-1;ia/=3;}
        for(int j=l+r-1;j>=l;--j){digits[j]=ib%3-1;ib/=3;}
    }
};
extern "C" void* create(int m,const I* data){return new Search(m,data);}
extern "C" void solve(void* p,const I* w,int* out){static_cast<Search*>(p)->solve(w,out);}
extern "C" void destroy(void* p){delete static_cast<Search*>(p);}
