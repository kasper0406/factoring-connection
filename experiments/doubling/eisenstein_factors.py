import json
from flint import fmpz
from eisenstein_search import hurwitz,bern_num
from experiment import HERE

if __name__=='__main__':
    G=hurwitz(1200)
    data=json.loads((HERE/'eisenstein_search.json').read_text())
    results=[]
    for p,hs in data['eligible']:
        residual=abs(int(G[p-1].p))
        for s in range(2,10000):
            while residual%s==0:
                residual//=s
        factors=fmpz(residual).factor_smooth(32)
        for q,e in factors:
            if q>p and q%4==1 and (q-1)%(p-1) in hs and q.is_prime() and bern_num(p-1)%int(q):
                results.append(dict(p=p,q=str(q)))
    (HERE/'eisenstein_factors.json').write_text(json.dumps(results,indent=2)+'\n')
    print('Additional hits',results)
