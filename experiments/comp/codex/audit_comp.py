"""Direct all-row check of both rank enclosures on two cells, without strip pruning."""
import csv
import itertools
import json
from flint import arb
from test_comp import DATA, HERE, SCALE, endpoint, geometry, rank

def audit():
    comp=[json.loads(s) for s in (HERE/'cells_comp.jsonl').read_text().splitlines()]
    selected={(c['bits'],c['m'],c['p']) for c in comp if c['bits']==60 and c['m']==c['p'] in (3,4)}
    assert len(selected)==2
    records={k:{} for k in selected}
    for r in csv.DictReader((DATA/'embeddings.csv').open()):
        k=tuple(int(r[c]) for c in ('bits','m','p'))
        if k in selected:
            records[k][r['side'],int(r['root'])]=r
    for l in (DATA/'cells.jsonl').read_text().splitlines():
        c=json.loads(l)
        k=tuple(c[j] for j in ('bits','m','p'))
        if k not in selected:
            continue
        sites,*_=geometry(c,records[k])
        m=c['m']
        rec=next(r for r in comp if tuple(r[j] for j in ('bits','m','p'))==k)
        for lab,w in [('half',arb(m)/2),('one',arb(m)),('two',arb(2*m))]:
            radii=[(endpoint(SCALE*(w-loss).exp()-2*m*arb(2).sqrt()),
                    endpoint(SCALE*(w-loss).exp()+2*m*arb(2).sqrt(),True)) for _,loss in sites]
            bases=[[],[]]
            for row in itertools.product((-1,0,1),repeat=2*m):
                if next((s for s in row if s),0)<=0:
                    continue
                dist=[sum(sum(x*y for x,y in zip(row,v[i::2]))**2 for i in (0,1)) for v,_ in sites]
                for j,basis in enumerate(bases):
                    if len(basis)<m and all(d<=r[j]**2 for d,r in zip(dist,radii)) and rank(basis+[row])>len(basis):
                        basis.append(row)
            assert [len(b) for b in bases]==[u[0] for u in rec[lab]]
            print(*k,lab,*[len(b) for b in bases],flush=True)
if __name__=='__main__':
    audit()
