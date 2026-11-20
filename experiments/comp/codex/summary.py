import json
from pathlib import Path

here=Path(__file__).resolve().parent
def lines(name):
    return [json.loads(s) for s in (here/name).read_text().splitlines()]
cells=lines('cells_comp.jsonl')
pencils=lines('pencils.jsonl')
assert len(cells)==195 and len(pencils)==390
assert all(r[w][0]==r[w][1] for r in cells for w in ('half','one','two'))

def certified(r,w):
    rank,one,t=r[w][0]
    return rank==r['m'] or (2*rank>=r['m'] and one and t)

def report(name,s):
    print(name, len(s),sum(r['Q_at_quarter'] for r in s),sum(r['Q_nonpos'] for r in s),
          *[sum(certified(r,w) for r in s) for w in ('half','one','two')],
          *[sum(r[w][0][0]==r['m'] for r in s) for w in ('half','one','two')])
print('group cells Q>=X/4 Q<=0 rank_cert[w=m/2,m,2m] full_rank[w=m/2,m,2m]')
report('all',cells)
report('loss premise',[r for r in cells if r['prem_count']])
report('H<=4X',[r for r in cells if r['prem_height']])
report('both premise cuts',[r for r in cells if r['prem_count'] and r['prem_height']])
for m in range(3,9):
    report('m='+str(m),[r for r in cells if r['m']==m])
for k in ('H_over_X','Q_over_X'):
    print(k,'outer range in billionths',min(r[k][0] for r in cells),max(r[k][1] for r in cells))
print('Cells with an affine unit',sum(r['affine_min']==1 for r in cells))
print('Q<=0 plus eligible affine witness below X and capped',
      sum(r['Q_nonpos'] and r['affine_below_X'] and r['affine_capped'] for r in cells))
print('Most favorable affine norm',min(r['affine_min'] for r in cells))
print('Pencil rows: coverage, leading height >=X, bound >=X, leading passed but bound failed, mixed factors present')
for name,rows in [('eligible',[r for r in pencils if r['unit']]),('ineligible',[r for r in pencils if not r['unit']])]:
    print(name,len(rows),sum(r['lead_X'] for r in rows),sum(r['bound_X'] for r in rows),
          sum(r['lead_X'] and not r['bound_X'] for r in rows),sum(r['mixed_factors']>0 for r in rows))
bycell={}
for r in pencils:
    bycell.setdefault((r['bits'],r['m'],r['p']),[]).append(r)
assert all(len(rs)==2 for rs in bycell.values())
print('Pencil cells passing, empty',sum(any(r['unit'] for r in rs) and all(r['bound_X'] for r in rs if r['unit']) for rs in bycell.values()),sum(not any(r['unit'] for r in rs) for rs in bycell.values()))
print('Gamma upper bound (billionths)',max(r['Gamma'][1] for r in pencils))
