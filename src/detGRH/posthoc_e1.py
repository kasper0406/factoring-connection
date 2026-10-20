# POST HOC diagnostics for detGRH E1 (not pre-registered). Run: PYTHONPATH=. .venv/bin/python src/detGRH/posthoc_e1.py
import json
from src.detGRH import core, diag
cfg=json.load(open('experiments/configs/detGRH_e1.json'))
tot=own=own_triv=other=other_nd=0; signodd=0; nrows=0
for cls in ['blum','mixed','both1','tri']:
    for seed in range(4):
        bits=32
        N,fs=diag.make_instance(cls,bits,seed,min_factor=core.fb_bound(1<<bits))
        fb=core.FactorBase.make(N)
        rel=core.collect(N,fb,K=cfg['K'],X=fb.B,T=cfg['T'],surplus=cfg['surplus'])
        nrows+=len(rel.kernel); signodd+=sum(v[0]%2 for v in rel.kernel)
        basis=core.F2Basis(rel.kernel,fb.n)
        # index of each seed's kernel row
        idx=[i for i,s in enumerate(rel.kernel_src) if s=='seed']
        for k,(r,R,eR) in enumerate(rel.squares):
            i=idx[k] if k<len(idx) else None
            combo=basis.express(core.bits_of(eR))
            if combo is None: continue
            y=core.canonical_root(fb,basis,eR); chi=y*pow(r,-1,N)%N
            tot+=1
            if combo==(1<<i):
                own+=1; own_triv+= chi in (1,N-1)
            else:
                other+=1; other_nd+= chi not in (1,N-1)
print('seeds',tot,'canonical combo = own relation only:',own,'(chi trivial in',own_triv,')','other:',other,'nondiag among other:',other_nd, round(other_nd/max(other,1),3))
print('kernel rows with odd -1 coordinate:',signodd,'of',nrows)
