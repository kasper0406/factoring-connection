"""Descriptive full/uncertain results for the ordinary box and q lattices separately."""
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from threadpoolctl import threadpool_limits
from rbase import YS
from rfit import HEADS,NAMES
from rstats import data,block_arrays,hierarchical,loss,write


def render():
    rows=[]
    for l,ycut in enumerate(YS):
        path=Path(f'results/A_{l}'); d=data(path/'test_data.npz')
        for fam in NAMES:
            z=data(path/f'test_{fam}.npz')
            for h,head in enumerate(HEADS):
                for geometry,gm in [('U_box',d['block']>=32),('Q_lattice',d['block']<32)]:
                    for label,sel in [('full',z['domain']&gm),('uncertain',z['domain']&gm&d['ok'][:,h])]:
                        if not sel.any():continue
                        ys=d['y'][sel,h];n=d['n'][sel];b=d['block'][sel];pb=d['p'][sel,h];pm=z['pred'][sel,h]
                        base=loss(ys,pb,h<3);aug=loss(ys,pm,h<3);gain=base-aug
                        av,bt,_=hierarchical(block_arrays(np.column_stack((gain,base,aug,ys)),n,b))
                        auc=[]
                        if h<3:
                            for k in np.unique(n):
                                ix=n==k
                                if np.unique(ys[ix]).size==2:auc.append(roc_auc_score(ys[ix],pm[ix])-roc_auc_score(ys[ix],pb[ix]))
                        rows.append(dict(ycut=ycut,family=fam,head=head,geometry=geometry,domain=label,rows=int(sel.sum()),
                            hits=int(ys.sum()) if h<3 else None,base_loss=av[1],model_loss=av[2],gain=av[0],lo=np.quantile(bt[:,0],.005),
                            hi=np.quantile(bt[:,0],.995),delta_AUC=float(np.mean(auc)) if auc else None,auc_eligible_N=len(auc)))
    write('results/by_regime.csv',rows)


if __name__=='__main__':
    with threadpool_limits(limits=1):render()
