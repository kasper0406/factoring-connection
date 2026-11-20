"""Compact crosswalk of registered gains and family-native sensitivity."""
import csv
from pathlib import Path
from rstats import write
from rfit import NAMES
from rbase import YS

def read(path):
    with open(path,newline='') as f:return list(csv.DictReader(f))

def main():
    orig=read('results/A_families.csv');native=read('results/native_sensitivity.csv')
    gates=read('results/native_gates_0.csv')+read('results/native_gates_1.csv'); rows=[]
    for fam in NAMES[:8]:
        for y in map(str,YS):
            r=next(t for t in orig if t['ycut']==y and t['family']==fam and t['head']=='joint')
            row=dict(y=y,family=fam,eligible_rows=r['rows'],hits=r['hits'],joint_gain=r['gain'],lo99=r['ci_lo'],
                     hi99=r['ci_hi'],WY_validation=r['val_WY'],Holm_test=r['Holm'],discovery=r['declared'])
            for h in ('joint','algebraic'):
                q=next(t for t in native if t['ycut']==y and t['family']==fam and t['head']==h)
                v,t=[g for g in gates if g['ycut']==y and g['family']==fam and g['head']==h]
                row[h+'_2x_test_power']=q['conditional_test_power']
                row[h+'_2x_observed_test_pass']=q['detected']
                row[h+'_2x_all_gates']=bool(q['detected']=='True' and float(v['val_or_test_WY'])<=.025 and
                    float(t['cell_p'])<.01 and float(t['per_hit_lo'])>.02)
            rows.append(row)
    write('results/family_summary.csv',rows)
if __name__=='__main__': main()
