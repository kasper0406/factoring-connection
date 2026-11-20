"""Publish compact family-level tables with equal weight for moduli and geometries."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import t


def read(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def numeric(v):
    return float(v) if v not in ("",None) else np.nan


def interval(v):
    keep=np.asarray([a for a in v if np.isfinite(a)])
    if not len(keep):
        return [0,None,None,None,None]
    mean=float(keep.mean())
    if len(keep)==1:
        return [1,mean,None,None,None]
    se=float(keep.std(ddof=1)/np.sqrt(len(keep)))
    half=t.ppf(.975,len(keep)-1)*se
    return [len(keep),mean,se,mean-half,mean+half]


def table(path,cols,rows):
    with open(path, "w", newline="") as f:
        wr=csv.writer(f)
        wr.writerow(cols)
        wr.writerows(rows)


def summarize(root):
    table_rows, ledger_rows=[],[]
    for directory in sorted(root.glob("results/*")):
        if not (directory/"heldout_summary.json").exists():
            continue
        summary=json.loads((directory/"heldout_summary.json").read_text())
        mode,locked,matched=summary["mode"],summary["locked"],summary["matched"]
        rows=read(directory/"by_modulus_geometry.csv")
        lookup={(r["N_id"],r["geometry"],r["head"],r["model"]): r for r in rows}
        for d in ("3","4"):
            ids=sorted({r["N_id"] for r in rows if r["degree"]==d}, key=int)
            for head in ("joint","rational","algebraic"):
                for model in sorted({r["model"] for r in rows}):
                    capacity=model.split("_")[0]+"_base"
                    refs=["absolute"]
                    if model!="logistic":
                        refs+=["logistic"]
                        if model!=capacity: refs.append(capacity)
                        if model.endswith("history"): refs.append(model.replace("history","geometry"))
                    for ref in refs:
                        for metric in ("loss","auc"):
                            for label,geometries in [("both",["box","special_q"]),("box",["box"]),("special_q",["special_q"])]:
                                vals=[]
                                for k in ids:
                                    v=np.array([numeric(lookup[(k,g,head,model)][metric]) for g in geometries])
                                    if ref!="absolute":
                                        base=np.array([numeric(lookup[(k,g,head,ref)][metric]) for g in geometries])
                                        v=(base-v) if metric=="loss" else (v-base)
                                    vals.append(np.mean(v))
                                table_rows.append([mode,d,head,model,ref,label,metric,*interval(vals)])
        ledger=read(directory/"yield_ledger.csv")
        for d in ("all","3","4"):
            for g in ("both","box","special_q"):
                for b in range(-1,4):
                    subset=[r for r in ledger if (d=="all" or r["degree"]==d) and
                            (g=="both" or r["geometry"]==g) and int(r["u_bin"])==b]
                    n,pos,budget,base,control,model=(sum(int(r[k]) for r in subset)
                        for k in ["rows","positives","budget","logistic",matched,locked])
                    ledger_rows.append([mode,d,g,b,n,pos,budget,base,control,model,
                                        model/base if base else None,model/control if control else None])
    table(root/"results"/"family_metrics.csv",["experiment","degree","head","model","comparator","geometry",
          "metric","eligible_moduli","estimate","SE","CI95_low","CI95_high"],table_rows)
    table(root/"results"/"pooled_yield.csv",["experiment","degree","geometry","u_bin","candidates","positives",
          "checks","logistic_relations","matched_relations","selected_relations","selected_over_logistic",
          "selected_over_matched"],ledger_rows)


if __name__=="__main__":
    import sys
    for path in sys.argv[1:] or [".","power_cohort"]:
        summarize(Path(path))
