"""Train and lock all three experiments without loading final test data."""
import gc
import hashlib
import json
import os
import time
from pathlib import Path
import warnings

import joblib
import numpy as np
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
from experiment import *
from settings import FIT_LIMIT

warnings.filterwarnings("ignore", message="Stochastic Optimizer: Maximum iterations",category=ConvergenceWarning)


def fit(mode, real_teacher=None):
    start = time.monotonic()
    train = load_split("fit",mode,real_teacher)
    cal = load_split("cal",mode,real_teacher)
    sel = load_split("select",mode,real_teacher)
    # Planted control modifies only the joint outcome; its two size regressions are identical.
    size = real_teacher[0] if mode.startswith("planted") else CrossSize().fit(train)
    rt, rc, rs = size.scores(train,True), size.scores(cal), size.scores(sel)
    params = calibration_for(cal.y,rc)
    pcal, psel = calibrated(rc,params), calibrated(rs,params)
    take = np.sort(np.random.default_rng(SEED).choice(len(train.n),FIT_LIMIT,replace=False))
    xt, yt = make_x(train,rt,take), train.y[take]
    xc, xs = make_x(cal,rc), make_x(sel,rs)
    del train, rt
    base_loss = grouped_mean(loss(sel.y[:,0], psel[:,0]), sel)
    score, fitted = {"logistic":base_loss}, {}
    for form in FORMS:
        t = time.monotonic()
        model = Flexible(form).fit(xt,yt).calibrate(xc,cal.y,pcal)
        p = model.predict(xs,psel)
        score[form] = grouped_mean(loss(sel.y[:,0],p[:,0]), sel)
        fitted[form] = model
        print(mode,form,round(time.monotonic()-t,1),score[form],flush=True)
    margins = {f:min(base_loss, score[f.split("_")[0]+"_base"]) - score[f] for f in EXTRAS}
    winner = max(EXTRAS,key=margins.get)
    root = Path("results") / mode
    root.mkdir(exist_ok=True)
    joblib.dump((size,params,fitted),root/"models.joblib",compress=2)
    lock = dict(mode=mode, selected=winner, selection_joint_loss=score, selection_min_gain=margins,
                weights={f:m.weights.tolist() for f,m in fitted.items()},size_calibration=params.tolist(),
                seconds=time.monotonic()-start, model_sha256=hashlib.sha256((root/"models.joblib").read_bytes()).hexdigest(),
                prereg_sha256=hashlib.sha256(Path("PREREGISTRATION.md").read_bytes()).hexdigest())
    dump(root/"selection_lock.json",lock)
    return size,params


if __name__ == "__main__":
    start, cpu = time.monotonic(), time.process_time()
    with threadpool_limits(limits=1, user_api="blas"), threadpool_limits(limits=2,user_api="openmp"):
        prepare_edges()
        teacher = fit("true")
        gc.collect()
        fit("shuffled")
        gc.collect()
        fit("planted", teacher)
    dump("logs/training_time.json",dict(wall=time.monotonic()-start,cpu=time.process_time()-cpu))
