"""Run the pre-announced stronger sensitivity control if the original misses."""
import json
import time
from pathlib import Path

import joblib
from threadpoolctl import threadpool_limits

from evaluate import evaluate
from fit_models import fit
from experiment import dump

if __name__ == "__main__":
    if not json.loads(Path("results/planted/heldout_summary.json").read_text())["crosses_rule"]:
        wall,cpu=time.monotonic(),time.process_time()
        with threadpool_limits(limits=1,user_api="blas"),threadpool_limits(limits=2,user_api="openmp"):
            size,cal,_=joblib.load("results/true/models.joblib")
            fit("planted_strong",(size,cal))
            evaluate("planted_strong",(size,cal))
        dump("logs/strong_control_time.json",dict(wall=time.monotonic()-wall,cpu=time.process_time()-cpu))
