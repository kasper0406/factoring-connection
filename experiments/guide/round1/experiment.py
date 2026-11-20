"""Feature allowlist, split loaders, and estimators used by fit and final testing."""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logit
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

from settings import BOUND, PER_FAMILY, POWER, REFS, SEED
EPS = 1e-6
HEADS = ["joint", "rational", "algebraic"]
MODES = ["true", "shuffled", "planted"]
FORMS = ["boost_base", "boost_geometry", "boost_history", "mlp_base", "mlp_geometry", "mlp_history"]
EXTRAS = [key for key in FORMS if not key.endswith("_base")]
SIZES = [27, 33, 39]  # both log sizes + 22 valuations + 3 out-of-fold baseline logits; 6 geometry, 6 history


def manifest():
    return json.loads(Path("data/manifest.json").read_text())


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def clip(p):
    return np.clip(p, EPS, 1-EPS)


def loss(y, p):
    p = clip(np.asarray(p, dtype=np.float64))
    return -np.where(y, np.log(p), np.log1p(-p))


def geometry(coords):
    a, b, i, g = coords.T
    x, t, h = i // 32, (b-1)//32, a//256
    # Stable integer hash of coarse cells. The salt is fixed before any labels enter.
    key = (x.astype(np.int64)+17) + 43*t + 1723*g
    key = ((key ^ (key >> 16))*0x45d9f3b) & 0xffffffff
    key = ((key ^ (key >> 16))*0x45d9f3b) & 0xffffffff
    bit = ((key ^ (key >> 16)) & 1)
    return np.column_stack((g, x, t, h, (a<0).astype(int), bit)).astype(np.float32)


def guard_table(limit):
    """Sieve determinants directly; table[d] says no forbidden prime divides d."""
    prime = np.ones(BOUND+1, dtype=np.bool_)
    prime[:2] = False
    for p in range(2, int(BOUND**0.5)+1):
        if prime[p]:
            prime[p*p::p] = False
    allowed = np.ones(limit+1, dtype=np.bool_)
    for p in np.flatnonzero(prime & (np.arange(BOUND+1)>31)):
        allowed[p::p] = False
    allowed[0] = False
    return allowed


def history_edges(coords, table):
    a, b = coords[REFS:, :2].astype(np.int64).T
    active = np.arange(len(a))
    graph = np.empty((len(a),16), dtype=np.uint16 if POWER else np.uint8)
    counts = np.zeros(len(a),dtype=np.int32)
    for j, (ra,rb) in enumerate(coords[:REFS,:2]):
        if not len(active):
            break
        det = np.abs(a[active]*rb-b[active]*ra)
        accepted = active[table[det]]
        graph[accepted, counts[accepted]] = j
        counts[accepted] += 1
        active = active[counts[active] < 16]
    assert np.all(counts == 16)
    return graph


def prepare_edges():
    t = time.monotonic()
    # |a| <= 31*224 + 15*b; donors b<=768, targets b<=512.
    table = guard_table((31*224+15*512)*768 + (31*224+15*768)*512)
    for item in manifest():
        k = item["id"]
        with np.load(f"data/n{k:02}.npz") as z:
            graph = history_edges(z["coords"], table)
        np.savez_compressed(f"data/edges{k:02}.npz", indices=graph)
    dump("logs/history_construction.json", {"seconds":time.monotonic()-t})


@dataclass
class Batch:
    base: np.ndarray
    geo: np.ndarray
    hist: np.ndarray
    y: np.ndarray
    n: np.ndarray
    oracle: np.ndarray


def load_split(split, mode="true", real_teacher=None):
    parts = [[] for _ in range(6)]
    global_outcomes = {}
    if mode == "global_shuffle":
        ids = [item["id"] for item in manifest() if item["split"] == split]
        pairs = []
        for k in ids:
            with np.load(f"data/n{k:02}.npz") as z:
                pairs.append(z["sides"])
        rng = np.random.default_rng(SEED + 12941 + ["fit","cal","select","test"].index(split))
        references = rng.permutation(np.concatenate([p[:REFS] for p in pairs]))
        targets = rng.permutation(np.concatenate([p[REFS:] for p in pairs]))
        offset = 0
        for j,(k,pair) in enumerate(zip(ids,pairs)):
            length = len(pair)-REFS
            global_outcomes[k] = np.vstack((references[j*REFS:(j+1)*REFS],targets[offset:offset+length]))
            offset += length
    for it in manifest():
        if it["split"] != split:
            continue
        k = it["id"]
        with np.load(f"data/n{k:02}.npz") as z:
            coords, sides = z["coords"], z["sides"]
            words, val = z["cofactor_words"], z["valuation"]
        with np.load(f"data/edges{k:02}.npz") as z:
            edges = z["indices"]
        size = words[:,:,0].astype(np.float64) + words[:,:,1].astype(np.float64)*2.0**64
        base = np.column_stack((np.log(size)/np.log(BOUND), val.reshape(len(val),22))).astype(np.float32)
        geo = geometry(coords)
        if mode == "global_shuffle":
            sides = global_outcomes[k]
        y = np.column_stack((sides.prod(axis=1), sides)).astype(np.uint8)
        rng = np.random.default_rng(SEED+193*k+(7907 if mode=="shuffled" else 7909))
        oracle = np.zeros(len(y), dtype=np.float32)
        if mode == "shuffled":
            for isref in [True, False]:
                for g in [0,1]:
                    ix = np.flatnonzero(((np.arange(len(y))<REFS)==isref) & (coords[:,3]==g))
                    y[ix] = y[rng.permutation(ix)]
        elif mode in ("planted", "planted_strong"):
            # Fit points and their references receive predictions from a teacher excluding their N.
            teacher, calibration = real_teacher
            prediction = calibrated(teacher.predict(base, k if split=="fit" else None), calibration)[:,0]
            amplitude = 0.70 if mode == "planted_strong" else 0.35
            oracle = expit(logit(clip(prediction)) + amplitude*(2*geo[:,-1]-1))
            y[:,0] = rng.random(len(y)) < oracle
        ref_y = y[:REFS][edges]
        hist = np.column_stack((ref_y[:,:4].mean(axis=1), ref_y.mean(axis=1))).astype(np.float32)
        for out, v in zip(parts, (base[REFS:], geo[REFS:], hist, y[REFS:],
                                  np.full(len(edges), k, dtype=np.int16), oracle[REFS:])):
            out.append(v)
    return Batch(*(np.concatenate(p) for p in parts))


def side_design(x, side):
    u = x[:,side]
    knots = np.array([1,1.2,1.4,1.6,1.8,2,2.25,2.5,2.75,3,3.5,4,4.5,5], dtype=np.float64)
    return np.column_stack((x, u <= 1, np.maximum(0, u[:,None]-knots)))


class SizeModel:
    def __init__(self):
        self.models = []

    def fit(self, x, y):
        for side in [0,1]:
            model = LogisticRegression(C=2, solver="newton-cholesky", max_iter=50)
            model.fit(side_design(x,side), y[:,1+side])
            self.models.append(model)
        return self

    def predict(self, x):
        p = [m.predict_proba(side_design(x,j))[:,1] for j,m in enumerate(self.models)]
        return np.column_stack((p[0]*p[1], *p))


class CrossSize:
    def fit(self, batch):
        self.folds = {}
        # Each fold has two fit moduli of each degree, selected by within-family index.
        self.fold_for = {int(k): (int(k)%PER_FAMILY)%4 for k in np.unique(batch.n)}
        assignment = np.array([self.fold_for[int(k)] for k in batch.n], dtype=np.uint8)
        for f in range(4):
            keep = assignment != f
            self.folds[f] = SizeModel().fit(batch.base[keep], batch.y[keep])
        self.full = SizeModel().fit(batch.base, batch.y)
        return self

    def predict(self, x, fit_n=None):
        model = self.full if fit_n is None else self.folds[self.fold_for[int(fit_n)]]
        return model.predict(x)

    def scores(self, batch, oof=False):
        if not oof:
            return self.predict(batch.base)
        out = np.empty((len(batch.base),3))
        for k in np.unique(batch.n):
            mask = batch.n==k
            out[mask] = self.predict(batch.base[mask], k)
        return out


def sigmoid_fit(y,p):
    z = logit(clip(p)).astype(np.float64)
    def objective(v):
        a = v[0] + v[1]*z
        resid = expit(a)-y
        return np.mean(np.logaddexp(0,a)-y*a), np.array([np.mean(resid),np.mean(z*resid)])
    opt = minimize(objective, [0,1], jac=True, method="L-BFGS-B", bounds=[(-25,25),(0,12)],
                   options={"gtol":1e-9,"ftol":1e-13})
    if not opt.success:
        raise ValueError(opt.message)
    return opt.x


def calibration_for(y, p):
    return np.stack([sigmoid_fit(y[:,j],p[:,j]) for j in range(3)])


def calibrated(p, params):
    return clip(expit(params[:,0] + params[:,1]*logit(clip(p))))


def make_x(batch, raw, take=None):
    ix = slice(None) if take is None else take
    return np.column_stack((batch.base[ix], logit(clip(raw[ix])).astype(np.float32),
                            batch.geo[ix], batch.hist[ix]))


def dimensions(form):
    return SIZES[["base","geometry","history"].index(form.split("_")[1])]


class Flexible:
    def __init__(self, form):
        self.form = form
        self.dim = dimensions(form)

    def fit(self, x, y):
        z = x[:,:self.dim]
        if self.form.startswith("boost"):
            self.models = []
            for j in range(3):
                model = HistGradientBoostingClassifier(max_iter=85, max_leaf_nodes=15,
                    min_samples_leaf=80, learning_rate=.08, l2_regularization=3,
                    early_stopping=False, random_state=SEED)
                model.fit(z,y[:,j])
                self.models.append(model)
        else:
            self.scaler = StandardScaler().fit(z)
            self.model = MLPClassifier(hidden_layer_sizes=(40,20), batch_size=512,
                learning_rate_init=.001, alpha=.5, max_iter=18, tol=0,
                n_iter_no_change=20, random_state=SEED)
            self.model.fit(self.scaler.transform(z), y)
        return self

    def raw(self,x):
        z = x[:,:self.dim]
        if self.form.startswith("boost"):
            return np.column_stack([m.predict_proba(z)[:,1] for m in self.models])
        return self.model.predict_proba(self.scaler.transform(z)).astype(np.float64)

    def calibrate(self,x,y,p0):
        raw = self.raw(x)
        self.cal = calibration_for(y,raw)
        new = calibrated(raw, self.cal)
        self.weights = []
        for j in range(3):
            def objective(w):
                return loss(y[:,j], (1-w)*p0[:,j]+w*new[:,j]).mean()
            opt = minimize_scalar(objective,bounds=(0,1),method="bounded")
            self.weights.append(min([0, float(opt.x), 1], key=objective))
        self.weights = np.array(self.weights)
        return self

    def predict(self,x,p0):
        return clip(p0 + self.weights*(calibrated(self.raw(x),self.cal)-p0))

    def joint(self,x,p0):
        if self.weights[0] == 0:
            return p0
        z=x[:,:self.dim]
        if self.form.startswith("boost"):
            raw = self.models[0].predict_proba(z)[:,1]
        else:
            raw = self.model.predict_proba(self.scaler.transform(z))[:,0].astype(np.float64)
        new = clip(expit(self.cal[0,0]+self.cal[0,1]*logit(clip(raw))))
        return clip(p0+self.weights[0]*(new-p0))


def grouped_mean(values, batch):
    return float(np.mean([values[(batch.n==k)&(batch.geo[:,0]==g)].mean()
                          for k in np.unique(batch.n) for g in (0,1)]))
