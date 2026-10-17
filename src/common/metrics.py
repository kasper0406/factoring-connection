"""Cost counters and JSONL metrics records (shared infrastructure, owner: infraA).

A single `Metrics` object is threaded through ring / baseline / candidate code.
It holds
  * `counts`  : additive event counters, e.g. ``rn_mul``, ``rn_add``, ``gcd``,
                ``inv``, ``poly_mul`` ...  (use ``inc``)
  * `maxima`  : running maxima, e.g. ``max_operand_bits``, ``max_degree``,
                ``stored_coeffs``, ``circuit_nodes``, ``ext_degree`` (use ``hi``)
  * `fields`  : free-form scalar fields (result, seed, timings ...) (use ``set``)

Records are serialised with `to_record()` and appended to JSONL files with
`append_jsonl`.  Everything is plain Python ints / strs so output is exact.
"""
from __future__ import annotations

import json
import os
import time
from typing import Any, Dict, Optional


class Metrics:
    __slots__ = ("counts", "maxima", "fields", "_t0", "_timers")

    def __init__(self, **fields: Any) -> None:
        self.counts: Dict[str, int] = {}
        self.maxima: Dict[str, int] = {}
        self.fields: Dict[str, Any] = dict(fields)
        self._t0 = time.perf_counter()
        self._timers: Dict[str, float] = {}

    # --- counters -------------------------------------------------------
    def inc(self, key: str, amount: int = 1) -> None:
        self.counts[key] = self.counts.get(key, 0) + amount

    def hi(self, key: str, value: int) -> None:
        """Record a running maximum."""
        if value > self.maxima.get(key, -1):
            self.maxima[key] = value

    def set(self, key: str, value: Any) -> None:
        self.fields[key] = value

    def get(self, key: str, default: int = 0) -> int:
        if key in self.counts:
            return self.counts[key]
        return self.maxima.get(key, default)

    def merge(self, other: "Metrics", prefix: str = "") -> None:
        for k, v in other.counts.items():
            self.inc(prefix + k, v)
        for k, v in other.maxima.items():
            self.hi(prefix + k, v)

    # --- timing ---------------------------------------------------------
    def elapsed(self) -> float:
        return time.perf_counter() - self._t0

    def start(self, name: str) -> None:
        self._timers[name] = time.perf_counter()

    def stop(self, name: str) -> float:
        dt = time.perf_counter() - self._timers.pop(name)
        self.fields["time_" + name] = self.fields.get("time_" + name, 0.0) + dt
        return dt

    # --- output ---------------------------------------------------------
    def to_record(self, **extra: Any) -> Dict[str, Any]:
        rec: Dict[str, Any] = dict(self.fields)
        rec.update(extra)
        rec["counts"] = dict(sorted(self.counts.items()))
        rec["maxima"] = dict(sorted(self.maxima.items()))
        rec.setdefault("wall_s", round(self.elapsed(), 6))
        return rec

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Metrics(counts={self.counts}, maxima={self.maxima}, fields={self.fields})"


def _default(o: Any) -> Any:
    if hasattr(o, "__int__"):
        return int(o)
    return str(o)


def append_jsonl(path: str, record: Dict[str, Any]) -> None:
    """Append one JSON record (exact ints preserved) to `path`, creating dirs."""
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "a") as fh:
        fh.write(json.dumps(record, default=_default, sort_keys=False) + "\n")


def read_jsonl(path: str):
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


class NullMetrics(Metrics):
    """Metrics sink that ignores everything (for speed-critical inner loops)."""

    def inc(self, key: str, amount: int = 1) -> None:  # noqa: D401
        pass

    def hi(self, key: str, value: int) -> None:
        pass


def ensure(m: Optional[Metrics]) -> Metrics:
    return m if m is not None else NullMetrics()
