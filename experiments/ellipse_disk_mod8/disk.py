"""Symbolic disk over coordinate and quotient bits. Count assignments in Z/8Z."""

from collections import Counter
from math import isqrt
import warnings
from dd import cudd


def ring_count(manager, root, variables, work, bits):
    memo = {}

    def edge(node, level):
        if node == manager.false:
            return 0
        if node == manager.true:
            gap = variables - level
            return 0 if gap >= bits else 1 << gap
        gap = node.level - level
        if gap >= bits:
            return 0
        key = int(node)
        if key not in memo:
            pos, lo, hi = manager.succ(node)
            if node.negated:
                lo, hi = ~lo, ~hi
            memo[key] = (edge(lo, pos + 1) + edge(hi, pos + 1)) & ((1 << bits) - 1)
        return (memo[key] << gap) & ((1 << bits) - 1)

    result = edge(root, 0)
    work["count_nodes"] += len(memo)
    return result


def count_residue(manager, planes, value):
    f = manager.true
    for bit, plane in enumerate(planes):
        f &= plane if (value >> bit) & 1 else ~plane
    return f


def forget_digit(manager, states, name, work, parity):
    """Sum out a coordinate bit at fixed carry in the selected coefficient ring."""
    if parity:
        updated = {}
        for (c, mult), f in states.items():
            left = manager.let({name: False}, f)
            right = manager.let({name: True}, f)
            if (left & right) != manager.false:
                work["overlapping_cofactors"] += 1
            merged = manager.apply("xor", left, right)
            if merged != manager.false:
                updated[c, 1] = merged
        return updated
    by_carry = {}
    for (carry, mult), predicate in states.items():
        planes = by_carry.setdefault(carry, [manager.false] * 3)
        for bit in range(3):
            if (mult >> bit) & 1:
                planes[bit] |= predicate
    updated = {}
    for carry, planes in by_carry.items():
        left = [manager.let({name: False}, f) for f in planes]
        right = [manager.let({name: True}, f) for f in planes]
        overflow = manager.false
        total = []
        for l, r in zip(left, right):
            half_sum = manager.apply("xor", l, r)
            total.append(manager.apply("xor", half_sum, overflow))
            overflow = (l & r) | (overflow & half_sum)
        if overflow != manager.false:
            work["overflow_carry_classes"] += 1
        for mult in range(1, 8):
            f = count_residue(manager, total, mult)
            if f != manager.false:
                updated[carry, mult] = f
    return updated


class DiskCounter:
    """One symbolic relation, all B up to cap. Bounds ensure exact integer equality."""

    def __init__(
        self,
        br,
        order="lsb",
        limit=None,
        early=False,
        parity=False,
        project_quotient=False,
    ):
        if order not in ("lsb", "msb", "grouped"):
            raise ValueError("unknown variable order")
        self.br = br
        self.work = Counter()
        self.parity = parity
        self.early = early
        self.bdd = m = cudd.BDD()
        m.configure(reordering=False)
        self.limit = br.cap if limit is None else limit
        self.project_quotient = project_quotient
        if not 1 <= self.limit <= br.cap:
            raise ValueError("invalid disk bound")
        top = 4 * br.n * self.limit
        sizes = dict(
            w=isqrt(top).bit_length(),
            y=isqrt(top // br.delta).bit_length(),
            j=self.limit.bit_length(),
        )
        names = [f"{v}{i}" for v in sizes for i in range(sizes[v])]
        if order != "grouped":
            names.sort(key=lambda name: int(name[1:]) * (1 if order == "lsb" else -1))
        m.declare(*names)
        self.variables = len(names)
        vb = {v: [m.var(f"{v}{i}") for i in range(n)] for v, n in sizes.items()}
        self.y = vb["y"]
        self.j = vb["j"]
        # Beyond this width every summand is zero and the final signed carry is exact.
        bits = (
            max(
                (
                    2 * ((1 << sizes["w"]) - 1) ** 2 + br.delta * (1 << sizes["y"]) ** 2
                ).bit_length(),
                (4 * br.n * (1 << sizes["j"])).bit_length(),
            )
            + 1
        )
        columns = [[] for _ in range(bits)]

        def term(c, node, offset):
            sign = 1 if c > 0 else -1
            c = abs(c)
            i = 0
            while c:
                if c & 1:
                    columns[offset + i].append((sign, node))
                c >>= 1
                i += 1

        for v, c in [("w", 1), ("y", br.delta)]:
            for i, x in enumerate(vb[v]):
                term(c, x, 2 * i)
                for k in range(i):
                    term(c, x & vb[v][k], i + k + 1)
        for i, x in enumerate(vb["j"]):
            term(-4 * br.n, x, i)
        last_use = {}
        for index, terms in enumerate(columns):
            for _, condition in terms:
                for name in condition.support:
                    if name.startswith("w") or (
                        project_quotient and name.startswith("j")
                    ):
                        last_use[name] = index

        def positive(digits):
            f = m.false
            for digit in digits:
                f |= digit
            return f

        valid = positive(self.y) & positive(vb["w"]) & self.leq(self.limit)
        states = {(0, 1): valid}
        for index, terms in enumerate(columns):
            for sign, node in terms:
                fresh = {}
                for (c, mult), f in states.items():
                    for idx, g in (
                        ((c, mult), f & ~node),
                        ((c + sign, mult), f & node),
                    ):
                        if g != m.false:
                            fresh[idx] = g if idx not in fresh else fresh[idx] | g
                states = fresh
                self.work["max_carry_states"] = max(
                    self.work["max_carry_states"], len(states)
                )
            states = {
                (c // 2, mult): f for (c, mult), f in states.items() if c % 2 == 0
            }
            if early:
                for name in [n for n, last in last_use.items() if last == index]:
                    states = forget_digit(m, states, name, self.work, parity)
                    self.work[name[0] + "_bits_eliminated"] += 1
        self.masses = {mult: f for (c, mult), f in states.items() if c == 0}
        # Excluded coordinates have already been summed out and must count once.
        # Pin them to zero solely for the final traversal over the original order.
        if early:
            pinned = m.true
            for digit in vb["w"] + (vb["j"] if project_quotient else []):
                pinned &= ~digit
            self.masses = {mult: f & pinned for mult, f in self.masses.items()}
        self.work["final_relation_nodes"] = sum(len(f) for f in self.masses.values())
        self.work["variables"] = self.variables

    def leq(self, B):
        f = self.bdd.true
        for i, x in enumerate(self.j):
            f = (~x | f) if (B >> i) & 1 else (~x & f)
        return f

    def count(self, B, bit=None):
        bits = () if bit is None else bit if isinstance(bit, tuple) else (bit,)
        if not 0 <= B <= self.limit or any(i < 0 for i in bits):
            raise ValueError("outside query range")
        self.work["queries"] += 1
        if B == 0 or any(i >= len(self.y) for i in bits):
            return 0
        if self.project_quotient and B != self.limit:
            raise ValueError("quotient bits were projected after fixing this cutoff")
        condition = self.bdd.true if self.project_quotient else self.leq(B)
        for i in bits:
            condition &= self.y[i]
        return sum(
            mult
            * ring_count(
                self.bdd,
                f & condition,
                self.variables,
                self.work,
                1 if self.parity else 3,
            )
            for mult, f in self.masses.items()
        ) & (1 if self.parity else 7)

    def metrics(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            stats = self.bdd.statistics()
        return dict(
            self.work,
            peak_live_nodes=stats["peak_live_nodes"],
            cache_lookups=int(stats["cache_lookups"]),
        )


class AdaptiveDisk:
    """Build only growing disk ranges actually used by the search."""

    def __init__(self, br, order="lsb", early=False, parity=False):
        self.br = br
        self.order = order
        self.parity = parity
        self.early = early
        self.current = None
        self.completed = []

    def count(self, B, bit=None):
        if not 0 <= B <= self.br.cap:
            raise ValueError("outside query range")
        if B == 0:
            return 0
        if self.current is None or self.current.limit < B:
            target = B
            if self.current is not None:
                target = max(target, 2 * self.current.limit)
                self.completed.append(self.current.metrics())
                self.current = None
            self.current = DiskCounter(
                self.br,
                self.order,
                min(target, self.br.cap),
                early=self.early,
                parity=self.parity,
            )
        return self.current.count(B, bit)

    def metrics(self):
        builds = [*self.completed, self.current.metrics()]
        return dict(
            builds=len(builds),
            max_disk=self.current.limit,
            max_peak_live_nodes=max(b["peak_live_nodes"] for b in builds),
            total_cache_lookups=sum(b["cache_lookups"] for b in builds),
            final_relation_nodes=builds[-1]["final_relation_nodes"],
            overflow_carry_classes=sum(
                b.get("overflow_carry_classes", 0) for b in builds
            ),
            overlapping_cofactors=sum(
                b.get("overlapping_cofactors", 0) for b in builds
            ),
            count_nodes=sum(b["count_nodes"] for b in builds),
            query_calls_excluding_zero=sum(b["queries"] for b in builds),
        )


class FixedDisks:
    """Eliminate w and j at a fixed cutoff. Charge a separate diagram per cutoff."""

    def __init__(self, br):
        self.br = br
        self.disks = {}

    def count(self, B, bit=None):
        if B == 0:
            return 0
        if B not in self.disks:
            self.disks[B] = DiskCounter(
                self.br, limit=B, early=True, parity=True, project_quotient=True
            )
        return self.disks[B].count(B, bit)

    def metrics(self):
        ms = [d.metrics() for d in self.disks.values()]
        return dict(
            builds=len(ms),
            max_peak_live_nodes=max(m["peak_live_nodes"] for m in ms),
            total_cache_lookups=sum(m["cache_lookups"] for m in ms),
            sum_final_nodes=sum(m["final_relation_nodes"] for m in ms),
            overlapping_cofactors=sum(m.get("overlapping_cofactors", 0) for m in ms),
        )
