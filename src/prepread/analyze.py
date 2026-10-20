"""Summaries of prepread raw data (prints tables used in notes/prepread_results.md §6)."""
from __future__ import annotations
import json, os, statistics, sys
from collections import defaultdict

RAW = os.path.join(os.path.dirname(__file__), '..', '..', 'experiments', 'raw')


def load(name):
    path = os.path.join(RAW, f'prepread_{name}.jsonl')
    return [json.loads(l) for l in open(path)] if os.path.exists(path) else []


def e1(tag='e1'):
    rows = load(tag)
    base = {r['p']: r for r in rows if r.get('kind') == 'baseline'}
    mem = [r for r in rows if r.get('kind') == 'member']
    print(f'== {tag}: primes', sorted(base), 'B_p', {p: round(b['B_p'], 4) for p, b in sorted(base.items())})
    fam = defaultdict(list)
    for r in mem:
        fam[r['family']].append(r)
    for f, rs in sorted(fam.items()):
        nd = [r for r in rs if r['nondeg']]
        dg = [r for r in rs if not r['nondeg']]
        mx_nd = max((r['ratio'] for r in nd), default=None)
        mx_dg = max((r['ratio'] for r in dg), default=None)
        hi = [r for r in rs if r.get('width') is not None and r['width'] >= (r['p'] ** 0.5)]
        print(f, 'n', len(rs), 'nondeg', len(nd), 'max ratio nondeg', mx_nd and round(mx_nd, 2),
              'max ratio degen', mx_dg and round(mx_dg, 2),
              'median ratio (width>=sqrt p)', hi and round(statistics.median(r['ratio'] for r in hi), 2),
              'max ratio (width>=sqrt p)', hi and round(max(r['ratio'] for r in hi), 2), 'n_hi', len(hi))
    weil = [r for r in mem if 'weil_ok' in r]
    print('H2 Weil checks', len(weil), 'violations', sum(not r['weil_ok'] for r in weil),
          'max M/bound', weil and round(max(r['M_allnonzero'] / r['weil_bound'] for r in weil), 3))
    h1 = [r for r in mem if r['family'] in ('F1', 'F2', 'F3') and r['nondeg']]
    print('H1 nondeg F1-F3', len(h1), 'violations (ratio>2)', sum(r['ratio'] > 2 for r in h1))
    allh = [r for r in mem if r['family'] in ('F1', 'F2', 'F3', 'E1b')]
    print('all non-control ratio max', round(max(r['ratio'] for r in allh), 3),
          'kill (ratio>=10)', sum(r['ratio'] >= 10 for r in allh))
    for r in mem:
        if r['family'] == 'F4':
            print('F4', r['p'], round(r['M'], 4), round(r['ratio'], 1))
    top = sorted(allh, key=lambda r: -r['ratio'])[:8]
    for r in top:
        print('top', r['name'], r['p'], round(r['ratio'], 3), 'image', r['image'], 'width', r.get('width'))


if __name__ == '__main__':
    for a in sys.argv[1:]:
        globals()[a.split(':')[0]](*a.split(':')[1:])
