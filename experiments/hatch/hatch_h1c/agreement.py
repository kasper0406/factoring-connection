"""Descriptive follow-up. Read the recorded fixtures; emit a pooled list curve and optional PDF."""
import argparse
import bisect
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from sympy import primerange

from experiment import V, full_lattice_sieve, interval_sieve, norm


def curves():
    # Extra resolution below the smallest preregistered threshold of J/4.
    thresholds = sorted({10**(-3 + i/80) for i in range(241)} | {0.25,0.5,1.0})
    totals = defaultdict(lambda: [0] * len(thresholds))
    for line in Path('runs.jsonl').read_text().splitlines():
        case = json.loads(line)
        start = case['interval']['X']
        p1,p2 = (list(primerange(y+1)) for y in case['bounds'])
        rem,_,_ = interval_sieve(start,p1)
        form = case['m2']['form']
        _,left,bad,_ = full_lattice_sieve(form,p2)
        parts = {'interval': [(start+i)//r for i,r in enumerate(rem)],
                 'm2': [norm(form,*z)//r//bad.get(z,1) for z,r in left.items()]}
        assert len(parts['interval']) == V
        for name,values in parts.items():
            J = case[name]['curve']['J']
            scores = sorted(math.log2(v)/J for v in values)
            counts = [len(scores)-bisect.bisect_left(scores,t) for t in thresholds]
            assert [counts[thresholds.index(t)] for t in (0.25,0.5,1.0)] == case[name]['curve']['johnson_counts']
            out = totals[name,case['bits']]
            for i,count in enumerate(counts):
                out[i] += count
    with Path('agreement.csv').open('w') as out:
        writer = csv.writer(out)
        writer.writerow(['domain','bits','agreement_over_J','pooled_list_size'])
        for (name,bits),counts in totals.items():
            writer.writerows((name,bits,t,k) for t,k in zip(thresholds,counts))
    return thresholds,totals


def plot(thresholds,totals):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(8.8,3.7),layout='constrained')
    colors = ('#176B70','#92400E','#6653A1','#922954')
    for ax,name,title in zip(axes,['interval','m2'],['Integer intervals','Primitive quadratic norms']):
        for bits,col in zip((40,64,80,96),colors):
            ax.step(thresholds,totals[name,bits],where='pre',label=f'{bits} bits',color=col,lw=1.6)
        ax.axvline(1,ls='--',color='0.45',lw=1)
        ax.set(xscale='log',yscale='symlog',ylim=(-0.15,300000),xlim=(0.001,1.2),
               xlabel='Agreement / Johnson benchmark',title=title)
        ax.set_yticks([0,1,10,100,1000,10000,100000],['0','1','10','100','1k','10k','100k'])
        ax.grid(which='major',axis='both',alpha=.16)
    axes[0].set_ylabel('List size across six fixtures')
    axes[1].legend(frameon=False,fontsize=8)
    fig.savefig('agreement.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot',action='store_true')
    args = parser.parse_args()
    result = curves()
    if args.plot:
        plot(*result)
