"""COLL-8 analysis: the ordering across seven cells.

Rates are normalised to spikes per second because PY 1 and AB/PD 2 were run
at a 30 s horizon (COLL-6, COLL-5) and the five new cells at 35 s.
"""
import glob
import itertools
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C5 = os.path.join(HERE, '..', 'coll5')
C6 = os.path.join(HERE, '..', 'coll6')


def pair_from(files_off, files_on, secs):
    off, on = {}, {}
    for f in files_off:
        for r in json.load(open(f))['rows']:
            off[r['seed']] = r
    for f in files_on:
        for r in json.load(open(f))['rows']:
            on[r['seed']] = r
    s = sorted(set(off) & set(on))
    n = len(s)
    saved = sum(1 for x in s if not off[x]['alive'] and on[x]['alive'])
    killed = sum(1 for x in s if off[x]['alive'] and not on[x]['alive'])
    dead = [x for x in s if not off[x]['alive']]
    live = [x for x in s if off[x]['alive']]
    rate = float(np.mean([off[x]['total_spikes'] for x in s])) / secs
    return dict(n=n, saved=saved, killed=killed, diff=(saved - killed) / n,
                rate=rate, n_dead=len(dead), n_live=len(live),
                kill_healthy=(killed / len(live)) if live else float('nan'),
                rescue=(saved / len(dead)) if dead else float('nan'))


cells = {}
for c, slug in (('AB/PD 1', 'ABPD1'), ('AB/PD 5', 'ABPD5'), ('AB/PD 3', 'ABPD3')):
    cells[c] = pair_from([f'{HERE}/coll8a_{slug}_off.json'],
                         [f'{HERE}/coll8a_{slug}_on.json'], 35.0)
for c, slug in (('LP 3', 'LP3'), ('AB/PD 4', 'ABPD4')):
    cells[c] = pair_from([f'{HERE}/coll8b_{slug}_off.json'],
                         [f'{HERE}/coll8b_{slug}_on.json'], 35.0)
cells['PY 1'] = pair_from(sorted(glob.glob(f'{C6}/coll6_lot*_off.json')),
                          sorted(glob.glob(f'{C6}/coll6_lot*_on.json')), 30.0)
cells['AB/PD 2'] = pair_from([f'{C5}/coll5_g0.045_off.json'],
                             [f'{C5}/coll5_g0.045_on.json'], 30.0)

order = sorted(cells, key=lambda c: cells[c]['rate'])
print("COLL-8  seven cells, ordered by the OFF arm's own firing rate at g_ie 0.045")
print("=" * 86)
print(f"{'cell':>9} {'spikes/s':>9} {'n':>4} {'OFF dead':>9} {'saved':>6} "
      f"{'killed':>7} {'diff':>8} {'kill|healthy':>13} {'rescue|dead':>12}")
for c in order:
    d = cells[c]
    rs = f"{d['saved']}/{d['n_dead']}" if d['n_dead'] else '   --'
    print(f"{c:>9} {d['rate']:>9.1f} {d['n']:>4} {d['n_dead']:>9} {d['saved']:>6} "
          f"{d['killed']:>7} {d['diff']:>+8.3f} {d['kill_healthy']:>13.3f} {rs:>12}")


def spearman(x, y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def exact_p(x, y):
    n = len(x)
    obs = abs(spearman(x, y))
    base = list(range(n))
    tot = cnt = 0
    for p in itertools.permutations(base):
        tot += 1
        cnt += abs(spearman(base, list(p))) >= obs - 1e-9
    return obs, cnt / tot, tot


rates = [cells[c]['rate'] for c in order]
for lbl, key in (('diff = saved - killed (ceiling moves with OFF-death rate)', 'diff'),
                 ('kill rate on healthy wirings (no such ceiling)', 'kill_healthy')):
    vals = [cells[c][key] for c in order]
    rho, p, tot = exact_p(rates, vals)
    sgn = spearman(rates, vals)
    print(f"\n  {lbl}")
    print(f"    Spearman rho = {sgn:+.3f}   |rho| = {rho:.3f}   "
          f"exact two-sided p = {p:.4f}  ({tot} permutations)")
    print(f"    registered bar |rho| >= 0.786: "
          f"{'MET' if rho >= 0.786 else 'NOT MET'}")
