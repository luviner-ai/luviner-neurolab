"""COLL-8b: the two cells that bracket the declared crossing.

Registered in BRAIN-SIMULATION.md, `## Run: COLL-8 (registration)`, commit
6d2a958, before any arm ran. The per-cell sign is not buyable -- 98 wirings
near the crossing -- so the claim is the Spearman ordering across seven cells
and each cell needs a stable estimate rather than a significant one.

LP 3 and AB/PD 4 sit either side of the crossing declared at 6,596 spikes,
where the sign is nearest zero and hardest to estimate, so they get 16 wirings
against 8a's 10 -- the orchestrator's own mitigation, in the direction they
proposed. A rank-based trend test is unaffected by the unequal n.

LP 3 is also the conductance set used as the INHIBITORY population in every
dish here, so in this row alone the E and I cells are the same set. Declared
in the registration, not discovered afterwards.

Written out rather than sed-derived from coll6.py. Two rows in a row lost
something a parent script had; the rule recorded at PHASE-T4 and again at
COLL-6 step B is that a derived script needs its consumers checked, and the
cheapest way to check them is to write them.
"""

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
sys.path.insert(0, os.path.join(HERE, '..', 'dish'))
sys.path.insert(0, os.path.join(HERE, '..', 'coll4'))
sys.path.insert(0, os.path.join(HERE, '..', 'coll5'))

from price import build
from stgbatch import BatchedSTGDish
from luviner.biophysics import stg
from stg_pm import N_EXC, burst_windows

CELLS = ['LP 3', 'AB/PD 4']
G_IE = 0.045
SEEDS = list(range(30, 44))          # 14 wirings per cell
# 16 was priced from AB/PD 2's 1.095 s/s constant, which does not
# transfer across cells: COLL-8a measured 84.7 s per wiring against
# the 76.6 priced, so 16 would have been 46.7 min against a 45 cap.
# The SE of (saved - killed)/n moves 0.177 -> 0.189, which a
# rank-based trend test does not notice.
T_MS = 35000.0
CHUNK_MS = 5000.0
Q = 24
BUDGET = (N_EXC - 1) * 0.15 * 0.015
STOP_MIN = 44.0
TAG = 'coll8b'
T0 = time.process_time()
W0 = time.time()


def check_clock(tag):
    cpu = (time.process_time() - T0) / 60.0
    wall = (time.time() - W0) / 60.0
    if wall > 1.5 * cpu + 2.0:
        print(f"    [clock] wall {wall:.1f} vs cpu {cpu:.1f} at {tag}", flush=True)
        os.system("pmset -g log 2>/dev/null | grep -i 'Entering Sleep' | tail -2")
    if cpu > STOP_MIN:
        raise SystemExit(f"STOP RULE: {cpu:.1f} min of process time at {tag}")
    return cpu, wall


def slug(cell):
    return cell.replace('/', '').replace(' ', '')


def gate(cell, seeds=(30, 31, 32), ms=2000.0):
    """Batched == separate, bit for bit, for THIS cell at the arms' level."""
    singles = [build(cell, s, g_ie=G_IE) for s in seeds]
    refs = [build(cell, s, g_ie=G_IE) for s in seeds]
    b = BatchedSTGDish(singles, list(stg.PRINZ_TABLE2[cell]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE, budget=BUDGET)
    b.run(ms, learn=True)
    for r in refs:
        r.run(ms, learn=True)
    same, dw = True, 0.0
    for k, r in enumerate(refs):
        for i in range(N_EXC):
            x = np.asarray(b.spikes_of(k, i), dtype=float)
            y = np.asarray(r.spikes(i), dtype=float)
            if x.shape != y.shape or not np.array_equal(x, y):
                same = False
        dw = max(dw, float(np.max(np.abs(b.weights_of(k) - r.weights))))
    return same, dw


def lot(cell, learn):
    tag = f"{slug(cell)}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'{TAG}_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t = time.process_time()
    dishes = [build(cell, s, g_ie=G_IE) for s in SEEDS]
    b = BatchedSTGDish(dishes, list(stg.PRINZ_TABLE2[cell]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE, budget=BUDGET)
    trace = [[] for _ in SEEDS]
    for c in range(int(round(T_MS / CHUNK_MS))):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn)
        for k in range(len(SEEDS)):
            spk = [np.asarray(b.spikes_of(k, i), dtype=float)
                   for i in range(N_EXC)]
            act = sum(1 for s in spk if ((s > lo) & (s <= lo + CHUNK_MS)).any())
            trace[k].append(dict(t_s=(lo + CHUNK_MS) / 1000.0, cells=act))
    rows = []
    for k, seed in enumerate(SEEDS):
        spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        rows.append(dict(cell=cell, g_ie=G_IE, learn=learn, seed=seed,
                         cells_tail=tail, alive=bool(tail >= Q),
                         alive_strict=bool(trace[k][-1]['cells'] >= Q),
                         bursts_tail=len(wins),
                         total_spikes=sum(int(s.size) for s in spk),
                         trace=trace[k]))
    res = dict(cell=cell, learn=learn, rows=rows,
               cpu_s=time.process_time() - t)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, wall = check_clock(tag)
    print(f"  {tag}: {sum(1 for r in rows if r['alive']):2d}/{len(SEEDS)} alive,"
          f" {res['cpu_s'] / 60:.1f} min cpu  [total {cpu:.1f} cpu / {wall:.1f}"
          f" wall of {STOP_MIN}]", flush=True)
    return res


def table(off, on, key='alive'):
    seeds = [s for s in sorted(off) if s in on]
    saved = sum(1 for s in seeds if not off[s][key] and on[s][key])
    killed = sum(1 for s in seeds if off[s][key] and not on[s][key])
    dd = sum(1 for s in seeds if not off[s][key] and not on[s][key])
    aa = len(seeds) - saved - killed - dd
    n = len(seeds)
    return dict(n=n, aa=aa, killed=killed, saved=saved, dd=dd,
                diff=(saved - killed) / n if n else float('nan'))


def main():
    print(f"COLL-8b: {CELLS}, {len(SEEDS)} wirings each, g_ie {G_IE}, "
          f"{T_MS / 1000:.0f} s\n", flush=True)
    print("  gates, per cell, at the arms' own level", flush=True)
    for cell in CELLS:
        done = all(os.path.exists(os.path.join(HERE, f'{TAG}_{slug(cell)}_{a}.json'))
                   for a in ('off', 'on'))
        if done:
            print(f"    {cell:>9}: both lots on disk, gate already paid", flush=True)
            continue
        same, dw = gate(cell)
        print(f"    {cell:>9}: spikes identical {same}, max |dw| {dw:.3e} "
              f"-> {'PASS' if same and dw == 0.0 else 'FAIL'}", flush=True)
        if not (same and dw == 0.0):
            raise SystemExit(f"GATE FAILED on {cell}; nothing else runs.")
    print(f"  [gates done at {(time.process_time() - T0) / 60:.1f} min]\n",
          flush=True)

    out = {}
    for cell in CELLS:
        off = {r['seed']: r for r in lot(cell, False)['rows']}
        on = {r['seed']: r for r in lot(cell, True)['rows']}
        out[cell] = (off, on)
        print(flush=True)

    print("=" * 74)
    print("COLL-8b  saved - killed on the two cells bracketing the crossing")
    print("=" * 74 + "\n")
    print(f"    {'cell':>9} {'OFF spikes':>11} {'n':>3} {'saved':>6} "
          f"{'killed':>7} {'diff':>7} {'strict diff':>12}")
    summ = {}
    for cell in CELLS:
        off, on = out[cell]
        t = table(off, on)
        ts = table(off, on, key='alive_strict')
        rate = float(np.mean([off[s]['total_spikes'] for s in off]))
        summ[cell] = dict(rate=rate, **t, strict_diff=ts['diff'])
        print(f"    {cell:>9} {rate:>11,.0f} {t['n']:>3} {t['saved']:>6} "
              f"{t['killed']:>7} {t['diff']:>+7.3f} {ts['diff']:>+12.3f}")
    json.dump(summ, open(os.path.join(HERE, f'{TAG}_summary.json'), 'w'),
              indent=1, default=float)

    print("\n  the pre-registered ordering used the COLL-4 scan at g_ie 0.030;")
    print("  the OFF rates above are this row's own, at the arms' level, free.")
    scan = {'LP 3': 6276, 'AB/PD 4': 8016}
    pre = sorted(CELLS, key=lambda c: scan[c])
    now = sorted(CELLS, key=lambda c: summ[c]['rate'])
    print(f"    pre-registered order: {pre}")
    print(f"    measured order:       {now}")
    print(f"    agree: {pre == now}")
    print(f"\n  total {(time.process_time() - T0) / 60:.1f} min cpu, "
          f"{(time.time() - W0) / 60:.1f} min wall of {STOP_MIN}")


if __name__ == '__main__':
    main()
