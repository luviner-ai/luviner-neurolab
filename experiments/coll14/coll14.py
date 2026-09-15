"""COLL-14: PY 1 and AB/PD 2 re-run at 35 s, so the ordering has one horizon.

Registered in BRAIN-SIMULATION.md at 97f1254 before any lot, with the
consequence for PAPER-2's abstract of each outcome written there.

Only the horizon changes: same seeds, same wirings, same rule, same budget,
same g_ie, same quorum. Both detectors are recorded per wiring, which also
closes the outstanding AB/PD 2 strict-detector recheck (referee point D7).

LOTS is the unit of splitting: 14a runs the first two PY 1 lots, 14b the third
plus AB/PD 2. Each lot is written to disk before the clock check, so a stop
rule costs a lot and never a table.
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

G_IE = 0.045
T_MS, CHUNK_MS, Q = 35000.0, 5000.0, 24
BUDGET = (N_EXC - 1) * 0.15 * 0.015

HALF = sys.argv[1] if len(sys.argv) > 1 else 'a'
if HALF == 'a':
    LOTS = [('PY 1', list(range(30, 42))), ('PY 1', list(range(42, 54)))]
    STOP_MIN = 32.0
else:
    LOTS = [('PY 1', list(range(54, 66))), ('AB/PD 2', list(range(30, 42)))]
    STOP_MIN = 32.0

T0 = time.process_time()
W0 = time.time()


def check_clock(tag):
    cpu = (time.process_time() - T0) / 60.0
    wall = (time.time() - W0) / 60.0
    if wall > 1.5 * cpu + 2.0:
        print(f"    [clock] wall {wall:.1f} vs cpu {cpu:.1f} at {tag}", flush=True)
        os.system("pmset -g log 2>/dev/null | grep -i 'Entering Sleep' | tail -2")
    if cpu > STOP_MIN:
        raise SystemExit(f"STOP RULE: {cpu:.1f} min at {tag}")
    return cpu, wall


def slug(cell, seeds):
    return f"{cell.replace('/', '').replace(' ', '')}_{seeds[0]}"


def lot(cell, seeds, learn):
    tag = f"{slug(cell, seeds)}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll14_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t = time.process_time()
    dishes = [build(cell, s, g_ie=G_IE) for s in seeds]
    b = BatchedSTGDish(dishes, list(stg.PRINZ_TABLE2[cell]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE, budget=BUDGET)
    trace = [[] for _ in seeds]
    for c in range(int(round(T_MS / CHUNK_MS))):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn)
        for k in range(len(seeds)):
            spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
            trace[k].append(sum(1 for s in spk
                                if ((s > lo) & (s <= lo + CHUNK_MS)).any()))
    rows = []
    for k, seed in enumerate(seeds):
        spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        rows.append(dict(cell=cell, g_ie=G_IE, learn=learn, seed=seed,
                         horizon_s=T_MS / 1000.0, cells_tail=tail,
                         alive=bool(tail >= Q),
                         alive_strict=bool(trace[k][-1] >= Q),
                         n_bursts=len(wins),
                         total_spikes=sum(int(s.size) for s in spk),
                         trace=trace[k]))
    res = dict(cell=cell, learn=learn, seeds=seeds, rows=rows,
               cpu_s=time.process_time() - t)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, _ = check_clock(tag)
    print(f"  {tag}: {sum(1 for r in rows if r['alive'])}/{len(seeds)} alive "
          f"(strict {sum(1 for r in rows if r['alive_strict'])}), "
          f"{res['cpu_s'] / 60:.1f} min  [total {cpu:.1f} of {STOP_MIN}]",
          flush=True)
    return res


def main():
    print(f"COLL-14{HALF}: {[(c, s[0]) for c, s in LOTS]} at {T_MS/1000:.0f} s, "
          f"g_ie {G_IE}\n", flush=True)
    for cell, seeds in LOTS:
        for learn in (False, True):
            lot(cell, seeds, learn)
        print(flush=True)
    print(f"  total {(time.process_time() - T0) / 60:.1f} min cpu, "
          f"{(time.time() - W0) / 60:.1f} min wall of {STOP_MIN}", flush=True)


if __name__ == '__main__':
    main()
