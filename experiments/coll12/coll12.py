"""COLL-12: does the two-sign result survive initialisation and network size?

Registered in BRAIN-SIMULATION.md at 8ac4db4 before any lot, with the price
and the power arithmetic. The falsifier is a CROSSING between the two cells
within a condition, which is a two-sample question and needs ten
opportunities per cell rather than the eighteen a one-sample sign test would.

    families   A  w0 ~ w_rec * (1 + 0.2 * N(0,1)), clipped at 0 -- the
                  initialisation every previous row used
               B  w0 ~ Uniform(0, 2 * w_rec) -- same mean, so the homeostatic
                  budget starts identical; far wider spread
    sizes      48 excitatory + 12 inhibitory, and 24 + 6
    cells      PY 1 at its own drive (its control collapses unaided) and
               AB/PD 1 driven to -0.05, where COLL-9 measured 12/12

An arm yielding fewer than 3 collapsed controls is UNDERPOWERED and is not
reported as a rate; that rule is standing across this programme.
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
sys.path.insert(0, os.path.join(HERE, '..', 'coll5'))

from stgbatch import BatchedSTGDish
from luviner.biophysics import stg
from luviner.biophysics.stg_dish import STGDish, STGDriver
from stg_pm import CONN, W_REC, burst_windows

G_IE = 0.045
T_MS, CHUNK_MS = 35000.0, 5000.0
FAMILY = sys.argv[1] if len(sys.argv) > 1 else 'A'
STOP_MIN = 41.0
# cell -> (injected current, wirings)  : the fragility-matched pair
CELLS = [('PY 1', 0.0, 11), ('AB/PD 1', -0.05, 10)]
SIZES = [48, 24]
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


def build_sized(cell, seed, n_exc, family):
    """The dish of every previous row, with size and initialisation as
    parameters. Family B replaces w0 before the dish is batched, which is
    where BatchedSTGDish reads it."""
    n_inh = n_exc // 4
    budget = (n_exc - 1) * CONN * W_REC
    g_e, g_i = list(stg.PRINZ_TABLE2[cell]), list(stg.PRINZ_TABLE2['LP 3'])

    def factory(n):
        cells = ([stg.STGNeuron(list(g_e), name=f'E{i}') for i in range(n_exc)]
                 + [stg.STGNeuron(list(g_i), name=f'I{i}') for i in range(n_inh)])
        return STGDriver(stg.STGNetwork(cells, []))

    d = STGDish(factory, n_exc=n_exc, n_inh=n_inh, connectivity=CONN,
                w_rec=W_REC, g_ie=G_IE, budget=budget,
                depress=(0.20, 400.0), seed=seed)
    if family == 'B':
        rng = np.random.RandomState(10_000 + seed)
        d.w0 = np.clip(rng.uniform(0.0, 2.0 * W_REC, d.w0.size), 0.0, None)
    return d, budget


def lot(cell, extra, seeds, n_exc, family, learn):
    tag = f"{cell.replace('/','').replace(' ','')}_n{n_exc}_{family}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll12_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t = time.process_time()
    dishes, budget = [], None
    for s in seeds:
        d, budget = build_sized(cell, s, n_exc, family)
        dishes.append(d)
    b = BatchedSTGDish(dishes, list(stg.PRINZ_TABLE2[cell]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE, budget=budget)
    e = np.zeros(b.n)
    if extra:
        n1 = b.n_exc1 + b.n_inh1
        for k in range(len(seeds)):
            e[np.arange(b.n_exc1) + k * n1] = extra
    q = n_exc // 2
    trace = [[] for _ in seeds]
    for c in range(int(round(T_MS / CHUNK_MS))):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn, extra=(e if extra else None))
        for k in range(len(seeds)):
            spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(n_exc)]
            trace[k].append(sum(1 for s in spk
                                if ((s > lo) & (s <= lo + CHUNK_MS)).any()))
    rows = []
    for k, seed in enumerate(seeds):
        spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(n_exc)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, n_exc)
        rows.append(dict(cell=cell, extra=extra, family=family, n_exc=n_exc,
                         quorum=q, learn=learn, seed=seed, cells_tail=tail,
                         alive=bool(tail >= q),
                         alive_strict=bool(trace[k][-1] >= q),
                         n_bursts=len(wins),
                         rate_hz=sum(int(s.size) for s in spk) / (T_MS / 1000.0),
                         trace=trace[k]))
    res = dict(cell=cell, family=family, n_exc=n_exc, learn=learn,
               rows=rows, cpu_s=time.process_time() - t)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, _ = check_clock(tag)
    print(f"  {tag}: {sum(1 for r in rows if r['alive'])}/{len(seeds)} alive "
          f"(Q={q}), mean {np.mean([r['rate_hz'] for r in rows]):.0f} sp/s, "
          f"{res['cpu_s']/60:.1f} min  [total {cpu:.1f} of {STOP_MIN}]", flush=True)
    return res


def main():
    print(f"COLL-12, initialisation family {FAMILY}: sizes {SIZES}, "
          f"cells {[c for c,_,_ in CELLS]}, g_ie {G_IE}\n", flush=True)
    for n_exc in SIZES:
        for cell, extra, nw in CELLS:
            seeds = list(range(30, 30 + nw))
            for learn in (False, True):
                lot(cell, extra, seeds, n_exc, FAMILY, learn)
        print(flush=True)
    print(f"  total {(time.process_time()-T0)/60:.1f} min cpu, "
          f"{(time.time()-W0)/60:.1f} min wall of {STOP_MIN}", flush=True)


if __name__ == '__main__':
    main()
