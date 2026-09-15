"""COLL-11a: which single conductance carries the rescue?

Registered in BRAIN-SIMULATION.md at 319f060 before any arm. Arm 0 is the
full swap -- PY 1 by construction -- and is the positive control: if it does
not rescue under this protocol, nothing else in the row can be read.

An arm with fewer than 3 dead-OFF wirings is UNDERPOWERED for rescue and is
not reported as a null. Declared in the registration, not discovered here.
"""

import json
import os
import sys
import time

import numpy as np
from math import comb

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

BASE, TARGET = 'AB/PD 1', 'PY 1'
NAMES = stg.CONDUCTANCE_NAMES
G_IE, EXTRA = 0.045, -0.03
SEEDS = list(range(30, 37))          # 7 wirings
T_MS, CHUNK_MS, Q = 35000.0, 5000.0, 24
BUDGET = (N_EXC - 1) * 0.15 * 0.015
ARMS = [('full', None), ('H', 6), ('KCa', 4), ('CaS', 2)]
STOP_MIN = 41.0
T0 = time.process_time()
W0 = time.time()


def vec(idx):
    """AB/PD 1 with conductance `idx` taken from PY 1; None means all of them."""
    g = list(stg.PRINZ_TABLE2[BASE])
    t = list(stg.PRINZ_TABLE2[TARGET])
    if idx is None:
        return list(t)
    g[idx] = t[idx]
    return g


def check_clock(tag):
    cpu = (time.process_time() - T0) / 60.0
    wall = (time.time() - W0) / 60.0
    if wall > 1.5 * cpu + 2.0:
        print(f"    [clock] wall {wall:.1f} vs cpu {cpu:.1f} at {tag}", flush=True)
        os.system("pmset -g log 2>/dev/null | grep -i 'Entering Sleep' | tail -2")
    if cpu > STOP_MIN:
        raise SystemExit(f"STOP RULE: {cpu:.1f} min at {tag}")
    return cpu, wall


def make(seeds, g_e):
    dishes = [build(BASE, s, g_ie=G_IE) for s in seeds]
    b = BatchedSTGDish(dishes, list(g_e), list(stg.PRINZ_TABLE2['LP 3']),
                       g_ie=G_IE, budget=BUDGET)
    e = np.zeros(b.n)
    n1 = b.n_exc1 + b.n_inh1
    for k in range(len(seeds)):
        e[np.arange(b.n_exc1) + k * n1] = EXTRA
    return b, e


def spikes_of_run(g_e, seeds, learn, ms):
    b, e = make(seeds, g_e)
    b.run(ms, learn=learn, extra=e)
    return [[np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
            for k in range(len(seeds))]


def lot(name, g_e, learn):
    tag = f"{name}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll11a_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t = time.process_time()
    b, e = make(SEEDS, g_e)
    trace = [[] for _ in SEEDS]
    for c in range(int(round(T_MS / CHUNK_MS))):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn, extra=e)
        for k in range(len(SEEDS)):
            spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
            trace[k].append(sum(1 for s in spk
                                if ((s > lo) & (s <= lo + CHUNK_MS)).any()))
    rows = []
    for k, seed in enumerate(SEEDS):
        spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        rows.append(dict(arm=name, learn=learn, seed=seed, cells_tail=tail,
                         alive=bool(tail >= Q),
                         alive_strict=bool(trace[k][-1] >= Q),
                         n_bursts=len(wins),
                         rate_hz=sum(int(s.size) for s in spk) / (T_MS / 1000.0),
                         trace=trace[k]))
    res = dict(arm=name, learn=learn, rows=rows, cpu_s=time.process_time() - t)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, _ = check_clock(tag)
    print(f"  {tag}: {sum(1 for r in rows if r['alive'])}/{len(SEEDS)} alive, "
          f"mean {np.mean([r['rate_hz'] for r in rows]):.0f} sp/s, "
          f"{res['cpu_s'] / 60:.1f} min  [total {cpu:.1f} of {STOP_MIN}]", flush=True)
    return res


def fisher_ge(k, m, c=0, d=21):
    R, N = k + c, m + d
    return sum(comb(m, i) * comb(d, R - i) / comb(N, R)
               for i in range(k, min(m, R) + 1))


def main():
    print(f"COLL-11a: {BASE} -> {TARGET}, arms {[a for a, _ in ARMS]}, "
          f"{len(SEEDS)} wirings, extra {EXTRA}\n", flush=True)

    # ---- identity gate: the full swap IS PY 1 -----------------------------
    same_vec = vec(None) == list(stg.PRINZ_TABLE2[TARGET])
    print(f"  gate identity: full-swap vector == PY 1: {same_vec}", flush=True)
    a = spikes_of_run(vec(None), [30], True, 2000.0)
    b2 = BatchedSTGDish([build(TARGET, 30, g_ie=G_IE)],
                        list(stg.PRINZ_TABLE2[TARGET]),
                        list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE,
                        budget=BUDGET)
    e2 = np.zeros(b2.n)
    e2[np.arange(b2.n_exc1)] = EXTRA
    b2.run(2000.0, learn=True, extra=e2)
    ok_id = all(np.array_equal(a[0][i],
                               np.asarray(b2.spikes_of(0, i), dtype=float))
                for i in range(N_EXC))
    print(f"  gate identity: full swap == a directly built PY 1 dish: {ok_id}",
          flush=True)

    # ---- null-swap gate: CaT is identical, so it must change nothing ------
    x = spikes_of_run(vec(1), [30], True, 2000.0)
    y = spikes_of_run(list(stg.PRINZ_TABLE2[BASE]), [30], True, 2000.0)
    ok_null = all(np.array_equal(x[0][i], y[0][i]) for i in range(N_EXC))
    print(f"  gate null-swap (CaT): identical to unswapped {BASE}: {ok_null}",
          flush=True)
    if not (same_vec and ok_id and ok_null):
        raise SystemExit("GATE FAILED; nothing else runs.")
    print(f"  [gates at {(time.process_time() - T0) / 60:.1f} min]\n", flush=True)

    out = {}
    for name, idx in ARMS:
        g_e = vec(idx)
        off = {r['seed']: r for r in lot(name, g_e, False)['rows']}
        on = {r['seed']: r for r in lot(name, g_e, True)['rows']}
        out[name] = (off, on)
        print(flush=True)

    print("=" * 78)
    print("COLL-11a  rescue by single-conductance swap, against COLL-9's 0/21")
    print("=" * 78 + "\n")
    print(f"    {'arm':>6} {'sp/s OFF':>9} {'OFF dead':>9} {'rescued':>9} "
          f"{'killed':>7} {'Fisher vs 0/21':>15}")
    summ = {}
    for name, idx in ARMS:
        if name not in out:
            continue
        off, on = out[name]
        s = sorted(off)
        dead = [x for x in s if not off[x]['alive']]
        live = [x for x in s if off[x]['alive']]
        k = sum(1 for x in dead if on[x]['alive'])
        killed = sum(1 for x in live if not on[x]['alive'])
        rate = float(np.mean([off[x]['rate_hz'] for x in s]))
        p = fisher_ge(k, len(dead)) if dead else float('nan')
        under = len(dead) < 3
        summ[name] = dict(rate=rate, n_dead=len(dead), rescued=k, killed=killed,
                          p=p, underpowered=under)
        flag = '  UNDERPOWERED' if under else ''
        print(f"    {name:>6} {rate:>9.0f} {len(dead)}/{len(s):<7} "
              f"{k}/{len(dead) if dead else 0:<7} {killed:>7} "
              f"{p:>15.4f}{flag}")

    print("\n  the positive control")
    if 'full' in summ:
        f = summ['full']
        if f['underpowered']:
            print("      full swap UNDERPOWERED -- the row cannot be read")
        elif f['rescued'] > 0 and f['p'] < 0.05:
            print(f"      full swap rescues {f['rescued']}/{f['n_dead']}, "
                  f"p = {f['p']:.4f} -- the protocol CAN show a rescue, "
                  f"the swaps are readable")
        else:
            print(f"      full swap rescues {f['rescued']}/{f['n_dead']} -- "
                  f"**THE ROW IS VOID**: the protocol cannot show a rescue "
                  f"even with PY 1's own conductances")

    print("\n  P27 (mine, H) and P28 (mine, KCa) against luviner-55 (burst duration: CaS or KCa)")
    for name in ('H', 'KCa', 'CaS'):
        if name in summ:
            d = summ[name]
            v = ('UNDERPOWERED' if d['underpowered']
                 else ('RESCUES' if d['rescued'] > 0 and d['p'] < 0.05 else 'null'))
            print(f"      {name:>4}: {d['rescued']}/{d['n_dead']} rescued, "
                  f"p = {d['p']:.4f}  -> {v}")
    winners = [n for n, d in summ.items()
               if n != 'full' and not d['underpowered']
               and d['rescued'] > 0 and d['p'] < 0.05]
    print(f"\n      arms that rescue: {winners if winners else 'NONE in this half'}")
    if not winners and 'full' in summ and summ['full']['rescued'] > 0:
        print("      -> so far consistent with: no SINGLE conductance carries it")

    json.dump(summ, open(os.path.join(HERE, 'coll11a_summary.json'), 'w'),
              indent=1, default=float)
    print(f"\n  total {(time.process_time() - T0) / 60:.1f} min cpu, "
          f"{(time.time() - W0) / 60:.1f} min wall of {STOP_MIN}")


if __name__ == '__main__':
    main()
