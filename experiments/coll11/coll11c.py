"""COLL-11c: H and leak, each at its own dying drive.

Registered in BRAIN-SIMULATION.md as `## Run: COLL-11c (registration)`,
before any arm, with 7.4 min of drive probing reported there.

These two swaps RAISE the rate and so lost their controls at the -0.03 the
other six arms ran at. Each runs here at its own drive, chosen from the probe
to sit below COLL-9's measured death band: H at -0.30, leak at -0.14.

QUALIFYING RULE: an OFF lot must yield at least 4 of 7 dead or its ON arm is
not run and the arm is reported unreachable. Fewer than 3 dead is
UNDERPOWERED and is not a null.

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
G_IE = 0.045
SEEDS = list(range(30, 37))          # 7 wirings
T_MS, CHUNK_MS, Q = 35000.0, 5000.0, 24
BUDGET = (N_EXC - 1) * 0.15 * 0.015
ARMS = [('H', 6, -0.30), ('leak', 7, -0.14)]   # per-arm drive
STOP_MIN = 21.5
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


def make(seeds, g_e, extra):
    dishes = [build(BASE, s, g_ie=G_IE) for s in seeds]
    b = BatchedSTGDish(dishes, list(g_e), list(stg.PRINZ_TABLE2['LP 3']),
                       g_ie=G_IE, budget=BUDGET)
    e = np.zeros(b.n)
    n1 = b.n_exc1 + b.n_inh1
    for k in range(len(seeds)):
        e[np.arange(b.n_exc1) + k * n1] = extra
    return b, e


def spikes_of_run(g_e, seeds, learn, ms, extra):
    b, e = make(seeds, g_e, extra)
    b.run(ms, learn=learn, extra=e)
    return [[np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
            for k in range(len(seeds))]


def lot(name, g_e, learn, extra):
    tag = f"{name}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll11c_{tag}.json')
    # NO legacy alias here. coll11b reused coll11a's CaS_off because it was
    # the same condition; in this row H and leak run at DIFFERENT drives
    # (-0.30 and -0.14 against -0.03), so reusing 11a's lots would silently
    # answer a different question. Third consumer of this class to bite.
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t = time.process_time()
    b, e = make(SEEDS, g_e, extra)
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
        rows.append(dict(arm=name, extra=extra, learn=learn, seed=seed,
                         cells_tail=tail,
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
    print(f"COLL-11c: {BASE} -> {TARGET}, arms {[(a, e) for a, _, e in ARMS]}, "
          f"{len(SEEDS)} wirings\n", flush=True)

    print(f"  gates were paid in COLL-11a (identity, null-swap, batch); "
          f"the swap machinery is unchanged here\n", flush=True)

    out = {}
    for name, idx, extra in ARMS:
        g_e = vec(idx)
        offlot = lot(name, g_e, False, extra)
        off = {r['seed']: r for r in offlot['rows']}
        n_dead = sum(1 for r in offlot['rows'] if not r['alive'])
        if n_dead < 4:
            print(f"  {name}: OFF gives {n_dead}/7 dead, below the registered 4/7 "
                  f"-> ON arm NOT run, arm reported unreachable\n", flush=True)
            out[name] = (off, None)
            continue
        print(f"  {name}: OFF gives {n_dead}/7 dead, qualifies", flush=True)
        on = {r['seed']: r for r in lot(name, g_e, True, extra)['rows']}
        out[name] = (off, on)
        print(flush=True)

    print("=" * 78)
    print("COLL-11c  rescue by single-conductance swap, against COLL-9's 0/21")
    print("=" * 78 + "\n")
    print(f"    {'arm':>6} {'sp/s OFF':>9} {'OFF dead':>9} {'rescued':>9} "
          f"{'killed':>7} {'Fisher vs 0/21':>15}")
    summ = {}
    for name, idx, extra in ARMS:
        if name not in out:
            continue
        off, on = out[name]
        if on is None:
            sl = sorted(off)
            nd = sum(1 for x in sl if not off[x]['alive'])
            print(f"    {name:>6} {np.mean([off[x]['rate_hz'] for x in sl]):>9.0f} "
                  f"{nd}/{len(sl):<7} {'unreachable -- ON arm not run':>30}")
            continue
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

    # ---- P29 needs all six swaps, so 11a's arms are read from disk -------
    import glob
    prev = {}
    for f in (glob.glob(os.path.join(HERE, 'coll11a_*.json'))
              + glob.glob(os.path.join(HERE, 'coll11b_*.json'))):
        if 'summary' in f:
            continue
        d = json.load(open(f))
        prev.setdefault(d['arm'], {})[d['learn']] = {r['seed']: r
                                                     for r in d['rows']}
    for name, lots in prev.items():
        if name in summ or len(lots) < 2:
            continue
        off, on = lots[False], lots[True]
        sl = sorted(off)
        dead = [x for x in sl if not off[x]['alive']]
        live = [x for x in sl if off[x]['alive']]
        k = sum(1 for x in dead if on[x]['alive'])
        summ[name] = dict(rate=float(np.mean([off[x]['rate_hz'] for x in sl])),
                          n_dead=len(dead), rescued=k,
                          killed=sum(1 for x in live if not on[x]['alive']),
                          p=fisher_ge(k, len(dead)) if dead else float('nan'),
                          underpowered=len(dead) < 3,
                          source=('COLL-11a' if 'coll11a' in f else 'COLL-11b'))

    print("\n  the whole single-swap grid, 11a and 11b together")
    print(f"    {'arm':>6} {'OFF sp/s':>9} {'OFF dead':>9} {'rescued':>9} "
          f"{'Fisher':>8}  where")
    for name in ('full', 'Na', 'CaS', 'KCa', 'Kd', 'H', 'leak'):
        if name not in summ:
            print(f"    {name:>6} {'not run':>9}")
            continue
        d = summ[name]
        src = d.get('source', 'COLL-11b')
        flag = ' UNDERPOWERED' if d['underpowered'] else ''
        print(f"    {name:>6} {d['rate']:>9.0f} {d['n_dead']:>9} "
              f"{d['rescued']:>9} {d['p']:>8.4f}  {src}{flag}")

    print("\n  P29  no single swap other than the full one reaches p < 0.05")
    singles = {n: d for n, d in summ.items()
               if n != 'full' and not d['underpowered']}
    hits = [n for n, d in singles.items() if d['rescued'] > 0 and d['p'] < 0.05]
    print(f"      testable single swaps: {sorted(singles)}")
    print(f"      reaching p < 0.05: {hits if hits else 'NONE'}")
    if not hits and singles:
        print("      P29 HOLDS so far -- rescuability is not carried by any "
              "single conductance tested")
    elif hits:
        print(f"      P29 REFUTED -- {hits} carries it")
    missing = [n for n in ('Na', 'CaS', 'KCa', 'Kd', 'H', 'leak')
               if n not in summ or summ[n]['underpowered']]
    if missing:
        print(f"      NOT settled: {missing} untested or underpowered; "
              f"no verdict is written for them")

    json.dump(summ, open(os.path.join(HERE, 'coll11c_summary.json'), 'w'),
              indent=1, default=float)
    print(f"\n  total {(time.process_time() - T0) / 60:.1f} min cpu, "
          f"{(time.time() - W0) / 60:.1f} min wall of {STOP_MIN}")


if __name__ == '__main__':
    main()
