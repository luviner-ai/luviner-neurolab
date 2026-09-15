"""COLL-5: the fraction of unselected wirings that learning saves, by g_ie.

Registered in BRAIN-SIMULATION.md, `## Run: COLL-5 (registration and gate)`,
commit 016f409, before any arm ran. Stage 0 (batched == separate, bit for
bit) is `gate.py` and passed at 0.000e+00 on both arms and both levels.

Twelve unselected wirings, three inhibition levels, two arms, 30 s, batched
K = 12. Lot order puts the two REGISTERED levels before the bonus one, so a
stop rule can only cost the bonus.
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

from price import build
from stgbatch import BatchedSTGDish
from luviner.biophysics import stg
from stg_pm import N_EXC, burst_windows

CELL = 'AB/PD 2'
SEEDS = list(range(30, 42))
LEVELS = [0.030, 0.045, 0.060]      # registered two first, bonus last
T_MS = 30000.0
CHUNK_MS = 5000.0
Q = 24
BUDGET = (N_EXC - 1) * 0.15 * 0.015
STOP_MIN = 42.0
T0 = time.process_time()
W0 = time.time()


def check_clock(tag):
    cpu = (time.process_time() - T0) / 60.0
    wall = (time.time() - W0) / 60.0
    if wall > 1.5 * cpu + 2.0:
        print(f"    [clock] wall {wall:.1f} min vs cpu {cpu:.1f} at {tag} -- "
              f"checking pmset")
        os.system("pmset -g log 2>/dev/null | grep -i 'Entering Sleep' "
                  "| tail -2")
    if cpu > STOP_MIN:
        raise SystemExit(f"STOP RULE: {cpu:.1f} min of process time at {tag}")
    return cpu, wall


def lot(g_ie, learn):
    tag = f"g{g_ie:.3f}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll5_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped")
        return json.load(open(path))
    check_clock(tag)
    t, w = time.process_time(), time.time()
    dishes = [build(CELL, s, g_ie=g_ie) for s in SEEDS]
    b = BatchedSTGDish(dishes, list(stg.PRINZ_TABLE2[CELL]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=g_ie,
                       budget=BUDGET)
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
        spk = [np.asarray(b.spikes_of(k, i), dtype=float)
               for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        wt = b.weights_of(k)
        rows.append(dict(cell=CELL, g_ie=g_ie, learn=learn, seed=seed,
                         cells_tail=tail, alive=bool(tail >= Q),
                         bursts_tail=len(wins),
                         total_spikes=sum(int(s.size) for s in spk),
                         tot_w=float(np.mean(np.bincount(
                             dishes[k].post_e, weights=wt, minlength=N_EXC))),
                         budget=BUDGET, trace=trace[k]))
    n_alive = sum(1 for r in rows if r['alive'])
    res = dict(g_ie=g_ie, learn=learn, rows=rows,
               cpu_s=time.process_time() - t, wall_s=time.time() - w)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, wall = check_clock(tag)
    print(f"  {tag}: {n_alive:2d}/12 alive, {res['cpu_s'] / 60:.1f} min cpu "
          f"[total {cpu:.1f} cpu / {wall:.1f} wall of {STOP_MIN}]")
    return res


def wilson(k, n, z=1.96):
    """Wilson interval, which behaves at 0 and n where the normal one does not."""
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    print(f"COLL-5: {CELL}, {len(SEEDS)} unselected wirings, "
          f"levels {LEVELS}, {T_MS / 1000:.0f} s, batched K={len(SEEDS)}\n")
    out = {}
    for g in LEVELS:
        for learn in (False, True):
            out[(g, learn)] = lot(g, learn)
        print()

    print("=" * 74)
    print("The 2x2 tables, per level")
    print("=" * 74)
    summ = {}
    for g in LEVELS:
        off = {r['seed']: r['alive'] for r in out[(g, False)]['rows']}
        on = {r['seed']: r['alive'] for r in out[(g, True)]['rows']}
        aa = sum(1 for s in SEEDS if off[s] and on[s])
        ad = sum(1 for s in SEEDS if off[s] and not on[s])      # killed
        da = sum(1 for s in SEEDS if not off[s] and on[s])      # saved
        dd = sum(1 for s in SEEDS if not off[s] and not on[s])
        n = len(SEEDS)
        killed, saved = ad / n, da / n
        lo_k, hi_k = wilson(ad, n)
        lo_s, hi_s = wilson(da, n)
        off_dead = (da + dd) / n
        summ[g] = dict(alive_alive=aa, killed=ad, saved=da, dead_dead=dd,
                       frac_killed=killed, frac_saved=saved,
                       diff=saved - killed, off_dead=off_dead,
                       ci_killed=[lo_k, hi_k], ci_saved=[lo_s, hi_s])
        print(f"\n  g_ie = {g:.3f}          ON alive   ON dead")
        print(f"    OFF alive              {aa:>6}    {ad:>6}   <- killed")
        print(f"    OFF dead               {da:>6}    {dd:>6}")
        print(f"                           ^ saved")
        print(f"    saved  {da}/12 = {saved:.3f}  (Wilson 95% "
              f"[{lo_s:.3f}, {hi_s:.3f}])")
        print(f"    killed {ad}/12 = {killed:.3f}  (Wilson 95% "
              f"[{lo_k:.3f}, {hi_k:.3f}])")
        print(f"    saved - killed = {saved - killed:+.3f}")
        print(f"    OFF arm death rate = {off_dead:.3f}  "
              f"(the ceiling on `saved`)")

    print("\n" + "=" * 74)
    print("saved - killed against g_ie")
    print("=" * 74 + "\n")
    print(f"    {'g_ie':>6} {'saved':>7} {'killed':>7} {'diff':>8} "
          f"{'OFF dead':>9}")
    for g in LEVELS:
        s = summ[g]
        print(f"    {g:>6.3f} {s['frac_saved']:>7.3f} {s['frac_killed']:>7.3f} "
              f"{s['diff']:>+8.3f} {s['off_dead']:>9.3f}")
    flip = any(summ[g]['diff'] > 0 for g in LEVELS)
    print(f"\n    orchestrator (sign flips, positive somewhere): "
          f"{'HOLDS' if flip else 'REFUTED'}")
    print(f"    mine (no positive value at any level): "
          f"{'HOLDS' if not flip else 'REFUTED'}")
    bound_ok = all(summ[g]['frac_saved'] <= summ[g]['off_dead'] + 1e-12
                   for g in LEVELS)
    print(f"    the bound (saved <= OFF death rate): "
          f"{'holds, as it must' if bound_ok else 'VIOLATED -- instrument bug'}")
    json.dump({str(g): summ[g] for g in summ},
              open(os.path.join(HERE, 'coll5_summary.json'), 'w'), indent=1)
    print(f"\n  total {(time.process_time() - T0) / 60:.1f} min cpu, "
          f"{(time.time() - W0) / 60:.1f} min wall")


if __name__ == '__main__':
    main()
