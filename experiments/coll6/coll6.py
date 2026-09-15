"""COLL-6: does the occasion-conditioned rescue survive a frequent occasion?

Registered in BRAIN-SIMULATION.md, `## Run: COLL-6 (registration)`, commit
033aecb, before any wiring ran. Stage 0 (batched == separate, bit for bit,
on PY 1 at the arms' own level) is `gate.py` and passed at 0.000e+00.

One level, 36 unselected wirings, three lots of K = 12 -- COLL-5's measured
lot size, so its cost constant applies unchanged. Lot order pairs the arms
so a stop rule costs whole wiring-pairs and never half a 2x2 table, with the
OFF gate lot first.
"""

import json
import os
import sys
import time
from math import comb

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

CELL = 'PY 1'
G_IE = 0.045
LOTS = [list(range(30, 42)), list(range(42, 54)), list(range(54, 66))]
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
              f"checking pmset", flush=True)
        os.system("pmset -g log 2>/dev/null | grep -i 'Entering Sleep' | tail -2")
    if cpu > STOP_MIN:
        raise SystemExit(f"STOP RULE: {cpu:.1f} min of process time at {tag}")
    return cpu, wall


def lot(seeds, learn, li):
    tag = f"lot{li}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll6_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t, w = time.process_time(), time.time()
    dishes = [build(CELL, s, g_ie=G_IE) for s in seeds]
    b = BatchedSTGDish(dishes, list(stg.PRINZ_TABLE2[CELL]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE,
                       budget=BUDGET)
    trace = [[] for _ in seeds]
    for c in range(int(round(T_MS / CHUNK_MS))):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn)
        for k in range(len(seeds)):
            spk = [np.asarray(b.spikes_of(k, i), dtype=float)
                   for i in range(N_EXC)]
            act = sum(1 for s in spk if ((s > lo) & (s <= lo + CHUNK_MS)).any())
            trace[k].append(dict(t_s=(lo + CHUNK_MS) / 1000.0, cells=act))
    rows = []
    for k, seed in enumerate(seeds):
        spk = [np.asarray(b.spikes_of(k, i), dtype=float)
               for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        wt = b.weights_of(k)
        # first chunk that falls under quorum -- the graded severity P13 uses
        under = [d['t_s'] for d in trace[k] if d['cells'] < Q]
        rows.append(dict(cell=CELL, g_ie=G_IE, learn=learn, seed=seed,
                         cells_tail=tail, alive=bool(tail >= Q),
                         t_under_s=(under[0] if under else None),
                         bursts_tail=len(wins),
                         total_spikes=sum(int(s.size) for s in spk),
                         tot_w=float(np.mean(np.bincount(
                             dishes[k].post_e, weights=wt, minlength=N_EXC))),
                         budget=BUDGET, trace=trace[k]))
    n_alive = sum(1 for r in rows if r['alive'])
    res = dict(seeds=seeds, learn=learn, rows=rows,
               cpu_s=time.process_time() - t, wall_s=time.time() - w)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, wall = check_clock(tag)
    print(f"  {tag}: {n_alive:2d}/{len(seeds)} alive, "
          f"{res['cpu_s'] / 60:.1f} min cpu "
          f"[total {cpu:.1f} cpu / {wall:.1f} wall of {STOP_MIN}]", flush=True)
    return res


def wilson(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def binom_p_ge(k, n, p0=0.5):
    """One-sided exact binomial: P(X >= k | p = p0)."""
    return sum(comb(n, i) * p0 ** i * (1 - p0) ** (n - i) for i in range(k, n + 1))


def main():
    print(f"COLL-6: {CELL} at g_ie {G_IE}, {sum(len(l) for l in LOTS)} "
          f"unselected wirings, {T_MS / 1000:.0f} s, three lots of "
          f"K={len(LOTS[0])}\n", flush=True)
    off, on = {}, {}
    for li, seeds in enumerate(LOTS, 1):
        for learn, store in ((False, off), (True, on)):
            r = lot(seeds, learn, li)
            for row in r['rows']:
                store[row['seed']] = row
        print(flush=True)

    seeds = [s for s in sorted(off) if s in on]
    aa = [s for s in seeds if off[s]['alive'] and on[s]['alive']]
    ad = [s for s in seeds if off[s]['alive'] and not on[s]['alive']]
    da = [s for s in seeds if not off[s]['alive'] and on[s]['alive']]
    dd = [s for s in seeds if not off[s]['alive'] and not on[s]['alive']]
    n = len(seeds)
    n_dead_off, n_alive_off = len(da) + len(dd), len(aa) + len(ad)

    print("=" * 74)
    print(f"COLL-6  the 2x2 table, {n} wirings at g_ie {G_IE}")
    print("=" * 74)
    print(f"\n                    ON alive   ON dead")
    print(f"    OFF alive       {len(aa):8d}  {len(ad):8d}")
    print(f"    OFF dead        {len(da):8d}  {len(dd):8d}")
    print(f"\n    OFF-death fraction {n_dead_off}/{n} = "
          f"{n_dead_off / n if n else float('nan'):.3f}")

    if n_dead_off < 12:
        print(f"\n  *** UNDERPOWERED: {n_dead_off} dead-OFF cases against the "
              f"12 registered as the floor.")
        print(f"      The interval is reported; no wirings are added.")

    print("\n" + "-" * 74)
    print("The registered quantity: P(ON alive | OFF dead)")
    if n_dead_off:
        k = len(da)
        p = k / n_dead_off
        lo, hi = wilson(k, n_dead_off)
        pv = binom_p_ge(k, n_dead_off, 0.5)
        rej = pv <= 0.05
        print(f"      {k}/{n_dead_off} = {p:.3f}   Wilson 95% [{lo:.3f}, {hi:.3f}]")
        print(f"      exact one-sided binomial against p <= 0.5: p = {pv:.4f}"
              f"  -> {'rejects' if rej else 'does NOT reject'}")
        if n_dead_off < 12:
            verdict = 'UNDERPOWERED'
        elif rej and p >= 0.8:
            verdict = 'CONFIRMED -- rescue is real and high'
        elif rej:
            verdict = ('MIDDLE -- rescue is real and a majority; ">= 0.8" is '
                       'NOT established')
        else:
            verdict = 'FALSIFIED -- not established above chance'
        print(f"\n      REGISTERED VERDICT: {verdict}")
        print(f"      P12 (mine, point estimate in [0.45, 0.70]): "
              f"{'HOLDS' if 0.45 <= p <= 0.70 else 'REFUTED'}")
        print(f"      luviner-55 (point >= 0.8 and Wilson lower >= 0.6): "
              f"{'HOLDS' if p >= 0.8 and lo >= 0.6 else 'REFUTED'}")

    print("\nThe companion: P(ON dead | OFF alive)")
    if n_alive_off:
        k2 = len(ad)
        p2 = k2 / n_alive_off
        lo2, hi2 = wilson(k2, n_alive_off)
        print(f"      {k2}/{n_alive_off} = {p2:.3f}   Wilson 95% "
              f"[{lo2:.3f}, {hi2:.3f}]")
        print(f"      P14 (both of us, in [0.2, 0.5]): "
              f"{'HOLDS' if 0.2 <= p2 <= 0.5 else 'REFUTED'}")

    print("\n" + "-" * 74)
    print("P13  does severity grade the rescue? (free, from the OFF trace)")
    dead = [s for s in seeds if not off[s]['alive']]
    timed = [(off[s]['t_under_s'] if off[s]['t_under_s'] is not None else 1e9, s)
             for s in dead]
    timed.sort()
    if len(timed) >= 6:
        t3 = len(timed) // 3
        early = [s for _, s in timed[:t3]]
        late = [s for _, s in timed[-t3:]]
        re_ = sum(1 for s in early if on[s]['alive']) / len(early)
        rl = sum(1 for s in late if on[s]['alive']) / len(late)
        print(f"      earliest-dying third (n={len(early)}): rescued {re_:.3f}")
        print(f"      latest-dying third   (n={len(late)}): rescued {rl:.3f}")
        print(f"      difference {rl - re_:+.3f} "
              f"(predicted the early third at least 0.20 lower)")
        print(f"      P13: {'HOLDS' if rl - re_ >= 0.20 else 'REFUTED'}")
    else:
        print(f"      only {len(timed)} dead-OFF wirings; terciles need 6. "
              f"Not evaluated.")

    json.dump(dict(n=n, alive_alive=len(aa), killed=len(ad), saved=len(da),
                   dead_dead=len(dd), seeds_saved=da, seeds_dead_dead=dd,
                   seeds_killed=ad),
              open(os.path.join(HERE, 'coll6_summary.json'), 'w'), indent=1)
    cpu, wall = (time.process_time() - T0) / 60, (time.time() - W0) / 60
    print(f"\n  total {cpu:.1f} min cpu, {wall:.1f} min wall of {STOP_MIN}")


if __name__ == '__main__':
    main()
