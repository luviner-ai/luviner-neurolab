"""COLL-9: AB/PD 1 driven down into COLL-8's gap.

Registered in BRAIN-SIMULATION.md, `## Run: COLL-9 (registration)`, commit
6536bfd, before any arm ran. The calibration there found the drive to be
bistable rather than gradable, so the settings are the two reachable states
inside the gap, and the primary statistic is the OFF-death rate rather than
the kill rate on healthy wirings, which is 0.000 at both ends of the gap.

Written out rather than derived. `extra` is a parameter of the shipped
`STGDish.run` that no previous row used, so the gate covers it.
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

CELL = 'AB/PD 1'
G_IE = 0.045
SETTINGS = [-0.05, -0.03]            # ~82 spikes/s, and the bimodal branch point
SEEDS = list(range(30, 42))          # 12 wirings per setting
T_MS = 35000.0
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
        print(f"    [clock] wall {wall:.1f} vs cpu {cpu:.1f} at {tag}", flush=True)
        os.system("pmset -g log 2>/dev/null | grep -i 'Entering Sleep' | tail -2")
    if cpu > STOP_MIN:
        raise SystemExit(f"STOP RULE: {cpu:.1f} min of process time at {tag}")
    return cpu, wall


def slug(ex):
    return f"m{abs(ex):.2f}".replace('.', '')


def make(seeds, ex):
    dishes = [build(CELL, s, g_ie=G_IE) for s in seeds]
    b = BatchedSTGDish(dishes, list(stg.PRINZ_TABLE2[CELL]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=G_IE, budget=BUDGET)
    e = np.zeros(b.n)
    n1 = b.n_exc1 + b.n_inh1
    for k in range(len(seeds)):
        e[np.arange(b.n_exc1) + k * n1] = ex
    return b, e


def gate(ex, seeds=(30, 31, 32), ms=2000.0):
    """Bit for bit WITH `extra` -- the parameter, not just the batching."""
    b, e = make(list(seeds), ex)
    refs = [build(CELL, s, g_ie=G_IE) for s in seeds]
    b.run(ms, learn=True, extra=e)
    single = np.full(refs[0].n, 0.0)
    single[np.arange(refs[0].n_exc)] = ex
    for r in refs:
        r.run(ms, learn=True, extra=single)
    same, dw = True, 0.0
    for k, r in enumerate(refs):
        for i in range(N_EXC):
            x = np.asarray(b.spikes_of(k, i), dtype=float)
            y = np.asarray(r.spikes(i), dtype=float)
            if x.shape != y.shape or not np.array_equal(x, y):
                same = False
        dw = max(dw, float(np.max(np.abs(b.weights_of(k) - r.weights))))
    return same, dw


def lot(ex, learn):
    tag = f"{slug(ex)}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll9_{tag}.json')
    if os.path.exists(path):
        print(f"  {tag}: on disk, skipped", flush=True)
        return json.load(open(path))
    check_clock(tag)
    t = time.process_time()
    b, e = make(SEEDS, ex)
    trace = [[] for _ in SEEDS]
    for c in range(int(round(T_MS / CHUNK_MS))):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn, extra=e)
        for k in range(len(SEEDS)):
            spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
            act = sum(1 for s in spk if ((s > lo) & (s <= lo + CHUNK_MS)).any())
            trace[k].append(dict(t_s=(lo + CHUNK_MS) / 1000.0, cells=act))
    rows = []
    for k, seed in enumerate(SEEDS):
        spk = [np.asarray(b.spikes_of(k, i), dtype=float) for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        # burst structure -- the risk luviner-55 declared, reported per wiring
        durs = [float(b1 - a1) for a1, b1 in wins] if wins else []
        per = ([float(wins[i + 1][0] - wins[i][0]) for i in range(len(wins) - 1)]
               if len(wins) > 1 else [])
        rows.append(dict(cell=CELL, extra=ex, g_ie=G_IE, learn=learn, seed=seed,
                         cells_tail=tail, alive=bool(tail >= Q),
                         alive_strict=bool(trace[k][-1]['cells'] >= Q),
                         n_bursts=len(wins),
                         burst_ms=float(np.mean(durs)) if durs else None,
                         period_ms=float(np.mean(per)) if per else None,
                         duty=(float(np.mean(durs) / np.mean(per))
                               if durs and per else None),
                         rate_hz=sum(int(s.size) for s in spk) / (T_MS / 1000.0),
                         total_spikes=sum(int(s.size) for s in spk),
                         trace=trace[k]))
    res = dict(extra=ex, learn=learn, rows=rows, cpu_s=time.process_time() - t)
    json.dump(res, open(path, 'w'), indent=1)
    cpu, wall = check_clock(tag)
    print(f"  {tag}: {sum(1 for r in rows if r['alive']):2d}/{len(SEEDS)} alive, "
          f"mean {np.mean([r['rate_hz'] for r in rows]):.0f} sp/s, "
          f"{res['cpu_s'] / 60:.1f} min cpu  [total {cpu:.1f} of {STOP_MIN}]",
          flush=True)
    return res


def wilson(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    print(f"COLL-9: {CELL} at g_ie {G_IE}, extra {SETTINGS}, "
          f"{len(SEEDS)} wirings each, {T_MS / 1000:.0f} s\n", flush=True)
    for ex in SETTINGS:
        done = all(os.path.exists(os.path.join(HERE, f'coll9_{slug(ex)}_{a}.json'))
                   for a in ('off', 'on'))
        if done:
            print(f"    extra {ex}: both lots on disk, gate already paid", flush=True)
            continue
        same, dw = gate(ex)
        print(f"    gate extra {ex}: identical {same}, max |dw| {dw:.3e} -> "
              f"{'PASS' if same and dw == 0.0 else 'FAIL'}", flush=True)
        if not (same and dw == 0.0):
            raise SystemExit(f"GATE FAILED at extra {ex}; nothing else runs.")
    print(f"  [gates done at {(time.process_time() - T0) / 60:.1f} min]\n", flush=True)

    out = {}
    for ex in SETTINGS:
        off = {r['seed']: r for r in lot(ex, False)['rows']}
        on = {r['seed']: r for r in lot(ex, True)['rows']}
        out[ex] = (off, on)
        print(flush=True)

    print("=" * 78)
    print("COLL-9  does driving AB/PD 1 into the gap make it fragile?")
    print("=" * 78 + "\n")
    print(f"    {'extra':>6} {'OFF sp/s':>9} {'n':>3} {'OFF dead':>9} "
          f"{'Wilson 95%':>18} {'saved':>6} {'killed':>7}")
    summ = {}
    rows_all = []
    for ex in SETTINGS:
        off, on = out[ex]
        s = sorted(off)
        dead = [x for x in s if not off[x]['alive']]
        live = [x for x in s if off[x]['alive']]
        saved = sum(1 for x in dead if on[x]['alive'])
        killed = sum(1 for x in live if not on[x]['alive'])
        lo, hi = wilson(len(dead), len(s))
        rate = float(np.mean([off[x]['rate_hz'] for x in s]))
        summ[ex] = dict(rate=rate, n=len(s), n_dead=len(dead), saved=saved,
                        killed=killed, ci=[lo, hi])
        print(f"    {ex:>6.2f} {rate:>9.0f} {len(s):>3} "
              f"{len(dead)}/{len(s):<7} [{lo:.3f}, {hi:.3f}]{'':>4} "
              f"{saved:>6} {killed:>7}")
        rows_all += [(off[x]['rate_hz'], not off[x]['alive']) for x in s]

    print("\n  P17  OFF-death above 0.25 at the ~82 sp/s setting")
    ex0 = SETTINGS[0]
    f0 = summ[ex0]['n_dead'] / summ[ex0]['n']
    print(f"      extra {ex0}: OFF-death {f0:.3f}")
    print(f"      P17: {'HOLDS' if f0 > 0.25 else 'REFUTED'}")
    print(f"      F17 (0 dead -> rate is a proxy for something else): "
          f"{'FIRES' if summ[ex0]['n_dead'] == 0 else 'does not fire'}")

    print("\n  P18  rescue where the opportunity exists")
    for ex in SETTINGS:
        d = summ[ex]
        if d['n_dead']:
            lo, hi = wilson(d['saved'], d['n_dead'])
            print(f"      extra {ex:>6.2f}: {d['saved']}/{d['n_dead']} = "
                  f"{d['saved'] / d['n_dead']:.3f}  Wilson [{lo:.3f}, {hi:.3f}]"
                  f"   (PY 1 was 33/33)")
        else:
            print(f"      extra {ex:>6.2f}: no dead OFF wirings, not defined")

    print("\n  P19  within the bimodal setting, does the low branch die more?")
    ex1 = SETTINGS[1]
    off1 = out[ex1][0]
    rates = sorted(off1[x]['rate_hz'] for x in off1)
    cut = (min(rates) + max(rates)) / 2
    low = [x for x in off1 if off1[x]['rate_hz'] < cut]
    high = [x for x in off1 if off1[x]['rate_hz'] >= cut]
    print(f"      branches split at {cut:.0f} sp/s: "
          f"low n={len(low)} (mean {np.mean([off1[x]['rate_hz'] for x in low]):.0f}), "
          f"high n={len(high)} (mean {np.mean([off1[x]['rate_hz'] for x in high]):.0f})"
          if low and high else "      no split: all wirings on one branch")
    if low and high:
        dl = sum(1 for x in low if not off1[x]['alive']) / len(low)
        dh = sum(1 for x in high if not off1[x]['alive']) / len(high)
        print(f"      OFF-death low {dl:.3f} vs high {dh:.3f}, difference {dl - dh:+.3f}")
        print(f"      P19: {'HOLDS' if dl > dh else 'REFUTED'}")

    print("\n  burst structure, the declared risk -- OFF arm")
    print(f"    {'extra':>6} {'sp/s':>7} {'bursts':>7} {'burst ms':>9} "
          f"{'period ms':>10} {'duty':>6}")
    for ex in SETTINGS:
        off = out[ex][0]
        g = lambda k: [off[x][k] for x in off if off[x][k] is not None]
        print(f"    {ex:>6.2f} {np.mean([off[x]['rate_hz'] for x in off]):>7.0f} "
              f"{np.mean([off[x]['n_bursts'] for x in off]):>7.1f} "
              f"{(np.mean(g('burst_ms')) if g('burst_ms') else float('nan')):>9.1f} "
              f"{(np.mean(g('period_ms')) if g('period_ms') else float('nan')):>10.1f} "
              f"{(np.mean(g('duty')) if g('duty') else float('nan')):>6.3f}")
    print("    COLL-8a AB/PD 1 at extra 0.00: 483 sp/s, OFF 0/10 dead, kill 0/10")

    json.dump(summ, open(os.path.join(HERE, 'coll9_summary.json'), 'w'),
              indent=1, default=float)
    print(f"\n  total {(time.process_time() - T0) / 60:.1f} min cpu, "
          f"{(time.time() - W0) / 60:.1f} min wall of {STOP_MIN}")


if __name__ == '__main__':
    main()
