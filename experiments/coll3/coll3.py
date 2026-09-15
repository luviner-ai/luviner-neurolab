"""COLL-3: does learning kill networks on cells that are not the STG burster?

Registered in BRAIN-SIMULATION.md, `## Run: COLL-3 (registration and
economic gate)`, commit 234901f, before any arm ran.

Protocol as `## Run: COLL-2`: two arms, STDP on and off, same wiring and the
same homeostasis in both, on four wirings per cell model. Detector: the
cell-count quorum of RESCORE-1, Q = 24 of 48 in the final 10 s, with the
burst clause dropped and the burst count reported anyway (declared in the
registration).

Two cell models, one constructor call apart:
  cortical   NeuronPopulation.cortical  -- what SiliconDish ships
  HH squid   NeuronPopulation(kinetics='squid') -- textbook Hodgkin-Huxley

Drive is recalibrated per cell by the rule fixed in the registration,
because at the dish's own 0.90 neither cell fires at all.
"""

import contextlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
sys.path.insert(0, os.path.join(HERE, '..', 'dish'))

from cells import CorticalDish, HHDish, N_EXC, N_INH, BUDGET, build
from luviner.biophysics import DishBatch
from luviner.biophysics.population import NeuronPopulation
from stg_pm import burst_windows

SEEDS = [30, 31, 32, 33]
T_MS = 60000.0
CHUNK_MS = 5000.0
Q = 24
GRID = {'cortical': [1.0, 1.1, 1.2, 1.3, 1.5],
        'hh': [3.0, 4.0, 5.0, 6.0, 8.0]}
CELLS = {'cortical': CorticalDish, 'hh': HHDish}
CAL_MS = 2000.0
RATE_BAR = 20.0
STOP_MIN = 40.0
T0 = time.process_time()
W0 = time.time()


def budget_check(tag):
    used = (time.process_time() - T0) / 60.0
    if used > STOP_MIN:
        raise SystemExit(f"STOP RULE: {used:.1f} min of process time at {tag}")
    return used


@contextlib.contextmanager
def squid_population():
    """Make `DishBatch` build a squid population instead of a cortical one.

    `DishBatch._build` hardcodes `NeuronPopulation.cortical`, the same way
    `SiliconDish._build` does. Rather than duplicate its 70-line body, the
    factory is swapped for the duration of the constructor. `gate_batch`
    proves the result matches a single HH dish bit for bit, which is
    `dish-batch`'s own acceptance criterion re-run on the new path.
    """
    orig = NeuronPopulation.cortical
    NeuronPopulation.cortical = classmethod(
        lambda cls, n, **kw: cls(n, kinetics='squid'))
    try:
        yield
    finally:
        NeuronPopulation.cortical = orig


def make_batch(cell, dishes):
    if cell == 'hh':
        with squid_population():
            return DishBatch(dishes)
    return DishBatch(dishes)


def sp(pop, i):
    return np.asarray(pop.spike_times[i], dtype=float)


# ------------------------------------------------------------------ stage A
def calibrate():
    print("=" * 72)
    print("Stage A  drive calibration, registered rule: lowest drive on the")
    print("         grid with 48/48 cells firing and >= 20 Hz/cell, STDP off")
    print("=" * 72)
    out = {}
    for cell, cls in CELLS.items():
        print(f"\n  {cell}")
        chosen, rows = None, []
        for dr in GRID[cell]:
            budget_check(f"calibrate/{cell}/{dr}")
            d = build(cls, SEEDS[0], drive=dr, plastic=False)
            d.run(CAL_MS, learn=False)
            ns = sum(int(sp(d.pop, i).size) for i in range(N_EXC))
            nc = sum(1 for i in range(N_EXC) if sp(d.pop, i).size > 0)
            hz = ns / (N_EXC * CAL_MS / 1000.0)
            rows.append(dict(drive=dr, spikes=ns, cells=nc, hz=hz))
            ok = (nc == N_EXC and hz >= RATE_BAR)
            print(f"    drive {dr:5.2f}   {nc:2d}/48 cells   {hz:6.1f} Hz/cell"
                  f"   {'<- chosen' if ok and chosen is None else ''}")
            if ok and chosen is None:
                chosen = dr
        if chosen is None:
            raise SystemExit(f"NO DRIVE on the {cell} grid meets the rule; "
                             f"the row stops here rather than widening it")
        out[cell] = dict(chosen=chosen, scan=rows)
    print()
    for cell in CELLS:
        print(f"    {cell:>9}: drive {out[cell]['chosen']}")
    json.dump(out, open(os.path.join(HERE, 'calibration.json'), 'w'), indent=1)
    return {c: out[c]['chosen'] for c in out}


# ------------------------------------------------------------------ stage B
def gate_batch(drives):
    """dish-batch's acceptance criterion, re-run on both cell paths."""
    print("\n" + "=" * 72)
    print("Stage B (gate)  batched == separate, on both cell models")
    print("=" * 72 + "\n")
    out = {}
    for cell, cls in CELLS.items():
        budget_check(f"gate/{cell}")
        dr = drives[cell]
        singles = [build(cls, s, drive=dr) for s in SEEDS[:2]]
        refs = [build(cls, s, drive=dr) for s in SEEDS[:2]]
        b = make_batch(cell, singles)
        b.run(2000.0, learn=True)
        for r in refs:
            r.run(2000.0, learn=True)
        same, dw = True, 0.0
        for k, r in enumerate(refs):
            gi = b.cells(k)
            for i in range(N_EXC):
                x, y = sp(b.pop, int(gi[i])), sp(r.pop, i)
                if x.shape != y.shape or not np.array_equal(x, y):
                    same = False
            dw = max(dw, float(np.max(np.abs(b.weights(k) - r.conn.weights))))
        out[cell] = dict(spikes_identical=bool(same), max_dw=dw)
        print(f"    {cell:>9}  2 s, K=2, STDP on:  spikes identical {same}, "
              f"max |dw| {dw:.3e}")
    return out


def lot(cell, learn, drive):
    """One (cell, arm): four wirings batched, T_MS ms, result to disk."""
    tag = f"{cell}_{'on' if learn else 'off'}"
    path = os.path.join(HERE, f'coll3_{tag}.json')
    if os.path.exists(path):
        print(f"    {tag}: already on disk, skipped")
        return json.load(open(path))
    budget_check(f"lot/{tag}")
    t0, w0 = time.process_time(), time.time()
    dishes = [build(CELLS[cell], s, drive=drive, plastic=learn)
              for s in SEEDS]
    b = make_batch(cell, dishes)
    n_chunks = int(round(T_MS / CHUNK_MS))
    trace = [[] for _ in SEEDS]
    for c in range(n_chunks):
        lo = c * CHUNK_MS
        b.run(CHUNK_MS, learn=learn)
        for k in range(len(SEEDS)):
            gi = b.cells(k)
            spk = [sp(b.pop, int(gi[i])) for i in range(N_EXC)]
            act = sum(1 for s in spk if ((s > lo) & (s <= lo + CHUNK_MS)).any())
            cnt = sum(int(((s > lo) & (s <= lo + CHUNK_MS)).sum()) for s in spk)
            V = b.pop.V[gi[:N_EXC]]
            w = b.weights(k)
            tot = float(np.mean(np.bincount(
                dishes[k].post, weights=w, minlength=N_EXC)))
            trace[k].append(dict(t_s=(lo + CHUNK_MS) / 1000.0, cells=act,
                                 spikes=cnt, mean_V=float(np.mean(V)),
                                 tot_w=tot,
                                 frac_min=float(np.mean(w <= 1e-9))))
    rows = []
    for k, seed in enumerate(SEEDS):
        gi = b.cells(k)
        spk = [sp(b.pop, int(gi[i])) for i in range(N_EXC)]
        tail = sum(1 for s in spk if (s > T_MS - 10000.0).any())
        wins, _ = burst_windows(spk, T_MS - 30000.0, T_MS, N_EXC)
        w = b.weights(k)
        rows.append(dict(
            cell=cell, learn=learn, seed=seed, drive=drive, t_ms=T_MS,
            cells_tail=tail, alive=bool(tail >= Q), bursts_tail=len(wins),
            mean_V_tail=trace[k][-1]['mean_V'],
            tot_w_final=trace[k][-1]['tot_w'], budget=BUDGET,
            frac_at_min=trace[k][-1]['frac_min'],
            counts=[int(s.size) for s in spk], trace=trace[k]))
        print(f"    {tag} seed {seed}: {tail:2d}/48 cells, "
              f"{'ALIVE' if tail >= Q else 'DEAD '}, {len(wins):2d} bursts, "
              f"V {trace[k][-1]['mean_V']:7.2f} mV, "
              f"sum_w {trace[k][-1]['tot_w']:.4f} vs budget {BUDGET}")
    res = dict(cell=cell, learn=learn, drive=drive, rows=rows,
               cpu_s=time.process_time() - t0, wall_s=time.time() - w0)
    json.dump(res, open(path, 'w'), indent=1)
    print(f"    {tag}: {res['cpu_s'] / 60:.1f} min cpu, "
          f"{res['wall_s'] / 60:.1f} min wall -> {os.path.basename(path)}")
    return res


def main():
    # Calibration and the batch gate are cached: they are deterministic and
    # already cost 2.9 min of the 45 authorised. Re-running them to relaunch
    # the lots would spend that twice.
    cal_p = os.path.join(HERE, 'calibration.json')
    gate_p = os.path.join(HERE, 'gate.json')
    if os.path.exists(cal_p):
        cal = json.load(open(cal_p))
        drives = {c: cal[c]['chosen'] for c in cal}
        print(f"Stage A  cached: " + ", ".join(
            f"{c} drive {drives[c]}" for c in drives))
    else:
        drives = calibrate()
    if os.path.exists(gate_p):
        print(f"Stage B  cached: {json.load(open(gate_p))}")
    else:
        gate = gate_batch(drives)
        json.dump(gate, open(gate_p, 'w'), indent=1)
    print(f"\n  after calibration and gate: "
          f"{(time.process_time() - T0) / 60:.1f} min of 40\n")
    print("=" * 72)
    print(f"Stage C  four lots, K=4 batched, {T_MS / 1000:.0f} s each")
    print("=" * 72 + "\n")
    for cell in CELLS:
        for learn in (False, True):
            lot(cell, learn, drives[cell])
            print(f"      elapsed {(time.process_time() - T0) / 60:.1f} min "
                  f"of 40\n")
    print(f"  total: {(time.process_time() - T0) / 60:.1f} min process, "
          f"{(time.time() - W0) / 60:.1f} min wall")


if __name__ == '__main__':
    main()
