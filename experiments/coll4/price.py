"""COLL-4 economic gate and P1b, before any arm.

The row asks for "an endogenous burster that is not the STG cell". The
first thing to establish is what COLL-2's cell actually was, because the
answer changes the design: `stg_pm.VALIDATION` points at
`experiments/pyloric/validation.json` and `susc1.build` reads
`cell_conductances_mS_cm2['ABPD']` for the 48 excitatory cells --
(100, 2.5, 6, 50, 5, 100, 0.01, 0), which is **`PRINZ_TABLE2['AB/PD 2']`
exactly**. COLL-2's dish IS the pyloric pacemaker cell.

So `pyloric.py` is not a different instrument; it is the same one. What IS
different, and available for free, are the OTHER pacemaker variants --
Prinz 2004 Table 2 lists five AB/PD models, all with published intrinsic
burst periods -- plus the follower cells LP and PY.

P1b is therefore run as a SCAN: which cell models, at the dish's own
sub-rheobase drive with STDP off, produce a collective burst rhythm?
Only those can carry the two arms.
"""

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
sys.path.insert(0, os.path.join(HERE, '..', 'dish'))

from luviner.biophysics import stg
from luviner.biophysics.stg_dish import STGDish, STGDriver
from stg_pm import N_EXC, N_INH, CONN, W_REC, burst_windows

BUDGET = (N_EXC - 1) * CONN * W_REC
G_IE = 0.030
PROBE_MS = 10000.0
CANDIDATES = ['AB/PD 2', 'AB/PD 1', 'AB/PD 3', 'AB/PD 4', 'AB/PD 5',
              'LP 3', 'PY 1']


def build(cell_name, seed, g_ie=G_IE, plastic=True):
    """The COLL-2 dish with the excitatory conductance set as a parameter.

    Inhibitory cells stay 'LP 3', exactly as `susc1.build` has them, so the
    only thing that changes between arms of this scan is the E cell.
    """
    g_e = list(stg.PRINZ_TABLE2[cell_name])
    g_i = list(stg.PRINZ_TABLE2['LP 3'])

    def factory(n):
        cells = ([stg.STGNeuron(list(g_e), name=f'E{i}') for i in range(N_EXC)]
                 + [stg.STGNeuron(list(g_i), name=f'I{i}')
                    for i in range(N_INH)])
        return STGDriver(stg.STGNetwork(cells, []))

    return STGDish(factory, n_exc=N_EXC, n_inh=N_INH, connectivity=CONN,
                   w_rec=W_REC, g_ie=g_ie, budget=BUDGET,
                   depress=(0.20, 400.0), seed=seed)


def main():
    t0 = time.process_time()
    print("COLL-4 gate: P1b as a scan, plus the price\n")
    print(f"  COLL-2's excitatory cell, read from the file it used:")
    v = json.load(open(os.path.join(
        HERE, '..', 'pyloric', 'validation.json')))
    abpd = tuple(v['model']['cell_conductances_mS_cm2']['ABPD'])
    for name, g in stg.PRINZ_TABLE2.items():
        if tuple(g) == abpd:
            print(f"    validation.json ABPD == PRINZ_TABLE2['{name}']")
    print()
    print(f"  {PROBE_MS / 1000:.0f} s, STDP off, drive as the dish has it, "
          f"g_ie = {G_IE}\n")
    print(f"    {'cell':>9} {'cells firing':>13} {'spikes':>9} "
          f"{'net bursts':>11} {'per 10 s':>9} {'P1b':>6} {'cpu':>7}")
    out = {}
    for name in CANDIDATES:
        t = time.process_time()
        d = build(name, 30)
        d.run(PROBE_MS, learn=False)   # COLL-2's own off arm: learn=False
        sp = [np.asarray(d.spikes(i), dtype=float) for i in range(N_EXC)]
        nc = sum(1 for s in sp if s.size > 0)
        ns = sum(int(s.size) for s in sp)
        wins, _ = burst_windows(sp, 0.0, PROBE_MS, N_EXC)
        rate = len(wins) / (PROBE_MS / 10000.0)
        cpu = time.process_time() - t
        ok = (nc == N_EXC and rate >= 1.0)
        out[name] = dict(cells=nc, spikes=ns, bursts=len(wins),
                         per_10s=rate, p1b=bool(ok), cpu_s=cpu,
                         s_per_s=cpu / (PROBE_MS / 1000.0))
        print(f"    {name:>9} {nc:>9}/48 {ns:>9,d} {len(wins):>11d} "
              f"{rate:>9.1f} {'MET' if ok else 'no':>6} "
              f"{cpu / (PROBE_MS / 1000.0):>6.2f}s/s")
    cost = float(np.median([out[n]['s_per_s'] for n in out]))
    print(f"\n  median cost {cost:.2f} s cpu / s sim, one dish, unbatched")
    print(f"  one 60 s run: {cost * 60 / 60:.1f} min")
    print(f"  2 arms x 4 wirings x 60 s, one cell model: "
          f"{8 * cost * 60 / 60:.1f} min")
    print(f"  gate spent: {(time.process_time() - t0) / 60:.1f} min")
    out['_price'] = dict(median_s_per_s=cost,
                         gate_cpu_min=(time.process_time() - t0) / 60.0)
    json.dump(out, open(os.path.join(HERE, 'gate.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
