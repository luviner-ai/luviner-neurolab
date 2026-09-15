"""COLL-6 stage 0: the same acceptance criterion, at PY 1 and the arms' own level.

COLL-5 gated on AB/PD 2. PY 1 is a different conductance set, so the gate is
re-run here rather than inherited -- COLL-4's lesson was that a gate at one
setting says nothing about another.

Batched and separate must agree BIT FOR BIT -- a bug reads as gross
disagreement, not a plausible number -- on spikes and on weights, with the
rule on and with it off. Nothing in COLL-5 runs until this passes.
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
sys.path.insert(0, os.path.join(HERE, '..', 'dish'))
sys.path.insert(0, os.path.join(HERE, '..', 'coll4'))
sys.path.insert(0, os.path.join(HERE, '..', 'coll5'))   # BatchedSTGDish, built in COLL-5

from price import build            # COLL-4's builder, unchanged
from stgbatch import BatchedSTGDish
from luviner.biophysics import stg
from stg_pm import N_EXC

CELL = 'PY 1'
BUDGET = (48 - 1) * 0.15 * 0.015
MS = 2000.0


def run_gate(g_ie, learn, seeds=(30, 31, 32)):
    singles = [build(CELL, s, g_ie=g_ie) for s in seeds]
    refs = [build(CELL, s, g_ie=g_ie) for s in seeds]
    b = BatchedSTGDish(singles, list(stg.PRINZ_TABLE2[CELL]),
                       list(stg.PRINZ_TABLE2['LP 3']), g_ie=g_ie,
                       budget=BUDGET)
    t = time.process_time()
    b.run(MS, learn=learn)
    cpu_b = time.process_time() - t
    t = time.process_time()
    for r in refs:
        r.run(MS, learn=learn)
    cpu_s = time.process_time() - t

    same, dw = True, 0.0
    for k, r in enumerate(refs):
        for i in range(N_EXC):
            x = np.asarray(b.spikes_of(k, i), dtype=float)
            y = np.asarray(r.spikes(i), dtype=float)
            if x.shape != y.shape or not np.array_equal(x, y):
                same = False
        dw = max(dw, float(np.max(np.abs(b.weights_of(k) - r.weights))))
    return same, dw, cpu_b, cpu_s


if __name__ == '__main__':
    print("COLL-6 stage 0: batched == separate, bit for bit\n")
    print(f"    {'g_ie':>6} {'arm':>5} {'spikes identical':>17} "
          f"{'max |dw|':>10} {'cpu batched':>12} {'cpu separate':>13} "
          f"{'speedup':>8}")
    ok = True
    for g in (0.045,):        # the arms' own level, and only it
        for learn in (False, True):
            s, dw, cb, cs = run_gate(g, learn)
            ok = ok and s and dw == 0.0
            print(f"    {g:>6.3f} {'on' if learn else 'off':>5} "
                  f"{str(s):>17} {dw:>10.3e} {cb:>11.1f}s {cs:>12.1f}s "
                  f"{cs / cb:>7.2f}x")
    print(f"\n  gate: {'PASS' if ok else 'FAIL -- nothing else runs'}")
