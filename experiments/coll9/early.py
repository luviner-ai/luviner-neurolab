"""COLL-9 control: the early rate, which cannot depend on survival.

Registered in BRAIN-SIMULATION.md, `## Run: COLL-9 (registration addendum --
the early-rate control)`, before any number existed. Same seeds, same
deterministic integration, 5 s: this is the opening of the run already on
disk, taken before any network has died.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from coll9 import CELL, G_IE, SEEDS, SETTINGS, make, slug   # noqa: E402
from stg_pm import N_EXC                                     # noqa: E402

MS = 5000.0
T0 = time.process_time()
out = {}
for ex in SETTINGS:
    b, e = make(SEEDS, ex)
    b.run(MS, learn=False, extra=e)
    rates = []
    for k in range(len(SEEDS)):
        n = sum(int(np.asarray(b.spikes_of(k, i)).size) for i in range(N_EXC))
        rates.append(n / (MS / 1000.0))
    out[ex] = dict(zip(SEEDS, rates))
    print(f"  extra {ex:>6.2f}: early rates " +
          ' '.join(f'{r:.0f}' for r in rates), flush=True)
json.dump({str(k): v for k, v in out.items()},
          open(os.path.join(HERE, 'coll9_early.json'), 'w'), indent=1)

print("\n  P20  does the EARLY rate predict death?")
for ex in SETTINGS:
    off = {r['seed']: r for r in
           json.load(open(os.path.join(HERE, f'coll9_{slug(ex)}_off.json')))['rows']}
    dead = [out[ex][s] for s in SEEDS if not off[s]['alive']]
    live = [out[ex][s] for s in SEEDS if off[s]['alive']]
    if dead and live:
        md, ml = float(np.mean(dead)), float(np.mean(live))
        print(f"    extra {ex:>6.2f}: died n={len(dead)} early {md:.0f} sp/s, "
              f"lived n={len(live)} early {ml:.0f} sp/s, ratio {ml / md:.2f}x")
        print(f"      P20 (>= 2x separation): "
              f"{'HOLDS' if ml / md >= 2.0 or md / ml >= 2.0 else 'REFUTED'}")
    else:
        print(f"    extra {ex:>6.2f}: all {'died' if not live else 'lived'}, "
              f"early rate {np.mean(dead or live):.0f} sp/s -- no contrast")
print(f"\n  control cost {(time.process_time() - T0) / 60:.1f} min")
