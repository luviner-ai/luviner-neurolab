"""PM-3: the replication the post-hoc split demands, with regime controlled.

`## Run: stg-P/M` scored (a) refuted and (b) one clause of two: the
strong-minus-weak `P/M` difference was positive in 10/10 wirings but its
between-wiring spread slightly exceeded its mean, so the registered bar
was not cleared. The ten wirings turned out to be two regimes separable
by burst count, each of which passes the spread clause alone -- **but
that split was noticed after seeing the result**, so it explained
nothing and rescued nothing.

This run makes the regime an independent variable **assigned before the
outcome is measured**, which is the only thing that can turn that
observation into evidence.

Step 1, here: the SCREEN
------------------------
Regime is an emergent property of the wiring, not a knob. So it is
classified from a short run that measures burst rate and nothing else,
and the assignment is written to disk **before** any `P/M` is computed
for that wiring.

The screen has to be validated before it can be trusted, and there is
free ground truth: seeds 0-9 were measured at 60 s in `stg-P/M`, where
seeds **0 and 8** gave 17 bursts and the other eight gave 22.

**Registered before running the screen: a 15 s screen reproduces that
assignment for all ten seeds.** If it does not, the screen is not a
screen and the regime cannot be pre-assigned -- in which case PM-3 as
designed is not runnable and that is reported instead of worked around.

Cost: 3.7 s wall per biological second, so 10 wirings x 15 s = 9.2 min,
one batch under the line.

Step 2, once the screen is validated: sizing
--------------------------------------------
The 17-burst regime was 2 of 10. Reaching ten wirings in the rarer
regime needs roughly fifty screened. **That estimate gets stated and
approved before it is spent, not after** -- the screen's measured
frequency replaces the guess.

The registered predictions for step 3 (from the board, unchanged):

  (1) within regime, strong-minus-weak positive >= 8/10 AND exceeding
      the within-regime spread
  (2) the causal-pair excess positive >= 9/10
  (3) the `P/M` difference tracks the causal-fraction difference across
      wirings

**All three, or the accumulation reading dies. No partial credit** --
the two-clause lesson from stg-P/M, applied deliberately.

The wiring is persisted explicitly this time. In stg-P/M it was
recoverable from the seed and I verified that, but only because `w_in`
happened to be a witness. Persist the thing, not a way to rebuild it.

Usage:  /tmp/venv/bin/python experiments/dish/stg_pm3.py --screen 0 9
"""

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
sys.path.insert(0, HERE)

from luviner.biophysics import stg
from luviner.biophysics.stg_dish import STGDish, STGDriver
from stg_pm import burst_windows, N_EXC, N_INH, CONN, W_REC, G_IE, VALIDATION

SCREEN_SECONDS = 15.0
# ground truth from stg-P/M at 60 s: these seeds gave 17 bursts, the rest 22
KNOWN_17 = {0, 8}


def screen_one(seed, seconds=SCREEN_SECONDS):
    sets = json.load(open(VALIDATION))['model']['cell_conductances_mS_cm2']
    g_e, g_i = sets['ABPD'], sets['LP']
    budget = (N_EXC - 1) * CONN * W_REC

    def factory(n):
        cells = ([stg.STGNeuron(list(g_e), name=f'E{i}') for i in range(N_EXC)]
                 + [stg.STGNeuron(list(g_i), name=f'I{i}')
                    for i in range(N_INH)])
        return STGDriver(stg.STGNetwork(cells, []))

    d = STGDish(factory, n_exc=N_EXC, n_inh=N_INH, connectivity=CONN,
                w_rec=W_REC, g_ie=G_IE, budget=budget,
                depress=(0.20, 400.0), seed=seed)
    d.run(seconds * 1000.0, learn=True)
    t_lo, t_hi = seconds * 1000.0 * 0.5, seconds * 1000.0
    spikes = [d.spikes(i, t_lo=t_lo) for i in range(N_EXC)]
    wins, widths = burst_windows(spikes, t_lo, t_hi, N_EXC)
    half = (t_hi - t_lo) / 1000.0
    return {'seed': seed, 'n_bursts': len(wins),
            'burst_hz': len(wins) / half,
            'dur_ms': float(widths.mean()) if widths.size else 0.0,
            # persist the wiring EXPLICITLY -- not a way to rebuild it
            'pre': d.pre_e.tolist(), 'post': d.post_e.tolist(),
            'capped': int(d.capped.sum())}


def main():
    i = sys.argv.index('--screen')
    lo, hi = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    seeds = list(range(lo, hi + 1))
    est = len(seeds) * SCREEN_SECONDS * 3.7
    print(f"screen, seeds {lo}-{hi}: {len(seeds)} wirings x "
          f"{SCREEN_SECONDS:.0f} s\nestimated {est / 60:.1f} min", flush=True)
    if est > 15 * 60:
        print("over the line for one launch -- smaller batches")
        return
    t0 = time.time()
    print(f"\n{'seed':>4} {'bursts':>7} {'burst Hz':>9} {'dur ms':>7}  regime")
    rows = []
    for s in seeds:
        r = screen_one(s)
        # the classification rule, fixed before the numbers: the two
        # regimes differ ~1.4x in rate, so the midpoint of the 60 s
        # rates (0.57 and 0.73 Hz) is the boundary
        r['regime'] = 'low' if r['burst_hz'] < 0.65 else 'high'
        rows.append(r)
        print(f"{s:>4} {r['n_bursts']:>7} {r['burst_hz']:>9.2f} "
              f"{r['dur_ms']:>7.0f}  {r['regime']}", flush=True)

    known = [r for r in rows if r['seed'] <= 9]
    if known:
        agree = all((r['regime'] == 'low') == (r['seed'] in KNOWN_17)
                    for r in known)
        got_low = sorted(r['seed'] for r in known if r['regime'] == 'low')
        print(f"\nvalidation against the 60 s ground truth:")
        print(f"  screen says low-rate: {got_low}")
        print(f"  60 s said 17-burst:   {sorted(KNOWN_17)}")
        print(f"  -> {'SCREEN VALIDATED' if agree else 'SCREEN FAILED'}")
        if not agree:
            print("  the regime cannot be pre-assigned from a 15 s run, so "
                  "PM-3 as designed is not runnable at this screen length")
    n_low = sum(r['regime'] == 'low' for r in rows)
    print(f"\nlow-rate frequency: {n_low}/{len(rows)}")
    print(f"{(time.time() - t0) / 60:.1f} min actual (estimated {est / 60:.1f})")
    p = os.path.join(HERE, f'stg_pm3_screen_{lo}_{hi}.json')
    json.dump(rows, open(p, 'w'))
    print("written", p)


if __name__ == '__main__':
    main()
