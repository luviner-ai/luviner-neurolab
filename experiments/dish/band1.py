"""BAND-1: is the silent absorbing state reachable from a BAND?

DUR-1 and DUR-2 each killed more wirings in the middle of a lever's range
than at either end:

    g_ie  0.030  0/5      CaS  1.0  0/5
    g_ie  0.060  2/5      CaS  0.8  4/5
    g_ie  0.100  0/5      CaS  0.6  1/5

Two unrelated parameters -- one synaptic, one intracellular. Either the
absorbing state is reachable from a band, or 2/5 and 4/5 out of five is
small-n noise that happened twice.

Registered
----------
  P1  on this finer grid at ten wirings per level, death probability is
      NON-MONOTONE with an INTERIOR MAXIMUM, for BOTH levers: some
      interior level kills strictly more than both the lowest and the
      highest level tested.

Falsifier: monotone OR flat death across the finer grid -> the band was
small-n noise from five wirings, BAND-1 is wrong, the row closes. In
those words.

Order and gate, fixed before the numbers
----------------------------------------
`CaS` runs first: it carries the stronger signal, and "it came first
chronologically" is not a reason.

  * CaS monotone or flat -> the falsifier fires and `g_ie` DOES NOT RUN.
    If the strong evidence dissolves, the weak evidence cannot carry the
    row alone.
  * CaS confirms the band -> `g_ie` runs ANYWAY. P1 asks for both levers;
    a band on one lever is a different result from the registered one,
    and confirming CaS does not licence stopping early and claiming P1.

The grid step
-------------
The standing rule says a grid coarser than the transition steps over it,
and here **the band is the transition**. The step is justified by being
ANCHORED on a point already known to be interior -- 0.8 for CaS, 0.060
for g_ie -- so the grid cannot miss the band entirely whatever its width.
A band narrower than the step appears as one killing level with living
neighbours: still non-monotone, still an answer to P1. Resolving the
band's width is not the question.

The detector, validated before use
----------------------------------
Death is fast: the one wiring examined in detail fired for 2,821 ms and
never again in 57 s. So 10 s should classify, at a quarter of 60 s. That
rests on ONE observation, so `--stage0` re-runs the five `g_ie = 0.060`
wirings at 10 s and must reproduce the 60-second classification **5/5**.
A 4/5 is a redesign, not an interpretation.

A late death misread as alive biases every level's rate down uniformly,
which can mask a band but cannot manufacture one.

    dead := spiked during the first half, and produced no spike in the
            final 5 s. Distinguishes death from never having started.
"""

import glob
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
sys.path.insert(0, HERE)

from luviner.biophysics import stg
from luviner.biophysics.stg_dish import STGDish, STGDriver
from stg_pm import N_EXC, N_INH, CONN, W_REC, VALIDATION

# 35 s / 10 s tail: the probe found CaS deaths at 22,135 ms, where the
# g_ie deaths stage 0 validated on were at 2,821 ms. A 10 s window could
# never have seen one. Deaths later than 25 s are still missed, which
# biases every level's rate uniformly downward -- it can mask a band but
# cannot invent one.
T_MS = 35000.0
QUIET_TAIL_MS = 10000.0
CAS_LEVELS = (1.0, 0.9, 0.8, 0.6)
# Four levels, arranged like the CaS arm so the two are comparable
# without clauses: extreme / shoulder / known-interior / far extreme.
GIE_LEVELS = (0.030, 0.045, 0.060, 0.090)
CAS_INDEX = 2
KNOWN_60S = {30: 'dead', 31: 'alive', 32: 'alive', 33: 'alive', 34: 'dead'}


def run_once(seed, g_ie, cas_mult):
    sets = json.load(open(VALIDATION))['model']['cell_conductances_mS_cm2']
    g_e, g_i = list(sets['ABPD']), sets['LP']
    g_e[CAS_INDEX] = sets['ABPD'][CAS_INDEX] * cas_mult

    def factory(n):
        cells = ([stg.STGNeuron(list(g_e), name=f'E{i}') for i in range(N_EXC)]
                 + [stg.STGNeuron(list(g_i), name=f'I{i}')
                    for i in range(N_INH)])
        return STGDriver(stg.STGNetwork(cells, []))

    d = STGDish(factory, n_exc=N_EXC, n_inh=N_INH, connectivity=CONN,
                w_rec=W_REC, g_ie=g_ie, budget=(N_EXC - 1) * CONN * W_REC,
                depress=(0.20, 400.0), seed=seed)
    d.run(T_MS, learn=True)
    sp = [d.spikes(i) for i in range(N_EXC)]
    n_early = sum(int(np.sum(s < T_MS / 2)) for s in sp)
    n_tail = sum(int(np.sum(s >= T_MS - QUIET_TAIL_MS)) for s in sp)
    return {'seed': seed, 'g_ie': g_ie, 'cas_mult': cas_mult,
            'n_early': n_early, 'n_tail': n_tail,
            'started': n_early > 0,
            'dead': bool(n_early > 0 and n_tail == 0)}


def stage0():
    print("stage 0: does a 10 s run reproduce the 60 s classification?\n")
    print(f"{'seed':>5} {'early':>7} {'tail':>6} {'10 s says':>10} "
          f"{'60 s said':>10} {'':>6}")
    t0 = time.time()
    ok = 0
    for s in sorted(KNOWN_60S):
        r = run_once(s, 0.060, 1.0)
        says = 'dead' if r['dead'] else ('alive' if r['started']
                                         else 'never started')
        match = says == KNOWN_60S[s]
        ok += match
        print(f"{s:>5} {r['n_early']:>7} {r['n_tail']:>6} {says:>10} "
              f"{KNOWN_60S[s]:>10} {'ok' if match else 'MISMATCH':>6}",
              flush=True)
    print(f"\n  {ok}/5 agree.  {'DETECTOR VALIDATED -- the map may run.' if ok == 5 else 'NOT 5/5 -- the map is REDESIGNED, not reinterpreted.'}")
    print(f"  {(time.time() - t0) / 60:.1f} min")
    return ok == 5


def main():
    if '--stage0' in sys.argv:
        stage0()
        return
    arm = sys.argv[sys.argv.index('--arm') + 1]
    i = sys.argv.index('--seeds')
    lo, hi = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    seeds = list(range(lo, hi + 1))
    levels = CAS_LEVELS if arm == 'cas' else GIE_LEVELS
    est = len(seeds) * len(levels) * (T_MS / 1000.0) * 3.7
    print(f"BAND-1 {arm} arm, seeds {seeds} x {levels}\n"
          f"estimated {est / 60:.1f} min", flush=True)
    t0 = time.time()
    print(f"\n{'seed':>5} {'level':>7} {'early':>7} {'tail':>6} {'verdict':>14}")
    rows = []
    for s in seeds:
        for L in levels:
            r = (run_once(s, 0.03, L) if arm == 'cas'
                 else run_once(s, L, 1.0))
            rows.append(r)
            v = ('DEAD' if r['dead'] else
                 ('alive' if r['started'] else 'never started'))
            print(f"{s:>5} {L:>7.3f} {r['n_early']:>7} {r['n_tail']:>6} "
                  f"{v:>14}", flush=True)
    p = os.path.join(HERE, f'band1_{arm}_s{lo}_{hi}.json')
    json.dump(rows, open(p, 'w'))
    print(f"\n{(time.time() - t0) / 60:.1f} min actual\nwritten {p}")


if __name__ == '__main__':
    main()
