"""DUR-1: does burst duration control sorting strength, or did it look that way?

`## Run: PM-4` found rho(magnitude, burst duration) = -0.47 at n = 2 --
shorter bursts sorting harder -- and recorded it post-hoc with its own
precondition attached: duration must be fixed as a PREDICTOR before
anyone looks. This is that test.

The lever, and why it is not the one the brief named
----------------------------------------------------
N8b's `U` lever was measured on the CORTICAL dish, where duration is
recruitment-limited. This is the STG dish, where `## Run: STIM-1`
established the rhythm is intrinsically generated and the synapses do not
recruit at rest. `U` is a synaptic knob on a substrate whose duration is
not synaptic; it would likely move nothing, and a null from an inert
lever cannot be told apart from a null that refutes the hypothesis.

The lever used instead is `g_ie`, from this line's own measured band:
0.03 gives 437 ms and 0.10 gives 364 ms **at an identical 0.73 Hz**.
Inhibition truncates a burst whose length is set intrinsically.

Registered before the numbers
-----------------------------
  P1  within each wiring, the causal-count excess is LARGER at
      g_ie = 0.10 (shorter bursts) than at 0.03, in >= 4 of 5 wirings
  P2  the same for strong-minus-weak P/M, in >= 4 of 5
  P3  across the five wirings at g_ie = 0.06, duration correlates
      NEGATIVELY with sorting strength, rho <= -0.4

**The direction is counterintuitive and is stated in advance**: more
inhibition, shorter bursts, STRONGER sorting -- backwards from "more time
to accumulate". The mechanism that would explain it is named now so it
cannot be invented afterwards: PM-4 measured sorting as a COUNTING effect
whose exponential weighting added nothing, so it is carried by the causal
FRACTION rather than the count, and a tighter burst orders its spikes
more sharply.

Covariate plan, named in advance
--------------------------------
Rate is recorded per run. If mean rate differs between levels by more
than 10%, rate becomes a covariate: each level's runs split at the median
rate and the direction must hold inside both strata. Above 25%, the
experimental arm is void and only P3 is scored.

Duration is recorded per run. **If a level does not actually move
duration in a given wiring, that wiring contributes no evidence about
duration** and is reported as such rather than counted.

Falsifiers
----------
  * P1 and P2 both fail -> the anecdote dies as registered, and duration
    does not come back without a different preparation.
  * P1/P2 pass but P3 fails -> the intervention moved sorting while
    natural duration did not track it. That is evidence for `g_ie`, NOT
    for duration, and is reported in those words.
  * the guard fails at any level in any wiring -> that cell is void, not
    reinterpreted.

Nothing is pooled across wirings: `## Run: PM-3` found a regime
bifurcation, and the experimental arm is paired within wiring precisely
so it cannot contaminate the comparison.

Cost: 3 levels x 60 s x 3.7 s/s = 11.1 min per wiring, one wiring per
batch, five batches.
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
from luviner.biophysics.dish import rule_with_integral
from luviner.biophysics.stg_dish import STGDish, STGDriver
from stg_pm import (MAX_LAG, N_EXC, N_INH, CONN, SECONDS, TAU_MINUS,
                    TAU_PLUS, VALIDATION, W_REC, burst_windows, spearman)

LEVELS = (0.03, 0.06, 0.10)
OBS_LEVEL = 0.06


def one(seed, g_ie):
    sets = json.load(open(VALIDATION))['model']['cell_conductances_mS_cm2']
    g_e, g_i = sets['ABPD'], sets['LP']

    def factory(n):
        cells = ([stg.STGNeuron(list(g_e), name=f'E{i}') for i in range(N_EXC)]
                 + [stg.STGNeuron(list(g_i), name=f'I{i}')
                    for i in range(N_INH)])
        return STGDriver(stg.STGNetwork(cells, []))

    d = STGDish(factory, n_exc=N_EXC, n_inh=N_INH, connectivity=CONN,
                w_rec=W_REC, g_ie=g_ie, budget=(N_EXC - 1) * CONN * W_REC,
                depress=(0.20, 400.0), seed=seed)
    d.run(SECONDS * 1000.0, learn=True)

    t_lo, t_hi = SECONDS * 1000.0 * 0.5, SECONDS * 1000.0
    spikes = [d.spikes(i, t_lo=t_lo) for i in range(N_EXC)]
    windows, widths = burst_windows(spikes, t_lo, t_hi, N_EXC)
    w = d.conn.weights
    strong = w > np.median(w)

    P = {True: 0.0, False: 0.0}
    M = {True: 0.0, False: 0.0}
    causal = {True: 0, False: 0}
    anti = {True: 0, False: 0}
    for a, b in windows:
        cache = {}

        def sp(c):
            if c not in cache:
                s = spikes[c]
                cache[c] = s[(s >= a) & (s <= b)]
            return cache[c]
        for j in range(w.size):
            pre_s, post_s = sp(d.pre_e[j]), sp(d.post_e[j])
            if pre_s.size == 0 or post_s.size == 0:
                continue
            lag = (post_s[None, :] - pre_s[:, None]).ravel()
            lag = lag[np.abs(lag) <= MAX_LAG]
            if lag.size == 0:
                continue
            s = bool(strong[j])
            pos, neg = lag[lag >= 0], lag[lag < 0]
            causal[s] += pos.size
            anti[s] += neg.size
            P[s] += float(np.sum(np.exp(-pos / TAU_PLUS)))
            M[s] += float(np.sum(np.exp(neg / TAU_MINUS)))

    pm = {k: (P[k] / M[k] if M[k] > 0 else float('nan')) for k in (True, False)}
    cf = {k: (causal[k] / (causal[k] + anti[k])
              if causal[k] + anti[k] else float('nan'))
          for k in (True, False)}
    n_burst = len(windows)
    return {
        'seed': seed, 'g_ie': g_ie,
        'n_bursts': n_burst,
        'rate_hz': n_burst / (SECONDS / 2.0),
        'duration_ms': float(widths.mean()) if widths.size else 0.0,
        'guard_ok': bool(n_burst >= 3 and (widths.max() <= 0.2 * (t_hi - t_lo)
                                           if widths.size else False)),
        'pm_diff': pm[True] - pm[False],
        'cf_diff': cf[True] - cf[False],
        'capped': int(d.capped.sum()),
        'pre': d.pre_e.tolist(), 'post': d.post_e.tolist(),
        'weights': w.tolist(), 'strong': strong.tolist(),
        'windows': [list(x) for x in windows],
        'spikes': [s.tolist() for s in spikes],
    }


def score():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'dur1_s*.json'))):
        rows += json.load(open(f))
    seeds = sorted({r['seed'] for r in rows})
    if len(seeds) < 5:
        print(f"only {len(seeds)}/5 wirings -- not scoring a partial set")
        return
    by = {(r['seed'], r['g_ie']): r for r in rows}

    print(f"\n{'seed':>5} {'g_ie':>6} {'dur ms':>8} {'rate Hz':>8} "
          f"{'causal diff':>12} {'P/M diff':>10} {'guard':>6}")
    for s in seeds:
        for g in LEVELS:
            r = by[(s, g)]
            print(f"{s:>5} {g:>6.2f} {r['duration_ms']:>8.0f} "
                  f"{r['rate_hz']:>8.2f} {r['cf_diff']:>+12.5f} "
                  f"{r['pm_diff']:>+10.4f} {str(r['guard_ok']):>6}")
        print()

    void = [(s, g) for s in seeds for g in LEVELS if not by[(s, g)]['guard_ok']]
    if void:
        print(f"  VOID cells (guard failed): {void}")
    # did the lever move duration in each wiring?
    inert = [s for s in seeds
             if abs(by[(s, 0.10)]['duration_ms'] - by[(s, 0.03)]['duration_ms'])
             < 0.05 * max(by[(s, 0.03)]['duration_ms'], 1e-9)]
    if inert:
        print(f"  LEVER INERT (<5% duration change): wirings {inert} "
              f"contribute no evidence about duration")

    # A void CELL is (wiring, level), which is what the registration says.
    # An earlier version of this scorer required all three levels to pass
    # and so dropped wirings 30 and 34 from a comparison that never touches
    # their void level -- turning 5/5 into 3/3 and printing that the
    # anecdote had died. P1 and P2 contrast 0.03 against 0.10 only.
    usable = [s for s in seeds if s not in inert
              and by[(s, 0.03)]['guard_ok'] and by[(s, 0.10)]['guard_ok']]
    obs = [s for s in seeds if by[(s, OBS_LEVEL)]['guard_ok']]
    print(f"\n  wirings usable for P1/P2 (0.03 and 0.10 valid): {usable}")
    print(f"  wirings usable for P3 (0.06 valid): {obs}")

    rates = {g: np.array([by[(s, g)]['rate_hz'] for s in usable])
             for g in (0.03, 0.10)}
    spread = ((max(r.mean() for r in rates.values())
               - min(r.mean() for r in rates.values()))
              / max(1e-9, np.mean([r.mean() for r in rates.values()])))
    print(f"  rate spread between levels: {100 * spread:.1f}%  "
          + ("(>25%: experimental arm VOID, only P3 scores)" if spread > 0.25
             else "(>10%: rate enters as a covariate)" if spread > 0.10
             else "(rate matched, no covariate needed)"))

    print(f"\n=== scored ===")
    if spread > 0.25:
        print("  P1, P2 VOID -- rate was not held; only P3 is scored.")
        p1 = p2 = None
    else:
        p1 = sum(1 for s in usable
                 if by[(s, 0.10)]['cf_diff'] > by[(s, 0.03)]['cf_diff'])
        p2 = sum(1 for s in usable
                 if by[(s, 0.10)]['pm_diff'] > by[(s, 0.03)]['pm_diff'])
        print(f"  P1  causal excess larger at g_ie=0.10   {p1}/{len(usable)}   "
              f"{'MET' if p1 >= 4 else 'NOT MET'}")
        print(f"  P2  P/M diff larger at g_ie=0.10        {p2}/{len(usable)}   "
              f"{'MET' if p2 >= 4 else 'NOT MET'}")
    dur = np.array([by[(s, OBS_LEVEL)]['duration_ms'] for s in obs])
    srt = np.array([by[(s, OBS_LEVEL)]['cf_diff'] for s in obs])
    span = (dur.max() - dur.min()) / dur.mean() if dur.size else 0.0
    if dur.size < 4 or span < 0.05:
        print(f"  P3  UNINFORMATIVE, not failed: n={dur.size} wirings, "
              f"duration spans {dur.min():.0f}-{dur.max():.0f} ms "
              f"({100 * span:.1f}%)")
        print(f"      The observational arm needs natural duration variation "
              f"and there is none\n      at this level. rho over these points "
              f"is {spearman(dur, srt):+.3f} and means nothing.")
        rho = None
    else:
        rho = spearman(dur, srt)
        print(f"  P3  rho(duration, sorting) at g_ie=0.06  {rho:+.3f}   "
              f"{'MET' if rho <= -0.4 else 'NOT MET'}")

    if p1 is not None and p1 >= 4 and p2 >= 4 and rho is not None and rho > -0.4:
        print(f"\n  The intervention moved sorting and natural duration did")
        print(f"  not track it. That is evidence for g_ie, NOT for duration.")
    elif p1 is not None and p1 < 4 and p2 < 4:
        print(f"\n  The anecdote dies as registered. Duration does not come")
        print(f"  back without a different preparation.")


def main():
    if '--score' in sys.argv:
        score()
        return
    i = sys.argv.index('--seeds')
    lo, hi = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    seeds = list(range(lo, hi + 1))
    est = len(seeds) * len(LEVELS) * SECONDS * 3.7
    rule = rule_with_integral(-0.134)
    print(f"DUR-1 seeds {seeds} x g_ie {LEVELS}\n"
          f"estimated {est / 60:.1f} min   neutral ratio "
          f"{rule.A_minus / rule.A_plus:.4f}", flush=True)
    if est > 15 * 60:
        print("over the line -- smaller batches")
        return
    t0 = time.time()
    print(f"\n{'seed':>5} {'g_ie':>6} {'dur ms':>8} {'rate Hz':>8} "
          f"{'causal diff':>12} {'P/M diff':>10} {'guard':>6}")
    rows = []
    for s in seeds:
        for g in LEVELS:
            r = one(s, g)
            rows.append(r)
            print(f"{s:>5} {g:>6.2f} {r['duration_ms']:>8.0f} "
                  f"{r['rate_hz']:>8.2f} {r['cf_diff']:>+12.5f} "
                  f"{r['pm_diff']:>+10.4f} {str(r['guard_ok']):>6}", flush=True)
    p = os.path.join(HERE, f'dur1_s{lo}_{hi}.json')
    json.dump(rows, open(p, 'w'))
    print(f"\n{(time.time() - t0) / 60:.1f} min actual\nwritten {p}")


if __name__ == '__main__':
    main()
