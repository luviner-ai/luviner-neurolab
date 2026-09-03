"""Where a plasticity rule is neutral, and why it is not at zero area.

`STDP.window_integral` is the net change per pair *when the lags are
spread evenly over the window*. A stimulated network does not spread them
evenly: it concentrates them near zero, where the window is
potentiation-dominated because its causal side is taller. So a rule with
zero area still potentiates, and `rule_with_integral(0)` -- the obvious
ablation, and the one this document registered as the DishBrain
falsifier -- is not an ablation at all.

Neutrality is scale-free and is a *ratio*, not an area:

    A_minus / A_plus = P / M

with P and M the two kernels of the lag distribution the network itself
produces. The ladder below is therefore run in the ratio, not in the
integral: the neutral integral moves with A_plus, the neutral ratio does
not.

Registered before the measurement: the realised mean weight change
crosses zero at the ratio P/M computed from a frozen dish's lags, and
P/M is below tau_plus/tau_minus = 0.499^-1 ... below the value that would
make the window integral zero, so the neutral rule is net-depressing by
area.

A first version of this ladder ran 24 stimuli and measured the wrong
thing: the weights lost 63% of their initial value and 3-4% sat on a
rail, so the ladder was reading the clipping, not the drift. Here the
protocol is short and A_plus small enough that nothing reaches a bound,
which is checked and reported rather than assumed.

Falsified if the crossing is at the ratio that zeroes the area
(tau_plus/tau_minus = 0.499), or if it does not track the prediction.

Usage:  python experiments/dish/neutral_point.py
"""

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from luviner.biophysics import (SiliconDish, STDP, rule_with_expected_change,
                                lag_kernels, observed_lags)

STIM_SITES = (0, 1, 2)
N_PULSES, ISI = 9, 300.0
A_PLUS = 0.003
TAU_PLUS, TAU_MINUS = 16.8, 33.7
RATIOS = (0.30, 0.42, 0.50, 0.55, 0.62, 0.75, 0.95)
SEEDS = (0, 1, 2, 3, 4)


def protocol(dish, learn=True):
    t0 = dish.t
    for k in range(N_PULSES):
        dish.pulse(STIM_SITES[k % len(STIM_SITES)], follow=ISI - 1.0,
                   learn=learn)
    return t0, dish.t


def predicted(seed):
    """P/M from a frozen dish, whose lags cannot depend on the rule."""
    d = SiliconDish(seed=seed, budget=None, plastic=False)
    d.run(100.0, learn=False)
    t0, t1 = protocol(d, learn=False)
    lags = observed_lags(d, t0, t1)
    P, M = lag_kernels(lags, TAU_PLUS, TAU_MINUS)
    return P / M, lags.size, float(np.mean(lags > 0))


def measured(seed, ratio):
    rule = STDP(A_plus=A_PLUS, A_minus=A_PLUS * ratio, tau_plus=TAU_PLUS,
                tau_minus=TAU_MINUS, w_min=0.0, w_max=0.06, mu=0.0)
    d = SiliconDish(seed=seed, budget=None, rule=rule)
    d.run(100.0)
    before = d.weights.copy()
    protocol(d)
    dw = d.weights - before
    rail = float(np.mean((d.weights <= 1e-9)
                         | (d.weights >= rule.w_max - 1e-9)))
    return float(dw.mean()), rail, float(np.abs(dw).mean() / before.mean())


def main():
    out = {'ratios': list(RATIOS), 'seeds': list(SEEDS), 'A_plus': A_PLUS}
    preds = []
    for s in SEEDS:
        pm, n, causal = predicted(s)
        preds.append(pm)
        print(f"seed {s}: {n} pairs, {100 * causal:.1f}% causal -> "
              f"predicted neutral ratio P/M = {pm:.4f}", flush=True)
    out['predicted_ratio'] = preds
    print(f"prediction: P/M = {np.mean(preds):.4f} "
          f"+/- {np.std(preds, ddof=1):.4f}   "
          f"(zero-area ratio is tau_plus/tau_minus = "
          f"{TAU_PLUS / TAU_MINUS:.4f})\n", flush=True)

    table = np.zeros((len(RATIOS), len(SEEDS)))
    rails = np.zeros_like(table)
    rel = np.zeros_like(table)
    for i, r in enumerate(RATIOS):
        for j, s in enumerate(SEEDS):
            table[i, j], rails[i, j], rel[i, j] = measured(s, r)
        print(f"ratio {r:.2f} (integral "
              f"{A_PLUS * (TAU_PLUS - TAU_MINUS * r):+.4f}): "
              f"mean dw {table[i].mean():+.7f} +/- "
              f"{table[i].std(ddof=1):.7f}   rails "
              f"{100 * rails[i].mean():.1f}%   |dw|/w0 "
              f"{100 * rel[i].mean():.1f}%", flush=True)
    out['measured'] = table.tolist()
    out['at_rails'] = rails.tolist()
    out['relative_move'] = rel.tolist()

    m = table.mean(axis=1)
    x = np.array(RATIOS)
    a, b = np.polyfit(x, m, 1)
    cross = -b / a
    print(f"\nmeasured drift = {a:.6f}*ratio {b:+.6f} -> zero at "
          f"ratio {cross:.4f}")
    print(f"predicted from the lags alone : {np.mean(preds):.4f}")
    print(f"zero-area ratio               : {TAU_PLUS / TAU_MINUS:.4f}")
    out['crossing_measured'] = float(cross)
    out['crossing_predicted'] = float(np.mean(preds))
    out['zero_area_ratio'] = TAU_PLUS / TAU_MINUS
    path = os.path.join(os.path.dirname(__file__), 'neutral_point.json')
    json.dump(out, open(path, 'w'))
    print("written", path)


if __name__ == '__main__':
    t = time.time()
    main()
    print(f"{time.time() - t:.0f}s")
