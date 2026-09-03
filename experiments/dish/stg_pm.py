"""Do long, structured bursts read the weight order? N8's instrument, better bursts.

`## Run: N8` measured `P/M` by weight class on 18-46 ms events and got
`+0.0041 +/- 0.0432` -- a tenth of the between-seed spread, with both
classes on the same side of the neutral point in 3 of 5 seeds. Bursts in
that preparation were "approximately a no-op": they left the weight order
where they found it.

The open question was whether that is a property of *bursts* or of
*those* bursts. `## Run: dish-slow2c` built a preparation whose bursts
carry **22.5-24.8 spikes each** over 444-492 ms, in 10 of 10 wirings,
without latching. **If leaders exist anywhere in this model class, they
exist here**, and a null has no "events too brief" escape left.

Two registered claims, hierarchically ordered
----------------------------------------------
`P/M` is an outcome. The mechanism it depends on is that firing ORDER
within a burst tracks a cell's incoming weight: a cell driven harder
should fire earlier. If order does not track weight, `P/M` cannot sort
and we learn *why* rather than only *that*.

**(a) MECHANISM, registered prediction: Spearman rho < 0 between a
cell's median within-burst first-spike latency and its total incoming
weight, in >= 8/10 wirings.** Negative because more weight should mean
earlier firing, hence shorter latency.

**(b) OUTCOME, registered prediction: strong-minus-weak `P/M` positive
in >= 8/10 wirings, with the effect exceeding the between-wiring
spread.** Prior: **45%**. The bursts now have the structure N8's lacked,
but synchronised bursting may impose an order set by intrinsic phase
rather than by synaptic weight -- which is what (a) exists to tell us.

**The binding, fixed before any number exists.** (b) alone scores the
consolidation claim. **(a) can never rescue a null in (b).** Each is
scored separately and reported whatever it says. This is a double-
barrelled registration precisely so that the second barrel cannot be
substituted for the first after the fact.

**The decision tree, also fixed now:**

    (a) fails, (b) fails   the mechanism is absent -- consolidation
                           fails, with its reason attached
    (a) holds, (b) fails   order exists but STDP does not convert it,
                           which points at within-burst timing statistics
    (b) holds              the positive

The instrument, unchanged from N8
----------------------------------
Classes split at the median conductance per wiring, equal by
construction. Burst windows by N8's own definition -- 10 ms detection
bins, 20% active fraction, 50 ms gap merge, then full width at 20% of
each burst's own peak at 1 ms resolution. Lags inside those windows
only, `|lag| <= 250 ms`. The dividing line is the rule's own
`A_minus/A_plus = 0.8299` -- N7's correction -- not zero.

Everything is persisted: spikes, weights, classes, burst windows and
latencies in ONE file, so a third question about these bursts costs no
re-run. `## Run: N8` stored its analyses and not its data and paid for
it; `stg_pilot.py` repeated that miss and this does not.

Cost: 3.7 s wall per biological second measured with plasticity on, so
10 wirings x 60 s is ~37 min, run as three batches under the line.

Usage:  /tmp/venv/bin/python experiments/dish/stg_pm.py --seeds 0 3
"""

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))

from luviner.biophysics import stg
from luviner.biophysics.dish import rule_with_integral
from luviner.biophysics.stg_dish import STGDish, STGDriver

VALIDATION = os.path.join(ROOT, 'experiments', 'pyloric', 'validation.json')
N_EXC, N_INH, CONN = 48, 12, 0.15
W_REC, G_IE, SECONDS = 0.015, 0.06, 60.0
BIN_MS, ACTIVE_FRAC, GAP_MS = 10.0, 0.20, 50.0
FINE_MS, PEAK_FRAC, MAX_LAG = 1.0, 0.20, 250.0
TAU_PLUS, TAU_MINUS = 16.8, 33.7


def burst_windows(spikes_by_cell, t_lo, t_hi, n_exc):
    """N8's definition, mirrored: coarse detection, fine width at 20%."""
    sp = np.concatenate([s for s in spikes_by_cell] + [np.zeros(0)])
    sp = sp[(sp > t_lo) & (sp <= t_hi)]
    coarse = np.arange(t_lo, t_hi + BIN_MS, BIN_MS)
    counts = np.histogram(sp, bins=coarse)[0]
    active = counts >= ACTIVE_FRAC * n_exc
    gap_bins = int(round(GAP_MS / BIN_MS))
    runs, run_gap, inside = [], gap_bins, False
    for i, a in enumerate(active):
        if a:
            if not inside and run_gap >= gap_bins:
                runs.append([i, i])
            elif runs:
                runs[-1][1] = i
            inside, run_gap = True, 0
        else:
            run_gap += 1
            inside = False
    fine = np.arange(t_lo, t_hi + FINE_MS, FINE_MS)
    fcounts = np.histogram(sp, bins=fine)[0]
    windows, widths = [], []
    for a, b in runs:
        lo = max(0, int((coarse[a] - t_lo - 50.0) / FINE_MS))
        hi = min(fcounts.size, int((coarse[b + 1] - t_lo + 50.0) / FINE_MS))
        seg = fcounts[lo:hi]
        if seg.size == 0 or seg.max() == 0:
            continue
        over = np.flatnonzero(seg >= PEAK_FRAC * seg.max())
        widths.append((over[-1] - over[0] + 1) * FINE_MS)
        windows.append((t_lo + (lo + over[0]) * FINE_MS,
                        t_lo + (lo + over[-1] + 1) * FINE_MS))
    return windows, np.array(widths)


def spearman(x, y):
    """Rank correlation, no scipy in this project."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return float('nan')
    rx = np.argsort(np.argsort(x[ok])).astype(float)
    ry = np.argsort(np.argsort(y[ok])).astype(float)
    if rx.std() == 0 or ry.std() == 0:
        return float('nan')
    return float(np.corrcoef(rx, ry)[0, 1])


def one(seed):
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
    d.run(SECONDS * 1000.0, learn=True)

    t_lo, t_hi = SECONDS * 1000.0 * 0.5, SECONDS * 1000.0
    spikes = [d.spikes(i, t_lo=t_lo) for i in range(N_EXC)]
    windows, widths = burst_windows(spikes, t_lo, t_hi, N_EXC)
    w = d.conn.weights
    strong = w > np.median(w)

    # ---- (a) does firing order track incoming weight? -------------------
    w_in = np.zeros(N_EXC)
    np.add.at(w_in, d.post_e, w)
    lat = [[] for _ in range(N_EXC)]
    for a, b in windows:
        for c in range(N_EXC):
            s = spikes[c][(spikes[c] >= a) & (spikes[c] <= b)]
            if s.size:
                lat[c].append(float(s[0] - a))
    med_lat = np.array([np.median(v) if v else np.nan for v in lat])
    rho = spearman(med_lat, w_in)

    # ---- (b) P/M by class, N8's instrument ------------------------------
    P = {True: 0.0, False: 0.0}
    M = {True: 0.0, False: 0.0}
    n_lag = {True: 0, False: 0}
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
            n_lag[s] += lag.size
            P[s] += float(np.sum(np.exp(-lag[lag >= 0] / TAU_PLUS)))
            M[s] += float(np.sum(np.exp(lag[lag < 0] / TAU_MINUS)))
    pm = {k: (P[k] / M[k] if M[k] > 0 else float('nan'))
          for k in (True, False)}

    return {
        'seed': seed, 'n_bursts': len(windows),
        'dur_ms': float(widths.mean()) if widths.size else 0.0,
        'rho_latency_weight': rho,
        'pm_strong': pm[True], 'pm_weak': pm[False],
        'pm_diff': pm[True] - pm[False],
        'n_lags_strong': n_lag[True], 'n_lags_weak': n_lag[False],
        'capped': int(d.capped.sum()),
        # everything, in one file
        'spikes': [s.tolist() for s in spikes],
        'weights': w.tolist(), 'strong': strong.tolist(),
        'windows': [list(x) for x in windows],
        'w_in': w_in.tolist(), 'median_latency': med_lat.tolist(),
    }


def main():
    i = sys.argv.index('--seeds')
    lo, hi = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    seeds = list(range(lo, hi + 1))
    rule = rule_with_integral(-0.134)
    ratio = rule.A_minus / rule.A_plus
    est = len(seeds) * SECONDS * 3.7
    print(f"P/M by class, seeds {lo}-{hi}: {len(seeds)} wirings x "
          f"{SECONDS:.0f} s\nestimated {est / 60:.1f} min   "
          f"neutral ratio A-/A+ = {ratio:.4f}", flush=True)
    if est > 15 * 60:
        print("over the line for one launch -- use smaller batches")
        return
    t0 = time.time()
    print(f"\n{'seed':>4} {'bursts':>7} {'dur ms':>7} {'rho':>7} "
          f"{'P/M strong':>11} {'P/M weak':>9} {'diff':>9} {'lags':>10}")
    rows = []
    for s in seeds:
        r = one(s)
        rows.append(r)
        print(f"{s:>4} {r['n_bursts']:>7} {r['dur_ms']:>7.0f} "
              f"{r['rho_latency_weight']:>+7.3f} {r['pm_strong']:>11.4f} "
              f"{r['pm_weak']:>9.4f} {r['pm_diff']:>+9.4f} "
              f"{r['n_lags_strong'] + r['n_lags_weak']:>10}", flush=True)
    p = os.path.join(HERE, f'stg_pm_s{lo}_{hi}.json')
    json.dump(rows, open(p, 'w'))
    print(f"\n{(time.time() - t0) / 60:.1f} min actual")
    print("written", p, f"({os.path.getsize(p) / 1e6:.1f} MB, "
          f"spikes/weights/classes/windows included)")


if __name__ == '__main__':
    main()
