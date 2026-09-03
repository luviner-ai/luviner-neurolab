"""PM-4: the accumulation reading, direction-only, on data that did not exist.

`## Run: stg-P/M` found the strong-minus-weak `P/M` difference positive
in 10/10 wirings but failed its own second clause: the between-wiring
spread slightly exceeded the mean. `## Run: PM-3` then established WHY,
before this registration was written: the two regimes are two outcomes of
the same learning process, a wiring can convert between them mid-run, and
so the variance is intrinsic to the dynamics rather than a nuisance
factor any stratification could block.

**That is the entire justification for dropping the magnitude clause,
and the order matters.** The clause is removed because its failure is
explained by a measured mechanism, established beforehand, on a
prediction registered in advance — not because it was inconvenient. A
bar may change between registrations; it may not change after seeing
which way it decides. The compensation is that this runs on **wirings
that do not yet exist** and is **all-or-nothing**.

Registered, on ten FRESH wirings (seeds 10-19, never simulated)
----------------------------------------------------------------
    (1) strong-minus-weak `P/M` positive in >= 9/10
    (2) causal-pair excess on strong synapses positive in >= 9/10
    (3) the `P/M` difference tracks the causal-fraction difference
        across wirings, Spearman rho >= 0.6

**All three, or the accumulation reading dies.** No partial credit, no
clause reported without the others, no seed excluded.

(2) and (3) were *observations* in stg-P/M, found after scoring. Here
they are predictions on unseen data for the first time. That is what
makes a pass mean something and a miss fatal.

**Reported but NEVER scored**, so the descriptive record cannot become a
retroactive criterion: the magnitude of every difference, the
between-wiring spread, and each wiring's burst count as its regime
marker. Visible, never load-bearing.

The wiring is persisted explicitly, and spikes, weights, classes,
windows and causal counts all go to one file.

Cost: 3.7 s wall per biological second, 10 wirings x 60 s = 37 min in
batches under the line.

Usage:  /tmp/venv/bin/python experiments/dish/stg_pm4.py --seeds 10 13
        /tmp/venv/bin/python experiments/dish/stg_pm4.py --score
"""

import glob
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
from luviner.biophysics.dish import rule_with_integral
from luviner.biophysics.stg_dish import STGDish, STGDriver
from stg_pm import (burst_windows, spearman, N_EXC, N_INH, CONN, W_REC,
                    G_IE, SECONDS, MAX_LAG, TAU_PLUS, TAU_MINUS, VALIDATION)


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

    return {
        'seed': seed, 'n_bursts': len(windows),
        'dur_ms': float(widths.mean()) if widths.size else 0.0,
        'pm_strong': pm[True], 'pm_weak': pm[False],
        'pm_diff': pm[True] - pm[False],
        'cf_strong': cf[True], 'cf_weak': cf[False],
        'cf_diff': cf[True] - cf[False],
        'capped': int(d.capped.sum()),
        # persisted explicitly: the wiring, not a way to rebuild it
        'pre': d.pre_e.tolist(), 'post': d.post_e.tolist(),
        'weights': w.tolist(), 'strong': strong.tolist(),
        'windows': [list(x) for x in windows],
        'spikes': [s.tolist() for s in spikes],
    }


def score():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'stg_pm4_s*.json'))):
        rows += json.load(open(f))
    rows.sort(key=lambda r: r['seed'])
    if len(rows) < 10:
        print(f"only {len(rows)}/10 wirings present -- not scoring a "
              f"partial set")
        return
    pm = np.array([r['pm_diff'] for r in rows])
    cf = np.array([r['cf_diff'] for r in rows])
    print(f"\n{'seed':>4} {'bursts':>7} {'P/M diff':>10} {'causal diff':>12}")
    for r in rows:
        print(f"{r['seed']:>4} {r['n_bursts']:>7} {r['pm_diff']:>+10.4f} "
              f"{r['cf_diff']:>+12.5f}")
    rho = spearman(pm, cf)
    p1, p2 = int((pm > 0).sum()), int((cf > 0).sum())
    print(f"\n=== SCORED, all-or-nothing ===")
    print(f"  (1) P/M diff positive        {p1}/10   "
          f"{'MET' if p1 >= 9 else 'NOT MET'}   (>= 9 required)")
    print(f"  (2) causal excess positive   {p2}/10   "
          f"{'MET' if p2 >= 9 else 'NOT MET'}   (>= 9 required)")
    print(f"  (3) Spearman rho(P/M, causal) {rho:+.3f}   "
          f"{'MET' if rho >= 0.6 else 'NOT MET'}   (>= 0.60 required)")
    ok = p1 >= 9 and p2 >= 9 and rho >= 0.6
    print(f"\n  -> accumulation reading {'ESTABLISHED' if ok else 'DEAD'}")
    print(f"\n--- descriptive only, NEVER scored ---")
    print(f"  P/M magnitude {pm.mean():+.4f} +/- {pm.std(ddof=1):.4f}   "
          f"effect/spread {pm.mean() / pm.std(ddof=1):.2f}")
    nb = np.array([r['n_bursts'] for r in rows])
    print(f"  burst counts {sorted(nb.tolist())}")


def main():
    if '--score' in sys.argv:
        score()
        return
    i = sys.argv.index('--seeds')
    lo, hi = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    seeds = list(range(lo, hi + 1))
    assert lo >= 10, "seeds 0-9 are not fresh; PM-4 runs on unseen wirings"
    est = len(seeds) * SECONDS * 3.7
    rule = rule_with_integral(-0.134)
    print(f"PM-4, FRESH seeds {lo}-{hi}: {len(seeds)} wirings x "
          f"{SECONDS:.0f} s\nestimated {est / 60:.1f} min   "
          f"neutral ratio {rule.A_minus / rule.A_plus:.4f}", flush=True)
    if est > 15 * 60:
        print("over the line for one launch -- smaller batches")
        return
    t0 = time.time()
    print(f"\n{'seed':>4} {'bursts':>7} {'dur':>6} {'P/M diff':>10} "
          f"{'causal diff':>12} {'capped':>7}")
    rows = []
    for s in seeds:
        r = one(s)
        rows.append(r)
        print(f"{s:>4} {r['n_bursts']:>7} {r['dur_ms']:>6.0f} "
              f"{r['pm_diff']:>+10.4f} {r['cf_diff']:>+12.5f} "
              f"{r['capped']:>7}", flush=True)
    p = os.path.join(HERE, f'stg_pm4_s{lo}_{hi}.json')
    json.dump(rows, open(p, 'w'))
    print(f"\n{(time.time() - t0) / 60:.1f} min actual")
    print("written", p, f"({os.path.getsize(p) / 1e6:.1f} MB)")


if __name__ == '__main__':
    main()
