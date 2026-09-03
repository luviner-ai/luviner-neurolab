"""DUR-2: triangulate the duration result with a lever that touches no synapse.

`## Run: DUR-1` shortened bursts with `g_ie` and sorting rose 6-10x in
5/5 wirings. But `g_ie` is inhibition, and scoping said plainly that no
knob in this preparation moves duration alone. This arm shortens the
burst from **inside the cell** -- `g_CaS`, from the validated Prinz set --
so it carries a different confound, and the two arms triangulate.

Why `g_CaS` and not `g_KCa`: measured on an isolated cell, `g_KCa x2.0`
gives 156 ms but drags rate +43%, while `g_CaS x0.6` gives 209 ms at
unchanged rate. `g_KCa` would have smuggled rate back in.

The rate profiles differ, and that is the point
-----------------------------------------------
`g_ie` held rate nearly fixed (8.6% spread) while cutting duration 41%.
`g_CaS` in the network cuts duration ~56% while dropping rate ~25%. **If
sorting tracked rate rather than duration, the two arms would disagree.**
A triangulating lever that shared the confound would buy nothing.

The ceiling, registered before the run
--------------------------------------
Spikes per burst falls with duration under BOTH levers -- 1095 -> 191
here, 26,736 -> 6,458 in DUR-1 wiring 30. They are collinear in both
arms, so this design separates duration from **inhibition** and from
**rate**, and NOT from **spikes per burst**. A pass licenses "sorting
tracks burst duration, or the spike count that falls with it". That is
the honest maximum and it is written before the numbers exist.

Registered
----------
  P1  within each wiring, sorting is stronger at CaS x0.6 than at x1.0,
      in >= 4 of 5
  P2  pooling both arms against REALIZED duration, rho <= -0.5, **and
      the relationship holds within each arm separately** -- lever
      identity is a registered covariate, so a pooled correlation made of
      two clusters at different offsets cannot pass
  P3  within the intrinsic arm, where rate falls WITH duration, sorting
      still rises as duration falls -- rate cannot account for it

Falsifiers
----------
  * the intrinsic arm does not strengthen sorting -> **"duration" dies,
    and DUR-1's result belongs to inhibition.** In those words.
  * P2 passes pooled but fails within an arm -> the correlation is two
    clusters, not a relationship, reported as an artifact
  * sorting tracks spikes/burst better than realized duration -> nothing
    is separated; reported as unresolved, not as support for either
  * **guard failure voids that CELL, not that wiring** -- DUR-1's scorer
    voided the wiring and turned 5/5 into a printed refutation. Written
    into the registration so the scorer cannot repeat it.

Baseline is free: DUR-1's `g_ie = 0.03` runs ARE `CaS x1.0`, and all five
wirings were alive there. Cost is 2 new levels x 5 wirings = 37 min.
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
from stg_pm import (MAX_LAG, N_EXC, N_INH, CONN, SECONDS, TAU_MINUS,
                    TAU_PLUS, VALIDATION, W_REC, burst_windows, spearman)

G_IE = 0.03                      # where DUR-1 had all five wirings alive
CAS_LEVELS = (0.8, 0.6)          # x1.0 comes free from DUR-1
CAS_INDEX = 2                    # (Na, CaT, CaS, A, KCa, Kd, H, leak)


def measure(d, seed, tag):
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
    n_spikes = int(sum(len(s) for s in spikes))
    return {
        'seed': seed, 'arm': tag, 'n_bursts': n_burst,
        'rate_hz': n_burst / (SECONDS / 2.0),
        'duration_ms': float(widths.mean()) if widths.size else 0.0,
        'spikes_per_burst': n_spikes / max(n_burst, 1),
        'guard_ok': bool(n_burst >= 3 and (widths.max() <= 0.2 * (t_hi - t_lo)
                                           if widths.size else False)),
        'pm_diff': pm[True] - pm[False], 'cf_diff': cf[True] - cf[False],
        'spikes': [s.tolist() for s in spikes],
        'weights': w.tolist(), 'pre': d.pre_e.tolist(),
        'post': d.post_e.tolist(),
    }


def one(seed, cas_mult):
    sets = json.load(open(VALIDATION))['model']['cell_conductances_mS_cm2']
    g_e, g_i = list(sets['ABPD']), sets['LP']
    g_e[CAS_INDEX] = sets['ABPD'][CAS_INDEX] * cas_mult

    def factory(n):
        cells = ([stg.STGNeuron(list(g_e), name=f'E{i}') for i in range(N_EXC)]
                 + [stg.STGNeuron(list(g_i), name=f'I{i}')
                    for i in range(N_INH)])
        return STGDriver(stg.STGNetwork(cells, []))

    d = STGDish(factory, n_exc=N_EXC, n_inh=N_INH, connectivity=CONN,
                w_rec=W_REC, g_ie=G_IE, budget=(N_EXC - 1) * CONN * W_REC,
                depress=(0.20, 400.0), seed=seed)
    d.run(SECONDS * 1000.0, learn=True)
    r = measure(d, seed, 'intrinsic')
    r['cas_mult'] = cas_mult
    return r


def _load():
    """DUR-1's g_ie=0.03 rows are the CaS x1.0 baseline, free."""
    base = []
    for f in sorted(glob.glob(os.path.join(HERE, 'dur1_s*.json'))):
        for r in json.load(open(f)):
            if abs(r['g_ie'] - 0.03) < 1e-9:
                sp = [np.array(s) for s in r['spikes']]
                base.append({'seed': r['seed'], 'arm': 'inhibitory-baseline',
                             'cas_mult': 1.0,
                             'duration_ms': r['duration_ms'],
                             'rate_hz': r['rate_hz'],
                             'spikes_per_burst':
                                 sum(len(s) for s in sp) / max(r['n_bursts'], 1),
                             'cf_diff': r['cf_diff'], 'pm_diff': r['pm_diff'],
                             'guard_ok': r['guard_ok']})
    new = []
    for f in sorted(glob.glob(os.path.join(HERE, 'dur2_s*.json'))):
        new += [{k: v for k, v in r.items()
                 if k not in ('spikes', 'weights', 'pre', 'post')}
                for r in json.load(open(f))]
    # the g_ie arm, for the pooled test: DUR-1's 0.03 and 0.10 cells
    gie = []
    for f in sorted(glob.glob(os.path.join(HERE, 'dur1_s*.json'))):
        for r in json.load(open(f)):
            if r['guard_ok']:
                sp = [np.array(s) for s in r['spikes']]
                gie.append({'seed': r['seed'], 'arm': 'inhibitory',
                            'duration_ms': r['duration_ms'],
                            'rate_hz': r['rate_hz'],
                            'spikes_per_burst':
                                sum(len(s) for s in sp) / max(r['n_bursts'], 1),
                            'cf_diff': r['cf_diff'], 'pm_diff': r['pm_diff'],
                            'guard_ok': True})
    return base, new, gie


def score():
    base, new, gie = _load()
    if len({r['seed'] for r in new}) < 5:
        print(f"only {len({r['seed'] for r in new})}/5 wirings in the "
              f"intrinsic arm -- not scoring a partial set")
        return
    b = {r['seed']: r for r in base if r['guard_ok']}
    n06 = {r['seed']: r for r in new
           if abs(r['cas_mult'] - 0.6) < 1e-9 and r['guard_ok']}
    n08 = {r['seed']: r for r in new
           if abs(r['cas_mult'] - 0.8) < 1e-9 and r['guard_ok']}

    print(f"\n{'seed':>5} {'arm':>12} {'CaS':>5} {'dur ms':>8} {'rate':>6} "
          f"{'spk/brst':>9} {'causal':>10} {'P/M':>9} {'guard':>6}")
    for s in sorted({r['seed'] for r in new} | set(b)):
        for lab, d in (('baseline', b.get(s)), ('CaS x0.8', n08.get(s)),
                       ('CaS x0.6', n06.get(s))):
            if d is None:
                print(f"{s:>5} {lab:>12} {'':>5} {'VOID or missing':>34}")
                continue
            print(f"{s:>5} {lab:>12} {d.get('cas_mult', 1.0):>5.1f} "
                  f"{d['duration_ms']:>8.0f} {d['rate_hz']:>6.2f} "
                  f"{d['spikes_per_burst']:>9.0f} {d['cf_diff']:>+10.5f} "
                  f"{d['pm_diff']:>+9.4f} {str(d['guard_ok']):>6}")
        print()

    paired = [s for s in b if s in n06]
    p1 = sum(1 for s in paired if n06[s]['cf_diff'] > b[s]['cf_diff'])
    print(f"=== scored ===")
    print(f"  P1  sorting stronger at CaS x0.6 than x1.0   {p1}/{len(paired)}"
          f"   {'MET' if p1 >= 4 else 'NOT MET'}")

    intr = [r for r in new if r['guard_ok']] + [b[s] for s in b]
    inh = [r for r in gie]
    pool = intr + inh
    rho_pool = spearman(np.array([r['duration_ms'] for r in pool]),
                        np.array([r['cf_diff'] for r in pool]))
    rho_i = spearman(np.array([r['duration_ms'] for r in intr]),
                     np.array([r['cf_diff'] for r in intr]))
    rho_h = spearman(np.array([r['duration_ms'] for r in inh]),
                     np.array([r['cf_diff'] for r in inh]))
    ok2 = rho_pool <= -0.5 and rho_i <= -0.5 and rho_h <= -0.5
    print(f"  P2  rho(duration, sorting)  pooled {rho_pool:+.3f}  "
          f"intrinsic {rho_i:+.3f}  inhibitory {rho_h:+.3f}   "
          f"{'MET' if ok2 else 'NOT MET'}")
    if rho_pool <= -0.5 and not (rho_i <= -0.5 and rho_h <= -0.5):
        print(f"      pooled passes, an arm does not -- two clusters, not a "
              f"relationship. Reported as an artifact.")

    ii = sorted(intr, key=lambda r: r['duration_ms'])
    p3 = ii[0]['cf_diff'] > ii[-1]['cf_diff']
    print(f"  P3  within the intrinsic arm, shortest beats longest   "
          f"{ii[0]['duration_ms']:.0f} ms {ii[0]['cf_diff']:+.5f} vs "
          f"{ii[-1]['duration_ms']:.0f} ms {ii[-1]['cf_diff']:+.5f}   "
          f"{'MET' if p3 else 'NOT MET'}")

    rho_spk = spearman(np.array([r['spikes_per_burst'] for r in pool]),
                       np.array([r['cf_diff'] for r in pool]))
    print(f"\n  registered tie-check: rho(spikes/burst, sorting) "
          f"{rho_spk:+.3f} vs rho(duration, sorting) {rho_pool:+.3f}")
    if abs(rho_spk) >= abs(rho_pool):
        print(f"  Spike count tracks sorting at least as well as duration. "
              f"Nothing is\n  separated: reported as UNRESOLVED, not as "
              f"support for either.")
    if p1 < 4:
        print(f"\n  The intrinsic arm did not strengthen sorting. DURATION "
              f"DIES, and\n  DUR-1's result belongs to inhibition.")


def main():
    if '--score' in sys.argv:
        score()
        return
    i = sys.argv.index('--seeds')
    lo, hi = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    seeds = list(range(lo, hi + 1))
    est = len(seeds) * len(CAS_LEVELS) * SECONDS * 3.7
    print(f"DUR-2 intrinsic arm, seeds {seeds} x CaS {CAS_LEVELS} "
          f"at g_ie={G_IE}\nestimated {est / 60:.1f} min", flush=True)
    if est > 15 * 60:
        print("over the line -- smaller batches")
        return
    t0 = time.time()
    print(f"\n{'seed':>5} {'CaS':>5} {'dur ms':>8} {'rate':>6} {'spk/brst':>9} "
          f"{'causal':>10} {'P/M':>9} {'guard':>6}")
    rows = []
    for s in seeds:
        for m in CAS_LEVELS:
            r = one(s, m)
            rows.append(r)
            print(f"{s:>5} {m:>5.1f} {r['duration_ms']:>8.0f} "
                  f"{r['rate_hz']:>6.2f} {r['spikes_per_burst']:>9.0f} "
                  f"{r['cf_diff']:>+10.5f} {r['pm_diff']:>+9.4f} "
                  f"{str(r['guard_ok']):>6}", flush=True)
    p = os.path.join(HERE, f'dur2_s{lo}_{hi}.json')
    json.dump(rows, open(p, 'w'))
    print(f"\n{(time.time() - t0) / 60:.1f} min actual\nwritten {p}")


if __name__ == '__main__':
    main()
