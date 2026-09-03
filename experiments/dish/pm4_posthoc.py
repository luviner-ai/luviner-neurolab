"""PM-4 post-hoc: continuous weight, and the kernel ablation.

    [EXPLORATORY, POST-HOC -- NOT REGISTERED]

Neither analysis here was in PM-4's registration. PM-4 registered a MEDIAN
SPLIT into strong and weak halves and scored the direction of the gap; the
continuous-weight analysis below replaces that split with a graded predictor
CHOSEN AFTER the results were seen, and the 25/75 quantile contrast is
likewise a post-hoc choice. They are reported as exploratory throughout, and
neither can promote itself into evidence at the registered bar.

No simulation is run. Everything is recomputed from the spikes, burst
windows, wiring and weights persisted by `stg_pm4.py`.

The first thing this file does is REPRODUCE the persisted `pm_diff` and
`cf_diff` from those raw fields. If the recomputation disagreed with what
the original run wrote, every number below would be measured with a
different instrument than the one that produced the registered result --
which is precisely how MEA-1's first pass went wrong.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stg_pm import MAX_LAG, TAU_PLUS, TAU_MINUS      # noqa: E402


def cohort():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'stg_pm4_s*.json'))):
        rows += json.load(open(f))
    return sorted(rows, key=lambda r: r['seed'])


def per_synapse(r):
    """Per-synapse causal/anti counts and kernel-weighted P/M.

    The lag convention is copied from `stg_pm4.measure`, including the
    detail that lag == 0 counts as causal (`pos = lag[lag >= 0]`). Changing
    that would silently shift every fraction.
    """
    spikes = [np.asarray(s) for s in r['spikes']]
    pre, post = np.asarray(r['pre']), np.asarray(r['post'])
    n = pre.size
    caus = np.zeros(n); anti = np.zeros(n)
    P = np.zeros(n); M = np.zeros(n)
    for a, b in r['windows']:
        cache = {}

        def sp(c):
            if c not in cache:
                s = spikes[c]
                cache[c] = s[(s >= a) & (s <= b)]
            return cache[c]
        for j in range(n):
            pre_s, post_s = sp(pre[j]), sp(post[j])
            if pre_s.size == 0 or post_s.size == 0:
                continue
            lag = (post_s[None, :] - pre_s[:, None]).ravel()
            lag = lag[np.abs(lag) <= MAX_LAG]
            if lag.size == 0:
                continue
            pos, neg = lag[lag >= 0], lag[lag < 0]
            caus[j] += pos.size; anti[j] += neg.size
            P[j] += np.exp(-pos / TAU_PLUS).sum()
            M[j] += np.exp(neg / TAU_MINUS).sum()
    return caus, anti, P, M


def rank(x):
    return np.argsort(np.argsort(np.asarray(x, float)))


def spearman(x, y):
    return float(np.corrcoef(rank(x), rank(y))[0, 1])


def main():
    rows = cohort()
    out = []
    print("reproduction check -- recomputed vs persisted, per wiring")
    print(f"{'seed':>5} {'cf_diff (rec)':>14} {'cf_diff (persisted)':>20} "
          f"{'pm_diff (rec)':>14} {'pm_diff (persisted)':>20}")
    for r in rows:
        caus, anti, P, M = per_synapse(r)
        st = np.asarray(r['strong'], bool)
        cf = {k: caus[m].sum() / (caus[m].sum() + anti[m].sum())
              for k, m in ((True, st), (False, ~st))}
        pm = {k: P[m].sum() / M[m].sum() for k, m in ((True, st), (False, ~st))}
        cfd, pmd = cf[True] - cf[False], pm[True] - pm[False]
        print(f"{r['seed']:5d} {cfd:14.8f} {r['cf_diff']:20.8f} "
              f"{pmd:14.6f} {r['pm_diff']:20.6f}")
        out.append({'seed': r['seed'], 'caus': caus, 'anti': anti, 'P': P,
                    'M': M, 'w': np.asarray(r['weights']), 'strong': st,
                    'cfd': cfd, 'pmd': pmd,
                    'cfd_persisted': r['cf_diff'], 'pmd_persisted': r['pm_diff']})
    dmax = max(abs(o['cfd'] - o['cfd_persisted']) for o in out)
    print(f"\nlargest cf_diff discrepancy: {dmax:.2e}"
          f"  -> {'SAME INSTRUMENT' if dmax < 1e-9 else 'DIFFERENT INSTRUMENT -- STOP'}")
    if dmax >= 1e-9:
        sys.exit(1)

    # ---- (a) continuous weight [EXPLORATORY, POST-HOC] ------------------
    print("\n(a) weight as a continuous predictor  [EXPLORATORY, POST-HOC]")
    print(f"{'seed':>5} {'n syn':>6} {'rho(w, causal frac)':>21} "
          f"{'q25 frac':>10} {'q75 frac':>10} {'q75-q25 (pp)':>13}")
    rhos, gaps = [], []
    for o in out:
        tot = o['caus'] + o['anti']
        m = tot > 0
        w, cf = o['w'][m], o['caus'][m] / tot[m]
        rho = spearman(w, cf)
        lo, hi = np.percentile(w, 25), np.percentile(w, 75)
        # the 25/75 contrast is a post-hoc choice, declared as such
        g = cf[w >= hi].mean() - cf[w <= lo].mean()
        rhos.append(rho); gaps.append(g)
        print(f"{o['seed']:5d} {m.sum():6d} {rho:21.3f} "
              f"{cf[w <= lo].mean():10.5f} {cf[w >= hi].mean():10.5f} {g * 100:13.3f}")
    rhos, gaps = np.array(rhos), np.array(gaps)
    print(f"  mean rho {rhos.mean():+.3f} (sd {rhos.std(ddof=1):.3f}), "
          f"positive in {(rhos > 0).sum()}/{len(rhos)}")
    print(f"  mean q75-q25 gap {gaps.mean() * 100:+.3f} pp, "
          f"positive in {(gaps > 0).sum()}/{len(gaps)}")
    print(f"  registered median split, for comparison: "
          f"mean cf_diff {np.mean([o['cfd'] for o in out]) * 100:+.3f} pp, "
          f"positive in {sum(1 for o in out if o['cfd'] > 0)}/{len(out)}")

    # ---- (b) kernel ablation [EXPLORATORY, POST-HOC] --------------------
    print("\n(b) flat count vs exponential kernel  [EXPLORATORY, POST-HOC]")
    cfd = np.array([o['cfd'] for o in out]); pmd = np.array([o['pmd'] for o in out])
    ratio = pmd / cfd
    print(f"  RECORD's published quantity: CV of pm_diff/cf_diff across wirings"
          f" = {ratio.std(ddof=1) / ratio.mean():.3f}  (record: 0.13)")
    print(f"  rho(cf_diff, pm_diff) across wirings = {spearman(cfd, pmd):+.3f}")
    print(f"  sign agreement: flat {int((cfd > 0).sum())}/10, "
          f"exponential {int((pmd > 0).sum())}/10, "
          f"same wirings: {bool(np.all((cfd > 0) == (pmd > 0)))}")
    ps = []
    for o in out:
        m = (o['caus'] + o['anti']) > 0
        m &= o['M'] > 0
        ps.append(spearman(o['caus'][m] / (o['caus'][m] + o['anti'][m]),
                           o['P'][m] / o['M'][m]))
    ps = np.array(ps)
    print(f"  ADDITIONAL check, explicitly a different quantity from the"
          f" record's:\n    per-synapse rho(count fraction, P/M), mean over"
          f" wirings = {ps.mean():+.3f} (min {ps.min():+.3f})")
    np.savez(os.path.join(HERE, 'pm4_posthoc.npz'),
             rho_weight=rhos, gap_2575=gaps, cfd=cfd, pmd=pmd,
             rho_kernel=ps,
             w=np.concatenate([o['w'] for o in out]),
             cf=np.concatenate([o['caus'] / np.maximum(o['caus'] + o['anti'], 1)
                                for o in out]),
             tot=np.concatenate([o['caus'] + o['anti'] for o in out]),
             seed=np.concatenate([np.full(o['w'].size, o['seed']) for o in out]))
    print("\nwrote pm4_posthoc.npz")


if __name__ == '__main__':
    main()
