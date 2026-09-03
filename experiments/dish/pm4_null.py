"""(d) A per-synapse burst-shuffle null for the carrier IN SIMULATION.

    [EXPLORATORY, POST-HOC -- NOT REGISTERED]

Sections 3.1-3.4 scored the causal-pair excess against the strong/weak
contrast but never against a surrogate; section 3.5 (MEA-1) did. That
asymmetry is a real gap and this closes it, using the SAME null MEA-1
registered before its download:

  A's spikes from burst i are paired against B's from burst j, i != j.
  Spike count, firing rate and each unit's typical position within a burst
  are preserved; only same-burst coordination is destroyed.

The exact mean over all i != j pairings needs no sampling -- pool the
onset-aligned spikes without the gap and every (i, j) combination is counted
at once, then remove the i == j diagonal, which is the observed count. The
derangements are run only for the spread.

The null also subtracts any static phase shift a synapse imposes on every
burst, so the excess measured here is a LOWER BOUND on same-burst
coordination: it can miss true signal but cannot manufacture it.

No simulation. Recomputed from what stg_pm4.py persisted.
"""
import json
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stg_pm import MAX_LAG                                    # noqa: E402

R_DERANGE = 50
GAP = 1e6


def _count(A, B, W):
    if A.size == 0 or B.size == 0:
        return 0, 0
    c = np.searchsorted(B, A + W, 'right') - np.searchsorted(B, A, 'left')
    d = np.searchsorted(A, B + W, 'right') - np.searchsorted(A, B, 'left')
    return int(c.sum()), int(d.sum())


def _offset(sc, ks, perm=None):
    src = ks if perm is None else perm
    xs = [sc[k] + i * GAP for i, k in enumerate(src) if sc[k].size]
    return np.concatenate(xs) if xs else np.empty(0)


def one(r, rng):
    spikes = [np.asarray(s) for s in r['spikes']]
    pre, post = np.asarray(r['pre']), np.asarray(r['post'])
    strong = np.asarray(r['strong'], bool)
    wins = r['windows']
    ks = list(range(len(wins)))
    S = {c: [np.sort(spikes[c][(spikes[c] >= a) & (spikes[c] <= b)] - a)
             for a, b in wins]
         for c in set(pre.tolist() + post.tolist())}
    obs, null, sd, keep = [], [], [], []
    for j in range(pre.size):
        A, B = _offset(S[pre[j]], ks), _offset(S[post[j]], ks)
        co, ao = _count(A, B, MAX_LAG)
        if co + ao == 0:
            continue
        Ap = np.sort(np.concatenate([S[pre[j]][k] for k in ks]))
        Bp = np.sort(np.concatenate([S[post[j]][k] for k in ks]))
        ct, at = _count(Ap, Bp, MAX_LAG)
        cn, an = ct - co, at - ao
        if cn + an <= 0:
            continue
        fr = []
        for _ in range(R_DERANGE):
            p = list(ks)
            for _ in range(50):
                rng.shuffle(p)
                if all(x != y for x, y in zip(p, ks)):
                    break
            c2, a2 = _count(A, _offset(S[post[j]], ks, p), MAX_LAG)
            if c2 + a2:
                fr.append(c2 / (c2 + a2))
        obs.append(co / (co + ao)); null.append(cn / (cn + an))
        sd.append(float(np.std(fr)) if fr else np.nan); keep.append(j)
    obs, null, sd = np.array(obs), np.array(null), np.array(sd)
    st = strong[np.array(keep)]
    exc = obs - null
    return {'seed': r['seed'], 'n': int(obs.size),
            'excess_strong': float(exc[st].mean()),
            'excess_weak': float(exc[~st].mean()),
            'excess_diff': float(exc[st].mean() - exc[~st].mean()),
            'raw_diff': float(obs[st].mean() - obs[~st].mean()),
            'null_diff': float(null[st].mean() - null[~st].mean()),
            'mean_null_sd': float(np.nanmedian(sd)),
            'z_diff': float((exc[st].mean() - exc[~st].mean()) /
                            (np.nanmedian(sd) * np.sqrt(1 / st.sum()
                                                        + 1 / (~st).sum())))}


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'stg_pm4_s*.json'))):
        rows += json.load(open(f))
    rows.sort(key=lambda r: r['seed'])
    rng = np.random.default_rng(0)
    out = [one(r, rng) for r in rows]
    print("(d) per-synapse burst-shuffle null  [EXPLORATORY, POST-HOC]")
    print(f"{'seed':>5} {'n syn':>6} {'raw diff':>10} {'null diff':>10} "
          f"{'excess-over-null':>17} {'z':>7}")
    for o in out:
        print(f"{o['seed']:5d} {o['n']:6d} {o['raw_diff'] * 100:+10.4f} "
              f"{o['null_diff'] * 100:+10.4f} {o['excess_diff'] * 100:+17.4f} "
              f"{o['z_diff']:+7.1f}")
    e = np.array([o['excess_diff'] for o in out])
    raw = np.array([o['raw_diff'] for o in out])
    nul = np.array([o['null_diff'] for o in out])
    print(f"\n  raw strong-weak excess       mean {raw.mean() * 100:+.4f} pp, "
          f"positive {int((raw > 0).sum())}/10")
    print(f"  the null's own strong-weak   mean {nul.mean() * 100:+.4f} pp, "
          f"positive {int((nul > 0).sum())}/10")
    print(f"  excess OVER the null         mean {e.mean() * 100:+.4f} pp, "
          f"positive {int((e > 0).sum())}/10")
    print(f"  registered value, for scale: +0.346 pp (mean cf_diff)")
    json.dump(out, open(os.path.join(HERE, 'pm4_null.json'), 'w'))


if __name__ == '__main__':
    main()
