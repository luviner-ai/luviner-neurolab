"""Bootstrap CI on the DUR-2 tie-check gap, de-duplicated pool.

    [EXPLORATORY, POST-HOC -- NOT REGISTERED]

The registered scorer pools `intr + inh`, and the CaS x1.0 baseline rows ARE
DUR-1's g_ie = 0.03 rows, so 5 of its 23 points are exact duplicates of 5
others. This runs on the 18 unique cells, the only defensible pool.

The two correlations are DEPENDENT -- computed on the same rows, against
predictors that themselves correlate at rho = +0.93 -- so the interval comes
from resampling the CELLS, which carries that dependence through, not from
any closed form for independent correlations.

A note on sign, because it is the kind of thing that becomes an erratum.
Both correlations are negative. The signed difference rho_spike - rho_dur is
POSITIVE when duration is the stronger predictor, while the registered
tie-clause compares MAGNITUDES and fires when |rho_spike| >= |rho_dur|. Both
are reported below with the clause's own verdict spelled out.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dur2 import _load                                       # noqa: E402

N_BOOT = 10000


def spearman(x, y):
    if np.unique(x).size < 3 or np.unique(y).size < 3:
        return np.nan
    return float(np.corrcoef(np.argsort(np.argsort(x)),
                             np.argsort(np.argsort(y)))[0, 1])


def pool_unique():
    base, new, gie = _load()
    b = {r['seed']: r for r in base if r['guard_ok']}
    pool = [r for r in new if r['guard_ok']] + [b[s] for s in b] + list(gie)
    seen, out = set(), []
    for r in pool:
        k = (round(r['duration_ms'], 6), round(r['cf_diff'], 9),
             round(r['spikes_per_burst'], 4))
        if k not in seen:
            seen.add(k); out.append(r)
    return pool, out


def main():
    full, uniq = pool_unique()
    d = np.array([r['duration_ms'] for r in uniq])
    s = np.array([r['spikes_per_burst'] for r in uniq])
    y = np.array([r['cf_diff'] for r in uniq])
    n = d.size
    rd, rs = spearman(d, y), spearman(s, y)
    print(f"pool: {len(full)} rows as the scorer builds it, {n} unique\n")
    print(f"  rho(duration, sorting)        {rd:+.4f}")
    print(f"  rho(spikes/burst, sorting)    {rs:+.4f}")
    print(f"  rho(duration, spikes/burst)   {spearman(d, s):+.4f}"
          f"   <- the two predictors are nearly the same variable")
    print(f"  signed difference rho_spk - rho_dur   {rs - rd:+.4f}")
    print(f"  magnitude difference |rho_spk| - |rho_dur|   {abs(rs) - abs(rd):+.4f}")

    rng = np.random.default_rng(0)
    diff, mag = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        a, b_ = spearman(d[i], y[i]), spearman(s[i], y[i])
        if np.isnan(a) or np.isnan(b_):
            continue
        diff.append(b_ - a); mag.append(abs(b_) - abs(a))
    diff, mag = np.array(diff), np.array(mag)
    lo, hi = np.percentile(diff, [2.5, 97.5])
    mlo, mhi = np.percentile(mag, [2.5, 97.5])
    print(f"\n  bootstrap over cells, {len(diff)} usable of {N_BOOT} resamples")
    print(f"    signed    rho_spk - rho_dur      {rs - rd:+.4f}"
          f"   95% CI [{lo:+.4f}, {hi:+.4f}]")
    print(f"    magnitude |rho_spk| - |rho_dur|  {abs(rs) - abs(rd):+.4f}"
          f"   95% CI [{mlo:+.4f}, {mhi:+.4f}]")
    print(f"    fraction of resamples where the tie-clause would fire "
          f"(|rho_spk| >= |rho_dur|): {float((mag >= 0).mean()):.3f}")
    verdict = ("INDISTINGUISHABLE -- the interval spans zero"
               if lo <= 0 <= hi else "separated")
    print(f"\n  {verdict}")


if __name__ == '__main__':
    main()
