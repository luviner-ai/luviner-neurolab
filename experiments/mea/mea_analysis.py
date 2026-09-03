"""MEA-1: the consolidation prediction, tested on real cortical tissue.

Data: Wagenaar, Pine & Potter (2006), BMC Neuroscience 7:11 -- dissociated
rat cortical cultures on 59-electrode MEAs, spontaneous activity. Public
archive; `manifest.tsv` carries the exact URLs and checksums. The raw files
are not vendored here.

What this measures
------------------
`## Run: PM-4` established that the consolidation carrier in the model dish
is an excess of causally ordered spike pairs -- a pure count -- among the
stronger half of synapses. Section 5 of the preprint turned that into a
prediction for MEA cultures. This is that test.

Three things make it a test rather than a tautology, and each was registered
before the download:

1. SPLIT-HALF BY BURST. On an MEA there are no weights: strength must be
   estimated from spikes, and a cross-correlogram peak is essentially a
   causal pair count -- so scoring both on the same events would make the
   prediction true by construction. Strength is estimated on odd bursts and
   the excess scored on even, then swapped and averaged. The spikes that
   define a pair's class are never the spikes that score it.

2. A PER-PAIR NULL, and it is not jitter. Class assignment correlates with
   firing rate, so a pair's chance floor is not 0.5. Jitter cannot supply
   that floor here: the statistic lives at burst scale, so jitter wide
   enough to matter makes both trains uniform and returns exactly 0.5 for
   every pair -- collapsing into the nominal null it was meant to replace --
   while jitter small compared to the window leaves the data untouched.
   Instead each pair is scored against a BURST-SHUFFLE null: A's spikes from
   burst i against B's from burst j, i != j. That preserves spike count,
   rate, and each unit's typical position within a burst, and destroys only
   same-burst coordination. The null also subtracts any synapse-induced
   static phase shift, so the measured excess is a LOWER BOUND on same-burst
   coordination -- it can miss true signal but cannot manufacture it.

3. NO INHERITED WINDOW. The model's 250 ms was chosen for ~450 ms model
   bursts. Here the window is half the median burst duration OF THAT
   RECORDING, which comes out 31-332 ms. Inheriting 250 ms would have
   counted pairs across a window several times the structure it sits in.

Burst detection is the dataset's own published method (SIMMUX), not this
board's detector, and burst DURATION is measured the way its authors measure
it -- on the ASDR profile, between the 20% points around the peak. That
distinction is not cosmetic: taking the burstlet extent as the duration gave
a 1713 ms median at 25 div against the paper's published <200 ms, while the
ASDR measurement gives 186 ms, agreeing with a figure other people
published on the same data.

The first 300 s of every recording is discarded. The paper reports that
moving a dish into the rig raises SYNCHRONY for about five minutes almost
without raising firing rate -- a transient in precisely the quantity here.
"""
import json
import os
import sys
import time

import numpy as np

DROP_HEAD_S = 300.0     # the paper's own mechanical-perturbation transient
BIN_S = 0.001
SIGMA_S = 0.010         # the paper's 10 ms Gaussian
CCG_LAG = 0.025         # strength-proxy window, fixed pre-registration
R_DERANGE = 50          # burst-shuffle repetitions, for the null spread
GAP = 100.0             # keeps bursts from interacting in the aligned frame
MIN_BURSTS = 8          # registered void guard
MIN_SPIKES_PER_UNIT = 100


# --------------------------------------------------------------- loading
def load(path):
    a = np.loadtxt(path)
    t, ch = a[:, 0], a[:, 1].astype(int)
    keep = (ch >= 0) & (ch <= 59) & (t >= DROP_HEAD_S)   # 60 = stimulus marker
    return t[keep], ch[keep]


# ------------------------------------------------- SIMMUX burst detection
def burstlets(t, ch):
    """The dataset's own definition, quoted from its Methods: >=4 spikes on
    one electrode with every ISI below 1/4 of that electrode's inverse mean
    rate, or below 100 ms if that rate is under 10 Hz."""
    out = []
    span = t.max() - t.min()
    for c in np.unique(ch):
        ts = np.sort(t[ch == c])
        if ts.size < 4:
            continue
        rate = ts.size / span
        thr = 0.100 if rate < 10.0 else 0.25 / rate
        ok = np.diff(ts) < thr
        i = 0
        while i < ok.size:
            if not ok[i]:
                i += 1
                continue
            j = i
            while j < ok.size and ok[j]:
                j += 1
            if (j - i + 1) >= 4:
                out.append((ts[i], ts[j], c))
            i = j + 1
    return sorted(out)


def merge(bl):
    """'Any group of burstlets across several electrodes that overlapped in
    time was considered a burst.'"""
    bursts = []
    for s, e, c in bl:
        if bursts and s <= bursts[-1][1]:
            bursts[-1][1] = max(bursts[-1][1], e)
            bursts[-1][2].add(c)
        else:
            bursts.append([s, e, {c}])
    return bursts


def _gauss(x, sigma_bins):
    n = int(4 * sigma_bins)
    k = np.exp(-0.5 * (np.arange(-n, n + 1) / sigma_bins) ** 2)
    return np.convolve(x, k / k.sum(), mode='same')


def windows(t, bursts, baseline):
    """One window per SIMMUX region, spanning the 20% points either side of
    that region's dominant ASDR peak -- the paper's own convention.

    An earlier version split on every threshold crossing instead. It agreed
    with this function on the one recording used to validate it (186 vs
    185 ms) and disagreed by 10x on the rest, fragmenting one burst into ten
    and collapsing the lag window from ~90 ms to ~13 ms. That run is void.
    Validate the code that runs, not a sibling of it.
    """
    out = []
    for s, e, _ in bursts:
        lo, hi = s - 0.5, e + 0.5
        nb = int((hi - lo) / BIN_S) + 1
        idx = ((t[(t >= lo) & (t < hi)] - lo) / BIN_S).astype(int)
        a = _gauss(np.bincount(idx, minlength=nb).astype(float) / BIN_S,
                   SIGMA_S / BIN_S)
        pk = int(a.argmax())
        thr20 = baseline + 0.20 * (a[pk] - baseline)
        L = pk - int(np.argmax(a[pk::-1] < thr20)) if (a[:pk + 1] < thr20).any() else 0
        R = pk + int(np.argmax(a[pk:] < thr20)) if (a[pk:] < thr20).any() else nb - 1
        if R > L:
            out.append((lo + L * BIN_S, lo + R * BIN_S))
    return sorted(out)


# ------------------------------------------------------- pair statistics
def per_burst(t, ch, wins, nch=60):
    S = [[None] * len(wins) for _ in range(nch)]
    for k, (s, e) in enumerate(wins):
        m = (t >= s) & (t <= e)
        tt, cc = t[m] - s, ch[m]
        for c in range(nch):
            S[c][k] = np.sort(tt[cc == c])
    return S


def _count(A, B, W):
    """(#pairs with 0 < b-a <= W, #pairs with 0 < a-b <= W). Counted, never
    enumerated -- which is why a 10^5-pair statistic costs O(n log n)."""
    if A.size == 0 or B.size == 0:
        return 0, 0
    c = np.searchsorted(B, A + W, 'right') - np.searchsorted(B, A, 'right')
    d = np.searchsorted(A, B + W, 'right') - np.searchsorted(A, B, 'right')
    return int(c.sum()), int(d.sum())


def _offset(S_c, ks, perm=None):
    src = ks if perm is None else perm
    xs = [S_c[k] + i * GAP for i, k in enumerate(src) if S_c[k].size]
    return np.concatenate(xs) if xs else np.empty(0)


def causal(S, a, b, ks, W, rng):
    """Observed causal fraction, the exact burst-shuffle null mean, and the
    null's spread.

    The exact mean over ALL burst pairings i != j needs no sampling: pool the
    onset-aligned spikes without the gap and every (i, j) combination is
    counted at once, then remove the i == j diagonal, which is the observed
    count. The derangements are run only for the spread.
    """
    A, B = _offset(S[a], ks), _offset(S[b], ks)
    co, ao = _count(A, B, W)
    if co + ao == 0:
        return None
    Ap = np.sort(np.concatenate([S[a][k] for k in ks] or [np.empty(0)]))
    Bp = np.sort(np.concatenate([S[b][k] for k in ks] or [np.empty(0)]))
    ct, at = _count(Ap, Bp, W)
    cn, an = ct - co, at - ao
    if cn + an <= 0:
        return None
    null_mean = cn / (cn + an)
    fr = []
    for _ in range(R_DERANGE):
        p = list(ks)
        for _ in range(50):
            rng.shuffle(p)
            if all(x != y for x, y in zip(p, ks)):
                break
        c2, a2 = _count(A, _offset(S[b], ks, p), W)
        if c2 + a2:
            fr.append(c2 / (c2 + a2))
    return co / (co + ao), null_mean, (float(np.std(fr)) if fr else float('nan'))


def ccg_peak(S, a, b, ks, rate_a, rate_b):
    """Strength proxy: coincidences within 25 ms inside burst windows,
    normalised by the product of the two units' in-burst counts."""
    A, B = _offset(S[a], ks), _offset(S[b], ks)
    if A.size == 0 or B.size == 0 or rate_a <= 0 or rate_b <= 0:
        return 0.0
    c, d = _count(A, B, CCG_LAG)
    return (c + d) / (rate_a * rate_b)


# ------------------------------------------------------------------- run
def run(path, div, culture):
    T = {}
    t0 = time.time(); t, ch = load(path); T['load'] = time.time() - t0

    t0 = time.time()
    asdr = np.bincount(((t - t.min()) / 1.0).astype(int)) / 1.0
    baseline = float(np.median(asdr))
    regions = [x for x in merge(burstlets(t, ch)) if len(x[2]) >= 5]
    wins = windows(t, regions, baseline)
    T['detect'] = time.time() - t0

    base = os.path.basename(path)
    if len(wins) < MIN_BURSTS:
        return {'file': base, 'div': div, 'culture': culture,
                'void': f'fewer than {MIN_BURSTS} bursts', 'n_bursts': len(wins)}

    active = [c for c in range(60) if (ch == c).sum() >= MIN_SPIKES_PER_UNIT]
    S = per_burst(t, ch, wins)
    K = len(wins)
    durs = np.array([e - s for s, e in wins])
    W = float(np.median(durs) / 2.0)          # per recording, never inherited
    rng = np.random.default_rng(0)
    counts = {c: sum(S[c][k].size for k in range(K)) for c in active}
    odd = [k for k in range(K) if k % 2 == 1]
    even = [k for k in range(K) if k % 2 == 0]

    t0 = time.time()
    pairs, lat_rho = [], []
    for half, (est, sco) in enumerate(((odd, even), (even, odd))):
        st = {}
        for i, a in enumerate(active):
            for b in active[i + 1:]:
                st[(a, b)] = ccg_peak(S, a, b, est, counts[a], counts[b])
        med = float(np.median(list(st.values())))
        for (a, b), s in st.items():
            r = causal(S, a, b, sco, W, rng)
            if r is None or not np.isfinite(r[1]):
                continue
            pairs.append({'half': half, 'a': a, 'b': b,
                          'strong': bool(s > med), 'strength': s,
                          'obs': r[0], 'null': r[1], 'null_sd': r[2],
                          'excess': r[0] - r[1]})
        # P3 keeps the same split-half discipline: strength from the
        # estimation half, latency from the scoring half.
        inc = {u: float(np.mean([v for (x, y), v in st.items() if u in (x, y)]))
               for u in active}
        lat = {}
        for u in active:
            fs = [S[u][k][0] for k in sco if S[u][k].size]
            if len(fs) >= 5:
                lat[u] = float(np.median(fs))
        us = [u for u in active if u in lat]
        if len(us) >= 8:
            rx = np.argsort(np.argsort([lat[u] for u in us]))
            ry = np.argsort(np.argsort([inc[u] for u in us]))
            lat_rho.append(float(np.corrcoef(rx, ry)[0, 1]))
    T['pairs'] = time.time() - t0

    strong = np.array([p['excess'] for p in pairs if p['strong']])
    weak = np.array([p['excess'] for p in pairs if not p['strong']])
    sd = np.array([p['null_sd'] for p in pairs if np.isfinite(p['null_sd'])])
    return {
        'file': base, 'div': div, 'culture': culture,
        'n_bursts': K, 'W_ms': W * 1e3,
        'median_dur_ms': float(np.median(durs) * 1e3),
        'active_units': len(active),
        'mean_rate_hz': float(t.size / (t.max() - t.min()) / len(active)),
        'strong_excess': float(strong.mean()), 'weak_excess': float(weak.mean()),
        'sorting': float(strong.mean() - weak.mean()),
        'n_strong': int(strong.size), 'n_weak': int(weak.size),
        # the per-pair null spread, propagated to the recording level: the
        # chance floor against which the sorting statistic must be read
        'null_sd_median': float(np.median(sd)),
        'sorting_null_sd': float(np.sqrt(np.mean(sd ** 2) *
                                         (1 / max(strong.size, 1) +
                                          1 / max(weak.size, 1)))),
        'latency_rho': float(np.mean(lat_rho)) if lat_rho else None,
        'timing_s': T, 'pairs': pairs,
    }


if __name__ == '__main__':
    path, div, culture, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    r = run(path, div, culture)
    with open(out, 'w') as f:
        json.dump(r, f)
    print(f"{r['file']:20s} {sum(r.get('timing_s', {}).values()):6.1f}s"
          f"  bursts {r.get('n_bursts', 0):5d}"
          + ("  VOID" if 'void' in r else f"  sorting {r['sorting']:+.5f}"),
          flush=True)
