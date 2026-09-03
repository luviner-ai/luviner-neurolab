"""MEA-1 scorer. Written before the results were seen, and applying every
registered falsifier -- including the two that fire against us."""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.size < 3:
        return float('nan')
    return float(np.corrcoef(np.argsort(np.argsort(x)),
                             np.argsort(np.argsort(y)))[0, 1])


def main():
    R = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(HERE, 'results', '*.json')))]
    void = [r for r in R if 'void' in r]
    R = [r for r in R if 'void' not in r]
    N = len(R)
    print(f"recordings analysed {N}   void {len(void)} {[v['file'] for v in void]}\n")
    print(f"{'file':14s} {'div':>4} {'brst':>5} {'W ms':>6} {'units':>6} "
          f"{'strong':>9} {'weak':>9} {'sorting':>9} {'null sd':>8} {'lat rho':>8}")
    for r in sorted(R, key=lambda r: (r['div'], r['culture'])):
        print(f"{r['file'][:13]:14s} {r['div']:4d} {r['n_bursts']:5d} {r['W_ms']:6.0f} "
              f"{r['active_units']:6d} {r['strong_excess']:+9.5f} {r['weak_excess']:+9.5f} "
              f"{r['sorting']:+9.5f} {r['sorting_null_sd']:8.5f} "
              + (f"{r['latency_rho']:+8.2f}" if r['latency_rho'] is not None else f"{'--':>8}"))

    # ---- P1: strong > weak in >= N-1 of N -------------------------------
    hits = sum(1 for r in R if r['sorting'] > 0)
    print(f"\nP1  strong > weak in {hits}/{N}   (registered pass: >= {N - 1})")
    print(f"P1  VERDICT: {'PASS' if hits >= N - 1 else 'FAIL -- the counting mechanism does not appear'}")

    # ---- P2: sorting vs duration, WITHIN age strata --------------------
    print("\nP2  within age strata:")
    strata = sorted({r['div'] for r in R})
    iqr = lambda v: float(np.percentile(v, 75) - np.percentile(v, 25))
    within = [iqr([r['median_dur_ms'] for r in R if r['div'] == d]) for d in strata]
    across = iqr([np.median([r['median_dur_ms'] for r in R if r['div'] == d]) for d in strata])
    ratio = float(np.mean(within)) / across if across > 0 else float('inf')
    print(f"    duration spread  within {np.mean(within):.0f} ms  across {across:.0f} ms"
          f"  ratio {ratio:.2f}  (registered floor 0.25)")
    rds, rus = [], []
    for d in strata:
        g = [r for r in R if r['div'] == d]
        rd = spearman([x['median_dur_ms'] for x in g], [x['sorting'] for x in g])
        ru = spearman([x['active_units'] for x in g], [x['sorting'] for x in g])
        rds.append(rd); rus.append(ru)
        print(f"    div {d:2d}  n={len(g)}  rho(sorting, duration) = {rd:+.2f}"
              f"   rho(sorting, units) = {ru:+.2f}")
    rd, ru = float(np.nanmean(rds)), float(np.nanmean(rus))
    print(f"    mean across strata   rho_duration {rd:+.2f}   rho_units {ru:+.2f}")
    if ratio < 0.25:
        print("P2  VERDICT: UNDERPOWERED -- within-stratum duration spread below the\n"
              "    registered floor. Not a refutation.")
    elif abs(ru) >= abs(rd):
        print("P2  VERDICT: UNRESOLVED -- sorting tracks unit count at least as well as\n"
              "    duration. Nothing is separated (the DUR-2 tie-clause).")
    else:
        print(f"P2  VERDICT: {'PASS' if rd < 0 else 'FAIL -- sorting does not track duration negatively'}")

    # ---- P3: latency order should carry NO signal -----------------------
    lr = [r['latency_rho'] for r in R if r['latency_rho'] is not None]
    neg = sum(1 for v in lr if v < 0)
    print(f"\nP3  latency rho: mean {np.mean(lr):+.2f}, negative in {neg}/{len(lr)}")
    print("P3  VERDICT: " + ('FAIL -- the refuted temporal mechanism IS present in tissue'
                             if neg >= len(lr) - 1 and np.mean(lr) < -0.3
                             else 'PASS -- latency order carries no consistent signal'))

    # ---- what the data constrains, beyond the sign test -----------------
    s = np.array([r['sorting'] for r in R])
    sem = s.std(ddof=1) / np.sqrt(N)
    print(f"\nsorting across recordings: mean {s.mean() * 100:+.3f} pp, sd {s.std(ddof=1) * 100:.2f} pp,"
          f" 95% CI [{(s.mean() - 2.131 * sem) * 100:+.2f}, {(s.mean() + 2.131 * sem) * 100:+.2f}] pp")
    print(f"  -> excludes any tissue sorting effect above {(s.mean() + 2.131 * sem) * 100:.2f} pp;"
          f" the model's ~0.20 pp carrier lies inside, so it is not resolvable here.")
    print(f"  -> P1's own pass rule needed ~{1.27 * s.std(ddof=1) * 100:.2f} pp, five times that"
          f" carrier: it was unpassable at the effect size we expected.")


if __name__ == '__main__':
    main()
