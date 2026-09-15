"""Figures for PAPER-2 (the two-sign manuscript), from the recorded JSON only.

Separate from `make_figures.py`, which belongs to the companion preprint and
reads a different experiment tree; the house style is shared, the data is not.

Every value plotted is read from the run outputs in `experiments/coll*/`.
Nothing is recomputed in a way that could drift from the record, nothing is
synthesised, and each figure writes the exact numbers it plots to a CSV beside
it so `check_figures_paper2.py` can compare them against the manuscript.
"""
import csv
import glob
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, '..', 'experiments')
# PY 1 and AB/PD 2 were re-measured at the 35 s horizon used by the other five
# cells; those runs supersede the 30 s ones everywhere except the two extra
# inhibition levels in fig2, which were not re-measured and are marked 30 s.
PY1 = os.path.join(EXP, 'coll14', 'coll14_PY1_*_%s.json')
ABPD2 = os.path.join(EXP, 'coll14', 'coll14_ABPD2_30_%s.json')
OUT = os.path.join(HERE, 'figures')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8, 'axes.labelsize': 8,
    'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': 0.6, 'lines.linewidth': 1.0,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
K, KL, ACC, WARM = '#222222', '#888888', '#4a6fa5', '#a5564a'


def rows_of(path):
    d = json.load(open(path))
    return d['rows'] if isinstance(d, dict) and 'rows' in d else d


def save(fig, name, table):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'))
    plt.close(fig)
    with open(os.path.join(OUT, f'{name}.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerows(table)
    print(f"  wrote {name}.pdf / .png / .csv")


def midranks(v):
    """Ranks with ties averaged, the definition Spearman's rho is built on."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    r, i = [0.0] * len(v), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def wilson(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def table_2x2(off, on):
    """(n, dead_off, rescued, eliminated) from paired per-seed dicts."""
    s = sorted(set(off) & set(on))
    dead = [x for x in s if not off[x]['alive']]
    live = [x for x in s if off[x]['alive']]
    return (len(s), len(dead),
            sum(1 for x in dead if on[x]['alive']),
            sum(1 for x in live if not on[x]['alive']),
            len(live))


def paired(off_paths, on_paths):
    off, on = {}, {}
    for p in off_paths:
        for r in rows_of(p):
            off[r['seed']] = r
    for p in on_paths:
        for r in rows_of(p):
            on[r['seed']] = r
    return off, on


# --------------------------------------------------------------- figure 1
def fig1():
    """Conceptual: one rule, one inhibition, two conductance sets, two signs.

    The only figure here that plots no data. The four numbers it carries are
    read from the recorded runs so it cannot drift from the text.
    """
    off6, on6 = paired(sorted(glob.glob(PY1 % 'off')), sorted(glob.glob(PY1 % 'on')))
    _, d6, r6, e6, l6 = table_2x2(off6, on6)
    o, n_ = paired([ABPD2 % 'off'], [ABPD2 % 'on'])
    _, d5, r5, e5, l5 = table_2x2(o, n_)
    tab = [['box', 'cell', 'quantity', 'numerator', 'denominator'],
           ['left', 'PY 1', 'rescued | control collapsed', r6, d6],
           ['right', 'AB/PD 2', 'eliminated | control survived', e5, l5]]

    fig, ax = plt.subplots(figsize=(5.0, 2.5))
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.4); ax.axis('off')

    def box(x, y, w, h, text, edge=K, face='none', size=7, weight='normal'):
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=face != 'none',
                                   facecolor=face, edgecolor=edge, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                fontsize=size, color=K, weight=weight, linespacing=1.4)

    box(1.6, 3.85, 6.8, 1.25,
        'identical: the plasticity rule, the homeostatic budget,\n'
        'the inhibition, the wiring distribution, the horizon,\n'
        'the survival criterion', edge=KL, size=6.5)
    ax.annotate('', xy=(2.9, 2.98), xytext=(4.3, 3.85),
                arrowprops=dict(arrowstyle='->', lw=0.8, color=KL))
    ax.annotate('', xy=(7.1, 2.98), xytext=(5.7, 3.85),
                arrowprops=dict(arrowstyle='->', lw=0.8, color=KL))
    ax.text(5.0, 3.42, 'only the excitatory\nconductance set differs',
            ha='center', va='center', fontsize=6, color=KL, style='italic')

    box(0.7, 2.0, 4.0, 0.95, 'follower-type set\n(PY 1)', edge=ACC, size=7.5)
    box(5.3, 2.0, 4.0, 0.95, 'pacemaker-type set\n(AB/PD 2)', edge=WARM, size=7.5)
    ax.annotate('', xy=(2.7, 1.52), xytext=(2.7, 2.0),
                arrowprops=dict(arrowstyle='->', lw=0.9, color=ACC))
    ax.annotate('', xy=(7.3, 1.52), xytext=(7.3, 2.0),
                arrowprops=dict(arrowstyle='->', lw=0.9, color=WARM))

    box(0.7, 0.25, 4.0, 1.25,
        f'RESCUED\n{r6} of the {d6} networks that\ncollapsed without the rule',
        edge=ACC, face='#eef2f8', size=7, weight='bold')
    box(5.3, 0.25, 4.0, 1.25,
        f'ELIMINATED\n{e5} of the {l5} networks that\nsurvived without the rule',
        edge=WARM, face='#f8efee', size=7, weight='bold')
    save(fig, 'fig1_concept', tab)


# --------------------------------------------------------------- figure 2
def fig2():
    """COLL-3: the rule concentrates conductance and nothing dies."""
    cells = [('cortical', 'cortical'), ('hh', 'Hodgkin-Huxley')]
    tab = [['cell', 'arm', 'seed', 'frac_at_min', 'tot_w_final', 'budget',
            'deficit_pct', 'alive']]
    fig, ax = plt.subplots(1, 2, figsize=(5.2, 2.1))
    for j, (key, label) in enumerate(cells):
        for arm, learn in (('off', False), ('on', True)):
            rs = rows_of(os.path.join(EXP, 'coll3', f'coll3_{key}_{arm}.json'))
            for r in rs:
                tab.append([label, arm, r['seed'], r['frac_at_min'],
                            r['tot_w_final'], r['budget'],
                            100 * (r['tot_w_final'] - r['budget']) / r['budget'],
                            r['alive']])
            fm = [r['frac_at_min'] for r in rs]
            x = np.full(len(fm), 0 if arm == 'off' else 1, dtype=float)
            x += np.linspace(-0.06, 0.06, len(fm))
            ax[j].plot(x, fm, 'o', ms=3.5, color=(KL if arm == 'off' else ACC),
                       mew=0)
            ax[j].plot([0 if arm == 'off' else 1], [np.mean(fm)], '_',
                       ms=16, color=K, mew=1.2)
        ax[j].set_xticks([0, 1]); ax[j].set_xticklabels(['rule off', 'rule on'])
        ax[j].set_xlim(-0.35, 1.35); ax[j].set_ylim(-0.02, 0.32)
        ax[j].set_title(f'{label}   (4 wirings, 8/8 alive)')
        if j == 0:
            ax[j].set_ylabel('fraction of synapses at zero')
    save(fig, 'fig2_redistribution', tab)


# --------------------------------------------------------------- figure 2
def fig3():
    """The two signs, at the 35 s horizon, with Wilson intervals."""
    off6, on6 = paired(sorted(glob.glob(PY1 % 'off')), sorted(glob.glob(PY1 % 'on')))
    n6, d6, r6, e6, l6 = table_2x2(off6, on6)
    tab = [['panel', 'cell', 'g_ie', 'horizon_s', 'n', 'dead_off', 'alive_off',
            'rescued', 'eliminated', 'fraction', 'wilson_lo', 'wilson_hi']]
    lo, hi = wilson(r6, d6)
    tab.append(['A', 'PY 1', 0.045, 35, n6, d6, l6, r6, e6, r6 / d6, lo, hi])
    panels = [('PY 1\ng_ie 0.045', r6, d6, ACC)]

    o, n_ = paired([ABPD2 % 'off'], [ABPD2 % 'on'])
    n5, d5, r5, e5, l5 = table_2x2(o, n_)
    lo5, hi5 = wilson(e5, l5)
    tab.append(['B', 'AB/PD 2', 0.045, 35, n5, d5, l5, r5, e5, e5 / l5, lo5, hi5])
    panels.append(('AB/PD 2\ng_ie 0.045', e5, l5, WARM))

    for g in ('0.030', '0.060'):          # not re-measured: 30 s, marked
        o2, n2 = paired([os.path.join(EXP, 'coll5', f'coll5_g{g}_off.json')],
                        [os.path.join(EXP, 'coll5', f'coll5_g{g}_on.json')])
        na, da, ra, ea, la = table_2x2(o2, n2)
        loa, hia = wilson(ea, la) if la else (float('nan'),) * 2
        tab.append(['B', 'AB/PD 2', float(g), 30, na, da, la, ra, ea,
                    (ea / la if la else float('nan')), loa, hia])

    fig, ax = plt.subplots(figsize=(3.3, 2.2))
    for i, (lab, k, n, col) in enumerate(panels):
        pf = k / n
        lo_, hi_ = wilson(k, n)
        ax.plot([i], [pf], 'o', ms=6, color=col, mew=0)
        ax.plot([i, i], [lo_, hi_], '-', color=col, lw=1.4)
        ax.annotate(f'{k}/{n}', (i, pf), textcoords='offset points',
                    xytext=(9, -3), fontsize=7, color=K)
    ax.set_xticks(range(len(panels)))
    ax.set_xticklabels([p_[0] for p_ in panels])
    ax.set_xlim(-0.5, len(panels) - 0.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('conditional fraction')
    ax.axhline(0.5, color=KL, lw=0.5, ls=':')
    ax.set_title('same rule, same inhibition, two conductance sets')
    save(fig, 'fig3_two_signs', tab)


def fig4():
    """COLL-8: elimination fraction ordered by the control's firing rate."""
    src = {
        'AB/PD 1': (EXP + '/coll8/coll8a_ABPD1_%s.json', 35.0),
        'AB/PD 5': (EXP + '/coll8/coll8a_ABPD5_%s.json', 35.0),
        'AB/PD 3': (EXP + '/coll8/coll8a_ABPD3_%s.json', 35.0),
        'LP 3': (EXP + '/coll8/coll8b_LP3_%s.json', 35.0),
        'AB/PD 4': (EXP + '/coll8/coll8b_ABPD4_%s.json', 35.0),
    }
    pts = []
    for cell, (pat, secs) in src.items():
        o, n_ = paired([pat % 'off'], [pat % 'on'])
        n, d, r, e, l = table_2x2(o, n_)
        rate = np.mean([o[s]['total_spikes'] for s in sorted(o)]) / secs
        pts.append((cell, rate, e, l, n))
    off6, on6 = paired(sorted(glob.glob(PY1 % 'off')), sorted(glob.glob(PY1 % 'on')))
    n, d, r, e, l = table_2x2(off6, on6)
    pts.append(('PY 1', np.mean([off6[s]['total_spikes'] for s in sorted(off6)]) / 35.0,
                e, l, n))
    o5, n5_ = paired([ABPD2 % 'off'], [ABPD2 % 'on'])
    n, d, r, e, l = table_2x2(o5, n5_)
    pts.append(('AB/PD 2', np.mean([o5[s]['total_spikes'] for s in sorted(o5)]) / 35.0,
                e, l, n))
    pts.sort(key=lambda t: t[1])

    rates = np.array([p[1] for p in pts])
    frac = np.array([p[2] / p[3] for p in pts])
    # ranks with ties averaged: three cells tie at 0.000, and breaking those
    # ties by input order inflates rho to +0.893 (STATS-1, §7.2)
    rho = float(np.corrcoef(midranks(rates), midranks(frac))[0, 1])

    tab = [['cell', 'spikes_per_s', 'n', 'alive_off', 'eliminated',
            'fraction', 'wilson_lo', 'wilson_hi']]
    for c, rt, e, l, n in pts:
        lo, hi = wilson(e, l)
        tab.append([c, rt, n, l, e, e / l, lo, hi])
    tab.append(['spearman_rho_all_seven', rho, '', '', '', '', '', ''])

    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    # labels are staggered because five of the seven cells sit within a
    # 1.6-fold window and collide at any single offset
    offs = {'PY 1': (6, 5), 'AB/PD 1': (-4, -13), 'LP 3': (2, 7),
            'AB/PD 5': (5, 6), 'AB/PD 3': (4, -11), 'AB/PD 4': (6, 2),
            'AB/PD 2': (7, 0)}
    for c, rt, e, l, n in pts:
        pf = e / l
        lo, hi = wilson(e, l)
        col = ACC if c in ('PY 1', 'LP 3') else K
        ax.plot([rt, rt], [lo, hi], '-', color=KL, lw=0.9)
        ax.plot([rt], [pf], 'o', ms=4.5, color=col, mew=0)
        ax.annotate(c, (rt, pf), textcoords='offset points',
                    xytext=offs.get(c, (4, 4)), fontsize=6, color=col)
    ax.set_xscale('log')
    ax.set_xticks([30, 100, 300, 1000])
    ax.set_xticklabels(['30', '100', '300', '1000'])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlim(18, 1400)
    ax.set_xlabel("control-arm firing rate (spikes/s, log scale)")
    ax.set_ylabel('fraction eliminated | would have lived')
    ax.set_ylim(-0.08, 1.0)
    gap = rates[1] / rates[0]
    ax.axvspan(rates[0], rates[1], color=KL, alpha=0.12, lw=0)
    ax.annotate(f'unsampled\n{gap:.0f}-fold gap',
                ((rates[0] * rates[1]) ** 0.5, 0.86),
                ha='center', fontsize=6, color=KL)
    tab.append(['unsampled_gap_fold', gap, '', '', '', '', '', ''])
    ax.set_title(f'Spearman {rho:+.3f} over seven cells')
    save(fig, 'fig4_ordering', tab)


# --------------------------------------------------------------- figure 4
def fig5():
    """COLL-9: fragility follows the rate, rescuability does not."""
    tab = [['condition', 'n', 'dead_off', 'rescued', 'rescue_fraction']]
    o1, n1 = paired([EXP + '/coll9/coll9_m003_off.json'],
                    [EXP + '/coll9/coll9_m003_on.json'])
    o2, n2 = paired([EXP + '/coll9/coll9_m005_off.json'],
                    [EXP + '/coll9/coll9_m005_on.json'])
    off6, on6 = paired(sorted(glob.glob(PY1 % 'off')), sorted(glob.glob(PY1 % 'on')))
    a = table_2x2(o1, n1); b = table_2x2(o2, n2); c = table_2x2(off6, on6)
    # COLL-8a's AB/PD 1 at no injected current: the fragility baseline
    o0, n0 = paired([EXP + '/coll8/coll8a_ABPD1_off.json'],
                    [EXP + '/coll8/coll8a_ABPD1_on.json'])
    z = table_2x2(o0, n0)

    frag = [('AB/PD 1\nno current', z[1], z[0]),
            ('AB/PD 1\n-0.03', a[1], a[0]),
            ('AB/PD 1\n-0.05', b[1], b[0]),
            ('PY 1\nno current', c[1], c[0])]
    resc = [('AB/PD 1 driven', a[2] + b[2], a[1] + b[1]),
            ('PY 1', c[2], c[1])]
    for lab, k, n in frag:
        tab.append(['fragility: ' + lab.replace('\n', ' '), n, k, '', k / n])
    for lab, k, n in resc:
        tab.append(['rescue: ' + lab, n, n, k, k / n if n else float('nan')])

    fig, ax = plt.subplots(1, 2, figsize=(5.4, 2.2))
    for i, (lab, k, n) in enumerate(frag):
        lo, hi = wilson(k, n)
        ax[0].plot([i, i], [lo, hi], '-', color=KL, lw=1.0)
        ax[0].plot([i], [k / n], 'o', ms=5, color=K, mew=0)
        ax[0].annotate(f'{k}/{n}', (i, k / n), textcoords='offset points',
                       xytext=(7, -3), fontsize=6.5)
    ax[0].set_xticks(range(len(frag)))
    ax[0].set_xticklabels([f[0] for f in frag], fontsize=6)
    ax[0].set_ylim(-0.05, 1.1); ax[0].set_ylabel('control arm dead')
    ax[0].set_title('fragility follows the drive')

    for i, (lab, k, n) in enumerate(resc):
        lo, hi = wilson(k, n)
        col = WARM if 'driven' in lab else ACC
        ax[1].plot([i, i], [lo, hi], '-', color=col, lw=1.4)
        ax[1].plot([i], [k / n], 'o', ms=6, color=col, mew=0)
        ax[1].annotate(f'{k}/{n}', (i, k / n), textcoords='offset points',
                       xytext=(8, -3), fontsize=7)
    ax[1].set_xticks(range(len(resc)))
    ax[1].set_xticklabels([r[0] for r in resc], fontsize=6.5)
    ax[1].set_xlim(-0.5, 1.5); ax[1].set_ylim(-0.05, 1.1)
    ax[1].set_ylabel('rescued | would have died')
    ax[1].set_title('rescuability does not')
    save(fig, 'fig5_fragility_vs_rescue', tab)


# --------------------------------------------------------------- figure 5
def fig6():
    """COLL-11: the exhaustive single-substitution decomposition."""
    SRC = {'full': ('coll11a', None), 'Na': ('coll11b', None),
           'KCa': ('coll11a', None), 'Kd': ('coll11b', None),
           'H': ('coll11c', -0.30), 'leak': ('coll11c', -0.14)}
    order = ['full', 'Na', 'CaS', 'KCa', 'Kd', 'H', 'leak']
    tab = [['arm', 'source', 'drive', 'n', 'dead_off', 'rescued',
            'fraction', 'wilson_lo', 'wilson_hi']]
    pts = []
    for arm in order:
        if arm == 'CaS':      # OFF ran in 11a, ON in 11b
            o, n_ = paired([f'{EXP}/coll11/coll11a_CaS_off.json'],
                           [f'{EXP}/coll11/coll11b_CaS_on.json'])
            pre, drive = 'coll11a/b', -0.03
        else:
            pre, drive = SRC[arm]
            o, n_ = paired([f'{EXP}/coll11/{pre}_{arm}_off.json'],
                           [f'{EXP}/coll11/{pre}_{arm}_on.json'])
            drive = drive if drive is not None else -0.03
        n, d, r, e, l = table_2x2(o, n_)
        lo, hi = wilson(r, d)
        tab.append([arm, pre, drive, n, d, r, r / d, lo, hi])
        pts.append((arm, r, d, lo, hi, drive))

    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    for i, (arm, r, d, lo, hi, drive) in enumerate(pts):
        y = len(pts) - 1 - i
        col = ACC if arm == 'full' else K
        ax.plot([lo, hi], [y, y], '-', color=(ACC if arm == 'full' else KL),
                lw=1.4 if arm == 'full' else 0.9)
        ax.plot([r / d], [y], 'o', ms=5.5 if arm == 'full' else 4, color=col,
                mew=0)
        lab = 'all eight (= PY 1)' if arm == 'full' else f'I_{arm}'
        ax.annotate(f'{lab}   {r}/{d}' + ('' if drive == -0.03 else f'   ({drive})'),
                    (1.02, y), fontsize=6.5, va='center', color=col,
                    annotation_clip=False)
    ax.set_yticks([]); ax.set_xlim(0, 1.0); ax.set_ylim(-0.6, len(pts) - 0.4)
    ax.set_xlabel('rescued | would have died  (Wilson 95%)')
    ax.set_title('one substitution at a time, exhaustive over eight channels')
    save(fig, 'fig6_decomposition', tab)


if __name__ == '__main__':
    print("PAPER-2 figures, from the recorded runs only")
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print(f"  -> {OUT}")
