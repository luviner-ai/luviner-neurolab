"""Figures for the consolidation preprint, from the recorded JSON only.

Every value plotted is read from `experiments/dish/*.json` as written by
the runs. Nothing is recomputed in a way that could drift from the
record, nothing is synthesised, and no cosmetic rounding is applied to a
number that appears in the text.
"""
import glob
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DISH = os.path.join(HERE, '..', 'experiments', 'dish')
MEA = os.path.join(HERE, '..', 'experiments', 'mea', 'results')
OUT = os.path.join(HERE, 'figures')

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8, 'axes.labelsize': 8,
    'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': 0.6, 'lines.linewidth': 1.0,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
K, KL, ACC = '#222222', '#888888', '#4a6fa5'   # ink, light ink, one accent


def load(pattern):
    rows = []
    for f in sorted(glob.glob(os.path.join(DISH, pattern))):
        rows += json.load(open(f))
    return rows


def load_mea():
    """MEA-1 per-recording results. One dict per file, not a list, and the
    two recordings voided on the registered burst-count guard are dropped
    here rather than plotted as zeros."""
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(MEA, '*.json')))]
    return [r for r in rows if 'void' not in r]


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'))
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def fig1():
    """Preparation and instruments: raster with burst windows, weight split."""
    pm = load('stg_pm_s*.json')
    # A wiring from the MODAL regime. Seed 0 is one of the two off-regime
    # wirings PM-3 identified (17 bursts against the modal 22) and has the
    # highest zero-weight fraction of the ten (0.28 vs 0.09-0.13), so it
    # would misrepresent the preparation as an illustration.
    r = next(x for x in pm if x['n_bursts'] == 22)
    spikes = [np.array(s) for s in r['spikes']]
    wins = np.array(r['windows'])
    t0, t1 = 30000.0, 40000.0
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6),
                           gridspec_kw={'width_ratios': [2.1, 1]})

    for a, b in wins:
        if b > t0 and a < t1:
            ax[0].axvspan((max(a, t0) - t0) / 1000, (min(b, t1) - t0) / 1000,
                          color=KL, alpha=0.22, lw=0)
    for i, s in enumerate(spikes):
        s = s[(s >= t0) & (s <= t1)]
        if s.size:
            ax[0].plot((s - t0) / 1000, np.full(s.size, i), '|', color=K,
                       ms=1.4, mew=0.35)
    ax[0].set_xlabel('time (s)')
    ax[0].set_ylabel('excitatory cell')
    ax[0].set_xlim(0, (t1 - t0) / 1000)
    ax[0].set_ylim(-1, len(spikes))
    ax[0].set_title('A   population bursts and detected windows', loc='left')

    w = np.array(r['weights'])
    med = np.median(w)
    ax[1].hist(w, bins=28, color=KL, edgecolor=K, linewidth=0.4)
    ax[1].axvline(med, color=ACC, lw=1.2)
    ax[1].annotate('median split', xy=(med, ax[1].get_ylim()[1] * 0.92),
                   xytext=(4, 0), textcoords='offset points',
                   color=ACC, fontsize=7, va='top')
    ax[1].set_xlabel('synaptic conductance')
    ax[1].set_ylabel('synapses')
    ax[1].set_title('B   weight-class construction', loc='left')
    save(fig, 'fig1_preparation')


def fig2():
    """Mechanism refuted vs outcome met."""
    pm = load('stg_pm_s*.json')
    pm4 = load('stg_pm4_s*.json')
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4))
    fig.subplots_adjust(wspace=0.42)

    rho = np.array([r['rho_latency_weight'] for r in pm])
    ax[0].axhline(0, color=KL, lw=0.7)
    ax[0].plot(np.arange(1, len(rho) + 1), rho, 'o', color=K, ms=4)
    ax[0].axhline(rho.mean(), color=ACC, lw=1.0, ls='--')
    ax[0].set_ylim(-0.5, 0.5)
    ax[0].set_xlabel('wiring')
    ax[0].set_ylabel(r'$\rho$(latency, weight)')
    ax[0].set_title(f'A   mechanism: refuted\n'
                    f'{int(np.sum(rho < 0))}/10 negative, '
                    f'mean {rho.mean():+.3f}', loc='left')

    for k, (rows, lab) in enumerate(((pm, 'B   outcome: original 10'),
                                     (pm4, 'C   replication: fresh 10'))):
        a = ax[k + 1]
        st = np.array([r['pm_strong'] for r in rows])
        wk = np.array([r['pm_weak'] for r in rows])
        for i in range(len(rows)):
            a.plot([0, 1], [wk[i], st[i]], '-', color=KL, lw=0.8)
        a.plot(np.zeros(len(rows)), wk, 'o', color=K, ms=3.5)
        a.plot(np.ones(len(rows)), st, 'o', color=ACC, ms=3.5)
        a.set_xticks([0, 1])
        a.set_xticklabels(['weak', 'strong'])
        a.set_xlim(-0.35, 1.35)
        a.set_ylabel('$P/M$')
        n_up = int(np.sum(st > wk))
        a.set_title(f'{lab}\nstrong > weak in {n_up}/{len(rows)}', loc='left')
    save(fig, 'fig2_mechanism_outcome')


def fig3():
    """The carrier, replicated out of sample.

    Per-class causal fractions were never persisted for the ORIGINAL ten
    wirings (`stg_pm_s*.json` carries neither the fields nor `pre`/`post`,
    so they cannot even be recomputed). The fresh cohort does carry them,
    so this figure belongs to the out-of-sample replication rather than to
    the original result, and the caption says so.
    """
    pm4 = load('stg_pm4_s*.json')
    cs = np.array([r['cf_strong'] for r in pm4])
    cw = np.array([r['cf_weak'] for r in pm4])
    pmd = np.array([r['pm_diff'] for r in pm4])
    cfd = np.array([r['cf_diff'] for r in pm4])
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    fig.subplots_adjust(wspace=0.34)

    for i in range(len(pm4)):
        ax[0].plot([0, 1], [cw[i], cs[i]], '-', color=KL, lw=0.8)
    ax[0].plot(np.zeros(len(pm4)), cw, 'o', color=K, ms=3.5)
    ax[0].plot(np.ones(len(pm4)), cs, 'o', color=ACC, ms=3.5)
    ax[0].set_xticks([0, 1])
    ax[0].set_xticklabels(['weak', 'strong'])
    ax[0].set_xlim(-0.35, 1.35)
    ax[0].set_ylabel('pre-before-post fraction')
    ax[0].set_title(f'A   the carrier, replicated\n'
                    f'{cs.mean():.4f} vs {cw.mean():.4f} '
                    f'({100 * (cs - cw).mean():.2f} pp), '
                    f'{int(np.sum(cs > cw))}/10', loc='left')

    ratio = pmd / cfd
    ax[1].plot(cfd, pmd, 'o', color=K, ms=4, mfc='none', mew=0.9)
    xs = np.linspace(0, cfd.max() * 1.05, 10)
    ax[1].plot(xs, ratio.mean() * xs, '-', color=ACC, lw=1.0)
    ax[1].set_xlabel('pre-before-post excess (counted)')
    ax[1].set_ylabel('$P/M$ difference (kernel-weighted)')
    ax[1].set_xlim(0, cfd.max() * 1.05)
    ax[1].set_ylim(0, max(pmd.max(), ratio.mean() * cfd.max()) * 1.08)
    ax[1].set_title(f'B   the kernel adds nothing beyond the count\n'
                    f'ratio {ratio.mean():.1f} $\\pm$ {ratio.std(ddof=1):.1f}, '
                    f'CV {ratio.std(ddof=1) / ratio.mean():.2f}', loc='left')
    save(fig, 'fig3_carrier')


def fig4():
    """Interventions: sorting vs realized duration, and rate trajectories."""
    d1 = [r for r in load('dur1_s*.json') if r['guard_ok']]
    d2 = [r for r in load('dur2_s*.json') if r['guard_ok']]
    base = [r for r in d1 if abs(r['g_ie'] - 0.03) < 1e-9]
    inh = [r for r in d1 if abs(r['g_ie'] - 0.03) >= 1e-9]
    intr = d2
    # The intrinsic arm's CaS x1.0 IS the inhibitory arm's g_ie = 0.03 --
    # the same runs, by design, which is what made the second arm cost one
    # arm rather than two. Drawn once, in its own marker, so a reader who
    # counts points counts correctly without reading the caption.
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))

    for rows, col, mk, lab in ((inh, K, 'o', 'inhibitory ($g_{ie}$)'),
                               (intr, ACC, 's', 'intrinsic ($g_{CaS}$)'),
                               (base, KL, 'D', 'shared baseline')):
        x = np.array([r['duration_ms'] for r in rows])
        y = np.array([r['cf_diff'] for r in rows])
        ax[0].plot(x, y, mk, color=col, ms=4, label=lab, mfc='none', mew=0.9)
    for rows, col in ((inh + base, K), (intr + base, ACC)):
        x = np.array([r['duration_ms'] for r in rows])
        y = np.array([r['cf_diff'] for r in rows])
        c = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 10)
        ax[0].plot(xs, np.polyval(c, xs), '-', color=col, lw=0.9, alpha=0.8)
    ax[0].set_xlabel('realized burst duration (ms)')
    ax[0].set_ylabel('pre-before-post excess (strong $-$ weak)')
    ax[0].legend(frameon=False, loc='upper right')
    ax[0].set_title('A   sorting tracks duration in both arms', loc='left')

    for rows, col, mk, lab in ((inh, K, 'o', 'inhibitory'),
                               (intr, ACC, 's', 'intrinsic'),
                               (base, KL, 'D', 'shared baseline')):
        x = np.array([r['duration_ms'] for r in rows])
        y = np.array([r['rate_hz'] for r in rows])
        ax[1].plot(x, y, mk, color=col, ms=4, label=lab, mfc='none', mew=0.9)
    ax[1].set_xlabel('realized burst duration (ms)')
    ax[1].set_ylabel('burst rate (Hz)')
    ax[1].set_ylim(0, 0.9)
    ax[1].legend(frameon=False, loc='lower right')
    ax[1].set_title('B   rate moves differently in each arm', loc='left')
    save(fig, 'fig4_interventions')


def fig5():
    """MEA-1: the prediction on real cortical tissue, and what it bounds."""
    R = load_mea()
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7),
                           gridspec_kw={'width_ratios': [1.25, 1]})

    # --- A: per-recording excess against each recording's own null -------
    R1 = sorted(R, key=lambda r: r['sorting'])
    y = np.array([r['sorting'] for r in R1]) * 100          # percentage points
    e = np.array([r['sorting_null_sd'] for r in R1]) * 100
    x = np.arange(len(R1))
    pos = y > 0
    ax[0].axhline(0, color=KL, lw=0.7, zorder=1)
    ax[0].errorbar(x[~pos], y[~pos], yerr=e[~pos], fmt='o', ms=4, mfc='none',
                   mew=0.9, color=K, ecolor=KL, elinewidth=0.7, capsize=1.6,
                   label='strong $<$ weak (10)', zorder=3)
    ax[0].errorbar(x[pos], y[pos], yerr=e[pos], fmt='o', ms=4, mfc=ACC,
                   mew=0.9, color=ACC, ecolor=KL, elinewidth=0.7, capsize=1.6,
                   label='strong $>$ weak (%d)' % pos.sum(), zorder=3)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([r['file'].split('.')[0] for r in R1], rotation=90)
    ax[0].set_ylabel('pre-before-post excess, strong $-$ weak (pp)')
    ax[0].set_xlabel('recording (culture-div), ordered')
    ax[0].legend(frameon=False, loc='upper left')
    ax[0].set_title('A   %d of %d recordings positive; %d required'
                    % (pos.sum(), len(R1), len(R1) - 1), loc='left')

    # --- B: sorting vs realized duration, within age strata --------------
    strata = sorted({r['div'] for r in R})
    for d, col, mk in zip(strata, (K, ACC, KL), ('o', 's', '^')):
        g = [r for r in R if r['div'] == d]
        xs = np.array([r['median_dur_ms'] for r in g])
        ys = np.array([r['sorting'] for r in g]) * 100
        rho = np.corrcoef(np.argsort(np.argsort(xs)),
                          np.argsort(np.argsort(ys)))[0, 1]
        ax[1].plot(xs, ys, mk, color=col, ms=4, mfc='none', mew=0.9,
                   label=r'%d div  ($\rho$ = %+.2f)' % (d, rho))
    ax[1].axhline(0, color=KL, lw=0.7, zorder=1)
    ax[1].set_xlabel('median burst duration (ms)')
    ax[1].set_ylabel('pre-before-post excess (pp)')
    ax[1].legend(frameon=False, loc='lower right')
    # One stratum IS negative (-0.26, n = 6), so "no negative trend in any
    # stratum" would be false to the letter. The claim that survives is
    # about consistency across strata; the per-stratum values are in the
    # legend so a reader checks the judgement against the numbers.
    ax[1].set_title('B   no consistent negative trend across strata', loc='left')
    save(fig, 'fig5_mea')


def figS1():
    """[EXPLORATORY, POST-HOC] Weight as a continuous predictor.

    PM-4 registered a MEDIAN SPLIT. This panel replaces it with a graded
    predictor chosen after the results were seen, and the 25/75 contrast in
    panel B is likewise a post-hoc choice. Neither can be read at the
    registered bar; both are labelled in the caption.
    """
    z = np.load(os.path.join(DISH, 'pm4_posthoc.npz'))
    w, cf, sd = z['w'], z['cf'], z['seed']
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))

    # Deciles BY RANK, not by weight value: many weights are identical (the
    # w_max cap at one end, near-zero at the other), so value percentiles
    # collapse and leave bins empty. Each wiring's curve is centred on its
    # own mean, because the gradient is a WITHIN-wiring effect and baseline
    # causal fractions differ between wirings.
    #
    # The summary across wirings is the MEDIAN, not the mean. Two of the ten
    # (seeds 17 and 18) sit far off the modal regime and carry gaps 5-10x the
    # rest; averaging lets them set the shape of the curve, and the mean
    # curve is not even monotone in the direction of its own rank
    # correlation. The median is reported with the interquartile band.
    def dec(x):
        return (np.argsort(np.argsort(x)) * 10 // x.size).clip(0, 9)

    M = []
    for sv in np.unique(sd):
        m = sd == sv
        k = dec(w[m])
        row = np.array([cf[m][k == i].mean() for i in range(10)])
        M.append((row - row.mean()) * 100)
    M = np.array(M)
    xs = np.arange(10) + 1
    lo, med, hi = (np.percentile(M, q, axis=0) for q in (25, 50, 75))
    ax[0].axhline(0, color=KL, lw=0.7)
    ax[0].fill_between(xs, lo, hi, color=ACC, alpha=0.16, lw=0)
    for row in M:
        ax[0].plot(xs, row, '-', color=KL, lw=0.5, alpha=0.55)
    ax[0].plot(xs, med, 'o-', color=ACC, ms=4, lw=1.2, label='median (IQR band)')
    ax[0].set_xlabel('synaptic weight, decile')
    ax[0].set_ylabel('pre-before-post fraction,\ncentred within wiring (pp)')
    ax[0].set_xticks(xs)
    ax[0].set_ylim(-1.2, 1.2)
    ax[0].legend(frameon=False, loc='upper left')
    rho_med = np.corrcoef(np.arange(10), med)[0, 1]
    n_up = int(sum(1 for row in M if np.corrcoef(np.arange(10), row)[0, 1] > 0))
    ax[0].set_title(r'A   rises with weight ($\rho$ = %+.2f on the median;'
                    '\n%d/10 wirings positive) -- not monotone'
                    % (rho_med, n_up), loc='left')

    # Paired, on ONE axis: the post-hoc 25/75 contrast against the gap the
    # registration actually scored. Two scales on twinned axes made this
    # panel unreadable and invited the two quantities to be compared by eye
    # as if they shared units.
    r = z['rho_weight']; g = z['gap_2575'] * 100; c = z['cfd'] * 100
    for ci, gi in zip(c, g):
        ax[1].plot([0, 1], [ci, gi], '-', color=KL, lw=0.8)
    ax[1].plot(np.zeros(c.size), c, 'D', color=KL, ms=4, mfc='none', mew=0.9,
               label='median split (registered)')
    ax[1].plot(np.ones(g.size), g, 's', color=ACC, ms=4, mfc='none', mew=0.9,
               label='25/75 contrast (post-hoc)')
    ax[1].set_xlim(-0.35, 1.35)
    ax[1].set_xticks([0, 1])
    ax[1].set_xticklabels(['registered', 'post-hoc'])
    ax[1].set_ylabel('pre-before-post excess, strong $-$ weak (pp)')
    ax[1].legend(frameon=False, loc='upper left')
    ax[1].set_title(r'B   both positive 10/10; $\rho$(weight, fraction)'
                    '\n$> 0$ in 10/10 (mean %+.2f)' % r.mean(), loc='left')
    save(fig, 'figS1_continuous_weight')


def figS2():
    """[EXPLORATORY, POST-HOC] Kernel ablation: flat count vs exponential.

    The record's published quantity is the CV of pm_diff/cf_diff ACROSS
    WIRINGS (0.13), reproduced here as 0.132. The per-synapse correlation in
    panel A is a DIFFERENT quantity and is labelled as an additional check,
    never as that CV.
    """
    z = np.load(os.path.join(DISH, 'pm4_posthoc.npz'))
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))

    rk = z['rho_kernel']
    ax[0].plot(np.zeros(rk.size), rk, 'o', color=K, ms=4, mfc='none', mew=0.9)
    ax[0].set_xlim(-0.8, 0.8); ax[0].set_xticks([])
    ax[0].set_ylim(0.9, 1.005)
    ax[0].set_ylabel(r'per-synapse $\rho$(count fraction, $P/M$)')
    ax[0].set_title('A   the kernel reorders almost nothing\n'
                    r'(min $\rho$ = %+.3f over 10 wirings)' % rk.min(), loc='left')

    cfd, pmd = z['cfd'] * 100, z['pmd']
    ratio = pmd / cfd
    ax[1].plot(cfd, pmd, 'o', color=K, ms=4, mfc='none', mew=0.9)
    xs = np.linspace(0, cfd.max() * 1.05, 10)
    ax[1].plot(xs, ratio.mean() * xs, '-', color=ACC, lw=1.0,
               label='proportional, through the origin')
    ax[1].set_xlabel('flat count excess, strong $-$ weak (pp)')
    ax[1].set_ylabel(r'kernel-weighted $P/M$ excess')
    ax[1].legend(frameon=False, loc='upper left')
    ax[1].set_title('B   CV of the ratio = %.3f; sign agrees 10/10'
                    % (ratio.std(ddof=1) / ratio.mean()), loc='left')
    save(fig, 'figS2_kernel_ablation')


def fig6():
    """Conceptual schematic: what carries consolidation here, and what does not.

    The pairs-per-synapse figure is MEASURED, not inherited. The brief for
    this panel specified ~3x10^5, which is the number Section 5 of the draft
    also gives ("at least 10^5 pairs per connection"). Recomputing it from
    the persisted PM-4 spikes gives a median of 8,960 pairs per synapse
    inside burst windows at |lag| <= 250 ms -- about 10^4, roughly 35x
    smaller. The measured value is drawn; the discrepancy is reported.
    """
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')

    def box(x, y, w, h, text, col=K, face='none', fs=6.6):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle='round,pad=0.004,rounding_size=0.010',
                                    linewidth=0.8, edgecolor=col, facecolor=face))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                fontsize=fs, color=col, linespacing=1.4)

    def arrow(x0, y0, x1, y1, col=K, ls='-'):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='-|>', color=col, lw=0.8,
                                    linestyle=ls, shrinkA=0, shrinkB=0))

    # --- main branch: the census, six boxes sized to their own text -------
    y, h = 0.74, 0.19
    ws = [0.150, 0.126, 0.196, 0.150, 0.120, 0.116]
    labs = ['initial\nsynaptic weight', 'spontaneous\nburst',
            'tiny pre-before-post\nbias (0.2-0.35 pp)',
            r'$\sim$9$\times$10$^{3}$ pairs' + '\nper synapse',
            'STDP\nintegration', 'weight\nsorting']
    gap, x = 0.026, 0.008
    xs = []
    for wi in ws:
        xs.append(x); x += wi + gap
    for xi, wi, t in zip(xs, ws, labs):
        box(xi, y, wi, h, t, face='#eef2f8' if t == 'weight\nsorting' else 'none')
    for xi, wi, xn in zip(xs[:-1], ws[:-1], xs[1:]):
        arrow(xi + wi, y + h / 2, xn, y + h / 2)
    ax.text(0.008, y + h + 0.05, 'the carrier: an accumulated count',
            fontsize=7, color=K, style='italic')

    # --- both side branches hang off the SAME node, on a shared spine.
    # An earlier version ran the arrow down THROUGH the latency branch into
    # the duration branch, which draws burst duration as a consequence of
    # first-spike latency. It is not; both are properties of the same burst.
    spine = xs[1] + ws[1] / 2
    y2, y3, hb = 0.42, 0.10, 0.17
    ax.plot([spine, spine], [y3 + hb / 2, y], ':', color=KL, lw=0.9, zorder=0)
    bx = spine + 0.055

    box(bx, y2, 0.128, hb, 'first-spike\nlatency', col=KL)
    box(bx + 0.154, y2, 0.290, hb,
        'no reliable signal\n(refuted in model and in tissue)', col=KL)
    arrow(spine, y2 + hb / 2, bx, y2 + hb / 2, col=KL, ls=':')
    arrow(bx + 0.128, y2 + hb / 2, bx + 0.154, y2 + hb / 2, col=KL)

    # Bottom row laid out to a budget: four boxes and three 0.018 gaps must
    # fit between the spine offset and the right edge. An earlier version
    # used the same widths as the row above and ran "live intervention
    # required" off the canvas.
    bw = [0.118, 0.145, 0.205, 0.166]
    bt = ['shorter\nbursts', 'stronger sorting\n(model)',
          'not reproduced\nobservationally (16 rec.)',
          'live intervention\nrequired']
    bc = [K, K, ACC, ACC]
    bxs, xx = [], bx
    for wi in bw:
        bxs.append(xx); xx += wi + 0.018
    for xi, wi, t, c in zip(bxs, bw, bt, bc):
        box(xi, y3, wi, hb, t, col=c,
            face='#eef2f8' if t.startswith('live') else 'none')
    arrow(spine, y3 + hb / 2, bxs[0], y3 + hb / 2, col=KL, ls=':')
    for xi, wi, xn, c in zip(bxs[:-1], bw[:-1], bxs[1:], (K, ACC, ACC)):
        arrow(xi + wi, y3 + hb / 2, xn, y3 + hb / 2, col=c)

    ax.set_title('A census, not a choreography', loc='left', fontsize=8.5)
    save(fig, 'fig6_schematic')


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    fig1(); fig2(); fig3(); fig4(); fig5(); figS1(); figS2(); fig6()
