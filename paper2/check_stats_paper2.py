"""Every statistic PAPER-2 quotes, recomputed from the run outputs.

`check_figures_paper2.py` joins the figures to the text. It cannot catch a
wrong p-value, because no p is plotted: two p quoted from the one-sided tail
survived five review rounds, two outside models and a referee, and were found
by recomputation rather than by rereading. This script is the recomputation,
kept so that a reader can run it instead of trusting it.

Counts, conditional probabilities, Wilson intervals and every test are derived
here from `experiments/coll*/**.json` with this file's own loader. It never
reads the figure CSVs: those are written by `make_figures_paper2.py`, and
checking the text against them would check one derivation against itself.

Tests are exact and implemented here rather than imported, so the script runs
on the deposit's declared dependencies (numpy and matplotlib) and in fact needs
neither: Fisher by enumeration of the hypergeometric over both tails, binomial
by summation, Spearman by full enumeration of the permutation null, McNemar as
the exact binomial on discordant pairs.

A mismatch is an error, not a warning, and the rule registered before this
script existed is that a failure corrects the manuscript and never the
calculation. What the script cannot reach is listed by name at the end: a
checker silent about its blind spots manufactures the assurance that let those
two p stand.
"""
import glob
import json
import os
import re
import sys
from itertools import permutations
from math import comb, floor, log10, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, '..', 'experiments')
DOC = os.path.join(HERE, 'preprint2.md')

HORIZON = 35.0          # seconds; the horizon of every run quoted in §2.2-§2.6
PY1 = os.path.join(EXP, 'coll14', 'coll14_PY1_*_%s.json')
ABPD2 = os.path.join(EXP, 'coll14', 'coll14_ABPD2_30_%s.json')


# --------------------------------------------------------------- exact tests
def hypergeom(a, b, c, d):
    """P(X = x) over the 2x2 tables with the observed margins."""
    n, r1, c1 = a + b + c + d, a + b, a + c
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return {x: comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
            for x in range(lo, hi + 1)}


def fisher(a, b, c, d):
    """Two-sided Fisher exact: every table no more probable than the observed."""
    p = hypergeom(a, b, c, d)
    return min(1.0, sum(v for v in p.values() if v <= p[a] * (1 + 1e-9)))


def fisher_one_sided(a, b, c, d):
    return sum(v for x, v in hypergeom(a, b, c, d).items() if x >= a)


def binom_ge(k, n, p=0.5):
    """Exact binomial, P(X >= k) — the test against p <= 0.5 in the Methods."""
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def mcnemar(b, c):
    """Exact McNemar: two-sided binomial on the discordant pairs."""
    n, k = b + c, min(b, c)
    return min(1.0, 2 * sum(comb(n, i) * 0.5 ** n for i in range(k + 1)))


def wilson(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    p, d = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def midranks(v):
    """Ranks with ties averaged — the definition Spearman's rho is built on."""
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


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / sqrt(sxx * syy)


def spearman(x, y):
    return pearson(midranks(x), midranks(y))


def spearman_perm(x, y):
    """rho and its permutation null by full enumeration (n! orderings of y)."""
    rho, rx = spearman(x, y), midranks(x)
    ge = two = total = 0
    for perm in permutations(y):
        r = pearson(rx, midranks(list(perm)))
        total += 1
        ge += r >= rho - 1e-12
        two += abs(r) >= abs(rho) - 1e-12
    return rho, ge / total, two / total, total


# ------------------------------------------------------------- loading runs
def rows_of(path):
    d = json.load(open(path))
    return d['rows'] if isinstance(d, dict) and 'rows' in d else d


def paired(off_glob, on_glob):
    """Per-seed control and learning rows, keyed by seed."""
    off, on = {}, {}
    for pat, into in ((off_glob, off), (on_glob, on)):
        paths = sorted(glob.glob(pat)) if '*' in pat else [pat]
        if not paths:
            raise SystemExit(f"check_stats: no run files matched {pat}")
        for p in paths:
            for r in rows_of(p):
                into[r['seed']] = r
    return off, on


def table(off, on):
    """(n, dead_off, rescued, alive_off, eliminated) over the shared seeds."""
    seeds = sorted(set(off) & set(on))
    dead = [s for s in seeds if not off[s]['alive']]
    live = [s for s in seeds if off[s]['alive']]
    return (len(seeds), len(dead), sum(1 for s in dead if on[s]['alive']),
            len(live), sum(1 for s in live if not on[s]['alive']))


def spikes_per_s(off):
    return sum(off[s]['total_spikes'] for s in off) / len(off) / HORIZON


# ------------------------------------------------------------- the assertions
doc = open(DOC, encoding='utf-8').read()
# The manuscript is hard-wrapped, so every pattern is matched against a copy
# with runs of whitespace flattened to one space: a check that breaks when a
# sentence rewraps is testing the wrapping, not the number.
flat = re.sub(r'\s+', ' ', doc)
fails, uncovered = [], []
checks = 0

# the manuscript spells small numbers out in prose and prints them in tables
WORD = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five',
        6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten',
        11: 'eleven', 12: 'twelve'}


def cnt(k):
    """A count as either digits or the word the manuscript may spell it as."""
    return f'(?:{k}|{WORD[k]})' if k in WORD else str(k)


def says(pattern, label):
    """The manuscript must contain `pattern`, which carries a computed value."""
    global checks
    checks += 1
    if not re.search(pattern, flat):
        fails.append(f"{label}: the manuscript does not contain /{pattern}/")


def sci(p, digits=1):
    """A p-value formatted the way the manuscript writes it."""
    e = floor(log10(p))
    m = round(p / 10 ** e, digits)
    if m >= 10:
        m, e = m / 10, e + 1
    return rf"{m:.{digits}f} × 10<sup>−{-e}</sup>"


def esc(s):
    return re.escape(s)


# ---- 2.2, PY 1: the rescue conditional ----------------------------------
off, on = paired(PY1 % 'off', PY1 % 'on')
n, dead, rescued, live, killed = table(off, on)
on_dead_of_dead = dead - rescued
lo, hi = wilson(rescued, dead)
p_binom = binom_ge(rescued, dead)

says(rf'collapses in {dead} of {n}\b', '2.2 PY 1 control mortality')
says(rf'\|\s*\*\*rule off, dead\*\*\s*\|\s*\*\*{rescued}\*\*\s*\|\s*{on_dead_of_dead}\s*\|',
     '2.2 PY 1 2x2 bottom row')
says(rf'\|\s*\*\*rule off, alive\*\*\s*\|\s*{live}\s*\|\s*{killed}\s*\|',
     '2.2 PY 1 2x2 top row')
says(rf'{rescued}/{dead} = {rescued / dead:.3f}', '2.2 PY 1 conditional')
says(rf'\[{lo:.3f}, {hi:.3f}\]', '2.2 PY 1 Wilson interval')
says(esc(sci(p_binom)), '2.2 PY 1 exact binomial against p <= 0.5')
says(r'the rule eliminates none' if killed == 0 else r'\Znever',
     '2.2 PY 1 eliminates none of the survivors')
says(rf'{cnt(live)} wirings that would have survived',
     '2.2 PY 1 number of survivors in prose')
agree = all(r['alive'] == r['alive_strict'] for r in list(off.values()) + list(on.values()))
says(r'the two detectors agree on every wiring' if agree else r'\Znever',
     '2.2 PY 1 the two detectors agree')

# ---- 2.2, AB/PD 2: the elimination conditional ---------------------------
o2, n2 = paired(ABPD2 % 'off', ABPD2 % 'on')
_, dead2, rescued2, live2, killed2 = table(o2, n2)
lo2, hi2 = wilson(killed2, live2)

says(rf'survives in {live2} of {live2 + dead2}\b', '2.2 AB/PD 2 control survival')
says(rf'eliminates {killed2} of those {live2}\*\*', '2.2 AB/PD 2 eliminated')
says(rf'\({killed2 / live2:.3f}, Wilson', '2.2 AB/PD 2 elimination fraction')
says(rf'Wilson \[{lo2:.3f}, {hi2:.3f}\]', '2.2 AB/PD 2 Wilson interval')
says(r'rescuing the single wiring that died' if (dead2, rescued2) == (1, 1) else r'\Znever',
     '2.2 AB/PD 2 the one control death is rescued')
agree2 = all(r['alive'] == r['alive_strict'] for r in list(o2.values()) + list(n2.values()))
says(r'On both cells, at 35 s, the two' if agree2 else r'\Znever',
     '2.2 AB/PD 2 the two detectors agree')

# ---- 2.2, the other two inhibition levels (COLL-5) -----------------------
C5 = os.path.join(EXP, 'coll5', 'coll5_g%s_%s.json')
o6, n6 = paired(C5 % ('0.060', 'off'), C5 % ('0.060', 'on'))
_, _, _, live6, killed6 = table(o6, n6)
says(rf'eliminates {killed6} of {live6}\b', '2.2 AB/PD 2 at g_ie 0.060')
o3, n3 = paired(C5 % ('0.030', 'off'), C5 % ('0.030', 'on'))
_, dead3, _, _, killed3 = table(o3, n3)
says(r'nothing dies on either arm' if (dead3, killed3) == (0, 0) else r'\Znever',
     '2.2 nothing dies at g_ie 0.030')

# ---- 2.3: the seven-cell ordering ---------------------------------------
C8 = os.path.join(EXP, 'coll8', 'coll8%s_%s_%%s.json')
SEVEN = [('AB/PD 1', C8 % ('a', 'ABPD1')), ('AB/PD 5', C8 % ('a', 'ABPD5')),
         ('AB/PD 3', C8 % ('a', 'ABPD3')), ('LP 3', C8 % ('b', 'LP3')),
         ('AB/PD 4', C8 % ('b', 'ABPD4'))]
cells = []
for name, pat in SEVEN:
    o, n_ = paired(pat % 'off', pat % 'on')
    tot, _, _, alive, elim = table(o, n_)
    cells.append((name, spikes_per_s(o), tot, alive, elim))
cells.append(('PY 1', spikes_per_s(off), n, live, killed))
cells.append(('AB/PD 2', spikes_per_s(o2), live2 + dead2, live2, killed2))
cells.sort(key=lambda c: c[1])

for name, rate, tot, alive, elim in cells:
    row = rf'\| {esc(name)} \| {rate:.1f} \| {tot} \| {elim / alive:.3f} \|'
    says(row, f'2.3 table row for {name}')

rates = [c[1] for c in cells]
fracs = [c[4] / c[3] for c in cells]
rho7, one7, two7, n_perm7 = spearman_perm(rates, fracs)
says(rf'Spearman \+{rho7:.3f}, p = {two7:.4f}\*\*', '2.3 Spearman over seven cells')
says(rf'\({n_perm7:,} *\n?permutations at seven cells\)|{n_perm7:,}\s*\n?permutations',
     '2.3 the permutation count in the Methods')

six = [c for c in cells if c[0] != 'PY 1']
rho6, one6, two6, _ = spearman_perm([c[1] for c in six], [c[4] / c[3] for c in six])
says(rf'Spearman \+{rho6:.3f} at p = {two6:.3f} over six cells',
     '2.3 Spearman excluding PY 1')

p_mcn = mcnemar(killed2, rescued2)
says(rf'{cnt(killed2)} eliminated against {cnt(rescued2)} rescued is {cnt(killed2 + rescued2)} discordant wirings, p = {p_mcn:.3f}',
     '2.3 McNemar on AB/PD 2')

# ---- 2.4 and 2.5: the within-cell drive manipulation ---------------------
C9 = os.path.join(EXP, 'coll9', 'coll9_m%s_%s.json')
nat_o, nat_n = paired(C8 % ('a', 'ABPD1') % 'off', C8 % ('a', 'ABPD1') % 'on')
_, nat_dead, _, nat_live, nat_elim = table(nat_o, nat_n)
d3_o, d3_n = paired(C9 % ('003', 'off'), C9 % ('003', 'on'))
_, dead_003, resc_003, live_003, elim_003 = table(d3_o, d3_n)
d5_o, d5_n = paired(C9 % ('005', 'off'), C9 % ('005', 'on'))
_, dead_005, resc_005, live_005, elim_005 = table(d5_o, d5_n)

says(rf'\*\*{nat_elim / nat_live:.3f} at its natural rate of {spikes_per_s(nat_o):.0f} spikes/s to {elim_003 / live_003:.3f} at the driven setting\*\*',
     '2.4 elimination moves from 0.000 to 1.000')
says(rf'eliminated all {cnt(live_003)} wirings that would otherwise have survived'
     if elim_003 == live_003 else r'\Znever', '2.4 all survivors eliminated')

says(rf'no injected current \| ~{spikes_per_s(nat_o):.0f} spikes/s \| {nat_dead}/{nat_dead + nat_live} \|',
     '2.5 table, AB/PD 1 at its natural drive')
says(rf'−0\.03 µA/cm² \| — \| {dead_003}/{dead_003 + live_003} \| {resc_003}/{dead_003} \|',
     '2.5 table, AB/PD 1 at -0.03')
says(rf'\| {dead_005}/{dead_005 + live_005} \| {resc_005}/{dead_005} \|',
     '2.5 table, AB/PD 1 at -0.05')
says(rf'\*\*{rescued}/{dead}\*\* \|', '2.5 table, PY 1')
says(rf'{nat_dead} of {nat_dead + nat_live} controls collapse at baseline and {dead_005} of {dead_005 + live_005} at the driven setting',
     '2.5 fragility follows the rate')

early = json.load(open(os.path.join(EXP, 'coll9', 'coll9_early.json')))
e5 = early['-0.05']
says(rf'The {sum(e5.values()) / len(e5):.0f} spikes/s above is measured over the first 5 s',
     '2.5 the early rate of the driven arm')
whole = sum(r['rate_hz'] for r in d5_o.values()) / len(d5_o)
says(rf'averages \*\*{whole:.0f} spikes/s and produces zero network bursts\*\*'
     if all(r['n_bursts'] == 0 for r in d5_o.values()) else r'\Znever',
     '2.5 the whole-run rate and the absence of bursts')

driven_resc = resc_003 + resc_005
driven_opp = dead_003 + dead_005
p_matched = fisher(rescued, dead - rescued, driven_resc, driven_opp - driven_resc)
says(rf'\*\*{driven_resc} of {driven_opp}\*\* \(\*\*p = {esc(sci(p_matched))}',
     '2.5 Fisher, PY 1 against the driven cell at matched fragility')

# ---- 2.5: the non-circular early-rate control ---------------------------
e3 = early['-0.03']
surv = {str(s) for s in d3_o if d3_o[s]['alive']}
ranked = sorted(e3, key=lambda s: -e3[s])
separated = set(ranked[:len(surv)]) == surv
p_sep = 1 / comb(len(e3), len(surv))
says(rf'the {cnt(len(surv))} wirings that survived are exactly the {cnt(len(surv))} with the highest early rates'
     if separated else r'\Znever', '2.5 complete rank separation')
says(rf'1/C\({len(e3)},{len(surv)}\) = {p_sep:.4f}', '2.5 the rank-separation p')

# ---- 2.6: the substitution grid -----------------------------------------
C11 = os.path.join(EXP, 'coll11', '%s_%s.json')
ARMS = [('all eight (= PY 1)', C11 % ('coll11a_full', 'off'), C11 % ('coll11a_full', 'on')),
        ('I_Na, fast sodium', C11 % ('coll11b_Na', 'off'), C11 % ('coll11b_Na', 'on')),
        ('I_CaS, slow calcium', C11 % ('coll11a_CaS', 'off'), C11 % ('coll11b_CaS', 'on')),
        ('I_KCa, calcium-activated potassium', C11 % ('coll11a_KCa', 'off'), C11 % ('coll11a_KCa', 'on')),
        ('I_Kd, delayed-rectifier potassium', C11 % ('coll11b_Kd', 'off'), C11 % ('coll11b_Kd', 'on')),
        ('I_H, hyperpolarisation-activated inward', C11 % ('coll11c_H', 'off'), C11 % ('coll11c_H', 'on')),
        ('I_leak, leak', C11 % ('coll11c_leak', 'off'), C11 % ('coll11c_leak', 'on'))]
subs = []
for label, off_p, on_p in ARMS:
    o, n_ = paired(off_p, on_p)
    _, d, r, _, _ = table(o, n_)
    rate = sum(x['rate_hz'] for x in o.values()) / len(o)
    p = fisher(r, d - r, driven_resc, driven_opp - driven_resc)
    subs.append((label, d, r, rate, p))
    B = r'\*?\*?'       # the complete substitution's row is bold throughout
    says(rf'\| {B}{esc(label)}{B} \| −?\d[\d.]* \| {B}{d}/7{B} \| {B}{r}/{d}{B} \|',
         f'2.6 table counts for {label}')
    cell = rf'{B}< 10<sup>−4</sup>{B}' if p < 1e-4 else rf'{p:.3f}'
    says(rf'\| {B}{esc(label)}{B} \|[^|]*\|[^|]*\|[^|]*\| {cell} \|',
         f'2.6 Fisher against 0/21 for {label}')

pooled_r = sum(s[2] for s in subs[1:])
pooled_d = sum(s[1] for s in subs[1:])
full_r, full_d = subs[0][2], subs[0][1]
p_pool = fisher(full_r, full_d - full_r, pooled_r, pooled_d - pooled_r)
says(rf'pooled rescue {pooled_r} of {pooled_d}; the complete substitution rescues {full_r} of {full_d}; p = {esc(sci(p_pool))}',
     '2.6 the pooled comparison')
says(rf'control rate is {subs[0][3]:.0f} spikes/s', '2.6 the complete arm rate')

matched = [s for s in subs[1:] if abs(s[3] - subs[0][3]) < 10]
says(rf'I_Na at {matched[0][3]:.0f} and I_KCa at {matched[1][3]:.0f} spikes/s',
     '2.6 the two rate-matched arms')
p_na = fisher(full_r, full_d - full_r, matched[0][2], matched[0][1] - matched[0][2])
p_kca = fisher(full_r, full_d - full_r, matched[1][2], matched[1][1] - matched[1][2])
says(rf'\*\*p = {p_na:.4f} and p = {p_kca:.3f}\*\*', '2.6 the two direct comparisons')
unmatched = [s[3] for s in subs[1:] if s not in matched]
says(r'\(' + ', '.join(f'{r:.0f}' for r in unmatched[:-1]) + rf' and {unmatched[-1]:.0f} spikes/s\)',
     '2.6 the four arms that are not rate-matched')

# ---- 4: the limitation that carries PY 1's own elimination estimate ------
lo_py, hi_py = wilson(killed, live)
says(rf'elimination estimate is \*\*{killed} of {live}\*\*, Wilson \[{lo_py:.3f}, {hi_py:.3f}\]',
     '4 PY 1 elimination estimate and its interval')

# ---- 7.2: the withdrawn arm that the minimum-denominator rule stopped ----
lo_w, on_w = paired(C11 % ('coll11b_leak', 'off'), C11 % ('coll11b_leak', 'on'))
_, d_w, r_w, _, _ = table(lo_w, on_w)
p_w = fisher_one_sided(r_w, d_w - r_w, driven_resc, driven_opp - driven_resc)
says(rf'scoring {r_w} of {d_w} at p = {p_w:.4f}', '7.2 the withdrawn underpowered arm')
says(rf'that arm is {subs[6][2]} of {subs[6][1]}\b' if subs[6][0].startswith('I_leak') else r'\Znever',
     '7.2 the same arm at a proper denominator')

# ---- 6: the robustness grid ---------------------------------------------
C12 = os.path.join(EXP, 'coll12', 'coll12_%s_n%d_%s_%s.json')
for size in (48, 24):
    for fam in ('A', 'B'):
        counts = {}
        for cell in ('PY1', 'ABPD1'):
            o, n_ = paired(C12 % (cell, size, fam, 'off'), C12 % (cell, size, fam, 'on'))
            _, d, r, _, _ = table(o, n_)
            counts[cell] = (d, r, sum(x['rate_hz'] for x in o.values()) / len(o))
        (dp, rp, ratep), (da, ra, _) = counts['PY1'], counts['ABPD1']
        label = rf'family {fam}, {size} excitatory cells'
        if dp and da:
            p12 = fisher(rp, dp - rp, ra, da - ra)
            says(rf'\| {label} \| {rp}/{dp} = {rp / dp:.3f} \| {ra}/{da} \| \*\*no\*\*, p = {esc(sci(p12))} \|',
                 f'6 row for {label}')
        else:
            says(rf'\| {label} \| — \| {ra}/{da} \| \*\*not testable\*\* \|',
                 f'6 row for {label}')
        says(rf'{ratep:.0f}', f'6 control-arm rate for family {fam} at {size} cells')
says(r'does not collapse at all, in either initialisation family \(0 of 11\)',
     '6 the size-24 control arm does not collapse')

# --------------------------------------------------------------- blind spots
uncovered += [
    "the acceptance gates (batched vs unbatched, identity swap): they compare "
    "spike times, not statistics, and `coll5/gate.py` and `coll6/gate.py` run "
    "them",
    "the COLL-3 redistribution percentages in §2.1: checked against the "
    "figures by check_figures_paper2.py, which reads the same runs",
    "the execution record and lot timings in Supplementary Note S1: process "
    "times, not statistics",
    "every count in §2.2-§2.6 whose denominator is a design parameter (seeds "
    "requested, horizon, quorum): those are inputs, not results",
]

print(f"{checks} checks against {os.path.basename(DOC)}")
print(f"not covered by this script ({len(uncovered)}):")
for u in uncovered:
    print(f"  - {u}")
if fails:
    print(f"\n{len(fails)} MISMATCH(ES):")
    for f in fails:
        print(f"  - {f}")
    sys.exit(1)
print("\nevery statistic in the manuscript matches the deposited runs")
