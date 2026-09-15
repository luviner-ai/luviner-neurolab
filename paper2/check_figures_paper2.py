"""Every number a PAPER-2 figure plots, checked against the manuscript text.

The figures are generated from the run JSON; the manuscript is written by
hand. This script is the join between them. It parses the numbers out of
`preprint2.md` where they appear and compares them to the CSVs that
`make_figures_paper2.py` writes beside each figure.

A mismatch is an error, not a warning. Run it after editing either side.
"""
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figures')
DOC = os.path.join(HERE, 'preprint2.md')


def load(name):
    with open(os.path.join(FIG, f'{name}.csv')) as fh:
        return list(csv.DictReader(fh))


def num(rows, where, col):
    """Single value from the row matching the `where` dict."""
    hits = [r for r in rows if all(r[k] == v for k, v in where.items())]
    if len(hits) != 1:
        raise SystemExit(f"check: {where} matched {len(hits)} rows, expected 1")
    return float(hits[0][col])


def mean(rows, where, col):
    hits = [float(r[col]) for r in rows
            if all(r[k] == v for k, v in where.items())]
    if not hits:
        raise SystemExit(f"check: {where} matched nothing")
    return sum(hits) / len(hits)


doc = open(DOC, encoding='utf-8').read()
fails, checks = [], 0


def says(pattern, label):
    """The manuscript must contain `pattern` (a regex)."""
    global checks
    checks += 1
    if not re.search(pattern, doc):
        fails.append(f"{label}: manuscript does not contain /{pattern}/")


def close(a, b, tol, label):
    global checks
    checks += 1
    if abs(a - b) > tol:
        fails.append(f"{label}: figure {a:.6g} vs text {b:.6g} "
                     f"(tolerance {tol:g})")


# ---- figure 1: redistribution ------------------------------------------
f1 = load('fig2_redistribution')
cort = mean(f1, {'cell': 'cortical', 'arm': 'on'}, 'frac_at_min')
hh = mean(f1, {'cell': 'Hodgkin-Huxley', 'arm': 'on'}, 'frac_at_min')
cort_d = mean(f1, {'cell': 'cortical', 'arm': 'on'}, 'deficit_pct')
hh_d = mean(f1, {'cell': 'Hodgkin-Huxley', 'arm': 'on'}, 'deficit_pct')
close(cort, 0.222, 0.001, 'fig6 cortical fraction at zero')
close(hh, 0.153, 0.001, 'fig6 HH fraction at zero')
close(cort_d, -2.15, 0.02, 'fig6 cortical budget deficit')
close(hh_d, -10.15, 0.02, 'fig6 HH budget deficit')
says(r'\*\*22% on the cortical model\*\*', 'fig6 cortical % in text')
says(r'15%', 'fig6 HH % in text')
says(r'2\.15%', 'fig6 cortical deficit in text')
says(r'10\.15%', 'fig6 HH deficit in text')
says(r'0 of 4 wirings', 'fig6 survival in text')

# ---- figure 2: the two signs -------------------------------------------
f2 = load('fig3_two_signs')
close(num(f2, {'panel': 'A'}, 'rescued'), 31, 0, 'fig6 PY 1 rescued')
close(num(f2, {'panel': 'A'}, 'dead_off'), 33, 0, 'fig6 PY 1 opportunities')
close(num(f2, {'panel': 'A'}, 'wilson_lo'), 0.804, 0.0005, 'fig6 PY 1 CI low')
close(num(f2, {'panel': 'B', 'g_ie': '0.045'}, 'eliminated'), 6, 0,
      'fig6 AB/PD 2 eliminated')
close(num(f2, {'panel': 'B', 'g_ie': '0.045'}, 'alive_off'), 11, 0,
      'fig6 AB/PD 2 survivors')
close(num(f2, {'panel': 'B', 'g_ie': '0.06'}, 'eliminated'), 4, 0,
      'fig6 AB/PD 2 at 0.060')
close(num(f2, {'panel': 'B', 'g_ie': '0.03'}, 'dead_off'), 0, 0,
      'fig6 empty table at 0.030')
says(r'31/33 = 0\.939', 'fig6 PY 1 in text')
says(r'\[0\.804, 0\.983\]', 'fig6 CI in text')
says(r'eliminates 6 of those 11', 'fig6 AB/PD 2 in text')
says(r'eliminates 4 of', 'fig6 0.060 in text')

# ---- figure 3: the ordering --------------------------------------------
f3 = load('fig4_ordering')
rho = num(f3, {'cell': 'spearman_rho_all_seven'}, 'spikes_per_s')
close(rho, 0.852, 0.0005, 'fig6 Spearman')
for cell, rate, frac in (('PY 1', 25.6, 0.000), ('AB/PD 1', 482.7, 0.000),
                         ('LP 3', 517.7, 0.000), ('AB/PD 5', 536.2, 0.300),
                         ('AB/PD 3', 549.9, 0.200), ('AB/PD 4', 610.8, 0.273),
                         ('AB/PD 2', 756.2, 0.545)):
    close(num(f3, {'cell': cell}, 'spikes_per_s'), rate, 0.05,
          f'fig6 {cell} rate')
    close(num(f3, {'cell': cell}, 'fraction'), frac, 0.0005,
          f'fig6 {cell} elimination')
    says(re.escape(f'{rate}'), f'fig6 {cell} rate in text')
says(r'\+0\.852', 'fig6 Spearman in text')
gap = num(f3, {'cell': 'unsampled_gap_fold'}, 'spikes_per_s')
close(gap, 18.9, 0.1, 'fig6 unsampled gap')
says(r'nineteenfold gap', 'fig6 gap wording in text')
says(r'p = 0\.0123', 'fig6 exact p in text')

# ---- figure 4: fragility vs rescue -------------------------------------
f4 = load('fig5_fragility_vs_rescue')
close(num(f4, {'condition': 'fragility: AB/PD 1 no current'}, 'dead_off'), 0, 0,
      'fig6 baseline fragility')
close(num(f4, {'condition': 'fragility: AB/PD 1 -0.05'}, 'dead_off'), 12, 0,
      'fig6 driven fragility')
close(num(f4, {'condition': 'rescue: AB/PD 1 driven'}, 'rescued'), 0, 0,
      'fig6 driven rescue')
close(num(f4, {'condition': 'rescue: AB/PD 1 driven'}, 'n'), 21, 0,
      'fig6 driven opportunities')
close(num(f4, {'condition': 'rescue: PY 1'}, 'rescued'), 31, 0,
      'fig6 PY 1 rescue')
says(r'0/10', 'fig6 baseline in text')
says(r'12/12', 'fig6 driven in text')
says(r'\*\*0 of 21\*\*', 'fig6 driven rescue in text')

# ---- figure 5: the decomposition ---------------------------------------
f5 = load('fig6_decomposition')
pool_k = sum(num(f5, {'arm': a}, 'rescued')
             for a in ('Na', 'CaS', 'KCa', 'Kd', 'H', 'leak'))
pool_n = sum(num(f5, {'arm': a}, 'dead_off')
             for a in ('Na', 'CaS', 'KCa', 'Kd', 'H', 'leak'))
close(pool_k, 2, 0, 'fig6 pooled single-substitution rescues')
close(pool_n, 42, 0, 'fig6 pooled opportunities')
close(num(f5, {'arm': 'full'}, 'rescued'), 6, 0, 'fig6 full rescues')
close(num(f5, {'arm': 'full'}, 'dead_off'), 6, 0, 'fig6 full opportunities')
close(num(f5, {'arm': 'H'}, 'drive'), -0.30, 1e-9, 'fig6 H drive')
close(num(f5, {'arm': 'leak'}, 'drive'), -0.14, 1e-9, 'fig6 leak drive')
says(r'2 of 42', 'fig6 pooled in text')
says(r'6 of 6', 'fig6 full in text')
says(r'−0\.30|-0\.30', 'fig6 H drive in text')

h2 = {r['horizon_s'] for r in load('fig3_two_signs') if r['g_ie'] == '0.045'}
checks += 1
if h2 != {'35'}:
    fails.append(f'fig6 main panels are not all at 35 s: {h2}')

f0 = load('fig1_concept')
close(num(f0, {'box': 'left'}, 'numerator'), 31, 0, 'fig1 concept rescued')
close(num(f0, {'box': 'left'}, 'denominator'), 33, 0, 'fig1 concept opportunities')
close(num(f0, {'box': 'right'}, 'numerator'), 6, 0, 'fig1 concept eliminated')
close(num(f0, {'box': 'right'}, 'denominator'), 11, 0, 'fig1 concept survivors')

print(f"{checks} checks against {os.path.basename(DOC)}")
if fails:
    print(f"\n{len(fails)} MISMATCH(ES):")
    for f in fails:
        print(f"  - {f}")
    sys.exit(1)
print("all figure numbers match the manuscript")
