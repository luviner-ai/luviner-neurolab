# Paper 2 — the same plasticity rule rescues networks built from one conductance set and destabilises networks built from another

Manuscript, figures and the registration record for the second preprint from
this line of work. The simulator it runs on is the repository's `src/luviner/`
and the runs are in `../experiments/`; nothing is duplicated here.

| file | what it is |
|---|---|
| `preprint2.md` | the manuscript |
| `PAPER-2-SUPPLEMENTARY-S1.md` | audit tables, execution record, acceptance gates, and every wiring at both survival horizons |
| `LEARNING-REDISTRIBUTES.md` | the pre-registration ledger: for each result, the commit that registered the design and falsifiers before the numbers existed |
| `figures/` | the six figures as PDF and PNG, each with the CSV of exactly the values it plots |
| `make_figures_paper2.py` | regenerates all six from the run outputs |
| `check_figures_paper2.py` | 68 assertions joining every plotted number to the manuscript text |
| `check_stats_paper2.py` | 71 assertions recomputing every statistic from the runs in `../experiments/` |

## Reproducing the figures

From this directory, with Python 3.12, numpy 2.5.3 and matplotlib 3.11.1:

    python make_figures_paper2.py     # writes figures/fig1..fig6
    python check_figures_paper2.py    # exits non-zero on any mismatch
    python check_stats_paper2.py      # the same, for every statistic

The figure checker cannot catch a wrong p-value, because no p is plotted. The
stats checker recomputes every count, interval and test from the run JSON —
never from the figure CSVs, which would check one derivation against itself —
and lists at the end what it does not reach. Between them they found a rank
correlation whose ties were broken by input order and two p quoted from the
one-sided tail; §7.2 of the preprint records both.

The checker reads the CSV each figure writes beside itself and compares those
values against the manuscript. It is included because it is what caught three
number errors and one figure whose annotation disagreed with the text during
revision: a reader can re-run the check rather than take it on trust.

## Re-running the simulations

The run scripts are in `../experiments/coll*/`, each self-contained, each
writing JSON beside itself and skipping a lot whose JSON already exists. From
the directory of a script:

    cd ../experiments/coll6  && python coll6.py     # the rescue conditional
    cd ../experiments/coll8  && python coll8a.py    # five of the seven cells
    cd ../experiments/coll9  && python coll9.py     # the within-cell manipulation
    cd ../experiments/coll11 && python coll11a.py   # the substitution grid
    cd ../experiments/coll14 && python coll14.py a  # both endpoints at 35 s
    cd ../experiments/coll12 && python coll12.py A  # the robustness test

Delete the JSON first to reproduce it rather than skip it. Runs are
deterministic: the same seeds give the same spike times, which is what the
batched-versus-unbatched acceptance gates in `coll5/gate.py` and
`coll6/gate.py` check, bit for bit.

## Companion

The first preprint from the same simulator is in `../paper/`.
