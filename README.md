# A pre-registered computational neuroscience lab, in a git repository

This repository holds the complete materials for one study: the code, the
recorded results, the paper, and the laboratory notebook that was written
while the work was done.

The study asks whether spontaneous population bursts reinforce synaptic
weight differences that already exist, in a conductance-based model of a
bursting neuronal culture — 48 excitatory and 12 inhibitory
stomatogastric-ganglion model neurons with sparse recurrent connectivity,
short-term depression and pair-based STDP. The answer is that they do, and
that the statistic carrying it is not the one we expected: not first-spike
order, but a very small excess of pre-before-post spike *pairs* — 0.2–0.35
percentage points, accumulated over about 10⁴ pairs per synapse in 60
simulated seconds. Shortening bursts amplifies the effect roughly tenfold.

We then pre-registered the corresponding prediction for real tissue and
tested it on 16 archived recordings of developing cortical cultures
(Wagenaar, Pine & Potter, 2006). **It was not detected.** That negative
result is here in full, with the exact bound it supports.

What is unusual is the method, not the scale. Every prediction, every
falsifier and every stopping rule was committed to a version-controlled
record *before* the run it governs, by a founder working with AI agents on
a laptop. `LAB-RECORD.md` is the curated export of that notebook, each
section carrying the hash of the commit that introduced it. Corrections are
appended rather than applied, so several sections end by contradicting
their own opening — including three published numbers that a second,
independent recomputation caught after the fact.

## What was found, and how far it goes

| Finding | Status |
|---|---|
| Stronger synapses receive greater net potentiation across bursts | **Replicated** — 10/10 wirings in an out-of-sample cohort. The pre-registered *magnitude* criterion was **not** met; only the direction |
| The carrier is a pre-before-post **pair count**, not kernel-weighted timing | **Established in this model** — the exponential kernel adds nothing beyond the count (ratio CV 0.13) |
| First-spike latency order carries the sorting | **Refuted**, in the model and again in tissue |
| Shortening bursts strengthens sorting (~10×) | **Pre-registered intervention, met** in 5/5 wirings, and reproduced by a second, physiologically distinct lever |
| Burst *duration* rather than per-burst spike count is the operative variable | **Not established.** The two cannot be separated here: they correlate at ρ = +0.93 and the gap between them is indistinguishable from zero |
| Both levers show a silent absorbing state at *intermediate* values | **Non-monotone in both**, on the same seven wirings — so not two independent confirmations |
| The same strength-ordered excess appears in cortical cultures | **Not detected.** A bounded null: the 16 recordings exclude effects above +0.32 pp but cannot resolve the model's ~0.2 pp |
| The predicted negative relationship with burst duration appears in tissue | **Refuted**, and adequately powered under the pre-registered criterion |

All interventional evidence is computational. The study does **not**
establish that the mechanism operates in living cortical networks. The
decisive test is an intervention on cultured networks with an
independently validated measure of synaptic strength, which this work does
not contain.

## Reproducing

```sh
./reproduce.sh          # venv, dependencies, regenerate all seven figures
```

The figures are rebuilt from the JSON results committed here, not from a
fresh simulation — the recorded numbers are the artifact. To re-run the
simulations themselves, the scripts in `experiments/dish/` are the ones
that produced those files; each carries its registration in its docstring.

For the tissue analysis, `experiments/mea/fetch.sh` downloads the public
recordings and verifies them against `manifest.tsv`; the raw data belongs
to its authors and is not redistributed here. `run.sh` then reproduces the
analysis in about nine minutes on one core.

## Layout

    paper/          the preprint, its figures, and the scripts that build both
    src/            the simulation library (biophysics only)
    experiments/    the registered runs and the JSON they wrote
    LAB-RECORD.md   the notebook: predictions, falsifiers, verdicts, corrections

## Citing

DOI (concept, always the latest release): [10.5281/zenodo.22287792](https://doi.org/10.5281/zenodo.22287792); this release (v1.0.1): [10.5281/zenodo.22287793](https://doi.org/10.5281/zenodo.22287793).

Preprint: [10.64898/2026.09.03.749248](https://doi.org/10.64898/2026.09.03.749248).

If you use the data or the code, please cite the preprint (see
`CITATION.cff`) and, for the culture recordings, Wagenaar, Pine & Potter
(2006), *BMC Neuroscience* 7:11.

## Licence

Code under MIT (`LICENSE`); the paper text and figures under CC-BY 4.0
(`LICENSE-CC-BY`).

More about the authors' other work: <https://luviner.com>
