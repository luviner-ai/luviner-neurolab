"""
Population-level analysis: rates, rhythms, and synchrony.

A circuit does things a single neuron does not, and those things are
visible only in the aggregate. These are the measurements for that:
the population firing rate over time, its spectrum, and whether an
apparent rhythm survives a surrogate control.

The surrogate test matters more than it looks. A firing-rate signal
estimated over a short window has a ragged spectrum, and its largest
peak can sit well above the mean power purely by chance. Reporting that
peak as an oscillation is the easiest mistake to make here, so
`is_oscillatory` compares against spike trains shuffled to destroy
synchrony while preserving spike count.
"""

import numpy as np


def population_rate(spike_times, duration, bin_ms=1.0, n_neurons=None):
    """Population firing rate in Hz per neuron, binned over time.

    `spike_times` is a list of per-neuron spike time arrays, as returned
    by Network.simulate or NeuronPopulation.simulate.
    """
    n = len(spike_times) if n_neurons is None else n_neurons
    edges = np.arange(0.0, duration + bin_ms, bin_ms)
    present = [s for s in spike_times if len(s)]
    allsp = np.concatenate(present) if present else np.zeros(0)
    counts, _ = np.histogram(allsp, bins=edges)
    return edges[:-1], counts / (n * bin_ms / 1000.0)


def rate_spectrum(rate, bin_ms=1.0, f_lo=5.0, f_hi=200.0):
    """Power spectrum of a population rate signal.

    Returns (frequencies, power) restricted to [f_lo, f_hi]. The mean is
    removed and a Hann window applied, so the DC component and edge
    discontinuities do not dominate.
    """
    r = np.asarray(rate, dtype=float)
    r = r - r.mean()
    if len(r) < 8:
        return np.zeros(0), np.zeros(0)
    power = np.abs(np.fft.rfft(r * np.hanning(len(r)))) ** 2
    freq = np.fft.rfftfreq(len(r), d=bin_ms / 1000.0)
    band = (freq >= f_lo) & (freq <= f_hi)
    return freq[band], power[band]


def spectral_peak(rate, bin_ms=1.0, f_lo=5.0, f_hi=200.0):
    """Return (peak frequency in Hz, peak-to-mean power ratio).

    The ratio is a coherence proxy: how far the strongest rhythm stands
    above the rest of the band. It is meaningful only when compared
    against a surrogate -- see `is_oscillatory`.
    """
    f, p = rate_spectrum(rate, bin_ms, f_lo, f_hi)
    if len(f) == 0 or p.mean() == 0:
        return 0.0, 0.0
    return float(f[np.argmax(p)]), float(p.max() / p.mean())


def _ratio_of(spikes, duration, bin_ms, n):
    _, rate = population_rate(spikes, duration, bin_ms, n)
    return spectral_peak(rate, bin_ms)


def _uniform_surrogate(spike_times, duration, rng):
    """Redraw every spike time uniformly. Destroys everything: synchrony
    between cells AND the interval structure within each cell."""
    n_spikes = sum(len(s) for s in spike_times)
    return [rng.uniform(0.0, duration, size=n_spikes)]


def _shift_surrogate(spike_times, duration, rng):
    """Circularly shift each spike train by a random amount.

    Every neuron keeps its own spike count and its complete interval
    structure; only the alignment BETWEEN neurons is destroyed. This is
    the control that separates a genuine population rhythm from cells
    that merely happen to fire at similar rates.
    """
    out = []
    for s in spike_times:
        s = np.asarray(s, dtype=float)
        if len(s) == 0:
            out.append(s)
            continue
        out.append(np.sort((s + rng.uniform(0.0, duration)) % duration))
    return out


def is_oscillatory(spike_times, duration, bin_ms=1.0, n_surrogates=20,
                   n_neurons=None, seed=0, z_threshold=3.0,
                   surrogate='shift'):
    """Test a population rhythm against surrogate spike trains.

    Two surrogates are available, and they answer different questions.

    `surrogate='uniform'` redraws every spike time uniformly. It only
    asks whether the peak beats noise. A peak can pass this test purely
    because the cells fire at similar rates, with no synchrony at all --
    the interval structure of each train is destroyed along with the
    synchrony, so both contribute to the difference.

    `surrogate='shift'` (the default) circularly shifts each train by a
    random offset. Each neuron keeps its exact spike count and interval
    distribution; only the alignment between neurons is broken. Passing
    this means the rhythm lives in the relationship between cells, which
    is what "population oscillation" should mean.

    Returns the observed peak and power ratio, the surrogate
    distribution, the z score, and a verdict.
    """
    n = len(spike_times) if n_neurons is None else n_neurons
    freq, ratio = _ratio_of(spike_times, duration, bin_ms, n)

    present = [np.asarray(s) for s in spike_times if len(s)]
    if not present:
        return {'peak_hz': 0.0, 'ratio': 0.0, 'z': 0.0,
                'surrogate_mean': 0.0, 'surrogate_std': 0.0,
                'surrogate': surrogate, 'oscillatory': False}

    make = {'uniform': _uniform_surrogate,
            'shift': _shift_surrogate}[surrogate]
    rng = np.random.RandomState(seed)
    ratios = []
    for _ in range(n_surrogates):
        fake = make(present, duration, rng)
        ratios.append(_ratio_of(fake, duration, bin_ms, n)[1])
    ratios = np.array(ratios)
    sd = ratios.std()
    z = (ratio - ratios.mean()) / sd if sd > 0 else 0.0
    return {'peak_hz': freq, 'ratio': ratio, 'z': float(z),
            'surrogate_mean': float(ratios.mean()),
            'surrogate_std': float(sd),
            'surrogate': surrogate,
            'oscillatory': bool(z > z_threshold)}


def synchrony_index(V):
    """Golomb-Rinzel synchrony measure from membrane potentials.

    Variance of the population-averaged voltage over the mean of the
    per-neuron variances. 0 means asynchronous, 1 means fully
    synchronous. `V` has shape (n_steps, n_neurons).
    """
    V = np.asarray(V, dtype=float)
    mean_trace = V.mean(axis=1)
    num = mean_trace.var()
    den = V.var(axis=0).mean()
    return float(num / den) if den > 0 else 0.0
