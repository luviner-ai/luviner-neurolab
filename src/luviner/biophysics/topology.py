"""
Connectivity structure.

Up to here every circuit has been wired at random: pick a source, pick a
target, connect with fixed probability. That produces rhythm, because
rhythm needs only recurrent excitation and a brake. It cannot produce
function, because function needs to know *which* neuron connects to
*which*.

The generators here build structured connectivity, and they all return
the same thing -- arrays of (pre, post, weight) -- so a structured
network and a random one can be run through identical machinery and
compared. `shuffle_targets` is the control that makes such a comparison
mean something: it destroys structure while preserving the exact
multiset of weights and the in-degree of every neuron, so the only
difference is where the synapses land.

The ring architecture implements the classic continuous attractor for
spatial working memory: neurons are laid out on a circle, excitation is
local (a Gaussian in angular distance) and inhibition is global. A brief
input creates a localized hump of activity that outlasts it, and the
position of the hump encodes what the input was.

Reference:
    Compte A, Brunel N, Goldman-Rakic PS, Wang XJ (2000). "Synaptic
    mechanisms and network dynamics underlying spatial working memory in
    a cortical network model." Cereb Cortex 10:910-923.
    Ben-Yishai R, Bar-Or RL, Sompolinsky H (1995). PNAS 92:3844-3848.
"""

import numpy as np


def ring_positions(n):
    """Angular position of each neuron on a ring, in radians [0, 2pi)."""
    return np.arange(n) * 2.0 * np.pi / n


def angular_distance(a, b):
    """Shortest angular separation, in radians [0, pi]."""
    d = np.abs(a - b) % (2.0 * np.pi)
    return np.minimum(d, 2.0 * np.pi - d)


def ring_weights(n_pre, n_post=None, sigma_deg=43.2, w_min=0.0, w_max=1.0,
                 self_connections=False, threshold=0.05):
    """Gaussian ring connectivity.

    Weight falls off with angular distance between neurons:

        w(d) = w_min + (w_max - w_min) * exp(-d^2 / (2 sigma^2))

    Returns (pre, post, weight) arrays with weights below `threshold`
    (relative to w_max) dropped, so a broad profile does not cost N^2
    synapses.

    sigma_deg defaults to the footprint used by Compte et al. (2000).
    """
    n_post = n_pre if n_post is None else n_post
    theta_pre = ring_positions(n_pre)
    theta_post = ring_positions(n_post)
    d = angular_distance(theta_pre[:, None], theta_post[None, :])
    sigma = np.radians(sigma_deg)
    w = w_min + (w_max - w_min) * np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    if not self_connections and n_pre == n_post:
        np.fill_diagonal(w, 0.0)
    keep = w > threshold * w_max
    pre, post = np.nonzero(keep)
    return pre, post, w[keep]


def lateral_weights(n_pre, n_post, sigma_deg=45.0, invert=False, w_max=1.0,
                    threshold=0.12):
    """Distance-dependent weights between two ring populations.

    With `invert=False` the profile is Gaussian: strong nearby, weak far.
    With `invert=True` it is the complement -- weak nearby, strong far --
    which is what turns inhibition from a global brake into a competitive
    one.

    The distinction matters more than the amount. Inhibition applied
    uniformly subtracts the same quantity from every unit, so it
    compresses the differences between them; measured here, raising a
    global inhibitory conductance sevenfold left the ratio between a
    strong and a weak response unchanged at 0.50. Lateral inhibition,
    where each region suppresses the *others*, drives the same ratio to
    0.06 -- selection rather than attenuation, and the winner ends up
    firing faster than it would with no competition at all, because it is
    the one unit its own inhibition spares.

    This is the circuit-level counterpart of the dendritic result: what
    matters is not how much input arrives but where it goes.
    """
    theta_pre = ring_positions(n_pre)
    theta_post = ring_positions(n_post)
    d = angular_distance(theta_pre[:, None], theta_post[None, :])
    g = np.exp(-(d ** 2) / (2.0 * np.radians(sigma_deg) ** 2))
    w = (1.0 - g) * w_max if invert else g * w_max
    keep = w > threshold * w_max
    pre, post = np.nonzero(keep)
    return pre, post, w[keep]


def uniform_weights(pre_ids, post_ids, weight, probability=1.0, seed=0,
                    self_connections=False):
    """All-to-all or random connectivity with one shared weight."""
    rng = np.random.RandomState(seed)
    pre, post = np.meshgrid(np.asarray(pre_ids), np.asarray(post_ids),
                            indexing='ij')
    pre, post = pre.ravel(), post.ravel()
    keep = np.ones(len(pre), dtype=bool)
    if not self_connections:
        keep &= pre != post
    if probability < 1.0:
        keep &= rng.rand(len(pre)) < probability
    pre, post = pre[keep], post[keep]
    return pre, post, np.full(len(pre), float(weight))


def shuffle_targets(pre, post, weight, seed=0):
    """Randomly permute where synapses land, keeping everything else.

    Each neuron keeps its exact out-degree and its exact set of outgoing
    weights; only the identity of the targets is permuted. Total synaptic
    drive into the network is unchanged. This is the control for any
    claim that a result comes from structure rather than from the amount
    of connectivity.
    """
    rng = np.random.RandomState(seed)
    return np.asarray(pre), rng.permutation(np.asarray(post)), np.asarray(weight)


# ---- reading out a bump --------------------------------------------------

def population_vector(rates, positions=None):
    """Direction and strength of activity on a ring.

    Returns (angle in radians, coherence in [0, 1]). Coherence near 0
    means activity is spread all round the ring; near 1 means it is
    concentrated at one place. This is how the network's "memory" is
    read out: the angle says what it is holding.
    """
    rates = np.asarray(rates, dtype=float)
    n = len(rates)
    theta = ring_positions(n) if positions is None else np.asarray(positions)
    total = rates.sum()
    if total <= 0:
        return 0.0, 0.0
    z = (rates * np.exp(1j * theta)).sum() / total
    return float(np.angle(z) % (2.0 * np.pi)), float(np.abs(z))


def bump_width(rates, positions=None):
    """Circular standard deviation of the activity profile, in degrees.

    A tight bump gives a small number; activity spread evenly round the
    ring approaches the uniform-distribution value (~81 degrees).
    """
    _, coherence = population_vector(rates, positions)
    if coherence <= 0:
        return float('nan')
    coherence = min(coherence, 1.0 - 1e-12)
    return float(np.degrees(np.sqrt(-2.0 * np.log(coherence))))


def rates_in_window(spike_times, t_start, t_end, n_neurons=None):
    """Per-neuron firing rate over a time window, in Hz."""
    n = len(spike_times) if n_neurons is None else n_neurons
    span = (t_end - t_start) / 1000.0
    out = np.zeros(n)
    for i, s in enumerate(spike_times[:n]):
        s = np.asarray(s)
        out[i] = ((s >= t_start) & (s < t_end)).sum() / span
    return out
