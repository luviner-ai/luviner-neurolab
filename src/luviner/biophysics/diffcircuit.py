"""
Differentiable circuits: solving for a circuit instead of tuning it.

Every circuit in this package was found by hand. Pick conductances, run
it, look at the output, adjust, repeat. The feed-forward inhibition
result took a dozen sweeps; the competition result took several more.
That is how circuits are built everywhere, because no simulator can
differentiate through a network of conductance-based neurons.

This module can. The state is a trajectory of ordinary differential
equations with no discontinuities anywhere, so a gradient runs from any
scalar measured on the output all the way back to every conductance in
the circuit, and the question inverts: instead of "what does this
circuit do?", ask **"what circuit does this?"**

The one obstacle is the spike. A conventional synapse is driven by
discrete presynaptic spikes, and a discrete event has no derivative.
The way round it is not a numerical trick but a real synapse type:
**graded release**. Transmitter output rises smoothly with presynaptic
voltage rather than waiting for a threshold crossing, and synapses that
work this way are found in the retina, in many invertebrate circuits,
and throughout C. elegans -- which has no action potentials at all.

    release(V_pre) = sigmoid((V_pre - V_half) / k)

So the differentiable circuit is not an approximation of the spiking
one. It is a different, equally biological, signalling regime -- and the
one that happens to admit gradients.

Reference for graded synaptic transmission:
    Juusola M, French AS, Uusitalo RO, Weckstrom M (1996). "Information
        processing by graded-potential transmission through tonically
        active synapses." Trends Neurosci 19:292-297.
    Liu Q, Hollopeter G, Jorgensen EM (2009). "Graded synaptic
        transmission at the C. elegans neuromuscular junction."
        PNAS 106:10823-10828.
"""

import numpy as np

from . import autodiff as ad
from .autodiff import Var


class GradedSynapse:
    """A synapse whose transmitter release varies smoothly with voltage.

    Parameters
    ----------
    g_max : peak conductance, stored in log space so gradient descent
        cannot make it negative.
    E_rev : reversal potential. Excitatory near 0 mV, inhibitory below
        rest.
    V_half, k : the release curve — half-maximal presynaptic voltage and
        its steepness.
    tau : how fast the postsynaptic conductance follows release.
    """

    def __init__(self, g_max=0.1, E_rev=0.0, V_half=-40.0, k=5.0, tau=3.0,
                 trainable=True, train_threshold=False):
        self.log_g = Var(np.log(max(g_max, 1e-6)), requires_grad=trainable)
        # The release threshold is a Var too, so an optimizer can move
        # *when* a synapse starts transmitting as well as how strongly.
        # It has to be used directly inside release(); keeping a plain
        # float alongside it silently detaches the gradient.
        self.V_half = Var(V_half, requires_grad=train_threshold)
        self.E_rev = E_rev
        self.k = k
        self.tau = tau

    @property
    def parameters(self):
        return [p for p in (self.log_g, self.V_half) if p.requires_grad]

    @property
    def threshold(self):
        return float(self.V_half.value)

    @property
    def g_max(self):
        return float(np.exp(self.log_g.value))

    def release(self, V_pre):
        """Fraction of maximal transmitter release, in [0, 1]."""
        return ad.sigmoid((V_pre - self.V_half) * (1.0 / self.k))

    def clamp(self, lo=-70.0, hi=-15.0):
        """Keep the release threshold inside a physiological range."""
        if self.V_half.requires_grad:
            self.V_half.value = np.clip(self.V_half.value, lo, hi)
        return self

    def current(self, s, V_post):
        """Synaptic current given the conductance state `s`."""
        return ad.exp(self.log_g) * s * (V_post - self.E_rev)


class DifferentiableCircuit:
    """A small network of conductance-based cells, end-to-end differentiable.

    Cells use the cortical (Traub-Miles) kinetics. Each connection is a
    GradedSynapse with its own first-order conductance state.
    """

    def __init__(self, n_neurons, g_Na=50.0, g_K=5.0, g_L=0.1,
                 E_Na=50.0, E_K=-90.0, E_L=-70.0, C_m=1.0, V_T=-63.0):
        self.n = n_neurons
        self.g_Na, self.g_K, self.g_L = g_Na, g_K, g_L
        self.E_Na, self.E_K, self.E_L = E_Na, E_K, E_L
        self.C_m, self.V_T = C_m, V_T
        self.connections = []          # (pre, post, GradedSynapse)

    def connect(self, pre, post, synapse):
        if not (0 <= pre < self.n and 0 <= post < self.n):
            raise IndexError("neuron index out of range")
        self.connections.append((pre, post, synapse))
        return synapse

    @property
    def parameters(self):
        return [p for _, _, s in self.connections for p in s.parameters]

    # ---- kinetics --------------------------------------------------------

    def _rates(self, V):
        x = V - self.V_T
        return (1.28 * ad.x_over_1_minus_exp_neg_x((x - 13.0) * (1.0 / 4.0)),
                1.4 * ad.x_over_1_minus_exp_neg_x(-(x - 40.0) * (1.0 / 5.0)),
                0.128 * ad.exp(-(x - 17.0) * (1.0 / 18.0)),
                4.0 / (1.0 + ad.exp(-(x - 40.0) * (1.0 / 5.0))),
                0.16 * ad.x_over_1_minus_exp_neg_x((x - 15.0) * (1.0 / 5.0)),
                0.5 * ad.exp(-(x - 10.0) * (1.0 / 40.0)))

    def initial_state(self):
        V = Var(np.full(self.n, self.E_L))
        am, bm, ah, bh, an, bn = self._rates(V)
        m = am / (am + bm)
        h = ah / (ah + bh)
        n = an / (an + bn)
        s = [Var(np.zeros(1)) for _ in self.connections]
        return V, m, h, n, s

    # ---- dynamics --------------------------------------------------------

    def _derivatives(self, V, m, h, n, s, I_ext):
        I_Na = self.g_Na * (m ** 3) * h * (V - self.E_Na)
        I_K = self.g_K * (n ** 4) * (V - self.E_K)
        I_L = self.g_L * (V - self.E_L)

        # synaptic currents, accumulated per postsynaptic cell
        I_syn = [Var(np.zeros(1)) for _ in range(self.n)]
        ds = []
        for (pre, post, syn), state in zip(self.connections, s):
            I_syn[post] = I_syn[post] + syn.current(state, V[post:post + 1])
            ds.append((syn.release(V[pre:pre + 1]) - state) * (1.0 / syn.tau))
        I_total = ad.stack([c.sum() for c in I_syn])

        am, bm, ah, bh, an, bn = self._rates(V)
        return ((I_ext - I_Na - I_K - I_L - I_total) * (1.0 / self.C_m),
                am * (1.0 - m) - bm * m,
                ah * (1.0 - h) - bh * h,
                an * (1.0 - n) - bn * n,
                ds)

    def simulate(self, I_ext, duration, dt=0.05):
        """Integrate the circuit with forward Euler; returns V over time.

        Euler rather than RK4: the tape holds every intermediate, and a
        four-stage method quadruples it. At dt <= 0.05 ms with graded
        synapses the difference is small compared with the parameter
        changes an optimizer makes.
        """
        V, m, h, n, s = self.initial_state()
        I_ext = np.broadcast_to(np.asarray(I_ext, float), (self.n,))
        I = Var(I_ext)
        trace = []
        for _ in range(int(round(duration / dt))):
            trace.append(V)
            dV, dm, dh, dn, ds = self._derivatives(V, m, h, n, s, I)
            V = V + dt * dV
            m = m + dt * dm
            h = h + dt * dh
            n = n + dt * dn
            s = [a + dt * b for a, b in zip(s, ds)]
        return ad.stack(trace)


def fully_connected(n_neurons, g_init=0.05, excitatory_fraction=0.5,
                    V_half=-50.0, k=4.0, tau=5.0, seed=0, types=None):
    """A circuit with every possible connection, all weights trainable.

    Rather than wiring a hypothesis and tuning it, this hands the
    optimizer the whole space and lets it decide which connections
    matter. A synapse whose conductance is driven to zero has been
    removed by the search; what survives is the structure the function
    required.

    Each connection is assigned an excitatory or inhibitory reversal
    potential at construction — sign is a property of the presynaptic
    cell in biology, not something a single synapse flips — so the
    optimizer chooses strengths, not signs.
    """
    rng = np.random.RandomState(seed)
    circuit = DifferentiableCircuit(n_neurons)
    # `types` fixes which cells are excitatory, so repeated searches
    # differ only in weight initialization. Left random, every run gets a
    # different cell composition and the structures found cannot be
    # compared.
    if types is None:
        types = [rng.rand() < excitatory_fraction for _ in range(n_neurons)]
    for pre in range(n_neurons):
        excitatory = bool(types[pre])
        for post in range(n_neurons):
            if pre == post:
                continue
            circuit.connect(pre, post, GradedSynapse(
                g_max=g_init * (0.5 + rng.rand()),
                E_rev=0.0 if excitatory else -80.0,
                V_half=V_half, k=k, tau=tau, train_threshold=True))
    return circuit


def step_response(circuit, neuron, amplitude, duration=120.0, dt=0.05,
                  onset_fraction=0.25):
    """Response of one cell to a current step, split into early and late.

    Returns (early, late) mean voltages after the step arrives. A unit
    that adapts has early > late; one that does not has early == late.
    """
    n_steps = int(round(duration / dt))
    onset = int(n_steps * onset_fraction)
    V, m, h, n, s = circuit.initial_state()
    trace = []
    off = Var(np.zeros(circuit.n))
    on = Var(np.full(circuit.n, amplitude))
    for i in range(n_steps):
        trace.append(V)
        I = on if i >= onset else off
        dV, dm, dh, dn, ds = circuit._derivatives(V, m, h, n, s, I)
        V = V + dt * dV
        m = m + dt * dm
        h = h + dt * dh
        n = n + dt * dn
        s = [a + dt * b for a, b in zip(s, ds)]
    V_trace = ad.stack(trace)
    early_lo = onset + int(0.05 * n_steps)
    early_hi = onset + int(0.20 * n_steps)
    late_lo = onset + int(0.50 * n_steps)
    return (V_trace[early_lo:early_hi, neuron].mean(),
            V_trace[late_lo:, neuron].mean())


def sequence_response(circuit, readout, inputs, duration=90.0, dt=0.03,
                      amplitude=3.0, pulse_ms=12.0, measure_ms=40.0):
    """Response of one cell to timed pulses delivered to several cells.

    `inputs` maps neuron index -> onset time in ms. The same two pulses
    delivered in opposite orders is the simplest test of whether a
    circuit can tell temporal order at all — a property no feed-forward
    summation has, since summation is commutative.

    Returns the peak depolarization of the readout above its baseline.
    """
    n_steps = int(round(duration / dt))
    V, m, h, n, s = circuit.initial_state()
    on_off = np.zeros((n_steps, circuit.n))
    for neuron, onset in inputs.items():
        lo = int(onset / dt)
        hi = int((onset + pulse_ms) / dt)
        on_off[lo:hi, neuron] = amplitude
    trace = []
    for i in range(n_steps):
        trace.append(V)
        I = Var(on_off[i])
        dV, dm, dh, dn, ds = circuit._derivatives(V, m, h, n, s, I)
        V = V + dt * dV
        m = m + dt * dm
        h = h + dt * dh
        n = n + dt * dn
        s = [a + dt * b for a, b in zip(s, ds)]
    V_trace = ad.stack(trace)
    baseline = float(V_trace.value[:int(10.0 / dt), readout].mean())

    # Measure inside a FIXED window after the last pulse, not over the
    # whole trace. A soft maximum taken over everything is an average in
    # disguise, so its value shifts with the total duration: measured on
    # the same circuit, 70 ms gave +0.0128 and 90 ms gave -0.0053 — the
    # sign flipped because the shorter run truncated the response being
    # measured. Anchoring the window to the stimulus makes the result
    # independent of how long the simulation happens to run.
    last_onset = max(inputs.values())
    lo = int(last_onset / dt)
    hi = min(int((last_onset + measure_ms) / dt), n_steps)
    scaled = (V_trace[lo:hi, readout] - baseline) * (1.0 / 8.0)
    return ad.log(ad.exp(scaled).mean()) * 8.0


def mean_depolarization(V_trace, neuron, skip_fraction=0.4):
    """Average voltage of one cell over the latter part of a run.

    A smooth, differentiable stand-in for firing rate: a cell driven
    harder sits more depolarized on average, whether or not one counts
    its spikes.
    """
    start = int(V_trace.value.shape[0] * skip_fraction)
    return V_trace[start:, neuron].mean()
