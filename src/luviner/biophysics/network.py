"""
Networks of biophysical neurons connected by chemical synapses.

The network integrates every neuron and every synapse on a shared clock.
Per step:

  1. Each synapse's current is evaluated against the *present* postsynaptic
     voltage and subtracted from that neuron's input current.
  2. Every neuron advances one RK4 step.
  3. Every synapse advances one exponential step.
  4. Spikes are detected as upward crossings of a threshold and pushed
     into the axons of the outgoing synapses, where they wait out the
     conduction delay.

Synaptic coupling is evaluated once per step rather than inside the RK4
stages, so coupling is first-order accurate in dt even though each cell's
own dynamics are fourth-order. This is what NEURON and Brian do as well;
it is accurate because synaptic conductances change on a ~ms timescale
while dt is ~0.01 ms. Halving dt must not change the result appreciably,
which the test suite checks directly.
"""

import numpy as np

from .hodgkin_huxley import HodgkinHuxleyNeuron


class Network:
    """A set of neurons wired by conductance-based synapses.

    Parameters
    ----------
    neurons : list of HodgkinHuxleyNeuron, or an int (that many defaults)
    spike_threshold : mV, upward crossing counts as a spike
    """

    def __init__(self, neurons, spike_threshold=0.0):
        if isinstance(neurons, int):
            neurons = [HodgkinHuxleyNeuron() for _ in range(neurons)]
        self.neurons = list(neurons)
        self.spike_threshold = spike_threshold
        self.connections = []   # (pre_idx, post_idx, synapse)
        self.reset()

    @property
    def n_neurons(self):
        return len(self.neurons)

    def connect(self, pre, post, synapse):
        """Wire pre -> post. Returns the synapse for convenience.

        A neuron may synapse onto itself (autapse) and any pair may be
        connected more than once, as in real tissue.
        """
        if not (0 <= pre < self.n_neurons and 0 <= post < self.n_neurons):
            raise IndexError("neuron index out of range")
        self.connections.append((pre, post, synapse))
        return synapse

    def reset(self, V0=None):
        """Put every neuron at its resting state and clear all synapses."""
        for n in self.neurons:
            n.reset(n.resting_potential() if V0 is None else V0)
        for _, _, s in self.connections:
            s.reset()
        self._t = 0.0
        self._prev_V = np.array([n.V for n in self.neurons])
        self.spike_times = [[] for _ in range(self.n_neurons)]
        return self

    # ---- dynamics --------------------------------------------------------

    def synaptic_currents(self):
        """Total synaptic current into each neuron, uA/cm^2."""
        I = np.zeros(self.n_neurons)
        for _, post, syn in self.connections:
            I[post] += syn.current(self.neurons[post].V)
        return I

    def step(self, I_ext=None, dt=0.01):
        """Advance the whole network one timestep."""
        n = self.n_neurons
        I_ext = np.zeros(n) if I_ext is None else np.broadcast_to(
            np.asarray(I_ext, dtype=float), (n,))

        # Synaptic current opposes the membrane like any ionic current.
        I_syn = self.synaptic_currents()
        for i, neuron in enumerate(self.neurons):
            neuron.step(I_ext[i] - I_syn[i], dt)
        for _, _, syn in self.connections:
            syn.step(dt)

        self._t += dt
        V = np.array([nn.V for nn in self.neurons])
        fired = (self._prev_V <= self.spike_threshold) & (V > self.spike_threshold)
        if fired.any():
            for i in np.flatnonzero(fired):
                self.spike_times[i].append(self._t)
            for pre, _, syn in self.connections:
                if fired[pre]:
                    syn.receive_spike(self._t)
        self._prev_V = V
        return V

    def simulate(self, duration, dt=0.01, I_ext=None, record_synapses=False):
        """Run the network and record every membrane potential.

        I_ext may be None, a scalar, a per-neuron array, a callable
        t -> array, or an array of shape (n_steps, n_neurons).
        """
        n_steps = int(round(duration / dt))
        n = self.n_neurons

        if I_ext is None:
            current_at = lambda i: None
        elif callable(I_ext):
            current_at = lambda i: I_ext(i * dt)
        else:
            arr = np.asarray(I_ext, dtype=float)
            if arr.ndim == 2:
                current_at = lambda i: arr[i]
            else:
                current_at = lambda i: arr

        t = np.arange(n_steps) * dt
        V = np.empty((n_steps, n))
        G = (np.empty((n_steps, len(self.connections)))
             if record_synapses else None)

        for i in range(n_steps):
            V[i] = [nn.V for nn in self.neurons]
            if record_synapses:
                G[i] = [s.conductance(self.neurons[p].V)
                        for _, p, s in self.connections]
            self.step(current_at(i), dt)

        out = {'t': t, 'V': V, 'dt': dt,
               'spikes': [np.array(s) for s in self.spike_times]}
        if record_synapses:
            out['g_syn'] = G
        return out


# ---- analysis helpers ----------------------------------------------------

def psp_amplitude(V, V_baseline):
    """Peak deflection of a postsynaptic potential from baseline, mV.

    Signed: positive for an EPSP, negative for an IPSP.
    """
    V = np.asarray(V)
    up = V.max() - V_baseline
    down = V.min() - V_baseline
    return up if abs(up) >= abs(down) else down


def transmission_delay(pre_spikes, post_spikes):
    """Latency from each presynaptic spike to the next postsynaptic one, ms.

    Returns an empty array if the postsynaptic cell never followed.
    """
    pre_spikes = np.asarray(pre_spikes)
    post_spikes = np.asarray(post_spikes)
    lags = []
    for tp in pre_spikes:
        later = post_spikes[post_spikes > tp]
        if len(later):
            lags.append(later[0] - tp)
    return np.array(lags)
