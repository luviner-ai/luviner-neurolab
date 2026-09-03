"""
Synaptic plasticity: connections that change with use.

Every synapse so far has had a fixed conductance. That is enough for a
circuit to have dynamics, and even to hold a memory for a second in
persistent activity, but it is not enough for a circuit to *learn*
anything: switch the network off and everything it ever saw is gone.

Spike-timing dependent plasticity is the mechanism by which real
synapses change. What matters is not that two neurons were active, but
in which order and how close together:

    pre fires shortly BEFORE post   ->  strengthen  (LTP)
    pre fires shortly AFTER post    ->  weaken      (LTD)
    far apart in time               ->  no change

That asymmetry is causal, not correlational: it strengthens connections
where the presynaptic cell plausibly helped cause the postsynaptic
spike, and weakens the ones that arrived too late to have contributed.

Implemented with the standard trace formulation rather than by storing
spike pairs. Each synapse carries two decaying traces:

    dx_pre/dt  = -x_pre / tau_plus     x_pre  += 1 on each presynaptic spike
    dy_post/dt = -y_post / tau_minus   y_post += 1 on each postsynaptic spike

    on a presynaptic spike:   w -= A_minus * y_post
    on a postsynaptic spike:  w += A_plus  * x_pre

This is mathematically equivalent to the all-to-all pair rule and costs
O(1) per synapse per step, so it scales with the vectorized population.

Weight dependence: with `mu=0` updates are additive and weights must be
clipped, which drives them to the rails and makes the distribution
bimodal. With `mu=1` they are multiplicative -- a strong synapse
potentiates less -- which keeps a graded, unimodal distribution. Real
synapses sit somewhere between; `mu` selects where.

Two homeostatic processes are provided, and they answer different
questions. `SynapticScaling` holds a neuron's *firing rate* at a set
point; `WeightNormalization` holds the *total conductance* arriving at
it. The second is what keeps selectivity and gain from trading against
each other, which is the obstacle to building a usable cell assembly.

STDP on its own has no set point. It is driven by correlations and
knows nothing about how much a neuron is firing, so a cell that starts
slightly too excitable potentiates its inputs, fires more, and
potentiates them further. Real cortex does not do this, because a second
and much slower process watches the firing rate and scales every
incoming synapse of a neuron up or down together (Turrigiano et al.
1998). Scaling is *multiplicative*, which is the property that matters:
it changes the total drive a neuron receives while preserving the ratios
between its inputs, so it stabilises the rate without erasing what the
Hebbian rule learned. `SynapticScaling` implements it.

References:
    Turrigiano GG, Leslie KR, Desai NS, Rutherford LC, Nelson SB (1998).
        "Activity-dependent scaling of quantal amplitude in neocortical
        neurons." Nature 391:892-896.
    Bi GQ, Poo MM (1998). "Synaptic modifications in cultured hippocampal
        neurons: dependence on spike timing, synaptic strength, and
        postsynaptic cell type." J Neurosci 18:10464-10472.
    Song S, Miller KD, Abbott LF (2000). "Competitive Hebbian learning
        through spike-timing-dependent synaptic plasticity."
        Nat Neurosci 3:919-926.
    Morrison A, Diesmann M, Gerstner W (2008). "Phenomenological models
        of synaptic plasticity based on spike timing."
        Biol Cybern 98:459-478.
"""

import numpy as np


class STDP:
    """Pair-based spike-timing dependent plasticity with weight bounds.

    Parameters
    ----------
    A_plus, A_minus : learning rates for potentiation and depression, as
        a fraction of w_max. LTD is conventionally a little weaker per
        pair but acts over a slightly longer window, which makes the rule
        net-depressing for uncorrelated activity and so keeps runaway
        potentiation in check.
    tau_plus, tau_minus : ms, the widths of the two sides of the timing
        window. ~20 ms matches the hippocampal measurements.
    w_min, w_max : bounds on the conductance.
    mu : weight dependence. 0 = additive (bimodal weights),
         1 = multiplicative (graded weights).
    """

    def __init__(self, A_plus=0.008, A_minus=0.0088,
                 tau_plus=16.8, tau_minus=33.7,
                 w_min=0.0, w_max=1.0, mu=0.0):
        self.A_plus = A_plus
        self.A_minus = A_minus
        self.tau_plus = tau_plus
        self.tau_minus = tau_minus
        self.w_min = w_min
        self.w_max = w_max
        self.mu = mu

    # ---- the timing window ----------------------------------------------

    def window(self, dt_ms):
        """Weight change for a single pre/post pair separated by dt_ms.

        dt_ms > 0 means presynaptic first (causal, potentiating).
        Returned as a fraction of w_max, for comparison with the
        experimental curves which are plotted as percent change.
        """
        dt = np.asarray(dt_ms, dtype=float)
        return np.where(dt >= 0,
                        self.A_plus * np.exp(-dt / self.tau_plus),
                        -self.A_minus * np.exp(dt / self.tau_minus))

    def window_integral(self):
        """Net weight change per pair when the timing is broadly spread.

        The area under the window, `A_plus * tau_plus - A_minus *
        tau_minus`. It is the quantity that decides what the rule does to
        a group of cells that fire together but not in a fixed order:
        their relative spike times land on both sides of zero, so what
        survives is the integral, not the peak.

        Negative means the rule *dissolves* co-active groups rather than
        binding them -- which is the correct default for stability, and
        the reason a net-depressing window cannot by itself build a cell
        assembly, however long it is trained.
        """
        return self.A_plus * self.tau_plus - self.A_minus * self.tau_minus

    def expected_change(self, lags_ms):
        """Mean weight change per pair over an observed set of lags.

        `window_integral` assumes the lags are spread evenly. When they
        are not -- a protocol with a systematic order, or a network whose
        jitter is narrower than the window -- this is the honest version:
        evaluate the window on the lags the simulation actually produced.
        """
        lags = np.asarray(lags_ms, dtype=float)
        if lags.size == 0:
            return 0.0
        return float(np.mean(self.window(lags)))

    # ---- bounded updates -------------------------------------------------

    def _potentiation_scale(self, w):
        if self.mu == 0.0:
            return np.ones_like(w)
        return ((self.w_max - w) / (self.w_max - self.w_min)) ** self.mu

    def _depression_scale(self, w):
        if self.mu == 0.0:
            return np.ones_like(w)
        return ((w - self.w_min) / (self.w_max - self.w_min)) ** self.mu

    def potentiate(self, w, trace_pre):
        """Apply LTP driven by the presynaptic trace, in place-safe form."""
        dw = self.A_plus * self.w_max * trace_pre * self._potentiation_scale(w)
        return np.clip(w + dw, self.w_min, self.w_max)

    def depress(self, w, trace_post):
        """Apply LTD driven by the postsynaptic trace."""
        dw = self.A_minus * self.w_max * trace_post * self._depression_scale(w)
        return np.clip(w - dw, self.w_min, self.w_max)


class SynapticScaling:
    """Homeostatic scaling of the weights arriving at each neuron.

    The neuron measures its own firing rate with a slow trace and scales
    all of its incoming weights by a common factor:

        dw/dt = w * (r_target - r_post) / (r_target * tau_scale)

    A silent neuron therefore grows its inputs with time constant
    `tau_scale`, and one firing at twice the target shrinks them at the
    same rate. Because the factor is common to every synapse of that
    neuron, the *relative* weights are untouched: whatever contrast STDP
    has built survives, while the total is pulled to a set point.

    That combination is what makes Hebbian learning competitive. With a
    bounded total, a synapse that potentiates does so at the expense of
    the others on the same cell -- which is the mechanism that turns
    correlated input into a selective assembly rather than into uniform
    saturation.

    Parameters
    ----------
    target_rate : Hz, the set point.
    tau_rate : ms, the window over which the neuron estimates its own
        rate. Long compared with the interspike interval, short compared
        with `tau_scale`.
    tau_scale : ms, how fast the scaling acts. In cortex this is hours;
        here it is set to the timescale of the experiment, which is the
        usual compression in a simulation.
    w_min, w_max : bounds. Scaling is a separate process from the Hebbian
        rule and keeps its own, deliberately wider, ceiling.
    """

    def __init__(self, target_rate=5.0, tau_rate=200.0, tau_scale=500.0,
                 w_min=0.0, w_max=np.inf):
        if target_rate <= 0:
            raise ValueError("target_rate must be positive")
        self.target_rate = float(target_rate)
        self.tau_rate = float(tau_rate)
        self.tau_scale = float(tau_scale)
        self.w_min = w_min
        self.w_max = w_max
        self.rate = None

    def reset(self, n_neurons):
        """Start every neuron at the set point, so scaling begins neutral."""
        self.rate = np.full(int(n_neurons), self.target_rate)
        return self

    def step(self, weights, post, fired, dt):
        """Advance the rate estimate and scale `weights` accordingly.

        `weights` are the per-connection conductances, `post` their
        postsynaptic neuron indices, `fired` the current spike vector.
        """
        if self.rate is None:
            self.reset(len(fired))
        # Rate trace in Hz: each spike adds 1000/tau, so steady firing at
        # f Hz settles at f.
        self.rate *= np.exp(-dt / self.tau_rate)
        self.rate[fired] += 1000.0 / self.tau_rate
        if not len(weights):
            return weights
        error = (self.target_rate - self.rate[post]) / self.target_rate
        return np.clip(weights * (1.0 + dt * error / self.tau_scale),
                       self.w_min, self.w_max)


class WeightNormalization:
    """Heterosynaptic competition: a fixed total input budget per neuron.

    Scaling to a *rate* set point ties two things together that need to
    be separate. Raising the target rate is the only way to make the
    recurrent weights large enough to recruit a cell, but a higher target
    also makes the whole network more active during training, and once
    everything is co-active there is no longer anything for the Hebbian
    rule to distinguish. Selectivity and gain end up fighting each other.

    Normalising the *total* incoming conductance separates them. The
    budget fixes how much drive a neuron can receive; the Hebbian rule
    decides how that budget is divided. Contrast then costs nothing in
    gain, and gain costs nothing in contrast.

        w  ->  w * (budget / sum of that neuron's incoming weights)

    applied gradually with time constant `tau`, so it is a process the
    network runs rather than a projection imposed on it each step.

    Reference:
        Royer S, Pare D (2003). "Conservation of total synaptic weight
            through balanced synaptic depression and potentiation."
            Nature 422:518-522.
        Miller KD, MacKay DJC (1994). "The role of constraints in Hebbian
            learning." Neural Comput 6:100-126.
    """

    def __init__(self, budget=0.1, tau=100.0, w_min=0.0, w_max=np.inf):
        if budget <= 0:
            raise ValueError("budget must be positive")
        self.budget = float(budget)
        self.tau = float(tau)
        self.w_min = w_min
        self.w_max = w_max

    def reset(self, n_neurons):
        self.n_neurons = int(n_neurons)
        return self

    def totals(self, weights, post):
        """Incoming weight sum per postsynaptic neuron."""
        total = np.zeros(self.n_neurons)
        np.add.at(total, post, weights)
        return total

    def step(self, weights, post, fired, dt):
        if not len(weights):
            return weights
        total = self.totals(weights, post)
        factor = np.where(total > 0, self.budget / np.maximum(total, 1e-12), 1.0)
        return np.clip(weights * (1.0 + (dt / self.tau) * (factor[post] - 1.0)),
                       self.w_min, self.w_max)


class PlasticConnections:
    """A set of synapses whose weights follow an STDP rule.

    Holds per-connection weights and the pre/post traces, and updates
    them from spike vectors. Designed to sit alongside a
    NeuronPopulation: call `step` once per timestep with the boolean
    spike vector, and read `weights` to get the current conductances.
    """

    def __init__(self, pre, post, weights, rule=None, n_neurons=None,
                 homeostasis=None):
        self.pre = np.asarray(pre, dtype=int)
        self.post = np.asarray(post, dtype=int)
        self.weights = np.asarray(weights, dtype=float).copy()
        if not (len(self.pre) == len(self.post) == len(self.weights)):
            raise ValueError("pre, post and weights must have equal length")
        self.rule = rule or STDP()
        self.homeostasis = homeostasis
        n = (max(self.pre.max(), self.post.max()) + 1 if len(self.pre)
             else 0) if n_neurons is None else n_neurons
        self.n_neurons = int(n)
        self.reset_traces()

    def reset_traces(self):
        """Clear the eligibility traces, leaving weights untouched."""
        self.trace_pre = np.zeros(self.n_neurons)
        self.trace_post = np.zeros(self.n_neurons)
        if self.homeostasis is not None:
            self.homeostasis.reset(self.n_neurons)
        return self

    def step(self, fired, dt=0.05, learn=True):
        """Advance traces by dt and apply STDP for the spikes in `fired`.

        `fired` is a boolean array over neurons. Order matters and
        follows the standard convention: depression is evaluated on a
        presynaptic spike against the existing postsynaptic trace, and
        potentiation on a postsynaptic spike against the existing
        presynaptic trace, both before the traces are incremented by the
        current spikes. A pre and post spike in the *same* timestep
        therefore changes nothing: neither trace has been incremented
        yet, so both terms are zero. Exactly simultaneous firing is not a
        causal pair in either direction, and the rule treats it as one.
        """
        fired = np.asarray(fired, dtype=bool)
        if learn and len(self.weights):
            pre_fired = fired[self.pre]
            post_fired = fired[self.post]
            if post_fired.any():
                w = self.weights
                new = self.rule.potentiate(w, self.trace_pre[self.pre])
                self.weights = np.where(post_fired, new, w)
            if pre_fired.any():
                w = self.weights
                new = self.rule.depress(w, self.trace_post[self.post])
                self.weights = np.where(pre_fired, new, w)

        self.trace_pre *= np.exp(-dt / self.rule.tau_plus)
        self.trace_post *= np.exp(-dt / self.rule.tau_minus)
        self.trace_pre[fired] += 1.0
        self.trace_post[fired] += 1.0

        # Homeostasis runs on every step, learning or not: it is a
        # property of the cell, not of the training protocol.
        if self.homeostasis is not None:
            self.weights = self.homeostasis.step(self.weights, self.post,
                                                 fired, dt)
        return self

    # ---- readout ---------------------------------------------------------

    def weight_matrix(self):
        """Dense (n_neurons, n_neurons) view of the current weights."""
        W = np.zeros((self.n_neurons, self.n_neurons))
        np.add.at(W, (self.pre, self.post), self.weights)
        return W

    def mean_weight(self):
        return float(self.weights.mean()) if len(self.weights) else 0.0

    def bimodality(self):
        """Fraction of weights sitting within 10% of either bound.

        Additive STDP drives weights to the rails; multiplicative STDP
        keeps them graded. This number distinguishes the two regimes.
        """
        if not len(self.weights):
            return 0.0
        span = self.rule.w_max - self.rule.w_min
        low = self.weights < self.rule.w_min + 0.1 * span
        high = self.weights > self.rule.w_max - 0.1 * span
        return float((low | high).mean())


def pair_protocol(rule, dt_ms, n_pairs=60, pair_interval=1000.0, dt=0.1,
                  w0=0.5):
    """Simulate the classic pairing experiment on a single synapse.

    Repeatedly fires pre and post separated by `dt_ms` (positive = pre
    first), spaced far enough apart that pairs do not interact, and
    returns the final weight as a fraction of the initial one -- the
    quantity plotted in the experimental papers.
    """
    conn = PlasticConnections([0], [1], [w0], rule=rule, n_neurons=2)
    steps_between = int(pair_interval / dt)
    lag = int(round(abs(dt_ms) / dt))
    first, second = (0, 1) if dt_ms >= 0 else (1, 0)
    for _ in range(n_pairs):
        spike = np.zeros(2, dtype=bool)
        spike[first] = True
        conn.step(spike, dt)
        for _ in range(max(lag - 1, 0)):
            conn.step(np.zeros(2, dtype=bool), dt)
        spike = np.zeros(2, dtype=bool)
        spike[second] = True
        conn.step(spike, dt)
        for _ in range(steps_between):
            conn.step(np.zeros(2, dtype=bool), dt)
    return conn.weights[0] / w0
