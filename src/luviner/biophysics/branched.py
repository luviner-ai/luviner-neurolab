"""
Vectorized populations of branched neurons, and plasticity that is gated
by the branch rather than by the soma.

`dendrite.py` and `plasticity.py` have never been combined, and the gap
between them is where the interesting claim lives. The dendritic result
was negative *at the soma*: an active branch computes the same function
as a passive one with rescaled synapses, once you only look at whether
the cell spiked. The plasticity result was obtained on point neurons,
which have no "where". Put them together and the readout changes -- a
synapse that learns from *its own branch's* depolarisation rather than
from the somatic spike is reading a variable the point neuron does not
have.

`BranchedPopulation` is `MultiCompartmentNeuron` at network scale, in the
one geometry the claim needs: a star, `B` passive branches on a
spike-generating soma, held as (n_cells, 1 + B) arrays. The object-based
class stays the reference; the tests check the two agree.

Two things are deliberately not general here. Every cell has the same
number of branches, and every branch is one compartment. A cable of
several segments per branch is what `build_star(comps_per_branch=...)`
does and is not needed to ask whether learning gated locally partitions
better than learning gated at the soma.

References:
    Losonczy A, Makara JK, Magee JC (2008). "Compartmentalized dendritic
        plasticity and input feature storage in neurons." Nature 452:436.
    Kastellakis G, Silva AJ, Poirazi P (2016). "Linking memories across
        time via neuronal and dendritic overlaps in model neurons."
        Cell Rep 17:1491-1509.
    Poirazi P, Brannon T, Mel BW (2003). "Pyramidal neuron as two-layer
        neural network." Neuron 37:989-999.
"""

import numpy as np

from .hodgkin_huxley import _x_over_1_minus_exp_neg_x as _xo
from .plasticity import PlasticConnections


class BranchedPopulation:
    """`n` cells, each a soma plus `n_branches` passive dendrites.

    Compartment 0 of every cell is the soma and carries the
    spike-generating conductances; compartments 1..B are the branches and
    are passive by default, which is the arrangement that matters -- the
    branch gets its own local voltage while the output spike stays where
    it belongs.

    State is held as (n_cells, 1 + n_branches). Kinetics are Traub-Miles
    referred to V_T, the same Class I cortical cell as
    `NeuronPopulation.cortical`, so a branched network and a point network
    can be compared without changing the cell.
    """

    def __init__(self, n_cells, n_branches=2, g_axial=0.12,
                 g_Na=50.0, g_K=5.0, g_L=0.1,
                 E_Na=50.0, E_K=-90.0, E_L=-70.0, C_m=1.0,
                 V_T=-63.0, temperature=6.3, spike_threshold=0.0):
        self.n = int(n_cells)
        self.B = int(n_branches)
        self.C = self.B + 1
        shape = (self.n, self.C)
        self.g_axial = float(g_axial)
        self.V_T = V_T
        self.C_m = np.full(shape, float(C_m))
        self.g_L = np.full(shape, float(g_L))
        self.E_L = np.full(shape, float(E_L))
        self.E_Na = np.full(shape, float(E_Na))
        self.E_K = np.full(shape, float(E_K))
        # Only the soma spikes.
        self.g_Na = np.zeros(shape)
        self.g_K = np.zeros(shape)
        self.g_Na[:, 0] = g_Na
        self.g_K[:, 0] = g_K
        self.phi = 3.0 ** ((float(temperature) - 6.3) / 10.0)
        self.spike_threshold = float(spike_threshold)

        self._pre, self._post, self._branch = [], [], []
        self._g_max, self._E_rev, self._tau_r, self._tau_d = [], [], [], []
        self._delay, self._mg, self._slices = [], [], []
        self._built = False
        self.reset()

    # ---- construction ----------------------------------------------------

    def connect(self, pre, post, branch, synapse, weights=None):
        """Wire pre -> (post, branch).

        `branch` is 0 for the soma and 1..B for a dendrite. It is the
        argument the point-neuron version does not have.
        """
        pre = np.atleast_1d(np.asarray(pre, dtype=int))
        post = np.atleast_1d(np.asarray(post, dtype=int))
        branch = np.atleast_1d(np.asarray(branch, dtype=int))
        pre, post, branch = np.broadcast_arrays(pre, post, branch)
        pre, post, branch = pre.ravel(), post.ravel(), branch.ravel()
        if pre.size == 0:
            return self
        if pre.min() < 0 or pre.max() >= self.n or post.min() < 0 \
                or post.max() >= self.n:
            raise IndexError("neuron index out of range")
        if branch.min() < 0 or branch.max() >= self.C:
            raise IndexError("branch index out of range")
        k = pre.size
        g = (np.full(k, synapse.g_max) if weights is None
             else np.broadcast_to(np.asarray(weights, float), (k,)).copy())
        start = sum(len(a) for a in self._g_max)
        self._slices.append(slice(start, start + k))
        self._pre.append(pre)
        self._post.append(post)
        self._branch.append(branch)
        self._g_max.append(g)
        self._E_rev.append(np.full(k, synapse.E_rev))
        self._tau_r.append(np.full(k, synapse.tau_rise))
        self._tau_d.append(np.full(k, synapse.tau_decay))
        self._delay.append(np.full(k, synapse.delay))
        self._mg.append(np.full(k, synapse.mg_conc if synapse.mg_block
                                else 0.0))
        self._built = False
        return self

    def connection_slice(self, index=-1):
        return self._slices[index]

    def _build(self, dt):
        cat = lambda p: (np.concatenate(p) if p else np.zeros(0, dtype=float))
        self.pre = (np.concatenate(self._pre) if self._pre
                    else np.zeros(0, dtype=int))
        self.post = (np.concatenate(self._post) if self._post
                     else np.zeros(0, dtype=int))
        self.branch = (np.concatenate(self._branch) if self._branch
                       else np.zeros(0, dtype=int))
        self.g_max = cat(self._g_max)
        self.E_rev = cat(self._E_rev)
        self.tau_r = cat(self._tau_r)
        self.tau_d = cat(self._tau_d)
        self.mg_conc = cat(self._mg)
        delay = cat(self._delay)
        self.n_syn = self.pre.size
        if self.n_syn:
            t_peak = ((self.tau_r * self.tau_d) / (self.tau_d - self.tau_r)
                      * np.log(self.tau_d / self.tau_r))
            self.norm = 1.0 / (np.exp(-t_peak / self.tau_d)
                               - np.exp(-t_peak / self.tau_r))
            self.decay_r = np.exp(-dt / self.tau_r)
            self.decay_d = np.exp(-dt / self.tau_d)
            self.has_mg = self.mg_conc > 0.0
            self.delay_steps = np.maximum(
                np.round(delay / dt).astype(int) - 1, 0)
            self._buf_len = int(self.delay_steps.max()) + 2
        else:
            self.norm = self.decay_r = self.decay_d = cat([])
            self.has_mg = np.zeros(0, dtype=bool)
            self.delay_steps = np.zeros(0, dtype=int)
            self._buf_len = 2
        self.A = np.zeros(self.n_syn)
        self.Bsyn = np.zeros(self.n_syn)
        self._spike_buf = np.zeros((self._buf_len, self.n), dtype=bool)
        self._ptr = 0
        self._dt_built = dt
        self._built = True

    # ---- state -----------------------------------------------------------

    def reset(self, V0=None):
        V0 = self.E_L[0, 0] if V0 is None else V0
        self.V = np.full((self.n, self.C), float(V0))
        self.m, self.h, self.nn = self._steady_state(self.V)
        self.t = 0.0
        self.spike_times = [[] for _ in range(self.n)]
        self.last_fired = np.zeros(self.n, dtype=bool)
        if self._built:
            self.A[:] = 0.0
            self.Bsyn[:] = 0.0
            self._spike_buf[:] = False
            self._ptr = 0
        return self

    def _rates(self, V):
        x = V - self.V_T
        return (1.28 * _xo((x - 13.0) / 4.0),
                1.4 * _xo(-(x - 40.0) / 5.0),
                0.128 * np.exp(-(x - 17.0) / 18.0),
                4.0 / (1.0 + np.exp(-(x - 40.0) / 5.0)),
                0.16 * _xo((x - 15.0) / 5.0),
                0.5 * np.exp(-(x - 10.0) / 40.0))

    def _steady_state(self, V):
        am, bm, ah, bh, an, bn = self._rates(V)
        return am / (am + bm), ah / (ah + bh), an / (an + bn)

    # ---- dynamics --------------------------------------------------------

    def axial_currents(self, V):
        """Star coupling: every branch exchanges with the soma only.

        Sign convention follows `MultiCompartmentNeuron.axial_currents` --
        the flow is out of the child and into the parent, and both appear
        in the derivative with a minus sign.
        """
        I = np.zeros_like(V)
        flow = self.g_axial * (V[:, 1:] - V[:, :1])
        I[:, 1:] = flow
        I[:, 0] = -flow.sum(axis=1)
        return I

    def synaptic_currents(self, V):
        I = np.zeros_like(V)
        if not self.n_syn:
            return I
        g = self.g_max * self.norm * (self.Bsyn - self.A)
        V_post = V[self.post, self.branch]
        if self.has_mg.any():
            unblock = 1.0 / (1.0 + np.exp(-0.062 * V_post)
                             * (self.mg_conc / 3.57))
            g = np.where(self.has_mg, g * unblock, g)
        np.add.at(I, (self.post, self.branch), g * (V_post - self.E_rev))
        return I

    def _derivatives(self, V, m, h, n, I_ext):
        I_Na = self.g_Na * m ** 3 * h * (V - self.E_Na)
        I_K = self.g_K * n ** 4 * (V - self.E_K)
        I_L = self.g_L * (V - self.E_L)
        I_ax = self.axial_currents(V)
        I_syn = self.synaptic_currents(V)
        am, bm, ah, bh, an, bn = self._rates(V)
        p = self.phi
        return ((I_ext - I_Na - I_K - I_L - I_syn - I_ax) / self.C_m,
                p * (am * (1.0 - m) - bm * m),
                p * (ah * (1.0 - h) - bh * h),
                p * (an * (1.0 - n) - bn * n))

    def step(self, I_ext=0.0, dt=0.01):
        if not self._built or self._dt_built != dt:
            self._build(dt)
        I = np.broadcast_to(np.asarray(I_ext, dtype=float), (self.n, self.C))
        base = (self.V, self.m, self.h, self.nn)

        def deriv(k=None, frac=0.0):
            st = base if k is None else tuple(
                s + frac * dt * kk for s, kk in zip(base, k))
            return self._derivatives(st[0], st[1], st[2], st[3], I)

        k1 = deriv()
        k2 = deriv(k1, 0.5)
        k3 = deriv(k2, 0.5)
        k4 = deriv(k3, 1.0)
        out = [s + (dt / 6.0) * (a + 2.0 * b + 2.0 * c + d)
               for s, a, b, c, d in zip(base, k1, k2, k3, k4)]
        prev_soma = self.V[:, 0]
        self.V = out[0]
        self.m = np.clip(out[1], 0.0, 1.0)
        self.h = np.clip(out[2], 0.0, 1.0)
        self.nn = np.clip(out[3], 0.0, 1.0)

        if self.n_syn:
            arrived = self._spike_buf[
                (self._ptr - self.delay_steps) % self._buf_len, self.pre]
            self.A = (self.A + arrived) * self.decay_r
            self.Bsyn = (self.Bsyn + arrived) * self.decay_d

        self.t += dt
        fired = (prev_soma <= self.spike_threshold) & \
                (self.V[:, 0] > self.spike_threshold)
        self.last_fired = fired
        self._ptr = (self._ptr + 1) % self._buf_len
        self._spike_buf[self._ptr] = fired
        if fired.any():
            for i in np.flatnonzero(fired):
                self.spike_times[i].append(self.t)
        return self.V

    # ---- what the branch knows -------------------------------------------

    def branch_voltage(self):
        """(n_cells, n_branches) -- the local variable a soma cannot see."""
        return self.V[:, 1:]

    def branch_active(self, threshold=-40.0):
        """Which branches are depolarised past the plateau threshold.

        This is the gate. A synapse whose branch is above threshold is in
        a position to change; one whose branch is not, is not -- whatever
        the soma did.
        """
        return self.V[:, 1:] >= float(threshold)


class BranchGatedConnections(PlasticConnections):
    """STDP whose updates are permitted by the postsynaptic *branch*.

    The pair rule is unchanged -- same traces, same window, same bounds.
    What changes is who is allowed to use it: a synapse may move only
    while its own branch is depolarised past the plateau threshold. A
    somatic spike no longer licenses every synapse on the cell, which is
    the whole difference between this and the point-neuron version.

    The mechanism the claim rests on: assembly A's synapses cluster on
    one branch and B's on another, and A's activity during B's training
    cannot touch B's synapses, because the gate that permits learning is
    local. Interference is removed by geometry rather than by inhibition.

    `plateau` is the branch voltage at which learning is permitted. It is
    a threshold on a real variable, not a free knob: the NMDA plateau in
    `dendrite.py` reaches 28-57 mV, and the branch sits near E_L = -70 mV
    when nothing is happening to it.
    """

    def __init__(self, pre, post, branch, weights, plateau=-40.0,
                 n_branches=None, gate_homeostasis=None, **kw):
        super().__init__(pre, post, weights, **kw)
        self.branch = np.asarray(branch, dtype=int)
        if self.branch.size != self.pre.size:
            raise ValueError("branch must have one entry per connection")
        self.gate_homeostasis = gate_homeostasis
        if gate_homeostasis is not None:
            if n_branches is None:
                raise ValueError("regulating the gate needs n_branches")
            # One threshold per branch, so each can find its own set point.
            self.plateau = np.full((self.n_neurons, int(n_branches)),
                                   float(plateau))
            gate_homeostasis.reset(self.plateau.shape)
        else:
            self.plateau = float(plateau)

    def _threshold(self):
        """The threshold as a per-connection vector."""
        if np.isscalar(self.plateau):
            return self.plateau
        return self.plateau[self.post, self.branch - 1]

    def gate(self, V):
        """Per-connection: is this synapse's own branch depolarised?"""
        return V[self.post, self.branch] >= self._threshold()

    def step(self, fired, V=None, dt=0.05, learn=True):
        """As `PlasticConnections.step`, with the branch gate applied.

        `V` is the (n_cells, 1 + n_branches) voltage array. Passing None
        removes the gate, which makes this class fall back exactly on its
        parent and gives the ungated control for free.
        """
        fired = np.asarray(fired, dtype=bool)
        if learn and len(self.weights):
            g = np.ones(len(self.weights), dtype=bool) if V is None \
                else self.gate(V)
            pre_fired = fired[self.pre] & g
            post_fired = fired[self.post] & g
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

        if self.homeostasis is not None:
            self.weights = self.homeostasis.step(self.weights, self.post,
                                                 fired, dt)
        if self.gate_homeostasis is not None and V is not None:
            self.plateau = self.gate_homeostasis.step(self.plateau,
                                                      V[:, 1:], dt)
        return self

    def gated_fraction(self, V):
        """How much of the time the gate is open, for diagnosis.

        A gate that is always open is not gating, and a gate that is never
        open is not learning. Either reads as a null for reasons that have
        nothing to do with the claim, so it is measured rather than
        assumed.
        """
        return float(self.gate(V).mean()) if len(self.weights) else 0.0


class BranchNormalization:
    """A conductance budget per *branch* rather than per cell.

    `WeightNormalization` holds the total conductance arriving at a
    neuron. That is the right constraint when the neuron is a point, and
    the wrong one as soon as learning is gated locally: the gate opens on
    a plateau, potentiation raises the branch's conductance, the plateau
    latches, and the gate never closes again. Measured in
    experiments/branches/stream_separation.py, that runaway takes one
    block -- segregation falls from 0.021 to zero and stays there, with
    the gate open on every step from the first.

    A gate that cannot close is not a gate, and holding the budget per
    branch does close it -- measured, the open fraction falls from 1.000
    to 0.833.

    **It is the wrong constraint for making branches specialise, and the
    reason is worth keeping.** Normalising within a branch drives every
    branch to the same total, so two branches receiving the same input
    end up identical: the homeostasis erases the very asymmetry that
    clustering consists of. For specialisation the budget belongs to the
    *cell* (`WeightNormalization`), which lets one branch grow at
    another's expense, while the branch gate decides which one. Use this
    class when the runaway is the problem, not when the segregation is.

        w -> w * (budget / that branch's incoming total)

    applied with time constant `tau`, so it is a process the cell runs
    rather than a projection imposed each step.
    """

    def __init__(self, budget=0.12, tau=100.0, w_min=0.0, w_max=np.inf):
        if budget <= 0:
            raise ValueError("budget must be positive")
        self.budget = float(budget)
        self.tau = float(tau)
        self.w_min = w_min
        self.w_max = w_max
        self.branch = None
        self.n_branches = None

    def bind(self, branch, n_branches):
        """Told which branch each connection lands on, and how many there
        are. `PlasticConnections` only knows `post`, so the compartment
        has to be supplied from outside."""
        self.branch = np.asarray(branch, dtype=int)
        self.n_branches = int(n_branches)
        return self

    def reset(self, n_neurons):
        self.n_neurons = int(n_neurons)
        return self

    def _key(self, post):
        return post * (self.n_branches + 1) + self.branch

    def step(self, weights, post, fired, dt):
        if not len(weights) or self.branch is None:
            return weights
        size = self.n_neurons * (self.n_branches + 1)
        total = np.zeros(size)
        key = self._key(post)
        np.add.at(total, key, weights)
        factor = np.where(total > 0, self.budget / np.maximum(total, 1e-12),
                          1.0)
        return np.clip(weights * (1.0 + (dt / self.tau) * (factor[key] - 1.0)),
                       self.w_min, self.w_max)


class GateHomeostasis:
    """Regulates the plateau threshold so a branch's gate keeps gating.

    The gate has two useless states and they are the ones it reaches on
    its own. Shut, nothing learns. Latched open, every synapse on the
    cell is licensed at once and the gate carries no information about
    which branch did what -- and because a plateau is regenerative,
    potentiation drives it to latch.

    Measured in experiments/branches/stream_separation.py: across 40 runs
    in which the gate opened at all, the duty cycle ran 0.50 to 1.00, and
    the single run that separated its two input streams is the one with
    the lowest duty cycle. That is one event, and one event out of forty
    is a lead rather than a result -- but it is a lead with a mechanism
    and a cheap test, which is to hold the duty cycle where it is
    informative and see whether separation becomes reliable.

    So: give the gate its own set point. Each branch tracks the fraction
    of time it spends above threshold and moves its threshold to hold
    that fraction at `target`, slowly enough to be a process the cell
    runs rather than a constraint imposed each step. Raising the
    threshold when the branch is open too often is the direct analogue of
    what `WeightNormalization` does for total conductance, applied to the
    variable that decides *when learning is allowed* instead of to the
    weights themselves.
    """

    def __init__(self, target=0.5, tau=200.0, gain=40.0, tau_trace=100.0,
                 lo=-70.0, hi=-10.0):
        if not 0.0 < target < 1.0:
            raise ValueError("target must be a fraction strictly inside (0,1)")
        self.target = float(target)
        self.tau = float(tau)
        self.gain = float(gain)
        self.tau_trace = float(tau_trace)
        self.lo, self.hi = float(lo), float(hi)
        self.open_fraction = None

    def reset(self, shape):
        """`shape` is (n_cells, n_branches)."""
        self.open_fraction = np.full(shape, self.target)
        return self

    def step(self, plateau, V_branch, dt):
        """Advance the trace and return the updated per-branch threshold.

        `plateau` and `V_branch` are both (n_cells, n_branches).
        """
        if self.open_fraction is None:
            self.reset(V_branch.shape)
        is_open = (V_branch >= plateau).astype(float)
        a = np.exp(-dt / self.tau_trace)
        self.open_fraction = a * self.open_fraction + (1.0 - a) * is_open
        err = self.open_fraction - self.target
        return np.clip(plateau + (dt / self.tau) * self.gain * err,
                       self.lo, self.hi)
