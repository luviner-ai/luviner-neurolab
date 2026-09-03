"""
Multi-compartment neurons: computation in the dendrites.

Every neuron so far has been a point: one voltage, one set of channels,
inputs summing into a single node. That is a useful abstraction and it is
what most network models use, but it discards the part of the cell where
most of the machinery is. Around 90% of a pyramidal neuron's membrane is
dendrite, and synapses land there, not on the soma.

The difference is not cosmetic. A point neuron sums its inputs and then
applies one nonlinearity. A neuron with dendrites applies a nonlinearity
*in each branch* and then sums the results, which is a two-layer network
rather than a one-layer one. Two inputs arriving on the same branch
interact; the same two inputs on different branches do not. That is a
computation a point neuron cannot express at all, and it is why a single
pyramidal cell has been argued to be equivalent to a small network
(Poirazi, Brannon & Mel 2003).

Compartments are coupled by axial conductance, from cable theory:

    I_axial = g_axial * (V_neighbour - V_this)

with g_axial set by the geometry of the connecting segment. Each
compartment carries its own channels, so a dendritic branch can have NMDA
receptors and voltage-gated calcium while the soma has the sodium and
potassium that make the output spike.

The tree is stored as a parent array -- compartment i's parent is
`parent[i]`, with the soma at index 0 -- which is enough to describe any
branching structure without cycles, and lets the axial currents be
computed in one vectorized pass.

Compartments may also carry any channel from `channels.py`. That is what
turns a dendrite from a passive cable into an active one: a branch with
its own Ca2+ or K+ conductance regenerates or shapes the signal locally
instead of merely attenuating it, so the branch nonlinearity becomes a
spike rather than a bend. Each channel keeps its own state per
compartment, and intracellular calcium is tracked per compartment,
because a local dendritic Ca2+ transient is local.

References:
    Rall W (1959). "Branching dendritic trees and motoneuron membrane
        resistivity." Exp Neurol 1:491-527.
    Poirazi P, Brannon T, Mel BW (2003). "Pyramidal neuron as two-layer
        neural network." Neuron 37:989-999.
    Schiller J, Major G, Koester HJ, Schiller Y (2000). "NMDA spikes in
        basal dendrites of cortical pyramidal neurons."
        Nature 404:285-289.
"""

import numpy as np

from .hodgkin_huxley import HodgkinHuxleyNeuron as _H


class Compartment:
    """One electrical segment of a neuron.

    Parameters
    ----------
    parent : index of the compartment this one attaches to; None for the
        soma.
    g_axial : coupling conductance to the parent, mS/cm^2. Larger means
        more electrically continuous; smaller means more isolated, which
        is what lets a thin branch behave as its own computational unit.
    g_Na, g_K : spike-generating conductances. Zero makes the compartment
        passive, which is the right default for a distal dendrite.
    channels : IonChannel instances from `channels.py` carried by this
        compartment. A channel object holds parameters only -- the state
        lives in the neuron -- so one instance may be shared between
        compartments without them interfering.
    """

    def __init__(self, parent=None, g_axial=1.0, g_Na=0.0, g_K=0.0,
                 g_L=0.1, E_Na=50.0, E_K=-90.0, E_L=-70.0, C_m=1.0,
                 label=None, channels=()):
        self.parent = parent
        self.g_axial = g_axial
        self.g_Na, self.g_K, self.g_L = g_Na, g_K, g_L
        self.E_Na, self.E_K, self.E_L = E_Na, E_K, E_L
        self.C_m = C_m
        self.channels = list(channels)
        self.label = label or ('soma' if parent is None else 'dendrite')


class MultiCompartmentNeuron:
    """A neuron built from coupled compartments.

    Synaptic input is delivered per-compartment, which is the point: where
    an input lands changes what it does.
    """

    def __init__(self, compartments, V_T=-63.0, temperature=6.3,
                 calcium=None):
        self.comps = list(compartments)
        self.n = len(self.comps)
        if self.comps[0].parent is not None:
            raise ValueError("compartment 0 must be the soma (parent=None)")
        for i, c in enumerate(self.comps[1:], 1):
            if c.parent is None or not (0 <= c.parent < i):
                raise ValueError(
                    f"compartment {i} must attach to an earlier one")
        arr = lambda f: np.array([getattr(c, f) for c in self.comps], float)
        self.g_Na, self.g_K, self.g_L = arr('g_Na'), arr('g_K'), arr('g_L')
        self.E_Na, self.E_K, self.E_L = arr('E_Na'), arr('E_K'), arr('E_L')
        self.C_m = arr('C_m')
        self.g_axial = arr('g_axial')
        self.parent = np.array([-1 if c.parent is None else c.parent
                                for c in self.comps], int)
        self.V_T = V_T
        self.phi = 3.0 ** ((temperature - 6.3) / 10.0)

        # Extra channels, flattened to (compartment, channel) pairs. Their
        # gating state is held in one flat vector so the RK4 step treats it
        # exactly like V, m, h and n.
        self.channels = [(i, c) for i, comp in enumerate(self.comps)
                         for c in comp.channels]
        self._n_extra = sum(len(c.state_names) for _, c in self.channels)
        self._track_ca = any(c.needs_calcium or c.carries_calcium
                             for _, c in self.channels)
        if self._track_ca:
            from .channels import CalciumDynamics
            self.calcium = calcium or CalciumDynamics()
        else:
            self.calcium = calcium

        self.synapses = []          # (compartment, Synapse-like params)
        self.reset()

    # ---- kinetics (Traub-Miles, referred to V_T) -------------------------

    def _rates(self, V):
        from .hodgkin_huxley import _x_over_1_minus_exp_neg_x as xo
        x = V - self.V_T
        return (1.28 * xo((x - 13.0) / 4.0),
                1.4 * xo(-(x - 40.0) / 5.0),
                0.128 * np.exp(-(x - 17.0) / 18.0),
                4.0 / (1.0 + np.exp(-(x - 40.0) / 5.0)),
                0.16 * xo((x - 15.0) / 5.0),
                0.5 * np.exp(-(x - 10.0) / 40.0))

    def reset(self, V0=None):
        V0 = self.E_L.copy() if V0 is None else np.broadcast_to(
            np.asarray(V0, float), (self.n,)).astype(float).copy()
        self.V = V0
        am, bm, ah, bh, an, bn = self._rates(self.V)
        self.m, self.h, self.nn = am / (am + bm), ah / (ah + bh), an / (an + bn)
        self.extra = np.array(
            [float(v) for i, c in self.channels
             for v in c.steady_state(self.V[i])])
        self.Ca = (np.full(self.n, self.calcium.Ca_rest)
                   if self._track_ca else None)
        self.t = 0.0
        for syn in self.synapses:
            syn['A'] = 0.0
            syn['B'] = 0.0
        return self

    # ---- synapses --------------------------------------------------------

    def add_synapse(self, compartment, synapse):
        """Attach a synapse to a named compartment.

        Location matters: the same synapse on a thin distal branch and on
        the soma produce different computations, which is the whole point
        of having compartments.
        """
        if not (0 <= compartment < self.n):
            raise IndexError("no such compartment")
        t_peak = ((synapse.tau_rise * synapse.tau_decay)
                  / (synapse.tau_decay - synapse.tau_rise)
                  * np.log(synapse.tau_decay / synapse.tau_rise))
        norm = 1.0 / (np.exp(-t_peak / synapse.tau_decay)
                      - np.exp(-t_peak / synapse.tau_rise))
        rec = {'comp': compartment, 'g_max': synapse.g_max,
               'E_rev': synapse.E_rev, 'tau_r': synapse.tau_rise,
               'tau_d': synapse.tau_decay, 'norm': norm,
               'mg': synapse.mg_conc if synapse.mg_block else 0.0,
               'A': 0.0, 'B': 0.0}
        self.synapses.append(rec)
        return len(self.synapses) - 1

    def activate(self, index, n_quanta=1.0):
        """Deliver a presynaptic spike to synapse `index`."""
        self.synapses[index]['A'] += n_quanta
        self.synapses[index]['B'] += n_quanta
        return self

    def synaptic_currents(self, V):
        I = np.zeros(self.n)
        for s in self.synapses:
            g = s['g_max'] * s['norm'] * (s['B'] - s['A'])
            if s['mg'] > 0:
                Vp = V[s['comp']]
                g = g / (1.0 + np.exp(-0.062 * Vp) * (s['mg'] / 3.57))
            I[s['comp']] += g * (V[s['comp']] - s['E_rev'])
        return I

    # ---- extra channels --------------------------------------------------

    def _split_extra(self, extra):
        """Slice the flat gating vector into one tuple per channel."""
        out, i = [], 0
        for _, c in self.channels:
            k = len(c.state_names)
            out.append(tuple(extra[i:i + k]))
            i += k
        return out

    def channel_currents(self, V=None, extra=None, Ca=None):
        """Per-compartment current from each extra channel.

        Returned as a list of (compartment, name, current), because the
        same channel type may sit on several branches and the whole point
        is which branch it is on.
        """
        V = self.V if V is None else V
        extra = self.extra if extra is None else extra
        Ca = self.Ca if Ca is None else Ca
        return [(i, c.name,
                 float(c.current(V[i], st, None if Ca is None else Ca[i])))
                for (i, c), st in zip(self.channels, self._split_extra(extra))]

    def _extra_currents(self, V, extra, Ca):
        """Total extra-channel current per compartment, its calcium part,
        and the gating derivatives, in one pass."""
        I_extra = np.zeros(self.n)
        I_Ca = np.zeros(self.n)
        d_extra = []
        for (i, c), st in zip(self.channels, self._split_extra(extra)):
            Ca_i = None if Ca is None else Ca[i]
            I_c = c.current(V[i], st, Ca_i)
            I_extra[i] += I_c
            if c.carries_calcium:
                I_Ca[i] += I_c
            d_extra.extend(c.derivatives(V[i], st, Ca_i))
        return I_extra, I_Ca, np.array(d_extra, dtype=float)

    # ---- axial coupling --------------------------------------------------

    def axial_currents(self, V):
        """Current flowing between compartments, uA/cm^2.

        Each compartment exchanges with its parent; by Kirchhoff the
        parent receives the opposite of what the child sends, so one pass
        over the children accounts for both directions.
        """
        I = np.zeros(self.n)
        child = np.arange(1, self.n)
        par = self.parent[1:]
        flow = self.g_axial[1:] * (V[child] - V[par])
        np.add.at(I, par, -flow)     # into the parent
        np.add.at(I, child, flow)    # out of the child
        return I

    # ---- dynamics --------------------------------------------------------

    def _derivatives(self, state, I_ext):
        V, m, h, n, extra = state[0], state[1], state[2], state[3], state[4]
        Ca = state[5] if self._track_ca else None
        I_Na = self.g_Na * m ** 3 * h * (V - self.E_Na)
        I_K = self.g_K * n ** 4 * (V - self.E_K)
        I_L = self.g_L * (V - self.E_L)
        I_ax = self.axial_currents(V)
        I_syn = self.synaptic_currents(V)
        I_extra, I_Ca, d_extra = self._extra_currents(V, extra, Ca)
        am, bm, ah, bh, an, bn = self._rates(V)
        p = self.phi
        out = [(I_ext - I_Na - I_K - I_L - I_syn - I_ax - I_extra) / self.C_m,
               p * (am * (1 - m) - bm * m),
               p * (ah * (1 - h) - bh * h),
               p * (an * (1 - n) - bn * n),
               d_extra]
        if self._track_ca:
            out.append(self.calcium.derivative(Ca, I_Ca))
        return out

    def step(self, I_ext=0.0, dt=0.01):
        I_ext = np.broadcast_to(np.asarray(I_ext, float), (self.n,))
        s = [self.V, self.m, self.h, self.nn, self.extra]
        if self._track_ca:
            s.append(self.Ca)
        k1 = self._derivatives(s, I_ext)
        k2 = self._derivatives([a + 0.5 * dt * b for a, b in zip(s, k1)], I_ext)
        k3 = self._derivatives([a + 0.5 * dt * b for a, b in zip(s, k2)], I_ext)
        k4 = self._derivatives([a + dt * b for a, b in zip(s, k3)], I_ext)
        out = [a + (dt / 6.0) * (b + 2 * c + 2 * d + e)
               for a, b, c, d, e in zip(s, k1, k2, k3, k4)]
        self.V = out[0]
        self.m = np.clip(out[1], 0.0, 1.0)
        self.h = np.clip(out[2], 0.0, 1.0)
        self.nn = np.clip(out[3], 0.0, 1.0)
        self.extra = np.clip(out[4], 0.0, 1.0)
        if self._track_ca:
            self.Ca = np.maximum(out[5], 0.0)
        for syn in self.synapses:
            syn['A'] *= np.exp(-dt / syn['tau_r'])
            syn['B'] *= np.exp(-dt / syn['tau_d'])
        self.t += dt
        return self.V

    def simulate(self, duration, dt=0.01, I_ext=None, events=()):
        """Run the cell. `events` is a list of (time_ms, synapse_index)."""
        n_steps = int(round(duration / dt))
        queue = sorted(events, key=lambda e: e[0])
        qi = 0
        V = np.empty((n_steps, self.n))
        for i in range(n_steps):
            t = i * dt
            while qi < len(queue) and queue[qi][0] <= t:
                self.activate(queue[qi][1])
                qi += 1
            V[i] = self.V
            self.step(0.0 if I_ext is None else I_ext, dt)
        return {'t': np.arange(n_steps) * dt, 'V': V, 'dt': dt}


def build_star(n_branches=2, comps_per_branch=1, g_axial=0.6,
               g_Na_soma=50.0, g_K_soma=5.0, g_L=0.1, dendrite_active=False,
               dendrite_channels=()):
    """A soma with `n_branches` dendrites attached to it.

    The simplest structure that can show the thing that matters: two
    inputs on one branch interact, the same two inputs on different
    branches do not.

    `dendrite_channels` puts extra conductances from `channels.py` on
    every dendritic compartment and none on the soma, which is the
    arrangement that matters: the branch gets its own active machinery
    while the output spike stays where it belongs.
    """
    comps = [Compartment(None, g_axial=0.0, g_Na=g_Na_soma, g_K=g_K_soma,
                         g_L=g_L, label='soma')]
    for b in range(n_branches):
        parent = 0
        for k in range(comps_per_branch):
            comps.append(Compartment(
                parent, g_axial=g_axial,
                g_Na=(15.0 if dendrite_active else 0.0),
                g_K=(3.0 if dendrite_active else 0.0),
                g_L=g_L, channels=dendrite_channels,
                label=f'branch{b}_seg{k}'))
            parent = len(comps) - 1
    return MultiCompartmentNeuron(comps)
