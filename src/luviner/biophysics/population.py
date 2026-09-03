"""
Vectorized populations of Hodgkin-Huxley neurons.

`Network` in network.py holds one Python object per neuron and steps them
in a loop. That is readable and it is what the test suite validates
against, but it pays interpreter overhead per neuron per timestep, so it
stops being usable somewhere around a hundred cells.

`NeuronPopulation` holds the same state as flat arrays and advances every
neuron and every synapse with whole-array operations. The equations are
identical -- the test suite checks the two implementations agree to
floating-point tolerance -- but the cost per step no longer scales with
interpreter work.

Synapses are stored per-connection rather than as a dense matrix, so a
sparsely wired network costs what its connections cost, not N^2.
Axonal delays use a circular spike buffer, so an arbitrary per-connection
delay costs one gather.
"""

import numpy as np

from .hodgkin_huxley import (HodgkinHuxleyNeuron,
                              _x_over_1_minus_exp_neg_x as _xo)


class NeuronPopulation:
    """N Hodgkin-Huxley neurons plus their synapses, held as arrays.

    Per-neuron parameters may be given as scalars (shared) or as arrays
    of length n_neurons (heterogeneous populations).
    """

    # Cortical defaults (Traub-Miles kinetics, Pospischil conductances).
    CORTICAL = dict(g_Na=50.0, g_K=5.0, g_L=0.1,
                    E_Na=50.0, E_K=-90.0, E_L=-70.0)

    @classmethod
    def cortical(cls, n_neurons, **kw):
        """Build a population of Class I cortical cells.

        Same cell as CorticalNeuron: fires continuously from near zero
        rather than jumping to a ~55 Hz floor, so it can adapt instead of
        switching off. Extra keywords override the cortical defaults.
        """
        params = dict(cls.CORTICAL)
        params.update(kw)
        return cls(n_neurons, kinetics='cortical', **params)

    def __init__(self, n_neurons, C_m=1.0, g_Na=120.0, g_K=36.0, g_L=0.3,
                 E_Na=50.0, E_K=-77.0, E_L=-54.387, temperature=6.3,
                 spike_threshold=0.0, kinetics='squid', V_T=-63.0, g_M=0.0):
        if kinetics not in ('squid', 'cortical'):
            raise ValueError("kinetics must be 'squid' or 'cortical'")
        self.n = int(n_neurons)
        self.kinetics = kinetics
        b = lambda x: np.broadcast_to(np.asarray(x, dtype=float), (self.n,)).copy()
        self.V_T = b(V_T)
        # Per-neuron muscarinic K+ conductance; zero disables it.
        self.g_M = b(g_M)
        self._has_M = bool(np.any(self.g_M > 0.0))
        self.C_m = b(C_m)
        self.g_Na = b(g_Na)
        self.g_K = b(g_K)
        self.g_L = b(g_L)
        self.E_Na = b(E_Na)
        self.E_K = b(E_K)
        self.E_L = b(E_L)
        self.phi = 3.0 ** ((b(temperature) - 6.3) / 10.0)
        self.spike_threshold = spike_threshold

        # Per-connection synapse arrays, grown by connect().
        self._pre = []
        self._post = []
        self._g_max = []
        self._E_rev = []
        self._tau_r = []
        self._tau_d = []
        self._delay = []
        self._mg = []
        self._slices = []
        self._built = False
        # Optional short-term plasticity over the synapse array. None is
        # the model this package has always simulated.
        self.stp = None

        self.reset()

    # ---- construction ----------------------------------------------------

    def connect(self, pre, post, synapse, weights=None):
        """Wire pre -> post using a Synapse instance as the parameter set.

        The Synapse object is read for its parameters only; its own state
        is not used, so one prototype can be reused for many connections.

        `weights` gives a per-connection peak conductance, overriding the
        prototype's g_max. Structured connectivity -- a distance-dependent
        profile, a learned weight matrix -- is expressed this way.
        """
        pre = np.atleast_1d(np.asarray(pre, dtype=int))
        post = np.atleast_1d(np.asarray(post, dtype=int))
        pre, post = np.broadcast_arrays(pre, post)
        if pre.size == 0:
            return self
        if pre.min() < 0 or pre.max() >= self.n or post.min() < 0 or post.max() >= self.n:
            raise IndexError("neuron index out of range")
        k = pre.size
        if weights is None:
            g = np.full(k, synapse.g_max)
        else:
            g = np.broadcast_to(np.asarray(weights, dtype=float), (k,)).copy()
        self._pre.append(pre.ravel())
        self._post.append(post.ravel())
        self._g_max.append(g)
        # Record where this batch lands in the flat arrays, so a caller
        # holding a plastic rule can write updated weights back into it.
        start = sum(len(a) for a in self._g_max[:-1])
        self._slices.append(slice(start, start + k))
        self._E_rev.append(np.full(k, synapse.E_rev))
        self._tau_r.append(np.full(k, synapse.tau_rise))
        self._tau_d.append(np.full(k, synapse.tau_decay))
        self._delay.append(np.full(k, synapse.delay))
        self._mg.append(np.full(k, synapse.mg_conc if synapse.mg_block else 0.0))
        self._built = False
        return self

    def connection_slice(self, index=-1):
        """Where the synapses from one connect() call live in `g_max`.

        Lets an external learning rule write updated conductances back:
        `pop.g_max[pop.connection_slice(0)] = new_weights`.
        """
        return self._slices[index]

    def _build(self, dt):
        """Flatten the connection lists and precompute step constants."""
        cat = lambda parts: (np.concatenate(parts) if parts
                             else np.zeros(0, dtype=float))
        self.pre = (np.concatenate(self._pre) if self._pre
                    else np.zeros(0, dtype=int))
        self.post = (np.concatenate(self._post) if self._post
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
            # Circular buffer deep enough for the longest axon.
            # The buffer is read before the write pointer advances, so it
            # already carries one step of latency; subtract it so a spike
            # lands exactly `delay` after it was detected, matching the
            # reference implementation. Zero delay still costs one step,
            # since a spike cannot act on the step that detected it.
            self.delay_steps = np.maximum(
                np.round(delay / dt).astype(int) - 1, 0)
            self._buf_len = int(self.delay_steps.max()) + 2
        else:
            self.norm = self.decay_r = self.decay_d = cat([])
            self.has_mg = np.zeros(0, dtype=bool)
            self.delay_steps = np.zeros(0, dtype=int)
            self._buf_len = 2

        self.A = np.zeros(self.n_syn)
        self.B = np.zeros(self.n_syn)
        if self.stp is not None:
            self.stp.reset(self.n_syn, dt)
        self._spike_buf = np.zeros((self._buf_len, self.n), dtype=bool)
        self._ptr = 0
        self._dt_built = dt
        self._built = True

    # ---- state -----------------------------------------------------------

    def reset(self, V0=None):
        """Rest every neuron and clear all synaptic state."""
        if V0 is None:
            V0 = self.resting_potential()
        self.V = np.broadcast_to(np.asarray(V0, dtype=float), (self.n,)).copy()
        self.m, self.h, self.nn = self._steady_state(self.V)
        self.w = self._w_inf(self.V) if self._has_M else None
        self.t = 0.0
        self.spike_times = [[] for _ in range(self.n)]
        self.last_fired = np.zeros(self.n, dtype=bool)
        self._clamped = np.zeros(self.n, dtype=bool)
        if self._built:
            self.A[:] = 0.0
            self.B[:] = 0.0
            self._spike_buf[:] = False
            self._ptr = 0
        return self

    def resting_potential(self):
        """Per-neuron resting potential, solved by bisection."""
        lo = np.full(self.n, -90.0)
        hi = np.full(self.n, -50.0)
        f = lambda V: self._net_ionic(V)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            take_hi = f(lo) * f(mid) <= 0
            hi = np.where(take_hi, mid, hi)
            lo = np.where(take_hi, lo, mid)
        return 0.5 * (lo + hi)

    def _net_ionic(self, V):
        m, h, n = self._steady_state(V)
        I = (self.g_Na * m ** 3 * h * (V - self.E_Na)
             + self.g_K * n ** 4 * (V - self.E_K)
             + self.g_L * (V - self.E_L))
        if self._has_M:
            I = I + self.g_M * self._w_inf(V) * (V - self.E_K)
        return I

    def _rates(self, V):
        """Return (alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n).

        Dispatches on the kinetics family. Both are fully vectorized;
        the cortical set is written relative to the per-neuron threshold
        V_T, so a population may hold cells of differing excitability.
        """
        H = HodgkinHuxleyNeuron
        if self.kinetics == 'squid':
            return (H.alpha_m(V), H.beta_m(V), H.alpha_h(V),
                    H.beta_h(V), H.alpha_n(V), H.beta_n(V))
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

    @staticmethod
    def _w_inf(V):
        """I_M activation at equilibrium."""
        return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))

    @staticmethod
    def _tau_w(V, tau_max=1000.0):
        return tau_max / (3.3 * np.exp((V + 35.0) / 20.0)
                          + np.exp(-(V + 35.0) / 20.0))

    def voltage_clamp(self, idx, V_hold):
        """Hold the given neurons at V_hold; their gates keep evolving."""
        self._clamped[idx] = True
        self.V[idx] = V_hold
        return self

    def release_clamp(self, idx=None):
        self._clamped[:] = False if idx is None else self._clamped[idx]
        if idx is not None:
            self._clamped[idx] = False
        return self

    def block_sodium(self, idx, fraction=1.0):
        self.g_Na[idx] = 120.0 * (1.0 - fraction)
        return self

    # ---- dynamics --------------------------------------------------------

    def _derivatives(self, V, m, h, n, I_ext, w=None):
        I_Na = self.g_Na * m ** 3 * h * (V - self.E_Na)
        I_K = self.g_K * n ** 4 * (V - self.E_K)
        I_L = self.g_L * (V - self.E_L)
        I_M = (self.g_M * w * (V - self.E_K)) if self._has_M else 0.0
        am, bm, ah, bh, an, bn = self._rates(V)
        p = self.phi
        out = ((I_ext - I_Na - I_K - I_L - I_M) / self.C_m,
               p * (am * (1.0 - m) - bm * m),
               p * (ah * (1.0 - h) - bh * h),
               p * (an * (1.0 - n) - bn * n))
        if self._has_M:
            out = out + ((self._w_inf(V) - w) / self._tau_w(V),)
        return out

    def synaptic_currents(self):
        """Total synaptic current into each neuron, uA/cm^2."""
        I = np.zeros(self.n)
        if not self.n_syn:
            return I
        g = self.g_max * self.norm * (self.B - self.A)
        V_post = self.V[self.post]
        if self.has_mg.any():
            unblock = 1.0 / (1.0 + np.exp(-0.062 * V_post)
                             * (self.mg_conc / 3.57))
            g = np.where(self.has_mg, g * unblock, g)
        np.add.at(I, self.post, g * (V_post - self.E_rev))
        return I

    def step(self, I_ext=0.0, dt=0.01):
        """Advance the whole population one timestep (RK4)."""
        if not self._built or self._dt_built != dt:
            self._build(dt)

        I = np.broadcast_to(np.asarray(I_ext, dtype=float), (self.n,))
        I = I - self.synaptic_currents()

        base = (self.V, self.m, self.h, self.nn)
        if self._has_M:
            base = base + (self.w,)
        def deriv_at(k=None, frac=0.0):
            """Derivatives at the base state, or at a partial RK4 step."""
            st = base if k is None else tuple(
                s + frac * dt * kk for s, kk in zip(base, k))
            return self._derivatives(st[0], st[1], st[2], st[3], I,
                                     st[4] if self._has_M else None)

        k1 = deriv_at()
        k2 = deriv_at(k1, 0.5)
        k3 = deriv_at(k2, 0.5)
        k4 = deriv_at(k3, 1.0)
        out = [s + (dt / 6.0) * (a + 2.0 * b + 2.0 * c + d)
               for s, a, b, c, d in zip(base, k1, k2, k3, k4)]

        V_new = np.where(self._clamped, self.V, out[0])
        self.V = V_new
        self.m = np.clip(out[1], 0.0, 1.0)
        self.h = np.clip(out[2], 0.0, 1.0)
        self.nn = np.clip(out[3], 0.0, 1.0)
        if self._has_M:
            self.w = np.clip(out[4], 0.0, 1.0)

        # Synaptic conductances decay exactly; arriving spikes kick them.
        if self.n_syn:
            arrived = self._spike_buf[
                (self._ptr - self.delay_steps) % self._buf_len, self.pre]
            # With no short-term plasticity the kick is the boolean itself,
            # which is what it has always been; `stp` scales it by the
            # available resource instead. The None branch is left as an
            # exact identity so every existing result is unaffected.
            kick = arrived if self.stp is None else self.stp.release(arrived, dt)
            self.A = (self.A + kick) * self.decay_r
            self.B = (self.B + kick) * self.decay_d

        prev = base[0]        # membrane potential before this step
        self.t += dt
        fired = (prev <= self.spike_threshold) & (self.V > self.spike_threshold)
        # Exposed so a learning rule can be driven from the same step
        # without reconstructing it from spike_times.
        self.last_fired = fired
        self._ptr = (self._ptr + 1) % self._buf_len
        self._spike_buf[self._ptr] = fired
        if fired.any():
            for i in np.flatnonzero(fired):
                self.spike_times[i].append(self.t)
        return self.V

    def simulate(self, duration, dt=0.01, I_ext=None, record=True):
        """Run and optionally record every membrane potential."""
        n_steps = int(round(duration / dt))
        if I_ext is None:
            at = lambda i: 0.0
        elif callable(I_ext):
            at = I_ext
        else:
            arr = np.asarray(I_ext, dtype=float)
            at = (lambda i: arr[i]) if arr.ndim == 2 else (lambda i: arr)

        V = np.empty((n_steps, self.n)) if record else None
        for i in range(n_steps):
            if record:
                V[i] = self.V
            self.step(at(i), dt)
        return {'t': np.arange(n_steps) * dt, 'V': V, 'dt': dt,
                'spikes': [np.array(s) for s in self.spike_times]}
