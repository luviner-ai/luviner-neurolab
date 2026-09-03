"""
Networks integrated with per-neuron adaptive steps.

A single neuron benefits enormously from an adaptive step: it is
quiescent most of the time, so the step can grow between spikes. A
network appears not to. Whenever any cell spikes the step must shrink to
resolve it, and with N cells firing at 10 Hz the mean gap between spikes
*somewhere* in the network is 100/N ms. Past ten thousand neurons that
gap is smaller than the step a spike needs, and a globally adaptive
method degenerates into a fixed one. This is why large simulations --
NEST, and the cortical models run on Fugaku -- use a fixed step.

The way round it is the axonal delay. Every synapse has one, and if the
smallest in the network is `min_delay`, then no neuron can influence any
other in less than that. Each cell may therefore integrate freely for a
whole `min_delay` window, at whatever step size its own dynamics call
for, and only at the window boundary does anyone need to exchange
spikes.

Crucially the decoupling is exact rather than approximate. Within a
window a cell's synaptic conductances are a known analytic function of
time: they decay exponentially from their value at the window boundary,
and every spike that will arrive during the window was emitted at least
min_delay ago, so its arrival time is already fixed. The cell can
evaluate its synaptic current at any instant its integrator asks for,
including the intermediate stages of a Runge-Kutta step.

Spike times are recovered by interpolating the threshold crossing rather
than being rounded to a step boundary, so they stay accurate even where
the local step is large.

What this buys is measured in `stats`: the number of derivative
evaluations, which is the work that would have to be done in any
language. The Python wall-clock time is dominated by interpreter
overhead and is not the figure of merit.

The synchronisation-by-minimum-delay idea is standard in parallel
simulators, where it is used to avoid communication. Here it is used for
something else: to let each neuron keep its own step size.
"""

import numpy as np

from .adaptive import AdaptiveIntegrator
from .hodgkin_huxley import HodgkinHuxleyNeuron as _H


class _SynapseGroup:
    """Per-connection synaptic state with analytic within-window decay."""

    def __init__(self, pre, post, g_max, E_rev, tau_rise, tau_decay, delay,
                 mg_conc):
        self.pre = np.asarray(pre, dtype=int)
        self.post = np.asarray(post, dtype=int)
        self.g_max = np.asarray(g_max, dtype=float)
        self.E_rev = np.asarray(E_rev, dtype=float)
        self.tau_r = np.asarray(tau_rise, dtype=float)
        self.tau_d = np.asarray(tau_decay, dtype=float)
        self.delay = np.asarray(delay, dtype=float)
        self.mg_conc = np.asarray(mg_conc, dtype=float)
        self.has_mg = self.mg_conc > 0
        t_peak = ((self.tau_r * self.tau_d) / (self.tau_d - self.tau_r)
                  * np.log(self.tau_d / self.tau_r))
        self.norm = 1.0 / (np.exp(-t_peak / self.tau_d)
                           - np.exp(-t_peak / self.tau_r))
        self.A = np.zeros(len(self.pre))
        self.B = np.zeros(len(self.pre))


class _PopulationIntegrator(AdaptiveIntegrator):
    """Bogacki-Shampine over an (n, 4) population state.

    Only the error norm differs from the per-neuron integrator: each
    neuron's own RMS over its four variables, maximised over neurons, so
    an accepted step satisfies EVERY neuron's tolerance individually.
    The step arithmetic and the PI controller are inherited rather than
    re-implemented, so a comparison against the per-neuron path isolates
    vectorisation instead of confounding it with a second controller.
    """

    def error_norm(self, y, y3, y2):
        scale = self.atol + self.rtol * np.maximum(np.abs(y), np.abs(y3))
        per_neuron = np.sqrt(np.mean(((y3 - y2) / scale) ** 2, axis=1))
        return float(np.max(per_neuron))


class AdaptiveNetwork:
    """Network of Hodgkin-Huxley cells, each on its own adaptive step.

    Parameters
    ----------
    neurons : list of HodgkinHuxleyNeuron (or CorticalNeuron)
    min_delay : ms. Must be no larger than the smallest axonal delay in
        the network; it sets the synchronisation window.
    rtol, atol : per-neuron integration tolerances.
    """

    def __init__(self, neurons, min_delay=1.0, rtol=1e-4, atol=1e-6,
                 spike_threshold=0.0, dt_max=None):
        self.neurons = list(neurons)
        self.n = len(self.neurons)
        self.min_delay = float(min_delay)
        self.spike_threshold = spike_threshold
        self.rtol, self.atol = rtol, atol
        self.dt_max = dt_max if dt_max is not None else min_delay
        self._pre, self._post, self._params = [], [], []
        self.syn = None
        self.reset()

    # ---- construction ----------------------------------------------------

    def connect(self, pre, post, synapse):
        pre = np.atleast_1d(np.asarray(pre, dtype=int))
        post = np.atleast_1d(np.asarray(post, dtype=int))
        pre, post = np.broadcast_arrays(pre, post)
        if synapse.delay < self.min_delay - 1e-12:
            raise ValueError(
                f"synapse delay {synapse.delay} ms is shorter than "
                f"min_delay {self.min_delay} ms; the window would not be "
                "safe to decouple")
        k = pre.size
        self._pre.append(pre.ravel())
        self._post.append(post.ravel())
        self._params.append((np.full(k, synapse.g_max),
                             np.full(k, synapse.E_rev),
                             np.full(k, synapse.tau_rise),
                             np.full(k, synapse.tau_decay),
                             np.full(k, synapse.delay),
                             np.full(k, synapse.mg_conc if synapse.mg_block
                                     else 0.0)))
        self.syn = None
        return self

    def _build(self):
        if not self._pre:
            self.syn = _SynapseGroup([], [], [], [], [1.0], [2.0], [], [])
            return
        cols = list(zip(*self._params))
        self.syn = _SynapseGroup(np.concatenate(self._pre),
                                 np.concatenate(self._post),
                                 *[np.concatenate(c) for c in cols])
        # Connection indices grouped by postsynaptic cell, so a neuron's
        # own synaptic current can be evaluated without touching anyone
        # else's. See `_synaptic_current_for`.
        self._by_post = [np.flatnonzero(self.syn.post == i)
                         for i in range(self.n)]

    # ---- state -----------------------------------------------------------

    def reset(self):
        self.state = np.empty((self.n, 4))
        for i, cell in enumerate(self.neurons):
            Vr = cell.resting_potential()
            cell.reset(Vr)
            self.state[i] = (cell.V, cell.m, cell.h, cell.n)
        self.t = 0.0
        self.spike_times = [[] for _ in range(self.n)]
        self._pending = []            # (arrival_time, connection_index)
        if self.syn is not None:
            self.syn.A[:] = 0.0
            self.syn.B[:] = 0.0
        self.n_evaluations = 0
        self.n_vector_evaluations = 0
        self.n_rejected_shared = 0
        self.n_steps = 0
        self.dt_samples = []
        return self

    @property
    def stats(self):
        return {'evaluations': self.n_evaluations, 'steps': self.n_steps,
                'mean_dt': float(np.mean(self.dt_samples))
                if self.dt_samples else 0.0,
                'windows': int(round(self.t / self.min_delay))}

    # ---- within-window synaptic current ---------------------------------

    def _window_arrivals(self, t0, t1):
        """Split pending spike arrivals into those inside this window and
        those still in the future."""
        inside, future = [], []
        for arrival, idx in self._pending:
            (inside if arrival < t1 else future).append((arrival, idx))
        self._pending = future
        return inside

    def _synaptic_current_fn(self, t0, arrivals):
        """Build f(t, V_all) -> synaptic current per neuron, valid on the
        window starting at t0.

        Conductances decay analytically from their boundary values, plus
        one exponential per arrival that has already occurred by t. No
        state is advanced, so the integrator may evaluate at any t in the
        window, in any order, as Runge-Kutta stages require.
        """
        syn = self.syn
        if not len(syn.pre):
            return lambda t, V: np.zeros(self.n)
        by_conn = {}
        for arrival, idx in arrivals:
            by_conn.setdefault(idx, []).append(arrival)

        def current(t, V_all):
            dtau = t - t0
            A = syn.A * np.exp(-dtau / syn.tau_r)
            B = syn.B * np.exp(-dtau / syn.tau_d)
            for idx, times in by_conn.items():
                for a in times:
                    if a <= t:
                        A[idx] += np.exp(-(t - a) / syn.tau_r[idx])
                        B[idx] += np.exp(-(t - a) / syn.tau_d[idx])
            g = syn.g_max * syn.norm * (B - A)
            V_post = V_all[syn.post]
            if syn.has_mg.any():
                unblock = 1.0 / (1.0 + np.exp(-0.062 * V_post)
                                 * (syn.mg_conc / 3.57))
                g = np.where(syn.has_mg, g * unblock, g)
            I = np.zeros(self.n)
            np.add.at(I, syn.post, g * (V_post - syn.E_rev))
            return I
        return current

    def _arrivals_by_connection(self, arrivals):
        by_conn = {}
        for arrival, idx in arrivals:
            by_conn.setdefault(idx, []).append(arrival)
        return by_conn

    def _synaptic_current_for(self, i, t0, by_conn):
        """f(t, V_i) -> the synaptic current into neuron `i` alone.

        The same arithmetic as `_synaptic_current_fn`, restricted to the
        connections whose postsynaptic cell is `i`.

        **This is not an approximation.** Every term in that function
        depends on `V_post` and on nothing else, so for synapses onto `i`
        the only voltage that enters is `V_i` -- the whole-network vector
        was built, scattered into, and then indexed at `[i]` by the
        caller, with the other `n-1` entries discarded. Restricting the
        evaluation turns an `O(n_syn)` cost per neuron per Runge-Kutta
        stage into `O(fan-in)`, and removes the need for the frozen
        voltage snapshot entirely.
        """
        syn = self.syn
        if not len(syn.pre):
            return lambda t, V_i: 0.0
        idx = self._by_post[i]
        if idx.size == 0:
            return lambda t, V_i: 0.0
        A0, B0 = syn.A[idx], syn.B[idx]
        g_max, norm = syn.g_max[idx], syn.norm[idx]
        tau_r, tau_d = syn.tau_r[idx], syn.tau_d[idx]
        E_rev, has_mg, mg = syn.E_rev[idx], syn.has_mg[idx], syn.mg_conc[idx]
        any_mg = bool(has_mg.any())
        where = {int(c): int(j) for j, c in enumerate(idx)}
        local = [(where[c], by_conn[c]) for c in by_conn if c in where]

        def current(t, V_i):
            dtau = t - t0
            A = A0 * np.exp(-dtau / tau_r)
            B = B0 * np.exp(-dtau / tau_d)
            for j, times in local:
                for a in times:
                    if a <= t:
                        A[j] += np.exp(-(t - a) / tau_r[j])
                        B[j] += np.exp(-(t - a) / tau_d[j])
            g = g_max * norm * (B - A)
            if any_mg:
                unblock = 1.0 / (1.0 + np.exp(-0.062 * V_i) * (mg / 3.57))
                g = np.where(has_mg, g * unblock, g)
            return float(np.sum(g * (V_i - E_rev)))
        return current

    # ---- vectorised population path (T1e) --------------------------------
    #
    # The per-neuron path above pays one Python call and one small numpy
    # call per neuron per Runge-Kutta stage. T1c made that cost flat in
    # `n`; it did not make it small. This path integrates the whole
    # population as one (n, 4) state under a single shared step, which is
    # what `## Run: T1d` measured an upper bound for before it was
    # written.
    #
    # Nothing above this comment is modified: the existing `simulate` is
    # left exactly as it is so that it remains the reference the vector
    # path is checked against.

    def _cell_arrays(self):
        """Per-cell parameters gathered once, so cells may be heterogeneous."""
        c = self.neurons
        g = lambda name: np.array([getattr(x, name) for x in c], dtype=float)
        return {k: g(k) for k in ('C_m', 'g_Na', 'g_K', 'g_L',
                                  'E_Na', 'E_K', 'E_L', 'phi')}

    def _vector_synaptic_current(self, t0, arrivals):
        """f(t, V) -> synaptic current for every neuron at once.

        Same arithmetic as `_synaptic_current_for`, evaluated across all
        connections and scattered onto postsynaptic cells. The per-neuron
        version's Python loop over arrivals is replaced by a masked
        scatter, which is the whole point of the path.
        """
        syn = self.syn
        n = self.n
        if not len(syn.pre):
            zero = np.zeros(n)
            return lambda t, V: zero
        post = syn.post
        A0, B0 = syn.A.copy(), syn.B.copy()
        g_max, norm = syn.g_max, syn.norm
        tau_r, tau_d = syn.tau_r, syn.tau_d
        E_rev, has_mg, mg = syn.E_rev, syn.has_mg, syn.mg_conc
        any_mg = bool(has_mg.any())
        if arrivals:
            arr_t = np.array([a for a, _ in arrivals], dtype=float)
            arr_i = np.array([i for _, i in arrivals], dtype=int)
        else:
            arr_t = np.empty(0)
            arr_i = np.empty(0, dtype=int)

        def current(t, V):
            dtau = t - t0
            A = A0 * np.exp(-dtau / tau_r)
            B = B0 * np.exp(-dtau / tau_d)
            if arr_t.size:
                m = arr_t <= t
                if m.any():
                    ii = arr_i[m]
                    d = t - arr_t[m]
                    np.add.at(A, ii, np.exp(-d / tau_r[ii]))
                    np.add.at(B, ii, np.exp(-d / tau_d[ii]))
            g = g_max * norm * (B - A)
            Vp = V[post]
            if any_mg:
                unblock = 1.0 / (1.0 + np.exp(-0.062 * Vp) * (mg / 3.57))
                g = np.where(has_mg, g * unblock, g)
            I = np.zeros(n)
            np.add.at(I, post, g * (Vp - E_rev))
            return I
        return current

    def _vector_derivatives(self, syn_current, I_ext, cells=None):
        """f(t, Y) for the whole population; Y is (n, 4) as [V, m, h, n]."""
        c = cells if cells is not None else self._cell_arrays()
        H = self.neurons[0]
        C_m, phi = c['C_m'], c['phi']
        g_Na, g_K, g_L = c['g_Na'], c['g_K'], c['g_L']
        E_Na, E_K, E_L = c['E_Na'], c['E_K'], c['E_L']

        def f(t, Y):
            V, m, h, nn = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3]
            self.n_vector_evaluations += 1
            I_syn = syn_current(t, V)
            I_Na = g_Na * m ** 3 * h * (V - E_Na)
            I_K = g_K * nn ** 4 * (V - E_K)
            I_L = g_L * (V - E_L)
            out = np.empty_like(Y)
            out[:, 0] = (I_ext - I_syn - I_Na - I_K - I_L) / C_m
            out[:, 1] = phi * (H.alpha_m(V) * (1.0 - m) - H.beta_m(V) * m)
            out[:, 2] = phi * (H.alpha_h(V) * (1.0 - h) - H.beta_h(V) * h)
            out[:, 3] = phi * (H.alpha_n(V) * (1.0 - nn) - H.beta_n(V) * nn)
            return out
        return f

    def simulate_shared(self, duration, I_ext=None, dt0=0.01):
        """Integrate the population under one step shared by all neurons.

        The step is governed by the stiffest neuron at that moment: the
        error norm is each neuron's own RMS over its four state
        variables, maximised over neurons. **Not** an RMS over all 4n
        entries, which would divide one spiking cell's error by sqrt(n)
        and buy a larger step by quietly loosening the tolerance every
        other measurement on this board was taken at.
        """
        if self.syn is None:
            self._build()
        I_ext = (np.zeros(self.n) if I_ext is None
                 else np.broadcast_to(np.asarray(I_ext, dtype=float),
                                      (self.n,)).astype(float))
        cells = self._cell_arrays()
        integ = _PopulationIntegrator(rtol=self.rtol, atol=self.atol,
                                      dt_max=self.dt_max)
        n_windows = int(np.ceil(duration / self.min_delay))

        for w in range(n_windows):
            t0 = self.t
            t1 = min(t0 + self.min_delay, duration)
            if t1 <= t0:
                break
            arrivals = self._window_arrivals(t0, t1)
            f = self._vector_derivatives(
                self._vector_synaptic_current(t0, arrivals), I_ext, cells)
            resume = integ.last_dt
            ts, ys = integ.integrate(f, self.state, (t0, t1),
                                     dt0=(resume if resume else dt0))
            fired = []
            for i in range(self.n):
                V = ys[:, i, 0]
                above = V > self.spike_threshold
                for c in np.flatnonzero(~above[:-1] & above[1:]):
                    v0, v1 = V[c], V[c + 1]
                    frac = (self.spike_threshold - v0) / (v1 - v0)
                    t_spike = float(ts[c] + frac * (ts[c + 1] - ts[c]))
                    self.spike_times[i].append(t_spike)
                    fired.append((i, t_spike))
            self.state = ys[-1]
            self.n_steps += integ.stats['accepted']
            self.n_rejected_shared += integ.n_rejected
            self.dt_samples.extend(integ.dt_history)
            carry = integ.last_dt
            integ.reset_stats()
            integ.last_dt = carry

            self._advance_synapses(t0, t1, arrivals)
            if len(self.syn.pre):
                for i, t_spike in fired:
                    for idx in np.flatnonzero(self.syn.pre == i):
                        self._pending.append(
                            (t_spike + self.syn.delay[idx], int(idx)))
            self.t = t1

        return {'spikes': [np.array(s) for s in self.spike_times],
                't_end': self.t}

    def _advance_synapses(self, t0, t1, arrivals):
        """Move the synaptic state to the end of the window."""
        syn = self.syn
        if not len(syn.pre):
            return
        w = t1 - t0
        syn.A *= np.exp(-w / syn.tau_r)
        syn.B *= np.exp(-w / syn.tau_d)
        for arrival, idx in arrivals:
            rem = t1 - arrival
            syn.A[idx] += np.exp(-rem / syn.tau_r[idx])
            syn.B[idx] += np.exp(-rem / syn.tau_d[idx])

    # ---- integration -----------------------------------------------------

    def _derivatives_for(self, i, syn_current, I_ext_i):
        cell = self.neurons[i]

        def f(t, y):
            V, m, h, n = y
            self.n_evaluations += 1
            I_syn = syn_current(t, V)
            I_Na = cell.g_Na * m ** 3 * h * (V - cell.E_Na)
            I_K = cell.g_K * n ** 4 * (V - cell.E_K)
            I_L = cell.g_L * (V - cell.E_L)
            p = cell.phi
            return np.array([
                (I_ext_i - I_syn - I_Na - I_K - I_L) / cell.C_m,
                p * (cell.alpha_m(V) * (1.0 - m) - cell.beta_m(V) * m),
                p * (cell.alpha_h(V) * (1.0 - h) - cell.beta_h(V) * h),
                p * (cell.alpha_n(V) * (1.0 - n) - cell.beta_n(V) * n),
            ])
        return f

    def _V_snapshot_with(self, i, V):
        snap = self._V_frozen
        snap[i] = V
        return snap

    def simulate(self, duration, I_ext=None, dt0=0.01, record_dt=None):
        """Run the network for `duration` ms.

        Returns spike times per neuron, plus a voltage trace on a uniform
        grid if `record_dt` is given.
        """
        if self.syn is None:
            self._build()
        I_ext = (np.zeros(self.n) if I_ext is None
                 else np.broadcast_to(np.asarray(I_ext, dtype=float),
                                      (self.n,)))
        n_windows = int(np.ceil(duration / self.min_delay))
        integrators = [AdaptiveIntegrator(rtol=self.rtol, atol=self.atol,
                                          dt_max=self.dt_max)
                       for _ in range(self.n)]
        recorded_t, recorded_V = [], []

        for w in range(n_windows):
            t0 = self.t
            t1 = min(t0 + self.min_delay, duration)
            if t1 <= t0:
                break
            arrivals = self._window_arrivals(t0, t1)
            by_conn = self._arrivals_by_connection(arrivals)

            fired = []
            for i in range(self.n):
                f = self._derivatives_for(
                    i, self._synaptic_current_for(i, t0, by_conn), I_ext[i])
                resume = integrators[i].last_dt
                ts, ys = integrators[i].integrate(
                    f, self.state[i], (t0, t1),
                    dt0=(resume if resume else dt0))
                # threshold crossings, interpolated
                V = ys[:, 0]
                above = V > self.spike_threshold
                cross = np.flatnonzero(~above[:-1] & above[1:])
                for c in cross:
                    v0, v1 = V[c], V[c + 1]
                    frac = (self.spike_threshold - v0) / (v1 - v0)
                    t_spike = ts[c] + frac * (ts[c + 1] - ts[c])
                    self.spike_times[i].append(float(t_spike))
                    fired.append((i, float(t_spike)))
                self.state[i] = ys[-1]
                if record_dt is not None:
                    recorded_t.append(ts)
                    recorded_V.append(ys[:, 0])

            for integ in integrators:
                self.n_steps += integ.stats['accepted']
                self.dt_samples.extend(integ.dt_history)
                carry = integ.last_dt
                integ.reset_stats()
                integ.last_dt = carry       # survives the stats reset

            self._advance_synapses(t0, t1, arrivals)
            # queue the spikes emitted in this window
            if len(self.syn.pre):
                for i, t_spike in fired:
                    for idx in np.flatnonzero(self.syn.pre == i):
                        self._pending.append(
                            (t_spike + self.syn.delay[idx], int(idx)))
            self.t = t1

        out = {'spikes': [np.array(s) for s in self.spike_times],
               't_end': self.t}
        if record_dt is not None and recorded_t:
            grid = np.arange(0.0, self.t, record_dt)
            V = np.full((len(grid), self.n), np.nan)
            for i in range(self.n):
                ts = np.concatenate(recorded_t[i::self.n])
                vs = np.concatenate(recorded_V[i::self.n])
                order = np.argsort(ts)
                V[:, i] = np.interp(grid, ts[order], vs[order])
            out['t'], out['V'] = grid, V
        return out
