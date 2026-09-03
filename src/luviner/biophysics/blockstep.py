"""
Block timesteps: adaptive stepping that vectorizes.

Per-neuron adaptive stepping wins about 10x in work (adaptive_network.py),
and the win holds as the network grows. But it cannot be vectorized: each
cell is on its own step, so each needs its own integrator, and the Python
loop over cells costs far more than the arithmetic it saves.

The trick, borrowed from N-body astrophysics where the same problem
appears with orbits of wildly different periods, is to stop giving every
body its own arbitrary step and instead quantize steps to powers of two:

    dt_i = dt_base * 2^k_i

A window of 2^K base steps is then divided so that at sub-step number i,
exactly those neurons whose level k satisfies (i mod 2^k == 0) are due for
an update. Everything at the same level shares a step size, so a level is
one vectorized NumPy operation over an array of neurons rather than a loop.

The synchronization is free rather than enforced: because the steps are
powers of two, cells at every level land on common boundaries
automatically, and a cell at level 3 is always up to date whenever a cell
at level 0 needs to read it.

The saving is the same as before -- quiescent cells sit at high levels and
are touched rarely, spiking cells drop to level 0 -- but now the work at
each level is vectorized. In a sparse cortical network most cells are at
the top level on any given sub-step, so most sub-steps update only a
handful of neurons.

Level selection is driven by the rate of change of the membrane potential,
which is what actually limits the step in a Hodgkin-Huxley cell: during an
upstroke dV/dt reaches ~300 mV/ms and the step must be small; at rest it
is essentially zero. A cell may drop several levels at once when it starts
to move, but rises one level at a time, which prevents oscillation between
levels.

Reference for the scheme:
    Aarseth SJ (2003). Gravitational N-Body Simulations. Cambridge.
    (block time steps, ch. 2)
"""

import numpy as np

from .hodgkin_huxley import HodgkinHuxleyNeuron as _H


class BlockTimestepNetwork:
    """Vectorized network with power-of-two adaptive steps per neuron.

    Parameters
    ----------
    n_neurons : number of cells.
    kinetics : 'squid' or 'cortical', matching NeuronPopulation.
    dt_base : the finest step, in ms. Cells resolving a spike run here.
    n_levels : how many doublings above dt_base are allowed. The coarsest
        step is dt_base * 2^(n_levels-1), capped by min_delay.
    dv_tol : the target voltage change per step, in mV. This is what
        selects a cell's level; smaller means finer steps everywhere.
    min_delay : synchronization window, ms. Must not exceed the smallest
        axonal delay in the network.
    """

    def __init__(self, n_neurons, kinetics='cortical', dt_base=0.01,
                 n_levels=7, dv_tol=0.25, min_delay=1.0,
                 spike_threshold=0.0, **params):
        from .population import NeuronPopulation
        self.pop = (NeuronPopulation.cortical(n_neurons, **params)
                    if kinetics == 'cortical'
                    else NeuronPopulation(n_neurons, **params))
        self.n = n_neurons
        self.dt_base = float(dt_base)
        self.min_delay = float(min_delay)
        self.dv_tol = float(dv_tol)
        self.spike_threshold = spike_threshold

        # The coarsest level may not overshoot the synchronization window.
        max_by_delay = int(np.floor(np.log2(max(min_delay / dt_base, 1.0))))
        self.n_levels = int(max(1, min(n_levels, max_by_delay + 1)))
        self.steps_per_window = int(round(min_delay / dt_base))
        self.reset()

    # ---- construction ----------------------------------------------------

    def connect(self, pre, post, synapse, weights=None):
        self.pop.connect(pre, post, synapse, weights=weights)
        return self

    def reset(self):
        self.pop.reset()
        self.level = np.full(self.n, self.n_levels - 1, dtype=int)
        self.t = 0.0
        self.spike_times = [[] for _ in range(self.n)]
        self.n_updates = 0          # neuron-updates actually performed
        self.n_substeps = 0
        self._level_census = np.zeros(self.n_levels, dtype=np.int64)
        return self

    @property
    def stats(self):
        """Work done, and where the neurons were sitting.

        `updates` counts neuron-updates: the quantity a fixed-step scheme
        would pay n_neurons * n_substeps for.
        """
        fixed = self.n * self.n_substeps
        census = self._level_census.astype(float)
        return {'updates': int(self.n_updates), 'substeps': self.n_substeps,
                'fixed_equivalent': int(fixed),
                'reduction': (fixed / self.n_updates) if self.n_updates else 0.0,
                'level_occupancy': (census / census.sum()).tolist()
                if census.sum() else [],
                'mean_level': float((census * np.arange(self.n_levels)).sum()
                                    / census.sum()) if census.sum() else 0.0}

    # ---- level assignment ------------------------------------------------

    def _choose_levels(self, dVdt):
        """Pick each neuron's level from how fast its voltage is moving.

        A cell may fall any distance at once -- a spike starts abruptly
        and must be caught immediately -- but climbs only one level per
        window, so a cell that has just finished spiking does not jump
        straight back to the coarsest step and miss its own afterpotential.
        """
        # With exponential Euler the gates are exact, so the step is no
        # longer capped by their time constants -- only by how fast the
        # coefficients (i.e. the voltage) are moving.
        speed = np.abs(dVdt) + 1e-12
        ideal_dt = self.dv_tol / speed
        want = np.floor(np.log2(np.maximum(ideal_dt / self.dt_base, 1.0)))
        want = np.clip(want, 0, self.n_levels - 1).astype(int)
        rising = want > self.level
        self.level = np.where(rising, self.level + 1, want)
        return self.level

    # ---- integration -----------------------------------------------------

    def _step_subset(self, idx, dt):
        """Advance the neurons in `idx` by dt with exponential Euler.

        Hodgkin-Huxley is stiff: the sodium activation gate has tau ~ 0.05
        ms even when the voltage is not moving at all, so an explicit
        method is limited by that gate at every step regardless of what V
        is doing. That, not the spike itself, is what pins conventional
        simulations to a small step.

        But every equation here has the same form, dx/dt = (x_inf - x)/tau,
        whose solution over a step of constant coefficients is exact:

            x(t + dt) = x_inf + (x - x_inf) * exp(-dt / tau)

        For the gates x_inf and tau come straight from the rate functions.
        For the voltage the same rearrangement applies, with the total
        conductance setting tau_V and the conductance-weighted mean of the
        reversal potentials setting V_inf.

        This is unconditionally stable and exact at the frozen-coefficient
        limit, so the step is limited by how fast the coefficients change
        rather than by the stiffest time constant. It is the method NEURON
        calls `cnexp`.

        Vectorized over the subset: a handful of NumPy calls regardless of
        how many neurons are in it.
        """
        pop = self.pop
        V, m, h, n = pop.V[idx], pop.m[idx], pop.h[idx], pop.nn[idx]
        I = self._I_total[idx]
        p = pop.phi[idx]

        am, bm, ah, bh, an, bn = self._rates_subset(V, idx)
        # gates: exact update at frozen V
        for rate_a, rate_b, arr in ((am, bm, 'm'), (ah, bh, 'h'), (an, bn, 'n')):
            tot = (rate_a + rate_b) * p
            inf = rate_a / (rate_a + rate_b)
            decay = np.exp(-dt * tot)
            cur = {'m': m, 'h': h, 'n': n}[arr]
            new_val = inf + (cur - inf) * decay
            if arr == 'm':
                m_new = new_val
            elif arr == 'h':
                h_new = new_val
            else:
                n_new = new_val

        # voltage: same form, with conductances at the midpoint of the gate
        # update so the step stays second-order in the slow variable
        mm = 0.5 * (m + m_new)
        hh = 0.5 * (h + h_new)
        nn_ = 0.5 * (n + n_new)
        g_Na = pop.g_Na[idx] * mm ** 3 * hh
        g_K = pop.g_K[idx] * nn_ ** 4
        g_L = pop.g_L[idx]
        g_tot = g_Na + g_K + g_L
        V_inf = ((g_Na * pop.E_Na[idx] + g_K * pop.E_K[idx]
                  + g_L * pop.E_L[idx] + I) / g_tot)
        tau_V = pop.C_m[idx] / g_tot
        V_new = V_inf + (V - V_inf) * np.exp(-dt / tau_V)

        pop.V[idx] = V_new
        pop.m[idx] = np.clip(m_new, 0.0, 1.0)
        pop.h[idx] = np.clip(h_new, 0.0, 1.0)
        pop.nn[idx] = np.clip(n_new, 0.0, 1.0)
        self.n_updates += len(idx)
        return (V_inf - V) / tau_V      # dV/dt at the start of the step

    def _rates_subset(self, V, idx):
        pop = self.pop
        if pop.kinetics == 'squid':
            return (_H.alpha_m(V), _H.beta_m(V), _H.alpha_h(V),
                    _H.beta_h(V), _H.alpha_n(V), _H.beta_n(V))
        from .population import _xo
        x = V - pop.V_T[idx]
        return (1.28 * _xo((x - 13.0) / 4.0),
                1.4 * _xo(-(x - 40.0) / 5.0),
                0.128 * np.exp(-(x - 17.0) / 18.0),
                4.0 / (1.0 + np.exp(-(x - 40.0) / 5.0)),
                0.16 * _xo((x - 15.0) / 5.0),
                0.5 * np.exp(-(x - 10.0) / 40.0))

    def simulate(self, duration, I_ext=None, record_spikes=True):
        """Run the network. Synapses advance on the base step."""
        pop = self.pop
        if not pop._built or pop._dt_built != self.dt_base:
            pop._build(self.dt_base)
        I_ext = (np.zeros(self.n) if I_ext is None
                 else np.broadcast_to(np.asarray(I_ext, float), (self.n,)))

        n_windows = int(np.ceil(duration / self.min_delay))
        dVdt = np.zeros(self.n)

        for _ in range(n_windows):
            self._choose_levels(dVdt)
            np.add.at(self._level_census, self.level, 1)

            # Group neurons by level once per window rather than at every
            # sub-step. The membership cannot change inside a window, so
            # recomputing it 2^K times was pure overhead -- and with the
            # arithmetic now reduced several-fold, that overhead was what
            # the wall clock was actually spending its time on.
            by_level = [np.flatnonzero(self.level == k)
                        for k in range(self.n_levels)]
            has_syn = bool(pop.n_syn)

            for i in range(self.steps_per_window):
                self._I_total = (I_ext - pop.synaptic_currents() if has_syn
                                 else I_ext)
                due_levels = [k for k in range(self.n_levels)
                              if len(by_level[k]) and (i & ((1 << k) - 1)) == 0]
                if due_levels:
                    fired_idx = []
                    for k in due_levels:
                        sel = by_level[k]
                        V_before = pop.V[sel]
                        d = self._step_subset(sel, self.dt_base * (1 << k))
                        dVdt[sel] = d
                        crossed = ((V_before <= self.spike_threshold)
                                   & (pop.V[sel] > self.spike_threshold))
                        if crossed.any():
                            fired_idx.append(sel[crossed])
                    if fired_idx:
                        fired_all = np.concatenate(fired_idx)
                        if record_spikes:
                            t_now = self.t + i * self.dt_base
                            for j in fired_all:
                                self.spike_times[j].append(t_now)
                        self._deliver(fired_all)
                if has_syn:
                    self._advance_synapses()
                self.n_substeps += 1
            self.t += self.min_delay

        return {'spikes': [np.array(s) for s in self.spike_times],
                't_end': self.t}

    def _deliver(self, fired_idx):
        pop = self.pop
        if not pop.n_syn or not len(fired_idx):
            return
        mark = np.zeros(self.n, dtype=bool)
        mark[fired_idx] = True
        pop._ptr = (pop._ptr + 1) % pop._buf_len
        pop._spike_buf[pop._ptr] = mark
        self._delivered = True

    def _advance_synapses(self):
        pop = self.pop
        if not pop.n_syn:
            return
        if not getattr(self, '_delivered', False):
            pop._ptr = (pop._ptr + 1) % pop._buf_len
            pop._spike_buf[pop._ptr] = False
        self._delivered = False
        arrived = pop._spike_buf[
            (pop._ptr - pop.delay_steps) % pop._buf_len, pop.pre]
        pop.A = (pop.A + arrived) * pop.decay_r
        pop.B = (pop.B + arrived) * pop.decay_d
