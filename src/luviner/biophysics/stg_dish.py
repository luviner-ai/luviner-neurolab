"""STG bursting cells wired with the dish's plastic synapses.

Three mechanisms failed to give this project culture-length bursts:
short-term depression (`## Run: N8b`), a voltage-gated slow K+ current
and a calcium-gated one (`## Run: dish-slow`, `## Run: dish-slow2a`).
The wall stands even in an isolated cortical cell. The pyloric lane's
AB/PD cell clears every criterion by construction -- 0.687 Hz, 494 ms
self-terminating bursts, 30 spikes each -- because it is *bistable with a
slow variable*, which is what bursting actually requires and what a
monostable cortical cell cannot be talked into.

This module puts those cells in the dish's synaptic setting: spike-
mediated excitatory synapses with the dish's own AMPA/NMDA kinetics,
short-term depression, STDP under a homeostatic budget, and structured
inhibition. The question it exists to answer is whether the long bursts
**survive the network without latching**.

**These are not cortical cells.** This preparation tests the mechanism
question -- do long structured events read the weight order -- and drops
the culture-fidelity claim explicitly. That is a deliberate trade, not an
oversight.

How it couples to the STG lane
------------------------------
The STG network integrates in *conductance* form (exponential Euler on
`g_tot` and `gE_tot`) and exposes a per-step API with an `I_ext` hook.
Synaptic current is therefore injected as a current computed from the
voltage at the **start** of the step, rather than folded into `g_tot`
where a conductance-based synapse properly belongs. That is a
first-order approximation and `tests/test_stg_dish.py` measures its
error against a halved timestep rather than assuming it is small.

The driver contract, so this module can be built and tested before the
other lane's commit lands:

    driver.n                  number of cells
    driver.V                  (n,) mV, current voltages
    driver.step(dt, I_ext)    advance one step, return (n,) bool of spikes

Anything satisfying that works, including the fake in the tests.
`STGDriver` below adapts the real `stg.STGNetwork` to it.

**Their step returns events, not a boolean, and the events are one step
late**: a peak at sample `i` is only confirmed while computing `i + 1`,
which is the price of causal peak detection. The lateness is uniform
across cells, so pair *timing differences* -- which is all STDP reads --
are unaffected; only the absolute alignment to the stimulus shifts by
`dt`. Their peak times are used for the recorded spike train, which is
more accurate than counting steps here.
"""

import numpy as np

from .plasticity import PlasticConnections, WeightNormalization
from .shortterm import ShortTermDepression
from .synapse import Synapse


class SpikeSynapses:
    """Spike-driven two-state conductances, the dish's kinetics exactly.

    The same `(A, B)` pair, `norm`, circular delay buffer and short-term
    depression hook as `NeuronPopulation`, lifted out so they can be
    driven by any spike source. Kept deliberately identical: if this
    diverged from the dish's synapse the comparison between preparations
    would be between two synapse models rather than two cell models.
    """

    def __init__(self, pre, post, kinds, weights, dt, stp=None):
        self.pre = np.asarray(pre, dtype=int)
        self.post = np.asarray(post, dtype=int)
        self.n_syn = self.pre.size
        self.weights = np.asarray(weights, dtype=float).copy()
        g = [k.g_max for k in kinds]
        self.g_max = np.asarray(g, dtype=float)
        self.E_rev = np.array([k.E_rev for k in kinds], dtype=float)
        self.tau_r = np.array([k.tau_rise for k in kinds], dtype=float)
        self.tau_d = np.array([k.tau_decay for k in kinds], dtype=float)
        self.mg = np.array([k.mg_conc if k.mg_block else 0.0 for k in kinds],
                           dtype=float)
        self.has_mg = self.mg > 0.0
        delay = np.array([k.delay for k in kinds], dtype=float)
        t_peak = ((self.tau_r * self.tau_d) / (self.tau_d - self.tau_r)
                  * np.log(self.tau_d / self.tau_r))
        self.norm = 1.0 / (np.exp(-t_peak / self.tau_d)
                           - np.exp(-t_peak / self.tau_r))
        self.stp = stp
        self.reset(dt)
        self.delay_steps = np.maximum(
            np.round(delay / dt).astype(int) - 1, 0)
        self._buf_len = int(self.delay_steps.max()) + 2 if self.n_syn else 2
        self._buf = np.zeros((self._buf_len, int(self.post.max()) + 1
                              if self.n_syn else 1), dtype=bool)
        self._ptr = 0

    def reset(self, dt):
        self.dt = float(dt)
        self.decay_r = np.exp(-dt / self.tau_r)
        self.decay_d = np.exp(-dt / self.tau_d)
        self.A = np.zeros(self.n_syn)
        self.B = np.zeros(self.n_syn)
        if self.stp is not None:
            self.stp.reset(self.n_syn, dt)
        return self

    def currents(self, V, n_cells):
        """Synaptic current into each cell, uA/cm^2, sink convention.

        Positive means current leaving the cell, matching
        `NeuronPopulation.synaptic_currents`. The caller negates it to
        feed an `I_ext` hook.
        """
        I = np.zeros(n_cells)
        if not self.n_syn:
            return I
        g = self.weights * self.g_max * self.norm * (self.B - self.A)
        V_post = V[self.post]
        if self.has_mg.any():
            unblock = 1.0 / (1.0 + np.exp(-0.062 * V_post) * (self.mg / 3.57))
            g = np.where(self.has_mg, g * unblock, g)
        np.add.at(I, self.post, g * (V_post - self.E_rev))
        return I

    def advance(self, fired):
        """Decay the conductances and deliver the spikes arriving now."""
        if not self.n_syn:
            return
        self._ptr = (self._ptr + 1) % self._buf_len
        self._buf[self._ptr] = fired
        rows = (self._ptr - self.delay_steps) % self._buf_len
        # The BOOLEAN array, not a float cast: `ShortTermDepression.release`
        # masks with `arrived & self.mask`, and `population.py` passes the
        # raw booleans. Casting here diverged from the dish's synapse in
        # exactly the path the kinetics-equality test did not cover,
        # because that test ran with `stp=None`.
        arrived = self._buf[rows, self.pre]
        kick = arrived if self.stp is None else self.stp.release(arrived,
                                                                 self.dt)
        self.A = (self.A + kick) * self.decay_r
        self.B = (self.B + kick) * self.decay_d


class STGDriver:
    """`stg.STGNetwork` behind the driver contract.

    Adapts their signature to this module's, rather than the other way
    round: their file, their signature. `reset_state` must have been
    called -- their docstring names driving with `step()` as exactly this
    use case.
    """

    def __init__(self, net, V0=-65.0, rng=None, v0_jitter=0.0):
        self.net = net
        self.n = int(net.n)
        net.reset_state(V0=V0, rng=rng, v0_jitter=v0_jitter)
        self.events = []

    @property
    def V(self):
        return self.net.V

    def step(self, dt, I_ext=0.0):
        self.events = self.net.step(dt, I_ext)
        fired = np.zeros(self.n, dtype=bool)
        for cell, _t in self.events:
            fired[cell] = True
        return fired


def burst_events(spike_times, window_ms, gap_ms=50.0, min_spikes=3):
    """Bursts from one cell's spike train, with the saturation guard.

    The guard is re-specified by purpose rather than by a duty threshold:
    **at least three discrete events in the window, and no single event
    spanning more than 20% of it.** Latching -- one event covering the
    whole window -- fails it; a clean 1.46 s limit cycle at 34% duty
    passes. The earlier `duty < 0.2` form rejected physiological AB/PD
    bursting, whose published window is 0.3-0.4.

    A burst is at least `min_spikes` spikes: two spikes once is not a
    bursting regime, which is the error `dish-slow2a` caught in itself.

    Returns (starts, durations, sizes, ok).
    """
    sp = np.asarray(spike_times, dtype=float)
    if sp.size == 0:
        return np.zeros(0), np.zeros(0), np.zeros(0, dtype=int), False
    brk = np.flatnonzero(np.diff(sp) >= gap_ms)
    i0 = np.concatenate([[0], brk + 1])
    i1 = np.concatenate([brk, [sp.size - 1]])
    starts, durs = sp[i0], sp[i1] - sp[i0]
    sizes = i1 - i0 + 1
    real = sizes >= min_spikes
    ok = bool(real.sum() >= 3
              and (durs[real].max() <= 0.20 * window_ms if real.any()
                   else False))
    return starts[real], durs[real], sizes[real], ok


class STGDish:
    """Excitatory STG cells, dish synapses, structured inhibition.

    Only the E->E synapses are plastic, as in `SiliconDish`. E->I and
    I->E are fixed conductances, and the inhibition is global by default
    -- the arrangement `## Run: dish-2` pinned.

    **The synaptic conductance scale is NOT yet measured for these
    cells and the defaults here are the dish's, not the STG lane's.**
    `w_rec = 0.015` and the rule's `w_max = 0.06` are right for cortical
    cells with `g_L = 0.1`; STG cells carry maximal conductances two to
    three orders larger. Two library defaults have already been wrong for
    this project by 5x (`I_M`) and 20-60x (`I_KCa`), both caught only by
    measuring. A conductance sweep on one coupled pair is the required
    step before the pilot, and `budget` must satisfy
    `budget <= w_max * min_fan_in` or homeostasis cannot reach it -- the
    clip wins silently and every weight sits at the rail.
    """

    def __init__(self, driver_factory, n_exc=10, n_inh=3, connectivity=0.15,
                 w_rec=0.015, g_ei=0.3, g_ie=3.0, nmda_ratio=2.0,
                 budget=None, tau_norm=200.0, rule=None, dt=0.05,
                 depress=(0.20, 400.0), seed=0):
        self.n_exc, self.n_inh = int(n_exc), int(n_inh)
        self.n = self.n_exc + self.n_inh
        self.dt = float(dt)
        self.seed = int(seed)
        self.driver = driver_factory(self.n)
        if self.driver.n != self.n:
            raise ValueError(f"driver has {self.driver.n} cells, expected "
                             f"{self.n}")

        rng = np.random.RandomState(seed)
        mask = rng.rand(self.n_exc, self.n_exc) < connectivity
        np.fill_diagonal(mask, False)
        pre_e, post_e = (a.copy() for a in np.nonzero(mask))
        w0 = np.clip(w_rec * (1.0 + 0.2 * rng.randn(pre_e.size)), 0.0, None)
        self.pre_e, self.post_e, self.w0 = pre_e, post_e, w0

        exc, inh = np.arange(self.n_exc), self.n_exc + np.arange(self.n_inh)
        pre_ei, post_ei = (x.ravel() for x in
                           np.meshgrid(exc, inh, indexing='ij'))
        pre_ie, post_ie = (x.ravel() for x in
                           np.meshgrid(inh, exc, indexing='ij'))

        pre = np.concatenate([pre_e, pre_e, pre_ei, pre_ie])
        post = np.concatenate([post_e, post_e, post_ei, post_ie])
        kinds = ([Synapse.ampa(g_max=1.0, delay=1.0)] * pre_e.size
                 + [Synapse.nmda(g_max=1.0, delay=1.0)] * pre_e.size
                 + [Synapse.ampa(g_max=1.0, delay=1.0)] * pre_ei.size
                 + [Synapse.gaba_a(g_max=1.0, delay=1.0, E_rev=-75.0)]
                 * pre_ie.size)
        weights = np.concatenate([w0, w0 * nmda_ratio,
                                  np.full(pre_ei.size, g_ei),
                                  np.full(pre_ie.size, g_ie)])
        self._slice_e = slice(0, pre_e.size)
        self._slice_nmda = slice(pre_e.size, 2 * pre_e.size)

        stp = None
        if depress is not None:
            m = np.zeros(pre.size, dtype=bool)
            m[self._slice_e] = True
            m[self._slice_nmda] = True
            stp = ShortTermDepression(U=depress[0], tau_rec=depress[1],
                                      mask=m)
        self.syn = SpikeSynapses(pre, post, kinds, weights, self.dt, stp=stp)

        from .dish import rule_with_integral
        homeo = (WeightNormalization(budget=budget, tau=tau_norm,
                                     w_max=(rule or rule_with_integral(
                                         -0.134)).w_max)
                 if budget is not None else None)
        self.conn = PlasticConnections(
            pre_e, post_e, w0.copy(),
            rule=rule if rule is not None else rule_with_integral(-0.134),
            n_neurons=self.n, homeostasis=homeo)

        # Cells that cannot reach the budget at ANY weight -- `fan_in *
        # w_max < budget` -- are a wiring fact, not a dynamical one, and
        # `## Run: phase-link` established the precedent: exclude them by
        # a criterion computable from the wiring alone, before any data
        # is looked at, rather than let the clip act silently.
        #
        # The first version RAISED here, which was too strong: at 48 cells
        # and connectivity 0.15 one wiring in ten draws a cell with
        # fan-in 1, and refusing the whole wiring over 2 cells of 48
        # discards a usable experiment. It is recorded instead.
        self.capped = np.zeros(self.n_exc, dtype=bool)
        if budget is not None:
            k_in = np.zeros(self.n_exc)
            np.add.at(k_in, post_e, 1.0)
            has_input = k_in > 0
            self.capped = has_input & (self.conn.rule.w_max * k_in < budget)
            if self.capped.all():
                raise ValueError(
                    f"budget {budget} is above w_max * fan-in for EVERY "
                    f"cell (max reachable "
                    f"{self.conn.rule.w_max * k_in.max():.4g}); homeostasis "
                    f"would clip at the rail everywhere. Lower the budget, "
                    f"raise w_max, or connect more densely.")

        self.t = 0.0
        self.spike_times = [[] for _ in range(self.n)]

    # ---- running ---------------------------------------------------------

    def _push_weights(self):
        w = self.conn.weights
        self.syn.weights[self._slice_e] = w
        self.syn.weights[self._slice_nmda] = w * 2.0

    def run(self, duration, learn=True, extra=None):
        """Advance the dish.

        `extra` adds an external current, either a fixed per-cell array
        or a callable of the current time returning one. The callable
        form is what closed-loop stimulation needs: it is invoked before
        each step, so it may read `spike_times` up to the previous step
        and decide what to inject now. It must not write to the dish.
        """
        n_steps = int(round(duration / self.dt))
        live = callable(extra)
        for _ in range(n_steps):
            I_syn = self.syn.currents(self.driver.V, self.n)
            I_ext = -I_syn                      # the hook is depolarizing
            if live:
                I_ext = I_ext + extra(self.t)
            elif extra is not None:
                I_ext = I_ext + extra
            fired = np.asarray(self.driver.step(self.dt, I_ext), dtype=bool)
            self.t += self.dt
            events = getattr(self.driver, 'events', None)
            if events:                      # their peak times, not step counts
                for cell, t_pk in events:
                    self.spike_times[cell].append(float(t_pk))
            else:
                for i in np.flatnonzero(fired):
                    self.spike_times[i].append(self.t)
            self.syn.advance(fired)
            self.conn.step(fired, dt=self.dt, learn=learn)
            if learn or self.conn.homeostasis is not None:
                self._push_weights()
        return self

    # ---- readout ---------------------------------------------------------

    @property
    def weights(self):
        return self.conn.weights

    def spikes(self, i, t_lo=0.0, t_hi=None):
        s = np.asarray(self.spike_times[i], dtype=float)
        t_hi = self.t if t_hi is None else t_hi
        return s[(s > t_lo) & (s <= t_hi)]

    def bursts(self, i, t_lo=0.0, t_hi=None, **kw):
        t_hi = self.t if t_hi is None else t_hi
        return burst_events(self.spikes(i, t_lo, t_hi), t_hi - t_lo, **kw)
