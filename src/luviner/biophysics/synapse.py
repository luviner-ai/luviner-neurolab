"""
Conductance-based chemical synapses.

A synapse is not a weight. In a real neuron a presynaptic spike opens
ligand-gated ion channels on the postsynaptic membrane, and the current
that flows depends on how far the postsynaptic voltage sits from that
channel's reversal potential:

    I_syn = g_syn(t) * B(V_post) * (V_post - E_rev)

This has consequences an artificial "weight" does not reproduce:

  - The same synapse delivers less current as the postsynaptic cell
    depolarizes toward E_rev, and reverses sign beyond it.
  - Inhibition can be *shunting*: GABA_A sits near rest, so it injects
    almost no current but raises membrane conductance, dividing the
    effect of every other input. A negative weight cannot do this.
  - NMDA is gated by voltage as well as glutamate (Mg2+ block), making
    it a coincidence detector rather than a linear summator.

Conductance kinetics use the two-state dual-exponential form:

    g(t) = g_max * norm * (B - A),  dA/dt = -A/tau_rise,
                                    dB/dt = -B/tau_decay

with A, B each incremented by 1 on spike arrival. Both equations are
linear, so they are integrated exactly with an exponential step rather
than numerically -- unconditionally stable and cheaper than RK4.

Kinetic parameters are experimental values from voltage-clamp recordings
at cortical and hippocampal synapses.

References:
    Jahr CE, Stevens CF (1990). "Voltage dependence of NMDA-activated
        macroscopic conductances predicted by single-channel kinetics."
        J Neurosci 10(9):3178-3182.   [the Mg2+ block equation]
    Destexhe A, Mainen ZF, Sejnowski TJ (1994). "An efficient method for
        computing synaptic conductances." Neural Comput 6:14-18.
    Roth A, van Rossum MCW (2009). "Modeling synapses." In: Computational
        Modeling Methods for Neuroscientists, MIT Press.
"""

import numpy as np


class Synapse:
    """A conductance-based chemical synapse with axonal delay.

    Parameters
    ----------
    g_max : peak conductance, mS/cm^2
    E_rev : reversal potential, mV. Sets whether the synapse is
            excitatory (E_rev far above rest) or inhibitory (at/below).
    tau_rise, tau_decay : ms
    delay : axonal conduction + transmitter release latency, ms
    mg_block : apply the Jahr-Stevens Mg2+ block (NMDA only)
    mg_conc : extracellular [Mg2+], mM
    """

    def __init__(self, g_max, E_rev, tau_rise, tau_decay,
                 delay=1.0, mg_block=False, mg_conc=1.0, label=None):
        if tau_rise >= tau_decay:
            raise ValueError("tau_rise must be smaller than tau_decay")
        self.g_max = g_max
        self.E_rev = E_rev
        self.tau_rise = tau_rise
        self.tau_decay = tau_decay
        self.delay = delay
        self.mg_block = mg_block
        self.mg_conc = mg_conc
        self.label = label or 'synapse'

        # Normalize so that a single spike produces a peak conductance of
        # exactly g_max, independent of the rise/decay pair chosen.
        t_peak = ((tau_rise * tau_decay) / (tau_decay - tau_rise)
                  * np.log(tau_decay / tau_rise))
        self.t_peak = t_peak
        self._norm = 1.0 / (np.exp(-t_peak / tau_decay)
                            - np.exp(-t_peak / tau_rise))

        self.reset()

    # ---- presets ---------------------------------------------------------

    @classmethod
    def ampa(cls, g_max=0.05, delay=1.0):
        """Fast excitatory glutamatergic. Carries the bulk of cortical
        excitation; decay ~2 ms means it barely summates at low rates."""
        return cls(g_max, E_rev=0.0, tau_rise=0.4, tau_decay=2.0,
                   delay=delay, label='AMPA')

    @classmethod
    def nmda(cls, g_max=0.05, delay=1.0, mg_conc=1.0):
        """Slow excitatory glutamatergic with voltage-dependent Mg2+
        block. Needs glutamate AND a depolarized postsynaptic cell, so it
        fires only on coincidence -- the substrate of Hebbian plasticity."""
        return cls(g_max, E_rev=0.0, tau_rise=2.0, tau_decay=100.0,
                   delay=delay, mg_block=True, mg_conc=mg_conc, label='NMDA')

    @classmethod
    def gaba_a(cls, g_max=0.05, delay=1.0, E_rev=-70.0):
        """Fast inhibitory, chloride-mediated. E_rev sits near rest, so
        its main effect is shunting rather than hyperpolarizing."""
        return cls(g_max, E_rev=E_rev, tau_rise=0.5, tau_decay=6.0,
                   delay=delay, label='GABA_A')

    @classmethod
    def gaba_b(cls, g_max=0.02, delay=1.0):
        """Slow inhibitory, K+-mediated, strongly hyperpolarizing."""
        return cls(g_max, E_rev=-95.0, tau_rise=25.0, tau_decay=200.0,
                   delay=delay, label='GABA_B')

    # ---- state -----------------------------------------------------------

    def reset(self):
        self.A = 0.0          # fast (rise) state
        self.B = 0.0          # slow (decay) state
        self._pending = []    # spike arrival times still in the axon
        self._t = 0.0
        return self

    # ---- dynamics --------------------------------------------------------

    def receive_spike(self, t_now):
        """Register a presynaptic spike; it arrives after `delay`."""
        self._pending.append(t_now + self.delay)

    def mg_unblock(self, V_post):
        """Fraction of NMDA channels not blocked by Mg2+ (Jahr-Stevens).

        ~6% at rest, ~78% at 0 mV: the channel is essentially shut unless
        the postsynaptic cell is already depolarized.
        """
        if not self.mg_block:
            return 1.0
        return 1.0 / (1.0 + np.exp(-0.062 * V_post) * (self.mg_conc / 3.57))

    def conductance(self, V_post=None):
        """Current synaptic conductance in mS/cm^2 (Mg block included)."""
        g = self.g_max * self._norm * (self.B - self.A)
        if self.mg_block and V_post is not None:
            g = g * self.mg_unblock(V_post)
        return g

    def current(self, V_post):
        """Synaptic current in uA/cm^2. Positive = outward = hyperpolarizing,
        matching the sign convention of the ionic currents in the neuron."""
        return self.conductance(V_post) * (V_post - self.E_rev)

    def step(self, dt):
        """Advance the conductance by dt. Exact exponential integration."""
        self._t += dt
        # Deliver any spikes whose axonal delay has elapsed.
        #
        # The comparison carries half a step of tolerance. Both _t and the
        # arrival time are sums of dt in floating point, so an exact "<="
        # lands on either side of the boundary depending on accumulated
        # rounding, and the same synapse then delivers one spike after
        # `delay` and the next after `delay + dt`. Rounding to the nearest
        # step removes that jitter.
        if self._pending:
            cutoff = self._t + 0.5 * dt
            still = []
            for t_arrive in self._pending:
                if t_arrive <= cutoff:
                    self.A += 1.0
                    self.B += 1.0
                else:
                    still.append(t_arrive)
            self._pending = still
        self.A *= np.exp(-dt / self.tau_rise)
        self.B *= np.exp(-dt / self.tau_decay)
