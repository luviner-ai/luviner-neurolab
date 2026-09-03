"""
Modulation of a circuit's intrinsic properties, as an alternative to
modulating its synapses.

The stomatogastric ganglion runs different rhythms out of the same thirty
neurons under different neuromodulators, and the modulators act mostly on
*ion channels* -- time constants, thresholds, intrinsic currents -- not on
synaptic weights. Every multi-task mechanism in machine learning does the
opposite: FiLM, adapters and hypernetworks all modulate weights.

The two are not equally expensive. Intrinsic parameters are order `N` per
modulator state; synaptic weights are order `N^2`. If the intrinsic route
reaches the same behaviours, the multifunctionality per parameter is
higher by a factor of order `N`, and that is a claim about what a
network's task-switching mechanism should touch.

`ModulatedCircuit` makes the intrinsic parameters per-neuron and
differentiable, so the same gradient that designs a circuit by its
synapses can design it by its channels instead. Two are exposed, and they
are the two neuromodulators actually reach:

    V_T   the spike threshold, which the Traub-Miles rate functions are
          all written relative to, so moving it shifts the whole
          spike-generating machinery of that cell together;
    g_L   the leak, which sets the membrane time constant tau = C_m/g_L.

Both are held positive-definite where they must be, and a modulator state
is just a pair of vectors that can be swapped in between runs.

References:
    Marder E (2012). "Neuromodulation of neuronal circuits: back to the
        future." Neuron 76:1-11.
    Marder E, Goaillard JM (2006). "Variability, compensation and
        homeostasis in neuron and network function."
        Nat Rev Neurosci 7:563-574.
    Perez FJ, Vidal-Gadea AG (2019) and the wider graded-transmission
        literature cited in `diffcircuit.py`.
"""

import numpy as np

from . import autodiff as ad
from .autodiff import Var
from .diffcircuit import DifferentiableCircuit


class ModulatedCircuit(DifferentiableCircuit):
    """A differentiable circuit whose intrinsic parameters are per-neuron.

    `V_T` and `g_L` become length-`n` `Var`s. Everything else -- the
    kinetics, the graded synapses, the integrator -- is inherited
    unchanged, so a circuit modulated intrinsically and the same circuit
    modulated synaptically differ in exactly one thing.
    """

    G_L_MIN, G_L_MAX = 0.02, 0.6
    V_T_MIN, V_T_MAX = -75.0, -45.0

    def __init__(self, n_neurons, **kw):
        V_T0 = kw.pop('V_T', -63.0)
        g_L0 = kw.pop('g_L', 0.1)
        # The base class assigns scalars to self.V_T and self.g_L, which
        # are read-only properties here. Give it somewhere harmless to
        # write, then install the real parameters over the top.
        object.__setattr__(self, '_base_scalars', {})
        super().__init__(n_neurons, V_T=V_T0, g_L=g_L0, **kw)
        span = self.V_T_MAX - self.V_T_MIN
        u = np.clip((float(V_T0) - self.V_T_MIN) / span, 1e-4, 1 - 1e-4)
        self.v_t_raw = Var(np.full(n_neurons, float(np.log(u / (1 - u)))))
        self.log_g_L = Var(np.full(n_neurons, float(np.log(g_L0))))

    # ---- the intrinsics, as the equations see them -----------------------

    @property
    def V_T(self):
        span = self.V_T_MAX - self.V_T_MIN
        return self.V_T_MIN + span * ad.sigmoid(self.v_t_raw)

    @V_T.setter
    def V_T(self, value):
        # Swallowed: the base class sets a scalar in __init__, before the
        # real parameter exists. Nothing else writes here.
        self._base_scalars['V_T'] = value

    @property
    def g_L(self):
        return ad.exp(self.log_g_L)

    @g_L.setter
    def g_L(self, value):
        self._base_scalars['g_L'] = value

    # ---- the two parameter sets, kept separable --------------------------

    @property
    def intrinsic_parameters(self):
        """Order N per modulator state."""
        return [self.v_t_raw, self.log_g_L]

    @property
    def synaptic_parameters(self):
        """Order N^2 per modulator state.

        The conductance is held as `log_g`, which is what the optimizer
        moves and what keeps it positive without a clamp.
        """
        return [s.log_g for _, _, s in self.connections]

    @property
    def parameters(self):
        return self.intrinsic_parameters + self.synaptic_parameters

    def clamp_intrinsic(self):
        """Keep the leak in the range a cell can actually have.

        The threshold needs no clamp -- its coordinate is bounded by
        construction. The leak does: without it the optimizer walks it to
        zero, which is not a modulator state but a broken cell, since the
        membrane time constant diverges and the forward pass stops
        meaning anything.
        """
        self.log_g_L.value = np.clip(self.log_g_L.value,
                                     np.log(self.G_L_MIN),
                                     np.log(self.G_L_MAX))
        return self

    def clamp_synaptic(self, lo=1e-4, hi=0.4):
        for _, _, s in self.connections:
            s.log_g.value = np.clip(s.log_g.value, np.log(lo), np.log(hi))
        return self

    # ---- modulator states ------------------------------------------------

    def intrinsic_state(self):
        return (self.v_t_raw.value.copy(), self.log_g_L.value.copy())

    def load_intrinsic(self, state):
        self.v_t_raw.value = np.array(state[0], dtype=float)
        self.log_g_L.value = np.array(state[1], dtype=float)
        return self

    def readable_intrinsics(self):
        """The two vectors in the units a physiologist would quote."""
        span = self.V_T_MAX - self.V_T_MIN
        u = 1.0 / (1.0 + np.exp(-self.v_t_raw.value))
        return self.V_T_MIN + span * u, np.exp(self.log_g_L.value)

    def synaptic_state(self):
        return [np.array(s.log_g.value, dtype=float, copy=True)
                for _, _, s in self.connections]

    def load_synaptic(self, state):
        for (_, _, s), v in zip(self.connections, state):
            s.log_g.value = np.array(v, dtype=float)
        return self


def modulated_fully_connected(n_neurons, g_init=0.05,
                              excitatory_fraction=0.5, V_half=-50.0, k=4.0,
                              tau=5.0, seed=0, types=None, **kw):
    """`fully_connected`, but on a `ModulatedCircuit`.

    Duplicated rather than parameterised because the original belongs to
    another line of work and this one has to stay swappable against it.
    """
    from .diffcircuit import GradedSynapse
    rng = np.random.RandomState(seed)
    if types is None:
        types = [rng.rand() < excitatory_fraction for _ in range(n_neurons)]
    c = ModulatedCircuit(n_neurons, **kw)
    for pre in range(n_neurons):
        for post in range(n_neurons):
            if pre == post:
                continue
            E_rev = 0.0 if types[pre] else -80.0
            c.connect(pre, post,
                      GradedSynapse(g_max=g_init * (1.0 + 0.1 * rng.randn()),
                                    E_rev=E_rev, V_half=V_half, k=k, tau=tau))
    c.types = types
    return c
