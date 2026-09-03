"""
Additional voltage- and calcium-gated channels.

The Hodgkin-Huxley neuron carries three currents: fast Na+, delayed
rectifier K+, and leak. That is enough for an action potential and it is
what the squid giant axon has, but it is not what a cortical neuron has.
Its consequences are visible in the base model: firing starts abruptly at
~55 Hz (Class II), and the rate never drops during a sustained stimulus.
Real cortical cells start firing from arbitrarily low rates and slow down
as they fire.

The channels here are what account for that difference. Each one is an
independent gating system that contributes a current and carries its own
state variables; the neuron integrates them alongside V, m, h and n in
the same RK4 step.

    PotassiumA    I_A     transient K+, inactivating. Delays the first
                          spike after hyperpolarization, and converts
                          Class II to Class I excitability.
    PotassiumM    I_M     muscarinic K+, slow, non-inactivating. Builds
                          up over hundreds of ms -> spike frequency
                          adaptation.
    CalciumT      I_T     low-threshold Ca2+. Source of Ca2+ influx and
                          of post-inhibitory rebound bursts.
    PotassiumCa   I_KCa   K+ gated by intracellular Ca2+, not voltage.
                          Reads the cell's own recent activity ->
                          adaptation and the slow afterhyperpolarization.
    HCurrent      I_h     mixed cation current activated BY hyper-
                          polarization. Produces the voltage sag and
                          makes the cell a resonator.

Parameters follow the published minimal cortical models. Where a paper
gives a range, the value here is the one that reproduces the reference
behaviour, and the test suite checks that behaviour rather than the
parameter.

References:
    Connor JA, Stevens CF (1971). J Physiol 213:31-53.   [I_A]
    Pospischil M et al. (2008). "Minimal Hodgkin-Huxley type models for
        different classes of cortical and thalamic neurons."
        Biol Cybern 99:427-441.                          [I_M, I_T]
    Destexhe A, Babloyantz A, Sejnowski TJ (1993). Biophys J 65:1538.
                                                          [I_h, Ca2+]
    Traub RD, Miles R (1991). Neuronal Networks of the Hippocampus.
                                                          [I_KCa]
"""

import numpy as np


def _sig(x):
    """Numerically safe logistic, for arguments of either sign."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    e = np.exp(x[~pos])
    out[~pos] = e / (1.0 + e)
    return out if out.ndim else float(out)


class IonChannel:
    """Base class. A channel owns gating variables and a current.

    `state_names` lists the gating variables in the order the neuron
    stores them. `needs_calcium` marks channels that read intracellular
    [Ca2+]; `carries_calcium` marks those that raise it.
    """

    name = 'channel'
    state_names = ()
    needs_calcium = False
    carries_calcium = False

    def steady_state(self, V):
        """Gating values at equilibrium for a held voltage."""
        raise NotImplementedError

    def current(self, V, state, Ca=None):
        """Current in uA/cm^2, outward positive."""
        raise NotImplementedError

    def derivatives(self, V, state, Ca=None):
        """d(state)/dt for each gating variable, in 1/ms."""
        raise NotImplementedError


class PotassiumA(IonChannel):
    """Transient A-type K+ current (Connor-Stevens).

    Activates fast and inactivates slowly. At rest it is largely
    inactivated; hyperpolarization removes the inactivation, so a cell
    coming out of inhibition fires only after I_A has decayed. That delay
    is the signature to test for.
    """

    name = 'I_A'
    state_names = ('a', 'b')

    def __init__(self, g_max=47.7, E_K=-77.0):
        self.g_max = g_max
        self.E_K = E_K

    @staticmethod
    def a_inf(V):
        return (0.0761 * np.exp((V + 94.22) / 31.84)
                / (1.0 + np.exp((V + 1.17) / 28.93))) ** (1.0 / 3.0)

    @staticmethod
    def tau_a(V):
        return 0.3632 + 1.158 / (1.0 + np.exp((V + 55.96) / 20.12))

    @staticmethod
    def b_inf(V):
        return (1.0 / (1.0 + np.exp((V + 53.3) / 14.54))) ** 4.0

    @staticmethod
    def tau_b(V):
        return 1.24 + 2.678 / (1.0 + np.exp((V + 50.0) / 16.027))

    def steady_state(self, V):
        return self.a_inf(V), self.b_inf(V)

    def current(self, V, state, Ca=None):
        a, b = state
        return self.g_max * a ** 3 * b * (V - self.E_K)

    def derivatives(self, V, state, Ca=None):
        a, b = state
        return ((self.a_inf(V) - a) / self.tau_a(V),
                (self.b_inf(V) - b) / self.tau_b(V))


class PotassiumM(IonChannel):
    """Muscarinic slow K+ current (Pospischil et al. 2008).

    Non-inactivating and slow (tau up to ~1 s), so it accumulates while
    the cell fires and progressively opposes further firing. This is the
    canonical mechanism of spike frequency adaptation.
    """

    name = 'I_M'
    state_names = ('w',)

    def __init__(self, g_max=0.5, E_K=-77.0, tau_max=1000.0):
        self.g_max = g_max
        self.E_K = E_K
        self.tau_max = tau_max

    @staticmethod
    def w_inf(V):
        return _sig((V + 35.0) / 10.0)

    def tau_w(self, V):
        return self.tau_max / (3.3 * np.exp((V + 35.0) / 20.0)
                               + np.exp(-(V + 35.0) / 20.0))

    def steady_state(self, V):
        return (self.w_inf(V),)

    def current(self, V, state, Ca=None):
        return self.g_max * state[0] * (V - self.E_K)

    def derivatives(self, V, state, Ca=None):
        return ((self.w_inf(V) - state[0]) / self.tau_w(V),)


class CalciumT(IonChannel):
    """Low-threshold T-type Ca2+ current (Pospischil et al. 2008).

    De-inactivates when hyperpolarized, so releasing a cell from
    inhibition triggers a rebound burst. Also the calcium source that
    feeds I_KCa.
    """

    name = 'I_T'
    state_names = ('s', 'u')
    carries_calcium = True

    def __init__(self, g_max=0.4, E_Ca=120.0):
        self.g_max = g_max
        self.E_Ca = E_Ca

    @staticmethod
    def s_inf(V):
        return _sig((V + 57.0) / 6.2)

    @staticmethod
    def u_inf(V):
        return _sig(-(V + 81.0) / 4.0)

    @staticmethod
    def tau_s(V):
        return 0.612 + 1.0 / (np.exp(-(V + 132.0) / 16.7)
                              + np.exp((V + 16.8) / 18.2))

    @staticmethod
    def tau_u(V):
        return np.where(V < -80.0,
                        np.exp((V + 467.0) / 66.6),
                        28.0 + np.exp(-(V + 22.0) / 10.5))

    def steady_state(self, V):
        return self.s_inf(V), self.u_inf(V)

    def current(self, V, state, Ca=None):
        s, u = state
        return self.g_max * s ** 2 * u * (V - self.E_Ca)

    def derivatives(self, V, state, Ca=None):
        s, u = state
        return ((self.s_inf(V) - s) / self.tau_s(V),
                (self.u_inf(V) - u) / self.tau_u(V))


class PotassiumCa(IonChannel):
    """Ca2+-activated K+ current (Traub-Miles form).

    Gated by intracellular calcium rather than by voltage, so it reads
    the cell's own recent spiking history. Produces the slow after-
    hyperpolarization and adaptation on a slower timescale than I_M.
    """

    name = 'I_KCa'
    state_names = ()
    needs_calcium = True

    def __init__(self, g_max=1.0, E_K=-77.0, K_d=0.0005):
        self.g_max = g_max
        self.E_K = E_K
        self.K_d = K_d          # mM; half activation

    def steady_state(self, V):
        return ()

    def current(self, V, state, Ca=None):
        activation = Ca / (Ca + self.K_d)
        return self.g_max * activation * (V - self.E_K)

    def derivatives(self, V, state, Ca=None):
        return ()


class HCurrent(IonChannel):
    """Hyperpolarization-activated mixed cation current I_h.

    Runs backwards relative to every other channel: hyperpolarization
    opens it, and because E_h sits above rest it then depolarizes the
    cell. The result is the voltage sag during a hyperpolarizing step and
    a rebound on release.
    """

    name = 'I_h'
    state_names = ('r',)

    def __init__(self, g_max=0.05, E_h=-43.0):
        self.g_max = g_max
        self.E_h = E_h

    @staticmethod
    def r_inf(V):
        return _sig(-(V + 75.0) / 5.5)

    @staticmethod
    def tau_r(V):
        return 1.0 / (np.exp(-14.59 - 0.086 * V)
                      + np.exp(-1.87 + 0.0701 * V))

    def steady_state(self, V):
        return (self.r_inf(V),)

    def current(self, V, state, Ca=None):
        return self.g_max * state[0] * (V - self.E_h)

    def derivatives(self, V, state, Ca=None):
        return ((self.r_inf(V) - state[0]) / self.tau_r(V),)


class DendriticCaAP(IonChannel):
    """Dendritic calcium action potential with an amplitude that falls.

    Every other channel in this file conducts more as the membrane
    depolarizes further. This one does not. Its activation rises with
    voltage and then an inactivation term overtakes it, so the current
    peaks at an intermediate depolarization and *decreases* beyond it.

    That non-monotonicity is the point. A unit whose output rises and
    then falls with input strength can answer "yes" to a moderate input
    and "no" to a strong one, which is exclusive-or. A perceptron cannot
    compute XOR at all -- it is the canonical example of a function
    requiring more than one layer -- so a dendritic branch that can do it
    is not a summing unit with a threshold.

    Gidon et al. found exactly this in the dendrites of human cortical
    layer 2/3 pyramidal neurons: a calcium-mediated spike whose amplitude
    is maximal near threshold and declines with stronger stimulation,
    supporting XOR within a single branch.

    Modelled phenomenologically rather than from channel kinetics, since
    the underlying conductances are not fully identified:

        g = g_max * m(V) * (1 - n(V))

    with m a sigmoid opening at `V_half_act` and n a sigmoid closing at
    the higher `V_half_inact`. Between the two the branch conducts; above
    both it shuts down again.

    Reference:
        Gidon A, Zolnik TA, Fidzinski P, Bolduan F, Papoutsi A, Poirazi P,
        Holtkamp M, Vida I, Larkum ME (2020). "Dendritic action
        potentials and computation in human layer 2/3 cortical neurons."
        Science 367:83-87.
    """

    name = 'I_dCaAP'
    state_names = ('c',)

    def __init__(self, g_max=1.2, E_rev=70.0, V_half_act=-40.0,
                 k_act=4.0, V_half_inact=-25.0, k_inact=4.0, tau=2.0):
        self.g_max = g_max
        self.E_rev = E_rev
        self.V_half_act = V_half_act
        self.k_act = k_act
        self.V_half_inact = V_half_inact
        self.k_inact = k_inact
        self.tau = tau

    def activation(self, V):
        """The non-monotonic conductance profile, as a function of V.

        Rises through `V_half_act`, peaks, and is shut off again by the
        inactivation sigmoid centred on `V_half_inact`.
        """
        m = _sig((np.asarray(V, dtype=float) - self.V_half_act) / self.k_act)
        n = _sig((np.asarray(V, dtype=float) - self.V_half_inact)
                 / self.k_inact)
        return m * (1.0 - n)

    def steady_state(self, V):
        return (self.activation(V),)

    def current(self, V, state, Ca=None):
        return self.g_max * state[0] * (V - self.E_rev)

    def derivatives(self, V, state, Ca=None):
        return ((self.activation(V) - state[0]) / self.tau,)


class CalciumDynamics:
    """Intracellular Ca2+ concentration, in mM.

        d[Ca]/dt = -k * I_Ca - ([Ca] - Ca_rest) / tau

    Calcium entering through Ca2+ channels raises the concentration in a
    thin shell under the membrane; pumps clear it with time constant tau.
    The concentration is what gates I_KCa, which is how a cell's firing
    history feeds back onto its excitability.
    """

    def __init__(self, Ca_rest=1e-4, tau=200.0, k=0.003):
        self.Ca_rest = Ca_rest
        self.tau = tau
        self.k = k

    def derivative(self, Ca, I_Ca):
        return -self.k * I_Ca - (Ca - self.Ca_rest) / self.tau
