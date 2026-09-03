"""
STG model neuron and graded synapses — the substrate of the pyloric circuit.

Single-compartment stomatogastric ganglion (STG) neuron of Prinz,
Billimoria & Marder (2003). Eight currents:

    I_Na    fast sodium                     m^3 h
    I_CaT   fast transient calcium          m^3 h
    I_CaS   slow transient calcium          m^3 h
    I_A     transient potassium             m^3 h
    I_KCa   calcium-dependent potassium     m^4     (m_inf reads [Ca])
    I_Kd    delayed rectifier potassium     m^4
    I_H     hyperpolarization-activated     m
    I_leak  passive leak

Gating kinetics are the voltage-clamp fits to lobster STG neurons
(Turrigiano, LeMasson & Marder 1995, as formulated by Liu, Golowasch,
Abbott & Marder 1998) exactly as tabulated in Prinz et al. 2003, Table 1.
All gates follow tau(V) dx/dt = x_inf(V) - x.

The intracellular calcium pool (Prinz 2003, from Liu 1998):

    tau_Ca dCa/dt = -f * (I_CaT + I_CaS) - Ca + Ca0

with tau_Ca = 200 ms, f = 14.96 uM/nA, Ca0 = 0.05 uM, and the calcium
reversal potential recomputed every step from the Nernst equation with
[Ca]_ext = 3 mM at T = 283 K (10 C, the standard STG bath temperature).

Chemical synapses are GRADED (Prinz, Bucher & Marder 2004, Methods
p. 1351, after Abbott & Marder, "Modeling Small Networks", 1998):
transmitter release is a continuous sigmoidal function of presynaptic
voltage, not a spike event —

    I_s = g_s * s * (V_post - E_s)
    tau_s(V_pre) ds/dt = s_inf(V_pre) - s
    s_inf = 1 / (1 + exp((V_th - V_pre)/Delta)),  tau_s = (1 - s_inf)/k_minus

with V_th = -35 mV, Delta = 5 mV for both transmitter types;
glutamatergic: E_s = -70 mV, k_minus = 1/40 ms^-1 (fast);
cholinergic:   E_s = -80 mV, k_minus = 1/100 ms^-1 (slow).

Units (house convention, and Prinz's own):
    V mV, t ms, [Ca] uM, conductances mS/cm^2, currents uA/cm^2.
    Membrane area 0.628e-3 cm^2, C_m = 1 uF/cm^2 (total C = 0.628 nF),
    so absolute synaptic strengths in nS convert to densities via
    g [mS/cm^2] = g [nS] * 1e-6 / 0.628e-3.

Integration is the exponential Euler method (Dayan & Abbott 2001), the
method Prinz et al. used, at their dt = 0.05 ms. Cross-checked line by
line against the published simulator of Prinz 2004 as reimplemented by
the Macke lab (github.com/mackelab/pyloric, MIT; used in Deistler,
Macke & Goncalves, PNAS 2022).

References:
    Prinz AA, Billimoria CP, Marder E (2003). "Alternative to hand-tuning
        conductance-based models: construction and analysis of databases
        of model neurons." J Neurophysiol 90:3998-4015.  [cell, Table 1+2]
    Prinz AA, Bucher D, Marder E (2004). "Similar network activity from
        disparate circuit parameters." Nat Neurosci 7:1345-1352.
        [synapse model p.1351, cell conductance sets Table 2]
    Liu Z, Golowasch J, Abbott LF, Marder E (1998). "A model neuron with
        activity-dependent conductances regulated by multiple calcium
        sensors." J Neurosci 18:2309-2320.  [channel kinetics, Ca pool]
    Turrigiano G, LeMasson G, Marder E (1995). J Neurosci 15:3640-3652.
        [the underlying voltage-clamp measurements]
"""

import numpy as np

# ---- physical constants (Prinz 2003/2004) --------------------------------

E_NA = 50.0        # mV
E_K = -80.0        # mV
E_H = -20.0        # mV
E_LEAK = -50.0     # mV
C_M = 1.0          # uF/cm^2
AREA = 0.628e-3    # cm^2  (Prinz 2003: 0.628e-3 cm^2, C = 0.628 nF)

TAU_CA = 200.0     # ms
CA_0 = 0.05        # uM   resting/steady calcium
CA_EXT = 3000.0    # uM   extracellular calcium (3 mM)
# f = 14.96 uM/nA converts calcium current to concentration change.
# In current-density units: I[nA] = I[uA/cm^2] * AREA * 1e3, so
# f_density = 14.96 * 0.628 = 9.3949 uM per (uA/cm^2).
F_CA = 14.96 * AREA * 1e3   # uM / (uA/cm^2)
# Nernst prefactor RT/zF for divalent Ca2+ at T = 283 K (10 C):
RT_2F = 8.314462 * 283.0 / (2.0 * 96485.332) * 1e3   # mV

# nS -> mS/cm^2 on the standard membrane area
NS_TO_MS_CM2 = 1e-6 / AREA


def _sig(V, half, slope):
    """1 / (1 + exp((V + half)/slope)) — the STG Boltzmann form.

    slope < 0 gives an activation curve, slope > 0 an inactivation curve,
    matching the sign convention of Prinz 2003 Table 1.
    """
    return 1.0 / (1.0 + np.exp((V + half) / slope))


# ---- channel kinetics — Prinz et al. 2003, Table 1 -----------------------
# Each function returns (x_inf, tau_x[ms]) for one gate at voltage V (mV).

def m_na(V):
    return _sig(V, 25.5, -5.29), 2.64 - 2.52 * _sig(V, 120.0, -25.0)


def h_na(V):
    tau = (1.34 * _sig(V, 62.9, -10.0)) * (1.5 + _sig(V, 34.9, 3.6))
    return _sig(V, 48.9, 5.18), tau


def m_cat(V):
    return _sig(V, 27.1, -7.2), 43.4 - 42.6 * _sig(V, 68.1, -20.5)


def h_cat(V):
    return _sig(V, 32.1, 5.5), 210.0 - 179.6 * _sig(V, 55.0, -16.9)


def m_cas(V):
    tau = 2.8 + 14.0 / (np.exp((V + 27.0) / 10.0) + np.exp((V + 70.0) / -13.0))
    return _sig(V, 33.0, -8.1), tau


def h_cas(V):
    tau = 120.0 + 300.0 / (np.exp((V + 55.0) / 9.0) + np.exp((V + 65.0) / -16.0))
    return _sig(V, 60.0, 6.2), tau


def m_a(V):
    return _sig(V, 27.2, -8.7), 23.2 - 20.8 * _sig(V, 32.9, -15.2)


def h_a(V):
    return _sig(V, 56.9, 4.9), 77.2 - 58.4 * _sig(V, 38.9, -26.5)


def m_kca(V, Ca):
    """The KCa activation reads intracellular calcium: the m_inf Boltzmann
    is scaled by Ca/(Ca + 3 uM), so the channel is shut at resting [Ca]."""
    inf = (Ca / (Ca + 3.0)) * _sig(V, 28.3, -12.6)
    return inf, 180.6 - 150.2 * _sig(V, 46.0, -22.7)


def m_kd(V):
    return _sig(V, 12.3, -11.8), 14.4 - 12.8 * _sig(V, 28.3, -19.2)


def m_h(V):
    tau = 2.0 / (np.exp(-14.59 - 0.086 * V) + np.exp(-1.87 + 0.0701 * V))
    return _sig(V, 75.0, 5.5), tau


# ---- the model cell -------------------------------------------------------

# Column order used throughout for conductance vectors.
CONDUCTANCE_NAMES = ('Na', 'CaT', 'CaS', 'A', 'KCa', 'Kd', 'H', 'leak')

# Maximal conductance densities (mS/cm^2) of the 16 model neurons used in
# the pyloric network database — Prinz, Bucher & Marder 2004, Table 2.
# Verbatim:                Na    CaT  CaS   A   KCa   Kd    H     leak
PRINZ_TABLE2 = {
    'AB/PD 1': (400.0,  2.5, 6.0, 50.0, 10.0, 100.0, 0.01, 0.00),
    'AB/PD 2': (100.0,  2.5, 6.0, 50.0,  5.0, 100.0, 0.01, 0.00),
    'AB/PD 3': (200.0,  2.5, 4.0, 50.0,  5.0,  50.0, 0.01, 0.00),
    'AB/PD 4': (200.0,  5.0, 4.0, 40.0,  5.0, 125.0, 0.01, 0.00),
    'AB/PD 5': (300.0,  2.5, 2.0, 10.0,  5.0, 125.0, 0.01, 0.00),
    'LP 1':    (100.0,  0.0, 8.0, 40.0,  5.0,  75.0, 0.05, 0.02),
    'LP 2':    (100.0,  0.0, 6.0, 30.0,  5.0,  50.0, 0.05, 0.02),
    'LP 3':    (100.0,  0.0, 10.0, 50.0, 5.0, 100.0, 0.00, 0.03),
    'LP 4':    (100.0,  0.0, 4.0, 20.0,  0.0,  25.0, 0.05, 0.03),
    'LP 5':    (100.0,  0.0, 6.0, 30.0,  0.0,  50.0, 0.03, 0.02),
    'PY 1':    (100.0,  2.5, 2.0, 50.0,  0.0, 125.0, 0.05, 0.01),
    'PY 2':    (200.0,  7.5, 0.0, 50.0,  0.0,  75.0, 0.05, 0.00),
    'PY 3':    (200.0, 10.0, 0.0, 50.0,  0.0, 100.0, 0.03, 0.00),
    'PY 4':    (400.0,  2.5, 2.0, 50.0,  0.0,  75.0, 0.05, 0.00),
    'PY 5':    (500.0,  2.5, 2.0, 40.0,  0.0, 125.0, 0.01, 0.03),
    'PY 6':    (500.0,  2.5, 2.0, 40.0,  0.0, 125.0, 0.00, 0.02),
}

# Intrinsic burst periods of the five pacemaker models when isolated,
# as published in Prinz 2004, Methods p. 1351.
ABPD_INTRINSIC_PERIOD_S = {
    'AB/PD 1': 1.46, 'AB/PD 2': 1.49, 'AB/PD 3': 1.58,
    'AB/PD 4': 1.61, 'AB/PD 5': 1.64,
}


class STGNeuron:
    """One single-compartment STG model cell: a bag of eight maximal
    conductances. The kinetics are shared by every cell; the conductance
    set is what distinguishes a pacemaker from a follower (Prinz 2003).
    """

    def __init__(self, g, name='STG'):
        g = np.asarray(g, dtype=float)
        if g.shape != (8,):
            raise ValueError("g must be 8 conductances: %s"
                             % (CONDUCTANCE_NAMES,))
        self.g = g
        self.name = name

    @classmethod
    def from_table(cls, name):
        """A cell from Prinz 2004, Table 2 — e.g. 'AB/PD 2', 'LP 4', 'PY 1'."""
        return cls(PRINZ_TABLE2[name], name=name)

    def simulate(self, duration, dt=0.05, V0=-65.0, I_ext=0.0):
        """Simulate the isolated cell. Thin wrapper over STGNetwork."""
        net = STGNetwork([self])
        return net.simulate(duration, dt=dt, V0=V0, I_ext=I_ext)


class GradedSynapse:
    """Graded inhibitory chemical synapse (Prinz 2004, Methods p. 1351).

    Continuous transmitter release: s tracks a sigmoid of presynaptic
    voltage with a release-dependent time constant. No spike detection
    is involved — subthreshold slow-wave depolarization releases
    transmitter, which is how the biological pyloric synapses work.

    g is given in nS (the unit Prinz scanned: 0/1/3/10/30/100 nS) and
    converted to a conductance density on the standard membrane area.
    """

    V_TH = -35.0    # mV, half-activation of release
    DELTA = 5.0     # mV, slope of the release sigmoid

    def __init__(self, pre, post, g_nS, E_s, k_minus, label='syn'):
        self.pre = int(pre)
        self.post = int(post)
        self.g_nS = float(g_nS)
        self.g = float(g_nS) * NS_TO_MS_CM2   # mS/cm^2
        self.E_s = float(E_s)
        self.k_minus = float(k_minus)         # 1/ms
        self.label = label

    @classmethod
    def glutamatergic(cls, pre, post, g_nS, label='glut'):
        """Fast graded inhibition (AB, LP, PY are glutamatergic).
        E_s = -70 mV, k_minus = 1/40 ms^-1 (Prinz 2004)."""
        return cls(pre, post, g_nS, E_s=-70.0, k_minus=1.0 / 40.0, label=label)

    @classmethod
    def cholinergic(cls, pre, post, g_nS, label='chol'):
        """Slow graded inhibition (PD is cholinergic).
        E_s = -80 mV, k_minus = 1/100 ms^-1 (Prinz 2004)."""
        return cls(pre, post, g_nS, E_s=-80.0, k_minus=1.0 / 100.0, label=label)

    def s_inf(self, V_pre):
        return 1.0 / (1.0 + np.exp((self.V_TH - V_pre) / self.DELTA))

    def tau_s(self, V_pre):
        return (1.0 - self.s_inf(V_pre)) / self.k_minus


class STGNetwork:
    """A set of STG cells coupled by graded synapses, integrated together.

    Exponential Euler throughout (the integrator Prinz used): V is
    advanced against the instantaneous conductance-weighted equilibrium,
    every gate against its own (x_inf, tau) at the previous voltage, and
    the synaptic activation s against (s_inf, tau_s) of the previous
    presynaptic voltage. All updates read previous-step state only.
    """

    def __init__(self, cells, synapses=()):
        self.cells = list(cells)
        self.synapses = list(synapses)
        self.n = len(self.cells)
        self.G = np.array([c.g for c in self.cells])   # (n, 8)
        for s in self.synapses:
            if not (0 <= s.pre < self.n and 0 <= s.post < self.n):
                raise ValueError("synapse endpoints out of range")
        self._pre = np.array([sy.pre for sy in self.synapses], dtype=int)
        self._post = np.array([sy.post for sy in self.synapses], dtype=int)
        self._g_syn = np.array([sy.g for sy in self.synapses])
        self._E_syn = np.array([sy.E_s for sy in self.synapses])
        self._k_min = np.array([sy.k_minus for sy in self.synapses])

    def steady_state(self, V, Ca=CA_0):
        """All 11 gate values at equilibrium for voltage V, calcium Ca."""
        return np.array([
            m_na(V)[0], m_cat(V)[0], m_cas(V)[0], m_a(V)[0],
            m_kca(V, Ca)[0], m_kd(V)[0], m_h(V)[0],
            h_na(V)[0], h_cat(V)[0], h_cas(V)[0], h_a(V)[0],
        ])

    def reset_state(self, V0=-65.0, rng=None, v0_jitter=0.0):
        """Prepare step-wise integration: V at V0 (optionally jittered),
        calcium at rest, gates at their steady state, synapses fully
        deactivated. Called by simulate(); call it directly when driving
        the network with step() — e.g. a plasticity rule that needs
        per-step access to voltages and spikes.

        V0 : scalar or per-cell array, mV. With v0_jitter > 0, each cell
             starts at V0 + U(-jitter, +jitter) drawn from `rng` — used
             to check that a rhythm is an attractor, not an initial
             condition.
        """
        n = self.n
        V = np.full(n, V0, dtype=float) if np.ndim(V0) == 0 \
            else np.asarray(V0, dtype=float).copy()
        if v0_jitter > 0.0:
            rng = rng or np.random.default_rng()
            V = V + rng.uniform(-v0_jitter, v0_jitter, size=n)
        self.V = V
        self.Ca = np.full(n, CA_0)
        self.gates = self.steady_state(V, self.Ca)     # (11, n)
        self.s = np.zeros(len(self.synapses))
        self._V_prev = None       # previous sample, for peak detection
        self._i = 0               # index of the current sample self.V
        self._last_peak = np.full(n, -np.inf)   # last raw peak time, ms
        return self

    def step(self, dt, I_ext=0.0):
        """Advance the whole network one exponential-Euler step.

        I_ext : scalar or per-cell (n,) array, uA/cm^2 — external
            current added to this step's voltage update. This is the
            hook for a caller that computes its own synaptic or stimulus
            currents each step: convert them to current density and pass
            them here (1 nA on the standard membrane area is
            1e-3 / 0.628e-3 = 1.592 uA/cm^2).

        Returns the spike events confirmed by this step as a list of
        (cell_index, peak_time_ms). A spike is a local maximum of V
        above SPIKE_THRESHOLD (same criterion as detect_spikes), so a
        peak at sample i is only confirmed while computing sample i+1 —
        events arrive one step late, the price of causal peak detection.

        All updates read previous-step state only, so the step is
        order-independent and simulate() built on it is bit-identical
        to a monolithic loop.
        """
        V, Ca, gates, s = self.V, self.Ca, self.gates, self.s
        n, m = self.n, len(self.synapses)
        gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gL = self.G.T

        (mNa, mCaT, mCaS, mA, mKCa, mKd, mH,
         hNa, hCaT, hCaS, hA) = gates

        # Channel conductances at previous-step gates, mS/cm^2
        cNa = gNa * mNa ** 3 * hNa
        cCaT = gCaT * mCaT ** 3 * hCaT
        cCaS = gCaS * mCaS ** 3 * hCaS
        cA = gA * mA ** 3 * hA
        cKCa = gKCa * mKCa ** 4
        cKd = gKd * mKd ** 4
        cH = gH * mH

        E_Ca = RT_2F * np.log(CA_EXT / Ca)             # Nernst, per cell

        # Synaptic conductance aggregated onto each postsynaptic cell
        g_sum = np.zeros(n)
        gE_sum = np.zeros(n)
        if m:
            gs = self._g_syn * s
            np.add.at(g_sum, self._post, gs)
            np.add.at(gE_sum, self._post, gs * self._E_syn)

        # Exponential Euler on V: dV/dt = (gE_tot - g_tot V + I)/C
        I = np.broadcast_to(np.asarray(I_ext, dtype=float), (n,))
        g_tot = (g_sum + cNa + cCaT + cCaS + cA + cKCa + cKd + cH + gL)
        gE_tot = (gE_sum + cNa * E_NA + (cCaT + cCaS) * E_Ca
                  + (cA + cKCa + cKd) * E_K + cH * E_H + gL * E_LEAK
                  + I)
        V_inf = gE_tot / g_tot
        V_new = V_inf + (V - V_inf) * np.exp(-dt * g_tot / C_M)

        # Calcium pool: exponential update toward Ca_inf(I_Ca(V_prev))
        I_Ca = (cCaT + cCaS) * (V - E_Ca)              # uA/cm^2
        Ca_inf = CA_0 - F_CA * I_Ca
        Ca_new = Ca_inf + (Ca - Ca_inf) * np.exp(-dt / TAU_CA)

        # Gates: exponential update toward x_inf(V_prev)
        new_gates = np.empty_like(gates)
        for j, (fn, needs_ca) in enumerate(_GATE_FNS):
            inf, tau = (fn(V, Ca) if needs_ca else fn(V))
            new_gates[j] = inf + (gates[j] - inf) * np.exp(-dt / tau)

        # Graded synapses: exponential update toward s_inf(V_pre_prev)
        if m:
            sinf = 1.0 / (1.0 + np.exp((GradedSynapse.V_TH - V[self._pre])
                                       / GradedSynapse.DELTA))
            tau_s = (1.0 - sinf) / self._k_min
            # tau_s -> 0 when release saturates; the exponential step
            # goes to s_inf exactly, no truncation needed.
            with np.errstate(divide='ignore'):
                self.s = sinf + (s - sinf) * np.exp(-dt / tau_s)

        # Spike events: peak at the current sample, confirmed by V_new.
        # Pure comparisons — the numerics above are untouched.
        events = []
        if self._V_prev is not None:
            peaks = ((V > SPIKE_THRESHOLD) & (V - self._V_prev > 0.0)
                     & (V_new - V <= 0.0))
            if peaks.any():
                t_pk = self._i * dt
                for c in np.where(peaks)[0]:
                    if t_pk - self._last_peak[c] >= SPIKE_MIN_SEPARATION:
                        events.append((int(c), t_pk))
                    self._last_peak[c] = t_pk

        self._V_prev = V
        self.V = V_new
        self.Ca = np.maximum(Ca_new, 1e-9)
        self.gates = new_gates
        self._i += 1
        return events

    def simulate(self, duration, dt=0.05, V0=-65.0, I_ext=0.0,
                 record_s=False, rng=None, v0_jitter=0.0):
        """Integrate for `duration` ms. A thin recording loop over
        reset_state() + step().

        V0, rng, v0_jitter : see reset_state().
        I_ext : scalar, per-cell array, or (n_steps, n) array, uA/cm^2.

        Returns dict with 't' (ms), 'V' (n_steps, n), 'Ca' (n_steps, n),
        and 's' (n_steps, n_syn) if record_s. Sample i is the state
        BEFORE step i, matching t = i * dt.
        """
        n, dt = self.n, float(dt)
        n_steps = int(round(duration / dt))

        if np.ndim(I_ext) == 2:
            I_of = np.asarray(I_ext, dtype=float)
            current_at = lambda i: I_of[i]
        else:
            I_const = np.broadcast_to(np.asarray(I_ext, float), (n,)).copy()
            current_at = lambda i: I_const

        self.reset_state(V0=V0, rng=rng, v0_jitter=v0_jitter)
        m = len(self.synapses)
        t_out = np.arange(n_steps) * dt
        V_out = np.empty((n_steps, n))
        Ca_out = np.empty((n_steps, n))
        s_out = np.empty((n_steps, m)) if record_s else None

        for i in range(n_steps):
            V_out[i] = self.V
            Ca_out[i] = self.Ca
            if record_s:
                s_out[i] = self.s
            self.step(dt, current_at(i))

        out = {'t': t_out, 'V': V_out, 'Ca': Ca_out, 'dt': dt}
        if record_s:
            out['s'] = s_out
        return out


_GATE_FNS = (
    (m_na, False), (m_cat, False), (m_cas, False), (m_a, False),
    (m_kca, True), (m_kd, False), (m_h, False),
    (h_na, False), (h_cat, False), (h_cas, False), (h_a, False),
)


SPIKE_THRESHOLD = -10.0        # mV — STG spike criterion (see below)
SPIKE_MIN_SEPARATION = 3.0     # ms — merge window for double maxima


def detect_spikes(t, V, threshold=SPIKE_THRESHOLD,
                  min_separation=SPIKE_MIN_SEPARATION):
    """Spike times as local maxima of V above `threshold` (mV).

    STG spikes ride on a slow wave that can itself cross a fixed level,
    so upward-crossing detection (the HH convention) merges plateau and
    spikes. Local-maximum detection above -10 mV is the convention of the
    published pyloric analyses (mackelab/pyloric summary statistics).
    Maxima closer than `min_separation` ms are merged (keep the first).
    """
    V = np.asarray(V)
    peak = ((V[1:-1] > threshold)
            & (np.diff(V[:-1]) > 0) & (np.diff(V[1:]) <= 0))
    times = t[np.where(peak)[0] + 1]
    if len(times) < 2:
        return times
    keep = np.concatenate(([True], np.diff(times) >= min_separation))
    return times[keep]
