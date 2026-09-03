"""
Hodgkin-Huxley neuron — real biophysical simulation.

This is NOT a neural network engine. It does not learn, it does not
classify. It integrates the conductance-based membrane equations that
Hodgkin and Huxley fitted to squid giant axon voltage-clamp data in 1952
(Nobel Prize 1963). Parameters are measured physical quantities, not
trained weights.

Membrane equation (current balance across the lipid bilayer):

    C_m dV/dt = I_ext - I_Na - I_K - I_L

    I_Na = g_Na * m^3 * h * (V - E_Na)     fast sodium, 3 activation
                                            gates + 1 inactivation gate
    I_K  = g_K  * n^4 * (V - E_K)          delayed rectifier potassium
    I_L  = g_L  * (V - E_L)                passive leak

Each gating variable x in {m, h, n} is a probability in [0, 1] obeying
first-order kinetics between an open and a closed conformation:

    dx/dt = alpha_x(V) * (1 - x) - beta_x(V) * x

The rate functions alpha/beta are empirical fits to voltage-clamp
measurements. m activates in ~0.1 ms (spike upstroke), h inactivates in
~1 ms (spike termination + refractory period), n activates in ~5 ms
(repolarization). The separation of these three timescales is what
produces an action potential.

Units follow the original paper:
    V     mV
    t     ms
    C_m   uF/cm^2
    g     mS/cm^2
    I     uA/cm^2

Temperature: rate constants are scaled by a Q10 factor. The canonical
parameter set was measured at 6.3 C.

Reference:
    Hodgkin AL, Huxley AF (1952). "A quantitative description of membrane
    current and its application to conduction and excitation in nerve."
    J Physiol 117(4):500-544.
"""

import numpy as np


def _x_over_1_minus_exp_neg_x(x):
    """Compute x / (1 - exp(-x)), removable singularity at x = 0.

    The HH alpha_m and alpha_n rate functions are 0/0 at V = -40 and
    V = -55 respectively. The limit is 1. Naive evaluation returns NaN
    and silently destroys the simulation, so switch to the Taylor
    expansion 1 + x/2 + x^2/12 near zero.
    """
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-6
    safe = np.where(small, 1.0, x)
    return np.where(small,
                    1.0 + x / 2.0 + x * x / 12.0,
                    safe / (1.0 - np.exp(-safe)))


class HodgkinHuxleyNeuron:
    """Single-compartment conductance-based neuron.

    The whole cell is treated as one isopotential patch of membrane:
    no dendrites, no spatial extent, no axonal propagation. This is the
    standard "point neuron" abstraction and is the correct starting
    point before adding compartments.
    """

    # Canonical squid giant axon parameters (HH 1952, 6.3 C)
    def __init__(self,
                 C_m=1.0,
                 g_Na=120.0, g_K=36.0, g_L=0.3,
                 E_Na=50.0, E_K=-77.0, E_L=-54.387,
                 temperature=6.3, channels=(), calcium=None):
        self.C_m = C_m
        self.g_Na = g_Na
        self.g_K = g_K
        self.g_L = g_L
        self.E_Na = E_Na
        self.E_K = E_K
        self.E_L = E_L
        self.temperature = temperature

        # Q10 = 3 for channel gating kinetics; conductances are unscaled.
        self.phi = 3.0 ** ((temperature - 6.3) / 10.0)

        # Optional extra conductances (I_A, I_M, I_T, I_KCa, I_h ...).
        # With none of them the neuron is exactly the 1952 model, so the
        # base behaviour and its tests are untouched.
        self.channels = list(channels)
        self._n_extra = sum(len(c.state_names) for c in self.channels)
        self._needs_ca = any(c.needs_calcium for c in self.channels)
        self._has_ca_source = any(c.carries_calcium for c in self.channels)
        self._track_ca = self._needs_ca or self._has_ca_source
        if self._track_ca:
            from .channels import CalciumDynamics
            self.calcium = calcium or CalciumDynamics()
        else:
            self.calcium = calcium

        self._clamped = False
        self._clamp_current = 0.0
        self.reset()

    # ---- gating kinetics -------------------------------------------------

    @staticmethod
    def alpha_m(V):
        return 1.0 * _x_over_1_minus_exp_neg_x((V + 40.0) / 10.0)

    @staticmethod
    def beta_m(V):
        return 4.0 * np.exp(-(V + 65.0) / 18.0)

    @staticmethod
    def alpha_h(V):
        return 0.07 * np.exp(-(V + 65.0) / 20.0)

    @staticmethod
    def beta_h(V):
        return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))

    @staticmethod
    def alpha_n(V):
        return 0.1 * _x_over_1_minus_exp_neg_x((V + 55.0) / 10.0)

    @staticmethod
    def beta_n(V):
        return 0.125 * np.exp(-(V + 65.0) / 80.0)

    def steady_state(self, V):
        """Return (m_inf, h_inf, n_inf) — gate values at equilibrium for V."""
        am, bm = self.alpha_m(V), self.beta_m(V)
        ah, bh = self.alpha_h(V), self.beta_h(V)
        an, bn = self.alpha_n(V), self.beta_n(V)
        return am / (am + bm), ah / (ah + bh), an / (an + bn)

    def time_constants(self, V):
        """Return (tau_m, tau_h, tau_n) in ms — how fast each gate responds.

        The ~50x spread between tau_m and tau_n at rest is the timescale
        separation that makes the action potential possible.
        """
        f = self.phi
        return (1.0 / (f * (self.alpha_m(V) + self.beta_m(V))),
                1.0 / (f * (self.alpha_h(V) + self.beta_h(V))),
                1.0 / (f * (self.alpha_n(V) + self.beta_n(V))))

    # ---- state -----------------------------------------------------------

    def reset(self, V0=-65.0):
        """Set the neuron to rest: gates at their steady state for V0."""
        self.V = float(V0)
        self.m, self.h, self.n = self.steady_state(V0)
        self.m, self.h, self.n = float(self.m), float(self.h), float(self.n)
        self.extra = np.array(
            [float(v) for c in self.channels for v in c.steady_state(V0)])
        self.Ca = (self.calcium.Ca_rest if self._track_ca else None)
        self._clamped = False
        self._clamp_current = 0.0
        return self

    # ---- experimental manipulations --------------------------------------

    def voltage_clamp(self, V_hold):
        """Hold the membrane at V_hold, as a patch-clamp amplifier does.

        Gating variables keep evolving; only V is pinned. This is the
        instrument Hodgkin and Huxley used to measure the conductances in
        the first place, and it is the only way to observe a synaptic
        conductance without the cell's own spikes contaminating it.
        """
        self._clamped = True
        self.V = float(V_hold)
        return self

    def release_clamp(self):
        """Return to current-clamp (free-running) mode."""
        self._clamped = False
        return self

    @property
    def clamp_current(self):
        """Current the amplifier had to inject to hold V, uA/cm^2.

        Sign convention matches the recorded trace: this is the command
        current, equal and opposite to the net membrane current.
        """
        return self._clamp_current

    def block_sodium(self, fraction=1.0):
        """Block a fraction of Na+ channels — tetrodotoxin in software.

        Fully blocking abolishes spikes while leaving synaptic input
        intact, which is how subthreshold PSPs are isolated experimentally.
        """
        self.g_Na = 120.0 * (1.0 - fraction)
        return self

    def block_potassium(self, fraction=1.0):
        """Block a fraction of delayed-rectifier K+ channels (TEA)."""
        self.g_K = 36.0 * (1.0 - fraction)
        return self

    def resting_potential(self, tol=1e-9, max_iter=200):
        """Solve for the true resting potential (net ionic current = 0).

        Bisection on the steady-state I-V curve. The textbook value of
        -65 mV is approximate; the exact root depends on E_L.
        """
        def net_current(V):
            m, h, n = self.steady_state(V)
            I = (self.g_Na * m ** 3 * h * (V - self.E_Na)
                 + self.g_K * n ** 4 * (V - self.E_K)
                 + self.g_L * (V - self.E_L))
            Ca = self.calcium.Ca_rest if self._track_ca else None
            for c in self.channels:
                I += c.current(V, c.steady_state(V), Ca)
            return I

        lo, hi = -90.0, -50.0
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            if net_current(lo) * net_current(mid) <= 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < tol:
                break
        return 0.5 * (lo + hi)

    # ---- dynamics --------------------------------------------------------

    def currents(self, V, m, h, n):
        """Return (I_Na, I_K, I_L) in uA/cm^2 — the ionic current breakdown."""
        return (self.g_Na * m ** 3 * h * (V - self.E_Na),
                self.g_K * n ** 4 * (V - self.E_K),
                self.g_L * (V - self.E_L))

    def _split_extra(self, extra):
        """Slice the flat extra-state vector into one tuple per channel."""
        out, i = [], 0
        for c in self.channels:
            k = len(c.state_names)
            out.append(tuple(extra[i:i + k]))
            i += k
        return out

    def channel_currents(self, V=None, extra=None, Ca=None):
        """Current contributed by each extra channel, keyed by name."""
        V = self.V if V is None else V
        extra = self.extra if extra is None else extra
        Ca = self.Ca if Ca is None else Ca
        return {c.name: c.current(V, st, Ca)
                for c, st in zip(self.channels, self._split_extra(extra))}

    def _derivatives(self, state, I_ext):
        V, m, h, n = state[0], state[1], state[2], state[3]
        I_Na, I_K, I_L = self.currents(V, m, h, n)
        f = self.phi
        base = [
            None,   # dV/dt, filled in once the channel currents are known
            f * (self.alpha_m(V) * (1.0 - m) - self.beta_m(V) * m),
            f * (self.alpha_h(V) * (1.0 - h) - self.beta_h(V) * h),
            f * (self.alpha_n(V) * (1.0 - n) - self.beta_n(V) * n),
        ]
        if not self.channels:
            base[0] = (I_ext - I_Na - I_K - I_L) / self.C_m
            return np.array(base)

        extra = state[4:4 + self._n_extra]
        Ca = state[4 + self._n_extra] if self._track_ca else None
        I_extra, d_extra, I_Ca = 0.0, [], 0.0
        for c, st in zip(self.channels, self._split_extra(extra)):
            I_c = c.current(V, st, Ca)
            I_extra += I_c
            if c.carries_calcium:
                I_Ca += I_c
            d_extra.extend(c.derivatives(V, st, Ca))
        base[0] = (I_ext - I_Na - I_K - I_L - I_extra) / self.C_m

        out = base + list(d_extra)
        if self._track_ca:
            out.append(self.calcium.derivative(Ca, I_Ca))
        return np.array(out)

    def step(self, I_ext=0.0, dt=0.01):
        """Advance one timestep with RK4.

        dt = 0.01 ms is the standard choice. The sodium activation gate
        has tau_m ~ 0.05 ms during the upstroke, so larger steps distort
        spike amplitude; forward Euler goes unstable above ~0.04 ms.
        """
        s = np.concatenate(([self.V, self.m, self.h, self.n], self.extra,
                            [self.Ca] if self._track_ca else []))
        k1 = self._derivatives(s, I_ext)
        k2 = self._derivatives(s + 0.5 * dt * k1, I_ext)
        k3 = self._derivatives(s + 0.5 * dt * k2, I_ext)
        k4 = self._derivatives(s + dt * k3, I_ext)
        s = s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        if getattr(self, '_clamped', False):
            # The amplifier holds V and supplies whatever current that
            # takes. Record it; the gates still advance.
            I_Na, I_K, I_L = self.currents(self.V, self.m, self.h, self.n)
            self._clamp_current = I_Na + I_K + I_L - I_ext
            s[0] = self.V

        # Gates are probabilities. RK4 can overshoot by ~1e-12 at the
        # rails; clip so the invariant holds exactly.
        self.V = float(s[0])
        self.m = float(np.clip(s[1], 0.0, 1.0))
        self.h = float(np.clip(s[2], 0.0, 1.0))
        self.n = float(np.clip(s[3], 0.0, 1.0))
        if self._n_extra:
            self.extra = np.clip(s[4:4 + self._n_extra], 0.0, 1.0)
        if self._track_ca:
            # A concentration cannot go negative.
            self.Ca = float(max(s[4 + self._n_extra], 0.0))
        return self.V

    def simulate(self, I_ext, duration=None, dt=0.01, record_currents=False):
        """Run the neuron and record its trajectory.

        I_ext may be a constant (uA/cm^2), a callable t -> current, or an
        array of per-step currents. Returns a dict of numpy arrays.
        """
        if callable(I_ext):
            n_steps = int(round(duration / dt))
            current_at = lambda i: I_ext(i * dt)
        elif np.ndim(I_ext) > 0:
            I_arr = np.asarray(I_ext, dtype=float)
            n_steps = len(I_arr)
            current_at = lambda i: I_arr[i]
        else:
            n_steps = int(round(duration / dt))
            current_at = lambda i: I_ext

        t = np.empty(n_steps)
        V = np.empty(n_steps)
        m = np.empty(n_steps)
        h = np.empty(n_steps)
        n = np.empty(n_steps)
        I_rec = np.empty((n_steps, 3)) if record_currents else None

        for i in range(n_steps):
            I = current_at(i)
            t[i], V[i], m[i], h[i], n[i] = i * dt, self.V, self.m, self.h, self.n
            if record_currents:
                I_rec[i] = self.currents(self.V, self.m, self.h, self.n)
            self.step(I, dt)

        out = {'t': t, 'V': V, 'm': m, 'h': h, 'n': n, 'dt': dt}
        if record_currents:
            out['I_Na'] = I_rec[:, 0]
            out['I_K'] = I_rec[:, 1]
            out['I_L'] = I_rec[:, 2]
        return out


# ---- spike train analysis ------------------------------------------------

def detect_spikes(t, V, threshold=0.0):
    """Return spike times (ms) as upward crossings of `threshold`.

    0 mV is the conventional detection level for HH: it sits far above
    any subthreshold oscillation and well below the ~+40 mV peak, so the
    count is insensitive to the exact value.
    """
    above = V > threshold
    crossings = np.where(~above[:-1] & above[1:])[0]
    return t[crossings + 1]


def firing_rate(t, V, threshold=0.0, skip_transient=20.0):
    """Steady-state firing rate in Hz, ignoring the onset transient.

    Measured from the mean inter-spike interval rather than a spike count
    so the result does not depend on where the window happens to end.
    """
    spikes = detect_spikes(t, V, threshold)
    spikes = spikes[spikes >= skip_transient]
    if len(spikes) < 2:
        return 0.0
    return 1000.0 / np.mean(np.diff(spikes))
