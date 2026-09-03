"""
Adaptive-step integration.

A Hodgkin-Huxley neuron spends almost all of its time doing nothing. It
rests, and then for about one millisecond in every hundred it fires. But
a fixed-step integrator runs at the speed of the fastest event in the
system at all times: dt must be ~0.01-0.05 ms because the sodium gate
has tau_m ~ 0.1 ms during the upstroke, so a cell firing at 10 Hz gets
20,000 updates per second to produce 10 events.

That is 0.05% efficiency, and it is the reason detailed simulations need
supercomputers. The cost is not the biology; it is the method.

Embedded Runge-Kutta fixes it. Two solutions of different order are
computed from the same function evaluations; their difference estimates
the local error, and the step size is adjusted to keep that error at a
requested tolerance. Where the trajectory is flat the step grows until
accuracy stops improving; during a spike it shrinks automatically. No
one has to know in advance when the spikes will happen.

Bogacki-Shampine (3,2) is used here: four stages, third-order solution
with an embedded second-order estimate, and the First Same As Last
property so a rejected step costs nothing extra and an accepted one
reuses its last evaluation as the next one's first. Three effective
function evaluations per step against RK4's four -- so the method is
cheaper per step as well as needing fewer of them.

Accuracy is controlled by `rtol`/`atol`, not by a step size chosen by
hand, which also makes it far harder to get a silently wrong answer from
a step that was slightly too large.

Reference:
    Bogacki P, Shampine LF (1989). "A 3(2) pair of Runge-Kutta formulas."
    Appl Math Lett 2:321-325.
"""

import numpy as np


class AdaptiveIntegrator:
    """Bogacki-Shampine 3(2) with PI step-size control.

    Parameters
    ----------
    rtol, atol : relative and absolute error tolerances per state
        variable. Voltage is in mV and gates are in [0, 1], so atol is
        applied per-component via `scale`.
    dt_min, dt_max : bounds on the step, in ms.
    safety : factor applied to the ideal step, below 1 so that a step
        computed as exactly on tolerance is not immediately rejected.
    """

    def __init__(self, rtol=1e-4, atol=1e-6, dt_min=1e-5, dt_max=1.0,
                 safety=0.9, min_scale=0.2, max_scale=5.0):
        self.rtol = rtol
        self.atol = atol
        self.dt_min = dt_min
        self.dt_max = dt_max
        self.safety = safety
        self.min_scale = min_scale
        self.max_scale = max_scale
        self.reset_stats()

    def reset_stats(self):
        self.n_accepted = 0
        self.n_rejected = 0
        self.n_evaluations = 0
        self.dt_history = []
        self.last_dt = None

    @property
    def stats(self):
        total = self.n_accepted + self.n_rejected
        return {'accepted': self.n_accepted, 'rejected': self.n_rejected,
                'evaluations': self.n_evaluations,
                'reject_rate': (self.n_rejected / total) if total else 0.0,
                'mean_dt': float(np.mean(self.dt_history))
                if self.dt_history else 0.0,
                'min_dt': float(np.min(self.dt_history))
                if self.dt_history else 0.0,
                'max_dt': float(np.max(self.dt_history))
                if self.dt_history else 0.0}

    # ---- one step --------------------------------------------------------

    def try_step(self, f, t, y, dt, k1=None):
        """Attempt one step. Returns (y3, error_norm, k_last, k1_used).

        `f(t, y)` returns dy/dt. `k1` is the derivative at (t, y), reused
        from the previous accepted step under the FSAL property.
        """
        if k1 is None:
            k1 = f(t, y)
            self.n_evaluations += 1
        k2 = f(t + 0.5 * dt, y + dt * 0.5 * k1)
        k3 = f(t + 0.75 * dt, y + dt * 0.75 * k2)
        self.n_evaluations += 2

        # third-order solution
        y3 = y + dt * (2.0 / 9.0 * k1 + 1.0 / 3.0 * k2 + 4.0 / 9.0 * k3)
        k4 = f(t + dt, y3)
        self.n_evaluations += 1
        # embedded second-order solution
        y2 = y + dt * (7.0 / 24.0 * k1 + 0.25 * k2
                       + 1.0 / 3.0 * k3 + 0.125 * k4)

        return y3, self.error_norm(y, y3, y2), k4

    def error_norm(self, y, y3, y2):
        """RMS of the embedded error, scaled per component.

        Factored out so a population integrator can impose a different
        norm over a matrix state without duplicating the Bogacki-Shampine
        arithmetic. For a single neuron's `y = [V, m, h, n]` this is
        unchanged from the inline form it replaces.
        """
        scale = self.atol + self.rtol * np.maximum(np.abs(y), np.abs(y3))
        return float(np.sqrt(np.mean(((y3 - y2) / scale) ** 2)))

    def integrate(self, f, y0, t_span, dt0=0.01, max_steps=10_000_000,
                  dense_dt=None):
        """Integrate f from t_span[0] to t_span[1].

        Returns (times, states). If `dense_dt` is given, the solution is
        additionally sampled on a uniform grid by linear interpolation,
        which is what a caller wanting a regular trace should use --
        the internal steps are deliberately irregular.
        """
        t0, t1 = t_span
        t, y = float(t0), np.asarray(y0, dtype=float).copy()
        dt = min(dt0, self.dt_max)
        times, states = [t], [y.copy()]
        k1 = None
        steps = 0

        while t < t1 and steps < max_steps:
            dt = min(dt, t1 - t)
            y_new, err, k_last = self.try_step(f, t, y, dt, k1)

            if err <= 1.0 or dt <= self.dt_min:
                t += dt
                y = y_new
                k1 = k_last          # FSAL: reuse as next k1
                self.n_accepted += 1
                self.dt_history.append(dt)
                times.append(t)
                states.append(y.copy())
            else:
                self.n_rejected += 1
                k1 = None            # step rejected; k1 no longer valid

            factor = self.safety * (1.0 / max(err, 1e-12)) ** (1.0 / 3.0)
            factor = np.clip(factor, self.min_scale, self.max_scale)
            dt = float(np.clip(dt * factor, self.dt_min, self.dt_max))
            steps += 1

        # Remember where the controller had settled. A caller integrating
        # in consecutive windows should resume from here rather than from
        # a cold dt0: restarting small throws away the large step the
        # controller worked up to, and in a network that is most of the
        # saving.
        self.last_dt = dt
        times = np.array(times)
        states = np.array(states)
        if dense_dt is None:
            return times, states
        grid = np.arange(t0, t1, dense_dt)
        dense = np.empty((len(grid), states.shape[1]))
        for j in range(states.shape[1]):
            dense[:, j] = np.interp(grid, times, states[:, j])
        return grid, dense


def hh_derivatives(neuron, I_of_t):
    """Wrap a HodgkinHuxleyNeuron as an f(t, y) for the integrator."""
    from .hodgkin_huxley import HodgkinHuxleyNeuron as H

    def f(t, y):
        V, m, h, n = y[0], y[1], y[2], y[3]
        I_Na = neuron.g_Na * m ** 3 * h * (V - neuron.E_Na)
        I_K = neuron.g_K * n ** 4 * (V - neuron.E_K)
        I_L = neuron.g_L * (V - neuron.E_L)
        p = neuron.phi
        I = I_of_t(t)
        return np.array([
            (I - I_Na - I_K - I_L) / neuron.C_m,
            p * (H.alpha_m(V) * (1.0 - m) - H.beta_m(V) * m),
            p * (H.alpha_h(V) * (1.0 - h) - H.beta_h(V) * h),
            p * (H.alpha_n(V) * (1.0 - n) - H.beta_n(V) * n),
        ])
    return f
